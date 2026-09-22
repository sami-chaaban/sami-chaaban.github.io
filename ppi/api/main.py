"""
FastAPI backend for the protein-protein interaction demo.

This backend exposes routes for listing chains, analyzing interfaces, and
returning reports. The analysis logic uses PDBe Arpeggio to derive
atom-level interaction contacts.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Optional
import asyncio
import codecs
import gc
import gzip
import hashlib
import json
import os
import re
import select
import signal
import subprocess
import tempfile
import threading
import time
from pathlib import Path
import sys
from urllib import error as urlerror
from urllib import parse as urlparse
from urllib import request as urlrequest

from fastapi import FastAPI, HTTPException, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from starlette.middleware.gzip import GZipMiddleware

from .analysis import (
    cache_key,
    fetch_mmcif,
    list_chains,
)
from .cache import ReportCache
from .analysis_worker import (
    AnalysisBusy, AnalysisDisconnected, BoundedAnalysisRunner,
    ReportDelivery, analyze_and_serialize, analyze_stored_and_serialize, prepare_delivery,
)
from .structure_store import StructureStore, StructureGone, StructureCapacity
from .explain import explain_report
from .models import (
    AnalyzeRequest,
    ChainsRequest,
    ChapiMeshRequest,
    ExplainRequest,
    LocalCompanionJsonRequest,
    RibbonRequest,
)


def _env_positive_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        value = int(str(raw).strip())
        return value if value > 0 else default
    except Exception:
        return default


def _env_nonnegative_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        value = int(str(raw).strip())
        return value if value >= 0 else default
    except Exception:
        return default


def _env_enabled(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    token = str(raw).strip().lower()
    if token in {"1", "true", "yes", "on"}:
        return True
    if token in {"0", "false", "no", "off"}:
        return False
    return default


DEMO_CHAINS = {
    "4hhb": ["A", "B"],
    "1a3n": ["A", "B"],
}

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_URL_TEMPLATE = "https://rest.uniprot.org/uniprotkb/{accession}.json"
PD_BE_SUMMARY_URL_TEMPLATE = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id}"
PD_BE_MOLECULES_URL_TEMPLATE = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id}"
RCSB_ENTRY_URL_TEMPLATE = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
UNIPROT_SEARCH_FETCH_SIZE = 15
PDB_SUMMARY_TTL_SECONDS = 60 * 60 * 24
PDB_SUMMARY_MAX_WORKERS = 8
UNIPROT_SEARCH_FIELDS = ",".join(
    [
        "accession",
        "id",
        "gene_names",
        "organism_name",
        "length",
        "cc_function",
        "xref_pdb",
        "xref_alphafolddb",
    ]
)
FIRST_SENTENCE_REGEX = re.compile(r"(?<=[.!?])\s+")
PDB_SUMMARY_CACHE: dict[str, tuple[float, dict[str, Any]]] = {}
PDB_SUMMARY_CACHE_LOCK = threading.Lock()
ALPHAFOLD_REFERENCE = {
    "referenceAuthors": "Jumper et al.",
    "referenceJournal": "Nature",
    "referenceYear": "2021",
}


app = FastAPI(title="PPI Demo API", description="Backend for protein interface demo")

# Allow the browser front-end to call the API from any origin. In production,
# restrict this to your known domains.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)
# Mesh responses contain megabytes of repeated numeric JSON. Level 1 retains
# every coordinate while reducing transfer size without expensive compression.
app.add_middleware(GZipMiddleware, minimum_size=1024, compresslevel=1)

REPORT_CACHE_LOW_MEMORY = _env_enabled("CHAPI_LOW_MEMORY_MODE", _env_enabled("RENDER", False))
# Shared admission includes mesh preprocessing and scientific workers, not just
# native Coot. A completed analysis releases admission only after worker exit.
HEAVY_WORK_GATE = threading.BoundedSemaphore(_env_positive_int('HEAVY_WORKERS', 1))
structure_store = StructureStore(
    max_bytes=_env_positive_int('STRUCTURE_STORE_MAX_BYTES', 256 * 1024**2),
    max_entries=_env_positive_int('STRUCTURE_STORE_MAX_ENTRIES', 8),
    ttl_seconds=_env_positive_int('STRUCTURE_STORE_TTL_SECONDS', 1800),
    max_upload_bytes=_env_positive_int('STRUCTURE_UPLOAD_MAX_BYTES', 64 * 1024**2),
)
STRUCTURE_UPLOAD_SLOTS = asyncio.Semaphore(2)
cache = ReportCache(max_bytes=_env_nonnegative_int(
    "ANALYSIS_REPORT_CACHE_MAX_BYTES", (32 if REPORT_CACHE_LOW_MEMORY else 128) * 1024 * 1024,
))
report_store = ReportCache(
    ttl_seconds=_env_positive_int("REPORT_STORE_TTL_SECONDS", 60 * 60 * 6),
    max_entries=_env_positive_int("REPORT_STORE_MAX_ENTRIES", 64),
    max_bytes=_env_nonnegative_int("REPORT_STORE_MAX_BYTES", (16 if REPORT_CACHE_LOW_MEMORY else 64) * 1024 * 1024),
)
analysis_runner = BoundedAnalysisRunner(
    workers=_env_positive_int("ANALYZE_WORKERS", 1),
    max_pending=_env_positive_int("ANALYZE_MAX_PENDING", 8),
    recycle_workers=_env_enabled('ANALYZE_RECYCLE_WORKERS', True),
    execution_gate=HEAVY_WORK_GATE,
)

CHAPI_PYTHON = os.environ.get("CHAPI_PYTHON") or sys.executable
CHAPI_PREFIX = os.environ.get("COOT_PREFIX") or os.environ.get("CONDA_PREFIX") or ""
CHAPI_BRIDGE = Path(__file__).with_name("chapi_bridge.py")
CHAPI_LOW_MEMORY_MODE = _env_enabled("CHAPI_LOW_MEMORY_MODE", _env_enabled("RENDER", False))
STRUCTURE_TEXT_CACHE = ReportCache(
    ttl_seconds=_env_positive_int("STRUCTURE_TEXT_CACHE_TTL_SECONDS", 60 * 60 * 6),
    max_entries=_env_positive_int("STRUCTURE_TEXT_CACHE_MAX_ENTRIES", 1 if CHAPI_LOW_MEMORY_MODE else 24),
)
CHAPI_MESH_CACHE = ReportCache(
    ttl_seconds=_env_positive_int("CHAPI_MESH_CACHE_TTL_SECONDS", 60 * 30),
    max_entries=_env_positive_int("CHAPI_MESH_CACHE_MAX_ENTRIES", 1 if CHAPI_LOW_MEMORY_MODE else 8),
)
CHAPI_MESH_CACHE_MAX_BYTES = _env_nonnegative_int(
    "CHAPI_MESH_CACHE_MAX_BYTES",
    0 if CHAPI_LOW_MEMORY_MODE else 64 * 1024 * 1024,
)
CHAPI_MESH_INFLIGHT: dict[str, threading.Event] = {}
CHAPI_MESH_INFLIGHT_LOCK = threading.Lock()
CHAPI_PERSISTENT_WORKER = _env_enabled("CHAPI_PERSISTENT_WORKER", not CHAPI_LOW_MEMORY_MODE)
CHAPI_WORKER_WARMUP_AT_STARTUP = _env_enabled("CHAPI_WORKER_WARMUP_AT_STARTUP", not CHAPI_LOW_MEMORY_MODE)
CHAPI_WORKER_TIMEOUT_SECONDS = max(10, _env_positive_int("CHAPI_WORKER_TIMEOUT_SECONDS", 180))
CHAPI_SPLIT_CHAIN_BATCH_LIMIT = _env_nonnegative_int(
    "CHAPI_SPLIT_CHAIN_BATCH_LIMIT",
    1 if CHAPI_LOW_MEMORY_MODE else 0,
)
CHAPI_WORKER_PROC: Optional[subprocess.Popen] = None
CHAPI_WORKER_LOCK = threading.Lock()
# The HTTP handler runs off the event loop, but native jobs stay serialized to
# retain the existing memory bound, including in one-shot/low-memory mode.
CHAPI_MESH_EXECUTION_LOCK = threading.Lock()
# Remember successful reader recovery without retaining coordinate text.
CHAPI_MMDB_READER_CACHE = ReportCache(ttl_seconds=60 * 30, max_entries=64)


@app.get("/")
async def root():
    """Basic root endpoint so hitting the API base URL is informative."""
    return {
        "name": "PPI Demo API",
        "status": "ok",
        "routes": [
            "/chains",
            "/analyze",
            "/chapi-mesh",
            "/explain",
            "/health",
            "/local-companion-json",
            "/protein-search",
            "/protein-structures/{accession}",
        ],
    }


@app.get("/health")
async def health():
    return {"status": "ok"}


@app.post('/structures')
async def upload_structure(request: Request, format: str = 'mmcif'):
    if format not in {'mmcif', 'pdb'}:
        raise HTTPException(status_code=400, detail='Unsupported structure format')
    length = request.headers.get('content-length')
    if length and length.isdigit() and int(length) > structure_store.max_upload_bytes:
        raise HTTPException(status_code=413, detail='Structure exceeds the upload limit')
    async with STRUCTURE_UPLOAD_SLOTS:
        path = structure_store.upload_path()
        size = 0
        digest = hashlib.sha256((format + '\n').encode())
        decoder = codecs.getincrementaldecoder('utf-8')()
        try:
            with path.open('wb') as stream:
                async for chunk in request.stream():
                    size += len(chunk)
                    if size > structure_store.max_upload_bytes:
                        raise HTTPException(status_code=413, detail='Structure exceeds the upload limit')
                    decoder.decode(chunk)
                    digest.update(chunk)
                    await asyncio.to_thread(stream.write, chunk)
                decoder.decode(b'', final=True)
            if not size:
                raise HTTPException(status_code=400, detail='Structure upload is empty')
            entry = await asyncio.to_thread(structure_store.register, path, format, digest.hexdigest())
            return {'structureId': entry.identifier, 'format': format,
                    'expiresIn': structure_store.ttl_seconds}
        except UnicodeDecodeError as exc:
            raise HTTPException(status_code=400, detail='Coordinates must be UTF-8 text') from exc
        except StructureCapacity as exc:
            raise HTTPException(status_code=503, detail=str(exc), headers={'Retry-After': '1'}) from exc
        finally:
            path.unlink(missing_ok=True)


@app.delete('/structures/{structure_id}', status_code=204)
def delete_structure(structure_id: str):
    structure_store.delete(structure_id)
    return Response(status_code=204)


def _lease_structure(identifier):
    try:
        return structure_store.acquire(identifier)
    except StructureGone as exc:
        raise HTTPException(status_code=410, detail=str(exc)) from exc


@app.post("/local-companion-json")
async def local_companion_json(request: LocalCompanionJsonRequest):
    """Return a same-directory prediction JSON for a local structure path."""
    raw_structure_path = str(request.structurePath or "").strip()
    if not raw_structure_path:
        raise HTTPException(status_code=400, detail="structurePath is required")

    try:
        structure_path = Path(raw_structure_path).expanduser().resolve(strict=False)
    except Exception as exc:
        raise HTTPException(status_code=400, detail="Invalid structurePath") from exc

    if not structure_path.is_absolute() or not structure_path.is_file():
        raise HTTPException(status_code=404, detail="Structure path not found")

    directory = structure_path.parent
    candidate_names = request.candidateNames or []
    for raw_name in candidate_names[:12]:
        name = Path(str(raw_name or "").strip()).name
        if not name or Path(name).suffix.lower() != ".json":
            continue
        candidate_path = (directory / name).resolve(strict=False)
        try:
            candidate_path.relative_to(directory)
        except ValueError:
            continue
        if not candidate_path.is_file():
            continue
        try:
            text = candidate_path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            text = candidate_path.read_text(encoding="utf-8-sig")
        except OSError as exc:
            raise HTTPException(status_code=500, detail=f"Unable to read companion JSON: {exc}") from exc
        return Response(
            content=text,
            media_type="application/json",
            headers={"X-Roami-Local-Prediction-Filename": candidate_path.name},
        )

    raise HTTPException(status_code=404, detail="Companion prediction JSON not found")


@app.on_event("startup")
async def startup_chapi_worker() -> None:
    if not CHAPI_PERSISTENT_WORKER or not CHAPI_WORKER_WARMUP_AT_STARTUP:
        return
    try:
        with CHAPI_WORKER_LOCK:
            _start_chapi_worker_locked()
    except Exception:
        # Keep startup resilient: worker failures should fall back to one-shot bridge calls.
        pass


@app.on_event("shutdown")
async def shutdown_chapi_worker() -> None:
    analysis_runner.shutdown()
    with CHAPI_WORKER_LOCK:
        _stop_chapi_worker_locked()


def collapse_whitespace(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def dedupe_keep_order(values: list[str]) -> list[str]:
    output: list[str] = []
    seen: set[str] = set()
    for raw in values:
        value = collapse_whitespace(raw)
        if not value:
            continue
        token = value.upper()
        if token in seen:
            continue
        seen.add(token)
        output.append(value)
    return output


def parse_publication_year(raw: Any) -> str:
    if isinstance(raw, (int, float)) and not isinstance(raw, bool):
        year = int(raw)
        return str(year) if year > 0 else ""
    token = collapse_whitespace(raw)
    if not token:
        return ""
    match = re.search(r"\b(19|20)\d{2}\b", token)
    if match:
        return match.group(0)
    if token.isdigit():
        year = int(token)
        return str(year) if year > 0 else ""
    return ""


def normalize_author_list(raw: Any) -> str:
    values: list[str] = []
    if isinstance(raw, str):
        token = collapse_whitespace(raw)
        if token:
            values.append(token)
    elif isinstance(raw, list):
        for row in raw:
            if isinstance(row, str):
                token = collapse_whitespace(row)
                if token:
                    values.append(token)
                continue
            if not isinstance(row, dict):
                continue
            token = collapse_whitespace(row.get("name") or row.get("fullName") or row.get("value"))
            if token:
                values.append(token)
    elif isinstance(raw, dict):
        token = collapse_whitespace(raw.get("name") or raw.get("fullName") or raw.get("value"))
        if token:
            values.append(token)
    return ", ".join(dedupe_keep_order(values))


def normalize_reference_fields(raw: Any) -> dict[str, str]:
    if not isinstance(raw, dict):
        return {
            "referenceAuthors": "",
            "referenceJournal": "",
            "referenceYear": "",
        }
    return {
        "referenceAuthors": normalize_author_list(
            raw.get("referenceAuthors")
            or raw.get("authors")
            or raw.get("authorList")
            or raw.get("author_list")
            or raw.get("entry_authors")
            or raw.get("entryAuthors")
            or raw.get("rcsb_authors")
            or raw.get("rcsbAuthors")
        ),
        "referenceJournal": collapse_whitespace(
            raw.get("referenceJournal")
            or raw.get("journal")
            or raw.get("journal_abbrev")
            or raw.get("journalAbbrev")
            or raw.get("journal_name")
            or raw.get("journalName")
            or raw.get("journal_full")
            or raw.get("rcsb_journal_abbrev")
            or raw.get("rcsbJournalAbbrev")
        ),
        "referenceYear": parse_publication_year(
            raw.get("referenceYear")
            or raw.get("year")
            or raw.get("publication_year")
            or raw.get("publicationYear")
            or raw.get("journal_year")
            or raw.get("journalYear")
        ),
    }


ORGANISM_VALUE_KEYS = (
    "organism",
    "organism_name",
    "organism_scientific_name",
    "scientific_name",
    "ncbi_scientific_name",
    "pdbx_organism_scientific",
    "pdbx_gene_src_scientific_name",
    "gene_src_scientific_name",
)

ORGANISM_CONTAINER_KEYS = (
    "source",
    "sources",
    "entity_src_nat",
    "entity_src_gen",
    "pdbx_entity_src_syn",
    "rcsb_entity_source_organism",
    "polymer_entities",
    "entity",
)


def _collect_organism_tokens(raw: Any, output: list[str]) -> None:
    if raw is None:
        return
    if isinstance(raw, str):
        token = collapse_whitespace(raw)
        if token and token.upper() not in {"-", "—", "N/A", "NA", "UNKNOWN"}:
            output.append(token)
        return
    if isinstance(raw, list):
        for row in raw:
            _collect_organism_tokens(row, output)
        return
    if not isinstance(raw, dict):
        return
    for key in ORGANISM_VALUE_KEYS:
        if key in raw:
            _collect_organism_tokens(raw.get(key), output)
    for key in ORGANISM_CONTAINER_KEYS:
        if key in raw:
            _collect_organism_tokens(raw.get(key), output)


def normalize_organism_list(raw: Any) -> list[str]:
    tokens: list[str] = []
    _collect_organism_tokens(raw, tokens)
    return dedupe_keep_order(tokens)


def merge_organism_lists(primary: Any, fallback: Any) -> list[str]:
    merged = normalize_organism_list(primary)
    merged.extend(normalize_organism_list(fallback))
    return dedupe_keep_order(merged)


def merge_reference_fields(primary: Any, fallback: Any) -> dict[str, str]:
    left = normalize_reference_fields(primary)
    right = normalize_reference_fields(fallback)
    return {
        "referenceAuthors": left["referenceAuthors"] or right["referenceAuthors"],
        "referenceJournal": left["referenceJournal"] or right["referenceJournal"],
        "referenceYear": left["referenceYear"] or right["referenceYear"],
    }


def has_complete_reference_fields(raw: Any) -> bool:
    fields = normalize_reference_fields(raw)
    return bool(
        fields["referenceAuthors"]
        and fields["referenceJournal"]
        and fields["referenceYear"]
    )


def first_sentence(text: str, max_chars: int = 180) -> str:
    normalized = collapse_whitespace(text)
    if not normalized:
        return ""
    pieces = FIRST_SENTENCE_REGEX.split(normalized, maxsplit=1)
    candidate = pieces[0] if pieces else normalized
    if len(candidate) <= max_chars:
        return candidate
    return f"{candidate[: max_chars - 1].rstrip()}…"


def fetch_json(url: str, timeout: float = 18.0) -> dict[str, Any]:
    req = urlrequest.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "ppi-demo-api/1.0",
        },
    )
    try:
        with urlrequest.urlopen(req, timeout=timeout) as response:
            body = response.read().decode("utf-8")
    except urlerror.HTTPError as exc:
        detail = ""
        try:
            detail = collapse_whitespace(exc.read().decode("utf-8", errors="ignore"))[:180]
        except Exception:
            detail = ""
        suffix = f" ({detail})" if detail else ""
        raise RuntimeError(f"HTTP {exc.code}{suffix}") from exc
    except urlerror.URLError as exc:
        raise RuntimeError(f"Network error: {exc.reason}") from exc
    except TimeoutError as exc:
        raise RuntimeError("Request timed out") from exc
    try:
        payload = json.loads(body)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Invalid JSON response: {exc}") from exc
    if not isinstance(payload, dict):
        raise RuntimeError("Unexpected response payload")
    return payload


def normalize_uniprot_query(query: str) -> str:
    text = collapse_whitespace(query).rstrip("*")
    return text


def build_uniprot_query_variants(query: str) -> list[str]:
    normalized = normalize_uniprot_query(query)
    if not normalized:
        return []
    variants = [normalized]
    if "-" in normalized or "_" in normalized:
        collapsed_separators = re.sub(r"[-_]+", "", normalized)
        if collapsed_separators and collapsed_separators.upper() != normalized.upper():
            variants.append(collapsed_separators)
    return dedupe_keep_order(variants)


def extract_organism_name(entry: dict[str, Any]) -> str:
    organism = entry.get("organism")
    if isinstance(organism, dict):
        direct = collapse_whitespace(
            organism.get("scientificName")
            or organism.get("commonName")
            or organism.get("taxonName")
        )
        if direct:
            return direct
    return collapse_whitespace(entry.get("organism_name") or entry.get("organismName"))


def extract_length(entry: dict[str, Any]) -> Optional[int]:
    sequence = entry.get("sequence")
    if isinstance(sequence, dict):
        value = sequence.get("length")
        if isinstance(value, int):
            return value
        if isinstance(value, str) and value.strip().isdigit():
            return int(value.strip())
    value = entry.get("length")
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.strip().isdigit():
        return int(value.strip())
    return None


def extract_function_text(entry: dict[str, Any]) -> str:
    comments = entry.get("comments")
    if isinstance(comments, list):
        for comment in comments:
            if not isinstance(comment, dict):
                continue
            comment_type = collapse_whitespace(comment.get("commentType")).upper()
            if comment_type != "FUNCTION":
                continue
            texts = comment.get("texts")
            if isinstance(texts, list):
                for text_obj in texts:
                    if not isinstance(text_obj, dict):
                        continue
                    value = collapse_whitespace(text_obj.get("value"))
                    if value:
                        return value
    fallback = entry.get("cc_function")
    if isinstance(fallback, str):
        return fallback
    if isinstance(fallback, list):
        for row in fallback:
            if isinstance(row, str) and row.strip():
                return row
            if isinstance(row, dict):
                value = collapse_whitespace(row.get("value") or row.get("text"))
                if value:
                    return value
    return ""


def extract_gene_names(entry: dict[str, Any]) -> list[str]:
    candidates: list[str] = []
    genes = entry.get("genes")
    if isinstance(genes, list):
        for gene in genes:
            if not isinstance(gene, dict):
                continue
            gene_name = gene.get("geneName")
            if isinstance(gene_name, dict):
                value = collapse_whitespace(gene_name.get("value"))
                if value:
                    candidates.append(value)
            synonyms = gene.get("synonyms")
            if isinstance(synonyms, list):
                for synonym in synonyms:
                    if not isinstance(synonym, dict):
                        continue
                    value = collapse_whitespace(synonym.get("value"))
                    if value:
                        candidates.append(value)
    fallback = entry.get("gene_names") or entry.get("geneNames")
    if isinstance(fallback, str):
        candidates.extend(part.strip() for part in fallback.split() if part.strip())
    elif isinstance(fallback, list):
        for row in fallback:
            if isinstance(row, str):
                candidates.append(row)
            elif isinstance(row, dict):
                candidates.append(collapse_whitespace(row.get("value") or row.get("name")))
    return dedupe_keep_order(candidates)


def extract_cross_reference_ids(
    entry: dict[str, Any], database_name: str, fallback_key: Optional[str] = None
) -> list[str]:
    wanted = database_name.strip().upper()
    output: list[str] = []
    refs = entry.get("uniProtKBCrossReferences")
    if isinstance(refs, list):
        for ref in refs:
            if not isinstance(ref, dict):
                continue
            db = collapse_whitespace(ref.get("database")).upper()
            if db != wanted:
                continue
            ref_id = collapse_whitespace(ref.get("id"))
            if ref_id:
                output.append(ref_id)
    if not output and fallback_key:
        fallback = entry.get(fallback_key)
        if isinstance(fallback, list):
            for row in fallback:
                if isinstance(row, str):
                    output.append(row)
                    continue
                if not isinstance(row, dict):
                    continue
                ref_id = collapse_whitespace(row.get("id") or row.get("value") or row.get("name"))
                if ref_id:
                    output.append(ref_id)
        elif isinstance(fallback, str):
            output.extend(part.strip() for part in re.split(r"[;, ]+", fallback) if part.strip())
    return dedupe_keep_order(output)


def parse_uniprot_entry_summary(entry: dict[str, Any]) -> dict[str, Any]:
    accession = collapse_whitespace(entry.get("primaryAccession") or entry.get("accession"))
    entry_name = collapse_whitespace(
        entry.get("uniProtkbId") or entry.get("id") or entry.get("entryName")
    )
    gene_names = extract_gene_names(entry)
    organism = extract_organism_name(entry)
    length = extract_length(entry)
    function_snippet = first_sentence(extract_function_text(entry))
    pdb_ids = extract_cross_reference_ids(entry, "PDB", fallback_key="xref_pdb")
    alphafold_ids = extract_cross_reference_ids(
        entry, "AlphaFoldDB", fallback_key="xref_alphafolddb"
    )
    return {
        "accession": accession,
        "entryName": entry_name,
        "geneNames": gene_names,
        "geneName": gene_names[0] if gene_names else entry_name,
        "organism": organism,
        "length": length,
        "functionSnippet": function_snippet,
        "structures": {
            "pdbCount": len(pdb_ids),
            "hasAlphaFold": bool(alphafold_ids),
            "alphaFoldCount": len(alphafold_ids),
            "pdbIds": pdb_ids,
            "alphaFoldIds": alphafold_ids,
        },
    }


def extract_resolution_value(raw: Any) -> Optional[float]:
    if isinstance(raw, (int, float)):
        value = float(raw)
        return value if value > 0 else None
    if isinstance(raw, str):
        text = raw.strip()
        if not text:
            return None
        try:
            value = float(text)
        except ValueError:
            match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)
            if not match:
                return None
            try:
                value = float(match.group(0))
            except ValueError:
                return None
        return value if value > 0 else None
    if isinstance(raw, dict):
        for key in (
            "resolution",
            "resolution_combined",
            "ls_d_res_high",
            "d_resolution_high",
            "highest_res",
            "value",
        ):
            if key not in raw:
                continue
            value = extract_resolution_value(raw.get(key))
            if value is not None:
                return value
    if isinstance(raw, list):
        for item in raw:
            value = extract_resolution_value(item)
            if value is not None:
                return value
    return None


def extract_model_count(raw: Any) -> Optional[int]:
    if isinstance(raw, bool):
        return None
    if isinstance(raw, int):
        return raw if raw > 0 else None
    if isinstance(raw, float):
        value = int(raw)
        return value if value > 0 else None
    if isinstance(raw, str):
        text = raw.strip()
        if not text:
            return None
        if text.isdigit():
            value = int(text)
            return value if value > 0 else None
        match = re.search(r"\d+", text)
        if not match:
            return None
        value = int(match.group(0))
        return value if value > 0 else None
    if isinstance(raw, dict):
        for key in (
            "deposited_model_count",
            "number_of_models",
            "model_count",
            "conformers_submitted_total_number",
            "value",
        ):
            if key not in raw:
                continue
            value = extract_model_count(raw.get(key))
            if value is not None:
                return value
    if isinstance(raw, list):
        for item in raw:
            value = extract_model_count(item)
            if value is not None:
                return value
    return None


def default_pdb_summary(pdb_id: str) -> dict[str, Any]:
    normalized = collapse_whitespace(pdb_id).upper()
    return {
        "pdbId": normalized,
        "title": "—",
        "method": "—",
        "resolution": None,
        "modelCount": None,
        "organisms": [],
        "referenceAuthors": "",
        "referenceJournal": "",
        "referenceYear": "",
    }


def merge_pdb_summary(primary: dict[str, Any], fallback: dict[str, Any]) -> dict[str, Any]:
    merged_reference = merge_reference_fields(primary, fallback)
    merged = {
        "pdbId": collapse_whitespace(primary.get("pdbId") or fallback.get("pdbId")).upper(),
        "title": collapse_whitespace(primary.get("title")) or "—",
        "method": collapse_whitespace(primary.get("method")) or "—",
        "resolution": extract_resolution_value(primary.get("resolution")),
        "modelCount": extract_model_count(primary.get("modelCount")),
        "organisms": merge_organism_lists(primary.get("organisms"), fallback.get("organisms")),
        "referenceAuthors": merged_reference["referenceAuthors"],
        "referenceJournal": merged_reference["referenceJournal"],
        "referenceYear": merged_reference["referenceYear"],
    }
    fallback_title = collapse_whitespace(fallback.get("title")) or "—"
    fallback_method = collapse_whitespace(fallback.get("method")) or "—"
    fallback_resolution = extract_resolution_value(fallback.get("resolution"))
    fallback_model_count = extract_model_count(fallback.get("modelCount"))
    if merged["title"] == "—" and fallback_title != "—":
        merged["title"] = fallback_title
    if merged["method"] == "—" and fallback_method != "—":
        merged["method"] = fallback_method
    if merged["resolution"] is None and fallback_resolution is not None:
        merged["resolution"] = fallback_resolution
    if merged["modelCount"] is None and fallback_model_count is not None:
        merged["modelCount"] = fallback_model_count
    return merged


def parse_pdbe_summary_payload(pdb_id: str, payload: dict[str, Any]) -> dict[str, Any]:
    key_lower = pdb_id.lower()
    key_upper = pdb_id.upper()
    raw_entry = payload.get(key_lower) or payload.get(key_upper) or payload.get(pdb_id)
    entry: dict[str, Any] = {}
    if isinstance(raw_entry, list) and raw_entry:
        first_row = raw_entry[0]
        if isinstance(first_row, dict):
            entry = first_row
    elif isinstance(raw_entry, dict):
        entry = raw_entry
    title = collapse_whitespace(entry.get("title"))
    method_raw = entry.get("experimental_method")
    if isinstance(method_raw, list):
        methods = [collapse_whitespace(row) for row in method_raw if collapse_whitespace(row)]
        method = ", ".join(methods) if methods else "—"
    else:
        method = collapse_whitespace(method_raw) or "—"
    resolution_value = extract_resolution_value(entry.get("resolution"))
    if resolution_value is None:
        resolution_value = extract_resolution_value(entry.get("resolution_combined"))
    if resolution_value is None:
        resolution_value = extract_resolution_value(entry.get("resolution_high"))
    if resolution_value is None:
        resolution_value = extract_resolution_value(entry.get("highest_resolution"))
    model_count = extract_model_count(
        entry.get("number_of_models")
        or entry.get("model_count")
        or entry.get("deposited_model_count")
    )
    organisms = normalize_organism_list(
        {
            "source": entry.get("source"),
            "sources": entry.get("sources"),
            "organism": entry.get("organism"),
            "organism_name": entry.get("organism_name"),
            "organism_scientific_name": entry.get("organism_scientific_name"),
        }
    )
    citation_rows: list[dict[str, Any]] = []
    for key in ("citation", "citations"):
        value = entry.get(key)
        if isinstance(value, list):
            citation_rows.extend(row for row in value if isinstance(row, dict))
    if isinstance(entry.get("primary_citation"), dict):
        citation_rows.append(entry["primary_citation"])
    reference = normalize_reference_fields(
        {
            "authors": entry.get("entry_authors") or entry.get("entryAuthors") or entry.get("authors")
        }
    )
    for citation in citation_rows:
        parsed = normalize_reference_fields(
            {
                "authors": citation.get("author_list")
                or citation.get("authorList")
                or citation.get("authors")
                or citation.get("rcsb_authors"),
                "journal": citation.get("journal_abbrev")
                or citation.get("journalAbbrev")
                or citation.get("journal_full")
                or citation.get("journal_name")
                or citation.get("journal"),
                "year": citation.get("year")
                or citation.get("publication_year")
                or citation.get("journal_year"),
            }
        )
        if not reference["referenceAuthors"] and parsed["referenceAuthors"]:
            reference["referenceAuthors"] = parsed["referenceAuthors"]
        if not reference["referenceJournal"] and parsed["referenceJournal"]:
            reference["referenceJournal"] = parsed["referenceJournal"]
        if not reference["referenceYear"] and parsed["referenceYear"]:
            reference["referenceYear"] = parsed["referenceYear"]
    return {
        "pdbId": key_upper,
        "title": title or "—",
        "method": method,
        "resolution": resolution_value,
        "modelCount": model_count,
        "organisms": organisms,
        "referenceAuthors": reference["referenceAuthors"],
        "referenceJournal": reference["referenceJournal"],
        "referenceYear": reference["referenceYear"],
    }


def parse_pdbe_molecules_organisms_payload(pdb_id: str, payload: dict[str, Any]) -> list[str]:
    key_lower = pdb_id.lower()
    key_upper = pdb_id.upper()
    raw_entry = payload.get(key_lower) or payload.get(key_upper) or payload.get(pdb_id)
    return normalize_organism_list(raw_entry)


def parse_rcsb_entry_payload(pdb_id: str, payload: dict[str, Any]) -> dict[str, Any]:
    key_upper = collapse_whitespace(pdb_id).upper()
    struct_block = payload.get("struct")
    title = collapse_whitespace(struct_block.get("title")) if isinstance(struct_block, dict) else ""
    methods: list[str] = []
    exptl = payload.get("exptl")
    if isinstance(exptl, list):
        for row in exptl:
            if not isinstance(row, dict):
                continue
            method_value = collapse_whitespace(row.get("method"))
            if method_value:
                methods.append(method_value)
    if not methods:
        entry_info = payload.get("rcsb_entry_info")
        if isinstance(entry_info, dict):
            raw_method = entry_info.get("experimental_method")
            if isinstance(raw_method, list):
                methods.extend(
                    collapse_whitespace(item)
                    for item in raw_method
                    if collapse_whitespace(item)
                )
            else:
                method_value = collapse_whitespace(raw_method)
                if method_value:
                    methods.append(method_value)
    method = ", ".join(dedupe_keep_order(methods)) if methods else "—"
    entry_info = payload.get("rcsb_entry_info")
    resolution_value = extract_resolution_value(
        entry_info.get("resolution_combined") if isinstance(entry_info, dict) else None
    )
    if resolution_value is None:
        resolution_value = extract_resolution_value(payload.get("refine"))
    if resolution_value is None:
        resolution_value = extract_resolution_value(payload.get("em_3d_reconstruction"))
    model_count = extract_model_count(
        entry_info.get("deposited_model_count") if isinstance(entry_info, dict) else None
    )
    if model_count is None:
        model_count = extract_model_count(payload.get("pdbx_nmr_ensemble"))
    organisms = normalize_organism_list(
        {
            "rcsb_entity_source_organism": payload.get("rcsb_entity_source_organism"),
            "entity_src_nat": payload.get("entity_src_nat"),
            "entity_src_gen": payload.get("entity_src_gen"),
            "pdbx_entity_src_syn": payload.get("pdbx_entity_src_syn"),
            "polymer_entities": payload.get("polymer_entities"),
        }
    )
    citation_rows: list[dict[str, Any]] = []
    if isinstance(payload.get("rcsb_primary_citation"), dict):
        citation_rows.append(payload["rcsb_primary_citation"])
    citation_value = payload.get("citation")
    if isinstance(citation_value, list):
        citation_rows.extend(row for row in citation_value if isinstance(row, dict))
    elif isinstance(citation_value, dict):
        citation_rows.append(citation_value)
    reference = {
        "referenceAuthors": "",
        "referenceJournal": "",
        "referenceYear": "",
    }
    for citation in citation_rows:
        parsed = normalize_reference_fields(
            {
                "authors": citation.get("rcsb_authors")
                or citation.get("authors")
                or citation.get("author_list"),
                "journal": citation.get("rcsb_journal_abbrev")
                or citation.get("journal_abbrev")
                or citation.get("journal_name")
                or citation.get("journal_full")
                or citation.get("journal"),
                "year": citation.get("year")
                or citation.get("publication_year")
                or citation.get("journal_year"),
            }
        )
        if not reference["referenceAuthors"] and parsed["referenceAuthors"]:
            reference["referenceAuthors"] = parsed["referenceAuthors"]
        if not reference["referenceJournal"] and parsed["referenceJournal"]:
            reference["referenceJournal"] = parsed["referenceJournal"]
        if not reference["referenceYear"] and parsed["referenceYear"]:
            reference["referenceYear"] = parsed["referenceYear"]
    return {
        "pdbId": key_upper,
        "title": title or "—",
        "method": method,
        "resolution": resolution_value,
        "modelCount": model_count,
        "organisms": organisms,
        "referenceAuthors": reference["referenceAuthors"],
        "referenceJournal": reference["referenceJournal"],
        "referenceYear": reference["referenceYear"],
    }


def fetch_pdbe_molecule_organisms(pdb_id: str) -> list[str]:
    normalized = collapse_whitespace(pdb_id).upper()
    if not normalized:
        return []
    url = PD_BE_MOLECULES_URL_TEMPLATE.format(pdb_id=normalized.lower())
    try:
        payload = fetch_json(url)
    except Exception:
        return []
    if not isinstance(payload, dict):
        return []
    return parse_pdbe_molecules_organisms_payload(normalized, payload)


def fetch_pdbe_summary(pdb_id: str) -> dict[str, Any]:
    normalized = collapse_whitespace(pdb_id).upper()
    if not normalized:
        return default_pdb_summary("")
    now = time.time()
    with PDB_SUMMARY_CACHE_LOCK:
        cached = PDB_SUMMARY_CACHE.get(normalized)
        if cached and (now - cached[0]) < PDB_SUMMARY_TTL_SECONDS:
            cached_summary = cached[1]
            cached_method = collapse_whitespace(cached_summary.get("method")).upper()
            cached_resolution = extract_resolution_value(cached_summary.get("resolution"))
            cached_title = collapse_whitespace(cached_summary.get("title"))
            cached_organisms = normalize_organism_list(cached_summary.get("organisms"))
            cache_is_missing_critical_fields = (
                not cached_title
                or cached_method in ("", "—")
                or not cached_organisms
                or not has_complete_reference_fields(cached_summary)
                or (
                    cached_resolution is None
                    and (
                        "X-RAY" in cached_method
                        or "ELECTRON" in cached_method
                        or "CRYO-EM" in cached_method
                    )
                )
            )
            if not cache_is_missing_critical_fields:
                return cached_summary
    summary = default_pdb_summary(normalized)
    url = PD_BE_SUMMARY_URL_TEMPLATE.format(pdb_id=normalized.lower())
    try:
        payload = fetch_json(url)
        summary = parse_pdbe_summary_payload(normalized, payload)
    except Exception:
        pass
    if not normalize_organism_list(summary.get("organisms")):
        molecule_organisms = fetch_pdbe_molecule_organisms(normalized)
        if molecule_organisms:
            summary["organisms"] = merge_organism_lists(summary.get("organisms"), molecule_organisms)
    # PDBe summary frequently omits resolution/method/citation details. Fill from RCSB.
    if (
        summary["resolution"] is None
        or summary["method"] == "—"
        or summary["title"] == "—"
        or not normalize_organism_list(summary.get("organisms"))
        or not has_complete_reference_fields(summary)
    ):
        rcsb_url = RCSB_ENTRY_URL_TEMPLATE.format(pdb_id=normalized.upper())
        try:
            rcsb_payload = fetch_json(rcsb_url)
            rcsb_summary = parse_rcsb_entry_payload(normalized, rcsb_payload)
            summary = merge_pdb_summary(summary, rcsb_summary)
        except Exception:
            pass
    with PDB_SUMMARY_CACHE_LOCK:
        PDB_SUMMARY_CACHE[normalized] = (now, summary)
    return summary


def fetch_pdb_summaries_parallel(pdb_ids: list[str]) -> list[dict[str, Any]]:
    ordered_ids = dedupe_keep_order([collapse_whitespace(pdb_id).upper() for pdb_id in pdb_ids])
    if not ordered_ids:
        return []
    workers = max(1, min(PDB_SUMMARY_MAX_WORKERS, len(ordered_ids)))
    results: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        future_map = {pool.submit(fetch_pdbe_summary, pdb_id): pdb_id for pdb_id in ordered_ids}
        for future in as_completed(future_map):
            pdb_id = future_map[future]
            try:
                results[pdb_id] = future.result()
            except Exception:
                results[pdb_id] = default_pdb_summary(pdb_id)
    return [
        results.get(pdb_id, default_pdb_summary(pdb_id))
        for pdb_id in ordered_ids
    ]


def annotate_pdb_loadability(summary: dict[str, Any]) -> dict[str, Any]:
    model_count = extract_model_count(summary.get("modelCount"))
    loadable = True
    reason = ""
    # Viewer cannot reliably load multi-model ensembles.
    if model_count is not None and model_count > 1:
        loadable = False
        reason = f"Cannot load (ensemble; {model_count} models)"
    out = dict(summary)
    out["modelCount"] = model_count
    out["isEnsemble"] = bool(model_count is not None and model_count > 1)
    out["loadable"] = loadable
    if reason:
        out["loadabilityReason"] = reason
    return out


def fetch_uniprot_search_variant(query: str, reviewed: bool = False) -> dict[str, Any]:
    normalized = normalize_uniprot_query(query)
    if not normalized:
        return {"results": []}
    built_query = f"{normalized}*"
    if reviewed:
        built_query = f"({built_query}) AND reviewed:true"
    params = {
        "query": built_query,
        "format": "json",
        "size": str(UNIPROT_SEARCH_FETCH_SIZE),
        "fields": UNIPROT_SEARCH_FIELDS,
    }
    url = f"{UNIPROT_SEARCH_URL}?{urlparse.urlencode(params)}"
    return fetch_json(url)


def fetch_uniprot_search(query: str, reviewed: bool = False) -> dict[str, Any]:
    variants = build_uniprot_query_variants(query)
    if not variants:
        return {"results": []}
    merged_results: list[dict[str, Any]] = []
    seen_accessions: set[str] = set()
    had_success = False
    last_error: Optional[Exception] = None
    for variant in variants:
        try:
            payload = fetch_uniprot_search_variant(variant, reviewed=reviewed)
        except Exception as exc:
            last_error = exc
            continue
        had_success = True
        raw_results = payload.get("results")
        if not isinstance(raw_results, list):
            continue
        for row in raw_results:
            if not isinstance(row, dict):
                continue
            accession = collapse_whitespace(row.get("primaryAccession") or row.get("accession")).upper()
            if accession and accession in seen_accessions:
                continue
            if accession:
                seen_accessions.add(accession)
            merged_results.append(row)
    if had_success:
        return {"results": merged_results}
    if last_error is not None:
        raise last_error
    return {"results": []}


def fetch_uniprot_entry(accession: str) -> dict[str, Any]:
    normalized = collapse_whitespace(accession).upper()
    if not normalized:
        raise RuntimeError("No UniProt accession provided")
    url = UNIPROT_ENTRY_URL_TEMPLATE.format(accession=urlparse.quote(normalized, safe=""))
    return fetch_json(url)


@app.get("/protein-search")
def protein_search(query: str, reviewed: bool = False):
    normalized = normalize_uniprot_query(query)
    if not normalized:
        return {"query": "", "items": [], "resultCount": 0, "fetchSize": UNIPROT_SEARCH_FETCH_SIZE}
    try:
        payload = fetch_uniprot_search(normalized, reviewed=reviewed)
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"UniProt search failed: {exc}") from exc
    raw_results = payload.get("results")
    items: list[dict[str, Any]] = []
    if isinstance(raw_results, list):
        for row in raw_results:
            if not isinstance(row, dict):
                continue
            parsed = parse_uniprot_entry_summary(row)
            if parsed.get("accession"):
                items.append(parsed)
    return {
        "query": normalized,
        "items": items,
        "resultCount": len(items),
        "fetchSize": UNIPROT_SEARCH_FETCH_SIZE,
    }


@app.get("/protein-structures/{accession}")
def protein_structures(accession: str):
    normalized = collapse_whitespace(accession).upper()
    if not normalized:
        raise HTTPException(status_code=400, detail="accession is required")
    try:
        entry = fetch_uniprot_entry(normalized)
    except Exception as exc:
        message = str(exc)
        if "HTTP 404" in message:
            raise HTTPException(status_code=404, detail=f"UniProt entry not found: {normalized}") from exc
        raise HTTPException(status_code=502, detail=f"UniProt entry fetch failed: {exc}") from exc

    summary = parse_uniprot_entry_summary(entry)
    pdb_ids = summary["structures"]["pdbIds"]
    alphafold_ids = summary["structures"]["alphaFoldIds"]
    pdb_entries_all = fetch_pdb_summaries_parallel(pdb_ids)
    pdb_entries = [annotate_pdb_loadability(entry) for entry in pdb_entries_all]
    alphafold_entries = [
        {
            "alphafoldId": af_id,
            **ALPHAFOLD_REFERENCE,
        }
        for af_id in dedupe_keep_order(alphafold_ids)
    ]
    return {
        "accession": summary.get("accession") or normalized,
        "entryName": summary.get("entryName") or "",
        "geneNames": summary.get("geneNames") or [],
        "geneName": summary.get("geneName") or "",
        "organism": summary.get("organism") or "",
        "length": summary.get("length"),
        "functionSnippet": summary.get("functionSnippet") or "",
        "pdb": pdb_entries,
        "alphafold": alphafold_entries,
    }


def resolve_mmcif(pdb_id: Optional[str], mmcif_text: Optional[str]) -> Optional[str]:
    if mmcif_text:
        return mmcif_text
    if not pdb_id:
        return None
    normalized = collapse_whitespace(pdb_id).lower()
    if not normalized:
        return None
    cache_key_text = f"mmcif:{normalized}"
    if not CHAPI_LOW_MEMORY_MODE:
        cached = STRUCTURE_TEXT_CACHE.get(cache_key_text)
        if isinstance(cached, str) and cached:
            return cached
    text = fetch_mmcif(normalized)
    if text and not CHAPI_LOW_MEMORY_MODE:
        STRUCTURE_TEXT_CACHE.set(cache_key_text, text)
    return text


def write_temp_structure(text: str, suffix: str) -> Path:
    temp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    temp.write(text.encode("utf-8"))
    temp.flush()
    temp.close()
    return Path(temp.name)


def build_chapi_mesh_cache_key(payload: dict, source_key: Optional[str] = None) -> str:
    raw_chain_ids = payload.get("chainIds")
    normalized_chain_ids: list[str] = []
    if isinstance(raw_chain_ids, list):
        seen_chain_ids: set[str] = set()
        for chain in raw_chain_ids:
            token = str(chain or "").strip()
            if not token or token in seen_chain_ids:
                continue
            seen_chain_ids.add(token)
            normalized_chain_ids.append(token)
        normalized_chain_ids.sort()
    rep_options = {
        'outputFormat': payload.get('outputFormat', 'json'),
        "format": payload.get("format"),
        "representation": payload.get("representation"),
        "mode": payload.get("mode"),
        "againstDarkBackground": bool(payload.get("againstDarkBackground", False)),
        "bondWidth": payload.get("bondWidth"),
        "atomRadiusToBondWidthRatio": payload.get("atomRadiusToBondWidthRatio"),
        "smoothnessFactor": payload.get("smoothnessFactor"),
        "nonDrawCids": payload.get("nonDrawCids"),
        "carbonColor": payload.get("carbonColor"),
        "cid": payload.get("cid"),
        "colourScheme": payload.get("colourScheme"),
        "style": payload.get("style"),
        "secondaryStructureUsage": payload.get("secondaryStructureUsage"),
        "splitByChain": bool(payload.get("splitByChain", False)),
        "chainIds": normalized_chain_ids,
    }
    options_json = json.dumps(rep_options, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    if source_key:
        structure_part = source_key
    else:
        text = payload.get("text") or ""
        fmt = str(payload.get("format") or "")
        text_digest = hashlib.blake2b(text.encode("utf-8"), digest_size=20).hexdigest()
        structure_part = f"{fmt}:{len(text)}:{text_digest}"
    digest_src = f"{structure_part}|{options_json}"
    return hashlib.blake2b(digest_src.encode("utf-8"), digest_size=20).hexdigest()


def chapi_cid_selects_single_chain(cid: Any) -> bool:
    token = str(cid or "").strip()
    if not token.startswith("//") or token in {"//", "///"}:
        return False
    chain_token = token[2:].split("/", 1)[0].strip()
    return bool(chain_token and not any(sep in chain_token for sep in (",", ";", "|", " ")))


def chapi_single_chain_from_cid(cid: Any) -> Optional[str]:
    if not chapi_cid_selects_single_chain(cid):
        return None
    return str(cid or "").strip()[2:].split("/", 1)[0].strip() or None


def pdb_line_chain_id(line: str) -> str:
    return line[21].strip() if len(line) > 21 else ""


def filter_pdb_text_to_single_chain(text: str, chain_id: str) -> tuple[str, int]:
    wanted = str(chain_id or "").strip()
    if not text or not wanted:
        return "", 0
    kept: list[str] = []
    atom_rows = 0
    for line in text.splitlines():
        record = line[:6].strip().upper()
        if record in {"ATOM", "HETATM", "ANISOU", "TER"}:
            if pdb_line_chain_id(line) != wanted:
                continue
            if record in {"ATOM", "HETATM"}:
                atom_rows += 1
        elif record in {"CONECT", "MASTER"}:
            continue
        kept.append(line)
    if atom_rows <= 0:
        return "", 0
    if not kept or kept[-1].strip().upper() != "END":
        kept.append("END")
    return "\n".join(kept) + "\n", atom_rows


def filter_mmcif_text_to_single_chain(text: str, chain_id: str) -> tuple[str, int]:
    """Filter actual CIF rows, independently of physical line wrapping/packing.

    Gemmi preserves all other categories, quoting and missing-value tokens. The
    same conformer/model selection is subsequently applied by the structure reader.
    """
    wanted = str(chain_id or "").strip()
    if not text or not wanted:
        return "", 0
    import gemmi
    document = gemmi.cif.read_string(text)
    matched_atom_rows = 0
    for block in document:
        table = block.find_mmcif_category("_atom_site.")
        if not table:
            continue
        tags = [str(tag).lower() for tag in table.tags]
        indices = [tags.index(tag) for tag in ("_atom_site.auth_asym_id", "_atom_site.label_asym_id") if tag in tags]
        if not indices:
            return "", 0
        columns: list[list[str]] = [[] for _ in tags]
        for row in table:
            if wanted in {gemmi.cif.as_string(row[column]) for column in indices}:
                matched_atom_rows += 1
                for column, value in zip(columns, row):
                    column.append(value)
        # Replace in one pass: repeated row erasure shifts the underlying array
        # and becomes quadratic for large assemblies.
        table.loop.set_all_values(columns)
    if not matched_atom_rows:
        return "", 0
    return document.as_string(), matched_atom_rows


def filter_structure_text_to_single_chain(text: str, fmt: str, chain_id: str) -> tuple[str, int]:
    if fmt == "pdb":
        return filter_pdb_text_to_single_chain(text, chain_id)
    if fmt == "mmcif":
        return filter_mmcif_text_to_single_chain(text, chain_id)
    return "", 0


def _build_chapi_bridge_env() -> dict:
    env = os.environ.copy()
    if CHAPI_PREFIX:
        env.setdefault("COOT_PREFIX", CHAPI_PREFIX)
    # Keep a locally rebuilt Coot library scoped to the mesh worker. The API
    # process and other tools using the same conda environment keep their loader
    # settings. An explicitly empty CHAPI_LIBRARY_DIR disables the local overlay.
    native_library = {
        "darwin": "libcootmoleculestotriangles.1.1.dylib",
        "linux": "libcootmoleculestotriangles.so.1.1",
    }.get(sys.platform)
    library_dir = env.get("CHAPI_LIBRARY_DIR")
    if library_dir is None and native_library:
        local_dir = CHAPI_BRIDGE.parent.parent / ".coot-ribbon" / "lib"
        if (local_dir / native_library).is_file():
            library_dir = str(local_dir)
    if library_dir:
        resolved_dir = Path(library_dir).expanduser().resolve()
        if not native_library or not (resolved_dir / native_library).is_file():
            raise RuntimeError("CHAPI_LIBRARY_DIR must contain the matching Coot 1.1.20 mesh library.")
        loader_key = "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
        inherited = env.get(loader_key, "")
        env[loader_key] = str(resolved_dir) + (os.pathsep + inherited if inherited else "")
    env.setdefault("PYTHONUNBUFFERED", "1")
    return env


class ChapiWorkerTransportError(RuntimeError):
    def __init__(self, message: str, *, retryable: bool = False):
        super().__init__(message)
        self.retryable = retryable


class ChapiNativeProcessError(RuntimeError):
    def __init__(self, returncode: int, detail: str = ""):
        self.returncode = returncode
        message = f"chapi bridge failed (exit code {returncode})"
        super().__init__(f"{message}: {detail}" if detail else message)


def _chapi_worker_native_exit(proc: subprocess.Popen) -> Optional[ChapiNativeProcessError]:
    returncode = proc.poll()
    if returncode is None:
        try:
            returncode = proc.wait(timeout=1.0)
        except subprocess.TimeoutExpired:
            return None
    if returncode:
        return ChapiNativeProcessError(returncode)
    return None


def classify_chapi_mesh_error(exc: Exception) -> tuple[str, str]:
    message = str(exc or "").strip() or "CHAPI mesh generation failed."
    normalized = message.lower()

    if "pdbid, pdbtext, or mmciftext is required" in normalized:
        return "CHAPI-REQ-001", message
    if "chapi_bridge.py not found" in normalized:
        return "CHAPI-CONFIG-001", message
    if "failed to fetch mmcif" in normalized:
        return "CHAPI-UPSTREAM-001", message
    if any(
        token in normalized
        for token in (
            "out of memory",
            "memoryerror",
            "exit code -9",
            "signal 9",
            "killed",
        )
    ):
        return "CHAPI-OOM-001", message
    if "timed out" in normalized:
        return "CHAPI-WORKER-001", message
    if any(
        token in normalized
        for token in (
            "stdout is unavailable",
            "stdin is unavailable",
            "closed output unexpectedly",
            "communication failed",
            "worker response",
            "unknown chapi worker status",
        )
    ):
        return "CHAPI-WORKER-002", message
    if any(
        token in normalized
        for token in (
            "invalid json",
            "empty json",
            "bridge failed",
            "failed to read structure in chapi bridge",
        )
    ):
        return "CHAPI-BRIDGE-001", message
    return "CHAPI-MESH-001", message


def _stop_chapi_worker_locked() -> None:
    global CHAPI_WORKER_PROC
    proc = CHAPI_WORKER_PROC
    CHAPI_WORKER_PROC = None
    if not proc:
        return
    try:
        if proc.stdin:
            try:
                proc.stdin.close()
            except Exception:
                pass
        proc.terminate()
        proc.wait(timeout=2.0)
    except Exception:
        try:
            proc.kill()
            proc.wait(timeout=1.0)
        except Exception:
            pass
    finally:
        for stream in (proc.stdout, proc.stderr):
            if stream:
                stream.close()


def _start_chapi_worker_locked() -> subprocess.Popen:
    global CHAPI_WORKER_PROC
    if CHAPI_WORKER_PROC and CHAPI_WORKER_PROC.poll() is None:
        return CHAPI_WORKER_PROC
    _stop_chapi_worker_locked()
    if not CHAPI_BRIDGE.exists():
        raise RuntimeError("chapi_bridge.py not found in api directory.")
    proc = subprocess.Popen(
        [CHAPI_PYTHON, str(CHAPI_BRIDGE), "--server"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        env=_build_chapi_bridge_env(),
    )
    CHAPI_WORKER_PROC = proc
    return proc


def _read_chapi_worker_line(proc: subprocess.Popen, timeout_seconds: int) -> bytes:
    stdout = proc.stdout
    if stdout is None:
        raise ChapiWorkerTransportError("CHAPI worker stdout is unavailable.")
    ready, _, _ = select.select([stdout], [], [], float(timeout_seconds))
    if not ready:
        raise ChapiWorkerTransportError(f"CHAPI worker timed out after {timeout_seconds}s.")
    line = stdout.readline()
    if not line:
        native_error = _chapi_worker_native_exit(proc)
        if native_error:
            raise native_error
        raise ChapiWorkerTransportError("CHAPI worker closed output unexpectedly.")
    return line


def _run_chapi_mesh_worker(payload: dict) -> bytes:
    request_line = (
        json.dumps({"payload": payload}, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
        + b"\n"
    )
    with CHAPI_WORKER_LOCK:
        try:
            proc = _start_chapi_worker_locked()
        except OSError as exc:
            raise ChapiWorkerTransportError(str(exc), retryable=True) from exc
        stdin = proc.stdin
        if stdin is None:
            _stop_chapi_worker_locked()
            raise ChapiWorkerTransportError("CHAPI worker stdin is unavailable.", retryable=True)
        try:
            stdin.write(request_line)
            stdin.flush()
            line = _read_chapi_worker_line(proc, CHAPI_WORKER_TIMEOUT_SECONDS)
            if line.startswith(b'BINARY\t'):
                try:
                    length = int(line.split(b'\t', 1)[1].strip())
                except ValueError as exc:
                    raise ChapiWorkerTransportError('Invalid binary worker response length') from exc
                if length < 12 or length > 512 * 1024**2:
                    raise ChapiWorkerTransportError('Invalid binary worker response length')
                body = _read_chapi_worker_bytes(proc, length, CHAPI_WORKER_TIMEOUT_SECONDS)
                if not body.startswith(b'ROAMIM01'):
                    raise ChapiWorkerTransportError('Invalid binary worker response')
                return body
        except (ChapiWorkerTransportError, ChapiNativeProcessError):
            _stop_chapi_worker_locked()
            raise
        except Exception as exc:
            native_error = _chapi_worker_native_exit(proc)
            _stop_chapi_worker_locked()
            if native_error:
                raise native_error from exc
            raise ChapiWorkerTransportError(str(exc) or "CHAPI worker communication failed.") from exc

    line = line.rstrip(b"\r\n")
    status, sep, body = line.partition(b"\t")
    if not sep:
        preview = line[:200].decode("utf-8", errors="replace")
        raise ChapiWorkerTransportError(f"Invalid CHAPI worker response: {preview}")
    if status == b"OK":
        stripped = body.lstrip()
        if not stripped.startswith((b"{", b"[")):
            preview = body[:200].decode("utf-8", errors="replace")
            raise ChapiWorkerTransportError(f"Invalid JSON from CHAPI worker: {preview}")
        return body
    if status == b"ERR":
        message = ""
        try:
            parsed = json.loads(body.decode("utf-8", errors="replace"))
            if isinstance(parsed, dict):
                message = str(parsed.get("error") or "").strip()
        except Exception:
            message = ""
        if not message:
            message = body.decode("utf-8", errors="replace").strip() or "CHAPI worker error"
        raise RuntimeError(message)
    preview = line[:200].decode("utf-8", errors="replace")
    raise ChapiWorkerTransportError(f"Unknown CHAPI worker status: {preview}")


def _read_chapi_worker_bytes(proc, length, timeout_seconds):
    """Read a framed body with a deadline, including bytes buffered by readline."""
    stdout = proc.stdout
    descriptor = stdout.fileno()
    deadline = time.monotonic() + timeout_seconds
    output = bytearray()
    pipe_ready = False
    os.set_blocking(descriptor, False)
    try:
        while len(output) < length:
            chunk = stdout.read1(min(65536, length - len(output)))
            if chunk:
                output.extend(chunk)
                pipe_ready = False
                continue
            if pipe_ready or proc.poll() is not None:
                native_error = _chapi_worker_native_exit(proc)
                if native_error:
                    raise native_error
                raise ChapiWorkerTransportError('CHAPI worker closed binary output unexpectedly')
            remaining = deadline - time.monotonic()
            if remaining <= 0 or not select.select([stdout], [], [], remaining)[0]:
                raise ChapiWorkerTransportError(f'CHAPI worker timed out after {timeout_seconds}s.')
            pipe_ready = True
        return bytes(output)
    finally:
        os.set_blocking(descriptor, True)


def _run_chapi_mesh_subprocess(payload: dict) -> bytes:
    if not CHAPI_BRIDGE.exists():
        raise RuntimeError("chapi_bridge.py not found in api directory.")
    payload_bytes = json.dumps(payload, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    proc = subprocess.run(
        [CHAPI_PYTHON, str(CHAPI_BRIDGE)],
        input=payload_bytes,
        capture_output=True,
        env=_build_chapi_bridge_env(),
        timeout=CHAPI_WORKER_TIMEOUT_SECONDS,
    )
    if proc.returncode != 0:
        stderr_text = proc.stderr.decode("utf-8", errors="replace").strip()
        stdout_text = proc.stdout.decode("utf-8", errors="replace").strip()
        detail = stderr_text or stdout_text
        raise ChapiNativeProcessError(proc.returncode, detail)
    raw = proc.stdout or b""
    if not raw:
        raise RuntimeError("Empty JSON from chapi bridge.")
    stripped = raw.lstrip()
    if not stripped.startswith((b"{", b"[", b'ROAMIM01')):
        preview = raw[:160].decode("utf-8", errors="replace")
        raise RuntimeError(f"Invalid JSON from chapi bridge: {preview}")
    return raw


def _dispatch_chapi_mesh(payload: dict) -> bytes:
    if CHAPI_PERSISTENT_WORKER:
        try:
            return _run_chapi_mesh_worker(payload)
        except ChapiWorkerTransportError as exc:
            # Only retry failures before a request could reach the worker.
            # Timeouts, native exits and malformed responses are not safe to
            # repeat with the same reader.
            if not exc.retryable:
                raise
    return _run_chapi_mesh_subprocess(payload)


def _chapi_reader_source_key(payload: dict) -> str:
    text = str(payload.get("text") or "")
    digest = hashlib.blake2b(text.encode("utf-8"), digest_size=20).hexdigest()
    return f"{payload.get('format') or ''}:{digest}"


def run_chapi_mesh(payload: dict) -> bytes:
    with CHAPI_MESH_EXECUTION_LOCK:
        source_key = payload.get("_reader_source_key") or _chapi_reader_source_key(payload)
        use_mmdb = payload.get("_use_mmdb_reader") is True or bool(CHAPI_MMDB_READER_CACHE.get(source_key))
        attempt = dict(payload, _use_mmdb_reader=True) if use_mmdb else payload
        try:
            return _dispatch_chapi_mesh(attempt)
        except ChapiNativeProcessError as exc:
            # Coot's Gemmi-to-MMDB sequence transfer can segfault on valid
            # structures (e.g. 6VXX). Retry once with its alternative reader;
            # never treat OOM, a timeout or an ordinary parse error as this bug.
            if exc.returncode != -signal.SIGSEGV or use_mmdb:
                raise
            output = _dispatch_chapi_mesh(dict(payload, _use_mmdb_reader=True))
            CHAPI_MMDB_READER_CACHE.set(source_key, True)
            return output


def maybe_cache_chapi_mesh(cache_key: str, output: bytes) -> None:
    if CHAPI_MESH_CACHE_MAX_BYTES <= 0:
        return
    if len(output) <= CHAPI_MESH_CACHE_MAX_BYTES:
        CHAPI_MESH_CACHE.set(cache_key, output)


def run_chapi_mesh_cached(payload: dict, cache_key: str) -> bytes:
    if CHAPI_MESH_CACHE_MAX_BYTES > 0:
        cached = CHAPI_MESH_CACHE.get(cache_key)
        if isinstance(cached, (bytes, bytearray)):
            return bytes(cached)

    leader = False
    wait_event: Optional[threading.Event] = None
    with CHAPI_MESH_INFLIGHT_LOCK:
        if CHAPI_MESH_CACHE_MAX_BYTES > 0:
            cached = CHAPI_MESH_CACHE.get(cache_key)
            if isinstance(cached, (bytes, bytearray)):
                return bytes(cached)
        wait_event = CHAPI_MESH_INFLIGHT.get(cache_key)
        if wait_event is None:
            wait_event = threading.Event()
            CHAPI_MESH_INFLIGHT[cache_key] = wait_event
            leader = True

    if leader:
        try:
            output = run_chapi_mesh(payload)
            maybe_cache_chapi_mesh(cache_key, output)
            return output
        finally:
            with CHAPI_MESH_INFLIGHT_LOCK:
                done_event = CHAPI_MESH_INFLIGHT.pop(cache_key, None)
                if done_event:
                    done_event.set()

    if wait_event:
        wait_event.wait(timeout=125.0)
    if CHAPI_MESH_CACHE_MAX_BYTES > 0:
        cached = CHAPI_MESH_CACHE.get(cache_key)
        if isinstance(cached, (bytes, bytearray)):
            return bytes(cached)

    output = run_chapi_mesh(payload)
    maybe_cache_chapi_mesh(cache_key, output)
    return output


def stub_report(pdb_id: str, chain_a: str, chain_b: str) -> dict:
    contact_template = {
        "residueA": {"chain": chain_a, "resName": "LEU", "seq": "42", "atom": "CD1"},
        "residueB": {"chain": chain_b, "resName": "VAL", "seq": "10", "atom": "CG1"},
        "distance": 3.8,
        "type": "hydrophobic",
        "category": "hydrophobic",
        "atomKeyA": f"{chain_a}:42:LEU:CD1",
        "atomKeyB": f"{chain_b}:10:VAL:CG1",
        "pairKey": (
            f"{chain_a}:42:LEU:CD1|{chain_b}:10:VAL:CG1"
            if f"{chain_a}:42:LEU:CD1" <= f"{chain_b}:10:VAL:CG1"
            else f"{chain_b}:10:VAL:CG1|{chain_a}:42:LEU:CD1"
        ),
        "arpeggio": {
            "type": "atom-atom",
            "terms": ["HYDROPHOBIC", "VDW"],
            "distance": 3.8,
            "interactingEntities": "INTER",
        },
        "asserted": {
            "family": "hydrophobic",
            "confidence": "medium",
            "evidence": ["nonpolar_pair", "distance_ok"],
        },
    }
    contacts = {
        "hydrogen_bonds": [],
        "polar_contacts": [],
        "base_pairing": [],
        "salt_bridges": [],
        "hydrophobic": [contact_template],
        "metal_coordination": [],
        "pi_pi": [],
        "pi_cation": [],
        "aromatic_packing": [],
        "other": [],
    }
    per_residue = {
        f"{chain_a}:42": {
            "chain": chain_a,
            "resName": "LEU",
            "seq": "42",
            "hydrophobic": 1,
            "hbond": 0,
            "polar_contact": 0,
            "base_pairing": 0,
            "salt_bridge": 0,
            "metal_coordination": 0,
            "pi_pi": 0,
            "pi_cation": 0,
            "aromatic_packing": 0,
            "vdw": 0,
            "clash": 0,
            "other": 0,
            "total": 1,
        },
        f"{chain_b}:10": {
            "chain": chain_b,
            "resName": "VAL",
            "seq": "10",
            "hydrophobic": 1,
            "hbond": 0,
            "polar_contact": 0,
            "base_pairing": 0,
            "salt_bridge": 0,
            "metal_coordination": 0,
            "pi_pi": 0,
            "pi_cation": 0,
            "aromatic_packing": 0,
            "vdw": 0,
            "clash": 0,
            "other": 0,
            "total": 1,
        },
    }
    return {
        "pdbId": pdb_id,
        "chainA": chain_a,
        "chainB": chain_b,
        "analysisVersion": "stub",
        "contacts": contacts,
        "perResidue": per_residue,
        "interfaceArea": None,
        "buriedFraction": None,
        "contactingResidueFraction": {chain_a: 0.5, chain_b: 0.5},
        "approxDeltaG": None,
        "meta": {
            "engine": "stub",
            "analysisVersion": "stub",
            "classifier": "plausibility+assertion:v1",
            "note": "Stub report returned because structure could not be fetched.",
        },
    }


@app.get("/chains")
def get_chains(pdbId: str):
    """Return chain identifiers for a given PDB entry."""
    pdb_lower = pdbId.strip().lower()
    try:
        mmcif_text = resolve_mmcif(pdb_lower, None)
        if not mmcif_text:
            raise RuntimeError("No structure data")
        chains, aliases = list_chains(mmcif_text)
        return {"chains": chains, "aliases": aliases.label_to_auth}
    except Exception:
        demo = DEMO_CHAINS.get(pdb_lower)
        if demo:
            return {"chains": demo, "aliases": {}}
        raise HTTPException(status_code=404, detail="Unable to load structure")


@app.post("/chains")
def post_chains(request: ChainsRequest):
    if not request.mmcifText:
        raise HTTPException(status_code=400, detail="mmcifText is required")
    chains, aliases = list_chains(request.mmcifText)
    return {"chains": chains, "aliases": aliases.label_to_auth}


def _resolve_analysis_source(request: AnalyzeRequest):
    pdb_id = (request.pdbId or "").strip().lower() or None
    if request.pdbText:
        text, fmt = request.pdbText, "pdb"
    elif request.mmcifText:
        text, fmt = request.mmcifText, "mmcif"
    else:
        try:
            text, fmt = resolve_mmcif(pdb_id, None), "mmcif"
        except Exception:
            text, fmt = None, "mmcif"
    key = cache_key(
        pdb_id, f"{fmt}\n{text or ''}", request.chainA, request.chainB,
        request.mode or "all", focus_residue=(request.focusResidue or "").strip() or None,
    )
    return pdb_id, text, fmt, key


def _accepts_gzip(request: Optional[Request]) -> bool:
    header = str(getattr(request, "headers", {}).get("accept-encoding", ""))
    qualities = {}
    for item in header.lower().split(','):
        parts = [part.strip() for part in item.split(';')]
        if not parts[0]:
            continue
        quality = 1.0
        for parameter in parts[1:]:
            if parameter.startswith('q='):
                try:
                    quality = float(parameter[2:])
                except ValueError:
                    quality = 0.0
        qualities[parts[0]] = quality
    return qualities.get('gzip', qualities.get('*', 0.0)) > 0.0


async def _canonical_response(compressed: bytes, request: Optional[Request]) -> Response:
    headers = {"Vary": "Accept-Encoding"}
    if _accepts_gzip(request):
        headers["Content-Encoding"] = "gzip"
        content = compressed
    else:
        # A diagnostic export can be tens of MB; keep decompression off the loop.
        content = await asyncio.to_thread(gzip.decompress, compressed)
    return Response(content=content, media_type="application/json", headers=headers)


async def _delivery_response(delivery: ReportDelivery, include_diagnostics: bool,
                             request: Optional[Request] = None) -> Response:
    if include_diagnostics:
        return await _canonical_response(delivery.canonical_gzip, request)
    return Response(content=delivery.display, media_type="application/json")


@app.post("/analyze")
async def analyze(request: AnalyzeRequest, http_request: Request):
    lease = None
    if request.structureId and not request.pdbText and not request.mmcifText:
        # Acquisition transfers ownership synchronously. Cancelling a to_thread
        # waiter here could lose a lease obtained after the caller disappeared.
        lease = _lease_structure(request.structureId)
    try:
        return await _analyze_request(request, http_request, lease)
    finally:
        # The runner takes a second lease for queued/running work. This request
        # lease covers cache lookup and remains safe to release on disconnect.
        if lease:
            lease.release()


async def _analyze_request(request, http_request, lease=None):
    # Source fetching/hashing and JSON encoding are also substantial for large
    # inputs, so neither runs on the HTTP event loop.
    if lease:
        entry = lease.entry
        pdb_id = (request.pdbId or '').strip().lower() or None
        structure_text, structure_format = str(entry.path), entry.format
        key = cache_key(entry.digest, None, request.chainA, request.chainB,
                        request.mode or 'all', focus_residue=(request.focusResidue or '').strip() or None)
    else:
        pdb_id, structure_text, structure_format, key = await asyncio.to_thread(_resolve_analysis_source, request)
    cached = cache.get(key)
    if cached is not None:
        report_store.set(cached.report_id, cached.canonical_gzip)
        return await _delivery_response(cached, request.includeDiagnostics, http_request)
    if not structure_text:
        if pdb_id and pdb_id in DEMO_CHAINS:
            delivery = await asyncio.to_thread(prepare_delivery, stub_report(pdb_id, request.chainA, request.chainB))
        else:
            raise HTTPException(status_code=404, detail="Structure not available")
    else:
        worker_lease = _lease_structure(request.structureId) if lease else None
        try:
            delivery = await analysis_runner.run(
                analyze_stored_and_serialize if lease else analyze_and_serialize,
                structure_text, request.chainA, request.chainB,
                request.mode or "all", structure_format,
                (request.focusResidue or "").strip() or None, pdb_id,
                request=http_request,
                on_complete=worker_lease.release if worker_lease else None,
            )
        except AnalysisBusy as exc:
            raise HTTPException(status_code=503, detail=str(exc), headers={"Retry-After": "1"}) from exc
        except AnalysisDisconnected as exc:
            raise HTTPException(status_code=499, detail=str(exc)) from exc
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        except RuntimeError as exc:
            raise HTTPException(status_code=500, detail=str(exc)) from exc
    report_store.set(delivery.report_id, delivery.canonical_gzip)
    cache.set(key, delivery)
    return await _delivery_response(delivery, request.includeDiagnostics, http_request)


@app.get("/report/{report_id}")
async def get_report(report_id: str, http_request: Request = None):
    report = report_store.get(report_id)
    if not report:
        raise HTTPException(status_code=404, detail="Report not found")
    return await _canonical_response(report, http_request)


@app.post("/explain")
async def explain(request: ExplainRequest):
    narrative = await asyncio.to_thread(explain_report, request.report, request.images, request.notes)
    return {"narrative": narrative}


@app.post("/ribbon")
def ribbon(request: RibbonRequest):
    try:
        from .ribbon_backend import structure_to_ribbon_json
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"Ribbon backend not available: {exc}") from exc

    pdb_id = (request.pdbId or "").strip().lower() or None
    step = request.step or 0.35
    text: Optional[str] = None
    suffix = None

    if request.pdbText:
        text = request.pdbText
        suffix = ".pdb"
    elif request.mmcifText:
        text = request.mmcifText
        suffix = ".cif"
    elif pdb_id:
        text = resolve_mmcif(pdb_id, None)
        suffix = ".cif"

    if not text or not suffix:
        raise HTTPException(status_code=400, detail="pdbId, pdbText, or mmcifText is required")

    temp_path = None
    try:
        temp_path = write_temp_structure(text, suffix)
        data = structure_to_ribbon_json(str(temp_path), step=step)
        return data
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    finally:
        if temp_path and temp_path.exists():
            temp_path.unlink(missing_ok=True)


@app.post("/chapi-mesh")
def chapi_mesh(request: ChapiMeshRequest):
    # Admission starts before parsing/filtering, which used to overlap other
    # native work and leave allocations in the long-lived API process.
    with HEAVY_WORK_GATE:
        lease = (_lease_structure(request.structureId)
                 if request.structureId and not request.pdbText and not request.mmcifText else None)
        try:
            return _chapi_mesh_request(request, lease.entry if lease else None)
        finally:
            if lease:
                lease.release()


def _chapi_mesh_request(request: ChapiMeshRequest, stored=None):
    pdb_id = (request.pdbId or "").strip().lower() or None
    use_pdb_id_source_key = bool(pdb_id and not request.pdbText and not request.mmcifText)
    request_chain_ids: Optional[list[str]] = None
    low_memory_single_chain_id: Optional[str] = None
    if isinstance(request.chainIds, list):
        seen_chain_ids: set[str] = set()
        cleaned_chain_ids: list[str] = []
        for chain in request.chainIds:
            token = str(chain or "").strip()
            if not token or token in seen_chain_ids:
                continue
            seen_chain_ids.add(token)
            cleaned_chain_ids.append(token)
        if cleaned_chain_ids:
            request_chain_ids = cleaned_chain_ids
            if len(cleaned_chain_ids) == 1:
                low_memory_single_chain_id = cleaned_chain_ids[0]
    text: Optional[str] = None
    fmt = None
    if (
        CHAPI_SPLIT_CHAIN_BATCH_LIMIT > 0
        and request.representation in {"ribbon", "surface"}
        and request.splitByChain
    ):
        single_chain_cid = chapi_single_chain_from_cid(request.cid)
        if low_memory_single_chain_id is None and single_chain_cid:
            low_memory_single_chain_id = single_chain_cid
        if request_chain_ids is None and not single_chain_cid:
            raise HTTPException(
                status_code=413,
                detail={
                    "message": (
                        "CHAPI low-memory mode requires ribbon/surface requests "
                        "to be chunked with explicit chainIds."
                    ),
                    "errorCode": "CHAPI-LOWMEM-001",
                },
            )
        if request_chain_ids is not None and len(request_chain_ids) > CHAPI_SPLIT_CHAIN_BATCH_LIMIT:
            raise HTTPException(
                status_code=413,
                detail={
                    "message": (
                        "CHAPI low-memory mode rejected a multi-chain ribbon/surface batch; "
                        "request one chain at a time."
                    ),
                    "errorCode": "CHAPI-LOWMEM-002",
                },
            )

    should_filter_low_memory_single_chain = bool(
        CHAPI_LOW_MEMORY_MODE
        and request.representation in {"ribbon", "surface"}
        and request.splitByChain
        and low_memory_single_chain_id
    )

    source_path = None
    if stored:
        source_path, fmt = stored.path, stored.format
    elif request.pdbText:
        text = request.pdbText
        fmt = "pdb"
    elif request.mmcifText:
        text = request.mmcifText
        fmt = "mmcif"
    elif pdb_id:
        text = resolve_mmcif(pdb_id, None)
        fmt = "mmcif"

    if (not text and not source_path) or not fmt:
        raise HTTPException(
            status_code=400,
            detail={
                "message": "pdbId, pdbText, or mmcifText is required",
                "errorCode": "CHAPI-REQ-001",
            },
        )

    # Use the unfiltered input identity so a successful reader recovery is
    # shared by subsequent low-memory requests for other chains of this file.
    reader_source_key = (stored.digest if stored else
                         _chapi_reader_source_key({"text": text, "format": fmt}))
    if stored and should_filter_low_memory_single_chain:
        try:
            source_path, atom_rows = structure_store.chain_path(
                stored, low_memory_single_chain_id, CHAPI_PYTHON,
                _build_chapi_bridge_env(), CHAPI_WORKER_TIMEOUT_SECONDS,
            )
        except StructureCapacity as exc:
            raise HTTPException(status_code=503, detail=str(exc)) from exc
        except (ValueError, subprocess.TimeoutExpired) as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc
        if not source_path or not atom_rows:
            return Response(content='{"meshType":"chains","meshes":[],"chainIds":[],"status":"empty-chain-filter"}',
                            media_type='application/json')
    elif should_filter_low_memory_single_chain and low_memory_single_chain_id:
        original_text_length = len(text)
        filtered_text, filtered_atom_rows = filter_structure_text_to_single_chain(
            text,
            fmt,
            low_memory_single_chain_id,
        )
        if filtered_atom_rows <= 0 or not filtered_text:
            empty_payload = {
                "meshType": "chains",
                "meshes": [],
                "chainIds": [],
                "status": "empty-chain-filter",
            }
            return Response(
                content=json.dumps(empty_payload, separators=(",", ":"), ensure_ascii=False),
                media_type="application/json",
            )
        text = filtered_text
        request_chain_ids = [low_memory_single_chain_id]
        if fmt == "pdb":
            request.pdbText = None
        elif fmt == "mmcif":
            request.mmcifText = None
        if original_text_length > len(filtered_text):
            gc.collect()

    payload = {
        "text": text,
        "sourcePath": str(source_path) if source_path else None,
        "outputFormat": request.outputFormat,
        "format": fmt,
        "_reader_source_key": reader_source_key,
        "representation": request.representation,
        "mode": request.mode,
        "againstDarkBackground": request.againstDarkBackground,
        "bondWidth": request.bondWidth,
        "atomRadiusToBondWidthRatio": request.atomRadiusToBondWidthRatio,
        "smoothnessFactor": request.smoothnessFactor,
        "nonDrawCids": request.nonDrawCids,
        "carbonColor": request.carbonColor,
        "cid": request.cid,
        "colourScheme": request.colourScheme,
        "style": request.style,
        "secondaryStructureUsage": request.secondaryStructureUsage,
        "splitByChain": request.splitByChain,
        "chainIds": request_chain_ids,
    }

    try:
        source_key = stored.digest if stored else (f"pdbid:{pdb_id}" if use_pdb_id_source_key and pdb_id else None)
        mesh_cache_key = build_chapi_mesh_cache_key(payload, source_key=source_key)
        mesh_json = run_chapi_mesh_cached(payload, mesh_cache_key)
        return Response(content=mesh_json, media_type=(
            'application/vnd.roami.mesh' if mesh_json.startswith(b'ROAMIM01') else 'application/json'))
    except HTTPException:
        raise
    except Exception as exc:
        error_code, message = classify_chapi_mesh_error(exc)
        raise HTTPException(
            status_code=500,
            detail={
                "message": message,
                "errorCode": error_code,
            },
        ) from exc
    finally:
        if stored and source_path and source_path != stored.path:
            structure_store.release_chain(stored, source_path)


@app.get("/image/{report_id}/{view}")
async def get_image(report_id: str, view: str):
    raise HTTPException(status_code=501, detail="Image rendering not implemented")
