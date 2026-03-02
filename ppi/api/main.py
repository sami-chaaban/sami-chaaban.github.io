"""
FastAPI backend for the protein-protein interaction demo.

This backend exposes routes for listing chains, analyzing interfaces, and
returning reports. The analysis logic uses a lightweight, deterministic
"toy" engine so the demo runs without heavy scientific dependencies.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Optional
import json
import os
import re
import subprocess
import tempfile
import threading
import time
from pathlib import Path
import sys
import uuid
from urllib import error as urlerror
from urllib import parse as urlparse
from urllib import request as urlrequest

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

from .analysis import (
    analyze_interface,
    cache_key,
    fetch_mmcif,
    list_chains,
)
from .cache import ReportCache
from .explain import explain_report
from .models import AnalyzeRequest, ChainsRequest, ExplainRequest, RibbonRequest, ChapiMeshRequest


DEMO_CHAINS = {
    "4hhb": ["A", "B"],
    "1a3n": ["A", "B"],
}

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_URL_TEMPLATE = "https://rest.uniprot.org/uniprotkb/{accession}.json"
PD_BE_SUMMARY_URL_TEMPLATE = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/{pdb_id}"
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

cache = ReportCache()
report_store: dict = {}

CHAPI_PYTHON = os.environ.get("CHAPI_PYTHON") or sys.executable
CHAPI_PREFIX = os.environ.get("COOT_PREFIX") or os.environ.get("CONDA_PREFIX") or ""
CHAPI_BRIDGE = Path(__file__).with_name("chapi_bridge.py")


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
            "/protein-search",
            "/protein-structures/{accession}",
        ],
    }


@app.get("/health")
async def health():
    return {"status": "ok"}


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
        "referenceAuthors": reference["referenceAuthors"],
        "referenceJournal": reference["referenceJournal"],
        "referenceYear": reference["referenceYear"],
    }


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
        "referenceAuthors": reference["referenceAuthors"],
        "referenceJournal": reference["referenceJournal"],
        "referenceYear": reference["referenceYear"],
    }


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
            cache_is_missing_critical_fields = (
                not cached_title
                or cached_method in ("", "—")
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
    # PDBe summary frequently omits resolution/method/citation details. Fill from RCSB.
    if (
        summary["resolution"] is None
        or summary["method"] == "—"
        or summary["title"] == "—"
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


def fetch_uniprot_search(query: str, reviewed: bool = False) -> dict[str, Any]:
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


def fetch_uniprot_entry(accession: str) -> dict[str, Any]:
    normalized = collapse_whitespace(accession).upper()
    if not normalized:
        raise RuntimeError("No UniProt accession provided")
    url = UNIPROT_ENTRY_URL_TEMPLATE.format(accession=urlparse.quote(normalized, safe=""))
    return fetch_json(url)


@app.get("/protein-search")
async def protein_search(query: str, reviewed: bool = False):
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
async def protein_structures(accession: str):
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
    return fetch_mmcif(pdb_id)


def write_temp_structure(text: str, suffix: str) -> Path:
    temp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    temp.write(text.encode("utf-8"))
    temp.flush()
    temp.close()
    return Path(temp.name)


def run_chapi_mesh(payload: dict) -> dict:
    if not CHAPI_BRIDGE.exists():
        raise RuntimeError("chapi_bridge.py not found in api directory.")
    env = os.environ.copy()
    if CHAPI_PREFIX:
        env.setdefault("COOT_PREFIX", CHAPI_PREFIX)
    env.setdefault("PYTHONUNBUFFERED", "1")
    proc = subprocess.run(
        [CHAPI_PYTHON, str(CHAPI_BRIDGE)],
        input=json.dumps(payload),
        text=True,
        capture_output=True,
        env=env,
        timeout=120,
    )
    if proc.returncode != 0:
        detail = proc.stderr.strip() or proc.stdout.strip() or "chapi bridge failed"
        raise RuntimeError(detail)
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Invalid JSON from chapi bridge: {exc}") from exc


def stub_report(pdb_id: str, chain_a: str, chain_b: str) -> dict:
    contact_template = {
        "residueA": {"chain": chain_a, "resName": "LEU", "seq": "42", "atom": "CD1"},
        "residueB": {"chain": chain_b, "resName": "VAL", "seq": "10", "atom": "CG1"},
        "distance": 3.8,
        "type": "hydrophobic",
    }
    contacts = {
        "hydrogen_bonds": [],
        "salt_bridges": [],
        "hydrophobic": [contact_template],
        "metal_coordination": [],
        "pi_pi": [],
        "pi_cation": [],
        "other": [],
    }
    per_residue = {
        f"{chain_a}:42": {
            "chain": chain_a,
            "resName": "LEU",
            "seq": "42",
            "hydrophobic": 1,
            "hbond": 0,
            "salt_bridge": 0,
            "metal_coordination": 0,
            "pi_pi": 0,
            "pi_cation": 0,
            "other": 0,
            "total": 1,
        },
        f"{chain_b}:10": {
            "chain": chain_b,
            "resName": "VAL",
            "seq": "10",
            "hydrophobic": 1,
            "hbond": 0,
            "salt_bridge": 0,
            "metal_coordination": 0,
            "pi_pi": 0,
            "pi_cation": 0,
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
        "buriedFraction": {chain_a: 0.5, chain_b: 0.5},
        "approxDeltaG": None,
        "meta": {
            "engine": "stub",
            "note": "Stub report returned because structure could not be fetched.",
        },
    }


@app.get("/chains")
async def get_chains(pdbId: str):
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
async def post_chains(request: ChainsRequest):
    if not request.mmcifText:
        raise HTTPException(status_code=400, detail="mmcifText is required")
    chains, aliases = list_chains(request.mmcifText)
    return {"chains": chains, "aliases": aliases.label_to_auth}


@app.post("/analyze")
async def analyze(request: AnalyzeRequest):
    pdb_id = (request.pdbId or "").strip().lower() or None
    chain_a = request.chainA
    chain_b = request.chainB
    mode = request.mode or "all"

    structure_text: Optional[str] = None
    structure_format: Optional[str] = None

    if request.pdbText:
        structure_text = request.pdbText
        structure_format = "pdb"
    elif request.mmcifText:
        structure_text = request.mmcifText
        structure_format = "mmcif"
    else:
        try:
            structure_text = resolve_mmcif(pdb_id, None)
            structure_format = "mmcif"
        except Exception:
            structure_text = None
            structure_format = None

    if not structure_text:
        if pdb_id and pdb_id in DEMO_CHAINS:
            report = stub_report(pdb_id, chain_a, chain_b)
        else:
            raise HTTPException(status_code=404, detail="Structure not available")
    else:
        cache_source = (
            structure_text
            if pdb_id
            else f"{structure_format or 'unknown'}\n{structure_text}"
        )
        key = cache_key(pdb_id, cache_source, chain_a, chain_b, mode)
        cached = cache.get(key)
        if cached:
            return cached
        try:
            report = analyze_interface(
                structure_text,
                chain_a,
                chain_b,
                mode,
                structure_format=structure_format or "mmcif",
            )
            if pdb_id:
                report["pdbId"] = pdb_id
        except ValueError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    report_id = uuid.uuid4().hex[:10]
    report["reportId"] = report_id
    report_store[report_id] = report
    cache_source = (
        structure_text
        if pdb_id
        else f"{structure_format or 'unknown'}\n{structure_text or ''}"
    )
    cache.set(cache_key(pdb_id, cache_source, chain_a, chain_b, mode), report)
    return report


@app.get("/report/{report_id}")
async def get_report(report_id: str):
    report = report_store.get(report_id)
    if not report:
        raise HTTPException(status_code=404, detail="Report not found")
    return report


@app.post("/explain")
async def explain(request: ExplainRequest):
    narrative = explain_report(request.report, request.images, request.notes)
    return {"narrative": narrative}


@app.post("/ribbon")
async def ribbon(request: RibbonRequest):
    try:
        from ribbon_backend import structure_to_ribbon_json
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
async def chapi_mesh(request: ChapiMeshRequest):
    pdb_id = (request.pdbId or "").strip().lower() or None
    text: Optional[str] = None
    fmt = None

    if request.pdbText:
        text = request.pdbText
        fmt = "pdb"
    elif request.mmcifText:
        text = request.mmcifText
        fmt = "mmcif"
    elif pdb_id:
        text = resolve_mmcif(pdb_id, None)
        fmt = "mmcif"

    if not text or not fmt:
        raise HTTPException(status_code=400, detail="pdbId, pdbText, or mmcifText is required")

    payload = {
        "text": text,
        "format": fmt,
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
    }

    try:
        return run_chapi_mesh(payload)
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc


@app.get("/image/{report_id}/{view}")
async def get_image(report_id: str, view: str):
    raise HTTPException(status_code=501, detail="Image rendering not implemented")
