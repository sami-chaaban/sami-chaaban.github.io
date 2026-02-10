"""
FastAPI backend for the protein-protein interaction demo.

This backend exposes routes for listing chains, analyzing interfaces, and
returning reports. The analysis logic uses a lightweight, deterministic
"toy" engine so the demo runs without heavy scientific dependencies.
"""

from __future__ import annotations

from typing import Optional
import json
import os
import subprocess
import tempfile
from pathlib import Path
import sys
import uuid

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

from .analysis import analyze_interface, cache_key, fetch_mmcif, list_chains
from .cache import ReportCache
from .explain import explain_report
from .models import AnalyzeRequest, ChainsRequest, ExplainRequest, RibbonRequest, ChapiMeshRequest


DEMO_CHAINS = {
    "4hhb": ["A", "B"],
    "1a3n": ["A", "B"],
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
        "routes": ["/chains", "/analyze", "/chapi-mesh", "/explain", "/health"],
    }


@app.get("/health")
async def health():
    return {"status": "ok"}


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
