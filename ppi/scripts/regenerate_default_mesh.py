#!/usr/bin/env python3
"""Rebuild bundled 9GNQ with Roami's CHAPI worker and validate before replacing it."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
from Bio.PDB import MMCIFParser


SITE = Path(__file__).resolve().parents[2]
CHAINS = ["B", "A", "K"]
EXPECTED_STARTS = {"A:1", "A:47", "B:1", "K:8", "K:199"}


def audit(mesh_data: dict, structure: Path) -> dict:
    model = next(MMCIFParser(QUIET=True).get_structure("9GNQ", str(structure)).get_models())
    failures, results, observed_starts = [], [], set()
    meshes = mesh_data.get("meshes", [])
    if mesh_data.get("chainIds") != CHAINS or [m.get("chainId") for m in meshes] != CHAINS:
        failures.append("Expected exactly the original B, A, K chain order")
    for mesh in meshes:
        cid = mesh["chainId"]
        positions = np.asarray(mesh["positions"], dtype=float).reshape(-1, 3)
        normals = np.asarray(mesh["normals"], dtype=float).reshape(-1, 3)
        colors = np.asarray(mesh["colors"], dtype=float).reshape(-1, 4)
        indices = np.asarray(mesh["indices"]).reshape(-1, 3)
        valid = bool(
            len(positions) and len(indices)
            and mesh["vertexCount"] == len(positions) == len(normals) == len(colors)
            and mesh["triangleCount"] == len(indices)
            and all(np.isfinite(values).all() for values in (positions, normals, colors, indices))
            and np.issubdtype(indices.dtype, np.integer)
            and (indices >= 0).all() and (indices < len(positions)).all()
        )
        if not valid:
            failures.append(f"{cid}: invalid geometry or non-finite values")
        segments = []
        for residue in model[cid]:
            if residue.id[0].strip() or "CA" not in residue:
                continue
            if not segments or np.linalg.norm(residue["CA"].coord - segments[-1][-1]["CA"].coord) > 4.1:
                segments.append([])
            segments[-1].append(residue)
        distances = {}
        for segment in segments:
            if len(segment) < 3:
                continue
            for endpoint, residue in (("start", segment[0]), ("end", segment[-1])):
                label = f"{cid}:{residue.id[1]}{residue.id[2].strip()}"
                distance = float(np.linalg.norm(positions - residue["CA"].coord, axis=1).min()) if len(positions) else float("inf")
                distances[label] = {"endpoint": endpoint, "distance_angstrom": round(distance, 6)}
                if endpoint == "start":
                    observed_starts.add(label)
                if not np.isfinite(distance) or distance > 0.6:
                    failures.append(f"{label}: terminal CA is {distance:.6f} Angstrom from the mesh")
        results.append({"chain": cid, "vertices": len(positions), "triangles": len(indices),
                        "valid_geometry": valid, "termini": distances})
    if observed_starts != EXPECTED_STARTS:
        failures.append(f"Unexpected segment starts: {sorted(observed_starts)}")
    return {"passed": not failures, "max_terminal_distance_angstrom": 0.6,
            "failures": failures, "results": results}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=SITE / "public/ppi/data/9GNQ.chapi.json")
    parser.add_argument("--report", type=Path, default=SITE / "roami-tests/coot-ribbon-termini/default-mesh-validation.json")
    args = parser.parse_args()
    structure = SITE / "public/ppi/data/9GNQ.cif"
    lock = json.loads((SITE / "ppi/coot-ribbon-source.json").read_text())
    previous_bytes = args.output.read_bytes() if args.output.exists() else None
    previous = audit(json.loads(previous_bytes), structure) if previous_bytes else None

    # A one-shot worker avoids keeping a native subprocess alive after generation.
    os.environ["CHAPI_PERSISTENT_WORKER"] = "0"
    sys.path.insert(0, str(SITE / "ppi"))
    from api.main import run_chapi_mesh

    options = {"representation": "ribbon", "cid": "//", "colourScheme": "colorRampChainsScheme",
               "style": "Ribbon", "secondaryStructureUsage": 2, "splitByChain": True, "chainIds": CHAINS}
    data = json.loads(run_chapi_mesh({"text": structure.read_text(), "format": "mmcif", **options}))
    validation = audit(data, structure)
    if not validation["passed"]:
        raise RuntimeError("Generated mesh failed validation; original asset retained:\n" + "\n".join(validation["failures"]))
    data["generation"] = {
        "structureSha256": hashlib.sha256(structure.read_bytes()).hexdigest(),
        "configuredCootSourceCommit": lock["commit"],
        "options": options,
    }
    encoded = (json.dumps(data, separators=(",", ":"), allow_nan=False) + "\n").encode()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    staged = args.output.with_suffix(args.output.suffix + ".tmp")
    staged.write_bytes(encoded)
    staged.replace(args.output)
    report = {"asset": args.output.name, "sha256": hashlib.sha256(encoded).hexdigest(),
              "generation": data["generation"], "validation": validation,
              "previous_sha256": hashlib.sha256(previous_bytes).hexdigest() if previous_bytes else None,
              "previous_validation": previous}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"Rebuilt {args.output.name}: {len(data['meshes'])} chains; all five starts and segment ends pass.")
    print(f"Validation report: {args.report}")


if __name__ == "__main__":
    main()
