#!/usr/bin/env python3
"""Compare complete scientific outputs with stock and build-patched Arpeggio.

Runs each side in a fresh process; no runtime monkeypatching or contact pruning.
Only contact ordering is canonicalized before hashing. Input fixtures are local.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time


PPI = Path(__file__).resolve().parents[1]


def encode(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def digest(value):
    return hashlib.sha256(encode(value).encode()).hexdigest()


def capture(package_root: Path, fixture: Path, chain_a: str, chain_b: str):
    sys.path.insert(0, str(PPI))
    sys.path.insert(0, str(package_root))
    from api import analysis
    started = time.perf_counter()
    report = analysis.analyze_interface(fixture.read_text(), chain_a, chain_b)
    elapsed = time.perf_counter() - started
    raw = next(iter(analysis.ARPEGGIO_CONTACTS_CACHE._entries.values())).value
    canonical = {**report, "contacts": {bucket: sorted(records, key=encode) for bucket, records in report["contacts"].items()}}
    from api.analysis_worker import prepare_delivery
    delivery = prepare_delivery(report)
    return {
        "seconds": elapsed, "rawCount": len(raw), "rawSha256": digest(sorted(raw, key=encode)),
        "reportSha256": digest(canonical), "canonicalBytes": len(delivery.canonical),
        "displayBytes": len(delivery.display),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overlay", type=Path, default=PPI / ".arpeggio" / "lib")
    parser.add_argument("--fixture-dir", type=Path, default=PPI.parent / "roami-tests" / "audit-2026-09-20" / "chemistry")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--capture", nargs=4, metavar=("PACKAGE_ROOT", "FIXTURE", "CHAIN_A", "CHAIN_B"))
    args = parser.parse_args()
    if args.capture:
        root, fixture, a, b = args.capture
        print(json.dumps(capture(Path(root), Path(fixture), a, b)))
        return
    stock = Path(importlib.util.find_spec("arpeggio").origin).parent.parent
    if stock.resolve() == args.overlay.resolve():
        raise RuntimeError("Run this verifier from the stock environment, without the overlay on PYTHONPATH")
    results = []
    for pdb_id, a, b in [("1bna", "A", "B"), ("1cll", "A", "A"), ("1a6m", "A", "A")]:
        fixture = args.fixture_dir / f"{pdb_id}.cif"
        comparisons = {}
        # Arpeggio serializes a few internal sets as lists; use the same hash seed
        # on both sides instead of changing those lists during comparison.
        child_env = {**os.environ, "PYTHONHASHSEED": "0"}
        for label, root in [("stock", stock), ("patched", args.overlay)]:
            completed = subprocess.run(
                [sys.executable, str(Path(__file__).resolve()), "--capture", str(root), str(fixture), a, b],
                capture_output=True, text=True, timeout=120, check=True,
                env=child_env,
            )
            comparisons[label] = json.loads(completed.stdout.strip().splitlines()[-1])
        before, after = comparisons["stock"], comparisons["patched"]
        result = {
            "pdbId": pdb_id, **comparisons, "speedup": before["seconds"] / after["seconds"],
            "rawEqual": before["rawSha256"] == after["rawSha256"],
            "reportEqual": before["reportSha256"] == after["reportSha256"],
        }
        results.append(result)
        print(json.dumps(result), flush=True)
        if not result["rawEqual"] or not result["reportEqual"]:
            raise RuntimeError(f"Scientific equivalence failed for {pdb_id}")
    if args.output:
        args.output.write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
