#!/usr/bin/env python3
"""Apply the reviewed selection-set hoist to pinned PDBe Arpeggio 1.4.4.

Docker uses --installed at build time. Local startup uses --overlay DIR, which
copies the package into a project-owned import directory without changing the
shared conda environment. No methods are replaced in a running interpreter.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import shutil


VERSION = "1.4.4"
ORIGINAL_SHA256 = "2d867a59a2778d6979fbf200951171a684ab842cdf758830261f4b6802001ca8"
OLD = "        for atom_bgn, atom_end in self.ns.search_all(interacting_cutoff):\n\n            selection_set = set(self.selection)\n"
NEW = "        selection_set = set(self.selection)\n        for atom_bgn, atom_end in self.ns.search_all(interacting_cutoff):\n"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def patched_source(source: bytes) -> bytes:
    if digest(source) == ORIGINAL_SHA256:
        if source.count(OLD.encode()) != 1:
            raise RuntimeError("Pinned Arpeggio source does not contain exactly one reviewed patch site")
        return source.replace(OLD.encode(), NEW.encode(), 1)
    # Permit repeat builds only if reversing this one exact edit restores the
    # audited upstream source; unrelated source changes must fail closed.
    if source.count(NEW.encode()) == 1 and digest(source.replace(NEW.encode(), OLD.encode(), 1)) == ORIGINAL_SHA256:
        return source
    raise RuntimeError("Arpeggio source hash differs from the reviewed 1.4.4 dependency; refusing to patch")


def build(overlay: Path | None = None) -> dict:
    version = importlib.metadata.version("pdbe-arpeggio")
    if version != VERSION:
        raise RuntimeError(f"Expected pdbe-arpeggio=={VERSION}, found {version}")
    spec = importlib.util.find_spec("arpeggio")
    if spec is None or spec.origin is None:
        raise RuntimeError("Cannot locate installed arpeggio package")
    source_package = Path(spec.origin).parent.resolve()
    source = (source_package / "core" / "interactions.py").read_bytes()
    patched = patched_source(source)
    target_package = source_package
    if overlay is not None:
        target_package = overlay.resolve() / "arpeggio"
        if target_package != source_package:
            shutil.copytree(source_package, target_package, dirs_exist_ok=True, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    target = target_package / "core" / "interactions.py"
    target.write_bytes(patched)
    # Never allow stale copied bytecode to shadow the patched source.
    for bytecode in (target.parent / "__pycache__").glob("interactions.*.pyc"):
        bytecode.unlink()
    manifest = {
        "package": "pdbe-arpeggio", "version": VERSION,
        "patch": "selection-set-hoist-v1", "upstreamSha256": ORIGINAL_SHA256,
        "patchedSha256": digest(patched),
    }
    manifest_path = target_package.parent / "roami-arpeggio-build.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument("--overlay", type=Path)
    target.add_argument("--installed", action="store_true")
    args = parser.parse_args()
    print(json.dumps(build(args.overlay), indent=2))
