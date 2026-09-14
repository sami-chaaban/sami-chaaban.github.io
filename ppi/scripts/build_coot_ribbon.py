#!/usr/bin/env python3
"""Build and verify Roami's pinned Coot ribbon library without updating chapi."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request


LIBRARY = "libcootmoleculestotriangles.so.1.1"


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def coot_package(prefix: Path) -> dict:
    records = list((prefix / "conda-meta").glob("coot-headless-*.json"))
    if len(records) != 1:
        raise RuntimeError(f"Expected exactly one coot-headless package in {prefix}")
    package = json.loads(records[0].read_text())
    return {key: package[key] for key in ("name", "version", "build", "build_number", "subdir")}


def check_package(prefix: Path, lock: dict) -> dict:
    package = coot_package(prefix)
    if package["version"] != lock["coot_version"] or package["build_number"] != lock["coot_build_number"]:
        raise RuntimeError(f"Coot package does not match source lock: {package}")
    return package


def extract_source(archive: Path, destination: Path) -> Path:
    """Extract regular files/directories only; reject traversal and special files."""
    with tarfile.open(archive, "r:gz") as source:
        for member in source.getmembers():
            target = destination / member.name
            if not target.resolve().is_relative_to(destination.resolve()):
                raise RuntimeError(f"Unsafe source archive path: {member.name}")
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                target.parent.mkdir(parents=True, exist_ok=True)
                with source.extractfile(member) as incoming, target.open("wb") as outgoing:
                    shutil.copyfileobj(incoming, outgoing)
                target.chmod(member.mode & 0o777)
            else:
                raise RuntimeError(f"Unexpected non-regular source archive entry: {member.name}")
    roots = list(destination.iterdir())
    if len(roots) != 1 or not (roots[0] / "roami/CMakeLists.txt").is_file():
        raise RuntimeError("Source archive must contain one Coot checkout with the Roami build driver")
    return roots[0]


def run(arguments: list[str], **kwargs) -> None:
    subprocess.run(arguments, check=True, **kwargs)


def verify(prefix: Path, manifest_path: Path) -> None:
    """Verify the actual shared object imported by CHAPI in this process."""
    manifest = json.loads(manifest_path.read_text())
    if coot_package(prefix) != manifest["coot_package"]:
        raise RuntimeError("Installed Coot package differs from the build environment")
    library = prefix / "lib" / LIBRARY
    if digest(library) != manifest["library_sha256"]:
        raise RuntimeError("Installed ribbon library does not match the build manifest")
    import coot_headless_api  # noqa: F401: importing forces the native dependencies to load

    maps = Path("/proc/self/maps")
    if not maps.is_file():
        raise RuntimeError("Shared-library verification requires Linux /proc/self/maps")
    loaded = {
        line.split(maxsplit=5)[5].strip()
        for line in maps.read_text().splitlines()
        if "libcootmoleculestotriangles" in line
    }
    if loaded != {str(library.resolve())}:
        raise RuntimeError(f"CHAPI loaded an unexpected ribbon library: {sorted(loaded)}")
    print(f"Verified CHAPI loaded {library} ({manifest['source']['commit']})")


def build(args: argparse.Namespace) -> None:
    lock = json.loads(args.lock.read_text())
    if lock["repository"] != "https://github.com/sami-chaaban/coot":
        raise RuntimeError("Unexpected Coot source repository")
    if not re.fullmatch(r"[0-9a-f]{40}", lock["commit"]):
        raise RuntimeError("The source lock needs an immutable 40-character Git commit")
    if not re.fullmatch(r"[0-9a-f]{64}", lock["archive_sha256"]):
        raise RuntimeError("The source lock needs a SHA256 digest of the source archive")
    package = check_package(args.prefix, lock)
    if not package["subdir"].startswith("linux-"):
        raise RuntimeError("This backend build targets Linux; use Coot's roami driver directly on macOS")
    args.output.mkdir(parents=True, exist_ok=True)
    archive_url = f"{lock['repository']}/archive/{lock['commit']}.tar.gz"
    with tempfile.TemporaryDirectory(prefix="roami-coot-") as work_name:
        work = Path(work_name)
        archive = work / "coot.tar.gz"
        request = urllib.request.Request(archive_url, headers={"User-Agent": "Roami-Coot-Build/1"})
        with urllib.request.urlopen(request, timeout=120) as response, archive.open("wb") as target:
            shutil.copyfileobj(response, target)
        if digest(archive) != lock["archive_sha256"]:
            raise RuntimeError("Coot source archive failed SHA256 verification")
        source_dir = work / "source"
        source_dir.mkdir()
        source = extract_source(archive, source_dir)
        build_dir = work / "build"
        run([
            "cmake", "-S", str(source / "roami"), "-B", str(build_dir), "-G", "Ninja",
            "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_C_COMPILER=/usr/bin/gcc", "-DCMAKE_CXX_COMPILER=/usr/bin/g++",
            f"-DCOOT_DEPENDENCY_PREFIX={args.prefix}", f"-DCOOT_BOOST_INCLUDE_DIR={args.boost_include}",
            "-DCMAKE_BUILD_WITH_INSTALL_RPATH=ON",
        ])
        run(["cmake", "--build", str(build_dir), "--parallel", str(args.jobs)])
        run(["ctest", "--test-dir", str(build_dir), "--output-on-failure"])
        built_library = build_dir / LIBRARY
        if not built_library.is_file():
            raise RuntimeError(f"Coot build did not produce {LIBRARY}")
        output_library = args.output / LIBRARY
        shutil.copy2(built_library, output_library)
        manifest = {
            "source": lock,
            "archive_url": archive_url,
            "coot_package": package,
            "library": LIBRARY,
            "library_sha256": digest(output_library),
        }
        manifest_path = args.output / "coot-ribbon-build.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        license_dir = args.output / "licenses"
        license_dir.mkdir()
        for filename in ("COPYING", "COPYING.LESSERv3"):
            shutil.copy2(source / filename, license_dir / filename)
        # Replace only in the disposable builder, then test through CHAPI. The
        # final stage performs its own load verification after copying the file.
        installed_library = args.prefix / "lib" / LIBRARY
        if not installed_library.is_file():
            raise RuntimeError(f"Expected existing Coot shared library: {installed_library}")
        installed_library.unlink()
        shutil.copy2(output_library, installed_library)
        child_env = os.environ.copy()
        child_env.update({
            "COOT_PREFIX": str(args.prefix),
            "COOT_REFMAC_LIB_DIR": str(args.prefix / "share/coot/lib"),
            "COOT_MONOMER_LIB_DIR": str(args.prefix / "share/coot/lib/data/monomers"),
            "CLIBD_MON": str(args.prefix / "share/coot/lib/data/monomers"),
            "COOT_STANDARD_RESIDUES": str(args.prefix / "share/coot/standard-residues.pdb"),
            "SYMINFO": str(args.prefix / "share/coot/data/syminfo.lib"),
        })
        python = str(args.prefix / "bin/python")
        run([python, str(Path(__file__).resolve()), "verify", "--prefix", str(args.prefix), "--manifest", str(manifest_path)], env=child_env)
        run([python, str(source / "roami/mesh_regression.py"), "--structure", str(args.structure), "--output", str(args.output / "mesh-regression.json")], env=child_env)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    builder = commands.add_parser("build")
    builder.add_argument("--lock", type=Path, required=True)
    builder.add_argument("--prefix", type=Path, required=True)
    builder.add_argument("--boost-include", type=Path, required=True)
    builder.add_argument("--structure", type=Path, required=True)
    builder.add_argument("--output", type=Path, required=True)
    builder.add_argument("--jobs", type=int, default=2)
    verifier = commands.add_parser("verify")
    verifier.add_argument("--prefix", type=Path, required=True)
    verifier.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "build":
        if args.jobs < 1:
            parser.error("--jobs must be positive")
        build(args)
    else:
        verify(args.prefix, args.manifest)


if __name__ == "__main__":
    main()
