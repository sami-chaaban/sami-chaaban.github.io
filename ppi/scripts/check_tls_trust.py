#!/usr/bin/env python3
"""Check this Python environment's patched CA packages and optional HTTPS trust."""

import argparse
import json
from pathlib import Path
import re
import ssl
import sys
from urllib import error, request


MINIMUM_VERSION = (2026, 10, 19)
CONTROL_URL = "https://www.python.org/"
COMPROMISED_URLS = (
    "https://dv-g3b-emudhra.we-have-a-problem.com/",
    "https://dv-g3b-aaa-emudhra.we-have-a-problem.com/",
)


def require_patched(name: str, version: str) -> None:
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        raise RuntimeError(f"Cannot validate {name} version: {version!r}")
    if tuple(map(int, version.split("."))) < MINIMUM_VERSION:
        raise RuntimeError(f"{name} {version} is older than required 2026.10.19")


def check_packages() -> str:
    print(f"Python: {sys.version.split()[0]} ({sys.executable})")
    print(f"Prefix: {sys.prefix}")
    print(f"Default TLS paths: {ssl.get_default_verify_paths()}")
    try:
        import certifi
    except ImportError as exc:
        raise RuntimeError("certifi is not installed in this Python environment") from exc
    print(f"certifi: {certifi.__version__}; bundle: {certifi.where()}")
    require_patched("certifi", certifi.__version__)
    try:
        from pip._vendor import certifi as pip_certifi
    except ImportError as exc:
        raise RuntimeError("Cannot validate pip's bundled certifi") from exc
    print(f"pip certifi: {pip_certifi.__version__}; bundle: {pip_certifi.where()}")
    require_patched("pip's bundled certifi", pip_certifi.__version__)
    records = []
    for path in (Path(sys.prefix) / "conda-meta").glob("ca-certificates-*.json"):
        record = json.loads(path.read_text())
        if record.get("name") == "ca-certificates":
            records.append(record)
    if len(records) != 1:
        raise RuntimeError("Expected one ca-certificates package in this prefix's conda-meta")
    version = records[0]["version"]
    print(f"Conda ca-certificates: {version}")
    require_patched("Conda ca-certificates", version)
    return certifi.where()


def classify_error(exc: Exception) -> str:
    # URLError normally wraps SSL failures in reason. Strings and unrelated
    # connection failures are never evidence that a certificate was rejected.
    cause = exc.reason if isinstance(exc, error.URLError) else exc
    if isinstance(cause, ssl.SSLCertVerificationError):
        return "certificate_rejected"
    return "request_failed"


class NoRedirect(request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        # A different host's certificate cannot prove this host was rejected.
        return None


def probe(url: str, context: ssl.SSLContext) -> tuple[str, str]:
    opener = request.build_opener(request.HTTPSHandler(context=context), NoRedirect())
    try:
        with opener.open(url, timeout=20) as response:
            status = response.status
            return ("verified" if 200 <= status < 300 else "request_failed", f"HTTP {status}")
    except error.HTTPError as exc:
        return "request_failed", f"HTTP {exc.code}"
    except (error.URLError, OSError) as exc:
        return classify_error(exc), str(exc)


def check_network(certifi_bundle: str) -> None:
    failures = []
    contexts = (
        ("default OpenSSL", ssl.create_default_context()),
        ("certifi", ssl.create_default_context(cafile=certifi_bundle)),
    )
    for label, context in contexts:
        for url, expected in ((CONTROL_URL, "verified"), *(
            (url, "certificate_rejected") for url in COMPROMISED_URLS
        )):
            outcome, detail = probe(url, context)
            print(f"{label}: {url}: {outcome} ({detail})")
            if outcome != expected:
                failures.append(f"{label}: {url}: expected {expected}, got {outcome}")
    if failures:
        raise RuntimeError("HTTPS trust checks failed:\n" + "\n".join(failures))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", action="store_true", help="also check live HTTPS trust")
    args = parser.parse_args()
    try:
        bundle = check_packages()
        if args.network:
            check_network(bundle)
    except (RuntimeError, OSError, ValueError, KeyError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1
    print("PASS: patched CA packages" + (" and HTTPS trust" if args.network else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
