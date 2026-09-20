# Analysis runtime

Start the backend with `CHAPI_PYTHON=/path/to/chapi/bin/python ./start.sh` from this directory. Local startup builds the reviewed Arpeggio patch into the ignored `.arpeggio/lib` package overlay and adds it to `PYTHONPATH`; the shared Python environment remains unchanged. Docker applies the same patch while building the image. Direct `uvicorn` invocation bypasses local overlay preparation unless its `PYTHONPATH` has been configured explicitly.

The dependency is pinned to `pdbe-arpeggio==1.4.4` and `openbabel-wheel==3.1.1.23`. `scripts/patch_arpeggio.py` checks the exact reviewed upstream source hash and applies only the selection-set hoist. It accepts an already patched source only if reversing that edit restores the original hash. No running Python methods are replaced.

Analysis runs in spawned worker processes. `ANALYZE_WORKERS` defaults to 1; `ANALYZE_MAX_PENDING` defaults to 8, including active requests. Excess admission returns HTTP 503 with `Retry-After: 1`. Disconnected waiting requests are discarded before native execution. A job already executing finishes in its existing process and retains its worker slot until completion, even if its client disconnects. One worker is the recommended memory bound for the existing low-memory deployment. Coot mesh jobs retain their separate existing serial execution bound.

`POST /analyze` includes diagnostic contacts by default for compatibility. Setting `includeDiagnostics: false` omits only records explicitly flagged `debugOnly`, `debug_only`, or `asserted.debugOnly`; low-confidence scientific contacts are retained. The complete canonical result remains at `GET /report/{reportId}`. A compact response's `meta.diagnostics` records omission counts by bucket and its canonical report URL. Canonical and compact JSON are encoded off the event loop; the cache is independent of projection choice. Reports expire according to `REPORT_STORE_TTL_SECONDS` (default six hours); another cached analyze response refreshes access to that canonical report.

`contactingResidueFraction` describes the fraction of the contacting residues assigned to each chain. `buriedFraction`, `interfaceArea`, and `approxDeltaG` are null when those physical quantities have not been computed.

Run scheduling, source-identity, delivery, CIF, explanation and build-patch regressions:

```sh
PYTHONPATH=. /path/to/chapi/bin/python -m unittest api.test_api_audit -v
```

Validate the actual built package against stock Arpeggio on the local DNA, calmodulin and heme fixtures:

```sh
/path/to/chapi/bin/python scripts/patch_arpeggio.py --overlay .arpeggio/lib
/path/to/chapi/bin/python scripts/verify_arpeggio_equivalence.py --output /tmp/arpeggio-equivalence.json
```

The verifier starts fresh interpreters for each side and compares all raw contacts and the complete Roami report, canonicalizing contact ordering only. This establishes output equivalence for those fixtures, not independent chemical ground truth.
