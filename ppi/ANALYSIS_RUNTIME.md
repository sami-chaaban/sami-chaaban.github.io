# Analysis runtime

Start the backend with `CHAPI_PYTHON=/path/to/chapi/bin/python ./start.sh` from this directory. Local startup builds the reviewed Arpeggio patch into the ignored `.arpeggio/lib` package overlay and adds it to `PYTHONPATH`; the shared Python environment remains unchanged. Docker applies the same patch while building the image. Direct `uvicorn` invocation bypasses local overlay preparation unless its `PYTHONPATH` has been configured explicitly.

The dependency is pinned to `pdbe-arpeggio==1.4.4` and `openbabel-wheel==3.1.1.23`. `scripts/patch_arpeggio.py` checks the exact reviewed upstream source hash and applies only the selection-set hoist. It accepts an already patched source only if reversing that edit restores the original hash. No running Python methods are replaced.

Analysis runs in spawned worker processes. `ANALYZE_WORKERS` defaults to 1; `ANALYZE_MAX_PENDING` defaults to 8, including active requests. Excess admission returns HTTP 503 with `Retry-After: 1`. Disconnected waiting requests are discarded before native execution. A job already executing finishes in its existing process and retains its worker slot until process exit, even if its client disconnects.

`ANALYZE_RECYCLE_WORKERS=1` is the default: each uncached analysis gets a fresh process, releasing scientific caches and native allocator memory on exit. This works on Python 3.10 without `max_tasks_per_child`. Set it to `0` only when deliberately trading retained memory for warm worker reuse. `HEAVY_WORKERS` defaults to 1 and shares admission between analyses and mesh preprocessing/generation. Increasing `ANALYZE_WORKERS` alone does not increase simultaneous heavy work. Cached analysis responses bypass native execution.

## Coordinate and mesh delivery

The browser uploads each actual coordinate source once with `POST /structures?format=mmcif` (or `pdb`), using a raw UTF-8 text body. The response contains an opaque `structureId`, `format`, and `expiresIn`. Subsequent `/chapi-mesh` and `/analyze` requests can supply that handle instead of inline text. Optional `pdbId` is report metadata; a handle never fetches replacement coordinates by accession. Inline source fields remain supported and retain precedence for compatibility.

Uploads live in private temporary files. Defaults: `STRUCTURE_UPLOAD_MAX_BYTES=67108864` (64 MiB per upload), `STRUCTURE_STORE_MAX_BYTES=268435456` (256 MiB including prepared chains), `STRUCTURE_STORE_MAX_ENTRIES=8`, and `STRUCTURE_STORE_TTL_SECONDS=1800` (30 minutes idle). Eviction/expiry is checked on access. Active jobs hold leases so files survive eviction requests and client disconnects until consumers exit. `DELETE /structures/{structureId}` releases an unused structure; unknown/expired handles return 410 on use. The browser shares one replacement upload and retries once, and releases handles when replacing data. API restarts naturally invalidate old handles.

In low-memory mode, the first ribbon/surface request prepares shared metadata and chain atom fragments in a disposable subprocess. Gemmi parses mmCIF once, retaining auth/label chain matching, all models, raw CIF token quoting, and non-atom metadata. Later chain requests assemble one temporary native input from those files, then remove it after rendering; metadata is stored once and the API process never constructs the full parsed structure. Input and derived files are bounded on disk rather than held in Python caches.

`/chapi-mesh` defaults to JSON for existing clients. `outputFormat: "binary"` selects `application/vnd.roami.mesh`. The `ROAMIM01` envelope contains an eight-byte magic, little-endian uint32 JSON-metadata length, UTF-8 metadata padded to a four-byte boundary, and packed little-endian Float32 position/normal/RGBA and Uint32 index buffers. Array descriptors contain byte offsets relative to the packed section, element counts, and types. Native vectors pack directly without intermediate Python numeric lists. Browser typed-array views share the response buffer; geometry is unchanged at its existing Float32 rendering precision. Gzip remains negotiated for either representation. Frontends fall back to inline JSON when an older backend rejects `/structures` with 404/405.

Worker recycling bounds retained memory, not the peak of a single scientific job. Large full-structure Arpeggio analyses can still exceed a small container's resources; these delivery changes do not alter their scientific scope or impose a new analysis timeout.

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
