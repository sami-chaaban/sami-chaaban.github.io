# Structure loading and Arpeggio input compatibility

The frontend batches up to eight ribbon chains per request instead of repeatedly uploading and parsing the complete structure for every chain. It checks that every requested chain is present and restores the original chain order. A backend that returns a Coot low-memory batch-limit response is retried one chain at a time; that limit is remembered for the page session. Large-structure mode retains single-chain requests. Network failures stop the remaining batches, and canceled loads cannot render a partial model.

The backend negotiates lossless gzip compression, avoids a redundant Coot dictionary initialization, and reads each native vertex vector once during JSON conversion. Mesh generation runs outside the HTTP event loop. A lock still permits only one native job at a time, preserving the previous memory bound.

Testing the larger public entry 6VXX also exposed a native crash in Coot's default Gemmi reader. A confirmed SIGSEGV now triggers one retry with Coot's MMDB reader. Successful recovery is remembered by input hash for later chain requests. Timeouts, memory exhaustion, and ordinary parse errors are not retried. Both persistent-worker and one-shot execution rendered all three 6VXX protein chains with finite geometry and valid indices; subsequent requests used the remembered reader directly.

For the reported 8,660-atom, nine-chain CIF, the browser now completes loading with two successful mesh requests. The unchanged mesh JSON is 10,634,259 bytes and transfers as approximately 2,251,859 gzip bytes. Local uncached mesh delivery measured about 1.98 seconds for the original nine requests versus 0.80 seconds for batches of eight, before the final vector-conversion improvement. These are local mesh-processing/delivery measurements, not end-to-end hosted-page timings. Actual ribbon and bond-selection payloads were compared before and after the backend performance changes and were identical.

The analysis error `PDBe Arpeggio failed: 'type'` occurred before contact calculation because the exported CIF lacked `_chem_comp.type`. Roami now completes the metadata needed by Arpeggio in a temporary copy. Existing classifications and all coordinates, atom identifiers and chain labels are retained. Standard residue types and unambiguous polymer metadata inform missing types; unfamiliar components remain unknown. The original input file is never rewritten. The analysis cache version was updated so old reports are not reused.

The reported input completed analysis through the local HTTP API: Ca/Da returned 209 contacts in 1.90 seconds. Independent direct checks also passed for Ha/Hb (4,278 contacts) and the public 9GNQ A/B interface (2,303 contacts).

Run frontend checks from the site root with `npm run test:roami`. Run backend checks using the CHAPI environment with `python -m unittest discover -s ppi/api -p 'test_*.py'` and `PYTHONPATH=ppi`. Setting `ROAMI_CIF_REGRESSION_FILE` to a local CIF enables the additional file regression without adding its contents to the repository.

Final validation passed all 16 frontend tests and all 25 backend tests, with the optional local-file regressions enabled and no skips. A browser smoke test loaded the reported file with two successful mesh requests and no browser errors.

Publish the frontend changes and rebuild/restart the backend together. These changes are local pending the user's GitHub Desktop publishing workflow.
