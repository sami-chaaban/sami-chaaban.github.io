// Keep source identity tied to the actual coordinates, their format and API.
// A few entries accommodate original, ribbon-proxy and analysis-proxy sources.
export function createStructureTransport({ fetchImpl = (...args) => fetch(...args),
  getApiBase = () => '', maxEntries = 3, now = () => Date.now() } = {}) {
  const entries = [];
  const unsupportedApis = new Set();
  const abortError = () => new DOMException('Structure request was cancelled.', 'AbortError');
  const checkAbort = signal => { if (signal?.aborted) throw signal.reason || abortError(); };

  function waitFor(promise, signal) {
    checkAbort(signal);
    if (!signal) return promise;
    return new Promise((resolve, reject) => {
      const cancel = () => reject(signal.reason || abortError());
      signal.addEventListener('abort', cancel, { once: true });
      promise.then(resolve, reject).finally(() => signal.removeEventListener('abort', cancel));
    });
  }

  function deleteWhenUnused(entry) {
    if (!entry.retired || entry.users || entry.pending || !entry.structureId) return;
    const id = entry.structureId;
    entry.structureId = '';
    // Expiry is the backstop if the page closes or the network is unavailable.
    Promise.resolve().then(() => fetchImpl(`${entry.api}/structures/${encodeURIComponent(id)}`, {
      method: 'DELETE', keepalive: true,
    })).catch(() => {});
  }

  function retire(entry) {
    entry.retired = true;
    deleteWhenUnused(entry);
  }

  function clear() {
    for (const entry of entries.splice(0)) retire(entry);
  }

  function entryFor(api, source) {
    const found = entries.findIndex(entry => entry.api === api && entry.format === source.format && entry.text === source.text);
    if (found >= 0) {
      const [entry] = entries.splice(found, 1);
      entries.push(entry);
      return entry;
    }
    const entry = { ...source, api, structureId: '', expiresAt: 0, pending: null,
      error: null, users: 0, retired: false };
    entries.push(entry);
    while (entries.length > Math.max(1, maxEntries)) retire(entries.shift());
    return entry;
  }

  async function upload(entry) {
    if (entry.retired) throw abortError();
    if (unsupportedApis.has(entry.api)) return null;
    if (entry.error) throw entry.error;
    if (entry.structureId && now() < entry.expiresAt) return entry.structureId;
    if (entry.pending) return entry.pending;
    entry.pending = (async () => {
      const response = await fetchImpl(`${entry.api}/structures?format=${entry.format}`, {
        method: 'POST', headers: { 'Content-Type': 'text/plain; charset=utf-8', Accept: 'application/json' },
        body: entry.text,
      });
      if (response.status === 404 || response.status === 405) {
        unsupportedApis.add(entry.api);
        return null;
      }
      if (!response.ok) {
        let message = `Structure upload failed (${response.status}).`;
        try {
          const error = await response.json();
          const detail = error.detail;
          message = (typeof detail === 'string' ? detail : detail?.message) || error.message || message;
        } catch (_) { /* Keep the HTTP diagnostic for non-JSON errors. */ }
        const error = new Error(message);
        error.status = response.status;
        throw error;
      }
      const body = await response.json();
      if (!body || typeof body.structureId !== 'string' || !body.structureId.trim() ||
          body.format !== entry.format || !Number.isFinite(body.expiresIn) || body.expiresIn <= 0) {
        throw new Error('Structure upload returned an invalid handle.');
      }
      entry.structureId = body.structureId;
      entry.expiresAt = now() + body.expiresIn * 1000;
      return entry.structureId;
    })().catch(error => {
      // A failed upload is shared by this source's callers. Never turn failure
      // into dozens of full-coordinate uploads or an original-PDB refetch.
      error.structureUploadFailed = true;
      entry.error = error;
      throw error;
    }).finally(() => {
      entry.pending = null;
      deleteWhenUnused(entry);
    });
    return entry.pending;
  }

  async function post(path, payload, options = {}) {
    const signal = options.signal;
    checkAbort(signal);
    const api = String(getApiBase() || '').replace(/\/+$/, '');
    const source = typeof payload.pdbText === 'string' && payload.pdbText.length
      ? { text: payload.pdbText, format: 'pdb' }
      : typeof payload.mmcifText === 'string' && payload.mmcifText.length
        ? { text: payload.mmcifText, format: 'mmcif' } : null;
    const entry = source && !unsupportedApis.has(api) ? entryFor(api, source) : null;
    if (entry) entry.users += 1;
    try {
      for (let attempt = 0; attempt < 2; attempt += 1) {
        const structureId = entry ? await waitFor(upload(entry), signal) : null;
        checkAbort(signal);
        if (entry?.retired) throw abortError();
        const body = { ...payload };
        if (structureId) {
          delete body.pdbText; delete body.mmcifText;
          // pdbId is descriptive metadata. The backend resolves structureId
          // first, preserving displayed coordinates and report/export labels.
          body.structureId = structureId;
        }
        const headers = { 'Content-Type': 'application/json', ...options.headers };
        if (unsupportedApis.has(api)) {
          // Older backends keep the existing inline JSON mesh API.
          delete body.outputFormat;
          headers.Accept = 'application/json';
        }
        const response = await fetchImpl(`${api}${path}`, {
          method: 'POST', headers, signal, body: JSON.stringify(body),
        });
        if (response.status !== 410 || !structureId || attempt) return response;
        // Multiple requests may discover expiry simultaneously. Invalidate only
        // the handle this response used, then share the replacement upload.
        if (entry.structureId === structureId) {
          entry.structureId = '';
          entry.expiresAt = 0;
        }
        await response.body?.cancel();
      }
    } finally {
      if (entry) { entry.users -= 1; deleteWhenUnused(entry); }
    }
  }

  return { post, clear };
}
