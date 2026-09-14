import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';

const html = readFileSync(new URL('../public/ppi/index.html', import.meta.url), 'utf8');
const names = new Set(['buildChapiRequestOptions', 'fetchChapiMeshJson',
  'normalizeChapiChunkChainIds', 'extractChapiChunkMeshes', 'fetchChunkedChapiRibbonMeshData', 'parseApiErrorPayload']);
const functions = [...html.matchAll(/^      (?:async )?function ([^(]+)[^\n]*\n[\s\S]*?^      \}/gm)]
  .filter(match => names.has(match[1])).map(match => match[0]).join('\n');
const constants = ['CHAPI_RIBBON_CHUNK_MIN_CHAIN_COUNT', 'CHAPI_RIBBON_BATCH_MAX_CHAIN_COUNT']
  .map(name => html.match(new RegExp(`^      const ${name} = .+;$`, 'm'))[0]).join('\n');
const chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'];
const meshFor = chainId => ({ chainId, vertexCount: 3, positions: [0, 0, 0, 1, 0, 0, 0, 1, 0], indices: [0, 1, 2] });
const success = payload => new Response(JSON.stringify({
  meshType: 'chains', meshes: [...payload.chainIds].reverse().map(meshFor),
}), { status: 200 });

function harness(respond = success) {
  const requests = [];
  const context = vm.createContext({
    Response, console, API_BASE: 'http://local.test', state: { largeStructureMode: false },
    collectChapiRequestChainIds: () => chains,
    getActiveRenderPreset: () => ({ theme: 'dark' }), getSceneTheme: () => ({ id: 'dark' }),
    toErrorMessage: error => String(error?.message || error),
    appendLoadDiagnostic: (entries, text) => entries?.push(text),
    appendRibbonLoadErrorCode: () => {}, appendRibbonBackendErrorCode: () => {},
    RIBBON_LOAD_ERROR_CODES: { REQUEST_FAILED: 'request', BACKEND_HTTP: 'http', INVALID_PAYLOAD: 'json', INCOMPATIBLE_MESH: 'incomplete' },
    fetch: async (_url, options) => {
      const payload = JSON.parse(options.body);
      requests.push(payload);
      return respond(payload, requests.length);
    },
  });
  vm.runInContext(`${constants}\nlet chapiRibbonBatchLimit = CHAPI_RIBBON_BATCH_MAX_CHAIN_COUNT;\n${functions}`, context);
  const load = async (source = { mmcifText: 'data_synthetic\n#\n' }, options = {}) => {
    context.source = source;
    context.options = { chainIds: chains, ...options };
    return vm.runInContext('fetchChunkedChapiRibbonMeshData(source, options)', context);
  };
  return { load, requests, context };
}

for (const source of [{ mmcifText: 'data_synthetic\n#\n' }, { pdbId: '6VXX' }]) {
  test(`batches ${Object.keys(source)[0]} loads without changing chain order or mesh settings`, async () => {
    const { load, requests } = harness();
    const result = await load(source);
    assert.deepEqual(requests.map(request => request.chainIds.length), [8, 1]);
    assert.deepEqual(Array.from(result.meshes, mesh => mesh.chainId), chains);
    assert.deepEqual(Array.from(result.expectedChainIds), chains);
    for (const request of requests) {
      assert.equal(request.secondaryStructureUsage, 2);
      assert.equal(request.style, 'Ribbon');
      assert.equal(request.colourScheme, 'colorRampChainsScheme');
      for (const key of Object.keys(source)) assert.equal(request[key], source[key]);
    }
  });
}

test('learns a low-memory backend limit and retries each chain exactly once', async () => {
  const { load, requests } = harness(payload => payload.chainIds.length > 1
    ? new Response(JSON.stringify({ detail: { errorCode: 'CHAPI-LOWMEM-002', message: 'One chain per request' } }), { status: 413 })
    : success(payload));
  assert.equal((await load()).complete, true);
  assert.deepEqual(requests.map(request => request.chainIds.length), [8, ...chains.map(() => 1)]);
  requests.length = 0;
  assert.equal((await load()).complete, true);
  assert.deepEqual(requests.map(request => request.chainIds.length), chains.map(() => 1));
});

test('retains one-chain requests for the large-structure memory safeguard', async () => {
  const { load, requests, context } = harness();
  context.state.largeStructureMode = true;
  assert.equal((await load()).complete, true);
  assert.deepEqual(requests.map(request => request.chainIds.length), chains.map(() => 1));
});

test('refuses incomplete batches without rendering a partial model', async () => {
  const { load, requests } = harness(() => new Response(JSON.stringify({ meshes: [meshFor('A')] })));
  assert.equal(await load(), null);
  assert.equal(requests.length, 1);
});

test('stops after a network failure instead of re-uploading the file for every chain', async () => {
  const { load, requests } = harness(() => { throw new Error('Backend unavailable'); });
  assert.equal(await load(), null);
  assert.equal(requests.length, 1);
});

test('does not retry unrelated HTTP 413 failures as single-chain requests', async () => {
  const { load, requests } = harness(() => new Response('Upload too large', { status: 413 }));
  assert.equal(await load(), null);
  assert.equal(requests.length, 1);
});

test('a newer structure load cancels the remaining batches', async () => {
  let stale = false;
  const { load, requests } = harness(payload => { stale = true; return success(payload); });
  assert.equal(await load(undefined, { isStaleLoad: () => stale }), null);
  assert.equal(requests.length, 1);
});
