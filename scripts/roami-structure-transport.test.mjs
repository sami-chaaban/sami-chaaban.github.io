import test from 'node:test';
import assert from 'node:assert/strict';
import { createStructureTransport } from '../public/ppi/structure_transport.js';

const json = (data, status = 200) => new Response(JSON.stringify(data), { status });
const handle = (id, format = 'mmcif') => json({ structureId: id, format, expiresIn: 60 });
const flush = () => new Promise(resolve => setImmediate(resolve));

function harness(respond, options = {}) {
  const calls = [];
  let api = 'https://example.test/api';
  const client = createStructureTransport({
    ...options, getApiBase: () => api,
    fetchImpl: async (url, init) => {
      const call = { url, ...init };
      calls.push(call);
      return respond(call, calls);
    },
  });
  return { ...client, calls, setApi: value => { api = value; } };
}

test('one raw upload is shared by concurrent mesh and analysis requests using the actual source', async () => {
  const client = harness(({ url }) => url.includes('/structures?') ? handle('same-source') : json({ ok: true }));
  const source = { mmcifText: 'data_actual\nATOM 1 C A 0 1 2\n', pdbId: '7z8g' };
  await Promise.all([
    client.post('/chapi-mesh', { ...source, chainIds: ['A'], outputFormat: 'binary' }),
    client.post('/analyze', { ...source, chainA: 'A', chainB: 'B' }),
  ]);
  assert.equal(client.calls.length, 3);
  assert.equal(client.calls[0].body, source.mmcifText);
  assert.match(client.calls[0].headers['Content-Type'], /^text\/plain/);
  for (const call of client.calls.slice(1)) {
    const payload = JSON.parse(call.body);
    assert.equal(payload.structureId, 'same-source');
    assert.equal(payload.mmcifText, undefined);
    assert.equal(payload.pdbId, '7z8g'); // Label metadata survives; structureId identifies the coordinates.
  }
});

test('identity includes changed coordinates, format, and API, never merely PDB ID', async () => {
  let uploads = 0;
  const client = harness(({ url }) => url.includes('/structures?')
    ? handle(`source-${++uploads}`, new URL(url).searchParams.get('format')) : json({}));
  const request = source => client.post('/chapi-mesh', { ...source, pdbId: '7z8g' });
  await request({ mmcifText: 'original' });
  await request({ mmcifText: 'original' });
  await request({ mmcifText: 'edited' });
  await request({ pdbText: 'edited' });
  client.setApi('https://second.test');
  await request({ mmcifText: 'edited' });
  assert.equal(uploads, 4);
  assert.deepEqual(client.calls.filter(c => c.url.endsWith('/chapi-mesh')).map(c => JSON.parse(c.body).structureId),
    ['source-1', 'source-1', 'source-2', 'source-3', 'source-4']);
});

test('simultaneous expired-handle responses share exactly one replacement upload', async () => {
  let uploads = 0;
  const client = harness(({ url, body }) => {
    if (url.includes('/structures?')) return handle(`id-${++uploads}`);
    return JSON.parse(body).structureId === 'id-1' ? json({ detail: 'Expired' }, 410) : json({ ok: true });
  });
  const responses = await Promise.all([
    client.post('/chapi-mesh', { mmcifText: 'coordinates', chainIds: ['A'] }),
    client.post('/chapi-mesh', { mmcifText: 'coordinates', chainIds: ['B'] }),
  ]);
  assert.equal(uploads, 2);
  assert.deepEqual(responses.map(r => r.status), [200, 200]);
  const payloads = client.calls.filter(c => c.url.endsWith('/chapi-mesh')).map(c => JSON.parse(c.body));
  assert.deepEqual(payloads.map(p => p.structureId), ['id-1', 'id-1', 'id-2', 'id-2']);
  assert.ok(payloads.every(p => !p.mmcifText && !p.pdbId));
});

test('a second 410 is surfaced without an upload or request loop', async () => {
  let uploads = 0;
  const client = harness(({ url }) => url.includes('/structures?') ? handle(`id-${++uploads}`) : json({}, 410));
  assert.equal((await client.post('/analyze', { mmcifText: 'coordinates' })).status, 410);
  assert.equal(uploads, 2);
  assert.equal(client.calls.filter(c => c.url.endsWith('/analyze')).length, 2);
});

test('the advertised expiry renews the source before sending another request', async () => {
  let time = 0, uploads = 0;
  const client = harness(({ url }) => url.includes('/structures?') ? handle(`id-${++uploads}`) : json({}), { now: () => time });
  await client.post('/chapi-mesh', { mmcifText: 'coordinates' });
  time = 60_000;
  await client.post('/chapi-mesh', { mmcifText: 'coordinates' });
  assert.equal(uploads, 2);
});

for (const status of [404, 405]) {
  test(`older backend ${status} is probed once then uses original inline JSON`, async () => {
    const client = harness(({ url }) => url.includes('/structures?') ? json({}, status) : json({ ok: true }));
    const payload = { mmcifText: 'user-edited-coordinates', outputFormat: 'binary' };
    await client.post('/chapi-mesh', payload, { headers: { Accept: 'application/vnd.roami.mesh' } });
    await client.post('/chapi-mesh', payload);
    assert.equal(client.calls.filter(c => c.url.includes('/structures?')).length, 1);
    for (const call of client.calls.filter(c => c.url.endsWith('/chapi-mesh'))) {
      assert.deepEqual(JSON.parse(call.body), { mmcifText: payload.mmcifText });
      assert.equal(call.headers.Accept, 'application/json');
    }
    client.setApi('https://upgraded.test');
    await client.post('/chapi-mesh', payload);
    assert.equal(client.calls.filter(c => c.url.includes('/structures?')).length, 2);
  });
}

for (const failure of [() => json({ detail: 'Upload too large' }, 413), () => { throw new Error('Network down'); }]) {
  test(`upload failure is shared and never turns into repeated inline fallbacks (${failure.toString()})`, async () => {
    const client = harness(failure);
    const payload = { mmcifText: 'coordinates' };
    const results = await Promise.allSettled([client.post('/chapi-mesh', payload), client.post('/analyze', payload)]);
    assert.ok(results.every(r => r.status === 'rejected' && r.reason.structureUploadFailed));
    await assert.rejects(client.post('/chapi-mesh', payload));
    assert.equal(client.calls.length, 1);
    client.clear();
    await assert.rejects(client.post('/chapi-mesh', payload));
    assert.equal(client.calls.length, 2);
  });
}

test('aborting one upload waiter leaves the shared upload usable by another', async () => {
  let finishUpload;
  const client = harness(({ url }) => url.includes('/structures?')
    ? new Promise(resolve => { finishUpload = () => resolve(handle('id')); }) : json({}));
  const controller = new AbortController();
  const cancelled = client.post('/analyze', { mmcifText: 'coordinates' }, { signal: controller.signal });
  const surviving = client.post('/chapi-mesh', { mmcifText: 'coordinates' });
  controller.abort();
  await assert.rejects(cancelled, { name: 'AbortError' });
  finishUpload();
  assert.equal((await surviving).status, 200);
  assert.equal(client.calls.length, 2);
});

test('unload releases a handle after in-flight requests finish, and uploads that finish after unload are deleted', async () => {
  let finishRequest;
  const client = harness(({ url, method }) => url.includes('/structures?') ? handle('id')
    : method === 'DELETE' ? json({}) : new Promise(resolve => { finishRequest = resolve; }));
  const pending = client.post('/analyze', { mmcifText: 'coordinates' });
  await flush();
  client.clear();
  assert.equal(client.calls.filter(c => c.method === 'DELETE').length, 0);
  finishRequest(json({})); await pending; await flush();
  assert.equal(client.calls.filter(c => c.method === 'DELETE').length, 1);

  let finishUpload;
  const duringUpload = harness(({ url }) => url.includes('/structures?')
    ? new Promise(resolve => { finishUpload = resolve; }) : json({}));
  const obsolete = duringUpload.post('/chapi-mesh', { mmcifText: 'old' });
  duringUpload.clear(); finishUpload(handle('old-id'));
  await assert.rejects(obsolete, { name: 'AbortError' }); await flush();
  assert.equal(duringUpload.calls.at(-1).method, 'DELETE');
  assert.match(duringUpload.calls.at(-1).url, /old-id$/);
});

test('the bounded source cache retires an unused old handle', async () => {
  let uploads = 0;
  const client = harness(({ url }) => url.includes('/structures?') ? handle(`id-${++uploads}`) : json({}), { maxEntries: 1 });
  await client.post('/chapi-mesh', { mmcifText: 'one' });
  await client.post('/chapi-mesh', { mmcifText: 'two' }); await flush();
  assert.equal(client.calls.filter(c => c.method === 'DELETE').length, 1);
  assert.match(client.calls.find(c => c.method === 'DELETE').url, /id-1$/);
});
