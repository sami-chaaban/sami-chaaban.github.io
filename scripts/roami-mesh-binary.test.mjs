import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import * as THREE from '../public/ppi/vendor/three/three.module.min.js';
import { decodeMeshBinary, MESH_BINARY_CONTENT_TYPE } from '../public/ppi/mesh_binary.js';

const source = {
  chainId: 'α-chain', vertexCount: 3, triangleCount: 1,
  positions: [0, 0, 0, 1, 0, 0, 0, 1, 0],
  normals: [0, 0, 1, 0, 0, 1, 0, 0, 1],
  colors: [1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1], indices: [0, 1, 2],
};

// Independent encoder fixture; production encoding lives in the Python bridge.
function fixture({ chained = true, mutate = () => {}, root = null } = {}) {
  const fields = ['positions', 'normals', 'colors', 'indices'];
  let offset = 0;
  const mesh = { ...source };
  for (const field of fields) {
    mesh[field] = { offset, count: source[field].length, type: field === 'indices' ? 'uint32' : 'float32' };
    offset += source[field].length * 4;
  }
  mutate(mesh);
  const metadata = root || (chained ? { meshType: 'chains', meshes: [mesh], chainIds: [mesh.chainId] } : mesh);
  const text = new TextEncoder().encode(JSON.stringify(metadata));
  const dataStart = Math.ceil((12 + text.length) / 4) * 4;
  const buffer = new ArrayBuffer(dataStart + offset);
  new Uint8Array(buffer).set(new TextEncoder().encode('ROAMIM01'));
  const view = new DataView(buffer);
  view.setUint32(8, text.length, true);
  new Uint8Array(buffer, 12, text.length).set(text);
  let cursor = dataStart;
  for (const field of fields) {
    for (const value of source[field]) {
      view[field === 'indices' ? 'setUint32' : 'setFloat32'](cursor, value, true);
      cursor += 4;
    }
  }
  return buffer;
}

for (const chained of [true, false]) {
  test(`decodes ${chained ? 'chain list' : 'single mesh'} as shared typed views preserving values`, () => {
    const buffer = fixture({ chained });
    const payload = decodeMeshBinary(buffer);
    const mesh = chained ? payload.meshes[0] : payload;
    assert.equal(mesh.chainId, source.chainId);
    assert.equal(mesh.vertexCount, source.vertexCount);
    for (const field of ['positions', 'normals', 'colors', 'indices']) {
      assert.ok(mesh[field] instanceof (field === 'indices' ? Uint32Array : Float32Array));
      assert.equal(mesh[field].buffer, buffer);
      assert.deepEqual([...mesh[field]], source[field]);
    }
  });
}

test('empty chain payload remains a valid empty result', () => {
  assert.deepEqual(decodeMeshBinary(fixture({ root: { meshType: 'chains', meshes: [], chainIds: [], status: 'empty' } })),
    { meshType: 'chains', meshes: [], chainIds: [], status: 'empty' });
});

test('rejects truncated header, metadata, and numeric buffers', () => {
  const buffer = fixture();
  for (const length of [0, 8, 11, 12, 24, buffer.byteLength - 1]) {
    assert.throws(() => decodeMeshBinary(buffer.slice(0, length)), /Truncated|descriptor/);
  }
  const badMagic = buffer.slice(0); new Uint8Array(badMagic)[0] = 0;
  assert.throws(() => decodeMeshBinary(badMagic), /magic/);
});

for (const [label, mutate] of [
  ['negative offset', mesh => { mesh.positions.offset = -4; }],
  ['unaligned offset', mesh => { mesh.positions.offset = 1; }],
  ['out-of-bounds offset', mesh => { mesh.positions.offset = 2 ** 32; }],
  ['negative count', mesh => { mesh.positions.count = -1; }],
  ['fractional count', mesh => { mesh.positions.count = 1.5; }],
  ['unsafe count', mesh => { mesh.positions.count = Number.MAX_SAFE_INTEGER + 1; }],
  ['unknown type', mesh => { mesh.positions.type = 'float64'; }],
  ['incorrect vertex count', mesh => { mesh.vertexCount = 20; }],
  ['incorrect triangle count', mesh => { mesh.triangleCount = 2; }],
  ['missing position array', mesh => { delete mesh.positions; }],
  ['missing triangle indices', mesh => { delete mesh.indices; }],
  ['missing triangle count', mesh => { delete mesh.triangleCount; }],
  ['missing normals descriptor', mesh => { delete mesh.normals; }],
  ['missing colors descriptor', mesh => { delete mesh.colors; }],
]) {
  test(`rejects ${label}`, () => assert.throws(() => decodeMeshBinary(fixture({ mutate })), /Invalid|Inconsistent/));
}

test('rejects triangle references outside the position array', () => {
  const buffer = fixture();
  new DataView(buffer).setUint32(buffer.byteLength - 4, 99, true);
  assert.throws(() => decodeMeshBinary(buffer), /out of range/);
});

const html = readFileSync(new URL('../public/ppi/index.html', import.meta.url), 'utf8');
test('production geometry builder accepts typed indices and keeps cached positions unchanged', () => {
  const mesh = decodeMeshBinary(fixture()).meshes[0];
  const definition = html.match(/^      function buildChapiGeometry[^\n]+\n[\s\S]*?^      \}/m)[0];
  const context = vm.createContext({ THREE, state: { largeStructureMode: true, renderSettings: { ribbonInflate: .12 } },
    LARGE_GEOMETRY_VERTEX_THRESHOLD: 10000, detectChapiColorScale: () => 1, mesh });
  vm.runInContext(definition, context);
  const geometry = vm.runInContext('buildChapiGeometry(mesh)', context);
  assert.ok(geometry.index instanceof THREE.BufferAttribute);
  assert.equal(geometry.index.array, mesh.indices);
  assert.equal(geometry.index.count, 3);
  assert.equal(geometry.getAttribute('position').count, 3);
  assert.deepEqual([...mesh.positions], source.positions);
  assert.ok(geometry.getAttribute('position').array[2] > 0);
  geometry.dispose();
});

test('production mesh loader negotiates binary and still accepts legacy JSON', async () => {
  const definition = html.match(/^      async function fetchChapiMeshJson[^\n]+\n[\s\S]*?^      \}/m)[0];
  for (const binary of [true, false]) {
    const calls = [];
    const context = vm.createContext({ API_BASE: 'https://backend.test', decodeMeshBinary, MESH_BINARY_CONTENT_TYPE,
      structureTransport: { post: async (...args) => {
        calls.push(args);
        return binary ? new Response(fixture(), { headers: { 'Content-Type': MESH_BINARY_CONTENT_TYPE } })
          : new Response(JSON.stringify(source), { headers: { 'Content-Type': 'application/json' } });
      } },
      appendLoadDiagnostic: () => {}, appendRibbonLoadErrorCode: () => {}, toErrorMessage: err => err.message,
      RIBBON_LOAD_ERROR_CODES: {},
    });
    vm.runInContext(definition, context);
    const result = await vm.runInContext("fetchChapiMeshJson({mmcifText:'coordinates'})", context);
    assert.ok(result);
    assert.equal(calls[0][1].outputFormat, 'binary');
    assert.equal(calls[0][2].headers.Accept, MESH_BINARY_CONTENT_TYPE);
    assert.equal(binary ? result.meshes[0].positions[3] : result.positions[3], 1);
  }
});
