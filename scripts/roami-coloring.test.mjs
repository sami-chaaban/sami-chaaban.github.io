import assert from 'node:assert/strict';
import test from 'node:test';
import { harness } from './roami-test-harness.mjs';

test('nearest atom sampling is exact across empty cells, cell boundaries and distant large structures', () => {
  const { run } = harness();
  assert.equal(run(`sampleStructureSpatialValueAtPoint(createStructureSpatialValueSampler([{x:-3.19,y:-3.19,z:-3.19,value:10},{x:6.4,y:3.19,z:3.19,value:90}],100),3.19,3.19,3.19)`), 90);
  const errors = run(`(() => {
    const samples=Array.from({length:2001},(_,i)=>({x:(i%19)*2.7,y:(i%37)*5.6,z:(i%17)*11.7,value:i}));
    const sampler=createStructureSpatialValueSampler(samples,0);
    const points=Array.from({length:40},(_,i)=>({x:i*11.33-100,y:i*7.81-120,z:i*21.73-150}));
    return points.filter(p=>{
      const d=s=>Math.hypot(s.x-p.x,s.y-p.y,s.z-p.z);
      const nearest=samples.reduce((a,b)=>d(a)<d(b)?a:b);
      const actual=samples[sampleStructureSpatialValueAtPoint(sampler,p.x,p.y,p.z)];
      return Math.abs(d(actual)-d(nearest))>1e-8;
    }).length;
  })()`);
  assert.equal(errors, 0);
});

test('missing metrics remain absent and use a missing-value color', () => {
  const { run } = harness();
  assert.equal(run('computeStructureBFactorStats([{bFactor:null},{bFactor:undefined},{bFactor:""}]).hasData'), false);
  assert.equal(run('computeStructureBFactorStats([{bFactor:null},{bFactor:0},{bFactor:100}]).mean'), 50);
  assert.equal(run(`setColorFromStructureMetricValue(null).getHexString()`), '8ea4be');
  assert.equal(run(`(() => {
    state.atoms=[{chain:'A',x:0,y:0,z:0,bFactor:null},{chain:'A',x:10,y:0,z:0,bFactor:100}];state.structureHasBFactor=true;
    return sampleChainBFactorAtPoint(getChainBFactorSampler('A'),0,0,0);
  })()`), null);
});

test('B-factor uses the observed range while pLDDT retains 0–100', () => {
  const { run } = harness();
  const result = run(`(() => {
    state.structureBFactorMin=20;state.structureBFactorMax=200;
    const experimental=[100,150,200].map(v=>setColorFromStructureMetricValue(v).getHexString());
    const range=getStructureMetricRange();state.structureMetricKind='plddt';
    const plddt=getStructureMetricRange();
    return {experimental,range,plddt};
  })()`);
  assert.equal(new Set(result.experimental).size, 3);
  assert.deepEqual(JSON.parse(JSON.stringify(result.range)), {min:20,max:200});
  assert.deepEqual(JSON.parse(JSON.stringify(result.plddt)), {min:0,max:100});
});

test('ribbon color modes reuse vertex associations and retain independent color buffers', () => {
  const { run } = harness();
  const result = run(`(() => {
    const spatial=createStructureSpatialValueSampler([{x:0,y:0,z:0},{x:10,y:0,z:0}],0);
    const metric={spatial,samples:[{value:20},{value:90}]},rainbow={spatial,samples:[{value:0},{value:1}]};
    const mesh=new THREE.Mesh(new THREE.BufferGeometry(),new THREE.MeshBasicMaterial());
    mesh.geometry.setAttribute('position',new THREE.Float32BufferAttribute([0,0,0,10,0,0],3));
    applyBFactorColorsToMeshGeometry(mesh,'A',metric);
    const mapping=structureVertexColorCache.get(mesh.geometry).indices,colors=mesh.geometry.getAttribute('color');
    applyRainbowColorsToMeshGeometry(mesh,'A',rainbow);
    const reused=mapping===structureVertexColorCache.get(mesh.geometry).indices;
    applyBFactorColorsToMeshGeometry(mesh,'A',metric);
    const restored=colors===mesh.geometry.getAttribute('color');
    mesh.geometry.getAttribute('position').setX(0,9);mesh.geometry.getAttribute('position').needsUpdate=true;
    applyBFactorColorsToMeshGeometry(mesh,'A',metric);
    return [reused,restored,mapping!==structureVertexColorCache.get(mesh.geometry).indices,structureVertexColorCache.get(mesh.geometry).indices[0]];
  })()`);
  assert.deepEqual(Array.from(result), [true, true, true, 1]);
});

test('atom mode recolors instances in place and preserves element colors', () => {
  const { run } = harness();
  const result = run(`(() => {
    state.atoms=[{chain:'A',element:'C',bFactor:20,x:0,y:0,z:0},{chain:'A',element:'O',bFactor:90,x:1,y:0,z:0}];
    state.structureHasBFactor=true;state.structureColorMode='bfactor';
    const mesh=new THREE.InstancedMesh(new THREE.SphereGeometry(1),new THREE.MeshBasicMaterial(),2);
    mesh.userData.structureColorAtoms=state.atoms;
    const geometry=mesh.geometry,matrix=mesh.instanceMatrix;
    applyStructureInstanceColors(mesh,'A',new THREE.Color('#123456'));
    const metric=mesh.instanceColor;const oxygen=new THREE.Color();mesh.getColorAt(1,oxygen);
    state.structureColorMode='chains';applyStructureInstanceColors(mesh,'A',new THREE.Color('#123456'));
    const chain=new THREE.Color();mesh.getColorAt(0,chain);
    state.structureColorMode='bfactor';applyStructureInstanceColors(mesh,'A',new THREE.Color('#123456'));
    return [geometry===mesh.geometry,matrix===mesh.instanceMatrix,metric===mesh.instanceColor,oxygen.getHexString(),chain.getHexString()];
  })()`);
  assert.deepEqual(Array.from(result), [true, true, true, 'e53935', '123456']);
});
