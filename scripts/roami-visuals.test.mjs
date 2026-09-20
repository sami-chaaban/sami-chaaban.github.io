import assert from 'node:assert/strict';
import test from 'node:test';
import { harness } from './roami-test-harness.mjs';
import { buildCartoonGroup } from '../public/ppi/ribbon_geometry.js';

function geometryHarness() {
  const h = harness();
  h.run(`
    function fixtureAtom(atomName,resName,chain,seq,x,y,z) { return {atomName,element:atomName.startsWith('N')?'N':'C',resName,chain,seq,resKey:chain+':'+seq,x,y,z}; }
    function fixtureRing(names,resName,chain,seq,cx,cy,z,radius=1.4) { return names.map((name,i)=>fixtureAtom(name,resName,chain,seq,cx+Math.cos(2*Math.PI*i/names.length)*radius,cy+Math.sin(2*Math.PI*i/names.length)*radius,z)); }
    function fixtureAtoms(atoms) {
      state.atoms=atoms;state.residueAtoms=new Map();state.loadId++;
      atoms.forEach((a,i)=>{const list=state.residueAtoms.get(a.resKey)||[];list.push(i);state.residueAtoms.set(a.resKey,list);});
    }
  `);
  return h;
}

test('protonated histidine can be the cation opposite an aromatic ring; neutral HIS cannot', () => {
  const { run } = geometryHarness();
  const result = run(`(() => {
    fixtureAtoms([...fixtureRing(['CG','CD1','CE1','CZ','CE2','CD2'],'PHE','A','1',0,0,0), ...fixtureRing(['CG','ND1','CE1','NE2','CD2'],'HIP','B','2',0,0,3.6,1.2)]);
    const contact={source:'arpeggio',type:'pi_cation',residueA:{chain:'A',seq:'1',resName:'PHE',atom:'CZ'},residueB:{chain:'B',seq:'2',resName:'HIP',atom:'ND1'}};
    const protonated=resolvePiCationParticipants(contact);
    const profile=resolveContactSideChemicalProfile(contact,'B');
    const charged={...contact,chargeB:1,residueB:{...contact.residueB,resName:'HIS'}};
    const neutral={...contact,residueB:{...contact.residueB,resName:'HIS'}};
    return [protonated?.cationSide,profile.canAccept,!!resolvePiCationParticipants(charged),!!resolvePiCationParticipants(neutral)];
  })()`);
  assert.deepEqual(Array.from(result), ['B', false, true, false]);
});

test('fused-ring visualization matches explicit membership and reported centroid distance', () => {
  const { run } = geometryHarness();
  const result = run(`(() => {
    const trp=fixtureRing(['CE2','CD2','CE3','CZ3','CH2','CZ2'],'TRP','T','10',0,0,0);
    let x=trp[1].x,y=trp[1].y,dx=trp[1].x-trp[0].x,dy=trp[1].y-trp[0].y;
    for(const name of ['CG','CD1','NE1']) { const a=-2*Math.PI/5,ndx=dx*Math.cos(a)-dy*Math.sin(a),ndy=dx*Math.sin(a)+dy*Math.cos(a);dx=ndx;dy=ndy;x+=dx;y+=dy;trp.push(fixtureAtom(name,'TRP','T','10',x,y,0)); }
    const names=['CG','CD1','NE1','CE2','CD2'];
    const ring=computeRingDescriptorFromAtomSet(trp.filter(a=>names.includes(a.atomName)));
    const partner=fixtureRing(['CG','CD1','CE1','CZ','CE2','CD2'],'PHE','P','20',ring.centroid.x,ring.centroid.y,3.6);
    fixtureAtoms([...trp,...partner]);
    const contact={source:'arpeggio',type:'pi_pi',distance:3.6,ringKeyA:'T:10|'+names.join(','),ring:{centroid_distance:3.6},residueA:{chain:'T',seq:'10',resName:'TRP',atom:'CD2'},residueB:{chain:'P',seq:'20',resName:'PHE',atom:'CZ'}};
    const explicit=resolvePiPiParticipants(contact);
    const backend={...contact,ringKeyA:'T:10:TRP:indole:opaque',ringAtomNamesA:names};
    const exported=resolvePiPiParticipants(backend);
    const allNames=trp.map(a=>a.atomName);
    const system=ringDescriptorsForContactSide({...contact,ringAtomNamesA:allNames},'A',contact.residueA)[0];
    return [preferredRingDescriptorsForInteraction(contact.residueA).map(r=>r.atomNames.length).sort().join(','),explicit?.ringA.centroid.distanceTo(ring.centroid),explicit?.distance,exported?.ringA.centroid.distanceTo(ring.centroid),system.atomNames.length];
  })()`);
  assert.equal(result[0], '5,6');
  assert.ok(result[1] < 1e-10);
  assert.ok(Math.abs(result[2] - 3.6) < 1e-10);
  assert.ok(result[3] < 1e-10);
  assert.equal(result[4], 9, 'an exported full-system descriptor must remain a full-system descriptor');
});

test('soft clipping composes with Basic and Physical shaders without conditional worldPosition', () => {
  const { run } = harness();
  const shaders = run(`['basic','physical'].map(kind=>{
    const material=kind==='basic'?new THREE.MeshBasicMaterial():new THREE.MeshPhysicalMaterial();
    const shader={vertexShader:THREE.ShaderLib[kind].vertexShader,fragmentShader:THREE.ShaderLib[kind].fragmentShader,uniforms:{}};
    ensureSoftClipFadeOnMaterial(material);material.onBeforeCompile(shader);
    return shader.vertexShader;
  })`);
  for (const shader of shaders) {
    assert.ok(!shader.includes('vSoftClipWorldPos = worldPosition.xyz'));
    assert.ok(shader.includes('vec4(transformed, 1.0)'));
    assert.ok(shader.includes('roamiWorldPosition = batchingMatrix * roamiWorldPosition'));
    assert.ok(shader.includes('roamiWorldPosition = instanceMatrix * roamiWorldPosition'));
    assert.ok(shader.includes('vSoftClipWorldPos = (modelMatrix * roamiWorldPosition).xyz'));
  }
});

test('beta-strand sweep includes all four sides of the cross-section', () => {
  const samples = {pos:[0,0,0,0,0,1,0,0,2,0,0,3], n:[1,0,0,1,0,0,1,0,0,1,0,0],b:[0,1,0,0,1,0,0,1,0,0,1,0],w:[1,1,1,1],h:[0.3,0.3,0.3,0.3],ss:['E','E','E','E']};
  const result = buildCartoonGroup({chains:[{id:'A',segments:[{samples}]}]});
  const core = result.userData.chainData.get('A').coreMeshes[0];
  assert.equal(core.geometry.index.count, 3 * 4 * 2 * 3, '3 longitudinal spans × 4 sides × 2 triangles × 3 indices');
});

test('material cloning never serializes GPU shader graphs or copies installed-hook flags', () => {
  const { run } = harness();
  const result = run(`(() => {
    const material=new THREE.MeshBasicMaterial();ensureSoftClipFadeOnMaterial(material);
    setMaterialRuntimeValue(material,'shader',{texture:{toJSON(){throw new Error('GPU graph serialized');}}});
    const copy=material.clone();
    const staleFlag=Boolean(copy.userData.softClipInstalled);
    ensureSoftClipFadeOnMaterial(copy);
    return [staleFlag,Boolean(copy.userData.softClipInstalled),Object.hasOwn(copy.userData,'shader')];
  })()`);
  assert.deepEqual(Array.from(result), [false, true, false]);
});
