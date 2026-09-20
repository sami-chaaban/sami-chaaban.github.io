import assert from 'node:assert/strict';
import test from 'node:test';
import { harness, html } from './roami-test-harness.mjs';

test('PAE preserves direction and raw large values while missing entries stay missing', () => {
  const { run } = harness();
  const result = run(`(() => {
    const normalized=normalizeAlphaFoldPaeMatrix([[0,null,43],[12,0,''],[7,18,0]]);
    state.alphaFoldPae={matrix:normalized.matrix,axis:[{chain:'A',seq:'1'},{chain:'B',seq:'2'},{chain:'B',seq:'3'}]};
    const missing=resolveAlphaFoldPaePoint(0,1),forward=resolveAlphaFoldPaePoint(0,2),reverse=resolveAlphaFoldPaePoint(2,0);
    return [missing.paeAngstrom,forward.paeAngstrom,reverse.paeAngstrom,forward.alignmentFrameLabel,forward.assessedPositionLabel,
      getAlphaFoldPaeDisplayMax(),normalized.observedMax,formatPaePanelValue(null),formatPaePanelValue(0),
      [null,undefined,'',false,[],{}].every(v=>Number.isNaN(parseAlphaFoldPaeNumber(v)))];
  })()`);
  assert.deepEqual(Array.from(result), [null,43,7,'A1','B3',30,43,'NA','0.00',true]);
});

test('PAE heatmap distinguishes missing from zero and saturates only the display color', () => {
  const { run, context } = harness();
  context.ImageData = class {
    constructor(width, height) { this.data = new Uint8ClampedArray(width * height * 4); }
  };
  run(html.match(/^      const PAE_COLORMAP_START_RGB = [\s\S]*?\n      \}\);/m)[0]);
  const result = run(`(() => {
    const pae={matrix:[[0,null],[12,43]]};
    const pixels=createPaePanelHeatmapImageData(pae,2,2).data;
    return [Array.from(pixels.slice(0,4)),Array.from(pixels.slice(4,8)),Array.from(pixels.slice(12,16)),pae.matrix[1][1]];
  })()`);
  assert.notDeepEqual(Array.from(result[0]), Array.from(result[1]));
  assert.deepEqual(Array.from(result[1]), [17,24,34,255]);
  assert.deepEqual(Array.from(result[2]), [255,255,255,255]);
  assert.equal(result[3],43);
});

test('PAE focus metadata separates error from physical separation and never invents missing confidence', () => {
  const { run } = harness();
  const result = run(`(() => {
    buildPaeFocusContactResidue=key=>({chain:key[0],seq:'1',resName:'ALA',atom:'CA'});
    getResidueSidechainCenterOfMass=key=>new THREE.Vector3(key==='A:1'?0:12,0,0);
    const known=buildSyntheticPaePairContact('A:1','B:1',{paeAngstrom:2.5});
    const missing=buildSyntheticPaePairContact('A:1','B:1',{paeAngstrom:null});
    return [known.distance,known.paeAngstrom,known.paeAlignmentFrame,known.paeAssessedPosition,
      formatContactDistanceLabel(known),missing.distance,missing.paeAngstrom,formatContactDistanceLabel(missing)];
  })()`);
  assert.deepEqual(Array.from(result), [12,2.5,'A:1','B:1','PAE 2.50 Å',12,null,'PAE unavailable']);
});

test('native pLDDT imports retain the 0–100 scale and original experimental values', () => {
  const { run } = harness();
  const result = run(`(() => {
    state.atoms=[{bFactor:20},{bFactor:30},{bFactor:40}];
    applyLocalPlddtToAtomValues([0,1,null]);
    const first=state.atoms.map(a=>[a.bFactor,a.originalBFactor,a.plddt]);
    applyLocalPlddtToAtomValues([0.5,null,90]);
    const second=state.atoms.map(a=>[a.bFactor,a.originalBFactor,a.plddt]);
    const invalid=normalizeLocalPredictionNumericArray([null,'',false,[],{},-1,101]);
    return {first,second,invalid};
  })()`);
  assert.deepEqual(JSON.parse(JSON.stringify(result)), {
    first:[[0,20,0],[1,30,1],[null,40,null]],
    second:[[0.5,20,0.5],[null,30,null],[90,40,90]],
    invalid:null,
  });
});

test('residue pLDDT imports clear scores for missing and unmapped atoms without losing B-factors', () => {
  const { run } = harness();
  const result = run(`(() => {
    state.atoms=[{bFactor:20},{bFactor:30},{bFactor:40},{bFactor:50}];
    state.residueAtoms=new Map([['A:1',[0,1]],['A:2',[2]],['B:1',[3]]]);
    const applied=applyLocalPlddtToResidueValues([85,null],[{resKey:'A:1'},{resKey:'A:2'}]);
    const before=state.atoms.map(a=>[a.bFactor,a.originalBFactor]);
    const failed=applyLocalPlddtToResidueValues([70],[{resKey:'Z:99'}]);
    return {applied:applied.appliedAtomCount,before,failed,after:state.atoms.map(a=>[a.bFactor,a.originalBFactor])};
  })()`);
  const normalized = JSON.parse(JSON.stringify(result));
  assert.equal(normalized.applied,2);
  assert.deepEqual(normalized.before,[[85,20],[85,30],[null,40],[null,50]]);
  assert.equal(normalized.failed,null);
  assert.deepEqual(normalized.after,normalized.before);
});

test('metric educational panels follow experimental method rather than any electron keyword', () => {
  const { run } = harness();
  const result = run(`(() => {
    state.structureColorMode='bfactor';state.structureHasBFactor=true;
    const methods=['X-RAY DIFFRACTION','ELECTRON MICROSCOPY','CRYO-EM','ELECTRON CRYSTALLOGRAPHY','SOLUTION NMR','X-RAY DIFFRACTION; ELECTRON MICROSCOPY',''];
    const modes=methods.map(method=>{state.structureExperimentMethod=method;return getStructureMetricLegendInfoMode();});
    state.structureMetricKind='plddt';modes.push(getStructureMetricLegendInfoMode());
    return modes;
  })()`);
  assert.deepEqual(Array.from(result), ['bfactor','bfactor-em','bfactor-em','','','','','plddt']);
});
