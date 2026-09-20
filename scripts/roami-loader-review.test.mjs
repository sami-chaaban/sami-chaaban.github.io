import test from 'node:test';
import assert from 'node:assert/strict';
import vm from 'node:vm';
import { parseStructureRecords, serializeStructureRecords } from '../public/ppi/structure_parser.js';
import { harness, html } from './roami-test-harness.mjs';

const atom = overrides => ({chain:'A',seq:'42',atomName:'CA',resName:'ALA',element:'C',x:1,y:2,z:3,bFactor:20,occupancy:1,...overrides});

test('Coot proxy uses lossless CIF for more than 62 chains and unrepresentable PDB residue identifiers', () => {
  const {run,context} = harness();
  Object.assign(context,{parseStructureRecords,serializeStructureRecords,CHAPI_PROXY_CHAIN_SYMBOLS:'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'});
  context.state.structureFormat='mmcif';context.state.structureText=serializeStructureRecords([atom()]);
  context.state.atoms=Array.from({length:70},(_,i)=>atom({chain:`copy_${i}`,seq:'10001AB',x:i+0.123456789}));
  const proxy=run('buildChapiProxyFromExpandedAtoms(state.atoms)');
  assert.ok(proxy.mmcifText);assert.equal(proxy.pdbText,undefined);
  const parsed=parseStructureRecords(proxy.mmcifText);
  assert.equal(parsed.records.length,70);assert.equal(parsed.records[69].seq,'10001AB');
  assert.equal(parsed.records[69].x,context.state.atoms[69].x);
  assert.equal(proxy.chainForward.get('copy_69'),'copy_69');
  context.state.atoms=[atom({seq:'10000'})];
  assert.ok(run('buildChapiProxyFromExpandedAtoms(state.atoms).mmcifText'));
});

test('synthetic Coot payload cannot silently fall back to the original PDB entry', () => {
  const {run,context} = harness();
  context.state.chapiProxyMmcifText='data_displayed';context.state.chapiProxyPdbText='';
  assert.deepEqual(JSON.parse(JSON.stringify(run("getChapiStructurePayload({pdbId:'1abc',mmcifText:'data_original'})"))),{mmcifText:'data_displayed'});
  context.state.chapiProxyMmcifText='';context.state.chapiProxyPdbText='ATOM_DISPLAYED';
  assert.deepEqual(JSON.parse(JSON.stringify(run("getChapiStructurePayload({pdbId:'1abc'})"))),{pdbText:'ATOM_DISPLAYED'});
  assert.match(html,/const canTryRemoteIdFirst = Boolean\(normalizedPdbId && !state\.chapiProxyPdbText && !state\.chapiProxyMmcifText\)/);
});

test('cancelled parser workers cannot overwrite caches even without a load revision', async () => {
  const {run,context} = harness();const workers=[];
  context.parseStructureRecords=parseStructureRecords;
  context.Worker=class {constructor(){workers.push(this);}postMessage(){}terminate(){this.terminated=true;}};
  const first=run("prepareStructureCoordinateData('first','pdb')");
  const second=run("prepareStructureCoordinateData('second','pdb')");
  assert.equal(await first,null);
  workers[1].onmessage({data:{data:{records:[]}}});await second;
  workers[0].onmessage({data:{data:{records:[atom()]}}});
  assert.equal(run('getStructureCoordinateData.cache.text'),'second');
  assert.equal(workers[0].terminated,true);
});

test('failed older load stops reporting errors when a new load starts during restoration', async () => {
  const fn=html.match(/^      async function loadStructureFromFile[^\n]+\n[\s\S]*?^      \}/m)[0];
  let current=0;const history=[];let beginRestore,finishRestore;
  const restoring=new Promise(resolve=>{beginRestore=resolve;});
  const restored=new Promise(resolve=>{finishRestore=resolve;});
  const noop=()=>{};const state={viewer:{}};
  const ctx=vm.createContext({state,cancelActiveAnalysisRun:noop,nextLoadId:()=>++current,captureStructureLoadSnapshot:()=>({}),
    logDebug:value=>history.push(['debug',value]),setStructureDbInfo:noop,setViewerStructureTitle:noop,setStatus:noop,
    appendLoadDiagnostic:noop,inferFormat:()=> 'pdb',isCurrentLoad:id=>id===current,getLoadDiagnostics:()=>[],
    setStartupInteractionLocked:noop,setStartupRevealCurtainOpacity:noop,attachLoadDiagnostics:error=>error,
    restoreStructureLoadSnapshot:async()=>{beginRestore();await restored;},formatLoadDebugMessage:()=> 'error',
    setViewerPlaceholder:()=>history.push(['placeholder']),showRibbonLoadErrorStatus:()=>history.push(['error'])});
  vm.runInContext(fn,ctx);
  const old=ctx.loadStructureFromFile({name:'failed.pdb',text:async()=>{throw new Error('bad input');}});
  await restoring;current++;
  const previous=history.length;finishRestore();await old;
  assert.equal(history.length,previous);
});
