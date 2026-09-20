import test from 'node:test';
import assert from 'node:assert/strict';
import vm from 'node:vm';
import { parseCif, parseStructureRecords, serializeStructureRecords, selectConformers } from '../public/ppi/structure_parser.js';
import { harness, html } from './roami-test-harness.mjs';

const atom = (overrides={}) => ({chain:'A',authChain:'A',labelChain:'A',modelToken:'1',seq:'42',atomName:'N',resName:'LIG',element:'N',x:0,y:1,z:2,bFactor:null,occupancy:1,altLoc:'',...overrides});
const snapshot = value => JSON.parse(JSON.stringify(value));
function parserHarness() { const h = harness(); Object.assign(h.context, {parseStructureRecords,serializeStructureRecords}); return h; }

test('CIF rows are independent of physical lines and retain quoted/multiline metadata', () => {
  const source = serializeStructureRecords([atom(),atom({atomName:'O',element:'O',seq:'42A'})]);
  const tagsEnd = source.indexOf('_atom_site.pdbx_PDB_model_num')+'_atom_site.pdbx_PDB_model_num'.length;
  const tokens = source.slice(tagsEnd).trim().replace(/\n#$/, '').split(/\s+/);
  const packed = source.slice(0,tagsEnd)+'\n'+tokens.map((t,i)=>t+(i%7===0?'\n':' ')).join('')+'\n#\n';
  assert.deepEqual(parseStructureRecords(packed).records,parseStructureRecords(source).records);
  const metadata = parseCif('data_x\n_entity.pdbx_description\n;A long\nquoted description\n;\nloop_\n_chem_comp.id\n_chem_comp.name\nLIG "iron ligand" XYZ \'a "quoted" label\'\n#\n');
  assert.equal(metadata.items.get('_entity.pdbx_description'),'A long\nquoted description');
  assert.deepEqual(metadata.tables[0].rows,[['LIG','iron ligand'],['XYZ','a "quoted" label']]);
});

test('conformer selection retains B-only residues, shared atoms, and one coherent occupancy winner', () => {
  const records=[atom(),atom({atomName:'O',altLoc:'A',occupancy:0.3,x:1}),atom({atomName:'O',altLoc:'B',occupancy:0.7,x:2}),atom({atomName:'C',altLoc:'B',occupancy:0.7}),atom({seq:'43',altLoc:'B'}),atom({seq:'44',altLoc:'B',occupancy:0.5}),atom({seq:'44',altLoc:'A',occupancy:0.5})];
  assert.deepEqual(selectConformers(records).map(a=>[a.seq,a.atomName,a.altLoc,a.x]),[['42','N','',0],['42','O','B',2],['42','C','B',0],['43','N','B',0],['44','N','A',0]]);
});

test('displayed assembly payload retains transformed coordinates, insertion codes, elements and >62 chains', () => {
  const {run,context}=parserHarness();
  context.state.atoms=Array.from({length:70},(_,i)=>atom({chain:`copy_${i}`,seq:'42A',x:i*10.123456,element:i%2?'CL':'FE',atomName:i%2?'CL1':'FE1'}));
  Object.assign(context.state,{assemblyExpanded:true,structureFormat:'mmcif',structureText:serializeStructureRecords([atom()]),backendChainAliases:new Map([['copy_69','A']]),structureChainAliases:new Map([['copy_69','A']])});
  const result=run("createAnalyzePayload('copy_0','copy_69')");
  assert.equal(result.payload.chainA,'copy_0');assert.equal(result.payload.chainB,'copy_69');
  const parsed=parseStructureRecords(result.payload.mmcifText);
  assert.equal(new Set(parsed.records.map(a=>a.chain)).size,70);
  assert.equal(parsed.records[69].x,context.state.atoms[69].x);
  assert.equal(parsed.records[69].seq,'42A');assert.equal(parsed.records[69].element,'CL');
  assert.equal(parsed.records[0].element,'FE');assert.equal(result.payload.includeDiagnostics,false);
  assert.equal(run("chainsEquivalentForAnalysis('copy_69','A')"),false);
});

test('CIF explicit bonds never join distinct insertion-code residues', () => {
  const {run,context}=parserHarness();
  context.input=serializeStructureRecords([atom({atomName:'X1'}),atom({atomName:'X2',x:1}),atom({seq:'42A',atomName:'X1',x:20}),atom({seq:'42A',atomName:'X2',x:21})])+`loop_\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_auth_asym_id\n_struct_conn.ptnr1_auth_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.pdbx_ptnr1_PDB_ins_code\n_struct_conn.ptnr2_auth_asym_id\n_struct_conn.ptnr2_auth_seq_id\n_struct_conn.ptnr2_label_atom_id\n_struct_conn.pdbx_ptnr2_PDB_ins_code\ncovale A 42 X1 A A 42 X2 A\n#\n`;
  const parsed=run('parseMmcifAtoms(input)');
  assert.deepEqual(snapshot(parsed.explicitBonds),[[2,3]]);
});

test('a stale file read cannot replace a newer structure or clear its status', async () => {
  const fn=html.match(/^      async function loadStructureFromFile[^\n]+\n[\s\S]*?^      \}/m)[0];
  const state={viewer:{},loaded:null};let current=0;const history=[];const noop=()=>{};
  const ctx=vm.createContext({state,cancelActiveAnalysisRun:noop,nextLoadId:()=>++current,captureStructureLoadSnapshot:()=>({}),logDebug:noop,
    setStructureDbInfo:noop,setViewerStructureTitle:noop,setStatus:message=>history.push(['status',message]),appendLoadDiagnostic:noop,
    inferFormat:()=> 'pdb',prepareStructureCoordinateData:async()=>{},clearViewer:async()=>{history.push(['clear']);state.loaded=null;},
    storeStructureAtoms:text=>{state.loaded=text;state.atoms=[];state.structureData={chains:new Map()};history.push(['store',text]);},
    tryLoadCompanionPredictionJsonForStructure:async()=>{},runWithMiniLoadRenderSuppressed:fn=>fn(),tryLoadChapiMeshFromBackend:async(_p,o)=>o.loadId===current,
    isCurrentLoad:id=>id===current,runStartupRevealCameraAnimation:async()=>true,applyViewerStyle:async()=>{},updatePaePanel:noop,updateViewerStructureDbInfoUI:noop,flashStatus:noop,attachLoadDiagnostics:e=>e});
  vm.runInContext(fn,ctx);let finishRead;
  const pending=new Promise(resolve=>{finishRead=resolve;});
  const first=ctx.loadStructureFromFile({name:'slow.pdb',text:()=>pending});
  await ctx.loadStructureFromFile({name:'new.pdb',text:async()=> 'new'});
  const after=history.length;finishRead('old');await first;
  assert.equal(state.loaded,'new');assert.equal(state.structureLabel,'new.pdb');assert.equal(history.length,after);
});

test('superseded parser worker is terminated and cannot overwrite the coordinate cache', async () => {
  const {run,context}=parserHarness();const workers=[];
  class Worker { constructor(){workers.push(this);} postMessage(payload){this.payload=payload;} terminate(){this.terminated=true;} }
  context.Worker=Worker;context.state.loadId=1;
  const first=run("prepareStructureCoordinateData('first','pdb',1)");
  run('nextLoadId()');assert.equal(await first,null);assert.equal(workers[0].terminated,true);
  const second=run("prepareStructureCoordinateData('second','pdb',2)");
  workers[1].onmessage({data:{data:{records:[atom()]}}});await second;
  workers[0].onmessage({data:{data:{records:[]}}});
  assert.equal(run('getStructureCoordinateData.cache.text'),'second');
});
