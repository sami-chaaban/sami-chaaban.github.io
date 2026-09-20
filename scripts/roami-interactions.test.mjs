import assert from 'node:assert/strict';
import test from 'node:test';
import { harness } from './roami-test-harness.mjs';

test('halogen precedence suppresses only the same atom pair, in either endpoint order', () => {
  const { run } = harness();
  assert.deepEqual(Array.from(run(`(() => {
    const halogen = { type:'halogen_bond', residueA:{chain:'L',seq:'1',resName:'LIG',atom:'CL1'}, residueB:{chain:'A',seq:'2',resName:'MET',atom:'O'} };
    const distinct = { type:'hydrophobic', residueA:{...halogen.residueA,atom:'C9'}, residueB:{...halogen.residueB,atom:'CE'} };
    const same = {...halogen,type:'hydrophobic'};
    const reverse = {...same,residueA:same.residueB,residueB:same.residueA};
    const atomless = {...distinct,residueA:{...distinct.residueA,atom:''}};
    const contacts = {halogen_bonds:[halogen],hydrophobic:[distinct,same,reverse]};
    return [distinct,same,reverse,atomless].map(c=>shouldSuppressHydrophobicOrPackingByHalogenPrecedence(c,contacts));
  })()`)), [false, true, true, true]);
});

test('water, metal and ligand atom names cannot create a polymer backbone endpoint', () => {
  const { run } = harness();
  assert.deepEqual(Array.from(run(`[
    ['SER','OG','HOH','O'], ['ASP','OD1','CA','CA'], ['ASN','ND2','LIG','N'],
    ['SER','N','HOH','O'], ['A',"O2'",'LIG','N']
  ].map(([ra,aa,rb,ab])=>classifyContactAnatomy({type:'hbond',residueA:{chain:'A',seq:'1',resName:ra,atom:aa},residueB:{chain:'B',seq:'2',resName:rb,atom:ab}}))`)),
  ['sidechain', 'sidechain', 'sidechain', 'backbone', 'backbone']);
});

test('normalizers retain every distinct cation-ring pair and direct clashes', async () => {
  const { run } = harness();
  run(`globalThis.payload = {
    pi_cation:[1,2,3].map(n=>({source:'arpeggio',type:'pi_cation',distance:3+n/10,ringPairKey:'ring:'+n,residueA:{chain:'A',seq:'1',resName:'TRP',atom:'CD1'},residueB:{chain:'B',seq:String(n),resName:'LYS',atom:'NZ'}})),
    clash:[{source:'arpeggio',type:'clash',distance:1.8,residueA:{chain:'A',seq:'1',resName:'ALA',atom:'CB'},residueB:{chain:'B',seq:'1',resName:'ALA',atom:'CB'}}]
  };`);
  const sync = run('normalizeContactsPayload(payload)');
  const asyncResult = await run('normalizeContactsPayloadAsync(payload)');
  for (const result of [sync, asyncResult]) {
    assert.equal(result.pi_cation.length, 3);
    assert.equal(result.clash.length, 1);
  }
});

test('category index caches counts and refreshes for debug, structure and report changes', () => {
  const { run } = harness();
  const result = run(`(() => {
    const contact=(n,debugOnly=false)=>({type:'hbond',debugOnly,distance:2.9,residueA:{chain:'A',seq:String(n),resName:'SER',atom:'OG'},residueB:{chain:'B',seq:String(n),resName:'ASP',atom:'OD1'}});
    const contacts={hydrogen_bonds:[contact(1),contact(2,true)]};
    const a=contactsForMode(contacts,'hbond');
    const reused=a===contactsForMode(contacts,'hbond');
    window.__PPI_INTERACTION_DEBUG_MODE=true;
    const debug=contactsForMode(contacts,'hbond').length;
    window.__PPI_INTERACTION_DEBUG_MODE=false;
    contacts.hydrogen_bonds.push(contact(3));
    const appended=contactsForMode(contacts,'hbond').length;
    const before=contactsForMode(contacts,'hbond');
    state.structureRevision++;
    const revised=before!==contactsForMode(contacts,'hbond');
    contacts.hydrogen_bonds=[contact(4)];
    return [a.length,reused,debug,appended,revised,computeInteractionModeCounts({contacts}).hbond];
  })()`);
  assert.deepEqual(Array.from(result), [1, true, 2, 2, true, 1]);
});

test('base-pair metadata is visible without replacing the specific H-bond family', () => {
  const { run } = harness();
  assert.deepEqual(Array.from(run(`(() => {
    const contact={type:'hbond',basePair:{canonical:true},distance:2.9,residueA:{chain:'A',seq:'1',resName:'A',atom:'N6'},residueB:{chain:'B',seq:'1',resName:'U',atom:'O4'}};
    const contacts={hydrogen_bonds:[contact]};
    return [contactsForMode(contacts,'hbond').length,contactsForMode(contacts,'base_pairing').length,getContactCategory(contact)];
  })()`)), [1, 1, 'hbond']);
});

test('visual budgets keep canonical contacts intact and always include the focused cation contact', () => {
  const { run } = harness();
  const result = run(`(() => {
    const contacts=[1,2,3,4].map(n=>({source:'arpeggio',type:'pi_cation',distance:3+n/10,ringKeyA:'A:1:TRP:site-'+(n===4?'two':'one'),ringPairKey:'pair:'+n,residueA:{chain:'A',seq:'1',resName:'TRP',atom:'CD1'},residueB:{chain:'B',seq:String(n),resName:'LYS',atom:'NZ'}}));
    const selected=selectContactsForVisualBudget(contacts,contacts[2]);
    const normalized=normalizeContactsPayload({pi_cation:contacts});
    window.__PPI_INTERACTION_DEBUG_MODE=true;
    return [contacts.length,normalized.pi_cation.length,selected.contacts.length,selected.hiddenCount,selected.contacts.includes(contacts[2]),selected.contacts.includes(contacts[3]),selectContactsForVisualBudget(contacts).contacts.length];
  })()`);
  assert.deepEqual(Array.from(result), [4, 4, 3, 1, true, true, 4]);
});

test('precedence caches refresh after a contact bucket changes', () => {
  const { run } = harness();
  assert.deepEqual(Array.from(run(`(() => {
    const a={chain:'L',seq:'1',resName:'LIG',atom:'CL1'}, b={chain:'A',seq:'2',resName:'MET',atom:'O'};
    const contact={type:'hydrophobic',residueA:a,residueB:b};
    const contacts={hydrophobic:[contact],halogen_bonds:[]};
    const before=contactsForMode(contacts,'hydrophobic').length;
    contacts.halogen_bonds.push({...contact,type:'halogen_bond'});
    return [before,contactsForMode(contacts,'hydrophobic').length];
  })()`)), [1, 0]);
});

test('diagnostic reports use separate cache variants and toggling reloads the current scope', async () => {
  const { run } = harness();
  const compact = run("buildAnalyzePairReportCacheKey('A','B')");
  run('window.__PPI_INTERACTION_DEBUG_MODE = true;');
  assert.notEqual(compact, run("buildAnalyzePairReportCacheKey('A','B')"));
  assert.equal(compact, run("buildAnalyzePairReportCacheKey('A','B',{includeDiagnostics:false})"));
  run(`state.report={chainA:'A',chainB:'B',meta:{scope:'residue',residue:'A:42'}}; handleAnalyze=async (options)=>{globalThis.requestedScope=options;};`);
  await run('setInteractionDebugMode(false)');
  assert.equal(run('requestedScope.target'), 'residue');
  assert.equal(run('requestedScope.residueKey'), 'A:42');
});
