import assert from 'node:assert/strict';
import test from 'node:test';
import { readFileSync } from 'node:fs';
import { harness } from './roami-test-harness.mjs';

function fixture() {
  const ra={chain:'A',seq:'123',resName:'LEU'},rb={chain:'B',seq:'456',resName:'PHE'};
  const a={id:'leu:CD1',atomName:'CD1',coordinates:[0,0,0],residue:ra,element:'C'};
  const bs=['CG','CD1','CE1'].map((name,i)=>({id:`phe:${name}`,atomName:name,coordinates:[3.6,i*.3,0],residue:rb,element:'C'}));
  const features={leu:{id:'leu',kind:'group',label:'nonpolar side chain',residue:ra,atoms:[a],centroid:[0,0,0]},phe:{id:'phe',kind:'ring',label:'aromatic ring',residue:rb,atoms:bs,centroid:[3.6,.3,0]}};
  const raw=bs.map((b,i)=>({type:'hydrophobic',source:'test_geometry',atomKeyA:a.id,atomKeyB:b.id,residueA:{...ra,atom:a.atomName},residueB:{...rb,atom:b.atomName},distance:Math.hypot(...b.coordinates),displayGrouping:{groupId:'patch:1',featureIds:['leu','phe'],rule:'nonpolar_feature_patch'},semantics:{version:1,identity:`raw:${i}`,family:'hydrophobic',directionality:'symmetric',participants:[{side:'A',role:'contact_atom',site:{id:a.id,kind:'atom',residue:ra,atoms:[a]}},{side:'B',role:'contact_atom',site:{id:b.id,kind:'atom',residue:rb,atoms:[b]}}],geometry:{distance:{value:Math.hypot(...b.coordinates),kind:'atom_pair',unit:'angstrom'}},evidence:{level:'chemically_supported',supportingRecords:1}}}));
  const group={id:'patch:1',family:'hydrophobic',rule:'nonpolar_feature_patch',featureIds:['leu','phe'],representativeContactId:'raw:0',contactIds:raw.map(c=>c.semantics.identity),ambiguityFlags:[],supportingAtomContacts:bs.map((b,i)=>({id:`support:${i}`,canonicalContactId:`raw:${i}`,atomA:a,atomB:b,distance:raw[i].distance,source:'test_geometry',typing:{A:{atomTypes:['hydrophobe']},B:{atomTypes:['aromatic','hydrophobe']}},provenance:['input_coordinates']})),counts:{canonicalContacts:3,supportingAtomContacts:3,chemicalUnits:1,displayObjects:1},geometry:{closestAtomDistance:3.6,centroidDistance:Math.hypot(3.6,.3)},renderGeometry:{points:[[0,0,0],[3.6,.3,0]],schematic:true}};
  return {contacts:{hydrophobic:raw},displayGroups:{version:1,features,groups:[group]},meta:{interactionSemanticsVersion:1}};
}

test('three raw atom contacts survive JSON while both normalization paths display and count one unit',async()=>{
  const {context,run}=harness();context.report=fixture();const original=JSON.stringify(context.report);
  for(const normalized of [run('normalizeAnalysisReport(report)'),await run('normalizeAnalysisReportAsync(report)')]) {
    context.normalized=normalized;assert.equal(normalized.contacts.hydrophobic.length,3);assert.equal(normalized.perResidue['A:123'].total,1);
    const list=run("contactsForMode(normalized.contacts,'hydrophobic')");assert.equal(list.length,1);context.grouped=list[0];
    assert.deepEqual(JSON.parse(JSON.stringify(run('interactionDisplayCounts([grouped])'))),{chemicalUnits:1,canonicalContacts:3,supportingAtomContacts:3,unknownSupport:0,displayObjects:1});
    const supports=run('getContactSupportingContacts(grouped)');assert.equal(supports.length,3);assert.equal(new Set(supports.map(s=>s.atomKeyB)).size,3);
    for(const child of supports){context.child=child;assert.equal(run('getContactDisplayGroup(child)'),null);assert.match(run('formatContactGeometryLabel(child)'),/Atom-pair separation/);}
    assert.equal(run('getInteractionContactIndex(normalized.contacts).canonical.length'),3);
    assert.equal(run("contactsForMode(normalized.contacts,'hydrophobic')===contactsForMode(normalized.contacts,'hydrophobic')"),true);
  }
  assert.equal(JSON.stringify(context.report),original,'canonical input and metadata not mutated');
});

test('filtering, merged focus reports and proxy remapping retain groups and atom details',async()=>{
  const {context,run}=harness();context.report=fixture();context.normalized=run('normalizeAnalysisReport(report)');
  for(const filtered of [run("filterContactsForResidue(normalized.contacts,'A:123')"),await run("filterContactsForResidueAsync(normalized.contacts,'B:456')")]){
    context.filtered=filtered;assert.equal(run("contactsForMode(filtered,'hydrophobic').length"),1);assert.equal(run("getContactSupportingContacts(contactsForMode(filtered,'hydrophobic')[0]).length"),3);
  }
  context.focused=await run("buildResidueFocusedReport([report,report],'A:123',[['A','B']])");assert.equal(context.focused.contacts.hydrophobic.length,3);assert.equal(run("contactsForMode(focused.contacts,'hydrophobic').length"),1);assert.equal(run("getContactSupportingContacts(contactsForMode(focused.contacts,'hydrophobic')[0]).length"),3);
  context.reverse=new Map([['A','chain-long-A'],['B','chain-long-B']]);
  for(const mapped of [run('remapAnalyzeReportProxyChains(report,reverse)'),await run('remapAnalyzeReportProxyChainsAsync(report,reverse)'),await run('remapAnalyzeReportProxyChainsContactsOnlyAsync(report,reverse)')]){
    assert.equal(mapped.displayGroups.features.leu.residue.chain,'chain-long-A');assert.equal(mapped.displayGroups.groups[0].supportingAtomContacts[0].atomB.residue.chain,'chain-long-B');assert.equal(mapped.displayGroups.groups[0].id,'patch:1');context.mapped=mapped;assert.equal(run("contactsForMode(normalizeAnalysisReport(mapped).contacts,'hydrophobic').length"),1);
  }
});

test('subset projection retains only scoped support and multiplicity never changes confidence',()=>{
  const {context,run}=harness();context.report=fixture();context.normalized=run('normalizeAnalysisReport(report)');context.subset={hydrophobic:[context.normalized.contacts.hydrophobic[1]]};run('inheritInteractionDisplayGroups(normalized.contacts,subset)');context.grouped=run("contactsForMode(subset,'hydrophobic')[0]");
  assert.equal(context.grouped.__displayGroup.counts.supportingAtomContacts,1);assert.equal(run('getContactSupportingContacts(grouped).length'),1);
  const before=run("resolveContactStrengthNormalized(grouped,'hydrophobic')");context.grouped.__displayGroup.counts.supportingAtomContacts=300;assert.equal(run("resolveContactStrengthNormalized(grouped,'hydrophobic')"),before);
  context.raw=context.normalized.contacts.hydrophobic[0];assert.equal(run('getContactDisplayGroup(raw)'),null);assert.equal(run('allInteractionContacts(normalized.contacts).length'),3);
});

test('context attached to a nonrepresentative member remains a view of the same display unit',()=>{
  const {context,run}=harness();context.report=fixture();
  context.report.contacts.hydrophobic[1].basePair={supportCount:1};
  context.report.displayGroups.groups[0].basePairContexts={'raw:1':{supportCount:1}};
  context.normalized=run('normalizeAnalysisReport(report)');
  assert.equal(run("contactsForMode(normalized.contacts,'base_pairing').length"),1);
  assert.equal(run("contactsForMode(normalized.contacts,'hydrophobic').length"),1);
  assert.equal(run('getInteractionContactIndex(normalized.contacts).all.length'),1);
});

test('actual Python producer agrees with frontend grouping, roles, supports and per-family object counts', async () => {
  const {context,run}=harness();
  context.report=JSON.parse(readFileSync(new URL('../roami-tests/display-grouping-2026-09-20/canonical-fixture.json',import.meta.url)));
  const before=new Map(Object.values(context.report.contacts).flat().map(c=>[c.semantics.identity,c]));
  context.normalized=await run('normalizeAnalysisReportAsync(report)');
  const raw=Object.values(context.normalized.contacts).flat();
  assert.equal(raw.length,before.size);
  for(const contact of raw) assert.deepEqual(JSON.parse(JSON.stringify(contact.semantics)),JSON.parse(JSON.stringify(before.get(contact.semantics.identity).semantics)));
  const groups=run('getInteractionContactIndex(normalized.contacts).all');
  assert.equal(groups.length,context.report.displayGroups.groups.filter(g=>!g.debugOnly).length);
  const expected=new Map(context.report.displayGroups.groups.map(g=>[g.id,g]));
  for(const contact of groups){
    context.contact=contact;
    const g=run('getContactDisplayGroup(contact)'), source=expected.get(g.id);
    assert.ok(source);assert.equal(g.counts.chemicalUnits,1);
    assert.equal(run('getContactSupportingContacts(contact).length'),source.supportingAtomContacts.length);
    assert.equal(g.counts.supportingAtomContacts,source.counts.supportingAtomContacts);
    for(const support of run('getContactSupportingContacts(contact)')){
      context.support=support; assert.equal(run('getContactDisplayGroup(support)'),null);
      assert.match(run('formatContactTypeLabel(support)'),/Supporting atom contact/);
    }
  }
  assert.equal(run("contactsForMode(normalized.contacts,'metal').length"),2,'two donors remain two coordination edges');
  assert.ok(groups.some(c=>c.__displayGroup.family==='hydrophobic' && c.__displayGroup.counts.supportingAtomContacts>=3));
  context.window.__PPI_INTERACTION_DEBUG_MODE=true;
  assert.equal(run('getInteractionContactIndex(normalized.contacts).all.length'),context.report.displayGroups.groups.length);
});
