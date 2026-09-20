import assert from 'node:assert/strict';
import test from 'node:test';
import { harness, html } from './roami-test-harness.mjs';

const atomSite = (id, name, coordinates, resName = 'LIG', chain = 'A', seq = '1', element = 'N') => ({
  kind: 'atom', id, residue: { chain, seq, resName },
  atoms: [{ id, atomName: name, element, coordinates, modelId: 1, altloc: '' }],
});
const participant = (side, role, site) => ({ side, role, site, roleProvenance: ['arpeggio_atom_typing'] });
const contact = (family, a, b, kind = 'atom_pair', direction = null) => ({
  residueA: { ...a.site.residue, atom: a.site.atoms[0]?.atomName },
  residueB: { ...b.site.residue, atom: b.site.atoms[0]?.atomName },
  distance: 999,
  semantics: { version: 1, family, participants: [a, b],
    directionality: direction ? 'directional' : ['salt_bridge','metal_coordination','pi_cation'].includes(family) ? 'role_specific' : 'symmetric',
    direction: direction || { from: null, to: null, certainty: 'not_applicable' },
    geometry: { distance: { value: 3.33, kind, unit: 'angstrom', source: 'input_coordinates' }, measurements: [] },
    evidence: { level: 'candidate', chemicalCompatibility: 'supported', geometrySupport: 'partial', chemicalState: 'assumed', ambiguityFlags: [], biologicalInterpretation: 'not_evaluated' },
  },
});
const ringSite = (id, size, z = 0) => ({
  kind: 'ring', id, residue: { chain: 'R', seq: '5', resName: 'LIG' }, centroid: [0,0,z], normal: [0,0,1],
  atoms: Array.from({ length: size }, (_, i) => ({ id: `${id}:${i}`, atomName: `C${i+1}`, element: 'C', coordinates: [Math.cos(i*2*Math.PI/size),Math.sin(i*2*Math.PI/size),z] })),
});

test('directional labels and arrows follow canonical donor roles even when endpoints are reversed', () => {
  const { run, context } = harness();
  context.contact = contact('hbond',
    participant('A','acceptor',atomSite('glu:OE1','OE1',[3.33,0,0],'GLU','A','157','O')),
    participant('B','donor',atomSite('ser:OG','OG',[0,0,0],'SER','B','423','O')),
    'donor_acceptor',{ from:'B',to:'A',certainty:'inferred' });
  assert.match(run('formatContactSemanticPartners(contact)'), /B 423 SER OG donor → A 157 GLU OE1 acceptor.*direction inferred/);
  assert.equal(run('formatContactDistanceLabel(contact)'), 'Donor–acceptor heavy-atom distance: 3.33 Å');
  assert.equal(run('formatContactTypeLabel(contact)'), 'hydrogen bond (candidate)');
  assert.deepEqual(Array.from(run('buildContactSemanticArrow(contact).userData.semanticStart')), [0,0,0]);
  assert.deepEqual(Array.from(run('buildContactSemanticArrow(contact).userData.semanticEnd')), [3.33,0,0]);
  context.contact.semantics.direction.certainty = 'ambiguous';
  context.contact.semantics.participants.forEach(p => { p.role = 'donor_or_acceptor'; });
  context.contact.donorSide = 'B'; // A legacy hint must not override canonical uncertainty.
  context.contact.asserted = { subtype:'hbond_confirmed' };
  assert.doesNotMatch(run('formatContactSemanticPartners(contact)'), /→/);
  assert.equal(run('buildContactSemanticArrow(contact)'), null);
});

test('backbone and halogen direction comes from metadata, never atom names or a confirmed subtype', () => {
  const { run, context } = harness();
  context.contact = contact('halogen_bond',participant('A','acceptor',atomSite('a','O',[0,0,0])),participant('B','halogen_donor',atomSite('b','CL',[0,3,0],'LIG','B','1','Cl')),'atom_pair',{from:'B',to:'A',certainty:'certain'});
  assert.match(run('formatContactSemanticPartners(contact)'), /CL halogen donor →.*O acceptor/);
  assert.equal(run('buildContactSemanticArrow(contact).userData.semanticDirection.from'), 'B');
  context.contact = contact('hbond',participant('A','donor',atomSite('a','N',[0,0,0],'ALA')),participant('B','acceptor',atomSite('b','O',[3,0,0],'ALA','B','2','O')),'donor_acceptor',{from:'A',to:'B',certainty:'certain'});
  context.contact.semantics.geometry.hydrogen = { id:'h',atomName:'H',element:'H',coordinates:[1,0,0],source:'input_hydrogens',inferred:false };
  assert.match(run('formatContactSemanticEvidence(contact)'), /H: input hydrogens/);
  const legacy = structuredClone(context.contact); delete legacy.semantics;
  legacy.donorSide = 'A'; legacy.asserted = { family:'hbond',subtype:'hbond_confirmed' };
  context.contact = legacy;
  assert.equal(run('getContactSemanticDirection(contact)'), null);
  assert.doesNotMatch(run('formatContactSemanticPartners(contact)'), /→/);
  assert.doesNotMatch(run('formatContactTypeLabel(contact)'), /confirmed|strong/);
});

test('metal display preserves distinct donor atoms and does not invent charge or H-bond arrows', () => {
  const { run, context } = harness();
  const metal = participant('B','metal_center',{...atomSite('zn','ZN',[0,0,0],'ZN','M','1','Zn'),kind:'metal'});
  context.first = contact('metal_coordination',participant('A','coordinating_atom',atomSite('oe1','OE1',[2,0,0],'GLU')),metal,'metal_donor');
  context.second = contact('metal_coordination',participant('A','coordinating_atom',atomSite('oe2','OE2',[0,2,0],'GLU')),metal,'metal_donor');
  assert.match(run('formatContactSemanticPartners(first)'), /ZN metal center —.*OE1 coordinating atom/);
  assert.doesNotMatch(run('formatContactSemanticPartners(first)'), /\+2|→/);
  assert.deepEqual(Array.from(run('resolveMetalCoordinationParticipants(first).ligandPoint.toArray()')), [2,0,0]);
  assert.deepEqual(Array.from(run('resolveMetalCoordinationParticipants(second).ligandPoint.toArray()')), [0,2,0]);
  assert.equal(run('buildContactSemanticArrow(first)'), null);
});

test('salt-bridge colors follow canonical endpoint polarity and never infer the opposite sign', () => {
  const { run, context } = harness();
  context.contact = contact('salt_bridge',participant('A','positive_site',atomSite('p','NZ',[0,0,0],'LYS')),participant('B','negative_site',atomSite('n','OE1',[0,3,0],'GLU','B','2')),'charge_site');
  const colors = () => Array.from(run(`(() => {
    const group=buildElectrostaticBeadBridge(new THREE.Vector3(0,0,0),new THREE.Vector3(0,3,0),{contact,colorPlus:'#ff0000',colorMinus:'#0000ff'});
    return group.children[0].children.filter(m=>m.geometry?.type==='SphereGeometry').slice(0,2).map(m=>m.material.color.getHexString());
  })()`));
  assert.deepEqual(colors(), ['ff0000','0000ff']);
  assert.match(run('formatContactSemanticPartners(contact)'), /NZ \(\+\) ↔.*OE1 \(−\)/);
  context.contact.semantics.participants[1].role = 'unresolved';
  assert.notEqual(colors()[1], '0000ff');
  assert.equal(run('buildContactSemanticArrow(contact)'), null);
});

test('cation–π and π–π render the supplied physical sites including distinct fused-ring membership', () => {
  const { run, context } = harness();
  context.contact = contact('pi_cation',participant('A','aromatic_ring',ringSite('fused-five',5)),participant('B','cation',atomSite('cation','N9',[0,0,4],'LIG','C','1')),'cation_centroid');
  assert.match(run('formatContactSemanticPartners(contact)'), /N9 cation ↔.*ring \[C1, C2, C3, C4, C5\]/);
  assert.equal(run('resolvePiCationParticipants(contact).ring.siteId'), 'fused-five');
  assert.equal(run('resolvePiCationParticipants(contact).ring.atoms.length'),5);
  assert.equal(run('buildContactSemanticArrow(contact)'),null);
  context.contact = contact('pi_pi',participant('A','aromatic_ring',ringSite('fused-five',5)),participant('B','aromatic_ring',ringSite('fused-six',6,3.5)),'ring_centroid');
  const result = run('resolvePiPiParticipants(contact)');
  assert.equal(result.ringA.atoms.length,5);
  assert.equal(result.ringB.atoms.length,6);
  assert.equal(result.distance,3.5);
  assert.equal(run('formatContactDistanceLabel(contact)'), 'Ring centroid distance: 3.33 Å'); // Reported metric, not schematic recomputation.
  assert.equal(run('buildContactSemanticArrow(contact)'),null);
});

test('symmetric contacts stay symmetric and incomplete semantics cannot be treated as canonical evidence', () => {
  const { run, context } = harness();
  for (const family of ['hydrophobic','clash','polar_contact']) {
    context.contact = contact(family,participant('A','contact_atom',atomSite('a','C1',[0,0,0])),participant('B','contact_atom',atomSite('b','C2',[2,0,0],'LIG','B','1')));
    assert.match(run('formatContactSemanticPartners(contact)'), /↔/);
    assert.equal(run('buildContactSemanticArrow(contact)'),null);
    assert.doesNotMatch(run('formatContactTypeLabel(contact)'), /strong|confirmed|stabiliz|favorable/);
  }
  for (const participants of [[],[null,null],[{},{}],[{side:'A',site:{}},{side:'A',site:{}}]]) {
    context.contact.semantics.participants = participants;
    assert.equal(run('getContactSemantics(contact)'),null);
  }
});

test('all canonical distance kinds preserve measured values, units and missingness without schematic substitutions', () => {
  const { run, context } = harness();
  context.contact = contact('hbond',participant('A','donor',atomSite('a','N',[0,0,0])),participant('B','acceptor',atomSite('b','O',[99,0,0],'LIG','B','2')));
  const labels = {
    donor_acceptor:'Donor–acceptor heavy-atom distance',metal_donor:'Metal–donor distance',
    ring_centroid:'Ring centroid distance',closest_atom:'Closest atom distance',charge_site:'Charged-site atom separation',
    cation_centroid:'Cation-atom–ring centroid distance',atom_pair:'Atom-pair separation',
  };
  for (const [kind,label] of Object.entries(labels)) {
    context.contact.semantics.geometry.distance.kind = kind;
    assert.equal(run('formatContactDistanceLabel(contact)'),`${label}: 3.33 Å`);
  }
  context.contact.semantics.geometry.distance.value = null;
  assert.equal(run('formatContactDistanceLabel(contact)'),'Distance unavailable');
  context.contact.semantics.geometry.measurements = [
    {kind:'donor_hydrogen_acceptor_angle',value:152.4,unit:'degree',source:'input_hydrogens'},
    {kind:'vdw_overlap',value:0,unit:'angstrom',source:'vdw_radii'},
    {kind:'missing_angle',value:null,unit:'degree'},
  ];
  assert.match(run('formatContactSemanticMeasurements(contact)'),/152.40 °.*0.00 Å/);
  assert.doesNotMatch(run('formatContactSemanticMeasurements(contact)'),/missing angle/);
  context.contact.energy = -90;
  const before = run('resolveContactStrengthNormalized(contact)');
  context.contact.energy = 100;
  assert.equal(run('resolveContactStrengthNormalized(contact)'),before);
  context.contact.basePair = {annotation:'watson_crick_candidate'};
  assert.match(run('formatContactSemanticEvidence(contact)'),/Base-pair context: watson crick candidate/);
  assert.match(run('formatContactTypeLabel(contact)'),/^hydrogen bond/);
});

test('charged groups connect measured endpoint atoms and retain the full chemical group identity', () => {
  const { run, context } = harness();
  const group = { ...atomSite('group','NZ',[9,9,9],'LYS'),kind:'group',centroid:[8,8,8],
    contactAtom:{id:'contact',atomName:'NZ',element:'N',coordinates:[1,2,3]} };
  context.site = group;
  assert.deepEqual(Array.from(run('resolveSemanticSitePoint(site).toArray()')),[1,2,3]);
  context.contact = contact('salt_bridge',participant('A','positive_site',group),participant('B','negative_site',atomSite('n','OE1',[0,0,0])),'charge_site');
  assert.match(run('formatContactSemanticPartners(contact)'),/group \[NZ\] \(contact atom NZ\)/);
  delete context.contact.semantics;
  assert.equal(run('buildMetalCoordinationBridge(contact)'),null);
  assert.equal(run('buildPiCationBridge(contact)'),null);
});

test('canonical list rows skip source-strength heuristics and keep clashes out of the debug tier', () => {
  const { run, context } = harness();
  for (const name of ['formatContactWeaknessBadge','classifyNeighborInteractionTier']) {
    const definition = html.match(new RegExp(`^        const ${name} = [\\s\\S]*?^        };`, 'm'));
    assert.ok(definition,`missing actual list helper ${name}`);
    run(definition[0]);
  }
  run(`resolveBackendWeaknessLabel=()=>{throw new Error('legacy strength heuristic reached');};
    isContactWeakOutsideTypicalRange=()=>{throw new Error('legacy distance heuristic reached');};`);
  context.contact = contact('clash',participant('A','contact_atom',atomSite('a','N',[0,0,0])),participant('B','contact_atom',atomSite('b','O',[1,0,0])));
  context.contact.strength = 'weak';
  assert.equal(run('formatContactWeaknessBadge(contact)'),'');
  assert.equal(run('classifyNeighborInteractionTier(contact)'),'secondary');
  context.contact.semantics.evidence.level = 'geometrically_supported';
  assert.equal(run('classifyNeighborInteractionTier(contact)'),'primary');
});

test('actual H-bond dash rendering follows donor order and ambiguous roles produce no arrow', () => {
  const { run, context } = harness();
  context.contact = contact('hbond',participant('A','acceptor',atomSite('o','O',[0,0,0],'GLU','A','1','O')),participant('B','donor',atomSite('n','N',[3,0,0],'ALA','B','2','N')),'donor_acceptor',{from:'B',to:'A',certainty:'inferred'});
  run(`state.atoms=[{chain:'A',seq:'1',resName:'GLU',atomName:'O',element:'O',x:0,y:0,z:0},{chain:'B',seq:'2',resName:'ALA',atomName:'N',element:'N',x:3,y:0,z:0}];
    state.visualTuning={hbond:{dash:0.2,gap:0.1,radius:0.05,opacity:0.8}};
    isHbondLineRenderableContact=()=>true;isAtomRenderMode=()=>true;isResidueChemBackboneAtomPreviewActive=()=>false;
    findAtomIndex=(chain,seq,name)=>state.atoms.findIndex(a=>a.chain===chain&&a.seq===seq&&a.atomName===name);`);
  const first = run(`(() => {
    const group=buildHbondLineGroup(new Set(['A:1:O','B:2:N']),null,[contact]);
    const arrow=group.children.find(mesh=>mesh.userData.semanticDirection);
    const dash=group.children.find(mesh=>mesh.userData.isHbondDash);
    return [arrow.userData.semanticStart,arrow.userData.semanticEnd,new THREE.Vector3(0,1,0).applyQuaternion(dash.quaternion).x];
  })()`);
  assert.deepEqual(Array.from(first[0]),[3,0,0]);
  assert.deepEqual(Array.from(first[1]),[0,0,0]);
  assert.ok(first[2] < -0.99);
  context.contact.semantics.direction.certainty = 'ambiguous';
  assert.equal(run(`buildHbondLineGroup(new Set(['A:1:O','B:2:N']),null,[contact]).children.some(mesh=>mesh.userData.semanticDirection)`),false);
});
