import assert from 'node:assert/strict';
import test from 'node:test';
import { harness } from './roami-test-harness.mjs';

const site = (chain, seq, atoms, centroid) => ({
  id: `${chain}:${seq}:${atoms.map(atom => atom.atomName).join(',')}`,
  kind: 'group', label: 'nonpolar carbon segment', residue: { chain, seq, resName: 'LEU' }, atoms, centroid,
});
const atom = (id, atomName, coordinates) => ({ id, atomName, coordinates, element: 'C', modelId: '1', altloc: '' });
const fixture = () => {
  const a = site('A', '1', [atom('a1','CB',[0,0,0]),atom('a2','CG',[0,2,0])], [0,1,0]);
  const b = site('B', '2', [atom('b1','CB',[3.5,0,0]),atom('b2','CG',[3.5,2,0])], [3.5,1,0]);
  const raw = (index) => ({
    id: `raw:${index}`, residueA: { ...a.residue, atom: a.atoms[index].atomName },
    residueB: { ...b.residue, atom: b.atoms[index].atomName }, distance: 3.5,
    semantics: { version: 1, identity: `hydrophobic:${index}`, family: 'hydrophobic', directionality: 'symmetric',
      direction: { from: null, to: null, certainty: 'not_applicable' },
      participants: [{ side:'A',role:'contact_atom',site:{...a,kind:'atom',atoms:[a.atoms[index]],centroid:undefined} },
        { side:'B',role:'contact_atom',site:{...b,kind:'atom',atoms:[b.atoms[index]],centroid:undefined} }],
      geometry: { distance:{value:3.5,kind:'atom_pair',unit:'angstrom',source:'input_coordinates'},measurements:[] },
      evidence:{level:'chemically_supported',chemicalCompatibility:'supported',geometrySupport:'supported',chemicalState:'input',ambiguityFlags:[]},
    },
  });
  const supports = [raw(0), raw(1)];
  return { ...supports[0], __supportingContacts: supports, __canonicalContacts: supports,
    __displayGroup: { id:'hydrophobic:patch1',family:'hydrophobic',rule:'connected_nonpolar_feature_pair_patch',
      features:[a,b], counts:{canonicalContacts:2,supportingAtomContacts:2,chemicalUnits:1,displayObjects:1},
      geometry:{closestAtomDistance:3.5,centroidDistance:3.5},ambiguityFlags:[] },
  };
};

test('feature labels separate one natural unit, raw support count, and measured geometry', () => {
  const { context, run } = harness();
  context.contact = fixture();
  assert.equal(run('formatContactNaturalUnitCount([contact])'), '1 nonpolar feature contact');
  assert.match(run('formatContactSemanticPartners(contact)'), /CB, CG.*nonpolar carbon segment.*↔.*CB, CG/);
  assert.equal(run('formatContactGeometryLabel(contact)'), 'Closest supporting atom distance: 3.50 Å');
  assert.match(run('formatContactDisplayGroupDetails(contact)'), /One nonpolar feature contact; 2 supporting atom contacts/);
  assert.match(run('formatContactDisplayGroupDetails(contact)'), /not independent measures of interaction strength/);
  context.contact.__displayGroup.geometry.centroidDistance = 99;
  assert.equal(run('formatContactGeometryLabel(contact)'), 'Closest supporting atom distance: 3.50 Å');
  context.contact.__displayGroup.geometry.closestAtomDistance = null;
  assert.equal(run('formatContactGeometryLabel(contact)'), 'Closest supporting atom distance unavailable');
  context.contact.__displayGroup.counts.supportingAtomContacts = null;
  assert.match(run('formatContactDisplayGroupDetails(contact)'), /total supporting atom-contact count unavailable; 2 exported contacts available for inspection/);
  context.contact.semantics.geometry.distance = {kind:'atom_ring_centroid',value:4.2,unit:'angstrom',source:'arpeggio_atom_plane'};
  assert.equal(run('formatContactGeometryLabel(contact)'), 'Atom–ring centroid distance: 4.20 Å');
  assert.equal(run('resolveContactDistanceValue(contact)'),4.2);
  assert.doesNotMatch(run('formatContactGeometryLabel(contact)'), /Closest.*4.20/);
  context.contact.__supportingContacts = [];
  assert.match(run('formatContactDisplayGroupDetails(contact)'), /supporting atom contacts not exported/);
  context.contact.semantics.geometry.distance = {kind:'atom_pair',value:4.2,unit:'angstrom'};
  assert.equal(run('Number.isNaN(resolveContactDistanceValue(contact))'),true,'a missing group closest distance is not replaced by a representative atom distance');
});

test('each family has a natural unit and raw supporting records never inherit the group', () => {
  const { context, run } = harness();
  context.contact = fixture();
  const expected = {hydrophobic:'nonpolar feature contact',salt_bridge:'charge-group pair',hbond:'donor–acceptor pair',
    halogen_bond:'halogen donor–acceptor pair',metal_coordination:'metal–donor contact',pi_pi:'ring pair',
    pi_cation:'cation–ring pair',water_pi:'water–ring pair',aromatic_packing:'ring pair',clash:'clashing atom pair',
    polar_contact:'polar atom pair',other:'atom pair'};
  for (const [family, label] of Object.entries(expected)) {
    context.contact.__displayGroup.family = family;
    assert.equal(run('contactNaturalUnitLabel(contact)'), label);
  }
  context.contact.__supportingAtomContact = true;
  assert.equal(run('getContactDisplayGroup(contact)'), null);
  assert.equal(run('contactNaturalUnitLabel(contact)'), 'supporting atom contact');
  assert.equal(run('formatContactGeometryLabel(contact)'), 'Atom-pair separation: 3.50 Å');
});

test('one feature visual uses schematic centers and does not scale strength by raw multiplicity', () => {
  const { context, run } = harness();
  context.contact = fixture();
  context.SULFUR_HYDROPHOBIC_TINT = '#ffff00';
  context.TUNING_DEFAULTS = { forcefield:{hydrophobicOpacity:0.4,hydrophobicRadius:1.2} };
  run(`buildForcefieldBridge = (a,b,options) => {
    const object = new THREE.Group();
    object.userData.options = options;
    object.userData.endpoints = [a.toArray(),b.toArray()];
    return object;
  };`);
  const draw = () => run(`(() => {
    const group = new THREE.Group();
    addContactCategoryVisuals(group,contact,null,TUNING_DEFAULTS.forcefield);
    return group;
  })()`);
  const initial = draw();
  assert.equal(initial.children.length, 1);
  assert.equal(initial.children[0].userData.chemicalUnits, 1);
  assert.deepEqual(Array.from(initial.children[0].userData.endpoints, point => Array.from(point)), [[0,1,0],[3.5,1,0]]);
  const firstOptions = initial.children[0].userData.options;
  context.contact.__displayGroup.counts.supportingAtomContacts = 20;
  context.contact.__displayGroup.counts.canonicalContacts = 20;
  const repeated = draw();
  assert.equal(repeated.children.length, 1);
  for (const key of ['opacity','fieldRadius','motionScale','pulseScale']) {
    assert.equal(repeated.children[0].userData.options[key], firstOptions[key]);
  }
  context.contact.__displayGroup.renderGeometry = {
    kind:'hydrophobic_patch',points:[[0,12,0],[3.5,12,0]],schematic:true,source:'supporting_atom_patch_centroids',
  };
  const patch = draw();
  assert.deepEqual(Array.from(patch.children[0].userData.endpoints, point => Array.from(point)), [[0,12,0],[3.5,12,0]]);
  assert.equal(run('formatContactGeometryLabel(contact)'), 'Closest supporting atom distance: 3.50 Å');
});

test('three atom contacts against one Phe ring produce one actual animated feature object', () => {
  const { context, run } = harness();
  const grouped = fixture();
  const names = ['CG','CD1','CE1','CZ','CE2','CD2'];
  const ringAtoms = names.map((name, index) => atom(`phe:${name}`, name,
    [Math.cos(index * Math.PI / 3), Math.sin(index * Math.PI / 3), 0]));
  const ring = { id:'PHE:phenyl',kind:'ring',label:'phenyl ring',residue:{chain:'A',seq:'30',resName:'PHE'},
    atoms:ringAtoms,centroid:[0,0,0],normal:[0,0,1] };
  const supports = ringAtoms.slice(0,3).map((ringAtom, index) => {
    const raw = structuredClone(grouped.__supportingContacts[0]);
    const partner = atom(`tail:${index}`, `C${index+1}`, [ringAtom.coordinates[0],ringAtom.coordinates[1],3.5]);
    raw.id = `phe-support:${index}`;
    raw.residueA = {...ring.residue,atom:ringAtom.atomName};
    raw.residueB = {chain:'B',seq:'8',resName:'LIG',atom:partner.atomName};
    raw.semantics.identity = raw.id;
    raw.semantics.participants[0].site = {...ring,kind:'atom',atoms:[ringAtom],centroid:undefined};
    raw.semantics.participants[1].site = {id:partner.id,kind:'atom',residue:raw.residueB,atoms:[partner]};
    return raw;
  });
  const tail = {id:'LIG:tail',kind:'group',label:'carbon segment',residue:{chain:'B',seq:'8',resName:'LIG'},
    atoms:supports.map(raw => raw.semantics.participants[1].site.atoms[0]),centroid:[0,0,3.5]};
  context.contact = {...supports[0],__supportingContacts:supports,__canonicalContacts:supports,
    __displayGroup:{...grouped.__displayGroup,features:[ring,tail],
      counts:{canonicalContacts:3,supportingAtomContacts:3,chemicalUnits:1,displayObjects:1}}};
  context.SULFUR_HYDROPHOBIC_TINT = '#f0da67';
  context.TUNING_DEFAULTS = {forcefield:{hydrophobicOpacity:0.4,hydrophobicRadius:1.2}};
  const draw = () => run(`(() => {
    const group=new THREE.Group(); addContactCategoryVisuals(group,contact,null,TUNING_DEFAULTS.forcefield); return group;
  })()`);
  const mesh = draw();
  assert.equal(mesh.children.length,1);
  assert.equal(mesh.children[0].userData.chemicalUnits,1);
  assert.equal(typeof mesh.children[0].userData.dynamicEffect.update,'function');
  assert.match(run('formatContactDisplayGroupDetails(contact)'),/3 supporting atom contacts/);
  assert.match(run('formatContactSemanticPartners(contact)'),/A 30 PHE ring \[CG, CD1, CE1, CZ, CE2, CD2\]/);
  const materials = [];
  mesh.traverse(node => { if(node.material) materials.push(node.material); });
  assert.ok(materials.length > 0, 'actual forcefield meshes were built');
  mesh.children[0].userData.dynamicEffect.update(1);
  const opacity = materials.map(material => material.opacity);
  context.contact.__displayGroup.counts.supportingAtomContacts = 30;
  const second = draw();
  second.children[0].userData.dynamicEffect.update(1);
  const repeatedOpacity = [];
  second.traverse(node => {if(node.material) repeatedOpacity.push(node.material.opacity);});
  assert.deepEqual(repeatedOpacity,opacity,'actual animated material opacity ignores raw support multiplicity');
});

test('supporting records use one actual atom-pair mesh without inherited ring or chemical arrows', () => {
  const { context, run } = harness();
  context.contact = {...fixture().__supportingContacts[0],__supportingAtomContact:true};
  for (const family of ['pi_pi','pi_cation','salt_bridge','hbond','halogen_bond','metal_coordination']) {
    context.contact.semantics.family = family;
    context.contact.semantics.directionality = 'directional';
    context.contact.semantics.direction = {from:'A',to:'B',certainty:'certain'};
    context.contact.semantics.participants[0].role = family === 'halogen_bond' ? 'halogen_donor' : 'donor';
    context.contact.semantics.participants[1].role = 'acceptor';
    const group = run(`(() => {
      const group = new THREE.Group(); addContactCategoryVisuals(group,contact,null,{}); return group;
    })()`);
    assert.equal(group.children.length,1,family);
    assert.equal(group.children[0].geometry.type,'CylinderGeometry');
    assert.equal(group.children[0].geometry.parameters.height,3.5);
    assert.equal(group.children[0].userData.supportingAtomContact,true);
    assert.deepEqual(Array.from(group.children[0].userData.atomPairEndpoints,point=>Array.from(point)),[[0,0,0],[3.5,0,0]]);
    assert.equal(run('formatContactTypeLabel(contact)'),'Supporting atom contact');
    assert.equal(run('formatContactGeometryLabel(contact)'),'Atom-pair separation: 3.50 Å');
    assert.doesNotMatch(run('formatContactSemanticPartners(contact)'),/ring|donor|acceptor|→/);
    assert.equal(run('getContactSemanticDirection(contact)'),null);
    assert.equal(run('isHbondLineRenderableContact(contact)'),false);
  }
});

test('group and supporting atom focus target exact local endpoints while retaining the full context data', () => {
  const { context, run } = harness();
  context.contact = fixture();
  run(`state.viewer = {camera:{position:new THREE.Vector3(0,0,100)},controls:{offset:new THREE.Vector3(0,0,40),up:new THREE.Vector3(0,1,0),radius:40}};
    chooseLeastOccludedFocusDirection = (center, direction) => direction;
    applySequenceBiasToViewState = value => value;
    getSequenceCenterBiasPx = () => 0;
    data = {coreA:[new THREE.Vector3(-200,0,0),new THREE.Vector3(200,0,0)],coreB:[new THREE.Vector3(0,0,0)],
      points:[new THREE.Vector3(-200,0,0),new THREE.Vector3(200,0,0)],
      centerA:new THREE.Vector3(-100,0,0),centerB:new THREE.Vector3(100,0,0),
      contactPointA:new THREE.Vector3(0,0,0),contactPointB:new THREE.Vector3(3.5,0,0)};`);
  const initial = run('computeContactViewState(contact,{data,pinned:true})');
  assert.deepEqual(Array.from(initial.target.toArray()), [1.75,0,0]);
  assert.ok(initial.offset.length() < 40, 'remote atoms do not force a whole-ligand zoom');
  context.contact.__supportingAtomContact = true;
  run('data.contactPointA.set(80,2,0); data.contactPointB.set(83.5,2,0);');
  const second = run('computeContactViewState(contact,{data,pinned:true})');
  assert.deepEqual(Array.from(second.target.toArray()), [81.75,2,0]);
  assert.equal(run('data.points.length'),2);
  assert.equal(run('data.coreA[0].x'),-200,'complete structural context is not pruned');
});

test('focus cache preserves ordered endpoints and refreshes revised display geometry', () => {
  const { context, run } = harness();
  context.contact = fixture();
  context.state.atoms = context.contact.__displayGroup.features.flatMap(feature => feature.atoms.map(member => ({
    ...member,...feature.residue,resKey:`${feature.residue.chain}:${feature.residue.seq}`,
    x:member.coordinates[0],y:member.coordinates[1],z:member.coordinates[2],
  })));
  context.state.contactCache = new Map();
  context.FOCUS_NEIGHBOR_SIDECHAIN_RADIUS = 8;
  run(`getResidueAtomIndices = key => key === 'A:1' ? [0,1] : [2,3];
    findNeighborSidechainAtomsForDisplay = () => [];`);
  const first = run('getContactAtomData(contact)');
  context.reversed = structuredClone(context.contact);
  [context.reversed.residueA,context.reversed.residueB] = [context.reversed.residueB,context.reversed.residueA];
  context.reversed.__displayGroup.features.reverse();
  context.reversed.semantics.participants.forEach(participant => {participant.side = participant.side === 'A' ? 'B' : 'A';});
  const reversed = run('getContactAtomData(reversed)');
  assert.notEqual(first,reversed);
  assert.deepEqual(Array.from(first.contactPointA.toArray()),[0,1,0]);
  assert.deepEqual(Array.from(reversed.contactPointA.toArray()),[3.5,1,0]);
  context.contact.__displayGroup.renderGeometry = {schematic:true,points:[[0,5,0],[3.5,5,0]]};
  context.state.reportRevision = 2;
  assert.deepEqual(Array.from(run('getContactAtomData(contact).contactPointA.toArray()')),[0,5,0]);
});

class Element {
  constructor(tag = 'div') { this.tagName = tag; this.children = []; this.style = {}; this.dataset = {}; this.events = {}; this.className = ''; this._text = ''; this._html = ''; }
  appendChild(node) { this.children.push(node); return node; }
  addEventListener(type, callback) { this.events[type] = callback; }
  fire(type, event = {}) { this.events[type]?.(event); }
  set innerHTML(value) { this._html = value; this.children = []; }
  get innerHTML() { return this._html; }
  set textContent(value) { this._text = value; }
  get textContent() { return this._text + this._html.replace(/<[^>]*>/g, '') + this.children.map(child => child.textContent).join(''); }
  get classList() {
    const element = this;
    return { contains: token => element.className.split(' ').includes(token), toggle(token) {
      const values = new Set(element.className.split(' ').filter(Boolean));
      const add = !values.has(token); if (add) values.add(token); else values.delete(token);
      element.className = [...values].join(' '); return add;
    } };
  }
}
const descendants = (element) => element.children.flatMap(child => [child, ...descendants(child)]);

test('contact-rich rows omit zero-contact focus placeholders and label the listed contacts', () => {
  const { context, run } = harness();
  const container = new Element();
  context.document = { createElement: tag => new Element(tag), getElementById: () => container };
  run(`updateActiveListItem = () => {};
    focusResidue = (key, pinned) => { state.testFocusedResidue = key; state.testPinned = pinned; };
    renderHotspots({ perResidue: {
      'Db:459': { chain:'Db', seq:'459', resName:'LYS', total:0 },
      'Db:457': { chain:'Db', seq:'457', resName:'SER', total:2 }
    }});`);
  assert.equal(container.children.length, 1);
  assert.match(container.textContent, /Db 457 SER/);
  assert.match(container.textContent, /2 listed contacts/);
  assert.doesNotMatch(container.textContent, /459|interaction units/);
  container.children[0].fire('click');
  assert.equal(context.state.testFocusedResidue, 'Db:457');
  assert.equal(context.state.testPinned, true);
  run("renderHotspots({perResidue:{'Db:459':{chain:'Db',seq:'459',total:0}}});");
  assert.match(container.textContent, /No contact-rich residues identified/);
});

test('group row lazily expands every supporting atom contact into independently selectable rows', () => {
  const { context, run } = harness();
  const container = new Element();
  context.document = { createElement: tag => new Element(tag), getElementById: () => container };
  context.contact = fixture();
  context.state.contactMap = new Map();
  context.state.collapsedInteractionGroups = new Set();
  run(`contactsForMode = () => [contact];
    classifyContactAnatomy = () => 'sidechain';
    getResidueScopeKeyFromReport = () => '';
    resolvePrimaryResidueKeyForContact = () => 'A:1';
    resolveResidueLabelForContact = () => 'A 1 LEU';
    updateActiveListItem = () => {};
    focusContact = (id, pinned) => { state.selectedTestContact = state.contactMap.get(id); state.testPinned = pinned; };
    renderContactList('hydrophobic', { contacts:{} }, 'testList');`);
  const find = className => descendants(container).filter(node => node.className.split(' ').includes(className));
  assert.equal(find('list-item').length, 0, 'residue bucket starts lazy');
  find('contact-protein-toggle')[0].fire('click');
  assert.equal(find('list-item').length, 1, 'one primary chemical feature');
  const details = find('contact-support-details')[0];
  assert.equal(details.children[0].textContent, '2 supporting atom contacts');
  details.open = true; details.fire('toggle');
  const rows = find('list-item');
  assert.equal(rows.length, 3, 'one feature row plus two raw support rows');
  rows[1].fire('click');
  assert.equal(context.state.selectedTestContact.id, 'raw:0');
  assert.equal(context.state.selectedTestContact.__supportingAtomContact, true);
  assert.equal(context.state.testPinned, true);
  rows[2].fire('click');
  assert.equal(context.state.selectedTestContact.id, 'raw:1');
  assert.equal(context.contact.__supportingContacts[0].__supportingAtomContact, undefined, 'source records stay raw');
  details.open = false; details.fire('toggle'); details.open = true; details.fire('toggle');
  assert.equal(find('list-item').length, 3, 'reopening does not duplicate rows');
  assert.match(container.textContent, /1 interaction · 2 supporting atom contacts/);
});
