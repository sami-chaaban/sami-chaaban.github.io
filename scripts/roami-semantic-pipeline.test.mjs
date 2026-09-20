import assert from 'node:assert/strict';
import test from 'node:test';
import { readFileSync } from 'node:fs';
import { harness } from './roami-test-harness.mjs';

function contact(family, atomA = 'OG', atomB = 'OE1', options = {}) {
  const residueA = { chain: 'A', seq: '423', resName: 'SER', atom: atomA };
  const residueB = { chain: 'B', seq: '157', resName: 'GLU', atom: atomB };
  const roles = options.roles || (family === 'hbond' ? ['donor', 'acceptor'] : ['contact_atom', 'contact_atom']);
  const participants = [residueA, residueB].map((residue, index) => ({
    side: index === 0 ? 'A' : 'B', role: roles[index], roleProvenance: ['arpeggio_atom_typing'],
    site: { kind: options.kinds?.[index] || 'atom', id: options.sites?.[index] || `${residue.chain}:${residue.seq}:${residue.atom}`,
      residue: { chain: residue.chain, seq: residue.seq, resName: residue.resName },
      atoms: [{ id: `${residue.chain}:${residue.seq}:${residue.atom}`, atomName: residue.atom, element: 'O', coordinates: [index * 2.9, 0, 0] }] },
  }));
  return { type: family, source: 'pdbe-arpeggio', residueA, residueB, distance: 99,
    asserted: { family, confidence: 'medium', evidence: ['typed'] },
    semantics: { version: 1, family,
      directionality: family === 'hbond' ? 'directional' : options.directionality || 'symmetric', participants,
      direction: family === 'hbond' ? { from: 'A', to: 'B', certainty: 'inferred' } : { from: null, to: null, certainty: 'not_applicable' },
      geometry: { distance: { value: 2.9, kind: family === 'hbond' ? 'donor_acceptor' : 'atom_pair', unit: 'angstrom', source: 'model_coordinates' }, measurements: [] },
      evidence: { level: 'candidate', chemicalCompatibility: 'supported', geometrySupport: 'partial', chemicalState: 'assumed', ambiguityFlags: ['hydrogen_missing'], biologicalInterpretation: 'not_evaluated' },
    },
  };
}

function reversed(original) {
  const copy = structuredClone(original);
  [copy.residueA, copy.residueB] = [copy.residueB, copy.residueA];
  copy.semantics.participants.reverse();
  for (const participant of copy.semantics.participants) participant.side = participant.side === 'A' ? 'B' : 'A';
  for (const key of ['from', 'to']) {
    const side = copy.semantics.direction[key];
    if (side) copy.semantics.direction[key] = side === 'A' ? 'B' : 'A';
  }
  return copy;
}

test('canonical normalization removes reversal duplicates and retains atom edges and distinct families', async () => {
  const { run, context } = harness();
  const hbond = contact('hbond');
  const hydrophobic = contact('hydrophobic', 'CB', 'CG');
  context.payload = { hydrogen_bonds: [hbond, reversed(hbond)],
    hydrophobic: [hydrophobic, reversed(hydrophobic), contact('hydrophobic', 'CB', 'CD')],
    salt_bridges: [contact('salt_bridge', 'OG', 'OE1', { roles: ['positive_site', 'negative_site'], directionality: 'role_specific' })],
    other: [contact('clash', 'C', 'O')],
  };
  for (const normalized of [run('normalizeContactsPayload(payload)'), await run('normalizeContactsPayloadAsync(payload)')]) {
    context.normalized = normalized;
    assert.equal(normalized.hydrogen_bonds.length, 1);
    assert.equal(normalized.hydrophobic.length, 2);
    assert.equal(normalized.salt_bridges.length, 1);
    assert.equal(normalized.clash.length, 1);
    assert.deepEqual(Array.from(run("['hbond','hydrophobic','electrostatic','clash'].map(mode=>contactsForMode(normalized,mode).length)")), [1, 2, 1, 1]);
    assert.equal(run('shouldSuppressContactByPrecedence(normalized.hydrophobic[0], normalized)'), false);
  }
  assert.equal(hbond.semantics.participants[0].role, 'donor');
  assert.equal(hbond.__groupSize, undefined, 'normalization must not mutate original API records');
});

test('canonical bidentate coordination and separate fused rings survive normalization and category lists', async () => {
  const { run, context } = harness();
  context.payload = {
    metal_coordination: ['OE1', 'OE2'].map(atom => contact('metal_coordination', 'CA', atom, { roles: ['metal_center', 'coordinating_atom'], kinds: ['metal', 'atom'], directionality: 'role_specific' })),
    pi_cation: ['five', 'six'].map(ring => contact('pi_cation', 'NZ', 'CD2', { roles: ['cation', 'aromatic_ring'], kinds: ['atom', 'ring'], sites: ['LYS:NZ', `TRP:ring:${ring}`], directionality: 'role_specific' })),
    pi_pi: ['five', 'six'].map(ring => contact('pi_pi', 'CD2', 'CG', { roles: ['aromatic_ring', 'aromatic_ring'], kinds: ['ring', 'ring'], sites: [`TRP:ring:${ring}`, 'PHE:ring:six'] })),
  };
  const sync = run('normalizeContactsPayload(payload)');
  const asyncOutput = await run('normalizeContactsPayloadAsync(payload)');
  for (const result of [sync, asyncOutput]) {
    assert.equal(result.metal_coordination.length, 2);
    assert.equal(result.pi_cation.length, 2);
    assert.equal(result.pi_pi.length, 2);
    context.normalized = result;
    assert.equal(run("contactsForMode(normalized, 'aromatic').length"), 4);
    assert.equal(run("contactsForMode(normalized, 'metal').length"), 2);
  }
});

test('duplicate evidence and independent base-pair context are retained through JSON and normalization', () => {
  const { run, context } = harness();
  const first = contact('hbond');
  const other = reversed(first);
  other.basePair = { isCanonicalAtomPattern: false, supportCount: 2 };
  other.semantics.participants.find(p => p.role === 'donor').roleProvenance = ['supplied_chemical_state'];
  other.semantics.geometry.measurements = [{ kind: 'donor_acceptor', value: 2.9, unit: 'angstrom', source: 'coordinate_check' }];
  other.asserted.evidence = ['second_engine_record'];
  context.payload = JSON.parse(JSON.stringify({ hydrogen_bonds: [first, other] }));
  const result = run('normalizeContactsPayload(payload).hydrogen_bonds[0]');
  assert.deepEqual(Array.from(result.semantics.participants.find(p => p.role === 'donor').roleProvenance), ['arpeggio_atom_typing', 'supplied_chemical_state']);
  assert.equal(result.semantics.geometry.measurements[0].source, 'coordinate_check');
  assert.equal(result.basePair.supportCount, 2);
  context.normalized = { hydrogen_bonds: [result] };
  assert.equal(run("contactsForMode(normalized, 'hbond').length"), 1);
  assert.equal(run("contactsForMode(normalized, 'base_pairing').length"), 1);
});

test('geometry values and family come from canonical semantics without inferred replacements', () => {
  const { run, context } = harness();
  context.contact = contact('hbond');
  context.contact.type = 'hydrophobic';
  context.contact.asserted.family = 'other';
  assert.equal(run('getContactCategory(contact)'), 'hbond');
  assert.equal(run('resolveContactDistanceValue(contact)'), 2.9);
  context.contact.semantics.geometry.distance.value = null;
  assert.equal(run('Number.isNaN(resolveContactDistanceValue(contact))'), true);
  assert.equal(run('getContactDistanceValue(contact)'), Infinity);
});

test('chain alias remapping preserves nested role sites and opaque identity in sync and async paths', async () => {
  const { run, context } = harness();
  const original = contact('hbond');
  original.semantics.identity = 'opaque-physical-edge';
  context.payload = { hydrogen_bonds: [original] };
  run("globalThis.reverseMap = new Map([['A','protein-long'],['B','ligand-long']]);");
  const outputs = [run('remapAnalyzeContactsPayload(payload, reverseMap)'), await run('remapAnalyzeContactsPayloadAsync(payload, reverseMap)')];
  for (const output of outputs) {
    const result = output.hydrogen_bonds[0];
    assert.equal(result.residueA.chain, 'protein-long');
    assert.equal(result.semantics.participants[0].site.residue.chain, 'protein-long');
    assert.equal(result.semantics.participants[1].site.residue.chain, 'ligand-long');
    assert.equal(result.semantics.identity, 'opaque-physical-edge');
    assert.deepEqual(JSON.parse(JSON.stringify(result.semantics.direction)), { from: 'A', to: 'B', certainty: 'inferred' });
  }
  assert.equal(original.semantics.participants[0].site.residue.chain, 'A');
});

test('canonical asynchronous normalization honors cancellation rather than delivering stale roles', async () => {
  const { run, context } = harness();
  context.payload = { hydrogen_bonds: [contact('hbond')] };
  await assert.rejects(run('normalizeContactsPayloadAsync(payload, {signal:{aborted:true}})'), { name: 'AbortError' });
});

test('canonical empty reports never acquire frontend-inferred chemistry in focused views', async () => {
  const { run, context } = harness();
  context.report = { contacts: {}, meta: { interactionSemanticsVersion: 1, engine: 'typed-backend' } };
  run(`
    buildIntraResidueMetalCoordinationContacts = () => { throw new Error('inferred metal'); };
    inferHydrophobicContactsForFocusedResidueAsync = () => { throw new Error('inferred hydrophobic'); };
    inferNucleobaseHydrogenBondContactsForAnalyzePairAsync = () => { throw new Error('inferred base pair'); };
    getResidueStubAtom = () => null;
  `);
  const focused = await run("buildResidueFocusedReport([report], 'A:423', [['A','B']])");
  assert.equal(focused.meta.interactionSemanticsVersion, 1);
  assert.equal(Object.values(focused.contacts).flat().length, 0);
  assert.equal(await run("appendNucleobaseHydrogenBondFallbackToAnalyzeReportAsync(report,'A','B')"), context.report);
});

test('base-pair views retain canonical atom edges beyond the legacy residue-pair display budget', () => {
  const { run, context } = harness();
  context.payload = { hydrogen_bonds: Array.from({ length: 9 }, (_, index) => ({
    ...contact('hbond', `N${index}`, `O${index}`), basePair: { isCanonicalAtomPattern: false },
  })) };
  assert.equal(run("contactsForMode(normalizeContactsPayload(payload), 'base_pairing').length"), 9);
});

test('an intra-residue edge contributes once to residue totals in both normalization paths', async () => {
  const { run, context } = harness();
  const intra = contact('metal_coordination', 'ZN', 'N1', { roles: ['metal_center','coordinating_atom'] });
  intra.residueB = { ...intra.residueA, atom: 'N1' };
  context.payload = { metal_coordination: [intra] };
  const sync = run('derivePerResidueFromContacts(payload)');
  const asyncOutput = await run('derivePerResidueFromContactsAsync(payload)');
  assert.equal(sync['A:423'].total, 1);
  assert.equal(asyncOutput['A:423'].total, 1);
});

const producedReport = JSON.parse(readFileSync(new URL('../roami-tests/interaction-semantics-2026-09-20/canonical-fixture.json', import.meta.url), 'utf8'));

test('production Python report retains every canonical identity, role and metric through JavaScript normalization', async () => {
  const { run, context } = harness();
  context.report = structuredClone(producedReport);
  const originals = Object.values(producedReport.contacts).flat();
  assert.equal(originals.length, 16);
  assert.equal(new Set(originals.map(c => c.semantics.identity)).size, originals.length);
  for (const result of [run('normalizeAnalysisReport(report)'), await run('normalizeAnalysisReportAsync(report)')]) {
    const delivered = Object.values(result.contacts).flat();
    assert.equal(delivered.length, originals.length);
    for (const original of originals) {
      const actual = delivered.find(c => c.semantics.identity === original.semantics.identity);
      assert.ok(actual, original.semantics.identity);
      assert.deepEqual(JSON.parse(JSON.stringify(actual.semantics)), JSON.parse(JSON.stringify(original.semantics)));
      context.contact = actual;
      assert.equal(run('getContactCategory(contact)'), original.semantics.family);
      assert.equal(run('resolveContactDistanceValue(contact)'), original.semantics.geometry.distance.value);
    }
  }
});

test('production Python role records agree with frontend arrows, ring features and panel counts', () => {
  const { run, context } = harness();
  context.report = structuredClone(producedReport);
  run('globalThis.normalized = normalizeAnalysisReport(report);');
  assert.deepEqual(Array.from(run("['hbond','electrostatic','halogen','hydrophobic','metal','aromatic','clash'].map(mode => contactsForMode(normalized.contacts,mode).length)")), [5,1,1,2,2,4,1]);
  for (const contact of Object.values(producedReport.contacts).flat()) {
    context.contact = contact;
    const s = contact.semantics;
    const arrow = run('buildContactSemanticArrow(contact)');
    const directional = s.directionality === 'directional' && ['certain','inferred'].includes(s.direction.certainty);
    assert.equal(Boolean(arrow), directional, `${s.family}: ${s.direction.certainty}`);
    const label = run('formatContactSemanticPartners(contact)');
    assert.equal(label.includes('→'), directional);
    assert.doesNotMatch(run('formatContactTypeLabel(contact)'), /confirmed|strong|stabiliz|important|favorable/);
    if (s.family === 'pi_pi') {
      const rings = run('resolvePiPiParticipants(contact)');
      assert.ok(rings);
      assert.equal(rings.ringA.siteId,s.participants.find(p=>p.side==='A').site.id);
      assert.equal(rings.ringB.siteId,s.participants.find(p=>p.side==='B').site.id);
    }
    if (s.family === 'pi_cation') {
      assert.equal(run('resolvePiCationParticipants(contact).ring.siteId'),s.participants.find(p=>p.role==='aromatic_ring').site.id);
    }
  }
});

test('duplicate source states survive merge and side-specific ring metrics keep the correct partner', () => {
  const { run, context } = harness();
  const first = contact('pi_pi', 'CG', 'CZ', {roles:['aromatic_ring','aromatic_ring'], kinds:['ring','ring']});
  first.semantics.geometry.measurements = [{kind:'ring_plane_separation_A', value:2, unit:'angstrom'}];
  first.semantics.participants[0].charge = {value:0,source:'first',inferred:true};
  const second = reversed(first);
  second.semantics.geometry.measurements = [{kind:'ring_plane_separation_A', value:3, unit:'angstrom'}];
  second.semantics.participants[1].charge = {value:1,source:'second',inferred:true};
  context.payload = {pi_pi:[first,second]};
  const result = run('normalizeContactsPayload(payload).pi_pi[0].semantics');
  assert.equal(result.evidence.level,'ambiguous');
  assert.ok(result.evidence.ambiguityFlags.includes('charge_state_observations_differ'));
  assert.equal(result.sourceObservations.length,2);
  assert.deepEqual(JSON.parse(JSON.stringify(result.geometry.measurements)),[
    {kind:'ring_plane_separation_A',value:2,unit:'angstrom'},
    {kind:'ring_plane_separation_B',value:3,unit:'angstrom'},
  ]);
  assert.equal(first.semantics.sourceObservations,undefined);
});

test('base-pair context does not move a chemically unresolved polar edge into the H-bond panel', () => {
  const {run,context} = harness();
  context.payload = {polar_contacts:[{...contact('polar_contact'),basePair:{annotation:'paired_residue_context'}}],hydrogen_bonds:[contact('hbond','N','O')]};
  run('globalThis.normalized = normalizeContactsPayload(payload);');
  assert.equal(run("collectContactsForPanelMode(normalized, 'hbond').length"),1);
  assert.equal(run("contactsForMode(normalized, 'base_pairing')[0].semantics.family"),'polar_contact');
  assert.equal(run("contactsForMode(normalized, 'polar').length"),1);
  run("state.interactionVisibility = {base_pairing:false};");
  assert.equal(run("contactPassesInteractionVisibility(normalized.polar_contacts[0], '', false, {contacts:normalized})"),false);
  assert.equal(run("contactPassesInteractionVisibility(normalized.hydrogen_bonds[0], '', false, {contacts:normalized})"),true);
  const participant=context.payload.polar_contacts[0].semantics.participants[0];
  participant.role='unresolved';
  participant.roleAlternatives=[{role:'donor',source:'arpeggio_openbabel_atom_types'},{role:'acceptor',source:'standard_residue_template'}];
  assert.match(run('formatContactSemanticEvidence(payload.polar_contacts[0])'),/role alternatives: donor \(arpeggio openbabel atom types\) or acceptor \(standard residue template\)/);
  assert.doesNotMatch(run('formatContactSemanticPartners(payload.polar_contacts[0])'),/→/);
});
