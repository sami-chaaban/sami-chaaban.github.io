import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import { parseStructureRecords } from '../public/ppi/structure_parser.js';

const html = readFileSync(new URL('../public/ppi/index.html', import.meta.url), 'utf8');
const functions = [...html.matchAll(/^      (?:async )?function [^\n]+\n[\s\S]*?^      \}/gm)].map(match => match[0]).join('\n');
function parserContext(input) {
  let parseCount = 0;
  const context = vm.createContext({ input, parseStructureRecords: (...args) => { parseCount++; return parseStructureRecords(...args); } });
  vm.runInContext(functions, context);
  return { run: expression => JSON.parse(vm.runInContext(`JSON.stringify(${expression})`, context)), count: () => parseCount };
}

test('biological organism labels exclude recombinant expression hosts', () => {
  const { run } = parserContext({
    rcsb_entity_source_organism: [{ ncbi_scientific_name: 'Homo sapiens' }],
    rcsb_entity_host_organism: [{ ncbi_scientific_name: 'Escherichia coli' }],
    hosts: [{ scientific_name: 'Spodoptera frugiperda' }],
    polymer_entities: [{ entity_src_nat: [{ pdbx_organism_scientific: 'Sus scrofa' }] }],
    entity_src_gen: [{ pdbx_gene_src_scientific_name: 'Homo sapiens', pdbx_host_org_scientific_name: 'Escherichia coli' }],
  });
  assert.deepEqual(run('normalizeProteinOrganismList(input)'), ['Homo sapiens', 'Sus scrofa']);
  assert.deepEqual(run('parseRcsbEntryPayloadForProtein("9GNQ", input).organisms'), ['Homo sapiens', 'Sus scrofa']);
  assert.deepEqual(run('normalizeProteinOrganismList({ host: input.rcsb_entity_host_organism })'), []);
});

const fixture = `data_metadata
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 1 0 0 0 0 1 0 0 0 0 1 0 2 1 0 0
3.8 0 1 0 0 0 0 1 0
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '(1,2)' A 1 1 B
loop_
_entity.id
_entity.pdbx_description
1
;First protein
with multiline description
;
2 'Second protein'
loop_
_struct_asym.id
_struct_asym.entity_id
A 1 B 2
_entity_poly.entity_id 1
_entity_poly.pdbx_strand_id 'X,Y'
loop_
_exptl.method
'X-RAY DIFFRACTION' 'SOLUTION NMR'
_pdbx_helical_symmetry.rotation_per_n_subunits 360
_pdbx_helical_symmetry.rise_per_n_subunits 60
_pdbx_helical_symmetry.number_of_subunits 12
loop_
_struct_conf.conf_type_id
_struct_conf.beg_auth_asym_id
_struct_conf.end_auth_asym_id
_struct_conf.beg_auth_seq_id
_struct_conf.end_auth_seq_id
_struct_conf.pdbx_beg_PDB_ins_code
_struct_conf.pdbx_end_PDB_ins_code
HELX_P A A 42 44 A B HELX_P B B
10 20 ? ?
loop_
_chem_comp.id
_chem_comp.name
LIG
;Long ligand
chemical name
;
FE 'IRON ION'
`;

test('assembly transformations retain packed and wrapped operation rows', () => {
  const { run } = parserContext(fixture);
  const assembly = run(`(() => {const parsed=parseMmcifAssemblyMetadata(input);return {expanded:Boolean(parsed),chains:[...parsed.chainTransforms]};})()`);
  assert.equal(assembly.expanded, true);
  assert.equal(assembly.chains[0][0], 'A');
  assert.equal(assembly.chains[0][1].length, 2);
  assert.equal(assembly.chains[0][1][1].operations[0].v1, 3.8);
  assert.equal(assembly.chains[1][1][0].isIdentity, true);
});

test('all metadata consumers share one token parse and preserve scalar/multiline values', () => {
  const { run, count } = parserContext(fixture);
  assert.deepEqual(run('[...parseMmcifChainDescriptions(input)]'), [
    ['A', 'First protein with multiline description'], ['B', 'Second protein'],
    ['X', 'First protein with multiline description'], ['Y', 'First protein with multiline description'],
  ]);
  assert.deepEqual(run('parseMmcifExperimentMethods(input)'), ['X-RAY DIFFRACTION', 'SOLUTION NMR']);
  assert.deepEqual(run('[...parseMmcifResidueFullNames(input)]'), [['LIG', 'Long ligand chemical name'], ['FE', 'IRON ION']]);
  assert.deepEqual(run('parseMmcifSecondaryRanges(input)'), [
    { type: 'H', chain: 'A', startSeq: '42A', endSeq: '44B' },
    { type: 'H', chain: 'B', startSeq: '10', endSeq: '20' },
  ]);
  const helical = run('parseMmcifHelicalSymmetryMetadata(input)');
  assert.equal(helical.stepAngleDeg, 30);
  assert.equal(helical.stepRise, 5);
  assert.equal(helical.source, 'scalar');
  run('Boolean(parseMmcifAssemblyMetadata(input))');
  assert.equal(count(), 1);
});

test('metadata tag case does not change assembly, chemistry labels or secondary ranges', () => {
  const baseline = parserContext(fixture);
  const upper = parserContext(fixture.replace(/^_\S+|^loop_$/gm, token => token.toUpperCase()));
  for (const expression of ['[...parseMmcifAssemblyMetadata(input).chainTransforms]', 'parseMmcifHelicalSymmetryMetadata(input)', '[...parseMmcifChainDescriptions(input)]', 'parseMmcifExperimentMethods(input)', 'parseMmcifSecondaryRanges(input)', '[...parseMmcifResidueFullNames(input)]']) {
    assert.deepEqual(upper.run(expression), baseline.run(expression));
  }
});

test('single-record categories are accepted without artificial loops', () => {
  const { run } = parserContext(`data_scalar
_entity.id 1
_entity.pdbx_description 'A protein'
_struct_asym.id A
_struct_asym.entity_id 1
_chem_comp.id LIG
_chem_comp.name 'A ligand'
_exptl.method 'ELECTRON MICROSCOPY'
_struct_sheet_range.beg_label_asym_id A
_struct_sheet_range.end_label_asym_id A
_struct_sheet_range.beg_label_seq_id 3
_struct_sheet_range.end_label_seq_id 7
`);
  assert.deepEqual(run('[...parseMmcifChainDescriptions(input)]'), [['A', 'A protein']]);
  assert.deepEqual(run('[...parseMmcifResidueFullNames(input)]'), [['LIG', 'A ligand']]);
  assert.deepEqual(run('parseMmcifExperimentMethods(input)'), ['ELECTRON MICROSCOPY']);
  assert.deepEqual(run('parseMmcifSecondaryRanges(input)'), [{type:'E',chain:'A',startSeq:'3',endSeq:'7'}]);
});
