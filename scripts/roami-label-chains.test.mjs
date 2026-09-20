import test from 'node:test';
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { parseStructureRecords, serializeStructureRecords, selectConformers } from '../public/ppi/structure_parser.js';
import { harness } from './roami-test-harness.mjs';

const tags=['group_PDB','id','type_symbol','label_atom_id','label_alt_id','label_comp_id','label_asym_id','label_entity_id','label_seq_id','pdbx_PDB_ins_code','Cartn_x','Cartn_y','Cartn_z','occupancy','B_iso_or_equiv','auth_seq_id','auth_comp_id','auth_asym_id','auth_atom_id','pdbx_PDB_model_num'];
const header='data_collision\nloop_\n'+tags.map(tag=>'_atom_site.'+tag).join('\n')+'\n';
const row=(id,label,auth,seq,name,x,labelName=name)=>`ATOM ${id} C ${labelName} . ALA ${label} 1 1 . ${x} 0 0 1 20 ${seq} ALA ${auth} ${name} 1`;
const fixture=header+[
  row(1,'A','A',42,'X1',0),row(2,'A','A',42,'X2',1),
  row(3,'B','A',42,'X1',20),row(4,'B','A',42,'X2',21),
  row(5,'Z','B',7,'CA',40),
].join('\n')+'\n#\n';
const setup=text=>{const h=harness();Object.assign(h.context,{parseStructureRecords,serializeStructureRecords,input:text});return h;};

test('colliding label chains preserve every atom and use the backend canonical namespace', () => {
  const {run}=setup(fixture);
  const parsed=run('parseMmcifAtoms(input)');
  assert.equal(parsed.atoms.length,5);
  assert.deepEqual([...parsed.residueAtoms.keys()],['A[A]:42','A[B]:42','B:7']);
  assert.equal(parsed.hasLabelChainCollisions,true);
  const backend=spawnSync(process.env.ROAMI_PYTHON || 'python3',['-c',
    'import sys,json;sys.path.insert(0,sys.argv[1]);from api.analysis import parse_mmcif_atoms; atoms,_=parse_mmcif_atoms(sys.stdin.read());print(json.dumps([[a.chain_id,a.res_seq,a.atom_name,a.x] for a in atoms]))',
    fileURLToPath(new URL('../ppi',import.meta.url))],{input:fixture,encoding:'utf8'});
  assert.equal(backend.status,0,backend.stderr);
  const frontend=parsed.atoms.map(atom=>[atom.chain,atom.seq,atom.atomName,atom.x]);
  assert.deepEqual(JSON.parse(JSON.stringify(frontend)),JSON.parse(backend.stdout));
});

test('disjoint author residue identities retain their ordinary author-chain display', () => {
  const source=header+[row(1,'X','A',42,'CA',0),row(2,'Y','A',43,'CA',20)].join('\n')+'\n#\n';
  const parsed=parseStructureRecords(source);
  assert.equal(parsed.hasLabelChainCollisions,false);
  assert.deepEqual(parsed.records.map(atom=>atom.chain),['A','A']);
  const collision=source.replace('20 0 0 1 20 43 ALA A','20 0 0 1 20 42 ALA A');
  assert.deepEqual(parseStructureRecords(collision).records.map(atom=>atom.chain),['X','Y']);
});

test('label-aware conformer selection cannot delete atoms in a different molecule', () => {
  const atoms=['X','Y'].map((labelChain,index)=>({chain:'A',authChain:'A',labelChain,modelToken:'1',seq:'1',atomName:'CA',altLoc:'',occupancy:1,x:index*20}));
  assert.equal(selectConformers(atoms).length,2);
});

test('struct_conn label identifiers disambiguate shared author endpoints without inventing bonds', () => {
  const connection=`loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
covale A 42 B 1 X1 A 42 B 1 X2
#\n`;
  const ligandFixture=fixture.replaceAll('ALA','LIG');
  const {run}=setup(ligandFixture+connection);
  assert.deepEqual(JSON.parse(JSON.stringify(run('parseMmcifAtoms(input).explicitBonds'))),[[2,3]]);
  const ambiguous=setup(ligandFixture+connection.replace('covale A 42 B 1 X1 A 42 B 1 X2','covale A 42 ? ? X1 A 42 ? ? X2'));
  assert.deepEqual(JSON.parse(JSON.stringify(ambiguous.run('parseMmcifAtoms(input).explicitBonds'))),[]);
  const renamed=header+[row(1,'A','A',42,'X1',0,'L1'),row(2,'A','A',42,'X2',1,'L2'),row(3,'B','A',42,'X1',20,'L1'),row(4,'B','A',42,'X2',21,'L2')].join('\n')+'\n#\n';
  const differentAtomNames=setup(renamed.replaceAll('ALA','LIG')+connection.replaceAll('X1','L1').replaceAll('X2','L2'));
  assert.deepEqual(JSON.parse(JSON.stringify(differentAtomNames.run('parseMmcifAtoms(input).explicitBonds'))),[[2,3]]);
});

test('analysis proxy preserves canonical label-chain endpoints and original coordinates', () => {
  const {run,context}=setup(fixture);
  const parsed=run('parseMmcifAtoms(input)');
  Object.assign(context.state,{atoms:parsed.atoms,hasLabelChainCollisions:true,structureText:fixture,structureFormat:'mmcif',
    structureChainAliases:parsed.chainAliases,backendChainAliases:parsed.backendChainAliases});
  const request=run("createAnalyzePayload('A[A]','A[B]')");
  assert.equal(request.payload.chainA,'A[A]');assert.equal(request.payload.chainB,'A[B]');
  assert.deepEqual(parseStructureRecords(request.payload.mmcifText).records.map(a=>[a.chain,a.seq,a.atomName,a.x]),
    JSON.parse(JSON.stringify(parsed.atoms.map(a=>[a.chain,a.seq,a.atomName,a.x]))));
});

test('assembly operations use label identity when a label also names a different author chain', () => {
  const {run,context}=setup(fixture);
  const parsed=run('parseMmcifAtoms(input)');
  context.parsed=parsed;
  context.assembly={chainTransforms:new Map([['B',[{label:'2',isIdentity:false,operations:[]}]]])};
  const normalized=run('normalizeMmcifAssemblyTransformsForParsedStructure(assembly,parsed)');
  assert.deepEqual([...normalized.normalizedAssemblyTransforms.keys()],['A[B]']);
});

test('secondary-structure ranges resolve their label chain without crossing colliding residues', () => {
  const metadata=`_struct_conf.conf_type_id HELX_P
_struct_conf.beg_auth_asym_id A
_struct_conf.end_auth_asym_id A
_struct_conf.beg_label_asym_id B
_struct_conf.end_label_asym_id B
_struct_conf.beg_auth_seq_id 42
_struct_conf.end_auth_seq_id 42
`;
  const {run}=setup(fixture+metadata);
  assert.deepEqual(JSON.parse(JSON.stringify(run('parseMmcifSecondaryRanges(input)'))),[{type:'H',chain:'A[B]',startSeq:'42',endSeq:'42'}]);
});

test('chain descriptions follow label identity even when the displayed name occupies another author namespace', () => {
  const metadata=`loop_
_entity.id
_entity.pdbx_description
1 First 2 Second 3 Third
loop_
_struct_asym.id
_struct_asym.entity_id
A 1 B 2 Z 3
#\n`;
  const {run}=setup(fixture+metadata);
  const descriptions=run('parseMmcif(input).chainDescriptions');
  assert.equal(descriptions.get('A[A]'),'First');
  assert.equal(descriptions.get('A[B]'),'Second');
  assert.equal(descriptions.get('B'),'Third');
});
