"""Chemical role/geometry identity acceptance tests through the real report builder.

Engine observations are controlled here; separate chemistry tests validate upstream
Open Babel typing and deposited 1CLL/1A6M/1BNA structures. No UI inference is used.
"""
import copy
import json
import math
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from api import analysis as a
from test_chemistry_accuracy import atom, raw_pair, classify


def reverse(raw):
    raw = copy.deepcopy(raw)
    raw['bgn'], raw['end'] = raw['end'], raw['bgn']
    if raw.get('hbond_donor_side') in {'A','B'}:
        raw['hbond_donor_side'] = 'B' if raw['hbond_donor_side']=='A' else 'A'
    return raw


def report(atoms, rows):
    aliases = a._identity_chain_aliases({at.chain_id for at in atoms})
    with patch.object(a, '_parse_structure_cached', return_value=(atoms, aliases)), patch.object(a, '_run_arpeggio_contacts', return_value=rows):
        return a.analyze_interface('controlled semantic fixture', 'A', 'B')


def records(result):
    return [r for bucket in result['contacts'].values() for r in bucket]


def roles(sem):
    return {p['side']:p['role'] for p in sem['participants']}


def ring(chain, seq, residue='PHE', center=(0,0,0)):
    names = ['CG','CD1','CE1','CZ','CE2','CD2']
    return [atom(chain,seq,residue,name,'C',(center[0]+1.4*math.cos(i*math.pi/3),center[1]+1.4*math.sin(i*math.pi/3),center[2])) for i,name in enumerate(names)]


def runtime_fixture():
    atoms, rows = [], []
    def add(at, raw):
        atoms.extend(at); rows.append(raw)
    ser=atom('A',1,'SER','OG','O',(0,0,0)); glu=atom('B',1,'GLU','OE1','O',(2.8,0,0)); h=atom('A',1,'SER','HG','H',(0.96,0,0))
    supplied=raw_pair(ser,glu,['hbond','polar']); add([ser,glu,h],supplied); rows.append(reverse(supplied))
    ser2=atom('A',2,'SER','OG','O',(0,10,0)); glu2=atom('B',2,'GLU','OE1','O',(2.8,10,0)); add([ser2,glu2],raw_pair(ser2,glu2,['polar']))
    n=atom('A',3,'ALA','N','N',(0,20,0)); o=atom('B',3,'GLY','O','O',(2.9,20,0)); add([n,o],raw_pair(n,o,['polar']))
    zn=atom('A',4,'ZN','ZN','ZN',(0,30,0)); oe1=atom('B',4,'GLU','OE1','O',(2.1,30,0)); oe2=atom('B',4,'GLU','OE2','O',(0,32.2,0))
    add([zn,oe1,oe2],raw_pair(zn,oe1,['metal_complex'])); rows.append(raw_pair(zn,oe2,['metal_complex']))
    lys=atom('A',5,'LYS','NZ','N',(0,40,0)); asp=atom('B',5,'ASP','OD1','O',(3.0,40,0)); ionic=raw_pair(lys,asp,['ionic']); ionic['bgn'].update(formal_charge=1,formal_charge_source='input_formal_charge'); ionic['end'].update(formal_charge=-1,formal_charge_source='input_formal_charge'); add([lys,asp],ionic)
    br=atom('A',6,'LIG','BR1','BR',(0,50,0)); anchor=atom('A',6,'LIG','C1','C',(-1.8,50,0)); acceptor=atom('B',6,'ASN','OD1','O',(3.3,50,0)); xb=raw_pair(br,acceptor,['xbond']); xb['bgn']['bonded_atoms']=[dict(auth_atom_id='C1',type_symbol='C',coordinates=[-1.8,50,0])]; xb['end']['atom_types']=['hbond acceptor']; add([br,anchor,acceptor],xb)
    ra,rb=ring('A',7,center=(0,60,0)),ring('B',7,center=(0,60,3.5)); add(ra+rb,raw_pair(ra[0],rb[0],['aromatic']))
    rp=ring('B',8,center=(0,70,0)); cat=atom('A',8,'LYS','NZ','N',(0,70,4)); cp=raw_pair(cat,rp[0],['cationpi']); cp['bgn']['formal_charge']=1; add([cat]+rp,cp)
    c1=atom('A',9,'LEU','CD1','C',(0,80,0)); c2=atom('A',9,'LEU','CD2','C',(0,80,1)); d1=atom('B',9,'VAL','CG1','C',(3.7,80,0)); add([c1,c2,d1],raw_pair(c1,d1,['hydrophobic'])); rows.append(raw_pair(c2,d1,['hydrophobic']))
    cc=atom('A',10,'ALA','CB','C',(0,90,0)); dd=atom('B',10,'ALA','CB','C',(2,90,0)); add([cc,dd],raw_pair(cc,dd,['vdw_clash']))
    oh=atom('A',11,'A',"O2'",'O',(0,100,0)); po=atom('B',11,'U','OP1','O',(2.8,100,0)); add([oh,po],raw_pair(oh,po,['polar']))
    his=atom('A',12,'HIS','ND1','N',(0,110,0)); other=atom('B',12,'SER','OG','O',(2.8,110,0)); add([his,other],raw_pair(his,other,['hbond'],hbond_angle=170))
    # Two physical indole rings, from planar fused pentagon/hexagon vertices.
    positions={'CD2':(0,0),'CE2':(0,1.4),'CG':(-1.33,-0.43),'CD1':(-2.15,0.7),'NE1':(-1.33,1.83),'CZ2':(1.21,2.1),'CH2':(2.42,1.4),'CZ3':(2.42,0),'CE3':(1.21,-0.7)}
    trp=[atom('A',13,'TRP',name,'N' if name=='NE1' else 'C',(xy[0],xy[1]+120,0)) for name,xy in positions.items()]
    counterpart=ring('B',13,center=(0,120,3.5)); atoms.extend(trp+counterpart)
    rows += [raw_pair(next(at for at in trp if at.atom_name==name),counterpart[0],['aromatic']) for name in ('CD1','CZ3')]
    return report(atoms,rows)


class InteractionSemanticTests(unittest.TestCase):
    def test_ser_glu_and_backbone_roles_are_chemical_not_row_order(self):
        for res,name,element,accres,accname in [('SER','OG','O','GLU','OE1'),('ALA','N','N','GLY','O')]:
            d=atom('A',1,res,name,element,(0,0,0)); acc=atom('B',2,accres,accname,'O',(2.8,0,0)); raw=raw_pair(d,acc,['polar'])
            sem=classify(raw,[d,acc])['semantics']; flipped=classify(reverse(raw),[d,acc])['semantics']
            self.assertEqual(roles(sem),{'A':'donor','B':'acceptor'}); self.assertEqual(roles(flipped),{'A':'acceptor','B':'donor'})
            self.assertEqual(sem['identity'],flipped['identity']); self.assertEqual(sem['direction']['certainty'],'inferred')
            self.assertEqual(sem['geometry']['distance']['kind'],'donor_acceptor')

    def test_deposited_hydrogen_geometry_and_provenance_survive_report(self):
        d=atom('A',1,'SER','OG','O',(0,0,0)); h=atom('A',1,'SER','HG','H',(.96,0,0)); acc=atom('B',2,'GLU','OE1','O',(2.8,0,0)); sem=records(report([d,h,acc],[raw_pair(d,acc,['hbond'])]))[0]['semantics']
        self.assertEqual(sem['geometry']['hydrogen']['atomName'],'HG'); self.assertFalse(sem['geometry']['hydrogen']['inferred'])
        self.assertEqual(sem['evidence']['level'],'geometrically_supported'); self.assertEqual(sem['direction']['certainty'],'certain')
        self.assertEqual(sem['evidence']['biologicalInterpretation'],'not_evaluated')

    def test_angle_without_actual_hydrogen_cannot_confirm(self):
        d=atom('A',1,'SER','OG','O',(0,0,0)); acc=atom('B',2,'GLU','OE1','O',(2.8,0,0))
        result=classify(raw_pair(d,acc,['hbond'],hbond_angle=180),[d,acc]); sem=result['semantics']
        self.assertEqual(sem['evidence']['level'],'candidate'); self.assertNotIn('hydrogen',sem['geometry']); self.assertEqual(result['subtype'],'hbond_candidate')

    def test_modelled_hydrogen_is_never_deposited_evidence(self):
        d=atom('A',1,'SER','OG','O',(0,0,0)); acc=atom('B',2,'GLU','OE1','O',(2.8,0,0)); raw=raw_pair(d,acc,['hbond'],hbond_donor_side='A',hbond_angle=180,hbond_geometry_source='arpeggio_generated_hydrogens',hbond_hydrogen={'coordinates':[1,0,0]})
        sem=classify(raw,[d,acc])['semantics']; self.assertTrue(sem['geometry']['hydrogen']['inferred']); self.assertEqual(sem['evidence']['level'],'candidate')

    def test_bad_actual_hydrogen_geometry_is_retained_without_support(self):
        d=atom('A',1,'SER','OG','O',(0,0,0)); h=atom('A',1,'SER','HG','H',(-.96,0,0)); acc=atom('B',2,'GLU','OE1','O',(2.8,0,0)); sem=classify(raw_pair(d,acc,['hbond']),[d,h,acc])['semantics']
        self.assertEqual(sem['evidence']['level'],'candidate'); self.assertIn('hydrogen_angle_not_supportive',sem['evidence']['ambiguityFlags'])

    def test_dual_role_histidine_is_explicitly_ambiguous(self):
        d=atom('A',1,'HIS','ND1','N',(0,0,0)); acc=atom('B',2,'SER','OG','O',(2.8,0,0)); sem=classify(raw_pair(d,acc,['hbond'],hbond_angle=170),[d,acc])['semantics']
        self.assertEqual(sem['direction']['certainty'],'ambiguous'); self.assertIsNone(sem['direction']['from']); self.assertEqual(len(sem['direction']['alternatives']),2)
        self.assertIn('histidine_protonation_or_tautomer_uncertain',sem['evidence']['ambiguityFlags'])

    def test_explicit_histidine_tautomers_have_distinct_roles(self):
        self.assertEqual(a._partner_hbond_roles({},'HID','ND1','N'),(True,False,True))
        self.assertEqual(a._partner_hbond_roles({},'HID','NE2','N'),(False,True,True))
        self.assertEqual(a._partner_hbond_roles({},'HIP','NE2','N'),(True,False,True))

    def test_halogen_graph_anchor_and_angle_are_preserved(self):
        report_=runtime_fixture(); sem=report_['contacts']['halogen_bonds'][0]['semantics']
        self.assertEqual(roles(sem),{'A':'halogen_donor','B':'acceptor'}); self.assertEqual(sem['geometry']['donorAnchor']['atomName'],'C1'); self.assertEqual(sem['geometry']['donorAnchor']['source'],'arpeggio_bond_graph')
        self.assertEqual(next(m['value'] for m in sem['geometry']['measurements'] if m['kind']=='donor_anchor_halogen_acceptor_angle'),180)

    def test_unknown_metal_oxidation_and_bidentate_atoms_remain_distinct(self):
        rows=runtime_fixture()['contacts']['metal_coordination']; self.assertEqual(len(rows),2)
        self.assertEqual({r['semantics']['participants'][1]['site']['atoms'][0]['atomName'] for r in rows},{'OE1','OE2'})
        for r in rows:
            self.assertIsNone(r['semantics']['participants'][0]['charge']['value']); self.assertEqual(r['asserted']['denticity'],2)

    def test_explicit_neutral_lysine_cannot_make_salt_bridge(self):
        lys=atom('A',1,'LYS','NZ','N',(0,0,0)); glu=atom('B',2,'GLU','OE1','O',(3,0,0)); raw=raw_pair(lys,glu,['ionic']);raw['bgn']['formal_charge']=0;raw['end']['formal_charge']=-1
        self.assertNotEqual(classify(raw,[lys,glu])['family'],'salt_bridge'); self.assertEqual(a._semantic_charge(raw['bgn'],{'resName':'LYS','atom':'NZ','element':'N'})['value'],0)

    def test_ligand_formal_charges_support_ionic_roles(self):
        n=atom('A',1,'QAM','N1','N',(0,0,0)); o=atom('B',2,'ACT','O1','O',(3,0,0)); raw=raw_pair(n,o,['ionic']);raw['bgn']['formal_charge']=1;raw['end']['formal_charge']=-1
        result=classify(raw,[n,o]);self.assertEqual(result['family'],'salt_bridge'); self.assertEqual(roles(result['semantics']),{'A':'positive_site','B':'negative_site'})

    def test_salt_far_apart_and_same_charge_are_not_asserted(self):
        for separation,right_charge in [(8,-1),(3,1)]:
            n=atom('A',1,'LYS','NZ','N',(0,0,0)); o=atom('B',2,'GLU','OE1','O',(separation,0,0)); raw=raw_pair(n,o,['ionic']);raw['bgn']['formal_charge']=1;raw['end']['formal_charge']=right_charge
            self.assertNotEqual(classify(raw,[n,o])['family'],'salt_bridge')

    def test_cation_pi_and_pi_pi_distances_have_physical_sites(self):
        result=runtime_fixture(); cp=result['contacts']['pi_cation'][0]['semantics']; self.assertEqual(cp['geometry']['distance']['kind'],'cation_centroid');self.assertEqual(cp['geometry']['distance']['value'],4)
        self.assertEqual(roles(cp),{'A':'cation','B':'aromatic_ring'}); self.assertEqual(len(cp['participants'][1]['site']['atoms']),6)
        pp=next(r['semantics'] for r in result['contacts']['pi_pi'] if r['residueA']['seq']=='7');self.assertEqual(pp['geometry']['distance']['kind'],'ring_centroid');self.assertEqual(pp['geometry']['distance']['value'],3.5)

    def test_fused_ring_identities_and_normals_are_not_collapsed(self):
        rows=[r for r in records(runtime_fixture()) if r['residueA']['seq']=='13']; sites={p['site']['id']:p['site'] for r in rows for p in r['semantics']['participants'] if p['side']=='A'}
        self.assertEqual(sorted(len(site['atoms']) for site in sites.values()),[5,6]); self.assertTrue(all(len(site['normal'])==3 for site in sites.values()))

    def test_symmetric_contacts_and_real_atom_pairs_survive(self):
        result=runtime_fixture(); hyd=result['contacts']['hydrophobic'];self.assertEqual(len(hyd),2)
        for record in hyd+[r for r in result['contacts']['other'] if r['type']=='clash']:
            self.assertEqual(record['semantics']['directionality'],'symmetric');self.assertEqual(record['semantics']['direction']['certainty'],'not_applicable')
        clash=[r for r in result['contacts']['other'] if r['type']=='clash'][0]['semantics'];self.assertGreater(next(m['value'] for m in clash['geometry']['measurements'] if m['kind']=='vdw_overlap'),0)

    def test_hbond_distinct_acceptor_oxygens_and_families_survive(self):
        d=atom('A',1,'SER','OG','O',(0,0,0)); x=atom('B',2,'GLU','OE1','O',(2.8,0,0)); y=atom('B',2,'GLU','OE2','O',(0,2.8,0)); rows=records(report([d,x,y],[raw_pair(d,x,['hbond']),raw_pair(d,y,['hbond'])]));self.assertEqual(len(rows),2)
        one=rows[0]['semantics'];other=dict(one,family='salt_bridge');self.assertNotEqual(a._semantic_identity(one),a._semantic_identity(other))

    def test_rna_hydroxyl_phosphate_roles_and_base_context_are_separate(self):
        row=next(r for r in records(runtime_fixture()) if r['residueA']['seq']=='11');self.assertEqual(roles(row['semantics']),{'A':'donor','B':'acceptor'});self.assertEqual(row['semantics']['family'],'hbond')

    def test_true_reversals_merge_but_source_evidence_and_context_survive(self):
        d=atom('A',1,'SER','OG','O',(0,0,0));acc=atom('B',2,'GLU','OE1','O',(2.8,0,0));raw=raw_pair(d,acc,['polar']); second=reverse(raw);second['contact']=['hbond']
        rows=records(report([d,acc],[raw,second]));self.assertEqual(len(rows),1);sem=rows[0]['semantics'];self.assertEqual(sem['evidence']['supportingRecords'],2);self.assertTrue(all('participants' in item and 'direction' in item for item in sem['sourceObservations']))

    def test_model_conformer_and_actual_hydrogen_are_identity_components(self):
        d=atom('A',1,'SER','OG','O',(0,0,0),model_id='2',altloc='B');acc=atom('B',2,'GLU','OE1','O',(2.8,0,0),model_id='2',altloc='B');sem=classify(raw_pair(d,acc,['hbond']),[d,acc])['semantics'];self.assertIn('@model=2@alt=B',sem['identity'])
        left=copy.deepcopy(sem);right=copy.deepcopy(sem);left['geometry']['hydrogen']={'id':'H1'};right['geometry']['hydrogen']={'id':'H2'};self.assertNotEqual(a._semantic_identity(left),a._semantic_identity(right))

    def test_group_charge_retains_neutral_resonance_atom_without_neutral_group(self):
        lys=atom('A',1,'LYS','NZ','N',(0,0,0));x=atom('B',2,'GLU','OE1','O',(3,0,0));y=atom('B',2,'GLU','OE2','O',(3,.8,0))
        rows=[]
        for acceptor,charge in ((x,0),(y,-1)):
            raw=raw_pair(lys,acceptor,['ionic']);raw['bgn']['formal_charge']=1;raw['end']['formal_charge']=charge;rows.append(raw)
        actual=report([lys,x,y],rows)['contacts']['salt_bridges'];self.assertEqual(len(actual),1)
        sem=actual[0]['semantics'];negative=next(p for p in sem['participants'] if p['role']=='negative_site')
        self.assertEqual(negative['site']['kind'],'group');self.assertEqual(negative['charge']['value'],-1);self.assertEqual(negative['charge']['atomValue'],0)
        self.assertNotIn('charge_state_observations_differ',sem['evidence']['ambiguityFlags'])
        changed=next(o for o in sem['sourceObservations'] if o.get('sourceEndpoints'))
        self.assertEqual(changed['sourceEndpoints'][1]['atom']['atomName'],'OE2');self.assertEqual(changed['sourceEndpoints'][1]['charge']['value'],-1)

    def test_ring_coordinates_alone_do_not_support_diagnostic_geometry(self):
        left,right=ring('A',1),ring('B',2,center=(0,0,10));idx=a.ResidueAtomIndex(left+right)
        payload=lambda at:dict(a.residue_payload(at),element=at.element)
        for family in ('pi_pi','aromatic_proximal'):
            sem=a._build_interaction_semantics(raw_pair(left[0],right[0],['aromatic']),{'family':family},payload(left[0]),payload(right[0]),idx)
            self.assertEqual(sem['evidence']['level'],'candidate');self.assertEqual(sem['evidence']['geometrySupport'],'partial')

    def test_neutral_cation_pi_and_incompatible_ligand_hbond_are_rejected(self):
        aromatic=ring('B',2);n=atom('A',1,'LYS','NZ','N',(0,0,4));raw=raw_pair(n,aromatic[0],['cationpi']);raw['bgn']['formal_charge']=0
        self.assertNotEqual(classify(raw,[n]+aromatic)['family'],'pi_cation')
        acceptor=atom('B',3,'GLU','OE1','O',(2.8,0,4));raw=raw_pair(n,acceptor,['hbond']);raw['bgn']['atom_types']=['hbond acceptor'];raw['end']['atom_types']=['hbond acceptor']
        self.assertNotEqual(classify(raw,[n,acceptor])['family'],'hbond')

    def test_nested_chain_remapping_changes_identity_once_and_keeps_source_chains(self):
        row=runtime_fixture()['contacts']['hydrogen_bonds'][0]
        aliases=a._identity_chain_aliases({'A','B'});aliases.canonical_to_auth.update(A='B',B='C')
        mapped=a._remap_contact_record_for_external(row,aliases);sem=mapped['semantics']
        self.assertEqual([p['site']['residue']['chain'] for p in sem['participants']],['B','C'])
        self.assertTrue(sem['participants'][0]['site']['atoms'][0]['id'].startswith('B:'))
        self.assertTrue(sem['geometry']['hydrogen']['id'].startswith('B:'))
        self.assertEqual({p['site']['residue']['chain'] for o in sem['sourceObservations'] for p in o['participants']},{'B','C'})
        self.assertEqual(sem['identity'],a._semantic_identity(sem));self.assertTrue(row['semantics']['participants'][0]['site']['id'].startswith('A:'))

    def test_duplicate_merge_retains_original_provenance_snapshots(self):
        row=runtime_fixture()['contacts']['salt_bridges'][0];other=copy.deepcopy(row)
        other['semantics']['participants'][0]['roleProvenance']=['extra_independent_typing']
        merged=a._merge_semantic_contact_records(row,other)
        self.assertIn('extra_independent_typing',merged['semantics']['participants'][0]['roleProvenance'])
        original=merged['semantics']['sourceObservations'][0]['participants'][0]['roleProvenance']
        self.assertNotIn('extra_independent_typing',original)

    def test_supplied_observation_types_precede_cached_chemistry(self):
        at=atom('A',1,'LIG','N1','N',(0,0,0));idx=a.ResidueAtomIndex([at]);idx.chemical_nodes[('A','1','N1')]={'auth_atom_id':'N1','formal_charge':-1}
        supplied={'auth_atom_id':'N1','formal_charge':1};self.assertEqual(a._semantic_node(dict(a.residue_payload(at),element='N'),supplied,idx)['formal_charge'],1)

    def test_supplied_and_perceived_charge_state_are_distinguished(self):
        sem=runtime_fixture()['contacts']['salt_bridges'][0]['semantics'];self.assertEqual(sem['evidence']['chemicalState'],'input');self.assertNotIn('chemical_state_assumed',sem['evidence']['ambiguityFlags'])
        n=atom('A',1,'QAM','N1','N',(0,0,0));o=atom('B',2,'ACT','O1','O',(3,0,0));raw=raw_pair(n,o,['ionic']);raw['bgn']['formal_charge']=1;raw['end']['formal_charge']=-1
        self.assertEqual(classify(raw,[n,o])['semantics']['evidence']['chemicalState'],'modelled')

    def test_cation_pi_in_ring_plane_is_not_geometrically_supported(self):
        aromatic=ring('B',2);n=atom('A',1,'LYS','NZ','N',(4,0,0));raw=raw_pair(n,aromatic[0],['cationpi']);raw['bgn']['formal_charge']=1
        sem=a._build_interaction_semantics(raw,{'family':'pi_cation'},dict(a.residue_payload(n),element='N'),dict(a.residue_payload(aromatic[0]),element='C'),a.ResidueAtomIndex([n]+aromatic))
        self.assertEqual(sem['evidence']['level'],'candidate');self.assertEqual(next(m['value'] for m in sem['geometry']['measurements'] if m['kind']=='cation_ring_normal_angle'),90)

    def test_canonical_schema_survives_json_and_empty_reports(self):
        result=json.loads(json.dumps(runtime_fixture(),allow_nan=False));self.assertEqual(result['meta']['interactionSemanticsVersion'],1)
        self.assertTrue(all(r['semantics']['version']==1 and r['semantics']['family']==r['type'] for r in records(result)))
        d=atom('A',1,'ALA','CA','C',(0,0,0));acc=atom('B',2,'ALA','CA','C',(20,0,0));self.assertEqual(report([d,acc],[])['meta']['interactionSemanticsVersion'],1)


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--write-fixture':
        target=Path(__file__).resolve().parents[2]/'roami-tests/interaction-semantics-2026-09-20/canonical-fixture.json'; target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(runtime_fixture(),indent=2,allow_nan=False)+'\n');print(target)
    else:
        unittest.main()
