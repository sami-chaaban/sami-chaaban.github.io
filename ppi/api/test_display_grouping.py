"""Display grouping acceptance: chemistry identities, locality, evidence and counts."""
import copy
import json
import math
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from api import analysis as a
from api.display_grouping import attach_display_groups, FeatureCatalog
from api.analysis_worker import prepare_delivery
from test_chemistry_accuracy import atom, raw_pair
from test_interaction_semantics import report, records, ring, reverse, runtime_fixture


def hydro_fixture():
    left = [atom('A',1,'LEU','CD1','C',(0,0,0)),atom('A',1,'LEU','CD2','C',(0,1.5,0))]
    right = ring('B',1,center=(0,0,3.8))
    rows = [raw_pair(left[0],right[i],['hydrophobic']) for i in (0,1,2)]
    rows += [raw_pair(left[1],right[1],['hydrophobic']),reverse(rows[0])]
    return report(left+right,rows)


def ligand_chain(atoms):
    nodes = {}
    for i,at in enumerate(atoms):
        nodes[at.atom_name] = {'atom_types':['hydrophobe'], 'bonded_atoms':[
            {'auth_atom_id':atoms[j].atom_name,'type_symbol':'C','same_residue':True,'bond_order':1,'aromatic':False}
            for j in (i-1,i+1) if 0<=j<len(atoms)]}
    return nodes


def grouping_fixture():
    """Controlled engine inputs passed through the production report builder."""
    atoms,rows=[],[]
    left=[atom('A',1,'LEU','CD1','C',(0,0,0)),atom('A',1,'LEU','CD2','C',(0,1.5,0))]
    right=ring('B',1,center=(0,0,3.8)); atoms+=left+right
    rows+=[raw_pair(left[0],right[i],['hydrophobic']) for i in (0,1,2)]
    rows+=[raw_pair(left[1],right[1],['hydrophobic']),reverse(rows[0])]
    lig=[atom('A',10,'LIG',f'C{i+1}','C',(20+(i//2)*10+(i%2)*1.5,0,0)) for i in range(4)]
    partners=[atom('B',10,'LEU','CD1','C',(20,0,3.8)),atom('B',10,'LEU','CD2','C',(30,0,3.8))]
    nodes={**ligand_chain(lig[:2]),**ligand_chain(lig[2:])}; atoms+=lig+partners
    for i,at in enumerate(lig):
        raw=raw_pair(at,partners[i//2],['hydrophobic']); raw['bgn'].update(nodes[at.atom_name]); rows.append(raw)
    positions={'CD2':(0,0),'CE2':(0,1.4),'CG':(-1.33,-0.43),'CD1':(-2.15,0.7),'NE1':(-1.33,1.83),'CZ2':(1.21,2.1),'CH2':(2.42,1.4),'CZ3':(2.42,0),'CE3':(1.21,-0.7)}
    trp=[atom('A',20,'TRP',name,'N' if name=='NE1' else 'C',(xy[0],xy[1]+60,0)) for name,xy in positions.items()]
    counterpart=ring('B',20,center=(0,60,3.5)); atoms+=trp+counterpart
    rows += [raw_pair(next(at for at in trp if at.atom_name==name),counterpart[0],['aromatic']) for name in ('CD1','CZ3')]
    met=atom('A',30,'MET','CG','C',(0,90,0)); phenyl=ring('B',30,center=(0,90,4.06)); atoms += [met]+phenyl
    plane=raw_pair(met,phenyl[0],['carbonpi']); plane.update(type='atom-plane',distance=4.06); plane['end']['auth_atom_id']=','.join(at.atom_name for at in phenyl); rows.append(plane)
    def add_pair(left,right,terms,extra=()):
        atoms.extend([left,right,*extra]); raw=raw_pair(left,right,terms); rows.append(raw); return raw
    add_pair(atom('A',40,'SER','OG','O',(0,120,0)),atom('B',40,'GLU','OE1','O',(2.8,120,0)),['hbond'],[atom('A',40,'SER','HG','H',(.96,120,0))])
    ionic=add_pair(atom('A',50,'LYS','NZ','N',(0,150,0)),atom('B',50,'ASP','OD1','O',(3,150,0)),['ionic'])
    ionic['bgn'].update(formal_charge=1,formal_charge_source='input_formal_charge'); ionic['end'].update(formal_charge=0,formal_charge_source='input_formal_charge')
    second_oxygen=atom('B',50,'ASP','OD2','O',(3,150.8,0)); atoms.append(second_oxygen)
    second_ionic=copy.deepcopy(ionic); second_ionic['end']['auth_atom_id']='OD2'; second_ionic['end']['formal_charge']=-1; second_ionic['distance']=math.sqrt(9+.64); rows.append(second_ionic)
    zn=atom('A',60,'ZN','ZN','ZN',(0,180,0)); oe1=atom('B',60,'GLU','OE1','O',(2.1,180,0)); oe2=atom('B',60,'GLU','OE2','O',(0,182.2,0)); atoms += [zn,oe1,oe2]; rows += [raw_pair(zn,o,['metal_complex']) for o in (oe1,oe2)]
    xb=add_pair(atom('A',70,'LIG','BR1','BR',(0,210,0)),atom('B',70,'ASN','OD1','O',(3.3,210,0)),['xbond'],[atom('A',70,'LIG','C1','C',(-1.8,210,0))])
    xb['bgn']['bonded_atoms']=[{'auth_atom_id':'C1','type_symbol':'C','coordinates':[-1.8,210,0]}]; xb['end']['atom_types']=['hbond acceptor']
    add_pair(atom('A',80,'ALA','CB','C',(0,240,0)),atom('B',80,'ALA','CB','C',(2,240,0)),['vdw_clash'])
    add_pair(atom('A',90,'HIS','ND1','N',(0,270,0)),atom('B',90,'SER','OG','O',(2.8,270,0)),['polar'])
    cat_ring=ring('B',100,center=(0,300,0)); guan=[atom('A',100,'ARG',name,'N',(x,300,4)) for name,x in (('NH1',0),('NH2',.8),('NE',-.8))]; atoms+=cat_ring+guan
    for nitrogen in guan[:2]:
        raw=raw_pair(nitrogen,cat_ring[0],['cationpi']); raw['bgn'].update(formal_charge=1,formal_charge_source='input_formal_charge'); rows.append(raw)
    carbonyl=atom('A',110,'ALA','C','C',(0,330,0)); phenyl=ring('B',110,center=(0,330,4.59)); atoms += [carbonyl]+phenyl
    plane=raw_pair(carbonyl,phenyl[0],['carbonpi']); plane.update(type='atom-plane',distance=4.59); plane['end']['auth_atom_id']=','.join(at.atom_name for at in phenyl); rows.append(plane)
    polar=add_pair(atom('A',120,'ASN','OD1','O',(0,360,0)),atom('B',120,'GLU','OE1','O',(2.9,360,0)),['polar'])
    polar['bgn']['atom_types']=['hbond acceptor']; polar['end']['atom_types']=['hbond acceptor']
    return report(atoms,rows)


class DisplayGroupingTests(unittest.TestCase):
    def test_standard_aromatic_hydrophobic_features_use_physical_membership(self):
        for residue,contact_name,expected in [('PHE','CG',6),('TYR','CG',6),('HIS','CG',5),('TRP','CD1',5),('TRP','CZ3',6)]:
            names=sorted(a._residue_ring_atom_names(residue)); atoms=[atom('A',1,residue,name,'N' if name.startswith('N') else 'C',(math.cos(i),math.sin(i),0)) for i,name in enumerate(names)]
            index=a.ResidueAtomIndex(atoms); payload={'chain':'A','seq':'1','resName':residue,'atom':contact_name}
            participant={'side':'A','role':'contact_atom','site':a._semantic_site(payload,{},index)}
            feature=FeatureCatalog(index,a._identity_chain_aliases({'A'})).for_participant(participant,'hydrophobic')
            with self.subTest(residue=residue,atom=contact_name):
                self.assertEqual(feature['kind'],'ring'); self.assertEqual(len(feature['atoms']),expected)
        payload['atom']='CD2'; participant['site']=a._semantic_site(payload,{},index)
        shared=FeatureCatalog(index,a._identity_chain_aliases({'A'})).for_participant(participant,'hydrophobic')
        self.assertEqual(shared['kind'],'atom'); self.assertIn('shared_fused_ring_atom',shared['ambiguityFlags'])

    def test_two_perceived_ligand_rings_remain_two_hydrophobic_features(self):
        left=atom('A',1,'LEU','CD1','C',(0,0,3.8)); rings=[]; index_nodes={}
        for prefix,cx in (('C',-3),('D',3)):
            members=[atom('B',1,'LIG',f'{prefix}{i+1}','C',(cx+1.4*math.cos(i*math.pi/3),1.4*math.sin(i*math.pi/3),0)) for i in range(6)]
            rings+=members
            for member in members: index_nodes[('B','1',member.atom_name)]={'atom_types':['hydrophobe','aromatic'],'aromatic_rings':[sorted(at.atom_name for at in members)]}
        index=a.ResidueAtomIndex([left]+rings); index.chemical_nodes.update(index_nodes); rows=[]
        for target in (rings[0],rings[9]):
            ra={'chain':'A','seq':'1','resName':'LEU','atom':'CD1'}; rb={'chain':'B','seq':'1','resName':'LIG','atom':target.atom_name}
            raw=raw_pair(left,target,['hydrophobic']); raw['end'].update(index_nodes[('B','1',target.atom_name)])
            sem=a._build_interaction_semantics(raw,{'family':'hydrophobic'},ra,rb,index)
            rows.append({'semantics':sem,'arpeggio':{'type':'atom-atom'},'atomKeyA':a._build_atom_key_from_payload(ra),'atomKeyB':a._build_atom_key_from_payload(rb)})
        result=attach_display_groups({'contacts':{'hydrophobic':rows}},index,a._identity_chain_aliases({'A','B'}))
        self.assertEqual(len(result['displayGroups']['groups']),2)
        ring_features=[feature for feature in result['displayGroups']['features'].values() if feature['kind']=='ring']
        self.assertEqual(len(ring_features),2)
        self.assertTrue(all('arpeggio_openbabel_ring_graph' in feature['provenance'] for feature in ring_features))

    def test_alias_ring_site_ids_share_physical_display_unit_only_within_family(self):
        first=copy.deepcopy(runtime_fixture()['contacts']['pi_pi'][0]); second=copy.deepcopy(first)
        second['semantics']['participants'][0]['site']['id'] += ':opaque-alias'
        second['semantics']['identity']=a._semantic_identity(second['semantics'])
        other=copy.deepcopy(second); other['semantics']['family']='aromatic_packing'; other['semantics']['identity']=a._semantic_identity(other['semantics'])
        atoms=ring('A',7,center=(0,60,0))+ring('B',7,center=(0,60,3.5))
        result=attach_display_groups({'contacts':{'pi_pi':[first,second],'aromatic_packing':[other]}},a.ResidueAtomIndex(atoms),a._identity_chain_aliases({'A','B'}))
        groups=result['displayGroups']['groups']; self.assertEqual(len(groups),2)
        pi=next(g for g in groups if g['family']=='pi_pi')
        self.assertEqual(pi['counts']['canonicalContacts'],2); self.assertEqual(pi['counts']['supportingAtomContacts'],1)
        self.assertEqual(len(result['contacts']['pi_pi']),2)

    def test_salt_group_keeps_both_source_atom_pairs_and_dedupes_reversal(self):
        lys=atom('A',1,'LYS','NZ','N',(0,0,0)); oxygens=[atom('B',2,'GLU',name,'O',(3,y,0)) for name,y in (('OE1',0),('OE2',.8))]; rows=[]
        for oxygen,charge in zip(oxygens,(0,-1)):
            raw=raw_pair(lys,oxygen,['ionic']); raw['bgn']['formal_charge']=1; raw['end']['formal_charge']=charge; rows.append(raw)
        rows.append(reverse(rows[1])); result=report([lys]+oxygens,rows)
        self.assertEqual(len(result['contacts']['salt_bridges']),1)
        group=next(g for g in result['displayGroups']['groups'] if g['family']=='salt_bridge')
        self.assertEqual(group['counts'],{'canonicalContacts':1,'supportingAtomContacts':2,'chemicalUnits':1,'displayObjects':1})
        self.assertEqual({s['atomB']['atomName'] for s in group['supportingAtomContacts']},{'OE1','OE2'})
        self.assertEqual(sorted(s['distance'] for s in group['supportingAtomContacts']),[3.0,round(math.sqrt(9+.64),6)])

    def test_guanidinium_cation_ring_is_one_unit_with_all_atomic_observations(self):
        ring_atoms=ring('B',1)
        nitrogen=[atom('A',1,'ARG','NH1','N',(0,0,4)),atom('A',1,'ARG','NH2','N',(.8,0,4)),atom('A',1,'ARG','NE','N',(-.8,0,4))]
        rows=[]
        for at in nitrogen[:2]:
            raw=raw_pair(at,ring_atoms[0],['cationpi']); raw['bgn'].update(formal_charge=1,formal_charge_source='input_formal_charge'); rows.append(raw)
        result=report(nitrogen+ring_atoms,rows); groups=[g for g in result['displayGroups']['groups'] if g['family']=='pi_cation']
        self.assertEqual(len(result['contacts']['pi_cation']),2); self.assertEqual(len(groups),1)
        self.assertEqual(groups[0]['counts'],{'canonicalContacts':2,'supportingAtomContacts':2,'chemicalUnits':1,'displayObjects':1})
        feature=result['displayGroups']['features'][groups[0]['featureIds'][0]]
        self.assertEqual(feature['kind'],'group'); self.assertEqual({at['atomName'] for at in feature['atoms']},{'NH1','NH2','NE'})
        self.assertNotIn('renderGeometry',groups[0])

    def test_distinct_cations_on_same_ring_remain_distinct(self):
        rings=ring('B',1); left=[atom('A',seq,'LYS','NZ','N',(offset,0,4)) for seq,offset in ((1,0),(2,.8))]; rows=[]
        for at in left:
            raw=raw_pair(at,rings[0],['cationpi']); raw['bgn']['formal_charge']=1; rows.append(raw)
        result=report(left+rings,rows)
        self.assertEqual(len([g for g in result['displayGroups']['groups'] if g['family']=='pi_cation']),2)

    def test_native_ring_observation_joins_unique_atom_support_patch_losslessly(self):
        left=atom('A',1,'LEU','CD1','C',(0,0,0)); rings=ring('B',1,center=(0,0,3.8)); index=a.ResidueAtomIndex([left]+rings)
        result=report([left]+rings,[raw_pair(left,rings[i],['hydrophobic']) for i in (0,1,2)])
        raw=raw_pair(left,rings[0],['carbonpi']); raw.update(type='atom-plane',distance=3.8); raw['end']['auth_atom_id']=[at.atom_name for at in rings]
        ra={'chain':'A','seq':'1','resName':'LEU','atom':'CD1'}; rb={'chain':'B','seq':'1','resName':'PHE','atom':'CD1'}
        sem=a._build_interaction_semantics(raw,{'family':'hydrophobic'},ra,rb,index)
        # Real mmCIF remapping preserves an opaque ring hash from internal chain
        # IDs while support atoms already carry author-chain IDs.
        sem['participants'][1]['site']['id']='B:1:PHE:phenyl:engine_internal_chain_hash@model=1@alt=.'
        sem['identity']=a._semantic_identity(sem)
        native={'semantics':sem,'arpeggio':{'type':'atom-plane'},'atomKeyA':a._build_atom_key_from_payload(ra),'atomKeyB':a._build_atom_key_from_payload(rb)}
        result['contacts']['hydrophobic'].append(native)
        attach_display_groups(result,index,a._identity_chain_aliases({'A','B'}))
        groups=result['displayGroups']['groups']; self.assertEqual(len(groups),1)
        self.assertEqual(groups[0]['counts']['canonicalContacts'],4); self.assertIsNone(groups[0]['counts']['supportingAtomContacts'])
        self.assertEqual(len(groups[0]['supportingAtomContacts']),3)
        ring_feature=result['displayGroups']['features'][groups[0]['featureIds'][1]]
        self.assertIn('engine_internal_chain_hash',ring_feature['id'])
        self.assertGreaterEqual(len(ring_feature['sourceSiteIds']),2)
        self.assertEqual(sem['geometry']['distance']['kind'],'atom_ring_centroid')

    def test_feature_observation_matching_two_patches_remains_ambiguous(self):
        left=[atom('A',1,'LEU',name,'C',(x,0,0)) for name,x in (('CD1',-2),('CD2',2),('CB',0))]
        rings=ring('B',1,center=(0,0,3.8)); index=a.ResidueAtomIndex(left+rings); rows=[]
        for i,at in enumerate(left):
            ra={'chain':'A','seq':'1','resName':'LEU','atom':at.atom_name}; rb={'chain':'B','seq':'1','resName':'PHE','atom':'CG'}
            raw=raw_pair(at,rings[0],['hydrophobic'])
            if i==2:
                raw['type']='atom-plane'; raw['end']['auth_atom_id']=[at.atom_name for at in rings]
            sem=a._build_interaction_semantics(raw,{'family':'hydrophobic'},ra,rb,index)
            rows.append({'semantics':sem,'arpeggio':{'type':raw['type']},'atomKeyA':a._build_atom_key_from_payload(ra),'atomKeyB':a._build_atom_key_from_payload(rb)})
        result=attach_display_groups({'contacts':{'hydrophobic':rows}},index,a._identity_chain_aliases({'A','B'}))
        self.assertEqual(len(result['displayGroups']['groups']),3)
        self.assertEqual(sum('source_feature_matches_multiple_patches' in group['ambiguityFlags'] for group in result['displayGroups']['groups']),1)

    def test_same_pair_stronger_family_retains_raw_hydrophobic_term(self):
        left=atom('A',1,'ALA','CB','C',(0,0,0)); right=atom('B',1,'VAL','CG1','C',(2,0,0))
        result=report([left,right],[raw_pair(left,right,['hydrophobic','vdw_clash'])])
        self.assertEqual(len(records(result)),1); self.assertEqual(records(result)[0]['type'],'clash')
        self.assertIn('hydrophobic',[term.lower() for term in records(result)[0]['arpeggio']['terms']])
        self.assertEqual(result['displayGroups']['groups'][0]['family'],'clash')

    def test_native_atom_plane_preserves_ring_metric_in_every_family(self):
        left=atom('A',145,'MET','CG','C',(0,0,0)); right=ring('B',141,center=(0,0,4.06)); index=a.ResidueAtomIndex([left]+right)
        raw=raw_pair(left,right[0],['carbonpi']); raw['type']='atom-plane'; raw['distance']=4.06
        raw['end']['auth_atom_id']=','.join(at.atom_name for at in right)
        ra={'chain':'A','seq':'145','resName':'MET','atom':'CG'}; rb={'chain':'B','seq':'141','resName':'PHE','atom':'CD1'}
        for family in ('hydrophobic','packing_contact','other'):
            sem=a._build_interaction_semantics(raw,{'family':family},ra,rb,index)
            self.assertEqual(sem['family'],family); self.assertEqual(sem['participants'][1]['site']['kind'],'ring')
            self.assertEqual(len(sem['participants'][1]['site']['atoms']),6)
            self.assertEqual(sem['geometry']['distance']['kind'],'atom_ring_centroid'); self.assertAlmostEqual(sem['geometry']['distance']['value'],4.06)
            record={'semantics':sem,'arpeggio':{'type':'atom-plane'},'atomKeyA':a._build_atom_key_from_payload(ra),'atomKeyB':a._build_atom_key_from_payload(rb)}
            result=attach_display_groups({'contacts':{'other':[record]}},index,a._identity_chain_aliases({'A','B'}))
            group=result['displayGroups']['groups'][0]
            self.assertIsNone(group['counts']['supportingAtomContacts']); self.assertEqual(group['supportingAtomContacts'],[])
            self.assertEqual(result['displayGroups']['features'][group['featureIds'][1]]['kind'],'ring')

    def test_one_carbon_three_ring_atoms_one_unit_preserves_all(self):
        left=atom('A',1,'LEU','CD1','C',(0,0,0)); right=ring('B',1,center=(0,0,3.8))
        result=report([left]+right,[raw_pair(left,right[i],['hydrophobic']) for i in (0,1,2)])
        groups=result['displayGroups']['groups']; self.assertEqual(len(groups),1)
        self.assertEqual(groups[0]['counts'],{'canonicalContacts':3,'supportingAtomContacts':3,'chemicalUnits':1,'displayObjects':1})
        self.assertEqual(len(records(result)),3); self.assertEqual(len(groups[0]['supportingAtomContacts']),3)
        self.assertTrue(groups[0]['renderGeometry']['schematic'])

    def test_multiple_carbons_ring_and_reverse_keep_unique_support(self):
        result=hydro_fixture(); groups=result['displayGroups']['groups']
        self.assertEqual(len(groups),1); self.assertEqual(groups[0]['counts']['supportingAtomContacts'],4)
        self.assertEqual(len(records(result)),4)
        self.assertEqual(sum(r['semantics']['evidence'].get('supportingRecords',1) for r in records(result)),5)
        for support in groups[0]['supportingAtomContacts']:
            self.assertEqual(support['atomA']['residue']['chain'],'A'); self.assertEqual(support['atomB']['residue']['chain'],'B')
            self.assertAlmostEqual(support['distance'],math.dist(support['atomA']['coordinates'],support['atomB']['coordinates']),places=5)

    def test_ring_identity_is_reused_from_canonical_aromatic_site(self):
        left=atom('A',1,'LEU','CD1','C',(0,0,0)); right=ring('B',1,center=(0,0,3.8))
        result=report([left]+right,[raw_pair(left,right[0],['hydrophobic'])])
        expected=a._semantic_ring_site({'chain':'B','seq':'1','resName':'PHE','atom':'CG'},{},a.ResidueAtomIndex(right))['id']
        self.assertIn(expected,result['displayGroups']['features'])

    def test_separate_ligand_nonpolar_groups_do_not_merge_across_linker(self):
        left=[atom('A',1,'LIG',f'C{i+1}','C',(i*1.5,0,0)) for i in range(2)]
        left += [atom('A',1,'LIG',f'C{i+3}','C',(10+i*1.5,0,0)) for i in range(2)]
        right=[atom('B',1,'LEU','CD1','C',(0,0,3.8)),atom('B',1,'LEU','CD2','C',(10,0,3.8))]
        nodes={**ligand_chain(left[:2]),**ligand_chain(left[2:])}
        rows=[]
        for at,partner in zip((left[0],left[2]),right):
            raw=raw_pair(at,partner,['hydrophobic']); raw['bgn'].update(nodes[at.atom_name]); rows.append(raw)
        # Graph neighbors not themselves contacts are available from all raw rows.
        for at,partner in zip((left[1],left[3]),right):
            raw=raw_pair(at,partner,['hydrophobic']); raw['bgn'].update(nodes[at.atom_name]); rows.append(raw)
        result=report(left+right,rows)
        self.assertEqual(len(result['displayGroups']['groups']),2)
        self.assertEqual(sorted(g['counts']['supportingAtomContacts'] for g in result['displayGroups']['groups']),[2,2])

    def test_long_connected_ligand_splits_complete_link_patches(self):
        left=[atom('A',1,'LIG',f'C{i}','C',(i*1.5,0,0)) for i in range(7)]
        right=[atom('B',1,'LIG',f'D{i}','C',(i*1.5,0,3.8)) for i in range(7)]
        ln,rn=ligand_chain(left),ligand_chain(right)
        rows=[]
        for at,bt in zip(left,right):
            raw=raw_pair(at,bt,['hydrophobic']); raw['bgn'].update(ln[at.atom_name]); raw['end'].update(rn[bt.atom_name]); rows.append(raw)
        result=report(left+right,rows); groups=result['displayGroups']['groups']
        self.assertGreater(len(groups),1)
        self.assertEqual(sum(g['counts']['supportingAtomContacts'] for g in groups),7)
        centers=[]
        for group in groups:
            centers.append(group['renderGeometry']['points'][0])
            for side in ('atomA','atomB'):
                xyz=[s[side]['coordinates'] for s in group['supportingAtomContacts']]
                self.assertTrue(all(math.dist(p,q)<=3.5 for p in xyz for q in xyz))
        self.assertEqual(len({tuple(center) for center in centers}),len(centers))

    def test_crossing_connectors_cannot_join_by_midpoint(self):
        left=[atom('A',1,'LEU','CD1','C',(-3,0,0)),atom('A',1,'LEU','CD2','C',(3,0,0))]
        right=[atom('B',1,'LEU','CD1','C',(3,0,0)),atom('B',1,'LEU','CD2','C',(-3,0,0))]
        # Use direct canonical semantic builder to test policy independent of
        # engine distance cutoffs (these long lines share exactly one midpoint).
        source=hydro_fixture(); templates=records(source)[:2]
        index=a.ResidueAtomIndex(left+right); rows=[]
        for i in range(2):
            row=copy.deepcopy(templates[i]); ra={'chain':'A','seq':'1','resName':'LEU','atom':left[i].atom_name}; rb={'chain':'B','seq':'1','resName':'LEU','atom':right[i].atom_name}
            raw=raw_pair(left[i],right[i],['hydrophobic']); row['semantics']=a._build_interaction_semantics(raw,{'family':'hydrophobic'},ra,rb,index)
            row['atomKeyA']=a._build_atom_key_from_payload(ra); row['atomKeyB']=a._build_atom_key_from_payload(rb); rows.append(row)
        result=attach_display_groups({'contacts':{'hydrophobic':rows}},index,a._identity_chain_aliases({'A','B'}))
        self.assertEqual(len(result['displayGroups']['groups']),2)

    def test_uncertain_ligand_topology_remains_atom_features(self):
        left=[atom('A',1,'LIG','C1','C',(0,0,0)),atom('A',1,'LIG','C2','C',(1.5,0,0))]; right=[atom('B',1,'LEU','CD1','C',(0,0,3.8))]
        rows=[raw_pair(at,right[0],['hydrophobic']) for at in left]
        for row in rows: row['bgn']['atom_types']=['hydrophobe']
        result=report(left+right,rows)
        self.assertEqual(len(result['displayGroups']['groups']),2)
        self.assertTrue(all('nonpolar_feature_unresolved' in g['ambiguityFlags'] for g in result['displayGroups']['groups']))

    def test_fused_physical_rings_remain_separate(self):
        result=runtime_fixture(); pi=[g for g in result['displayGroups']['groups'] if g['family']=='pi_pi' and any(result['displayGroups']['features'][f]['residue']['resName']=='TRP' for f in g['featureIds'])]
        self.assertEqual(len(pi),2)
        ring_ids={f for g in pi for f in g['featureIds'] if result['displayGroups']['features'][f]['residue']['resName']=='TRP'}
        self.assertEqual(len(ring_ids),2)

    def test_grouping_never_changes_evidence_or_strength(self):
        result=hydro_fixture()
        before=[{k:copy.deepcopy(v) for k,v in r.items() if k!='displayGrouping'} for r in records(result)]
        atoms=[atom('A',1,'LEU','CD1','C',(0,0,0)),atom('A',1,'LEU','CD2','C',(0,1.5,0))]+ring('B',1,center=(0,0,3.8))
        attach_display_groups(result,a.ResidueAtomIndex(atoms),a._identity_chain_aliases({'A','B'}))
        self.assertEqual(before,[{k:v for k,v in r.items() if k!='displayGrouping'} for r in records(result)])
        self.assertTrue(all('confidence' not in g and 'energy' not in g and 'strength' not in g for g in result['displayGroups']['groups']))

    def test_edge_families_stay_distinct_and_base_context_not_extra_group(self):
        result=runtime_fixture(); groups=result['displayGroups']['groups']
        for family in ('hbond','metal_coordination','halogen_bond','clash','polar_contact'):
            grouped=[g for g in groups if g['family']==family]
            actual=[r for r in records(result) if r['semantics']['family']==family]
            self.assertEqual(len(grouped),len(actual))
            self.assertTrue(all(g['counts']['canonicalContacts']==1 for g in grouped))
        self.assertEqual(sum(g['counts']['canonicalContacts'] for g in groups),len(records(result)))
        self.assertFalse(any(g['family']=='base_pairing' for g in groups))

    def test_feature_only_ring_observation_does_not_invent_atom_support(self):
        result=runtime_fixture(); row=copy.deepcopy(result['contacts']['pi_pi'][0]); row['arpeggio']['type']='ring-ring'; row['semantics'].pop('sourceObservations',None)
        row['semantics'].pop('sourceEndpoints',None); row['atomKeyA']='A:7:PHE:CG,CD1,CE1,CZ,CE2,CD2'; row['atomKeyB']='B:7:PHE:CG,CD1,CE1,CZ,CE2,CD2'
        out=attach_display_groups({'contacts':{'pi_pi':[row]}},a.ResidueAtomIndex([]),a._identity_chain_aliases({'A','B'}))
        group=out['displayGroups']['groups'][0]
        self.assertIsNone(group['counts']['supportingAtomContacts']); self.assertEqual(group['supportingAtomContacts'],[])
        self.assertIsNone(group['geometry']['closestAtomDistance']); self.assertIn('atom_support_not_exported',group['ambiguityFlags'])

    def test_diagnostic_delivery_drops_groups_without_mutating_canonical(self):
        result=hydro_fixture(); row=copy.deepcopy(records(result)[0]); row['debugOnly']=True; row['semantics']['identity']+=':diagnostic'; result['contacts']['hydrophobic'].append(row)
        atoms=[atom('A',1,'LEU','CD1','C',(0,0,0)),atom('A',1,'LEU','CD2','C',(0,1.5,0))]+ring('B',1,center=(0,0,3.8))
        attach_display_groups(result,a.ResidueAtomIndex(atoms),a._identity_chain_aliases({'A','B'}))
        delivery=prepare_delivery(result); canonical=json.loads(delivery.canonical); display=json.loads(delivery.display)
        self.assertEqual(len(canonical['displayGroups']['groups']),2); self.assertEqual(len(display['displayGroups']['groups']),1)
        self.assertFalse(display['displayGroups']['groups'][0]['debugOnly'])

    def test_empty_report_has_grouping_version(self):
        left=atom('A',1,'ALA','CB','C',(0,0,0)); right=atom('B',1,'ALA','CB','C',(20,0,0)); result=report([left,right],[])
        self.assertEqual(result['displayGroups']['groups'],[]); self.assertEqual(result['meta']['displayGroupsVersion'],1)

    def test_order_and_reversal_do_not_change_group_or_support_ids(self):
        left=atom('A',1,'LEU','CD1','C',(0,0,0)); right=ring('B',1,center=(0,0,3.8)); rows=[raw_pair(left,right[i],['hydrophobic']) for i in (0,1,2)]
        first=report([left]+right,rows); second=report([left]+right,[reverse(row) for row in rows[::-1]])
        g1,g2=first['displayGroups']['groups'][0],second['displayGroups']['groups'][0]
        self.assertEqual(g1['id'],g2['id']); self.assertEqual({s['id'] for s in g1['supportingAtomContacts']},{s['id'] for s in g2['supportingAtomContacts']})


if __name__=='__main__':
    if '--write-fixture' in sys.argv:
        path=Path(__file__).resolve().parents[2]/'roami-tests/display-grouping-2026-09-20'
        path.mkdir(parents=True,exist_ok=True)
        (path/'canonical-fixture.json').write_text(json.dumps(grouping_fixture(),indent=2))
        print(path/'canonical-fixture.json')
    else:
        unittest.main()
