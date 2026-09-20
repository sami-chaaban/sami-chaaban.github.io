"""Base pairing stays residue context in real DNA engine/report output."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from api import analysis as a

FIXTURE=Path(__file__).resolve().parents[2]/'roami-tests/audit-2026-09-20/chemistry/1bna.cif'


@unittest.skipIf(a.InteractionComplex is None or not FIXTURE.exists(),'Arpeggio or deposited DNA fixture unavailable')
class DepositedBasePairAtomSemanticsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=a.analyze_interface(FIXTURE.read_text(),'A','B')
        cls.rows=[r for group in cls.report['contacts'].values() for r in group]

    def matching(self,left,right):
        return [r for r in self.rows if {(r[side]['resName'],r[side]['atom']) for side in ('residueA','residueB')}=={left,right} and r.get('basePair')]

    def test_atom_contacts_never_replace_chemistry_with_base_pair_family(self):
        self.assertEqual(self.report['contacts']['base_pairing'],[])
        self.assertFalse(any(r['type']=='base_pairing' or r['semantics']['family']=='base_pairing' for r in self.rows))
        contextual=[r for r in self.rows if r.get('basePair')]
        self.assertGreater(len(contextual),0)
        self.assertTrue(any(r['basePair']['isCanonicalAtomPattern'] for r in contextual))
        self.assertTrue(all(p['role']!='base' for r in contextual for p in r['semantics']['participants']))

    def test_engine_donor_donor_gc_contact_is_ambiguous_polar_not_hbond(self):
        rows=self.matching(('DC','N3'),('DG','N1'));self.assertGreater(len(rows),0)
        for row in rows:
            self.assertEqual(row['type'],'polar_contact');sem=row['semantics']
            self.assertEqual(sem['evidence']['level'],'ambiguous');self.assertEqual(sem['direction']['certainty'],'ambiguous');self.assertIsNone(sem['direction']['from'])
            c=next(p for p in sem['participants'] if p['site']['residue']['resName']=='DC')
            self.assertEqual(c['role'],'unresolved');self.assertIn('hbond donor',c['atomTypes'])
            self.assertEqual({(v['role'],v['source']) for v in c['roleAlternatives']},{('donor','arpeggio_openbabel_atom_types'),('acceptor','standard_residue_template')})
            self.assertTrue(row['basePair']['isCanonicalAtomPattern']);self.assertIn('nucleobase_typing_template_conflict',sem['evidence']['ambiguityFlags'])

    def test_missing_typed_amino_donor_is_not_reconstructed_from_base_name(self):
        rows=self.matching(('DC','N4'),('DG','O6'));self.assertGreater(len(rows),0)
        for row in rows:
            self.assertEqual(row['type'],'polar_contact');sem=row['semantics']
            c=next(p for p in sem['participants'] if p['site']['residue']['resName']=='DC')
            self.assertNotIn('hbond donor',c['atomTypes']);self.assertEqual(c['role'],'unresolved')
            self.assertIn({'role':'donor','source':'standard_residue_template'},c['roleAlternatives'])
            self.assertEqual(sem['evidence']['level'],'ambiguous')

    def test_compatible_typed_gc_atoms_keep_hbond_roles_and_pair_context(self):
        rows=self.matching(('DC','O2'),('DG','N1'));self.assertGreater(len(rows),0)
        supported=[r for r in rows if r['type']=='hbond'];self.assertGreater(len(supported),0)
        for row in supported:
            self.assertEqual({p['role'] for p in row['semantics']['participants']},{'donor','acceptor'})
            self.assertEqual(row['semantics']['geometry']['distance']['kind'],'donor_acceptor')
            self.assertEqual(row['basePair']['family'],'GC')

    def test_canonical_n2_o2_pair_uses_supplied_types_not_residue_identity(self):
        from test_interaction_semantics import atom, raw_pair, report, records
        donor=atom('A',1,'DG','N2','N',(0,0,0));acceptor=atom('B',2,'DC','O2','O',(2.8,0,0));raw=raw_pair(donor,acceptor,['vdw_clash'])
        raw['bgn']['atom_types']=['hbond donor'];raw['end']['atom_types']=['hbond acceptor']
        result=records(report([donor,acceptor],[raw]))[0]
        self.assertEqual(result['type'],'hbond');self.assertTrue(result['basePair']['isCanonicalAtomPattern'])
        self.assertEqual({p['role'] for p in result['semantics']['participants']},{'donor','acceptor'})
        raw['bgn']['atom_types']=['pos ionisable']
        result=records(report([donor,acceptor],[raw]))[0]
        self.assertEqual(result['type'],'polar_contact');self.assertTrue(result['basePair']['isCanonicalAtomPattern'])
        self.assertEqual(result['semantics']['evidence']['level'],'ambiguous')

    def test_base_pair_mode_filters_context_without_rewriting_atom_family(self):
        filtered=a.filter_contacts_by_mode(self.report['contacts'],'base_pair')
        rows=[r for group in filtered.values() for r in group]
        self.assertEqual({r['semantics']['identity'] for r in rows},{r['semantics']['identity'] for r in self.rows if r.get('basePair')})
        self.assertTrue(all(r['type']==r['semantics']['family'] and r['type']!='base_pairing' for r in rows))
        self.assertTrue(any(r['type']=='hbond' for r in rows));self.assertTrue(any(r['type']=='polar_contact' for r in rows))


if __name__=='__main__':unittest.main()
