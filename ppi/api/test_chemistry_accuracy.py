"""Regression tests for coordinate identity, chemical evidence, and pruning."""
from dataclasses import replace
import gzip
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from api import analysis as a

CAPTURES = Path(__file__).resolve().parents[2] / 'roami-tests/audit-2026-09-20/chemistry'


def atom(chain, seq, residue, name, element, xyz, **kwargs):
    return a.AtomRecord(chain, chain, chain, residue, str(seq), name, element, *xyz, **kwargs)


def raw_pair(left, right, terms, **kwargs):
    def node(at):
        return dict(auth_asym_id=at.chain_id, auth_seq_id=at.res_seq,
                    label_comp_id=at.res_name, auth_atom_id=at.atom_name)
    return dict(type='atom-atom', distance=a.distance(left, right),
                contact=terms, bgn=node(left), end=node(right), **kwargs)


def classify(raw, atoms, **kwargs):
    aliases = a._identity_chain_aliases({at.chain_id for at in atoms})
    left = a._build_residue_payload_from_arpeggio_partner(raw['bgn'], aliases)
    right = a._build_residue_payload_from_arpeggio_partner(raw['end'], aliases)
    return a._assert_interaction(raw, left, right, aliases, a.ResidueAtomIndex(atoms), **kwargs)


class ChemistryAccuracyTests(unittest.TestCase):
    def test_coordinate_elements_restore_halogen_and_metal_chemistry(self):
        oxygen = atom('B', 2, 'ASN', 'OD1', 'O', (3.3, 0, 0))
        for element in ('CL', 'BR'):
            halogen = atom('A', 1, 'LIG', element + '1', element, (0, 0, 0))
            anchor = atom('A', 1, 'LIG', 'C1', 'C', (-1.8, 0, 0))
            self.assertEqual(classify(raw_pair(halogen, oxygen, ['xbond']), [halogen, oxygen, anchor])['family'], 'halogen_bond')
        iron = atom('A', 1, 'HEM', 'FE', 'FE', (0, 0, 0))
        histidine = atom('B', 2, 'HIS', 'NE2', 'N', (2.06, 0, 0))
        self.assertEqual(classify(raw_pair(iron, histidine, ['metal_complex', 'covalent']), [iron, histidine])['family'], 'metal_coordination')

    def test_unknown_ligand_nitrogen_does_not_invent_hbond(self):
        nitrogen = atom('A', 1, 'PYR', 'N1', 'N', (0, 0, 0))
        oxygen = atom('B', 2, 'ASN', 'OD1', 'O', (2.8, 0, 0))
        self.assertNotEqual(classify(raw_pair(nitrogen, oxygen, ['proximal']), [nitrogen, oxygen])['family'], 'hbond')

    def test_ligand_types_distinguish_alcohol_pyridine_amide_and_quaternary_n(self):
        oxygen = atom('B', 2, 'ASN', 'OD1', 'O', (2.8, 0, 0))
        donor = atom('B', 2, 'ASN', 'ND2', 'N', (2.8, 0, 0))
        cases = [('EOH', 'O', 'O', ['hbond donor', 'hbond acceptor'], oxygen, True),
                 ('PYR', 'N1', 'N', ['hbond acceptor'], oxygen, False),
                 ('AMD', 'N1', 'N', ['hbond donor'], donor, False),
                 ('QAM', 'N1', 'N', ['pos ionisable'], donor, False)]
        for residue, name, element, roles, partner, expected in cases:
            with self.subTest(residue=residue):
                ligand = atom('A', 1, residue, name, element, (0, 0, 0))
                raw = raw_pair(ligand, partner, ['hbond', 'polar'], hbond_angle=170)
                raw['bgn']['atom_types'] = roles
                self.assertEqual(classify(raw, [ligand, partner])['family'] == 'hbond', expected)

    @unittest.skipIf(a.InteractionComplex is None or a._openbabel is None, 'Arpeggio/Open Babel unavailable')
    def test_upstream_smarts_validates_ligand_role_fixtures(self):
        from arpeggio.core import config
        ob = a._openbabel
        a.InteractionComplex.address_ambiguities(None)
        for smiles, index, expected in [('CCO', 3, (True, True)), ('n1ccccc1', 1, (False, True)),
                                        ('CC(=O)N', 4, (True, False)), ('C[N+](C)(C)C', 2, (False, False)),
                                        ('[NH4+]', 1, (True, False))]:
            with self.subTest(smiles=smiles):
                converter, molecule = ob.OBConversion(), ob.OBMol()
                converter.SetInFormat('smi')
                self.assertTrue(converter.ReadString(molecule, smiles))
                roles = []
                for role in ('hbond donor', 'hbond acceptor'):
                    for expression in config.ATOM_TYPES[role].values():
                        pattern = ob.OBSmartsPattern()
                        self.assertTrue(pattern.Init(expression))
                        pattern.Match(molecule)
                        if any(index in list(match) for match in pattern.GetMapList()):
                            roles.append(role)
                            break
                donor, acceptor, known = a._partner_hbond_roles({'atom_types': roles}, 'LIG', 'X', 'N')
                self.assertEqual((donor, acceptor), expected)
                self.assertTrue(known)

    def test_explicit_ligand_hbond_is_retained_when_older_export_lacks_roles(self):
        alcohol = atom('A', 1, 'EOH', 'O', 'O', (0, 0, 0))
        oxygen = atom('B', 2, 'ASN', 'OD1', 'O', (2.8, 0, 0))
        result = classify(raw_pair(alcohol, oxygen, ['hbond', 'polar']), [alcohol, oxygen])
        self.assertEqual(result['family'], 'hbond')
        self.assertIn('upstream_hbond_typing_evidence', result['evidence'])

    def test_hydroxyl_axis_neither_rejects_upstream_bond_nor_confirms_candidate(self):
        donor = atom('A', 1, 'SER', 'OG', 'O', (0, 0, 0))
        anchor = atom('A', 1, 'SER', 'CB', 'C', (1.43, 0, 0))
        for degrees in (109.5, 180):
            acceptor = atom('B', 2, 'ASN', 'OD1', 'O', (2.8 * math.cos(math.radians(degrees)), 2.8 * math.sin(math.radians(degrees)), 0))
            for terms in (['hbond', 'polar'], ['polar']):
                result = classify(raw_pair(donor, acceptor, terms), [donor, anchor, acceptor])
                self.assertEqual(result['family'], 'hbond')
                self.assertEqual(result['subtype'], 'hbond_candidate')
                self.assertNotEqual(result['confidence'], 'high')

    def test_rna_hydroxyl_can_donate_to_phosphate_but_ether_cannot(self):
        phosphate = atom('B', 20, 'G', 'OP1', 'O', (2.8, 0, 0))
        for name, allowed in [("O2'", True), ("O4'", False)]:
            sugar = atom('A', 1, 'A', name, 'O', (0, 0, 0))
            result = classify(raw_pair(sugar, phosphate, ['hbond', 'polar'], hbond_angle=170), [sugar, phosphate])
            self.assertEqual(result['family'] == 'hbond', allowed)

    def test_basepair_identity_is_orthogonal_to_hbond_term(self):
        left = atom('A', 1, 'DA', 'N6', 'N', (0, 0, 0))
        right = atom('B', 10, 'DT', 'O4', 'O', (3, 0, 0))
        key = a._unordered_residue_pair_key(a.residue_payload(left), a.residue_payload(right), prefix='basepair_support:')
        stats = {key: dict(supportCount=2, mutualBestMatch=True, coplanaritySupported=True)}
        for term in ('proximal', 'hbond', 'polar'):
            result = classify(raw_pair(left, right, [term]), [left, right], base_pair_pair_stats=stats)
            self.assertTrue(result['basePair']['isCanonicalAtomPattern'])
            self.assertEqual(result['basePair']['family'], 'AT')

    def test_engine_export_retains_chemical_types_and_real_hydrogen_geometry(self):
        class BioAtom:
            def __init__(self, name, element, xyz, roles, hydrogen=()):
                self.name, self.element, self.coord = name, element, xyz
                self.atom_types, self.h_coords = roles, hydrogen
                self.formal_charge, self.vdw_radius = 0, 1.52
            def get_altloc(self):
                return 'B'
        donor = BioAtom('O', 'O', (0, 0, 0), {'hbond donor'}, [(1, 0, 0)])
        acceptor = BioAtom('OD1', 'O', (2.8, 0, 0), {'hbond acceptor'})
        raw = dict(type='atom-atom', bgn={'auth_atom_id': 'O'}, end={'auth_atom_id': 'OD1'})
        engine = SimpleNamespace(atom_contacts=[SimpleNamespace(bgn_atom=donor, end_atom=acceptor)], input_has_hydrogens=False)
        result = a._enrich_arpeggio_contacts([raw], engine)[0]
        self.assertEqual(result['bgn']['atom_types'], ['hbond donor'])
        self.assertEqual(result['hbond_angle'], 180)
        self.assertEqual(result['hbond_geometry_source'], 'arpeggio_generated_hydrogens')
        self.assertEqual(result['hbond_donor_side'], 'A')

    def test_generated_hydrogen_geometry_is_not_labeled_confirmed(self):
        donor = atom('A', 1, 'SER', 'OG', 'O', (0, 0, 0))
        acceptor = atom('B', 2, 'ASN', 'OD1', 'O', (2.8, 0, 0))
        for source in ('arpeggio_generated_hydrogens', 'arpeggio_completed_hydrogens'):
            raw = raw_pair(donor, acceptor, ['hbond'], hbond_angle=180, hbond_geometry_source=source)
            result = classify(raw, [donor, acceptor])
            self.assertEqual(result['subtype'], 'hbond_candidate')
            self.assertNotEqual(result['confidence'], 'high')
            self.assertIn(source, result['evidence'])

    def test_selected_altloc_letters_on_different_residues_are_compatible(self):
        donor = atom('A', 1, 'SER', 'OG', 'O', (0, 0, 0), altloc='A')
        acceptor = atom('B', 2, 'ASN', 'OD1', 'O', (2.8, 0, 0), altloc='B')
        raw = raw_pair(donor, acceptor, ['hbond'])
        raw['bgn']['label_alt_id'], raw['end']['label_alt_id'] = 'A', 'B'
        self.assertEqual(classify(raw, [donor, acceptor])['family'], 'hbond')

    def test_fused_ring_membership_matches_exported_physical_ring(self):
        names = ['N9', 'C8', 'N7', 'C5', 'C4', 'N1', 'C2', 'N3', 'C6']
        atoms = [atom(chain, 1, 'DA', name, name[0], (i % 3, i // 3, z))
                 for chain, z in [('A', 0), ('B', 3.5)] for i, name in enumerate(names)]
        index = a.ResidueAtomIndex(atoms)
        left = dict(chain='A', seq='1', resName='DA', atom='N9,C8,N7,C5,C4')
        right = dict(chain='B', seq='1', resName='DA', atom='N9,C8,N7,C5,C4')
        descriptors = a._residue_ring_descriptors(left, index)
        self.assertEqual(sorted(len(row['atom_names']) for row in descriptors), [5, 6])
        ring = a._resolve_aromatic_ring_site_keys(left, right, index)
        self.assertEqual(set(ring['ringAtomNamesA']), {'N9', 'C8', 'N7', 'C5', 'C4'})
        self.assertEqual(ring['ringAtomNamesA'], ring['ringAtomNamesB'])

    def test_aromatic_output_preserves_distinct_ring_sites_beyond_render_budget(self):
        left, right = dict(chain='A', seq='1', resName='LIG', atom='C1'), dict(chain='B', seq='2', resName='LYS', atom='NZ')
        entries = []
        for i, key in enumerate(('ring1', 'ring2', 'ring3', 'ring1')):
            record = dict(residueA=left, residueB=right, ringPairKey=key, distance=3.5,
                          asserted=dict(family='pi_cation', confidence='medium'))
            entries.append(dict(family='pi_cation', residueA=left, residueB=right, record=record, rank=i))
        selected = a._select_aromatic_records_for_output(entries)
        self.assertEqual({entry['record']['ringPairKey'] for entry in selected}, {'ring1', 'ring2', 'ring3'})
        self.assertEqual(len(selected), 3)

    def test_ring_features_computed_once_per_residue_and_index(self):
        coordinates = [(math.cos(i * math.pi / 3), math.sin(i * math.pi / 3), 0) for i in range(6)]
        atoms = [atom('A', 1, 'LIG', 'C' + str(i), 'C', xyz) for i, xyz in enumerate(coordinates)]
        index = a.ResidueAtomIndex(atoms)
        with patch.object(a, '_residue_ring_descriptors_uncached', wraps=a._residue_ring_descriptors_uncached) as compute:
            first = a._residue_ring_descriptors(a.residue_payload(atoms[0]), index)
            second = a._residue_ring_descriptors(a.residue_payload(atoms[0]), index)
            self.assertIs(first, second)
            self.assertEqual(compute.call_count, 1)


class CoordinatePolicyTests(unittest.TestCase):
    def test_mean_occupancy_selects_one_whole_residue_conformer_and_first_model(self):
        base = atom('A', 1, 'SER', 'CA', 'C', (0, 0, 0))
        atoms = [base,
                 replace(base, atom_name='CB', altloc='A', occupancy=.3),
                 replace(base, atom_name='OG', element='O', altloc='A', occupancy=.3),
                 replace(base, atom_name='CB', altloc='B', occupancy=.8),
                 replace(base, atom_name='OG', element='O', altloc='B', occupancy=.8),
                 replace(base, atom_name='N', element='N', altloc='A', occupancy=.3),
                 replace(base, model_id='2', x=99)]
        selected = a._select_model_conformers(atoms)
        self.assertEqual([(at.atom_name, at.altloc) for at in selected], [('CA', ''), ('CB', 'B'), ('OG', 'B')])
        self.assertEqual({at.model_id for at in selected}, {'1'})

    def test_ties_prefer_a_then_one_then_lexical_and_shared_atom(self):
        base = atom('A', 1, 'SER', 'CA', 'C', (0, 0, 0))
        atoms = [replace(base, altloc=alt) for alt in ('C', '1', 'A')] + [base]
        selected = a._select_model_conformers(atoms)
        self.assertEqual(selected, [base])
        selected = a._select_model_conformers(atoms[:-1])
        self.assertEqual(selected[0].altloc, 'A')

    def test_uploaded_coordinates_override_pdb_id_in_cache_identity(self):
        self.assertNotEqual(a.cache_key('1ABC', 'model1', 'A', 'B', 'all'), a.cache_key('1ABC', 'model2', 'A', 'B', 'all'))

    @unittest.skipIf(a._gemmi is None, 'Gemmi unavailable')
    def test_engine_coordinate_copy_matches_first_model_parser(self):
        text = (Path(__file__).resolve().parents[2] / 'scripts/fixtures/mmcif-case.cif').read_text()
        prepared = a._prepare_arpeggio_mmcif_text(text)
        self.assertEqual(a.parse_mmcif_atoms(text)[0], a.parse_mmcif_atoms(prepared)[0])
        models = a._gemmi.cif.read_string(prepared).sole_block().find_values('_atom_site.pdbx_PDB_model_num')
        self.assertEqual(set(models), {'1'})


@unittest.skipUnless((CAPTURES / '1cll-raw.json.gz').exists(), 'audit fixture captures unavailable')
class DepositedStructureTests(unittest.TestCase):
    def reclassify(self, name, chain_a, chain_b):
        raw = json.loads(gzip.decompress((CAPTURES / f'{name}-raw.json.gz').read_bytes()))
        with patch.object(a, '_run_arpeggio_contacts', return_value=raw):
            return a.analyze_interface((CAPTURES / f'{name}.cif').read_text(), chain_a, chain_b)

    @unittest.skipIf(a.InteractionComplex is None or a._openbabel is None, 'Arpeggio/Open Babel unavailable')
    def test_partial_hydrogen_calmodulin_roles_preserve_atoms_without_invented_geometry(self):
        import tempfile
        ob = a._openbabel
        text = a._prepare_arpeggio_mmcif_text((CAPTURES / '1cll.cif').read_text())
        with tempfile.NamedTemporaryFile('w', suffix='.cif') as handle:
            handle.write(text)
            handle.flush()
            engine = a.InteractionComplex(handle.name)
            before = {at.GetId(): (at.GetAtomicNum(), at.GetFormalCharge(), at.GetX(), at.GetY(), at.GetZ()) for at in ob.OBMolAtomIter(engine.ob_mol)}
            self.assertTrue(engine.input_has_hydrogens)
            self.assertGreater(a._perceive_partial_hydrogen_roles(engine), 0)
            after = {at.GetId(): (at.GetAtomicNum(), at.GetFormalCharge(), at.GetX(), at.GetY(), at.GetZ()) for at in ob.OBMolAtomIter(engine.ob_mol)}
            self.assertTrue(all(after.get(key) == value for key, value in before.items()))
            self.assertEqual(engine._roami_hydrogen_typing_source, 'implicit_hydrogen_valence')
            self.assertEqual(set(after), set(before))
            for name in ('structure_checks', 'address_ambiguities', 'initialize'):
                getattr(engine, name)()
            self.assertGreater(sum('hbond donor' in at.atom_types for at in engine.s_atoms), 100)
            self.assertLessEqual(sum(len(at.h_coords) for at in engine.s_atoms), 1)
            repeated = a.InteractionComplex(handle.name)
            a._perceive_partial_hydrogen_roles(repeated)
            for name in ('structure_checks', 'address_ambiguities', 'initialize'):
                getattr(repeated, name)()
            def fingerprint(complex_):
                return [(at.get_full_id(), tuple(sorted(at.atom_types)),
                         tuple(tuple(float(value) for value in xyz) for xyz in at.h_coords))
                        for at in complex_.s_atoms]
            self.assertEqual(fingerprint(engine), fingerprint(repeated))

    def test_calmodulin_preserves_all_28_coordination_edges_and_four_bidentate_glutamates(self):
        report = self.reclassify('1cll', 'A', 'A')
        metals = report['contacts']['metal_coordination']
        self.assertEqual(len(metals), 28)
        bidentate = [row for row in metals if row['asserted']['denticity'] == 2]
        self.assertEqual(len(bidentate), 8)
        self.assertTrue(all(row['asserted']['coordinationDonorAtoms'] == ['OE1', 'OE2'] for row in bidentate))
        self.assertFalse(any(row['asserted'].get('ambiguousDonorAtoms') for row in metals))
        self.assertIsNone(report['buriedFraction'])
        self.assertIn('contactingResidueFraction', report)

    def test_myoglobin_recovers_heme_iron_coordination(self):
        report = self.reclassify('1a6m', 'A', 'A')
        metal = report['contacts']['metal_coordination']
        pairs = [{row['residueA']['atom'], row['residueB']['atom']} for row in metal]
        self.assertIn({'NE2', 'FE'}, pairs)
        self.assertIn({'O1', 'FE'}, pairs)
        self.assertTrue(all(row['metalElement'] == 'FE' for row in metal))

    def test_dna_has_no_noncanonical_edge_asserted_from_residue_support_alone(self):
        report = self.reclassify('1bna', 'A', 'B')
        rows = [row for group in report['contacts'].values() for row in group]
        noncanonical = [row for row in rows if row['asserted'].get('subtype') == 'base_pair_noncanonical']
        self.assertTrue(all(set(row['arpeggioContact']) & {'HBOND', 'WEAK_HBOND'} or row['arpeggio'].get('hbond_angle', 0) >= 150 for row in noncanonical))
        self.assertTrue(any(row.get('basePair', {}).get('isCanonicalAtomPattern') for row in rows))


if __name__ == '__main__':
    unittest.main()
