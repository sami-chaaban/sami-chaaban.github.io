import sys
import unittest
import os
import re
from pathlib import Path


sys.path.insert(0, "my-site/ppi")
from api import analysis


MINIMAL_AMBIGUOUS_CHAIN_MMCIF = """data_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 S SG . CYS D 1 33 ? 0.000 0.000 0.000 1.00 10.00 33 CYS Y SG 1
HETATM 2 ZN ZN . ZN J 2 . ? 2.330 0.000 0.000 1.00 10.00 501 ZN Y ZN 1
ATOM 3 C CA . GLY E 3 1 ? 5.000 0.000 0.000 1.00 10.00 1 GLY J CA 1
#
"""


class AnalysisRegressionTests(unittest.TestCase):
    def test_metal_complex_contact_is_not_demoted_to_covalent_or_clash(self) -> None:
        aliases = analysis._identity_chain_aliases({"Y"})
        raw_contact = {
            "type": "atom-atom",
            "contact": ["COVALENT", "METAL_COMPLEX"],
            "distance": 2.33,
            "bgn": {
                "auth_asym_id": "Y",
                "auth_seq_id": 33,
                "label_comp_id": "CYS",
                "auth_atom_id": "SG",
                "type_symbol": "S",
            },
            "end": {
                "auth_asym_id": "Y",
                "auth_seq_id": 501,
                "label_comp_id": "ZN",
                "auth_atom_id": "ZN",
                "type_symbol": "ZN",
            },
        }
        residue_a = {"chain": "Y", "resName": "CYS", "seq": "33", "atom": "SG", "element": "S"}
        residue_b = {"chain": "Y", "resName": "ZN", "seq": "501", "atom": "ZN", "element": "ZN"}

        asserted = analysis._assert_interaction(
            raw_contact,
            residue_a,
            residue_b,
            aliases,
            residue_atoms_index={},
        )

        self.assertEqual(asserted.get("family"), "metal_coordination")
        self.assertEqual(asserted.get("confidence"), "high")

    def test_ambiguous_auth_chain_resolves_to_disambiguated_internal_chain_ids(self) -> None:
        chains, aliases = analysis.list_chains(MINIMAL_AMBIGUOUS_CHAIN_MMCIF)

        self.assertIn("Y", chains)
        self.assertEqual(aliases.resolve("Y"), {"D", "Y[J]"})
        self.assertNotIn("D", aliases.label_to_auth)

        zinc_payload = analysis._build_residue_payload_from_arpeggio_partner(
            {
                "auth_asym_id": "Y",
                "auth_seq_id": 501,
                "label_comp_id": "ZN",
                "auth_atom_id": "ZN",
            },
            aliases,
        )
        cysteine_payload = analysis._build_residue_payload_from_arpeggio_partner(
            {
                "auth_asym_id": "Y",
                "auth_seq_id": 33,
                "label_comp_id": "CYS",
                "auth_atom_id": "SG",
            },
            aliases,
        )

        self.assertIsNotNone(zinc_payload)
        self.assertIsNotNone(cysteine_payload)
        self.assertEqual(zinc_payload["chain"], "Y[J]")
        self.assertEqual(cysteine_payload["chain"], "D")

    def test_external_report_payloads_are_remapped_back_to_auth_chain_ids(self) -> None:
        _, aliases = analysis.list_chains(MINIMAL_AMBIGUOUS_CHAIN_MMCIF)
        contacts = {
            "metal_coordination": [
                {
                    "residueA": {"chain": "Y[J]", "resName": "ZN", "seq": "501", "atom": "ZN"},
                    "residueB": {"chain": "D", "resName": "CYS", "seq": "33", "atom": "SG"},
                    "atomKeyA": "Y[J]:501:ZN",
                    "atomKeyB": "D:33:SG",
                    "pairKey": "D:33:SG|Y[J]:501:ZN",
                    "ringKeyA": "Y[J]:501:ring",
                    "ringKeyB": "D:33:ring",
                    "ringPairKey": "D:33:ring|Y[J]:501:ring",
                    "asserted": {
                        "family": "metal_coordination",
                        "ringKeyA": "Y[J]:501:ring",
                        "ringKeyB": "D:33:ring",
                        "ringPairKey": "D:33:ring|Y[J]:501:ring",
                    },
                }
            ]
        }
        per_residue = {
            "Y[J]:501": {
                "chain": "Y[J]",
                "resName": "ZN",
                "seq": "501",
                "metal_coordination": 1,
                "total": 1,
            },
            "D:33": {
                "chain": "D",
                "resName": "CYS",
                "seq": "33",
                "metal_coordination": 1,
                "total": 1,
            },
        }

        (
            chain_a,
            chain_b,
            remapped_contacts,
            remapped_per_residue,
            _buried_fraction,
            remapped_meta,
        ) = analysis._remap_report_to_external_chains(
            chain_a="Y[J]",
            chain_b="D",
            contacts=contacts,
            per_residue=per_residue,
            buried_fraction={"Y[J]": 0.5, "D": 0.5},
            meta={"focusResidue": "Y[J]:501"},
            aliases=aliases,
        )

        self.assertEqual(chain_a, "Y")
        self.assertEqual(chain_b, "Y")
        self.assertEqual(remapped_meta["focusResidue"], "Y:501")
        self.assertIn("Y:501", remapped_per_residue)
        self.assertIn("Y:33", remapped_per_residue)
        remapped_contact = remapped_contacts["metal_coordination"][0]
        self.assertEqual(remapped_contact["residueA"]["chain"], "Y")
        self.assertEqual(remapped_contact["residueB"]["chain"], "Y")
        self.assertEqual(remapped_contact["atomKeyA"], "Y:501:ZN")
        self.assertEqual(remapped_contact["atomKeyB"], "Y:33:SG")
        self.assertEqual(remapped_contact["pairKey"], "Y:33:SG|Y:501:ZN")
        self.assertEqual(remapped_contact["ringKeyA"], "Y:501:ring")
        self.assertEqual(remapped_contact["ringKeyB"], "Y:33:ring")
        self.assertEqual(remapped_contact["ringPairKey"], "Y:33:ring|Y:501:ring")


class MmcifCaseRegressionTests(unittest.TestCase):
    fixture = (Path(__file__).resolve().parents[2] / "scripts/fixtures/mmcif-case.cif").read_text()

    def test_coordinate_header_case_does_not_discard_atoms(self) -> None:
        canonical = analysis.parse_mmcif_atoms(self.fixture)
        self.assertEqual(len(canonical[0]), 7)
        coot_text = self.fixture.replace("_atom_site.Cartn_", "_atom_site.cartn_")
        self.assertEqual(analysis.parse_mmcif_atoms(coot_text), canonical)

    def test_tag_and_keyword_case_preserves_coordinates_names_and_chain_ids(self) -> None:
        canonical = analysis.parse_mmcif_atoms(self.fixture)
        transforms = [str.lower, str.upper, str.swapcase]
        for transform in transforms:
            with self.subTest(case=transform.__name__):
                text = re.sub(r"^_\S+|^loop_$", lambda m: transform(m.group()), self.fixture, flags=re.M)
                parsed = analysis.parse_mmcif_atoms(text)
                self.assertEqual(parsed, canonical)
                self.assertEqual({a.chain_auth for a in parsed[0]}, {"Da", "DA", "rA"})
                self.assertEqual([a.atom_name for a in parsed[0]], ["N", "CA", "CA", "C4'", "O5'", "CA", "CA"])

    def test_uppercase_loop_boundary_does_not_become_atom_data(self) -> None:
        text = re.sub(r"^_\S+|^loop_$", lambda m: m.group().upper(), self.fixture, flags=re.M)
        text += "LOOP_\n_other.id\n_other.name\n1 unrelated\n"
        self.assertEqual(analysis.parse_mmcif_atoms(text), analysis.parse_mmcif_atoms(self.fixture))

    def test_missing_coordinate_column_is_not_replaced_with_another_column(self) -> None:
        text = self.fixture.replace("_atom_site.Cartn_x", "_atom_site.unused_x")
        self.assertEqual(analysis.parse_mmcif_atoms(text)[0], [])

    @unittest.skipUnless(os.environ.get("ROAMI_CIF_REGRESSION_FILE"), "No local regression file specified")
    def test_local_regression_file(self) -> None:
        text = Path(os.environ["ROAMI_CIF_REGRESSION_FILE"]).read_text()
        atoms, aliases = analysis.parse_mmcif_atoms(text)
        self.assertEqual(len(atoms), 8660)
        self.assertEqual(len({a.chain_auth for a in atoms}), 9)
        self.assertEqual({a.atom_name for a in atoms if a.chain_auth == "Da" and a.res_seq == "710"},
                         {"O", "N", "CA", "C", "CB", "CG", "CD1", "CD2"})


if __name__ == "__main__":
    unittest.main()
