"""Narrative consumers must preserve roles without inventing chemistry."""
import copy
import json
import unittest
from unittest.mock import patch

from api.analysis_worker import prepare_delivery
from api.explain import contact_counts, describe_contact, explain_report


def participant(side, role, name, atom, *, kind="atom", atoms=None):
    residue = {"chain": side, "seq": "423" if side == "A" else "157", "resName": name}
    names = atoms or [atom]
    return {
        "side": side,
        "role": role,
        "roleProvenance": ["arpeggio_openbabel_atom_types"],
        "site": {
            "kind": kind, "id": f"{side}:{name}:{','.join(names)}", "residue": residue,
            "atoms": [{"id": f"{side}:{name}:{value}", "atomName": value, "element": "O", "modelId": "1", "altloc": ""} for value in names],
        },
    }


def contact(family="hbond", *, roles=("donor", "acceptor"), directionality="directional", distance_kind="donor_acceptor", level="candidate"):
    partners = [participant("A", roles[0], "SER", "OG"), participant("B", roles[1], "GLU", "OE1")]
    return {
        "residueA": {**partners[0]["site"]["residue"], "atom": "OG"},
        "residueB": {**partners[1]["site"]["residue"], "atom": "OE1"},
        "asserted": {"family": family},
        "semantics": {
            "version": 1, "family": family, "directionality": directionality,
            "participants": partners,
            "direction": {"from": "A" if directionality == "directional" else None, "to": "B" if directionality == "directional" else None, "certainty": "inferred" if directionality == "directional" else "not_applicable"},
            "geometry": {"distance": {"value": 3.33, "kind": distance_kind, "unit": "angstrom", "source": "model_coordinates"}, "measurements": []},
            "evidence": {"level": level, "chemicalCompatibility": "supported", "geometrySupport": "partial", "chemicalState": "assumed", "ambiguityFlags": [], "biologicalInterpretation": "not_evaluated"},
        },
    }


class RoleAwareExplanationTests(unittest.TestCase):
    def test_candidate_uses_preserved_roles_and_heavy_atom_metric(self):
        text = describe_contact(contact())
        self.assertIn("SER A:423 OG donor → GLU B:157 OE1 acceptor", text)
        self.assertIn("(inferred direction)", text)
        self.assertIn("donor–acceptor heavy-atom distance 3.33 Å", text)
        self.assertIn("(candidate)", text)
        self.assertIn("role sources: arpeggio openbabel atom types", text)
        self.assertNotIn("confirmed", text)
        self.assertNotIn("stabilizing", text)

    def test_reversed_endpoint_order_keeps_donor_first(self):
        row = contact()
        for p in row["semantics"]["participants"]:
            p["side"] = "B" if p["side"] == "A" else "A"
        row["semantics"]["direction"].update({"from": "B", "to": "A"})
        text = describe_contact(row)
        self.assertIn("SER A:423 OG donor → GLU B:157 OE1 acceptor", text)

    def test_ambiguous_histidine_chemistry_does_not_get_an_arrow(self):
        row = contact(roles=("donor_or_acceptor", "donor_or_acceptor"), level="ambiguous")
        row["semantics"]["participants"][0]["site"]["residue"]["resName"] = "HIS"
        row["semantics"]["direction"].update({"from": None, "to": None, "certainty": "ambiguous"})
        row["semantics"]["evidence"]["ambiguityFlags"] = ["histidine_tautomer_unknown"]
        text = describe_contact(row)
        self.assertNotIn("→", text)
        self.assertIn("direction unresolved", text)
        self.assertIn("histidine tautomer unknown", text)

    def test_invalid_or_unresolved_direction_metadata_cannot_force_arrow(self):
        row = contact(roles=("unresolved", "unresolved"))
        self.assertNotIn("→", describe_contact(row))
        row = contact()
        row["semantics"]["direction"]["to"] = "A"
        self.assertNotIn("→", describe_contact(row))

    def test_supplied_vs_modeled_hydrogens_have_distinct_provenance(self):
        for source, level in (("input_hydrogens", "geometrically_supported"), ("arpeggio_generated_hydrogens", "candidate")):
            with self.subTest(source=source):
                row = contact(level=level)
                row["semantics"]["geometry"]["hydrogen"] = {"id": "H:1", "atomName": "HG", "source": source, "inferred": source != "input_hydrogens"}
                row["semantics"]["geometry"]["measurements"] = [{"kind": "donor_h_acceptor_angle", "value": 165, "unit": "degree", "source": source}]
                text = describe_contact(row)
                self.assertIn(source.replace("_", " "), text)
                self.assertIn("donor–H–acceptor angle 165 °", text)
                self.assertIn(level.replace("_", " "), text)
                self.assertNotIn("confirmed", text)

    def test_halogen_names_donor_anchor_and_actual_angle(self):
        row = contact("halogen_bond", roles=("halogen_donor", "acceptor"), level="geometrically_supported")
        row["semantics"]["geometry"]["donorAnchor"] = {"atomName": "C1", "source": "arpeggio_bond_graph"}
        row["semantics"]["geometry"]["measurements"] = [{"kind": "halogen_angle", "value": 160, "unit": "degree", "source": "model_coordinates"}]
        text = describe_contact(row)
        self.assertIn("halogen donor →", text)
        self.assertIn("donor anchor C1 (source: arpeggio bond graph)", text)
        self.assertIn("halogen···acceptor angle 160 °", text)

    def test_metal_roles_do_not_invent_oxidation_state_or_direction(self):
        row = contact("metal_coordination", roles=("metal_center", "coordinating_atom"), directionality="role_specific", distance_kind="metal_donor")
        row["semantics"]["participants"][0] = participant("A", "metal_center", "ZN", "ZN", kind="metal")
        text = describe_contact(row)
        self.assertIn("ZN metal center —", text)
        self.assertIn("OE1 coordinating atom", text)
        self.assertIn("metal–donor distance", text)
        self.assertNotIn("→", text)
        self.assertNotIn("2+", text)
        self.assertNotIn("charge", text)

    def test_charge_provenance_and_zero_value_are_preserved(self):
        row = contact("salt_bridge", roles=("positive_site", "negative_site"), directionality="role_specific", distance_kind="charge_site")
        row["semantics"]["participants"][0]["charge"] = {"value": 1, "source": "standard_residue_charge_template", "inferred": True}
        row["semantics"]["participants"][1]["charge"] = {"value": -1, "source": "arpeggio_openbabel_formal_charge", "inferred": True}
        text = describe_contact(row)
        self.assertIn("positive site (charge +1, inferred", text)
        self.assertIn("negative site (charge -1, inferred", text)
        self.assertIn("↔", text)
        self.assertIn("charged-site separation", text)
        row["semantics"]["participants"][0]["charge"]["value"] = 0
        self.assertIn("charge 0", describe_contact(row))

    def test_charge_group_names_the_actual_distance_endpoint(self):
        row = contact("salt_bridge", roles=("positive_site", "negative_site"), directionality="role_specific", distance_kind="charge_site")
        negative = participant("B", "negative_site", "GLU", "OE1", kind="group", atoms=["OE1", "OE2"])
        negative["site"]["contactAtom"] = {"atomName": "OE1", "id": "B:GLU:OE1"}
        row["semantics"]["participants"][1] = negative
        text = describe_contact(row)
        self.assertIn("[OE1, OE2] group, negative site (distance endpoint: OE1)", text)
        self.assertIn("charged-site separation", text)
        self.assertNotIn("centroid", text)

    def test_inferred_charge_sign_does_not_invent_a_formal_magnitude(self):
        row = contact("salt_bridge", roles=("positive_site", "negative_site"), directionality="role_specific", distance_kind="charge_site")
        row["semantics"]["participants"][0]["charge"] = {"value": None, "sign": 1, "source": "standard_residue_charge_template", "inferred": True}
        text = describe_contact(row)
        self.assertIn("positive charge sign inferred, magnitude unspecified; standard residue charge template", text)
        self.assertNotIn("charge +1", text)

    def test_ring_identity_and_distance_are_not_replaced_by_arbitrary_atoms(self):
        row = contact("pi_cation", roles=("cation", "aromatic_ring"), directionality="role_specific", distance_kind="cation_centroid")
        row["semantics"]["participants"][1] = participant("B", "aromatic_ring", "TRP", "CD2", kind="ring", atoms=["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"])
        text = describe_contact(row)
        self.assertIn("[CD2, CE2, CE3, CZ2, CZ3, CH2] aromatic ring", text)
        self.assertIn("cation–ring centroid distance", text)
        self.assertNotIn("→", text)
        row["semantics"]["family"] = "pi_pi"
        row["semantics"]["participants"][0] = participant("A", "aromatic_ring", "TRP", "CG", kind="ring", atoms=["CG", "CD1", "NE1", "CE2", "CD2"])
        row["semantics"]["geometry"]["distance"]["kind"] = "ring_centroid"
        self.assertIn("[CG, CD1, NE1, CE2, CD2] aromatic ring ↔", describe_contact(row))

    def test_site_charge_and_bond_association_provenance_remain_distinct(self):
        row = contact("pi_cation", roles=("cation", "aromatic_ring"), directionality="role_specific", distance_kind="cation_centroid")
        row["semantics"]["participants"][0]["charge"] = {"value": 1, "atomValue": 0, "scope": "group", "source": "formal_charge_sum_of_site_atoms", "inferred": True}
        row["semantics"]["participants"][1]["site"]["provenance"] = "coordinate_ring_inference"
        text = describe_contact(row)
        self.assertIn("group charge +1, inferred; formal charge sum of site atoms", text)
        self.assertIn("distance-endpoint atomic charge 0", text)
        self.assertIn("site sources: coordinate ring inference", text)
        hbond = contact()
        hbond["semantics"]["geometry"]["hydrogen"] = {"atomName": "HG", "source": "input_hydrogens", "bondSource": "coordinate_bond_inference"}
        self.assertIn("hydrogen HG (source: input hydrogens; donor association: coordinate bond inference)", describe_contact(hbond))

    def test_symmetric_contacts_and_clash_overlap_are_not_directional(self):
        for family in ("hydrophobic", "clash"):
            row = contact(family, roles=("contact_atom", "contact_atom"), directionality="symmetric", distance_kind="atom_pair")
            if family == "clash":
                row["semantics"]["geometry"]["measurements"] = [{"kind": "vdw_overlap", "value": 0.8, "unit": "angstrom", "source": "vdw_radii"}]
            text = describe_contact(row)
            self.assertIn("↔", text)
            self.assertNotIn("→", text)
            if family == "clash":
                self.assertIn("van der Waals overlap 0.8 Å", text)

    def test_missing_measurement_metadata_never_relabels_legacy_distance(self):
        row = contact("pi_pi", directionality="symmetric")
        row["distance"] = 4.2
        row["semantics"]["geometry"] = {}
        text = describe_contact(row)
        self.assertIn("measurement metadata unavailable", text)
        self.assertNotIn("4.2", text)
        self.assertNotIn("centroid distance", text)

    def test_legacy_record_keeps_count_but_does_not_invent_roles(self):
        row = contact()
        del row["semantics"]
        row["arpeggio"] = {"hbondDonorSide": "A"}
        text = describe_contact(row, "hydrogen_bonds")
        self.assertIn("canonical role and measurement metadata unavailable", text)
        self.assertNotIn("→", text)
        self.assertNotIn("OG donor", text)
        self.assertEqual(contact_counts({"hydrogen_bonds": [row]}), {"hydrogen_bonds": 1})

    def test_report_counts_semantic_family_and_omits_diagnostic_examples(self):
        row = contact()
        row["asserted"]["family"] = "other"
        row["residueA"] = {"chain": "LEGACY", "seq": "999"}
        hidden = copy.deepcopy(row)
        hidden["debugOnly"] = True
        hidden["semantics"]["participants"][0]["site"]["residue"]["resName"] = "HIDDEN"
        report = {"contacts": {"other": [row, hidden]}}
        text = explain_report(report)
        self.assertIn("1 hydrogen-bond assignment", text)
        self.assertIn("Canonical evidence levels: 1 candidate", text)
        self.assertIn("A:423 (1 assignment)", text)
        self.assertNotIn("HIDDEN", text)
        self.assertNotIn("LEGACY:999", text)
        self.assertEqual(contact_counts(report["contacts"]), {"hydrogen_bonds": 1})
        row["semantics"]["family"] = "polar_proximal"
        row["asserted"]["family"] = "hbond"
        self.assertEqual(contact_counts({"hydrogen_bonds": [row]}), {"polar_proximal": 1})
        self.assertTrue(describe_contact(row, "hydrogen_bonds").startswith("polar proximal"))

    def test_base_pair_metadata_is_context_not_duplicate_family_count(self):
        row = contact()
        row["basePair"] = {"annotation": "watson_crick_candidate"}
        text = explain_report({"contacts": {"hydrogen_bonds": [row]}})
        self.assertIn("1 hydrogen-bond assignment", text)
        self.assertIn("Base-pair context accompanies 1 assignment(s)", text)
        self.assertNotIn("1 base-pairing contact", text)

    def test_base_pair_context_does_not_promote_conflicting_polar_roles_to_hbond(self):
        row = contact("polar_contact", roles=("unresolved", "contact_atom"), directionality="symmetric", distance_kind="atom_pair", level="ambiguous")
        row["basePair"] = {"annotation": "watson_crick_candidate", "family": "G-C"}
        row["asserted"]["family"] = "base_pairing"
        participant_a = row["semantics"]["participants"][0]
        participant_a["site"]["residue"]["resName"] = "DG"
        participant_a["site"]["atoms"][0]["atomName"] = "N1"
        participant_a["roleAlternatives"] = [
            {"role": "acceptor", "source": "arpeggio_openbabel_atom_types"},
            {"role": "donor", "source": "standard_residue_template"},
        ]
        row["semantics"]["evidence"]["ambiguityFlags"] = ["nucleobase_typing_template_conflict"]
        row["semantics"]["direction"]["certainty"] = "ambiguous"
        text = explain_report({"contacts": {"base_pairing": [row]}})
        self.assertEqual(contact_counts({"base_pairing": [row]}), {"polar_contacts": 1})
        self.assertIn("1 polar contact", text)
        self.assertIn("Base-pair context accompanies 1 assignment(s)", text)
        self.assertIn("role alternatives: acceptor from arpeggio openbabel atom types; donor from standard residue template", text)
        self.assertIn("nucleobase typing template conflict", text)
        self.assertNotIn("→", text)
        self.assertNotIn("1 hydrogen-bond assignment", text)
        self.assertNotIn("1 base-pairing contact", text)

    def test_examples_are_bounded_without_serializing_every_contact(self):
        rows = [contact() for _ in range(250)]
        with patch("api.explain.describe_contact", wraps=describe_contact) as describe:
            text = explain_report({"contacts": {"hydrogen_bonds": rows}})
        self.assertEqual(describe.call_count, 1)
        self.assertIn("250 hydrogen-bond assignments", text)
        self.assertEqual(text.count("SER A:423 OG donor →"), 1)

    def test_api_json_roundtrip_preserves_explanation_and_does_not_mutate_record(self):
        row = contact()
        before = copy.deepcopy(row)
        self.assertEqual(describe_contact(row), describe_contact(json.loads(json.dumps(row))))
        self.assertEqual(row, before)

    def test_canonical_and_compact_api_delivery_preserve_surviving_roles(self):
        row = contact()
        hidden = copy.deepcopy(row)
        hidden["debugOnly"] = True
        delivery = prepare_delivery({"contacts": {"hydrogen_bonds": [row, hidden]}})
        canonical, compact = json.loads(delivery.canonical), json.loads(delivery.display)
        self.assertEqual(len(canonical["contacts"]["hydrogen_bonds"]), 2)
        self.assertEqual(len(compact["contacts"]["hydrogen_bonds"]), 1)
        for report in (canonical, compact):
            retained = report["contacts"]["hydrogen_bonds"][0]
            self.assertEqual(retained["semantics"], row["semantics"])
            self.assertEqual(describe_contact(retained), describe_contact(row))
        self.assertEqual(contact_counts(canonical["contacts"]), contact_counts(compact["contacts"]))


if __name__ == "__main__":
    unittest.main()
