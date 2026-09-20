"""Explanation counting units must follow producer grouping, not raw density."""
import copy
import unittest
from unittest.mock import patch

from api.explain import describe_contact, display_group_summary, explain_report
from api.test_explain_semantics import contact


def row(identity, family="hydrophobic", *, debug=False):
    result = contact(family, roles=("contact_atom", "contact_atom"), directionality="symmetric", distance_kind="atom_pair")
    result["semantics"]["identity"] = identity
    if family == "hydrophobic":
        suffix = "".join(value for value in identity if value.isdigit())
        index = int(suffix or "0")
        for side, residue_name, atom_name in (("A", "LEU", ["CG", "CD1", "CD2"][(index // 3) % 3]), ("B", "VAL", ["CB", "CG1", "CG2"][index % 3])):
            participant = next(p for p in result["semantics"]["participants"] if p["side"] == side)
            participant["site"]["residue"]["resName"] = residue_name
            participant["site"]["atoms"][0].update(id=f"{side}:{residue_name}:{atom_name}", atomName=atom_name, element="C")
    if debug:
        result["debugOnly"] = True
    return result


def group(identity, rows, *, family=None, known=True, support_ids=None):
    ids = [value["semantics"]["identity"] for value in rows]
    supports = [{"id": value, "canonicalContactId": ids[index % len(ids)], "distance": 3.2} for index, value in enumerate(support_ids or ids)] if known else []
    return {"id": identity, "family": family or rows[0]["semantics"]["family"], "rule": "nonpolar_feature_local_patch",
            "featureIds": ["featureA", "featureB"], "contactIds": ids, "representativeContactId": ids[0],
            "counts": {"canonicalContacts": len(ids), "chemicalUnits": 1, "displayObjects": 1, "supportingAtomContacts": len(supports) if known else None},
            "supportingAtomContacts": supports, "geometry": {"closestAtomDistance": 3.2, "meanAtomDistance": 3.6, "medianAtomDistance": 3.5, "centroidDistance": 3.8},
            "renderGeometry": {"kind": "hydrophobic_patch", "schematic": True}, "ambiguityFlags": [] if known else ["atom_support_not_exported"]}


def report(rows, groups):
    return {"contacts": {"other": rows}, "displayGroups": {"version": 1, "groups": groups, "features": {}}}


class GroupedExplanationTests(unittest.TestCase):
    def test_dense_hydrophobic_support_is_one_unit_and_one_available_object(self):
        rows = [row(f"pair{i}") for i in range(3)]
        rows[0]["semantics"]["evidence"]["supportingRecords"] = 90
        data = report(rows, [group("patch", rows)])
        text = explain_report(data)
        self.assertIn("1 assigned interaction unit represented by 1 available display object", text)
        self.assertIn("1 nonpolar contact region (3 known supporting atom contacts)", text)
        self.assertIn("covers 3 of 3 retained canonical contact records", text)
        self.assertIn("A:423 (1 assigned unit)", text)
        self.assertNotIn("A:423 (3", text)
        self.assertNotIn("90", text)
        self.assertIn("closest supporting atom distance 3.2 Å", text)
        self.assertIn("schematic feature-center separation 3.8 Å", text)
        self.assertIn("not the number of connectors actually drawn", text)

    def test_distinct_nonpolar_features_on_same_residues_remain_two_units(self):
        rows = [row("tail1"), row("tail2")]
        data = report(rows, [group("patch1", rows[:1]), group("patch2", rows[1:])])
        self.assertEqual(display_group_summary(data)["chemicalUnits"], 2)
        self.assertIn("2 nonpolar contact regions", explain_report(data))

    def test_support_identity_is_deduplicated_and_diagnostics_do_not_leak(self):
        rows = [row("pair1"), row("pair2"), row("hidden", debug=True)]
        groups = [group("p1", rows[:1], support_ids=["same-reversed-pair"]), group("p2", rows[1:2], support_ids=["same-reversed-pair"]), group("hidden", rows[2:])]
        groups[-1]["debugOnly"] = True
        result = display_group_summary(report(rows, groups))
        self.assertEqual(result["chemicalUnits"], 2)
        self.assertEqual(result["knownSupportingAtomContacts"], 1)
        self.assertEqual(result["canonicalContacts"], 2)
        self.assertEqual(result["coveredContacts"], 2)

    def test_ring_feature_observation_does_not_invent_atom_pairs(self):
        rows = [row("ringpair", "pi_pi")]
        data = report(rows, [group("ringgroup", rows, known=False)])
        data["displayGroups"]["groups"][0]["geometry"] = {"closestAtomDistance": None}
        result = display_group_summary(data)
        self.assertEqual(result["unknownSupportGroups"], 1)
        text = explain_report(data)
        self.assertIn("1 ring-pair assignment", text)
        self.assertIn("supporting atom-contact count unavailable", text)
        self.assertIn("closest supporting atom distance unavailable", text)
        self.assertIn("missing counts are not zero", text)

    def test_two_metal_donors_are_two_edges_even_with_one_metal(self):
        rows = [row("metalO1", "metal_coordination"), row("metalO2", "metal_coordination")]
        text = explain_report(report(rows, [group("edge1", rows[:1]), group("edge2", rows[1:])]))
        self.assertIn("2 metal–donor edges", text)
        self.assertIn("A:423 (2 assigned units)", text)

    def test_native_atom_ring_contact_has_centroid_metric_without_fake_atom_support(self):
        rows = [row("native-atom-ring")]
        rows[0]["semantics"]["geometry"]["distance"].update(kind="atom_ring_centroid", value=4.06)
        ring = rows[0]["semantics"]["participants"][1]
        ring["role"] = "aromatic_ring"
        ring["site"]["kind"] = "ring"
        data = report(rows, [group("atom-ring", rows, known=False)])
        data["displayGroups"]["groups"][0]["geometry"] = {"closestAtomDistance": None}
        text = explain_report(data)
        self.assertIn("atom–ring centroid distance 4.06 Å", text)
        self.assertIn("closest supporting atom distance unavailable", text)
        self.assertNotIn("atom-pair distance 4.06", text)

    def test_missing_or_partial_grouping_is_not_a_complete_unit_count(self):
        rows = [row("covered"), row("unmapped")]
        data = report(rows, [group("patch", rows[:1])])
        text = explain_report(data)
        self.assertIn("covers 1 of 2 retained canonical contact records", text)
        self.assertIn("Ungrouped records remain", text)
        self.assertIn("A:423 (2 assignments)", text)
        del data["displayGroups"]
        self.assertIn("family-specific chemical-unit and display-object counts are unavailable", explain_report(data))

    def test_stale_scope_counts_and_overlapping_context_cannot_inflate_units(self):
        rows = [row("actual")]
        first = group("primary", rows)
        context = copy.deepcopy(first)
        context["id"] = "overlapping-context"
        stale = copy.deepcopy(first)
        stale["id"] = "old-scope"
        stale["contactIds"] = ["not-in-report"]
        data = report(rows, [first, context, stale])
        self.assertEqual(display_group_summary(data)["chemicalUnits"], 1)

    def test_base_pair_context_does_not_add_another_unit(self):
        rows = [row("polar", "polar_contact")]
        rows[0]["basePair"] = {"family": "G-C"}
        text = explain_report(report(rows, [group("polar-group", rows)]))
        self.assertIn("1 polar assignment", text)
        self.assertIn("Base-pair context accompanies 1 assignment(s)", text)
        self.assertNotIn("1 base-pairing contact", text)

    def test_examples_remain_bounded_and_input_records_are_not_mutated(self):
        rows = [row(f"pair{i}") for i in range(100)]
        data = report(rows, [group(f"patch{i}", [value]) for i, value in enumerate(rows)])
        original = copy.deepcopy(data)
        with patch("api.explain.describe_contact", wraps=describe_contact) as describe:
            text = explain_report(data)
        self.assertEqual(describe.call_count, 1)
        self.assertIn("100 nonpolar contact regions", text)
        self.assertEqual(data, original)


if __name__ == "__main__":
    unittest.main()
