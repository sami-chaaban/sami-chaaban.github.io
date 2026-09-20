from __future__ import annotations

import math
from typing import List, Optional

from .analysis_worker import is_diagnostic_contact


CONTACT_LABELS = {
    "hydrophobic": ("hydrophobic contact", "hydrophobic contacts"),
    "hydrogen_bonds": ("hydrogen-bond assignment", "hydrogen-bond assignments"),
    "polar_contacts": ("polar contact", "polar contacts"),
    "packing_contact": ("packing contact", "packing contacts"),
    "salt_bridges": ("salt bridge", "salt bridges"),
    "halogen_bonds": ("halogen-bond assignment", "halogen-bond assignments"),
    "metal_coordination": ("metal-coordination contact", "metal-coordination contacts"),
    "pi_pi": ("π–π contact", "π–π contacts"),
    "pi_cation": ("cation–π contact", "cation–π contacts"),
    "aromatic_packing": ("aromatic-packing contact", "aromatic-packing contacts"),
    "base_pairing": ("base-pairing contact", "base-pairing contacts"),
    "clash": ("steric clash", "steric clashes"),
    "other": ("other contact", "other contacts"),
}
FAMILY_BUCKETS = {
    "hbond": "hydrogen_bonds", "hydrogen_bond": "hydrogen_bonds",
    "polar_contact": "polar_contacts", "salt_bridge": "salt_bridges",
    "halogen_bond": "halogen_bonds",
}
EVIDENCE_LABELS = {
    "candidate": "candidate",
    "chemically_supported": "chemically supported",
    "geometrically_supported": "geometrically supported",
    "ambiguous": "ambiguous",
}
DISTANCE_LABELS = {
    "donor_acceptor": "donor–acceptor heavy-atom distance",
    "metal_donor": "metal–donor distance",
    "ring_centroid": "ring centroid distance",
    "closest_atom": "closest atom distance",
    "charge_site": "charged-site separation",
    "cation_centroid": "cation–ring centroid distance",
    "atom_ring_centroid": "atom–ring centroid distance",
    "atom_pair": "atom-pair distance",
}
MEASUREMENT_LABELS = {
    **DISTANCE_LABELS,
    "hbond_angle": "donor–H–acceptor angle",
    "donor_h_acceptor_angle": "donor–H–acceptor angle",
    "donor_hydrogen_acceptor_angle": "donor–H–acceptor angle",
    "halogen_angle": "donor-anchor–halogen···acceptor angle",
    "donor_anchor_halogen_acceptor_angle": "donor-anchor–halogen···acceptor angle",
    "ring_normal_angle": "ring-normal angle",
    "normal_angle": "ring-normal angle",
    "ring_lateral_offset": "ring lateral offset",
    "lateral_offset": "lateral offset",
    "ring_interplanar_distance": "ring interplanar distance",
    "interplanar_distance": "ring interplanar distance",
    "vdw_overlap": "van der Waals overlap",
    "vdw_reference": "van der Waals reference separation",
    "vdw_sum": "sum of van der Waals radii",
    "vdw_radius_sum": "sum of van der Waals radii",
    "cation_plane_offset": "cation-to-ring-plane distance",
    "cation_lateral_offset": "cation lateral offset",
    "cation_ring_normal_angle": "cation–centroid vector angle from ring normal",
    "ring_plane_separation_A": "separation from ring A plane",
    "ring_plane_separation_B": "separation from ring B plane",
    "ring_lateral_offset_A": "lateral offset in ring A plane",
    "ring_lateral_offset_B": "lateral offset in ring B plane",
    "hydrogen_acceptor": "H–acceptor distance",
}
ROLE_LABELS = {
    "donor": "donor", "acceptor": "acceptor",
    "donor_or_acceptor": "donor/acceptor unresolved",
    "halogen_donor": "halogen donor",
    "metal_center": "metal center", "coordinating_atom": "coordinating atom",
    "positive_site": "positive site", "negative_site": "negative site",
    "cation": "cationic site", "aromatic_ring": "aromatic ring",
    "contact_atom": "contacting atom", "base": "base",
    "unresolved": "role unresolved",
}
GROUP_LABELS = {
    "hydrophobic": ("nonpolar contact region", "nonpolar contact regions"),
    "hydrogen_bonds": ("hydrogen-bond assignment", "hydrogen-bond assignments"),
    "salt_bridges": ("charged-feature pair", "charged-feature pairs"),
    "halogen_bonds": ("halogen-bond assignment", "halogen-bond assignments"),
    "metal_coordination": ("metal–donor edge", "metal–donor edges"),
    "pi_pi": ("ring-pair assignment", "ring-pair assignments"),
    "pi_cation": ("cation–ring assignment", "cation–ring assignments"),
    "aromatic_packing": ("aromatic contact region", "aromatic contact regions"),
    "polar_contacts": ("polar assignment", "polar assignments"),
    "packing_contact": ("packing assignment", "packing assignments"),
    "clash": ("atom-pair clash", "atom-pair clashes"),
}


def canonical_semantics(contact: dict) -> Optional[dict]:
    """Consumers never recover chemical roles from endpoint names or buckets."""
    value = contact.get("semantics")
    return value if isinstance(value, dict) and value.get("version") == 1 else None


def _human(value) -> str:
    return str(value or "unknown").replace("_", " ")


def _number(value) -> Optional[str]:
    if isinstance(value, bool):
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(parsed):
        return None
    return f"{parsed:.3f}".rstrip("0").rstrip(".") if parsed else "0"


def _participant_label(participant: dict) -> str:
    site = participant.get("site") or {}
    residue = site.get("residue") or {}
    chain, seq = residue.get("chain"), residue.get("seq")
    residue_id = f"{chain}:{seq}" if chain is not None and seq is not None else "unidentified residue"
    residue_name = str(residue.get("resName") or "").strip()
    prefix = " ".join(value for value in (residue_name, residue_id) if value)
    atoms = site.get("atoms") or []
    atom_names = [str(atom.get("atomName")) for atom in atoms if isinstance(atom, dict) and atom.get("atomName")]
    role = ROLE_LABELS.get(participant.get("role"), _human(participant.get("role")))
    kind = site.get("kind")
    if kind in {"ring", "group"}:
        feature = f"[{', '.join(atom_names)}]" if atom_names else str(site.get("id") or f"unidentified {kind}")
        role = role if kind in role else f"{kind}, {role}"
    else:
        feature = ", ".join(atom_names) or str(site.get("id") or "unidentified atom")
    result = f"{prefix} {feature} {role}"
    alternatives = [
        f"{ROLE_LABELS.get(item.get('role'), _human(item.get('role')))} from {_human(item.get('source'))}"
        for item in participant.get("roleAlternatives") or []
        if isinstance(item, dict) and item.get("role")
    ]
    if alternatives:
        result += f" (role alternatives: {'; '.join(alternatives)})"
    contact_atom = site.get("contactAtom")
    if kind == "group" and isinstance(contact_atom, dict) and contact_atom.get("atomName"):
        result += f" (distance endpoint: {contact_atom['atomName']})"
    charge = participant.get("charge")
    if isinstance(charge, dict):
        value = _number(charge.get("value"))
        if value is not None:
            sign = "+" if float(value) > 0 else ""
            status = "inferred" if charge.get("inferred") else "reported"
            scope = "group charge" if charge.get("scope") == "group" else "charge"
            result += f" ({scope} {sign}{value}, {status}; {_human(charge.get('source'))})"
            atom_value = _number(charge.get("atomValue"))
            if charge.get("scope") == "group" and atom_value is not None:
                atom_sign = "+" if float(atom_value) > 0 else ""
                result += f" (distance-endpoint atomic charge {atom_sign}{atom_value})"
        elif charge.get("sign") in {-1, 1}:
            sign = "positive" if charge["sign"] == 1 else "negative"
            status = "inferred" if charge.get("inferred") else "reported"
            result += f" ({sign} charge sign {status}, magnitude unspecified; {_human(charge.get('source'))})"
    return result


def _geometry_text(geometry: dict) -> str:
    descriptions = []
    distance = geometry.get("distance")
    if isinstance(distance, dict):
        value = _number(distance.get("value"))
        if value is not None and float(value) >= 0:
            kind = distance.get("kind")
            label = DISTANCE_LABELS.get(kind, f"reported {_human(kind)} metric")
            unit = "Å" if distance.get("unit") == "angstrom" else _human(distance.get("unit"))
            descriptions.append(f"{label} {value} {unit} (source: {_human(distance.get('source'))})")
    for item in geometry.get("measurements") or []:
        if not isinstance(item, dict):
            continue
        value = _number(item.get("value"))
        if value is None:
            continue
        label = MEASUREMENT_LABELS.get(item.get("kind"), _human(item.get("kind")))
        unit = {"angstrom": "Å", "degree": "°", "degrees": "°"}.get(item.get("unit"), _human(item.get("unit")))
        descriptions.append(f"{label} {value} {unit} (source: {_human(item.get('source'))})")
    for field, label in (("hydrogen", "hydrogen"), ("donorAnchor", "donor anchor")):
        atom = geometry.get(field)
        if not isinstance(atom, dict):
            continue
        identity = atom.get("atomName") or atom.get("id") or "identity unavailable"
        source = _human(atom.get("source"))
        bond_source = f"; donor association: {_human(atom['bondSource'])}" if atom.get("bondSource") else ""
        descriptions.append(f"{label} {identity} (source: {source}{bond_source})")
    return "; ".join(descriptions) if descriptions else "measurement metadata unavailable"


def describe_contact(contact: dict, bucket: str = "other") -> str:
    """Describe one canonical contact without inventing roles, state or geometry.

    Suitable for a future focused explanation as well as bounded report examples.
    Legacy labels remain available but never become chemical-role evidence.
    """
    family = contact_family(bucket, contact)
    family_label = CONTACT_LABELS.get(family, (family.replace("_", " "),))[0]
    semantics = canonical_semantics(contact)
    if semantics is None:
        return f"{family_label}: canonical role and measurement metadata unavailable; no direction or chemical state inferred from the legacy label."
    evidence = semantics.get("evidence") or {}
    level = EVIDENCE_LABELS.get(evidence.get("level"), "support unspecified")
    participants = {
        p.get("side"): p for p in semantics.get("participants") or []
        if isinstance(p, dict) and p.get("side") in {"A", "B"}
    }
    direction = semantics.get("direction") or {}
    from_side, to_side = direction.get("from"), direction.get("to")
    directed = (
        semantics.get("directionality") == "directional"
        and direction.get("certainty") in {"certain", "inferred"}
        and from_side in participants and to_side in participants and from_side != to_side
        and participants[from_side].get("role") in {"donor", "halogen_donor"}
        and participants[to_side].get("role") == "acceptor"
    )
    ordered = [participants[from_side], participants[to_side]] if directed else [participants[side] for side in ("A", "B") if side in participants]
    if directed:
        separator = " → "
    elif semantics.get("directionality") == "directional":
        separator = " — "
    elif any(p.get("role") == "metal_center" for p in ordered):
        separator = " — "
    else:
        separator = " ↔ "
    partners = separator.join(_participant_label(p) for p in ordered) or "participant metadata unavailable"
    if semantics.get("directionality") == "directional":
        partners += " (inferred direction)" if directed and direction.get("certainty") == "inferred" else ""
        if not directed:
            partners += " (direction unresolved)"
    provenance = sorted({
        _human(source) for participant in ordered for source in participant.get("roleProvenance") or []
        if isinstance(source, str) and source
    })
    evidence_text = (
        f"Chemical compatibility: {_human(evidence.get('chemicalCompatibility'))}; "
        f"chemical state: {_human(evidence.get('chemicalState'))}; "
        f"geometry support: {_human(evidence.get('geometrySupport'))}"
    )
    if provenance:
        evidence_text += f"; role sources: {', '.join(provenance)}"
    site_sources = sorted({_human(p['site']['provenance']) for p in ordered if isinstance(p.get('site'), dict) and p['site'].get('provenance')})
    if site_sources:
        evidence_text += f"; site sources: {', '.join(site_sources)}"
    ambiguity = [_human(flag) for flag in evidence.get("ambiguityFlags") or [] if isinstance(flag, str)]
    if ambiguity:
        evidence_text += f"; ambiguity: {', '.join(ambiguity)}"
    return f"{family_label} ({level}): {partners}. Model geometry: {_geometry_text(semantics.get('geometry') or {})}. {evidence_text}. Biological and energetic interpretation not evaluated."


def non_diagnostic_contacts(contacts: dict):
    """Use the same records for the narrative counts and residue ranking."""
    for bucket, records in contacts.items():
        if not isinstance(records, list):
            continue
        for record in records:
            if isinstance(record, dict) and not is_diagnostic_contact(record):
                yield bucket, record


def contact_family(bucket: str, contact: dict) -> str:
    semantics = canonical_semantics(contact)
    if semantics:
        semantic_family = FAMILY_BUCKETS.get(semantics.get("family"), semantics.get("family"))
        if isinstance(semantic_family, str) and semantic_family.strip():
            return semantic_family.strip()
    asserted = contact.get("asserted")
    family = str(asserted.get("family") or "").strip().lower() if isinstance(asserted, dict) else ""
    family = FAMILY_BUCKETS.get(family, family)
    # Canonical reports can store an asserted clash in the broad 'other' bucket.
    return family if family in CONTACT_LABELS else FAMILY_BUCKETS.get(bucket, bucket)


def _count(value) -> Optional[int]:
    return value if isinstance(value, int) and not isinstance(value, bool) and value >= 0 else None


def display_group_summary(report: dict) -> Optional[dict]:
    """Count declared grouping units without confusing source observations with atoms.

    Grouping is produced by the chemistry pipeline; this consumer never clusters
    contacts or guesses atom support for a ring/feature observation.
    """
    payload = report.get("displayGroups")
    if not isinstance(payload, dict) or payload.get("version") != 1 or not isinstance(payload.get("groups"), list):
        return None
    records = {}
    retained_count = 0
    for bucket, record in non_diagnostic_contacts(report.get("contacts") or {}):
        retained_count += 1
        semantics = canonical_semantics(record)
        identity = semantics.get("identity") if semantics else None
        if identity:
            records[identity] = (bucket, record)
    families, groups, covered, group_ids, support_ids = {}, [], set(), set(), set()
    for group in payload["groups"]:
        if not isinstance(group, dict) or group.get("debugOnly") or group.get("debug_only"):
            continue
        identity = group.get("id")
        members = group.get("contactIds")
        counts = group.get("counts") or {}
        units, objects = _count(counts.get("chemicalUnits")), _count(counts.get("displayObjects"))
        if not identity or identity in group_ids or not isinstance(members, list) or not members or units is None or objects is None:
            continue
        # A stale/mixed scope must not contribute hidden contacts or cached totals.
        if any(member not in records for member in members):
            continue
        family = FAMILY_BUCKETS.get(group.get("family"), group.get("family"))
        if not isinstance(family, str) or not family or any(contact_family(*records[member]) != family for member in members):
            continue
        # Primary groups partition canonical records; overlapping context views
        # must not create another counted chemical unit.
        if covered.intersection(members):
            continue
        group_ids.add(identity)
        covered.update(members)
        groups.append(group)
        values = families.setdefault(family, {"chemicalUnits": 0, "displayObjects": 0, "canonicalContacts": 0,
                                             "supportIds": set(), "unknownSupportGroups": 0, "example": group})
        values["chemicalUnits"] += units
        values["displayObjects"] += objects
        values["canonicalContacts"] += len(set(members))
        atom_contacts = group.get("supportingAtomContacts")
        declared_support = _count(counts.get("supportingAtomContacts"))
        ids = {item.get("id") for item in atom_contacts if isinstance(item, dict) and isinstance(item.get("id"), str) and item["id"]} if isinstance(atom_contacts, list) else set()
        values["supportIds"].update(ids)
        support_ids.update(ids)
        if declared_support is None or declared_support != len(ids):
            values["unknownSupportGroups"] += 1
    return {"families": families, "groups": groups, "records": records, "coveredContacts": len(covered),
            "canonicalContacts": retained_count, "chemicalUnits": sum(v["chemicalUnits"] for v in families.values()),
            "displayObjects": sum(v["displayObjects"] for v in families.values()),
            "knownSupportingAtomContacts": len(support_ids),
            "unknownSupportGroups": sum(v["unknownSupportGroups"] for v in families.values())}


def _group_example(group: dict, records: dict) -> str:
    family = FAMILY_BUCKETS.get(group.get("family"), group.get("family"))
    label = GROUP_LABELS.get(family, (f"{_human(family)} unit",))[0]
    count = _count((group.get("counts") or {}).get("supportingAtomContacts"))
    ids = {item.get("id") for item in group.get("supportingAtomContacts") or [] if isinstance(item, dict) and isinstance(item.get("id"), str) and item["id"]}
    if count != len(ids):
        count = None
    support = f"{count} supporting atom contact(s)" if count is not None else "supporting atom-contact count unavailable"
    text = f"{label}: {support}; grouping rule {_human(group.get('rule'))}"
    geometry = group.get("geometry") or {}
    distances = []
    for field, name in (("closestAtomDistance", "closest supporting atom distance"), ("meanAtomDistance", "mean supporting atom distance"), ("medianAtomDistance", "median supporting atom distance")):
        value = _number(geometry.get(field))
        if value is not None:
            distances.append(f"{name} {value} Å")
        elif field == "closestAtomDistance":
            distances.append("closest supporting atom distance unavailable")
    centroid = _number(geometry.get("centroidDistance"))
    if centroid is not None:
        label = "schematic feature-center separation" if (group.get("renderGeometry") or {}).get("schematic") else "feature-centroid separation"
        distances.append(f"{label} {centroid} Å")
    if distances:
        text += "; " + "; ".join(distances)
    flags = [_human(flag) for flag in group.get("ambiguityFlags") or [] if isinstance(flag, str)]
    if flags:
        text += "; grouping ambiguity: " + ", ".join(flags)
    representative = group.get("representativeContactId")
    if representative not in group.get("contactIds", []):
        representative = None
    if representative in records:
        bucket, record = records[representative]
        text += ". Representative retained assignment: " + describe_contact(record, bucket)
    return text + ("" if text.endswith(".") else ".")


def _top_group_residues(grouping: dict, limit: int = 3) -> List[dict]:
    ranked = {}
    for group in grouping["groups"]:
        seen = set()
        for identity in group["contactIds"]:
            _, record = grouping["records"][identity]
            semantics = canonical_semantics(record) or {}
            for participant in semantics.get("participants") or []:
                residue = (participant.get("site") or {}).get("residue") or {}
                chain, seq = residue.get("chain"), residue.get("seq")
                if chain is None or seq is None or str(chain) == "" or str(seq) == "":
                    continue
                key = f"{chain}:{seq}"
                if key in seen:
                    continue
                seen.add(key)
                entry = ranked.setdefault(key, {"id": key, "total": 0})
                entry["total"] += group["counts"]["chemicalUnits"]
    return sorted(ranked.values(), key=lambda item: (-item["total"], item["id"]))[:limit]


def explain_report(report: dict, images: Optional[List[str]] = None, notes: Optional[str] = None) -> str:
    images = images or []
    contacts = report.get("contacts") or {}
    counts = contact_counts(contacts)
    grouping = display_group_summary(report)
    chain_a = report.get("chainA", "A")
    chain_b = report.get("chainB", "B")

    def counted_label(bucket, count):
        labels = CONTACT_LABELS.get(bucket, (bucket.replace('_', ' ') + ' contact', bucket.replace('_', ' ') + ' contacts'))
        return f"{count} {labels[0 if count == 1 else 1]}"

    contact_text = ", ".join(counted_label(bucket, count) for bucket, count in counts.items() if count)
    scope = f"within chain {chain_a}" if chain_a == chain_b else f"between chains {chain_a} and {chain_b}"
    summary = f"The analysis {scope} " + (
        f"includes {contact_text}." if contact_text else "contains no non-diagnostic contact assignments under the current criteria."
    )
    grouping_text = ""
    if grouping is not None and (grouping["groups"] or grouping["canonicalContacts"] == 0):
        family_parts = []
        for family, values in grouping["families"].items():
            units = values["chemicalUnits"]
            labels = GROUP_LABELS.get(family, (f"{_human(family)} unit", f"{_human(family)} units"))
            support_count = len(values["supportIds"])
            atom_support = f"{support_count} known supporting atom {'contact' if support_count == 1 else 'contacts'}"
            if values["unknownSupportGroups"]:
                unknown = f"atom-support count unavailable for {values['unknownSupportGroups']} group(s)"
                atom_support = atom_support + "; " + unknown if support_count else unknown
            family_parts.append(f"{units} {labels[0 if units == 1 else 1]} ({atom_support})")
        summary = (
            f"The analysis {scope} has {grouping['chemicalUnits']} assigned interaction {'unit' if grouping['chemicalUnits'] == 1 else 'units'} "
            f"represented by {grouping['displayObjects']} available display {'object' if grouping['displayObjects'] == 1 else 'objects'}"
            + (": " + ", ".join(family_parts) if family_parts else "") + "."
        )
        grouping_text = (
            f"The grouping covers {grouping['coveredContacts']} of {grouping['canonicalContacts']} retained canonical contact records "
            f"and contains {grouping['knownSupportingAtomContacts']} known supporting atom contacts. "
            "Chemical units use family-specific features or atom assignments; display objects are available groups, "
            "not the number of connectors actually drawn in a particular view. "
        )
        if grouping["unknownSupportGroups"]:
            grouping_text += f"Atom-contact support is not fully exported for {grouping['unknownSupportGroups']} group(s); missing counts are not zero. "
        if grouping["coveredContacts"] < grouping["canonicalContacts"]:
            grouping_text += "Ungrouped records remain in the canonical counts; their chemical-unit and display counts are unavailable. "
        if contact_text:
            grouping_text += "Underlying canonical records: " + contact_text + ". "
        grouping_text += "Multiple atom contacts can describe one nonpolar region; neither their number nor the number of groups measures independent bonds or binding energy."
    elif counts:
        grouping_text = "These counts are retained canonical contact records; family-specific chemical-unit and display-object counts are unavailable in this report."
    interpretation = (
        "These assignments describe geometric and chemical compatibility in the supplied model; "
        "evidence levels summarize the available support, not calibrated probabilities or binding energies. "
        "Hydrogen-bond assignments depend on donor/acceptor typing, protonation and hydrogen placement; "
        "missing hydrogens or uncertain ligand chemistry limit certainty. "
        "Contact counts and short distances alone do not establish net stabilization or binding affinity."
    )

    evidence_counts = {}
    examples = {}
    missing_semantics = 0
    base_pair_context_count = 0
    for bucket, contact in non_diagnostic_contacts(contacts):
        semantics = canonical_semantics(contact)
        if semantics:
            evidence = semantics.get("evidence") or {}
            level = EVIDENCE_LABELS.get(evidence.get("level"), "support unspecified")
            evidence_counts[level] = evidence_counts.get(level, 0) + 1
            family = contact_family(bucket, contact)
            if family not in examples and not (grouping and grouping["groups"]):
                examples[family] = describe_contact(contact, bucket)
        else:
            missing_semantics += 1
        if isinstance(contact.get("basePair") or (contact.get("asserted") or {}).get("basePair"), dict):
            base_pair_context_count += 1
    semantic_parts = []
    if evidence_counts:
        semantic_parts.append("Canonical evidence levels: " + ", ".join(f"{number} {label}" for label, number in evidence_counts.items()) + " (counts of canonical records).")
    if missing_semantics:
        semantic_parts.append(f"Canonical role metadata is unavailable for {missing_semantics} assignment(s); no donor, charge or ring roles have been reconstructed from their names.")
    if base_pair_context_count:
        semantic_parts.append(f"Base-pair context accompanies {base_pair_context_count} assignment(s); it is separate from their atom-level families and is not a count of unique base pairs.")
    if grouping and grouping["groups"]:
        semantic_parts.append("Representative group examples (one per grouped family, selected for display rather than energetic importance): " + " ".join(
            _group_example(values["example"], grouping["records"]) for values in grouping["families"].values()
        ))
    elif examples:
        semantic_parts.append("Representative role and geometry examples (one per retained family, not an energetic ranking): " + " ".join(examples.values()))
    semantics_text = " ".join(semantic_parts)

    grouped_ranking = grouping is not None and grouping["groups"] and grouping["coveredContacts"] == grouping["canonicalContacts"]
    hotspots = _top_group_residues(grouping) if grouped_ranking else top_residues(report.get("perResidue", {}), limit=3, contacts=contacts)
    hotspot_text = ""
    if hotspots:
        unit = "assigned unit" if grouped_ranking else "assignment"
        hotspot_list = ", ".join(f"{item['id']} ({item['total']} {unit if item['total'] == 1 else unit + 's'})" for item in hotspots)
        basis = "assigned chemical units" if grouped_ranking else "retained canonical assignments"
        hotspot_text = f"The most connected residues by {basis} are {hotspot_list}; these counts do not identify energetic hotspots."

    image_text = ""
    if images:
        image_text = f"{len(images)} image reference(s) were supplied; this text-only explanation has not inspected them."
    note_text = f"User note: {notes.strip()}" if notes else ""
    return " ".join(part for part in [summary, grouping_text, interpretation, semantics_text, hotspot_text, image_text, note_text] if part)


def contact_counts(contacts: dict) -> dict:
    counts = {}
    for bucket, contact in non_diagnostic_contacts(contacts):
        family = contact_family(bucket, contact)
        counts[family] = counts.get(family, 0) + 1
    return counts


def top_residues(per_residue: dict, limit: int = 3, *, contacts: Optional[dict] = None) -> List[dict]:
    if contacts is not None:
        ranked = {}
        for _, contact in non_diagnostic_contacts(contacts):
            seen = set()
            semantics = canonical_semantics(contact)
            residues = [
                (participant.get("site") or {}).get("residue")
                for participant in semantics.get("participants") or []
                if isinstance(participant, dict)
            ] if semantics else [contact.get(side) for side in ("residueA", "residueB")]
            for residue in residues:
                if not isinstance(residue, dict):
                    continue
                chain, seq = residue.get("chain"), residue.get("seq")
                if chain is None or seq is None or str(chain) == "" or str(seq) == "":
                    continue
                key = f"{chain}:{seq}"
                if key in seen:
                    continue
                seen.add(key)
                entry = ranked.setdefault(key, {"id": key, "total": 0, "resName": residue.get("resName", "")})
                entry["total"] += 1
        entries = ranked.values()
    else:
        entries = (
            {"id": key, "total": value.get("total", 0), "resName": value.get("resName", "")}
            for key, value in per_residue.items() if isinstance(value, dict)
        )
    return sorted(entries, key=lambda item: (-item["total"], item["id"]))[:limit]
