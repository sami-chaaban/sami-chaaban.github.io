from __future__ import annotations

from typing import List, Optional


def explain_report(report: dict, images: Optional[List[str]] = None, notes: Optional[str] = None) -> str:
    images = images or []
    counts = contact_counts(report.get("contacts", {}))
    total = sum(counts.values())
    chain_a = report.get("chainA", "A")
    chain_b = report.get("chainB", "B")

    summary = (
        f"The interface between chain {chain_a} and chain {chain_b} includes "
        f"{counts.get('hydrophobic', 0)} hydrophobic contacts, "
        f"{counts.get('hydrogen_bonds', 0)} hydrogen bonds, "
        f"{counts.get('salt_bridges', 0)} salt bridges, "
        f"{counts.get('metal_coordination', 0)} metal coord. contacts, "
        f"{counts.get('pi_pi', 0)} pi-pi contacts, and "
        f"{counts.get('pi_cation', 0)} pi-cation contacts."
    )

    hotspots = top_residues(report.get("perResidue", {}), limit=3)
    hotspot_text = ""
    if hotspots:
        hotspot_list = ", ".join(
            f"{item['id']} ({item['total']} contacts)" for item in hotspots
        )
        hotspot_text = f"Key hotspots by contact density are {hotspot_list}."

    image_text = ""
    if images:
        image_text = (
            "In image 1 you can see the overall interface footprint; "
            "later images highlight the strongest contact clusters."
        )

    note_text = ""
    if notes:
        note_text = f"User note: {notes.strip()}"

    parts = [summary, hotspot_text, image_text, note_text]
    return " ".join(part for part in parts if part)


def contact_counts(contacts: dict) -> dict:
    return {key: len(value) for key, value in contacts.items()}


def top_residues(per_residue: dict, limit: int = 3) -> List[dict]:
    ranked = sorted(
        (
            {
                "id": key,
                "total": value.get("total", 0),
                "resName": value.get("resName", ""),
            }
            for key, value in per_residue.items()
        ),
        key=lambda item: item["total"],
        reverse=True,
    )
    return ranked[:limit]
