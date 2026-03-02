from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple
import hashlib
import math
import re
import shlex
import urllib.request

TOOL_VERSION = "toy-1.1"
MODEL_SERVER_URL = (
    "https://www.ebi.ac.uk/pdbe/model-server/v1/{pdb_id}/full"
    "?encoding=cif&data_source=pdb-h"
)

HYDROPHOBIC_RESIDUES = {
    "ALA",
    "VAL",
    "ILE",
    "LEU",
    "MET",
    "PHE",
    "TRP",
    "PRO",
    "TYR",
}
AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}
ACIDIC_RESIDUES = {"ASP", "GLU"}
BASIC_RESIDUES = {"LYS", "ARG", "HIS"}
STANDARD_AMINO_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "SEC", "PYL", "MSE",
}
NUCLEIC_RESIDUES = {
    "A", "C", "G", "U", "I",
    "DA", "DC", "DG", "DT", "DI", "DU",
    "ADE", "CYT", "GUA", "THY", "URA",
    "AMP", "CMP", "GMP", "TMP", "UMP",
}
POLYMER_RESIDUES = STANDARD_AMINO_RESIDUES | NUCLEIC_RESIDUES
WATER_RESIDUES = {"HOH", "WAT", "H2O"}
AROMATIC_CANDIDATE_ELEMENTS = {"C", "N", "O", "S"}
AROMATIC_BOND_MIN = 1.18
AROMATIC_BOND_MAX = 1.82
AROMATIC_PLANAR_RMSD_MAX = 0.34
AROMATIC_LOCAL_RING_PLANAR_RMSD_MAX = 0.46
METAL_ELEMENTS = {
    "LI", "NA", "K", "RB", "CS",
    "BE", "MG", "CA", "SR", "BA",
    "AL", "GA", "IN", "TL",
    "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
    "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD",
    "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG",
    "PB", "BI",
}
METAL_DONOR_ELEMENTS = {"O", "N", "S", "SE"}
METAL_COORDINATION_CUTOFF = {
    "LI": 2.35,
    "NA": 2.95,
    "K": 3.45,
    "RB": 3.6,
    "CS": 3.75,
    # Mg coordination can be undercalled in medium-resolution models with a strict 2.65 A gate.
    # Use a slightly wider bound so bidentate phosphate contacts are retained.
    "MG": 2.85,
    "CA": 3.05,
    "SR": 3.25,
    "BA": 3.45,
    "MN": 2.85,
    "FE": 2.85,
    "CO": 2.8,
    "NI": 2.75,
    "CU": 2.8,
    "ZN": 2.8,
    "CD": 2.95,
    "HG": 3.0,
    "AL": 2.45,
}


@dataclass(frozen=True)
class AtomRecord:
    chain_id: str
    chain_label: str
    res_name: str
    res_seq: str
    atom_name: str
    element: str
    x: float
    y: float
    z: float

    @property
    def residue_key(self) -> str:
        return f"{self.chain_id}:{self.res_seq}"


@dataclass(frozen=True)
class ChainAliases:
    label_to_auth: Dict[str, str]

    def normalize(self, chain_id: str) -> str:
        mapped = self.label_to_auth.get(chain_id, chain_id)
        # Keep uppercase/lowercase chain IDs distinct (e.g. "A" vs "a").
        if (
            isinstance(mapped, str)
            and isinstance(chain_id, str)
            and mapped != chain_id
            and mapped.lower() == chain_id.lower()
        ):
            return chain_id
        return mapped


@dataclass(frozen=True)
class ResidueProperties:
    aromatic: bool
    acidic: bool
    basic: bool
    hydrophobic: bool
    non_polymer: bool


def _distance_sq(atom_a: AtomRecord, atom_b: AtomRecord) -> float:
    dx = atom_a.x - atom_b.x
    dy = atom_a.y - atom_b.y
    dz = atom_a.z - atom_b.z
    return dx * dx + dy * dy + dz * dz


def _residue_element_counts(atoms: List[AtomRecord]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for atom in atoms:
        element = (atom.element or "").upper()
        if not element or element == "H":
            continue
        counts[element] = counts.get(element, 0) + 1
    return counts


def _best_plane_normal(points: List[Tuple[float, float, float]]) -> Optional[Tuple[float, float, float]]:
    if len(points) < 3:
        return None
    best = (0.0, 0.0, 0.0)
    best_norm_sq = 0.0
    for i in range(len(points) - 2):
        p0 = points[i]
        for j in range(i + 1, len(points) - 1):
            p1 = points[j]
            v1x = p1[0] - p0[0]
            v1y = p1[1] - p0[1]
            v1z = p1[2] - p0[2]
            for k in range(j + 1, len(points)):
                p2 = points[k]
                v2x = p2[0] - p0[0]
                v2y = p2[1] - p0[1]
                v2z = p2[2] - p0[2]
                cx = v1y * v2z - v1z * v2y
                cy = v1z * v2x - v1x * v2z
                cz = v1x * v2y - v1y * v2x
                norm_sq = cx * cx + cy * cy + cz * cz
                if norm_sq > best_norm_sq:
                    best = (cx, cy, cz)
                    best_norm_sq = norm_sq
    if best_norm_sq < 1e-10:
        return None
    inv_len = 1.0 / math.sqrt(best_norm_sq)
    return (best[0] * inv_len, best[1] * inv_len, best[2] * inv_len)


def _is_planar_component(atoms: List[AtomRecord], max_rmsd: float) -> bool:
    if len(atoms) < 3:
        return False
    points = [(atom.x, atom.y, atom.z) for atom in atoms]
    normal = _best_plane_normal(points)
    if not normal:
        return False
    cx = sum(p[0] for p in points) / len(points)
    cy = sum(p[1] for p in points) / len(points)
    cz = sum(p[2] for p in points) / len(points)
    rmsd_sq = 0.0
    for px, py, pz in points:
        d = abs((px - cx) * normal[0] + (py - cy) * normal[1] + (pz - cz) * normal[2])
        rmsd_sq += d * d
    rmsd = math.sqrt(rmsd_sq / len(points))
    return rmsd <= max_rmsd


def _is_likely_aromatic_nonpolymer(atoms: List[AtomRecord]) -> bool:
    candidates = [atom for atom in atoms if (atom.element or "").upper() in AROMATIC_CANDIDATE_ELEMENTS]
    if len(candidates) < 5 or len(candidates) > 80:
        return False

    min_d_sq = AROMATIC_BOND_MIN * AROMATIC_BOND_MIN
    max_d_sq = AROMATIC_BOND_MAX * AROMATIC_BOND_MAX
    adjacency: List[List[int]] = [[] for _ in candidates]
    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            d_sq = _distance_sq(candidates[i], candidates[j])
            if min_d_sq <= d_sq <= max_d_sq:
                adjacency[i].append(j)
                adjacency[j].append(i)

    def _is_aromatic_component(component_indices: List[int]) -> bool:
        size = len(component_indices)
        if size < 5 or size > 28:
            return False
        comp_set = set(component_indices)
        edge_count = 0
        for idx in component_indices:
            edge_count += sum(1 for nbr in adjacency[idx] if nbr in comp_set)
        edge_count //= 2
        if edge_count < size:
            return False
        avg_degree = (2.0 * edge_count) / size
        if avg_degree < 1.45:
            return False
        component_atoms = [candidates[idx] for idx in component_indices]
        return _is_planar_component(component_atoms, AROMATIC_PLANAR_RMSD_MAX)

    def _contains_planar_local_ring(component_indices: List[int]) -> bool:
        # Large ligands can contain aromatic sub-rings embedded in a much bigger
        # non-planar graph (e.g. taxol-like scaffolds). Detect local 5-7 atom
        # cycles instead of requiring the whole component to be aromatic-like.
        if len(component_indices) < 5:
            return False
        comp_set = set(component_indices)
        neighbors: Dict[int, List[int]] = {
            idx: [nbr for nbr in adjacency[idx] if nbr in comp_set]
            for idx in comp_set
        }
        min_cycle = 5
        max_cycle = 7

        for start in sorted(comp_set):
            stack: List[Tuple[int, List[int], set[int]]] = [(start, [start], {start})]
            while stack:
                current, path, visited_nodes = stack.pop()
                path_len = len(path)
                for nbr in neighbors.get(current, []):
                    if nbr == start:
                        if min_cycle <= path_len <= max_cycle:
                            ring_atoms = [candidates[idx] for idx in path]
                            if _is_planar_component(ring_atoms, AROMATIC_LOCAL_RING_PLANAR_RMSD_MAX):
                                return True
                        continue
                    if nbr in visited_nodes:
                        continue
                    # Cheap canonicalization to reduce equivalent traversals.
                    if nbr < start:
                        continue
                    if path_len >= max_cycle:
                        continue
                    next_path = path + [nbr]
                    next_visited = set(visited_nodes)
                    next_visited.add(nbr)
                    stack.append((nbr, next_path, next_visited))
        return False

    def _two_core(component_indices: List[int]) -> List[int]:
        # Strip leaf-like substituents so aromatic ring cores in large non-polymers
        # (e.g. nucleotides with sugar/phosphate tails) are still detected.
        active = set(component_indices)
        if len(active) < 3:
            return []
        degree: Dict[int, int] = {
            idx: sum(1 for nbr in adjacency[idx] if nbr in active)
            for idx in active
        }
        queue = [idx for idx, d in degree.items() if d < 2]
        while queue:
            idx = queue.pop()
            if idx not in active:
                continue
            active.remove(idx)
            for nbr in adjacency[idx]:
                if nbr not in active:
                    continue
                degree[nbr] = degree.get(nbr, 0) - 1
                if degree[nbr] == 1:
                    queue.append(nbr)
        return sorted(active)

    visited = [False] * len(candidates)
    for start in range(len(candidates)):
        if visited[start]:
            continue
        stack = [start]
        visited[start] = True
        component: List[int] = []
        while stack:
            idx = stack.pop()
            component.append(idx)
            for nbr in adjacency[idx]:
                if not visited[nbr]:
                    visited[nbr] = True
                    stack.append(nbr)

        # Try both the raw component and its 2-core (ring-like kernel).
        if _is_aromatic_component(component):
            return True
        core = _two_core(component)
        if core and core != component and _is_aromatic_component(core):
            return True
        if _contains_planar_local_ring(component):
            return True
        if core and core != component and _contains_planar_local_ring(core):
            return True
    return False


def _is_likely_acidic_nonpolymer(atoms: List[AtomRecord], counts: Dict[str, int]) -> bool:
    oxygen = counts.get("O", 0)
    carbon = counts.get("C", 0)
    nitrogen = counts.get("N", 0)
    phosphorus = counts.get("P", 0)
    sulfur = counts.get("S", 0)
    if oxygen >= 3 and (phosphorus >= 1 or sulfur >= 1):
        return True
    if oxygen < 2 or carbon < 1 or nitrogen > 0:
        return False
    oxygens = [atom for atom in atoms if (atom.element or "").upper() == "O"]
    carbons = [atom for atom in atoms if (atom.element or "").upper() == "C"]
    short_co_pairs = 0
    for oxygen_atom in oxygens:
        if any(_distance_sq(oxygen_atom, carbon_atom) <= 1.38 * 1.38 for carbon_atom in carbons):
            short_co_pairs += 1
    return short_co_pairs >= 2


def _is_likely_basic_nonpolymer(res_name: str, counts: Dict[str, int]) -> bool:
    if "+" in res_name:
        return True
    nitrogen = counts.get("N", 0)
    if nitrogen == 0:
        return False
    oxygen = counts.get("O", 0)
    sulfur = counts.get("S", 0)
    phosphorus = counts.get("P", 0)
    heavy_atoms = sum(counts.values())
    if heavy_atoms <= 2:
        return True
    if nitrogen >= 1 and oxygen == 0 and sulfur == 0 and phosphorus == 0 and heavy_atoms <= 24:
        return True
    if nitrogen >= 2 and oxygen <= 1 and phosphorus == 0 and heavy_atoms <= 28:
        return True
    return False


def _is_likely_hydrophobic_nonpolymer(counts: Dict[str, int], aromatic: bool) -> bool:
    heavy_atoms = sum(counts.values())
    carbon = counts.get("C", 0)
    hetero = heavy_atoms - carbon
    # Keep aromatic ligands (e.g. nucleotide bases) eligible for hydrophobic
    # contacts through their carbon-rich ring systems even when they also carry
    # strongly polar substituents.
    if aromatic and carbon >= 5 and heavy_atoms <= 80:
        return True
    if aromatic and carbon >= 5 and hetero <= 3:
        return True
    return carbon >= 6 and hetero <= 2


def build_residue_properties(atoms: List[AtomRecord]) -> Dict[str, ResidueProperties]:
    residue_atoms: Dict[str, List[AtomRecord]] = {}
    for atom in atoms:
        residue_atoms.setdefault(atom.residue_key, []).append(atom)

    props: Dict[str, ResidueProperties] = {}
    for residue_key, residue_atom_list in residue_atoms.items():
        anchor = residue_atom_list[0]
        res_name = (anchor.res_name or "").upper()
        non_polymer = (
            res_name not in POLYMER_RESIDUES
            and res_name not in WATER_RESIDUES
            and res_name not in METAL_ELEMENTS
        )
        is_nucleic = res_name in NUCLEIC_RESIDUES
        aromatic = res_name in AROMATIC_RESIDUES or is_nucleic
        acidic = res_name in ACIDIC_RESIDUES
        basic = res_name in BASIC_RESIDUES
        hydrophobic = res_name in HYDROPHOBIC_RESIDUES or is_nucleic

        counts = _residue_element_counts(residue_atom_list)
        if non_polymer:
            if not aromatic:
                aromatic = _is_likely_aromatic_nonpolymer(residue_atom_list)
            if not acidic:
                acidic = _is_likely_acidic_nonpolymer(residue_atom_list, counts)
            if not basic:
                basic = _is_likely_basic_nonpolymer(res_name, counts)
            if not hydrophobic:
                hydrophobic = _is_likely_hydrophobic_nonpolymer(counts, aromatic)

        props[residue_key] = ResidueProperties(
            aromatic=aromatic,
            acidic=acidic,
            basic=basic,
            hydrophobic=hydrophobic,
            non_polymer=non_polymer,
        )
    return props


def fetch_mmcif(pdb_id: str) -> str:
    url = MODEL_SERVER_URL.format(pdb_id=pdb_id.lower())
    try:
        with urllib.request.urlopen(url, timeout=20) as response:
            data = response.read()
    except Exception as exc:
        raise RuntimeError(f"Failed to fetch mmCIF for {pdb_id}: {exc}") from exc
    return data.decode("utf-8")


def list_chains(mmcif_text: str) -> Tuple[List[str], ChainAliases]:
    atoms, aliases = parse_mmcif_atoms(mmcif_text)
    chains = sorted({atom.chain_id for atom in atoms})
    return chains, aliases


def analyze_interface(
    structure_text: str,
    chain_a: str,
    chain_b: str,
    mode: str = "all",
    structure_format: str = "mmcif",
) -> dict:
    fmt = (structure_format or "mmcif").lower()
    if fmt == "pdb":
        atoms, aliases = parse_pdb_atoms(structure_text)
    else:
        atoms, aliases = parse_mmcif_atoms(structure_text)
    chain_a = aliases.normalize(chain_a)
    chain_b = aliases.normalize(chain_b)
    intrachain = chain_a == chain_b

    atoms_a = [atom for atom in atoms if atom.chain_id == chain_a]
    atoms_b = [atom for atom in atoms if atom.chain_id == chain_b]
    if not atoms_a or not atoms_b:
        raise ValueError("No atoms found for one or both chains")

    contacts = {
        "hydrogen_bonds": [],
        "salt_bridges": [],
        "hydrophobic": [],
        "metal_coordination": [],
        "pi_pi": [],
        "pi_cation": [],
        "other": [],
    }
    per_residue: Dict[str, dict] = {}
    residue_props = build_residue_properties(atoms)

    def bump_residue(atom: AtomRecord, key: str) -> None:
        entry = per_residue.setdefault(
            atom.residue_key,
            {
                "chain": atom.chain_id,
                "resName": atom.res_name,
                "seq": atom.res_seq,
                "hydrophobic": 0,
                "hbond": 0,
                "salt_bridge": 0,
                "metal_coordination": 0,
                "pi_pi": 0,
                "pi_cation": 0,
                "other": 0,
                "total": 0,
            },
        )
        entry[key] += 1
        entry["total"] += 1

    max_distance = 6.0
    cell_size = 5.0
    grid = build_grid(atoms_b, cell_size)
    atom_index_b = {id(atom): idx for idx, atom in enumerate(atoms_b)} if intrachain else {}
    for atom_a in atoms_a:
        if atom_a.element == "H":
            continue
        idx_a = atom_index_b.get(id(atom_a), -1) if intrachain else -1
        for atom_b in iter_neighbor_atoms(grid, atom_a, cell_size):
            if atom_b.element == "H":
                continue
            same_residue_pair = False
            if intrachain:
                idx_b = atom_index_b.get(id(atom_b), -1)
                if idx_b <= idx_a:
                    continue
                same_residue_pair = atom_a.residue_key == atom_b.residue_key
                if same_residue_pair:
                    # Most same-residue contacts are covalent/internal geometry noise.
                    # Keep only potential metal-ligand pairs for same-residue cofactors
                    # (e.g. heme Fe-N coordination), then enforce exact class below.
                    if is_metal_atom(atom_a) == is_metal_atom(atom_b):
                        continue
            dist = distance(atom_a, atom_b)
            if dist > max_distance:
                continue
            contact_type, strength = classify_contact(
                atom_a,
                atom_b,
                dist,
                residue_props.get(atom_a.residue_key),
                residue_props.get(atom_b.residue_key),
            )
            if contact_type is None:
                continue
            if same_residue_pair and contact_type != "metal_coordination":
                continue
            if (
                intrachain
                and contact_type == "hydrogen_bonds"
                and _is_adjacent_backbone_oxygen_nitrogen_pair(atom_a, atom_b)
            ):
                continue
            if intrachain and dist <= 2.2 and contact_type != "metal_coordination":
                # Keep intrachain metal coordination (typically short) while
                # continuing to suppress near-covalent short contacts.
                continue
            record = {
                "residueA": residue_payload(atom_a),
                "residueB": residue_payload(atom_b),
                "distance": round(dist, 3),
                "type": contact_type,
            }
            if contact_type == "metal_coordination":
                metal_side, metal_element = resolve_metal_contact_side(atom_a, atom_b)
                if metal_side:
                    record["metalSide"] = metal_side
                if metal_element:
                    record["metalElement"] = metal_element
            if strength:
                record["strength"] = strength

            contacts[contact_type].append(record)
            if contact_type == "hydrogen_bonds":
                bump_residue(atom_a, "hbond")
                bump_residue(atom_b, "hbond")
            elif contact_type == "salt_bridges":
                bump_residue(atom_a, "salt_bridge")
                bump_residue(atom_b, "salt_bridge")
            elif contact_type == "metal_coordination":
                bump_residue(atom_a, "metal_coordination")
                bump_residue(atom_b, "metal_coordination")
            else:
                bump_residue(atom_a, contact_type)
                bump_residue(atom_b, contact_type)

    mode = (mode or "all").lower()
    if mode != "all":
        contacts = filter_contacts_by_mode(contacts, mode)

    interface_res_a = {key for key in per_residue if key.startswith(f"{chain_a}:")}
    interface_res_b = {key for key in per_residue if key.startswith(f"{chain_b}:")}
    buried_fraction = None
    if intrachain:
        total_res = len(interface_res_a)
        if total_res:
            buried_fraction = {chain_a: 1.0}
    else:
        total_res = len(interface_res_a) + len(interface_res_b)
        if total_res:
            buried_fraction = {
                chain_a: round(len(interface_res_a) / total_res, 3),
                chain_b: round(len(interface_res_b) / total_res, 3),
            }

    return {
        "chainA": chain_a,
        "chainB": chain_b,
        "analysisVersion": TOOL_VERSION,
        "contacts": contacts,
        "perResidue": per_residue,
        "interfaceArea": None,
        "buriedFraction": buried_fraction,
        "approxDeltaG": None,
        "meta": {
            "engine": "toy-distance",
            "scope": "intrachain" if intrachain else "interchain",
            "note": "Interface area and deltaG require PISA or similar tools.",
        },
    }


def filter_contacts_by_mode(contacts: dict, mode: str) -> dict:
    mode = mode.lower()
    if mode == "hydrophobic":
        return {"hydrophobic": contacts["hydrophobic"]}
    if mode in {"electrostatic", "ionic", "salt"}:
        return {"salt_bridges": contacts["salt_bridges"]}
    if mode in {"metal", "metal_coordination", "coordination"}:
        return {"metal_coordination": contacts["metal_coordination"]}
    if mode in {"hbond", "hbond_network", "hydrogen"}:
        return {"hydrogen_bonds": contacts["hydrogen_bonds"]}
    if mode in {"aromatic", "pi"}:
        return {
            "pi_pi": contacts["pi_pi"],
            "pi_cation": contacts["pi_cation"],
        }
    return contacts


def residue_payload(atom: AtomRecord) -> dict:
    return {
        "chain": atom.chain_id,
        "resName": atom.res_name,
        "seq": atom.res_seq,
        "atom": atom.atom_name,
    }


def _prop(props: Optional[ResidueProperties], name: str, fallback: bool = False) -> bool:
    if not props:
        return fallback
    return bool(getattr(props, name, fallback))


def _is_cationic_atom_like(atom: AtomRecord, residue_basic: bool) -> bool:
    element = (atom.element or "").upper()
    if residue_basic:
        return element in {"N", "C", "S"}
    return element == "N"


def _is_anionic_atom_like(atom: AtomRecord, residue_acidic: bool) -> bool:
    element = (atom.element or "").upper()
    if element not in {"O", "S"}:
        return False
    if residue_acidic:
        return True
    atom_name = (atom.atom_name or "").strip().upper()
    return atom_name.startswith("O") or atom_name.startswith("S")


def _parse_res_seq_index(res_seq: str) -> Optional[int]:
    token = str(res_seq or "").strip()
    if not token:
        return None
    match = re.match(r"^-?\d+", token)
    if not match:
        return None
    try:
        return int(match.group(0))
    except Exception:
        return None


def _is_adjacent_backbone_oxygen_nitrogen_pair(atom_a: AtomRecord, atom_b: AtomRecord) -> bool:
    if not atom_a or not atom_b:
        return False
    if atom_a.chain_id != atom_b.chain_id:
        return False
    if atom_a.res_name not in STANDARD_AMINO_RESIDUES or atom_b.res_name not in STANDARD_AMINO_RESIDUES:
        return False
    seq_a = _parse_res_seq_index(atom_a.res_seq)
    seq_b = _parse_res_seq_index(atom_b.res_seq)
    if seq_a is None or seq_b is None:
        return False
    if abs(seq_a - seq_b) != 1:
        return False

    atom_name_a = (atom_a.atom_name or "").strip().upper()
    atom_name_b = (atom_b.atom_name or "").strip().upper()
    # Exclude peptide-linkage O···N neighbors between consecutive residues.
    # These are bonded-geometry contacts, not noncovalent hydrogen bonds.
    return {atom_name_a, atom_name_b} == {"O", "N"}


def classify_contact(
    atom_a: AtomRecord,
    atom_b: AtomRecord,
    dist: float,
    props_a: Optional[ResidueProperties] = None,
    props_b: Optional[ResidueProperties] = None,
) -> Tuple[Optional[str], Optional[str]]:
    res_a = atom_a.res_name
    res_b = atom_b.res_name
    element_a = atom_a.element
    element_b = atom_b.element

    is_acidic_a = _prop(props_a, "acidic", res_a in ACIDIC_RESIDUES)
    is_acidic_b = _prop(props_b, "acidic", res_b in ACIDIC_RESIDUES)
    is_basic_a = _prop(props_a, "basic", res_a in BASIC_RESIDUES)
    is_basic_b = _prop(props_b, "basic", res_b in BASIC_RESIDUES)
    is_aromatic_a = _prop(props_a, "aromatic", res_a in AROMATIC_RESIDUES)
    is_aromatic_b = _prop(props_b, "aromatic", res_b in AROMATIC_RESIDUES)
    is_hydrophobic_a = _prop(props_a, "hydrophobic", res_a in HYDROPHOBIC_RESIDUES)
    is_hydrophobic_b = _prop(props_b, "hydrophobic", res_b in HYDROPHOBIC_RESIDUES)
    has_nonpolymer = _prop(props_a, "non_polymer", False) or _prop(props_b, "non_polymer", False)
    is_metal_a = is_metal_atom(atom_a)
    is_metal_b = is_metal_atom(atom_b)

    if is_metal_a != is_metal_b:
        metal_atom = atom_a if is_metal_a else atom_b
        ligand_atom = atom_b if is_metal_a else atom_a
        coordination_cutoff = metal_coordination_cutoff(metal_atom.element)
        if ligand_atom.element in METAL_DONOR_ELEMENTS and dist <= coordination_cutoff:
            strength = "strong" if dist <= coordination_cutoff - 0.35 else "weak"
            return "metal_coordination", strength

    if (
        ((is_acidic_a and is_basic_b) or (is_basic_a and is_acidic_b))
        and dist <= 4.0
        and {element_a, element_b} <= {"O", "N", "S"}
    ):
        return "salt_bridges", None

    if has_nonpolymer and dist <= 4.6:
        if (
            (is_acidic_a and _is_cationic_atom_like(atom_b, is_basic_b))
            or (is_acidic_b and _is_cationic_atom_like(atom_a, is_basic_a))
            or (is_basic_a and _is_anionic_atom_like(atom_b, is_acidic_b))
            or (is_basic_b and _is_anionic_atom_like(atom_a, is_acidic_a))
        ):
            return "salt_bridges", None

    if element_a in {"O", "N"} and element_b in {"O", "N"} and dist <= 3.5:
        strength = "strong" if dist <= 3.0 else "weak"
        return "hydrogen_bonds", strength

    if is_aromatic_a and is_aromatic_b and dist <= 5.0:
        return "pi_pi", None

    if (is_aromatic_a and is_basic_b) or (is_aromatic_b and is_basic_a):
        if dist <= 6.0:
            return "pi_cation", None

    if has_nonpolymer and dist <= 6.0:
        if (is_aromatic_a and _is_cationic_atom_like(atom_b, is_basic_b)) or (
            is_aromatic_b and _is_cationic_atom_like(atom_a, is_basic_a)
        ):
            return "pi_cation", None

    if (
        is_hydrophobic_a
        and is_hydrophobic_b
        and {element_a, element_b} <= {"C", "S"}
        and dist <= 4.5
    ):
        return "hydrophobic", None

    # Non-polymer-specific fallback: allow aromatic carbon-rich ligands to
    # contribute hydrophobic contacts even if the whole residue is not globally
    # classified hydrophobic due phosphate/charged substituents.
    if has_nonpolymer and {element_a, element_b} <= {"C", "S"} and dist <= 4.8:
        aromatic_hydrophobic_mix = (
            (is_aromatic_a and is_hydrophobic_b)
            or (is_aromatic_b and is_hydrophobic_a)
            or (is_aromatic_a and is_aromatic_b)
        )
        if aromatic_hydrophobic_mix:
            return "hydrophobic", None

    if dist <= 4.0:
        return "other", None

    return None, None


def parse_pdb_atoms(pdb_text: str) -> Tuple[List[AtomRecord], ChainAliases]:
    atoms: List[AtomRecord] = []
    for line in pdb_text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            alt_loc = line[16:17].strip()
            if alt_loc and alt_loc not in {"A", "1"}:
                continue
            chain_id = line[21:22].strip() or "A"
            res_name = line[17:20].strip().upper()
            res_seq = line[22:26].strip()
            ins_code = line[26:27].strip()
            if ins_code:
                res_seq = f"{res_seq}{ins_code}"
            atom_name = line[12:16].strip()
            element = line[76:78].strip().upper()
            if not element:
                element = guess_element(atom_name)
            if element == "H":
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except Exception:
            continue

        atoms.append(
            AtomRecord(
                chain_id=chain_id,
                chain_label=chain_id,
                res_name=res_name,
                res_seq=res_seq,
                atom_name=atom_name,
                element=element,
                x=x,
                y=y,
                z=z,
            )
        )

    # PDB has no label/auth chain alias distinction in this flow.
    return atoms, ChainAliases(label_to_auth={})


def is_metal_atom(atom: AtomRecord) -> bool:
    if not atom:
        return False
    element = (atom.element or "").upper()
    if element in METAL_ELEMENTS:
        return True
    # Some deposits encode ions using residue names more reliably than element symbols.
    res_name = (atom.res_name or "").upper().replace("+", "").replace("-", "")
    return res_name in METAL_ELEMENTS


def metal_coordination_cutoff(element: str) -> float:
    key = (element or "").upper()
    return METAL_COORDINATION_CUTOFF.get(key, 2.85)


def resolve_metal_contact_side(atom_a: AtomRecord, atom_b: AtomRecord) -> Tuple[Optional[str], Optional[str]]:
    is_metal_a = is_metal_atom(atom_a)
    is_metal_b = is_metal_atom(atom_b)
    if is_metal_a and not is_metal_b:
        return "A", (atom_a.element or "").upper() or (atom_a.res_name or "").upper()
    if is_metal_b and not is_metal_a:
        return "B", (atom_b.element or "").upper() or (atom_b.res_name or "").upper()
    return None, None


def build_grid(atoms: Iterable[AtomRecord], cell_size: float) -> Dict[Tuple[int, int, int], List[AtomRecord]]:
    grid: Dict[Tuple[int, int, int], List[AtomRecord]] = {}
    for atom in atoms:
        key = grid_key(atom, cell_size)
        grid.setdefault(key, []).append(atom)
    return grid


def grid_key(atom: AtomRecord, cell_size: float) -> Tuple[int, int, int]:
    return (
        int(math.floor(atom.x / cell_size)),
        int(math.floor(atom.y / cell_size)),
        int(math.floor(atom.z / cell_size)),
    )


def iter_neighbor_atoms(
    grid: Dict[Tuple[int, int, int], List[AtomRecord]],
    atom: AtomRecord,
    cell_size: float,
) -> Iterable[AtomRecord]:
    base = grid_key(atom, cell_size)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                key = (base[0] + dx, base[1] + dy, base[2] + dz)
                for neighbor in grid.get(key, []):
                    yield neighbor


def distance(atom_a: AtomRecord, atom_b: AtomRecord) -> float:
    dx = atom_a.x - atom_b.x
    dy = atom_a.y - atom_b.y
    dz = atom_a.z - atom_b.z
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def parse_mmcif_atoms(mmcif_text: str) -> Tuple[List[AtomRecord], ChainAliases]:
    lines = mmcif_text.splitlines()
    columns: List[str] = []
    data_rows: List[List[str]] = []
    in_atom_site = False

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "loop_":
            j = i + 1
            cols: List[str] = []
            while j < len(lines):
                col_line = lines[j].strip()
                if col_line.startswith("_atom_site."):
                    cols.append(col_line.split()[0])
                    j += 1
                    continue
                break
            if cols:
                columns = cols
                data_rows = list(iter_loop_rows(lines, j, len(columns)))
                in_atom_site = True
                break
            i = j
            continue
        i += 1

    if not in_atom_site:
        return [], ChainAliases(label_to_auth={})

    col_index = {col: idx for idx, col in enumerate(columns)}
    def idx(*names: str) -> Optional[int]:
        for name in names:
            if name in col_index:
                return col_index[name]
        return None

    chain_idx = idx("_atom_site.auth_asym_id", "_atom_site.label_asym_id")
    label_chain_idx = idx("_atom_site.label_asym_id")
    auth_res_seq_idx = idx("_atom_site.auth_seq_id")
    label_res_seq_idx = idx("_atom_site.label_seq_id")
    res_name_idx = idx("_atom_site.auth_comp_id", "_atom_site.label_comp_id")
    atom_name_idx = idx("_atom_site.auth_atom_id", "_atom_site.label_atom_id")
    element_idx = idx("_atom_site.type_symbol")
    x_idx = idx("_atom_site.Cartn_x")
    y_idx = idx("_atom_site.Cartn_y")
    z_idx = idx("_atom_site.Cartn_z")
    alt_idx = idx("_atom_site.label_alt_id", "_atom_site.pdbx_PDB_alt_id")
    ins_code_idx = idx("_atom_site.pdbx_PDB_ins_code")

    if (
        None in (chain_idx, res_name_idx, atom_name_idx, x_idx, y_idx, z_idx)
        or (auth_res_seq_idx is None and label_res_seq_idx is None)
    ):
        return [], ChainAliases(label_to_auth={})

    atoms: List[AtomRecord] = []
    label_to_auth: Dict[str, str] = {}

    for row in data_rows:
        try:
            chain_id = row[chain_idx]
            if chain_id in {"?", "."}:
                continue
            res_seq = ""
            if auth_res_seq_idx is not None:
                res_seq = str(row[auth_res_seq_idx]).strip()
            if (not res_seq or res_seq in {".", "?"}) and label_res_seq_idx is not None:
                res_seq = str(row[label_res_seq_idx]).strip()
            if not res_seq or res_seq in {".", "?"}:
                continue
            if ins_code_idx is not None:
                ins_code = str(row[ins_code_idx]).strip()
                if ins_code and ins_code not in {".", "?"}:
                    res_seq = f"{res_seq}{ins_code}"
            res_name = row[res_name_idx].upper()
            atom_name = row[atom_name_idx]
            element = row[element_idx].upper() if element_idx is not None else ""
            if element in {"", "?", "."}:
                element = guess_element(atom_name)
            if element == "H":
                continue
            if alt_idx is not None:
                alt_id = row[alt_idx]
                if alt_id not in {".", "?", "A", "1"}:
                    continue
            x = float(row[x_idx])
            y = float(row[y_idx])
            z = float(row[z_idx])
        except Exception:
            continue

        label_chain = row[label_chain_idx] if label_chain_idx is not None else chain_id
        if (
            label_chain not in {".", "?"}
            and chain_id not in {".", "?"}
            and label_chain != chain_id
            and str(label_chain).lower() != str(chain_id).lower()
        ):
            label_to_auth.setdefault(label_chain, chain_id)

        atoms.append(
            AtomRecord(
                chain_id=chain_id,
                chain_label=label_chain,
                res_name=res_name,
                res_seq=res_seq,
                atom_name=atom_name,
                element=element,
                x=x,
                y=y,
                z=z,
            )
        )

    return atoms, ChainAliases(label_to_auth=label_to_auth)


def iter_loop_rows(lines: List[str], start: int, columns: int) -> Iterable[List[str]]:
    tokens: List[str] = []
    i = start
    while i < len(lines):
        raw = lines[i]
        stripped = raw.strip()
        if not stripped:
            i += 1
            continue
        if stripped.startswith("#"):
            i += 1
            continue
        if stripped.startswith("_") or stripped.startswith("loop_") or stripped.startswith("data_"):
            break
        if stripped.startswith(";"):
            block, i = read_semicolon_block(lines, i)
            tokens.append(block)
        else:
            tokens.extend(shlex.split(raw, posix=True))
            i += 1
        while len(tokens) >= columns:
            row = tokens[:columns]
            tokens = tokens[columns:]
            yield row
    return


def read_semicolon_block(lines: List[str], start: int) -> Tuple[str, int]:
    block_lines = []
    i = start + 1
    while i < len(lines):
        line = lines[i]
        if line.startswith(";"):
            return "\n".join(block_lines), i + 1
        block_lines.append(line.rstrip("\n"))
        i += 1
    return "\n".join(block_lines), i


def guess_element(atom_name: str) -> str:
    name = atom_name.strip()
    if not name:
        return ""
    for char in name:
        if char.isalpha():
            return char.upper()
    return ""


def cache_key(pdb_id: Optional[str], mmcif_text: Optional[str], chain_a: str, chain_b: str, mode: str) -> str:
    if pdb_id:
        source = pdb_id.lower()
    else:
        digest = hashlib.sha256((mmcif_text or "").encode("utf-8")).hexdigest()
        source = digest[:12]
    return f"{source}:{chain_a}:{chain_b}:{mode}:{TOOL_VERSION}"
