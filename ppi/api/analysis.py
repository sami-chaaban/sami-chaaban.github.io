from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple
import hashlib
import math
import shlex
import urllib.request

TOOL_VERSION = "toy-0.4"
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
WATER_RESIDUES = {"HOH", "WAT", "H2O"}
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
    "MG": 2.65,
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
        return self.label_to_auth.get(chain_id, chain_id)


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
            if intrachain:
                idx_b = atom_index_b.get(id(atom_b), -1)
                if idx_b <= idx_a:
                    continue
                # Exclude same-residue atom pairs and covalent-neighbor distances.
                if atom_a.residue_key == atom_b.residue_key:
                    continue
            dist = distance(atom_a, atom_b)
            if dist > max_distance:
                continue
            contact_type, strength = classify_contact(atom_a, atom_b, dist)
            if contact_type is None:
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


def classify_contact(
    atom_a: AtomRecord, atom_b: AtomRecord, dist: float
) -> Tuple[Optional[str], Optional[str]]:
    res_a = atom_a.res_name
    res_b = atom_b.res_name
    element_a = atom_a.element
    element_b = atom_b.element

    is_acidic_a = res_a in ACIDIC_RESIDUES
    is_acidic_b = res_b in ACIDIC_RESIDUES
    is_basic_a = res_a in BASIC_RESIDUES
    is_basic_b = res_b in BASIC_RESIDUES
    is_aromatic_a = res_a in AROMATIC_RESIDUES
    is_aromatic_b = res_b in AROMATIC_RESIDUES
    is_hydrophobic_a = res_a in HYDROPHOBIC_RESIDUES
    is_hydrophobic_b = res_b in HYDROPHOBIC_RESIDUES
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

    if element_a in {"O", "N"} and element_b in {"O", "N"} and dist <= 3.5:
        strength = "strong" if dist <= 3.0 else "weak"
        return "hydrogen_bonds", strength

    if is_aromatic_a and is_aromatic_b and dist <= 5.0:
        return "pi_pi", None

    if (is_aromatic_a and is_basic_b) or (is_aromatic_b and is_basic_a):
        if dist <= 6.0:
            return "pi_cation", None

    if (
        is_hydrophobic_a
        and is_hydrophobic_b
        and {element_a, element_b} <= {"C", "S"}
        and dist <= 4.5
    ):
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
            if res_name in WATER_RESIDUES:
                continue
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
    res_seq_idx = idx("_atom_site.auth_seq_id", "_atom_site.label_seq_id")
    res_name_idx = idx("_atom_site.auth_comp_id", "_atom_site.label_comp_id")
    atom_name_idx = idx("_atom_site.auth_atom_id", "_atom_site.label_atom_id")
    element_idx = idx("_atom_site.type_symbol")
    x_idx = idx("_atom_site.Cartn_x")
    y_idx = idx("_atom_site.Cartn_y")
    z_idx = idx("_atom_site.Cartn_z")
    alt_idx = idx("_atom_site.label_alt_id", "_atom_site.pdbx_PDB_alt_id")

    if None in (chain_idx, res_seq_idx, res_name_idx, atom_name_idx, x_idx, y_idx, z_idx):
        return [], ChainAliases(label_to_auth={})

    atoms: List[AtomRecord] = []
    label_to_auth: Dict[str, str] = {}

    for row in data_rows:
        try:
            chain_id = row[chain_idx]
            if chain_id in {"?", "."}:
                continue
            res_seq = row[res_seq_idx]
            res_name = row[res_name_idx].upper()
            if res_name in WATER_RESIDUES:
                continue
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
        if label_chain not in {".", "?"} and chain_id not in {".", "?"}:
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
