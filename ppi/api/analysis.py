from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple
import hashlib
import copy
import json
import math
import os
import re
import tempfile
import urllib.request

from .cache import ReportCache

ARPEGGIO_IMPORT_ERROR: Optional[Exception] = None
try:
    from arpeggio.core import InteractionComplex
except Exception as exc:  # pragma: no cover - optional runtime dependency
    InteractionComplex = None
    ARPEGGIO_IMPORT_ERROR = exc

GEMMI_IMPORT_ERROR: Optional[Exception] = None
try:
    import gemmi as _gemmi  # type: ignore
except Exception as exc:  # pragma: no cover - optional runtime dependency
    _gemmi = None
    GEMMI_IMPORT_ERROR = exc

# Open Babel can emit repetitive aromatic kekulization warnings for some ligands.
# Keep stderr focused on actionable failures.
try:  # pragma: no cover - optional runtime dependency side effect
    from openbabel import openbabel as _openbabel  # type: ignore
except Exception:
    _openbabel = None
else:
    try:
        _openbabel.obErrorLog.SetOutputLevel(getattr(_openbabel, "obError", 0))
    except Exception:
        pass


TOOL_VERSION = "roami-assertion-2.4"
MODEL_SERVER_URL = (
    "https://www.ebi.ac.uk/pdbe/model-server/v1/{pdb_id}/full"
    "?encoding=cif&data_source=pdb-h"
)
ARPEGGIO_INTERACTING_CUTOFF = 6.0
ARPEGGIO_VDW_COMP = 0.1
ARPEGGIO_INTERFACE_SELECTION_CUTOFF = 8.0
ARPEGGIO_HYDROPHOBIC_MAX_DISTANCE = 4.6

CONTACT_BUCKET_TO_CATEGORY = {
    "hydrogen_bonds": "hbond",
    "polar_contacts": "polar_contact",
    "base_pairing": "base_pairing",
    "salt_bridges": "salt_bridge",
    "halogen_bonds": "halogen_bond",
    "hydrophobic": "hydrophobic",
    "metal_coordination": "metal_coordination",
    "pi_pi": "pi_pi",
    "pi_cation": "pi_cation",
    "aromatic_packing": "aromatic_packing",
    "other": "other",
}
CONTACT_CATEGORY_TO_BUCKET = {
    value: key for key, value in CONTACT_BUCKET_TO_CATEGORY.items()
}
CONTACT_CATEGORY_TO_PER_RESIDUE_KEY = {
    "hbond": "hbond",
    "polar_contact": "polar_contact",
    "polar_proximal": "other",
    "base_pairing": "base_pairing",
    "salt_bridge": "salt_bridge",
    "halogen_bond": "halogen_bond",
    "hydrophobic": "hydrophobic",
    "metal_coordination": "metal_coordination",
    "aromatic_packing": "aromatic_packing",
    "aromatic_proximal": "other",
    "pi_pi": "pi_pi",
    "pi_cation": "pi_cation",
    "packing_contact": "other",
    "proximal": "other",
    "vdw": "vdw",
    "clash": "clash",
    "invalid_contact": "other",
    "other": "other",
}
ASSERTED_FAMILY_TO_BUCKET = {
    "hbond": "hydrogen_bonds",
    "polar_contact": "polar_contacts",
    "polar_proximal": "other",
    "base_pairing": "base_pairing",
    "salt_bridge": "salt_bridges",
    "halogen_bond": "halogen_bonds",
    "hydrophobic": "hydrophobic",
    "metal_coordination": "metal_coordination",
    "aromatic_packing": "aromatic_packing",
    "aromatic_proximal": "other",
    "pi_pi": "pi_pi",
    "pi_cation": "pi_cation",
    "packing_contact": "other",
    "proximal": "other",
    "vdw": "other",
    "clash": "other",
    "invalid_contact": "other",
    "covalent_bond": "other",
    "other": "other",
}

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
POLAR_OR_CHARGED_SIDECHAIN_RESIDUES = {
    "ASP", "GLU", "ASN", "GLN", "ARG", "LYS", "HIS", "SER", "THR", "TYR", "CYS",
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
METAL_CONTACT_TERMS = {"METAL", "METAL_COMPLEX"}
HBOND_EXPLICIT_TERMS = {"HBOND", "WEAK_HBOND"}
HBOND_POLAR_FALLBACK_TERMS = {"POLAR", "WEAK_POLAR", "CARBONYL"}
HBOND_POLAR_CONTACT_TERMS = HBOND_EXPLICIT_TERMS | HBOND_POLAR_FALLBACK_TERMS
POLAR_CONTACT_MAX_DISTANCE = 3.8
HBOND_POLAR_FALLBACK_MAX_DISTANCE = POLAR_CONTACT_MAX_DISTANCE
HBOND_EXPLICIT_MAX_DISTANCE = 3.7
HBOND_CANDIDATE_MAX_DISTANCE = 3.6
HBOND_CANDIDATE_MEDIUM_CONFIDENCE_MAX_DISTANCE = 3.2
HBOND_STRONG_ANGLE_MIN = 150.0
HBOND_PROXY_ANGLE_MIN = 125.0
HBOND_PROXY_STRONG_ANGLE_MIN = 150.0
HBOND_PROXY_FAIL_ANGLE_MAX = 110.0
HBOND_UNUSUALLY_SHORT_DISTANCE = 2.6
HBOND_EXTREME_SHORT_DISTANCE = 2.4
HBOND_BORDERLINE_DISTANCE = 3.5
HBOND_HEAVY_MIN_DISTANCE = 1.45
HBOND_MAX_DISTANCE_OO = 3.55
HBOND_MAX_DISTANCE_NO = 3.75
HBOND_MAX_DISTANCE_NN = 3.5
HBOND_MAX_DISTANCE_CHALCOGEN = 3.7
HALOGEN_BOND_EXPLICIT_TERMS = {"HALOGEN", "XBOND", "X_BOND", "HALOGEN_BOND"}
HALOGEN_BOND_DONOR_ELEMENTS = {"CL", "BR", "I"}
HALOGEN_BOND_ACCEPTOR_ELEMENTS = {"N", "O", "S", "SE"}
HALOGEN_BOND_DISTANCE_CUTOFF_BY_ELEMENT = {
    "CL": 3.5,
    "BR": 3.7,
    "I": 3.9,
}
HALOGEN_BOND_DEFAULT_DISTANCE_CUTOFF = 3.8
HALOGEN_BOND_STRONG_ANGLE_MIN = 155.0
HALOGEN_BOND_MEDIUM_ANGLE_MIN = 145.0
HALOGEN_DONOR_BOND_MIN_DISTANCE = 1.2
HALOGEN_DONOR_BOND_MAX_DISTANCE = 2.3
SALT_BRIDGE_MAX_DISTANCE = 4.8
SALT_BRIDGE_CONFIDENT_DISTANCE = 4.2
HYDROPHOBIC_MIN_DISTANCE = 2.2
HYDROPHOBIC_MEDIUM_MAX_DISTANCE = 4.2
HYDROPHOBIC_MAX_ALLOWED_OVERLAP = 0.5
POLAR_PRECEDENCE_MAX_DISTANCE = 4.0
CLASH_SOFT_OVERLAP = 0.25
CLASH_HARD_OVERLAP = 0.55
CLASH_HARD_OVERLAP_POLAR = 0.7
SOFT_CLASH_PRECLASSIFY = 0.4
HARD_CLASH_PRECLASSIFY = 0.7
INVALID_NONBONDED_MIN_DISTANCE = 0.8
INVALID_NONBONDED_STRICT_DISTANCE = 1.0
PACKING_CONTACT_TOLERANCE = 0.5
PACKING_CONTACT_MAX_DISTANCE = 4.6
PROXIMAL_CONTACT_MAX_DISTANCE = 8.0
PACKING_ELIGIBLE_ELEMENTS = {"C", "S", "SE", "F", "CL", "BR", "I"}
BASE_PAIR_CANDIDATE_MAX_DISTANCE = 3.6
BASE_PAIR_MIN_POLAR_PAIR_SUPPORT = 2
BASE_PAIR_SINGLE_PAIR_LOW_CONFIDENCE_DISTANCE = 3.5
BASE_PAIR_RING_PLANE_MAX_NORMAL_ANGLE = 25.0
BASE_PAIR_RING_PLANE_MAX_INTERPLANAR_DISTANCE = 0.7
BASE_PAIR_RING_PLANE_MIN_LATERAL_OFFSET = 3.5
BASE_PAIR_RING_PLANE_MAX_LATERAL_OFFSET = 6.5
BASE_PAIR_DISTANCE_STRONG_MAX = 3.1
BASE_PAIR_DISTANCE_MEDIUM_MAX = 3.3
BASE_PAIR_SCORE_WEIGHT_SUPPORT = 0.4
BASE_PAIR_SCORE_WEIGHT_DISTANCE = 0.22
BASE_PAIR_SCORE_WEIGHT_ANGLE = 0.16
BASE_PAIR_SCORE_WEIGHT_MUTUAL_BEST = 0.12
BASE_PAIR_SCORE_WEIGHT_COPLANARITY = 0.1
AROMATIC_PACKING_MIN_DISTANCE = 3.2
AROMATIC_PACKING_MAX_DISTANCE = 4.2
PI_PI_MIN_CENTROID_DISTANCE = 3.3
PI_PI_MAX_CENTROID_DISTANCE = 6.2
PI_PI_MIN_INTERPLANAR_DISTANCE = 2.6
PI_PI_MAX_INTERPLANAR_DISTANCE = 4.3
PI_PI_MAX_LATERAL_OFFSET = 3.2
PI_PI_STACKED_MAX_NORMAL_ANGLE = 30.0
PI_PI_TSHAPED_MIN_NORMAL_ANGLE = 60.0
PI_PI_TSHAPED_MAX_LATERAL_OFFSET = 4.0
AROMATIC_ASSERTED_FAMILIES = {"pi_pi", "pi_cation", "aromatic_packing", "aromatic_proximal"}
AROMATIC_TOP_K_PER_RESIDUE_PAIR_BY_FAMILY = {
    "pi_pi": 2,
    "pi_cation": 2,
    "aromatic_packing": 3,
    "aromatic_proximal": 2,
}
COVALENT_DISTANCE_MAX_PO = 1.9
DEFAULT_VDW_RADIUS = 1.7
VDW_RADIUS_BY_ELEMENT = {
    "H": 1.2,
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.8,
    "S": 1.8,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "SE": 1.9,
    "ZN": 1.39,
    "MG": 1.73,
    "CA": 1.94,
    "NA": 2.27,
    "K": 2.75,
}
MIN_NONBONDED_DISTANCE_BY_ELEMENT_PAIR = {
    frozenset({"N", "N"}): 2.4,
    frozenset({"N", "O"}): 2.3,
    frozenset({"O", "O"}): 2.4,
    frozenset({"C", "C"}): 2.8,
    frozenset({"C", "N"}): 2.7,
    frozenset({"C", "O"}): 2.7,
    frozenset({"C", "CL"}): 3.0,
    frozenset({"C", "S"}): 2.9,
    frozenset({"CL", "CL"}): 3.3,
}
HYDROPHOBIC_MIN_DISTANCE_BY_ELEMENT_PAIR = {
    frozenset({"C", "C"}): 2.9,
    frozenset({"C", "CL"}): 3.0,
    frozenset({"CL", "CL"}): 3.3,
    frozenset({"C", "S"}): 2.9,
}
ALLOW_CRYSTAL_CONTACTS = False
ATOM_REUSE_CONFIDENCE_THRESHOLD_BY_BUCKET = {
    "hydrogen_bonds": 3,
    "halogen_bonds": 1,
    "hydrophobic": 4,
}
POLAR_CONTACT_ELEMENTS = {"N", "O", "S", "SE"}
PROTEIN_SIDECHAIN_CARBONYL_CARBON_BY_RESIDUE = {
    "ASP": {"CG"},
    "GLU": {"CD"},
    "ASN": {"CG"},
    "GLN": {"CD"},
    # Ambiguous residue aliases in some structures.
    "ASX": {"CG"},
    "GLX": {"CD"},
}
PROTEIN_CHARGED_GROUP_ASSOCIATED_CARBON_BY_RESIDUE = {
    # Guanidinium central carbon.
    "ARG": {"CZ"},
}
HBOND_DONOR_ANTECEDENT_BY_RESIDUE = {
    "SER": {"OG": "CB"},
    "THR": {"OG1": "CB"},
    "TYR": {"OH": "CZ"},
}
NUCLEOBASE_ATOMS_BY_FAMILY = {
    "A": {"N1", "C2", "N3", "C4", "C5", "C6", "N6", "N7", "C8", "N9"},
    "G": {"N1", "C2", "N2", "N3", "C4", "C5", "C6", "O6", "N7", "C8", "N9"},
    "C": {"N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"},
    "U": {"N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"},
    "T": {"N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6", "C7"},
    "I": {"N1", "C2", "O2", "N3", "C4", "O6", "C5", "C6", "N7", "C8", "N9"},
}
NUCLEOBASE_PAIRING_EDGE_ATOMS_BY_FAMILY = {
    "A": {"N1", "N6", "N7"},
    "G": {"N1", "N2", "O6"},
    "C": {"N3", "N4", "O2"},
    "U": {"N3", "O4", "O2"},
    "T": {"N3", "O4", "O2"},
    "I": {"N1", "O6"},
}
NUCLEOBASE_GLYCOSIDIC_ATOM_BY_FAMILY = {
    "A": "N9",
    "G": "N9",
    "C": "N1",
    "U": "N1",
    "T": "N1",
    "I": "N9",
}
NUCLEOTIDE_BACKBONE_ATOMS = {
    "P",
    "OP1",
    "OP2",
    "OP3",
    "O1P",
    "O2P",
    "O3P",
    "O5'",
    "C5'",
    "C4'",
    "O4'",
    "C3'",
    "O3'",
    "C2'",
    "O2'",
    "C1'",
    "O5*",
    "C5*",
    "C4*",
    "O4*",
    "C3*",
    "O3*",
    "C2*",
    "O2*",
    "C1*",
}
NUCLEOTIDE_LINKAGE_PHOSPHATE_ATOMS = {
    "P",
    "OP1",
    "OP2",
    "OP3",
    "O1P",
    "O2P",
    "O3P",
    "O5'",
    "O5*",
}
NUCLEOTIDE_LINKAGE_O3_ATOMS = {"O3'", "O3*"}
PI_NON_RING_SURFACE_TERMS = {"DONORPI", "CARBONPI", "HALOGENPI", "METSULPHURPI"}
AROMATIC_RING_ATOMS_BY_RESIDUE = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HID": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HIE": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HIP": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HSD": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HSE": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "HSP": {"CG", "ND1", "CD2", "CE1", "NE2"},
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
}
NUCLEOBASE_RING_ATOMS_BY_FAMILY = {
    "A": {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"},
    "G": {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"},
    "C": {"N1", "C2", "N3", "C4", "C5", "C6"},
    "U": {"N1", "C2", "N3", "C4", "C5", "C6"},
    "T": {"N1", "C2", "N3", "C4", "C5", "C6"},
    "I": {"N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"},
}
NUCLEIC_A_FAMILY = {"A", "RA", "DA", "ADE", "AMP", "ADP", "ATP", "ANP"}
NUCLEIC_G_FAMILY = {"G", "RG", "DG", "GUA", "GMP", "GDP", "GTP"}
NUCLEIC_C_FAMILY = {"C", "RC", "DC", "CYT", "CMP", "CDP", "CTP"}
NUCLEIC_U_FAMILY = {"U", "RU", "DU", "URA", "UMP", "UDP", "UTP"}
NUCLEIC_T_FAMILY = {"T", "RT", "DT", "THY", "TMP", "TDP", "TTP"}
NUCLEIC_I_FAMILY = {"I", "RI", "DI", "IMP", "IDP", "ITP"}
NUCLEIC_BASE_PAIR_FAMILY_BY_BASES = {
    frozenset({"A", "U"}): "AU",
    frozenset({"A", "T"}): "AT",
    frozenset({"G", "C"}): "GC",
}
CARBOXYLATE_ACCEPTOR_ATOMS_BY_RESIDUE = {
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
    "ASX": {"OD1", "OD2"},
    "GLX": {"OE1", "OE2"},
}
SALT_BRIDGE_CATION_SITE_ATOMS_BY_RESIDUE = {
    "ARG": ("NH1", "NH2", "NE"),
    "LYS": ("NZ",),
    "HIP": ("NE2", "ND1"),
    "HSP": ("NE2", "ND1"),
}
SALT_BRIDGE_ANION_SITE_ATOMS_BY_RESIDUE = {
    "ASP": ("OD1", "OD2"),
    "ASX": ("OD1", "OD2"),
    "GLU": ("OE1", "OE2"),
    "GLX": ("OE1", "OE2"),
}
PROTONATED_HISTIDINE_RESIDUE_NAMES = {"HIP", "HSP"}
SULFATE_LIKE_RESIDUE_NAMES = {"SO4", "HSO4", "SUL", "SULF"}
NUCLEIC_SUGAR_DONOR_OXYGENS = {"O2'", "O3'"}
NUCLEIC_SUGAR_ACCEPTOR_OXYGENS = {"O2'", "O3'", "O4'", "O5'"}
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
METAL_DEFAULT_COORDINATION_CUTOFF = 2.85


def _env_positive_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        parsed = int(str(raw).strip())
    except Exception:
        return default
    return parsed if parsed > 0 else default


PARSED_STRUCTURE_CACHE = ReportCache(
    ttl_seconds=_env_positive_int("ANALYZE_PARSED_STRUCTURE_CACHE_TTL_SECONDS", 60 * 60 * 6),
    max_entries=_env_positive_int("ANALYZE_PARSED_STRUCTURE_CACHE_MAX_ENTRIES", 12),
)
ARPEGGIO_CONTACTS_CACHE = ReportCache(
    ttl_seconds=_env_positive_int("ANALYZE_ARPEGGIO_CONTACT_CACHE_TTL_SECONDS", 60 * 60 * 6),
    max_entries=_env_positive_int("ANALYZE_ARPEGGIO_CONTACT_CACHE_MAX_ENTRIES", 64),
)


@dataclass(frozen=True)
class AtomRecord:
    chain_id: str
    chain_auth: str
    chain_label: str
    res_name: str
    res_seq: str
    atom_name: str
    element: str
    x: float
    y: float
    z: float
    altloc: str = ""
    occupancy: float = 1.0
    model_id: str = "1"

    @property
    def residue_key(self) -> str:
        return f"{self.chain_id}:{self.res_seq}"


class ResidueAtomIndex(dict):
    """Per-analysis immutable atom lookup and ring-feature caches."""
    def __init__(self, atoms=()):
        super().__init__()
        self.atom_lookup = {}
        self.ring_descriptors = {}
        self.ring_points = {}
        self.chemical_nodes = {}
        for atom in atoms:
            self.setdefault((atom.chain_id, atom.res_seq), []).append(atom)
            self.atom_lookup[(atom.chain_id, atom.res_seq, _normalize_atom_name(atom.atom_name))] = atom


@dataclass(frozen=True)
class ChainAliases:
    label_to_auth: Dict[str, str]
    auth_ids: Set[str]
    canonical_ids: Set[str]
    query_to_canonical: Dict[str, str]
    auth_to_canonical: Dict[str, Tuple[str, ...]]
    canonical_to_auth: Dict[str, str]
    auth_residue_to_canonical: Dict[Tuple[str, str, str], Tuple[str, ...]]
    auth_residue_atom_to_canonical: Dict[Tuple[str, str, str, str], Tuple[str, ...]]
    query_ids: Set[str]

    def normalize(self, chain_id: str) -> str:
        token = str(chain_id or "").strip()
        if not token:
            return token
        resolved = sorted(self.resolve(token))
        if len(resolved) == 1:
            return resolved[0]
        return token

    def resolve(self, chain_id: str) -> Set[str]:
        token = str(chain_id or "").strip()
        if not token:
            return set()
        if token in self.canonical_ids:
            return {token}
        mapped = str(self.query_to_canonical.get(token, "") or "").strip()
        if mapped:
            return {mapped}
        candidates = self.auth_to_canonical.get(token)
        if candidates:
            return set(candidates)
        if token in self.auth_ids:
            return {token}
        return set()

    def selection_ids(self, chain_id: str) -> Set[str]:
        resolved = self.resolve(chain_id)
        if resolved:
            return {
                str(self.canonical_to_auth.get(token, token) or "").strip() or token
                for token in resolved
            }
        token = str(chain_id or "").strip()
        if token in self.auth_ids:
            return {token}
        return set()

    def resolve_partner_chain(
        self,
        *,
        auth_chain: str = "",
        label_chain: str = "",
        res_seq: str = "",
        res_name: str = "",
        atom_name: str = "",
    ) -> str:
        label_token = str(label_chain or "").strip()
        if label_token and label_token not in {"?", "."}:
            resolved_label = sorted(self.resolve(label_token))
            if len(resolved_label) == 1:
                return resolved_label[0]

        auth_token = str(auth_chain or "").strip()
        seq_token = str(res_seq or "").strip()
        residue_token = str(res_name or "").strip().upper()
        atom_token = _normalize_atom_name(atom_name)
        if auth_token and auth_token not in {"?", "."}:
            if seq_token and residue_token and atom_token:
                atom_candidates = self.auth_residue_atom_to_canonical.get(
                    (auth_token, seq_token, residue_token, atom_token),
                    (),
                )
                if len(atom_candidates) == 1:
                    return atom_candidates[0]
            if seq_token and residue_token:
                residue_candidates = self.auth_residue_to_canonical.get(
                    (auth_token, seq_token, residue_token),
                    (),
                )
                if len(residue_candidates) == 1:
                    return residue_candidates[0]
            resolved_auth = sorted(self.resolve(auth_token))
            if len(resolved_auth) == 1:
                return resolved_auth[0]

        fallback = self.normalize(label_token or auth_token)
        return fallback or (label_token or auth_token)


def _identity_chain_aliases(chain_ids: Set[str]) -> ChainAliases:
    sanitized_ids = {
        str(chain_id or "").strip()
        for chain_id in chain_ids
        if str(chain_id or "").strip()
    }
    auth_to_canonical = {
        chain_id: (chain_id,)
        for chain_id in sorted(sanitized_ids)
    }
    canonical_to_auth = {
        chain_id: chain_id
        for chain_id in sorted(sanitized_ids)
    }
    return ChainAliases(
        label_to_auth={},
        auth_ids=set(sanitized_ids),
        canonical_ids=set(sanitized_ids),
        query_to_canonical={},
        auth_to_canonical=auth_to_canonical,
        canonical_to_auth=canonical_to_auth,
        auth_residue_to_canonical={},
        auth_residue_atom_to_canonical={},
        query_ids=set(sanitized_ids),
    )


PER_RESIDUE_COUNT_FIELDS = {
    "hydrophobic",
    "hbond",
    "polar_contact",
    "base_pairing",
    "salt_bridge",
    "halogen_bond",
    "metal_coordination",
    "pi_pi",
    "pi_cation",
    "aromatic_packing",
    "vdw",
    "clash",
    "other",
    "total",
}


def _external_chain_id(chain_id: str, aliases: ChainAliases) -> str:
    token = str(chain_id or "").strip()
    if not token:
        return token
    mapped = str(aliases.canonical_to_auth.get(token, "") or "").strip()
    return mapped or token


def _remap_prefixed_chain_key_for_external(key: object, aliases: ChainAliases) -> str:
    raw = str(key or "").strip()
    if not raw:
        return raw
    if ":" not in raw:
        return _external_chain_id(raw, aliases)
    prefix, suffix = raw.split(":", 1)
    mapped_prefix = _external_chain_id(prefix, aliases)
    if not mapped_prefix or mapped_prefix == prefix:
        return raw
    return f"{mapped_prefix}:{suffix}"


def _remap_residue_payload_for_external(residue: object, aliases: ChainAliases) -> object:
    if not isinstance(residue, dict):
        return residue
    chain = str(residue.get("chain") or "").strip()
    mapped_chain = _external_chain_id(chain, aliases)
    if not mapped_chain or mapped_chain == chain:
        return residue
    remapped = dict(residue)
    remapped["chain"] = mapped_chain
    return remapped


def _remap_contact_record_for_external(record: object, aliases: ChainAliases) -> object:
    if not isinstance(record, dict):
        return record
    residue_a = _remap_residue_payload_for_external(record.get("residueA"), aliases)
    residue_b = _remap_residue_payload_for_external(record.get("residueB"), aliases)
    remapped = dict(record)
    if residue_a is not record.get("residueA"):
        remapped["residueA"] = residue_a
    if residue_b is not record.get("residueB"):
        remapped["residueB"] = residue_b

    atom_key_a = _remap_prefixed_chain_key_for_external(remapped.get("atomKeyA"), aliases)
    atom_key_b = _remap_prefixed_chain_key_for_external(remapped.get("atomKeyB"), aliases)
    if atom_key_a:
        remapped["atomKeyA"] = atom_key_a
    if atom_key_b:
        remapped["atomKeyB"] = atom_key_b
    if atom_key_a and atom_key_b:
        remapped["pairKey"] = _build_unordered_pair_key(atom_key_a, atom_key_b)

    ring_key_a = _remap_prefixed_chain_key_for_external(remapped.get("ringKeyA"), aliases)
    ring_key_b = _remap_prefixed_chain_key_for_external(remapped.get("ringKeyB"), aliases)
    if ring_key_a:
        remapped["ringKeyA"] = ring_key_a
    if ring_key_b:
        remapped["ringKeyB"] = ring_key_b
    if ring_key_a and ring_key_b:
        remapped["ringPairKey"] = _build_unordered_pair_key(ring_key_a, ring_key_b)

    asserted = remapped.get("asserted")
    if isinstance(asserted, dict):
        asserted_remapped = dict(asserted)
        if ring_key_a:
            asserted_remapped["ringKeyA"] = ring_key_a
        if ring_key_b:
            asserted_remapped["ringKeyB"] = ring_key_b
        if ring_key_a and ring_key_b:
            asserted_remapped["ringPairKey"] = _build_unordered_pair_key(ring_key_a, ring_key_b)
        remapped["asserted"] = asserted_remapped

    if isinstance(record.get('semantics'), dict):
        semantics = copy.deepcopy(record['semantics'])
        replacements = {}
        for participant in semantics.get('participants', []):
            site = participant.get('site', {})
            residue = site.get('residue', {})
            old_chain = str(residue.get('chain') or '')
            new_chain = _external_chain_id(old_chain, aliases)
            if old_chain and new_chain and old_chain != new_chain:
                replacements[old_chain + ':'] = new_chain + ':'
        pattern = re.compile(r'(?<![A-Za-z0-9_])(?:' + '|'.join(re.escape(old) for old in sorted(replacements,key=len,reverse=True)) + ')') if replacements else None
        def remap_nested(value):
            if isinstance(value, dict):
                return {k: _external_chain_id(v, aliases) if k == 'chain' and isinstance(v, str) else remap_nested(v) for k,v in value.items()}
            if isinstance(value, list):
                return [remap_nested(v) for v in value]
            if isinstance(value, str) and pattern:
                return pattern.sub(lambda match: replacements[match.group(0)], value)
            return value
        semantics = remap_nested(semantics)
        semantics['identity'] = _semantic_identity(semantics)
        remapped['semantics'] = semantics
    return remapped


def _merge_external_per_residue_entry(target: dict, source: dict) -> dict:
    for field in PER_RESIDUE_COUNT_FIELDS:
        source_value = source.get(field)
        if isinstance(source_value, bool):
            continue
        if isinstance(source_value, (int, float)):
            target[field] = target.get(field, 0) + source_value
    for field in ("chain", "resName", "seq"):
        if not str(target.get(field) or "").strip():
            source_value = str(source.get(field) or "").strip()
            if source_value:
                target[field] = source_value
    return target


def _remap_per_residue_payload_for_external(
    per_residue: Dict[str, dict],
    aliases: ChainAliases,
) -> Dict[str, dict]:
    remapped: Dict[str, dict] = {}
    for residue_key, entry in per_residue.items():
        raw_key = str(residue_key or "").strip()
        mapped_key = _remap_prefixed_chain_key_for_external(raw_key, aliases)
        mapped_entry = dict(entry) if isinstance(entry, dict) else entry
        if isinstance(mapped_entry, dict):
            chain = str(mapped_entry.get("chain") or "").strip()
            mapped_chain = _external_chain_id(chain, aliases)
            if mapped_chain and mapped_chain != chain:
                mapped_entry["chain"] = mapped_chain
        if mapped_key in remapped and isinstance(mapped_entry, dict) and isinstance(remapped[mapped_key], dict):
            remapped[mapped_key] = _merge_external_per_residue_entry(remapped[mapped_key], mapped_entry)
            continue
        remapped[mapped_key] = mapped_entry
    return remapped


def _remap_focus_residue_for_external(focus_residue: str, aliases: ChainAliases) -> str:
    raw = str(focus_residue or "").strip()
    if not raw or ":" not in raw:
        return _external_chain_id(raw, aliases)
    chain, seq = raw.split(":", 1)
    mapped_chain = _external_chain_id(chain, aliases)
    if not mapped_chain or mapped_chain == chain:
        return raw
    return f"{mapped_chain}:{seq}"


def _remap_report_to_external_chains(
    *,
    chain_a: str,
    chain_b: str,
    contacts: Dict[str, List[dict]],
    per_residue: Dict[str, dict],
    buried_fraction: Optional[Dict[str, float]],
    meta: Dict[str, object],
    aliases: ChainAliases,
) -> Tuple[str, str, Dict[str, List[dict]], Dict[str, dict], Optional[Dict[str, float]], Dict[str, object]]:
    remapped_contacts = {
        bucket: [
            _remap_contact_record_for_external(record, aliases)
            for record in rows
        ]
        for bucket, rows in contacts.items()
    }
    remapped_per_residue = _remap_per_residue_payload_for_external(per_residue, aliases)
    remapped_buried_fraction = None
    if isinstance(buried_fraction, dict):
        remapped_buried_fraction = {}
        for chain_id, value in buried_fraction.items():
            remapped_chain = _external_chain_id(chain_id, aliases)
            remapped_buried_fraction[remapped_chain or chain_id] = value
    remapped_meta = dict(meta)
    focus_residue = str(remapped_meta.get("focusResidue") or "").strip()
    if focus_residue:
        remapped_meta["focusResidue"] = _remap_focus_residue_for_external(focus_residue, aliases)
    return (
        _external_chain_id(chain_a, aliases),
        _external_chain_id(chain_b, aliases),
        remapped_contacts,
        remapped_per_residue,
        remapped_buried_fraction,
        remapped_meta,
    )


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


def _build_aromatic_candidate_graph(
    atoms: List[AtomRecord],
) -> Tuple[List[AtomRecord], List[List[int]], List[List[int]]]:
    candidates = [atom for atom in atoms if (atom.element or "").upper() in AROMATIC_CANDIDATE_ELEMENTS]
    if len(candidates) < 5 or len(candidates) > 80:
        return candidates, [], []

    min_d_sq = AROMATIC_BOND_MIN * AROMATIC_BOND_MIN
    max_d_sq = AROMATIC_BOND_MAX * AROMATIC_BOND_MAX
    adjacency: List[List[int]] = [[] for _ in candidates]
    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            d_sq = _distance_sq(candidates[i], candidates[j])
            if min_d_sq <= d_sq <= max_d_sq:
                adjacency[i].append(j)
                adjacency[j].append(i)

    visited = [False] * len(candidates)
    components: List[List[int]] = []
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
        components.append(component)
    return candidates, adjacency, components


def _component_edge_count(component_indices: List[int], adjacency: List[List[int]]) -> int:
    comp_set = set(component_indices)
    edge_count = 0
    for idx in component_indices:
        edge_count += sum(1 for nbr in adjacency[idx] if nbr in comp_set)
    return edge_count // 2


def _is_aromatic_component_indices(
    *,
    candidates: List[AtomRecord],
    adjacency: List[List[int]],
    component_indices: List[int],
) -> bool:
    size = len(component_indices)
    if size < 5 or size > 28:
        return False
    edge_count = _component_edge_count(component_indices, adjacency)
    if edge_count < size:
        return False
    avg_degree = (2.0 * edge_count) / size
    if avg_degree < 1.45:
        return False
    component_atoms = [candidates[idx] for idx in component_indices]
    return _is_planar_component(component_atoms, AROMATIC_PLANAR_RMSD_MAX)


def _two_core_component_indices(
    component_indices: List[int],
    adjacency: List[List[int]],
) -> List[int]:
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


def _find_planar_local_ring_component_index_sets(
    *,
    candidates: List[AtomRecord],
    adjacency: List[List[int]],
    component_indices: List[int],
) -> List[List[int]]:
    # Large ligands can contain aromatic sub-rings embedded in a much bigger
    # non-planar graph (e.g. taxol-like scaffolds). Detect local 5-7 atom
    # cycles instead of requiring the whole component to be aromatic-like.
    if len(component_indices) < 5:
        return []
    comp_set = set(component_indices)
    neighbors: Dict[int, List[int]] = {
        idx: [nbr for nbr in adjacency[idx] if nbr in comp_set]
        for idx in comp_set
    }
    min_cycle = 5
    max_cycle = 7
    ring_keys: Dict[Tuple[int, ...], None] = {}

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
                            ring_key = tuple(sorted(path))
                            ring_keys[ring_key] = None
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
    ordered_keys = sorted(
        ring_keys.keys(),
        key=lambda key: (-len(key), key),
    )
    return [list(key) for key in ordered_keys]


def _find_planar_local_ring_component_indices(
    *,
    candidates: List[AtomRecord],
    adjacency: List[List[int]],
    component_indices: List[int],
) -> List[int]:
    ring_sets = _find_planar_local_ring_component_index_sets(
        candidates=candidates,
        adjacency=adjacency,
        component_indices=component_indices,
    )
    if not ring_sets:
        return []
    return ring_sets[0]


def _infer_nonpolymer_ring_candidate_index_sets(
    atoms: List[AtomRecord],
) -> Tuple[List[AtomRecord], List[List[int]]]:
    candidates, adjacency, components = _build_aromatic_candidate_graph(atoms)
    if not candidates or not adjacency or not components:
        return candidates, []

    ring_map: Dict[Tuple[int, ...], Tuple[int, int]] = {}

    def _consider(indices: List[int]) -> None:
        if not indices:
            return
        unique_indices = tuple(sorted(set(indices)))
        if len(unique_indices) < 5:
            return
        edge_count = _component_edge_count(list(unique_indices), adjacency)
        score = (len(unique_indices), edge_count)
        current = ring_map.get(unique_indices)
        if current is None or score > current:
            ring_map[unique_indices] = score

    for component in components:
        if _is_aromatic_component_indices(
            candidates=candidates,
            adjacency=adjacency,
            component_indices=component,
        ):
            _consider(component)

        core = _two_core_component_indices(component, adjacency)
        if core and core != component and _is_aromatic_component_indices(
            candidates=candidates,
            adjacency=adjacency,
            component_indices=core,
        ):
            _consider(core)

        local_rings = _find_planar_local_ring_component_index_sets(
            candidates=candidates,
            adjacency=adjacency,
            component_indices=component,
        )
        for local_ring in local_rings:
            _consider(local_ring)

        if core and core != component:
            local_rings_core = _find_planar_local_ring_component_index_sets(
                candidates=candidates,
                adjacency=adjacency,
                component_indices=core,
            )
            for local_ring_core in local_rings_core:
                _consider(local_ring_core)

    ordered = sorted(
        ring_map.items(),
        key=lambda item: (-item[1][0], -item[1][1], item[0]),
    )
    ring_sets = [list(indices) for indices, _ in ordered]
    return candidates, ring_sets


def _infer_nonpolymer_ring_candidates(
    atoms: List[AtomRecord],
) -> Tuple[List[AtomRecord], List[int]]:
    candidates, ring_sets = _infer_nonpolymer_ring_candidate_index_sets(atoms)
    if not ring_sets:
        return candidates, []
    return candidates, ring_sets[0]


def _is_likely_aromatic_nonpolymer(atoms: List[AtomRecord]) -> bool:
    _, ring_indices = _infer_nonpolymer_ring_candidates(atoms)
    return bool(ring_indices)


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
    chains = sorted(aliases.auth_ids)
    return chains, aliases


def _structure_digest(structure_text: str) -> str:
    return hashlib.sha256((structure_text or "").encode("utf-8")).hexdigest()[:16]


def _parse_structure_cached(structure_text: str, structure_format: str) -> Tuple[List[AtomRecord], ChainAliases]:
    fmt = str(structure_format or "mmcif").strip().lower()
    digest = _structure_digest(structure_text)
    key = f"{fmt}:{digest}"
    cached = PARSED_STRUCTURE_CACHE.get(key)
    if isinstance(cached, tuple) and len(cached) == 2:
        atoms = cached[0]
        aliases = cached[1]
        if isinstance(atoms, list) and isinstance(aliases, ChainAliases):
            return atoms, aliases
    if fmt == "pdb":
        parsed = parse_pdb_atoms(structure_text)
    else:
        parsed = parse_mmcif_atoms(structure_text)
    PARSED_STRUCTURE_CACHE.set(key, parsed)
    return parsed


def _residue_sort_key(res_seq: str) -> Tuple[int, str]:
    token = str(res_seq or "").strip()
    if not token:
        return (10**9, "")
    numeric = _parse_res_seq_index(token)
    if numeric is None:
        return (10**9, token)
    suffix = token[len(str(numeric)) :] if token.startswith(str(numeric)) else token
    return (numeric, suffix)


def _build_arpeggio_selection_for_chains(
    atoms: List[AtomRecord],
    chain_ids: Set[str],
) -> List[str]:
    residues = sorted(
        {
            (
                str(atom.chain_auth or atom.chain_id).strip(),
                str(atom.res_seq or "").strip(),
            )
            for atom in atoms
            if atom.chain_id in chain_ids
            and str(atom.chain_auth or atom.chain_id).strip()
            and atom.res_seq
        },
        key=lambda row: (row[0], _residue_sort_key(row[1])),
    )
    selection = [f"/{chain_id}/{res_seq}/" for chain_id, res_seq in residues if res_seq]
    if selection:
        return selection
    auth_chains = sorted(
        {
            str(atom.chain_auth or atom.chain_id).strip()
            for atom in atoms
            if atom.chain_id in chain_ids and str(atom.chain_auth or atom.chain_id).strip()
        }
    )
    return [f"/{chain_id}/" for chain_id in auth_chains]


def _parse_focus_residue_key(residue_key: object) -> Optional[Tuple[str, str]]:
    token = str(residue_key or "").strip()
    if not token:
        return None
    if ":" not in token:
        return None
    chain_id, res_seq = token.split(":", 1)
    chain_token = chain_id.strip()
    seq_token = res_seq.strip()
    if not chain_token or not seq_token:
        return None
    return chain_token, seq_token


def _spatial_cell_key(x: float, y: float, z: float, cell_size: float) -> Tuple[int, int, int]:
    return (
        int(math.floor(x / cell_size)),
        int(math.floor(y / cell_size)),
        int(math.floor(z / cell_size)),
    )


def _source_residues_within_cutoff(
    source_atoms: List[AtomRecord],
    target_atoms: List[AtomRecord],
    cutoff: float,
) -> Set[Tuple[str, str]]:
    if not source_atoms or not target_atoms:
        return set()
    cutoff_sq = float(cutoff) * float(cutoff)
    cell_size = max(1.0, float(cutoff))
    target_grid: Dict[Tuple[int, int, int], List[AtomRecord]] = {}
    for atom in target_atoms:
        key = _spatial_cell_key(atom.x, atom.y, atom.z, cell_size)
        target_grid.setdefault(key, []).append(atom)

    nearby_residues: Set[Tuple[str, str]] = set()
    for atom in source_atoms:
        selection_chain = str(atom.chain_auth or atom.chain_id).strip()
        res_seq = str(atom.res_seq or "").strip()
        if not selection_chain or not res_seq:
            continue
        base_key = _spatial_cell_key(atom.x, atom.y, atom.z, cell_size)
        found = False
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    bucket = (
                        base_key[0] + dx,
                        base_key[1] + dy,
                        base_key[2] + dz,
                    )
                    for target in target_grid.get(bucket, ()):
                        if _distance_sq(atom, target) <= cutoff_sq:
                            nearby_residues.add((selection_chain, res_seq))
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break
    return nearby_residues


def _build_arpeggio_selection(
    atoms: List[AtomRecord],
    chain_ids_a: Set[str],
    chain_ids_b: Set[str],
    focus_residue_candidates: Optional[List[Tuple[str, str]]] = None,
) -> Tuple[List[str], List[Tuple[str, str]]]:
    if focus_residue_candidates:
        selection: List[str] = []
        applied_focuses: List[Tuple[str, str]] = []
        seen_focuses: Set[Tuple[str, str]] = set()
        seen_selection_tokens: Set[Tuple[str, str]] = set()
        for focus_chain, focus_seq in focus_residue_candidates:
            chain_token = str(focus_chain or "").strip()
            seq_token = str(focus_seq or "").strip()
            key = (chain_token, seq_token)
            if not chain_token or not seq_token or key in seen_focuses:
                continue
            seen_focuses.add(key)
            matching_atoms = [
                atom
                for atom in atoms
                if atom.chain_id == chain_token and str(atom.res_seq or "").strip() == seq_token
            ]
            if not matching_atoms:
                continue
            for atom in matching_atoms:
                selection_chain = str(atom.chain_auth or atom.chain_id).strip()
                selection_key = (selection_chain, seq_token)
                if not selection_chain or selection_key in seen_selection_tokens:
                    continue
                seen_selection_tokens.add(selection_key)
                selection.append(f"/{selection_chain}/{seq_token}/")
            applied_focuses.append(key)
        if selection:
            return selection, applied_focuses

    if chain_ids_a != chain_ids_b:
        atoms_a = [atom for atom in atoms if atom.chain_id in chain_ids_a]
        atoms_b = [atom for atom in atoms if atom.chain_id in chain_ids_b]
        interface_residues = _source_residues_within_cutoff(
            atoms_a,
            atoms_b,
            ARPEGGIO_INTERFACE_SELECTION_CUTOFF,
        )
        if interface_residues:
            residues = sorted(
                interface_residues,
                key=lambda row: (row[0], _residue_sort_key(row[1])),
            )
            selection = [f"/{chain_id}/{res_seq}/" for chain_id, res_seq in residues if res_seq]
            if selection:
                return selection, []

    return _build_arpeggio_selection_for_chains(atoms, chain_ids_a), []


def _normalize_arpeggio_contact_terms(value: object) -> List[str]:
    if isinstance(value, str):
        text = value.strip().upper()
        return [text] if text else []
    if not isinstance(value, list):
        return []
    output: List[str] = []
    seen: set[str] = set()
    for row in value:
        token = str(row or "").strip().upper()
        if not token or token in seen:
            continue
        seen.add(token)
        output.append(token)
    return output


def _coerce_float(value: object) -> Optional[float]:
    try:
        numeric = float(value)
    except Exception:
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def _normalize_contact_category(value: object) -> str:
    token = str(value or "").strip().lower().replace("-", "_")
    if not token:
        return ""
    if token in {"hbond", "hydrogen_bond", "hydrogen_bonds"}:
        return "hbond"
    if token in {"salt_bridge", "salt_bridges", "electrostatic", "ionic"}:
        return "salt_bridge"
    if token in {"halogen_bond", "halogen_bonds", "halogen", "xbond", "x_bond"}:
        return "halogen_bond"
    if token in {"hydrophobic"}:
        return "hydrophobic"
    if token in {"metal_coordination", "metal"}:
        return "metal_coordination"
    if token in {"aromatic_packing", "pi_packing", "stacking_contact", "stacking_adjacent"}:
        return "aromatic_packing"
    if token in {"aromatic_proximal", "ring_neighbor", "ring_neighbour"}:
        return "aromatic_proximal"
    if token in {"pi_pi", "pi"}:
        return "pi_pi"
    if token in {"pi_cation", "cation_pi", "cationpi"}:
        return "pi_cation"
    if token in {"packing_contact", "packing", "vdw_close"}:
        return "packing_contact"
    if token in {"proximal", "neighborhood"}:
        return "proximal"
    if token in {"polar_contact", "polar_contacts", "polar"}:
        return "polar_contact"
    if token in {"polar_proximal", "polar_neighbor", "polar_neighbour"}:
        return "polar_proximal"
    if token in {"base_pairing", "base_pair", "base_pairs"}:
        return "base_pairing"
    if token in {"covalent_bond", "covalent"}:
        return "covalent_bond"
    if token in {"vdw"}:
        return "vdw"
    if token in {"clash", "vdw_clash"}:
        return "clash"
    if token in {"other"}:
        return "other"
    return ""


def _normalize_atom_name(atom_name: object) -> str:
    return str(atom_name or "").strip().upper().replace("*", "'")


def _is_protein_backbone_atom_name(atom_name: object) -> bool:
    name = _normalize_atom_name(atom_name)
    return name in {"N", "CA", "C", "O", "OXT", "OT1", "OT2"}


def _is_likely_carbonyl_carbon_for_hydrophobic(
    residue_name: object,
    atom_name: object,
) -> bool:
    res_name = str(residue_name or "").strip().upper()
    atom = _normalize_atom_name(atom_name)
    if not res_name or not atom:
        return False
    # Protein backbone carbonyl carbon.
    if res_name in STANDARD_AMINO_RESIDUES and atom == "C":
        return True
    sidechain_set = PROTEIN_SIDECHAIN_CARBONYL_CARBON_BY_RESIDUE.get(res_name)
    if sidechain_set and atom in sidechain_set:
        return True
    return False


def _is_charged_group_associated_carbon_for_hydrophobic(
    residue_name: object,
    atom_name: object,
) -> bool:
    res_name = str(residue_name or "").strip().upper()
    atom = _normalize_atom_name(atom_name)
    if not res_name or not atom:
        return False
    sidechain_set = PROTEIN_CHARGED_GROUP_ASSOCIATED_CARBON_BY_RESIDUE.get(res_name)
    if sidechain_set and atom in sidechain_set:
        return True
    return False


def _is_hydrophobic_contact_atom_candidate(
    element: str,
    atom_name: str,
    residue_name: str,
) -> bool:
    res_name = str(residue_name or "").strip().upper()
    if not res_name or res_name in WATER_RESIDUES:
        return False
    element_token = str(element or "").strip().upper()
    if not element_token:
        element_token = guess_element(str(atom_name or "")).upper()
    if element_token not in {"C", "S", "SE", "F", "CL", "BR", "I"}:
        return False
    # Carbonyl/carboxylamide carbons are polar-adjacent electrophilic centers and
    # should not be promoted as hydrophobic atom sites.
    if element_token == "C" and _is_likely_carbonyl_carbon_for_hydrophobic(res_name, atom_name):
        return False
    # Charged functional-group core carbons (e.g., ARG CZ) are better
    # interpreted as packing near charged sites than true hydrophobic contacts.
    if element_token == "C" and _is_charged_group_associated_carbon_for_hydrophobic(res_name, atom_name):
        return False
    # Treat protein backbone heavy atoms as structural scaffold, not hydrophobic contacts.
    if res_name in STANDARD_AMINO_RESIDUES and _is_protein_backbone_atom_name(atom_name):
        return False
    return True


def _is_polar_or_charged_sidechain_residue(residue_name: object) -> bool:
    res_name = str(residue_name or "").strip().upper()
    if not res_name or res_name not in STANDARD_AMINO_RESIDUES:
        return False
    return res_name in POLAR_OR_CHARGED_SIDECHAIN_RESIDUES


def _is_packing_contact_atom_candidate(
    element: str,
    atom_name: str,
    residue_name: str,
) -> bool:
    res_name = str(residue_name or "").strip().upper()
    if not res_name or res_name in WATER_RESIDUES:
        return False
    element_token = str(element or "").strip().upper()
    if not element_token:
        element_token = guess_element(str(atom_name or "")).upper()
    if element_token not in PACKING_ELIGIBLE_ELEMENTS:
        return False
    if element_token == "C" and _is_likely_carbonyl_carbon_for_hydrophobic(res_name, atom_name):
        return False
    return True


def _is_phosphate_oxygen_atom_name(atom_name: object) -> bool:
    name = _normalize_atom_name(atom_name)
    if not name:
        return False
    if name in {"OP1", "OP2", "OP3", "O1P", "O2P", "O3P"}:
        return True
    # Covers nucleotide/ligand phosphate oxygens such as O1A/O2A/O3A, O1B/O2B/O3B, O1G/O2G/O3G.
    if re.match(r"^O[123][A-Z]$", name):
        return True
    return False


def _is_sulfate_oxygen_site(
    res_name: object,
    atom_name: object,
    element: object = "",
) -> bool:
    residue = str(res_name or "").strip().upper().replace("+", "").replace("-", "")
    if residue not in SULFATE_LIKE_RESIDUE_NAMES:
        return False
    atom = _normalize_atom_name(atom_name)
    if not atom or not atom.startswith("O"):
        return False
    atom_element = str(element or "").strip().upper()
    if not atom_element:
        atom_element = guess_element(atom).upper()
    return atom_element == "O"


def _nucleic_base_family(res_name: object) -> str:
    residue = str(res_name or "").strip().upper()
    if residue in NUCLEIC_A_FAMILY:
        return "A"
    if residue in NUCLEIC_G_FAMILY:
        return "G"
    if residue in NUCLEIC_C_FAMILY:
        return "C"
    if residue in NUCLEIC_U_FAMILY:
        return "U"
    if residue in NUCLEIC_T_FAMILY:
        return "T"
    if residue in NUCLEIC_I_FAMILY:
        return "I"
    return ""


def _is_nucleobase_donor_nitrogen(res_name: object, atom_name: object) -> bool:
    atom = _normalize_atom_name(atom_name)
    family = _nucleic_base_family(res_name)
    if family == "A":
        return atom == "N6"
    if family == "G":
        return atom in {"N1", "N2"}
    if family == "C":
        return atom == "N4"
    if family in {"U", "T"}:
        return atom == "N3"
    if family == "I":
        return atom == "N1"
    return False


def _is_nucleobase_acceptor_nitrogen(res_name: object, atom_name: object) -> bool:
    atom = _normalize_atom_name(atom_name)
    family = _nucleic_base_family(res_name)
    if family == "A":
        return atom in {"N1", "N3", "N7"}
    if family == "G":
        return atom in {"N3", "N7"}
    if family == "C":
        return atom == "N3"
    if family == "I":
        return atom in {"N3", "N7"}
    return False


def _is_nucleobase_acceptor_oxygen(atom_name: object) -> bool:
    atom = _normalize_atom_name(atom_name)
    return atom in {"O2", "O4", "O6"}


def _is_nucleobase_atom(res_name: object, atom_name: object) -> bool:
    family = _nucleic_base_family(res_name)
    if not family:
        return False
    atom = _normalize_atom_name(atom_name)
    if not atom:
        return False
    allowed = NUCLEOBASE_ATOMS_BY_FAMILY.get(family)
    if not allowed:
        return False
    return atom in allowed


def _is_nucleobase_pairing_edge_atom(res_name: object, atom_name: object) -> bool:
    family = _nucleic_base_family(res_name)
    if not family:
        return False
    atom = _normalize_atom_name(atom_name)
    if not atom:
        return False
    allowed = NUCLEOBASE_PAIRING_EDGE_ATOMS_BY_FAMILY.get(family)
    if not allowed:
        return False
    return atom in allowed


def _is_nucleobase_glycosidic_atom(res_name: object, atom_name: object) -> bool:
    family = _nucleic_base_family(res_name)
    if not family:
        return False
    atom = _normalize_atom_name(atom_name)
    if not atom:
        return False
    glyco = str(NUCLEOBASE_GLYCOSIDIC_ATOM_BY_FAMILY.get(family) or "").strip().upper()
    if not glyco:
        return False
    return atom == glyco


def _is_nucleotide_backbone_atom_name(atom_name: object) -> bool:
    atom = _normalize_atom_name(atom_name)
    if not atom:
        return False
    if atom in NUCLEOTIDE_BACKBONE_ATOMS:
        return True
    if atom.startswith("OP"):
        return True
    if re.match(r"^O[123]P$", atom):
        return True
    return False


def _is_nucleic_backbone_oxygen_or_phosphate_site(
    *,
    res_name: object,
    atom_name: object,
    element: object,
) -> bool:
    family = _nucleic_base_family(res_name)
    if not family:
        return False
    atom = _normalize_atom_name(atom_name)
    if not atom:
        return False
    atom_element = str(element or "").strip().upper()
    if not atom_element:
        atom_element = guess_element(atom).upper()
    if atom_element == "P":
        return atom == "P"
    if atom_element != "O":
        return False
    if _is_phosphate_oxygen_atom_name(atom):
        return True
    if atom in NUCLEIC_SUGAR_ACCEPTOR_OXYGENS:
        return True
    return False


def _is_nucleic_phosphate_sugar_oxygen_pair(
    *,
    res_name_a: object,
    atom_name_a: object,
    element_a: object,
    res_name_b: object,
    atom_name_b: object,
    element_b: object,
) -> bool:
    family_a = _nucleic_base_family(res_name_a)
    family_b = _nucleic_base_family(res_name_b)
    if not family_a or not family_b:
        return False
    atom_a = _normalize_atom_name(atom_name_a)
    atom_b = _normalize_atom_name(atom_name_b)
    if not atom_a or not atom_b:
        return False
    element_token_a = str(element_a or "").strip().upper() or guess_element(atom_a).upper()
    element_token_b = str(element_b or "").strip().upper() or guess_element(atom_b).upper()
    if element_token_a != "O" or element_token_b != "O":
        return False
    phosphate_like_a = _is_phosphate_oxygen_atom_name(atom_a)
    phosphate_like_b = _is_phosphate_oxygen_atom_name(atom_b)
    sugar_like_a = atom_a in NUCLEIC_SUGAR_ACCEPTOR_OXYGENS
    sugar_like_b = atom_b in NUCLEIC_SUGAR_ACCEPTOR_OXYGENS
    return bool((phosphate_like_a and sugar_like_b) or (phosphate_like_b and sugar_like_a))


def _is_nucleic_backbone_oxygen_neighborhood_pair(
    *,
    res_name_a: object,
    atom_name_a: object,
    element_a: object,
    res_name_b: object,
    atom_name_b: object,
    element_b: object,
) -> bool:
    family_a = _nucleic_base_family(res_name_a)
    family_b = _nucleic_base_family(res_name_b)
    if not family_a or not family_b:
        return False
    atom_a = _normalize_atom_name(atom_name_a)
    atom_b = _normalize_atom_name(atom_name_b)
    if not atom_a or not atom_b:
        return False
    element_token_a = str(element_a or "").strip().upper() or guess_element(atom_a).upper()
    element_token_b = str(element_b or "").strip().upper() or guess_element(atom_b).upper()
    if element_token_a != "O" or element_token_b != "O":
        return False

    phosphate_like_a = _is_phosphate_oxygen_atom_name(atom_a)
    phosphate_like_b = _is_phosphate_oxygen_atom_name(atom_b)
    sugar_like_a = atom_a in NUCLEIC_SUGAR_ACCEPTOR_OXYGENS
    sugar_like_b = atom_b in NUCLEIC_SUGAR_ACCEPTOR_OXYGENS
    terminal_hydroxyl_like_a = atom_a in {"O3'", "O3*", "O5'", "O5*"}
    terminal_hydroxyl_like_b = atom_b in {"O3'", "O3*", "O5'", "O5*"}

    if phosphate_like_a and phosphate_like_b:
        return True
    if (phosphate_like_a and sugar_like_b) or (phosphate_like_b and sugar_like_a):
        return True
    if (phosphate_like_a and terminal_hydroxyl_like_b) or (phosphate_like_b and terminal_hydroxyl_like_a):
        return True
    return False


def _is_likely_adjacent_nucleotide_linkage_contact(
    *,
    residue_a: dict,
    residue_b: dict,
    base_family_a: str,
    base_family_b: str,
    atom_name_a: str,
    atom_name_b: str,
) -> bool:
    if not base_family_a or not base_family_b:
        return False
    if not isinstance(residue_a, dict) or not isinstance(residue_b, dict):
        return False
    chain_a = str(residue_a.get("chain") or "").strip()
    chain_b = str(residue_b.get("chain") or "").strip()
    if not chain_a or not chain_b or chain_a != chain_b:
        return False
    seq_a = _parse_res_seq_index(str(residue_a.get("seq") or "").strip())
    seq_b = _parse_res_seq_index(str(residue_b.get("seq") or "").strip())
    if seq_a is None or seq_b is None:
        return False
    if abs(seq_a - seq_b) > 1:
        return False

    atom_a = _normalize_atom_name(atom_name_a)
    atom_b = _normalize_atom_name(atom_name_b)
    if not atom_a or not atom_b:
        return False
    if not (_is_nucleotide_backbone_atom_name(atom_a) and _is_nucleotide_backbone_atom_name(atom_b)):
        return False

    a_phosphate_like = atom_a in NUCLEOTIDE_LINKAGE_PHOSPHATE_ATOMS
    b_phosphate_like = atom_b in NUCLEOTIDE_LINKAGE_PHOSPHATE_ATOMS
    a_o3 = atom_a in NUCLEOTIDE_LINKAGE_O3_ATOMS
    b_o3 = atom_b in NUCLEOTIDE_LINKAGE_O3_ATOMS
    if (a_phosphate_like and b_o3) or (b_phosphate_like and a_o3):
        return True
    # Extra conservative guard: hide short O3'/O5' neighboring linkage contacts.
    if (a_o3 and atom_b in {"O5'", "O5*"}) or (b_o3 and atom_a in {"O5'", "O5*"}):
        return True
    return False


def _is_sequence_adjacent_nucleotide_pair(
    *,
    residue_a: dict,
    residue_b: dict,
    base_family_a: str,
    base_family_b: str,
) -> bool:
    if not base_family_a or not base_family_b:
        return False
    if not isinstance(residue_a, dict) or not isinstance(residue_b, dict):
        return False
    chain_a = str(residue_a.get("chain") or "").strip()
    chain_b = str(residue_b.get("chain") or "").strip()
    if not chain_a or not chain_b or chain_a != chain_b:
        return False
    seq_a = _parse_res_seq_index(str(residue_a.get("seq") or "").strip())
    seq_b = _parse_res_seq_index(str(residue_b.get("seq") or "").strip())
    if seq_a is None or seq_b is None:
        return False
    return abs(seq_a - seq_b) == 1


def _residue_ring_atom_names(res_name: str) -> Set[str]:
    residue = str(res_name or "").strip().upper()
    if not residue:
        return set()
    aromatic_set = AROMATIC_RING_ATOMS_BY_RESIDUE.get(residue)
    if aromatic_set:
        return set(aromatic_set)
    family = _nucleic_base_family(residue)
    if not family:
        return set()
    return set(NUCLEOBASE_RING_ATOMS_BY_FAMILY.get(family) or ())


def _residue_ring_point_sets(residue, residue_atoms_index, *, nucleobase_only=False):
    cache = getattr(residue_atoms_index, "ring_points", None)
    key = (residue.get("chain"), residue.get("seq"), residue.get("resName"), nucleobase_only)
    if cache is not None and key in cache:
        return cache[key]
    result = _residue_ring_point_sets_uncached(residue, residue_atoms_index, nucleobase_only=nucleobase_only)
    if cache is not None:
        cache[key] = result
    return result


def _residue_ring_point_sets_uncached(
    residue: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
    *,
    nucleobase_only: bool = False,
) -> List[List[Tuple[float, float, float]]]:
    if not isinstance(residue, dict) or not residue_atoms_index:
        return []
    chain = str(residue.get("chain") or "").strip()
    seq = str(residue.get("seq") or "").strip()
    res_name = str(residue.get("resName") or "").strip().upper()
    if not chain or not seq or not res_name:
        return []
    atoms = residue_atoms_index.get((chain, seq), [])
    if not atoms:
        return []

    point_sets: List[List[Tuple[float, float, float]]] = []
    ring_names = _residue_ring_atom_names(res_name)
    if nucleobase_only and not _nucleic_base_family(res_name):
        return []
    if ring_names:
        points: List[Tuple[float, float, float]] = []
        for atom in atoms:
            atom_name = _normalize_atom_name(atom.atom_name)
            if atom_name not in ring_names:
                continue
            element = str(atom.element or "").strip().upper()
            if element == "H":
                continue
            points.append((atom.x, atom.y, atom.z))
        if points:
            point_sets.append(points)
    else:
        if nucleobase_only:
            return []
        is_nonpolymer_residue = bool(
            res_name
            and res_name not in POLYMER_RESIDUES
            and res_name not in WATER_RESIDUES
            and res_name not in METAL_ELEMENTS
        )
        if not is_nonpolymer_residue:
            return []
        ring_candidates, ring_index_sets = _infer_nonpolymer_ring_candidate_index_sets(atoms)
        if not ring_index_sets:
            return []
        for ring_indices in ring_index_sets:
            points: List[Tuple[float, float, float]] = []
            for idx in ring_indices:
                if idx < 0 or idx >= len(ring_candidates):
                    continue
                atom = ring_candidates[idx]
                element = str(atom.element or "").strip().upper()
                if element == "H":
                    continue
                points.append((atom.x, atom.y, atom.z))
            if len(points) >= 4:
                point_sets.append(points)
    return point_sets


def _build_atom_key_from_record(atom: AtomRecord) -> str:
    if not isinstance(atom, AtomRecord):
        return ""
    chain = str(atom.chain_id or "").strip()
    seq = str(atom.res_seq or "").strip()
    if not chain or not seq:
        return ""
    res_name = str(atom.res_name or "").strip().upper() or "UNK"
    atom_name = _normalize_atom_name(atom.atom_name) or "?"
    return f"{chain}:{seq}:{res_name}:{atom_name}"


def _ring_label_for_residue_name(res_name: str, descriptor_index: int = 0) -> str:
    residue = str(res_name or "").strip().upper()
    if residue == "HIS":
        return "imidazole"
    if residue in {"PHE", "TYR"}:
        return "phenyl"
    if residue == "TRP":
        return "indole"
    nucleic_family = _nucleic_base_family(residue)
    if nucleic_family in {"A", "G", "I"}:
        return "purine"
    if nucleic_family in {"C", "U", "T"}:
        return "pyrimidine"
    return f"ring_{max(0, descriptor_index) + 1}"


def _build_ring_descriptor_from_atoms(
    ring_atoms: List[AtomRecord],
    *,
    residue_name: str,
    descriptor_index: int = 0,
) -> Optional[dict]:
    if not ring_atoms:
        return None
    points: List[Tuple[float, float, float]] = []
    atom_names: List[str] = []
    atom_keys: List[str] = []
    for atom in ring_atoms:
        if not isinstance(atom, AtomRecord):
            continue
        element = str(atom.element or "").strip().upper()
        if element == "H":
            continue
        points.append((atom.x, atom.y, atom.z))
        atom_name = _normalize_atom_name(atom.atom_name)
        if atom_name:
            atom_names.append(atom_name)
        atom_key = _build_atom_key_from_record(atom)
        if atom_key:
            atom_keys.append(atom_key)
    if len(points) < 3:
        return None
    atom_names = sorted(set(atom_names))
    atom_keys = sorted(set(atom_keys))
    hash_source_tokens = atom_keys if atom_keys else atom_names
    if not hash_source_tokens:
        hash_source_tokens = [
            ",".join(f"{coord:.3f}" for coord in point)
            for point in points
        ]
    digest = hashlib.sha1("|".join(hash_source_tokens).encode("utf-8")).hexdigest()[:12]
    label = _ring_label_for_residue_name(residue_name, descriptor_index=descriptor_index)
    centroid = (
        sum(point[0] for point in points) / len(points),
        sum(point[1] for point in points) / len(points),
        sum(point[2] for point in points) / len(points),
    )
    return {
        "points": points,
        "atom_names": atom_names,
        "atom_keys": atom_keys,
        "hash": digest,
        "label": label,
        "centroid": centroid,
    }


def _residue_ring_descriptors(residue, residue_atoms_index, *, nucleobase_only=False):
    cache = getattr(residue_atoms_index, "ring_descriptors", None)
    key = (residue.get("chain"), residue.get("seq"), residue.get("resName"), nucleobase_only)
    if cache is not None and key in cache:
        return cache[key]
    result = _residue_ring_descriptors_uncached(residue, residue_atoms_index, nucleobase_only=nucleobase_only)
    if cache is not None:
        cache[key] = result
    return result


def _residue_ring_descriptors_uncached(
    residue: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
    *,
    nucleobase_only: bool = False,
) -> List[dict]:
    if not isinstance(residue, dict) or not residue_atoms_index:
        return []
    chain = str(residue.get("chain") or "").strip()
    seq = str(residue.get("seq") or "").strip()
    res_name = str(residue.get("resName") or "").strip().upper()
    if not chain or not seq or not res_name:
        return []
    atoms = residue_atoms_index.get((chain, seq), [])
    if not atoms:
        return []

    descriptors: List[dict] = []
    ring_names = _residue_ring_atom_names(res_name)
    if nucleobase_only and not _nucleic_base_family(res_name):
        return []
    if ring_names:
        ring_sets = [ring_names]
        if res_name == "TRP":
            ring_sets = [{"CG", "CD1", "NE1", "CE2", "CD2"}, {"CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"}]
        elif _nucleic_base_family(res_name) in {"A", "G", "I"}:
            ring_sets = [{"N9", "C8", "N7", "C5", "C4"}, {"N1", "C2", "N3", "C4", "C5", "C6"}]
        for index, names in enumerate(ring_sets):
            ring_atoms = [atom for atom in atoms if _normalize_atom_name(atom.atom_name) in names and atom.element != "H"]
            descriptor = _build_ring_descriptor_from_atoms(ring_atoms, residue_name=res_name, descriptor_index=index)
            if descriptor:
                descriptors.append(descriptor)
    else:
        if nucleobase_only:
            return []
        is_nonpolymer_residue = bool(
            res_name
            and res_name not in POLYMER_RESIDUES
            and res_name not in WATER_RESIDUES
            and res_name not in METAL_ELEMENTS
        )
        if not is_nonpolymer_residue:
            return []
        ring_candidates, ring_index_sets = _infer_nonpolymer_ring_candidate_index_sets(atoms)
        if not ring_index_sets:
            return []
        for ring_index, ring_indices in enumerate(ring_index_sets):
            ring_atoms: List[AtomRecord] = []
            for idx in ring_indices:
                if idx < 0 or idx >= len(ring_candidates):
                    continue
                atom = ring_candidates[idx]
                element = str(atom.element or "").strip().upper()
                if element == "H":
                    continue
                ring_atoms.append(atom)
            descriptor = _build_ring_descriptor_from_atoms(
                ring_atoms,
                residue_name=res_name,
                descriptor_index=ring_index,
            )
            if descriptor and len(descriptor.get("points") or []) >= 4:
                descriptors.append(descriptor)
    descriptors.sort(
        key=lambda item: (
            str(item.get("hash") or ""),
            str(item.get("label") or ""),
        )
    )
    return descriptors


def _ring_descriptor_contains_atom_name(descriptor: dict, atom_name: str) -> bool:
    if not isinstance(descriptor, dict):
        return False
    token = _normalize_atom_name(atom_name)
    if not token:
        return False
    atom_names = descriptor.get("atom_names")
    if not isinstance(atom_names, list):
        return False
    return token in atom_names


def _select_ring_descriptor_pair_by_contact(
    residue_a: dict,
    residue_b: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
    *,
    atom_name_a: str = "",
    atom_name_b: str = "",
) -> Optional[Tuple[dict, dict]]:
    descriptors_a = _residue_ring_descriptors(
        residue_a,
        residue_atoms_index,
        nucleobase_only=False,
    )
    descriptors_b = _residue_ring_descriptors(
        residue_b,
        residue_atoms_index,
        nucleobase_only=False,
    )
    if not descriptors_a or not descriptors_b:
        return None

    best_pair: Optional[Tuple[dict, dict]] = None
    best_hint_misses = math.inf
    best_min_sq = math.inf
    best_centroid_sq = math.inf
    preferred_atoms_a = {_normalize_atom_name(token) for token in str(residue_a.get("atom") or atom_name_a).split(",") if token}
    preferred_atoms_b = {_normalize_atom_name(token) for token in str(residue_b.get("atom") or atom_name_b).split(",") if token}
    for descriptor_a in descriptors_a:
        points_a = descriptor_a.get("points")
        if not isinstance(points_a, list) or not points_a:
            continue
        for descriptor_b in descriptors_b:
            points_b = descriptor_b.get("points")
            if not isinstance(points_b, list) or not points_b:
                continue
            min_sq = _min_distance_sq_between_point_sets(points_a, points_b)
            if not math.isfinite(min_sq):
                continue
            centroid_sq = _centroid_distance_sq_for_point_sets(points_a, points_b)
            hint_misses = len(preferred_atoms_a - set(descriptor_a.get("atom_names", [])))
            hint_misses += len(preferred_atoms_b - set(descriptor_b.get("atom_names", [])))
            if (
                hint_misses + 1e-9 < best_hint_misses
                or (
                    abs(hint_misses - best_hint_misses) <= 1e-9
                    and (
                        min_sq + 1e-9 < best_min_sq
                        or (
                            abs(min_sq - best_min_sq) <= 1e-9
                            and centroid_sq + 1e-9 < best_centroid_sq
                        )
                    )
                )
            ):
                best_pair = (descriptor_a, descriptor_b)
                best_hint_misses = float(hint_misses)
                best_min_sq = min_sq
                best_centroid_sq = centroid_sq
    return best_pair


def _compute_ring_metrics_from_point_sets(
    points_a: List[Tuple[float, float, float]],
    points_b: List[Tuple[float, float, float]],
) -> Dict[str, float]:
    if not points_a or not points_b:
        return {}
    normal_a = _best_plane_normal(points_a)
    normal_b = _best_plane_normal(points_b)
    if not normal_a or not normal_b:
        return {}
    centroid_a = (
        sum(point[0] for point in points_a) / len(points_a),
        sum(point[1] for point in points_a) / len(points_a),
        sum(point[2] for point in points_a) / len(points_a),
    )
    centroid_b = (
        sum(point[0] for point in points_b) / len(points_b),
        sum(point[1] for point in points_b) / len(points_b),
        sum(point[2] for point in points_b) / len(points_b),
    )
    dx = centroid_b[0] - centroid_a[0]
    dy = centroid_b[1] - centroid_a[1]
    dz = centroid_b[2] - centroid_a[2]
    centroid_distance = math.sqrt(dx * dx + dy * dy + dz * dz)
    dot_raw = (
        normal_a[0] * normal_b[0]
        + normal_a[1] * normal_b[1]
        + normal_a[2] * normal_b[2]
    )
    dot = abs(dot_raw)
    dot = max(0.0, min(1.0, dot))
    normal_angle = math.degrees(math.acos(dot))

    aligned_normal_b = normal_b
    if dot_raw < 0:
        aligned_normal_b = (-normal_b[0], -normal_b[1], -normal_b[2])
    avg_nx = normal_a[0] + aligned_normal_b[0]
    avg_ny = normal_a[1] + aligned_normal_b[1]
    avg_nz = normal_a[2] + aligned_normal_b[2]
    avg_norm_sq = avg_nx * avg_nx + avg_ny * avg_ny + avg_nz * avg_nz
    if avg_norm_sq < 1e-10:
        avg_nx, avg_ny, avg_nz = normal_a
        avg_norm_sq = avg_nx * avg_nx + avg_ny * avg_ny + avg_nz * avg_nz
    inv_avg_len = 1.0 / math.sqrt(avg_norm_sq)
    avg_nx *= inv_avg_len
    avg_ny *= inv_avg_len
    avg_nz *= inv_avg_len
    interplanar_distance = abs(dx * avg_nx + dy * avg_ny + dz * avg_nz)
    lateral_sq = max(0.0, centroid_distance * centroid_distance - interplanar_distance * interplanar_distance)
    lateral_offset = math.sqrt(lateral_sq)
    min_sq = _min_distance_sq_between_point_sets(points_a, points_b)
    min_atom_distance = math.sqrt(min_sq) if math.isfinite(min_sq) else None

    payload: Dict[str, float] = {
        "ring_centroid_distance": centroid_distance,
        "ring_normal_angle": normal_angle,
        "ring_interplanar_distance": interplanar_distance,
        "ring_lateral_offset": lateral_offset,
    }
    if min_atom_distance is not None:
        payload["ring_min_atom_distance"] = min_atom_distance
    return payload


def _compute_ring_geometry_metrics_for_contact(
    residue_a: dict,
    residue_b: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
    *,
    atom_name_a: str = "",
    atom_name_b: str = "",
) -> Dict[str, float]:
    ring_pair = _select_ring_descriptor_pair_by_contact(
        residue_a,
        residue_b,
        residue_atoms_index,
        atom_name_a=atom_name_a,
        atom_name_b=atom_name_b,
    )
    if not ring_pair:
        return _compute_ring_geometry_metrics(
            residue_a,
            residue_b,
            residue_atoms_index,
        )
    descriptor_a, descriptor_b = ring_pair
    points_a = descriptor_a.get("points") if isinstance(descriptor_a, dict) else None
    points_b = descriptor_b.get("points") if isinstance(descriptor_b, dict) else None
    if not isinstance(points_a, list) or not isinstance(points_b, list):
        return _compute_ring_geometry_metrics(
            residue_a,
            residue_b,
            residue_atoms_index,
        )
    return _compute_ring_metrics_from_point_sets(points_a, points_b)


def _is_nonpolymer_residue_name(res_name: str) -> bool:
    token = str(res_name or "").strip().upper()
    return bool(
        token
        and token not in POLYMER_RESIDUES
        and token not in WATER_RESIDUES
        and token not in METAL_ELEMENTS
    )


def _is_contact_atom_in_residue_ring(
    residue: dict,
    res_name: str,
    atom_name: str,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> bool:
    atom_token = _normalize_atom_name(atom_name)
    if not atom_token:
        return False
    ring_names = _residue_ring_atom_names(res_name)
    if ring_names:
        return atom_token in ring_names
    if not _is_nonpolymer_residue_name(res_name):
        return False
    descriptors = _residue_ring_descriptors(
        residue,
        residue_atoms_index,
        nucleobase_only=False,
    )
    for descriptor in descriptors:
        if _ring_descriptor_contains_atom_name(descriptor, atom_token):
            return True
    return False


def _residue_ring_points(
    residue: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
    *,
    nucleobase_only: bool = False,
) -> List[Tuple[float, float, float]]:
    point_sets = _residue_ring_point_sets(
        residue,
        residue_atoms_index,
        nucleobase_only=nucleobase_only,
    )
    if not point_sets:
        return []
    return point_sets[0]


def _min_distance_sq_between_point_sets(
    points_a: List[Tuple[float, float, float]],
    points_b: List[Tuple[float, float, float]],
) -> float:
    min_distance_sq = math.inf
    for ax, ay, az in points_a:
        for bx, by, bz in points_b:
            dx = ax - bx
            dy = ay - by
            dz = az - bz
            distance_sq = dx * dx + dy * dy + dz * dz
            if distance_sq < min_distance_sq:
                min_distance_sq = distance_sq
    return min_distance_sq


def _centroid_distance_sq_for_point_sets(
    points_a: List[Tuple[float, float, float]],
    points_b: List[Tuple[float, float, float]],
) -> float:
    if not points_a or not points_b:
        return math.inf
    centroid_a = (
        sum(point[0] for point in points_a) / len(points_a),
        sum(point[1] for point in points_a) / len(points_a),
        sum(point[2] for point in points_a) / len(points_a),
    )
    centroid_b = (
        sum(point[0] for point in points_b) / len(points_b),
        sum(point[1] for point in points_b) / len(points_b),
        sum(point[2] for point in points_b) / len(points_b),
    )
    dx = centroid_b[0] - centroid_a[0]
    dy = centroid_b[1] - centroid_a[1]
    dz = centroid_b[2] - centroid_a[2]
    return dx * dx + dy * dy + dz * dz


def _select_ring_point_set_pair_by_proximity(
    point_sets_a: List[List[Tuple[float, float, float]]],
    point_sets_b: List[List[Tuple[float, float, float]]],
) -> Optional[Tuple[List[Tuple[float, float, float]], List[Tuple[float, float, float]]]]:
    best_pair: Optional[Tuple[List[Tuple[float, float, float]], List[Tuple[float, float, float]]]] = None
    best_min_sq = math.inf
    best_centroid_sq = math.inf
    for points_a in point_sets_a:
        if not points_a:
            continue
        for points_b in point_sets_b:
            if not points_b:
                continue
            min_sq = _min_distance_sq_between_point_sets(points_a, points_b)
            if not math.isfinite(min_sq):
                continue
            centroid_sq = _centroid_distance_sq_for_point_sets(points_a, points_b)
            if (
                min_sq + 1e-9 < best_min_sq
                or (
                    abs(min_sq - best_min_sq) <= 1e-9
                    and centroid_sq + 1e-9 < best_centroid_sq
                )
            ):
                best_pair = (points_a, points_b)
                best_min_sq = min_sq
                best_centroid_sq = centroid_sq
    return best_pair


def _residue_ring_geometry(
    residue: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]:
    points = _residue_ring_points(
        residue,
        residue_atoms_index,
        nucleobase_only=False,
    )
    if len(points) < 4:
        return None
    normal = _best_plane_normal(points)
    if not normal:
        return None
    centroid = (
        sum(point[0] for point in points) / len(points),
        sum(point[1] for point in points) / len(points),
        sum(point[2] for point in points) / len(points),
    )
    return centroid, normal


def _compute_ring_geometry_metrics(
    residue_a: dict,
    residue_b: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Dict[str, float]:
    point_sets_a = _residue_ring_point_sets(
        residue_a,
        residue_atoms_index,
        nucleobase_only=False,
    )
    point_sets_b = _residue_ring_point_sets(
        residue_b,
        residue_atoms_index,
        nucleobase_only=False,
    )
    ring_pair = _select_ring_point_set_pair_by_proximity(point_sets_a, point_sets_b)
    if not ring_pair:
        return {}
    points_a, points_b = ring_pair
    normal_a = _best_plane_normal(points_a)
    normal_b = _best_plane_normal(points_b)
    if not normal_a or not normal_b:
        return {}
    centroid_a = (
        sum(point[0] for point in points_a) / len(points_a),
        sum(point[1] for point in points_a) / len(points_a),
        sum(point[2] for point in points_a) / len(points_a),
    )
    centroid_b = (
        sum(point[0] for point in points_b) / len(points_b),
        sum(point[1] for point in points_b) / len(points_b),
        sum(point[2] for point in points_b) / len(points_b),
    )

    dx = centroid_b[0] - centroid_a[0]
    dy = centroid_b[1] - centroid_a[1]
    dz = centroid_b[2] - centroid_a[2]
    centroid_distance = math.sqrt(dx * dx + dy * dy + dz * dz)

    dot_raw = (
        normal_a[0] * normal_b[0]
        + normal_a[1] * normal_b[1]
        + normal_a[2] * normal_b[2]
    )
    dot = abs(dot_raw)
    dot = max(0.0, min(1.0, dot))
    normal_angle = math.degrees(math.acos(dot))

    aligned_normal_b = normal_b
    if dot_raw < 0:
        aligned_normal_b = (-normal_b[0], -normal_b[1], -normal_b[2])
    avg_nx = normal_a[0] + aligned_normal_b[0]
    avg_ny = normal_a[1] + aligned_normal_b[1]
    avg_nz = normal_a[2] + aligned_normal_b[2]
    avg_norm_sq = avg_nx * avg_nx + avg_ny * avg_ny + avg_nz * avg_nz
    if avg_norm_sq < 1e-10:
        avg_nx, avg_ny, avg_nz = normal_a
        avg_norm_sq = avg_nx * avg_nx + avg_ny * avg_ny + avg_nz * avg_nz
    inv_avg_len = 1.0 / math.sqrt(avg_norm_sq)
    avg_nx *= inv_avg_len
    avg_ny *= inv_avg_len
    avg_nz *= inv_avg_len

    interplanar_distance = abs(dx * avg_nx + dy * avg_ny + dz * avg_nz)
    lateral_sq = max(0.0, centroid_distance * centroid_distance - interplanar_distance * interplanar_distance)
    lateral_offset = math.sqrt(lateral_sq)
    return {
        "ring_centroid_distance": centroid_distance,
        "ring_normal_angle": normal_angle,
        "ring_interplanar_distance": interplanar_distance,
        "ring_lateral_offset": lateral_offset,
    }


def _round_metric_or_none(value: Optional[float]) -> Optional[float]:
    number = _coerce_float(value)
    if number is None or not math.isfinite(number):
        return None
    return round(number, 3)


def _build_ring_metrics_payload(
    *,
    centroid_distance: Optional[float],
    min_atom_distance: Optional[float],
    interplanar_distance: Optional[float],
    lateral_offset: Optional[float],
    normal_angle: Optional[float],
) -> dict:
    return {
        "centroid_distance": _round_metric_or_none(centroid_distance),
        "min_atom_distance": _round_metric_or_none(min_atom_distance),
        "interplanar_distance": _round_metric_or_none(interplanar_distance),
        "lateral_offset": _round_metric_or_none(lateral_offset),
        "normal_angle": _round_metric_or_none(normal_angle),
    }


def _resolve_ring_display_distance_for_family(
    family: str,
    ring_payload: Optional[dict],
) -> Tuple[Optional[float], Optional[str]]:
    if not isinstance(ring_payload, dict):
        return None, None
    family_key = str(family or "").strip().lower()
    centroid = _coerce_float(ring_payload.get("centroid_distance"))
    min_atom = _coerce_float(ring_payload.get("min_atom_distance"))
    if family_key in {"pi_pi", "pi_cation", "water_pi"}:
        return centroid, "ring_centroid_distance"
    if family_key == "aromatic_packing":
        return min_atom, "ring_min_atom_distance"
    if family_key == "aromatic_proximal":
        if centroid is not None:
            return centroid, "ring_centroid_distance"
        return min_atom, "ring_min_atom_distance"
    return None, None


def _resolve_contact_atom_record_from_payload(
    residue: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Optional[AtomRecord]:
    if not isinstance(residue, dict) or not residue_atoms_index:
        return None
    chain = str(residue.get("chain") or "").strip()
    seq = str(residue.get("seq") or "").strip()
    atom_name = _normalize_atom_name(residue.get("atom"))
    if not chain or not seq or not atom_name:
        return None
    lookup = getattr(residue_atoms_index, "atom_lookup", None)
    if lookup is not None:
        return lookup.get((chain, seq, atom_name))
    atoms = residue_atoms_index.get((chain, seq), [])
    if not atoms:
        return None
    for atom in atoms:
        if _normalize_atom_name(atom.atom_name) == atom_name:
            return atom
    return None


def _distance_between_payload_atoms(
    residue_a: dict,
    residue_b: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Optional[float]:
    atom_a = _resolve_contact_atom_record_from_payload(residue_a, residue_atoms_index)
    atom_b = _resolve_contact_atom_record_from_payload(residue_b, residue_atoms_index)
    if atom_a is None or atom_b is None:
        return None
    dx = atom_a.x - atom_b.x
    dy = atom_a.y - atom_b.y
    dz = atom_a.z - atom_b.z
    distance = math.sqrt(dx * dx + dy * dy + dz * dz)
    if not math.isfinite(distance):
        return None
    return distance


def _compute_nucleobase_ring_min_atom_distance(
    residue_a: dict,
    residue_b: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Optional[float]:
    if not isinstance(residue_a, dict) or not isinstance(residue_b, dict):
        return None
    if not residue_atoms_index:
        return None

    points_a = _residue_ring_points(
        residue_a,
        residue_atoms_index,
        nucleobase_only=True,
    )
    points_b = _residue_ring_points(
        residue_b,
        residue_atoms_index,
        nucleobase_only=True,
    )
    if not points_a or not points_b:
        return None

    min_distance_sq = math.inf
    for ax, ay, az in points_a:
        for bx, by, bz in points_b:
            dx = ax - bx
            dy = ay - by
            dz = az - bz
            distance_sq = dx * dx + dy * dy + dz * dz
            if distance_sq < min_distance_sq:
                min_distance_sq = distance_sq
    if not math.isfinite(min_distance_sq):
        return None
    if min_distance_sq < 0:
        return None
    distance = math.sqrt(min_distance_sq)
    if not math.isfinite(distance):
        return None
    return distance


def _compute_ring_min_atom_distance(
    residue_a: dict,
    residue_b: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Optional[float]:
    if not isinstance(residue_a, dict) or not isinstance(residue_b, dict):
        return None
    if not residue_atoms_index:
        return None
    point_sets_a = _residue_ring_point_sets(
        residue_a,
        residue_atoms_index,
        nucleobase_only=False,
    )
    point_sets_b = _residue_ring_point_sets(
        residue_b,
        residue_atoms_index,
        nucleobase_only=False,
    )
    if not point_sets_a or not point_sets_b:
        return None
    min_distance_sq = math.inf
    for points_a in point_sets_a:
        if not points_a:
            continue
        for points_b in point_sets_b:
            if not points_b:
                continue
            distance_sq = _min_distance_sq_between_point_sets(points_a, points_b)
            if distance_sq < min_distance_sq:
                min_distance_sq = distance_sq
    if not math.isfinite(min_distance_sq) or min_distance_sq < 0:
        return None
    distance = math.sqrt(min_distance_sq)
    if not math.isfinite(distance):
        return None
    return distance


def _same_atom_payload_endpoint(residue_a: dict, residue_b: dict) -> bool:
    if not isinstance(residue_a, dict) or not isinstance(residue_b, dict):
        return False
    chain_a = str(residue_a.get("chain") or "").strip()
    chain_b = str(residue_b.get("chain") or "").strip()
    seq_a = str(residue_a.get("seq") or "").strip()
    seq_b = str(residue_b.get("seq") or "").strip()
    atom_a = _normalize_atom_name(residue_a.get("atom"))
    atom_b = _normalize_atom_name(residue_b.get("atom"))
    if not chain_a or not chain_b or not seq_a or not seq_b or not atom_a or not atom_b:
        return False
    return chain_a == chain_b and seq_a == seq_b and atom_a == atom_b


def _resolve_contact_distance_value(
    raw_contact: dict,
    residue_a: dict,
    residue_b: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Tuple[Optional[float], Optional[str]]:
    raw_distance = _coerce_float(raw_contact.get("distance"))
    if raw_distance is not None and raw_distance > 0.0:
        return raw_distance, None
    if raw_distance is not None and raw_distance <= 0.0 and _same_atom_payload_endpoint(residue_a, residue_b):
        return raw_distance, None
    recomputed = _distance_between_payload_atoms(residue_a, residue_b, residue_atoms_index)
    if recomputed is not None and recomputed > 0.0:
        return recomputed, "distance_recomputed_from_coordinates"
    return None, "distance_missing_or_invalid"


def _angle_between_vectors_degrees(
    vec_a: Tuple[float, float, float],
    vec_b: Tuple[float, float, float],
) -> Optional[float]:
    ax, ay, az = vec_a
    bx, by, bz = vec_b
    len_a = math.sqrt(ax * ax + ay * ay + az * az)
    len_b = math.sqrt(bx * bx + by * by + bz * bz)
    if len_a <= 1e-8 or len_b <= 1e-8:
        return None
    dot = (ax * bx + ay * by + az * bz) / (len_a * len_b)
    dot = max(-1.0, min(1.0, dot))
    angle = math.degrees(math.acos(dot))
    if not math.isfinite(angle):
        return None
    return angle


def _resolve_hbond_role_candidates(
    *,
    res_name_a: str,
    atom_name_a: str,
    element_a: str,
    res_name_b: str,
    atom_name_b: str,
    element_b: str,
) -> List[Tuple[str, str]]:
    a_donor = _is_hbond_donor_capable(
        res_name=res_name_a,
        atom_name=atom_name_a,
        element=element_a,
    )
    b_donor = _is_hbond_donor_capable(
        res_name=res_name_b,
        atom_name=atom_name_b,
        element=element_b,
    )
    a_acceptor = _is_hbond_acceptor_capable(
        res_name=res_name_a,
        atom_name=atom_name_a,
        element=element_a,
    )
    b_acceptor = _is_hbond_acceptor_capable(
        res_name=res_name_b,
        atom_name=atom_name_b,
        element=element_b,
    )
    roles: List[Tuple[str, str]] = []
    if a_donor and b_acceptor:
        roles.append(("A", "B"))
    if b_donor and a_acceptor:
        roles.append(("B", "A"))
    return roles


def _resolve_residue_atom_record_by_name(
    residue: dict,
    atom_name: str,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Optional[AtomRecord]:
    if not isinstance(residue, dict):
        return None
    atom = _normalize_atom_name(atom_name)
    if not atom:
        return None
    payload = dict(residue)
    payload["atom"] = atom
    return _resolve_contact_atom_record_from_payload(payload, residue_atoms_index)


def _resolve_halogen_donor_anchor_atom(
    donor_residue: dict,
    donor_atom_name: str,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Optional[AtomRecord]:
    donor_atom = _resolve_residue_atom_record_by_name(
        donor_residue,
        donor_atom_name,
        residue_atoms_index,
    )
    if donor_atom is None:
        return None
    donor_element = str(donor_atom.element or "").strip().upper() or guess_element(donor_atom.atom_name).upper()
    if donor_element not in HALOGEN_BOND_DONOR_ELEMENTS:
        return None
    node = _semantic_node(dict(donor_residue, atom=donor_atom_name), {}, residue_atoms_index)
    if 'bonded_atoms' in node:
        for neighbor in node['bonded_atoms']:
            if neighbor.get('type_symbol') == 'C':
                return _resolve_residue_atom_record_by_name(donor_residue, neighbor.get('auth_atom_id'), residue_atoms_index)
        return None
    candidate_atoms = _residue_atoms_from_payload(donor_residue, residue_atoms_index)
    best_anchor: Optional[AtomRecord] = None
    best_distance = math.inf
    donor_atom_name_token = _normalize_atom_name(donor_atom.atom_name)
    for candidate in candidate_atoms:
        candidate_name = _normalize_atom_name(candidate.atom_name)
        if not candidate_name or candidate_name == donor_atom_name_token:
            continue
        candidate_element = str(candidate.element or "").strip().upper() or guess_element(candidate_name).upper()
        if candidate_element != "C":
            continue
        d = math.sqrt(_distance_sq(donor_atom, candidate))
        if not math.isfinite(d):
            continue
        if d < HALOGEN_DONOR_BOND_MIN_DISTANCE or d > HALOGEN_DONOR_BOND_MAX_DISTANCE:
            continue
        if d < best_distance:
            best_distance = d
            best_anchor = candidate
    return best_anchor


def _compute_halogen_bond_angle(
    *,
    donor_residue: dict,
    donor_atom_name: str,
    acceptor_residue: dict,
    acceptor_atom_name: str,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Optional[float]:
    donor_atom = _resolve_residue_atom_record_by_name(
        donor_residue,
        donor_atom_name,
        residue_atoms_index,
    )
    acceptor_atom = _resolve_residue_atom_record_by_name(
        acceptor_residue,
        acceptor_atom_name,
        residue_atoms_index,
    )
    anchor_atom = _resolve_halogen_donor_anchor_atom(
        donor_residue,
        donor_atom_name,
        residue_atoms_index,
    )
    if donor_atom is None or acceptor_atom is None or anchor_atom is None:
        return None
    donor_to_anchor = (
        anchor_atom.x - donor_atom.x,
        anchor_atom.y - donor_atom.y,
        anchor_atom.z - donor_atom.z,
    )
    donor_to_acceptor = (
        acceptor_atom.x - donor_atom.x,
        acceptor_atom.y - donor_atom.y,
        acceptor_atom.z - donor_atom.z,
    )
    return _angle_between_vectors_degrees(donor_to_anchor, donor_to_acceptor)


def _resolve_hbond_donor_antecedent_atom_name(
    *,
    donor_res_name: str,
    donor_atom_name: str,
    donor_element: str,
) -> str:
    residue = str(donor_res_name or "").strip().upper()
    atom = _normalize_atom_name(donor_atom_name)
    element = str(donor_element or "").strip().upper() or guess_element(atom).upper()
    if element != "O":
        return ""
    antecedents = HBOND_DONOR_ANTECEDENT_BY_RESIDUE.get(residue) or {}
    return str(antecedents.get(atom) or "").strip().upper()


def _compute_hydroxyl_hbond_proxy_angle(
    *,
    donor_residue: dict,
    donor_res_name: str,
    donor_atom_name: str,
    donor_element: str,
    acceptor_residue: dict,
    acceptor_atom_name: str,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Optional[float]:
    antecedent_name = _resolve_hbond_donor_antecedent_atom_name(
        donor_res_name=donor_res_name,
        donor_atom_name=donor_atom_name,
        donor_element=donor_element,
    )
    if not antecedent_name:
        return None

    donor_atom = _resolve_residue_atom_record_by_name(
        donor_residue,
        donor_atom_name,
        residue_atoms_index,
    )
    acceptor_atom = _resolve_residue_atom_record_by_name(
        acceptor_residue,
        acceptor_atom_name,
        residue_atoms_index,
    )
    antecedent_atom = _resolve_residue_atom_record_by_name(
        donor_residue,
        antecedent_name,
        residue_atoms_index,
    )
    if donor_atom is None or acceptor_atom is None or antecedent_atom is None:
        return None

    # Diagnostic heavy-atom axis only: this is not a donor-H-acceptor angle.
    donor_to_inferred_h = (
        donor_atom.x - antecedent_atom.x,
        donor_atom.y - antecedent_atom.y,
        donor_atom.z - antecedent_atom.z,
    )
    donor_to_acceptor = (
        acceptor_atom.x - donor_atom.x,
        acceptor_atom.y - donor_atom.y,
        acceptor_atom.z - donor_atom.z,
    )
    donor_frame_angle = _angle_between_vectors_degrees(donor_to_inferred_h, donor_to_acceptor)
    if donor_frame_angle is None:
        return None
    # Approximate D-H...A from donor-frame vectors.
    proxy_dha = 180.0 - donor_frame_angle
    if not math.isfinite(proxy_dha):
        return None
    return max(0.0, min(180.0, proxy_dha))


def _compute_hbond_proxy_angle(
    *,
    residue_a: dict,
    residue_b: dict,
    res_name_a: str,
    atom_name_a: str,
    element_a: str,
    res_name_b: str,
    atom_name_b: str,
    element_b: str,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Tuple[Optional[float], Optional[str]]:
    roles = _resolve_hbond_role_candidates(
        res_name_a=res_name_a,
        atom_name_a=atom_name_a,
        element_a=element_a,
        res_name_b=res_name_b,
        atom_name_b=atom_name_b,
        element_b=element_b,
    )
    if not roles:
        return None, None

    best_angle: Optional[float] = None
    best_method: Optional[str] = None
    for donor_side, acceptor_side in roles:
        if donor_side == "A":
            angle = _compute_hydroxyl_hbond_proxy_angle(
                donor_residue=residue_a,
                donor_res_name=res_name_a,
                donor_atom_name=atom_name_a,
                donor_element=element_a,
                acceptor_residue=residue_b,
                acceptor_atom_name=atom_name_b,
                residue_atoms_index=residue_atoms_index,
            )
        else:
            angle = _compute_hydroxyl_hbond_proxy_angle(
                donor_residue=residue_b,
                donor_res_name=res_name_b,
                donor_atom_name=atom_name_b,
                donor_element=element_b,
                acceptor_residue=residue_a,
                acceptor_atom_name=atom_name_a,
                residue_atoms_index=residue_atoms_index,
            )
        if angle is None:
            continue
        if best_angle is None or angle > best_angle:
            best_angle = angle
            best_method = "heavy_atom_axis"

    return best_angle, best_method


def _is_hbond_donor_capable(
    *,
    res_name: str,
    atom_name: str,
    element: str,
) -> bool:
    residue = str(res_name or "").strip().upper()
    atom = _normalize_atom_name(atom_name)
    atom_element = str(element or "").strip().upper()
    if not atom_element:
        atom_element = guess_element(atom).upper()

    if atom_element == "N":
        if residue in STANDARD_AMINO_RESIDUES:
            if atom == "N":
                return residue != "PRO"
            if residue == "LYS":
                return atom == "NZ"
            if residue == "ARG":
                return atom in {"NE", "NH1", "NH2"}
            if residue == "HIS":
                return atom in {"ND1", "NE2"}
            if residue == "ASN":
                return atom == "ND2"
            if residue == "GLN":
                return atom == "NE2"
            if residue == "TRP":
                return atom == "NE1"
            return False
        if _nucleic_base_family(residue):
            return _is_nucleobase_donor_nitrogen(residue, atom)
        if residue in WATER_RESIDUES:
            return False
        # Ligand valence/protonation requires Arpeggio/Open Babel atom typing.
        return False

    if atom_element == "O":
        if residue in WATER_RESIDUES:
            return True
        if residue == "SER":
            return atom == "OG"
        if residue == "THR":
            return atom == "OG1"
        if residue == "TYR":
            return atom == "OH"
        if _nucleic_base_family(residue):
            return atom in NUCLEIC_SUGAR_DONOR_OXYGENS
        return False

    if atom_element in {"S", "SE"}:
        return (residue == "CYS" and atom == "SG") or (residue == "SEC" and atom == "SE")

    return False


def _is_hbond_acceptor_capable(
    *,
    res_name: str,
    atom_name: str,
    element: str,
) -> bool:
    residue = str(res_name or "").strip().upper()
    atom = _normalize_atom_name(atom_name)
    atom_element = str(element or "").strip().upper()
    if not atom_element:
        atom_element = guess_element(atom).upper()

    if atom_element == "O":
        if residue in STANDARD_AMINO_RESIDUES:
            if atom in {"O", "OXT", "OT1", "OT2"}:
                return True
            if residue == "ASP":
                return atom in {"OD1", "OD2"}
            if residue == "GLU":
                return atom in {"OE1", "OE2"}
            if residue == "ASN":
                return atom == "OD1"
            if residue == "GLN":
                return atom == "OE1"
            if residue == "SER":
                return atom == "OG"
            if residue == "THR":
                return atom == "OG1"
            if residue == "TYR":
                return atom == "OH"
            if residue == "CYS":
                return atom == "SG"
            if residue == "SEC":
                return atom == "SE"
            if residue == "MET":
                return atom == "SD"
            return atom.startswith("O")
        if _nucleic_base_family(residue):
            if _is_phosphate_oxygen_atom_name(atom):
                return True
            if _is_nucleobase_acceptor_oxygen(atom):
                return True
            return atom in NUCLEIC_SUGAR_ACCEPTOR_OXYGENS or atom.startswith("O")
        return residue in WATER_RESIDUES

    if atom_element == "N":
        if residue in STANDARD_AMINO_RESIDUES:
            if atom == "N":
                return False
            if residue == "HIS":
                return atom in {"ND1", "NE2"}
            if residue == "LYS":
                return atom != "NZ"
            if residue == "ARG":
                return atom not in {"NE", "NH1", "NH2"}
            if residue == "ASN":
                return atom != "ND2"
            if residue == "GLN":
                return atom != "NE2"
            if residue == "TRP":
                return atom != "NE1"
            return False
        if _nucleic_base_family(residue):
            return _is_nucleobase_acceptor_nitrogen(residue, atom)
        return False

    if atom_element in {"S", "SE"}:
        if residue == "MET":
            return atom == "SD"
        if residue == "CYS":
            return atom == "SG"
        if residue == "SEC":
            return atom == "SE"
        return False

    return False


def _is_weak_hbond_acceptor_site(
    *,
    res_name: str,
    atom_name: str,
    element: str,
) -> bool:
    residue = str(res_name or "").strip().upper()
    atom = _normalize_atom_name(atom_name)
    atom_element = str(element or "").strip().upper()
    if not atom_element:
        atom_element = guess_element(atom).upper()
    if atom_element != "O":
        return False
    # Sidechain hydroxyl oxygens are weak acceptors and require directional evidence.
    if residue == "SER":
        return atom == "OG"
    if residue == "THR":
        return atom == "OG1"
    if residue == "TYR":
        return atom == "OH"
    return False


def _is_histidine_protonation_dependent_donor_site(
    *,
    res_name: str,
    atom_name: str,
    element: str,
) -> bool:
    residue = str(res_name or "").strip().upper()
    atom = _normalize_atom_name(atom_name)
    atom_element = str(element or "").strip().upper()
    if not atom_element:
        atom_element = guess_element(atom).upper()
    return residue == "HIS" and atom_element == "N" and atom in {"ND1", "NE2"}


def _is_hbond_donor_acceptor_pair(
    *,
    res_name_a: str,
    atom_name_a: str,
    element_a: str,
    res_name_b: str,
    atom_name_b: str,
    element_b: str,
) -> bool:
    a_donor = _is_hbond_donor_capable(
        res_name=res_name_a,
        atom_name=atom_name_a,
        element=element_a,
    )
    b_donor = _is_hbond_donor_capable(
        res_name=res_name_b,
        atom_name=atom_name_b,
        element=element_b,
    )
    a_acceptor = _is_hbond_acceptor_capable(
        res_name=res_name_a,
        atom_name=atom_name_a,
        element=element_a,
    )
    b_acceptor = _is_hbond_acceptor_capable(
        res_name=res_name_b,
        atom_name=atom_name_b,
        element=element_b,
    )
    return (a_donor and b_acceptor) or (b_donor and a_acceptor)


def _is_canonical_base_pair_hbond_pair(
    *,
    res_name_a: str,
    atom_name_a: str,
    element_a: str,
    res_name_b: str,
    atom_name_b: str,
    element_b: str,
) -> bool:
    family_a = _nucleic_base_family(res_name_a)
    family_b = _nucleic_base_family(res_name_b)
    if not family_a or not family_b:
        return False

    atom_a = _normalize_atom_name(atom_name_a)
    atom_b = _normalize_atom_name(atom_name_b)
    elem_a = str(element_a or "").strip().upper() or guess_element(atom_a).upper()
    elem_b = str(element_b or "").strip().upper() or guess_element(atom_b).upper()
    if elem_a not in {"N", "O"} or elem_b not in {"N", "O"}:
        return False

    # Canonical Watson-Crick atom pairs:
    # A-U/T: A(N6-H)->U/T(O4), U/T(N3-H)->A(N1)
    # G-C:   G(N1-H)->C(N3), G(N2-H)->C(O2), C(N4-H)->G(O6)
    canonical_pairs = {
        ("A", "N6", "U", "O4"),
        ("U", "N3", "A", "N1"),
        ("A", "N6", "T", "O4"),
        ("T", "N3", "A", "N1"),
        ("G", "N1", "C", "N3"),
        ("G", "N2", "C", "O2"),
        ("C", "N4", "G", "O6"),
    }
    return (
        (family_a, atom_a, family_b, atom_b) in canonical_pairs
        or (family_b, atom_b, family_a, atom_a) in canonical_pairs
    )


def _is_salt_bridge_cation_site(
    res_name: str,
    atom_name: str,
    element: str,
) -> bool:
    residue = str(res_name or "").strip().upper()
    atom = _normalize_atom_name(atom_name)
    atom_element = str(element or "").strip().upper()
    if not atom_element:
        atom_element = guess_element(atom).upper()
    if atom_element in METAL_ELEMENTS or residue in METAL_ELEMENTS:
        return True
    if atom_element != "N":
        return False
    if residue in {"LYS"}:
        return atom == "NZ"
    if residue in {"ARG"}:
        return atom in {"NE", "NH1", "NH2"}
    if residue in PROTONATED_HISTIDINE_RESIDUE_NAMES:
        return atom in {"ND1", "NE2"}
    return False


def _is_salt_bridge_anion_site(
    res_name: str,
    atom_name: str,
    element: str,
) -> bool:
    residue = str(res_name or "").strip().upper()
    atom = _normalize_atom_name(atom_name)
    atom_element = str(element or "").strip().upper()
    if not atom_element:
        atom_element = guess_element(atom).upper()
    if atom_element != "O":
        return False
    if residue == "ASP":
        return atom in {"OD1", "OD2"}
    if residue == "GLU":
        return atom in {"OE1", "OE2"}
    if _is_phosphate_oxygen_atom_name(atom):
        return True
    if _is_sulfate_oxygen_site(residue, atom, atom_element):
        return True
    if residue not in POLYMER_RESIDUES and "-" in residue:
        return True
    return False


def _is_valid_salt_bridge_pair(
    *,
    res_name_a: str,
    atom_name_a: str,
    element_a: str,
    res_name_b: str,
    atom_name_b: str,
    element_b: str,
) -> bool:
    a_is_cation = _is_salt_bridge_cation_site(res_name_a, atom_name_a, element_a)
    b_is_cation = _is_salt_bridge_cation_site(res_name_b, atom_name_b, element_b)
    a_is_anion = _is_salt_bridge_anion_site(res_name_a, atom_name_a, element_a)
    b_is_anion = _is_salt_bridge_anion_site(res_name_b, atom_name_b, element_b)
    return (a_is_cation and b_is_anion) or (b_is_cation and a_is_anion)


def _residue_atoms_from_payload(
    residue: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> List[AtomRecord]:
    if not isinstance(residue, dict) or not residue_atoms_index:
        return []
    chain = str(residue.get("chain") or "").strip()
    seq = str(residue.get("seq") or "").strip()
    if not chain or not seq:
        return []
    atoms = residue_atoms_index.get((chain, seq), [])
    return list(atoms) if atoms else []


def _collect_salt_bridge_site_atoms(
    *,
    residue: dict,
    res_name: str,
    site_kind: str,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> List[AtomRecord]:
    atoms = _residue_atoms_from_payload(residue, residue_atoms_index)
    if not atoms:
        return []
    residue_name = str(res_name or "").strip().upper()
    if not residue_name:
        return []

    by_atom_name: Dict[str, AtomRecord] = {}
    for atom in atoms:
        atom_name = _normalize_atom_name(atom.atom_name)
        if atom_name and atom_name not in by_atom_name:
            by_atom_name[atom_name] = atom

    if site_kind == "cation":
        ordered = SALT_BRIDGE_CATION_SITE_ATOMS_BY_RESIDUE.get(residue_name) or ()
        selected: List[AtomRecord] = []
        for atom_name in ordered:
            atom = by_atom_name.get(atom_name)
            if atom is None:
                continue
            element = str(atom.element or "").strip().upper() or guess_element(atom_name).upper()
            if _is_salt_bridge_cation_site(residue_name, atom_name, element):
                selected.append(atom)
        if selected:
            return selected
        if residue_name not in POLYMER_RESIDUES and "+" in residue_name:
            return [
                atom
                for atom in atoms
                if _is_salt_bridge_cation_site(
                    residue_name,
                    atom.atom_name,
                    str(atom.element or "").strip().upper(),
                )
            ]
        return []

    if site_kind == "anion":
        ordered = SALT_BRIDGE_ANION_SITE_ATOMS_BY_RESIDUE.get(residue_name) or ()
        selected = []
        for atom_name in ordered:
            atom = by_atom_name.get(atom_name)
            if atom is None:
                continue
            element = str(atom.element or "").strip().upper() or guess_element(atom_name).upper()
            if _is_salt_bridge_anion_site(residue_name, atom_name, element):
                selected.append(atom)
        if selected:
            return selected

        phosphate_candidates = [
            atom
            for atom in atoms
            if _is_salt_bridge_anion_site(
                residue_name,
                atom.atom_name,
                str(atom.element or "").strip().upper(),
            )
            and _is_phosphate_oxygen_atom_name(atom.atom_name)
        ]
        if phosphate_candidates:
            return sorted(
                phosphate_candidates,
                key=lambda atom: _normalize_atom_name(atom.atom_name),
            )

        sulfate_candidates = [
            atom
            for atom in atoms
            if _is_sulfate_oxygen_site(
                residue_name,
                atom.atom_name,
                str(atom.element or "").strip().upper(),
            )
        ]
        if sulfate_candidates:
            return sorted(
                sulfate_candidates,
                key=lambda atom: _normalize_atom_name(atom.atom_name),
            )

        if residue_name not in POLYMER_RESIDUES and "-" in residue_name:
            return [
                atom
                for atom in atoms
                if _is_salt_bridge_anion_site(
                    residue_name,
                    atom.atom_name,
                    str(atom.element or "").strip().upper(),
                )
            ]
    return []


def _best_salt_bridge_site_pair(
    cation_atoms: List[AtomRecord],
    anion_atoms: List[AtomRecord],
) -> Optional[Tuple[AtomRecord, AtomRecord, float]]:
    if not cation_atoms or not anion_atoms:
        return None
    best_pair: Optional[Tuple[AtomRecord, AtomRecord, float]] = None
    best_distance_sq = math.inf
    best_tiebreak: Tuple[str, str] = ("", "")
    for cation_atom in cation_atoms:
        cation_name = _normalize_atom_name(cation_atom.atom_name)
        for anion_atom in anion_atoms:
            anion_name = _normalize_atom_name(anion_atom.atom_name)
            distance_sq = _distance_sq(cation_atom, anion_atom)
            if not math.isfinite(distance_sq):
                continue
            tiebreak = (cation_name, anion_name)
            if (
                distance_sq + 1e-10 < best_distance_sq
                or (
                    abs(distance_sq - best_distance_sq) <= 1e-10
                    and (best_pair is None or tiebreak < best_tiebreak)
                )
            ):
                best_distance_sq = distance_sq
                best_tiebreak = tiebreak
                best_pair = (cation_atom, anion_atom, math.sqrt(distance_sq))
    return best_pair


def _resolve_salt_bridge_endpoint_override(
    *,
    residue_a: dict,
    residue_b: dict,
    res_name_a: str,
    res_name_b: str,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Optional[dict]:
    if not residue_atoms_index:
        return None
    cation_atoms_a = _collect_salt_bridge_site_atoms(
        residue=residue_a,
        res_name=res_name_a,
        site_kind="cation",
        residue_atoms_index=residue_atoms_index,
    )
    anion_atoms_a = _collect_salt_bridge_site_atoms(
        residue=residue_a,
        res_name=res_name_a,
        site_kind="anion",
        residue_atoms_index=residue_atoms_index,
    )
    cation_atoms_b = _collect_salt_bridge_site_atoms(
        residue=residue_b,
        res_name=res_name_b,
        site_kind="cation",
        residue_atoms_index=residue_atoms_index,
    )
    anion_atoms_b = _collect_salt_bridge_site_atoms(
        residue=residue_b,
        res_name=res_name_b,
        site_kind="anion",
        residue_atoms_index=residue_atoms_index,
    )

    candidates: List[Tuple[str, str, AtomRecord, AtomRecord, float]] = []
    pair_a_to_b = _best_salt_bridge_site_pair(cation_atoms_a, anion_atoms_b)
    if pair_a_to_b is not None:
        cation_atom, anion_atom, best_distance = pair_a_to_b
        candidates.append(("A", "B", cation_atom, anion_atom, best_distance))
    pair_b_to_a = _best_salt_bridge_site_pair(cation_atoms_b, anion_atoms_a)
    if pair_b_to_a is not None:
        cation_atom, anion_atom, best_distance = pair_b_to_a
        candidates.append(("B", "A", cation_atom, anion_atom, best_distance))

    if not candidates:
        return None

    cation_side, anion_side, cation_atom, anion_atom, best_distance = min(
        candidates,
        key=lambda item: (
            item[4],
            _normalize_atom_name(item[2].atom_name),
            _normalize_atom_name(item[3].atom_name),
            item[0],
        ),
    )
    if cation_side == "A":
        atom_name_a = _normalize_atom_name(cation_atom.atom_name)
        element_a = str(cation_atom.element or "").strip().upper() or guess_element(atom_name_a).upper()
        atom_name_b = _normalize_atom_name(anion_atom.atom_name)
        element_b = str(anion_atom.element or "").strip().upper() or guess_element(atom_name_b).upper()
    else:
        atom_name_a = _normalize_atom_name(anion_atom.atom_name)
        element_a = str(anion_atom.element or "").strip().upper() or guess_element(atom_name_a).upper()
        atom_name_b = _normalize_atom_name(cation_atom.atom_name)
        element_b = str(cation_atom.element or "").strip().upper() or guess_element(atom_name_b).upper()

    if not atom_name_a or not atom_name_b:
        return None
    return {
        "atomNameA": atom_name_a,
        "elementA": element_a,
        "atomNameB": atom_name_b,
        "elementB": element_b,
        "distance": best_distance if math.isfinite(best_distance) else None,
        "cationSide": cation_side,
        "anionSide": anion_side,
    }


def _primary_contact_atom_name(atom_name: object) -> str:
    text = str(atom_name or "").strip()
    if "," in text:
        text = text.split(",", 1)[0].strip()
    return _normalize_atom_name(text)


def _hbond_distance_limit(element_a: str, element_b: str) -> float:
    elem_a = str(element_a or "").strip().upper()
    elem_b = str(element_b or "").strip().upper()
    if elem_a in {"S", "SE"} or elem_b in {"S", "SE"}:
        return HBOND_MAX_DISTANCE_CHALCOGEN
    pair = {elem_a, elem_b}
    if pair == {"O"}:
        return HBOND_MAX_DISTANCE_OO
    if pair == {"N"}:
        return HBOND_MAX_DISTANCE_NN
    if pair == {"N", "O"}:
        return HBOND_MAX_DISTANCE_NO
    return HBOND_POLAR_FALLBACK_MAX_DISTANCE


def _hbond_distance_is_within_limits(distance: Optional[float], element_a: str, element_b: str) -> bool:
    if distance is None:
        return False
    pair_floor = _pair_specific_min_nonbonded_distance(element_a, element_b)
    minimum_distance = HBOND_HEAVY_MIN_DISTANCE
    if pair_floor is not None:
        minimum_distance = max(minimum_distance, pair_floor)
    if distance < minimum_distance:
        return False
    return distance <= _hbond_distance_limit(element_a, element_b)


def _halogen_bond_distance_limit(halogen_element: str) -> float:
    token = str(halogen_element or "").strip().upper()
    if not token:
        return HALOGEN_BOND_DEFAULT_DISTANCE_CUTOFF
    return HALOGEN_BOND_DISTANCE_CUTOFF_BY_ELEMENT.get(token, HALOGEN_BOND_DEFAULT_DISTANCE_CUTOFF)


def _is_halogen_donor_atom(
    *,
    res_name: str,
    atom_name: str,
    element: str,
) -> bool:
    residue = str(res_name or "").strip().upper()
    atom = _normalize_atom_name(atom_name)
    atom_element = str(element or "").strip().upper()
    if not atom_element:
        atom_element = guess_element(atom).upper()
    if atom_element not in HALOGEN_BOND_DONOR_ELEMENTS:
        return False
    if residue in WATER_RESIDUES:
        return False
    if residue in METAL_ELEMENTS:
        return False
    return True


def _resolve_halogen_bond_context(
    *,
    res_name_a: str,
    atom_name_a: str,
    element_a: str,
    res_name_b: str,
    atom_name_b: str,
    element_b: str,
    terms: Set[str],
    distance: Optional[float],
) -> Optional[dict]:
    explicit_term = bool(terms.intersection(HALOGEN_BOND_EXPLICIT_TERMS))
    candidates: List[dict] = []
    for donor_side in ("A", "B"):
        if donor_side == "A":
            donor_res_name = res_name_a
            donor_atom_name = atom_name_a
            donor_element = element_a
            acceptor_res_name = res_name_b
            acceptor_atom_name = atom_name_b
            acceptor_element = element_b
            acceptor_side = "B"
        else:
            donor_res_name = res_name_b
            donor_atom_name = atom_name_b
            donor_element = element_b
            acceptor_res_name = res_name_a
            acceptor_atom_name = atom_name_a
            acceptor_element = element_a
            acceptor_side = "A"

        donor_valid = _is_halogen_donor_atom(
            res_name=donor_res_name,
            atom_name=donor_atom_name,
            element=donor_element,
        )
        if not donor_valid:
            continue
        acceptor_capable = _is_hbond_acceptor_capable(
            res_name=acceptor_res_name,
            atom_name=acceptor_atom_name,
            element=acceptor_element,
        )
        acceptor_valid = bool(
            acceptor_capable
            and str(acceptor_element or "").strip().upper() in HALOGEN_BOND_ACCEPTOR_ELEMENTS
        )
        distance_limit = _halogen_bond_distance_limit(donor_element)
        distance_ok = bool(
            distance is not None
            and distance >= HBOND_HEAVY_MIN_DISTANCE
            and distance <= distance_limit
        )
        candidates.append(
            {
                "explicit_term": explicit_term,
                "donor_side": donor_side,
                "acceptor_side": acceptor_side,
                "donor_element": str(donor_element or "").strip().upper(),
                "acceptor_element": str(acceptor_element or "").strip().upper(),
                "donor_atom_name": _normalize_atom_name(donor_atom_name),
                "acceptor_atom_name": _normalize_atom_name(acceptor_atom_name),
                "donor_valid": donor_valid,
                "acceptor_valid": acceptor_valid,
                "distance_limit": distance_limit,
                "distance_ok": distance_ok,
            }
        )

    if candidates:
        candidates.sort(
            key=lambda row: (
                0 if bool(row.get("acceptor_valid")) else 1,
                0 if bool(row.get("distance_ok")) else 1,
                float(row.get("distance_limit") or HALOGEN_BOND_DEFAULT_DISTANCE_CUTOFF),
            )
        )
        return candidates[0]

    if not explicit_term:
        return None
    # Explicit halogen-bond term present but endpoint typing failed.
    donor_like_side = ""
    donor_like_element = ""
    acceptor_like_side = ""
    acceptor_like_element = ""
    if str(element_a or "").strip().upper() in HALOGEN_BOND_DONOR_ELEMENTS:
        donor_like_side = "A"
        donor_like_element = str(element_a or "").strip().upper()
        acceptor_like_side = "B"
        acceptor_like_element = str(element_b or "").strip().upper()
    elif str(element_b or "").strip().upper() in HALOGEN_BOND_DONOR_ELEMENTS:
        donor_like_side = "B"
        donor_like_element = str(element_b or "").strip().upper()
        acceptor_like_side = "A"
        acceptor_like_element = str(element_a or "").strip().upper()
    distance_limit = _halogen_bond_distance_limit(donor_like_element)
    return {
        "explicit_term": True,
        "donor_side": donor_like_side,
        "acceptor_side": acceptor_like_side,
        "donor_element": donor_like_element,
        "acceptor_element": acceptor_like_element,
        "donor_atom_name": "",
        "acceptor_atom_name": "",
        "donor_valid": bool(donor_like_side),
        "acceptor_valid": bool(acceptor_like_element in HALOGEN_BOND_ACCEPTOR_ELEMENTS),
        "distance_limit": distance_limit,
        "distance_ok": bool(
            distance is not None
            and distance >= HBOND_HEAVY_MIN_DISTANCE
            and distance <= distance_limit
        ),
    }


def _vdw_radius(element: str) -> float:
    token = str(element or "").strip().upper()
    if not token:
        return DEFAULT_VDW_RADIUS
    return VDW_RADIUS_BY_ELEMENT.get(token, DEFAULT_VDW_RADIUS)


def _estimate_vdw_overlap(
    distance: Optional[float],
    element_a: str,
    element_b: str,
) -> Optional[float]:
    if distance is None or distance <= 0:
        return None
    radius_sum = _vdw_radius(element_a) + _vdw_radius(element_b)
    return radius_sum - float(distance)


def _pair_specific_min_nonbonded_distance(element_a: str, element_b: str) -> Optional[float]:
    elem_a = str(element_a or "").strip().upper()
    elem_b = str(element_b or "").strip().upper()
    if not elem_a or not elem_b:
        return None
    return MIN_NONBONDED_DISTANCE_BY_ELEMENT_PAIR.get(frozenset({elem_a, elem_b}))


def _hydrophobic_min_distance_for_pair(element_a: str, element_b: str) -> float:
    elem_a = str(element_a or "").strip().upper()
    elem_b = str(element_b or "").strip().upper()
    pair_min = HYDROPHOBIC_MIN_DISTANCE_BY_ELEMENT_PAIR.get(frozenset({elem_a, elem_b}))
    vdw_based_min = _vdw_radius(elem_a) + _vdw_radius(elem_b) - 0.5
    min_distance = max(HYDROPHOBIC_MIN_DISTANCE, vdw_based_min)
    if pair_min is not None:
        min_distance = max(min_distance, pair_min)
    return min_distance


def _clean_optional_token(value: object) -> str:
    token = str(value or "").strip()
    if token in {"", ".", "?"}:
        return ""
    return token


def _extract_altloc_family_from_node(node: dict) -> str:
    if not isinstance(node, dict):
        return ""
    for key in (
        "label_alt_id",
        "pdbx_PDB_alt_id",
        "alt_id",
        "altloc",
        "alt_loc",
    ):
        token = _clean_optional_token(node.get(key))
        if token:
            return token.upper()
    return ""


def _extract_model_id_from_node(node: dict) -> str:
    if not isinstance(node, dict):
        return ""
    for key in (
        "pdbx_PDB_model_num",
        "model_num",
        "model_id",
        "auth_model_id",
        "label_model_id",
        "model",
    ):
        token = _clean_optional_token(node.get(key))
        if token:
            return token
    return ""


def _extract_chain_id_from_node(node: dict) -> str:
    if not isinstance(node, dict):
        return ""
    token = str(node.get("auth_asym_id") or node.get("label_asym_id") or "").strip()
    if token in {".", "?"}:
        return ""
    return token


def _symmetry_token_is_non_identity(value: object) -> bool:
    token = _clean_optional_token(value)
    if not token:
        return False
    compact = token.upper().replace(" ", "")
    if compact in {"1", "1_555", "1555", "X,Y,Z"}:
        return False
    if re.fullmatch(r"\d+_\d{3}", compact):
        return compact != "1_555"
    if re.fullmatch(r"\d+", compact):
        return compact != "1"
    if compact == "X,Y,Z":
        return False
    if "X,Y,Z" in compact and ("+" in compact or "-" in compact):
        return True
    # If a symmetry token is present and not one of the identity encodings above,
    # treat it as crystal-generated contact and suppress by default.
    return True


def _contact_has_non_identity_symmetry(raw_contact: dict, node_a: dict, node_b: dict) -> bool:
    if not isinstance(raw_contact, dict):
        return False
    for key in (
        "symmetry",
        "symop",
        "symmetry_op",
        "pdbx_ptnr1_symmetry",
        "pdbx_ptnr2_symmetry",
        "ptnr1_symmetry",
        "ptnr2_symmetry",
    ):
        if _symmetry_token_is_non_identity(raw_contact.get(key)):
            return True
    for node in (node_a, node_b):
        if not isinstance(node, dict):
            continue
        for key in (
            "symmetry",
            "symop",
            "symmetry_op",
            "pdbx_ptnr_symmetry",
            "ptnr_symmetry",
        ):
            if _symmetry_token_is_non_identity(node.get(key)):
                return True
    return False


def _resolve_contact_identity_issue(
    *,
    raw_contact: dict,
    residue_a: dict,
    residue_b: dict,
    node_a: dict,
    node_b: dict,
    atom_name_a: str,
    atom_name_b: str,
    aliases: ChainAliases,
) -> Optional[dict]:
    atom_key_a = _build_atom_key_from_payload(residue_a)
    atom_key_b = _build_atom_key_from_payload(residue_b)
    same_atom_signature = bool(atom_key_a and atom_key_b and atom_key_a == atom_key_b)
    if same_atom_signature:
        return {
            "reason": "same_atom_signature",
            "evidence": ["invalid_same_atom_signature", "same_atom_normalized"],
        }

    chain_a = str(residue_a.get("chain") or "").strip()
    chain_b = str(residue_b.get("chain") or "").strip()
    seq_a = str(residue_a.get("seq") or "").strip()
    seq_b = str(residue_b.get("seq") or "").strip()
    same_residue = bool(chain_a and chain_b and seq_a and seq_b and chain_a == chain_b and seq_a == seq_b)
    same_atom_name = bool(atom_name_a and atom_name_b and atom_name_a == atom_name_b)

    altloc_a = _extract_altloc_family_from_node(node_a)
    altloc_b = _extract_altloc_family_from_node(node_b)
    if same_residue and altloc_a and altloc_b and altloc_a != altloc_b:
        return {
            "reason": "incompatible_altlocs",
            "evidence": ["invalid_altloc_incompatible", "altloc_conflict"],
        }
    if same_residue and same_atom_name and altloc_a and altloc_b and altloc_a == altloc_b:
        return {
            "reason": "same_residue_atom_altloc",
            "evidence": ["invalid_same_residue_atom_altloc", "same_atom_normalized"],
        }

    model_a = _extract_model_id_from_node(node_a)
    model_b = _extract_model_id_from_node(node_b)
    if model_a and model_b and model_a != model_b:
        return {
            "reason": "different_models",
            "evidence": ["invalid_model_mismatch", "model_conflict"],
        }

    raw_chain_a = _extract_chain_id_from_node(node_a)
    raw_chain_b = _extract_chain_id_from_node(node_b)
    if (
        raw_chain_a
        and raw_chain_b
        and raw_chain_a != raw_chain_b
        and same_atom_name
        and seq_a
        and seq_b
        and seq_a == seq_b
        and aliases.normalize(raw_chain_a) == aliases.normalize(raw_chain_b)
    ):
        return {
            "reason": "duplicate_mapping_after_alias",
            "evidence": ["invalid_duplicate_mapping_after_alias", "chain_alias_collision"],
        }

    if not ALLOW_CRYSTAL_CONTACTS and _contact_has_non_identity_symmetry(raw_contact, node_a, node_b):
        return {
            "reason": "symmetry_generated_contact",
            "evidence": ["invalid_symmetry_generated_contact"],
        }
    return None


def _resolve_impossible_contact_preclassification(
    *,
    distance: Optional[float],
    element_a: str,
    element_b: str,
    suspect_invalid_mapping: bool,
    metal_coordination_supported: bool = False,
) -> Optional[dict]:
    if distance is None:
        return None
    if metal_coordination_supported:
        return None
    vdw_overlap = _estimate_vdw_overlap(distance, element_a, element_b)
    pair_min_distance = _pair_specific_min_nonbonded_distance(element_a, element_b)
    evidence: List[str] = []
    if vdw_overlap is not None and vdw_overlap >= SOFT_CLASH_PRECLASSIFY:
        evidence.append("soft_vdw_overlap_preclassification")
    if pair_min_distance is not None:
        if distance < pair_min_distance:
            evidence.extend(
                [
                    "distance_below_pair_minimum",
                    "pair_specific_distance_impossible",
                ]
            )
            if distance < INVALID_NONBONDED_STRICT_DISTANCE or suspect_invalid_mapping:
                evidence.append("invalid_contact_preclassified")
                return {
                    "family": "invalid_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "pair_specific_minimum_distance_violation",
                }
            evidence.append("clash_preclassified_from_pair_minimum")
            return {
                "family": "clash",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "pair_specific_minimum_distance_violation",
            }

    if distance < INVALID_NONBONDED_MIN_DISTANCE:
        evidence.extend(["distance_below_absolute_minimum", "invalid_contact_preclassified"])
        return {
            "family": "invalid_contact",
            "confidence": "low",
            "evidence": evidence,
            "reason_dropped": "nonbonded_distance_below_absolute_minimum",
        }

    if vdw_overlap is not None and vdw_overlap >= HARD_CLASH_PRECLASSIFY:
        evidence.extend(
            [
                "hard_vdw_overlap_preclassified",
                "clash_preclassified_before_family_assignment",
            ]
        )
        return {
            "family": "clash",
            "confidence": "low",
            "evidence": evidence,
            "reason_dropped": "hard_vdw_overlap_preclassification",
        }
    return None


def _is_likely_covalent_nonbonded_false_positive(
    *,
    distance: Optional[float],
    element_a: str,
    atom_name_a: str,
    res_name_a: str,
    element_b: str,
    atom_name_b: str,
    res_name_b: str,
    terms: Set[str],
) -> bool:
    if distance is None or distance <= 0:
        return False
    if terms.intersection({"COVALENT", "COVALENT_BOND"}):
        metal_context = _resolve_metal_contact_context(
            res_name_a=res_name_a,
            atom_name_a=atom_name_a,
            element_a=element_a,
            res_name_b=res_name_b,
            atom_name_b=atom_name_b,
            element_b=element_b,
            terms=terms,
            distance=distance,
        )
        if metal_context and bool(metal_context.get("coordination_supported")):
            return False
        return True
    atom_a = _normalize_atom_name(atom_name_a)
    atom_b = _normalize_atom_name(atom_name_b)
    elem_a = str(element_a or "").strip().upper() or guess_element(atom_a).upper()
    elem_b = str(element_b or "").strip().upper() or guess_element(atom_b).upper()
    family_a = _nucleic_base_family(res_name_a)
    family_b = _nucleic_base_family(res_name_b)
    po_pair = {elem_a, elem_b} == {"P", "O"}
    if po_pair and distance <= COVALENT_DISTANCE_MAX_PO:
        if (
            family_a
            or family_b
            or _is_nucleotide_backbone_atom_name(atom_a)
            or _is_nucleotide_backbone_atom_name(atom_b)
            or _is_phosphate_oxygen_atom_name(atom_a)
            or _is_phosphate_oxygen_atom_name(atom_b)
        ):
            return True
    return False


def _extract_arpeggio_hbond_angle(raw_contact: dict) -> Optional[float]:
    if not isinstance(raw_contact, dict):
        return None
    for key in (
        "hbond_angle",
        "hb_angle",
        "hbAngle",
        "angle",
        "dha_angle",
        "DHA_angle",
        "angle_donor_h_acceptor",
    ):
        value = _coerce_float(raw_contact.get(key))
        if value is None:
            continue
        if value <= 0 or value > 180:
            continue
        return value
    return None


def _extract_ring_geometry_angle(raw_contact: dict) -> Optional[float]:
    if not isinstance(raw_contact, dict):
        return None
    for key in ("ring_normal_angle", "ring_angle", "pi_stack_angle"):
        value = _coerce_float(raw_contact.get(key))
        if value is None:
            continue
        if value < 0 or value > 180:
            continue
        return value
    return None


def _extract_ring_geometry_interplanar_distance(raw_contact: dict) -> Optional[float]:
    if not isinstance(raw_contact, dict):
        return None
    for key in (
        "ring_interplanar_distance",
        "interplanar_distance",
        "ring_plane_distance",
        "plane_distance",
    ):
        value = _coerce_float(raw_contact.get(key))
        if value is None:
            continue
        if value < 0:
            continue
        return value
    return None


def _extract_ring_geometry_lateral_offset(raw_contact: dict) -> Optional[float]:
    if not isinstance(raw_contact, dict):
        return None
    for key in (
        "ring_lateral_offset",
        "lateral_offset",
        "ring_offset",
    ):
        value = _coerce_float(raw_contact.get(key))
        if value is None:
            continue
        if value < 0:
            continue
        return value
    return None


def _has_ring_geometry_metrics(raw_contact: dict) -> bool:
    if not isinstance(raw_contact, dict):
        return False
    for key in (
        "ring_centroid_distance",
        "ring_normal_angle",
        "ring_angle",
        "pi_stack_angle",
        "ring_interplanar_distance",
        "interplanar_distance",
        "ring_plane_distance",
        "plane_distance",
        "ring_lateral_offset",
        "lateral_offset",
        "ring_offset",
    ):
        value = _coerce_float(raw_contact.get(key))
        if value is not None:
            return True
    return False


def _resolve_metal_contact_context(
    *,
    res_name_a: str,
    atom_name_a: str,
    element_a: str,
    res_name_b: str,
    atom_name_b: str,
    element_b: str,
    terms: Set[str],
    distance: Optional[float],
) -> Optional[dict]:
    res_name_a_clean = str(res_name_a or "").strip().upper().replace("+", "").replace("-", "")
    res_name_b_clean = str(res_name_b or "").strip().upper().replace("+", "").replace("-", "")
    metal_on_a = element_a in METAL_ELEMENTS or res_name_a_clean in METAL_ELEMENTS
    metal_on_b = element_b in METAL_ELEMENTS or res_name_b_clean in METAL_ELEMENTS
    if not (metal_on_a or metal_on_b):
        return None

    metal_side = ""
    metal_element = ""
    donor_element = ""
    donor_atom_name = ""

    def _resolved_metal_element(element: str, residue_name: str) -> str:
        element_token = str(element or "").strip().upper()
        residue_token = str(residue_name or "").strip().upper()
        if element_token in METAL_ELEMENTS:
            return element_token
        if residue_token in METAL_ELEMENTS:
            return residue_token
        return element_token or residue_token

    if metal_on_a and not metal_on_b:
        metal_side = "A"
        metal_element = _resolved_metal_element(element_a, res_name_a_clean)
        donor_element = element_b
        donor_atom_name = _primary_contact_atom_name(atom_name_b)
    elif metal_on_b and not metal_on_a:
        metal_side = "B"
        metal_element = _resolved_metal_element(element_b, res_name_b_clean)
        donor_element = element_a
        donor_atom_name = _primary_contact_atom_name(atom_name_a)
    else:
        metal_element = (
            _resolved_metal_element(element_a, res_name_a_clean)
            or _resolved_metal_element(element_b, res_name_b_clean)
        )

    explicit_term = bool(terms.intersection(METAL_CONTACT_TERMS))
    donor_valid = bool(donor_element in METAL_DONOR_ELEMENTS)
    cutoff = METAL_COORDINATION_CUTOFF.get(metal_element, METAL_DEFAULT_COORDINATION_CUTOFF)
    distance_ok = bool(distance is not None and distance <= cutoff)
    coordination_supported = bool(donor_valid and distance_ok)

    return {
        "metal_on_a": metal_on_a,
        "metal_on_b": metal_on_b,
        "metal_side": metal_side,
        "metal_element": metal_element,
        "donor_element": donor_element,
        "donor_atom_name": donor_atom_name,
        "donor_valid": donor_valid,
        "cutoff": cutoff,
        "distance_ok": distance_ok,
        "explicit_term": explicit_term,
        "coordination_supported": coordination_supported,
    }


def _build_atom_key_from_payload(residue: dict) -> str:
    if not isinstance(residue, dict):
        return ""
    chain = str(residue.get("chain") or "").strip()
    seq = str(residue.get("seq") or "").strip()
    if not chain or not seq:
        return ""
    res_name = str(residue.get("resName") or "").strip().upper() or "UNK"
    atom = _primary_contact_atom_name(residue.get("atom"))
    atom_token = atom or "?"
    return f"{chain}:{seq}:{res_name}:{atom_token}"


def _build_unordered_pair_key(left: str, right: str) -> str:
    a = str(left or "").strip()
    b = str(right or "").strip()
    if not a and not b:
        return ""
    if not a:
        return b
    if not b:
        return a
    return f"{a}|{b}" if a <= b else f"{b}|{a}"


def _base_pair_family(base_a: str, base_b: str) -> str:
    token_a = str(base_a or "").strip().upper()
    token_b = str(base_b or "").strip().upper()
    if not token_a or not token_b:
        return "other"
    return NUCLEIC_BASE_PAIR_FAMILY_BY_BASES.get(frozenset({token_a, token_b}), "other")


def _append_evidence(evidence: List[str], token: str) -> None:
    item = str(token or "").strip()
    if not item:
        return
    if item in evidence:
        return
    evidence.append(item)


def _clamp_unit_interval(value: float) -> float:
    if value <= 0.0:
        return 0.0
    if value >= 1.0:
        return 1.0
    return value


def _base_pair_support_score_component(support_count: int) -> float:
    if support_count >= BASE_PAIR_MIN_POLAR_PAIR_SUPPORT:
        return 1.0
    if support_count == 1:
        return 0.45
    return 0.0


def _base_pair_distance_score_component(distance: Optional[float]) -> float:
    if distance is None:
        return 0.0
    if distance <= BASE_PAIR_DISTANCE_STRONG_MAX:
        return 1.0
    if distance <= BASE_PAIR_DISTANCE_MEDIUM_MAX:
        return 0.82
    if distance <= BASE_PAIR_CANDIDATE_MAX_DISTANCE:
        return 0.62
    return 0.0


def _base_pair_angle_score_component(
    evaluated_count: int,
    passed_count: int,
    strong_count: int,
) -> float:
    if evaluated_count <= 0:
        return 0.0
    weighted_passes = max(0.0, float(passed_count) + 0.5 * float(strong_count))
    return _clamp_unit_interval(weighted_passes / 2.0)


def _base_pair_coplanarity_supported_from_metrics(
    ring_normal_angle: Optional[float],
    ring_interplanar_distance: Optional[float],
    ring_lateral_offset: Optional[float],
) -> bool:
    if ring_normal_angle is None:
        return False
    if ring_interplanar_distance is None:
        return False
    if ring_lateral_offset is None:
        return False
    if ring_normal_angle > BASE_PAIR_RING_PLANE_MAX_NORMAL_ANGLE:
        return False
    if ring_interplanar_distance > BASE_PAIR_RING_PLANE_MAX_INTERPLANAR_DISTANCE:
        return False
    if ring_lateral_offset < BASE_PAIR_RING_PLANE_MIN_LATERAL_OFFSET:
        return False
    if ring_lateral_offset > BASE_PAIR_RING_PLANE_MAX_LATERAL_OFFSET:
        return False
    return True


def _base_pair_score_total(components: Dict[str, float]) -> float:
    support_component = _clamp_unit_interval(float(components.get("support_count") or 0.0))
    distance_component = _clamp_unit_interval(float(components.get("best_distance") or 0.0))
    angle_component = _clamp_unit_interval(float(components.get("angles") or 0.0))
    mutual_component = _clamp_unit_interval(float(components.get("mutual_best_match") or 0.0))
    coplanarity_component = _clamp_unit_interval(float(components.get("coplanarity") or 0.0))
    total = (
        BASE_PAIR_SCORE_WEIGHT_SUPPORT * support_component
        + BASE_PAIR_SCORE_WEIGHT_DISTANCE * distance_component
        + BASE_PAIR_SCORE_WEIGHT_ANGLE * angle_component
        + BASE_PAIR_SCORE_WEIGHT_MUTUAL_BEST * mutual_component
        + BASE_PAIR_SCORE_WEIGHT_COPLANARITY * coplanarity_component
    )
    return _clamp_unit_interval(total)


def _base_pair_pair_rank_tuple(pair_stat: dict) -> Tuple[float, float, float, float, float, float, float]:
    support_count = int(pair_stat.get("supportCount") or 0)
    best_distance = _coerce_float(pair_stat.get("bestDistance"))
    angle_passed_count = int(pair_stat.get("anglePassedCount") or 0)
    angle_strong_count = int(pair_stat.get("angleStrongCount") or 0)
    coplanarity_supported = bool(pair_stat.get("coplanaritySupported"))
    canonical_match_count = int(pair_stat.get("canonicalTemplateMatches") or 0)
    return (
        float(support_count),
        _base_pair_distance_score_component(best_distance),
        -best_distance if best_distance is not None else -math.inf,
        float(angle_passed_count),
        float(angle_strong_count),
        1.0 if coplanarity_supported else 0.0,
        float(canonical_match_count),
    )


def _build_base_pair_score_payload(pair_stat: Optional[dict]) -> Optional[dict]:
    if not isinstance(pair_stat, dict):
        return None
    support_count = int(pair_stat.get("supportCount") or 0)
    if support_count <= 0:
        return None
    best_distance = _coerce_float(pair_stat.get("bestDistance"))
    angle_evaluated = int(pair_stat.get("angleEvaluatedCount") or 0)
    angle_passed = int(pair_stat.get("anglePassedCount") or 0)
    angle_strong = int(pair_stat.get("angleStrongCount") or 0)
    coplanarity_supported = bool(pair_stat.get("coplanaritySupported"))
    mutual_best_match = bool(pair_stat.get("mutualBestMatch"))
    components = dict(pair_stat.get("scoreComponents") or {})
    payload = {
        "total": round(_clamp_unit_interval(float(pair_stat.get("scoreTotal") or 0.0)), 3),
        "components": {
            "support_count": round(_clamp_unit_interval(float(components.get("support_count") or 0.0)), 3),
            "best_distance": round(_clamp_unit_interval(float(components.get("best_distance") or 0.0)), 3),
            "angles": round(_clamp_unit_interval(float(components.get("angles") or 0.0)), 3),
            "mutual_best_match": round(_clamp_unit_interval(float(components.get("mutual_best_match") or 0.0)), 3),
            "coplanarity": round(_clamp_unit_interval(float(components.get("coplanarity") or 0.0)), 3),
        },
        "supportCount": support_count,
        "angleEvaluatedCount": angle_evaluated,
        "anglePassedCount": angle_passed,
        "angleStrongCount": angle_strong,
        "mutualBestMatch": mutual_best_match,
        "coplanaritySupported": coplanarity_supported,
    }
    if best_distance is not None:
        payload["bestDistance"] = round(best_distance, 3)

    ring_normal_angle = _coerce_float(pair_stat.get("ringNormalAngle"))
    ring_interplanar_distance = _coerce_float(pair_stat.get("ringInterplanarDistance"))
    ring_lateral_offset = _coerce_float(pair_stat.get("ringLateralOffset"))
    ring_geometry = {}
    if ring_normal_angle is not None:
        ring_geometry["normalAngle"] = round(ring_normal_angle, 3)
    if ring_interplanar_distance is not None:
        ring_geometry["interplanarDistance"] = round(ring_interplanar_distance, 3)
    if ring_lateral_offset is not None:
        ring_geometry["lateralOffset"] = round(ring_lateral_offset, 3)
    if ring_geometry:
        payload["ringGeometry"] = ring_geometry

    return payload


def _classify_aromatic_context_contact(
    *,
    residue_a: dict,
    residue_b: dict,
    res_name_a: str,
    res_name_b: str,
    atom_name_a: str,
    atom_name_b: str,
    distance: Optional[float],
    terms: Set[str],
    evidence: List[str],
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Optional[dict]:
    if not residue_atoms_index:
        return None
    non_polymer_a = _is_nonpolymer_residue_name(res_name_a)
    non_polymer_b = _is_nonpolymer_residue_name(res_name_b)
    if not (non_polymer_a or non_polymer_b):
        return None
    ring_site_a = _is_contact_atom_in_residue_ring(
        residue_a,
        res_name_a,
        atom_name_a,
        residue_atoms_index,
    )
    ring_site_b = _is_contact_atom_in_residue_ring(
        residue_b,
        res_name_b,
        atom_name_b,
        residue_atoms_index,
    )
    if not (ring_site_a and ring_site_b):
        return None

    _append_evidence(evidence, "aromatic_context_possible")
    aromatic_context_ring_geometry = _compute_ring_geometry_metrics_for_contact(
        residue_a,
        residue_b,
        residue_atoms_index,
        atom_name_a=atom_name_a,
        atom_name_b=atom_name_b,
    )
    ring_angle = _coerce_float(aromatic_context_ring_geometry.get("ring_normal_angle"))
    ring_centroid_distance = _coerce_float(aromatic_context_ring_geometry.get("ring_centroid_distance"))
    ring_interplanar_distance = _coerce_float(aromatic_context_ring_geometry.get("ring_interplanar_distance"))
    ring_lateral_offset = _coerce_float(aromatic_context_ring_geometry.get("ring_lateral_offset"))
    ring_min_atom_distance = _coerce_float(aromatic_context_ring_geometry.get("ring_min_atom_distance"))
    aromatic_context_ring_payload = _build_ring_metrics_payload(
        centroid_distance=ring_centroid_distance,
        min_atom_distance=ring_min_atom_distance,
        interplanar_distance=ring_interplanar_distance,
        lateral_offset=ring_lateral_offset,
        normal_angle=ring_angle,
    )
    if (
        ring_angle is not None
        and ring_centroid_distance is not None
        and ring_interplanar_distance is not None
        and ring_lateral_offset is not None
    ):
        _append_evidence(evidence, "ring_geometry_available")
        _append_evidence(evidence, "ring_geometry_computed")
        centroid_distance_ok = bool(
            PI_PI_MIN_CENTROID_DISTANCE <= ring_centroid_distance <= PI_PI_MAX_CENTROID_DISTANCE
        )
        interplanar_distance_ok = bool(
            PI_PI_MIN_INTERPLANAR_DISTANCE <= ring_interplanar_distance <= PI_PI_MAX_INTERPLANAR_DISTANCE
        )
        lateral_offset_ok = bool(ring_lateral_offset <= PI_PI_MAX_LATERAL_OFFSET)
        if centroid_distance_ok:
            _append_evidence(evidence, "distance_ok_for_pi")
        else:
            _append_evidence(evidence, "distance_out_of_range_for_pi")
        if interplanar_distance_ok:
            _append_evidence(evidence, "interplanar_distance_ok_for_pi")
        else:
            _append_evidence(evidence, "interplanar_distance_out_of_range_for_pi")
        if lateral_offset_ok:
            _append_evidence(evidence, "lateral_offset_ok_for_pi")
        else:
            _append_evidence(evidence, "lateral_offset_out_of_range_for_pi")

        if centroid_distance_ok and interplanar_distance_ok and lateral_offset_ok:
            subtype = None
            if ring_angle <= PI_PI_STACKED_MAX_NORMAL_ANGLE:
                subtype = "pi_pi_stacked"
            elif (
                ring_angle >= PI_PI_TSHAPED_MIN_NORMAL_ANGLE
                and ring_lateral_offset <= PI_PI_TSHAPED_MAX_LATERAL_OFFSET
            ):
                subtype = "pi_pi_tshaped"
            _append_evidence(evidence, "aromatic_context_promoted_to_pi_pi")
            return {
                "family": "pi_pi",
                "subtype": subtype,
                "confidence": "medium",
                "evidence": evidence,
                "ring": aromatic_context_ring_payload,
            }
        _append_evidence(evidence, "aromatic_context_pi_geometry_not_supported")
    else:
        _append_evidence(evidence, "ring_geometry_missing")

    if ring_min_atom_distance is not None:
        _append_evidence(evidence, "ring_min_atom_distance_computed")
        if (
            ring_min_atom_distance >= AROMATIC_PACKING_MIN_DISTANCE
            and ring_min_atom_distance <= AROMATIC_PACKING_MAX_DISTANCE
        ):
            _append_evidence(evidence, "distance_ok_for_aromatic_packing")
            _append_evidence(evidence, "aromatic_reclassified_from_pi_context")
            result = {
                "family": "aromatic_packing",
                "confidence": "medium",
                "evidence": evidence,
                "reason_dropped": "aromatic_context_pi_geometry_not_supported",
                "ring": aromatic_context_ring_payload,
            }
            if distance is None or abs(ring_min_atom_distance - distance) > 1e-3:
                result["distanceOverride"] = ring_min_atom_distance
            return result
        _append_evidence(evidence, "distance_out_of_range_for_aromatic_packing")
        if ring_min_atom_distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
            _append_evidence(evidence, "distance_ok_for_aromatic_proximal")
            _append_evidence(evidence, "aromatic_reclassified_from_pi_context")
            result = {
                "family": "aromatic_proximal",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "aromatic_context_pi_geometry_not_supported",
                "debugOnly": True,
                "ring": aromatic_context_ring_payload,
            }
            if ring_centroid_distance is not None and (distance is None or abs(ring_centroid_distance - distance) > 1e-3):
                result["distanceOverride"] = ring_centroid_distance
            return result
        _append_evidence(evidence, "distance_out_of_range_for_aromatic_proximal")
        return None

    if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
        _append_evidence(evidence, "distance_ok_for_aromatic_proximal")
        _append_evidence(evidence, "aromatic_proximal_distance_fallback_contact_pair")
        result = {
            "family": "aromatic_proximal",
            "confidence": "low",
            "evidence": evidence,
            "reason_dropped": "aromatic_context_pi_geometry_not_supported",
            "debugOnly": True,
            "ring": aromatic_context_ring_payload,
        }
        return result
    _append_evidence(evidence, "distance_out_of_range_for_aromatic_proximal")
    return None


def _hydrate_contact_partner(residue, node, residue_atoms_index):
    result = dict(residue)
    if not isinstance(node, dict) or isinstance(node.get("auth_atom_id"), list) or "," in str(result.get("atom", "")):
        return result
    atom = _resolve_contact_atom_record_from_payload(result, residue_atoms_index)
    if atom is not None:
        result["element"] = atom.element
        result["altloc"] = atom.altloc
        result["modelId"] = atom.model_id
    return result


def _partner_hbond_roles(node, res_name, atom_name, element):
    types = node.get("atom_types")
    if isinstance(types, (list, tuple, set)):
        types = {str(value).strip().lower() for value in types}
        return "hbond donor" in types, "hbond acceptor" in types, True
    if res_name in {'HID', 'HSD', 'HIE', 'HSE', 'HIP', 'HSP'} and atom_name in {'ND1', 'NE2'}:
        protonated = {'ND1', 'NE2'} if res_name in {'HIP', 'HSP'} else {'ND1'} if res_name in {'HID', 'HSD'} else {'NE2'}
        return atom_name in protonated, atom_name not in protonated, True
    known = res_name in POLYMER_RESIDUES or res_name in WATER_RESIDUES
    return (
        _is_hbond_donor_capable(res_name=res_name, atom_name=atom_name, element=element),
        _is_hbond_acceptor_capable(res_name=res_name, atom_name=atom_name, element=element),
        known,
    )


def _semantic_atom(residue, node, residue_atoms_index, atom_name=None):
    """Stable coordinate identity; no roles are reconstructed by a consumer."""
    payload = dict(residue)
    if atom_name is not None:
        payload['atom'] = atom_name
    atom = _resolve_contact_atom_record_from_payload(payload, residue_atoms_index)
    name = _primary_contact_atom_name(payload.get('atom'))
    model = str(atom.model_id if atom is not None else payload.get('modelId') or node.get('model_id') or '')
    altloc = str(atom.altloc if atom is not None else payload.get('altloc') or node.get('label_alt_id') or '').strip()
    coords = [atom.x, atom.y, atom.z] if atom is not None else node.get('coordinates')
    coords = [round(float(v), 6) for v in coords] if isinstance(coords, (list, tuple)) and len(coords) == 3 and all(_coerce_float(v) is not None for v in coords) else None
    return {'id': f'{_build_atom_key_from_payload(payload)}@model={model or "?"}@alt={altloc or "."}',
            'atomName': name, 'element': str(atom.element if atom is not None else payload.get('element') or node.get('type_symbol') or '').upper(),
            'coordinates': coords, 'modelId': model or None, 'altloc': altloc}


def _semantic_node(residue, node, residue_atoms_index):
    name = _primary_contact_atom_name(residue.get('atom'))
    if name == _primary_contact_atom_name(node.get('auth_atom_id') or node.get('label_atom_id')) and any(k in node for k in ('atom_types','formal_charge','bonded_atoms')):
        return node
    cache = getattr(residue_atoms_index, 'chemical_nodes', {})
    cached = cache.get((residue.get('chain'), residue.get('seq'), name))
    if cached is not None:
        return cached
    if name == _primary_contact_atom_name(node.get('auth_atom_id') or node.get('label_atom_id')):
        return node
    return {}


def _semantic_charge(node, residue):
    formal = _coerce_float(node.get('formal_charge'))
    if formal is not None:
        return {'value': int(formal) if formal.is_integer() else formal,
                'source': str(node.get('formal_charge_source') or 'arpeggio_openbabel_formal_charge'),
                'inferred': node.get('formal_charge_source') != 'input_formal_charge'}
    types = set(node.get('atom_types') or ())
    pos, neg = 'pos ionisable' in types, 'neg ionisable' in types
    if pos != neg:
        return {'value': None, 'sign': 1 if pos else -1, 'source': 'arpeggio_ionisable_type', 'inferred': True}
    args = (residue.get('resName'), residue.get('atom'), residue.get('element'))
    pos, neg = _is_salt_bridge_cation_site(*args), _is_salt_bridge_anion_site(*args)
    if pos != neg and str(residue.get('element') or '').upper() not in METAL_ELEMENTS:
        return {'value': None, 'sign': 1 if pos else -1, 'source': 'standard_residue_charge_template', 'inferred': True}
    return {'value': None, 'source': 'unknown', 'inferred': True}


def _semantic_site_charge(node, residue, residue_atoms_index):
    charge = _semantic_charge(node, residue)
    template_sign = 1 if _is_salt_bridge_cation_site(residue.get('resName'),residue.get('atom'),residue.get('element')) else -1 if _is_salt_bridge_anion_site(residue.get('resName'),residue.get('atom'),residue.get('element')) else 0
    if template_sign:
        atoms = _collect_salt_bridge_site_atoms(residue=residue,res_name=residue.get('resName'),site_kind='cation' if template_sign > 0 else 'anion',residue_atoms_index=residue_atoms_index)
        charges = []
        for atom in atoms:
            payload = dict(residue, atom=atom.atom_name, element=atom.element)
            charges.append(_semantic_charge(_semantic_node(payload,node,residue_atoms_index),payload))
        if len(charges) > 1 and all(q['value'] is not None for q in charges):
            charge = {'value': sum(q['value'] for q in charges), 'source': 'formal_charge_sum_of_site_atoms',
                      'inferred': any(q['inferred'] for q in charges), 'atomValue': charge.get('value'), 'scope': 'group'}
    return charge


def _semantic_charge_sign(charge):
    value = charge.get('value')
    if value is not None:
        return 1 if value > 0 else -1 if value < 0 else 0
    return int(charge.get('sign') or 0)


def _semantic_hydrogen_evidence(raw, residues, nodes, roles, residue_atoms_index):
    candidates = []
    exported = raw.get('hbond_hydrogen')
    exported_side = str(raw.get('hbond_donor_side') or '').upper()
    if isinstance(exported, dict) and exported_side in ('A', 'B'):
        donor_idx = 0 if exported_side == 'A' else 1
        donor = _semantic_atom(residues[donor_idx], nodes[donor_idx], residue_atoms_index)
        h = dict(exported)
        h.setdefault('element', 'H')
        h.setdefault('modelId', donor['modelId'])
        h.setdefault('altloc', donor['altloc'])
        h.setdefault('atomName', None)
        source = str(raw.get('hbond_geometry_source') or h.get('source') or 'unknown')
        h['source'], h['inferred'] = source, source != 'input_hydrogens'
        if h.get('atomName'):
            h['id'] = _semantic_atom(residues[donor_idx], h, residue_atoms_index, h['atomName'])['id']
        else:
            xyz = ','.join(str(round(float(v), 6)) for v in h.get('coordinates') or [])
            h['id'] = donor['id'] + ':H@' + xyz
        acceptor = _semantic_atom(residues[1-donor_idx], nodes[1-donor_idx], residue_atoms_index)
        if h.get('coordinates') and donor.get('coordinates') and acceptor.get('coordinates'):
            dh = math.dist(donor['coordinates'], h['coordinates'])
            if 0.45 <= dh <= (1.65 if donor['element'] in {'S','SE'} else 1.3):
                angle = _angle_between_vectors_degrees(tuple(donor['coordinates'][i]-h['coordinates'][i] for i in range(3)),tuple(acceptor['coordinates'][i]-h['coordinates'][i] for i in range(3)))
                candidates.append((angle,donor_idx,h,math.dist(h['coordinates'],acceptor['coordinates'])))
    # Deposited H coordinates can provide evidence even for older engine exports.
    for donor_idx, acceptor_idx in roles:
        donor = _resolve_contact_atom_record_from_payload(residues[donor_idx], residue_atoms_index)
        acceptor = _resolve_contact_atom_record_from_payload(residues[acceptor_idx], residue_atoms_index)
        if donor is None or acceptor is None:
            continue
        for hydrogen in _residue_atoms_from_payload(residues[donor_idx], residue_atoms_index):
            if hydrogen.element not in {'H', 'D'}:
                continue
            # An H must belong to this donor, not just be any nearby residue H.
            neighbors = [at for at in _residue_atoms_from_payload(residues[donor_idx], residue_atoms_index) if at.element not in {'H', 'D'}]
            nearest = min(neighbors, key=lambda at: _distance_sq(at, hydrogen), default=None)
            if nearest is not donor or not 0.45 <= distance(donor, hydrogen) <= (1.65 if donor.element in {'S', 'SE'} else 1.3):
                continue
            angle = _angle_between_vectors_degrees((donor.x-hydrogen.x, donor.y-hydrogen.y, donor.z-hydrogen.z),
                                                    (acceptor.x-hydrogen.x, acceptor.y-hydrogen.y, acceptor.z-hydrogen.z))
            h = _semantic_atom(residues[donor_idx], {}, residue_atoms_index, hydrogen.atom_name)
            h.update(source='input_hydrogens', inferred=False, bondSource='coordinate_bond_inference')
            candidates.append((angle, donor_idx, h, distance(hydrogen, acceptor)))
    candidates = [row for row in candidates if _coerce_float(row[0]) is not None]
    return max(candidates, key=lambda row: (not row[2].get('inferred'), row[0], row[2]['id']), default=None)


def _semantic_site(residue, node, residue_atoms_index, kind='atom', names=None, site_id=None):
    names = list(names) if names else [_primary_contact_atom_name(residue.get('atom'))]
    atoms = [_semantic_atom(residue, node if len(names) == 1 else {}, residue_atoms_index, name) for name in sorted(set(names)) if name]
    identity = site_id or '|'.join(atom['id'] for atom in atoms)
    return {'kind': kind, 'id': identity, 'residue': {key: residue.get(key) for key in ('chain', 'seq', 'resName')}, 'atoms': atoms}


def _semantic_ring_site(residue, other, residue_atoms_index, names=None):
    descriptors = _residue_ring_descriptors(residue, residue_atoms_index)
    hints = set(names or str(residue.get('atom') or '').split(','))
    partner = _resolve_contact_atom_record_from_payload(other, residue_atoms_index)
    def rank(desc):
        center = desc.get('centroid') or (0, 0, 0)
        d2 = sum((center[i] - (partner.x, partner.y, partner.z)[i])**2 for i in range(3)) if partner else 0
        return len(hints - set(desc.get('atom_names') or [])), d2, desc.get('hash', '')
    descriptor = min(descriptors, key=rank, default=None)
    if descriptor is None:
        return _semantic_site(residue, {}, residue_atoms_index, 'ring', names=names or list(hints))
    site = _semantic_site(residue, {}, residue_atoms_index, 'ring', names=descriptor['atom_names'])
    scope = site['atoms'][0] if site['atoms'] else {}
    site['id'] = _ring_site_key_from_descriptor(residue, descriptor) + f'@model={scope.get("modelId") or "?"}@alt={scope.get("altloc") or "."}'
    site['centroid'] = [round(float(v), 6) for v in descriptor['centroid']]
    normal = _best_plane_normal(descriptor['points'])
    if normal:
        site['normal'] = [round(float(v), 6) for v in normal]
    site['provenance'] = 'standard_residue_ring_template' if _residue_ring_atom_names(residue.get('resName')) else 'coordinate_ring_inference'
    return site


def _semantic_identity(semantics):
    # Canonical reversal-invariant identity. Different features/H atoms remain distinct.
    sites = sorted(f'{p.get("role", "unresolved")}:{p.get("site", {}).get("id", "?")}' for p in semantics.get('participants', []))
    hydrogen = semantics.get('geometry', {}).get('hydrogen') or {}
    return 'sem1:' + str(semantics.get('family') or 'other') + ':' + '||'.join(sites) + (':H=' + str(hydrogen['id']) if hydrogen.get('id') else '')


def _build_interaction_semantics(raw, asserted, residue_a, residue_b, residue_atoms_index):
    """One chemical record for the API, explanations, labels and rendering.

    Geometry support is not experimental confirmation, and perceived atom types or
    formal charges are not measurements of the input structure's protonation state.
    """
    family = str(asserted.get('family') or 'other')
    residues = [dict(residue_a), dict(residue_b)]
    nodes = [raw.get('bgn') or {}, raw.get('end') or {}]
    for i, side in enumerate(('A', 'B')):
        if asserted.get('atomOverride' + side):
            residues[i]['atom'] = asserted['atomOverride' + side]
        if asserted.get('elementOverride' + side):
            residues[i]['element'] = asserted['elementOverride' + side]
        nodes[i] = _semantic_node(residues[i], nodes[i], residue_atoms_index)
    profiles = [_partner_hbond_roles(node, res.get('resName', ''), res.get('atom', ''), res.get('element', '')) for node, res in zip(nodes, residues)]
    participants = []
    flags = []
    for side, residue, node in zip(('A', 'B'), residues, nodes):
        typed = isinstance(node.get('atom_types'), (list, tuple, set))
        provenance = ['arpeggio_openbabel_atom_types'] if typed else ['standard_residue_template'] if residue.get('resName') in POLYMER_RESIDUES | WATER_RESIDUES else []
        if residue.get('resName') in {'HID','HSD','HIE','HSE','HIP','HSP'}:
            provenance.append('explicit_protonation_label')
        participants.append({'side': side, 'role': 'contact_atom', 'site': _semantic_site(residue, node, residue_atoms_index),
                             'roleProvenance': provenance, 'atomTypes': sorted(node.get('atom_types') or [])})
    value, distance_source = _resolve_contact_distance_value(raw, *residues, residue_atoms_index)
    override = _coerce_float(asserted.get('distanceOverride'))
    if override is not None:
        value, distance_source = override, 'asserted_site_coordinates'
    geometry = {'distance': {'value': _round_metric_or_none(value), 'kind': 'atom_pair', 'unit': 'angstrom',
                            'source': 'model_coordinates' if distance_source == 'distance_recomputed_from_coordinates' else 'engine_geometry'}, 'measurements': []}
    direction = {'from': None, 'to': None, 'certainty': 'not_applicable'}
    directionality = 'symmetric'
    chemical = 'unknown'
    geometric = 'partial' if value is not None else 'missing'
    state = 'modelled' if any('arpeggio_openbabel_atom_types' in p['roleProvenance'] for p in participants) else 'assumed'
    level = 'candidate'
    def measure(kind, number, unit='angstrom', source='model_coordinates'):
        number = _round_metric_or_none(number)
        if number is not None:
            geometry['measurements'].append({'kind': kind, 'value': number, 'unit': unit, 'source': source})
    roles = [(0, 1)] if profiles[0][0] and profiles[1][1] else []
    if profiles[1][0] and profiles[0][1]:
        roles.append((1, 0))
    if family == 'hbond':
        directionality = 'directional'
        geometry['distance']['kind'] = 'donor_acceptor'
        chemical = 'supported' if all(isinstance(n.get('atom_types'), (list, tuple, set)) for n in nodes) and roles else 'inferred' if roles else 'unknown'
        h_evidence = _semantic_hydrogen_evidence(raw, residues, nodes, roles, residue_atoms_index)
        selected = roles[0] if len(roles) == 1 else None
        if h_evidence and (h_evidence[1], 1-h_evidence[1]) in roles:
            angle, donor_idx, hydrogen, h_distance = h_evidence
            geometry['hydrogen'] = hydrogen
            measure('donor_hydrogen_acceptor_angle', angle, 'degree', hydrogen['source'])
            measure('hydrogen_acceptor', h_distance, source=hydrogen['source'])
            max_h_distance = 1.2 + _vdw_radius(residues[1-donor_idx].get('element','')) + ARPEGGIO_VDW_COMP
            if angle >= HBOND_STRONG_ANGLE_MIN and h_distance is not None and 0.7 <= h_distance <= max_h_distance:
                selected = (donor_idx, 1-donor_idx)
                if not hydrogen['inferred'] and _hbond_distance_is_within_limits(value, residues[0].get('element',''), residues[1].get('element','')):
                    geometric, level = 'supported', 'geometrically_supported'
                elif hydrogen['inferred']:
                    flags.append('hydrogen_position_modelled')
                else:
                    flags.append('hydrogen_donor_acceptor_distance_outside_criteria')
            else:
                if angle < HBOND_STRONG_ANGLE_MIN:
                    flags.append('hydrogen_angle_not_supportive')
                if h_distance is None or not 0.7 <= h_distance <= max_h_distance:
                    flags.append('hydrogen_acceptor_distance_outside_criteria')
        else:
            flags.append('hydrogen_missing')
            # A legacy angle with no associated H/provenance is an unverified metric.
            if _extract_arpeggio_hbond_angle(raw) is not None:
                measure('donor_hydrogen_acceptor_angle', _extract_arpeggio_hbond_angle(raw), 'degree', str(raw.get('hbond_geometry_source') or 'unverified_engine_angle'))
                flags.append('hydrogen_angle_provenance_incomplete')
        for i, (residue, node) in enumerate(zip(residues, nodes)):
            if residue.get('resName') == 'HIS' and residue.get('atom') in {'ND1','NE2'} and not (geometry.get('hydrogen') and selected and selected[0] == i and not geometry['hydrogen']['inferred']):
                flags.append('histidine_protonation_or_tautomer_uncertain')
        if selected:
            donor_idx, acceptor_idx = selected
            participants[donor_idx]['role'], participants[acceptor_idx]['role'] = 'donor', 'acceptor'
            certainty = 'certain' if geometry.get('hydrogen') and not geometry['hydrogen']['inferred'] and geometric == 'supported' else 'inferred'
            if 'histidine_protonation_or_tautomer_uncertain' in flags:
                certainty, level = 'ambiguous', 'ambiguous'
            direction.update({'from': ('A','B')[donor_idx], 'to': ('A','B')[acceptor_idx], 'certainty': certainty})
        else:
            for p in participants:
                p['role'] = 'donor_or_acceptor' if roles else 'unresolved'
            direction['certainty'] = 'ambiguous'
            direction['alternatives'] = [{'from': ('A','B')[d], 'to': ('A','B')[a]} for d,a in roles]
            flags.append('donor_acceptor_direction_ambiguous')
            level = 'ambiguous'
    elif family == 'polar_contact':
        # Preserve individual chemical roles even when they do not form a valid
        # H-bond pair; no directional bond is invented for incompatible sites.
        for participant, (donor, acceptor, known) in zip(participants, profiles):
            participant['role'] = 'donor_or_acceptor' if donor and acceptor else 'donor' if donor else 'acceptor' if acceptor else 'unresolved'
        chemical = 'inferred' if roles else 'unknown'
        if asserted.get('basePair'):
            direction['certainty'] = 'ambiguous'
            level = 'ambiguous'
            flags.append('base_pair_atom_hbond_assignment_unresolved')
    elif family == 'halogen_bond':
        directionality = 'directional'
        geometry['distance']['kind'] = 'donor_acceptor'
        donor_idx = next((i for i,r in enumerate(residues) if r.get('element') in HALOGEN_BOND_DONOR_ELEMENTS), None)
        if donor_idx is not None:
            acceptor_idx = 1-donor_idx
            participants[donor_idx]['role'], participants[acceptor_idx]['role'] = 'halogen_donor', 'acceptor'
            direction.update({'from': ('A','B')[donor_idx], 'to': ('A','B')[acceptor_idx], 'certainty':'inferred'})
            node = nodes[donor_idx]
            anchors = [n for n in node.get('bonded_atoms') or [] if n.get('type_symbol') == 'C']
            anchor = None
            if anchors:
                anchor = _semantic_atom(residues[donor_idx], anchors[0], residue_atoms_index, anchors[0].get('auth_atom_id'))
                anchor['source'] = 'arpeggio_bond_graph'
            else:
                at = _resolve_halogen_donor_anchor_atom(residues[donor_idx], residues[donor_idx].get('atom'), residue_atoms_index)
                if at is not None:
                    anchor = _semantic_atom(residues[donor_idx], {}, residue_atoms_index, at.atom_name)
                    anchor['source'] = 'coordinate_bond_inference'
                    flags.append('halogen_covalent_anchor_inferred')
            if anchor:
                geometry['donorAnchor'] = anchor
            angle = _coerce_float(asserted.get('halogenAngle'))
            if angle is None:
                angle = _compute_halogen_bond_angle(donor_residue=residues[donor_idx], donor_atom_name=residues[donor_idx].get('atom'), acceptor_residue=residues[acceptor_idx], acceptor_atom_name=residues[acceptor_idx].get('atom'), residue_atoms_index=residue_atoms_index)
            measure('donor_anchor_halogen_acceptor_angle', angle, 'degree', anchor['source'] if anchor else 'model_coordinates')
            compatible = profiles[acceptor_idx][1] or 'xbond acceptor' in nodes[acceptor_idx].get('atom_types', [])
            chemical = 'supported' if compatible and anchors else 'inferred' if compatible else 'unknown'
            if anchor and angle is not None and angle >= HALOGEN_BOND_MEDIUM_ANGLE_MIN and chemical != 'unknown':
                geometric, level = 'supported', 'geometrically_supported'
            else:
                flags.append('halogen_geometry_or_typing_incomplete')
        else:
            direction['certainty'] = 'ambiguous'
            flags.append('halogen_donor_unresolved')
    elif family == 'metal_coordination':
        directionality = 'role_specific'
        geometry['distance']['kind'] = 'metal_donor'
        metal_idx = next((i for i,r in enumerate(residues) if r.get('element') in METAL_ELEMENTS), None)
        if metal_idx is not None:
            donor_idx = 1-metal_idx
            participants[metal_idx]['role'], participants[metal_idx]['site']['kind'] = 'metal_center', 'metal'
            participants[metal_idx]['roleProvenance'] = ['element_identity']
            participants[metal_idx]['charge'] = _semantic_charge(nodes[metal_idx], residues[metal_idx])
            participants[donor_idx]['role'] = 'coordinating_atom'
            participants[metal_idx]['site']['coordinationSiteId'] = participants[metal_idx]['site']['id']
            chemical = 'supported' if profiles[donor_idx][1] and isinstance(nodes[donor_idx].get('atom_types'), (list,tuple,set)) else 'inferred' if profiles[donor_idx][1] else 'unknown'
            geometric = 'supported' if value is not None and value <= METAL_COORDINATION_CUTOFF.get(residues[metal_idx].get('element'), METAL_DEFAULT_COORDINATION_CUTOFF) else geometric
            level = 'geometrically_supported' if geometric == 'supported' and chemical != 'unknown' else 'candidate'
            if participants[metal_idx]['charge']['value'] is None:
                flags.append('metal_oxidation_state_unknown')
            if chemical == 'unknown':
                flags.append('coordination_donor_chemistry_unresolved')
    elif family in {'salt_bridge', 'pi_cation'}:
        directionality = 'role_specific'
        charges = [_semantic_site_charge(n,r,residue_atoms_index) for n,r in zip(nodes,residues)]
        signs = [_semantic_charge_sign(q) for q in charges]
        for i,q in enumerate(charges):
            participants[i]['charge'] = q
        if family == 'salt_bridge':
            geometry['distance']['kind'] = 'charge_site'
            for i, residue in enumerate(residues):
                sign = signs[i]
                if not sign and charges[i]['value'] == 0:
                    flags.append('explicit_neutral_site')
                participants[i]['role'] = 'positive_site' if sign > 0 else 'negative_site' if sign < 0 else 'unresolved'
                group_atoms = _collect_salt_bridge_site_atoms(residue=residue,res_name=residue.get('resName'),site_kind='cation' if sign > 0 else 'anion',residue_atoms_index=residue_atoms_index) if sign else []
                if len(group_atoms) > 1:
                    participants[i]['site'] = _semantic_site(residue,nodes[i],residue_atoms_index,'group',names=[a.atom_name for a in group_atoms])
                    participants[i]['site']['contactAtom'] = _semantic_atom(residue,nodes[i],residue_atoms_index)
            chemical = 'inferred' if {p['role'] for p in participants} == {'positive_site','negative_site'} else 'unknown'
            if chemical == 'inferred' and all(q['value'] is not None and q['value'] != 0 for q in charges):
                chemical = 'supported'
            if value is not None and value <= SALT_BRIDGE_MAX_DISTANCE and chemical != 'unknown':
                geometric = 'supported'
                level = 'chemically_supported' if chemical == 'supported' else 'candidate'
            flags.append('charge_state_not_experimentally_determined')
        else:
            cation_idx = next((i for i,sign in enumerate(signs) if sign > 0), None)
            if cation_idx is not None:
                ring_idx = 1-cation_idx
                participants[cation_idx]['role'] = 'cation'
                participants[ring_idx]['role'] = 'aromatic_ring'
                ring_names = _resolve_aromatic_ring_site_keys(*residues,residue_atoms_index).get('ringAtomNames' + ('A','B')[ring_idx])
                participants[ring_idx]['site'] = _semantic_ring_site(residues[ring_idx],residues[cation_idx],residue_atoms_index,ring_names)
                participants[ring_idx]['roleProvenance'] = ['ring_membership']
                center = participants[ring_idx]['site'].get('centroid')
                xyz = participants[cation_idx]['site']['atoms'][0].get('coordinates')
                if center and xyz:
                    vector = [xyz[i]-center[i] for i in range(3)]
                    value = math.sqrt(sum(v*v for v in vector))
                    geometry['distance'].update(value=round(value,3),kind='cation_centroid',source='model_coordinates')
                    normal = participants[ring_idx]['site'].get('normal')
                    if normal:
                        offset = abs(sum(vector[i]*normal[i] for i in range(3)))
                        measure('cation_plane_offset',offset)
                        measure('cation_lateral_offset',math.sqrt(max(0,value*value-offset*offset)))
                        measure('cation_ring_normal_angle',math.degrees(math.acos(max(0,min(1,offset/value)))) if value else None,'degree')
                else:
                    geometry['distance'].update(value=None,kind='cation_centroid')
                    flags.append('ring_coordinates_unavailable')
                chemical = 'supported' if signs[cation_idx] > 0 and charges[cation_idx]['value'] is not None else 'inferred'
                lateral = next((m['value'] for m in geometry['measurements'] if m['kind']=='cation_lateral_offset'), None)
                theta = next((m['value'] for m in geometry['measurements'] if m['kind']=='cation_ring_normal_angle'), None)
                # Match the upstream atom-aromatic 4.5 A / 30 degree normal-angle criteria.
                supported = bool(center and xyz and normal and PI_PI_MIN_CENTROID_DISTANCE <= value <= 4.5 and theta is not None and theta <= 30.0)
                geometric = 'supported' if supported else 'partial' if center and xyz else 'missing'
                level = 'geometrically_supported' if supported else 'candidate'
                if not supported:
                    flags.append('cation_pi_geometry_incomplete_or_outside_criteria')
            else:
                flags.append('cation_site_unresolved')
                geometry['distance'].update(value=None,kind='cation_centroid')
                level = 'ambiguous'
    elif family in {'pi_pi','aromatic_packing','aromatic_proximal'}:
        keys = _resolve_aromatic_ring_site_keys(*residues,residue_atoms_index)
        for i, side in enumerate(('A','B')):
            participants[i]['role'] = 'aromatic_ring'
            participants[i]['site'] = _semantic_ring_site(residues[i],residues[1-i],residue_atoms_index,keys.get('ringAtomNames'+side))
            participants[i]['roleProvenance'] = ['ring_membership']
        sa,sb = [p['site'] for p in participants]
        ca,cb = sa.get('centroid'),sb.get('centroid')
        if ca and cb:
            vector = [cb[i]-ca[i] for i in range(3)]
            centroid_distance = math.sqrt(sum(v*v for v in vector))
            measure('ring_centroid',centroid_distance)
            points_a, points_b = [a['coordinates'] for a in sa['atoms'] if a['coordinates']], [a['coordinates'] for a in sb['atoms'] if a['coordinates']]
            closest = min((math.dist(a,b) for a in points_a for b in points_b),default=None)
            measure('closest_atom',closest)
            normals = sa.get('normal'),sb.get('normal')
            if all(normals):
                angle = math.degrees(math.acos(max(0,min(1,abs(sum(normals[0][i]*normals[1][i] for i in range(3)))))))
                measure('ring_normal_angle',angle,'degree')
                for i,normal in enumerate(normals):
                    height=abs(sum(vector[j]*normal[j] for j in range(3)))
                    measure('ring_plane_separation_'+('A','B')[i],height)
                    measure('ring_lateral_offset_'+('A','B')[i],math.sqrt(max(0,centroid_distance**2-height**2)))
            geometry['distance'].update(value=_round_metric_or_none(closest if family=='aromatic_packing' else centroid_distance),kind='closest_atom' if family=='aromatic_packing' else 'ring_centroid',source='model_coordinates')
            # Support requires the family criteria, not merely computable centroids.
            ring_geometry = _compute_ring_metrics_from_point_sets(sa['atoms'] and points_a, sb['atoms'] and points_b)
            plane = _coerce_float(ring_geometry.get('ring_interplanar_distance'))
            lateral = _coerce_float(ring_geometry.get('ring_lateral_offset'))
            normal_angle = _coerce_float(ring_geometry.get('ring_normal_angle'))
            supported = (AROMATIC_PACKING_MIN_DISTANCE <= closest <= AROMATIC_PACKING_MAX_DISTANCE) if family == 'aromatic_packing' and closest is not None else bool(family == 'pi_pi' and PI_PI_MIN_CENTROID_DISTANCE <= centroid_distance <= PI_PI_MAX_CENTROID_DISTANCE and plane is not None and PI_PI_MIN_INTERPLANAR_DISTANCE <= plane <= PI_PI_MAX_INTERPLANAR_DISTANCE and lateral is not None and normal_angle is not None and ((normal_angle <= PI_PI_STACKED_MAX_NORMAL_ANGLE and lateral <= PI_PI_MAX_LATERAL_OFFSET) or (normal_angle >= PI_PI_TSHAPED_MIN_NORMAL_ANGLE and lateral <= PI_PI_TSHAPED_MAX_LATERAL_OFFSET)))
            geometric='supported' if supported else 'partial'
            if not supported:
                flags.append('aromatic_geometry_incomplete_or_outside_family_criteria')
            chemical='inferred' if 'coordinate_ring_inference' in {sa.get('provenance'),sb.get('provenance')} else 'supported'
            level='geometrically_supported' if supported else 'candidate'
        else:
            geometry['distance'].update(value=None,kind='ring_centroid')
            flags.append('ring_coordinates_unavailable')
    elif family == 'clash':
        geometric='supported' if value is not None else 'missing'
        radii=[]
        for node,residue in zip(nodes,residues):
            radii.append(_coerce_float(node.get('vdw_radius')) or _vdw_radius(residue.get('element','')))
        measure('vdw_radius_sum',sum(radii),source='arpeggio_openbabel_radii' if all(n.get('vdw_radius') is not None for n in nodes) else 'roami_vdw_reference')
        measure('vdw_overlap',sum(radii)-value if value is not None else None,source='vdw_reference_minus_separation')
        level='geometrically_supported' if geometric=='supported' else 'candidate'
    elif family in {'hydrophobic','packing_contact','vdw'}:
        chemical='supported' if all('hydrophobe' in n.get('atom_types',[]) for n in nodes) else 'inferred'
    elif family == 'base_pairing':
        for p in participants:
            p['role']='base'
        flags.append('base_pair_context_not_atom_bond_direction')
    # Engine atom–ring observations retain the ring even when their family is
    # generic/packing. The source distance is never an arbitrarily chosen atom pair.
    if str(raw.get('type') or '').lower() in {'atom-plane','plane-atom'} and family not in {'pi_cation','pi_pi','aromatic_packing','aromatic_proximal'}:
        ring_side = 1 if str(raw.get('type')).lower() == 'atom-plane' else 0
        raw_ring = raw.get(('bgn','end')[ring_side]) or {}
        names_raw = raw_ring.get('auth_atom_id') or raw_ring.get('label_atom_id') or ''
        names = [str(name).strip() for name in (names_raw if isinstance(names_raw,(list,tuple)) else str(names_raw).split(',')) if str(name).strip()]
        descriptors = _residue_ring_descriptors(residues[ring_side],residue_atoms_index)
        if names and any(set(names)==set(desc['atom_names']) for desc in descriptors):
            participants[ring_side]['site'] = _semantic_ring_site(residues[ring_side],residues[1-ring_side],residue_atoms_index,names)
        elif len(names)>=3:
            site = _semantic_site(residues[ring_side],{},residue_atoms_index,kind='ring',names=names)
            site['provenance'] = 'engine_ring_membership'
            points = [atom['coordinates'] for atom in site['atoms'] if atom.get('coordinates')]
            if len(points)==len(names):
                site['centroid'] = [round(sum(point[i] for point in points)/len(points),6) for i in range(3)]
                normal = _best_plane_normal(points)
                if normal:
                    site['normal'] = list(normal)
            participants[ring_side]['site'] = site
        else:
            site = participants[ring_side]['site']
            participants[ring_side]['site'] = {'kind':'ring','id':site['id']+':ring_membership_unknown','residue':site['residue'],'atoms':[],'provenance':'engine_ring_membership_unavailable'}
            flags.append('ring_membership_unavailable')
        participants[ring_side]['role'] = 'aromatic_ring'
        participants[ring_side]['roleProvenance'] = ['engine_atom_ring_observation','ring_membership']
        center = participants[ring_side]['site'].get('centroid')
        point = participants[1-ring_side]['site']['atoms'][0].get('coordinates')
        geometry['distance'].update(kind='atom_ring_centroid',value=_round_metric_or_none(math.dist(point,center)) if point and center else None,source='model_coordinates')
        flags.append('atom_support_not_exported')
    if asserted.get('basePair') and family in {'hbond','polar_contact'}:
        conflicts = False
        def role_name(donor, acceptor):
            return 'donor_or_acceptor' if donor and acceptor else 'donor' if donor else 'acceptor' if acceptor else 'unresolved'
        for i, (node, residue) in enumerate(zip(nodes, residues)):
            if not isinstance(node.get('atom_types'), (list,tuple,set)) or not _nucleic_base_family(residue.get('resName')):
                continue
            expected = _partner_hbond_roles({},residue.get('resName',''),residue.get('atom',''),residue.get('element',''))
            if tuple(profiles[i][:2]) != tuple(expected[:2]):
                conflicts = True
                participants[i]['roleAlternatives'] = [
                    {'role':role_name(*profiles[i][:2]),'source':'arpeggio_openbabel_atom_types'},
                    {'role':role_name(*expected[:2]),'source':'standard_residue_template'}]
                participants[i]['role'] = 'unresolved'
                participants[i]['roleProvenance'] = sorted(set(participants[i]['roleProvenance']) | {'standard_residue_template'})
        if conflicts:
            flags.append('nucleobase_typing_template_conflict')
            level, chemical = 'ambiguous', 'unknown'
            if direction.get('from') and direction.get('to'):
                direction['alternatives'] = [{'from':direction['from'],'to':direction['to'],'source':'arpeggio_openbabel_atom_types'}]
            direction.update({'from':None,'to':None,'certainty':'ambiguous'})
    if family in {'salt_bridge','pi_cation'}:
        required_charges = [p.get('charge', {}) for p in participants if p['role'] in {'positive_site','negative_site','cation'}]
        if required_charges and all(q.get('value') is not None for q in required_charges):
            state = 'modelled' if any(q.get('inferred', True) for q in required_charges) else 'input'
        elif required_charges:
            state = 'assumed'
        else:
            state = 'unknown'
    if chemical=='unknown':
        flags.append('chemical_compatibility_unresolved')
    if state=='assumed':
        flags.append('chemical_state_assumed')
    semantics={'version':1,'family':family,'directionality':directionality,'participants':participants,'direction':direction,
               'geometry':geometry,'evidence':{'level':level,'chemicalCompatibility':chemical,'geometrySupport':geometric,
               'chemicalState':state,'ambiguityFlags':sorted(set(flags)),'biologicalInterpretation':'not_evaluated'}}
    if any(asserted.get('atomOverride' + side) and asserted['atomOverride' + side] != residue.get('atom') for side,residue in zip(('A','B'),(residue_a,residue_b))):
        semantics['sourceEndpoints'] = [{'side':side,'atom':_semantic_atom(residue,node,residue_atoms_index),'atomTypes':sorted(node.get('atom_types') or []),'charge':_semantic_charge(node,residue)} for side,residue,node in zip(('A','B'),(residue_a,residue_b),(raw.get('bgn') or {},raw.get('end') or {}))]
    semantics['identity']=_semantic_identity(semantics)
    return semantics


def _merge_semantic_contact_records(first, second):
    ranks = {'candidate': 0, 'ambiguous': 0, 'chemically_supported': 1, 'geometrically_supported': 2}
    def rank(record):
        sem = record['semantics']
        return (ranks.get(sem['evidence']['level'], 0), sem['direction']['certainty'] == 'certain', -(record.get('distance') or math.inf))
    best, other = (second, first) if rank(second) > rank(first) else (first, second)
    best = copy.deepcopy(best)
    for key in ('basePair', 'canonicalBasePair'):
        if key not in best and key in other:
            best[key] = copy.deepcopy(other[key])
    for section, field in (('asserted', 'evidence'), ('arpeggio', 'terms')):
        best[section][field] = sorted(set(best.get(section, {}).get(field, [])) | set(other.get(section, {}).get(field, [])))
    evidence = best['semantics']['evidence']
    evidence['ambiguityFlags'] = sorted(set(evidence['ambiguityFlags']) | set(other['semantics']['evidence']['ambiguityFlags']))
    def observation(record):
        return copy.deepcopy({'atomPair': [record.get('atomKeyA'), record.get('atomKeyB')], 'source': record.get('source'), 'type': record.get('arpeggio',{}).get('type'),
                'terms': record.get('arpeggio', {}).get('terms', []), 'geometry': record['semantics']['geometry'],
                'sourceEndpoints':record['semantics'].get('sourceEndpoints'),
                'participants': [{**{k:v for k,v in p.items() if k != 'site'}, 'site': {k:v for k,v in p['site'].items() if k in {'id','kind','residue','contactAtom'}}} for p in record['semantics']['participants']], 'direction': record['semantics']['direction'],
                'evidence': {k:v for k,v in record['semantics']['evidence'].items() if k != 'supportingRecords'}})
    observations = best['semantics'].get('sourceObservations') or [observation(best)]
    observations += other['semantics'].get('sourceObservations') or [observation(other)]
    seen = set()
    best['semantics']['sourceObservations'] = []
    for item in observations:
        token = json.dumps(item, sort_keys=True)
        if token not in seen:
            seen.add(token)
            best['semantics']['sourceObservations'].append(item)
    for participant in best['semantics']['participants']:
        counterpart = next((p for p in other['semantics']['participants'] if p['role'] == participant['role'] and p['site']['id'] == participant['site']['id']), None)
        if counterpart:
            for key in ('roleProvenance', 'atomTypes'):
                participant[key] = sorted(set(participant.get(key, [])) | set(counterpart.get(key, [])))
            left_charge, right_charge = participant.get('charge') or {}, counterpart.get('charge') or {}
            left_value, right_value = left_charge.get('value'), right_charge.get('value')
            incompatible_charge = left_value is not None and right_value is not None and left_value != right_value
            if incompatible_charge:
                evidence['ambiguityFlags'] = sorted(set(evidence['ambiguityFlags']) | {'charge_state_observations_differ'})
                evidence['level'] = 'ambiguous'
    evidence['supportingRecords'] = sum(r['semantics']['evidence'].get('supportingRecords',1) for r in (first,second))
    return best


def _assert_interaction(raw, residue_a, residue_b, aliases, residue_atoms_index=None,
                        base_pair_support_counts=None, base_pair_pair_stats=None):
    residue_a = _hydrate_contact_partner(residue_a, raw.get("bgn"), residue_atoms_index)
    residue_b = _hydrate_contact_partner(residue_b, raw.get("end"), residue_atoms_index)
    if hasattr(residue_atoms_index, 'chemical_nodes'):
        for residue, node in ((residue_a, raw.get('bgn')), (residue_b, raw.get('end'))):
            if isinstance(node, dict):
                residue_atoms_index.chemical_nodes[(residue.get('chain'),residue.get('seq'),_primary_contact_atom_name(residue.get('atom')))] = node
    result = _assert_interaction_core(raw, residue_a, residue_b, aliases, residue_atoms_index,
                                     base_pair_support_counts, base_pair_pair_stats)
    # Pair identity is residue-level context, independent of atom-level H-bond evidence.
    if result.get("family") in {"hbond", "polar_contact"} and not result.get("basePair"):
        key = _unordered_residue_pair_key(residue_a, residue_b, prefix="basepair_support:")
        stats = (base_pair_pair_stats or {}).get(key, {})
        fa, fb = _nucleic_base_family(residue_a.get("resName")), _nucleic_base_family(residue_b.get("resName"))
        canonical = _is_canonical_base_pair_hbond_pair(
            res_name_a=residue_a.get("resName", ""), atom_name_a=residue_a.get("atom", ""), element_a=residue_a.get("element", ""),
            res_name_b=residue_b.get("resName", ""), atom_name_b=residue_b.get("atom", ""), element_b=residue_b.get("element", ""))
        edge_pair = _is_nucleobase_pairing_edge_atom(residue_a.get("resName"), residue_a.get("atom")) and _is_nucleobase_pairing_edge_atom(residue_b.get("resName"), residue_b.get("atom"))
        if fa and fb and edge_pair and stats.get("mutualBestMatch") and not _is_sequence_adjacent_nucleotide_pair(
                residue_a=residue_a, residue_b=residue_b, base_family_a=fa, base_family_b=fb):
            result["basePair"] = {"isCanonicalAtomPattern": canonical,
                "family": _base_pair_family(fa, fb), "annotation": "watson_crick_candidate" if canonical else "paired_residue_context",
                "supportingPolarPairs": int(stats.get("supportCount") or 0), "score": _build_base_pair_score_payload(stats)}
    if result.get('family') in {'invalid_contact', 'covalent_bond'}:
        return result
    result['semantics'] = _build_interaction_semantics(raw, result, residue_a, residue_b, residue_atoms_index)
    semantics = result['semantics']
    if result.get('family') == 'pi_cation' and 'cation_site_unresolved' in semantics['evidence']['ambiguityFlags']:
        result.update(family='proximal', confidence='low', reason_dropped='cation_charge_state_unresolved', debugOnly=True)
        result['semantics'] = _build_interaction_semantics(raw, result, residue_a, residue_b, residue_atoms_index)
    if result.get('family') == 'hbond':
        result['subtype'] = 'hbond_geometrically_supported' if semantics['evidence']['level'] == 'geometrically_supported' else 'hbond_candidate'
        if semantics['evidence']['level'] != 'geometrically_supported' and result.get('confidence') == 'high':
            result['confidence'] = 'medium'
    return result


def _assert_interaction_core(
    raw: dict,
    residue_a: dict,
    residue_b: dict,
    aliases: ChainAliases,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]] = None,
    base_pair_support_counts: Optional[Dict[str, int]] = None,
    base_pair_pair_stats: Optional[Dict[str, dict]] = None,
) -> dict:
    node_a = raw.get("bgn") if isinstance(raw.get("bgn"), dict) else {}
    node_b = raw.get("end") if isinstance(raw.get("end"), dict) else {}

    res_name_a = str(
        node_a.get("label_comp_id")
        or node_a.get("auth_comp_id")
        or residue_a.get("resName")
        or ""
    ).strip().upper()
    res_name_b = str(
        node_b.get("label_comp_id")
        or node_b.get("auth_comp_id")
        or residue_b.get("resName")
        or ""
    ).strip().upper()
    atom_name_a = _primary_contact_atom_name(
        node_a.get("auth_atom_id") or node_a.get("label_atom_id") or residue_a.get("atom")
    )
    atom_name_b = _primary_contact_atom_name(
        node_b.get("auth_atom_id") or node_b.get("label_atom_id") or residue_b.get("atom")
    )
    element_a = str(
        residue_a.get("element")
        or node_a.get("type_symbol")
        or guess_element(atom_name_a)
        or ""
    ).strip().upper()
    element_b = str(
        residue_b.get("element")
        or node_b.get("type_symbol")
        or guess_element(atom_name_b)
        or ""
    ).strip().upper()

    terms = set(_normalize_arpeggio_contact_terms(raw.get("contact")))
    distance, distance_resolution = _resolve_contact_distance_value(
        raw,
        residue_a,
        residue_b,
        residue_atoms_index,
    )

    evidence: List[str] = []
    if distance_resolution == "distance_recomputed_from_coordinates":
        _append_evidence(evidence, "distance_recomputed_from_coordinates")
    elif distance_resolution == "distance_missing_or_invalid":
        _append_evidence(evidence, "distance_missing_or_invalid")

    identity_issue = _resolve_contact_identity_issue(
        raw_contact=raw,
        residue_a=residue_a,
        residue_b=residue_b,
        node_a=node_a,
        node_b=node_b,
        atom_name_a=atom_name_a,
        atom_name_b=atom_name_b,
        aliases=aliases,
    )
    if identity_issue:
        for token in identity_issue.get("evidence") or []:
            _append_evidence(evidence, token)
        _append_evidence(evidence, "invalid_contact_preclassified")
        return {
            "family": "invalid_contact",
            "confidence": "low",
            "evidence": evidence,
            "reason_dropped": str(identity_issue.get("reason") or "identity_or_mapping_gate"),
            "debugOnly": True,
        }

    covalent_false_positive = _is_likely_covalent_nonbonded_false_positive(
        distance=distance,
        element_a=element_a,
        atom_name_a=atom_name_a,
        res_name_a=res_name_a,
        element_b=element_b,
        atom_name_b=atom_name_b,
        res_name_b=res_name_b,
        terms=terms,
    )
    if covalent_false_positive:
        _append_evidence(evidence, "covalent_neighbor_filtered")
        if distance is not None:
            _append_evidence(evidence, "distance_covalent_like")
        return {
            "family": "covalent_bond",
            "confidence": "high",
            "evidence": evidence,
            "reason_dropped": "covalent_neighbor",
            "excludeFromNoncovalent": True,
        }

    preclassification_metal_context = _resolve_metal_contact_context(
        res_name_a=res_name_a,
        atom_name_a=atom_name_a,
        element_a=element_a,
        res_name_b=res_name_b,
        atom_name_b=atom_name_b,
        element_b=element_b,
        terms=terms,
        distance=distance,
    )
    preclassified_impossible = _resolve_impossible_contact_preclassification(
        distance=distance,
        element_a=element_a,
        element_b=element_b,
        suspect_invalid_mapping=False,
        metal_coordination_supported=bool(
            preclassification_metal_context
            and preclassification_metal_context.get("coordination_supported")
        ),
    )
    if isinstance(preclassified_impossible, dict):
        for token in preclassified_impossible.get("evidence") or []:
            _append_evidence(evidence, token)
        return {
            "family": str(preclassified_impossible.get("family") or "clash"),
            "confidence": str(preclassified_impossible.get("confidence") or "low"),
            "evidence": evidence,
            "reason_dropped": str(preclassified_impossible.get("reason_dropped") or "preclassification_gate"),
            "debugOnly": str(preclassified_impossible.get("family") or "").strip().lower() == "invalid_contact",
        }

    interaction_type = str(raw.get("type") or "").strip().lower()
    plausible_category = _classify_arpeggio_contact(raw)

    chain_token_a = str(residue_a.get("chain") or "").strip()
    chain_token_b = str(residue_b.get("chain") or "").strip()
    seq_token_a = str(residue_a.get("seq") or "").strip()
    seq_token_b = str(residue_b.get("seq") or "").strip()
    same_residue_pair = bool(
        chain_token_a
        and chain_token_b
        and seq_token_a
        and seq_token_b
        and chain_token_a == chain_token_b
        and seq_token_a == seq_token_b
    )
    if same_residue_pair:
        same_residue_metal_context = bool(
            element_a in METAL_ELEMENTS
            or element_b in METAL_ELEMENTS
            or terms.intersection(METAL_CONTACT_TERMS)
            or plausible_category == "metal_coordination"
        )
        same_residue_aromatic_context = bool(
            plausible_category in {"pi_pi", "pi_cation"}
            or interaction_type in {
                "plane-plane",
                "group-group",
                "atom-plane",
                "plane-atom",
                "group-plane",
                "plane-group",
            }
            or terms.intersection({"AROMATIC", "CATIONPI"})
        )
        if not (same_residue_metal_context or same_residue_aromatic_context):
            _append_evidence(evidence, "same_residue_nonbonded_filtered")
            _append_evidence(evidence, "invalid_contact_preclassified")
            return {
                "family": "invalid_contact",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "same_residue_nonbonded_contact",
                "debugOnly": True,
            }

    has_explicit_hbond_term = bool(terms.intersection(HBOND_EXPLICIT_TERMS))
    has_polar_fallback_term = bool(terms.intersection(HBOND_POLAR_FALLBACK_TERMS))
    has_hbond_like_term = has_explicit_hbond_term or has_polar_fallback_term
    donor_capable_a, acceptor_capable_a, roles_known_a = _partner_hbond_roles(node_a, res_name_a, atom_name_a, element_a)
    donor_capable_b, acceptor_capable_b, roles_known_b = _partner_hbond_roles(node_b, res_name_b, atom_name_b, element_b)
    donor_acceptor_consistent = (donor_capable_a and acceptor_capable_b) or (donor_capable_b and acceptor_capable_a)
    if not (roles_known_a and roles_known_b):
        _append_evidence(evidence, "chemical_roles_partially_unknown")
        # An explicit engine H-bond includes atom typing and hydrogen geometry.
        # Do not discard it solely because an older contact export omitted roles.
        if has_explicit_hbond_term and element_a in POLAR_CONTACT_ELEMENTS and element_b in POLAR_CONTACT_ELEMENTS:
            donor_acceptor_consistent = True
            _append_evidence(evidence, "upstream_hbond_typing_evidence")
    elif isinstance(node_a.get("atom_types"), list) or isinstance(node_b.get("atom_types"), list):
        _append_evidence(evidence, "arpeggio_chemical_atom_types")
    if raw.get("hydrogen_typing_source"):
        _append_evidence(evidence, str(raw["hydrogen_typing_source"]))
    if has_explicit_hbond_term:
        _append_evidence(evidence, "upstream_hbond_geometry_supported")
    strict_hbond_distance_ok = _hbond_distance_is_within_limits(distance, element_a, element_b)
    hbond_assert_distance_ok = bool(
        distance is None or distance <= HBOND_EXPLICIT_MAX_DISTANCE
    )
    polar_fallback_distance_ok = bool(
        distance is not None and distance <= HBOND_POLAR_FALLBACK_MAX_DISTANCE
    )
    angle_value = _extract_arpeggio_hbond_angle(raw)
    angle_available = angle_value is not None
    angle_passed = bool(angle_available and angle_value >= HBOND_STRONG_ANGLE_MIN)
    upstream_angle_supported = bool(has_explicit_hbond_term and raw.get("hbond_geometry_source") and angle_available and angle_value >= 90.0)
    if raw.get("hbond_geometry_source"):
        _append_evidence(evidence, str(raw["hbond_geometry_source"]))
    angle_proxy_value: Optional[float] = None
    angle_proxy_method: Optional[str] = None
    # A heavy-atom hydroxyl axis cannot determine the freely rotating O-H direction.
    # Only an exported donor-H-acceptor angle is directional validation.
    angle_proxy_available = angle_proxy_value is not None
    angle_proxy_passed = bool(
        angle_proxy_available and angle_proxy_value is not None and angle_proxy_value >= HBOND_PROXY_ANGLE_MIN
    )
    angle_proxy_strong = bool(
        angle_proxy_available
        and angle_proxy_value is not None
        and angle_proxy_value >= HBOND_PROXY_STRONG_ANGLE_MIN
    )
    angle_proxy_failed_hard = bool(
        angle_proxy_available
        and angle_proxy_value is not None
        and angle_proxy_value < HBOND_PROXY_FAIL_ANGLE_MAX
    )
    angle_checked_failed = bool(angle_available and not angle_passed)
    angle_proxy_failed = bool(angle_proxy_available and not angle_proxy_passed)

    base_family_a = _nucleic_base_family(res_name_a)
    base_family_b = _nucleic_base_family(res_name_b)
    nucleotide_residue_pair = bool(base_family_a and base_family_b)
    base_atom_a = _is_nucleobase_atom(res_name_a, atom_name_a)
    base_atom_b = _is_nucleobase_atom(res_name_b, atom_name_b)
    nucleobase_pair = bool(nucleotide_residue_pair and base_atom_a and base_atom_b)
    base_pair_edge_atom_a = _is_nucleobase_pairing_edge_atom(res_name_a, atom_name_a)
    base_pair_edge_atom_b = _is_nucleobase_pairing_edge_atom(res_name_b, atom_name_b)
    glycosidic_atom_a = _is_nucleobase_glycosidic_atom(res_name_a, atom_name_a)
    glycosidic_atom_b = _is_nucleobase_glycosidic_atom(res_name_b, atom_name_b)
    base_pair_edge_eligible = bool(
        nucleobase_pair
        and base_pair_edge_atom_a
        and base_pair_edge_atom_b
        and not glycosidic_atom_a
        and not glycosidic_atom_b
    )
    nucleotide_backbone_contact = bool(
        nucleotide_residue_pair
        and not nucleobase_pair
        and (
            _is_nucleotide_backbone_atom_name(atom_name_a)
            or _is_nucleotide_backbone_atom_name(atom_name_b)
        )
    )
    nucleotide_backbone_op_pair = bool(
        _is_nucleic_backbone_oxygen_or_phosphate_site(
            res_name=res_name_a,
            atom_name=atom_name_a,
            element=element_a,
        )
        and _is_nucleic_backbone_oxygen_or_phosphate_site(
            res_name=res_name_b,
            atom_name=atom_name_b,
            element=element_b,
        )
    )
    phosphate_sugar_oxygen_pair = _is_nucleic_phosphate_sugar_oxygen_pair(
        res_name_a=res_name_a,
        atom_name_a=atom_name_a,
        element_a=element_a,
        res_name_b=res_name_b,
        atom_name_b=atom_name_b,
        element_b=element_b,
    )
    nucleotide_backbone_oxygen_neighborhood_pair = _is_nucleic_backbone_oxygen_neighborhood_pair(
        res_name_a=res_name_a,
        atom_name_a=atom_name_a,
        element_a=element_a,
        res_name_b=res_name_b,
        atom_name_b=atom_name_b,
        element_b=element_b,
    )
    adjacent_nucleotide_linkage_contact = _is_likely_adjacent_nucleotide_linkage_contact(
        residue_a=residue_a,
        residue_b=residue_b,
        base_family_a=base_family_a,
        base_family_b=base_family_b,
        atom_name_a=atom_name_a,
        atom_name_b=atom_name_b,
    )
    sequence_adjacent_nucleotide_pair = _is_sequence_adjacent_nucleotide_pair(
        residue_a=residue_a,
        residue_b=residue_b,
        base_family_a=base_family_a,
        base_family_b=base_family_b,
    )
    canonical_base_pair = _is_canonical_base_pair_hbond_pair(
        res_name_a=res_name_a,
        atom_name_a=atom_name_a,
        element_a=element_a,
        res_name_b=res_name_b,
        atom_name_b=atom_name_b,
        element_b=element_b,
    )
    base_pair_distance_ok = bool(
        distance is not None and distance <= BASE_PAIR_CANDIDATE_MAX_DISTANCE
    )
    canonical_base_pair_template_match = bool(canonical_base_pair and base_pair_distance_ok)
    weak_acceptor_site_a = _is_weak_hbond_acceptor_site(
        res_name=res_name_a,
        atom_name=atom_name_a,
        element=element_a,
    )
    weak_acceptor_site_b = _is_weak_hbond_acceptor_site(
        res_name=res_name_b,
        atom_name=atom_name_b,
        element=element_b,
    )
    weak_acceptor_involved = bool(
        (donor_capable_a and weak_acceptor_site_b and acceptor_capable_b)
        or (donor_capable_b and weak_acceptor_site_a and acceptor_capable_a)
    )
    polar_role_known_a = bool(donor_capable_a or acceptor_capable_a)
    polar_role_known_b = bool(donor_capable_b or acceptor_capable_b)
    both_polar_elements = bool(
        element_a in POLAR_CONTACT_ELEMENTS
        and element_b in POLAR_CONTACT_ELEMENTS
    )
    polar_candidate_pair = bool(
        donor_acceptor_consistent
        and both_polar_elements
        and distance is not None
        and distance <= HBOND_POLAR_FALLBACK_MAX_DISTANCE
    )
    histidine_donor_involved = bool(
        (
            _is_histidine_protonation_dependent_donor_site(
                res_name=res_name_a,
                atom_name=atom_name_a,
                element=element_a,
            )
            and donor_capable_a
            and acceptor_capable_b
        )
        or (
            _is_histidine_protonation_dependent_donor_site(
                res_name=res_name_b,
                atom_name=atom_name_b,
                element=element_b,
            )
            and donor_capable_b
            and acceptor_capable_a
        )
    )
    extreme_short_no_hbond_distance = bool(
        distance is not None
        and distance < HBOND_EXTREME_SHORT_DISTANCE
        and {element_a, element_b} == {"N", "O"}
    )
    hbond_borderline_distance = bool(
        distance is not None and distance >= HBOND_BORDERLINE_DISTANCE
    )
    hbond_candidate_distance_ok = bool(
        distance is not None and distance <= HBOND_CANDIDATE_MAX_DISTANCE
    )
    polar_site_a = bool(element_a in POLAR_CONTACT_ELEMENTS and polar_role_known_a)
    polar_site_b = bool(element_b in POLAR_CONTACT_ELEMENTS and polar_role_known_b)
    packing_rejected_for_polar_site = bool(polar_site_a or polar_site_b)
    packing_atom_a = _is_packing_contact_atom_candidate(element_a, atom_name_a, res_name_a)
    packing_atom_b = _is_packing_contact_atom_candidate(element_b, atom_name_b, res_name_b)
    packing_eligible_pair = bool(packing_atom_a and packing_atom_b)
    packing_limit = min(
        PACKING_CONTACT_MAX_DISTANCE,
        _vdw_radius(element_a) + _vdw_radius(element_b) + PACKING_CONTACT_TOLERANCE,
    )
    base_pair_candidate = bool(
        base_pair_edge_eligible
        and donor_acceptor_consistent
        and base_pair_distance_ok
    )
    base_pair_residue_key = _unordered_residue_pair_key(
        residue_a,
        residue_b,
        prefix="basepair_support:",
    )
    base_pair_pair_stat = (base_pair_pair_stats or {}).get(base_pair_residue_key)
    base_pair_score_payload = _build_base_pair_score_payload(base_pair_pair_stat)
    base_pair_mutual_best_match = False
    base_pair_coplanarity_supported = False
    if isinstance(base_pair_pair_stat, dict):
        base_pair_mutual_best_match = bool(base_pair_pair_stat.get("mutualBestMatch"))
        base_pair_coplanarity_supported = bool(base_pair_pair_stat.get("coplanaritySupported"))
    base_pair_primary_partner = bool(base_pair_mutual_best_match)
    base_pair_support_pair_count = int((base_pair_support_counts or {}).get(base_pair_residue_key, 0))
    if isinstance(base_pair_pair_stat, dict):
        base_pair_support_pair_count = int(
            base_pair_pair_stat.get("supportCount") or base_pair_support_pair_count
        )
    base_pair_multi_polar_support = base_pair_support_pair_count >= BASE_PAIR_MIN_POLAR_PAIR_SUPPORT
    base_pair_ring_plane_candidate = False
    base_pair_ring_geometry: Dict[str, float] = {}
    ring_normal_angle: Optional[float] = None
    ring_interplanar_distance: Optional[float] = None
    ring_lateral_offset: Optional[float] = None
    if isinstance(base_pair_pair_stat, dict):
        ring_normal_angle = _coerce_float(base_pair_pair_stat.get("ringNormalAngle"))
        ring_interplanar_distance = _coerce_float(base_pair_pair_stat.get("ringInterplanarDistance"))
        ring_lateral_offset = _coerce_float(base_pair_pair_stat.get("ringLateralOffset"))
        if ring_normal_angle is not None:
            base_pair_ring_geometry["ring_normal_angle"] = ring_normal_angle
        if ring_interplanar_distance is not None:
            base_pair_ring_geometry["ring_interplanar_distance"] = ring_interplanar_distance
        if ring_lateral_offset is not None:
            base_pair_ring_geometry["ring_lateral_offset"] = ring_lateral_offset
        base_pair_ring_plane_candidate = bool(base_pair_pair_stat.get("coplanaritySupported"))
    if nucleobase_pair and base_pair_candidate and not canonical_base_pair_template_match:
        if not base_pair_ring_geometry:
            base_pair_ring_geometry = _compute_ring_geometry_metrics(
                residue_a,
                residue_b,
                residue_atoms_index,
            )
            ring_normal_angle = _coerce_float(base_pair_ring_geometry.get("ring_normal_angle"))
            ring_interplanar_distance = _coerce_float(base_pair_ring_geometry.get("ring_interplanar_distance"))
            ring_lateral_offset = _coerce_float(base_pair_ring_geometry.get("ring_lateral_offset"))
        if (
            ring_normal_angle is not None
            and ring_interplanar_distance is not None
            and ring_lateral_offset is not None
            and ring_normal_angle <= BASE_PAIR_RING_PLANE_MAX_NORMAL_ANGLE
            and ring_interplanar_distance <= BASE_PAIR_RING_PLANE_MAX_INTERPLANAR_DISTANCE
            and BASE_PAIR_RING_PLANE_MIN_LATERAL_OFFSET
            <= ring_lateral_offset
            <= BASE_PAIR_RING_PLANE_MAX_LATERAL_OFFSET
        ):
            base_pair_ring_plane_candidate = True
    base_pair_blocked_by_sequence_adjacency = bool(
        nucleobase_pair and sequence_adjacent_nucleotide_pair
    )
    base_ring_carbon_pair = bool(
        nucleobase_pair
        and element_a == "C"
        and element_b == "C"
    )

    def _classify_polar_candidate_with_hbond_promotion() -> dict:
        _append_evidence(evidence, "polar_candidate_pair")
        _append_evidence(evidence, "donor_acceptor_consistent")
        _append_evidence(evidence, "distance_ok_for_polar")

        if angle_available:
            _append_evidence(evidence, "angle_checked")
            if angle_checked_failed:
                _append_evidence(evidence, "angle_failed")
                return {
                    "family": "polar_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "angle_failed",
                }
            _append_evidence(evidence, "angle_passed")
        elif angle_proxy_available:
            _append_evidence(evidence, "angle_proxy_checked")
            if angle_proxy_method:
                _append_evidence(evidence, f"angle_proxy_method_{angle_proxy_method}")
            if angle_proxy_failed:
                _append_evidence(evidence, "angle_proxy_failed")
                if angle_proxy_failed_hard:
                    _append_evidence(evidence, "angle_proxy_failed_hard")
                return {
                    "family": "polar_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "angle_failed",
                }
            _append_evidence(evidence, "angle_proxy_passed")
        else:
            _append_evidence(evidence, "angle_missing")
            _append_evidence(evidence, "angle_not_evaluated")

        if distance is not None and distance < HBOND_EXTREME_SHORT_DISTANCE:
            _append_evidence(evidence, "unusually_short_distance")
            _append_evidence(evidence, "suspected_geometry_issue")
            _append_evidence(evidence, "suspected_clash")
            if not (angle_passed or angle_proxy_passed):
                _append_evidence(evidence, "requires_angle_to_confirm")
                return {
                    "family": "polar_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "extreme_short_hbond_distance",
                }

        promotion_eligible = bool(
            donor_acceptor_consistent
            and both_polar_elements
            and strict_hbond_distance_ok
            and hbond_candidate_distance_ok
        )
        if promotion_eligible:
            _append_evidence(evidence, "distance_ok_for_hbond")
            _append_evidence(evidence, "promoted_from_polar_by_distance_and_roles")
            confidence = (
                "medium"
                if distance is not None
                and distance <= HBOND_CANDIDATE_MEDIUM_CONFIDENCE_MAX_DISTANCE
                else "low"
            )
            return {
                "family": "hbond",
                "subtype": "hbond_candidate",
                "confidence": confidence,
                "evidence": evidence,
            }

        return {
            "family": "polar_contact",
            "confidence": "low",
            "evidence": evidence,
        }

    def _classify_base_pair_candidate() -> Optional[dict]:
        canonical_base_pair_assertable = bool(
            canonical_base_pair_template_match
            and base_pair_primary_partner
            and not base_pair_blocked_by_sequence_adjacency
        )
        base_pair_asserted = bool(
            canonical_base_pair_assertable
            or (
                base_pair_multi_polar_support
                and (has_explicit_hbond_term or angle_passed)
                and base_pair_primary_partner
                and not base_pair_blocked_by_sequence_adjacency
            )
        )
        if not (nucleobase_pair and (canonical_base_pair_template_match or base_pair_candidate)):
            return None
        _append_evidence(evidence, "nucleobase_pair")
        if base_pair_blocked_by_sequence_adjacency:
            _append_evidence(evidence, "adjacent_nucleotide_pair_blocked")
            _append_evidence(evidence, "sequence_adjacent_blocked")
        if base_pair_edge_eligible:
            _append_evidence(evidence, "base_pair_edge_atoms_eligible")
        if canonical_base_pair_template_match:
            _append_evidence(evidence, "base_pair_atom_pattern_matched")
        elif canonical_base_pair:
            _append_evidence(evidence, "base_pair_template_distance_out_of_range")
        if base_pair_support_pair_count >= BASE_PAIR_MIN_POLAR_PAIR_SUPPORT:
            _append_evidence(evidence, "base_pair_support_pairs_ge_2")
        elif base_pair_support_pair_count == 1:
            _append_evidence(evidence, "base_pair_support_single_pair")
        else:
            _append_evidence(evidence, "base_pair_support_pairs_lt_2")
        if glycosidic_atom_a or glycosidic_atom_b:
            _append_evidence(evidence, "base_pair_glycosidic_atom_excluded")
        if donor_acceptor_consistent:
            _append_evidence(evidence, "donor_acceptor_consistent")
        if base_pair_distance_ok:
            _append_evidence(evidence, "distance_ok_for_base_pair")
        else:
            _append_evidence(evidence, "distance_out_of_range_for_base_pair")
        if has_polar_fallback_term:
            _append_evidence(evidence, "polar_term_present")
        if has_explicit_hbond_term:
            _append_evidence(evidence, "hbond_term_present")
        if base_pair_candidate:
            _append_evidence(evidence, "polar_candidate_pair")
        if base_pair_score_payload:
            _append_evidence(evidence, "base_pair_score_available")
            if bool(base_pair_score_payload.get("mutualBestMatch")):
                _append_evidence(evidence, "base_pair_mutual_best_match")
            else:
                _append_evidence(evidence, "base_pair_not_mutual_best_match")
            if bool(base_pair_score_payload.get("coplanaritySupported")):
                _append_evidence(evidence, "base_pair_coplanarity_supported")
            else:
                _append_evidence(evidence, "base_pair_coplanarity_not_supported")
            if base_pair_primary_partner:
                _append_evidence(evidence, "base_pair_primary_partner")
            else:
                _append_evidence(evidence, "base_pair_non_primary_partner")
        if base_pair_ring_geometry:
            _append_evidence(evidence, "base_pair_ring_geometry_available")
            if ring_normal_angle is not None:
                if ring_normal_angle <= BASE_PAIR_RING_PLANE_MAX_NORMAL_ANGLE:
                    _append_evidence(evidence, "base_pair_ring_normal_angle_ok")
                else:
                    _append_evidence(evidence, "base_pair_ring_normal_angle_out_of_range")
            if ring_interplanar_distance is not None:
                if ring_interplanar_distance <= BASE_PAIR_RING_PLANE_MAX_INTERPLANAR_DISTANCE:
                    _append_evidence(evidence, "base_pair_ring_interplanar_distance_ok")
                else:
                    _append_evidence(evidence, "base_pair_ring_interplanar_distance_out_of_range")
            if ring_lateral_offset is not None:
                if (
                    BASE_PAIR_RING_PLANE_MIN_LATERAL_OFFSET
                    <= ring_lateral_offset
                    <= BASE_PAIR_RING_PLANE_MAX_LATERAL_OFFSET
                ):
                    _append_evidence(evidence, "base_pair_ring_lateral_offset_ok")
                else:
                    _append_evidence(evidence, "base_pair_ring_lateral_offset_out_of_range")
            if base_pair_ring_plane_candidate:
                _append_evidence(evidence, "base_pair_ring_plane_geometry_supported")
            else:
                _append_evidence(evidence, "base_pair_ring_plane_geometry_not_supported")
        elif base_pair_candidate and not canonical_base_pair_template_match:
            _append_evidence(evidence, "base_pair_ring_geometry_missing")
        if angle_available:
            _append_evidence(evidence, "angle_checked")
            if angle_passed:
                _append_evidence(evidence, "angle_passed")
            else:
                _append_evidence(evidence, "angle_failed")
        else:
            _append_evidence(evidence, "angle_missing")

        if base_pair_asserted:
            confidence = "low"
            if canonical_base_pair_assertable:
                confidence = "medium"
                if has_explicit_hbond_term and strict_hbond_distance_ok and donor_acceptor_consistent and angle_passed:
                    confidence = "high"
            elif base_pair_multi_polar_support:
                confidence = "medium"
            base_pair_payload = {
                "isCanonicalAtomPattern": canonical_base_pair,
                "family": _base_pair_family(base_family_a, base_family_b),
                "annotation": (
                    "watson_crick_candidate"
                    if canonical_base_pair
                    else "noncanonical_candidate"
                ),
                "supportingPolarPairs": base_pair_support_pair_count,
            }
            if base_pair_score_payload:
                base_pair_payload["score"] = base_pair_score_payload

            # Pairing is residue context, never a replacement for atom chemistry.
            # Use actual endpoint typing; a canonical atom-name pattern cannot
            # override a donor/donor or otherwise incompatible engine assignment.
            complementary_roles = bool((donor_capable_a and acceptor_capable_b) or
                                       (donor_capable_b and acceptor_capable_a))
            hbond_compatible = bool(complementary_roles and strict_hbond_distance_ok and
                                    hbond_candidate_distance_ok and not angle_checked_failed)
            _append_evidence(evidence, "base_pair_context_retained_on_atom_contact")
            if not hbond_compatible:
                _append_evidence(evidence, "base_pair_atom_hbond_chemistry_or_geometry_unresolved")
            return {
                "family": "hbond" if hbond_compatible else "polar_contact",
                "subtype": "hbond_candidate" if hbond_compatible else "base_pair_polar_candidate",
                "confidence": "medium" if hbond_compatible else "low",
                "evidence": evidence,
                "basePair": base_pair_payload,
            }

        _append_evidence(evidence, "base_pair_assertion_not_met")
        _append_evidence(
            evidence,
            "base_pair_requires_template_or_multi_pair_support",
        )
        if base_pair_score_payload and not base_pair_primary_partner and not base_pair_blocked_by_sequence_adjacency:
            _append_evidence(evidence, "base_pair_non_primary_partner_suppressed")
            result = {
                "family": "polar_proximal",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "base_pair_non_primary_partner",
                "debugOnly": True,
            }
            result["basePair"] = {
                "isCanonicalAtomPattern": canonical_base_pair,
                "family": _base_pair_family(base_family_a, base_family_b),
                "annotation": (
                    "watson_crick_candidate"
                    if canonical_base_pair
                    else "noncanonical_candidate"
                ),
                "supportingPolarPairs": base_pair_support_pair_count,
                "score": base_pair_score_payload,
            }
            return result
        fallback_confidence = "low"
        if (
            distance is not None
            and distance < BASE_PAIR_SINGLE_PAIR_LOW_CONFIDENCE_DISTANCE
            and donor_acceptor_consistent
            and has_explicit_hbond_term
        ):
            fallback_confidence = "medium"
        result = {
            "family": "polar_contact",
            "confidence": fallback_confidence,
            "evidence": evidence,
            "reason_dropped": (
                "sequence_adjacent"
                if base_pair_blocked_by_sequence_adjacency
                else "base_pair_assertion_not_met"
            ),
        }
        if base_pair_score_payload:
            result["basePair"] = {
                "isCanonicalAtomPattern": canonical_base_pair,
                "family": _base_pair_family(base_family_a, base_family_b),
                "annotation": (
                    "watson_crick_candidate"
                    if canonical_base_pair
                    else "noncanonical_candidate"
                ),
                "supportingPolarPairs": base_pair_support_pair_count,
                "score": base_pair_score_payload,
            }
        return result
    if nucleotide_backbone_contact:
        _append_evidence(evidence, "nucleotide_backbone_contact")
    if adjacent_nucleotide_linkage_contact:
        _append_evidence(evidence, "adjacent_nucleotide_linkage_contact")
        if distance is not None and distance <= HBOND_POLAR_FALLBACK_MAX_DISTANCE:
            _append_evidence(evidence, "distance_ok_for_polar")
        elif distance is not None:
            _append_evidence(evidence, "distance_out_of_range_for_polar")
        _append_evidence(evidence, "suppressed_adjacent_backbone_linkage")
        return {
            "family": "polar_contact",
            "confidence": "low",
            "evidence": evidence,
            "reason_dropped": "adjacent_nucleotide_linkage_contact",
            "debugOnly": True,
        }

    metal_context = _resolve_metal_contact_context(
        res_name_a=res_name_a,
        atom_name_a=atom_name_a,
        element_a=element_a,
        res_name_b=res_name_b,
        atom_name_b=atom_name_b,
        element_b=element_b,
        terms=terms,
        distance=distance,
    )
    metal_involved = bool(metal_context)
    has_metal_terms = bool(metal_context and metal_context.get("explicit_term"))
    if metal_context and (
        has_metal_terms
        or plausible_category == "metal_coordination"
        or metal_involved
    ):
        metal_side = str(metal_context.get("metal_side") or "").strip()
        metal_element = str(metal_context.get("metal_element") or "").strip().upper()
        donor_element = str(metal_context.get("donor_element") or "").strip().upper()
        donor_valid = bool(metal_context.get("donor_valid"))
        distance_ok = bool(metal_context.get("distance_ok"))

        if has_metal_terms:
            _append_evidence(evidence, "metal_term_present")
        if metal_side:
            _append_evidence(evidence, "metal_partner_detected")

        if donor_valid:
            _append_evidence(evidence, "donor_element_valid")
        elif donor_element:
            _append_evidence(evidence, "donor_element_invalid")
        else:
            _append_evidence(evidence, "donor_element_unknown")

        if distance_ok:
            _append_evidence(evidence, "distance_ok_for_metal")
        else:
            _append_evidence(evidence, "distance_out_of_range_for_metal")
        _append_evidence(evidence, "coordination_number_unknown")

        if donor_valid and distance_ok:
            confidence = "high" if has_metal_terms else "medium"
            return {
                "family": "metal_coordination",
                "confidence": confidence,
                "evidence": evidence,
                "metalElement": metal_element,
                "metalSide": metal_side,
                "donorElement": donor_element,
            }

        _append_evidence(evidence, "metal_coordination_not_asserted")
        if not donor_valid:
            _append_evidence(evidence, "metal_donor_not_eligible")
        if distance is None:
            _append_evidence(evidence, "metal_distance_missing")
        elif not distance_ok:
            _append_evidence(evidence, "metal_distance_outside_coordination_window")

    if nucleotide_backbone_op_pair and not adjacent_nucleotide_linkage_contact:
        _append_evidence(evidence, "nucleotide_backbone_op_pair")
        if phosphate_sugar_oxygen_pair and not donor_acceptor_consistent:
            _append_evidence(evidence, "phosphate_sugar_oxygen_pair")
            if distance is not None and distance <= POLAR_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_polar")
                return {
                    "family": "polar_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "phosphate_sugar_oxygen_hbond_blocked",
                }
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                _append_evidence(evidence, "suppressed_backbone_proximity")
                return {
                    "family": "polar_proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "nucleotide_backbone_proximity_suppressed",
                    "debugOnly": True,
                }
        elif nucleotide_backbone_oxygen_neighborhood_pair and not donor_acceptor_consistent:
            _append_evidence(evidence, "phosphate_backbone_oxygen_neighborhood_pair")
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                _append_evidence(evidence, "suppressed_backbone_proximity")
                return {
                    "family": "polar_proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "phosphate_backbone_oxygen_neighborhood_suppressed",
                    "debugOnly": True,
                }
            if distance is not None:
                _append_evidence(evidence, "distance_out_of_range_for_proximal")
        elif (
            not donor_acceptor_consistent
            and distance is not None
            and distance <= PROXIMAL_CONTACT_MAX_DISTANCE
        ):
            _append_evidence(evidence, "suppressed_backbone_proximity")
            _append_evidence(evidence, "distance_ok_for_proximal")
            return {
                "family": "polar_proximal",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "nucleotide_backbone_proximity_suppressed",
                "debugOnly": True,
            }

    has_ionic_term = bool("IONIC" in terms or plausible_category == "salt_bridge")
    if has_ionic_term and not metal_involved:
        _append_evidence(evidence, "ionic_term_present")
        salt_atom_name_a = atom_name_a
        salt_atom_name_b = atom_name_b
        salt_element_a = element_a
        salt_element_b = element_b
        salt_distance = distance
        salt_endpoint_override = _resolve_salt_bridge_endpoint_override(
            residue_a=residue_a,
            residue_b=residue_b,
            res_name_a=res_name_a,
            res_name_b=res_name_b,
            residue_atoms_index=residue_atoms_index,
        )
        if isinstance(salt_endpoint_override, dict):
            override_atom_a = _primary_contact_atom_name(salt_endpoint_override.get("atomNameA"))
            override_atom_b = _primary_contact_atom_name(salt_endpoint_override.get("atomNameB"))
            override_element_a = str(salt_endpoint_override.get("elementA") or "").strip().upper()
            override_element_b = str(salt_endpoint_override.get("elementB") or "").strip().upper()
            override_distance = _coerce_float(salt_endpoint_override.get("distance"))
            if override_atom_a and override_atom_b:
                if override_atom_a != atom_name_a or override_atom_b != atom_name_b:
                    _append_evidence(evidence, "salt_endpoints_reassigned_to_charged_sites")
                salt_atom_name_a = override_atom_a
                salt_atom_name_b = override_atom_b
            if override_element_a:
                salt_element_a = override_element_a
            if override_element_b:
                salt_element_b = override_element_b
            if override_distance is not None:
                if distance is None or abs(override_distance - distance) > 1e-3:
                    _append_evidence(evidence, "salt_distance_from_charged_sites")
                salt_distance = override_distance

        salt_payload_a = dict(residue_a, atom=salt_atom_name_a, element=salt_element_a)
        salt_payload_b = dict(residue_b, atom=salt_atom_name_b, element=salt_element_b)
        charge_a = _semantic_site_charge(_semantic_node(salt_payload_a,node_a,residue_atoms_index),salt_payload_a,residue_atoms_index)
        charge_b = _semantic_site_charge(_semantic_node(salt_payload_b,node_b,residue_atoms_index),salt_payload_b,residue_atoms_index)
        sign_a, sign_b = _semantic_charge_sign(charge_a), _semantic_charge_sign(charge_b)
        cation_site_a, cation_site_b = sign_a > 0, sign_b > 0
        anion_site_a, anion_site_b = sign_a < 0, sign_b < 0
        cation_site_unambiguous = bool(
            (cation_site_a and not cation_site_b) or (cation_site_b and not cation_site_a)
        )
        anion_site_unambiguous = bool(
            (anion_site_a and not anion_site_b) or (anion_site_b and not anion_site_a)
        )
        valid_pair = bool((cation_site_a and anion_site_b) or (cation_site_b and anion_site_a))
        if valid_pair:
            _append_evidence(evidence, "cation_anion_pattern_matched")
        else:
            _append_evidence(evidence, "cation_anion_pattern_not_matched")
        if cation_site_unambiguous:
            _append_evidence(evidence, "cation_site_unambiguous")
        else:
            _append_evidence(evidence, "cation_site_ambiguous_or_unknown")
        if anion_site_unambiguous:
            _append_evidence(evidence, "anion_site_unambiguous")
        else:
            _append_evidence(evidence, "anion_site_ambiguous_or_unknown")

        his_like_names = {"HIS", "HID", "HIE", "HSD", "HSE"}
        his_like_site_a = bool(
            res_name_a in his_like_names and salt_atom_name_a in {"ND1", "NE2"} and salt_element_a == "N"
        )
        his_like_site_b = bool(
            res_name_b in his_like_names and salt_atom_name_b in {"ND1", "NE2"} and salt_element_b == "N"
        )
        if (his_like_site_a or his_like_site_b) and not (
            res_name_a in PROTONATED_HISTIDINE_RESIDUE_NAMES
            or res_name_b in PROTONATED_HISTIDINE_RESIDUE_NAMES
        ):
            _append_evidence(evidence, "histidine_protonation_unknown_for_salt")

        distance_ok = bool(salt_distance is not None and salt_distance <= SALT_BRIDGE_MAX_DISTANCE)
        distance_confident = bool(
            salt_distance is not None and salt_distance <= SALT_BRIDGE_CONFIDENT_DISTANCE
        )
        if distance_ok:
            _append_evidence(evidence, "distance_ok_for_salt")
            if distance_confident:
                _append_evidence(evidence, "distance_ok_for_salt_confident")
            else:
                _append_evidence(evidence, "distance_near_salt_upper_bound")
        else:
            _append_evidence(evidence, "distance_out_of_range_for_salt")

        is_polymer_a = res_name_a in POLYMER_RESIDUES
        is_polymer_b = res_name_b in POLYMER_RESIDUES
        phosphate_a = _is_phosphate_oxygen_atom_name(salt_atom_name_a) and salt_element_a == "O"
        phosphate_b = _is_phosphate_oxygen_atom_name(salt_atom_name_b) and salt_element_b == "O"
        protein_cation_a = res_name_a in STANDARD_AMINO_RESIDUES and _is_salt_bridge_cation_site(
            res_name_a, salt_atom_name_a, salt_element_a
        )
        protein_cation_b = res_name_b in STANDARD_AMINO_RESIDUES and _is_salt_bridge_cation_site(
            res_name_b, salt_atom_name_b, salt_element_b
        )
        has_polymer_charge_context = (
            (is_polymer_a and is_polymer_b)
            or (protein_cation_a and phosphate_b)
            or (protein_cation_b and phosphate_a)
        )
        non_polymer_endpoint = not is_polymer_a or not is_polymer_b
        if non_polymer_endpoint and "+" not in res_name_a and "-" not in res_name_a and "+" not in res_name_b and "-" not in res_name_b:
            _append_evidence(evidence, "ligand_charge_unknown")

        salt_assertion_ok = bool(
            cation_site_unambiguous
            and anion_site_unambiguous
            and valid_pair
            and distance_ok
        )
        if not salt_assertion_ok:
            _append_evidence(evidence, "salt_bridge_assertion_not_met")
            if has_hbond_like_term or polar_candidate_pair:
                _append_evidence(evidence, "salt_bridge_failed_fallback_to_directional_polar")
            else:
                _append_evidence(evidence, "ionic_term_demoted_to_polar")
                if (
                    salt_distance is not None
                    and salt_distance <= POLAR_CONTACT_MAX_DISTANCE
                    and salt_element_a in {"N", "O", "S", "SE"}
                    and salt_element_b in {"N", "O", "S", "SE"}
                ):
                    _append_evidence(evidence, "distance_ok_for_polar")
                result = {
                    "family": "polar_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "salt_bridge_constraints_not_met",
                }
                if salt_atom_name_a and salt_atom_name_a != atom_name_a:
                    result["atomOverrideA"] = salt_atom_name_a
                if salt_atom_name_b and salt_atom_name_b != atom_name_b:
                    result["atomOverrideB"] = salt_atom_name_b
                if salt_element_a and salt_element_a != element_a:
                    result["elementOverrideA"] = salt_element_a
                if salt_element_b and salt_element_b != element_b:
                    result["elementOverrideB"] = salt_element_b
                if salt_distance is not None and (distance is None or abs(salt_distance - distance) > 1e-3):
                    result["distanceOverride"] = salt_distance
                return result

        if salt_assertion_ok:
            confidence = "low"
            if distance_confident and has_polymer_charge_context:
                confidence = "high"
            elif distance_confident:
                confidence = "medium"
            result = {
                "family": "salt_bridge",
                "confidence": confidence,
                "evidence": evidence,
            }
            if salt_atom_name_a and salt_atom_name_a != atom_name_a:
                result["atomOverrideA"] = salt_atom_name_a
            if salt_atom_name_b and salt_atom_name_b != atom_name_b:
                result["atomOverrideB"] = salt_atom_name_b
            if salt_element_a and salt_element_a != element_a:
                result["elementOverrideA"] = salt_element_a
            if salt_element_b and salt_element_b != element_b:
                result["elementOverrideB"] = salt_element_b
            if salt_distance is not None and (distance is None or abs(salt_distance - distance) > 1e-3):
                result["distanceOverride"] = salt_distance
            return result

    halogen_context = _resolve_halogen_bond_context(
        res_name_a=res_name_a,
        atom_name_a=atom_name_a,
        element_a=element_a,
        res_name_b=res_name_b,
        atom_name_b=atom_name_b,
        element_b=element_b,
        terms=terms,
        distance=distance,
    )
    halogen_like = bool(
        halogen_context
        or plausible_category == "halogen_bond"
        or terms.intersection(HALOGEN_BOND_EXPLICIT_TERMS)
    )
    if halogen_like:
        explicit_halogen_term = bool(terms.intersection(HALOGEN_BOND_EXPLICIT_TERMS))
        if explicit_halogen_term:
            _append_evidence(evidence, "halogen_term_present")
        donor_side = str((halogen_context or {}).get("donor_side") or "").strip().upper()
        if not donor_side:
            if element_a in HALOGEN_BOND_DONOR_ELEMENTS and element_b in HALOGEN_BOND_ACCEPTOR_ELEMENTS:
                donor_side = "A"
            elif element_b in HALOGEN_BOND_DONOR_ELEMENTS and element_a in HALOGEN_BOND_ACCEPTOR_ELEMENTS:
                donor_side = "B"
        if donor_side == "A":
            donor_residue = residue_a
            donor_res_name = res_name_a
            donor_atom_name = atom_name_a
            donor_element = element_a
            acceptor_residue = residue_b
            acceptor_res_name = res_name_b
            acceptor_atom_name = atom_name_b
            acceptor_element = element_b
            acceptor_side = "B"
        elif donor_side == "B":
            donor_residue = residue_b
            donor_res_name = res_name_b
            donor_atom_name = atom_name_b
            donor_element = element_b
            acceptor_residue = residue_a
            acceptor_res_name = res_name_a
            acceptor_atom_name = atom_name_a
            acceptor_element = element_a
            acceptor_side = "A"
        else:
            donor_residue = {}
            donor_res_name = ""
            donor_atom_name = ""
            donor_element = ""
            acceptor_residue = {}
            acceptor_res_name = ""
            acceptor_atom_name = ""
            acceptor_element = ""
            acceptor_side = ""

        donor_element_token = str(donor_element or "").strip().upper()
        acceptor_element_token = str(acceptor_element or "").strip().upper()
        donor_valid = bool(donor_element_token in HALOGEN_BOND_DONOR_ELEMENTS)
        if donor_valid:
            _append_evidence(evidence, "halogen_donor_element_valid")
        elif donor_side:
            _append_evidence(evidence, "halogen_donor_element_invalid")
        acceptor_role = acceptor_capable_a if acceptor_side == "A" else acceptor_capable_b
        acceptor_roles_known = roles_known_a if acceptor_side == "A" else roles_known_b
        acceptor_valid = bool(
            acceptor_element_token in HALOGEN_BOND_ACCEPTOR_ELEMENTS
            and (acceptor_role or (explicit_halogen_term and not acceptor_roles_known))
        )
        if acceptor_element_token == "C":
            _append_evidence(evidence, "halogen_acceptor_carbon_forbidden")
            acceptor_valid = False
        elif acceptor_valid:
            _append_evidence(evidence, "halogen_acceptor_valid")
        elif acceptor_side:
            _append_evidence(evidence, "halogen_acceptor_invalid")

        distance_limit = _halogen_bond_distance_limit(donor_element_token)
        distance_ok = bool(
            distance is not None
            and distance >= HBOND_HEAVY_MIN_DISTANCE
            and distance <= distance_limit
        )
        if distance_ok:
            _append_evidence(evidence, "distance_ok_for_halogen")
        elif distance is not None:
            _append_evidence(evidence, "distance_out_of_range_for_halogen")
        halogen_overlap = _estimate_vdw_overlap(distance, donor_element_token, acceptor_element_token)
        halogen_overlap_limit = min(HYDROPHOBIC_MAX_ALLOWED_OVERLAP, SOFT_CLASH_PRECLASSIFY)
        halogen_overlap_ok = bool(
            halogen_overlap is None
            or halogen_overlap <= halogen_overlap_limit
        )
        if halogen_overlap is not None:
            if halogen_overlap_ok:
                _append_evidence(evidence, "halogen_overlap_within_limit")
            else:
                _append_evidence(evidence, "halogen_overlap_exceeds_limit")

        donor_anchor = _resolve_halogen_donor_anchor_atom(
            donor_residue,
            donor_atom_name,
            residue_atoms_index,
        )
        donor_bound_to_carbon = donor_anchor is not None
        if donor_bound_to_carbon:
            _append_evidence(evidence, "halogen_donor_bound_to_carbon")
        elif donor_side:
            _append_evidence(evidence, "halogen_donor_not_bound_to_carbon")

        halogen_angle = _compute_halogen_bond_angle(
            donor_residue=donor_residue,
            donor_atom_name=donor_atom_name,
            acceptor_residue=acceptor_residue,
            acceptor_atom_name=acceptor_atom_name,
            residue_atoms_index=residue_atoms_index,
        )
        angle_strong = bool(halogen_angle is not None and halogen_angle >= HALOGEN_BOND_STRONG_ANGLE_MIN)
        angle_medium = bool(halogen_angle is not None and halogen_angle >= HALOGEN_BOND_MEDIUM_ANGLE_MIN)
        if angle_strong:
            _append_evidence(evidence, "halogen_angle_strong")
        elif angle_medium:
            _append_evidence(evidence, "halogen_angle_medium")
        elif halogen_angle is not None:
            _append_evidence(evidence, "halogen_angle_failed")
        else:
            _append_evidence(evidence, "halogen_angle_missing")

        halogen_assertion_ok = bool(
            donor_valid
            and acceptor_valid
            and donor_bound_to_carbon
            and distance_ok
            and halogen_overlap_ok
            and angle_medium
        )
        if halogen_assertion_ok:
            confidence = "high" if angle_strong and explicit_halogen_term else "medium"
            subtype = "halogen_bond_strong" if angle_strong else "halogen_bond_candidate"
            result = {
                "family": "halogen_bond",
                "subtype": subtype,
                "confidence": confidence,
                "evidence": evidence,
            }
            if halogen_angle is not None:
                result["halogenAngle"] = round(halogen_angle, 3)
            return result

        if explicit_halogen_term:
            _append_evidence(evidence, "halogen_assertion_not_met")
            halogen_overlap_hard_conflict = bool(
                halogen_overlap is not None and not halogen_overlap_ok
            )
            halogen_identity_invalid = bool(
                (donor_side and not donor_valid)
                or (acceptor_side and not acceptor_valid)
                or (donor_side and donor_valid and not donor_bound_to_carbon)
            )
            halogen_geometry_weak = bool(
                donor_valid and acceptor_valid and donor_bound_to_carbon and distance_ok and not angle_medium
            )
            halogen_distance_or_geometry_poor = bool(
                donor_valid and acceptor_valid and donor_bound_to_carbon and not distance_ok
            )
            if halogen_overlap_hard_conflict:
                _append_evidence(evidence, "halogen_demoted_to_clash_by_overlap")
                return {
                    "family": "clash",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "halogen_overlap_too_high",
                }
            if (
                halogen_geometry_weak
                and distance is not None
                and distance <= POLAR_CONTACT_MAX_DISTANCE
            ):
                _append_evidence(evidence, "distance_ok_for_polar")
                _append_evidence(evidence, "halogen_demoted_to_polar_by_angle")
                return {
                    "family": "polar_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "halogen_angle_not_supported",
                }
            if (
                halogen_distance_or_geometry_poor
                and distance is not None
                and distance <= PROXIMAL_CONTACT_MAX_DISTANCE
            ):
                _append_evidence(evidence, "distance_ok_for_proximal")
                _append_evidence(evidence, "halogen_demoted_to_proximal_by_geometry")
                return {
                    "family": "proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "halogen_geometry_not_supported",
                    "debugOnly": True,
                }
            if halogen_identity_invalid:
                _append_evidence(evidence, "halogen_endpoint_identity_invalid")
                if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                    _append_evidence(evidence, "distance_ok_for_proximal")
                    return {
                        "family": "proximal",
                        "confidence": "low",
                        "evidence": evidence,
                        "reason_dropped": "halogen_endpoint_identity_invalid",
                        "debugOnly": True,
                    }
                return {
                    "family": "other",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "halogen_endpoint_identity_invalid",
                }
            if (
                distance is not None
                and distance <= POLAR_CONTACT_MAX_DISTANCE
                and (
                    (element_a in POLAR_CONTACT_ELEMENTS and element_b in POLAR_CONTACT_ELEMENTS)
                    or acceptor_valid
                )
            ):
                _append_evidence(evidence, "distance_ok_for_polar")
                _append_evidence(evidence, "halogen_demoted_to_polar")
                return {
                    "family": "polar_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "halogen_constraints_not_met",
                }
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                _append_evidence(evidence, "halogen_demoted_to_proximal")
                return {
                    "family": "proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "halogen_constraints_not_met",
                    "debugOnly": True,
                }

    pi_like = bool(
        plausible_category in {"pi_pi", "pi_cation"}
        or "CATIONPI" in terms
        or "AROMATIC" in terms
        or interaction_type in {"plane-plane", "group-group", "atom-plane", "plane-atom", "group-plane", "plane-group"}
    )
    if pi_like and not (has_hbond_like_term or polar_candidate_pair):
        pi_evidence = list(evidence)
        family = "pi_cation" if plausible_category == "pi_cation" or "CATIONPI" in terms else "pi_pi"
        ring_angle = _extract_ring_geometry_angle(raw)
        ring_centroid_distance = _coerce_float(raw.get("ring_centroid_distance"))
        ring_interplanar_distance = _extract_ring_geometry_interplanar_distance(raw)
        ring_lateral_offset = _extract_ring_geometry_lateral_offset(raw)
        ring_min_atom_distance = _coerce_float(raw.get("ring_min_atom_distance"))
        computed_ring_geometry = {}
        if (
            ring_angle is None
            or ring_centroid_distance is None
            or ring_interplanar_distance is None
            or ring_lateral_offset is None
        ):
            computed_ring_geometry = _compute_ring_geometry_metrics_for_contact(
                residue_a,
                residue_b,
                residue_atoms_index,
                atom_name_a=atom_name_a,
                atom_name_b=atom_name_b,
            )
            if ring_angle is None:
                ring_angle = _coerce_float(computed_ring_geometry.get("ring_normal_angle"))
            if ring_centroid_distance is None:
                ring_centroid_distance = _coerce_float(computed_ring_geometry.get("ring_centroid_distance"))
            if ring_interplanar_distance is None:
                ring_interplanar_distance = _coerce_float(computed_ring_geometry.get("ring_interplanar_distance"))
            if ring_lateral_offset is None:
                ring_lateral_offset = _coerce_float(computed_ring_geometry.get("ring_lateral_offset"))
        if ring_min_atom_distance is None:
            ring_min_atom_distance = _coerce_float(computed_ring_geometry.get("ring_min_atom_distance"))
        if ring_min_atom_distance is None:
            ring_min_atom_distance = _compute_ring_min_atom_distance(
                residue_a,
                residue_b,
                residue_atoms_index,
            )
            if ring_min_atom_distance is not None:
                _append_evidence(pi_evidence, "ring_min_atom_distance_computed")
            else:
                _append_evidence(pi_evidence, "ring_min_atom_distance_missing")
        ring_payload = _build_ring_metrics_payload(
            centroid_distance=ring_centroid_distance,
            min_atom_distance=ring_min_atom_distance,
            interplanar_distance=ring_interplanar_distance,
            lateral_offset=ring_lateral_offset,
            normal_angle=ring_angle,
        )

        if (
            ring_angle is not None
            or ring_centroid_distance is not None
            or ring_interplanar_distance is not None
            or ring_lateral_offset is not None
        ):
            _append_evidence(pi_evidence, "ring_geometry_available")
        else:
            _append_evidence(pi_evidence, "ring_geometry_missing")
        if computed_ring_geometry:
            _append_evidence(pi_evidence, "ring_geometry_computed")
        elif _has_ring_geometry_metrics(raw):
            _append_evidence(pi_evidence, "ring_geometry_reported_by_arpeggio")

        centroid_distance_ok = bool(
            ring_centroid_distance is not None
            and PI_PI_MIN_CENTROID_DISTANCE <= ring_centroid_distance <= PI_PI_MAX_CENTROID_DISTANCE
        )
        interplanar_distance_ok = bool(
            ring_interplanar_distance is not None
            and PI_PI_MIN_INTERPLANAR_DISTANCE <= ring_interplanar_distance <= PI_PI_MAX_INTERPLANAR_DISTANCE
        )
        lateral_offset_ok = bool(
            ring_lateral_offset is not None and ring_lateral_offset <= PI_PI_MAX_LATERAL_OFFSET
        )
        if centroid_distance_ok:
            _append_evidence(pi_evidence, "distance_ok_for_pi")
        elif ring_centroid_distance is not None:
            _append_evidence(pi_evidence, "distance_out_of_range_for_pi")
        elif distance is not None:
            _append_evidence(pi_evidence, "distance_missing_for_pi_geometry")

        if ring_interplanar_distance is not None:
            if interplanar_distance_ok:
                _append_evidence(pi_evidence, "interplanar_distance_ok_for_pi")
            else:
                _append_evidence(pi_evidence, "interplanar_distance_out_of_range_for_pi")
        if ring_lateral_offset is not None:
            if lateral_offset_ok:
                _append_evidence(pi_evidence, "lateral_offset_ok_for_pi")
            else:
                _append_evidence(pi_evidence, "lateral_offset_out_of_range_for_pi")

        subtype = None
        full_pi_geometry_available = bool(
            ring_angle is not None
            and ring_centroid_distance is not None
            and ring_interplanar_distance is not None
            and ring_lateral_offset is not None
        )
        if family == "pi_pi" and ring_angle is not None:
            if (
                ring_angle <= PI_PI_STACKED_MAX_NORMAL_ANGLE
                and centroid_distance_ok
                and interplanar_distance_ok
                and lateral_offset_ok
            ):
                subtype = "pi_pi_stacked"
                return {
                    "family": family,
                    "subtype": subtype,
                    "confidence": "medium",
                    "evidence": pi_evidence,
                    "ring": ring_payload,
                }
            elif (
                ring_angle >= PI_PI_TSHAPED_MIN_NORMAL_ANGLE
                and centroid_distance_ok
                and (ring_lateral_offset is None or ring_lateral_offset <= PI_PI_TSHAPED_MAX_LATERAL_OFFSET)
            ):
                subtype = "pi_pi_tshaped"
                return {
                    "family": family,
                    "subtype": subtype,
                    "confidence": "medium",
                    "evidence": pi_evidence,
                    "ring": ring_payload,
                }
            if full_pi_geometry_available:
                _append_evidence(pi_evidence, "pi_geometry_computed_but_failed_thresholds")
                if ring_min_atom_distance is not None:
                    if (
                        ring_min_atom_distance >= HYDROPHOBIC_MIN_DISTANCE
                        and ring_min_atom_distance <= AROMATIC_PACKING_MAX_DISTANCE
                    ):
                        _append_evidence(pi_evidence, "distance_ok_for_aromatic_packing")
                        _append_evidence(pi_evidence, "aromatic_reclassified_from_pi_geometry")
                        result = {
                            "family": "aromatic_packing",
                            "confidence": "medium",
                            "evidence": pi_evidence,
                            "reason_dropped": "pi_geometry_not_supported",
                            "ring": ring_payload,
                        }
                        if distance is None or abs(ring_min_atom_distance - distance) > 1e-3:
                            result["distanceOverride"] = ring_min_atom_distance
                        return result
                    _append_evidence(pi_evidence, "distance_out_of_range_for_aromatic_packing")
                    if ring_min_atom_distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                        _append_evidence(pi_evidence, "distance_ok_for_aromatic_proximal")
                        _append_evidence(pi_evidence, "aromatic_reclassified_from_pi_geometry")
                        result = {
                            "family": "aromatic_proximal",
                            "confidence": "low",
                            "evidence": pi_evidence,
                            "reason_dropped": "pi_geometry_not_supported",
                            "debugOnly": True,
                            "ring": ring_payload,
                        }
                        if ring_centroid_distance is not None and (distance is None or abs(ring_centroid_distance - distance) > 1e-3):
                            result["distanceOverride"] = ring_centroid_distance
                        return result
                    _append_evidence(pi_evidence, "distance_out_of_range_for_aromatic_proximal")
                elif distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                    _append_evidence(pi_evidence, "distance_ok_for_aromatic_proximal")
                    _append_evidence(pi_evidence, "aromatic_proximal_distance_fallback_contact_pair")
                    _append_evidence(pi_evidence, "aromatic_reclassified_from_pi_geometry")
                    return {
                        "family": "aromatic_proximal",
                        "confidence": "low",
                        "evidence": pi_evidence,
                        "reason_dropped": "pi_geometry_not_supported",
                        "debugOnly": True,
                        "ring": ring_payload,
                    }
            else:
                _append_evidence(pi_evidence, "pi_assertion_not_met_geometry_missing")
                return {
                    "family": family,
                    "subtype": None,
                    "confidence": "low",
                    "evidence": pi_evidence,
                    "ring": ring_payload,
                }
        elif family == "pi_cation":
            cation_geometry_supported = bool(
                centroid_distance_ok
                and (ring_lateral_offset is None or ring_lateral_offset <= PI_PI_TSHAPED_MAX_LATERAL_OFFSET)
            )
            if cation_geometry_supported:
                return {
                    "family": family,
                    "subtype": subtype,
                    "confidence": "medium",
                    "evidence": pi_evidence,
                    "ring": ring_payload,
                }
            if ring_centroid_distance is None:
                _append_evidence(pi_evidence, "pi_assertion_not_met_geometry_missing")
                return {
                    "family": family,
                    "subtype": subtype,
                    "confidence": "low",
                    "evidence": pi_evidence,
                    "ring": ring_payload,
                }
            _append_evidence(pi_evidence, "pi_geometry_not_supported_for_asserted_pi_cation")

    if has_hbond_like_term or plausible_category == "hbond":
        if phosphate_sugar_oxygen_pair and not donor_acceptor_consistent:
            _append_evidence(evidence, "phosphate_sugar_oxygen_pair")
            _append_evidence(evidence, "hbond_blocked_for_backbone_oxygen_pair")
            if distance is not None and distance <= POLAR_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_polar")
                return {
                    "family": "polar_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "phosphate_sugar_oxygen_hbond_blocked",
                }
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                return {
                    "family": "polar_proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "phosphate_sugar_oxygen_hbond_blocked",
                    "debugOnly": True,
                }
        elif nucleotide_backbone_oxygen_neighborhood_pair and not donor_acceptor_consistent:
            _append_evidence(evidence, "phosphate_backbone_oxygen_neighborhood_pair")
            _append_evidence(evidence, "hbond_blocked_for_backbone_oxygen_pair")
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                return {
                    "family": "polar_proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "phosphate_backbone_oxygen_neighborhood_suppressed",
                    "debugOnly": True,
                }
            if distance is not None:
                _append_evidence(evidence, "distance_out_of_range_for_proximal")
        if donor_acceptor_consistent and strict_hbond_distance_ok and hbond_assert_distance_ok:
            _append_evidence(evidence, "distance_ok_for_hbond")
            _append_evidence(evidence, "donor_acceptor_consistent")
            if histidine_donor_involved:
                _append_evidence(evidence, "protonation_dependent")
            if weak_acceptor_involved:
                _append_evidence(evidence, "weak_acceptor_site_present")
                if not angle_available and not angle_proxy_passed:
                    _append_evidence(evidence, "weak_acceptor_requires_geometry")
                    if distance is not None and distance <= HBOND_POLAR_FALLBACK_MAX_DISTANCE:
                        _append_evidence(evidence, "distance_ok_for_polar")
                    return {
                        "family": "polar_contact",
                        "confidence": "low",
                        "evidence": evidence,
                        "reason_dropped": "weak_acceptor_geometry_missing",
                    }
            confidence = "medium"
            subtype = "hbond_candidate"
            if angle_available:
                _append_evidence(evidence, "angle_checked")
                if angle_passed:
                    confidence = "high"
                    subtype = "hbond_confirmed"
                    if raw.get("hbond_geometry_source") in {"arpeggio_generated_hydrogens", "arpeggio_completed_hydrogens"}:
                        confidence = "medium"
                        subtype = "hbond_candidate"
                        _append_evidence(evidence, "modeled_hydrogen_geometry")
                    _append_evidence(evidence, "angle_passed")
                elif upstream_angle_supported:
                    _append_evidence(evidence, "upstream_hbond_angle_supported")
                    confidence = "medium"
                else:
                    _append_evidence(evidence, "angle_failed")
                    if distance is not None and distance <= HBOND_POLAR_FALLBACK_MAX_DISTANCE:
                        _append_evidence(evidence, "distance_ok_for_polar")
                    _append_evidence(evidence, "hbond_demoted_by_angle")
                    return {
                        "family": "polar_contact",
                        "confidence": "low",
                        "evidence": evidence,
                        "reason_dropped": "angle_failed",
                    }
            else:
                _append_evidence(evidence, "angle_missing")
                _append_evidence(evidence, "angle_not_evaluated")
                if distance is not None and distance < HBOND_UNUSUALLY_SHORT_DISTANCE:
                    _append_evidence(evidence, "unusually_short_distance")
                if angle_proxy_available:
                    _append_evidence(evidence, "angle_proxy_checked")
                    if angle_proxy_method:
                        _append_evidence(evidence, f"angle_proxy_method_{angle_proxy_method}")
                    if angle_proxy_passed:
                        _append_evidence(evidence, "angle_proxy_passed")
                        if angle_proxy_strong:
                            _append_evidence(evidence, "angle_proxy_strong")
                            confidence = "high"
                            subtype = "hbond_confirmed"
                        else:
                            subtype = "hbond_confirmed"
                    else:
                        _append_evidence(evidence, "angle_proxy_failed")
                        if angle_proxy_failed_hard:
                            _append_evidence(evidence, "angle_proxy_failed_hard")
                        if distance is not None and distance <= HBOND_POLAR_FALLBACK_MAX_DISTANCE:
                            _append_evidence(evidence, "distance_ok_for_polar")
                            _append_evidence(evidence, "hbond_demoted_by_angle_proxy")
                            return {
                                "family": "polar_contact",
                                "confidence": "low",
                                "evidence": evidence,
                                "reason_dropped": "angle_proxy_failed",
                            }
            if extreme_short_no_hbond_distance:
                _append_evidence(evidence, "suspected_geometry_issue")
                _append_evidence(evidence, "requires_angle_to_confirm")
                vdw_overlap = _estimate_vdw_overlap(distance, element_a, element_b)
                if vdw_overlap is not None and vdw_overlap > CLASH_SOFT_OVERLAP:
                    _append_evidence(evidence, "clash_soft_overlap")
                if not angle_available or not angle_passed:
                    if not angle_available:
                        _append_evidence(evidence, "extreme_short_distance_angle_missing")
                    else:
                        _append_evidence(evidence, "extreme_short_distance_angle_failed")
                    if vdw_overlap is not None and vdw_overlap > CLASH_HARD_OVERLAP_POLAR:
                        _append_evidence(evidence, "overlap_hard_for_clash")
                        _append_evidence(evidence, "overlap_ok_for_clash")
                        _append_evidence(evidence, "hbond_rejected_extreme_short_distance")
                        return {
                            "family": "clash",
                            "confidence": "low",
                            "evidence": evidence,
                            "reason_dropped": "extreme_short_hbond_distance",
                        }
                    if distance is not None and distance <= HBOND_POLAR_FALLBACK_MAX_DISTANCE:
                        _append_evidence(evidence, "distance_ok_for_polar")
                    _append_evidence(evidence, "hbond_rejected_extreme_short_distance")
                    return {
                        "family": "polar_contact",
                        "confidence": "low",
                        "evidence": evidence,
                        "reason_dropped": "extreme_short_hbond_distance",
                    }
                _append_evidence(evidence, "extreme_short_distance_angle_passed")
                _append_evidence(evidence, "hydrogen_placement_unverified")
                confidence = "medium"
                subtype = "hbond_candidate"
            if hbond_borderline_distance:
                _append_evidence(evidence, "distance_near_hbond_upper_bound")
                if not angle_passed and not angle_proxy_passed:
                    confidence = "low"
                    _append_evidence(evidence, "confidence_downgraded_distance_tier")
            if histidine_donor_involved and not angle_available:
                _append_evidence(evidence, "protonation_uncertain_without_angle")
                confidence = "low"
                subtype = "hbond_candidate"
            elif histidine_donor_involved and not (angle_passed or angle_proxy_passed):
                _append_evidence(evidence, "protonation_uncertain_angle_failed")
                confidence = "low"
                subtype = "hbond_candidate"
            return {
                "family": "hbond",
                "subtype": subtype,
                "confidence": confidence,
                "evidence": evidence,
            }

        if has_polar_fallback_term and polar_fallback_distance_ok:
            _append_evidence(evidence, "polar_term_present")
            if polar_candidate_pair:
                return _classify_polar_candidate_with_hbond_promotion()
            if both_polar_elements:
                _append_evidence(evidence, "both_atoms_polar_eligible")
            else:
                _append_evidence(evidence, "one_sided_polar_site")
            if polar_role_known_a and polar_role_known_b:
                _append_evidence(evidence, "polar_roles_known")
            elif polar_role_known_a or polar_role_known_b:
                _append_evidence(evidence, "polar_roles_partially_known")
            else:
                _append_evidence(evidence, "polar_roles_unknown")
            if donor_acceptor_consistent:
                _append_evidence(evidence, "donor_acceptor_consistent")
            elif polar_role_known_a or polar_role_known_b:
                _append_evidence(evidence, "donor_acceptor_noncomplementary")
            else:
                _append_evidence(evidence, "donor_acceptor_unknown")
            _append_evidence(evidence, "distance_ok_for_polar")
            if packing_eligible_pair and distance is not None and distance <= packing_limit:
                _append_evidence(evidence, "distance_ok_for_packing")
                _append_evidence(evidence, "polar_site_neighbor_packing")
                return {
                    "family": "packing_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "polar_term_noncomplementary",
                }
            _append_evidence(evidence, "distance_out_of_range_for_packing")
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                _append_evidence(evidence, "polar_site_neighbor_proximal")
                return {
                    "family": "polar_proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "polar_term_noncomplementary",
                    "debugOnly": True,
                }
            _append_evidence(evidence, "distance_out_of_range_for_proximal")
            return {
                "family": "other",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "polar_term_noncomplementary",
            }

        if has_explicit_hbond_term and (
            distance is None or distance <= HBOND_POLAR_FALLBACK_MAX_DISTANCE
        ):
            _append_evidence(evidence, "hbond_term_present")
            if polar_candidate_pair:
                return _classify_polar_candidate_with_hbond_promotion()
            if both_polar_elements:
                _append_evidence(evidence, "both_atoms_polar_eligible")
            else:
                _append_evidence(evidence, "one_sided_polar_site")
            if polar_role_known_a and polar_role_known_b:
                _append_evidence(evidence, "polar_roles_known")
            elif polar_role_known_a or polar_role_known_b:
                _append_evidence(evidence, "polar_roles_partially_known")
            else:
                _append_evidence(evidence, "polar_roles_unknown")
            if donor_acceptor_consistent:
                _append_evidence(evidence, "donor_acceptor_consistent")
            elif polar_role_known_a or polar_role_known_b:
                _append_evidence(evidence, "donor_acceptor_noncomplementary")
            else:
                _append_evidence(evidence, "donor_acceptor_unknown")
            _append_evidence(evidence, "angle_missing")
            if distance is not None:
                if distance <= HBOND_POLAR_FALLBACK_MAX_DISTANCE:
                    _append_evidence(evidence, "distance_ok_for_polar")
                else:
                    _append_evidence(evidence, "distance_out_of_range_for_polar")
            if distance is not None and packing_eligible_pair and distance <= packing_limit:
                _append_evidence(evidence, "distance_ok_for_packing")
                _append_evidence(evidence, "polar_site_neighbor_packing")
                return {
                    "family": "packing_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "hbond_term_noncomplementary",
                }
            if distance is not None and (not packing_eligible_pair or distance > packing_limit):
                _append_evidence(evidence, "distance_out_of_range_for_packing")
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                _append_evidence(evidence, "polar_site_neighbor_proximal")
                return {
                    "family": "polar_proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "hbond_term_noncomplementary",
                    "debugOnly": True,
                }
            if distance is not None:
                _append_evidence(evidence, "distance_out_of_range_for_proximal")
            return {
                "family": "other",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "hbond_term_noncomplementary",
            }

    base_pair_result = _classify_base_pair_candidate()
    if isinstance(base_pair_result, dict):
        return base_pair_result

    if polar_candidate_pair:
        return _classify_polar_candidate_with_hbond_promotion()

    if base_ring_carbon_pair:
        _append_evidence(evidence, "nucleobase_pair")
        _append_evidence(evidence, "base_ring_carbon_pair")
        ring_geometry = _compute_ring_geometry_metrics(
            residue_a,
            residue_b,
            residue_atoms_index,
        )
        ring_centroid_distance = _coerce_float(ring_geometry.get("ring_centroid_distance"))
        ring_interplanar_distance = _coerce_float(ring_geometry.get("ring_interplanar_distance"))
        ring_lateral_offset = _coerce_float(ring_geometry.get("ring_lateral_offset"))
        ring_normal_angle = _coerce_float(ring_geometry.get("ring_normal_angle"))
        ring_min_atom_distance = _compute_nucleobase_ring_min_atom_distance(
            residue_a,
            residue_b,
            residue_atoms_index,
        )
        if ring_min_atom_distance is not None:
            _append_evidence(evidence, "ring_min_atom_distance_computed")
        else:
            _append_evidence(evidence, "ring_min_atom_distance_missing")
        ring_payload = _build_ring_metrics_payload(
            centroid_distance=ring_centroid_distance,
            min_atom_distance=ring_min_atom_distance,
            interplanar_distance=ring_interplanar_distance,
            lateral_offset=ring_lateral_offset,
            normal_angle=ring_normal_angle,
        )

        aromatic_distance = ring_min_atom_distance
        if (
            aromatic_distance is not None
            and aromatic_distance >= AROMATIC_PACKING_MIN_DISTANCE
            and aromatic_distance <= AROMATIC_PACKING_MAX_DISTANCE
        ):
            _append_evidence(evidence, "distance_ok_for_aromatic_packing")
            if terms.intersection({"AROMATIC", "VDW", "PROXIMAL", "HYDROPHOBIC"}):
                _append_evidence(evidence, "aromatic_packing_term_context")
            result = {
                "family": "aromatic_packing",
                "confidence": "medium",
                "evidence": evidence,
                "ring": ring_payload,
            }
            if (
                distance is None
                or abs(aromatic_distance - distance) > 1e-3
            ):
                result["distanceOverride"] = aromatic_distance
            return result

        if aromatic_distance is not None:
            _append_evidence(evidence, "distance_out_of_range_for_aromatic_packing")
            if aromatic_distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_aromatic_proximal")
                if terms.intersection({"AROMATIC", "VDW", "PROXIMAL", "HYDROPHOBIC"}):
                    _append_evidence(evidence, "aromatic_proximal_term_context")
                result = {
                    "family": "aromatic_proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "aromatic_packing_distance_out_of_range",
                    "debugOnly": True,
                    "ring": ring_payload,
                }
                proximal_distance = ring_centroid_distance if ring_centroid_distance is not None else aromatic_distance
                if proximal_distance is not None and (distance is None or abs(proximal_distance - distance) > 1e-3):
                    result["distanceOverride"] = proximal_distance
                return result
            _append_evidence(evidence, "distance_out_of_range_for_aromatic_proximal")
        elif distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
            _append_evidence(evidence, "distance_ok_for_aromatic_proximal")
            _append_evidence(evidence, "aromatic_proximal_distance_fallback_contact_pair")
            if terms.intersection({"AROMATIC", "VDW", "PROXIMAL", "HYDROPHOBIC"}):
                _append_evidence(evidence, "aromatic_proximal_term_context")
            return {
                "family": "aromatic_proximal",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "aromatic_packing_distance_out_of_range",
                "debugOnly": True,
                "ring": ring_payload,
            }

    hydrophobic_atom_a = _is_hydrophobic_contact_atom_candidate(element_a, atom_name_a, res_name_a)
    hydrophobic_atom_b = _is_hydrophobic_contact_atom_candidate(element_b, atom_name_b, res_name_b)
    nonpolar_pair = hydrophobic_atom_a and hydrophobic_atom_b
    hydrophobic_pair_min_distance = _hydrophobic_min_distance_for_pair(element_a, element_b)
    hydrophobic_overlap = _estimate_vdw_overlap(distance, element_a, element_b)
    hydrophobic_overlap_ok = bool(
        hydrophobic_overlap is None or hydrophobic_overlap <= HYDROPHOBIC_MAX_ALLOWED_OVERLAP
    )
    hydrophobic_directional_block = bool(
        has_hbond_like_term
        or polar_candidate_pair
        or plausible_category in {"hbond", "halogen_bond"}
        or terms.intersection(HALOGEN_BOND_EXPLICIT_TERMS)
    )
    hydrophobic_distance_ok = bool(
        distance is not None
        and distance >= hydrophobic_pair_min_distance
        and distance <= ARPEGGIO_HYDROPHOBIC_MAX_DISTANCE
    )
    hydrophobic_like = bool(
        not base_ring_carbon_pair
        and not hydrophobic_directional_block
        and (
            plausible_category == "hydrophobic"
            or ("HYDROPHOBIC" in terms)
            or (terms.intersection({"PROXIMAL", "VDW", "VDW_CLASH"}) and nonpolar_pair)
        )
    )
    if hydrophobic_like:
        aromatic_context_result = _classify_aromatic_context_contact(
            residue_a=residue_a,
            residue_b=residue_b,
            res_name_a=res_name_a,
            res_name_b=res_name_b,
            atom_name_a=atom_name_a,
            atom_name_b=atom_name_b,
            distance=distance,
            terms=terms,
            evidence=evidence,
            residue_atoms_index=residue_atoms_index,
        )
        if isinstance(aromatic_context_result, dict):
            return aromatic_context_result

    if hydrophobic_like and nonpolar_pair:
        _append_evidence(evidence, "nonpolar_pair")
        if hydrophobic_overlap is not None:
            if hydrophobic_overlap_ok:
                _append_evidence(evidence, "hydrophobic_overlap_within_limit")
            else:
                _append_evidence(evidence, "hydrophobic_overlap_exceeds_limit")
        if not hydrophobic_overlap_ok:
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                _append_evidence(evidence, "hydrophobic_rejected_by_overlap")
                return {
                    "family": "proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "hydrophobic_overlap_too_high",
                    "debugOnly": True,
                }
            _append_evidence(evidence, "distance_out_of_range_for_proximal")
            return {
                "family": "other",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "hydrophobic_overlap_too_high",
            }
        non_polymer_a = bool(
            res_name_a
            and res_name_a not in POLYMER_RESIDUES
            and res_name_a not in WATER_RESIDUES
            and res_name_a not in METAL_ELEMENTS
        )
        non_polymer_b = bool(
            res_name_b
            and res_name_b not in POLYMER_RESIDUES
            and res_name_b not in WATER_RESIDUES
            and res_name_b not in METAL_ELEMENTS
        )
        aromatic_site_a = bool(
            res_name_a in AROMATIC_RESIDUES
            and element_a == "C"
            and not _is_protein_backbone_atom_name(atom_name_a)
        )
        aromatic_site_b = bool(
            res_name_b in AROMATIC_RESIDUES
            and element_b == "C"
            and not _is_protein_backbone_atom_name(atom_name_b)
        )
        if (aromatic_site_a and non_polymer_b) or (aromatic_site_b and non_polymer_a):
            _append_evidence(evidence, "aromatic_context_possible")
            aromatic_context_ring_geometry = _compute_ring_geometry_metrics_for_contact(
                residue_a,
                residue_b,
                residue_atoms_index,
                atom_name_a=atom_name_a,
                atom_name_b=atom_name_b,
            )
            ring_angle = _coerce_float(aromatic_context_ring_geometry.get("ring_normal_angle"))
            ring_centroid_distance = _coerce_float(aromatic_context_ring_geometry.get("ring_centroid_distance"))
            ring_interplanar_distance = _coerce_float(aromatic_context_ring_geometry.get("ring_interplanar_distance"))
            ring_lateral_offset = _coerce_float(aromatic_context_ring_geometry.get("ring_lateral_offset"))
            ring_min_atom_distance = _coerce_float(aromatic_context_ring_geometry.get("ring_min_atom_distance"))
            if ring_min_atom_distance is None:
                ring_min_atom_distance = _compute_ring_min_atom_distance(
                    residue_a,
                    residue_b,
                    residue_atoms_index,
                )
            aromatic_context_ring_payload = _build_ring_metrics_payload(
                centroid_distance=ring_centroid_distance,
                min_atom_distance=ring_min_atom_distance,
                interplanar_distance=ring_interplanar_distance,
                lateral_offset=ring_lateral_offset,
                normal_angle=ring_angle,
            )
            if (
                ring_angle is not None
                and ring_centroid_distance is not None
                and ring_interplanar_distance is not None
                and ring_lateral_offset is not None
            ):
                _append_evidence(evidence, "ring_geometry_available")
                _append_evidence(evidence, "ring_geometry_computed")
                centroid_distance_ok = bool(
                    PI_PI_MIN_CENTROID_DISTANCE <= ring_centroid_distance <= PI_PI_MAX_CENTROID_DISTANCE
                )
                interplanar_distance_ok = bool(
                    PI_PI_MIN_INTERPLANAR_DISTANCE <= ring_interplanar_distance <= PI_PI_MAX_INTERPLANAR_DISTANCE
                )
                lateral_offset_ok = bool(ring_lateral_offset <= PI_PI_MAX_LATERAL_OFFSET)
                if centroid_distance_ok:
                    _append_evidence(evidence, "distance_ok_for_pi")
                else:
                    _append_evidence(evidence, "distance_out_of_range_for_pi")
                if interplanar_distance_ok:
                    _append_evidence(evidence, "interplanar_distance_ok_for_pi")
                else:
                    _append_evidence(evidence, "interplanar_distance_out_of_range_for_pi")
                if lateral_offset_ok:
                    _append_evidence(evidence, "lateral_offset_ok_for_pi")
                else:
                    _append_evidence(evidence, "lateral_offset_out_of_range_for_pi")

                if centroid_distance_ok and interplanar_distance_ok and lateral_offset_ok:
                    subtype = None
                    if ring_angle <= PI_PI_STACKED_MAX_NORMAL_ANGLE:
                        subtype = "pi_pi_stacked"
                    elif (
                        ring_angle >= PI_PI_TSHAPED_MIN_NORMAL_ANGLE
                        and ring_lateral_offset <= PI_PI_TSHAPED_MAX_LATERAL_OFFSET
                    ):
                        subtype = "pi_pi_tshaped"
                    _append_evidence(evidence, "aromatic_context_promoted_to_pi_pi")
                    return {
                        "family": "pi_pi",
                        "subtype": subtype,
                        "confidence": "medium",
                        "evidence": evidence,
                        "ring": aromatic_context_ring_payload,
                    }
                _append_evidence(evidence, "aromatic_context_pi_geometry_not_supported")
                if ring_min_atom_distance is not None:
                    _append_evidence(evidence, "ring_min_atom_distance_computed")
                    if (
                        ring_min_atom_distance >= HYDROPHOBIC_MIN_DISTANCE
                        and ring_min_atom_distance <= AROMATIC_PACKING_MAX_DISTANCE
                    ):
                        _append_evidence(evidence, "distance_ok_for_aromatic_packing")
                        _append_evidence(evidence, "aromatic_reclassified_from_pi_context")
                        result = {
                            "family": "aromatic_packing",
                            "confidence": "medium",
                            "evidence": evidence,
                            "reason_dropped": "aromatic_context_pi_geometry_not_supported",
                            "ring": aromatic_context_ring_payload,
                        }
                        if distance is None or abs(ring_min_atom_distance - distance) > 1e-3:
                            result["distanceOverride"] = ring_min_atom_distance
                        return result
                    _append_evidence(evidence, "distance_out_of_range_for_aromatic_packing")
                    if ring_min_atom_distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                        _append_evidence(evidence, "distance_ok_for_aromatic_proximal")
                        _append_evidence(evidence, "aromatic_reclassified_from_pi_context")
                        result = {
                            "family": "aromatic_proximal",
                            "confidence": "low",
                            "evidence": evidence,
                            "reason_dropped": "aromatic_context_pi_geometry_not_supported",
                            "debugOnly": True,
                            "ring": aromatic_context_ring_payload,
                        }
                        if ring_centroid_distance is not None and (distance is None or abs(ring_centroid_distance - distance) > 1e-3):
                            result["distanceOverride"] = ring_centroid_distance
                        return result
                    _append_evidence(evidence, "distance_out_of_range_for_aromatic_proximal")
            else:
                _append_evidence(evidence, "ring_geometry_missing")
        if (
            _is_polar_or_charged_sidechain_residue(res_name_a)
            or _is_polar_or_charged_sidechain_residue(res_name_b)
        ):
            _append_evidence(evidence, "polar_residue_backbone_or_sidechain_involved")
        nonpolar_packing_semantic = bool(
            "aromatic_context_possible" in evidence
            or "polar_residue_backbone_or_sidechain_involved" in evidence
        )
        if hydrophobic_distance_ok:
            confidence = "low"
            if distance is not None and distance <= HYDROPHOBIC_MEDIUM_MAX_DISTANCE:
                confidence = (
                    "medium"
                    if "HYDROPHOBIC" in terms or plausible_category == "hydrophobic"
                    else "low"
                )
            else:
                _append_evidence(evidence, "distance_near_hydrophobic_upper_bound")
            if nonpolar_packing_semantic:
                _append_evidence(evidence, "distance_ok_for_packing")
                _append_evidence(evidence, "reclassified_nonpolar_packing")
                return {
                    "family": "packing_contact",
                    "confidence": confidence,
                    "evidence": evidence,
                    "reason_dropped": "nonpolar_packing_reclassified_from_hydrophobic",
                }
            _append_evidence(evidence, "distance_ok_for_hydrophobic")
            return {
                "family": "hydrophobic",
                "confidence": confidence,
                "evidence": evidence,
            }
        _append_evidence(evidence, "distance_out_of_range_for_hydrophobic")
        if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
            _append_evidence(evidence, "distance_ok_for_proximal")
            _append_evidence(evidence, "downgraded_from_hydrophobic")
            return {
                "family": "proximal",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "hydrophobic_distance_out_of_range",
                "debugOnly": True,
            }
        _append_evidence(evidence, "distance_out_of_range_for_proximal")
        return {
            "family": "other",
            "confidence": "low",
            "evidence": evidence,
            "reason_dropped": "hydrophobic_distance_out_of_range",
        }
    if hydrophobic_like and not nonpolar_pair:
        _append_evidence(evidence, "nonpolar_pair_rejected")
        if _is_likely_carbonyl_carbon_for_hydrophobic(res_name_a, atom_name_a) or _is_likely_carbonyl_carbon_for_hydrophobic(
            res_name_b, atom_name_b
        ):
            _append_evidence(evidence, "carbonyl_carbon_present")
        if _is_charged_group_associated_carbon_for_hydrophobic(
            res_name_a,
            atom_name_a,
        ) or _is_charged_group_associated_carbon_for_hydrophobic(
            res_name_b,
            atom_name_b,
        ):
            _append_evidence(evidence, "charged_group_core_carbon_present")
        if packing_rejected_for_polar_site:
            _append_evidence(evidence, "packing_rejected_polar_site")
            if both_polar_elements:
                _append_evidence(evidence, "both_atoms_polar_eligible")
            else:
                _append_evidence(evidence, "one_sided_polar_site")
            if polar_role_known_a and polar_role_known_b:
                _append_evidence(evidence, "polar_roles_known")
            elif polar_role_known_a or polar_role_known_b:
                _append_evidence(evidence, "polar_roles_partially_known")
            else:
                _append_evidence(evidence, "polar_roles_unknown")
            if donor_acceptor_consistent:
                _append_evidence(evidence, "donor_acceptor_consistent")
            elif polar_role_known_a or polar_role_known_b:
                _append_evidence(evidence, "donor_acceptor_noncomplementary")
            else:
                _append_evidence(evidence, "donor_acceptor_unknown")
            if polar_candidate_pair:
                _append_evidence(evidence, "distance_ok_for_polar")
                return {
                    "family": "polar_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "hydrophobic_nonpolar_pair_failed",
                }
            if distance is not None:
                _append_evidence(evidence, "distance_out_of_range_for_polar")
            if distance is not None and packing_eligible_pair and distance <= packing_limit:
                _append_evidence(evidence, "distance_ok_for_packing")
                _append_evidence(evidence, "polar_site_neighbor_packing")
                return {
                    "family": "packing_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "hydrophobic_nonpolar_pair_failed",
                }
            if distance is not None and (not packing_eligible_pair or distance > packing_limit):
                _append_evidence(evidence, "distance_out_of_range_for_packing")
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                _append_evidence(evidence, "polar_site_neighbor_proximal")
                return {
                    "family": "polar_proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "hydrophobic_nonpolar_pair_failed",
                    "debugOnly": True,
                }
            _append_evidence(evidence, "distance_out_of_range_for_proximal")
            return {
                "family": "other",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "hydrophobic_nonpolar_pair_failed",
            }
        if not packing_eligible_pair:
            _append_evidence(evidence, "packing_nonpolar_eligibility_failed")
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                _append_evidence(evidence, "hydrophobic_downgraded_to_proximal")
                return {
                    "family": "proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "hydrophobic_nonpolar_pair_failed",
                    "debugOnly": True,
                }
            _append_evidence(evidence, "distance_out_of_range_for_proximal")
            return {
                "family": "other",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "hydrophobic_nonpolar_pair_failed",
            }
        if distance is not None and distance <= packing_limit:
            _append_evidence(evidence, "distance_ok_for_packing")
            _append_evidence(evidence, "hydrophobic_downgraded_to_packing")
            return {
                "family": "packing_contact",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "hydrophobic_nonpolar_pair_failed",
            }
        if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
            _append_evidence(evidence, "distance_ok_for_proximal")
            _append_evidence(evidence, "hydrophobic_downgraded_to_proximal")
            return {
                "family": "proximal",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "hydrophobic_nonpolar_pair_failed",
                "debugOnly": True,
            }
        _append_evidence(evidence, "distance_out_of_range_for_packing")
        _append_evidence(evidence, "distance_out_of_range_for_proximal")
        return {
            "family": "other",
            "confidence": "low",
            "evidence": evidence,
            "reason_dropped": "hydrophobic_nonpolar_pair_failed",
        }

    if "VDW_CLASH" in terms:
        vdw_overlap = _estimate_vdw_overlap(distance, element_a, element_b)
        polar_precedence = bool(
            polar_candidate_pair
            or has_hbond_like_term
            or (
                nucleobase_pair
                and donor_acceptor_consistent
                and distance is not None
                and distance <= POLAR_PRECEDENCE_MAX_DISTANCE
            )
        )
        hard_overlap_cutoff = (
            CLASH_HARD_OVERLAP_POLAR if polar_precedence else CLASH_HARD_OVERLAP
        )
        if vdw_overlap is not None:
            if vdw_overlap > CLASH_SOFT_OVERLAP:
                _append_evidence(evidence, "clash_soft_overlap")
            if vdw_overlap > hard_overlap_cutoff:
                _append_evidence(evidence, "overlap_hard_for_clash")
                _append_evidence(evidence, "clash_term_present")
                if vdw_overlap is not None:
                    _append_evidence(evidence, "overlap_ok_for_clash")
                return {
                    "family": "clash",
                    "confidence": "low",
                    "evidence": evidence,
                }
        if polar_precedence:
            _append_evidence(evidence, "polar_precedence_applied")
        elif vdw_overlap is None:
            _append_evidence(evidence, "clash_overlap_unknown")
        else:
            _append_evidence(evidence, "tight_contact")
        if nucleobase_pair:
            _append_evidence(evidence, "nucleobase_pair")
            if (
                canonical_base_pair
                and distance is not None
                and distance <= BASE_PAIR_CANDIDATE_MAX_DISTANCE
            ):
                _append_evidence(evidence, "base_pair_atom_pattern_matched")
        if packing_rejected_for_polar_site:
            _append_evidence(evidence, "packing_rejected_polar_site")
            if both_polar_elements:
                _append_evidence(evidence, "both_atoms_polar_eligible")
            else:
                _append_evidence(evidence, "one_sided_polar_site")
            if polar_role_known_a and polar_role_known_b:
                _append_evidence(evidence, "polar_roles_known")
            elif polar_role_known_a or polar_role_known_b:
                _append_evidence(evidence, "polar_roles_partially_known")
            else:
                _append_evidence(evidence, "polar_roles_unknown")
            if donor_acceptor_consistent:
                _append_evidence(evidence, "donor_acceptor_consistent")
            elif polar_role_known_a or polar_role_known_b:
                _append_evidence(evidence, "donor_acceptor_noncomplementary")
            else:
                _append_evidence(evidence, "donor_acceptor_unknown")
            if polar_candidate_pair:
                _append_evidence(evidence, "distance_ok_for_polar")
                return {
                    "family": "polar_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "packing_rejected_polar_site",
                }
            if distance is not None:
                _append_evidence(evidence, "distance_out_of_range_for_polar")
            if distance is not None and packing_eligible_pair and distance <= packing_limit:
                _append_evidence(evidence, "distance_ok_for_packing")
                _append_evidence(evidence, "polar_site_neighbor_packing")
                return {
                    "family": "packing_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "packing_rejected_polar_site",
                }
            if distance is not None and (not packing_eligible_pair or distance > packing_limit):
                _append_evidence(evidence, "distance_out_of_range_for_packing")
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                _append_evidence(evidence, "polar_site_neighbor_proximal")
                return {
                    "family": "polar_proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "packing_rejected_polar_site",
                    "debugOnly": True,
                }
            _append_evidence(evidence, "distance_out_of_range_for_proximal")
            return {
                "family": "other",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "packing_rejected_polar_site",
            }
        if not packing_eligible_pair:
            _append_evidence(evidence, "packing_nonpolar_eligibility_failed")
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                _append_evidence(evidence, "proximal_term_present")
                return {
                    "family": "proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "packing_nonpolar_eligibility_failed",
                    "debugOnly": True,
                }
            _append_evidence(evidence, "distance_out_of_range_for_packing")
            return {
                "family": "other",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "packing_nonpolar_eligibility_failed",
            }
        if distance is not None and distance <= packing_limit:
            _append_evidence(evidence, "distance_ok_for_packing")
            _append_evidence(evidence, "vdw_term_present")
            return {
                "family": "packing_contact",
                "confidence": "low",
                "evidence": evidence,
            }
        if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
            _append_evidence(evidence, "distance_ok_for_proximal")
            _append_evidence(evidence, "proximal_term_present")
            return {
                "family": "proximal",
                "confidence": "low",
                "evidence": evidence,
                "debugOnly": True,
            }
        _append_evidence(evidence, "distance_out_of_range_for_packing")
        return {
            "family": "other",
            "confidence": "low",
            "evidence": evidence,
            "reason_dropped": "packing_distance_out_of_range",
        }

    if terms.intersection({"VDW", "PROXIMAL"}):
        if nucleobase_pair:
            _append_evidence(evidence, "nucleobase_pair")
            if (
                canonical_base_pair
                and distance is not None
                and distance <= BASE_PAIR_CANDIDATE_MAX_DISTANCE
            ):
                _append_evidence(evidence, "base_pair_atom_pattern_matched")
            if polar_candidate_pair:
                _append_evidence(evidence, "polar_candidate_pair")
        if packing_rejected_for_polar_site:
            _append_evidence(evidence, "packing_rejected_polar_site")
            if both_polar_elements:
                _append_evidence(evidence, "both_atoms_polar_eligible")
            else:
                _append_evidence(evidence, "one_sided_polar_site")
            if polar_role_known_a and polar_role_known_b:
                _append_evidence(evidence, "polar_roles_known")
            elif polar_role_known_a or polar_role_known_b:
                _append_evidence(evidence, "polar_roles_partially_known")
            else:
                _append_evidence(evidence, "polar_roles_unknown")
            if donor_acceptor_consistent:
                _append_evidence(evidence, "donor_acceptor_consistent")
            elif polar_role_known_a or polar_role_known_b:
                _append_evidence(evidence, "donor_acceptor_noncomplementary")
            else:
                _append_evidence(evidence, "donor_acceptor_unknown")
            if polar_candidate_pair:
                _append_evidence(evidence, "distance_ok_for_polar")
                return {
                    "family": "polar_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "packing_rejected_polar_site",
                }
            if distance is not None:
                _append_evidence(evidence, "distance_out_of_range_for_polar")
            if distance is not None and packing_eligible_pair and distance <= packing_limit:
                _append_evidence(evidence, "distance_ok_for_packing")
                _append_evidence(evidence, "polar_site_neighbor_packing")
                return {
                    "family": "packing_contact",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "packing_rejected_polar_site",
                }
            if distance is not None and (not packing_eligible_pair or distance > packing_limit):
                _append_evidence(evidence, "distance_out_of_range_for_packing")
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                _append_evidence(evidence, "polar_site_neighbor_proximal")
                return {
                    "family": "polar_proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "packing_rejected_polar_site",
                    "debugOnly": True,
                }
            _append_evidence(evidence, "distance_out_of_range_for_proximal")
            return {
                "family": "other",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "packing_rejected_polar_site",
            }
        if not packing_eligible_pair:
            _append_evidence(evidence, "packing_nonpolar_eligibility_failed")
            if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
                _append_evidence(evidence, "distance_ok_for_proximal")
                _append_evidence(evidence, "proximal_term_present")
                return {
                    "family": "proximal",
                    "confidence": "low",
                    "evidence": evidence,
                    "reason_dropped": "packing_nonpolar_eligibility_failed",
                    "debugOnly": True,
                }
            _append_evidence(evidence, "distance_out_of_range_for_packing")
            return {
                "family": "other",
                "confidence": "low",
                "evidence": evidence,
                "reason_dropped": "packing_nonpolar_eligibility_failed",
            }
        if distance is not None and distance <= packing_limit:
            _append_evidence(evidence, "distance_ok_for_packing")
            _append_evidence(evidence, "vdw_term_present")
            return {
                "family": "packing_contact",
                "confidence": "low",
                "evidence": evidence,
            }
        if distance is not None and distance <= PROXIMAL_CONTACT_MAX_DISTANCE:
            _append_evidence(evidence, "distance_ok_for_proximal")
            _append_evidence(evidence, "proximal_term_present")
            return {
                "family": "proximal",
                "confidence": "low",
                "evidence": evidence,
                "debugOnly": True,
            }
        _append_evidence(evidence, "distance_out_of_range_for_packing")
        return {
            "family": "other",
            "confidence": "low",
            "evidence": evidence,
            "reason_dropped": "packing_distance_out_of_range",
        }

    _append_evidence(evidence, "no_assertion_rule_matched")
    return {
        "family": "other",
        "confidence": "low",
        "evidence": evidence,
        "reason_dropped": "no_assertion_rule_matched",
    }


def _classify_arpeggio_contact(raw_contact: dict) -> str:
    interaction_type = str(raw_contact.get("type") or "").strip().lower()
    terms = set(_normalize_arpeggio_contact_terms(raw_contact.get("contact")))

    node_a = raw_contact.get("bgn") if isinstance(raw_contact.get("bgn"), dict) else {}
    node_b = raw_contact.get("end") if isinstance(raw_contact.get("end"), dict) else {}
    atom_name_a = str(node_a.get("auth_atom_id") or node_a.get("label_atom_id") or "").strip()
    atom_name_b = str(node_b.get("auth_atom_id") or node_b.get("label_atom_id") or "").strip()
    element_a = str(node_a.get("type_symbol") or "").strip().upper()
    element_b = str(node_b.get("type_symbol") or "").strip().upper()
    if not element_a:
        element_a = guess_element(atom_name_a).upper()
    if not element_b:
        element_b = guess_element(atom_name_b).upper()
    res_name_a = str(node_a.get("label_comp_id") or node_a.get("auth_comp_id") or "").strip().upper()
    res_name_b = str(node_b.get("label_comp_id") or node_b.get("auth_comp_id") or "").strip().upper()
    distance = _coerce_float(raw_contact.get("distance"))

    polar_elements = {"N", "O", "S", "SE"}
    is_hydrophobic_atom_a = _is_hydrophobic_contact_atom_candidate(
        element_a,
        atom_name_a,
        res_name_a,
    )
    is_hydrophobic_atom_b = _is_hydrophobic_contact_atom_candidate(
        element_b,
        atom_name_b,
        res_name_b,
    )
    hydrophobic_min_distance = _hydrophobic_min_distance_for_pair(element_a, element_b)
    hydrophobic_overlap = _estimate_vdw_overlap(distance, element_a, element_b)
    hydrophobic_overlap_ok = bool(
        hydrophobic_overlap is None or hydrophobic_overlap <= HYDROPHOBIC_MAX_ALLOWED_OVERLAP
    )
    hydrophobic_distance_ok = bool(
        distance is not None
        and distance >= hydrophobic_min_distance
        and distance <= ARPEGGIO_HYDROPHOBIC_MAX_DISTANCE
    )
    hydrophobic_pair_ok = bool(
        is_hydrophobic_atom_a
        and is_hydrophobic_atom_b
        and hydrophobic_distance_ok
        and hydrophobic_overlap_ok
    )

    metal_context = _resolve_metal_contact_context(
        res_name_a=res_name_a,
        atom_name_a=atom_name_a,
        element_a=element_a,
        res_name_b=res_name_b,
        atom_name_b=atom_name_b,
        element_b=element_b,
        terms=terms,
        distance=distance,
    )
    if metal_context:
        if bool(metal_context.get("coordination_supported")):
            return "metal_coordination"
    if "IONIC" in terms:
        if _is_valid_salt_bridge_pair(
            res_name_a=res_name_a,
            atom_name_a=atom_name_a,
            element_a=element_a,
            res_name_b=res_name_b,
            atom_name_b=atom_name_b,
            element_b=element_b,
        ):
            return "salt_bridge"
    has_valid_hbond_roles = _is_hbond_donor_acceptor_pair(
        res_name_a=res_name_a,
        atom_name_a=atom_name_a,
        element_a=element_a,
        res_name_b=res_name_b,
        atom_name_b=atom_name_b,
        element_b=element_b,
    )
    halogen_context = _resolve_halogen_bond_context(
        res_name_a=res_name_a,
        atom_name_a=atom_name_a,
        element_a=element_a,
        res_name_b=res_name_b,
        atom_name_b=atom_name_b,
        element_b=element_b,
        terms=terms,
        distance=distance,
    )
    if halogen_context:
        if bool(halogen_context.get("donor_valid")) and bool(halogen_context.get("acceptor_valid")):
            return "halogen_bond"
        if bool(halogen_context.get("explicit_term")):
            return "halogen_bond"
    if terms.intersection(HBOND_EXPLICIT_TERMS):
        if distance is None or distance <= HBOND_POLAR_FALLBACK_MAX_DISTANCE:
            if has_valid_hbond_roles:
                return "hbond"
    if "CATIONPI" in terms:
        return "pi_cation"

    if interaction_type in {"plane-plane", "group-group"}:
        return "pi_pi"
    if interaction_type in {"atom-plane", "plane-atom", "group-plane", "plane-group"}:
        if "CATIONPI" in terms:
            return "pi_cation"
        if terms.intersection(PI_NON_RING_SURFACE_TERMS) and hydrophobic_pair_ok:
            return "hydrophobic"
        if terms.intersection(PI_NON_RING_SURFACE_TERMS):
            return "other"

    if "AROMATIC" in terms and not terms.intersection(PI_NON_RING_SURFACE_TERMS):
        return "pi_pi"
    if terms.intersection(PI_NON_RING_SURFACE_TERMS) and hydrophobic_pair_ok:
        return "hydrophobic"
    if "HYDROPHOBIC" in terms and hydrophobic_pair_ok:
        return "hydrophobic"

    # Arpeggio frequently reports non-bonded close contacts as proximal/vdw without
    # a specific "hydrophobic" flag; map obvious nonpolar contacts into the existing UI bucket.
    if terms.intersection({"PROXIMAL", "VDW", "VDW_CLASH"}):
        if hydrophobic_pair_ok:
            return "hydrophobic"

    # Conservative fallback: only promote POLAR-like terms to hydrogen-bond when
    # donor/acceptor chemistry is consistent.
    if terms.intersection(HBOND_POLAR_FALLBACK_TERMS):
        if (
            distance is not None
            and distance <= HBOND_POLAR_FALLBACK_MAX_DISTANCE
            and element_a in polar_elements
            and element_b in polar_elements
            and has_valid_hbond_roles
        ):
            return "hbond"

    return "other"


def _format_arpeggio_res_seq(node: dict) -> str:
    auth_seq = str(node.get("auth_seq_id") or node.get("label_seq_id") or "").strip()
    if not auth_seq or auth_seq in {"?", "."}:
        auth_seq = "?"
    ins_code = str(node.get("pdbx_PDB_ins_code") or "").strip()
    if ins_code and ins_code not in {"?", "."}:
        auth_seq = f"{auth_seq}{ins_code}"
    return auth_seq


def _build_residue_payload_from_arpeggio_partner(
    node: dict,
    aliases: ChainAliases,
) -> Optional[dict]:
    if not isinstance(node, dict):
        return None
    res_name = str(node.get("label_comp_id") or node.get("auth_comp_id") or "").strip().upper()
    seq = _format_arpeggio_res_seq(node)
    atom_name_raw = node.get("auth_atom_id") or node.get("label_atom_id") or ""
    if isinstance(atom_name_raw, list):
        atom_name = ",".join(str(token or "").strip() for token in atom_name_raw if str(token or "").strip())
    else:
        atom_name = str(atom_name_raw or "").strip()
    chain = aliases.resolve_partner_chain(
        auth_chain=str(node.get("auth_asym_id") or "").strip(),
        label_chain=str(node.get("label_asym_id") or "").strip(),
        res_seq=seq,
        res_name=res_name,
        atom_name=atom_name,
    )
    if not chain or chain in {"?", "."}:
        return None
    payload = {
        "chain": chain,
        "resName": res_name,
        "seq": seq,
        "atom": atom_name,
    }
    element = str(node.get("type_symbol") or "").strip().upper()
    if element:
        payload["element"] = element
    return payload


def _residue_pair_identity_token(residue: dict) -> str:
    if not isinstance(residue, dict):
        return "?"
    chain = str(residue.get("chain") or "").strip()
    seq = str(residue.get("seq") or "").strip()
    if chain and seq:
        return f"{chain}:{seq}"
    # Fallback for malformed payloads where chain/seq are missing.
    res_name = str(residue.get("resName") or "").strip().upper()
    atom_name = str(residue.get("atom") or "").strip().upper()
    return f"{chain}:{res_name}:{seq}:{atom_name}"


def _unordered_residue_pair_key(residue_a: dict, residue_b: dict, prefix: str = "") -> str:
    token_a = _residue_pair_identity_token(residue_a)
    token_b = _residue_pair_identity_token(residue_b)
    if token_a <= token_b:
        return f"{prefix}{token_a}|{token_b}"
    return f"{prefix}{token_b}|{token_a}"


def _carboxylate_group_token(res_name: str, atom_name: str) -> str:
    residue = str(res_name or "").strip().upper()
    atom = _normalize_atom_name(atom_name)
    if not residue or not atom:
        return ""
    allowed = CARBOXYLATE_ACCEPTOR_ATOMS_BY_RESIDUE.get(residue)
    if not allowed or atom not in allowed:
        return ""
    if residue in {"ASP", "ASX"}:
        return "asp_carboxylate"
    if residue in {"GLU", "GLX"}:
        return "glu_carboxylate"
    return ""


def _resolve_hbond_carboxylate_dedupe_key(
    residue_a: dict,
    residue_b: dict,
) -> Optional[Tuple[str, str]]:
    if not isinstance(residue_a, dict) or not isinstance(residue_b, dict):
        return None
    res_name_a = str(residue_a.get("resName") or "").strip().upper()
    res_name_b = str(residue_b.get("resName") or "").strip().upper()
    atom_name_a = _primary_contact_atom_name(residue_a.get("atom"))
    atom_name_b = _primary_contact_atom_name(residue_b.get("atom"))
    if not res_name_a or not res_name_b or not atom_name_a or not atom_name_b:
        return None
    element_a = str(residue_a.get("element") or guess_element(atom_name_a) or "").strip().upper()
    element_b = str(residue_b.get("element") or guess_element(atom_name_b) or "").strip().upper()
    if not element_a or not element_b:
        return None

    a_donor = _is_hbond_donor_capable(
        res_name=res_name_a,
        atom_name=atom_name_a,
        element=element_a,
    )
    b_donor = _is_hbond_donor_capable(
        res_name=res_name_b,
        atom_name=atom_name_b,
        element=element_b,
    )
    a_acceptor = _is_hbond_acceptor_capable(
        res_name=res_name_a,
        atom_name=atom_name_a,
        element=element_a,
    )
    b_acceptor = _is_hbond_acceptor_capable(
        res_name=res_name_b,
        atom_name=atom_name_b,
        element=element_b,
    )

    donor_residue = None
    acceptor_residue = None
    acceptor_res_name = ""
    acceptor_atom_name = ""
    if a_donor and b_acceptor and not (b_donor and a_acceptor):
        donor_residue = residue_a
        acceptor_residue = residue_b
        acceptor_res_name = res_name_b
        acceptor_atom_name = atom_name_b
    elif b_donor and a_acceptor and not (a_donor and b_acceptor):
        donor_residue = residue_b
        acceptor_residue = residue_a
        acceptor_res_name = res_name_a
        acceptor_atom_name = atom_name_a
    else:
        return None

    group_token = _carboxylate_group_token(acceptor_res_name, acceptor_atom_name)
    if not group_token:
        return None
    donor_atom_key = _build_atom_key_from_payload(donor_residue)
    acceptor_residue_key = _residue_pair_identity_token(acceptor_residue)
    if not donor_atom_key or not acceptor_residue_key:
        return None
    dedupe_key = f"hbond_carboxylate:{donor_atom_key}|{acceptor_residue_key}|{group_token}"
    return dedupe_key, _normalize_atom_name(acceptor_atom_name)


def _unordered_residue_atom_pair_key(
    residue_a: dict,
    atom_name_a: str,
    residue_b: dict,
    atom_name_b: str,
    prefix: str = "",
) -> str:
    token_a = _residue_pair_identity_token(residue_a)
    token_b = _residue_pair_identity_token(residue_b)
    atom_a = _primary_contact_atom_name(atom_name_a) or "?"
    atom_b = _primary_contact_atom_name(atom_name_b) or "?"
    left = f"{token_a}:{atom_a}"
    right = f"{token_b}:{atom_b}"
    if left <= right:
        return f"{prefix}{left}|{right}"
    return f"{prefix}{right}|{left}"


def _contact_atom_feature_signature(
    *,
    res_name: str,
    atom_name: str,
    element: str,
) -> str:
    atom = _primary_contact_atom_name(atom_name)
    atom_element = str(element or "").strip().upper() or guess_element(atom).upper()
    donor = _is_hbond_donor_capable(
        res_name=res_name,
        atom_name=atom,
        element=atom_element,
    )
    acceptor = _is_hbond_acceptor_capable(
        res_name=res_name,
        atom_name=atom,
        element=atom_element,
    )
    hydrophobic = _is_hydrophobic_contact_atom_candidate(atom_element, atom, res_name)
    backbone = bool(_is_protein_backbone_atom_name(atom) or _is_nucleotide_backbone_atom_name(atom))
    return (
        f"{atom_element}:d{int(donor)}:a{int(acceptor)}:h{int(hydrophobic)}:b{int(backbone)}"
    )


def _build_preclassification_duplicate_key(raw: dict, residue_a: dict, residue_b: dict) -> str:
    if not isinstance(raw, dict) or not isinstance(residue_a, dict) or not isinstance(residue_b, dict):
        return ""
    node_a = raw.get("bgn") if isinstance(raw.get("bgn"), dict) else {}
    node_b = raw.get("end") if isinstance(raw.get("end"), dict) else {}

    res_name_a = str(
        node_a.get("label_comp_id")
        or node_a.get("auth_comp_id")
        or residue_a.get("resName")
        or ""
    ).strip().upper()
    res_name_b = str(
        node_b.get("label_comp_id")
        or node_b.get("auth_comp_id")
        or residue_b.get("resName")
        or ""
    ).strip().upper()
    atom_name_a = _primary_contact_atom_name(
        node_a.get("auth_atom_id") or node_a.get("label_atom_id") or residue_a.get("atom")
    )
    atom_name_b = _primary_contact_atom_name(
        node_b.get("auth_atom_id") or node_b.get("label_atom_id") or residue_b.get("atom")
    )
    element_a = str(
        node_a.get("type_symbol")
        or residue_a.get("element")
        or guess_element(atom_name_a)
        or ""
    ).strip().upper()
    element_b = str(
        node_b.get("type_symbol")
        or residue_b.get("element")
        or guess_element(atom_name_b)
        or ""
    ).strip().upper()

    residue_pair_key = _unordered_residue_pair_key(residue_a, residue_b, prefix="preclass_res:")
    atom_pair_key = _unordered_residue_atom_pair_key(
        residue_a,
        atom_name_a,
        residue_b,
        atom_name_b,
        prefix="preclass_atom:",
    )
    feature_a = _contact_atom_feature_signature(
        res_name=res_name_a,
        atom_name=atom_name_a,
        element=element_a,
    )
    feature_b = _contact_atom_feature_signature(
        res_name=res_name_b,
        atom_name=atom_name_b,
        element=element_b,
    )
    feature_pair_key = _build_unordered_pair_key(feature_a, feature_b)
    interaction_type = str(raw.get("type") or "").strip().lower()
    motif_terms = "|".join(sorted(_normalize_arpeggio_contact_terms(raw.get("contact"))))
    altloc_a = _extract_altloc_family_from_node(node_a)
    altloc_b = _extract_altloc_family_from_node(node_b)
    if altloc_a and altloc_b:
        if altloc_a == altloc_b:
            altloc_class = f"same:{altloc_a}"
        else:
            altloc_pair = sorted([altloc_a, altloc_b])
            altloc_class = f"incompatible:{altloc_pair[0]}|{altloc_pair[1]}"
    elif altloc_a or altloc_b:
        altloc_class = f"mixed:{altloc_a or altloc_b}"
    else:
        altloc_class = "none"
    symmetry_class = (
        "non_identity"
        if _contact_has_non_identity_symmetry(raw, node_a, node_b)
        else "identity"
    )
    role_candidates = _resolve_hbond_role_candidates(
        res_name_a=res_name_a,
        atom_name_a=atom_name_a,
        element_a=element_a,
        res_name_b=res_name_b,
        atom_name_b=atom_name_b,
        element_b=element_b,
    )
    if len(role_candidates) == 1:
        donor_acceptor_direction = f"{role_candidates[0][0]}>{role_candidates[0][1]}"
    elif len(role_candidates) > 1:
        donor_acceptor_direction = "both"
    else:
        donor_acceptor_direction = "none"
    ring_key_a = _fallback_ring_site_key(residue_a, atom_name_a)
    ring_key_b = _fallback_ring_site_key(residue_b, atom_name_b)
    ring_pair_key = _build_unordered_pair_key(ring_key_a, ring_key_b)
    return (
        f"{residue_pair_key}|{atom_pair_key}|"
        f"feature:{feature_pair_key}|type:{interaction_type}|terms:{motif_terms}|"
        f"symmetry:{symmetry_class}|altloc:{altloc_class}|"
        f"ring:{ring_pair_key}|donor_acceptor:{donor_acceptor_direction}"
    )


def _compute_base_pair_pair_stats(
    prepared_contacts: List[Tuple[dict, dict, dict]],
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]] = None,
) -> Dict[str, dict]:
    pair_stats: Dict[str, dict] = {}
    best_partner_by_residue: Dict[str, Tuple[str, Tuple[float, float, float, float, float, float, float]]] = {}

    def _update_best_partner(
        residue_token: str,
        partner_token: str,
        rank_tuple: Tuple[float, float, float, float, float, float, float],
    ) -> None:
        existing = best_partner_by_residue.get(residue_token)
        if existing is None or rank_tuple > existing[1]:
            best_partner_by_residue[residue_token] = (partner_token, rank_tuple)

    for contact_rank, (raw, residue_a, residue_b) in enumerate(prepared_contacts):
        if not isinstance(raw, dict) or not isinstance(residue_a, dict) or not isinstance(residue_b, dict):
            continue
        node_a = raw.get("bgn") if isinstance(raw.get("bgn"), dict) else {}
        node_b = raw.get("end") if isinstance(raw.get("end"), dict) else {}
        res_name_a = str(
            node_a.get("label_comp_id")
            or node_a.get("auth_comp_id")
            or residue_a.get("resName")
            or ""
        ).strip().upper()
        res_name_b = str(
            node_b.get("label_comp_id")
            or node_b.get("auth_comp_id")
            or residue_b.get("resName")
            or ""
        ).strip().upper()
        atom_name_a = _primary_contact_atom_name(
            node_a.get("auth_atom_id") or node_a.get("label_atom_id") or residue_a.get("atom")
        )
        atom_name_b = _primary_contact_atom_name(
            node_b.get("auth_atom_id") or node_b.get("label_atom_id") or residue_b.get("atom")
        )
        if not atom_name_a or not atom_name_b:
            continue
        base_family_a = _nucleic_base_family(res_name_a)
        base_family_b = _nucleic_base_family(res_name_b)
        if not base_family_a or not base_family_b:
            continue
        # Same-residue contacts are never valid base-pair support evidence.
        chain_token_a = str(residue_a.get("chain") or "").strip()
        chain_token_b = str(residue_b.get("chain") or "").strip()
        seq_token_a = str(residue_a.get("seq") or "").strip()
        seq_token_b = str(residue_b.get("seq") or "").strip()
        if (
            chain_token_a
            and chain_token_b
            and seq_token_a
            and seq_token_b
            and chain_token_a == chain_token_b
            and seq_token_a == seq_token_b
        ):
            continue
        # Sequence-adjacent nucleotides are linkage neighbors, not base-pair partners.
        if _is_sequence_adjacent_nucleotide_pair(
            residue_a=residue_a,
            residue_b=residue_b,
            base_family_a=base_family_a,
            base_family_b=base_family_b,
        ):
            continue
        if not (_is_nucleobase_atom(res_name_a, atom_name_a) and _is_nucleobase_atom(res_name_b, atom_name_b)):
            continue
        if _is_nucleobase_glycosidic_atom(res_name_a, atom_name_a) or _is_nucleobase_glycosidic_atom(
            res_name_b, atom_name_b
        ):
            continue
        if not (
            _is_nucleobase_pairing_edge_atom(res_name_a, atom_name_a)
            and _is_nucleobase_pairing_edge_atom(res_name_b, atom_name_b)
        ):
            continue
        element_a = str(
            node_a.get("type_symbol")
            or residue_a.get("element")
            or guess_element(atom_name_a)
            or ""
        ).strip().upper()
        element_b = str(
            node_b.get("type_symbol")
            or residue_b.get("element")
            or guess_element(atom_name_b)
            or ""
        ).strip().upper()
        if element_a not in POLAR_CONTACT_ELEMENTS or element_b not in POLAR_CONTACT_ELEMENTS:
            continue
        donor_acceptor_consistent = _is_hbond_donor_acceptor_pair(
            res_name_a=res_name_a,
            atom_name_a=atom_name_a,
            element_a=element_a,
            res_name_b=res_name_b,
            atom_name_b=atom_name_b,
            element_b=element_b,
        )
        if not donor_acceptor_consistent:
            continue
        distance, _ = _resolve_contact_distance_value(
            raw,
            residue_a,
            residue_b,
            residue_atoms_index,
        )
        if distance is None or distance > BASE_PAIR_CANDIDATE_MAX_DISTANCE:
            continue
        preclassified = _resolve_impossible_contact_preclassification(
            distance=distance,
            element_a=element_a,
            element_b=element_b,
            suspect_invalid_mapping=False,
        )
        if isinstance(preclassified, dict):
            continue

        residue_pair_key = _unordered_residue_pair_key(
            residue_a,
            residue_b,
            prefix="basepair_support:",
        )
        residue_token_a = _residue_pair_identity_token(residue_a)
        residue_token_b = _residue_pair_identity_token(residue_b)
        atom_pair_key = _unordered_residue_atom_pair_key(
            residue_a,
            atom_name_a,
            residue_b,
            atom_name_b,
            prefix="basepair_atom:",
        )
        pair_stat = pair_stats.get(residue_pair_key)
        if pair_stat is None:
            pair_stat = {
                "residueTokenA": residue_token_a,
                "residueTokenB": residue_token_b,
                "sampleResidueA": dict(residue_a),
                "sampleResidueB": dict(residue_b),
                "supportAtomPairs": set(),
                "supportCount": 0,
                "candidateCount": 0,
                "bestDistance": None,
                "angleEvaluatedCount": 0,
                "anglePassedCount": 0,
                "angleStrongCount": 0,
                "canonicalTemplateMatches": 0,
                "coplanaritySupported": False,
                "mutualBestMatch": False,
            }
            pair_stats[residue_pair_key] = pair_stat

        support_atom_pairs = pair_stat.get("supportAtomPairs")
        if not isinstance(support_atom_pairs, set):
            support_atom_pairs = set()
            pair_stat["supportAtomPairs"] = support_atom_pairs
        support_atom_pairs.add(atom_pair_key)
        pair_stat["supportCount"] = len(support_atom_pairs)
        pair_stat["candidateCount"] = int(pair_stat.get("candidateCount") or 0) + 1

        best_distance = _coerce_float(pair_stat.get("bestDistance"))
        if best_distance is None or distance < best_distance:
            pair_stat["bestDistance"] = distance

        canonical_template_match = _is_canonical_base_pair_hbond_pair(
            res_name_a=res_name_a,
            atom_name_a=atom_name_a,
            element_a=element_a,
            res_name_b=res_name_b,
            atom_name_b=atom_name_b,
            element_b=element_b,
        )
        if canonical_template_match:
            pair_stat["canonicalTemplateMatches"] = int(pair_stat.get("canonicalTemplateMatches") or 0) + 1

        angle_value = _extract_arpeggio_hbond_angle(raw)
        if angle_value is None:
            proxy_angle, _ = _compute_hbond_proxy_angle(
                residue_a=residue_a,
                residue_b=residue_b,
                res_name_a=res_name_a,
                atom_name_a=atom_name_a,
                element_a=element_a,
                res_name_b=res_name_b,
                atom_name_b=atom_name_b,
                element_b=element_b,
                residue_atoms_index=residue_atoms_index,
            )
            angle_value = proxy_angle
        if angle_value is not None:
            pair_stat["angleEvaluatedCount"] = int(pair_stat.get("angleEvaluatedCount") or 0) + 1
            if angle_value >= HBOND_PROXY_ANGLE_MIN:
                pair_stat["anglePassedCount"] = int(pair_stat.get("anglePassedCount") or 0) + 1
            if angle_value >= HBOND_PROXY_STRONG_ANGLE_MIN:
                pair_stat["angleStrongCount"] = int(pair_stat.get("angleStrongCount") or 0) + 1

    for pair_stat in pair_stats.values():
        sample_residue_a = pair_stat.get("sampleResidueA")
        sample_residue_b = pair_stat.get("sampleResidueB")
        ring_geometry = {}
        if isinstance(sample_residue_a, dict) and isinstance(sample_residue_b, dict):
            ring_geometry = _compute_ring_geometry_metrics(
                sample_residue_a,
                sample_residue_b,
                residue_atoms_index,
            )
        ring_normal_angle = _coerce_float(ring_geometry.get("ring_normal_angle"))
        ring_interplanar_distance = _coerce_float(ring_geometry.get("ring_interplanar_distance"))
        ring_lateral_offset = _coerce_float(ring_geometry.get("ring_lateral_offset"))
        coplanarity_supported = _base_pair_coplanarity_supported_from_metrics(
            ring_normal_angle,
            ring_interplanar_distance,
            ring_lateral_offset,
        )
        pair_stat["coplanaritySupported"] = coplanarity_supported
        pair_stat["ringNormalAngle"] = ring_normal_angle
        pair_stat["ringInterplanarDistance"] = ring_interplanar_distance
        pair_stat["ringLateralOffset"] = ring_lateral_offset

        support_count = int(pair_stat.get("supportCount") or 0)
        best_distance = _coerce_float(pair_stat.get("bestDistance"))
        angle_evaluated_count = int(pair_stat.get("angleEvaluatedCount") or 0)
        angle_passed_count = int(pair_stat.get("anglePassedCount") or 0)
        angle_strong_count = int(pair_stat.get("angleStrongCount") or 0)
        score_components = {
            "support_count": _base_pair_support_score_component(support_count),
            "best_distance": _base_pair_distance_score_component(best_distance),
            "angles": _base_pair_angle_score_component(
                angle_evaluated_count,
                angle_passed_count,
                angle_strong_count,
            ),
            "mutual_best_match": 0.0,
            "coplanarity": 1.0 if coplanarity_supported else 0.0,
        }
        pair_stat["scoreComponents"] = score_components
        pair_stat["scoreTotal"] = _base_pair_score_total(score_components)

    for pair_stat in pair_stats.values():
        residue_token_a = str(pair_stat.get("residueTokenA") or "").strip()
        residue_token_b = str(pair_stat.get("residueTokenB") or "").strip()
        if not residue_token_a or not residue_token_b:
            continue
        rank_tuple = _base_pair_pair_rank_tuple(pair_stat)
        _update_best_partner(residue_token_a, residue_token_b, rank_tuple)
        _update_best_partner(residue_token_b, residue_token_a, rank_tuple)

    for pair_stat in pair_stats.values():
        residue_token_a = str(pair_stat.get("residueTokenA") or "").strip()
        residue_token_b = str(pair_stat.get("residueTokenB") or "").strip()
        best_a = best_partner_by_residue.get(residue_token_a)
        best_b = best_partner_by_residue.get(residue_token_b)
        mutual_best_match = bool(
            best_a is not None
            and best_b is not None
            and best_a[0] == residue_token_b
            and best_b[0] == residue_token_a
        )
        pair_stat["mutualBestMatch"] = mutual_best_match
        components = dict(pair_stat.get("scoreComponents") or {})
        components["mutual_best_match"] = 1.0 if mutual_best_match else 0.0
        pair_stat["scoreComponents"] = components
        pair_stat["scoreTotal"] = _base_pair_score_total(components)

    return pair_stats


def _contact_distance_for_preference(contact: dict) -> float:
    if not isinstance(contact, dict):
        return math.inf
    numeric = _coerce_float(contact.get("distance"))
    if numeric is None:
        return math.inf
    return numeric


def _prefer_contact_by_shorter_distance(candidate: dict, current: dict) -> bool:
    if not isinstance(candidate, dict):
        return False
    if not isinstance(current, dict):
        return True
    cand_distance = _contact_distance_for_preference(candidate)
    curr_distance = _contact_distance_for_preference(current)
    if cand_distance + 1e-6 < curr_distance:
        return True
    if cand_distance - 1e-6 > curr_distance:
        return False
    cand_atom_a = str(candidate.get("residueA", {}).get("atom") or "").strip()
    cand_atom_b = str(candidate.get("residueB", {}).get("atom") or "").strip()
    curr_atom_a = str(current.get("residueA", {}).get("atom") or "").strip()
    curr_atom_b = str(current.get("residueB", {}).get("atom") or "").strip()
    cand_specificity = int(bool(cand_atom_a)) + int(bool(cand_atom_b))
    curr_specificity = int(bool(curr_atom_a)) + int(bool(curr_atom_b))
    return cand_specificity > curr_specificity


def _sanitize_ring_label_token(label: str) -> str:
    token = re.sub(r"[^A-Za-z0-9_]+", "_", str(label or "").strip().lower()).strip("_")
    return token or "ring"


def _fallback_ring_site_key(residue: dict, atom_name_hint: str = "") -> str:
    residue_token = _residue_pair_identity_token(residue)
    residue_name = str(residue.get("resName") or "").strip().upper() or "UNK"
    atom_name = _primary_contact_atom_name(atom_name_hint or residue.get("atom")) or "?"
    return f"{residue_token}:{residue_name}:ring_atom:{atom_name}"


def _ring_site_key_from_descriptor(residue: dict, descriptor: Optional[dict], atom_name_hint: str = "") -> str:
    if not isinstance(descriptor, dict):
        return _fallback_ring_site_key(residue, atom_name_hint)
    residue_token = _residue_pair_identity_token(residue)
    residue_name = str(residue.get("resName") or "").strip().upper() or "UNK"
    label = _sanitize_ring_label_token(str(descriptor.get("label") or "ring"))
    digest = str(descriptor.get("hash") or "").strip().lower() or "unknown"
    return f"{residue_token}:{residue_name}:{label}:{digest}"


def _resolve_aromatic_ring_site_keys(
    residue_a: dict,
    residue_b: dict,
    residue_atoms_index: Optional[Dict[Tuple[str, str], List[AtomRecord]]],
) -> Dict[str, str]:
    atom_hint_a = _primary_contact_atom_name(residue_a.get("atom")) if isinstance(residue_a, dict) else ""
    atom_hint_b = _primary_contact_atom_name(residue_b.get("atom")) if isinstance(residue_b, dict) else ""
    ring_pair = _select_ring_descriptor_pair_by_contact(
        residue_a,
        residue_b,
        residue_atoms_index,
        atom_name_a=atom_hint_a,
        atom_name_b=atom_hint_b,
    )
    def single_ring(residue, other):
        descriptors = _residue_ring_descriptors(residue, residue_atoms_index)
        if not descriptors:
            return None
        names = {_normalize_atom_name(name) for name in str(residue.get("atom") or "").split(",") if name}
        partner = _resolve_contact_atom_record_from_payload(other, residue_atoms_index)
        def rank(descriptor):
            misses = len(names - set(descriptor.get("atom_names", [])))
            centroid = descriptor.get("centroid", (0.0, 0.0, 0.0))
            distance = sum((centroid[i] - (partner.x, partner.y, partner.z)[i]) ** 2 for i in range(3)) if partner else 0.0
            return misses, distance, descriptor.get("hash", "")
        return min(descriptors, key=rank)
    descriptor_a = ring_pair[0] if ring_pair else single_ring(residue_a, residue_b)
    descriptor_b = ring_pair[1] if ring_pair else single_ring(residue_b, residue_a)
    ring_key_a = _ring_site_key_from_descriptor(residue_a, descriptor_a, atom_name_hint=atom_hint_a)
    ring_key_b = _ring_site_key_from_descriptor(residue_b, descriptor_b, atom_name_hint=atom_hint_b)
    ring_pair_key = _build_unordered_pair_key(ring_key_a, ring_key_b)
    return {
        "ringKeyA": ring_key_a,
        "ringKeyB": ring_key_b,
        "ringPairKey": ring_pair_key,
        "ringAtomNamesA": list((descriptor_a or {}).get("atom_names", [])),
        "ringAtomNamesB": list((descriptor_b or {}).get("atom_names", [])),
    }


def _aromatic_score_distance_to_target(
    value: Optional[float],
    *,
    target: float,
    tolerance: float,
    min_value: Optional[float] = None,
    max_value: Optional[float] = None,
) -> float:
    numeric = _coerce_float(value)
    if numeric is None:
        return 0.0
    if min_value is not None and numeric < min_value:
        return 0.0
    if max_value is not None and numeric > max_value:
        return 0.0
    window = max(1e-6, float(tolerance))
    return max(0.0, 1.0 - abs(numeric - float(target)) / window)


def _aromatic_record_confidence_bonus(record: dict) -> float:
    if not isinstance(record, dict):
        return 0.0
    asserted = record.get("asserted")
    if not isinstance(asserted, dict):
        return 0.0
    confidence = str(asserted.get("confidence") or "").strip().lower()
    if confidence == "high":
        return 0.15
    if confidence == "medium":
        return 0.08
    return 0.0


def _aromatic_record_ring_metric(record: dict, metric_key: str) -> Optional[float]:
    if not isinstance(record, dict):
        return None
    ring = record.get("ring")
    if isinstance(ring, dict):
        numeric = _coerce_float(ring.get(metric_key))
        if numeric is not None:
            return numeric
    asserted = record.get("asserted")
    if isinstance(asserted, dict):
        asserted_ring = asserted.get("ring")
        if isinstance(asserted_ring, dict):
            return _coerce_float(asserted_ring.get(metric_key))
    return None


def _aromatic_record_score(record: dict, family: str) -> float:
    family_key = str(family or "").strip().lower()
    centroid = _aromatic_record_ring_metric(record, "centroid_distance")
    min_atom = _aromatic_record_ring_metric(record, "min_atom_distance")
    interplanar = _aromatic_record_ring_metric(record, "interplanar_distance")
    lateral = _aromatic_record_ring_metric(record, "lateral_offset")
    normal = _aromatic_record_ring_metric(record, "normal_angle")
    confidence_bonus = _aromatic_record_confidence_bonus(record)

    if family_key == "pi_pi":
        centroid_score = _aromatic_score_distance_to_target(
            centroid,
            target=4.8,
            tolerance=1.9,
            min_value=PI_PI_MIN_CENTROID_DISTANCE,
            max_value=PI_PI_MAX_CENTROID_DISTANCE,
        )
        interplanar_score = _aromatic_score_distance_to_target(
            interplanar,
            target=3.4,
            tolerance=0.9,
            min_value=PI_PI_MIN_INTERPLANAR_DISTANCE,
            max_value=PI_PI_MAX_INTERPLANAR_DISTANCE,
        )
        lateral_score = 0.0
        if lateral is not None:
            lateral_score = max(0.0, 1.0 - max(0.0, lateral) / max(1e-6, PI_PI_MAX_LATERAL_OFFSET))
        angle_score = 0.0
        if normal is not None:
            stacked = max(0.0, 1.0 - normal / max(1e-6, PI_PI_STACKED_MAX_NORMAL_ANGLE))
            tshaped = max(0.0, 1.0 - abs(normal - 90.0) / 30.0)
            angle_score = max(stacked, tshaped)
        return (
            0.35 * centroid_score
            + 0.25 * interplanar_score
            + 0.2 * lateral_score
            + 0.2 * angle_score
            + confidence_bonus
        )

    if family_key == "pi_cation":
        centroid_score = _aromatic_score_distance_to_target(
            centroid,
            target=4.8,
            tolerance=2.1,
            min_value=PI_PI_MIN_CENTROID_DISTANCE,
            max_value=PI_PI_MAX_CENTROID_DISTANCE,
        )
        lateral_score = 0.0
        if lateral is not None:
            lateral_score = max(0.0, 1.0 - max(0.0, lateral) / max(1e-6, PI_PI_TSHAPED_MAX_LATERAL_OFFSET))
        return (0.65 * centroid_score + 0.35 * lateral_score + confidence_bonus)

    if family_key == "aromatic_packing":
        min_atom_score = _aromatic_score_distance_to_target(
            min_atom,
            target=3.6,
            tolerance=1.0,
            min_value=AROMATIC_PACKING_MIN_DISTANCE,
            max_value=AROMATIC_PACKING_MAX_DISTANCE,
        )
        return min_atom_score + confidence_bonus

    if family_key == "aromatic_proximal":
        centroid_score = _aromatic_score_distance_to_target(
            centroid,
            target=5.4,
            tolerance=2.6,
            min_value=AROMATIC_PACKING_MIN_DISTANCE,
            max_value=PROXIMAL_CONTACT_MAX_DISTANCE,
        )
        return 0.65 * centroid_score + confidence_bonus

    return confidence_bonus


def _aromatic_family_top_k(family: str) -> int:
    family_key = str(family or "").strip().lower()
    return max(1, int(AROMATIC_TOP_K_PER_RESIDUE_PAIR_BY_FAMILY.get(family_key, 1)))


def _prefer_aromatic_record(candidate: dict, current: dict, family: str) -> bool:
    if not isinstance(candidate, dict):
        return False
    if not isinstance(current, dict):
        return True
    family_key = str(family or "").strip().lower()
    candidate_score = _aromatic_record_score(candidate, family_key)
    current_score = _aromatic_record_score(current, family_key)
    if candidate_score > current_score + 1e-6:
        return True
    if candidate_score + 1e-6 < current_score:
        return False
    return _prefer_contact_by_shorter_distance(candidate, current)


def _select_aromatic_records_for_output(
    pending_records: List[dict],
) -> List[dict]:
    if not pending_records:
        return []
    grouped_by_residue_pair_and_family: Dict[str, List[dict]] = {}
    for entry in pending_records:
        if not isinstance(entry, dict):
            continue
        family = str(entry.get("family") or "").strip().lower()
        residue_a = entry.get("residueA")
        residue_b = entry.get("residueB")
        if family not in AROMATIC_ASSERTED_FAMILIES:
            continue
        if not isinstance(residue_a, dict) or not isinstance(residue_b, dict):
            continue
        residue_pair_key = _unordered_residue_pair_key(
            residue_a,
            residue_b,
            prefix=f"aromatic_cap:{family}:",
        )
        grouped_by_residue_pair_and_family.setdefault(residue_pair_key, []).append(entry)

    selected: List[dict] = []
    for entries in grouped_by_residue_pair_and_family.values():
        if not entries:
            continue
        family = str(entries[0].get("family") or "").strip().lower()
        deduped_by_site: Dict[str, dict] = {}
        for entry in entries:
            record = entry.get("record")
            if not isinstance(record, dict):
                continue
            asserted_payload = record.get("asserted")
            subtype = ""
            if isinstance(asserted_payload, dict):
                subtype = str(asserted_payload.get("subtype") or "").strip().lower()
            ring_pair_key = str(
                record.get("ringPairKey")
                or (asserted_payload.get("ringPairKey") if isinstance(asserted_payload, dict) else "")
                or record.get("pairKey")
                or ""
            ).strip()
            if not ring_pair_key:
                residue_a = entry.get("residueA") or {}
                residue_b = entry.get("residueB") or {}
                ring_pair_key = _unordered_residue_pair_key(
                    residue_a,
                    residue_b,
                    prefix="ringpair_fallback:",
                )
            dedupe_key = f"{family}:{subtype}:{ring_pair_key}"
            current = deduped_by_site.get(dedupe_key)
            if current is None:
                deduped_by_site[dedupe_key] = entry
                continue
            if _prefer_aromatic_record(record, current.get("record"), family):
                deduped_by_site[dedupe_key] = entry

        ranked = sorted(
            deduped_by_site.values(),
            key=lambda entry: (
                -_aromatic_record_score(entry.get("record") or {}, family),
                _contact_distance_for_preference(entry.get("record") or {}),
                -(int(bool(str((entry.get("record") or {}).get("residueA", {}).get("atom") or "").strip())) + int(bool(str((entry.get("record") or {}).get("residueB", {}).get("atom") or "").strip()))),
                int(entry.get("rank") or 0),
            ),
        )
        # Connector budgets belong to rendering; retain every distinct scientific ring site.
        selected.extend(ranked)
    return selected


def _contact_matches_chain_pair(
    residue_a: dict,
    residue_b: dict,
    chain_ids_a: Set[str],
    chain_ids_b: Set[str],
    intrachain: bool,
) -> bool:
    if not residue_a or not residue_b:
        return False
    contact_chain_a = str(residue_a.get("chain") or "").strip()
    contact_chain_b = str(residue_b.get("chain") or "").strip()
    if not contact_chain_a or not contact_chain_b:
        return False
    if intrachain:
        return contact_chain_a in chain_ids_a and contact_chain_b in chain_ids_a
    return (
        (contact_chain_a in chain_ids_a and contact_chain_b in chain_ids_b)
        or (contact_chain_a in chain_ids_b and contact_chain_b in chain_ids_a)
    )


def _new_per_residue_entry(residue: dict) -> dict:
    return {
        "chain": str(residue.get("chain") or "").strip(),
        "resName": str(residue.get("resName") or "").strip().upper(),
        "seq": str(residue.get("seq") or "").strip(),
        "hydrophobic": 0,
        "hbond": 0,
        "polar_contact": 0,
        "base_pairing": 0,
        "salt_bridge": 0,
        "halogen_bond": 0,
        "metal_coordination": 0,
        "pi_pi": 0,
        "pi_cation": 0,
        "aromatic_packing": 0,
        "vdw": 0,
        "clash": 0,
        "other": 0,
        "total": 0,
    }


def _bump_per_residue_from_payload(per_residue: Dict[str, dict], residue: dict, category: str) -> None:
    if not residue:
        return
    chain = str(residue.get("chain") or "").strip()
    seq = str(residue.get("seq") or "").strip()
    if not chain or not seq:
        return
    stat_key = CONTACT_CATEGORY_TO_PER_RESIDUE_KEY.get(category, "other")
    residue_key = f"{chain}:{seq}"
    entry = per_residue.setdefault(residue_key, _new_per_residue_entry(residue))
    if stat_key not in entry:
        stat_key = "other"
    entry[stat_key] += 1
    entry["total"] += 1


def _downgrade_confidence_level(level: str, steps: int = 1) -> str:
    order = ["low", "medium", "high"]
    token = str(level or "").strip().lower()
    if token not in order:
        token = "low"
    idx = order.index(token)
    next_idx = max(0, idx - max(1, int(steps)))
    return order[next_idx]


def _apply_atom_reuse_confidence_penalties(contacts: Dict[str, List[dict]]) -> None:
    if not isinstance(contacts, dict):
        return
    for bucket, threshold in ATOM_REUSE_CONFIDENCE_THRESHOLD_BY_BUCKET.items():
        rows = contacts.get(bucket)
        if not isinstance(rows, list) or not rows:
            continue
        atom_counts: Dict[str, int] = {}
        for record in rows:
            if not isinstance(record, dict):
                continue
            atom_key_a = str(record.get("atomKeyA") or "").strip()
            atom_key_b = str(record.get("atomKeyB") or "").strip()
            if atom_key_a:
                atom_counts[atom_key_a] = atom_counts.get(atom_key_a, 0) + 1
            if atom_key_b:
                atom_counts[atom_key_b] = atom_counts.get(atom_key_b, 0) + 1
        for record in rows:
            if not isinstance(record, dict):
                continue
            atom_key_a = str(record.get("atomKeyA") or "").strip()
            atom_key_b = str(record.get("atomKeyB") or "").strip()
            max_reuse = max(
                atom_counts.get(atom_key_a, 0) if atom_key_a else 0,
                atom_counts.get(atom_key_b, 0) if atom_key_b else 0,
            )
            if max_reuse <= threshold:
                continue
            asserted = record.get("asserted")
            if not isinstance(asserted, dict):
                continue
            current = str(asserted.get("confidence") or "low").strip().lower() or "low"
            penalty_steps = 2 if max_reuse >= threshold + 3 else 1
            asserted["confidence"] = _downgrade_confidence_level(current, penalty_steps)
            evidence = asserted.get("evidence")
            if not isinstance(evidence, list):
                evidence = []
                asserted["evidence"] = evidence
            _append_evidence(evidence, "atom_reuse_penalty")
            _append_evidence(evidence, f"atom_reuse_count_{max_reuse}")
            _append_evidence(evidence, f"atom_reuse_threshold_{threshold}")

    # Cross-family cleanup: if the same atom pair appears in multiple families,
    # downgrade weaker semantic families to reduce ambiguous table noise.
    family_strength_rank = {
        "metal_coordination": 0,
        "salt_bridge": 1,
        "halogen_bond": 2,
        "hbond": 3,
        "pi_cation": 4,
        "pi_pi": 5,
        "base_pairing": 6,
        "aromatic_packing": 7,
        "hydrophobic": 8,
        "polar_contact": 9,
        "packing_contact": 10,
        "proximal": 11,
        "polar_proximal": 12,
        "other": 13,
    }
    weaker_families = {
        "hydrophobic",
        "aromatic_packing",
        "polar_contact",
        "packing_contact",
        "proximal",
        "polar_proximal",
        "other",
    }
    pair_to_rows: Dict[str, List[dict]] = {}
    for bucket_rows in contacts.values():
        if not isinstance(bucket_rows, list):
            continue
        for record in bucket_rows:
            if not isinstance(record, dict):
                continue
            pair_key = str(record.get("pairKey") or "").strip()
            if not pair_key:
                atom_key_a = str(record.get("atomKeyA") or "").strip()
                atom_key_b = str(record.get("atomKeyB") or "").strip()
                pair_key = _build_unordered_pair_key(atom_key_a, atom_key_b)
            if not pair_key:
                continue
            pair_to_rows.setdefault(pair_key, []).append(record)
    for rows in pair_to_rows.values():
        family_rows: List[Tuple[str, dict]] = []
        for record in rows:
            if not isinstance(record, dict):
                continue
            family = str(record.get("type") or "").strip().lower()
            if not family:
                asserted = record.get("asserted")
                if isinstance(asserted, dict):
                    family = str(asserted.get("family") or "").strip().lower()
            if not family:
                family = "other"
            family_rows.append((family, record))
        unique_families = {family for family, _ in family_rows}
        if len(unique_families) <= 1:
            continue
        strongest_family = min(
            unique_families,
            key=lambda family: int(family_strength_rank.get(family, 99)),
        )
        strongest_rank = int(family_strength_rank.get(strongest_family, 99))
        for family, record in family_rows:
            if family == strongest_family:
                continue
            if family not in weaker_families:
                continue
            family_rank = int(family_strength_rank.get(family, 99))
            if family_rank <= strongest_rank:
                continue
            asserted = record.get("asserted")
            if not isinstance(asserted, dict):
                continue
            current = str(asserted.get("confidence") or "low").strip().lower() or "low"
            asserted["confidence"] = _downgrade_confidence_level(current, 1)
            evidence = asserted.get("evidence")
            if not isinstance(evidence, list):
                evidence = []
                asserted["evidence"] = evidence
            _append_evidence(evidence, "multi_family_competition_penalty")
            _append_evidence(evidence, f"competing_family_{strongest_family}")
            _append_evidence(evidence, f"suppressed_weaker_family_{family}")


def _arpeggio_contacts_cache_key(
    structure_text: str,
    structure_format: str,
    selection: List[str],
) -> str:
    selection_digest = hashlib.sha256("\n".join(selection).encode("utf-8")).hexdigest()[:16]
    return (
        f"{TOOL_VERSION}:{str(structure_format or '').strip().lower()}:"
        f"{_structure_digest(structure_text)}:{selection_digest}"
    )


def _convert_pdb_text_to_mmcif_text(pdb_text: str) -> str:
    if _gemmi is None:
        detail = str(GEMMI_IMPORT_ERROR) if GEMMI_IMPORT_ERROR else "unknown import error"
        raise RuntimeError(f"gemmi unavailable for PDB -> mmCIF conversion ({detail})")
    handle = tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False, encoding="utf-8")
    temp_path = handle.name
    try:
        handle.write(pdb_text)
        handle.flush()
        handle.close()
        structure = _gemmi.read_structure(temp_path)
        mmcif_doc = structure.make_mmcif_document()
        return mmcif_doc.as_string()
    finally:
        try:
            os.unlink(temp_path)
        except Exception:
            pass


def _filter_arpeggio_coordinate_rows(block):
    category = block.get_mmcif_category("_atom_site.")
    columns = {key.lower(): key for key in category}
    if not category:
        return False
    def value(index, *names, default=""):
        for name in names:
            key = columns.get(name.lower())
            if key is not None:
                token = _clean_optional_token(category[key][index])
                if token:
                    return token
        return default
    records = []
    for index in range(len(next(iter(category.values())))):
        chain = value(index, "auth_asym_id", "label_asym_id", default="A")
        seq = value(index, "auth_seq_id", "label_seq_id") + value(index, "pdbx_PDB_ins_code")
        records.append(AtomRecord(chain, chain, value(index, "label_asym_id", default=chain),
            value(index, "auth_comp_id", "label_comp_id"), seq,
            value(index, "auth_atom_id", "label_atom_id"), value(index, "type_symbol").upper(),
            float(index), 0.0, 0.0, value(index, "label_alt_id", "pdbx_PDB_alt_id"),
            _coordinate_occupancy(value(index, "occupancy")), value(index, "pdbx_PDB_model_num", default="1")))
    selected = _select_model_conformers(records)
    indexes = sorted(int(atom.x) for atom in selected)
    if len(indexes) == len(records):
        return False
    block.set_mmcif_category("_atom_site.", {key: [values[i] for i in indexes] for key, values in category.items()})
    return True


def _prepare_arpeggio_mmcif_text(structure_text: str) -> str:
    """Complete Arpeggio's component metadata in its temporary input copy.

    Model-building exports can omit the optional _chem_comp.type/name columns
    (or the whole category), but Arpeggio indexes them unconditionally. Preserve
    supplied classifications, coordinates and identifiers. Unknown components
    retain an unknown type rather than acquiring an invented chemical class.
    """
    if _gemmi is None:
        return structure_text
    document = _gemmi.cif.read_string(structure_text)
    block = document.sole_block()
    raw_components = block.get_mmcif_category("_chem_comp.")
    components = {key.lower(): values for key, values in raw_components.items()}
    changed = _filter_arpeggio_coordinate_rows(block)
    changed = changed or components.keys() != raw_components.keys()
    atom_components = [_gemmi.cif.as_string(value) for value in block.find_values("_atom_site.label_comp_id")]
    component_ids = list(components.get("id", []))
    if components and not component_ids:
        raise ValueError("The mmCIF _chem_comp category has no component identifiers")

    # An explicit polymer entity can identify an otherwise unfamiliar modified
    # residue. Require one consistent entity class before using that evidence.
    entity_poly = {key.lower(): values for key, values in block.get_mmcif_category("_entity_poly.").items()}
    entity_types = {}
    for entity_id, polymer_type in zip(entity_poly.get("entity_id", []), entity_poly.get("type", [])):
        token = str(polymer_type or "").lower()
        if token.startswith("polypeptide"):
            entity_types[entity_id] = "peptide linking"
        elif token == "polydeoxyribonucleotide":
            entity_types[entity_id] = "DNA linking"
        elif token == "polyribonucleotide":
            entity_types[entity_id] = "RNA linking"
        elif token.startswith("polysaccharide"):
            entity_types[entity_id] = "saccharide"
    component_entity_types: Dict[str, Set[Optional[str]]] = {}
    atom_entities = [_gemmi.cif.as_string(value) for value in block.find_values("_atom_site.label_entity_id")]
    for component_id, entity_id in zip(atom_components, atom_entities):
        component_entity_types.setdefault(component_id, set()).add(entity_types.get(entity_id))

    def inferred_type(component_id: str) -> Optional[str]:
        token = str(component_id).upper()
        if token in STANDARD_AMINO_RESIDUES:
            return "peptide linking"
        if token in {"DA", "DC", "DG", "DT", "DI", "DU"}:
            return "DNA linking"
        if token in {"A", "C", "G", "U", "I"}:
            return "RNA linking"
        if token in WATER_RESIDUES:
            return "non-polymer"
        candidates = component_entity_types.get(component_id, set())
        return next(iter(candidates)) if len(candidates) == 1 else None

    if "id" not in components:
        components["id"] = component_ids
    if "type" not in components:
        components["type"] = [inferred_type(component_id) for component_id in component_ids]
        changed = True
    if "name" not in components:
        components["name"] = [None] * len(component_ids)
        changed = True
    known_ids = set(component_ids)
    for component_id in atom_components:
        if not component_id or component_id in {".", "?"} or component_id in known_ids:
            continue
        known_ids.add(component_id)
        for key, values in components.items():
            value = None
            if key == "id":
                value = component_id
            elif key == "type":
                value = inferred_type(component_id)
            values.append(value)
        changed = True
    for index, component_id in enumerate(components["id"]):
        if not components["name"][index]:
            components["name"][index] = "WATER" if str(component_id).upper() in WATER_RESIDUES else component_id
            changed = True
    if not changed:
        return structure_text
    block.set_mmcif_category("_chem_comp.", components)
    return document.as_string()


def _perceive_partial_hydrogen_roles(interaction_complex):
    """Infer missing H valence without inventing underconstrained H positions.

    Open Babel's reader treats any explicit H as a fully hydrogenated input.
    Typical implicit valence restores donor typing on partially hydrogenated
    structures. Generating new coordinates here gives random water/hydroxyl
    orientations, so omitted H remain implicit and have no angle evidence.
    """
    if _openbabel is None or not getattr(interaction_complex, "input_has_hydrogens", False):
        return 0
    original = getattr(interaction_complex, "ob_mol", None)
    if original is None:
        return 0
    before = {atom.GetId(): (atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetX(), atom.GetY(), atom.GetZ())
              for atom in _openbabel.OBMolAtomIter(original)}
    perceived = _openbabel.OBMol(original)
    inferred = 0
    for atom in _openbabel.OBMolAtomIter(perceived):
        if atom.GetAtomicNum() not in {1, 0}:
            previous = atom.GetImplicitHCount()
            _openbabel.OBAtomAssignTypicalImplicitHydrogens(atom)
            inferred += max(0, atom.GetImplicitHCount() - previous)
    after = {atom.GetId(): (atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetX(), atom.GetY(), atom.GetZ())
             for atom in _openbabel.OBMolAtomIter(perceived)}
    if after != before:
        return 0
    if inferred:
        interaction_complex.ob_mol = perceived
        interaction_complex._roami_hydrogen_typing_source = "implicit_hydrogen_valence"
    return inferred


def _enrich_arpeggio_contacts(contacts, interaction_complex):
    """Export the actual engine chemical features, bond graph and hydrogen geometry.

    Coordinates/types are cached per engine atom; enrichment never protonates or
    modifies the engine molecule. Perceived charge is labelled as modelled.
    """
    atom_contacts = getattr(interaction_complex, "atom_contacts", ())
    payload_cache = {}
    hydrogen_cache = {}
    rings_by_atom = {}
    for ring in getattr(getattr(interaction_complex, 'biopython_str', None), 'rings', {}).values():
        ring_atoms = ring.get('atoms', [])
        names = sorted(str(atom.name) for atom in ring_atoms)
        if ring_atoms and all(atom.get_parent() is ring_atoms[0].get_parent() for atom in ring_atoms):
            for atom in ring_atoms:
                rings_by_atom.setdefault(id(atom), []).append(names)
    def atom_payload(atom):
        key = id(atom)
        if key in payload_cache:
            return payload_cache[key]
        payload = {"type_symbol": str(atom.element).upper(), "auth_atom_id": str(atom.name),
                   "coordinates": [round(float(v), 6) for v in atom.coord],
                   "label_alt_id": str(atom.get_altloc()).strip()}
        if hasattr(atom, "atom_types"):
            payload["atom_types"] = sorted(atom.atom_types)
        if hasattr(atom, "formal_charge"):
            payload["formal_charge"] = int(atom.formal_charge)
            payload["formal_charge_source"] = "arpeggio_openbabel_formal_charge"
        if hasattr(atom, "vdw_radius"):
            payload["vdw_radius"] = float(atom.vdw_radius)
        if key in rings_by_atom:
            payload['aromatic_rings'] = sorted(rings_by_atom[key])
        payload_cache[key] = payload
        hydrogens, anchors = [], []
        molecule = getattr(interaction_complex, "ob_mol", None)
        ob_id = getattr(interaction_complex, "bio_to_ob", {}).get(atom)
        if molecule is not None and ob_id is not None and _openbabel is not None:
            ob_atom = molecule.GetAtomById(ob_id)
            for neighbor in _openbabel.OBAtomAtomIter(ob_atom):
                bio_atom = getattr(interaction_complex, "ob_to_bio", {}).get(neighbor.GetId())
                coords = [round(float(v), 6) for v in (neighbor.GetX(), neighbor.GetY(), neighbor.GetZ())]
                if neighbor.GetAtomicNum() == 1:
                    hydrogens.append({"atomName": str(bio_atom.name) if bio_atom is not None else None,
                                      "coordinates": coords, "bondSource": "arpeggio_bond_graph"})
                elif bio_atom is not None:
                    bond = ob_atom.GetBond(neighbor)
                    anchors.append({"auth_atom_id": str(bio_atom.name), "type_symbol": str(bio_atom.element).upper(),
                                    "coordinates": coords, "label_alt_id": str(bio_atom.get_altloc()).strip(),
                                    "same_residue": bio_atom.get_parent() is atom.get_parent(),
                                    "bond_order": int(bond.GetBondOrder()) if bond else None,
                                    "aromatic": bool(bond.IsAromatic()) if bond else False})
        if anchors:
            payload['bonded_atoms'] = sorted(anchors, key=lambda item: item['auth_atom_id'])
        hydrogen_cache[key] = hydrogens
        return payload
    for row, contact in zip(contacts, atom_contacts):
        if row.get("type") != "atom-atom":
            break
        a, b = getattr(contact, "bgn_atom", None), getattr(contact, "end_atom", None)
        if a is None or b is None:
            continue
        if any(str(row.get(side, {}).get("auth_atom_id")) != str(atom.name) for side, atom in (("bgn", a), ("end", b))):
            continue
        for side, atom in (("bgn", a), ("end", b)):
            row[side].update(atom_payload(atom))
        if getattr(interaction_complex, "_roami_hydrogen_typing_source", None):
            row["hydrogen_typing_source"] = interaction_complex._roami_hydrogen_typing_source
        candidates = []
        for donor_side, donor, acceptor in (("A", a, b), ("B", b, a)):
            if "hbond donor" not in getattr(donor, "atom_types", ()) or "hbond acceptor" not in getattr(acceptor, "atom_types", ()):
                continue
            for hydrogen in getattr(donor, "h_coords", ()):
                hd = tuple(float(donor.coord[i] - hydrogen[i]) for i in range(3))
                ha = tuple(float(acceptor.coord[i] - hydrogen[i]) for i in range(3))
                h_distance = math.sqrt(sum(component * component for component in ha))
                maximum = 1.2 + float(getattr(acceptor, "vdw_radius", 1.7)) + ARPEGGIO_VDW_COMP
                if h_distance > maximum:
                    continue
                angle = _angle_between_vectors_degrees(hd, ha)
                if angle is None:
                    continue
                actual = next((h for h in hydrogen_cache[id(donor)] if all(abs(h['coordinates'][i]-float(hydrogen[i])) < 1e-4 for i in range(3))), None)
                h_payload = dict(actual) if actual else {"atomName": None, "coordinates": [round(float(v),6) for v in hydrogen]}
                source = "input_hydrogens" if actual and actual.get('atomName') else "arpeggio_generated_hydrogens"
                # Some test/adaptor engines expose H coordinates but no bond graph.
                if not actual and getattr(interaction_complex, "input_has_hydrogens", False):
                    source = 'input_hydrogens'
                h_payload.update(element='H', source=source, inferred=source != 'input_hydrogens')
                candidates.append((angle, donor_side, h_distance, h_payload))
        if candidates:
            angle, donor_side, h_distance, hydrogen = max(candidates, key=lambda c: (c[0], c[1], str(c[3]['coordinates'])))
            row["hbond_angle"] = round(angle, 3)
            row["hbond_donor_side"] = donor_side
            row["hbond_h_acceptor_distance"] = round(h_distance, 3)
            row["hbond_geometry_source"] = hydrogen['source']
            row['hbond_hydrogen'] = hydrogen
    return contacts


def _run_arpeggio_contacts(structure_text: str, structure_format: str, selection: List[str]) -> List[dict]:
    if InteractionComplex is None:
        detail = str(ARPEGGIO_IMPORT_ERROR) if ARPEGGIO_IMPORT_ERROR else "unknown import error"
        raise ValueError(
            "PDBe Arpeggio runtime unavailable. "
            "Ensure 'pdbe-arpeggio' and its Python dependencies are installed "
            f"(current import error: {detail})."
        )
    contacts_cache_key = _arpeggio_contacts_cache_key(structure_text, structure_format, selection)
    cached_contacts = ARPEGGIO_CONTACTS_CACHE.get(contacts_cache_key)
    if isinstance(cached_contacts, list):
        return cached_contacts

    def _run_once(text_payload: str, fmt: str) -> List[dict]:
        if str(fmt or "").strip().lower() == "pdb":
            text_payload = _convert_pdb_text_to_mmcif_text(text_payload)
        suffix = ".cif"
        text_payload = _prepare_arpeggio_mmcif_text(text_payload)
        handle = tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False, encoding="utf-8")
        path = handle.name
        try:
            handle.write(text_payload)
            handle.flush()
            handle.close()

            interaction_complex = InteractionComplex(path)
            _perceive_partial_hydrogen_roles(interaction_complex)
            for method_name in ("structure_checks", "address_ambiguities", "initialize"):
                method = getattr(interaction_complex, method_name, None)
                if callable(method):
                    method()
            run_arpeggio = getattr(interaction_complex, "run_arpeggio", None)
            if not callable(run_arpeggio):
                raise RuntimeError("PDBe Arpeggio runtime does not expose run_arpeggio().")
            try:
                run_arpeggio(
                    selection,
                    interacting_cutoff=ARPEGGIO_INTERACTING_CUTOFF,
                    vdw_comp=ARPEGGIO_VDW_COMP,
                    include_sequence_adjacent=False,
                )
            except TypeError:
                # Compatibility path for older signatures.
                run_arpeggio(selection)

            getter = getattr(interaction_complex, "get_contacts", None)
            if not callable(getter):
                raise RuntimeError("PDBe Arpeggio runtime does not expose get_contacts().")
            contacts = getter()
            if not isinstance(contacts, list):
                return []
            return _enrich_arpeggio_contacts([row for row in contacts if isinstance(row, dict)], interaction_complex)
        finally:
            try:
                os.unlink(path)
            except Exception:
                pass

    try:
        normalized = _run_once(structure_text, structure_format)
        ARPEGGIO_CONTACTS_CACHE.set(contacts_cache_key, normalized)
        return normalized
    except Exception as exc:
        format_token = str(structure_format or "").strip().lower()
        message = str(exc).lower()
        should_retry_as_mmcif = (
            format_token == "pdb" and
            ("expected block header" in message or "data_" in message)
        )
        if should_retry_as_mmcif:
            try:
                mmcif_text = _convert_pdb_text_to_mmcif_text(structure_text)
                normalized = _run_once(mmcif_text, "mmcif")
                ARPEGGIO_CONTACTS_CACHE.set(contacts_cache_key, normalized)
                return normalized
            except Exception as retry_exc:
                raise RuntimeError(
                    "PDBe Arpeggio failed for PDB input; mmCIF fallback retry also failed: "
                    f"{retry_exc}"
                ) from retry_exc
        if isinstance(exc, ValueError):
            raise
        raise RuntimeError(f"PDBe Arpeggio failed: {exc}") from exc


def analyze_interface(
    structure_text: str,
    chain_a: str,
    chain_b: str,
    mode: str = "all",
    structure_format: str = "mmcif",
    focus_residue: Optional[str] = None,
) -> dict:
    fmt = (structure_format or "mmcif").lower()
    atoms, aliases = _parse_structure_cached(structure_text, fmt)
    parsed_focus = _parse_focus_residue_key(focus_residue)
    requested_chain_a = str(chain_a or "").strip()
    requested_chain_b = str(chain_b or "").strip()
    resolved_chain_ids_a = aliases.resolve(requested_chain_a)
    resolved_chain_ids_b = aliases.resolve(requested_chain_b)
    chain_a = aliases.normalize(requested_chain_a)
    chain_b = aliases.normalize(requested_chain_b)
    focus_candidates: List[Tuple[str, str]] = []
    if parsed_focus:
        focus_chain_raw, focus_seq = parsed_focus
        seen_focus_candidates: Set[Tuple[str, str]] = set()
        for focus_chain_resolved in sorted(aliases.resolve(focus_chain_raw)):
            key = (focus_chain_resolved, focus_seq)
            if key in seen_focus_candidates:
                continue
            seen_focus_candidates.add(key)
            focus_candidates.append(key)
        focus_chain_normalized = aliases.normalize(focus_chain_raw)
        normalized_key = (focus_chain_normalized, focus_seq)
        if focus_chain_normalized and normalized_key not in seen_focus_candidates:
            focus_candidates.append(normalized_key)
    intrachain = bool(
        resolved_chain_ids_a
        and resolved_chain_ids_b
        and resolved_chain_ids_a == resolved_chain_ids_b
    )

    atoms_a = [atom for atom in atoms if atom.chain_id in resolved_chain_ids_a]
    atoms_b = [atom for atom in atoms if atom.chain_id in resolved_chain_ids_b]
    if not atoms_a or not atoms_b:
        raise ValueError("No atoms found for one or both chains")

    contacts = {bucket: [] for bucket in CONTACT_BUCKET_TO_CATEGORY}
    per_residue: Dict[str, dict] = {}
    hydrophobic_pair_to_index: Dict[str, int] = {}
    salt_pair_to_index: Dict[str, int] = {}
    metal_pair_to_index: Dict[str, int] = {}
    metal_pair_to_donor_atoms: Dict[str, Set[str]] = {}
    hbond_carboxylate_pair_to_index: Dict[str, int] = {}
    hbond_carboxylate_pair_to_acceptors: Dict[str, Set[str]] = {}
    aromatic_pending_records: List[dict] = []
    residue_atoms_index = ResidueAtomIndex(atoms)

    selection, applied_focuses = _build_arpeggio_selection(
        atoms,
        resolved_chain_ids_a,
        resolved_chain_ids_b,
        focus_residue_candidates=focus_candidates,
    )
    raw_contacts = _run_arpeggio_contacts(structure_text, fmt, selection)
    focus_match_keys: Set[Tuple[str, str]] = set()
    for focus_chain, focus_seq in applied_focuses:
        chain_token = str(focus_chain or "").strip()
        seq_token = str(focus_seq or "").strip()
        if not chain_token or not seq_token:
            continue
        focus_match_keys.add((chain_token, seq_token))
        for normalized_chain in sorted(aliases.resolve(chain_token)):
            if normalized_chain and normalized_chain != chain_token:
                focus_match_keys.add((normalized_chain, seq_token))

    prepared_contact_candidates: Dict[str, Tuple[float, int, Tuple[dict, dict, dict]]] = {}
    for raw_rank, raw in enumerate(raw_contacts):
        node_a = raw.get("bgn")
        node_b = raw.get("end")
        residue_a = _build_residue_payload_from_arpeggio_partner(node_a, aliases)
        residue_b = _build_residue_payload_from_arpeggio_partner(node_b, aliases)
        if not residue_a or not residue_b:
            continue
        residue_a = _hydrate_contact_partner(residue_a, node_a, residue_atoms_index)
        residue_b = _hydrate_contact_partner(residue_b, node_b, residue_atoms_index)
        for residue, node in ((residue_a, node_a), (residue_b, node_b)):
            if isinstance(node, dict):
                residue_atoms_index.chemical_nodes[(residue['chain'], residue['seq'], _primary_contact_atom_name(residue.get('atom')))] = node
        if not _contact_matches_chain_pair(
            residue_a,
            residue_b,
            resolved_chain_ids_a,
            resolved_chain_ids_b,
            intrachain,
        ):
            continue
        if focus_match_keys:
            residue_key_a = (
                str(residue_a.get("chain") or "").strip(),
                str(residue_a.get("seq") or "").strip(),
            )
            residue_key_b = (
                str(residue_b.get("chain") or "").strip(),
                str(residue_b.get("seq") or "").strip(),
            )
            focus_hit_a = residue_key_a in focus_match_keys
            focus_hit_b = residue_key_b in focus_match_keys
            if not (focus_hit_a or focus_hit_b):
                continue
        for residue, node in ((residue_a, node_a), (residue_b, node_b)):
            if isinstance(node, dict):
                residue_atoms_index.chemical_nodes[(residue['chain'], residue['seq'], _primary_contact_atom_name(residue.get('atom')))] = node
        # Preserve every distinct observation until roles/sites/H identities exist.
        duplicate_key = f"raw:{raw_rank}"
        candidate_distance, _ = _resolve_contact_distance_value(
            raw,
            residue_a,
            residue_b,
            residue_atoms_index,
        )
        candidate_distance_for_rank = candidate_distance if candidate_distance is not None else math.inf
        existing = prepared_contact_candidates.get(duplicate_key)
        if existing is None:
            prepared_contact_candidates[duplicate_key] = (
                candidate_distance_for_rank,
                raw_rank,
                (raw, residue_a, residue_b),
            )
            continue
        existing_distance, existing_rank, _ = existing
        if candidate_distance_for_rank + 1e-6 < existing_distance:
            prepared_contact_candidates[duplicate_key] = (
                candidate_distance_for_rank,
                raw_rank,
                (raw, residue_a, residue_b),
            )
            continue
        if abs(candidate_distance_for_rank - existing_distance) <= 1e-6 and raw_rank < existing_rank:
            prepared_contact_candidates[duplicate_key] = (
                candidate_distance_for_rank,
                raw_rank,
                (raw, residue_a, residue_b),
            )

    prepared_contacts: List[Tuple[dict, dict, dict]] = [
        entry[2]
        for entry in sorted(
            prepared_contact_candidates.values(),
            key=lambda row: (row[1], row[0]),
        )
    ]

    base_pair_pair_stats = _compute_base_pair_pair_stats(
        prepared_contacts,
        residue_atoms_index=residue_atoms_index,
    )
    base_pair_support_counts = {
        pair_key: int(pair_stat.get("supportCount") or 0)
        for pair_key, pair_stat in base_pair_pair_stats.items()
    }

    for contact_rank, (raw, residue_a, residue_b) in enumerate(prepared_contacts):
        asserted = _assert_interaction(
            raw,
            residue_a,
            residue_b,
            aliases,
            residue_atoms_index=residue_atoms_index,
            base_pair_support_counts=base_pair_support_counts,
            base_pair_pair_stats=base_pair_pair_stats,
        )
        if bool(asserted.get("excludeFromNoncovalent")):
            continue
        family = str(asserted.get("family") or "").strip().lower() or "other"
        if family in {"covalent_bond", "invalid_contact"}:
            continue
        if family not in ASSERTED_FAMILY_TO_BUCKET:
            family = "other"
        terms = _normalize_arpeggio_contact_terms(raw.get("contact"))
        raw_type = str(raw.get("type") or "").strip().lower()
        interacting_entities = str(raw.get("interacting_entities") or "").strip().upper()
        raw_distance_value = _coerce_float(raw.get("distance"))

        record_residue_a = dict(residue_a)
        record_residue_b = dict(residue_b)
        atom_override_a = _primary_contact_atom_name(asserted.get("atomOverrideA"))
        atom_override_b = _primary_contact_atom_name(asserted.get("atomOverrideB"))
        element_override_a = str(asserted.get("elementOverrideA") or "").strip().upper()
        element_override_b = str(asserted.get("elementOverrideB") or "").strip().upper()
        if atom_override_a:
            record_residue_a["atom"] = atom_override_a
        if atom_override_b:
            record_residue_b["atom"] = atom_override_b
        if element_override_a:
            record_residue_a["element"] = element_override_a
        if element_override_b:
            record_residue_b["element"] = element_override_b

        distance_value, distance_resolution = _resolve_contact_distance_value(
            raw,
            record_residue_a,
            record_residue_b,
            residue_atoms_index,
        )
        distance_override_value = _coerce_float(asserted.get("distanceOverride"))
        if distance_override_value is not None:
            if distance_override_value > 0.0 or _same_atom_payload_endpoint(record_residue_a, record_residue_b):
                distance_value = distance_override_value
                distance_resolution = "distance_override_from_asserted"

        asserted_ring = asserted.get("ring") if isinstance(asserted.get("ring"), dict) else {}
        ring_payload = _build_ring_metrics_payload(
            centroid_distance=_coerce_float(asserted_ring.get("centroid_distance") or asserted_ring.get("ring_centroid_distance")),
            min_atom_distance=_coerce_float(asserted_ring.get("min_atom_distance") or asserted_ring.get("ring_min_atom_distance")),
            interplanar_distance=_coerce_float(asserted_ring.get("interplanar_distance") or asserted_ring.get("ring_interplanar_distance")),
            lateral_offset=_coerce_float(asserted_ring.get("lateral_offset") or asserted_ring.get("ring_lateral_offset")),
            normal_angle=_coerce_float(asserted_ring.get("normal_angle") or asserted_ring.get("ring_normal_angle")),
        )
        ring_payload_has_any = any(
            _coerce_float(ring_payload.get(key)) is not None
            for key in ("centroid_distance", "min_atom_distance", "interplanar_distance", "lateral_offset", "normal_angle")
        )
        ring_family = family in {"pi_pi", "aromatic_packing", "aromatic_proximal"}
        if not ring_payload_has_any and ring_family:
            computed_ring = _compute_ring_geometry_metrics(
                record_residue_a,
                record_residue_b,
                residue_atoms_index,
            )
            ring_payload = _build_ring_metrics_payload(
                centroid_distance=_coerce_float(raw.get("ring_centroid_distance")) or _coerce_float(computed_ring.get("ring_centroid_distance")),
                min_atom_distance=_compute_ring_min_atom_distance(
                    record_residue_a,
                    record_residue_b,
                    residue_atoms_index,
                ),
                interplanar_distance=_extract_ring_geometry_interplanar_distance(raw) or _coerce_float(computed_ring.get("ring_interplanar_distance")),
                lateral_offset=_extract_ring_geometry_lateral_offset(raw) or _coerce_float(computed_ring.get("ring_lateral_offset")),
                normal_angle=_extract_ring_geometry_angle(raw) or _coerce_float(computed_ring.get("ring_normal_angle")),
            )
            ring_payload_has_any = any(
                _coerce_float(ring_payload.get(key)) is not None
                for key in ("centroid_distance", "min_atom_distance", "interplanar_distance", "lateral_offset", "normal_angle")
            )
        include_ring_payload = ring_family or ring_payload_has_any

        asserted_evidence = list(asserted.get("evidence") or [])
        if family == "aromatic_packing":
            ring_min_distance = _coerce_float(ring_payload.get("min_atom_distance")) if ring_payload_has_any else None
            if (
                ring_min_distance is None
                or ring_min_distance < AROMATIC_PACKING_MIN_DISTANCE
                or ring_min_distance > AROMATIC_PACKING_MAX_DISTANCE
            ):
                family = "aromatic_proximal"
                _append_evidence(asserted_evidence, "aromatic_packing_distance_guard_reclassified_to_proximal")
                asserted["confidence"] = "low"
                if not str(asserted.get("reason_dropped") or "").strip():
                    asserted["reason_dropped"] = "aromatic_packing_distance_out_of_range"
                asserted["debugOnly"] = True

        ring_display_distance, ring_distance_kind = _resolve_ring_display_distance_for_family(
            family,
            ring_payload if include_ring_payload else None,
        )
        if ring_display_distance is not None and (
            ring_display_distance > 0.0 or _same_atom_payload_endpoint(record_residue_a, record_residue_b)
        ):
            distance_value = ring_display_distance
            distance_resolution = "distance_override_from_asserted"

        category = family
        bucket = ASSERTED_FAMILY_TO_BUCKET.get(family, "other")

        arpeggio_layer = {
            "type": raw_type,
            "terms": terms,
            "interactingEntities": interacting_entities,
        }
        if raw_distance_value is not None:
            arpeggio_layer["distance"] = round(raw_distance_value, 3)
        for geometry_key in (
            "hbond_angle",
            "hb_angle",
            "angle",
            "ring_centroid_distance",
            "ring_normal_angle",
            "ring_interplanar_distance",
            "interplanar_distance",
            "ring_plane_distance",
            "plane_distance",
            "ring_lateral_offset",
            "lateral_offset",
            "ring_offset",
            "ring_angle",
            "pi_stack_angle",
        ):
            value = _coerce_float(raw.get(geometry_key))
            if value is not None:
                arpeggio_layer[geometry_key] = round(value, 3)

        if raw.get("hbond_geometry_source"):
            arpeggio_layer["hbondGeometrySource"] = raw["hbond_geometry_source"]
            arpeggio_layer["hbondDonorSide"] = raw.get("hbond_donor_side")
        asserted_payload = {
            "family": family,
            "confidence": str(asserted.get("confidence") or "low").strip().lower() or "low",
            "evidence": asserted_evidence,
        }
        if distance_resolution == "distance_override_from_asserted":
            _append_evidence(asserted_payload["evidence"], "distance_override_from_asserted")
        elif distance_resolution == "distance_recomputed_from_coordinates":
            _append_evidence(asserted_payload["evidence"], "distance_recomputed_from_coordinates")
        elif distance_resolution == "distance_missing_or_invalid":
            _append_evidence(asserted_payload["evidence"], "distance_missing_or_invalid")
        if ring_distance_kind:
            _append_evidence(asserted_payload["evidence"], f"display_distance_{ring_distance_kind}")
        subtype = str(asserted.get("subtype") or "").strip()
        if subtype:
            asserted_payload["subtype"] = subtype
        reason_dropped = str(asserted.get("reason_dropped") or "").strip()
        if reason_dropped:
            asserted_payload["reason_dropped"] = reason_dropped
        if include_ring_payload:
            asserted_payload["ring"] = ring_payload
            asserted_payload["displayDistanceKind"] = ring_distance_kind or ""
        debug_only = bool(asserted.get("debugOnly"))
        if debug_only:
            asserted_payload["debugOnly"] = True

        atom_key_a = _build_atom_key_from_payload(record_residue_a)
        atom_key_b = _build_atom_key_from_payload(record_residue_b)
        pair_key = _build_unordered_pair_key(atom_key_a, atom_key_b)
        record = {
            "residueA": record_residue_a,
            "residueB": record_residue_b,
            "type": family,
            "category": category,
            "source": "pdbe-arpeggio",
            "atomKeyA": atom_key_a,
            "atomKeyB": atom_key_b,
            "pairKey": pair_key,
            "asserted": asserted_payload,
            "arpeggio": arpeggio_layer,
            "arpeggioType": raw_type,
            "arpeggioContact": terms,
            "interactingEntities": interacting_entities,
        }
        semantic_assertion = dict(asserted, family=family)
        record['semantics'] = asserted['semantics'] if asserted.get('semantics', {}).get('family') == family else _build_interaction_semantics(raw, semantic_assertion, record_residue_a, record_residue_b, residue_atoms_index)
        canonical_distance = record['semantics']['geometry']['distance']
        if canonical_distance.get('value') is not None:
            distance_value = canonical_distance['value']
        if debug_only:
            record["debugOnly"] = True
        if include_ring_payload:
            record["ring"] = ring_payload
            if ring_distance_kind:
                record["displayDistanceKind"] = ring_distance_kind
        if distance_value is not None:
            record["distance"] = round(distance_value, 3)
        base_pair_info = asserted.get("basePair")
        if isinstance(base_pair_info, dict):
            record["basePair"] = base_pair_info
            record["canonicalBasePair"] = bool(base_pair_info.get("isCanonicalAtomPattern"))
        metal_side = str(asserted.get("metalSide") or "").strip()
        metal_element = str(asserted.get("metalElement") or "").strip().upper()
        donor_element = str(asserted.get("donorElement") or "").strip().upper()
        if metal_side:
            record["metalSide"] = metal_side
        if metal_element:
            record["metalElement"] = metal_element
        if donor_element:
            record["donorElement"] = donor_element

        if family in AROMATIC_ASSERTED_FAMILIES:
            ring_site_keys = _resolve_aromatic_ring_site_keys(
                record_residue_a,
                record_residue_b,
                residue_atoms_index,
            )
            ring_key_a = str(ring_site_keys.get("ringKeyA") or "").strip()
            ring_key_b = str(ring_site_keys.get("ringKeyB") or "").strip()
            ring_pair_key = str(ring_site_keys.get("ringPairKey") or "").strip()
            if ring_key_a:
                record["ringKeyA"] = ring_key_a
                asserted_payload["ringKeyA"] = ring_key_a
            if ring_key_b:
                record["ringKeyB"] = ring_key_b
                asserted_payload["ringKeyB"] = ring_key_b
            if ring_pair_key:
                record["ringPairKey"] = ring_pair_key
                asserted_payload["ringPairKey"] = ring_pair_key
            for name in ("ringAtomNamesA", "ringAtomNamesB"):
                record[name] = ring_site_keys.get(name, [])
                asserted_payload[name] = record[name]

        contacts[bucket].append(record)

    # Chemical identity is independent of A/B presentation order. Separate atoms,
    # physical rings, charge groups and actual hydrogens remain separate edges.
    for bucket, records in contacts.items():
        unique = {}
        for record in records:
            identity = record['semantics']['identity']
            if identity in unique:
                unique[identity] = _merge_semantic_contact_records(unique[identity], record)
            else:
                unique[identity] = record
        contacts[bucket] = list(unique.values())
        for record in contacts[bucket]:
            family = record['semantics']['family']
            for residue in (record['residueA'], record['residueB']):
                _bump_per_residue_from_payload(per_residue, residue, family)

    coordination_groups = {}
    for record in contacts["metal_coordination"]:
        side = record.get("metalSide")
        metal = record["residueA"] if side == "A" else record["residueB"]
        donor = record["residueB"] if side == "A" else record["residueA"]
        key = (_build_atom_key_from_payload(metal), _residue_pair_identity_token(donor))
        coordination_groups.setdefault(key, []).append((record, donor.get("atom")))
    for group in coordination_groups.values():
        donor_atoms = sorted({atom for _, atom in group if atom})
        for record, _ in group:
            record["asserted"]["denticity"] = len(donor_atoms)
            record["asserted"]["coordinationDonorAtoms"] = donor_atoms
            if len(donor_atoms) > 1:
                _append_evidence(record["asserted"]["evidence"], "multidentate_coordination")


    _apply_atom_reuse_confidence_penalties(contacts)

    mode = (mode or "all").lower()
    if mode != "all":
        contacts = filter_contacts_by_mode(contacts, mode)

    interface_res_a = {
        key
        for key in per_residue
        if str(key).split(":", 1)[0] in resolved_chain_ids_a
    }
    interface_res_b = {
        key
        for key in per_residue
        if str(key).split(":", 1)[0] in resolved_chain_ids_b
    }
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

    meta = {
        "engine": "pdbe-arpeggio",
        "analysisVersion": TOOL_VERSION,
        "interactionSemanticsVersion": 1,
        "scope": "intrachain" if intrachain else "interchain",
        "classifier": "plausibility+assertion:v2",
        "note": "Contacts preserve PDBe Arpeggio plausibility and include preclassification validity/clash gating with assertion/confidence layers.",
    }
    if focus_match_keys:
        if parsed_focus:
            meta["focusResidue"] = f"{parsed_focus[0]}:{parsed_focus[1]}"
        else:
            preferred = None
            for candidate in focus_match_keys:
                if candidate[0] in resolved_chain_ids_a:
                    preferred = candidate
                    break
            if preferred is None:
                preferred = next(iter(focus_match_keys))
            meta["focusResidue"] = f"{preferred[0]}:{preferred[1]}"

    (
        chain_a,
        chain_b,
        contacts,
        per_residue,
        buried_fraction,
        meta,
    ) = _remap_report_to_external_chains(
        chain_a=chain_a,
        chain_b=chain_b,
        contacts=contacts,
        per_residue=per_residue,
        buried_fraction=buried_fraction,
        meta=meta,
        aliases=aliases,
    )

    report = {
        "chainA": chain_a,
        "chainB": chain_b,
        "analysisVersion": TOOL_VERSION,
        "contacts": contacts,
        "perResidue": per_residue,
        "interfaceArea": None,
        "buriedFraction": None,
        "contactingResidueFraction": buried_fraction,
        "approxDeltaG": None,
        "meta": meta,
    }
    from .display_grouping import attach_display_groups
    return attach_display_groups(report, residue_atoms_index, aliases)


def filter_contacts_by_mode(contacts: dict, mode: str) -> dict:
    mode = mode.lower()
    def _get(bucket: str) -> List[dict]:
        value = contacts.get(bucket)
        return value if isinstance(value, list) else []

    if mode == "hydrophobic":
        return {"hydrophobic": _get("hydrophobic")}
    if mode in {"electrostatic", "ionic", "salt"}:
        return {"salt_bridges": _get("salt_bridges")}
    if mode in {"polar", "polar_contact", "polar_contacts"}:
        return {
            "polar_contacts": _get("polar_contacts"),
            "halogen_bonds": _get("halogen_bonds"),
        }
    if mode in {"base_pair", "base_pairs", "base_pairing"}:
        # Select contextual atom contacts without rewriting their chemical family.
        return {bucket: [record for record in rows if isinstance(record.get("basePair"), dict) or record.get("type") == "base_pairing"]
                for bucket, rows in contacts.items() if isinstance(rows, list)}
    if mode in {"metal", "metal_coordination", "coordination"}:
        return {"metal_coordination": _get("metal_coordination")}
    if mode in {"hbond", "hbond_network", "hydrogen"}:
        return {"hydrogen_bonds": _get("hydrogen_bonds")}
    if mode in {"halogen", "halogen_bond", "halogen_bonds", "xbond"}:
        return {"halogen_bonds": _get("halogen_bonds")}
    if mode in {"aromatic", "pi"}:
        return {
            "pi_pi": _get("pi_pi"),
            "pi_cation": _get("pi_cation"),
            "aromatic_packing": _get("aromatic_packing"),
        }
    if mode in {"other"}:
        return {"other": _get("other")}
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


def _select_model_conformers(atoms):
    """One coherent residue conformer in the first model; blank atoms are shared."""
    if not atoms:
        return []
    first_model = atoms[0].model_id
    groups = {}
    for atom in atoms:
        if atom.model_id != first_model:
            continue
        key = (atom.chain_id, atom.chain_label, atom.res_seq)
        groups.setdefault(key, []).append(atom)
    selected = []
    for group in groups.values():
        occupancies = {}
        for atom in group:
            if atom.altloc:
                occupancies.setdefault(atom.altloc, []).append(atom.occupancy)
        preferred = min(occupancies, key=lambda alt: (
            -sum(occupancies[alt]) / len(occupancies[alt]),
            0 if alt == "A" else 1 if alt == "1" else 2, alt)) if occupancies else ""
        by_name = {}
        for atom in group:
            if atom.altloc and atom.altloc != preferred:
                continue
            name = _normalize_atom_name(atom.atom_name)
            old = by_name.get(name)
            if old is None or (not atom.altloc, atom.occupancy) > (not old.altloc, old.occupancy):
                by_name[name] = atom
        selected.extend(by_name.values())
    retained = {id(atom) for atom in selected}
    return [atom for atom in atoms if id(atom) in retained]


def _coordinate_occupancy(value):
    parsed = _coerce_float(value)
    return parsed if parsed is not None and math.isfinite(parsed) else 1.0


def parse_pdb_atoms(pdb_text: str) -> Tuple[List[AtomRecord], ChainAliases]:
    atoms: List[AtomRecord] = []
    model_id = "1"
    for line in pdb_text.splitlines():
        if line.startswith("MODEL"):
            model_id = line[10:14].strip() or "1"
            continue
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        try:
            alt_loc = line[16:17].strip()
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
            occupancy = _coordinate_occupancy(line[54:60].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except Exception:
            continue

        atoms.append(
            AtomRecord(
                chain_id=chain_id,
                chain_auth=chain_id,
                chain_label=chain_id,
                res_name=res_name,
                res_seq=res_seq,
                atom_name=atom_name,
                element=element,
                x=x,
                y=y,
                z=z,
                altloc=alt_loc,
                occupancy=occupancy,
                model_id=model_id,
            )
        )

    atoms = [atom for atom in _select_model_conformers(atoms) if atom.element not in {"H", "D"}]
    # PDB has no label/auth chain alias distinction in this flow.
    auth_ids = {str(atom.chain_id or "").strip() for atom in atoms if str(atom.chain_id or "").strip()}
    return atoms, _identity_chain_aliases(auth_ids)


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
        if line.lower() == "loop_":
            j = i + 1
            cols: List[str] = []
            while j < len(lines):
                col_line = lines[j].strip()
                if col_line.lower().startswith("_atom_site."):
                    cols.append(col_line.split()[0].lower())
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
        return [], _identity_chain_aliases(set())

    # CIF tags are case-insensitive; atom names and chain IDs are not.
    col_index = {col: idx for idx, col in enumerate(columns)}
    def idx(*names: str) -> Optional[int]:
        for name in names:
            name = name.lower()
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
    occupancy_idx = idx("_atom_site.occupancy")
    model_idx = idx("_atom_site.pdbx_PDB_model_num")

    if (
        None in (chain_idx, res_name_idx, atom_name_idx, x_idx, y_idx, z_idx)
        or (auth_res_seq_idx is None and label_res_seq_idx is None)
    ):
        return [], _identity_chain_aliases(set())

    parsed_rows: List[dict] = []
    auth_to_labels: Dict[str, Set[str]] = {}

    for row in data_rows:
        try:
            auth_chain = str(row[chain_idx] or "").strip()
            if auth_chain in {"", "?", "."}:
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
            alt_id = _clean_optional_token(row[alt_idx]) if alt_idx is not None else ""
            model_id = _clean_optional_token(row[model_idx]) if model_idx is not None else "1"
            occupancy = _coordinate_occupancy(row[occupancy_idx]) if occupancy_idx is not None else 1.0
            x = float(row[x_idx])
            y = float(row[y_idx])
            z = float(row[z_idx])
        except Exception:
            continue

        label_chain = ""
        if label_chain_idx is not None:
            label_chain = str(row[label_chain_idx] or "").strip()
        if label_chain in {"", "?", "."}:
            label_chain = auth_chain

        parsed_rows.append(
            {
                "auth_chain": auth_chain,
                "label_chain": label_chain,
                "res_name": res_name,
                "res_seq": res_seq,
                "atom_name": atom_name,
                "element": element,
                "x": x,
                "y": y,
                "z": z,
                "altloc": alt_id,
                "occupancy": occupancy,
                "model_id": model_id or "1",
            }
        )
        auth_to_labels.setdefault(auth_chain, set()).add(label_chain)

    if parsed_rows:
        first_model = parsed_rows[0]["model_id"]
        parsed_rows = [row for row in parsed_rows if row["model_id"] == first_model]
        auth_to_labels = {}
        for row in parsed_rows:
            auth_to_labels.setdefault(row["auth_chain"], set()).add(row["label_chain"])

    ambiguous_auth_ids = {
        auth_chain
        for auth_chain, label_chains in auth_to_labels.items()
        if len({str(token or "").strip() for token in label_chains if str(token or "").strip()}) > 1
    }
    auth_ids = {str(auth_chain or "").strip() for auth_chain in auth_to_labels if str(auth_chain or "").strip()}
    canonical_chain_by_pair: Dict[Tuple[str, str], str] = {}
    reserved_canonical_ids: Set[str] = set(auth_ids)

    for auth_chain in sorted(auth_to_labels):
        auth_token = str(auth_chain or "").strip()
        label_tokens = sorted(
            {
                str(token or "").strip()
                for token in auth_to_labels.get(auth_chain, set())
                if str(token or "").strip()
            }
        )
        if not auth_token:
            continue
        if auth_token not in ambiguous_auth_ids:
            canonical_chain_by_pair[(auth_token, auth_token)] = auth_token
            for label_token in label_tokens:
                canonical_chain_by_pair[(auth_token, label_token)] = auth_token
            reserved_canonical_ids.add(auth_token)
            continue
        for label_token in label_tokens:
            preferred = label_token
            if not preferred or preferred in reserved_canonical_ids:
                preferred = f"{auth_token}[{label_token or auth_token}]"
            suffix = 2
            while preferred in reserved_canonical_ids:
                preferred = f"{auth_token}[{label_token or auth_token}.{suffix}]"
                suffix += 1
            canonical_chain_by_pair[(auth_token, label_token)] = preferred
            reserved_canonical_ids.add(preferred)

    atoms: List[AtomRecord] = []
    canonical_to_auth: Dict[str, str] = {}
    auth_to_canonical_sets: Dict[str, Set[str]] = {}
    auth_residue_to_canonical_sets: Dict[Tuple[str, str, str], Set[str]] = {}
    auth_residue_atom_to_canonical_sets: Dict[Tuple[str, str, str, str], Set[str]] = {}

    def _canonical_chain_id(auth_chain: str, label_chain: str) -> str:
        auth_token = str(auth_chain or "").strip()
        label_token = str(label_chain or "").strip()
        if not auth_token and not label_token:
            return ""
        if auth_token and label_token:
            mapped = str(canonical_chain_by_pair.get((auth_token, label_token), "") or "").strip()
            if mapped:
                return mapped
        if auth_token:
            mapped = str(canonical_chain_by_pair.get((auth_token, auth_token), "") or "").strip()
            if mapped:
                return mapped
        return auth_token or label_token

    for parsed in parsed_rows:
        auth_chain = str(parsed.get("auth_chain") or "").strip()
        label_chain = str(parsed.get("label_chain") or "").strip()
        chain_id = _canonical_chain_id(auth_chain, label_chain)
        if not chain_id:
            continue
        canonical_to_auth.setdefault(chain_id, auth_chain or chain_id)
        auth_to_canonical_sets.setdefault(auth_chain or chain_id, set()).add(chain_id)

        res_name = str(parsed.get("res_name") or "").strip().upper()
        res_seq = str(parsed.get("res_seq") or "").strip()
        atom_name = str(parsed.get("atom_name") or "").strip()
        atom_key = _normalize_atom_name(atom_name)
        residue_key = (auth_chain or chain_id, res_seq, res_name)
        auth_residue_to_canonical_sets.setdefault(residue_key, set()).add(chain_id)
        if atom_key:
            atom_resolution_key = (auth_chain or chain_id, res_seq, res_name, atom_key)
            auth_residue_atom_to_canonical_sets.setdefault(atom_resolution_key, set()).add(chain_id)

        atoms.append(
            AtomRecord(
                chain_id=chain_id,
                chain_auth=auth_chain or chain_id,
                chain_label=label_chain or chain_id,
                res_name=res_name,
                res_seq=res_seq,
                atom_name=atom_name,
                element=str(parsed.get("element") or "").strip().upper(),
                x=float(parsed.get("x")),
                y=float(parsed.get("y")),
                z=float(parsed.get("z")),
                altloc=parsed["altloc"],
                occupancy=parsed["occupancy"],
                model_id=parsed["model_id"],
            )
        )

    atoms = [atom for atom in _select_model_conformers(atoms) if atom.element not in {"H", "D"}]
    canonical_ids = {str(atom.chain_id or "").strip() for atom in atoms if str(atom.chain_id or "").strip()}
    auth_to_canonical = {
        auth_chain: tuple(sorted(chain_ids))
        for auth_chain, chain_ids in auth_to_canonical_sets.items()
        if auth_chain
    }
    auth_residue_to_canonical = {
        key: tuple(sorted(chain_ids))
        for key, chain_ids in auth_residue_to_canonical_sets.items()
        if key[0] and key[1] and key[2]
    }
    auth_residue_atom_to_canonical = {
        key: tuple(sorted(chain_ids))
        for key, chain_ids in auth_residue_atom_to_canonical_sets.items()
        if key[0] and key[1] and key[2] and key[3]
    }
    query_ids = set(canonical_ids)
    query_ids.update(auth_ids)
    query_to_canonical: Dict[str, str] = {}
    sanitized_aliases: Dict[str, str] = {}
    for parsed in parsed_rows:
        label_chain = str(parsed.get("label_chain") or "").strip()
        auth_chain = str(parsed.get("auth_chain") or "").strip()
        canonical_chain = _canonical_chain_id(auth_chain, label_chain)
        if not label_chain or not auth_chain:
            continue
        if label_chain != canonical_chain and label_chain not in auth_ids:
            query_to_canonical.setdefault(label_chain, canonical_chain)
        if (
            label_chain == canonical_chain
            or label_chain == auth_chain
            or str(label_chain).lower() == str(auth_chain).lower()
            or label_chain in auth_ids
        ):
            continue
        sanitized_aliases.setdefault(label_chain, auth_chain)
    return atoms, ChainAliases(
        label_to_auth=sanitized_aliases,
        auth_ids=auth_ids,
        canonical_ids=set(canonical_ids),
        query_to_canonical=query_to_canonical,
        auth_to_canonical=auth_to_canonical,
        canonical_to_auth=canonical_to_auth,
        auth_residue_to_canonical=auth_residue_to_canonical,
        auth_residue_atom_to_canonical=auth_residue_atom_to_canonical,
        query_ids=query_ids,
    )


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
        if stripped.startswith("_") or stripped.lower().startswith(("loop_", "data_", "save_", "stop_")):
            break
        if stripped.startswith(";"):
            block, i = read_semicolon_block(lines, i)
            tokens.append(block)
        else:
            tokens.extend(_split_mmcif_row(raw))
            i += 1
        while len(tokens) >= columns:
            row = tokens[:columns]
            tokens = tokens[columns:]
            yield row
    return


def _split_mmcif_row(raw: str) -> List[str]:
    """Split a single mmCIF loop row into data tokens.

    mmCIF uses whitespace-separated tokens, optionally quoted with single or
    double quotes. Unquoted tokens may legitimately contain apostrophes (e.g.
    atom names like O5'), which breaks shlex-based tokenization.
    """
    text = str(raw or "")
    n = len(text)
    i = 0
    out: List[str] = []
    while i < n:
        while i < n and text[i].isspace():
            i += 1
        if i >= n:
            break
        quote = text[i]
        if quote in {"'", '"'}:
            i += 1
            start = i
            while i < n:
                ch = text[i]
                if ch == quote and (i + 1 == n or text[i + 1].isspace()):
                    break
                i += 1
            out.append(text[start:i])
            if i < n and text[i] == quote:
                i += 1
            continue
        start = i
        while i < n and not text[i].isspace():
            i += 1
        out.append(text[start:i])
    return out


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


def cache_key(
    pdb_id: Optional[str],
    mmcif_text: Optional[str],
    chain_a: str,
    chain_b: str,
    mode: str,
    focus_residue: Optional[str] = None,
) -> str:
    if mmcif_text:
        source = hashlib.sha256(mmcif_text.encode("utf-8")).hexdigest()
    else:
        source = str(pdb_id or "").lower()
    focus_token = str(focus_residue or "").strip()
    return f"{source}:{chain_a}:{chain_b}:{mode}:{focus_token}:{TOOL_VERSION}"
