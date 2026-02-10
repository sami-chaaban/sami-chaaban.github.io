from __future__ import annotations

import json
import tempfile
import warnings
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.interpolate import CubicSpline

from Bio.PDB import PDBIO, PDBParser, MMCIFParser
from Bio.PDB.DSSP import DSSP


Vec3 = np.ndarray  # shape (3,)


def _norm(v: Vec3, eps: float = 1e-12) -> float:
    return float(np.linalg.norm(v) + eps)


def _normalize(v: Vec3, eps: float = 1e-12) -> Vec3:
    n = np.linalg.norm(v)
    if n < eps:
        return np.array([1.0, 0.0, 0.0], dtype=np.float32)
    return (v / n).astype(np.float32)


def _project_out(v: Vec3, axis: Vec3) -> Vec3:
    return v - axis * float(np.dot(v, axis))


def _cumulative_lengths(points: np.ndarray) -> np.ndarray:
    d = np.linalg.norm(points[1:] - points[:-1], axis=1)
    s = np.concatenate([[0.0], np.cumsum(d)])
    return s


@dataclass
class ResidueBackbone:
    res_id: Tuple[str, int, str]  # (chain_id, resseq, icode)
    ca: Vec3
    c: Optional[Vec3]
    o: Optional[Vec3]
    ss: str  # 'H', 'E', 'C'


def load_structure(path: str):
    if path.lower().endswith(".cif") or path.lower().endswith(".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", path)
    return structure


def assign_ss_dssp(model, path: str) -> Dict[Tuple[str, int, str], str]:
    ss_map: Dict[Tuple[str, int, str], str] = {}
    temp_path = None
    dssp_path = path
    file_type = None
    if path.lower().endswith((".cif", ".mmcif")):
        temp = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
        temp.close()
        io = PDBIO()
        io.set_structure(model)
        io.save(temp.name)
        temp_path = temp.name
        dssp_path = temp.name
        file_type = "PDB"

    try:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="parse error at line 1: This file does not seem to be an mmCIF file",
            )
            if file_type:
                try:
                    dssp = DSSP(model, dssp_path, "mkdssp", "Sander", file_type)
                except TypeError:
                    try:
                        dssp = DSSP(model, dssp_path, dssp="mkdssp", file_type=file_type)
                    except TypeError:
                        dssp = DSSP(model, dssp_path, dssp="mkdssp")
            else:
                dssp = DSSP(model, dssp_path, dssp="mkdssp")
        for (chain_id, res_id), d in dssp.property_dict.items():
            resseq = int(res_id[1])
            icode = str(res_id[2]).strip() if res_id[2] != " " else ""
            if isinstance(d, dict):
                raw = d.get("ss")
            elif isinstance(d, (tuple, list)) and len(d) > 2:
                raw = d[2]
            else:
                raw = None
            if raw is None:
                continue
            if raw in ("H", "G", "I"):
                ss = "H"
            elif raw in ("E", "B"):
                ss = "E"
            else:
                ss = "C"
            ss_map[(chain_id, resseq, icode)] = ss
    except Exception:
        pass
    finally:
        if temp_path:
            try:
                Path(temp_path).unlink(missing_ok=True)
            except Exception:
                pass
    return ss_map


def extract_backbone(model, ss_map: Dict[Tuple[str, int, str], str]) -> Dict[str, List[ResidueBackbone]]:
    chains: Dict[str, List[ResidueBackbone]] = {}
    for chain in model:
        cid = chain.id
        residues: List[ResidueBackbone] = []
        for res in chain:
            hetflag = res.id[0]
            if hetflag.strip():
                continue
            resseq = int(res.id[1])
            icode = str(res.id[2]).strip() if res.id[2] != " " else ""
            key = (cid, resseq, icode)

            if "CA" not in res:
                continue
            ca = np.array(res["CA"].get_coord(), dtype=np.float32)

            c = np.array(res["C"].get_coord(), dtype=np.float32) if "C" in res else None
            o = np.array(res["O"].get_coord(), dtype=np.float32) if "O" in res else None

            ss = ss_map.get(key, "C")
            residues.append(ResidueBackbone(key, ca, c, o, ss))
        if residues:
            chains[cid] = residues
    return chains


def split_continuous_segments(
    residues: List[ResidueBackbone], max_ca_gap: float = 4.5
) -> List[List[ResidueBackbone]]:
    segs: List[List[ResidueBackbone]] = []
    cur: List[ResidueBackbone] = []
    for r in residues:
        if not cur:
            cur.append(r)
            continue
        gap = _norm(r.ca - cur[-1].ca)
        if gap > max_ca_gap:
            if len(cur) >= 4:
                segs.append(cur)
            cur = [r]
        else:
            cur.append(r)
    if len(cur) >= 4:
        segs.append(cur)
    return segs


def catmull_like_spline(points: np.ndarray) -> Tuple[CubicSpline, CubicSpline, CubicSpline, np.ndarray]:
    s = _cumulative_lengths(points)
    if np.any(np.diff(s) < 1e-6):
        s = s + np.linspace(0, 1e-3, len(s))

    csx = CubicSpline(s, points[:, 0], bc_type="natural")
    csy = CubicSpline(s, points[:, 1], bc_type="natural")
    csz = CubicSpline(s, points[:, 2], bc_type="natural")
    return csx, csy, csz, s


def sample_spline(
    csx, csy, csz, s_max: float, step: float
) -> Tuple[np.ndarray, np.ndarray]:
    s_samples = np.arange(0.0, s_max, step, dtype=np.float32)
    x = csx(s_samples)
    y = csy(s_samples)
    z = csz(s_samples)
    pos = np.stack([x, y, z], axis=1).astype(np.float32)

    dx = csx(s_samples, 1)
    dy = csy(s_samples, 1)
    dz = csz(s_samples, 1)
    tan = np.stack([dx, dy, dz], axis=1).astype(np.float32)
    tan = np.array([_normalize(t) for t in tan], dtype=np.float32)
    return pos, tan


def nearest_residue_index_by_arclen(res_s: np.ndarray, sample_s: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(res_s, sample_s, side="left")
    idx = np.clip(idx, 0, len(res_s) - 1)
    idx2 = np.clip(idx - 1, 0, len(res_s) - 1)
    choose_left = np.abs(res_s[idx2] - sample_s) < np.abs(res_s[idx] - sample_s)
    idx = np.where(choose_left, idx2, idx)
    return idx.astype(np.int32)


def rotation_minimizing_frames(tan: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    K = tan.shape[0]
    n = np.zeros((K, 3), dtype=np.float32)
    b = np.zeros((K, 3), dtype=np.float32)

    t0 = tan[0]
    ref = np.array([0.0, 0.0, 1.0], dtype=np.float32)
    if abs(float(np.dot(t0, ref))) > 0.9:
        ref = np.array([0.0, 1.0, 0.0], dtype=np.float32)
    n0 = _normalize(np.cross(ref, t0))
    b0 = _normalize(np.cross(t0, n0))
    n[0], b[0] = n0, b0

    for i in range(1, K):
        t_prev = tan[i - 1]
        t_cur = tan[i]
        v = t_prev + t_cur
        if _norm(v) < 1e-6:
            n[i] = n[i - 1]
            b[i] = b[i - 1]
            continue
        v = _normalize(v)

        ni = n[i - 1] - 2.0 * float(np.dot(n[i - 1], v)) * v
        ni = _normalize(_project_out(ni, t_cur))
        bi = _normalize(np.cross(t_cur, ni))
        n[i], b[i] = ni, bi

    return n, b


def oxygen_guided_normals(
    tan: np.ndarray,
    residues: List[ResidueBackbone],
    res_idx_for_sample: np.ndarray,
    smooth_window: int = 7,
) -> Optional[np.ndarray]:
    K = tan.shape[0]
    n = np.zeros((K, 3), dtype=np.float32)

    ok = 0
    for i in range(K):
        r = residues[int(res_idx_for_sample[i])]
        if r.o is None:
            continue
        origin = r.c if r.c is not None else r.ca
        ovec = (r.o - origin).astype(np.float32)
        ni = _project_out(ovec, tan[i])
        if _norm(ni) < 1e-6:
            continue
        n[i] = _normalize(ni)
        ok += 1

    if ok < max(10, int(0.3 * K)):
        return None

    w = smooth_window
    pad = w // 2
    n_pad = np.pad(n, ((pad, pad), (0, 0)), mode="edge")
    n_s = np.zeros_like(n)
    for i in range(K):
        avg = np.mean(n_pad[i : i + w], axis=0)
        avg = _project_out(avg, tan[i])
        n_s[i] = _normalize(avg)

    return n_s


def build_samples_for_segment(
    residues: List[ResidueBackbone],
    step: float = 0.35,
    arrow_len_angstrom: float = 6.0,
) -> Dict:
    ca = np.stack([r.ca for r in residues], axis=0)
    res_s = _cumulative_lengths(ca)

    csx, csy, csz, s = catmull_like_spline(ca)
    s_max = float(s[-1])
    pos, tan = sample_spline(csx, csy, csz, s_max, step)
    K = pos.shape[0]
    sample_s = np.arange(0.0, K * step, step, dtype=np.float32)
    sample_s = np.clip(sample_s, 0.0, s_max)

    res_idx = nearest_residue_index_by_arclen(res_s, sample_s)

    n_oxy = oxygen_guided_normals(tan, residues, res_idx)
    if n_oxy is None:
        n, b = rotation_minimizing_frames(tan)
    else:
        n = n_oxy
        b = np.array([_normalize(np.cross(tan[i], n[i])) for i in range(K)], dtype=np.float32)

    ss = np.array([residues[int(j)].ss for j in res_idx], dtype="U1")

    w = np.zeros((K,), dtype=np.float32)
    h = np.zeros((K,), dtype=np.float32)
    for i in range(K):
        if ss[i] == "H":
            w[i], h[i] = 2.2, 0.60
        elif ss[i] == "E":
            w[i], h[i] = 2.4, 0.25
        else:
            w[i], h[i] = 0.85, 0.85

    isE = ss == "E"
    i = 0
    while i < K:
        if not isE[i]:
            i += 1
            continue
        j = i
        while j < K and isE[j]:
            j += 1
        seg_len = float((j - i) * step)
        taper_len = min(arrow_len_angstrom, seg_len)
        taper_samples = max(1, int(taper_len / step))
        start = max(i, j - taper_samples)
        for t_i in range(start, j):
            u = (t_i - start) / max(1, (j - start - 1))
            w[t_i] *= 1.0 + 0.8 * u
            h[t_i] *= 1.0 - 0.95 * u
        i = j

    def flat(a: np.ndarray) -> List[float]:
        return a.reshape(-1).astype(np.float32).tolist()

    return {
        "pos": flat(pos),
        "n": flat(n),
        "b": flat(b),
        "w": w.tolist(),
        "h": h.tolist(),
        "ss": ss.tolist(),
        "resIndex": res_idx.tolist(),
    }


def structure_to_ribbon_json(path: str, model_index: int = 0, step: float = 0.35) -> Dict:
    structure = load_structure(path)
    model = list(structure.get_models())[model_index]

    ss_map = assign_ss_dssp(model, path)
    chains = extract_backbone(model, ss_map)

    out = {"source": path, "step": step, "chains": []}

    for cid, residues in chains.items():
        segments = split_continuous_segments(residues)
        chain_out = {"id": cid, "segments": []}
        for seg in segments:
            samples = build_samples_for_segment(seg, step=step)
            chain_out["segments"].append({"samples": samples})
        out["chains"].append(chain_out)

    return out


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="PDB or mmCIF file")
    ap.add_argument("--step", type=float, default=0.35, help="sampling step in Å")
    ap.add_argument("--out", default="ribbon.json")
    args = ap.parse_args()

    data = structure_to_ribbon_json(args.input, step=args.step)
    with open(args.out, "w") as f:
        json.dump(data, f)
    print(f"Wrote {args.out}")
