from __future__ import annotations

import json
import os
import sys
import tempfile
from typing import Any, Dict, List, Optional, Tuple


def _suppress_stdio():
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    stdout_fd = os.dup(1)
    stderr_fd = os.dup(2)
    os.dup2(devnull_fd, 1)
    os.dup2(devnull_fd, 2)
    return stdout_fd, stderr_fd, devnull_fd


def _restore_stdio(stdout_fd: int, stderr_fd: int, devnull_fd: int) -> None:
    os.dup2(stdout_fd, 1)
    os.dup2(stderr_fd, 2)
    os.close(stdout_fd)
    os.close(stderr_fd)
    os.close(devnull_fd)


def _mesh_to_json(mesh) -> Dict[str, Any]:
    vertices = getattr(mesh, "vertices", None) or []
    triangles = getattr(mesh, "triangles", None) or []
    vertex_count = len(vertices)
    triangle_count = len(triangles)
    positions: List[float] = [float(v.pos[i]) for v in vertices for i in (0, 1, 2)]
    normals: List[float] = [float(v.normal[i]) for v in vertices for i in (0, 1, 2)]
    colors: List[float] = [float(v.color[i]) for v in vertices for i in (0, 1, 2, 3)]
    indices: List[int] = [int(tri.point_id[i]) for tri in triangles for i in (0, 1, 2)]

    return {
        "positions": positions,
        "normals": normals,
        "colors": colors,
        "indices": indices,
        "vertexCount": vertex_count,
        "triangleCount": triangle_count,
        "name": getattr(mesh, "name", ""),
        "status": getattr(mesh, "status", None),
    }


def _empty_mesh(status: str = "empty") -> Dict[str, Any]:
    return {
        "positions": [],
        "normals": [],
        "colors": [],
        "indices": [],
        "vertexCount": 0,
        "triangleCount": 0,
        "name": "",
        "status": status,
    }


def _gltf_component_size(component_type: int) -> int:
    if component_type == 5126:  # FLOAT
        return 4
    if component_type == 5125:  # UNSIGNED_INT
        return 4
    if component_type == 5123:  # UNSIGNED_SHORT
        return 2
    if component_type == 5121:  # UNSIGNED_BYTE
        return 1
    raise ValueError(f"Unsupported glTF component type {component_type}")


def _gltf_num_components(accessor_type: str) -> int:
    if accessor_type == "SCALAR":
        return 1
    if accessor_type == "VEC2":
        return 2
    if accessor_type == "VEC3":
        return 3
    if accessor_type == "VEC4":
        return 4
    raise ValueError(f"Unsupported glTF accessor type {accessor_type}")


def _read_accessor(
    accessor: dict,
    buffer_view: dict,
    bin_blob: bytes,
) -> List[float]:
    component_type = accessor["componentType"]
    count = accessor["count"]
    accessor_type = accessor["type"]
    num_components = _gltf_num_components(accessor_type)
    component_size = _gltf_component_size(component_type)
    normalized = bool(accessor.get("normalized"))

    view_offset = buffer_view.get("byteOffset", 0)
    accessor_offset = accessor.get("byteOffset", 0)
    stride = buffer_view.get("byteStride") or num_components * component_size
    start = view_offset + accessor_offset
    total = count * num_components
    out: List[float] = [0.0] * total

    import struct

    fmt_char = {
        5126: "f",  # FLOAT
        5125: "I",  # UNSIGNED_INT
        5123: "H",  # UNSIGNED_SHORT
        5121: "B",  # UNSIGNED_BYTE
    }.get(component_type)

    packed_stride = num_components * component_size
    if fmt_char and stride == packed_stride and count > 0:
        mv = memoryview(bin_blob)[start : start + count * stride]
        pack_fmt = "<" + (fmt_char * num_components)
        out_i = 0
        if normalized:
            scale = (
                1.0 / 255.0
                if component_type == 5121
                else 1.0 / 65535.0
                if component_type == 5123
                else 1.0 / 4294967295.0
                if component_type == 5125
                else 1.0
            )
            for vals in struct.iter_unpack(pack_fmt, mv):
                for value in vals:
                    out[out_i] = float(value) * scale
                    out_i += 1
        else:
            for vals in struct.iter_unpack(pack_fmt, mv):
                for value in vals:
                    out[out_i] = float(value)
                    out_i += 1
        return out

    out_i = 0
    for i in range(count):
        base = start + i * stride
        for c in range(num_components):
            offset = base + c * component_size
            if component_type == 5126:
                value = struct.unpack_from("<f", bin_blob, offset)[0]
            elif component_type == 5125:
                value = struct.unpack_from("<I", bin_blob, offset)[0]
            elif component_type == 5123:
                value = struct.unpack_from("<H", bin_blob, offset)[0]
            elif component_type == 5121:
                value = struct.unpack_from("<B", bin_blob, offset)[0]
            else:
                raise ValueError(f"Unsupported glTF component type {component_type}")
            if normalized:
                if component_type == 5121:
                    value = value / 255.0
                elif component_type == 5123:
                    value = value / 65535.0
                elif component_type == 5125:
                    value = value / 4294967295.0
            out[out_i] = float(value)
            out_i += 1
    if len(out) != total:
        raise ValueError("Failed to read glTF accessor")
    return out


def _parse_gltf_glb(path: str) -> Dict[str, Any]:
    with open(path, "rb") as fh:
        data = fh.read()
    if data[:4] != b"glTF":
        raise ValueError("Not a glTF binary file")
    length = int.from_bytes(data[8:12], "little")
    if length != len(data):
        data = data[:length]
    offset = 12
    json_chunk: Optional[bytes] = None
    bin_chunk: Optional[bytes] = None
    while offset + 8 <= len(data):
        chunk_len = int.from_bytes(data[offset : offset + 4], "little")
        chunk_type = int.from_bytes(data[offset + 4 : offset + 8], "little")
        offset += 8
        chunk = data[offset : offset + chunk_len]
        offset += chunk_len
        if chunk_type == 0x4E4F534A:  # JSON
            json_chunk = chunk
        elif chunk_type == 0x004E4942:  # BIN
            bin_chunk = chunk
    if not json_chunk or bin_chunk is None:
        raise ValueError("Missing glTF JSON or BIN chunk")
    gltf = json.loads(json_chunk.decode("utf-8"))
    meshes = gltf.get("meshes") or []
    if not meshes:
        raise ValueError("No meshes in glTF")
    buffer_views = gltf.get("bufferViews") or []
    accessors = gltf.get("accessors") or []
    primitives = meshes[0].get("primitives") or []
    if not primitives:
        raise ValueError("No primitives in glTF mesh")

    positions: List[float] = []
    normals: List[float] = []
    colors: List[float] = []
    indices: List[int] = []
    vertex_offset = 0

    for prim in primitives:
        attrs = prim.get("attributes") or {}
        pos_idx = attrs.get("POSITION")
        if pos_idx is None:
            continue
        pos_acc = accessors[pos_idx]
        pos_view = buffer_views[pos_acc["bufferView"]]
        pos = _read_accessor(pos_acc, pos_view, bin_chunk)
        positions.extend(pos)

        nrm_idx = attrs.get("NORMAL")
        if nrm_idx is not None:
            nrm_acc = accessors[nrm_idx]
            nrm_view = buffer_views[nrm_acc["bufferView"]]
            normals.extend(_read_accessor(nrm_acc, nrm_view, bin_chunk))

        col_idx = attrs.get("COLOR_0")
        if col_idx is not None:
            col_acc = accessors[col_idx]
            col_view = buffer_views[col_acc["bufferView"]]
            col = _read_accessor(col_acc, col_view, bin_chunk)
            if col_acc.get("type") == "VEC3":
                # pad alpha
                for i in range(0, len(col), 3):
                    colors.extend([col[i], col[i + 1], col[i + 2], 1.0])
            else:
                colors.extend(col)

        idx_idx = prim.get("indices")
        if idx_idx is not None:
            idx_acc = accessors[idx_idx]
            idx_view = buffer_views[idx_acc["bufferView"]]
            idx_vals = _read_accessor(idx_acc, idx_view, bin_chunk)
            indices.extend([int(i) + vertex_offset for i in idx_vals])
        else:
            count = int(len(pos) / 3)
            indices.extend([i + vertex_offset for i in range(count)])

        vertex_offset += int(len(pos) / 3)

    return {
        "positions": positions,
        "normals": normals,
        "colors": colors,
        "indices": indices,
        "vertexCount": int(len(positions) / 3),
        "triangleCount": int(len(indices) / 3),
        "name": meshes[0].get("name", ""),
        "status": None,
    }


def _write_temp_structure(text: str, suffix: str) -> str:
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    try:
        tmp.write(text.encode("utf-8"))
        tmp.flush()
    finally:
        tmp.close()
    return tmp.name


def _sanitize_pdb_for_chapi(text: str) -> str:
    if not text:
        return text
    lines = text.splitlines()
    kept = [line for line in lines if not line.startswith("SEQRES")]
    if len(kept) == len(lines):
        return text
    return "\n".join(kept) + "\n"


def _run_payload(payload: dict) -> Dict[str, Any]:
    if not isinstance(payload, dict):
        raise RuntimeError("Payload must be a JSON object.")
    text = payload.get("text")
    fmt = payload.get("format")
    if not text or fmt not in ("pdb", "mmcif"):
        raise RuntimeError("Payload must include text and format ('pdb' or 'mmcif').")

    rep = payload.get("representation", "bonds")
    split_by_chain = bool(payload.get("splitByChain", False))

    temp_path = None
    suppressed = None
    error = None
    result = None

    try:
        suppressed = _suppress_stdio()
        import coot_headless_api as ch  # type: ignore

        suffix = ".pdb" if fmt == "pdb" else ".cif"
        if fmt == "pdb":
            text = _sanitize_pdb_for_chapi(text)
        temp_path = _write_temp_structure(text, suffix)

        container = ch.molecules_container_t(False)
        try:
            container.geometry_init_standard()
        except Exception:
            pass
        if fmt == "pdb":
            imol = container.read_pdb(temp_path)
        else:
            imol = container.read_coordinates(temp_path)

        if imol is None or int(imol) < 0:
            raise RuntimeError("Failed to read structure in chapi bridge.")

        if rep == "bonds":
            mode = payload.get("mode", "COLOUR-BY-CHAIN-AND-DICTIONARY")
            against_dark = bool(payload.get("againstDarkBackground", False))
            bond_width = float(payload.get("bondWidth", 0.12))
            atom_ratio = float(payload.get("atomRadiusToBondWidthRatio", 1.0))
            smoothness = int(payload.get("smoothnessFactor", 2))
            mesh = container.get_bonds_mesh(
                int(imol),
                mode,
                against_dark,
                bond_width,
                atom_ratio,
                smoothness,
            )
        elif rep == "bonds-selection":
            selection_cid = payload.get("cid", "//")
            mode = payload.get("mode", "COLOUR-BY-CHAIN-AND-DICTIONARY")
            against_dark = bool(payload.get("againstDarkBackground", False))
            bond_width = float(payload.get("bondWidth", 0.12))
            atom_ratio = float(payload.get("atomRadiusToBondWidthRatio", 1.0))
            smoothness = int(payload.get("smoothnessFactor", 2))
            draw_h = bool(payload.get("drawHydrogenAtoms", False))
            draw_missing = bool(payload.get("drawMissingResidueLoops", False))
            non_draw_cids = payload.get("nonDrawCids") or []
            carbon_color = payload.get("carbonColor")
            glb_path = tempfile.NamedTemporaryFile(delete=False, suffix=".glb").name
            try:
                if non_draw_cids:
                    for cid in non_draw_cids:
                        if cid:
                            try:
                                container.add_to_non_drawn_bonds(int(imol), str(cid))
                            except Exception:
                                pass
                if carbon_color:
                    try:
                        container.set_use_bespoke_carbon_atom_colour(int(imol), True)
                        container.set_bespoke_carbon_atom_colour(int(imol), str(carbon_color))
                    except Exception:
                        pass
                container.export_model_molecule_as_gltf(
                    int(imol),
                    selection_cid,
                    mode,
                    against_dark,
                    bond_width,
                    atom_ratio,
                    smoothness,
                    draw_h,
                    draw_missing,
                    glb_path,
                )
                try:
                    if not os.path.exists(glb_path) or os.path.getsize(glb_path) < 20:
                        result = _empty_mesh("empty")
                    else:
                        result = _parse_gltf_glb(glb_path)
                except Exception:
                    result = _empty_mesh("empty")
            finally:
                try:
                    if non_draw_cids:
                        container.clear_non_drawn_bonds(int(imol))
                except Exception:
                    pass
                try:
                    if carbon_color:
                        container.set_use_bespoke_carbon_atom_colour(int(imol), False)
                except Exception:
                    pass
                try:
                    os.unlink(glb_path)
                except Exception:
                    pass
        elif rep in ("ribbon", "surface"):
            colour_scheme = payload.get("colourScheme", "Chain")
            style = payload.get("style", "Ribbon" if rep == "ribbon" else "MolecularSurface")
            ss_flag = int(payload.get("secondaryStructureUsage", 0))
            if split_by_chain:
                requested_chain_ids_raw = payload.get("chainIds")
                requested_chain_ids: List[str] = []
                if isinstance(requested_chain_ids_raw, list):
                    seen_requested: set[str] = set()
                    for chain in requested_chain_ids_raw:
                        chain_id = str(chain).strip()
                        if not chain_id or chain_id in seen_requested:
                            continue
                        seen_requested.add(chain_id)
                        requested_chain_ids.append(chain_id)
                requested_chain: Optional[str] = None
                requested_cid = str(payload.get("cid", "//") or "//").strip()
                if requested_cid.startswith("//") and requested_cid not in {"//", "///"}:
                    chain_token = requested_cid[2:].split("/", 1)[0].strip()
                    if chain_token and not any(sep in chain_token for sep in (",", ";", "|", " ")):
                        requested_chain = chain_token

                raw_chain_ids = (
                    requested_chain_ids
                    if requested_chain_ids
                    else (container.get_chains_in_model(int(imol)) or [])
                )
                chain_ids: List[str] = []
                seen_chain_ids: set[str] = set()
                for chain in raw_chain_ids:
                    chain_id = str(chain).strip()
                    if not chain_id or chain_id in seen_chain_ids:
                        continue
                    if requested_chain and chain_id != requested_chain:
                        continue
                    seen_chain_ids.add(chain_id)
                    chain_ids.append(chain_id)
                meshes = []
                for chain in chain_ids:
                    cid = f"//{chain}"
                    mesh = container.get_molecular_representation_mesh(
                        int(imol),
                        cid,
                        colour_scheme,
                        style,
                        ss_flag,
                    )
                    if getattr(mesh, "vertices", None):
                        mesh_json = _mesh_to_json(mesh)
                        if mesh_json.get("vertexCount", 0) > 0:
                            mesh_json["chainId"] = chain
                            meshes.append(mesh_json)
                if meshes:
                    result = {
                        "meshType": "chains",
                        "meshes": meshes,
                        "chainIds": [m["chainId"] for m in meshes],
                    }
                else:
                    cid = payload.get("cid", "//")
                    mesh = container.get_molecular_representation_mesh(
                        int(imol),
                        cid,
                        colour_scheme,
                        style,
                        ss_flag,
                    )
                    result = _mesh_to_json(mesh)
            else:
                cid = payload.get("cid", "//")
                mesh = container.get_molecular_representation_mesh(
                    int(imol),
                    cid,
                    colour_scheme,
                    style,
                    ss_flag,
                )
        else:
            raise RuntimeError(f"Unknown representation '{rep}'.")
        if result is None:
            result = _mesh_to_json(mesh)
    except Exception as exc:
        error = str(exc)
    finally:
        if temp_path:
            try:
                os.unlink(temp_path)
            except Exception:
                pass
        if suppressed:
            _restore_stdio(*suppressed)

    if error:
        raise RuntimeError(error)
    return result or _empty_mesh("empty")


def _write_server_response(ok: bool, body: str) -> None:
    prefix = "OK" if ok else "ERR"
    sys.stdout.write(f"{prefix}\t{body}\n")
    sys.stdout.flush()


def _run_server() -> int:
    for raw_line in sys.stdin:
        line = raw_line.strip()
        if not line:
            continue
        try:
            request_obj = json.loads(line)
        except Exception as exc:
            error_body = json.dumps({"error": f"Invalid JSON input: {exc}"}, separators=(",", ":"), ensure_ascii=False)
            _write_server_response(False, error_body)
            continue

        if isinstance(request_obj, dict) and request_obj.get("op") == "ping":
            _write_server_response(True, '{"pong":true}')
            continue

        payload = request_obj.get("payload") if isinstance(request_obj, dict) and "payload" in request_obj else request_obj
        try:
            result = _run_payload(payload)
            result_json = json.dumps(result, separators=(",", ":"), ensure_ascii=False)
            _write_server_response(True, result_json)
        except Exception as exc:
            error_body = json.dumps({"error": str(exc)}, separators=(",", ":"), ensure_ascii=False)
            _write_server_response(False, error_body)
    return 0


def main() -> int:
    if any(arg == "--server" for arg in sys.argv[1:]):
        return _run_server()

    raw = sys.stdin.read()
    if not raw.strip():
        sys.stderr.write("No input provided to chapi bridge.\n")
        return 2

    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as exc:
        sys.stderr.write(f"Invalid JSON input: {exc}\n")
        return 2

    try:
        result = _run_payload(payload)
    except Exception as exc:
        sys.stderr.write(str(exc) + "\n")
        return 1

    sys.stdout.write(json.dumps(result, separators=(",", ":"), ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
