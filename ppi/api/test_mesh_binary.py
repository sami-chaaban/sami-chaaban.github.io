"""Packed mesh wire-format and bridge regressions, without native Coot."""
from __future__ import annotations

import io
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from api import chapi_bridge
from api.mesh_binary import BinaryMesh, MAGIC


def native_mesh():
    return SimpleNamespace(
        vertices=[SimpleNamespace(pos=(1.25, -2.5, 3), normal=(0, 0, 1), color=(0.2, 0.4, 0.6, 1))
                  for _ in range(3)],
        triangles=[SimpleNamespace(point_id=(2, 1, 0))], name="chain α", status=1,
    )


def decode(body):
    assert body[:8] == MAGIC
    metadata_length = struct.unpack_from("<I", body, 8)[0]
    metadata = json.loads(body[12:12 + metadata_length])
    buffer_start = (12 + metadata_length + 3) & ~3

    def unpack(descriptor):
        code = "f" if descriptor["type"] == "float32" else "I"
        return struct.unpack_from("<" + code * descriptor["count"], body, buffer_start + descriptor["offset"])

    return metadata, unpack, buffer_start


def serialize(binary):
    stream = io.BytesIO()
    binary.write_to(stream)
    return stream.getvalue()


class BinaryEncodingTests(unittest.TestCase):
    def test_multiple_meshes_are_aligned_and_preserve_geometry_and_identity(self):
        binary = BinaryMesh()
        meshes = [binary.native_mesh(native_mesh()), binary.native_mesh(native_mesh())]
        for chain, mesh in zip(("A", "B"), meshes):
            mesh["chainId"] = chain
        binary.metadata = {"meshType": "chains", "meshes": meshes, "chainIds": ["A", "B"]}
        body = serialize(binary)
        metadata, unpack, data_start = decode(body)
        self.assertEqual(data_start % 4, 0)
        next_offset = 0
        for i, mesh in enumerate(metadata["meshes"]):
            self.assertEqual(mesh["name"], "chain α")
            self.assertEqual(mesh["chainId"], "AB"[i])
            self.assertEqual(mesh["status"], 1)
            self.assertEqual(mesh["vertexCount"], 3)
            self.assertEqual(mesh["triangleCount"], 1)
            self.assertEqual(unpack(mesh["positions"]), (1.25, -2.5, 3) * 3)
            self.assertEqual(unpack(mesh["indices"]), (2, 1, 0))
            for field in ("positions", "normals", "colors", "indices"):
                descriptor = mesh[field]
                self.assertEqual(descriptor["offset"], next_offset)
                next_offset += descriptor["count"] * 4
        self.assertEqual(data_start + next_offset, len(body))

    def test_nonfinite_and_unrepresentable_components_have_finite_defaults(self):
        mesh = native_mesh()
        mesh.vertices[0].pos = (float("nan"), float("inf"), 1e100)
        mesh.vertices[0].normal = (-float("inf"), None, "invalid")
        mesh.vertices[0].color = (float("nan"), float("inf"), 1e100, float("nan"))
        binary = BinaryMesh()
        binary.metadata = binary.native_mesh(mesh)
        metadata, unpack, _ = decode(serialize(binary))
        self.assertEqual(unpack(metadata["positions"])[:3], (0, 0, 0))
        self.assertEqual(unpack(metadata["normals"])[:3], (0, 0, 0))
        self.assertEqual(unpack(metadata["colors"])[:4], (0, 0, 0, 1))

    def test_invalid_triangle_indices_and_inconsistent_buffers_are_rejected(self):
        for index in (-1, 3, 2**32, 0.5):
            mesh = native_mesh()
            mesh.triangles[0].point_id = (0, 1, index)
            with self.subTest(index=index), self.assertRaises(ValueError):
                BinaryMesh().native_mesh(mesh)
        mesh = native_mesh()
        mesh.vertices[0].normal = (0, 1)
        with self.assertRaisesRegex(ValueError, "normal count"):
            BinaryMesh().native_mesh(mesh)

    def test_native_vector_properties_are_fetched_once(self):
        class Vertex:
            def __init__(self):
                self.reads = {"pos": 0, "normal": 0, "color": 0}

            def __getattr__(self, name):
                self.reads[name] += 1
                return (0, 0, 0, 1) if name == "color" else (0, 0, 1)

        mesh = native_mesh()
        mesh.vertices = [Vertex() for _ in range(3)]
        BinaryMesh().native_mesh(mesh)
        self.assertTrue(all(v.reads == {"pos": 1, "normal": 1, "color": 1} for v in mesh.vertices))

    def test_empty_chain_result_still_has_valid_envelope(self):
        binary = BinaryMesh()
        binary.metadata = {"meshType": "chains", "meshes": [], "chainIds": [], "status": "empty"}
        body = serialize(binary)
        metadata, _, data_start = decode(body)
        self.assertEqual(metadata, binary.metadata)
        self.assertEqual(data_start, len(body))

    def test_gltf_mesh_allows_missing_optional_normals_and_colors(self):
        mesh = chapi_bridge._mesh_to_json(native_mesh())
        mesh["normals"] = []
        mesh["colors"] = []
        binary = BinaryMesh()
        binary.metadata = binary.json_mesh(mesh)
        metadata, unpack, _ = decode(serialize(binary))
        self.assertEqual(unpack(metadata["normals"]), ())
        self.assertEqual(unpack(metadata["indices"]), (2, 1, 0))


class BridgeBinaryTests(unittest.TestCase):
    def test_native_representations_bypass_python_list_serialization(self):
        container = Mock()
        container.read_coordinates.return_value = 1
        container.get_molecular_representation_mesh.side_effect = lambda *_: native_mesh()
        container.get_bonds_mesh.side_effect = lambda *_: native_mesh()
        container.get_chains_in_model.return_value = ["A", "B"]
        module = SimpleNamespace(molecules_container_t=lambda *_: container)
        for representation, split in (("bonds", False), ("ribbon", False), ("surface", False), ("ribbon", True)):
            payload = {"text": "data_test", "format": "mmcif", "outputFormat": "binary",
                       "representation": representation, "splitByChain": split}
            with self.subTest(representation=representation, split=split), \
                    patch.dict(sys.modules, {"coot_headless_api": module}), \
                    patch.object(chapi_bridge, "_suppress_stdio", return_value=None), \
                    patch.object(chapi_bridge, "_mesh_to_json", side_effect=AssertionError("numeric list path used")):
                result = chapi_bridge._run_payload(payload)
            metadata, _, _ = decode(serialize(result))
            if split:
                self.assertEqual(metadata["chainIds"], ["A", "B"])
            else:
                self.assertEqual(metadata["triangleCount"], 1)

    def test_bonds_selection_packs_existing_gltf_result(self):
        container = Mock()
        container.read_coordinates.return_value = 1
        container.export_model_molecule_as_gltf.side_effect = lambda *args: Path(args[-1]).write_bytes(b"fake GLB" * 4)
        module = SimpleNamespace(molecules_container_t=lambda *_: container)
        gltf_mesh = chapi_bridge._mesh_to_json(native_mesh())
        with patch.dict(sys.modules, {"coot_headless_api": module}), \
                patch.object(chapi_bridge, "_suppress_stdio", return_value=None), \
                patch.object(chapi_bridge, "_parse_gltf_glb", return_value=gltf_mesh):
            result = chapi_bridge._run_payload({"text": "data_test", "format": "mmcif",
                "representation": "bonds-selection", "outputFormat": "binary"})
        metadata, unpack, _ = decode(serialize(result))
        self.assertEqual(unpack(metadata["indices"]), (2, 1, 0))
        self.assertEqual(unpack(metadata["positions"]), (1.25, -2.5, 3) * 3)

    def test_leased_source_file_is_read_directly_and_never_deleted(self):
        container = Mock()
        container.read_coordinates.return_value = 1
        container.get_bonds_mesh.return_value = native_mesh()
        module = SimpleNamespace(molecules_container_t=lambda *_: container)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "source.cif"
            path.write_text("data_test")
            with patch.dict(sys.modules, {"coot_headless_api": module}), \
                    patch.object(chapi_bridge, "_suppress_stdio", return_value=None), \
                    patch.object(chapi_bridge, "_write_temp_structure", side_effect=AssertionError("copied source")):
                chapi_bridge._run_payload({"sourcePath": str(path), "format": "mmcif", "outputFormat": "binary"})
            container.read_coordinates.assert_called_once_with(str(path))
            self.assertEqual(path.read_text(), "data_test")


FAKE_COOT = '''
from types import SimpleNamespace
class molecules_container_t:
    def __init__(self, *_): pass
    def read_coordinates(self, *_): return 0
    def get_bonds_mesh(self, *_):
        return SimpleNamespace(vertices=[SimpleNamespace(pos=(1,2,3),normal=(0,0,1),color=(1,0,0,1)) for _ in range(3)],
                               triangles=[SimpleNamespace(point_id=(0,1,2))],name="test",status=1)
'''


class BridgeCliTests(unittest.TestCase):
    def run_bridge(self, text, *, server=False):
        with tempfile.TemporaryDirectory() as directory:
            (Path(directory) / "coot_headless_api.py").write_text(FAKE_COOT)
            environment = dict(os.environ, PYTHONPATH=directory)
            return subprocess.run([sys.executable, chapi_bridge.__file__] + (["--server"] if server else []),
                                  input=text.encode(), capture_output=True, env=environment, timeout=10)

    def test_default_cli_remains_json_and_binary_cli_is_packed(self):
        payload = {"text": "data_test", "format": "mmcif"}
        original = self.run_bridge(json.dumps(payload))
        self.assertEqual(original.returncode, 0, original.stderr)
        self.assertIsInstance(json.loads(original.stdout)["positions"], list)
        binary = self.run_bridge(json.dumps(dict(payload, outputFormat="binary")))
        self.assertEqual(binary.returncode, 0, binary.stderr)
        metadata, unpack, _ = decode(binary.stdout)
        self.assertEqual(unpack(metadata["positions"]), tuple(json.loads(original.stdout)["positions"]))

    def test_bad_input_has_nonzero_exit_and_no_binary_success_body(self):
        result = self.run_bridge(json.dumps({"outputFormat": "binary"}))
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, b"")
        self.assertIn(b"Payload must include", result.stderr)

    def test_server_binary_length_preserves_following_json_and_error_frames(self):
        payload = {"text": "data_test", "format": "mmcif", "outputFormat": "binary"}
        lines = [json.dumps(payload), json.dumps({"op": "ping"}), "malformed",
                 json.dumps(dict(payload, outputFormat="json"))]
        result = self.run_bridge("\n".join(lines) + "\n", server=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        stream = io.BytesIO(result.stdout)
        prefix, length = stream.readline().split(b"\t")
        self.assertEqual(prefix, b"BINARY")
        metadata, _, _ = decode(stream.read(int(length)))
        self.assertEqual(metadata["triangleCount"], 1)
        self.assertEqual(stream.readline(), b'OK\t{"pong":true}\n')
        self.assertTrue(stream.readline().startswith(b'ERR\t{"error":"Invalid JSON input:'))
        prefix, body = stream.readline().split(b"\t", 1)
        self.assertEqual(prefix, b"OK")
        self.assertIsInstance(json.loads(body)["positions"], list)
        self.assertEqual(stream.read(), b"")


if __name__ == "__main__":
    unittest.main()
