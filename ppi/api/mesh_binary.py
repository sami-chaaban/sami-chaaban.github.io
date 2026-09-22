"""ROAMI mesh envelope with little-endian packed arrays and small JSON metadata.

Wire layout: b'ROAMIM01', uint32 metadata length, UTF-8 JSON, zero padding
to a four-byte boundary, then packed buffers. Descriptor offsets are relative
to the start of the packed buffers, and counts are numbers of elements.
"""
from __future__ import annotations

from array import array
import json
import math
import struct
import sys
from typing import Any, BinaryIO, Dict, Iterable


MAGIC = b"ROAMIM01"
FLOAT32_MAX = 3.4028234663852886e38


def _float32(value: Any, fallback: float = 0.0) -> float:
    try:
        number = float(value)
    except (ValueError, TypeError, OverflowError):
        return fallback
    return number if math.isfinite(number) and abs(number) <= FLOAT32_MAX else fallback


def _indices(values: Iterable[Any], vertex_count: int) -> Iterable[int]:
    for value in values:
        index = int(value)
        if index != value or index < 0 or index >= vertex_count or index > 0xFFFFFFFF:
            raise ValueError("Mesh triangle index is outside the vertex buffer.")
        yield index


class BinaryMesh:
    """Collect packed buffers without materializing Python numeric lists.

    Streaming each buffer avoids another full-size joined response allocation
    in the native worker. Arrays live only as long as this response.
    """

    def __init__(self) -> None:
        self.metadata: Dict[str, Any] = {}
        self._buffers: list[array] = []
        self._buffer_bytes = 0

    def _append(self, values: Iterable[Any], kind: str) -> Dict[str, Any]:
        packed = array("f" if kind == "float32" else "I", values)
        if packed.itemsize != 4:
            raise RuntimeError("Mesh binary output requires four-byte floats and unsigned ints.")
        if sys.byteorder != "little":
            packed.byteswap()
        descriptor = {"offset": self._buffer_bytes, "count": len(packed), "type": kind}
        self._buffers.append(packed)
        self._buffer_bytes += len(packed) * 4
        return descriptor

    def _check_shape(self, mesh: Dict[str, Any]) -> Dict[str, Any]:
        vertices = mesh["vertexCount"]
        triangles = mesh["triangleCount"]
        if mesh["positions"]["count"] != vertices * 3:
            raise ValueError("Mesh position count does not match its vertices.")
        if mesh["normals"]["count"] not in (0, vertices * 3):
            raise ValueError("Mesh normal count does not match its vertices.")
        if mesh["colors"]["count"] not in (0, vertices * 4):
            raise ValueError("Mesh color count does not match its vertices.")
        if mesh["indices"]["count"] != triangles * 3:
            raise ValueError("Mesh index count does not match its triangles.")
        return mesh

    def native_mesh(self, mesh: Any) -> Dict[str, Any]:
        vertices = getattr(mesh, "vertices", None) or []
        triangles = getattr(mesh, "triangles", None) or []
        vertex_count = len(vertices)
        # Each native vector property is fetched once, rather than once per
        # component. Generator -> array never constructs Python numeric lists.
        return self._check_shape({
            "positions": self._append((_float32(x) for v in vertices for x in v.pos), "float32"),
            "normals": self._append((_float32(x) for v in vertices for x in v.normal), "float32"),
            "colors": self._append((_float32(x, 1.0 if i == 3 else 0.0)
                                    for v in vertices for i, x in enumerate(v.color)), "float32"),
            "indices": self._append(_indices((x for t in triangles for x in t.point_id), vertex_count), "uint32"),
            "vertexCount": vertex_count,
            "triangleCount": len(triangles),
            "name": getattr(mesh, "name", ""),
            "status": getattr(mesh, "status", None),
        })

    def json_mesh(self, mesh: Dict[str, Any]) -> Dict[str, Any]:
        """Pack the existing bonds-selection glTF decoder's output."""
        result = dict(mesh)
        for field in ("positions", "normals", "colors"):
            result[field] = self._append((_float32(x, 1.0 if field == "colors" and i % 4 == 3 else 0.0)
                                          for i, x in enumerate(mesh[field])), "float32")
        result["indices"] = self._append(_indices(mesh["indices"], mesh["vertexCount"]), "uint32")
        return self._check_shape(result)

    def header(self) -> bytes:
        metadata = json.dumps(self.metadata, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode("utf-8")
        if len(metadata) > 0xFFFFFFFF:
            raise ValueError("Mesh metadata exceeds the binary format limit.")
        return MAGIC + struct.pack("<I", len(metadata)) + metadata + b"\0" * (-len(metadata) % 4)

    def write_to(self, stream: BinaryIO, *, framed: bool = False) -> None:
        header = self.header()
        if framed:
            stream.write(f"BINARY\t{len(header) + self._buffer_bytes}\n".encode("ascii"))
        stream.write(header)
        for packed in self._buffers:
            stream.write(memoryview(packed).cast("B"))
        stream.flush()
