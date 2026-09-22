export const MESH_BINARY_CONTENT_TYPE = 'application/vnd.roami.mesh';
const MAGIC = 'ROAMIM01';

// Metadata holds the existing mesh shape; array descriptors refer to the packed
// data section. Views share the response buffer instead of expanding JSON numbers.
export function decodeMeshBinary(buffer) {
  if (!(buffer instanceof ArrayBuffer) || buffer.byteLength < 12) {
    throw new Error('Truncated binary mesh header.');
  }
  const header = new DataView(buffer);
  for (let i = 0; i < MAGIC.length; i += 1) {
    if (header.getUint8(i) !== MAGIC.charCodeAt(i)) throw new Error('Invalid binary mesh magic.');
  }
  const metadataLength = header.getUint32(8, true);
  const dataOffset = Math.ceil((12 + metadataLength) / 4) * 4;
  if (dataOffset > buffer.byteLength) throw new Error('Truncated binary mesh metadata.');
  let metadata;
  try {
    metadata = JSON.parse(new TextDecoder('utf-8', { fatal: true }).decode(new Uint8Array(buffer, 12, metadataLength)));
  } catch (_) { throw new Error('Invalid binary mesh metadata.'); }
  if (!metadata || typeof metadata !== 'object' || Array.isArray(metadata)) {
    throw new Error('Invalid binary mesh payload.');
  }
  const littleEndian = new Uint8Array(new Uint32Array([1]).buffer)[0] === 1;
  function array(descriptor, field) {
    const expected = field === 'indices' ? 'uint32' : 'float32';
    if (!descriptor || typeof descriptor !== 'object' || descriptor.type !== expected ||
        !Number.isSafeInteger(descriptor.offset) || descriptor.offset < 0 || descriptor.offset % 4 ||
        !Number.isSafeInteger(descriptor.count) || descriptor.count < 0 ||
        descriptor.count > Math.floor((buffer.byteLength - dataOffset - descriptor.offset) / 4)) {
      throw new Error(`Invalid binary mesh ${field} descriptor.`);
    }
    const Type = expected === 'uint32' ? Uint32Array : Float32Array;
    const offset = dataOffset + descriptor.offset;
    if (littleEndian) return new Type(buffer, offset, descriptor.count);
    const values = new Type(descriptor.count);
    for (let i = 0; i < values.length; i += 1) {
      values[i] = header[expected === 'uint32' ? 'getUint32' : 'getFloat32'](offset + i * 4, true);
    }
    return values;
  }
  function mesh(value) {
    if (!value || typeof value !== 'object' || Array.isArray(value)) throw new Error('Invalid binary mesh entry.');
    for (const field of ['positions', 'normals', 'colors', 'indices']) {
      value[field] = array(value[field], field);
    }
    const vertices = value.positions?.length / 3;
    if (!Number.isSafeInteger(vertices) || vertices < 0 || value.vertexCount !== vertices ||
        (value.normals?.length && value.normals.length !== vertices * 3) ||
        (value.colors?.length && value.colors.length !== vertices * 4) ||
        (value.indices && (value.indices.length % 3 || value.triangleCount !== value.indices.length / 3))) {
      throw new Error('Inconsistent binary mesh array counts.');
    }
    if (value.indices) {
      for (const index of value.indices) if (index >= vertices) throw new Error('Binary mesh index is out of range.');
    }
    return value;
  }
  if (Object.hasOwn(metadata, 'meshes')) {
    if (!Array.isArray(metadata.meshes)) throw new Error('Invalid binary mesh list.');
    metadata.meshes.forEach(mesh);
    return metadata;
  }
  return mesh(metadata);
}
