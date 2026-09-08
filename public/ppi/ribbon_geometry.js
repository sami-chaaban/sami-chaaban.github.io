import * as THREE from './vendor/three/three.module.min.js';

function applySurfaceScatterMaterial(material, settings = {}) {
  if (!material) return;
  const defaults = {
    enabled: false,
    roughnessBoost: 0.22,
  };
  const scatter = {
    enabled: Boolean(settings.enabled),
    roughnessBoost: Number.isFinite(settings.roughnessBoost)
      ? settings.roughnessBoost
      : defaults.roughnessBoost,
  };
  material.userData = material.userData || {};
  material.userData.surfaceScatter = scatter;
  if (!scatter.enabled) {
    return;
  }
  const previous = material.onBeforeCompile;
  material.onBeforeCompile = (shader) => {
    if (previous) previous(shader);
    shader.uniforms.surfaceScatterEnabled = { value: scatter.enabled ? 1.0 : 0.0 };
    shader.uniforms.surfaceScatterStrength = { value: scatter.roughnessBoost };
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <common>',
      `#include <common>\n            uniform float surfaceScatterEnabled;\n            uniform float surfaceScatterStrength;`
    );
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <normal_fragment_maps>',
      `#include <normal_fragment_maps>\n            if (surfaceScatterEnabled > 0.5) {\n              vec3 dxy = max(abs(dFdx(normal)), abs(dFdy(normal)));\n              float geometryRoughness = max(max(dxy.x, dxy.y), dxy.z);\n              roughnessFactor = clamp(\n                max(roughnessFactor, geometryRoughness * max(surfaceScatterStrength, 0.0)),\n                0.045,\n                1.0\n              );\n            }`
    );
    material.userData.shader = shader;
  };
}

function applyHalftoneMaterial(material, settings = {}) {
  if (!material) return;
  const defaults = {
    enabled: false,
    scale: 7.5,
    strength: 0.9,
    angle: 0.55,
    softness: 0.1,
  };
  const halftone = {
    enabled: Boolean(settings.enabled),
    scale: Number.isFinite(settings.scale) ? settings.scale : defaults.scale,
    strength: Number.isFinite(settings.strength) ? settings.strength : defaults.strength,
    angle: Number.isFinite(settings.angle) ? settings.angle : defaults.angle,
    softness: Number.isFinite(settings.softness) ? settings.softness : defaults.softness,
  };
  material.userData = material.userData || {};
  material.userData.halftone = halftone;
  if (!halftone.enabled) {
    return;
  }
  const previous = material.onBeforeCompile;
  material.onBeforeCompile = (shader) => {
    if (previous) previous(shader);
    shader.uniforms.halftoneEnabled = { value: halftone.enabled ? 1.0 : 0.0 };
    shader.uniforms.halftoneScale = { value: halftone.scale };
    shader.uniforms.halftoneStrength = { value: halftone.strength };
    shader.uniforms.halftoneAngle = { value: halftone.angle };
    shader.uniforms.halftoneSoftness = { value: halftone.softness };
    shader.vertexShader = shader.vertexShader.replace(
      '#include <common>',
      `#include <common>\n            varying vec3 vWorldPos;\n            varying vec3 vWorldNormal;`
    );
    shader.vertexShader = shader.vertexShader.replace(
      '#include <beginnormal_vertex>',
      `#include <beginnormal_vertex>\n            vWorldNormal = normalize(mat3(modelMatrix) * objectNormal);`
    );
    shader.vertexShader = shader.vertexShader.replace(
      '#include <worldpos_vertex>',
      `#include <worldpos_vertex>\n            vWorldPos = worldPosition.xyz;`
    );
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <common>',
      `#include <common>\n            uniform float halftoneEnabled;\n            uniform float halftoneScale;\n            uniform float halftoneStrength;\n            uniform float halftoneAngle;\n            uniform float halftoneSoftness;\n            varying vec3 vWorldPos;\n            varying vec3 vWorldNormal;`
    );
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <output_fragment>',
      `if (halftoneEnabled > 0.5) {\n              vec3 n = normalize(vWorldNormal);\n              vec3 absN = abs(n);\n              vec2 proj;\n              if (absN.x > absN.y && absN.x > absN.z) {\n                proj = vWorldPos.zy;\n              } else if (absN.y > absN.z) {\n                proj = vWorldPos.xz;\n              } else {\n                proj = vWorldPos.xy;\n              }\n              vec3 halftoneLightDir = normalize(vec3(0.35, 0.55, 0.85));\n              float geomLit = clamp(dot(n, halftoneLightDir), 0.0, 1.0);\n              float geomShade = pow(1.0 - geomLit, 1.15);\n              float luma = clamp(dot(outgoingLight, vec3(0.2126, 0.7152, 0.0722)), 0.0, 1.0);\n              float lightShade = pow(1.0 - luma, 1.1);\n              float shade = clamp(mix(lightShade, geomShade, 0.65), 0.0, 1.0);\n              float strength = clamp(halftoneStrength, 0.0, 1.0);\n              float radius = shade * 0.5 * strength;\n              float softness = max(0.001, halftoneSoftness);\n              float angle = halftoneAngle;\n              float c = cos(angle);\n              float s = sin(angle);\n              mat2 rot = mat2(c, -s, s, c);\n              vec2 grid = (rot * proj) / max(0.0001, halftoneScale);\n              vec2 cell = fract(grid) - 0.5;\n              float dist = length(cell);\n              float dotMask = 1.0 - smoothstep(radius, radius + softness, dist);\n              dotMask *= shade * strength;\n              outgoingLight *= (1.0 - dotMask * 0.85);\n            }\n            #include <output_fragment>`
    );
    material.userData.shader = shader;
  };
}

function makeRibbonSection(segments = 18) {
  const pts = [];
  const steps = Math.max(6, segments);
  for (let i = 0; i <= steps; i += 1) {
    const t = Math.PI - (i / steps) * Math.PI;
    const x = Math.cos(t);
    const y = -0.5 + Math.sin(t);
    pts.push([x, y]);
  }
  return pts;
}

function makeStripSection() {
  return [
    [-1, -1],
    [1, -1],
    [1, 1],
    [-1, 1],
  ];
}

function sweepGeometry(samples, section2D, closedRing, widthScale, heightScale) {
  const pos = samples.pos;
  const nrmN = samples.n;
  const nrmB = samples.b;
  const wArr = samples.w;
  const hArr = samples.h;

  const K = wArr.length;
  const S = section2D.length;

  const vertices = new Float32Array(K * S * 3);
  const uvs = new Float32Array(K * S * 2);
  const wrap = closedRing ? S : S - 1;
  const quadCount = Math.max(0, K - 1) * Math.max(0, wrap);
  const indexCount = quadCount * 6;
  const indexArray =
    K * S > 65535
      ? new Uint32Array(indexCount)
      : new Uint16Array(indexCount);

  for (let k = 0; k < K; k += 1) {
    const px = pos[3 * k + 0];
    const py = pos[3 * k + 1];
    const pz = pos[3 * k + 2];
    const nx = nrmN[3 * k + 0];
    const ny = nrmN[3 * k + 1];
    const nz = nrmN[3 * k + 2];
    const bx = nrmB[3 * k + 0];
    const by = nrmB[3 * k + 1];
    const bz = nrmB[3 * k + 2];
    const w = wArr[k] * widthScale;
    const h = hArr[k] * heightScale;

    for (let s = 0; s < S; s += 1) {
      const u = section2D[s][0] * w;
      const v = section2D[s][1] * h;

      const x = px + u * nx + v * bx;
      const y = py + u * ny + v * by;
      const z = pz + u * nz + v * bz;

      const base = (k * S + s) * 3;
      vertices[base + 0] = x;
      vertices[base + 1] = y;
      vertices[base + 2] = z;

      const uvBase = (k * S + s) * 2;
      const uvU = S > 1 ? s / (S - 1) : 0;
      const uvV = K > 1 ? k / (K - 1) : 0;
      uvs[uvBase + 0] = uvU;
      uvs[uvBase + 1] = uvV;
    }
  }

  let iOff = 0;
  for (let k = 0; k < K - 1; k += 1) {
    for (let s = 0; s < wrap; s += 1) {
      const s2 = s + 1;
      const a = k * S + s;
      const b = k * S + (closedRing ? s2 % S : s2);
      const c = (k + 1) * S + s;
      const d = (k + 1) * S + (closedRing ? s2 % S : s2);
      indexArray[iOff + 0] = a;
      indexArray[iOff + 1] = c;
      indexArray[iOff + 2] = b;
      indexArray[iOff + 3] = b;
      indexArray[iOff + 4] = c;
      indexArray[iOff + 5] = d;
      iOff += 6;
    }
  }

  const geom = new THREE.BufferGeometry();
  geom.setAttribute('position', new THREE.BufferAttribute(vertices, 3));
  geom.setAttribute('uv', new THREE.BufferAttribute(uvs, 2));
  geom.setIndex(new THREE.BufferAttribute(indexArray, 1));
  geom.computeVertexNormals();
  return geom;
}

function buildSectionShape(section2D, width, height) {
  const shape = new THREE.Shape();
  if (!section2D.length) return shape;
  shape.moveTo(section2D[0][0] * width, section2D[0][1] * height);
  for (let i = 1; i < section2D.length; i += 1) {
    shape.lineTo(section2D[i][0] * width, section2D[i][1] * height);
  }
  shape.closePath();
  return shape;
}

function buildCapGeometry(section2D, center, normal, binormal, width, height) {
  const shape = buildSectionShape(section2D, width, height);
  const geom = new THREE.ShapeGeometry(shape, 8);
  const pos = geom.attributes.position;
  const temp = new THREE.Vector3();
  for (let i = 0; i < pos.count; i += 1) {
    const x = pos.getX(i);
    const y = pos.getY(i);
    temp.copy(center).addScaledVector(normal, x).addScaledVector(binormal, y);
    pos.setXYZ(i, temp.x, temp.y, temp.z);
  }
  geom.computeVertexNormals();
  return geom;
}

function sliceSamples(samples, start, end) {
  const sliceView = (arr, from, to) => (
    ArrayBuffer.isView(arr) && typeof arr.subarray === 'function'
      ? arr.subarray(from, to)
      : arr.slice(from, to)
  );
  return {
    pos: sliceView(samples.pos, 3 * start, 3 * end),
    n: sliceView(samples.n, 3 * start, 3 * end),
    b: sliceView(samples.b, 3 * start, 3 * end),
    w: sliceView(samples.w, start, end),
    h: sliceView(samples.h, start, end),
    ss: samples.ss.slice(start, end),
  };
}

function ensureTypedSamples(samples) {
  if (!samples) {
    return {
      pos: new Float32Array(0),
      n: new Float32Array(0),
      b: new Float32Array(0),
      w: new Float32Array(0),
      h: new Float32Array(0),
      ss: [],
    };
  }
  if (samples.__typedRibbonSamples) {
    return samples.__typedRibbonSamples;
  }
  const typed = {
    pos: samples.pos instanceof Float32Array ? samples.pos : new Float32Array(samples.pos || []),
    n: samples.n instanceof Float32Array ? samples.n : new Float32Array(samples.n || []),
    b: samples.b instanceof Float32Array ? samples.b : new Float32Array(samples.b || []),
    w: samples.w instanceof Float32Array ? samples.w : new Float32Array(samples.w || []),
    h: samples.h instanceof Float32Array ? samples.h : new Float32Array(samples.h || []),
    ss: Array.isArray(samples.ss) ? samples.ss : [],
  };
  if (Object.isExtensible(samples)) {
    Object.defineProperty(samples, '__typedRibbonSamples', {
      value: typed,
      enumerable: false,
      configurable: true,
    });
  }
  return typed;
}

function createOutlineShellGeometry(sourceGeometry, thickness = 0.08) {
  if (!sourceGeometry) return null;
  const geometry = sourceGeometry.clone();
  let normal = geometry.getAttribute('normal');
  if (!normal) {
    geometry.computeVertexNormals();
    normal = geometry.getAttribute('normal');
  }
  const position = geometry.getAttribute('position');
  if (!position || !normal || position.count !== normal.count) {
    return geometry;
  }
  const posArr = position.array;
  const nrmArr = normal.array;
  for (let i = 0; i < posArr.length; i += 3) {
    posArr[i + 0] += nrmArr[i + 0] * thickness;
    posArr[i + 1] += nrmArr[i + 1] * thickness;
    posArr[i + 2] += nrmArr[i + 2] * thickness;
  }
  position.needsUpdate = true;
  geometry.computeBoundingBox();
  geometry.computeBoundingSphere();
  return geometry;
}

export function buildCartoonGroup(ribbonJSON, options = {}) {
  const group = new THREE.Group();

  const helixSection = makeRibbonSection(18);
  const coilSection = makeRibbonSection(14);
  const strandSection = makeStripSection();

  const widthScale = 0.65 * (Number.isFinite(options.widthScale) ? options.widthScale : 1);
  const heightScale = 0.75 * (Number.isFinite(options.heightScale) ? options.heightScale : 1);
  const chainMaterials = new Map();
  const halftoneSettings = options.halftone || {};
  const surfaceScatterSettings = options.surfaceScatter || {};
  const chainData = new Map();
  const chainMeshes = [];
  const materialDefaults = {
    roughness: 0.18,
    metalness: 0.1,
    clearcoat: 0.78,
    clearcoatRoughness: 0.16,
    envMapIntensity: 1.1,
  };
  const palette = [
    '#6B8FB5',
    '#B0727E',
    '#6B9E7C',
    '#8A76B0',
    '#5E9BA5',
    '#A8745A',
    '#7D965B',
    '#6E7FB8',
  ];

  for (const [chainIndex, chain] of (ribbonJSON.chains || []).entries()) {
    const chainId = chain.id || String(chainIndex);
    const chainGroup = new THREE.Group();
    chainGroup.userData.chainId = chainId;
    group.add(chainGroup);
    if (!chainMaterials.has(chainId)) {
      const mat = new THREE.MeshPhysicalMaterial({
        roughness: materialDefaults.roughness,
        metalness: materialDefaults.metalness,
        clearcoat: materialDefaults.clearcoat,
        clearcoatRoughness: materialDefaults.clearcoatRoughness,
        envMapIntensity: materialDefaults.envMapIntensity,
      });
      mat.color.set(palette[chainIndex % palette.length]);
      mat.transparent = false;
      mat.opacity = 1.0;
      mat.side = THREE.FrontSide;
      mat.depthWrite = true;
      mat.dithering = true;
      applySurfaceScatterMaterial(mat, surfaceScatterSettings);
      applyHalftoneMaterial(mat, halftoneSettings);
      chainMaterials.set(chainId, mat);
    }
    const chainMat = chainMaterials.get(chainId);
    if (!chainData.has(chainId)) {
      const outlineMat = new THREE.MeshBasicMaterial({
        color: 0xffffff,
        transparent: true,
        opacity: 0.0,
        depthWrite: true,
        depthTest: true,
        fog: false,
        toneMapped: false,
        side: THREE.BackSide,
      });
      const capMaterial = chainMat.clone();
      capMaterial.side = THREE.DoubleSide;
      capMaterial.onBeforeCompile = null;
      applySurfaceScatterMaterial(capMaterial, surfaceScatterSettings);
      applyHalftoneMaterial(capMaterial, halftoneSettings);
      chainData.set(chainId, {
        id: chainId,
        material: chainMat,
        capMaterial,
        outlineMaterial: outlineMat,
        baseColor: chainMat.color.clone(),
        meshes: [],
        outlines: [],
        coreMeshes: [],
        capMeshes: [],
        group: chainGroup,
      });
    }
    const chainEntry = chainData.get(chainId);
    for (const seg of chain.segments || []) {
      const s = ensureTypedSamples(seg.samples);
      const K = s.w.length;
      let start = 0;
      const pad = 2;
      while (start < K) {
        const kind = s.ss[start];
        let end = start + 1;
        while (end < K && s.ss[end] === kind) end += 1;

        const runStart = start === 0 ? start : Math.max(0, start - pad);
        const runEnd = end === K ? end : Math.min(K, end + pad);

        if (runEnd - runStart >= 4) {
          const sub = sliceSamples(s, runStart, runEnd);

          let geom;
          let section = coilSection;
          let closedRing = true;
          if (kind === 'H') {
            section = helixSection;
            closedRing = true;
          } else if (kind === 'E') {
            section = strandSection;
            closedRing = false;
          }
          geom = sweepGeometry(sub, section, closedRing, widthScale, heightScale);

          const outlineGeom = createOutlineShellGeometry(geom, 0.08);
          if (outlineGeom) {
            const outline = new THREE.Mesh(outlineGeom, chainEntry.outlineMaterial);
            outline.userData.chainId = chainId;
            outline.userData.isOutline = true;
            outline.renderOrder = 3;
            outline.frustumCulled = false;
            outline.castShadow = false;
            outline.receiveShadow = false;
            chainGroup.add(outline);
            chainEntry.outlines.push(outline);
          }

          const mesh = new THREE.Mesh(geom, chainMat);
          mesh.userData.chainId = chainId;
          mesh.userData.isRibbonSurface = true;
          mesh.renderOrder = 1;
          mesh.frustumCulled = false;
          mesh.castShadow = true;
          mesh.receiveShadow = true;
          chainGroup.add(mesh);
          chainMeshes.push(mesh);
          if (chainEntry) {
            chainEntry.meshes.push(mesh);
            chainEntry.coreMeshes.push(mesh);
          }

          const startIdx = 0;
          const endIdx = sub.w.length - 1;
          const startCenter = new THREE.Vector3(sub.pos[0], sub.pos[1], sub.pos[2]);
          const endCenter = new THREE.Vector3(
            sub.pos[3 * endIdx],
            sub.pos[3 * endIdx + 1],
            sub.pos[3 * endIdx + 2]
          );
          const startNormal = new THREE.Vector3(sub.n[0], sub.n[1], sub.n[2]).normalize();
          const startBinormal = new THREE.Vector3(sub.b[0], sub.b[1], sub.b[2]).normalize();
          const endNormal = new THREE.Vector3(
            sub.n[3 * endIdx],
            sub.n[3 * endIdx + 1],
            sub.n[3 * endIdx + 2]
          ).normalize();
          const endBinormal = new THREE.Vector3(
            sub.b[3 * endIdx],
            sub.b[3 * endIdx + 1],
            sub.b[3 * endIdx + 2]
          ).normalize();
          const startW = sub.w[startIdx] * widthScale;
          const startH = sub.h[startIdx] * heightScale;
          const endW = sub.w[endIdx] * widthScale;
          const endH = sub.h[endIdx] * heightScale;
          const capGeomStart = buildCapGeometry(section, startCenter, startNormal, startBinormal, startW, startH);
          const capGeomEnd = buildCapGeometry(section, endCenter, endNormal, endBinormal, endW, endH);
          const capStart = new THREE.Mesh(capGeomStart, chainEntry.capMaterial);
          capStart.userData.chainId = chainId;
          capStart.userData.isRibbonCap = true;
          capStart.renderOrder = 1;
          capStart.frustumCulled = false;
          capStart.castShadow = true;
          capStart.receiveShadow = true;
          const capEnd = new THREE.Mesh(capGeomEnd, chainEntry.capMaterial);
          capEnd.userData.chainId = chainId;
          capEnd.userData.isRibbonCap = true;
          capEnd.renderOrder = 1;
          capEnd.frustumCulled = false;
          capEnd.castShadow = true;
          capEnd.receiveShadow = true;
          chainGroup.add(capStart);
          chainGroup.add(capEnd);
          chainEntry.meshes.push(capStart, capEnd);
          chainEntry.capMeshes.push(capStart, capEnd);
        }

        start = end;
      }
    }
  }

  group.userData.chainData = chainData;
  group.userData.chainMeshes = chainMeshes;
  return group;
}
