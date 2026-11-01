import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { depthAngle, depthBlend, keyframePriority, orbitCamera, warpMatrix } from '../src/lib/hero-depth.ts';

const camera = JSON.parse(readFileSync(new URL('../assets/hero-depth/camera.json', import.meta.url), 'utf8'));
const desktop = { left: 0, top: 0, width: camera.width, height: camera.height };
const mobile = { left: 555, top: 0, width: 960, height: camera.height };

// The vertex shader's mapping from a source pixel centre to target view pixels.
function project(view, matrix, x, y, inverseDepth) {
  const fx = x + 0.5 + view.left;
  const fy = y + 0.5 + view.top;
  const ndc = [fx / camera.width * 2 - 1, 1 - fy / camera.height * 2];
  const p = camera.projection;
  const point = [(ndc[0] + p[8]) / p[0], (ndc[1] + p[9]) / p[5], -1, inverseDepth];
  const clip = [0, 1, 2, 3].map(row => point.reduce((sum, value, column) => sum + matrix[column * 4 + row] * value, 0));
  return [(clip[0] / clip[3] + 1) / 2 * view.width - 0.5, (1 - clip[1] / clip[3]) / 2 * view.height - 0.5];
}

const near = camera.invMax;
const far = camera.invMax / 40;

test('a keyframe drawn at its own angle lands on its own pixels, cropped or not', () => {
  for (const view of [desktop, mobile]) {
    const matrix = warpMatrix(camera, view, 1.2, 1.2);
    for (const [x, y, inverse] of [[0, 0, near], [400, 300, far], [view.width - 1, view.height - 1, 0]]) {
      // Float32 matrices with scene-scale translations leave about a thousandth of a pixel.
      const [u, v] = project(view, matrix, x, y, inverse);
      assert.ok(Math.abs(u - x) < 1e-2 && Math.abs(v - y) < 1e-2, `${x},${y} moved to ${u},${v}`);
    }
  }
});

test('the mobile crop shows the same warp as the full frame, shifted by the crop', () => {
  const full = warpMatrix(camera, desktop, 1.2, 1.5);
  const cropped = warpMatrix(camera, mobile, 1.2, 1.5);
  for (const [x, y, inverse] of [[100, 200, near], [700, 50, far], [900, 900, 0]]) {
    const [fu, fv] = project(desktop, full, x + mobile.left, y, inverse);
    const [cu, cv] = project(mobile, cropped, x, y, inverse);
    assert.ok(Math.abs(fu - mobile.left - cu) < 1e-2 && Math.abs(fv - cv) < 1e-2);
  }
});

test('turning the camera moves near surfaces further than far ones', () => {
  const matrix = warpMatrix(camera, desktop, 1.2, 1.4);
  const nearShift = project(desktop, matrix, 1000, 400, near)[0] - 1000;
  const farShift = project(desktop, matrix, 1000, 400, far)[0] - 1000;
  assert.ok(Math.abs(nearShift - farShift) > 1, `parallax ${nearShift} vs ${farShift}`);
});

test('the orbit turns about world Y through the pivot', () => {
  const start = orbitCamera(camera, 0);
  const turned = orbitCamera(camera, 3.857);
  const distance = (m) => Math.hypot(m[12] - camera.pivot[0], m[14] - camera.pivot[2]);
  assert.ok(start.every((value, i) => Math.abs(value - camera.camera[i]) < 1e-4));
  assert.ok(Math.abs(turned[13] - start[13]) < 1e-3);
  assert.ok(Math.abs(distance(turned) - distance(start)) < 1e-3);
  assert.ok(Math.abs(turned[12] - start[12]) > 1);
});

const frames = [0, 4, 8].map(position => ({ position, angle: (8 - position) * camera.degreesPerFrame, color: `c${position}`, depth: `d${position}` }));

test('playback blends the loaded keyframes either side of a position', () => {
  assert.deepEqual(depthBlend(5, frames, () => true), { a: 1, b: 2, mix: 0.25 });
  assert.deepEqual(depthBlend(4, frames, () => true), { a: 1, b: 1, mix: 0 });
  assert.deepEqual(depthBlend(5, frames, index => index !== 2), { a: 1, b: 1, mix: 0 });
  assert.deepEqual(depthBlend(5, frames, index => index === 0), { a: 0, b: 0, mix: 0 });
  assert.equal(depthBlend(5, frames, () => false), undefined);
});

test('keyframes load nearest first and angles run backwards with position', () => {
  assert.deepEqual(keyframePriority(7, frames), [2, 1, 0]);
  assert.equal(depthAngle(8, frames, camera), 0);
  assert.ok(Math.abs(depthAngle(2, frames, camera) - 6 * camera.degreesPerFrame) < 1e-9);
});
