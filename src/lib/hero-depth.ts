// Camera maths for depth-warped hero playback. Matrices are column-major, as WebGL reads them.

export type DepthCamera = {
  width: number;
  height: number;
  degreesPerFrame: number;
  // Inverse depth of the nearest surface; depth images store inverse depth / invMax.
  invMax: number;
  pivot: number[];
  projection: number[];
  // Camera-to-world at source frame 0.
  camera: number[];
};
// The part of the full render a variant shows, in full-render pixels.
export type DepthView = { left: number; top: number; width: number; height: number };
export type DepthKeyframe = { position: number; angle: number; color: string; depth: string };
export type DepthVariant = { view: DepthView; frames: DepthKeyframe[] };
export type DepthHero = { camera: DepthCamera; desktop: DepthVariant; mobile: DepthVariant };

type Matrix = Float32Array;

export function multiply(a: ArrayLike<number>, b: ArrayLike<number>): Matrix {
  const out = new Float32Array(16);
  for (let column = 0; column < 4; column++) {
    for (let row = 0; row < 4; row++) {
      let sum = 0;
      for (let k = 0; k < 4; k++) sum += a[k * 4 + row] * b[column * 4 + k];
      out[column * 4 + row] = sum;
    }
  }
  return out;
}

function rigidInverse(m: ArrayLike<number>): Matrix {
  const out = new Float32Array(16);
  for (let row = 0; row < 3; row++) for (let column = 0; column < 3; column++) out[column * 4 + row] = m[row * 4 + column];
  for (let row = 0; row < 3; row++) out[12 + row] = -(out[row] * m[12] + out[4 + row] * m[13] + out[8 + row] * m[14]);
  out[15] = 1;
  return out;
}

const translation = ([x, y, z]: number[]) => new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, x, y, z, 1]);

// The Blender camera rides an Empty that turns about world Y through the pivot.
export function orbitCamera(camera: DepthCamera, angleDegrees: number): Matrix {
  const radians = angleDegrees * Math.PI / 180;
  const c = Math.cos(radians), s = Math.sin(radians);
  const turn = new Float32Array([c, 0, -s, 0, 0, 1, 0, 0, s, 0, c, 0, 0, 0, 0, 1]);
  const toPivot = translation(camera.pivot.map(value => -value));
  return multiply(translation(camera.pivot), multiply(turn, multiply(toPivot, camera.camera)));
}

// Clip transform for points in a keyframe's camera space, seen from the target angle
// and cropped to the variant's view. Only x and y change with the crop.
export function warpMatrix(camera: DepthCamera, view: DepthView, sourceAngle: number, targetAngle: number): Matrix {
  const { width, height } = camera;
  const crop = new Float32Array([
    width / view.width, 0, 0, 0,
    0, height / view.height, 0, 0,
    0, 0, 1, 0,
    (width - 2 * view.left) / view.width - 1, 1 - (height - 2 * view.top) / view.height, 0, 1,
  ]);
  const toTarget = multiply(rigidInverse(orbitCamera(camera, targetAngle)), orbitCamera(camera, sourceAngle));
  return multiply(crop, multiply(camera.projection, toTarget));
}

// Playback runs frames in reverse, so angle falls as position rises toward the last keyframe.
export function depthAngle(position: number, frames: DepthKeyframe[], camera: DepthCamera): number {
  const end = frames[frames.length - 1];
  return end.angle + (end.position - position) * camera.degreesPerFrame;
}

// Nearest loaded keyframes either side of a position. Any loaded keyframe can be warped
// to any angle, so one side is enough while the other downloads.
export function depthBlend(position: number, frames: DepthKeyframe[], ready: (index: number) => boolean) {
  let a = -1;
  let b = -1;
  for (let index = 0; index < frames.length; index++) {
    if (!ready(index)) continue;
    if (frames[index].position <= position) a = index;
    if (frames[index].position >= position && b < 0) b = index;
  }
  if (a < 0 && b < 0) return undefined;
  if (a < 0) a = b;
  if (b < 0) b = a;
  const span = frames[b].position - frames[a].position;
  return { a, b, mix: span > 0 ? (position - frames[a].position) / span : 0 };
}

export function keyframePriority(position: number, frames: DepthKeyframe[]): number[] {
  return frames.map((_, index) => index)
    .sort((a, b) => Math.abs(frames[a].position - position) - Math.abs(frames[b].position - position));
}
