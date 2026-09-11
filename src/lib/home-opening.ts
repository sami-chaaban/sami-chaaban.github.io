import { clamp, smoothstep } from './motion-math.ts';

// The text reveal keeps its own timing within the opening transition.
export function homeOpening(progress: number, reducedMotion = false) {
  const t = clamp(progress);
  const titleOpacity = 1 - smoothstep(0.03, 0.22, t);
  return {
    titleOpacity: reducedMotion ? 1 : titleOpacity,
    researchOpacity: reducedMotion ? 1 : smoothstep(0.26, 0.46, t),
  };
}

export function homeBackground(openingProgress: number, cellArrivalProgress: number, outroArrivalProgress: number, reducedMotion = false) {
  // Each arrival plays one third, holding between transitions while the content is read.
  const frameProgress = (clamp(openingProgress) + clamp(cellArrivalProgress) + clamp(outroArrivalProgress)) / 3;
  return {
    frameProgress: reducedMotion ? 0 : frameProgress,
    tintProgress: reducedMotion ? 1 : smoothstep(0, 1, frameProgress),
  };
}

export function homeFramePosition(progress: number, frameCount: number) {
  const last = Math.max(0, frameCount - 1);
  // Both reading pauses land on real images, even when thirds fall between frames.
  const first = Math.round(last / 3);
  const second = Math.round(last * 2 / 3);
  return clamp(progress * 3) * first
    + clamp(progress * 3 - 1) * (second - first)
    + clamp(progress * 3 - 2) * (last - second);
}
