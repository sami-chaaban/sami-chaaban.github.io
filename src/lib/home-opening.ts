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

export function homeBackground(openingProgress: number, reducedMotion = false) {
  return {
    // Complete playback at research arrival and hold that frame in later sections.
    frameProgress: reducedMotion ? 0 : clamp(openingProgress),
    // The correction and dark overlay share the same opening swipe.
    tintProgress: reducedMotion ? 1 : smoothstep(0, 1, openingProgress),
  };
}

export function homeFramePosition(progress: number, frameCount: number) {
  const last = Math.max(0, frameCount - 1);
  return clamp(progress) * last;
}
