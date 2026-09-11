import { clamp, smoothstep } from './motion-math.ts';

const animationEnd = 0.98;

// One reversible scroll timeline: title, molecular interlude, research arrival.
export function homeOpening(progress: number, reducedMotion = false) {
  const t = clamp(progress);
  const titleOpacity = 1 - smoothstep(0.03, 0.22, t);
  const frameProgress = clamp(progress / animationEnd);
  return {
    titleOpacity: reducedMotion ? 1 : titleOpacity,
    frameProgress: reducedMotion ? 0 : frameProgress,
    tintProgress: reducedMotion ? 1 : smoothstep(0, 1, frameProgress),
    researchOpacity: reducedMotion ? 1 : smoothstep(0.32, 0.52, t),
  };
}
