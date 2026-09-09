import { clamp, smoothstep } from './motion-math.ts';

const readingFadeEnd = 0.76;

// One reversible scroll timeline: title, molecular interlude, research arrival.
export function homeOpening(progress: number, reducedMotion = false) {
  const t = clamp(progress);
  const titleOpacity = 1 - smoothstep(0.03, 0.22, t);
  const frameProgress = clamp(progress / readingFadeEnd);
  const reveal = smoothstep(0.06, 0.22, t);
  const reading = smoothstep(0.54, readingFadeEnd, t);
  return {
    titleOpacity: reducedMotion ? 1 : titleOpacity,
    frameProgress: reducedMotion ? 0 : frameProgress,
    veilOpacity: reducedMotion ? 0.9 : 0.9 - reveal * 0.6 + reading * 0.05,
    shadeOpacity: reducedMotion ? 1 : 1 - reveal + reading,
    researchOpacity: reducedMotion ? 1 : smoothstep(0.32, 0.52, t),
  };
}
