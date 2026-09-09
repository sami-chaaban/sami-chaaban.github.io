import { smoothstep } from './motion-math.ts';

export function researchPinLayout(viewportHeight: number, contentHeight: number, reducedMotion = false) {
  const hold = reducedMotion ? 0 : viewportHeight * 0.5;
  return {
    // Tall text scrolls into view before pinning its lower edge to the viewport.
    top: reducedMotion ? 0 : Math.min(0, viewportHeight - contentHeight),
    height: contentHeight + hold,
    hold,
  };
}

export function researchFocusOpacity(progress: number, reducedMotion = false) {
  // Use the entire hold for the reveal, beginning as soon as the section pins.
  return reducedMotion ? 1 : smoothstep(0, 1, progress);
}
