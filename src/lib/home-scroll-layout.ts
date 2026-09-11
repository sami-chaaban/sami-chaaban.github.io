// Each stop describes the same point in the story before and after a resize:
// opening, research arrival/hold, division arrival/hold, profile arrival, and page end.
export function remapHomeScroll(position: number, previous: number[], next: number[]): number {
  if (position <= previous[0]) return next[0];
  const last = previous.length - 1;
  if (position >= previous[last] - 1) return next[last];
  for (let i = 1; i <= last; i++) {
    const travel = previous[i] - previous[i - 1];
    if (travel > 0 && position <= previous[i]) {
      const progress = (position - previous[i - 1]) / travel;
      return next[i - 1] + progress * (next[i] - next[i - 1]);
    }
  }
  return next[last];
}
