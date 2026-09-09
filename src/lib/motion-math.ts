export const clamp = (value: number, min = 0, max = 1) => Math.min(max, Math.max(min, value));

export function smoothstep(start: number, end: number, value: number): number {
  const amount = clamp((value - start) / (end - start));
  return amount * amount * (3 - 2 * amount);
}

// Equal elapsed time produces equal settling at 60 Hz and 120 Hz.
export function approach(current: number, target: number, seconds: number, timeConstant = 0.16): number {
  return current + (target - current) * (1 - Math.exp(-Math.max(0, seconds) / timeConstant));
}
