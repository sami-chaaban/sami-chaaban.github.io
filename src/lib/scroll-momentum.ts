import { approach, clamp } from './motion-math.ts';

// Ease toward the wheel's requested distance without adding extra travel.
export class ScrollMomentum {
  position = 0;
  private target = 0;

  get moving() { return this.position !== this.target; }

  push(delta: number, maximum: number) {
    // Reversing direction discards the remaining glide immediately.
    if (delta * (this.target - this.position) < 0) this.target = this.position;
    this.target = clamp(this.target + delta, 0, Math.max(0, maximum));
  }

  advance(seconds: number, maximum: number) {
    this.target = clamp(this.target, 0, Math.max(0, maximum));
    this.position = clamp(approach(this.position, this.target, seconds, 0.14), 0, Math.max(0, maximum));
    if (Math.abs(this.target - this.position) < 0.5) this.position = this.target;
    return this.position;
  }

  stop(position: number) {
    this.position = this.target = position;
  }
}
