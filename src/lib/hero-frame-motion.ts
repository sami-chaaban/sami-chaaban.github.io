import { approach } from './motion-math.ts';

// Follow continuous scrolling; snap only after scrolling ends, never at every image boundary.
export class HeroFrameMotion {
  position = 0;
  private from = 0;
  private destination = 0;
  private elapsed = 0;
  private duration = 0.18;
  private settling = false;

  get moving() { return this.position !== this.destination; }

  advance(target: number, seconds: number, scrolling = false): number {
    if (scrolling) {
      this.settling = false;
      this.destination = target;
      this.position = approach(this.position, target, seconds, 0.06);
      return this.position;
    }
    const destination = Math.round(target);
    if (!this.settling || destination !== this.destination) {
      this.settling = true;
      this.from = this.position;
      this.destination = destination;
      this.elapsed = 0;
    }
    this.elapsed = Math.min(this.duration, this.elapsed + Math.max(0, seconds));
    if (this.elapsed >= this.duration) this.position = this.destination;
    else {
      const progress = 1 - (1 - this.elapsed / this.duration) ** 3;
      this.position = this.from + (this.destination - this.from) * progress;
    }
    return this.position;
  }

  hold(index: number) {
    this.position = this.from = this.destination = Math.round(index);
    this.elapsed = this.duration;
    this.settling = false;
  }
}
