import { approach } from './motion-math.ts';

// Follow continuous playback; settle when its target stops, even if the page keeps scrolling.
export class HeroFrameMotion {
  position = 0;
  private from = 0;
  private destination = 0;
  private elapsed = 0;
  private duration = 0.18;
  private settling = false;
  private settleTarget = 0;
  private previousTarget = NaN;
  private stationaryFor = 0;

  get moving() { return this.position !== this.destination; }

  advance(target: number, seconds: number, scrolling = false): number {
    this.stationaryFor = target === this.previousTarget ? this.stationaryFor + Math.max(0, seconds) : 0;
    this.previousTarget = target;
    // Once sharpening starts, tiny momentum updates must not reopen the blend
    // or switch its winning frame. Half a source-frame of travel releases it;
    // measuring from the captured target avoids jitter at rounding boundaries.
    const keepChoice = this.settling && Math.abs(target - this.settleTarget) < 0.5;
    // Use the same 100 ms quiet window for a paused animation as for a stopped page.
    if (scrolling && this.stationaryFor < 0.1 && !keepChoice) {
      this.settling = false;
      this.destination = target;
      this.position = approach(this.position, target, seconds, 0.06);
      return this.position;
    }
    const destination = keepChoice ? this.destination : Math.round(target);
    if (!this.settling || destination !== this.destination) {
      this.settling = true;
      this.from = this.position;
      this.destination = destination;
      this.settleTarget = target;
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
