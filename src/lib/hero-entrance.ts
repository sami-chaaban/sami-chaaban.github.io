// The opening occupies the first part of the scroll timeline in both directions.
export function heroEntranceOpacity(progress: number, reducedMotion = false): number {
  if (reducedMotion) return progress > 0 ? 0 : 0.9;
  const amount = Math.min(1, Math.max(0, progress / 0.1));
  return 0.9 * (1 - amount * amount * (3 - 2 * amount));
}

// Either the first gesture or the opening image can arrive first.
export class HeroEntrance {
  private ready = false;
  private started = false;
  private pendingProgress?: number;
  private onStart: (progress: number) => void;

  constructor(onStart: (progress: number) => void) {
    this.onStart = onStart;
  }

  get hasStarted() {
    return this.started;
  }

  request(progress: number): boolean {
    if (this.started) return false;
    // Do not skip the opening if someone keeps scrolling while it loads.
    this.pendingProgress ??= progress;
    this.startIfReady();
    return true;
  }

  markReady() {
    this.ready = true;
    this.startIfReady();
  }

  private startIfReady() {
    if (!this.ready || this.started || this.pendingProgress === undefined) return;
    this.started = true;
    this.onStart(this.pendingProgress);
  }
}
