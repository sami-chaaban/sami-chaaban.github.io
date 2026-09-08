export type ConnectionHints = {
  saveData?: boolean;
  effectiveType?: string;
  downlink?: number;
};

export function useStaticHero(reducedMotion: boolean, connection?: ConnectionHints): boolean {
  return reducedMotion || Boolean(connection?.saveData)
    || ['slow-2g', '2g'].includes(connection?.effectiveType ?? '');
}

export function framePriority(position: number, count: number): number[] {
  return Array.from({ length: count }, (_, index) => index)
    .sort((a, b) => Math.abs(a - position) - Math.abs(b - position));
}

type LoaderOptions<T> = {
  sources: string[];
  load: (source: string, signal: AbortSignal) => Promise<T>;
  onReady: (index: number) => void;
  onFallback: () => void;
  release?: (frame: T) => void;
  concurrency?: number;
  timeoutMs?: number;
};

// A bounded queue that follows the scroll position and never exposes undecoded frames.
export class HeroLoader<T> {
  readonly ready = new Map<number, T>();
  private queue: number[] = [];
  private pending = new Map<number, AbortController>();
  private failed = new Set<number>();
  private paused = true;
  private stopped = false;
  private disposed = false;
  private options: LoaderOptions<T>;

  constructor(options: LoaderOptions<T>) {
    this.options = options;
  }

  request(priority: number[]) {
    if (this.stopped) return;
    this.queue = [...new Set(priority)].filter((index) => index >= 0 && index < this.options.sources.length);
    this.paused = false;
    this.pump();
  }

  pause() {
    this.paused = true;
    this.queue = [];
    this.pending.forEach((controller) => controller.abort());
  }

  stop() {
    if (this.stopped) return;
    this.stopped = true;
    this.pause();
    this.options.onFallback();
  }

  dispose() {
    this.disposed = true;
    this.stopped = true;
    this.pause();
    this.ready.forEach((frame) => this.options.release?.(frame));
    this.ready.clear();
  }

  private pump() {
    if (this.paused || this.stopped) return;
    while (this.pending.size < (this.options.concurrency ?? 2)) {
      const index = this.queue.find((candidate) => !this.ready.has(candidate) && !this.pending.has(candidate) && !this.failed.has(candidate));
      if (index === undefined) break;
      const controller = new AbortController();
      this.pending.set(index, controller);
      const timeout = setTimeout(() => {
        if (!controller.signal.aborted) this.stop();
      }, this.options.timeoutMs ?? 1800);

      Promise.resolve().then(() => this.options.load(this.options.sources[index], controller.signal)).then((frame) => {
        if (controller.signal.aborted || this.disposed) {
          this.options.release?.(frame);
          return;
        }
        this.ready.set(index, frame);
        this.options.onReady(index);
      }).catch(() => {
        if (controller.signal.aborted) return;
        this.failed.add(index);
        if (this.failed.size >= 2) this.stop();
      }).finally(() => {
        clearTimeout(timeout);
        this.pending.delete(index);
        this.pump();
      });
    }
  }
}

export type DecodedFrame = { url: string; image: HTMLImageElement; owned: boolean };

export async function loadDecodedFrame(source: string, signal: AbortSignal): Promise<DecodedFrame> {
  const response = await fetch(source, { signal, cache: 'force-cache' });
  if (!response.ok) throw new Error(`Hero image returned ${response.status}`);
  const url = URL.createObjectURL(await response.blob());
  try {
    const image = new Image();
    image.decoding = 'async';
    image.src = url;
    await image.decode();
    signal.throwIfAborted();
    return { url, image, owned: true };
  } catch (error) {
    URL.revokeObjectURL(url);
    throw error;
  }
}
