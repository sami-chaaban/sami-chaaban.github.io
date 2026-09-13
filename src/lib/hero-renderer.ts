export type HeroPaintFrame = { index: number; image: HTMLImageElement };

// Keep one painted surface alive instead of changing the src of visible images.
export class HeroRenderer {
  private canvas: HTMLCanvasElement;
  private context: CanvasRenderingContext2D | null;
  private mirror?: HTMLCanvasElement | null;
  private mirrorContext: CanvasRenderingContext2D | null;
  private base?: HeroPaintFrame;
  private upper?: HeroPaintFrame;
  private blend = 0;

  constructor(canvas: HTMLCanvasElement, mirror?: HTMLCanvasElement | null) {
    this.canvas = canvas;
    this.context = canvas.getContext('2d', { alpha: false });
    this.mirror = mirror;
    this.mirrorContext = mirror?.getContext('2d', { alpha: false }) ?? null;
  }

  get available() {
    return this.context !== null && (!this.mirror || this.mirrorContext !== null);
  }

  render(base: HeroPaintFrame, upper?: HeroPaintFrame, blend = 0): boolean {
    const context = this.context;
    const { naturalWidth: width, naturalHeight: height } = base.image;
    if (!this.available || !context || !width || !height) return false;
    if (upper && (!upper.image.naturalWidth || !upper.image.naturalHeight)) return false;

    // Resizing clears a canvas: only allocate when the source dimensions change,
    // and fill it in the same task before exposing it to the browser.
    if (this.canvas.width !== width) this.canvas.width = width;
    if (this.canvas.height !== height) this.canvas.height = height;
    const opacity = Math.min(1, Math.max(0, blend));
    context.globalAlpha = 1;
    context.drawImage(base.image, 0, 0, width, height);
    if (upper && opacity > 0) {
      context.globalAlpha = opacity;
      context.drawImage(upper.image, 0, 0, width, height);
      context.globalAlpha = 1;
    }
    // Copy the completed blend in the same paint. The masked correction must
    // show exactly this frame, including sharp holds, without another loader.
    if (this.mirror && this.mirrorContext) {
      if (this.mirror.width !== width) this.mirror.width = width;
      if (this.mirror.height !== height) this.mirror.height = height;
      this.mirrorContext.drawImage(this.canvas, 0, 0, width, height);
      this.mirror.style.opacity = '1';
    }
    this.base = base;
    this.upper = upper;
    this.blend = opacity;
    this.canvas.style.opacity = '1';
    return true;
  }

  hold(): number | undefined {
    const frame = this.upper && this.blend >= 0.5 ? this.upper : this.base;
    if (frame) this.render(frame);
    return frame?.index;
  }

  reset() {
    this.canvas.style.opacity = '0';
    if (this.mirror) this.mirror.style.opacity = '0';
    this.base = undefined;
    this.upper = undefined;
    this.blend = 0;
  }
}
