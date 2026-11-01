import { warpMatrix } from './hero-depth.ts';
import type { DepthCamera, DepthView } from './hero-depth.ts';

export type DepthPaintFrame = { angle: number; color: HTMLImageElement; depth: HTMLImageElement };

const WARP_VS = `#version 300 es
precision highp float;
precision highp int;
uniform highp sampler2D uDepth;
uniform mat4 uWarp;
uniform vec4 uProjection;
uniform vec4 uView;
uniform vec2 uSize;
uniform float uInvMax;
uniform int uMargin;
out vec2 vSource;
void main() {
  // One triangle strip per pixel row, overhanging the frame so moving edges stay covered.
  ivec2 cell = ivec2(gl_VertexID >> 1, gl_InstanceID + (gl_VertexID & 1)) - uMargin;
  float inverseDepth = texelFetch(uDepth, clamp(cell, ivec2(0), ivec2(uSize) - 1), 0).r * uInvMax;
  vSource = vec2(cell) + 0.5;
  vec2 full = vSource + uView.xy;
  vec2 ndc = vec2(full.x / uView.z * 2.0 - 1.0, 1.0 - full.y / uView.w * 2.0);
  // Homogeneous camera-space point with w = 1/depth, so empty background stays finite.
  vec4 clip = uWarp * vec4((ndc.x + uProjection.z) / uProjection.x, (ndc.y + uProjection.w) / uProjection.y, -1.0, inverseDepth);
  // Depth-test on inverse depth: linear in screen space and precise far away.
  float nearness = clamp(inverseDepth / (max(clip.w, 1e-6) * uInvMax), 0.0, 1.0);
  gl_Position = vec4(clip.xy, (1.0 - 2.0 * nearness) * clip.w, clip.w);
}`;

const WARP_FS = `#version 300 es
precision highp float;
uniform sampler2D uColor;
uniform vec2 uSize;
in vec2 vSource;
out vec4 outColor;
void main() {
  // Source pixels per screen pixel: about 1 on surfaces, near 0 on triangles that
  // bridge a depth edge. Those stretched sheets are only trusted as a last resort.
  vec2 dx = dFdx(vSource), dy = dFdy(vSource);
  float area = abs(dx.x * dy.y - dx.y * dy.x);
  outColor = vec4(texture(uColor, vSource / uSize).rgb, smoothstep(0.08, 0.6, area));
}`;

const BLEND_VS = `#version 300 es
out vec2 vUv;
void main() {
  vUv = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
  gl_Position = vec4(vUv * 2.0 - 1.0, 0.0, 1.0);
}`;

const BLEND_FS = `#version 300 es
precision highp float;
uniform sampler2D uWarpA, uWarpB;
uniform float uMix;
in vec2 vUv;
out vec4 outColor;
void main() {
  vec4 a = texture(uWarpA, vUv), b = texture(uWarpB, vUv);
  float wa = (1.0 - uMix) * a.a, wb = uMix * b.a;
  outColor = vec4(wa + wb > 1e-4 ? (a.rgb * wa + b.rgb * wb) / (wa + wb) : mix(a.rgb, b.rgb, uMix), 1.0);
}`;

const MARGIN = 32;

type Program = { program: WebGLProgram; uniforms: Record<string, WebGLUniformLocation | null> };
type Target = { texture: WebGLTexture; depth: WebGLRenderbuffer; framebuffer: WebGLFramebuffer; width: number; height: number };

function compile(gl: WebGL2RenderingContext, vertex: string, fragment: string): Program {
  const program = gl.createProgram()!;
  for (const [type, source] of [[gl.VERTEX_SHADER, vertex], [gl.FRAGMENT_SHADER, fragment]] as const) {
    const shader = gl.createShader(type)!;
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(shader) ?? 'shader');
    gl.attachShader(program, shader);
  }
  gl.linkProgram(program);
  if (!gl.getProgramParameter(program, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(program) ?? 'program');
  const uniforms: Program['uniforms'] = {};
  for (let i = 0; i < gl.getProgramParameter(program, gl.ACTIVE_UNIFORMS); i++) {
    const { name } = gl.getActiveUniform(program, i)!;
    uniforms[name] = gl.getUniformLocation(program, name);
  }
  return { program, uniforms };
}

// Re-projects Cycles keyframes through their Z pass into the camera for any angle on the
// orbit, so playback moves with true parallax instead of stepping between frames.
export class HeroDepthRenderer {
  private gl: WebGL2RenderingContext;
  private canvas: HTMLCanvasElement;
  private mirror?: HTMLCanvasElement | null;
  private mirrorContext: CanvasRenderingContext2D | null;
  private camera: DepthCamera;
  private view: DepthView;
  private warp: Program;
  private blend: Program;
  private vao: WebGLVertexArrayObject;
  private targets: Target[] = [];
  private textures = new Map<HTMLImageElement, WebGLTexture>();
  private lost = false;

  // Null when WebGL2 is unavailable, leaving the canvas free for a 2D context. A context
  // that exists but cannot compile keeps the canvas, so the page falls back to the poster.
  static create(canvas: HTMLCanvasElement, mirror: HTMLCanvasElement | null | undefined,
    camera: DepthCamera, view: DepthView, onLost: () => void): HeroDepthRenderer | null {
    const gl = canvas.getContext('webgl2', { alpha: false, antialias: false, depth: false, stencil: false });
    if (!gl) return null;
    try {
      return new HeroDepthRenderer(canvas, gl, mirror, camera, view, onLost);
    } catch {
      return null;
    }
  }

  private constructor(canvas: HTMLCanvasElement, gl: WebGL2RenderingContext, mirror: HTMLCanvasElement | null | undefined,
    camera: DepthCamera, view: DepthView, onLost: () => void) {
    this.canvas = canvas;
    this.gl = gl;
    this.mirror = mirror;
    this.mirrorContext = mirror?.getContext('2d', { alpha: false }) ?? null;
    this.camera = camera;
    this.view = view;
    this.warp = compile(gl, WARP_VS, WARP_FS);
    this.blend = compile(gl, BLEND_VS, BLEND_FS);
    this.vao = gl.createVertexArray()!;
    gl.pixelStorei(gl.UNPACK_COLORSPACE_CONVERSION_WEBGL, gl.NONE);
    canvas.addEventListener('webglcontextlost', () => {
      this.lost = true;
      this.textures.clear();
      this.reset();
      onLost();
    });
  }

  get available() {
    return !this.lost && (!this.mirror || this.mirrorContext !== null);
  }

  // Applied at the next render, so the last picture stays until new frames can replace it.
  setView(view: DepthView) {
    this.view = view;
    this.textures.forEach(texture => this.gl.deleteTexture(texture));
    this.textures.clear();
  }

  release(image: HTMLImageElement) {
    const texture = this.textures.get(image);
    if (texture) this.gl.deleteTexture(texture);
    this.textures.delete(image);
  }

  render(a: DepthPaintFrame, b: DepthPaintFrame, mix: number, angle: number): boolean {
    const gl = this.gl;
    if (!this.available || gl.isContextLost()) return false;
    const images = [a.color, a.depth, b.color, b.depth];
    if (images.some(image => !image.naturalWidth || !image.naturalHeight)) return false;
    const { width, height } = this.view;
    // Resizing clears the canvas: only allocate when the view changes, then fill it in this task.
    if (this.canvas.width !== width) this.canvas.width = width;
    if (this.canvas.height !== height) this.canvas.height = height;
    this.prepareTargets(width, height);
    this.draw(this.targets[0], a, angle);
    if (b.color !== a.color) this.draw(this.targets[1], b, angle);

    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, width, height);
    const { program, uniforms } = this.blend;
    gl.useProgram(program);
    this.bind(0, uniforms.uWarpA, this.targets[0].texture);
    this.bind(1, uniforms.uWarpB, this.targets[b.color !== a.color ? 1 : 0].texture);
    gl.uniform1f(uniforms.uMix, b.color !== a.color ? Math.min(1, Math.max(0, mix)) : 0);
    gl.bindVertexArray(this.vao);
    gl.drawArrays(gl.TRIANGLES, 0, 3);

    // Copy in the same task, before the WebGL buffer is presented and cleared.
    if (this.mirror && this.mirrorContext) {
      if (this.mirror.width !== width) this.mirror.width = width;
      if (this.mirror.height !== height) this.mirror.height = height;
      this.mirrorContext.drawImage(this.canvas, 0, 0, width, height);
      this.mirror.style.opacity = '1';
    }
    this.canvas.style.opacity = '1';
    return true;
  }

  // The presented picture persists between renders.
  hold(): number | undefined {
    return undefined;
  }

  reset() {
    this.canvas.style.opacity = '0';
    if (this.mirror) this.mirror.style.opacity = '0';
  }

  private draw(target: Target, frame: DepthPaintFrame, angle: number) {
    const gl = this.gl;
    const { width, height } = this.view;
    const { program, uniforms } = this.warp;
    gl.bindFramebuffer(gl.FRAMEBUFFER, target.framebuffer);
    gl.viewport(0, 0, width, height);
    gl.clearColor(0, 0, 0, 0);
    gl.clearDepth(1);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.enable(gl.DEPTH_TEST);
    // Empty background sits exactly on the far plane, which equals the cleared depth.
    gl.depthFunc(gl.LEQUAL);
    gl.useProgram(program);
    const p = this.camera.projection;
    gl.uniformMatrix4fv(uniforms.uWarp, false, warpMatrix(this.camera, this.view, frame.angle, angle));
    gl.uniform4f(uniforms.uProjection, p[0], p[5], p[8], p[9]);
    gl.uniform4f(uniforms.uView, this.view.left, this.view.top, this.camera.width, this.camera.height);
    gl.uniform2f(uniforms.uSize, width, height);
    gl.uniform1f(uniforms.uInvMax, this.camera.invMax);
    gl.uniform1i(uniforms.uMargin, MARGIN);
    this.bind(0, uniforms.uDepth, this.texture(frame.depth, true));
    this.bind(1, uniforms.uColor, this.texture(frame.color, false));
    gl.bindVertexArray(this.vao);
    gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, (width + 2 * MARGIN) * 2, height + 2 * MARGIN - 1);
    gl.disable(gl.DEPTH_TEST);
  }

  private bind(unit: number, location: WebGLUniformLocation | null, texture: WebGLTexture) {
    const gl = this.gl;
    gl.activeTexture(gl.TEXTURE0 + unit);
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.uniform1i(location, unit);
  }

  private texture(image: HTMLImageElement, depth: boolean): WebGLTexture {
    const gl = this.gl;
    let texture = this.textures.get(image);
    if (texture) return texture;
    texture = gl.createTexture()!;
    gl.bindTexture(gl.TEXTURE_2D, texture);
    if (depth) gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, gl.RED, gl.UNSIGNED_BYTE, image);
    else gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, gl.RGBA, gl.UNSIGNED_BYTE, image);
    // Depth values must never be filtered across an edge.
    const filter = depth ? gl.NEAREST : gl.LINEAR;
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, filter);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, filter);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    this.textures.set(image, texture);
    return texture;
  }

  private prepareTargets(width: number, height: number) {
    const gl = this.gl;
    if (this.targets.length && this.targets[0].width === width && this.targets[0].height === height) return;
    this.targets.forEach(({ texture, depth, framebuffer }) => {
      gl.deleteTexture(texture);
      gl.deleteRenderbuffer(depth);
      gl.deleteFramebuffer(framebuffer);
    });
    this.targets = [0, 1].map(() => {
      const texture = gl.createTexture()!;
      gl.bindTexture(gl.TEXTURE_2D, texture);
      gl.texStorage2D(gl.TEXTURE_2D, 1, gl.RGBA8, width, height);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      const depth = gl.createRenderbuffer()!;
      gl.bindRenderbuffer(gl.RENDERBUFFER, depth);
      gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT24, width, height);
      const framebuffer = gl.createFramebuffer()!;
      gl.bindFramebuffer(gl.FRAMEBUFFER, framebuffer);
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, texture, 0);
      gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, depth);
      return { texture, depth, framebuffer, width, height };
    });
  }
}
