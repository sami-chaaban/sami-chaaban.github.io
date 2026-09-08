import test from 'node:test';
import assert from 'node:assert/strict';
import { HeroRenderer } from '../src/lib/hero-renderer.ts';

function surface() {
  let pixel = 0;
  let allocations = 0;
  let width = 300;
  let height = 150;
  const context = {
    globalAlpha: 1,
    drawImage(image) {
      pixel = image.pixel * this.globalAlpha + pixel * (1 - this.globalAlpha);
    },
  };
  const canvas = {
    style: { opacity: '0' },
    get width() { return width; },
    set width(value) { width = value; pixel = 0; allocations++; },
    get height() { return height; },
    set height(value) { height = value; pixel = 0; allocations++; },
    getContext() { return context; },
  };
  return { canvas, context, get pixel() { return pixel; }, get allocations() { return allocations; } };
}

function frame(index, pixel) {
  return { index, image: {
    naturalWidth: 960,
    naturalHeight: 936,
    pixel,
    set src(_) { assert.fail('decoded images must not have their source replaced'); },
  } };
}

test('scrolling in either direction blends complete images without clearing the surface', () => {
  const output = surface();
  const renderer = new HeroRenderer(output.canvas);
  assert.equal(renderer.render(frame(0, 100), frame(1, 200), 0.25), true);
  assert.equal(output.pixel, 125);
  assert.equal(output.canvas.style.opacity, '1');
  const allocations = output.allocations;

  renderer.render(frame(8, 200), frame(9, 100), 0.75);
  assert.equal(output.pixel, 125);
  renderer.render(frame(1, 40), frame(2, 80), 0.5);
  assert.equal(output.pixel, 60);
  assert.equal(output.allocations, allocations);
  assert.equal(output.context.globalAlpha, 1);
});

test('an unavailable image leaves the last complete picture visible', () => {
  const output = surface();
  const renderer = new HeroRenderer(output.canvas);
  renderer.render(frame(0, 100));
  const missing = frame(1, 200);
  missing.image.naturalWidth = 0;
  assert.equal(renderer.render(missing), false);
  assert.equal(renderer.render(frame(2, 50), missing, 0.5), false);
  assert.equal(output.pixel, 100);
  assert.equal(output.canvas.style.opacity, '1');
});

test('slow-loading fallback holds the dominant frame sharply in either blend direction', () => {
  const output = surface();
  const renderer = new HeroRenderer(output.canvas);
  assert.equal(renderer.hold(), undefined);
  renderer.render(frame(3, 80), frame(4, 160), 0.8);
  assert.equal(renderer.hold(), 4);
  assert.equal(output.pixel, 160);
  renderer.render(frame(2, 40), frame(3, 80), 0.2);
  assert.equal(renderer.hold(), 2);
  assert.equal(output.pixel, 40);
  renderer.reset();
  assert.equal(output.canvas.style.opacity, '0');
  assert.equal(renderer.hold(), undefined);
});

test('canvas-unavailable browsers keep the responsive poster', () => {
  const canvas = { style: { opacity: '0' }, getContext() { return null; } };
  const renderer = new HeroRenderer(canvas);
  assert.equal(renderer.available, false);
  assert.equal(renderer.render(frame(0, 100)), false);
  assert.equal(canvas.style.opacity, '0');
});
