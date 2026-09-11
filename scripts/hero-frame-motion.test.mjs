import test from 'node:test';
import assert from 'node:assert/strict';
import { HeroFrameMotion } from '../src/lib/hero-frame-motion.ts';

test('scrolling follows fractional targets across image boundaries without rounding', () => {
  const motion = new HeroFrameMotion();
  const before = motion.advance(0.49, 1 / 60, true);
  const after = motion.advance(0.51, 1 / 60, true);
  assert.ok(before > 0 && before < 0.49);
  assert.ok(after > before && after < 0.51);
  assert.ok(after - before < 0.15);
});

test('scrolling down and back to the top repeatedly always restores frame zero', () => {
  const motion = new HeroFrameMotion();
  for (let cycle = 0; cycle < 5; cycle++) {
    for (let frame = 0; frame <= 18; frame += 0.25) motion.advance(frame, 1 / 60, true);
    for (let frame = 18; frame >= 0; frame -= 0.25) motion.advance(frame, 1 / 60, true);
    motion.advance(0, 0.18, false);
    assert.equal(motion.position, 0);
    assert.equal(motion.moving, false);
  }
});

test('fractional scroll positions finish on a whole frame within 180 ms', () => {
  for (const target of [0.25, 0.5, 4.33, 10.8, 18]) {
    const motion = new HeroFrameMotion();
    motion.advance(target, 0.09);
    motion.advance(target, 0.09);
    assert.equal(motion.position, Math.round(target));
    assert.equal(motion.moving, false);
  }
});

test('a paused animation settles sharply even while page scrolling continues', () => {
  for (const fps of [30, 60, 120]) {
    const motion = new HeroFrameMotion();
    motion.hold(13);
    for (let i = 0; i < Math.ceil(fps * 0.3); i++) motion.advance(14.5, 1 / fps, true);
    assert.equal(motion.position, 15);
    assert.equal(motion.moving, false);
    for (let i = 0; i < fps * 3; i++) {
      assert.equal(motion.advance(14.5, 1 / fps, true), 15);
    }
    assert.ok(motion.advance(15.5, 1 / fps, true) > 15);
  }
});

test('late momentum cannot interrupt sharpening or change the chosen frame at a rounding boundary', () => {
  for (const direction of [-1, 1]) {
    const motion = new HeroFrameMotion();
    const start = direction > 0 ? 4 : 5;
    motion.hold(start);
    motion.advance(4.5 - direction * 0.01, 0.1, true);
    let previous = motion.advance(4.5 - direction * 0.01, 0.09, false);
    for (let i = 0; i < 60; i++) {
      const target = 4.5 + (i % 2 ? -0.01 : 0.01);
      const position = motion.advance(target, 1 / 60, true);
      assert.ok(Math.abs(position - start) <= Math.abs(previous - start));
      previous = position;
    }
    assert.equal(motion.position, start);
    assert.equal(motion.moving, false);
  }
});

test('a finished sharp frame stays locked through residual drift and resumes smoothly in either direction', () => {
  for (const nextTarget of [3.9, 5.1]) {
    const motion = new HeroFrameMotion();
    motion.hold(4);
    motion.advance(4.49, 0.18, false);
    for (const target of [4.51, 4.48, 4.53]) {
      assert.equal(motion.advance(target, 1 / 60, true), 4);
      assert.equal(motion.moving, false);
    }
    assert.equal(motion.advance(nextTarget, 0, true), 4);
    const resumed = motion.advance(nextTarget, 1 / 60, true);
    assert.ok(resumed > Math.min(4, nextTarget) && resumed < Math.max(4, nextTarget));
    motion.advance(nextTarget, 0.18, false);
    assert.equal(motion.position, Math.round(nextTarget));
  }
});

test('transition timing agrees across display refresh rates', () => {
  const at100ms = fps => {
    const motion = new HeroFrameMotion();
    for (let i = 0; i < fps / 10; i++) motion.advance(12, 1 / fps);
    return motion.position;
  };
  assert.ok(Math.abs(at100ms(60) - at100ms(120)) < 1e-12);
  assert.ok(at100ms(60) > 0 && at100ms(60) < 12);
});

test('rapid reversal continues from the visible position and settles sharply', () => {
  const motion = new HeroFrameMotion();
  const turningPoint = motion.advance(18, 0.06);
  assert.equal(motion.advance(0, 0), turningPoint);
  assert.ok(motion.advance(0, 0.06) < turningPoint);
  motion.advance(0, 0.12);
  assert.equal(motion.position, 0);
  assert.equal(motion.moving, false);
});

test('a stalled image holds a complete frame, then resumes without a jump', () => {
  const motion = new HeroFrameMotion();
  motion.advance(10, 0.06);
  motion.hold(3);
  assert.equal(motion.position, 3);
  assert.equal(motion.moving, false);
  assert.equal(motion.advance(10, 0), 3);
  motion.advance(10, 0.18);
  assert.equal(motion.position, 10);
});
