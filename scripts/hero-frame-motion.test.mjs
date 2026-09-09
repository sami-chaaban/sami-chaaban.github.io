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
