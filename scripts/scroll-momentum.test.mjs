import test from 'node:test';
import assert from 'node:assert/strict';
import { ScrollMomentum } from '../src/lib/scroll-momentum.ts';

test('a wheel gesture continues after release with decreasing travel, then stops at the requested distance', () => {
  const motion = new ScrollMomentum();
  motion.push(240, 2000);
  let previous = 0;
  let lastStep = Infinity;
  for (let frame = 0; frame < 120; frame++) {
    const position = motion.advance(1 / 60, 2000);
    const step = position - previous;
    assert.ok(position >= previous && position <= 240);
    if (motion.moving) assert.ok(step <= lastStep);
    previous = position;
    lastStep = step;
  }
  assert.equal(motion.position, 240);
  assert.equal(motion.moving, false);
});

test('reversing direction cancels outstanding travel immediately', () => {
  const motion = new ScrollMomentum();
  motion.stop(500);
  motion.push(400, 2000);
  const turningPoint = motion.advance(0.07, 2000);
  motion.push(-80, 2000);
  assert.ok(motion.advance(1 / 60, 2000) < turningPoint);
  for (let i = 0; i < 120; i++) motion.advance(1 / 60, 2000);
  assert.equal(motion.position, turningPoint - 80);
});

test('edges, a shrinking page, and explicit interruption cannot leave residual momentum', () => {
  const motion = new ScrollMomentum();
  motion.push(-200, 1000);
  assert.equal(motion.moving, false);
  motion.stop(990);
  motion.push(200, 1000);
  for (let i = 0; i < 120; i++) motion.advance(1 / 60, 1000);
  assert.equal(motion.position, 1000);
  assert.equal(motion.advance(1 / 60, 600), 600);
  motion.push(-400, 600);
  motion.stop(450);
  assert.equal(motion.advance(1, 600), 450);
  assert.equal(motion.moving, false);
});

test('the glide has the same timing on 60 Hz and 120 Hz displays', () => {
  const travel = fps => {
    const motion = new ScrollMomentum();
    motion.push(300, 2000);
    for (let i = 0; i < fps / 2; i++) motion.advance(1 / fps, 2000);
    return motion.position;
  };
  assert.ok(Math.abs(travel(60) - travel(120)) < 1e-9);
});
