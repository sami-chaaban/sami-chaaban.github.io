import test from 'node:test';
import assert from 'node:assert/strict';
import { approach, smoothstep } from '../src/lib/motion-math.ts';
import { cellPath } from '../src/lib/home-cell.ts';

test('a half-second of motion settles equally across refresh rates', () => {
  const settle = fps => {
    let value = 0;
    for (let frame = 0; frame < fps / 2; frame++) value = approach(value, 1, 1 / fps);
    return value;
  };
  assert.ok(Math.abs(settle(60) - settle(120)) < 1e-12);
  assert.ok(Math.abs(settle(30) - settle(144)) < 1e-12);
});

test('reversing scroll stays bounded and responds without overshoot', () => {
  let position = 0;
  for (let i = 0; i < 10; i++) position = approach(position, 1, 1 / 60);
  const turningPoint = position;
  for (let i = 0; i < 30; i++) {
    position = approach(position, 0, 1 / 60);
    assert.ok(position >= 0 && position < turningPoint);
  }
  assert.ok(position < 0.04);
  assert.equal(approach(position, 1, 0), position);
  assert.equal(approach(position, 1, -1), position);
});

test('overscroll clamps fades and produces finite cell outlines at both ends', () => {
  assert.equal(smoothstep(0.1, 1.5, -1), 0);
  assert.equal(smoothstep(0.1, 1.5, 10), 1);
  assert.equal(cellPath(-1), cellPath(0));
  assert.equal(cellPath(2), cellPath(1));
  for (const progress of [0, 0.48, 0.72, 1]) {
    assert.doesNotMatch(cellPath(progress), /NaN|Infinity/);
  }
  assert.equal((cellPath(0).match(/M /g) ?? []).length, 1);
  assert.equal((cellPath(1).match(/M /g) ?? []).length, 2);
});
