import test from 'node:test';
import assert from 'node:assert/strict';
import { remapHomeScroll } from '../src/lib/home-scroll-layout.ts';

const portrait = [0, 675, 800, 1222, 2066, 2910, 3630, 3630];
const landscape = [0, 312, 312, 507, 897, 1287, 1677, 1710];

test('rotation keeps the same animation frame and pinned reading progress', () => {
  for (const index of [1, 3, 4, 5, 6]) {
    for (const progress of [0.1, 0.35, 0.7, 0.95]) {
      const original = portrait[index - 1] + progress * (portrait[index] - portrait[index - 1]);
      const rotated = remapHomeScroll(original, portrait, landscape);
      const actual = (rotated - landscape[index - 1]) / (landscape[index] - landscape[index - 1]);
      assert.ok(Math.abs(actual - progress) < 1e-12);
      assert.ok(Math.abs(remapHomeScroll(rotated, landscape, portrait) - original) < 1e-9);
    }
  }
});

test('rotation preserves top, bottom, and progress into the final profile', () => {
  assert.equal(remapHomeScroll(0, portrait, landscape), 0);
  assert.equal(remapHomeScroll(-20, portrait, landscape), 0);
  assert.equal(remapHomeScroll(3629.5, portrait, landscape), 1710);
  assert.equal(remapHomeScroll(4000, portrait, landscape), 1710);
  assert.equal(remapHomeScroll((2910 + 3630) / 2, portrait, landscape), (1287 + 1677) / 2);
});

test('the expanded arrival section has finite positions when rotating back from a compact layout', () => {
  for (const position of [0, 311, 312, 313, 507, 800, 897, 1287, 1710]) {
    const rotated = remapHomeScroll(position, landscape, portrait);
    assert.ok(Number.isFinite(rotated));
    assert.ok(rotated >= 0 && rotated <= 3630);
  }
});
