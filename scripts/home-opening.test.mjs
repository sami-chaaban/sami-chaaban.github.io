import test from 'node:test';
import assert from 'node:assert/strict';
import { homeOpening } from '../src/lib/home-opening.ts';

test('molecular playback begins on the first scroll movement while the title is still visible', () => {
  assert.equal(homeOpening(0).frameProgress, 0);
  const firstMovement = homeOpening(1 / 1000);
  assert.ok(firstMovement.frameProgress > 0);
  assert.ok(firstMovement.titleOpacity > 0);
});

test('research appears while molecular playback and tint expansion continue', () => {
  assert.equal(homeOpening(1).researchOpacity, 1);
  assert.equal(homeOpening(0.32).researchOpacity, 0);
  assert.ok(homeOpening(0.4).researchOpacity > 0);
  assert.equal(homeOpening(0.52).researchOpacity, 1);
  assert.ok(homeOpening(0.52).frameProgress < 1);
  assert.ok(homeOpening(0.52).tintProgress > 0);
  assert.ok(homeOpening(0.52).tintProgress < 1);
});

test('the tint covers the full screen at the last frame, with no late fade', () => {
  const before = homeOpening(0.97);
  const finished = homeOpening(0.98);
  assert.ok(before.frameProgress < 1);
  assert.ok(before.tintProgress < 1);
  assert.equal(finished.frameProgress, 1);
  assert.equal(finished.tintProgress, 1);
  assert.equal(homeOpening(1).frameProgress, 1);
});

test('the left-to-right tint expands gradually throughout playback without clearing again', () => {
  assert.equal(homeOpening(0).tintProgress, 0);
  const positions = [0, 0.001, 0.12, 0.22, 0.38, 0.54, 0.76, 0.97, 0.98];
  const expansion = positions.map(progress => homeOpening(progress).tintProgress);
  assert.ok(expansion.every((value, index) => index === 0 || value > expansion[index - 1]));
  assert.equal(homeOpening(0.49).tintProgress, 0.5);
  assert.equal(homeOpening(1).tintProgress, 1);
});

test('reversals, overscroll and jumps restore the exact opening without persistent state', () => {
  const positions = [0, 0.1, 0.3, 0.5, 0.9, 1];
  const forward = positions.map(progress => homeOpening(progress));
  const reverse = [...positions].reverse().map(progress => homeOpening(progress));
  assert.deepEqual(reverse, [...forward].reverse());
  assert.deepEqual(homeOpening(-0.5), homeOpening(0));
  assert.deepEqual(homeOpening(5), homeOpening(1));
});

test('reduced motion keeps the opening and research text available with a static image', () => {
  for (const progress of [0, 0.5, 1]) {
    const state = homeOpening(progress, true);
    assert.equal(state.titleOpacity, 1);
    assert.equal(state.researchOpacity, 1);
    assert.equal(state.frameProgress, 0);
    assert.equal(state.tintProgress, 1);
  }
});
