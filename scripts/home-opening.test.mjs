import test from 'node:test';
import assert from 'node:assert/strict';
import { homeOpening, homeBackground, homeFramePosition } from '../src/lib/home-opening.ts';

test('molecular playback begins on the first scroll movement while the title is still visible', () => {
  assert.equal(homeBackground(0).frameProgress, 0);
  const firstMovement = homeBackground(1 / 1000);
  assert.ok(firstMovement.frameProgress > 0);
  assert.ok(homeOpening(1 / 1000).titleOpacity > 0);
});

test('research appears while molecular playback and tint expansion continue', () => {
  assert.equal(homeOpening(1).researchOpacity, 1);
  assert.equal(homeOpening(0.26).researchOpacity, 0);
  assert.ok(homeOpening(0.32).researchOpacity > 0);
  assert.equal(homeOpening(0.46).researchOpacity, 1);
  assert.equal(homeBackground(0.46).frameProgress, 0.46);
  assert.ok(homeBackground(0.46).tintProgress > 0);
  assert.ok(homeBackground(0.46).tintProgress < 1);
});

test('the opening plays the entire sequence and later sections hold the final frame', () => {
  assert.equal(homeBackground(0.5).frameProgress, 0.5);
  assert.equal(homeFramePosition(homeBackground(0.5).frameProgress, 37), 18);
  for (const progress of [1, 1.25, 2, 3, 4, 5, 10]) {
    assert.equal(homeBackground(progress).frameProgress, 1);
    assert.equal(homeFramePosition(homeBackground(progress).frameProgress, 37), 36);
  }
});

test('playback reaches the exact final image and stays within the available sequence', () => {
  for (const count of [0, 1, 29, 30, 31, 35, 37]) {
    const last = Math.max(0, count - 1);
    assert.equal(homeFramePosition(-1, count), 0);
    assert.equal(homeFramePosition(0, count), 0);
    assert.equal(homeFramePosition(1, count), last);
    assert.equal(homeFramePosition(5, count), last);
    const beforeArrival = homeFramePosition(0.9999, count);
    assert.ok(beforeArrival <= last);
    assert.ok(last - beforeArrival < 0.01);
  }
});

test('the correction and overlay finish with the animation at research arrival', () => {
  const before = homeBackground(0.99);
  const finished = homeBackground(1);
  assert.ok(before.frameProgress < 1);
  assert.ok(before.tintProgress < 1);
  assert.equal(finished.frameProgress, 1);
  assert.equal(finished.tintProgress, 1);
  for (const progress of [1.1, 2, 3, 5]) {
    assert.deepEqual(homeBackground(progress), finished);
  }
});

test('the tint expands throughout the opening and stays complete in later sections', () => {
  assert.equal(homeBackground(0).tintProgress, 0);
  const expansion = [0, 0.25, 0.5, 0.75, 1].map(progress => homeBackground(progress).tintProgress);
  assert.ok(expansion.every((value, index) => index === 0 || value > expansion[index - 1]));
  assert.equal(homeBackground(0.5).tintProgress, 0.5);
  for (const progress of [1, 1.5, 2, 3, 4, 5]) {
    assert.equal(homeBackground(progress).tintProgress, 1);
  }
});

test('reversals, overscroll and jumps restore the exact opening without persistent state', () => {
  const positions = [0, 0.1, 0.3, 0.5, 0.9, 1, 2, 3, 5];
  const forward = positions.map(progress => homeOpening(progress));
  const reverse = [...positions].reverse().map(progress => homeOpening(progress));
  assert.deepEqual(reverse, [...forward].reverse());
  assert.deepEqual(homeOpening(-0.5), homeOpening(0));
  assert.deepEqual(homeOpening(5), homeOpening(1));
  const playback = positions.map(progress => homeBackground(progress));
  const rewind = [...positions].reverse().map(progress => homeBackground(progress));
  assert.deepEqual(rewind, [...playback].reverse());
  assert.deepEqual(homeBackground(-1), homeBackground(0));
});

test('reduced motion keeps the opening and research text available with a static image', () => {
  for (const progress of [0, 0.5, 1, 5]) {
    const state = homeOpening(progress, true);
    assert.equal(state.titleOpacity, 1);
    assert.equal(state.researchOpacity, 1);
    const background = homeBackground(progress, true);
    assert.equal(background.frameProgress, 0);
    assert.equal(background.tintProgress, 1);
  }
});
