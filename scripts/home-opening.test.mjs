import test from 'node:test';
import assert from 'node:assert/strict';
import { homeOpening, homeBackground, homeFramePosition } from '../src/lib/home-opening.ts';

test('molecular playback begins on the first scroll movement while the title is still visible', () => {
  assert.equal(homeBackground(0, 0, 0).frameProgress, 0);
  const firstMovement = homeBackground(1 / 1000, 0, 0);
  assert.ok(firstMovement.frameProgress > 0);
  assert.ok(homeOpening(1 / 1000).titleOpacity > 0);
});

test('research appears while molecular playback and tint expansion continue', () => {
  assert.equal(homeOpening(1).researchOpacity, 1);
  assert.equal(homeOpening(0.26).researchOpacity, 0);
  assert.ok(homeOpening(0.32).researchOpacity > 0);
  assert.equal(homeOpening(0.46).researchOpacity, 1);
  assert.ok(homeBackground(0.46, 0, 0).frameProgress < 1 / 3);
  assert.ok(homeBackground(0.46, 0, 0).tintProgress > 0);
  assert.ok(homeBackground(0.46, 0, 0).tintProgress < 1);
});

test('each transition plays one third with holds during research and cell division', () => {
  assert.equal(homeBackground(0.5, -1, -2).frameProgress, 1 / 6);
  assert.equal(homeBackground(1, 0.5, -1).frameProgress, 0.5);
  assert.equal(homeBackground(1, 1, 0.5).frameProgress, 5 / 6);
  for (const openingProgress of [1, 1.25, 2]) {
    for (const cellArrivalProgress of [-2, -0.5, 0]) {
      assert.equal(homeBackground(openingProgress, cellArrivalProgress, -1).frameProgress, 1 / 3);
    }
    for (const outroArrivalProgress of [-2, -0.5, 0]) {
      assert.equal(homeBackground(openingProgress, 2, outroArrivalProgress).frameProgress, 2 / 3);
    }
  }
});

test('both reading holds target real frames and all transitions join without reversing', () => {
  for (const count of [0, 1, 29, 30, 31]) {
    const last = Math.max(0, count - 1);
    assert.equal(homeFramePosition(0, count), 0);
    assert.equal(homeFramePosition(1, count), last);
    for (const progress of [1 / 3, 2 / 3]) {
      const frame = homeFramePosition(progress, count);
      assert.equal(frame, Math.round(last * progress));
      assert.ok(homeFramePosition(progress - 0.0001, count) <= frame);
      assert.ok(homeFramePosition(progress + 0.0001, count) >= frame);
      assert.ok(Math.abs(homeFramePosition(progress + 0.0001, count) - frame) < 0.01);
    }
  }
});

test('the image and tint finish together only when the profile arrives', () => {
  assert.equal(homeBackground(1, 1, 0).frameProgress, 2 / 3);
  const before = homeBackground(1, 1, 0.99);
  const finished = homeBackground(1, 1, 1);
  assert.ok(before.frameProgress < 1);
  assert.ok(before.tintProgress < 1);
  assert.equal(finished.frameProgress, 1);
  assert.equal(finished.tintProgress, 1);
  assert.deepEqual(homeBackground(2, 3, 4), finished);
});

test('the left-to-right tint expands gradually throughout playback without clearing again', () => {
  assert.equal(homeBackground(0, 0, 0).tintProgress, 0);
  const positions = [[0, 0, 0], [0.5, 0, 0], [1, 0, 0], [1, 0.5, 0], [1, 1, 0], [1, 1, 0.5], [1, 1, 1]];
  const expansion = positions.map(progress => homeBackground(...progress).tintProgress);
  assert.ok(expansion.every((value, index) => index === 0 || value > expansion[index - 1]));
  assert.equal(homeBackground(1, 0.5, 0).tintProgress, 0.5);
  assert.equal(homeBackground(1, 1, 1).tintProgress, 1);
});

test('reversals, overscroll and jumps restore the exact opening without persistent state', () => {
  const positions = [0, 0.1, 0.3, 0.5, 0.9, 1];
  const forward = positions.map(progress => homeOpening(progress));
  const reverse = [...positions].reverse().map(progress => homeOpening(progress));
  assert.deepEqual(reverse, [...forward].reverse());
  assert.deepEqual(homeOpening(-0.5), homeOpening(0));
  assert.deepEqual(homeOpening(5), homeOpening(1));
  const timeline = [[0, -2, -3], [0.4, -1, -2], [1, -0.5, -1], [1.2, 0, -1], [2, 0.5, -1], [3, 1, -0.5], [4, 2, 0], [5, 3, 0.5], [6, 4, 1]];
  const playback = timeline.map(progress => homeBackground(...progress));
  const rewind = [...timeline].reverse().map(progress => homeBackground(...progress));
  assert.deepEqual(rewind, [...playback].reverse());
  assert.deepEqual(homeBackground(-1, -1, -1), homeBackground(0, 0, 0));
});

test('reduced motion keeps the opening and research text available with a static image', () => {
  for (const progress of [0, 0.5, 1]) {
    const state = homeOpening(progress, true);
    assert.equal(state.titleOpacity, 1);
    assert.equal(state.researchOpacity, 1);
    const background = homeBackground(progress, progress, progress, true);
    assert.equal(background.frameProgress, 0);
    assert.equal(background.tintProgress, 1);
  }
});
