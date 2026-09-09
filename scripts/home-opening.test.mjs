import test from 'node:test';
import assert from 'node:assert/strict';
import { homeOpening } from '../src/lib/home-opening.ts';

test('molecular playback begins on the first scroll movement while the title is still visible', () => {
  assert.equal(homeOpening(0).frameProgress, 0);
  const firstMovement = homeOpening(1 / 1000);
  assert.ok(firstMovement.frameProgress > 0);
  assert.ok(firstMovement.titleOpacity > 0);
});

test('research appears earlier while molecular playback and the delayed shading continue', () => {
  assert.equal(homeOpening(1).researchOpacity, 1);
  assert.equal(homeOpening(0.32).researchOpacity, 0);
  assert.ok(homeOpening(0.4).researchOpacity > 0);
  assert.equal(homeOpening(0.52).researchOpacity, 1);
  assert.ok(homeOpening(0.52).frameProgress < 1);
  assert.ok(Math.abs(homeOpening(0.54).veilOpacity - 0.3) < 1e-12);
});

test('the last frame coincides with the reading filter finishing its return', () => {
  const before = homeOpening(0.97);
  const finished = homeOpening(0.98);
  assert.ok(before.frameProgress < 1);
  assert.ok(before.veilOpacity < finished.veilOpacity);
  assert.equal(finished.frameProgress, 1);
  assert.ok(Math.abs(finished.veilOpacity - 0.35) < 1e-12);
  assert.equal(finished.shadeOpacity, 1);
  assert.equal(homeOpening(1).frameProgress, 1);
});

test('the image is revealed during playback, then shaded again for the research text', () => {
  const opening = homeOpening(0);
  const interlude = homeOpening(0.45);
  const research = homeOpening(1);
  assert.equal(opening.veilOpacity, 0.35);
  assert.ok(interlude.veilOpacity < opening.veilOpacity);
  assert.ok(Math.abs(interlude.veilOpacity - 0.3) < 1e-12, 'filter is 70% transparent');
  assert.ok(homeOpening(0.12).veilOpacity < opening.veilOpacity, 'filter lifts early in the scroll');
  assert.ok(homeOpening(0.22).veilOpacity > interlude.veilOpacity, 'fade-out is still progressing at its old endpoint');
  assert.ok(Math.abs(homeOpening(0.38).veilOpacity - interlude.veilOpacity) < 1e-12);
  assert.equal(interlude.shadeOpacity, 0, 'no extra gradient obscures the animation');
  assert.ok(research.veilOpacity > interlude.veilOpacity);
  assert.ok(homeOpening(0.76).veilOpacity < research.veilOpacity, 'fade-in is still progressing at its old endpoint');
  assert.ok(Math.abs(research.veilOpacity - opening.veilOpacity) < 1e-12);
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
    assert.equal(state.veilOpacity, 0.35);
  }
});
