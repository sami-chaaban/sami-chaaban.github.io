import test from 'node:test';
import assert from 'node:assert/strict';
import { researchPinLayout, researchFocusOpacity } from '../src/lib/home-research.ts';

test('the final paragraph starts immediately and fades throughout the entire sticky hold', () => {
  for (const progress of [-1, 0]) assert.equal(researchFocusOpacity(progress), 0);
  const opacity = [0.01, 0.25, 0.5, 0.75, 0.99].map(progress => researchFocusOpacity(progress));
  assert.ok(opacity[0] > 0);
  assert.equal(opacity[2], 0.5);
  assert.ok(opacity.at(-1) < 1);
  assert.ok(opacity.every((value, index) => index === 0 || value > opacity[index - 1]));
  assert.equal(researchFocusOpacity(1), 1);
});

test('desktop and tall narrow-screen content have the same hold after all text can enter view', () => {
  for (const [viewport, content] of [[720, 720], [844, 1160], [568, 1300]]) {
    const layout = researchPinLayout(viewport, content);
    assert.equal(layout.hold, viewport * 0.5);
    assert.equal(layout.height - content, viewport * 0.5);
    assert.ok(layout.top <= 0);
    assert.equal(layout.top + content, viewport);
  }
});

test('scrolling backwards hides the paragraph at the same point without layout changes', () => {
  const progress = [0, 0.5, 0.68, 0.76, 0.84, 1];
  assert.deepEqual([...progress].reverse().map(value => researchFocusOpacity(value)),
    progress.map(value => researchFocusOpacity(value)).reverse());
});

test('reduced motion shows all text and adds no pinned scroll distance', () => {
  assert.equal(researchFocusOpacity(0, true), 1);
  assert.deepEqual(researchPinLayout(844, 1160, true), { top: 0, height: 1160, hold: 0 });
});
