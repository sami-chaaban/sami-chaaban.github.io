import test from 'node:test';
import assert from 'node:assert/strict';
import { HeroEntrance, heroEntranceOpacity } from '../src/lib/hero-entrance.ts';

test('the background fades from ten percent visibility and returns there when rewound', () => {
  const forward = [0, 0.025, 0.05, 0.1, 1].map((progress) => heroEntranceOpacity(progress));
  assert.equal(forward[0], 0.9);
  assert.equal(forward.at(-1), 0);
  assert(forward[1] > forward[2] && forward[2] > forward[3]);
  const backward = [1, 0.1, 0.05, 0.025, 0].map((progress) => heroEntranceOpacity(progress));
  assert.deepEqual(backward, [...forward].reverse());
  assert.equal(heroEntranceOpacity(-1), 0.9);
  assert.equal(heroEntranceOpacity(3), 0);
});

test('reduced motion still supports returning to the dim opening without a fade', () => {
  assert.equal(heroEntranceOpacity(0, true), 0.9);
  assert.equal(heroEntranceOpacity(0.01, true), 0);
  assert.equal(heroEntranceOpacity(0, true), 0.9);
});

test('preloading never reveals the page until the first gesture', () => {
  const starts = [];
  const entrance = new HeroEntrance((progress) => starts.push(progress));
  entrance.markReady();
  assert.equal(entrance.hasStarted, false);
  assert.deepEqual(starts, []);
  assert.equal(entrance.request(0.05), true);
  assert.equal(entrance.hasStarted, true);
  assert.deepEqual(starts, [0.05]);
});

test('early scrolling waits for the image without skipping ahead while it loads', () => {
  const starts = [];
  const entrance = new HeroEntrance((progress) => starts.push(progress));
  entrance.request(0.05);
  entrance.request(0.2);
  entrance.request(0.5);
  assert.equal(entrance.hasStarted, false);
  assert.deepEqual(starts, []);
  entrance.markReady();
  assert.deepEqual(starts, [0.05]);
});

test('keyboard entry at zero and repeated readiness signals reveal only once', () => {
  const starts = [];
  const entrance = new HeroEntrance((progress) => starts.push(progress));
  entrance.request(0);
  entrance.markReady();
  entrance.markReady();
  assert.equal(entrance.request(0.4), false);
  assert.deepEqual(starts, [0]);
});
