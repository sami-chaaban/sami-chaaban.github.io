import test from 'node:test';
import assert from 'node:assert/strict';
import { HeroLoader, framePriority, useStaticHero, loadDecodedFrame } from '../src/lib/hero-loader.ts';

const flush = () => new Promise((resolve) => setImmediate(resolve));
const sources = ['0', '1', '2', '3', '4'];

test('data saver, reduced motion and slow connections use the static hero', () => {
  assert.equal(useStaticHero(true), true);
  for (const hints of [{ saveData: true }, { effectiveType: 'slow-2g' }, { effectiveType: '2g' }]) {
    assert.equal(useStaticHero(false, hints), true);
  }
  assert.equal(useStaticHero(false), false);
  // Rough estimates alone must not disable animation; actual load deadlines handle slow links.
  assert.equal(useStaticHero(false, { effectiveType: '3g', downlink: 1.5 }), false);
  assert.equal(useStaticHero(false, { effectiveType: '4g', downlink: 10 }), false);
});

test('loading is deferred, bounded, and reprioritised after a rapid scroll', async () => {
  const started = [];
  const complete = new Map();
  const loader = new HeroLoader({ sources, load: (src) => {
    started.push(src);
    return new Promise((resolve) => complete.set(src, resolve));
  }, onReady() {}, onFallback() {} });
  assert.deepEqual(started, []);
  loader.request(framePriority(0, sources.length));
  await flush();
  assert.deepEqual(started, ['0', '1']);
  assert.equal(loader.ready.size, 0);
  loader.request(framePriority(4, sources.length));
  complete.get('0')('decoded-0');
  await flush();
  assert.deepEqual(started, ['0', '1', '4']);
  assert.equal(loader.ready.get(0), 'decoded-0');
  loader.dispose();
  complete.get('1')('decoded-1');
  complete.get('4')('decoded-4');
  await flush();
});

test('a slow image aborts remaining requests and retains the last good frame', async () => {
  const signals = [];
  let fallbacks = 0;
  const loader = new HeroLoader({ sources, timeoutMs: 15, load: (_, signal) => {
    signals.push(signal);
    return new Promise((_, reject) => signal.addEventListener('abort', () => reject(new Error('aborted')), { once: true }));
  }, onReady() {}, onFallback() { fallbacks++; } });
  loader.ready.set(0, 'poster');
  loader.request(framePriority(1, sources.length));
  await new Promise((resolve) => setTimeout(resolve, 40));
  assert.equal(fallbacks, 1);
  assert.equal(signals.length, 2);
  assert(signals.every((signal) => signal.aborted));
  assert.equal(loader.ready.get(0), 'poster');
  loader.request([4]);
  await flush();
  assert.equal(signals.length, 2);
  loader.dispose();
});

test('failed images never become ready and repeated errors stop background downloads', async () => {
  let fallbacks = 0;
  const loader = new HeroLoader({ sources, concurrency: 1, load: async () => { throw new Error('404'); },
    onReady() { assert.fail('failed frames must not be displayed'); }, onFallback() { fallbacks++; } });
  loader.ready.set(0, 'poster');
  loader.request([1, 2, 3, 4]);
  await flush();
  assert.equal(fallbacks, 1);
  assert.deepEqual([...loader.ready.values()], ['poster']);
  loader.dispose();
});

test('leaving the visible animation cancels pending work and scrolling back resumes it', async () => {
  const started = [];
  const loader = new HeroLoader({ sources, load: (source, signal) => {
    started.push(source);
    return new Promise((_, reject) => signal.addEventListener('abort', () => reject(new Error('aborted')), { once: true }));
  }, onReady() {}, onFallback() { assert.fail('intentional pauses are not failures'); } });
  loader.request([0, 1, 2]);
  await flush();
  loader.pause();
  await flush();
  assert.deepEqual(started, ['0', '1']);
  loader.request([4, 3]);
  await flush();
  assert.deepEqual(started, ['0', '1', '4', '3']);
  loader.dispose();
  await flush();
});

test('late results after disposal are released rather than displayed', async () => {
  let complete;
  const released = [];
  const loader = new HeroLoader({ sources: ['0'], load: () => new Promise((resolve) => { complete = resolve; }),
    onReady() { assert.fail('disposed frames must not be displayed'); }, onFallback() {}, release: (frame) => released.push(frame) });
  loader.request([0]);
  await flush();
  loader.dispose();
  complete('late-frame');
  await flush();
  assert.deepEqual(released, ['late-frame']);
  assert.equal(loader.ready.size, 0);
});

test('download completion alone does not expose a frame before decoding finishes', async (t) => {
  let finishDecode;
  let startedDecode;
  const decodingStarted = new Promise((resolve) => { startedDecode = resolve; });
  const originalImage = globalThis.Image;
  globalThis.Image = class {
    decode() { return new Promise((resolve) => { finishDecode = resolve; startedDecode(); }); }
  };
  t.after(() => { if (originalImage) globalThis.Image = originalImage; else delete globalThis.Image; });
  t.mock.method(globalThis, 'fetch', async () => new Response(new Blob(['frame']), { status: 200 }));
  let resolved = false;
  const promise = loadDecodedFrame('/test.webp', new AbortController().signal).then((frame) => { resolved = true; return frame; });
  await decodingStarted;
  assert.equal(resolved, false);
  finishDecode();
  const frame = await promise;
  assert.equal(resolved, true);
  assert(frame.url.startsWith('blob:'));
  URL.revokeObjectURL(frame.url);
});
