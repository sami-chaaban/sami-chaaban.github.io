import { cellPath } from './home-cell';
import { HeroLoader, framePriority, loadDecodedFrame, useStaticHero } from './hero-loader';
import type { ConnectionHints, DecodedFrame } from './hero-loader';
import { HeroRenderer } from './hero-renderer';
import { HeroFrameMotion } from './hero-frame-motion';
import { homeOpening } from './home-opening';
import { researchPinLayout, researchFocusOpacity } from './home-research';
import { approach, clamp, smoothstep } from './motion-math';

export function initFluidHome(home: HTMLElement) {
  const poster = home.querySelector<HTMLImageElement>('[data-scene-poster]');
  const canvas = home.querySelector<HTMLCanvasElement>('[data-scene-canvas]');
  const veil = home.querySelector<HTMLElement>('[data-scene-veil]');
  const shade = home.querySelector<HTMLElement>('[data-scene-shade]');
  const heroCopy = home.querySelector<HTMLElement>('[data-hero-copy]');
  const heroFooter = home.querySelector<HTMLElement>('[data-hero-footer]');
  const openingSequence = home.querySelector<HTMLElement>('[data-opening-sequence]');
  const openingStage = home.querySelector<HTMLElement>('[data-opening-stage]');
  const researchContent = home.querySelector<HTMLElement>('#our-research');
  const researchSequence = home.querySelector<HTMLElement>('[data-research-sequence]');
  const researchFocus = home.querySelector<HTMLElement>('[data-research-focus]');
  const researchCopy = home.querySelector<HTMLElement>('.inquiry-copy');
  const researchIntro = home.querySelector<HTMLElement>('.inquiry-intro');
  const division = home.querySelector<HTMLElement>('[data-division]');
  const divisionStage = home.querySelector<HTMLElement>('.division-stage');
  const path = home.querySelector<SVGPathElement>('[data-cell-path]');
  const progressBar = home.querySelector<HTMLElement>('[data-division-progress]');
  const cutCopy = home.querySelector<HTMLElement>('[data-division-cut]');
  const cellCopies = [...home.querySelectorAll<HTMLElement>('.division-half')];
  if (!poster || !canvas || !division) return;

  const reduced = window.matchMedia('(prefers-reduced-motion: reduce)');
  const portrait = window.matchMedia('(max-width: 780px) and (orientation: portrait)');
  const connection = (navigator as Navigator & { connection?: ConnectionHints & EventTarget }).connection;
  const frames = JSON.parse((portrait.matches ? home.dataset.mobileFrames : home.dataset.frames) ?? '[]') as string[];
  const renderer = new HeroRenderer(canvas);
  const frameMotion = new HeroFrameMotion();
  let staticBackground = useStaticHero(reduced.matches, connection) || !renderer.available;
  let posterReady = false;
  let disposed = false;
  let raf = 0;
  let previousTime = 0;
  let scrollY = Math.max(0, window.scrollY);
  let lastScrollTime = -Infinity;
  let viewportHeight = window.innerHeight;
  let openingTop = 0;
  let openingTravel = 1;
  let researchStart = 0;
  let researchTravel = 1;
  let researchCopyOffset = 0;
  let divisionTop = 0;
  let divisionTravel = 1;
  let cutEntrance = 18;
  let cell = 0;
  let lastRequestKey = '';
  let lastPaintKey = '';
  let lastCell = -1;

  const hold = () => {
    frameMotion.hold(renderer.hold() ?? 0);
    lastPaintKey = '';
  };
  const loader = new HeroLoader<DecodedFrame>({
    sources: frames,
    load: loadDecodedFrame,
    timeoutMs: 6000,
    onReady: () => { lastPaintKey = ''; wake(); },
    onFallback: () => { hold(); staticBackground = true; home.dataset.heroMode = 'static'; },
    release: frame => { if (frame.owned) URL.revokeObjectURL(frame.url); },
  });
  home.dataset.heroMode = staticBackground ? 'static' : 'loading';

  function measure() {
    viewportHeight = window.innerHeight;
    if (openingSequence && openingStage) {
      openingTop = openingSequence.getBoundingClientRect().top + window.scrollY;
      openingTravel = Math.max(1, openingSequence.offsetHeight - openingStage.offsetHeight);
    }
    if (researchSequence && researchContent) {
      const layout = researchPinLayout(viewportHeight, researchContent.offsetHeight, reduced.matches);
      const height = `${layout.height}px`;
      const top = `${layout.top}px`;
      if (researchSequence.style.getPropertyValue('--research-sequence-height') !== height) {
        researchSequence.style.setProperty('--research-sequence-height', height);
      }
      if (researchSequence.style.getPropertyValue('--research-pin-top') !== top) {
        researchSequence.style.setProperty('--research-pin-top', top);
      }
      researchStart = researchSequence.getBoundingClientRect().top + window.scrollY - layout.top;
      researchTravel = Math.max(1, layout.hold);
    }
    // Center only the visible introductory paragraphs until the final paragraph enters.
    researchCopyOffset = researchCopy && researchIntro && !reduced.matches && window.innerWidth > 1000
      ? Math.max(0, (researchCopy.offsetHeight - researchIntro.offsetHeight) / 2) : 0;
    const narrowDivision = window.matchMedia('(max-width: 1000px)').matches;
    cutEntrance = narrowDivision ? 8 : 18;
    let divisionPinTop = 0;
    if (narrowDivision && divisionStage && !reduced.matches) {
      const stageHeight = divisionStage.offsetHeight;
      divisionPinTop = Math.min(0, viewportHeight - stageHeight);
      division!.style.setProperty('--division-sequence-height', `${stageHeight + viewportHeight}px`);
      division!.style.setProperty('--division-pin-top', `${divisionPinTop}px`);
    } else {
      division!.style.removeProperty('--division-sequence-height');
      division!.style.removeProperty('--division-pin-top');
    }
    // Begin division only after the stage has reached its sticky resting point.
    divisionTop = division!.getBoundingClientRect().top + window.scrollY - divisionPinTop;
    divisionTravel = Math.max(viewportHeight * 0.45, division!.offsetHeight - viewportHeight);
    if (narrowDivision) divisionTravel = viewportHeight;
    lastCell = -1;
    wake();
  }

  function requestFrames(target: number) {
    if (!posterReady || staticBackground || !frames.length) return;
    const requestKey = String(target);
    if (lastRequestKey !== requestKey) {
      lastRequestKey = requestKey;
      loader.request(framePriority(target, frames.length));
    }
  }

  function paintFrame(position: number): boolean {
    if (!posterReady || staticBackground || !frames.length) return false;
    const lowerIndex = Math.floor(position);
    const upperIndex = Math.ceil(position);
    const lower = loader.ready.get(lowerIndex);
    const upper = loader.ready.get(upperIndex);
    // Never leave an incomplete blend visible while waiting for a download.
    if (!lower || !upper) { hold(); return false; }
    const blend = position - lowerIndex;
    const paintKey = `${lowerIndex}:${upperIndex}:${blend.toFixed(5)}`;
    if (paintKey === lastPaintKey) return true;
    const painted = renderer.render(
      { index: lowerIndex, image: lower.image },
      upperIndex !== lowerIndex ? { index: upperIndex, image: upper.image } : undefined,
      blend,
    );
    if (painted) lastPaintKey = paintKey;
    return painted;
  }

  function draw(time: number) {
    raf = 0;
    if (disposed || document.hidden) return;
    scrollY = Math.max(0, window.scrollY);
    const dt = previousTime ? Math.min((time - previousTime) / 1000, 0.1) : 1 / 60;
    previousTime = time;
    const opening = homeOpening((scrollY - openingTop) / openingTravel, reduced.matches);
    const targetCell = reduced.matches ? 1 : clamp((scrollY - divisionTop) / divisionTravel);
    cell = reduced.matches ? 1 : approach(cell, targetCell, dt, 0.12);
    if (Math.abs(cell - targetCell) < 0.00005) cell = targetCell;
    const scrolling = time - lastScrollTime < 100;
    const targetFrame = opening.frameProgress * Math.max(0, frames.length - 1);
    requestFrames(Math.round(targetFrame));
    const painted = paintFrame(frameMotion.advance(targetFrame, dt, scrolling));
    if (shade) shade.style.opacity = String(opening.shadeOpacity);
    if (veil) veil.style.opacity = String(opening.veilOpacity);
    if (heroCopy) {
      heroCopy.style.opacity = String(opening.titleOpacity);
      heroCopy.style.transform = reduced.matches ? '' : `translate3d(0, ${(1 - opening.titleOpacity) * -44}px, 0)`;
    }
    if (heroFooter) {
      heroFooter.style.opacity = String(opening.titleOpacity);
      heroFooter.inert = opening.titleOpacity < 0.1;
      heroFooter.setAttribute('aria-hidden', String(opening.titleOpacity < 0.1));
    }
    if (researchContent) {
      researchContent.style.opacity = String(opening.researchOpacity);
      researchContent.inert = opening.researchOpacity < 0.1;
      researchContent.setAttribute('aria-hidden', String(opening.researchOpacity < 0.1));
    }
    if (researchFocus) {
      const opacity = researchFocusOpacity((scrollY - researchStart) / researchTravel, reduced.matches);
      researchFocus.style.opacity = String(opacity);
      researchFocus.style.transform = reduced.matches ? '' : `translate3d(0, ${(1 - opacity) * 18}px, 0)`;
      researchCopy?.style.setProperty('--research-copy-offset', `${researchCopyOffset}px`);
      researchFocus.setAttribute('aria-hidden', String(opacity === 0));
    }
    if (lastCell !== cell) {
      path?.setAttribute('d', cellPath(cell));
      if (progressBar) progressBar.style.transform = `scaleX(${cell})`;
      if (cutCopy) {
        const opacity = reduced.matches ? 1 : smoothstep(0.35, 0.49, cell);
        cutCopy.style.opacity = String(opacity);
        cutCopy.style.transform = reduced.matches ? '' : `translate3d(0, ${cutEntrance * (1 - opacity)}px, 0)`;
        cutCopy.setAttribute('aria-hidden', String(opacity === 0));
      }
      lastCell = cell;
    }
    if (scrolling || cell !== targetCell || (painted && frameMotion.moving)) raf = requestAnimationFrame(draw);
    else previousTime = 0;
  }

  function wake() {
    if (!disposed && !document.hidden && !raf) raf = requestAnimationFrame(draw);
  }
  const onScroll = () => { lastScrollTime = performance.now(); wake(); };
  const onVisibility = () => {
    if (document.hidden) { loader.pause(); cancelAnimationFrame(raf); raf = 0; previousTime = 0; }
    else { lastRequestKey = ''; wake(); }
  };
  const onPreference = () => {
    if (useStaticHero(reduced.matches, connection)) loader.stop();
    measure();
  };
  const onPortraitChange = () => {
    loader.stop();
    void poster!.decode().then(() => renderer.reset()).catch(() => undefined);
    measure();
  };
  window.addEventListener('scroll', onScroll, { passive: true });
  window.addEventListener('resize', measure, { passive: true });
  reduced.addEventListener('change', onPreference);
  portrait.addEventListener('change', onPortraitChange);
  connection?.addEventListener('change', onPreference);
  document.addEventListener('visibilitychange', onVisibility);
  const observer = new ResizeObserver(measure);
  observer.observe(home);
  if (researchContent) observer.observe(researchContent);
  if (researchCopy) observer.observe(researchCopy);
  if (researchIntro) observer.observe(researchIntro);
  cellCopies.forEach(copy => observer.observe(copy));
  measure();
  void poster.decode().then(() => {
    if (disposed) return;
    posterReady = true;
    if (staticBackground) return;
    loader.ready.set(0, { url: poster.currentSrc || poster.src, image: poster, owned: false });
    home.dataset.heroMode = 'adaptive';
    wake();
  }).catch(() => { loader.stop(); });
  window.addEventListener('pageshow', event => {
    if (event.persisted) { scrollY = window.scrollY; lastRequestKey = ''; measure(); }
  });
  window.addEventListener('pagehide', event => {
    cancelAnimationFrame(raf); raf = 0; previousTime = 0;
    loader.pause();
    if (event.persisted) return;
    disposed = true;
    loader.dispose();
    observer.disconnect();
    window.removeEventListener('scroll', onScroll);
    window.removeEventListener('resize', measure);
    document.removeEventListener('visibilitychange', onVisibility);
    reduced.removeEventListener('change', onPreference);
    portrait.removeEventListener('change', onPortraitChange);
    connection?.removeEventListener('change', onPreference);
  });
}
