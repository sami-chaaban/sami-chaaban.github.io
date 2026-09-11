import { ScrollMomentum } from './scroll-momentum';

export function initDesktopScroll() {
  const finePointer = window.matchMedia('(hover: hover) and (pointer: fine)');
  const reduced = window.matchMedia('(prefers-reduced-motion: reduce)');
  const motion = new ScrollMomentum();
  let frame = 0;
  let previousTime = 0;
  let writtenPosition = window.scrollY;
  const maximum = () => Math.max(0, document.documentElement.scrollHeight - window.innerHeight);

  const stop = () => {
    cancelAnimationFrame(frame);
    frame = 0;
    previousTime = 0;
    motion.stop(window.scrollY);
  };

  const animate = (time: number) => {
    const seconds = previousTime ? Math.min((time - previousTime) / 1000, 0.05) : 1 / 60;
    previousTime = time;
    writtenPosition = motion.advance(seconds, maximum());
    window.scrollTo({ top: writtenPosition, behavior: 'instant' });
    writtenPosition = window.scrollY;
    if (motion.moving) frame = requestAnimationFrame(animate);
    else { frame = 0; previousTime = 0; }
  };

  const needsNativeScroll = (event: WheelEvent) => {
    for (const element of event.composedPath()) {
      if (!(element instanceof Element)) continue;
      if (element === document.body || element === document.documentElement) break;
      if (element.matches('input, textarea, select, [contenteditable]:not([contenteditable="false"]), [role="slider"], [role="spinbutton"], dialog, [role="dialog"], [data-native-scroll]')) return true;
      if (element.scrollHeight > element.clientHeight + 1 && /^(auto|scroll|overlay)$/.test(getComputedStyle(element).overflowY)) return true;
    }
    return false;
  };

  const onWheel = (event: WheelEvent) => {
    if (event.defaultPrevented || !event.cancelable || event.ctrlKey || event.metaKey || event.altKey || event.shiftKey
      || !event.deltaY || Math.abs(event.deltaX) > Math.abs(event.deltaY) || needsNativeScroll(event)) {
      stop();
      return;
    }
    const delta = event.deltaY * (event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? window.innerHeight : 1);
    if (!frame) motion.stop(window.scrollY);
    motion.push(delta, maximum());
    if (!motion.moving) return;
    event.preventDefault();
    if (!frame) frame = requestAnimationFrame(animate);
  };

  // Yield to anchors, keyboard navigation, scrollbar dragging and scroll restoration.
  const onScroll = () => {
    if (frame && Math.abs(window.scrollY - writtenPosition) > 1) stop();
  };
  const onVisibility = () => { if (document.hidden) stop(); };
  const configure = () => {
    stop();
    window.removeEventListener('wheel', onWheel);
    if (finePointer.matches && !reduced.matches) window.addEventListener('wheel', onWheel, { passive: false });
  };

  configure();
  finePointer.addEventListener('change', configure);
  reduced.addEventListener('change', configure);
  window.addEventListener('scroll', onScroll, { passive: true });
  window.addEventListener('resize', stop, { passive: true });
  window.addEventListener('pointerdown', stop, { passive: true, capture: true });
  window.addEventListener('touchstart', stop, { passive: true, capture: true });
  window.addEventListener('keydown', stop, true);
  window.addEventListener('click', stop, true);
  document.addEventListener('visibilitychange', onVisibility);
  window.addEventListener('pageshow', configure);
  window.addEventListener('pagehide', stop);

  return () => {
    stop();
    window.removeEventListener('wheel', onWheel);
    finePointer.removeEventListener('change', configure);
    reduced.removeEventListener('change', configure);
    window.removeEventListener('scroll', onScroll);
    window.removeEventListener('resize', stop);
    window.removeEventListener('pointerdown', stop, true);
    window.removeEventListener('touchstart', stop, true);
    window.removeEventListener('keydown', stop, true);
    window.removeEventListener('click', stop, true);
    document.removeEventListener('visibilitychange', onVisibility);
    window.removeEventListener('pageshow', configure);
    window.removeEventListener('pagehide', stop);
  };
}
