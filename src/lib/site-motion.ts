import { initDesktopScroll } from './desktop-scroll';

export function initSiteMotion() {
  const disposeScroll = initDesktopScroll();
  const reduced = window.matchMedia('(prefers-reduced-motion: reduce)');
  const finePointer = window.matchMedia('(hover: hover) and (pointer: fine)');
  const animations = new Set<Animation>();
  const pendingReveals = new Set<HTMLElement>();
  const ease = 'cubic-bezier(0.22, 1, 0.36, 1)';

  const configureScroll = () => {
    if (reduced.matches) {
      animations.forEach(animation => animation.cancel());
      animations.clear();
      pendingReveals.forEach(element => element.style.removeProperty('opacity'));
      pendingReveals.clear();
    }
  };
  configureScroll();
  reduced.addEventListener('change', configureScroll);

  // Content is visible by default, including with blocked scripts or failed animation.
  const reveals = new IntersectionObserver(entries => {
    for (const entry of entries) {
      if (!entry.isIntersecting) continue;
      reveals.unobserve(entry.target);
      const element = entry.target as HTMLElement;
      pendingReveals.delete(element);
      element.style.removeProperty('opacity');
      if (reduced.matches) continue;
      const animation = element.animate([
        { opacity: 0, transform: 'translate3d(0, 26px, 0)' },
        { opacity: 1, transform: 'translate3d(0, 0, 0)' },
      ], {
        duration: 950,
        delay: Math.min(400, Number(element.dataset.delay) || 0),
        easing: ease,
        fill: 'backwards',
      });
      animations.add(animation);
      void animation.finished.then(() => animations.delete(animation)).catch(() => animations.delete(animation));
    }
  }, { threshold: 0.08 });
  document.querySelectorAll<HTMLElement>('[data-motion]').forEach(element => {
    // Never hide already-visible content, including the opening when returning to the top.
    if (reduced.matches || element.getBoundingClientRect().top < window.innerHeight) return;
    element.style.opacity = '0';
    pendingReveals.add(element);
    reveals.observe(element);
  });

  const nav = document.querySelector<HTMLElement>('.site-nav');
  const active = nav?.querySelector<HTMLElement>('.nav-link.is-active');
  const highlight = nav?.querySelector<HTMLElement>('.nav-highlight');
  let current = active;
  let navFrame = 0;
  let navigating = false;
  function alignHighlight(animate = false) {
    if (!nav || !highlight || !nav.offsetWidth) return;
    // Initial placement and font/viewport corrections must never animate in.
    highlight.toggleAttribute('data-animate', animate && !reduced.matches);
    if (!current) { highlight.style.opacity = '0'; return; }
    highlight.style.width = `${current.offsetWidth}px`;
    highlight.style.height = `${current.offsetHeight}px`;
    highlight.style.transform = `translate3d(${current.offsetLeft}px, ${current.offsetTop}px, 0)`;
    highlight.style.opacity = '1';
    nav.dataset.sliding = '';
  }
  function scheduleNav() {
    cancelAnimationFrame(navFrame);
    navFrame = requestAnimationFrame(() => alignHighlight(false));
  }
  const targetLink = (event: Event) => event.target instanceof Element
    ? event.target.closest<HTMLElement>('.nav-link') : null;
  const onHover = (event: PointerEvent) => {
    if (!finePointer.matches || navigating) return;
    const link = targetLink(event);
    if (link) { current = link; alignHighlight(true); }
  };
  const onFocus = (event: FocusEvent) => {
    const link = targetLink(event);
    if (link && !navigating) { current = link; alignHighlight(true); }
  };
  const onLeave = () => {
    if (navigating) return;
    current = nav?.querySelector<HTMLElement>('.nav-link:focus-visible') ?? active;
    alignHighlight(true);
  };
  const onNavClick = (event: MouseEvent) => {
    if (event.button !== 0 || event.metaKey || event.ctrlKey || event.shiftKey || event.altKey || event.defaultPrevented) return;
    const link = targetLink(event);
    if (!link) return;
    navigating = true;
    current = link;
    alignHighlight(false);
  };
  const onFocusOut = () => requestAnimationFrame(onLeave);
  nav?.addEventListener('pointerover', onHover);
  nav?.addEventListener('pointerleave', onLeave);
  nav?.addEventListener('focusin', onFocus);
  nav?.addEventListener('focusout', onFocusOut);
  nav?.addEventListener('click', onNavClick);
  const navObserver = new ResizeObserver(scheduleNav);
  if (nav) navObserver.observe(nav);
  void document.fonts.ready.then(scheduleNav);
  scheduleNav();

  // Warm internal documents on deliberate pointer/focus intent; links stay native.
  const prefetched = new Set<string>();
  const connection = (navigator as Navigator & { connection?: { saveData?: boolean; effectiveType?: string } }).connection;
  const prefetch = (event: Event) => {
    if (connection?.saveData || /2g/.test(connection?.effectiveType ?? '')) return;
    if (!(event.target instanceof Element)) return;
    const anchor = event.target.closest<HTMLAnchorElement>('a[href]');
    if (!anchor || anchor.hasAttribute('download') || anchor.target) return;
    const url = new URL(anchor.href);
    if (url.origin !== location.origin || url.pathname === location.pathname || !url.pathname.endsWith('/')) return;
    url.hash = '';
    if (prefetched.has(url.href)) return;
    prefetched.add(url.href);
    const link = document.createElement('link');
    link.rel = 'prefetch';
    link.href = url.href;
    document.head.append(link);
  };
  document.addEventListener('pointerover', prefetch, { passive: true });
  document.addEventListener('focusin', prefetch);
  window.addEventListener('pageshow', event => {
    if (event.persisted) { navigating = false; current = active; configureScroll(); scheduleNav(); }
  });
  window.addEventListener('pagehide', event => {
    cancelAnimationFrame(navFrame);
    animations.forEach(animation => animation.cancel());
    animations.clear();
    pendingReveals.forEach(element => element.style.removeProperty('opacity'));
    pendingReveals.clear();
    reveals.disconnect();
    if (event.persisted) return;
    disposeScroll();
    navObserver.disconnect();
    reduced.removeEventListener('change', configureScroll);
    document.removeEventListener('pointerover', prefetch);
    document.removeEventListener('focusin', prefetch);
    nav?.removeEventListener('pointerover', onHover);
    nav?.removeEventListener('pointerleave', onLeave);
    nav?.removeEventListener('focusin', onFocus);
    nav?.removeEventListener('focusout', onFocusOut);
    nav?.removeEventListener('click', onNavClick);
  });
}
