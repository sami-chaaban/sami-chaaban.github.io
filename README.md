# sami-chaaban.github.io

Personal website built with Astro.

## Hero animation frames

Original numbered frames live in `assets/hero-animation/`, outside `public/` so
they are not copied into the deployed site. Keep this folder in the repository
because both `npm run dev` and `npm run build` use it to generate the mobile and
desktop WebP frames and manifest in `public/images/hero-optimized/`.

The September 9 fluidity trial uses native wheel, touch, and keyboard scrolling.
The opening is visible immediately. The molecular sequence starts on the first
scroll movement while the title fades out, and spans roughly 0.6 viewports of
scrolling with a 70% transparent filter and no extra gradient shading. The tint
lifts as the title leaves. The research text begins appearing at 32% of the opening
sequence, with the darker reading tint following at 54%. The last molecular frame
and the completed return of the reading filter both occur at 76%. “Explore our research”
is an initial anchor shortcut that fades out with the title. The research section
overlaps the end of the opening stage to avoid an empty intermediate screen,
and a cell diagram divides before the profile and links. The entire opening
reverses with native scrolling; reduced motion removes the extended sticky stage.
The research section then pins for half a viewport of scrolling. Its final “A central
focus…” paragraph fades in throughout that entire hold, beginning immediately,
with light easing between wheel ticks and a small upward movement. On desktop, the
initial right-hand paragraphs are centered against the heading and remain fixed
as the final paragraph enters from below. Tall layouts scroll fully into view before
pinning their lower edge; reduced motion shows the paragraph immediately without
adding the hold. The state before this change is saved in
`.design-trials/before-research-pinning-2026-09-09/`.
On narrow layouts, the entire cell fits inside the viewport and pins for one
viewport of continued scrolling. The second half starts revealing at 35% of this
hold. The SVG retains equal scale on both axes, so division ends in circles.
Portrait layouts stack the circles; landscape layouts place them side by side.
Body text remains 16px, with compact spacing and decorative chapter numbers
omitted on short screens. Reduced motion shows both halves without the hold.
The persistent canvas follows continuous scroll positions with elapsed-time easing,
then settles on one complete image after scrolling stops. It does not restart its
easing at each image boundary. While an image loads, the last complete frame stays
visible. Independent camera transforms and synthetic wheel smoothing are removed.
The original 19 images remain the
motion source; this is not a live 3D molecular model.

Shared motion adds staggered text entrances, a sliding navigation highlight,
internal page prefetch on pointer/focus intent, and softer native page transitions.
The shared `src/styles/site-theme.css` carries the toolbar geometry across pages
while preserving their original colours and layouts. Publications retains its full
background and responsive title placement; CV retains its light palette and unboxed
institution logos. Navigation highlighting
is placed immediately on load, font changes, resize, and navigation; only hover
and keyboard-focus changes animate it.
Content is visible without JavaScript. Reduced-motion users get native scrolling,
static imagery, and an unpinned division section. Slow connections retain a poster
or the last complete frame.

Run `npm run test:hero` for the existing loader/renderer and legacy entrance tests,
and `npm run test:motion` for refresh-rate independence, continuous frame movement,
sharp settling, and repeated scroll reversal. Browser checks are recorded in
`.design-trials/fluidity-browser-check-2026-09-09.md`.
The report includes a follow-up browser check of the opening interlude and the
full-hold research fade on desktop and narrow layouts. Browser access worked on
retry; the local preview server was restarted before checking. The interlude's
preceding source is saved in
`.design-trials/before-animation-interlude-2026-09-09/`.

The complete pre-trial source, assets, previous edits, and Git history are archived
outside this project in `../backups/before-fluidity-trial-2026-09-09-095705/`.
See that folder's `RESTORE.md` and verified SHA-256 manifest. This trial is local;
no deployment has been performed.

## Update the site logo

From `my-site`, run:

```sh
npm run set-logo -- ../logo5.png
```

Pass exactly one PNG path (quote paths containing spaces). Run `npm install`
first if dependencies are not installed. The script preserves transparency and
aspect ratio, adding transparent padding to make square icons.

The original is copied to `public/images/site-logo-source.png`. The script
regenerates `public/images/site-logo.png`, `public/icon-{16,32,180,192,512}.png`,
and the multi-size `public/favicon.ico`. These fixed filenames are already used
by the site and PPI viewer, so no template edits are needed when changing logos.

Run `npm run build` when preparing the updated site for deployment. Browsers may
cache favicons, so an old tab icon may persist until the cache is refreshed.
