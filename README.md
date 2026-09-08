# sami-chaaban.github.io

Personal website built with Astro.

## Hero animation frames

Original numbered frames live in `assets/hero-animation/`, outside `public/` so
they are not copied into the deployed site. Keep this folder in the repository
because both `npm run dev` and `npm run build` use it to generate the mobile and
desktop WebP frames and manifest in `public/images/hero-optimized/`.

The homepage keeps its title and navigation visible over the first frame at
10% visibility while frames preload. Scrolling brings the background to full
visibility and starts the sequence; the title and institute line fade while
drifting slightly left from the first scroll. Scrolling back restores the dim opening and title.
Decoded frames are blended on a persistent canvas to avoid empty image swaps.
Run `npm run test:hero` to check loading, rendering, and entrance behaviour.

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
