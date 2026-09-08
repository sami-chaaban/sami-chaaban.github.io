# sami-chaaban.github.io

Personal website built with Astro.

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
