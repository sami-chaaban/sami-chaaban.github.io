import { mkdir, readFile, readdir, stat, writeFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import { createHash } from 'node:crypto';
import sharp from 'sharp';

const root = new URL('../', import.meta.url);
const source = new URL('assets/hero-animation/', root);
const files = (await readdir(source)).filter((name) => /^\d+\.(png|jpe?g)$/i.test(name)).sort();
if (!files.length) throw new Error('No hero source frames found');

const mobileCropWidth = 960;
// Keep native pixel detail; discard only the sides outside a portrait viewport.
const mobileCrop = (width, height) => {
  const cropWidth = Math.min(mobileCropWidth, width);
  return { left: Math.round((width - cropWidth) / 2), top: 0, width: cropWidth, height };
};

const manifest = {};
for (const variant of ['mobile', 'desktop']) {
  const destination = new URL(`public/images/hero-optimized/${variant}/`, root);
  await mkdir(destination, { recursive: true });
  let bytes = 0;
  manifest[variant] = [];
  for (const name of files) {
    const input = new URL(name, source);
    const buffer = await readFile(input);
    const settings = { variant, quality: 86, mobileCropWidth, revision: 2 };
    const hash = createHash('sha256').update(buffer).update(JSON.stringify(settings)).digest('hex').slice(0, 12);
    const filename = `${name.replace(/\.[^.]+$/, '')}-${hash}.webp`;
    const output = new URL(filename, destination);
    const existing = await stat(output).catch(() => null);
    if (!existing) {
      let image = sharp(buffer);
      if (variant === 'mobile') {
        const { width, height } = await image.metadata();
        image = image.extract(mobileCrop(width, height));
      }
      await image.webp({ quality: settings.quality, effort: 4 }).toFile(fileURLToPath(output));
    }
    bytes += (await stat(output)).size;
    manifest[variant].push(`/images/hero-optimized/${variant}/${filename}`);
  }
  console.log(`Hero ${variant}: ${files.length} frames, ${(bytes / 1e6).toFixed(2)} MB`);
}

// Depth-warped playback (src/lib/hero-depth-renderer.ts) needs only keyframes: their
// colour frames above plus each one's inverse depth, re-projected in WebGL. The depth
// sources are rendered by the hero-depth scripts kept beside the Blender scene.
const depthSource = new URL('assets/hero-depth/', root);
const camera = JSON.parse(await readFile(new URL('camera.json', depthSource), 'utf8').catch(() => 'null'));
if (camera) {
  const keyframeStep = 4;
  const last = files.length - 1;
  const keyframes = [...new Set(Array.from({ length: Math.floor(last / keyframeStep) + 1 }, (_, i) => i * keyframeStep).concat(last))];
  manifest.depth = { keyframes, camera, views: {} };
  for (const variant of ['mobile', 'desktop']) {
    const destination = new URL(`public/images/hero-optimized/depth-${variant}/`, root);
    await mkdir(destination, { recursive: true });
    let bytes = 0;
    manifest.depth[variant] = [];
    manifest.depth.views[variant] = variant === 'mobile'
      ? mobileCrop(camera.width, camera.height) : { left: 0, top: 0, width: camera.width, height: camera.height };
    for (const index of keyframes) {
      const name = `${files[index].replace(/\.[^.]+$/, '')}.png`;
      const buffer = await readFile(new URL(name, depthSource));
      const settings = { variant, encoding: 'inverse-depth-8bit', mobileCropWidth, revision: 1 };
      const hash = createHash('sha256').update(buffer).update(JSON.stringify(settings)).digest('hex').slice(0, 12);
      const filename = `${name.replace(/\.[^.]+$/, '')}-${hash}.webp`;
      const output = new URL(filename, destination);
      const existing = await stat(output).catch(() => null);
      if (!existing) {
        // 8 bits of inverse depth move parallax by well under a pixel between keyframes.
        const { data, info } = await sharp(buffer).toColourspace('grey16').raw({ depth: 'ushort' })
          .toBuffer({ resolveWithObject: true });
        const depth = Buffer.alloc(info.width * info.height);
        // Raw samples come out in machine byte order, little-endian on supported hosts.
        for (let i = 0; i < depth.length; i++) depth[i] = Math.round(data.readUInt16LE(i * 2) / 257);
        let image = sharp(depth, { raw: { width: info.width, height: info.height, channels: 1 } });
        if (variant === 'mobile') image = image.extract(mobileCrop(info.width, info.height));
        // Lossless: compression artefacts would tear surfaces apart at silhouettes.
        await image.webp({ lossless: true, effort: 4 }).toFile(fileURLToPath(output));
      }
      bytes += (await stat(output)).size;
      manifest.depth[variant].push(`/images/hero-optimized/depth-${variant}/${filename}`);
    }
    console.log(`Hero depth ${variant}: ${keyframes.length} keyframes, ${(bytes / 1e6).toFixed(2)} MB`);
  }
}
const manifestPath = new URL('public/images/hero-optimized/manifest.json', root);
const serialized = `${JSON.stringify(manifest, null, 2)}\n`;
if (await readFile(manifestPath, 'utf8').catch(() => '') !== serialized) await writeFile(manifestPath, serialized);
