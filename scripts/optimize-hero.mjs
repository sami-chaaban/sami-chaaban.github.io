import { mkdir, readFile, readdir, stat, writeFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import { createHash } from 'node:crypto';
import sharp from 'sharp';

const root = new URL('../', import.meta.url);
const source = new URL('public/images/hero-animation/', root);
const files = (await readdir(source)).filter((name) => /^\d+\.(png|jpe?g)$/i.test(name)).sort();
if (!files.length) throw new Error('No hero source frames found');

const manifest = {};
for (const variant of ['mobile', 'desktop']) {
  const destination = new URL(`public/images/hero-optimized/${variant}/`, root);
  await mkdir(destination, { recursive: true });
  let bytes = 0;
  manifest[variant] = [];
  for (const name of files) {
    const input = new URL(name, source);
    const buffer = await readFile(input);
    const settings = { variant, quality: 86, mobileCropWidth: 960, revision: 2 };
    const hash = createHash('sha256').update(buffer).update(JSON.stringify(settings)).digest('hex').slice(0, 12);
    const filename = `${name.replace(/\.[^.]+$/, '')}-${hash}.webp`;
    const output = new URL(filename, destination);
    const existing = await stat(output).catch(() => null);
    if (!existing) {
      let image = sharp(buffer);
      if (variant === 'mobile') {
        // Keep native pixel detail; discard only the sides outside a portrait viewport.
        const { width, height } = await image.metadata();
        const cropWidth = Math.min(settings.mobileCropWidth, width);
        image = image.extract({ left: Math.round((width - cropWidth) / 2), top: 0, width: cropWidth, height });
      }
      await image.webp({ quality: settings.quality, effort: 4 }).toFile(fileURLToPath(output));
    }
    bytes += (await stat(output)).size;
    manifest[variant].push(`/images/hero-optimized/${variant}/${filename}`);
  }
  console.log(`Hero ${variant}: ${files.length} frames, ${(bytes / 1e6).toFixed(2)} MB`);
}
const manifestPath = new URL('public/images/hero-optimized/manifest.json', root);
const serialized = `${JSON.stringify(manifest, null, 2)}\n`;
if (await readFile(manifestPath, 'utf8').catch(() => '') !== serialized) await writeFile(manifestPath, serialized);
