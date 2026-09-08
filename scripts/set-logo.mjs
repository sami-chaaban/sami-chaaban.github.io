#!/usr/bin/env node
import { readFile, mkdir, writeFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import path from 'node:path';
import sharp from 'sharp';

const publicDir = fileURLToPath(new URL('../public/', import.meta.url));

// ICO directory entries point to embedded PNG images, preserving their alpha.
function makeIco(images) {
  const directory = Buffer.alloc(6 + 16 * images.length);
  directory.writeUInt16LE(1, 2);
  directory.writeUInt16LE(images.length, 4);
  let offset = directory.length;
  images.forEach(({ size, data }, index) => {
    const entry = 6 + index * 16;
    directory[entry] = size === 256 ? 0 : size;
    directory[entry + 1] = size === 256 ? 0 : size;
    directory.writeUInt16LE(1, entry + 4);
    directory.writeUInt16LE(32, entry + 6);
    directory.writeUInt32LE(data.length, entry + 8);
    directory.writeUInt32LE(offset, entry + 12);
    offset += data.length;
  });
  return Buffer.concat([directory, ...images.map(({ data }) => data)]);
}

async function main() {
  const args = process.argv.slice(2);
  if (args.length !== 1) throw new Error('Usage: npm run set-logo -- <path-to-logo.png>');
  // Read before writing so an existing generated logo can also be used as input.
  const source = await readFile(path.resolve(args[0]));
  if ((await sharp(source).metadata()).format !== 'png') {
    throw new Error('The input must be a PNG image.');
  }
  const sizes = [16, 32, 48, 64, 180, 192, 512];
  const images = await Promise.all(sizes.map(async (size) => ({
    size,
    data: await sharp(source).rotate().resize(size, size, {
      fit: 'contain',
      background: { r: 0, g: 0, b: 0, alpha: 0 },
    }).png().toBuffer(),
  })));
  const bySize = new Map(images.map(({ size, data }) => [size, data]));
  // Complete decoding and conversion before replacing any site assets.
  const outputs = new Map([
    ['images/site-logo-source.png', source],
    ['images/site-logo.png', bySize.get(512)],
    ...[16, 32, 180, 192, 512].map((size) => [`icon-${size}.png`, bySize.get(size)]),
    ['favicon.ico', makeIco(images.filter(({ size }) => size <= 64))],
  ]);
  await mkdir(path.join(publicDir, 'images'), { recursive: true });
  for (const [name, data] of outputs) {
    await writeFile(path.join(publicDir, name), data);
  }
  console.log('Updated site logo, PNG icons, and favicon.ico (transparency preserved).');
  console.log('Original saved as public/images/site-logo-source.png.');
}

main().catch((error) => {
  console.error(`Logo update failed: ${error.message}`);
  process.exitCode = 1;
});
