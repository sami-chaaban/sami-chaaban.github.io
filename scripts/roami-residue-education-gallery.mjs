import { mkdir, writeFile } from 'node:fs/promises';
import sharp from 'sharp';
import { components, residueHarness } from './roami-residue-education-harness.mjs';

const output = new URL('../roami-tests/education-audit-2026-09-20/residue-gallery/', import.meta.url);
await mkdir(output, { recursive: true });
const { run } = residueHarness();
const codes = Object.keys(components).sort((a, b) => a.length === b.length ? a.localeCompare(b) : b.length - a.length);
const rendered = [];
for (const code of codes) {
  const markup = run(`buildResidueDiagramSvg(loadResidueFixture(${JSON.stringify(code)}))`);
  const svg = markup.match(/<svg[\s\S]*?<\/svg>/)[0];
  await writeFile(new URL(`${code}.svg`, output), svg);
  const figure = await sharp(Buffer.from(svg)).resize(360, 430, { fit: 'contain', background: '#f6fbff' }).png().toBuffer();
  const header = Buffer.from(`<svg width="400" height="46"><rect width="400" height="46" fill="#f6fbff"/><text x="20" y="30" font-size="22" font-family="sans-serif">${code} · connectivity</text></svg>`);
  const tile = await sharp({ create: { width: 400, height: 490, channels: 4, background: '#f6fbff' } })
    .composite([{ input: header, left: 0, top: 0 }, { input: figure, left: 20, top: 46 }]).png().toBuffer();
  rendered.push({ code, markup, tile });
}
for (let page = 0; page * 9 < rendered.length; page += 1) {
  const subset = rendered.slice(page * 9, (page + 1) * 9);
  await sharp({ create: { width: 1200, height: Math.ceil(subset.length / 3) * 490, channels: 4, background: '#e7edf2' } })
    .composite(subset.map(({ tile }, i) => ({ input: tile, left: (i % 3) * 400, top: Math.floor(i / 3) * 490 })))
    .png().toFile(new URL(`page-${page + 1}.png`, output).pathname);
}
await writeFile(new URL('index.html', output), `<!doctype html><meta charset="utf-8"><title>Roami residue diagram audit</title>
<style>body{font:16px system-ui;background:#edf3f8;color:#203243;margin:24px}main{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:20px}article{background:#f6fbff;border:1px solid #ccd8e2;padding:16px;border-radius:10px}svg{width:100%;height:420px}.residue-chem-diagram-note{font-size:12px;color:#52687b}h2{margin:0 0 10px}</style>
<h1>Residue diagram audit</h1><p>Actual Roami renderer; official CCD ideal coordinates and heavy atoms. Connectivity view, without bond orders, formal charges, stereochemistry or H atoms. CCD terminal leaving atoms are present in these standalone components.</p>
<main>${rendered.map(({code, markup}) => `<article><h2>${code}</h2>${markup}</article>`).join('')}</main>`);
console.log(`Rendered ${rendered.length} components to ${output.pathname}`);
