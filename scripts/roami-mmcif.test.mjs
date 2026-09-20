import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import * as THREE from '../public/ppi/vendor/three/three.module.min.js';
import { parseStructureRecords, serializeStructureRecords } from '../public/ppi/structure_parser.js';

// Exercise the production inline parsers without starting the viewer or its API.
const html = readFileSync(new URL('../public/ppi/index.html', import.meta.url), 'utf8');
const functions = [...html.matchAll(/^      (?:async )?function [^\n]+\n[\s\S]*?^      \}/gm)]
  .map((match) => match[0]).join('\n');
const constants = ['RESIDUE_TO_ONE', 'PROTEIN_RESIDUE_NAMES', 'RESIDUE_TYPING_ALIASES', 'NUCLEIC_RESIDUE_NAMES']
  .map((name) => {
    const match = html.match(new RegExp(`^      const ${name}\\b[\\s\\S]*?(?=^      (?:const|function|async function) )`, 'm'));
    assert.ok(match, `Missing parser constant ${name}`);
    return match[0];
  }).join('\n');
const context = vm.createContext({ THREE, parseStructureRecords, serializeStructureRecords });
vm.runInContext(`${constants}\n${functions}`, context);
const fixture = readFileSync(new URL('./fixtures/mmcif-case.cif', import.meta.url), 'utf8');

function parse(text, options = {}) {
  context.input = text;
  context.options = options;
  // Convert VM collections to a plain snapshot for equality assertions.
  return JSON.parse(vm.runInContext(`JSON.stringify((() => {
    const parsed = parseMmcifAtoms(input, options);
    const ribbon = parseMmcif(input, options);
    return {
      atoms: parsed.atoms,
      residueKeys: [...parsed.residueAtoms.keys()],
      aliases: [...parsed.chainAliases],
      models: countMmcifAtomSiteDistinctModelNumbers(input),
      selectedModels: parsed.modelTokens,
      hasMultipleModels: parsed.hasMultipleModels,
      chains: [...ribbon.chains].map(([id, points]) => [id, points.map(p => p.toArray())]),
    };
  })())`, context));
}

function changeHeaderCase(text, transform) {
  return text.replace(/^_\S+|^loop_$/gm, transform);
}

const variants = {
  lowercase: (token) => token.toLowerCase(),
  uppercase: (token) => token.toUpperCase(),
  mixed: (token) => [...token].map((c, i) => i % 2 ? c.toUpperCase() : c.toLowerCase()).join(''),
};

test('canonical CIF retains atom values, case-distinct chains and model numbers', () => {
  const parsed = parse(fixture);
  assert.equal(parsed.atoms.length, 7);
  assert.equal(parsed.residueKeys.length, 5);
  assert.equal(parsed.chains.length, 5);
  assert.deepEqual(parsed.chains.map(([id]) => id), ['Da~m1', 'DA~m1', 'rA~m1', 'Da~m2', 'DA~m2']);
  assert.deepEqual(parsed.models, { hasModelColumn: true, distinctModelCount: 2, modelTokens: ['1', '2'] });
  assert.deepEqual(parsed.atoms.map(a => a.atomName), ['N', 'CA', 'CA', "C4'", "O5'", 'CA', 'CA']);
  assert.equal(parsed.atoms[0].bFactor, 98.5);
  assert.deepEqual([parsed.atoms[0].x, parsed.atoms[0].y, parsed.atoms[0].z], [1, 2, 3]);
});

test('Coot lowercase cartn_x/y/z headers load all atoms and ribbon chains', () => {
  const cootText = fixture.replaceAll('_atom_site.Cartn_', '_atom_site.cartn_');
  assert.deepEqual(parse(cootText), parse(fixture));
});

for (const [name, transform] of Object.entries(variants)) {
  test(`${name} atom-site tags and loop keywords preserve the full parsed structure`, () => {
    const text = changeHeaderCase(fixture, transform);
    assert.deepEqual(parse(text), parse(fixture));
    assert.deepEqual(parse(text, { selectedModelToken: '2' }), parse(fixture, { selectedModelToken: '2' }));
    assert.equal(parse(text, { selectedModelToken: '2' }).atoms.length, 2);
  });
}

test('label identifiers still work when author identifiers are absent', () => {
  const labelOnly = fixture.replace(/^_atom_site\.auth_/gm, '_unused.auth_');
  assert.equal(parse(labelOnly).atoms.length, 7);
  assert.deepEqual(parse(changeHeaderCase(labelOnly, variants.uppercase)), parse(labelOnly));
});

test('missing coordinates remain missing rather than reading a different column', () => {
  const missingX = fixture.replace('_atom_site.Cartn_x', '_unused.x');
  assert.equal(parse(missingX).atoms.length, 0);
  assert.equal(parse(missingX).chains.length, 0);
});

test('local regression CIFs load without changing their headers', { skip: !process.env.ROAMI_CIF_REGRESSION_FILE }, () => {
  const text = readFileSync(process.env.ROAMI_CIF_REGRESSION_FILE, 'utf8');
  const parsed = parse(text);
  assert.equal(parsed.atoms.length, 8660);
  assert.equal(parsed.residueKeys.length, 1102);
  assert.equal(parsed.chains.length, 9);
  assert.deepEqual(parse(changeHeaderCase(text, variants.uppercase)), parsed);
});
