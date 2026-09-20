import assert from 'node:assert/strict';
import test from 'node:test';
import { components, residueHarness } from './roami-residue-education-harness.mjs';

test('canonical heavy-atom connectivity agrees with independent CCD for all 20 amino acids and six nucleotide families', () => {
  const { run } = residueHarness();
  for (const [code, component] of Object.entries(components)) {
    if (['PSU', 'SEP', 'MSE'].includes(code)) continue;
    const names = new Set(component.atoms.map(atom => atom.atomName));
    const expected = component.bonds.map(bond => bond.atoms.join('|')).sort();
    const actual = Array.from(run(`(() => {
      const code = ${JSON.stringify(code)};
      return [...(CANONICAL_AMINO_INTRA_RESIDUE_BOND_NAME_PAIR_LOOKUP[code] ||
        CANONICAL_NUCLEIC_INTRA_RESIDUE_BOND_NAME_PAIR_LOOKUP[normalizeNucleobaseFamilyToken(code)])];
    })()`)).filter(pair => pair.split('|').every(name => names.has(name))).sort();
    assert.deepEqual(actual, expected, code);
  }
});

test('actual residue diagram bond construction recovers CCD topology with and without explicit bonds', () => {
  const { run } = residueHarness();
  for (const [code, component] of Object.entries(components)) {
    const expected = component.bonds.map(bond => bond.atoms.join('|')).sort();
    for (const explicit of [false, true]) {
      const actual = Array.from(run(`buildResidueDiagramBonds(loadResidueFixture(${JSON.stringify(code)}, ${explicit}))
        .map(pair => pair.map(atom => atom.atomName).sort().join('|')).sort()`));
      assert.deepEqual(actual, expected, `${code} explicit=${explicit}`);
    }
  }
});

test('ribose shortcut is rejected even for a distorted close approach or erroneous explicit bond', () => {
  const { run } = residueHarness();
  for (const code of ['A', 'G', 'C', 'U', 'DT', 'I']) {
    const result = run(`(() => {
      const meta = loadResidueFixture(${JSON.stringify(code)});
      const a = meta.atoms.find(atom => atom.atomName === "C1'");
      const b = meta.atoms.find(atom => atom.atomName === "C4'");
      b.x = a.x + 1.5; b.y = a.y; b.z = a.z;
      state.explicitAtomBonds = [[a.index, b.index]];
      return { explicitAllowed: shouldKeepExplicitCanonicalIntraResidueBond(a, b),
        bonds: buildResidueDiagramBonds(meta).map(pair => pair.map(atom => atom.atomName).sort().join('|')) };
    })()`);
    assert.equal(result.explicitAllowed, false, code);
    assert.ok(!result.bonds.includes("C1'|C4'"), code);
  }
});

test('residue diagrams render all CCD heavy atoms with finite coordinates and preserve base atom identities', () => {
  const { run } = residueHarness();
  for (const [code, component] of Object.entries(components)) {
    const svg = run(`buildResidueDiagramSvg(loadResidueFixture(${JSON.stringify(code)}))`);
    assert.doesNotMatch(svg, /NaN|Infinity|undefined/, code);
    assert.equal(new Set([...svg.matchAll(/data-chem-atom-id="([^"]+)"/g)].map(match => match[1])).size, component.atoms.length, code);
  }
  assert.equal(run(`formatResidueDiagramAtomLabel({atomName: 'N9'}, {resName:'DA'})`), 'N9');
  assert.equal(run(`formatResidueDiagramAtomLabel({atomName: "C1'"}, {resName:'A'})`), 'C1′');
  assert.doesNotMatch(run(`buildResidueDiagramSvg({resName:'LIG', atoms:[]})`), /<circle|<line/);
});

test('education controls resolve every standard residue and nucleotide aliases', () => {
  const { run } = residueHarness();
  for (const code of [...Object.keys(components), 'DA', 'DG', 'DC', 'DU', 'ATP', 'GTP', 'HIP', 'HIE', 'HID']) {
    assert.equal(run(`Boolean(RESIDUE_INTERFACE_TOOLTIP_CONTENT[getResidueChemTooltipMode(${JSON.stringify(code)})])`), true, code);
  }
});

test('modified aliases keep their own chemistry and pseudouridine C1′–C5 attachment', () => {
  const { run } = residueHarness();
  for (const code of ['PSU', 'SEP', 'MSE']) {
    assert.equal(run(`getResidueChemTooltipMode(${JSON.stringify(code)})`), 'residue-modified', code);
  }
  assert.equal(run("getResidueChemDescriptor('S', 'SEP')"), 'Modified residue');
  assert.equal(run("formatResidueDiagramAtomLabel({atomName:'SE', element:'SE'}, {resName:'MSE', isAminoAcid:true})"), 'Se');
  assert.equal(run("formatResidueDiagramAtomLabel({atomName:'O1P', element:'O'}, {resName:'SEP', isAminoAcid:true})"), 'O1P');
  for (const explicit of [true, false]) {
    const result = run(`(() => {
      const meta = loadResidueFixture('PSU', ${explicit});
      const bonds = buildResidueDiagramBonds(meta).map(pair => pair.map(atom => atom.atomName).sort().join('|'));
      return { bonds, markup: buildResidueDiagramSvg(meta) };
    })()`);
    assert.ok(result.bonds.includes("C1'|C5"));
    assert.ok(!result.bonds.includes("C1'|N1"));
    if (explicit) assert.doesNotMatch(result.markup, /Some bonds are approximate/);
    else assert.match(result.markup, /Some bonds are approximate/);
  }
});
