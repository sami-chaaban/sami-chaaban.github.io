import { readFileSync } from 'node:fs';
import { harness, html } from './roami-test-harness.mjs';

// Independent wwPDB CCD atom identities, ideal coordinates and bond topology.
// Source URLs and downloaded-definition hashes accompany the fixture.
export const components = JSON.parse(readFileSync(new URL(
  '../roami-tests/education-audit-2026-09-20/ccd/heavy-atom-fixtures.json', import.meta.url,
), 'utf8'));

export function residueHarness() {
  const result = harness();
  const constants = [
    'RESIDUE_INTERFACE_TOOLTIP_CONTENT', 'AMINO_ACID_REFERENCE',
    'CANONICAL_AMINO_HEAVY_ATOM_REQUIREMENTS', 'CANONICAL_AMINO_SIDECHAIN_BOND_SPECS',
    'CANONICAL_AMINO_INTRA_RESIDUE_BOND_SPECS', 'CANONICAL_AMINO_INTRA_RESIDUE_BOND_NAME_PAIR_LOOKUP',
    'CANONICAL_AMINO_HEAVY_ATOM_NAME_LOOKUP', 'CANONICAL_AROMATIC_RING_BOND_SPECS',
    'CANONICAL_AROMATIC_RING_BOND_NAME_PAIR_LOOKUP', 'CANONICAL_AROMATIC_RING_ATOM_NAME_LOOKUP',
    'CANONICAL_NUCLEIC_INTRA_RESIDUE_BOND_SPECS', 'CANONICAL_NUCLEIC_INTRA_RESIDUE_BOND_NAME_PAIR_LOOKUP',
    'CANONICAL_NUCLEIC_HEAVY_ATOM_NAME_LOOKUP', 'DIAGRAM_SULFUR_LIGHT_COLOR', 'SULFUR_ATOM_COLOR',
  ];
  constants.forEach((name) => {
    const match = html.match(new RegExp(`^      const ${name} =[\\s\\S]*?(?=^      (?:const|let|function|async function) )`, 'm'));
    if (!match) throw new Error(`Missing frontend constant ${name}`);
    result.run(match[0]);
  });
  result.run("state.sceneTheme = 'light'; state.explicitAtomBonds = [];");
  result.context.ccd = components;
  result.run(`function loadResidueFixture(code, explicit = false) {
    const component = ccd[code];
    state.atoms = component.atoms.map((atom, index) => ({
      ...atom, index, serial: index + 1, chain: 'A', seq: '1',
      resKey: 'A:1', resName: code,
    }));
    state.residueAtoms = new Map([['A:1', state.atoms.map((_, index) => index)]]);
    state.explicitAtomBonds = explicit ? component.bonds.map(({atoms}) =>
      atoms.map(name => state.atoms.findIndex(atom => atom.atomName === name))) : [];
    return { key: 'A:1', resName: code, atoms: state.atoms, isAminoAcid: isProteinResidueName(code) };
  }`);
  return result;
}
