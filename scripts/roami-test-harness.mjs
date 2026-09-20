import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import { performance } from 'node:perf_hooks';
import * as THREE from '../public/ppi/vendor/three/three.module.min.js';

export const html = readFileSync(new URL('../public/ppi/index.html', import.meta.url), 'utf8');
const definitions = [...html.matchAll(/^      (?:async )?function ([^(]+)[^\n]*\n[\s\S]*?^      \}/gm)]
  .map((match) => match[0]).join('\n');
const constants = [
  'RESIDUE_TO_ONE', 'PROTEIN_RESIDUE_NAMES', 'RESIDUE_TYPING_ALIASES', 'NUCLEIC_RESIDUE_NAMES',
  'AROMATIC_RESIDUE_NAMES', 'BASIC_RESIDUE_NAMES', 'WATER_RESIDUE_NAMES', 'ACIDIC_RESIDUE_NAMES',
  'HALIDE_ION_ELEMENTS', 'METAL_ELEMENTS', 'BASE_PAIRING_RESIDUE_PAIR_TOP_K',
  'NUCLEOBASE_OTHER_CONTACT_TOP_K', 'PACKING_SUPPRESSION_STRONG_FAMILIES',
  'PACKING_SUPPRESSION_INDEX_CACHE', 'PI_PI_RESIDUE_PAIR_INDEX_CACHE', 'HALOGEN_PRECEDENCE_INDEX_CACHE',
  'PANEL_MODE_CONTACT_SIGNATURE_CACHE', 'CONTACT_LIST_MODE_TO_CONTAINER_ID', 'INTERACTION_VISIBILITY_DEFAULTS',
  'contactCategoryCache', 'interactionContactIndexCache', 'interactionDisplayRegistryCache', 'aromaticRingDescriptorCache',
  'AROMATIC_RING_DESCRIPTOR_CACHE_MAX', 'PURINE_AROMATIC_RING_ATOM_NAMES',
  'PYRIMIDINE_AROMATIC_RING_ATOM_NAMES', 'AROMATIC_RING_ATOM_NAMES',
  'PI_CATION_MAX_CENTROID_DISTANCE', 'PI_CATION_MAX_PLANE_OFFSET', 'PI_CATION_LATERAL_OFFSET_RADIUS_SCALE',
  'PI_CATION_MIN_LATERAL_OFFSET_LIMIT', 'PI_CATION_MAX_LATERAL_OFFSET_LIMIT',
  'PI_PI_MIN_CENTROID_DISTANCE', 'PI_PI_MAX_CENTROID_DISTANCE', 'PI_PI_MIN_CLOSEST_ATOM_DISTANCE',
  'PI_PI_MAX_CLOSEST_ATOM_DISTANCE', 'PI_PI_MAX_PLANE_SEPARATION', 'PI_PI_STACKING_NORMAL_DOT_MIN',
  'PI_PI_T_SHAPED_NORMAL_DOT_MAX',
  'STRUCTURE_METRIC_COLOR_STOPS', 'STRUCTURE_RAINBOW_COLOR_STOPS', 'STRUCTURE_METRIC_COLORMAP_MAX_FRACTION',
  'structureVertexColorCache',
];

export function harness() {
  const context = vm.createContext({ THREE, performance, console, Map, Set, WeakMap,
    window: { __PPI_INTERACTION_DEBUG_MODE: false }, INTERACTION_DEBUG_MODE_ENABLED: false,
    state: { atoms: [], residueAtoms: new Map(), loadId: 1, structureRevision: 1 },
  });
  const run = (code) => vm.runInContext(code, context);
  constants.forEach((name) => {
    const match = html.match(new RegExp(`^      const ${name} = [\\s\\S]*?;(?= *[^\\n]*\\n)`, 'm'));
    if (!match) throw new Error(`Missing frontend constant ${name}`);
    run(match[0]);
  });
  run(definitions);
  // The scheduling primitive is separate from the normalization under test.
  run('maybeYieldAnalyzePreparation = async () => {};');
  return { run, context };
}
