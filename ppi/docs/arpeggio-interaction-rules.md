# Arpeggio Interaction Rules and Filters (Current Full Spec)

This document is the implementation-level specification for interaction classification and filtering in:

- `my-site/ppi/api/analysis.py` (backend generation, assertion, dedupe)
- `my-site/public/ppi/index.html` (frontend normalization, suppression, mode visibility)

Snapshot date from code state: 2026-03-06.

## 1) End-to-end Pipeline

### 1.1 Structure ingest and atom parsing

Backend accepts PDB/mmCIF text (or fetches mmCIF via PDBe model server when needed), then parses heavy atoms.

Atom-level ingest filters:

- Hydrogens are dropped.
- mmCIF loop tokenization is CIF-safe for nucleic atom names with apostrophes (for example `O5'`), so chain/sequence/atom fields are not misaligned during parsing.
- Altloc is filtered at parse time:
  - PDB: keep only blank, `A`, or `1`.
  - mmCIF: keep only `.`, `?`, `A`, or `1`.
- Chains are normalized via `ChainAliases` (label/auth remapping).

### 1.2 Contact candidate generation

`analyze_interface()` does the following:

1. Build residue atom index keyed by normalized `(chain, seq)`.
2. Build Arpeggio selection:
   - If focused residue is provided and valid, select that residue only.
   - Else if inter-chain, select chain-A residues within `8.0 A` of chain B.
   - Else fallback to chain-wide selection.
3. Run PDBe Arpeggio with:
   - `interacting_cutoff = 6.0 A`
   - `vdw_comp = 0.1`
   - sequence-adjacent inclusion disabled when supported.
4. Keep only contacts matching requested chain pair.
5. If residue focus is active, keep only contacts touching focus residue.

### 1.3 Pre-assertion duplicate suppression (raw candidate stage)

Before classification, candidates are deduped by a normalized signature combining:

- residue pair,
- atom pair,
- atom feature-pair signature (element + donor/acceptor/hydrophobic/backbone flags),
- Arpeggio interaction type,
- Arpeggio term motif,
- symmetry class (identity vs non-identity),
- altloc compatibility class,
- ring-site token,
- inferred donor/acceptor direction (`A->B`, `B->A`, both, or none).

If duplicates exist, shortest-distance contact wins (then earliest raw order tiebreak).

### 1.4 Residue-level base-pair prepass

`_compute_base_pair_pair_stats()` computes residue-pair support metrics before per-contact assertion.

A contact contributes only when all are true:

- both residues are nucleic bases,
- residues are not the same residue identity (same chain + sequence),
- residues are not sequence-adjacent,
- both atoms are nucleobase atoms,
- both are pairing-edge atoms,
- neither is glycosidic atom,
- both elements are polar (`N/O/S/SE`),
- donor/acceptor complementarity holds,
- distance `<= 3.6 A`,
- and the contact does not fail impossible-contact preclassification.

Computed per residue-pair statistics include support count, best distance, angle pass counts, canonical-template matches, coplanarity support, mutual-best partner, and score components.

### 1.5 Per-contact authoritative assertion

Every remaining contact goes through `_assert_interaction()` (the final family assignment).

Global precedence in code is:

1. Identity/conformer/model/symmetry gate (invalid-contact gate).
2. Covalent/near-covalent nonbonded artifact exclusion.
3. Impossible-contact preclassification gate (invalid/clash).
4. Same-residue nonbonded suppression (unless explicitly metal/aromatic-context exception).
5. Family assertions and demotions:
   - adjacent nucleotide linkage suppression,
   - metal coordination,
   - nucleic backbone O/P suppression,
   - salt bridge (evaluated before hbond; if salt assertion fails, candidate may continue to later directional halogen/hbond/polar evaluation),
   - halogen bond,
   - pi/aromatic,
   - hydrogen bond,
   - base-pairing branch (residue-aware; evaluated after directional atom-pair families),
   - nucleobase polar fallbacks,
   - hydrophobic / packing,
   - VDW clash branch,
   - VDW/proximal fallback,
   - `other` fallback.

### 1.6 Post-assertion backend processing

After assertion:

- Drop `excludeFromNoncovalent` and `covalent_bond`.
- Drop internal `invalid_contact` (hidden from API output).
- Apply atom/element/distance overrides from assertion.
- Build normalized `record` payload with asserted evidence.
- Bucket-specific dedupe/caps:
  - hydrophobic: one per residue pair (best distance/specificity wins),
  - salt: one per residue pair,
  - metal: dedupe by `(metal atom, donor residue)`; ambiguous donor atoms annotated,
  - hbond carboxylate oxygen ambiguity dedupe,
  - aromatic families: ring-site dedupe + top-k caps.
- Apply atom-reuse confidence penalties for crowded families:
  - hydrogen bonds: threshold 3,
  - halogen bonds: threshold 1,
  - hydrophobic: threshold 4.
- Apply optional API mode filter.

`meta.classifier` is now `plausibility+assertion:v2`.

## 2) Universal Gates (Applied Before Family Assignment)

### 2.1 Identity/conformer/model/symmetry gate

Classifies as internal `invalid_contact` if any of these hold:

- same normalized atom signature on both endpoints,
- incompatible altlocs,
- same residue + same atom + same altloc family,
- different model IDs,
- unresolved duplicate mapping after chain alias normalization,
- non-identity symmetry-generated contact (crystal contact), when `ALLOW_CRYSTAL_CONTACTS = False`.

Result is hidden from normal output (`debugOnly`, dropped before buckets).

### 2.2 Covalent-neighbor artifact gate

Classifies as `covalent_bond` (excluded from noncovalent output) for:

- explicit covalent terms,
- likely P–O covalent neighbors in nucleotide-like context at `<= 1.9 A`.

### 2.3 Same-residue nonbonded suppression gate

By default, same-residue nonbonded atom-pairs are filtered as internal `invalid_contact`.
These same-residue pairs are also excluded from residue-level base-pair support prepass.

Exceptions are intentionally narrow:

- explicit/likely metal-coordination context,
- explicit/likely aromatic/pi intramolecular context.

### 2.4 Impossible-contact hard-stop gate

Uses heavy-atom distance and VdW overlap before family logic:

- `overlap = vdw(A) + vdw(B) - d`
- pair-specific minimum nonbonded distances are enforced,
- absolute invalid floor at `0.8 A`,
- preclassification hard clash when overlap `>= 0.7 A`.

Outcomes:

- `invalid_contact` for impossible mapping-like distances,
- `clash` for hard steric impossibility.

This prevents impossible pairs from being rescued into hbond/base-pair/hydrophobic.
Operational split used in practice: `invalid_contact` is used when identity/mapping/conformer duplication is implicated, while `clash` is used when atoms are otherwise legitimate but sterically impossible.

## 3) Assertion Rules by Family

### 3.1 Base pairing (`base_pairing`)

Residue-pair-aware nucleobase branch with atom-level checks.
It is evaluated after salt/halogen/hbond atom-pair branches so directional atom chemistry has precedence.

Asserts when either:

- canonical Watson-Crick atom template + mutual-best partner + not sequence-adjacent,
- or multi-polar support (`>=2`) + mutual-best partner + not sequence-adjacent.

Otherwise demotes to `polar_contact` / debug `polar_proximal` with reasons.

### 3.2 Metal coordination (`metal_coordination`)

Requires metal endpoint + donor element in `{O,N,S,SE}` + distance within metal-specific cutoff.

Confidence:

- high with explicit metal term,
- otherwise medium.

### 3.3 Salt bridge (`salt_bridge`)

Triggered by ionic context (unless preempted by metal coordination), with endpoint reassignment to best charged sites.
Salt evaluation occurs before hbond evaluation and wins when ionic topology is unambiguous.

Requires:

- valid cation/anion topology,
- unambiguous cation and anion side,
- distance `<= 4.8 A`.

If salt constraints fail:

- fallback is `polar_contact` when no directional halogen/hbond/polar path is available,
- otherwise control falls through to later directional halogen/hbond/polar evaluation.

### 3.4 Halogen bond (`halogen_bond`)

Dedicated branch before hydrophobic:

- donor element must be one of `{Cl, Br, I}`,
- acceptor element must be one of `{O, N, S, Se}`,
- explicit `X···C` acceptor is forbidden,
- donor must be covalently bound to carbon (local geometry from coordinates),
- distance within donor-element cutoff,
- overlap veto: halogen assertion is blocked when VdW overlap exceeds soft-overlap/hydrophobic-overlap guardrails,
- directional angle threshold:
  - strong `>= 155°`,
  - medium `>= 145°`.

If explicit halogen term exists but constraints fail, demotion is failure-mode specific:

- chemically valid donor/acceptor + acceptable distance + weak angle -> `polar_contact`,
- distance too long / poor geometry -> debug `proximal`,
- overlap above halogen overlap guard -> `clash`,
- invalid donor/acceptor identity -> `other` or debug `proximal`.

### 3.5 Pi / aromatic (`pi_pi`, `pi_cation`, `aromatic_packing`, `aromatic_proximal`)

Uses reported or recomputed ring geometry:

- `pi_pi` requires centroid/interplanar/lateral/angle geometry support,
- `pi_cation` requires cation-ring geometry support,
- if strict pi geometry fails, may reclassify to `aromatic_packing`,
- then to debug `aromatic_proximal` if only weak proximal geometry remains.

### 3.6 Hydrogen bond and polar (`hbond`, `polar_contact`, `polar_proximal`)

Hbond requires donor/acceptor validity + strict distance constraints; angle and proxy-angle logic controls confidence and demotion.

Distance floor behavior:

- `1.45 A` remains a global corrupt-data floor,
- pair-specific nonbonded minima are applied before chemistry and are also respected by strict H-bond distance checks.

Special NA backbone tightening:

- phosphate/sugar oxygen pairs are blocked from hbond assertion,
- phosphate/sugar oxygen pairs demote to `polar_contact` at short polar-compatible distance and otherwise to debug `polar_proximal`,
- broader phosphate-backbone oxygen neighborhood pairs (including phosphate-phosphate O···O) demote to debug `polar_proximal` by default.

Extreme-short suspicious contacts are demoted or converted to clash based on overlap and angle evidence.

### 3.7 Hydrophobic and packing (`hydrophobic`, `packing_contact`, `proximal`)

Hydrophobic assertion requires:

- nonpolar atom eligibility,
- pair-specific minimum distance (`max(pair-min, vdw-sum-0.5, global minimum)`),
- distance within hydrophobic max (`<= 4.6 A`),
- VdW overlap not exceeding `0.5 A`,
- no stronger directional chemistry context (hbond/halogen-like).

If hydrophobic-like but not assertable, demotes through `packing_contact`, debug `proximal`, or `other`.

### 3.8 Clash and generic fallback (`clash`, `packing_contact`, `proximal`, `other`)

After earlier branches:

- explicit `VDW_CLASH` terms can still produce clash if overlap thresholds are exceeded,
- otherwise VDW/PROXIMAL terms map to packing/proximal/other according to chemistry and distance.

If nothing matches: `other`.

## 4) Thresholds and Constants

Key numeric values currently used:

- Polar contact max distance: `3.8 A`.
- Hbond explicit max: `3.7 A`.
- Hbond candidate max: `3.6 A`.
- Hbond strong angle: `150°`.
- Hbond heavy-atom minimum: `1.45 A`.
- Salt bridge max/confident: `4.8 / 4.2 A`.
- Hydrophobic max: `4.6 A`.
- Hydrophobic overlap max: `0.5 A`.
- Hard clash preclassification overlap: `0.7 A`.
- Soft clash preclassification overlap flag: `0.4 A`.
- Halogen overlap assertion veto uses soft/hydrophobic overlap guard (`min(0.4, 0.5) = 0.4 A` effective cap).
- Invalid absolute nonbonded minimum: `0.8 A`.
- Pair-specific nonbonded minima include:
  - `N···N >= 2.4 A`,
  - `N···O >= 2.3 A`,
  - `O···O >= 2.4 A`,
  - `C···C >= 2.8 A`,
  - `C···N >= 2.7 A`,
  - `C···O >= 2.7 A`,
  - `C···Cl >= 3.0 A`,
  - `C···S >= 2.9 A`,
  - `Cl···Cl >= 3.3 A`.
- For element pairs not explicitly listed above, preclassification falls back to overlap-driven clash/invalid gating plus the global absolute minimum floor (`0.8 A`), rather than a complete hardcoded pair table.
- Halogen donor distance limits:
  - `Cl: 3.5 A`,
  - `Br: 3.7 A`,
  - `I: 3.9 A`.
- Halogen angle thresholds:
  - strong `>= 155°`,
  - medium `>= 145°`.

## 5) Backend Output Buckets and Hidden Types

Public contact buckets:

- `hydrogen_bonds`
- `polar_contacts`
- `base_pairing`
- `salt_bridges`
- `halogen_bonds`
- `hydrophobic`
- `metal_coordination`
- `pi_pi`
- `pi_cation`
- `aromatic_packing`
- `other`

Internal/non-public families:

- `invalid_contact` (debug/internal; removed before output)
- `covalent_bond` (excluded)
- debug-style families usually routed to `other` bucket (`proximal`, `polar_proximal`, `aromatic_proximal`, `packing_contact`, etc.).

## 6) Backend Mode Filter (`filter_contacts_by_mode`)

Accepted mode aliases:

- hydrophobic mode: `hydrophobic`
- electrostatic mode: `electrostatic|ionic|salt` -> `salt_bridges`
- polar mode: `polar|polar_contact|polar_contacts` -> `polar_contacts + halogen_bonds`
- base-pair mode: `base_pair|base_pairs|base_pairing`
- metal mode: `metal|metal_coordination|coordination`
- hbond mode: `hbond|hbond_network|hydrogen`
- halogen mode: `halogen|halogen_bond|halogen_bonds|xbond`
- aromatic mode: `aromatic|pi` -> `pi_pi + pi_cation + aromatic_packing`
- other mode: `other`

## 7) Frontend Normalization and Filters

Frontend applies additional normalization and suppression before rendering.

### 7.1 Category normalization

Declared/asserted category tokens are normalized. `halogen_bond` is recognized as its own interaction mode (`halogen`) and is also included by Polar aggregate mode queries.

### 7.2 Contact source normalization and inference

Frontend can:

- normalize malformed bucket keys,
- split/rebucket contacts by category,
- infer aromatic non-polymer contacts when needed,
- for residue-focused analysis, start from nearby chain partners but broaden to full structure-chain partners when only self-chain is found (or when focused residue is non-polymer), so valid cross-chain contacts are not missed,
- canonicalize common modified amino-acid and nucleotide residue names for polymer typing (for example phospho-residues and common modified bases), with residue-atom signature fallback so ribbon/backbone handling remains consistent for modified polymers,
- classify polymer-vs-nonpolymer and backbone-vs-sidechain in rendering with residue-key-aware polymer inference (name aliases first, atom-signature fallback second), so modified residues do not leak backbone atoms into sidechain views and modified nucleotides remain ribbon-compatible,
- dedupe within mode-specific grouping keys.

### 7.3 Global suppression by precedence

`shouldSuppressContactByPrecedence()` suppresses weaker overlaps when stronger families exist for same residue pair, including aromatic-over-hydrophobic suppression in ring contexts.
Halogen precedence is explicit:

- `halogen_bond` suppresses weaker `hydrophobic` and `packing_contact` in the same atom-pair or residue-pair context,
- `halogen_bond` is not suppressed by generic `polar_contact`.
Suppression order for halogen precedence is deterministic: exact atom-pair match is checked first, then residue-pair fallback is applied.

### 7.4 Debug-only suppression

`debugOnly` contacts are hidden unless interaction debug mode is enabled.

### 7.5 Visibility and anatomy filters

Final visibility requires:

- mode visibility enabled,
- not precedence-suppressed,
- passes anatomy isolation (`sidechain` or `backbone`) when active,
- passes debug-mode gate,
- sidechain/backbone anatomy counts and submenu items are computed from the same visible interaction-mode set; `other` is excluded unless debug mode is enabled,
- excludes `other`-mode contacts from focused animations when debug mode is off,
- for focused interaction animations, passes active panel-mode gating (`summary` shows all visible families; specific panels animate only contacts that belong to that panel’s displayed set, with `hbond` panel also including base-pair subitems).

### 7.6 UI mode mapping (important)

Rendered interaction modes are:

- Hydrophobic
- H-bond
- Polar (includes `polar_contact` and also includes halogen contacts in aggregated mode queries)
- Halogen (dedicated)
- Base pairing
- Salt bridge
- Clash (hidden by default via mode visibility; still available explicitly)
- Metal coordination
- Aromatic
- Other

### 7.7 PAE Hover Overlay Behavior

When a residue is hovered/focused from the 3D model and mapped onto the PAE matrix, the panel draws:

- a diagonal-centered square marker at `(i, i)` for each mapped axis index,
- dashed horizontal and vertical guide lines extending from that square toward panel edges.

This replaces the previous full-width/full-height trace rectangles.

## 8) What is Specifically Hidden or Demoted

Explicitly hidden from normal output:

- identity/conformer/model/symmetry artifacts (`invalid_contact`),
- covalent-neighbor artifacts (`covalent_bond`),
- debug-only proximal-style contacts unless debug mode is enabled.

Systematic demotions happen when constraints fail:

- halogen -> `polar_contact` / `proximal` / `clash` / `other` (failure-mode dependent),
- hbond -> polar,
- base-pair candidate -> polar/polar_proximal,
- hydrophobic -> packing/proximal/other,
- aromatic strict family -> aromatic_proximal/other.

## 9) Practical Consequences of Current Rules

- Impossible self/duplicate/altloc/model/symmetry contacts are blocked before chemistry assignment.
- Hard overlaps are globally preclassified as clash (or invalid if mapping-like).
- Base-pair inference remains residue-pair-aware, and now runs after directional atom-pair families.
- Halogen bonds now have a dedicated assertion path (no longer forced to hydrophobic fallback).
- Short nonpolar artifacts are strongly reduced by pair-specific minima plus overlap gating.
