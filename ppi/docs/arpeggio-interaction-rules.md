# Arpeggio Interaction Rules and Filters (Current Full Spec)

This document is the implementation-level specification for interaction classification and filtering in:

- `my-site/ppi/api/analysis.py` (backend generation, assertion, dedupe)
- `my-site/public/ppi/index.html` (frontend normalization, suppression, mode visibility)

Snapshot date from code state: 2026-09-20.

The thresholds below are Roami's screening rules, not universal definitions of chemical bonds. Support levels describe the available chemical and geometric evidence; they are not calibrated probabilities, experimental confirmation or estimates of interaction energy. Contact counts and schematic visual effects do not measure binding affinity.

## 1) End-to-end Pipeline

### 1.1 Structure ingest and atom parsing

Backend accepts PDB/mmCIF text (or fetches mmCIF via PDBe model server when needed), then parses heavy atoms.

Atom-level ingest filters:

- The heavy-atom coordinate index excludes H/D. Supplied hydrogen coordinates remain available through the Arpeggio input and exported H-geometry metadata; excluding them from the heavy-atom index does not mean they are discarded as evidence.
- mmCIF loop tokenization is CIF-safe for nucleic atom names with apostrophes (for example `O5'`), so chain/sequence/atom fields are not misaligned during parsing.
- One nonblank conformer is selected per residue by highest mean occupancy, retaining shared blank-altloc atoms. Ties prefer `A`, then `1`, then lexical order. See coordinate and evidence consistency below for duplicate-atom and model handling.
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

### 1.3 Raw candidate retention

Each source observation is retained until classification has supplied chemical roles, physical site identities and hydrogen identity. No residue-pair or raw atom-pair cap removes an observation before that information exists. Equivalent canonical interactions are merged after assertion, preserving their provenance.

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
   - nucleobase pairing-context screening, retaining `hbond` only for compatible atom roles/geometry and otherwise a polar family,
   - adjacent nucleotide linkage suppression,
   - metal coordination,
   - nucleic backbone O/P suppression,
   - salt bridge (evaluated before hbond; if salt assertion fails, candidate may continue to later directional halogen/hbond/polar evaluation),
   - halogen bond,
   - pi/aromatic,
   - hydrogen bond,
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
- Build a normalized `record` with versioned `semantics`, including roles, physical sites, direction, typed geometry and evidence.
- Merge only matching canonical identities within a family. Identity includes role/site IDs, model/conformer identity and an actual hydrogen ID where available. Distinct hydrophobic atom contacts, charged features, carboxylate acceptors and physical rings remain distinct.
- Retain merged source observations and role/typing provenance; charge-state disagreement is explicit ambiguity.
- Preserve separate metal–donor edges and annotate `denticity`/`coordinationDonorAtoms` for simultaneous chelating atoms.
- Apply atom-reuse penalties to legacy `asserted.confidence` for crowded families; these do not replace `semantics.evidence`:
  - hydrogen bonds: threshold 3,
  - halogen bonds: threshold 1,
  - hydrophobic: threshold 4.
- Apply optional API mode filter.

`meta.classifier` is `plausibility+assertion:v2`, `analysisVersion` is `roami-assertion-2.4`, and `meta.interactionSemanticsVersion` is `1`. The shared role/evidence contract is documented in [interaction-semantics.md](interaction-semantics.md).

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
- `clash` for excessive overlap or a pair-specific nonbonded-distance violation in the model.

This prevents these flagged pairs from being reassigned as hbond/base-pair/hydrophobic. Supported metal-coordination contacts bypass this overlap preclassification because ordinary nonbonded radii are not appropriate for their coordination distances.
Operationally, `invalid_contact` identifies mapping-like or extreme-distance artifacts, while `clash` identifies suspicious nonbonded geometry. Radii are approximate: a clash flag is not proof of literal steric impossibility or of its cause. Roami's count is not the all-atom MolProbity clashscore.

## 3) Assertion Rules by Family

### 3.1 Base-pair context (`basePair`)

Residue-pair-aware nucleobase screening retains pairing context separately from atom-level chemistry. The early nucleobase branch emits `hbond` only when actual endpoint donor/acceptor roles and H-bond geometry/distance criteria are compatible; otherwise it emits a polar contact. It does not emit an atom-level `base_pairing` family.

Pairing-context support comes from either:

- canonical Watson-Crick atom template + mutual-best partner + not sequence-adjacent,
- or multi-polar support (`>=2`) + an explicit upstream H-bond term or a passing donor–H–acceptor angle + mutual-best partner + not sequence-adjacent.

Otherwise the atom contact can remain `polar_contact` / debug `polar_proximal` with reasons. The second branch requires both residue-pair support and the stated directional evidence. Neither context route overrides incompatible actual atom typing. When an engine role conflicts with a neutral nucleobase template, `roleAlternatives` preserves both sources, the affected role is unresolved, and the assignment/direction is ambiguous. These heuristics are not a complete Leontis–Westhof annotation and can miss unusual pairs or triples. `basePair` context is not a second atom contact or evidence for an H-bond by itself.

### 3.2 Metal coordination (`metal_coordination`)

Requires metal endpoint + donor element in `{O,N,S,SE}` + distance within metal-specific cutoff.
This element/distance screen identifies coordination candidates. It does not establish donor lone-pair availability, metal identity or oxidation state, coordination number, or the geometry of the complete shell.

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

Uses reported or recomputed ring geometry. In the main pi branch:

- Stacked `pi_pi` requires centroid separation 3.3–6.2 Å, interplanar separation 2.6–4.3 Å, lateral offset ≤3.2 Å and normal angle ≤30°.
- T-shaped `pi_pi` requires normal angle ≥60° and centroid separation 3.3–6.2 Å; lateral offset must be ≤4.0 Å when present. This branch does not enforce the stacked interplanar-distance window.
- Geometry-supported `pi_cation` uses the 3.3–6.2 Å centroid-distance window and lateral offset ≤4.0 Å when present.
- Incomplete ring geometry can retain a pi family at low support in some branches; a pi category therefore does not guarantee that every geometric check was performed. Complete geometry that fails pi criteria may instead yield `aromatic_packing`, then debug `aromatic_proximal` for broader proximity.

Additional aromatic-context branches can promote contacts with supporting ring geometry without assigning a stacked/T-shaped subtype. Inspect the subtype, evidence and ring metrics rather than treating every pi-family row as the same fully validated geometry.

### 3.6 Hydrogen bond and polar (`hbond`, `polar_contact`, `polar_proximal`)

H-bond classification requires complementary donor/acceptor roles and strict distance constraints. Arpeggio/Open Babel atom types carry ligand valence/protonation information; an unknown ligand nitrogen or oxygen is not assigned a role from its element alone. Explicit upstream H-bond evidence is retained when an older contact export lacks role metadata.

An exported donor–H–acceptor angle can support geometry; the report identifies whether hydrogens came from the input or were generated by Arpeggio. Generated hydrogen geometry is a modeled hypothesis, not experimental confirmation. `hbond_confirmed` is reserved for a supplied H-angle passing the strong angle threshold; generated or completed hydrogen positions remain `hbond_candidate`, at most medium confidence. Explicit Arpeggio H-bonds below that strong threshold remain candidates when their exported angle satisfies the engine criterion. A hydroxyl heavy-atom axis neither validates nor rejects an H-bond: hydroxyl rotation and its bent bond geometry make that axis an invalid substitute for the hydrogen position. Without hydrogen geometry, distance/role candidates remain at most medium confidence.

Distance floor behavior:

- `1.45 A` is the H-bond heavy-atom floor before any stricter pair-specific minimum; the separate universal invalid-contact floor is `0.8 A`,
- pair-specific nonbonded minima are applied before chemistry and are also respected by strict H-bond distance checks.

Special NA backbone tightening:

- chemically donor-capable RNA hydroxyls (including 2′-OH) can hydrogen-bond to phosphate acceptors with the same evidence requirements as other H-bonds,
- acceptor-only phosphate/ether oxygen pairs demote to `polar_contact` at short polar-compatible distance and otherwise to debug `polar_proximal`,
- broader acceptor-only phosphate-backbone oxygen neighborhoods (including phosphate-phosphate O···O) demote to debug `polar_proximal`,
- actual covalent and adjacent linkage contacts remain excluded.

Extreme-short suspicious contacts are demoted or converted to clash based on overlap and angle evidence.

### 3.7 Hydrophobic and packing (`hydrophobic`, `packing_contact`, `proximal`)

Hydrophobic assertion requires:

- nonpolar atom eligibility,
- pair-specific minimum distance (`max(pair-min, vdw-sum-0.5, global minimum)`),
- distance within hydrophobic max (`<= 4.6 A`),
- VdW overlap not exceeding `0.5 A`,
- no higher-priority directional chemistry context (hbond/halogen-like).

Classification priority is not an interaction-energy ranking. Nonpolar proximity alone does not measure solvent burial or the contribution of the hydrophobic effect.

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
- `base_pairing` (legacy compatibility bucket; current pairing context remains on atom-level contacts)
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

Frontend normalizes canonical records without reclassifying their chemical roles or pruning distinct identities. Legacy records without versioned semantics still use the compatibility inference/suppression paths described below.

### 7.1 Category normalization

`semantics.family` is authoritative for canonical records. Legacy declared/asserted category tokens are normalized with `asserted.family` taking precedence when present. `halogen_bond` is recognized as its own interaction mode (`halogen`) and is also included by Polar aggregate mode queries. A supplied `clash` bucket is retained; the backend may also deliver clash-classified records through `other`.

### 7.2 Contact source normalization and inference

For all reports, frontend normalization preserves record metadata, normalizes category/bucket names and applies presentation grouping. Canonical identities are deduplicated with their source observations retained; distinct chemical edges survive. Chain aliases remap nested participant residues without rebuilding chemical roles or opaque source identities. Canonical reports, including empty reports carrying `interactionSemanticsVersion: 1`, bypass local aromatic, nucleobase and intra-residue metal inference.

Legacy compatibility and structure-display paths can:

- normalize malformed bucket keys,
- split/rebucket contacts by category,
- infer aromatic non-polymer contacts when needed,
- for residue-focused analysis, start from nearby chain partners but broaden to full structure-chain partners when only self-chain is found (or when focused residue is non-polymer), so valid cross-chain contacts are not missed,
- canonicalize common modified amino-acid and nucleotide residue names for polymer typing (for example phospho-residues and common modified bases), with residue-atom signature fallback so ribbon/backbone handling remains consistent for modified polymers,
- classify polymer-vs-nonpolymer and backbone-vs-sidechain in rendering with residue-key-aware polymer inference (name aliases first, atom-signature fallback second), so modified residues do not leak backbone atoms into sidechain views and modified nucleotides remain ribbon-compatible,
- dedupe within mode-specific grouping keys.

### 7.3 Global suppression by precedence

For canonical contacts, `shouldSuppressContactByPrecedence()` returns false: the producer already applied family precedence, and a second heuristic must not remove a distinct chemical edge. For legacy contacts it removes weaker duplicate explanations from the displayed set. Generic packing uses stronger-family context; aromatic-over-hydrophobic suppression requires overlapping ring atoms, while non-ring aliphatic contacts remain independent.
Halogen precedence is explicit:

- `halogen_bond` suppresses weaker `hydrophobic` and `packing_contact` records for the same unordered atom pair,
- independent, fully identified atom pairs between those residues remain visible,
- residue-pair fallback applies only when the weaker record lacks an atom name on at least one endpoint,
- `halogen_bond` is not suppressed by generic `polar_contact`.

Exact atom-pair matching is checked first and works in either endpoint order. A failed exact match between fully identified endpoints does not trigger residue-wide suppression. These precedence filters are disabled in interaction debug mode.

### 7.4 Debug-only suppression

`debugOnly` contacts are hidden unless interaction debug mode is enabled. Normal `/analyze` requests use `includeDiagnostics: false`; the smaller response can omit diagnostic records. The canonical report remains available from `GET /report/{report_id}`, and `includeDiagnostics: true` requests the diagnostic variant directly.

The frontend caches compact and diagnostic variants separately. `?debug=true` enables diagnostics at startup. For a live transition, use `await window.setInteractionDebugMode(true)` or `false`; this reloads the existing chain or residue analysis scope. Assigning the legacy `window.__PPI_INTERACTION_DEBUG_MODE` flag alone cannot recover records absent from an already loaded compact response.

### 7.5 Visibility and anatomy filters

Final visibility requires:

- mode visibility enabled,
- not precedence-suppressed,
- passes anatomy isolation (`sidechain` or `backbone`) when active,
- passes debug-mode gate,
- sidechain/backbone anatomy counts and submenu items are computed from the same visible interaction-mode set; `other` is excluded unless debug mode is enabled,
- excludes `other`-mode contacts from focused animations when debug mode is off,
- for focused interaction animations, passes active panel-mode gating (`summary` shows all visible families; specific panels animate only contacts that belong to that panel’s displayed set, with `hbond` panel also including base-pair subitems).

Backbone anatomy requires a polymer endpoint with a backbone atom name. A water `O`, ligand `N`, or calcium `CA` cannot turn a protein side-chain contact into a backbone contact merely by sharing a familiar atom name. Protonated histidine may serve as the cation in a cation–pi interaction against an aromatic partner; neutral histidine is not assumed to have that role.

### 7.6 UI mode mapping (important)

Rendered interaction modes are:

- Hydrophobic
- H-bond
- Polar (includes `polar_contact` and also includes halogen contacts in aggregated mode queries)
- Halogen (dedicated)
- Base-pair context
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

### 7.8 Shared category and display contract

Backend `semantics` supplies the authoritative family, roles, physical sites, direction, measurement definitions and evidence; legacy fields remain compatibility data. API family-mode filtering projects the corresponding buckets; the frontend requests `mode: 'all'` and derives aggregate views from those records and their metadata. Frontend normalization can deduplicate equivalent identities and group records for display, but a rendering budget must not delete distinct chemical edges from the canonical contact set.

The UI modes overlap deliberately:

- Polar includes halogen contacts, which also appear in the dedicated Halogen view.
- `basePair` metadata records residue-pair context independently of the atom-level family. An H-bond carrying that metadata appears in both H-bond and Base-pair context views and keeps its H-bond label. A polar contact carrying context remains polar, is excluded from the H-bond view, and contributes only once to unique-contact totals.
- Sidechain and Backbone describe anatomy across interaction families; they are not additional chemical families.

Mode totals therefore must not be added together as a count of unique interactions. Counts, list membership, anatomy isolation, and active-panel animation gating use the same category queries. Their index is reused for an unchanged report and rebuilt when the structure, debug mode, or contact buckets change. Expanding a collapsed residue group creates its rows on demand without changing its contact count.

`report.displayGroups` adds family-specific chemical/display units without deleting the canonical evidence. Dense hydrophobic atom contacts are grouped by nonpolar features and local spatial patches; metal–donor, H-bond, halogen and clash assignments retain their atom-specific units. Supporting atom contacts, source observations, canonical assignments, available display objects and actually rendered connectors have separate definitions. See [contact-aggregation.md](contact-aggregation.md) for the grouping contract, 3.5 Å hydrophobic endpoint-patch policy, diagnostic scope and missing-support handling. Group/atom counts and minimum/mean/median distances do not measure independent bond energies or strength.

### 7.9 Ring rendering and support labels

Canonical ring participants carry the exact ring membership, centroid and normal supplied by the backend; rendered endpoints consume those sites. Frontend descriptors for legacy records use separate physical five- and six-membered rings for fused Trp/purine systems, with exported `ringAtomNamesA/B` taking precedence. Older reports identifying a whole aromatic system retain that exact membership rather than silently substituting a smaller ring.

Cation–pi context overlays display at most two connectors per identified ring for legibility. This cap is applied by `selectContactsForVisualBudget()` during drawing, not during report normalization. A specifically selected contact is prioritized, all distinct contacts remain listed, and debug mode bypasses the connector cap. The overlay records its omitted connector count in `group.userData.hiddenInteractionCount`; aromatic help text explains the cap.

Canonical labels display `candidate`, `chemically_supported`, `geometrically_supported` or `ambiguous`, with direction and provenance supplied by the record. Legacy high/medium/low labels and `hbond_confirmed` tokens remain compatibility data, not calibrated probabilities or experimental confirmation. Contact-count rankings are labeled “Contact-rich residues”; they do not establish energetic binding hotspots. Ribbon bridge positions are schematic surface projections, whereas reported distances retain their declared atom- or ring-based definitions.

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
- Hard overlaps are preclassified as clash (or invalid if mapping-like), with the supported-metal-coordination exception described above.
- Base-pair inference remains residue-pair-aware, and now runs after directional atom-pair families.
- Halogen bonds now have a dedicated assertion path (no longer forced to hydrophobic fallback).
- Short nonpolar artifacts are strongly reduced by pair-specific minima plus overlap gating.


## Coordinate and evidence consistency (20 September 2026)

The coordinate index supplies authoritative element identity before assertion and deduplication. This preserves two-letter elements such as Fe, Cl and Br when Arpeggio exports only atom names. Direct API analysis selects the first coordinate model. Within each residue, one nonblank conformer is selected by highest mean atom occupancy (missing occupancy defaults to 1); ties prefer A, then 1, then lexical order. Blank shared atoms are retained, and duplicate names prefer blank then higher occupancy then input order. The temporary Arpeggio input applies the same policy. Alternate labels on different residues are not intrinsically incompatible.

`basePair` records residue-pair context independently of the atom-level family, so an additional H-bond label cannot erase it. Noncanonical atom-edge assertions require directional evidence in addition to residue support; nearby diagonal contacts do not acquire a bond merely because their residues form a pair.

Ring descriptors are cached per analysis, and `ringAtomNamesA/B` communicate the exact selected physical ring membership for drawing. Fused Trp/purine systems have separate five- and six-membered descriptors; whole-base plane tests continue to use the base atom set.

`buriedFraction` is null because solvent-accessible surface burial is not calculated. `contactingResidueFraction` names the former count ratio accurately. Inline coordinates determine cache identity even when a PDB accession is also present.

For partially hydrogenated inputs, Open Babel perceives omitted hydrogen valence at the supplied formal-charge/protonation state. Every original atom ID, coordinate and charge must remain unchanged. Missing hydrogens stay implicit: their underconstrained orientations are not generated. These donor-role candidates carry `implicit_hydrogen_valence` evidence and require actual supplied H geometry for angle validation. This avoids treating the upstream any-hydrogen shortcut as evidence that all other atoms lack donor hydrogens, while keeping repeated analyses deterministic.
