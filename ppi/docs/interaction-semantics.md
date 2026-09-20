# Canonical interaction roles and evidence

Each current contact carries a versioned `semantics` object. This is the chemical contract consumed by normalization, deduplication, the API, the UI, explanations and rendering. `meta.interactionSemanticsVersion: 1` also identifies the contract on empty reports. Legacy `residueA`, `residueB`, `distance`, `asserted` and source-engine fields remain for compatibility; they must not replace the semantic roles or redefine a distance metric.

The contract separates five questions:

1. **Model geometry:** what distance, angle, centroid or overlap was calculated from the supplied coordinates?
2. **Chemical compatibility:** what typing supports the proposed donor, acceptor, charge or aromatic feature?
3. **Chemical state:** was a protonation/tautomer/charge state supplied, modeled, assumed or left unknown?
4. **Assignment:** which interaction family and categorical evidence level describe this record?
5. **Interpretation:** has biological importance or an energetic contribution actually been evaluated? The current answer is `not_evaluated`.

Here “observed geometry” means geometry present in the coordinate model. Coordinates can themselves be predictions, and a hydrogen added by a chemistry engine is a model rather than experimental evidence.

## Record shape

```text
contact.semantics = {
  version: 1,
  identity: string,
  family: string,
  directionality: "directional" | "role_specific" | "symmetric",
  participants: [
    {
      side: "A" | "B",
      role: string,
      site: {
        kind: "atom" | "ring" | "group" | "metal",
        id: string,
        residue: { chain, seq, resName },
        atoms: [{ id, atomName, element, coordinates, modelId, altloc }],
        contactAtom?: { id, atomName, element, coordinates, modelId, altloc },
        coordinationSiteId?: string,
        centroid?: [x, y, z],
        normal?: [x, y, z],
        provenance?: string
      },
      roleProvenance: [string, ...],
      roleAlternatives?: [{ role: string, source: string }, ...],
      atomTypes?: [string, ...],
      charge?: { value: number | null, sign?: -1 | 1, source: string, inferred: boolean,
                 scope?: "group", atomValue?: number }
    }, ...
  ],
  direction: {
    from: "A" | "B" | null,
    to: "A" | "B" | null,
    certainty: "certain" | "inferred" | "ambiguous" | "not_applicable",
    alternatives?: [...]
  },
  geometry: {
    distance: { value, kind, unit: "angstrom", source },
    measurements: [{ kind, value, unit, source }, ...],
    hydrogen?: { id, atomName, element, coordinates, modelId, altloc, source, inferred, bondSource? },
    donorAnchor?: { id, atomName, element, coordinates, modelId, altloc, source }
  },
  evidence: {
    level: "candidate" | "chemically_supported" | "geometrically_supported" | "ambiguous",
    chemicalCompatibility: "supported" | "inferred" | "unknown",
    geometrySupport: "supported" | "partial" | "missing",
    chemicalState: "input" | "modelled" | "assumed" | "unknown",
    ambiguityFlags: [string, ...],
    biologicalInterpretation: "not_evaluated",
    supportingRecords?: number
  },
  sourceEndpoints?: [{ side, atom, atomTypes, charge }, ...],
  sourceObservations?: [{ atomPair, source, terms, geometry, participants, direction, evidence, sourceEndpoints? }, ...]
}
```

Absent measurements and unknown state remain absent or null. They must not become zero. Conversely, a supplied formal charge of zero is a real value and must not be treated as missing. Numeric zero is not a substitute for an unknown metal oxidation state.

`side` ties a participant to the record's A/B endpoint ordering. It is not a chemical role. Reversing storage order must preserve site identity, assigned roles, direction, hydrogen identity and geometric meaning together.

## Roles and notation

| Family | Semantics | Participants | Display relationship |
| --- | --- | --- | --- |
| Hydrogen bond | Directional | `donor`, `acceptor`; unresolved alternatives can use `donor_or_acceptor` or `unresolved` | Donor → acceptor only for a resolved certain/inferred direction; visibly annotate inferred direction. |
| Halogen bond | Directional | `halogen_donor`, `acceptor`; donor anchor retained separately | Halogen donor → acceptor when supported. |
| Metal coordination | Role-specific | `metal_center`, `coordinating_atom` | Metal — coordinating atom; do not append an oxidation state without charge evidence. |
| Salt bridge | Role-specific | `positive_site`, `negative_site` | Positive site ↔ negative site, with charge provenance. |
| Cation–π | Role-specific | `cation`, `aromatic_ring` | Cation ↔ identified physical ring. |
| π–π | Symmetric | Two `aromatic_ring` sites | Ring ↔ ring. |
| Hydrophobic contact | Symmetric | `contact_atom` or an explicitly represented group | Atom/group ↔ atom/group. |
| Clash | Symmetric | Two `contact_atom` sites | Atom ↔ atom; separation and overlap describe a geometric conflict. |
| Base-pair context | Residue-pair metadata alongside the atom-level family | Preserve the atom-level donor/acceptor or unresolved roles | Context never replaces the atom-level family or supplies a donor arrow. |
| Generic polar/packing/proximity | Symmetric unless a supported family supplies more specific roles | Retained contact sites, or `unresolved` | Do not manufacture donor/acceptor or charge roles from proximity. |

Direction certainty describes the role assignment under the available chemical model. `certain` does not mean experimental confirmation of a bond, true protonation in solution, or known binding energetics. `inferred` must remain visible as an inference. An `ambiguous` direction has no single asserted donor arrow.

Common role-provenance values include `arpeggio_openbabel_atom_types`, `standard_residue_template`, `explicit_protonation_label`, `coordinate_hydrogen`, `engine_assignment`, `element_identity` and `ring_membership`. A standard template is a chemical assumption, not a residue-name shortcut to be repeated independently in the UI. The backend owns that inference and its uncertainty.

When actual nucleobase atom typing conflicts with the neutral-residue template, the affected participant is `unresolved` and `roleAlternatives` records each source's role. The assignment is ambiguous and has no asserted donor direction. A canonical G–C or A–T residue pattern cannot override incompatible actual atom roles. Current reports retain `basePair` as context on `hbond` or `polar_contact` records; they do not create a separate atom-level `base_pairing` family.

Charge provenance distinguishes perceived Open Babel formal charge, an engine ionisable type and a standard-residue charge template. When only a conditional charge sign is inferred, `value` stays null and `sign` records that conditional positive/negative role; it is not a formal charge magnitude of +1 or −1. `formal_charge_sum_of_site_atoms` with `scope: "group"` is the sum over a charged feature, while `atomValue` preserves the charge of the distance endpoint atom. These are distinct quantities. A perceived formal charge can still be inferred, because it need not have been explicitly supplied in the coordinate file. Metal identity alone supplies no oxidation state.

For salt bridges and cation–π, the required charged-site provenance sets `chemicalState`: supplied non-inferred formal charges are `input`, perceived formal charges are `modelled`, assumed charge signs are `assumed`, and unresolved required charged sites are `unknown`. Supplied charge is an input-model state, not experimental verification of the solution state.

## Measurement semantics

| `geometry.distance.kind` | Meaning |
| --- | --- |
| `donor_acceptor` | Donor–acceptor heavy-atom separation; not H–acceptor distance. |
| `metal_donor` | Metal-center to coordinating-atom separation. |
| `ring_centroid` | Separation between the two specified ring centroids. |
| `cation_centroid` | Cationic site to the specified ring centroid. |
| `atom_ring_centroid` | Specified atom to the selected physical ring centroid, including native engine atom–ring observations; not an atom-pair distance. |
| `charge_site` | Separation of the specified charged-site contact atoms. A group can retain the full charged feature while `site.contactAtom` identifies the atom used for this measurement. |
| `closest_atom` | Closest constituent-atom separation for the specified features. |
| `atom_pair` | Separation of the named atom pair. |

Other measured quantities have their own `kind`, `unit` and `source`; examples include donor–H–acceptor angle, donor-anchor–halogen···acceptor angle, ring-normal angle, interplanar distance, lateral offset, sum of van der Waals radii and overlap. Angles based on an actual supplied hydrogen and angles based on a generated hydrogen have distinct provenance. A heavy-atom hydroxyl proxy is never a measured donor–H–acceptor angle.

`cation_ring_normal_angle` measures the angle of the centroid-to-cation vector from the ring normal. The current cation–π support screen requires a centroid distance at most 4.5 Å and that angle at most 30°, as well as the implemented lower distance bound. A cation lying near the ring plane does not gain geometric support from centroid proximity alone.

`geometry.hydrogen.source` distinguishes `input_hydrogens` from `arpeggio_generated_hydrogens`. The hydrogen retains its identity, coordinates and inference flag; `bondSource` separately records how its donor association was obtained. `geometry.donorAnchor.source` distinguishes the engine bond graph from coordinate-based bond inference. Neither inferred hydrogens nor an inferred bond graph establish experimentally measured bonding.

Ring membership identifies the actual physical site. Five- and six-membered rings in fused systems remain distinct. A ring normal is an orientation vector, not a physical force. Missing ring coordinates cannot be replaced by an arbitrary residue atom and then reported as a centroid measurement. Surface-projected, offset or animated connector lengths are schematic and are never independent numerical measurements.

## Evidence levels

The backend assigns `evidence.level`; consumers display it rather than reconstructing a level from distance or legacy words such as `strong` and `confirmed`. The separate compatibility, geometry, state and ambiguity fields are retained even when they differ in support. No level implies a stabilizing energy or biological function.

- `candidate` preserves a plausible interaction that is not fully supported by the required chemistry/state/geometry. Heavy-atom proximity alone can support candidacy, not direct H-bond confirmation. Generated-hydrogen H-bond geometry remains a model-dependent candidate.
- `chemically_supported` identifies chemical role compatibility while retaining any incomplete geometric or state evidence. The current salt-bridge implementation uses it for compatible separation with opposite nonzero formal site charges; assumed charge signs remain candidates. It does not mean that the displayed separation measures energy.
- `geometrically_supported` identifies an assignment satisfying the family's implemented geometric checks under the recorded chemical model. For H-bonds, this requires supplied, non-proxy hydrogen geometry, a valid donor association, compatible heavy-atom and H–acceptor separations, and the required donor–H–acceptor angle. Generated-H geometry remains candidate. π families must pass their implemented separation, orientation and offset checks; available centroids alone are insufficient. Aromatic packing instead uses its closest-atom separation range.
- `ambiguous` explicitly marks unresolved chemically meaningful roles or alternatives. It does not erase measured distances or any compatibility evidence that is available.

These labels are categorical summaries rather than a universal ranked scale. In particular, a certain donor/acceptor direction can coexist with a candidate interaction, and a precisely calculated heavy-atom distance can coexist with unknown protonation. Consult the family-specific classification rules and recorded flags for the exact reason a level was assigned.

## Identity, normalization and pruning

`semantics.identity` is family-aware and constructed from the chemical roles and physical site IDs, including model/alternate-conformer identity and a relevant explicit hydrogen ID. Reversed storage order is the same interaction only if the roles and physical identities agree. Different families, donor atoms, charged features or physical rings remain distinct even if the residue pair is identical.

Normalization and chain-alias remapping must carry participant residues, atom/site IDs, direction endpoints, hydrogen/anchor metadata, geometry and provenance together. The backend external-chain projection transforms its internal site/atom IDs and rebuilds the semantic identity once in the external chain namespace. Subsequent frontend display aliases remap nested residue references while preserving those exported opaque identities. Consumers must not discard this object and reconstruct roles from atom names. A UI group is a presentation container, not a replacement canonical interaction. Connector limits may hide drawings for readability; they must not remove distinct canonical edges or their detail rows.

Equivalent observations merge role provenance and typing while retaining the source geometry, participants, direction and evidence in `sourceObservations`. Snapshot participants reference physical site IDs rather than repeating full ring atom coordinates; their source role, charge and typing remain available. If assertion replaces an atom endpoint with another member of the same charged group, `sourceEndpoints` retains the original atom and its typing/charge. Conflicting numerical charges for the same physical site are flagged as ambiguous; different atomic charges within one resonance-delocalized group do not by themselves imply different group charge. `supportingRecords` counts merged observations, not independent interaction energies or experimental replication.

Canonical versus diagnostic API projections preserve the same surviving contact records and their semantics. Diagnostic omission is explicit and separate from evidence level: a diagnostic proximity record can retain a precisely measured separation without meeting a specific interaction family's criteria. Base-pair context is additional metadata, not a duplicate H-bond or an independent count of unique paired residues.

Family-specific display aggregation is a separate layer documented in [contact-aggregation.md](contact-aggregation.md). `report.displayGroups` preserves canonical member identities and supporting atom evidence while describing natural feature/edge units. In particular, many hydrophobic atom contacts can support one local nonpolar region. Source-observation, canonical-contact, supporting-atom, assigned-unit and display-object counts have different meanings.

## Explanation contract

`api.explain.describe_contact()` consumes the versioned roles and measurement metadata directly. A legacy contact without supported semantics remains countable, but its direction, charge and distance type are reported as unavailable rather than guessed. A supplied partial record cannot cause a missing ring metric to be relabeled as a centroid distance.

`explain_report()` keeps non-diagnostic canonical record counts separate from available grouping counts. When versioned grouping metadata is complete, it leads with family-specific assigned units, unique known atom support and available display objects; residue rankings use assigned units. Legacy or incomplete data retain explicitly labeled record counts/rankings. Actual viewport connector counts are unavailable to the backend. It describes evidence-level coverage and shows one representative example per retained family; selection is not an energetic ranking. This bound avoids constructing thousands of detailed narratives for large reports. It does not discard any contact from the canonical report. A future focused explanation can use `describe_contact()` for any selected record.

## Scientific context

The evidence separation follows the distinction between structural screening and demonstrated bonding in the [IUPAC hydrogen-bond definition](https://publications.iupac.org/pac/83/8/1637/index.html) and [IUPAC halogen-bond definition](https://publications.iupac.org/pac/85/8/1711/index.html). Atom typing and geometric detection build on [Arpeggio's methods](https://pubmed.ncbi.nlm.nih.gov/27964945/). Physical ring/charged-site roles matter for [cation–π contacts](https://pubmed.ncbi.nlm.nih.gov/10449714/), and atom-level contacts should not replace the richer edge/orientation description of [RNA base pairing](https://pubmed.ncbi.nlm.nih.gov/11345429/). Metal identity and complete coordination-site validation require additional checks, as illustrated by [CheckMyMetal](https://pubmed.ncbi.nlm.nih.gov/24356774/).

These sources motivate the chemical distinctions. Roami's numerical thresholds and serialized field contract are implementation choices, documented here and in `arpeggio-interaction-rules.md`.
