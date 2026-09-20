# Contact evidence, chemical units and display groups

Roami preserves atom-level assignments while grouping dense depictions by family-specific chemical features. One local nonpolar region can have many supporting atom contacts. A metal coordination site can contain several distinct metal–donor edges. A base-pair context can contain both H-bond and polar assignments. These are different units and must not be added together as if each were an independent bond or an equal energetic contribution.

The [Arpeggio methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC5282402/) distinguish atom, atom–ring and ring–ring observations. [PLIP's methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC4489249/) and [official documentation](https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md?plain=1) illustrate reducing dense hydrophobic diagrams while retaining a declared screening policy. Roami's feature partitioning and spatial thresholds are its own reporting rules, not universal chemical definitions or energy models.

## Counting units

| Count | Meaning | Does not mean |
| --- | --- | --- |
| Source observations | Engine evidence records, including repeated/reversed observations where retained | Unique atom contacts or independent measurements |
| Canonical contacts | Retained versioned assignments after exact chemical-identity deduplication | A universal number of physical bonds |
| Supporting atom contacts | Unique exported atom-contact identities, including actual H identity where relevant | All possible atom combinations inside a ring or feature |
| Assigned chemical units | Family-specific feature/edge units under the grouping policy | Confirmed bonds, energetic strength or additive stabilization |
| Available display objects | Groups available for drawing and inspection | Connectors actually visible after view filters and render limits |
| Rendered objects | Objects actually drawn in the current viewport | A changed amount of underlying chemical evidence |

`semantics.evidence.supportingRecords` counts source observations. It must not be repurposed as supporting atom-contact count. Family views can overlap, and the same supporting atom identity can occur in multiple contexts; use identity-aware totals rather than summing all view labels. Base-pair context is metadata on existing groups, not an additional group contributing to a unique-contact total.

## Version 1 grouping contract

```text
report.displayGroups = {
  version: 1,
  groups: [{
    id, family, rule,
    featureIds: [featureIdA, featureIdB],
    contactIds: [canonicalSemanticsIdentity, ...],
    representativeContactId,
    debugOnly?: boolean,
    counts: {
      canonicalContacts: integer,
      supportingAtomContacts: integer | null,
      chemicalUnits: 1,
      displayObjects: 1
    },
    geometry: {
      closestAtomDistance: number | null,
      meanAtomDistance?: number | null,
      medianAtomDistance?: number | null,
      centroidDistance?: number | null
    },
    renderGeometry?: {
      kind: "hydrophobic_patch",
      points: [[x, y, z], [x, y, z]],
      schematic: true,
      source: "supporting_atom_patch_centroids"
    },
    ambiguityFlags: [string, ...],
    supportingAtomContacts: [{
      id, canonicalContactId,
      atomA: {id, atomName, element, coordinates, modelId, altloc, residue},
      atomB: {id, atomName, element, coordinates, modelId, altloc, residue},
      distance: number | null,
      source,
      typing: {A: {atomTypes, roleProvenance, charge?, role?, roleAlternatives?}, B: {...}},
      hydrogen?, provenance: [...]
    }, ...]
  }, ...],
  features: {featureId: {id, kind, label?, residue, atoms, centroid?, normal?, provenance?, sourceSiteIds?}, ...},
  policy: {hydrophobicEndpointDiameterAngstrom: 3.5, interpretation: "display_grouping_only"}
}
```

All distance fields are in ångströms. `contactIds` refer to preserved canonical `semantics.identity` values; contact-level `displayGrouping` metadata references the group instead of duplicating its evidence. Group keys canonicalize endpoint order, while contact-level feature references remain associated with their A/B participants. Diagnostic and ordinary records do not mix in one group.

`meta.displayGroupsVersion: 1` identifies this reporting layer. Optional `basePairContexts` on a group retains member-specific pairing metadata; it does not create an overlapping counted group.

Supporting atom identities are invariant to reversal. The original source atoms survive charged-group endpoint reassignment. An atom–ring or ring–ring source observation without exported atom support gives its group `supportingAtomContacts: null` in its counts and `atom_support_not_exported` uncertainty. The support list retains any known pairs from other group members and is empty when none were exported. Consumers must not invent the Cartesian product of ring atoms or treat unavailable support as zero.

A native hydrophobic/packing/generic atom–ring observation retains the physical ring participant and `atom_ring_centroid` distance. That metric must not be attached to an arbitrary ring atom or called its closest atom distance. If no supporting atom pair was exported for any group member, the group's closest supporting atom distance remains unavailable even though atom-to-ring-center geometry is known. A native feature observation can accompany a unique compatible same-feature patch; several compatible patches require explicit selection uncertainty rather than arbitrary merging.

## Family-specific units

| Family | Unit preserved for display | Required retained detail |
| --- | --- | --- |
| Hydrophobic | Local contact region between nonpolar features | Every supporting atom pair, feature members, typing/connectivity provenance and grouping rule |
| H-bond | Donor–acceptor assignment with actual H identity where present | Roles, direction/ambiguity, H geometry and alternatives; distinct donor/acceptor edges remain distinct |
| Halogen bond | Halogen-donor/acceptor edge | Covalent donor anchor, measured angle and role provenance |
| Salt bridge | Charged-feature pair | Constituent atoms, original supporting pairs, group charge and exact distance endpoint |
| Metal coordination | Metal–donor edge, optionally contained in a coordination-site view | Each donor edge; bidentate coordination retains two donors |
| Cation–π | Cation feature and physical ring | Cation site, ring membership and centroid/normal geometry |
| π–π | Physical ring pair | Both ring memberships, selected geometry and alternative-site uncertainty |
| Aromatic packing | Aromatic feature-pair region | Closest atom support; grouping does not promote packing to π stacking |
| Polar/packing/proximity | Atom-pair or native atom/ring assignment, optionally in a presentation container | Exact sites, distance kind, uncertainty and diagnostic status |
| Clash | Actual atom-pair overlap | Atom identities, radius reference, separation and overlap |
| Base-pair context | Context metadata alongside the atom-level units | All constituent assignments; no extra interaction count or automatic H-bond classification |

Nonpolar feature assignment uses chemical typing and connectivity rather than one feature for every ligand residue. Separate nonpolar segments divided by polar/charged linkers and distinct aromatic/aliphatic features remain identifiable. Missing topology or ambiguous feature membership needs an explicit conservative fallback.

Salt-bridge and cation–π grouping use the assigned charged feature, including known standard-residue groups such as arginine guanidinium, rather than creating a separate chemical unit for every atom belonging to that feature. The individual canonical atoms and their measurements remain inspectable. If charge-group membership is unresolved, an atomic charge-site fallback carries that uncertainty; the consumer does not assume the entire ligand is one charged group. Feature membership does not establish protonation, oxidation state or an energetic contribution.

Hydrophobic local patches use a complete-link endpoint-distance policy: supporting contacts in one patch have a maximum separation of 3.5 Å on each corresponding feature endpoint. Checking both endpoints prevents crossed connectors with similar midpoints from merging distant contact regions; complete linkage prevents a long chain of small spatial steps from merging remote patches. This threshold controls presentation granularity, not the existence or energy of a hydrophobic bond.

The representative is chosen by a declared geometric policy and stable identity tie-breaking. The shortest supporting atom distance is a measured model quantity. A patch-center separation is a schematic depiction metric; its line length and the mean/median supporting distances are not bond energies or strength rankings. Selecting a representative does not remove the other atom contacts.

Fused rings retain their physical memberships. The display catalog reuses canonical ring sites and unifies identical scoped atom membership across descriptor or chain-alias paths, retaining source site IDs as provenance. A different descriptor hash cannot create a second physical feature for the same ring. Different five-/six-membered memberships remain distinct. If an observation is compatible with multiple overlapping ring choices, preserve selection uncertainty and alternatives. Several ring depictions from the same aromatic system must not be described as additive independent bonds.

## Explanation and fallback behavior

`api.explain.display_group_summary()` consumes the producer's groups without reclustering. It checks that each counted group belongs to the retained, non-diagnostic canonical scope, deduplicates group/support identities and reports coverage. Its assigned-unit and available-object totals are separate from canonical-record counts and the count of unique known atom support.

When grouping metadata is complete, residue rankings count assigned units, so a dense hydrophobic region contributes one unit per involved residue. Incomplete or legacy reports retain explicitly labeled canonical-record rankings; missing grouping or atom-support coverage is stated rather than reconstructed. Examples remain bounded to one representative group per family.

The backend explanation cannot know the current viewport's actual rendered-object count. Visibility and connector budgets belong to the UI. Neither changing a view nor increasing the number of atom rows establishes stronger binding.
