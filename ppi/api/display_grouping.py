"""Lossless, family-specific display projection over canonical chemical records.

The 3.5 Å endpoint-diameter bound is a display policy, not a hydrophobic energy
criterion. This module never reclassifies, removes or strengthens a contact.
"""
import copy
import hashlib
import math
import statistics
from dataclasses import replace

PATCH_DIAMETER = 3.5
# Nonaromatic side-chain carbon segments. Polar linkers and carbonyl carbons
# deliberately break segments; absent atoms do not change known connectivity.
ALKYL_SEGMENTS = {
    'ALA': [('CB',)], 'VAL': [('CB','CG1','CG2')],
    'LEU': [('CB','CG','CD1','CD2')], 'ILE': [('CB','CG1','CG2','CD1')],
    'LYS': [('CB','CG','CD','CE')], 'ARG': [('CB','CG','CD')],
    'PRO': [('CB','CG','CD')], 'GLU': [('CB','CG')], 'GLN': [('CB','CG')],
    'ASP': [('CB',)], 'ASN': [('CB',)], 'MET': [('CB','CG'),('CE',)],
    'MSE': [('CB','CG'),('CE',)], 'PHE': [('CB',)], 'TYR': [('CB',)],
    'TRP': [('CB',)], 'HIS': [('CB',)], 'THR': [('CG2',)],
}
FAMILY_RULES = {
    'hbond': 'donor_hydrogen_acceptor_edge',
    'metal_coordination': 'metal_donor_edge',
    'halogen_bond': 'halogen_acceptor_edge',
    'clash': 'clashing_atom_pair', 'steric_clash': 'clashing_atom_pair',
    'salt_bridge': 'charge_group_pair',
    'pi_pi': 'physical_ring_pair', 'aromatic_packing': 'physical_ring_pair',
    'aromatic_proximal': 'physical_ring_pair', 'pi_cation': 'cation_physical_ring_pair',
    'polar_contact': 'typed_polar_atom_pair',
}


def _id(prefix, values):
    return prefix + ':' + hashlib.sha256('\x1f'.join(map(str, values)).encode()).hexdigest()[:24]


def _xyz(atom):
    xyz = atom.get('coordinates') if isinstance(atom, dict) else None
    return xyz if isinstance(xyz, (list, tuple)) and len(xyz) == 3 and all(isinstance(v,(int,float)) and math.isfinite(v) for v in xyz) else None


def _center(atoms):
    points = [_xyz(atom) for atom in atoms]
    return [round(sum(p[i] for p in points)/len(points),6) for i in range(3)] if points and all(points) else None


def _distance(a, b):
    return round(math.dist(a,b),6) if a is not None and b is not None else None


def _diagnostic(contact):
    return bool(contact.get('debugOnly') or contact.get('debug_only') or contact.get('asserted',{}).get('debugOnly'))


class FeatureCatalog:
    def __init__(self, index, aliases):
        from . import analysis as a
        self.a = a
        self.index = a.ResidueAtomIndex(replace(atom, chain_id=a._external_chain_id(atom.chain_id, aliases)) for atoms in index.values() for atom in atoms)
        self.index.chemical_nodes = {(a._external_chain_id(chain,aliases),seq,name):node for (chain,seq,name),node in index.chemical_nodes.items()}
        self.features = {}
        self.ring_membership = {}
        self.hydrophobic = {}
        self.atom_by_key = {}
        for atoms in self.index.values():
            for atom in atoms:
                residue = {'chain':atom.chain_id,'seq':atom.res_seq,'resName':atom.res_name,'atom':atom.atom_name}
                item = {**a._semantic_atom(residue,{},self.index), 'residue':{k:residue[k] for k in ('chain','seq','resName')}}
                self.atom_by_key[a._build_atom_key_from_payload(residue)] = item
                self.atom_by_key[item['id']] = item

    def register(self, feature):
        membership = tuple(sorted(atom['id'] for atom in feature.get('atoms',[]))) if feature.get('kind')=='ring' else ()
        if len(membership)>=3:
            existing_id = self.ring_membership.get(membership)
            if existing_id:
                existing = self.features[existing_id]
                existing['sourceSiteIds'] = sorted(set(existing.get('sourceSiteIds',[])) | {feature['id']})
                if feature.get('label') and not existing.get('label'):
                    existing['label'] = feature['label']
                return existing
            self.ring_membership[membership] = feature['id']
            feature = {**feature,'sourceSiteIds':[feature['id']]}
        if feature['id'] not in self.features:
            self.features[feature['id']] = copy.deepcopy(feature)
        return self.features[feature['id']]

    def atom_feature(self, atom, flags=()):
        site = {'id':atom['id'],'kind':'atom','label':atom['atomName'],'residue':atom['residue'],
                'atoms':[{k:v for k,v in atom.items() if k != 'residue'}], 'provenance':['coordinate_atom_identity']}
        if flags:
            site['ambiguityFlags'] = list(flags)
        return self.register(site)

    def make_feature(self, residue, names, kind, label, provenance):
        atoms = [self.atom_by_key.get(f"{residue['chain']}:{residue['seq']}:{residue['resName']}:{name}") for name in sorted(names)]
        atoms = [atom for atom in atoms if atom]
        if not atoms:
            return None
        feature = {'id':_id('display-feature', [kind]+[atom['id'] for atom in atoms]),'kind':kind,'label':label,
                   'residue':dict(residue),'atoms':[{k:v for k,v in atom.items() if k != 'residue'} for atom in atoms],
                   'centroid':_center(atoms),'provenance':[provenance]}
        if kind == 'ring':
            descriptor = next((item for item in self.a._residue_ring_descriptors(residue,self.index) if set(item['atom_names'])==set(names)),None)
            if descriptor is not None:
                identity = self.a._ring_site_key_from_descriptor(residue,descriptor)
                scope = atoms[0]
                feature['id'] = identity + f"@model={scope['modelId'] or '?'}@alt={scope['altloc'] or '.'}"
            normal = self.a._best_plane_normal([atom['coordinates'] for atom in atoms])
            if normal is not None:
                feature['normal'] = list(normal)
        return self.register(feature)

    def _hydro_residue(self, residue):
        key = (residue['chain'],residue['seq'])
        if key in self.hydrophobic:
            return self.hydrophobic[key]
        atoms = self.index.get(key,[])
        nodes = {atom.atom_name:self.index.chemical_nodes.get((*key,atom.atom_name),{}) for atom in atoms}
        mapping = {atom.atom_name:[] for atom in atoms}
        ring_sets = set()
        # Actual perceived physical rings for ligands; known standard rings can
        # also be recovered when a controlled/older exporter lacks bond features.
        for node in nodes.values():
            for names in node.get('aromatic_rings',[]):
                ring_sets.add(tuple(sorted(names)))
        standard_rings = self.a._residue_ring_atom_names(residue['resName'])
        if standard_rings:
            for descriptor in self.a._residue_ring_descriptors(residue,self.index):
                ring_sets.add(tuple(sorted(descriptor['atom_names'])))
        ring_names = set()
        for names in sorted(ring_sets):
            feature = self.make_feature(residue,names,'ring',f"{residue['resName']} ring ({', '.join(names)})",
                                        'arpeggio_openbabel_ring_graph' if any(list(names) in node.get('aromatic_rings',[]) for node in nodes.values()) else 'standard_residue_ring_template')
            if feature:
                for name in names:
                    if name in mapping:
                        mapping[name].append(feature)
                        ring_names.add(name)
        segments = ALKYL_SEGMENTS.get(residue['resName'])
        if segments is not None:
            for names in segments:
                # Explicit polar/charged typing vetoes nonpolar segment membership.
                selected = [name for name in names if name in mapping and name not in ring_names and
                            not (set(nodes.get(name,{}).get('atom_types',[])) & {'hbond donor','hbond acceptor','pos ionisable','neg ionisable'})]
                feature = self.make_feature(residue,selected,'group',f"{residue['resName']} alkyl segment ({', '.join(selected)})",'standard_residue_sidechain_connectivity')
                if feature:
                    for name in selected:
                        mapping[name].append(feature)
        else:
            # Only actual typed, graph-connected nonaromatic carbons. Never infer
            # an entire ligand feature from its residue name or spatial proximity.
            eligible = {atom.atom_name for atom in atoms if atom.element.upper()=='C' and atom.atom_name not in ring_names and
                        'hydrophobe' in nodes[atom.atom_name].get('atom_types',[]) and
                        not (set(nodes[atom.atom_name].get('atom_types',[])) & {'hbond donor','hbond acceptor','pos ionisable','neg ionisable'})}
            adjacency = {name:set() for name in eligible}
            for name in eligible:
                for neighbor in nodes[name].get('bonded_atoms',[]):
                    other = neighbor.get('auth_atom_id')
                    if other in eligible and neighbor.get('same_residue') is True and not neighbor.get('aromatic') and neighbor.get('bond_order',1)==1:
                        adjacency[name].add(other); adjacency[other].add(name)
            remaining = set(eligible)
            while remaining:
                component, stack = set(),[min(remaining)]
                while stack:
                    current = stack.pop()
                    if current in component:
                        continue
                    component.add(current); stack.extend(adjacency[current]-component)
                remaining -= component
                if len(component)>1:
                    feature = self.make_feature(residue,component,'group',f"Nonpolar carbon segment ({', '.join(sorted(component))})",'arpeggio_openbabel_bond_graph')
                    for name in component:
                        mapping[name].append(feature)
        self.hydrophobic[key] = mapping
        return mapping

    def for_participant(self, participant, family):
        site = participant['site']
        if family in {'pi_cation','salt_bridge'} and participant.get('role') in {'cation','positive_site','negative_site'}:
            sign = 'anion' if participant['role']=='negative_site' else 'cation'
            residue = dict(site['residue'])
            endpoint = site.get('contactAtom') or next(iter(site.get('atoms',[])),{})
            residue.update(atom=endpoint.get('atomName'),element=endpoint.get('element'))
            group_atoms = self.a._collect_salt_bridge_site_atoms(residue=residue,res_name=residue['resName'],site_kind=sign,residue_atoms_index=self.index)
            if len(group_atoms)>1:
                group = self.a._semantic_site(residue,{},self.index,'group',names=[atom.atom_name for atom in group_atoms])
                group.update(label=f"{residue['resName']} {'guanidinium group' if residue['resName']=='ARG' else 'charge group'}",provenance=['standard_residue_charge_group_membership'])
                return self.register(group)
            if site.get('kind')=='atom' and residue['resName'] not in self.a.POLYMER_RESIDUES and residue['resName'] not in self.a.METAL_ELEMENTS:
                return self.register({**site,'label':'Atomic charge site','ambiguityFlags':['charge_group_membership_unresolved']})
        if family != 'hydrophobic':
            return self.register(site)
        if site.get('kind') == 'ring':
            return self.register(site)
        residue = site['residue']
        atom = site.get('contactAtom') or next(iter(site.get('atoms',[])),None)
        if not atom:
            return self.register(site)
        matches = self._hydro_residue(residue).get(atom['atomName'],[])
        if len(matches)==1:
            return matches[0]
        return self.atom_feature({**atom,'residue':residue}, ['shared_fused_ring_atom'] if len(matches)>1 else ['nonpolar_feature_unresolved'])


def _support_for_contact(contact, catalog):
    sem = contact['semantics']
    observations = sem.get('sourceObservations') or [{
        'atomPair':[contact.get('atomKeyA'),contact.get('atomKeyB')],
        'type':contact.get('arpeggio',{}).get('type'), 'source':contact.get('source'),
        'participants':sem['participants'],'sourceEndpoints':sem.get('sourceEndpoints'),
        'geometry':sem['geometry'],
    }]
    supports, unknown = {}, False
    canonical_sides = {p['side']:p for p in sem['participants']}
    for observation in observations:
        endpoints = observation.get('sourceEndpoints')
        participants = {p['side']:p for p in observation.get('participants',sem['participants'])}
        resolved = []
        endpoint_typing = []
        if endpoints and len(endpoints)==2:
            by_side = {p['side']:p for p in endpoints}
            for side in ('A','B'):
                endpoint = by_side.get(side,{})
                atom = copy.deepcopy(endpoint.get('atom') or {})
                residue = participants.get(side,{}).get('site',{}).get('residue')
                if atom and residue:
                    atom['residue'] = copy.deepcopy(residue)
                resolved.append(atom if atom and residue else None)
                endpoint_typing.append({**participants.get(side,{}),**endpoint})
        else:
            for side, token in zip(('A','B'),observation.get('atomPair') or (None,None)):
                # Comma-separated ring membership is not an atom. Never generate
                # a cross-product or relabel a centroid distance as atom support.
                resolved.append(copy.deepcopy(catalog.atom_by_key.get(token)) if token and ',' not in token else None)
                endpoint_typing.append(participants.get(side,{}))
        if len(resolved)!=2 or not all(resolved):
            unknown = True
            continue
        source_type = observation.get('type') or contact.get('arpeggio',{}).get('type')
        if source_type not in (None,'atom-atom') and not endpoints:
            unknown = True
            continue
        # Source observations can be true reversals; orient to this contact's A/B.
        def belongs(atom, participant):
            site = participant['site']
            return atom['id'] in {at['id'] for at in site.get('atoms',[])} or atom['residue']==site.get('residue')
        original_a = participants.get('A',{}).get('site',{}).get('id')
        reversed_identity = original_a == canonical_sides['B']['site']['id'] and original_a != canonical_sides['A']['site']['id']
        if reversed_identity or (not belongs(resolved[0],canonical_sides['A']) and belongs(resolved[1],canonical_sides['A'])):
            resolved.reverse(); endpoint_typing.reverse()
        hydrogen = observation.get('geometry',{}).get('hydrogen')
        token = _id('atom-support', sorted(atom['id'] for atom in resolved)+[hydrogen.get('id','') if hydrogen else ''])
        support = {'id':token,'canonicalContactId':sem['identity'],'atomA':resolved[0],'atomB':resolved[1],
                   'distance':_distance(_xyz(resolved[0]),_xyz(resolved[1])),
                   'source':observation.get('source') or contact.get('source') or 'pdbe-arpeggio',
                   'typing':{side:{key:copy.deepcopy(value) for key,value in item.items() if key in {'atomTypes','roleProvenance','charge','role','roleAlternatives'}} for side,item in zip(('A','B'),endpoint_typing)},
                   'provenance':['source_atom_coordinates']}
        if hydrogen:
            support['hydrogen'] = copy.deepcopy(hydrogen)
        if token in supports:
            for side in ('A','B'):
                for key in ('atomTypes','roleProvenance'):
                    supports[token]['typing'][side][key] = sorted(set(supports[token]['typing'][side].get(key,[])) | set(support['typing'][side].get(key,[])))
        else:
            supports[token] = support
    return list(supports.values()), unknown


def _patch_points(item):
    points = [[],[]]
    for support in item['supports']:
        for i,key in enumerate(('atomA','atomB')):
            xyz = _xyz(support[key])
            if xyz is not None:
                points[i].append(xyz)
    if not all(points) and item['unknown']:
        # A source atom–ring observation can locate a patch without claiming that
        # any particular cross-feature atom pair was exported by the engine.
        participants = sorted(item['contact']['semantics']['participants'],key=lambda p:p['side'])
        points = [[xyz for atom in participant['site'].get('atoms',[]) if (xyz:=_xyz(atom)) is not None] for participant in participants]
    if item['featureIds'][0] > item['featureIds'][1]:
        points.reverse()
    return points


def _coherent(patch,item):
    points = _patch_points(item)
    if not all(points):
        return False
    if not all(all(math.dist(p,q)<=PATCH_DIAMETER+1e-9 for i,p in enumerate(side) for q in side[i+1:]) for side in points):
        return False
    for old in patch:
        old_points = _patch_points(old)
        if not all(old_points):
            return False
        for side in (0,1):
            if any(math.dist(p,q)>PATCH_DIAMETER+1e-9 for p in points[side] for q in old_points[side]):
                return False
    return True


def attach_display_groups(report, residue_atoms_index, aliases):
    catalog = FeatureCatalog(residue_atoms_index,aliases)
    # Canonical ring hashes can contain the engine's internal chain namespace.
    # Reuse an existing canonical site identity, joining equivalent physical rings
    # by exact scoped membership rather than recomputing that opaque digest.
    canonical_rings = {p['site']['id']:p['site'] for rows in report.get('contacts',{}).values() for contact in rows
                       for p in contact.get('semantics',{}).get('participants',[]) if p.get('site',{}).get('kind')=='ring'}
    for identity in sorted(canonical_rings):
        catalog.register(canonical_rings[identity])
    buckets = {}
    for rows in report.get('contacts',{}).values():
        for contact in rows:
            sem = contact.get('semantics') or {}
            if sem.get('version')!=1 or len(sem.get('participants',[]))!=2:
                continue
            participants = sorted(sem['participants'], key=lambda p:p['side'])
            features = [catalog.for_participant(p,sem['family']) for p in participants]
            feature_ids = [feature['id'] for feature in features]
            supports, unknown = _support_for_contact(contact,catalog)
            family = sem['family']
            rule = 'nonpolar_feature_pair_complete_link_patch' if family=='hydrophobic' else FAMILY_RULES.get(family,'individual_atom_pair')
            if rule=='individual_atom_pair' and any(feature['kind']=='ring' for feature in features):
                rule='atom_physical_ring_pair'
            # Non-hydrophobic canonical identities already encode the appropriate
            # physical ring/charge group or individual directional atom edge.
            grouped_features = family in {'pi_cation','salt_bridge','pi_pi','aromatic_packing','aromatic_proximal'} or rule=='atom_physical_ring_pair'
            natural = tuple(sorted(feature_ids)) if family=='hydrophobic' else tuple(sorted((p['role'],feature) for p,feature in zip(participants,feature_ids))) if grouped_features else sem['identity']
            key = (family,rule,natural,_diagnostic(contact))
            buckets.setdefault(key,[]).append({'contact':contact,'featureIds':feature_ids,'supports':supports,'unknown':unknown})
    groups = []
    for key,items in sorted(buckets.items(),key=lambda x:str(x[0])):
        family,rule,_,debug = key
        patches = []
        for item in sorted(items,key=lambda item:(item['unknown'],item['contact']['semantics']['identity'])):
            compatible = [p for p in patches if family!='hydrophobic' or _coherent(p,item)]
            target = compatible[0] if compatible and not (item['unknown'] and len(compatible)>1) else None
            if item['unknown'] and len(compatible)>1:
                item['ambiguousPatch'] = True
            if target is None:
                patches.append([item])
            else:
                target.append(item)
        for patch in patches:
            representative = min(patch,key=lambda item:(min((s['distance'] for s in item['supports'] if s['distance'] is not None),default=math.inf),item['contact']['semantics']['identity']))
            feature_ids = representative['featureIds']
            contacts = sorted(item['contact']['semantics']['identity'] for item in patch)
            group_id = _id('display-group',[family,rule,str(debug)]+sorted(feature_ids)+contacts)
            supports = {}
            flags = set()
            for item in patch:
                item['contact']['displayGrouping'] = {'groupId':group_id,'rule':rule,'featureIds':item['featureIds']}
                reverse = item['featureIds'] != feature_ids and item['featureIds'][::-1] == feature_ids
                for support in item['supports']:
                    support = copy.deepcopy(support)
                    if reverse:
                        support['atomA'],support['atomB'] = support['atomB'],support['atomA']
                        support['typing']['A'],support['typing']['B'] = support['typing']['B'],support['typing']['A']
                    supports.setdefault(support['id'],support)
                for feature_id in item['featureIds']:
                    flags.update(catalog.features[feature_id].get('ambiguityFlags',[]))
                if item['unknown']:
                    flags.add('atom_support_not_exported')
                if item.get('ambiguousPatch'):
                    flags.add('source_feature_matches_multiple_patches')
            ordered_supports = sorted(supports.values(),key=lambda s:s['id'])
            distances = [s['distance'] for s in ordered_supports if s['distance'] is not None]
            geometry = {'closestAtomDistance':min(distances) if distances else None,
                        'meanAtomDistance':round(statistics.mean(distances),6) if distances else None,
                        'medianAtomDistance':statistics.median(distances) if distances else None}
            group = {'id':group_id,'family':family,'rule':rule,'featureIds':feature_ids,'contactIds':contacts,
                     'representativeContactId':representative['contact']['semantics']['identity'],
                     'counts':{'canonicalContacts':len(patch),'supportingAtomContacts':None if 'atom_support_not_exported' in flags else len(supports),'chemicalUnits':1,'displayObjects':1},
                     'geometry':geometry,'ambiguityFlags':sorted(flags),'supportingAtomContacts':ordered_supports,'debugOnly':debug}
            if family=='hydrophobic' and supports:
                centers = []
                for side,feature_id in zip(('atomA','atomB'),feature_ids):
                    feature = catalog.features[feature_id]
                    unique = {s[side]['id']:s[side] for s in ordered_supports}
                    centers.append(feature.get('centroid') if feature['kind']=='ring' else _center(list(unique.values())))
                if all(centers):
                    group['renderGeometry'] = {'kind':'hydrophobic_patch','points':centers,'schematic':True,'source':'supporting_atom_patch_centroids'}
                    geometry['centroidDistance'] = _distance(*centers)
            contexts = {item['contact']['semantics']['identity']:copy.deepcopy(item['contact']['basePair']) for item in patch if item['contact'].get('basePair')}
            if contexts:
                group['basePairContexts'] = contexts
            groups.append(group)
    referenced = {feature for group in groups for feature in group['featureIds']}
    report['displayGroups'] = {'version':1,'groups':sorted(groups,key=lambda g:g['id']),
                               'features':{key:catalog.features[key] for key in sorted(referenced)},
                               'policy':{'hydrophobicEndpointDiameterAngstrom':PATCH_DIAMETER,'interpretation':'display_grouping_only'}}
    report.setdefault('meta',{})['displayGroupsVersion'] = 1
    return report


def project_display_groups(registry, contacts):
    """Drop diagnostic groups/features from delivery; canonical registry is intact."""
    identities = {row.get('semantics',{}).get('identity') for rows in contacts.values() for row in rows}
    groups = [group for group in registry.get('groups',[]) if group['contactIds'] and all(identity in identities for identity in group['contactIds'])]
    feature_ids = {feature for group in groups for feature in group['featureIds']}
    return {**registry,'groups':groups,'features':{key:value for key,value in registry.get('features',{}).items() if key in feature_ids}}
