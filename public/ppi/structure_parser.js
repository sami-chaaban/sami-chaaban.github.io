// Shared coordinate parsing for the viewer and its module worker. No DOM/Three.js dependencies.
const clean = value => value == null || value === '.' || value === '?' ? '' : String(value).trim();
const modelToken = value => { const token = clean(value); return /^[+-]?\d+$/.test(token) ? String(parseInt(token, 10)) : token; };
const numeric = value => { const token = clean(value); if (!token) return null; const n = Number(token); return Number.isFinite(n) ? n : null; };
const compareTokens = (a, b) => a.localeCompare(b, 'en', { numeric: true });

function* cifTokens(text) {
  let i = 0;
  while (i < text.length) {
    if (/\s/.test(text[i])) { i++; continue; }
    if (text[i] === '#') { while (i < text.length && text[i] !== '\n') i++; continue; }
    if (text[i] === ';' && (i === 0 || text[i - 1] === '\n')) {
      const start = ++i;
      const end = text.indexOf('\n;', i);
      if (end < 0) throw new Error('Unterminated CIF multiline value.');
      yield { value: text.slice(start, end).replace(/^\r?\n/, ''), quoted: true };
      i = end + 2;
      continue;
    }
    const quote = text[i] === "'" || text[i] === '"' ? text[i++] : '';
    const start = i;
    if (quote) {
      while (i < text.length && !(text[i] === quote && (i + 1 === text.length || /\s/.test(text[i + 1])))) i++;
      if (i === text.length) throw new Error('Unterminated CIF quoted value.');
      yield { value: text.slice(start, i++), quoted: true };
    } else {
      while (i < text.length && !/\s/.test(text[i])) i++;
      yield { value: text.slice(start, i), quoted: false };
    }
  }
}

export function parseCif(text, options = {}) {
  const iterator = cifTokens(String(text || ''));
  let current = iterator.next();
  const next = () => { const token = current.value; current = iterator.next(); return token; };
  const control = token => token && !token.quoted && (/^(?:loop_|stop_|global_|data_|save_)/i.test(token.value) || token.value.startsWith('_'));
  const tables = [];
  const items = new Map();
  while (!current.done) {
    const token = next();
    if (!token.quoted && token.value.toLowerCase() === 'loop_') {
      const columns = [];
      while (!current.done && !current.value.quoted && current.value.value.startsWith('_')) columns.push(next().value.toLowerCase());
      if (!columns.length) throw new Error('CIF loop has no columns.');
      const table = { columns, rows: [] };
      const atomTable = columns.some(name => name.startsWith('_atom_site.'));
      while (!current.done && !control(current.value)) {
        const row = [];
        for (let j = 0; j < columns.length; j++) {
          if (current.done || control(current.value)) throw new Error('Incomplete CIF loop row.');
          row.push(next().value);
        }
        if (options.onLoopRow) options.onLoopRow(columns, row);
        if (!atomTable || options.keepAtomRows !== false) table.rows.push(row);
      }
      tables.push(table);
    } else if (!token.quoted && token.value.startsWith('_')) {
      if (current.done || control(current.value)) throw new Error(`Missing CIF value for ${token.value}.`);
      items.set(token.value.toLowerCase(), next().value);
    }
  }
  return { tables, items };
}

function inferElement(name, residue, pdbField = '') {
  const atom = clean(name).toUpperCase().replace(/^\d+/, '');
  if (pdbField.startsWith(' ')) return atom.slice(0, 1);
  if (/^(CL|BR)/.test(atom)) return atom.slice(0, 2);
  const metals = new Set(['FE','ZN','MG','MN','CU','CO','NI','NA','CA','CD','HG','SR','SE']);
  const prefix = atom.slice(0, 2);
  if (metals.has(prefix) && (prefix !== 'CA' || clean(residue).toUpperCase() === 'CA')) return prefix;
  return atom.slice(0, 1);
}

export function selectConformers(records) {
  const groups = new Map();
  for (const atom of records) {
    const key = `${atom.modelToken}|${atom.chain}|${atom.labelChain || atom.chain}|${atom.seq}`;
    let group = groups.get(key);
    if (!group) { group = { atoms: [], alternatives: new Map() }; groups.set(key, group); }
    group.atoms.push(atom);
    if (atom.altLoc) {
      const score = group.alternatives.get(atom.altLoc) || { sum: 0, count: 0 };
      score.sum += atom.occupancy == null ? 1 : atom.occupancy;
      score.count++;
      group.alternatives.set(atom.altLoc, score);
    }
  }
  const retained = new Set();
  const rank = alt => alt === 'A' ? 0 : alt === '1' ? 1 : 2;
  for (const group of groups.values()) {
    const choices = [...group.alternatives].sort(([a, sa], [b, sb]) =>
      (sb.sum / sb.count - sa.sum / sa.count) || rank(a) - rank(b) || compareTokens(a, b));
    const chosen = choices.length ? choices[0][0] : '';
    const names = new Map();
    for (const atom of group.atoms) {
      if (atom.altLoc && atom.altLoc !== chosen) continue;
      const previous = names.get(atom.atomName);
      if (!previous || (previous.altLoc && !atom.altLoc) ||
          (previous.altLoc === atom.altLoc && (atom.occupancy ?? 1) > (previous.occupancy ?? 1))) names.set(atom.atomName, atom);
    }
    for (const atom of names.values()) retained.add(atom);
  }
  // File order is significant for sequence/token correspondence.
  return records.filter(atom => retained.has(atom));
}

function separateCollidingLabelChains(records) {
  const labelsByAuthor = new Map(), residueLabels = new Map(), collisions = new Set();
  for (const atom of records) {
    const author = atom.authChain || atom.chain, label = atom.labelChain || author;
    if (!labelsByAuthor.has(author)) labelsByAuthor.set(author, new Set());
    labelsByAuthor.get(author).add(label);
    const residue = `${atom.modelToken}|${author}|${atom.seq}`;
    if (residueLabels.has(residue) && residueLabels.get(residue) !== label) collisions.add(author);
    else residueLabels.set(residue, label);
  }
  if (!collisions.size) return false;
  // Allocate the same canonical namespace as the backend, including names
  // reserved by other author chains. Ordinary, noncolliding author groups keep
  // their familiar display name; only overlapping residue identities are split.
  const reserved = new Set(labelsByAuthor.keys()), canonical = new Map();
  for (const author of [...labelsByAuthor.keys()].sort()) {
    const labels = [...labelsByAuthor.get(author)].sort();
    if (labels.length < 2) continue;
    for (const label of labels) {
      let name = label;
      if (!name || reserved.has(name)) name = `${author}[${label || author}]`;
      let suffix = 2;
      while (reserved.has(name)) name = `${author}[${label || author}.${suffix++}]`;
      reserved.add(name);
      canonical.set(`${author}|${label}`, name);
    }
  }
  for (const atom of records) {
    const author = atom.authChain || atom.chain;
    if (collisions.has(author)) atom.chain = canonical.get(`${author}|${atom.labelChain || author}`) || atom.chain;
  }
  return true;
}

export function parseStructureRecords(text, format = 'mmcif') {
  const records = [];
  const tokens = new Set();
  let hasModelColumn = false;
  let hasExplicitModelRecords = false;
  let metadata = null;
  const connections = [];
  const push = atom => {
    if (!atom.seq || !atom.atomName || ![atom.x,atom.y,atom.z].every(Number.isFinite)) return;
    records.push(atom);
    if (atom.modelToken) tokens.add(atom.modelToken);
  };
  if (format === 'pdb') {
    let model = '1';
    let ordinal = 0;
    for (const line of String(text || '').split(/\r?\n/)) {
      const record = line.slice(0, 6).trim().toUpperCase();
      if (record === 'MODEL') { hasExplicitModelRecords = true; ordinal++; model = modelToken(line.slice(10).trim().split(/\s+/)[0]) || String(ordinal); tokens.add(model); continue; }
      if (record === 'ENDMDL') { model = ''; continue; }
      if (record === 'CONECT') { connections.push({ modelToken: model, serials: (line.slice(6).match(/\d+/g) || []).map(Number) }); continue; }
      if (record !== 'ATOM' && record !== 'HETATM') continue;
      const atomName = clean(line.slice(12,16));
      const resName = clean(line.slice(17,20));
      const chain = clean(line.slice(21,22)) || 'A';
      push({ atomName, resName, chain, authChain:chain, labelChain:chain,
        seq: clean(line.slice(22,26)) + clean(line.slice(26,27)),
        element: clean(line.slice(76,78)).toUpperCase() || inferElement(atomName,resName,line.slice(12,16)),
        x:numeric(line.slice(30,38)), y:numeric(line.slice(38,46)), z:numeric(line.slice(46,54)),
        bFactor:numeric(line.slice(60,66)), occupancy:numeric(line.slice(54,60)), altLoc:clean(line.slice(16,17)),
        serial:numeric(line.slice(6,11)), modelToken:model || '1', mmcifGroupPdb:record, mmcifLabelSeqId:'' });
    }
  } else {
    let lastColumns = null;
    let fields;
    metadata = parseCif(text, { keepAtomRows:false, onLoopRow(columns,row) {
      if (!columns[0].startsWith('_atom_site.')) return;
      if (columns !== lastColumns) { fields = new Map(columns.map((name,i)=>[name.slice(11),i])); lastColumns=columns; }
      const get = name => clean(row[fields.get(name)]);
      const authChain = get('auth_asym_id'), labelChain = get('label_asym_id');
      const atomName = get('auth_atom_id') || get('label_atom_id');
      const resName = get('auth_comp_id') || get('label_comp_id');
      const sequence = get('auth_seq_id') || get('label_seq_id');
      const model = modelToken(get('pdbx_pdb_model_num')) || '1';
      hasModelColumn ||= fields.has('pdbx_pdb_model_num');
      push({ atomName,resName,chain:authChain || labelChain || 'A',authChain,labelChain,
        authAtomName:get('auth_atom_id') || atomName,labelAtomName:get('label_atom_id') || atomName,
        seq:sequence ? sequence + get('pdbx_pdb_ins_code') : '',
        element:get('type_symbol').toUpperCase() || inferElement(atomName,resName),
        x:numeric(get('cartn_x')),y:numeric(get('cartn_y')),z:numeric(get('cartn_z')),
        bFactor:numeric(get('b_iso_or_equiv')), occupancy:numeric(get('occupancy')),
        altLoc:get('label_alt_id') || get('pdbx_pdb_alt_id'),modelToken:model,serial:numeric(get('id')),
        mmcifGroupPdb:get('group_pdb').toUpperCase(),mmcifLabelSeqId:get('label_seq_id'),mmcifEntityId:get('label_entity_id') });
    }});
  }
  const hasLabelChainCollisions = format !== 'pdb' && separateCollidingLabelChains(records);
  const selected = selectConformers(records);
  return { format, records:selected, modelTokens:[...tokens].sort(compareTokens),hasModelColumn,hasExplicitModelRecords,
    hasAlternateConformers:records.some(atom=>atom.altLoc),hasLabelChainCollisions, metadata,connections };
}

export function cifValue(value) {
  if (value == null || value === '') return '.';
  const text = String(value);
  if (!/[\s#'";]/.test(text) && !/^_|^(?:data_|loop_|stop_|save_|global_)/i.test(text)) return text;
  if (!text.includes('"') && !/[\r\n]/.test(text)) return `"${text}"`;
  if (!text.includes("'") && !/[\r\n]/.test(text)) return `'${text}'`;
  return `\n;${text}\n;\n`;
}

export function serializeStructureRecords(atoms, metadata = null) {
  const lines = ['data_roami_displayed', '#'];
  // Keep ligand chemistry dictionaries when creating transformed/selected-model input.
  for (const table of metadata?.tables || []) {
    if (!/^_chem_comp(?:\.|_atom\.|_bond\.)/.test(table.columns[0] || '')) continue;
    lines.push('loop_',...table.columns);
    for (const row of table.rows) lines.push(row.map(cifValue).join(' '));
    lines.push('#');
  }
  for (const [key,value] of metadata?.items || []) {
    if (/^_chem_comp(?:\.|_atom\.|_bond\.)/.test(key)) lines.push(key+' '+cifValue(value));
  }
  const tags = ['group_PDB','id','type_symbol','label_atom_id','label_alt_id','label_comp_id','label_asym_id',
    'label_entity_id','label_seq_id','pdbx_PDB_ins_code','Cartn_x','Cartn_y','Cartn_z','occupancy','B_iso_or_equiv',
    'auth_seq_id','auth_comp_id','auth_asym_id','auth_atom_id','pdbx_PDB_model_num'];
  lines.push('loop_',...tags.map(tag=>'_atom_site.'+tag));
  const entityIds = new Map();
  const residueIndices = new Map();
  for (const [index,atom] of atoms.entries()) {
    if (!entityIds.has(atom.chain)) { entityIds.set(atom.chain,String(entityIds.size+1)); residueIndices.set(atom.chain,new Map()); }
    const residues = residueIndices.get(atom.chain);
    if (!residues.has(atom.seq)) residues.set(atom.seq,String(residues.size+1));
    const match = String(atom.seq || '').match(/^([+-]?\d+)(.*)$/);
    const seq = match ? match[1] : String(atom.seq || '1');
    const ins = match ? match[2] : '';
    const values = [atom.mmcifGroupPdb || 'ATOM',index+1,atom.element,atom.atomName,'.',atom.resName,atom.chain,
      entityIds.get(atom.chain),residues.get(atom.seq),ins,atom.x,atom.y,atom.z,atom.occupancy ?? 1,atom.bFactor,
      seq,atom.resName,atom.chain,atom.atomName,1];
    lines.push(values.map(cifValue).join(' '));
  }
  lines.push('#');
  return lines.join('\n')+'\n';
}
