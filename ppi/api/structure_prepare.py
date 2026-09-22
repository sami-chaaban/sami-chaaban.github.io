"""Prepare shared metadata and chain rows once in a disposable parser process.

The manifest references fragments instead of duplicating large CIF metadata for
every chain. Gemmi raw tokens preserve quoting, multiline values and models.
"""
from __future__ import annotations

from collections import OrderedDict
import hashlib
import json
from pathlib import Path
import sys


class FragmentWriter:
    """Bound both disk writes and buffered open files while streaming atom rows."""

    def __init__(self, destination, max_bytes):
        self.destination = destination
        self.max_bytes = max_bytes
        self.written = 0
        self.handles = OrderedDict()

    def append(self, name, text):
        data = text.encode('utf-8')
        if self.written + len(data) > self.max_bytes:
            raise ValueError('Prepared structure exceeds the disk storage limit')
        stream = self.handles.pop(name, None)
        if stream is None:
            if len(self.handles) >= 32:
                _, oldest = self.handles.popitem(last=False)
                oldest.close()
            stream = (self.destination / name).open('ab', buffering=16 * 1024)
        self.handles[name] = stream
        stream.write(data)
        self.written += len(data)

    def close(self):
        while self.handles:
            _, stream = self.handles.popitem()
            stream.close()


def chain_token(chain):
    return hashlib.sha256(chain.encode('utf-8')).hexdigest()


def prepare_cif(source, writer):
    import gemmi
    document = gemmi.cif.read_file(str(source))
    blocks, all_chains = [], set()
    for index, block in enumerate(document):
        table = block.find_mmcif_category('_atom_site.')
        tags = list(table.tags) if table else []
        lowered = [tag.lower() for tag in tags]
        columns = [lowered.index(tag) for tag in
                   ('_atom_site.auth_asym_id', '_atom_site.label_asym_id') if tag in lowered]
        counts, fragments = {}, {}
        if table and not columns:
            raise ValueError('Structure atom table has no chain identifiers')
        for row in table or []:
            chains = {gemmi.cif.as_string(row[column]) for column in columns}
            tokens = list(row)
            # Semicolon fields, including the closing delimiter, need their own
            # physical lines. Other raw tokens retain Gemmi's original quoting.
            row_text = ('\n'.join(tokens) if any(t.startswith(';') for t in tokens)
                        else ' '.join(tokens)) + '\n'
            for chain in chains:
                fragment = f'{index}-{chain_token(chain)}.rows.cif'
                if chain not in counts:
                    writer.append(fragment, 'loop_\n' + '\n'.join(tags) + '\n')
                    fragments[chain] = fragment
                writer.append(fragment, row_text)
                counts[chain] = counts.get(chain, 0) + 1
                all_chains.add(chain)
        if table:
            table.erase()
        metadata = f'{index}.metadata.cif'
        writer.append(metadata, block.as_string() + '\n')
        blocks.append((metadata, fragments, counts))
    manifest = {}
    for chain in sorted(all_chains):
        parts, atoms = [], 0
        for metadata, fragments, counts in blocks:
            parts.append(metadata)
            if chain in fragments:
                parts.append(fragments[chain])
                atoms += counts[chain]
        manifest[chain] = {'parts': parts, 'atoms': atoms}
    return manifest


def prepare_pdb(source, writer):
    # Preserve record order, including MODEL/ENDMDL and metadata between atom
    # sections. Each shared metadata run is written once, independent of chains.
    segments, atom_counts = [], {}
    with Path(source).open() as stream:
        for original in stream:
            line = original.rstrip('\r\n') + '\n'
            record = line[:6].strip().upper()
            if record in {'CONECT', 'MASTER'}:
                continue
            coordinates = record in {'ATOM', 'HETATM', 'ANISOU', 'TER'}
            kind = 'coordinates' if coordinates else 'metadata'
            if not segments or segments[-1]['kind'] != kind:
                segments.append({'kind': kind, 'fragments': {}, 'last': {}})
            segment = segments[-1]
            if coordinates:
                chain = line[21:22].strip()
                name = f'{len(segments)-1}-{chain_token(chain)}.rows.pdb'
                segment['fragments'][chain] = name
                segment['last'][chain] = line
                if record in {'ATOM', 'HETATM'}:
                    atom_counts[chain] = atom_counts.get(chain, 0) + 1
            else:
                name = f'{len(segments)-1}.metadata.pdb'
                segment['metadata'] = name
                segment['last_line'] = line
            writer.append(name, line)
    manifest, needs_terminator = {}, False
    for chain, count in atom_counts.items():
        parts, last_line = [], ''
        for segment in segments:
            if segment['kind'] == 'metadata':
                parts.append(segment['metadata'])
                last_line = segment['last_line']
            elif chain in segment['fragments']:
                parts.append(segment['fragments'][chain])
                last_line = segment['last'][chain]
        if last_line.strip().upper() != 'END':
            parts.append('terminator.pdb')
            needs_terminator = True
        manifest[chain] = {'parts': parts, 'atoms': count}
    if needs_terminator:
        writer.append('terminator.pdb', 'END\n')
    return manifest


def prepare(source, fmt, destination, max_bytes):
    destination = Path(destination)
    destination.mkdir()
    writer = FragmentWriter(destination, max_bytes)
    try:
        if fmt == 'mmcif':
            manifest = prepare_cif(source, writer)
        elif fmt == 'pdb':
            manifest = prepare_pdb(source, writer)
        else:
            raise ValueError('Unsupported structure format')
        if not manifest:
            raise ValueError('Structure contains no atom rows')
        writer.append('manifest.json', json.dumps(manifest, separators=(',', ':')))
        return manifest
    finally:
        writer.close()


if __name__ == '__main__':
    try:
        prepare(*sys.argv[1:4], int(sys.argv[4]))
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
