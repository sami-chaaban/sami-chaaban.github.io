"""Bounded, expiring disk storage for uploaded coordinates and prepared chains.

Only opaque, unguessable handles cross the HTTP boundary. Leases keep files alive
until their native consumer exits, including after a browser disconnects.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass, field
import json
from pathlib import Path
import secrets
import shutil
import subprocess
import tempfile
import threading
import time


class StructureGone(ValueError):
    pass


class StructureCapacity(ValueError):
    pass


@dataclass
class StoredStructure:
    identifier: str
    directory: Path
    format: str
    digest: str
    size: int
    touched: float
    readers: int = 0
    deleted: bool = False
    chains: dict | None = None
    activefiles: dict[Path, int] = field(default_factory=dict)
    prepare_lock: threading.Lock = field(default_factory=threading.Lock)

    @property
    def path(self):
        return self.directory / ('source.cif' if self.format == 'mmcif' else 'source.pdb')


class StructureLease:
    def __init__(self, store, entry):
        self.store, self.entry, self.released = store, entry, False

    def release(self):
        with self.store.lock:
            if not self.released:
                self.released = True
                self.entry.readers -= 1
                if self.entry.deleted and not self.entry.readers and not self.entry.activefiles:
                    self.store._remove(self.entry.identifier)

    def __enter__(self):
        return self.entry

    def __exit__(self, *_):
        self.release()


class StructureStore:
    def __init__(self, max_bytes=256 * 1024**2, max_entries=8, ttl_seconds=1800,
                 max_upload_bytes=64 * 1024**2, root=None):
        self.max_bytes, self.max_entries = max_bytes, max_entries
        self.ttl_seconds, self.max_upload_bytes = ttl_seconds, max_upload_bytes
        self._temporary = tempfile.TemporaryDirectory(prefix='roami-structures-', dir=root)
        self.root = Path(self._temporary.name)
        self.entries = OrderedDict()
        self.lock = threading.RLock()

    def upload_path(self):
        descriptor, name = tempfile.mkstemp(prefix='upload-', dir=self.root)
        import os
        os.close(descriptor)
        return Path(name)

    def _remove(self, identifier):
        entry = self.entries.pop(identifier, None)
        if entry:
            shutil.rmtree(entry.directory, ignore_errors=True)

    def _prune(self):
        now = time.monotonic()
        for identifier, entry in list(self.entries.items()):
            if not entry.readers and not entry.activefiles and (entry.deleted or now - entry.touched > self.ttl_seconds):
                self._remove(identifier)

    def _make_room(self, extra_bytes, extra_entries=0, protected=None):
        if extra_bytes > self.max_bytes or extra_entries > self.max_entries:
            raise StructureCapacity('Structure exceeds the disk storage limit.')
        self._prune()
        while (sum(e.size for e in self.entries.values()) + extra_bytes > self.max_bytes
               or len(self.entries) + extra_entries > self.max_entries):
            victim = next((e for e in self.entries.values()
                           if not e.readers and not e.activefiles and e.identifier != protected), None)
            if victim is None:
                raise StructureCapacity('Structure storage is busy or full. Please retry shortly.')
            self._remove(victim.identifier)

    def register(self, path, fmt, digest):
        size = path.stat().st_size
        if not size or size > self.max_upload_bytes:
            raise StructureCapacity('Structure upload is empty or exceeds the upload limit.')
        with self.lock:
            self._make_room(size, 1)
            identifier = secrets.token_urlsafe(24)
            directory = self.root / identifier
            directory.mkdir(mode=0o700)
            entry = StoredStructure(identifier, directory, fmt, digest, size, time.monotonic())
            path.replace(entry.path)
            self.entries[identifier] = entry
            return entry

    def acquire(self, identifier):
        with self.lock:
            self._prune()
            entry = self.entries.get(identifier)
            if entry is None or entry.deleted:
                raise StructureGone('Structure expired. Upload the coordinates again.')
            entry.readers += 1
            entry.touched = time.monotonic()
            self.entries.move_to_end(identifier)
            return StructureLease(self, entry)

    def delete(self, identifier):
        with self.lock:
            entry = self.entries.get(identifier)
            if entry:
                entry.deleted = True
                if not entry.readers and not entry.activefiles:
                    self._remove(identifier)

    def chain_path(self, entry, chain_id, python, env, timeout=180):
        """Prepare shared fragments once; materialize a leased temporary chain.

        The caller must release_chain() after native consumption, before
        releasing its entry lease.
        """
        with entry.prepare_lock:
            if entry.chains is None:
                destination = entry.directory / 'chains'
                allowance = 0
                reserved = False
                try:
                    # Reserve before the child writes, so concurrent uploads or
                    # preparations cannot spend the same remaining allowance.
                    with self.lock:
                        self._make_room(0, protected=entry.identifier)
                        allowance = self.max_bytes - sum(e.size for e in self.entries.values())
                        entry.size += allowance
                        reserved = True
                    result = subprocess.run(
                        [python, str(Path(__file__).with_name('structure_prepare.py')),
                         str(entry.path), entry.format, str(destination), str(allowance)],
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env, timeout=timeout,
                    )
                    if result.returncode:
                        raise ValueError(result.stderr.decode('utf-8', errors='replace').strip()
                                         or 'Structure preparation failed')
                    manifest = json.loads((destination / 'manifest.json').read_text())
                    size = sum(p.stat().st_size for p in destination.iterdir() if p.is_file())
                    with self.lock:
                        if size > allowance:
                            raise StructureCapacity('Prepared structure exceeds the disk storage limit.')
                        entry.size += size - allowance
                        reserved = False
                        entry.chains = manifest
                except BaseException:
                    shutil.rmtree(destination, ignore_errors=True)
                    if reserved:
                        with self.lock:
                            entry.size -= allowance
                    raise
            info = entry.chains.get(chain_id)
            if info is None:
                return None, 0
            parts = [entry.directory / 'chains' / name for name in info['parts']]
            size = sum(part.stat().st_size for part in parts)
            suffix = '.cif' if entry.format == 'mmcif' else '.pdb'
            path = entry.directory / ('active-' + secrets.token_urlsafe(18) + suffix)
            with self.lock:
                self._make_room(size, protected=entry.identifier)
                entry.size += size
                entry.activefiles[path] = size
            try:
                with path.open('xb') as output:
                    for part in parts:
                        with part.open('rb') as fragment:
                            shutil.copyfileobj(fragment, output, length=65536)
                return path, info['atoms']
            except BaseException:
                self.release_chain(entry, path)
                raise

    def release_chain(self, entry, path):
        """Delete a materialized chain and return its reserved bytes exactly once."""
        path = Path(path)
        with self.lock:
            size = entry.activefiles.get(path)
            if size is None:
                return
            path.unlink(missing_ok=True)
            entry.activefiles.pop(path)
            entry.size -= size
            if entry.deleted and not entry.readers and not entry.activefiles:
                self._remove(entry.identifier)

    def close(self):
        self._temporary.cleanup()
