"""Stored-coordinate lifetime, lossless preparation and API source identity."""
from __future__ import annotations

import asyncio
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gemmi
from fastapi import HTTPException
from starlette.requests import Request
from api import main
from api.analysis_worker import BoundedAnalysisRunner, prepare_delivery
from api.cache import ReportCache
from api.models import AnalyzeRequest, ChapiMeshRequest
from api.structure_prepare import prepare
from api.structure_store import StructureCapacity, StructureGone, StructureStore


CIF = """data_first
_entry.id example
_struct.title 'A quoted title'
loop_
_atom_site.auth_asym_id
_atom_site.label_asym_id
_atom_site.auth_seq_id
_atom_site.label_atom_id
_atom_site.pdbx_PDB_model_num
_atom_site.details
'A chain' X 42 "O5'" 1 'one'
B 'A chain' 7 N 1
;multiline
detail
;
C C 8 CA 1 .
'A chain' X 42 "O5'" 2 ?
#
data_second
_entry.id second
_struct.title
;second block
title
;
loop_
_atom_site.AUTH_ASYM_ID
_atom_site.label_asym_id
_atom_site.auth_seq_id
_atom_site.label_atom_id
_atom_site.pdbx_PDB_model_num
Z X 4 CA 1
C C 5 CA 2
#
"""


def store_bytes(store, data, fmt='mmcif'):
    path = store.upload_path()
    path.write_bytes(data)
    try:
        return store.register(path, fmt, hashlib.sha256((fmt + '\n').encode() + data).hexdigest())
    finally:
        path.unlink(missing_ok=True)


def cif_semantics(text):
    output = []
    for block in gemmi.cif.read_string(text):
        categories = {}
        for item in block:
            if item.pair is not None:
                tag, value = item.pair
                categories[tag.lower()] = gemmi.cif.as_string(value)
            elif item.loop is not None and item.loop.values:
                width = len(item.loop.tags)
                for index, tag in enumerate(item.loop.tags):
                    categories[tag.lower()] = [gemmi.cif.as_string(value)
                                               for value in item.loop.values[index::width]]
        output.append((block.name, categories))
    return output


class StoreTests(unittest.TestCase):
    def make_store(self, **options):
        store = StructureStore(**options)
        self.addCleanup(store.close)
        return store

    def test_expiry_and_explicit_delete_preserve_active_lease_until_last_release(self):
        store = self.make_store(ttl_seconds=10)
        with patch('api.structure_store.time.monotonic', return_value=0):
            idle = store_bytes(store, b'idle')
            active = store_bytes(store, b'active')
            first, second = store.acquire(active.identifier), store.acquire(active.identifier)
        with patch('api.structure_store.time.monotonic', return_value=11):
            with self.assertRaises(StructureGone):
                store.acquire(idle.identifier)
            self.assertFalse(idle.path.exists())
            self.assertTrue(active.path.exists())
            store.delete(active.identifier)
            with self.assertRaises(StructureGone):
                store.acquire(active.identifier)
            first.release()
            first.release()
            self.assertTrue(active.path.exists())
            self.assertEqual(active.readers, 1)
            second.release()
            self.assertFalse(active.directory.exists())

    def test_capacity_evicts_unleased_lru_and_never_deletes_active_input(self):
        store = self.make_store(max_bytes=10, max_entries=2)
        first = store_bytes(store, b'1111')
        second = store_bytes(store, b'2222')
        lease = store.acquire(first.identifier)
        third = store_bytes(store, b'333333')
        self.assertTrue(first.path.exists())
        self.assertFalse(second.path.exists())
        third_lease = store.acquire(third.identifier)
        with self.assertRaises(StructureCapacity):
            store_bytes(store, b'4')
        self.assertTrue(first.path.exists())
        self.assertTrue(third.path.exists())
        lease.release()
        third_lease.release()

    def test_oversized_upload_does_not_evict_existing_structures(self):
        store = self.make_store(max_bytes=5, max_upload_bytes=20)
        existing = store_bytes(store, b'1234')
        with self.assertRaises(StructureCapacity):
            store_bytes(store, b'123456')
        self.assertTrue(existing.path.exists())

    def test_handles_cannot_be_used_as_filesystem_paths(self):
        store = self.make_store()
        entry = store_bytes(store, b'data_example')
        self.assertNotIn('/', entry.identifier)
        for identifier in ('../source.cif', str(entry.path), str(store.root), '../' + entry.identifier):
            with self.assertRaises(StructureGone):
                store.acquire(identifier)
        self.assertTrue(entry.path.exists())

    def test_chain_preparation_is_reused_and_counts_against_capacity(self):
        store = self.make_store(max_bytes=100000)
        entry = store_bytes(store, CIF.encode())
        with store.acquire(entry.identifier):
            path, atoms = store.chain_path(entry, 'A chain', sys.executable, os.environ.copy())
            self.assertEqual(atoms, 3)
            self.assertTrue(path.exists())
            self.assertGreater(entry.size, len(CIF.encode()))
            self.assertEqual(entry.size, sum(p.stat().st_size for p in entry.directory.rglob('*') if p.is_file()))
            self.assertIn(path, entry.activefiles)
            active_size = path.stat().st_size
            cached_size = entry.size - active_size
            self.store_release_chain(store, entry, path)
            self.assertEqual(entry.size, cached_size)
            with patch('api.structure_store.subprocess.run') as launch:
                second_path, second_atoms = store.chain_path(entry, 'A chain', sys.executable, os.environ.copy())
                launch.assert_not_called()
            self.assertNotEqual(second_path, path)
            self.assertEqual(second_atoms, atoms)
            self.store_release_chain(store, entry, second_path)
            self.assertEqual(entry.size, cached_size)
            self.assertEqual(store.chain_path(entry, '../../x', sys.executable, os.environ.copy()), (None, 0))

    def store_release_chain(self, store, entry, path):
        store.release_chain(entry, path)
        store.release_chain(entry, path)
        self.assertFalse(path.exists())
        self.assertNotIn(path, entry.activefiles)

    def test_failed_preparation_removes_partial_files_and_keeps_source(self):
        data = CIF.encode()
        store = self.make_store(max_bytes=len(data) + 10)
        entry = store_bytes(store, data)
        with store.acquire(entry.identifier):
            with self.assertRaises(ValueError):
                store.chain_path(entry, 'A chain', sys.executable, os.environ.copy())
        self.assertTrue(entry.path.exists())
        self.assertFalse((entry.directory / 'chains').exists())
        self.assertEqual(entry.size, len(data))
        self.assertIsNone(entry.chains)

    def test_metadata_is_shared_under_budget_and_repeated_chain_files_are_released(self):
        text = ('data_large_header\n_struct.title\n;' + 'metadata ' * 2000 + '\n;\n'
                + 'loop_\n_atom_site.auth_asym_id\n_atom_site.label_asym_id\n_atom_site.id\n'
                + ''.join(f'C{index} C{index} {index}\n' for index in range(40)))
        store = self.make_store(max_bytes=100000)
        entry = store_bytes(store, text.encode())
        with store.acquire(entry.identifier):
            retained = None
            for index in range(40):
                path, atoms = store.chain_path(entry, f'C{index}', sys.executable, os.environ.copy())
                self.assertEqual(atoms, 1)
                self.assertLessEqual(entry.size, store.max_bytes)
                store.release_chain(entry, path)
                if retained is None:
                    retained = entry.size
                self.assertEqual(entry.size, retained)
            self.assertEqual(len(list((entry.directory / 'chains').glob('*.metadata.cif'))), 1)
            self.assertEqual(entry.activefiles, {})
            self.assertLess(entry.size, 3 * len(text.encode()))

    def test_materialization_reserves_bytes_and_cleans_partial_copy_failure(self):
        store = self.make_store(max_bytes=100000)
        entry = store_bytes(store, CIF.encode())
        with store.acquire(entry.identifier):
            path, _ = store.chain_path(entry, 'X', sys.executable, os.environ.copy())
            active_size = path.stat().st_size
            store.release_chain(entry, path)
            retained = entry.size
            def failing_copy(fragment, output, length):
                self.assertEqual(entry.size, retained + active_size)
                self.assertEqual(length, 65536)
                output.write(fragment.read(5))
                raise OSError('simulated disk failure')
            with patch('api.structure_store.shutil.copyfileobj', side_effect=failing_copy):
                with self.assertRaisesRegex(OSError, 'simulated disk failure'):
                    store.chain_path(entry, 'X', sys.executable, os.environ.copy())
            self.assertEqual(entry.size, retained)
            self.assertEqual(entry.activefiles, {})
            self.assertFalse(list(entry.directory.glob('active-*')))
            store.max_bytes = retained + active_size - 1
            with self.assertRaises(StructureCapacity):
                store.chain_path(entry, 'X', sys.executable, os.environ.copy())
            self.assertEqual(entry.size, retained)
            self.assertEqual(entry.activefiles, {})

    def test_active_chain_file_survives_delete_until_explicit_release(self):
        store = self.make_store(max_bytes=100000)
        entry = store_bytes(store, CIF.encode())
        lease = store.acquire(entry.identifier)
        path, _ = store.chain_path(entry, 'X', sys.executable, os.environ.copy())
        store.delete(entry.identifier)
        lease.release()
        self.assertTrue(path.exists())
        with self.assertRaises(StructureGone):
            store.acquire(entry.identifier)
        store.release_chain(entry, path)
        self.assertFalse(entry.directory.exists())


class PrepareTests(unittest.TestCase):
    def test_actual_gemmi_matches_existing_filter_for_aliases_multiline_blocks_models(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'source.cif'
            source.write_text(CIF)
            output = Path(directory) / 'chains'
            manifest = prepare(source, 'mmcif', output, 100000)
            for chain in ('A chain', 'X', 'B', 'C', 'Z'):
                with self.subTest(chain=chain):
                    expected, count = main.filter_mmcif_text_to_single_chain(CIF, chain)
                    actual = ''.join((output / part).read_text() for part in manifest[chain]['parts'])
                    self.assertEqual(manifest[chain]['atoms'], count)
                    self.assertEqual(cif_semantics(actual), cif_semantics(expected))
                    self.assertEqual(len(gemmi.cif.read_string(actual)), 2)
            self.assertEqual(len(list(output.glob('*.metadata.cif'))), 2)

    def test_pdb_filter_parity_keeps_model_records_and_removes_stale_connectivity(self):
        def atom(record, number, chain):
            return f'{record:<6}{number:5d}  CA  ALA {chain}   1       1.000   2.000   3.000\n'
        pdb = ('HEADER    test\nMODEL        1\n' + atom('ATOM', 1, 'A') + atom('ANISOU', 1, 'A')
               + atom('ATOM', 2, 'B') + atom('TER', 3, 'B') + 'ENDMDL\nMODEL        2\n'
               + atom('HETATM', 4, 'A') + 'ENDMDL\nCONECT    1    2\nMASTER        4\n')
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'source.pdb'
            source.write_text(pdb)
            output = Path(directory) / 'chains'
            manifest = prepare(source, 'pdb', output, 100000)
            for chain in ('A', 'B'):
                expected, count = main.filter_pdb_text_to_single_chain(pdb, chain)
                actual = ''.join((output / part).read_text() for part in manifest[chain]['parts'])
                self.assertEqual(actual, expected)
                self.assertEqual(manifest[chain]['atoms'], count)


class Client:
    disconnected = False
    headers = {}

    async def is_disconnected(self):
        return self.disconnected


async def until(predicate):
    async def poll():
        while not predicate():
            await asyncio.sleep(0.005)
    await asyncio.wait_for(poll(), 5)


def stream_request(chunks, headers=()):
    remaining = iter(chunks)
    async def receive():
        try:
            body = next(remaining)
            return {'type': 'http.request', 'body': body, 'more_body': True}
        except StopIteration:
            return {'type': 'http.request', 'body': b'', 'more_body': False}
    return Request({'type': 'http', 'method': 'POST', 'path': '/structures', 'headers': list(headers)}, receive)


class StructureApiTests(unittest.IsolatedAsyncioTestCase):
    def setUp(self):
        self.store = StructureStore(max_upload_bytes=4096)
        self.addCleanup(self.store.close)
        self.context = ExitStack()
        self.addCleanup(self.context.close)
        self.context.enter_context(patch.object(main, 'structure_store', self.store))
        self.context.enter_context(patch.object(main, 'STRUCTURE_UPLOAD_SLOTS', asyncio.Semaphore(2)))
        self.context.enter_context(patch.object(main, 'cache', ReportCache()))
        self.context.enter_context(patch.object(main, 'report_store', ReportCache()))

    async def test_upload_preserves_split_utf8_and_content_identity(self):
        text = 'data_test\n_struct.title "α protein"\n'
        encoded = text.encode()
        split = encoded.index('α'.encode()) + 1
        response = await main.upload_structure(stream_request([encoded[:split], encoded[split:]]), 'mmcif')
        entry = self.store.entries[response['structureId']]
        self.assertEqual(entry.path.read_bytes(), encoded)
        self.assertEqual(entry.digest, hashlib.sha256(b'mmcif\n' + encoded).hexdigest())
        self.assertEqual(set(response), {'structureId', 'format', 'expiresIn'})
        self.assertFalse(list(self.store.root.glob('upload-*')))

    async def test_actual_upload_and_delete_routes_accept_raw_streamed_coordinates(self):
        async def call(method, path, body=(), query=b''):
            messages = []
            request = stream_request(body)
            scope = {'type': 'http', 'asgi': {'version': '3.0'}, 'http_version': '1.1',
                     'method': method, 'scheme': 'http', 'path': path, 'raw_path': path.encode(),
                     'query_string': query, 'root_path': '', 'headers': [(b'accept-encoding', b'identity')],
                     'server': ('test', 80), 'client': ('test', 1)}
            async def send(message):
                messages.append(message)
            await main.app(scope, request.receive, send)
            status = next(item['status'] for item in messages if item['type'] == 'http.response.start')
            data = b''.join(item.get('body', b'') for item in messages if item['type'] == 'http.response.body')
            return status, data
        status, body = await call('POST', '/structures', [b'HEAD', b'ER    raw\nEND\n'], b'format=pdb')
        self.assertEqual(status, 200)
        identifier = json.loads(body)['structureId']
        entry = self.store.entries[identifier]
        self.assertEqual(entry.format, 'pdb')
        self.assertEqual(entry.path.read_bytes(), b'HEADER    raw\nEND\n')
        status, body = await call('DELETE', '/structures/' + identifier)
        self.assertEqual((status, body), (204, b''))
        self.assertFalse(entry.directory.exists())

    async def test_invalid_streams_and_size_limits_leave_no_uploads_or_entries(self):
        self.store.max_upload_bytes = 4
        cases = [([], [], 400), ([b'abc', b'de'], [], 413),
                 ([b'abc'], [(b'content-length', b'5')], 413),
                 ([b'\xff'], [], 400), ([b'\xce'], [], 400)]
        for chunks, headers, status in cases:
            with self.subTest(chunks=chunks, headers=headers):
                with self.assertRaises(HTTPException) as error:
                    await main.upload_structure(stream_request(chunks, headers), 'mmcif')
                self.assertEqual(error.exception.status_code, status)
                self.assertFalse(list(self.store.root.iterdir()))
        with self.assertRaises(HTTPException) as error:
            await main.upload_structure(stream_request([b'x']), 'other')
        self.assertEqual(error.exception.status_code, 400)

    async def test_stored_and_inline_analyses_share_content_key_without_confusing_changed_coordinates(self):
        first = store_bytes(self.store, b'INPUT1', 'pdb')
        duplicate = store_bytes(self.store, b'INPUT1', 'pdb')
        changed = store_bytes(self.store, b'INPUT2', 'pdb')
        calls = []
        def native(path, *args):
            text = Path(path).read_text()
            calls.append(text)
            return prepare_delivery({'contacts': {}, 'input': text})
        runner = BoundedAnalysisRunner(executor=ThreadPoolExecutor(max_workers=1))
        self.addCleanup(runner.shutdown)
        with patch.object(main, 'analysis_runner', runner), patch.object(main, 'analyze_stored_and_serialize', native):
            responses = []
            for entry in (first, duplicate, changed):
                response = await main.analyze(AnalyzeRequest(structureId=entry.identifier, pdbId='1abc', chainA='A', chainB='B'), Client())
                responses.append(json.loads(response.body))
                self.assertEqual(entry.readers, 0)
            inline = await main.analyze(AnalyzeRequest(pdbText='INPUT1', pdbId='1abc', chainA='A', chainB='B'), Client())
        self.assertEqual(calls, ['INPUT1', 'INPUT2'])
        self.assertEqual(responses[0]['reportId'], responses[1]['reportId'])
        self.assertEqual(responses[0]['reportId'], json.loads(inline.body)['reportId'])
        self.assertNotEqual(responses[0]['reportId'], responses[2]['reportId'])

    async def test_deleted_analysis_input_remains_until_disconnected_worker_finishes(self):
        entry = store_bytes(self.store, b'INPUT', 'pdb')
        started, release = threading.Event(), threading.Event()
        def native(path, *args):
            started.set()
            release.wait(3)
            return prepare_delivery({'contacts': {}, 'input': Path(path).read_text()})
        runner = BoundedAnalysisRunner(executor=ThreadPoolExecutor(max_workers=1))
        self.addCleanup(runner.shutdown)
        client = Client()
        with patch.object(main, 'analysis_runner', runner), patch.object(main, 'analyze_stored_and_serialize', native):
            request = asyncio.create_task(main.analyze(AnalyzeRequest(structureId=entry.identifier, chainA='A', chainB='B'), client))
            await until(started.is_set)
            self.store.delete(entry.identifier)
            client.disconnected = True
            try:
                with self.assertRaises(HTTPException) as error:
                    await request
                self.assertEqual(error.exception.status_code, 499)
                self.assertEqual(entry.readers, 1)
                self.assertTrue(entry.path.exists())
            finally:
                release.set()
                await until(lambda: runner.pending == 0)
            self.assertFalse(entry.directory.exists())

    async def test_cancelled_analysis_waiting_for_gate_releases_both_input_leases(self):
        entry = store_bytes(self.store, b'INPUT', 'pdb')
        runner = BoundedAnalysisRunner(execution_gate=threading.Semaphore(0))
        with patch.object(main, 'analysis_runner', runner), patch.object(main, 'analyze_stored_and_serialize') as native:
            request = asyncio.create_task(main.analyze(
                AnalyzeRequest(structureId=entry.identifier, chainA='A', chainB='B'), Client()))
            await until(lambda: entry.readers == 2 and runner.slots.locked())
            request.cancel()
            with self.assertRaises(asyncio.CancelledError):
                await request
            self.assertEqual(entry.readers, 0)
            self.assertEqual(runner.pending, 0)
            native.assert_not_called()
            self.store.delete(entry.identifier)
            self.assertFalse(entry.directory.exists())

    async def test_stored_mesh_uses_prepared_path_and_original_digest(self):
        entry = store_bytes(self.store, CIF.encode())
        observed = []
        def mesh(payload, key):
            self.assertIsNone(payload['text'])
            path = Path(payload['sourcePath'])
            self.assertNotEqual(path, entry.path)
            expected, _ = main.filter_mmcif_text_to_single_chain(CIF, 'X')
            self.assertEqual(cif_semantics(path.read_text()), cif_semantics(expected))
            observed.append((payload['_reader_source_key'], key))
            return b'{"meshes":[]}'
        with patch.object(main, 'CHAPI_LOW_MEMORY_MODE', True), patch.object(main, 'CHAPI_SPLIT_CHAIN_BATCH_LIMIT', 1), patch.object(main, 'CHAPI_PYTHON', sys.executable), patch.object(main, '_build_chapi_bridge_env', return_value=os.environ.copy()), patch.object(main, 'run_chapi_mesh_cached', side_effect=mesh):
            request = ChapiMeshRequest(structureId=entry.identifier, representation='ribbon', splitByChain=True, chainIds=['X'])
            response = await asyncio.to_thread(main.chapi_mesh, request)
        self.assertEqual(response.status_code, 200)
        self.assertEqual(observed[0][0], entry.digest)
        self.assertEqual(entry.readers, 0)
        self.assertTrue(entry.path.exists())
        self.assertEqual(entry.activefiles, {})
        self.assertFalse(list(entry.directory.glob('active-*')))

    async def test_native_mesh_failure_releases_materialized_chain_and_input_lease(self):
        entry = store_bytes(self.store, CIF.encode())
        with patch.object(main, 'CHAPI_LOW_MEMORY_MODE', True), patch.object(main, 'CHAPI_SPLIT_CHAIN_BATCH_LIMIT', 1), patch.object(main, 'CHAPI_PYTHON', sys.executable), patch.object(main, '_build_chapi_bridge_env', return_value=os.environ.copy()), patch.object(main, 'run_chapi_mesh_cached', side_effect=RuntimeError('native failure')):
            request = ChapiMeshRequest(structureId=entry.identifier, representation='ribbon', splitByChain=True, chainIds=['X'])
            with self.assertRaises(HTTPException) as error:
                await asyncio.to_thread(main.chapi_mesh, request)
        self.assertEqual(error.exception.status_code, 500)
        self.assertEqual(entry.readers, 0)
        self.assertEqual(entry.activefiles, {})
        self.assertFalse(list(entry.directory.glob('active-*')))
        self.assertTrue(entry.path.exists())

    async def test_stored_mesh_cache_key_uses_content_not_handle_or_pdb_name(self):
        first = store_bytes(self.store, CIF.encode())
        duplicate = store_bytes(self.store, CIF.encode())
        changed = store_bytes(self.store, CIF.replace('42 "O5\'" 1', '43 "O5\'" 1').encode())
        keys = []
        def mesh(payload, key):
            keys.append(key)
            self.assertIsNone(payload['text'])
            self.assertTrue(Path(payload['sourcePath']).exists())
            return b'{"meshes":[]}'
        with patch.object(main, 'CHAPI_LOW_MEMORY_MODE', False), patch.object(main, 'run_chapi_mesh_cached', side_effect=mesh):
            for entry in (first, duplicate, changed):
                await asyncio.to_thread(main.chapi_mesh, ChapiMeshRequest(
                    structureId=entry.identifier, pdbId='7z8g', representation='bonds'))
        self.assertEqual(keys[0], keys[1])
        self.assertNotEqual(keys[0], keys[2])

    async def test_deleted_and_pathlike_handles_return_410(self):
        entry = store_bytes(self.store, CIF.encode())
        self.store.delete(entry.identifier)
        for identifier in (entry.identifier, str(entry.path), '../anything'):
            with self.assertRaises(HTTPException) as error:
                await main.analyze(AnalyzeRequest(structureId=identifier, chainA='A', chainB='B'), Client())
            self.assertEqual(error.exception.status_code, 410)
            with self.assertRaises(HTTPException) as error:
                main.chapi_mesh(ChapiMeshRequest(structureId=identifier))
            self.assertEqual(error.exception.status_code, 410)


if __name__ == '__main__':
    unittest.main()
