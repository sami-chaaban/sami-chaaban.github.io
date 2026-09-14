"""Mesh delivery regressions; no native Coot or network access required."""
from __future__ import annotations

import asyncio
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
import io
import signal
import subprocess
from pathlib import Path
import sys
import threading
import time
import unittest
from unittest.mock import Mock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from api import main


async def request(method: str, path: str, payload=None, encoding="identity"):
    body = json.dumps(payload).encode() if payload is not None else b""
    messages = []
    sent = False

    async def receive():
        nonlocal sent
        if not sent:
            sent = True
            return {"type": "http.request", "body": body, "more_body": False}
        return {"type": "http.disconnect"}

    async def send(message):
        messages.append(message)

    await main.app({
        "type": "http", "asgi": {"version": "3.0"}, "http_version": "1.1",
        "method": method, "scheme": "http", "path": path, "raw_path": path.encode(),
        "query_string": b"", "root_path": "", "server": ("localhost", 80),
        "client": ("127.0.0.1", 12345),
        "headers": [(b"content-type", b"application/json"), (b"accept-encoding", encoding.encode())],
    }, receive, send)
    start = next(message for message in messages if message["type"] == "http.response.start")
    content = b"".join(message.get("body", b"") for message in messages if message["type"] == "http.response.body")
    return start["status"], dict(start["headers"]), content


PAYLOAD = {"mmcifText": "data_test\n", "representation": "ribbon", "splitByChain": True, "chainIds": ["A"]}


class MeshDeliveryTests(unittest.IsolatedAsyncioTestCase):
    async def test_mesh_compression_is_lossless_and_negotiated(self):
        mesh_bytes = json.dumps({"positions": [1.23456, -7.891, 0.00001] * 10000}, separators=(",", ":")).encode()
        with patch.object(main, "run_chapi_mesh_cached", return_value=mesh_bytes):
            status, headers, body = await request("POST", "/chapi-mesh", PAYLOAD, "gzip")
            self.assertEqual(status, 200)
            self.assertEqual(headers[b"content-encoding"], b"gzip")
            self.assertIn(b"Accept-Encoding", headers[b"vary"])
            self.assertEqual(gzip.decompress(body), mesh_bytes)
            self.assertLess(len(body), len(mesh_bytes) // 2)
            status, headers, body = await request("POST", "/chapi-mesh", PAYLOAD)
            self.assertEqual(status, 200)
            self.assertNotIn(b"content-encoding", headers)
            self.assertEqual(body, mesh_bytes)

    async def test_slow_mesh_does_not_block_other_http_requests(self):
        loop = asyncio.get_running_loop()
        entered = asyncio.Event()
        release = threading.Event()

        def slow_mesh(*_):
            loop.call_soon_threadsafe(entered.set)
            if not release.wait(timeout=2):
                raise RuntimeError("Mesh handler blocked the event loop")
            return b'{"meshes":[]}'

        with patch.object(main, "run_chapi_mesh_cached", side_effect=slow_mesh):
            mesh_task = asyncio.create_task(request("POST", "/chapi-mesh", PAYLOAD))
            try:
                await asyncio.wait_for(entered.wait(), timeout=1)
                status, _, _ = await asyncio.wait_for(request("GET", "/"), timeout=0.5)
                self.assertEqual(status, 200)
            finally:
                release.set()
                status, _, _ = await mesh_task
            self.assertEqual(status, 200)

    async def test_low_memory_batch_limit_is_retained(self):
        with patch.object(main, "CHAPI_SPLIT_CHAIN_BATCH_LIMIT", 1), patch.object(main, "run_chapi_mesh_cached") as generate:
            status, _, body = await request("POST", "/chapi-mesh", dict(PAYLOAD, chainIds=["A", "B"]))
            self.assertEqual(status, 413)
            self.assertEqual(json.loads(body)["detail"]["errorCode"], "CHAPI-LOWMEM-002")
            generate.assert_not_called()


class NativeMemoryBoundTests(unittest.TestCase):
    def test_concurrent_requests_keep_native_jobs_serial(self):
        active = 0
        peak = 0
        counter_lock = threading.Lock()
        start = threading.Barrier(4)

        def native(*_):
            nonlocal active, peak
            with counter_lock:
                active += 1
                peak = max(peak, active)
            try:
                time.sleep(0.02)
                return b'{"meshes":[]}'
            finally:
                with counter_lock:
                    active -= 1

        def submit(_):
            start.wait(timeout=1)
            return main.run_chapi_mesh({})

        with patch.object(main, "CHAPI_PERSISTENT_WORKER", False), patch.object(main, "_run_chapi_mesh_subprocess", side_effect=native):
            with ThreadPoolExecutor(max_workers=4) as pool:
                result = list(pool.map(submit, range(4)))
        self.assertEqual(peak, 1)
        self.assertEqual(result, [b'{"meshes":[]}'] * 4)


class NativeReaderRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.cache_patch = patch.object(main, "CHAPI_MMDB_READER_CACHE", main.ReportCache(max_entries=4))
        self.cache_patch.start()
        self.addCleanup(self.cache_patch.stop)
        self.payload = {"text": "data_reader_test", "format": "mmcif", "chainIds": ["A"]}

    def test_successful_default_reader_is_unchanged(self):
        with patch.object(main, "_dispatch_chapi_mesh", return_value=b"mesh") as dispatch:
            self.assertEqual(main.run_chapi_mesh(self.payload), b"mesh")
        self.assertEqual(dispatch.call_args.args[0], self.payload)
        self.assertNotIn("_use_mmdb_reader", self.payload)

    def test_segfault_recovers_once_and_remembers_reader_for_other_chains(self):
        with patch.object(main, "_dispatch_chapi_mesh", side_effect=[main.ChapiNativeProcessError(-signal.SIGSEGV), b"recovered", b"next-chain"]) as dispatch:
            self.assertEqual(main.run_chapi_mesh(self.payload), b"recovered")
            self.assertEqual(main.run_chapi_mesh(dict(self.payload, chainIds=["B"])), b"next-chain")
        calls = [call.args[0] for call in dispatch.call_args_list]
        self.assertNotIn("_use_mmdb_reader", calls[0])
        self.assertTrue(calls[1]["_use_mmdb_reader"])
        self.assertTrue(calls[2]["_use_mmdb_reader"])
        self.assertNotIn("_use_mmdb_reader", self.payload)

    def test_mmdb_failure_does_not_loop(self):
        with patch.object(main, "_dispatch_chapi_mesh", side_effect=main.ChapiNativeProcessError(-signal.SIGSEGV)) as dispatch:
            with self.assertRaises(main.ChapiNativeProcessError):
                main.run_chapi_mesh(self.payload)
        self.assertEqual(dispatch.call_count, 2)
        self.assertIsNone(main.CHAPI_MMDB_READER_CACHE.get(main._chapi_reader_source_key(self.payload)))

    def test_oom_timeout_and_bad_input_are_not_reader_retries(self):
        failures = [
            main.ChapiNativeProcessError(-9), main.ChapiNativeProcessError(1),
            subprocess.TimeoutExpired("coot", 180),
            main.ChapiWorkerTransportError("CHAPI worker timed out after 180s."),
            RuntimeError("Failed to read structure in chapi bridge."),
        ]
        for failure in failures:
            with self.subTest(error=str(failure)), patch.object(main, "_dispatch_chapi_mesh", side_effect=failure) as dispatch:
                with self.assertRaises(type(failure)):
                    main.run_chapi_mesh(self.payload)
                self.assertEqual(dispatch.call_count, 1)

    def test_persistent_eof_preserves_native_signal_without_same_reader_retry(self):
        proc = Mock()
        proc.stdin = io.BytesIO()
        proc.stdout = io.BytesIO()
        proc.poll.return_value = -signal.SIGSEGV
        with patch.object(main, "_start_chapi_worker_locked", return_value=proc), patch.object(main, "_stop_chapi_worker_locked"), patch.object(main.select, "select", return_value=([proc.stdout], [], [])):
            with self.assertRaises(main.ChapiNativeProcessError) as error:
                main._run_chapi_mesh_worker(self.payload)
        self.assertEqual(error.exception.returncode, -signal.SIGSEGV)

    def test_one_shot_preserves_native_signal(self):
        result = subprocess.CompletedProcess(["coot"], -signal.SIGSEGV, b"", b"")
        with patch.object(main.subprocess, "run", return_value=result):
            with self.assertRaises(main.ChapiNativeProcessError) as error:
                main._run_chapi_mesh_subprocess(self.payload)
        self.assertEqual(error.exception.returncode, -signal.SIGSEGV)

    def test_worker_timeout_does_not_start_second_process(self):
        with patch.object(main, "CHAPI_PERSISTENT_WORKER", True), patch.object(main, "_run_chapi_mesh_worker", side_effect=main.ChapiWorkerTransportError("CHAPI worker timed out after 180s.")), patch.object(main, "_run_chapi_mesh_subprocess") as oneshot:
            with self.assertRaises(main.ChapiWorkerTransportError):
                main._dispatch_chapi_mesh(self.payload)
            oneshot.assert_not_called()


if __name__ == "__main__":
    unittest.main()
