"""Binary HTTP negotiation and persistent worker framing (no native libraries)."""
import asyncio
import gzip
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from api import main
from api.test_mesh_delivery import request, PAYLOAD


class BinaryDeliveryTests(unittest.IsolatedAsyncioTestCase):
    async def test_binary_http_is_compressed_without_numeric_json_conversion(self):
        body = b'ROAMIM01' + struct.pack('<I', 2) + b'{}\0\0' + bytes(8000)
        with patch.object(main, 'run_chapi_mesh_cached', return_value=body), \
             patch.object(main, 'CHAPI_LOW_MEMORY_MODE', False):
            status, headers, received = await request('POST', '/chapi-mesh',
                dict(PAYLOAD, outputFormat='binary'), 'gzip')
        self.assertEqual(status, 200)
        self.assertEqual(headers[b'content-type'], b'application/vnd.roami.mesh')
        self.assertEqual(gzip.decompress(received), body)

    async def test_binary_and_json_have_separate_cache_keys(self):
        payload = {'text': 'data_test', 'format': 'mmcif'}
        self.assertNotEqual(main.build_chapi_mesh_cache_key(payload),
                            main.build_chapi_mesh_cache_key(dict(payload, outputFormat='binary')))


class PersistentBinaryFramingTests(unittest.TestCase):
    def test_closed_pipe_is_detected_even_before_worker_exit(self):
        read_fd, write_fd = os.pipe()
        os.close(write_fd)
        with os.fdopen(read_fd, 'rb') as stdout:
            proc = SimpleNamespace(stdout=stdout, poll=lambda: None,
                                   wait=Mock(side_effect=subprocess.TimeoutExpired('worker', 1)))
            started = time.monotonic()
            with self.assertRaisesRegex(main.ChapiWorkerTransportError, 'closed binary output'):
                main._read_chapi_worker_bytes(proc, 100, 1)
            self.assertLess(time.monotonic() - started, 0.5)

    def test_prefetched_binary_payload_does_not_consume_next_response(self):
        # One write deliberately puts the header and a large body in the same
        # pipe, exercising bytes buffered by the existing line reader.
        with tempfile.TemporaryDirectory() as directory:
            bridge = Path(directory) / 'bridge.py'
            bridge.write_text('''import sys
for request in sys.stdin:
    body = b"ROAMIM01" + bytes(100001)
    sys.stdout.buffer.write(b"BINARY\\t" + str(len(body)).encode() + b"\\n" + body)
    sys.stdout.buffer.flush()
''')
            with patch.object(main, 'CHAPI_BRIDGE', bridge), \
                 patch.object(main, 'CHAPI_PYTHON', sys.executable):
                try:
                    first = main._run_chapi_mesh_worker({'outputFormat': 'binary'})
                    second = main._run_chapi_mesh_worker({'outputFormat': 'binary'})
                finally:
                    with main.CHAPI_WORKER_LOCK:
                        main._stop_chapi_worker_locked()
        self.assertEqual(first, b'ROAMIM01' + bytes(100001))
        self.assertEqual(second, first)

    def test_truncated_binary_frame_is_a_transport_error_and_worker_is_stopped(self):
        with tempfile.TemporaryDirectory() as directory:
            bridge = Path(directory) / 'bridge.py'
            bridge.write_text('''import sys
sys.stdin.readline()
sys.stdout.buffer.write(b"BINARY\\t100\\nROAMIM01")
sys.stdout.buffer.flush()
''')
            with patch.object(main, 'CHAPI_BRIDGE', bridge), \
                 patch.object(main, 'CHAPI_PYTHON', sys.executable):
                with self.assertRaises(main.ChapiWorkerTransportError):
                    main._run_chapi_mesh_worker({'outputFormat': 'binary'})
                self.assertIsNone(main.CHAPI_WORKER_PROC)


if __name__ == '__main__':
    unittest.main()
