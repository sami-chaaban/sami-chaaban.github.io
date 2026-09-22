"""Native worker admission, cleanup ownership, recycling and failure recovery."""
from __future__ import annotations

import asyncio
from concurrent.futures import Future, ThreadPoolExecutor
from concurrent.futures.process import BrokenProcessPool
import os
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from api import analysis_worker
from api.analysis_worker import AnalysisBusy, AnalysisDisconnected, BoundedAnalysisRunner


def worker_pid():
    return os.getpid()


def crash_worker():
    os._exit(7)


class Client:
    disconnected = False

    async def is_disconnected(self):
        return self.disconnected


class ControlledExecutor:
    def __init__(self):
        self.submitted = threading.Event()
        self.shutdown_started = threading.Event()
        self.allow_shutdown = threading.Event()
        self.future = Future()
        self.shutdown_calls = []

    def submit(self, function, *args):
        self.submitted.set()
        return self.future

    def shutdown(self, wait, cancel_futures):
        self.shutdown_calls.append((wait, cancel_futures, threading.current_thread().ident))
        self.shutdown_started.set()
        if not self.allow_shutdown.wait(3):
            raise AssertionError("Test did not permit executor shutdown")


async def until(predicate):
    async def poll():
        while not predicate():
            await asyncio.sleep(0.005)
    await asyncio.wait_for(poll(), timeout=5)


class WorkerLifecycleTests(unittest.IsolatedAsyncioTestCase):
    async def test_recycled_jobs_use_distinct_pids_and_return_after_process_exit(self):
        runner = BoundedAnalysisRunner(recycle_workers=True)
        callbacks = []
        try:
            first = await runner.run(worker_pid, on_complete=lambda: callbacks.append("first"))
            second = await runner.run(worker_pid, on_complete=lambda: callbacks.append("second"))
            self.assertNotEqual(first, second)
            self.assertNotIn(os.getpid(), (first, second))
            for pid in (first, second):
                with self.assertRaises(ProcessLookupError):
                    os.kill(pid, 0)
            self.assertEqual(callbacks, ["first", "second"])
            self.assertEqual(runner.pending, 0)
        finally:
            runner.shutdown()

    async def test_cancelled_request_holds_gate_slot_and_lease_through_worker_exit(self):
        executor = ControlledExecutor()
        gate = threading.Semaphore(1)
        runner = BoundedAnalysisRunner(recycle_workers=True, execution_gate=gate)
        releases = []
        with patch.object(analysis_worker, "ProcessPoolExecutor", return_value=executor):
            request = asyncio.create_task(runner.run(worker_pid, on_complete=lambda: releases.append(True)))
            await until(executor.submitted.is_set)
            executor.future.set_result(123)
            await until(executor.shutdown_started.is_set)
            request.cancel()
            with self.assertRaises(asyncio.CancelledError):
                await request
            try:
                self.assertFalse(gate.acquire(blocking=False))
                self.assertEqual(runner.pending, 1)
                self.assertEqual(releases, [])
                # Shutdown must be in another thread, leaving the loop responsive.
                self.assertNotEqual(executor.shutdown_calls[0][2], threading.current_thread().ident)
                self.assertEqual(executor.shutdown_calls[0][:2], (True, True))
            finally:
                executor.allow_shutdown.set()
                await until(lambda: runner.pending == 0)
            self.assertTrue(gate.acquire(blocking=False))
            gate.release()
            self.assertEqual(releases, [True])
            self.assertEqual(len(executor.shutdown_calls), 1)

    async def test_disconnected_native_job_finishes_before_lease_and_gate_release(self):
        executor = ControlledExecutor()
        gate = threading.Semaphore(1)
        runner = BoundedAnalysisRunner(recycle_workers=True, execution_gate=gate)
        client, releases = Client(), []
        with patch.object(analysis_worker, "ProcessPoolExecutor", return_value=executor):
            request = asyncio.create_task(runner.run(worker_pid, request=client,
                                                     on_complete=lambda: releases.append(True)))
            await until(executor.submitted.is_set)
            client.disconnected = True
            with self.assertRaises(AnalysisDisconnected):
                await request
            self.assertEqual(releases, [])
            executor.future.set_result(123)
            await until(executor.shutdown_started.is_set)
            try:
                self.assertFalse(gate.acquire(blocking=False))
                self.assertEqual(releases, [])
            finally:
                executor.allow_shutdown.set()
                await until(lambda: runner.pending == 0)
            self.assertEqual(releases, [True])

    async def test_next_job_cannot_start_during_previous_worker_shutdown(self):
        first_executor, second_executor = ControlledExecutor(), ControlledExecutor()
        gate = threading.Semaphore(1)
        runner = BoundedAnalysisRunner(workers=2, recycle_workers=True, execution_gate=gate)
        with patch.object(analysis_worker, "ProcessPoolExecutor", side_effect=[first_executor, second_executor]):
            first = asyncio.create_task(runner.run(worker_pid))
            await until(first_executor.submitted.is_set)
            second = asyncio.create_task(runner.run(worker_pid))
            first_executor.future.set_result(1)
            await until(first_executor.shutdown_started.is_set)
            try:
                await asyncio.sleep(0.06)
                self.assertFalse(second_executor.submitted.is_set())
                self.assertFalse(first.done())
            finally:
                first_executor.allow_shutdown.set()
                self.assertEqual(await first, 1)
                await until(second_executor.submitted.is_set)
                second_executor.future.set_result(2)
                second_executor.allow_shutdown.set()
                self.assertEqual(await second, 2)
            self.assertEqual(runner.pending, 0)

    async def test_direct_lifecycle_cancellation_still_waits_for_worker_exit(self):
        executor = ControlledExecutor()
        gate = threading.Semaphore(1)
        runner = BoundedAnalysisRunner(recycle_workers=True, execution_gate=gate)
        releases = []
        with patch.object(analysis_worker, "ProcessPoolExecutor", return_value=executor):
            request = asyncio.create_task(runner.run(worker_pid, on_complete=lambda: releases.append(True)))
            await until(executor.submitted.is_set)
            executor.future.set_result(123)
            await until(executor.shutdown_started.is_set)
            lifecycle = next(iter(runner._jobs))
            lifecycle.cancel()
            try:
                await asyncio.sleep(0.01)
                self.assertFalse(gate.acquire(blocking=False))
                self.assertEqual(releases, [])
            finally:
                executor.allow_shutdown.set()
                with self.assertRaises(asyncio.CancelledError):
                    await request
            self.assertEqual(releases, [True])
            self.assertEqual(runner.pending, 0)

    async def test_gate_wait_cancellation_releases_lease_without_submitting(self):
        gate = threading.Semaphore(0)
        executor = ThreadPoolExecutor(max_workers=1)
        runner = BoundedAnalysisRunner(executor=executor, execution_gate=gate)
        releases, calls = [], []
        request = asyncio.create_task(runner.run(lambda: calls.append(True),
                                                 on_complete=lambda: releases.append(True)))
        await until(lambda: runner.pending == 1 and runner.slots.locked())
        request.cancel()
        with self.assertRaises(asyncio.CancelledError):
            await request
        self.assertEqual((calls, releases, runner.pending), ([], [True], 0))
        gate.release()
        self.assertTrue(gate.acquire(blocking=False))
        runner.shutdown()

    async def test_disconnected_queue_and_queue_full_release_each_lease_once(self):
        gate = threading.Semaphore(0)
        runner = BoundedAnalysisRunner(max_pending=2, execution_gate=gate)
        releases = []
        first = asyncio.create_task(runner.run(worker_pid, on_complete=lambda: releases.append("first")))
        await until(lambda: runner.slots.locked())
        client = Client()
        second = asyncio.create_task(runner.run(worker_pid, request=client,
                                                on_complete=lambda: releases.append("second")))
        await until(lambda: runner.pending == 2)
        try:
            with self.assertRaises(AnalysisBusy):
                await runner.run(worker_pid, on_complete=lambda: releases.append("rejected"))
            client.disconnected = True
            with self.assertRaises(AnalysisDisconnected):
                await second
        finally:
            first.cancel()
            with self.assertRaises(asyncio.CancelledError):
                await first
        self.assertEqual(sorted(releases), ["first", "rejected", "second"])
        self.assertEqual(runner.pending, 0)

    async def test_persistent_pool_recovers_after_native_process_crash(self):
        gate = threading.Semaphore(1)
        runner = BoundedAnalysisRunner(execution_gate=gate)
        releases = []
        try:
            with self.assertRaises(BrokenProcessPool):
                await runner.run(crash_worker, on_complete=lambda: releases.append("failed"))
            self.assertIsNone(runner.executor)
            self.assertEqual(runner.pending, 0)
            result = await runner.run(worker_pid, on_complete=lambda: releases.append("recovered"))
            self.assertNotEqual(result, os.getpid())
            self.assertEqual(releases, ["failed", "recovered"])
            self.assertTrue(gate.acquire(blocking=False))
            gate.release()
        finally:
            runner.shutdown()

    async def test_injected_executor_is_reused_when_recycling_requested(self):
        executor = ThreadPoolExecutor(max_workers=1)
        runner = BoundedAnalysisRunner(executor=executor, recycle_workers=True)
        try:
            self.assertEqual(await runner.run(lambda: 1), 1)
            self.assertEqual(await runner.run(lambda: 2), 2)
            self.assertIs(runner.executor, executor)
        finally:
            runner.shutdown()


class StoredInputTests(unittest.TestCase):
    def test_stored_input_is_read_inside_worker_entrypoint(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "structure.cif"
            path.write_text("data_test\n#\n", encoding="utf-8")
            with patch.object(analysis_worker, "analyze_and_serialize", return_value="report") as analyze:
                result = analysis_worker.analyze_stored_and_serialize(
                    path, "A", "B", "all", "mmcif", "A:42", "7z8g")
            self.assertEqual(result, "report")
            analyze.assert_called_once_with("data_test\n#\n", "A", "B", "all", "mmcif", "A:42", "7z8g")


if __name__ == "__main__":
    unittest.main()
