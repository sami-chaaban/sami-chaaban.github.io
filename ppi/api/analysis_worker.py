"""Bounded scientific execution and lossless report delivery projections."""
from __future__ import annotations

import asyncio
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from dataclasses import dataclass
from functools import partial
import json
import gzip
import multiprocessing
from pathlib import Path
from typing import Callable
import uuid


class AnalysisBusy(RuntimeError):
    pass


class AnalysisDisconnected(RuntimeError):
    pass


class BoundedAnalysisRunner:
    """Admit jobs before submitting: the process pool never has a stale backlog.

    A disconnected/cancelled running job finishes in its existing worker. Recycled
    jobs retain their slot and shared gate until that process exits; waiting jobs
    can be discarded without ever invoking the scientific engine. on_complete
    releases caller-owned input resources after this lifecycle, including errors.
    """

    def __init__(self, workers: int = 1, max_pending: int = 8, executor=None,
                 recycle_workers: bool = False, execution_gate=None):
        self.workers = max(1, workers)
        self.max_pending = max(self.workers, max_pending)
        self.executor = executor
        self.recycle_workers = bool(recycle_workers)
        self.execution_gate = execution_gate
        # An injected executor retains its existing lifetime/ownership contract.
        self._injected_executor = executor is not None
        self._jobs = set()
        self._retirements = {}
        self._shutdown_requested = False
        self.slots = asyncio.Semaphore(self.workers)
        self.pending = 0

    async def run(self, function: Callable, *args, request=None, on_complete=None):
        if self.pending >= self.max_pending:
            try:
                raise AnalysisBusy("Analysis queue is full. Please retry shortly.")
            finally:
                if on_complete is not None:
                    on_complete()
        self.pending += 1
        acquired = False
        delegated = False
        gate_acquired = False
        admission = asyncio.create_task(self.slots.acquire())
        try:
            while not admission.done():
                await asyncio.wait({admission}, timeout=0.05)
                if request is not None and await request.is_disconnected():
                    raise AnalysisDisconnected("Analysis request disconnected")
            acquired = await admission
            if request is not None and await request.is_disconnected():
                raise AnalysisDisconnected("Analysis request disconnected")
            if self.execution_gate is not None:
                # Nonblocking acquisition avoids tying up the loop or a thread
                # whose eventual acquisition could outlive a cancelled waiter.
                while not self.execution_gate.acquire(blocking=False):
                    await asyncio.sleep(0.05)
                    if request is not None and await request.is_disconnected():
                        raise AnalysisDisconnected("Analysis request disconnected")
                gate_acquired = True
            if request is not None and await request.is_disconnected():
                raise AnalysisDisconnected("Analysis request disconnected")
            # This independent task owns the slot and gate from here through
            # native completion and (when recycling) full worker process exit.
            # Cancelling the HTTP request must never cancel that cleanup.
            result = asyncio.create_task(self._execute(function, args, gate_acquired, on_complete))
            delegated = True
            self._jobs.add(result)

            def finished(done):
                self._jobs.discard(done)
                if not done.cancelled():
                    done.exception()  # Consume exceptions from disconnected jobs.

            result.add_done_callback(finished)
            while not result.done():
                await asyncio.wait({result}, timeout=0.05)
                if request is not None and await request.is_disconnected():
                    raise AnalysisDisconnected("Analysis request disconnected")
            return result.result()
        finally:
            if not delegated:
                # acquire can complete while cancellation/disconnect is being checked.
                if not admission.done():
                    admission.cancel()
                    await asyncio.gather(admission, return_exceptions=True)
                elif not admission.cancelled() and admission.exception() is None:
                    acquired = bool(admission.result())
                if acquired:
                    self.slots.release()
                if gate_acquired:
                    self.execution_gate.release()
                self.pending -= 1
                if on_complete is not None:
                    on_complete()

    async def _retire_executor(self, executor):
        cleanup = self._retirements.get(executor)
        if cleanup is None:
            # Use an executor Future rather than a Task: event-loop shutdown
            # cancels Tasks, but must not cancel our worker-exit notification.
            cleanup = asyncio.get_running_loop().run_in_executor(
                None, partial(executor.shutdown, wait=True, cancel_futures=True),
            )
            self._retirements[executor] = cleanup
        cancelled = False
        try:
            # Even event-loop shutdown must not release the shared gate while a
            # native process is still exiting. Request cancellation never reaches
            # this task, but direct lifecycle cancellation can occur at shutdown.
            while True:
                try:
                    await asyncio.shield(cleanup)
                    break
                except asyncio.CancelledError:
                    cancelled = True
                    if cleanup.cancelled():
                        raise
            if cancelled:
                raise asyncio.CancelledError
        finally:
            if cleanup.done():
                self._retirements.pop(executor, None)

    async def _execute(self, function, args, gate_acquired, on_complete):
        executor = None
        recycle = self.recycle_workers and not self._injected_executor
        broken = False
        try:
            if recycle:
                # max_tasks_per_child is unavailable on Python 3.10. A dedicated
                # one-job pool releases native allocator state by process exit.
                executor = ProcessPoolExecutor(
                    max_workers=1, mp_context=multiprocessing.get_context("spawn"),
                )
            else:
                if self.executor is None:
                    self.executor = ProcessPoolExecutor(
                        max_workers=self.workers,
                        mp_context=multiprocessing.get_context("spawn"),
                    )
                executor = self.executor
            try:
                future = executor.submit(function, *args)
                return await asyncio.wrap_future(future)
            except BrokenProcessPool:
                broken = True
                if self.executor is executor:
                    self.executor = None
                raise
        finally:
            try:
                if executor is not None and (recycle or broken or self._shutdown_requested):
                    if self.executor is executor:
                        self.executor = None
                    await self._retire_executor(executor)
            finally:
                if gate_acquired:
                    self.execution_gate.release()
                self.slots.release()
                self.pending -= 1
                if on_complete is not None:
                    on_complete()

    def shutdown(self):
        self._shutdown_requested = True
        # Submitted jobs retain ownership until their asynchronous cleanup ends.
        # An idle persistent executor can start shutdown without blocking the loop.
        if self.executor is not None and not any(not job.done() for job in self._jobs):
            self.executor.shutdown(wait=False, cancel_futures=True)
            self.executor = None


@dataclass(frozen=True)
class ReportDelivery:
    report_id: str
    canonical_gzip: bytes
    display: bytes

    @property
    def canonical(self) -> bytes:
        """Compatibility view; never retain an uncompressed copy in the cache."""
        return gzip.decompress(self.canonical_gzip)

    @property
    def cache_size_bytes(self) -> int:
        return len(self.canonical_gzip) + len(self.display)


def is_diagnostic_contact(contact: dict) -> bool:
    asserted = contact.get("asserted") or {}
    return bool(contact.get("debugOnly") or contact.get("debug_only") or asserted.get("debugOnly"))


def prepare_delivery(report: dict) -> ReportDelivery:
    """Keep canonical records untouched; omit only explicitly diagnostic records."""
    report = dict(report)
    report_id = report.get("reportId") or uuid.uuid4().hex[:10]
    report["reportId"] = report_id
    contacts = report.get("contacts") or {}
    display_contacts = {
        bucket: [record for record in records if not is_diagnostic_contact(record)]
        for bucket, records in contacts.items()
    }
    omitted = {bucket: len(records) - len(display_contacts[bucket]) for bucket, records in contacts.items()}
    display = {
        **report,
        "contacts": display_contacts,
        "meta": {
            **(report.get("meta") or {}),
            "diagnostics": {
                "included": False,
                "omittedCount": sum(omitted.values()),
                "omittedByBucket": omitted,
                "reason": "Explicitly diagnostic contacts are available in the canonical report.",
                "canonicalReportUrl": f"/report/{report_id}",
            },
        },
    }
    if isinstance(report.get('displayGroups'), dict):
        from .display_grouping import project_display_groups
        display['displayGroups'] = project_display_groups(report['displayGroups'], display_contacts)
    encode = lambda value: json.dumps(value, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode("utf-8")
    # Compress in the worker before IPC. Most canonical rows are diagnostics,
    # whose repeated chemistry metadata compresses well without pruning evidence.
    canonical_gzip = gzip.compress(encode(report), compresslevel=1, mtime=0)
    return ReportDelivery(report_id, canonical_gzip, encode(display))


def analyze_and_serialize(structure_text, chain_a, chain_b, mode, structure_format, focus_residue, pdb_id):
    # Import in the spawned process: native libraries and their mutable state stay
    # isolated from the HTTP event loop and from concurrently running workers.
    from .analysis import analyze_interface
    report = analyze_interface(
        structure_text, chain_a, chain_b, mode,
        structure_format=structure_format, focus_residue=focus_residue,
    )
    if pdb_id:
        report["pdbId"] = pdb_id
    return prepare_delivery(report)


def analyze_stored_and_serialize(structure_path, chain_a, chain_b, mode, structure_format, focus_residue, pdb_id):
    """Load a server-owned structure only after admission to the native worker."""
    structure_text = Path(structure_path).read_text(encoding="utf-8")
    return analyze_and_serialize(
        structure_text, chain_a, chain_b, mode, structure_format, focus_residue, pdb_id,
    )
