"""Bounded scientific execution and lossless report delivery projections."""
from __future__ import annotations

import asyncio
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import json
import gzip
import multiprocessing
from typing import Callable
import uuid


class AnalysisBusy(RuntimeError):
    pass


class AnalysisDisconnected(RuntimeError):
    pass


class BoundedAnalysisRunner:
    """Admit jobs before submitting: the process pool never has a stale backlog.

    A disconnected/cancelled running job finishes in its existing worker, retaining
    its slot until completion. Killing native chemistry midway is unsafe; waiting
    jobs, however, can be discarded without ever invoking the scientific engine.
    """

    def __init__(self, workers: int = 1, max_pending: int = 8, executor=None):
        self.workers = max(1, workers)
        self.max_pending = max(self.workers, max_pending)
        self.executor = executor
        self.slots = asyncio.Semaphore(self.workers)
        self.pending = 0

    async def run(self, function: Callable, *args, request=None):
        if self.pending >= self.max_pending:
            raise AnalysisBusy("Analysis queue is full. Please retry shortly.")
        self.pending += 1
        acquired = False
        submitted = False
        admission = asyncio.create_task(self.slots.acquire())
        try:
            while not admission.done():
                await asyncio.wait({admission}, timeout=0.05)
                if request is not None and await request.is_disconnected():
                    raise AnalysisDisconnected("Analysis request disconnected")
            acquired = await admission
            if request is not None and await request.is_disconnected():
                raise AnalysisDisconnected("Analysis request disconnected")
            if self.executor is None:
                self.executor = ProcessPoolExecutor(
                    max_workers=self.workers,
                    mp_context=multiprocessing.get_context("spawn"),
                )
            loop = asyncio.get_running_loop()
            future = self.executor.submit(function, *args)
            submitted = True

            def release_slot(_):
                def release():
                    self.slots.release()
                    self.pending -= 1
                if not loop.is_closed():
                    loop.call_soon_threadsafe(release)

            future.add_done_callback(release_slot)
            result = asyncio.wrap_future(future)
            # Retrieve abandoned exceptions as well, without cancelling native jobs.
            result.add_done_callback(lambda done: None if done.cancelled() else done.exception())
            while not result.done():
                await asyncio.wait({result}, timeout=0.05)
                if request is not None and await request.is_disconnected():
                    raise AnalysisDisconnected("Analysis request disconnected")
            return result.result()
        finally:
            if not submitted:
                # acquire can complete while cancellation/disconnect is being checked.
                if not admission.done():
                    admission.cancel()
                    await asyncio.gather(admission, return_exceptions=True)
                elif not admission.cancelled() and admission.exception() is None:
                    acquired = bool(admission.result())
                if acquired:
                    self.slots.release()
                self.pending -= 1

    def shutdown(self):
        if self.executor is not None:
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
