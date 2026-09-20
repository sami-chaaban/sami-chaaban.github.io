"""API delivery/scheduling regressions from the September 2026 audit."""
from __future__ import annotations

import asyncio
from concurrent.futures import ThreadPoolExecutor
import importlib.util
import json
import os
from pathlib import Path
import sys
import threading
import time
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from api import analysis, main
from api.analysis_worker import AnalysisBusy, AnalysisDisconnected, BoundedAnalysisRunner, prepare_delivery
from api.cache import ReportCache
from api.explain import contact_counts, explain_report, top_residues
from api.models import AnalyzeRequest


def worker_pid():
    return os.getpid()


class Client:
    disconnected = False

    async def is_disconnected(self):
        return self.disconnected


class SchedulingTests(unittest.IsolatedAsyncioTestCase):
    async def test_production_runner_uses_another_process(self):
        runner = BoundedAnalysisRunner()
        try:
            self.assertNotEqual(await runner.run(worker_pid), os.getpid())
        finally:
            runner.shutdown()

    async def test_http_health_stays_responsive_and_native_concurrency_is_bounded(self):
        runner = BoundedAnalysisRunner(workers=1, executor=ThreadPoolExecutor(max_workers=4))
        active = peak = 0
        started = threading.Event()
        release = threading.Event()

        def compute(*args):
            nonlocal active, peak
            active += 1
            peak = max(peak, active)
            started.set()
            release.wait(2)
            active -= 1
            return prepare_delivery({"contacts": {}, "input": args[0]})

        with patch.object(main, "analysis_runner", runner), patch.object(main, "analyze_and_serialize", compute), patch.object(main, "cache", ReportCache()), patch.object(main, "report_store", ReportCache()):
            jobs = [asyncio.create_task(main.analyze(AnalyzeRequest(pdbText=f"INPUT{i}", chainA="A", chainB="B"), Client())) for i in range(2)]
            try:
                while not started.is_set():
                    await asyncio.sleep(0.005)
                self.assertEqual(await asyncio.wait_for(main.health(), timeout=0.1), {"status": "ok"})
                await asyncio.sleep(0.01)
                self.assertFalse(jobs[0].done())
            finally:
                release.set()
                await asyncio.gather(*jobs)
                runner.shutdown()
        self.assertEqual(peak, 1)

    async def test_disconnected_waiting_job_never_reaches_native_executor(self):
        runner = BoundedAnalysisRunner(workers=1, max_pending=2, executor=ThreadPoolExecutor(max_workers=2))
        started, release = threading.Event(), threading.Event()
        calls = []
        def native(value):
            calls.append(value)
            started.set()
            release.wait(2)
            return value
        first = asyncio.create_task(runner.run(native, "first"))
        while not started.is_set():
            await asyncio.sleep(0.005)
        client = Client()
        queued = asyncio.create_task(runner.run(native, "stale", request=client))
        await asyncio.sleep(0)
        with self.assertRaises(AnalysisBusy):
            await runner.run(native, "excess")
        client.disconnected = True
        try:
            with self.assertRaises(AnalysisDisconnected):
                await queued
        finally:
            release.set()
            await first
            runner.shutdown()
        self.assertEqual(calls, ["first"])
        self.assertEqual(runner.pending, 0)

    async def test_cancelling_running_request_keeps_native_slot_until_completion(self):
        runner = BoundedAnalysisRunner(workers=1, executor=ThreadPoolExecutor(max_workers=2))
        started, release = threading.Event(), threading.Event()
        calls = []
        def native(value):
            calls.append(value)
            started.set()
            release.wait(2)
            return value
        first = asyncio.create_task(runner.run(native, "first"))
        while not started.is_set():
            await asyncio.sleep(0.005)
        first.cancel()
        with self.assertRaises(asyncio.CancelledError):
            await first
        second = asyncio.create_task(runner.run(native, "second"))
        await asyncio.sleep(0.02)
        self.assertEqual(calls, ["first"])
        release.set()
        await second
        runner.shutdown()
        self.assertEqual(calls, ["first", "second"])


class ReportDeliveryTests(unittest.IsolatedAsyncioTestCase):
    async def test_projection_keeps_weak_visible_records_and_full_report_is_retrievable(self):
        original = {
            "contacts": {"other": [
                {"id": "visible", "confidence": "low"}, {"id": "a", "debugOnly": True},
                {"id": "b", "debug_only": True}, {"id": "c", "asserted": {"debugOnly": True}},
            ]}, "perResidue": {"A:1": {"total": 1}},
        }
        delivery = prepare_delivery(original)
        full, display = json.loads(delivery.canonical), json.loads(delivery.display)
        self.assertEqual(full["contacts"], original["contacts"])
        self.assertNotIn("reportId", original)
        self.assertEqual(display["contacts"]["other"], [{"id": "visible", "confidence": "low"}])
        self.assertEqual(display["meta"]["diagnostics"]["omittedByBucket"], {"other": 3})
        self.assertEqual(display["perResidue"], original["perResidue"])
        with patch.object(main, "report_store", ReportCache()):
            main.report_store.set(delivery.report_id, delivery.canonical_gzip)
            response = await main.get_report(delivery.report_id)
            self.assertEqual(json.loads(response.body), full)

    async def test_cache_uses_coordinates_even_with_pdb_id_and_projection_is_not_a_new_analysis(self):
        runner = BoundedAnalysisRunner(executor=ThreadPoolExecutor(max_workers=1))
        calls = []
        def compute(*args):
            calls.append(args[0])
            return prepare_delivery({"contacts": {"other": [{"debugOnly": True}]}, "input": args[0]})
        with patch.object(main, "analysis_runner", runner), patch.object(main, "analyze_and_serialize", compute), patch.object(main, "cache", ReportCache()), patch.object(main, "report_store", ReportCache()):
            first = await main.analyze(AnalyzeRequest(pdbId="1abc", pdbText="A", chainA="A", chainB="B"), Client())
            second = await main.analyze(AnalyzeRequest(pdbId="1abc", pdbText="B", chainA="A", chainB="B"), Client())
            compact = await main.analyze(AnalyzeRequest(pdbId="1abc", pdbText="B", chainA="A", chainB="B", includeDiagnostics=False), Client())
            self.assertTrue(all(entry.value.startswith(b"\x1f\x8b") for entry in main.report_store._entries.values()))
            self.assertTrue(all("canonical" not in vars(entry.value) for entry in main.cache._entries.values()))
            runner.shutdown()
        self.assertEqual(calls, ["A", "B"])
        self.assertNotEqual(json.loads(first.body)["reportId"], json.loads(second.body)["reportId"])
        self.assertEqual(json.loads(compact.body)["reportId"], json.loads(second.body)["reportId"])
        self.assertEqual(json.loads(compact.body)["contacts"]["other"], [])


class CifFilterTests(unittest.TestCase):
    def test_packed_wrapped_quoted_rows_preserve_insertion_codes_and_other_categories(self):
        text = """data_test
_entry.id test
loop_
_atom_site.auth_asym_id
_atom_site.label_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.label_atom_id
'A chain' A 42 A "O5'" B B 10 ? N
'A chain' A
42 B CA
#
_struct.title 'kept title'
"""
        filtered, count = main.filter_mmcif_text_to_single_chain(text, "A chain")
        block = analysis._gemmi.cif.read_string(filtered).sole_block()
        self.assertEqual(count, 2)
        self.assertEqual(list(block.find_values("_atom_site.pdbx_PDB_ins_code")), ["A", "B"])
        self.assertEqual(analysis._gemmi.cif.as_string(block.find_value("_struct.title")), "kept title")
        self.assertEqual(analysis._gemmi.cif.as_string(block.find_values("_atom_site.label_atom_id")[0]), "O5'")

    def test_legacy_ribbon_import_uses_package_layout(self):
        with patch("api.ribbon_backend.structure_to_ribbon_json", return_value={"chains": []}) as generate:
            self.assertEqual(main.ribbon(main.RibbonRequest(pdbText="END\n")), {"chains": []})
            self.assertTrue(generate.called)


class ExplanationTests(unittest.TestCase):
    def test_organism_metadata_separates_biological_sources_from_expression_hosts(self):
        payload = {
            "rcsb_entity_source_organism": [{"ncbi_scientific_name": "Homo sapiens"}],
            "rcsb_entity_host_organism": [{"ncbi_scientific_name": "Escherichia coli"}],
            "hosts": [{"scientific_name": "Spodoptera frugiperda"}],
            "polymer_entities": [{"entity_src_nat": [{"pdbx_organism_scientific": "Sus scrofa"}]}],
        }
        self.assertEqual(main.normalize_organism_list(payload), ["Homo sapiens", "Sus scrofa"])
        self.assertEqual(main.normalize_organism_list({"host": payload["rcsb_entity_host_organism"]}), [])

    def test_all_categories_and_uncertainty_are_explained_without_image_invention(self):
        contact = {"residueA": {"chain": "A", "seq": "1"}, "residueB": {"chain": "B", "seq": "2"}}
        narrative = explain_report({"contacts": {"polar_contacts": [contact], "halogen_bonds": [contact], "base_pairing": [contact]}, "perResidue": {"A:1": {"total": 100}}}, images=["arbitrary string"])
        for phrase in ("1 polar contact", "1 halogen-bond assignment", "1 base-pairing contact", "most connected residues", "A:1 (3 assignments)", "not calibrated probabilities", "has not inspected"):
            self.assertIn(phrase, narrative)
        self.assertNotIn("100 assignments", narrative)
        self.assertNotIn("strongest contact clusters", narrative)

    def test_asserted_clashes_and_diagnostics_are_explained_consistently(self):
        contact = {"asserted": {"family": "clash"}, "residueA": {"chain": "A", "seq": 0}, "residueB": {"chain": "B", "seq": "2A"}}
        diagnostic = {"debugOnly": True, "residueA": {"chain": "Z", "seq": 99}}
        contacts = {"other": [contact, diagnostic], "hydrogen_bonds": [{"asserted": {"debugOnly": True}}]}
        self.assertEqual(contact_counts(contacts), {"clash": 1})
        ranking = top_residues({"Z:99": {"total": 100}}, contacts=contacts)
        self.assertEqual([(r["id"], r["total"]) for r in ranking], [("A:0", 1), ("B:2A", 1)])
        narrative = explain_report({"chainA": "A", "chainB": "A", "contacts": contacts})
        self.assertIn("within chain A", narrative)
        self.assertIn("1 steric clash", narrative)
        self.assertNotIn("other contact", narrative)
        self.assertNotIn("Z:99", narrative)

    def test_same_residue_endpoint_counts_once_and_empty_results_do_not_invent_hotspots(self):
        residue = {"chain": "A", "seq": 1}
        contacts = {"other": [{"residueA": residue, "residueB": residue}]}
        self.assertEqual(top_residues({}, contacts=contacts)[0]["total"], 1)
        narrative = explain_report({"contacts": {"other": [{"debug_only": True}]}, "perResidue": {"A:1": {"total": 10}}})
        self.assertIn("no non-diagnostic contact assignments under the current criteria", narrative)
        self.assertNotIn("most connected residues", narrative)


class DependencyPatchTests(unittest.TestCase):
    def test_patch_is_idempotent_and_rejects_unreviewed_source(self):
        path = Path(__file__).resolve().parents[1] / "scripts" / "patch_arpeggio.py"
        spec = importlib.util.spec_from_file_location("patch_arpeggio", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        import arpeggio
        source = (Path(arpeggio.__file__).parent / "core" / "interactions.py").read_bytes()
        patched = module.patched_source(source)
        self.assertEqual(module.patched_source(patched), patched)
        with self.assertRaises(RuntimeError):
            module.patched_source(source + b"\n# unreviewed change\n")


if __name__ == "__main__":
    unittest.main()
