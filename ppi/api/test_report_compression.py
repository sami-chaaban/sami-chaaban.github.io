"""Lossless canonical delivery and bounded serialized-report cache regressions."""
import asyncio
import gzip
import json
import pickle
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from api import main
from api.analysis_worker import prepare_delivery
from api.cache import ReportCache


class CompressionTests(unittest.IsolatedAsyncioTestCase):
    def fixture(self):
        path = Path(__file__).resolve().parents[2] / 'roami-tests/interaction-semantics-2026-09-20/canonical-fixture.json'
        return json.loads(path.read_text())

    async def test_actual_semantics_roundtrip_and_ipc_keep_only_compressed_canonical(self):
        original=self.fixture();delivery=prepare_delivery(original);restored=pickle.loads(pickle.dumps(delivery))
        canonical=json.loads(restored.canonical)
        self.assertEqual(canonical['contacts'],original['contacts'])
        self.assertEqual(canonical['perResidue'],original['perResidue'])
        self.assertNotIn('canonical',vars(delivery));self.assertNotIn('canonical',vars(restored))
        self.assertTrue(restored.canonical_gzip.startswith(b'\x1f\x8b'))
        self.assertLess(len(restored.canonical_gzip),len(restored.canonical)/3)
        self.assertEqual(restored.cache_size_bytes,len(restored.canonical_gzip)+len(restored.display))
        delivery.canonical
        self.assertNotIn('canonical',vars(delivery))

    async def test_canonical_http_negotiation_preserves_json_and_skips_recompression(self):
        delivery=prepare_delivery(self.fixture())
        accepted=SimpleNamespace(headers={'accept-encoding':'gzip, deflate, br'})
        compressed=await main._delivery_response(delivery,True,accepted)
        self.assertIs(compressed.body,delivery.canonical_gzip)
        self.assertEqual(compressed.headers['content-encoding'],'gzip')
        self.assertEqual(compressed.headers['vary'],'Accept-Encoding')
        self.assertEqual(gzip.decompress(compressed.body),delivery.canonical)
        plain=await main._delivery_response(delivery,True,SimpleNamespace(headers={'accept-encoding':'gzip;q=0, identity'}))
        self.assertEqual(plain.body,delivery.canonical);self.assertNotIn('content-encoding',plain.headers)
        compact=await main._delivery_response(delivery,False,accepted)
        self.assertIs(compact.body,delivery.display)

    async def test_report_store_retains_gzip_and_uncompressed_reads_run_off_event_loop(self):
        delivery=prepare_delivery(self.fixture());cache=ReportCache(max_bytes=1_000_000)
        with patch.object(main,'report_store',cache),patch.object(main.asyncio,'to_thread',wraps=asyncio.to_thread) as offload:
            cache.set(delivery.report_id,delivery.canonical_gzip)
            response=await main.get_report(delivery.report_id)
            self.assertEqual(response.body,delivery.canonical);offload.assert_awaited_once()
            self.assertIs(cache.get(delivery.report_id),delivery.canonical_gzip)
            response=await main.get_report(delivery.report_id,SimpleNamespace(headers={'accept-encoding':'gzip'}))
            self.assertIs(response.body,delivery.canonical_gzip)
            self.assertEqual(offload.await_count,1)

    async def test_actual_asgi_report_route_serves_identical_json_with_and_without_gzip(self):
        delivery=prepare_delivery(self.fixture())
        with patch.object(main,'report_store',ReportCache(max_bytes=1_000_000)):
            main.report_store.set(delivery.report_id,delivery.canonical_gzip)
            for encoding in ('gzip','identity'):
                sent=[]
                async def receive():return {'type':'http.request','body':b'','more_body':False}
                async def send(message):sent.append(message)
                path='/report/'+delivery.report_id
                scope={'type':'http','asgi':{'version':'3.0'},'http_version':'1.1','method':'GET','scheme':'http','path':path,'raw_path':path.encode(),'query_string':b'','root_path':'','headers':[(b'accept-encoding',encoding.encode())],'server':('test',80),'client':('test',1)}
                await main.app(scope,receive,send)
                start=next(m for m in sent if m['type']=='http.response.start');headers=dict(start['headers']);body=b''.join(m.get('body',b'') for m in sent if m['type']=='http.response.body')
                self.assertEqual(start['status'],200)
                if encoding=='gzip':
                    self.assertEqual(headers[b'content-encoding'],b'gzip');self.assertEqual(body,delivery.canonical_gzip)
                    body=gzip.decompress(body)
                else:self.assertNotIn(b'content-encoding',headers)
                self.assertEqual(body,delivery.canonical)

    def test_explicit_gzip_rejection_overrides_wildcard(self):
        self.assertFalse(main._accepts_gzip(SimpleNamespace(headers={'accept-encoding':'gzip;q=0, *;q=1'})))
        self.assertTrue(main._accepts_gzip(SimpleNamespace(headers={'accept-encoding':'br, *;q=.5'})))
        self.assertFalse(main._accepts_gzip(None))


class CacheBudgetTests(unittest.TestCase):
    def test_evicts_until_byte_budget_fits_and_accounts_replacement(self):
        cache=ReportCache(max_entries=10,max_bytes=10)
        with patch('api.cache.time.time',side_effect=[1,2,3,4,5,6,7,8,9,10]):
            cache.set('a',b'aaaa');cache.set('b',b'bbbb');cache.set('c',b'ccccccc')
            self.assertIsNone(cache.get('a'));self.assertIsNone(cache.get('b'));self.assertEqual(cache.total_bytes,7)
            cache.set('c',b'xx');self.assertEqual(cache.total_bytes,2)
            cache.set('d',b'12345678');self.assertEqual(cache.total_bytes,10)
        self.assertEqual(set(cache._entries),{'c','d'})

    def test_oversized_entry_is_not_retained_or_allowed_to_leave_stale_replacement(self):
        cache=ReportCache(max_bytes=4);cache.set('a',b'1234');cache.set('a',b'12345')
        self.assertIsNone(cache.get('a'));self.assertEqual(cache.total_bytes,0)
        cache.set('b',b'abcde');self.assertEqual(cache.total_bytes,0)

    def test_ttl_removal_and_entry_limit_release_payload_budget(self):
        cache=ReportCache(ttl_seconds=10,max_entries=1,max_bytes=20)
        with patch('api.cache.time.time',return_value=0):cache.set('a',b'12345')
        with patch('api.cache.time.time',return_value=11):self.assertIsNone(cache.get('a'))
        self.assertEqual(cache.total_bytes,0)
        cache.set('b',b'123');cache.set('c',b'123456')
        self.assertEqual(cache.total_bytes,6);self.assertIsNone(cache.get('b'))

    def test_delivery_budget_counts_compressed_canonical_plus_display(self):
        delivery=prepare_delivery({'contacts':{'other':[{'debugOnly':True,'evidence':'x'*1000}]}})
        cache=ReportCache(max_bytes=delivery.cache_size_bytes)
        cache.set('a',delivery);self.assertIs(cache.get('a'),delivery)
        self.assertEqual(cache.total_bytes,len(delivery.canonical_gzip)+len(delivery.display))
        cache.set('b',delivery);self.assertIsNone(cache.get('a'));self.assertIs(cache.get('b'),delivery)

    def test_custom_sizer_and_zero_budget(self):
        cache=ReportCache(max_bytes=3,size_of=lambda item:len(item['payload']))
        cache.set('x',{'payload':b'1234'});self.assertIsNone(cache.get('x'))
        cache.set('x',{'payload':b'123'});self.assertEqual(cache.total_bytes,3)
        disabled=ReportCache(max_bytes=0);disabled.set('x',b'1');self.assertIsNone(disabled.get('x'))


if __name__=='__main__':unittest.main()
