"""Prevent availability failures from masquerading as safe certificate rejection."""

import socket
import ssl
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import MagicMock, patch
from urllib import error

import check_tls_trust as checker


class TrustProbeTests(unittest.TestCase):
    def test_only_certificate_verification_errors_count_as_rejection(self):
        rejection = ssl.SSLCertVerificationError(1, "certificate verify failed")
        for exc in (rejection, error.URLError(rejection)):
            with self.subTest(exc=exc):
                self.assertEqual(checker.classify_error(exc), "certificate_rejected")
        for exc in (
            error.URLError("certificate verify failed"),
            error.URLError(socket.gaierror("DNS failed")),
            TimeoutError("timed out"),
            ConnectionResetError("reset"),
            ssl.SSLError("TLS handshake failed"),
        ):
            with self.subTest(exc=exc):
                self.assertEqual(checker.classify_error(exc), "request_failed")

    def test_http_errors_do_not_count_as_rejection(self):
        for status in (301, 403, 500):
            exc = error.HTTPError("https://example.com/", status, "response", {}, None)
            with self.subTest(status=status), patch.object(checker.request, "build_opener") as factory:
                factory.return_value.open.side_effect = exc
                self.assertEqual(checker.probe("https://example.com/", ssl.create_default_context())[0], "request_failed")

    def test_successful_https_response_is_verified(self):
        response = MagicMock()
        response.__enter__.return_value.status = 200
        with patch.object(checker.request, "build_opener") as factory:
            factory.return_value.open.return_value = response
            self.assertEqual(checker.probe("https://example.com/", ssl.create_default_context())[0], "verified")

    def test_version_floor(self):
        for version in ("2026.10.19", "2026.10.23", "2027.1.1"):
            checker.require_patched("test", version)
        for version in ("2026.2.25", "2026.10.18", "unknown"):
            with self.subTest(version=version), self.assertRaises(RuntimeError):
                checker.require_patched("test", version)

    def test_stale_pip_bundle_fails_even_with_patched_environment_packages(self):
        certifi = SimpleNamespace(__version__="2026.10.19", where=lambda: "/test/ca.pem")
        bundled = SimpleNamespace(__version__="2026.6.17", where=lambda: "/test/pip-ca.pem")
        modules = {"certifi": certifi, "pip._vendor": SimpleNamespace(certifi=bundled)}
        with tempfile.TemporaryDirectory() as prefix:
            metadata = Path(prefix) / "conda-meta"
            metadata.mkdir()
            (metadata / "ca-certificates-2026.10.19-0.json").write_text(json.dumps({
                "name": "ca-certificates", "version": "2026.10.19",
            }))
            with patch.dict(checker.sys.modules, modules), patch.object(checker.sys, "prefix", prefix), patch("builtins.print"):
                with self.assertRaisesRegex(RuntimeError, "pip's bundled certifi"):
                    checker.check_packages()
                bundled.__version__ = "2026.10.19"
                self.assertEqual(checker.check_packages(), "/test/ca.pem")


if __name__ == "__main__":
    unittest.main()
