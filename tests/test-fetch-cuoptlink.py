"""Unit tests for fetch-cuoptlink.py.

Run with:
    python3 tests/test-fetch-cuoptlink.py -v

Both this test file and the script under test have hyphens in their
filenames, so this file is run directly rather than via `python3 -m
unittest <module>` (hyphenated names aren't importable as modules), and
the script under test is loaded via importlib instead of a normal
`import` statement.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import os
import sys
import tempfile
import unittest
from unittest import mock

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_MODULE_PATH = os.path.join(_REPO_ROOT, "fetch-cuoptlink.py")


def _load_module():
    spec = importlib.util.spec_from_file_location("fetch_cuoptlink", _MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


fc = _load_module()


CUOPT_RELEASE_YAML = """solverConfig:
  - cuopt:
      minVersion: 49
      fileType: 1001
      dictType: 0
      licCodes: 000102030405
      scriptName: gmscuopt.run
      executableName: gmscuopt.out
      modelTypes:
        - LP
        - RMIP
        - MIP
        - QCP
        - MIQCP
        - RMIQCP
"""


class TempGamsDirMixin:
    def setUp(self) -> None:
        self._tmpdir = tempfile.TemporaryDirectory()
        self.gams_dir = self._tmpdir.name
        self.addCleanup(self._tmpdir.cleanup)

    def write(self, name: str, content: str) -> str:
        path = os.path.join(self.gams_dir, name)
        with open(path, "w", encoding="utf-8") as file:
            file.write(content)
        return path

    def read(self, name: str) -> str:
        with open(os.path.join(self.gams_dir, name), encoding="utf-8") as file:
            return file.read()


class QuotedKeyNoDuplicateMergeTest(TempGamsDirMixin, unittest.TestCase):
    """Review finding #1: quoted `cuopt` keys must be recognized as already present."""

    def test_quoted_cuopt_key_prevents_duplicate_merge(self):
        original = (
            "solverConfig:\n"
            '  - "cuopt":\n'
            "      minVersion: 1\n"
        )
        self.write(fc.CONFIG_FILE, original)
        self.write(fc.CUOPT_CONFIG_FILE, CUOPT_RELEASE_YAML)

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            fc._merge_cuopt_config(self.gams_dir)

        self.assertIn("already exists", out.getvalue())
        # The existing config must be left untouched - no duplicate/second entry added.
        self.assertEqual(self.read(fc.CONFIG_FILE), original)
        # The release fragment is still consumed (matches existing behavior for a skipped merge).
        self.assertFalse(os.path.exists(os.path.join(self.gams_dir, fc.CUOPT_CONFIG_FILE)))

    def test_single_quoted_cuopt_key_also_detected(self):
        original = "solverConfig:\n  - 'cuopt':\n      minVersion: 1\n"
        self.write(fc.CONFIG_FILE, original)
        self.write(fc.CUOPT_CONFIG_FILE, CUOPT_RELEASE_YAML)

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            fc._merge_cuopt_config(self.gams_dir)

        self.assertIn("already exists", out.getvalue())
        self.assertEqual(self.read(fc.CONFIG_FILE), original)

    def test_dict_style_quoted_cuopt_key_also_detected(self):
        original = 'solverConfig:\n  "cuopt":\n    minVersion: 1\n'
        self.write(fc.CONFIG_FILE, original)
        self.write(fc.CUOPT_CONFIG_FILE, CUOPT_RELEASE_YAML)

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            fc._merge_cuopt_config(self.gams_dir)

        self.assertIn("already exists", out.getvalue())
        self.assertEqual(self.read(fc.CONFIG_FILE), original)


class MalformedYamlRejectionTest(TempGamsDirMixin, unittest.TestCase):
    """Review finding #2: malformed/unsupported YAML must never be merged into or corrupt gamsconfig.yaml."""

    def assert_rejected(self, original: str, force_fallback: bool = False) -> None:
        self.write(fc.CONFIG_FILE, original)
        self.write(fc.CUOPT_CONFIG_FILE, CUOPT_RELEASE_YAML)

        patch_ctx = mock.patch.object(fc, "_pyyaml", None) if force_fallback else contextlib.nullcontext()
        out = io.StringIO()
        with patch_ctx, contextlib.redirect_stdout(out):
            with self.assertRaises(SystemExit) as ctx:
                fc._merge_cuopt_config(self.gams_dir)

        self.assertEqual(ctx.exception.code, 1)
        # Depending on the engine and the exact input, rejection happens
        # either while parsing (real YAML syntax error / unsupported
        # construct) or afterwards, when the parsed document turns out not
        # to be a mapping of sections - either is an acceptable safe outcome.
        message = out.getvalue()
        self.assertTrue(
            "not a valid YAML file" in message or "must map configuration sections" in message,
            f"unexpected message: {message!r}",
        )
        # The file must be byte-for-byte unchanged: no cuopt block written into it.
        self.assertEqual(self.read(fc.CONFIG_FILE), original)

    # These cases are genuinely invalid YAML, so they must be rejected
    # whichever engine is in effect - real pyyaml when available, or the
    # dependency-free fallback parser otherwise.
    def test_unterminated_flow_sequence_is_rejected(self):
        self.assert_rejected("other:\n  malformed: [unterminated\n", force_fallback=False)

    def test_unterminated_flow_sequence_is_rejected_by_fallback(self):
        self.assert_rejected("other:\n  malformed: [unterminated\n", force_fallback=True)

    def test_tab_indentation_is_rejected(self):
        self.assert_rejected("other:\n\tmalformed: 1\n", force_fallback=False)

    def test_tab_indentation_is_rejected_by_fallback(self):
        self.assert_rejected("other:\n\tmalformed: 1\n", force_fallback=True)

    def test_unterminated_quote_is_rejected(self):
        self.assert_rejected('other: "unterminated\n', force_fallback=False)

    def test_unterminated_quote_is_rejected_by_fallback(self):
        self.assert_rejected('other: "unterminated\n', force_fallback=True)

    def test_scalar_root_document_is_rejected(self):
        self.assert_rejected("just a scalar document\n", force_fallback=False)

    def test_scalar_root_document_is_rejected_by_fallback(self):
        self.assert_rejected("just a scalar document\n", force_fallback=True)

    def test_flow_style_mapping_is_rejected_by_fallback_parser(self):
        # Flow collections ("{...}") are valid YAML but outside the
        # conservative subset the dependency-free fallback parser supports,
        # so it must refuse rather than guess. (When real pyyaml is
        # available it is used instead and correctly accepts this - see
        # PyyamlIntegrationTest below.)
        self.assert_rejected("other: {a: 1}\n", force_fallback=True)


@unittest.skipIf(fc._pyyaml is None, "pyyaml is not installed in this environment")
class PyyamlIntegrationTest(TempGamsDirMixin, unittest.TestCase):
    """When pyyaml is installed, `_yaml_safe_load` must use it instead of the fallback parser."""

    def test_yaml_safe_load_uses_pyyaml_when_available(self):
        with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as file:
            file.write("solverConfig:\n  cuopt:\n    minVersion: 1\n")
            path = file.name
        self.addCleanup(os.unlink, path)

        with mock.patch.object(fc._pyyaml, "safe_load", wraps=fc._pyyaml.safe_load) as spy:
            result = fc._yaml_safe_load(path)

        spy.assert_called_once()
        self.assertEqual(result, {"solverConfig": {"cuopt": {"minVersion": 1}}})

    def test_pyyaml_syntax_errors_are_translated_to_yamlerror(self):
        with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as file:
            file.write("other:\n  malformed: [unterminated\n")
            path = file.name
        self.addCleanup(os.unlink, path)

        with self.assertRaises(fc._YamlError):
            fc._yaml_safe_load(path)

    def test_flow_style_mapping_is_accepted_via_pyyaml(self):
        # Full YAML support (flow collections included) becomes available
        # automatically once pyyaml is installed, unlike the dependency-free
        # fallback parser exercised in MalformedYamlRejectionTest above.
        self.write(fc.CONFIG_FILE, "other: {a: 1}\n")
        self.write(fc.CUOPT_CONFIG_FILE, CUOPT_RELEASE_YAML)

        with contextlib.redirect_stdout(io.StringIO()):
            fc._merge_cuopt_config(self.gams_dir)

        merged = self.read(fc.CONFIG_FILE)
        self.assertIn("other: {a: 1}", merged)
        self.assertIn("cuopt:", merged)

    def test_quoted_key_merge_behavior_matches_across_engines(self):
        # The duplicate-registration fix (review finding #1) must hold
        # whichever engine parsed the file.
        original = 'solverConfig:\n  - "cuopt":\n      minVersion: 1\n'
        self.write(fc.CONFIG_FILE, original)
        self.write(fc.CUOPT_CONFIG_FILE, CUOPT_RELEASE_YAML)

        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            fc._merge_cuopt_config(self.gams_dir)

        self.assertIn("already exists", out.getvalue())
        self.assertEqual(self.read(fc.CONFIG_FILE), original)


class NormalMergeTest(TempGamsDirMixin, unittest.TestCase):
    """Supported configuration merge behavior must keep working."""

    def test_merge_into_fresh_directory_with_no_existing_config(self):
        self.write(fc.CUOPT_CONFIG_FILE, CUOPT_RELEASE_YAML)

        with contextlib.redirect_stdout(io.StringIO()):
            fc._merge_cuopt_config(self.gams_dir)

        merged = self.read(fc.CONFIG_FILE)
        self.assertIn(fc.MARKER_BEGIN, merged)
        self.assertIn(fc.MARKER_END, merged)
        self.assertIn("solverConfig:", merged)
        self.assertIn("cuopt:", merged)
        self.assertFalse(os.path.exists(os.path.join(self.gams_dir, fc.CUOPT_CONFIG_FILE)))

        # Parses back as valid, supported YAML with cuopt registered.
        config = fc._load_config(os.path.join(self.gams_dir, fc.CONFIG_FILE))
        solver_config = config["solverConfig"]
        self.assertTrue(any("cuopt" in item for item in solver_config))

    def test_merge_preserves_unrelated_existing_sections(self):
        original = (
            "commandLineParameters:\n"
            "  - LP:\n"
            "      value: HIGHS\n"
        )
        self.write(fc.CONFIG_FILE, original)
        self.write(fc.CUOPT_CONFIG_FILE, CUOPT_RELEASE_YAML)

        with contextlib.redirect_stdout(io.StringIO()):
            fc._merge_cuopt_config(self.gams_dir)

        merged = self.read(fc.CONFIG_FILE)
        self.assertIn("commandLineParameters:", merged)
        self.assertIn("cuopt:", merged)

    def test_merge_is_idempotent_on_second_run(self):
        self.write(fc.CUOPT_CONFIG_FILE, CUOPT_RELEASE_YAML)
        with contextlib.redirect_stdout(io.StringIO()):
            fc._merge_cuopt_config(self.gams_dir)

        self.write(fc.CUOPT_CONFIG_FILE, CUOPT_RELEASE_YAML)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            fc._merge_cuopt_config(self.gams_dir)

        self.assertIn("already exists", out.getvalue())
        merged = self.read(fc.CONFIG_FILE)
        self.assertEqual(merged.count(fc.MARKER_BEGIN), 1)

    def test_remove_after_merge_restores_previous_content(self):
        original = "commandLineParameters:\n  - LP:\n      value: HIGHS\n"
        self.write(fc.CONFIG_FILE, original)
        self.write(fc.CUOPT_CONFIG_FILE, CUOPT_RELEASE_YAML)

        with contextlib.redirect_stdout(io.StringIO()):
            fc._merge_cuopt_config(self.gams_dir)
        with contextlib.redirect_stdout(io.StringIO()):
            fc._remove_cuopt_from_config(self.gams_dir)

        self.assertEqual(self.read(fc.CONFIG_FILE), original)


class TimeoutErrorHandlingTest(unittest.TestCase):
    """Review finding #3: a direct TimeoutError must be handled cleanly, not crash."""

    def test_get_asset_urls_handles_direct_timeout_error(self):
        with mock.patch.object(fc.urllib.request, "urlopen", side_effect=TimeoutError("timed out")):
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                with self.assertRaises(SystemExit) as ctx:
                    fc._get_asset_urls(["some-asset.zip"])

        self.assertEqual(ctx.exception.code, 1)
        self.assertIn("Could not reach GitHub API", out.getvalue())
        self.assertIn("timed out", out.getvalue())

    def test_download_handles_direct_timeout_error(self):
        with mock.patch.object(fc.urllib.request, "urlopen", side_effect=TimeoutError("timed out")):
            out = io.StringIO()
            with tempfile.TemporaryDirectory() as tmp:
                target = os.path.join(tmp, "asset.zip")
                with contextlib.redirect_stdout(out):
                    with self.assertRaises(SystemExit) as ctx:
                        fc._download("https://example.invalid/asset.zip", target)

        self.assertEqual(ctx.exception.code, 1)
        self.assertIn("Could not download", out.getvalue())
        self.assertIn("timed out", out.getvalue())

    def test_get_asset_urls_still_reports_http_errors(self):
        import urllib.error

        http_error = urllib.error.HTTPError(
            url="https://api.github.com/x", code=404, msg="Not Found", hdrs=None,
            fp=io.BytesIO(b"not found"),
        )
        with mock.patch.object(fc.urllib.request, "urlopen", side_effect=http_error):
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                with self.assertRaises(SystemExit):
                    fc._get_asset_urls(["some-asset.zip"])

        self.assertIn("404", out.getvalue())


class CliCompatibilityTest(unittest.TestCase):
    """Basic CLI wiring must still work, without touching the network or a real GAMS install."""

    def test_install_help(self):
        parser = fc._build_parser()
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            with self.assertRaises(SystemExit) as ctx:
                parser.parse_args(["install", "--help"])
        self.assertEqual(ctx.exception.code, 0)
        text = out.getvalue()
        self.assertIn("--gams-dir", text)
        self.assertIn("--cuda-version", text)
        self.assertIn("--cuda-runtime", text)
        self.assertIn("--release", text)

    def test_uninstall_help(self):
        parser = fc._build_parser()
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            with self.assertRaises(SystemExit) as ctx:
                parser.parse_args(["uninstall", "--help"])
        self.assertEqual(ctx.exception.code, 0)
        self.assertIn("--gams-dir", out.getvalue())

    def test_no_command_triggers_interactive_install(self):
        with mock.patch.object(fc, "_run_interactive_install") as interactive:
            with mock.patch.object(sys, "argv", ["fetch-cuoptlink.py"]):
                fc.main()
        interactive.assert_called_once_with()

    def test_install_subcommand_dispatches_to_execute_install(self):
        with mock.patch.object(fc, "_execute_install") as execute:
            with mock.patch.object(
                sys, "argv", ["fetch-cuoptlink.py", "install", "-g", "/tmp/gamsdir"]
            ):
                fc.main()
        execute.assert_called_once_with("/tmp/gamsdir", None, False, "latest")


if __name__ == "__main__":
    unittest.main()
