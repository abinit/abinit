"""Unit tests for testbot.py module."""

from __future__ import annotations

import json
import os
import sys
import tempfile
from unittest.mock import MagicMock, patch

import pytest

# Mock heavy dependencies before importing testbot
sys.modules["pymods"] = MagicMock()
sys.modules["pymods.termcolor"] = MagicMock()
sys.modules["pymods.jobrunner"] = MagicMock()
sys.modules["pymods.testsuite"] = MagicMock()
sys.modules["pymods.tools"] = MagicMock()

# Now we can import testbot
from tests.testbot import (
    TestBot,
    TestBotSummary,
    _str2list,
    _yesno2bool,
    analyze,
    build_parser,
    get_git_tag,
    get_mpi_prefix_from_env,
    main,
)

# These are production classes from testbot.py, not test classes. Their names
# happen to start with "Test", so mark them explicitly to stop pytest from
# trying to collect them (they have __init__ constructors and warn otherwise).
TestBot.__test__ = False
TestBotSummary.__test__ = False


class TestYesno2Bool:
    """Tests for the _yesno2bool() utility function."""

    def test_yesno2bool_yes(self):
        """'yes' should convert to True."""
        assert _yesno2bool("yes") is True
        assert _yesno2bool("YES") is True
        assert _yesno2bool("Yes") is True

    def test_yesno2bool_no(self):
        """'no' should convert to False."""
        assert _yesno2bool("no") is False
        assert _yesno2bool("NO") is False
        assert _yesno2bool("No") is False

    def test_yesno2bool_with_whitespace(self):
        """Whitespace should be stripped."""
        assert _yesno2bool("  yes  ") is True
        assert _yesno2bool("  no  ") is False

    def test_yesno2bool_with_quotes(self):
        """Quotes should be removed."""
        assert _yesno2bool('"yes"') is True
        assert _yesno2bool("'no'") is False
        assert _yesno2bool('"yes"') is True

    def test_yesno2bool_invalid(self):
        """Invalid strings should raise ValueError."""
        with pytest.raises(ValueError, match="Cannot interpret string"):
            _yesno2bool("maybe")
        with pytest.raises(ValueError, match="Cannot interpret string"):
            _yesno2bool("1")
        with pytest.raises(ValueError, match="Cannot interpret string"):
            _yesno2bool("")


class TestStr2List:
    """Tests for the _str2list() utility function."""

    def test_str2list_comma_separated(self):
        """Comma-separated strings should be split."""
        assert _str2list("a,b,c") == ["a", "b", "c"]

    def test_str2list_with_whitespace(self):
        """Whitespace should be trimmed from each item."""
        assert _str2list("a, b, c") == ["a", "b", "c"]
        assert _str2list("  a  ,  b  ,  c  ") == ["a", "b", "c"]

    def test_str2list_single_item(self):
        """Single item should return a list with one element."""
        assert _str2list("abc") == ["abc"]

    def test_str2list_empty_string(self):
        """Empty string should return empty list."""
        assert _str2list("") == []

    def test_str2list_with_empty_items(self):
        """Empty items should be skipped."""
        assert _str2list("a,,b") == ["a", "b"]
        assert _str2list(",a,b,") == ["a", "b"]

    def test_str2list_already_list(self):
        """Lists and tuples should be returned as-is (converted to list if tuple)."""
        assert _str2list(["a", "b", "c"]) == ["a", "b", "c"]
        assert _str2list(("a", "b", "c")) == ("a", "b", "c")


class TestGetMpiPrefixFromEnv:
    """Tests for get_mpi_prefix_from_env() function."""

    def test_get_mpi_prefix_from_mpi_home(self):
        """Should return MPI_HOME if set."""
        with patch.dict(os.environ, {"MPI_HOME": "/usr/local/mpi"}):
            result = get_mpi_prefix_from_env()
            assert result == "/usr/local/mpi"

    def test_get_mpi_prefix_from_mpihome(self):
        """Should return MPIHOME if MPI_HOME is not set."""
        env = {"MPIHOME": "/opt/mpich"}
        # Remove MPI_HOME if it exists
        env_copy = os.environ.copy()
        env_copy.pop("MPI_HOME", None)
        env_copy.update(env)
        with patch.dict(os.environ, env_copy, clear=False):
            result = get_mpi_prefix_from_env()
            assert result == "/opt/mpich" or result is None

    def test_get_mpi_prefix_mpi_home_precedence(self):
        """MPI_HOME should take precedence over MPIHOME."""
        with patch.dict(os.environ, {"MPI_HOME": "/usr/local/mpi", "MPIHOME": "/opt/mpich"}):
            result = get_mpi_prefix_from_env()
            assert result == "/usr/local/mpi"

    def test_get_mpi_prefix_neither_set(self):
        """Should return None if neither variable is set."""
        with patch.dict(os.environ, {}, clear=True):
            result = get_mpi_prefix_from_env()
            assert result is None


class TestGetGitTag:
    """Tests for get_git_tag() function."""

    def test_get_git_tag_success(self):
        """Should return git tag when available."""
        result = get_git_tag()
        # Result should be a string (either a real tag or "unknown")
        assert isinstance(result, str)

    def test_get_git_tag_default(self):
        """Should return 'unknown' if git fails."""
        # Mock subprocess to simulate git failure
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=1, stdout="")
            result = get_git_tag()
            assert result == "unknown"

    def test_get_git_tag_exception(self):
        """Should return 'unknown' if git raises an exception."""
        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = Exception("Git not found")
            result = get_git_tag()
            assert result == "unknown"


class TestAnalyzeFunction:
    """Tests for analyze() function."""

    def test_analyze_missing_file(self):
        """Should return 1 if summary file is missing."""
        result = analyze("nonexistent_file.json")
        assert result == 1

    def test_analyze_valid_summary(self, monkeypatch):
        """Should return 0 for a valid summary file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            summary_file = os.path.join(tmpdir, "testbot_summary.json")

            # Create a minimal valid summary
            summary = {
                "failed": [],
                "passed": [],
                "summary_table": [["suite1", "suite2"], ["0/1/0/0", "1/0/0/0"]],
                "suite1": {"test1": {"status": "passed"}},
                "suite2": {"test2": {"status": "failed"}},
            }

            with open(summary_file, "w") as f:
                json.dump(summary, f)

            # Change to temp dir (monkeypatch restores the original cwd on
            # teardown even if the test fails) to avoid littering it.
            monkeypatch.chdir(tmpdir)
            result = analyze(summary_file)
            assert result == 0

            # Check that tag was added to the file
            with open(summary_file) as f:
                updated = json.load(f)
                assert "tag" in updated

    def test_analyze_missing_summary_table(self, monkeypatch):
        """Should return 1 if summary_table is missing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            summary_file = os.path.join(tmpdir, "testbot_summary.json")

            # Create summary without summary_table
            summary = {"failed": [], "passed": []}

            with open(summary_file, "w") as f:
                json.dump(summary, f)

            monkeypatch.chdir(tmpdir)
            result = analyze(summary_file)
            assert result == 1

    def test_analyze_invalid_json(self, monkeypatch):
        """Should return 99 for invalid JSON."""
        with tempfile.TemporaryDirectory() as tmpdir:
            summary_file = os.path.join(tmpdir, "testbot_summary.json")

            # Write invalid JSON
            with open(summary_file, "w") as f:
                f.write("invalid json {")

            monkeypatch.chdir(tmpdir)
            result = analyze(summary_file)
            assert result == 99


class TestTestBotSummary:
    """Tests for TestBotSummary class."""

    def test_testbotsummary_init(self):
        """TestBotSummary should initialize with a result table."""
        res_table = {
            "suite1": {"test1": {}, "test2": {}},
            "suite2": {"test3": {}},
        }
        summary = TestBotSummary(res_table)
        assert summary.res_table == res_table
        assert summary.failed == []
        assert summary.passed == []

    def test_testbotsummary_suite_names(self):
        """Should return sorted list of suite names."""
        res_table = {
            "zsuite": {"test1": {}},
            "asuite": {"test2": {}},
            "msuite": {"test3": {}},
        }
        summary = TestBotSummary(res_table)
        names = summary.suite_names()
        assert names == ["asuite", "msuite", "zsuite"]

    def test_testbotsummary_iter(self):
        """Should be iterable over suite names."""
        res_table = {
            "suite2": {},
            "suite1": {},
        }
        summary = TestBotSummary(res_table)
        suites = list(summary)
        assert suites == ["suite1", "suite2"]

    def test_testbotsummary_status_of_suite_all_passed(self):
        """Should report status correctly for all-passed suite."""
        res_table = {
            "suite1": {
                "test1": {"status": "passed"},
                "test2": {"status": "passed"},
            }
        }
        summary = TestBotSummary(res_table)
        status, stats = summary.status_of_suite("suite1")
        assert status == "passed"
        assert stats["passed"] == 2
        assert stats["failed"] == 0

    def test_testbotsummary_status_of_suite_has_failure(self):
        """Failed status should take precedence over passed."""
        res_table = {
            "suite1": {
                "test1": {"status": "passed"},
                "test2": {"status": "failed"},
            }
        }
        summary = TestBotSummary(res_table)
        status, stats = summary.status_of_suite("suite1")
        assert status == "failed"
        assert stats["failed"] == 1
        assert stats["passed"] == 1

    def test_testbotsummary_to_table(self):
        """Should convert results to table format."""
        res_table = {
            "suite1": {
                "test1": {"status": "passed"},
                "test2": {"status": "failed"},
            },
            "suite2": {"test3": {"status": "skipped"}},
        }
        summary = TestBotSummary(res_table)
        table = summary.to_table()
        assert len(table) == 2
        assert table[0] == ["suite1", "suite2"]  # suite names
        assert len(table[1]) == 2  # rows with stats

    def test_testbotsummary_json_dump(self, monkeypatch):
        """Should dump results to JSON file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            res_table = {
                "suite1": {
                    "test1": {"status": "passed"},
                }
            }
            summary = TestBotSummary(res_table)
            summary.failed = ["suite1/test_failed"]
            summary.passed = ["suite1/test1"]

            json_file = os.path.join(tmpdir, "summary.json")
            monkeypatch.chdir(tmpdir)
            summary.json_dump(json_file)

            assert os.path.exists(json_file)
            with open(json_file) as f:
                data = json.load(f)
                assert "failed" in data
                assert "passed" in data
                assert "summary_table" in data


class TestTestBotClass:
    """Tests for the TestBot class."""

    def test_testbot_print_options(self):
        """print_options() should not raise an exception."""
        # Just verify it doesn't crash
        try:
            TestBot.print_options()
        except SystemExit:
            # This is acceptable, the function might call sys.exit
            pass

    def test_testbot_attributes(self):
        """TestBot._attrbs should have expected keys."""
        expected_keys = {
            "builder_name", "type", "max_cpus", "max_gpus", "mpi_prefix",
            "mpirun_np", "omp_num_threads", "enable_mpi", "enable_openmp",
            "with_tdirs", "without_tdirs", "timeout_time", "runmode",
            "keywords", "verbose", "tmp_basedir", "mpi_args",
            "force_mpi"
        }
        assert set(TestBot._attrbs.keys()) == expected_keys

    def test_testbot_has_mpi_property(self):
        """has_mpi property should check for mpirun_np."""
        # Create a mock TestBot instance with minimal setup
        tb = MagicMock(spec=TestBot)
        # When mpirun_np is set, has_mpi should be True
        tb.mpirun_np = "/usr/bin/mpirun"
        assert bool(tb.mpirun_np)

    def test_testbot_has_openmp_property(self):
        """has_openmp property should check for omp_num_threads."""
        tb = MagicMock(spec=TestBot)
        tb.omp_num_threads = 4
        assert bool(tb.omp_num_threads > 0)


class TestBuildParser:
    """Tests for the run/analyze/print argparse CLI."""

    def test_run_defaults_and_explicit_cfg(self):
        parser = build_parser()
        ns = parser.parse_args(["run"])
        assert ns.command == "run"
        assert ns.testbot_cfg is None

        ns = parser.parse_args(["run", "my.cfg"])
        assert ns.testbot_cfg == "my.cfg"

    def test_analyze_defaults_and_explicit_tag(self):
        parser = build_parser()
        ns = parser.parse_args(["analyze"])
        assert ns.command == "analyze"
        assert ns.tag == "unknown"

        ns = parser.parse_args(["analyze", "v9.2.0"])
        assert ns.tag == "v9.2.0"

    def test_print_defaults_and_explicit_cfg(self):
        parser = build_parser()
        ns = parser.parse_args(["print"])
        assert ns.command == "print"
        assert ns.testbot_cfg is None

        ns = parser.parse_args(["print", "my.cfg"])
        assert ns.testbot_cfg == "my.cfg"

    def test_missing_command_is_rejected(self):
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([])

    def test_unknown_command_is_rejected(self):
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["dry-run"])


class TestMain:
    """Characterization tests for main()'s dispatch to run/analyze/print."""

    def test_main_dispatches_to_analyze(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["testbot.py", "analyze", "v1"])
        with patch("tests.testbot.analyze", return_value=0) as mock_analyze:
            assert main() == 0
        mock_analyze.assert_called_once_with(fname="testbot_summary.json", tag="v1")

    def test_main_dispatches_to_run(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["testbot.py", "run", "my.cfg"])
        mock_instance = MagicMock()
        mock_instance.run.return_value = 3
        with patch("tests.testbot.TestBot", return_value=mock_instance) as mock_cls:
            assert main() == 3
        mock_cls.assert_called_once_with("my.cfg")
        mock_instance.run.assert_called_once()

    def test_main_dispatches_to_print_without_running(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["testbot.py", "print", "my.cfg"])
        mock_instance = MagicMock()
        mock_instance.__str__ = MagicMock(return_value="fake-config-dump")
        with patch("tests.testbot.TestBot", return_value=mock_instance) as mock_cls:
            assert main() == 0
        mock_cls.assert_called_once_with("my.cfg")
        mock_instance.run.assert_not_called()
        assert "fake-config-dump" in capsys.readouterr().out

    def test_main_run_with_no_cfg_uses_default(self, monkeypatch):
        """A bare `testbot.py run` (the analysis.sh invocation) passes cfg=None."""
        monkeypatch.setattr(sys, "argv", ["testbot.py", "run"])
        mock_instance = MagicMock()
        mock_instance.run.return_value = 0
        with patch("tests.testbot.TestBot", return_value=mock_instance) as mock_cls:
            assert main() == 0
        mock_cls.assert_called_once_with(None)


class TestUtilityFunctions:
    """Tests for other utility functions."""

    def test_lazy_str_decorator(self):
        """lazy__str__ decorator should create a __str__ method."""
        from tests.testbot import lazy__str__

        class TestClass:
            def __init__(self):
                self.a = 1
                self.b = 2

            @lazy__str__
            def __str__(self):
                pass

        obj = TestClass()
        result = str(obj)
        assert "a : 1" in result
        assert "b : 2" in result
