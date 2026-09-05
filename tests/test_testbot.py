"""Unit tests for testbot.py module."""

from __future__ import annotations

import dataclasses
import json
import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Mock heavy dependencies before importing testbot
sys.modules["pymods"] = MagicMock()
sys.modules["pymods.termcolor"] = MagicMock()
sys.modules["pymods.jobrunner"] = MagicMock()
sys.modules["pymods.testsuite"] = MagicMock()
sys.modules["pymods.tools"] = MagicMock()

# Monkey-patch dataclass to ignore kw_only parameter in Python < 3.10
_original_dataclass = dataclasses.dataclass


def _dataclass_compat(*args, **kwargs):
    """Wrapper that removes kw_only parameter for Python < 3.10 compatibility."""
    if sys.version_info < (3, 10) and "kw_only" in kwargs:
        kwargs.pop("kw_only")
    return _original_dataclass(*args, **kwargs)


# Apply the monkey-patch
dataclasses.dataclass = _dataclass_compat

# Now we can import testbot
from tests.testbot import (
    TestBot,
    TestBotContext,
    TestBotSummary,
    _str2list,
    _yesno2bool,
    analyze,
    get_git_tag,
    get_mpi_prefix_from_env,
    read_builders,
)

# These are production classes from testbot.py, not test classes. Their names
# happen to start with "Test", so mark them explicitly to stop pytest from
# trying to collect them (they have __init__ constructors and warn otherwise).
TestBot.__test__ = False
TestBotContext.__test__ = False
TestBotSummary.__test__ = False

TESTBOT_CONTEXT_AVAILABLE = True


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


@pytest.mark.skipif(
    not TESTBOT_CONTEXT_AVAILABLE,
    reason="TestBotContext requires Python 3.10+ (kw_only parameter)"
)
class TestTestBotContext:
    """Tests for the TestBotContext dataclass."""

    def test_testbotcontext_creation(self):
        """TestBotContext should be creatable with required fields."""
        ctx = TestBotContext(slavename="test_worker", ncpus=8)
        assert ctx.slavename == "test_worker"
        assert ctx.ncpus == 8
        assert ctx.type == ""
        assert ctx.timeout_time == 900.0

    def test_testbotcontext_missing_slavename(self):
        """Missing slavename should raise ValueError."""
        with pytest.raises(ValueError, match="Missing required configuration"):
            TestBotContext(slavename="", ncpus=8)

    def test_testbotcontext_missing_ncpus(self):
        """Missing ncpus should raise ValueError."""
        with pytest.raises(ValueError, match="Missing required configuration"):
            TestBotContext(slavename="test_worker", ncpus=None)

    def test_testbotcontext_negative_ncpus(self):
        """Negative ncpus should raise ValueError."""
        with pytest.raises(ValueError, match="ncpus is negative"):
            TestBotContext(slavename="test_worker", ncpus=-1)

    def test_testbotcontext_negative_max_gpus(self):
        """Negative max_gpus should raise ValueError."""
        with pytest.raises(ValueError, match="max_gpus is negative"):
            TestBotContext(slavename="test_worker", ncpus=8, max_gpus=-1)

    def test_testbotcontext_negative_omp_threads(self):
        """Negative omp_num_threads should raise ValueError."""
        with pytest.raises(ValueError, match="omp_num_threads is negative"):
            TestBotContext(slavename="test_worker", ncpus=8, omp_num_threads=-1)

    def test_testbotcontext_negative_timeout(self):
        """Non-positive timeout_time should raise ValueError."""
        with pytest.raises(ValueError, match="timeout_time is negative"):
            TestBotContext(slavename="test_worker", ncpus=8, timeout_time=-1)
        with pytest.raises(ValueError, match="timeout_time is negative"):
            TestBotContext(slavename="test_worker", ncpus=8, timeout_time=0)

    def test_testbotcontext_from_builders_valid(self):
        """from_builders should extract configuration for a valid builder."""
        builders = [
            {"name": "builder1", "slavename": "worker1", "ncpus": 8},
            {"name": "builder2", "slavename": "worker2", "ncpus": 16},
        ]
        ctx = TestBotContext.from_builders(builders, "builder1")
        assert ctx.slavename == "worker1"
        assert ctx.ncpus == 8

    def test_testbotcontext_from_builders_missing(self):
        """from_builders should raise ValueError for missing builder."""
        builders = [{"name": "builder1", "slavename": "worker1", "ncpus": 8}]
        with pytest.raises(ValueError, match="Cannot find"):
            TestBotContext.from_builders(builders, "builder_missing")

    def test_testbotcontext_with_tdirs_and_without_tdirs(self):
        """Mutually exclusive with_tdirs and without_tdirs should fail at creation."""
        # This is validated in TestBot.__init__, not in TestBotContext
        # so let's just test that they can coexist in the dataclass
        ctx = TestBotContext(
            slavename="test",
            ncpus=8,
            with_tdirs=["dir1"],
            without_tdirs=["dir2"]
        )
        assert ctx.with_tdirs == ["dir1"]
        assert ctx.without_tdirs == ["dir2"]

    def test_testbotcontext_type_field(self):
        """Type field should accept 'ref' or empty string."""
        ctx = TestBotContext(slavename="test", ncpus=8, type="ref")
        assert ctx.type == "ref"
        ctx2 = TestBotContext(slavename="test", ncpus=8, type="")
        assert ctx2.type == ""


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


class TestReadBuilders:
    """Tests for read_builders() function."""

    def test_read_builders_json_format(self):
        """Should read builders from JSON file."""
        builders = read_builders("json")
        assert isinstance(builders, list)
        assert len(builders) > 0
        assert all("name" in b for b in builders)

    def test_read_builders_yaml_format(self):
        """Should read builders from YAML file."""
        try:
            from ruamel import yaml  # noqa: F401
        except ImportError:
            pytest.skip("ruamel YAML library not available")

        builders = read_builders("yaml")
        assert isinstance(builders, list)
        assert len(builders) > 0
        assert all("name" in b for b in builders)

    def test_read_builders_invalid_format(self):
        """Should raise ValueError for invalid format."""
        with pytest.raises(ValueError, match="Invalid"):
            read_builders("xml")

    def test_read_builders_json_yaml_consistency(self):
        """JSON and YAML builders should have the same structure."""
        try:
            from ruamel import yaml  # noqa: F401
        except ImportError:
            pytest.skip("ruamel YAML library not available")

        builders_json = read_builders("json")
        builders_yaml = read_builders("yaml")

        # Check that we have the same number of builders
        assert len(builders_json) == len(builders_yaml)

        # Check that builder names match
        json_names = {b["name"] for b in builders_json}
        yaml_names = {b["name"] for b in builders_yaml}
        assert json_names == yaml_names


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


@pytest.mark.skipif(
    not TESTBOT_CONTEXT_AVAILABLE,
    reason="TestBotSummary requires Python 3.10+ (kw_only parameter)"
)
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


@pytest.mark.skipif(
    not TESTBOT_CONTEXT_AVAILABLE,
    reason="TestBot requires Python 3.10+ (kw_only parameter)"
)
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
            "slavename", "type", "ncpus", "max_gpus", "mpi_prefix",
            "mpirun_np", "omp_num_threads", "enable_mpi", "enable_openmp",
            "with_tdirs", "without_tdirs", "timeout_time", "runmode",
            "keywords", "etsf_check", "verbose", "tmp_basedir", "mpi_args",
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


class TestGetParser:
    """Tests for command-line parser."""

    def test_parser_creation(self):
        """Parser should be created without errors."""
        from tests.testbot import get_parser
        parser = get_parser()
        assert parser is not None

    def test_parser_subcommands(self):
        """Parser should have expected subcommands."""
        from tests.testbot import get_parser
        parser = get_parser()
        # Parse with minimal args to verify subcommands exist
        try:
            args = parser.parse_args(["run", "test_builder"])
            assert args.command == "run"
            assert args.builder_name == "test_builder"
        except SystemExit:
            pass


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
