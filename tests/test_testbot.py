"""Unit tests for testbot.py module."""

from __future__ import annotations

import dataclasses
import json
import os
import sys
import tempfile
import warnings
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
    BenchmarkResult,
    TestBot,
    TestBotSummary,
    TestRunSummary,
    analyze,
    benchmark,
    build_parser,
    default_mpirun_np,
    get_git_tag,
    get_mpi_prefix_from_env,
    main,
)

# These are production classes from testbot.py, not test classes. Their names
# happen to start with "Test", so mark them explicitly to stop pytest from
# trying to collect them (they have __init__ constructors and warn otherwise).
TestBot.__test__ = False
BenchmarkResult.__test__ = False
TestRunSummary.__test__ = False
TestBotSummary.__test__ = False


class TestGetMpiPrefixFromEnv:
    """Tests for get_mpi_prefix_from_env() function."""

    def test_get_mpi_prefix_from_mpi_home(self):
        """Should return MPI_HOME if set."""
        with patch.dict(os.environ, {"MPI_HOME": "/usr/local/mpi"}):
            result = get_mpi_prefix_from_env()
            assert result == "/usr/local/mpi"

    def test_get_mpi_prefix_from_mpihome(self):
        """Should return MPIHOME if MPI_HOME is not set.

        Regression test: this used to build a *copy* of os.environ with
        MPI_HOME popped out, then apply it via patch.dict(..., clear=False)
        -- but clear=False only adds/overwrites keys, it never removes ones
        merely absent from the given dict, so a real MPI_HOME already set on
        the host (common on an HPC worker) leaked straight through and beat
        the injected MPIHOME. clear=True (the pattern every other test in
        this class already uses) makes the environment fully deterministic.
        """
        with patch.dict(os.environ, {"MPIHOME": "/opt/mpich"}, clear=True):
            result = get_mpi_prefix_from_env()
            assert result == "/opt/mpich"

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

    def test_get_mpi_prefix_from_ebrootopenmpi(self):
        """Should fall back to EBROOTOPENMPI, joined with 'bin', if MPI_HOME/MPIHOME unset."""
        with patch.dict(os.environ, {"EBROOTOPENMPI": "/eb/OpenMPI/4.1.6"}, clear=True):
            assert get_mpi_prefix_from_env() == "/eb/OpenMPI/4.1.6/bin"

    def test_get_mpi_prefix_ebroot_precedence_order(self):
        """EBROOTOPENMPI should win over EBROOTMPICH/EBROOTIMPI if multiple are set."""
        env = {"EBROOTOPENMPI": "/eb/ompi", "EBROOTMPICH": "/eb/mpich", "EBROOTIMPI": "/eb/impi"}
        with patch.dict(os.environ, env, clear=True):
            assert get_mpi_prefix_from_env() == "/eb/ompi/bin"

    def test_get_mpi_prefix_mpi_home_beats_ebroot(self):
        """MPI_HOME/MPIHOME must take precedence over any EBROOT* variable."""
        env = {"MPI_HOME": "/usr/local/mpi", "EBROOTMPICH": "/eb/mpich"}
        with patch.dict(os.environ, env, clear=True):
            assert get_mpi_prefix_from_env() == "/usr/local/mpi"


class TestDefaultMpirunNp:
    """Tests for default_mpirun_np(), the fallback when mpirun_np is unset."""

    @staticmethod
    def _make_launcher(directory, name="mpiexec"):
        """Create an executable stub named `name` in `directory` and return its path."""
        os.makedirs(directory, exist_ok=True)
        path = os.path.join(directory, name)
        with open(path, "w") as fh:
            fh.write("#!/bin/sh\n")
        os.chmod(path, 0o700)
        return path

    def test_prefers_launcher_under_mpi_home_prefix(self, tmp_path):
        """MPI_HOME-style prefix: the launcher in <prefix>/bin must win over $PATH."""
        exe = self._make_launcher(os.path.join(str(tmp_path), "bin"))
        assert default_mpirun_np(str(tmp_path)) == f"{exe} -n"

    def test_accepts_ebroot_style_prefix_already_ending_in_bin(self, tmp_path):
        """get_mpi_prefix_from_env returns <root>/bin for EBROOT*; that must resolve too."""
        exe = self._make_launcher(str(tmp_path))
        assert default_mpirun_np(str(tmp_path)) == f"{exe} -n"

    def test_ignores_a_non_executable_candidate(self, tmp_path):
        """A non-executable file named mpiexec must not be picked as the launcher."""
        bindir = tmp_path / "bin"
        bindir.mkdir()
        (bindir / "mpiexec").write_text("not executable\n")
        assert default_mpirun_np(str(tmp_path)) == "mpiexec -n"

    def test_falls_back_to_path_lookup_without_a_prefix(self, monkeypatch):
        """With no usable prefix, fall back to whichever launcher $PATH provides."""
        import tests.testbot as tb_module

        monkeypatch.setattr(
            tb_module.shutil, "which",
            lambda name: None if name == "mpiexec" else "/usr/bin/mpirun",
        )
        assert default_mpirun_np("") == "mpirun -n"

    def test_returns_standard_spelling_when_nothing_is_found(self, monkeypatch):
        """No prefix and nothing in $PATH: still emit a real command name.

        A bogus "mpiexec: command not found" is far easier to diagnose than the
        mangled command line an empty mpirun_np used to produce.
        """
        import tests.testbot as tb_module

        monkeypatch.setattr(tb_module.shutil, "which", lambda name: None)
        assert default_mpirun_np(None) == "mpiexec -n"


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

    def test_analyze_valid_summary(self, monkeypatch, capsys):
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

            output = capsys.readouterr().out
            update_message = (
                f"Updating {summary_file}: adding the Git tag to the merged test results."
            )
            assert update_message in output
            assert (
                "Writing testbot_analysis.html: suite-level status and timing totals "
                "for the results page."
            ) in output

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

    def test_testbotsummary_json_dump(self, monkeypatch, capsys):
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

            assert (
                f"Writing {json_file}: merged test results across all MPI and OpenMP "
                "execution configurations."
            ) in capsys.readouterr().out

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
        """TestBot's dataclass fields should have the expected config keys."""
        expected_keys = {
            "builder_name", "type", "max_cpus", "max_gpus",
            "mpirun_np", "omp_num_threads", "has_mpi",
            "with_tdirs", "without_tdirs", "timeout_time", "runmode",
            "keywords", "verbose", "tmp_basedir", "mpi_args",
            "force_mpi", "dont_exclude_builders", "with_parametrized"
        }
        config_fields = {f.name for f in dataclasses.fields(TestBot) if f.init}
        assert config_fields == expected_keys

    def test_testbot_has_openmp_property(self):
        """has_openmp property should check for omp_num_threads."""
        tb = MagicMock(spec=TestBot)
        tb.omp_num_threads = 4
        assert bool(tb.omp_num_threads > 0)

    def test_write_run_summaries_creates_json_and_html(self, tmp_path, monkeypatch, capsys):
        """Per-run output should be machine-readable and link to detailed reports."""
        monkeypatch.chdir(tmp_path)
        testbot = TestBot.__new__(TestBot)
        testbot.run_summaries = [
            TestRunSummary(
                mpi_nprocs=4,
                omp_nthreads=2,
                py_nprocs=3,
                runmode="static",
                workdir_name="TestBot_MPI4_OMP2",
                nfailed=1,
                npassed=2,
                nsucceeded=7,
                nskipped=4,
                ndisabled=1,
                nexecuted=10,
            )
        ]

        testbot.write_run_summaries()

        assert capsys.readouterr().out.splitlines() == [
            "Writing testbot_runs.json: per-execution test counts and parallel configuration "
            "for machine processing.",
            "Writing testbot_runs.html: per-execution test counts and links to detailed reports "
            "for the results page.",
        ]

        data = json.loads((tmp_path / "testbot_runs.json").read_text())
        assert data == [testbot.run_summaries[0].as_dict()]
        report = (tmp_path / "testbot_runs.html").read_text()
        assert 'class="testbot-runs"' in report
        assert "TestBot_MPI4_OMP2" in report
        assert 'href="TestBot_MPI4_OMP2/"' in report
        assert "<td>10</td><td>1</td><td>2</td><td>7</td>" in report

        # Report must be the second column (right after Configuration), not
        # trailing at the end -- readers scanning left to right want the link
        # to the detailed report before wading through the pass/fail counts.
        assert "<th>Configuration</th><th>Report</th><th>MPI</th>" in report
        assert '<td>TestBot_MPI4_OMP2</td><td><a class="run-report-link"' in report

    def test_write_run_summaries_escapes_workdir(self, tmp_path, monkeypatch):
        """Generated report paths must not permit HTML injection."""
        monkeypatch.chdir(tmp_path)
        testbot = TestBot.__new__(TestBot)
        testbot.run_summaries = [
            TestRunSummary(
                mpi_nprocs=1,
                omp_nthreads=1,
                py_nprocs=1,
                runmode="static",
                workdir_name='TestBot_<script>"',
                nfailed=0,
                npassed=0,
                nsucceeded=1,
                nskipped=0,
                ndisabled=0,
                nexecuted=1,
            )
        ]

        testbot.write_run_summaries()

        report = (tmp_path / "testbot_runs.html").read_text()
        assert "<script>" not in report
        assert "TestBot_&lt;script&gt;&quot;" in report

    def test_write_run_summaries_flags_failed_runs_for_highlighting(self, tmp_path, monkeypatch):
        """A run with nfailed > 0 must get "suite-failed" on its <tr>, reusing
        the same row_class convention to_table() already uses for
        testbot_analysis.html -- so the results page's existing CSS for that
        class (a light, readable red tint) highlights this table's failed
        rows too, with no new styling needed.
        """
        monkeypatch.chdir(tmp_path)
        testbot = TestBot.__new__(TestBot)
        testbot.run_summaries = [
            TestRunSummary(
                mpi_nprocs=1, omp_nthreads=1, py_nprocs=1, runmode="static",
                workdir_name="TestBot_failed", nfailed=3, npassed=0,
                nsucceeded=0, nskipped=0, ndisabled=0, nexecuted=3,
            ),
            TestRunSummary(
                mpi_nprocs=1, omp_nthreads=1, py_nprocs=1, runmode="static",
                workdir_name="TestBot_ok", nfailed=0, npassed=3,
                nsucceeded=0, nskipped=0, ndisabled=0, nexecuted=3,
            ),
        ]

        testbot.write_run_summaries()

        report = (tmp_path / "testbot_runs.html").read_text()
        assert '<tr class="suite-failed"><td>TestBot_failed</td>' in report
        assert '<tr class="suite-ok"><td>TestBot_ok</td>' in report


class TestTestBotFromJson:
    """End-to-end characterization tests for TestBot.from_json()."""

    def _mock_environment(self, monkeypatch, *, defined_cppvars=(), has_timeout=False):
        """Mock the module-level names TestBot.__post_init__ relies on."""
        import tests.testbot as tb_module

        mock_database = MagicMock()
        mock_database.init_result_table.return_value = {}
        monkeypatch.setattr(tb_module.abitests, "get_database", lambda: mock_database)

        mock_build_env = MagicMock()
        mock_build_env.defined_cppvars = list(defined_cppvars)
        mock_build_env.has_bin.return_value = has_timeout
        mock_build_env.path_of_bin.return_value = "/usr/bin/timeout"
        monkeypatch.setattr(tb_module, "BuildEnvironment", lambda *a, **kw: mock_build_env)

        return mock_build_env

    def test_from_json_builds_a_real_instance(self, tmp_path, monkeypatch):
        """Feeding a real testbot.json through from_json() must produce a working TestBot."""
        self._mock_environment(monkeypatch, defined_cppvars=[])
        testbot_json = tmp_path / "testbot.json"
        testbot_json.write_text(json.dumps({
            "builder_name": "eos_gnu_13.2_serial",
            "max_cpus": 4,
            "with_tdirs": ["v1", "v2"],
            "force_mpi": True,
            "has_mpi": False,
        }))

        testbot = TestBot.from_json(str(testbot_json))

        assert testbot.builder_name == "eos_gnu_13.2_serial"
        assert testbot.max_cpus == 4
        assert testbot.with_tdirs == ["v1", "v2"]
        assert testbot.force_mpi is True
        assert testbot.has_mpi is False
        # Defaults for everything not in the JSON payload.
        assert testbot.max_gpus == 0
        assert testbot.without_tdirs == []

    def test_from_json_with_and_without_tdirs_merges_with_warning(self, tmp_path, monkeypatch):
        """Both filters set (e.g. a Force build's with_tdirs landing on a
        builder whose static config already sets without_tdirs) must warn
        and merge -- with_tdirs wins, minus anything without_tdirs excludes --
        rather than raise.
        """
        self._mock_environment(monkeypatch, defined_cppvars=[])
        testbot_json = tmp_path / "testbot.json"
        testbot_json.write_text(json.dumps({
            "builder_name": "b",
            "max_cpus": 2,
            "has_mpi": False,
            "with_tdirs": ["v1", "v2"],
            "without_tdirs": ["v2", "hpc_gpu_omp"],
        }))

        with pytest.warns(UserWarning, match="mutually exclusive"):
            testbot = TestBot.from_json(str(testbot_json))

        assert testbot.with_tdirs == ["v1"]
        assert testbot.without_tdirs == []

    def test_from_json_with_and_without_tdirs_fully_overlapping_runs_everything(self, tmp_path, monkeypatch):
        """If without_tdirs excludes every entry in with_tdirs, the merged
        with_tdirs is empty -- a second warning, and both fields end up
        empty so run() falls back to its own "no filter" default (run all
        suites) instead of a with_tdirs=[] that would (per the plain
        `if self.with_tdirs:` check in run()) look identical to "not set".
        """
        self._mock_environment(monkeypatch, defined_cppvars=[])
        testbot_json = tmp_path / "testbot.json"
        testbot_json.write_text(json.dumps({
            "builder_name": "b",
            "max_cpus": 2,
            "has_mpi": False,
            "with_tdirs": ["v1", "v2"],
            "without_tdirs": ["v1", "v2"],
        }))

        with pytest.warns(UserWarning, match="falling back to running all test suites"):
            testbot = TestBot.from_json(str(testbot_json))

        assert testbot.with_tdirs == []
        assert testbot.without_tdirs == []

    def test_from_json_missing_mandatory_key_raises(self, tmp_path):
        """A testbot.json missing builder_name/max_cpus must raise a clear ValueError."""
        testbot_json = tmp_path / "testbot.json"
        testbot_json.write_text(json.dumps({"builder_name": "only_this"}))

        with pytest.raises(ValueError, match="Mandatory option max_cpus is not declared"):
            TestBot.from_json(str(testbot_json))

    def test_from_json_ignores_unknown_keys(self, tmp_path, monkeypatch):
        """Extra keys in the JSON payload must be silently ignored, not raise."""
        self._mock_environment(monkeypatch, defined_cppvars=["HAVE_MPI"])
        testbot_json = tmp_path / "testbot.json"
        testbot_json.write_text(json.dumps({
            "builder_name": "b",
            "max_cpus": 2,
            "some_future_field": "unused",
        }))

        testbot = TestBot.from_json(str(testbot_json))
        assert testbot.builder_name == "b"
        assert not hasattr(testbot, "some_future_field")

    def test_from_json_gpu_mismatch_raises(self, tmp_path, monkeypatch):
        """max_gpus > 0 on a build without HAVE_GPU must clamp to 0 (with a warning)."""
        self._mock_environment(monkeypatch, defined_cppvars=[])
        testbot_json = tmp_path / "testbot.json"
        testbot_json.write_text(json.dumps({"builder_name": "b", "max_cpus": 2, "max_gpus": 2, "has_mpi": False}))

        with pytest.warns(UserWarning, match="not compiled with GPU support"):
            testbot = TestBot.from_json(str(testbot_json))
        assert testbot.max_gpus == 0

    def test_from_json_mpi_true_but_build_lacks_it_raises(self, tmp_path, monkeypatch):
        """has_mpi=True (default) but build lacks HAVE_MPI must raise ValueError."""
        self._mock_environment(monkeypatch, defined_cppvars=[])
        testbot_json = tmp_path / "testbot.json"
        testbot_json.write_text(json.dumps({"builder_name": "b", "max_cpus": 2}))
        with pytest.raises(ValueError, match="declares has_mpi=True"):
            TestBot.from_json(str(testbot_json))

    def test_from_json_defaults_mpirun_np_when_builder_omits_it(self, tmp_path, monkeypatch):
        """An MPI builder with no mpirun_np must get a working launcher, with a warning.

        Regression test for alps_gnu_14.2_cov: mpirun_np is optional in
        builders.yaml, and the empty value reached JobRunner, which still took
        its MPI branch and built a command line starting with the process count
        ("/bin/timeout: failed to run command '2'", retcode 127 on every np > 1
        test).
        """
        self._mock_environment(monkeypatch, defined_cppvars=["HAVE_MPI"])
        testbot_json = tmp_path / "testbot.json"
        testbot_json.write_text(json.dumps({"builder_name": "alps_gnu_14.2_cov", "max_cpus": 8}))

        with pytest.warns(UserWarning, match="no mpirun_np"):
            testbot = TestBot.from_json(str(testbot_json))

        assert testbot.mpirun_np
        assert testbot.mpirun_np.endswith(" -n")

    def test_from_json_keeps_an_explicit_mpirun_np(self, tmp_path, monkeypatch):
        """A builder that pins mpirun_np must keep it verbatim, with no warning."""
        self._mock_environment(monkeypatch, defined_cppvars=["HAVE_MPI"])
        testbot_json = tmp_path / "testbot.json"
        testbot_json.write_text(json.dumps({
            "builder_name": "b", "max_cpus": 2, "mpirun_np": "/opt/mpich/bin/mpiexec -n",
        }))

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            testbot = TestBot.from_json(str(testbot_json))

        assert testbot.mpirun_np == "/opt/mpich/bin/mpiexec -n"
        # Other warnings (e.g. a missing timeout binary) are unrelated and fine.
        assert not [w for w in caught if "mpirun_np" in str(w.message)]

    def test_from_json_leaves_mpirun_np_alone_for_serial_builders(self, tmp_path, monkeypatch):
        """has_mpi=False builders never launch under MPI, so nothing is filled in."""
        self._mock_environment(monkeypatch, defined_cppvars=[])
        testbot_json = tmp_path / "testbot.json"
        testbot_json.write_text(json.dumps({"builder_name": "b", "max_cpus": 2, "has_mpi": False}))

        testbot = TestBot.from_json(str(testbot_json))
        assert testbot.mpirun_np == ""

    def test_from_json_mpi_false_but_build_has_it_raises(self, tmp_path, monkeypatch):
        """has_mpi=False but build has HAVE_MPI must raise ValueError."""
        self._mock_environment(monkeypatch, defined_cppvars=["HAVE_MPI"])
        testbot_json = tmp_path / "testbot.json"
        testbot_json.write_text(json.dumps({"builder_name": "b", "max_cpus": 2, "has_mpi": False}))
        with pytest.raises(ValueError, match="declares has_mpi=False"):
            TestBot.from_json(str(testbot_json))


class TestWorkdirReplacement:
    """Tests for the opt-in stale work-directory removal behavior."""

    @staticmethod
    def make_bot(remove_existing: bool) -> TestBot:
        bot = object.__new__(TestBot)
        bot.remove_existing_workdirs = remove_existing
        return bot

    def test_default_still_refuses_existing_directory(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "TestBot_MPI1").mkdir()
        with pytest.raises(RuntimeError, match="already exists"):
            self.make_bot(False)._create_workdir("TestBot_MPI1")

    def test_rf_removes_and_recreates_existing_directory(self, tmp_path, monkeypatch, capsys):
        monkeypatch.chdir(tmp_path)
        workdir = tmp_path / "TestBot_MPI1"
        workdir.mkdir()
        (workdir / "stale.txt").write_text("old")

        self.make_bot(True)._create_workdir("TestBot_MPI1")

        assert workdir.is_dir()
        assert not (workdir / "stale.txt").exists()
        assert "Removing existing TestBot work directory" in capsys.readouterr().out

    @pytest.mark.parametrize("existing_kind", ["file", "symlink"])
    def test_rf_refuses_files_and_symlinks(self, tmp_path, monkeypatch, existing_kind):
        monkeypatch.chdir(tmp_path)
        target = tmp_path / "TestBot_MPI1"
        if existing_kind == "file":
            target.write_text("not a directory")
        else:
            real_dir = tmp_path / "real"
            real_dir.mkdir()
            target.symlink_to(real_dir, target_is_directory=True)

        with pytest.raises(RuntimeError, match="Refusing to remove"):
            self.make_bot(True)._create_workdir("TestBot_MPI1")


class TestBuildParser:
    """Tests for the run/analyze/print argparse CLI."""

    def test_run_defaults_and_explicit_json(self):
        parser = build_parser()
        ns = parser.parse_args(["run"])
        assert ns.command == "run"
        assert ns.testbot_json is None

        ns = parser.parse_args(["run", "my.json"])
        assert ns.testbot_json == "my.json"
        assert ns.remove_existing_workdirs is False

        ns = parser.parse_args(["run", "-rf", "my.json"])
        assert ns.remove_existing_workdirs is True

    def test_analyze_defaults_and_explicit_tag(self):
        parser = build_parser()
        ns = parser.parse_args(["analyze"])
        assert ns.command == "analyze"
        assert ns.tag == "unknown"

        ns = parser.parse_args(["analyze", "v9.2.0"])
        assert ns.tag == "v9.2.0"

    def test_benchmark_options(self):
        parser = build_parser()
        ns = parser.parse_args(["benchmark", "my.json", "--py-nprocs", "1", "4", "8", "--profile"])
        assert ns.command == "benchmark"
        assert ns.testbot_json == "my.json"
        assert ns.py_nprocs == [1, 4, 8]
        assert ns.profile is True
        assert ns.remove_existing_workdirs is False

        ns = parser.parse_args(
            ["benchmark", "my.json", "--py-nprocs", "1", "4", "-rf"]
        )
        assert ns.remove_existing_workdirs is True

    def test_print_defaults_and_explicit_json(self):
        parser = build_parser()
        ns = parser.parse_args(["print"])
        assert ns.command == "print"
        assert ns.testbot_json is None

        ns = parser.parse_args(["print", "my.json"])
        assert ns.testbot_json == "my.json"

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
        monkeypatch.setattr(sys, "argv", ["testbot.py", "run", "my.json"])
        mock_instance = MagicMock()
        mock_instance.run.return_value = 3
        with patch("tests.testbot.TestBot.from_json", return_value=mock_instance) as mock_from_json:
            assert main() == 3
        mock_from_json.assert_called_once_with("my.json")
        mock_instance.run.assert_called_once()
        assert mock_instance.remove_existing_workdirs is False

    def test_main_propagates_remove_existing_workdirs(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["testbot.py", "run", "-rf", "my.json"])
        mock_instance = MagicMock()
        mock_instance.run.return_value = 0
        with patch("tests.testbot.TestBot.from_json", return_value=mock_instance):
            assert main() == 0
        assert mock_instance.remove_existing_workdirs is True

    def test_main_dispatches_to_print_without_running(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["testbot.py", "print", "my.json"])
        mock_instance = MagicMock()
        mock_instance.__str__ = MagicMock(return_value="fake-config-dump")
        with patch("tests.testbot.TestBot.from_json", return_value=mock_instance) as mock_from_json:
            assert main() == 0
        mock_from_json.assert_called_once_with("my.json")
        mock_instance.run.assert_not_called()
        assert "fake-config-dump" in capsys.readouterr().out

    def test_main_run_with_no_json_uses_default(self, monkeypatch):
        """A bare `testbot.py run` (the analysis.sh invocation) passes path=None."""
        monkeypatch.setattr(sys, "argv", ["testbot.py", "run"])
        mock_instance = MagicMock()
        mock_instance.run.return_value = 0
        with patch("tests.testbot.TestBot.from_json", return_value=mock_instance) as mock_from_json:
            assert main() == 0
        mock_from_json.assert_called_once_with(None)

    def test_main_dispatches_to_benchmark(self, monkeypatch):
        monkeypatch.setattr(
            sys, "argv", ["testbot.py", "benchmark", "my.json", "--py-nprocs", "1", "4"]
        )
        with patch("tests.testbot.benchmark", return_value=2) as mock_benchmark:
            assert main() == 2
        mock_benchmark.assert_called_once_with("my.json", [1, 4], False)

    def test_main_rejects_invalid_benchmark_cpu_values(self, monkeypatch):
        monkeypatch.setattr(
            sys, "argv", ["testbot.py", "benchmark", "my.json", "--py-nprocs", "0"]
        )
        with pytest.raises(SystemExit, match="positive integers"):
            main()


def test_benchmark_writes_summary(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    bots = []

    def make_bot(_path):
        bot = MagicMock()
        bot.run.return_value = 0
        bot.run_summaries = []
        bots.append(bot)
        return bot

    monkeypatch.setattr("tests.testbot.TestBot.from_json", make_bot)
    times = iter([10.0, 14.0, 20.0, 28.0])
    monkeypatch.setattr("tests.testbot.time.perf_counter", lambda: next(times))

    assert benchmark("testbot.json", [1, 4]) == 0

    data = json.loads((tmp_path / "testbot_benchmark.json").read_text())
    assert data == [
        {
            "py_nprocs": 1,
            "wall_time": 4.0,
            "speedup": 1.0,
            "parallel_efficiency": 1.0,
            "nexecuted": 0,
            "tests_per_second": 0.0,
            "returncode": 0,
        },
        {
            "py_nprocs": 4,
            "wall_time": 8.0,
            "speedup": 0.5,
            "parallel_efficiency": 0.125,
            "nexecuted": 0,
            "tests_per_second": 0.0,
            "returncode": 0,
        },
    ]
    assert bots[0].py_nprocs_override == 1
    assert bots[1].py_nprocs_override == 4
    assert bots[0].workdir_prefix == "Benchmark_PY1_"
    assert bots[1].workdir_prefix == "Benchmark_PY4_"
    output = capsys.readouterr().out
    assert "Benchmark summary" in output
    assert "efficiency" in output
    assert "12.5%" in output


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
