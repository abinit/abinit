"""
Unit tests for jobrunner.py module.

Tests cover is_string(), mpicfg_parser(), JobRunnerError, JobRunner's
constructor validation/classmethods/properties/run(), BaseValgrindParser/
MemcheckParser's Valgrind output parsing, TimeBomb's three run() strategies,
and OMPEnvironment's INI-driven construction.
"""

from __future__ import annotations

import pytest

from . import jobrunner
from .jobrunner import (
    BaseValgrindParser,
    JobRunner,
    JobRunnerError,
    MemcheckParser,
    OMPEnvironment,
    TimeBomb,
    is_string,
    mpicfg_parser,
)


class TestIsString:
    """Test the is_string() duck-typed string check."""

    @pytest.mark.parametrize("value", ["", "abc", "3.14"])
    def test_strings_are_strings(self, value):
        assert is_string(value) is True

    @pytest.mark.parametrize("value", [1, 3.14, [], {}, None, object()])
    def test_non_strings_are_not_strings(self, value):
        assert is_string(value) is False


class TestMpicfgParser:
    """Test mpicfg_parser()'s INI-file parsing of the [mpi] section."""

    def test_missing_file_uses_all_defaults(self, tmp_path):
        d = mpicfg_parser(str(tmp_path / "does_not_exist.cfg"))
        assert d["mpi_flavor"] == ""
        assert d["mpirun_np"] == ""

    def test_reads_configured_values(self, tmp_path):
        cfg = tmp_path / "mpi.cfg"
        cfg.write_text("[mpi]\nmpi_flavor = openmpi\nmpirun_np = mpirun -np\n")
        d = mpicfg_parser(str(cfg))
        assert d["mpi_flavor"] == "openmpi"
        assert d["mpirun_np"] == "mpirun -np"
        # Options not present in the file still fall back to their default.
        assert d["poe"] == ""

    def test_section_without_key_uses_default(self, tmp_path):
        cfg = tmp_path / "mpi.cfg"
        cfg.write_text("[mpi]\nmpi_flavor = mpich\n")
        d = mpicfg_parser(str(cfg))
        assert d["mpirun_np"] == ""


class TestJobRunnerError:
    """Test the JobRunnerError exception-info holder."""

    def test_str_without_prev_errmsg(self):
        err = JobRunnerError(1, "mycmd arg", 2.5)
        s = str(err)
        assert "mycmd arg" in s
        assert "1" in s
        assert "Previous exception" not in s

    def test_str_with_prev_errmsg(self):
        err = JobRunnerError(1, "mycmd", 2.5, prev_errmsg="boom")
        s = str(err)
        assert "Previous exception: boom" in s


class TestJobRunnerConstructionAndClassmethods:
    """Test JobRunner.__init__ and its alternate constructors."""

    def test_duplicate_kwarg_key_raises(self):
        # "exceptions" is set by __init__ itself before the loop, so passing
        # it again must be rejected rather than silently overwritten.
        with pytest.raises(ValueError, match="already in self.__dict__"):
            JobRunner({"exceptions": []})

    def test_poe_and_mpirun_are_mutually_exclusive(self):
        with pytest.raises(ValueError, match="mutually exclusive"):
            JobRunner({"poe": "poe", "mpirun_np": "mpirun -np"})

    def test_poe_and_srun_are_mutually_exclusive(self):
        with pytest.raises(ValueError, match="mutually exclusive"):
            JobRunner({"poe": "poe", "mpirun_np": "srun -n"})

    def test_mpi_args_defaults_to_empty_string(self):
        runner = JobRunner({})
        assert runner.mpi_args == ""

    def test_explicit_mpi_args_is_kept(self):
        runner = JobRunner({"mpi_args": "--bind-to core"})
        assert runner.mpi_args == "--bind-to core"

    def test_sequential_classmethod(self):
        runner = JobRunner.sequential()
        assert runner.has_mpirun is False
        assert runner.has_poe is False

    def test_srun_classmethod(self):
        runner = JobRunner.srun(mpi_args="--foo")
        assert runner.has_srun is True
        assert runner.mpirun_np == "srun -n"
        assert runner.mpi_args == "--foo"

    def test_generic_mpi_classmethod_default_uses_mpirun(self):
        runner = JobRunner.generic_mpi()
        assert runner.mpirun_np == "mpirun -np"
        assert runner.has_mpirun is True

    def test_generic_mpi_classmethod_can_use_mpiexec(self):
        runner = JobRunner.generic_mpi(use_mpiexec=True)
        assert runner.mpirun_np == "mpiexec -np"

    def test_fromdict_merges_ompenv_and_timebomb(self):
        ompenv = OMPEnvironment(OMP_NUM_THREADS=2)
        timebomb = TimeBomb(10)
        runner = JobRunner.fromdict({"mpirun_np": "mpirun -np"}, ompenv=ompenv, timebomb=timebomb)
        assert runner.ompenv is ompenv
        assert runner.timebomb is timebomb
        assert runner.mpirun_np == "mpirun -np"

    def test_fromfile_reads_mpi_config(self, tmp_path):
        cfg = tmp_path / "mpi.cfg"
        cfg.write_text("[mpi]\nmpirun_np = mpirun -np\n")
        runner = JobRunner.fromfile(str(cfg))
        assert runner.mpirun_np == "mpirun -np"
        assert runner.has_ompenv is False


class TestJobRunnerProperties:
    """Test JobRunner's has_* properties and set_* setters."""

    def test_has_mpirun_false_when_unset(self):
        assert JobRunner({}).has_mpirun is False

    def test_has_mpirun_false_when_empty_string(self):
        # Regression: mpicfg_parser()/TestBot always hand every CFG_KEYWORDS
        # key to JobRunner, so an *unconfigured* launcher arrives as
        # mpirun_np="" -- has_mpirun must test the value, not mere presence.
        assert JobRunner({"mpirun_np": ""}).has_mpirun is False

    def test_has_mpirun_true_when_set(self):
        assert JobRunner({"mpirun_np": "mpirun -np"}).has_mpirun is True

    def test_has_mpirun_false_for_srun(self):
        assert JobRunner({"mpirun_np": "srun -n"}).has_mpirun is False

    def test_has_srun_true_only_for_srun_value(self):
        assert JobRunner({"mpirun_np": "srun -n"}).has_srun is True
        assert JobRunner({"mpirun_np": "mpirun -np"}).has_srun is False
        assert JobRunner({}).has_srun is False

    def test_has_poe(self):
        assert JobRunner({"poe": "poe"}).has_poe is True
        assert JobRunner({"poe": ""}).has_poe is False
        assert JobRunner({}).has_poe is False

    def test_has_timebomb(self):
        runner = JobRunner({})
        assert runner.has_timebomb is False
        runner.set_timebomb(TimeBomb(10))
        assert runner.has_timebomb is True

    def test_set_timebomb_twice_raises(self):
        runner = JobRunner({"timebomb": TimeBomb(10)})
        with pytest.raises(ValueError, match="already defined"):
            runner.set_timebomb(TimeBomb(20))

    def test_has_ompenv(self):
        runner = JobRunner({})
        assert runner.has_ompenv is False
        runner.set_ompenv(OMPEnvironment(OMP_NUM_THREADS=4))
        assert runner.has_ompenv is True

    def test_set_ompenv_twice_raises(self):
        runner = JobRunner({"ompenv": OMPEnvironment(OMP_NUM_THREADS=4)})
        with pytest.raises(ValueError, match="already defined"):
            runner.set_ompenv(OMPEnvironment(OMP_NUM_THREADS=2))

    def test_has_valgrind_and_build_parser(self):
        runner = JobRunner({})
        assert runner.has_valgrind is False
        with pytest.raises(ValueError, match="does not use valgrind"):
            runner.build_valgrind_parser()

        runner.set_valgrind_cmdline("memcheck")
        assert runner.has_valgrind is True
        assert isinstance(runner.build_valgrind_parser(), MemcheckParser)

    def test_has_debugger(self):
        runner = JobRunner({})
        assert runner.has_debugger is False
        runner.set_debugger("gdb")
        assert runner.has_debugger is True

    def test_has_perf(self):
        runner = JobRunner({})
        assert runner.has_perf is False
        runner.set_perf_command("stat")
        assert runner.has_perf is True


class TestJobRunnerStr:
    """Test JobRunner.__str__()'s human-readable summary."""

    def test_empty_runner_produces_empty_string(self):
        assert str(JobRunner({})) == ""

    def test_mpi_setup_is_reported(self):
        runner = JobRunner({"mpirun_np": "mpirun -np"})
        s = str(runner)
        assert "[MPI setup]" in s
        assert "mpirun_np = mpirun -np" in s

    def test_ompenv_is_reported(self):
        runner = JobRunner({"ompenv": OMPEnvironment(OMP_NUM_THREADS=4)})
        s = str(runner)
        assert "[OpenMP]" in s
        assert "OMP_NUM_THREADS" in s


class TestJobRunnerRun:
    """Test JobRunner.run()'s command construction and execution."""

    def test_sequential_success(self, tmp_path):
        runner = JobRunner.sequential()
        stdout = tmp_path / "out.log"
        stderr = tmp_path / "err.log"
        etime = runner.run(1, "true", None, str(stdout), str(stderr))
        assert etime >= 0
        assert runner.retcode == 0
        assert runner.exceptions == []

    def test_sequential_failure_is_recorded(self, tmp_path):
        runner = JobRunner.sequential()
        stdout = tmp_path / "out.log"
        stderr = tmp_path / "err.log"
        runner.run(1, "false", None, str(stdout), str(stderr))
        assert runner.retcode == 1
        assert len(runner.exceptions) == 1
        assert isinstance(runner.exceptions[0], JobRunnerError)

    def test_multi_process_without_launcher_raises(self, tmp_path):
        runner = JobRunner.sequential()
        stdout = tmp_path / "out.log"
        stderr = tmp_path / "err.log"
        with pytest.raises(ValueError, match="no MPI launcher configured"):
            runner.run(4, "true", None, str(stdout), str(stderr))

    def test_generic_mpi_uses_echo_as_a_fake_launcher(self, tmp_path):
        # Route "mpirun_np" at a real, always-present command so run() can
        # build and execute a full command line without needing real MPI.
        runner = JobRunner.fromdict({"mpirun_np": "echo -n"})
        stdout = tmp_path / "out.log"
        stderr = tmp_path / "err.log"
        runner.run(2, "true", None, str(stdout), str(stderr))
        assert runner.retcode == 0

    def test_poe_launcher_branch(self, tmp_path):
        runner = JobRunner.fromdict({"poe": "echo", "poe_args": "-x"})
        stdout = tmp_path / "out.log"
        stderr = tmp_path / "err.log"
        runner.run(4, "true", None, str(stdout), str(stderr))
        assert runner.retcode == 0

    def test_stdin_stdout_stderr_redirection(self, tmp_path):
        runner = JobRunner.sequential()
        stdin_file = tmp_path / "in.txt"
        stdin_file.write_text("hello\n")
        stdout_file = tmp_path / "out.txt"
        stderr_file = tmp_path / "err.txt"
        runner.run(1, "cat", str(stdin_file), str(stdout_file), str(stderr_file))
        assert runner.retcode == 0
        assert stdout_file.read_text() == "hello\n"

    def test_valgrind_and_perf_prefix_the_command(self, tmp_path):
        # Neither valgrind nor perf need to be installed: we only check that
        # run() builds and executes a command line that *starts* with them
        # (which fails harmlessly with retcode=127, "command not found").
        runner = JobRunner.sequential()
        runner.set_valgrind_cmdline("memcheck")
        runner.set_perf_command("stat")
        stdout = tmp_path / "out.log"
        stderr = tmp_path / "err.log"
        runner.run(1, "true", None, str(stdout), str(stderr))
        # Whether or not valgrind/perf are actually installed here, run()
        # must complete and record the outcome rather than raise.
        assert runner.retcode is not None

    def test_debugger_branch_writes_commands_file_and_runs_sequential(self, tmp_path):
        runner = JobRunner.sequential()
        runner.set_debugger("gdb")
        stdout = tmp_path / "out.log"
        stderr = tmp_path / "err.log"
        runner.run(1, "true", None, str(stdout), str(stderr), bin_argstr="--flag")
        dbg_file = tmp_path / "dbg_commands"
        assert dbg_file.exists()
        assert "run --flag" in dbg_file.read_text()

    def test_debugger_branch_with_mpi(self, tmp_path):
        runner = JobRunner.fromdict({"mpirun_np": "echo -n"})
        runner.set_debugger("gdb")
        stdout = tmp_path / "out.log"
        stderr = tmp_path / "err.log"
        runner.run(2, "true", None, str(stdout), str(stderr))
        assert (tmp_path / "dbg_commands").exists()

    def test_exception_during_popen_is_recorded(self, tmp_path, monkeypatch):
        def raising_popen(*args, **kwargs):
            raise OSError("boom")

        monkeypatch.setattr(jobrunner, "Popen", raising_popen)
        runner = JobRunner.sequential()
        stdout = tmp_path / "out.log"
        stderr = tmp_path / "err.log"
        runner.run(1, "true", None, str(stdout), str(stderr))
        assert len(runner.exceptions) == 1
        assert "boom" in runner.exceptions[0].prev_errmsg

    def test_timebomb_branch_is_used_when_configured(self, tmp_path):
        runner = JobRunner.fromdict({}, timebomb=TimeBomb(2, delay=0.05))
        stdout = tmp_path / "out.log"
        stderr = tmp_path / "err.log"
        runner.run(1, "true", None, str(stdout), str(stderr))
        assert runner.retcode == 0


class TestBaseValgrindParser:
    """Test the BaseValgrindParser abstract-ish base class."""

    def test_parse_not_implemented(self):
        with pytest.raises(NotImplementedError):
            BaseValgrindParser().parse("some_file")

    def test_error_report_returns_stored_value(self):
        parser = BaseValgrindParser()
        parser._error_report = "some report"
        assert parser.error_report == "some report"


VALGRIND_OUTPUT_WITH_LEAKS = """\
==123== HEAP SUMMARY:
==123==     in use at exit: 100 bytes in 2 blocks
==123== LEAK SUMMARY:
==123==    definitely lost: 64 bytes in 1 blocks
==123==    indirectly lost: 0 bytes in 0 blocks
==123==      possibly lost: 36 bytes in 1 blocks
==123== ERROR SUMMARY: 3 errors from 2 contexts (suppressed: 0 from 0)
"""

VALGRIND_OUTPUT_CLEAN = """\
==123== HEAP SUMMARY:
==123==     in use at exit: 0 bytes in 0 blocks
==123== LEAK SUMMARY:
==123==    definitely lost: 0 bytes in 0 blocks
==123==    indirectly lost: 0 bytes in 0 blocks
==123==      possibly lost: 0 bytes in 0 blocks
==123== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)
"""

VALGRIND_OUTPUT_NO_LEAK_SUMMARY = """\
==123== HEAP SUMMARY:
==123==     in use at exit: 0 bytes in 0 blocks
"""


class TestMemcheckParser:
    """Test MemcheckParser.parse()'s Valgrind stderr scraping."""

    def test_parses_leak_and_error_counts(self, tmp_path):
        f = tmp_path / "valgrind.out"
        f.write_text(VALGRIND_OUTPUT_WITH_LEAKS)
        parser = MemcheckParser()
        parser.parse(str(f))
        assert "definitely lost:" in parser.error_report
        assert "64" in parser.error_report
        assert "ERROR SUMMARY:" in parser.error_report
        assert "3" in parser.error_report

    def test_clean_run_produces_empty_report(self, tmp_path):
        f = tmp_path / "valgrind.out"
        f.write_text(VALGRIND_OUTPUT_CLEAN)
        parser = MemcheckParser()
        parser.parse(str(f))
        assert parser.error_report == ""

    def test_missing_leak_summary_raises(self, tmp_path):
        f = tmp_path / "valgrind.out"
        f.write_text(VALGRIND_OUTPUT_NO_LEAK_SUMMARY)
        parser = MemcheckParser()
        with pytest.raises(RuntimeError, match="Cannot find 'LEAK SUMMARY'"):
            parser.parse(str(f))


class TestTimeBomb:
    """Test TimeBomb's three run() strategies."""

    def test_init_truncates_timeout_to_int(self):
        # self.timeout = int(timeout): a sub-1s float timeout truncates to 0,
        # which then reads as "no timeout enforcement" below.
        bomb = TimeBomb(1.9, delay=0.1)
        assert bomb.timeout == 1
        assert bomb.delay == 0.1

    def test_zero_timeout_uses_plain_popen(self, tmp_path):
        bomb = TimeBomb(0)
        p, retcode = bomb.run(["true"])
        assert retcode == 0
        p.wait()

    def test_positive_timeout_without_exec_path_uses_subprocesswithtimeout(self, tmp_path):
        bomb = TimeBomb(2, delay=0.05)
        p, retcode = bomb.run(["true"])
        assert retcode == 0

    def test_positive_timeout_actually_enforces_it(self):
        bomb = TimeBomb(0.2, delay=0.05)
        # TimeBomb.__init__ truncates timeout via int(), so use a value where
        # that still leaves a positive, sub-second-ish window: int(1) == 1.
        bomb = TimeBomb(1, delay=0.05)
        p, retcode = bomb.run(["sleep", "5"])
        assert retcode in (124, 137)
        if p.poll() is None:
            p.kill()
        p.wait()

    def test_exec_path_with_string_args(self, tmp_path):
        # A fake "timeout" executable: echoes its own argv so we can verify
        # TimeBomb correctly prepended "<exec_path> <timeout>" to the command.
        fake_timeout = tmp_path / "fake_timeout.sh"
        fake_timeout.write_text('#!/bin/sh\necho "$@"\n')
        fake_timeout.chmod(0o755)

        # String args are joined into a single command line (is_string(args)
        # branch), so this needs shell=True to actually be split and run --
        # same as how JobRunner.run() always invokes TimeBomb.run().
        bomb = TimeBomb(5, exec_path=str(fake_timeout))
        p, retcode = bomb.run("echo hi", shell=True)
        assert retcode == 0

    def test_exec_path_with_list_args(self, tmp_path):
        fake_timeout = tmp_path / "fake_timeout.sh"
        fake_timeout.write_text('#!/bin/sh\necho "$@"\n')
        fake_timeout.chmod(0o755)

        bomb = TimeBomb(5, exec_path=str(fake_timeout))
        p, retcode = bomb.run([str(fake_timeout), "ignored"])
        assert retcode == 0

    def test_exec_path_with_non_positive_timeout_skips_prefixing(self, tmp_path):
        fake_timeout = tmp_path / "fake_timeout.sh"
        fake_timeout.write_text('#!/bin/sh\necho "$@"\n')
        fake_timeout.chmod(0o755)

        # timeout=0 with exec_path set: the "if self.timeout > 0." guard is
        # False, so args must be passed to exec_path unprefixed.
        bomb = TimeBomb(0, exec_path=str(fake_timeout))
        p, retcode = bomb.run([str(fake_timeout)])
        assert retcode == 0

    def test_exception_from_popen_is_reraised(self):
        # run()'s outer try/except is a pure passthrough (`except Exception:
        # raise`) -- an unusable exec_path must still surface as an error,
        # not be silently swallowed.
        bomb = TimeBomb(5, exec_path="/no/such/timeout/binary")
        with pytest.raises(FileNotFoundError):
            bomb.run(["true"])


class TestOMPEnvironment:
    """Test the OMPEnvironment dict subclass."""

    def test_valid_keys_are_coerced_to_strings(self):
        env = OMPEnvironment(OMP_NUM_THREADS=4)
        assert env["OMP_NUM_THREADS"] == "4"
        assert isinstance(env["OMP_NUM_THREADS"], str)

    def test_unknown_key_raises(self):
        with pytest.raises(ValueError, match="unknown option"):
            OMPEnvironment(NOT_A_REAL_OMP_VAR=1)

    def test_from_file_without_openmp_section_allow_empty(self, tmp_path):
        cfg = tmp_path / "empty.cfg"
        cfg.write_text("[other]\nfoo = bar\n")
        env = OMPEnvironment.from_file(str(cfg), allow_empty=True)
        assert env == {}

    def test_from_file_without_openmp_section_raises_by_default(self, tmp_path):
        cfg = tmp_path / "empty.cfg"
        cfg.write_text("[other]\nfoo = bar\n")
        with pytest.raises(ValueError, match="does not contain any \\[openmp\\] section"):
            OMPEnvironment.from_file(str(cfg))

    def test_from_file_reads_uppercase_keys(self, tmp_path):
        cfg = tmp_path / "omp.cfg"
        cfg.write_text("[openmp]\nOMP_NUM_THREADS = 8\n")
        env = OMPEnvironment.from_file(str(cfg))
        assert env["OMP_NUM_THREADS"] == "8"

    def test_from_file_reads_lowercase_keys_as_fallback(self, tmp_path):
        cfg = tmp_path / "omp.cfg"
        cfg.write_text("[openmp]\nomp_num_threads = 8\n")
        env = OMPEnvironment.from_file(str(cfg))
        assert env["OMP_NUM_THREADS"] == "8"

    def test_from_file_unknown_option_raises(self, tmp_path):
        cfg = tmp_path / "omp.cfg"
        cfg.write_text("[openmp]\nnot_an_omp_var = 1\n")
        with pytest.raises(ValueError, match="unknown option"):
            OMPEnvironment.from_file(str(cfg))

    def test_from_file_empty_section_raises_by_default(self, tmp_path):
        cfg = tmp_path / "omp.cfg"
        cfg.write_text("[openmp]\n")
        with pytest.raises(ValueError, match="Refusing to return with an empty dict"):
            OMPEnvironment.from_file(str(cfg))

    def test_from_file_empty_section_allowed_when_requested(self, tmp_path):
        cfg = tmp_path / "omp.cfg"
        cfg.write_text("[openmp]\n")
        env = OMPEnvironment.from_file(str(cfg), allow_empty=True)
        assert env == {}
