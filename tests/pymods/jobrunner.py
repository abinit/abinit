"""
Job execution and management tools for ABINIT tests.

This module provides a unified interface for running ABINIT binaries in various
environments, including sequential, MPI (OpenMPI, MPICH, srun, poe), and OpenMP.
It also includes facilities for enforcing timeouts, profiling with Valgrind or
perf, and managing environment variables.
"""
from __future__ import annotations

import os
import shlex
import sys
import time
from collections.abc import Callable
from configparser import NoOptionError
from subprocess import Popen
from typing import Any, cast

from .subprocesswithtimeout import SubProcessWithTimeout

try:
    from configparser import SafeConfigParser
except ImportError:
    from configparser import ConfigParser as SafeConfigParser

import logging

logger = logging.getLogger(__name__)

__version__ = "0.1"
__author__ = "Matteo Giantomassi"

__all__ = [
    "JobRunner",
    "OMPEnvironment",
    "TimeBomb",
]


CFG_KEYWORDS: dict[str, tuple[type, str, str, str]] = {
# keyword             : (parser, default value i.e. NO MPI, section, description)
# [mpi]
"mpi_flavor"         : (str, "", "mpi", "Name of the MPI implementation e.g. openmpi, mpich2 ..."),
"mpi_version"        : (str, "", "mpi", "Version of the MPI implementation e.g. 1.2.4"),
"mpi_prefix"         : (str, "", "mpi", "Top directory of the MPI library. e.g. /shared/openmpi-ifc10"),
"mpirun_np"          : (str, "", "mpi", "String specifying how to execute a binary with N processors"),
#"mpirun_extra_args"  : (str, "", "mpi", "Options passed after the 'mpirun -np 3' command"),
#"np_option"         : (str, "", "-np", "")
"poe"                : (str, "", "mpi", "poe location"),
"poe_args"           : (str, "", "mpi", "arguments passed to poe"),
#"info"              : (str, "", "mpi", "String with optional information"),
}


def is_string(s: Any) -> bool:
    """
    Check if the input is a string-like object.

    Uses duck typing (concatenation with a str succeeds) rather than
    `isinstance(s, str)`, so any object supporting `+` with a string
    (e.g. a str subclass) is accepted too.

    Args:
        s: The object to check.

    Returns:
        bool: True if s is a string, False otherwise.
    """
    try:
        s + "hello"
        return True
    except TypeError:
        return False


def mpicfg_parser(fname: str, defaults: dict[str, str] | None = None) -> dict[str, str]:
    """
    Parse a configuration file (INI format) for MPI options.

    Args:
        fname: Path to the configuration file.
        defaults: Default values for the parser.

    Returns:
        dict: A dictionary containing the parsed MPI options.

    Raises:
        ValueError: If a keyword's value cannot be converted by its declared parser type.
    """
    logger.debug(f"Parsing [MPI] section in file: {fname}")

    parser = SafeConfigParser(defaults)
    parser.read(fname)

    # Read variables needed to handle the job.
    d = {}

    for key, tup in CFG_KEYWORDS.items():
        line_parser = tup[0]
        section = tup[2]
        if section in parser.sections():
            try:
                d[key] = parser.get(section, key)
            except NoOptionError:
                # Section exists but option is not specified. Use default value.
                d[key] = tup[1]
        else:
            # Section does not exist. Use default value.
            d[key] = tup[1]

        # Process the line
        try:
            d[key] = line_parser(d[key])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Wrong line: key = {key} d[key] = {d[key]}") from exc

    return d


class JobRunnerError:
    """
    Exception-like object to store information about job execution failures.

    Inherits from `object`, not `Exception`, because pickling fails if
    JobRunnerError inherits from Exception. This should not be a serious
    problem in practice since JobRunner exceptions are never caught elsewhere
    -- this object is mainly used to store info about the exception in
    JobRunner.exceptions (see JobRunner.run()).
    """

    def __init__(self, return_code: int, cmd: str, run_etime: float, prev_errmsg: str | None = None) -> None:
        """
        Args:
            return_code: Return code of the subprocess
            cmd: Command executed
            run_etime: Elapsed-time
            prev_errmsg: Previous error message.
        """
        self.return_code = return_code
        self.cmd = cmd
        self.run_etime = run_etime
        self.prev_errmsg = prev_errmsg

    def __str__(self) -> str:
        string = f"Command {self.cmd}\n returned exit_code: {self.return_code}\n"
        if self.prev_errmsg:
            string += f"Previous exception: {self.prev_errmsg}"

        return string


class JobRunner:
    """
    Manages the execution of jobs in an MPI or sequential environment.

    This class provides a unified interface for running binaries with support for
    MPI (via mpirun, srun, or poe), OpenMP, and diagnostic tools like Valgrind or perf.

    It abstracts the complexities of different MPI launchers and environment
    configurations, providing resource management and timeout enforcement.
    """

    @classmethod
    def fromdict(cls, kwargs: dict[str, Any], ompenv: OMPEnvironment | None = None, timebomb: TimeBomb | None = None) -> JobRunner:
        """
        Create a JobRunner instance from a dictionary of options.

        Args:
            kwargs: Dictionary containing runner options.
            ompenv: Optional OMPEnvironment instance.
            timebomb: Optional TimeBomb instance for enforcing timeouts.

        Returns:
            JobRunner: A new instance configured with the provided options.
        """
        d = dict(ompenv=ompenv, timebomb=timebomb)
        d.update(kwargs)

        return cls(d)

    @classmethod
    def fromfile(cls, fname: str, timebomb: TimeBomb | None = None) -> JobRunner:
        """
        Create a JobRunner instance from an INI configuration file.

        Args:
            fname: Path to the configuration file.
            timebomb: Optional TimeBomb instance.

        Returns:
            JobRunner: A new instance configured from the file.
        """
        d: dict[str, Any] = mpicfg_parser(fname)
        d["ompenv"] = OMPEnvironment.from_file(fname, allow_empty=True)
        d["timebomb"] = timebomb

        return cls(d)

    @classmethod
    def sequential(cls, ompenv: OMPEnvironment | None = None, timebomb: TimeBomb | None = None) -> JobRunner:
        """
        Create a JobRunner for sequential (non-MPI) execution.

        Args:
            ompenv: Optional OMPEnvironment.
            timebomb: Optional TimeBomb.

        Returns:
            JobRunner: A sequential runner instance.
        """
        return cls(dict(ompenv=ompenv, timebomb=timebomb))

    @classmethod
    def srun(cls, ompenv: OMPEnvironment | None = None, timebomb: TimeBomb | None = None, mpi_args: str = "") -> JobRunner:
        """
        Create a JobRunner configured for Slurm's `srun`.

        Args:
            ompenv: Optional OMPEnvironment.
            timebomb: Optional TimeBomb.
            mpi_args: Extra arguments for the mpirun command.

        Returns:
            JobRunner: A runner instance configured for Slurm.
        """
        d = dict(ompenv=ompenv, timebomb=timebomb, mpi_args=mpi_args)
        d["mpirun_np"] = "srun -n"
        return cls(d)

    @classmethod
    def generic_mpi(cls, ompenv: OMPEnvironment | None = None, use_mpiexec: bool = False, mpi_args: str = "",
                     timebomb: TimeBomb | None = None) -> JobRunner:
        """
        Create a JobRunner for generic MPI execution (mpirun or mpiexec).

        Args:
            ompenv: Optional OMPEnvironment.
            use_mpiexec: If True, use `mpiexec` instead of `mpirun`.
            mpi_args: Extra arguments for the MPI launcher.
            timebomb: Optional TimeBomb.

        Returns:
            JobRunner: A generic MPI runner instance.
        """
        # It should work, provided that the shell environment is properly defined.
        d = dict(ompenv=ompenv, timebomb=timebomb, mpi_args=mpi_args)

        if use_mpiexec:
            d["mpirun_np"] = "mpiexec -np"
        else:
            d["mpirun_np"] = "mpirun -np"

        return cls(d)

    def __init__(self, dic: dict[str, Any]) -> None:
        """
        Initialize the JobRunner.

        Args:
            dic (dict): Dictionary of configuration options.

        Raises:
            ValueError: If poe and (mpirun or srun) are both specified.
        """
        self.exceptions: list[JobRunnerError] = []

        for k, v in dic.items():
            if k not in self.__dict__:
                self.__dict__[k] = v
            else:
                raise ValueError(f"key {k} is already in self.__dict__, cannot overwrite")

        if "mpi_args" not in dic:
            self.mpi_args = ""

        if self.has_poe and (self.has_mpirun or self.has_srun):
            raise ValueError("poe and (mpirun||srun) are mutually exclusive")

    def __str__(self) -> str:
        """
        Return a string representation of the job runner's configuration.

        Returns:
            str: The configuration summary.
        """
        string = ""
        for key in CFG_KEYWORDS:
            attr = getattr(self, str(key), None)
            if attr: string += f"{key} = {attr}\n"

        if string:
            string = "[MPI setup]\n" + string

        if self.has_ompenv:
            string += "[OpenMP]\n" + str(self.ompenv)

        return string

    def set_timebomb(self, timebomb: TimeBomb) -> None:
        """
        Set the timebomb for the runner.

        Args:
            timebomb (TimeBomb): The TimeBomb instance.

        Raises:
            ValueError: If a timebomb is already defined.
        """
        if self.has_timebomb:
            raise ValueError("timebomb is already defined")
        self.timebomb = timebomb

    def set_valgrind_cmdline(self, cmdline: str) -> None:
        """
        Set the command line options to be passed to VALGRIND.

        Args:
            cmdline (str): The command line options.
        """
        self.valgrind_cmdline = cmdline

    @property
    def has_valgrind(self) -> bool:
        """True if we are running the code with VALGRIND."""
        return hasattr(self, "valgrind_cmdline")

    def build_valgrind_parser(self) -> MemcheckParser:
        """
        Build and return a parser for Valgrind output.

        Returns:
            MemcheckParser: The initialized parser.

        Raises:
            ValueError: If Valgrind is not enabled for this runner.
        """
        if not self.has_valgrind: raise ValueError("Runner does not use valgrind!")
        return MemcheckParser()

    def set_debugger(self, debugger: str) -> None:
        """
        Set the debugger executable.

        Args:
            debugger (str): Path to the debugger executable.
        """
        self.debugger = debugger

    @property
    def has_debugger(self) -> bool:
        """True if we are running the executable under the control of a debugger."""
        return hasattr(self, "debugger")

    def set_perf_command(self, perf_command: str) -> None:
        """
        Set the perf command to be used.

        Args:
            perf_command (str): The perf command string.
        """
        self.perf_command = perf_command

    @property
    def has_perf(self) -> bool:
        """True if we are profiling the run with perf."""
        return hasattr(self, "perf_command")

    @property
    def has_srun(self) -> bool:
        """True if we are running with Slurm srun"""
        return hasattr(self, "mpirun_np") and self.mpirun_np == "srun -n"

    @property
    def has_mpirun(self) -> bool:
        """
        True if we are running a MPI job with mpirun.

        Tests the attribute's *value*, not just its presence: both
        mpicfg_parser and TestBot hand us every CFG_KEYWORDS key, so an
        unconfigured launcher arrives as mpirun_np="". Reporting True for that
        made run() emit [<empty>, nprocs, ..., bin_path, ...], i.e. a command
        line whose first token was the process count ("failed to run command
        '2'", retcode 127) instead of failing outright.
        """
        return bool(getattr(self, "mpirun_np", "")) and self.mpirun_np != "srun -n"

    @property
    def has_poe(self) -> bool:
        """True if are using IBM poe for MPI executions."""
        return hasattr(self, "poe") and bool(self.poe)

    @property
    def has_timebomb(self) -> bool:
        """
        True if we are running the job under the control of
        an application that will enforce a timeout.
        """
        return hasattr(self, "timebomb") and bool(self.timebomb)

    @property
    def has_ompenv(self) -> bool:
        """True if we are using OpenMP."""
        return hasattr(self, "ompenv") and bool(self.ompenv)

    def set_ompenv(self, ompenv: OMPEnvironment) -> None:
        """
        Set the value of the OpenMP environmental variables.

        Args:
            ompenv (OMPEnvironment): The OMP environment instance.

        Raises:
            ValueError: If an OMP environment is already defined.
        """
        if self.has_ompenv:
            raise ValueError("ompenv is already defined")
        self.ompenv = ompenv

    def run(self, mpi_nprocs: int, bin_path: str, stdin_fname: str | None, stdout_fname: str, stderr_fname: str,
            bin_argstr: str = "", cwd: str | None = None) -> float:
        """
        Execute the job.

        Args:
            mpi_nprocs: Number of MPI processes to launch.
            bin_path: Path to the executable.
            stdin_fname: Optional path to the input file.
            stdout_fname: Optional path to the output file.
            stderr_fname: Optional path to the error file.
            bin_argstr: Extra command-line arguments for the binary.
            cwd: Optional working directory for the execution.

        Returns:
            float: Elapsed time of the execution (in seconds).

        Raises:
            ValueError: If mpi_nprocs != 1 but no MPI launcher (mpirun_np/poe)
                is configured on this runner.

        Note:
            Execution failures (non-zero return code, or an exception while
            launching the subprocess) are not raised -- they are recorded as
            JobRunnerError instances in self.exceptions, and self.retcode
            carries the raw return code. Callers must inspect self.exceptions
            (or self.retcode) after this returns to detect failure.
        """
        env = os.environ.copy()
        if self.has_ompenv: env.update(self.ompenv)

        # Build valgrind command line
        valcmd = ""
        if self.has_valgrind:
            valcmd = f"valgrind --tool={self.valgrind_cmdline} "

        # Perf command
        perf_cmd = ""
        if self.has_perf:
            perf_cmd = f"perf {self.perf_command} "

        # Quoted because these are genuine single filesystem paths (unlike
        # e.g. mpi_args/bin_argstr, which are meant to expand to multiple
        # shell tokens) -- an unquoted path containing a space would
        # otherwise silently split into multiple shell arguments.
        bin_path_q = shlex.quote(bin_path)
        stdin = f" < {shlex.quote(stdin_fname)} " if stdin_fname else ""
        stdout = f" > {shlex.quote(stdout_fname)} " if stdout_fname else ""
        stderr = f" 2> {shlex.quote(stderr_fname)} " if stderr_fname else ""

        if self.has_mpirun or self.has_srun:
            mpirun_np = cast("str", getattr(self, "mpirun_np", ""))
            args = [perf_cmd, mpirun_np, str(mpi_nprocs), f" {self.mpi_args} ",
                    valcmd, bin_path_q, bin_argstr, stdin, stdout, stderr]

        elif self.has_poe:
            # example ${poe} abinit ${poe_args} -procs 4
            # no support for valgrind, debugger, bin_argstr or perf here since poe uses a weird syntax for command line options.
            poe = cast("str", getattr(self, "poe", ""))
            poe_args = cast("str", getattr(self, "poe_args", ""))
            args = [poe, bin_path_q, poe_args, " -procs "+ str(mpi_nprocs),
                    stdin, stdout, stderr]
        else:
            if mpi_nprocs != 1:
                raise ValueError(
                    f"Cannot run with mpi_nprocs={mpi_nprocs}: this JobRunner has no MPI "
                    "launcher configured (mpirun_np and poe are both empty). Set mpirun_np "
                    "(e.g. 'mpiexec -n') in the builder configuration or the [mpi] section "
                    "of the config file."
                )
            args = [perf_cmd, valcmd, bin_path_q, bin_argstr, stdin, stdout, stderr]

        if self.has_debugger:
            # Use completely different syntax if we are running under the control of gdb.
            # mpirun -np 2 xterm -e gdb fftprof --command=dbg_file

            # Get the working directory (warning: I assume that stderr_fname is an absolute path).
            workdir = os.path.dirname(stderr_fname)

            dbg_filepath = os.path.join(workdir, "dbg_commands")
            dbg_filepath_q = shlex.quote(dbg_filepath)

            with open(dbg_filepath, "w") as fh:
                fh.write(f"run {bin_argstr} {stdin}") # Use dbg syntax

            if self.has_mpirun or self.has_srun:
                mpirun_np = cast("str", getattr(self, "mpirun_np", ""))
                args = [mpirun_np, str(mpi_nprocs), "xterm -e gdb", bin_path_q, f"--command={dbg_filepath_q}"]
            else:
                args = ["gdb", bin_path_q, f"--command={dbg_filepath_q}"]

        cmd = " ".join(args)
        #print(cmd)

        #if self.has_valgrind: print(f"Invoking valgrind:\n {cmd}")
        logger.debug(f"About to execute command:\n{cmd}")

        start_time = time.time()
        self.retcode = -1

        try:
            if self.has_timebomb:
                p, self.retcode = self.timebomb.run(cmd, shell=True, cwd=cwd, env=env)
            else:
                p = Popen(cmd, shell=True, cwd=cwd, env=env)
                self.retcode = p.wait()

            run_etime = time.time() - start_time

            if self.retcode != 0:
                exc = JobRunnerError(self.retcode, " ".join(args), run_etime)
                logger.debug(str(exc))
                self.exceptions.append(exc)

        except:
            run_etime = time.time() - start_time
            prev_errmsg = str(sys.exc_info()[1])
            exc = JobRunnerError(self.retcode, " ".join(args), run_etime, prev_errmsg=prev_errmsg)
            self.exceptions.append(exc)

        return run_etime


class BaseValgrindParser:
    """
    Abstract base class for Valgrind output parsers.

    Subclasses must implement the `parse(filename)` method.
    """
    _error_report: str

    # I really miss python 2.6 abc and context managers but must be compatible with py 2.4
    def parse(self, filename: str) -> None:
        raise NotImplementedError("You cannot call the base class")

    @property
    def error_report(self) -> str:
        return self._error_report


class MemcheckParser(BaseValgrindParser):
    """Parser for Valgrind Memcheck tool output."""
    #==3851== HEAP SUMMARY:
    #==3851==     in use at exit: 25,149 bytes in 13 blocks
    #==3851==   total heap usage: 841 allocs, 828 frees, 579,777,815 bytes allocated
    #==3851==
    #==3851== LEAK SUMMARY:
    #==3851==    definitely lost: 0 bytes in 0 blocks
    #==3851==    indirectly lost: 0 bytes in 0 blocks
    #==3851==      possibly lost: 0 bytes in 0 blocks
    #==3851==    still reachable: 25,061 bytes in 12 blocks
    #==3851==         suppressed: 88 bytes in 1 blocks
    #==3851== Rerun with --leak-check=full to see details of leaked memory
    #==3851==
    #==3851== For counts of detected and suppressed errors, rerun with: -v
    #==3851== Use --track-origins=yes to see where uninitialised values come from
    #==3851== ERROR SUMMARY: 10000000 errors from 60 contexts (suppressed: 0 from 0)

    def parse(self, filename: str) -> None:
        """
        Parse the Memcheck output file and store the error report.

        Args:
            filename (str): Path to the output file.

        Raises:
            RuntimeError: If 'LEAK SUMMARY' section is missing.
            ValueError: If a required key is not found in the line.
        """

        def fragile_parser(key: str, string: str) -> int:
            """
            Extract number from a line in the form: key number ignored_tokens

            "Fragile" because this scrapes Valgrind's human-readable text
            output by fixed key position rather than a structured format
            (e.g. XML); any change to that output layout can silently break
            this.
            """
            start = line.find(key)
            if start == -1: raise ValueError(f"Cannot find key {key} in string {string}")
            bytes_lost = int(string[start + len(key):].split(maxsplit=1)[0])
            return bytes_lost

        lost_bytes = 0

        # A `with` block (rather than a bare open()/close()) so the file is
        # closed even when 'LEAK SUMMARY' is missing or fragile_parser()
        # raises below -- both used to leak the handle, since the matching
        # fh.close() at the end was only ever reached on the success path.
        with open(filename) as fh:
            for line in fh:
                if "LEAK SUMMARY:" in line: break
            else:
                raise RuntimeError("Cannot find 'LEAK SUMMARY' section in valgrind stderr file")

            keys = [
                "definitely lost:",
                "indirectly lost:",
                "possibly lost:",
            ]

            # Inspect the next len(keys) line (memleak section)
            errors = {}
            for key, line in zip(keys, fh):
                bytes = fragile_parser(key, line)
                if bytes:
                    errors[key] = bytes

            # Get total number of errors.
            key = "ERROR SUMMARY:"
            for line in fh:
                if key in line:
                    num_errors = fragile_parser(key, line)
                    if num_errors: errors[key] = num_errors

        self._error_report = ""
        if errors: self._error_report = str(errors)


class TimeBomb:
    """Enforces execution timeouts on subprocesses."""

    def __init__(self, timeout: float, delay: float = 0.05, exec_path: str | None = None) -> None:
        """
        Initialize the TimeBomb object.

        Args:
            timeout: Timeout in seconds.
            delay: Delay between checks.
            exec_path: Path to the timeout executable.
        """
        self.timeout = int(timeout)
        self.delay = float(delay)
        self.exec_path = exec_path

    def run(self, args: str | list[str],
            bufsize: int = 0, executable: str | None = None, stdin: Any = None, stdout: Any = None, stderr: Any = None, preexec_fn: Callable[[], Any] | None = None,
            close_fds: bool = False, shell: bool = False, cwd: str | None = None, env: dict[str, str] | None = None, universal_newlines: bool = False, startupinfo: Any = None, creationflags: int = 0) -> tuple[Popen[Any], int]:
        """
        Execute a command with the configured timeout.

        Supports the same interface as subprocess.Popen. Dispatches to one of
        three execution paths depending on configuration:

        1. `exec_path` is set and `timeout > 0`: wraps `args` with the
           external timeout executable (e.g. `timeout <seconds> <args>`) and
           runs it via a plain `Popen`.
        2. `exec_path` is not set and `timeout > 0`: falls back to
           `SubProcessWithTimeout`, which enforces the timeout itself (no
           external timeout executable needed).
        3. `timeout <= 0` (regardless of `exec_path`): no timeout is
           enforced at all -- runs `args` via a plain `Popen`.

        Returns:
            tuple: (subprocess.Popen object, return_code)
        """

        if self.exec_path:
            # timeout exec is available.
            if self.timeout > 0.:
                logger.debug(f"Using timeout function: {self.exec_path}")
                if is_string(args):
                    args = " ".join([self.exec_path, str(self.timeout), cast("str", args)])
                else:
                    args = [self.exec_path, str(self.timeout)] + cast("list[str]", args)

            p = Popen(args,
                      bufsize=bufsize, executable=executable, stdin=stdin, stdout=stdout, stderr=stderr, preexec_fn=preexec_fn,
                      close_fds=close_fds, shell=shell, cwd=cwd, env=env, universal_newlines=universal_newlines, startupinfo=startupinfo,
                      creationflags=creationflags)

            ret_code = p.wait()

        # timeout exec is NOT available.
        elif self.timeout > 0.0:
            logger.debug(f"Using SubprocesswithTimeout and timeout_time: {self.timeout}")
            timeout_proc = SubProcessWithTimeout(self.timeout, delay=self.delay)

            p_temp, ret_code = timeout_proc.run(args,
                bufsize=bufsize, executable=executable, stdin=stdin, stdout=stdout, stderr=stderr, preexec_fn=preexec_fn,
                close_fds=close_fds, shell=shell, cwd=cwd, env=env, universal_newlines=universal_newlines, startupinfo=startupinfo,
                creationflags=creationflags)
            p = cast("Popen[Any]", p_temp)
        else:
            logger.debug("Using Popen (no timeout_time)")
            p = Popen(args,
                      bufsize=bufsize, executable=executable, stdin=stdin, stdout=stdout, stderr=stderr, preexec_fn=preexec_fn,
                      close_fds=close_fds, shell=shell, cwd=cwd, env=env, universal_newlines=universal_newlines, startupinfo=startupinfo,
                      creationflags=creationflags)

            ret_code = p.wait()

        return p, ret_code



class OMPEnvironment(dict):
    """
    Dictionary-like object storing OpenMP environment variables.

    Supports validation of OpenMP-standard keys and initialization from INI files.
    """
    _keys = [
       "OMP_SCHEDULE",
       "OMP_NUM_THREADS",
       "OMP_DYNAMIC",
       "OMP_PROC_BIND",
       "OMP_NESTED",
       "OMP_STACKSIZE",
       "OMP_WAIT_POLICY",
       "OMP_MAX_ACTIVE_LEVELS",
       "OMP_THREAD_LIMIT",
    ]

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        """
        Constructor method inherited from dictionary:

        >>> OMPEnvironment(OMP_NUM_THREADS=1)
        {'OMP_NUM_THREADS': '1'}

        To create an instance from the INI file fname, use:
           OMPEnvironment.from_file(fname)
        """
        self.update(*args, **kwargs)

        err_msg = ""
        for key, value in self.items():
            self[key] = str(value)
            if key not in OMPEnvironment._keys:
                err_msg += f"unknown option {key}"
        if err_msg: raise ValueError(err_msg)

    @classmethod
    def from_file(cls, fname: str, allow_empty: bool = False) -> OMPEnvironment:
        """Initialize the object from file."""
        parser = SafeConfigParser()
        parser.read(fname)

        inst = OMPEnvironment()

        # Consistency check. Note that we only check if the option name is correct,
        # we do not check whether the value is correct or not.
        if "openmp" not in parser.sections():
            if not allow_empty:
                raise ValueError(f"{fname} does not contain any [openmp] section")
            return inst

        err_msg = ""
        for key in parser.options("openmp"):
            if key.upper() not in OMPEnvironment._keys:
                err_msg += f"unknown option {key}, maybe a typo"
        if err_msg:
            raise ValueError(err_msg)

        for key in OMPEnvironment._keys:
            try:
                inst[key] = str(parser.get("openmp", key))
            except NoOptionError:
                try:
                    inst[key] = str(parser.get("openmp", key.lower()))
                except NoOptionError:
                    pass

        if not allow_empty and not inst:
            raise ValueError("Refusing to return with an empty dict")

        return inst
