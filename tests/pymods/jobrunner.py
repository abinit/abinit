"""
Job execution and management tools for ABINIT tests.

This module provides a unified interface for running ABINIT binaries in various
environments, including sequential, MPI (OpenMPI, MPICH, srun, poe), and OpenMP.
It also includes facilities for enforcing timeouts, profiling with Valgrind or
perf, and managing environment variables.
"""
import os
import sys
import time

from configparser import NoOptionError
from subprocess import Popen

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


CFG_KEYWORDS = {
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


def is_string(s):
    """
    Check if the input is a string-like object.

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


def mpicfg_parser(fname, defaults=None):
    """
    Parse a configuration file (INI format) for MPI options.

    Args:
        fname: Path to the configuration file.
        defaults: Default values for the parser.

    Returns:
        dict: A dictionary containing the parsed MPI options.
    """
    logger.debug("Parsing [MPI] section in file : " + str(fname))

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
        except:
            raise ValueError("Wrong line: key = " + str(key) + " d[key] = " + str(d[key]) )

    return d

#NB: Pickle fail if JobRunnerError inherits from Exception
# so we inherit from object. This should represents a serious problem
# because we never catch JobRunner Exceptions. This object is mainly
# used to store info about the exception in JobRunner exceptions (see run method)

class JobRunnerError:
    """Exception-like object to store information about job execution failures."""

    def __init__(self, return_code, cmd, run_etime, prev_errmsg=None):
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

    def __str__(self):
        string = "Command %s\n returned exit_code: %s\n" % (self.cmd, self.return_code)
        if self.prev_errmsg:
            string += "Previous exception: %s" % self.prev_errmsg

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
    def fromdict(cls, kwargs, ompenv=None, timebomb=None):
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
    def fromfile(cls, fname, timebomb=None):
        """
        Create a JobRunner instance from an INI configuration file.

        Args:
            fname: Path to the configuration file.
            timebomb: Optional TimeBomb instance.

        Returns:
            JobRunner: A new instance configured from the file.
        """
        d = mpicfg_parser(fname)
        d["ompenv"] = OMPEnvironment.from_file(fname, allow_empty=True)
        d["timebomb"] = timebomb

        return cls(d)

    @classmethod
    def sequential(cls, ompenv=None, timebomb=None):
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
    def srun(cls, ompenv=None, timebomb=None, mpi_args=""):
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
    def generic_mpi(cls, ompenv=None, use_mpiexec=False, mpi_args="", timebomb=None):
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

    def __init__(self, dic):
        """
        Initialize the JobRunner.

        Args:
            dic (dict): Dictionary of configuration options.

        Raises:
            ValueError: If poe and (mpirun or srun) are both specified.
        """
        self.exceptions = []

        for k, v in dic.items():
            if k not in self.__dict__:
                self.__dict__[k] = v
            else:
                raise ValueError("key %s is already in self.__dict__, cannot overwrite" % k)

        if "mpi_args" not in dic:
            self.mpi_args = ""

        if self.has_poe and (self.has_mpirun or self.has_srun):
            raise ValueError("poe and (mpirun||srun) are mutually exclusive")

    def __str__(self):
        """
        Return a string representation of the job runner's configuration.

        Returns:
            str: The configuration summary.
        """
        string = ""
        for key in CFG_KEYWORDS:
            attr = getattr(self, str(key), None)
            if attr: string += "%s = %s\n" % (key, attr)

        if string:
            string = "[MPI setup]\n" + string

        if self.has_ompenv:
            string += "[OpenMP]\n" + str(self.ompenv)

        return string

    def set_timebomb(self, timebomb):
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

    def set_valgrind_cmdline(self, cmdline):
        """
        Set the command line options to be passed to VALGRIND.

        Args:
            cmdline (str): The command line options.
        """
        self.valgrind_cmdline = cmdline

    @property
    def has_valgrind(self):
        """True if we are running the code with VALGRIND."""
        return hasattr(self, "valgrind_cmdline")

    def build_valgrind_parser(self):
        """
        Build and return a parser for Valgrind output.

        Returns:
            MemcheckParser: The initialized parser.

        Raises:
            ValueError: If Valgrind is not enabled for this runner.
        """
        if not self.has_valgrind: raise ValueError("Runner does not use valgrind!")
        return MemcheckParser()

    def set_debugger(self, debugger):
        """
        Set the debugger executable.

        Args:
            debugger (str): Path to the debugger executable.
        """
        self.debugger = debugger

    @property
    def has_debugger(self):
        """True if we are running the executable under the control of a debugger."""
        return hasattr(self, "debugger")

    def set_perf_command(self, perf_command):
        """
        Set the perf command to be used.

        Args:
            perf_command (str): The perf command string.
        """
        self.perf_command = perf_command

    @property
    def has_perf(self):
        """True if we are profiling the run with perf."""
        return hasattr(self, "perf_command")

    @property
    def has_srun(self):
        """True if we are running with Slurm srun"""
        return hasattr(self, "mpirun_np") and self.mpirun_np == "srun -n"

    @property
    def has_mpirun(self):
        """True if we are running a MPI job with mpirun"""
        return hasattr(self, "mpirun_np") and self.mpirun_np != "srun -n"


    @property
    def has_poe(self):
        """True if are using IBM poe for MPI executions."""
        return hasattr(self, "poe") and bool(self.poe)

    @property
    def has_timebomb(self):
        """
        True if we are running the job under the control of
        an application that will enforce a timeout.
        """
        return hasattr(self, "timebomb") and bool(self.timebomb)

    @property
    def has_ompenv(self):
        """True if we are using OpenMP."""
        return hasattr(self, "ompenv") and bool(self.ompenv)

    def set_ompenv(self, ompenv):
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

    def run(self, mpi_nprocs, bin_path, stdin_fname, stdout_fname, stderr_fname, bin_argstr="", cwd=None):
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
        """
        env = os.environ.copy()
        if self.has_ompenv: env.update(self.ompenv)

        # Build valgrind command line
        valcmd = ""
        if self.has_valgrind:
            valcmd = "valgrind --tool=%s " % self.valgrind_cmdline

        # Perf command
        perf_cmd = ""
        if self.has_perf:
            perf_cmd = "perf %s " % self.perf_command

        stdin = " < %s " % stdin_fname if stdin_fname else ""
        stdout = " > %s " % stdout_fname if stdout_fname else ""
        stderr = " 2> %s " % stderr_fname if stderr_fname else ""

        if self.has_mpirun or self.has_srun:
            args = [perf_cmd, self.mpirun_np, str(mpi_nprocs), " %s " % self.mpi_args,
                    valcmd, bin_path, bin_argstr, stdin, stdout, stderr]

        elif self.has_poe:
            # example ${poe} abinit ${poe_args} -procs 4
            # no support for valgrind, debugger, bin_argstr or perf here since poe uses a weird syntax for command line options.
            args = [self.poe, bin_path, self.poe_args, " -procs "+ str(mpi_nprocs),
                    stdin, stdout, stderr]
        else:
            assert mpi_nprocs == 1
            args = [perf_cmd, valcmd, bin_path, bin_argstr, stdin, stdout, stderr]

        if self.has_debugger:
            # Use completely different syntax if we are running under the control of gdb.
            # mpirun -np 2 xterm -e gdb fftprof --command=dbg_file

            # Get the working directory (warning: I assume that stderr_fname is an absolute path).
            workdir = os.path.dirname(stderr_fname)

            dbg_filepath = os.path.join(workdir, "dbg_commands")

            with open(dbg_filepath, "w") as fh:
                fh.write("run %s %s" % (bin_argstr, stdin)) # Use dbg syntax

            if self.has_mpirun or self.has_srun:
                args = [self.mpirun_np, str(mpi_nprocs), "xterm -e gdb", bin_path, "--command=%s" % dbg_filepath]
            else:
                args = ["gdb", bin_path, "--command=%s" % dbg_filepath]

        cmd = " ".join(args)
        #print(cmd)

        #if self.has_valgrind: print("Invoking valgrind:\n %s" % cmd)
        logger.debug("About to execute command:\n" + cmd)

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
    # I really miss python 2.6 abc and context managers but must be compatible with py 2.4
    def parse(self, filename):
        raise NotImplementedError("You cannot call the base class")

    @property
    def error_report(self):
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

    def parse(self, filename):
        """
        Parse the Memcheck output file.

        Args:
            filename (str): Path to the output file.

        Returns:
            dict: Dictionary of detected errors.

        Raises:
            RuntimeError: If 'LEAK SUMMARY' section is missing.
            ValueError: If a required key is not found in the line.
        """

        def fragile_parser(key, string):
            """
            Extract number from a line in the form: key number ignored_tokens
            """
            start = line.find(key)
            if start == -1: raise ValueError("Cannot find key %s in string %s" % (key, string))
            bytes_lost = int(string[start + len(key):].split()[0])
            return bytes_lost

        lost_bytes = 0
        fh = open(filename)

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

        fh.close()

        self._error_report = ""
        if errors: self._error_report = str(errors)


class TimeBomb:
    """Enforces execution timeouts on subprocesses."""

    def __init__(self, timeout, delay=.05, exec_path=None):
        """
        Initialize the TimeBomb object.

        Args:
            timeout (int): Timeout in seconds.
            delay (float, optional): Delay between checks.
            exec_path (str, optional): Path to the timeout executable.
        """
        self.timeout = int(timeout)
        self.delay = float(delay)
        self.exec_path = exec_path

    def run(self, args,
            bufsize=0, executable=None, stdin=None, stdout=None, stderr=None, preexec_fn=None,
            close_fds=False, shell=False, cwd=None, env=None, universal_newlines=False, startupinfo=None, creationflags=0):
        """
        Execute a command with the configured timeout.

        Supports the same interface as subprocess.Popen.

        Returns:
            tuple: (subprocess.Popen object, return_code)
        """
        try:

            if self.exec_path:
                #
                # timeout exec is available.
                #
                if self.timeout > 0.:
                    logger.debug("Using timeout function: " + self.exec_path)
                    if is_string(args):
                        args = " ".join([self.exec_path, str(self.timeout), args])
                    else:
                        args = [self.exec_path, str(self.timeout)] + args

                p = Popen(args,
                          bufsize=bufsize, executable=executable, stdin=stdin, stdout=stdout, stderr=stderr, preexec_fn=preexec_fn,
                          close_fds=close_fds, shell=shell, cwd=cwd, env=env, universal_newlines=universal_newlines, startupinfo=startupinfo,
                          creationflags=creationflags)

                ret_code = p.wait()

            #
            # timeout exec is NOT available.
            #
            elif self.timeout > 0.0:
                logger.debug("Using SubprocesswithTimeout and timeout_time : "+str(self.timeout))
                p = SubProcessWithTimeout(self.timeout, delay=self.delay)

                p, ret_code = p.run(args,
                    bufsize=bufsize, executable=executable, stdin=stdin, stdout=stdout, stderr=stderr, preexec_fn=preexec_fn,
                    close_fds=close_fds, shell=shell, cwd=cwd, env=env, universal_newlines=universal_newlines, startupinfo=startupinfo,
                    creationflags=creationflags)
            else:
                logger.debug("Using Popen (no timeout_time)")
                p = Popen(args,
                          bufsize=bufsize, executable=executable, stdin=stdin, stdout=stdout, stderr=stderr, preexec_fn=preexec_fn,
                          close_fds=close_fds, shell=shell, cwd=cwd, env=env, universal_newlines=universal_newlines, startupinfo=startupinfo,
                          creationflags=creationflags)

                ret_code = p.wait()

            return p, ret_code

        except Exception:
            raise


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
       "OMP_STACKSIZE",
       "OMP_PROC_BIND",
    ]

    def __init__(self, *args, **kwargs):
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
                err_msg += "unknown option %s" % key
        if err_msg: raise ValueError(err_msg)

    @classmethod
    def from_file(cls, fname, allow_empty=False):
        """Initialize the object from file."""
        parser = SafeConfigParser()
        parser.read(fname)

        inst = OMPEnvironment()

        # Consistency check. Note that we only check if the option name is correct,
        # we do not check whether the value is correct or not.
        if "openmp" not in parser.sections():
            if not allow_empty:
                raise ValueError("%s does not contain any [openmp] section" % fname)
            return inst

        err_msg = ""
        for key in parser.options("openmp"):
            if key.upper() not in OMPEnvironment._keys:
                err_msg += "unknown option %s, maybe a typo" % key
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
