#!/usr/bin/env python
from __future__ import print_function, division, absolute_import  #, unicode_literals

import sys
import os
import platform

from os.path import join as pj, abspath as absp, basename
from socket import gethostname
from warnings import warn

import logging
logger = logging.getLogger(__name__)

try:
    from ConfigParser import SafeConfigParser, NoOptionError
except ImportError:  # The ConfigParser module has been renamed to configparser in Python 3
    from configparser import ConfigParser as SafeConfigParser, NoOptionError

# Add the directory [...]/abinit/tests to $PYTHONPATH
pack_dir, x = os.path.split(absp(__file__))
pack_dir, x = os.path.split(pack_dir)
sys.path.insert(0, pack_dir)
import tests

__version__ = "0.3"
__author__ = "Matteo Giantomassi"

_my_name = basename(__file__) + "-" + __version__

abitests = tests.abitests
abenv = tests.abenv

from pymods import termcolor
from pymods.testsuite import BuildEnvironment
from pymods.jobrunner import TimeBomb, JobRunner, OMPEnvironment

# Helper functions


def lazy__str__(func):
    "Lazy decorator for __str__ methods"
    def oncall(*args, **kwargs):
        self = args[0]
        return "\n".join(str(k) + " : " + str(v) for (k, v) in self.__dict__.items())
    return oncall


def _yesno2bool(string):
    string = string.lower().strip().replace('"', "").replace("'", "")
    if string == "yes":
        return True
    elif string == "no":
        return False
    else:
        raise ValueError("Cannot interpret string: %s" % string)


def _str2list(string):
    return [s.strip() for s in string.split(",") if s]


class TestBot(object):
    """
    This object drives the execution of the abinit automatic tests:

      1) Read setup options from the file testbot.cfg uploaded by the slave master on the slave.
      2) Initialize the jobrunner and other objects used to run the tests.
      3) run the tests (see the run method), and return the number of tests that failed

    Step 1-2 are done in the creation method.
    """
    _attrbs = {
      # name         -->   (default, parser)
      "slavename"        : (None, str),
      "type"             : ("",   str),  # "ref" if this slave is the reference slave, e.g. testf
      "ncpus"            : (None, int),  # Max number of CPUs that can be used by TestBot
      "mpi_prefix"       : ("",   str),  # MPI runner.
      "mpirun_np"        : ("",   str),  # String used to execute `mpirun -n#NUM`
      "omp_num_threads"  : (0,    int),  # Number of OpenMp threads to be used
      "enable_mpi"       : (None, _yesno2bool),
      "enable_openmp"    : (None, _yesno2bool),
      "poe"              : ("", str),
      "poe_args"         : ("", str),
      "with_tdirs"       : ("", _str2list),   # List of abinit test subsuites to include
      "without_tdirs"    : ("", _str2list),   # List of abinit test subsuites to exclude
      "timeout_time"     : (900, float),      # Timeout time in seconds.
      "cygwin_dir"       : (None, str),
      "runmode"          : ("static", str),
      "keywords"         : ("", _str2list),   # String with the keywords that should be selected/ignored
      "etsf_check"       : ("no", _yesno2bool),  # yes to activate validation of netcdf files produced by Abinit.
    }

    def __init__(self, testbot_cfg=None):

        # 1) Read the options specified in the testbot configuration file.
        if testbot_cfg is None:
            basedir, x = os.path.split(absp(__file__))
            testbot_cfg = pj(basedir, "testbot.cfg")

        parser = SafeConfigParser()
        parser.read(testbot_cfg)

        attrs2read = [
            "slavename",
            "type",
            "ncpus",
            "omp_num_threads",
            "mpi_prefix",
            "mpirun_np",
            "poe",
            "poe_args",
            "with_tdirs",
            "without_tdirs",
            "timeout_time",
            "cygwin_dir",
            "runmode",
            "keywords",
            "etsf_check",
        ]

        for attr in attrs2read:
            default, parse = TestBot._attrbs[attr]
            try:
                value = parser.get("testbot", attr)
            except NoOptionError:
                value = default

            if value is None:
                # Write out the cfg file and raise
                for section in parser.sections():
                    print("[" + section + "]")
                    for opt in parser.options(section):
                        print(opt + " = " + parser.get(section, opt))
                err_msg = "Mandatory option %s is not declared" % attr
                raise ValueError(err_msg)

            self.__dict__[attr] = parse(value)

        if self.with_tdirs and self.without_tdirs:
            err_msg = "with_tdirs and without_tdirs attribute are mutually exclusive"
            raise ValueError(err_msg)

        system, node, release, version, machine, processor = platform.uname()
        print("Running on %s -- slave %s -- system %s -- ncpus %s -- Python %s -- %s" % (
            gethostname(), self.slavename, system, self.ncpus, platform.python_version(), _my_name))

        # Set the logger level.
        # loglevel is bound to the string value obtained from the command line argument.
        # Convert to upper case to allow the user to specify --loglevel=DEBUG or --loglevel=debug
        # numeric_level = getattr(logging, options.loglevel.upper(), None)
        numeric_level = getattr(logging, "ERROR", None)

        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: %s' % numeric_level)
        logging.basicConfig(level=numeric_level)

        # Read testfarm configuration file.
        build_examples = abenv.apath_of(pj("config", "specs", "testfarm.conf"))

        parser = SafeConfigParser()
        parser.read(build_examples)

        if self.slavename not in parser.sections():
            err_msg = "%s is not a valid buildbot slave." % self.slavename
            raise ValueError(err_msg)

        # TODO
        # Consistency check
        # d = self.__dict__
        # for attr, (default, parse) in TestBot._attrbs.items():
        #   try:
        #     d[attr] = parser.get(self.slavename, attr)
        #   except NoOptionError:
        #     print "option %s is not declared" % attr
        #     #raise ValueError(err_msg)
        #     d[attr] = default
        #     if default is None:
        #       err_msg = "option %s is not declared" % attr
        #       raise ValueError(err_msg)
        #   d[attr] = parse(d[attr])

        # 2 )Initialize the jobrunner.
        self.build_env = build_env = BuildEnvironment(os.curdir, cygwin_instdir=self.__dict__["cygwin_dir"])
        self.build_env.set_buildbot_builder(self.slavename)

        if build_env.has_bin("timeout") and self.timeout_time > 0:
            # We can run executables under the control of timeout.c
            timeout_path = build_env.path_of_bin("timeout")
            timebomb = TimeBomb(self.timeout_time, exec_path=timeout_path)
        else:
            warn("Cannot find timeout executable at: %s" % build_env.path_of_bin("timeout"))
            timebomb = TimeBomb(self.timeout_time)

        print("Initalizing JobRunner for sequential runs.")
        self.seq_runner = JobRunner.sequential(timebomb=timebomb)
        print(self.seq_runner)

        if self.has_mpi:
            print("Initalizing MPI JobRunner from self.__dict__")
            self.mpi_runner = JobRunner.fromdict(self.__dict__, timebomb=timebomb)
            print(self.mpi_runner)

        if self.omp_num_threads > 0:
            print("Initalizing OMP environment with omp_num_threads %d " % self.omp_num_threads)
            omp_env = OMPEnvironment(OMP_NUM_THREADS=self.omp_num_threads)
            self.seq_runner.set_ompenv(omp_env)
            if self.has_mpi:
                self.mpi_runner.set_ompenv(omp_env)

        self.targz_fnames = []
        print(self)

        # Initalize the table to store the final results.
        # The table include all the abinit tests (also those that will be skipped)
        # values are initialized with None.
        # FIXME
        database = abitests.get_database()
        res_table = database.init_result_table()
        self.summary = TestBotSummary(res_table)
        # print(self.summary)

    @lazy__str__
    def __str__(self): pass

    @property
    def has_mpi(self):
        """True if we have the MPI runner"""
        return bool(self.mpirun_np) or bool(self.poe)

    @property
    def has_openmp(self):
        """True if the tests must be executed with OpenMP threads."""
        return self.omp_num_threads > 0

    def run_tests_with_np(self, mpi_nprocs, suite_args=None, runmode="static"):
        """
        Run the tests specified by suite_args, using mpi_nprocs MPI processors
        Returns: (nfailed, npassed, nexecuted)
        """
        omp_nthreads = max(self.omp_num_threads, 1)
        py_nprocs = self.ncpus // (mpi_nprocs * omp_nthreads)
        assert py_nprocs > 0

        test_suite = abitests.select_tests(suite_args, keys=self.keywords, regenerate=False)

        # Create workdir.
        workdir = "TestBot_MPI" + str(mpi_nprocs)
        if self.has_openmp:
            workdir += "_OMP" + str(omp_nthreads)

        if os.path.exists(workdir):
            raise RuntimeError("%s already exists!" % workdir)
        else:
            os.mkdir(workdir)

        # Run the tests.
        if self.has_openmp:
            msg = "Running ntests = %s, MPI_nprocs = %s, OMP_nthreads %s, py_nprocs = %s..." % (
                  test_suite.full_length, mpi_nprocs, omp_nthreads, py_nprocs)
        else:
            msg = "Running ntests = %s, MPI_nprocs = %s, py_nprocs = %s..." % (
                  test_suite.full_length, mpi_nprocs, py_nprocs)
        print(msg)

        runner = self.seq_runner
        if mpi_nprocs > 1:
            runner = self.mpi_runner

        results = test_suite.run_tests(
            self.build_env, workdir, runner, make_html_diff=1,
            nprocs=mpi_nprocs, py_nprocs=py_nprocs, runmode=self.runmode)
        # Cannot use this option on the test farm because hdf5 is not thread-safe.
        # See https://www.hdfgroup.org/hdf5-quest.html#tsafe
        # etsf_check=self.etsf_check)

        if results is None:
            print("Test suite is empty, returning 0 0 0 ")
            return 0, 0, 0

        # Store the results in the summary table,
        # taking into account that an input file might be executed multiple times
        # with a different environment (MPI, OMP ...)
        run_info = {}
        # run_info = [self.build_env, workdir, runner, mpi_nprocs, py_nprocs]
        self.summary.merge_results(test_suite, run_info)

        # Push the location of the tarball file
        self.targz_fnames.append(results.targz_fname)

        return results.nfailed, results.npassed, results.nexecuted

    def run(self):
        """
        Run all the automatic tests depending on the environment and the options
        specified in the testbot configuration file.
        Return the number of failing tests (+ no. passed tests if this is the reference slave).
        """
        # If with_tdirs and without_tdirs are not given => execute all tests.
        # else create a list of strings with the suites to (execute|exclude).
        suite_args = None

        # XG130410 Crude hack, to avoid paral and mpiio test directories
        # in case of enable_mpi=no in config/specs/testfarm.conf
        if not self.has_mpi:
            suite_args = "paral- mpiio-".split()

        if self.with_tdirs:
            suite_args = self.with_tdirs
        elif self.without_tdirs:
            # Append "-" to the string to signal that the suite should be excluded.
            suite_args = "- ".join([s for s in self.without_tdirs]) + "-"
            suite_args = suite_args.split()

        # Run tests with 1 processor
        # print(suite_args)
        # runmode = "static"
        # nexecuted = 0

        if self.runmode == "static":
            # Old mode: run all available tests with 1 MPI node here,
            # then use MPI mode in mp_suites (paral/mpiio)
            np_list = [2, 4, 10, 24]
            nfailed, npassed, nexecuted = self.run_tests_with_np(1, suite_args=suite_args, runmode=self.runmode)
        else:
            # New mode: run all available tests with 2 MPI procs here,
            # then enter the mp_suites with [4, 10, 24] CPUs
            np_list = [4, 10, 24]
            nfailed, npassed, nexecuted = self.run_tests_with_np(2, suite_args=suite_args, runmode=self.runmode)

        if self.has_mpi:
            # Run the parallel tests in the multi-parallel suites.
            mp_suites = abitests.multi_parallel_suites()
            mp_names = [suite.name for suite in mp_suites]

            suite_args = mp_names

            # Prune dirs.
            if self.with_tdirs:
                suite_args = [s for s in self.with_tdirs if s in mp_names]

            elif self.without_tdirs:
                suite_args = [s for s in mp_names if s not in self.without_tdirs]

            for np in np_list:
                if np > self.ncpus or not suite_args:
                    continue
                print("Running multi-parallel tests with %s MPI processors, suite_args %s" % (np, suite_args))
                para_nfailed, para_npassed, para_nexec = self.run_tests_with_np(np, runmode="static", suite_args=suite_args)

                # Accumulate the counters.
                nfailed += para_nfailed
                npassed += para_npassed
                nexecuted += para_nexec

        # TODO Collect tarball files.
        # for fname in self.targz_fnames:
        #     print("got targz file: ", fname)

        # Dump the table with the final results.
        # The buildbot master will upload it and the info will be used
        # to generate the HTML summary table
        # for suite_name in self.summary_table:
        #   print suite_name, self.summary_table.suite_status(suite_name)
        from pymods.tools import pprint_table
        table = self.summary.totable()
        pprint_table(table)

        self.summary.json_dump("testbot_summary.json")

        # Empty list of tests (usually due to the use of with, without options)
        # Create file to signal this condition and return 0
        if nexecuted == 0:
            print("No file found")
            with open("__emptylist__", "w") as fh:
                fh.write("nfailed = %d, npassed = %d, nexecuted %d" % (nfailed, npassed, nexecuted))
                return 0

        # Status error depends on the builder type.
        if self.type == "ref":
            return nfailed + npassed  # Reference slave (all the tests must pass).
        else:
            return nfailed


class TestBotSummary(object):
    """Stores the final results of the tests performed by TestBot."""

    _possible_status = ["failed", "passed", "succeeded", "skipped"]

    def __init__(self, res_table):
        self.res_table = res_table
        self.failed = []
        self.passed = []

    @lazy__str__
    def __str__(self): pass

    def _min_status(self, items):
        indices = [self._possible_status.index(item) for item in items]
        return self._possible_status[min(indices)]

    def totable(self):
        """
        Return a table (list of lists whose elements are strings).
        The table has the format:
            - [ suite_name1, suite_name2, ... ]
            - [ stats1,      stats2,      ... ]

        where stats is given by "nfail/npass/nsucces"
        """
        def stats2string(stats):
            "helper function that returns (nfail/npass/nsucc/nskipped)"
            return "/".join(str(stats[k]) for k in self._possible_status)

        table = [self.suite_names()]
        row = []
        for suite_name in self:
            status, stats = self.status_of_suite(suite_name)
            row.append(stats2string(stats))
        table.append(row)

        return table

    def merge_results(self, test_suite, run_info):
        # assert test_suite._executed
        self.run_info = run_info

        for test in test_suite:
            d = self.res_table[test.suite_name][test.id]

            if "status" not in d:
                # Entry is not initialized. Save the status of this test.
                d["status"] = test.status
                d["number_of_runs"] = 1
                d["run_etime"] = test.run_etime
                d["tot_etime"] = test.tot_etime
            else:
                # Handle the case where we have executed the same test but with a different setup.
                # Ordering is: failed < passed < succeeded
                # hence a test is marked as "failed" if at least one run failed.
                d["status"] = self._min_status([test.status, d["status"]])
                d["number_of_runs"] += 1
                d["run_etime"] += test.run_etime
                d["tot_etime"] += test.tot_etime

        # Store the full_id of the tests that (failed|passed).
        for test in test_suite.failed_tests():
            self.failed.append(test.full_id)

        for test in test_suite.passed_tests():
            self.passed.append(test.full_id)

    def __iter__(self):
        """Iterate over the suite names in alphabetical order."""
        for suite_name in self.suite_names():
            yield suite_name

    def suite_names(self):
        """List of suite names in alphabetical order."""
        return sorted(list(self.res_table.keys()))

    def status_of_suite(self, suite_name):
        """
        suite_name : string with the name of the suite.
        return (suite_status, stats) where
          suite_status is one of the possile status in `_possible_status`
          stats is a dictionary : {failed:1, passed:2, succeeded:0, skipped:0}
        """

        # Initialize stats setting the keys to 0
        stats = dict.fromkeys(self._possible_status, 0)

        for test_id, d in self.res_table[suite_name].items():
            if "status" not in d:
                d["status"] = "skipped"
            stats[d["status"]] += 1

        for status in self._possible_status:
            if stats[status] > 0:
                suite_status = status
                break
        else:
            # Ignore error if suite is empty else raise.
            if not self.res_table[suite_name]:
                suite_status = "succeeded"
            else:
                raise RuntimeError("[%s] Wrong list of status values!" % suite_name)

        return suite_status, stats

    def json_dump(self, fname):
        """
        Save self.res_table, self.failed, self.passed in json format.
        The file will be transferred from the slave to the buildbot master.
        """
        def as_ascii(s):
            # return str(s.encode("ascii", "ignore"))
            return str(s)

        d = {}
        # list of strings with the name of the tests that (failed|passed)
        d[as_ascii("failed")] = self.failed
        d[as_ascii("passed")] = self.passed
        d[as_ascii("summary_table")] = self.totable()

        for suite_name in self:
            suite_name = as_ascii(suite_name)
            print("---", suite_name, self.res_table[suite_name])
            if suite_name in d:
                # raise KeyError("Cannot overwrite key %s" % suite_name)
                print("Warning: About to overwrite key %s" % suite_name)
            d[suite_name] = self.res_table[suite_name]

        import json
        with open(fname, "wt") as fh:
            json.dump(d, fh)


def main():
    # Configuration file (hardcoded or from command line)
    testbot_cfg = None
    if len(sys.argv) > 1:
        testbot_cfg = sys.argv[1]

    # Disable colors
    termcolor.enable(False)

    return TestBot(testbot_cfg).run()


if __name__ == "__main__":
    sys.exit(main())
