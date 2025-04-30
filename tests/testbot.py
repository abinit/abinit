#!/usr/bin/env python
from __future__ import print_function, division, absolute_import  #, unicode_literals

import sys
import os
# Set ABI_PSPDIR env variable to point to the absolute path of Pspdir
os.environ["ABI_PSPDIR"] = os.path.abspath(os.path.join(os.path.dirname(__file__), "Pspdir"))
import platform
import shutil
import tempfile
import json

from os.path import join as pj, abspath as absp, basename
from socket import gethostname
from warnings import warn

import logging
logger = logging.getLogger(__name__)

try:
    from ConfigParser import SafeConfigParser, NoOptionError
except ImportError:
    # The ConfigParser module has been renamed to configparser in Python 3
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
from pymods.tools import pprint_table


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
    raise ValueError("Cannot interpret string: %s" % string)


def _str2list(string):
    return [s.strip() for s in string.split(",") if s]


class TestBot(object):
    """
    This object drives the execution of the abinit automatic tests:

      1) Read setup options from the file testbot.cfg uploaded by the master on the worker.
      2) Initialize the job_runner and other objects used to run the tests.
      3) run the tests (see run method), and return the number of tests that failed.

    Step 1-2 are performed in the creation method.
    """
    _attrbs = {
      # name           --> (default, parser, info)
      # If default is None, the option must be specified.
      "slavename"        : (None, str, "Name of buildbot worker"),
      "type"             : ("",   str, "'ref' if this worker is the reference worker where all tests should pass"),
      # TODO: ncpus should be replaced by max_cpus for clarity reasons
      "ncpus"            : (None, int, "Max number of CPUs that can be used by TestBot"),
      "max_gpus"         : (0,    int, "Max number of GPUs that can be used by TestBot"),
      "mpi_prefix"       : ("",   str, "MPI runner"),
      "mpirun_np"        : ("",   str, "String used to execute `mpirun -n#NUM`"),
      "omp_num_threads"  : (0,    int, "Number of OpenMP threads. 0 if OpenMP should not be used"),
      "enable_mpi"       : (None, _yesno2bool, "yes if MPI is activated else no."),
      "enable_openmp"    : (None, _yesno2bool, "yes if OpenMP is activated else no."),
      "poe"              : ("", str, "This option is deprecated"),
      "poe_args"         : ("", str, "This option is deprecated"),
      "with_tdirs"       : ("", _str2list, "List of subsuites to include."),
      "without_tdirs"    : ("", _str2list, "List of subsuites to exclude."),
      "timeout_time"     : (900, float, "Timeout time in seconds."),
      "cygwin_dir"       : ("", str, "This option is deprecated"),
      "runmode"          : ("static", str, "'static to run all tests with 1 MPI proc and use np > 1 only for multiparallel tests'"),
      "keywords"         : ("", _str2list, "String with the keywords that should be selected/ignored."),
      "etsf_check"       : ("no", _yesno2bool, "yes to activate the validation of the netcdf files produced by Abinit."),
      "verbose"          : (0,    int, "Verbosity level"),
      "tmp_basedir"      : ("", str, "temporary folder where the tests will be executed and copied back"),
      "mpi_args"         : ("", str, "args passed to the mpi command"),
      "force_mpi"        : ("no", _yesno2bool, "force usage of of mpirun_np prefix"),
    }

    @classmethod
    def print_options(cls):
        """
        Print the different options supported in testbot.cfg with a
        brief description and the default value.
        """
        for opt_name, t in cls._attrbs.items():
            default, parser, info = t
            print("#", info)
            print(opt_name, "=", default)

        print("# NB If default is None, the option must be specified.")

    def __init__(self, testbot_cfg=None):

        # Read the options specified in the testbot configuration file.
        if testbot_cfg is None:
            basedir, _ = os.path.split(absp(__file__))
            testbot_cfg = pj(basedir, "testbot.cfg")

        parser = SafeConfigParser()
        parser.read(testbot_cfg)

        attrs2read = [
            "slavename",
            "type",
            "ncpus",
            "mpi_prefix",
            "mpirun_np",
            "omp_num_threads",
            "poe",
            "poe_args",
            "with_tdirs",
            "without_tdirs",
            "timeout_time",
            "cygwin_dir",
            "runmode",
            "keywords",
            "etsf_check",
            "verbose",
            "tmp_basedir",
            "mpi_args",
            "force_mpi",
        ]

        for attr in attrs2read:
            default, parse, info = TestBot._attrbs[attr]
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
                raise ValueError("Mandatory option %s is not declared" % attr)

            self.__dict__[attr] = parse(value)

        if self.with_tdirs and self.without_tdirs:
            raise ValueError("with_tdirs and without_tdirs attribute are mutually exclusive")

        # TODO: ncpus should be replaced by max_cpus for clarity reasons
        self.max_cpus = self.ncpus

        system, node, release, version, machine, processor = platform.uname()
        print("Running on %s -- worker %s -- system %s -- max_cpus %s -- Python %s -- %s" % (
              gethostname(), self.slavename, system, self.max_cpus, platform.python_version(), _my_name))

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
            raise ValueError("%s is not a valid buildbot worker." % self.slavename)

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

        # 2 )Initialize the job_runner.
        self.build_env = build_env = BuildEnvironment(os.curdir, cygwin_instdir=self.__dict__["cygwin_dir"])
        self.build_env.set_buildbot_builder(self.slavename)

        # TODO: These parameters should be passed to testbot.cfg
        from tests.pymods.devtools import number_of_cpus, number_of_gpus
        #max_cpus = max(1, number_of_cpus())
        if "HAVE_GPU" not in self.build_env.defined_cppvars:
            self.max_gpus = 0
        else:
            self.max_gpus = max(0, number_of_gpus())

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

        # Initialize the table to store the final results.
        # The table include all the abinit tests (also those that will be skipped)
        # values are initialized with None.
        # FIXME
        database = abitests.get_database()
        res_table = database.init_result_table()
        self.summary = TestBotSummary(res_table)
        #print(self.summary)

    def __str__(self):
        """String representation."""
        lines = []
        app = lines.append
        for attr_name, t in self._attrbs.items():
            default, parser, info = t
            value = getattr(self, attr_name, "undefined")
            app("%s = %s" % (attr_name, value))

        return "\n".join(lines)

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
        Run the tests specified by suite_args, using mpi_nprocs MPI processors.

        Returns: (nfailed, npassed, nexecuted)
        """
        # Compute number of python processes, note that self.omp_num_threads might be zero.
        py_nprocs = self.max_cpus // (mpi_nprocs * max(self.omp_num_threads, 1))
        if py_nprocs < 1:
            raise RuntimeError("py_nprocs = %s" % py_nprocs)

        test_suite = abitests.select_tests(suite_args, keys=self.keywords, regenerate=False)

        # Create workdir.
        workdir_name = "TestBot_MPI" + str(mpi_nprocs)
        if self.has_openmp:
            workdir_name += "_OMP" + str(self.omp_num_threads)

        if os.path.exists(workdir_name):
            raise RuntimeError("%s already exists!" % workdir_name)
        else:
            os.mkdir(workdir_name)

        if self.tmp_basedir:
            workdir = os.path.join(tempfile.mkdtemp(dir=self.tmp_basedir), workdir_name)
        else:
            workdir = workdir_name

        # Run the tests.
        if self.has_openmp:
            msg = "Running ntests = %s, MPI_nprocs = %s, OMP_nthreads %s, Max GPUs: %s, py_nprocs = %s..." % (
                  test_suite.full_length, mpi_nprocs, self.omp_num_threads, self.max_gpus, py_nprocs)
        else:
            msg = "Running ntests = %s, MPI_nprocs = %s, Max GPUs: %s py_nprocs = %s..." % (
                  test_suite.full_length, mpi_nprocs, self.max_gpus, py_nprocs)
        print(msg)

        job_runner = self.seq_runner
        if mpi_nprocs > 1 or self.force_mpi:
            job_runner = self.mpi_runner

        results = test_suite.run_tests(self.build_env, workdir, job_runner,
                                       mpi_nprocs=mpi_nprocs,
                                       omp_nthreads=self.omp_num_threads,
                                       max_cpus=self.max_cpus,
                                       max_gpus=self.max_gpus,
                                       py_nprocs=py_nprocs,
                                       runmode=self.runmode,
                                       verbose=self.verbose,
                                       make_html_diff=1)

        # Cannot use this option on the test farm because hdf5 is not thread/process-safe.
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

        if self.tmp_basedir:
            for fn in ["results.tar.gz", "suite_report.html"]:
                try:
                    shutil.copy2(os.path.join(workdir, fn), workdir_name)
                except:
                    print("Could not copy back file ", fn)

        return results.nfailed, results.npassed, results.nexecuted

    def run(self):
        """
        Run all the automatic tests depending on the environment and the options
        specified in the testbot.cfg configuration file.
        Return the number of failing tests (+ no. passed tests if this is the reference worker).
        """
        # If with_tdirs and without_tdirs are not given => execute all tests.
        # else create a list of strings with the suites that should be executed|excluded.
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

        if self.runmode == "static":
            # Old mode: run all available tests with 1 MPI proc here,
            # then use MPI mode in mp_suites (paral/mpiio)
            np_list = [2, 4, 10, 24, 64]
            #nfailed, npassed, nexecuted = 0, 0 , 0
            nfailed, npassed, nexecuted = self.run_tests_with_np(1, suite_args=suite_args, runmode=self.runmode)
        else:
            # New mode: run all available tests with 2 MPI procs here,
            # then enter the mp_suites with np_list CPUs
            np_list = [4, 10, 24, 64]
            nfailed, npassed, nexecuted = self.run_tests_with_np(2, suite_args=suite_args, runmode=self.runmode)

        if self.has_mpi:
            # Run the parallel tests in the multi-parallel suites.
            mp_suites = abitests.multi_parallel_suites()
            suite_args = [suite.name for suite in mp_suites]

            # Prune dirs.
            if self.with_tdirs:
                suite_args = [s for s in self.with_tdirs if s in suite_args]

            elif self.without_tdirs:
                suite_args = [s for s in suite_args if s not in self.without_tdirs]

            for np in np_list:
                if np > self.max_cpus or not suite_args:
                    print("Skipping tests with np: %s as max_cpus: %s" % (np, self.max_cpus))
                    continue

                print("Running multi-parallel tests with %s MPI processors, suite_args %s" % (np, suite_args))
                para_nfailed, para_npassed, para_nexec = self.run_tests_with_np(np, runmode="static", suite_args=suite_args)

                # Accumulate the counters.
                nfailed += para_nfailed
                npassed += para_npassed
                nexecuted += para_nexec

        table = self.summary.to_table()
        pprint_table(table)

        self.summary.json_dump("testbot_summary.json")

        # Empty list of tests (usually due to the use of with, without options)
        # Create file to signal this condition and return 0
        if nexecuted == 0:
            print("No file found")
            with open("__emptylist__", "wt") as fh:
                fh.write("nfailed = %d, npassed = %d, nexecuted = %d" % (nfailed, npassed, nexecuted))
                return 0

        # The status error depends on the builder type.
        if self.type == "ref":
            # Reference worker --> all the tests must pass.
            return nfailed + npassed
        else:
            return nfailed

    def finalize(self):
        """
        This piece of code has been extracted from analysis9
        """
        fname = "testbot_summary.json"
        with open(fname, "rt") as data_file:
           d = json.load(data_file)

        # FIXME What is this?
        d['tag'] = sys.argv[1]

        with open(fname, 'wt') as data_file:
           json.dump(d, data_file)

        try:
            tests_status = dict(zip(d["summary_table"][0],d["summary_table"][1]))

            dashline = "=========================================================================="
            print( dashline )
            print(     "          Serie   #failed   #passed  #succes  #skip  |   #CPU      #WALL")
            print(dashline)
            rtime = 0.0
            ttime = 0.0
            paral = ''
            mpiio = ''
            for t, s in sorted(tests_status.items()):
                kt = False
                for i in d[t].keys():
                   if  d[t][i]['status'] != "skipped":
                      kt = True
                      rtime += d[t][i]['run_etime']
                      ttime += d[t][i]['tot_etime']
                if kt:
                     temp = ''.join(['%5s   |' % l for l in  s.split('/') ])
                     temp = '%15s | %10s %7.1f  | %7.1f' % (t,temp,rtime,ttime)
                     if t == 'mpiio':
                        mpiio = temp
                     elif t == 'paral':
                        paral = temp
                     else:
                        print(temp)
                rtime = ttime = 0.0

            print(dashline)
            putline = 0
            if paral != '':
                print(paral)
                putline=1
            if mpiio != '':
                print(mpiio)
                putline=1
            if putline == 1:
                print(dashline)
        except:
            print("no results")
            sys.exit(99)


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

    def to_table(self):
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
        rows = []
        for suite_name in self:
            status, stats = self.status_of_suite(suite_name)
            rows.append(stats2string(stats))
        table.append(rows)

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
        Args:
            suite_name: string with the name of the suite.

        return: (suite_status, stats) where
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
        The file will be transferred from the buildbot worker to the buildbot master.
        """
        # list of strings with the name of the tests that (failed|passed)
        d = {}
        d["failed"] = self.failed
        d["passed"] = self.passed
        d["summary_table"] = self.to_table()

        for suite_name in self:
            print("---", suite_name, self.res_table[suite_name])
            if suite_name in d:
                #raise KeyError("Cannot overwrite key %s" % suite_name)
                print("Warning: About to overwrite key %s" % suite_name)
            d[suite_name] = self.res_table[suite_name]


        with open(fname, "wt") as fh:
            json.dump(d, fh)


def main():
    if "--help" in sys.argv or "-h" in sys.argv:
        # Print help and exit.
        TestBot.print_options()
        return 0

    # Configuration file (hardcoded or from command line)
    testbot_cfg = None
    if len(sys.argv) > 1:
        testbot_cfg = sys.argv[1]
        print("Reading testbot.cfg configuration file from: ", testbot_cfg)

    # Disable colors
    termcolor.enable(False)

    testbot = TestBot(testbot_cfg)
    if "--dry-run" in sys.argv or "-d" in sys.argv:
        print("Running in dry-run mode, will return immediately.")
        print(testbot)
        return 0

    return testbot.run()


if __name__ == "__main__":
    #import multiprocessing
    #multiprocessing.set_start_method("fork")  # Ensure compatibility on macOS/Linux
    sys.exit(main())
