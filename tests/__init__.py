from __future__ import annotations

"""
The ABINIT test suite package.

This package provides classes and utilities for managing the ABINIT automatic
test suite, including environment discovery, test selection, and database
management.
"""

import json
import logging
import os
import pickle
import platform
import re
import sys
from collections.abc import Callable, Iterable, Iterator
from io import StringIO
from pprint import pprint
from socket import gethostname
from typing import Any

from tests.pymods.devtools import FileLock
from tests.pymods.termcolor import cprint
from tests.pymods.testsuite import AbinitTestSuite, ChainOfTests

logger = logging.getLogger(__name__)


__all__ = []

# Dictionary mapping test keywords to string with human-readable description.
_json_path = os.path.join(os.path.dirname(__file__), "known_keywords.json")
with open(_json_path, encoding="utf-8") as _fh:
    KNOWN_KEYWORDS = json.load(_fh)


class AbinitEnvironment:
    """
    Environmental context and path manager for the ABINIT source tree.

    Provides information about the system (hostname, username) and helper
    methods to construct absolute paths relative to the ABINIT root directory.
    """

    def __init__(self):
        self.uname = platform.uname()
        self.hostname = gethostname()
        try:
            self.username = os.getlogin()
        except Exception:
            self.username = "No_username"

        filepath = os.path.abspath(__file__)
        # Follows symlink (if any)
        filepath = os.path.realpath(filepath)

        # Paths of important dirs.
        self.tests_dir, tail = os.path.split(filepath)
        self.home_dir, tail = os.path.split(self.tests_dir)

        self.src_dir = os.path.join(self.home_dir, "src")
        self.psps_dir = os.path.join(self.tests_dir, "Pspdir")
        self.fldiff_path = os.path.join(self.tests_dir, "Scripts", "fldiff.pl")

    def __str__(self):
        """
        Return a summary of the environment.

        Returns:
            str: Environment details.
        """
        return "\n".join([str(k) + " : " + str(v) for (k, v) in self.__dict__.items()])

    def apath_of(self, *p: str) -> str:
        """
        Compute the absolute path of a subpath relative to the ABINIT root.

        Args:
            *p: One or more path components.

        Returns:
            str: The joined absolute path.
        """
        return os.path.join(self.home_dir, *p)

    def isbuild(self) -> bool:
        """bool: True if the code has been built in the current directory."""
        configh_path = os.path.join(self.home_dir, "config.h")
        abinit_path = os.path.join(self.home_dir, "src", "98_main", "abinit")
        return os.path.isfile(configh_path) and os.path.isfile(abinit_path)

    def touch_srcfiles(self, patterns: list[str]) -> None:
        """
        Touch (update timestamp) all source files containing any of the patterns.

        Args:
            patterns: List of string patterns to search for in source files.
        """

        def touch(fname: str) -> None:
            """
            Python touch
            See also http://stackoverflow.com/questions/1158076/implement-touch-using-python
            for a race-condition free version based on py3k
            """
            try:
                os.utime(fname, None)
            except Exception:
                open(fname, "a").close()

        top = os.path.join(abenv.home_dir, "src")
        print("Walking through directory:", top)
        print("Touching all files containing:", str(patterns))

        for root, dirs, files in os.walk(self.src_dir):
            for path in files:
                if not path.lower().endswith(".f90"):
                    continue
                path = os.path.join(root, path)
                with open(path) as fh:
                    # print(path)
                    for line in fh:
                        if any(p in line for p in patterns):
                            print("Touching %s" % path)
                            touch(path)
                            break

    def start_watching_sources(self) -> dict[str, os.stat_result]:
        """
        Initialize a snapshot of the source files' statistics for modification tracking.

        Returns:
            dict: Mapping of file paths to their `os.stat` results.
        """

        def is_source(fname: str) -> bool:
            # NB: `.finc` include files are ignored on purpose because
            # make is not tracking them. one should change the Fortran file
            # that includes .finc to force recompilation.
            _, ext = os.path.splitext(fname)
            return ext.lower() in (".f", ".f90", ".c")  # ".finc"

        self.srcpath_stat = {}
        for root, dirs, files in os.walk(self.src_dir):
            for fname in files:
                if not is_source(fname):
                    continue
                path = os.path.join(root, fname)
                self.srcpath_stat[path] = os.stat(path)

        return self.srcpath_stat

    def changed_sources(self) -> list[str]:
        """
        Detect and return a list of source files modified since the last check.

        Returns:
            list: Paths of modified files.
        """
        changed = []
        for path, old_stat in self.srcpath_stat.items():
            now_stat = os.stat(path)
            if now_stat.st_mtime != old_stat.st_mtime:
                changed.append(path)
                self.srcpath_stat[path] = now_stat
        return changed


abenv = AbinitEnvironment()

database_path = os.path.join(abenv.tests_dir, "test_suite.cpkl")


_tsuite_dirs = [
    "atompaw",
    "atdep",
    "bigdft",
    "bigdft_paral",
    "built-in",
    # "cpu",      This directory is disabled
    "etsf_io",
    "fast",
    "gwr_suite",
    "gwpt_suite",
    "psml",
    "gpu",
    "libxc",
    "mpiio",
    "paral",
    # "hpc",
    "hpc_gpu_omp",
    "gpu_omp",
    "gpu_kokkos",
    "rttddft_suite",
    "seq",
    "tutoatdep",
    "tutomultibinit",
    "tutoparal",
    "tutoplugs",
    "tutorespfn",
    "tutorial",
    "unitary",
    "v1",
    "v2",
    "v3",
    "v4",
    "v5",
    "v6",
    "v67mbpt",
    "v7",
    "v8",
    "v9",
    "v10",
    "vdwxc",
    "wannier90",
]

_tsuite_dirs.sort()
_tsuite_dirs = tuple(
    [os.path.join(abenv.tests_dir, dir_name) for dir_name in _tsuite_dirs]
)


def load_mod(filepath: str) -> Any:
    """
    Dynamically load a Python module from a file path.

    Supports both older `imp` and modern `importlib` mechanisms.

    Args:
        filepath: Absolute path to the .py file.

    Returns:
        ModuleType: The loaded module.
    """
    try:
        import imp

        return imp.load_source(filepath, filepath)
    except ModuleNotFoundError:
        from importlib.machinery import SourceFileLoader

        return SourceFileLoader(filepath, filepath).load_module()


class Suite:
    """
    Representation of an ABINIT test suite.

    A suite corresponds to a directory containing an `__init__.py` file and
    an `Input/` subdirectory. It manages keywords, required CPP variables,
    and subsuite logic.
    """

    def __init__(self, suite_path: str):
        suite_path = os.path.abspath(suite_path)

        self.suite_path = os.path.abspath(suite_path)
        self.name = os.path.basename(suite_path)

        module_name = os.path.join(suite_path, "__init__.py")
        module = load_mod(os.path.join(suite_path, "__init__.py"))

        self.keywords = set(module.keywords)
        self.need_cpp_vars = set(module.need_cpp_vars)

        # Divide tests into (active|disabled).
        self.inp_paths = [p for p in module.inp_files if not p.startswith("-")]
        self.disabled_inp_paths = [p[1:] for p in module.inp_files if p.startswith("-")]

        # Use absolute paths
        self.inp_paths = [os.path.join(suite_path, "Input", p) for p in self.inp_paths]
        self.disabled_inp_paths = [
            os.path.join(suite_path, "Input", p) for p in self.disabled_inp_paths
        ]

        # True if the suite contains tests that should be executed with different numbers of MPI processes
        self.is_multi_parallel = False
        if hasattr(module, "is_multi_parallel"):
            self.is_multi_parallel = module.is_multi_parallel

        self.subsuites = {}
        if hasattr(module, "subsuites"):
            subsuite_names = module.subsuites
            for k in subsuite_names:
                self.subsuites[k] = []

            for name in subsuite_names:
                pattern = re.compile("-?t" + name + r"_\d+\.abi")
                for inp in module.inp_files:
                    if pattern.match(inp):
                        # print(inp, "--> subsuite: ", name)
                        inp_path = os.path.join(suite_path, "Input", inp)
                        self.subsuites[name].append(inp_path)

            nfound = sum([len(paths) for paths in self.subsuites.values()])
            if nfound != len(module.inp_files):
                err_msg = (
                    "At least one input_file does not belong to a subsuite, nfound = %s, __init__.nfiles = %s\n"
                    % (nfound, len(module.inp_files))
                )

                for inp in module.inp_files:
                    for paths in self.subsuites.values():
                        fnames = [os.path.basename(p) for p in paths]
                        if inp in fnames:
                            break
                    else:
                        err_msg += "%s not found\n" % inp

                raise ValueError(err_msg)

            # Remove disabled tests (if any).
            for sub_name, paths in self.subsuites.items():
                inp_paths = [p for p in paths if not p.endswith("-")]
                self.subsuites[sub_name] = inp_paths

    def has_subsuite(self, subsuite_name: str) -> bool:
        """bool: True if the suite contains a subsuite with the given name."""
        return subsuite_name in self.subsuites

    def inputs_of_subsuite(self, subsuite_name: str) -> list[str]:
        """
        Return the absolute paths of the input files in a subsuite.

        Args:
            subsuite_name: Name of the subsuite.

        Returns:
            list: List of absolute paths to .abi files.
        """
        return self.subsuites[subsuite_name]

    def __str__(self):
        return "\n".join([str(k) + " : " + str(v) for k, v in self.__dict__.items()])


class AbinitTestsDatabase(dict):
    """
    Database of test instances, indexed by suite name.

    Inherits from `dict`. Each value is typically a list of test objects
    (AbinitTestInfo or similar). Providing facilities to query and iterate
    over the entire test collection.
    """

    def __init__(self, suites: dict[str, Suite]):
        dict.__init__(self)
        self._suites = suites

    def iter_tests(self) -> Iterator[Any]:
        """
        Ordered iterator over all tests in the database.

        Yields:
            AbinitTest: The next test object.
        """
        for suite in self.values():
            for test in suite:
                yield test

    @property
    def suite_names(self) -> list[str]:
        """List of suite names"""
        return [suite.name for suite in self._suites.values()]

    @property
    def authors_snames(self) -> set[str]:
        """List of authors' second names extracted from the tests"""
        all_snames = []
        for test in self.iter_tests():
            all_snames.extend(test._authors_snames)

        return set(all_snames)

    def test_chains(self) -> list[ChainOfTests]:
        """
        Return a list with all the chained tests.

        Returns:
            list[ChainOfTests]: List of tests that are part of a chain.
        """
        return [t for t in self.iter_tests() if isinstance(t, ChainOfTests)]

    def tests_with_variables(self, ivars: Iterable[str]) -> list[Any]:
        """
        Return a list with all the tests that contain the input variables ivars.

        Args:
            ivars: List of input variable names to search for.

        Returns:
            list: List of test objects containing those variables.
        """
        return [t for t in self.iter_tests() if t.has_variables(ivars)]

    def add_test_suite(self, suite_name: str, test_suite: AbinitTestSuite) -> None:
        """
        Add an executed or selected test suite to the database.

        Args:
            suite_name: Key used to index the suite.
            test_suite: The test suite object to add.

        Raises:
            ValueError: If the suite name is invalid or already exists.
        """
        if suite_name not in self.suite_names:
            raise ValueError("%s is not a valid suite name" % suite_name)

        if suite_name in self:
            raise ValueError("%s is already in the database" % suite_name)

        self[suite_name] = test_suite

    def init_result_table(self) -> dict[str, dict[str, dict]]:
        """
        Initialize a nested dictionary used to store test results.

        Structure: `res_table[suite_name][test_id] = {}`.
        Ensures all possible tests have an entry for Buildbot output aggregation.

        Returns:
            dict: The initialized nested results table.

        Raises:
            ValueError: If duplicate test IDs are found within the same suite.
        """
        res_table = {}
        for suite_name in self.suite_names:
            res_table[suite_name] = {}

        for suite_name, suite in self.items():
            for test in suite:
                # test.id should be unique inside a suite. Check it once again.
                if test.id in res_table[suite_name]:
                    other = res_table[suite_name][test.id]
                    print("test\n:", test)
                    print("other:\n", other)
                    raise ValueError(
                        "Replicated test.id %s in suite %s" % (test.id, suite_name)
                    )

                res_table[suite_name][test.id] = {}

        return res_table

    def get_test_suite(
        self,
        suite_name: str,
        subsuite_name: str | None = None,
        slice_obj: slice | None = None,
    ) -> AbinitTestSuite:
        """
        Retrieve a selection of tests from a suite.

        Args:
            suite_name: Name of the suite to query.
            subsuite_name: Optional subsuite filter.
            slice_obj: Optional slice to select a range of tests.

        Returns:
            AbinitTestSuite: An object containing the requested tests.
        """
        test_suite = self[suite_name]

        if subsuite_name is not None:
            # Build the tests in the subsuite from the input files.
            suite = self._suites[suite_name]

            if not suite.has_subsuite(subsuite_name):
                raise ValueError(
                    "suite %s does not have subsuite %s" % (suite_name, subsuite_name)
                )

            sub_inputs = suite.inputs_of_subsuite(subsuite_name)

            abenv = test_suite.abenv
            test_suite = AbinitTestSuite(
                abenv,
                inp_files=sub_inputs,
                keywords=suite.keywords,
                need_cpp_vars=suite.need_cpp_vars,
            )

        if slice_obj is None:
            return test_suite
        logger.debug("will slice test_suite with slice_obj= %s " % slice_obj)
        return test_suite[slice_obj]

    def find_unknown_wrong_keywords(self) -> tuple[set[str], set[str]]:
        """
        Identify keywords in `TEST_INFO` sections that are undocumented or malformed.

        Returns:
            tuple: (unknown_keywords: set, malformed_keywords: set)
        """
        unknowns, wrong = set(), set()
        for suite_name, suite in self.items():
            for test in suite:
                for key in test.keywords:
                    if key not in KNOWN_KEYWORDS:
                        unknowns.add(key)
                    if " " in key:
                        wrong.add(key)

        return unknowns, wrong

    def find_stale_or_lost_inputs(self) -> str:
        """
        Verify that all input files in `Input/` directories are referenced by tests.

        Returns:
            str: Error report with unreferenced or double-referenced input files.
        """
        # Build the list of files that are tested.
        err = StringIO()

        for suite_name, suite in self.items():
            # List all the files in suite_name/Input (exclude hidden files or vim backup files).
            inp_dir = abenv.apath_of("tests", suite_name, "Input")

            listdir = [f for f in os.listdir(inp_dir) if not exclude_path(f)]
            inp_fnames = [os.path.join(inp_dir, f) for f in listdir]

            # Mapping inp_fname --> number of tests using it.
            inp2test = dict().fromkeys(inp_fnames, 0)

            for test in suite:
                for ius in test.inputs_used:
                    if ius not in inp2test:
                        raise ValueError(
                            "Input [%s] [%s] does not appear in Input2keys!"
                            % (suite_name, ius)
                        )
                    inp2test[ius] += 1

            def remove_file(fname):
                # XG130810 : When the report.in files in abirules/Input/report.in  buildsys/Input/report.in
                # will have been suppressed, one might replace the next line by the simpler :
                # return fname.endswith(".files")
                return fname.endswith(".files") or os.path.basename(fname) in [
                    "report.in"
                ]

            keys = []
            for fname, ntimes in inp2test.items():
                if not remove_file(fname):
                    keys.append(fname)
                else:
                    ntest = ntimes
                    assert ntest == 0

            # inp2test = {k: inp2test[k] for k in keys} # requires py2.7
            inp2test = dict([(k, inp2test[k]) for k in keys])

            # At this point inp2test should be >= 1.
            for fname, ntimes in inp2test.items():
                # if ntimes != 1:
                if ntimes == 0:
                    err.write(
                        "Input file %s is used %s time(s)\n" % (path2str(fname), ntimes)
                    )

        return err.getvalue()

    def find_stale_or_lost_refs(self) -> str:
        """
        Verify that all reference files in `Refs/` directories are tracked by at least one test.

        Returns:
            str: Error report with unreferenced reference files.
        """
        # Build the list of files that are tested.
        err = StringIO()

        for suite_name, suite in self.items():
            # List all the files in suite_name/Refs (ignore hidden files).
            ref_dir = abenv.apath_of("tests", suite_name, "Refs")

            if not os.path.exists(ref_dir):
                err.write("%s does not exist\n" % ref_dir)
                continue

            listdir = [f for f in os.listdir(ref_dir) if not exclude_path(f)]
            # use absolute path.
            ref_fnames = [os.path.join(ref_dir, f) for f in listdir]

            # Mapping ref_fname --> number of tests using it.
            ref2test = dict().fromkeys(ref_fnames, 0)

            for test in suite:
                files_to_test = [
                    os.path.join(ref_dir, f.name) for f in test.files_to_test
                ]

                for o in files_to_test:
                    # FIXME due to out --> stdout replacement
                    if o.endswith(".stdout"):
                        o = o[:-7] + ".out"
                    if o not in ref2test:
                        err.write("files_to_test %s does not appear in Refs!\n" % o)
                    else:
                        ref2test[o] += 1

            # At this point ref2test should contain only ones.
            for ref_fname, ntimes in ref2test.items():
                if ntimes != 1:
                    err.write(
                        "Reference file %s is tested %s time(s)\n"
                        % (path2str(ref_fname), ntimes)
                    )

        return err.getvalue()

    def check_testinfo_options(self) -> str:
        """
        Test the presence of important options in the TEST_INFO section of each test.
        """

        def check_options_in_test(test: Any) -> dict[str, str]:
            recommended_opts = [
                "keywords",
                "description",
                "authors",
                "max_nprocs",
            ]

            d = {}
            for opt in recommended_opts:
                try:
                    value = getattr(test, opt)
                    if not value:
                        d[opt] = "EMPTY"
                except AttributeError:
                    d[opt] = "MISSING"

            return d

        lines = []
        app = lines.append
        for suite_name, suite in self.items():
            for test in suite:
                if not isinstance(test, ChainOfTests):
                    wrong_options = check_options_in_test(test)
                    for opt, stat in wrong_options.items():
                        app("%s: option %s is %s" % (test.full_id, opt, stat))
                else:
                    # print("In test chain %s" % test.full_id)
                    for t in test:
                        wrong_options = check_options_in_test(t)
                        for opt, stat in wrong_options.items():
                            app("%s: option %s is %s" % (t.full_id, opt, stat))

        return "\n".join(lines)


def exclude_path(p: str) -> bool:
    p = os.path.basename(p)
    if p.startswith(".") or p.endswith("~"):
        return True
    return False


def path2str(path: str) -> str:
    head, fname = os.path.split(path)
    head, x = os.path.split(head)
    _, dirname = os.path.split(head)

    return "[" + dirname + "][" + fname + "]"


class AbinitTests:
    """
    High-level manager for the entire collection of ABINIT automatic tests.

    This class handles basic suite discovery, indexing, and high-level
    database construction and selection tasks.
    """

    def __init__(self):
        self.suite_names = tuple([os.path.basename(d) for d in _tsuite_dirs])
        self.suite_paths = tuple(
            [os.path.join(abenv.tests_dir, d) for d in _tsuite_dirs]
        )

        self._suites = dict()
        for suite_name, suite_path in self.walk_suites():
            self._suites[suite_name] = Suite(suite_path)

        # Check suite_names and subsuite_names
        all_subsuite_names = self.all_subsuite_names
        for suite_name, suite_path in self.walk_suites():
            if suite_name in all_subsuite_names:
                print("Found suite and subsuite with the same name: %s" % suite_name)

    def walk_suites(self) -> Iterator[tuple[str, str]]:
        """
        Iterator over all registered suites.

        Yields:
            tuple: (suite_name, suite_path)
        """
        return zip(self.suite_names, self.suite_paths)

    def __str__(self):
        return "\n".join([str(k) + " : " + str(v) for (k, v) in self.__dict__.items()])

    def get_suite(self, suite_name: str) -> Suite:
        """
        Get the `Suite` object for a given name.

        Args:
            suite_name: Name of the suite to retrieve.

        Returns:
            Suite: The suite object.
        """
        return self._suites[suite_name]

    @property
    def suites(self) -> Iterable[Suite]:
        return self._suites.values()

    def multi_parallel_suites(self) -> list[Suite]:
        """
        List of all suites containing multi-parallel tests.

        Returns:
            list[Suite]: List of suites where `is_multi_parallel` is True.
        """
        return [s for s in self.suites if s.is_multi_parallel]

    def suite_of_subsuite(self, subsuite_name: str) -> Suite:
        """
        Find the parent suite containing a specific subsuite.

        Args:
            subsuite_name: Name of the subsuite to search for.

        Returns:
            Suite: The parent suite object.

        Raises:
            ValueError: If the subsuite name is not registered.
        """
        for suite in self.suites:
            if suite.has_subsuite(subsuite_name):
                return suite

        raise ValueError("subsuite %s not found" % subsuite_name)

    @property
    def all_subsuite_names(self) -> list[str]:
        """List with the names of all the registered subsuites."""
        all_subnames = []
        for suite in self.suites:
            all_subnames.extend(suite.subsuites.keys())

        if len(all_subnames) != len(set(all_subnames)):
            raise RuntimeError(
                "The suite/subsuite name must be unique\n"
                "Please change the name of the suite/subsuite"
            )

        return all_subnames

    def keywords_of_suite(self, suite_name: str) -> set[str]:
        return self._suites[suite_name].keywords

    def cpp_vars_of_suite(self, suite_name: str) -> set[str]:
        return self._suites[suite_name].need_cpp_vars

    # def get_all_need_cppvars(self):
    #    all_need_cppvars = set()

    #    for suite_name in self.suite_names:
    #        cpp_vars = self.cpp_vars_of_suite(suite_name)
    #        if cpp_vars: print(cpp_vars)
    #        all_need_cppvars = all_need_cppvars.union(cpp_vars)

    #    database = self.build_database(with_disabled=False)
    #    for test in database:
    #        print(test)

    #    return all_need_cppvars

    def inputs_of_suite(self, suite_name: str, active: bool = True) -> list[str]:
        """
        Get the list of input files for a specific suite.

        Args:
            suite_name (str): Name of the suite.
            active (bool): If True, return active tests. Otherwise, return disabled ones.

        Returns:
            list: Paths to the input files.
        """
        if active:
            return self._suites[suite_name].inp_paths
        return self._suites[suite_name].disabled_inp_paths

    def build_database(self, with_disabled: bool = False) -> AbinitTestsDatabase:
        """
        Build a comprehensive tests database from the filesystem.

        Args:
            with_disabled: If True, include tests marked as disabled in the database.

        Returns:
            AbinitTestsDatabase: The fully populated database instance.
        """
        database = AbinitTestsDatabase(self._suites)

        for suite_name in self.suite_names:
            inp_files = self.inputs_of_suite(suite_name, active=True)
            if with_disabled:
                inp_files.extend(self.inputs_of_suite(suite_name, active=False))

            test_suite = AbinitTestSuite(
                abenv,
                inp_files=inp_files,
                keywords=self.keywords_of_suite(suite_name),
                need_cpp_vars=self.cpp_vars_of_suite(suite_name),
            )

            database.add_test_suite(suite_name, test_suite)

        return database

    def get_database(
        self, regenerate: bool = False, with_pickle: bool = False
    ) -> AbinitTestsDatabase:
        """
        Retrieve the tests database, optionally loading from a pickle cache.

        Args:
            regenerate: If True, force recalculation of the database instead of
                loading from cache.
            with_pickle: If True, save the database to a pickle file after generation.

        Returns:
            AbinitTestsDatabase: The tests database instance.
        """
        if regenerate or not os.path.exists(database_path):
            cprint("Regenerating database...", "yellow")
            database = self.build_database()

            # Save the database in the cpickle file.
            # Use file locking mechanism to prevent IO from other processes.
            if with_pickle:
                print("Saving database to %s" % database_path)
                with FileLock(database_path), open(database_path, "wb") as fh:
                    pickle.dump(database, fh, protocol=-1)

        else:
            cprint("Loading database from: %s" % database_path, "yellow")

            # Read the database from the cpickle file.
            # Use file locking mechanism to prevent IO from other processes.
            with FileLock(database_path), open(database_path, "rb") as fh:
                database = pickle.load(fh)

        return database

    def _suite_args_parser(
        self, args: list[str] | None = None
    ) -> dict[tuple[str, str | None], list[slice]]:
        """
        Parse script arguments. Return a mapping suite_name --> [slice objects]
        Three forms are possible
        0)                               Run all tests.
        1) v2[34:35] v1[12:] v5 v6[:45]  Select slices in the suites
        2) v4- v5-                       Exclude suites
        """

        # Mapping (suite_name, subsuite_name) --> slice_obj
        def all_tests() -> dict[tuple[str, str | None], list[slice]]:
            tuples = [(name, None) for name in self.suite_names]
            return dict.fromkeys(
                tuples,
                [
                    slice(None),
                ],
            )

        if args is None or not args:
            # Run all tests.
            return all_tests()

        args = [s.replace(" ", "") for s in args]

        exclude_mode = any(arg.endswith("-") for arg in args)
        if exclude_mode:
            d = all_tests()
            for arg in args:
                if arg.endswith("-"):
                    arg = arg[:-1]
                    if arg in self.suite_names:
                        d.pop((arg, None))  # v4- --> Remove v4
                    else:
                        # TODO
                        raise NotImplementedError(
                            "exclude_mode does not support subsuites"
                        )
                        suite = self.suite_of_subsuite(arg)
                        d.pop((suite.name, arg))  # gw1- --> skip tutorial/gw1

        else:
            re_slice = re.compile(r"^\[(\d*):(\d*)\]$")
            re_single = re.compile(r"^\[(\d*)\]$")

            d = {}
            for arg in args:
                start_stop = slice(None)
                idx = arg.find("[")
                if idx != -1:
                    arg, string = arg[:idx], arg[idx:]
                    match = re_slice.search(string)
                    if match:
                        start, stop = match.group(1), match.group(2)
                        # if not start: start = 1
                        if not start:
                            start = 0
                        if not stop:
                            stop = None
                        if stop is None:
                            start_stop = slice(int(start), stop)
                        else:
                            start_stop = slice(int(start), int(stop))
                    else:
                        match = re_single.search(string)
                        if match:
                            start = int(match.group(1))
                            start_stop = slice(start, start + 1)
                        else:
                            raise ValueError("Wrong or unknown argument: %s" % arg)

                if arg in self.suite_names:
                    tp = (arg, None)
                elif arg in self.all_subsuite_names:
                    suite = self.suite_of_subsuite(arg)
                    tp = (suite.name, arg)
                else:
                    raise ValueError(
                        "Wrong (suite_name|subsuite_name): `%s`. Did you remove the initial `t`?"
                        % arg
                    )

                if tp not in d:
                    d[tp] = [start_stop]
                else:
                    d[tp].append(start_stop)

        return d

    def select_tests(
        self,
        suite_args: list[str] | None,
        regenerate: bool = False,
        keys: list[str] | None = None,
        authors: list[str] | None = None,
        ivars: list[str] | None = None,
        with_pickle: bool = True,
        flat_list: bool = False,
    ) -> AbinitTestSuite | list[Any]:
        """
        Construct a test suite based on selection arguments and filters.

        Args:
            suite_args: Arguments specifying suites and slices (e.g., "v1[1:10]").
            regenerate: If True, regenerate the tests database.
            keys: List of keywords to filter by.
            authors: List of authors to filter by.
            ivars: List of input variables to filter by.
            with_pickle: If True, use the pickle cache for the database.
            flat_list: If True, return a flat list of tests instead of chained ones.

        Returns:
            AbinitTestSuite: The suite of selected tests.
        """
        tests_todo = self._suite_args_parser(suite_args)

        # Load the full database.
        database = self.get_database(regenerate=regenerate, with_pickle=with_pickle)

        # Extract the tests to run as specified by suite_args i.e by the string "v1[1:4] v3 ..."
        # TODO waiting for changes in the naming scheme
        # suites_without_slicing = ["tutoparal", "paral", "mpiio", "built-in", "seq"]
        # suites_without_slicing = ["tutoparal", "mpiio", "built-in", "seq"]
        # suites_without_slicing = ["tutoparal", "built-in",]
        suites_without_slicing = [
            "built-in",
        ]

        tests = AbinitTestSuite(abenv, test_list=[])

        # FIXME Not the sorting algorithm we would like to have!
        tuples = sorted(tests_todo.keys())

        for t in tuples:
            suite_name, subsuite_name = t
            for slice_obj in tests_todo[t]:
                # print("Extracting suite_name: %s, subsuite_name: %s, slice_obj: %s" % (suite_name, subsuite_name, slice_obj))

                # FIXME
                if suite_name in suites_without_slicing:
                    slice_obj = None

                tests = tests + database.get_test_suite(
                    suite_name, subsuite_name=subsuite_name, slice_obj=slice_obj
                )

        if keys or authors or ivars:
            # Create new suite whose tests contain the specified keywords.
            with_keys, exclude_keys, with_authors, exclude_authors = 4 * (None,)

            if keys:
                with_keys = [k for k in keys if not k.endswith("-")]
                exclude_keys = [k[:-1] for k in keys if k.endswith("-")]
                print(
                    "Extracting tests with keywords = %s, without keywords %s"
                    % (with_keys, exclude_keys)
                )

                if "VASP" in with_keys:
                    from pymods.tools import ascii_wasp

                    print(ascii_wasp())
                    print("Maybe you meant ABINIT!")
                    sys.exit(1)

                if "PGI" in with_keys:
                    from pymods.tools import ascii_scream

                    print(ascii_scream())
                    print(
                        "I really can't imagine how PGI could pass the ABINIT test suite!"
                    )
                    sys.exit(1)

            if authors:
                with_authors = [a for a in authors if not a.endswith("-")]
                exclude_authors = [a[:-1] for a in authors if a.endswith("-")]
                print(
                    "Extracting tests with authors = %s, without authors %s"
                    % (with_authors, exclude_authors)
                )

            tests = tests.select_tests(
                with_keys=with_keys,
                exclude_keys=exclude_keys,
                with_authors=with_authors,
                exclude_authors=exclude_authors,
                ivars=ivars,
            )
        if not flat_list:
            return tests

        # Build flat list of tests.
        flat = []
        for t in tests:
            if isinstance(t, ChainOfTests):
                # DO NOT use isinstance to check if ChainOfTests but rely on duck typing.
                # if hasattr(t, "tests"):
                flat.extend(t.tests)
            else:
                flat.append(t)
        return flat

    def generate_html_listoftests(self) -> None:
        """Generate the ListOfTests files"""
        database = self.get_database(regenerate=True)

        for suite_name, suite_path in self.walk_suites():
            suite = database.get_test_suite(suite_name)

            fname = os.path.join(suite_path, "ListOfTests.html")
            print("Writing ListOfTests HTML file: ", fname)
            with open(fname, "w") as fh:
                fh.write(suite.make_listoftests(width=160, html=True))

            fname = os.path.join(suite_path, "ListOfTests.txt")
            print("Writing ListOfTests text file: ", fname)
            with open(fname, "w") as fh:
                fh.write(suite.make_listoftests(width=100, html=False))

    def show_info(self, verbose: int = 0) -> None:
        """
        Print info on the test suite.
        """
        table = [["Suite", "# Activated Tests", "# Disabled Tests"]]
        for suite_name in self.suite_names:
            active_tests = self.inputs_of_suite(suite_name, active=True)
            disabled_tests = self.inputs_of_suite(suite_name, active=False)
            table.append([suite_name, str(len(active_tests)), str(len(disabled_tests))])
        from tests.pymods.tools import pprint_table

        print()
        pprint_table(table)
        print()

        print(8 * "=" + " KEYWORDS " + 8 * "=")
        width = max([len(k) for k in KNOWN_KEYWORDS]) + 5
        for skey in sorted(KNOWN_KEYWORDS.keys()):
            info = KNOWN_KEYWORDS[skey]
            print(skey.ljust(width), info)
        print()

        if verbose:
            for suite_name in self.suite_names:
                suite = self.get_suite(suite_name)
                if not suite.need_cpp_vars:
                    continue
                print("suite:", suite_name, "needs CPP variables:", suite.need_cpp_vars)

            # if verbose > 1:
            #    # list authors
            #    database = self.get_database(regenerate=True)
            #    pprint(database.authors_snames)

        # TODO: add this test to check_test_suite
        # chains = database.test_chains()
        # for c in chains:
        #   string, nlinks = c.info_on_chain()
        #   if nlinks == 0:
        #       print(15 * "*" + " Warning: found 0 explicit links " + 15 * "*")
        # print(string)

        print("\nUse verbose > 0 to print more info.")

        # This to print tests with keywords
        """
        database = self.get_database(regenerate=False)
        key = 'linear electro-optical coefficient'
        key = "UJDET"

        print("Reporting tests with key:", key)
        for suite_name, tests in database.items():
            #print(suite_name)
            for test in tests:
                if key in test.keywords:
                    try:
                        print(test, test.inp_fname)
                    except AttributeError:
                        print(test, [t.inp_fname for t in test])
        """


abitests = AbinitTests()
