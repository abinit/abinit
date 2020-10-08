from __future__ import print_function, division, absolute_import #, unicode_literals

import sys
import os
import imp
import platform
import re

from tests.pymods.termcolor import cprint

# Handle py2, py3k differences.
py2 = sys.version_info[0] <= 2
if py2:
    import cPickle as pickle
    from cStringIO import StringIO
else:
    import pickle
    from io import StringIO

from socket import gethostname
from pprint import pprint

from tests.pymods.testsuite import ChainOfTests, AbinitTestSuite
#from tests.pymods.tools import pprint_table
from tests.pymods.devtools import FileLock

import logging
logger = logging.getLogger(__name__)


__all__ = [
]


class AbinitEnvironment(object):
    """
    Container with information on the abinit source tree.
    Provide helper functions to construct the absolute path of directories.
    """
    def __init__(self):
        self.uname = platform.uname()
        self.hostname = gethostname()
        try:
            self.username = os.getlogin()
        except:
            self.username = "No_username"

        filepath = os.path.abspath(__file__)
        # Follows symlink (if any)
        filepath = os.path.realpath(filepath)

        # Paths of important dirs.
        self.tests_dir, tail = os.path.split(filepath)
        self.home_dir, tail = os.path.split(self.tests_dir)

        self.src_dir = os.path.join(self.home_dir, "src")
        self.psps_dir = os.path.join(self.tests_dir, "Psps_for_tests")
        self.fldiff_path = os.path.join(self.tests_dir, "Scripts", "fldiff.pl")

    def __str__(self):
        return "\n".join([str(k) + " : " + str(v) for (k, v) in self.__dict__.items()])

    def apath_of(self, *p):
        """
        Return the absolute path of p where p is one or more pathname components
        relative to the top level directory of the package.
        """
        return os.path.join(self.home_dir, *p)

    def isbuild(self):
        """True if we have built the code in the top level directory."""
        configh_path = os.path.join(self.home_dir, "config.h")
        abinit_path = os.path.join(self.home_dir, "src", "98_main", "abinit")
        return os.path.isfile(configh_path) and os.path.isfile(abinit_path)

    def touch_srcfiles(self, patterns):
        """
        Touch all Abinit source files containing the specified list of `patterns`.
        """
        def touch(fname):
            """
            Python touch
            See also http://stackoverflow.com/questions/1158076/implement-touch-using-python
            for a race-condition free version based on py3k
            """
            try:
                os.utime(fname, None)
            except:
                open(fname, 'a').close()

        top = os.path.join(abenv.home_dir, "src")
        print("Walking through directory:", top)
        print("Touching all files containing:", str(patterns))

        for root, dirs, files in os.walk(self.src_dir):
            for path in files:
                if not path.lower().endswith(".f90"): continue
                path = os.path.join(root, path)
                with open(path, "rt") as fh:
                    #print(path)
                    for line in fh:
                        if any(p in line for p in patterns):
                            print("Touching %s" % path)
                            touch(path)
                            break

    def start_watching_sources(self):
        def is_source(fname):
            # NB: `.finc` include files are ignored on purpose because
            # make is not tracking them. one should change the Fortran file
            # that includes .finc to force recompilation.
            _, ext = os.path.splitext(fname)
            return ext.lower() in (".f", ".f90", ".c") # ".finc"

        self.srcpath_stat = {}
        for root, dirs, files in os.walk(self.src_dir):
            for fname in files:
                if not is_source(fname): continue
                path = os.path.join(root, fname)
                self.srcpath_stat[path] = os.stat(path)

        return self.srcpath_stat

    def changed_sources(self):
        """Return list of files that have been modified."""
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
    #"abirules",
    "atompaw",
    "bigdft",
    "bigdft_paral",
    #"buildsys",
    "built-in",
    #"cpu",      This directory is disabled
    "etsf_io",
    "fast",
    "psml",
    "gpu",
    "libxc",
    "mpiio",
    "paral",
    #"hpc",
    #"physics",
    "seq",
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
    "vdwxc",
    "wannier90",
]

_tsuite_dirs.sort()
_tsuite_dirs = tuple([os.path.join(abenv.tests_dir, dir_name) for dir_name in _tsuite_dirs])


class Suite(object):
    """Information on one test suite"""
    def __init__(self, suite_path):

        suite_path = os.path.abspath(suite_path)

        self.suite_path = os.path.abspath(suite_path)
        self.name = os.path.basename(suite_path)

        module_name = os.path.join(suite_path, "__init__.py")
        module = imp.load_source(module_name, os.path.join(suite_path, "__init__.py") )

        self.keywords = set(module.keywords)
        self.need_cpp_vars = set(module.need_cpp_vars)

        # Divide tests into (active|disabled).
        self.inp_paths = [p for p in module.inp_files if not p.startswith("-")]
        self.disabled_inp_paths = [p[1:] for p in module.inp_files if p.startswith("-")]

        # Use absolute paths
        self.inp_paths = [os.path.join(suite_path, "Input", p) for p in self.inp_paths]
        self.disabled_inp_paths = [os.path.join(suite_path, "Input", p) for p in self.disabled_inp_paths]

        # True if the suite contains tests that should be executed with
        # different numbers of MPI processes
        self.is_multi_parallel = False
        if hasattr(module, "is_multi_parallel"):
            self.is_multi_parallel = module.is_multi_parallel

        self.subsuites = {}
        if hasattr(module, "subsuites"):
            subsuite_names = module.subsuites
            for k in subsuite_names: self.subsuites[k] = []

            for name in subsuite_names:
                pattern = re.compile("-?t" + name + "_\d+\.in")
                for inp in module.inp_files:
                    if pattern.match(inp):
                        #print(inp, "--> subsuite: ", name)
                        inp_path = os.path.join(suite_path, "Input", inp)
                        self.subsuites[name].append(inp_path)

            nfound = sum([len(paths) for paths in self.subsuites.values()])
            if nfound != len(module.inp_files):
                err_msg = ("At least one input_file does not belong to a subsuite, nfound = %s, __init__.nfiles = %s\n"
                           % (nfound, len(module.inp_files)))

                for inp in module.inp_files:
                    for paths in self.subsuites.values():
                        fnames = [os.path.basename(p) for p in paths]
                        if inp in fnames: break
                    else:
                        err_msg += "%s not found\n" % inp

                raise ValueError(err_msg)

            # Remove disabled tests (if any).
            for sub_name, paths in self.subsuites.items():
                inp_paths = [p for p in paths if not p.endswith("-")]
                self.subsuites[sub_name] = inp_paths

    def has_subsuite(self, subsuite_name):
        return subsuite_name in self.subsuites

    def inputs_of_subsuite(self, subsuite_name):
        """Return the absolut path of the input files in the subsuite."""
        return self.subsuites[subsuite_name]

    def __str__(self):
        return "\n".join([str(k) + " : " + str(v) for k, v in self.__dict__.items()])


class AbinitTestsDatabase(dict):
    """Database of tests indexed by the name of the abinit suite"""

    def __init__(self, suites):
        dict.__init__(self)
        self._suites = suites

    def iter_tests(self):
        """Iterate over all the tests."""
        for suite in self.values():
            for test in suite:
                yield test

    @property
    def suite_names(self):
        """List of suite names"""
        return [suite.name for suite in self._suites.values()]

    @property
    def authors_snames(self):
        """List of authors' second names extracted from the tests"""
        all_snames = []
        for test in self.iter_tests():
            all_snames.extend(test._authors_snames)

        return set(all_snames)

    def test_chains(self):
        """Return a list with all the chained tests."""
        return [t for t in self.iter_tests() if isinstance(t, ChainOfTests)]

    def tests_with_variables(self, ivars):
        """
        Return a list with all the tests that contain the input variables ivars
        """
        return [t for t in self.iter_tests() if t.has_variables(ivars)]

    def add_test_suite(self, suite_name, test_suite):
        """Add test_suite to the database using the key suite_name."""
        if suite_name not in self.suite_names:
            raise ValueError("%s is not a valid suite name" % suite_name)

        if suite_name in self:
            raise ValueError("%s is already in the database" % suite_name)

        self[suite_name] = test_suite

    def init_result_table(self):
        """
        Initialize a nested dictionary indexed by [suite_name][test.id].
        The dictionary contains ALL the tests present in the database
        and is used in testbot.py to store the results of the tests executed
        by the buildbot slave. Values are initialized with {}.
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
                    raise ValueError("Replicated test.id %s in suite %s" % (test.id, suite_name))

                res_table[suite_name][test.id] = {}

        return res_table

    def get_test_suite(self, suite_name, subsuite_name=None, slice_obj=None):
        """
        Return a list of tests belonging to the suite suite_name.
        If subsuite_name is not None, only the tests in suite_name/subsuite_name are returned.
        slice_obj is a slice object that can be used to specify the initial and final number of the test.
        if slice_obj is None, all the tests in suite_name/subsuite_name are returned.
        """
        test_suite = self[suite_name]

        if subsuite_name is not None:
            # Build the tests in the subsuite from the input files.
            suite = self._suites[suite_name]

            if not suite.has_subsuite(subsuite_name):
                raise ValueError("suite %s does not have subsuite %s" % (suite_name, subsuite_name))

            sub_inputs = suite.inputs_of_subsuite(subsuite_name)

            abenv = test_suite.abenv
            test_suite = AbinitTestSuite(abenv,
                                         inp_files=sub_inputs,
                                         keywords=suite.keywords,
                                         need_cpp_vars=suite.need_cpp_vars)

        if slice_obj is None:
            return test_suite
        else:
            logger.debug("will slice test_suite with slice_obj= %s " % slice_obj)
            return test_suite[slice_obj]

    def find_unknown_wrong_keywords(self):
        """
        Check whether TEST_INFO keywords have been documented.
        Returns set with the unknown keywords.
        """
        unknowns, wrong = set(), set()
        for suite_name, suite in self.items():
            for test in suite:
                for key in test.keywords:
                    #if not key: continue
                    if key not in KNOWN_KEYWORDS:
                        unknowns.add(key)
                    if " " in key:
                        wrong.add(key)

        return unknowns, wrong

    def find_stale_or_lost_inputs(self):
        """
        Check whether all reference files located in the Refs directories
        are tested namely that if they appear in the files_to_test field.

        Returns:
            String with the errors, empty string if OK.
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
                        raise ValueError("Input [%s] [%s] does not appear in Input2keys!" % (suite_name, ius))
                    else:
                        inp2test[ius] += 1

            def remove_file(fname):
                # XG130810 : When the report.in files in abirules/Input/report.in  buildsys/Input/report.in
                # will have been suppressed, one might replace the next line by the simpler :
                # return fname.endswith(".files")
                return fname.endswith(".files") or os.path.basename(fname) in ["report.in"]

            keys = []
            for fname, ntimes in inp2test.items():
                if not remove_file(fname):
                    keys.append(fname)
                else:
                    ntest = inp2test[fname]
                    assert ntest == 0
            #inp2test = {k: inp2test[k] for k in keys} # requires py2.7
            inp2test = dict([(k, inp2test[k]) for k in keys])

            # At this point inp2test should be >= 1.
            for fname, ntimes in inp2test.items():
                #if ntimes != 1:
                if ntimes == 0:
                    err.write("Input file %s is used %s time(s)\n" % (path2str(fname), ntimes))

        return err.getvalue()

    def find_stale_or_lost_refs(self):
        """
        Check whether all reference files located in the Refs directories
        are tested namely that if they compare in the files_to_test field.

        Returns:
            String with the errors, empty string if OK
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
            ref_fnames = [os.path.join(ref_dir, f) for f in listdir] # use absolute path.

            # Mapping ref_fname --> number of tests using it.
            ref2test = dict().fromkeys(ref_fnames, 0)

            for test in suite:
                files_to_test = [os.path.join(ref_dir, f.name) for f in test.files_to_test]

                for o in files_to_test:
                    # FIXME due to out --> stdout replacement
                    if o.endswith(".stdout"): o = o[:-7] + ".out"
                    if o not in ref2test:
                        err.write("files_to_test %s does not appear in Refs!\n" % o)
                    else:
                        ref2test[o] += 1

            # At this point ref2test should contain only ones.
            for ref_fname, ntimes in ref2test.items():
                if ntimes != 1:
                    err.write("Ref file %s is tested %s time(s)\n" % (path2str(ref_fname), ntimes))

        return err.getvalue()

    def check_testinfo_options(self):
        """
        Test the presence of important options in the TEST_INFO section of each test.
        """
        def check_options_in_test(test):
            recommended_opts = [
                "keywords",
                "description",
                "authors",
                "max_nprocs",]

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
                    #print("In test chain %s" % test.full_id)
                    for t in test:
                        wrong_options = check_options_in_test(t)
                        for opt, stat in wrong_options.items():
                            app("%s: option %s is %s" % (t.full_id, opt, stat))

        return "\n".join(lines)


def exclude_path(p):
    p = os.path.basename(p)
    if p.startswith(".") or p.endswith("~"): return True
    return False


def path2str(path):
    head, fname = os.path.split(path)
    head, x = os.path.split(head)
    _, dirname = os.path.split(head)

    return "["+dirname+"]["+fname+"]"


class AbinitTests(object):
    """
    Object describing the collection of automatic tests
    """
    def __init__(self):
        self.suite_names = tuple([os.path.basename(d) for d in _tsuite_dirs])
        self.suite_paths = tuple([os.path.join(abenv.tests_dir, d) for d in _tsuite_dirs])

        self._suites = dict()
        for (suite_name, suite_path) in self.walk_suites():
            self._suites[suite_name] = Suite(suite_path)

        # Check suite_names and subsuite_names
        all_subsuite_names = self.all_subsuite_names
        for suite_name, suite_path in self.walk_suites():
            if suite_name in all_subsuite_names:
                print("Found suite and subsuite with the same name: %s" % suite_name)

    def walk_suites(self):
        """return list of (suite_name, suite_paths)"""
        return zip(self.suite_names, self.suite_paths)

    def __str__(self):
        return "\n".join([str(k) + " : " + str(v) for (k,v) in self.__dict__.items()])

    def get_suite(self, suite_name):
        return self._suites[suite_name]

    @property
    def suites(self):
        return self._suites.values()

    def multi_parallel_suites(self):
        suites = self.suites
        return [s for s in suites if s.is_multi_parallel]

    def suite_of_subsuite(self, subsuite_name):
        """Return the suite corresponding to the given subsuite."""
        for suite in self.suites:
            if suite.has_subsuite(subsuite_name): return suite
        else:
            raise ValueError("subsuite %s not found" % subsuite_name)

    @property
    def all_subsuite_names(self):
        """List with the names of all the registered subsuites."""
        all_subnames = []
        for suite in self.suites:
            all_subnames.extend(suite.subsuites.keys())

        if len(all_subnames) != len(set(all_subnames)):
            raise RuntimeError("The suite/subsuite name must be unique\n" +
                "Please change the name of the suite/subsuite")

        return all_subnames

    def keywords_of_suite(self, suite_name):
        return self._suites[suite_name].keywords

    def cpp_vars_of_suite(self, suite_name):
        return self._suites[suite_name].need_cpp_vars

    def inputs_of_suite(self, suite_name, active=True):
        if active:
            return self._suites[suite_name].inp_paths
        else:
            return self._suites[suite_name].disabled_inp_paths

    def build_database(self, with_disabled=False):
        """
        Return an instance of TestsDatabase containing all the ABINIT automatic tests.
        with_disabled specifies whether disabled tests should be included in the database.
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
                need_cpp_vars=self.cpp_vars_of_suite(suite_name))

            database.add_test_suite(suite_name, test_suite)

        return database

    def get_database(self, regenerate=False, with_pickle=False):
        """
        Return an instance of TestsDatabase initialized from an external pickle file.

        Args:
            regenerate:
                True to force the regeneration of the database
                and the writing of a new pickle file.
            with_pickle:
                Save the generated database in pickle format.
        """
        if regenerate or not os.path.exists(database_path):
            cprint("Regenerating database...", "yellow")
            database = self.build_database()

            # Save the database in the cpickle file.
            # Use file locking mechanism to prevent IO from other processes.
            if with_pickle:
                print("Saving database to %s" % database_path)
                lock = FileLock(database_path)

                lock.acquire()

                with open(database_path, "wb") as fh:
                    pickle.dump(database, fh, protocol=-1)

                lock.release()

        else:
            cprint("Loading database from: %s" % database_path, "yellow")

            # Read the database from the cpickle file.
            # Use file locking mechanism to prevent IO from other processes.
            lock = FileLock(database_path)

            lock.acquire()

            with open(database_path, "rb") as fh:
                database = pickle.load(fh)

            lock.release()

        return database

    def _suite_args_parser(self, args=None):
        """
        Parse script arguments. Return a mapping suite_name --> [slice objects]
        Three forms are possible
        0)                               Run all tests.
        1) v2[34:35] v1[12:] v5 v6[:45]  Select slices in the suites
        2) v4- v5-                       Exclude suites
        """
        # Mapping (suite_name, subsuite_name) --> slice_obj
        def all_tests():
            tuples = [(name, None) for name in self.suite_names]
            return dict.fromkeys(tuples, [slice(None),])

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
                        raise NotImplementedError("exclude_mode does not support subsuites")
                        suite = self.suite_of_subsuite(arg)
                        d.pop((suite.name, arg))  # gw1- --> skip tutorial/gw1

        else:
            import re
            re_slice = re.compile("^\[(\d*):(\d*)\]$")
            re_single = re.compile("^\[(\d*)\]$")

            d = {}
            for arg in args:
                start_stop = slice(None)
                idx = arg.find("[")
                if idx != -1:
                    arg, string = arg[:idx], arg[idx:]
                    match = re_slice.search(string)
                    if match:
                        start, stop = match.group(1), match.group(2)
                        #if not start: start = 1
                        if not start: start = 0
                        if not stop: stop = None
                        if stop is None:
                            start_stop = slice(int(start), stop)
                        else:
                            start_stop = slice(int(start), int(stop))
                    else:
                        match = re_single.search(string)
                        if match:
                            start = int(match.group(1))
                            start_stop = slice(start, start+1)
                        else:
                            raise ValueError("Wrong or unknown argument: %s" % arg)

                if arg in self.suite_names:
                    tp = (arg, None)
                elif arg in self.all_subsuite_names:
                    suite = self.suite_of_subsuite(arg)
                    tp = (suite.name, arg)
                else:
                    raise ValueError("Wrong (suite_name|subsuite_name) : %s" % arg)

                if tp not in d:
                    d[tp] = [start_stop]
                else:
                    d[tp].append(start_stop)

        return d

    def select_tests(self, suite_args, regenerate=False, keys=None,
                     authors=None, ivars=None, with_pickle=True, flat_list=False):
        """
        Main entry point for client code. Return an instance of AbinitTestSuite

        Args:
            with_pickle: Save the generated database in pickle format.
            flat_list: True if a flat list of tests should be returned instead
                of a list that can contain ChainOfTests objects.
        """
        tests_todo = self._suite_args_parser(suite_args)

        # Load the full database.
        database = self.get_database(regenerate=regenerate, with_pickle=with_pickle)

        # Extract the tests to run as specified by suite_args i.e by the string "v1[1:4] v3 ..."
        # TODO waiting for changes in the naming scheme
        #suites_without_slicing = ["tutoparal", "paral", "mpiio", "built-in", "seq"]
        #suites_without_slicing = ["tutoparal", "mpiio", "built-in", "seq"]
        #suites_without_slicing = ["tutoparal", "built-in",]
        suites_without_slicing = ["built-in",]

        tests = AbinitTestSuite(abenv, test_list=[])

        # FIXME Not the sorting algorithm we would like to have!
        tuples = sorted(tests_todo.keys())

        for t in tuples:
            suite_name, subsuite_name = t
            for slice_obj in tests_todo[t]:
                #print("Extracting suite_name: %s, subsuite_name: %s, slice_obj: %s" % (suite_name, subsuite_name, slice_obj))

                # FIXME
                if suite_name in suites_without_slicing: slice_obj = None

                tests = tests + database.get_test_suite(suite_name, subsuite_name=subsuite_name, slice_obj=slice_obj)

        if keys or authors or ivars:
            # Create new suite whose tests contain the specified keywords.
            with_keys, exclude_keys, with_authors, exclude_authors = 4 * (None,)

            if keys:
                with_keys = [k for k in keys if not k.endswith("-")]
                exclude_keys = [k[:-1] for k in keys if k.endswith("-")]
                print("Extracting tests with keywords = %s, without keywords %s" % (with_keys, exclude_keys))

                if "VASP" in with_keys:
                    from pymods.tools import ascii_wasp
                    print(ascii_wasp())
                    print("Maybe you meant ABINIT!")
                    sys.exit(1)

                if "PGI" in with_keys:
                    from pymods.tools import ascii_scream
                    print(ascii_scream())
                    print("I really can't imagine how PGI could pass the ABINIT test suite!")
                    sys.exit(1)

            if authors:
                with_authors = [a for a in authors if not a.endswith("-")]
                exclude_authors = [a[:-1] for a in authors if a.endswith("-")]
                print("Extracting tests with authors = %s, without authors %s" % (with_authors, exclude_authors))

            tests = tests.select_tests(with_keys= with_keys,
                                       exclude_keys=exclude_keys,
                                       with_authors=with_authors,
                                       exclude_authors=exclude_authors,
                                       ivars=ivars
                                       )
        if not flat_list:
            return tests
        else:
            # Build flat list of tests.
            flat = []
            for t in tests:
                if isinstance(t, ChainOfTests):
                # DO NOT use isinstance to check if ChainOfTests but rely on duck typing.
                #if hasattr(t, "tests"):
                    flat.extend(t.tests)
                else:
                    flat.append(t)
            return flat

    def generate_html_listoftests(self):
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

    def show_info(self):
        table = [["Suite",  "# Activated Tests", "# Disabled Tests"]]
        for suite_name in self.suite_names:
            active_tests = self.inputs_of_suite(suite_name, active=True)
            disabled_tests = self.inputs_of_suite(suite_name, active=False)
            table.append([suite_name, str(len(active_tests)), str(len(disabled_tests))])

        #pprint_table(table)

        skeys = sorted(KNOWN_KEYWORDS.keys())
        print(8 * "=" + " KEYWORDS " + 8 * "=")
        width = max( [len(k) for k in KNOWN_KEYWORDS] ) + 5
        for skey in skeys:
            desc = KNOWN_KEYWORDS[skey]
            print(skey.ljust(width), desc)

        #for suite_name in self.suite_names:
        #  suite = self.get_suite(suite_name)
        #  print(suite_name, suite.need_cpp_vars)

        #for subsuite_name in self.all_subsuite_names:
        #  print(subsuite_name)

        #for suite in self.multi_parallel_suites():
        #  print(suite.name)
        #sys.exit(0)

        # list authors
        database = self.get_database(regenerate=True)
        pprint(database.authors_snames)

        # TODO: add this test to chec_test_suite
        #chains = database.test_chains()
        #for c in chains:
        #    string, nlinks = c.info_on_chain()
        #    if nlinks == 0:
        #        print(15 * "*" + " Warning: found 0 explicit links " + 15 * "*")
        #print(string)

abitests = AbinitTests()


# Dictionary mapping test keywords to string with human-readable description.
KNOWN_KEYWORDS = {
    # Main executables
    "abinit": "Test abinit code.",
    "anaddb": "Test anaddb code",
    "atompaw": "Atompaw tests",
    "band2eps": "Test band2eps code",
    "cut3d": "Test cut3d code",
    "mrgscr": "Test mrgscr code",
    "mrggkk": "Test mrggkk code",
    "mrgddb": "Test mrgddb code",
    'mrgdv': "Test mrgdv code",
    "optic": "Test optic code",
    "ujdet": "Test ujdet code",
    "aim": "Test aim code",
    "conducti": "Test conducti code",
    "fftprof": "Test fftprof code and low-level FFT routines",
    "wannier90": "Test the interface with Wannier90",
    "macroave": "Test macroave code",
    "bigdft": "The the interface with Bigdft",
    # Keywords describing the test.
    "NC": "Calculations with norm-conserving pseudos",
    "PAW": "Calculations with PAW datasets",
    "GW": "GW calculations",
    "GWGamma": "GW calculations with vertex corrections",
    "BSE": "Bethe-Salpeter calculations",
    "DDK": "DDK calculations",
    "EPH": "Electron-phonon calculations",
    "DFTU": "DFT+U calculations",
    "LEXX": "Local exact exchange calculations",
    "DMFT": "Dynamical mean-field theory",
    "TDDFT": "Time-dependent Density Functional Theory",
    "STRING": "String method",
    "CML": "Chemical Markup Language",
    "WanT": "Interface with Want code",
    "XML": "XML output (test prtxml)",
    "CIF": "Tests producing Crystallographic Information Files",
    "positron": "Positron calculations",
    "DFPT": "Density functional perturbation theory",
    "SOC": "Spin-orbit coupling",
    "EFG": "Electric field gradient",
    "IMAGES": "Parallelization over images",
    "PARTIAL_DOS": "Partial DOS",
    "RTA": "Relaxation-time approximation",
    "DOS": "electronic DOS calculations",
    "STS": "STS calculations",
    "CTQMC": "CTQMC method for DMFT",
    "ETSF_IO": "Tests using ETSF-IO library",
    "netcdf": "Tests using netcdf",
    "WVL": "Wavelets",
    "PIMD": "Path integral molecular dynamics",
    "HF": "Calculations including Hartree-Fock term in DFT",
    "psp8": "Tests using pseudos in PSP8 format",
    "MPI_FFT": "Unit tests for MPI-FFT",
    "PBE0": "PBE0 calculations",
    "cRPA": "RPA for correlated electrons.",
    "FAILS_IFMPI": "Tests failing if MPI is used",
    'NVT': "(N,V,T) enseble",
    'ELASTIC': "Calculations of elastic constants",
    'INTERNAL_STRAIN': "Calculations of internal strain",
    'DFT-D3(BJ)': "DFT-D3 dispersion correction",
    'DFT-D3': "DFT-D3 dispersion correction",
    '3-BODY_TERM': "DFT-D3 dispersion correction",
    'GWLS': "GW with the Sternheimer approach",
    'Projected_Wannier': "Projected Wannier functions",
    'VDW': "van der Wall interaction",
    'LOBSTER': "Interface with LOBSTER code",
    'RELAXATION': "Structural relaxations",
    'magnetic_constraint': "Tests employing magnetic constraints",
    "FOLD2BLOCH": "Fold2Bloch tests.",
}
