from __future__ import print_function, division, absolute_import  # , unicode_literals

import sys
import os
import time
import shutil
import textwrap
import platform
import tarfile
import re
import warnings
from base64 import b64encode

from socket import gethostname
from subprocess import Popen, PIPE
from multiprocessing import Process, Queue, Lock
from threading import Thread
from contextlib import contextmanager

# Handle py2, py3k differences.
py2 = sys.version_info[0] <= 2
if py2:
    import cPickle as pickle
    from StringIO import StringIO
    from ConfigParser import SafeConfigParser, NoOptionError, ParsingError as CPError
    from Queue import Empty as EmptyQueueError
else:
    import pickle
    from io import StringIO
    from configparser import ConfigParser, ParsingError as CPError
    from queue import Empty as EmptyQueueError

from collections import OrderedDict
from .jobrunner import TimeBomb
from .tools import (RestrictedShell, unzip, tail_file, pprint_table, Patcher, Editor)
from .xyaptu import xcopier
from .devtools import NoErrorFileLock, makeunique
from .memprof import AbimemFile
from .termcolor import cprint
from .fldiff import Differ as FlDiffer

import logging
logger = logging.getLogger(__name__)

__version__ = "0.5"
__author__ = "Matteo Giantomassi"

__all__ = [
    "BuildEnvironment",
    "AbinitTestSuite",
]

fldebug = 'FLDIFF_DEBUG' in os.environ and os.environ['FLDIFF_DEBUG']

_MY_NAME = os.path.basename(__file__)[:-3] + "-" + __version__


# Helper functions and tools


def my_getlogin():
    """
    Returns the user logged in to the controlling terminal of the process.

    https://stackoverflow.com/questions/4399617/python-os-getlogin-problem
    """
    username = "No_username"
    if hasattr(os, 'getlogin'):
        try:
            username = os.getlogin()
        except Exception: # FileNotFoundError
            try:
                import pwd
                getlogin = lambda: pwd.getpwuid(os.getuid())[0]
                username = getlogin()
            except Exception:
                username = "No_username_tried_pwd"

    return username


@makeunique
def genid():
    '''
    Produce a random sequence 12 bytes represented as 16 ascii characters.
    The decorator ensure that output is different at each call.
    '''
    return b64encode(os.urandom(12)).decode('ascii')


def html_colorize_text(string, code):
    return "<FONT COLOR='%s'>%s</FONT>" % (code, string)


_status2htmlcolor = {
    "succeeded": lambda string: html_colorize_text(string, 'Green'),
    "passed": lambda string: html_colorize_text(string, 'DeepSkyBlue'),
    "failed": lambda string: html_colorize_text(string, 'Red'),
    "disabled": lambda string: html_colorize_text(string, 'Cyan'),
    "skipped": lambda string: html_colorize_text(string, 'Cyan'),
}


def status2html(status):
    """Convert test status in a colored HTML string."""
    return _status2htmlcolor[status](status)


def sec2str(seconds):
    """Convert seconds to string."""
    return "%.2f" % seconds


def str2html(string, end="<br>"):
    """Returns a HTML string."""
    lines = string.splitlines()
    return "<br>".join(lines) + end


def args2htmltr(*args):
    string = ""
    for arg in args:
        string += "<td>" + str(arg) + "</td>"
    return string


def html_link(string, href=None):
    """Create a HTML link from a string. Use href as link of href is not None."""
    if href is not None:
        return "<a href='%s'>%s</a>" % (href, string)
    else:
        return "<a href='%s'>%s</a>" % (string, string)


def is_string(s):
    """True is s is a string (duck typying test)"""
    try:
        s + "hello"
        return True
    except TypeError:
        return False


def has_exts(path, exts):
    """True if path ends with extensions exts"""
    root, ext = os.path.splitext(path)
    if is_string(exts):
        return ext == exts
    else:
        return ext in exts


def lazy__str__(func):
    """Lazy decorator for __str__ methods"""
    def oncall(*args, **kwargs):
        self = args[0]
        return "\n".join(str(k) + " : " + str(v) for (k, v) in self.__dict__.items())
    return oncall


# Helper functions for performing IO

def lazy_read(fname):
    if not py2:
        with open(fname, "rt", encoding="utf-8") as fh:
            return fh.read()

    else:
        with open(fname, "rt") as fh:
            return fh.read()


def lazy_readlines(fname):
    if not py2:
        with open(fname, "rt", encoding="utf-8") as fh:
            return fh.readlines()
    else:
        with open(fname, "rt") as fh:
            return fh.readlines()


def rmrf(top, exclude_paths=None):
    """
    Recursively remove all files and directories contained in directory top.

    Args:
        exclude_paths:
            list with the absolute paths that should be preserved

    Returns the list of files and the directories that have been removed.
    """
    exc_paths = []
    if exclude_paths is not None:
        if is_string(exclude_paths):
            exc_paths = [exclude_paths]
        else:
            exc_paths = exclude_paths

    removed = []
    for (root, dirs, files) in os.walk(top):
        for f in files:
            file_path = os.path.join(root, f)
            if file_path not in exc_paths:
                os.unlink(file_path)
                removed.append(file_path)
        for d in dirs:
            dir_path = os.path.join(root, d)
            if dir_path not in exc_paths:
                shutil.rmtree(dir_path)
                removed.append(dir_path)

    return removed


def find_abortfile(workdir):
    """
    Return the absolute path of the MPIABORTFILE file produced by (Abinit|Abinit_with_libpaw)
    Empty string if file is not present.

    Args:
        workdir: Working directory of the test.

    .. Note::

        __LIBPAW_MPIABORFILE__ is produced if abinit uses libpaw and when we die inside libpaw.
    """
    for s in ("__ABI_MPIABORTFILE__", "__LIBPAW_MPIABORFILE__"):
        path = os.path.join(workdir, s)
        if os.path.exists(s):
            return path
    return ""


def read_yaml_errmsg(path):
    """
    Extract the YAML error message from file `path`.
    Returns string with message, empty string if message is not found.

    The Yaml error message is in the form:

    --- !ERROR
    src_file: m_io_screening.F90
    src_line: 648
    message: |
        Unsupported value of iomode
    ...
    """
    errlines, inerr = [], 0

    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("---") and ("ERROR" in line or "BUG" in line):
                inerr = 1

            if inerr:
                errlines.append(line)
                if line.startswith("..."):
                    break

    return "".join(errlines)


def extract_errinfo_from_files(workdir):
    """
    Extract information from the files produced by the code when we run tests in debug mode.

    Return:
        String with the content of the files. Empty string if no debug file is found.
    """
    registered_exts = {".flun", ".mocc"}
    errinfo = []

    for path in os.listdir(workdir):
        _, ext = os.path.splitext(path)
        if ext not in registered_exts:
            continue
        with open(os.path.join(workdir, path), "rt") as fh:
            errinfo.append(" ")
            errinfo.append("From file: %s " % path)
            errinfo.extend(l.strip() for l in fh)

    return "\n".join(errinfo)


class FileToTest(object):
    """This object contains information on the output file that will be analyzed by fldiff"""
    #  atr_name,   default, conversion function. None designes mandatory attributes.
    _attrbs = [
        ("name", None, str),
        ("tolnlines", 0, int),    # fldiff tolerances
        ("tolabs", 0, float),
        ("tolrel", 0, float),
        ("fld_options", "", str),     # options passed to fldiff.
        ("fldiff_fname", "", str),
        ("hdiff_fname", "", str),
        ("diff_fname", "", str),
        ("use_yaml", "no", str),
        ("verbose_report", "no", str),
        # ("pydiff_fname","",str),
    ]

    def __init__(self, dic):

        for atr in FileToTest._attrbs:
            atr_name = atr[0]
            default = atr[1]
            f = atr[2]
            value = dic.get(atr_name, default)
            if value is None:
                raise ValueError("%s must be defined" % atr_name)

            value = f(value)
            if hasattr(value, "strip"):
                value = value.strip()
            self.__dict__[atr_name] = value

        # Postprocess fld_options
        self.fld_options = self.fld_options.split()
        for opt in self.fld_options:
            if not opt.startswith("-"):
                raise ValueError("Wrong fldiff option: %s" % opt)

        self.has_line_count_error = False
        self.do_html_diff = False

    @lazy__str__
    def __str__(self): pass

    def compare(self, fldiff_path, ref_dir, workdir, yaml_test, timebomb=None,
                outf=sys.stdout):
        """
        Use fldiff_path to compare the reference file located in ref_dir with
        the output file located in workdir. Results are written to stream outf.
        """
        ref_fname = os.path.abspath(os.path.join(ref_dir, self.name))
        # FIXME Hack due to the stdout-out ambiguity
        if not os.path.exists(ref_fname) and ref_fname.endswith(".stdout"):
            ref_fname = ref_fname[:-7] + ".out"
        out_fname = os.path.abspath(os.path.join(workdir, self.name))

        opts = {
            'label': self.name,
            'ignore': True,
            'ignoreP': True,
            'debug': fldebug,
        }

        if '-medium' in self.fld_options:
            opts['tolerance'] = 1.01e-8
        elif '-easy' in self.fld_options:
            opts['tolerance'] = 1.01e-5
        elif '-ridiculous' in self.fld_options:
            opts['tolerance'] = 1.01e-2

        if '-include' in self.fld_options:
            opts['ignore'] = False

        if '-includeP' in self.fld_options:
            opts['ignoreP'] = False

        if self.verbose_report == 'yes':
            opts['verbose'] = True

        if self.use_yaml not in ('yes', 'no', 'only'):
            # raise ParameterError
            pass

        if self.use_yaml == 'yes':
            opts['use_yaml'] = True
            opts['use_fl'] = True
        elif self.use_yaml == 'only':
            opts['use_yaml'] = True
            opts['use_fl'] = False
        else:  # self.use_yaml == 'no':
            opts['use_yaml'] = False
            opts['use_fl'] = True

        differ = FlDiffer(yaml_test=yaml_test, **opts)

        def make_diff():
            result = differ.diff(ref_fname, out_fname)
            result.dump_details(outf)

            return (result.passed_within_tols(
                self.tolnlines, self.tolabs, self.tolrel
            ), result.has_line_count_error())

        if fldebug:
            # fail on first error and output the traceback
            (isok, status, msg), has_line_count_error = make_diff()
        else:
            try:
                (isok, status, msg), has_line_count_error = make_diff()
            except Exception as e:
                warnings.warn(('[{}] Something went wrong with this test:\n'
                               '{}: {}\n').format(self.name, type(e).__name__,
                                                  str(e)))
                isok, status = False, 'failed'
                msg = 'internal error:\n{}: {}'.format(type(e).__name__,
                                                       str(e))
                has_line_count_error = False
        msg += ' [file={}]'.format(os.path.basename(ref_fname))

        # Save comparison results.
        self.fld_isok = isok
        self.fld_status = status
        self.fld_msg = msg
        self.has_line_count_error = has_line_count_error

        return isok, status, msg

# Parsers used for the different TEST_INFO options


def _str2filestotest(string):
    """
    Parse the files_to_test section.
    Returns a tuple of `FileToTest` objects.
    """
    if not string:
        return []

    if ";" in string:
        file_lines = [s for s in string.split(";") if s.strip()]
    else:
        file_lines = [string]

    files_to_test = []
    for line in file_lines:
        tokens = line.split(",")
        d = {"name": tokens[0]}
        for tok in tokens[1:]:
            k, v = [s.strip() for s in tok.split("=")]
            if k in d:
                err_msg = "Found multiple occurences of keyword %s" % k
                raise AbinitTestInfoParserError(err_msg)
            d[k] = v
        files_to_test.append(FileToTest(d))

    return tuple(files_to_test)


def _str2list(string):    return [s.strip() for s in string.split(",") if s]
def _str2intlist(string): return [int(item) for item in _str2list(string)]
def _str2set(string):     return {s.strip() for s in string.split(",") if s}
def _str2cmds(string):    return [s.strip() for s in string.split(";") if s]


def _str2bool(string):
    string = string.strip().lower()
    return string == "yes"


# TEST_INFO specifications
TESTCNF_KEYWORDS = {
    # keyword        : (parser, default, section, description)
    # [setup]
    "executable"     : (str       , None , "setup", "Name of the executable e.g. abinit"),
    #"use_files_file" : (_str2bool , "yes" , "setup", "Pass files file to executable (legacy mode)"),
    "use_files_file" : (_str2bool , "no" , "setup", "Pass files file to executable (legacy mode)"),
    "exec_args"      : (str       , ""   , "setup", "Arguments passed to executable on the command line."),
    "test_chain"     : (_str2list , ""   , "setup", "Defines a ChainOfTest i.e. a list of tests that are connected together."),
    "need_cpp_vars"  : (_str2set  , ""   , "setup", "CPP variables that must be defined in config.h in order to enable the test."),
    "exclude_hosts"  : (_str2list , ""   , "setup", "The test is not executed if we are running on a slave that matches compiler@hostname"),
    "exclude_builders": (_str2list, ""   , "setup", "The test is not executed if we are using a builder whose name is in the list"),
    "input_prefix"   : (str       , ""   , "setup", "Prefix for input files (used for the ABINIT files file)"),
    "output_prefix"  : (str       , ""   , "setup", "Prefix for output files (used for the ABINIT files file)"),
    "expected_failure": (_str2bool, "no" , "setup", "yes if the subprocess executing executable is expected to fail (retcode != 0) (default: no)"),
    "input_ddb"      : (str       , ""   , "setup", "The input DDB file read by anaddb"),
    "input_gkk"      : (str       , ""   , "setup", "The input GKK file read by anaddb"),
    "system_xml"     : (str       , ""   , "setup","The system.xml file read by multibinit"),
    "coeff_xml"      : (str       , ""   , "setup","The coeff.xml file read by multibinit"),
    "md_hist"        : (str       , ""   , "setup","The hist file file read by multibinit"),
    "no_check"       : (_str2bool , "no" , "setup","Explicitly do not check any files"),
    # [files]
    "files_to_test"  : (_str2filestotest, "", "files", "List with the output files that are be compared with the reference results. Format:\n" +
                                                       "\t file_name, tolnlines = int, tolabs = float, tolrel = float [,fld_options = -medium]\n" +
                                                       "\t    tolnlines: the tolerance on the number of differing lines\n" +
                                                       "\t    tolabs:the tolerance on the absolute error\n" +
                                                       "\t    tolrel: tolerance on the relative error\n" +
                                                       "\t    fld_options: options passed to fldiff.pl (optional).\n" +
                                                       "\t    Multiple files are separated by ; e.g.\n" +
                                                       "\t    foo.out, tolnlines = 2, tolabs = 0.1, tolrel = 1.0e-01;\n" +
                                                       "\t    bar.out, tolnlines = 4, tolabs = 0.0, tolrel = 1.0e-01"
                                                           ),
    "psp_files"      : (_str2list,        "", "files", "List of pseudopotential files (located in the Psps_for_tests directory)."),
    "extra_inputs"   : (_str2list,        "", "files", "List of extra input files."),
    # [shell]
    "pre_commands"   : (_str2cmds, "", "shell", "List of commands to execute before starting the test"),
    "post_commands"  : (_str2cmds, "", "shell", "List of commands to execute after the test is completed"),
    # [paral_info]
    "max_nprocs"     : (int      ,  1 , "paral_info", "Maximum number of MPI processors (1 for sequential run)"),
    "nprocs_to_test" : (_str2intlist, "","paral_info","List with the number of MPI processes that should be used for the test"),
    "exclude_nprocs" : (_str2intlist, "","paral_info","List with the number of MPI processes that should not be used for the test"),
    # [extra_info]
    "authors"         : (_str2set , "Unknown"                 , "extra_info", "Author(s) of the test"),
    "keywords"       : (_str2set , ""                         , "extra_info", "List of keywords associated to the test"),
    "description"    : (str      , "No description available",  "extra_info", "String containing extra information on the test"),
    "topics"         : (_str2list, "",  "extra_info", "Topics associated to the test"),
    "references"     : (_str2list, "",  "extra_info", "List of references to papers or other articles"),
    "file"           : (str, "", "yaml_test", "File path to the YAML config file relative to the input file."),
    "yaml"           : (str, "", "yaml_test", "Raw YAML config for quick config."),
}

# TESTCNF_SECTIONS = set( [ TESTCNF_KEYWORDS[k][2] for k in TESTCNF_KEYWORDS ] )

# This extra list is hardcoded in order to have a fixed order of the sections in doc_testcfn_format.
# OrderedDict have been introduced in python2.7 sigh!
TESTCNF_SECTIONS = {
    "setup",
    "files",
    "shell",
    "paral_info",
    "extra_info",
    "yaml_test",
}

# consistency check.
for key, tup in TESTCNF_KEYWORDS.items():
    if tup[2] not in TESTCNF_SECTIONS:
        raise ValueError("Please add the new section %s to TESTCNF_SECTIONS" % tup[2])


def line_starts_with_section_or_option(string):
    """True if string start with a TEST_INFO section or option."""
    from re import compile
    re_ncpu = compile(r"^NCPU_(\d+)$")
    s = string.strip()
    idx = s.find("=")
    if idx == -1:  # might be a section.
        if s.startswith("[") and s.endswith("]"):
            if s[1:-1] in TESTCNF_SECTIONS:
                return 1  # [files]...
            if re_ncpu.search(s[1:-1]):
                return 1    # [NCPU_1] ...
    else:
        if s[:idx].strip() in TESTCNF_KEYWORDS:
            return 2

    return 0


def doc_testcnf_format(fh=sys.stdout):
    """Automatic documentation of the TEST_INFO sections and related options."""
    def writen(string):
        fh.write(string + "\n")

    writen("Automatic documentation of the TEST_INFO sections and options.")

    for section in TESTCNF_SECTIONS:
        writen("\n[" + section + "]")
        for key in TESTCNF_KEYWORDS:
            tup = TESTCNF_KEYWORDS[key]
            if section == tup[2]:
                # line_parser = tup[0]
                default = tup[1]
                if default is None:
                    default = "Mandatory"
                desc = tup[3]
                if default:
                    msg = "%s =  %s (DEFAULT: %s)" % (key, desc, default)
                else:
                    msg = "%s =  %s" % (key, desc)
                writen(msg)


class AbinitTestInfo(object):
    """Container storing the options specified in the TEST_INFO section."""
    def __init__(self, dct):
        for k, v in dct.items():
            self.__dict__[k] = v

        # if self.nprocs_to_test and self.test_chain:
        #     raise TestInfoParserError("test_chain and nprocs_to_test are mutually exclusive")

        # Add the executable name to the list of keywords.
        self.add_keywords([self.executable])

    @lazy__str__
    def __str__(self): pass

    def add_cpp_vars(self, need_cpp_vars):
        """Add new set of CPP variables."""
        self.need_cpp_vars = self.need_cpp_vars.union(need_cpp_vars)

    def add_keywords(self, keywords):
        """Add new set of keywords."""
        self.keywords = self.keywords.union(keywords)

    def make_test_id(self):
        """
        Generate the string with the test identifier
        A special treatment is used for the multi-parallel tests.
        In this case, the test_id is constructed by appending the string _MPI#
        where # is the number of MPI processors.
        """
        # FIXME Assumes inp_fname is in the form name.in
        test_id = os.path.basename(self.inp_fname).split(".")[0]
        if self.ismulti_parallel:
            test_id += "_MPI%d" % self.max_nprocs
        return test_id

    @property
    def ismulti_parallel(self):
        """True is this is a multi-parallel test."""
        return self._ismulti_paral


class AbinitTestInfoParserError(Exception):
    """Error class raised by the parse"""


class AbinitTestInfoParser(object):
    """This object parses the TEST_INFO section that describes the test."""
    Error = AbinitTestInfoParserError

    def __init__(self, inp_fname, defaults=None):
        """
        Args:
            inp_fname: test input file
            defaults: default values passed to the INI parser.
        """
        logger.info("Parsing TEST_INFO section from input file : " + str(inp_fname))

        self.inp_fname = os.path.abspath(inp_fname)
        self.inp_dir, x = os.path.split(self.inp_fname)

        SENTINEL = '#%%'
        HEADER = "<BEGIN TEST_INFO>\n"
        FOOTER = "<END TEST_INFO>\n"

        lines = lazy_readlines(inp_fname)
        lines = [l.replace(SENTINEL, "", 1)
                 for l in lines if l.startswith(SENTINEL)]

        try:
            start, stop = lines.index(HEADER), lines.index(FOOTER)
        except ValueError:
            raise self.Error("{} does not contain any valid testcnf section!"
                             .format(inp_fname))

        # Keep only test section lines and remove one space at the begining
        lines = [line[1:] if line.startswith(' ') else line
                 for i, line in enumerate(lines) if start < i < stop]

        if not lines:
            raise self.Error("%s does not contain any valid testcnf section!" % inp_fname)

        # Interface in python 3 is richer so we rebuilt part of it
        if py2:
            class MySafeConfigParser(SafeConfigParser):
                """Wrap the get method of SafeConfigParser to disable the interpolation of raw_options."""
                raw_options = {"description"}

                def get(self, section, option, raw=False, vars=None):
                    if option in self.raw_options and section == TESTCNF_KEYWORDS[option][2]:
                        logger.debug("Disabling interpolation for section = %s, option = %s" % (section, option))
                        return SafeConfigParser.get(self, section, option, raw=True, vars=vars)
                    else:
                        return SafeConfigParser.get(self, section, option, raw, vars)

                def read_string(self, string, source='<string>'):
                    s = StringIO(string)
                    SafeConfigParser.readfp(self, s, filename=source)

            self.parser = MySafeConfigParser(defaults)
        else:
            self.parser = ConfigParser(defaults, interpolation=None)

        try:
            self.parser.read_string("".join(lines), source=inp_fname)
        except CPError as exc:
            cprint("Exception while parsing: %s\n%s" % (inp_fname, exc), "red")
            for l in lines:
                print(l, end="")
            cprint("A common problem is inappropriate indentation. The rules is"
                   ": do not indent options more than section title, indent "
                   "lines that belong to the option above.", "red")
            raise exc

        # Consistency check
        opt = "test_chain"
        section = TESTCNF_KEYWORDS[opt][2]
        pars = TESTCNF_KEYWORDS[opt][0]

        if self.parser.has_option(section, opt):
            string = self.parser.get(section, opt)
            chain = pars(string)
            ones = [chain.count(value) for value in chain]
            if sum(ones) != len(ones):
                err_msg = "%s : test_chain contains repeated tests %s" % (inp_fname, string)
                raise self.Error(err_msg)

    def generate_testinfo_nprocs(self, nprocs):
        """Returns a record with the variables needed to handle the job with nprocs."""
        d = {}

        d['yaml_test'] = self.yaml_test()

        # First read and parse the global options.
        for key in TESTCNF_KEYWORDS:
            tup = TESTCNF_KEYWORDS[key]
            line_parser = tup[0]
            section = tup[2]

            if section == 'yaml_test':
                # special case: handle this separatly
                continue
            elif section in self.parser.sections() and self.parser.has_option(
                section, key
            ):
                d[key] = self.parser.get(section, key)
            else:
                d[key] = tup[1]  # Section does not exist. Use default value.

            # Process the line
            try:
                d[key] = line_parser(d[key])
            except Exception as exc:
                err_msg = ("In file: %s\nWrong line:\n key = %s, d[key] = %s\n"
                           "%s: %s") % (self.inp_fname, key, d[key],
                                        type(exc).__name__, str(exc))
                raise self.Error(err_msg)

        # At this point info contains the parsed global values.
        # Now check if this is a parallel test and, in case, overwrite the values
        # using those reported in the [CPU_nprocs] sections.
        # Set also the value of info._ismulti_paral so that we know how to create the test id
        if not d['nprocs_to_test']:
            assert nprocs == 1
            d['_ismulti_paral'] = False
        else:
            logger.debug("multi parallel case")
            if nprocs not in d['nprocs_to_test']:
                err_msg = "in file: %s. nprocs = %s > not in nprocs_to_test = %s" % (self.inp_fname, nprocs, d['nprocs_to_test'])
                raise self.Error(err_msg)

            if nprocs > d['max_nprocs']:
                if hasattr(self, 'max_nprocs'):
                    err_msg = "in file: %s. nprocs = %s > max_nprocs = %s" % (self.inp_fname, nprocs, self.max_nprocs)
                else:
                    err_msg = "in file: %s\nmax_nprocs is not defined" % self.inp_fname

                raise self.Error(err_msg)

            # Redefine variables related to the number of CPUs.
            d['_ismulti_paral'] = True
            d['nprocs_to_test'] = [nprocs]
            d['max_nprocs'] = nprocs

            d['exclude_nprocs'] = list(range(1, nprocs))
            # print(self.inp_fname, nprocs, d['exclude_nprocs'])

            ncpu_section = "NCPU_" + str(nprocs)
            if not self.parser.has_section(ncpu_section):
                err_msg = "Cannot find section %s in %s" % (ncpu_section, self.inp_fname)
                raise self.Error(err_msg)

            for key in self.parser.options(ncpu_section):
                if key in self.parser.defaults():
                    continue
                opt = self.parser.get(ncpu_section, key)
                tup = TESTCNF_KEYWORDS[key]
                line_parser = tup[0]

                # Process the line and replace the global value.
                try:
                    d[key] = line_parser(opt)
                except Exception as exc:
                    err_msg = ("In file: %s\nWrong line:\n"
                               " key = %s, d[key] = %s\n %s: %s") % (
                                   self.inp_fname, key, d[key],
                                   type(exc).__name__, str(exc)
                    )
                    raise self.Error(err_msg)

                # print(self.inp_fname, d["max_nprocs"])

        # Add the name of the input file.
        d['inp_fname'] = self.inp_fname

        return AbinitTestInfo(d)

    @property
    def nprocs_to_test(self):
        """List with the number of MPI processors to be tested."""
        key = "nprocs_to_test"
        opt_parser = TESTCNF_KEYWORDS[key][0]
        default = TESTCNF_KEYWORDS[key][1]
        section = TESTCNF_KEYWORDS[key][2]

        if self.parser.has_option(section, key):
            opt = self.parser.get(section, key)
        else:
            opt = default

        return opt_parser(opt)

    @property
    def is_testchain(self):
        """True if this is a chain of tests"""
        opt = "test_chain"
        section = TESTCNF_KEYWORDS[opt][2]
        return self.parser.has_option(section, opt)

    def chain_inputs(self):
        """Return a list with the path of the input files belonging to the test chain"""
        assert self.is_testchain
        opt = "test_chain"
        section = TESTCNF_KEYWORDS[opt][2]
        parse = TESTCNF_KEYWORDS[opt][0]

        fnames = parse(self.parser.get(section, opt))
        return [os.path.join(self.inp_dir, fname) for fname in fnames]

    def yaml_test(self):
        sec_name = 'yaml_test'
        ytest = {}

        if self.parser.has_section(sec_name):
            scalar_key = ['file', 'yaml']
            for key in scalar_key:
                if self.parser.has_option(sec_name, key):
                    ytest[key] = self.parser.get(sec_name, key)

            if 'file' in ytest:
                val = ytest['file']
                base = os.path.realpath(os.path.dirname(self.inp_fname))
                ytest['file'] = os.path.join(base, val)

        return ytest


def find_top_build_tree(start_path, with_abinit=True, ntrials=10):
    """
    Returns the absolute path of the ABINIT build tree.
    Assume start_path is within the build tree.

    Raises:
        `RuntimeError` if build tree is not found after ntrials attempts.
    """
    abs_path = os.path.abspath(start_path)

    for _ in range(ntrials):
        config_h = os.path.join(abs_path, "config.h")
        abinit_bin = os.path.join(abs_path, "src", "98_main", "abinit")
        # Check if we are in the top of the ABINIT source tree
        if with_abinit:
            found = os.path.isfile(config_h) and os.path.isfile(abinit_bin)
        else:
            found = os.path.isfile(config_h)

        if found:
            return abs_path
        else:
            abs_path, _ = os.path.split(abs_path)

    raise RuntimeError("Cannot find the ABINIT build tree after %s trials" % ntrials)


class Compiler(object):
    """
    Base class for C,Fortran,C++ compilers.
    Usually instantiated through the class method from_defined_cpp_vars.
    """
    def __init__(self, name, version=None):
        self.name = name
        self.version = version

    def __str__(self):
        return "%s: %s %s" % (type(self).__name__, self.name, self.version)

    @classmethod
    def from_defined_cpp_vars(cls, defined_cpp_vars):
        for var in defined_cpp_vars:
            # TODO: version may be useful but it's not reported in config.h
            if var in cls._KNOWN_CPP_VARS:
                # Build the name of the compiler.
                name = var.lower().split("_")[1]
                if name == "gnu":
                    name = "gfortran"
                if name == "pathscale":
                    name = "psc"
                return cls(name=name, version=None)
        else:
            err_msg = "Cannot detect the name of the %s\n. Defined CPP vars: %s " % (cls.__name__, str(defined_cpp_vars))
            raise RuntimeError(err_msg)


class FortranCompiler(Compiler):
    """
    Store information on the Fortran compiler used to build abinit.
    """
    # CPP variables used in config.h
    _KNOWN_CPP_VARS = [
        "FC_ABSOFT",
        "FC_FUJITSU",
        "FC_G95",
        "FC_GNU",
        "FC_HITACHI",
        "FC_IBM",
        "FC_INTEL",
        "FC_MIPSPRO",
        "FC_NAG",
        "FC_OPEN64",
        "FC_PATHSCALE",
        "FC_PGI",
        "FC_SUN",
    ]


class CPreProcessorError(Exception):
    """Errors raised by `CPreProcessors`"""


class CPreProcessor(object):
    """Pre-process source code with ANSI CPP."""
    Error = CPreProcessorError

    def __init__(self, includes=None, opts=None, bin="cpp", verbose=0):
        self.includes = ["."]
        if includes is not None:
            self.includes = includes
        self.opts = ["-DHAVE_CONFIG_H"]
        if opts is not None:
            self.opts = opts
        self.bin, self.verbose = bin, verbose

    def process_file(self, filepath, remove_lhash=True):
        """
        Read source from filepath, call CPP wit the includes and the
        options passed to the constructor.

        Returns:
            preprocessed text.
        """
        if self.bin is None:
            # No pre-processing, return raw string.
            with open(filepath, "r") as f:
                return f.read()

        cmd = [self.bin]
        if self.opts:
            cmd += self.opts
        cmd += ["-ansi"]
        if self.includes:
            cmd += ["-I" + inc for inc in self.includes]
        cmd += [filepath]
        cmd = " ".join(cmd)
        if self.verbose:
            print(cmd)

        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()

        if p.returncode:
            raise self.Error("C-preprocessor returned %d\n stderr:\n%s" % (p.returncode, stderr))

        # Remove leading hash symbols added by CPP
        if not remove_lhash:
            return stdout
        else:
            return "\n".join(str(l) for l in stdout.splitlines() if not l.startswith("#"))


class FortranBacktrace(object):
    def __init__(self, text):
        self.text = text
        self.trace = []
        self.parse()

    def __str__(self):
        return str(self.trace)

    def parse(self):
        raise NotImplementedError("parse method must be implemented by the subclass")

    def locate_srcfile(self, base_name):
        top = find_top_build_tree(start_path=".", with_abinit=True)
        top = os.path.join(top, "src")

        for dirpath, dirnames, filenames in os.walk(top):
            if base_name in filenames:
                apath = os.path.join(dirpath, base_name)
                return apath
        else:
            cprint("Cannot find file: %s" % base_name, "red")
            return None

    def edit_source(self, editor=None):
        if not self.trace:
            return

        if editor is None:
            editor = Editor()
        src_file, lineno = self.trace[0]
        src_file = self.locate_srcfile(src_file)

        return editor.edit_file(src_file, lineno=lineno)


class NagBacktrace(FortranBacktrace):

    def parse(self):
        # Example
        #
        # Runtime Error: opernl4a_cpp.f90, line 871: INTEGER(int32) overflow for 2146435072 * 3
        # Program terminated by fatal error
        # opernl4a_cpp.f90, line 871: Error occurred in OPERNL4A
        if not self.text:
            return

        # MAGIC = "Program terminated by fatal error"
        # for i, line in enumerate(self.text):
        #    if MAGIC in line: break
        # else:
        #     return

        re_nagline = re.compile(r"(\w+\.f90), line (\d+): (.+)")

        for line in self.text:
            m = re_nagline.match(line)
            if not m:
                continue
            src_file, lineno = m.group(1), m.group(2)
            self.trace.append((src_file, int(lineno)))


class BuildEnvironment(object):
    """Store information on the build environment."""

    def __init__(self, build_dir, cygwin_instdir=None):
        """
        Args:
            build_dir: Path to the top level directory of the build.
            cygwin_instdir: Installation directory of cygwin. Defaults to '/cygwin'
        """
        # Try to figure out the top level directory of the build tree.
        try:
            build_dir = find_top_build_tree(build_dir)
        except Exception as e:
            raise e

        self.uname = platform.uname()
        self.hostname = gethostname().split(".")[0]
        self.username = my_getlogin()

        self.build_dir = os.path.abspath(build_dir)
        self.configh_path = os.path.join(self.build_dir, "config.h")
        self.binary_dir = os.path.join(self.build_dir, "src", "98_main")

        self._cygwin_instdir = ""
        if cygwin_instdir is not None:
            self._cygwin_instdir = cygwin_instdir

        # Binaries that are not located in src/98_main
        self._external_bins = {
            #"atompaw": os.path.join(self.build_dir, "fallbacks", "exports", "bin", "atompaw-abinit"),
            "atompaw": os.path.join(self.build_dir, "src", "98_main", "atompaw"),
            "timeout": os.path.join(self.build_dir, "tests", "Timeout", "timeout"),
        }

        # Check if this is a valid ABINIT build tree.
        if not (os.path.isfile(self.configh_path) and os.path.isfile(self.path_of_bin("abinit"))):
            raise ValueError("%s is not a valid ABINIT build tree." % self.build_dir)

        # Get the list of CPP variables defined in the build.
        self.defined_cppvars = parse_configh_file(self.configh_path)

        # Get info on the compilers
        self.fortran_compiler = FortranCompiler.from_defined_cpp_vars(self.defined_cppvars)
        # print(self.fortran_compiler)
        # if not self.has_bin("timeout"): print("Cannot find timeout executable!")

        self.buildbot_builder = None

    @lazy__str__
    def __str__(self): pass

    def issrctree(self):
        """True if this is a source tree."""
        configac_path = os.path.join(self.build_dir, "configure.ac")
        abinitF90_path = os.path.join(self.build_dir, "src", "98_main", "abinit.F90")

        return os.path.isfile(configac_path) and os.path.isfile(abinitF90_path)

    def iscygwin(self):
        """True if we are running under CYGWIN"""
        return "CYGWIN" in self.uname[0].upper()

    def _addext(self, string):
        """Append .exe extension, needed for cygwin"""
        if self.iscygwin(): string += ".exe"
        return string

    def path_of_bin(self, bin_name, try_syspath=True):
        """Return the absolute path of bin_name."""
        if bin_name in self._external_bins:
            bin_path = self._external_bins[bin_name]
        else:
            bin_path = os.path.join(self.binary_dir, bin_name)  # It's in src/98_main

        bin_path = self._addext(bin_path)

        # Handle external bins that are installed system wide (such as atompaw on woopy)
        if bin_name in self._external_bins and not os.path.isfile(bin_path):
            if not try_syspath:
                return ""
            # Search it in PATH.
            paths = os.getenv("PATH").split(os.pathsep)
            for p in paths:
                bin_path = os.path.join(p, bin_name)
                if os.path.isfile(bin_path):
                    break
            else:
                # err_msg = ("Cannot find path of bin_name %s, neither in the build directory nor in PATH %s" %
                #           (bin_name, paths))
                # warnings.warn(err_msg)
                bin_path = ""

        return bin_path

    def cygwin_path_of_bin(self, bin_name):
        """
        Mangle the name of the executable. Needed for Windows
        when we have to call an executable that is not located
        within the CYGWIN filesystem (aka $Win$ application).
        """
        path = self.path_of_bin(bin_name)
        if self.iscygwin():
            path = self._cygwin_instdir + path
        return path

    def has_bin(self, bin_name, try_syspath=True):
        """True if binary bin_name is present in the build."""
        return os.path.isfile(self.path_of_bin(bin_name, try_syspath=try_syspath))

    def set_buildbot_builder(self, builder):
        """
        Set the name of the buildbot builder.
        Used to skip tests defining `exclude_builders` in the TEST_INFO_SECTION
        """
        self.buildbot_builder = builder


def parse_configh_file(fname):
    """
    Parse the configuration file config.h,
    Returns a list with the CCP variables that are #defined.

    Note:
      Not very robust. It does not handle instructions such as:

      #ifdef HAVE_FOO
      #  define HAVE_BAR 1
      #endif

    Handling this case would require a real preprocessing with CPP and then the parsing.
    Not easy to implement in a portable way especially on IBM machines with XLF.
    """
    with open(fname, "rt") as fh:
        defined_cppvars = {}
        for l in fh:
            l = l.lstrip()
            if l.startswith("#define "):
                tokens = l.split()
                varname = tokens[1]
                if len(tokens) >= 3:
                    value = tokens[2]
                    defined_cppvars[varname] = value

        return defined_cppvars


def input_file_has_vars(fname, ivars, comment="#", mode="any"):
    """
    Primitive parser that searches for the occurrence of input variables in the input file fname

    Args:
        fname:
            Input file
        ivars:
            dictionary whose keys are strings with the input variables to search.
            ivar[varname] can be either None or an integer
            if ivar[varname] is None, we have a match if varname is present
            if ivar[varname] is int, we have a match if varname is present and it has value int

        mode: "all" or "any"

        return:
            (bool, d)
            bool is True is the input file contains the specified variables
            d is a dictionary with the matching lines (empty dict if no occurence).
    """
    # This algorithm is not very robust as it assumes that the variable and the line
    # are placed on the same line.
    with open(fname, "rt") as fh:
        lines = []
        for line in fh:
            line = line.lower().strip()
            idx = line.find(comment)
            if idx != -1:
                line = line[:idx]
            lines.append(line)

    matches = {}
    for k in ivars:
        matches[k] = []

    items = ivars.items()

    re_ivars = {}
    for varname in ivars:
        re_ivars[varname] = re.compile(varname + r"\d*\s*(\d+)\s*")

    nfound = 0
    for line in lines:
        for varname, varvalue in items:
            re_match = re_ivars[varname].match(line)
            # print("match")
            if varvalue is None and varname in line:
                nfound += 1
                matches[varname].append(line)
            elif re_match:
                num = int(re_match.group(1))
                if num == int(varvalue):
                    # print line
                    matches[varname].append(line)
                    nfound += 1

    if nfound == 0:
        return False, {}

    if mode == "all":
        return all(bool(v) for v in matches.values()), matches
    elif mode == "any":
        return any(bool(v) for v in matches.values()), matches
    else:
        raise ValueError("Wrong mode %s" % mode)


def make_abitest_from_input(inp_fname, abenv, keywords=None, need_cpp_vars=None, with_np=1):
    """
    Factory function to generate a Test object from the input file inp_fname
    """
    inp_fname = os.path.abspath(inp_fname)

    parser = AbinitTestInfoParser(inp_fname)

    nprocs_to_test = parser.nprocs_to_test
    ntests = len(nprocs_to_test)

    if ntests == 0:
        nprocs_to_test = [1]
        ntests = 1

    test_info = parser.generate_testinfo_nprocs(with_np)

    # Add global cpp variables.
    test_info.add_cpp_vars(need_cpp_vars)

    # Add global keywords.
    test_info.add_keywords(keywords)

    # Single test with np processors.
    # Istanciate the appropriate subclass depending on the name of the executable. Default is BaseTest.
    cls = exec2class(test_info.executable)

    return cls(test_info, abenv)


def make_abitests_from_inputs(input_fnames, abenv, keywords=None, need_cpp_vars=None):
    """
    Factory function. Return a list of tests generated from the TEST_INFO section reported
    in the input files inp_fnames.
    """
    if is_string(input_fnames):
        input_fnames = [input_fnames]

    inp_fnames = [os.path.abspath(p) for p in input_fnames]

    out_tests = []

    while inp_fnames:
        inp_fname = inp_fnames.pop(0)

        # print("inp_fname", inp_fname)
        parser = AbinitTestInfoParser(inp_fname)
        nprocs_to_test = parser.nprocs_to_test

        if len(nprocs_to_test) == 0:
            nprocs_to_test = [1]

        if not parser.is_testchain:
            # No dependency --> generate a list of test by changing the number np of MPI processors.
            for np in nprocs_to_test:
                test_info = parser.generate_testinfo_nprocs(np)

                test_info.add_cpp_vars(need_cpp_vars)  # Add global cpp variables.
                test_info.add_keywords(keywords)       # Add global keywords.

                # Istanciate the appropriate subclass depending on the name of the executable. Default is BaseTest.
                cls = exec2class(test_info.executable)
                out_tests.append(cls(test_info, abenv))

        else:
            logger.info("got chain input %s" % inp_fname)
            # print(parser.chain_inputs())

            # Build the test chain with np nprocessors.
            for np in nprocs_to_test:
                tchain_list = []
                for cht_fname in parser.chain_inputs():
                    t = make_abitest_from_input(cht_fname, abenv, keywords=keywords, need_cpp_vars=need_cpp_vars, with_np=np)
                    tchain_list.append(t)

                if not tchain_list:
                    raise RuntimeError("tchain_list is empty, inp_fname %s" % inp_fname)

                out_tests.append(ChainOfTests(tchain_list))

            # Remove the input files of the chain
            for s in parser.chain_inputs()[1:]:
                try:
                    idx = inp_fnames.index(s)
                except ValueError:
                    raise RuntimeError("%s not found in inp_fnames" % inp_fnames)

                inp_fnames.pop(idx)

    return out_tests


class NotALock:
    '''
    NOP context manager
    '''
    def __enter__(self):
        pass

    def __exit__(self, *args):
        pass


class BaseTestError(Exception):
    """Base Error class raised by Test objects"""


class BaseTest(object):
    """
    Base class describing a single test. Tests associated to other executables should
    sublcass BaseTest and redefine the method make_stdin.
    Then change exec2cls so that the appropriate instance is returned.
    """
    Error = BaseTestError

    # Possible status of the test.
    _possible_status = ["failed", "passed", "succeeded", "skipped", "disabled"]

    def __init__(self, test_info, abenv):
        logger.info("Initializing BaseTest from inp_fname: ", test_info.inp_fname)

        self._rid = genid()

        self.inp_fname = os.path.abspath(test_info.inp_fname)
        self.abenv = abenv
        self.id = test_info.make_test_id()  # The test identifier (takes into account the multi_parallel case)
        self.nprocs = 1  # Start with 1 MPI process.

        # FIXME Assumes inp_fname is in the form tests/suite_name/Input/name.in
        suite_name = os.path.dirname(self.inp_fname)
        suite_name = os.path.dirname(suite_name)

        self.suite_name = os.path.basename(suite_name)
        self.ref_dir = abenv.apath_of("tests", suite_name, "Refs")
        self.inp_dir = abenv.apath_of("tests", suite_name, "Input")

        self._executed = False
        self._status = None
        self._isok = None
        self.stdout_fname = None
        self._print_lock = NotALock()
        self.exec_error = False

        self.had_timeout = False
        self.force_skip = False

        if os.path.basename(self.inp_fname).startswith("-"):
            self._status = "disabled"

        # Initial list of local files that should not be removed.
        self._files_to_keep = []

        # Default values.
        self.make_html_diff = 0   # 0 => Do not produce diff files in HTML format
                                  # 1 => Produced HTML diff but only if test failed
                                  # 2 => Produce HTML diff independently of the final status

        self.sub_timeout = 30     # Timeout for subprocesses (in seconds)

        self.erase_files = 2      # 0 => Keep all files.
                                  # 1 => Remove files but only if the test passes or succeeds
                                  # 2 => Remove files even when the test fail.

        # Incorporate the attributes of test_info in self.
        err_msg = ""
        for k in test_info.__dict__:
            if k in self.__dict__ and test_info.__dict__[k] != self.__dict__[k]:
                err_msg += "Cannot overwrite key %s\n" % k
                # print(test_info.__dict__[k],  self.__dict__[k])

        if err_msg:
            raise self.Error(err_msg)

        self.__dict__.update(test_info.__dict__)

        if self.no_check:
            self.files_to_test = []
        elif not self.files_to_test:  # no file to test
            raise ValueError(
                self.full_id + 'This test have no files_to_test attribute.'
                ' It is forbidden unless you had "no_check = yes" to its'
                ' [setup] section in test configuration.'
            )

        # Save authors' second names to speed up the search.
        # Well, let's hope that we don't have authors with the same second name!
        second_names = []
        for string in self.authors:
            idx = string.rfind(".")
            f, s = ("", string)
            if idx != -1:
                try:
                    f, s = string[:idx + 1], string[idx + 2:]
                except IndexError:
                    raise ValueError("Wrong author(s) name")

            if not f and s and s != "Unknown":
                print("author(s) first name is missing in file %s, string = %s " % (self.full_id, string))

            second_names.append(s)

        self._authors_snames = set(second_names)

    def __repr__(self):
        return self.full_id

    def __str__(self):
        return repr(self)

    def stdin_readlines(self):
        return lazy_readlines(self.stdin_fname)

    def stdin_read(self):
        return lazy_read(self.stdin_fname)

    def stdout_readlines(self):
        return lazy_readlines(self.stdout_fname)

    def stdout_read(self):
        return lazy_read(self.stdout_fname)

    def stderr_readlines(self):
        return lazy_readlines(self.stderr_fname)

    def stderr_read(self):
        return lazy_read(self.stderr_fname)

    def cprint(self, msg='', color=None):
        with self._print_lock:
            if color is not None:
                cprint(msg, color)
            else:
                print(msg)

    @property
    def has_empty_stderr(self):
        return not bool(self.stderr_read())

    @property
    def full_id(self):
        """Full identifier of the test."""
        return "[%s][%s][np=%s]" % (self.suite_name, self.id, self.nprocs)

    @property
    def bin_path(self):
        """The absolute path of the executable needed to run the test."""
        return self.build_env.path_of_bin(self.executable)

    def cpkl_dump(self, protocol=-1):
        """Save the instance in a pickle file"""
        self.cpkl_fname = os.path.join(self.workdir, self.id + ".cpkl")

        with open(self.cpkl_fname, "wb") as fh:
            pickle.dump(self, fh, protocol=protocol)
            self.keep_files(self.cpkl_fname)

    def has_keywords(self, keywords, mode="any"):
        """
        True if test has keywords
        mode == "all" --> check if all keywords are present
        mode == "any" --> check if at least one keyword is present
        """
        if mode == "all":
            return set(keywords).issubset(self.keywords)
        elif mode == "any":
            return set(keywords).intersection(self.keywords)
        else:
            raise ValueError("wrong mode %s" % mode)

    def has_authors(self, authors, mode="any"):
        """
        True if test has authors

        mode == "all" --> check if all authors are present
        mode == "any" --> check if at least one author is present
        """
        if mode == "all":
            return set(authors).issubset(self._authors_snames)
        elif mode == "any":
            return set(authors).intersection(self._authors_snames)
        else:
            raise ValueError("wrong mode %s" % mode)

    def get_varname_set(self):
        """
        Return set of variables used by this test.
        Mainly used to check if all variables in the doc are documented/tested.

        .. note:

            Dataset index (if any) is removed.
        """
        # See abio.abivars.AbinitInputParser
        import io
        lines = []
        with io.open(self.inp_fname, "rt", encoding="utf-8") as fh:
            for line in fh:
                line.strip()
                # Remove comments from lines.
                i = line.find("#")
                if i != -1:
                    line = line[:i]
                i = line.find("!")
                if i != -1:
                    line = line[:i]
                if line:
                    lines.append(line)

        vnames = []
        # 1) Build string of the form "var1 value1 var2 value2"
        tokens = " ".join(lines).split()
        for pos, tok in enumerate(tokens):
            if tok[0].isalpha():
                # TODO
                # Either new variable, string defining the unit or operator e.g. sqrt
                # if is_abiunit(tok) or tok in ABI_OPERATORS or "?" in tok:
                #    continue

                # Have new variable
                if tok[-1].isdigit():
                    # and "?" not in tok:
                    # Handle dataset index.
                    # l = []
                    for i, c in enumerate(tok[::-1]):
                        if c.isalpha():
                            break
                        # l.append(c)
                    else:
                        raise ValueError("Cannot find dataset index in token: %s" % tok)
                    tok = tok[:len(tok) - i]
                    # l.reverse()
                    # print("tok", tok, l)
                    # tok = l
                vnames.append(tok)

        # print(vnames)
        return set(v.lower() for v in vnames)

    def has_variables(self, ivars, mode="any"):
        """True if test has the input variables ivars (dict {varname:varvalue})"""
        found, d = input_file_has_vars(self.inp_fname, ivars, mode=mode)
        return found

    def edit_input(self, editor=None):
        """
        Call editor to edit the input file of the test.
        A default editor is provided if editor is None (use $EDITOR shell variable)
        """
        if editor is None:
            editor = Editor()
        try:
            editor.edit_file(self.inp_fname)
        except Exception as e:
            raise e

    def listoftests(self, width=100, html=True, abslink=True):
        string = self.description.lstrip()
        if self.references:
            string += "References:\n" + "\n".join(self.references)
        string = textwrap.dedent(string)
        string = textwrap.fill(string, width=width)
        if not html:
            return self.full_id + ":\n" + string
        else:
            if abslink:
                link = html_link(self.full_id, self.inp_fname)
            else:
                # Use relative path so that we can upload the HTML file on
                # the buildbot master and browse the pages.
                link = html_link(self.full_id, os.path.basename(self.inp_fname))
            string = link + "<br>" + string.replace("\n", "<br>") + "\n"
        return string

    def make_stdin(self):
        """
        Generate the standard input of the test.
        The base implementation writes the content of inp_fname to stdin.
        Subclasses should redefine this method according to their needs.
        """
        t_stdin = StringIO()
        with open(self.inp_fname, "rt") as fh:
            t_stdin.writelines(fh)

        return t_stdin.getvalue()

    def get_pseudo_paths(self, dir_and_names=False):
        """
        Return list of absolut paths for pseudos.
        If `dir_and_names` is True, the function returns (dirname, basenames)
        where dirname is the common directory and basenames is a list of basenames in dirname.
        If a common directory cannot be found, dirname is set to None and basename is a list of absolute paths.
        """
        # Path to the pseudopotential files.
        # 1) pp files are searched in psps_dir first then in workdir.
        psp_paths = [os.path.join(self.abenv.psps_dir, pname) for pname in self.psp_files]

        for i, psp in enumerate(psp_paths):
            if not os.path.isfile(psp):
                pname = os.path.join(self.workdir, os.path.basename(psp))
                if os.path.isfile(pname):
                    # Use local pseudo. This is needed for atompaw tests.
                    psp_paths[i] = pname
                else:
                    err_msg = "Cannot find pp file %s, neither in Psps_for_tests nor in self.workdir" % pname
                    self.exceptions.append(self.Error(err_msg))

        if not dir_and_names:
            return psp_paths

        dirnames = [os.path.dirname(p) for p in psp_paths]
        basenames = [os.path.basename(p) for p in psp_paths]
        dirname = None
        if all(d == dirnames[0] for d in dirnames): dirname = dirnames[0]
        if dirname is not None:
            return dirname, basenames
        else:
            return None, psp_paths

    def get_extra_inputs(self):
        """Copy extra inputs from inp_dir to workdir."""
        # First copy the main input file (useful for debugging the test)
        # Avoid raising exceptions as python processes do not handle them correctly.
        try:
            src = self.inp_fname
            dest = os.path.join(self.workdir, os.path.basename(self.inp_fname))
            shutil.copy(src, dest)
            self.keep_files(dest)  # Do not remove it after the test.
        except Exception:
            self.exceptions.append(self.Error("copying %s => %s" % (src, dest)))

        for extra in self.extra_inputs:
            src = os.path.join(self.inp_dir, extra)
            dest = os.path.join(self.workdir, extra)

            if not os.path.isfile(src):
                self.exceptions.append(self.Error("%s: no such file" % src))
                continue

            shutil.copy(src, dest)
            if dest.endswith(".gz"):  # Decompress the file
                unzip(dest)
                dest = dest[:-3]
            # self.keep_files(dest)  # Do not remove dest after the test.

    @property
    def inputs_used(self):
        """List with the input files used by the test."""
        inputs = [self.inp_fname] + [os.path.join(self.inp_dir, f) for f in self.extra_inputs]

        # Add files appearing in the shell sections.
        for cmd_str in (self.pre_commands + self.post_commands):
            if cmd_str.startswith("iw_"):
                tokens = cmd_str.split()
                inp = os.path.join(self.inp_dir, tokens[1])
                inputs.append(inp)

        return inputs

    @property
    def status(self):
        """The status of the test"""
        if self._status is None:
            if self.no_check:
                self._status = "succeeded"
            else:
                all_fldstats = {f.fld_status for f in self.files_to_test}
                if "failed" in all_fldstats:
                    self._status = "failed"
                elif "passed" in all_fldstats:
                    self._status = "passed"
                else:
                    assert all_fldstats == {"succeeded"}, (
                        "Unexpected test status: {}".format(all_fldstats))
                    self._status = "succeeded"

        return self._status

    @property
    def isok(self):
        """Return true if test is OK (test passed and not python exceptions."""
        if self._isok is None:
            self._isok = self.fld_isok and not self.exceptions
        return self._isok

    @property
    def files_to_keep(self):
        """List with the files that should not be erased once the test completed"""
        return self._files_to_keep

    def keep_files(self, files):
        """Add files to the list of paths that should not be erased"""
        if is_string(files):
            self._files_to_keep.append(files)
        else:
            self._files_to_keep.extend(files)

    def compute_nprocs(self, build_env, nprocs, runmode):
        """
        Compute the number of MPI processes that can be used for the test from the initial guess nprocs

        Return: (nprocs, string)

        where nprocs = 0 if the test cannot be executed.
        string contains a human-readable message explaining the reason why the test will be skipped.

        A test cannot be executed if:

          1) It requires CPP variables that are not defined in the build.
          2) The user asks for more MPI nodes than max_nprocs (this value is reported in the TEST_INFO section).
          3) We have a multiparallel test (e.g. paral/tA.in) and nprocs is not in in nprocs_to_test
          4) nprocs is in exclude_nprocs
        """
        # !HAVE_FOO --> HAVE_FOO should not be present.
        errors = []
        eapp = errors.append
        for var in self.need_cpp_vars:
            if not var.startswith("!") and var not in build_env.defined_cppvars:
                eapp("Build environment does not define the CPP variable %s" % var)
            elif var[1:] in build_env.defined_cppvars:
                eapp("Build environment defines the CPP variable %s" % var[1:])

        # Remove this check to run the entire test suite in parallel
        # runmode ="dynamic"

        if runmode == "static":
            if nprocs > self.max_nprocs:
                eapp("nprocs: %s > max_nprocs: %s" % (nprocs, self.max_nprocs))

        elif runmode == "dynamic":
            # Will select the minimum between max_nprocs and nprocs
            pass

        else:
            raise ValueError("Wrong runmode %s" % runmode)

        if self.nprocs_to_test and nprocs != self.nprocs_to_test[0]:
            eapp("nprocs: %s != nprocs_to_test: %s" % (nprocs, self.nprocs_to_test[0]))

        if nprocs in self.exclude_nprocs:
            eapp("nprocs: %s in exclude_nprocs: %s" % (nprocs, self.exclude_nprocs))

        if self.force_skip:
            eapp("forced to be skipped by the chain of test.")

        err_msg = "\n".join(errors)
        if err_msg:
            real_nprocs = 0
        else:
            real_nprocs = min(self.max_nprocs, nprocs)

        # if err_msg: print(err_msg)
        return real_nprocs, err_msg

    def skip_host(self):
        """
        Return True if the test should be skipped since we are running on a banned host.
        """
        compilers, slaves = [], []

        for s in self.exclude_hosts:
            compiler, host = None, s
            if "@" in s:
                compiler, host = s.split("@")
            else:
                # TODO: validate TEST_INFO at the level of the parser.
                warnings.warn("Wrong string %s in exclude_hosts" % s)

            compilers.append(compiler)
            slaves.append(host)

        # Find the slave and compare the name of the compiler.
        try:
            idx = slaves.index(self.build_env.hostname)
        except ValueError:
            return False

        return compilers[idx] == self.build_env.fortran_compiler.name

    def skip_buildbot_builder(self):
        """
        Return True if the test should be skipped since we are running on a banned builder.
        """
        if getattr(self.build_env, "buildbot_builder", None) is None:
            return False
        for builder in self.exclude_builders:
            if any(c in builder for c in "*?![]{}"):
                # Interpret builder as regex.
                m = re.compile(builder)
                if m.match(self.build_env.buildbot_builder):
                    return True
            else:
                if builder == self.build_env.buildbot_builder:
                    return True
        return False

    def run(self, build_env, runner, workdir, print_lock=None, nprocs=1,
            runmode="static", **kwargs):
        """
        Run the test with nprocs MPI nodes in the build environment build_env using the `JobRunner` runner.
        Results are produced in directory workdir. kwargs is used to pass additional options

        ================  ====================================================================
        kwargs            Meaning
        ================  ====================================================================
        pedantic           Mark tests as failed if stderr is not empty.
        erase_files        0 => Keep all files produced by the test
                           1 => Remove files but only if the test passed or succeeded.
                           2 => Remove files even if the test failed.
                           default=2
        make_html_diff     True to produce diff in HTML format. Default: False.
        sub_timeout        Timeout for subprocesses.
        abimem_check       True if abimem.mocc files should be analyzes for possible errors.
                           Requires HAVE_MEM_PROFILE and `call abimem_init(2)` in main.
                           Default: False
        etsf_check         True if netcdf files should be validated. Requires netcdf4.
                           Default: False
        ================  ====================================================================

        .. warning:
            This method must be thread-safe, DO NOT change build_env or runner.
        """
        import copy
        runner = copy.deepcopy(runner)
        start_time = time.time()

        if print_lock is not None:
            self._print_lock = print_lock

        workdir = os.path.abspath(workdir)
        if not os.path.exists(workdir): os.mkdir(workdir)
        self.workdir = workdir

        self.build_env = build_env

        self.exceptions = []
        self.fld_isok = True  # False if at least one file comparison fails.

        # Extract options from kwargs
        self.pedantic = kwargs.get("pedantic", False)
        self.erase_files = kwargs.get("erase_files", self.erase_files)
        self.make_html_diff = kwargs.get("make_html_diff", self.make_html_diff)
        self.sub_timeout = kwargs.get("sub_timeout", self.sub_timeout)

        timeout = self.sub_timeout
        if self.build_env.has_bin("timeout") and timeout > 0.0:
            exec_path = self.build_env.path_of_bin("timeout")
            self.timebomb = TimeBomb(timeout, delay=0.05, exec_path=exec_path)
        else:
            self.timebomb = TimeBomb(timeout, delay=0.05)

        status2txtcolor = {
            "succeeded": "green",
            "passed": "blue",
            "failed": "red",
            "disabled": "cyan",
            "skipped": "cyan",
        }

        # Check whether the test can be executed.
        can_run = True
        if self._status == "disabled":
            msg = self.full_id + ": Disabled"
            can_run = False
            self.cprint(msg, status2txtcolor[self._status])

        # Here we get the number of MPI nodes for test.
        self.nprocs, self.skip_msg = self.compute_nprocs(self.build_env, nprocs, runmode=runmode)

        if self.skip_msg:
            self._status = "skipped"
            msg = self.full_id + ": Skipped."
            self.cprint(msg, status2txtcolor[self._status])
            for l in self.skip_msg.splitlines():
                self.cprint("\t" + l, status2txtcolor[self._status])
            self.cprint()
            can_run = False

        if self.skip_host():
            self._status = "skipped"
            msg = self.full_id + ": Skipped: this hostname has been excluded."
            self.cprint(msg, status2txtcolor[self._status])
            can_run = False

        if self.skip_buildbot_builder():
            self._status = "skipped"
            msg = self.full_id + ": Skipped: this buildbot builder has been excluded."
            self.cprint(msg, status2txtcolor[self._status])
            can_run = False

        self.run_etime = 0.0

        if can_run:
            # Execute pre_commands in workdir.
            rshell = RestrictedShell(self.inp_dir, self.workdir, self.abenv.psps_dir)

            for cmd_str in self.pre_commands:
                rshell.execute(cmd_str)

            if rshell.exceptions:
                self.exceptions.extend(rshell.exceptions)
                rshell.empty_exceptions()

            # Copy extra inputs in workdir (if any).
            self.get_extra_inputs()

            # Create stdin file in the workdir.
            self.stdin_fname = os.path.join(self.workdir, self.id + ".stdin")
            self.stdout_fname = os.path.join(self.workdir, self.id + ".stdout")
            self.stderr_fname = os.path.join(self.workdir, self.id + ".stderr")

            # Run the code (run_etime is the wall time spent to execute the test)
            # FIXME: Add support for more executables
            use_files_file = self.use_files_file
            if self.executable not in ("abinit", "anaddb"):
                use_files_file = True

            if use_files_file:
                self.keep_files([self.stdin_fname, self.stdout_fname, self.stderr_fname])
                # Legacy mode: create files file and invoke exec with syntax: `abinit < run.files`
                with open(self.stdin_fname, "wt") as fh:
                    fh.writelines(self.make_stdin())
                stdin_fname = self.stdin_fname
                bin_argstr = self.exec_args

            else:
                # New CLI mode: invoke exec with syntax `abinit run.abi`
                # stdin_fname won't be created
                self.keep_files([self.stdout_fname, self.stderr_fname])
                stdin_fname = ""
                self.prepare_new_cli_invokation()

                path = os.path.join(self.workdir, os.path.basename(self.inp_fname))
                bin_argstr = path + self.exec_args
                print("Using .abi mode with bin_argstr", bin_argstr)

            self.run_etime = runner.run(self.nprocs, self.bin_path,
                                        stdin_fname, self.stdout_fname, self.stderr_fname,
                                        bin_argstr=bin_argstr, cwd=self.workdir)

            # Save exceptions (if any).
            if runner.exceptions:
                self.exec_error = True
                self.exceptions.extend(runner.exceptions)
                if not self.expected_failure:
                    for exc in runner.exceptions:
                        self.cprint(exc)

            # Execute post_commands in workdir.
            for cmd_str in self.post_commands:
                rshell.execute(cmd_str)

            # Save exceptions (if any).
            if rshell.exceptions:
                self.exceptions.extend(rshell.exceptions)
                rshell.empty_exceptions()

            # Check final results:
            # 1) use fldiff to compare ref and output files.
            # 2) fldiff stdout is redirected to fldiff_fname.
            for f in self.files_to_test:
                fldiff_fname = os.path.join(self.workdir, f.name + ".fldiff")
                self.keep_files(fldiff_fname)

                with open(fldiff_fname, "w") as fh:
                    f.fldiff_fname = fldiff_fname

                    isok, status, msg = f.compare(self.abenv.fldiff_path, self.ref_dir, self.workdir,
                                                  yaml_test=self.yaml_test, timebomb=self.timebomb, outf=fh)
                self.keep_files(os.path.join(self.workdir, f.name))
                self.fld_isok = self.fld_isok and isok

                if not self.exec_error and f.has_line_count_error:
                    f.do_html_diff = True

                msg = ": ".join([self.full_id, msg])
                self.cprint(msg, status2txtcolor[status])

            # Check if the test is expected to fail.
            if runner.retcode == 124:
                self._status = "failed"
                self.had_timeout = True
                msg = self.full_id + "test has reached timeout and has been killed (SIGTERM)."
                self.cprint(msg, status2txtcolor["failed"])

            elif runner.retcode == 137:
                self._status = "failed"
                self.had_timeout = True
                msg = self.full_id + "test has reached timeout and has been killed (SIGKILL)."
                self.cprint(msg, status2txtcolor["failed"])

            elif runner.retcode != 0 and not self.expected_failure:
                self._status = "failed"
                msg = (self.full_id + "Test was not expected to fail but subprocesses returned %s" % runner.retcode)
                self.cprint(msg, status2txtcolor["failed"])

            # If pedantic, stderr must be empty unless the test is expected to fail!
            if self.pedantic and not self.expected_failure:
                try:
                    errout = self.stderr_read()
                    if errout:
                        # TODO: Not very clean, I should introduce a new status and a setter method.
                        self._status = "failed"
                except Exception as exc:
                    self.exceptions.append(exc)

            # Check stderr for presence of valgrind errors.
            if runner.has_valgrind:
                try:
                    # Build a parser from the command line options and parse the stderr.
                    parser = runner.build_valgrind_parser()
                    parser.parse(self.stderr_fname)

                    if parser.error_report:
                        # TODO: Not very clean, I should introduce a new status and a setter method.
                        self._status = "failed"
                        msg = " ".join([self.full_id, "VALGRIND ERROR:", parser.error_report])
                        self.cprint(msg, status2txtcolor["failed"])

                except Exception as exc:
                    self.exceptions.append(exc)

            if self.status == "failed":
                # Print the first line of the stderr if it's not empty.
                # Look also for MPIABORTFILE
                try:
                    errout = self.stderr_read()
                    if errout:
                        self.cprint(errout, status2txtcolor["failed"])

                    # Extract YAML error message from ABORTFILE or stdout.
                    abort_file = os.path.join(self.workdir, "__ABI_MPIABORTFILE__")
                    if os.path.exists(abort_file):
                        with open(abort_file, "rt") as f:
                            self.cprint(12 * "=" + " ABI_MPIABORTFILE " + 12 * "=")
                            self.cprint(f.read(), status2txtcolor["failed"])
                            f.close()
                    else:
                        yamlerr = read_yaml_errmsg(self.stdout_fname)
                        if yamlerr:
                            self.cprint("YAML Error found in the stdout of: " + repr(self))
                            self.cprint(yamlerr, status2txtcolor["failed"])
                        else:
                            self.cprint("No YAML Error found in: " + repr(self))

                except Exception as exc:
                    self.exceptions.append(exc)

            if kwargs.get("abimem_check", False):
                paths = [os.path.join(self.workdir, f) for f in os.listdir(self.workdir)
                         if f.startswith("abimem") and f.endswith(".mocc")]
                self.cprint("Found %s abimem files" % len(paths))
                # abimem_retcode = 0
                for path in paths:
                    memfile = AbimemFile(path)
                    memfile.find_memleaks()
                    #if rc: parser.show_errors()
                    #abimem_retcode += rc

            # if False and kwargs.get("etsf_check", False):
            if kwargs.get("etsf_check", False):
                # Mark the test as failed and create a custom Exception
                # developers will have to inspect the xreport file for the full list of errors.
                try:
                    from . import etsf_specs as etsf
                except ImportError:
                    etsf = None

                errmsg = ""
                if etsf is None:
                    errmsg = "etsf_check is activated but netcdf4 module is not available"
                    nc_retcode = 1

                else:
                    nc_retcode = 0
                    all_errors = []
                    for p in os.listdir(self.workdir):
                        if p.endswith(".nc"):
                            path = os.path.join(self.workdir, p)
                            elist = []
                            # elist += etsf.validate_vars(path))
                            elist += etsf.validate_ncfile(path)

                            if elist:
                                all_errors.append(elist)
                                self.cprint("%s [FAILED]" % p, "red")
                            else:
                                self.cprint("%s [OK]" % p, "green")

                    nc_retcode = len(all_errors)

                    if nc_retcode != 0:
                        errmsg = ("Setting status to failed because nc_retcode=%s\n"
                                  "The netcdf files produced by this tests either is not consistent with the etsf specs.\n"
                                  "or it has not been registered in ~abinit/tests/pymods/etsf_specs.py\n"
                                  "Please, control the errors messages in the xreport file produced by buildbot."
                                  ) % nc_retcode

                if nc_retcode != 0:
                    # TODO: Not very clean, I should introduce a new status and a setter method.
                    self._status = "failed"
                    # Store the exception and continue.
                    self.exceptions.append(Exception(errmsg))
                    self.cprint(errmsg)
                else:
                    self.cprint("netcdf validation [OK]", "green")

        self._executed = True
        self.tot_etime = time.time() - start_time

    def results_load(self, d):
        """
        Load the run results from a run in a different process.
        """
        self._status = d['status']
        self.stdout_fname = d['stdout']
        self._files_to_keep = d['files_to_keep']
        self.tot_etime = d['tot_etime']
        self.run_etime = d['run_etime']
        self._executed = d['executed']
        self._isok = d['isok']
        self.exec_error = d['exec_error']
        self.workdir = d['workdir']

    def results_dump(self, skipped_info=False):
        """
        Dump the run results to pass it to a different process
        """
        return {
            'id': self._rid,
            'status': self.status,
            'stdout': self.stdout_fname,
            'files_to_keep': self.files_to_keep,
            'tot_etime': self.tot_etime,
            'run_etime': self.run_etime,
            'executed': self._executed,
            'exec_error': self.exec_error,
            'isok': self.isok,
            'workdir': self.workdir,
        }

    @property
    def executed(self):
        return self._executed

    def clean_workdir(self, other_test_files=None):
        """Remove the files produced in self.workdir."""
        assert self._executed
        if not os.path.exists(self.workdir) or self.erase_files == 0:
            return

        save_files = self._files_to_keep[:]
        if other_test_files is not None:
            save_files += other_test_files

        # Add harcoded list of files
        hard_files = ["perf.data", "__ABI_MPIABORTFILE__"]
        save_files += [os.path.join(self.workdir, f) for f in hard_files]

        # List of file extensions to be preserved.
        keep_exts = [".flun", ".mocc"]

        if (self.erase_files == 1 and self.isok) or self.erase_files == 2:
            entries = [os.path.join(self.workdir, e) for e in os.listdir(self.workdir)]
            for entry in entries:
                if entry in save_files: continue
                _, ext = os.path.splitext(entry)
                if ext in keep_exts: continue
                if os.path.isfile(entry):
                    try:
                        os.remove(entry)
                    except OSError:
                        pass
                else:
                    raise NotImplementedError("Found directory: %s in workdir!!" % entry)

    def patch(self, patcher=None):
        """
        Patch the output files of the test with the specified patcher.
        A default patcher is provided if patcher is None (use $PATCHER shell variable)
        """
        assert self._executed
        from tests.pymods import Patcher
        for f in self.files_to_test:
            ref_fname = os.path.abspath(os.path.join(self.ref_dir, f.name))
            out_fname = os.path.abspath(os.path.join(self.workdir, f.name))
            raise NotImplementedError("patcher should be tested")
            Patcher(patcher).patch(out_fname, ref_fname)

    def make_html_diff_files(self):
        """Generate and write diff files in HTML format."""
        assert self._executed
        if self.make_html_diff == 0 or self._status in {"disabled", "skipped"}:
            return

        diffpy = self.abenv.apath_of("tests", "pymods", "diff.py")

        for f in self.files_to_test:
            if not f.do_html_diff and self.make_html_diff == 1:
                continue

            ref_fname = os.path.abspath(os.path.join(self.ref_dir, f.name))

            if not os.path.isfile(ref_fname) and ref_fname.endswith(".stdout"):
                ref_fname = ref_fname[:-7] + ".out"  # FIXME Hack due to the stdout-out ambiguity

            out_fname = os.path.abspath(os.path.join(self.workdir, f.name))

            # Check whether output and ref file exist.
            out_exists = os.path.isfile(out_fname)
            ref_exists = os.path.isfile(ref_fname)

            hdiff_fname = os.path.abspath(os.path.join(self.workdir, f.name + ".diff.html"))

            f.hdiff_fname = hdiff_fname

            x, ext = os.path.splitext(f.name)
            safe_hdiff = ext in {".out", ".stdout"}  # Create HTML diff file only for these files

            if ref_exists and out_exists and safe_hdiff:
                out_opt = "-m"
                # out_opt = "-t"   # For simple HTML table. (can get stuck)
                # args = ["python", diffpy, out_opt, "-f " + hdiff_fname, out_fname, ref_fname ]
                args = [diffpy, out_opt, "-j",  "-f " + hdiff_fname, out_fname, ref_fname]
                cmd = " ".join(args)
                # print("Diff", cmd)

                p, ret_code = self.timebomb.run(cmd, shell=True, cwd=self.workdir)

                if ret_code != 0:
                    err_msg = "Timeout error (%s s) while executing %s, retcode = %s" % (
                        self.timebomb.timeout, str(args), ret_code
                    )
                    self.exceptions.append(self.Error(err_msg))
                else:
                    self.keep_files(hdiff_fname)

    def make_txt_diff_files(self):
        """Generate and write diff files in txt format."""
        assert self._executed
        if self._status in {"disabled", "skipped"}:
            return

        # print(self._status)
        # if self._status not in {"failed", "passed"}: return
        diffpy = self.abenv.apath_of("tests", "pymods", "diff.py")

        for f in self.files_to_test:
            # print(f, f.fld_isok)
            if f.fld_isok: continue

            ref_fname = os.path.abspath(os.path.join(self.ref_dir, f.name))

            if not os.path.isfile(ref_fname) and ref_fname.endswith(".stdout"):
                ref_fname = ref_fname[:-7] + ".out"  # FIXME Hack due to the stdout-out ambiguity

            out_fname = os.path.abspath(os.path.join(self.workdir, f.name))

            # Check whether output and ref file exist.
            out_exists = os.path.isfile(out_fname)
            ref_exists = os.path.isfile(ref_fname)

            diff_fname = os.path.abspath(os.path.join(self.workdir, f.name + ".diff"))

            f.diff_fname = diff_fname

            x, ext = os.path.splitext(f.name)

            if ref_exists and out_exists:
                # n is for ndiff format, c for context, u for unified
                # for out_opt in ["-n", "-c"]:
                # out_opt = "-n"
                # out_opt = "-c"
                out_opt = "-u"
                args = [diffpy, out_opt, "-j", "-f " + diff_fname, out_fname,
                        ref_fname]
                cmd = " ".join(args)

                (p, ret_code) = self.timebomb.run(cmd, shell=True, cwd=self.workdir)

                if ret_code != 0:
                    err_msg = "Timeout error (%s s) while executing %s, retcode = %s" % (
                        self.timebomb.timeout, str(args), ret_code)
                    self.exceptions.append(self.Error(err_msg))
                else:
                    self.keep_files(diff_fname)

    def write_html_report(self, fh=None, oc="oc"):
        """Write the HTML file summarizing the results of the test."""
        assert self._executed

        close_fh = False
        if fh is None:
            close_fh = True
            html_report = os.path.join(self.workdir, "test_report.html")
            fh = open(html_report, "wt")

        self.keep_files(fh.name)

        self.make_html_diff_files()
        self.make_txt_diff_files()

        # Try to read stdout, stderr and the abort_file produced by Abinit in parallel
        # Ignore errors (fock takes years to flush the stdout)
        # stdout_text, stderr_text = 2*("",)
        nlast = 120
        stderr_text, stdout_text, abiabort_text = 3 * (" ",)
        abort_file = find_abortfile(self.workdir)
        # self.fld_isok = False
        errinfo_text = " "
        # print("fld_isok:", self.fld_isok)
        if not self.fld_isok or self.status == "failed":
            try:
                stderr_text = str2html(self.stderr_read())
                stdout_text = str2html(tail_file(self.stdout_fname, nlast))
                abiabort_text = "No __ABI_MPIABORTFILE__ found"

                if abort_file:
                    with open(abort_file, "rt") as f:
                        abiabort_text = (
                            12 * "=" + os.path.basename(abort_file)
                            + 12 * "=" + 2 * "\n" + str(f.read())
                        )

            except Exception as exc:
                s = "Exception while trying to get info from stderr, stdout and __ABI_MPIABORTFILE\n" + str(exc)
                stderr_text, stdout_text, abiabort_text = 3 * (s,)

            # Look for extra info on the error in selected files produced by the code.
            try:
                errinfo_text = str2html(extract_errinfo_from_files(self.workdir))
            except Exception as exc:
                errinfo_text = "Exception while trying to get error info from extra files\n" + str(exc)

        ##################################################
        # Document Name Space that serves as the substitution
        # namespace for instantiating a doc template.
        username = my_getlogin()

        DNS = {
            "self": self,
            "page_title": "page_title",
            "user_name": username,
            "hostname": gethostname(),
            "Headings": ['File_to_test', 'Status', 'fld_output', 'fld_options', 'txt_diff', 'html_diff'],
            "nlast": nlast,
            "stderr_text": stderr_text,
            "stdout_text": stdout_text,
            "abiabort_text": abiabort_text,
            "errinfo_text": errinfo_text,
            # Functions and modules available in the template.
            "time": time,
            "pj": os.path.join,
            "basename": os.path.basename,
            "str2html": str2html,
            "sec2str": sec2str,
            "args2htmltr": args2htmltr,
            "html_link": html_link,
            "status2html": status2html
        }

        header = """
        <html>
         <head><title>$page_title</title></head>
         <body bgcolor="#FFFFFF" text="#000000">
        """

        if self.status in {"skipped", "disabled"}:
            if self.status == "skipped":
                template = str2html(self.skip_msg)
            else:
                template = "This test has been disabled!"
        else:
            template = """
              <hr>
              <h1>Results of test ${self.full_id}</h1>
                 MPI nprocs =  ${self.nprocs},
                 run_etime = ${sec2str(self.run_etime)} s,
                 tot_etime = ${sec2str(self.tot_etime)} s
               <br>
               ${html_link("stdin",  basename(self.stdin_fname))},
               ${html_link("stdout", basename(self.stdout_fname))},
               ${html_link("stderr", basename(self.stderr_fname))}
              <p>
              <table width="100%" border="0" cellspacing="0" cellpadding="2">
                <tr valign="top" align="left">
                <py-open code = "for h in Headings:"> </py-open>
                  <th>${h}</th>
                <py-close/>
                </tr>
                <py-open>for idx, f in enumerate(self.files_to_test):</py-open>
                 <tr valign="top" align="left">
                  <py-line code = "fld_link = html_link(basename(f.fldiff_fname))"/>
                  <py-line code = "txt_diff_link = html_link(basename(f.diff_fname))"/>
                  <py-line code = "html_diff_link = html_link(basename(f.hdiff_fname))"/>
                  <py-line code = "tab_row = args2htmltr(f.name, status2html(f.fld_status), fld_link, f.fld_options, txt_diff_link, html_diff_link)"/>
                  ${tab_row}
                 </tr>
                <py-close/>
              </table>

              <py-open>for idx, f in enumerate(self.files_to_test):</py-open>
                <py-open code="if f.fld_status != 'succeeded':"/>
                <p> ${f.name} ${f.fld_msg} </p>
              <py-close/>

              <py-open code="if self.status == "failed":"/>
                <py-open code="if self.exceptions:"/>
                  <hr><p>
                  <h1>Exceptions raised at run-time:</h1>
                  <py-open code="for idx, e in enumerate(self.exceptions):"/>
                    <p> $idx) ${str2html(str(e))}</p>
                  <py-close/>
                  <br>
                <py-close/>
                <hr><p>
                <h1>Standard Error of test ${self.id}:</h1>
                  ${stderr_text}
                <hr><p>
                <h1>__MPIABORTFILE__ of test ${self.id}:</h1>
                  ${abiabort_text}
                <hr><p>
                <h1>Info extracted from debug files produced by ${self.id}:</h1>
                  ${errinfo_text}
                <hr><p>
                <h1>Standard output of test ${self.id} (last ${nlast} lines):</h1>
                  ${stdout_text}
                <br>
              <py-close/>
              <p>
              <h3>Extra Information</h3>
              <py-line code = "authors = ', '.join(a for a in self.authors)" />
              <p>Authors = ${authors}</p>
              <py-line code = "keys = ', '.join(k for k in self.keywords)" />
              <p>Keywords = ${keys}</p>
              <p>${self.listoftests(abslink=False)}</p>
            """

        footer = """
          <hr>
          Automatically generated by %s on %s. Logged on as %s@%s
          Python version: %s
          <hr>
          </body>
          </html> """ % (_MY_NAME, time.asctime(), username, gethostname(), platform.python_version())

        if "o" in oc:
            template = header + template
        if "c" in oc:
            template += footer

        # Set a file-like object to template
        template_stream = StringIO(template)

        # Initialise an xyaptu xcopier, and call xcopy
        xcp = xcopier(DNS, ouf=fh)
        xcp.xcopy(template_stream)

        if close_fh:
            fh.close()

    def _get_one_backtrace(self):
        return NagBacktrace(self.stderr_readlines())

    def get_backtraces(self):
        return [self._get_one_backtrace()]

#############################################################################################################
# Subclasses needed to handle the different executables
#############################################################################################################


class AbinitTest(BaseTest):
    """
    Class for Abinit tests. Redefine the make_stdin method of BaseTest,
    provides `prepare_new_cli_invokation`
    """
    def make_stdin(self):
        t_stdin = StringIO()

        # Use the basename instead of the absolute path because the input has been already copied
        # and we might want to change it especially if we are debugging the code
        inp_fname = self.inp_fname
        t_stdin.write(os.path.basename(inp_fname) + "\n")
        t_stdin.write(self.id + ".out" + "\n")

        # Prefix for input/output/temporary files
        i_prefix = self.input_prefix if self.input_prefix else self.id + "i"
        o_prefix = self.output_prefix if self.output_prefix else self.id + "o"
        # FIXME: Use t prefix and change iofn
        t_prefix = self.id  # + "t"

        t_stdin.writelines(l + "\n" for l in [i_prefix, o_prefix, t_prefix])

        # Path to the pseudopotential files.
        # 1) pp files are searched in psps_dir first then in workdir.
        psp_paths = self.get_pseudo_paths()

        t_stdin.writelines(p + "\n" for p in psp_paths)

        return t_stdin.getvalue()

    def prepare_new_cli_invokation(self):
        """Perform operations required to execute test with new CLI."""
        # Need to add pseudopotential info to input.
        with open(self.inp_fname, "rt") as fh:
            line = fh.read()

        extra = ["# Added by runtests.py"]
        app = extra.append
        app('output_file = "%s"' % (self.id + ".out"))

        # This is needed for ATOMPAW as the pseudo will be generated at runtime.
        #dirname, pp_names = self.get_pseudo_paths(dir_and_names=True)
        #if dirname is not None:
        #    app('pp_dirpath = "%s"' % (dirname))
        #    app('pseudos = "%s"' % (",".join(pp_names)))
        #else:
        #    app('pseudos = "%s"' % (",\n".join(pp_names)))

        # This is to check whether the parser supports "long strings"
        #app('pseudos = "%s"' % (", ".join(self.get_pseudo_paths())))

        # This is to check pseudos with an env variable defining directory
        app('pp_dirpath = "$ABI_PSPDIR"')
        app('pseudos = "%s"' % (", ".join(os.path.basename(p) for p in self.get_pseudo_paths())))

        # Prefix for input/output/temporary files
        i_prefix = self.input_prefix if self.input_prefix else self.id + "i"
        o_prefix = self.output_prefix if self.output_prefix else self.id + "o"
        # FIXME: Use temp prefix and change iofn
        t_prefix = self.id  # + "t"

        app('indata_prefix "%s"' % i_prefix)
        app('outdata_prefix "%s"' % o_prefix)
        app('tmpdata_prefix "%s"' % t_prefix)
        app("# end runtests.py section\n\n")

        path = os.path.join(self.workdir, os.path.basename(self.inp_fname))
        with open(path, "wt") as fh:
            fh.write("\n".join(extra) + line)


class AnaddbTest(BaseTest):
    """
    Class for Anaddb tests. Redefine the make_stdin method of BaseTest
    provides `prepare_new_cli_invokation`
    """

    def get_ddb_path(self):
        """Return the path to the input DDB file."""
        iddb_fname = self.id + ".ddb.in"
        if self.input_ddb:
            # Use output DDB of a previous run.
            iddb_fname = os.path.join(self.workdir, self.input_ddb)
            if not os.path.isfile(iddb_fname):
                self.exceptions.append(self.Error("%s no such DDB file: " % iddb_fname))
        return iddb_fname

    def get_gkk_path(self):
        """Return the path to the input GKK file for EPH calculations."""
        input_gkk = self.id + ".gkk"
        if self.input_gkk:
            input_gkk = os.path.join(self.workdir, self.input_gkk)  # Use output GKK of a previous run.
            if not os.path.isfile(input_gkk):
                self.exceptions.append(self.Error("%s no such GKK file: " % input_gkk))

        if not os.path.isfile(input_gkk): input_gkk = ""
        return input_gkk

    def get_ddk_path(self):
        """Return the path to the input DKK file for EPH calculations."""
        input_ddk = self.id + ".ddk"
        if not os.path.isfile(input_ddk):
            # Try in input directory:
            # FIXME: Someone has to rewrite the treatment of the anaddb files file
            input_ddk = os.path.join(self.inp_dir, input_ddk)
        if not os.path.isfile(input_ddk): input_dkk = ""
        return input_ddk

    def make_stdin(self):
        t_stdin = StringIO()

        t_stdin.write(self.inp_fname + "\n")         # 1) formatted input file
        t_stdin.write(self.id + ".out" + "\n")       # 2) formatted output file e.g. t13.out
        t_stdin.write(self.get_ddb_path() + "\n")    # 3) input derivative database e.g. t13.ddb.in
        t_stdin.write(self.id + ".md" + "\n")        # 4) output molecular dynamics e.g. t13.md
        t_stdin.write(self.get_gkk_path() + "\n")    # 5) input elphon matrix elements  (GKK file) :
        t_stdin.write(self.id + "\n")                # 6) base name for elphon output files e.g. t13
        t_stdin.write(self.get_ddk_path() + "\n")    # 7) file containing ddk filenames for elphon/transport:

        return t_stdin.getvalue()

    def prepare_new_cli_invokation(self):
        """Perform operations required to execute test with new CLI."""
        # Need to add extra variables depending on calculation type.
        with open(self.inp_fname, "rt") as fh:
            line = fh.read()

        extra = ["# Added by runtests.py"]
        app = extra.append
        app('ddb_path = "%s"' % (self.get_ddb_path()))
        app('output_file = "%s"' % (self.id + ".out"))

        # EPH stuff
        gkk_path = self.get_gkk_path()
        if gkk_path: app('gkk_path = "%s"' % gkk_path)
        ddk_path = self.get_ddk_path()
        if ddk_path: app('ddk_path = "%s"' % ddk_path)

        if gkk_path or ddk_path:
            # EPH calculation
            app('eph_prefix = "%s"' % self.id)

        app("# end runtests.py section\n\n")

        path = os.path.join(self.workdir, os.path.basename(self.inp_fname))
        with open(path, "wt") as fh:
            fh.write("\n".join(extra) + line)


class MultibinitTest(BaseTest):
    """
    Class for Multibinit tests. Redefine the make_stdin method of BaseTest
    """
    def make_stdin(self):
        t_stdin = StringIO()

        t_stdin.write(self.inp_fname + "\n")         # 1) formatted input file
        t_stdin.write(self.id + ".out" + "\n")       # 2) formatted output file e.g. t13.out

        if self.input_ddb:
            iddb_fname = os.path.join(self.inp_dir, self.input_ddb)
            if not os.path.isfile(iddb_fname):
                self.exceptions.append(self.Error("%s no such DDB file: " % iddb_fname))
            t_stdin.write(iddb_fname + "\n")         # 3) input derivative database e.g. ddb.in
        else:
            if self.system_xml:
                sys_xml_fname = os.path.join(self.inp_dir, self.system_xml)
                if not os.path.isfile(sys_xml_fname):
                    self.exceptions.append(self.Error("%s no such XML file: " % sys_xml_fname))
                t_stdin.write(sys_xml_fname + "\n")  # 3) input for system.xml XML
            else:
                self.exceptions.append(self.Error("%s no file available for the system"))

        if self.coeff_xml:
            coeffxml_fname = os.path.join(self.inp_dir, self.coeff_xml)
            if not os.path.isfile(coeffxml_fname):
                self.exceptions.append(self.Error("%s no such XML file for coeffs: " % coeffxml_fname))

            t_stdin.write(coeffxml_fname + "\n")  # 4) input for coefficients
        else:
            coeffxml_fname = "no"
            t_stdin.write(coeffxml_fname + "\n")

        if self.md_hist:
            md_hist_fname = os.path.join(self.inp_dir, self.md_hist)
            if not os.path.isfile(md_hist_fname):
                self.exceptions.append(self.Error("%s no such XML file for coeffs: " % md_hist_fname))

            t_stdin.write(md_hist_fname + "\n")  # 5) input for coefficients
        else:
            md_hist_fname = "no"
            t_stdin.write(md_hist_fname + "\n")

        return t_stdin.getvalue()


class TdepTest(BaseTest):
    """
    Class for TDEP tests. Redefine the make_stdin method of BaseTest
    """
    def make_stdin(self):
        t_stdin = StringIO()

        inp_fname = os.path.basename(self.inp_fname)
        t_stdin.write(inp_fname + "\n")              # 1) formatted input file

        md_hist_fname = os.path.join(self.inp_dir, self.md_hist)
        if not os.path.isfile(md_hist_fname):
            self.exceptions.append(self.Error("%s no such hist file: " % md_hist_fname))

        t_stdin.write(md_hist_fname + "\n")
        t_stdin.write(self.id + "\n")       # 2) formatted output file e.g. t13.out

        return t_stdin.getvalue()


class AimTest(BaseTest):
    """
    Class for Aim tests. Redefine the make_stdin method of BaseTest
    """
    def make_stdin(self):
        t_stdin = StringIO()

        t_stdin.write(self.inp_fname + "\n")   # formatted input file e.g. .../Input/t57.in

        iden_fname = self.id + "i_DEN"
        t_stdin.write(iden_fname + "\n")  # input density  e.g. t57i_DEN
        t_stdin.write(self.id + "\n")     # t57

        # Path to the pseudopotential files.
        psp_paths = [os.path.join(self.abenv.psps_dir, pname) for pname in self.psp_files]

        t_stdin.writelines(p + "\n" for p in psp_paths)

        return t_stdin.getvalue()


class ConductiTest(BaseTest):
    """
    Class for Conducti tests. Redefine the make_stdin method of BaseTest
    """
    def make_stdin(self):
        t_stdin = StringIO()
        t_stdin.write(self.inp_fname + "\n")  # formatted input file e.g. .../Input/t57.in
        t_stdin.write(self.id + "\n")         # will be used as the prefix of the log file names e.g. t57

        return t_stdin.getvalue()


class OpticTest(BaseTest):
    """
    Class for Optic tests. Redefine the make_stdin method of BaseTest
    """
    def make_stdin(self):
        t_stdin = StringIO()

        t_stdin.write(self.inp_fname + "\n")  # optic input file e.g. .../Input/t57.in
        t_stdin.write(self.id + ".out\n")     # Output. e.g t57.out
        t_stdin.write(self.id + "\n")         # Used as suffix to diff and prefix to log file names,
                                               # and also for roots for temporaries

        return t_stdin.getvalue()


class Band2epsTest(BaseTest):
    """How to waste lines of code just to test a F90 code that can be implemented with a few python commands!"""

    def make_stdin(self):
        t_stdin = StringIO()

        t_stdin.write(self.inp_fname + "\n")   # input file e.g. .../Input/t51.in
        t_stdin.write(self.id + ".out.eps\n")  # output file e.g. t51.out.eps

        inp_freq = os.path.join(self.inp_dir, self.id + ".in_freq")
        t_stdin.write(inp_freq + "\n")          # input freq file e.g Input/t51.in_freq

        inp_displ = os.path.join(self.inp_dir, self.id + ".in_displ")
        if not os.path.isfile(inp_displ): inp_displ = "no"
        t_stdin.write(inp_displ + "\n")         # input displ file e.g Input/t51.in_displ

        return t_stdin.getvalue()


class AtompawTest(BaseTest):
    """
    Class for Atompaw tests. Redefine the methods clean_workdir and bin_path provided by BaseTest
    """
    def clean_workdir(self, other_test_files=None):
        """Keep all atompaw output files."""

    @property
    def bin_path(self):
        """atompaw is not located in src/98_main"""
        return self.build_env.path_of_bin("atompaw")


def exec2class(exec_name):
    """
    Return the test class associated to the executable. Default is BaseTest.
    """
    return {
        "abinit": AbinitTest,
        "anaddb": AnaddbTest,
        "aim": AimTest,
        "conducti": ConductiTest,
        "atompaw": AtompawTest,
        "band2eps": Band2epsTest,
        "optic": OpticTest,
        "multibinit": MultibinitTest,
        "tdep": TdepTest,
    }.get(exec_name, BaseTest)


class ChainOfTests(object):
    """
    A list of tests that should be executed together due to inter-dependencies.
    It provides the same interface as the one given by BaseTest
    """
    Error = BaseTestError

    def __init__(self, tests):
        self.tests = tuple(t for t in tests)

        self.inp_dir = tests[0].inp_dir
        self.suite_name = tests[0].suite_name

        # Consistency check.
        self._rid = genid()
        for t in tests:
            if self.inp_dir != t.inp_dir or self.suite_name != t.suite_name:
                raise self.Error("All tests should be located in the same directory")
            self._rid += ':' + t._rid

        all_keys = [t.keywords for t in self.tests]
        self.keywords = set()
        for ks in all_keys:
            self.keywords = self.keywords.union(ks)

        all_cpp_vars = [t.need_cpp_vars for t in self.tests]
        self.need_cpp_vars = set()
        for vs in all_cpp_vars:
            self.need_cpp_vars = self.need_cpp_vars.union(vs)

        self._priv_executed = None
        self._status = None
        self._tot_etime = None
        self._run_etime = None
        self._isok = None
        self._files_to_keep = []

    def __len__(self):
        return len(self.tests)

    def __str__(self):
        return "\n".join(str(t) for t in self)

    def __iter__(self):
        for t in self.tests:
            yield t

    def info_on_chain(self):
        attr_names = ["extra_inputs", "pre_commands", "post_commands"]
        string = "Info on chain: %s\n" % self.full_id

        nlinks = 0
        for test in self:
            string += test.full_id + "executable " + test.executable + ":\n"
            for (attr, value) in test.__dict__.items():
                if (value and (attr in attr_names or attr.startswith("input_")
                               or attr.startswith("output_"))):
                    string += "  %s = %s\n" % (attr, value)
                    nlinks += 1

        return string, nlinks

    # A lot of boilerplate code!
    # See the doc strings of BaseTest
    @property
    def id(self):
        return "-".join(test.id for test in self)

    @property
    def full_id(self):
        return "[{}][{}]".format(self.suite_name, self.id)

    @property
    def max_nprocs(self):
        return max(test.max_nprocs for test in self)

    @property
    def _executed(self):
        if self._priv_executed is None:
            self._priv_executed = all(test._executed for test in self)
        return self._priv_executed

    @property
    def ref_dir(self):
        ref_dirs = [test.ref_dir for test in self]
        assert all(d == ref_dirs[0] for d in ref_dirs)
        return ref_dirs[0]

    def listoftests(self, width=100, html=True, abslink=True):
        string = ""
        if not html:
            string += "\n".join(test.listoftests(width, html, abslink) for test in self)
            string = self.full_id + ":\n" + string
        else:
            string += "<br>".join(test.listoftests(width, html, abslink) for test in self)
            string = "Test Chain " + self.full_id + ":<br>" + string
        return string

    @property
    def files_to_test(self):
        files = []
        for test in self:
            files.extend(test.files_to_test)
        return files

    @property
    def extra_inputs(self):
        extra_inputs = []
        for test in self:
            extra_inputs.extend(test.extra_inputs)
        return extra_inputs

    @property
    def inputs_used(self):
        inputs = []
        for test in self:
            inputs.extend(test.inputs_used)
        return inputs

    @property
    def run_etime(self):
        if self._run_etime is None:
            self._run_etime = sum(test.run_etime for test in self)
        return self._run_etime

    @property
    def tot_etime(self):
        if self._tot_etime is None:
            self._tot_etime = sum(test.tot_etime for test in self)
        return self._tot_etime

    @property
    def isok(self):
        if self._isok is None:
            self._isok = all(test.isok for test in self)
        return self._isok

    @property
    def exceptions(self):
        excs = []
        for test in self:
            excs.extend(test.exceptions)

        return excs

    @property
    def status(self):
        if self._status is None:
            _stats = {test.status for test in self}
            if "disabled" in _stats or "skipped" in _stats:
                if len(_stats) > 1:  # it must be {'skipped'} or {'disabled'}
                    self._status = 'failed'
                else:
                    self._status = _stats.pop()
            else:
                all_fldstats = {f.fld_status for f in self.files_to_test}

                if "failed" in all_fldstats:
                    self._status = "failed"
                elif "passed" in all_fldstats:
                    self._status = "passed"
                elif all_fldstats != {"succeeded"}:
                    print(self)
                    print("WARNING, expecting {'succeeded'} but got\n%s" % str(all_fldstats))
                    self._status = "failed"
                else:
                    self._status = "succeeded"

        return self._status

    def keep_files(self, files):
        if is_string(files):
            self._files_to_keep.append(files)
        else:
            self._files_to_keep.extend(files)

    @property
    def files_to_keep(self):
        # The files produced by the individual tests.
        files_of_tests = []
        for test in self:
            files_of_tests.extend(test.files_to_keep)

        # Add the files produced by self.
        self._files_to_keep += files_of_tests
        return self._files_to_keep

    def cpkl_dump(self, protocol=-1):
        self.cpkl_fname = os.path.join(self.workdir, self.id + ".cpkl")
        with open(self.cpkl_fname, "wb") as fh:
            pickle.dump(self, fh, protocol=protocol)
            self.files_to_keep.append(self.cpkl_fname)

    def has_keywords(self, keywords, mode="any"):
        if mode == "all":
            return set(keywords).issubset(self.keywords)
        elif mode == "any":
            return set(keywords).intersection(self.keywords)
        else:
            raise ValueError("wrong mode %s" % mode)

    def has_variables(self, ivars):
        for test in self:
            matches = test.has_variables(ivars)
            if matches:
                return matches

        return []

    def edit_input(self, editor=None):
        if editor is None:
            editor = Editor()

        for test in self:
            try:
                test.edit_input(editor=editor)
            except Exception as e:
                raise e

    @property
    def _authors_snames(self):
        snames = set()
        for test in self:
            snames = snames.union(test._authors_snames)
        return snames

    def has_authors(self, authors, mode="any"):
        # return set(authors).issubset(self._authors_snames)
        if mode == "all":
            return set(authors).issubset(self._authors_snames)
        elif mode == "any":
            return set(authors).intersection(self._authors_snames)
        else:
            raise ValueError("wrong mode %s" % mode)

    def write_html_report(self):
        html_report = os.path.join(self.workdir, "test_report.html")
        with open(html_report, "wt") as fh:
            for idx, test in enumerate(self):
                oc = ""
                if idx == 0:
                    oc += "o"
                if idx == (len(self) - 1):
                    oc += "c"
                test.write_html_report(fh=fh, oc=oc)

    def run(self, build_env, runner, workdir, nprocs=1, **kwargs):

        workdir = os.path.abspath(workdir)
        if not os.path.exists(workdir):
            os.mkdir(workdir)
        self.workdir = workdir

        fail_all = False

        for test in self:
            if fail_all:
                test.force_skip = True
            test.run(build_env, runner, workdir=self.workdir, nprocs=nprocs, **kwargs)
            if test.had_timeout:
                fail_all = True

    def results_load(self, d):
        """
        Load the run results from a run in a different process.
        """
        self._status = d['status']
        self._files_to_keep = d['files_to_keep']
        self._tot_etime = d['tot_etime']
        self._run_etime = d['run_etime']
        self._priv_executed = d['executed']
        self._isok = d['isok']
        self.workdir = d['workdir']

    def results_dump(self):
        """
        Dump the run results to pass it to a different process
        """
        return {
            'id': self._rid,
            'status': self.status,
            'files_to_keep': self.files_to_keep,
            'tot_etime': self.tot_etime,
            'run_etime': self.run_etime,
            'executed': self._executed,
            'isok': self.isok,
            'workdir': self.workdir
        }

    def clean_workdir(self, other_test_files=None):
        for test in self:
            test.clean_workdir(other_test_files=self.files_to_keep)

    def patch(self, patcher=None):
        for test in self:
            test.patch(patcher)

    def get_backtraces(self):
        return [test._get_one_backtrace() for test in self]


class AbinitTestSuite(object):
    """
    List of BaseTest instances. Provide methods to:

    1) select subset of tests according to keywords, authors, numbers
    2) run tests in parallel with python processes
    3) analyze the final results
    """
    def __init__(self, abenv, inp_files=None, test_list=None, keywords=None, need_cpp_vars=None):

        # One and only one should be provided
        if (inp_files is None) == (test_list is None):
            raise ValueError(
                "One and only one of inp_file and test_list is expected but"
                " found: {} and {}".format(inp_files, test_list)
            )

        self._executed = False
        self._kill_me = False
        self.abenv = abenv
        self.exceptions = []
        self._processes = []

        if inp_files is not None:
            self.tests = make_abitests_from_inputs(
                inp_files, abenv,
                keywords=keywords, need_cpp_vars=need_cpp_vars
            )

        elif test_list is not None:
            assert keywords is None, ("keywords argument is not expected with test_list")
            assert need_cpp_vars is None, ("need_cpp_vars argument is not expected with test_list.")
            self.tests = tuple(test_list)

    def __str__(self):
        return "\n".join(str(t) for t in self.tests)

    def __add__(self, other):
        test_list = [t for t in self] + [t for t in other]
        return self.__class__(self.abenv, test_list=test_list)

    def __len__(self):
        return len(self.tests)

    def __iter__(self):
        for t in self.tests:
            yield t

    def __getitem__(self, key):
        """Called by self[key]."""
        # FIXME: this won't work for tutorial, paral and other test suites.
        if isinstance(key, slice):
            return self.__getslice(key)
        else:
            raise NotImplementedError("__getitem__ expects a slice instance")

    def __getslice(self, slice):
        start = slice.start
        if start is None:
            start = 0
        stop = slice.stop
        if stop is None:
            stop = 10000  # Not very elegant, but cannot use len(self) since indices are not contiguous
        assert slice.step is None, ("Slices with steps (e.g. [1:4:2]) are not  supported.")

        # Rules for the test id:
        # Simple case: t01, tgw1_1
        # test chain (no MPI): t81-t82-t83-t84, tudet_1-tudet_2-tudet_3
        # multi-parallel tests:  t74_MPI2, t51_MPI1-t52_MPI1-t53_MPI1, tdfpt_01_MPI2 ...

        test_list = []
        for test in self:
            # print("ID",test.id)
            # extract the ID of the first test (if test_chain)
            tokens = test.id.split("-")
            assert tokens[0][0] == "t"  # Assume first character is "t"
            num = tokens[0][1:]

            if "_MPI" in test.id:
                # Handle multi-parallel tests.
                # print(test.id)
                idx = test.id.find("_MPI")
                tok = test.id[1:idx]
                # print(tok)
                idx = tok.rfind("_")
                if idx != -1:
                    # Handle tdfpt_01_MPI2 ...
                    # FIXME: this will fail if _OMP2_MPI2
                    tok = tok[idx + 1:]
                try:
                    num = int(tok)
                except ValueError:
                    raise ValueError("Cannot convert %s to integer" % tok)

            else:
                # Simple case or test_chain
                idx = num.rfind("_")
                if idx != -1:
                    num = int(num[idx + 1:])

            num = int(num)
            if num in range(start, stop):
                # print "got", test.id
                test_list.append(test)

        return self.__class__(self.abenv, test_list=test_list)

    @property
    def full_length(self):
        return sum(getattr(test, "__len__", lambda: 1)() for test in self)

    @property
    def run_etime(self):
        """Total elapsed time i.e. the wall-time spent in the sub-processes e.g. abinit)"""
        assert self._executed
        return sum(test.run_etime for test in self)

    @property
    def keywords(self):
        keys = []
        for test in self:
            keys.extend(test.keywords)
        return set(keys)

    def has_keywords(self, keywords):
        return set(keywords).issubset(self.keywords)

    @property
    def need_cpp_vars(self):
        cpp_vars = []
        for test in self:
            cpp_vars.extend(test.need_cpp_vars)
        return set(cpp_vars)

    def on_refslave(self):
        """True if we are running on a reference slave e.g. abiref."""
        try:
            return self._on_ref_slave
        except AttributeError:
            return False

    def set_on_refslave(self, value=True):
        """Attribute setter"""
        self._on_ref_slave = bool(value)

    def all_exceptions(self):
        """Return my exceptions + test exceptions."""
        all_excs = self.exceptions
        for test in self:
            all_excs.extend(test.exceptions)
        return all_excs

    def cpkl_dump(self, protocol=-1):
        self.cpkl_fname = os.path.join(self.workdir, "test_suite.cpkl")
        with open(self.cpkl_fname, "wb") as fh:
            pickle.dump(self, fh, protocol=protocol)

    def _tests_with_status(self, status):
        assert status in BaseTest._possible_status
        # assert self._executed
        return [test for test in self if test.status == status]

    def succeeded_tests(self): return self._tests_with_status("succeeded")
    def passed_tests(self):    return self._tests_with_status("passed")
    def failed_tests(self):    return self._tests_with_status("failed")
    def skipped_tests(self):   return self._tests_with_status("skipped")
    def disabled_tests(self):  return self._tests_with_status("disabled")

    @property
    def targz_fname(self):
        """
        Location of the tarball file with the results in HTML format
        Returns None if the tarball has not been created.
        """
        return getattr(self, '_targz_fname', None)

    def create_targz_results(self):
        """Create the tarball file results.tar.gz in the working directory."""
        assert self._executed
        exclude_exts = [".cpkl", ".py", "pyc", ]

        self._targz_fname = None
        ofname = os.path.join(self.workdir, "results.tar.gz")

        # The most delicate part here is the treatment of the exceptions
        # since the test might not have produced the reference files
        # we want to copy to the server. If something goes wrong, we simply
        # register the exception and we continue the execution.
        try:
            targz = tarfile.open(ofname, "w:gz")

            for test in self:
                # Don't try to collect files for tests that are disabled or skipped.
                if test.status in {"disabled", "skipped"}:
                    continue

                files = set(test.files_to_keep)
                save_files = {f for f in files if not has_exts(f, exclude_exts)}
                # print(save_files)

                # Store stdout files only if the test failed.
                important_status = {"failed", }

                # Special treatment for reference machines
                if self.on_refslave:
                    important_status = {"passed", "failed", }

                if test.status not in important_status:
                    if isinstance(test, ChainOfTests):
                        for t in test:
                            # print "Removing Test Chain", t.stdout_fname
                            save_files.discard(t.stdout_fname)
                    else:
                        # print "Removing", test.stdout_fname
                        save_files.discard(test.stdout_fname)

                for p in save_files:
                    # if not os.path.exists(os.path.join(self.workdir, p)): continue
                    # /foo/bar/suite_workdir/test_workdir/file --> test_workdir/t01/file
                    rpath = os.path.relpath(p, start=self.workdir)
                    # arcname = str(rpath.encode("ascii", "ignore"))
                    arcname = str(rpath)
                    try:
                        # print("targz.add: adding:", p," with arcname ", arcname)
                        # print(type(p), type(arcname))
                        targz.add(p, arcname=arcname)
                    except Exception as exc:
                        # Handle the case in which the output file has not been produced.
                        warnings.warn("exception while adding %s to tarball:\n%s" % (p, exc))
                        self.exceptions.append(exc)

            targz.close()

            # Save the name of the tarball file.
            self._targz_fname = ofname

        except Exception as exc:
            warnings.warn("exception while creating tarball file: %s" % str(exc))
            self.exceptions.append(exc)

    def sanity_check(self):
        all_full_ids = [test.full_id for test in self]
        if len(all_full_ids) != len(set(all_full_ids)):
            raise ValueError("Cannot have more than two tests with the same full_id")

    def start_workers(self, nprocs, runner):
        """
        Start nprocs new processes that will get tests from a queue and run
        them with runner and put the result of runner in a output queue.
        Return the task/input queue (to be closed only) and the results/output queue.
        """

        def worker(qin, qout, print_lock, thread_mode=False):
            done = {'type': 'proc_done'}
            all_done = False
            try:
                while not all_done and not (thread_mode and self._kill_me):
                    test = qin.get(block=True, timeout=2)
                    if test is None:  # reached the end
                        all_done = True
                    else:
                        qout.put(runner(test, print_lock=print_lock))
            except EmptyQueueError:
                # If that happen it is a probably a bug
                done['error'] = RuntimeError('Task queue is unexpectedly empty.')
            except Exception as e:
                # Any other error is reported
                done['error'] = e
                try:
                    done['task'] = test.full_id
                except (AttributeError, NameError):
                    pass
            finally:
                qout.put(done)

        print_lock = Lock()
        task_q = Queue()
        res_q = Queue()

        for test in self:
            # fill the queue
            task_q.put(test)

        for _ in range(nprocs):
            # one end signal for each worker
            task_q.put(None)

        for i in range(nprocs - 1):
            # create and start subprocesses
            p = Process(target=worker, args=(task_q, res_q, print_lock))
            self._processes.append(p)
            p.start()

        # Add the worker as a thread of the main process
        t = Thread(target=worker, args=(task_q, res_q, print_lock, True))
        # make it daemon so it will die if the main process is interupted early
        t.daemon = True

        t.start()

        return task_q, res_q

    def wait_loop(self, nprocs, ntasks, timeout, queue):
        '''
        Wait for all tests to be done by workers. Receives tests results from
        queue and update the local tests objects.
        '''
        results = {}
        proc_running, task_remaining = nprocs, ntasks
        try:
            while proc_running > 0:
                msg = queue.get(block=True, timeout=(
                    1 + 2 * task_remaining * timeout / proc_running
                ))
                if msg['type'] == 'proc_done':
                    proc_running -= 1
                    if 'error' in msg:
                        e = msg['error']
                        if 'task' in msg:
                            task_remaining -= 1
                            warnings.warn(
                                'Error append in a worker on test {}:\n{}: {}'
                                .format(msg['task'], type(e).__name__, e)
                            )
                        else:
                            warnings.warn(
                                'Error append in a worker:\n{}: {}'
                                .format(type(e).__name__, e)
                            )

                    logger.info("{} worker(s) remaining for {} tasks."
                                .format(proc_running, task_remaining))
                elif msg['type'] == 'result':
                    results[msg['id']] = msg
                    task_remaining -= 1

        except KeyboardInterrupt:
            self.terminate()
            raise KeyboardInterrupt()

        except EmptyQueueError:
            warnings.warn(
                ("Workers have been hanging until timeout. There were {} procs"
                 " working on {} tasks.").format(proc_running, task_remaining)
            )
            self.terminate()
            return None

        return results

    def run_tests(self, build_env, workdir, runner, nprocs=1, py_nprocs=1,
                  runmode="static", **kwargs):
        """
        Execute the list of tests (main entry point for client code)

        Args:
            build_env: `BuildEnv` instance with info on the build environment.
            workdir: Working directory (string)
            runner: `JobRunner` instance
            nprocs: number of MPI processes to use for a single test.
            py_nprocs: number of py_nprocs for tests
        """
        self.sanity_check()

        if len(self) == 0:
            warnings.warn("len(self) == 0")
            return

        workdir = os.path.abspath(workdir)
        if not os.path.exists(workdir):
            os.mkdir(workdir)

        self.workdir = workdir

        # Acquire the lock file.
        self.lock = NoErrorFileLock(os.path.join(workdir, "__run_tests_lock__"),
                                    timeout=3)

        with self.lock as locked:
            # aquire the global file lock
            if not locked:
                msg = (
                    "Timeout occured while trying to acquire lock in:\n\t{}\n"
                    "Perhaps a previous run did not exit cleanly or another "
                    "process is running in the same directory.\n If you are"
                    "sure no other process is in execution, remove the "
                    "directory with `rm -rf` and rerun.\n"
                ).format(self.workdir)

                cprint(msg, "red")
                return

            # Remove all stale files present in workdir (except the lock!)
            rmrf(self.workdir, exclude_paths=self.lock.lockfile)

            self.nprocs = nprocs
            self.py_nprocs = py_nprocs

            def run_and_check_test(test, print_lock=None):
                """Helper function to execute the test. Must be thread-safe."""

                testdir = os.path.abspath(os.path.join(self.workdir, test.suite_name + "_" + test.id))

                # Run the test
                test.run(build_env, runner, testdir, print_lock=print_lock,
                         nprocs=nprocs, runmode=runmode, **kwargs)

                # Write HTML summary
                test.write_html_report()

                # Remove useless files in workdir.
                test.clean_workdir()

                d = test.results_dump()
                d['type'] = 'result'
                return d

            ##############################
            # And now let's run the tests
            ##############################
            start_time = time.time()

            if py_nprocs == 1:
                logger.info("Sequential version")
                for test in self:
                    # discard the return value because tests are directly modified
                    run_and_check_test(test)

            elif py_nprocs > 1:
                logger.info("Parallel version with py_nprocs = %s" % py_nprocs)

                task_q, res_q = self.start_workers(py_nprocs, run_and_check_test)

                timeout_1test = float(runner.timebomb.timeout)
                if timeout_1test <= 0.1:
                    timeout_1test = 240.

                # Wait for all tests to be done gathering results
                results = self.wait_loop(py_nprocs, len(self.tests),
                                         timeout_1test, res_q)

                # remove this to let python garbage collect processes and avoid
                # Pickle to complain (it does not accept processes for security reasons)
                self._processes = []
                task_q.close()
                res_q.close()

                # update local tests instances with the results of their running in
                # a remote process
                for test in self.tests:
                    if test._rid not in results:
                        # This error will only happen if there is a bug
                        raise RuntimeError((
                            "I did not get the results of the test {}. It"
                            " means that something fishy happen in the worker."
                        ).format(test.full_id))
                    test.results_load(results[test._rid])

            # Run completed.
            self._executed = True

            # Collect HTML files in a tarball
            self.create_targz_results()

            nsucc = len(self.succeeded_tests())
            npass = len(self.passed_tests())
            nfail = len(self.failed_tests())
            nskip = len(self.skipped_tests())
            ndisa = len(self.disabled_tests())

            self.tot_etime = time.time() - start_time

            # Print summary table.
            stats_suite = {}
            for test in self:
                if test.suite_name not in stats_suite:
                    d = dict.fromkeys(BaseTest._possible_status, 0)
                    d["run_etime"] = 0.0
                    d["tot_etime"] = 0.0
                    stats_suite[test.suite_name] = d

                stats_suite[test.suite_name][test.status] += 1
                stats_suite[test.suite_name]["run_etime"] += test.run_etime
                stats_suite[test.suite_name]["tot_etime"] += test.tot_etime

            suite_names = sorted(stats_suite.keys())

            times = ["run_etime", "tot_etime"]

            table = [["Suite"] + BaseTest._possible_status + times]
            for suite_name in suite_names:
                stats = stats_suite[suite_name]
                row = [suite_name] + [str(stats[s]) for s in BaseTest._possible_status] + ["%.2f" % stats[s] for s in times]
                table.append(row)

            print("")
            pprint_table(table)
            print("")

            executed = [t for t in self if t.status != "skipped"]
            if executed:
                mean_etime = sum(test.run_etime for test in executed) / len(executed)
                dev_etime = (sum((test.run_etime - mean_etime)**2 for test in executed) / len(executed))**0.5

                cprint("Completed in %.2f [s]. Average time for test=%.2f [s], stdev=%.2f [s]" % (
                    self.tot_etime, mean_etime, dev_etime), "yellow"
                )

                msg = "Summary: failed=%s, succeeded=%s, passed=%s, skipped=%s, disabled=%s" % (
                    nfail, nsucc, npass, nskip, ndisa)

                if nfail:
                    cprint(msg, "red", attrs=['underline'])
                else:
                    cprint(msg, "green")

                # Print outliers
                if False and dev_etime > 0.0:
                    for test in self:
                        if abs(test.run_etime) > 0.0 and abs(test.run_etime - mean_etime) > 2 * dev_etime:
                            print("%s has run_etime %.2f s" % (test.full_id, test.run_etime))

            with open(os.path.join(self.workdir, "results.txt"), "wt") as fh:
                pprint_table(table, out=fh)

            username = my_getlogin()

            # Create the HTML index.
            DNS = {
                "self": self,
                "runner": runner,
                "user_name": username,
                "hostname": gethostname(),
                "test_headings": ['ID', 'Status', 'run_etime (s)', 'tot_etime (s)'],
                "suite_headings": ['failed', 'passed', 'succeeded', 'skipped', 'disabled'],
                # Functions and modules available in the template.
                "time": time,
                "pj": os.path.join,
                "basename": os.path.basename,
                "str2html": str2html,
                "sec2str": sec2str,
                "args2htmltr": args2htmltr,
                "html_link": html_link,
                "status2html": status2html,
            }

            fname = os.path.join(self.workdir, "suite_report.html")
            fh = open(fname, "w")

            header = """
            <html>
            <head><title>Suite Summary</title></head>
            <body bgcolor="#FFFFFF" text="#000000">
            <hr>
            <h1>Suite Summary</h1>
                <table width="100%" border="0" cellspacing="0" cellpadding="2">
                <tr valign="top" align="left">
                <py-open code = "for h in suite_headings:"> </py-open>
                <th>${status2html(h)}</th>
                <py-close/>
                </tr>
                <tr valign="top" align="left">
                <py-open code = "for h in suite_headings:"> </py-open>
                <td> ${len(self._tests_with_status(h))} </td>
                <py-close/>
                </tr>
                </table>
                <p>
                tot_etime = ${sec2str(self.tot_etime)} <br>
                run_etime = ${sec2str(self.run_etime)} <br>
                no_pyprocs = ${self.py_nprocs} <br>
                no_MPI = ${self.nprocs} <br>
                ${str2html(str(runner))}
            <hr>
            """

            table = """
            <p>
            <h1>Test Results</h1>
            <table width="100%" border="0" cellspacing="0" cellpadding="2">
                <tr valign="top" align="left">
                <py-open code = "for h in test_headings:"> </py-open>
                <th>$h</th>
                <py-close/>
                </tr>
            """

            for status in BaseTest._possible_status:
                table += self._pyhtml_table_section(status)

            table += "</table>"

            footer = """
            <hr>
            <h1>Suite Info</h1>
                <py-line code = "keys = ', '.join(self.keywords)" />
                <p>Keywords = ${keys}</p>
                <py-line code = "cpp_vars = ', '.join(self.need_cpp_vars)"/>
                <p>Required CPP variables = ${cpp_vars}</p>
            <hr>
                Automatically generated by %s on %s. Logged on as %s@%s
            <hr>
            </body>
            </html> """ % (_MY_NAME, time.asctime(), username, gethostname())

            template = header + table + footer

            template_stream = StringIO(template)

            # Initialise an xyaptu xcopier, and call xcopy
            xcp = xcopier(DNS, ouf=fh)
            xcp.xcopy(template_stream)
            fh.close()

        return Results(self)

    def terminate(self):
        '''
        Kill all workers
        '''
        for p in self._processes:
            p.terminate()
        self._kill_me = True
        self._processes = []

    @staticmethod
    def _pyhtml_table_section(status):
        # ['ID', 'Status', 'run_etime', 'tot_etime'],
        string = """
           <py-open code="for test in self.%s_tests():"/>
            <py-line code = "report_link = pj(basename(test.workdir),'test_report.html') " />
            <tr valign="top" align="left">
             <td> ${html_link(test.full_id, report_link)}</td>
             <td> ${status2html(test.status)} </td>
             <td> ${sec2str(test.run_etime)} </td>
             <td> ${sec2str(test.tot_etime)} </td>
            </tr>
           <py-close/>
           """ % status
        return string

    def patch(self, patcher=None):
        """
        Patch the output files of the test with the specified patcher.
        A default patcher is provided if patcher is None (use $PATCHER shell variable)
        """
        for test in self:
            test.patch(patcher)

    def select_tests(self, with_keys=None, exclude_keys=None, with_authors=None, exclude_authors=None,
                     ivars=None, mode="any"):
        """
        Extract the subset of tests matching the given conditions.

        Returns:
            `AbinitTestSuite` instance
        """
        test_list = [test for test in self]

        if with_keys:
            test_list = [test for test in test_list if test.has_keywords(with_keys, mode=mode)]

        if exclude_keys:
            test_list = [test for test in test_list if not test.has_keywords(exclude_keys, mode=mode)]

        if with_authors:
            test_list = [test for test in test_list if test.has_authors(with_authors, mode=mode)]

        if exclude_authors:
            test_list = [test for test in test_list if not test.has_authors(exclude_authors, mode=mode)]

        if ivars:
            test_list = [test for test in test_list if test.has_variables(ivars)]

        return AbinitTestSuite(self.abenv, test_list=test_list)

    def make_listoftests(self, width=100, html=True):
        """Create the ListOfTests files."""
        if not html:
            return "\n\n".join(test.listoftests(width, html) for test in self)
        else:
            header = """
             <html>
             <head><title>"LIST OF TESTS" FILE</title></head>
             <body bgcolor="#FFFFFF" text="#000000">
             <!-- Automatically generated by %s on %s. ****DO NOT EDIT**** -->""" % (_MY_NAME, time.asctime())

            body = "<hr>".join(test.listoftests(width, html) for test in self)

            footer = """
              <hr>
               Automatically generated by %s on %s.
              <hr>
              </body>
              </html>""" % (_MY_NAME, time.asctime())

            return header + body + footer


class Results(object):
    """Stores the final results."""
    def __init__(self, test_suite):
        # assert test_suite._executed
        self.test_suite = test_suite
        self.failed_tests = test_suite.failed_tests()
        self.passed_tests = test_suite.passed_tests()
        self.succeeded_tests = test_suite.succeeded_tests()
        self.skipped_tests = test_suite.skipped_tests()
        self.disabled_tests = test_suite.disabled_tests()
        self.targz_fname = test_suite.targz_fname

    @lazy__str__
    def __str__(self): pass

    def tests_with_status(self, status):
        return {
            "succeeded": self.succeeded_tests,
            "passed": self.passed_tests,
            "failed": self.failed_tests,
            "disabled": self.disabled_tests,
            "skipped": self.skipped_tests,
            "all": [test for test in self.test_suite]
        }[status]

    @property
    def nfailed(self):
        """Number of failures"""
        return len(self.failed_tests)

    @property
    def npassed(self):
        """Number of tests marked as passed."""
        return len(self.passed_tests)

    @property
    def nexecuted(self):
        """Number of tests executed."""
        n = 0
        for test in self.test_suite:
            if isinstance(test, ChainOfTests):
                n += len([t for t in test if t.executed])
            else:
                if test.executed:
                    n += 1
        return n

    def outref_files(self, status):
        """
        Return (out_files, ref_files)
        where out and ref are list with the output files and the reference
        files of the tests with the given status.
        """
        out_files, ref_files = [], []
        for test in (self.tests_with_status(status)):
            for f in test.files_to_test:
                # if status != "all" and f.status != status: continue

                out_files.append(os.path.join(test.workdir, f.name))
                ref_fname = os.path.join(test.ref_dir, f.name)
                # FIXME Hack due to the ambiguity stdout, out!
                if not os.path.exists(ref_fname) and ref_fname.endswith(".stdout"):
                    ref_fname = ref_fname[:-7] + ".out"
                ref_files.append(ref_fname)

        return out_files, ref_files

    def in_files(self, status):
        """List with the input files of the tests with the given status."""
        in_files = []
        for test in (self.tests_with_status(status)):
            if isinstance(test, ChainOfTests):
                in_files.extend(t.inp_fname for t in test)
            else:
                in_files.append(test.inp_fname)

        return in_files

    def patch_refs(self, status="failed"):
        """Patch the reference files of the tests with the specified status."""
        out_files, ref_files = self.outref_files(status=status)
        # for r, o in zip(out_files, ref_files): print("reference: %s, output %s" % (r, o))

        return Patcher().patch_files(out_files, ref_files)

    def edit_inputs(self, status="failed"):
        """Edit the input files of the tests with the specified status."""
        in_files = self.in_files(status=status)
        # for r, o in zip(out_files, ref_files):
        #     print("reference: %s, output %s" % (r, o))

        return Editor().edit_files(in_files)

    def inspect_stderrs(self, status="failed"):
        """Open the stderr of the tests with the give status in `Editor`."""
        return Editor().edit_files(self.stderr_files(status))

    def stderr_files(self, status="failed"):
        """List of non-empty error files of the tests with the specified status."""
        # Loop over the tests, open the stderr to see if it's empty ot not
        # and add it to the list.
        err_files = []
        for test in self.tests_with_status(status):
            if isinstance(test, ChainOfTests):
                es = [t.stderr_fname for t in test if not t.has_empty_stderr]
                if es:
                    err_files.extend(es)
            else:
                if not test.has_empty_stderr:
                    err_files.append(test.stderr_fname)

        return err_files

    def cpkl_dump(self, cpkl_fname, protocol=-1):
        """Save the object in pickle format."""
        with open(cpkl_fname, "wb") as fh:
            pickle.dump(self, fh, protocol=protocol)


if __name__ == "__main__":
    # Automatic documentation of the TEST_INFO options.
    doc_testcnf_format()
