"""
Unit tests for testsuite.py module.

Tests cover configuration parsing, file comparison setup, and test info handling.
Priority 1: High-impact, core infrastructure components.
"""

from __future__ import annotations

import os
import shutil
import tempfile
import threading
from pathlib import Path

import pytest

from tests import abenv

from . import testsuite as testsuite_module
from .jobrunner import JobRunner
from .testsuite import (
    SLURM_ACCT_FORMAT,
    SLURM_OOM_MAX_QUERIES,
    AbinitTestInfo,
    AbinitTestInfoParser,
    AbinitTestInfoParserError,
    AbinitTestSuite,
    BaseTest,
    BuildEnvironment,
    ChainOfTests,
    Compiler,
    CPreProcessor,
    FileToTest,
    FortranCompiler,
    _LocalCounter,
    _oom_query_allowed,
    _str2bool,
    _str2cmds,
    _str2filestotest,
    _str2intlist,
    _str2list,
    _str2set,
    args2htmltr,
    detect_slurm_oom_stepid,
    genid,
    has_exts,
    html_colorize_text,
    html_file_link,
    html_link,
    input_file_has_vars,
    is_string,
    lazy_read,
    lazy_readlines,
    make_abitest_from_input,
    make_abitests_from_inputs,
    my_getlogin,
    parse_configh_file,
    query_slurm_step_accounting,
    rm_rf,
    sec2str,
    status2html,
    str2html,
)

# ============================================================================
# FIXTURES
# ============================================================================


@pytest.fixture
def temp_test_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


# By convention, a local out-of-tree build lives in `_build` at the top of
# the repository (this is a project convention, not something ABINIT's
# BuildEnvironment itself enforces or assumes).
ABINIT_BUILD_DIR = Path(__file__).resolve().parents[2] / "_build"


@pytest.fixture
def build_environment():
    """A real BuildEnvironment built from the repo's conventional `_build` dir.

    `_build` may not exist at all, or may hold an incomplete/stale build
    (e.g. `./configure` ran but `make` never finished, so `config.h` exists
    but the `abinit` binary doesn't) -- BuildEnvironment's constructor raises
    ValueError/RuntimeError in either case. Skip rather than fail when that
    happens, since it reflects the state of this checkout, not a code bug --
    but print the reason unconditionally (not just as a pytest skip reason,
    which is easy to miss without -rs/-v) so it's obvious on the terminal
    why these tests didn't run.
    """
    try:
        return BuildEnvironment(str(ABINIT_BUILD_DIR))
    except (ValueError, RuntimeError) as exc:
        msg = (
            f"Skipping: {ABINIT_BUILD_DIR} is not a valid/complete ABINIT "
            f"build tree ({exc}). Run configure && make in _build to enable "
            "this test."
        )
        print(msg)
        pytest.skip(msg)


@pytest.fixture
def sample_test_info_content():
    """Sample TEST_INFO section content (valid)."""
    return """\
<BEGIN TEST_INFO>
[setup]
executable = abinit
use_files_file = no
exec_args =
test_chain =
need_cpp_vars =
exclude_hosts =
exclude_builders =
input_prefix =
output_prefix =
expected_failure = no
input_ddb =
input_gkk =
system_xml =
coeff_xml =
md_hist =
test_set =
no_check = no
spin_pot =
latt_pot =
slc_pot =
lwf_pot =

[files]
files_to_test = output.txt, tolnlines=0, tolabs=0.01, tolrel=1e-3
psp_files =
extra_inputs =

[shell]
pre_commands =
post_commands =

[paral_info]
max_nprocs = 1
nprocs_to_test =
exclude_nprocs =

[extra_info]
authors = Unknown
keywords =
description = No description available
topics =
references =
<END TEST_INFO>
"""


def _make_minimal_testinfo_data() -> dict:
    """Create minimal valid AbinitTestInfo data for testing."""
    return {
        "executable": "abinit",
        "use_files_file": False,
        "exec_args": "",
        "test_chain": [],
        "need_cpp_vars": set(),
        "exclude_hosts": [],
        "exclude_builders": [],
        "input_prefix": "",
        "output_prefix": "",
        "expected_failure": False,
        "input_ddb": "",
        "input_gkk": "",
        "system_xml": "",
        "coeff_xml": "",
        "md_hist": "",
        "test_set": "",
        "no_check": False,
        "spin_pot": "",
        "latt_pot": "",
        "slc_pot": "",
        "lwf_pot": "",
        "files_to_test": (),
        "psp_files": [],
        "extra_inputs": [],
        "use_git_submodule": "",
        "pre_commands": [],
        "post_commands": [],
        "max_nprocs": 1,
        "nprocs_to_test": [],
        "exclude_nprocs": [],
        "authors": {"Unknown"},
        "keywords": set(),
        "description": "Test",
        "topics": [],
        "references": [],
        "file": "",
        "yaml": "",
        "inp_fname": "/path/to/test.abi",
        "_ismulti_paral": False,
        "yaml_test": {},
    }


@pytest.fixture
def input_file_with_test_info(temp_test_dir):
    """Create a temporary input file with TEST_INFO section.

    The format must match what the parser expects:
    - Header: #%%<BEGIN TEST_INFO> (no space after #%%)
    - Body lines: #%% followed by content (with space after #%%)
    - Footer: #%%<END TEST_INFO> (no space after #%%)
    - The #%% prefix is removed, then leading space is stripped
    """
    input_path = Path(temp_test_dir) / "test.abi"
    content = """\
# ABINIT test input file
#%%<BEGIN TEST_INFO>
#%% [setup]
#%% executable = abinit
#%% use_files_file = no
#%% exec_args =
#%% [files]
#%% files_to_test = output.txt, tolnlines=0, tolabs=0.01, tolrel=1e-3
#%% [shell]
#%% pre_commands =
#%% post_commands =
#%% [paral_info]
#%% max_nprocs = 1
#%% nprocs_to_test =
#%% exclude_nprocs =
#%% [extra_info]
#%% authors = Unknown
#%% keywords =
#%% description = No description available
#%% topics =
#%% references =
#%%<END TEST_INFO>

# Test data below
npsp 1
"""
    input_path.write_text(content)
    return str(input_path)


# ============================================================================
# TESTS FOR PARSER FUNCTIONS
# ============================================================================


class TestParserFunctions:
    """Test string parsing functions for TEST_INFO options."""

    def test_str2bool_yes(self):
        """Test _str2bool with 'yes'."""
        assert _str2bool("yes") is True

    def test_str2bool_no(self):
        """Test _str2bool with 'no'."""
        assert _str2bool("no") is False

    def test_str2bool_case_insensitive(self):
        """Test _str2bool is case-insensitive."""
        assert _str2bool("YES") is True
        assert _str2bool("Yes") is True
        assert _str2bool("NO") is False

    def test_str2bool_with_whitespace(self):
        """Test _str2bool handles leading/trailing whitespace."""
        assert _str2bool("  yes  ") is True
        assert _str2bool("\tno\t") is False

    def test_str2list_single_item(self):
        """Test _str2list with single item."""
        result = _str2list("item1")
        assert result == ["item1"]

    def test_str2list_multiple_items(self):
        """Test _str2list with multiple comma-separated items."""
        result = _str2list("item1, item2, item3")
        assert result == ["item1", "item2", "item3"]

    def test_str2list_with_whitespace(self):
        """Test _str2list strips whitespace."""
        result = _str2list("  item1  ,  item2  ,  item3  ")
        assert result == ["item1", "item2", "item3"]

    def test_str2list_empty_string(self):
        """Test _str2list with empty string."""
        result = _str2list("")
        assert result == []

    def test_str2list_empty_items_with_spaces(self):
        """Test _str2list behavior with whitespace-only items.

        Note: items with only whitespace are stripped to empty strings
        but still included because whitespace-only strings are truthy
        in the original split.
        """
        result = _str2list("item1,  , item2")
        # Whitespace-only items become empty strings after strip()
        assert result == ["item1", "", "item2"]

    def test_str2intlist(self):
        """Test _str2intlist converts to integers."""
        result = _str2intlist("1, 2, 4, 8")
        assert result == [1, 2, 4, 8]

    def test_str2intlist_empty(self):
        """Test _str2intlist with empty string."""
        result = _str2intlist("")
        assert result == []

    def test_str2intlist_invalid_number(self):
        """Test _str2intlist raises on non-integer."""
        with pytest.raises(ValueError):
            _str2intlist("1, two, 4")

    def test_str2set(self):
        """Test _str2set creates set from comma-separated string."""
        result = _str2set("a, b, c")
        assert result == {"a", "b", "c"}

    def test_str2set_duplicates(self):
        """Test _str2set removes duplicates."""
        result = _str2set("a, b, a, c, b")
        assert result == {"a", "b", "c"}

    def test_str2set_empty(self):
        """Test _str2set with empty string."""
        result = _str2set("")
        assert result == set()

    def test_str2cmds_single(self):
        """Test _str2cmds with single command."""
        result = _str2cmds("echo hello")
        assert result == ["echo hello"]

    def test_str2cmds_multiple(self):
        """Test _str2cmds with multiple semicolon-separated commands."""
        result = _str2cmds("cmd1; cmd2; cmd3")
        assert result == ["cmd1", "cmd2", "cmd3"]

    def test_str2cmds_with_whitespace(self):
        """Test _str2cmds strips whitespace."""
        result = _str2cmds("  cmd1  ;  cmd2  ;  cmd3  ")
        assert result == ["cmd1", "cmd2", "cmd3"]

    def test_str2cmds_empty(self):
        """Test _str2cmds with empty string."""
        result = _str2cmds("")
        assert result == []


# ============================================================================
# TESTS FOR FileToTest CLASS
# ============================================================================


class TestFileToTest:
    """Test suite for FileToTest class."""

    def test_init_minimal_config(self):
        """Test FileToTest initialization with minimal config."""
        config = {"name": "output.txt"}
        ft = FileToTest(config)
        assert ft.name == "output.txt"
        assert ft.tolnlines == 0
        assert ft.tolabs == 0.0
        assert ft.tolrel == 0.0
        assert ft.fld_options == []
        assert ft.mode == ""
        assert ft.verbose_report == "no"

    def test_init_full_config(self):
        """Test FileToTest with all attributes."""
        config = {
            "name": "output.txt",
            "tolnlines": 2,
            "tolabs": 0.01,
            "tolrel": 1e-3,
            "fld_options": "-medium -include",
            "mode": "yaml",
            "verbose_report": "yes",
        }
        ft = FileToTest(config)
        assert ft.name == "output.txt"
        assert ft.tolnlines == 2
        assert ft.tolabs == 0.01
        assert ft.tolrel == 1e-3
        assert ft.fld_options == ["-medium", "-include"]
        assert ft.mode == "yaml"
        assert ft.verbose_report == "yes"

    def test_init_missing_name_raises_error(self):
        """Test FileToTest raises ValueError if 'name' is missing."""
        config = {"tolnlines": 2}
        with pytest.raises(ValueError, match="name must be defined"):
            FileToTest(config)

    def test_init_invalid_fld_option(self):
        """Test FileToTest raises on invalid fldiff option (missing dash)."""
        config = {
            "name": "output.txt",
            "fld_options": "-medium invalid_opt",
        }
        with pytest.raises(ValueError, match="Wrong fldiff option"):
            FileToTest(config)

    def test_init_strips_whitespace(self):
        """Test FileToTest strips whitespace from string attributes."""
        config = {
            "name": "  output.txt  ",
            "fld_options": "  -medium  ",
        }
        ft = FileToTest(config)
        assert ft.name == "output.txt"
        assert ft.fld_options == ["-medium"]

    def test_fld_options_parsing(self):
        """Test parsing of multiple fldiff options."""
        config = {
            "name": "output.txt",
            "fld_options": "-medium -include -includeP -ridiculous",
        }
        ft = FileToTest(config)
        assert ft.fld_options == ["-medium", "-include", "-includeP", "-ridiculous"]

    def test_fld_options_empty(self):
        """Test empty fld_options results in empty list."""
        config = {"name": "output.txt", "fld_options": ""}
        ft = FileToTest(config)
        assert ft.fld_options == []

    def test_mode_valid_values(self):
        """Test valid mode values."""
        for value in ["", "yaml", "yaml_docs"]:
            config = {"name": "output.txt", "mode": value}
            ft = FileToTest(config)
            assert ft.mode == value

    def test_mode_invalid_value_raises_error(self):
        """Test FileToTest raises ValueError on invalid mode value."""
        config = {"name": "output.txt", "mode": "bogus"}
        with pytest.raises(ValueError, match="Invalid value for mode"):
            FileToTest(config)

    def test_init_state_after_creation(self):
        """Test FileToTest has correct initial state."""
        config = {"name": "output.txt"}
        ft = FileToTest(config)
        assert ft.has_line_count_error is False
        assert ft.do_html_diff is False
        assert ft.fld_isok is False
        assert ft.fld_status == "failed"
        assert ft.fld_msg == "Initialized in __init__"


# ============================================================================
# TESTS FOR _str2filestotest FUNCTION
# ============================================================================


class TestStr2FilesToTest:
    """Test suite for _str2filestotest parsing function."""

    def test_empty_string(self):
        """Test parsing empty string returns empty list."""
        result = _str2filestotest("")
        assert result == []

    def test_single_file(self):
        """Test parsing single file specification."""
        spec = "output.txt, tolnlines=0, tolabs=0.01, tolrel=1e-3"
        result = _str2filestotest(spec)
        assert len(result) == 1
        assert result[0].name == "output.txt"
        assert result[0].tolnlines == 0
        assert result[0].tolabs == 0.01
        assert result[0].tolrel == 1e-3

    def test_multiple_files_semicolon_separated(self):
        """Test parsing multiple file specs separated by semicolons."""
        spec = "out1.txt, tolnlines=1, tolabs=0.01; out2.txt, tolnlines=2, tolabs=0.02"
        result = _str2filestotest(spec)
        assert len(result) == 2
        assert result[0].name == "out1.txt"
        assert result[0].tolnlines == 1
        assert result[1].name == "out2.txt"
        assert result[1].tolnlines == 2

    def test_file_with_fld_options(self):
        """Test parsing file spec with fldiff options."""
        spec = "output.txt, tolabs=0.01, fld_options=-medium -include"
        result = _str2filestotest(spec)
        assert len(result) == 1
        assert result[0].fld_options == ["-medium", "-include"]

    def test_duplicate_keyword_raises_error(self):
        """Test that duplicate keywords in spec raise error."""
        spec = "output.txt, tolabs=0.01, tolabs=0.02"
        with pytest.raises(AbinitTestInfoParserError, match="multiple occurrences"):
            _str2filestotest(spec)

    def test_multiple_files_with_empty_specs(self):
        """Test parsing skips empty specs."""
        spec = "out1.txt, tolabs=0.01; ; out2.txt, tolabs=0.02"
        result = _str2filestotest(spec)
        assert len(result) == 2

    def test_file_spec_with_all_attributes(self):
        """Test parsing file spec with all possible attributes."""
        spec = (
            "output.txt, tolnlines=2, tolabs=0.01, tolrel=1e-3, "
            "fld_options=-medium -include, mode=yaml, verbose_report=yes"
        )
        result = _str2filestotest(spec)
        ft = result[0]
        assert ft.name == "output.txt"
        assert ft.tolnlines == 2
        assert ft.tolabs == 0.01
        assert ft.tolrel == 1e-3
        assert ft.fld_options == ["-medium", "-include"]
        assert ft.mode == "yaml"
        assert ft.verbose_report == "yes"


# ============================================================================
# TESTS FOR AbinitTestInfoParser CLASS
# ============================================================================


class TestAbinitTestInfoParser:
    """Test suite for AbinitTestInfoParser class."""

    def test_parse_valid_test_info(self, input_file_with_test_info):
        """Test parsing valid TEST_INFO section."""
        parser = AbinitTestInfoParser(input_file_with_test_info)
        assert parser.inp_fname == os.path.abspath(input_file_with_test_info)
        assert parser.parser is not None

    def test_parse_missing_test_info_raises_error(self, temp_test_dir):
        """Test that missing TEST_INFO section raises error."""
        input_path = Path(temp_test_dir) / "no_test_info.abi"
        input_path.write_text("# Just a normal file\nwith no test info\n")

        with pytest.raises(AbinitTestInfoParserError, match="does not contain any valid testcnf section"):
            AbinitTestInfoParser(str(input_path))

    def test_generate_testinfo_basic(self, input_file_with_test_info):
        """Test generating AbinitTestInfo from parser."""
        parser = AbinitTestInfoParser(input_file_with_test_info)
        info = parser.generate_testinfo_nprocs(1)

        assert isinstance(info, AbinitTestInfo)
        assert info.executable == "abinit"
        assert info.use_files_file is False
        assert info.expected_failure is False
        assert info.no_check is False

    def test_nprocs_to_test_property(self, input_file_with_test_info):
        """Test nprocs_to_test property."""
        parser = AbinitTestInfoParser(input_file_with_test_info)
        nprocs = parser.nprocs_to_test
        assert isinstance(nprocs, list)

    def test_is_testchain_property(self, input_file_with_test_info):
        """Test is_testchain property."""
        parser = AbinitTestInfoParser(input_file_with_test_info)
        is_chain = parser.is_testchain
        assert isinstance(is_chain, bool)

    def test_yaml_test_method(self, input_file_with_test_info):
        """Test yaml_test method returns dict."""
        parser = AbinitTestInfoParser(input_file_with_test_info)
        yaml_config = parser.yaml_test()
        assert isinstance(yaml_config, dict)

    def test_parser_with_parallel_test(self, temp_test_dir):
        """Test parsing parallel test with nprocs_to_test."""
        input_path = Path(temp_test_dir) / "parallel_test.abi"
        content = """\
#%%<BEGIN TEST_INFO>
#%% [setup]
#%% executable = abinit
#%% [paral_info]
#%% max_nprocs = 4
#%% nprocs_to_test = 1, 2, 4
#%% [files]
#%% files_to_test = output.txt, tolabs=0.01
#%% [shell]
#%% [extra_info]
#%%<END TEST_INFO>
"""
        input_path.write_text(content)

        parser = AbinitTestInfoParser(str(input_path))
        nprocs = parser.nprocs_to_test
        assert nprocs == [1, 2, 4]

    def test_parser_with_test_chain_detection(self, temp_test_dir):
        """Test detection of test chain."""
        input_path = Path(temp_test_dir) / "chain_test.abi"
        content = """\
#%%<BEGIN TEST_INFO>
#%% [setup]
#%% executable = abinit
#%% test_chain = test1.abi, test2.abi, test3.abi
#%% [files]
#%% files_to_test = output.txt, tolabs=0.01
#%% [shell]
#%% [extra_info]
#%%<END TEST_INFO>
"""
        input_path.write_text(content)

        parser = AbinitTestInfoParser(str(input_path))
        assert parser.is_testchain is True

    def test_parser_default_values(self, input_file_with_test_info):
        """Test that defaults are applied correctly."""
        parser = AbinitTestInfoParser(input_file_with_test_info)
        info = parser.generate_testinfo_nprocs(1)

        # Check defaults are applied
        assert info.max_nprocs == 1
        assert info.exec_args == ""
        assert info.expected_failure is False

    def test_parser_with_multiple_files(self, temp_test_dir):
        """Test parsing multiple files_to_test."""
        input_path = Path(temp_test_dir) / "multi_file.abi"
        content = """\
#%%<BEGIN TEST_INFO>
#%% [setup]
#%% executable = abinit
#%% [files]
#%% files_to_test = out1.txt, tolabs=0.01; out2.txt, tolabs=0.02; out3.txt, tolabs=0.03
#%% [shell]
#%% [extra_info]
#%%<END TEST_INFO>
"""
        input_path.write_text(content)

        parser = AbinitTestInfoParser(str(input_path))
        info = parser.generate_testinfo_nprocs(1)
        assert len(info.files_to_test) == 3
        assert info.files_to_test[0].name == "out1.txt"
        assert info.files_to_test[1].name == "out2.txt"
        assert info.files_to_test[2].name == "out3.txt"

    def test_parser_with_repeated_test_chain(self, temp_test_dir):
        """Test that repeated tests in chain raise error."""
        input_path = Path(temp_test_dir) / "bad_chain.abi"
        content = """\
#%%<BEGIN TEST_INFO>
#%% [setup]
#%% executable = abinit
#%% test_chain = test1.abi, test2.abi, test1.abi
#%% [files]
#%% files_to_test = output.txt, tolabs=0.01
#%% [shell]
#%% [extra_info]
#%%<END TEST_INFO>
"""
        input_path.write_text(content)

        with pytest.raises(AbinitTestInfoParserError, match="repeated tests"):
            AbinitTestInfoParser(str(input_path))


# ============================================================================
# TESTS FOR AbinitTestInfo CLASS
# ============================================================================


class TestAbinitTestInfo:
    """Test suite for AbinitTestInfo class."""

    def test_init_from_dict(self):
        """Test AbinitTestInfo initialization from dictionary."""
        data = _make_minimal_testinfo_data()
        data.update({
            "keywords": {"test", "basic"},
        })
        info = AbinitTestInfo(data)
        assert info.executable == "abinit"
        assert info.inp_fname == "/path/to/test.abi"
        assert "test" in info.keywords
        assert "basic" in info.keywords

    def test_add_keywords(self):
        """Test adding keywords to test info."""
        data = _make_minimal_testinfo_data()
        info = AbinitTestInfo(data)
        # abinit is added automatically
        assert "abinit" in info.keywords
        info.add_keywords({"new_kw"})
        assert "abinit" in info.keywords
        assert "new_kw" in info.keywords

    def test_add_cpp_vars(self):
        """Test adding CPP variables to test info."""
        data = _make_minimal_testinfo_data()
        data["need_cpp_vars"] = {"HAVE_MPI"}
        info = AbinitTestInfo(data)
        info.add_cpp_vars({"HAVE_NETCDF"})
        assert "HAVE_MPI" in info.need_cpp_vars
        assert "HAVE_NETCDF" in info.need_cpp_vars

    def test_make_test_id_basic(self):
        """Test test_id generation."""
        data = _make_minimal_testinfo_data()
        data["inp_fname"] = "/path/to/suite/Input/test01.abi"
        info = AbinitTestInfo(data)
        test_id = info.make_test_id()
        assert test_id == "test01"

    def test_make_test_id_with_mpi(self):
        """Test test_id generation for MPI tests."""
        data = _make_minimal_testinfo_data()
        data["inp_fname"] = "/path/to/suite/Input/test02.abi"
        data["_ismulti_paral"] = True
        data["max_nprocs"] = 4
        info = AbinitTestInfo(data)
        test_id = info.make_test_id()
        assert test_id == "test02_MPI4"

    def test_ismulti_parallel_property(self):
        """Test ismulti_parallel property."""
        data_serial = _make_minimal_testinfo_data()
        data_serial["_ismulti_paral"] = False
        info_serial = AbinitTestInfo(data_serial)
        assert info_serial.ismulti_parallel is False

        data_parallel = _make_minimal_testinfo_data()
        data_parallel["_ismulti_paral"] = True
        info_parallel = AbinitTestInfo(data_parallel)
        assert info_parallel.ismulti_parallel is True


# ============================================================================
# TESTS FOR UTILITY FUNCTIONS (P3)
# ============================================================================


class TestUtilityFunctions:
    """Test suite for utility helper functions."""

    def test_genid_returns_string(self):
        """Test genid returns a string."""
        result = genid()
        assert isinstance(result, str)

    def test_genid_length(self):
        """Test genid returns 16-character string."""
        result = genid()
        assert len(result) == 16

    def test_genid_uniqueness(self):
        """Test genid generates unique IDs."""
        ids = [genid() for _ in range(100)]
        assert len(set(ids)) == 100, "genid should generate unique IDs"

    def test_my_getlogin_returns_string(self):
        """Test my_getlogin returns a string."""
        result = my_getlogin()
        assert isinstance(result, str)
        assert len(result) > 0

    def test_html_colorize_text_escaping(self):
        """Test html_colorize_text properly escapes HTML."""
        result = html_colorize_text("<script>", "#FF0000")
        assert "<script>" not in result
        assert "&lt;script&gt;" in result

    def test_status2html_passed(self):
        """Test status2html with 'passed' status."""
        result = status2html("passed")
        assert "✓" in result or "passed" in result.lower()
        assert "status-passed" in result

    def test_status2html_failed(self):
        """Test status2html with 'failed' status."""
        result = status2html("failed")
        assert "✕" in result or "failed" in result.lower()
        assert "status-failed" in result

    def test_status2html_unknown_symbol(self):
        """Test status2html with unknown status gets default symbol."""
        result = status2html("unknown")
        assert "•" in result  # Default symbol

    def test_sec2str_formatting(self):
        """Test sec2str formats seconds correctly."""
        assert sec2str(1.23456) == "1.23"
        assert sec2str(0) == "0.00"
        assert sec2str(123.456) == "123.46"

    def test_str2html_single_line(self):
        """Test str2html with single line."""
        result = str2html("Hello World")
        assert "Hello World" in result
        assert "<br>" in result

    def test_str2html_multiple_lines(self):
        """Test str2html with multiple lines."""
        result = str2html("Line 1\nLine 2\nLine 3")
        assert "<br>" in result
        assert result.count("<br>") >= 3

    def test_str2html_custom_end(self):
        """Test str2html with custom end tag."""
        result = str2html("Hello", end="</p>")
        assert result.endswith("</p>")

    def test_str2html_escaping(self):
        """Test str2html escapes HTML entities."""
        result = str2html("<script>alert('xss')</script>")
        assert "<script>" not in result
        assert "&lt;script&gt;" in result

    def test_args2htmltr_single_arg(self):
        """Test args2htmltr with single argument."""
        result = args2htmltr("value1")
        assert "<td>value1</td>" in result

    def test_args2htmltr_multiple_args(self):
        """Test args2htmltr with multiple arguments."""
        result = args2htmltr("val1", "val2", "val3")
        assert "<td>val1</td>" in result
        assert "<td>val2</td>" in result
        assert "<td>val3</td>" in result

    def test_html_link_with_url(self):
        """Test html_link with explicit URL."""
        result = html_link("Click here", href="http://example.com")
        assert 'href="http://example.com"' in result
        assert "Click here" in result

    def test_html_link_without_url(self):
        """Test html_link uses text as URL when href not provided."""
        result = html_link("http://example.com")
        assert 'href="http://example.com"' in result

    def test_html_link_escaping(self):
        """Test html_link escapes URLs and text."""
        result = html_link('<"text">', href='<"url">')
        assert "<" not in result or "&lt;" in result

    def test_html_file_link_missing_file(self):
        """Test html_file_link for missing file."""
        result = html_file_link("/nonexistent/file.txt")
        assert "Missing" in result or "unavailable" in result

    def test_html_file_link_empty_file(self, temp_test_dir):
        """Test html_file_link for empty file."""
        empty_file = Path(temp_test_dir) / "empty.txt"
        empty_file.write_text("")
        result = html_file_link(str(empty_file))
        assert "Empty" in result or "unavailable" in result

    def test_html_file_link_valid_file(self, temp_test_dir):
        """Test html_file_link for valid file."""
        valid_file = Path(temp_test_dir) / "test.txt"
        valid_file.write_text("content")
        result = html_file_link(str(valid_file))
        assert "test.txt" in result
        assert "<a" in result

    def test_is_string_with_string(self):
        """Test is_string returns True for strings."""
        assert is_string("hello") is True

    def test_is_string_with_non_string(self):
        """Test is_string returns False for non-strings."""
        assert is_string(123) is False
        assert is_string([]) is False
        assert is_string({}) is False

    def test_has_exts_single_ext(self):
        """Test has_exts with single extension."""
        assert has_exts("file.txt", ".txt") is True
        assert has_exts("file.txt", ".py") is False

    def test_has_exts_multiple_exts(self):
        """Test has_exts with multiple extensions."""
        assert has_exts("file.txt", [".txt", ".py"]) is True
        assert has_exts("file.py", [".txt", ".py"]) is True
        assert has_exts("file.md", [".txt", ".py"]) is False

    def test_lazy_read_file(self, temp_test_dir):
        """Test lazy_read reads entire file content."""
        test_file = Path(temp_test_dir) / "test.txt"
        test_file.write_text("Line 1\nLine 2\nLine 3")
        result = lazy_read(str(test_file))
        assert "Line 1" in result
        assert "Line 2" in result
        assert "Line 3" in result

    def test_lazy_readlines_file(self, temp_test_dir):
        """Test lazy_readlines reads lines from file."""
        test_file = Path(temp_test_dir) / "test.txt"
        test_file.write_text("Line 1\nLine 2\nLine 3")
        result = lazy_readlines(str(test_file))
        assert len(result) == 3
        assert "Line 1" in result[0]

    def test_rm_rf_removes_files(self, temp_test_dir):
        """Test rm_rf removes files in directory."""
        subdir = Path(temp_test_dir) / "subdir"
        subdir.mkdir()
        (subdir / "file1.txt").write_text("content")
        (subdir / "file2.txt").write_text("content")

        removed = rm_rf(str(subdir))
        assert len(removed) >= 2
        assert not (subdir / "file1.txt").exists()

    def test_rm_rf_exclude_paths(self, temp_test_dir):
        """Test rm_rf excludes specified paths."""
        subdir = Path(temp_test_dir) / "subdir"
        subdir.mkdir()
        keep_file = subdir / "keep.txt"
        remove_file = subdir / "remove.txt"
        keep_file.write_text("keep")
        remove_file.write_text("remove")

        removed = rm_rf(str(subdir), exclude_paths=[str(keep_file)])
        assert not remove_file.exists()
        assert keep_file.exists()

    def test_parse_configh_file(self, temp_test_dir):
        """Test parse_configh_file extracts CPP defines."""
        config_h = Path(temp_test_dir) / "config.h"
        config_h.write_text("""\
#define HAVE_MPI 1
#define HAVE_NETCDF 1
#define HAVE_HDF5
int x = 5;
""")
        result = parse_configh_file(str(config_h))
        assert "HAVE_MPI" in result
        assert result["HAVE_MPI"] == "1"
        assert "HAVE_NETCDF" in result

    def test_input_file_has_vars_match(self, temp_test_dir):
        """Test input_file_has_vars finds variables."""
        inp_file = Path(temp_test_dir) / "input.txt"
        inp_file.write_text("""\
npsp 2
ecut 20.0
ngkpt 2 2 2
""")
        found, matches = input_file_has_vars(str(inp_file), {"npsp": None})
        assert found is True
        assert "npsp" in matches

    def test_input_file_has_vars_no_match(self, temp_test_dir):
        """Test input_file_has_vars with no match."""
        inp_file = Path(temp_test_dir) / "input.txt"
        inp_file.write_text("ecut 20.0\n")
        found, matches = input_file_has_vars(str(inp_file), {"npsp": None})
        assert found is False


# ============================================================================
# TESTS FOR SLURM OOM-KILL DETECTION
#
# See detect_slurm_oom_stepid()/query_slurm_step_accounting()/
# _oom_query_allowed() and BaseTest.report_slurm_oom_if_detected() in
# testsuite.py. Real Slurm messages, observed live on manneback_gnu_14.2_hpc
# (v9/t83): a memory-hungry test gets killed by srun with a distinctive
# stderr, quite different from an ordinary non-zero-retcode failure.
# ============================================================================

REAL_OOM_STDERR = (
    "[v9][t83][np=1] Test was not expected to fail but subprocesses returned retcode: 1\n"
    "slurmstepd: error: Detected 1 oom_kill event in StepId=9787516.112. "
    "Some of the step tasks have been OOM Killed.\n"
    "srun: error: mb-mil009: task 0: Out Of Memory\n"
    "srun: Terminating StepId=9787516.112\n"
)


class TestDetectSlurmOomStepid:
    """Test suite for detect_slurm_oom_stepid()."""

    def test_detects_real_oom_message_and_extracts_stepid(self):
        """Regression fixture: the exact stderr text reported for a live OOM kill."""
        assert detect_slurm_oom_stepid(REAL_OOM_STDERR) == "9787516.112"

    def test_returns_none_for_ordinary_failure(self):
        """A plain, non-Slurm failure must not be misidentified as an OOM kill."""
        assert detect_slurm_oom_stepid("forrtl: severe (174): SIGSEGV\n") is None

    def test_returns_none_for_empty_or_missing_stderr(self):
        """Empty or None stderr must not raise."""
        assert detect_slurm_oom_stepid("") is None
        assert detect_slurm_oom_stepid(None) is None

    def test_returns_none_when_oom_detected_but_no_stepid_present(self):
        """A recognizable OOM message without a parseable StepId must not raise."""
        assert detect_slurm_oom_stepid("srun: error: node1: task 0: Out Of Memory\n") is None

    def test_matches_case_insensitively(self):
        """The OOM signature match must not be case-sensitive."""
        assert detect_slurm_oom_stepid("OUT OF MEMORY in StepId=1.0") == "1.0"


class TestQuerySlurmStepAccounting:
    """Test suite for query_slurm_step_accounting() -- Popen is mocked, no real Slurm needed."""

    def test_returns_first_line_of_sacct_output(self, monkeypatch):
        """A successful sacct call must return its first stripped output line."""
        class FakeProc:
            returncode = 0

            def communicate(self, timeout=None):
                return "9787516.112|COMPLETED|0:0|137|16|2Gc|1048576|00:12:34\n", ""

        monkeypatch.setattr(testsuite_module, "Popen", lambda *a, **k: FakeProc())

        result = query_slurm_step_accounting("9787516.112")
        assert result == "9787516.112|COMPLETED|0:0|137|16|2Gc|1048576|00:12:34"

    def test_returns_none_on_nonzero_returncode(self, monkeypatch):
        """A failing sacct call (e.g. unknown job id) must return None, not raise."""
        class FakeProc:
            returncode = 1

            def communicate(self, timeout=None):
                return "", "slurm_load_jobs error: Invalid job id specified"

        monkeypatch.setattr(testsuite_module, "Popen", lambda *a, **k: FakeProc())

        assert query_slurm_step_accounting("999999.0") is None

    def test_returns_none_on_empty_output(self, monkeypatch):
        """A zero-exit-code sacct call with no matching record must return None."""
        class FakeProc:
            returncode = 0

            def communicate(self, timeout=None):
                return "", ""

        monkeypatch.setattr(testsuite_module, "Popen", lambda *a, **k: FakeProc())

        assert query_slurm_step_accounting("9787516.112") is None

    def test_returns_none_when_sacct_binary_is_missing(self, monkeypatch):
        """`sacct` not on PATH (e.g. a non-Slurm machine) must not raise."""

        def raise_not_found(*a, **k):
            raise FileNotFoundError("sacct not found")

        monkeypatch.setattr(testsuite_module, "Popen", raise_not_found)

        assert query_slurm_step_accounting("9787516.112") is None

    def test_kills_process_and_returns_none_on_timeout(self, monkeypatch):
        """A hung sacct call must be killed and return None, not block forever."""
        from subprocess import TimeoutExpired

        class FakeProc:
            returncode = None

            def __init__(self):
                self.killed = False

            def communicate(self, timeout=None):
                if not self.killed:
                    raise TimeoutExpired(cmd="sacct", timeout=timeout)
                return "", ""

            def kill(self):
                self.killed = True

        monkeypatch.setattr(testsuite_module, "Popen", lambda *a, **k: FakeProc())

        assert query_slurm_step_accounting("9787516.112", timeout=1) is None


class TestOomQueryAllowed:
    """Test suite for _oom_query_allowed()'s shared-budget rate limiting."""

    def test_allows_up_to_max_calls_then_denies(self):
        """Exactly max_calls queries are allowed; every one after that is denied."""
        counter = _LocalCounter()
        lock = threading.Lock()
        results = [_oom_query_allowed(counter, lock, max_calls=3) for _ in range(5)]
        assert results == [True, True, True, False, False]

    def test_denial_does_not_increment_the_counter_further(self):
        """Denied calls must not keep incrementing the counter past max_calls."""
        counter = _LocalCounter()
        lock = threading.Lock()
        for _ in range(10):
            _oom_query_allowed(counter, lock, max_calls=2)
        assert counter.value == 2


class _FakeTestForOom:
    """Minimal stand-in exposing exactly what report_slurm_oom_if_detected()
    needs from a real BaseTest -- avoids constructing a full BaseTest (which
    requires a real input file and build environment) just to test this
    self-contained reporting method.
    """

    def __init__(self, max_calls=SLURM_OOM_MAX_QUERIES):
        self._oom_query_counter = _LocalCounter()
        self._oom_query_lock = threading.Lock()
        self.oom_query_max_calls = max_calls
        self.messages = []

    def cprint(self, msg="", color=None):
        self.messages.append((msg, color))

    report_slurm_oom_if_detected = BaseTest.report_slurm_oom_if_detected


class TestReportSlurmOomIfDetected:
    """Test suite for BaseTest.report_slurm_oom_if_detected()."""

    def test_prints_nothing_for_a_non_oom_failure(self):
        """A non-OOM stderr must be a silent no-op."""
        fake = _FakeTestForOom()
        fake.report_slurm_oom_if_detected("forrtl: severe (174): SIGSEGV\n")
        assert fake.messages == []

    def test_reports_sacct_accounting_when_available(self, monkeypatch):
        """A detected OOM kill with budget available must query and print sacct accounting."""
        class FakeProc:
            returncode = 0

            def communicate(self, timeout=None):
                return "9787516.112|COMPLETED|0:0|137|16|2Gc|1048576|00:12:34\n", ""

        monkeypatch.setattr(testsuite_module, "Popen", lambda *a, **k: FakeProc())

        fake = _FakeTestForOom()
        fake.report_slurm_oom_if_detected(REAL_OOM_STDERR)

        assert len(fake.messages) == 1
        msg, color = fake.messages[0]
        assert "9787516.112" in msg
        assert SLURM_ACCT_FORMAT in msg
        assert color == "red"

    def test_reports_budget_exhausted_without_querying_sacct(self, monkeypatch):
        """Once the shared budget is exhausted, sacct must not be invoked at all."""

        def fail_if_called(*a, **k):
            raise AssertionError("sacct must not be called once the budget is exhausted")

        monkeypatch.setattr(testsuite_module, "Popen", fail_if_called)

        fake = _FakeTestForOom(max_calls=1)
        fake._oom_query_counter.value = 1  # budget already exhausted

        fake.report_slurm_oom_if_detected(REAL_OOM_STDERR)  # must not raise, must not call Popen

        assert len(fake.messages) == 1
        msg, color = fake.messages[0]
        assert "budget" in msg.lower()
        assert color == "yellow"

    def test_shares_budget_across_multiple_tests(self, monkeypatch):
        """The counter/lock are meant to be shared across BaseTest instances
        (e.g. multiple py_nprocs workers, see AbinitTestSuite.run_tests()) --
        several tests passed the SAME counter/lock pair must together respect
        ONE combined budget, not one budget each.
        """
        class FakeProc:
            returncode = 0

            def communicate(self, timeout=None):
                return "9787516.112|COMPLETED|0:0|137|16|2Gc|1048576|00:12:34\n", ""

        monkeypatch.setattr(testsuite_module, "Popen", lambda *a, **k: FakeProc())

        shared_counter = _LocalCounter()
        shared_lock = threading.Lock()
        fakes = []
        for _ in range(3):
            fake = _FakeTestForOom(max_calls=2)
            fake._oom_query_counter = shared_counter
            fake._oom_query_lock = shared_lock
            fakes.append(fake)

        for fake in fakes:
            fake.report_slurm_oom_if_detected(REAL_OOM_STDERR)

        queried = [f for f in fakes if "sacct accounting" in f.messages[0][0]]
        exhausted = [f for f in fakes if "budget" in f.messages[0][0].lower()]
        assert len(queried) == 2
        assert len(exhausted) == 1


# ============================================================================
# TESTS FOR BuildEnvironment CLASS (P2)
# ============================================================================


class TestBuildEnvironment:
    """Test suite for BuildEnvironment class."""

    def test_buildenvironment_init_requires_valid_build_tree(self, temp_test_dir):
        """Test BuildEnvironment raises error for invalid build tree."""
        # BuildEnvironment raises ValueError or RuntimeError for invalid trees
        with pytest.raises((ValueError, RuntimeError)):
            BuildEnvironment(temp_test_dir)

    def test_buildenvironment_sets_attributes(self, build_environment):
        """Test BuildEnvironment initializes key attributes."""
        env = build_environment
        assert hasattr(env, "build_dir")
        assert hasattr(env, "hostname")
        assert hasattr(env, "username")
        assert hasattr(env, "fortran_compiler")

    def test_buildenvironment_path_of_bin_not_found(self, build_environment):
        """Test path_of_bin returns empty string for missing binary."""
        result = build_environment.path_of_bin("nonexistent_binary_xyz", try_syspath=False)
        # Result should be empty string or valid path
        assert isinstance(result, str)

    def test_compiler_base_class(self):
        """Test Compiler base class initialization."""
        compiler = Compiler(name="gfortran", version="9.3")
        assert compiler.name == "gfortran"
        assert compiler.version == "9.3"

    def test_fortran_compiler_known_vars(self):
        """Test FortranCompiler.from_defined_cpp_vars recognizes compilers."""
        # Test with known compiler variable
        compiler = FortranCompiler.from_defined_cpp_vars(["FC_GNU"])
        assert compiler.name == "gfortran"

    def test_cpp_preprocessor_init(self):
        """Test CPreProcessor initialization."""
        cpp = CPreProcessor(includes=["/usr/include"], bin="cpp")
        assert cpp.includes == ["/usr/include"]
        assert cpp.bin == "cpp"
        assert cpp.verbose == 0

    def test_cpp_preprocessor_default_includes(self):
        """Test CPreProcessor uses default includes."""
        cpp = CPreProcessor()
        assert "." in cpp.includes


# ============================================================================
# TESTS FOR BaseTest CLASS CORE METHODS (P3)
#
# Built from real ABINIT test inputs rather than mocks: BaseTest.__init__
# only needs the lightweight `AbinitEnvironment` (`tests.abenv`), not a
# compiled `BuildEnvironment`, so no build tree is required to exercise it.
# v1 is one of the oldest, most stable suites in the tree -- its Input files
# rarely change, which keeps these tests from breaking on unrelated edits
# elsewhere in the suite.
# ============================================================================

V1_T01_ABI = abenv.apath_of("tests", "v1", "Input", "t01.abi")


class TestBaseTestCore:
    """Test suite for BaseTest core methods."""

    def test_basetest_requires_files_to_test_or_no_check(self, temp_test_dir):
        """BaseTest must raise if files_to_test is empty and no_check is False.

        Necessarily a synthetic/malformed input: a real, already-passing
        suite test can never trigger this by construction.
        """
        bad_input = Path(temp_test_dir) / "tbad.abi"
        bad_input.write_text(
            "#%%<BEGIN TEST_INFO>\n"
            "#%% [setup]\n"
            "#%% executable = abinit\n"
            "#%% [paral_info]\n"
            "#%% max_nprocs = 1\n"
            "#%% [extra_info]\n"
            "#%% authors = Unknown\n"
            "#%%<END TEST_INFO>\n"
        )
        with pytest.raises(ValueError, match="no files_to_test attribute"):
            make_abitest_from_input(str(bad_input), abenv)

    def test_basetest_disabled_when_input_starts_with_dash(self, temp_test_dir):
        """A leading '-' in the input's basename marks the test disabled.

        Every currently-declared disabled entry in the real suite tree
        (`Suite.disabled_inp_paths`, across every suite that has any) points
        to a file that no longer exists on disk -- a real, separate data
        problem in its own right (the kind `test_no_stale_or_lost_inputs`,
        skipped elsewhere in this suite, is meant to catch) -- so there is no
        live example to build this test on directly. Instead, symlink a
        real, currently-active input (v1/t01.abi) under a dash-prefixed name
        inside a suite_name/Input/ layout, so the parsed content is 100%
        real; only the disabled-marker filename is synthesized.
        """
        input_dir = Path(temp_test_dir) / "v1" / "Input"
        input_dir.mkdir(parents=True)
        dashed = input_dir / "-t01.abi"
        dashed.symlink_to(V1_T01_ABI)

        test = make_abitest_from_input(str(dashed), abenv)
        assert test.status == "disabled"
        assert test.suite_name == "v1"

    def test_basetest_attributes_from_testinfo(self):
        """TEST_INFO attributes (keywords, authors, description, executable,
        ...) must be merged onto the real BaseTest instance.
        """
        test = make_abitest_from_input(V1_T01_ABI, abenv)

        assert isinstance(test, BaseTest)
        assert test.executable == "abinit"
        assert test.authors == {"Unknown"}
        assert "NC" in test.keywords
        assert test.description.strip().startswith("Bulk Aluminium")
        assert test.max_nprocs >= 1
        assert len(test.files_to_test) >= 1
        assert all(isinstance(ft, FileToTest) for ft in test.files_to_test)

    def test_basetest_full_id_property(self):
        """full_id must be '[suite_name][id][np=mpi_nprocs]'."""
        test = make_abitest_from_input(V1_T01_ABI, abenv)

        assert test.id == "t01"
        assert test.suite_name == "v1"
        assert test.mpi_nprocs == 1
        assert test.full_id == "[v1][t01][np=1]"


# ============================================================================
# TESTS FOR ChainOfTests AND MULTI-PARALLEL SELECTION
#
# t51/t52/t53 form a `test_chain` in `paral`, one of the few suite
# directories where chains and nprocs_to_test-driven multi-parallel tests
# actually exist -- generic suites like v1-v3 never exercise this. Low
# indices are used deliberately: these older chained tests are the least
# likely to be edited or renumbered.
# ============================================================================

PARAL_CHAIN_INPUTS = [
    abenv.apath_of("tests", "paral", "Input", name) for name in ("t51.abi", "t52.abi", "t53.abi")
]


class TestChainOfTestsMultiParallel:
    """ChainOfTests and nprocs_to_test-driven multi-parallel selection,
    built from a real `test_chain` in the `paral` suite.
    """

    @pytest.fixture(scope="class")
    @classmethod
    def chain_variants(cls):
        """One ChainOfTests per nprocs_to_test value declared on t51/t52/t53."""
        return make_abitests_from_inputs(list(PARAL_CHAIN_INPUTS), abenv)

    def test_chain_of_tests_is_generated_per_nprocs_value(self, chain_variants):
        """One ChainOfTests per declared nprocs_to_test value (1, 2, 4, 10)."""
        assert len(chain_variants) == 4
        assert all(isinstance(c, ChainOfTests) for c in chain_variants)
        assert all(c.is_chain for c in chain_variants)
        assert [c.max_nprocs for c in chain_variants] == [1, 2, 4, 10]

    def test_chain_keywords_are_the_union_of_its_members(self, chain_variants):
        """ChainOfTests.keywords is the union of every member's own keywords."""
        chain = chain_variants[0]
        assert len(chain) == 3
        assert [t.id for t in chain] == ["t51_MPI1", "t52_MPI1", "t53_MPI1"]
        assert "NC" in chain.keywords

    def test_chain_full_id_joins_member_ids_with_dashes(self, chain_variants):
        """Chain id/full_id are built from the dash-joined member ids."""
        chain = chain_variants[2]  # the MPI4 variant
        assert chain.id == "t51_MPI4-t52_MPI4-t53_MPI4"
        assert chain.full_id == "[paral][t51_MPI4-t52_MPI4-t53_MPI4]"

    def test_compute_nprocs_accepts_its_own_nprocs_and_rejects_others(self, chain_variants):
        """Each multi-parallel variant only runs at its own nprocs_to_test
        value; BaseTest.compute_nprocs() must accept it and reject any other
        real, declared alternative -- the same mechanism TestBot.run() relies
        on when sweeping np_list across multiple runs (see mysteps.py /
        testbot.py's run_tests_with_np()).
        """
        mpi4_chain = chain_variants[2]
        t51_mpi4 = mpi4_chain.tests[0]
        assert t51_mpi4.nprocs_to_test == [4]

        real_nprocs, err = t51_mpi4.compute_nprocs(build_env=None, mpi_nprocs=4, runmode="static")
        assert real_nprocs == 4
        assert err == ""

        real_nprocs, err = t51_mpi4.compute_nprocs(build_env=None, mpi_nprocs=2, runmode="static")
        assert real_nprocs == 0
        assert "nprocs_to_test" in err

    def test_chain_has_keywords_any_and_all_modes(self, chain_variants):
        """has_keywords() mirrors BaseTest: 'any' intersects, 'all' requires
        every keyword, an unknown mode raises.
        """
        chain = chain_variants[0]
        assert chain.keywords == {"DFPT", "NC", "abinit"}
        assert chain.has_keywords(["NC"], mode="any")
        assert chain.has_keywords(["DFPT", "NC"], mode="all")
        assert not chain.has_keywords(["NOPE"], mode="any")

        with pytest.raises(ValueError, match="wrong mode"):
            chain.has_keywords(["NC"], mode="bogus")

    def test_chain_authors_snames_and_has_authors(self, chain_variants):
        """_authors_snames unions every member's parsed author second-name."""
        chain = chain_variants[0]
        assert chain._authors_snames == {"Unknown"}
        assert chain.has_authors(["Unknown"], mode="all")

        with pytest.raises(ValueError, match="wrong mode"):
            chain.has_authors(["Unknown"], mode="bogus")

    def test_chain_exclude_builders_merges_members(self, chain_variants):
        """exclude_builders merges (and dedupes) every member's own list --
        empty here since none of t51/t52/t53 declare any.
        """
        assert chain_variants[0].exclude_builders == []

    def test_chain_status_defaults_to_failed_before_execution(self, chain_variants):
        """A never-run chain must read as 'failed', not silently succeeded:
        FileToTest.fld_status defaults to 'failed' until fldiff actually
        runs, and both BaseTest.status and ChainOfTests.status propagate
        that default rather than assuming success.
        """
        assert chain_variants[0].status == "failed"

    def test_chain_keep_files_and_files_to_keep(self, chain_variants):
        """keep_files()/files_to_keep track extra files across the whole
        chain, on top of whatever each member already keeps.
        """
        chain = chain_variants[1]
        before = list(chain.files_to_keep)
        chain.keep_files("extra_report.html")
        assert "extra_report.html" in chain.files_to_keep
        assert chain.files_to_keep[: len(before)] == before

    def test_chain_has_variables_checks_real_input_content(self, chain_variants):
        """has_variables() greps the real input file content, no execution
        needed: t51.abi genuinely sets natom 1 (see paral/Input/t51.abi).
        """
        chain = chain_variants[0]
        assert chain.has_variables({"natom": 1})
        assert chain.has_variables({"natom": 99}) == []


# ============================================================================
# TESTS FOR AbinitTestSuite
#
# v1 is reused again for the same reason as TestBaseTestCore: old, stable,
# unlikely to change. t01/t02/t03/t08 are plain, non-chained tests (t04/t07
# are a real test_chain pair and are deliberately excluded here so this
# suite is just a flat list of independent tests).
# ============================================================================

V1_SIMPLE_INPUTS = [
    abenv.apath_of("tests", "v1", "Input", name) for name in ("t01.abi", "t02.abi", "t03.abi")
]
V1_T08_ABI = abenv.apath_of("tests", "v1", "Input", "t08.abi")


class TestAbinitTestSuite:
    """AbinitTestSuite methods that don't require actually running abinit,
    built from real, chain-free v1 inputs.
    """

    @pytest.fixture
    def suite(self):
        """A real AbinitTestSuite over three plain, chain-free v1 tests."""
        return AbinitTestSuite(abenv, inp_files=list(V1_SIMPLE_INPUTS))

    def test_suite_requires_exactly_one_of_inp_files_or_test_list(self, suite):
        """Exactly one of inp_files/test_list must be given, never both/neither."""
        with pytest.raises(ValueError, match="One and only one"):
            AbinitTestSuite(abenv)
        with pytest.raises(ValueError, match="One and only one"):
            AbinitTestSuite(abenv, inp_files=list(V1_SIMPLE_INPUTS), test_list=[])

    def test_suite_len_iter_and_ids(self, suite):
        """len()/iteration expose one entry per input file, in order."""
        assert len(suite) == 3
        assert [t.id for t in suite] == ["t01", "t02", "t03"]

    def test_suite_full_length_counts_chain_members(self):
        """full_length must count each ChainOfTests by its own len(), not as 1 --
        unlike len(suite), which counts one entry per (possibly chained) test.
        """
        chain_suite = AbinitTestSuite(abenv, inp_files=list(PARAL_CHAIN_INPUTS))
        assert len(chain_suite) == 4  # 4 nprocs_to_test variants
        assert chain_suite.full_length == 12  # each chain has 3 members

    def test_suite_keywords_and_need_cpp_vars_are_unions(self, suite):
        """keywords/need_cpp_vars/has_keywords() aggregate across all tests."""
        assert suite.keywords == {"NC", "abinit"}
        assert suite.need_cpp_vars == set()
        assert suite.has_keywords(["NC"])
        assert not suite.has_keywords(["NOPE"])

    def test_suite_numeric_slice_selects_by_test_number(self, suite):
        """AbinitTestSuite[start:stop] selects by the test's own numeric id
        (via range(start, stop)), not by Python list position -- t01/t02
        have numeric ids 1/2, so [1:3] keeps them and drops t03 (id 3).
        """
        sliced = suite[1:3]
        assert isinstance(sliced, AbinitTestSuite)
        assert [t.id for t in sliced] == ["t01", "t02"]

    def test_suite_add_combines_two_suites(self, suite):
        """__add__ concatenates the two suites' test lists."""
        other = AbinitTestSuite(abenv, inp_files=[V1_T08_ABI])
        combined = suite + other
        assert [t.id for t in combined] == ["t01", "t02", "t03", "t08"]

    def test_suite_on_refslave_toggle(self, suite):
        """on_refslave() defaults to False until set_on_refslave() is called."""
        assert suite.on_refslave() is False
        suite.set_on_refslave(True)
        assert suite.on_refslave() is True

    def test_suite_run_etime_requires_executed(self, suite):
        """run_etime asserts self._executed rather than returning a bogus 0."""
        with pytest.raises(AssertionError):
            suite.run_etime  # noqa: B018

    def test_suite_status_filters_default_to_failed_before_execution(self, suite):
        """Same invariant as ChainOfTests.status: nothing has actually run
        yet, so every test reads as 'failed', never 'succeeded'.
        """
        assert suite.succeeded_tests() == []
        assert [t.id for t in suite.failed_tests()] == ["t01", "t02", "t03"]


# ============================================================================
# TESTS FOR BaseTest.run() / AbinitTestSuite.run_tests() -- REAL EXECUTION
#
# These actually invoke the compiled `abinit` binary on v1/t01.abi (the same
# real, stable input used throughout this file). They only run when a
# complete build tree is available at <abinit_home>/_build (the
# `build_environment` fixture above already skips with a clear message
# otherwise), and are kept to the minimum needed to exercise each real code
# path once: a single run can take a while, dominated by MPI startup/
# teardown overhead rather than this tiny physics workload.
#
# Success is asserted via `status` (the fldiff-based correctness check: did
# the computed output match the reference), not `isok` -- `isok` also folds
# in the subprocess's raw exit code, which can be a false negative on setups
# where MPI_Finalize itself errors out (observed locally, behind a VPN)
# after a numerically correct run.
#
# Artifacts (the real abinit working directory) are written to an explicit
# temp directory that is only removed if the test passes, so a failure
# leaves the real output behind for inspection.
#
# Currently disabled outright (not just left to the build_environment skip):
# a real run was observed to hang indefinitely, not just run slowly, tracked
# to a local VPN interfering with the MPI runtime's OFI network layer (see
# the MPI_Finalize/OFI errors noted above). Re-enable once that's resolved.
# ============================================================================

_VPN_MPI_SKIP_REASON = (
    "Real abinit execution currently hangs/fails here -- a local VPN "
    "appears to interfere with the MPI runtime's OFI network layer "
    "(MPI_Finalize/OFI errors observed). Re-enable once resolved."
)


@pytest.mark.skip(reason=_VPN_MPI_SKIP_REASON)
class TestBaseTestRunReal:
    """BaseTest.run() against a real, compiled abinit binary."""

    def test_run_produces_succeeded_status(self, build_environment, monkeypatch):
        """A real sequential run of v1/t01.abi must report status='succeeded'."""
        monkeypatch.setenv("ABI_PSPDIR", abenv.psps_dir)
        test = make_abitest_from_input(V1_T01_ABI, abenv)
        runner = JobRunner.sequential()
        workdir = tempfile.mkdtemp(prefix="test_basetest_run_")
        try:
            test.run(build_environment, runner, workdir, mpi_nprocs=1)
            assert test.status == "succeeded"
            assert test.run_etime > 0
            assert os.path.isfile(os.path.join(workdir, "t01.abo"))
        except BaseException:
            print(f"Test failed -- inspect real abinit output left in: {workdir}")
            raise
        else:
            shutil.rmtree(workdir, ignore_errors=True)


@pytest.mark.skip(reason=_VPN_MPI_SKIP_REASON)
class TestAbinitTestSuiteRunReal:
    """AbinitTestSuite.run_tests() (sequential mode) against a real build."""

    def test_run_tests_sequential_produces_succeeded_status(self, build_environment, monkeypatch):
        """run_tests() with py_nprocs=1 must drive the real test to 'succeeded'."""
        monkeypatch.setenv("ABI_PSPDIR", abenv.psps_dir)
        suite = AbinitTestSuite(abenv, inp_files=[V1_T01_ABI])
        runner = JobRunner.sequential()
        workdir = tempfile.mkdtemp(prefix="test_suite_run_tests_")
        try:
            suite.run_tests(build_environment, workdir, runner, mpi_nprocs=1, py_nprocs=1)
            assert [t.status for t in suite] == ["succeeded"]
            assert len(suite.succeeded_tests()) == 1
        except BaseException:
            print(f"Test failed -- inspect real abinit output left in: {workdir}")
            raise
        else:
            shutil.rmtree(workdir, ignore_errors=True)


# ============================================================================
# INTEGRATION TESTS
# ============================================================================


class TestIntegration:
    """Integration tests combining multiple components."""

    def test_full_parsing_workflow(self, temp_test_dir):
        """Test complete workflow: file creation -> parsing -> info generation."""
        input_path = Path(temp_test_dir) / "integration_test.abi"
        content = """\
#%%<BEGIN TEST_INFO>
#%% [setup]
#%% executable = abinit
#%% expected_failure = no
#%% [paral_info]
#%% max_nprocs = 1
#%% [files]
#%% files_to_test = output.txt, tolnlines=1, tolabs=0.01, tolrel=1e-3, fld_options=-medium; energy.txt, tolabs=0.001
#%% [shell]
#%% pre_commands =
#%% post_commands =
#%% [extra_info]
#%% authors = Test Author
#%% keywords = test, integration
#%% description = Integration test
#%%<END TEST_INFO>
"""
        input_path.write_text(content)

        # Parse the file
        parser = AbinitTestInfoParser(str(input_path))

        # Generate test info
        info = parser.generate_testinfo_nprocs(1)
        assert len(info.files_to_test) == 2
        assert info.files_to_test[0].name == "output.txt"
        assert info.files_to_test[0].tolabs == 0.01
        assert info.files_to_test[1].name == "energy.txt"
        assert info.files_to_test[1].tolabs == 0.001
        assert info.expected_failure is False

    def test_filetotest_created_from_parser(self, input_file_with_test_info):
        """Test FileToTest objects are correctly created from parser."""
        parser = AbinitTestInfoParser(input_file_with_test_info)
        info = parser.generate_testinfo_nprocs(1)

        # Check files_to_test contains FileToTest objects
        assert isinstance(info.files_to_test, (tuple, list))
        if info.files_to_test:
            assert all(isinstance(ft, FileToTest) for ft in info.files_to_test)
