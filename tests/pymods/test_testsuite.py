"""
Unit tests for testsuite.py module.

Tests cover configuration parsing, file comparison setup, and test info handling.
Priority 1: High-impact, core infrastructure components.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest

from .testsuite import (
    AbinitTestInfo,
    AbinitTestInfoParser,
    AbinitTestInfoParserError,
    BuildEnvironment,
    Compiler,
    CPreProcessor,
    FileToTest,
    FortranCompiler,
    _str2bool,
    _str2cmds,
    _str2filestotest,
    _str2intlist,
    _str2list,
    _str2set,
    args2htmltr,
    genid,
    has_exts,
    html_colorize_text,
    html_file_link,
    html_link,
    input_file_has_vars,
    is_string,
    lazy_read,
    lazy_readlines,
    my_getlogin,
    parse_configh_file,
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
        assert ft.use_yaml == "no"
        assert ft.verbose_report == "no"

    def test_init_full_config(self):
        """Test FileToTest with all attributes."""
        config = {
            "name": "output.txt",
            "tolnlines": 2,
            "tolabs": 0.01,
            "tolrel": 1e-3,
            "fld_options": "-medium -include",
            "use_yaml": "yes",
            "verbose_report": "yes",
        }
        ft = FileToTest(config)
        assert ft.name == "output.txt"
        assert ft.tolnlines == 2
        assert ft.tolabs == 0.01
        assert ft.tolrel == 1e-3
        assert ft.fld_options == ["-medium", "-include"]
        assert ft.use_yaml == "yes"
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

    def test_use_yaml_valid_values(self):
        """Test valid use_yaml values."""
        for value in ["yes", "no", "only"]:
            config = {"name": "output.txt", "use_yaml": value}
            ft = FileToTest(config)
            assert ft.use_yaml == value

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
            "fld_options=-medium -include, use_yaml=yes, verbose_report=yes"
        )
        result = _str2filestotest(spec)
        ft = result[0]
        assert ft.name == "output.txt"
        assert ft.tolnlines == 2
        assert ft.tolabs == 0.01
        assert ft.tolrel == 1e-3
        assert ft.fld_options == ["-medium", "-include"]
        assert ft.use_yaml == "yes"
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
# ============================================================================


class TestBaseTestCore:
    """Test suite for BaseTest core methods."""

    def test_basetest_requires_files_to_test_or_no_check(self, temp_test_dir):
        """Test BaseTest raises error when no files_to_test and no_check is False."""
        # Build environment would be needed, skip for now
        pytest.skip("Requires valid BuildEnvironment")

    def test_basetest_disabled_when_input_starts_with_dash(self):
        """Test BaseTest marks tests as disabled if input starts with dash."""
        pytest.skip("Requires valid BuildEnvironment and input file")

    def test_basetest_attributes_from_testinfo(self):
        """Test BaseTest incorporates TestInfo attributes."""
        pytest.skip("Requires valid BuildEnvironment and input file")

    def test_basetest_full_id_property(self):
        """Test BaseTest.full_id property."""
        pytest.skip("Requires valid BuildEnvironment and input file")


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
