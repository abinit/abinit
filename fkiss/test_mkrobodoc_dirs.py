"""
Unit tests for mkrobodoc_dirs.py module.
"""

from __future__ import annotations

import sys

import pytest

from .mkrobodoc_dirs import WildCard, is_string, list_strings, main, mkrobodoc_files, robodoc_dheader


class TestIsString:
    def test_string_is_a_string(self):
        assert is_string("hello") is True

    @pytest.mark.parametrize("value", [1, 3.14, [], {}, None])
    def test_non_strings_are_not_strings(self, value):
        assert is_string(value) is False


class TestListStrings:
    def test_single_string_becomes_a_list(self):
        assert list_strings("A single string") == ["A single string"]

    def test_list_of_strings_is_unchanged(self):
        assert list_strings(["a", "b", "c"]) == ["a", "b", "c"]


class TestWildCard:
    def test_default_pattern_matches_everything(self):
        w = WildCard("")
        assert w.match("anything.txt") is True

    def test_filter_selects_matching_names(self):
        w = WildCard("*.nc|*.pdf")
        assert w.filter(["foo.nc", "bar.pdf", "hello.txt"]) == ["foo.nc", "bar.pdf"]

    def test_filter_accepts_a_single_string(self):
        w = WildCard("*.nc")
        assert w.filter("foo.nc") == ["foo.nc"]

    def test_match_true_and_false(self):
        w = WildCard("*.F90")
        assert w.match("m_crystal.F90") is True
        assert w.match("m_crystal.py") is False

    def test_str(self):
        w = WildCard("*.F90|*.finc")
        assert "WildCard" in str(w)
        assert w.pats == ["*.F90", "*.finc"]


class TestRobodocDheader:
    def test_header_contains_basename_and_sections(self):
        header = robodoc_dheader("/some/path/70_gw")
        assert "70_gw" in header
        assert "!!****d* ABINIT/70_gw" in header
        assert "!! CHILDREN" in header


class TestMkrobodocFiles:
    """Test mkrobodoc_files()'s directory walk and robodoc-file rewriting.

    os.walk() visits the given top directory itself too, and mkrobodoc_files()
    checks *every* visited directory for a robodoc file -- including top.
    Real usage always walks from "src" itself, which is why "src" is in
    EXCLUDE_BASEDIRS (it's never expected to have its own robodoc file); so
    every fixture below nests builder directories under a "src" top dir,
    matching real usage, instead of using tmp_path (an arbitrary pytest-
    generated name) directly as the walk root.
    """

    @pytest.fixture
    def src(self, tmp_path):
        d = tmp_path / "src"
        d.mkdir()
        return d

    def test_existing_robodoc_dir_gets_source_files_listed(self, src):
        d = src / "70_gw"
        d.mkdir()
        (d / "_70_gw_").write_text(
            "!!****d* ABINIT/70_gw\n!! NAME\n!! 70_gw\n!! CHILDREN\n"
        )
        (d / "m_b.F90").write_text("module m_b\nend module m_b\n")
        (d / "m_a.F90").write_text("module m_a\nend module m_a\n")
        (d / "notes.txt").write_text("not a source file\n")

        n_wrong = mkrobodoc_files(str(src))

        content = (d / "_70_gw_").read_text()
        # Sorted alphabetically, only F90/finc files listed, txt excluded.
        assert content.index("m_a.F90") < content.index("m_b.F90")
        assert "notes.txt" not in content
        assert content.rstrip().endswith("!!***")
        assert n_wrong == 0

    def test_dir_without_source_files_is_left_untouched(self, src):
        d = src / "70_gw"
        d.mkdir()
        original = "!!****d* ABINIT/70_gw\n!! CHILDREN\n"
        (d / "_70_gw_").write_text(original)

        mkrobodoc_files(str(src))

        assert (d / "_70_gw_").read_text() == original

    def test_missing_children_marker_raises(self, src):
        d = src / "70_gw"
        d.mkdir()
        (d / "_70_gw_").write_text("!!****d* ABINIT/70_gw\n!! NAME\n")
        (d / "m_a.F90").write_text("module m_a\nend module m_a\n")

        with pytest.raises(ValueError, match="CHILDREN"):
            mkrobodoc_files(str(src))

    def test_dir_without_robodoc_file_is_reported_as_wrong(self, src):
        d = src / "not_an_abinit_dir"
        d.mkdir()
        (d / "m_a.F90").write_text("module m_a\nend module m_a\n")

        n_wrong = mkrobodoc_files(str(src))

        assert n_wrong == 1
        # A fresh robodoc file must have been created for it.
        created = d / "_not_an_abinit_dir_"
        assert created.is_file()
        assert "!! CHILDREN" in created.read_text()

    def test_excluded_basedirs_are_not_reported_as_wrong(self, src):
        # "src" (the walk root here) is itself in EXCLUDE_BASEDIRS.
        n_wrong = mkrobodoc_files(str(src))

        assert n_wrong == 0
        assert not (src / "_src_").exists()


class TestMain:
    """Test the main() CLI entry point."""

    def test_main_delegates_to_mkrobodoc_files(self, monkeypatch, tmp_path):
        src = tmp_path / "src"
        src.mkdir()
        monkeypatch.setattr(sys, "argv", ["mkrobodoc_dirs.py", str(src)])
        assert main() == 0

    def test_main_without_argument_raises(self, monkeypatch):
        monkeypatch.setattr(sys, "argv", ["mkrobodoc_dirs.py"])
        with pytest.raises(ValueError, match="Top level directory must be specified"):
            main()
