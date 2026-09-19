"""
Unit tests for diff.py module.

Tests cover the Abinit-specific junk heuristics (abinit_line_junk,
abinit_char_junk) and the main() CLI entry point across its four output
formats (context, unified, ndiff, html) and its -j/-f flags.
"""

from __future__ import annotations

import sys

import pytest

from .diff import abinit_char_junk, abinit_line_junk, main


class TestAbinitLineJunk:
    """Test the line-level junk heuristic used to synchronize diffs."""

    def test_minus_prefixed_line_is_junk(self):
        assert abinit_line_junk("- ignored line\n") is True

    def test_plus_prefixed_line_is_junk(self):
        assert abinit_line_junk("+ ignored line\n") is True

    def test_whitespace_only_line_is_junk(self):
        assert abinit_line_junk("   \n") is True

    def test_regular_line_is_not_junk(self):
        assert abinit_line_junk(" Here are some numbers : 0.4546\n") is False

    def test_empty_string_is_not_junk(self):
        # "".isspace() is False (unlike a whitespace-only string).
        assert abinit_line_junk("") is False


class TestAbinitCharJunk:
    """Test the character-level junk heuristic used by ndiff's charjunk."""

    @pytest.mark.parametrize("c", [" ", "\t", "\n", "0", "5", "9"])
    def test_whitespace_and_digits_are_junk(self, c):
        assert abinit_char_junk(c) is True

    @pytest.mark.parametrize("c", ["a", "Z", ".", "-", "="])
    def test_letters_and_punctuation_are_not_junk(self, c):
        assert abinit_char_junk(c) is False


@pytest.fixture
def file_pair(tmp_path):
    """Create a reference/compared file pair with one differing line."""
    fromfile = tmp_path / "ref.out"
    tofile = tmp_path / "new.out"
    fromfile.write_text("line one\nline two\nline three\n")
    tofile.write_text("line one\nline TWO\nline three\n")
    return fromfile, tofile


class TestMain:
    """Test the main() CLI entry point across its output formats."""

    def test_context_diff_is_default(self, monkeypatch, capsys, file_pair):
        fromfile, tofile = file_pair
        monkeypatch.setattr(sys, "argv", ["diff.py", str(fromfile), str(tofile)])
        with pytest.raises(SystemExit) as exc_info:
            main()
        assert exc_info.value.code == 0
        out = capsys.readouterr().out
        assert "***" in out  # context diff marker
        assert "line TWO" in out

    def test_unified_diff(self, monkeypatch, capsys, file_pair):
        fromfile, tofile = file_pair
        monkeypatch.setattr(sys, "argv", ["diff.py", "-u", str(fromfile), str(tofile)])
        with pytest.raises(SystemExit):
            main()
        out = capsys.readouterr().out
        assert out.startswith("---")
        assert "+line TWO" in out

    def test_ndiff(self, monkeypatch, capsys, file_pair):
        fromfile, tofile = file_pair
        monkeypatch.setattr(sys, "argv", ["diff.py", "-n", str(fromfile), str(tofile)])
        with pytest.raises(SystemExit):
            main()
        out = capsys.readouterr().out
        assert "- line two\n" in out
        assert "+ line TWO\n" in out

    def test_ndiff_with_abinit_junk_flag(self, monkeypatch, capsys, tmp_path):
        # abinit_char_junk() treats digits/whitespace as junk, so ndiff should
        # not flag a purely-numeric difference as an inline character change.
        fromfile = tmp_path / "ref.out"
        tofile = tmp_path / "new.out"
        fromfile.write_text("value 0.4546\n")
        tofile.write_text("value 0.4547\n")
        monkeypatch.setattr(
            sys, "argv", ["diff.py", "-n", "-j", str(fromfile), str(tofile)]
        )
        with pytest.raises(SystemExit):
            main()
        out = capsys.readouterr().out
        assert "value 0.4546" in out
        assert "value 0.4547" in out

    def test_html_side_by_side_file(self, monkeypatch, capsys, file_pair):
        fromfile, tofile = file_pair
        monkeypatch.setattr(sys, "argv", ["diff.py", "-m", str(fromfile), str(tofile)])
        with pytest.raises(SystemExit):
            main()
        out = capsys.readouterr().out
        assert "<html>" in out.lower()
        # HtmlDiff renders spaces as &nbsp;, so "line TWO" never appears literally.
        assert "line&nbsp;TWO" in out

    def test_html_table(self, monkeypatch, capsys, file_pair):
        fromfile, tofile = file_pair
        monkeypatch.setattr(sys, "argv", ["diff.py", "-t", str(fromfile), str(tofile)])
        with pytest.raises(SystemExit):
            main()
        out = capsys.readouterr().out
        assert "<table" in out.lower()

    def test_output_written_to_file(self, monkeypatch, capsys, tmp_path, file_pair):
        fromfile, tofile = file_pair
        outfile = tmp_path / "diff.txt"
        monkeypatch.setattr(
            sys, "argv", ["diff.py", "-u", "-f", str(outfile), str(fromfile), str(tofile)]
        )
        with pytest.raises(SystemExit):
            main()
        # Nothing goes to stdout when -f/--file is given.
        assert capsys.readouterr().out == ""
        content = outfile.read_text()
        assert "+line TWO" in content

    def test_context_lines_option_is_honored(self, monkeypatch, capsys, tmp_path):
        fromfile = tmp_path / "ref.out"
        tofile = tmp_path / "new.out"
        fromfile.write_text("".join(f"line {i}\n" for i in range(10)))
        tofile.write_text("".join(f"line {i}\n" for i in range(10)).replace("line 5", "line FIVE"))
        monkeypatch.setattr(
            sys, "argv", ["diff.py", "-u", "-l", "1", str(fromfile), str(tofile)]
        )
        with pytest.raises(SystemExit):
            main()
        out = capsys.readouterr().out
        # With only 1 line of context, lines far from the change (e.g. "line 0")
        # must not appear, while immediate neighbours ("line 4"/"line 6") do.
        assert "line 0\n" not in out
        assert " line 4\n" in out
        assert " line 6\n" in out
