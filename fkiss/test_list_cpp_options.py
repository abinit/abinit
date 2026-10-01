"""
Unit tests for list_cpp_options.py module.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from .list_cpp_options import list_cpp_options

_SCRIPT_PATH = Path(__file__).with_name("list_cpp_options.py")


class TestListCppOptions:
    """Test list_cpp_options()'s CPP-directive scanning."""

    def test_counts_directives_across_fortran_files(self, tmp_path, capsys):
        (tmp_path / "a.F90").write_text(
            "#ifdef HAVE_MPI\n"
            "  call foo()\n"
            "#endif\n"
            "#ifdef HAVE_MPI\n"
            "  call bar()\n"
            "#endif\n"
        )
        (tmp_path / "b.f90").write_text("#ifdef HAVE_FFTW3\n  call baz()\n#endif\n")

        retcode = list_cpp_options(str(tmp_path))

        assert retcode == 0
        out = capsys.readouterr().out
        assert "HAVE_MPI" in out
        assert "HAVE_FFTW3" in out

    def test_non_fortran_files_are_ignored(self, tmp_path, capsys):
        (tmp_path / "notes.txt").write_text("#ifdef SHOULD_NOT_APPEAR\n")

        retcode = list_cpp_options(str(tmp_path))

        assert retcode == 0
        assert "SHOULD_NOT_APPEAR" not in capsys.readouterr().out

    def test_non_cpp_lines_are_ignored(self, tmp_path, capsys):
        (tmp_path / "a.F90").write_text("  ! not a cpp line: ifdef HAVE_MPI\n")

        retcode = list_cpp_options(str(tmp_path))

        assert retcode == 0
        assert "HAVE_MPI" not in capsys.readouterr().out

    def test_recurses_into_subdirectories(self, tmp_path, capsys):
        sub = tmp_path / "sub"
        sub.mkdir()
        (sub / "nested.F90").write_text("#ifdef HAVE_NETCDF\n")

        list_cpp_options(str(tmp_path))

        assert "HAVE_NETCDF" in capsys.readouterr().out

    def test_empty_tree_prints_header_only(self, tmp_path, capsys):
        retcode = list_cpp_options(str(tmp_path))
        assert retcode == 0
        out = capsys.readouterr().out
        assert "Option" in out
        assert "Occurrences" in out


class TestMainScript:
    """Test the __main__ CLI entry point runs as a real script.

    Regression test: `sys.exit(main(top))` referenced an undefined `main`
    (only `list_cpp_options` is defined here), so running this file
    directly (`python list_cpp_options.py <dir>`) always raised NameError.
    """

    def test_running_as_a_script_works(self, tmp_path):
        (tmp_path / "a.F90").write_text("#ifdef HAVE_MPI\n  call foo()\n#endif\n")

        result = subprocess.run(
            [sys.executable, str(_SCRIPT_PATH), str(tmp_path)],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert "HAVE_MPI" in result.stdout
