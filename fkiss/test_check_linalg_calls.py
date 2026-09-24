"""
Unit tests for check_linalg_calls.py module.
"""

from __future__ import annotations

from .check_linalg_calls import main


class TestMain:
    """Test main()'s BLAS/LAPACK/BLACS/ScaLAPACK call-counting report."""

    def test_counts_calls_by_library(self, tmp_path, capsys):
        (tmp_path / "a.F90").write_text(
            "  call caxpy(n, a, x, incx, y, incy)\n"
            "  call caxpy(n, a, x, incx, y, incy)\n"
            "  call cbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu)\n"
            "  call cbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu)\n"
            "  call blacs_abort(icontxt, errcode)\n"
            "  call blacs_abort(icontxt, errcode)\n"
            "  call pslaed3(a, b)\n"
            "  call pslaed3(a, b)\n"
        )

        main(str(tmp_path))

        out = capsys.readouterr().out
        assert "BLAS use:" in out
        assert "caxpy" in out and "2" in out
        assert "LAPACK use:" in out
        assert "cbdsqr" in out
        assert "BLACS use:" in out
        assert "blacs_abort" in out
        assert "ScaLAPACK use:" in out
        assert "pslaed3" in out

    def test_call_not_at_start_of_line_is_ignored(self, tmp_path, capsys):
        (tmp_path / "a.F90").write_text("  if (ierr /= 0) call caxpy(n, a, x, incx, y, incy)\n")

        main(str(tmp_path))

        out = capsys.readouterr().out
        # "call" must be the first token on the line to be recognized.
        assert "caxpy" not in out

    def test_unknown_subroutine_is_not_counted(self, tmp_path, capsys):
        (tmp_path / "a.F90").write_text("  call my_own_routine(a, b)\n")

        main(str(tmp_path))

        out = capsys.readouterr().out
        assert "my_own_routine" not in out

    def test_blank_and_short_lines_do_not_crash(self, tmp_path, capsys):
        (tmp_path / "a.F90").write_text("\n   \ncall\n  call caxpy(n)\n")

        main(str(tmp_path))

        assert "caxpy" in capsys.readouterr().out

    def test_non_fortran_files_are_ignored(self, tmp_path, capsys):
        (tmp_path / "notes.txt").write_text("call caxpy(n)\n")

        main(str(tmp_path))

        assert "caxpy" not in capsys.readouterr().out

    def test_empty_tree_reports_zero_percent_use(self, tmp_path, capsys):
        main(str(tmp_path))

        out = capsys.readouterr().out
        assert "BLAS use: 0.0%" in out
        assert "LAPACK use: 0.0%" in out
        assert "BLACS use: 0.0%" in out
        assert "ScaLAPACK use: 0.0%" in out
