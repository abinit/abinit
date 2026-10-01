"""
Unit tests for abo_file_analysis.py module.

Tests cover AboFileAnalysis.extract()'s dataset/iteration parsing (every line
pattern it recognizes) and AboFileAnalysis.compare_with()'s tolerance-based
comparison logic.
"""

from __future__ import annotations

import pytest

from .abo_file_analysis import AboDataset, AboFileAnalysis

ABO_CONTENT = """\
Some preamble line, before any dataset.
== DATASET  1 ==================================================
 meta: {optdriver: 0, }
 At Broyd/MD step   5, gradients are converged
 At SCF step   8, etot is converged
== DATASET  2 ==================================================
 meta: {optdriver: 3, }
 At SCF step   6,        nres2 [something]
is converged :  diff(etot_el-etot_pos)= 1.0e-12
 nstep=   12 was not enough SCF cycles to converge;
 ntime=   10 was not enough Broyd/MD steps to converge gradients
== END DATASET(S) ===============================================
"""


@pytest.fixture
def abo_file(tmp_path):
    def _make(content: str, name: str = "run.abo"):
        p = tmp_path / name
        p.write_text(content)
        return str(p)

    return _make


class TestExtract:
    """Test AboFileAnalysis.extract()'s line-pattern recognition."""

    def test_extract_creates_one_dataset_per_marker(self, abo_file):
        analysis = AboFileAnalysis(abo_file(ABO_CONTENT), option="iterations")
        assert len(analysis.dtsets) == 2
        assert [d.number for d in analysis.dtsets] == [1, 2]

    def test_optdriver_is_read(self, abo_file):
        analysis = AboFileAnalysis(abo_file(ABO_CONTENT), option="iterations")
        assert analysis.dtsets[0].optdriver == 0
        assert analysis.dtsets[1].optdriver == 3

    def test_scf_converged_on_same_line(self, abo_file):
        analysis = AboFileAnalysis(abo_file(ABO_CONTENT), option="iterations")
        assert analysis.dtsets[0].SCF_niter == [8]

    def test_md_converged_same_line(self, abo_file):
        analysis = AboFileAnalysis(abo_file(ABO_CONTENT), option="iterations")
        assert analysis.dtsets[0].MD_niter == 5

    def test_scf_converged_message_on_next_line(self, abo_file):
        analysis = AboFileAnalysis(abo_file(ABO_CONTENT), option="iterations")
        assert 6 in analysis.dtsets[1].SCF_niter

    def test_scf_not_enough_cycles(self, abo_file):
        analysis = AboFileAnalysis(abo_file(ABO_CONTENT), option="iterations")
        assert 12 in analysis.dtsets[1].SCF_niter

    def test_md_not_enough_steps(self, abo_file):
        analysis = AboFileAnalysis(abo_file(ABO_CONTENT), option="iterations")
        assert analysis.dtsets[1].MD_niter == 10

    def test_without_iterations_option_leaves_counters_at_defaults(self, abo_file):
        # option="" as the *extract* option (via the "iterations" in option
        # guard) -- not to be confused with the empty-string *constructor*
        # option, which skips extraction entirely (see next test).
        analysis = AboFileAnalysis(abo_file(ABO_CONTENT), option="dummy")
        assert analysis.dtsets[0].MD_niter is None
        assert analysis.dtsets[0].SCF_niter == []

    def test_empty_constructor_option_skips_extraction_entirely(self, abo_file):
        analysis = AboFileAnalysis(abo_file(ABO_CONTENT), option="")
        assert not hasattr(analysis, "dtsets")

    def test_single_dataset_without_explicit_end_marker(self, abo_file):
        # extract() only appends the in-progress dataset to dataset_list when
        # it hits "== END DATASET(S)" or the *next* "== DATASET" marker --
        # a file with exactly one dataset and no end marker is a real (if
        # unusual) edge case worth locking down.
        content = "== DATASET  1 ==========\n meta: {optdriver: 1, }\n"
        analysis = AboFileAnalysis(abo_file(content), option="iterations")
        assert analysis.dtsets == []


class TestCompareWith:
    """Test AboFileAnalysis.compare_with()'s tolerance-based comparison."""

    def _make_analysis(self, abo_file, md_niter: int, scf_niter: list[int], name: str) -> AboFileAnalysis:
        lines = ["== DATASET  1 ==========\n", " meta: {optdriver: 0, }\n"]
        lines.append(f" At Broyd/MD step   {md_niter}, gradients are converged\n")
        for n in scf_niter:
            lines.append(f" At SCF step   {n}, etot is converged\n")
            # A separator with no "converged" substring: extract() also
            # matches "At SCF step" lines whose *next* line says "converged"
            # (to catch a converged message wrapped onto a second line), so
            # two "same-line converged" SCF entries placed back to back would
            # otherwise double-count the first one via that second check.
            lines.append(" -- next iteration --\n")
        lines.append("== END DATASET(S) ======\n")
        return AboFileAnalysis(abo_file("".join(lines), name=name), option="iterations")

    def test_identical_files_succeed_with_zero_tolerance(self, abo_file):
        a = self._make_analysis(abo_file, 5, [8], "a.abo")
        b = self._make_analysis(abo_file, 5, [8], "b.abo")
        status, err_msg, err_msg_short = a.compare_with(b, option="iterations")
        assert status == "succeeded"
        assert err_msg == ""
        assert err_msg_short == ""

    def test_none_other_file_raises_valueerror(self, abo_file):
        a = self._make_analysis(abo_file, 5, [8], "a.abo")
        with pytest.raises(ValueError, match="no abo file provided"):
            a.compare_with(None, option="iterations")

    def test_mismatched_dataset_count_raises_valueerror(self, abo_file):
        a = self._make_analysis(abo_file, 5, [8], "a.abo")
        b = AboFileAnalysis(
            abo_file(
                "== DATASET  1 ==========\n"
                " meta: {optdriver: 0, }\n"
                "== DATASET  2 ==========\n"
                " meta: {optdriver: 0, }\n"
                "== END DATASET(S) ======\n",
                name="b.abo",
            ),
            option="iterations",
        )
        with pytest.raises(ValueError, match="different dataset numbers"):
            a.compare_with(b, option="iterations")

    def test_md_niter_within_small_tolerance_succeeds(self, abo_file):
        # MD_niter=5 (<=8, so tol_small applies): ceil(5*1.20) = 6.
        a = self._make_analysis(abo_file, 5, [8], "a.abo")
        b = self._make_analysis(abo_file, 6, [8], "b.abo")
        status, _, _ = a.compare_with(b, option="iterations", percent_allowed_small=20)
        assert status == "succeeded"

    def test_md_niter_beyond_small_tolerance_fails(self, abo_file):
        a = self._make_analysis(abo_file, 5, [8], "a.abo")
        b = self._make_analysis(abo_file, 7, [8], "b.abo")
        status, err_msg, err_msg_short = a.compare_with(b, option="iterations", percent_allowed_small=20)
        assert status == "failed"
        assert "MD/relax iterations differs" in err_msg
        assert "MD/relax cycle" in err_msg_short

    def test_md_niter_uses_large_tolerance_above_threshold(self, abo_file):
        # MD_niter=10 (>8, so tol_large applies, tol_small is ignored).
        a = self._make_analysis(abo_file, 10, [8], "a.abo")
        b = self._make_analysis(abo_file, 11, [8], "b.abo")
        # 10% of 10 rounds up to 1, so 11 is within tolerance under tol_large
        # but would fail under tol_small (which is 0 here, unused for this case).
        status, _, _ = a.compare_with(
            b, option="iterations", percent_allowed_small=0, percent_allowed_large=10
        )
        assert status == "succeeded"

    def test_scf_niter_single_cycle_beyond_tolerance_fails(self, abo_file):
        a = self._make_analysis(abo_file, 5, [8], "a.abo")
        b = self._make_analysis(abo_file, 5, [15], "b.abo")
        status, err_msg, err_msg_short = a.compare_with(b, option="iterations", percent_allowed_small=20)
        assert status == "failed"
        assert "[non-]SCF iterations differs" in err_msg
        assert "SCF_iter" in err_msg_short

    def test_scf_niter_multiple_cycles_reports_cycle_number(self, abo_file):
        a = self._make_analysis(abo_file, 5, [8, 8], "a.abo")
        b = self._make_analysis(abo_file, 5, [8, 20], "b.abo")
        status, err_msg, err_msg_short = a.compare_with(b, option="iterations", percent_allowed_small=20)
        assert status == "failed"
        assert "MD/relax cycle 2" in err_msg
        assert "MD/relax cycle 2, SCF_iter" in err_msg_short

    def test_comparison_skipped_when_iterations_not_requested(self, abo_file):
        a = self._make_analysis(abo_file, 5, [8], "a.abo")
        b = self._make_analysis(abo_file, 999, [999], "b.abo")
        # option requested here doesn't include "iterations".
        status, err_msg, _ = a.compare_with(b, option="dummy")
        assert status == "succeeded"
        assert err_msg == ""


class TestAboDataset:
    """Test the plain AboDataset data holder."""

    def test_defaults(self):
        d = AboDataset(3)
        assert d.number == 3
        assert d.MD_niter is None
        assert d.SCF_niter == []
