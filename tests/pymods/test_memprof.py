"""
Unit tests for memprof.py module.

Tests cover Entry parsing/derived properties, AbimemFile's log-file parsing
and analysis methods (small/large allocs, zero-sized, weird pointers, peaks,
hotspots, memleak detection), Heap/Stack containers, and light smoke coverage
of the matplotlib-based plot_* methods.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import pytest

from .memprof import AbimemFile, Entry, Heap, Stack, entries_to_dataframe


@pytest.fixture(autouse=True)
def _close_figures():
    """Close all matplotlib figures after each test.

    Without this, figures accumulate across tests in this file (plot_*()
    never closes what it creates) and stale cross-figure GridSpec state can
    make an unrelated *later* test's plt.figure()/add_subplot() call raise
    inside numpy/matplotlib internals -- each plot_* test passes in isolation
    but not after enough other figures have piled up in the same process.
    """
    yield
    import matplotlib.pyplot as plt

    plt.close("all")


def make_line(vname: str, action: str, ptr: int, size: int, file: str, line: int, tot_memory: int) -> str:
    """Build one raw abimem-log line matching Entry.from_line()'s format:
    a fixed 59-char-wide variable-name field, followed by whitespace-
    separated action/ptr/size/file/line/tot_memory tokens.
    """
    prefix = vname.ljust(59)
    return f"{prefix}{action} {ptr} {size} {file} {line} {tot_memory}\n"


class TestEntry:
    """Test the Entry namedtuple's parsing and derived (lazy) properties."""

    def test_from_line_parses_fields(self):
        line = make_line("myvar", "A", 4096, 800, "foo.F90", 10, 800)
        entry = Entry.from_line(line)
        assert entry.vname == "myvar"
        assert entry.action == "A"
        assert entry.ptr == 4096
        assert entry.size == 800
        assert entry.file == "foo.F90"
        assert entry.line == 10
        assert entry.tot_memory == 800

    def test_from_line_strips_internal_spaces_from_vname(self):
        # vname occupies line[:59]; internal spaces (padding artifacts) must
        # be removed, not just leading/trailing ones.
        line = make_line("my var", "A", 1, 8, "f.F90", 1, 8)
        entry = Entry.from_line(line)
        assert entry.vname == "myvar"

    def test_size_mb_positive_for_alloc(self):
        entry = Entry.from_line(make_line("v", "A", 1, 8 * 1024 ** 2, "f.F90", 1, 8 * 1024 ** 2))
        assert entry.size_mb == pytest.approx(1.0)

    def test_size_mb_negative_for_free(self):
        entry = Entry.from_line(make_line("v", "D", 1, 8 * 1024 ** 2, "f.F90", 1, 0))
        assert entry.size_mb == pytest.approx(-1.0)

    def test_tot_memory_mb(self):
        entry = Entry.from_line(make_line("v", "A", 1, 8, "f.F90", 1, 16 * 1024 ** 2))
        assert entry.tot_memory_mb == pytest.approx(2.0)

    def test_isalloc_isfree(self):
        alloc = Entry.from_line(make_line("v", "A", 1, 8, "f.F90", 1, 8))
        free = Entry.from_line(make_line("v", "D", 1, 8, "f.F90", 1, 0))
        assert alloc.isalloc is True
        assert alloc.isfree is False
        assert free.isalloc is False
        assert free.isfree is True

    def test_iszerosized(self):
        zero = Entry.from_line(make_line("v", "A", 1, 0, "f.F90", 1, 0))
        nonzero = Entry.from_line(make_line("v", "A", 1, 8, "f.F90", 1, 8))
        assert zero.iszerosized is True
        assert nonzero.iszerosized is False

    def test_locus_and_hash_and_eq(self):
        e1 = Entry.from_line(make_line("v", "A", 1, 8, "f.F90", 10, 8))
        e2 = Entry.from_line(make_line("v", "A", 2, 8, "f.F90", 10, 16))
        # locus does not include ptr/tot_memory, so these compare equal.
        assert e1.locus == e2.locus
        assert e1 == e2
        assert hash(e1) == hash(e2)
        assert not e1.__neq__(e2)

    def test_repr_and_to_repr(self):
        entry = Entry.from_line(make_line("v", "A", 255, 8, "f.F90", 1, 8))
        assert "0xff" in repr(entry)
        assert "addr=" in repr(entry)
        assert "addr=" not in entry.to_repr(with_addr=False)

    def test_frees_onheap_never_matches_a_free_against_its_own_allocation(self):
        # BUG (documented, not fixed here; also unused anywhere in the
        # module -- find_memleaks() reimplements the same check inline
        # instead of calling this method): the guard reads
        # `if (not self.isfree) or other.isalloc: return False`, which
        # returns False whenever `other` *is* an allocation -- i.e. for
        # every real "does this free deallocate that alloc" call. The `or`
        # is almost certainly missing a `not` on `other.isalloc`.
        alloc = Entry.from_line(make_line("v", "A", 1, 800, "f.F90", 1, 800))
        free = Entry.from_line(make_line("v", "D", 1, -800, "f.F90", 1, 0))
        assert free.frees_onheap(alloc) is False
        assert alloc.frees_onheap(free) is False  # alloc is not a free

    def test_frees_onheap_true_only_between_two_cancelling_frees(self):
        # The only inputs that satisfy the (buggy) guard above: self is a
        # free *and* other is also not an allocation.
        free1 = Entry.from_line(make_line("v", "D", 1, -800, "f.F90", 1, 0))
        free2 = Entry.from_line(make_line("v", "D", 2, 800, "f.F90", 1, 800))
        assert free2.frees_onheap(free1) is True

    def test_frees_onheap_size_mismatch(self):
        free1 = Entry.from_line(make_line("v", "D", 1, -800, "f.F90", 1, 0))
        free2 = Entry.from_line(make_line("v", "D", 2, 400, "f.F90", 1, 400))
        assert free2.frees_onheap(free1) is False

    def test_frees_onstack_can_never_match_a_free_to_its_own_allocation(self):
        # BUG (documented, not fixed here): locus embeds the entry's own
        # action letter ("A:..." vs "D:..."), so a free's locus can never
        # equal the locus of the allocation it deallocates -- even for the
        # exact same vname/file/line. frees_onstack() (and transitively
        # find_memleaks()'s stack-reconciliation branch) can therefore never
        # actually match a free against its own allocation.
        alloc = Entry.from_line(make_line("v", "A", 1, 800, "f.F90", 1, 800))
        free_same_line = Entry.from_line(make_line("v", "D", 1, -800, "f.F90", 1, 0))
        assert alloc.locus != free_same_line.locus
        assert free_same_line.frees_onstack(alloc) is False

    def test_frees_onstack_matches_two_free_entries_with_identical_locus(self):
        # The only way frees_onstack() can return True in practice: two "D"
        # entries sharing the exact same locus (so parked under the same
        # stack[] key in find_memleaks()) whose sizes cancel.
        free1 = Entry.from_line(make_line("v", "D", 1, -800, "f.F90", 1, 0))
        free2 = Entry.from_line(make_line("v", "D", 2, 800, "f.F90", 1, 800))
        assert free1.locus == free2.locus
        assert free2.frees_onstack(free1) is True

    def test_frees_onstack_still_requires_matching_line(self):
        free1 = Entry.from_line(make_line("v", "D", 1, -800, "f.F90", 1, 0))
        free2 = Entry.from_line(make_line("v", "D", 2, 800, "f.F90", 2, 800))
        assert free2.frees_onstack(free1) is False


class TestEntriesToDataframe:
    """Test entries_to_dataframe()."""

    def test_converts_entries_to_dataframe_rows(self):
        entries = [
            Entry.from_line(make_line("v1", "A", 1, 8, "f.F90", 1, 8)),
            Entry.from_line(make_line("v2", "A", 2, 16, "f.F90", 2, 24)),
        ]
        df = entries_to_dataframe(entries)
        assert len(df) == 2
        assert list(df["vname"]) == ["v1", "v2"]
        assert "size_mb" in df.columns


ABIMEM_CONTENT = "".join(
    [
        "# header line, must be skipped\n",
        # Clean alloc/free pair on var_a.
        make_line("var_a", "A", 1000, 800, "mod_a.F90", 10, 800),
        make_line("var_a", "D", 1000, -800, "mod_a.F90", 20, 0),
        # A large, never-freed allocation (a genuine "leak" left on the heap).
        make_line("var_leak", "A", 2000, 20 * 8 * 1024 * 1024, "mod_b.F90", 5, 20 * 8 * 1024 * 1024),
        # A small allocation (below find_small_allocs()'s default threshold).
        make_line("var_small", "A", 3000, 80, "mod_c.F90", 1, 20 * 8 * 1024 * 1024 + 80),
        # A zero-sized alloc/free pair.
        make_line("var_zero", "A", 4000, 0, "mod_d.F90", 1, 20 * 8 * 1024 * 1024 + 80),
        make_line("var_zero", "D", 4000, 0, "mod_d.F90", 2, 20 * 8 * 1024 * 1024 + 80),
        # A weird (non-positive) pointer.
        make_line("var_weird", "A", 0, 8, "mod_e.F90", 1, 20 * 8 * 1024 * 1024 + 88),
        # A free with no matching prior allocation (a "reallocation").
        make_line("var_realloc", "D", 5000, -8, "mod_f.F90", 1, 20 * 8 * 1024 * 1024 + 80),
        # Two entries at the very same locus (repeated calls in a loop),
        # to exercise accumulated_entries()/get_intense_dataframe().
        make_line("var_loop", "A", 6000, 8, "mod_g.F90", 1, 20 * 8 * 1024 * 1024 + 88),
        make_line("var_loop", "A", 6000, 8, "mod_g.F90", 1, 20 * 8 * 1024 * 1024 + 96),
    ]
)


@pytest.fixture
def abimem_file(tmp_path):
    p = tmp_path / "abimem.mocc"
    p.write_text(ABIMEM_CONTENT)
    return AbimemFile(str(p))


class TestAbimemFileParsing:
    """Test AbimemFile.all_entries and accumulated_entries."""

    def test_all_entries_skips_header_and_parses_every_data_line(self, abimem_file):
        entries = abimem_file.all_entries
        assert len(entries) == 10
        assert all(isinstance(e, Entry) for e in entries)

    def test_all_entries_reports_parse_errors(self, tmp_path, capsys):
        bad_file = tmp_path / "bad.mocc"
        bad_file.write_text(make_line("v", "A", 1, 8, "f.F90", 1, 8) + "not a valid abimem line\n")
        analysis = AbimemFile(str(bad_file))
        with pytest.raises(Exception):
            analysis.all_entries
        assert "Error while parsing lineno" in capsys.readouterr().out

    def test_accumulated_entries_sums_same_locus(self, abimem_file):
        acc = abimem_file.accumulated_entries
        loop_entry = next(e for e in acc if e.vname == "var_loop")
        assert loop_entry.size == 16  # 8 + 8
        assert loop_entry.tot_memory == max(
            20 * 8 * 1024 * 1024 + 88, 20 * 8 * 1024 * 1024 + 96
        )


class TestAbimemFileAnalysis:
    """Test AbimemFile's alloc-analysis helper methods."""

    def test_find_small_allocs(self, abimem_file):
        smalls = abimem_file.find_small_allocs(nbits=160 * 8)
        assert {e.vname for e in smalls} == {"var_a", "var_small", "var_weird", "var_loop"}

    def test_find_large_allocs(self, abimem_file):
        larges = abimem_file.find_large_allocs(nbits=10 * 8 * 1024 * 1024)
        assert [e.vname for e in larges] == ["var_leak"]

    def test_find_zerosized(self, abimem_file):
        zeros = abimem_file.find_zerosized()
        assert {e.vname for e in zeros} == {"var_zero"}

    def test_find_zerosized_as_dataframe(self, abimem_file):
        df = abimem_file.find_zerosized(as_dataframe=True)
        assert list(df["vname"]) == ["var_zero", "var_zero"]

    def test_find_weird_ptrs_reports_non_positive_pointers(self, abimem_file, capsys):
        weird = abimem_file.find_weird_ptrs()
        assert {e.vname for e in weird} == {"var_weird"}
        assert "Found 1 weird entries" in capsys.readouterr().out

    def test_find_weird_ptrs_reports_none_found(self, tmp_path, capsys):
        clean = tmp_path / "clean.mocc"
        clean.write_text(make_line("v", "A", 1, 8, "f.F90", 1, 8))
        assert AbimemFile(str(clean)).find_weird_ptrs() == []
        assert "No weird entries found" in capsys.readouterr().out

    def test_get_intense_dataframe_groups_by_locus(self, abimem_file):
        df = abimem_file.get_intense_dataframe()
        loop_locus = next(e for e in abimem_file.all_entries if e.vname == "var_loop").locus
        assert df.loc[loop_locus, "ncalls"] == 2

    def test_dataframe_and_dataframe_accumulated(self, abimem_file):
        assert len(abimem_file.dataframe) == 10
        assert len(abimem_file.dataframe_accumulated) == 9  # var_loop's pair merged into 1

    def test_get_peaks_returns_largest_allocations_first(self, abimem_file):
        peaks = abimem_file.get_peaks(maxlen=5)
        vnames = [e.vname for e in peaks]
        assert vnames[0] == "var_leak"  # the single largest allocation

    def test_get_peaks_as_dataframe(self, abimem_file):
        df = abimem_file.get_peaks(maxlen=5, as_dataframe=True)
        assert "var_leak" in list(df["vname"])

    def test_get_hotspots_dataframe(self, abimem_file):
        df = abimem_file.get_hotspots_dataframe(accumulated=False)
        assert "mod_b.F90" in df.index
        assert df.loc["mod_b.F90", "malloc_mb"] == pytest.approx(20.0)


class TestFindMemleaks:
    """Test AbimemFile.find_memleaks()'s heap/stack/reallocation bookkeeping."""

    def test_clean_pair_leaves_nothing_behind(self, tmp_path):
        content = (
            make_line("v", "A", 1, 800, "f.F90", 1, 800)
            + make_line("v", "D", 1, -800, "f.F90", 2, 0)
        )
        analysis = AbimemFile(str(_write(tmp_path, content)))
        assert analysis.find_memleaks() == 0

    def test_unmatched_allocation_is_a_leak(self, tmp_path):
        content = make_line("v", "A", 1, 800, "f.F90", 1, 800)
        analysis = AbimemFile(str(_write(tmp_path, content)))
        assert analysis.find_memleaks() == 1

    def test_zero_sized_entries_are_ignored(self, tmp_path):
        content = make_line("v", "A", 1, 0, "f.F90", 1, 0)
        analysis = AbimemFile(str(_write(tmp_path, content)))
        assert analysis.find_memleaks() == 0

    def test_free_without_prior_allocation_is_a_realloc(self, tmp_path):
        content = make_line("v", "D", 1, -800, "f.F90", 1, 0)
        analysis = AbimemFile(str(_write(tmp_path, content)))
        assert analysis.find_memleaks() == 1

    def test_ptr_reused_ends_up_parked_on_the_stack_unreconciled(self, tmp_path):
        # Same ptr allocated twice in a row (simulating a compiler
        # reallocation): the second alloc can't cleanly pop the first from
        # the heap, so it's parked on the stack under its own "A:..." locus
        # key; a later free at that ptr that doesn't cancel the *original*
        # heap entry either likewise gets parked, under a "D:..." key.
        # Per frees_onstack()'s locus-includes-action quirk (see
        # TestEntry.test_frees_onstack_can_never_match_a_free_to_its_own_allocation),
        # these two parked entries can never reconcile with each other, so
        # both linger: heap keeps the untouched original alloc (1), and the
        # stack keeps both parked entries under 2 distinct locus keys (2).
        content = (
            make_line("v", "A", 1, 800, "f.F90", 1, 800)
            + make_line("v", "A", 1, 400, "f.F90", 1, 1200)
            + make_line("v", "D", 1, -400, "f.F90", 1, 800)
        )
        analysis = AbimemFile(str(_write(tmp_path, content)))
        assert analysis.find_memleaks() == 3

    def test_verbose_mode_prints_diagnostics(self, tmp_path, capsys):
        content = (
            make_line("v", "A", 1, 800, "f.F90", 1, 800)
            + make_line("v", "A", 1, 400, "f.F90", 1, 1200)
        )
        analysis = AbimemFile(str(_write(tmp_path, content)))
        analysis.find_memleaks(verbose=1)
        assert "WARNING:" in capsys.readouterr().out

    def test_verbose_mode_reports_reallocs(self, tmp_path, capsys):
        content = make_line("v", "D", 1, -800, "f.F90", 1, 0)
        analysis = AbimemFile(str(_write(tmp_path, content)))
        analysis.find_memleaks(verbose=1)
        assert "Possible reallocations:" in capsys.readouterr().out


def _write(tmp_path, content: str):
    p = tmp_path / "abimem.mocc"
    p.write_text(content)
    return p


class TestHeapAndStack:
    """Test the Heap/Stack dict subclasses."""

    def test_heap_show_empty(self, capsys):
        Heap().show()
        assert "HEAP OF LEN 0" in capsys.readouterr().out

    def test_heap_show_nonempty(self, capsys):
        heap = Heap()
        heap[1] = ["fake entry"]
        heap.show()
        out = capsys.readouterr().out
        assert "HEAP OF LEN 1" in out

    def test_heap_pop_alloc_non_free_entry_returns_zero(self):
        alloc = Entry.from_line(make_line("v", "A", 1, 8, "f.F90", 1, 8))
        assert Heap().pop_alloc(alloc) == 0

    def test_heap_pop_alloc_missing_ptr_returns_zero(self):
        free = Entry.from_line(make_line("v", "D", 1, 8, "f.F90", 1, 0))
        assert Heap().pop_alloc(free) == 0

    def test_heap_pop_alloc_raises_on_nonempty_match(self):
        # BUG (documented, not fixed here; also unused anywhere in the
        # module): `for i, olde in elist:` iterates the *Entry* namedtuples
        # themselves (7 fields each), trying to unpack each one into just
        # two names (i, olde) -- it should be `enumerate(elist)`, matching
        # the equivalent loop in find_memleaks(). The one case this method
        # exists for (a non-empty elist) always raises instead of matching.
        alloc = Entry.from_line(make_line("v", "A", 1, 800, "f.F90", 1, 800))
        free = Entry.from_line(make_line("v", "D", 1, -800, "f.F90", 1, 0))
        heap = Heap()
        heap[1] = [alloc]
        with pytest.raises(ValueError, match="too many values to unpack"):
            heap.pop_alloc(free)

    def test_stack_show_empty(self, capsys):
        Stack().show()
        assert "STACK OF LEN 0" in capsys.readouterr().out

    def test_stack_show_nonempty(self, capsys):
        stack = Stack()
        stack["locus"] = ["fake entry"]
        stack.show()
        out = capsys.readouterr().out
        assert "STACK OF LEN 1" in out


class TestPlotMethods:
    """Light smoke tests for the matplotlib-based plot_* methods."""

    def test_plot_memory_usage_returns_a_figure(self, abimem_file):
        fig = abimem_file.plot_memory_usage(show=False)
        assert fig is not None

    def test_plot_peaks_returns_a_figure(self, abimem_file):
        fig = abimem_file.plot_peaks(show=False, maxlen=5)
        assert fig is not None

    def test_plot_hist_returns_a_figure(self, abimem_file):
        fig = abimem_file.plot_hist(show=False)
        assert fig is not None

    def test_to_string_and_str(self, abimem_file):
        s = str(abimem_file)
        assert "var_leak" in s
        assert abimem_file.to_string() == s

    def test_expose_shows_the_memory_usage_plot(self, abimem_file):
        # AbimemFile.expose() itself isn't pragma-excluded (only the
        # MplExpose class it delegates to is); under the Agg backend
        # plt.show() is a no-op, so this just exercises the wiring.
        abimem_file.expose()
