"""
Unit tests for tools.py module.

Tests cover the small file/shell utility functions (patch, unzip, touch,
tail_file, which), the numeric-parsing helpers, RestrictedShell's sandboxed
command execution, StringColorizer/stream_has_colours, prompt()/
user_wants_to_exit(), Editor, Patcher, pprint_table, the ascii_* banners, and
the lazy_property descriptor.
"""

from __future__ import annotations

import gzip
import io
import os
import stat

import pytest

from . import tools
from .tools import (
    Editor,
    Patcher,
    PatcherError,
    RestrictedShell,
    RShellError,
    StringColorizer,
    ascii_abinit,
    ascii_scream,
    ascii_wasp,
    lazy_property,
    patch,
    pprint_table,
    prompt,
    stream_has_colours,
    tail_file,
    touch,
    unzip,
    user_wants_to_exit,
    which,
)


class TestPatch:
    """Test the patch() Unix diff/patch wrapper."""

    def test_copies_fromfile_over_tofile_and_keeps_backup(self, tmp_path):
        fromfile = tmp_path / "from.txt"
        tofile = tmp_path / "to.txt"
        fromfile.write_text("new content\n")
        tofile.write_text("old content\n")

        retcode = patch(str(fromfile), str(tofile))

        assert retcode == 0
        assert tofile.read_text() == "new content\n"
        assert (tmp_path / "to.txt.orig").read_text() == "old content\n"

    def test_reverts_to_backup_when_copy_fails(self, tmp_path):
        # A nonexistent fromfile makes the "cp" step fail, which must trigger
        # the revert-from-backup path instead of leaving tofile corrupted.
        fromfile = tmp_path / "does_not_exist.txt"
        tofile = tmp_path / "to.txt"
        tofile.write_text("old content\n")

        with pytest.warns(UserWarning, match="reverting to original file"):
            retcode = patch(str(fromfile), str(tofile))

        assert retcode != 0
        assert tofile.read_text() == "old content\n"


class TestUnzip:
    """Test the unzip() gzip-decompression helper."""

    def test_decompresses_to_default_destination(self, tmp_path):
        gz_path = tmp_path / "data.txt.gz"
        with gzip.open(gz_path, "wb") as fh:
            fh.write(b"hello world")

        unzip(str(gz_path))

        assert (tmp_path / "data.txt").read_bytes() == b"hello world"

    def test_decompresses_to_explicit_destination(self, tmp_path):
        gz_path = tmp_path / "data.txt.gz"
        dest = tmp_path / "explicit.txt"
        with gzip.open(gz_path, "wb") as fh:
            fh.write(b"payload")

        unzip(str(gz_path), dest=str(dest))

        assert dest.read_bytes() == b"payload"

    def test_rejects_non_gz_filename(self):
        with pytest.raises(ValueError, match=r"should end with \.gz"):
            unzip("archive.tar")


class TestTouch:
    """Test the touch() Unix `touch` emulation."""

    def test_creates_missing_file(self, tmp_path):
        fname = tmp_path / "new_file.txt"
        assert not fname.exists()
        touch(str(fname))
        assert fname.exists()

    def test_updates_mtime_of_existing_file(self, tmp_path):
        fname = tmp_path / "existing.txt"
        fname.write_text("data")
        touch(str(fname), times=(1000000.0, 2000000.0))
        stat_result = os.stat(fname)
        assert stat_result.st_atime == 1000000.0
        assert stat_result.st_mtime == 2000000.0


class TestTailFile:
    """Test the tail_file() Unix `tail` emulation."""

    def test_returns_last_n_lines_as_string(self, tmp_path):
        fname = tmp_path / "log.txt"
        fname.write_text("\n".join(f"line {i}" for i in range(10)) + "\n")

        out = tail_file(str(fname), 3)

        assert out == "line 7\nline 8\nline 9\n"

    def test_returns_last_n_lines_as_list(self, tmp_path):
        fname = tmp_path / "log.txt"
        fname.write_text("\n".join(f"line {i}" for i in range(10)) + "\n")

        out = tail_file(str(fname), 2, aslist=True)

        assert out == ["line 8\n", "line 9\n"]

    def test_raises_runtimeerror_on_missing_file(self, tmp_path):
        with pytest.raises(RuntimeError):
            tail_file(str(tmp_path / "missing.txt"), 3)


class TestWhich:
    """Test the which() PATH-lookup helper."""

    def test_finds_executable_by_name(self):
        # "true" is a standard POSIX utility, present on every Unix PATH.
        assert which("true") is not None

    def test_returns_none_for_unknown_program(self):
        assert which("this_program_does_not_exist_xyz") is None

    def test_absolute_path_to_executable_is_returned_as_is(self, tmp_path):
        script = tmp_path / "myscript.sh"
        script.write_text("#!/bin/sh\necho hi\n")
        script.chmod(script.stat().st_mode | stat.S_IEXEC)

        assert which(str(script)) == str(script)

    def test_absolute_path_to_non_executable_returns_none(self, tmp_path):
        script = tmp_path / "not_executable.sh"
        script.write_text("echo hi\n")
        script.chmod(0o644)

        assert which(str(script)) is None


class TestRestrictedShell:
    """Test RestrictedShell's sandboxed cp/mv/touch command execution."""

    @pytest.fixture
    def shell(self, tmp_path):
        inp_dir = tmp_path / "inp"
        workdir = tmp_path / "work"
        psps_dir = tmp_path / "psps"
        for d in (inp_dir, workdir, psps_dir):
            d.mkdir()
        return RestrictedShell(str(inp_dir), str(workdir), str(psps_dir)), inp_dir, workdir, psps_dir

    def test_touch_creates_file_in_workdir(self, shell):
        rshell, inp_dir, workdir, psps_dir = shell
        rshell.execute("w_touch out.txt")
        assert (workdir / "out.txt").exists()
        assert rshell.exceptions == []

    def test_copy_from_input_dir_to_workdir(self, shell):
        # The two-letter prefix ("iw") pairs positionally with the two args
        # via zip(pre_s, args) -- "iw_cp src dest" means "copy src (from i)
        # to dest (in w)", *not* a per-argument "i_"/"w_" prefix on each token.
        rshell, inp_dir, workdir, psps_dir = shell
        (inp_dir / "src.txt").write_text("data")
        rshell.execute("iw_cp src.txt dest.txt")
        assert (workdir / "dest.txt").read_text() == "data"

    def test_move_from_psps_dir_to_workdir(self, shell):
        rshell, inp_dir, workdir, psps_dir = shell
        (psps_dir / "pseudo.psp8").write_text("psp data")
        rshell.execute("pw_mv pseudo.psp8 pseudo.psp8")
        assert (workdir / "pseudo.psp8").read_text() == "psp data"
        assert not (psps_dir / "pseudo.psp8").exists()

    def test_wrong_argument_count_is_recorded_as_exception(self, shell):
        rshell, *_ = shell
        result = rshell.execute("w_touch a b c")
        assert result is None
        assert len(rshell.exceptions) == 1
        assert isinstance(rshell.exceptions[0], RShellError)

    def test_unparseable_command_is_recorded_as_exception(self, shell):
        rshell, *_ = shell
        result = rshell.execute("nonsense")
        assert result is None
        assert len(rshell.exceptions) == 1

    def test_unknown_command_key_is_recorded_as_exception(self, shell):
        rshell, *_ = shell
        result = rshell.execute("w_rm somefile")
        assert result is None
        assert len(rshell.exceptions) == 1

    def test_cp_falls_back_to_etsf_nc_variant(self, shell):
        rshell, inp_dir, workdir, psps_dir = shell
        (inp_dir / "out_DEN-etsf.nc").write_text("netcdf data")
        rshell.execute("iw_cp out_DEN out_DEN")
        assert (workdir / "out_DEN-etsf.nc").read_text() == "netcdf data"

    def test_cp_falls_back_to_plain_nc_variant(self, shell):
        rshell, inp_dir, workdir, psps_dir = shell
        (inp_dir / "out_DEN.nc").write_text("netcdf data")
        rshell.execute("iw_cp out_DEN out_DEN")
        assert (workdir / "out_DEN.nc").read_text() == "netcdf data"

    def test_cp_missing_source_and_no_netcdf_variant_is_recorded(self, shell):
        rshell, *_ = shell
        rshell.execute("iw_cp missing.txt missing.txt")
        assert len(rshell.exceptions) == 1

    def test_empty_exceptions_clears_recorded_errors(self, shell):
        rshell, *_ = shell
        rshell.execute("nonsense")
        assert rshell.exceptions
        rshell.empty_exceptions()
        assert rshell.exceptions == []

    def test_more_than_two_args_is_recorded_as_not_implemented(self, shell, monkeypatch):
        # None of the real _key2command entries (cp/mv: 2 args, touch: 1)
        # can reach the nargs>2 branch; register a fake one to exercise it.
        rshell, *_ = shell
        monkeypatch.setitem(RestrictedShell._key2command, "triple", (lambda a, b, c: None, 3))
        rshell.execute("iww_triple a b c")
        assert len(rshell.exceptions) == 1
        assert "too large" in str(rshell.exceptions[0])


class TestStreamHasColours:
    """Test stream_has_colours()'s TTY/curses detection heuristic."""

    def test_stream_without_isatty_returns_false(self):
        class NoIsAtty:
            pass

        assert stream_has_colours(NoIsAtty()) is False

    def test_non_tty_stream_returns_false(self):
        class NotATty:
            def isatty(self):
                return False

        assert stream_has_colours(NotATty()) is False

    def test_tty_stream_with_curses_failure_returns_false(self, monkeypatch):
        class Tty:
            def isatty(self):
                return True

        def broken_setupterm(*args, **kwargs):
            raise Exception("no terminfo")

        monkeypatch.setattr("curses.setupterm", broken_setupterm)
        assert stream_has_colours(Tty()) is False

    def test_tty_stream_with_curses_success_and_enough_colors_returns_true(self, monkeypatch):
        class Tty:
            def isatty(self):
                return True

        monkeypatch.setattr("curses.setupterm", lambda *a, **k: None)
        monkeypatch.setattr("curses.tigetnum", lambda name: 256)
        assert stream_has_colours(Tty()) is True


class TestStringColorizer:
    """Test StringColorizer's __call__ colorization."""

    class _FakeStream:
        def isatty(self):
            return False  # deterministic: never claims TTY colour support

    def test_no_colour_support_returns_plain_string(self):
        colorizer = StringColorizer(self._FakeStream())
        assert colorizer.has_colours is False
        assert colorizer("hello", "red") == "hello"

    def test_colour_support_wraps_known_colour(self, monkeypatch):
        monkeypatch.setattr(tools, "stream_has_colours", lambda stream: True)
        colorizer = StringColorizer(self._FakeStream())
        assert colorizer.has_colours is True
        assert colorizer("hello", "red") == "\x1b[01;31mhello\x1b[00m"

    def test_colour_support_with_default_colour_is_unchanged(self, monkeypatch):
        monkeypatch.setattr(tools, "stream_has_colours", lambda stream: True)
        colorizer = StringColorizer(self._FakeStream())
        assert colorizer("hello", "default") == "hello"

    def test_colour_support_with_unknown_colour_is_unchanged(self, monkeypatch):
        monkeypatch.setattr(tools, "stream_has_colours", lambda stream: True)
        colorizer = StringColorizer(self._FakeStream())
        assert colorizer("hello", "not_a_real_colour") == "hello"


class TestPromptAndUserWantsToExit:
    """Test prompt() and user_wants_to_exit()."""

    def test_prompt_delegates_to_input(self, monkeypatch):
        monkeypatch.setattr("builtins.input", lambda question: "42")
        assert prompt("Enter a number: ") == "42"

    @pytest.mark.parametrize("answer", ["n", "N", "no", "NO", " n "])
    def test_user_wants_to_exit_true_for_negative_answers(self, monkeypatch, answer):
        monkeypatch.setattr(tools, "prompt", lambda question: answer)
        assert user_wants_to_exit() is True

    @pytest.mark.parametrize("answer", ["y", "yes", "", "whatever"])
    def test_user_wants_to_exit_false_for_other_answers(self, monkeypatch, answer):
        monkeypatch.setattr(tools, "prompt", lambda question: answer)
        assert user_wants_to_exit() is False

    def test_user_wants_to_exit_true_on_eof(self, monkeypatch):
        def raise_eof(question):
            raise EOFError

        monkeypatch.setattr(tools, "prompt", raise_eof)
        assert user_wants_to_exit() is True


class TestEditor:
    """Test the Editor class."""

    def test_default_editor_from_env(self, monkeypatch):
        monkeypatch.setenv("EDITOR", "nano")
        assert Editor().editor == "nano"

    def test_default_editor_fallback_to_vi(self, monkeypatch):
        monkeypatch.delenv("EDITOR", raising=False)
        assert Editor().editor == "vi"

    def test_explicit_editor_overrides_env(self):
        assert Editor("emacs").editor == "emacs"

    def test_edit_file_without_lineno(self, monkeypatch):
        # edit_file() does its own `from subprocess import call` *inside* the
        # method body, which rebinds the local name straight from the
        # subprocess module -- patching tools.call would not be seen here.
        calls = []
        monkeypatch.setattr("subprocess.call", lambda args: calls.append(args) or 0)
        retcode = Editor("myeditor").edit_file("foo.txt")
        assert retcode == 0
        assert calls == [["myeditor", "foo.txt"]]

    def test_edit_file_with_lineno_appends_plus_line_syntax(self, monkeypatch):
        calls = []
        monkeypatch.setattr("subprocess.call", lambda args: calls.append(args) or 0)
        Editor("vi").edit_file("foo.txt", lineno=42)
        assert calls == [["vi", "foo.txt", "+42"]]

    def test_edit_file_warns_on_nonzero_retcode(self, monkeypatch):
        monkeypatch.setattr("subprocess.call", lambda args: 1)
        with pytest.warns(UserWarning, match="Error while trying to edit file"):
            retcode = Editor("vi").edit_file("foo.txt")
        assert retcode == 1

    def test_edit_files_stops_early_when_user_declines_to_continue(self, monkeypatch):
        edited = []
        monkeypatch.setattr(
            Editor, "edit_file", lambda self, fname, lineno=None: edited.append(fname) or 0
        )
        monkeypatch.setattr(tools, "user_wants_to_exit", lambda: True)
        Editor("vi").edit_files(["a.txt", "b.txt", "c.txt"])
        assert edited == ["a.txt"]

    def test_edit_files_without_asking_for_exit_edits_everything(self, monkeypatch):
        edited = []
        monkeypatch.setattr(
            Editor, "edit_file", lambda self, fname, lineno=None: edited.append(fname) or 0
        )
        Editor("vi").edit_files(["a.txt", "b.txt"], ask_for_exit=False)
        assert edited == ["a.txt", "b.txt"]


class TestPatcher:
    """Test the Patcher class."""

    def test_unknown_patcher_name_raises(self, monkeypatch):
        monkeypatch.setattr(tools, "which", lambda program: "/usr/bin/whatever")
        with pytest.raises(ValueError, match="is not supported"):
            Patcher("not_a_real_patcher")

    def test_missing_executable_raises(self, monkeypatch):
        monkeypatch.setattr(tools, "which", lambda program: None)
        with pytest.raises(ValueError, match="Cannot find executable"):
            Patcher("patch")

    def test_default_patcher_from_env(self, monkeypatch):
        monkeypatch.setenv("PATCHER", "kdiff3")
        monkeypatch.setattr(tools, "which", lambda program: "/usr/bin/kdiff3")
        assert Patcher().patcher == "kdiff3"

    def test_is_interactive(self, monkeypatch):
        monkeypatch.setattr(tools, "which", lambda program: "/usr/bin/whatever")
        assert Patcher("vimdiff").is_interactive is True
        assert Patcher("patch").is_interactive is False

    def test_patch_with_auto_patcher_delegates_to_module_patch(self, monkeypatch, tmp_path):
        monkeypatch.setattr(tools, "which", lambda program: "/usr/bin/patch")
        calls = []
        monkeypatch.setattr(tools, "patch", lambda fromfile, tofile: calls.append((fromfile, tofile)) or 0)
        patcher = Patcher("patch")
        assert patcher.patch("a", "b") == 0
        assert calls == [("a", "b")]

    def test_patch_wraps_exception_from_auto_patcher(self, monkeypatch):
        monkeypatch.setattr(tools, "which", lambda program: "/usr/bin/patch")

        def raising_patch(fromfile, tofile):
            raise RuntimeError("boom")

        monkeypatch.setattr(tools, "patch", raising_patch)
        patcher = Patcher("patch")
        with pytest.raises(PatcherError):
            patcher.patch("a", "b")

    def test_patch_with_interactive_patcher_calls_subprocess(self, monkeypatch):
        monkeypatch.setattr(tools, "which", lambda program: "/usr/bin/vimdiff")
        calls = []
        monkeypatch.setattr(tools, "call", lambda args: calls.append(args) or 0)
        patcher = Patcher("vimdiff")
        assert patcher.patch("a", "b") == 0
        assert calls == [["vimdiff", "a", "b"]]

    def test_patch_wraps_exception_from_interactive_patcher(self, monkeypatch):
        monkeypatch.setattr(tools, "which", lambda program: "/usr/bin/vimdiff")

        def raising_call(args):
            raise RuntimeError("boom")

        monkeypatch.setattr(tools, "call", raising_call)
        patcher = Patcher("vimdiff")
        with pytest.raises(PatcherError):
            patcher.patch("a", "b")

    def test_patch_files_with_no_files_is_a_noop(self, monkeypatch):
        monkeypatch.setattr(tools, "which", lambda program: "/usr/bin/patch")
        patcher = Patcher("patch")
        assert patcher.patch_files([], []) == 0

    def test_patch_files_mismatched_lengths_raises(self, monkeypatch):
        monkeypatch.setattr(tools, "which", lambda program: "/usr/bin/patch")
        patcher = Patcher("patch")
        with pytest.raises(AssertionError):
            patcher.patch_files(["a"], ["a", "b"])

    def test_patch_files_auto_mode_asks_for_confirmation_and_aborts_on_no(self, monkeypatch):
        monkeypatch.setattr(tools, "which", lambda program: "/usr/bin/patch")
        monkeypatch.setattr(tools, "prompt", lambda question: "n")
        patched = []
        monkeypatch.setattr(Patcher, "patch", lambda self, f, t: patched.append((f, t)) or 0)
        patcher = Patcher("patch")
        assert patcher.patch_files(["a"], ["b"]) == 0
        assert patched == []

    def test_patch_files_auto_mode_confirmed_patches_all_with_backup(self, monkeypatch, tmp_path):
        monkeypatch.setattr(tools, "which", lambda program: "/usr/bin/patch")
        monkeypatch.setattr(tools, "prompt", lambda question: "y")
        patched = []
        monkeypatch.setattr(Patcher, "patch", lambda self, f, t: patched.append((f, t)) or 0)

        to1 = tmp_path / "to1.txt"
        to2 = tmp_path / "to2.txt"
        to1.write_text("orig1")
        to2.write_text("orig2")

        patcher = Patcher("patch")
        status = patcher.patch_files(["from1", "from2"], [str(to1), str(to2)])

        assert status == 0
        assert patched == [("from1", str(to1)), ("from2", str(to2))]
        assert (tmp_path / "to1.txt.orig").read_text() == "orig1"
        assert (tmp_path / "to2.txt.orig").read_text() == "orig2"

    def test_patch_files_auto_mode_stops_on_first_failure(self, monkeypatch, tmp_path):
        monkeypatch.setattr(tools, "which", lambda program: "/usr/bin/patch")
        monkeypatch.setattr(tools, "prompt", lambda question: "y")
        patched = []

        def fake_patch(self, f, t):
            patched.append((f, t))
            return 1  # failure

        monkeypatch.setattr(Patcher, "patch", fake_patch)

        to1 = tmp_path / "to1.txt"
        to2 = tmp_path / "to2.txt"
        to1.write_text("orig1")
        to2.write_text("orig2")

        patcher = Patcher("patch")
        status = patcher.patch_files(["from1", "from2"], [str(to1), str(to2)])

        assert status == 1
        # Must stop after the first failure, not attempt the second file.
        assert patched == [("from1", str(to1))]

    def test_patch_files_interactive_mode_stops_early_when_user_declines(self, monkeypatch):
        monkeypatch.setattr(tools, "which", lambda program: "/usr/bin/vimdiff")
        monkeypatch.setattr(tools, "user_wants_to_exit", lambda: True)
        patched = []
        monkeypatch.setattr(Patcher, "patch", lambda self, f, t: patched.append((f, t)) or 0)

        patcher = Patcher("vimdiff")
        patcher.patch_files(["a", "b", "c"], ["x", "y", "z"])

        assert patched == [("a", "x")]


class TestPprintTable:
    """Test pprint_table()'s column-aligned text output."""

    def test_basic_alignment(self):
        out = io.StringIO()
        pprint_table([["a", "1"], ["bb", "22"]], out=out)
        lines = out.getvalue().splitlines()
        # Column 0 is ljust'd to (max width + 1); column 1 is rjust'd to
        # (max width + 2) -- max widths here are 2 ("bb") and 2 ("22").
        assert lines[0] == "a  " + "   1"
        assert lines[1] == "bb " + "  22"

    def test_rstrip_removes_trailing_whitespace_from_cells(self):
        out = io.StringIO()
        pprint_table([["a  ", "1"], ["b", "2"]], out=out, rstrip=True)
        lines = out.getvalue().splitlines()
        # After rstrip, both first-column entries are a single character wide.
        assert lines[0] == "a " + "  1"
        assert lines[1] == "b " + "  2"


class TestAsciiArt:
    """Trivial smoke tests for the ascii_* banner functions."""

    def test_ascii_wasp_returns_nonempty_string(self):
        assert isinstance(ascii_wasp(), str)
        assert ascii_wasp().strip()

    def test_ascii_scream_returns_nonempty_string(self):
        assert isinstance(ascii_scream(), str)
        assert ascii_scream().strip()

    def test_ascii_abinit_returns_nonempty_string(self):
        assert isinstance(ascii_abinit(), str)
        assert ascii_abinit().strip()


class TestLazyProperty:
    """Test the lazy_property descriptor."""

    def test_value_is_computed_once_and_cached(self):
        calls = []

        class Foo:
            @lazy_property
            def value(self):
                calls.append(1)
                return 42

        foo = Foo()
        assert foo.value == 42
        assert foo.value == 42
        assert len(calls) == 1
        assert foo.__dict__["value"] == 42

    def test_class_level_access_returns_the_descriptor(self):
        class Foo:
            @lazy_property
            def value(self):
                return 42

        assert isinstance(Foo.__dict__["value"], lazy_property)
        assert isinstance(Foo.value, lazy_property)

    def test_dunder_style_name_is_mangled(self):
        class Foo:
            @lazy_property
            def __secret(self):
                return "hidden"

        foo = Foo()
        # __get__ mangles a leading-dunder, non-dunder-suffixed name to
        # _ClassName__attr before caching it in the instance __dict__.
        assert foo._Foo__secret == "hidden"
        assert "_Foo__secret" in foo.__dict__

    def test_invalidate_removes_cached_value(self):
        class Foo:
            @lazy_property
            def value(self):
                return 42

        foo = Foo()
        assert foo.value == 42
        lazy_property.invalidate(foo, "value")
        assert "value" not in foo.__dict__

    def test_invalidate_on_non_lazy_attribute_raises(self):
        class Foo:
            value = 42

        foo = Foo()
        with pytest.raises(AttributeError, match="is not a"):
            lazy_property.invalidate(foo, "value")

    def test_invalidate_on_object_without_dict_raises(self):
        class Foo:
            __slots__ = ()

            @lazy_property
            def value(self):
                return 42

        with pytest.raises(AttributeError, match="has no attribute '__dict__'"):
            lazy_property.invalidate(Foo(), "value")

    def test_get_on_object_without_dict_raises(self):
        # Distinct from the invalidate() check above: this exercises
        # __get__'s own "no __dict__" guard, hit on first attribute access
        # rather than on an explicit invalidate() call.
        class Foo:
            __slots__ = ()

            @lazy_property
            def value(self):
                return 42

        with pytest.raises(AttributeError, match="has no attribute '__dict__'"):
            Foo().value

    def test_invalidate_on_dunder_style_name_is_mangled(self):
        class Foo:
            @lazy_property
            def __secret(self):
                return "hidden"

        foo = Foo()
        assert foo._Foo__secret == "hidden"
        lazy_property.invalidate(foo, "__secret")
        assert "_Foo__secret" not in foo.__dict__
