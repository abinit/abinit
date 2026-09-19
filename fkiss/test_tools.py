"""
Unit tests for tools.py module.
"""

from __future__ import annotations

import io
import pickle
import stat

import pandas as pd
import pytest

from . import tools
from .tools import Editor, NotebookWriter, lazy_property, pprint_table, print_dataframe, prompt, user_wants_to_exit, which


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
        assert foo._Foo__secret == "hidden"
        assert "_Foo__secret" in foo.__dict__

    def test_get_on_object_without_dict_raises(self):
        class Foo:
            __slots__ = ()

            @lazy_property
            def value(self):
                return 42

        with pytest.raises(AttributeError, match="has no attribute '__dict__'"):
            Foo().value

    def test_invalidate_removes_cached_value(self):
        class Foo:
            @lazy_property
            def value(self):
                return 42

        foo = Foo()
        assert foo.value == 42
        lazy_property.invalidate(foo, "value")
        assert "value" not in foo.__dict__

    def test_invalidate_on_dunder_style_name_is_mangled(self):
        class Foo:
            @lazy_property
            def __secret(self):
                return "hidden"

        foo = Foo()
        assert foo._Foo__secret == "hidden"
        lazy_property.invalidate(foo, "__secret")
        assert "_Foo__secret" not in foo.__dict__

    def test_invalidate_on_non_lazy_attribute_raises(self):
        class Foo:
            value = 42

        with pytest.raises(AttributeError, match="is not a"):
            lazy_property.invalidate(Foo(), "value")

    def test_invalidate_on_object_without_dict_raises(self):
        class Foo:
            __slots__ = ()

            @lazy_property
            def value(self):
                return 42

        with pytest.raises(AttributeError, match="has no attribute '__dict__'"):
            lazy_property.invalidate(Foo(), "value")


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
        assert lines[0] == "a " + "  1"
        assert lines[1] == "b " + "  2"


class TestWhich:
    """Test the which() PATH-lookup helper."""

    def test_finds_executable_by_name(self):
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


class TestPrintDataframe:
    """Test print_dataframe()'s formatted output."""

    def test_returns_string_when_file_is_string(self):
        df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
        out = print_dataframe(df, file="string")
        assert isinstance(out, str)
        assert "a" in out and "b" in out

    def test_title_is_included(self):
        df = pd.DataFrame({"a": [1, 2]})
        out = print_dataframe(df, title="My Title", file="string")
        assert "My Title" in out

    def test_writes_to_a_real_file_object(self):
        df = pd.DataFrame({"a": [1, 2]})
        buf = io.StringIO()
        result = print_dataframe(df, file=buf)
        assert result is None
        assert "a" in buf.getvalue()

    def test_sortby_reorders_rows(self):
        df = pd.DataFrame({"a": [3, 1, 2]})
        out = print_dataframe(df, sortby="a", file="string")
        # Sorted ascending: 1 must appear before 2, which must appear before 3.
        assert out.index("1") < out.index("2") < out.index("3")

    def test_sortby_ignored_when_column_missing(self):
        df = pd.DataFrame({"a": [3, 1, 2]})
        # Must not raise even though "missing" isn't a column.
        out = print_dataframe(df, sortby="missing", file="string")
        assert isinstance(out, str)

    def test_display_html_returns_ipython_html_object(self):
        # IPython is a soft/optional dependency (only display="html" needs
        # it) -- skip rather than fail where it isn't installed, same
        # convention as TestNotebookWriter's abipy/nbformat skips below.
        pytest.importorskip("IPython.core.display", exc_type=ImportError)
        from IPython.core.display import HTML

        df = pd.DataFrame({"a": [1, 2]})
        result = print_dataframe(df, display="html")
        assert isinstance(result, HTML)
        assert "<table" in result.data or "<div" in result.data


class ConcreteNotebookWriter(NotebookWriter):
    """Minimal concrete subclass for exercising the NotebookWriter mixin."""

    def write_notebook(self, nbpath=None):
        nbformat, nbv, nb = self.get_nbformat_nbv_nb(title="Test notebook")
        return self._write_nb_nbpath(nb, nbpath)

    def yield_figs(self, **kwargs):
        return iter([])


class TestNotebookWriter:
    """Test the NotebookWriter mixin's pickle and notebook-building helpers."""

    def test_pickle_roundtrip_with_explicit_path(self, tmp_path):
        obj = ConcreteNotebookWriter()
        path = str(tmp_path / "obj.pickle")
        result_path = obj.pickle_dump(path)
        assert result_path == path

        loaded = ConcreteNotebookWriter.pickle_load(path)
        assert isinstance(loaded, ConcreteNotebookWriter)

    def test_pickle_dump_creates_tempfile_when_no_path_given(self):
        obj = ConcreteNotebookWriter()
        path = obj.pickle_dump()
        try:
            with open(path, "rb") as fh:
                assert isinstance(pickle.load(fh), ConcreteNotebookWriter)
        finally:
            import os

            os.remove(path)

    def test_get_nbformat_nbv(self):
        # nbformat is a soft/optional dependency of fkiss (only the
        # notebook-building helpers need it) -- skip rather than fail where
        # it isn't installed, same convention as the abipy skips above.
        pytest.importorskip("nbformat", exc_type=ImportError)
        obj = ConcreteNotebookWriter()
        nbformat, nbv = obj.get_nbformat_nbv()
        assert nbv.__name__.endswith("v4")

    def test_get_nbformat_nbv_nb_adds_title_cell(self):
        pytest.importorskip("nbformat", exc_type=ImportError)
        obj = ConcreteNotebookWriter()
        nbformat, nbv, nb = obj.get_nbformat_nbv_nb(title="Hello")
        assert any("Hello" in cell.get("source", "") for cell in nb.cells)

    def test_write_notebook_creates_a_valid_ipynb_file(self, tmp_path):
        pytest.importorskip("nbformat", exc_type=ImportError)
        obj = ConcreteNotebookWriter()
        nbpath = str(tmp_path / "test.ipynb")
        result = obj.write_notebook(nbpath=nbpath)
        assert result == nbpath

        import nbformat

        nb = nbformat.read(nbpath, as_version=4)
        assert any("Test notebook" in cell.get("source", "") for cell in nb.cells)

    def test_write_notebook_defaults_to_a_tempfile(self):
        pytest.importorskip("nbformat", exc_type=ImportError)
        obj = ConcreteNotebookWriter()
        nbpath = obj.write_notebook()
        try:
            assert nbpath.endswith(".ipynb")
            assert "abinb_" in nbpath
        finally:
            import os

            os.remove(nbpath)

    def test_expose_with_no_figures_does_not_raise(self):
        # Regression test: expose() used to import the name `MplExpose` from
        # abipy.tools.plotting, but abipy had renamed it to `MplExposer`, so
        # every call raised ImportError instead of showing figures; it also
        # passed slide_mode as the slide_timeout kwarg by copy-paste mistake.
        #
        # abipy is a soft/optional dependency of fkiss (only expose() needs
        # it), so skip rather than fail where it isn't installed.
        pytest.importorskip("abipy.tools.plotting", exc_type=ImportError)
        import matplotlib

        matplotlib.use("Agg")
        obj = ConcreteNotebookWriter()
        obj.expose()

    def test_expose_passes_slide_timeout_through_correctly(self, monkeypatch):
        pytest.importorskip("abipy.tools.plotting", exc_type=ImportError)
        captured_kwargs = {}

        class FakeExposer:
            def __init__(self, **kwargs):
                captured_kwargs.update(kwargs)

            def __enter__(self):
                return lambda obj: None

            def __exit__(self, *exc_info):
                return False

        monkeypatch.setattr("abipy.tools.plotting.MplExposer", FakeExposer)
        obj = ConcreteNotebookWriter()
        obj.expose(slide_mode=True, slide_timeout=30)

        assert captured_kwargs["slide_mode"] is True
        assert captured_kwargs["slide_timeout"] == 30
