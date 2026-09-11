"""
Unit tests for termcolor.py module.

Tests cover the enable/ison toggle, colored()/colored_map() text formatting
(color, background, attributes), cprint()/cprint_map() printing, stream
capability detection, and terminal size detection.
"""

from __future__ import annotations

import pytest

from . import termcolor
from .termcolor import (
    ATTRIBUTES,
    COLORS,
    HIGHLIGHTS,
    RESET,
    cprint,
    cprint_map,
    colored,
    colored_map,
    enable,
    get_terminal_size,
    ison,
    stream_has_colours,
)


@pytest.fixture(autouse=True)
def _reset_state(monkeypatch):
    """Ensure colorization is on and not disabled via env var for every test."""
    monkeypatch.delenv("ANSI_COLORS_DISABLED", raising=False)
    enable(True)
    yield
    enable(True)


class TestEnableIson:
    """Test the enable()/ison() global on/off switch."""

    def test_default_is_on(self):
        assert ison() is True

    def test_enable_false_turns_off(self):
        enable(False)
        assert ison() is False

    def test_enable_true_turns_back_on(self):
        enable(False)
        enable(True)
        assert ison() is True


class TestColored:
    """Test the colored() text-formatting function."""

    def test_plain_text_when_no_color_requested(self):
        # No color/on_color/attrs: text is still wrapped with RESET when ison().
        assert colored("hello") == "hello" + RESET

    def test_color_wraps_text_with_color_code(self):
        result = colored("hello", color="red")
        assert result == "\033[%dmhello%s" % (COLORS["red"], RESET)

    def test_on_color_wraps_text_with_highlight_code(self):
        result = colored("hello", on_color="on_blue")
        assert result == "\033[%dmhello%s" % (HIGHLIGHTS["on_blue"], RESET)

    def test_single_attribute(self):
        result = colored("hello", attrs=["bold"])
        assert result == "\033[%dmhello%s" % (ATTRIBUTES["bold"], RESET)

    def test_color_on_color_and_attrs_combine_in_order(self):
        result = colored("hi", color="green", on_color="on_red", attrs=["bold", "underline"])
        expected = "hi"
        expected = "\033[%dm%s" % (COLORS["green"], expected)
        expected = "\033[%dm%s" % (HIGHLIGHTS["on_red"], expected)
        expected = "\033[%dm%s" % (ATTRIBUTES["bold"], expected)
        expected = "\033[%dm%s" % (ATTRIBUTES["underline"], expected)
        expected += RESET
        assert result == expected

    def test_unknown_color_raises_keyerror(self):
        with pytest.raises(KeyError):
            colored("hello", color="not_a_color")

    def test_disabled_via_enable_returns_plain_text(self):
        enable(False)
        assert colored("hello", color="red", attrs=["bold"]) == "hello"

    def test_disabled_via_env_var_returns_plain_text(self, monkeypatch):
        monkeypatch.setenv("ANSI_COLORS_DISABLED", "1")
        assert colored("hello", color="red") == "hello"


class TestColoredMap:
    """Test colored_map(), which colorizes tokens found in a larger string."""

    def test_replaces_matching_token_with_color(self):
        result = colored_map("foo bar", {"bar": "green"})
        assert result == "foo " + colored("bar", color="green")

    def test_dict_value_supports_color_and_on_color(self):
        result = colored_map("foo bar", {"bar": {"color": "green", "on_color": "on_red"}})
        assert result == "foo " + colored("bar", color="green", on_color="on_red")

    def test_no_matching_token_leaves_text_untouched(self):
        assert colored_map("foo bar", {"baz": "green"}) == "foo bar"

    def test_disabled_returns_original_text_untouched(self):
        enable(False)
        # colored_map() checks the module-level switch itself (not colored()'s
        # ANSI_COLORS_DISABLED path), so no escape codes at all -- not even a
        # trailing RESET -- should appear.
        text = "foo bar"
        assert colored_map(text, {"bar": "green"}) == text


class TestCprint:
    """Test cprint()/cprint_map() actually print the colorized text."""

    def test_cprint_writes_colored_text(self, capsys):
        cprint("hello", "red")
        captured = capsys.readouterr()
        assert captured.out == colored("hello", "red") + "\n"

    def test_cprint_passes_through_print_kwargs(self, capsys):
        cprint("hello", "red", end="")
        captured = capsys.readouterr()
        assert captured.out == colored("hello", "red")

    def test_cprint_map_writes_colorized_text(self, capsys):
        cprint_map("Hello world", {"Hello": "red"})
        captured = capsys.readouterr()
        assert captured.out == colored_map("Hello world", {"Hello": "red"}) + "\n"

    def test_cprint_retries_without_flush_on_typeerror(self, monkeypatch, capsys):
        # Exercises the "flush is not supported by py2.7" fallback: the first
        # print() call is made to fail exactly once (as a real TypeError
        # would for an unsupported kwarg), and the retry must still print
        # the colorized text with `flush` dropped, not raise or print nothing.
        calls = []
        real_print = print

        def flaky_print(*args, **kwargs):
            calls.append(kwargs)
            if len(calls) == 1:
                raise TypeError("flush not supported")
            real_print(*args, **kwargs)

        monkeypatch.setattr(termcolor, "print", flaky_print, raising=False)
        cprint("hello", "red", flush=True)
        captured = capsys.readouterr()
        assert captured.out == colored("hello", "red") + "\n"
        assert "flush" not in calls[1]

    def test_cprint_map_retries_without_flush_on_typeerror(self, monkeypatch, capsys):
        calls = []
        real_print = print

        def flaky_print(*args, **kwargs):
            calls.append(kwargs)
            if len(calls) == 1:
                raise TypeError("flush not supported")
            real_print(*args, **kwargs)

        monkeypatch.setattr(termcolor, "print", flaky_print, raising=False)
        cprint_map("Hello world", {"Hello": "red"}, flush=True)
        captured = capsys.readouterr()
        assert captured.out == colored_map("Hello world", {"Hello": "red"}) + "\n"
        assert "flush" not in calls[1]


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

        # Force the "guess false in case of error" branch deterministically:
        # whether curses.setupterm() itself succeeds depends on the $TERM
        # of whatever environment runs this test, so don't rely on that.
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

    def test_tty_stream_with_too_few_colors_returns_false(self, monkeypatch):
        class Tty:
            def isatty(self):
                return True

        monkeypatch.setattr("curses.setupterm", lambda *a, **k: None)
        monkeypatch.setattr("curses.tigetnum", lambda name: 2)
        assert stream_has_colours(Tty()) is False


class TestGetTerminalSize:
    """Test get_terminal_size()'s cascading fallbacks."""

    def test_uses_stty_size_when_available(self, monkeypatch):
        class FakePopen:
            def read(self):
                return "40 100\n"

        monkeypatch.setattr(termcolor.os, "popen", lambda cmd, mode="r": FakePopen())
        assert get_terminal_size() == (40, 100)

    def test_falls_back_to_env_vars_when_all_probes_fail(self, monkeypatch):
        def broken_popen(cmd, mode="r"):
            raise OSError("no stty")

        monkeypatch.setattr(termcolor.os, "popen", broken_popen)

        # ioctl_GWINSZ() is a nested function re-importing fcntl/termios on
        # every call, so patching the real fcntl.ioctl (already loaded in
        # sys.modules) also fails it for fds 0/1/2 *and* the ctermid fd --
        # this must not depend on whether fd 0/1/2 happen to be a real TTY
        # in whatever environment runs this test.
        def broken_ioctl(fd, request, buf):
            raise OSError("not a tty")

        monkeypatch.setattr("fcntl.ioctl", broken_ioctl)

        def broken_ctermid():
            raise OSError("no ctermid")

        monkeypatch.setattr(termcolor.os, "ctermid", broken_ctermid)
        monkeypatch.setenv("LINES", "50")
        monkeypatch.setenv("COLUMNS", "120")
        assert get_terminal_size() == (50, 120)

    def test_falls_back_to_ioctl_when_stty_fails(self, monkeypatch):
        import struct

        def broken_popen(cmd, mode="r"):
            raise OSError("no stty")

        monkeypatch.setattr(termcolor.os, "popen", broken_popen)
        monkeypatch.setattr("fcntl.ioctl", lambda fd, request, buf: struct.pack("hh", 33, 111))
        assert get_terminal_size() == (33, 111)
