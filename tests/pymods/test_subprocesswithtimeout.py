"""
Unit tests for subprocesswithtimeout.py module.

Tests cover input validation, a process that completes well before its
timeout, and a process that has to be killed once its timeout expires.
The module's own embedded unittest.TestCase (never collected by pytest,
since it lives outside a test_*.py file) is superseded by this file.
"""

from __future__ import annotations

import errno
import os
import time

import pytest

from . import subprocesswithtimeout as swt
from .subprocesswithtimeout import SubProcessWithTimeout


class TestInit:
    """Test SubProcessWithTimeout.__init__()'s argument validation."""

    def test_valid_timeout_and_delay(self):
        proc = SubProcessWithTimeout(timeout=2, delay=0.5)
        assert proc.timeout == 2.0
        assert proc.delay == 0.5

    def test_accepts_string_numbers(self):
        # __init__ coerces via float(), so numeric strings must also work.
        proc = SubProcessWithTimeout(timeout="2", delay="0.1")
        assert proc.timeout == 2.0
        assert proc.delay == 0.1

    def test_delay_greater_than_timeout_raises(self):
        with pytest.raises(ValueError, match="delay and timeout must be positive"):
            SubProcessWithTimeout(timeout=1, delay=2)

    def test_zero_delay_raises(self):
        with pytest.raises(ValueError, match="delay and timeout must be positive"):
            SubProcessWithTimeout(timeout=1, delay=0)

    def test_negative_timeout_raises(self):
        with pytest.raises(ValueError, match="delay and timeout must be positive"):
            SubProcessWithTimeout(timeout=-1, delay=0.1)

    def test_zero_timeout_raises(self):
        with pytest.raises(ValueError, match="delay and timeout must be positive"):
            SubProcessWithTimeout(timeout=0, delay=0.1)


class TestRun:
    """Test SubProcessWithTimeout.run()'s two outcomes: finishes vs. killed."""

    def test_process_completing_before_timeout_returns_its_own_exit_code(self):
        # "true" always exits 0 almost instantly, well under a 2s timeout.
        proc, retcode = SubProcessWithTimeout(timeout=2, delay=0.05).run(["true"])
        assert retcode == 0
        assert proc.returncode == 0

    def test_process_completing_before_timeout_propagates_nonzero_exit_code(self):
        proc, retcode = SubProcessWithTimeout(timeout=2, delay=0.05).run(["false"])
        assert retcode == 1

    def test_process_exceeding_timeout_is_killed(self):
        start = time.time()
        proc, retcode = SubProcessWithTimeout(timeout=0.2, delay=0.05).run(["sleep", "5"])
        elapsed = time.time() - start
        # 124 (SIGTERM) or 137 (SIGKILL): whichever fires depends on whether
        # the process group signal (see _wait_testcomplete) actually reaches
        # the child in this environment, but either way run() must return
        # promptly rather than waiting out the full "sleep 5".
        assert retcode in (124, 137)
        assert elapsed < 4

        # Don't leak the child process/zombie into the rest of the test session.
        if proc.poll() is None:
            proc.kill()
        proc.wait()

    def test_process_killed_by_real_sigterm_returns_124(self):
        # Plain subprocess.Popen() does not make the child a process-group
        # leader, so os.kill(-pid, ...) (targeting a *group*) silently no-ops
        # with ESRCH in the test above -- returning 137 regardless of whether
        # SIGTERM actually did anything. preexec_fn=os.setsid makes the child
        # its own session/group leader so the real SIGTERM reaches it, and
        # plain "sleep" dies on SIGTERM -- exercising the intended "died
        # cleanly on SIGTERM, no SIGKILL needed" return path.
        start = time.time()
        proc, retcode = SubProcessWithTimeout(timeout=0.2, delay=0.05).run(
            ["sleep", "5"], preexec_fn=os.setsid
        )
        elapsed = time.time() - start
        assert retcode == 124
        assert elapsed < 4
        proc.wait()

    def test_unexpected_oserror_from_sigterm_is_not_swallowed(self, monkeypatch):
        # _wait_testcomplete() deliberately swallows ESRCH (the process
        # already exited on its own) but must re-raise anything else -- a
        # real permissions error must not look like a harmless timeout.
        def fake_kill(pid, sig):
            raise OSError(errno.EPERM, "no permission")

        monkeypatch.setattr(swt.os, "kill", fake_kill)
        runner = SubProcessWithTimeout(timeout=0.2, delay=0.05)
        try:
            with pytest.raises(OSError, match="no permission"):
                runner.run(["sleep", "5"])
        finally:
            # monkeypatch.setattr(swt.os, ...) patches the real os module, so
            # subprocess.Popen.kill()'s own internal os.kill() call would hit
            # our fake_kill too unless we restore the real one first.
            monkeypatch.undo()
            # self.proc is set before _wait_testcomplete() ever runs, so it's
            # reachable here even though run() raised instead of returning.
            runner.proc.kill()
            runner.proc.wait()
