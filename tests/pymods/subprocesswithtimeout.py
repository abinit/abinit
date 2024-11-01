from __future__ import print_function, division, absolute_import #, unicode_literals

import time
import os
import errno
import signal
import subprocess


class TimeoutError(Exception):
    """Exceptions raised by SubProcessWithTimeout."""


class SubProcessWithTimeout(object):
    """
    Based on from http://stackoverflow.com/questions/3876886/timeout-a-subprocess?rq=1
    """
    Error = TimeoutError

    def __init__(self, timeout, delay=.05):
        self.timeout = float(timeout)
        self.delay = float(delay)

        if self.delay > self.timeout or self.delay <= 0 or self.timeout <= 0:
            raise ValueError("delay and timeout must be positive with delay <= timeout")

    def run(self, args,
            bufsize=0, executable=None, stdin=None, stdout=None, stderr=None, preexec_fn=None,
            close_fds=False, shell=False, cwd=None, env=None, universal_newlines=False, startupinfo=None, creationflags=0):
        """Same interface as Popen"""

        self.proc = subprocess.Popen(args,
                                     bufsize, executable, stdin, stdout, stderr, preexec_fn,
                                     close_fds, shell, cwd, env, universal_newlines, startupinfo, creationflags)

        return_code = self._wait_testcomplete()
        return self.proc, return_code

    def _wait_testcomplete(self):
        start = time.time()
        while (time.time()-start) < self.timeout:
            if self.proc.poll() is not None:  # 0 just means successful exit
                return self.proc.returncode
            else:
                time.sleep(self.delay)
        # The process may exit between the time we check and the
        # time we send the signal.
        try:
            os.kill(-self.proc.pid, signal.SIGTERM)
        except OSError as e:
            # If it's not because the process no longer exists, something weird is wrong.
            if e.errno != errno.ESRCH:
                raise e
        time.sleep(1)

        if self.proc.poll() is None:  # Still hasn't exited.
            try:
                os.kill(-self.proc.pid, signal.SIGKILL)
            except OSError as e:
                if e.errno != errno.ESRCH:
                    raise e
            return 137  # timeout return code for SIGKILL
        else:
            return 124  # timeout return code for SIGTERM


#############################################################################################################
# Unit tests
import unittest


class TestSubProcessWithTimeout(unittest.TestCase):
    def test_with_sleep(self):
        """"Testing if sleep 5 raises TimeoutError"""
        proc, retcode = SubProcessWithTimeout(1).run(["sleep", "5"])
        self.assertEqual(retcode, 124)


if __name__ == "__main__":
    unittest.main()
