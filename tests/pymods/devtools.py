"""
General-purpose developer tools and utilities for the ABINIT test suite.
Includes system CPU/GPU detection, cross-platform file locking, and function decorators.
"""
from __future__ import annotations

import errno
import os
import shutil
import subprocess
import time
from functools import wraps


def number_of_cpus() -> int:
    """
    Detect the number of physical or virtual CPUs on the system.

    Returns:
        int: Number of CPUs detected, or -1 if detection fails.
    """
    import os
    import re
    import subprocess

    # Python 2.6+
    #try:
    #    import multiprocessing
    #    return multiprocessing.cpu_count()
    #except (ImportError, NotImplementedError):
    #    pass

    # POSIX
    try:
        res = int(os.sysconf("SC_NPROCESSORS_ONLN"))
        if res > 0: return res
    except (AttributeError, ValueError):
        pass

    # Windows
    try:
        res = int(os.environ["NUMBER_OF_PROCESSORS"])
        if res > 0: return res
    except (KeyError, ValueError):
        pass

    # jython
    try:
        from java.lang import Runtime
        runtime = Runtime.getRuntime()
        res = runtime.availableProcessors()
        if res > 0: return res
    except ImportError:
        pass

    # BSD
    try:
        sysctl = subprocess.Popen(["sysctl", "-n", "hw.ncpu"], stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)
        if res > 0: return res
    except (OSError, ValueError):
        pass

    # Linux
    try:
        res = open("/proc/cpuinfo").read().count("processor\t:")
        if res > 0: return res
    except OSError:
        pass

    # Solaris
    try:
        pseudoDevices = os.listdir("/devices/pseudo/")
        expr = re.compile("^cpuid@[0-9]+$")
        res = 0
        for pd in pseudoDevices:
            if expr.match(pd) is not None:
                res += 1
        if res > 0: return res
    except OSError:
        pass

    # Other UNIXes (heuristic)
    try:
        try:
            dmesg = open("/var/run/dmesg.boot").read()
        except OSError:
            dmesgProcess = subprocess.Popen(["dmesg"], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]

        res = 0
        while "\ncpu" + str(res) + ":" in dmesg:
            res += 1

        if res > 0: return res
    except OSError:
        pass

    return -1
    #raise Exception('Cannot determine number of CPUs on this system')

def number_of_gpus() -> int:
    """
    Detect the number of GPUs using vendor-specific tools (`nvidia-smi` or `roc-smi`).

    Returns:
        int: Number of GPUs detected, or 0 if none are available.
    """
    # Look for NVIDIA GPU first, then AMD GPU...
    nvidia_cmd = ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"]
    amdgpu_cmd = ["roc-smi", "--listgpu"]

    num_gpus = 0
    for gpu_cmd in [nvidia_cmd, amdgpu_cmd]:
        if shutil.which(gpu_cmd[0]) is None:
            continue

        try:
            result = subprocess.run(gpu_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
            # The text argument was introduced in Python 3.7 as an alias for universal_newlines=True.
            #result = subprocess.run(gpu_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            # Check if command failed (meaning it exists)
            if result.returncode != 0:
                print(f"Error while executing {gpu_cmd[1]}:\n{result.stderr}")
                num_gpus = 0

            # Command was successful, count the lines (one per GPU) and exit
            gpu_lines = result.stdout.strip().split("\n")
            num_gpus = len(gpu_lines)
            break

        except FileNotFoundError:
            # Command doesn't exist, continue the loop and try another
            continue

    return num_gpus


class FileLockException(Exception):
    """Exception raised by FileLock."""


class FileLock:
    """
    A cross-platform file locking mechanism with context manager support.

    This class implements a simple advisory lock by creating a '.lock' file.
    It supports use as a context manager for easy acquisition and release.
    Wait times and delays can be configured to handle lock contention.
    """
    # Create an alias for compatibility
    Error = FileLockException

    def __init__(self, file_name, timeout=10, delay=.05):
        """
        Initialize the file lock.

        Args:
            file_name: Name of the file to lock.
            timeout: Maximum time (in seconds) to wait for the lock.
            delay: Delay (in seconds) between successive lock attempts.
        """
        self.file_name = file_name
        self.lockfile = os.path.abspath(file_name) + ".lock"
        self.timeout = float(timeout)
        self.delay = float(delay)
        self.is_locked = False

        if (self.delay > self.timeout or
            self.delay   <= 0 or
            self.timeout <= 0):
            err_msg = "delay and timeout must be positive with delay <= timeout"
            raise ValueError(err_msg)

    @classmethod
    def FakeLock(cls, file_name: str, timeout: float = 10, delay: float = .05) -> FileLock:
        """
        Create a lock object that does nothing (monkey-patched acquire/release).

        Args:
            file_name: Path to the target file.
            timeout: Timeout for the lock attempt.
            delay: Interval between attempts.

        Returns:
            FileLock: A fake lock instance.
        """
        fake = cls(file_name, timeout=timeout, delay=delay)

        def nop():
            """Monkey patch."""

        fake.acquire = nop
        fake.release = nop
        return fake

    def acquire(self):
        """
        Acquire the lock.

        Retries every `delay` seconds until the lock is acquired or `timeout`
        is reached.

        Raises:
            FileLockException: If the lock cannot be acquired within the timeout.
        """
        start_time = time.time()
        while True:
            try:
                self.fd = os.open(self.lockfile, os.O_CREAT|os.O_EXCL|os.O_RDWR)
                break
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
                if (time.time() - start_time) >= self.timeout:
                    raise FileLockException("Timeout occurred.")
                time.sleep(self.delay)

        self.is_locked = True

    def release(self):
        """
        Release the lock by deleting the lock file.

        This is called automatically when using the context manager.
        """
        if self.is_locked:
            os.close(self.fd)
            os.unlink(self.lockfile)
            self.is_locked = False

    def __enter__(self):
        """
        Enter the runtime context related to this object.

        Automatically acquires the lock.

        Returns:
            FileLock: The locked instance.
        """
        if not self.is_locked: self.acquire()
        return self

    def __exit__(self, type, value, traceback):
        """
        Exit the runtime context related to this object.

        Automatically releases the lock.

        Args:
            type: Exception type.
            value: Exception value.
            traceback: Exception traceback.
        """
        if self.is_locked: self.release()

    def __del__(self):
        """
        Destructor to ensure the lock file is released when the instance is deleted.
        """
        self.release()


class NoErrorFileLock(FileLock):
    """
    A file locker that suppresses `FileLockException` during context entry.

    Returns True if the lock was acquired, False otherwise.
    """

    def __enter__(self):
        """
        Enter the runtime context and attempt to acquire the lock.

        Returns:
            bool: True if the lock was successfully acquired, False otherwise.
        """
        try:
            self.acquire()
        except self.Error:
            return False
        else:
            return True


def makeunique(gen):
    """
    Decorator that ensures a generator produces unique outputs by caching them.

    This is useful for generators that might produce duplicate items (e.g., random
    name generators) when unique values are required.

    Args:
        gen (callable): The generator function to wrap.

    Returns:
        callable: A wrapped generator that filters out duplicate values.
    """
    cache = set()

    @wraps(gen)
    def generator(*args):
        s = gen(*args)
        while s in cache:
            s = gen(*args)
        cache.add(s)
        return s

    return generator
