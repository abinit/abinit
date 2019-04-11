from __future__ import print_function, division, absolute_import #, unicode_literals

import os
import time
import errno
from functools import wraps


def number_of_cpus():
    """
    Number of virtual or physical CPUs on this system, i.e.
    user/real as output by time(1) when called with an optimally scaling userspace-only program
    Return -1 if ncpus cannot be detected
    taken from:
    http://stackoverflow.com/questions/1006289/how-to-find-out-the-number-of-cpus-in-python
    """
    import os, re, subprocess

    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError, NotImplementedError):
        pass

    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))
        if res > 0: return res
    except (AttributeError, ValueError):
        pass

    # Windows
    try:
        res = int(os.environ['NUMBER_OF_PROCESSORS'])
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
        sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'], stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)
        if res > 0: return res
    except (OSError, ValueError):
        pass

    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')
        if res > 0: return res
    except IOError:
        pass

    # Solaris
    try:
        pseudoDevices = os.listdir('/devices/pseudo/')
        expr = re.compile('^cpuid@[0-9]+$')
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
            dmesg = open('/var/run/dmesg.boot').read()
        except IOError:
            dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]

        res = 0
        while '\ncpu' + str(res) + ':' in dmesg:
            res += 1

        if res > 0: return res
    except OSError:
        pass

    return -1
    #raise Exception('Cannot determine number of CPUs on this system')


class FileLockException(Exception):
    """Exception raised by FileLock."""


class FileLock(object):
    """ A file locking mechanism that has context-manager support so
        you can use it in a with statement. This should be relatively cross
        compatible as it doesn't rely on msvcrt or fcntl for the locking.
        Taken from http://www.evanfosmark.com/2009/01/cross-platform-file-locking-support-in-python/
    """
    Error = FileLockException

    def __init__(self, file_name, timeout=10, delay=.05):
        """ Prepare the file locker. Specify the file to lock and optionally
            the maximum timeout and the delay between each attempt to lock.
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
    def FakeLock(cls, file_name, timeout=10, delay=.05):
        """Returns a fake lock file."""
        fake = cls(file_name, timeout=timeout, delay=delay)

        def nop():
            """Monkey patch."""

        fake.acquire = nop
        fake.release = nop
        return fake

    def acquire(self):
        """ Acquire the lock, if possible. If the lock is in use, it check again
            every `wait` seconds. It does this until it either gets the lock or
            exceeds `timeout` number of seconds, in which case it throws
            an exception.
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
                    raise FileLockException("Timeout occured.")
                time.sleep(self.delay)

        self.is_locked = True

    def release(self):
        """ Get rid of the lock by deleting the lockfile.
            When working in a `with` statement, this gets automatically
            called at the end.
        """
        if self.is_locked:
            os.close(self.fd)
            os.unlink(self.lockfile)
            self.is_locked = False

    def __enter__(self):
        """ Activated when used in the with statement.
            Should automatically acquire a lock to be used in the with block.
        """
        if not self.is_locked: self.acquire()
        return self

    def __exit__(self, type, value, traceback):
        """ Activated at the end of the with statement.
            It automatically releases the lock if it isn't locked.
        """
        if self.is_locked: self.release()

    def __del__(self):
        """ Make sure that the FileLock instance doesn't leave a lockfile
            lying around.
        """
        self.release()


class NoErrorFileLock(FileLock):
    '''
    A file locker that never raise a FileLockErrorin call of __enter__ but
    return a boolean to tell wether the lock.
    '''

    def __enter__(self):
        try:
            self.acquire()
        except self.Error:
            return False
        else:
            return True


def makeunique(gen):
    '''
    gen have to be random enought not to produce too often the same thing
    '''
    cache = set()

    @wraps(gen)
    def generator(*args):
        s = gen(*args)
        while s in cache:
            s = gen(*args)
        cache.add(s)
        return s

    return generator
