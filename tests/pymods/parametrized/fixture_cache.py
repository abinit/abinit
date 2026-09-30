"""
Process-safe cache of fixture builds (plan 8.4).

A fixture is built at most once per cache directory and key, even when several processes or threads
(runtests.py -j N) need it at the same time:

- an exclusive flock on <root>/<entry>.lock serializes the builders of the same entry
  (flock, not lockf: lockf locks are per process and would not exclude threads);
- the build runs in a private directory <root>/<entry>.tmp-<pid>-<thread>, and status.json is written
  before the directory is renamed atomically to <root>/<entry>, so a half-built entry is never visible;
- a failed build is recorded as such and is not retried with the same key.

The lock is released by the kernel if the holder dies. Only the standard library is used:
this module is imported by the pickled test objects (see paramtests.py).
"""
from __future__ import annotations

import fcntl
import json
import os
import shutil
import stat
import threading
import time
from typing import Any, Callable

STATUS_FILE = "status.json"


def entry_name(name: str, key: str) -> str:
    return f"{name}-{key[:16]}"


def read_status(path: str) -> dict[str, Any] | None:
    try:
        with open(os.path.join(path, STATUS_FILE)) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def make_read_only(path: str) -> None:
    """Remove the write permission of the files in `path`, so that a user of the fixture cannot modify them."""
    for dirpath, _dirnames, filenames in os.walk(path):
        for f in filenames:
            p = os.path.join(dirpath, f)
            if os.path.islink(p):
                continue
            mode = os.stat(p).st_mode
            os.chmod(p, mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))


def ensure(root: str, name: str, key: str,
           build: Callable[[str], tuple]) -> tuple[str, dict[str, Any]]:
    """
    Return (directory of the entry, status) for the fixture `name` with cache key `key`, building it if needed.

    build(tmpdir) runs the fixture inside tmpdir and returns (status, message) or (status, message, extra),
    where status is "ok", "failed" or "skipped" and extra is a dict stored in status.json.
    It is called at most once per (root, name, key) among all the processes sharing `root`.
    """
    os.makedirs(root, exist_ok=True)
    entry = entry_name(name, key)
    final = os.path.join(root, entry)

    status = read_status(final)
    if status is not None:
        return final, status

    with open(os.path.join(root, entry + ".lock"), "a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            status = read_status(final)
            if status is not None:              # built by another process while we waited
                return final, status

            tmp = os.path.join(root, f"{entry}.tmp-{os.getpid()}-{threading.get_ident()}")
            shutil.rmtree(tmp, ignore_errors=True)
            os.makedirs(tmp)
            start = time.time()
            extra: dict[str, Any] = {}
            try:
                state, message, *rest = build(tmp)
                extra = rest[0] if rest else {}
            except Exception as exc:            # a bug in the build must not leave the entry undefined
                state, message = "failed", f"{type(exc).__name__}: {exc}"
            status = {"name": name, "key": key, "status": state, "message": message,
                      "build_time": round(time.time() - start, 3), "pid": os.getpid(), **extra}
            with open(os.path.join(tmp, STATUS_FILE), "w") as fh:
                json.dump(status, fh, indent=2)
            if state == "ok":
                make_read_only(tmp)
            if os.path.exists(final):           # leftover of a crash after the rename: never trusted
                shutil.rmtree(final, ignore_errors=True)
            os.rename(tmp, final)
            return final, status
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def remove_all(root: str) -> None:
    """Remove the whole cache (run_params.py --rebuild-fixtures). Entries are read-only: restore the permission."""
    if not os.path.isdir(root):
        return
    for dirpath, dirnames, filenames in os.walk(root):
        os.chmod(dirpath, 0o755)
        for f in filenames:
            p = os.path.join(dirpath, f)
            if not os.path.islink(p):
                os.chmod(p, 0o644)
    shutil.rmtree(root)


__all__ = ["STATUS_FILE", "ensure", "entry_name", "make_read_only", "read_status", "remove_all"]
