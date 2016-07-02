"""MPI distribution of work load"""
import sys, traceback
from contextlib import contextmanager
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

i_am_master = bool(rank == 0)

@contextmanager
def mpi_abort_if_exception():
    """Terminate all mpi process if an exception is raised."""
    try:
        yield
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_traceback)
        MPI.COMM_WORLD.Abort(1)

def mpi_watch(f):
    """Decorator. Terminate all mpi process if an exception is raised."""
    def g(*args, **kwargs):
        with mpi_abort_if_exception():
            return f(*args, **kwargs)
    return g

def master_only(f):
    """Decorator. Let a function be executed only by master."""
    def g(*args, **kwargs):
        if i_am_master:
            with mpi_abort_if_exception():
                return f(*args, **kwargs)
        return
    return g

