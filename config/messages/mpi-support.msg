MPI support is disabled by default because it makes the
configure script abort when the build environment is not
properly set-up.

You may activate it using the --enable-mpi=\\"yes\\"
option, provided that the PATH and LD_LIBRARY_PATH
environment variables are properly set-up. If this is not
the case, the configure script will abort with the following
error message:

    configure: error: cannot run C compiled programs.
    If you meant to cross compile, use \\`--host\\'.
    See \\`config.log\\' for more details.
