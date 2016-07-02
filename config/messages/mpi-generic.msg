You are using a generic MPI installation. In this case, for
consistency reasons, both sequential and parallel binaries
have to be built with the same compiler - the MPI-capable
one - which means that the sequential binaries will depend
on the MPI libraries, even though they do not make use of
them.

This, in turn, means that you will not be able to run ABINIT
binaries unless you properly set the LD_LIBRARY_PATH
environment variable. The configure script just did it now,
in order to perform its task seamlessly. Please read the
\\"Environment variables\\" section in the chapter dedicated to
the configure script of the ABINIT Build System Guide that
you will find in the ~abinit/doc/build/ subdirectory of the
source tree to make the proper settings more permanent.

Another solution might be to pass the \\"-static\\" flag
to the C, C++ and Fortran linkers respectively through the
\\"--with-cc-ld-optflags\\", \\"--with-cxx-ld-optflags\\" and
\\"--with-fc-ld-optflags\\" options of the configure script, but
this cannot be done on all platforms.
