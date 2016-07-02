Your compiler will work properly only if the LD_LIBRARY_PATH
and INTEL_LICENSE_FILE environment variables are correctly
set-up. Please make sure that this is the case before
building ABINIT.

In most cases it means sourcing a script in the Intel
compiler install directory, e.g.:

    . /opt/intel/fc/9.1/bin/ifortvars.sh

for a Bourne-like shell \\(do not forget the leading dot\\), or

    source /opt/intel/fc/9.1/bin/ifortvars.csh

for a C-like shell.

If you do not know what this is all about, please refer to
the Intel documentation or ask your system administrator.

Please note that the same holds for the Intel C and C++
compilers.
