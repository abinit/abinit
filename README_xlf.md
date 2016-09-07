# Building Abinit 6 with the IBM XL Fortran compiler

## Troubleshooting

If the compiler complains about invalid syntax related to `config.h`,
the most likely origin of the problem is an incorrect default setting of the
preprocessor options in the config file of xlf, explicitely asking to
leave C comments in the preprocessed file.

This is solved by editing `xlf.cfg`, replacing:

    cppoptions = -C

by:

    cppoptions = -P

If you do not have write permissions on this file, then copy it into
your home directory before modifying it. It is usually located in `/etc`
or `/etc/opt/ibmcmp/xlf/<version>`.

You will also need either to reconfigure Abinit the following way:

    ../configure FC="xlf -F /path/to/xlf.cfg"

or to force `make` to use the modified compiler:

    make FC="xlf -F /path/to/xlf.cfg"

Of course, do not hesitate to replace `xlf` by any other name if needed,
e.g. `mpif90`.

If `make` succeed, you should probably consider asking your system
administrator to edit the system-wide config file or create another one.
