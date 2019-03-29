You have decided to use the Abinit version of linear
algebra libraries \\(BLAS, LAPACK\\), which is a random mix of
routines coming from official Netlib versions ranging from
3.0 \\(released June 30, 1999\\) to 3.2 \\(released November 18,
2008\\). They were downloaded from 'http://www.netlib.org/'.

The current version of Netlib\\'s libraries is 3.4.0 \\(released
November 11, 2011\\). Since version 3.2, more than 60
significant bugs - including memory leaks - have been fixed,
various performance improvements have been made available,
and routine interfaces have been refactored.

Moreover, the Abinit version of the libraries rely on a
heavily modified build system and will be compiled with
build parameters optimized for Abinit, not for BLAS or
LAPACK. This might prove somewhat suboptimal or inadequate
on some computer architectures. Please note as well that the
Abinit test suite does not provide specific tests for linear
algebra.

By using the libraries provided by Abinit, you thus accept
to experience poor calculation performance, unlikely but 
possible random crashes without error messages, and possibly 
inaccurate results. However, if this is not what you 
intended originally, you may instead select one of the 
various optimized linear algebra implementations supported 
by Abinit.

For production-grade calculations, you will find detailed
documentation on how to configure Abinit for optimized
linear algebra libraries in the \\"Linear algebra
connector\\" section of Abinit\\'s config file template:
'~abinit/doc/build/config-template.ac8'.

Do not hesitate to contact your system administrator in
case of doubt or if you need linear algebra libraries
installed into your local high-performance computing
environment.
