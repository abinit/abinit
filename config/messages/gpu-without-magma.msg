GPU support has been enabled but the MAGMA library is not
used. Calculations using the GPU will thus show poor
performance.

We highly recommend you to build Abinit with MAGMA support,
by including \\"magma\\" in the linear algebra flavor, e.g.
--with-linalg-flavor=\\"netlib+magma\\".

The minimal required version for MAGMA is 1.1.0.

For more details about the MAGMA library, please consult
http://icl.cs.utk.edu/magma/.

Please note that we do not provide support for the
installation of MAGMA.
