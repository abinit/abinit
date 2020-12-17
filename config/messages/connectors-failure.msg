The build system of Abinit has not been able to detect
one or more requested external packages, probably because
you did not provide enough information. This means
that the build cannot proceed, and the configure script
will thus abort with an error.

To address this issue, you may:

  \\* make sure that the prerequisites are available and enabled, e.g.
    --with-dft-flavor=\\"libxc+bigdft\\"\\ or
    --with-trio-flavor=\\"netcdf\\";

  \\* use the --with-\\*-bins, --with-\\*-incs, and
    --with-\\*-libs options of configure to provide
    additional information\\; please see the section related
    to connectors in
    \\~abinit/doc/build/config-template.ac9 for available
    options and details\\;

  \\* disable the corresponding extensions, by using the
    associated --with-\\*-flavor=\\"none\\" options of
    configure\\; please see the section associated to
    connectors in \\~abinit/doc/build/config-template.ac9 for
    details.

