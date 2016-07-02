The build system of Abinit has not been able to download
one or more required external packages. An unstable or
absent connection to internet of this computer is the most
probable culprit. This means that the build cannot
proceed, and the configure script will thus abort now with
an error.

To address this issue, you may:

  \\* use pre-compiled versions of the external packages,
    through the --with-\\*-bins, --with-\\*-incs, and
    --with-\\*-libs options of configure\\; please see
    the section associated to connectors in
    \\~abinit/doc/build/config-template.ac for available
    options and details\\;

  \\* manually upload the source tarballs to a dedicated
    directory, either \\~/.abinit/tarballs/ \\(the default\\) or
    any directory of your choice\\; please see the section
    associated to fallbacks in
    \\~abinit/doc/build/config-template.ac for details\\;

  \\* disable the corresponding extensions, by using the
    following options of configure:
      --with-dft-flavor=\\"none\\"\\;
      --with-trio-flavor=\\"none\\"\\;
    please see the section associated to connectors in
    \\~abinit/doc/build/config-template.ac for details.

