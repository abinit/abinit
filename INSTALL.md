# Installation instructions specific to ABINIT

Please report problems or comments on the ABINIT [forum](https://forum.abinit.org).

## Quick install, starting from the tar.gz file:

Download, gunzip and untar the file abinit-x.y.z.tar.gz with:

    tar -zxvf abinit-x.y.z.tar.gz

where x.y.z is to be replaced by the actual numbers that characterize the version of ABINIT.

cd abinit-x.y.z and run the configuration script with:

./configure --prefix=~/local/

where `prefix` specifies the directory in which binaries and documentation will be installed.

Use:

    ./make -j4

to compile the code with e.g. four thread
(set this number according the number of CPUs available on your machine.

Optionally, install the package inside the `prefix` directory with:

./make install

For troubleshooting, more elaborate installation instructions are contained
in the [doc/install\_notes](doc/install_notes) directory of this package
(take the html file that corresponds to the version that you are using), or
can also be found on-line at <https://www.abinit.org/installation> for the latest production version.

# GNU AUTOTOOLS SUPPORT

This version of ABINIT includes full support for the GNU Autotools. If you
make changes to the source code and/or the documentation, you may need
recent test versions of GNU build tools to regenerate the intermediate
files (e.g. Makefile.am, Makefile.in, config.h.in).

The following versions were used to generate the intermediate files in
this distribution:

 * GNU Autoconf 2.69 <ftp://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz>

 * GNU Automake 1.15 <ftp://ftp.gnu.org/gnu/automake/automake-1.15.tar.gz>

 * GNU Libtool 2.2.4 <ftp://ftp.gnu.org/gnu/libtool/libtool-2.2.4.tar.gz>

 * GNU M4 1.4 <ftp://ftp.gnu.org/gnu/m4/m4-1.4.4.tar.gz>

 * DocBook-utils 0.6.13 <http://sources.redhat.com/docbook-tools/>

Any question regarding the GNU Autotools support inside ABINIT should be
directed to Yann Pouillon <yann.pouillon@ehu.es>.
