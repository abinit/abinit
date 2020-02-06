#!/bin/sh
#
## Copyright (C) 2017-2019 Yann Pouillon <devops@materialsevolution.es>

# Stop at first error encountered
set -e

# Check that we are in the right directory
if test ! -s "./configure.ac" -o ! -s "src/m_paw_atom.F90"; then
  echo "This is not a LibPAW source tree - aborting now"
  exit 1
fi

# Make sure the directory tree is writable
chmod -R u+w .

# Create possibly missing directories
mkdir -p config/gnu config/m4

# Generate makefiles for Automake
# FIXME: lists of source files are not automatically updated
echo "[pawbuild]   Generating makefiles for Automake (FIXME: update source lists)"
cp config/makefiles/doc.am doc/Makefile.am
cp config/makefiles/libxc.am libxc/Makefile.am
cp config/makefiles/src.am src/Makefile.am
cp config/makefiles/top.am Makefile.am

# Generate libtool scripts
echo "[pawbuild]   Generating libtool scripts"
my_libtoolize="libtoolize"
${my_libtoolize} --version >/dev/null 2>&1
if test "${?}" != "0"; then 
  my_libtoolize="glibtoolize"
fi
${my_libtoolize} --version >/dev/null 2>&1
if test "${?}" != "0"; then 
  echo "Error: could not find a working version of libtoolize" >&2
  exit 1
fi
${my_libtoolize} --automake --copy --force

# Generate M4 includes
echo "[pawbuild]   Generating aclocal.m4"
aclocal -I config/m4

# Generate configure auxiliary files
echo "[pawbuild]   Generating config.h.in"
autoheader

# Generate configure
echo "[pawbuild]   Generating configure script"
autoconf

# Generate makefile inputs
# Do not use "automake --force-missing", as it overwrites the INSTALL file.
echo "[pawbuild]   Generating Makefile.in for each directory"
automake --add-missing --copy
