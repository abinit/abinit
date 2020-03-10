#!/bin/sh
#
# Copyright (C) 2010-2014 ABINIT group (Yann Pouillon)
#
# This file is part of the Abinit software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

# Stop at first error encountered
set -e

# Check that we are in the right directory
if test ! -s "./configure.ac" -o ! -s "config/specs/fallbacks.conf"; then
  echo "[fbkbuild]   This is not an Abinit Fallbacks source tree - aborting now" >&2
  exit 1
fi

# Create possibly missing directories
mkdir -p config/gnu config/m4

# Generate input files for configure
echo "[fbkbuild]   Generating configure inputs"
./config/scripts/make-config-dumper

# Generate M4 macros
echo "[fbkbuild]   Generating M4 macros"
./config/scripts/make-macros-environment
./config/scripts/make-macros-fallbacks
./config/scripts/make-macros-info
./config/scripts/make-macros-options
./config/scripts/make-macros-patches

# Generate source files
echo "[fbkbuild]   Generating source files"
./config/scripts/make-fallbacks-config-in
./config/scripts/make-install-symlinks-in

# Generate makefiles
echo "[fbkbuild]   Generating makefiles"
./config/scripts/make-makefiles-fallbacks

# Generate M4 includes
echo "[fbkbuild]   Generating aclocal.m4"
aclocal -I config/m4

# Generate configure auxiliary files
#echo "[fbkbuild]   Generating config.h.in"
#autoheader

# Generate configure
echo "[fbkbuild]   Generating configure script"
autoconf

# Generate libtool scripts
#echo "[fbkbuild]   Generating libtool scripts"
#my_libtoolize="libtoolize"
#${my_libtoolize} --version >/dev/null 2>&1
#if test "${?}" != "0"; then 
#  my_libtoolize="glibtoolize"
#fi
#${my_libtoolize} --version >/dev/null 2>&1
#if test "${?}" != "0"; then 
#  echo "[fbkbuild]   Error: could not find a working version of libtoolize" >&2
#  exit 1
#fi
#${my_libtoolize} --automake --copy --force

# Generate makefile inputs
# Do not use "automake --force-missing", as it overwrites the INSTALL file.
echo "[fbkbuild]   Generating Makefile.in for each directory"
automake --add-missing --copy
