#!/bin/sh
#
# Copyright (C) 2011-2024 ABINIT group (Yann Pouillon)
#
# This file is part of the Abinit Documentation software package. For license
# information, please see the COPYING file in the top-level directory of the
# source distribution.
#

# Stop at first error encountered
set -e

# Check that we are in the right directory
if test ! -s "./configure.ac" -o ! -s "config/specs/documents.conf"; then
  echo "[docbuild]   This is not an Abinit Documentation source tree - aborting now" >&2
  exit 1
fi

# Create possibly missing directories
mkdir -p config/gnu

# Generate M4 macros
echo "[docbuild]   Generating M4 macros"
./config/scripts/make-macros-doc

# Generate makefiles
echo "[docbuild]   Generating makefiles"
./config/scripts/make-makefiles-doc

# Generate M4 includes
echo "[docbuild]   Generating aclocal.m4"
aclocal -I config/m4

# Generate configure auxiliary files
#echo "[docbuild]   Generating config.h.in"
#autoheader

# Generate configure
echo "[docbuild]   Generating configure script"
autoconf

# Generate libtool scripts
#echo "[docbuild]   Generating libtool scripts"
#my_libtoolize="libtoolize"
#${my_libtoolize} --version >/dev/null 2>&1
#if test "${?}" != "0"; then 
#  my_libtoolize="glibtoolize"
#fi
#${my_libtoolize} --version >/dev/null 2>&1
#if test "${?}" != "0"; then 
#  echo "[docbuild]   Error: could not find a working version of libtoolize" >&2
#  exit 1
#fi
#${my_libtoolize} --automake --copy --force

# Generate makefile inputs
# Do not use "automake --force-missing", as it overwrites the INSTALL file.
echo "[docbuild]   Generating Makefile.in for each directory"
automake --add-missing --copy
