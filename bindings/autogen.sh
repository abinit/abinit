#!/bin/sh
#
# Copyright (C) 2014-2020 ABINIT group (Yann Pouillon)
#
# This file is part of the Abinit Bindings software package. For license
# information, please see the COPYING file in the top-level directory of the
# source distribution.
#

# Stop at first error encountered
set -e

# Check that we are in the right directory
if test ! -s "./configure.ac" -o ! -s "config/specs/bindings.conf"; then
  echo "[bndbuild]   This is not an Abinit Bindings source tree - aborting now" >&2
  exit 1
fi

# Create possibly missing directories
mkdir -p config/gnu

# Generate M4 macros
echo "[bndbuild]   Generating M4 macros"
./config/scripts/make-macros-bindings

# Generate makefiles
echo "[bndbuild]   Generating makefiles"
./config/scripts/make-makefiles-bindings

# Generate M4 includes
echo "[bndbuild]   Generating aclocal.m4"
aclocal -I config/m4

# Generate configure auxiliary files
#echo "[bndbuild]   Generating config.h.in"
#autoheader

# Generate configure
echo "[bndbuild]   Generating configure script"
autoconf

# Generate libtool scripts
#echo "[bndbuild]   Generating libtool scripts"
#my_libtoolize="libtoolize"
#${my_libtoolize} --version >/dev/null 2>&1
#if test "${?}" != "0"; then 
#  my_libtoolize="glibtoolize"
#fi
#${my_libtoolize} --version >/dev/null 2>&1
#if test "${?}" != "0"; then 
#  echo "[bndbuild]   Error: could not find a working version of libtoolize" >&2
#  exit 1
#fi
#${my_libtoolize} --automake --copy --force

# Generate makefile inputs
# Do not use "automake --force-missing", as it overwrites the INSTALL file.
echo "[bndbuild]   Generating Makefile.in for each directory"
automake --add-missing --copy
