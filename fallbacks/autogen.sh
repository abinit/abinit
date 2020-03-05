#!/bin/sh
#
# Copyright (C) 2011-2020 ABINIT group (Yann Pouillon)
#
# This file is part of the Abinit Test Suite software package. For license
# information, please see the COPYING file in the top-level directory of the
# source distribution.
#

# Stop at first error encountered
set -e

# Check that we are in the right directory
if test ! -s "./build-abinit-fallbacks.sh.in"; then
  echo "[fbkbuild]   This is not an Abinit fallbacks source tree - aborting now" >&2
  exit 1
fi

