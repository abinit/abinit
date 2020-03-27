#!/bin/bash
#
# Copyright (C) 2011-2020 ABINIT Group (Jean-Michel Beuken)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

sed -i -e 's/\t\$.FCFLAGS. \\//' src/98_main/Makefile
