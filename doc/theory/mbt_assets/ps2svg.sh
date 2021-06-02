#!/bin/bash

# you need gs-common, pstoedit and skencil to
# get this script working
export BASENAME="`basename $1`";

# Outline fonts
eps2eps -dNOCACHE ${BASENAME}.ps ${BASENAME}2.ps

# Fix bounding box
ps2epsi  ${BASENAME}2.ps  ${BASENAME}.ps
rm ${BASENAME}2.ps

# Convert to Sketch
pstoedit -f sk  ${BASENAME}.ps  ${BASENAME}.sk

# Convert to SVG
skconvert  ${BASENAME}.sk  ${BASENAME}.svg

