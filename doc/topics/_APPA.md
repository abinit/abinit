---
description: How to use the APPA post-processing tool for the analysis of molecular dynamics output files (trajectories)
authors: SS, XG, YG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to use the APPA post-processing tool for the analysis of molecular dynamics
output files (trajectories) with the ABINIT package.

## Introduction

APPA is a graphical (or text) post-processor for ABINIT dedicated to the
analysis of Molecular Dynamics output files (trajectories).

Output files from large scale molecular dynamics simulations are not easy to
manage. All the data can be extracted from output files in NetCDF or
ACSCII(text) formats. APPA is a graphical tool written in python (PyQt4
library) able to read ABINIT output, easy to use in order to follow a
trajectory from molecular dynamics simulation (total energy, pressure...). It
is also possible to calculate velocity autocorrelation function from the
simulation, radial pair distribution, vibrational density of states, etc....
If Python is able to handle efficiently these kinds of output files, this one
is not suited to compute some quantities coming from large-scale MD
simulation. That?s why APPA also uses Fortran, and thanks to the F2PY library,
the connection between these both languages is possible.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

