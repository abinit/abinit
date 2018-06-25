---
description: How to set parameters for a PAW calculation
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to set parameters for a PAW calculation with the ABINIT package.

## Introduction

The PAW atomic data can be used with plane waves as well as with wavelets.
Specificities of PAW for use with planewaves are presented here. See
[[topic:Wavelets]] for its use with wavelets.

The way the PAW method is implemented with planewaves in ABINIT is described
in [[cite:Torrent2008]].  
The use of PAW atomic data (equivalent to pseudopotential file for the norm-
conserving case) automatically launch a PAW calculation. ABINIT is provided
with the JTH [[cite:Jollet2014]] PAW atomic data table on the ABINIT web site.  
To perform a standard PAW calculation, the input file is the same than for a
norm-conserving one, except that the variable [[pawecutdg]] must be specified
(see below). In the case the input variable [[accuracy]] is used, the input
variable [[pawecutdg]] is automatically used.  
Some physical functionalities are available only in the PAW framework: DFT+U,
DMFT, local exact exchange,...



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* The tutorial on the use of PAW ([[tutorial:paw1|PAW1]]) presents the Projector-Augmented Wave method, implemented in ABINIT as an alternative to norm-conserving pseudopotentials, with a sizeable accuracy and CPU time advantage.
* The tutorial on the generation of PAW atomic data files ([[tutorial:paw2|PAW2]]) presents the generation of atomic data for use with the PAW method. Prerequisite: PAW1.
* The tutorial on the validation of a PAW atomic datafile ([[tutorial:paw3|PAW3]]) demonstrates how to test a generated PAW dataset using ABINIT, against the ELK all-electron code, for diamond and magnesium. Prerequisite: PAW1 and PAW2.

