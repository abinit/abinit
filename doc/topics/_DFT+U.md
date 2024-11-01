---
description: How to perform a DFT+U calculation
authors: BAmadon
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform a DFT+U calculation with the ABINIT package.

## Introduction

This feature is available only in PAW. The DFT+U framework is described in
[[cite:Anisimov1991]] and [[cite:Liechtenstein1995]]. In ABINIT, the DFT+U
approximation is implemented inside the PAW atomic spheres only. Two choices
of double counting are provided: the Full Localized limit and the Around Mean
Field approximation. Our implementation is described in [[cite:Amadon2008a]].
It follows the main lines of [[cite:Bengone2000]]. See also
[[cite:Czyzyk1994]]. Forces and stress are implemented. For details on
keywords ([[lpawu]], [[upawu]], [[jpawu]], [[usedmatpu]], [[dmatpuopt]],
[[dmatudiag]]) see keyword [[usepawu]] in input variables.

In both the output and log files, we can find:

- The DFT+U contribution of energy which is contained inside the PAW
  Spherical terms in the output file.

- The Decomposition of the DFT+U energy is given (Interaction energy, Double
  counting term, and sum of the two) in the log file.

- The orbital density matrix ($n_{m,m'}^{\sigma}$), also called occupation
matrix (corresponding to Eq.(9) of [[cite:Bengone2000]] and Eq.(1) of
[[cite:Liechtenstein1995]], see also [[cite:Amadon2008a]] and variable
[[dmatpuopt]]) is also given for each atom in the basis of real spherical
harmonics. It is given at each SCF step in the log file: one can thus check
the convergency of the calculation.

Consistency between total energy and forces in DFT+U have been checked.

The implementation of DFT+U in ABINIT allows also to impose a starting density
matrix in order to compare the energy of various electronic configuration (see
keywords [[usedmatpu]] and [[dmatpawu]]).



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:dftu|The tutorial on DFT+U]] shows how to perform a DFT+U calculation using ABINIT, and will lead to compute the projected DOS of NiO. Prerequisite: PAW1.

