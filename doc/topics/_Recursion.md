---
description: How to perform orbital-free calculations
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform orbital-free calculations with the ABINIT package.

## Introduction

It is possible to use Thomas-Fermi kinetic functional (explicit functional of
the density) or Thomas-Fermi-Weizsacker kinetic functional (with Gradient
Corrections) instead of Kohn-Sham kinetic energy functional (implicit
functional of the density through Kohn-Sham wavefunctions).  
See [[cite:Perrot1979]].  
The Recursion Method may be used in order to compute electronic density,
entropy, Fermi energy and eigenvalues energy. This method computes the density
without computing any orbital, is efficient at high temperature, with a
efficient parallelization (almost perfect scalability).

At present, it only works for local pseudopotentials, severely restricting the
use of this method.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

