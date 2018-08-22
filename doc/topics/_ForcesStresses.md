---
description: How to to tune the computation of forces and stresses
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to to tune the computation of forces and stresses with the ABINIT package.

## Introduction

Hellman-Feynman forces are computed from an analytical formula, and
corresponds exactly to the limit of finite differences of energy for
infinitesimally small atomic displacements when the ground-state calculation
is at convergence. This feature is available for all the cases where the total
energy can be computed. A correction for non-converged cases allows to get
accurate forces with less converged wavefunctions than without it. The
decomposition of the forces in their different components can be provided.

Stress can also be computed. This feature is available for all the cases where
the total energy can be computed (except wavelets). The decomposition of the
stresses in their different components can be provided. A smearing scheme
applied to the kinetic energy [[ecutsm]] allows one to get smooth energy
curves as a function of lattice parameters and angles. A target stress can be
given by the user ([[strtarget]]), the geometry optimization algorithm will
try to find the primitive cell and atomic positions that deliver that target stress.

The computation of forces and stresses is optional, see [[optforces]] and
[[optstress]]. They are used to define SCF stopping criteria ([[toldff]],
[[tolrff]]) or geometry optimization stopping criteria ([[tolmxf]]). For the
geometry optimization, combined cell shape and atomic position optimization
need a conversion scale, set by [[strprecon]].



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:base1|The tutorial 1]] deals with the H2 molecule: get the total energy, the electronic energies, the charge density, the bond length, the atomisation energy 
* [[tutorial:base2|The tutorial 2]] deals again with the H2 molecule: convergence studies, LDA versus GGA 
* [[tutorial:base3|The tutorial 3]] deals with crystalline silicon (an insulator): the definition of a k-point grid, the smearing of the cut-off energy, the computation of a band structure, and again, convergence studies ...
* [[tutorial:base4|The tutorial 4]] deals with crystalline aluminum (a metal), and its surface: occupation numbers, smearing the Fermi-Dirac distribution, the surface energy, and again, convergence studies ...

