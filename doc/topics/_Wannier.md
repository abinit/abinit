---
description: How to perform Wannier functions calculation
authors: BAmadon
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform Wannier functions calculation with the ABINIT package.

## Introduction

There are two ways to obtain Wannier functions with ABINIT:

* The first one is to use an internal implementation of Projected Local Orbital Wannier functions [[cite:Amadon2008]], [[cite:Amadon2015]] that can be activated with the input variable [[plowan_compute]]. The variables [[plowan_bandi]] and [[plowan_bandf]] specifies the Kohn Sham bands that are used to built Wannier functions, whereas variables [[plowan_natom]],[[plowan_iatom]],[[plowan_nbl]], [[plowan_lcalc]], and [[plowan_projcalc]] specify the atoms, angular momentum and projectors corresponding to the projected local Orbital Wannier functions. 

In order to do Wannier interpolation or analysis of hopping integrals, real
space Wannier functions can be built using the variables [[plowan_realspace]]
and [[plowan_nt]] as well as [[plowan_it]]. Real space Hamiltonian for
analysis is given in latex file whose name ending is "wanniereigen" and the
band structure is given in the output and log file.

* The second one is to use the external code wannier90 ([www.wannier.org](http://www.wannier.org)) 
  that calculates Maximally Localized Wannier Functions (MLWF). 
  After a ground state calculation the Wannier90 code will obtain the MLWFs requiring just two ingredients: 

  * The overlaps between the cell periodic part of the Bloch states.
  * As a starting guess the projection of the Bloch states onto trial localized orbitals.

What ABINIT do is to take the Bloch functions from a ground state calculation
and compute these two ingredients. Then, Wannier90 is run. Wannier90 is
included as a library and ABINIT and the process is automatic, so that in a
single run you can do both the ground state calculation and the computation of
MLWFs. The input variables [[prtwant]], [[w90iniprj]] and [[w90prtunk]] are
related to the use of the wannier90 librairy.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:wannier90|The tutorial on Wannier90]] deals with the Wannier90 library to obtain Maximally Localized Wannier Functions.

