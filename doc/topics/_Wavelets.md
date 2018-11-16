---
description: How to perform calculations on a wavelet basis
authors: MT
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform calculations on a wavelet basis with the ABINIT package.

## Introduction

A wavelet basis (instead of a plane wave basis) can be used in ABINIT. With a
wavelet basis, one can perform basic static DFT calculations with selected
norm-conserving pseudopotentials (HGH or GTH pseudopotentials
[[cite:Genovese2008]]), but also with PAW atomic data [[cite:Rangel2016]]).
Available also : the finite size corrections to the total energy, restart on
wavefunctions following the ETSF norm and geometry relaxation using BFGS.
Molecular dynamic is also available for test purposes.

However, DFPT or excited-state calculations (except Δ-SCF) cannot be
performed.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:paral_gswvl|Parallelism for ground-state calculations, with wavelets]] presents the parallelism of ABINIT, when wavelets are used as a basis function instead of planewaves, for the computation of total energy, density, and ground state properties

