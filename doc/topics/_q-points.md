---
description: How to set parameters related to the phonon wavevectors (q-points) in DFPT calculations
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to set parameters related to the phonon wavevectors (q-points) in DFPT
calculations with the ABINIT package.

## Introduction

Like the electronic wavefunctions, the collective atomic displacements that
are eigenmodes of the corresponding periodic Hamiltonian can be characterized
by a wavevector, denoted q-point.

In ABINIT, DFPT calculations for one dataset are done for one specific
q-point, that must be specified. In the simplest case, the user gives the
corresponding q-point for each dataset, setting [[nqpt]]=1 and specifying the
corresponding single [[qpt]]. However, very often, it is needed to run
calculations for dozens or hundreds of q-points. Hence, the following
mechanism has been set: the use can specify a set of q points, using input
variables similar to the k-points, and then, for each dataset, the number of
the q-point in the set is indicated thanks to [[iqpt]]. This applies to the
generation of q-point grids as well as to q-point paths to produce phonon band
structures.

The input variables for specifying q-points in ANADDB are specified in
[[topic:Phonons]] and [[topic:PhononBands]].



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

