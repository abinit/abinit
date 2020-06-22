---
description: How to perform a Tdep calculation
authors: YB, FBottin
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform thermodynamic, elastic and transport properties calculations including explicit temperature effects with the ABINIT package.  

User guide: [[pdf:a-TDEP_Guide| a-TDEP guide]]  
Theory: [[pdf:a-TDEP_Paper|a-TDEP paper]]

## Introduction

The Temperature Dependent Effective Potential (TDEP) method
has been developped by O. Hellman *et al.* [[cite:Hellman2011]],
[[cite:Hellman2013]], [[cite:Hellman2013a]] in 2011 and the |a-TDEP| implementation
in ABINIT has been performed and used for the first time in 2015 by
J. Bouchet and F. Bottin [[cite:Bouchet2015]], [[cite:Bouchet2017]].

The capture of thermal effects in solid state physic is a long standing
issue and several stand-alone or post-process computational codes are
available. Using different theoretical frameworks, they propose to provide
some thermodynamic quantities involving the so called anharmonic effects.
|a-TDEP| calculation can produce almost all the temperature-dependent
thermodynamic quantities you want, from a single *ab initio*
molecular dynamic (AIMD) trajectory and by means of a Graphical User
Interface (GUI) very easy to use ([[https://github.com/abinit/abiout|AGATE]]).

The original TDEP method [[cite:Hellman2011]] is implemented in ABINIT.
In particular, various algorithms can be used to obtain the Interatomic Force Constants (IFC).
The 2nd-order (and soon 3rd-order) IFCs are produced self-consistently using a least-square
method fitting the AIMD forces on a model Hamiltonian function of the displacements.  
Numerous thermodynamic quantities can be computed starting from the
2nd order IFCs. The 1st one is the phonon spectra, from which a large
number of other quantities flow : internal energy, entropy, free energy, specific heat...
The elastic constants and other usual elastic moduli (the bulk,
shear and Young moduli) can be also produced at this level. Using the 3rd
order IFCs, we could extract the Gruneisen parameter, the thermal
expansion, the sound velocities... and in particular, how to take into account
the anisotropy of the system within.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

