---
description: How to perform a Tdep calculation
authors: YB, FBottin
---
<!--- This  file has been generated automatically from the corresponding _* source file. DO NOT EDIT. Edit the source file instead. -->

This page gives hints on how to perform thermodynamic, elastic and transport properties calculations including explicit temperature effects with the ABINIT package.  

User guide: [[pdf:TDEP_Guide| TDEP guide]]  
Theory: [[pdf:TDEP_Paper|TDEP paper]]

## Introduction

The Temperature Dependent Effective Potential (TDEP) method
has been developped by O. Hellman *et al.* [[cite:Hellman2011]],
[[cite:Hellman2013]], [[cite:Hellman2013a]] in 2011 and the |aTDEP| implementation
in ABINIT has been performed and used for the first time in 2015 by
J. Bouchet and F. Bottin [[cite:Bouchet2015]], [[cite:Bouchet2017]].

The capture of thermal effects in solid state physic is a long standing
issue and several stand-alone or post-process computational codes are
available. Using different theoretical frameworks, they propose to provide
some thermodynamic quantities involving the so called anharmonic effects.
|aTDEP| calculation can produce almost all the temperature-dependent
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

*basic:*

- [[tdep:amu]]  Atomic masses in Mass Units
- [[tdep:angle]]  ANGLE alpha
- [[tdep:brav]]  BRAVais
- [[tdep:multiplicity]]  MULTIPLICITY
- [[tdep:natom]]  NATOM
- [[tdep:natom_unitcell]]  NATOM in the UNITCELL
- [[tdep:nstep_max]]  NSTEP at MAX
- [[tdep:nstep_min]]  NSTEP at MIN
- [[tdep:ntypat]]  NTYPAT
- [[tdep:rcut]]  Radius CUToff
- [[tdep:rprimd]]  RPRIMD
- [[tdep:temperature]]  TEMPERATURE
- [[tdep:typat]]  TYPAT
- [[tdep:typat_unitcell]]  TYPAT in the UNITCELL
- [[tdep:xred_unitcell]]  XRED in the UNITCELL
 
*expert:*

- [[tdep:bzpath]]  Brillouin Zone PATH
- [[tdep:dosdeltae]]  DOS delta Energy
- [[tdep:enunit]]  ENergy UNIT
- [[tdep:ngqpt1]]  Number of Grid points for Q PoinTs generation (coarse)
- [[tdep:ngqpt2]]  Number of Grid points for Q PoinTs generation (fine)
- [[tdep:order]]  ORDER for the IFC
- [[tdep:slice]]  SLICE
- [[tdep:use_ideal_positions]]  USE IDEAL POSITIONS
 

## Selected Input Files

*v8:*

- [[tests/v8/Input/t37.in]]
 

