---
description: How to perform a PIMD calculation
authors: GG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform a PIMD calculation with the ABINIT package.

## Introduction

Path-Integral Molecular Dynamics (PIMD) is a technique allowing to simulate
the quantum fluctuations of the nuclei at thermodynamic equilibrium
[[cite:Marx1996]]. It is implemented in ABINIT in the NVT ensemble since v7.8.2.

In the Path-Integral formalism of quantum statistical mechanics, the (quantum)
nuclei are replaced by a set of images (beads) treated by means of classical
mechanics, and interacting with each other through a specific effective
potential. In the limit of an infinite number of beads, the quantum system and
this many-beads classicle system have the same partition function, and thus
the same static observables. In PIMD, the classical system of beads is
simulated by standard Molecular Dynamics. The PIMD equations of motion are
integrated by using the Verlet algorithm. At each time step, a ground state
DFT calculation is performed for each image. PIMD can be used with any XC
functional and works in the PAW framework as well as in the norm-conserving
pseudopotential (NCPP) case.

PIMD in ABINIT follows the set of numerical schemes developed by several
authors in the 90's [[cite:Marx1996]], [[cite:Tuckerman1996]]. PIMD in the
canonical ensemble needs specific thermostats to ensure that the trajectories
are ergodic: the Nose-Hoover chains are implemented, as well as the Langevin
thermostat (controlled by the value of [[imgmov]]). Also, it is possible to
use coordinate transformations (staging and normal mode), that are controlled
by the keyword [[pitransform]]. In standard equilibrium PIMD, only static
observables are relevant (quantum time-correlation functions are not
accessible): the masses associated to the motion of the beads are controlled
by the keyword [[pimass]], whereas the true masses of the atoms are given by
[[amu]]. The values given in [[pimass]] are used to fix the so-called
fictitious masses [[cite:Marx1996]]. In the case where a coordinate
transformation is used, the fictitious masses are automatically fixed in the
code to match the so-called staging masses or normal mode masses. The number
of time steps of the trajectory is controlled by [[ntimimage]], the initial
and thermostat temperature by [[mdtemp]]. Except if specified, the images in
the initial configuration are assumed to be at the same position, and a random
distribution of velocities is applied (governed by [[mdtemp]]) to start the
dynamics.

At each time step, ABINIT delivers in the output file:

* (i) information about the ground state DFT calculation of the ground state for each image
* (ii) the instantaneous temperature, the instantaneous energy as given by the primitive and virial estimators, and the pressure tensor as given by the primitive estimator.

The number of images (keyword [[nimage]]) is associated to a specific
parallelization level (keyword [[npimage]]).

PIMD has been used with ABINIT to reproduce the large isotope effect on the
phase transition between phase I and phase II of dense hydrogen
[[cite:Geneste2012]], and also some aspects of diffusion at low and room
temperature in proton-conducting oxides for fuel cells [[cite:Geneste2015]].
PIMD in the NPT ensemble is not available yet.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

