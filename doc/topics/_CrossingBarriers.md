---
description: How to calculate crossing barriers
authors: GG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to calculate barriers for crossings with the ABINIT package.

## Introduction

The knowledge of geometries at which crossings between two electronic states happen,
with minimal energy, or geometries at which the energy difference between the ground state and the excited
state is small, and the energy is still low, plays an important role
in the study of non-radiative transitions.

It is possible to formulate the search for such geometries in terms of minimisation
of a functional that is the linear combination of the energy of the two states at the same geometry, with Lagrange multipliers (Ref. to be given).
This is also related with a simple approach to Ensemble DFT: just make a linear combination of the DFT energies, the XC correlation
energy being not computed with a single common density, but from each density separately.

In ABINIT, with [[imgmov]]==6, it is possible to deal with such 
linear combination of systems with the same geometry, but differing occupation factors [[occ]].
It is possible to find the geometry at which the resulting energy is minimal, for a given value of the mixing factors [[mixesimgf]].
Set [[nimage]]=2, and set the occupation numbers for image 1 to the ground-state occupations, and for image 2 to the excited-state occupations.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

