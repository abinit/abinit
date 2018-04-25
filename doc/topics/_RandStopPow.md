---
description: How to perform random stopping power calculation
authors: FB
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform random stopping power calculation with the ABINIT package.

## Introduction

The slowing down of a swift charged particle inside condensed matter has been
a subject of intense interest since the advent of quantum-mechanics. The
Lindhard formula [[cite:Lindhard1954]] that gives the polarizability of the
free electron gas has been developed specifically for this purpose. The
kinetic energy lost by the impinging particle by unit of path length is named
the stopping power. For large velocities, the stopping power is dominated by
its electronic contribution: the arriving particle induces electronic
excitations in the target. These electronic excitations in the target can be
related to the inverse dielectric function ε-1( **q** ,ω) provided that linear
response theory is valid.

As a consequence, the electronic stopping power randomized over all the
possible impact parameters reads

S( **v** ) = (4π Z2/N **q** Ω| **v** |)∑ **q** ∑ **G** Im{- ε-1[ **q** ,
**v.** ( **q** + **G** )]} ( **v.** ( **q** + **G** )/| **q** + **G** |2),

where Z and **v** are respectively the charge and the velocity of the
impinging particle,  Ω is the unit cell volume, N **q** is the number of **q**
-points in the first Brillouin zone, and **G** are reciprocal lattice vectors.

Apart from an overall factor of 2, this equation is identical to the formula
published [[cite:Campillo1998]].

The GW module of ABINIT gives access to the full inverse dielectric function
for a grid of frequencies ω. Then, the implementation of the above equation is
a post-processing employing a spline interpolation of the inverse dielectric
function in order to evaluate it at ω= **v.** ( **q** + **G** ). The energy
cutoff on **G** is governed by the [[ecuteps]], as in the GW module. The
integer [[npvel]] and the cartesian vector [[pvelmax]] control the
discretization of the particle velocity.

Note that the absolute convergence of the random electronic stopping power is
a delicate matter that generally requires thousands of empty states together
with large values of the energy cutoff.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

