---
description: How to specify the unit cell
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to specify the unit cell with the ABINIT package.

## Introduction

ABINIT needs three dimensioned non-coplanar vectors, forming the unit cell, to
set up the real space lattice.

An initial set of three vectors, specified in real space by [[rprim]] or as
unit vectors with angles [[angdeg]], are dimensioned in a second step using
scaling factors as specified by [[acell]] or by rescaling their cartesian
coordinates, as specified by [[scalecart]]. Internally, only the final result,
[[rprimd]] matters. The most detailed explanation can be found by looking at
[[rprim]].

Note that ABINIT expects the mixed product of the three vectors (R1xR2).R3 to
be positive. If it is not the case, exchange two of them.

More information about the way the real space lattice, the reciprocal lattice,
and symmetries are defined in ABINIT can be found [[pdf:geometry|here]].

**Smart symmetriser**

ABINIT has also a smart symmetriser capability, when [[spgroup]]!=0 and
[[brvltt]]=-1. In this case, the CONVENTIONAL unit cell must be input through
the usual input variables [[rprim]], [[angdeg]], [[acell]] and/or
[[scalecart]]. ABINIT will fold the conventional unit to the primitive cell,
and also generate all the nuclei positions from the irreducible ones. See
[[topic:SmartSymm]].



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:base4|Fourth basic tutorial]] Determination of the surface energy of aluminum (100): changing the orientation of the unit cell.

