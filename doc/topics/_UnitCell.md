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
be positive. If it is not the case, exchange two of them with the associated reduced coordinates.
More information about the way the real space lattice, the reciprocal lattice,
and symmetries are defined in ABINIT can be found [[pdf:geometry|here]].

Also note that that Abinit space group routines uses by default strict tolerances for the 
recognition of the symmetry operations. 
This means that lattice vectors and atomic positions must be give with enough figures 
so that the code can detect the correct space group.
This is especially true in the case of hexagonal or rhombohedral lattices.
Remember that one can specify rational numbers with the syntax:

    xred 1/3 1/3 1/3

instead of the less precise:

    xred 0.33333 0.33333 0.33333


Using [[acell]] and [[angdeg]] instead of [[rprimd]] may solve possible issues 
with the space group recognition.

If your input parameters correspond to a high-symmetry structure but the numerical values at hand
are *noisy*, you may want to increase the value of [[tolsym]] in the input file  
so that Abinit will resymmetrize automatically the input parameters.

Finally, one can use the [[structure]] input variable to initialize the crystalline geometry 
from an external file.
For instance, it is possible to read the crystalline structure from an external netcdf file or 
other formats such as POSCAR without having to specify the value of 
[[natom]], [[ntypat]], [[typat]] and [[znucl]].


**Smart symmetriser**

ABINIT has also a smart symmetriser capability, when [[spgroup]]!=0 and
[[brvltt]]=-1. In this case, the CONVENTIONAL unit cell must be input through
the usual input variables [[rprim]], [[angdeg]], [[acell]] and/or
[[scalecart]]. ABINIT will fold the conventional unit to the primitive cell,
and also generate all the nuclei positions from the irreducible ones. 
See [[topic:SmartSymm]].


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:base4|Fourth basic tutorial]] Determination of the surface energy of aluminum (100): changing the orientation of the unit cell.

