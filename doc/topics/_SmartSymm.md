---
description: How to use the symetry information to build the system from the irreducible part of the primitive cell
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to use the symetry information to build the system from the irreducible part of
the primitive cell with the ABINIT package.

## Introduction

Sometimes, the user knows the space group of the system, the conventional cell
vectors, as well as the positions of atoms in the asymmetric (irreducible)
part of the cell. From such data, ABINIT can generate the usual primitive cell
vectors, as well as the coordinates of all the atoms in this cell.

This is activated if [[spgroup]]!=0 and [[brvltt]]=-1. The user needs to
specify the number of atoms to be read [[natrd]] in the asymmetric part of the
cell, as well as the expected number of atoms [[natom]] in the primitive cell.
Some additional information on the axes orientation [[spgaxor]] and the cell
origin [[spgorig]] might also have to be given. 
See [[help:spacegroup|the space group help file]]

The specification of a magnetic space group is even possible
(antiferromagnetic). See [[spgroupma]] and [[genafm]].

See also the [[topic:UnitCell]].


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

