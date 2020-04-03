---
description: How to to specify a crystal, with atomic positions and symmetries
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to to specify a crystal, with atomic positions and symmetries with the ABINIT package.

## Introduction

In addition to the [[topic:UnitCell|Specification of the unit cell]] and
[[topic:AtomTypes|Atom types]], ABINIT must know the number of atoms inside
the cell, their type, and position. This is described by [[natom]], [[typat]]
and one of [[xred]] or [[xcart]].

ABINIT can automatically detect the Bravais lattice and space group, and
generate symmetries (e.g. [[nsym]], [[symrel]], [[tnons]]), from the primitive
cell and the position of atoms (provided they are not too inaccurate, see
[[tolsym]]). For this purpose, in the magnetic case, ABINIT will also take
into account the input atomic spin, through the knowledge of [[spinat]].

Alternatively, ABINIT can start from the specification of symmetries (either
from [[spgroup]] or from the list of symmetries -
[[nsym]], [[symrel]], [[tnons]]) and generate the atomic positions from the
asymmetric (irreducible) part of the primitive cell. This is described in the
[[topic:SmartSymm|Smart Symmetrizer]] topic.

ABINIT can treat antiferromagnetic symmetry operations, see [[symafm]].

In ABINIT, a database with the 230 spatial groups of symmetry (see
[[spgroup]]) and the 1191 Shubnikov anti-ferromagnetic space groups is present
(see also [[spgroupma]] and [[genafm]]).

There is also a (non-graphical) atom manipulator in ABINIT, see [[topic:AtomManipulator]].

ABINIT can read XYZ files, see [[xyzfile]].

Atomic positions can also be generated at random, see [[random_atpos]].

Details about the way the crystal structure is defined in ABINIT can be found [[theory:geometry|here]].

!!! tip

    If |AbiPy| in installed on your machine, you can use the |abistruct| script
    to automate several operations related to crystalline structures.
    Further details about the python API are available in the |AbipyStructureNb|.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:base1|The tutorial 1]] deals with the H2 molecule : get the total energy, the electronic energies, the charge density, the bond length, the atomisation energy 
* [[tutorial:base2|The tutorial 2]] deals again with the H2 molecule: convergence studies, LDA versus GGA 
* [[tutorial:base3|The tutorial 3]] deals with crystalline silicon (an insulator): the definition of a k-point grid, the smearing of the cut-off energy, the computation of a band structure, and again, convergence studies ...
* [[tutorial:base4|The tutorial 4]] deals with crystalline aluminum (a metal), and its surface: occupation numbers, smearing the Fermi-Dirac distribution, the surface energy, and again, convergence studies ...

