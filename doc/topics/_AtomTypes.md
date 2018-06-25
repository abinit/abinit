---
description: How to specify the types of atoms that form the system
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to specify the types of atoms that form the system.

## Introduction

ABINIT needs to know the different types of atoms that form the system.
The atoms assembled in a molecule or a solid are physically specified by their
nuclear charge (and their mass for dynamical properties).

However, in a pseudopotential or PAW approach, the knowledge of the nuclear
charge does not define the potential felt by the electron, only the atomic
data file (pseudopotential or PAW) defines it. Thus, in addition to the number
of types of atoms [[ntypat]], and their nuclear charge [[znucl]], ABINIT
requires to know the pseudopotential/PAW to use for each type of atom. The
latters are specified in the [[help:abinit#intro1|"files" file]]. Unless
alchemical potentials are used (see later), the number of pseudopotentials to
be read, [[npsp]], is the same as [[ntypat]]. Note that one cannot mix norm-
conserving pseudopotentials with PAW atomic data files in a single ABINIT run,
even for different datasets. One has to stick either to norm-conserving
pseudopotentials or to PAW.  
More on the pseudos/PAW in [[topic:PseudosPAW]].

ABINIT also has a default table of the atomic masses, but this can be
superceded by specifying [[amu]].

**Alchemical potentials**

For norm-conserving pseudopotentials, ABINIT can mix the pseudopotentials, as
described in [[https://wiki.abinit.org/doku.php?id=developers:pseudos|the ABINIT wiki]], 
to create so-called "alchemical potentials", see [[mixalch]].  
In this case, the number of pseudopotentials to be given, [[npsp]], will
usually be larger than the number of types of atoms, [[ntypat]]. Using
alchemical potentials makes sense to treat alloys in which similar ions are
present, and whose specific chemical properties are not crucial for the
property of interest. Usually it is done only for isovalent species, and ions
of quite similar radii. It is a reasonable interpolation technique for the
electronic properties.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:ffield|Tutorial on polarization and finite electric fields]]. Polarization, and responses to finite electric fields for AlAs. In the present topic, it is an example of the definition of several atom types ... 

