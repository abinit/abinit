---
description: How to compute the polarisation and take into account a finite homogeneous electric field
authors: JZ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to compute the polarisation and take into account a finite homogeneous electric
field with the ABINIT package.

## Introduction

The effect of an homogeneous static electric field on an insulator may be
treated in ABINIT from two perspectives. One is perturbative, and yields the
susceptibility in the form of the second derivative of the total energy with
respect to the electric field, at zero field strength (see [[topic:DFPT]]).

ABINIT can also be used to compute the effect of an electric field of finite
amplitude, using techniques from the Modern Theory of Polarization
[[cite:Resta1994]],[[cite:Nunes2001]],[[cite:Souza2002]]. The latter is based
on the notion of "Berry phase". In this approach, the total energy to minimize
includes the contribution due to the interaction of the external electric
field with the material polarization **P** Tot, as follows:

E = E0 \- Ω **P** Tot **.E** , where E0 is the usual ground state energy
obtained from Kohn-Sham DFT in the absence of the external field **E** , **P**
Tot is the polarization, made up of an ionic contribution and an electronic
contribution, and Ω the volume of the unit cell.

Some details of the implementation of The Modern Theory of Polarization in
ABINIT are given in [[cite:Gonze2016|the 2016 ABINIT publication]].

In the NCPP case, the electric field has no additional contribution to the
Hellmann-Feynman forces, because the electronic states do not depend
explicitly on ionic position [[cite:Souza2002]]. In the PAW case however, as
the projectors do depend on ion location, an additional force and additional
stresses terms arise [[cite:Zwanziger2012]].

The generalisation to fixed D-field or fixed reduced fields are also
available, as described in M. Stengel, N.A. Spaldin and D. Vanderbilt, Nat. Phys. 5,304 (2009).

The polarization and finite electric field calculation in ABINIT is accessed
through the variables [[berryopt]] and [[efield]]. In addition, displacement
fields and mixed boundary conditions (a mix of electric field and displacement
field) can be computed as well.


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* The [[tutorial:ffield]|tutorial on polarization and finite electric field deals with the computation of the polarization of an insulator (e.g. ferroelectric, or dielectric material) thanks to the Berry phase approach, and also presents the computation of materials properties in the presence of a finite electric field (also thanks to the Berry phase approach).
