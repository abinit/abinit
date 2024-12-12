---
description: How to calculate the effective Coulomb interaction
authors: BAmadon
---
<!--- This  file has been generated automatically from the corresponding _* source file. DO NOT EDIT. Edit the source file instead. -->

This page gives hints on how to calculate the effective Coulomb interaction with the ABINIT package.

## Introduction

DFT+U as well as DFT+DMFT requires as input values the effective Coulomb
interaction. Two ways to compute them are available in ABINIT.

Firstly, the constrained Random Phase Approximation [[cite:Aryasetiawan2004]]
[[ucrpa]] allows one to take into account the screening of the Coulomb
interaction between correlated electrons, by non-interacting electrons. For
non-entangled bands ([[ucrpa]]= 1), the bands excluded from the polarisability
can be specified either by a band index ([[ucrpa_bands]]) or an energy window
([[ucrpa_window]]) [[cite:Amadon2014]].

For entangled bands ([[ucrpa]]= 2}), the scheme used in ABINIT
[[cite:Shih2012]], [[cite:Sakuma2013]],[[cite:Amadon2014]] uses a band and
k-point dependent weight to define the polarisability, using Wannier orbitals
as correlated orbitals.

This method is well adapted to compute the effective interaction for the same
orbitals used in DFT+DMFT. To use the same orbitals as in DFT+U, the Wannier
functions can be ajusted such that the bare interaction is close to the bare
interaction of atomic orbitals as used in DFT+ _U_ (see tutorial).

Secondly, a linear response method [[cite:Cococcioni2005]] is implemented. The
implementation is not yet in production. The implementation in ABINIT takes
into account the truncated atomic orbitals from PAW and therefore differs from
the original work [[cite:Cococcioni2005]] treating full atomic orbitals. In
particular, considerably higher effective values for U are found.



## Related Input Variables

*compulsory:*

- [[abinit:ucrpa]]  calculation of the screened interaction U with the Constrained RPA method
 
*basic:*

- [[abinit:ucrpa_bands]]  For the calculation of U with the Constrained RPA method, gives correlated BANDS
- [[abinit:ucrpa_window]]  For the calculation of U with the Constrained RPA method, gives energy WINDOW
 

## Selected Input Files

*v7:*

- [[tests/v7/Input/t23.abi]]
- [[tests/v7/Input/t24.abi]]
- [[tests/v7/Input/t25.abi]]
- [[tests/v7/Input/t78.abi]]
- [[tests/v7/Input/t79.abi]]
 

## Tutorials

* The [[tutorial:ucalc_crpa|tutorial on the calculation of effective interactions U and J by the cRPA method]] shows how to determine the U value with the constrained Random Phase Approximation [[cite:Aryasetiawan2004]] using projected Wannier orbitals. Prerequisite: DFT+U.
* The [[tutorial:lruj|tutorial on the computation of U from linear response]] for DFT+U shows how to determine the U value with the linear response method [[cite:Cococcioni2005]], to be used in the DFT+U approach. Prerequisite: DFT+U.

