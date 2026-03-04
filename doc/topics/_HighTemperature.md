---
description: How to run high temperature calculations
authors: AB
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to run high temperature calculations with the ABINIT package.

## Introduction


## Orbital free molecular dynamics (Recursion method)
It is possible to use Thomas-Fermi kinetic functional (explicit functional of the density) or Thomas-Fermi-Weizsacker kinetic functional (with gradient corrections) instead of Kohn-Sham kinetic energy functional (implicit functional of the density through Kohn-Sham wavefunctions). See [[cite:Perrot1979]]. The Recursion method may be used in order to compute electronic density, entropy, Fermi energy and eigenvalues energy. This method computes the density without computing any orbital, is efficient at high temperature, with a efficient parallelization (almost perfect scalability).
At present, it only works for local pseudopotentials, severely restricting the use of this method.

## Extended first-principles molecular dynamics (ExtFPMD model)
Extended First-Principles Molecular Dynamics model (Ext. FPMD), also known as extended DFT, allows to perform high temperature simulations from few Kelvins to thousands of eVs, by drastically reducing the needed number of bands for high temperature simulations [[cite:Blanchet2020]]. The implementation and usage is described in [[cite:Blanchet2022]]. ExtFPMD calculation is enabled by setting [[useextfpmd]]=1 (see the input variable description for further information).

High energy orbitals are replaced with pure single plane waves description based on the homogeneous electron gas model. This model can be set to its non-relativistic or relativistic version (controlled with the input variable [[extfpmd_rel]]). Bands from 1 to [[nband]] are treated with the complete plane waves basis set as usual, and the rest of occupied bands from [[nband]] to the infinity are treated with the homogeneous electron gas model (analytic). Contributions to the electron density, energy, entropy, stresses, number of electrons and chemical potential are computed automatically after enabling the model with [[useextfpmd]]. Conventional convergency studies are still needed to get accurate results (especially on the parameter [[nband]]).

To enhance the SCF cycle convergence, it is recommended to set a band buffer with [[nbdbuf]].

A somewhat refined version of the ExtFPMD model aims at replacing the missing orbitals with a contribution for each point of the real space grid. This version is known as Hybrid Kohn-Sham + Thomas-Fermi scheme. As evaluation of Fermi-Dirac integrals are required for each grid point, this version can be significantly slower with large unit cells. Hybrid Kohn-Sham + Thomas-Fermi scheme is available with [[useextfpmd]]=10.

Numerical integration of the generalized Fermi-Dirac incomplete integrals are implemented with the method described in [[cite:Aparicio1998]].

Contributions to the number of electrons and to the energy are explicitly shown in the *_GSR.nc* output file with key **nelect_extfpmd** and **e_extfpmd** (**edc_extfpmd** for the double counting term). The energy shift (resulting from the constant background potential) is also printed in the *_GSR.nc* output file with the key **extfpmd_eshift**.

## Thermal exchange-correlation functionals

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

