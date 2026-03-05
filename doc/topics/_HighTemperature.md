---
description: How to run high temperature calculations
authors: AB
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to run high temperature calculations with the ABINIT package.

## Introduction
The Mermin finite-temperature formulation of Kohn-Sham DFT has proven to be successful in handling the complexity of the warm dense matter regime, for electronic temperatures close to the Fermi temperature of the material TF = EF /kB and even above it. Unfortunately, plane-wave-based DFT is rather limited to low temperatures because of the orbital wall [[cite:Blanchet2020]]. At high temperatures, the Fermi-Dirac distribution imposes to consider a large number of weakly occupied states in the dense high energy continuum. Computing electronic properties, such as the equation of state (EOS), for temperatures higher than T F is then completely unreachable, even on the largest supercomputers.

Moreover, running calculations at high temperatures requires custom crafted pseudopotentials or atomic datas in the PAW case. Indeed, since core electrons will be thermally activated, their energy will shift from their cold value, and the core orbitals will also be thermally ionized (not fully occupied anymore). In the high temperature limit, all electrons must be unfrozen from the core. Crafting atomic data suitable for high temperature calculations is achievable with the ATOMPAW generator [[topic:PseudosPAW]].

## Orbital free molecular dynamics (Recursion method)
Among the different methods aiming at bypassing the orbital wall, the orbital free methods completely supress the orbitals by approximating the kinetic energy as a functional of the density. In ABINIT, it's possible to use Thomas-Fermi kinetic functional (explicit functional of the density) or Thomas-Fermi-Weizsacker kinetic functional (with gradient corrections) instead of Kohn-Sham kinetic energy functional (implicit functional of the density through Kohn-Sham wavefunctions), through the recursion method by setting the input variable [[tfkinfunc]] (see [[cite:Perrot1979]]).

The Recursion method may be used in order to compute electronic density, entropy, Fermi energy and eigenvalues energy. This method computes the density without computing any orbital, is efficient at high temperature, with a efficient parallelization (almost perfect scalability). Since this kind of calculation is orbital free, the effects of the ionization of the inner electron shells at high temperature is not described. Furthermore, the kinetic energy density functional approximation is kown not to reproduce the material properties at low temperatures.

At present, the recursion method only works for local pseudopotentials, severely restricting the use of this method.

## Extended first-principles molecular dynamics (ExtFPMD model)
Extended First-Principles Molecular Dynamics model (Ext. FPMD), also known as extended DFT, allows to perform high temperature simulations from few Kelvins to thousands of eVs, where the plasma is completely ionized, by drastically reducing the number of bands required for high temperature simulations, bypassing the orbital wall [[cite:Blanchet2020]]. The implementation and usage is described in [[cite:Blanchet2022]]. ExtFPMD calculation is enabled by setting [[useextfpmd]]=1 (see the input variable description for further information). Unlike orbital free methods, the extended DFT model reproduces ionization shell effects, because core orbitals are still described with Kohn-Sham wavefunctions.

High energy orbitals are replaced with pure single plane waves description based on the homogeneous electron gas model. This model can be set to its non-relativistic or relativistic version (controlled with the input variable [[extfpmd_rel]]). Bands from 1 to [[nband]] are treated with the complete plane waves basis set (Kohn-Sham wave functions), and the rest of occupied bands from [[nband]] to the infinity are treated with the homogeneous electron gas model (analytic). Contributions to the electron density, energy, entropy, stresses, number of electrons and chemical potential are computed automatically after enabling the model with [[useextfpmd]]. Conventional convergency studies are still needed to get accurate results (especially on the parameter [[nband]]). To enhance the SCF cycle convergence, it is recommended to set a band buffer with [[nbdbuf]].

A somewhat refined version of the ExtFPMD model aims at replacing the missing orbitals with a contribution for each point of the real space grid, instead of an homogeneous gas. This version is known as Hybrid Kohn-Sham + Thomas-Fermi scheme. As evaluation of Fermi-Dirac integrals are required for each grid point, this version can be significantly slower with large unit cells. Hybrid Kohn-Sham + Thomas-Fermi scheme is available with [[useextfpmd]]=10.

Numerical integration of the generalized Fermi-Dirac incomplete integrals are implemented with the method described in [[cite:Aparicio1998]]. Contributions to the number of electrons, the energy and the entropy are explicitly shown in the *_GSR.nc* output file with key **nelect_extfpmd**, **e_extfpmd** (**edc_extfpmd** for the double counting term), and **entropy_extfpmd**. The energy shift (resulting from the constant background potential) is also printed in the *_GSR.nc* output file with the key **extfpmd_eshift**. ExtFPMD contributions are not shown when computing the density of states [[prtdos]].

## Relaxed core method
Section incoming...

## Thermal exchange-correlation functionals
Section incoming...

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

