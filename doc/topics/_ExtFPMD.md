---
description: How to enable the Extended FPMD method for high temperature simulations
authors: A. Blanchet
---
<!--- This is the source file for this topics. Can be edited. -->
This page gives hints on how to perform calculation with the extended FPMD model.

## Introduction

Extended First-Principles Molecular Dynamics (Ext. FPMD) method allows to perform high temperature simulations from few Kelvins to thousands of eVs, by drastically reducing the needed number of bands for high temperature simulations [[cite:Blanchet2020]]. The implementation and usage is described in [[cite:Blanchet2022]].

High energy orbitals are replaced with pure single plane waves description based on the Fermi gas model. Bands from 1 to [[nband]] are treated with the complete plane waves basis set as usual, and the rest of occupied bands from [[nband]] to the infinity are treated with the Fermi gas model. Contributions to the energy, entropy, stresses, number of electrons and chemical potential are computed automatically after enabling the model with variable [[useextfpmd]]. Conventional convergency studies are still needed to get accurate results.

Contributions to the number of electrons and to the energy are explicitly shown in the *_GSR.nc* output file with key **nelect_extfpmd** and **e_extfpmd** (**edc_extfpmd** for the double counting term). The energy shift (resulting from the constant background potential) is also printed in the *_GSR.nc* output file with key **shiftfactor_extfpmd**.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}
