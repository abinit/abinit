---
description: How to generically perform DFPT calculations
authors: MT
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to generically perform DFPT calculations with the ABINIT package.

## Introduction

Density-Functional Perturbation Theory (DFPT) allows one to address a large
variety of physical observables. Many properties of interest can be computed
directly from the derivatives of the energy, without the use of finite
differences: phonons modes, elastic tensors, effective charges, dielectric
tensors, etc... Even non-linear properties can be computed, like the Raman
intensities (for the latter, see [[topic:nonlinear]])..

A DFPT calculation workflow is conducted as follows:

* Run a Ground-State calculation in order to extract the Kohn-Sham pseudo wave-functions; these must be extremely well converged.
* If necessary, e.g., for the application of the derivative of the Hamiltonian with respect to an electric field, determine the derivatives of the wave functions with respect to the wave vector **k** , and keep them in a file. The keyword [[rfddk]] is used to perform this type of calculation.
* Compute the 2nd-order derivative matrix (i.e., 2nd derivatives of the energy with respect to different perturbations λ). This can be done thanks to the keywords [[rfphon]] (λ=atomic displacement), [[rfstrs]] (λ=strain), [[rfelfd]] (λ=electric field) or [[rfmagn]] (λ=magnetic field). 
* Launch the anaddb tool (distributed with ABINIT) to analyse the derivative database and compute relaxed tensors and thermodynamical properties.

Note that for PAW calculation, when performing the post-processing with
anaddb, it is recommended to include all the keywords enforcing the sum rules
(acoustic sum and charge neutrality). Indeed the PAW formalism involves, for
each atom, the calculation of a large number of real space integrals, whose
numerical effect may be to break the translational invariance.

Thanks to the locality provided by PAW partial wave basis, it is possible to
perform response function calculations for correlated electron materials. The
DFT+U formalism is usable without any restriction for the PAW+DFPT
calculations.

All the tutorials dedicated to response functions can be followed both with
norm-conserving pseudopotentials and with PAW atomic datasets.

DFPT in ABINIT is implemented for non-magnetic, collinear as well as non-collinear systems [[nspden]]=1, 2 as well as 4.
However, the treatment of the strain perturbation is not yet implemented with [[nspden]]=4 (non-collinear systems).

More detailed explanations to perform a response calculation are given in the [[help:respfn]].


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:rf1|The tutorial Response-Function 1 (RF1)]] presents the basics of DFPT calculations within ABINIT. The example given is the study of dynamical and dielectric properties of AlAs (an insulator): phonons at Gamma, dielectric constant, Born effective charges, LO-TO splitting, phonons in the whole Brillouin zone. The creation of the "Derivative Data Base" (DDB) is presented.
* [[tutorial:rf2|The tutorial Response-Function 2 (RF2)]] presents the analysis of the DDBs that have been introduced in the preceeding tutorial RF1. The computation of the interatomic forces and the computation of thermodynamical properties is an outcome of this tutorial.
* [[tutorial:elastic|The tutorial on the elastic properties]] presents the computation with respect to the strain perturbation and its responses: elastic constants, piezoelectricity.
* [[tutorial:paral_dfpt|Parallelism of response-function calculations]]. Additional information to use the DFPT in parallel.
* [[tutorial:nlo|Tutorial on static non-linear properties]]. Electronic non-linear susceptibility, non-resonant Raman tensor, electro-optic effect.

