---
description: How to perform a GW calculation, including self-consistency
authors: MG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform a GW calculation, including self-consistency with the ABINIT package.

## Introduction

DFT performs reasonably well for the determination of structural properties,
but fails to predict accurate band gaps. A more rigorous framework for the
description of excited states is provided by many-body perturbation theory
(MBPT) [[cite:Fetter1971]], [[cite:Abrikosov1975]], based on the Green's
functions formalism and the concept of quasi-particles [[cite:Onida2002]].

Within MBPT, one can calculate the quasi-particle (QP) energies, E, and
amplitudes, Ψ, by solving a nonlinear equation involving the non-Hermitian,
nonlocal and frequency dependent self-energy operator Σ.

This equation goes beyond the mean-field approximation of independent KS
particles as it accounts for the dynamic many-body effects in the electron-
electron interaction.

Details about the GW implementation in ABINIT can be found [[theory:mbt|here]]

A typical GW calculation consists of two different steps (following a DFT
calculation): first the screened interaction ε-1 is calculated and stored on
disk ([[optdriver]]=3), then the KS band structure and W are used to evaluate
the matrix elements of Σ, finally obtaining the QP corrections
([[optdriver]]=4).

The computation of the screened interaction is described in
[[topic:Susceptibility]], while the computation of the self-energy is
described in [[topic:SelfEnergy]]. The frequency meshes, used e.g. for
integration along the real and imaginary axes are described in
[[topic:FrequencyMeshMBPT]].



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:gw1]] The first tutorial on GW (GW1) deals with the computation of the quasi-particle band gap of Silicon (semiconductor), in the GW approximation (much better than the Kohn-Sham LDA band structure), with a plasmon-pole model. 
* [[tutorial:gw2]] The second tutorial on GW (GW2) deals with the computation of the quasi-particle band structure of Aluminum, in the GW approximation (so, much better than the Kohn-Sham LDA band structure) without using the plasmon-pole model. 
* [[tutorial:paral_mbt|Parallelism of Many-Body Perturbation calculations (GW)]] allows to speed up the calculation of accurate electronic structures (quasi-particle band structure, including many-body effects).

