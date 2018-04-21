---
description: How compute the frequency-dependent susceptibility matrix, and related screened interaction matrix, and inverse dielectric marix
authors: MG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to compute the frequency-dependent susceptibility matrix, and related screened
interaction matrix, and inverse dielectric marix with the ABINIT package.

## Introduction

In the independent-particle approximation, the frequency-dependent
susceptibility matrix, and related screened interaction matrix, and inverse
dielectric matrix can be computed.

This can be done on top of eigenfunctions and eigenvalues obtained from Kohn-
Sham, generalized Kohn-Sham (e.g. hybrid functionals), as well as self-
consistent quasiparticle methodology in the (generalized) Kohn-Sham basis.

This is a prerequisite to many-body perturbation theory calculations, see
[[topic:GW]] and [[topic:BSE]], to which we refer.

The frequency meshes, used e.g. for integration along the real and imaginary
axes, on which the susceptibility matrices (and related matrices) have to be
computed are described in [[topic:FrequencyMeshMBPT]].



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

