---
authors: FJ
description: TDDFT Abinit topic
---

This page gives hints on how to perform time-dependent density-functional theory calculations of neutral
excitation energies with the ABINIT package.

## Introduction

For finite systems (atoms and molecules), excited states can be computed
within TDDFT (Casida approach - only norm-conserving pseudopotentials). See
the explanations given in the [[lesson:tddft]] lesson of tutorial. The
[[iscf]] input variable must be set to -1.

In the non-spin-polarized case, spin-singlet as well as spin-triplet
excitations are computed. Spin-polarized case is also available.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[lesson:tddft|The lesson on TDDFT]] deals with the computation of the excitation spectrum of finite systems, thanks to the Time-Dependent Density Functional Theory approach, in the Cassida formalism.

