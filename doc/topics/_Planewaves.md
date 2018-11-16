---
description: How to perform numerically precise calculations with planewaves or projector- augmented waves and pseudopotentials
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform numerically precise calculations with planewaves or projector-
augmented waves and pseudopotentials with the ABINIT package.

## Introduction

The numerical precision of the calculations depends on many settings, among
which the definition of a basis set is likely the most important. With
planewaves, there is one single parameter, [[ecut]] that governs the
completeness of the basis set.

The wavefunction, density, potentials are represented in both reciprocal space
(plane waves) and real space, on a homogeneous grid of points. The
transformation from reciprocal space to real space and vice-versa is made
thanks to the Fast Fourier Transform (FFT) algorithm. With norm-conserving
pseudopotential, [[ecut]] is also the main parameter to define the real space
FFT grid, In PAW, the sampling for such quantities is governed by a more
independent variable, [[pawecutdg]]. More precise tuning might be done by
using [[boxcutmin]] and [[ngfft]].

Avoiding discontinuity issues with changing the size of the planewave basis
set is made possible thanks to [[ecutsm]].

The [[accuracy]] variable enables to tune the accuracy of a calculation by
setting automatically up to seventeen variables.

Many more parameters govern a PAW computation than a norm-conserving
pseudopotential calculation. They are described in a specific page
[[topic:PAW]]. For the settings related to wavelets, see [[topic:Wavelets]].



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:base2|The tutorial 2]] deals again with the H2 molecule: convergence studies, LDA versus GGA 
* [[tutorial:base3|The tutorial 3]] deals with crystalline silicon (an insulator): the definition of a k-point grid, the smearing of the cut-off energy, the computation of a band structure, and again, convergence studies ...
* The first tutorial on the [[tutorial:paw1|the projector-augmented wave]] technique.

