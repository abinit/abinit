---
description: How to compute phonon bands, density of states, interatomic force constants, sound velocity ...
authors: MT
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to compute phonon bands, density of states, interatomic force constants, sound
velocity ... with the ABINIT package.

## Introduction

The Fourier transformation of the phonon dynamical matrices generates
interatomic force constants in real space, as explained in
[[cite:Gonze1997a]]. Backtransforming to reciprocal space gives the Fourier
interpolation of the initial phonon band structure. After such Fourier
interpolation, the DOS can be produced (see [[cite:Lee1995]]), the phonon
eigenenergies plotted along lines, the slope of the energy versus cristalline
momentum evaluated (to give sound velocity).

The two-phonon sum and difference spectra can also be obtained, see [[anaddb:dossum]].

For the related computation of temperature-dependent properties, see [[topic:Temperature]].


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:rf2|The tutorial Response-Function 2 (RF2)]] presents the analysis of the DDBs that have been introduced in the [[tutorial:rf1]]. The computation of the interatomic forces and the computation of thermodynamical properties is an outcome of this tutorial.

