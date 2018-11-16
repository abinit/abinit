---
description: How to compute vibrational free energy, entropy, specific heat, thermal expansion, as well as atomic temperature factors
authors: XG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to compute vibrational free energy, entropy, specific heat, thermal expansion, as
well as atomic temperature factors with the ABINIT package.

## Introduction

When the phonon band structure and corresponding eigenvectors are known over
the whole Brillouin Zone, thanks fo Fourier interpolation (see
[[topic:PhononBands]]), integrals can be performed, allowing to obtain a
wealth of properties, like free energy, entropy, specific heat, as well as
atomic temperature factors. For applications of this technique, see
[[cite:Lee1995]].

Moreover, knowing such information for different volumes allows one to compute
the thermal expansion, see [[anaddb:gruns_ddbs]].

The input variables needed to perform the interpolation over the Brillouin
Zone are described in [[topic:PhononBands]] and are not listed again in the
present topic.


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:rf2|The tutorial Response-Function 2 (RF2)]] presents the analysis of the DDBs that have been introduced in the preceeding tutorial RF1. The computation of the interatomic forces and the computation of thermodynamical properties is an outcome of this tutorial.

