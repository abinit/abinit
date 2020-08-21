---
description: How to analyze the densities and potentials
authors: SS, XG, YG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to analyze the densities and potentials with the ABINIT package.

## Introduction

All the files that have the density and potential format, see
[[topic:printing]] can be analyzed with the "Cut3D" postprocessor. In
particular, it can produce two-dimensional cuts (or one-dimensional cuts)
through the three-dimensional data, suitable for later visualisation using
e.g. [[topic:Abipy]]. It can perform the Hirshfeld computation of atomic
charges. It can analyse the charge contained in an atomic sphere, and
determine the angular momentum projected charge (l=0 to 4) contained in that
sphere. (only available for norm-conserving pseudopotentials)

See the [[help:cut3d]], as well as the [[tutorial:cut3d]].

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:cut3d]] explains the use and input parameters needed for the "Cut 3-Dimensional files" post-processor of the ABINIT package

