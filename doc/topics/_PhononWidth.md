---
description: How to compute the linewidth (or lifetime) of phonons, due to the electron-phonon interaction
authors: MV
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to compute the linewidth (or lifetime) of phonons, due to the electron-phonon
interaction with the ABINIT package.

## Introduction

This topic concerns metals only.

After generating a GKK file (see [[topic:ElPhonInt]]), the Electron-Phonon
Coupling (EPC) analysis is performed in anaddb, setting [[anaddb:elphflag]]
variable to 1. Most of the procedure is automatic, but can be lengthy if a
large number of k-points is being used. The [[anaddb:nqpath]] and
[[anaddb:qpath]] variables must be set, specifying a path in reciprocal space.
anaddb generates files containing the phonon linewidths (suffixed `_LWD`) and
frequencies ωqj (suffixed `_BST`) along [[anaddb:qpath]]. One can calculate the
nesting function n(q) = ∑kii' δ(εk,i) δ(εk+q,i') by setting [[anaddb:prtnest]]
to 1 (output to `_NEST`).


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* The tutorial on the [[tutorial:eph|electron-phonon interaction]] presents the use of the utility MRGKK and ANADDB to examine the electron-phonon interaction and the subsequent calculation of superconductivity temperature (for bulk systems).

