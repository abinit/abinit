---
description: How to to perform a Δ-SCF calculation of neutral excitations
authors: FJ
---

This page gives hints on how to to perform a Δ-SCF calculation of neutral excitations with the ABINIT package.

## Introduction

Although formally not justified, difference in total energy using constrained
occupation numbers sometimes results in surprisingly good agreement with
experimental results, for neutral excitations, in molecules or doped solids.
See e.g. [[cite:Jia2017]].

The manual specification of the occupation numbers is needed in this case.
This is accomplished with the input variable [[occ]], coupled with
[[occopt]]=0 for the homogeneous occupation of a band throughout the Brillouin
Zone, or with [[occopt]]=2 for the specific occupation of a state for a
selected k-wavevector. For very big cells, both should be equivalent, and
[[occopt]]=0 is easier to use.


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

