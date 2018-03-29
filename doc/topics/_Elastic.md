---
description: How to compute elastic, piezoelectric and internal strain tensors from DFPT
authors: MT
---

This page gives hints on how to compute elastic, piezoelectric and internal strain tensors 
from DFPT with the ABINIT package.

## Introduction

The DFPT theory to compute the elastic, piezoelectric and internal strain
tensors is presented in [[cite:Hamann2005]].

See the [[topic:DFPT|generalities about DFPT]] as well as the tutorial
[[lesson:elastic]]. For the preliminary runs of ABINIT, see [[rfstrs]], while
for ANADDB, see [[anaddb:elaflag]] for the elastic tensor, and
[[anaddb:piezoflag]] for the piezoelectric tensor.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[lesson:elastic|The lesson on the elastic properties]] presents the computation with respect to the strain perturbation and its responses: elastic constants, piezoelectricity.
* [[lesson:paral_dfpt|Parallelism of response-function calculations]]. Additional information to use the DFPT in parallel.

