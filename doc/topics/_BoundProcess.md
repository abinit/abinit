---
description: How to bound a model in multibinit
authors: AM,ACGC
---

This page gives hints on how to bound a model in multibinit.

## Introduction

When a model is build, it can contains phonons instabilities and may produce a divergence during the dynamics.
To avoid such behaviour,  multibinit includes an automatic procedure to bound the model.

Multibinit will first generate a set of coefficients with a power range defined by [[multibinit:bound_rangePower]] and keep only the coefficients with even power.
Then, the procedure is similar to the fit process with the constrains to only keep positive coefficients. Next, the bounding process will add the new coefficients to the model up to [[multibinit:bound_maxCoeff]] checking if the model is bound in each step of the process.

## Tutorials
The [[tutorial:lattice_model|First lesson on Multibinit]] explains how to build a lattice model and to run a dynamics.
  
## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

