---
description: How to fit the anharmonic part of a model in multibinit
authors: AM,ACGC
---

This page gives hints on how to fit the anharmonic part of a model in multibinit.

## Introduction

The fit process implemented in multibinit is based on [[cite:Escorihuela-Sayalero2017]].
The fit process of multibinit contains two main options:

* Generate the list of symmetry-adapted terms
* Select the best coefficients in the list with the fit process
  
In the first option, multibinit will generate a set of coefficients only if [[multibinit:fit_generateCoeff]] is set to one. This generation is mainly parametrized by [[multibinit:fit_rangePower]] and [[multibinit:fit_cutoff]]. You can avoid the generation by providing a list of coefficients with the model_anharmonic.XML file (see [[help:multibinit]]).


Then, the fit process will select the coefficients one by one up to [[multibinit:fit_ncoeff]] according to the procedure details in [[cite:Escorihuela-Sayalero2017]]. This process requires to provide the training_set_HIST.nc list file (see [[help:multibinit]])
  

## Tutorials
The [[tutorial:lattice_model|First lesson on Multibinit]] explains how to build a lattice model and to run a dynamics.
    
## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

