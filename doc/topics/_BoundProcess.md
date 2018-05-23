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

## Example

This is an example of how to bound a model with multibinit:

* General flags:

        bound_model    = 2
        bound_maxCoeff = 10  

* Flags for the generation of the coefficients:
  
        bound_rangePower = 6 8
        bound_cutoff     = 10 
        bound_anhaStrain = 0
        bound_SPCoupling = 1
        bound_anhaStrain = 0
  
* Flags for the molecular dynamics during the bounding process:

        bound_cell  = 10 10 10
        bound_step  = 5000  
        bound_temp  = 500
  
The previous flags will activate the option two of the bound process ([[multibinit:bound_model]]=2) and try to add up to 10 anharmonic coefficients ([[multibinit:bound_maxCoeff]]=10). The list of all the possible coefficients will be generated with a range of power from 3 to 4 ([[multibinit:bound_rangePower]]=6 8) and a cut off of 10 bohr ([[multibinit:bound_cutoff]]=10). Moreover, the strain-phonon counpling is activated ([[multibinit:bound_SPCoupling]]=1) but you will not fit the anharmonic strain ([[multibinit:bound_anhaStrain]]=0). This option of the bounding process will add one by one the additional coefficients and test if the model is bound by running a molecular dynamics on a 10x10x10 supercell  ([[multibinit:bound_cell]]=10 10 10) at 500K  ([[multibinit:bound_temp]]=500)  over 5000 time steps ([[multibinit:bound_step]]=5000).

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

