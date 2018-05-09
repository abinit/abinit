---
description: How to bound a model in multibinit
authors: AM
---

This page gives hints on how to fit bound a model in multibinit.

## Introduction

When an anharmonic part is fitted with multibinit, the resulting model is not bound.

The bound process can be activated with two options:
  
The first option will generate all the possible combinaisons of coefficients from 1 to [[multibinit:bound_maxCoeff]]. Some constrains are imposed during the generation and the fit of the coefficients, they have to be positive and with even power. Finaly, the code will try all the possible combinaisons and try to find a bounded model.


The second option will generate a set of coefficients with a power range defined by [[multibinit:bound_rangePower]] and keep only the coefficients with even power. Then the procedure is similar to the fit process with the constrains to only keep positive coefficients. The bound process will select the coefficients one by one up to [[multibinit:bound_maxCoeff]] and try if the model is bound at each step of the process. 

## Example

This is an example of how to bound a model with multibinit:

* General flags:

        bound_coeff  = 2
        bound_maxCoeff = 10  

* Flag for the generation of the coefficients:
  
        bound_rangePower     = 6 8
        bound_cutoff         = 10 
        bound_anhaStrain     = 0

* Flag for the molecular dynamics during the bound process:

        bound_cell  = 10 10 10
        bound_step  = 5000  
        bound_temp  = 500
  
The previous flags will activate the option two of the bound process ([[multibinit:bound_coeff]]=2) and try to add up to 10 anharmonic coefficients ([[multibinit:bound_maxCoeff]]=10). The list of all the possible coefficients will be generated with a power range of 3 to 4 ([[multibinit:bound_rangePower]]=6 8) and a cut off of 10 bohr ([[multibinit:bound_cutoff]]=10). This option of the bound process will add one by one the additional coefficients and test if the model is bound by running a molecular dynamics on 10 10 10 supercell  ([[multibinit:bound_cell]]=10 10 10) at 500K  ([[multibinit:bound_temp]]=500)  over 5000 time step ([[multibinit:bound_step]]=5000).

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

