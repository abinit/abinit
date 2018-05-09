---
description: How to fit the anharmonic part of a model in multibinit
authors: AM
---

This page gives hints on how to fit the anharmonic part of a model in multibinit.

## Introduction

The fit process implemented in multibinit is based on this [[cite:Escorihuela-Sayalero2017|paper]].
The fit process of multibinit contains two main options:

* Generate the list of symetried adapted term
* Select the best coefficients in the list with the fit process
  
In the first option, multibinit will generate a set of coefficients only if [[multibinit:fit_generateCoeff]] is set to one. This generation is mainly parametrized by [[multibinit:fit_rangePower]] and [[multibinit:fit_cutoff]]. By the way, you can also provided a list of coefficients with the model_anharmonic.MXL (see [[help:multibinit]]).


Then the fit process will select the coefficients one by one up to [[multibinit:fit_ncoeff]] according to the procedure details in the [[cite:Escorihuela-Sayalero2017|paper]]. This process requires to providea training_set_HIST.nc file (see [[help:multibinit]])
  
## Example

This is an example of how to fit a model with multibinit:

* General flags:

        fit_coeff  = 1     
        fit_ncoeff = 10  

* Flag for the generation of the coefficients:
  
        fit_generateCoeff  = 1
        fit_rangePower     = 3 4 
        fit_cutoff         = 10 
        fit_SPCoupling     = 1 
        fit_anhaStrain     = 0


The previous flags will activate the fit process ([[multibinit:fit_coeff]]=1) and fit 10 anharmonic coefficients ([[multibinit:fit_ncoeff]]=10). The list of all the possible coefficients will be generated ([[multibinit:fit_generateCoeff]]=1) with a power range of 3 to 4 ([[multibinit:fit_rangePower]]=3 4) and a cut off of 10 bohr ([[multibinit:fit_cutoff]]=10). Moreover, the strain-phonon counpling is activated ([[multibinit:fit_SPCoupling]]=1) but you will not fit anharmonic strain ([[multibinit:fit_anhaStrain]]=0).

### Additional flags for the convergence of the fit process:

  These next flags can be used for adding a criteria to stop the fit process. For example if you set [[multibinit:fit_ncoeff]] to 100, the fit process can be stop before to reach the 100 coefficients if one of these tolerance flags:

* Tolerance on the Mean Standard Deviation of the Energy in (meV/atm):

        fit_tolMSDE  = 1E-01
  
* Tolerance on the Mean Standard Deviation of the Stresses in (eV^2/A^2):

        fit_tolMSDS  = 3E-06

* Tolerance on the Mean Standard Deviation of the Forces (eV^2/A^2):

        fit_tolMSDF  = 3E-06

* Tolerance on the Mean Standard Deviation of the Forces+Stresses (eV^2/A^2):

        fit_tolMSDFS = 3E-06
  

### Advanced user flags:

It is possible to disable the storage in memory of the datas for the fit, thus you allow multibinit to generate more possible coefficients but the computation time will increase. To disable the storage of the data, add this next flag:

          fit_initializeData = 1
  
The next flags can be used to explicitly fix or ban some coefficients during the fit process.
These flags can be very useful in many applications.

Here an example of using these flags:

* After a run, the fit process can selected these specific:

          12 => (Sr-Ti)^3(eta_1)^1
          452 => (Sr-O1)^3(eta_2)^1
          6543 => (Sr-O2)^3(eta_3)^1
          2434 => (Ti-O1)^3(eta_4)^1
          5234 => (Ti-O2)^3(eta_5)^1

* So now you want to add more coefficients to you model but not restart from scratch, you can add in the input:

          fit_nfixcoeff = 5
          fit_fixcoeff  = 12 452 6543 2434 5234

* Moreover, if you want to exclude a specific coefficient you can add:
  
          fit_nbancoeff = 1 
          fit_bancoeff  = 5234



***WARNING*** Sometimes, depending on the code used to generate the training set, the stresses can be expresse with a -1 factor (not the case in ABINIT). So you need to add in the input:

          fit_ts_option = 1

    
## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

