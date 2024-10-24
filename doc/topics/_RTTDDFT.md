---
description: How to perform Real-Time TDDFT calculations
authors: FBrieuc
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform a real-time time-dependent DFT (RT-TDDFT) calculation with the ABINIT package.  

## Introduction

The main goal of TDDFT is to compute the electronic response to an external time-dependent perturbation.
In order to do so, real-time TDDFT numerically integrates the time-dependent Kohn-Sham equations in _real-time_ and thus gives access to the time evolution of the electronic density directly.
Similarly to the widely used linear-response TDDFT approach (see [[topics:tddft]]), this method can be used to compute the electronic response in the linear regime giving access among other things to transport coefficients, such as the electrical conductivity, and optical properties, such as the absorption spectrum.
Moreover, RT-TDDFT is not restricted to the linear-response regime, and can thus be used to study the response to intense and/or rapidly varying perturbations. 
This method is thus particularly suited to investigate the non-equilibrium electron dynamics following an intense excitation such as the response to 
high intensity lasers.
RT-TDDFT has been successfully used to study a wide range of phenomena including high-harmonics generation, electron stopping power, core electron excitations etc. [[cite:]]
A detailed description of TDDFT including real-time propagation schemes can be found for instance in the book of C. Ullrich [[cite:Ullrich2012]].

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* See [[tutorial:rttddft|Tutorial on real-time TDDFT]], to learn more on using real-time TDDFT with ABINIT.
