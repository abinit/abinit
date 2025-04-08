---
description: How to perform Real-Time TDDFT calculations
authors: FBrieuc
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives some information on how to perform a real-time time-dependent DFT (RT-TDDFT) calculation with the ABINIT package.

!!! warning

    RT-TDDFT is under active development and should thus be used with caution!

## Introduction
The goal of TDDFT is usually to describe the electronic response to an external time-dependent perturbation.
To do so, real-time TDDFT numerically integrates the time-dependent Kohn-Sham (TDKS) equations in _real-time_ and
thus gives access to the time evolution of the electronic density directly.
Similarly to the widely used linear-response TDDFT approach (see [[topic:TDDFT]]), this method can be used
to compute the electronic response in the linear regime giving access among others to transport coefficients,
such as the electrical conductivity, and optical properties, such as the absorption spectrum.
Moreover, RT-TDDFT is not restricted to the linear regime, and can thus be used to study the response to
intense perturbations.
This method is thus particularly suited to investigate the non-equilibrium electron dynamics following an
intense excitation such as the response to high intensity lasers.
It has been successfully used to study different phenomena including high-harmonics generation,
electron stopping power, core electron excitations etc.
A detailed description of TDDFT including real-time propagation schemes can be found for instance
in the book of C. Ullrich [[cite:Ullrich2011]].

## Implementation in ABINIT
ABINIT implements RT-TDDFT in the so-called adiabatic approximation using the standard XC functionals
developed for ground state calculations. The implementation works with LDA and GGA functionals.
*It has not yet been tested for other types (meta-GGAs, hybrids) and is thus most probably not compatible for now.*
*Moreover, the implementation has not been yet tested on all possible cases such as magnetic cases or
with spin-orbit coupling and thus should not be used in theses cases for now.*
The integration of the TDKS equations works with both norm-conserving pseudopotentials and
the projector augmented wave (PAW) method [[topic:PseudosPAW]].

As of now only the simple Exponential Rule (ER) and the Exponential Mid-point Rule (EMR) propagators are
implemented which seem to be sufficiently stable in most cases (see for instance [[cite:Castro2004]]
for more information on propagators for RT-TDDFT). It is possible to apply an external impulse electric
field in order to compute the associated response functions, so typically the conductivity and the
dielectric function.
The [[tutorial:rttddft|Tutorial on real-time TDDFT]] describes how to run such calculations to compute the
dielectric function of Diamond.
*Note that the application of an external electric field is only possible in PAW for now.*

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* See [[tutorial:rttddft|Tutorial on real-time TDDFT]], to learn more on using real-time TDDFT with ABINIT.
