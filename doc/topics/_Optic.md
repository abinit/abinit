---
description: How to compute linear and non-linear optical properties in the independent-particle approximation
authors: SS, XG, YG 
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to compute linear and non-linear optical properties 
in the independent-particle approximation with the Abinit package.

## Introduction

Optical and non-linear optical properties can be computed with different
levels of approximation.

The simplest (and fastest) approach relies on the independent-particle
approximation (IPA): the electrons are supposed independent of each other when
reacting to the optical perturbation (even if the initial computation of the
band structure includes interactions in a mean-field sense, like with DFT).
This approximation is also referred to as a "Sum-Over-States" approach (SOS).
This neglects all electron-hole interaction (so no excitonic effects), but
might provide meaningful results in many case, sometimes even quantitatively.
A first problem is linked with the erroneous band gap of the material, but
this can be corrected by a scissor approximation, see [[scissor@optic]].

In Abinit one can either work in the IPA (see below), or take into account the
excitonic effects, see [[topic:BSE]].

In the Abinit package, there are two different utilities to compute optical
responses in the independent-particle approximation: [[help:optic]] and
conducti [[cite:Mazevet2010]]. They have been developed independently of each other, and thus
overlap significantly. The first one computes the linear and non-
linear optical properties as a function of the frequency. It provides the
optical dielectric tensor, the second-harmonic generation (SHG) as well as the
optical rectification tensor (or electro-optic tensor) - without the
contribution from the nuclear displacements. For the further inclusion of the
contribution from nuclear displacements, see [[topic:nonlinear]].

The second utility "conducti" has more capabilities but only at the linear level,
providing the electronic conductivity, dielectric tensor, index of refraction,
reflectivity, absorption, the thermal conductivity, and the thermopower
(electron transport, high temperature, Kubo-Greenwood formalism), the real as well
as imaginary parts.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* See [[tutorial:optic|The tutorial on Optic]], the utility that allows to obtain the
frequency dependent linear optical dielectric function and the frequency
dependent second order nonlinear optical susceptibility, in the simple
"Sum-Over-States" approximation.

