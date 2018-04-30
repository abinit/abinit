---
description: How to calculate the RPA correlation energy
authors: FB
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to calculate the RPA correlation energy with the ABINIT package.

## Introduction

In the adiabatic-connection fluctuation-dissipation framework, the correlation
energy of an electronic system can be related to the density-density
correlation function, also known as the reducible polarizability. When further
neglecting the exchange-correlation contribution to the polarizability, one
obtains the celebrated random-phase approximation (RPA) correlation energy.
This expression for the correlation energy can alternatively be derived from
many-body perturbation theory. In this context, the RPA correlation energy
corresponds to the GW total energy.

The RPA correlation energy can be expressed as an integral function of the
dielectric matrix (see [[cite:Gonze2016]]). The integral over the frequencies
is performed along the imaginary axis, where the integrand function is very
smooth. Only a few sampling frequencies are then necessary. In ABINIT, the RPA
correlation energy is triggered by setting the keyword [[gwrpacorr]] to 1.

The RPA correlation energy is a post-processed quantity from the GW module of
ABINIT, which takes care of evaluating the dielectric matrix for several
imaginary frequencies.

The RPA correlation has been shown to capture the weak van der Waals
interactions [[cite:Lebegue2010]] and to drastically improve defect formation
energies [[cite:Bruneval2012]].

The convergence versus empty states and energy cutoff is generally very slow.

It requires a careful convergence study. The situation can be improved with
the use of an extrapolation scheme ([[cite:Bruneval2008]], [[cite:Harl2010]]).



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

