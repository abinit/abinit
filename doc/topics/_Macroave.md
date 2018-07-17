---
description: How to perform macroscopic averages of the densities and potentials
authors: XG
---
<!--- This is the source file for this topics. Can be edited. -->

This page explains how to perform macroscopic averages of the densities and potentials with the ABINIT package.

## Introduction

The MACROAVE program implements the macroscopic average technique,
introduced by A. Baldereschi et al [[cite:Baldereschi1988]].
This powerful method relates microscopic quantities, typical outputs of first-principles codes,
with macroscopic magnitudes, needed to perform electrostatic analysis.
Within this methodology, one is able to wash out all the
wiggles of the rapidly-varying functions of position (resembling the underlying atomic structure) 
of the microscopic quantities,
blowing up only the macroscopic features.
It can be used to compute band offsets, work functions, effective
charges, and high frequency dielectric constants, among others.

See also [[cite:Colombo1991]] as well as [[cite:Shaltaf2008]], the latter for an application using 
the MACROAVE program delivered with ABINIT, see Fig.1 .

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

