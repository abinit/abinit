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

To use the MACROAVE program, prepare an input file named "macroave.in" . No other name is admitted.
The examples below should be renamed "macroave.in". 
Issue simply "macroave", without any argument. Two files will be generated, whose names
are constructed from the name of the density or potential file in the third line of the input macroave.in file,
let us take t41_DEN as an example.
The first one will have been appended with .PAV , which is an abbreviation for Planar AVerage.
Only the x-y average has been done to generate this file. From t41_DEN, one generates t41_DEN.PAV.
The second one will have been appended with .MAV , which is an abbreviation for Macroscopic AVerage.
The x-y average has been followed by the macroscopic average to generate this file. From t41_DEN, one generates t41_DEN.MAV.
Note that the planar average is always done on the x-y plane, and the macroscopic average in the z direction.

Note that macroave works in reduced coordinates, so it will not encounter problems with non-orthogonal primitive vectors.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

