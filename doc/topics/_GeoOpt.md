---
description: How to perform a geometry optimization
authors: GG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform a geometry optimization with the ABINIT package.

## Introduction

Different algorithms (Broyden; modified Broyden; Verlet with sudden stop of
atoms, FIRE) allows to find the equilibrium configuration of the nuclei, for which
the forces vanish, see [[ionmov]]. The cell parameters can also be optimized
concurently with the atomic positions [[optcell]], possibly with a state of
given stress as a target, [[strtarget]].

Specified lattice parameters, or angles, or atomic positions, can be kept
fixed if needed, see [[topic:GeoConstraints]].

A genetic algorithm has been coded, for global optimisation. See [[ga_rules]]. Not in production yet. 

The results of the structural relaxation are saved in the HIST.nc file.
This file can be used to restart the calculation or to analyze the results at the end of calculation.
For further details about the |AbiPy| API please consult the |HistFileNb|.

When performing cell parameter optimization, the default setting of ABINIT will not permit increase
of the lattice parameters from the initial ones. One needs to set [[dilatmx]] to some value larger
than the default 1.0 . If the starting lattice parameters are thought to be accurate and should
not increase by more than 5%, set [[dilatmx]] to 1.05 , which will lead to a moderate increase
of CPU time by about 15% . If the starting lattice parameters are just a guess, perform a first run
with [[chkdilatmx]]=0 . This will yield an estimate of the correct lattice parameters, but such estimate are 
not the numerically exact ones for the [[ecut]] value. 
You are done if you do not need accurate geometry estimation anyhow. Otherwise,
after this first run, perform a second run by restarting from the newly estimated value, 
with [[chkdilatmx]]=1 and [[dilatmx]] to 1.05 .


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

