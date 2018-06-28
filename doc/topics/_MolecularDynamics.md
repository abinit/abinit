---
description: How to perform a molecular dynamics calculation
authors: GG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform a molecular dynamics calculation with the ABINIT package.

## Introduction

Three molecular dynamics algorithm (Numerov, Verlet, Blanes and Moanes) allow
to perform simulations in real (simulated) time, see [[ionmov]]. The
displacement of atoms may be computed according to Newton's law, or by adding
a friction force to it. Nose-Hoover thermostat is available with Verlet
algorithm. Langevin dynamics is also available.

Specified lattice parameters, or angles, or atomic positions, can be kept
fixed if needed, see [[topic:GeoConstraints]].

The trajectories can be analyzed thanks to the [[topic:APPA|APPA postprocessor]].


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:paral_moldyn]] Parallelism for molecular dynamics calculations

