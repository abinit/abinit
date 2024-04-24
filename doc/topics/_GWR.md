---
description: How to perform a GWR calculation 
authors: MG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform GWR calculations

## Introduction

A cubic scaling real-space imaginary-time algorithm for GW and RPA is available.
See the theory in [[cite:Liu2016]] and related references. 
An overview is available in the tutorial [[tutorial:gwr_intro|An overview of the GWR code]].

This implementation relies on the minimax time-frequency grids
available in the GreenX library [[cite:Azizi2023]].
At present, only norm-conserving pseudopotentials can be used. The implementation is restricted to non-magnetic materials,
and without spin-orbit coupling.
Still, different types of flows and algorithms (including self-consistency) are available, see [[gwr_task]].

Activate it using [[optdriver]]=6, and specify [[gwr_task]].

NOTE: GWR code is under active development and not yet ready for production runs.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:gwr_intro|An overview of the GWR code]]. Covers the motivation, requirements, formalism, workflow.
* [[tutorial:gwr1|First tutorial on GWR]]. Still under development at the time of writing. Quasi-particle band structure of silicon in the GW approximation.
