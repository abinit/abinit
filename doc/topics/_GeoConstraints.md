---
description: How to constrain the geometry of the system in geometry optimization, molecular dynamics or searches
authors: GG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to constaint the geometry of the system in geometry optimization, molecular
dynamics or searches with the ABINIT package.

## Introduction

There are two mechanisms to put constraints on the atom positions in ABINIT.
They can be used in [[topic:GeoOpt|geometry optimization]],
[[topic:MolecularDynamics|molecular dynamics]] (including PIMD) or other
geometry algorithms (e.g. [[topic:TransPath|transition path searches]]).

The simplest one (entry point [[iatfix]]) simply define a set of atoms that
are fixed, either entirely, or only along one of the directions (for the
latter, see the warning in [[iatfix]])

A more complex one, but also much more powerful, allows to place constraints
on linear combinations of atomic positions. Thanks to such constraint, the
mean position of two atoms (or a fragment, like a molecule) can be fixed, or
constrained to stay within an arbitrary plane. One can thus also sample
different mean positions. See a complete description in [[wtatcon]].



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

