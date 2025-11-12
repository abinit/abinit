---
description: How to take into account an external magnetic field
authors: EB, JWZ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to take into account an external magnetic field with the ABINIT package.

## Introduction

An applied external magnetic field has several types of interactions
with a system of electrons and nuclei.  In particular, it couples to
the electronic spins (Zeeman term) as well as to the electronic
orbital motion (orbital term).  Input variables to deal with these two
situations are decoupled in ABINIT.  [[zeemanfield]] gives access to
the Zeeman coupling.  The coupling to the electronic orbital motion
requires a more elaborate calculation, and is triggered by the input
variable [[orbmag]]. It is implemented following the theory outlined
in [[cite:Zwanziger2023]]. For insulators, this calculation would
typically be used in conjunction with an imposed nuclear magnetic
dipole, see [[nucdipmom]], in order to compute the chemical shielding
as observed in NMR.

Input variables for both types of calculation are listed below, as
well as relevant test files and the [[tutorial:nuc|tutorial on
properties at the nucleus]].

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:nuc|The tutorial on properties at the nucleus]] describes the use of [[orbmag]]
and [[nucdipmom]] to compute the chemical shielding at a nuclear site.

