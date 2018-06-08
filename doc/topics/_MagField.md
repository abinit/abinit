---
description: How to take into account an external magnetic field
authors: EB
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to take into account an external magnetic field with the ABINIT package.

## Introduction

An applied external magnetic field has several types of interactions with a system of electrons and nuclei.
In particular, it couples to the electronic spins (Zeeman term) as well as to the electronic orbital motion (orbital term).
Input variables to deal with these two situations are decoupled in ABINIT.
[[zeemanfield]] gives access to the Zeeman coupling, and is in production.
The coupling to the electronic orbital motion is much more delicate to treat, and is currently in development.
Chern number calculation is available.
Both types of input variables are listed below.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

