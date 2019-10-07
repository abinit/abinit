---
description: How to perform calculation within constrained DFT
authors: EB and XG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform calculation with constrained DFT (atomic charge, atomic magnetic moments) with the ABINIT package.

## Introduction

Constrained Density Functional Theory imposes constraints on the charge density and magnetic moments. Usually
integrals of the charge density and magnetization (real-space functions) inside spheres are constrained to user-defined
values. This is described in e.g. [[cite:Kaduk2012]] or [[cite:Ma2015]].

The algorithm implemented in ABINIT is a clear improvement of the algorithm reported in [[cite:Kaduk2012]], 
initially reported in [[cite:Wu2005]], in that there is no double-loop cycle, but a direct (and thus faster) determined of the density or magnetization under constraint. 
A description of the ABINIT algorithm should be published sooner or later.
It is also an improvement on the algorithm presented in [[cite:Ma2015]], based on a penalty function (it is NOT a Lagrage multiplier approach, unlike claimed by these authors) also implemented in ABINIT,
see [[topic:MagMom]], in that it imposes to arbitrary numerical precision the constraint, instead of an approximate one with a tunable accuracy under the control of
[[magcon_lambda]],
and is not subject to instabilities that have been observed when [[magcon_lambda]] becomes larger and larger.

THIS TOPIC WILL BE FURTHER ENLARGED AFTER PUBLICATION OF THE ABOVE_MENTIONED ARTICLE.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

