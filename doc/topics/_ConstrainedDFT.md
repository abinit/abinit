---
description: How to perform calculation within constrained DFT
authors: EB and XG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform calculation with constrained DFT (atomic charge, atomic magnetic moments)
with the ABINIT package.

## Introduction

Constrained Density Functional Theory (cDFT) can be used to impose constraints on the atomic charges,
and on the atomic magnetic moments.
Usually integrals of the charge density and magnetization (real-space functions) inside spheres
are constrained to user-defined values.
This is described in e.g. [[cite:Kaduk2012]] or [[cite:Ma2015]].

Two algorithms are implemented in ABINIT.
The most recent one, see [[cite:Gonze2022]], is an improvement with respect to the algorithms reported
in both [[cite:Kaduk2012]] or [[cite:Ma2015]] papers.
It is activated thanks to the [[constraint_kind]] input variable.
The algorithm in [[cite:Kaduk2012]], initially reported in [[cite:Wu2005]], implements a double-loop cycle,
which is avoided in the present implementation.
The algorithm presented in [[cite:Ma2015]] is based on a penalty function (so this algorithm is NOT a Lagrange multiplier approach, unlike claimed by these authors).
The latter is also implemented in ABINIT, but only for the atomic magnetic moments, see [[topic:MagMom]].
The recent algorithm is also superior to the latter, as the recent algorithm enforces the constraint
to arbitrary numerical precision at convergence,
instead of imposing it approximately with a tunable accuracy under the control of [[magcon_lambda]].
The present algorithm is also not subject to instabilities that have been observed when [[magcon_lambda]] becomes larger and larger
in the [[cite:Wu2005]] algorithm.
Still, for GGA functionals, the recent algorithms is not unconditionally stable, and this is currently under analysis.

In the recent algorithm, ABINIT implements forces as well as stresses in cDFT.
Also, derivatives of the total energy with respect to the constraint are delivered.
For four types of contraints (the atomic charge constraint, the vector atomic magnetization constraint,
the atomic magnetization length constraint and the atomic magnetization axis constraint),
the derivative is computed with respect to the value of the constraint defined directly by the user,
while for the magnetization direction constraint, the derivative is evaluated with respect to the change of angle.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}
