---
description: How to perform an effective mass calculation
authors: JLaflamme
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform an effective mass calculation with the ABINIT package.

## Introduction

The direct estimation of effective masses from DFT band curvature using
[[topic:DFPT]] has been implemented within the linear response part of ABINIT
[[cite:Laflamme2016]]. This method avoids the use of finite differences to
estimate these masses, which eliminates the associated numerical noise and
convergence study. To compute the effective masses, one has to set the keyword
[[efmas]] to 1 within a calculation of the derivative of ground-state
wavefunctions with respect to wavevector ([[rfelfd]] = 1 or 2). The effective
masses will then be computed for all k-points and bands present in the
calculation. One can optionally specify the range of bands to be treated for
each k-point with the keyword [[efmas_bands]].

An additional feature of the effective mass implementation is the correct
treatment of degenerate bands. Indeed, the concept of effective mass breaks
down at degenerate band extrema since it is no longer possible to describe
band curvature using a tensor [[cite:Luttinger1955]], [[cite:Mecholsky2014]].
However, using the concept of ``transport equivalent effective mass''
[[cite:Mecholsky2014]] and its adaptation to the **k.p** framework, the
implementation is able to provide the user with effective mass tensors which,
while not describing the band curvature, describe accurately the contribution
of the individual bands to transport properties.

The implementation supports both NCPP and PAW schemes.

Spin-polarized systems ([[nspden]] = 2) as well as spinors ([[nspinor]] = 2)
can be treated, although the spin-orbit interaction can only be treated in the
PAW case.

The treatment of degeneracies is limited to the extremal points of the band
structure (which are the most relevant in any case).

By the way, the first derivative of the eigenenergies is also computed and
printed during a d/dk calculation, and corresponds to the electronic velocity.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

