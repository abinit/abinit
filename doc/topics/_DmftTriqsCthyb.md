---
description: How to perform a DFT+DMFT calculation with TRIQS/CT-HYB
authors: BAmadon
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives some information on how to perform a DFT+DMFT calculation with
the interface between ABINIT and TRIQS/CT-HYB.

## Introduction

DFT+DMFT is a method made to better describe correlated systems - like
$3d$ transition metals or $4f$ rare earths - where some electrons stay pretty
localized and show strong correlation effects.

Unlike regular DFT, which is based on the electronic density, DFT+DMFT works with
the Green's function. This means it can describe one-particle excitations and gives you
access to the spectral function directly.

In short, DFT+DMFT maps the complicated many-body problem of electron-electron interactions
onto a simpler Anderson impurity model. In this model, the "impurity" represents the correlated
orbitals on a single lattice site, while the self-consistent "bath" represents all the other
electrons - either on different sites or in other orbitals. The local interactions between the
impurity electrons are treated exactly, while the bath's influence is captured through a
dynamical mean field.

So, running a DFT+DMFT calculation means solving this impurity model, which is handled by an
impurity solver - the most computationally demanding part of the process. ABINIT already includes
its own solvers, but they use the density-density approximation, while TRIQS/CT-HYB is numerically
exact.

We have built an interface between ABINIT and TRIQS/CT-HYB, so you can run fully charge self-consistent
DFT+DMFT calculations on real materials. In this setup, TRIQS/CT-HYB acts as an external library that
solves the impurity problem, while ABINIT handles everything else. You can activate this interface by
setting [[dmft_solv]] to 6 or 7.

It is worth noting that this interface does not simply differ from ABINIT's internal DMFT in the
impurity solver - there are also important changes in the formulas, algorithms, and even some
input variables used in the self-consistent loop.

Our interface was specifically designed to compute the Baym-Kadanoff functional and implements a
stationary version of DFT+DMFT. It includes features like the exact double counting formula and
an analytical evaluation of the high-frequency moments of the Green's function - all of which are
not available in ABINIT's internal DMFT implementation.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:dmft_triqs|The tutorial on DFT+DMFT with TRIQS/CT-HYB]] shows how to use the interface
  and perform a self-consistent DFT+DMFT calculation on Fe.

