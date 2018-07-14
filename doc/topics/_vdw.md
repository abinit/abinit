---
description: How to use Van der Waals functionals
authors: YPouillon, BVanTroeye
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to use Van der Waals functionals with the ABINIT package.

## Introduction

It is well known that long range correlations responsible of van der Waals
interactions are out of reach for both LDA and GGA approximations to the
exchange-correlation energy in DFT. In recent years several methods have been
devised to include such interactions, which can be grouped into two
strategies, namely _ad hoc_ methods and self-consistent approaches. Currently
ABINIT can perform calculations based on either the DFT-D methods or the vdW-
WF methods, as described later, both belonging to the first group.

A fully customizable implementation of the vdW-DF method [[cite:Dion2004]], a
self-consistent approach, and an adaptation of the strategy followed by
G.Roman-Perez _et al._ [[cite:Romanperez2009]] to the case of ABINIT are under
development. It will offer around 25 ajustable parameters and be delivered
with graphical tools to help users assess the quality of their kernels. It
does not only aim at performing production calculations with vdW-DF, but also
at helping researchers who develop new density functionals optimised for
systems requiring van-der-Waals interactions.

The DFT-D methods have been implemented inside ABINIT, namely DFT-D2
[[cite:Grimme2006]], DFT-D3 [[cite:Grimme2010]] and DFT-D3(BJ)
[[cite:Grimme2011]]. In these cases, pair-wise terms (and 3-body corrections
for DFT-D3 and DFT-D3(BJ)) are added to the DFT energy, which are independent
of the electronic density, in order to mimic the vdW interactions. The
implementation includes the contributions of these methods to forces and
stresses, in view of geometry optimization, as well as to first-order response
functions like dynamical matrices, clamped elastic constants and internal
strain coupling parameters.

To activate DFT-D dispersion correction, two keywords are in use: [[vdw_xc]] =
5/6/7 to choose between DFT-D2, DFT-D3 and DFT-D3(BJ), and [[vdw_tol]], to
control the inclusion of largely distant pairs (those giving a contribution
below [[vdw_tol]] are ignored). It is also possible to include 3-body
corrections [[cite:Grimme2010]] (for ground-state only) with the keyword
[[vdw_tol_3bt]], which also controls the tolerance over this term.

Methods based on maximally localized Wannier functions (MLWFs) to calculate
vdW energy corrections have also been implemented in ABINIT. In this case the
pair-wise terms come from contributions of pairs of MLWFs rather than from
atoms. Among the implemented methods in ABINIT it is found vdW-WF1
[[cite:Silvestrelli2008]], [[cite:Silvestrelli2009]] vdW-WF2
[[cite:Ambrosetti2012]] and vdW-QHO-WF [[cite:Silvestrelli2013]]. A full
description of the implementation of vdW-WF1 is reported in
[[cite:Espejo2012]].

Selection of one of these 3 methods is achieved by using [[vdw_xc]]=10/11/14
respectivelly. Since vdW-WF1 and vdW-WF2 methods are approximations for the
dispersion energy of non overlapping electronic densities, it is necessary to
define the interacting fragments of the system whose dispersion energy is
going to be calculated. The latter is achieved by using the input variables
[[vdw_nfrag]] and [[vdw_typfrag]] to define the number of interacting
fragments in the unit cell and to assign each atom to a fragment. A given MLWF
belongs to the same fragment as its closer atom. The need for defining the
interacting fragments is overridden in the vdW-QHO-WF, for which these input
variables are not used. When dealing with periodic systems the input variable
[[vdw_supercell]] controls the number of neighbor unit cells that will be
included in the calculation. Each one of the 3 components of [[vdw_supercell]]
indicates the maximum number of cells along both positive or negative
directions of the corresponding primitive vector. This is useful for studying
the spacial convergency of the vdW energy. It should be noticed that the user
must set the variables associated to the calculation of MLWFs and that the
resulting vdW energies strongly depend on the obtained Wannier functions.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

