---
description: How to set parameters related to the electronic wavevectors (k-points)
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to set parameters related to the electronic wavevectors (k-points) with the ABINIT package.

## Introduction

Since ABINIT is based on periodic boundary conditions, every wavefunction is
characterized by a wavevector, usually denoted k-point.
Any list of k-point can be specified, thanks to the keywords [[nkpt]] and [[kpt]].

Still, k-points are used in two different contexts in the vast majority of cases:

  * the sampling of the Brillouin Zone, with the goal to produce integrated quantities 
    (e.g. the charge density, the electronic energy, the electronic DOS ...) that are numerically precise;
  * or the specific computation of wavefunctions and eigenenergies e.g. to get an electronic band structure. 

In the first case, please complete the present topic by reading
[[topic:BandOcc]], while in the second case, please read [[topic:ElecBandStructure]].

In the first case, the Brillouin zone must be sampled adequately, with grids
that, in general will be homogeneous distributions of k-points throughout the
Brillouin Zone (e.g. Monkhorst-Pack grids, or their generalisations). 
For such grids, see [[ngkpt]], [[nshiftk]], [[shiftk]] or even the more general [[kptrlatt]]. 
A list of interesting k point sets can be generated automatically, including a measure 
of their accuracy in term of integration within the Brillouin Zone, see [[prtkpt]], [[kptrlen]]. 
For metals, a joint convergence study on [[tsmear]] and the k-point grid is important.

For the definition of a path of k-points, see [[topic:ElecBandStructure]].  

More detailed explanation concerning the convergence with respect to the
k-point sampling. The number of k-points to be used for this sampling, in the
full Brillouin zone, is inversely proportional to the unit cell volume, but
may also vary a lot from system to system. As a rule of thumb, a system with a
large band gap will need few k-points, while metals will need lot of k-points
to produce converged results. For large systems, the inverse scale with
respect to the unit cell volume is unfortunately stopped because at least one
k-point must be used. The effective number of k-points to be used will be
strongly influenced by the symmetries of the system, since only the
irreducible part of the Brillouin zone must be sampled. 
Moreover the time-reversal symmetry (k equivalent to -k) can be used for ground-state
calculations, to reduce sometimes even further the portion of the Brillouin
zone to be sampled. 
The number of k points to be used in a calculation is named [[nkpt]]. 

There is another way to take advantage of the time-reversal
symmetry, in the specific case of k-points that are invariant under k => -k ,
or are sent to another vector distant of the original one by some vector of
the reciprocal lattice. See below for more explanation about the advantages of
using these k-points.  

As a rule of thumb, for homogeneous systems, a reasonable accuracy may be
reached when the product of the number of atoms by the number of k-points in
the full Brillouin zone is on the order of 50 or larger, for wide gap
insulators, on the order of 250 for small gap semiconductors like Si, and
beyond 500 for metals, depending on the value of the input variable [[tsmear]]. 
As soon as there is some vacuum in the system, the product [[natom]] * [[nkpt]] can be
much smaller than this (for an isolated molecule in a sufficiently large
supercell, one k-point is enough).


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:base3|The tutorial 3]] deals with crystalline silicon (an insulator): the definition of a k-point grid, the smearing of the cut-off energy, the computation of a band structure, and again, convergence studies ...

