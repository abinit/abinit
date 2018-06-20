---
description: How to compute the electronic self-energy (due to electron-electron interaction)
authors: MG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to compute the electronic self-energy (due to electron-electron interaction) with the ABINIT package.

## Introduction

In principle, the exact self-energy can be obtained by solving self-consistently 
the set of coupled integro-differential equations proposed by
Hedin [[cite:Hedin1965]]. The fundamental building blocks of Hedin's equations
are, besides $\Sigma(1,2)$, the Green's function of the interacting many-body system,
$G(1,2)$, the Green's function of an appropriate non-interacting system, $\Go(1,2)$,
and the irreducible polarizability, $\tchi(1,2)$, which, through the inverse
dielectric matrix $\ee^{-1}(1,2)$, re-normalizes the static Coulomb potential,
resulting in the dynamical screened interaction $W(1,2)$. Finally, the vertex
function $\Gamma(1,2,3)$ describes the interactions between virtual holes and
electrons.

A typical self-energy calculation combines a quasi-particle band structure
with a screened interaction and possibly a vertex correction to the QP
corrections ([[optdriver]]=4).

In the frequency domain, the GW self-energy $\Sigma(\omega)$ can be evaluated in ABINIT
with two different, more effective, techniques:

* integration with a plasmon-pole model (PPM)
* integration with contour deformation (CD).

In the former case, the frequency dependence of $\ee^{-1}(\omega)$, is modeled with a
simple analytic form, and the frequency convolution is carried out analytically.
In the latter approach, the integral is evaluated numerically extending the
functions in the complex plane in order have a smoother integrand.

Four different plasmon pole models (PPMs) are available in ABINIT. The choice
of the particular PPM to be used is controlled by the variable [[ppmodel]].
The first two options ([[ppmodel]] = 1, 2) refer to approximations employed in
the pioneering implementations of the GW formalism: the plasmon-pole models of
Godby-Needs [[cite:Godby1989]] (GN) and Hybertsen and Louie [[cite:Hybertsen1986]] (HL).

The contour deformation technique is activated by setting the input variable
[[gwcalctyp]] to 2. The integration along the imaginary axis requires the
calculation of $\ee^{-1}(\omega)$, for purely imaginary frequencies. 
The frequency mesh for the quadrature is governed by the input variable [[nfreqim]], and can be very
coarse since the integrands is very smooth in this region.

The evaluation of the residue of the poles requires the calculation of $\ee^{-1}(\omega)$
on a fine mesh along the real axis. This regular mesh, sampling the interval
[0, +∞], is defined by the two input variables [[nfreqre]] and [[freqremax]].

The CD approach requires many evaluations of $\ee^{-1}(\omega)$ and can therefore be
computationally highly demanding. On the other hand, it is the preferred
approach for calculating the QP correction of low-lying states. Moreover, it
is the only technique available in ABINIT to compute the imaginary part of
$\Sigma(\omega)$ and the spectral function $A(\omega)$.

It is possible to disable the full computation, and actually do an Hartree-Fock, 
screened exchange, COHSEX or hybrid functional calculation.
The calculation is done in a precomputed basis set, that can be Kohn-Sham
(e.g. PBE) or generalized Kohn-Sham (e.g. HSE06).

As vertex corrections, the bootstrap kernel and others can be included in the
self-consistent W. 
The Faleev method ([[cite:Faleev2004]]), is implemented.

Convergence over the number of unoccupied band is much improved with respect
to usual implementations of GW, thanks to the "extrapolar" method.

The frequency meshes, used e.g. for integration along the real and imaginary
axes are described in [[topic:FrequencyMeshMBPT]].


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* The first tutorial on GW ([[tutorial:gw1|GW1]]) deals with the computation of the quasi-particle band gap of Silicon (semiconductor), in the GW approximation (much better than the Kohn-Sham LDA band structure), with a plasmon-pole model. 
* The second tutorial on GW ([[tutorial:gw1|GW2]]) deals with the computation of the quasi-particle band structure of Aluminum, in the GW approximation (so, much better than the Kohn-Sham LDA band structure) without using the plasmon-pole model. 
* [[tutorial:paral_mbt|The tutorial on Parallelism of Many-Body Perturbation calculations (GW)]] allows to speed up the calculation of accurate electronic structures (quasi-particle band structure, including many-body effects).

