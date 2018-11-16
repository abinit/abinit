---
description: How to set parameters related to the exchange and correlation functionals
authors: FJ
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to set parameters related to the exchange and correlation functionals with the ABINIT package.

## Introduction

Total energy computation in ABINIT is done according to Density Functional
Theory (DFT). Although formally exact, an approximate exchange-correlation
(XC) functional must be chosen. This is governed by the input variable
[[ixc]].  
However, the pseudopotentials (or PAW data sets) are constructed for one
specific XC functional. If [[ixc]] is not specified, ABINIT will simply take
the [[ixc]] of the given pseudopoential(s) - hoping they are coherent with
each others. One introduces an error by using a pseudopotential generated with
an XC functional that is not the same as the one explicitly specified by
[[ixc]]. However, ABINIT will nevertheless do the calculation.

Many exchange-correlation functionals are available (see the list in the
description of [[ixc]]), through two different implementations: one is the
native ABINIT implementation, the other is the ETSF library of XC functionals
(LibXC is a plug-in to ABINIT). In the native ABINIT set, most of the
important local approximations (LDA) are available, including the Perdew-
Zunger one. Two different local spin density (LSD) are available, including
the Perdew Wang 92, and one due to M. Teter. The Perdew-Burke-Ernzerhof, the
revPBE, the RPBE and the HCTH GGAs (spin unpolarized as well as polarized) are
also available.  
In the LibXC 2.0 library, as interfaced with ABINIT, there are 24 functional
forms of the 3D LDA type, and over 80 functional forms of the GGA type. They
can be used with norm-conserving pseudopotentials as well as PAW atomic data.
Also, some metaGGA can be used with ABINIT (norm-conserving case only). For
response-function type calculations, the native ABINIT LDA and GGA kernels can
be used as well as the LibXC ones.  

#### **Hybrid functionals:**

  
The exchange can also be computed on the basis of the Fock expression (exact
exchange), and the correlation can be computed on the basis of the RPA
approximation (see the GW section). [[topic:Hybrids|Hybrid functionals]]
calculations (HSE06, PBE0, B3LYP) can be performed. The implementation of the
exact exchange, correlation and hybrid does not deliver the forces and
stresses at present, at the exception of forces with norm-conserving
pseudopotentials.

#### **Local exact exchange:**

When [[useexexch]]=1, the hybrid functional PBE0 is used in PAW, inside PAW
spheres only, and only for correlated orbitals given by [[lexexch]]. To change
the ratio of exact exchange, see also [[exchmix]]. The implementation of local
exact exchange in ABINIT is provided in [[cite:Jollet2009]]. See useful input
variables [[exchmix]], [[lexexch]] and [[useexexch]].  
  

#### **Van der Waals functionals:**

  
[[topic:vdw|Several Van der Waals functionals]] are available: Grimme (D2, D3,
D3(Becke-Johnson)), Silvestrelli.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* The [[tutorial:base2]] deals with the H2 molecule: convergence studies, LDA versus GGA 

