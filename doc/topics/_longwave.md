---
description: How to compute spatial dispersion properties with the longwave DFPT approach. 
authors: MR and MS
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to compute spatial dispersion properties 
(e.g. flexoelectric tensor or dynamical quadrupoles) with the longwave DFPT
driver of the ABINIT package.

## Introduction
In condensed-matter physics, spatial dispersion refers to the dependence of many material
properties on the wavevector **q** at which they are probed, or equivalently on the gradients of the
external field (electric, magnetic, strain...) and/or the response in real space. Remarkable
examples of such gradient effects include the natural optical rotation, [[cite:Belinicher1980]] whereby some crystals
rotate the polarization plane of the transmitted light, or the flexoelectric tensor, [[cite:Zubko2013]] which
describes the polarization response to a gradient of applied strain.

Since ABINIT v9.0.2, the calculation of several spatial dispersion quantities is accessible
via the longwave driver. The implementation, detailed in [[cite:Romero2020]], follows the formalism developed in 
[[cite:Royo2019]] that adapts the classic longwave method of Born and Huang [[cite:Born1954]] with the modern tools of 
the DFPT. Technically, the driver computes analytical third-order energy derivatives (readily converted to physical 
spatial dispersion quantities) with respect to two of the *standard* 
perturbations (atomic displacements, electric field and strain) and to the wavevector **q**. 

At present ( |today| ), ABINIT enables the calculation of four spatial dispersion tensors required 
in order to build all the contributions to the bulk flexoelectric tensor following the prescriptions exposed in [[cite:Stengel2013]]. 
These include, the clamped-ion flexoelectric tensor (a purely electronic contribution) and the first real-space moment of three other 
tensors (entering the mixed and lattice mediated contributions): the polarization 
response to an atomic displacement, the interatomic force constants and the piezoelectric force-response
tensor. The implementation also provides access, as a by-product of the flexoelectric formalism, to the dynamical 
quadrupoles, which can be considered as the spatial-dispersion counterpart of the Born effective 
charges (dynamical dipoles). After execution, the longwave routines generate a third-order derivative 
database that is subsequently used by ANADDB either to compute and print the different contributions 
to the flexoelectric tensor, or to consider quadrupolar fields within the Fourier interpolation procedure of 
the dynamical matrix [[cite:Royo2020]].

The underlying theory of the long-wave DFPT approach [[cite:Royo2019]] has been developed for its application on 
time-reversal symmetric insulating crystals only. Therefore, the usage of the longwave driver is restricted to materials
of this kind. An extension of the theoretical framework and the ABINIT implementation to magnetic insulators 
and/or metals will be hopefully pursued in the future. 

Regarding the flexoelectric tensor that ABINIT provides, a few remarks are in order. First, recall that 
this is the **bulk** flexoelectric tensor and that a surface counterpart is still missing in order to obtain 
the total flexoelectric tensor of a system.[[cite:Stengel2016]] Even though the implementation can be applied
to slabs or low-dimensional systems (such as 2D materials), via a supercell approach, the outcome of such a 
calculation will not directly 
produce the total (i.e., bulk+surface) flexoelectric tensor, neither the surface contribution can be readily estimated 
from it. How to obtain the total flexoelectric tensor within the current capabilities of ABINIT will be the 
subject of a future publication. 

The interested user must be likewise aware of the physical ambiguities existing in the definition of the bulk 
flexoelectric tensor which inherently affect the longwave driver. One of them precludes its usage to obtain the 
flexoelectric tensor of non-centrosymmetric (i.e., piezoelectric) materials (see section VII.c of [[cite:Stengel2013]]). 
The other one is related with the dependence of the bulk flexoelectric coefficients on the choice of an arbitrary 
reference energy. The average electrostatic potential has been taken as the energy reference 
within the longwave driver. Nonetheless, as illustrated in section IV.d of [[cite:Stengel2016a]], other choices might 
lead to quantitatively and qualitatively different outcomes. On the other hand, such ambiguity is well known to disappear
when the surface-specific part is accounted for, as done e.g. in [[cite:Stengel2014]].

The longwave implementation is still under heavy development. To date it requires the use of norm-conserving 
pseudopotentials without XC nonlinear core corrections and is limited to LDA XC functionals. 
The use of spherical harmonics for the nonlocal projectors is mandatory through the option [[useylm]]=1.   

The following steps are required to perform a longwave DFPT calculation of the bulk flexoelectric tensor
(see, e.g., tests [[test:lw_1]] to [[test:lw_3]]):

* Perform ground-state calculation.
* Perform ddk and d2_dkdk response function calculations.
* Perform response function calculations at **q** =Γ to atomic displacements, electric field and strain perturbations, 
including the option [[prepalw]]=1.
* Perform a longwave DFTP calculation of third-order energy derivatives ([[optdriver]]=10 and [[lw_flexo]]=1).
* Use MRGDDB to merge 1st, 2nd and 3rd order DDB files.
* Run ANADDB with [[anaddb:flexoflag]]=1.  

The following steps are required to perform a phonon dispersion calculation including quadrupolar fields in the
nonanalytical part of the dynamical matrix (see, e.g., tests [[test:lw_4]] to [[test:lw_6]]):

* Perform ground-state calculation.
* Perform ddk and d2_dkdk response function calculations.
* Perform response function calculations at **q** =Γ to atomic displacements and electric field, 
including the option [[prepalw]]=1.
* Perform a longwave DFTP calculation of third-order energy derivatives ([[optdriver]]=10 and [[lw_qdrpl]]=1).
* Perform response function calculations to atomic displacements at finite **q** (coarse grid). 
* Use MRGDDB to merge 2nd and 3rd order DDB files.
* Run a phonon dispersion calculation of ANADDB including [[anaddb:dipquad]]=1 and/or [[anaddb:quadquad]]=1.  


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

A tutorial is in preparation.

