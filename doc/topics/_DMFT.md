---
description: How to perform a DMFT calculation
authors: BAmadon
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform a DMFT calculation with the ABINIT package.

## Introduction

DFT fails to describe the ground state and/or the excited states such as many
lanthanides, actinides or transition metals. Indeed, exchange correlation
functionals are not (yet) able to describe the strong repulsive Coulomb
interactions occurring among electrons in partly filled localized d or f
orbitals.

A way to improve the description of strongly correlated systems is to
explicitly include the strong repulsive Coulomb interactions in the
Hamiltonian. Solving it in the static mean field approximation, gives the
DFT+U method ([[cite:Anisimov1991]], [[cite:Liechtenstein1995]]), implemented
in ABINIT [[cite:Amadon2008a]]. The Dynamical Mean Field Theory
[[cite:Georges1996]] (DMFT), goes beyond, by solving exactly the local
correlations for an atom in an effective field (i.e., an Anderson model). The
effective field reproduces the effect of the surrounding correlated atoms and
is thus self-consistently related to the solution of the Anderson model
[[cite:Georges1996]].

The combination of DFT with DMFT [[cite:Georges2004]], [[cite:Kotliar2006]]
([[usedmft]]= 1) relies on:

* The definition of correlated orbitals. In ABINIT, we use Wannier functions built using projected local orbitals 
  [[cite:Amadon2008]]. Wannier functions are unitarily related to a selected set of Kohn Sham (KS) wavefunctions, 
  specified in ABINIT by band indices [[dmftbandi]] and [[dmftbandf]]. As empty bands are necessary 
  to build Wannier functions, it is required in DMFT calculations that the KS Hamiltonian is correctly diagonalized: 
  use high values for [[nnsclo]] and [[nline]]. In order to make a first rough estimation of the orbital character 
  of KS bands and choose the band index, the band structure with highlighted atomic orbital character 
  (so called *fatbands* ) can be plotted, using the [[pawfatbnd]] variable. 
  Band structures obtained from projected orbitals Wannier functions can also be plotted using 
  [[plowan_compute]] and related variables. 

* The choice of the screened Coulomb interaction U ([[upawu]]) and J ([[jpawu]]). Note that up to version 7.10.5 (but not in later versions) [[jpawu]]= 0 is required if the density matrix in the correlated subspace is not diagonal.

* The choice of the double counting correction [[cite:Amadon2012]]. The current default choice in ABINIT is [[dmft_dc]]= 1 which corresponds to the full localized limit.

* The method of resolution of the Anderson model. In ABINIT, it can be the Hubbard I method 
  [[cite:Amadon2012]] ([[dmft_solv]]= 2), the Continuous time Quantum Monte Carlo (CTQMC) method 
  [[cite:Gull2011]],[[cite:Bieder2014]] ([[dmft_solv]]= 5) or the static mean field method ([[dmft_solv]]= 1), 
  equivalent to usual DFT+U [[cite:Amadon2012]]).

The practical solution of the DFT+DMFT scheme is usually presented as a double
loop over, first, the local Green's function, and second the electronic local
density [[cite:Amadon2012]]. The number of iterations of the two loops are
determined by [[dmft_iter]] and [[nstep]]. However, in the general case, the
most efficient way to carry out fully consistent DFT+DMFT calculations is to
keep only the loop governed by [[nstep]], while [[dmft_iter]]=1
[[cite:Bieder2014]], [[dmft_rslf]]= 1 (to read the self-energy file at each
step of the DFT loop) and [[prtden]]= -1 (to be able to restart the
calculation of each step of the DFT loop from the density file). Lastly, one
linear and one logarithmic grid are used for Matsubara frequencies
[[cite:Kotliar2006]] determined by [[dmft_nwli]] and [[dmft_nwlo]] (Typical
values are 10000 and 100, but convergence should be studied). More information
can be obtained in the log file by setting [[pawprtvol]]=3.

The main output of the calculations are the imaginary time Green's function ,
from which spectral functions can be obtained using an external maximum
entropy code [[cite:Bergeron2016]], self-energies, from which quasiparticle
renormalization weight can be extracted, the density matrix of correlated
orbitals, and the internal energies [[cite:Amadon2006]]. The electronic
entropic contribution to the free energy can also be obtained using
[[dmft_entropy]] and [[dmft_nlambda]].

The efficient CTQMC code in ABINIT, which is the most time consuming part of
DMFT, uses the hybridization expansion [[cite:Werner2006]], [[cite:Gull2011]]
with a _density-density_ multiorbital interaction [[cite:Gull2011]]. Moreover,
the hybridization function [[cite:Gull2011]] is assumed to be diagonal in the
orbital (or flavor) index. This is exact for cubic symmetry without spin orbit
coupling but, in general, one should always check that the off-diagonal terms
are much smaller than the diagonal ones. A link to the exact rotationally
invariant interaction CTQMC code of the TRIQS library is also available using
[[dmft_solv]]=7.

As the CTQMC solver uses a Fourier transform, the time grid [[dmftqmc_l]] in
imaginary space must be chosen so that the Nyquist frequency, defined by
π[[dmftqmc_l]] [[tsmear]], is around 2 or 3 Ha. A convergence study should
be performed on this variable. Moreover, the number of imaginary frequencies
([[dmft_nwlo]]) has to be set to at least twice the value of [[dmftqmc_l]].
Typical numbers of steps for the thermalization ([[dmftqmc_therm]]) and for
the Monte carlo runs ([[dmftqmc_n]]) are 106 and 109 respectively. The random
number generator can be initialized with the variable [[dmftqmc_seed]].
Several other variables are available. [[dmftctqmc_order]] gives a histogram
of the perturbation orders during the simulation, [[dmftctqmc_gmove]]
customizes the global move tries (mainly useful for systems with high/low spin
configurations), and [[dmftctqmc_meas]] sets the frequency of measurement of
quantities.



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:dmft|The tutorial on DFT+DMFT]] shows how to perform a DFT+DMFT calculation on SrVO3 using projected Wannier functions. Prerequisite: DFT+U.

