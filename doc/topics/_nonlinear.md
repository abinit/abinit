---
description: How to compute Raman intensity, and the related electro-optic coefficients
authors: RC, XG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to compute Raman intensity, and the related electro-optic 
coefficients with the ABINIT package.

## Introduction

In Raman experiments, the incident light, usually a polarized or unpolarized
laser, is scattered by the sample, and the energy as well as polarization of
the outgoing light is measured. A Raman spectrum, presenting the energy of the
outgoing photons, will consist of rather well-defined peaks, around an elastic peak.

At the lowest order of the theory, the dominant mechanism is the absorption or
emission of a phonon by a photon. The energy of the absorbed or emitted phonon
corresponds to the energy difference between the outgoing and incident
photons. Thus, even more straightforwardly than the IR spectrum, a Raman
spectrum is directly related to the energy of phonons at the Brillouin-zone
center: when the zero of the energy scale is set at the incident light energy,
the absolute value of the energy of the peaks corresponds to the energy of the
phonons.

The computation of phonon energies is presented in [[topic:Phonons]]. Raman
intensities due to one-phonon emission or absorption are not linked to second-
order derivatives of the total energy, but, within the adiabatic
approximation, to derivative of the dielectric phonon with respect to atomic
displacements. Moreover, when the frequency of the incident light (usually in
the 1.5 eV to 2.5 eV range) is small with respect to the band gap (e.g. for
gaps larger than 4 eV), the static approximation can be made, in which the
Raman intensity will be linked to the third-order derivative of the total
energy with respect (twice) to an homogeneous electric field and (once) with
respect to atomic displacements. Thus, DFPT can be used, see below. For the
case in which the incident light frequency is not negligible with respect to
the gap, the DFPT cannot be used, but, if the adiabatic approximation can be
used (valid when the phonon frequency is much smaller than the gap, and 
also, as a consequence, features
smaller than the largest phonon frequency are not resolved in the Raman
spectrum), one can compute the Raman intensities thanks to finite differences
of dielectric function, see [[cite:Gillet2013]]. For the two-phonon Raman
spectrum, see [[cite:Gillet2017]].

Both the derivatives of the linear electronic dielectric susceptibilities with
respect to atomic displacements and the non-linear electronic dielectric
susceptibilities required to evaluate the Raman intensities are thus non-
linear responses. 

In the ABINIT implementation, they are computed within the
density functional perturbation theory.
A first formalism (PEAD) is described in [[cite:Veithen2005]].
Thanks to the 2n+1 theorem, this formulation only requires the knowledge of
the ground-state and first-order changes in the wavefunctions,
but some quantity is evaluated thanks to a finite-difference in k-space. 
It is implemented only for NC pseudopotentials, within LDA.
There is another formalism, based on 2nd-order Sternheimer equation, 
recently made available (L. Baguet, to be published).
The latter is implemented both for NC pseudopotentials and for PAW, again only for LDA.
Both work for [[nsppol]]=1 or 2.

The PEAD non-linear response formalism has been successfully applied to a large
variety of systems. We have so far studied the Raman spectra of ferroelectric
oxides ( BaTiO3 and PbTiO3 [[cite:Hermet2009]]), different minerals under
pressure conditions characteristic to the interior of the Earth
[[cite:Caracas2007a]] or molecular solids under extreme conditions
[[cite:Caracas2008]]. The computation of the non-linear optical
susceptibilities has also been applied to several polar dielectrics
[[cite:Caracas2007]].

As a by-product of the calculation of the Raman tensor and non-linear optical
coefficients, it is also possible to determine directly within ABINIT the
electro-optic (EO) coefficients rijγ (Pockels effect) which describe the
change of optical dielectric tensor in a (quasi-)static electric field through
the following expression [[cite:Veithen2005]]: Δ(ε-1)ij=∑γ=1,3 rijγΕγ

The clamped (zero strain) EO coefficients include an electronic and an ionic
contribution directly accessible within ABINIT. The unclamped EO coefficients
include an additional piezoelectric contribution which must be computed
separately from the knowledge of the elasto-optic and piezoelectric strain
coefficients. This formalism was for instance applied to different
ferroelectric ABO3 compounds [[cite:Veithen2005a]].


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

[[tutorial:nlo|The tutorial on static non-linear properties]] presents the
computation of responses beyond the linear order, within Density-Functional
Perturbation Theory (beyond the simple Sum-Over-State approximation): Raman
scattering efficiencies (non-resonant case), non-linear electronic
susceptibility, electro-optic effect. Comparison with the finite field
technique (combining the computation of linear response functions with finite
difference calculations), is also provided.

