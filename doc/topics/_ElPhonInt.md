---
description: How to compute the matrix elements of the electron-phonon interaction
authors: MV
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to compute the matrix elements of the electron-phonon interaction with the ABINIT package.

## Introduction

The theory and details of the implementation are described in [[cite:Gonze2009]] and [[cite:Gonze2016]].

Basic calculations of electron-phonon interaction in ABINIT: one performs a
normal ground state, then DFPT phonon calculations (using [[rfphon]], with
added keywords [[prepgkk]] and [[prtgkk]], which saves the matrix elements to
files suffixed GKK. The main change in this respect is that [[prtgkk]] now
disables the use of symmetry in reducing q-points and perturbations. This
avoids ambiguities in wave function phases due to band degeneracies. The
resulting GKK files are merged using the mrggkk utility, and processed by anaddb.

With the implementation of phonons in PAW DFPT, the electron phonon coupling
is also available in PAW, though this has not yet been tested extensively. The
input variables for electron-phonon coupling in anaddb are described in
[[cite:Gonze2009]] and [[cite:Gonze2016]].

Some details about the calculation of electron-phonon quantities in ABINIT and
ANADDB can be found [[pdf:elphon_manual.pdf|here]].

Subsequently, the GKK file is used to compute many quantities, as explained in
[[topic:PhononWidth]], [[topic:TDepES]] and [[topic:ElPhonTransport]].

A brand new ABINIT driver, focusing on the treatment of electron-phonon
interaction is under heavy development. Most of the input variables for experts,
with [[optdriver]]==7 are related to this development. It is operational
as of v9.2, although the documentation is not yet fully upgraded.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:rf1|The tutorial Response-Function 1 (RF1)]] presents the basics of DFPT calculations within ABINIT. The example given is the study of dynamical and dielectric properties of AlAs (an insulator): phonons at Gamma, dielectric constant, Born effective charges, LO-TO splitting, phonons in the whole Brillouin zone. The creation of the "Derivative Data Base" (DDB) is presented.

* (Legacy implementation) [[tutorial:eph|The tutorial on the electron-phonon interaction]] presents the use of the utility MRGKK and ANADDB to examine the electron-phonon interaction and the subsequent calculation of superconductivity temperature (for bulk systems).
Also there is a tutorial for [[tutorial:tdepes|the temperature dependence of the electronic structure]]:.

* (New implementation) Three tutorials for the new procedure are available at [[tutorial:eph_intro|an overview of the EPH code]], 
[[tutorial:eph4zpr|Zero-point renormalization of the band gap and temperature-dependent band gaps]], and
[[tutorial:eph4mob|Phonon-limited mobility]]
:.
