---
description: How to compute phonon frequencies and modes, IR and Raman spectra, Born effective charges, IR reflectivity ...
authors: MT
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to compute phonon frequencies and modes, IR and Raman spectra, Born effective
charges, IR reflectivity with the ABINIT package.

## Introduction

The computation of the second-order derivative of the total energy with
respect to atomic displacements at an arbitrary wavevector, using
[[topic:DFPT]], opens the possibility to compute the dynamical matrix at that
wavevector, and hence, to compute the phonon eigenfrequency and
eigendisplacements. When the wavevector is (0,0,0), usually denoted as the
Gamma point, the combination of the atomic displacements and electric field
type perturbations opens also the access to Born effective charges, electronic
(for frequencies lower than the electronic band gap) dielectric constants, and
then, to infra-red reflectivity of materials (in the infinite lifetime
approximation). See [[cite:Gonze1997a]] for the presentation of the theory
with DFPT.

In ABINIT, with one dataset for a fixed wavevector (see [[topic:q-points]]),
one can compute all such second-order derivatives. ABINIT will already perform
some post-processing treatment of the second-order derivatives (e.g.
computation of the dynamical matrix, and corresponding eigenenergies and
eigendisplacements), although the most extended post-processing treatment is
provided by ANADDB. Thus, there is some overlap of the two executables, with
some common input variables. Usually, the action of an input variable with the
same name in the two executables is very similar, although there are some
input variables that govern more options in ANADDB then in ABINIT, because of
the previously mentioned difference in capabilities. In the database of input
variables, the input variables related to ABINIT or ANADDB are clearly
distinguished.

The band-by-band decomposition of the Born effective charge tensors can be
computed thanks to [[prtbbb]]. The related localization tensor (see
[[cite:Veithen2002]] can also be computed.

Phonon calculations are arbitrary q-points can be done under finite electric
field ([[topic:Berry]]).

It will be the easiest to discover the capabilities of these two executables
through the [[tutorial:rf1]] of the tutorial.

See [[topic:DFPT]] for the general information about DFPT, [[topic:q-points]]
for the specification of q-points, and [[topic:PhononBands]] for the
computation of full phonon bands.

!!! important

    More than 1500 phonon band structures for insulators, computed with ABINIT, are now available 
    on the [Materials Project web site](https://materialsproject.org), accompanied with derived 
    thermodynamic quantities, Born effective charges, and dielectric tensor [[cite:Petretto2018a]].
    The DDB file can be downloaded automatically with |AbiPy| starting from the materials project
    identifier. For futher information, please consult the |DdbFileNb|.

## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* [[tutorial:rf1|The tutorial Response-Function 1 (RF1)]] presents the basics of DFPT calculations within ABINIT. The example given is the study of dynamical and dielectric properties of AlAs (an insulator): phonons at Gamma, dielectric constant, Born effective charges, LO-TO splitting, phonons in the whole Brillouin zone. The creation of the "Derivative Data Base" (DDB) is presented.

