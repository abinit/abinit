---
description: How to compute transport properties that are determined by the electron-phonon interaction (electrical resistivity, superconductivity, thermal conductivity)
authors: MV
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to compute transport properties that are determined by the electron-phonon
interaction (electrical resistivity, superconductivity, thermal conductivity) with the ABINIT package.

## Introduction

Warning : this topic concerns metals only.

The calculation of bulk transport quantities (electrical and thermal
resistivities - the part that is determined by the electron-phonon
interaction) is possible using anaddb. Analogous quantities are obtained from
the conducti post-processor, but due to electron-electron scattering, instead
of electron-phonon.

A preliminary calculation of the derivatives of the wavefunctions with respect
to k-vector must be carried out. After generating a GKK file (see
[[topic:ElPhonInt]]), the Electron-Phonon Coupling (EPC) analysis is performed
in anaddb, setting [[anaddb:elphflag]] variable to 1. Most of the procedure is
automatic, but can be lengthy if a large number of k-points is being used.

While the legacy implementation of the transport properties in ABINIT is quite stable,
there is a new implementation under heavy development.
Most of the present information relates to the legacy implementation,
although some also relates to the
most recent procedure that relies on [[optdriver]]=7. The documentation
of the new procedure is given mostly by the related tutorials (introduction and mobility), see below.
Another tutorial for the new procedure for superconductivity calculations is still under development.


For the superconductivity calculations (legacy implementation), The electron-phonon interaction is
interpolated in reciprocal space, then integrated over the Fermi surface to
give the Eliashberg function. Several quadrature methods are available. The
default ([[anaddb:telphint]]=1) is to use Gaussian weighting, with a width
[[anaddb:elphsmear]]. Another option is the improved tetrahedron
[[cite:Bloechl1994a]] ([[anaddb:telphint]]=0). Finally
([[anaddb:telphint]]=2), one can integrate a given set of electron bands,
between [[anaddb:ep_b_max]] and [[anaddb:ep_b_min]]. The resulting integrated
quantities are the Eliashberg function (in a file suffixed `_A2F`), and the EPC
strength λ which is printed in the main output file.

The transport calculation is turned on by setting [[anaddb:ifltransport]] to 1
in anaddb. The transport quantities depend on the Fermi velocity for each
band, and the electronic band-dependence of the matrix elements must be
preserved before integration, by setting [[anaddb:ep_keepbands]] to 1. This
increases the memory used, by the square of the number of bands crossing EF.
The results are the transport Eliashberg function (in file `_A2F_TR`), the
electrical resistivity (in file `_RHO`), and the thermal conductivity (in file `_WTH`).

It is also possible to consider the temperature dependence of the Fermi
energy: cubic spline interpolation ([[anaddb:ep_nspline]]) enables to linearly
interpolate the transport arrays and reduce the memory usage. Besides setting
the Fermi level with [[anaddb:elph_fermie]] (in Hartree), it is also possible
to specify the extra electrons per unit cell, (i.e., the doping concentration
often expressed in cm-3) with [[anaddb:ep_extrael]].

Some details about the calculation of electron-phonon quantities in ABINIT and
ANADDB can be found [[pdf:elphon_manual|here]].


## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

## Tutorials

* (Legacy approach) [[tutorial:eph|The tutorial on the electron-phonon interaction]] presents the use of the utility MRGKK and ANADDB to examine the electron-phonon interaction and the subsequent calculation of superconductivity temperature (for bulk systems).

* (New procedure) Two tutorials are available at [[tutorial:eph_intro|an overview of the EPH code]], and
at [[tutorial:eph4mob|Phonon-limited mobility]]:.
