---
authors: MG
---

# An overview of the EPH code

This page provides a quick introduction to the new EPH driver integrated with the ABINIT executable.
We discuss important technical details related to the implementation and the associated input variables.
The drawbacks/advantages with respect to the ANADDB implementation are also discussed.

## Why a new EPH code?

First of all, let's try to answer the question: 
*why did we decide to implement a new code for electron-phonon calculations?*

As you may know, e-ph calculations have been available through the ANADDB executable for a long time.
This ANADDB-based implementation is essentially a post-processing of the e-ph matrix elements
computed at the end of the DFPT run. This approach, however, presents advantages as well as drawbacks.
On the one hand, most of the work required to compute e-ph matrix elements is implemented directly by the DFPT code.
This means that advanced features such as PAW, SOC, non-collinear magnetism etc are readily available in the ANADDB
version once support in the DFPT code is implemented.
On the other hand, this post-processing approach implies that the number of $\kk/\qq$-points in the e-ph matrix elements
is automatically fixed at the level of the DFPT run.

In other words, if you want to compute phonon-limited mobilities with e.g. a 90x90x90 $\kk$- and $\qq$-mesh,
you need to perform DFPT calculations with the same sampling thus rendering the computation quite heavy.
In principle, it is possible to use tricks such as a linear interpolation to densify the sampling
inside ANADDB, yet in order to get a decent interpolation we usually need initial BZ meshes that are significantly 
denser than the ones needed to converge the DFPT part alone.

Electrons, phonons and e-ph properties, indeed, present completely different convergence rates.
In silicon, for instance, a 9x9x9 mesh both for phonons and electrons is enough to converge
the electron density and the vibrational spectrum.
On the contrary, phonon-limited mobilities require e.g. a 45×45×45 k-grid and a 90×90×90 q-grid
to reach a 5% relative error [[cite:Brunin2020]].
Roughly speaking, an explicit computation of phonons with a 90×90×90 $\qq$-mesh in Si
requires ~20000 x 3 * [[natom]] DFPT calculations
so you can easily get an idea of the cost of a fully ab-initio DFPT evaluation.

The EPH code bypasses this bottleneck by interpolating the DFPT potentials in $\qq$-space
while Bloch states are computed explicitly on arbitrarily dense $\kk$-meshes with a NSCF run.
As a net result, the phonon and electron problems are now partly decoupled and can be converged separately.
For further information about the difference between the two approaches, see also [[cite:Gonze2019]] and [[cite:Brunin2020]].

<!--
Features available in anaddb that are not yet supported by EPH.

At the time of writing (today), the following features are not supported by EPH:

* PAW calculations
* Non-local part applied with [[useyml]] = 1
* Spin-orbit coupling
* Non-collinear magnetism ([[nspinor]] 2 with [[nspden]] 4

Crystalline symmetries are used throughout the code in order to reduce the number of $\kk$- and $\qq$-points
that must be explicitly included in the integrals.
To achieve good parallel efficiently, the most CPU demanding parts are parallelized with MPI employing 
a distribution schemes over $\qq$-points, perturbations
and bands (the band level is available only when computing the full self-energy).

[[istwfk]]


Note that all these capabilities are integrated directly in ABINIT.
This implementation (henceforth refered to as the **EPH code**) significantly differs from the one available in ANADDB:
the anaddb version acts as a direct post-processing of the e-ph matrix elements computed in the DFPT part
whereas the EPH code interfaced with ABINIT computes the e-ph matrix elements on the fly using
the GS WFK and the DFPT potentials stored in the DVDB file.
In a nutshell, the EPH code is more scalable and flexible as the $\qq$-sampling can be easily changed
at runtime while the anaddb implementation can easily support advanced features such as PAW as most of the
work is already done at the end of the DFPT calculation.
Electron-phonon (EPH) calculations have been available in ABINIT for a long time with the help
provided by the ANADDB tool designed as a post-processing step of the EPH matrix elements
computed at the end of the DFPT calculation.
On the one hand, this approach was relatively easy to implement as most of the work,
in particular the computation of the EPH matrix elements, was already performed by the DFPT code.
On the other hand, the resulting implementation was too rigid
as several important dimensions such as the number of $\kk$-points, $\qq$-points and bands
in the EPH matrix elements had to be fixed at the level of the DFPT calculation.
Performing convergence studies with respect to the $\kk$-point sampling, for instance,
required performing new (and more expensive) DFPT calculations with denser $\kk$-meshes.
Similarly, convergence studies for the $\qq$-points required additional DFPT computations, possibly
on meshes that were multiples of the initial sampling so to reuse the $\qq$-points computed previously.
To address these limitations, ABINIT v8 provides a new driver
explicitly designed to compute the EPH matrix elements and related physical properties.
A different philosophy is used, in which EPH matrix elements are computed directly starting from the basic ingredients, namely,
the GS wavefunctions stored in the WFK file, and the first-order change of the Kohn-Sham (KS) potential produced by the DFPT code.
This approach allows for more flexibility because electron and phonon calculations are now partly decoupled:
the $\kk$-mesh can be densified by performing non-self-consistent calculations,
thus bypassing the DFPT part, and interpolation schemes for the linear-response in $\qq$-space can be readily implemented.
Unlike the previous algorithms implemented in ANADDB, the new driver is directly interfaced with  the ABINIT executable.
This means that important ANADDB variables related to the computation and diagonalization
of the dynamical matrix such as [[asr]] and [[dipdip]] have been added to the ABINIT input file as well.
-->

A typical EPH workflow with arrows denoting dependencies between the different steps
is schematically represented in the below figure:

![](eph_intro_assets/eph_workflow.png){: style="height:400px;width:400px"}

<!--
Each box type represents a different kind of calculation.
self-consistent and non-self-consistent calculations to obtain wavefunctions, atomic perturbations
with respect to specific atoms and directions.
The DFPT potentials and the blocks of the dynamical matrix are merged and stored 
in two distinct files at the end of the DFPT part.
-->

The brown boxes represent standard DFPT calculations done with relatively coarse $\kk$- and $\qq$-meshes.
These calculations produce DDB files with dynamical matrix elements,
and POT files with the local part % (H + XC + vloc)
of the first-order change of the KS potential (referred to as the DFPT potential below).

A new utility, *mrgdv*, has been added to merge the DFPT potentials
in a single ``Derivative of V($\rr$) DataBase'' **DVDB** file, while
the partial DDB files can be merged, as in previous versions, with **mrgddb**.
The EPH driver (blue box) receives as input the total DDB and the DBDB as well as a GS WFK file that may
have been produced with a different $\kk$-mesh (or even with a different number of bands)
These ingredients are then used to compute the EPH matrix elements and associated physical properties.
%The $\kk$-mesh in the WFK file and the $\qq$-mesh in the DVDB file must be commensurate

The EPH calculation is activated by [[optdriver]] = 7 while
[[eph_task]] defines the physical properties to be computed.
Internally, the code starts by reading the DDB file to construct the interatomic force constants (IFCs) in $\RR$-space.
<!--
This part is very similar to what is done in ANADDB.
-->
Other external files (WFK, DVDB) may be read depending on the value of [[eph_task]].
At this point, the code computes phonon bands and phonon DOS.
Finally a specialized routine is invoked depending on [[eph_task]].

In this section, we mainy focus on the parts common to the different subdrivers:

* computation of vibrational properties via Fourier interpolation of the dynamical matrix
* Fourier interpolation of the DFPT potentials.

The usage of the different subdrivers is discussed in more detail in the specialized lessons.

## Phonon bands and DOS with EPH

<!--
ABINIT users know that in order to interpolate a phonon band structure,
one should first merge the partial DDB files with *mrgddb* and then use ANADDB.
-->
Since phonon frequencies and displacements are needed for e-ph calculations, it's not surprising to see
that some of the ANADDB features are now integrated in the EPH code as well.
In many cases, EPH uses the same name as in ANADDB especially for important variables
such as [[dipdip]], [[asr]], and [[chneut]]
There are however some differences with respect to the ANADDB interface.
More specifically, the name of the DDB file is specified by
[[getddb_filepath]] whereas the $\qq$-mesh associated to the DDB file
is given by [[ddb_ngqpt]].
<!-- list of IBZ $\qq$-points for which the DFPT calculations have been performed -->
These two variables are mandatory when performing EPH calculations.

!!! importat

    Note that in Abinit9 the default values of [[dipdip]], [[asr]], and [[chneut]] have been changed.

### Variables for phonon DOS

By default, the EPH code computes the phonon DOS by interpolating the IFCs on the *dense* $\qq$-mesh
specified by [[ph_ngqpt]].
The step of the (linear) frequency mesh is governed by [[ph_wstep]],
and the linear tetrahedron method is used by default.
The Gaussian method can be activated via [[prtphdos]] with [[ph_smear]] defining the Gaussian smearing.
The final results are stored in the PHDOS.nc file (same format at the one produced by ANADDB).
The computation of the PHDOS can be disabled by setting [[prtphdos]] = 0

### Variables for phonon band structure

In a typical EPH run, the computation of the phonon band structure is activated by default.
The $\qq$-path is specified in terms of [[ph_nqpath]] vertices listed in the [[ph_qpath]] array
while [[ph_ndivsm]] defines the number of divisions used to sample the smallest segment.
The computation of the phonon band structure can be deactivated by setting [[prtphbands]] = 0.
The final results are stored in the PHBST.nc file (same format at the one produced by ANADDB).

## Fourier interpolation of the DFPT potential

<!--
E-PH properties are rather sensitive to the BZ sampling and some sort of interpolation scheme
is needed to avoid explicit computations on very dense k/q meshes.
The approach used in the EPH code consists in interpolating the DFPT potentials in $\qq$-space
The advantage of such approach is that the interpolation in $\qq$-space in relatively easy
to implement without any use inte

In this document, we mainly focus on the input variables governing the interpolation of the DPFT potentials $\qq$-space.
In the other EPH tutorials, we discuss how to use smart tricks to reduce the number
of $\kk$-points that must be treated explicitly.
As the calculation of the DFPT potentials represents a significant fraction of the overall computational time,
especially when compared with the non-self-consistent computation of the WFK file,
the new EPH driver allows the user to densify the $\qq$-mesh for phonons using
-->

In EPH we use the Fourier interpolation proposed in [[cite:Eiguren2008]]
to obtain the scattering potentials at arbitrary $\qq$-points
In this method, one computes the Fourier transform

\begin{equation}
	\label{eq:dfpt_pot_realspace}
    W_{\kappa\alpha}(\rr,\RR) = \dfrac{1}{N_\qq} \sum_\qq e^{-i\qq\cdot(\RR - \rr)}\,
    \partial_{\kappa\alpha\qq}{v^{\text{scf}}}(\rr),
\end{equation}

where the sum is over the $\qq$-points belonging to the ab-initio [[ddb_ngqpt]] grid used for the DFPT calculation
and $\partial_{\kappa\alpha\qq}{v^{\text{scf}}}$ represents the (lattice-periodic) first order derivative
of the local part of the KS potential associated to atom $\kappa$ along the cartesian direction $\alpha$.
Once $W_{\kappa\alpha}(\rr,\RR)$ is known, one can interpolate the potential at an arbitrary point $\tilde{\qq}$
using the inverse transform

\begin{equation}
	\label{eq:dfpt_pot_interpolation}
    \partial v^{scf}_{\tilde\qq\kappa\alpha}(\rr) \approx \sum_\RR e^{i\tilde{\qq}\cdot(\RR - \rr)} W_{\kappa\alpha}(\rr,\RR).
\end{equation}

where the sum is over the lattice vectors inside the Born-von Karman supercell.
[[dvdb_rspace_cell]]

Note that the accuracy of the interpolation depends on the localization in $\RR$-space of $W$.
This means that the Born-von Karman supercell corresponding to [[ddb_ngqpt]] grid should be large
enough to capture the spatial decay of $W_{\kappa\alpha}(\rr,\RR)$ as a function of $\RR$.
It is worth stressing that the same consideration holds for Wannier-based approaches indipendently on the degree of localization
of the (maximally localized) Wannier functions.

In metals, $W$ is expected to be short-ranged provided we ignore Friedel oscilations associated to Kohn anomalies.
On the contrary, a special numerical treatment is needed in semiconductors and insulators.
<!--
the long-range behaviour of the DFPT potential in polar semiconductors and insulators is therefore problematic
and  is needed to enforce localization.
-->
To handle the long-range part, we use an approach that is similar in spirit to the one employed
for the Fourier interpolation of the dynamical matrix [[cite:Gonze1997]]
As discussed in [[cite:Verdi2015]], and [[cite:Giustino2017]],
the long-range part associated to the displacement of atom $\kappa$ along the cartesian direction $\alpha$ can be modeled with

\begin{equation}
   \label{eq:v1_long_range}
    V^{\mathcal{L}}_{\kappa\alpha\qq}(\rr) = i \dfrac{4\pi}{\Omega} \sum_{\GG \neq -\qq}
    \dfrac{(\qG)_\beta\cdot {\bm Z}^*_{\kappa\beta\alpha}\,
    e^{i (\qG) \cdot (\rr - {\bm{\tau}}_{\kappa})}} {(\qG) \cdot {\bm{\varepsilon}}^\infty \cdot (\qG)},
\end{equation}

where ${\bm{\tau}}_\kappa$ is the atom position, $\Omega$ is the volume of the unit cell,
$\bm{Z}^*$ and ${\bm{\varepsilon}}^\infty$ are the Born effective charge tensor and the dielectric tensor, respectively, and summation over the cartesian directions $\beta$ is implied.

Inspired by the approach used for the Fourier interpolation of the dynamical matrix in polar materials~\cite{Gonze1997},
we subtract the long-range part from the DFPT potentials before computing Eq.~\eqref{eq:dfpt_pot_realspace}.
%thus making the real-space representation amenable to Fourier interpolation.
%The non-analytical part (Eq.~\eqref{eq:v1_long_range}) is then restored back to Eq.~\eqref{eq:dfpt_pot_interpolation}
%when interpolating the potential at $\tilde{\qq}$.

The e-ph matrix elements $\gkkp$ are given by

\begin{equation}
\gkkp = \langle \psi_{m\kk+\qq}|\Delta_\qnu V^\KS|\psi_{n\kk} \rangle,
\label{eq:elphon_mel}
\end{equation}

with

$\Delta_\qnu V^\KS$ the first-order variation of the self-consistent KS potential that can be computed with DFPT

[[ddb_ngqpt]]
[[dvdb_add_lr]]
[[dvdb_qdamp]]
[[getdvdb_filepath]]
