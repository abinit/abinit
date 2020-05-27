---
authors: MG, GB
---

# An overview of the EPH code

This page provides a quick introduction to the new EPH driver integrated with the ABINIT executable.
We discuss important technical details related to the implementation and the associated input variables.
The drawbacks/advantages with respect to the implementation available in ANADDB are also discussed.

## Why a new EPH code?

First of all, let's try to answer the question:
*why did we decide to implement a new code for electron-phonon calculations?*

As you may know, e-ph calculations have been available through the ANADDB executable for a long time.
This ANADDB-based implementation is essentially a post-processing of the e-ph matrix elements
computed by the DFPT code. This approach presents advantages as well as drawbacks.
On the one hand, most of the work required to compute e-ph matrix elements is implemented directly by the DFPT routines.
This means that e-ph calculations with advanced features such as PAW, SOC, non-collinear magnetism *etc*
are readily available in the ANADDB version once support in the DFPT part is implemented.
On the other hand, this post-processing approach implies that the number of $\kk/\qq$-points in the e-ph matrix elements
is automatically fixed at the level of the DFPT run.
In other words, if you want to compute phonon-limited mobilities with e.g. a 90×90×90 $\kk$- and $\qq$-mesh,
you need to perform DFPT calculations with the same sampling thus rendering the computation quite heavy.
In principle, it is possible to use tricks such as a linear interpolation to densify the sampling
inside ANADDB, but in order to get a decent interpolation one usually needs initial BZ meshes
that are significantly denser than the ones needed to converge the DFPT part alone.

As a matter of fact, electrons, phonons and e-ph properties present completely different convergence rates.
In silicon, for instance, a 9×9×9 mesh both for phonons and electrons is enough to converge
the electron density and the vibrational spectrum [[cite:Petretto2018]].
On the contrary, the phonon-limited mobility requires a 45×45×45 $\kk$-grid and a 90×90×90 $\qq$-grid
to reach a 5% relative error [[cite:Brunin2020]].
Roughly speaking, an explicit computation of phonons with a 90×90×90 $\qq$-mesh in Si
requires around 20000 × 3 × [[natom]] DFPT calculations
so you can easily get an idea of the cost of a fully ab-initio evaluation of the e-ph matrix elements
for all these wavevectors.

The EPH code bypasses this bottleneck by **interpolating the DFPT potentials** in $\qq$-space
while Bloch states are computed non-self-consistently on arbitrarily dense $\kk$-meshes.
As a net result, the three problems (electrons, phonons and electron-phonon) are now partly
decoupled and can be converged separately.
Keep in mind, however, that the fact that one can easily densify the sampling in the EPH code does not 
mean one can use under-converged values for the GS/DFPT part. 
Indeed, the quality of the ingredients used for the interpolation depends on the initial $\kk$- and $\qq$-meshes.
For further information about the difference between EPH and ANADDB, see also [[cite:Gonze2019]].
Further details about the EPH implementation are available in [[cite:Brunin2020]].

## Structure of an EPH calculation

A typical EPH workflow with arrows denoting dependencies between the different steps
is schematically represented in the figure below:

![](eph_intro_assets/eph_workflow.png){: style="height:400px;width:400px"}

The brown boxes are standard DFPT calculations done with relatively coarse $\kk$- and $\qq$-meshes (for instance, 9×9×9 in Si).
Each DFPT run produces a (partial) DDB file with a portion of the full dynamical matrix
as well as POT files with the first-order derivative of the KS potential (referred to as the DFPT potential below).
The partial POT files are merged with the **mrgdv** utility to produce a
single **DVDB file** (Derivatives of V($\rr$) DataBase).
As usual, the partial DDB files are merged with **mrgddb**.

The EPH driver (blue box) receives as input the total DDB and the DVDB as well as a GS WFK file that is usually
produced with a different $\kk$-mesh (in some cases, even with a different number of bands).
These ingredients are then used to compute the EPH matrix elements and associated physical properties.
<!--
The $\kk$-mesh in the WFK file must be commensurate with the $\qq$-mesh in the DVDB file
-->

The EPH calculation is activated by [[optdriver]] = 7 while
[[eph_task]] defines the physical properties to be computed.
To read the external files, one specifies the path on the file system with the three variables:
[[getwfk_filepath]], [[getddb_filepath]] and [[getdvdb_filepath]].

Internally, the code starts by reading the DDB file to construct the interatomic
force constants (IFCs) in the real-space supercell ($\RR$-space).
<!-- Other external files (WFK, DVDB) may be read depending on the value of [[eph_task]]. -->
Then, the EPH code computes phonon bands and DOS.
Finally a specialized routine is invoked depending on the value of [[eph_task]].

The following physical properties can be computed:

* Imaginary part of ph-e self-energy in metals ([[eph_task]] 1)

    * phonon linewidths induced by e-ph
    * Eliashberg function
    * superconducting properties within the isotropic Migdal-Eliashberg formalism
<!--
    * transport properties in metals with the LOVA approximation.
-->

* Real and imaginary parts of the e-ph self-energy ([[eph_task]] 4)

    * zero-point renormalization of the band gap
    * correction of the QP energies due to e-ph scattering
    * spectral function
    * Eliashberg functions

* Imaginary part of e-ph self-energy only ([[eph_task]] -4)

    * e-ph scattering rates
    * phonon-limited carrier mobility, electrical conductivity and Seebeck coefficient
    * phonon-limited carrier mean free path and relaxation times
    * all the calculations above can be done as a function of temperature and doping, for nonpolar and polar materials.

Other more advanced options are available...

<!--
Crystalline symmetries are used throughout the code in order to reduce the number
of $\kk$- and $\qq$-points that must be explicitly treated.
To achieve good parallel efficiently, the most CPU demanding parts are parallelized with MPI employing
a distribution schemes over $\qq$-points, perturbations
and bands (the band level is available only when computing the full self-energy).
[[istwfk]] [[kptopt]]
Features available in ANADDB that are not yet supported by EPH
-->

At the time of writing (|today|), the following features are not supported by EPH:

* PAW calculations
* Spin-orbit coupling
* Non-collinear magnetism ([[nspinor]] = 2 and [[nspden]] = 4)
* Non-local part of the pseudopotential applied with [[useylm]] = 1

In this introduction, we mainly focus on the parts common to the different sub-drivers:

1. Computation of vibrational properties via Fourier interpolation of the dynamical matrix.
2. Fourier interpolation of the DFPT potentials.

The use of the different sub-drivers is discussed in more detail in the specialized lessons.

## Phonon bands and DOS with EPH

Since phonon frequencies and displacements are needed for e-ph calculations, it is not surprising to see
that some of the ANADDB features related to the treatment of the dynamical matrix
are integrated in the EPH code as well.
In many cases, EPH uses the same name as in ANADDB especially for important variables
such as [[dipdip]], [[asr]], and [[chneut]].
There are however some differences with respect to the ANADDB interface.
More specifically, in EPH the name of the DDB file is specified by
[[getddb_filepath]] whereas the $\qq$-mesh associated to the DDB file is given by [[ddb_ngqpt]].
<!-- list of IBZ $\qq$-points for which the DFPT calculations have been performed -->
These two variables are **mandatory** when performing EPH calculations.

!!! importat

    Note that in Abinit9 the default values of [[dipdip]], [[asr]], and [[chneut]] have been changed.
    In particular, the acoustic sum rule and charge neutrality are always enforced by the EPH code.
    It is responsability of the user to check whether the breaking of these sum rules (always
    present due to numerical inaccuracies) are reasonable.
    By the same token, make sure that no vibrational instabilty is present in the phonon spectrum before
    embarking on EPH calculations.

### Variables for phonon DOS

By default, the EPH code computes the phonon DOS by interpolating the IFCs on the *dense* $\qq$-mesh
specified by [[ph_ngqpt]].
The default $\qq$ grid is 20×20×20, you may want to increase this value for more accurate results. 
The step of the (linear) frequency mesh is governed by [[ph_wstep]].
The linear tetrahedron integration method [[cite:Bloechl1994]] is used by default.
The Gaussian method can be activated via [[ph_intmeth]] with [[ph_smear]] defining the Gaussian smearing (in Hartree units by default).
The final results are stored in the PHDOS.nc file (same format at the one produced by ANADDB).
The computation of the PHDOS can be disabled by setting [[prtphdos]] = 0.

### Variables for phonon band structure

The computation of the phonon band structure is activated automatically, provided the input file
defines the high-symmetry $\qq$-path in terms of [[ph_nqpath]] vertices listed in the [[ph_qpath]] array.
[[ph_ndivsm]] defines the number of divisions used to sample the smallest segment.
The computation of the phonon band structure can be deactivated by setting [[prtphbands]] = 0.
The final results are stored in the PHBST.nc file (same format at the one produced by ANADDB).

## E-ph matrix elements

The e-ph matrix elements $\gkq$ are defined by

\begin{equation}
\gkq = \langle \psi_{m\kk+\qq}|\Delta_\qnu V^\KS|\psi_{n\kk} \rangle,
\label{eq:elphon_mel}
\end{equation}

where $\psi_{n\kk}$ is the KS Bloch state and $\Delta_\qnu V^\KS$ is the first-order variation of the
self-consistent KS potential induced by the phonon mode $\qnu$.
The scattering potential can be expresses as:

\begin{equation}
    \Delta_{\qq\nu} V^\KS(\rr) = e^{i\qq\cdot\rr} \Delta_{\qq\nu} v^\KS(\rr).
\label{eq:delta_vks}
\end{equation}

where $\Delta_{\qq\nu} v^\KS(\rr)$ is a lattice periodic function [[cite:Giustino2017]].
Note that ABINIT computes the response to an atomic perturbation defined by
the three variables [[qpt]], [[rfdir]] and [[rfatpol]] when [[rfphon]] is set to 1.
The connection between the phonon representation and the atomic perturbation employed by ABNIIT is given by
the equation

\begin{equation}
    \Delta_{\qq\nu} v^\KS(\rr) =
    \dfrac{1}{\sqrt{2\wqnu}} \sum_{\kappa\alpha}
    \dfrac{{e}_{\kappa\alpha,\nu}(\qq)}{\sqrt{M_\kappa}}\,\partial_{\kappa\alpha,\qq} v^\KS(\rr),
\end{equation}

where ${e}_{\kappa\alpha,\nu}(\qq)$ is the $\alpha$-th Cartesian component of the phonon eigenvector
for the atom $\kappa$ in the unit cell, $M_\kappa$ its atomic mass and
$\partial_{\kappa\alpha,\qq} v^\KS(\rr)$ is the first-order derivative of the KS potential
that can be obtained with DFPT by solving **self-consistently** a system of Sternheimer equations for a given
$(\kappa\alpha, \qq)$ perturbation [[cite:Gonze1997]] [[cite:Baroni2001]].

The DVDB file stores $\partial_{\kappa\alpha,\qq} v^\KS(\rr)$
for all the $\qq$-points in the IBZ and all the irreducible atomic perturbations.
In a more rigorous way, we should say that the DVDB file stores the local part of the DFPT potential
(variation of the Hartree + XC + local part of the pseudo)
but this is a rather technical point discussed in more detail in [[cite:Brunin2020]] that is not relevant
for the present discussion so we do not elaborate more on this.

<!--
In a pseudopotential-based implementation the KS potential is given by

\begin{equation}
    \begin{split}
    V^\KS(\rr,\rr') = \underbrace{\Bigl[ \VH[n](\rr)  + \Vxc[n](\rr) + \Vloc(\rr) \Bigr] }_{\Vscf(\rr)} & \\
    \times \delta(\rr-\rr') + \Vnl(\rr,\rr')
    \end{split}
\end{equation}

and consists of contributions from the Hartree part ($\VH$), the exchange-correlation (XC) potential ($\Vxc$),
and the bare pseudopotential term that, in turn, consists of the local ($\Vloc$) and non-local ($\Vnl$) parts~\cite{Payne1992}.
Following the internal \abinit convention, we group the Hartree, XC
and local terms in a single potential, $\Vscf$, although only the first two terms are computed self-consistently.
The lattice-periodic part of the first-order derivative of the KS potential thus reads

\begin{equation}
    \begin{split}
    \partial_{\kappa\alpha,\qq} v^\KS = \underbrace{ \Bigl [
        \partial_{\kappa\alpha,\qq} \vH +
        \partial_{\kappa\alpha,\qq} \vxc +
        \partial_{\kappa\alpha,\qq} \vloc \Bigr ]}_{\partial_{\kappa\alpha,\qq} \vscf} & \\
    +\, \partial_{\kappa\alpha,\qq} \vnl.
    \end{split}
    \label{eq:dvKS}
\end{equation}
-->

!!! important

    The number of irreducible atomic perturbations depends on the $\qq$-point and the number of symmetries [[nsym]].
    Actually only a particular subset of the point group symmetries of the unperturbed crystal can be
    exploited at the DFPT level.
    Fortunately, you don't have to worry about these technical details as symmetries are fully supported in EPH.
    Just run the DFPT calculation for the irreducible perturbations in the IBZ as usual.
    Also the computation of the WFK can be limited to the IBZ.
    EPH will employ symmetries to reconstruct the missing terms at runtime.

## Fourier interpolation of the DFPT potentials

The EPH code employs the Fourier interpolation proposed in [[cite:Eiguren2008]]
to obtain the scattering potentials at arbitrary $\qq$-points.
In this interpolation technique, one uses the DFPT potential stored in the DVDB file
to compute the $\qq \rightarrow \RR$ Fourier transform

\begin{equation}
	\label{eq:dfpt_pot_realspace}
    W_{\kappa\alpha}(\rr,\RR) = \dfrac{1}{N_\qq} \sum_\qq e^{-i\qq\cdot(\RR - \rr)}\,
    \partial_{\kappa\alpha\qq}{v^{\text{scf}}}(\rr),
\end{equation}

where the sum is over the BZ $\qq$-points belonging to the ab-initio [[ddb_ngqpt]] grid.
<!--
and $\partial_{\kappa\alpha\qq}{v^{\text{scf}}}$ represents the (lattice-periodic) first order derivative
of the local part of the KS potential associated to atom $\kappa$ along the Cartesian direction $\alpha$.
-->
Once $W_{\kappa\alpha}(\rr,\RR)$ is known, it is possible to interpolate the potential
at an arbitrary point $\tilde{\qq}$ using the inverse transform

\begin{equation}
	\label{eq:dfpt_pot_interpolation}
    \partial v^{scf}_{\tilde\qq\kappa\alpha}(\rr) \approx \sum_\RR e^{i\tilde{\qq}\cdot(\RR - \rr)} W_{\kappa\alpha}(\rr,\RR).
\end{equation}

where the sum is over the lattice vectors inside the Born-von Karman supercell.
The algorithm used to define the $\RR$ points of the supercell with the
corresponding weights is specified by [[dvdb_rspace_cell]].

The accuracy of the interpolation depends on the localization in $\RR$-space of $W_{\kappa\alpha}$.
This means that the Born-von Karman supercell corresponding to the [[ddb_ngqpt]] grid should be large
enough to capture the spatial decay of $W_{\kappa\alpha}(\rr,\RR)$ as a function of $\RR$.
As a consequence, [[ddb_ngqpt]] should be subject to convergence studies.



In metals, $W_{\kappa\alpha}$ is expected to be short-ranged provided we ignore possible Kohn anomalies.
On the contrary, a special numerical treatment is needed in semiconductors and insulators due to the presence of
long-ranged (LR) **dipolar** and **quadrupolar** fields in $\RR$-space.
These LR terms determine a non-analytic behaviour of the scattering potentials
in the long-wavelength limit $\qq \rightarrow 0$ [[cite:Vogl1976]].
To handle the LR part, EPH uses an approach that is similar in spirit to the one employed
for the Fourier interpolation of the dynamical matrix [[cite:Gonze1997]].
The idea is relatively simple.
One subtracts the LR part from the DFPT potentials before computing Eq. \eqref{eq:dfpt_pot_realspace}.
thus making the real-space representation amenable to Fourier interpolation.
The non-analytical part
<!-- (Eq.\eqref{eq:v1_long_range}) -->
is then restored back to Eq. \eqref{eq:dfpt_pot_interpolation} when interpolating the potential at $\tilde{\qq}$.

<!--
the long-range part associated to the displacement of atom $\kappa$ along the cartesian direction $\alpha$ can be modeled with
-->
In polar materials, the leading term is given by the dipolar field [[cite:Verdi2015]], [[cite:Sjakste2015]], [[cite:Giustino2017]],

\begin{equation}
   \label{eq:v1_long_range}
    V^{\mathcal{L}}_{\kappa\alpha\qq}(\rr) = i \dfrac{4\pi}{\Omega} \sum_{\GG \neq -\qq}
    \dfrac{(\qG)_\beta\cdot {\bm Z}^*_{\kappa\alpha,\beta}\,
    e^{i (\qG) \cdot (\rr - {\bm{\tau}}_{\kappa})}} {(\qG) \cdot {\bm{\varepsilon}}^\infty \cdot (\qG)},
\end{equation}

where ${\bm{\tau}}_\kappa$ is the atom position, $\Omega$ is the volume of the unit cell,
$\bm{Z}^*$ and ${\bm{\varepsilon}}^\infty$ are the Born effective charge tensor
and the dielectric tensor, respectively, and summation over the cartesian directions $\beta$ is implied.
This term diverges as $1/q$ for $\qq \rightarrow 0$ but the singularity is integrable in 3D systems.
$\bm{Z}^*$ and ${\bm{\varepsilon}}^\infty$ are read automatically from the DDB file (if present) so we
**strongly recommend** to compute these quantities with DFPT in order to prepare an EPH calculation
in semiconductors.

In non-polar materials, the Born effective charges are zero but the scattering potentials are still non-analytic
due to presence of jump discontinuities.
As discussed in [[cite:Brunin2020]] the non-analytic behaviour is fully captured by using:
<!-- 
In this case the leading term is associated to the dynamical quadrupoles [[cite:Royo2019]].
The expression for the LR model including both dipole and quadrupole terms reads:
-->

\begin{equation}
    V^{\mathcal{L}}_{\kappa\alpha,\qq}(\rr) =
    %\partial_{\kappa\alpha,\qq} v^{\Lc}(\rr) =
    \frac{4\pi}{\Omega} \sum_{\GG\neq\mathbf{-q}}
    \frac{ i(\qG)_\beta
    Z^*_{\kappa\alpha,\beta}
    -
    (\qG)_\beta
    ({Z^*_{\kappa\alpha,\beta}\cal{Q}}v^{\text{Hxc},\cal{E}_{\gamma}}(\rr)- \frac{1}{2}Q_{\kappa\alpha}^{\beta\gamma})(\qG)_\gamma
    }
    {(\qG) \cdot {\bm{\varepsilon}}^\infty \cdot (\qG)}
    e^{i (\qG) \cdot (\rr - {\bm{\tau}}_{\kappa})}
\label{eq:LRpot}
\end{equation}

TODO: Discuss more the integration with the DFPT part.

In the practical implementation, each Fourier component is multiplied by the Gaussian $e^{-\frac{|\qG|^2}{4\alpha}}$
in order to damp the short range components that are not relevant for the definition of the LR model.
The $\alpha$ parameter can be specified via the [[dvdb_qdamp]] input variable.
The default value worked well in the majority of the systems investigated so far yet
this parameter should be subject to convergence tests.
For testing purposes, it is possible to deactivate the treatment of the LR part by setting [[dvdb_add_lr]] to 0.
It can also be set to -1 if one is interested only in the short-range part of the e-ph scattering potentials.

### On the importance of the initial DFPT mesh

At this point it is worth commenting about the importance of the initial DFPT $\qq$-mesh.
The Fourier interpolation implicitly assumes that the signal in $\RR$-space decays quickly hence
the quality of the *interpolated* phonon frequencies and of the *interpolated* DFPT potentials,
between the ab-initio points depends on the spacing of the initial $\qq$-mesh that
in turns defines the size of the Born-von-Karman supercell.
In other words, the denser the DFPT mesh, the bigger the real-space supercell and the better the interpolation.
<!--
In semiconductors the atomic displacement induces dynamical dipoles and quadrupoles at the level of the density
that will generate long-range scattering potentials.
These potentials affect the behaviour of the e-ph matrix elements for $\qq \rightarrow 0$ and the
$\qq$-mesh must be dense enough to capture the full strenght of the coupling.
A more detailed discussion can be found in [[cite:Brunin2020]], [[cite:Verdi2015]] and [[cite:Sjakste2015]].
-->

From a more practical point of view, this implies that one should always monitor the convergence of the
physical properties with respect to the initial DFPT $\qq$-mesh.
The LR model implemented in ABINIT facilitates the convergence as the non-analytic behaviour for
$\qq \rightarrow 0$ is properly described yet the Fourier interpolation can introduce oscillations
between the *ab-initio* $\qq$-points and these oscillations may affect the quality of the physical results.
<!--
Note that the same consideration holds for Wannier-based approaches
independently on the degree of localization of the (maximally-localized) Wannier functions.
-->

## Tricks to accelerate the computation and reduce the memory requirements

Each sub-driver implements tricks to accelerate the calculation and reduce the memory requirements.
Here we focus on the techniques that are common to the different EPH sub-drivers.
Additional tricks specific to [[eph_task]] are discussed in more detail in the associated lesson.

First of all, note that the memory requirements for the $W_{\kappa\alpha}(\rr,\RR)$ array 
scales as [[nfft]] × product([[ddb_ngqpt]]). 
This quantity should be multiplied by 3 * [[natom]] if each MPI process stores all the perturbations in memory.
The MPI parallelism over perturbations (see [[eph_np_pqbks]]) allows one to decrease this prefactor from 
*3 × natom* to *3 × natom / nprocs_pert*.

Also, the total number of $\rr$-points ([[nfft]]) plays an important role both at the level of memory 
as well as the level of wall-time.
To optimize this part, one can decrease the value [[boxcutmin]] in the EPH calculation
to a value smaller than 2 e.g. 1.5 or the more aggressive 1.1.
You are not obliged to run the GS/DFPT part with the same [[boxcutmin]].
An exact representation of densities/potentials in $\GG$-space is obtained with [[boxcutmin]] = 2,
but we found that using a value of 1.1 does not change the result
but allows one to decrease the cost of the calculation and the memory by a factor ~8.

A significant fraction of the wall-time in EPH is spent for performing the FFTs required to apply $H^1$.
The use of single precision in the FFT routines allows one to decrease the computational cost without losing precision. 
This trick is activated by setting [[mixprec]] = 1 (support for FFTW3 or DFPT-MKL is required to take advantage 
of mixed precision FFTs).

!!! important

    The [[boxcutmin]] and [[mixprec]] tricks **are not activated by default** 
    because users are supposed to perform preliminary tests
    to make sure the quality of the results is not affected by these options.

We terminate the discussion with another trick that is not directly related to the EPH code but 
to the DFPT computation.
Since EPH does not need the first order change of the wavefunctions (1WF files)
we suggest to avoid the output of these files by using [[prtwf]] = -1 in the DFPT part
as these files are quite large and the overall space on disk scales as **nqpt × 3 × natom**.
When [[prtwf]] is set to -1, the DFPT code writes the 1WF only if the DFPT SCF cycle is not converged 
so that one can still restart from these wavefunctions if really needed 
(restarting a DFPT run from the 1WF file is more effective than restarting from the first order density).

<!--
On the one hand, this approach was relatively easy to implement as most of the work,
in particular the computation of the EPH matrix elements, was already performed by the DFPT code.
On the other hand, the resulting implementation was too rigid
as several important dimensions such as the number of $\kk$-points, $\qq$-points and bands
in the EPH matrix elements had to be fixed at the level of the DFPT calculation.
Performing convergence studies with respect to the $\kk$-point sampling, for instance,
required performing new (and more expensive) DFPT calculations with denser $\kk$-meshes.
Similarly, convergence studies for the $\qq$-points required additional DFPT computations, possibly
on meshes that were multiples of the initial sampling so to reuse the $\qq$-points computed previously.
A different philosophy is used, in which EPH matrix elements are computed directly starting from the basic ingredients, namely,
the GS wavefunctions stored in the WFK file, and the first-order change of the Kohn-Sham (KS) potential produced by the DFPT code.
This approach allows for more flexibility because electron and phonon calculations are now partly decoupled:
the $\kk$-mesh can be densified by performing non-self-consistent calculations,
thus bypassing the DFPT part, and interpolation schemes for the linear-response in $\qq$-space can be readily implemented.
Unlike the previous algorithms implemented in ANADDB, the new driver is directly interfaced with  the ABINIT executable.
-->
