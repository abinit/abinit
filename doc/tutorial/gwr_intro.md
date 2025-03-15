---
authors: MG
---

# An overview of the GWR code

This page offers a concise introduction to the new GWR driver of ABINIT.
We outline the technical details of its implementation, the relevant input variables,
and compare its advantages and limitations against the conventional GW approach,
which is formulated in Fourier space and the real-frequency domain.
In the following, we will refer to this traditional approach as the **conventional** or **legacy** GW code of ABINIT.

## Why a new GW code?

In the legacy GW approach, one uses the Lehmann representation of the Green's function in the frequency domain
to compute the irreducible polarizability, as explained in the [MBPT notes](/doc/theory/mbt).
Then the self-energy matrix elements in the KS representation are computed via an expensive convolution
in the frequency domain [[cite:Golze2019]] — by default using the plasmon-pole approximation [[cite:Giantomassi2011]],
which significantly accelerates calculations but introduces approximations and prevents direct access to the spectral function $A(\omega)$.
The conventional GW code exhibits quartic scaling with the number of atoms in the unit cell
and quadratic scaling with the number of $\kk$-points.

In contrast, the GWR code achieves cubic scaling with [[natom]] and linear scaling in the number of $\kk$-points
by computing the polarizability and the self-energy in the real-space supercell associated to the $\kk$-mesh.
FFTs are used to transform quantities from the supercell representation to Fourier space whenever needed.
For instance the equation for W (a convolution between $\ee^-1$ and the bare Coulomb interaction is best solved in Fourier space.
As concerns the frequency dependence, GWR evaluates the $G$, $W$ and $\Sigma$ along the imaginary axis using minimax meshes
that are adaptive grids that place points non-uniformly, concentrating them where they are most needed,
leading to higher accuracy with fewer points.
This approach was proposed for the first time in ... and then implemented
Citations relevant to the minimax mesh [[cite:Azizi2023]], [[cite:Azizi2024]]

followed by an analytic continuation (AC) to the real-frequency axis.
The AC approach enables access to the full frequency dependence of $\Sigma$ and $A$
at a substantially reduced computational cost, though the accuracy of the results now depends on the effectiveness of the AC step.

The frequency dependence of the self-energy $\Sigma(\omega)$ requires numerical integration over a range of energies.
Traditionally, these calculations use uniform energy grids or plasmon-pole approximations.
However, uniform grids become computationally expensive because they require many points
to accurately capture sharp features in in the Green's function $G(\omega)$, and the screened interaction $W(\omega)$.

The conventional GW code exhibits quartic scaling with the number of $\qq$-points, whereas GWR achieves linear scaling.

The two codes strongly differ also at the level of the MPI parallelization.
In the legacy GW code, MPI parallelization is available over the [[nband]] states; however,
key data structures, such as the screened interaction $W$, are not MPI-distributed.
As a result, the maximum number of usable MPI processes is limited by [[nband]],
and the workload becomes imbalanced when the number of MPI processes does not evenly divide [[nband]].
More critically, self-energy calculations in the legacy GW code are highly memory-intensive,
as they require storing both the wavefunctions (whose memory footprint scales with the number of MPI processes)
and $W$ (a non-scalable portion).
This memory requirement becomes particularly problematic for large [[npweps]] values or calculations
beyond the plasmon-pole approximation, where the full $W$ matrix must be stored for multiple frequencies.
Consequently, conventional GW calculations can be prohibitively demanding in terms of memory, especially for large systems
or systems with few symmetries.

In contrast, the GWR code distributes most data structures across MPI processes,
which enables handling larger systems when sufficient computing nodes are available.
However, this distribution comes at the cost of increased MPI communication in certain parts of the algorithm.
More specifically, GWR leverages Parallel BLAS (PBLAS) to efficiently distribute the memory required
for storing the Green’s functions and $W$, significantly improving scalability compared to the legacy implementation.

* [[ngkpt]] 4 4 4
* [[nshiftk]] 1
* [[shiftk]] 0 0 0

!!! important

    At the time of writing, the following features are **not yet supported** in GWR:

    * PAW method
    * Metallic systems as the our minimax meshes assume systems with an energy gap
    * Spinor wave-functions ([[nspinor]] = 2)
    * Temperature effects at the electronic level are not taken into account as we work with the T = 0 formalism.
    * Only $\Gamma$-centered $\kk$-meshes are supported in GWR, e.g.:

Select the task to be performed when [[optdriver]] == 6 i.e. GWR code.
while [[gwr_task]] defines the task to be performed.

## Formalism

The zero-temperature Green's function in imaginary-time is given by:

\begin{equation}
G(\rr, \rr', i\tau) =
\Theta(\tau) \Gove(\rr, \rr', i\tau) +
\Theta(-\tau) \Gund (\rr, \rr', i\tau)
\end{equation}

with $\Theta$ the Heaviside step-function, and

\begin{equation}
\Gove(\rr, \rr', i\tau) =
-\sum_n^{\text{unocc}} \psi_n(\rr)\psi_n^*(\rr') e^{-\ee_n\tau}
\qquad (\tau > 0)
\end{equation}

\begin{equation}
\Gund(\rr, \rr', i\tau) =
\sum_n^{\text{occ}} \psi_n(\rr)\psi_n^*(\rr') e^{-\ee_n\tau}
\qquad (\tau < 0).
\end{equation}

!!! important

    For simplicity, in all the equations we assume a spin-unpolarized semiconductor with scalar wavefunctions
    (i.e. [[nsppol]] = 1 with [[nspinor]] = 1) and the Fermi level $\mu$ is set to zero.
    At the level of the implementation, this means that the initial KS eigenvalues produced by the NSCF part of ABINIT
    are shifted with respect to the value of $\mu$ at zero temperature that is now located at mid-gap.

The GWR code constructs the Green's function from the KS wavefunctions and eigenvalues stored
in the WFK file specified via [[getwfk_filepath]] or [[getwfk]] in multi-dataset mode.
This WFK file is usually produced by performing a NSCF calculation including empty states and the GWR driver
provides a specialized option to perform a **direct diagonalization** of the KS Hamiltonian in parallel
with Scalapack (see section below).
When computing the KS Green's function $G^0$, the number of bands included
in the sum over states is controlled by [[nband]].
Clearly, it does not make any sense to ask for more bands than the ones available in the WFK file.

Note that GWR also needs the GS density produced by a previous GS SCF run.
The location of the density file can be specified via [[getden_filepath]] or [[getden]].

The imaginary axis is sampled using a minimax mesh with [[gwr_ntau]] points.
The other piece of information required for the selection of the minimax mesh
is the ratio between the **fundamental** gap and the maximum transition energy i.e.
the differerence between the highest KS eigenvalue for the empty states that, in turns, depends
on the value of [[nband]] and the energy of the lowest occupied state.

The irreducible polarizability is computed in real space in the supercell using

\begin{equation}
\chi(\rr,\RR', i\tau) = G(\rr, \RR', i\tau) G^*(\rr, \RR', -i\tau)
\end{equation}

and then immediately transformed to Fourier space.

\begin{equation}
\chi_\kk(\bg, \bg') =
\sum_{\rr\in\mcC} e^{-i(\kk+\bg)\rr} \chi(\rr, \GG' = \kk+\bg')
\end{equation}

The cutoff energy for the polarizability is given by [[ecuteps]]
while [[ecutsigx]] defines the number of $\gg$-vectors for the exchange part of the self-energy
that is computed using the standard summation over occupied states:

\begin{equation}\label{eq:Sigma_x}
\Sigma_x(\rr_1,\rr_2)= -\sum_\kk^\BZ
\sum_\nu^\text{occ} \Psi_{n\kk}(\rr_1){\Psi^\*_{n\kk}}(\rr_2)\,v(\rr_1,\rr_2)
\end{equation}

As concerns the treatment of the long-wavelenght limit $\qq \rightarrow 0$ in the polarizability,
we have the following input variables:

[[inclvkb]],
[[gw_qlwl]],
[[gwr_max_hwtene]]

## Real-space vs convolutions in the BZ

So far we have discussed the GWR equations used to compute the polarizability
and the self-energy in the real-space supercell.
This is the recommended approach if one needs to compute QP corrections for all the $\kk$-points in the IBZ.
This is the typical scenario for self-consistent calculations or if one needs to perform some kind of interpolation
of the QP results to obtain e.g. a band structure along a high-symmetry $\kk$-path.

It should be noted, however, that in several applications one is mainly interested in the QP corrections
at the band edges that are usually located at high-symmetry $\kk$-points.
In this case, it is more advantageous to use an alternative formulation that evaluates the matrix elements of
$\Sigma_\kk$ in terms of convolutions in the BZ according to

where $G_\qq(\rr,\rr')$ and $W_\qq(\rr,\rr')$ are now defined in the unit cell.

In this formalism, one can take advantage of the symmetries of the system to reduce the BZ summation to
the irreducible wedge defined by the little group of the $\kk$$ point.
In the best case scenario, both the CBM and the VBM are located at the $\Gamma$ point; hence the BZ summation
can be replaced by a much faster symmetrized sum over the wavevectors of the IBZ.


## GWR workflow for QP energies

A typical GWR workflow with arrows denoting dependencies between the different steps
is schematically represented in the figure below:

![](eph_intro_assets/eph_workflow.png){: style="height:400px;width:400px"}

In order to perform a standard one-shot GW calculation, one has to:

  1. Run a converged GS calculation to obtain the density.

  2. Perform a NSCF run to compute the KS eigenvalues and eigenfunctions
     including several empty states. Note that, unlike standard band structure calculations,
     here the KS states *must* be computed on a regular grid of **k**-points.

  3. Use [[optdriver]] = 3 to compute the independent-particle susceptibility $\chi^0$ on a regular grid of
     **q**-points, for at least two frequencies (usually, $\omega=0$ and a purely imaginary
     frequency - of the order of the plasmon frequency, a dozen of eV).

The Green's function is constructed from the KS wavefunctions and eigenvalues stored in the WFK file.

The following physical properties can be computed:

* Imaginary part of ph-e self-energy in metals (**eph_task 1**) that gives access to:

    * Phonon linewidths induced by e-ph coupling

* Real and imaginary parts of the e-ph self-energy (**eph_task 4**) that gives access to:

    * Zero-point renormalization of the band gap

## Tricks to accelerate the computation and reduce the memory requirements

    [[gwr_boxcutmin]]

## Self-consistency with GWR

The conventional GW code supports different kinds of self-consistency including the update
of the wavefunctions and the QPSCGW method proposed in [[cite:Faleev2004]].
The GWR code, on the contrary, only supports self-consistency on energies.

Performing energy-only self-consistent calculations with GWR is much easier than in the conventional GW code
as there is no need to chain different screening/sigma calculations together as the self-consistent loop
and the stopping criterion are implemented directly at the Fortran level.
To perform an energy-only self-consistent calculation with GWR, one has to specify the kind of self-consistency
via [[gwr_task]], the maximum number of iterations with [[gwr_nstep]] and the stopping criterion with [[gwr_tolqpe]].
Note that, at present, it is not possible to restart a self-consistent calculation from a previous checkpoint.

### Requirements

The GWR code requires an ABINIT build with Scalapack enabled.
Moreover, a significant fraction of the computing time is spent in performing FFTs thus we **strongly**
recommend to use vendor-optimized FFT libraries such as MKL-DFTI or FFTW3 instead
of the internal FFT version shipped with ABINIT.

Note that single-precision is the default mode as in the conventional $GW$ code.
To run computations in double-precision, one has to configure with  `--enable-gw-dpc="yes"` when
the command line interface is used or `enable_gw_dpc="yes"` when `--with-config-file=FILE` is used
to specify the configuration options via an external FILE.

### MPI parallelization

The GWR code employs a 4D MPI Cartesian grid to distribute both workload and memory over
collinear spins, points of the minimax mesh, $(\gg, \gg')$ components, and $\kk$-points.

The user can specify the number of MPI-processes for the different dimensions using the input variable [[gwr_np_kgts]],
although this is completely optional as GWR is able to build the MPI Cartesian grid at runtime on the basis
of the number of MPI processes available.

There are, however, some basic rules that is worth keeping in mind if good parallel performance is wanted.
Ideally the total number of MPI processes should be a multiple of [[gwr_ntau]] * [[nsppol]] as the parallelism
over minimax points and spin is rather efficient (few communications required).
On the contrary, the parallelism over $\gg$-vectors and $\kk$-points is much more network intensive, although
these two leves allow one to decrease the memory requirements.
