---
authors: MG
---

# An overview of the GWR code

This page provides a quick introduction to the new GWR driver of ABINIT
We discuss the technical details related to the implementation and the associated input variables.
The drawbacks/advantages with respect to the legacy GW implementation 
in Fourier-space and frequency domain are also discussed.

## Why a new GWR code?

The conventional GW algorithm has quartic scaling with the number of atoms whereas GWR scales cubically.
The legacy GW code 

Select the task to be performed when [[optdriver]] == 6 i.e. GWR code.
while [[gwr_task]] defines the task

### Requirements

* Scalapack
* Optmized FFT libraries (FFTW3 or MKL-DFTI)

Discuss single and double precision version. Single-precision is the default

## Formalism

The zero-temperature Green's function in the imaginary-time domain is given by:

\begin{equation}
G(\rr, \rr', i\tau) = 
\Theta(\tau) \Gove(\rr, \rr', i\tau) +
\Theta(-\tau) \Gund (\rr, \rr', i\tau) 
\end{equation}

with $\Theta$ the Heaviside step-function and 

\begin{equation}
\Gove(\rr, \rr', i\tau) =
-\sum_n^{\text{unocc}} \psi_n(\rr)\psi_n^*(\rr') e^{-\ee_n\tau}
\qquad (\tau > 0)
\end{equation}

\begin{equation}
\Gund(\rr, \rr', i\tau) = 
\sum_n^{\text{occ}} \psi_n(\rr)\psi_n^*(\rr') e^{-\ee_n\tau}
\qquad (\tau < 0)
\end{equation}

For simplicity we have assumed a spin-unpolarized semiconductor with scalar wavefunctions
and the Fermi level $\mu$ is set to zero.

The GWR code constructs the Green's function from the KS wavefunctions and eigenvalues stored 
in the WFK file specified via [[getwfk_filepath]] ([[getwfk]] in multi-dataset mode).
This WFK file is usually produced by performing a NSCF run with many empty states and the GWR driver
provides a specialized option to perform a direct diagonalization of the KS Hamiltonian with Scalapack.
(see section below)
When computing $G$, the number of bands included in the sum over states is controlled by [[nband]].
Clearly, it does not make any sense to ask for more bands than the ones available in the WFK file. 

The imaginary axis is sampled using a minimax mesh with [[gwr_ntau]] points.
The other piece of information required for the selection of the minimax mesh 
is the ratio between the fundamental gap and the maximum transition energy i.e. 
the differerence between the highest KS eigenvalue for empty states that clearly depends on [[nband]] and
the lowest occupied state.

The k-mesh must be k-centered e.g.

* [[ngkpt]] 4 4 4
* [[nshiftk]] 1 
* [[shiftk]] 0 0 0 

The cutoff energy for the polarizability is given by [[ecuteps]]
while [[ecutsigx]] defines the number of g-vectors for the exchange part of the self-energy.

Note that GWR needs the GS density produced during the GS run.
This file can be read via [[getden_filepath]] (recommended) or, alternatively 
with [[getden]] when using the multi dataset mode.

\begin{equation}
\chi(\rr,\RR', t) = G(\rr,\RR', i\tau) G^*(\rr,\RR', -i\tau)
\end{equation}

Treatment of the long-wavelenght limit

[[inclvkb]], [[gw_qlwl]], [[gwr_max_hwtene]]

## GWR workflow for QP energies

A typical GWR workflow with arrows denoting dependencies between the different steps
is schematically represented in the figure below:

![](eph_intro_assets/eph_workflow.png){: style="height:400px;width:400px"}

In order to perform a standard one-shot GW calculation one has to:

  1. Run a converged Ground State calculation to obtain the self-consistent density.

  2. Perform a non self-consistent run to compute the KS eigenvalues and the eigenfunctions
     including several empty states. Note that, unlike standard band structure calculations,
     here the KS states *must* be computed on a regular grid of **k**-points.

  3. Use [[optdriver]] = 3 to compute the independent-particle susceptibility $\chi^0$ on a regular grid of
     **q**-points, for at least two frequencies (usually, $\omega=0$ and a purely imaginary
     frequency - of the order of the plasmon frequency, a dozen of eV).

The Green's function is constructed from the KS wavefunctions and eigenvalues stored in the WFK file
hence the 

The following physical properties can be computed:

* Imaginary part of ph-e self-energy in metals (**eph_task 1**) that gives access to:

    * Phonon linewidths induced by e-ph coupling

* Real and imaginary parts of the e-ph self-energy (**eph_task 4**) that gives access to:

    * Zero-point renormalization of the band gap

At the time of writing, the following features are **not yet supported** in GWR:

* PAW calculations
* Spinors ([[nspinor]] = 2)

## Tricks to accelerate the computation and reduce the memory requirements
