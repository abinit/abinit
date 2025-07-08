---
authors: MG
---

# An overview of the GWR code

This page offers a concise introduction to the new GWR driver of ABINIT.
We outline the technical details of its implementation, the relevant input variables,
and compare its advantages and limitations against the conventional GW approach formulated
in Fourier space and real-frequency domain.
In the following, we will refer to this traditional approach as the **conventional** or **standard** GW code of ABINIT.

## Why a new GW code?

In the standard GW code, one uses the Lehmann representation of the Green's function in the real-frequency domain
to compute the irreducible polarizability, as explained in the [MBPT notes](/theory/mbt).
Then the self-energy matrix elements in the KS representation are computed via an expensive convolution
in the frequency domain [[cite:Golze2019]] — by default using the plasmon-pole approximation [[cite:Giantomassi2011]],
which significantly accelerates calculations but introduces approximations and prevents direct access
to the spectral function $A_\nk(\omega)$.
The conventional GW code exhibits quartic scaling with the number of atoms in the unit cell
and quadratic scaling with the number of $\kk$-points.
Moreover the computation of the polarizability with the Adler-Wiser expression represents
the most CPU-demanding part due to the double sum over conduction and valence states.

In contrast, the GWR code achieves cubic scaling with [[natom]] and linear scaling in the number of $\kk$-points
by computing the polarizability and the self-energy in the real-space supercell associated to the $\kk$-mesh [[cite:Liu2016]].
FFTs are used to transform quantities from the supercell to Fourier space whenever needed.
For instance, $W$ (a convolution between $\ee^-1$ and the bare Coulomb interaction) is
best solved in Fourier space and frequency domain whereas the polarizability and the self-energy are best evaluated
in real-space and (imaginary) time.
As concerns the frequency dependence, GWR evaluates $G$, $W$ and $\Sigma$ along the imaginary axis
using minimax meshes that are adaptive grids that place points non-uniformly, concentrating them where they are most needed,
leading to higher accuracy with fewer points [[cite:Kaltak2014]].
Inhomogeneous sine/cosine transforms are then used to from imaginary time to imaginary frequency with quadratic scaling
in the number of imaginary points. For this reason, it is very important to keep the number of points as small possible
without spoiling accuracy [[cite:Liu2016]].
The minimax meshes used in ABINIT are provided by the GreenX library ([[cite:Azizi2023]], [[cite:Azizi2024]]).

In the last step of a GWR calculation, the self-energy in imaginary time is transformed to imaginary frequency
followed by an analytic continuation (AC) to the real-frequency axis where the QP corrections are obtained
by solving the linearized version of the QP equation.

\begin{equation} \label{eq:implicit_QP_energy}
\ee^\QP = \ee^\KS + Z \langle\Psi^\KS|\Sigma(\ee^\KS)-\vxc|\Psi^\KS\rangle.
\end{equation}

where

\begin{equation} \label{eq:Z_factor}
Z \equiv \left[ 1- \langle \Psi^\KS|
\PDER{\Sigma}{\ee^\KS}|\Psi^\KS\rangle \right]^{-1}
\end{equation}

is the so-called renormalization factor where the derivative of the self-energy is computed at the KS energy.

The AC enables access to the full frequency dependence of $\Sigma_\nk$ and $A_\nk$
at a substantially reduced computational cost with respect to e.g. the countour deformation method used
in the conventional approach in which several points are needed to accurately capture the sharp features
in the Green's function $G(\omega)$, and the screened interaction $W(\omega)$ along the real-frequency axis,
The price to pay is that in GWR the accuracy of the results now depends on the stability of the AC procedure.

!!! important

    Although GWR exhibits better scaling with [[natom]], it is worth noticing that the standard implementation
    makes better use of symmetries [[cite:Aulbur2001]] to reduce the number of points
    in the BZ summation (see [[symchi]] and [[symsigma]])
    Consequently, the standard code is still competitive and can be faster than GWR
    in the case of high symmetric systems.

## Input variables

To enter the GWR driver of ABINIT, one has to use [[optdriver]] = 6,
and then select the computation to be performed via [[gwr_task]].

In the first step, you will very likely use [[gwr_task]] = "HDIAGO" or "HDIAGO_FULL"
to perform a **direct diagonalization** with ScaLAPACK or ELPA in order to generate a WFK with empty states and cutoff [[ecut]].
This feature can also be used if you want to use other many-body codes that are able to read ABINIT’s WFK file
Please take some time to read the documentation of [[gwr_task]] to see which calculations are supported.

The $\kk$-points and the band range for the QP corrections in $\Sigma_\nk$ can be specified in different ways.
In the explicit approach, one uses [[nkptgw]] to specify the number of $\kk$-points in $\Sigma_\nk$,
In this case, [[kptgw]] determines the list of $\kk$-points in reduced coordinates, while [[bdgw]] specifies the band range
for each $\kk$-point.
Alternatively, one can use [[gw_qprange]] to let ABINIT automatically select the wavevectors and the band range.
This is particularly useful for self-consistent calculations or if all the QP corrections in the IBZ are wanted.

As concerns the computation of the spectral function $A_\nk(\ww)$, [[nfreqsp]], [[freqspmin]], [[freqspmax]] can be used
to specify the frequency mesh along the real-axis.
It is worth mentioning that in the conventional GW code,
the computational cost increases quickly with [[nfreqsp]] while in the GWR algorithm
this part is very cheap as $A_\nk(\ww)$ is easily obtained with the Pade' approximant.

IF for some reasons, you decide to use iterative eigenvalue solvers to generate the WFK file,
remember to set [[istwfk]] = 1

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



In Fourier space, one obtains:

\begin{equation}
\Gove_\kk(\bg, \bg', i\tau) =
-\sum_n^{\text{unocc}} u_\nk(\bg) u_\nk^*(\bg') e^{-\ee_\nk\tau}
\qquad (\tau > 0)
\end{equation}

where the number of PWs is controlled by [[ecutwfn]] that by default is equal to [[ecut]].
These matrices are stored in memory for the entire duration of the calculation.
The $\bg, \bg'$ components are distributed inside the `g` level of [[gwr_np_kgts]].
Only the $\kk$-points in the IBZ are stored in memory and distributed inside the $\kk$-level
as the code can use symmetries to reconstruct the Green's function in the full on-the-fly.
Finally, the $\tau$ points are distributed inside the `t` level and collinear spins are distributed inside the `s` level.

!!! important

    For simplicity, in all the equations we assume a spin-unpolarized semiconductor with scalar wavefunctions
    (i.e. [[nsppol]] = 1 with [[nspinor]] = 1) and the Fermi level $\mu$ is set to zero.
    At the level of the implementation, this means that the initial KS eigenvalues produced by the NSCF part of ABINIT
    are shifted with respect to the value of $\mu$ at zero temperature that is now located at mid-gap.

The GWR code constructs the Green's function from the KS wavefunctions and eigenvalues stored
in the WFK file specified via [[getwfk_filepath]] ([[getwfk]] if you really insist in using multi-dataset mode).
This WFK file is usually produced by performing a NSCF calculation including empty states and the GWR driver
provides a specialized option to perform a **direct diagonalization** of the KS Hamiltonian in parallel
with Scalapack.

When computing the KS Green's function $G^0$, the number of bands is controlled by [[nband]].
Clearly, it does not make any sense to ask for more bands than the ones available in the WFK file.
Note that GWR also needs the GS density produced by a previous GS SCF run in order to compute
the matrix elements of $v_\xc[n]$ between KS states.
The location of the density file can be specified via [[getden_filepath]]
([[getden]] if you really insist in using multi-dataset mode).

The imaginary axis is sampled using a minimax mesh with [[gwr_ntau]] points.
The other piece of information required for the selection of the minimax mesh
is the ratio between the **fundamental** gap and the maximum transition energy i.e.
the differerence between the highest KS eigenvalue for the empty states that, in turns, depends
on the value of [[nband]] and the energy of the lowest occupied state.

The irreducible polarizability is computed in the real-space supercell using

\begin{equation}
\chi(\rr,\RR', i\tau) = G(\rr, \RR', i\tau) G^*(\rr, \RR', -i\tau)
\end{equation}

and then immediately transformed to Fourier space to obtain $\chi_\qq(\bg, \bg', i\tau)$.
Here, we use $\rr$ to indicate points in the unit cell, and $\RR$ to denote points in the supercell.

The cutoff energy for the polarizability is given by [[ecuteps]].
This is an important parameter that should be subject to convergence studies.
Only the $\qq$-points in the IBZ are stored in memory and distributed inside the k-level of [[gwr_np_kgts]].

The FFT mesh is defined according to [[gwr_boxcutmin]].
This parameter has a **big impact** on wall-time and memory requirements.
One usually start with a coarse FFT mesh e.g. [[gwr_boxcutmin]] 1.0 and
increase it during the last steps of the convergence study.

As discussed in [[cite:Baroni1986]], the treatment of the long-wavelenght limit $\qq \rightarrow 0$
in the polarizability, requires the inclusion of the commutator of the non-local part of the Hamiltonian
with the position operator.
See also [this section](/theory/bse/#5-matrix-elements-of-the-dipole-operator) of the Bethe-Salpeter notes.
By default the commutator of the non-local part is always included.
The input variable [[inclvkb]] can be used to deactivate it, if needed.

<!--
[[gw_nqlwl]],
[[gw_qlwl]],
[[gwr_max_hwtene]]
-->

The exchange part of $\Sigma$ is computed using the standard sum over occupied states
similarly to what is done in the conventional code:


<!--
\begin{equation}\label{eq:Sigma_x}
\Sigma_x(\rr_1,\rr_2)= -\sum_\kk^\BZ
\sum_\nu^\text{occ} \Psi_{n\kk}(\rr_1){\Psi^\*_{n\kk}}(\rr_2)\,v(\rr_1,\rr_2)
\end{equation}
-->

\begin{equation}
\label{eq:diag_mat_el_Sigma_x}
 \langle b_1\kk|\Sigma_x|b_1\kk\rangle =
 -\dfrac{4\pi}{V} \sum_{\nu}^{\text{occ}} \sum_\qq^\BZ \sum_{\GG}
  \dfrac{|M_\GG^{b_1b_1}(\kk,\qq)|^2}{|\qpG|^2}.
\end{equation}



The number of $\bg$-vectors in the sum is defined by [[ecutsigx]].
Note that computing $\Sigma_x$ is much cheaper than $\Sigma_c(i \tau)$.

## Real-space vs convolutions in the BZ

So far we have discussed the GWR equations used to compute the polarizability
and the self-energy in the real-space supercell.
This is the recommended approach if one needs to compute QP corrections for all the $\kk$-points in the IBZ as
we have linear scaling in the number of wave vectors.
This is the typical scenario for self-consistent GWR calculations or if one needs to perform some kind of interpolation
of the QP corrections to obtain a band structure along a high-symmetry $\kk$-path starting from the values
of $\Sigma_\nk$ for all $\kk$-points in the IBZ.

It should be noted, however, that in several applications one is mainly interested in the QP corrections
at the band edges that are usually located at high-symmetry $\kk$-points.
In this case, it is more advantageous to use an alternative formulation in which $\Sigma_\nk$
is expressed in terms of convolutions in the BZ according to:

\begin{equation}
\Sigma_\kk(\rr, \rr',i\tau) =
\sum_\qq
{G}_{\kq}(\rr,\rr',i\tau)
{W}_\qq(\rr',\rr, i\tau)
\end{equation}

where $G_\qq(\rr,\rr')$ and $W_\qq(\rr,\rr')$ are now quantities defined in the unit cell.

<!--
\begin{equation}
\tilde \chi_\qq(\rr, \rr',i\tau) =
\sum_\kk
{\tilde G}_{\kq}(\rr,\rr',i\tau)
{\tilde G}_\kk(\rr',\rr,-i\tau)
\end{equation}

Let us compute the RPA polarizabilty in real space using the mixed-space approach (the time dependence is not shown for the sake of simplicity)

\begin{equation}
G(\rr,\rr',i\tau) =
\sum_\kk e^{i\kk(\rr-\rr')}
{\tilde G}_\kk(\rr,\rr',i\tau)
\end{equation}

\begin{equation}
\chi(\rr,\rr',i\tau) =
\sum_\kk e^{i\kk(\rr - \rr')} {\tilde G}_\kk(\rr,\rr',i\tau)
\sum_{\kk'} e^{i\kk'(\rr' - \rr)} {\tilde G}_{\kk'}(\rr',\rr,-i\tau) =
\sum_\kk
{\tilde G}_{\kq}(\rr,\rr',i\tau)
{\tilde G}_\kk(\rr',\rr,-i\tau)
\end{equation}

hence

\begin{equation}
\tilde \chi_\qq(\rr, \rr',i\tau) =
\sum_\kk
{\tilde G}_{\kq}(\rr,\rr',i\tau)
{\tilde G}_\kk(\rr',\rr,-i\tau)
\end{equation}

The advantage is that FFTs are performed in the unit cell and symmetries are easier to exploit.
The disadvantage is that products in real space become convolutions in $\kk$-space, hence the computation is quadratic in
the number of $\qq$-points with a prefactor that depends on the symmetries of the system.
-->

The advantage of this formulation is that FFTs are performed in the unit cell
and one can take advantage of the symmetries of the system
to reduce the BZ integration to the irreducible wedge defined by the little group of the $\kk$ point in $\Sigma_\nk$.
In the best case scenario of a direct band gap semiconductor, both the CBM and the VBM are located at the $\Gamma$ point;
hence the BZ integral can be replaced by a much faster symmetrized integral over the wavevectors of the IBZ.
Clearly this approach is not the most efficient one if all the $\kk$-points are wanted
as the convolutions in the BZ lead to quadratic scaling with the number of $\kk$-points.

The input variable [[gwr_sigma_algo]] allows one to select the algorithm to be used.
[[gwr_chi_algo]] has a similar meaning but only for the polarizability .
This option renders the computation of the polarizability slower although it requires less memory
as only two-point functions in the unit cell need to be stored.

## Self-consistency with GWR

The conventional GW code supports different kinds of self-consistency, including the update
of the wavefunctions and the QPSCGW method proposed in [[cite:Faleev2004]].
On the contrary, at the time of writing, the GWR code only supports self-consistency on energies.

On the bright side, performing energy-only self-consistent calculations with GWR is much easier than in the conventional GW code
as there is no need to chain different screening/sigma calculations together: the self-consistent loop
and the stopping criterion are implemented directly at the Fortran level.
To perform an energy-only self-consistent calculation with GWR, one has to specify the kind of self-consistency
via [[gwr_task]], the maximum number of iterations with [[gwr_nstep]] and the stopping criterion with [[gwr_tolqpe]].
Note that, at present, it is not possible to restart a self-consistent calculation from a previous checkpoint file.

### MPI parallelization and performance

The two codes strongly differ at the level of the MPI parallelization.
In the standard GW code, MPI parallelization is available over the [[nband]] states and [[nsppol]]; however,
key data structures, such as $W$, are not MPI-distributed.
As a result, the maximum number of usable MPI processes is limited by [[nband]],
and the workload becomes imbalanced when the number of MPI processes does not evenly divide [[nband]].
More critically, $\Sigma$ calculations in the standard GW code are highly memory-intensive,
as they require storing both the wavefunctions (whose memory footprint scales with the number of MPI processes)
as well as $W$ (non-scalable portion).
This memory requirement becomes particularly problematic for large [[npweps]] values or calculations
beyond the plasmon-pole approximation, where the full $W$ matrix must be stored for multiple frequencies.
Consequently, conventional GW calculations can be prohibitively demanding in terms of memory,
especially for large systems or systems with few symmetries.

In contrast, the GWR code distributes most data structures across MPI processes,
which enables handling larger systems when sufficient computing nodes are available.
However, this distribution comes at the cost of increased MPI communication in certain parts of the algorithm.
More specifically, GWR leverages Parallel BLAS (PBLAS) to efficiently distribute the memory required
for storing the Green’s functions and $W$, significantly improving scalability compared to the standard implementation.

The GWR code employs a 4D MPI Cartesian grid to distribute both workload and memory over
collinear spins ([[nsppol]] = 2), points of the minimax mesh, $(\gg, \gg')$ components, and $\kk$-points.
The user can specify the number of MPI-processes for the different dimensions using the input variable [[gwr_np_kgts]],
although this is completely optional as GWR is able to build the MPI Cartesian grid at runtime on the basis
of the number of MPI processes available, the size of the problem and the memory allocated.

There are, however, some basic rules that are worth keeping in mind if decent parallel performance is wanted.
Ideally, the total number of MPI processes should divide (possibly a multiple of) [[gwr_ntau]] * [[nsppol]] as the parallelism
over minimax points and spins is rather efficient (few communications required).
On the contrary, the parallelism over $\gg$-vectors and $\kk$-points is much more network intensive,
although these two additional levels allow one to decrease the memory requirements required to store the two-point functions.

In case out-of-memory problems occur, we recommend increasing as much as possible the number of CPUs
used to distribute the matrices with PBLAS (g-level).
Once the memory usage for $\chi$ and $W$ has dropped to an acceptable level, one can start increasing
the number of MPI processes for \tau-parallelism.

<!--
TODO: Additional trick to be tested/documented in more detail:
[[gwr_ucsc_batch]]
[[gwr_rpa_ncut]]
[[gwr_boxcutmin]]
[[gwr_max_hwtene]]
gwr_regterm
-->

### Requirements

The GWR code requires an ABINIT build with Scalapack enabled.
Moreover, since a significant fraction of the computing time is spent in performing FFTs,
we **strongly** recommend to use vendor-optimized FFT libraries such as MKL-DFTI or FFTW3 instead
of the internal FFT version shipped with ABINIT ([[fftalg]] should be 312 or 512, this is done automatically
if an external FFT library is found by the build system).

If you are using the MKL library by intel, the configuration is relatively easy
as MKL provides all the required libraries. An example of `.ac` configuration file is reported below:

```
# BLAS/LAPACK with MKL
with_linalg_flavor="mkl"
LINALG_FCFLAGS="-I${MKLROOT}/include"
LINALG_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl"

# FFT from MKL
with_fft_flavor="dfti"
```

Depending on your version of the intel compiler and MKL, you may need to customize the list of libraries given in this example.
Please use the [link-line-advisor](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html)
provided by intel and select the option: "Link with Intel® oneMKL libraries explicitly" to get the value of `LINALG_LIBS`.

Note that single-precision is the default mode as in the conventional $GW$ code.
To run computations in double-precision, one has to configure the package with  `--enable-gw-dpc="yes"`
when the command line interface is used or `enable_gw_dpc="yes"` when `--with-config-file=FILE` is used
to pass the configuration options to `configure` via an external FILE.

To configure with ELPA + MKL, set the ELPAROOT env variable with :

```
export ELPAROOT="PATH_TO_ELPA_INSTALLATION_DIR"
```

and then

```
with_linalg_flavor="mkl+elpa"
LINALG_FCFLAGS="-I\${MKLROOT}/include -I\${ELPAROOT}/modules"
LINALG_LIBS="-L\${ELPAROOT}/lib -lelpa -L${MKLROOT}/lib -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl"
```

At the end of configure, you should obtain:

```
  * FFT flavor        : dfti (libs: user-defined)
  * LINALG flavor     : mkl+elpa (libs: user-defined)
  * SCALAPACK enabled : yes
  * ELPA enabled      : yes
```

!!! tip

    To run all the tests associated to the GWR code with `NUM` MPI processes
    and `PYTASKS` python multiprocessing, use:

    ```
    runtests.py -k GWR -n NUM -j PYTASKS
    ```
