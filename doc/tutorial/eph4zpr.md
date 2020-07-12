---
authors: MG
---

# Zero-point renormalization of the band gap and temperature-dependent band structures

This tutorial explains how to compute the electron self-energy due to phonons, evaluate the zero-point
renormalization (ZPR) of the band gap, obtain spectral functions and temperature-dependent band structures.
We start with a very brief overview of the many-body formalism in the context of the electron-phonon (e-ph) interaction,
then we discuss how to compute the e-ph self-energy and perform typical convergence studies.

It is assumed the user has already completed the two tutorials [RF1](rf1) and [RF2](rf2),
and that he/she is familiar with the calculation of ground state and response properties,
in particular phonons, Born effective charges and dielectric tensor.
It goes without saying that one should have read the [introduction page for the EPH code](eph_intro)
before running these examples.

This lesson should take about 1.5 hour.

## Formalism

The electron-phonon self-energy, $\Sigma^{\text{e-ph}}$, describes the renormalization of the single-particle
energy due to the interaction with phonons.
This term should be added to the electron-electron (e-e) self-energy $\Sigma^{\text{e-e}}$
that encodes many-body effects induced by the Coulomb interaction beyond the Hartree potential.
The e-e contribution can be estimated using, for instance, the $GW$ approximation but
in this tutorial we are mainly interested in $\Sigma^{\text{e-ph}}$ and its temperature dependence.

In semiconductors and insulators, indeed, most of temperature dependence of the electronic properties
at low T originates from the e-ph interaction and the thermal expansion of the unit cell.
Corrections due to $\Sigma^{\text{e-e}}$ are obviously important as KS gaps computed at the LDA/GGA level
are systematically underestimated with respect to experiments but the temperature dependence of $\Sigma^{\text{e-e}}$
is rather small as long as the fundamental gap is much larger than $kT$.

In state-of-the-art *ab-initio* methods, the e-ph coupling is described
within DFT by expanding the KS effective potential in the nuclear displacements,
and the vibrational properties are obtained with DFPT [[cite:Gonze1997]], [[cite:Baroni2001]].
The e-ph self-energy consists of two terms: Fan-Migdal (FM) and the Debye-Waller (DW) [[cite:Giustino2017]]

$$
\Sigma^\eph(\ww) = \Sigma^\FM(\ww) + \Sigma^{\DW}.
$$

The diagonal matrix elements of the FM self-energy in the KS representation are given by

\begin{equation}
\begin{split}
    \Sigma^\FM_{n\kk}(\omega,\ef,T) =
                & \sum_{m,\nu} \int_\BZ \frac{d\qq}{\Omega_\BZ} |\gkq|^2 \\
                & \times \left[
                    \frac{n_\qnu(T) + f_{m\kk+\qq}(\ef,T)}
                         {\omega - \emkq  + \wqnu + i \eta} \right.\\
                & \left. +
                    \frac{n_\qnu(T) + 1 - f_{m\kk+\qq}(\ef,T)}
                         {\omega - \emkq  - \wqnu + i \eta} \right] ,
\end{split}
\label{eq:fan_selfen}
\end{equation}

where $f_{m\kk+\qq}(\ef,T)$ and $n_\qnu(T)$ are the Fermi-Dirac and Bose-Einstein occupation functions
with $T$ the temperature and $\ef$ the Fermi level.
<!--
For the sake of simplicity, the temperature and Fermi level are considered as parameters, and the dependence
on $T$ and $\ef$ will be omitted in the following.
-->
The integration is performed over the $\qq$-points in the BZ of volume $\Omega_\BZ$ and $\eta$
is a positive real infinitesimal.

!!! important

    From a mathematical point of view, one should take the limit $\eta \rightarrow 0^+$.
    At the level of the implementation, the infinitesimal $\eta$ is replaced by a (small) 
    finite value given by the [[zcut]] variable that should be subject to convergence studies.
    More specifically, one should monitor the convergence of the physical properties of interest
    as a function of [[zcut]] and $\qq$-point sampling similarly to what is done 
    in metals for [[tsmear]] and [[ngkpt]].

<!--
First-principles calculations of the EPH self-energy are therefore crucial to understand the temperature-dependence
of band gaps, including the correction due to zero-point motion, as well as for computing phonon-limited mobilities
within the Boltzmann transport equation.
In ABINIT v9, it will be possible to compute the EPH self-energy in the Kohn-Sham representation using the EPH matrix elements.
The code employs optimized algorithms to compute either the full self-energy
(needed for QP corrections and spectral functions)
When computing the full self-energy, it is possible to reduce the number of empty states required for convergence
by using the first-order wavefunctions obtained by solving the relevant Sternheimer equation.
-->

The **static DW term** involves the second order derivative of the KS potential with respect to the nuclear displacements.
State-of-the-art implementations approximate the DW contribution with

\begin{equation}
\label{eq:dw_selfen}
\Sigma_{n\kk}^{\DW} = \sum_{\qq\nu m} (2 n_{\qq\nu} + 1) \dfrac{g_{mn\nu}^{2,DW}(\kk, \qq)}{\ee_{n\kk} - \ee_{m\kk}},
\end{equation}

where $g_{mn\nu}^{2,\DW}(\kk,\qq)$ is an effective matrix element that, within the **rigid-ion approximation**,
can be expressed in terms of the $\gkq$ matrix elements using the invariance of the QP energies 
under infinitesimal translation [[cite:Giustino2017]].

At the level of the implementation, the number of bands in the two expressions is defined by [[nband]]
while the $\qq$-mesh is specified by [[eph_ngqpt_fine]].
The potentials are interpolated if [[eph_ngqpt_fine]] differs from [[ddb_ngqpt]].

!!! tips

    The EPH code takes advantage of symmetries to reduce the BZ integration to an appropriate 
    irreducible wedge, $\text{IBZ}_k$, defined by the little group of the $\kk$-point.
    Calculations for high-symmetry $\kk$-points such as $\Gamma$ are therefore much faster as there 
    are more symmetries that leave the $\kk$-point unchanged within a reciprocal lattice vector.

    This symmetrization procedure is activated by default and can be deactivated by setting [[symsigma]]
    to zero for testing purposes.
    Note that when [[symsigma]] is set to 1, the code performs a final average of the QP results
    within each degenerate subspace.
    As a consequence, accidental degeneracies won't be removed when [[symsigma]] = 1.

<!--
Both the FM and the DW term converge slowly with respect to the number of empty states and the $\qq$-sampling.
Note that, in principle, one should sum an infinite number of states (read from the external WFK file)
and that the convergence is rather slow.
Note also that QP energy differences usually converge faster than QP energies (the same behaviour is observed in $GW$).
The EPH code implements numerical tricks to accelerate the convergence with respect to [[nband]] that are discussed
in more detail in the next sections.
Also, the convergence with the $\qq$-sampling is rather slow and delicate.
This is especially true in polar materials due to the divergence of the polar e-ph matrix elements
in the $\qq \rightarrow 0$ limit [[cite:Vogl1976]].
-->

### Quasi-particle corrections due to e-ph coupling

Strictly speaking, the quasi-particle (QP) excitations are defined by the solution(s)
in the complex plane of the equation $$ z = \ee_\nk + \Sigma_\nk^{\text{e-ph}}(z) $$
provided that the non-diagonal components of the self-energy can be neglected.
In practice, the problem is usually simplified by seeking approximated solutions along the real axis
following two different approaches: **on-the-mass-shell** and **linearized QP equation**.

In the on-the-mass-shell approximation, the QP energy is given by the real part
of the self-energy evaluated at the bare KS eigenvalue:

\begin{equation}
\ee^\QP_\nk = \ee_\nk + \Re\, \Sigma_\nk^{\text{e-ph}}(\ee_\nk).
\end{equation}

This approach is equivalent to standard time-dependent Rayleigh-Schrodinger perturbation theory.
In the linearized QP equation, on the contrary, the self-energy is Taylor-expanded around 
the KS eigenvalue and the QP correction is obtained using

\begin{equation}
  \ee^\QP_\nk = \ee_\nk
      + Z_\nk\,\Re  \Sigma^\text{e-ph}_\nk(\ee_\nk)
\end{equation}

with the renormalization factor $Z_\nk$ given by

\begin{equation}
  Z_\nk= \left(1 - \Re\left[ \frac{\partial\Sigma^\text{e-ph}_{\nk}}{\partial\ee}\right]\Bigg|_{\ee=\ee_\nk} \right)^{-1}.
\end{equation}

Both approaches are implemented in ABINIT although it should be noted that, according to recent works,
the on-the-mass-shell approach provides results that are closer to those obtained
with more advanced techniques based on the cumulant expansion [[cite:Nery2018]].

!!! important

    The EPH code can compute QP corrections only for $\nk$ states that are present in the input WFK file
    (a similar requirement is present in the $GW$ code as well).
    As a consequence, the $\kk$-mesh ([[ngkpt]], [[nshiftk]], [[shiftk]]) should be chosen carefully 
    especially tf the band edge is not located at an high-symmetry $\kk$-point.

There are different approaches one can use to specify the set of $\nk$ states in $\Sigma_{\nk}$.
Each approach has pros and cons.

The most direct way consists in specifying explicitly the $\kk$-points and the band range 
using the three variables: [[nkptgw]], [[kptgw]], [[bdgw]]
To compute the correction for the VBM at $\Gamma$ in a non-magnetic semiconductor
with 8 electrons per unit cell, one would use:

```sh
nkptgw 1
kptgw  0 0 0  # [3, nkptgw] array
bdgw   4 4    # [2, nkptgw] arary giving the initial and the last band index
              # for each nkptgw k-point
```

as the index of the valence band is given by 8 / 2 = 4.
When [[symsigma]] is set to 1 (default), the code may decide to enlarge the band range 
so that all degenerate states for a particular $\kk$-point are included in the calculation.

Alternatively, one can use [[gw_qprange]] or [[sigma_erange]]

### Spectral function and Eliashberg functions

For the spectral function, we use [[nfreqsp]], [[freqspmax]]

<!--
Further details concerning the implementation are given in Ref.[[cite:Gonze2019]].
In order to accelerate the convergence with the number of empty states, ABINIT replaces the contributions
given by the high-energy states above a certain band index $M$ with the solution
of a non-self-consistent Sternheimer equation in which only the first $M$ states are required.
The methodology, proposed in [[cite:Gonze2011]], is based on a quasi-static approximation
in which the phonon frequencies in the denominator of Eq.~(\ref{eq:fan_selfen}) are neglected and
the frequency dependence of $\Sigma$ is approximated with the value computed at $\omega = \enk$.
This approximation is justified when the bands above $M$ are sufficiently high in energy with respect
to the $n\kk$ states that must be corrected.
The code can compute the spectral functions
-->

## Typical workflow

A typical workflow for ZPR and T-dependent calculations requires the following steps:

1. **GS calculation** to obtain the WFK and the DEN file.
   Remember to set [[prtpot]] to 1 to produce the file with the KS potential required by the Sternheimer method.

2. **DFPT calculations** for all the IBZ $\qq$-points corresponding to the *ab-initio* [[ddb_ngqpt]] mesh
   used to perform the Fourier interpolation of the dynamical matrix and of the DFPT potentials.
   Remember to compute epsinf, BECS and Q*

3. **NSCF computation** of a WFK file on a dense $\kk$-mesh containing the wavevectors
   for which phonon-induced QP corrections are wanted. This part uses the DEN file produced in step #1.
   Remember to include enough empty states so that it is possible to perform convergence studies wrt [[nband]].

4. **Merge the partial DDB and POT files** with *mrgddb* and *mrgdvdb*, respectively

5. Use the full DDB file, the DVDB file and the WFK file obtained in step #3 to perform ZPR calculations
   with [[eph_task]] 4.

## Calculation of the ZPR for XXX

In this tutorial, we prefer to focus on the usage of the EPH code hence
we will be using pre-computed DDB and POT1 files obtained with a xxx $\qq$-mesh to bypass the DFPT computation
We also provide a DEN.nc file that can be used to start the NSCF calculations required
to generate the WKF file.

First of all, let's merge the DDB files using the following input file:

Now we can merge the DFPT potential with the *mrgdv* tool and the following input file:

For our first NSCF calculation, we use a xxx $\kk$-mesh and XXX bands.
The number of bands is sufficiently large so that we can perform initial convergence studies.
We also perform a NSCF calculation on a high-symmetry $\kk$-path to locate the position of the KS band edges
as these are the states we want to correct.

The input file is ...

We use [[getden_filepath]] to read the DEN.nc file instead of [[getden]] or [[irdden]] variables.

TODO: Mention ncdump and the input_string

In this particular case the CBM and the VMB ...

TODO: use and explain [[structure]]

<!--
!!! important

    Note that the EPH code can compute QP corrections only for $\kk$-points belonging to the $\kk$-mesh
    associated to WFK file.
    In some cases, it is not possible to generate a $\kk$-mesh containing the band edges.
    In this case, you will have to select the closest wavevectors in the grid.
-->

## Our first ZPR calculation

For our first ZPR we use a minimalistic input file that allows us to discuss
the basic input variables and the organization of the main output file.

Run the code, as usual, using:

```sh
abinit ...
```

Let's now discuss the meaning of the different variables in more details:

We use [[optdriver]] 7 and [[eph_task]] 4 to activate the computation of the full self-energy (real + imaginary part).
The paths to the external files are specified by [[getddb_filepath]], [[getwfk_filepath]], and [[getdvdb_filepath]].

[[ngkpt]] [[nshiftk]] and [[shitfk]]
[[ddb_ngqpt]] is set to XXX as this is the $\qq$-mesh we used in the DFPT part to generate the DDB and DVDB file.
The integations in $\qq$-space is done with the [[eph_ngqpt_fine]].
As [[eph_ngqpt_fine]] differs from [[ddb_ngqpt]], the code will automatically activate the interpolation of the DFPT potentials
as discussed in [introduction page for the EPH code](eph_intro).

The $\qq$-space integration is done [[eph_intmeth]]
[[zcut]]

Run the code with:

Let's have a look at the main output file...

!!! important

    To compute the imaginary part of $\Sigma^{\text{e-ph}}$ at the KS, we strongly recommend to use
    [[eph_task]] -4 as this option activates several important optimizations that are not possible
    when the full self-energy is wanted.
    Note, however, that [[eph_task]] -4 is not able to provide the full frequency dependence, only the
    values at the KS eigenvalue.
    The computation of spectral functions and Eliashberg functions therefore requires [[eph_task]] +4.

### Convergence wrt nband

At this point it should be not so difficult to write an input file to perform ZPR
calculations for different values of [[nband]] for fixed $\qq$-mesh and $\kk$ wavevectors.
It's just a matter of adding

```sh
ndtset
nband:
nband+

getwk1  0 # Cancel default
getwfk -1
```

Note how we use [[getwfk]] -1 to concatenate the different calculations to speedup a bit the calculation.
Obviously in a supercomputing center it's more advantageous to split the calculations and run them with independent jobs.
The price to pay is that we end up with a sequential input file.

An example of input file is given in

Run the calculation using

grep something gives

So XXX bands are needed to convergence the ZPR.

!!! important

    The multidatase syntax is quite handy when running small calculations but you can get a much better
    speedup if you split the calculation using different input files as these jobs are independent
    and can be executed in parallel.
    Note also that restart capabilities (see [[eph_restart]]) won't work in multidataset mode.

### How to reduce the number of bands with the Sternheimer method

In this section, we discuss how to take advantange of the Sternheimer equation to accelerate the convergence with [[nband]].
To activate the Sternheimer method, set [[eph_stern]] to 1
and use the [[getpot_filepath]] input variable to specify the external file storing the GS KS potential.
The Sternheimer equation is solved non-self-consistently using max [[nline]] NSCF iterations.
and the solver stops when the solution is converged within [[tolwfr]].
<!--
Default values are provided for these two variable are provided.
You might need to adjust the values for more complicated systems.
-->

!!! important

    The big advange of the Sternheimer method is that we don't need to compute WFK files with many
    empty bands to converge the self-energy.
    This means that one can use the computing power to densify the $\kk$-mesh while keeping the number of
    empty states at a reasonable level.
    Producing a WFK file with 1000 bands and a $100^3$ $\kk$-mesh is indeed way more expensive than
    performing the same computation with, let's say, 50 bands.

    Note however that the cost of the Sternheimer method quickly increases with [[nband]]
    due to the orthogonalization process. This means that a ZPR calculation with 300 bands
    **without** the Sternheimer method is much faster than the same computation done with [[eph_stern]] 1.
    As a matter of fact, one uses the Sternheimer method so that we don't need 300 bands to convergence the results.

### Convergence of the ZPR wrt the q-mesh

### How to compute the spectral function

### MPI parallelism and memory requirements

There are five different MPI levels that can be used to distribute the workload
and the most memory-demanding data structures.
By default, the code tries to reach some compromise between memory requirements and time to solution
by activating the parallelism over $\qq$-points if no other input is provided by the user.
You can however specify manually the MPI distribution across the five different levels
by using [[eph_np_pqbks]] (a list of 5 integers).
The product of these five numbers **must be equal** to the total number of MPI processes.
The first number gives the number of processes for the parallelization over perturbations.
The allowed value range between 1 and 3 × [[natom]], and should be a divisor
of 3 × [[natom]] to distribute the work equally.
The higher this number, the lower the memory requirements at the price of increased MPI communication.
The second number determines the parallelization over the $\qq$-points in the IBZ.
This parallelization level allows one to decrease both the computational time as well as memory although
it's not always possible to distribute the load equally among the processes.
The parallelization over bands is usually not relevant for mobility computations
as only a few states close to the VBM or CBM are considered.
It is however useful when the real part of the self-energy is needed.
The MPI parallelism over $\kk$-points and spins is very efficient
but it requires HDF5 with MPI-IO support and memory does not scale.
Use these additional levels if memory requirements are under control
and you want to boost the calculation.

There's an important difference with respect to [[eph_task]] -4 that is worth discussing in detail.
When computing the imaginary part at the KS energy for transport properties,
the EPH code is able to filter both $\kk$- and $\qq$-points so that only the relevant states around the band edge
are stored in memory.
Unfortunately, tn the case of full self-energy calculations, this filtering algorithm is not possible and each MPI process needs
to read and store all the KS wavefunctions in the IBZ so that one can compute e-ph matrix elements connecting $\kk$ to $\kq$.
In other words, the IBZ is not MPI-distributed and this leads to a significant increase in the memory requirements, especially
for dense meshes.
Fortunately, the code is able to MPI-distribute bands hence the memory required for this part will scale as [[nband]] / **np_band**.


If you are not using the Sternheimer method, it's important to use enough MPI proceeses for the band level
in order to decrease the memory requires to allocate the wavefunctions in $\GG$-space as

To recap:

1. Use enough *np_band* MPI processes to decrease the memory for the wavefunctions.
   Ideally, *np_band* should divide [[nband]] to distribute the work equally.
   Note that the band parallelism is beneficial also when the Sternheimer method is used
   as the Sternheimer solver will operate on MPI-distributed bands.
   Perhaps the parallel efficiency won't be perfect but memory for the wavefunctions will continue to scale.

2. Once the memory for the wavefunctions reaches a reasonable amount, activate the parallelism
   over perturbations in order to decrease the memory for $W(\rr, \RR, \text{3 natom})$.
   For better efficiency, *eph_nperts* should divide 3 * [[natom]].
   Usins *np_perts* = [[natom]] is usually a reasonable choice.

3. If the memory for the wavefuntions and $W(\rr, \RR, \text{3 natom})$ is under control,
   you may want to activate the $\qq$-point parallelism
   to speedup the calculation and decrease the memory required by the cache.

4. Finally, use the $\kk$-point parallelism if there are enough CPUs available to boost the calculation.
   Obviously, this kind of parallelism makes sense if you are computing QP corrections for sevaral $\kk$ points
   and the number of *np_kpt* MPI processes should be adjusted accordingly,
   Keep in mind, however, that each $\kk$-point has a different computational cost so load imbalance is expected.


Last but not least, remeber that setting [[boxcutmin]] to a value smaller than 2 will lead to a significant memory reduction.

### Estimating the ZPR with a generalized Fr\"ohlich model

Last but not least, one can estimate the correction to the zero-point renormalization
of the band gap in polar materials using a generalized Fr\"ohlich model based on
*ab initio** effective masses computed with DFPT [[cite:Laflamme2016]]
The formalism is detailed in XXX.
An example input file is available in [[test:v7_88]].
