---
authors: MG
---

# Zero-point renormalization of the band gap and temperature-dependent band gaps

This tutorial explains how to obtain the electron self-energy due to phonons, compute the zero-point
renormalization (ZPR) of the band gap as well as temperature-dependent band gaps and spectral functions.
We start with a very brief overview of the many-body formalism in the context of the electron-phonon (e-ph) interaction,
then we discuss how to evaluate the e-ph self-energy and perform typical convergence studies using MgO as example.

It is assumed the user has already completed the two tutorials [RF1](rf1) and [RF2](rf2),
and that he/she is familiar with the calculation of ground state and response properties,
in particular phonons, Born effective charges and dielectric tensor.
It goes without saying that one should have read the [introduction page for the EPH code](eph_intro)
before running these examples.

This lesson should take about 1.5 hour.

## Formalism

The electron-phonon self-energy, $\Sigma^{\text{e-ph}}$, describes the renormalization of the
charged quasi-particle excitation due to the interaction with phonons.
This term should be added to the electron-electron (e-e) self-energy $\Sigma^{\text{e-e}}$
that encodes many-body effects induced by the Coulomb interaction beyond the Hartree potential.
The e-e contribution can be estimated using, for instance, the $GW$ approximation but
in this tutorial we are mainly interested in $\Sigma^{\text{e-ph}}$ and its temperature dependence.

In semiconductors and insulators, indeed, most of temperature dependence of the electronic properties
at low T originates from the e-ph interaction
(and the thermal expansion of the unit cell that, however, is not treated in this lesson).
Corrections due to $\Sigma^{\text{e-e}}$ are obviously important as it is well known that KS gaps computed
at the LDA/GGA level are systematically underestimated with respect to experiments
but the temperature dependence of $\Sigma^{\text{e-e}}$
is rather small as long as $kT$ is much smaller than the fundamental gap.

In state-of-the-art *ab-initio* methods, the e-ph coupling is described
within DFT by expanding the KS effective potential in the nuclear displacements,
and the vibrational properties are obtained with DFPT [[cite:Gonze1997]], [[cite:Baroni2001]].
The e-ph self-energy consists of two terms: the **frequency-dependent Fan-Migdal** (FM) self-energy
and the **static and Hermitian Debye-Waller** (DW) part [[cite:Giustino2017]]

$$ \Sigma^\eph(\ww, T) = \Sigma^\FM(\ww, T) + \Sigma^{\DW}(T). $$

The diagonal matrix elements of the FM self-energy in the KS representation are given by

\begin{equation}
\begin{split}
    \Sigma^\FM_{n\kk}(\omega, T) =
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
with $T$ the temperature and $\ef$ the Fermi level that in turns depends on T and the number of electron per unit cell.
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
    Note also that [[zcut]] should be of the same order as the phonon frequency so around 0.01, 0.001 eV.

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
\Sigma_{n\kk}^{\DW}(T) = \sum_{\qq\nu m} (2 n_{\qq\nu}(T) + 1) \dfrac{g_{mn\nu}^{2,DW}(\kk, \qq)}{\ee_{n\kk} - \ee_{m\kk}},
\end{equation}

where $g_{mn\nu}^{2,\DW}(\kk,\qq)$ is an effective matrix element that, within the **rigid-ion approximation**,
can be expressed in terms of the $\gkq$ matrix elements using the invariance of the QP energies
under infinitesimal translation [[cite:Giustino2017]].

At the level of the implementation, the number of bands in the two sums is defined by [[nband]]
while the $\qq$-mesh is specified by [[eph_ngqpt_fine]] (or [[ddb_ngqpt]] if DFPT potentials should not be interpolated).
The list of temperatures in Kelvin is defined by the [[tmesh]] input variable.
For the sake of simplicity, the dependence on $T$ will be omitted in the following.
Keep in mind, however, that all the equations in which $\Sigma$ appears have an additional dependence on
the physical temperature T.

!!! important

    The EPH code takes advantage of (spatial and time-reversal) symmetries to reduce the BZ integration
    to an appropriate irreducible wedge, $\text{IBZ}_k$, defined by the little group of the $\kk$-point i.e.
    the set of point group operations of the crystal that leave the $\kk$-point invariant
    within a reciprocal lattice vector.
    Calculations at high-symmetry $\kk$-points such as $\Gamma$ are therefore much faster as there
    are more symmetries that can be exploited.

    This symmetrization procedure is activated by default and can be disabled by setting [[symsigma]]
    to 0 for testing purposes.
    Note that when [[symsigma]] is set to 1, the code performs a final average of the QP results
    within each degenerate subspace.
    As a consequence, **accidental degeneracies won't be removed when the sum is performed in the $\text{IBZ}_k$**.

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

### Quasi-particle corrections due to e-ph coupling

Strictly speaking, the quasi-particle (QP) excitations are defined by the solution(s)
in the complex plane of the equation $$ z = \ee_\nk + \Sigma_\nk^{\text{e-ph}}(z) $$
provided that the non-diagonal components of the self-energy can be neglected.
In practice, the problem is usually simplified by seeking approximated solutions along the real axis
following two different approaches: **on-the-mass-shell** and **linearized QP equation**.

In the on-the-mass-shell approximation, the QP energy is given by the real part
of the self-energy evaluated at the bare KS eigenvalue:

$$ \ee^\QP_\nk = \ee_\nk + \Re\, \Sigma_\nk^{\text{e-ph}}(\ee_\nk). $$

This approach is equivalent to a (thermal average of) standard time-dependent Rayleigh-Schrodinger perturbation theory.
In the linearized QP equation, on the contrary, the self-energy is Taylor-expanded around
the KS eigenvalue and the QP correction is obtained using

$$ \ee^\QP_\nk = \ee_\nk + Z_\nk\,\Re  \Sigma^\text{e-ph}_\nk(\ee_\nk) $$

with the renormalization factor $Z_\nk$ given by

$$
  Z_\nk= \left(1 - \Re\left[ \frac{\partial\Sigma^\text{e-ph}_{\nk}}{\partial\ee}\right]\Bigg|_{\ee=\ee_\nk} \right)^{-1}.
$$

Both approaches are implemented in ABINIT although it should be noted that, according to recent works,
the on-the-mass-shell approach provides results that are closer to those obtained
with more advanced techniques based on the cumulant expansion [[cite:Nery2018]].

!!! important

    The EPH code can compute QP corrections only for $\nk$ states that are present in the input WFK file
    (a similar requirement is present in the $GW$ code as well).
    As a consequence, the $\kk$-mesh ([[ngkpt]], [[nshiftk]], [[shiftk]]) for the WFK file
    should be chosen carefully especially if the band edge is not located at an high-symmetry $\kk$-point.

There are different approaches one can use to specify the set of $\nk$ states in $\Sigma_{\nk}$.
Each approach has pros and cons.

The most direct way consists in specifying explicitly the $\kk$-points and the band range
using the three variables: [[nkptgw]], [[kptgw]], [[bdgw]]
To compute the correction for the VBM/CBM at $\Gamma$ in silicon
(a non-magnetic semiconductor with 8 valence electrons per unit cell), one would use:

```sh
nkptgw 1
kptgw  0 0 0  # [3, nkptgw] array
bdgw   4 5    # [2, nkptgw] arary giving the initial and the last band index
              # for each nkptgw k-point
```

as the index of the valence band is given by 8 / 2 = 4.
Obviously this input file will only provide the ZPR of the optical gap as Si has an indirect bandgap.

!!! important

    When [[symsigma]] is set to 1 (default), the code may decide to enlarge the initial value of [[bdgw]]
    so that all degenerate states for that particular $\kk$-point are included in the calculation.

Alternatively, one can use [[gw_qprange]] or [[sigma_erange]]
Note that [[gw_qprange]] is mainly used to compute all the corrections for the occupied states plus
some conduction states while [[sigma_erange]] is usually employed for transport calculations with [[eph_task]] = -4.

### Spectral function and Eliashberg functions

The spectral function is defined by:

$$
A_\nk(\ww) = -\dfrac{1}{\pi} \dfrac{\Im \Sigma_\nk(\ww)} {(\ww - \ee_\nk - \Re \Sigma_\nk(\ww)) ^ 2 + \Im \Sigma_\nk(\ww) ^ 2}
$$

The computation of the spectral function requires the specification of [[nfreqsp]] and [[freqspmax]].
The results are stored in the SIGMAPH.nc file for each $\nk$ state specified by the user.

## Typical workflow for ZPR

A typical workflow for ZPR calculations requires the following steps:

1. **GS calculation** to obtain the WFK and the DEN file.
   The $\kk$-mesh should be dense enough to converge electronic and vibrational properties.
   Remember to set [[prtpot]] to 1 to produce the file with the KS potential required for the Sternheimer method.

2. **DFPT calculations** for all the IBZ $\qq$-points corresponding to the *ab-initio* [[ddb_ngqpt]] mesh
   used to perform the Fourier interpolation of the dynamical matrix and of the DFPT potentials.
   In the simplest case, the DFPT run uses the WFK produced in step #1
   provided the $\qq$-mesh is a submesh of the GS $\kk$-mesh.
   Remember to compute $\epsilon^{\infty}$, $Z^*$ and $Q^*$.

3. **NSCF computation** of a WFK file on a much denser $\kk$-mesh containing the wavevectors
   for which phonon-induced QP corrections are wanted. This job will use the DEN file produced in step #1.
   Remember to include enough empty states so that it is possible to perform
   convergence studies wrt [[nband]] afterwards.

4. **Merge the partial DDB and POT files** with *mrgddb* and *mrgdvdb*, respectively

5. **Start from the full DDB file, the DVDB file and the WFK file** obtained in step #3 to perform ZPR calculations
   with [[eph_task]] 4.

## Getting started

In this tutorial, we prefer to focus on the usage of the EPH code hence
we will be using **pre-computed DDB and DFPT POT files** to bypass the DFPT part.
We also provide a **precomputed DEN.nc file** that can be used to perform the NSCF calculations required
to generate the WKF file and the file with GS KS potential required for the Sternheimer equation.

If git is installed on your machine, one can easily fetch the entire repository with:

```sh
git clone ...
```

Alternatively, use *wget*:

```sh
wget
```

or *curl*:

```sh
curl
```

or simply copy the tarball by clicking the "download button" in the github interface.
Note **that the directory with the input files must be located in the same working directory**
in which you will be executing the tutorial.

The |AbiPy| script used to perform the prelimiary calculations is available here.
Note that AbiPy does not supports MultiDatasets so each directory corresponds to a single calculation.
In particular, all the DFPT task are in

To produce these files, we used the experimental parameters for hexagonal $MgB_2$ (a = 5.8317 and c/a= 1.1419)
and norm-conserving pseudopotentials with an energy cutoff [[ecut]] of 60 Ry.
All the calculations have been performed with a 40x40x40 [[ngkpt]] Gamma-centered grid for electrons,
and the Gaussian smearing [[occopt]] with [[tsmear]].
The DFPT computations have been done for a set of XXX irreducible $\qq$-points
corresponding to a $\Gamma$-centered 6x6x6 mesh.
This is the |AbiPy| script used to automate the GS + DFPT calculation:

### How to recover useful info from the output files

If what follows, we pretend we know nothing about MgO so that we can explain how to use
abitk and netcdf files to inspect the results of the previous calculations.
and use this piece of info to continue our calculations.

The input file of the GS run is stored in the DEN.nc file (*input_string* nc variable)
and one can easily access it with the *ncdump* utility

    ncdump -v input_string flow_zpr_mgo/w0/t0/outdata/out_DEN.nc

    FOOBAR

To print to terminal the crystalline structure, use:

```md
abitk crystal_print flow_zpr_mgo/w0/t0/outdata/out_DEN.nc

 ==== Info on the Cryst% object ====
 Real(R)+Recip(G) space primitive vectors, cartesian coordinates (Bohr,Bohr^-1):
 R(1)=  0.0000000  4.0182362  4.0182362  G(1)= -0.1244327  0.1244327  0.1244327
 R(2)=  4.0182362  0.0000000  4.0182362  G(2)=  0.1244327 -0.1244327  0.1244327
 R(3)=  4.0182362  4.0182362  0.0000000  G(3)=  0.1244327  0.1244327 -0.1244327
 Unit cell volume ucvol=  1.2975866E+02 bohr^3
 Angles (23,13,12)=  6.00000000E+01  6.00000000E+01  6.00000000E+01 degrees
 Time-reversal symmetry is present
 Reduced atomic positions [iatom, xred, symbol]:
    1)    0.0000000  0.0000000  0.0000000  Mg
    2)    0.5000000  0.5000000  0.5000000   O
```

<!--
To print the same info in a format that we can directly reuse in the ABINIT input file, use:
$ abitk crystal_print_abivars   flow_zpr_mgo/w0/t0/outdata/out_DEN.nc
-->

Since we want to compute the renormalization of the band gap due to phonons, it is also useful
to have a look at the gaps computed from the $\kk$-mesh used to generate the DEN file:

```md
abitk ebands_gaps flow_zpr_mgo/w0/t0/outdata/out_DEN.nc

  >>>> For spin  1
   Minimum direct gap =   4.45 (eV), located at k-point     : [ 0.0000E+00,  0.0000E+00,  0.0000E+00]
   Fundamental gap    =   4.45 (eV), Top of valence bands at: [ 0.0000E+00,  0.0000E+00,  0.0000E+00]
                                     Bottom of conduction at: [ 0.0000E+00,  0.0000E+00,  0.0000E+00]
   Valence Max:      4.51 (eV) at: [ 0.0000E+00,  0.0000E+00,  0.0000E+00]
   Conduction min:   8.96 (eV) at: [ 0.0000E+00,  0.0000E+00,  0.0000E+00]
```

!!! warning

    Our values for the gaps are consistent with the results for MgO given on the
    [materialsproject](https://materialsproject.org/materials/mp-1265/).
    Remember, however, that the values and the positions of the gaps may vary
    (in somes cases even significantly) depending on the kind of $\kk$-sampling employed.

    In this case, abitk reports the gaps computed from a $\kk$-mesh as the DEN file can only be produced
    by a SCF calculation that requires an BZ-mesh.
    **The results for MgO are OK simply because the CBM/VBM are at the $\Gamma$ point and this point
    belongs to our $\kk$-mesh**.
    Other systems (e.g. Si) may have the CBM/VBM at wavevectors that are not easily captured with a homogeneous mesh.
    **The most reliable approach to find the location of the CBM/VBM is to perform a band structure calculation
    on a high-symmetry $\kk$-path.**


To plot the band structure

```sh
abiopen.py flow_zpr_mgo/w0/t1/outdata/out_GSR.nc -e
```

If you don't remember the FFT mesh used to generate the DEN file, use

```sh
abitk hdr_print MgO_eph_zpr/flow_zpr_mgo/w0/t0/outdata/out_DEN.nc

===============================================================================
 ECHO of part of the ABINIT file header

 First record :
.codvsn,headform,fform = 9.1.5      80   52

 Second record :
 bantot,intxc,ixc,natom  =    96     0    11     2
 ngfft(1:3),nkpt         =    30    30    30     8
 nspden,nspinor          =     1     1
 nsppol,nsym,npsp,ntypat =     1    48     2     2
 occopt,pertcase,usepaw  =     1     0     0
 ecut,ecutdg,ecutsm      =  3.0000000000E+01  3.0000000000E+01  0.0000000000E+00
 ecut_eff                =  3.0000000000E+01
 qptn(1:3)               =  0.0000000000E+00  0.0000000000E+00  0.0000000000E+00
 rprimd(1:3,1)           =  0.0000000000E+00  4.0182361526E+00  4.0182361526E+00
 rprimd(1:3,2)           =  4.0182361526E+00  0.0000000000E+00  4.0182361526E+00
 rprimd(1:3,3)           =  4.0182361526E+00  4.0182361526E+00  0.0000000000E+00
 stmbias,tphysel,tsmear  =  0.0000000000E+00  0.0000000000E+00  1.0000000000E-02

 The header contain   4 additional records.
```

to print the header. Use `--prtvol 1` to output more records.

### Executing mrgddb

First of all, let's merge the partial DDB files with

```sh
mrgddb < teph4zpr_1.in
```

and the following input file:

{% dialog tests/tutorespfn/Input/teph4zpr_1.in %}

that lists the **relative paths** of the different partial DDB files and the syntax:

### Executing mrgdv

Now we can merge the DFPT potential with the *mrgdv* tool using the following input file:

{% dialog tests/tutorespfn/Input/teph4zpr_2.in %}

and the command

    mrgdvdb <  teph4zpr_2.in


!!! tip

    The number at the end of the POT file corresponds to the (*idir*, *ipert*) pertubation for that
    particular $\qq$-point.

    ```fortran
    percase = idir + ipert
    ```

    *idir* species the reduced direction ([1, 2, 3])
    whereas *ipert* specifies the kind of perturnation ([1, natom] if atomic displacement)
    All DFPT POT files with 1 <= index <= 3 x [[natom]] correspond to atomic pertubations.

In the output file produced by mrgdv 

{% dialog tests/tutorespfn/Refs/teph4zpr_2.stdout %}

there is a section for each $\qq$-point with the list of atomic perturbations 
that have been included in the database.

```md
 qpoint: [ 0.0000E+00,  0.0000E+00,  0.0000E+00] is present in the DVDB file
 The list of irreducible perturbations for this q vector is:
    1)  idir= 1, ipert=   1, type=independent, found=Yes
    2)  idir= 2, ipert=   1, type=symmetric, found=No
    3)  idir= 3, ipert=   1, type=symmetric, found=No
    4)  idir= 1, ipert=   2, type=independent, found=Yes
    5)  idir= 2, ipert=   2, type=symmetric, found=No
    6)  idir= 3, ipert=   2, type=symmetric, found=No
```

**symmetric** means that that particular *(idir, ipert)* can be reconstructed by symmetry from the **independent** entries.
All the independent entries should be present and you should get the following message:
at the end of the output file, 

```md
 All the independent perturbations are available
 Done
```

this indicates that our DVDB database is complete in the sense that the EPH code will be able
to reconstruct by symmetry all the 3 [[natom]] perturbations for each $\qq$-qpoint.

!!! warning

    If you don't get this message, it means that your DVDB cannot be used by the EPH code.
    In this case, check carefully your input files and the list of files that have been merged.

## Computing the WFK files

At this point we have all the ingredients required for the phonons and the DFPT potential and we
can finally start to generate the WFK files.

For our first NSCF calculation, we use a xxx $\kk$-mesh and XXX bands.
The number of bands is sufficiently large so that we can perform initial convergence studies.
We also perform a NSCF calculation on a high-symmetry $\kk$-path to locate
the position of the KS band edges as these are the states we want to correct.

We use [[getden_filepath]] to read the DEN.nc file instead of [[getden]] or [[irdden]].

Note that in all the input files of the tutorial, we will be using the [[structure]] variable
to initialize the unit cell from an external input file so that we don't need
to repeat it over and over again.

```sh
 structure = "abifile:MgO_eph_zpr/flow_zpr_mgo/w0/t0/outdata/out_DEN.nc"
```

```
 natom 2
 acell    1.0    1.0    1.0
 rprim
    0.0000000000    4.0182361526    4.0182361526
    4.0182361526    0.0000000000    4.0182361526
    4.0182361526    4.0182361526    0.0000000000
 xred_symbols
    0.0000000000    0.0000000000    0.0000000000 Mg
    0.5000000000    0.5000000000    0.5000000000 O
```

The list of pseudopotential given in the [[pseudos]] variable must be consistent
so the pseudo from Mg comes first as Mg is the first site in `xred_symbols`.

The input file is ...

{% dialog tests/tutorespfn/Input/teph4zpr_3.in %}


```sh
 nband 100
 nbdbuf 10
 ecut 35.0

 ngkpt 4 4 4
 nshiftk 1
 shiftk 0 0 0

 nstep 150
 tolwfr 1e-15
 iscf -2          # NSCF run from DEN.nc

 ngfft 30 30 30   # Important: Must use same FFT mesh as the one in DEN.nc

 getden_filepath "MgO_eph_zpr/flow_zpr_mgo/w0/t0/outdata/out_DEN.nc"
```

Note the use of [[nbdbuf]].

!!! important

    For mobility calculations, we have seen that it is possible
    to reduce significantly the cost of the WFK computation
    by restricting the NSCF calculation to the $\kk$-points inside the electron (hole) pockets.
    Unfortunately, this trick is not possible when computing the full self-energy.
    On the other hand, ZPR calculations can take advange of the Sternheimer method to reduce the number of [[nband]].

## Our first ZPR calculation

For our first ZPR calculation, we use a very minimalistic input file that allows us
to discuss the most important input variables and the info reported in the main output file.

First of all, run the code using:

```sh
abinit teph4zpr_4.in > teph4zpr_4.log 2> err &
```

with the following input file:

{% dialog tests/tutorespfn/Input/teph4zpr_4.in %}

!!! tip

    You may want to run the examples in parallel with MPI using e.g.

        mpirun -n 2 abinit teph4zpr_4.in > teph4zpr_4.log 2> err &

    to speed up the calculations of the tutorials.
    The EPH code will automatically distribute the workload using a predefined distribution scheme
    (not necessarily the most efficient one in terms of parallel efficiency though).
    In the last part of the tutorial, we explain how to specify a particular
    MPI distribution scheme with [[eph_np_pqbks]].

Let's now discuss the meaning of the different variables in more details.

We use [[optdriver]] 7 to enter the EPH code while [[eph_task]] 4 activates 
the computation of the full self-energy (real + imaginary part).
The paths to the external files (DDB, WFK, DVDB) are specified 
by the three variables [[getddb_filepath]], [[getwfk_filepath]], and [[getdvdb_filepath]]:

```sh
getddb_filepath "teph4zpr_1_DDB"
ddb_ngq pt 4 4 4          # The code expects to find in the DDB
                          # all the IBZ q-points corresponding to a 4x4x4 q-mesh

getdvdb_filepath "teph4zpr_2_DVDB"
getwfk_filepath "teph4zpr_3o_WFK"  # 4x4x4 k-mesh with 70 bands
```

The mesh for electrons ([[ngkpt]], [[nshiftk]] and [[shiftk]]) is the one used to generate the input WFK file.
[[ddb_ngqpt]] is set to 4x4x4 as this is the $\qq$-mesh we used in the DFPT part to generate the DDB and DVDB file.
but the integration in $\qq$-space is performed with the [[eph_ngqpt_fine]] mesh.
As [[eph_ngqpt_fine]] differs from [[ddb_ngqpt]], the code will automatically activate the interpolation of the DFPT potentials
as discussed in [introduction page for the EPH code](eph_intro).
The $\qq$-space integration is defined by [[eph_intmeth]] [[zcut]]

We can now have a look at the main output file.

{% dialog tests/tutorespfn/Refs/teph4zpr_4.out %}

First of all, we have a section that summarizes the most important parameters:

```md
 Number of bands in e-ph self-energy sum: 30
 From bsum_start: 1 to bsum_stop: 30
 Symsigma: 1 Timrev: 1
 Imaginary shift in the denominator (zcut): 0.010 [eV]
 Method for q-space integration:  Standard quadrature
 Both Real and Imaginary part of Sigma will be computed.
 Number of frequencies along the real axis: 0 , Step: 0.000 [eV]
 Number of frequency in generalized Eliashberg functions: 0
 Number of temperatures: 1 From: 0.000000E+00 to 0.000000E+00 [K]
 Ab-initio q-mesh from DDB file: [4, 4, 4]
 Q-mesh used for self-energy integration [ngqpt]: [8, 8, 8]
 Number of q-points in the IBZ: 29
 asr: 1 chneut: 1
 dipdip: 1 symdynmat: 1
 No special treatment of Frohlich divergence in gkq for q --> 0
 Number of k-points for self-energy corrections: 1
 List of K-points for self-energy corrections:
   1     1  [ 0.0000E+00,  0.0000E+00,  0.0000E+00]   6    9
```

Note how [[asr]] and [[chneut]] and [[dipdip]] are automaticall activated

Then we find another section defining how to different dimensions are distributed across the MPI processes:

```md
 === MPI parallelism ===
P Allocating and summing bands from my_bsum_start: 1 up to my_bsum_stop: 30
P Number of CPUs for parallelism over perturbations: 1
P Number of perturbations treated by this CPU: 6
P Number of CPUs for parallelism over q-points: 1
P Number of q-points in the IBZ treated by this proc: 29 of 29
P Number of CPUs for parallelism over bands: 1
P Number of CPUs for parallelism over spins: 1
P Number of CPUs for parallelism over k-points: 1
P Number of k-point in Sigma_nk treated by this proc: 1 of 1

 Cannot find eph_ngqpt_fine q-points in DVDB --> Activating Fourier interpolation.
```

In this case we are running in sequential but things will start to change here if you start to
use more than one core. [[eph_np_pqbks]]

Finally we have the section with the QP results:

```md
================================================================================
 Final results in eV.
 Notations:
     eKS: Kohn-Sham energy. eQP: quasi-particle energy.
     eQP - eKS: Difference between the QP and the KS energy.
     SE1(eKS): Real part of the self-energy computed at the KS energy, SE2 for imaginary part.
     Z(eKS): Renormalization factor.
     FAN: Real part of the Fan term at eKS. DW: Debye-Waller term.
     DeKS: KS energy difference between this band and band-1, DeQP same meaning but for eQP.
     OTMS: On-the-mass-shell approximation with eQP ~= eKS + Sigma(omega=eKS)
     TAU(eKS): Lifetime in femtoseconds computed at the KS energy.
     mu_e: Fermi level for given (T, nelect)


K-point: [ 0.0000E+00,  0.0000E+00,  0.0000E+00], T:    0.0 [K], mu_e:    7.568
   B    eKS     eQP    eQP-eKS   SE1(eKS)  SE2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP
   6   4.474    4.572    0.098    0.152   -0.006    0.644   -0.265    0.417    0.000    0.000
   7   4.474    4.572    0.098    0.152   -0.006    0.644   -0.265    0.417    0.000    0.000
   8   4.474    4.572    0.098    0.152   -0.006    0.644   -0.265    0.417    0.000    0.000
   9   8.959    8.857   -0.102   -0.106   -0.000    0.963    0.107   -0.213    4.485    4.285

 KS gap:    4.485 (assuming bval:8 ==> bcond:9)
 QP gap:    4.285 (OTMS:    4.227)
 QP_gap - KS_gap:   -0.200 (OTMS:   -0.258)
```

!!! important

    To compute the imaginary part of $\Sigma^{\text{e-ph}}$ at the KS energy, we **strongly recommend** to use
    [[eph_task]] -4 as this option activates several important optimizations that are not possible
    when the full self-energy is wanted.
    Note, however, that [[eph_task]] -4 is not able to provide the full frequency dependence, only the
    value of the imaginary part at the KS eigenvalue.
    The computation of spectral functions and Eliashberg functions therefore requires [[eph_task]] +4.


Since we have computed the phonon band structure with ([[ph_ndivsm]], [[ph_nqpath]], [[ph_qpath]]
and the phonon DOS with [[ph_ngqpt]]

```sh
abiopen.py teph4zpr_4o_PHBST.nc -e
```

### Exercise:

 Change the input file to use [[mixprec]] 1 and [[boxcutmin]] 1.1.
 Rerun the calculation and compare with the previous results.
 Do you see significant differences? What about the wall-time?

## Convergence wrt nband

At this point it should be not so difficult to write an input file to perform ZPR
calculations for different values of [[nband]] for fixed $\qq$-mesh and $\kk$ wavevectors.
It's just a matter of adding:

```sh
 ndtset 3
 nband: 20
 nband+ 20
```

An example of input file is given in

{% dialog tests/tutorespfn/Input/teph4zpr_5.in %}

Run the calculation using

grep something gives

{% dialog tests/tutorespfn/Refs/teph4zpr_5.out %}

```sh
abicomp.py sigeph teph4zpr_5o_DS*_SIGEPH.nc -e
```

The results are quite disappointing in the sense that the QP corrections do not converge at all!
In part this is due to the coarse nband 300 should be enough if you use a much denser $\kk$-mesh.
Note however that this typical of many-body approaches based on sum over states in the sense 
that QP corrections convergence slowly with respect to the number of bands.

!!! important

    The multidataset syntax is quite handy when running small calculations but you can get a much better
    speedup if you split the calculation using different input files as these jobs are independent
    and can be executed in parallel.
    Note also that restart capabilities (see [[eph_restart]]) won't work in multidataset mode.

## How to reduce the number of bands with the Sternheimer method

In this section, we discuss how to use the Sternheimer method to accelerate the convergence with [[nband]].

To activate the Sternheimer approach, we need to set [[eph_stern]] to 1
and use the [[getpot_filepath]] input variable to specify the external file with the GS KS potential
This file is produced at the end of the GS calculations provided we set [[prtpot]] to 1 in the input file.

The Sternheimer equation is solved non-self-consistently using max [[nline]] NSCF iterations
and the solver stops when the first-order wavefunction is converged within [[tolwfr]].
Default values for these two variables are provided if they are not specified by the user in the input file.
The code will abort if the algorithm cannot converge the solution.

In brief, we need to add the following section to our initial EPH input:

```sh
eph_stern 1
getpot_filepath "MgO_eph_zpr/flow_zpr_mgo/w0/t0/outdata/out_POT.nc"

# nline 100 tolwfr 1e-16 # Default values
```

An example of input file is provided in:

{% dialog tests/tutorespfn/Input/teph4zpr_6.in %}

Run the calculation, as usual, with:

    abinit teph4zpr_6.in > teph4zpr_6.log 2> err &

{% dialog tests/tutorespfn/Refs/teph4zpr_6.out %}

Let's extract the results:

```
grep OTMS teph4zpr_6.out
```

```
abicomp.py sigeph teph4zpr_6o_DS*_SIGEPH.nc -e
```

Now the convergence is much better and the ZPR valu is converged with 1 meV for nband ??
This is the value one should obtain when summing all the bands up to [[mpw]]

Exercise:

    Generate a new WFK with nband ~ mpw and run a ZPR calculation without the Sternheimer method.
    Compare the sum-over-states results with Sternheimer.
    Why it is not a good idea to use nband > mpw?


!!! important

    The big advange of the Sternheimer method is that we don't need to compute WFK files with many
    empty bands to converge the self-energy.
    This means that one can use the computing power to densify the $\kk$-mesh while keeping the number of
    empty states at a reasonable level.
    Producing a WFK file with 1000 bands and a $100^3$ $\kk$-mesh is indeed way more expensive than
    performing the same computation with, let's say, 50 bands.

    Note however that the cost of the Sternheimer method quickly increases with [[nband]]
    due to the orthogonalization process. This means that a ZPR calculation with 300 bands (occ + empty)
    **without** the Sternheimer method is much faster than the same computation done with [[eph_stern]] 1.
    As a matter of fact, one uses the Sternheimer method so that we don't need 300 bands to convergence the results.

## Convergence of the ZPR wrt the q-mesh

In the previous sections, we found that XXX bands with the Sternheimer method are enough to converge.
In this part of the tutorial, we perform a convergence study wrt to the q-sampling

In the third input *teph4zpr_7.in*, we have computed WFK files on different $\kk-$-meshes 
and a relatively small number of empty states (XXX).
We can finally use these extra WFK files to perform ZPR calculations

```sh
 ndtset 3

 ngkpt1 4 4 4;    eph_ngqpt_fine1 4 4 4;    getwfk_filepath1 "teph4zpr_3o_DS1_WFK"
 ngkpt2 8 8 8;    eph_ngqpt_fine2 8 8 8;    getwfk_filepath2 "teph4zpr_3o_DS2_WFK"
 ngkpt3 12 12 12; eph_ngqpt_fine3 12 12 12; getwfk_filepath3 "teph4zpr_3o_DS3_WFK"
```

{% dialog tests/tutorespfn/Input/teph4zpr_7.in %}

Run the calculation, as usual, with:

    abinit teph4zpr_7.in > teph4zpr_7.log 2> err &

Extracting the ZPR from the output file 

{% dialog tests/tutorespfn/Refs/teph4zpr_7.out %}

as a function of the number of k-points in the IBZ gives:

Obviously we are still far from convergence: one should test denser $\kk$-meshes and, last but not least,
monitor the behaviour for smaller [[zcut]].
Unfortunately, these additional convergence tests cannot be covered in this tutorial
and they are left as extra excercise.

<!--
### How to compute the spectral function

The maximum [[nfreqsp]] and [[freqspmax]].
-->

### MPI parallelism and memory requirements

<!--
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
-->

There is an important difference with respect to [[eph_task]] -4 that is worth discussing in detail.
When computing the imaginary part at the KS energy for transport properties,
the EPH code is able to filter both $\kk$- and $\qq$-points so that only the relevant states 
around the band edge are stored in memory.
Unfortunately, in the case of full self-energy calculations, this filtering algorithm 
is not possible and each MPI process needs to read and store all the KS wavefunctions in the IBZ 
so that one can compute e-ph matrix elements connecting $\kk$ to $\kq$.
In other words, the wavefunctions in the IBZ are not MPI-distributed and this leads to a significant 
increase in the memory requirements, especially for dense meshes.
Fortunately, the code is able to MPI-distribute bands hence the memory required 
for the orbitals will scale as [[nband]] / **np_band**.

<!--
If you are not using the Sternheimer method, it's important to use enough MPI processes for the band level
in order to decrease the memory required to allocate the wavefunctions in $\GG$-space as
-->

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
   Obviously, this kind of parallelism makes sense if you are computing QP corrections for several $\kk$ points
   and the number of *np_kpt* MPI processes should be adjusted accordingly,
   Keep in mind, however, that each $\kk$-point has a different computational cost so load imbalance is expected.

Finally, remember that setting [[boxcutmin]] to a value smaller than 2 (e.g. 1.1)
will lead to a significant memory reduction.

### Estimating the ZPR with a generalized Fr\"ohlich model

Last but not least, one can estimate the correction to the ZPR in polar materials
using a generalized Fr\"ohlich model based on *ab initio** effective masses
computed with DFPT [[cite:Laflamme2016]]
The formalism is detailed in XXX.
An example input file is available in [[test:v7_88]].
