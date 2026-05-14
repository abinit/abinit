---
authors: MG
---

# Zero-point renormalization of the band gap and temperature-dependent band gaps with GWPT

This tutorial is similar to [eph4zpr](../tutorial/eph4zpr.md).
Also in this lesson, we compute the electron self-energy due to phonons, obtain the zero-point
renormalization (ZPR) of the band gap and temperature-dependent band energies within the harmonic approximation.
The main difference with respect to [eph4zpr](../tutorial/eph4zpr.md), is that, in this lesson,
the e-ph matrix elements are computed within the GWPT formalism [[cite:Li2019]].

It is assumed that the user has already completed the two tutorials [RF1](../tutorial/rf1.md) and [RF2](../tutorial/rf2.md),
and that they are familiar with the calculation of ground state (GS) and response properties,
in particular phonons, Born effective charges (BECs) and the high-frequency dielectric tensor.

The user should have read the [introduction tutorial for the EPH code](../tutorial/eph_intro.md),
and (very important) the description of the [gstore-based approach](../tutorial/gstore.md), before running these examples.

Also, you are kindly invited to read the [first GW tutorial](../tutorial/gw1.md) if you are not familiar
with the GW implementation in Abinit.
A brief description of the formalism and of the equations implemented in the
code can be found in the [[theory:mbt|GW_notes]].

Note that in this tutorial, we will be using the conventional GW implementation in Fourier space and frequency-domain
with quartic scaling with the number of atoms.
A more efficient algorithm with cubic scaling is provided by the GWR code described in the [GWR1 tutorial](../tutorial/gwr1.md).
Unfortunately, at the time of writing, the GWR code is not yet interfaced with GWPT.

This lesson should take about 2.0 hours.

## Formalism

In the GWPT method [[cite:Li2019]], the e-ph matrix elements are computed by replacing the
first-order change of the KS Hamiltonian induced by a phonon with the variation of the $GW$ self-energy:

$$
g_{m n \nu \mathbf{k}\mathbf{q}}^{GW}(\varepsilon) = g_{m n \nu \mathbf{k} \mathbf{q}}^{\mathrm{KS}} +
\langle\psi_{m \mathbf{k}+\mathbf{q}} | \Delta_{\nu\mathbf{q}} \Sigma^{\mathrm{el}}(\varepsilon) | \psi_{n\mathbf{k}}\rangle
 - \langle \psi_{m \mathbf{k} + \mathbf{q}} | \Delta_{\nu\mathbf{q}} V^{\mathrm{xc}}[\rho^\mathrm{v}] | \psi_{n\mathbf{k}} \rangle
$$

with

\begin{align}\label{eq:sigma_g}
 \langle   \psi_{m \mathbf{k} + \mathbf{q}}  | \partial_{\kappa\alpha\mathbf{q}} \Sigma^{\mathrm{el}}(\varepsilon ) | \psi_{n\mathbf{k}}\rangle =&
  \frac{i}{2 \pi}\sum_{n'\mathbf{G G}^{\prime}} \int_\mathrm{BZ} \frac{\mathrm{d}\mathbf{p}}{\Omega}\langle\psi_{m \mathbf{k}+\mathbf{q}} | e^{i(\mathbf{p}+\mathbf{G}) \cdot \mathbf{r}} | \partial_{\kappa\alpha\mathbf{q}} \psi_{n^{\prime} \mathbf{k}-\mathbf{p}} \rangle \langle\psi_{n^{\prime} \mathbf{k}-\mathbf{p}} |e^{-i\left(\mathbf{p}+\mathbf{G}^{\prime}\right) \cdot \mathbf{r}^{\prime}}| \psi_{n \mathbf{k}}\rangle \widetilde{W}_{n'\mathbf{k-p,GG'}}(\varepsilon) \nonumber \\
   &+ \langle\psi_{m \mathbf{k}+\mathbf{q}} | e^{i(\mathbf{p}+\mathbf{G}) \cdot \mathbf{r}} | \psi_{n^{\prime} \mathbf{k}+\mathbf{q}-\mathbf{p}}
\rangle \langle\partial_{\kappa\alpha-\mathbf{q}} \psi_{n^{\prime} \mathbf{k}+\mathbf{q}-\mathbf{p}} | e^{-i (\mathbf{p}+\mathbf{G}^{\prime} ) \cdot \mathbf{r}^{\prime}} | \psi_{n \mathbf{k}}\rangle  \widetilde{W}_{n'\mathbf{k+q-p,GG'}}(\varepsilon).
\end{align}


and $\widetilde{W}$ is the frequency convolution of the screened Coulomb interaction:

\begin{align}
    \widetilde{W}_{n\mathbf{k-p,GG'}}(\varepsilon) \equiv & \int  \frac{\mathrm{d}\varepsilon' \, W_{\mathbf{p,GG'}}(\varepsilon')e^{i\eta\varepsilon'}}{\varepsilon-\varepsilon_{n\mathbf{k-p}}+ \varepsilon' +i\eta \sign(\varepsilon_{n\mathbf{k-p}}-\varepsilon^{\mathrm{F}})} \\
W_{\mathbf{p,GG'}}(\varepsilon') =& \frac{1}{V} \int_V \mathrm{d} \mathbf{r} \mathrm{d} \mathbf{r}^{\prime} \, W(\mathbf{r}, \mathbf{r}'; \varepsilon')\nonumber \\
    & \times  e^{-i(\mathbf{p}+\mathbf{G}) \cdot \mathbf{r}}  e^{i(\mathbf{p}+\mathbf{G}^{\prime}) \cdot \mathbf{r}^{\prime}},
    \label{eq:screeningGG}
\end{align}


Besides the summations over empty states $n'$, the equation requires the knowledge of the screened interaction $W$
as well as the (full) first-order derivative of the KS states
$\partial_{\kappa\alpha\mathbf{q}} \psi_{n^{\prime}}$
induced by an atomic displacement of atom $\kappa$ along direction $\alpha$ modulated by the wavevector $\mathbf{q}$.
Both terms can be computed by Abinit.

The screened interaction $W$ is computed in terms of a sum of states with [[optdriver]] 3, and the results
are stored in the SCR file.
The first-order derivative of the KS states, on the contrary, is computed on the fly by the GWPT subdriver
by solving a non-self-consistent (NSCF) Sternheimer equation as explained in the sections below.

Please note that ZPR computations at the GWPT level are still a field of active research,
especially in polar materials where additional long-range (LR) terms of many-body nature appear in the e-ph matrix elements.
In this tutorial, we won't be able to converge the calculation; hence, we mainly focus on explaining the different
steps involved and the input parameters affecting the quality of the calculation and the predictive power.
Hopefully, additional techniques and algorithmic improvements will be made available in forthcoming Abinit versions
in order to accelerate GWPT calculations without spoiling accuracy (stay tuned).

## Getting started

[TUTORIAL_README]

Before beginning, you might consider working in a different subdirectory as for the other tutorials.
Why not create Work_eph4zpr_gwpt in $ABI_TESTS/tutorespfn/Input?

```sh
cd $ABI_TESTS/tutorespfn/Input
mkdir Work_eph4zpr_gwpt
cd Work_eph4zpr_gwpt
```

In this tutorial, we prefer to focus on the use of the GWPT subdriver of the EPH code hence
we will be using **pre-computed** DDB and DFPT POT and DEN files to bypass the DFPT part.
We also provide a GS DEN.nc file to initialize the NSCF calculations
and a GS POT file with the KS potential required to solve the NSCF Sternheimer equation.

If *git* is installed on your machine, one can easily fetch the entire repository with:

```sh
git clone https://github.com/abinit/MgO_eph_zpr.git
```

Alternatively, use *wget*:

```sh
wget https://github.com/abinit/MgO_eph_zpr/archive/master.zip
```

or *curl*:

```sh
curl -L https://github.com/abinit/MgO_eph_zpr/archive/master.zip -o master.zip
```

or simply copy the tarball by clicking the "download button" available in the github web page,
unzip the file and rename the directory with:

```sh
unzip master.zip
mv MgO_eph_zpr-master MgO_eph_zpr_gwpt
```

!!! warning

    The directory with the precomputed files must be located in the same working directory
    in which you will be executing the tutorial and must be named `MgO_eph_zpr_gwpt`.


## Merging partial DDB, DFPT POT and DEN files

The first steps are similar to the ones used in the [eph4zpr tutorial](../tutorial/eph4zpr.md).
First, we merge the partial DDB files containing the dynamical matrix
evaluated for different $\mathbf{q}$-points and atomic perturbations by executing:

```sh
mrgddb < teph4zpr_gwpt_1.abi
```

with the following input file:

{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_1.abi %}

that lists the **relative paths** of the **partial DDB files** in the `MgO_eph_zpr` directory.

Then we merge the DFPT potential with the *mrgdv* tool using the command.

```sh
mrgdv < teph4zpr_gwpt_2.abi
```

with the following input file:

{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_2.abi %}

that lists the relative paths of the partial DFPT POT files in the MgO_eph_zpr directory.

Note that, in the case of GWPT, we also need to merge the files with the first-order derivatives of the density,
as the GWPT code needs to compute:

$$
\langle\psi_{m \mathbf{k}+\mathbf{q}} | \Delta_{\nu\mathbf{q}} \Sigma^{\mathrm{el}}(\varepsilon) | \psi_{n\mathbf{k}}\rangle
 - \langle \psi_{m \mathbf{k} + \mathbf{q}} | \Delta_{\nu\mathbf{q}} V^{\mathrm{xc}}[\rho^\mathrm{v}] | \psi_{n\mathbf{k}} \rangle
$$

This is done by executing

```sh
mrgdv < teph4zpr_gwpt_3.abi
```

with the following input file:

{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_3.abi %}


!!! Important

    The first-order densities are not written by default when performing DFPT calculations.
    Remember to use [[prtden]] 1 during the DFPT calculation.


## Computing the WFK files with empty states

At this point, we need to generate a WFK file with empty bands by performing an NSCF KS calculation
starting from a well-converged GS  density.
This WFK file will then be used to compute $W$, the $G_0W_0$ self-energy, and the GWPT matrix elements.

Let us start the NSCF calculation immediately by issuing:

```sh
mpirun -n 4 abinit teph4zpr_gwpt_4.abi > teph4zpr_4.log 2> err &
```

with the input file given by:

{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_4.abi %}

Here, we use a 4x4x4 $\Gamma$-centered $\kk$-mesh and 110 bands.
Note the use of [[getden_filepath]] to read the DEN.nc file instead of [[getden]] or [[irdden]].

!!! important

    This calculation is relatively fast.
    If you want to use more processors, remember to set [[autoparal]] to 1 in the input file.
    Note, however, that [[autoparal]] is only available in the GS part.

At this point, it is worth commenting about the use of [[nbdbuf]].
As mentioned in the documentation, **the highest energy states require more iterations to converge**.
To avoid wasting precious computing time, we use a buffer that is ~10% of [[nband]].
This trick significantly reduces the wall-time as the NSCF calculation completes
only when the first [[nband]] - [[nbdbuf]] states are converged within [[tolwfr]].
Obviously, one should not use the last [[nbdbuf]] states in the subsequent EPH calculation.

Another point worth mentioning here is that the $\kk$-mesh used in this step **strictly defines** the
electronic wavevectors $\kk$ that can be probed when computing the corrections to the QP energies,
both in the $GW$ part and in the $GWPT$ code.
In other words, if your CBM (VBM) is located at $\kk_0$, the $\kk$-mesh **must contain** $\kk_0$, as the present
implementation is not yet able to exploit interpolation techniques such as Wannierization to densify the $\kk$-mesh.

Before embarking on expensive GWPT, please make sure to understand very well the position of the KS band edges
by computing the KS band structure on a high-symmetry $\kk$-path.
In the case of MgO, the scenario is optimal as we are dealing with a direct-gap semiconductor with $\kk_0 = \Gamma$.
ZPR computations for the **fundamental band gap** in indirect band-gap semiconductors such as diamond
are significantly more challenging, as it is not always possible to find a $\kk$ mesh that can capture
the position of both the CBM and the VBM at the same time.

[HDIAGO_README]


## Computing the screened interaction W and QP corrections with one-shot GW

In this section, we use the WFK file generated in the previous section to compute the RPA polarizability,
the screened interaction $W$ and, finally, the QP corrections with one-shot GW.

To compute the screening, execute e.g.:

```sh
mpirun -n 4 abinit teph4zpr_gwpt_5.abi > teph4zpr_5.log 2> err &
```

with the input file given by:

{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_5.abi %}

!!! tip

    Since GW is more expensive than KS, you may want to use more MPI processes to speedup the calculation.
    Note that there is no need to specify how to distribute the work among processes as
    both the screening driver and the sigma driver are parallelized over bands.
    Clearly, one cannot use more MPI processes than [[nband]].

    An additional level of parallelism over $\qq$-points is available in the screening part.
    One can indeed split the computation of $W$ over $\qq$-points using [[nqptdm]] and [[qptdm]]
    and then merge the partial results with the `mrgscr` tool.


Here we generate a SCR file with only two frequencies: the static limit $\omega=0$
and $\omega=i\omega_p$, with $\omega_p$ being the Drude plasmon frequency computed from the number of (valence) electrons.
This is the **default behavior** of the screening driver (see also [[[ppmfrq]]).
The SCR file will then be used to construct the plasmon-pole model (PPM) when computing the self-energy
with [[optdriver]] 4 in the next step (see alsp [[ppmodel]])

To compute the $G_0W_0$ self-energy, execute e.g.:

```sh
mpirun -n 4 abinit teph4zpr_gwpt_6.abi > teph4zpr_6.log 2> err &
```

with the input file given by:

{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_6.abi %}

Let's now have a look at the QP results reported in the output file:

{% dialog tests/tutorespfn/Refs/teph4zpr_gwpt_6.abo %}

The KS bands gaps are

```
 === KS Band Gaps ===
  >>>> For spin  1
   Minimum direct gap =   4.4790 [eV], located at k-point      :   0.0000  0.0000  0.0000
   Fundamental gap    =   4.4790 [eV], Top of valence bands at :   0.0000  0.0000  0.0000
                                       Bottom of conduction at :   0.0000  0.0000  0.0000
```

and the QP corrections are summarized in this section:

```
 matrix elements of self-energy operator (all in [eV])

 Perturbative Calculation

--- !SelfEnergy_ee
iteration_state: {dtset: 1, }
kpoint     : [   0.000,    0.000,    0.000, ]
spin       : 1
KS_gap     :    4.479
QP_gap     :    6.528
Delta_QP_KS:    2.049
data: !SigmaeeData |
     Band     E0 <VxcDFT>   SigX SigC(E0)      Z dSigC/dE  Sig(E)    E-E0       E
        1 -68.031 -34.373 -58.114  13.281   0.273  -2.665 -37.227  -2.854 -70.885
        2 -34.759 -32.313 -45.678   6.515   0.808  -0.238 -37.845  -5.533 -40.292
        3 -34.759 -32.313 -45.678   6.515   0.808  -0.238 -37.845  -5.533 -40.292
        4 -34.759 -32.313 -45.678   6.515   0.808  -0.238 -37.845  -5.533 -40.292
        5 -12.564 -18.528 -30.572   9.400   0.667  -0.500 -20.290  -1.762 -14.326
        6   4.490 -20.036 -25.441   4.181   0.803  -0.246 -21.019  -0.983   3.507
        7   4.490 -20.036 -25.441   4.181   0.803  -0.246 -21.019  -0.983   3.507
        8   4.490 -20.036 -25.441   4.181   0.803  -0.246 -21.019  -0.983   3.507
        9   8.969 -13.280  -8.623  -3.404   0.851  -0.175 -12.214   1.066  10.035
       10  20.266  -9.448  -3.582  -4.476   0.861  -0.162  -8.251   1.197  21.462
       11  20.266  -9.448  -3.582  -4.476   0.861  -0.162  -8.251   1.197  21.462
       12  20.266  -9.448  -3.582  -4.476   0.861  -0.162  -8.251   1.197  21.462
...
```

For the meaning of the different columns, please consult the [GWR1 tutorial](../tutorial/gwr1.md).
The experimental gap of MgO is 7.67 eV.
A well converged G0W0 calculation should give 7.25 eV while our calculation gives 6.528.

At this stage, one should perform convergence studies for [[nband]], [[ecuteps]], [[ecutsigx]],
and [[ngkpt]] to ensure that the GW results are reasonably well converged.
The converged parameters can then be reused in the subsequent GWPT calculation.

For the sake of conciseness and performance reasons, these convergence studies are omitted here.

## Computing e-ph matrix elements with GWPT

At this point, we can finally start our first GWPT calculation by issuing:

```sh
mpirun -n 4 abinit teph4zpr_7.abi > teph4zpr_7.log 2> err &
```

!!! tip

    Feel free to use more MPI processes here as GWPT are expensive.
    The code will do its best to efficiently distribute the workload for the given number of MPI processes.
    If finer control is needed, please consult the documentation of [[gwpt_np_wpqbks]].

While the calculation is running, let's discuss the input file in more detail:

{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_7.abi %}

To activate the computation of the GWPT matrix elements, we use the
two variables [[optdriver]] and [[eph_task]]

```
optdriver 7  # Enter EPH driver.
eph_task 17  # GWPT computation.
```

Since we need a GSTORE for the ZPR only for specific $\kk$-points, we use
[[gstore_kfilter]] = "qprange" to select only the $\kk$-points associated with the band edges
and [[gstore_use_lgk]] = 1 to restrict the $\qq$-points to the IBZ_k.
These options are crucial to reduce the computational cost of the GWPT run.

Note that GWPT requires several external files in input.
The KS states are read from the WFK file via [[getwfk_filepath]].
This file defines the list of $\kk$-points in the e-ph matrix elements.
The value of [[ngkpt]], [[nshiftk]] and [[shiftk]] **must be consistent** with the ones used to generate the WFK file.
The number of bands in the $n'$ sum is given by [[nband]].
Clearly this value **cannot be greater** than the number of bands stored in the WFK file.

The DFPT KS potentials are read from [[getdvdb_filepath]], while the DFPT KS densities
are taken from [[getdrhodb_filepath]].
These two files define the list of $\qq$-points in the e-ph matrix elements.
This $\qq$-mesh must be identical to, or a submesh of, the $\kk$-mesh associated with the WFK file.

It is worth mentioning that it is possible to densify the $\qq$-mesh by using [[eph_ngqpt_fine]].
In a typical scenario, one generates a WFK file on a $\kk$-mesh much denser than the one used in the DFPT part,
and then use the Fourier interpolation of the DFPT potentials to reach a $\qq$-mesh that is equal or half the $kk$-mesh.
Further details on the interpolation of the DFPT scattering potentials are available on the [eph_intro page](eph_intro.md).

The screening is read from the SCR file specified with [[getscr_filepath]].
The cutoff energy in $W$ is given by [[ecuteps]], while [[ecutsigx]] defines
the cutoff energy for the exchange part of the self-energy.
Clearly, [[ecuteps]] **cannot be larger** than the value used in the previous screening calculation.
The SCR file defines the $\pp$-mesh for the integration over the transferred momenta in Eq.
This $\pp$-mesh must be identical to, or a submesh of, the $\kk$-mesh associated with the WFK file.
No interpolation in $\pp$-space is possible at present.
In a typical scenario, one works with a fixed reasonably-dense $\pp$-mesh, e.g. 6x6x6, uses a WFK file
on a $\kk$-mesh that is a multiple of 6x6x6 and uses [[eph_ngqpt_fine]] to reach the same BZ sampling as the
one used for electrons by using the interpolation of the DFPT scattering potentials.

Finally, the GWPT code needs to read the GS KS potential from the file specified with [[getpot_filepath]].
This file is used to build the GS Hamiltonian required by the Sternheimer solver.
The GS POT is produced at the end of the GS SCF cycle by setting [[prtpot]] to 1 (note that the default is 0).

The first-order derivative of the KS wavefunctions due to an atomic perturbation is computed on-the-fly
by solving the NSCF Sternheimer equation for each $n'$ band, atomic displacement, and quasi-momentum transfer.
There are two variables controlling the NSCF cycle:
[[nline]] defines the maximum number of NSCF iterations while [[tolwfr]] defines the stopping criterion.

!!! important

    [[tolwfr]] has a significant impact on the computational cost as more NSCF iterations of the
    Sternheimer equation are needed to reach small residuals.

The treatment of the frequency dependence in the GWPT matrix elements is governed by [[gwpt_wmode]].
By default, the GWPT matrix elements are computed at the energy of the incoming state $\varepsilon_\nk$.

Other variables worth mentioning here are [[zcut]] and [[elph2_imagden]].

!!! Important

    Techniques beyond the PPM are presently not available in the GWPT code.
    Also, only [[ppmodel]] 1 and 2 are presently supported.
    By default, we use the Godby-Needs model.


Now let's have a look at the output results:

{% dialog tests/tutorespfn/Refs/teph4zpr_gwpt_7.abo %}


## Computing the ZPR with GWPT e-ph matrix elements

In this section, we finally compute the ZPR of MgO at the GWPT level using the results produced previously.
As you will see the calculation is much faster than the previous step
as this run is essentially a post-processing of the data stored in the GSTORE file.

The price to pay is that there are few parameters that can be changed as this level.
In other words, the quality of the ZPR obtained here mainly depends on the density of the $\qq$- and $\pp$-meshes
used in the previous section, the number of empty states ($n'$ index) and the different cutoff energies.
In order to perform convergence studies, one should go back to the previous step, increase the relevant parameters
and monitor how the ZPR is affected by these settings.
This is left as an exercise to the reader...

After this preamble, let us start the calculation by issuing:

```sh
mpirun -n 4 abinit teph4zpr_gwpt_8.abi > teph4zpr_8.log 2> err &
```

with the following input file:

{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_8.abi %}

Let us now discuss the most relevant input variables:
To activate the computation of the e-ph self-energy from a GSTORE file we use
[[optdriver]], [[eph_task]] and [[getgstore_filepath]]

```
optdriver 7   # Enter EPH driver.
eph_task 24   # ZPR from GSTORE.
getgstore_filepath "teph4zpr_gwpt_7o_GSTORE.nc"
```

The temperature mesh in Kelvin is defined by [[tmesh]].
The imaginary shift in the denominator of the self-energy is given by [[zcut]].

We also use [[eph_stern]] 1 and [[getpot_filepath]] to activate the Sternheimer approach
to account for the contribution of the bands beyond the active space defined by [[nband]].

TODO: One should check that nband is consistent with nb_kq

!!! Important

  The Sternheimer method is exact if one is interested
  in the on-the-mass-shell corrections at the KS level in the adiabatic approximation.
  In the case of GWPT calculations, one assumes that the GWPT matrix elements are
  very close to the KS ones when $m$ > [[nband]].

Now let's have a look at the final results reported in the main output file:

{% dialog tests/tutorespfn/Refs/teph4zpr_gwpt_8.abo %}

One can use [[gstore_gname]] to select the kind of e-ph matrix elements that should
be read from the GSTORE.
In order to compute the ZPR with KS matrix elements, use [[gstore_gname]] = "gvals_ks"

!!! tip

  [[eph_ahc_type]] 0 can be used to use the adiabatic version of the Allen-Heine-Cardona equation.
  This is the version that should be used when comparing GWPT ZPR with finite-difference $G_0 W_0$
  calculations performed at fixed screening ignoring contributions beyond the rigid-ion approximation
  commonly used to deal with the Debye-Waller term.
