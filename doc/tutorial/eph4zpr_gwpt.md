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
and the description of the [gstore-based approach](../tutorial/gstore.md), before running these examples.

Also, you are kindly invited to read the [first GW tutorial](../tutorial/gw1.md) if you are not familiar
with the GW implementation in Abinit.
A brief description of the formalism and of the equations implemented in the
code can be found in the [[theory:mbt|GW_notes]].

This lesson should take about 2.0 hours.

## Formalism

In the GWPT method [[cite:Li2019]], the e-ph matrix elements are computed by replacing the
first-order change of the KS Hamiltonian due to a phonon with the variation of the $GW$ self-energy.

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


and  $\widetilde{W}$ is the frequency convolution of the screened Coulomb interaction that ensures the energy conservation of the plasmons:

\begin{align}
    \widetilde{W}_{n\mathbf{k-p,GG'}}(\varepsilon) \equiv & \int  \frac{\mathrm{d}\varepsilon' \, W_{\mathbf{p,GG'}}(\varepsilon')e^{i\eta\varepsilon'}}{\varepsilon-\varepsilon_{n\mathbf{k-p}}+ \varepsilon' +i\eta \sign(\varepsilon_{n\mathbf{k-p}}-\varepsilon^{\mathrm{F}})} \\
W_{\mathbf{p,GG'}}(\varepsilon') =& \frac{1}{V} \int_V \mathrm{d} \mathbf{r} \mathrm{d} \mathbf{r}^{\prime} \, W(\mathbf{r}, \mathbf{r}'; \varepsilon')\nonumber \\
    & \times  e^{-i(\mathbf{p}+\mathbf{G}) \cdot \mathbf{r}}  e^{i(\mathbf{p}+\mathbf{G}^{\prime}) \cdot \mathbf{r}^{\prime}},
    \label{eq:screeningGG}
\end{align}


Besides the summations over empty states $n'$, the equation requires the knowledge of the screened interaction $W$
as well as the (full) first-order derivative of the KS states
$\partial_{\kappa\alpha\mathbf{q}} \psi_{n^{\prime}}$
due to an atomic displacement of atom $\kappa$ along direction $\alpha$ modulated by the wavevector $\qq$.
Both terms can be computed by Abinit.
The screened interaction $W$ is computed in terms of a sum of states with [[optdriver]] 3, and the results
are stored in the SCR file.
The first-order derivative of the KS states, on the contrary, are computed on the fly by the GWPT subdriver
by solving a non-self-consistent (NSCF) Sternheimer equation as explained in the sections below.

Please note that ZPR computations at the GWPT level are still a field of active research,
especially in polar materials where additional long-range (LR) terms of many-body carachter appear in the e-ph matrix
elements.
In this tutorial, we won't be able to converge the calculation so we mainly focus on explaining the different
steps involved and the input parameters affecting the quality of the calculation and the predictive power.

<!--
A typical workflow for ZPR-GWPT requires the same step as the ones
[needed for KS-ZPR](../tutorial/eph_intro.md#typical_workflow_for_zpr)
plus additional computations for the screened interaction $W$.
-->

## Getting started

[TUTORIAL_README]

Before beginning, you might consider working in a different subdirectory as for the other tutorials.
Why not create Work_eph4zpr_gwpt in $ABI_TESTS/tutorespfn/Input?

```sh
cd $ABI_TESTS/tutorespfn/Input
mkdir Work_eph4zpr_gwpt
cd Work_eph4zpr_gwpt
```

In this tutorial, we prefer to focus on the use of the EPH code hence
we will be using **pre-computed** DDB and DFPT POT and DEN files to bypass the DFPT part.
We also provide a GS DEN.nc file to initialize the NSCF calculations
and a POT file with the GS KS potential required to solve the Sternheimer equation.

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

The first steps are similar to the ones used in the [eph4zpr tutorial](../tutorial/eph4zpr.md)
First of all, let's merge the partial DDB files with the command

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

Note that for GWPT we also need to merge the files with the first-order change of the density,
as we need to compute:

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
    Remember to use [[prtden]] 1 in the DFPT calculations.


## Computing the WFK files with empty states

At this point, we need to generate a WFK file with empty bands by performing a NSCF KS calculation
starting from a well-converged ground-state density.
This WFK file will then be used to compute W, the GW self-energy, and the GWPT matrix elements.

You may now run the NSCF calculation by issuing:

```sh
mpirun -n 4 abinit teph4zpr_gwpt_4.abi > teph4zpr_4.log 2> err &
```

with the input file given by:

{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_4.abi %}

Here, we use a 4x4x4 $\Gamma$-centered $\kk$-mesh and 110 bands.
Note the use of [[getden_filepath]] to read the DEN.nc file instead of [[getden]] or [[irdden]].

At this point, it is worth commenting about the use of [[nbdbuf]].
As mentioned in the documentation, **the highest energy states require more iterations to converge**.
To avoid wasting precious computing time, we use a buffer that is ~10% of [[nband]].
This trick significantly reduces the wall-time as the NSCF calculation completes
only when the first [[nband]] - [[nbdbuf]] states are converged within [[tolwfr]].
Obviously, one should not use the last [[nbdbuf]] states in the subsequent EPH calculation.

[HDIAGO_README]


## Computing the screened interaction W and QP corrections with one-shot GW

In this section, we use the WFK file generated in the previous section to compute the RPA polarizability
and the screened interaction $W$.
Note that here we generate a screening file with only two frequencies (default behaviour).
The SCR file will then be used to construct the plasmon-pole model (PPM) when computing
the self-energy using [[ppmodel]].

At this stage, one should perform convergence studies for [[nband]], [[ecuteps]], [[ecutsigx]],
and [[ngkpt]] to ensure that the GW results are reasonably well converged.
The converged parameters can then be reused in the subsequent GWPT calculation.

For the sake of conciseness, these convergence studies are omitted here, and reasonably
converged values are assumed in what follows.


To compute the SCR, execute e.g.:

```sh
mpirun -n 4 abinit teph4zpr_gwpt_5.abi > teph4zpr_5.log 2> err &
```

with the input file given by:

{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_5.abi %}


To compute the GW self-energy, execute e.g.:

```sh
mpirun -n 4 abinit teph4zpr_gwpt_6.abi > teph4zpr_6.log 2> err &
```

with the input file given by:

{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_6.abi %}

{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_6.abo %}


## Computing e-ph matrix elements with GWPT

GWPT computations require several external files in input.
In what follows, we describe all the files step by step and discuss the connection with
the parameters appearing in the equations.

You may want to start immediately the computation by issuing:

```sh
mpirun -n 4 abinit teph4zpr_7.abi > teph4zpr_7.log 2> err &
```

with the following input file:

{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_7.abi %}

that produces:

{% dialog tests/tutorespfn/Refs/teph4zpr_gwpt_7.abo %}

To activate the computation of the GWPT matrix elements, we use the two variables:
[[optdriver]] and [[eph_task]]

```
optdriver 7  # Enter EPH driver.
eph_task 17  # GWPT computation.
```

Since we plan to compute a GSTORE for the ZPR of the band gap, we use
[[gstore_kfilter]] = "qprange" to select only the $\kk$-points associated with the band edges
and [[gstore_use_lgk]] = 1 to restrict the $\qq$-points to the IBZ_k.
These options are crucial to reduce the computational cost of the GWPT part.

The KS states are read from the WFK file via [[getwfk_filepath]].
This file defines the list of $\kk$-points in the e-ph matrix elements.
The value of [[ngkpt]], [[nshiftk]] and [[shiftk]] **must be consistent** with the ones used to generate the WFK file.
The number of bands in the $n'$ sum is given by [[nband]].
Clearly this value **cannot be greater** than the number of bands stored in the WFK file.

The screening is read from the SCR file specified with [[getscr_filepath]].
The cutoff energy in $W$ is given by [[ecuteps]], while [[ecutsigx]] defines
the cutoff energy for the exchange part of the self-energy.
Note that [[ecuteps]] **cannot be larger** than the value used in the screening calculation.
The SCR file defines the $\pp$-mesh for the integration over transferred momenta in Eq.
This $\pp$-mesh must be identical to, or a submesh of, the $\kk$-mesh associated with the WFK file.
No interpolation in $\pp$-space is possible at this level.

The DFPT KS potentials are read from [[getdvdb_filepath]], while the DFPT KS densities
are taken from [[getdrhodb_filepath]].
These two files define the list of $\qq$-points in the GWPT matrix elements.
This $\qq$-mesh must be identical to, or a submesh of, the $\kk$-mesh associated with the WFK file.
Note that, in this case, it is possible to densify the $\qq$-mesh by using [[eph_ngqpt_fine]].
In a typical scenario, one generates a WFK file on a $\kk$-mesh much denser than the one used in the DFPT part
and then use the Fourier interpolation of the DFPT potentials to reach a $\qq$-mesh that is equal or half the $kk$-mesh.

Finally, the GWPT code needs to read the GS KS potential from the file specified with [[getpot_filepath]].
This file is produced at the end of the GS SCF cycle by setting [[prtpot]] to 1 (note that the default if 0).

The first order derivative of the KS wavefunctions due to an atomic perturbation is computed on-the-fly
by solving the NSCF Sternheimer equation.

There are two variables controlling the NSCF cycle:
[[nline]] defines the maximum number of iterations while [[tolwfr]] defines the stopping criterion.

!!! important

    [[tolwfr]] has a big impact on the computational cost as more NSCF iterations of the
    Sternheimer equation are needed to reach small residuals.

The treatment of the frequency dependence in the GWPT matrix elements is governed by [[gwpt_wmode]].
By default, the GWPT matrix elements are computed at the energy of the incoming state $\varepsilon_\nk$.

Other variables worth mentioning here are:

[[zcut]]
[[elph2_imagden]]

Also, [[ppmodel]] defines the kind of plasmon-pole approximation in the GWPT equation
By default we use the Godby-Needs model

!!! Important

    Techniques beyond the PPM are presently not available in the GWPT code.
    Also, only [[ppmodel]] 1 and 2 are presently supported.


## Computing the ZPR with GWPT

In this section, we can finally compute the ZPR of MgO at the GWPT level using the results produced previously.
It is important to understand that at this level of the calculations
there are few parameters that can be changed as this run is essentially
a post-processing of the e-ph matrix elements stored in the GSTORE file.

Start the calculation by issuing:

```sh
mpirun -n 4 abinit teph4zpr_gwpt_8.abi > teph4zpr_8.log 2> err &
```

with the following input file:


{% dialog tests/tutorespfn/Input/teph4zpr_gwpt_8.abi %}

```
optdriver 7
eph_task 24

gstore_gname "gvals"      # Use GWPT e-ph matrix elements from GSTORE (default)
eph_stern 1               # Activate Sternheimer to compute contribution given by states above nband
```

The temperature mesh is defined by [[tmesh]].
The imaginary shift in the denominator of the self-energy is given by [[zcut]].

!!! tip

    One can use [[gstore_gname]] to select the kind of e-ph matrix elements that should
    be read from the GSTORE.
    In order to compute the ZPR with KS matrix elements, use [[gstore_gname]] = "gvals_ks"


{% dialog tests/tutorespfn/Refs/teph4zpr_gwpt_8.abo %}









Finally, let us mention that [[eph_ahc_type]] 0 can be used
to use the adiabatic version of the Allen-Heine-Cardona equation to compute the ZPR.
This is the version that should be used when comparing with finite-difference GW calculations
at fixed screening.



