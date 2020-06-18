---
authors: GB, MG
---

# Phonon-limited mobility in AlAs

This tutorial shows how to compute phonon-limited carrier mobilities in semiconductors within
the relaxation time approximation.
It is assumed the user has already completed the two tutorials [RF1](rf1) and [RF2](rf2),
and that he/she is familiar with the calculation of ground state and response properties,
in particular phonons, Born effective charges and dielectric tensor.
It goes without saying that you should have read the [introduction page for the EPH code](eph_intro) 
before running these examples.

This lesson should take about 1 hour.

## Formalism

Before starting, it is worth summarizing the most important equations implemented in the code.
For a more detailed description of the ABINIT implementation, please consult [[cite:Brunin2020]].

Our goal is to find an approximated solution to the linearized
Boltzmann transport equation (BTE) [[cite:Ashcroft1976]] within the relaxation time approximation.
If what follows, we will be working within the so-called
self-energy relaxation time approximation (SERTA) [[cite:Giustino2017]].
The SERTA is more accurate than the constant relaxation time approximation (CRTA) as the
microscopic e-ph scattering is now included thus leading to linewidths that depend on the band index $n$
and the wavevector $\kk$.
Keep in mind, however, that the SERTA is still an approximation and that a more rigorous approach would require
to find the iterative solution of the BTE and/or the inclusion of many-body effects at different levels.
For a review of the different possible approaches see [[cite:Ponce2020]].

In the SERTA, the transport linewidth is given by
the imaginary part of the electron-phonon (e-ph) self-energy evaluated at the KS energy.
Only the Fan-Migdal (FM) part contributes to the linewidth as the Debye-Waller is Hermitian.
The linewidth of the electron state $n\kk$ due to the scattering with phonons are obtained from

\begin{equation}
\begin{split}
    \lim_{\eta \rightarrow 0^+} & \Im\{\Sigma^\FM_{n\kk}(\enk)\} =
                \pi \sum_{m,\nu} \int_\BZ \frac{d\qq}{\Omega_\BZ} |\gkkp|^2\\
                & \times \left[ (n_\qnu + f_{m\kk+\qq})
                                \delta(\enk - \emkq  + \wqnu) \right.\\
                & \left. + (n_\qnu + 1 - f_{m\kk+\qq})
                                \delta(\enk - \emkq  - \wqnu ) \right]
\end{split}
\label{eq:imagfanks_selfen}
\end{equation}

In this equation, $\nu$ is the phonon mode, $m$ the final electron state (after the scattering),
$n_\qnu$ is the Bose-Einstein occupation number, $f_{m\kk+\qq}$ is the Fermi-Dirac occupation function,
$\enk$ is the energy of the electron state, $\wqnu$ is the phonon frequency for the phonon wavevector $\qq$.
The integration is performed over the BZ for the phonon wavectors and $\gkkp$ is the e-ph matrix element.
The electron lifetime $\tau_{n\mathbf{k}}$ is inversely proportional to the linewidth:

\begin{align}
\frac{1}{\tau_{n\kk}} =
    2 \lim_{\eta \rightarrow 0^+} \Im\{\Sigma^\FM_{n\kk}(\enk)\}.
\label{eq:fanlifetime}
\end{align}

!!! important

    Note that this formalism does not take into account possibile contributions to the lifetime given by
    other scattering processes such as defects, ionized impurities in doped semiconductors, grain boundary scattering.
    These effects may be relevant depending on the system and/or the temperature under investigation
    but they are not treated in this tutorial.

The generalized transport coefficients are given by [[cite:Madsen2018]]

\begin{equation}
    \mcL^{(m)}_{\alpha\beta} =
    - \sum_n \int \frac{d\kk}{\Omega_\BZ} \vnka \vnkb\, \tau_{n\kk} 
    (\enk-\ef)^m
    \left.\frac{\partial f}{\partial\varepsilon}\right|_{\enk}
    \label{eq:transport_lc}
\end{equation}

where $\vnka$ is the $\alpha$-th Cartesian component of the matrix element the velocity operator

$$\vnk = \PDER{\enk}{\kk} = \langle \nk | \dfrac{\partial{H}}{\partial{\kk}} | \nk \rangle.$$

that can be computed with DFPT.
The generalized transport coefficients can be used to obtain different transport properties such as
the electrical conductivity, Peltier and Seebeck coefficients, and charge carrier contribution to the
thermal conductivity tensors [[cite:Madsen2018]].
The electrical conductivity tensor is given by

\begin{equation}
    \sigma_{\alpha\beta} =
    \frac{1}{\Omega} \mcL_{\alpha\beta}^{(0)} \label{eq:transport_sigma}
\end{equation}

and can be divided into hole and electron contributions

\begin{equation}
\sigma = n_e \mu_{e} + n_h \mu_{h} \label{eq:mobility}
\end{equation}

where $n_e$ and $n_h$ are the electron and hole concentrations in the conduction and valence bands respectively,
and $\mu_e$ and $\mu_h$ are the electron and hole mobilities,
which can be obtained by selecting the conduction or valences states $n$ in Eq. \eqref{eq:transport_lc}.

For electrons,

\begin{equation}
\begin{split}
    n_e = \sum_{n\in \text{CB}} \int \dfrac{d\kk}{\Omega_\BZ} f_{n\kk}, \\
    \mu_e = \dfrac{1}{n_e \Omega}\, \mcL_{n\in \text{CB}}^{(0)}
\end{split}
\end{equation}

where $n\in\text{CB}$ denotes states in the conduction bands. Similar expressions hold for holes.
At zero total carrier concentration, the Fermi level $\ef$ is located inside the band gap so that $n_e = n_h$.

A typical computation of mobilities requires different steps that are summarized 
in the [introduction page for the EPH code](eph_intro).
Here we only describe the e-ph related part, i.e the blue-box in the workflow presented in the previous page.
For this purpose, we use [[eph_task]] -4 to obtain only the imaginary part of the SE and explain other important
aspects in this tutorial.

<!--
  * Ground state computation to obtain the DEN and the WFK file
  * DFPT calculation for a set of $\qq$-points in the IBZ associated to a homogeneous mesh.
    Each DFPT run reads the WFK file produced in the previous step and produces **partial** DDB and POT1 files.
  * Merging of the partial DDB and POT1 files with *mrgddb* and *mrgdv*, respectively
  * Computation of the GS wavefunctions on a $\kk$-mesh that is much denser than the one
    used for the DFPT part as e-ph properties converge much slower than phonons.
  * Interpolation of the e-ph scattering potentials in $\qq$-space and computation of the lifetimes
    for the relevant $\kk$-points contributing to the mobility
  * Computation of transport properties

These steps can be summarized by the following graph:

![](eph4mob_assets/workflow.png ){: style="height:500px;width:400px"}
-->

All the results of the calculation are saved in netcdf format,
while the log and output files are used to report selected quantities, mainly for testing purposes.
Post-processing and visualisation tools are **not covered** in this tutorial.
Powerful tools based on python and matplotlib are provided by AbiPy.
See e.g. the README of [AbiPy](https://github.com/abinit/abipy)
and the [AbiPy tutorials](https://github.com/abinit/abitutorials).

[TUTORIAL_README]

## Ground state and phonons of fcc AlAs

*Before beginning, you might consider creating a different subdirectory to work in.
Why not create Work_eph4mob ?*

The file *teph4mob_1.in* is the input file for the first step (GS + DFPT perturbations).
Copy it to the working directory with:

```sh
cd $ABI_TUTORESPFN/Input
mkdir Work_eph4mob
cd Work_eph4mob
cp ../teph4mob_1.in .
```

{% dialog tests/tutorespfn/Input/teph4mob_1.in %}

!!! important

    Since ABINIT v9, the use of the files file is deprecated. The input file will be used instead,
    with an additional input variable listing the pseudopotentials.
    You should edit *teph4mob_1.in* to include the pseudopotentials, see [[pseudos]].
    This variable gives the different paths to the pseudopotentials.
    Use, for instance,

    ```sh
    pseudos "../../../Psps_for_tests/13al.981214.fhi, ../../../Psps_for_tests/33as.pspnc"
    ```

    to specify the relative path to the pseudos if you are following the conventions employed in this tutorials.
    You may need to adapt the paths depending on your *$ABI_PSPDIR*.

This tutorial starts with the DFPT calculation for all the $\qq$-points in the IBZ.
This step might be quite time-consuming so you can immediately start the job in background with:

```sh
abinit teph4mob_1.in > teph4mob_1.log 2> err &
```

The calculation is done for AlAs, the same crystalline material as in the first two DFPT tutorials.
<!-- Many input parameters are also quite similar. -->
<!-- as the EPH code can use symmetries can be used to reduce the number of atomic perturbations. -->
For more details about this first step, please refer to the first and second tutorials on DFPT.

!!! important

    Since AlAs is a **polar semiconductor**, we need to compute with DFPT the Born effective charges
    as well and the static dielectric tensor.
    These quantities are then used to treat the long-range (LR) part of the dynamical matrix in
    the Fourier interpolation of the phonon frequencies.
    As discussed in the EPH introduction. these quantities are also needed for the Fourier 
    interpolation of the DFPT potentials.

<!--
Only the (partial) DDB and POT files produced at the end of the DFPT run
are needed to perform e-ph calculation.
The files containing the first order wavefunctions (*1WF*) due to an atomic perturbation are not needed.
-->

## Merging the derivative databases and potentials

Once the DFPT calculation is completed, use *mrgddb* to merge the eight partial DDB files
corresponding to datasets 3-10 of *teph4mob_1*.
These partial DDBs contain the dynamical matrices for the
8 $\qq$-points in the IBZ, as well as the dielectric tensor and the Born effective charges.
Name the new DDB file *teph4mob_2_DDB*.

File *\$ABI_TUTORESPFN/Input/teph4mob_2.in* is an example of input file for *mrgddb*.

{% dialog tests/tutorespfn/Input/teph4mob_2.in %}

Copy the file in the *Work_eph4mob* directory, and run *mrgddb* using:

```sh
mrgddb < teph4mob_2.in
```

!!! tip

    Alternatively, one can specify the name of the output DDB and the list of DDB files
    to be merged via command line arguments.
    This approach is quite handy especially if used in conjuction with shell globbing and the "star" syntax:

    ```sh
    mrgddb teph4mob_2_DDB teph4mob_1o_DS*_DDB
    ```

    Use **mrgddb --help** to access the documentation.

Now use the *mrgdv* tool to merge the 29 DFPT POT files corresponding to datasets 3-10 of *teph4mob_1*.
Name the new file *teph4mob_3_DVDB*.

File *\$ABI_TUTORESPFN/Input/teph4mob_3.in* is an example of input file for *mrgdv*.

{% dialog tests/tutorespfn/Input/teph4mob_3.in %}

You can copy it in the *Work_eph4mob* directory, and then merge the files with:

```sh
mrgdv < teph4mob_3.in
```

!!! tip

    Alternatively, one can use the command line.

    ```sh
    mrgdv merge teph4mob_3_DVDB teph4mob_1o_DS*_POT*
    ```

    Use **mrgdv --help** to access the documentation.

We now have all the phonon-related files needed to compute the mobility.
The DDB will be used to Fourier interpolate the phonon frequencies on an **arbitrarily** dense $\qq$-mesh while
the DVDB will be used to Fourier interpolate the DFPT scattering potentials [[cite:Brunin2020]].
The only ingredient left is the WFK file with the GS wavefunctions on the dense $\kk$-mesh.

!!! important

    In real computations, you should always compute the electronic band structure along a $\kk$-path
    to have a qualitative understanding of the band dispersion, the position of the band edges,
    and the value of the band gap(s).
    Note also that there are several parts of the EPH code in which it is assumed that no vibrational instability
    is present so you should always look at the phonon spectrum computed by the code.
    Don't expect to get meaningful results if imaginary phonon frequencies (a.k.a **negative frequencies**) are present.

<!--
In principle, the code can deal with small instabilties for the acoustic modes around $\Gamma$
but this does not mean you should trust your results.
-->


## Calculation of the dense WFK file

In order to compute transport properties, we need a $\kk$-mesh that is dense enough to
sample the electron (hole) pockets.
<!--
to solve the Boltzmann Transport Equation (BTE) in the SERTA on a relatively dense $\kk$-mesh to have enough wavevectors
inside the electron (hole) pockets.
You will compute the electron lifetimes and group velocities on this dense mesh.
We will therefore need the wavefunctions on this dense mesh.
-->
Note that in this tutorial, a single dense $\kk$-mesh is used.
However, the mobility strongly depends on this mesh and a convergence study should be performed by 
increasing the $\kk$-mesh density,
as well as the $\qq$-mesh used for the integration of the imaginary part of the self-energy.
This study is explained later and left to the user as excercise.
It should be clear even that systems with small effective masses (e.g GaAs)
require denser homogeneous $\kk$-meshes and are more difficult to converge.

The computation of the dense WFK file is similar to a NSCF band structure computation.
The main difference is that we need wavefunctions on a $\kk$-mesh instead of a $\kk$-path
because these wavevectors will be used to evaluate integrals in the BZ.

The file *\$ABI_TUTORESPFN/Input/teph4mob_4.in* is an example of such computation.

{% dialog tests/tutorespfn/Input/teph4mob_4.in %}

It consists of two parts: the first one (dataset 1) computes the GS wavefunctions,
and the second one (datasets 2-3) computes the dense WFK that will be used to evaluate the mobility.
We also compute a denser WFK file that will be used with the double-grid method explained later.

As we want to compute the mobility of electrons in the conduction band, 
we need to consider conduction bands in the computation of the WFK ([[nband]] = 8).

Copy the file in the *Work_eph4mob* directory, and run ABINIT:

```sh
abinit teph4mob_4.in > teph4mob_4.log 2> err &
```

## Calculation of the mobility

The computation of the mobility requires different convergence studies.
We will explain them and their need in the following.
Let us first explain the different parameters in a standard mobility computation.
The file *\$ABI_TUTORESPFN/Input/teph4mob_5.in* is an example of such computation.

{% dialog tests/tutorespfn/Input/teph4mob_5.in %}

First of all, we need to read the WFK, the DDB and the DVDB files obtained previously.
Since it is not possible to run a mobility calculation with a single input file
with datasets, we use strings to specify the path to the different files:

```sh
getwfk_filepath "teph4mob_4o_DS2_WFK"
getddb_filepath "teph4mob_2_DDB"
getdvdb_filepath "teph4mob_3_DVDB"
```

Now copy the input file in the *Work_eph4mob* directory, and run the code with:

```sh
abinit teph4mob_5.in > teph4mob_5.log 2> err &
```

The job should take $\sim$15 seconds on a recent CPU.

Let's discuss the meaning of the e-ph variables in more details:

* [[optdriver]] 7 is required to activate the EPH driver

* [[eph_task]] -4 tells ABINIT that we only need the imaginary part
  of the e-ph self-energy at the KS energy, which directly gives the electron lifetimes.

* The homogeneous mesh corresponding to the WFK file is specified by [[ngkpt]] 24 24 24.
  The code will stop with an error if [[ngkpt]] is not the same as the one corresponding to the 
  dense WFK file.

* [[occopt]] 3 is required to correctly compute the
  location of the Fermi level using the Fermi-Dirac occupation function as we are dealing with the
  physical temperature and not a Fictitious broadening for integration purposes.

* [[ddb_ngqpt]] defined the initial grid used for the DFPT computation (4×4×4 in this example)

* [[eph_ngqpt_fine]] defines the dense $\qq$-mesh onto which the scattering potentials are interpolated
  and the e-ph matrix elements are computed. 


!!! warning

    [[ngkpt]] and [[eph_ngqpt_fine]] should be commensurate.
    More specifically, [[ngkpt]] must be a multiple of the $\qq$-mesh ([[eph_ngqpt_fine]])
    because the WFK should contain all the $\kk$- and $\qq$-points.
    In most cases, [[ngkpt]] = [[eph_ngqpt_fine]]. 
    It is however possible to use fewer $\qq$-points.
    Note that [[ngkpt]] does not necessarily correspond to the $\kk$-mesh used for the computation of
    transport quantities, see the following discussion.
    
<!--In this tutorial, we will use the same dense $\kk$- and $\qq$-meshes.
As a rule of thumb, a $\qq$-mesh twice as dense in each direction as the $\kk$-mesh,
is needed to achieve fast convergence of the integrals [[cite:Brunin2020]].
In this case, [[ngkpt]] = [[eph_ngqpt_fine]], but the use of [[sigma_ngkpt]]
allows to downsample the $\kk$-mesh used for the integration and it should be set to half
the values of [[ngkpt]].
An example will be given in this tutorial.
Possible exceptions are systems with very small effective masses (e.g. GaAs) in which
a very dense $\kk$-sampling is needed to sample the electron (hole) pocket.
In this case, using the same sampling for electrons and phonons may be enough to converge.-->

We will use the tetrahedron integration method [[cite:Blochl1994]] to obtain the lifetimes (integration
over the $\qq$-mesh). This allows to efficiently filter out the $\qq$-points that do not contribute
to the lifetimes. Indeed, only a small fraction of the $\qq$-points belonging to the $\qq$-mesh
ensure energy and momentum conservation for a given $\kk$-point.
All the other $\qq$-points do not need to be considered and can be filtered out.
The use of the tetrahedron method is automatically activated when [[eph_task]] is set to -4.
It is possible to change this behaviour by using [[eph_intmeth]] albeit not recommended 
as the calculation will become significantly slower.

The list of temperatures for which the mobility is computed is specified by [[tmesh]].
The carrier concentration is specified by [[eph_doping]] is |electron_charge| / cm^3 units so negative for electron-doping,
positive for hole doping.
To obtain results that are representative of the intrinsic mobility,
we suggest to use a very small number, for instance $10^{15}$ to $10^{18}$ electrons per cm$^3$.

!!! tip

    The computational cost increases with the number of temperatures.
    For the initial convergence studies, we suggest to start from a relatively small number
    of temperatures covering the region of interest. 
    The T-mesh can be densified aftwerwards when  converged parameters are found.
    Note that the convergence might be different depending on the temperature.
    Low temperatures are more difficult to converge.

The [[sigma_erange]] variable defines the energy window, below the VBM and above the
CBM where the lifetimes will be computed.
Since the BTE contains a derivative of the Fermi-Dirac occupation function centered on the Fermi level,
it is possible to filter the $\kk$-points that will contribute to the mobility and compute
the lifetimes for these $\kk$-points only. Indeed, the value of the derivative decreases rapidly
as we go further from the Fermi level. Only the states close to the band edges contribute.
This additional filtering technique allows one to compute only a few percents of all the lifetimes.
This variable should be subject to a convergence study, as explained in the next section.
<!--
A value of 0.2 eV represents a good starting point for further analysis.
Hopefully, in the next version this parameter will be automatically computed by the code.
-->

We can now examine the log file in detail.
After the standard output of the input variables,
the code reports the different parameters for the long-range potentials: the Born effective charges,
the dielectric constant, and the quadrupolar tensor.
Make sure to have all of them in order to have an
accurate interpolation of the scattering potentials, see discussion in [[cite:Brunin2020]].

*Note: for the moment (|today|), the computation of the quadrupole is not available in the public version,
and you should have the following in the log file:*

```sh
Have dielectric tensor: yes
Have Born effective charges: yes
Have quadrupoles: no
```

The code then outputs different information. The first of them is the location of the Fermi level
that will be used to compute the lifetimes. You can check that it is far enough from the band
edges so that the computed mobility can be considered as intrinsic.

```sh
Valence Maximum:   2.3564 (eV) at: [ 0.0000E+00,  0.0000E+00,  0.0000E+00]
Conduction minimum:   3.5254 (eV) at: [ 5.0000E-01,  5.0000E-01,  0.0000E+00]
Fermi level:  3.0878 (eV)
```

ABINIT also finds the list of $\kk$-point belonging to the dense mesh that
are located within the energy window:

```sh
Found 3 k-points within erange:  0.000  0.150  (eV)
min(nbcalc_ks): 1 MAX(nbcalc_ks): 1
```

Over the initial 413 $\kk$-points in the IBZ, only 3 will be computed!
ABINIT then reads the WFK file and interpolates the potentials to compute the e-ph matrix elements.
The use of the tetrahedron method allows for the filtering of the $\qq$-points:

```sh
qpoints_oracle: calculation of tau_nk will need: 39 q-points in the IBZ. (nqibz_eff / nqibz):   9.4 [%]
```

Again, this leads to a significant speed-up of the computation.
Once this is done, the code starts looping over the 3 $\kk$-points for which the lifetimes are needed.

```sh
Computing self-energy matrix elements for k-point: [ 4.5833E-01,  4.5833E-01,  0.0000E+00] [ 1 / 3 ]
```

You can find various information for each $\kk$-point, such as:

* the total number of $\qq$-points belonging to the little group of $\kk$ (called IBZ(k) in the code),
  i.e. the irreducible $\qq$-points for this specific $\kk$-point
* the number of $\qq$-point in the IBZ(k) contributing to the imaginary part of $\Sigma_\nk$ 
  (in most cases, this number will be much smaller than the total number of $\qq$-points in the IBZ(k))
* the wall-time each step takes.

Finally, we have the results for the lifetimes (TAU) in the *teph4mob_5.out* file:

```sh
K-point: [ 4.5833E-01,  4.5833E-01,  0.0000E+00], T=    5.0 [K]
    B    eKS    SE2(eKS)  TAU(eKS)  DeKS
	5   3.573    0.000  31553.2    0.000
```

Only the first temperature is printed in the output file, but all the results can be found in the SIGEPH.nc file.

!!! tip

    With |AbiPy|, one can easily have access to all the data of the computation. For instance, one can plot the
    electron linewidths:

    ```sh
    abiopen.py teph4mob_5o_SIGEPH.nc --expose
    ```

![](eph4mob_assets/linewidths.png)

At the end of the *.out* and *.log* files, the mobility is printed:

```sh
Temperature [K]             e/h density [cm^-3]          e/h mobility [cm^2/Vs]
            5.00        0.23E+17        0.00E+00            0.00            0.00
           64.00        0.23E+17        0.00E+00           38.72            0.00
          123.00        0.23E+17        0.00E+00          346.73            0.00
          182.00        0.23E+17        0.00E+00          424.55            0.00
          241.00        0.23E+17        0.00E+00          420.01            0.00
          300.00        0.23E+17        0.00E+00          364.31            0.00
```

The temperature is first given then the electron and hole densities followed by electron and hole mobilities.
In this computation, we consider only electrons and this explains why the values for holes are zero.
Note that the transport driver is automatically executed after the EPH run.
You can run the transport driver in standalone mode by setting [[eph_task]] 7, 
provided you already have the lifetimes in a SIGEPH.nc file,
<!-- This task can be performed only in serial and is very fast. -->

Now that you know how to obtain the mobility in a semiconductor for given k- and q-meshes,
we can give more details about convergence and additional tricks that can be used
to decrease the computational cost.

### Convergence w.r.t. the energy range

The first convergence study to be performed is to determine the energy range around the band edges
to be used for the computation of $\tau_{n\kk}$.
We can do that by performing the same mobility computation (same $\kk$- and $\qq$-mesh)
for an increasing energy window specified by [[sigma_erange]].

!!! important

    The code can compute both electron and hole mobilities in a single run 
    but this is not the recommended procedure as the $\qq$-point filtering is expected to be less efficient.
    Moreover electrons and holes may require a different $\kk$-sampling to convergence depending on the dispersion
    of the bands. So we suggest to run the calculations with different input files.

The file *$\$ABI_TUTORESPFN/Input/teph4mob_6.in* is an example of such computation.

{% dialog tests/tutorespfn/Input/teph4mob_6.in %}

Copy the input file in the *Work_eph4mob* directory, and run ABINIT:

```sh
abinit teph4mob_6.in > teph4mob_6.log 2> err &
```

This run should take a few minutes.

We can now analyze the variation of the mobility with respect to [[sigma_erange]].
Once enough $\kk$-points are included, the mobility should not vary, and computing $\kk$-points further in the bands will
only increase the computation time.
Note that you should perform this convergence study with a $\kk$-mesh that is already dense enough to
capture the band dispersion correctly.
In this case, we are using a 24×24×24 mesh, which is not very dense for such computations.
This means that, when increasing [[sigma_erange]], sometimes
no additional $\kk$-point is considered.
It is the case here for the first 3 datasets (3 $\kk$-points), and the last two datasets (6 $\kk$-points).
If a finer mesh was used, the number of $\kk$-points would have increased in a more monotonic way.

<!-- 
TODO: Discuss how to dope the system! 
Other quantities (Seebeck etc may have a different convergence behaviour
-->

### Convergence w.r.t. the k- and q-meshes

Once the energy window has been set, we can start to converge the mobility with respect to the
dense $\kk$- and $\qq$-meshes.
As a rule of thumb, a $\qq$-mesh twice as dense in each direction as the $\kk$-mesh,
is needed to achieve fast convergence of the integrals [[cite:Brunin2020]].
In this case, [[ngkpt]] = [[eph_ngqpt_fine]], but the use of [[sigma_ngkpt]]
allows to downsample the $\kk$-mesh used for the integration and it should be set to half
the values of [[ngkpt]].
Possible exceptions are systems with very small effective masses (e.g. GaAs) in which
a very dense $\kk$-sampling is needed to sample the electron (hole) pocket.
In this case, using the same sampling for electrons and phonons may be enough to converge.

The previous computations used 24×24×24 $\kk$- and $\qq$-meshes. This is far from convergence.
In silicon, for instance, a 45×45×45 $\kk$-mesh and 90×90×90 $\qq$-mesh are needed
to reach convergence within 5%.

In order to compute the mobility with a $\qq$-mesh twice as dense as the $\kk$-mesh, 
there are two possible approaches.
Let us take the previous example of silicon.

1. Run a computation with:

       * [[ngkpt]] 90 90 90,
       * [[eph_ngqpt_fine]] 90 90 90,
       * [[sigma_ngkpt]] 45 45 45.

   Using [[sigma_ngkpt]] will select the $\kk$-points belonging to the 45×45×45 mesh, but each lifetime will be computed
   with a 90×90×90 $\qq$-mesh.

2. Run a computation with:

      * [[ngkpt]] 90 90 90,
      * [[eph_ngqpt_fine]] 90 90 90,
      * [[sigma_ngkpt]] 90 90 90.

   In this way, you have the mobility with 90×90×90 $\kk$- and $\qq$-meshes. You can then run again the transport driver only,
   by setting [[eph_task]] 7, and setting [[transport_ngkpt]] 45 45 45. This will downsample the $\kk$-mesh used
   in the computation of the mobility. This second option has the advantage that it delivers two mobilities
   (useful for convergence studies), but it requires more computational time if only the 45×45×45/90×90×90 results are needed.

You can run again the previous input files by densifying the different meshes.
For the densest grids, you might need to run with multiple MPI processes.
You should obtain something like this, for $T$ = 300 K:

![](eph4mob_assets/teph4mob_conv.png)


Note that in order to get sensible results, one should use a denser DFPT $\qq$-mesh (around 9×9×9),
and a larger energy cutoff for the planewave basis set ([[ecut]]).
The inputs of this tutorial have been tuned so that the computations are quite fast,
but they are quite far from convergence.
In real studies, you should perform convergence tests to find suitable parameters.

### Double-grid technique

Another possibility to improve the results without increasing the computation time significantly
is the double-grid technique.
In this case, a coarse sampling is used for the $\kk$-mesh and the $\qq$-mesh for the e-ph matrix elements,
but a finer mesh is used for the phonon absorption-emission terms.
This technique allows one to better capture these processes, while computing the matrix elements on a coarser mesh.
The efficiency of this method depends
on the polar divergence of the matrix elements: if this divergence is very difficult to capture numerically,
the coarse $\qq$-mesh for the e-ph matrix elements will have to be denser.

The double-grid technique requires a second WFK file, containing the KS eigenvalues on the fine mesh.
You can specify the path to the fine WFK file using [[getwfkfine_filepath]]:

```sh
getwfkfine_filepath "teph4mob_4o_DS3_WFK"
```

The file *$\$ABI_TUTORESPFN/Input/teph4mob_7.in* is an example of such computation.

{% dialog tests/tutorespfn/Input/teph4mob_7.in %}

Copy the input file in the *Work_eph4mob* directory, and run ABINIT:

```sh
abinit teph4mob_7.in > teph4mob_7.log 2> err &
```

In the log file, you will now find information about the double-grid method:

```sh
coarse:                24          24          24
fine:                  48          48          48
```

The mobility obtained, at 300 K, is 157.59 cm$^2$/V/s.
Using a 48×48×48 $\qq$-mesh for the matrix elements as well would give 96.38.
The result is indeed improved, since using a 24×24×24 mesh for everything gives 364.31.
You can also use a finer mesh, but always a multiple of the initial coarse mesh
(in this case, 72×72×72, 96×96×96, etc).
However, we found that there is very little use to go beyond a mesh three times as dense as the coarse one.
Using a 72×72×72 fine mesh for the energies gives a mobility of 149.87 cm$^2$/V/s,
and a 96×96×96 mesh leads to 146.24 cm$^2$/V/s: the improvement is indeed rather limited.


### In-place restart

All the results of the calculation are stored in a single SIGEPH.nc file
for all the $\kk$-points (and spins) considered.
The list of $\kk$-points is initialized at the beginning of the calculation and an internal table
in the netcdf file stores the status of each $\kk$-point (whether it has been computed or not).
This means that calculations that get killed by the resource manager due to time limit can use
the SIGEPH file to perform an automatic in-place restart.
Just set [[eph_restart]] to 1 in the input file and rerun the job (there is no harm in setting it 
from the start).

!!! important

    Restart capabilities are not supported in multi-dataset mode!

### Transport calculation from SIGEPH.nc

The routine that computes carrier mobilites is automatically invoked when [[eph_task]] -4 is used
and a *TRANSPORT.nc* file with the final results is produced.
There are however cases in which one would like to compute mobilities from an already existing
SIGEPH.nc file without performing a full self-energy calculation.
In this case, one can use [[eph_task]] 7 and specify the name of the SIGEPH.nc file with [[getsigeph_filepath]].
The advanced input variable [[transport_ngkpt]] can be use to downsample the $\kk$-mesh used to evaluate the mobility integral.

### MPI parallelism and memory requirements

There are five different MPI levels that can be used to distribute the workload
and the most memory-demanding data structures.
By default, the code tries to reach some compromise between memory requirements and time to solution
by activating the parallelism over $\qq$-points if no other input is provided from the user.
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
Use these additional levels if the memory requirements are under control
and you want to boost the calculation.

### How to reduce memory requirements

As mentioned above, the memory should scale with the number of MPI processors used for the $\qq$-point and
the perturbation distribution.
However, there might be tricky systems in which you start to experience memory shortage that
prevents you from running with several MPI processes.
This problems should show up for very dense $\kk$/$\qq$ meshes.
As a rule of thumb, calculations with meshes denser than e.g 200x200x200 start to be very memory demanding
and become much slower because several algorithms and tables related to the BZ sampling will start to dominate.

The code uses an internal cache to store the DFPT potentials in the dense IBZ.
The size of the cache is defined by [[dvdb_qcache_mb]] whose default value is 1024 Mb.
The size of this buffer decreases with the number of MPI processese used for the qpt-level.
You can save some space by decreasing this value at the price of a global slow down of the calculation.

The code allocates a relatively small buffer to store the Bloch states involved in transport but unfortunately
the $\kk$-points are not easy to distribute with MPI.
To reduce the size of this part, one may opt for an internal buffer in single precision.
This option is enabled by using `enable_gw_dpc="no"` at configure time (note that this is the default behaviour).

If these tricks do not solve your problem, consider using OpenMP threads.
The code is not optimized for OpenMP but a few threads can be useful to avoid replicating memory at the MPI level.
As a rule of thumb 2-4 OpenMP threads should be OK provided you link with threaded FFT and BLAS libraries.

Last but not least, do not use datasets: large arrays allocated for $\kk$-points and the size depends on [[ndtset]].
Never ever use multiple datasets for big EPH calculations. You have been warned!

### How to compute only the $\kk$-points close to the band edges

<!-- part of the discussion can be moved to the eph_intro as SKW will be used also in phgamma -->

As we have already seen in the previous sections, a relatively small number of $\kk$-points
close to the band edges is usually sufficient to converge mobilities.
Yet, in the NSCF run, we computed a WFK file for all the $\kk$-points of the dense IBZ
hence we spent a lot of resources to compute and store states that are not
needed for phonon-limited mobilities.

In principle, it is possible to restrict the NSCF calculation to the relevant $\kk$-points
provided we have a cheap and good-enough method to predict whether the wavevector
is inside the energy window without solving the KS eigevalue problem exactly.
For example, one can use the star-function interpolation by Shankland-Koelling-Wood (SKW)
[[cite:Shankland1971]], [[cite:Koelling1986]], [[cite:Madsen2006]], [[cite:Madsen2018]]
which requires as input a set of eigenvalues in the IBZ and a single parameter defining the basis set for the interpolation.
There are several technical problems that should be addressed at the level of the internal implementation but
the idea is relatively simple and goes as follows:

1. Compute the KS eigenvalues on a relatively coarse $\kk$-mesh in the IBZ
2. Use this *coarse* WFK file to interpolate the eigenvalues on a much denser $\kk$-mesh specified by the user.
3. Find the wavevectors of the dense mesh inside an energy window specified by the user and
   store the list of $\kk$-points in a external file.
4. Use this external file to run a NSCF calculation only for these $\kk$-points.
   At the end of the NSCF job, ABINIT will produce a **customized** WFK file on the dense mesh that
   can be used to run calculations for phonon-limited mobilities

An example will help clarify.
Suppose we have computed a WFK file with a NSCF run using a 16x16x16 $\kk$-mesh (let's call it *161616_WFK*)
and we want to compute mobilites with the much denser 64x64x64 $\kk$-mesh.
In this case, use [[optdriver]] = 8 with [[wfk_task]] = "wfk_kpts_erange" to read
the WFK file specified by [[getwfk_filepath]], find the wavevectors belonging to the [[sigma_ngkpt]] $\kk$-mesh
inside the energy window defined by [[sigma_erange]] and produce a *KERANGE* file.
The parameters defining the SKW interpolation are specified by [[einterp]].
A typical input file for this step will look like:

```sh
optdriver 8
wfk_task "wfk_kpts_erange"
getwfk_filepath "161616_WFK"

# Define fine k-mesh for interpolation
sigma_ngkpt   64 64 64
sigma_nshiftk 1
sigma_shiftk  0 0 0

sigma_erange 0.0 0.2 eV   # Select kpts in fine mesh within this energy window.
einterp 1 5 0 0           # Parameters for star-function interpolation
```

This step produces a *KERANGE.nc* file (let's call it *out_KERANGE*) that can be used
via the [[getkerange_filepath]] variable as a starting point to perfom a NSCF run with the following variables:

```sh
getkerange_filepath "out_KERANGE.nc"
getden_filepath "161616_DEN"
getwfk_filepath "161616_WFK"    # Init GS wavefunctions from this file (optional)
iscf  -2
tolwfr 1e-18
kptopt 0                        # Important

# These variables must be consistent with the values of
# sigma_ngkpt, sigma_shiftk settings used in the previous step
ngkpt    64 64 64
nshiftk  1
shiftk   0.0 0.0 0.0
```

This calculation will produce a *customized* Fortran WFK file with [[ngkpt]] = 64 64 64 in which
only the states listed in the KERANGE file have been computed.
This WKF file can then be used in the EPH code to compute mobilites.
<!-- to compute mobilities with an energy window that shall not be greater than the one specifed in step 3.-->

For further examples see [[test:v9_57]], and [[test:v9_61]].
Also note that the two tests cannot be executed in multidataset mode with a single input file.
Note, however, that the quality of the interpolation depends on the initial coarse $\kk$-mesh.
So we recommended, to look at the interpolant.
It is also a good idea to use an energy window that is larger than the one employed for the mobility.
<!--as well as the position of the band edges-->
