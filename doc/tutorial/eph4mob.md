---
authors: GB, MG
---

# Phonon-limited mobility in AlAs

This tutorial shows how to compute phonon-limited carrier mobilities in semiconductors within
the relaxation time approximation.
It is assumed the user has already completed the Tutorials [RF1](rf1) and [RF2](rf2),
and that he/she is familiar with the calculation of ground state and response properties, 
in particular phonons, Born effective charges and dielectric tensor.

This lesson should take about 1 hour. 

## Formalism

Before starting, it is worth summarizing the most important equations implemented in the code.
For a more detailed description of the ABINIT implementation, please consult [[cite:Brunin2020]].

In this article, we focus on the solution of the linearized Boltzmann transport formulation [[cite:Ashcroft1976]]
within the relaxation time approximation (RTA).

If what follows, we will be working within the 
self-energy relaxation time approximation (SERTA) [[cite:Giustino2017]].
In the $\eta \rightarrow 0^+$ limit, 
the imaginary part of the electron-phonon (eph) self-energy evaluated at the KS energy is given by 

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

and corresponds to the linewidth of the electron state $n\kk$ due to the scattering with phonons.
The electron lifetime $\tau_{n\mathbf{k}}$ is inversely proportional to the linewidth:

\begin{align}
\frac{1}{\tau_{n\kk}} =
    2 \lim_{\eta \rightarrow 0^+} \Im\{\Sigma^\FM_{n\kk}(\enk)\}.
\label{eq:fanlifetime}
\end{align}

The generalized transport coefficients are given by [[cite:Madsen2018]]

\begin{equation}
    \begin{split}
    \mcL^{(m)}_{\alpha\beta} = 
    - \sum_n \int \frac{d\kk}{\Omega_\BZ} \vnka \vnkb\, \tau_{n\kk} & \\
    \times (\enk-\ef)^m 
    \left.\frac{\partial f}{\partial\varepsilon}\right|_{\enk}
    \label{eq:transport_lc}
    \end{split}
\end{equation}

where $\vnka$ is the $\alpha$-th component of the matrix element $\vnk$ of the electron velocity operator. 

The generalized transport coefficients can be used to obtain different transport properties such as
the electrical conductivity, Peltier and Seebeck coefficients, and charge carrier contribution to the
thermal conductivity tensors [[cite:Madsen2018]]
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
which can be obtained by selecting the conduction or valences states $n$ in Eq.~\eqref{eq:transport_lc}.

For electrons,

\begin{equation}
\begin{split}
    n_e = \sum_{n\in \text{CB}} \int \dfrac{d\kk}{\Omega_\BZ} f_{n\kk}, \\
    \mu_e = \dfrac{1}{n_e \Omega}\, \mcL_{n\in \text{CB}}^{(0)}
\end{split}
\end{equation}

where $n\in\text{CB}$ denotes states in the conduction bands. Similar expressions hold for holes.
At zero total carrier concentration, the Fermi level $\ef$ is located inside the band gap so that $n_e = n_h$.

A typical computation requires different steps:

  * Ground state computation to obtain the DEN and the WFK files
  * DFPT calculation for a set of $\qq$-points in the IBZ associated to a homogeneous mesh
  * Merging of the partial DDB and POT files
  * Computation of the GS wavefunctions on a dense $\kk$-mesh
  * Interpolation of the e-ph scattering potentials and computation of the electron lifetimes
    for the relevant $\kk$-points contributing to the mobility
  * Computation of transport properties

These steps can be summarized by the following graph:

![](eph4mob_assets/workflow.png )


!!! important

    Note that all these capabilities are integrated directly in ABINIT.
    This implementation (henceforth refered to as the **EPH code**) differs from the one available in ANADDB: 
    the anaddb version acts as a direct post-processing of the e-ph matrix elements computed in the DFPT part 
    whereas the EPH code interfaced with ABINIT computes the e-ph matrix elements directly using 
    the GS WFK and the DFPT potentials stored in the DVDB file.
    In a nutshell, the EPH code is more scalable and flexible as the $\qq$-sampling can be easily changed 
    at run time while the anaddb implementation easily supports advanced features such as PAW.
    For further information about the difference between the two approaches, see [[cite:Gonze2019]]

Note that all the results of the calculation are saved in netcdf format,
while the log and output files are used to report selected quantities mainly for testing purposes.
Post-processing and visualisation tools are NOT covered in this tutorial.
Powerful tools based on python and matplotlib are provided by Abipy
See the README of [Abipy](https://github.com/abinit/abipy)
and the [Abipy tutorials](https://github.com/abinit/abitutorials).

[TUTORIAL_README]

## Ground state and phonons of fcc AlAs

*Before beginning, you might consider creating a different subdirectory to work in. 
Why not create Work_eph4mob ?*

The file *teph4mob_1.in* is the input file for the first step (GS + DFPT perturbations).
You can copy it to the working directory with:

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

This tutorial starts with the DFPT calculation for all the $\qq$-points in the IBZ that might be quite time-consuming.
You can immediately start the job in background with:

```sh
abinit teph4mob_1.in > teph4mob_1.log 2> err &
```

The calculation is done for AlAs, the same crystalline material as in the first two DFPT tutorials.
Many input parameters are also quite similar. 
For more details about this first step, please refer to the first and second tutorials on DFPT.

Note this additional important points.
Since AlAs is a **polar semiconductor**, we need to compute with DFPT the Born effective charges 
as well and the static dielectric tensor
These quantities are then used to treat the long-range part of the dynamical matrix in 
the Fourier interpolation of the phonon frequencies. 
We will see that these quantities are also needed in the Fourier interpolation of the DFPT potential.


Only the (partial) DDB and POT files produced at the end of the DFPT run
are needed to perform e-ph calculation.
The files containing the first order wavefunctions (*1WF*) due to an atomic perturbatio are not needed.
As these files are quite larger and the overall space on disk scales as nq * 3 * natom, we suggest
to avoid the output of the DFPT wavefunctions by using [[prtwf]] = -1
In this case, the DFPT code writens the 1WF only if the DFPT SCF cycle is not converged so that one can restart
from it if needed.

## Merging the derivative databases and potentials

Once the DFPT calculation is completed, use *mrgddb* to merge the eight partial DDB files 
corresponding to datasets 3-10 of *teph4mob_1*.
These partial DDBs contain the dynamical matrices for the
8 $\qq$-points in the IBZ, as well as the dielectric tensor and the Born effective charges. 
Name the new DDB file *teph4mob_2_DDB*.

File *\$ABI_TUTORESPFN/Input/teph4mob_2.in* is an example of input file for *mrgddb*.

{% dialog tests/tutorespfn/Input/teph4mob_2.in %}

You can copy it in the *Work_eph4mob* directory, and run *mrgddb* using:

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

    Alternatively, one can the command line. 

    ```sh
    mrgdv merge teph4mob_3_DVDB teph4mob_1o_DS*_POT*
    ```

    Use **mrgdv --help** to access the documentation.

We now have all the phonon-related files needed to compute the mobility.
The DDB will be used to Fourier interpolate the phonon frequencies on an **arbitrarily** dense $\qq$-mesh while
the DVDB will be used to Fourier interpolate the DFPT scattering potentials [[cite:Brunin2020]].
The only ingredient left is the set of GS wavefunctions on the dense $\kk$-mesh.

## Calculation of the dense WFK file

In order to compute transport properties, we need to solve the Boltzmann
Transport Equation (BTE) in the SERTA on a relatively dense $\kk$-mesh to have enough wavevectors
inside the electron (hole) pockets.
<!--
You will compute the electron lifetimes and group velocities on this dense mesh. 
-->
We will therefore need the wavefunctions on this dense mesh. 
Note that in this tutorial, a single dense $\kk$-mesh will be used. 
However, the value for the mobility strongly depends on this
mesh and a convergence study should be performed by increasing the $\kk$-mesh density,
as well as the $\qq$-mesh used for the integration of the self-energy whose imaginary part
defines the electron lifetime.
This study is explained later and left to the user
but it should be clear even at this point that systems with small effective masses (e.g GaAs) will 
require denser $\kk$-meshes and will be more difficult to converge.

The computation of the dense WFK file is similar to a NSCF band structure computation.
The main difference is that we need wavefunctions on a $\kk$-mesh instead of a $\kk$-path
because these wavevectors will be used to evaluate integrals in the BZ.

The file *\$ABI_TUTORESPFN/Input/teph4mob_4.in* is an example of such computation.

{% dialog tests/tutorespfn/Input/teph4mob_4.in %}

It consists of two parts: the first one (dataset 1) computes the GS wavefunctions,
and the second one (datasets 2-3) computes the dense WFK that will be used to evaluate the mobility.
We also compute a denser WFK file that will be used with the double-grid method explained later.

We want to compute the mobility of electrons in the conduction band, therefore
we need to consider conduction bands in the computation of the WFK ([[nband]] = 8).

Copy the file in the *Work_eph4mob* directory, and run ABINIT:

```sh
abinit teph4mob_4.in > teph4mob_4.log 2> err &
```

*Note: do not forget to add the [[pseudos]] variable to this input file !*

## Calculation of the mobility

The computation of the mobility requires different convergence studies. 
We will explain them and their need in the following. 
Let us first explain the different parameters in a standard mobility computation.
The file *\$ABI_TUTORESPFN/Input/teph4mob_5.in* is an example of such computation.

{% dialog tests/tutorespfn/Input/teph4mob_5.in %}

We will need the WFK, the DDB and the DVDB files obtained previously.
One can specify the paths to these files using strings instead of integers:

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

Let's discuss the meaning of the e-ph variables in more details.

[[optdriver]] 7 is required to activate the e-ph driver while
[[eph_task]] -4 tells ABINIT that we only need the imaginary part 
of the e-ph self-energy, which directly gives the electron lifetimes, see [[cite:Brunin2020]].

The dense $\kk$-mesh is specified by [[ngkpt]] = 24 24 24. 
It should be the same mesh as the one used for the dense WFK computation. 
[[occopt]] 3 is required to correctly compute the
location of the Fermi level, using the Fermi-Dirac occupation function.
The initial grid used for the phonons (used in the DFPT computation) was a 4x4x4. 
These are the values that we need to specify for [[ddb_ngqpt]].
The dense $\qq$-mesh onto which the scattering potentials are interpolated and the 
e-ph matrix elements are computed is given by [[eph_ngqpt_fine]]. 
Note that [[eph_ngqpt_fine]] and [[ngkpt]] should be commensurated.

In this tutorial, we will use the same dense $\kk$- and $\qq$-meshes.
As a rule of thumb, a $\qq$-mesh twice as dense in each direction as the $\kk$-mesh, 
is needed to achieve fast convergence of the integrals.
Possible exceptions are systems with very small effective masses (e.g. GaAs) in which 
a very dense $\kk$-sampling is needed to to sample the electron (hole) pocket.
In this case, using the same sampling for electrons and phonons may be enough.

We will use the tetrahedron integration method to obtain the lifetimes (integration
over the $\qq$-mesh). This allows to efficiently filter out the $\qq$-points that do not contribute
to the lifetimes. Indeed, only a small fraction of the $\qq$-points belonging to the $\qq$-mesh
ensure energy and momentum conservation for a given $\kk$-point. 
All the other $\qq$-points do not need to be considered and can be filtered out.
The use of the tetrahedron method is automatically activated when only 
the imaginary part is wanted. 
It is possible to change this behaviour by using [[eph_intmeth]] albeit not recommended.

The list of temperatures for which the mobility is computed is specified by [[tmesh]].
The carrier concentration is deduced from the number of extra electrons in the unit cell,
specified by [[eph_extrael]]. 
To obtain results that are representative of the intrinsic mobility, 
we suggest to use a very small number, for instance $10^{15}$ to $10^{18}$ electrons per cm$^3$. 

The [[sigma_erange]] variable describes the energy window, below the VBM and above the
CBM where the lifetimes will be computed.
Since the BTE contains a derivative of the Fermi-Dirac occupation function centered on the Fermi level,
it is possible to filter the $\kk$-points that will contribute to the mobility and compute
the lifetimes for these $\kk$-points only. Indeed, the value of the derivative decreases rapidly
as we go further from the Fermi level. Only the states close to the band edges contribute.
This additional filtering technique allows one to compute only a few percents of all the lifetimes.
This variable should be subject to a convergence study, as explained in the next section.

Different tricks have been implemented to accelerate the computation.
The use of single precision in the FFT routines allows one to decrease the computational cost
without losing precision. This trick is activated by setting [[mixprec]] = 1.
Another trick to decrease the memory requirement is to decrease [[boxcutmin]] 
to a value smaller than 2 e.g. 1.5 or the more aggressive 1.1.
An exact representation of densities/potentials in $\GG$-space is obtained with [[boxcutmin]] = 2, 
but we found that using a value of 1.1 does not change the result 
but allows one to decrease the cost of the calculation and the memory by a factor ~8.
These tricks **are not activated by default** because users are supposed to perform preliminary tests
to make sure the quality of the results is not affected by these options.

If performance is really of concern, you can also try to set [[eph_mrta]] to 0.
By default, the code computes transport lifetimes both with the SERTA and the MRTA.
The MRTA requires the computation of the group velocities at $\kk$ and $\kk+\qq$. 
This part is relatively fast yet it does not come for free.
If you know in advance that you don't need MRTA results, it is possible to gain some speedup by disabling this part.

We can now examine the log file in detail.
After the standard output of the input variables,
the code reports the different parameters for the long-range potentials: the Born effective charges,
the dielectric constant, and the quadrupolar tensor. 
Make sure to have all of them in order to have an
accurate interpolation of the scattering potentials, see discussion in [[cite:Brunin2020]].

*Note: for the moment, the computation of the quadrupole is not available in the public version, 
and you should have the following in the log file:*

```sh
Have dielectric tensor: yes
Have Born effective charges: yes
Have quadrupoles: no
```

The code then outputs different information. The first of them is the location of the Fermi level
that will be used to compute the lifetimes. You can check that it is far enough from the band
edges so that the computed mobility will be intrinsic.

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

Over the initial 413 k-points in the IBZ, only 3 will be computed!
ABINIT then reads the WFK file and interpolates the potentials to compute the e-ph matrix elements.
The use of the tetrahedron method allows for the filtering of the q-points:

```sh
qpoints_oracle: calculation of tau_nk will need: 39 q-points in the IBZ. (nqibz_eff / nqibz):   9.4 [%]
```

Again, this leads to a significant speed-up of the computation.
Once this is done, the code starts looping over the 3 $\kk$-points for which the lifetimes are needed.

```sh
Computing self-energy matrix elements for k-point: [ 4.5833E-01,  4.5833E-01,  0.0000E+00] [ 1 / 3 ]
```

You can find various information for each $\kk$-point, such as the number of $\qq$-points 
belonging to the little group of $\kk$ (called IBZ(k) in the code), 
the number of $\qq$-point contributing to the imaginary part, the wall-time each step takes etc.
Finally, we have the results for the lifetimes (TAU) in the *teph4mob_5.out* file:

```sh
K-point: [ 4.5833E-01,  4.5833E-01,  0.0000E+00], T=    5.0 [K]
    B    eKS    SE2(eKS)  TAU(eKS)  DeKS
	5   3.573    0.000  31553.2    0.000
```

Only the first temperature is printed in the output file, but all the information can be found in the SIGEPH.nc file.

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
          64.00        0.23E+17        0.00E+00           38.37            0.00
         123.00        0.23E+17        0.00E+00          345.31            0.00
         182.00        0.23E+17        0.00E+00          423.32            0.00
         241.00        0.23E+17        0.00E+00          418.67            0.00
         300.00        0.23E+17        0.00E+00          363.11            0.00
```

The temperature is first given then the electron and hole densities followed by electron and hole mobilities. 
In this computation, we consider only electrons, so the values for holes are zero. 
Note that the transport driver is automatically executed after the e-ph calculation.
You can also run only the transport driver, provided you already have the lifetimes in a SIGEPH.nc file, by setting
[[eph_task]] = 7. 
This task can be performed only in serial and is very fast.

Now that you know how to obtain the mobility in a semiconductor for given k- and q-meshes, 
we can give more details about convergence and additional tricks that can be used 
to decrease the computational cost.

### Convergence w.r.t. the energy range

The first convergence study to be performed is to determine the energy range to be used for our computations. 
We can do that by performing the same mobility computation (same k- and q-meshes)
for an increasing energy window. 

The file *$\$ABI_TUTORESPFN/Input/teph4mob_6.in* is an example of such computation.

{% dialog tests/tutorespfn/Input/teph4mob_6.in %}

Copy the input file in the *Work_eph4mob* directory, and run ABINIT:

```sh
abinit teph4mob_6.in > teph4mob_6.log 2> err &
```

This run should take a few minutes.

We can now analyze the variation of the mobility with respect to [[sigma_erange]]. 
Once enough $\kk$-points are included, the mobility should not vary, and computing more $\kk$-points will
only increase the computation time.
Note that you should perform this convergence study with a $\kk$-mesh that is already dense enough to 
capture the band dispersion correctly. 
In this case, we are using a 24x24x24 mesh, which is not very dense for such computations. 
This means that, when increasing [[sigma_erange]], sometimes
no additional $\kk$-points are considered. 
It is the case here for the first 3 datasets (3 k-points), and the last two datasets (6 k-points). 
If a finer mesh was used, the number of $\kk$-points would have increased in a more monotonic way.

### Convergence w.r.t. the k- and q-meshes

Once the energy window has been set, we can start to converge the mobility with respect to the 
dense $\kk$- and $\qq$-meshes. 
The previous computations used 24x24x24 meshes. This is far from convergence.
In silicon, for instance, a 45x45x45 $\kk$-mesh and 90x90x90 $\qq$-mesh are needed 
to have results converged within 5%. 

In order to compute the mobility with a $\qq$-mesh twice as dense as the $\kk$-mesh, there are two possibilities. 
Let us take the previous example of silicon.

  * Run a computation with [[ngkpt]] = 90 90 90, [[eph_ngqpt_fine]] = 90 90 90, and [[sigma_ngkpt]] = 45 45 45.
Using [[sigma_ngkpt]] will select the $\kk$-points belonging to the 45x45x45 mesh, but each lifetime will be computed
with a 90x90x90 q-mesh.

  * Run a computation with [[ngkpt]] = 90 90 90, [[eph_ngqpt_fine]] = 90 90 90 and [[sigma_ngkpt]] = 90 90 90.
In this way, you have the mobility with 90x90x90 $\kk$- and $\qq$-meshes. You can then run again the transport driver only,
by setting [[eph_task]] = 7, and setting [[transport_ngkpt]] = 45 45 45. This will downsample the $\kk$-mesh used
in the computation of the mobility. This second option has the advantage that it delivers two mobilities 
(useful for convergence studies), but it requires more computational time.

You can run again the previous input files by densifying the different meshes.
For the densest grids, you might need to run with multiple MPI processes.
You should obtain something like this, for T = 300 K:

![](eph4mob_assets/teph4mob_conv.png)

Note that in order to get sensible results, one should use a denser DFPT $\qq$-mesh (around 9x9x9), 
and a larger energy cutoff for the planewave basis set ([[ecut]]). 
The inputs of this tutorial have been tuned so that the computations are quite fast,
but they are quite far from convergence. 
In real life, you should perform convergence studies to find suitable parameters.

### Double-grid technique

Another possibility to improve the results without increasing the computation time significantly is to use
the double-grid technique.
In this case, a coarse mesh is used for the $\kk$-mesh and the $\qq$-mesh for the e-ph matrix elements,
but a finer mesh is used for the phonon absorption-emission terms. 
This technique allows one to better capture these processes, while computing the matrix elements on a coarser mesh. 
The efficiency of this method depends
on the polar divergence of the matrix elements: if this divergence is very difficult to capture numerically,
the coarse $\qq$-mesh for the e-ph matrix elements will have to be dense.

The double-grid technique requires a second WFK file, containing the KS eigenvalues on the dense mesh. 
You can specify the path to the dense WFK file using [[getwfkfine_filepath]]:

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

	 coarse:                24          24          24
	 dense:                 48          48          48

The mobility obtained, at 300 K, is 158.01 cm$^2$/V/s. 
Using a 48x48x48 $\qq$-mesh for the matrix elements as well would give 96.09. 
The result is indeed improved, since using a 24x24x24 mesh for everything gives 363.11. 
You can also use a denser fine mesh, but always a multiple of the initial coarse mesh 
(in this case, 72x72x72, 96x96x96, etc). 
However, we found that there is very little use to go beyond a mesh three times as dense as the coarse one.
Using a 72x72x72 fine mesh for the energies gives a mobility of 149.87 cm$^2$/V/s,
and a 96x96x96 mesh leads to 146.24 cm$^2$/V/s: the improvement is indeed rather limited.

### On the importance of the initial DFPT mesh

At this point it is worth commenting about the importance of the initial DFPT $\qq$-mesh.
The Fourier interpolation implicitly assumes that the signal in $\RR$-space decays quickly hence
the quality of the *interpolated* phonon frequencies and of the *interpolated* DFPT potentials, 
between the ab-initio points depends on the spacing of the initial $\qq$-mesh that 
in turns defines the size of the Born-von-Karman supercell.

A more detailed discussion can be found in [[cite:Brunin2020]], [[cite:Verdi2015]] and [[cite:Sjakste2015]].

## Additional tricks

TODO: Discuss [[dipdip]], [[asr]], [[chneut]]

### MPI parallelism and memory requirements

There are five different MPI levels that can be used to distribute the workload 
and the most memory-demanding data structures.
By default, the code tries to reach some compromise between memory requirements and time to solution. 
by activating the parallelism over $\qq$-points if no other input is provided from the user.
You can however specify manually the MPI distribution across the five different levels
by using [[eph_np_pqbks]], a list of 5 integers. 
The product of these five numbers **must be equal** to the total number of MPI processes.
The first number gives the number of processes for the parallelization over perturbations. 
The allowed value range between 1 and 3 x [[natom]], and should be a divisor 
of 3 x [[natom]] to distribute the work equally.
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

TODO: Discuss [[dvdb_qcache_mb]]

### In-place restart

All the results of the calculation are stored in a single SIGEPH.nc file
for all the $\kk$-points (and spins) considered.
The list of $\kk$-points is initialized at the beginning of the run and an internal table 
stores the status of the k-point (computed or not).
This means that calculations that get killed by the resource manager due to time limit can use 
the netcdf file to perform an automatic in-place restart.
Just set [[eph_restart]] to 1 in the input file and rerun the job.

### Transport calculation from SIGEPH.nc 

TODO: [[transport_ngkpt]]

### How to compute only the k-points close to the CBM/VBM

TODO: [[getkerange_filepath]]
