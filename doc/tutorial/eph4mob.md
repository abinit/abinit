---
authors: GB, MG
---

# Electron-phonon tutorial 2

## Phonon-limited mobility in AlAs

This tutorial aims at showing how to get the phonon-limited (intrinsic) mobility in semiconductors.
Such a computation implies different steps:

  * ground state computation
  * DFPT calculation on a homogeneous q-grid
  * merging the DDB and POT files
  * getting the wavefunction on a dense k-mesh
  * electron-phonon interpolation of the scattering potentials to compute the electron lifetimes
  * computation of transport properties

They can be summarized by the following:

![](eph4mob_assets/workflow.png)

Anaddb will not be used in this tutorial. All the capabilities are integrated directly in ABINIT.
Note that in order to analyze the results in more details, you will need to have ABINIT compiled with netCDF.
Only the final results are printed in the log and output files.

[TUTORIAL_README]

This tutorial should take about 1 hour. It is strongly recommended to read [[cite:Brunin2020]]
before running the tutorial in order to understand it correctly.

Visualisation tools are NOT covered in this tutorial.
Powerful visualisation procedures have been developed in the Abipy context,
relying on matplotlib. See the README of [Abipy](https://github.com/abinit/abipy)
and the [Abipy tutorials](https://github.com/abinit/abitutorials).

## 1 Calculation of the ground state and phonon structure of fcc AlAs

*Before beginning, you might consider making a different subdirectory to work
in. Why not create Work_eph4mob ?*

It is presumed that the user has already followed the Tutorials [RF1](rf1) and [RF2](rf2),
and understands the calculation of ground state and response (phonons using density-functional
perturbation theory (DFPT)) properties with ABINIT.

The file *teph4mob_1.files* lists the file names and root names for the first run (GS+perturbations).
You can copy it to the working directory. You can also copy the file *teph4mob_1.in* to your working directory.

```sh
cd $ABI_TUTORESPFN/Input
mkdir Work_eph4mob
cd Work_eph4mob
cp ../teph4mob_1.file .
cp ../teph4mob_1.in .
```

{% dialog tests/tutorespfn/Input/teph4mob_1.in %}

This tutorial starts by the generation of a database, that might be quite time-consuming.
You can immediately start this run. 

    abinit < teph4mob_1.files > teph4mob_1.log 2> err &

The calculation is done for AlAs, the same crystalline material as for the first two tutorials on DFPT.
Many input parameters are also quite similar. For more details about this first step, please refer to
the first and second tutorials on DFPT.

## 2 Merging the derivative databases and potentials

Use MRGDDB to create the merged DDB from the eight DDB's corresponding to
datasets 3 to 10 of the teph4mob_1 job, containing the dynamical matrices for the
8 q points, as well as the response to the electric field (dielectric tensor
and Born effective charges). Name the new DDB *teph4mob_2_DDB*.
File *\$ABI_TUTORESPFN/Input/teph4mob_2.in* is an example of input file for MRGDDB.

{% dialog tests/tutorespfn/Input/teph4mob_2.in %}

You can copy it in the *Work_eph4mob* directory, and run the merge as follows:

    mrgddb < teph4mob_2.in

Use MRGDV to create the merged potentials from the 29 POT files corresponding to
datasets 3 to 10 of the teph4mob_1 job. Name the new file *teph4mob_3_DVDB*.
File *\$ABI_TUTORESPFN/Input/teph4mob_3.in* is an example of input file for MRGDV.

{% dialog tests/tutorespfn/Input/teph4mob_3.in %}

You can copy it in the *Work_eph4mob* directory, and run the merge as follows:

    mrgdv < teph4mob_3.in

You now have the phonon-related files that are needed for the computation of the mobility.
The DDB will be used to interpolate the phonon frequencies on the dense q-mesh, and
the DVDB file will be used to Fourier interpolate the scattering potentials and 
compute the e-ph matrix elements that we will need to compute the lifetimes. The only 
ingredient left is the wavefunction on the dense mesh, which is the next step.

## 3 Calculation of the dense WFK file

In order to compute transport properties, you need to solve the Boltzmann
Transport Equation (BTE) on a relatively dense k-mesh. You will compute the
electron lifetimes and velocities on this dense mesh. You will therefore need
the wavefunction on this dense mesh. Note that in this tutorial, a single dense
mesh will be used. However, the value for the mobility strongly depends on this
mesh and a convergence study should be performed by increasing the k-mesh density,
as well as the q-mesh used to obtain the lifetimes.
This study is explained later and left to the user.


The computation of the dense WFK file is similar to a band structure computation.
The file *$\$ABI_TUTORESPFN/Input/teph4mob_4.in* is an example of such computation.

{% dialog tests/tutorespfn/Input/teph4mob_4.in %}

It consists of two parts: the first one (dataset 1) computes the ground state wavefunction,
and the second one (datasets 2-3) computes the dense WFK that will be used to compute the transport.
We also compute a denser WFK file that will be used with the double-grid method explained later.
We want to compute the mobility of electrons in the conduction band, therefore
we need to consider conduction bands in the computation of the WFK ([[nband]] = 8).
You can copy it in the *Work_eph4mob* directory, and run ABINIT:

    abinit < teph4mob_4.files > teph4mob_4.log 2> err &

## 4 Calculation of the mobility

The computation of the mobility requires different convergence studies. We will
explain them and their need in the following. 
Let us first explain the different parameters in a simple mobility computation.
The file *$\$ABI_TUTORESPFN/Input/teph4mob_5.in* is an example of such computation.

{% dialog tests/tutorespfn/Input/teph4mob_5.in %}

You will need the wavefunction, DDB and DVDB files obtained previously.
You can give the paths to these file using:

    getwfk_path "teph4mob_4o_DS2_WFK"
    getddb_path "teph4mob_2_DDB"
    getdvdb_path "teph4mob_3_DVDB"

You can copy the input file in the *Work_eph4mob* directory, and run ABINIT:

    abinit < teph4mob_5.files > teph4mob_5.log 2> err &

This run should take $\sim$15 seconds on a recent CPU. 

[[optdriver]] = 7 is required to specify that we will run an e-ph computation.
[[eph_task]] = -4 tells ABINIT that we need only the imaginary part of the Fan-Migdal
self-energy, which directly gives the electron lifetimes, see [[cite:Brunin2020]].

The dense k-mesh is specified by [[ngkpt]] = 24 24 24. It should be the same mesh as the
one used for the dense WFK computation. [[occopt]] = 3 is required to correctly compute the
location of the Fermi level, using the Fermi-Dirac occupation function.
The initial grid used for the phonon (the DFPT computation) was a 4x4x4. These are the
values that we need to specify for [[ddb_ngqpt]].
The dense q-mesh onto which the scattering potentials are interpolated and the 
e-ph matrix elements are computed is given by
[[eph_ngqpt_fine]]. [[eph_ngqpt_fine]] and [[ngkpt]] should be commensurate.
We strongly advice to use a q-mesh twice as dense in each direction as the k-mesh.
In this tutorial, we will use the same dense k- and q-meshes.

Different tricks have been developed to accelerate the computation.
The use of single precision in some sums allows to decrease the memory
without losing precision. This trick is activated by setting [[mixprec]] = 1.
Another trick to decrease the memory requirement is to decrease [[boxcutmin]].
An exact representation is obtained with [[boxcutmin]] = 2, but we found that
using a value of 1.1 does not change the result (in the materials we tested; 
this might not be true in general and requires testing from the user) but allows
to decrease the memory by a factor close to 8. 

We will use the tetrahedron integration method to obtain the lifetimes (integration
over the q-mesh). This allows to efficiently filter out the q-points that do not contribute
to the lifetimes. Indeed, only a small fraction of the q-points belonging to the q-mesh
ensure energy and momentum conservation for a given k-point. All the other q-points do not
need to be considered and can be filtered out.
The use of the tetrahedron
method is automatically activated (see [[eph_intmeth]]).  

The temperatures at which the mobility is computed are given by [[tmesh]].
The carrier concentration is deduced from the number of extra electrons in the unit cell,
given by [[eph_extrael]]. You should use a very small number, corresponding to
$10^{15}$ to $10^{18}$ electrons per cm$^3$. 
The [[sigma_erange]] variable describes the energy window, below the VBM and above the
CBM, within which the lifetimes will be computed.
Indeed, since the BTE contains a derivative of the Fermi-Dirac occupation function centered on the Fermi level,
it is possible to filter the k-points that will contribute to the mobility and compute
the lifetimes for these k-points only. Indeed, the value of the derivative decreases rapidly
as we go further from the Fermi level. Only the states close to the band edges contribute.
This additional filtering method
allows one to compute only a few percents of all the lifetimes.
This variable should be subject to a 
convergence study, as explained in the next section.

You can now examine the log file. Let us describe it in detail here.
After the standard printing of the input variables,
the code reports the different parameters for the long-range potentials: the Born effective charges,
the dielectric constant, and the quadrupolar tensor. Make sure to have all of them in order to have a 
precise interpolation of the scattering potentials [[cite:Brunin2020]] !

*Note: for the moment, the computation of the quadrupole is not openly available, and you should have
the following in the log file:*

	 Have dielectric tensor: yes
	 Have Born effective charges: yes
 	 Have quadrupoles: no

The code then gives different information. The first of them is the location of the Fermi level
that will be used to compute the lifetimes. You can check that it is far enough from the band
edges so that the computed mobility will be intrinsic.

	 Valence Maximum:   2.3564 (eV) at: [ 0.0000E+00,  0.0000E+00,  0.0000E+00]
	 Conduction minimum:   3.5254 (eV) at: [ 5.0000E-01,  5.0000E-01,  0.0000E+00]
	 Fermi level:  3.0878 (eV)

ABINIT also finds the list of k-point belonging to the dense mesh that have states within the
energy window:

	 Found 3 k-points within erange:  0.000  0.150  (eV)
	 min(nbcalc_ks): 1 MAX(nbcalc_ks): 1

Over the 413 k-points, only 3 will be computed !
ABINIT then reads the WFK file and interpolates the potentials to compute the matrix elements.
The use of the tetrahedron method allows for the filtering of the q-points:

	 qpoints_oracle: calculation of tau_nk will need: 39 q-points in the IBZ. (nqibz_eff / nqibz):   9.4 [%]

Again, this leads to a speed-up of the computation.
Once this is done, the code starts looping over the 3 k-points for which the lifetimes are needed.

	 Computing self-energy matrix elements for k-point: [ 4.5833E-01,  4.5833E-01,  0.0000E+00] [ 1 / 3 ]

You can find various information for each k-point, such as the number of q-points belonging to the little group of k
(IBZ(k)), the number of q-point contributing, the time every step takes,... You can find the value of the lifetimes in the 
teph4mob_5.out file:

	K-point: [ 4.5833E-01,  4.5833E-01,  0.0000E+00], T=    5.0 [K]
	   B    eKS    SE2(eKS)  TAU(eKS)  DeKS
	   5   3.573    0.000  34495.8    0.000

Only the first temperature is printed in the output file, but all the information can be found in the SIGEPH.nc file.

!!! tip

    With |AbiPy|, one can easily have access to all the data of the computation. For instance, one can plot the
    electron linewidths:

        abiopen.py teph4mob_5o_SIGEPH.nc --expose

![](eph4mob_assets/linewidths.png)

At the end of the .out and .log files, the mobility is printed:

	Transport calculation results
	 Temperature [K]             e/h density [cm^-3]          e/h mobility [cm^2/Vs]
	            5.00        0.23E+17        0.00E+00            0.00            0.00
	           64.00        0.23E+17        0.00E+00           42.80            0.00
	          123.00        0.23E+17        0.00E+00          363.44            0.00
	          182.00        0.23E+17        0.00E+00          433.51            0.00
	          241.00        0.23E+17        0.00E+00          425.82            0.00
	          300.00        0.23E+17        0.00E+00          368.46            0.00

The temperature is first given then the electron and hole densities, and electron and hole mobilities. 
In this computation, we consider only the electrons, so the values for the holes are 0. 
Note that the transport driver is automatically run after the e-ph driver.
You can also run only the transport driver, if you have the lifetimes in a SIGEPH.nc file, by setting
[[eph_task]] = 7. This task can be performed only in serial and is very fast.


Now that you know how to obtain the mobility in a semiconductor for given k- and q-meshes, we can
give more details about convergence and tricks that can be used to decrease the computational cost.

### Convergence w.r.t. the energy range

The first convergence study to be performed is to determine the energy range that we will use
for our computations. We can do that by performing the same mobility computation (same k- and q-meshes)
for an increasing energy window. 
The file *$\$ABI_TUTORESPFN/Input/teph4mob_6.in* is an example of such computation.

{% dialog tests/tutorespfn/Input/teph4mob_6.in %}

You can copy the input file in the *Work_eph4mob* directory, and run ABINIT:

    abinit < teph4mob_6.files > teph4mob_6.log 2> err &

This run should take a few minutes.
You can then analyze the variation of the mobility with the value of [[sigma_erange]]. Once
enough k-points are included, the mobility should not vary, and computing more k-points will
only increase the computation time.
Note that you should perform this convergence study with a k-mesh already dense enough to 
capture the band structure correctly. In this case, we use a 24x24x24 mesh, which is not
very dense for such computations. This means that, when increasing [[sigma_erange]], sometimes
no additional k-points are considered. It is the case here for the first 3 datasets (3 k-points),
and the last two datasets (6 k-points). If a finer mesh was user, the number of k-points would have
increased in a more monotonic way.

### Convergence w.r.t. the k- and q-meshes

Once the energy window has been set, you should converge the mobility with respect to the 
dense k- and q-meshes. The previous computations used 24x24x24 meshes. This is far from convergence.
For instance, in silicon, a 45x45x45 k-mesh and 90x90x90 q-mesh are needed to have a result converged
within 5%. In order to compute the mobility with a q-mesh twice as dense as the k-mesh, there are two
possibilties. Lets us take the previous example of silicon.

  * Run a computation with [[ngkpt]] = 90 90 90, [[eph_ngqpt_fine]] = 90 90 90, and [[sigma_ngkpt]] = 45 45 45.
Using [[sigma_ngkpt]] will select the k-points belonging to the 45x45x45 mesh, but each lifetime will be computed
with a 90x90x90 q-mesh.
  * Run a computation with [[ngkpt]] = 90 90 90, [[eph_ngqpt_fine]] = 90 90 90 and [[sigma_ngkpt]] = 90 90 90.
In this way, you have the mobility with 90x90x90 k- and q-meshes. You can then run again the transport driver only,
by setting [[eph_task]] = 7, and setting [[transport_ngkpt]] = 45 45 45. This will downsample the k-mesh used
in the computation of the mobility. This second option has the advantage that it delivers two mobilities 
(useful for convergence studies), but it
requires more computational time.

You can run again the previous input files, by densifying the different meshes. 
You should obtain something like this, for T = 300 K:

![](eph4mob_assets/teph4mob_conv.png)

Note that in order to get truly sensible results, one should use a denser DFPT q-mesh (around 9x9x9), and a larger cut-off
for the planewave basis set ([[ecut]]). The inputs in this tutorial have been made so that the computations are quite fast,
but they are not converged. You should find the suitable parameters for your own work.

### Use of the double-grid technique

Another possibility to improve the results without increasing the computation time significantly is to use
a double-grid technique. A coarse mesh will be used for the k-mesh and the q-mesh for the e-ph matrix elements,
but a finer mesh will be used for the phonon absorption-emission terms. This technique helps to better capture
these processes, while computing the matrix elements on a coarser mesh. The efficiency of this method will depend
on the polar divergence of the matrix elements: if this divergence is very difficult to capture numerically,
the coarse q-mesh for the e-ph matrix elements will have to be dense.

The double-grid technique requires a second WFK file, containing the dense mesh. You can specify the path to the dense
WFK file using [[getwfkfine_path]]:

     getwfkfine_path "teph4mob_4o_DS3_WFK"

The file *$\$ABI_TUTORESPFN/Input/teph4mob_7.in* is an example of such computation.

{% dialog tests/tutorespfn/Input/teph4mob_7.in %}

You can copy the input file in the *Work_eph4mob* directory, and run ABINIT:

    abinit < teph4mob_7.files > teph4mob_7.log 2> err &

In the log file, you will now find information about the double-grid method:

	 coarse:                24          24          24
	 dense:                 48          48          48

The mobility obtained, at 300 K, is 160.34 cm$^2$/V/s. Using a 48x48x48 q-mesh for the matrix elements as well would give
97.10. The result is indeed improved, since using a 24x24x24 mesh for everything gives 368.46. 
You can also use a denser fine mesh, but always a multiple of the initial coarse mesh (in this case,
72x72x72, 96x96x96, etc). However, we found that there is very little use to go beyond a mesh three times
as dense as the coarse mesh.

### Note on the different levels of MPI parallelism

There are different MPI levels that are used to distribute the workload and the memory-demanding data structures.
By default, the code uses a compromise between memory and time consumption. Depending on your needs, you can however
tune it using the variabie [[eph_np_pqbks]], given by an array of 5 values. The product of the five numbers that you
give has to be equal to the number of MPI processes of your computation.
The first number concerns the number of processes for the parallelization over perturbation. This ranges between 1 and 3 x [[natom]],
and should be a divisor of 3 x [[natom]] in order to be efficient. The higher this number, the lower the memory requirements.
The second number determines the parallelization over the q-points. This parallelization level is efficient to decrease the
computational time. 
The parallelization over bands is usually not relevant in mobility computations, since only states close to the VBM or CBM are
considered. It is however useful when the real-part of the self-energy is needed.
The MPI parallelism over k-points and spins is efficient at the level of the wall-time
but it requires HDF5 + MPI-IO support and memory does not scale. Use these additional levels if the memory requirements
are under control and you need to boost the calculation. Note also that in this case the output results are written to
different text files, only the SIGEPH.nc file will contains all the k-points and spins.

