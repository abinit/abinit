---
authors: MG
---

# First tutorial on GWR (GW in real-space and imaginary time)

## The quasi-particle band structure of Silicon in the one-shot GW approximation.

This tutorial aims at showing how to calculate self-energy corrections to the
DFT Kohn-Sham (KS) eigenvalues in the one-shot GW approximation using the GWR code

The user should be familiarized with the four basic tutorials of ABINIT,
see the [tutorial home page](/tutorial),
and is strongly encouraged to read the [introduction to the GWR code](/tutorial/gwr_intro)
before running these examples.

This tutorial should take about 1.5 hours.

[TUTORIAL_README]

### Ground-state and band structure

*Before beginning, you might consider creating a different subdirectory to work in.
Why not create Work_gwr?*

The file *tgwr_1.abi* is the input file for the first step:
a SCF run followed by a KS NSCF band structure calculation along a high-symmetry $\kk$-path.
Copy it to the working directory with:

```
mkdir Work_gwr
cd Work_gwr
cp ../tgwr_1.abi .
```

You may want to immediately start the job in background with:

```
mpirun -n 1 tgwr_1.abi > tgwr_1.log 2> err &
```

so that we have some time to discuss the input while ABINIT is running.

{% dialog tests/tutorial/Input/tgwr_1.abi %}

!!! tip

    All the input files of this tutorial can be executed in parallel by just
    increasing the value of n in `mpirun -n`.
    Clearly the size of the problem is very small so do not expect
    great parallel performance if you start to use dozens of MPI processes.

The first dataset produces the KS density file that is then used to compute the band structure in the second dataset.
Since all the input files of this tutorial use the same crystalline structure, pseudos
and cutoff energy [[ecut]] for the wavefunctions,
we declare these variables in an external file that will be **included** in all the other input files
using the syntax:

```
# Include geometry and pseudos
include "gwr_include.abi"
```

If you open the include file:

{% dialog tests/tutorial/Input/gwr_include.abi %}

you will notice that we are using norm-conserving (NC) pseudos taken from
the standard scalar-relativistic table of the [PseudoDojo](https://www.pseudo-dojo.org/).

!!! important

    For GW calculations, we strongly recommend using NC pseudos from the **stringent table**
    as these pseudos have closed-shells treated as valence states.
    This is important for a correct description of the matrix elements of the exchange part of the self-energy
    that are rather sensitive to the overlap between the wavefunctions [[cite:Setten2018]].
    In the special case of silicon, there is no difference between the standard version and the stringent version of the pseudo,
    but if we consider e.g. Ga, you will notice that the stringent version includes all the 3spd states in valence
    besides the outermost 4sp electrons whereas the standard pseudo for Ga designed for GS calculations includes only the 3d states.
    Since we are using scalar-relativistic (SR) pseudos, it is not possible to treat spin-orbit coupling (SOC)
    and [[nspinor]] must be set to 1 (default value).
    To include SOC with NC, one should use fully-relativistic (FR) pseudos, use [[nspinor]] 2 in all the input files
    and then set the appropriate value of [[nspden]] to 4 or 1 depending on whether the system is magnetic or not.

    For efficiency reasons, these examples use underconverged parameters:
    the cutoff energy [[ecut]] is set to XXX that is slightly smaller
    than the recommended value reported in the PseudoDojo table

<!--
One final comment on MPI parallelism:
In the input, we set [[paral_kgb]] to 1 to enable k-point, band, and FFT parallelization,
and [[autoparal]] to 1 to allow ABINIT to determine an optimal distribution, with [[npfft]] fixed at 1.
-->

At this point, the calculation should have completed, and we can have a look
at the KS band structure to locate the position of the band edges.

!!! tip

    If |AbiPy| is installed on your machine, you can use the |abiopen| script
    with the `--expose` option to visualize the band structure from the GSR.nc file:

        abiopen.py tgwr_1o_DS2_GSR.nc --expose

    ![](base1_assets/abiopen_tgwr1_1.png)

    To print the results to terminal, use:

        abiopen.py tgwr_1o_DS2_GSR.nc -p


```
============================== Electronic Bands ==============================
Number of electrons: 8.0, Fermi level: 4.369 (eV)
nsppol: 1, nkpt: 101, mband: 12, nspinor: 1, nspden: 1
smearing scheme: none (occopt 1), tsmear_eV: 0.272, tsmear Kelvin: 3157.7
Direct gap:
    Energy: 2.564 (eV)
    Initial state: spin: 0, kpt: $\Gamma$ [+0.000, +0.000, +0.000], band: 3, eig: 4.369, occ: 2.000
    Final state:   spin: 0, kpt: $\Gamma$ [+0.000, +0.000, +0.000], band: 4, eig: 6.933, occ: 0.000
Fundamental gap:
    Energy: 0.590 (eV)
    Initial state: spin: 0, kpt: $\Gamma$ [+0.000, +0.000, +0.000], band: 3, eig: 4.369, occ: 2.000
    Final state:   spin: 0, kpt: [+0.429, +0.000, +0.429], band: 4, eig: 4.959, occ: 0.000
Bandwidth: 11.976 (eV)
Valence maximum located at kpt index 0:
    spin: 0, kpt: $\Gamma$ [+0.000, +0.000, +0.000], band: 3, eig: 4.369, occ: 2.000
Conduction minimum located at kpt index 12:
    spin: 0, kpt: [+0.429, +0.000, +0.429], band: 4, eig: 4.959, occ: 0.000
```

Silicon is an indirect band gap semiconductor:
the CBM is located at the $\Gamma$ point while the VMB is located at ~[+0.429, +0.000, +0.429],
At the PBE level, the direct gap is 2.564 while the fundamental band gap is ~0.590.
Both values strongly underestimate the experimental results that are ~3.4 eV and ~1.12 eV, respectively.

!!! important

    Similarly to the conventional GW code, also GWR can compute QP corrections only
    for the $\kk$-points belonging to the $\kk$-mesh used for generating the WFK file.
    Before running GW calculations is always a good idea to analyze carefully the KS band
    structure to understand the location of the band edges and then select
    the most appropriate $\kk$-mesh. for the GW computation.

    Also, we will use the GSR.nc file with the KS energies along the path
    to perform an interpolation of the GWR results in the last part of this tutorial.


### Generation of the WFK file with empty states

In the second input file, we use the density file computed previously (tgwr\_1o\_DS1\_DEN)
to generate three WFK files with three different $\Gamma$-centered $\kk$-meshes ([[ngkpt]] = 2x2x2 and 4x4x4, 6x6x6).

```
ndtset 3

ngkpt1   2 2 2
ngkpt2   4 4 4
ngkpt3   6 6 6
```

This allows us to perform convergence studies with respect to the BZ sampling in the last part of this tutorial.
Let us recall that shifted $\kk$-meshes are not supported by GWR.
In another words, [[nshiftk]] must be set to 1 with [[shiftk]] = 0 0 0 when producing the WFK file.
The density file, on the contrary, can be generated with shifted $\kk$-meshes, as usual.

Let's immediately start the job in background with:

```
mpirun -n 1 tgwr_2.abi > tgwr_2.log 2> err &
```

{% dialog tests/tutorial/Input/tgwr_2.abi %}

Note that here we are using [[gwr_task]] = "HDIAGO" to perform a **direct diagonalization**
of the KS Hamiltonian constructed from the input DEN file.
This procedure differs from the one used in the other GW tutorials in which the WFK file
is generated by performing an **iterative diagonalization** in which only the application of the Hamiltonian is required.
The reason is that the direct diagonalization outperforms iterative methods when many empty states are required,
especially if one can take advantage of ScalaPack to distribute the KS Hamiltonian matrix.

Here, we request 250 bands,
Let's recall that in the previous GW tutorial, [[nband]] = 100 was considered converged within 30 meV
but with GWR we can afford more bands since [[nband]] enters into play only during the initial construction
of the Green's function

Clearly, when studying new systems, the converged [[nband]] value is not known beforehand.
Therefore, it is important to plan ahead and choose a reasonably large number of bands
to avoid regenerating the WFK file multiple times just to increase [[nband]].

Note that [[paral_kgb]] = 1 is only available in the ground-state (GS) part.
The GWR code employs its own distribution scheme, which depends on the value of [[gwr_task]].
When "HDIAGO" is used, the distribution is handled automatically at runtime, and the user has no control over it.
Please refer to the note below for guidance on choosing an appropriate number of MPI processes
to ensure an efficient workload distribution.

!!! important

    The direct diagonalization is MPI-parallelized across three different levels:
    collinear spin $\sigma$ (not used here), $\kk$-points in the IBZ and
    Scalapack distribution of the $H^\sigma_\kk(\bg,\bg')$ matrix.

    ABINIT will try to find an "optimal" distribution of the workload at runtime, yet there are a couple
    of things worth keeping in mind when choosing the number of MPI processes for this step.
    Ideally the total number of cores should be a multiple of [[nkpt]] * [[nsppol]] to avoid load imbalance.

    In order to compute **all** the eigenvectors of the KS Hamiltonian, one can use gwr_task "HDIAGO_FULL".
    In this case the value of [[nband]] is automatically set to the total number of plawewaves
    for that particular $\kk$-point.
    No stopping criterion such as [[tolwfr]] or number of iterations [[nstep]] are required when Scalapack is used.

    Again, multi-datasets are **strongly discouraged** if you care about performance.


### Our first GWR calculation

For our first GWR run, we use a minimalistic input file that
performs a GWR calculation using the DEN and the WFK file produced previously
and different minimax meshes.

First of all, you may want to start immediately the computation by issuing:

```
mpirun -n 1 tgwr_3.abi > tgwr_3.log 2> err &
```

with the following input file:

{% dialog tests/tutorial/Input/tgwr_3.abi %}

This input contains some variables whose meaning is the same as in the conventional GW code,
and other variables that are specific to the GWR code whose name starts with `gwr_`.

We use [[optdriver]] 6 to enter the GWR code and [[gwr_task]] activates a one-shot GW calculation.
To reduce the wall-time, we use a minimax mesh with [[gwr_ntau]] = 6 points.
This is the minimum number of points that can be used.
Most likely, six points are not sufficient, but the convergence study for [[gwr_ntau]] is postponed to the next sections.

[[getden_filepath]] specifies the density file used to compute $v_{xc}[n](\rr)$,
while [[getwfk_filepath]] specifies the WFK file with empty states used to build the Green's function.
Keep in mind that the $\kk$-mesh specified in the input via [[ngkpt]], [[nshiftk]] and [[shiftk]] must
agree with the one found in the WFK file else the code will abort.

[[gwr_boxcutmin]] is set to 1.1 in order to accelerate the calculation and reduce memory requirements.
This is another parameter that should be subject to carefully convergence studies as systems with localized electrons
such as 3d or 4f electrons may require larger values (denser FFT meshes).
Please take some time to read the variable description of [[gwr_boxcutmin]] before proceeding.

Now, let us turn our attention to the variables that are also used in the conventional GW code.
The cutoff of the polarizability and $W$ is defined by [[ecuteps]] as in the conventional GW code.
The cutoff for the exchange part of the self-energy is given by [[ecutsigx]].
Theoretically, we need [[ecutsigx]] = 4 [[ecut]] to have an exact treatment of $\Sigma_x$.
The computation of the exchange is relatively fast as only occupied states are involved.

The $\kk$-points and the band range for the QP corrections are set expliclty via [[nkptgw]], [[kptgw]] and [[bdgw]].
Alternatively, one can use [[gw_qprange]].

For the spectral function, we have the following variables [[nfreqsp]], [[freqspmax]]
TODO: Mention [[inclvkb]] == 2, [[symsigma]]

We can now have a look at the main output file:

{% dialog tests/tutorial/Refs/tgwr_3.abo %}

First of all, we have a section that summarizes the most important GWR parameters in a Yaml document:

```yaml
--- !GWR_params
iteration_state: {dtset: 1, }
gwr_task: G0W0
nband: 100
ntau: 6
ngkpt: [2, 2, 2, ]
ngqpt: [2, 2, 2, ]
chi_algo: supercell
sigma_algo: 'BZ-convolutions'
nkibz: 3
nqibz: 3
inclvkb: 2
q0: [  1.00000000E-05,   2.00000000E-05,   3.00000000E-05, ]
gw_icutcoul: 6
green_mpw: 412
tchi_mpw: 190
g_ngfft: [12, 12, 12, 12, 12, 12, ]
gwr_boxcutmin:   1.10000000E+00
P gwr_np_kgts: [1, 1, 1, 1, ]
P np_kibz: [1, 1, 1, ]
P np_qibz: [1, 1, 1, ]
min_transition_energy_eV:   2.20247744E-02
max_transition_energy_eV:   3.84181232E+00
eratio:   1.74431404E+02
ft_max_err_t2w_cos:   3.39269547E-02
ft_max_err_w2t_cos:   4.27492521E-03
ft_max_err_t2w_sin:   1.00974389E+00
cosft_duality_error:   4.08855664E-04
Minimax imaginary tau/omega mesh: !Tabular | # tau, weight(tau), omega, weight(omega)
    1 1.16987E-01   3.20191E-01   8.71836E-03   1.84660E-02
    2 8.20298E-01   1.25063E+00   3.30978E-02   3.36974E-02
    3 3.19454E+00   3.95048E+00   8.75746E-02   8.49211E-02
    4 1.01456E+01   1.10178E+01   2.40282E-01   2.55119E-01
    5 2.85210E+01   2.82127E+01   7.50121E-01   9.22624E-01
    6 7.50861E+01   7.27575E+01   2.97572E+00   4.74388E+00
...
```

Some of the entries in this dictionary have a direct correspondence with ABINIT variables and
we don't discuss them further.

`nkibz`: Number of $\kk$-points in the IBZ
`qkibz`: Number of $\qq$-points in the IBZ
`green_mpw`: maximum number of PWs for the Green's function.
`tchi_mpw`: maximum number of PWs for the polarizability.
`g_ngfft` is the FFT mesh in the unit cell used for the different MBPT quantities.

`green_mpw` is computed from [[ecut]], while `tchi_mpw` is defined by [[ecuteps]]
This two number are directly related to the memory footprint as they define of the size of the matrices that
must be stored in memory.

TODO: Describe other entries

Finally, we have the most importan section with the QP results in eV units:

```
================================================================================
 QP results (energies in eV)
 Notations:
     E0: Kohn-Sham energy
     <VxcDFT>: Matrix elements of Vxc[n_val] without non-linear core correction (if any)
     SigX: Matrix elements of Sigma_x
     SigC(E0): Matrix elements of Sigma_c at E0
     Z: Renormalization factor
     E-E0: Difference between the QP and the KS energy.
     E-Eprev: Difference between QP energy at iteration i and i-1
     E: Quasi-particle energy
     Occ(E): Occupancy of QP state



--- !GWR_SelfEnergy_ee
iteration_state: {dtset: 1, }
kpoint     : [   0.000,    0.000,    0.000, ]
spin       : 1
gwr_scf_iteration: 1
gwr_task   : G0W0
QP_VBM_band: 4
QP_CBM_band: 5
KS_gap     :    2.433
QP_gap     :    3.401
Delta_QP_KS:    0.969
data: !Tabular |
     Band       E0 <VxcDFT>     SigX SigC(E0)        Z     E-E0  E-Eprev        E   Occ(E)
        2   -0.300  -11.415  -13.532    1.594    0.860   -0.447   -0.447   -0.746    2.000
        3   -0.300  -11.415  -13.532    1.594    0.860   -0.447   -0.447   -0.746    2.000
        4   -0.300  -11.415  -13.532    1.594    0.860   -0.447   -0.447   -0.746    2.000
        5    2.133   -9.958   -4.937   -4.393    0.834    0.522    0.522    2.655    0.000
        6    2.133   -9.958   -4.937   -4.393    0.834    0.522    0.522    2.655    0.000
        7    2.133   -9.958   -4.937   -4.393    0.834    0.522    0.522    2.655    0.000
...
```

The meaning of the different columns should be self-explanatory.

As usual, we can use `abiopen.py` with the `-p` (`--print`) option to have a summary
of the most important results printed to screen:
Try, for instance:

```
abiopen.py tgwr_3o_DS1_GWR.nc -p
```

The last section printed to the screen gives ous the QP direct gaps:

```
=============================== GWR parameters ===============================
gwr_task: G0W0
Number of k-points in Sigma_{nk}: 2
Number of bands included in e-e self-energy sum: 100
ecuteps: 6.0
ecutsigx: 20.0
ecut: 8.0
gwr_boxcutmin: 1.1

============================ QP direct gaps in eV ============================
            kpoint     kname  ks_dirgaps  qpz0_dirgaps  spin
0  [0.0, 0.0, 0.0]  $\Gamma$    2.564094      3.663357     0
1  [0.5, 0.0, 0.0]         L    2.802929      3.848154     0
```

The calculation has produced three different text files:

tgwr\_3o\_SIGC\_IT:
Diagonal elements of $\Sigma_c(i \tau)$ in atomic units

tgwr\_3o\_SIGXC\_IW:
Diagonal elements of $\Sigma_\xc(i \omega)$ in eV units

tgwr\_3o\_SIGXC\_RW:
Diagonal elements of $\Sigma_\xc(\omega)$ in eV units and spectral function $A_\nk(\omega)$

Finally, we have a netcdf file named tgwr\_3o\_DS1_GWR.nc storing the same data
in binary format that can be easily post-processed with AbiPy using
[this script][https://abinit.github.io/abipy/gallery/plot_gwr.html].
Please take some time to read the script and understand how this post-processing tool works before proceeding.


### Extracting useful info from the GWR log file

To extract the memory allocated for the most memory demanding arrays, one can use:

```
grep "<< MEM" log"
```

To extract the wall-time and cpu-time for the most important sections, one can use:

```
grep "<< TIME" log"
```

To obtain the wall-time and cpu-time required by the different datasets, use:

```
grep "dataset:" log | grep "<< TIME"
```

### Convergence study HOWTO

As discussed in [[cite:Setten2017]], the convergence studies for
the $\kk$-mesh, [[nband], and the cutoff energies can be decoupled.
This means that one can start with a relatively coarse $\kk$-mesh to determine
the converged values of [[nband]], [[ecuteps]], and [[ecutsigx]],
then fix these values, and refine the BZ sampling only at the end.

The recommended procedure for converging GWR gaps is therefore as follows:

1) Initial step

- Select the k-points where QP gaps are wanted.
  Usually the VBM and the CBM so that one can use [[gwr_sigma_algo]] 2
- Fix the [[ngkpt]] $\kk$-mesh in the WFK file to a resonable value and produce "enough" [[nband]] states
- Set an initial value for [[gwr_ntau]] in the GWR run.
- Perform convergence studies on [[nband]], [[ecuteps]], and [[ecutsigx]].

Be aware that [[gwr_ntau]] = 6 may be too small, potentially leading to unphysical results
such as renormalization factors $Z$ below 0.6.
If this occurs, increase [[gwr_ntau]] in order to fix the problem.
If the number of [[nband]] states in the WFK file is not large enough,
go back to point 1) and generate a new WFK with more bands else proceeed with the next step.

2) Convergence wrt [[gwr_boxcutmin]]

- After determining converged values for [[nband]], [[ecuteps]], and [[ecutsigx]], begin refining [[gwr_boxcutmin]].
- Increase it gradually to control memory usage and CPU time, which increase rapidly with this parameter.

3) Convergence wrt [[gwr_ntau]]

Once the results are converged with respect to [[gwr_boxcutmin]], you may start to increase [[gwr_ntau]]
while adjusting the number of MPI processes accordingly (you are not using multidatasets, right?)

4) Convergence wrt [[ngkpt]]

Finally, refine the BZ sampling to ensure full convergence.

Once a good setup have been found, one can use the same parameters to compute the QP corrections in the IBZ
using [[gwr_sigma_algo]] 2 and [[gw_qprange]] = `-NUM` to have a `GWR.nc` file that can be used to
perform an interpolation of the GW band structure as discussed in the last part of this tutorial.

In the next sections, we explain how to perform these convergence studies and how to use AbiPy to analyze the results.
Note that we will not provide ready-to-use input files.
Your task is therefore to modify tgwr\_3.abi, run the calculation (preferably in parallel), and then analyze the results.

Note that, due to cancellations of errors, QP gaps that are differences between QP energies
are usually much easier to convergence than QP values.
Fortunately, absolute values are important only in rather specialized studies such as Work Function
or Band Alignment in Heterostructures.
In this tutorial, we only focus on gaps, and we aim to achieve an overall convergence of the QP gaps 0.01 eV (10 meV).
As a consequence we will try to reach a convergence of 2 meV

## Convergence wrt nband and ecuteps

To perform a double convergence study in [[nband]] and [[ecuteps]],
we can use [[udtset]], and add these variables to tgwr\_3.abi:

```sh
ndtset 16    # NB: ndtset must be equal to udtset[0] * udtset[1]
udtset 4 4

#inner loop: increase nband. The number of iterations is given by udtset[0]
nband:? 100    nband+? 50

#outer loop: increase ecuteps. The number of iterations is given by udtset[1]
ecuteps?: 6     ecuteps?+  2
```

If we analyze the wall-time required by each dataset, we observe
that, at variance with the conventional GW code,
the values of [[nband]] and [[ecuteps]] have little impact of the computational cost.

```
grep "dataset:" log | grep "<< TIME"
 dataset: 1 , wall:  8.84 [s] , cpu:  8.83 [s] <<< TIME
 dataset: 2 , wall:  9.04 [s] , cpu:  9.03 [s] <<< TIME
 dataset: 3 , wall:  9.29 [s] , cpu:  9.24 [s] <<< TIME
 dataset: 4 , wall:  9.99 [s] , cpu:  9.98 [s] <<< TIME
 dataset: 5 , wall:  8.80 [s] , cpu:  8.79 [s] <<< TIME
 dataset: 6 , wall:  9.19 [s] , cpu:  9.15 [s] <<< TIME
 dataset: 7 , wall:  9.45 [s] , cpu:  9.45 [s] <<< TIME
 dataset: 8 , wall: 10.32 [s] , cpu: 10.28 [s] <<< TIME
 dataset: 9 , wall:  8.81 [s] , cpu:  8.81 [s] <<< TIME
 dataset: 10 , wall:  9.28 [s] , cpu:  9.22 [s] <<< TIME
 dataset: 11 , wall:  9.60 [s] , cpu:  9.59 [s] <<< TIME
 dataset: 12 , wall: 10.19 [s] , cpu: 10.17 [s] <<< TIME
 dataset: 13 , wall:  8.95 [s] , cpu:  8.94 [s] <<< TIME
 dataset: 14 , wall:  9.30 [s] , cpu:  9.29 [s] <<< TIME
 dataset: 15 , wall:  9.60 [s] , cpu:  9.60 [s] <<< TIME
 dataset: 16 , wall: 10.29 [s] , cpu: 10.28 [s] <<< TIME
```

As concerns [[nband]] this is expected, as this parameter enters into play
only during the computation of the Green's function.
For [[ecuteps]] this is due to the FFT algorithm.
Of course, we are not saying that the wall-time does not depend on these two parameters.
Increasing [[ecuteps]] definitely has a greater effect in terms of wall-time.

Once the calculation is finished, we can use the GwrRobot to load the list of GWR.nc files
and plot the convergence of the direct QP gaps with the following piece of code:

```python
#!/usr/bin/env python

udtset = [4, 4]  # Should be equal to the grid used in the ABINIT input file

filepaths = [f"tgwr_3o_DS{i}{j}_GWR.nc"
             for i in range(1, udtset[0]+1) for j in range(1, udtset[1]+1)]

from abipy.electrons.gwr import GwrRobot
robot = GwrRobot.from_files(filepaths)

# Convergence window in eV
abs_conv = 0.002

robot.plot_qpgaps_convergence_new(x="nband", y="qpz0_dirgaps",
                                  abs_conv=abs_conv, hue="ecuteps")

robot.plot_qpgaps_convergence_new(x="ecuteps", y="qpz0_dirgaps",
                                  abs_conv=abs_conv, hue="nband")
```

Save the script in a file, let’s say `conv_nband_ecuteps.py`, in the same
directory where we have launched the calculation, make the file executable with:

```
chmod u+x conv_nband_ecuteps.py
```

and then run it with:

```
./conv_nband_ecuteps.py
```

In the [first GW tutorial](/tutorial/gw1), we have already performed convergence studies,
and [[nband]] = 100 was found to give results converged within 30 mev, which is fair to compare with experimental accuracy.
Also ecuteps = 6.0 can be considered converged within 10 meV.
We will not repeat these convergence studies here, we just use these values for our calculation so that
we can focus on the analysis of the output file and the post-processing of the results.

Since we will continue using `tgwr_3.abi` to perform other convergence studies,
it is a good idea to move the previously generated GWR.nc files and the python script
to a new directory to avoid overlaps between different calculations.
You can, for example, use the following list of shell commands:

```
mkdir conv_nband_ecuteps
mv conv_nband_ecuteps.py conv_nband_ecuteps
mv tgwr_3o_* log conv_nband_ecuteps
```

## Convergence wrt gwr_boxcutmin

At this point edit `tgwr_3.abi`, replace the values of [[nband]] and [[ecuteps]] found earlier,
and define a new one-dimensional multidataset in which we increase [[gwr_boxcutmin]].

```
ndtset 5
gwr_boxcutmin: 1.1
gwr_boxcutmin+ 0.2
```

Note that different values of [[gwr_boxcutmin]] can lead to the same value of `g_ngfft`
as not all the number of points along the three directions of the FFT box are supported by the FFT library.

```
grep QP_gap tgwr_3.abo
````

```
grep g_nfft tgwr_3.abo
```

```
grep "dataset:" log | grep "<< TIME"
```

Once the calculation is finished, save the script reported below in a file,
let’s say `conv_boxcut.py`, in the same directory where we launched the calculation,
make it executable and execute it as usual


```python
#!/usr/bin/env python
filepaths = [f"tgwr_3o_DS{i}_GWR.nc" for i in range(1, 6)]

from abipy.electrons.gwr import GwrRobot
robot = GwrRobot.from_files(filepaths)

abs_conv = 0.01
robot.plot_qpgaps_convergence_new(x="gwr_boxcutmin", y="qpz0_dirgaps",
                                  abs_conv=abs_conv)

#robot.plot_qpgaps_convergence(sortby=sortby, abs_conv=abs_conv)
```

You should obtain the following figure:

![](base1_assets/abiopen_tgwr1_1.png)

In the case of silicon, the results are rather insensitive to the density
of the FFT mesh but keep in mind that other systems may behave differently.
We therefore fix [[gwr_boxcutmin]] to 1.1 and proceed with the other convergence studies.

## Convergence wrt gwr_ntau

Change `tgwr_3.abi` to perform a convergence study wrt [[gwr_ntau]] by adding the following section

```
ndtset      15
gwr_ntau:   6         # Number of imaginary-time points
gwr_ntau+   2
```

Once the calculation is finished, save the script reported below in a file,
let’s say `conv_ntau.py`, in the same directory where we launched the calculation,
make it executable and execute it as usual


```python
#!/usr/bin/env python
filepaths = [f"tgwr_3o_DS{i}_GWR.nc" for i in range(1, 16)]

from abipy.electrons.gwr import GwrRobot
robot = GwrRobot.from_files(filepaths)

def sortby(gwr):
    r"""$N_\tau$"""
    return gwr.params["gwr_ntau"]


abs_conv = 0.01
robot.plot_qpgaps_convergence_new(x="gwr_ntau", y="qpz0_dirgaps",
                                  abs_conv=abs_conv)

robot.plot_qpgaps_convergence(sortby=sortby, abs_conv=abs_conv)
```

You should get the following figure:

![](base1_assets/abiopen_tgwr1_1.png)

Clearly 6 minimax points are not enough to convergence.
Also the convergence with [[gwr_ntau]] is not variational, and there are meshes that perform better than others.
This study shows that we should use XXX minimax points in order to reach an accuracy of 0.01 eV (`abs_conv`).
Unfortunately, this would render the calculations more expensive, especially for a tutorial,
therefore we opt to continue with XXX minimax points for the next calculations


```
mkdir conv_ntau
mv conv_ntau.py conv_ntau
mv tgwr_3o_* log conv_ntau
```

### GWR calculations with different BZ meshes.

In this part of the tutorial, we run GWR calculations with different BZ meshes.

You may want to start immediately the computation by issuing:

```
mpirun -n 1 tgwr_4.abi > tgwr_4.log 2> err &
```

with the following input file:

{% dialog tests/tutorial/Input/tgwr_4.abi %}

The output file is reported here for your convenience:

{% dialog tests/tutorial/Refs/tgwr_4.abo %}

To analyze multiple GWR.nc files, we can use the `GwrRobot`, and the following python script:

```python
#!/usr/bin/env python
from abipy.electrons.gwr import GwrRobot

filepaths = [
    "tgwr_3o_GWR.nc",
    "tgwr_4o_GWR.nc",
]

robot = GwrRobot.from_files(filepaths)
#print(robot)

kpoint = (0, 0, 0)
df_dirgaps = robot.get_dirgaps_dataframe(kpoint=kpoint)
print(df_dirgaps)
```

To conclude, our best values for the direct QP gaps are ...

```
# WFK file with empty states. Must be consisten with ngkpt, shiftk
ndtset  3
ngkpt1   2 2 2
getwfk_filepath1 "tgwr_2o_DS1_WFK"

ngkpt2   4 4 4
getwfk_filepath2 "tgwr_2o_DS2_WFK"

ngkpt3   6 6 6
getwfk_filepath3 "tgwr_2o_DS3_WFK"
```

You should get the following figure:

![](base1_assets/abiopen_tgwr1_1.png)

### How to compute an interpolated band structure with GWR

In this last part of the tutorial, we discuss how to interpolate the QP corrections
along an arbitrary $\kk$-path using the star-function method discussed in
[this section](tutorial/eph_intro/#star-function-interpolation-of-the-ks-eigenvalues) of the EPH introduction.
This method is less accurate the e.g. the Wannier interpolation and may have problems in the presence
of band-crossings, but it has the big advantage of being much easier to use and with mininal user intervention.

First of all, we need to compute the QP corrections for **all** the $\kk$-points in the IBZ.
This is done in the following input:

{% dialog tests/tutorial/Input/tgwr_5.abi %}

Note the usage of  [[gw_qprange]] = `-NUM` and [[gwr_sigma_algo]] 1 to activate
the supercell algorithm for $\Sigma_\nk$ (the most efficient algorithm in this particular case).

To keep the wall-time at a reasonable level, we use a WFK file with 2x2x2 $\kk$-mesh,
but it is clear that a more accurate intepolation would require denser $\kk$-meshes.

Run the job with:

```
mpirun -n 1 tgwr_5.abi > tgwr_5.log 2> err &
```

The output file is reported here for your convenience:

{% dialog tests/tutorial/Refs/tgwr_5.abo %}

Now use the following AbiPy script to read the GWR results and interpolate the QP corrections.
In this case, we pass the KS band structure stored in `tgwr_1o_DS2_GSR.nc`.
This is the recommended approach as AbiPy will interpolate the QP corrections rather that the QP energies.
The interpolated QP corrections will then be added to the KS energies
thus improving the stability of the interpolation as QP corrections are usually smoother than QP energies.

```python
#!/usr/bin/env python
from abipy.electrons.gwr import GwrFile

gwr = GwrFile("tgwr_5o_GWR.nc")
r = gwr.interpolate(ks_ebands_kpath="tgwr_1o_DS2_GSR.nc")

from abipy.electrons.ebands import ElectronBandsPlotter
plotter = ElectronBandsPlotter()
plotter.add_ebands("KS bands", r.ks_ebands_kpath)
plotter.add_ebands("GW bands", r.qp_ebands_kpath)

#r.qp_ebands_kpath.plot(with_gaps=True)

plotter.combiplot()
#plotter.gridplot()
```

You should get the following figure:

![](base1_assets/abiopen_tgwr1_1.png)

### Analyzing the convergence with the minimax mesh

The script below shows how to use the `GWRobot` to perform a convergence study
for the QP direct gaps with respect to the number of minimax points.

```python
#!/usr/bin/env python
from abipy.electrons.gwr import GwrRobot
filepaths = [f"tgwr_3o_DS{num}_GWR.nc" for num in range(1, 15)]

robot = GwrRobot.from_files(filepaths)

kpoint = (0, 0, 0)
df_dirgaps = robot.get_dirgaps_dataframe(kpoint=kpoint)
print(df_dirgaps)

# This function is used to sort the GWR.nc files by the number of minimax points
# The doc string will be used as xlabel.

def sortby(gwr):
    r"""$N_\tau$"""
    return gwr.params["gwr_ntau"]

robot.plot_qpgaps_convergence(sortby=sortby, abs_conv=0.01)
```

If we analyze the results of a more expensive run done with multiple minimax meshes
and then plot the results with AbiPy using the same script but a different `filepaths`,
we obtain the following results:

Clearly 6 minimax points are not enough to convergence.
Also the convergence with [[gwr_ntau]] is not variational, and there are meshes that perform better than others.
This study shows that we should use XXX minimax points in order to reach an accuracy of 0.01 eV (`abs_conv`).
Unfortunately, this would render the calculations more expensive, especially for a tutorial,
therefore we opt to continue with XXX points.

The script below shows how to analyze a single GWR.nc file.

```python
!#/usr/bin/env python

from abipy.electrons.gwr import GwrFile
gwr = GwrFile("tgwr_3o_GWR.nc")
print(gwr)

spin = 0
kpoint = (0, 0, 0)
df_dirgaps = gwr.get_dirgaps_dataframe()
print(df_dirgaps)

df = gwr.get_dataframe_sk(spin=spin, kpoint=kpoint)
print(df)

gwr.plot_sigma_imag_axis(kpoint=kpoint)

gwr.plot_sigma_real_axis(kpoint=kpoint)

gwr.plot_spectral_function(kpoint=kpoint)

#gwr.plot_qps_vs_e0(with_fields="all")
```
