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

This tutorial should take about 2 hours.

[TUTORIAL_README]

### Ground-state and band structure

*Before beginning, you might consider creating a different subdirectory to work in. Why not create Work_gwr?*

The file *tgwr_1.abi* is the input file for the first step:
a SCF run followed by a KS NSCF band structure calculation along a high-symmetry $\kk$-path.
Copy it to the working directory with:

```sh
mkdir Work_gwr
cd Work_gwr
cp ../tgwr_1.abi .
```

You may want to immediately start the job in background with:

```sh
mpirun -n 1 tgwr_1.abi > tgwr_1.log 2> err &
```

so that we have some time to discuss the input while ABINIT is running.

{% dialog tests/tutorial/Input/tgwr_1.abi %}

!!! tip

    All the input files of this tutorial can be executed in parallel by just
    changing increasing the value of n in `mpirun -n`.
    Clearly the size of the problem is very small so do not expect great performance
    if you start to use dozens of MPI processes.

The first dataset produces the KS density file that is then used to compute the band structure in the second dataset.
Since all the input files of this tutorial use the same crystalline structure, pseudos and [[ecut]],
we declare these variables in an external file that will be included in all the other input files
using the syntax:

```
# Include geometry and pseudos
include "gwr_include.abi"
```

If you open the include file:

{% dialog tests/tutorial/Input/gwr_include.abi %}

you will notice that we are using norm-conserving (NC) pseudos taken from
the standard scalar-relativistic table of the [PseudoDojo](https://www.pseudo-dojo.org/).

For GW calculations, we strongly recommend using NC pseudos from the **stringent table**
as these pseudos have closed-shells treated as valence states.
This is important for a correct description of the matrix elements of the exchange part of the self-energy
that are rather sensitive to the overlap between the wavefunctions [[cite:Setten2018]].
In the special case of silicon, there is no difference between the standard version and the stringent version of the pseudo,
but if we consider e.g. Ga, you will notice that the stringent version includes all the 3spd states in valence
besides the outermost 4sp electrons whereas the standard pseudo for Ga designed for GS calculations includes only the 3d states.
Since we are using scalar-relativistic (SR) pseudos, it is not possible to treat spin-orbit coupling (SOC)
and [[nspinor]] must be set to 1 (default value).
To include SOC with NC, one should use fully-relativistic (FR) pseudos, set [[nspinor]] to 2 in all the input files
and then select the appropriate value of [[nspden]] (4 or 1) depending on whether the system is magnetic or not.

For efficiency reasons, in these examples we use underconverged parameters:
the cutoff energy [[ecut]] is set to XXX that is slightly smaller that the recommended value reported in the PseudoDojo table
Also, the density is computed with a 2x2x2 $\Gamma$-centered $\kk$-mesh.
In real life, one should compute the KS density with a much denser $\kk$-mesh so that we have an accurate $n(\rr)$
to compute KS states and eigenvalues for the GWR part.

One final comment on MPI parallelism:
In the input, we set [[paral_kgb]] to 1 to enable k-point, band, and FFT parallelization,
and [[autoparal]] to 1 to allow ABINIT to determine an optimal distribution, with [[npfft]] fixed at 1.

At this point, the calculation should have completed, and we can have a look
at the electronic band structure to understand the position of the band edges.

!!! tip

    If |AbiPy| is installed on your machine, you can use the |abiopen| script
    with the `--expose` option to visualize the band structure from the GSR.nc file:

        abiopen.py tgwr_1o_DS2_GSR.nc --expose

    ![](base1_assets/abiopen_tgwr1_1.png)

    To print the results to terminal, use:

        abiopen.py tgwr_1o_DS2_GSR.nc -p

Silicon is an indirect band gap semiconductor:
the CBM is located at the $\Gamma$ point
while the VMB is at the  ...

The PBE KS fundamental band gap is XXX that is strongly underestimated with respect to experiment.

!!! important

    Similarly to the conventional GW code, also GWR can compute QP corrections only
    for the $\kk$-points belonging to the $\kk$-mesh used for generating the WFK file.
    Before running GW calculations is always a good idea to analyze carefully the KS band
    structure to understand the location of the band edges and select the most appropriate $\kk$-mesh.


### Generation of the WFK file with empty states

In the second input file, we use the density file computed previously (tgwr\_1o\_DS1\_DEN)
to generate two WFK files with two different $\Gamma$-centered $\kk$-meshes ([[ngkpt]] = 2x2x2 and 4x4x4).
This allows us to perform convergence studies with respect to the BZ sampling in the last part of this tutorial.
Let us recall that shifted $\kk$-meshes are not supported by GWR.

Let's start the job in background with:

```sh
mpirun -n 1 tgwr_2.abi > tgwr_2.log 2> err &
```

{% dialog tests/tutorial/Input/tgwr_2.abi %}

Note that here we are using [[gwr_task]] = "HDIAGO" in order to perform a **direct diagonalization**
of the KS Hamiltonian constructed from the input DEN file.
This procedure differs from the one used in the other GW tutorials in which the WFK file
is generated by performing an **iterative diagonalization** in which only the application of the Hamiltonian is required.
The reason is that the direct diagonalization can outperform iterative methods when many empty states are required,
especially if one can take advantage of ScalaPack to distribute the KS Hamiltonian matrix.

Here, we request 100 states, as in the previous GW tutorial, where [[nband]] = 100
was considered converged within 30 meV.
Clearly, when studying new systems, the converged nband value is not known beforehand.
Therefore, it is important to plan ahead and choose a reasonably large number of states
to avoid regenerating the WFK file multiple times just to increase nband

Note that [[paral_kgb]] = 1 is only available in the ground-state (GS) part.
The GWR code employs its own distribution scheme, which depends on the value of [[gwr_task]].
When ‘HDIAGO’ is used, the distribution is handled automatically at runtime,
and the user has no control over it.
Refer to the note below for guidance on choosing an appropriate number of MPI processes
to ensure an efficient workload distribution.


!!! important

    The direct diagonalization is MPI-parallelized across three different levels:
    collinear spin $\sigma$ (not used here), $\kk$-points in the IBZ and scalapack distribution of the $H^\sigma_\kk(\bg,\bg')$ matrix.
    Abinit will try to find an "optimal" distribution of the workload at runtime, yet there are a couple
    of things worth keeping in mind when choosing the number of MPI processes for this step.
    Ideally the total number of cores should be a multiple of [[nkpt]] * [[nsppol]] to avoid load imbalance.

    In order to compute **all** the eigenvectors of the KS Hamiltonian, one can use gwr_task "HDIAGO_FULL".
    In this case the value of [[nband]] is automatically set to the total number of plawewaves
    for that particular $\kk$-point.
    No stopping criterion such as [[tolwfr]] or number of iterations [[nstep]] are required when Scalapack is used.
    Again, multi-datasets are **strongly discouraged** if you care about performance.


### Our first GWR calculation

For our first GWR run, we use a minimalistic input file that
performs a GWR calculation using the DEN and the WFK file produced previously.
First of all, you may want to start immediately the computation by issuing:

```sh
mpirun -n 1 tgwr_3.abi > tgwr_3.log 2> err &
```

with the following input file:

{% dialog tests/tutorial/Input/tgwr_3.abi %}

This input file contains some variables whose meaning is the same as in the conventional GW code,
and other variables specific to GWR.

We use [[optdriver]] 6 to enter the GWR code and [[gwr_task]] activates a one-shot GW calculation.
We ask for a minimax mesh with [[gwr_ntau]] = 6 points, the minimum number of points one can use.
Clearly, this parameter should be subject to convergence studies.
[[getden_filepath]] specifies the density file used to compute $v_{xc}[n](\rr)$,
while [[getwfk_filepath]] specifies the WFK file with empty states used to build the Green's function.
Keep in mind that the $\kk$-mesh specified in the input via [[ngkpt]], [[nshiftk]] and [[shiftk]] must
agree with the one found in the WFK file else the code will abort.

[[gwr_boxcutmin]] is set to 1.1 in order to accelerate the calculation and reduce memory requirements.
This is another parameter that should be subject to convergence studies as systems with localized electrons
such as 3d or 4f electrons may require larger values (denser FFT meshes).
Please take some time to read the variable description of [[gwr_boxcutmin]] before proceeding.

Now, let us turn our attention to the variables that are also used in the conventional GW code.
In the [first GW tutorial](/tutorial/gw1), we have already performed convergence studies,
and [[nband]] = 100 was found to give results converged within 30 mev, which is fair to compare with experimental accuracy.
Also ecuteps = 6.0 can be considered converged within 10 meV.
We will not repeat these convergence studies here, we just use these values for our calculation so that
we can focus on the analysis of the output file and the post-processing of the results.

The cutoff of the polarizability and $W$ is defined by [[ecuteps]] as in the conventional GW code.
The cutoff for the exchange part of the self-energy is defined by [[ecutsigx]].
Theoretically, we need [[ecutsigx]] = 4 [[ecut]] to have an exact treatment of $\Sigma_x$.
The computation of the exchange is relatively fast as only occupied states are involved.

The $\kk$-points and the band range for the QP corrections are set expliclty [[nkptgw]], [[kptgw]] and [[bdgw]]
Alternatively, one can use [[gw_qprange]]

For the spectral function, we have the following variables [[nfreqsp]], [[freqspmax]]

TODO: Mention [[inclvkb]] == 2, [[symsigma]]

We can now have a look at the main output file:

{% dialog tests/tutorial/Refs/tgwr_3.abo %}

First of all, we have a section that summarizes the most important GWR parameters
are summarized in this Yaml document

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

Finally, we have the section with the QP results in eV units

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

The job has produced three text files:

tgwr\_3o\_SIGC\_IT:
Diagonal elements of $\Sigma_c(i \tau)$ in atomic units

tgwr\_3o\_SIGXC\_IW:
Diagonal elements of $\Sigma_\xc(i \omega)$ in eV units

tgwr\_3o\_SIGXC\_RW:
Diagonal elements of $\Sigma_\xc(\omega)$ in eV units and spectral function $A_\nk(\omega)$

Finally, we have a netcdf file named tgwr\_3o\_GWR.nc that stores the same data
in binary format and that can be easily post-processed with AbiPy as discussed in the next section.

### Analyzing the GWR.nc file with AbiPy

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

### Convergence study HOWTO

As discussed in [[cite:Setten2017]], the convergence studies for
the $\kk$-mesh, [[nband], and the cutoff energies can be decoupled.
This means that one can start with a relatively coarse $\kk$-mesh to determine
the converged values of [[nband]], [[ecuteps]], and [[ecutsigx]],
then fix these values, and refine the BZ sampling.
The recommended procedure for converging GWR calculations is as follows:

1) Initial convergence tests

- Fix the [[ngkpt]] $\kk$-mesh in the WFK file to a resonable value and produce "enough" [[nband]] states
- Set an initial value for [[gwr_ntau]] in the GWR run.
- Perform convergence studies on [[nband]], [[ecuteps]], and [[ecutsigx]].

Be aware that [[gwr_ntau]] = 6 may be too small, potentially leading to unphysical results
such as renormalization factors $Z$ below 0.6.
If this occurs, increase [[gwr_ntau]] in order to fix the problem.
If the number of [[nband]] states in the WFK file is not enough,
go back to point 1) and produce a new WFK with more bands else proceeed with the next step.

2) Convergence wrt [[gwr_boxcutmin]]

- After determining converged values for [[nband]], [[ecuteps]], and [[ecutsigx]], begin refining [[gwr_boxcutmin]].
- Increase it gradually (e.g., 1.1, 1.2, 1.3) to control memory usage and CPU time, which increase rapidly with this parameter.

Once the results are converged with respect to [[gwr_boxcutmin]], you may start to increase [[gwr_ntau]]
while adjusting the number of MPI processes accordingly (you are not using multidatasets, right?)
Finally, refine the BZ sampling to ensure full convergence.
Note that, due to cancellations of errors, QP gaps (differences between QP energies)
are much easier to convergence than the QP values.

### GWR calculations with two different BZ meshes.

In this part of the tutorial, we run two GWR calculations with two different BZ meshes.
This allows us to simulate a typical convergence study and discuss how to use
AbiPy and pandas dataframes to analyze the results.

You may want to start immediately the computation by issuing:

```sh
mpirun -n 1 tgwr_4.abi > tgwr_4.log 2> err &
```

with the following input file:

{% dialog tests/tutorial/Input/tgwr_4.abi %}

The output file is reported here for your convenience:

{% dialog tests/tutorial/Refs/tgwr_4.abo %}

To analyze multiple GWR.nc files we can use the `GwrRobot`, and the following python script

```python
#!/usr/bin/env python
from abipy.electrons.gwr import GwrRobot

filepaths = [
    "tgwr_3o_GWR.nc",
    "tgwr_4o_GWR.nc",
]

robot = GwrRobot.from_files(filepaths)
#print(robot)

df_dirgaps = robot.get_dirgaps_dataframe()
print(df_dirgaps)

kpoint = (0, 0, 0)
spin = 0
```

### How to compute an interpolated band structure with GWR

In this last part of the tutorial, we discuss how to interpolate the QP corrections
along an arbitrary $\kk$-path using the star-function method discussed in
[this section](tutorial/eph_intro/#star-function-interpolation-of-the-ks-eigenvalues) of the EPH introduction.

First of all, we need to compute QP corrections for all the $\kk$-points in the IBZ.
This is what is done in the following input:

{% dialog tests/tutorial/Input/tgwr_5.abi %}

Note the use [[gw_qprange]] = `-NUM` to compute QP corrections for all the $\kk$-points in the IBZ.
and [[gwr_sigma_algo]] 1 to use the supercell algorithm for $\Sigma_\nk$
(the most efficient one in this particular case).

Run the calculation with:

```sh
mpirun -n 1 tgwr_5.abi > tgwr_5.log 2> err &
```

The output file is reported here for your convenience:

{% dialog tests/tutorial/Refs/tgwr_5.abo %}

Now use the following AbiPy script to read the GWR results and interpolate the QP corrections

```python
#!/usr/bin/env python
from abipy.electrons.gwr import GwrFile

gwr = GwrFile("tgwr_5o_GWR.nc")
r = gwr.intepolate(ebands_kpath="tgwr_1o_DS2_GSR.nc")

r.ks_ebands_kpath.plot()

from abipy.electrons.ebands import ElectronBandsPlotter
plotter = ElectronBandsPlotter()
#plotter.add_ebands("KS bands", "tgwr_1o_DS2_GSR.nc")
plotter.add_ebands("KS bands", r.ks_ebands_kpath)
plotter.add_ebands("GW bands", r.qp_ebands_kpath)

plotter.gridplot()
plotter.combiplot()
```
