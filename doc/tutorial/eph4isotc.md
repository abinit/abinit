---
authors: MG
---

# Superconducting properties within the isotropic Eliashberg formalism

This tutorial explains how to compute phonon linewidths in metals due to the electron-phonon interaction (e-ph)
and how to use the McMillan equation to estimate the (conventional) superconducting critical temperature $T_c$
within the isotropic Eliashberg formalism.
We start by presenting the basic equations implemented in the code and their connection with the ABINIT variables.
Then we discuss how to run isotropic $T_c$-calculations and how to perform typical convergence studies
using hexagonal MgB$_2$ as example.
For a more complete theoretical introduction see [[cite:Giustino2017]] and references therein.

It is assumed the user has already completed the two tutorials [RF1](rf1) and [RF2](rf2),
and that he/she is familiar with the calculation of ground state and vibrational properties **in metals**.
The user should have read the [fourth lesson on Al](base4) as well
as the [introduction page for the EPH code](eph_intro) before running these examples.

This lesson should take about 1.5 hour.

## Formalism and connection with the implementation

Due to the interaction with electrons, phonons acquire a **finite lifetime**
given by the imaginary part of the phonon-electron self-energy $\Pi_\qnu$.
In the so-called double-delta approximation, the linewidth $\gamma_{\qq\nu}$
of the ${\qq\nu}$ phonon is given by:

$$
\gamma_{\qq\nu} = 2\pi \ww_{\qq\nu} \sum_{mn\kk} |g_{mn\nu}(\kk, \qq)|^2
\delta(\ee_{\kpq m} -\ee_F) \delta(\ee_{\kk n} -\ee_F)
$$

where $\ww_{\qq\nu}$ is the phonon frequency,
the sum over the electron wavevector $\kk$ runs over the full BZ, $\ee_F$ is the Fermi level
and $g_{mn\nu}(\kk, \qq)$ are the e-ph matrix elements discussed in the [EPH introduction](eph_intro).

For a given phonon wavevector $\qq$, the double Dirac delta restricts the BZ integration to
transitions between $\kk$ and $\kq$ electronic states on the Fermi surface (FS).
Converging the double-delta integral therefore requires very dense $\kk$-meshes in order
to sample enough states around the FS.
The convergence rate, indeed, is expected to be much slower than the one required by the electron DOS at $\ee_F$:

$$
g(\ee_F) = \sum_{n\kk} \delta(\ee_F - \ee_{n\kk})
$$

in which a single Dirac delta is involved.

!!! note

    It is worth stressing that there is another important contribution to the phonon lifetimes induced by
    **non-harmonic** terms in the Taylor expansion of the Born-Oppenheimer energy surface around the equilibrium point.
    In the many-body language these non-harmonic leads to **phonon-phonon scattering processes** that can give
    a substantial contribution.
    In the rest of the tutorial, however, non-harmonic terms will be ignored and we will be mainly focusing
    on the computation of the imaginary part of the phonon self-energy $\Pi$ in the harmonic approximation.

At the level of the implementation, the integration of the double delta can be performed
either with the **tetrahedron** scheme or by replacing the Dirac distribution with a **Gaussian** of finite width.
The integration algorithm is defined by [[eph_intmeth]].
If the Gaussian method is selected ([[eph_intmeth]] == 1), one can choose between a constant broadening $\sigma$
specified in Hartree via [[eph_fsmear]] or an adaptive scheme (activated when [[eph_fsmear]] < 0)
in which a state-dependent broadening $\sigma_\nk$ is automatically computed from
the electron group velocities $v_\nk$ [[cite:Li2015]].

!!! important

    The tetrahedron method is more accurate and does not require any broadening parameter.
    Note, however, that in the present implementation the computational cost of the double
    delta with the tetrahedron method
    quickly increases with the size of the $\kk$-mesh so the adaptive Gaussian scheme represents a valid alternative,
    especially when a dense $\kk$-sampling is used.

The value of the Fermi level, $\ee_F$, is automatically computed from the KS eigenvalues stored
in the input WFK file according to the two input variables [[occopt]] and [[tsmear]].
These parameters are **usually equal** to the ones used for the GS/DFPT calculation
although it is possible to change the value of $\ee_F$ during an EPH calculation using
three (mutually exclusive) input variables: [[eph_fermie]], [[eph_extrael]] and [[eph_doping]].

The sum over bands can be restricted to Bloch states within an energy window of thickness [[eph_fsewin]]
around the Fermi level.
The $\kk$-mesh is defined by the input variables [[ngkpt]], [[nshiftk]] and [[shiftk]].
These parameters **must be equal** to the ones used to generate the input WFK file.
Convergence studies are performed by generating different WFK files on $\kk$-meshes of increasing density.

The code computes $\gamma_{\qq\nu}$ for each $\qq$-point in the IBZ associated to
the [[eph_ngqpt_fine]] $\qq$-mesh.
The DFPT potentials and phonon frequencies are interpolated starting from the [[ddb_ngqpt]] $\qq$-mesh
used to generate the input DVDB/DDB files.
If [[eph_ngqpt_fine]] is not specified, the [[ddb_ngqpt]] $\qq$-mesh is used and no interpolation
of the potentials in $\qq$-space in required.
Once the phonon linewidths $\gamma_{\qq\nu}$ are known in the IBZ,
the code computes the Eliashberg function defined by:

$$
\alpha^2F(\ww) = \dfrac{1}{N_F} \sum_{\qq\nu} \dfrac{\gamma_{\qq\nu}}{\ww_{\qq\nu}} \delta(\ww - \ww_{\qq \nu})
$$

where $N_F$ is the density of states (DOS) per spin at the Fermi level.
$\alpha^2F(\ww)$ gives the strength by which a phonon of energy $\ww$ scatters electronic
states on the FS (remember that ABINIT uses atomic units hence $\hbar = 1$).
This quantity is accessible in experiments and experience has shown that
$\alpha^2F(\ww)$ is qualitatively similar to the phonon DOS $F(\ww)$:

$$
F(\ww) = \sum_{\qq\nu} \delta(\ww - \ww_{\qq \nu})
$$

This is not surprising as the equation for $\alpha^2F(\ww)$ resembles the one for $F(\ww)$,
except for the weighting factor $\frac{\gamma_{\qq\nu}}{\ww_{\qq\nu}}$.

The technique used to compute $\alpha^2F(\ww)$ is defined by [[ph_intmeth]] (note the `ph_` prefix).
Both the Gaussian ([[ph_intmeth]] = 1 with [[ph_smear]]) and
the linear tetrahedron method [[ph_intmeth]] = 2, default) are implemented for the $\qq$-space integration
The total e-ph coupling strength $\lambda$ is defined as the first inverse moment of $\alpha^2F(\ww)$:

$$
\lambda = \int \dfrac{\alpha^2F(\ww)}{\ww}\dd\ww = \sum_{\qq\nu} \lambda_{\qq\nu}
$$

where we have introduced the mode dependent coupling strength:
<!-- For spin unpolarized systems: -->

$$
\lambda_{\qq\nu} = \dfrac{\gamma_{\qq\nu}}{\pi N_F \ww_{\qq\nu}^2}
$$

In principle, $T_c$ can be obtained by solving the isotropic Eliashberg equations
but many applications bypass the solution of these equations and
**estimate** $T_c$ using the McMillan equation:

$$
T_c = \dfrac{\ww_{log}}{1.2} \exp \Biggl [ \dfrac{-1.04 (1 + \lambda)}{\lambda ( 1 - 0.62 \mu^*) - \mu^*} \Biggr ]
$$

where $\mu^*$ describes the (screened) e-e interaction and
$\ww_{\text{log}}$ is the *logarithmic* average of the phonon frequencies defined by:

$$
\ww_{\text{log}} = \exp \Biggl [ \dfrac{2}{\lambda} \int \dfrac{\alpha^2F(\ww)}{\ww}\log(\ww)\dd\ww \Biggr ]
$$

In pratical applications, $\mu^*$ is treated as an **external parameter**.
The default value of [[eph_mustar]] used in the code is 0.1.

<!--
that usually ranges between .. and ,,

Nesting factor

\begin{equation}
    N(\qq) = \sum_{mn\kk} \delta(\ee_{\kpq m}) \delta(\ee_{\kk n})
\end{equation}

Implementation details:

\begin{equation}
    \gamma_{\qq\nu} = 2\pi \ww_{\qq\nu} \sum_{pp'} d_{\qq p}^* \tilde\gamma_{p p'}(\qq) d_{\qq p'}^*
\end{equation}

where

\begin{equation}
    \tilde\gamma_{p p'}(\qq) =
    \sum_{mn\kk} g_{mn,p}(\kk, \qq)^*  g_{mn,p'}(\kk, \qq)  \delta(\ee_{\kpq m}) \delta(\ee_{\kk n})
\end{equation}

$\tilde\gamma_{p p'}(\qq)$ has the same symmetries as the dynamical matrix.
The code computes $\tilde\gamma$ for all q-points in the IBZ, each matrix is re-symmetrized
if symdynmat == 1 so that degeneracies at high-symmetry q-points are correctly reproduced.
The matrix elements for q-points in the full BZ are then obtained by rotating the initial set of q-points in the BZ.
and Fourier transformed to real space with:

\begin{equation}
    \tilde\gamma_{p p'}(\RR) = \sum_\qq e^{i\qq\cdot\RR} \tilde\gamma_{p p'}(\qq)
\end{equation}

At this point, it is possible to interpolate the matrices via Fourier interpolation.
A similar approach is used for the dynamical matrix.
Note that we interpolate the matrices in this representation instead of the phonon mode representation
to avoid numerical instabilities introduce by band crossings.
-->

## Getting started

[TUTORIAL_READMEV9]

Before beginning, you might consider to work in a different subdirectory as for the other tutorials.
Why not create Work_eph4isotc in $ABI_TESTS/tutorespfn/Input?

```sh
cd $ABI_TESTS/tutorespfn/Input
mkdir Work_eph4isotc
cd Work_eph4isotc
```

In this lesson we prefer to focus on e-ph calculations and the associated convergence studies.
For this reason, we rely on **pre-computed** DEN.nc, DDB and DFPT potentials
to bypass both the GS and the DFPT part.
The DEN.nc file will be used to perform NSCF computations on arbitrarily dense $\kk$-meshes while the
DFPT POT.nc files will be merged with the *mrgdv* utility to produce the DVDB database required by the EPH code.

Note that these files **are not shipped** with the official ABINIT tarball as they are relatively large in size.
In order to run the examples of this tutorial, you need to download the files from
[this github repository](https://github.com/abinit/MgB2_eph4isotc).

If git is installed on your machine, one can easily fetch the entire repository with:

```sh
git clone https://github.com/abinit/MgB2_eph4isotc.git
```

Alternatively, use *wget*:

```sh
wget https://github.com/abinit/MgB2_eph4isotc/archive/master.zip
```

or *curl*:

```sh
curl -L https://github.com/abinit/MgB2_eph4isotc/archive/master.zip -o master.zip
```

or simply copy the tarball by clicking the "download button" in the github web page,
unzip the file and rename the directory with:

```sh
unzip master.zip
mv MgB2_eph4isotc-master MgB2_eph4isotc
```

!!! warning

    The directory with the precomputed files **must be located inside *Work_eph4isotc*.
    and must be named `MgB2_eph4isotc`.
    The reason is that all the input files and examples of this tutorial read data from external files
    specified in terms of **relative paths**.


The |AbiPy| script used to perform the GS + DFPT steps is available
[here](https://github.com/abinit/MgB2_eph4isotc/blob/main/run_mgb2_phonons.py).
Note that several parameters have been tuned in order to reach a reasonable **compromise between accuracy
and computational cost** so do not expect the results obtained at the end of the lesson to be fully converged.

We use norm-conserving pseudopotentials with a cutoff energy [[ecut]] of 38 Ha
and the experimental parameters for hexagonal MgB$_2$ (a = 5.8317 and c/a= 1.1416).
All the calculations are performed with a 12x12x12 [[ngkpt]] Gamma-centered $\kk$-grid for electrons (**too coarse**),
and the Marzari smearing ([[occopt]] == 4) with [[tsmear]] = 0.02 Ha.
The DFPT computations is done for 12 irreducible $\qq$-points corresponding
to a $\Gamma$-centered 4x4x4 $\qq$-mesh (again, **too coarse** as we will see in the next sections).
<!--
It is clear that, in more realistic applications, one should monitor the convergence of
lattice parameters and vibrational properties wrt to the $\kk$-mesh and the [[tsmear]] broadening
before embarking on EPH calculations.
This is especially true in metals in which a correct description of the FS is crucial.
In this tutorial, we are trying to find some kind of compromise between accuracy and computational cost.
-->


## Merging the DFPT potentials

To merge the DFPT potential files, copy the first input file with:

```sh
cp $ABI_TESTS//tutorespfn/Input/teph4isotc_1.abi .
```

and execute the *mrgdv* tool using:

```sh
mrgdv < teph4isotc_1.abi
```

The first line in *teph4isotc_1.abi* specifies the name of the output file, followed by the
number of partial DFPT POT files we want to merge and the full list of files:

{% dialog tests/tutorespfn/Input/teph4isotc_1.abi %}

This rather fast step produces the **teph4isotc_1_DVDB** file that will be used in the next examples.
Executing:

```sh
mrgdv info teph4isotc_1_DVDB
```

shows that all the independent atomic perturbations have been computed.

```md
 The list of irreducible perturbations for this q vector is:
    1)  idir= 1, ipert=   1, type=independent, found=Yes
    2)  idir= 2, ipert=   1, type=symmetric, found=No
    3)  idir= 3, ipert=   1, type=independent, found=Yes
    4)  idir= 1, ipert=   2, type=independent, found=Yes
    5)  idir= 2, ipert=   2, type=symmetric, found=No
    6)  idir= 3, ipert=   2, type=independent, found=Yes
    7)  idir= 1, ipert=   3, type=symmetric, found=No
    8)  idir= 2, ipert=   3, type=symmetric, found=No
    9)  idir= 3, ipert=   3, type=symmetric, found=No

 All the independent perturbations are available
 Done
```

## Analyzing electronic and vibrational properties

Before proceeding with the e-ph calculation, it is worth spending some time to analyze
in more detail the electronic and vibrational properties of MgB$_2$.
We will be using the |AbiPy| scripts to post-process the data stored in the precomputed netcdf files.

First of all, let's use |abiopen| to plot the electronic band structure stored in the GSR file produced
by the NSCF calculation (the second task of the first Work i.e. *w0/t1*):

```sh
abiopen.py MgB2_eph4isotc/flow_mgb2_phonons/w0/t1/outdata/out_GSR.nc -e
```

that produces:

![](eph4isotc_assets/MgB2_ebands.png)

The plot shows that there are three bands crossing the Fermi level.
TODO: Discussion about sigma, pi bands and fatbands. Add link to AbiPy fatbands plots

!!! tip

    The band structure is plotted along the high-symmetry $\kk$-path automatically computed by AbiPy.
    To get the explicit lists of high-symmetry $\kk$-point, use the |abistruct| script:

    ```sh
    abistruct.py kpath mgb2_DEN.nc

    ...

    # K-path in reduced coordinates:
     ndivsm 10
     kptopt -11
     kptbounds
        +0.00000  +0.00000  +0.00000  # $\Gamma$
        +0.50000  +0.00000  +0.00000  # M
        +0.33333  +0.33333  +0.00000  # K
        +0.00000  +0.00000  +0.00000  # $\Gamma$
        +0.00000  +0.00000  +0.50000  # A
        +0.50000  +0.00000  +0.50000  # L
        +0.33333  +0.33333  +0.50000  # H
        +0.00000  +0.00000  +0.50000  # A
        +0.50000  +0.00000  +0.50000  # L
        +0.50000  +0.00000  +0.00000  # M
        +0.33333  +0.33333  +0.00000  # K
        +0.33333  +0.33333  +0.50000  # H
    ```

Note that the phonon linewidths have a *geometrical contribution* due to Fermi surface since $\gamma_\qnu$
is expected to be large in correspondence of $\qq$ wave vectors connecting two portions of the FS.
Strictly speaking this is true only if the e-ph matrix elements are constant.
I real materials the amplitude of $g_{mn\nu}(\kk, \qq)$ is not constant
and this may enhance/suppress the value $\gamma_\qnu$ for particular modes.
Yet visualizing the FS is rather useful when discussing e-ph properties in metals.
For this reason, it is useful to have a look at the FS with an external

Graphical tools for the visualization of the FS usually require an external file with
electronic energies in the full BZ whereas ab-initio codes usually take advantage of symmetries
to compute $\ee_\nk$ in the IBZ only.
To produce a BXSF that can be used to visualize the FS with e.g. |xcrysden|,
use the `ebands_bxsf` command of the *abitk* utility located in *src/98_main*
and provide a *GSR.nc* (*WFK.nc*) file with energies computed on a $\kk$-mesh in the IBZ:

```sh
abitk ebands_bxsf MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_GSR.nc
```

This command reconstructs the KS energy in the full BZ by symmetry and produces the `out_GSR.nc_BXSF` file
that can be open with xcrysden using:

```sh
xcrysden --bxsf out_GSR.nc_BXSF
```

Other tools such as |pyprocar| or |fermisurfer| can read data in the BXSF format as well.
For a rather minimalistic matplotlib-based approach, use can also use the |abiview| script with the *fs* command:

```sh
abiview.py fs MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_GSR.nc
```

to visualize the FS in the **unit cell** of the reciprocal lattice.

!!! tip

    The BXSF file can be produced at the end of the GS calculations
    or inside the EPH code by setting [[prtfsurf]] to 1 but *abitk* is quite handy if you already
    have a file and you don't want to write a full Abinit input file and rerun the calculation
    from scratch.


Also the DOS in the region around the Fermi level plays a very important role when discussing
superconducting properties.
Actually this is one of the quantities that should be subject to convergence studies with respect to the $\kk$-grid
before embarking on DFPT calculations.
To compute the DOS with the tetrahedron method, use the `ebands_dos` command of *abitk*

```sh
abitk ebands_dos MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_GSR.nc
```

to read the eigenvalues in the IBZ from a `GSR.n` file (`WFK.nc` files are supported as well).
This command produces a text file named that can be visualized with:

```sh
foobar:
```

We now focus on the vibrational properties of MgB$_2$.
In principle, one can compute phonon frequencies either with *anaddb* or *abinit*.
However, for many applications, it is much easier
to automate the entire process by invoking *anaddb* via the |abiview| script
To compute the phonon band structure using the DDB file produced on the 4x4x4 $\qq$-mesh, use

```sh
abiview.py ddb MgB2_eph4isotc/flow_mgb2_phonons/w1/outdata/out_DDB
```

that produces:

![](eph4isotc_assets/MgB2_phbands.png)

The plot shows that the low-energy vibrations (below 40 meV) are mainly of Mg character.
The vibrational spectrum seems OK (no vibrational instability
the acoustic modes that tend to zero for $|q| \rightarrow 0$ linearly (also because the acoustic sum-rule
is automatically enforced).
Yet this does not mean that our results are fully converged.

In metals, the interatomic force-constants are short-ranged yet this does not guarantee that
a Fourier interpolation done starting from a 4x4x4 $\qq$-mesh can capture all the fine details.
This is clearly seen if we compare the phonon bands computed with a 4x4x4 and a 6x6x6 ab-iniiot $\qq$-mesh.
Also in this case, we can automate the process by using the `ddb` command of the |abicomp| script that takes
in input an arbitrary list of DDB files, calls *anaddb* for all the DDB files and finally compare the results:

```sh
abicomp.py ddb \
   MgB2_eph4isotc/flow_mgb2_phonons/w1/outdata/out_DDB \
   MgB2_eph4isotc/666q_DDB -e
```

![](eph4isotc_assets/abicomp_phbands.png)

The figure reveals that the phonon spectrum interpolated from the 4x4x4 $\qq$-mesh
underestimates the maximum phonon frequency.
Other differences are visible around ~.
Note also that phonon frequencies in metals are also quite sensitive to the $\kk$-mesh and the [[tsmear]] smearing
In real life one should perform an accurate convergence study ...

!!! note

    The *666q_DDB* file was produced with the same AbiPy script by just changing the value of *ngqpt*.
    The reason why we don't provide files obtained with a 6x6x6 sampling is that the size of the
    git repository including all the DFPT potentials would be around 200 Mb.

    Let's continue using the 4x4x4 DDB file but we should take this difference into account when comparing
    our results with previous works.

## Our first computation of the isotropic Tc

For our first example, we use a relatively simple input file that allows us to introduce
the most important variables and the organization of the main output file:
Copy the input file in the working directory and execute it using:

```sh
abinit teph4isotc_2.abi > log 2> err
```

If you prefer, you can run it in parallel with e.g two MPI processes using:

```sh
mpirun -n 2 abinit teph4isotc_2.abi > log 2> err
```

without having to introduce any input variable for the MPI parallelization
as the EPH code can automatically distribute the workload.
Further information concerning the MPI version are given in the last part of the tutorial.

We now discuss the meaning of the different variables in more detail.

{% dialog tests/tutorespfn/Input/teph4isotc_2.abi %}

To activate the computation of $\gamma_{\qq\nu}$ in metals, we use [[optdriver]] = 7 and [[eph_task]] = 1.
The location of the DDB, DVDB and WFK files is specified via
[[getddb_filepath]] [[getdvdb_filepath]] [[getwfk_filepath]], respectively.

```sh
getddb_filepath  "MgB2_eph4isotc/flow_mgb2_phonons/w1/outdata/out_DDB"

getwfk_filepath  "MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_WFK.nc"

getdvdb_filepath "teph4isotc_1_DVDB"
```

The DDB and the WFK files are read from the git repository while for the DVDB we use
the file produced by *mrgdv* in the previous section:
Note also the use of the new input variable [[structure]] (added in Abinit v9)
with the *abifile* prefix to read the crystalline structure from an external file

```sh
structure "abifile:MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_DEN.nc"
```

so that we do can avoid repeating the unit cell in every input file.

Next, we have the variables governing the EPH calculation:

```sh
ddb_ngqpt 4 4 4
dipdip 0         # No treatment of the dipole-dipole part. OK for metals

eph_intmeth 1
eph_fsmear -1.0
eph_fsewin 0.3 eV
```

In this first example, we prefer not to interpolate the DFPT potentials so only [[ddb_ngqpt]] is specified.
[[dipdip]] is set to zero as we are dealing with a metal.
The $\kk$-point integration is done with the adaptive Gaussian smearing
[[eph_intmeth]], [[eph_fsewin]] 0.3 eV [[eph_intmeth]] 1 [[eph_fsmear]] -1

Then we activate two tricks just to make the calculation faster:

```sh
mixprec 1
boxcutmin 1.1
```

!!! tip

    As explained in the documentation, using [[mixprec]] = 1 and [[boxcutmin]] = 1.1
    should not have significant effects on the final results yet these are not the default values
    as users are supposed to compare the results with/without this trick before running production calculations.

<!--

TODO: discuss ph_nqpath and ph_qpath
In this particular example, there's no difference between a calculation done with *dipdip* 0 or 1 as our
DDB does not contain $\ee^\infty$ and the Born-effective charges.
Finally, we have the specification of the $\kk$-mesh must in terms of [[ngkpt]], [[nshiftk]] and [[shiftk]]
that must be equal to the one used to generate the input WFK file:

```sh
ngkpt   12 12 12
nshiftk 1
shiftk  0.0 0.0 0.0
```

In this section, we use the interface provided by the EPH code to
compute vibrational properties using the pre-generated DDB file instead of the *anaddb* tool.
To compute vibrational properties without performing a full e-ph calculation,
use [[optdriver]] = 7 with [[eph_task]] = 0.
A typical template looks like:

```sh
optdriver 7
eph_task 0

# DDB file
getddb_filepath "test_submodule/mgb2_121212k_0.01tsmear_DDB"
ddb_ngqpt 8 8 8

dipdip 0  # No treatment of the dipole-dipole part. OK for metals

# phonon DOS with q-mesh (default tetrahedron method)
ph_ngqpt
#ph_intmeth

# phonon bands along high-symmetry q-path
ph_nqpath
ph_qpath
ph_ndivsm
```

[[getddb_filepath]] specifies the path to the external DDB file while [[ddb_ngqpt]] defines the
*ab-initio* $\qq$-mesh used to generate the DDB file.
The phonon DOS will be computed using the (dense) [[ph_ngqpt]] $\qq$-mesh and [[ph_intmeth]].
The high-symmetry $\qq$-path for the phonon band structure is specified with:
[[ph_nqpath]], [[ph_qpath]], and [[ph_ndivsm]].

!!! important

    [[dipdip]] can be set to zero as we are dealing with a metal and therefore the IFCs are short-ranged.
    Note that the default value of [[dipdip]] is designed for polar semiconductors so we recommend
    to override the default behaviour when performing calculations with [[eph_task]] = 1.

Since Abinit supports multidatasets, unlike anaddb, it is easy to define an input file to compute
the phonon DOS with multiple $\qq$-meshes.
This simple test allows us to get an initial (very qualitative) estimate of the $\qq$-sampling
required to convergence the  Eliashberg function as $\alpha^2F(\ww)$ is essentially a weighted phonon DOS.

The input file xxx, shows how to perform such a test with $\qq$-meshes of increasing density.

To analyse the results one can extract the data from the ... PHDOS.nc files

This is what you should get:

!!! important

    Multiple datasets may be handy when running small calculations as in this case.
    Remember, however, that the EPH code is not designed to be used with multiple datasets.
    In a nutshell, try to split your input files as much as possible so that they can be executed in parallel
    with different number of CPUS and different amount of memory.

Remember to discuss k-mesh and [[tsmear]] at the DFPT level.

getden_filepath "MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_DEN.nc"
-->

We now discuss in more detail the output file produced by the code.
After the standard section with info on the unit cell and the pseudopotentials used,
we find the output of the electronic DOS:

```sh
 Linear tetrahedron method.
 Mesh step:  10.0 (meV) with npts: 2757
 From emin:  -5.1 to emax:  22.5 (eV)
 Number of k-points in the IBZ: 133
 Fermi level:   7.63392108E+00 (eV)
 Total electron DOS at Fermi level in states/eV:   7.36154553E-01
 Total number of electrons at eF:    8.0

- Writing electron DOS to file: teph4isotc_2o_DS1_EDOS
```

Then the code outputs some basic info concerning the Fermi surface
and the integration method for the double delta:

```sh
 ==== Fermi surface info ====
 FS integration done with adaptive gaussian method
 Total number of k-points in the full mesh: 1728
 For spin: 1
    Number of BZ k-points close to the Fermi surface: 291 [ 16.8 %]
    Maximum number of bands crossing the Fermi level: 3
    min band: 3
    Max band: 5
```

Then, for each $\qq$-point in the IBZ, the code outputs the values of $\ww_\qnu$, $\gamma_\qnu$ and $\lambda_\qnu$:

```md
 q-point =    0.000000E+00    0.000000E+00    0.000000E+00
 Mode number    Frequency (Ha)  Linewidth (Ha)  Lambda(q,n)
    1        0.000000E+00    0.000000E+00    0.000000E+00
    2        0.000000E+00    0.000000E+00    0.000000E+00
    3        0.000000E+00    0.000000E+00    0.000000E+00
    4        1.444402E-03    4.898804E-07    3.731164E-03
    5        1.444402E-03    4.898804E-07    3.731164E-03
    6        1.674636E-03    2.955247E-15    1.674492E-11
    7        2.426616E-03    1.025439E-03    2.767185E+00
    8        2.426616E-03    1.025439E-03    2.767185E+00
    9        3.127043E-03    3.508418E-08    5.701304E-05
```

The linewidths of the acoustic mode at $\Gamma$ are zero.
The default value of [[asr]] is automatically set to 1

```sh
abiview.py ddb_asr MgB2_eph4isotc/flow_mgb2_phonons/w1/outdata/out_DDB
```

Finally, we have the total value of the isotropic $\lambda$:

```sh
 lambda=   0.4286
```

The calculation has produced the following output files:

```sh
$ ls teph4isotc_2o_DS1*

teph4isotc_2o_DS1_A2F.nc         teph4isotc_2o_DS1_NOINTP_A2FW    teph4isotc_2o_DS1_PHDOS.nc
teph4isotc_2o_DS1_A2FW           teph4isotc_2o_DS1_NOINTP_PH_A2FW teph4isotc_2o_DS1_PHGAMMA
teph4isotc_2o_DS1_EBANDS.agr     teph4isotc_2o_DS1_PHBANDS.agr    teph4isotc_2o_DS1_PH_A2FW
teph4isotc_2o_DS1_EDOS           teph4isotc_2o_DS1_PHBST.nc
```

The *A2F.nc* netcdf file stores all the results of the calculation in a format that can
can be visualized with |abiopen|:

```sh
abiopen.py teph4isotc_2o_DS1_A2F.nc -e --seaborn
```

![](eph4isotc_assets/skw_vs_abinitio_ebands.png)

Exercise:

 Since our WFK is defined on a 12x12x12 $\kk$-mesh it is very easy to activate the interpolation
 of the DFPT potentials to go from a 6x6x6 to a 12x12x12 $\qq$-mesh while keeping the same WFK file.
 Use [[eph_ngqpt_fine]] in the previous input file to densify the $\qq$-mesh for phonons and compare the
 results with those obtained with a 6x6x6 sampling.
 Use also [[mixprec]] 1 and [[boxcutmin]] 1.1 to accelerate the calculation and compare the results.
 You may want to run the calculation in parallel with mpirun.

## Preparing the convergence study wrt the k-mesh

<!--
Remember that in the EPH code the $\qq$-mesh can be changed at will thanks to the Fourier interpolation
of the dynamical matrix and of the DFPT potentials whereas the $\kk$-mesh must correspond
to the one used to generate the input WKF file.

Before running calculations with difference meshes, it is a good idea
to sketch the convergence study we want to perform.
Let's assume we want to test the following configurations:

A 6x6x6 $\qq$-mesh (without interpolation) and the following $\kk$-meshes:

    12x12x12, 18x18x18, 24x24x24 (36x36x36)

A 12x12x12 $\qq$-mesh (with interpolation) and the following $\kk$-meshes:

    12x12x12, 24x24x24, (36x36x36)

The pattern at this point should clear: we need to perform NSCF computations
for four different Gamma-centered $\kk$-meshes: 12x12x12, 18x18x18, 24x24x24, and (48x48x48).
The computational cost of the NSC run increases quickly with the size of the $\kk$-mesh but as already
stressed we only needed electron states in a appropriate energy window around the Fermi level.

To speed the computation of the WFK files we will be using the same trick based
on the star-function interpolation used in the mobility tutorial.
-->

Our goal is to perform calculations of $\gamma_\qnu$ and $\lambda_\qnu$ with different $\qq/\kk$-meshes
in order to analyze the convergence behaviour of $\alpha^2F(\ww)$ and $\lambda$.
For example ...

As the $\kk$-mesh must be a multiple of the $\qq$-mesh, we need to generate **different WFK files**
to prepare our convergence studies.
This preliminary step becomes quite CPU-consuming if very dense $\kk$-meshes are needed.
Fortunately, we can speed up this part since the computation of $\gamma_\qnu$, due to the double delta,
requires the knowledge of Bloch states only inside a relatively small energy window around $\ee_F$.

Similarly to what is done in the [eph4mob tutorial](eph4mob), we can take advantage of the star-function interpolation to
find the $\kk$ wavevectors whose energy is inside the [[sigma_erange]] energy window **around the Fermi level**
and then perform a NSCF calculations restricted to the relevant wavevectors.
This is the approach we will be using in this tutorial.

First of all, we recommend to test whether the SKW interpolation can reproduce the ab-initio results
with reasonable accuracy.
A possible approach consists in comparing the ab-initio band structure with the interpolated one
using the `skw_compare` command of *abitk*:

```sh
abitk skw_compare \
    MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_GSR.nc \
    MgB2_eph4isotc/flow_mgb2_phonons/w0/t1/outdata/out_GSR.nc  --is-metal
```

The fist argument is a GSR.nc (WFK.nc) file with the energies in the IBZ.
used to build the star-function interpolant.
The quality of the interpolant will obviously depend on the density of this $\kk$-mesh.
The second file contains the ab-initio eigenvalues along a $\kk$-path computed with a standard NSCF run.

!!! Note

    Note the use of the `--is-metal` option else *abitk* will complain that it cannot compute the gaps!

*skw_compare* produces two netcdf files (*abinitio_EBANDS.nc* and *skw_EBANDS.nc*) that we can compared  
using the |abicomp| script:

```sh
abicomp.py ebands abinitio_EBANDS.nc skw_EBANDS.nc -p combiplot
```

that produces:

![](eph4isotc_assets/skw_vs_abinitio_ebands.png)

The figure shows that the SKW interpolant nicely reproduces the ab-initio dispersion around the Fermi level.
Discrepancies between the ab-initio results and the interpolated values are visible
in correspondence of **band-crossings**.
This is expected since SKW is a Fourier-based approach and band-crossing leads to a non-analytic behaviour
of the signal that cannot be faithfully reproduced with a finite number of Fourier components.
Fortunately, these crossings are relatively far from the Fermi level and do not enter into play
in the computation of superconducting properties.

!!! tip

    Should the fit be problematic, use the command line options ... to improve the results
    and take note of the SKW parameters
    as the same values should be used when calling Abinit with [[einterp]].
    The worst case scenario of band crossing close to the Fermi level may be handled by enlarging the
    [[sigma_erange]] energy window.


At this point, we are confident that the SKW interpolation is OK and we can use it to locate
the $\kk$-wavevectors of a much denser $\kk$-mesh whose energy is inside the ([[sigma_erange]]) window
around the Fermi level.
This is what is done in *teph4isotc_4.abi*:

{% dialog tests/tutorespfn/Input/teph4isotc_3.abi %}

The most important section of the input file is reported below:

```sh
optdriver 8
wfk_task "wfk_kpts_erange"
getwfk_filepath  "MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_WFK.nc"

# Define fine k-mesh for the SKW interpolation
sigma_ngkpt   24 24 24
sigma_nshiftk 1
sigma_shiftk  0 0 0

sigma_erange -0.3 -0.3 eV  # Select kpts in the fine mesh within this energy window.
einterp 1 5 0 0            # Parameters for star-function interpolation (default values)
```

The first section ([[optdriver]] and [[wfk_task]]) activates the computation of the KERANGE file
with ab-initio energies taken from [[getwfk_filepath]].
The three variables [[sigma_ngkpt]], [[sigma_nshiftk]] and [[sigma_shiftk]] define the final dense $\kk$-mesh
In this case, we are using the default values of [[einterp]].
You may need to change these parameters

!!! important

    When dealing with metals, both entries in [[sigma_erange]] must be **negative** so that the energy window is
    refered to the Fermi level and not to the CBM/VBM.
    Use a value for the window that is reasonably large in order to account for possible oscillations
    and/or inaccuracies of the SKW interpolation.

that produces the *KERANGE.nc* file.
Once we have the *KERANGE.nc* file, we can use it to perform a NSCF calculation to generate
a customized WFK file on the dense $\kk$-mesh to save significant computing time and space of disk.

Finally, we can use these files to perform a NSCF calculation with *teph4isotc_4.abi* input:
with the [[getkerange_filepath]]

{% dialog tests/tutorespfn/Input/teph4isotc_4.abi %}

to perform a NSCF calculation with [[kptopt]] 0 to produce a new WFK file on the dense $\kk$-mesh.

```sh
iscf  -2
tolwfr 1e-18
kptopt 0                        # Important!
```

!!! note

    Note how we use [[getden_filepath]] and the syntax:

    ```sh
    getden_filepath "MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_DEN.nc"
    ```

    to start the NSCF run from the **precomputed** GS DEN file.

## Convergence study wrt to the k/q-mesh

At this point, we can run EPH calculations with denser $\kk$-meshes as discussed in the next section.
We will be using settings similar to the ones used in **teph4isotc_4.abi** except for
the use of [[eph_ngqpt_fine]] and [[ngkpt]]

{% dialog tests/tutorespfn/Input/teph4isotc_5.abi %}

## A more accurate converged calculation

## Notes on the MPI parallelism

EPH calculations done with [[eph_task]] = 1 support 5 different levels of MPI parallelism
and the number of MPI processes for each level can be specified via [[eph_np_pqbks]].
This variable is optional in the sense that whatever number of MPI processes you use, EPH will try to select
a reasonable distribution of the workload.
The distribution, however, may not be optimal as EPH tries to minimize memory requirements by focusing on the
perturbation/k-point parallelism.

As usual, MPI algorithms are quite efficient if the distribution of the workload is done at a very high-level.
In the case of $T_c$ calculations, the outermost loop is over the $\qq$-points in the IBZ hence the highest speedup
is achieved when most of the CPUs are used for the $\qq$-point parallelism.
Note, however, that this kind of MPI distribution does not distribute the wavefunctions and the scattering potentials.
