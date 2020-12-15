---
authors: MG
---

# Superconducting properties within the isotropic Eliashberg formalism

This tutorial explains how to compute electron-induced phonon linewidths in metals
and how to use the McMillan equation to estimate the superconducting critical temperature $T_c$ within
the isotropic Eliashberg formalism.
We start by presenting the basic equations implemented in the code and the connection with the ABINIT variables.
Then we discuss how to run isotropic $T_c$-calculations and how to perform typical convergence studies
using hexagonal $MgB_2$ as example.

It is assumed the user has already completed the two tutorials [RF1](rf1) and [RF2](rf2),
and that he/she is familiar with the calculation of ground state and vibrational properties in metals.
The user should have read the [fourth lesson on Al](base4) as well
as the [introduction page for the EPH code](eph_intro) before running these examples.

This lesson should take about 1.5 hour.

## Formalism and connection with the implementation

Due to the interaction with electrons, phonons acquire a finite lifetime
given by the imaginary part of the phonon-electron self-energy $\Pi_\qnu$.
In the so-called double-delta approximation, the phonon linewidth $\gamma_{\qq\nu}$ due
to the interaction of the ${\qq\nu}$ phonon with electrons is given by

$$
\gamma_{\qq\nu} = 2\pi \ww_{\qq\nu} \sum_{mn\kk} |g_{mn\nu}(\kk, \qq)|^2
\delta(\ee_{\kpq m} -\ee_F) \delta(\ee_{\kk n} -\ee_F)
$$

where $\ww_{\qq\nu}$ is the phonon frequency,
the sum over the electron wavevector $\kk$ runs over the full BZ, $\ee_F$ is the Fermi level
and $g_{mn\nu}(\kk, \qq)$ are the e-ph matrix elements discussed in the [EPH introduction](eph_intro).
For a given phonon wavevector $\qq$, the double delta restricts the BZ integration to
transitions between $\kk$ and $\kq$ electronic states on the Fermi surface (FS).
Converging the double-delta integral therefore requires very dense $\kk$-meshes in order
to sample enough states around the FS.
The convergence rate is indeed much slower than the one required by the electron DOS at $\ee_F$:

$$
g(\ee_F) = \sum_{n\kk} \delta(\ee_F - \ee_{n\kk})
$$

in which a single Dirac delta is involved.

!!! note

  Note that there is also another important contribution to the phonon lifetimes due to non-harmonic
  terms in the expansion of the Born-Oppenheimer energy surface around the equilibrium point.
  In the many-body language these non-harmonic leads to phonon-phonon scattering that can give a substantial
  contribution to the phonon linewidths.
  In the rest of the tutorial, however, non-harmonic terms will be ignored and we will be mainly focusing
  on the computation of the imaginary part of $\Pi$ in the harmonic approximation.

At the level of the implementation, the integration of the double delta can be performed
either with the tetrahedron scheme or by replacing the Dirac delta with a Gaussian of finite width.
The integration algorithm is defined by the [[eph_intmeth]] input variable.
In the case of Gaussian method, one can use a fixed broadening ([[eph_fsmear]] > 0)
or an adaptive scheme (activated by using a negative value for [[eph_fsmear]]) in which the broadening
is automatically computed from the electron group velocities [[cite:Li2015]].

!!! important

    The tetrahedron method is more accurate and does not require any broadening parameter.
    Note, however, that in the present implementation the computational cost of the double
    delta with the tetrahedron method
    quickly increases with the size of the $\kk$-mesh so the adaptive Gaussian scheme represents a valid alternative,
    especially when a dense $\kk$-sampling is used.

The value of the Fermi level, $\ee_F$, is automatically computed from the KS eigenvalues stored
in the input WFK file according to the two input variables [[occopt]] and [[tsmear]].
These parameters are usually equal to the ones used for the GS/DFPT calculation.
However, it is possible to change the value of $\ee_F$ at the EPH level using three (mutually exclusive) input variables:
[[eph_fermie]], [[eph_extrael]] and [[eph_doping]].

The summation can be limited to states within an energy window around the Fermi level specified via [[eph_fsewin]].

The $\kk$-mesh used by the EPH code is defined by the input variables [[ngkpt]], [[nshiftk]] and [[shiftk]]
Note, that the $\kk$-sampling  must correspond to the mesh used to generated the input WFK file.
Convergence studies are performed by generating different WFK files on $\kk$-meshes of increasing density.

The code computes $\gamma_{\qq\nu}$ for each $\qq$-point in the IBZ associated to
an arbitrary $\qq$-mesh specified by the user.
By default, the EPH code uses the [[ddb_ngqpt]] $\qq$-mesh corresponding to the DDB file
(that is assumed to be equal to the one used to generate the DVDB file).
In this case, all the DFPT scattering potentials are available and no interpolation in $qq$-space in required.
To increase the $qq$-sampling, one simply specifies [[eph_ngqpt_fine]] in the input file.
In this case, the code employs the Fourier interpolation to obtain the scattering potentials.

Once the phonon linewidths $\gamma_{\qq\nu}$ are known in the IBZ, EPH computes the Eliashberg function defined by:

$$
\alpha^2F(\ww) = -\dfrac{1}{N_F} \sum_{\qq\nu} \dfrac{\gamma_{\qq\nu}}{\ww_{\qq\nu}} \delta(\ww - \ww_{\qq \nu})
$$

where $N_F$ is the density of states (DOS) per spin at the Fermi level.
$\alpha^2F(\ww)$ gives the strength by which a phonon of energy $\ww$ scatters electronic
states on the FS (remember that Abinit uses atomic units hence $\hbar = 1$).
This quantity is accessible in experiments and experience has shown that
$\alpha^2F(\ww)$ is qualitatively similar to the phonon DOS $F(\ww)$:

$$
F(\ww) = \sum_{\qq\nu} \delta(\ww - \ww_{\qq \nu})
$$

This is not surprising as the equation for $\alpha^2F(\ww)$ resembles the one for the phonon DOS $F(\ww)$:
except for the weighting factor $\frac{\gamma_{\qq\nu}}{\ww_{\qq\nu}}$.
The technique used to compute $\alpha^2F(\ww)$ is defined by the two variables [[ph_intmeth]] and [[ph_smear]].
By default, the code uses the tetrahedron method for the $\qq$-space integration.

The total e-ph coupling strength $\lambda$ is defined as the first inverse moment of $\alpha^2F(\ww)$:

$$
\lambda = \int \dfrac{\alpha^2F(\ww)}{\ww}\dd\ww = \sum_{\qq\nu} \lambda_{\qq\nu}
$$

where we have introduced the mode dependent coupling strength:
<!-- For spin unpolarized systems: -->

$$
\lambda_{\qq\nu} = \dfrac{\gamma_{\qq\nu}}{\pi N_F \ww_{\qq\nu}^2}
$$

Finally, the isotropic superconducting temperature $T_c$ can be estimated using the McMillan equation:

$$
T_c = \dfrac{\ww_{log}}{1.2} \exp \Biggl [ \dfrac{-1.04 (1 + \lambda)}{\lambda ( 1 - 0.62 \mu^*) - \mu^*} \Biggr ]
$$

where $\mu^*$ is a semi-empirical parameter that describes the (screened) e-e interaction while
$\ww_{\text{log}}$ is the *logarithmic* average of the phonon frequencies defined by:

$$
\ww_{\text{log}} = \exp \Biggl [ \dfrac{2}{\lambda} \int \dfrac{\alpha^2F(\ww)}{\ww}\log(\ww)\dd\ww \Biggr ]
$$

!!! important

    The code computes $\gamma_{\qq\nu}$ for all the $\qq$-points in the IBZ associated
    to a homogeneous $\qq$-mesh as these quantities are then used to evalute integrals in $\qq$-space.
    [[ddb_ngqpt]] mesh that is the $\qq$-mesh used in the DFPT computations.
    In this case, the code does not performy any kind of interpolation in $\qq$-space.
    A (much denser) $\qq$-mesh can be specified via [[eph_ngqpt_fine]].
    In this case, the code interpolates the DFPT potentials on the dense mesh at runtime.
    EPH computes the electron DOS using the eigenvalues stored in the WFK file.

<!--
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

    The directory with the precomputed files must be located in the same working directory
    in which you will be executing the tutorial and must be named `MgB2_eph4isotc`.


The |AbiPy| script used to executed the DFPT part is available
[here](https://github.com/abinit/MgB2_eph4isotc/blob/main/run_mgb2_phonons.py).
Note that several parameters have been tuned to reach a reasonable **compromise between accuracy
and computational cost** so do not expect the results obtained at the end of the lesson to be fully converged.
More specifically, we use norm-conserving pseudopotentials with a cutoff energy [[ecut]] of 38 Ha
We use the experimental parameters for hexagonal $MgB_2$ (a = 5.8317 and c/a= 1.1416)
All the calculations have been performed with a 12x12x12 [[ngkpt]] Gamma-centered grid for electrons,
and the Marzari smearing [[occopt]] with [[tsmear]].
The DFPT computations is done for XXX irreducible $\qq$-points corresponding
to a $\Gamma$-centered 4x4x4 $\qq$-mesh (again, too coarse).

<!--
The input file of the GS run is also stored in the DEN.nc file and one can easily access it with the
*ncdump* utility

It is clear that, in more realistic applications, one should monitor the convergence of
lattice parameters and vibrational properties wrt to the $\kk$-mesh and the [[tsmear]] broadening
before embarking on EPH calculations.
This is especially true in metals in which a correct description of the FS is crucial.
In this tutorial, we are trying to find some kind of compromise between accuracy and computational cost.
-->


## Merging the DFPT potentials

To merge the DFPT potential files, execute the *mrgdv* tool using:

```sh
mrgdv < teph4isotc_1.abi
```

with the following input file that lists all the partial DFPT POT files already computed:

{% dialog tests/tutorespfn/Input/teph4isotc_1.abi %}

This (rather fast) step produces the **teph4isotc_1_DVDB** file that will be used in the next examples.
Executing:

```sh
mrgdv info teph4isotc_1_DVDB
```

shows that all the independent phonon perturbations are available

<!--
As mentioned in the [introduction page for the EPH code](eph_intro), the DVDB file is needed to reconstruct
the scattering potentials in the full BZ for all the 3 x [[natom]] atomic perturbations.
With this file, one can obtain $\Delta_\qnu V^\KS$ for all the $\qq$-points belonging to the initial 6x6x6 DFPT
$\qq$-mesh or even a much denser $\qq$-sampling when the Fourier interpolation of the potentials is exploited.
Still the computation of the e-ph matrix elements in the approach implemented in EPH requires the explicit
knowledge of Bloch states at $\kk$ and $\kq$.
This is the problem we will try to address in the next section.

At this point, we have all the ingredients required to compute phonon linewidths and $T_c$.
As mentioned in the introductions $\gamma_{\qq\nu}$ converge slowly wrt to the $\kk$-sampling
so we need to monitor the convergence wrt the number of wavevectors for electrons.

We will try to get an initial rough estimate by looking at the BZ sampling required to
properly describe quantities that are relatively easy to compute:

1. The electronic DOS at the Fermi level
2. The phonon DOS.

The first test gives us an idea of the $\kk$-mesh whereas the second test gives us an idea
of the $\qq$-mesh.

Let's start from the electron DOS.
-->

## Analyzing electronic and vibrational properties

Before proceeding with the e-ph calculation, it is worth spending some time to analyze 
in more detail the electron and phonon band structures computed by AbiPy.


```sh
abiopen.py MgB2_eph4isotc/flow_mgb2_phonons/w0/t1/outdata/out_GSR.nc -e
```

![](eph4isotc_assets/MgB2_ebands.png)


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

To compute the electron DOS with the tetrahedron method,
one can use the *abitk* utility located in *src/98_main* and the syntax:

```sh
abitk ebands_dos MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_GSR.nc 
```

to read the eigenvalues in the IBZ from a `GSR.n` file (`WFK.nc` files are supported as well).
We can also produce a BXSF file with the `ebands_bxsf` command:

```sh
abitk ebands_bxsf MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_GSR.nc 
```

to produce a file `out_GSR.nc_BXSF` that can be used to visualize the Fermi surface with Xcrysden e.g.

```
xcrysden --bxsf out_GSR.nc_BXSF
```

!!! tip

    The same feature is available inside Abinit if [[prtfsurf]] is set to 1.

<!--
```
abinit
```

{% dialog tests/tutorespfn/Input/teph4isotc_2.abi %}

All the NSCF calculations will start from the **pre-computed DEN.nc file** via the [[getden_filepath]] input variable.


!!! important

    Note how we use [[getwfk]] = -1 to read the WFK of the previous dataset to accelerate the calculation.
    Abinit, indeed, can initialize the wavefunctions from a previous WFK file even if the $\kk$-mesh is different.
    In our example this trick is beneficial as we are using a single input with multiple datasets.
    Keep in mind, however, that this kind of algorithm is intrinsically sequential in the sense that
    you need to complete dataset $n$ before moving to step $n+1$.
    If enough computing capabilities are avaiable, it is more efficient to perform NSCF calculations with different
    $\kk$-grids indipendently without using [[ndtset]] > 1.
    Consider that the workload is proportional to [[nkpt]] so running multiple datasets with the same number of
    MPI processes is not necessarily the most efficient way.
    Obviously you are free to use the multidataset philosophy whenever it makes calculations easier to handle
    but keep in mind that the EPH code is not designed with this idea in mind and that you will start
    to experience the multidaset inefficiency when the dimension of the problem increases.
-->

## Vibrational properties

Let's call anaddb via anaddb to compute the phonon band structure using a 4x4x4

```sh
abiview.py ddb MgB2_eph4isotc/flow_mgb2_phonons/w1/outdata/out_DDB
```

![](eph4isotc_assets/MgB2_phbands.png)

The fact that 4x4x4 is too coarse is clearly seen if we compare
The reason why we don't provide files obtained with a 6x6x6 is that the size of repo with the DFPT potentials 
would be around 200 Mb

```sh
abicomp.py ddb MgB2_eph4isotc/flow_mgb2_phonons/w1/outdata/out_DDB MgB2_eph4isotc/666q_DDB -e
```

![](eph4isotc_assets/abicomp_phbands.png)


<!--
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
    Note that the default value of [[dipdip]] is designed for polar semiconductors so we recommended
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

Remember to discuss k-mesh and tsmear at the DFPT level.
-->

## Our first computation of the isotropic Tc

For our first example, we use a relatively simple input file that allows us to introduce
the most important variables and the organization of the main output files:

{% dialog tests/tutorespfn/Input/teph4isotc_2.abi %}

To activate the computation of $\gamma_{\qq\nu}$ in metals, we use [[optdriver]] = 7 and [[eph_task]] = 1.
The location of the DDB, DVDB and WFK files is specified via
[[getddb_filepath]] [[getdvdb_filepath]] [[getwfk_filepath]], respectively:

```sh
getddb_filepath  "MgB2_eph4isotc/flow_mgb2_phonons/w1/outdata/out_DDB"

getwfk_filepath  "MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_WFK.nc"

getdvdb_filepath "teph4isotc_1_DVDB"
```

The DDB and the WFK are taken from the git repository while for the DVDB we specify 
the file produced by *mrgdv* in the previous step:
Note the usage of the new input variable [[structure]] (added in Abinit v9) to read the crystalline structure from
an external file in order to avoid repeating the unit cell in each input file.

```sh
structure "abifile:MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_DEN.nc"

getden_filepath "MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_DEN.nc"
```

Next, we have the variables governing the EPH calculation:

[[ddb_ngpqt]], [[dipdip]]
[[eph_intmeth]]
[[eph_fsewin]] 0.3 eV
[[eph_intmeth]] 1
[[eph_fsmear]] -1

``
ddb_ngqpt 4 4 4
dipdip 0         # No treatment of the dipole-dipole part. OK for metals

eph_intmeth 1
eph_fsmear -1.0
eph_fsewin 0.3 eV
```

In this example, we are using the adaptive Gaussian smearing

Tricks to make EPH calculations faster.
[[mixprec]] 1 and [[boxcutmin]] 1.1

```
mixprec 1
boxcutmin 1.1
```

Finally, the k-mesh must in terms of [[ngkpt]], [[nshift]] and [[shiftk]] 
that must be consistent with the sampling used to generate the WFK file:

```
ngkpt   12 12 12
nshiftk 1
shiftk  0.0 0.0 0.0
```

We now discuss in more detail the output files produced by the code.
TODO: Discuss output, netcdf files and provide AbiPy examples.

Exercise:
 
 Since our WFK is defined on a 12x12x12 $\kk$-mesh it is very easy to activate the interpolation
 of the DFPT potentials to go from a 6x6x6 to a 12x12x12 $\qq$-mesh while keeping the same WFK file.
 Use [[eph_ngqpt_fine]] in the previous input file to densify the $\qq$-mesh for phonons and compare the 
 results with those obtained with 6x6x6.
 Use also mixprec 1 and boxcutmin 1.1 to accelerate the calculation and compare the results.

## Preparing the convergence study wrt the k-mesh

Our goal is to perform calculations of $\gamma_{\qq\nu}$ with different $\qq/\kk$-meshes to
analyze the convergence behaviour of $\lambda$ and $\alpha^2F(\ww)$.
Remember that in the EPH code the $\qq$-mesh can be changed at will thanks to the Fourier interpolation
of the dynamical matrix and of the DFPT potentials whereas the $\kk$-mesh must correspond
to the one used to generate the input WKF file.
As the $\kk$-mesh must be a multiple of the $\qq$-mesh, we need to generate different WFK files
in order to perform our convergence studies.

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
To speed the computation of the WFK files we will be using the same trick based on the star-function interpolation
used in the Mobility tutorial.

TODO
How to reduce the number of k-points to be computed in the NSCF run
As we have seen, Eliashberg calculations require the knowledge of Bloch states
inside a relatively small energy window around $\ee_F$.
The NSCF computation of the WFK files with dense $\kk$-sampling [[sigma_erange]]

```
abitk skw_compare MgB2_eph4isotc/flow_mgb2_phonons/w0/t0/outdata/out_GSR.nc MgB2_eph4isotc/flow_mgb2_phonons/w0/t1/outdata/out_GSR.nc  --is-metal
```

```
abicomp.py ebands abinitio_EBANDS.nc skw_EBANDS.nc -p combiplot
```

![](eph4isotc_assets/skw_vs_abinitio_ebands.png)

The figure shows that the SKW interpolation nicely reproduces the band dispersion around the Fermi level.
Discrepancies between the ab-initio energies and the interpolated ones are observed 
in correspondence of band-crossings.
This is somehow expected since SKW is a Fourier-based interpolation scheme but in our fortunately the crossings
are relatively far from the Fermi level.
At this point, we are confident that the SKW interpolation works as expected and we can use it to locate
the $\kk$-wavevectors of a much denser $\kk$-mesh that are located inside an energy window [[sigma_erange]]
to generate the *KERANGE.nc* file

```
optdriver 8
wfk_task "wfk_kpts_erange"
getwfk_filepath "161616_WFK"

# Define fine k-mesh for the SKW interpolation
sigma_ngkpt   32 32 32
sigma_nshiftk 1
sigma_shiftk  0 0 0

sigma_erange -0.2 -0.2 eV  # Select kpts in the fine mesh within this energy window.
einterp 1 5 0 0            # Parameters for star-function interpolation
```

!!! important

    When dealing with metals, sigma_erange must be negative so that the energy window is
    refered to the Fermi level and not to the CBM/VBM.
    Use a value for the window that is reasonably large in order to account for possible oscillations
    and/or inaccuracies of the SKW interpolation.

Once we have the KERANGE.nc file, we can use it to perform a NSCF calculation to generate a customized WFK file
on the dense $\kk$-mesh.
This is what is done in ...


{% dialog tests/tutorespfn/Input/teph4isotc_3.abi %}

that produces the following (customized) WFK files:

- A
- B

At this point, we can run EPH calculations with denser $\kk$- and $\qq$-meshes.

## Convergence study wrt to the k/q-mesh

[[eph_ngqpt_fine]]

{% dialog tests/tutorespfn/Input/teph4isotc_3.abi %}

## Notes on the MPI parallelism

The EPH code supports 5 different levels of MPI parallelism and the number of MPI processes for each level
can be specified via [[eph_np_pqbks]].
This variable is optional in the sense that whatever number of MPI processes you use, EPH will try to select
a reasonable distribution of the workload.
The distribution, however, may not be optimal as EPH tries to minimize memory requirements by focusing on the
perturbation/k-point parallelism.
As usual, MPI algorithms are quite efficient if the distribution of the workload is done at a very high-level.
In the case of $T_c$ calculations, the outermost loop is over the $\qq$-points in the IBZ hence the highest speedup
is achieved when most of the CPUs are the used for the $\qq$-point parallelism.
Note, however, that this kind of MPI distribution does not distribute the wavefunctions and the scattering potentials.
