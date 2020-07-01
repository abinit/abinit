---
authors: MG
---

# Superconducting properties within the isotropic Eliashberg formalism

This tutorial explains how to compute electron-induced phonon linewidths in metals
and how to use the McMillan equation to estimate the superconducting critical temperature $T_c$ within
the isotropic Eliashberg formalism.
We start by presenting the equations implemented in the code and their connection with the ABINIT input variables.
Then we discuss how to run isotropic $T_c$-calculations and how to perform typical convergence studies
for hexagonal $MgB_2$.

It is assumed the user has already completed the two tutorials [RF1](rf1) and [RF2](rf2),
and that he/she is familiar with the calculation of ground state and vibrational properties in metals.
It goes without saying that one should have read the
[fourth lesson on Al](base4) as well as the [introduction page for the EPH code](eph_intro)
before running these examples.

This lesson should take about 1 hour.

## Formalism and connection with the implementation

In the so-called double-delta approximation, the phonon linewidth $\gamma_{\qq\nu}$ due 
to the interaction of the ${\qq\nu}$ phonon with electrons is given by

\begin{equation}
    \gamma_{\qq\nu} = 2\pi \ww_{\qq\nu} \sum_{mn\kk} |g_{mn\nu}(\kk, \qq)|^2
    \delta(\ee_{\kpq m} -\ee_F) \delta(\ee_{\kk n} -\ee_F)
\end{equation}

where $\ww_{\qq\nu}$ is the phonon frequency,
the sum over the electron wavevector $\kk$ runs over the full BZ, $\ee_F$ is the Fermi level at T = 0
and $g_{mn\nu}(\kk, \qq)$ are the e-ph matrix elements discussed in the [EPH introduction](eph_intro).
For a given phonon wavevector $\qq$, the double delta restricts the BZ integration to
transitions between $\kk$ and $\kq$ electron states on the Fermi surface (FS).
Converging the double-delta integral therefore requires very dense $\kk$-meshes in order to 
capture enough states on the FS.
The convergence rate is indeed much slower that the one required by the electron DOS at $\ee_F$:

\begin{equation}
    g(\ee_F) = \sum_{n\kk} \delta(\ee_F - \ee_{n\kk})
\end{equation}

in which a single Dirac delta is involved.

At the level of the implementation, the integration of the double delta can be performed either with the tetrahedron scheme
or by replacing the Dirac delta with a Gaussian function of finite width.
The integration algorithm is defined by the [[eph_intmeth]] input variable.
In the case of Gaussian method, one can use a fixed broadening ([[eph_fsmear]] > 0)
or an adaptive scheme ([[eph_fsmear]] < 0) in which the broadening is automatically computed from 
the electron group velocity [[cite:Li2015]].
The tetrahedron method is more accurate and does not require any broadening parameter. 
Note, however, that in the present implementation the computational cost of the double delta with the tetrahedron method
quickly increases with the size of the $\kk$-mesh so the adaptive Gaussian scheme represents a valid alternative,
especially when a dense $\kk$-sampling is used.

The value of the Fermi energy, $\ee_F$, is automatically computed from the KS eigenvalues stored 
in the input WFK file according to the two input variables [[occopt]] and [[tsmear]].
These parameters are usually equal to the ones used for the GS/DFPT calculation.
However, it is possible to change the value of $\ee_F$ at the EPH level using three (mutually exclusive) input variables:
[[eph_fermie]], [[eph_extrael]] and [[eph_doping]].

<!--
From a numerical point of view, the summation can be limited to states within an energy window around
the Fermi level that can be specified via [[eph_fsewin]].


!!! important

    The $\kk$-mesh used by the EPH code is defined by the input variables [[ngkpt]], [[nshiftk]] and [[shiftk]]
    Note, however, that the $\kk$-sampling  must correspond to the mesh used to generated the input WFK file.
    Convergence studies are performed by generating different WFK files on $\kk$-meshes of increasing density.
-->

The code computes $\gamma_{\qq\nu}$ for each $\qq$-point in the IBZ associated to 
a $\qq$-mesh that can be changed by the user.
By default, the code uses the [[ddb_ngqpt]] $qq$-mesh corresponding to the DDB file (assumed to be equal to the one
used to generate the DVDB file). 
In this case, all the DFPT scattering potentials are available and no interpolation in $qq$-space in required.
To increase the $qq$-sampling, one simply specifies [[eph_ngqpt_fine]] in the input file while 
In this case, the code employs the Fourier interpolation to obtain the scattering potentials.

Once the phonon linewidths $\gamma_{\qq\nu}$ are known in the IBZ, EPH computes the Eliashberg function defined by:

\begin{equation}
    \alpha^2F(\ww) = -\dfrac{1}{N_F} \sum_{\qq\nu} \dfrac{\gamma_{\qq\nu}}{\ww_{\qq\nu}} \delta(\ww - \ww_{\qq \nu})
\end{equation}

where $N_F$ is the density of states (DOS) per spin at the Fermi level.
$\alpha^2F(\ww)$ gives the strength by which a phonon of energy $\ww$ scatters electronic
states on the FS (remember that Abinit uses atomic units hence $\hbar = 1$).
This quantity is accessible in experiments and experience has shown that
$\alpha^2F(\ww)$ is qualitatively similar to the phonon DOS $F(\ww)$:

\begin{equation}
    F(\ww) = \sum_{\qq\nu} \delta(\ww - \ww_{\qq \nu})
\end{equation}

This is not surprising as the equation for $\alpha^2F(\ww)$ resembles the one for the phonon DOS $F(\ww)$:
except for the weighting factor $\frac{\gamma_{\qq\nu}}{\ww_{\qq\nu}}$.
The technique used to compute $\alpha^2F(\ww)$ is defined by the two variables [[ph_intmeth]] and [[ph_smear]].
By default, the code uses the tetrahedron method for the $\qq$-space integration. 

The total e-ph coupling strength $\lambda$ is defined as the first inverse moment of $\alpha^2F(\ww)$:

\begin{equation}
    \lambda = \int \dfrac{\alpha^2F(\ww)}{\ww}\dd\ww = \sum_{\qq\nu} \lambda_{\qq\nu}
\end{equation}

where we have introduced the mode dependent coupling strength:
<!-- For spin unpolarized systems: -->

\begin{equation}
    \lambda_{\qq\nu} = \dfrac{\gamma_{\qq\nu}}{\pi N_F \ww_{\qq\nu}^2}
\end{equation}

Finally, the isotropic superconducting temperature $T_c$ can be estimated using the McMillan expression:

\begin{equation}
    T_c = \dfrac{\ww_{log}}{1.2} \exp \Biggl [
        \dfrac{-1.04 (1 + \lambda)}{\lambda ( 1 - 0.62 \mu^*) - \mu^*}
    \Biggr ]
\end{equation}

where $\mu^*$ is a semi-empirical variable that descrives the (screened) e-e intercation while
$\ww_{\text{log}}$ is  the *logarithmic* average of the phonon frequencies given by:

\begin{equation}
    \ww_{\text{log}} = \exp \Biggl [ \dfrac{2}{\lambda} \int \dfrac{\alpha^2F(\ww)}{\ww}\log(\ww)\dd\ww \Biggr ]
\end{equation}

!!! important

    The code computes $\gamma_{\qq\nu}$ for all the $\qq$-points in the IBZ associated to a homogeneous mesh
    as these quantities are then used to evalute integrals in $\qq$-space.
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

## Preliminary steps

In this tutorial, we prefer to focus on e-ph calculations and the associcated convergence studies.
For this reason, we rely on **pre-computed DEN.nc, DDB and DFPT POT1.nc files** to bypass the DFPT part.
The DEN.nc file will be used to perform NSCF computations on arbitrarily dense $\kk$-meshes while the
POT1.nc files will be merged with the *mrgdv* utility to produce the DVDB database of scattering potentials.

Note that these files are not shipped with the official ABINIT tarball as they are relatively 
large in size.
In order to run the examples of this tutorial, you need to download these files from an external github repository.

If git is installed on your machine, one can easily fetch the repository with:

```sh
git clone ...
```

Alternatively one can use *wget*:

```sh
wget 
```

or *curl*:

```sh
curl
```

or simply copy the tarball by clicking the "download button" in the github interface.
Note that the directory with the input files must be located in the same working directory as the one in which
you will be executing the tutorial.

The input file of the GS run is also stored in the DEN.nc file and one can easily access it with the 
*ncdump* utility

!!! info

    To produce these files, we used the experimental parameters for hexagonal $MgB_2$ (a = 5.8317 and c/a= 1.1419)
    and norm-conserving pseudopotentials with an energy cutoff [[ecut]] of 60 Ry.
    All the calculations have been performed with a 40x40x40 [[ngkpt]] Gamma-centered grid for electrons,
    and the Gaussian smearing [[occopt]] with [[tsmear]].
    The DFPT computations have been done for a set of XXX irreducible $\qq$-points
    corresponding to a $\Gamma$-centered 6x6x6 mesh.
    This is the |AbiPy| script used to automate the GS + DFPT calculation:

<!--
!!! important

    It is clear that, in more realistic applications, one should monitor the convergence of
    lattice parameters and vibrational properties wrt to the $\kk$-mesh and the [[tsmear]] broadening
    before embarking on EPH calculations.
    This is especially true in metals in which a correct description of the FS is crucial.
    In this tutorial, we are trying to find some kind of compromise between accuracy and computational cost.
-->

To merge the POT1 files, execute the *mrgdv* tool using:

```sh
mrgdv <
```

with the following input file that lists all the partial POT1 files already computed for you:

    XXX

As mentioned in the [introduction page for the EPH code](eph_intro), the DVDB file is needed to recostruct
the scattering potentials in the full BZ for all the 3 x [[natom]] atomic perturbations.
With this file, one can obtain $\Delta_\qnu V^\KS$ for all the $\qq$-points belonging to the initial 6x6x6 DFPT
$\qq$-mesh or even a much denser $\qq$-sampling when the Fourier interpolation of the potentials is exploited.
Still the computation of the e-ph matrix elements in the approach implemented in EPH requires the explicit
knowledge of Bloch states at $\kk$ and $\kq$.
This is the problem we will address in the next section.

<!--
At this point, we have all the ingredients required to compute phonon linewidths and $T_c$.
As mentioned in the introductions $\gamma_{\qq\nu}$ coverge slowly wrt to the $\kk$-sampling
so we need to monitor the convergence wrt the number of wavevectors for electrons.

We will try to get an initial rough estimate by looking at the BZ sampling required to
properly descrive quantities that are relatively easy to compute:

1. The electronic DOS at the Fermi level
2. The phonon DOS.

The first test gives us an idea of the $\kk$-mesh whereas the second test gives us an idea
of the $\qq$-mesh.

Let's start from the electron DOS...
-->

### Battle plan

Our goal is to perform calculations of $\gamma_{\qq\nu}$ with different $\qq/\kk$-meshes to
analyse the convergence behaviour of $\lambda$ and $\alpha^2F(\ww)$.
Remember that in the EPH code the $\qq$-mesh can be changed at will thanks to the Fourier interpolation
of the dynamical matrix and of the DFPT potentials whereas the $\kk$-mesh must correspond 
to the one used to generate the input WKF file.
As the $\kk$-mesh must be a multiple of the $\qq$-mesh, we need to generate different WFK files 
in order to perform our convergence studies.

Before running calculations with difference meshes, it is a good idea
to sketch the kind of convergence study we want to perform.
Let's assume we want to test the following configurations:

A 6x6x6 $\qq$-mesh (without interpolation) and the following $\kk$-meshes:

    12x12x12, 18x18x18, 24x24x24 (36x36x36)

A 12x12x12 $\qq$-mesh (with interpolation) and the following $\kk$-meshes:

    12x12x12, 24x24x24, (36x36x36)

The pattern at this point should evident: we need to perform NSCF computations
for four different Gamma-centered $\kk$-meshes: 12x12x12, 18x18x18, 24x24x24, and (48x48x48).


## Electronic properties

This is what is done in the XXX input file.


Run the calculation with as it will take time using:

```
abinit
```

All the NSCF calculations will start from the **pre-computed DEN.nc file** via the [[getden_filepath]] input variable.
Note the usage of the new input variable [[structure]] (added in Abinit v9) to read the crystalline structure from 
an external file in order to avoid repeating the unit cell in each input file.

```sh
# Read MgB2 structure from external file
structure "abivars:test_submodule/mgb2.ucell"

pp_dirpath "$ABI_PSPDIR"
pseudos "12mg.pspnc, 05b.soft_tm"
```

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


```sh
abistruct.py kpath mgb2_DEN.nc
```

To compute the electron DOS, one can use the *abitk* (ABInit ToolKit) utility located in *src/98_main*

```sh
abitk ebands_dos
```

We can also produce a BXSF file that can be used to visualize the Fermi surface with [[prtfsurf]]

```sh
abitk ebands_bxsf
```

To summarize:

### Additional exercise:

First of all, we have to select the $\kk$-path for the NSCF run by using
the [[kptbounds]] and [[ndivsm]] input variables.

## Vibrational properties

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

Since Abinit supports multidatases, unlike anaddb, it's easy to define an input file to compute 
the phonon DOS with multiple $\qq$-meshes.
This simple test allows us to get an preliminary estimate of the $\qq$-sampling required to convergence the 
Eliashberg function as $\alpha^2F(\ww)$ is essentially a weighted phonon DOS.

The input file xxx, shows how to perform such a test with $\qq$-meshes of increasing density.

To analyse the results one can extract the data from the ... PHDOS.nc files

This is what you should get:

!!! Important

    Multiple datasets may be handy when running small calculations as in this case.
    Remember, however, that the EPH code is not designed to be used with multiple datasets.
    In a nutshell, try to split your input files as much as possible so that they can be executed in parallel
    with different number of CPUS and different amount of memory.

Remember to discuss k-mesh and tsmear at the DFPT level.

## Our first isotropic Tc computation

Now we can finally run our first Eliashberg calculation.
We start with a relatively simple input that allows us to introduce the most important variables
and discuss the most important sections of the main output file.

To compute $\gamma_{\qq\nu}$ in metals, one has to use [[optdriver]] 7 and [[eph_task]] 1.
As usual, the location of the DDB, DVDB and WFK files is given by
[[getddb_filepath]] [[getdvdb_filepath]] [[getwfk_filepath]], respectively.

[[eph_intmeth]]

Let's now discuss in more detail the output files produced by the code.

## Increasing the q-mesh

[[eph_ngqpt_fine]]

## Increasing the k-mesh

## Notes on the MPI parallelism

The EPH code supports 5 different levels of MPI parallelism and the number of MPI processes for each level 
can be specified via [[eph_np_pqbks]].
This variable is optional in the sense that whatever number of MPI processes you use, EPH will try to select
a reasonable distribution of the workload.
The distribution, however, may not be optimal as EPH tries to minimize memory requirements by focusing on the 
perturbation/k-point parallelism.
As usual, MPI algorithms are quite efficent if the distribution of the workload is done at a very high-level.
In the case of $T_c$ calculations, the outermost loop is over the $\qq$-points in the IBZ hence the highest speedup 
is achieved when most of the CPUs are the used for the $\qq$-point parallelism.
Note, however, that this kind of MPI distribution does not distribute the wavefunctions and the scattering potentials.

## How to reduce the number of k-points to be computed in the NSCF run

TODO
As we have seen, Eliashberg calculations require the knowledge of Bloch states 
inside a relatively small energy window around  $\ee_F$.
The NSCF computation of the WFK files with dense $\kk$-sampling
[[sigma_erange]]
