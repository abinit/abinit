---
authors: MG
---

# Superconducting properties within the isotropic Eliashberg formalism

This tutorial explains how to compute phonon linewidths in metals
and use the McMillan equation to estimate the superconducting critical temperature $T_c$ within 
the isotropic Eliashberg formalism.
We start by presenting the basic equations implemented in the code and 
the connection with the ABINIT input variables.
Then we discuss how to run EPH calculations and perform typical convergence studies using $MgB_2$ as example.

It is assumed the user has already completed the two tutorials [RF1](rf1) and [RF2](rf2),
and that he/she is familiar with the calculation of ground state and vibrational properties in metals.
It goes without saying that you should have read the 
[fourth lesson on Al](base4) as well as the [introduction page for the EPH code](eph_intro) 
before running these examples.

This lesson should take about 1 hour.

## Formalism and connection to the implementation

The phonon linewidths are defined by

\begin{equation}
    \gamma_{\qq\nu} = 2\pi \omega_{\qq\nu} \sum_{mn\kk} |g_{mn\nu}(\kk, \qq)|^2 
    \delta(\ee_{\kpq m} -\ee_F) \delta(\ee_{\kk n} -\ee_F)
\end{equation}

<!--
where, in order to simply the notation, the Fermi level $\ee_F$ has been set to zero.
so that $\delta(\ee_\kk -\ef)$ can be replaced by $\delta(\ee_\kk)$.
-->
For each $\qq$-point, the double delta $\delta(\ee_{\kpq m}) \delta(\ee_{\kk n})$
restricts the BZ integration to transitions between $\kk$ and $\kq$ states belonging to the Fermi surface (FS).
From a numerical standpoint, accurate phonon linewidths require a very dense $\kk$-sampling to capture enough 
transitions on the FS.
The convergence rate is expected to be slower that the one required 
for the DOS at $\ee_F$ in which only a single delta is involved.
<!--
As a consequence a very dense sampling is needed to obtain converged results with a $\kk$-mesh
that is usually denser than the one needed to converge 
-->
As concerns the implementation,
the integration of the double delta can be performed either with the tetrahedron scheme 
or by replacing the Dirac delta with a Gaussian method of finite broadening
(see [[eph_intmeth]] and [[eph_fsmear]]).
The value of the Fermi level is automatically taken from the input WFK file but it is possible to change the value 
of $\ee_F$ at runtime using either [[eph_fermie]] or [[eph_extrael]].
<!--
From a numerical point of view, the summation can be limited to states within an energy window around 
the Fermi level that can be specified via [[eph_fsewin]].
-->

!!! important

    At the level of the implementation, the $\kk$-mesh used by the EPH code is defined by the set of 
    input variables [[ngkpt]], [[nshiftk]] and [[shiftk]] 
    Note, however, that the $\kk$-sampling  must correspond to the mesh used to generated the input WFK file. 
    Convergence studies are performed by computing different WFK files on $\kk$-meshes of increasing density.

Once $\gamma_{\qq\nu}$ are know for a set of IBZ $\qq$-points, one can
obtain the Eliashberg function defined by:

\begin{equation}
    \alpha^2F(\omega) = -\dfrac{1}{N_F} \sum_{\qq\nu} \dfrac{\gamma_{\qq\nu}}{\omega_{\qq\nu}} \delta(\ww - \ww_{\qq \nu})
\end{equation}

where $N_F$ is the density of states (DOS) per spin evaluted at the Fermi level.
$\alpha^2F(\omega)$ gives the strength by which a phonon of energy $\ww$ scatters electronic states on the FS.
This quantity is accessible in experiments and experience has shown that, in many cases, 
$\alpha^2F$ shows the same trend as the phonon DOS.
As a matter of fact, the equation for $\alpha^2F$ is similar to the one for the phonon DOS $F(\omega)$:

\begin{equation}
    F(\ww) = \sum_{\qq\nu} \delta(\ww - \ww_{\qq \nu})
\end{equation}

except for the weighting factors $\frac{\gamma_{\qq\nu}}{\omega_{\qq\nu}}$ that is responsible 
for possible deviation between $\alpha^2F(\ww)$ and $F(\ww)$.

The relevant variables are: [[ph_intmeth]] [[ph_smear]] .

!!! important

    The code computes $\gamma_{\qq\nu}$ for all the $\qq$-points in the IBZ associated to a homogeneous mesh
    as these quantities are then used to evalute integrals in $\qq$-space.
    [[ddb_ngqpt]] mesh that is the $\qq$-mesh used in the DFPT computations.
    In this case, the code does not performy any kind of interpolation in $\qq$-space.
    A (much denser) $\qq$-mesh can be specified via [[eph_ngqpt_fine]].
    In this case, the code interpolates the DFPT potentials on the dense mesh at runtime.
    EPH computes the electron DOS using the eigenvalues stored in the WFK file.

The total e-ph coupling strength $\lambda$ is defined as the first inverse moment of $\alpha^2F(\ww)$:

\begin{equation}
    \lambda = \int \dfrac{\alpha^2F(\omega)}{\omega}\dd\omega = \sum_{\qq\nu} \lambda_{\qq\nu}
\end{equation}

where we have introduced the mode dependent coupling strength:
<!-- For spin unpolarized systems: -->

\begin{equation}
    \lambda_{\qq\nu} = \dfrac{\gamma_{\qq\nu}}{\pi N_F \ww_{\qq\nu}^2}
\end{equation}

Finally, the isotropic superconducting temperature can be estimated using the McMillan expression:

\begin{equation}
    T_c = \dfrac{\omega_{log}}{1.2} \exp \Biggl [ 
        \dfrac{-1.04 (1 + \lambda)}{\lambda ( 1 - 0.62 \mu^*) - \mu^*} 
    \Biggr ]
\end{equation}

where $\mu^*$ is a semi-empirical variable and 
the *logarithmic* average of the phonon frequencies, $\omega_{\text{log}}$, is given by:

\begin{equation}
    \omega_{\text{log}} = \exp \Biggl [ \dfrac{2}{\lambda} \int \dfrac{\alpha^2F(\omega)}{\omega}\log(\omega)\dd\omega \Biggr ]
\end{equation}

<!--
Nesting factor

\begin{equation}
    N(\qq) = \sum_{mn\kk} \delta(\ee_{\kpq m}) \delta(\ee_{\kk n})
\end{equation}


Implementation details:

\begin{equation}
    \gamma_{\qq\nu} = 2\pi \omega_{\qq\nu} \sum_{pp'} d_{\qq p}^* \tilde\gamma_{p p'}(\qq) d_{\qq p'}^*
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

In this tutorial, we prefer to focus on e-ph calculations rather that trying to compute 
all the ingredients from scratch as done in other tutorials.
For this reason, we will rely on **pre-computed** DEN.nc, DDB and POT1.nc files to bypass an expensive DFPT computation.
The DEN.nc file are used to perform NSCF computations on dense $\kk$-meshes while the 
POT1.nc files must be merged with the *mrgdv* utility to produce the DVDB database required by the EPH code.

To produce these files, we used the experimental structure parameters for hexagonal MgB2 (a = 5.8317 and c/a= 1.1419) 
and norm-conserving pseudopotentials with an energy cutoff [[ecut]] of 60 Ry. 
All the calculations have been performed with a 40x40x40 grid for electrons k [[ngkpt]], 
and the Gaussian smearing [[occopt]] with [[tsmear]].
The DFPT computations have been done for a set of XXX irreducible $\qq$-points 
corresponding to a $\Gamma$-centered 6x6x6 mesh.

TODO: Mention *ncdump* and the input file 

!!! important

    It is clear that, in more realistic applications, one should monitor the convergence of 
    lattice parameters and vibrational properties wrt to the $\kk$-mesh and the [[tsmear]] broadening 
    before embarking on EPH calculations.
    This is especially true in metals in which a correct description of the FS is crucial.
    In this tutorial, we are trying to find some kind of compromise between computational cost, disk-space 


To merge the POT1 files, run the *mrgdv* tool using:

```sh
mrgdv < 
```

with the following input file that lists all the partial POT1 files we have already computed for you:

    XXX


As mentioned in the [introduction page for the EPH code](eph_intro), the DVDB file is used to recostruct
the scattering potential in the full BZ. 
With this file, one can obtain $\Delta_\qnu V^\KS$ for all the $\qq$-points belonging to the initial 6x6x6 DFPT 
$\qq$-mesh or a much denser $\qq$-sampling if the Fourier interpolation of the potentials is exploited.
Still the computation of the e-ph matrix elements in the EPH approach requires the explicit 
knowledge of Bloch states and this is the problem that will be addressed in the next section.

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

## Electronic properties

Our goal is to perform calculations of $\gamma_{\qq\nu}$ with different $\qq/\kk$-meshes to 
analyse the convergence behaviour.
Remember, however, that in EPH the $\qq$-mesh can be changed at will thanks to the Fourier interpolation
whereas the $\kk$-mesh is bound to the input WKF file.
So it is always a good idea to sketch the kind of convergence study you want to perform before running calculations.

Let's assume we want to test the following configurations:

A 6x6x6 $\qq$-mesh (without interpolation) and the following $\kk$-meshes:

    12x12x12, 18x18x18, 24x24x24 (36x36x36)

A 12x12x12 $\qq$-mesh (with interpolation) and the following $\kk$-meshes:

    12x12x12, 24x24x24, (36x36x36)

The pattern at this point should evident: we need to perform NSCF computations
for four different Gamma-centered $\kk$-meshes: 12x12x12, 18x18x18, 24x24x24, and (48x48x48).
This is what is done in the XXX input file.

Run the calculation with as it will take time using:

```
abinit
```

All the NSCF calculations will start from the pre-computed DEN.nc file via the [[getden_filepath]] input variable.
We also use the new input variable [[structure]] (added in Abinit v9) to avoid 
repeating structral information in each input file.

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

In this section, we compute vibrational properties using the pre-generated DDB file.
Instead of using *anaddb* to Fourier interpolate the dynamical matrix, we employ 
the interface provided by the EPH code.

The input file xxx, computes the phonon DOS using different $\qq$-meshes
and band structures 

To compute vibrational properties without performing a full e-ph calculation, 
use: [[optdriver]] = 7 with [[eph_task]] = 0

[[getddb_filepath]] 
[[ddb_ngqpt]]

[[ph_ngqpt]]
[[ph_intmeth]]

[[ph_nqpath]]
[[ph_qpath]]
[[ph_ndivsm]]

[[dipdip]]

In a nustshell, one should add:

```sh
optdriver 7
eph_task 0 

ph_nqpath
ph_qpath
ph_ndivsm

ph_ngqpt
ph_intmeth
```

To summarize:
Remember to discuss k-mesh and tsmear at the DFPT level.

## Our first EPH computation 

Now we can finally run our first EPH calculation.
We start with a relatively simple input that allows us to introduce the most important variables
and discuss the output file. 

To compute $\gamma_{\qq\nu}$ in metals, one has to use [[optdriver]] 7 and [[eph_task]] 1.
As usual, the location of the DDB, DVDB and WFK files is given by 
[[getddb_filepath]] [[getdvdb_filepath]] [[getwfk_filepath]], respectively.

[[eph_intmeth]]

Let's now discuss in more detail the output files produced by the code.

## Increasing the q-mesh

[[eph_nqgpt_fine]]

## Notes on the MPI parallelism

[[eph_np_pqbks]]

## How to reduce the number of k-points.
