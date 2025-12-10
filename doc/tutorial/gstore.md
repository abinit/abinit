---
authors: MG
---

## GSTORE.nc file

In the initial implementation of the EPH code, the e-ph matrix elements were computed on the fly while evaluating the integrals
defining the physical properties of interest.
This approach has the advantage that ABINIT can automatically handle several key tasks, such as applying symmetry operations
to reduce the number of $\kk$ and $\qq$ points to the appropriate irreducible Brillouin zone,
or automatically filtering the bands in transport calculations.

However, this strategy also has important drawbacks.
The e-ph matrix elements must be recomputed from scratch every time a new physical quantity is evaluated.
More critically, there exist algorithms in which the same set of e-ph matrix elements is required multiple times.
A notable example is the VarPEq algorithm: here an external SCF loop is present, and at each iteration
the code must evaluate terms that depend on a fixed set of e-ph matrix elements.

To overcome this limitation, starting from version ??, ABINIT now provides the capability to precompute the $\gkq$
matrix elements using a dedicated EPH sub-driver that is activated using:

[[optdriver]] 7   # Enter EPH code.
[[eph_task]] 11   # GSTORE computation

It is important to understand, however, that the user is now responsible for specifying
how the $\kk$-mesh and $\qq$-mesh should be sampled and how symmetries should be applied to reduce the number of matrix elements.
All variables controlling the GSTORE computation start with the `gstore_` prefix.
Default values are provided that are generally well suited for standard electronic-structure workflows,
but in many situations you may need to customize or override the default behavior depending on your specific use case.
This guide aims to help you understand how to select the appropriate options.

The first step is to specify whether the $\kk$-points or the $\qq$-points should be restricted to the IBZ.
This is controlled by the variables [[gstore_kzone]] and [[gstore_qzone]].
The default behavior is:

[[gstore_kzone]] = "ibz"
[[gstore_qzone]] = "bz"

These settings are OK if you need to compute electronic properties such as the electron self-energy
required for ZPR or electronic transport calculations.
For the phonon self-energy, on the other hand, on should override the default behaviour using

[[gstore_kzone]] = "bz"
[[gstore_qzone]] = "ibz"


!!! important

    The combination [[gstore_kzone]] = "ibz" with [[gstore_qzone]] = "ibz" is not allowed.
    One usually restricts one wavevector to the IBZ while the other wavevector covers the full BZ.
    Using the BZ for both $\kk$ and $\qq$ is usually used for testing purposes (much slower).
    and it is not recommended for production runs unless you know that the post-processing step
    of the GSTORE does not support symmetries.


An additional reduction of the wavevectors can be achieved with the two
mutually exclusive variables [[gstore_use_lgq]] and [[gstore_use_lgk]].
In some cases, the integration over the BZ can indeed be restricted by symmetry to the irreducible wedge defined by the "external" wavevector.
The following examples will help clarify this point.

The electron self-energy Sigma_\nk is defined by an integration over $\qq$-points in the full BZ,
but one can use the symmetries of the little group of $\kk$ to restrict

$$
Sigma_\nk = \int_BZ d\qq = \int_{IBZ_\kk} w^\kk(q) [...]
$$

In this case, one can use

[[gstore_kzone]] "ibz"
[[gstore_qzone]] "bz"
[[gstore_use_lgk]] 1   # Default is 0


For phonon properties, one should use

[[gstore_kzone]] "bz"
[[gstore_qzone]] "ibz"
[[gstore_use_lgq]] 1   # Default is 0

!!! important

    Not all the e-ph calculations are compatible with the little group filtering.
    Please check the documentation and/or run small test calculations before firing big calculations.


Now we turn to the problem of selecting the bands that enter the e-ph matrix elements.
Several options are available, each tailored to simplify a different type of calculation.
Let us begin with the default behavior.
If no specific option is provided in the input file, ABINIT computes **all** matrix elements with $m$ and $n$ ranging from 1 up to [[nband]].
Clearly, this is rarely what you actually want: not all these transitions are needed to compute the final physical properties.
However, ABINIT cannot (yet) read your mind, so you must **explicitly** specify the band ranges in the input file.

The most basic variable is [[gstore_brange]], which defines the range of the $m$ and $n$ indices
(for each spin channel when [[nsppol]] = 2).
[[gstore_brange]] gives you full control over the bands to include, but it is not always
the most convenient option — especially when the relevant contributions to the physical properties come
from transitions located within an energy window around the Fermi level (as in metals)
or from windows starting at the band edges in semiconductors.
In this case, it is much easier to filter bands using an energy range defined by [[gstore_erange]].

!!! important

    [[gstore_erange]] is not compatible with [[gstore_brange]].


Finally, the [[gstore_kfilter]] variable allows you to apply an additional level of filtering directly on the electronic states.
As before, the most appropriate choice for this option depends strongly on the specific physical property you intend to compute.
There are, indeed, several possibile values TO BE DESCRIBED

We conclude this guide by providing examples of recommended settings for different classes of physical properties.
Note that not all gstore_ variables are explicitly included in these examples,
as we rely on the default behavior whenever appropriate.

For computing the ZPR of the fundamental/direct band gap:

[[gstore_kfilter]] "qprange"  # Compute g(k,q) only for |nk> at the band edges
[[gstore_use_lgk]] 1          # Only q-points in the IBZ_k
[[nband]]                     # Bands for the m index (from 1 up to nband)

## MPI parallelism in gstore computation

The GSTORE computation with [[eph_task]] 11 is parallelized over five different MPI levels.
The user can specify manually the MPI grid using [[eph_np_pqbks]].
In this case the product of the MPI processors along the different dimensions must be equal to the
total number of MPI processes allocated by the user, else the code will stop as idle processes are not supported.

If [[eph_np_pqbks]] is not specified in the input, the code will generate the MPI grid automatically
using the total number of MPI processors and the basic dimensions of the job computed at runtime.

If you decide to enforce your MPI grid with [[eph_np_pqbks]], take into account the following.
The parallelization levels over collinear spins, $\kk$-points and $\qq-points$ are the most efficient ones
but the the number of processors for $\kk$ or $\qq$ points should be adjusted according to the values
of [[gstore_kzone]], [[gstore_qzone]].
To reduce load imbalace, one should use less processors for the wavevector that is being restricted to the IBZ

The parallelism over perturbations should be activated only when the previous three MPI levels start to saturate.
Note that the parallelism over bands is not supported in GSTORE computation.

TODO: GWPT and [[gwpt_np_wpqbks]]
Also [[boxcutmin]] to accelerate computations.

!!! important

    The output of the GSTORE file requires a netcdf library with MPI-IO support.

## How to densify the q-mesh

By default, the e-ph matrix elements are computed using the coarse ab-initio $\qq$-mesh given by [[ddb_ngqpt]].
This is the $\qq$-mesh used in the DFPT calculation.

To densify the $\qq$-mesh, use [[eph_ngqpt_fine]] but remember that
The $\qq$-mesh must be identical to, or a submesh of, the $\kk$-mesh associated with the WFK file.

Further details on the interpolation of the DFPT scattering potentials are available in this section.

## Restarting a GSTORE computation

If your GSTORE calculation has been killed due to e.g. timeout limit,
you can always restart the computation by re-rurring the same input with the
addition of

[[getgstore_filepath]]  "out_GSTORE.nc"

where "out_GSTORE.nc" is the name of the output GSTORE file produced by the calculation
that was interrupted.

## How to compute physical properties from a GSTORE file

So far we have discussed how to generate a GSTORE.nc file.
Now we explain how to read the e-ph matrix elements from file and use them
to compute physical properties.

!!! critical

    Do not use a WFK file different from the one used to generate the GSTORE file.
    The (complex) e-ph matrix elements stored in the GSTORE depend on the gauge
    of the wavefunctions in the WFK file.


Reading a GSTORE file is very easy, use [[getgstore_filepath]] and then select the appropriate
value of [[eph_task]] to perform the post-processing step

[[optdriver]] 7         # Enter EPH code.
[[getgstore_filepath]]  "teph4zpr_10o_DS1_GSTORE.nc"

[[eph_task]] 24         # SIGMAPH from GSTORE

Other options or files may be needed depending on [[eph_task]].
Please consult the documentation or the available tutorials.

The GSTORE.nc stores additional quantities such as phonon frequencies and eigenvectors for all the $\qq$-points in the IBZ,
and additional metadata such as, for instance, a table that specifies whether
all the entries for a particular $\qq$-points have been computed.
This table is used to implement the automatic restart of the computation if the job is killed due to the time-limit.
In our implementation, we are also able to filter the set of $\kk$- and $\qq$-points as well as the set of
$m$ and $n$ bands in the e-ph matrix elements.
The kind of filtering technique that should be used depends on the application in mind.
For metals, for instance, one is usually interested in the e-ph matrix elements only for bands inside
an energy window around the Fermi level.
Moreover one can select only those $\kk$ and $\qq$ for which there is at least one electronic transition
from $\kk$ to $\kk+\qq$ inside the energy window.
%For the computation of the e-ph induced renormalization of the electronic states and the ZPR, on the other hand, one is usually interested in the corrections at the band edges. In this case, one can compute the GWPT matrix elements only for the $\nk$-states of the CBM/VBM, and then evaluate the coupling for all the $\qq$-points in the irreducible zone defined by the little group of $\kk$.
%For the ZPR we have to include a large number of empty states associated to the $m$ index and this clearly increases significantly the computational cost.
%If we assume, however, that the GWPT matrix elements do not differ significantly from the KS ones, one can use the Sternheimer method to account for the contribution to the sum beyond the active space.
