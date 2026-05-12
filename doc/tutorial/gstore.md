---
authors: MG
---

# The GSTORE.nc file

In the initial implementation of the EPH code, the e-ph matrix elements were computed on the fly
while evaluating the integrals defining the physical properties of interest.
This approach has the advantage that ABINIT can automatically handle several key tasks, such as
applying symmetry operations to reduce the number of $\kk$ and $\qq$ points to the
appropriate irreducible Brillouin zone (IBZ), or automatically filtering the bands in transport calculations.

However, this strategy also has important drawbacks since
the e-ph matrix elements must be recomputed from scratch every time a new physical quantity is evaluated.
More critically, there exist algorithms in which the same set of e-ph matrix elements is required multiple times.
A notable example is the [VarPEq algorithm](./eph4vpq.md): here an external SCF loop is present, and at each iteration
the code must evaluate terms that depend on a fixed set of e-ph matrix elements.

To overcome this limitation, ABINIT now provides the capability to precompute the $\gkq$
matrix elements using a dedicated EPH sub-driver that is activated using
[[optdriver]] and [[eph_task]]:

```
optdriver 7   # Enter EPH code.
eph_task 11   # GSTORE computation
```

It is important to understand, however, that the user is now **responsible** for specifying
how the $\kk$-mesh and $\qq$-mesh should be sampled and how symmetries should be applied
to reduce the number of matrix elements.
Default values are provided that are generally well suited for the computation of electronic properties,
but in many situations you may need to customize or override the default behavior depending on your specific use case.
This guide aims to help you understand how to select the appropriate options.

!!! important

    All variables controlling the GSTORE computation start with the `gstore_` prefix.


The first step is to specify whether the $\kk$-points or the $\qq$-points in $\gkq$ should be restricted to the
IBZ or the full Brillouin zone (BZ).
This is controlled by the variables [[gstore_kzone]] and [[gstore_qzone]].
The default behavior is:

```
gstore_kzone = "ibz"
gstore_qzone = "bz"
```

These settings are OK if you need to compute electronic properties such as the electron self-energy
$\Sigma_\nk$ required for the ZPR or electronic transport calculations.
For the phonon self-energy $\Pi_\qnu$, on the other hand, on should override the default behaviour using

```
gstore_kzone = "bz"
gstore_qzone = "ibz"
```

!!! important

    The combination [[gstore_kzone]] = "ibz" with [[gstore_qzone]] = "ibz" **is not allowed**.
    One usually restricts one wavevector to the IBZ while the other wavevector covers the full BZ.
    Using the BZ for both $\kk$ and $\qq$ is usually used for testing purposes (much slower).
    and it is not recommended for production runs unless you know that the post-processing step
    of the GSTORE does not support symmetries.


An additional reduction of the number of wavevectors can be achieved with the two
mutually exclusive variables [[gstore_use_lgq]] and [[gstore_use_lgk]].
In some cases, the integration over the BZ in the post-processing step can indeed be restricted
by symmetry to the irreducible wedge defined by the little group of the "external" wavevector ($\kk$ or $\qq$
We use the notation IBZ_k to denote the the irreducible wedge defined by the little group of $\kk$
and IBZ_q for the irrecudible wedge defined by the little group of $\qq$.
The following examples will help clarify this point.

The electron self-energy $\Sigma_\nk$ is defined by an integration over $\qq$-points in the full BZ,
but one can use the symmetries of the little group of $\kk$ to restrict the integration to a smaller zone.
Schematically:

$$
\Sigma_\nk = \int_{BZ} d\qq = \int_{IBZ_\kk} w^\kk(\qq) [...]
$$

In this case, one can use

```
gstore_kzone "ibz"
gstore_qzone "bz"
gstore_use_lgk 1   # Default is 0
```


For phonon properties, one should use

```
gstore_kzone "bz"
gstore_qzone "ibz"
gstore_use_lgq 1   # Default is 0
```

!!! important

    Not all the e-ph calculations are compatible with the little group filtering.
    Please check the documentation and/or run small test calculations before firing big calculations.


Now we turn to the problem of selecting the bands that enter the e-ph matrix elements.
Several options are available, each tailored to simplify a different type of calculation.
Let us begin with the default behavior.

If no specific option is provided in the input file, ABINIT computes **all** matrix elements
with $m$ and $n$ ranging from 1 up to [[nband]].
Clearly, this is **rarely** what you actually want: not all these transitions are needed to compute the final physical properties.
However, ABINIT cannot (yet) read your mind, so you must **explicitly specify** the band ranges in the input file.

The most basic variable is [[gstore_brange]], which defines the range of **both** the $m$ and $n$ indices
for each spin [[nsppol]] channel.

[[gstore_brange]] gives you full control over the bands to include and should be included when you need
all the e-ph matrix elements connecting all $\kk$ and $\kq$ states inside a band range as, for instance,
in polaron calculations.

In the case of metals or transport properties in semiconductors, the relevant contributions to the physical properties come
from transitions located within an energy window around the Fermi level (as in metals)
or from windows starting at the band edges in semiconductors.
In this case, it is much easier to filter bands automatically using an input energy range defined by [[gstore_erange]].

!!! important

    [[gstore_erange]] is not compatible with [[gstore_brange]].


There are however other applications in which we want to filter $\kk$-points and have a different range of bands
for the $\psi_\nk$ and $\psi_\mkq$ states entering the e-ph matrix elements.

In ZPR calculations of the (fundamental) band gap, for instance, the $\kk$-points and the $n$ index can be restricted
to the band edges that are automatically detected from the KS energies of the WFK file.
The $m$ index, on the contrary, should cover a much large band range to account for empty states in the summation
while the $\qq$-points should cover the full BZ or an appropriate irreducible wedge as discussed below.

Finally, the [[gstore_kfilter]] variable allows you to apply an additional level of filtering directly on the electronic states.
As before, the most appropriate choice for this option depends strongly on the specific physical property you intend to compute.
There are, indeed, several possibile values.
Please refer to the documentation of [[gstore_kfilter]] for further details.

We conclude this guide by providing examples of recommended settings for different classes of physical properties.
Note that not all gstore_ variables are explicitly included in these examples,
as we rely on the default behavior whenever appropriate.

For computing the ZPR of the fundamental/direct band gap, use:

```
gstore_kfilter "qprange"  # Compute g(k,q) only for |nk> at the band edges
gstore_use_lgk 1          # Only q-points in the IBZ_k
gstore_brange 1 12        # Range for the m index (last index cannot be greated than nband)
nband         12
```

!!! important

    Here, by band edges we refer to the $\kk$-points in the WFK file at which the
    conduction band minimum (CBM) and valence band maximum (VBM) are found.

    These $\kk$-points do not necessarily coincide with the true band extrema,
    as the latter may not lie on the chosen $\kk$-mesh. A typical example is
    silicon, where the CBM is located along the Γ–X direction at $k \approx 0.85\, \frac{2\pi}{a}$.

    If a more accurate description of the true band edges is required,
    generate a WFK file using a shifted $\kk$-mesh via [[nshiftk]] and [[shiftk]].


If you want to have full control on the list $\kk$-points and bands that should be considered
for the $|n\kk\rangle$ states, **remove** the "kfilter" option,
and use [[nkptgw]], [[kptgw]] and [[bdgw]] as in the example below:

```
gstore_use_lgk 1  # Only q-points in the IBZ_k
nkptgw 2
kptgw
0   0 0
0.8 0 0

bdgw
1 5               # Range for n index.
1 5

gstore_brange 1 12  # Range for the m index
nband           12  # Last index cannot be greater than nband
```


## MPI parallelism in gstore computation

There are two EPH subdrivers capable of generating a GSTORE.nc file.
[[eph_task]] = 11 computes the e-ph matrix elements at the KS level, whereas
[[eph_task]] = 17 employs the more expensive GWPT formalism [[cite:Li2019]] and produces a GSTORE.nc file
containing both GWPT and KS e-ph matrix elements.

[[eph_task]] 11 is parallelized over five different MPI levels.
The user can specify manually the MPI grid using [[eph_np_pqbks]].
In this case, the product of the MPI processors along the different dimensions must be equal to the
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

In the case of GWPT calculations ([[eph_task]] == 17) the MPI grid is defined by [[gwpt_np_wpqbks]]

!!! important

    The output of the GSTORE file requires a netcdf library with MPI-IO support.

## How to densify the q-mesh

By default, the e-ph matrix elements are computed using the coarse ab-initio $\qq$-mesh given by [[ddb_ngqpt]].
This is the $\qq$-mesh used in the DFPT calculation.

To densify the $\qq$-mesh, use [[eph_ngqpt_fine]].
Note, however, that the $\qq$-mesh must be identical to, or a submesh of, the $\kk$-mesh associated with the WFK file.
Also, be sure to compute the Born effective charges in polar materials,
and the dynamical quadrupoles in order to
properly describe the long-range part of the DFT scattering potentials and obtain a reliable Fourier interpolation.

Further details on the interpolation of the DFPT scattering potentials are available the [eph_intro page](eph_intro.md).

## Restarting a GSTORE computation

If your GSTORE calculation is interrupted, e.g. due to a timeout limit,
you can restart the computation simply by rerunning the same input.
Restart capabilities for GSTORE calculations are enabled by default, but they can be disabled
by setting [[eph_restart]] to 0.

<!--
with the addition of
[[getgstore_filepath]]  "out_GSTORE.nc"

where "out_GSTORE.nc" is the name of the output GSTORE file produced by the calculation
that was interrupted.
-->

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

```
optdriver 7         # Enter EPH code.
eph_task 24         # SIGMAPH from GSTORE

getgstore_filepath  "teph4zpr_10o_DS1_GSTORE.nc"
```


Other options or files may be needed depending on the value of [[eph_task]].
Please consult the documentation and the available tutorials.
