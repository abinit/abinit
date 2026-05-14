---
authors: MG
---

# The GSTORE.nc file

In the initial implementation of the EPH code, the e-ph matrix elements were computed on the fly
while evaluating the integrals defining the physical properties of interest.
This approach has the advantage that ABINIT can automatically handle several key tasks, such as
applying symmetry operations to reduce the number of $\kk$ and $\qq$ points to the
appropriate irreducible Brillouin zone (IBZ), or automatically filtering the bands in transport calculations.

However, this strategy has significant drawbacks, as the e-ph matrix elements must be recomputed
every time a different physical quantity is evaluated.
More critically, some algorithms require the same set of e-ph matrix elements multiple times.
A notable example is the [VarPEq algorithm](./eph4vpq.md), which involves an external SCF loop where,
at each iteration, the code must evaluate terms that depend on a fixed set of e-ph matrix elements.

To overcome this limitation, ABINIT now provides the capability to precompute the $\gkq$
matrix elements using a dedicated EPH sub-driver that is activated using [[optdriver]] and [[eph_task]]:

```
optdriver 7   # Enter EPH code.
eph_task 11   # GSTORE computation
```

It is important to understand, however, that the user is now **responsible**
for specifying how the $\kk$-mesh and $\qq$-mesh should be sampled,
and how symmetries should be applied to reduce the number of matrix elements.
Default values are provided that are generally well-suited for electronic properties,
but you may need to customize these settings depending on your specific use case.
This guide explains how to select the most appropriate options.

!!! important

    All variables controlling the GSTORE computation start with the `gstore_` prefix.


The first step is to specify whether the $\kk$-points or the $\qq$-points in $\gkq$ should
be restricted to the IBZ or the full Brillouin zone (BZ).
This is controlled by the variables [[gstore_kzone]] and [[gstore_qzone]].
The default behavior is:

```
gstore_kzone = "ibz"
gstore_qzone = "bz"
```

These settings are appropriate if you need to compute electronic properties such as the
electron self-energy $\Sigma_\nk$ required for ZPR or electronic transport calculations.
For the phonon self-energy $\Pi_\qnu$, on the other hand, you should override the default behavior using:

```
gstore_kzone = "bz"
gstore_qzone = "ibz"
```

!!! important

    The combination [[gstore_kzone]] = "ibz" with [[gstore_qzone]] = "ibz" **is not allowed**.
    One usually restricts one wavevector to the IBZ while the other covers the full BZ.
    Using the BZ for both $\kk$ and $\qq$ is mainly for testing purposes (it is much slower)
    and is not recommended for production runs unless you are certain
    that the post-processing step does not support symmetries.


An additional reduction of the number of wavevectors can be achieved with the two
mutually exclusive variables [[gstore_use_lgq]] and [[gstore_use_lgk]].

In some cases, the integration over the BZ during post-processing can be restricted
by symmetry to the irreducible wedge defined by the little group of the "external" wavevector ($\kk$ or $\qq$).
We use the notation IBZ_k to denote the irreducible wedge defined by the little group of $\kk$,
and IBZ_q for the irreducible wedge defined by the little group of $\qq$. The following examples clarify this approach.

The electron self-energy $\Sigma_\nk$, for instance, is defined by an integration over $\qq$-points in the full BZ,
but one can use the symmetries of the little group of $\kk$ to restrict the integration to a smaller zone.
Schematically:

$$
\Sigma_\nk = \int_{BZ} d\qq [...] = \int_{IBZ_\kk} d\qq w^\kk(\qq) [...]
$$

In this case, one can use

```
gstore_kzone "ibz"
gstore_qzone "bz"
gstore_use_lgk 1   # Default is 0
```

For phonon properties, on the contrary, one should use

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

If no specific option is provided, ABINIT computes **all** matrix elements with $m$ and $n$
ranging from 1 up to [[nband]].
This is **rarely** ideal, as not all transitions are needed for the final physical properties.
Therefore, you should **explicitly specify** the band ranges in the input file.

The most basic variable is [[gstore_brange]], which defines the range of **both** the $m$ and $n$ indices
for each **collinear** spin channel ([[nsppol]] = 2).

[[gstore_brange]] gives you full control over the bands to include and should be used when you need
all the e-ph matrix elements connecting all $\kk$ and $\kq$ states inside a band range as, for instance,
in polaron calculations.

In the case of metals or transport properties in semiconductors, the relevant contributions
come from transitions within an energy window around either the Fermi level (for metals)
or the band edges (for semiconductors).
In these cases, it is easier to filter bands automatically using an energy range defined by [[gstore_erange]].

!!! important

    [[gstore_erange]] is not compatible with [[gstore_brange]].


There are, however, other applications where we want to filter $\kk$-points and use different band ranges
for the $\psi_\nk$ and $\psi_\mkq$ states entering the e-ph matrix elements.

In ZPR calculations of the (fundamental) band gap, for instance, the $\kk$-points and the $n$ index can be restricted
to the band edges that are automatically detected from the KS energies of the WFK file.
The $m$ index, on the contrary, should cover a much larger band range to account for empty states in the summation,
while the $\qq$-points should cover the full BZ or an appropriate irreducible wedge as discussed below.

Finally, the [[gstore_kfilter]] variable allows you to apply an additional level of
filtering directly to the electronic states.
The most appropriate choice depends heavily on the specific physical property being computed;
indeed, several values are possible. Please refer to the [[gstore_kfilter]] documentation for more details.

Finally, the [[gstore_with_vk]] variable allows you to include electronic group velocities
(and optionally off-diagonal velocity matrix elements) in the GSTORE file.
This is particularly useful for transport calculations or when the velocity gauge must be consistent
with the wavefunctions used for the e-ph matrix elements.

The following examples provide recommended settings for different classes of physical properties.
Note that not all gstore_ variables are explicitly included in these examples,
as we rely on the default behavior whenever appropriate.

For computing the ZPR of the fundamental/direct band gap, use:

```
gstore_kfilter "qprange"  # Compute g(k,q) only for |nk> at the band edges
gstore_use_lgk 1          # Only q-points in the IBZ_k
gstore_brange 1 12        # Range for the m index (last index cannot be greater than nband)
nband         12
```

!!! important

    By "band edges," we refer to the $\kk$-points in the WFK file where the conduction
    band minimum (CBM) and valence band maximum (VBM) are found.

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


## Gstore computation and MPI parallelism

There are two EPH subdrivers capable of generating a GSTORE.nc file:
[[eph_task]] = 11 computes the e-ph matrix elements at the KS level, whereas
[[eph_task]] = 17 employs the more expensive GWPT formalism [[cite:Li2019]] and produces
a GSTORE.nc file containing both GWPT and KS e-ph matrix elements.

[[eph_task]] 11 is parallelized over five different MPI levels.
The user can manually specify the MPI grid using [[eph_np_pqbks]].
In this case, the product of the MPI processes along the different dimensions must equal
the total number of MPI processes; otherwise, the code will stop because idle processes are not supported.
If [[eph_np_pqbks]] is not specified, the code generates the MPI grid automatically using
the total number of MPI processes and the dimensions of the job at runtime.

If you decide to enforce your MPI grid with [[eph_np_pqbks]], consider the following:
parallelization over collinear spins, $\kk$-points, and $\qq$-points is the most efficient,
but the number of processes for $\kk$ or $\qq$ points should be adjusted according to [[gstore_kzone]] and [[gstore_qzone]].
To reduce load imbalance, use fewer processes for the wavevector restricted to the IBZ. Parallelism
over perturbations should be activated only when the previous three MPI levels saturate.
Note that parallelism over bands is not supported in the current version.

In the case of GWPT calculations ([[eph_task]] == 17) the MPI grid is defined by [[gwpt_np_wpqbks]]

!!! important

    The output of the GSTORE file requires a netcdf library with MPI-IO support.

## How to densify the q-mesh

By default, the e-ph matrix elements are computed using the coarse ab-initio $\qq$-mesh given by [[ddb_ngqpt]].
This is the $\qq$-mesh used in the DFPT calculation.

To densify the $\qq$-mesh, use [[eph_ngqpt_fine]].
Note, however, that the $\qq$-mesh must be identical to, or a submesh of, the electronic $\kk$-mesh.
Also, ensure you compute Born effective charges and dynamical quadrupoles to properly describe
the long-range part of the DFT scattering potentials and obtain reliable Fourier interpolation.
Further details can be found in the [EPH intro tutorial](eph_intro.md).

## Restarting a GSTORE computation

If your GSTORE calculation is interrupted, e.g. due to a timeout limit,
you can restart the computation simply by rerunning the same input.
Restart capabilities for GSTORE calculations are enabled by default, but they can be disabled
by setting [[eph_restart]] to 0.

## How to compute physical properties from a GSTORE file

So far we have discussed how to generate a GSTORE.nc file.
Now we explain how to read the e-ph matrix elements from file and use them
to compute physical properties.

!!! critical

    Do not use a WFK file different from the one used to generate the GSTORE file.
    The (complex) e-ph matrix elements stored in the GSTORE depend on the gauge
    of the wavefunctions in the WFK file.


Reading a GSTORE file is straightforward: use [[getgstore_filepath]] and then
select the appropriate [[eph_task]] to perform the post-processing step.
To compute the ZPR from GSTORE, use e.g.:

```
optdriver 7         # Enter EPH code.
eph_task 24         # SIGMAPH from GSTORE

getgstore_filepath  "teph4zpr_10o_DS1_GSTORE.nc"
```

Other options or files may be needed depending on the value of [[eph_task]].
Please consult the documentation and the available tutorials.
