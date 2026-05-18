---
authors: MG
---

# The GSTORE.nc file

In the initial implementation of the EPH code, electron-phonon (e-ph) matrix elements were computed on the fly during the evaluation of integrals for various physical properties. This approach allowed ABINIT to automatically handle tasks such as applying symmetry operations to reduce the number of $\kk$ and $\qq$ points to the irreducible Brillouin zone (IBZ) or filtering bands in transport calculations.

However, this strategy had significant drawbacks, as matrix elements had to be recomputed every time a different physical quantity was evaluated. Furthermore, certain algorithms—such as the [VarPEq algorithm](./eph4vpq.md)—require the same set of e-ph matrix elements multiple times across an external SCF loop.

To overcome these limitations, ABINIT now allows you to precompute $\gkq$ matrix elements using a dedicated EPH sub-driver activated via [[optdriver]] and [[eph_task]]:
```
optdriver 7   # Enter EPH code.
eph_task 11   # GSTORE computation
```

When using this approach, the user is responsible for specifying how the $\kk$-mesh and $\qq$-mesh should be sampled and how symmetries should be applied to reduce the total number of matrix elements. While default values are generally suitable for electronic properties, you may need to customize these settings for specific use cases. This guide explains how to select the most appropriate options.

!!! important
    All variables controlling the GSTORE computation start with the `gstore_` prefix.

The first step is to specify whether the $\kk$-points or $\qq$-points in $\gkq$ should be restricted to the IBZ or cover the full Brillouin zone (BZ). This is controlled by [[gstore_kzone]] and [[gstore_qzone]]. The default behavior is:
```
gstore_kzone = "ibz"
gstore_qzone = "bz"
```
These settings are appropriate for computing electronic properties, such as the electron self-energy $\Sigma_\nk$ required for ZPR or electronic transport calculations. For the phonon self-energy $\Pi_\qnu$, however, you should override these defaults using:
```
gstore_kzone = "bz"
gstore_qzone = "ibz"
```

!!! important
    The combination [[gstore_kzone]] = "ibz" and [[gstore_qzone]] = "ibz" is not allowed. Typically, one wavevector is restricted to the IBZ while the other covers the full BZ. Using the full BZ for both $\kk$ and $\qq$ is primarily for testing (as it is much slower) and is not recommended for production unless you are certain the post-processing step does not support symmetries.

You can achieve a further reduction in the number of wavevectors using the mutually exclusive variables [[gstore_use_lgq]] and [[gstore_use_lgk]]. In some cases, BZ integration during post-processing can be restricted by symmetry to the irreducible wedge defined by the little group of the "external" wavevector ($\kk$ or $\qq$). We denote these wedges as $IBZ_\kk$ and $IBZ_\qq$, respectively.

For instance, the electron self-energy $\Sigma_\nk$ is defined by an integration over $\qq$-points in the full BZ, but symmetries of the little group of $\kk$ can be used to restrict the integration to a smaller zone:
$$
\Sigma_\nk = \int_{BZ} d\qq [...] = \int_{IBZ_\kk} d\qq w^\kk(\qq) [...]
$$
In this case, use:
```
gstore_kzone "ibz"
gstore_qzone "bz"
gstore_use_lgk 1   # Default is 0
```

For phonon properties, use:
```
gstore_kzone "bz"
gstore_qzone "ibz"
gstore_use_lgq 1   # Default is 0
```

!!! important
    Not all e-ph calculations are compatible with little group filtering. Please verify the documentation or run small test calculations before launching large-scale jobs.

Next, you must select the bands that enter the e-ph matrix elements. Several options are available to simplify different types of calculations. By default, if no specific option is provided, ABINIT computes all matrix elements with $m$ and $n$ indices ranging from 1 to [[nband]]. This is rarely ideal, as not all transitions are required for final physical properties. You should explicitly specify the band ranges in the input file.

The variable [[gstore_brange]] defines the range for both $m$ and $n$ indices for each collinear spin channel ([[nsppol]] = 2). It provides full control over which bands to include and is useful when you need all matrix elements connecting states within a specific range, such as in polaron calculations.

For metals or transport properties in semiconductors, the relevant contributions usually come from transitions within an energy window around the Fermi level or band edges. In these cases, it is more efficient to filter bands automatically using an energy range defined by [[gstore_erange]].

!!! important
    [[gstore_erange]] is not compatible with [[gstore_brange]].

Other applications may require filtering $\kk$-points or using different band ranges for the $\psi_\nk$ and $\psi_\mkq$ states. In ZPR calculations of the fundamental band gap, for example, the $\kk$-points and $n$ index can be restricted to the band edges (automatically detected from the KS energies). The $m$ index, however, should cover a much larger range to account for empty states in the summation, while $\qq$-points cover the full BZ or an appropriate irreducible wedge.

The [[gstore_kfilter]] variable allows you to apply additional filtering directly to electronic states. The most appropriate choice depends on the specific physical property being computed; please refer to the [[gstore_kfilter]] documentation for details. Finally, [[gstore_with_vk]] allows you to include electronic group velocities (and optionally off-diagonal velocity matrix elements) in the GSTORE file, which is useful for transport calculations or for maintaining gauge consistency.

The following examples provide recommended settings for different properties. We rely on default behavior whenever appropriate.
For computing the ZPR of the fundamental/direct band gap:
```
gstore_kfilter "qprange"  # Compute g(k,q) only for |nk> at the band edges
gstore_use_lgk 1          # Only q-points in the IBZ_k
gstore_brange 1 12        # Range for the m index (last index cannot exceed nband)
nband         12
```

!!! important
    By "band edges," we refer to the $\kk$-points in the WFK file where the conduction band minimum (CBM) and valence band maximum (VBM) are found. These points do not necessarily coincide with the true band extrema if the latter do not lie on the chosen $\kk$-mesh (e.g., in Silicon). If a more accurate description of the true band edges is required, generate a WFK file using a shifted $\kk$-mesh via [[nshiftk]] and [[shiftk]].

To manually control the list of $\kk$-points and bands for the $|n\kk\rangle$ states, remove the "kfilter" option and use [[nkptgw]], [[kptgw]], and [[bdgw]]:
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
nband           12  # Last index cannot exceed nband
```

## Gstore computation and MPI parallelism
Two EPH sub-drivers can generate a `GSTORE.nc` file: [[eph_task]] = 11 computes matrix elements at the KS level, while [[eph_task]] = 17 uses the more expensive GWPT formalism [[cite:Li2019]] to produce a file containing both GWPT and KS matrix elements.

[[eph_task]] 11 is parallelized over five different MPI levels. You can manually specify the MPI grid using [[eph_np_pqbks]]. The product of MPI processes along all dimensions must equal the total number of processes, as idle processes are not supported. If [[eph_np_pqbks]] is not specified, the code generates the grid automatically at runtime.

When enforcing your own grid with [[eph_np_pqbks]], consider that parallelization over collinear spins, $\kk$-points, and $\qq$-points is most efficient. The number of processes for $\kk$ or $\qq$ should be adjusted based on [[gstore_kzone]] and [[gstore_qzone]]. To reduce load imbalance, use fewer processes for the wavevector restricted to the IBZ. Parallelism over perturbations should be activated only when the first three levels are saturated. Note that band parallelism is not supported in the current version.

For GWPT calculations ([[eph_task]] == 17), the MPI grid is defined by [[gwpt_np_wpqbks]].

!!! important
    Writing the GSTORE file requires a NetCDF library with MPI-IO support.

## Densifying the q-mesh
By default, e-ph matrix elements are computed on the coarse ab initio $\qq$-mesh defined by [[ddb_ngqpt]]. To densify this mesh, use [[eph_ngqpt_fine]]. The fine $\qq$-mesh must be identical to or a submesh of the electronic $\kk$-mesh. Ensure you compute Born effective charges and dynamical quadrupoles to properly describe the long-range part of the DFT potentials and achieve reliable Fourier interpolation. Refer to the [EPH intro tutorial](eph_intro.md) for more details.

## Restarting a GSTORE computation
If a GSTORE calculation is interrupted (e.g., due to a timeout), you can restart it simply by rerunning the same input. Restart capabilities are enabled by default but can be disabled by setting [[eph_restart]] to 0.

## Computing physical properties from a GSTORE file
Reading e-ph matrix elements from a file to compute physical properties is straightforward: use [[getgstore_filepath]] and select the appropriate [[eph_task]] for post-processing.

!!! critical
    Do not use a WFK file different from the one used to generate the GSTORE file. The complex e-ph matrix elements depend on the specific gauge of the wavefunctions in the WFK file.

To compute the ZPR from GSTORE, for example:
```
optdriver 7         # Enter EPH code.
eph_task 24         # SIGMAPH from GSTORE

getgstore_filepath  "teph4zpr_10o_DS1_GSTORE.nc"
```
Other options or files may be required depending on the value of [[eph_task]]. Please consult the documentation and available tutorials.
