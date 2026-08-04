---
description: How to compute self-trapped electron and hole polarons and their hopping paths
authors: VV
---
<!--- This is the source file for this topic. Can be edited. -->

This page summarizes how to compute self-trapped electron and hole polarons with the Variational Polaron Equations
(VarPEq) implemented in the ABINIT EPH driver.

## Introduction

A polaron is a quasiparticle formed when a charge carrier couples to the lattice vibrations of a crystal.
In the weak-coupling regime, the carrier remains delocalized and the interaction mainly renormalizes its band energy and
effective mass; see [[topic:TDepES]] and [[topic:ElPhonInt]].

In the strong-coupling regime, the carrier can create a lattice distortion that traps it in a localized state.
This process is known as **self-trapping** or **autolocalization**.
Its stability is characterized by the polaron binding energy: the total-energy difference between distorted and
undistorted systems containing $N\pm1$ electrons, where $N$ is the electron count of the neutral system.
A negative binding energy indicates a stable self-trapped state.

ABINIT treats this problem with the VarPEq formalism [[cite:Vasilchenko2022]], [[cite:Vasilchenko2025]], closely related
to the reciprocal-space polaron equations of [[cite:Sio2019]].
The binding energy is minimized with respect to coefficients describing the carrier wavefunction in a basis of Bloch
states and the lattice distortion in a basis of phonon modes.
The explicit variational functional and its gradients are presented in [[tutorial:eph4vpq]].

The reciprocal-space sampling defines the Born--von Karman (BvK) supercell hosting the polaron.
For example, a uniform $4\times4\times4$ mesh represents a $4\times4\times4$ supercell.
The usual workflow employs equal, uniform, $\Gamma$-centered $\mathbf{k}$- and $\mathbf{q}$-meshes, although commensurate
meshes are sufficient in principle.

## Calculation workflow

A VarPEq calculation uses the EPH driver ([[optdriver]] = 7) and is normally organized as follows:

1. Perform a ground-state calculation to obtain the density and wavefunctions.

2. Compute phonons and first-order potentials with [[topic:DFPT]].
   For polar materials, also compute the high-frequency dielectric tensor, Born effective charges, and, when required,
   dynamical quadrupoles for the long-range interpolation.

3. Merge the partial DDB and first-order potential files with `mrgddb` and `mrgdv`, respectively.

4. Perform an NSCF calculation on the $\mathbf{k}$-mesh representing the target BvK supercell.
   Include all valence or conduction bands that may participate in the polaron, including degenerate band-edge states.

5. Use [[eph_task]] = 11 to produce a GSTORE.nc file with electron-phonon matrix elements on the target mesh.
   Select the relevant band manifold with [[gstore_brange]] or [[gstore_erange]].

6. Use [[eph_task]] = 13 to optimize the polaron.
   The calculation reads the WFK, DDB, and GSTORE.nc files via [[getwfk_filepath]], [[getddb_filepath]], and
   [[getgstore_filepath]], and writes the solution to VPQ.nc.

7. Use [[eph_task]] = -13 and [[getvpq_filepath]] to reconstruct the charge distribution and atomic displacements as
   XSF files.

!!! important

    A GSTORE for an electron polaron must contain only the selected conduction bands, while a GSTORE for a hole polaron
    must contain only the selected valence bands.
    Mixing the two manifolds leads to incorrect results.

!!! warning

    Producing GSTORE.nc requires parallel I/O support in ABINIT and its HDF5 and NetCDF dependencies.
    GSTORE files may also require substantial memory and disk space.

## Main controls and convergence

The mandatory [[vpq_pkind]] variable selects an `"electron"` or `"hole"` polaron.
The initial carrier distribution is controlled by [[vpq_aseed]], with parameters such as [[vpq_gpr_energy]],
[[vpq_gpr_length]], or [[vpq_atloc]], depending on the selected seed.

The solver stops when the gradient residual reaches [[vpq_tolgrs]] or after [[vpq_nstep]] iterations.
Because the energy landscape can contain several local minima, different initial seeds should be tested.
Additional solutions can be searched for with [[vpq_nstates]]; [[vpq_nstep_ort]] and [[vpq_translate]] control how
previous solutions and their translationally equivalent images are excluded during this search.

An existing VPQ.nc solution can be restarted with [[eph_restart]] or interpolated onto a compatible mesh with
[[vpq_interp]].
This is useful for following the same localized state through increasingly dense meshes.
Use [[vpq_select]] to continue a particular state from a file containing several solutions.

The binding energy and spatial localization must be converged with respect to the ground-state and DFPT parameters, the
electronic and phonon band ranges, the optimization tolerance, and the BvK supercell size.
In polar materials, [[eph_frohl_ntheta]] and [[vpq_avg_g]] can account for the long-range contribution around $\Gamma$
and improve mesh convergence.
Calculations on several supercell sizes can also be extrapolated to the infinite-size limit.

## Output and post-processing

VPQ.nc contains the optimization history, energy decomposition, optimized electronic and vibrational coefficients, and
the metadata required for restart and post-processing.
|AbiPy| can inspect this file, plot convergence histories, analyze band and phonon-mode contributions, and compare results
from different meshes.

With [[eph_task]] = -13, ABINIT generates XSF files containing the polaron charge distribution, distorted structure, and
atomic-displacement vectors.
For large BvK supercells, [[vpq_mesh_fact]] reduces the real-space grid and memory requirement, while [[vpq_trvec]]
translates the solution before visualization.

## Polaron hopping

ABINIT can optimize a minimum-energy path between known polaron states with the simplified string method
[[cite:Weinan2007]], [[cite:Vasilchenko2026]].
Set [[vpq_mode]] to `"hopping"`, provide the initial state with [[vpq_hop_from_filepath]], and either read the final state
from [[vpq_hop_to_filepath]] or generate it by translation with [[vpq_hop_vec]].
In this mode, [[vpq_nstates]] specifies the number of images along the path, while [[vpq_hop_nstep]],
[[vpq_hop_tolgrs]], and [[vpq_hop_ts]] control its optimization.

## Related Input Variables

*compulsory:*

- [[abinit:vpq_pkind]]  Variational Polaron eQuations: Polaron KIND
 
*basic:*

- [[abinit:getvpq]]  GET the VPQ.nc from dataset
- [[abinit:getvpq_filepath]]  GET the VPQ.nc from FILEPATH
- [[abinit:vpq_aseed]]  Variational Polaron eQuations: A_nk-coefficients SEED
- [[abinit:vpq_atloc]]  Variational Polaron eQuations: ATom LOCalization site
- [[abinit:vpq_avg_g]]  Variational Polaron eQuations: AVeraGe matrix-elements at Gamma
- [[abinit:vpq_gpr_energy]]  Variational Polaron eQuations: Gaussian PaRameters -- electronic ENERGY
- [[abinit:vpq_gpr_length]]  Variational Polaron eQuations: Gaussian PaRameters -- localization LENGTH
- [[abinit:vpq_hop_from_filepath]]  VPQ.nc, HOPping FROM: FILEPATH
- [[abinit:vpq_hop_nstep]]  Variational Polaron eQuations, HOPping: Number of iteration STEPs
- [[abinit:vpq_hop_tolgrs]]  Variational Polaron eQuations, HOPping: TOLerance on the Gradient ReSidual
- [[abinit:vpq_hop_ts]]  Variational Polaron eQuations, HOPping: Time Step
- [[abinit:vpq_hop_vec]]  Variational Polaron eQuations: HOPping VECtor
- [[abinit:vpq_interp]]  Variational Polaron eQuations: INTERPolation
- [[abinit:vpq_mesh_fact]]  Variational Polaron eQuations: SCALE MESH for polaron wavefunction
- [[abinit:vpq_mode]]  Variational Polaron eQuations: MODE
- [[abinit:vpq_nstates]]  Variational Polaron eQuations: Number of polaronic STATES
- [[abinit:vpq_nstep]]  Variational Polaron eQuations: Number of iteration STEPs
- [[abinit:vpq_nstep_ort]]  Variational Polaron eQuations: Number of STEPs with ORThogonalisation
- [[abinit:vpq_select]]  Variational Polaron eQuations: SELECT polaronic state
- [[abinit:vpq_tolgrs]]  Variational Polaron eQuations: TOLerance on the Gradient ReSidual
- [[abinit:vpq_translate]]  Variational Polaron eQuations: TRANSLATE solutions
- [[abinit:vpq_trvec]]  Variational Polaron eQuations: TRanslation VECtor
 
*expert:*

- [[abinit:vpq_efilter]]  Variational Polaron eQuations: Energy FILTER
- [[abinit:vpq_hop_from_ip]]  Variational Polaron eQuations: HOPping FROM Ith Polaron state
- [[abinit:vpq_hop_from_site]]  Variational Polaron eQuations: HOPping FROM SITE
- [[abinit:vpq_hop_to_filepath]]  VPQ.nc, HOPping TO: FILEPATH
- [[abinit:vpq_hop_to_ip]]  Variational Polaron eQuations: HOPping TO Ith Polaron state
- [[abinit:vpq_hop_to_site]]  Variational Polaron eQuations: HOPping TO SITE
- [[abinit:vpq_mix_fact]]  Variational Polaron eQuations: MIXing FACTor
 

## Selected Input Files

No input file associated to this topic.

## Tutorials

* [[tutorial:eph4vpq|Self-trapped polarons and Variational Polaron Equations]] presents the full formalism and LiF examples
  for small and large polarons, multiple solutions, finite-size extrapolation, visualization, and hopping paths.
