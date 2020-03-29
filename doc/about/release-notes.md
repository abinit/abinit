## v9.0
  
Version 9.0, released on March 29, 2020.
List of changes with respect to version 8.10.

Many thanks to the contributors to the ABINIT project between
October 2018 and March 2020. These release notes
are relative to modifications/improvements of ABINIT v9.0 with respect to v8.10
(merge requests up to, and including, MR613 are taken into account,
as well as MR636).

The list of contributors includes:
B. Amadon, L. Baguet, J.-M. Beuken, J. Bieder, J. Bouchet, E. Bousquet, F. Bruneval, G. Brunin, Wei Chen, 
J.-B. Charraud, Ph. Ghosez, M. Giantomassi, O. Gingras, X. Gonze, F. Goudreault, 
B. Guster, G. Hautier, Xu He, N. Helbig, F. Jollet,
H. Miranda, F. Naccarato, R. Outerov, G. Petretto, N. Pike, Y. Pouillon, F. Ricci, M. Royo,
M. Schmitt, M. Stengel, M. Torrent, J. Van Bever, M. Verstraete, J. Zwanziger.

It is worth to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases.
This might take some time ...

Xavier

### **A.** Important remarks and warnings.

**A.1** At the occasion of the switch from ABINITv8 to ABINITv9, many improvements of the formats and content of files written
    by ABINIT have been made, so the backward compatibility of ABINITv9 may be broken. 
    The present ABINITv9.0 is NOT to be considered a production version. It is a beta release, allowing developers to get feedback
    from the users. Many features will work correctly, of course. Still, beginners are advised
    to stick to ABINITv8.10.3 except if ABINITv8.10.3 is not appropriate (or not working) for them.

In particular: 

1. The build system relies on new `.ac9` files (see [B.5](#v9.0.B.5)), superceeding the v8 `.ac` files.
   A bash script (`upgrade-build-config-file.sh`) located in the top level directory of the package can be used
   to convert from the old `.ac`format to `.ac9`.
2. The build system of ABINITv9 does not build anymore the hard dependencies (Linalg, NetCDF4, HDF5, LibXC, ...), 
as this was not sustainable (see [B.5](#v9.0.B.5)) and nowadays most users install themselves prerequired libraries.
3. The main ABINIT output file now contains sections written in YAML (sometimes replacing text sections, sometimes adding information).
    This means that some user-developed parsing tools might not work anymore, and should be adapted to the new ABINITv9 output file (see [B.8](#v9.0.B.8)). Note that the YAML output is still under development and modifications may appear in the next versions. A python API to extract the results of the calculation will be provided when the implementation is finalized.
4. Several default values have been changed, see [A.3](#v9.0.A.3).


**A.2** 
A new account of the ABINIT effort has been published [[cite:Gonze2020]], and provides description
of several new features.
A version of this paper, that is not formatted for Computer Phys. Comm., is also 
[available](https://www.abinit.org/sites/default/files/ABINIT20.pdf).
The licence allows the authors to put it on the Web.
Other specific publications are mentioned in the [Suggested acknowledgment page](../theory/acknowledgments).

<a name="v9.0.A.3"></a>
**A.3**  The default values of the following ABINIT input variables have been changed:
    [[ixcrot]], [[chneut]], [[ntime]], [[symsigma]], [[prtkden]].

**A.4** The initialization of the wavefunctions when [[paral_kgb]]=1 and [[nspinor]]=2 has been changed, since the previous one could prevent the code to converge.
    By M Torrent (MR 562).

**A.5** The input variable xangst has been disabled. Use xcart instead, and specify the unit, namely Angstrom.

* * *

### **B.** Most noticeable achievements

**B.1** Electron-phonon interaction (mobilities in the self-energy relaxation time approximation, 
temperature-dependent electronic band structures including the zero-point renormalization, etc.)

The new capabilities of ABINITv9 related to electron-phonon calculations are described 
fully in the Sec. 3.3.2 of [[cite:Gonze2020]], as follows.

>   In abinit v9, it is possible to compute the EPH self-energy
>   in the Kohn–Sham representation using the EPH matrix
>   elements. The code employs optimized algorithms to compute
>   either the full self-energy (needed for QP corrections and spectral
>   functions) or just the imaginary part that is then used to evaluate
>   mobilities within the self-energy relaxation time approximation
>   (SERTA). The computation of the mobility is fully
>   integrated inside abinit, and is an automatic output of the
>   computation of the imaginary part of the self-energy, bypassing
>   the need to post-process results. When computing the full self-energy,
>   it is possible to reduce the number of empty states
>   required for convergence by using the first-order wavefunctions
>   obtained by solving the relevant Sternheimer equation. 
>
>   In the case of lifetime computations, the code takes advantage of the
>   tetrahedron method to filter contributing q-points, a double-grid
>   integration technique to accelerate the convergence at marginal
>   additional computational cost, and samples the relevant regions
>   in the Brillouin zone contributing to transport properties thus
>   leading to a significant reduction of the computational effort.
>   Crystalline symmetries are used throughout the code in order to
>   reduce the number of k- and q-points that must be explicitly
>   included in the integrals. To achieve good parallel efficiently, the
>   most CPU demanding parts are parallelized with MPI employing a
>   distribution schemes over k/q-points, perturbations and bands (the
>   band level is available only when computing the full self-energy).

Moreover, the interpolation of the DFPT potential, described in Sec. 3.3.1 of [[cite:Gonze2020]] is fully operational,
with many tests provided.

List of tests: [[test:v9_50]], [[test:v9_61]] and [[test:v8_44]].

New input variables: [[dvdb_qcache_mb]], 
[[eph_phrange]], [[eph_tols_idelta]], [[eph_ecutosc]], [[eph_restart]], 
[[eph_stern]], [[eph_use_ftinterp]],
[[getdvdb]], [[getdvdb_filepath]], [[getkerange_filepath]],
[[irddvdb]], [[prteliash]], [[sigma_bsum_range]], [[sigma_erange]],
[[sigma_ngkpt]], [[sigma_nshiftk]], [[sigma_shiftk]], [[symv1scf]].

Note that the new EPH processing unit of ABINIT [[optdriver]]=7 has a different implementation than the one implemented in anaddb.
A new set of tutorials are in preparation and they will be made available in the forthcoming versions. 
For further details about the implementation, please consult this [preprint](https://arxiv.org/abs/2002.00630).

By G. Brunin, H. Miranda, M. Giantomassi, G.-M. Rignanese, G. Hautier.


**B.2** Flexoelectricity and dynamical quadrupoles

ABINITv9 allows one to compute flexoelectricity as well as dynamical quadrupoles, as described
in the Sec. V. D of [[cite:Romero2020]], with underlying theory and test calculations
presented in [[cite:Royo2019]]. The user is advised to consult these references.

At the practical level, see [[cite:Romero2020]]: 

>   In this way, both perturbations are generalized to finite q, as
>   is already the case for atomic displacements. This enables us
>   to carry out an analytical third order derivative of the energy
>   with respect to two of the standard perturbations, and to the
>   momentum q, which directly provides the sought-after spatial
>   dispersion tensors. Remarkably, by virtue of the 2n+1 theorem9,
>   the third-order energies are computed in one shot using
>   precalculated first-order response functions to the standard
>   perturbations, without the necessity of self-consistently computing
>   any response function to a perturbation gradient. After
>   execution, the long-wave DFPT routines generate a derivative
>   database that is subsequently used by post-processing tools
>   implemented in ANADDB to compute and print the different
>   contributions to the FxE tensor.

>   The dynamical quadrupoles are the spatial dispersion counterparts of the Born effective charges,
>   and can be used in lattice dynamics calculations
>   to improve the prevalent dipole-dipole treatment of the
>   long-range interactions. The ANADDB routines that carry
>   out the process of interpolating the dynamical matrix following
>   Ref. 34 have been adapted to incorporate the dipolequadrupole
>   and quadrupole-quadrupole electrostatic interactions
>   derived in Ref. 102. This new functionality results in a
>   faster convergence of the phonon bands calculation with respect
>   to the density of q points and, in some materials, represents
>   the only route to obtain the correct sound velocities.

>   Currently, the implementation is restricted to the use of
>   norm-conserving pseudopotentials without non-linear core corrections, and the LDA
>   functional.

A tutorial is in preparation, with tests [[tests:lw_1]] to [[tests:lw_2]].

New input variables have been defined: lw_flexo, lwqdvpl, prepalw, and d3e_pert2_strs. 
At present (v9.0.2) they are not yet documented, although explained and tested in the above-mentioned tutorial.
This will be completed for v9.2.

M. Royo, M. Stengel, M. Giantomassi (MR 618).


**B.2** DFT+DMFT

The new capabilities of ABINITv9 related to DFT+DMFT calculations are described
fully in the Sec. 3.7 of [[cite:Gonze2020]], as follows.

>   The DFT+DMFT parallelism was improved for large
>   systems. In particular, it is now possible to parallelize the calculation
>   on both k-points and bands/g-vectors by using the input
>   variable [[paral_kgb]] = 1 and related input variables.
>
>   Two new approaches to CT-QMC have been added to
>   solve the AIM. In the first one, the density–density CT-QMC code
>   available in abinit [[cite:Gonze2016]], [[cite:Bieder2014]] was generalized in order to take into
>   account off-diagonal elements of the hybridization function. This
>   implementation is activated with the input variable [[dmft_solv]]
>   = 8. Spin–orbit coupling calculations are possible, but using a
>   real valued imaginary time hybridization function. This solver was
>   used in Refs. [[cite:Amadon2015]], [[cite:Amadon2016]].
>
>   In the second approach, we use the Toolbox for Research
>   on Interacting Quantum System (TRIQS)library [[cite:Parcollet2015]], which is an
>   open-source project that provides a framework for many-body
>   quantum physics and more specifically for strongly-correlated
>   electronic systems. TRIQS provides an open source implementation
>   of the continuous-time hybridization expansion quantum
>   impurity solver (CT-HYB) [[cite:Seth2016]], considered a state-of-the art
>   solver for multi-orbital AIM. An interface between abinit and
>   the impurity solver TRIQS/CT-HYB is now available and will make
>   use of the independent progress made by the TRIQS library. <...>

Also, the DMFT k-resolved spectral function is available (MR 529, 490).

List of tests: [[test:paral_84]], [[test:paral_86]], [[test:paral_99]], [[test:v8_01]].
New input variables: [[dmft_charge_prec]] and [[dmft_kspectral_func]] (test to be provided for the latter). 
Also [[dmft_occnd_imag]], but only for keeping backward compatibility for tests.

By T. Cavignac, B. Amadon and O. Gingras.


**B.3** Spin model within Multibinit 

The new capabilities of Multibinit within ABINITv9 are described
fully in the Sec. 4.1 of [[cite:Gonze2020]]. See also Sec. [D.1](#v9.0.D.1).
In particular, a spin model, described specifically in Sec. 4.1.2 of [[cite:Gonze2020]], is available, as follows.

>Multibinit implements the most commonly used model for spin systems, 
>via the Heisenberg Hamiltonian including magnetic exchange and Dzyaloshinskii Moriya interactions. 
>Single ion anisotropy and dipole–dipole interactions are also included, 
>and all terms bear a very strong similarity to the quadratic part of the lattice model Hamiltonian.
>A number of open source spin dynamics codes already exist, such as UPPASD, VAMPIR, OOMF; 
>the distinguishing features of multibinit are the integration with abinit,
>to fit parameters, and the simultaneous dynamics with other
>degrees of freedom (in particular using the inter-atomic force constants). 

A tutorial for the multibinit spin model has been written, [[tutorial:spin_model]].

Many new input variables are present.
Not all these new input variables are present in automatic tests, though, in this beta-release. 
These "non-tested" input variables are indicated below with 'NT'.
This will be completed for the production version v9.2 .
List of tests in addition to those of the tutorial: [[test:v8_16]], [[test:v8_23]], [[test:v9_81]], and [[test:v9_82]].

New input variables tested and documented: 
[[spin_calc_thermo_obs@multibinit|spin_calc_thermo_obs]], 
[[spin_damping@multibinit|spin_damping]], 
[[spin_init_orientation@multibinit|spin_init_orientation]], 
[[spin_init_qpoint@multibinit|spin_init_qpoint]],
[[spin_init_rotate_axis@multibinit|spin_init_rotate_axis]], 
[[spin_init_state@multibinit|spin_init_state]], 
[[spin_ntime_pre@multibinit|spin_ntime_pre]], 
[[spin_projection_qpoint@multibinit|spin_projection_qpoint]], 
[[spin_sia_add@multibinit|spin_sia_add]], 
[[spin_sia_k1amp@multibinit|spin_sia_k1amp]], 
[[spin_sia_k1dir@multibinit|spin_sia_k1dir]], 
[[spin_temperature_start@multibinit|spin_temperature_start]], 
[[spin_temperature_end@multibinit|spin_temperature_end]],
[[spin_temperature_nstep@multibinit|spin_temperature_nstep]], 
[[spin_var_temperature@multibinit|spin_var_temperature]], 
[[spin_write_traj@multibinit|spin_write_traj]].
Additionnal new input variables:
[[slc_coupling@multibinit|slc_coupling]](NT),
[[spin_calc_correlation_obs@multibinit|spin_calc_correlation_obs]] (NT and not documented), 
[[spin_calc_traj_obs@multibinit|spin_calc_traj_obs]] (NT and not documented), 
[[spin_projection_qpoint@multibinit|spin_projection_qpoint]] (NT). 

By Xu He, N. Helbig, J. Bieder, Ph. Ghosez, M. Verstraete


**B.4** Constrained DFT

Constrained Density-Functional Theory (see [[topic:ConstrainedDFT]]) is available,
with a new algorithm allowing to impose the constraints to arbitrary precision,
whether it relates to the charge, magnetization, magnetization direction, or magnetisation size,
or a combination thereof for different atoms. The constraints are smeared spherical integrals
with ajustable sphere radius, centered on atoms. The algorithms has been demonstrated for norm-conserving pseudopotentials
as well as PAW. Forces and derivatives with respect to the constraints 
are available (i.e. magnetic torque for the non-collinear spin case).
Stresses are still to be coded, will be available in ABINITv9.2.

New tests: v8#24-29, v8#95-97 and v9#1-3.
New input variables: [[chrgat]], [[constraint_kind]], [[ratsm]].

By X. Gonze.

<a name="v9.0.B.5"></a>
**B.5** Large modifications of the build system 

The build system relies on new <hostname>.ac9 files, superceeding the v8 <hostname>.ac files.
Fully documented example files can be found in doc/build/config-examples.
The build system of ABINITv9 does not build anymore the (hard and soft) dependencies (Linalg, NetCDF4, HDF, LibXC, Wannier90, ...), as this was not sustainable.
Three libraries are now mandatory: linalg, NetCDF4/HDF5 and LibXC. Failing to link to them will prevent building ABINIT.
The other libraries are optional, there will only be a warning if they are not available.
If the user does not provide the path to these libraries,
the build system will try to find them in the "usual" directories, and inform the user that it has done so.
The build system also can make suggestions to the user, to complete its *.ac9 file.

By Y. Pouillon and JM Beuken

**B.6** New command line interface

There is a new (**recommended**) command line interface to run ABINIT, without the "files" file.
The new syntax is:

    abinit run.abi

or

    abinit run.abi > run.log 2> run.err &

where `run.abi` is the Abinit input file that now provides all the information related to pseudos 
and the prefixes that were previously passed via the "files" file. For comparison, the old syntax is

    abinit < run.files > run.log 2> run.err &      ! This is the old syntax

A file extension for the input file is highly recommended (in this example we use `.abi`)
as by default the parser will use the string before the file extension as root to build the prefixes 
for the input/output/temporary files.

The user can specify the name of the main output file thanks to the [[output_file]] input variable,
the list of pseudopotentials thanks to the [[pseudos]] input variable and the directory where 
all pseudos are located with [[pp_dirpath]].
The prefix for other input, output or temporary files can be specified with [[indata_prefix]], [[outdata_prefix]] and 
[[tmpdata_prefix]], respectively. 
A default set of prefixes computed from the basename of the input file is used if
these variables are not specified in the input.

For some examples, see tests [[test:v8_90]], [[test:v7_45]], and [[test:v5_54]]. See also [[topic:Control]]. 

A similar command line interface can also be used for the anaddb code.
In this case, the relevant variables are: 
[[output_file@anaddb]], [[ddb_filepath@anaddb]], [[ddk_filepath@anaddb]], [[gkk_filepath@anaddb]], [[eph_prefix@anaddb]].

The new syntax is:

    anaddb run.in > run.log 2> run.err &

See tests [[test:v8_52]] for a standard analysis of the DDB file and 
[[test:v7_94]] for the (old implementation) of electron-phonon calculations in anaddb.

!!! important

    The old "files file" interface is still operational although deprecated and will be **REMOVED** in Abinit v10.


By M. Giantomassi (MR 586).


**B.7** Reading strings from the input file

A new mechanism to read strings enclosed between **double quotation marks** from the input file has been activated. 
So, many new input keywords are reading strings as data, and, often, can be used alternatively to similar input keywords
that were expecting numerical values such as the `get*` and `ird*` variables.
The goal is to encourange a new approach for performing ABINIT calculations in which multiple datasets 
and `get*` variables are replaced by indipendent input files that are connected together via file paths.

List of new input variables that rely on this feature:

- [[getddb_filepath]], an alternative to [[getddb]] or [[irdddb]], see test [[test:v9_60]]
- [[getden_filepath]], an alternative to [[getden]] or [[irdden]], see test [[test:v8_36]] and [[test:v8_41]]
- [[getscr_filepath]], an alternative to [[getscr]] or [[irdscr]], see test [[test:v67mbpt_51]]
- [[getwfkfine_filepath]], an alternative to [[getwfkfine]] or [[irdwfkfine]], see tests [[test:v9_55]], [[test:v9_56]]
- [[getwfk_filepath]], an alternative to [[getwfk]] or [[irdwfk]], see test [[test:v9_60]]
- [[getwfq_filepath]], an alternative to [[getwfq]] or [[irdwfq]], see test [[test:v7_98]]
- [[getkerange_filepath]], see test [[test:v9_60]]
- [[getpot_filepath]], see test [[test:v8_44]]
- [[pseudos]], [[indata_prefix]], [[outdata_prefix]], [[tmpdata_prefix]], [[output_file]], see test [[test:v9_04]]
- [[pp_dirpath]]: cannot be tested  EXPLICITLY because system dependent but used by runtests.py when generating the input file. 
- [[output_file@anaddb]], [[ddb_filepath@anaddb]], [[gkk_filepath@anaddb]], [[eph_prefix@anaddb]].
  See tests [[test:v8_52]] for a standard analysis of the DDB file and 
  [[test:v7_94]] for the (old implementation) of electron-phonon calculations in anaddb.

By M. Giantomassi

<a name="v9.0.B.8"></a>
**B.8** YAML sections in the output file 

YAML sections are now generated in the output file, sometimes replacing text sections, sometime providing new information.
At present there is a YAML section for the components of the total energy, the GS results including forces and stresses as well as a YAML section for GW calculations, and some YAML sections giving information about the iteration status.
<!--
Example of tests: paral#86, v67mbpt#2. See the input variable use_yaml (TO BE DOCUMENTED).
-->
At the occasion of the development of this capability, and its adaptation to the test farm, the
perl script fldiff.pl has been replaced by a Python version.
See related information in Sec. 5.5 of [[cite:Gonze2020]].

By T. Cavignac, M. Giantomassi, GM Rignanese, X Gonze.


**B.9** New approach to define crystalline structures in the Abinit input

The new variable [[structure]] can be used to initialize the lattice vectors 
and the atomic positions **from an external file**.
Variables such as [[natom]], [[ntypat]], [[typat]] and [[znucl]] are automatically initialized
and need not to be specified in the ABINIT input.
At present, the code can read netcdf files produced by ABINIT (`GSR.nc`, `WFK.nc`, `DEN.nc`, `HIST.nc`)
and POSCAR files in VASP-5 format.
See the documentation for the syntax and limitations.

By M. Giantomassi.


**B.10** New capabilities of abipy and abiflows 

The abipy and abiflows projects have been significantly extended.
See Sec. 6 of [[cite:Gonze2020]], as well as the [gallery of plotting scripts](http://abinit.github.io/abipy/gallery/index.html) &nbsp;
and the [gallery of abipy workflows](http://abinit.github.io/abipy/flow_gallery/index.html) &nbsp;.

By M. Giantomassi, G. Petretto, F. Naccarato.


* * *

### **C.** Changes for the developers (also compilers)

**C.1** A python script to help ABINIT developers and development.

The new python script abisrc.py located in the top directory of the ABINIT package has been developed.
It has superceded abilint.py in the makemake procedure. 

Try

    ./abisrc.py --help

then follow the suggestions, to get info about files, directories, interfaces, to visualize
dependencies of the ABINIT subroutines, etc.

Note that there are dependencies of abisrc.py, to be installed prior being able to use some of its capabilities.
Use:

    pip install -r requirements.txt --user

to install the dependencies in user mode.

By M. Giantomassi.

**C.2** New characteristics of input variables

In the description of input variables (e.g. files abimkdocs/variables_abinit.py), a new field `added_in_version` has been introduced,
for example,

     added_in_version="9.0.0"

or, for variables introduced prior to v9,

     added_in_version="before_v9"

**C.2** Test farm: new and obsolete bots

* Bots introduced in the test farm: atlas_gnu_9.1_openmpi, buda2_gnu_8.2_cuda, cronos2_gnu_7.4_paral.
* Bots upgraded: abiref_gnu_5.3_* to abiref_gnu_9.2_* ; abiref_nag_6.2_openmpi to
    abiref_nag_7.0_openmpi; bob_gnu_5.3_openmp to bob_gnu_7.5_openmp; buda2_gnu_8.1_mpich3 
    to buda2_gnu_8.2_mpich3; graphene_gnu_6.4_macports to graphene_gnu_9.2_macports; max2_gnu_5.3_* to max2_gnu_6.5_*; ubu_gnu_5.3_openmpi to ubu_gnu_9.2_openmpi.
* Bots removed: atlas_gnu_7.2_fb.ac (no nightly test of the fallbacks anymore), cronos_gnu_5.3_paral(remplaced by cronos2), inca_gnu_6.3_py3k (inca too old), tikal_gnu_* (tikal too old).
Working also on a cronos-cronos2 cluster.

By JM Beuken

**C.3** Supported compilers

* gfort (GNU) compiler: v9 newly supported, v4 obsolete
* ifort (INTEL) compiler: v19 newly supported.
* NAG 7.0 instead of 6.2

By JM Beuken

**C.4** Unitary ttransposer#1 . Test of the transposer for linear algebra to KGB parallelisation.

By J. Bieder.

**C.5** Linkchecker has been reenabled, only for internal link checking.

By JM Beuken (MR 513).

**C.6** Enable the generation of a HTML side-by-side diff on the test farm when fldiff fails with a line count error and it was not caused by a crash of Abinit. The diff algorithms uses a specialized heuristic to improve line synchronization and prevent weird matching.

By Th. Cavignac (MR 526)

**C.7** Split of the source tree (ongoing).

In order to improve modularity, the source tree must be split in two parts, one for low-level routines, largely independent of ABINIT,
and one for more specific routines to ABINIT. The low-level routines should become a separate library, with its own build system and make. 
At present the low-level library have been moved out of src, inside the shared/common/src directory.
See related information in Sec. 5.4 of [[cite:Gonze2020]].

**C.8** New FFT specifications for the build system
See https://gitlab.abinit.org/pouillon/abinit/-/issues/33 .

By Y. Pouillon (MR 619)

* * *

<a name="v9.0.D.1"></a>
### **D.**  Other changes (or on-going developments, not yet finalized)

**D.1**  Miscellaneous improvements of Multibinit (lattice part)

Miscellaneous improvements have been made to the lattice part of Multibinit.
See the new input variables below, also see the Sec. 4.1.1 of [[cite:Gonze2020]].

New tests: [[test:v8_38]], [[test:v8_94]], [[test:v8_98]], [[test:v8_99]], [[test:v9_83]],  [[test:v9_84]], and [[test:v9_85]]. 
New input variables are listed below.
Not all these new input variables are present in automatic tests, though, in this beta-release.
These "non-tested" input variables are indicated below with 'NT'.
Also, not all of these are documented, or only partly documented (e.g. variable type, acronym, default, but no more).
This will be completed for the production version v9.2 .

- [[analyze_anh_pot@multibinit|analyze_anh_pot]] 
- fit_SPC_maxS@multibinit (NT)
- fit_iatom@multibinit [[test:paral_81]], [[test:paral_82]], [[test:v8_13]], [[test:v8_14]], but not documented. 
- latt_compressibility@multibinit, NOT TESTED, NOT DOCUMENTED
- [[latt_friction@multibinit|latt_friction]] 
- latt_mask@multibinit, NOT TESTED, NOT DOCUMENTED
- [[latt_taup@multibinit|latt_taup]] Despite the link, NOT TESTED, only partly DOCUMENTED
- [[latt_taut@multibinit|latt_taut]] Despite the link only partly DOCUMENTED
- [[opt_coeff@multibinit|opt_coeff]]
- [[opt_effpot@multibinit|opt_effpot]]
- [[opt_ncoeff@multibinit|opt_ncoeff]]
- prt_names, NOT TESTED, NOT DOCUMENTED
- [[test_effpot@multibinit|test_effpot]] v8#98

By M. Schmitt, Xu He, F. Ricci, M. Verstraete, Ph. Ghosez


**D.2** Miscellaneous improvements in the Chern number and orbital magnetization calculations,
including parallelization over k points of the Chern number calculation.
New test [[test:v8_39]].

By J. Zwanziger (MR 469, 500, 545, 588)

**D.3** Calculation of Debye-Waller tensor. New test [[test:v8_58]].

By M. Giantomassi


**D.4** Test of linear electro-optical coefficient, tutorespfn toptic#5.

By N. Pike (MR 575).

**D.5** NCPP Wavefunction mixing with Variational Energy
and minor improvements to prepare PAW+Hybrid variational energy.
New test [[test:v7_73]], simple system for testing Hartree-Fock and the SCF algorithms.

By X. Gonze (MR 434, 444, 445).

**D.6** New weight distribution of the Fourier components of the force constants.
Test tolerance in the new integration weights, tests [[test:v8_52]], [[test:v8_53]], [[test:v8_54]].

By H. Miranda and M. Giantomassi

**D.7** Test calculation of velocity matrix elements (DDK) with
 optdriver 8 and [[wfk_task]] "wfk_ddk”, see [[test:v8_59]].

By M. Giantomassi

**D.8** Upgraded [[tutorial:paral_gspw]], new version of auto paral (with threads)

By M. Torrent (MR502).

**D.9** Test wannier90 interface with [[nsppol]]=2 and [[nspden]]=2, [[test:wannier90_04]].

By Xu He

**D.10** Mixed precision for FFT transforms. New input variable [[mixprec]]
see [[test:v8_44]], [[test:v9_57]], [[test:v9_60]], and [[test:v9_61]].

From M. Giantomassi (MR491).

**D.11** Multibinit has been interfaced with scale-up,
https://www.secondprinciples.unican.es

By Marcus Schmitt, Jordan Bieder, Matthieu Verstraete and Philippe Ghosez

**D.12** The following units are now also allowed in input files:

- S Sec Second  for the ABINIT input file;
- nm (for nanometer)  got the ABINIT and ANADDB input files.

**D.13** TDEP utility:
added [[guide:tdep|A-TDep user guide]],
[[topic:Tdep|TDep topic]], and corresponding input variable documentation.
References: [[pdf:TDEP_Paper|TDEP paper]].
Also, see Sec. 4.2 of [[cite:Gonze2020]].

By F. Bottin, J. Bouchet, J. Bieder (MR491,422).

**D.14** Improvements of NLO calculations.
Optimize memory allocation. In particular for [[usepead]] = 1.
Write Raman susceptibilities to netcdf in anaddb.

By G. Petretto  (MR 599).

**D.15** MetaGGA + PAW in progress.

By M. Torrent, JB Charraud (MR 587, 558).

**D.16** Implementation of an alternate version of MPI reductions valid when the number of data exceeds 2^32.

By M. Torrent (MR 587)

**D.17** Build system and src support for llvm (clang/flang) and ARM (armflang/armclang/armpl).

By M Torrent (MR 571)

**D.18** Improvement of core WF reading (optics + electron-positron)

By M Torrent (MR 557)

**D.19** Fix bootstrap kernel convergence issues

By Wei Chen (MR 546)

**D.20** Upgrade versions of libraries used with ABINIT.
Upgrade atompaw to 4.1.0.6. Upgrade Libxc to 4+.

By M. Torrent, JM Beuken (MR 532, 470, 465, 441)

**D.21** Write yaml file for fatbands (phonopy format) with TDEP

By J. Bieder (MR510)

**D.22** Write polarization vectors in GSR.nc.

By H. Miranda (MR 462).

**D.23** Updated user interface of Raman_spec.py . Added one section to tutorial NLO to describe use of Raman_spec.py.

By N. Pike (MR 460, MR 581)

**D.24** XML format for core wavefunctions

By F. Jollet (MR 423)

**D.25** Wavefunction prediction for molecular dynamics.

By F. Jollet (MR 412)

**D.26** Added a preview for the toptic_4.in files in the optic tutorial.

By F. Goudreault (MR 408)

**D.27** Improve ELPA detection; make abinit compatible with ELPA2019

By M. Torrent (MR 626)

**D.28** Upgrade of [[tutorial:base3]] and [[tutorial:basepar]].

By X. Gonze (MR628)

**D.29** New input variable [[prtprocar]], see test [[test:v5_40]].

By M. Verstraete (MR630)

**D.30** The ABINIT input variable [[supercell_latt]] is documented. 
See also [[test:v8_94]].

By X. Gonze

**D.31** New implementation of the cRPA (to compute U and J in DFT+U or DFT+DMFT) 
that can consider any subset of orbitals in the system.
Compilation was a bit long on certain compilers, so new keywords have been added to the build system:
enable_crpa_optim and enable_crpa_no_optim

By R. Outerov and B. Amadon (MR622).

**D.32** Work on mGGA+PAW.
All internal and unit tests are OK. The implementation seems close to the end. Still need to perform physical tests.

By M. Torrent and J.-B. Charraud (MR 625).

**D.33** Work on refactoring the Coulomb interaction part of ABINIT.

By B. Guster, M. GIantomassi and X. Gonze (MR 627&633).

**D.34** Miscellaneous additional bug fixes and improvements of documentation.
L. Baguet, JM Beuken, J. Bieder, E. Bousquet, F. Bruneval, T. Cavignac, M. Giantomassi, X. Gonze, F. Jollet, N. Pike, Y Pouillon, M. Torrent, J. Van Bever, M. Verstraete, Xu He.


* * *

## v8.10

Version 8.10, released on October 15, 2018.
List of changes with respect to version 8.8.

Many thanks to the contributors to the ABINIT project between
April 2018 and October 2018. These release notes
are relative to modifications/improvements of ABINIT v8.10 with respect to v8.8.
The merge request #408 is the first MR not reported in these release notes. Then, #410-411, #413-416 have also been included.

The list of contributors includes:
B. Amadon, G. Antonius, L. Baguet, J.-M. Beuken, J. Bieder, E. Bousquet, F. Bruneval, Wei Chen, M. Cote, 
J. Denier, G. Geneste, Ph. Ghosez, M. Giantomassi, O. Gingras, X. Gonze, F. Goudreault, B. Guster, Xu He, Y. Jia, F. Jollet, 
A. Lherbier, A. Martin, H. Miranda, F. Naccarato, G. Petretto, N. Pike,
S. Ponce, Y. Pouillon, S. Prokhorenko, F. Ricci, M. Torrent, M. van Setten, B. Van Troeye, M. Verstraete, J. Zwanziger.

It is worth to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases.
This might take some time ...

Xavier

### A. Warnings and important remarks

A.1 The correct definition of the temperature has been implemented in the isokinetic algorithm [[ionmov]]=12.

A.2 The multibinit input variable "supercell" has been renamed "supercell_latt".

A.3 Tests v8#100-107 have been moved to v8#81-88.
    Tests paral#100-103 have been moved to paral#80-83

* * *

### B. Most noticeable achievements

B.1 The computation of the Raman intensity in DFPT with PAW is now possible (it was only available with norm conserving psps previously).
    This is based on the second-order Sternheimer equation for the derivative
    with respect to an electric field.
    See tests [[test:v8_81]] to [[test:v8_89]].  
    By L. Baguet and M. Torrent.

B.2 There has been a large effort to clean the documentation, that was needed and made possible thanks to the recent move to the mkdocs system.
    In particular, many new hyperlinks have been created in dozens of markdown files, 
    references have been centralized in the [[theory:bibliography]].  
    By B. Amadon, G. Antonius, L. Baguet, J. Bieder, E. Bousquet, F. Bruneval, Wei Chen, M. Cote, G. Geneste,
    M. Giantomassi, X. Gonze, Xu He, F. Jollet, A. Lherbier, H. Miranda, F. Naccarato, G. Petretto, N. Pike, 
    S. Ponce, Y. Pouillon, M. Torrent, M. van Setten, B. Van Troeye, M. Verstraete, J. Zwanziger.

B.3 The multibinit application (for second-principles calculations) has considerably progressed.
    Documentation has thus been set up: "topics" have been written, as well as a tutorial,
    in addition to the already existing input variable documentation and test cases.
    See [[topic:LatticeModel]], [[topic:BoundingProcess]], [[topic:FitProcess]] and [[topic:DynamicsMultibinit]],
    that are hub to the relevant tutorial, input variables and test cases 
    (e.g. [[lesson:lattice_model]],[[test:v8_15]], [[test:v8_16]]...).  
    By A. Martin, in collaboration with Fabio Ricci and Ph. Ghosez

B.4 Several new options are available for the [[ionmov]] input variable governing ionic dynamic or geometry optimization:

* [[ionmov]]=15 for the FIRE algorithm, [[test:v8_17]];
* for [[ionmov]]=12, isokinetic ensemble, the fixing of atomic positions is now allowed, [[test:v8_21]] and [[test:v8_22]].

Also, the documentation for the hybrid Monte Carlo algorithm has been improved, see [[test:v8_34]], and the new input variables [[hmcsst]] and [[hmctt]].

By Xu He, S. Prokhorenko and X. Gonze.

B.5 The linear combination of images is now allowed, with the new value for input variable [[imgmov]]=6,
    and mixing factor given by [[mixesimgf]].
    In this case, the total energy and forces are also assembled as a linear combination, and the geometry is optimized using 
    algorithms selected with the usual [[ionmov]] input variable.
    See test [[test:v8_20]].
    The wavefunctions from the previous itimimage value (see [[ntimimage]] input variables) can be stored,
    using the new input variable [[imgwfstor]]. This allows saving CPU time at the expense of memory, in all
    the image based algorithms.  
By X. Gonze, testing by Y. Jia.

B.6 Tutorial [[tutorial:nuc]]
    has now a section for the computation of the isomer shift (Mossbauer spectroscopy) based on Fermi contact interaction.  
By J. Zwanziger.

B.7 The Frohlich model is now implemented in the electron-phonon part of ABINIT, [[optdriver]]=7.
    The Frohlich average of effective masses is computed with the DFPT computation of effective masses, see [[test:v8_56]].
    Also, the zero-point renormalization of the band extrema is computed using a general formula valid for isotropic and anisotropic
    solids, as well as for non-degenerate or degenerate extrema, see [[eph_frohlichm]] and [[test:v8_57]].  
By X. Gonze.

* * *

### C. Changes for the developers (also compilers)

C.1 All F90 ABINIT sources are now inside modules.
    Makemake now aborts if F90 procedures outside modules.  
    By M. Giantomassi, with some help by J. Zwanziger, M. Torrent and B. Amadon.

C.2 Prepared the removal of the bindings subsystem.  
    By Y. Pouillon and M. Torrent.

C.2 New [Howto for developers](../developers/developers_howto) (variables, mkparents, robodoc, test_suite).  
    By M. Giantomassi.

C.3 New [Howto for the test suite](../developers/testsuite_howto).    
    By M. Giantomassi. 

* * *

### D.  Other changes (or on-going developments, not yet finalized)

D.1 New input variable [[prtkbff]].
    This input variable activates the output of the Kleynman-Bylander form factors in the netcdf WFK file produced at the end of the ground-state calculation. 
    The form factors are needed to compute the matrix elements of the commutator [Vnl, r] of the non-local part of the (NC) pseudopotentials. 
    This WFK file can therefore be used to perform optical and/or many-body calculations with external codes such as DP/EXC and Yambo. The option is ignored if PAW.  
By H. Miranda and M. Giantomassi.

D.2. A new routine to calculate and write the DDK matrix elements in a EVK.nc files for norm-conserving pseudo-potentials (and nspinor == 1),
    from a WFK file using the commutator [H, r] used to calculate the matrix elements for chi.
    These files are in the same format and contain the same information as the EVK.nc files (former DDK.nc) produced by the DFPT part of the code. 
    This routine is parallelized over k-points and bands and requires a much smaller memory footprint than the DFPT code.  
    By H. Miranda and M. Giantomassi.

D.3 New input variable [[slk_rankpp]].
    This variable controls how the number of processes to be used in Scalapack diagonalization algorithm: [[np_slk]] will be calculated according to this value.  
    By J. Bieder.

D.4 New input variables [[prtefmas]], [[irdefmas]] and [[getefmas]], to deal with the effective masses, e.g. to allow feeding computed effective masses
    from one dataset to another one.   
    By X. Gonze.

D.5 The spin dynamics has been implemented in multibinit.   
    By He Xu.

D.6 New value for input variable [[usepawu]]=4.
    The FLL double counting is used. However, and in comparison to usepaw=1, the calculation is done without polarization in the exchange correlation functional.  
By B. Amadon.

D.7 New extensive testing of the DFPT+PAW+GGA, see [[test:v8_51]].
    However, [[pawxcdev]]=0 is still needed.  
    By M. Torrent.

D.8 New input variable [[prtfull1wf]] to print the full perturbed wavefunctions.   
    By G. Antonius.

D.9 Allows one to suppress completion of the 2DTE for vanishing elements using symmetries only, 
    negative values of [[rfmeth]].  
    By X. Gonze.
    
D.10 Interface with TRIQS_2.0.  
    By O. Gingras.

D.11 Allows to define different occupations for different spins, when [[occopt]]=0 and [[nsppol]]=2.
     By X. Gonze.

D.12 NEB is now tested (at last), [[test:v6_29]].   
     By X. Gonze.

D.13 Many |abitutorials| have been improved.   
     By X. Gonze.
 
D.14 Work on orbital magnetization.   
     By J. Zwanziger.

D.15 Elastic and piezoelectric tensors in netcdf format.  
     By M. Giantomassi

D.16 Allow [[boxcutmin]] for DFPT case, after check by H. Miranda.

D.17 ABIWAN.nc files with Wannier output that can be analyzed wit AbiPy. 
     See the |AbiwanFileNb| for futher information.  
     By M. Giantomassi.

D.18 The [[topic:Macroave]] has been created.  
     By X. Gonze.

D.19. Include CTQMC Code using non diagonal hybridization function. See [[test:paral_83]].  
     By B. Amadon, J. Denier and J. Bieder. 

D.20. Finalize scalapack/elpa integration in lobpcg.
      Use ONLY if no openmp. Warning: scalapack is not thread-safe in general.  
      By J. Bieder.

D.21 Automatic test [[test:v8_37]] for TDep application.  
    By J. Bieder.

D.22 Python scripts to calculate physical properties that can be derived from the elastic tensor 
     (which is the result of an anaddb calculation).  
     By N. Pike.

D.23 Miscellaneous additional bug fixes and improvements of documentation.  
     By B. Amadon, G. Antonius, L Baguet, J.-M. Beuken, J. Bieder, F. Goudreault, F. Jollet, H. Miranda, 
     F. Nacaratto, N. Pike, M. Torrent, M. Verstraete.

* * *

## v8.8

Version 8.8, released on April 28, 2018.
List of changes with respect to version 8.6.

Many thanks to the contributors to the ABINIT project between
November 2017 and April 2018. These release notes
are relative to modifications/improvements of ABINIT v8.8 with respect to v8.6.
The merge request #285 is the first MR not reported in these release notes.

The list of contributors includes:
B. Amadon, G. Antonius, L. Baguet, J.-M. Beuken, J. Bieder, F. Bottin, Y. Bouchet, E. Bousquet, W. Chen, 
C. Espejo, Ph. Ghosez, M. Giantomassi, X. Gonze, F. Jollet, A. Martin,
H. Miranda, G. Petretto, N. Pike, Y. Pouillon, S. Prokhorenko, F. Ricci, 
G.-M. Rignanese, M. Torrent, M. Verstraete, J. Zwanziger

It is worth to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases ...
This might take some time ...

Xavier

### A. Warnings and important remarks

A.1 Due to the availability of new input variables, some obsolete input variables have been suppressed:

* input variable cgtyphf, see now [[fockoptmix]];
* input variable gwfockmix, see now [[hyb_mixing]] and [[hyb_mixing_sr]].

A.2 The algorithm used for computing the weights for the phonon band structure interpolation in ANADDB
    in the case [[brav@anaddb]] = 1 has changed. See B.5.

A.3 Tests v7#71-72 have been moved to v7#76-77.

A.4 Replaced *tribes* by *relevance* in doc/topics.

A.5 Replace `EELF` file extension by `ELF`. See [[prtelf]].
    By Guido Petretto.

A.6 Definition of a maximal value for dilatmx, at 1.15, than can be bypassed by setting chkdilatmx=0.
    This is to prevent users slowing down ABINIT too much inadvertantly.

* * *

### B. Most noticeable achievements

B.1 The whole ABINIT documentation has been placed under the control of [mkdocs](http://www.mkdocs.org/),
    and most files have been translated to markdown 
    (well, there are still a few remaining files not placed in this system, but these are quite few).
    The capabilities developed for [v8.6](#v86)
    (Topics -B.1-, central bibliography -B.2-, restructured and searcheable list of input variables -B.3-,
    frequency statistics -B.4-) have been maintained and consolidated.
    More documentation is available for developers than before. The whole system is better integrated and easier to maintain.
    The appearance is also new. 
    The work on documentation is nearly complete, still not all bibliographical references of the doc have been entered in this
    central bibliographic database.
    Entry point: see the new header of any ABINIT documentation file e.g. the [new user's guide](..).
    By M. Giantomassi, with some help from X. Gonze.

B.2 The DFPT has been extended to non-collinear systems ([[nspden]] = 4), with or without spin-orbit coupling,
    for the ddk, electric field and atomic displacement perturbations,
    as well as for the Zeeman magnetic field perturbation (see B.3). 
    See tests from [[test:v8_66]] to [[test:v8_80]]. For experts, see the new input variable [[ixcrot]].
    By F. Ricci, S. Prokhorenko, M. Verstraete, M. Torrent and E. Bousquet.
    
B.3 DFPT can now treat the magnetic field perturbation (Zeeman interaction - magnetic field couple to the spin).
    See the input variable [[rfmagn]], as well as tests [[test:v8_66]] to [[test:v8_70]].
    The new input variable [[tim1rev]] has been introduced, to allow treating perturbations with non-zero q wavevectors.
    By S. Prokhorenko and E. Bousquet.

B.4 The python library [AbiPy](https://github.com/abinit/abipy), for launching ABINIT (also in high-throughput mode)
    and the jupyter notebook based tutorials, nicknamed [abitutorials](https://github.com/abinit/abitutorials) are now sufficiently 
    mature to be advertised. They have been used at the ICTP electron-phonon doctoral school in March 2018.
    Feedback on AbiPy and abitutorials is welcome. 
    By M. Giantomassi.

B.5 A new algorithm (Wigner-Seitz cell based) for computing the weights for the phonon band structure interpolation in ANADDB
    has been implemented. It has replaced the old algorithm in case [[brav@anaddb]] = 1. 
    The old algorithm is still available for back-compatibility purposes, now corresponding to [[brav@anaddb]] = -1, 
    see [[test:v7_93]], although there is no real reason for using it. 
    The new algorithm is very general, respect better the symmetries, and should even supercede 
    the use of other values of [[brav@anaddb]].
    By G. Petretto following discussions with GM Rignanese and XGonze, and tests by Henrique Pereira Miranda.

B.6 The Chern number can be computed, in the norm-conserving case as well as in the PAW case.
    See the theory in [[cite:Ceresoli2006]].
    Associated input variable: [[orbmag]]. Associated test [[test:v8_33]].
    Nuclear magnetic dipole moment code has been improved for efficiency. In particular,
    this improvement is due to converted nucdipmom_k to complex type and explicit BLAS call. 
    Tutorial [[tutorial:nuc|nuc]] is nightly tested.
    By J. Zwanziger ([[tutorial:nuc|nuc]] testing by X. Gonze).

* * *

###C. Changes for the developers (also compilers)

C.1 Add support for NAG 6.2, use netcdf4 and hdf5 with NAG 6.2.
    New builders: abiref_nag_6.2_openmpi and atlas_gnu_7.2_fb

C.2 New version of ELPA module (2011 -> 2017 compatible)
    By M. Torrent

C.3 Replaced http://www.abinit.org by https://www.abinit.org everywhere in the doc. 
    By JM Beuken

C.4 Added fake -n/--no-split option to makemake. 
    This little change will let the test farm keep on working with all branches while the source tree is split into common + core.
    By Y. Pouillon

C.5 Upgrade Abinit to PSML API 1.1.
    By Y. Pouillon

C.6 Allowed for free-form link flags statements in configure options
    By Y. Pouillon

C.7 Intel 18.0 is now officially supported (one bot using it is present in the test farm).
    The situation is not perfect though, as  in mrgscr and mrgdv, the ADVANCE='NO' specification for the write instruction does not work, but
   it works in simple programs. Thus all tests (and their chain of tests) that rely on mrgscr and mrgdv have been disabled for this compiler.
   Namely, v3#87-91, vv67mbpt#37-39, v7#86-88, v8#41-44, v8#63
   Also, the reading of WFK files using MPIIO is not correct, for tests mpiio#26 and 62.

* * *

###D.  Other changes (or on-going developments, not yet finalized)

D.1 Implementation of the LDA-1/2 methodology (see the announcement B.10 of v8.6): [[test:v8_32]] has been provided.
    By F. Jollet.

D.2 Numerous progresses have been made related to the hybrid functionals (although hybrid functionals
    are not yet in production).

* Stresses are now computed correctly.

* The downsampling of the Brillouin Zone to build the Fock operator has been implemented and tested.
See the input variable [[fockdownsampling]] as well as tests [[test:libxc_72]] and [[test:paral_09]]. 

* The B3LYP functional has been implemented, [[ixc]] = -402.

* The new input variables [[hyb_mixing]], [[hyb_mixing_sr]], [[hyb_range_dft]], and [[hyb_range_fock]]
give sufficient flexibility in the PBE0 and HSE family of functionals. 

* GW calculations can now start on top of hybrid functional calculations.

* At variance, there is also now more flexibility to run hybrid calculations using the GW infrastructure 
([[gwcalctyp]] = 5, 15, 25) by the definition of the [[ixc_sigma]] input variable.

* There has been also important work concerning the self-consistency, although this work is not finalized yet
(one reason why hybrid functionals are not yet in production).
The self-consistency at fixed ACE operator can take advantage of an auxiliary XC functional to decrease
the number of inner iterations, see [[fockoptmix]], [[auxc_ixc]] and [[auxc_scal]], while for the outer loop,
in which the ACE operator is upgraded, the wavefunction mixing has been implemented ([[fockoptmix]] and [[wfmix]]).

See the new tests v7#67-72 libxc#44, 45, 72, 73, 74, 
and also the updated tests v4#86, 87, v67mbpt#09, v7#65, libxc#41, 42, 43, paral#09.
By X. Gonze and F. Jollet, with help by M. Torrent.

D.3 The [[tutorial:tdepes|tutorial on temperature-dependence of the electronic structure]] has been upgraded, and carefully tested.
    See all tests in `tutorespfn/tdepes*`.
    By X. Gonze and M. Giantomassi

D.4 Output of interpolated density in the MPI-IO case is now tested, [[test:mpiio_26]] and [[test:mpiio_27]].

D.5 Ongoing work on the multibinit project.
    New input variables fit_nfixcoeff, fit_fixcoeff, fix_generateTerm, 
    see [[test:v8_13]] and [[test:v8_14]].
    New input variable dipdip_prt, see [[test:v8_06]], as well as tests [[test:paral_96]] to paral[102].
    New generator for the polynomial coefficients, debug strain for the fit process, add tolerance in the fit process,
    add the plot of the comparison between model and DFT.
    By A. Martin, M. Verstraete and Ph. Ghosez.

D.6 Adjustment of tutorial tutoparal ucrpa, see test tutoparal#tucrpa_4.
    By B. Amadon

D.7 The ddk file is now available in netCDF format (lightweight version without first-order wavefunctions), 
    and test with the optic post-processor has been set up.
    See the new [[test:v7_49]]. 
    By M. Giantomassi

D.8 Continued development of the electron-phonon [[optdriver]] = 7 module of ABINIT.
    New input variable [[tmesh]], defining a linear mesh of temperatures, see tests [[test:v8_44]] and [[test:v8_45]].
    Also, debugging and improvement of doc.
    By M. Giantomassi

D.9 Added netcdf output of phonons for full grid, not just band structure. Only in tetrahedron prtdos 2 case.
    By M. Verstraete

D.10 On-going development: main executable `tdep`, for the TDEP algorithm, by Hellman and coworkers.
     See [[src:98_main/tdep.F90]], as well as directory 80_tdep. 
     No automatic tests provided yet, no documentation as well ...
     By F. Bottin, J. Bouchet, J. Bieder.

D.11 Capability to print perturbed vxc potential in response function calculations.
     By G. Antonius

D.12 On-going modularization of all source F90 files, to get rid off abilint.
     By M. Giantomassi

D.13 On-going improvements in the doc, to benefit from the new processing capabilities,
     like central bibliography, matjax, etc ...
     By M. Giantomassi, X. Gonze

D.14 Post-processing script for Raman calculations (script/post-processing/Raman_spec.py).
     Reads the anaddb output file and extracts the Raman tensor and then calculates 
     the Raman spectra as a function of frequency at a user-defined temperature.  
     Additionally, the script will automatically extract the dielectric tensor as a function of frequency if it is available.
     By N. Pike

D.15 Refactoring for DFPT+NON_COLL: first version for PAW
     By M. Torrent

D.16 Fix of several bugs in constrained magnetization calculations [[magconon]].
     By E. Bousquet

D.17 Wrong sign of the derivative of spherical harmonics for f orbitals.
     Can lead to problems in BSE and GW calculations if pseudos with explicit f-projectors are used.
     Detected and corrected by Henrique Pereira Miranda.
     See [[gitsha:87617d261b5905e368081af4b899b3ddd7ec83fe]] for the corrected expressions.

D.18 Add possibility to do DFT+U calculations without spin polarization 
     in the exchange and correlation functional: the spin polarization thus only comes from the U and J terms. 
     Can be used with [[usepawu]] = 4, but still under tests.
     By B. Amadon

D.19 GW is now available with [[nspinor]] = 2 with or without spin-orbit coupling,
    and with [[nspden]] = 1 or 4 (collinear or non-collinear spin-magnetisation).
    Implemented only in the norm-conserving case. Still under testing.
    See tests from [[test:v8_90]] to [[test:v8_93]]. 
    By M. Giantomassi.

D.20 Miscellaneous additional bug fixes and improvements of documentation by:
     L. Baguet, W. Chen, C. Espejo, M. Giantomassi, Y. Pouillon, M. Torrent, J. Zwanziger.


* * *
 
## v8.6 

Many thanks to the contributors to the ABINIT project between
May 2017 and October 2017. These release notes
are relative to modifications/improvements of ABINITv8.6 with respect to v8.4.

The list of contributors includes:
B. Amadon, G. Antonius, L. Baguet, J.-M. Beuken, J. Bieder, F. Bottin, Y. Bouchet, E. Bousquet,
M. Giantomassi, O. Gingras, Ph. Ghosez, M. Giantomassi, X. Gonze, F. Jollet, J. Junquera, A. Martin,
F. Naccarato, G. Petretto, N. Pike, Y. Pouillon, S. Prokhorenko, M. Torrent, M. Verstraete, J. Wiktor, J. Zwanziger

It is worth to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases.
This might take some time ...

Xavier

* * *

Version 8.6, released on November 3, 2017.

List of changes with respect to version 8.4 .

* * *

###A. Warnings and important remarks

A.1 The interface between ABINIT and TRIQS has been changed, such that TRIQS 1.4 is used now instead of TRIQS 1.3.
    No backward compatibility with TRIQS 1.3 has been coded, sorry. See D.4.

A.2 Some changes of names:

* input variable gwls_sternheimer_kmax has become [[gwls_stern_kmax]];
* input variable gwls_dielectric_model has become [[gwls_diel_model]];
* input variable prt_effpot has become prt_model (multibinit input variable);
* input variable effmass has become [[effmass_free]];
* tutorial tlda has become tdftu;

    Also, the input variable gwls_second_model_parameter, not used, has been suppressed.

A.3 The definition of Hund's coupling J as computed within cRPA has changed: it now uses the same convention as
    the convention used (for the variable [[jpawu]]) in DFT+U and DFT+DMFT in ABINIT (automatic tests, and tutorial
    are accordingly changed).
    By B. Amadon

* * *

###B.  Most noticeable achievements

B.1 The whole ABINIT documentation has been significantly improved by the introduction of Topics, replacing the
    previous "lists of ABINIT features". All the capabilities of ABINIT are now presented in about 70 short topic Web pages.
    Those topic web pages usually have:
    - a brief introduction;
    - the list of related tutorials -if any-;
    - the list of related input variables (ordered according to their importance for the topics -compulsory, basic, useful or expert-);
    - possibly example input files;
    - list of references.
    Entry point: see the new header of any ABINIT documentation file (e.g. the [new user's guide](..) )
    By F. Jollet and X. Gonze (also tests/fixes by B. Amadon, M. Torrent).

B.2 A central [[theory:bibliography]] database abiref.bib has been created, and linked to the
    above-mentioned topics (B.1) but also to other parts of the ABINIT documentation (e.g. input variable list,
    the tutorials, the theory documents, the acknowledgments ...).
    More than 200 bibliographical references are present.  Not all bibliographical references of the doc have been entered in this
    central bibliographic database, though.
    By X. Gonze and F. Jollet.

B.3 The list of input variables has been restructured, and is now searchable.
    The input variables for anaddb, aim and optic have been included in the database.
    By J. Bieder, X. Gonze and F. Jollet.

B.4 The frequency of usage of each input variable (in the set of automatic tests) is now automatically
    computed, and mentioned in the documentation for this input variable. Examples input files are also now mentioned in the documentation.
    The input files for the automatic tests can now be directly accessed on the Web, as well as the reference files for the tutorials.
    By. X. Gonze and F. Jollet.

B.5 Several important developments related to electron-phonon matrix element computation have been made.
    The Frohlich interpolation procedure for electron-phonon matrix elements, as explained in PRL 115, 176401 (2015), has been implemented.
    The long-range part of the phonon coupling potential is modeled with the Born effective charges and the dielectric tensor.
    This long-range part is substracted from the potential before the Fourier interpolation then added after the interpolation.
    The resulting potential is in much better agreement with the full calculation, as can be verified from the el-ph matrix elements.
    Also, a functionality has been added in the eph driver ([[eph_task]]=5) to only interpolate the phonon potential onto a fine q-point grid.
    The interpolation is performed one perturbation at a time, and is thus more memory efficient than the previous procedures.
    There is additional testing of the new "driver" optdrive=7 specifically dealing with electron-phonon
    related computations (including zero-point renormalisation), especially the interpolation.
    The symmetries have been fixed.
    See new tests [[test:v8_61]]-[[test:v8_65]].
    By  G. Antonius and M. Giantomassi.

B.6 ABINIT can now read pseudopotentials in the PSML 1.1 format, as described in https://arxiv.org/abs/1707.08938. This XML-based format is
    produced by ONCVPSP 3.2 and 3.3, as well as SIESTA's ATOM 4.2, and allows to perform calculations with the exact same pseudopotential files
    in both ABINIT and SIESTA. See the new directory ~abinit/tests/psml, tests [[test:psml_01]] to [[test:psml_14]].
    Note: patches are provided at https://launchpad.net/pspgenpatch to enable PSML output in ONCVPSP.
    By Y. Pouillon, M. Verstraete, J. Junquera and A. Garcia.

B.7 ABINIT is now interfaced with Libxc 3.0. The interface with Libxc 4.0 is in preparation.
    Tests [[test:libxc_06]], [[test:libxc_07]], [[test:libxc_17]], [[test:libxc_18]], [[test:libxc_20]], [[test:libxc_21]] 
    have been modified, because some functionals of libxc v2.0 have changed category in v3.0.
    By M. Torrent.

B.8 A new tutorial, called [[tutorial:positron|Electron-positron annihilation]] has been created.
    By J. Wiktor and M. Torrent.

B.9 The new input variable [[chkdilatmx]] has been introduced, to allow expert users to make
    ABINIT bypass the stopping criterion related to dilatmx. In practice, if the condition related
    to dilatmx is not met, ABINIT continues, and delivers an (approximate) optimized geometry and energy,
    that might be used by external drivers like e.g. USPEX to continue the search for global optimized structures.
    See input variable [[chkdilatmx]], and test [[test:v3_42]].
    By X. Gonze.

B.10 Implementation of the LDA-1/2 methodology.
     Tests to be provided.
     By F. Jollet.

* * *

###C. Changes for the developers (also compilers)

C.1 There are large changes of the procedure to document ABINIT, for most of the documentation files.
    The HTML files are now produced from YAML files, under the control of the script ~abinit/doc/generate_doc.py .
    The documentation that describes this procedure is available on the ABINIT wiki, at https://wiki.abinit.org/doku.php?id=developers:generate_doc .
    This is directly linked to the modifications in the doc presented in B1-B4.
    In particular, the Dokuwiki syntax is used for the hyperlinks.
    By F. Jollet and X. Gonze.

* * *

###D.  Other changes (or on-going developments, not yet finalized)

D.1 The "Adaptively Compressed Operator" approach for the fast application of the Fock operator has been implemented,
    and replaces the traditional way to apply the Fock operator in hybrid functionals (e.g. HSE06, PBE0, ...).
    Hybrid functionals are not yet to be considered in production, though (see D.6), but likely for ABINITv8.8.
    See tests libxc 51, 52, 53, 67, 68, 69, 70, 71, and also v7#65, 66, 70.
    By F. Jollet and X. Gonze.

D.2 A set of 6 input files and accompanying references (from 32 to 2048 procs), for benchmarking high-performance computing
    is available in the new directory ~abinit/tests/hpc .
    Not yet tested automatically, but this future capability is prepared.
    By M. Torrent.

D.3 The tutorial on the temperature-dependent electronic structure has been imported from the ABINIT wiki to the
    usual location ~abinit/doc/tutorial and suppressed from the Wiki. However, it is not yet operational.
    Work is also going on on tutorial fold2bloch.
    By X. Gonze.

D.4 Interfacing with TRIQS 1.4 (instead of 1.3).
    By O. Gingras, B. Amadon and J.-M. Beuken.

D.5 Anaddb can now interpolate and print out the DDB onto an arbitrary set of q-point.
    The same procedure was used to produce the phonon band structure.
    Now, with the input variable prtddb, anaddb will produce both the _DDB file and the _DDB.nc files, the latter being separated for each q-point.
    By G. Antonius.

D.6 On-going work on hybrid functionals: speed-up of the SCF loop, computation of stresses,
    joint computation of forces and stresses, downsampling the wavevectors.
    By X. Gonze and F. Jollet.

D.7 On-going work on the implementation of the TDEP algorithm (temperature dependent sampling).
    By J. Bieder, F. Bottin and Y. Bouchet.

D.8 Replacements of http:// by https:// in many documentation files.
    By J.M. Beuken.

D.9 Test of non-magnetic LDA+U and LDA+U+SO.
    See the new test v5#16
    By M. Torrent.

D.10 Make LDA+U and local EX-exchange compatible with nspden=1/nspinor=2
     By M. Torrent.

D.11 Test of the Velocity Verlet algorithm ionmov=24
     See the new test v8#13
     By S. Prokhorenko.

D.12 Make thermally occupied supercell of a given size, with input variable thermal_supercell.
     See test v8#46
     By M. Giantomassi

D.13 Write dielectric tensor to anaddb.nc when only perturbations w.r.t. electric field are present; Test for nlflag=2,3
     Test the computation of the nonlinear coefficients and first change of dielectric tensor.
     See test v8#47-50
     By F. Naccarato.

D.14 Ongoing work: Raman intensities, in the PAW case, using DFPT.
     By L. Baguet and M. Torrent.

D.15 Ongoing work on the multibinit project.
     New hist storage.
     New tests paral#101-102, to test the anharmonic part.
     Rationalization of supercell treatment with multibinit
     By A. Martin, M. Verstraete and Ph. Ghosez.

D.16 On-going work on the extension of DFPT within non-collinear magnetism.
     By F. Ricci, S. Prokhorenko, M. Verstraete, M. Torrent and E. Bousquet.

D.17 On-going work on DFPT with magnetic field perturbation (Zeeman field).
     By S. Prokhorenko and E. Bousquet.

D.18 Begin transport epc calculations within the eph part of the code.
     By M. Verstraete.

D.19 Improvements for the reading and initialization of density (esp. nspden=4).
     By M. Torrent.

D.20 New interface to the build system and new ac8 config file format.
     By Y. Pouillon.

D.21 Store fold2bloch results in NetCDF format.
     Add new option to cut3d to convert DEN/POT from Fortran to netcdf
     By M. Giantomassi.

D.22 Add mdtemp in the _HIST file. This can be useful (mandatory) for post processing MD/PIMD
     Add imgmov in the _HIST file. Convenient to know if it is PIMD or NEB/string for postprocessing
     By J. Bieder.

D.23 Use inversion symmetry if nspden == 4 and NC
     By M. Giantomassi.

D.24 Update elastic tutorial
     By J. Zwanziger.

D.25 Add LO-TO terms to netcdf files
     By M. Giantomassi.

D.26 Numerous miscellaneous additional bug fixes and improvements of documentation by:
     G. Antonius, J. Bieder, M. Giantomassi, F. Jollet,
     G. Petretto, N. Pike, Y. Pouillon, M. Verstraete, M. Torrent.

* * *

## v8.4

Many thanks to the contributors to the ABINIT project between
January 2017 and May 2017. 
These release notes are relative to modifications/improvements of ABINITv8.4 with respect to v8.2.

The list of contributors includes:
F. Altvater, G. Antonius, L. Baguet, J.-M. Beuken, J. Bieder, E. Bousquet, 
W. Chen, G. Geneste, M. Giantomassi, Y. Gillet, X. Gonze, F. Jollet, A. Martin, 
F. Naccarato, G. Petretto, S. Prokhorenko, F. Ricci, M. Torrent, M. Verstraete, J. Zwanziger

It is worth to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases ...
This might take some time ...

Xavier

* * *

Version 8.4, released on June 2, 2017.

List of changes with respect to version 8.2 .

* * *

A.  Warnings and important remarks

A.1 The Github ABINIT Web page, at https://github.com/abinit, allows one to
    access the mirror of the repository of ABINIT, as well as some other projects
    related to ABINIT, like AbiPy or the pseudo-dojo.

A.2 The content of the ABINIT Web portal, specifically for the pages related to the presentation of 
    ABINIT (http://www.abinit.org/about/what-is-abinit), has been upgraded.  

* * *

B.  Most noticeable achievements

B.1 Implementation of algorithms to interpolate the electronic band structure,
    based either on "star functions" of on "B-splines" (as alternatives to Wannier function interpolation).
    See the input variables [[einterp]], [[nkpath]], and [[prtebands]], and tests 
    old syntax: `Tlibxc#42, Tv8#04` replaced by 
    [[tests/libxc/Input/t41.in]], [[test:v8_04]]. 
    Work by M. Giantomassi

B.2 The Fock mixing factor for the HSE hybrid functional can be tuned thanks to the input variable gwfockmix  .
    (Warning: this is for the GW-type approach to electronic structure only, not for total energies)
    See test Tlibxc#43 .
    Work by W. Chen.

B.3 Implementation of spatially-varying chemical potential, for each atomic species.
    See the input variables [[chempot]] and [[nzchempot]], and tests Tv8#30 and 31.
    Work by X. Gonze

B.4 Linear geometrical constraints can now be imposed on PIMD runs.
    This allows ABINIT to use the blue moon sampling technique.
    See the input variable pimd_contraint, and test Tv8#05
    The use of history (HIST) files is now possible with images.
    Work by G. Geneste and M. Torrent.

B.5 Computation of linear reponse (optics executable as well as BSE part of ABINIT)
    in the case of temperature-dependent electronic structure.
    See tests Tv67mbpt#50-53.
    WARNING: This capability of ABINIT has not been fully tested. However, the
    basic tests for non-spin-polarized simple semiconductors are OK.
    As usual, use at your own risk.
    Work by Y. Gillet and M. Giantomassi.

B.6 Computation of Gruneisen parameters by finite differences, within ANADDB.
    See test v8#45
    Work by M. Giantomassi

B.7 New LOBPCG implementation [[wfoptalg]] = 114.
    This is the default when [[paral_kgb]] == 1.
    Performances are equivalent in standard use (MPI alone) and much better with openmp+multithreaded linalg. 
    It also allows to have only one block for large system and to reduce memory copies.
    This version has been developed keeping in mind the next generation of HPC.
    Work by J. Bieder

B.8 New algorithms for the displacement of nuclei ([[ionmov]]):
    - Hybrid Monte Carlo (HMC) predictor (ionmov=25)
    - Velocity Verlet (VV) NVE molecular dynamics predictor (ionmov=24)
    See Tv8#12 for ionmov=24. 
    Work by S. Prokhorenko

B.9 Refactoring of ANADDB for the production of DOS and other thermodynamic quantities,
    also the mean square displacement and mean square velocity.
    The DOS is obtained using usual DOS methods, instead of the histogram method,
    and converges much faster with [[anaddb:ng2qpt]]. Activated with [[anaddb:prtdos]] 1 or 2 in anaddb.
    Work by M. Verstraete.

* * *

C. Changes for the developers (also compilers)

C.1 Management of the test farm: the new bot ubu_intel_17_openmpi 
    has been activated, so that Intel 17 is now supported.
    Also, replacement of shiva_gnu_6.3_py3k by inca_gnu_6.3_py3k,
    update of graphene (MacPorts) to gcc6.3 + scalapack.
    By J.M. Beuken.

* * *

D.  Other changes (or on-going developments, not yet finalized).

D.1 The printing of potentials (e.g. [[prtvxc]] ...) works now for DFPT.
    By M. Verstraete.

D.2 New input variable [[prtphbands]].
    See tests v7#88 and v8#46.
    By M. Giantomassi

D.3 In case the ANADDB interpolation of phonon frequencies produces
    erroneously a negative slope around gamma, the new [[anaddb:nsphere]] = -1 possibility
    allows ANADDB to select a set of IFCs that does not lead to such non-physical behaviour.
    See test v8#46.
    Work by M. Giantomassi.

D.4 Added new value for [[anaddb:nlflag]] = 3 that computes only the non-linear susceptibility.
    By F. Naccarato.

D.5 Computation of forces is now possible in the Wavelets + PAW case.
    By M. Torrent.

D.6 Ongoing work concerning the new "driver" [[optdriver]] = 7 specifically dealing with electron-phonon 
    related computations (including zero-point renormalisation). 
    See new tests v8#41 to 44.
    Work by M. Giantomassi.

D.7 Ongoing work: Raman intensities, in the PAW case, using DFPT. 
    By L. Baguet and M. Torrent

D.8 Ongoing work related to non-collinear DFPT.
    Definition of the new input variable [[rfmagn]], as well as Tv8#20 .
    DFPT with static B-field (Zeeman term only) works for q=(0 0 0).
    Adjustmanet of dfpt corrections for metals to non-collinear case.
    Bug fix for GS calculations with finite magnetic field (case of collinear spins).
    By S. Prokhorenko, F. Ricci, and E. Bousquet

D.9 Ongoing work on the multibinit project.
    NPT simulations are possible.
    New tests v8#06 and paral#100, to check the correct IO of XML files.
    Note that tests Tv7#120-124 have been renumbered Tv8#07-11.
    Work by A. Martin

D.10 New test paral#62, to make sure [[paral_kgb]] = 0 works when there are some idle procs.
     By M. Giantomassi

D.11 Ongoing work on forces and stresses for hybrid functionals.
     By F. Jollet 

D.12 Ongoing work on the electron-phonon postprocessor ElectronPhononCoupling.
     Fix bug affecting MPI runs. 
     By G. Antonius.

D.13 Improvement of zheevd using MAGMA (from a msg on the forum)
     By M. Torrent

D.14 Implemented KSS.nc output with netcdf primitives
     By M. Giantomassi

D.15 Concerning  the Fourier interpolation of the phonon band structure,
     inside ANADDB, work by G Petretto:
     - Updates in the calculation of sound velocity, 
     - Small modifications for [[anaddb:nlflag]] == 3 and added some quantities to the anaddb netcdf file,
     - Fix a bug in the implementation of the new weights

D.16 Concerning Path Integral Molecular Dynamics with Quantum Thermal Bath:
     allow restart from history file.
     By M. Torrent

D.17 Refactored the computation of the electric dipole moment.
     By M. Torrent

D.18 The energy width for the bands in the Boltztrap intrans file is now automatically set.
     Previously a constant 0.4 Ryd, and the user should have checked by hand if
     the value was sufficient. Should be considered a bug fix.
     By M Verstraete

D.19 Numerous miscellaneous additional bug fixes and improvements of documentation by:
     F. Altvater, G. Antonius, J. Bieder, M. Giantomassi, F. Jollet, 
     G. Petretto, M. Verstraete, M. Torrent, J. Zwanziger. 

* * *

## v8.2
    
Many thanks to the contributors to the ABINIT project between
June 2016 and January 2017. These release notes
are relative to modifications/improvements of ABINITv8.2 with respect to v8.0.

Moreover, most of them are also described in the Computer Physics Communications 2016 ABINIT paper, 
doi:10.1016/j.cpc.2016.04.003
    
The list of contributors includes:
B. Amadon, G. Antonius, L. Baguet, J.-M. Beuken, J. Bieder, E. Bousquet, F. Bruneval,
W. Chen, M. Giantomassi, Y. Gillet, X. Gonze, G. Petretto, F. Jollet, A. Martin,
V. Planes, Y. Pouillon, T. Rangel, F. Ricci, M. Torrent, M. Verstraete

It is worth to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases ...
This might take some time ...

Xavier

* * *

Version 8.2, released on February 16, 2017.

List of changes with respect to version 8.0 .

* * *

A.  WARNINGS AND IMPORTANT REMARKS

A.0 The 2016 article by the ABINIT group is now mentioned in the acknowledgments:
    "Recent developments in the ABINIT software package. 
    Computer. Phys. Communications 205, 106 (2016)".
    See http://www.abinit.org/doc/helpfiles/for-v8.2/users/acknowledgments.html, as well as
    the notice at the end of ABINIT runs.
   
A.1 [[inclvkb]] 1 has been removed. Now the possible values are either 0 or 2

A.2 The default strategy of [[so_psp]] has been changed 
    (see the description of the input variable so_psp).

* * *

B.  Most noticeable achievements

B.1 Implementation of the Limited-memory Broyden-Fletcher-Goldfarb-Shanno (LBFGS) 
    minimization algorithm.  Activate this algorithm using [[ionmov]] = 22. 
    From the tests that have been run, this algorithm can be much better
    than the native implementation of BFGS in ABINIT when one approaches convergence, 
    perhaps because of better treatment of numerical details. 
    This algorithm might become the default in ABINIT, if the better behaviour is confirmed.
    Test case: v8#02 .
    The working routines were based on the original implementation of J. Nocedal
    available on netlib.org.  They have been reshaped and translated into modern fortran, 
    then interfaced to ABINIT by F. Bruneval (sources in 45_geomoptim/m_lbfgs.F90).

B.2 A new tutorial is available, on the calculation of the effective interactions U and J 
    using constrained Random Phase Approximation (cRPA) for DFT+DMFT (or DFT+U) calculations.
    See doc/tutorial/_ucalc_crpa.md as well as the automatic tests tutorial/tucrpa#1-5 .
    This tutorial was prepared by B. Amadon.

B.3 Implementation of temperature-dependent spectral functions 
    (electronic spectral function, with electron-phonon interactions),
    as well as real part of phonon self-energy (Pi) for gapped systems.
    Also, automatic test for spectral function, v7#89, and improved documentation.
    Work by G. Antonius.

B.4 The RPA one-shot bootstrap fxc kernel has been implemented for GW calculations (gwgamma=-8).
    See Rigamonti et al PRL 114, 146402 (2014) and Berger PRL 115, 137402 (2015).
    The test v67mbpt#36 has been updated.
    Work by W. Chen.

* * *

C. Changes for the developers (also compilers)

C.1 The version control system that is used for the development of ABINIT has been changed: 
    the whole ABINIT project
    has been ported from bzr to git.
    Work by J.-M. Beuken, Y. Pouillon, M. Giantomassi, X. Gonze, 
    with discussions with many developers.

C.2 New versions of Fortran compilers have been integrated in the test farm:
    - intel 16.0
    - gnu 6.1 and 6.2
    - IBM xlf compiler 14.1
    - NAG 5.3
    Corresponding examples are available in doc/build/config-examples.
    On the contrary, g95 is not tested anymore.
    Work by J.-M. Beuken

C.3 The v8 directory for tests has been initialized.
    By Matteo Giantomassi.

C.4 Python 3 >= 3.4 is now supported in build system scripts 
    (compatibility with py2_ version >= is maintained).
    By Matteo Giantomassi.

* * *

D.  Other changes
(or on-going developments, not yet finalized).

D.1 The main executable "multibinit" has been created.
    Its goal is to perform "second-principles calculations", building model Hamiltonians
    using the data provided by the DDB (or other info from ABINIT).
    Tests v7#120-124 (should be renamed) as well as paral#95-98.
    Work by A. Martin.

D.2 A new "driver" within ABINIT has been defined, specifically dealing with electron-phonon 
    related computations (including zero-point renormalisation). 
    Set optdriver=7 . 
    New input variables: ddb_shiftq, eph_task, eph_transport, prtphdos, prtphsurf.
    See tests v7#88 and 89.
    Work by M. Giantomassi and G. Antonius.

D.3 The generation of k-point meshes with kptrlatt and shiftk is now tested.
    See test v8#03 .
    Work by M. Giantomassi

D.4 As a follow-up of the Achievement B3 in the release notes of ABINITv8.0 (DMFT + TRIQS),
    new input variables have been defined for DMFT: dmft_tolfreq and dmftctqmc_triqs_nleg.
    Automatic tests have been set-up, tests v8#01 and paral#99.
    Work by B. Amadon and V. Planes

D.5 More systematic tests of the IO in parallel (different files) have been set up,
    in the norm-conserving, PAW and PAW + spin-orbit cases.
    See tests mpiio#26, 27 and 28.
    Also, the case [[paral_kgb]] = 0 is now tested with idle processors. See test paral#62.
    Work by M. Giantomassi

D.6 The load balancing for the repartition of the plane waves among procs is now monitored,
    and a dynamical equilibration is made possible thanks to the new input variable pw_unbal_thresh.
    See test mpiio#26.
    Work by M. Torrent

D.7 The capability to output spin-resolved DOS in case of spin-orbit NC pseudopotential
    calculations has been checked, and a test has been set-up (test v7#17).
    Work by M. Giantomassi

D.8 Files generated by ABINIT, and used by BOLTZTRAP are now tested.
    See v6#11.
    Work by M. Giantomassi

D.9 Unit tests (fftprof) have been set up for the use of the MKL-DFTI routines: 
    unitary#tfftmkl_03 and 04.
    Work by M. Giantomassi

D.10 Ongoing work: Raman intensities, in the PAW case, using DFPT. 
     By L. Baguet and M. Torrent

D.11 prtvcbm working with all parallelizations.
     Tests mpiio 26:28 have been reactivated on 4 processors.
     By M. Torrent

D.12 On going work related to non-collinear DFPT.
     By F. Ricci, S. Prokhorenko, and E. Bousquet

D.13 MBPT: support for the commutator [Vnl r] in the case of NC pseudos 
     with more than one projector per l-channel has been added. 
     Tests in v67mbpt[40] (GW run with psp8 files).
     SOC is not yet available, though.
     The KSS file continues to use the old implementation to maintain backward compatibility 
     hence outkss won't produce the file if multiple-projectors are detected
     inclvkb 1 has been removed. Now the possible values are either 0 or 2
     By M Giantomassi

D.14 For Hirshfeld and Bader: doc and warning.
     For Hirshfeld charges the output was unclear: the density integral (electrons only) 
     was called the hirshfeld charge, as opposed to the net one. 
     For Bader there was no check that the core charge file 
     corresponded to the pseudopotential used. 
     Now checks at least that it integrates to znucl-zion, which is a first step. 
     Important as the default fc files on the web do not have semicore electrons, 
     whereas many of the psp8 will.
     Contribution by M. Verstraete.
     
D.15 Updated acknowledgments, including the 2016 paper.
     By X. Gonze

D.16 Ongoing work on forces and stresses for hybrid functionals.
     By F. Jollet 

D.17 Ongoing work concerning weights for the Fourier interpolation inside ANADDB.
     By G. Petretto

D.18 Numerous miscellaneous additional bug fixes 
     (to the sources, as well as to the build system, including patches for the fallbacks), 
     and improvements of documentation by:
     G. Antonius, L. Baguet, J. Bieder, F. Bruneval,
     M. Giantomassi, Y. Gillet, G. Petretto, Y. Pouillon,
     M. Verstraete, M. Torrent (in particular, for DFPT+PAW).

* * *
