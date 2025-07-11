## v10.4

Version 10.4, released on May 15, 2025.
<!-- Release notes updated on June 4, 2025. -->

These release notes
are relative to modifications/improvements of ABINIT v10.4 with respect to v10.2.
<!-- Initially, beta-release v10.4.1, merge requests from MR1063, 1068, 1077, 1079, 1081, then 1083 to 1175, 
EXCEPT MR1087, 1088, 1090, 1091, 1092, 1095, 1097, 1100, 1104, 1106, 1114, 1115, that had alreacy been included in v10.2,
and EXCEPT MR1160, 1165, 1168, 1170, 1172 and 1173  that will be included in v10.6 .  
For release v10.3, also MR1176 to 1182 are included, as well as MR1184-1186, MR1188 and MR1190. MR1183 and MR1187 will be included in v10.6.
-->
<!-- (with some late contributions until XXX 2024). -->

Many thanks to the contributors to the ABINIT project between October 2024 and May 2025 !
The list of contributors includes:

G. Allemand, B. Amadon, G. Antonius, M. Azizi, L. Baguet, S. Bandyopadhyay, L. Bastogne, 
R. Béjaud, J.-M. Beuken, A. Blanchet, F. Bottin, J. Boust, J. Bouquiaux, F. Brieuc, 
E. Castiel, Q. Delacroix, F. Gendron, G. Geneste, M. Giantomassi, 
O. Gingras, F. Gomez-Ortiz, X. Gonze, R. Gupta, V. Karasiev,
I.-M. Lygatsika, O. Mattelaer, M. Mignolet, G. Petretto,
S. Ponce, S. Rostami, M. Royo, M. Sarraute, A. Sasani, Le Shu, M. Torrent, V. Trinquet, Hsiaoyi Tsai, 
V. Vasilchenko, M. Verstraete, Xu He, Huazhang Zhang, J. Zwanziger.

It is worthwhile to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases.
This might take some time ...

Xavier

* * *

### **A.** Remarks and warnings.

**A.1** Remark: forthcoming ABINIT paper

ABINIT will be the subject of a forthcoming paper, that will contain the brief description of the new features of ABINIT since 
the 2020 previous papers [[cite:Gonze2020]] and [[cite:Romero2020]]. This will provide complementary information to the present release notes.
So, when it will be available, the user is strongly suggested to rely on both the present release notes and this paper.

**A.2** Downgrading of [[ionmov]], in favor of [[geoopt]] and [[moldyn]] 

The variable [[ionmov]] (used to define geometry optimization and molecular dynamics) is now considered to be a developer variable. 
Two new input variables should be preferred to drive these calculatioins: [[geoopt]] and [[moldyn]]. 
For example, [[ionmov]] 2 corresponds now to [[geoopt]] "bfgs" .

All documentation and tutorials include these changes (there are many input files now using these). 
For backward compatibility reasons, [[ionmov]] can still be used. 
Multibinit has not been modified. 

See [[topic:MolecularDynamics]] and [[topic:GeoOpt]] for more informations.

By Augustin Blanchet, Romuald Béjaud, Fabien Brieuc, Grégory Geneste (MR1085). See also point B.5 .


**A.3** Renaming of test suites

* The directory tests/gwr has been renamed tests/gwr_suite
* The directory tests/gwpt has been renamed tests/gwpt_suite

This is to avoid confusion with the tutorial bearing the same name.
Accordingly, the new suite of tests for real-time time-dependent DFT has been named tests/rttddft_suite .

By M. Giantomassi and X. Gonze.

* * *


### **B.** Most noticeable developments

**B.1** Variational Polaron Equations

The variational polaron equation implementation in ABINIT is ready for production.
It allows one to compute polaron formation energies and associated quantities, in the
strong-coupling approximation. See [[cite:Vasilchenko2022]], and a forthcoming Vasilchenko2025 paper (reference to be added when published).
Activate it using [[optdriver]]=7 and [[eph_task]]=13 (see also [[eph_task]]=-13).
The corresponding [[tutorial:eph4vpq]] is available.
In this tutorial, there are 10 tests, [[test:eph4vpq_1]] to [[test:eph4vpq_10]].

Related input variables: [[getvpq]], [[getvpq_filepath]],[[vpq_aseed]], [[vpq_avg_g]], [[vpq_gpr_energy]], [[vpq_gpr_length]],
[[vpq_interp]], 
[[vpq_mix_fact]],
[[vpq_mesh_fact]],
[[vpq_nstates]],
[[vpq_nstep]],
[[vpq_nstep_ort]],
[[vpq_select]],
[[vpq_pkind]],
[[vpq_tolgrs]],
[[vpq_translate]],
[[vpq_trvec]].
The previously introduced input variables varpeq_<...> (or containing the string varpeq) have been renamed with the prefix vpq_ (or the string vpq) .
 
By V. Vasilchenko, with help from M. Giantomassi, X. Gonze and S. Ponce (MR1126, 1161, 1164, 1166)


**B.2** Real-time Time-Dependent Density-Functional-Theory

The real-time time-dependent density-functional theory implementation in ABINIT is ready for production.
It includes the calculation of the current density response to an impulse electric field and associated response functions
(conductivity, dielectric function). Different time-propagators are coded, see [[td_propagator]].
Activate it using [[optdriver]]=9. Control the number of time steps using [[ntime]], and the time step value using [[dtele]].
The corresponding [[tutorial:rttddft]] is available.

Related input variables : [[prtcurrent]], [[td_propagator]], [[dtele]], [[td_ef_type]], [[td_ef_induced_vecpot]], [[td_ef_tzero]], [[td_ef_pol]]
The new input units "fs" and "as" have been introduced (for "femtoseconds" and "attoseconds").
Related tests : [[test:rttddft_suite_01]] to  [[test:rttddft_suite_06]] as well as the tests of the tutorial, [[test:rttddft_1]] to [[test:rttddft_4]].

By F. Brieuc (MR1137)

 
**B.3** Electron-phonon interaction : thermal transport, cumulant and miscellaneous improvements.

Thanks to the treatment of thermal gradients in addition to electric fields, 
the calculation of Seebeck, electronic thermal conductivity and Peltier coefficients is available both in the Relexation-Time Approximation 
(RTA - [[eph_task]]=7)
and in the Iterative Boltzmann Transport Equation (IBTE - [[eph_task]]=8).
The implementation of chemical potential varying with temperature in the transport coefficient calculations has been fixed (in the RTA).
Related tests for RTA: [[test:v9_60]], [[test:v9_61]], [[test:eph4mob_5]], [[test:eph4mob_6]], [[test:eph4mob_7]].
Related tests for RTA and IBTE: [[test:v9_65]].

The calculation of mobility using Dyson-Migdal spectral functions 
(that can be generated by different means, including the cumulant approach - [[eph_task]]=9)
has also been implemented. 
Related [[test:v10_104]] to [[test:v10_107]]

Support for [[useylm]]=1 has been implemented in EPH routines.

There is a new [[eph_task]]=18 to compute e-ph matrix elements along a k- or q-path. See [[test:eph4zpr_9]] 
as well as the new input variables [[eph_path_brange]], [[eph_fix_korq]], and [[eph_fix_wavevec]].

By G. Allemand, R. Gupta, M. Giantomassi, M. Verstraete (MR1156, 1162, 1102)


**B.4** Multibinit 

Multibinit has been improved in many different aspects.

* Variable Maximum N-body Term Generation: Users can now generate terms with a variable maximum number of bodies using the fit_max_nbody@multibinit parameter, providing increased flexibility in model construction.
* Weighted Training Set: To improve fitting accuracy for structures near minima and saddle points, the [[fit_weight_T@multibinit]] input variable has been added. This allows users to apply a Boltzmann-distribution-like weighting to training samples, reducing the influence of high-energy, high-force, and high-stress structures. 
* Enhanced Bounding Strategy: A new bounding strategy ([[bound_option@multibinit]]=1) is implemented, restricting bounding to terms with the same "nbody" count. The previous scheme remains available via [[bound_option@multibinit]]=2. 
* Homogeneous Electric Fields: Support for homogeneous electric fields has been added through the [[efield@multibinit]] and [[efield_lambda@multibinit]] input variables. 
* Inhomogeneous Spatiotemporal Electric Fields: Implementation of inhomogeneous electric fields, varying in time and space, is now available. 
This feature has yielded successful results, as documented in arXiv:2501.17057 and arXiv:2411.16395. 
The following new input variables control this functionality: [[efield_period@multibinit]], 
[[efield_gmean@multibinit]], [[efield_gvel@multibinit]], [[efield_type@multibinit]], [[efield_phase@multibinit]], [[efield_sigma@multibinit]], [[efield_background@multibinit]], and [[nefield@multibinit]]. 

There have been also several bug Fixes:

* Symmetry-Equivalent Term Completeness: Corrected an issue where the number of symmetry-equivalent terms was incomplete. 
Dipole-Dipole Interaction Box Consistency: Resolved an inconsistency between the dipole-dipole interaction box and the simulation box when distances in the DDB exceeded the simulation box dimensions.
* Stress Correction Term: Fixed an error in the stress correction term. 
* Spin Monte Carlo: fix a bug which results in crash in spin Monte Carlo.

The input variables fit_ncoeff_per_cycle and fit_drop_rate have been disabled.

See [[test:v9_210]] to [[test:v9_215]] for the inhomogeneous spatiotemporal electric fields, [[test:v9_216]] for Weighting of training set, [[test:v8_104]] for the maximum number of bodies.

By Xu He, Louis Bastogne, Alireza Sasani, Fernando Gomez-Ortiz, Huazhang Zhang, Philippe Ghosez, with inspiration from Subhadeep Bandyopadhyay (MR1146 gives more details on each contribution, including the author or inspirator)


**B.5** New molecular dynamics algorithms, new features in the _HIST file, PIMD changes.

New molecular dynamics algorithms

* A new Langevin algorithm, [[moldyn]]="nvt_langevin" ([[ionmov]]=16), is now available for molecular dynamics. It was originally used for PIMD simulations ([[imgmov]]=9). "nvt_langevin" directly calls the PIMD algorithm, automatically setting [[nimage]]=1.
* A new molecular dynamics algorithm, [[moldyn]]="npt_langevin", is introduced for use only with [[optcell]]=2, incorporating friction coefficient for the thermostat ([[friction]]) and for the barostat ([[frictionbar]]).
Since both algorithms are non-deterministic, the restart of a molecular dynamics run controlled by [[restartxf]] is done only by reading the _HIST.nc file, bypassing the prediction algorithms.

These algorithms are tested in [[test:v10_20]] and [[test:v10_21]].

New features in the _HIST.nc file

* The output file _HIST.nc now includes the evolution of [[vel_cell]] (the VELocity of the CELL parameters). 
Previously, this variable was only included in _HIST.nc file for PIMD simulations (non-zero [[imgmov]]).
* The Number of Degrees Of Freedom (ndof) is now printed as a non-evolving integer scalar in the _HIST.nc file. 
Previously, it was not possible to determine the relationship between the kinetic energy and ionic temperature only from the _HIST.nc file. 
Post-processing tools may include this change, computing the ionic temperature as: T = 2Ek/ndof. 
The internal variable ndof is initialized as ndof=3N for non-contrained molecular dynamics algorithms. 
It is properly set to 3N-4 for moldyn="isokin" and 3N-3 for some other contrained algorithms. 
It also accounts for fixed atoms. However, it should be used with caution, as it is not tested for all algorithms.

PIMD changes

* The friction coefficient for the PIMD langevin algorithm has been switched from the [[vis]] input variable to the [[friction]] input variable.
* In the NPT Langevin algorithm, atomic velocities were previously expressed in reduced coordinates, leading to incorrect kinetic energy values the _HIST.nc file. 
This issue has now been corrected.

By Augustin Blanchet, Romuald Béjaud, Fabien Brieuc, Grégory Geneste (MR1085) - See also point A.2 .


**B.6** Support for entropy for thermal exchange-correlation functionals.

At finite temperature, the exchange-correlation free energy is defined as:
$F_\mathrm{xc} = E_\mathrm{xc} - TS_\mathrm{xc}$.

ABINIT was initially designed for zero-temperature exchange-correlation (xc) functionals, where $S_\mathrm{xc} \equiv 0$ and $F_\mathrm{xc} \equiv E_\mathrm{xc}$. However, this is not the case for finite-temperature xc functionals.
To obtain the correct total internal energy and total entropy, both needed for accurate equation-of-state data, $S_\mathrm{xc}$ must be accessible when summing energies after the SCF cycle.
Otherwise, only properties derived from the free energy, such as stresses, pressure, and forces, will be correctly computed. This implementation ensures that $S_\mathrm{xc}$ is properly accounted for.

In a previous release, ABINIT was adapted to communicate the electronic temperature to LibXC. However, the issue is that while LibXC provides the exchange-correlation free energy, it does not yet return the entropy contribution.
We have implemented native exchange-correlation free energy functionals to ensure accurate equation-of-state data.
We have introduced two new native thermal xc functionals:

* [[ixc]] = 51: corrKSDT (LDA based) [[cite:Karasiev2014]].
* [[ixc]] = 60: KDT16 (GGA based) [[cite:Karasiev2018]].

These functionals are not compatible with spin-polarized calculations yet. The documentation includes them, as well as a warning while using LibXC thermal functionals.
The already implemented Ichimaru thermal LDA functional ([[ixc]]=50) now sets xc entropy the same way as KSDT and KDT16 functionals.

See [[test:v10_22]] and [[test:v10_23]].

By A. Blanchet, V. V. Karasiev and M. Torrent (MR1121, MR1135).


**B.7** Computation of stress in mGGA

mGGA exchange-correlation functionals for the energy (e.g. R2SCAN01 and the other mGGA functionals available in libxc) are now 
available with the stresses, allowing to use these for cell optimization. 

See [[test:v9_95]].

By M. Torrent (MR1099)


**B.8** Orbital magnetization, chemical shielding, and electric field gradient : mGGA, abipy integration, GPU.

* The mGGA has been implemented for the DDK DFPT perturbation, enabling mGGA calculations of the orbital magnetization, that has been tested with R2SCAN, both with and without spinors.
* ZORA terms coupling electron spin to nuclear magnetic dipole moments are now computed. This accomplishes accurate chemical shielding calculations in the presence of significant relativistic effects.
* Orbital magnetism and electric field gradient have been added to the GSR netcdf file. Together with previous ORBMAG.nc additions, this completes abinit work for integration with abipy developments. 
* GPU porting : DDK and orbital magnetism response with nuclear dipole moment and [[nspinor]]= 2 are available.

See [[test:v10_41]], [[test:v10_42]] and [[test:gpu_omp_31]].

By J. Zwanziger with help from M. Giantomassi (MR1101, 1111, 1119, 1128)



**B.9** GPU porting of DMFT,  mGGA and spin-orbit and miscellaneous GPU improvements.

The global GPU porting of ABINIT using recent libraries/compilers, started three years ago, has been continued.
In the previous releases 10.0 and 10.2, ground-state calculations [[optdriver]]=0 including Fock operator, 
and Density-Functional Perturbation Theory [[optdriver]]=1 had been made available.

In the present release 10.4, DMFT and mGGA have been ported to GPU, also using OpenMP. Also, the previous implementations have been generalized to the spin-orbit case.
Bugs have been fixed in the previous implementations, slowdowns have been eliminated, and some documentation has been added.
The implementation is working with compiler NVHPC 24.9 openmpi (in addition to NVHPC 23.9).

The new input variable [[gpu_thread_limit]] has been introduced, allowing to avoid to use a high number of threads in some OpenMP regions of the code (NOT TESTED).

Numerous new tests are available, [[test:gpu_omp_30]] to [[test:gpu_omp_46]], including checks for PAW and norm-conserving formalism and nspden 4.
Also, [[test:gpu_06]] has been added for GPU+PAW+non-collinear magnetism.

See the description of the GPU possibilities of ABINIT in the documentation, input variable [[gpu_option]]=2 for the OpenMP capabilities.

By M. Sarraute, with help from L. Baguet and M. Torrent (MR 1120, 1123, 1125, 1129, 1134, 1136, 1142, 1144, 1163)



**B.10** aTDEP : Command-Line Interface, file format, and compatibility with Abipy.

The aTDEP utility has been upgraded to match the standard of most executables in the ABINIT package.

* New syntax : one can specify the input file as a command line argument, e.g.  atdep run.abi
* Three new input variables were introduced to replace the files file : [[output_file@atdep]], [[outdata_prefix@atdep]] and [[indata_prefix@atdep]].
* The IO name conventions are now in accordance with other executables inside the ABINIT package.
* The standard output is used instead of forcing the name 'atdep.log'.
* The input file admits free format. One can give input variables in any order, spread them over several lines, etc.
* The atdep utility now outputs a DDB file.
* New input variable [[debug_mode@atdep]].
* There is a final summary YAML block at the end of atdep standard output.

Moreover, aTDEP has been debugged and made more portable (no exclusion for selected workers of the ABINIT test farm). There are more tests 
of aTEDP features, e.g for initial atoms translation to the same unit cell and for rotation of the crystal to the conventional orientation.
The variable [[bzpath@atdep]] is no longer of mixed type (int + str), thanks to [[special_qpt@atdep]].
The optional variable [[znucl@atdep]] has been added so that it is included in the output DDB.

See [[test:atdep_05]], [[test:atdep_05]], [[test:atdep_38]]-[[test:atdep_40]].

By G. Antonius with some contribution by F. Bottin (MR 1105, 1112, 1113, 1117, 1145, 1155)


**B.11** Variable cell NEB 

Generalized solid-state NEB algorithm has been implemented.
New input variable : [[neb_cell_algo]].
See [[test:v10_83]].

By Q. Delacroix and M. Torrent (MR1149)


**B.12** Dynamical quadrupoles for standard pseudopotentials (with non-linear XC core correction).

Previously, dynamical quadrupoles could be computed only using norm-conserving pseudopotentials WITHOUT non-linear XC core correction.
However, the state-of-the-art pseudopotentials include this feature (e.g. pseudo-dojo).
This implementation has been done.

See [[test:v10_43]].

By M. Royo (MR1154).



**B.13** GWR (real-space imaginary-time GW) algorithm : tutorial and miscellaneous improvements.

* The [[tutorial:gwr1]] has been finalized (see also the introductory [[tutorial:gwr_intro]]). 
* [[ecutwfn]] can now be used for GWR. 

For the GWR, use [[optdriver]]=6, and specify [[gwr_task]]. 

See [[test:gwr_1]] to [[test:gwr_5]].

By M. Giantomassi (MR1147, MR1154) 


**B.14** Degree of filtering and its oracle in Chebyshev algorithm

The [[mdeg_filter]] input variable for spectrum filtering algorithms has been defined.
Moreover, it can be automatically computed:
[[chebfi_oracle]], [[oracle_factor]] and [[oracle_min_occ]] can be used to activate and control the 
computation of [[mdeg_filter]] parameter from the residuals. 

See [[test:v10_07]] to [[test:v10_09]] and [[test:paral_34]].

Also, the residuals in chebfi and lobpcg ([[wfoptalg]]==114 and 111) have been corrected, which solves the issue #168 (closed).
The timing in chebfi has been corrected.

By L. Baguet (MR1103) and M. Torrent (MR1098)


**B.15** Abipy : Automated computation of luminescence spectra

A set of Python modules integrated into the AbiPy framework now facilitates the automated computation of phonon-resolved luminescence spectra for defects in inorganic solids.  
The workflow automates key steps in a classical computational process, from the initial ∆SCF DFT calculations with constrained occupations to simulate the excited state, 
to the generation of defect phonon modes in large supercells, and finally the calculation of luminescence spectra based on Huang-Rhys theory, using the generating function formalism.

The workflow is available in the GitHub abipy repo. 
See for example <https://abinit.github.io/abipy/flow_gallery/run_lumi_Eu_doped_SLA.html>, 
and <https://abinit.github.io/abipy/flow_gallery/run_lumi_NV_center.html>.

By J. Bouquiaux, with help from M. Giantomassi (see the forthcoming ABINIT2025 paper and a forthcoming JOSS paper).


**B.16** Abipy : Automated computation of thermal expansion (volumetric as well as anisotropic)

Abipy workflows for thermal expansion calculations are available, within the Quasi-Harmonic Approximation (QHA) for materials with diverse crystallographic symmetries,
with efficient Zero Static Internal Stress Approximation (ZSISA) and its volume-constrained variant, v-ZSISA. 

The v-ZSISA workflow is available in the Github abipy repo,
<https://abinit.github.io/abipy/flow_gallery/run_qha_vzsisa.html> . The ZSISA workflow will soon be available.

By S. Rostami, with help from M. Giantomassi and X. Gonze (see forthcoming ABINIT2025 paper).


* * *


### **C.** Changes for the developers or for the installation of ABINIT (including information about compilers)


**C.1** Remote servers for the test farm.

The possibility to include remote servers in the ABINIT test farm has been developed.
This needed some modifications to the testbot.py script to handle a few more options.

By G. Petretto, with tests by JM Beuken (MR1169)

**C.2** Resource-aware version of runtests.py 

Now, runtests.py and testbot.py are aware of the maximum number of available CPUs/GPUs.
The Python processes share two counters tracking the number of CPUs/GPUs in use and
wait for sufficient resources to become available before launching the test.

By M. Giantomassi and M. Sarraute (MR1116)


**C.2** Improve autotools & CMake build systems 

* Correct a logic error for cmake variable ABINIT_ENABLE_MPI_INTERFACES_BUGFIX
* GPU markers, TRIQS, levmar
* Completion of CMakeList.txt
* Fix buffer overflow in m_vtorho
* Miscellaneous fixes

By M. Mignolet (MR1081), M. Torrent (MR1089, 1151), M. Sarraute (MR1096), M. Verstraete (MR1110)

**C.3** Add documentation on how to debug and profile MPI jobs

In the developer HowTo, added:

* New guide on how to use NVTX markers on CPU (devtool recently enabled since !1151 (merged)) and profile MPI jobs using NVIDIA Nsight.
* New tip on the gdb debugger section on how to use gdb for MPI jobs.

By I. Lygatsika (MR1158)

**C.4** Change HDF5 detection when using "with_hdf5=path/to/hdf5". Same for netcdf and netcdf-fortran

The situation was the following :

When giving a path in the ac9 file ("with_hdf5=path/to/hdf5"), the build system kind of ignores the command line if some environment variables are found (HDF5_CPPFLAGS,HDF5_CFLAGS,HDF5_FFLAGS,HDF5_FCFLAGS,HDF5_LDFLAGS,HDF5_LIBS). In that case the build system uses only the environment variables, whatever the path given in the ac9 file. This feature leads to confusing situations, especially when environment variables are set by the system, not the developer himself.

The new policy is the following :

The build system follows the command line "with_hdf5=path/to/hdf5", whatever the context. If some environment variables are found, a warning is printed and the developer is advised to set "with_hdf5=yes" if he or she wants to use environment variables.

The same change is done for netcdf and netcdf-fortran.
In addition, for netcdf (netcdf-fortran), fcflags and libs are filled with the information given by nc-config (nf-config), which is a binary found in "path/to/netcdf/bin/".
This merge improves hdf5 and netcdf detection on our supercomputers, and works on the testfarm.

By L. Baguet (MR1079).

**C.5** Add netcdf_fortran as dependency for 42_parser and 68_dmft   

The fortran file in 42_parser and 68_dmft both have "use netcdf"
but netcdf_fortran is not marked as dependency which does not seem consistent. This has been fixed.

By O. Mattelaer (MR1140)

**C.6** Intel debug

Corrections to improve the use of debug options from the intel compiler (still 2020 version) :

At first there were crashes at compiling time : the compiler does not like operations on arrays inside a matmul call.
Then many of warnings of temporary arrays have been removed, sometimes allocating explicitly an array (when the subarray is not contiguous, or very small), sometimes using a pointer (when it is contiguous). 
There are still many warnings, but for GS computations with lobpcg or chebfi and 
PAW pseudos the amount of warnings is limited, so one can run big computations without generating gigabytes of output every minute.

The debug options used are : -traceback -check all -check uninit -ftrapuv -debug all -warn declarations
All these changes do not affect the performances at all.

By L. Baguet (MR1139)

**C.7** OneAPI
Adding mpiifx and mpiicx before mpiifort and mpiicc in order to better support the new versions of oneAPI.

By Hsiaoyi Tsai (MR1114)

* * *

### **D.**  Other developments (possibly not yet finalized), other new tests, new input variables, new tutorials.


**D.1** Progresses in Dynamical Mean Field Theory (DMFT) 

There have been many improvements in the DMFT part of ABINIT, not yet though to be considered in production yet, because of lack of tests and/or documentation. This will be for ABINITv10.6 ... 

* There has been a massive optimization of the code, as well as a cleaning.
* New functionalities have been added (new double counting, possibility to work with an arbitrary orbital read from file...). List to be provided ...
* One can run dmft susceptibility with correlated and non correlated atoms.
* A new option to run DFT+DMFT has been added : one can use Wannier90 to construct localized orbitals and call the impurity solver using the python invocation scheme (still to be done: going back to DFT from the impurity solver).
* Several routines have been added in preparation for a (working) interface with TRIQS. About thirty new input variables have been defined, with the prefix dmft_triqs_XXXX (to be tested and documented).
* Two extra components for the energy are printed.
* New input variables (to be tested and documented) : dmft_fermi_step, dmft_nominal, dmft_orbital, dmft_prt_maxent, dmft_prtwan, dmft_shiftself, dmft_test, dmft_use_all_bands, dmft_use_full_chipsi, dmft_wanrad, dmft_x2my2d.

Input variables [[dmftctqmc_mov]], [[dmftctqmc_order]] and [[dmftctqmc_triqs_nleg]] are documented, but not tested.
There has also been miscellaneous DMFT fixes and improvements

By F. Castiel, F. Gendron, O. Gingras, B. Amadon (MR1063, 1077, 1084, 1107, 1148 give more details)


**D.2** Improvements and bug fixes in the implementation of [[cprj_in_memory]]

* Fixed a  memory issue in xg_tools for systems with high number of plane waves (solve the problem presented in merge !1058 (merged)).
* Some MPI communications outside the diago were too long for systems with [[nsppol]]=2, now it's fine.

As a consequence of the first point, xg_tools now uses OpenMP "collapse" directive when appropriate (as discussed with @torrent).

Other minor corrections :

* Prediction of memory footprint for chebfi with gemm_nonlop is now activated only if gemm_nonlop_use_gemm is activated
* Timing of mkinvovl in vtorho was producing negative values of "other" timing in vtorho, so it is suppressed
* Correction to reduce the numerical noise on stress and forces when comparing PAW calculations done with different compilers and/or different machines.
* Change the design in xg_nonlop, which solves several issues :
no more warnings "temporary array was created" for some cases (seen on ubuntu and the testfarm), and
no more crash for other (rare) cases.
* Also set [[cprj_in_memory]] to 0 if optdriver is not GS.

See [[test:paral_19]], [[test:paral_20]], [[test:v10_19]], [[test:v10_24]], [[test:v10_25]].

By L. Baguet (MR1067, 1108, 1133)

**D.3** Test cases for [[dilatmx]] with chebfi or lobpcg

See [[test:paral_33]] and [[test:paral_34]].

By L. Baguet (MR1122)

**D.4** Change [[nbdbuf]] for all testcases with chebfi2 to clearly show that band residuals of occupied bands are as good as with lobpcg.

By L. Baguet (MR1122)


**D.5** Math operation beautification (from Abidev24 Hackathon).

The goal was mainly to improve the readability of the code through the use of Fortran intrinsic math functions and vectorization capabilities.

By M. Mignolet, F. Gomez-Ortiz and I. Lygatsika (MR1083)


**D.6** Developments and bug fixes in extended First-Principles Molecular Dynamics (high temperature).

The extended First-Principles Molecular Dynamics model now allows values of [[occopt]] different from 3 by setting the electronic temperature during the initialization (metallic occupations are still required).

By A. Blanchet (MR1135).

**D.7** MRGDDB polishing.

Allow for single DDB file in command line interface for DDB conversion.
Remove unwanted output file fort.7.

By G. Antonius (MR1132). 


**D.8** Spinat check

Reestablish the test on [[spinat]] and update input and ref files [[test:v5_16]]

By Le Shu (MR 1093)


**D.9.** Dynamical allocation of the input string.

Now, INPUT_STRING is allocated dynamically instead of using strlen=2000000.
This reduces the size of the NetCDF files by ~2 MB and also decreases the memory footprint
per node when all cores on the node are used as MPI processes.

By M. Giantomassi (MR1109)

**D.10.** Optic utility : antiresonant terms and band decomposition 

Anti-resonant contributions were missing in the Optic calculations of second-harmonic generation.
Now they can be included. Use  do_antiresonant = .FALSE. parameter in the input file of optics (meaning that the antiresonant approximation is NOT done - this is confusing and should be changed).

Add option to print the contribution of each bands-triplet to the SHG at a given frequency.

By V. Trinquet (MR1143)


**D.11** Various improvements to Conducti in PAW.

Among others, the code can now handle different atomic types in the unit cell and compute the reflectivity at any angle.
A complete userguide is now available, and some bugs were fixed.

By J. Boust (MR1152)

* * *

### **E.**  Bug fixes, not yet mentioned in the previous sections of these release notes.

**E.1** The [[prtevk]] variable was not read from input file. Now fixed.

By M. Torrent (MR1141)


**E.2** The syntax for $(<.current_version) might not work properly in some older versions of sh, making $(cat .current_version) a more stable alternative.

By Hsiaoyi Tsai (MR1114)


**E.3** Using libyaml.a instead -lyaml if bigdft using internal yaml 

In fallbacks, bigdft will install its internal yaml library if the system does not have the yaml library. However, this is not shown in the dynamic linking during the execution of abinit; using static linking can avoid the issue of missing libraries.

By Hsiaoyi Tsai (MR1138)


**E.4** Fixing langevin_pimd moldyn algorithm

The recently implemented Langevin PIMD was encountering NaN values when Abinit was compiled with specific optimization options. This issue was fixed by initializing arrays.

By A. Blanchet (MR1094)


**E.5.** Tolerances have been cleaned in tests/Input/*.abi files.

By M. Azizi (MR1086)


**E.6** Miscellaneous bug fixes and cleaning

* Correct [[wfoptalg]] initialization for GPU (no more [[wfoptalg]]=14 for [[gpu_option]]>1)
* Correct [[test:v10_11]] test adding more bands so LOBPCG runs are correct
* Include [[nbdbuf]] in the check of the last occupancy in vtorho (and correct the exit condition for [[nsppol]]=2)
* xmpi_distrib_2D: fix edge case where the code would try to distribute n points over N>n processes
* m_krank: add support for non-diagonal [[kptrlatt]] in krank_from_kptrlatt()

By L. Baguet (MR1127), M. Mignolet (MR1118)

* * *








## v10.2

Version 10.2, released on October 10, 2024.
List of changes with respect to version 10.0.
<!-- Release notes updated on Feb 15, 2025. -->

Many thanks to the contributors to the ABINIT project between
March 2024 and October 2024.
<!-- (with some late contributions/fixes until Feb 2025). -->
These release notes
are relative to modifications/improvements of ABINIT v10.2 with respect to v10.0, with also some late fixes to v10.0, after March 2024,
that had not yet been documented.
<!-- Initially, merge requests from MR986 up to and including MR1062 in the beta-release v10.2.1. Also, MR1034, MR1064 to MR1067, MR1069 to MR1072, MR1074 are taken into account in v10.2.3, publicly released. For v10.2.5, MR1075, 1078 and 1087 are also included. For v10.2.7, MR1088, 1090, 1091, 1092, 1095, 1097, 1100, 1104, 1106, 1114, 1115 are also included.-->

The list of contributors includes:
G. Antonius, M. Azizi, L. Baguet, J.-M. Beuken, O. Bistoni, A. Blanchet, F. Bottin, F. Bruneval, Siyu Chen, F. Gendron, M. Giantomassi, X. Gonze,
P. Kestener, L. MacEnulty, M. Mignolet, C. Paillard,
S. Ponce, M. Royo, M. Sarraute, M. Torrent, V. Vasilchenko, M. Verstraete, Xu He, J. Zwanziger

It is worthwhile to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases.
This might take some time ...

Xavier

* * *

### **A.** Remarks and warnings.

**A.1** Change [[nline]] default to 6 for ChebFi algorithm .

By M. Torrent (MR1069)

**A.2** Renaming tests/Psps_for_tests to tests/Pspdir

By X. Gonze (fc5345534a and subsequent commits)

* * *


### **B.** Most noticeable developments

**B.1** GPU porting of the DFPT driver of ABINIT

The global GPU porting of ABINIT using recent libraries/compilers, started three years ago, has been continued.
In the previous release 10.0, two implementations (OpenMP or KOKKOS+CUDA) for ground-state calculations [[optdriver]]=0
had been made available.

In the present release 10.2, the Density-Functional Perturbation Theory [[optdriver]]=1
has been ported, using OpenMP, along with the following perturbations:
- phonons ([[rfphon]])
- electric field ([[rfelfd]])
- strains ([[rfstrs]])
Numerous tests are available, [[test:gpu_omp_11]] to [[test:gpu_omp_25]]

See the description of the GPU possibilities of ABINIT in the documentation, input variable [[gpu_option]]=2 for the OpenMP capabilities.
For the description of the modifications of ABINIT, see the Merge Requests (MR) below.

By M. Sarraute and M. Torrent (MR 1027, MR1055, MR1059, MR1071)


**B.2** GPU porting of Fock calculations

Similarly, the Fock calculation has been ported to GPU, using OpenMP.
Tests are available, [[test:gpu_omp_26]] to [[test:gpu_omp_29]].
See the description of the GPU possibilities of ABINIT in the documentation, input variable [[gpu_option]]=2.

For the description of the modifications of ABINIT, see the MR below.

By M. Sarraute and M. Torrent (MR 1062)


**B.3** Interface to coupled-cluster CC4S calculations

The file needed as input for computations with the CC4S package, <https://manuals.cc4s.org/user-manual>,
allowing to perform coupled-cluster calculations (e.g. CCSD), perturbative triples (and more),
can be produced using [[optdriver]]=6 and [[gwr_task]]="CC4S".
This functionality was available in v10.0. Now it can be considered in production, as it has been tested reasonably and debugged.
Moreover, there is a new option [[gwr_task]] = "CC4S_FROM_WFK" to read orbitals from file instead of using direct diagonalization.
The support for [[istwfk]]=2 has been added.
See tests: [[test:gwr_suite_07]], [[test:gwr_suite_08]], [[test:gwr_suite_09]], [[test:gwr_suite_10]] and [[test:gwr_suite_11]].

By M. Giantomassi (MR 1006) and F. Bruneval (MR1024).

**B.4** Speed-up of the calculation (PAW as well as norm-conserving), thanks to cprj_in_memory.

The coefficients of the wavefunctions giving the in-sphere contribution, denoted "cprj", can now be stored in memory,
thanks to the input variable [[cprj_in_memory]].
LOBPCG ([[wfoptalg]]=114) and Chebyshev filtering ([[wfoptalg]]=111) can be accelerated using this option.
There are many tests of this functionality, [[test:v10_10]] to [[test:v10_19]], [[test:paral_35]] to [[test:paral_40]],
 [[test:paral_48]] and  [[test:paral_49]], and  [[test:v9_71]] to  [[test:v9_74]].
However, as usual, the users are advised to check for themselves the results from turning on this new implementation (still mentioned as being "in development").

Many parts are now computed with matrix-matrix multiplications, done with BLAS routines ("gemm").
As a consequence, the efficiency of the implementation strongly depends on the BLAS library, which has to be multi-threaded if used with OpenMP.

Available only for the ground state. See the other limitations in the documentation of [[cprj_in_memory]].

It is available for both PAW and norm-conserving (NC) pseudopotentials, but with some constraints for the latter.
Indeed, NC can be used only with [[useylm]]=1 (set automatically if [[cprj_in_memory]]=1).
Also [[nspinor]]=2 does not work properly (but works well for PAW). This case is forbidden at the moment (test paral_50 is not activated).

At last: ABINIT also contains a small fix of [[tolwfr_diago]], which was not used in lobpcg ([[tolwfr]] was used instead).

In the case of Chebyshev filtering, two implementations are available, controlled by [[xg_nonlop_option]] (see [[test:paral_35]], [[test:paral_38]]).
In practice, only the default option has good performance.

By L. Baguet (MR 1058, MR1067, MR1070, MR1072)

**B.5** Wannier90 with spinors, and other Wannier90 improvements.

There has been several improvements of the Wannier90 usage together with ABINIT: enable w90 output with spinor wave functions, improvements to parallel calls,
general cleaning of mlwf routines.
See  [[test:wannier90_05]] and  [[test:wannier90_14]].

By M. Verstraete  and He Xu (MR1005).

**B.6** New developments and bug fixes in extended First-Principles Molecular Dynamics (high temperature).

Now ABINIT can compute extFPMD high energy contributions with a discrete sum on the orbitals.
This is activated with [[useextfpmd]]=5, and controlled with the new input variable [[extfpmd_nband]].
This allows ABINIT to compute smooth DOS between Kohn-Sham eingenstates and extended ideal electronic states.
See  [[test:v9_92]].

Also, [[useextfpmd]]=4 has been fixed, by computing the band shift. Now [[useextfpmd]]=1 and [[useextfpmd]]=4 show very similar results (not exactly equal though).

The value of extfpmd_eshift has been added in the header of wave function files (made optional in order not to break post processing tools): this allows for a faster restart of a calculation from the wavefunctions.

By A. Blanchet (MR1013).

* * *


### **C.** Changes for the developers or for the installation of ABINIT (including information about compilers)


**C.1** The "MINIMAL" tag has been added to a selected set of automatic tests (about 20 abi files).
They cover the basic functionalities of ABINIT, and the detection of failures is modified and made more robust.
There should not be false failures due to tight tolerances.
So, this should be the ideal set of tests to be executed after a new installation.
Use "runtests.py -k MINIMAL". Also, "make check-local" and "make check-am" trigger this suite.

By M. Verstrate, M. Torrent and X. Gonze (MR1035, MR1041)

**C.2** Fix support for new compilers: ifort / ifx oneAPI 2024 ( replace some implicit loop by an explicit one ) and improve detection of cray compilers.
In addition, oneapi and newer intel compilers refuse the old syntax -mkl=cluster, moved to -qmkl. One can presume that all future compilations will be done with this.
So, this has been modified in config/m4/sd_math_linalg_core.m4 .

By J.-M. Beuken (MR995) and M. Verstraete (MR1007)

**C.3** The build system has been modified in order to include cmake files in 'make dist'. Some modifications needed to compile with cmake have been made.

By M. Torrent (MR986, MR1032)

**C.4** There is an embryo documentation for cmake build. Also, the need to modify CMakeLists.txt is not mentioned in doc/developers/developers_howto.md .

By P. Kestener (MR994), and a comment from Maxime Mignolet.

**C.5** Manage GPU markers with compilation flag rather than input parameters.

By M. Sarraute (MR989)

**C.6** Improve NVPL lib support in build system

By M. Torrent (MR1069)

**C.7** Upgrade of mksite.py (abimkdocs).

Make mksite.py compatible with mkdocs-material==9.5.43, mkdocs==1.6.1 and python3.12

By M. Giantomassi (MR1078), with checking by M. Azizi of possible glitches in the web pages.

* * *


### **D.**  Other developments (possibly not yet finalized), other new tests, new input variables, new tutorial.

**D.1** The work on the ground-state GPU porting of ABINIT is continuing.

Concerning eigensolvers on AMD GPU:
By eigensolver routines, one refers to LAPACK DSYGVD/ZHEGVD and DSYEVD/ZHEEVD, which are provided on GPU by AMD HipSOLVER and NVIDIA CuSOLVER under specific names.
As of today, DSYGVD and DSYEVD are significantly slower on AMD GPU, and usually account for most of execution time in Adastra GPU partition.
But somehow, the complex variants ZHEGVD and ZHEEVD are incredibly slow, which lead to introduce a workaround,
use ZHEGVJ instead of ZHEGVD: J stands for Jacobi and it takes two extra parameter for tolerance and max step.
The maximum step number was raised from 1 to 100, as numerical issues would sometime occur otherwise. Performance impact seems light.
At varianes, ZHEEVD is used on CPU: AMD HIP provides a ZHEEVJ it was not possible to make it work.
Running ZHEEVD on CPU is slow... but faster than AMD HIP ZHEEVD

Concerning GPU-aware MPI:
This concept consists of providing GPU buffers addresses to MPI routines. If GPU-direct is enabled, data will transmit between GPU interconnect, aka NVlink (NVIDIA) or InfinityFabric (AMD), UAlink (next-gen AMD/Intel). In ABINIT, it can be enabled using enable_mpi_gpu_aware=yes  in Autotools buildsys (CMake buildsys uses autodetect).
This usually bring performance improvements, if GPU-direct is enabled, works fine on NVIDIA A100... but works erratically on AMD MI250.
It is only used to perform a MPI_SUM in xg_RayleighRitz, which is a performance bottleneck on GPU.
Anyway, this section was broken for some usecases ([[istwfk]]==2, AMD GPU) in xg submod so it has been fixed.

In addition, there were the following fixes:
fix unprotected OpenMP TARGET directives;
fix compilation with NVHPC 23.11-24.7 by switching mkrho optionals positions.

By M. Sarraute (MR1056, MR1060, MR1061)


**D.2** Tests of the GPU (KOKKOS+CUDA) porting of ABINIT are available,  [[test:gpu_kokkos_01]] and [[test:gpu_kokkos_02]], nonlop with BLAS, chebycheff.
See input variable [[gpu_option]]=3.

By M. Sarraute and M. Torrent (b3b2ce14ca1)


**D.3** New tests of the GPU (OpenMP) porting of ABINIT are available, [[test:gpu_omp_04]] to [[test:gpu_omp_09]]. See input variable [[gpu_option]]=2.
Bug fix:  [[istwfk]]==2 when using CHEBFI as diago algorithm.
Small fixes around ELPA+GPU, Cray and GPU-aware MPI have been done.

By M. Sarraute and M. Torrent (MR1002, MR1004)


**D.4** New developments and tests related to the chemical shielding, orbital magnetization, nuclear dipole moment.

PEAD wavefunctions have been added to orbital magnetism capabilities:
Can call [[orbmag]] and [[berryopt]]=-2 together, to compute orbital magnetic moment using PEAD derivative wavefunctions (instead of DFPT-based DDK wavefunctions).

Spatial symmetries ([[kptopt]]=4) have been added to nuclear dipole moment in ground state.
kptopt 4 now allows for systems with nuclear dipole moments ([[nucdipmom]]).
Usually not much symmetry is left once a single nuclear dipole moment has been added, but there may be some.
See [[test:v10_40]].

The requirement for m_orbmag routine to have [[pawxcdev]] 0 has been lifted. The default 1 works fine also, and the flexibility in pawxcdev makes the m_orbmag routine compatible with more features (mgga related, hybrid functional related, etc)

The fock_getghc code has been corrected to allow for cplex_dij = 2 case; necessary for use with [[nucdipmom]]. Now for example pbe0 calculations run together with nuclear dipole moments, allowing NMR shieldings to be calculated in the hybrid functional case.

By J.Zwanziger (MR993, MR1016, MR1051, MR1052)


**D.5** Miscellaneous developments for ground state.

This includes the rework of [[istwfk]]=2 with xg_tools (both CPU and GPU) and the extension of [[istwfk]]>2 to lobpcg2 and chebfi2
(but works only with DFTI - the FFT from MKL).
Tested in [[test:v10_03]], [[test:v10_04]], [[test:paral_44]] and [[test:paral_45]] for DFTI.
Same inputs are tested without DFTI in [[test:v10_05]], [[test:v10_06]], [[test:paral_46]] and [[test:paral_47]], but they use only [[istwfk]]=1 and 2.
The Chebicheff filtering algorithm has been generalized to PAW with Spin-Orbit Coupling.
See [[test:paral_44]] to [[test:paral_47]], as well as [[test:v10_03]] to [[test:v10_06]].

By L. Baguet (MR 1014,1021,1022)

**D.6**
Authorize Hartree-Fock with time-reversal symmetry when a single k-point is used.

By F. Bruneval (MR1038)

**D.7** Migration of a set of files previously present in the Dokuwiki of ABINIT, to the inner documentation of ABINIT. The dokuwiki has been suppressed,
as all its content is migrated either to the main ABINIT package, or to the Web site.

By M. Azizi (MR1046)

**D.8** Write magmoms to GSR.nc file in order to re-train chgnet.

By M. Giantomassi (MR1047)

**D.9** Add documentation for [[eph_task]]=14, and small correction to release notes.

By M. Mignolet (MR990)

**D.10**
Add new input variables [[getabiwan]] and [[getabiwan_filepath]].
Define the name and path of the ABIWAN.nc file with the Wannier rotation matrices produced by ABINIT
when wannier90 in invoked in library mode.
This is documented, but there is no test available yet as this feature is still under development.

By M. Giantomassi (ac10911a)


**D.11**
Add new input variables [[getgwan]] and [[getgwan_filepath]].
Define the name and path of the GWAN.nc file with the e-ph matrix elements in the Wannier representation
when wannier90 in invoked in library mode.
This is documented, but there is no test available yet as this feature is still under development.

By M. Giantomassi (5c6ff934)


**D.12** Add new input variables [[getdrhodb]], [[getdrhodb_filepath]], [[irddrhodb]].
to read the database with the first-order densities for GWPT calculations.
There is no test available yet (except the fake one `test:gwpt_04`), as this feature is still under development.

By Siyu Chen (867b11f9) and documentation of input variable by M. Giantomassi (MR1066)


**D.13** The post-process mode of Wannier90 with wfk_task "wannier" is now tested. Upgraded [[test:w90_7]] for that purpose.

By M. Giantomassi (5c5ff934b4c) on top of HeXu initial [[test:w90_7]] (63fb2a6be73).


**D.14** Very preliminary work to support wannier-based interpolation of e-ph matrix elements.

By M. Giantomassi (MR1047)

**D.15** Acceleration of the 'atdep' test series and normalization of eigenvectors.

By F. Bottin (MR1048)


**D.16** GWPT implementation is on-going. See [[eph_task]]=17. Testing has been initialized, see test:gwpt_01 to test:gwpt_04
New input variable [[gwpt_np_wpqbks]].
Note that for the time being, gwpt is commented in tests/\_\_init\_\_.py

By M. Giantomassi (5b9bb0e88ab, 306ccc446a2, c03fa967bde)

**D.17** New test for mGGA Norm conserving, [[test:v9_95]].
Copper atom isolated in a box.
Tests the reading of the upf2 potential, and inputs related to the model core kinetic energy density, then their usage.
[[nc_xccc_gspace]] 0 calculates model core contribution in real space, and is much cleaner (0 kinE and density outside NLCC radius).
Density (+kinE) mixing works much better than potential mixing: in the vacuum the mGGA potential is very sensitive.

By M. Vertraete  (43214101824)

**D.18** DMFT improvements.
Fixing size name of susceptibility files in DMFT.
Include local charge susceptibility in 62_ctqmc and results for tau = beta.

By F. Gendron (MR987 and MR1030)

**D.19** Add documentation for [[eph_task]]=14, and small correction to release notes.

By M. Mignolet (MR990)

**D.20** The variational polaron equation methodology
Ongoing implementation of the methodology described in [[cite:Vasilchenko2022]], activated with [[optdriver]]=7 and [[eph_task]]=13.
Current commit containts intermediate developments & variables, which are used for testing purposes only and are to be modified or removed later on.
Final release of the framework is expected in the next version of the code.
Related input variables: varpeq_aseed, varpeq_erange, varpeq_gau_params, varpeq_interpolate, varpeq_orth, varpeq_nstep, varpeq_pkind, varpeq_tolgrs, varpeq_pc_nupdate, varpeq_pc_factor, getvarpeq, getvarpeq_filepath.

By V. Vasilchenko, with help from M. Giantomassi (MR 1047, MR1065).

**D.21** A new tutorial, on the computation of dynamical quadrupoles, has been created.
Also, the elph legacy tutorial has been updated.

By S. Ponce, with contribution from M. Royo (MR1034).

**D.22**

Clickable version of the tutorial flowchart

By M. Azizi (MR1076, backported to ABINITv10.2)


* * *

### **E.**  Bug fixes, not yet mentioned in the previous sections of these release notes.


**E.1** Fix bugs in the recognition of symmetries (trigonal and hexagonal groups).

Work based on a user's report of a Wyckoff position error in space groups 151 and 153.
Update symsghexa and update relevant tests. Provide a new [[test:v10_02]] that tests all Wyckoff positions in groups 143-167, the trigonal groups.
See also [[test:v10_40]] and [[test:v3_98]].
By J.Zwanziger and M. Verstraete (MR1019, MR1025)


**E.2** Fix an error when using ABI_GPU_LEGACY with xg_tools ([[wfoptalg]]=114 or 111).

By L. Baguet (MR988)


**E.3** Fix: velocity was not read from hist file with [[restartxf]].

By He Xu (MR991)

**E.4** Multibinit fixes: multibinit spin monte carlo did not run; multibinit did not read more than one single-ion anisotropy from xml file.

By He Xu (MR992)


**E.5** Correct some small defects in the GPU version (test residue).

By M. Torrent (MR996)

**E.6** Pawrhoij fixing: Correct indexing in write/read of pawrhoij datastructures to POT and DEN files when qphase=2;
Correct indexing in storing of pawrhoij datastructures after rotation in the symmetrization subroutine when qphase=2 .

By Oliviero Bistoni (MR998)

**E.7**
Fix race condition in inwffil (cg_envlop kernel) occurring when ABINIT is compiled with -O0  with OpenMP parallelism enabled.
-O0  might be used for debugging and it is always handy when code numerically behaves the same regardless of optimization flags.

By M. Sarraute (MR1023)

**E.8**
Fix wrong start and count in [[eph_task]]=-2 when [[nsppol]] = 2.
Create indata/tmpdata directories if dirs do not exist.
Set OMP_NUM_THREADS to 1 if env var is undefined

By M. Giantomassi (MR124)

**E.9** Fixed a nasty bug in m_fftcore.F90 . The remainder != 1 check must be performed outside the loop. Otherwise, for a valid number of nproc_fft the code throws an error and stops.

By V. Vasilchenko (MR1028)

**E.10** Fix:Virtual Crystal Approximation and [[useylm]]=1.
In 64_psp/m_pspini.F90, filling the ffspl array when [[useylm]]=1 lead in certain cases to segmentation fault due iln set to ilmn,
thus running beyond ffspl dimension psps%lnmax.
The filling of ffspl has been revised, as well as the indexation in psps%indlmn (see lines 514 - 554, where most of the changes were made).
Similar index problems in 66_nonlocal/m_mkffnl.F90, that has thus been modified.
Now, one should be able to run non-linear response DFPT calculations with VCA, such as the electro-optic response.

By C. Paillard (MR1033)

**E.11**
A bug had been introduced in Abinit v9.8 when updating the format of _EIG2D file to make them real DDB files that can be merged.
In the case of [[ieig2rf]] > 2, only the Sternheimer part of the EIG2D should be written to file. This was achieved using the
outbsd routine by copying the eig2nkq array into eig2nkq_tmp in the eig2tot subroutine.
When the outbsd routine was replaced by ddb%set_d2eig_reshape(), the eig2nkq was passed instead of the correct eig2nkq_tmp one,
resulting in a bug --> the fully EIG2D was printed to file instead of the partial Sternheimer only one.
See documentation of [[ieig2rf]] for explanation.
This bug is corrected here. In addition, the coding was not optimal. The content of the former array eig2nkq_tmp could be printed directly in m_eig2d/eig2stern,
instead of inside m_dfpt_looppert, avoiding the allocation of this large array, and mirroring the same printing done in eig2tot.

Fix EIGR2D files for the case [[ieig2rf]]=5.
Move writing of EIGR2D files in eig2stern.
Add a test comparing an EIGR2D file against reference ([[test:eph_tdep_legacy_5).

The electron-phonon legacy scripts have been put back in the package, and have been updated to work with Python 3.6 and to support the new EIGR2D.nc DDB file format.

By S. Ponce (MR1036, MR1045) and G. Antonius (MR1043), with help from X. Gonze.

**E.12**
AbiDev 2024 hackathon tutorial fixes, specifically to the lruj tutorial, the eph_legacy and the DFPT2 modifications.
LRUJ tutorial fixes and acknowledgements changes.

By L. MacEnulty (MR1044)


**E.13** Fix regression in eph_driver when DDB file contains dynamical quadrupoles.

By M. Giantomassi (MR1047)


**E.14** Fix portability issues in GWR code.
Fix SIGSEGV when calling PZLATRA, likely due to ABI mismatch between gcc and intel scalapack
in functions returning complex results. Replace PZLATRA calls with home-made implementation.
Fix wrong index in slk_get_trace.
Allocate dummy_vec scalapack matrix before calling heevm and here the compiler is too picky because dummy_vec is never referenced when computing eigenvalues only!

By M. Giantomassi (MR1049, MR1050)


**E.15** Fix internal mapping of k-points when [[kptopt]] = 2 or 3.

m_kpts:kpts_map was always assuming spatial symmetries for mapping one set of k-points on another. This was causing gstore_compute to crash when trying to use a _WFK file generated with [[kptopt]]=3, and with [[gstore_kzone]]="bz".

By M. Mignolet (MR1057)

**E.16** New NVHPC version 24.9 dropped recently and was failing to compile (again). Switching loop order simply fixes the thing.

By M. Sarraute (MR1064)

**E.17** Minor fixes: fix compilation with gcc8.5 regarding contiguous pointers; Fix gpu[1], by fixing condition restricting GPU to DFPT, only with OpenMP; Reduce implicit copies when nspinor==1 in GEMM nonlop.

By M. Sarraute (MR1074)

**E.18** Several fixes

Bug 1: random errors on bob (only bot with threads) for some test cases (in fact: cases with bandpp>1).
It was due to calls of goedecker fft (non thread-safe) by multiple threads in the routine multithreaded_getghc.
Bug 2: error on eos for paral[45]
It was due to uninitialized variables in m_slk module.
Bug 3: random compiler errors with cray. Modifications done to compile and run with fortran cray compiler.

By L. Baguet (MR1087)

**E.19** Miscellaneous bug fixes and cleaning

By He Xu (MR1009), by M. Verstraete (MR1025)

* * *


## v10.0

Version 10.0, released on March 18, 2024.
List of changes with respect to version 9.10.
<!-- Release notes updated on Oct 2, 2024. -->

Many thanks to the contributors to the ABINIT project between
May 2023 and March 2024.
<!-- (with some late contributions until 11 Sep 2024). -->
These release notes
are relative to modifications/improvements of ABINIT v10.0 with respect to v9.10, with also some late fixes to v10.0, after March 2024,
that had not yet been documented.
<!-- TO BE CHANGED Merge requests up to and including MR984. Also, take into account: MR985 (backported from v10.1), MR999 to MR1003, MR1008, MR1011, MR1015, MR1018, MR1020, MR1021, MR1023, MR1031, MR1037, MR1039, MR1040, MR1053, MR1054 -->

The list of contributors includes:
G. Antonius, L. Baguet, J.-M. Beuken, F. Bottin, J. Bouchet, J. Bouquiaux, F. Brieuc, A. Donkov, F. Gendron, M. Giantomassi, X. Gonze,
F. Goudreault, B. Guster, P. Kestener, M. Mignolet,
S. Ponce, M. Sarraute, M. Torrent, V. Vasilchenko, M. Verstraete, J. Zwanziger

It is worthwhile to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases.
This might take some time ...

Xavier

* * *

### **A.** Important remarks and warnings.

**A.1** The input variable rfasr, used in DFPT, in case of phonon perturbations and electric field perturbations has been replaced
by [[asr]] and [[chneut]]. The latter variables have been used already for some time for a more detailed imposition of the acoustic sum rule
and the charge neutrality sum rule in the electron-phonon part of ABINIT, and even for a longer time
in the ANADDB utility. Actually, rfasr was used to initialize [[asr]] and [[chneut]] internally.
The default values for [[asr]] and [[chneut]], namely 1, are however not the same as the previous default value of rfasr, namely 0.
Thus a large fraction of the reference files of the test have been upgraded to the new default. However, in a sizeable number of reference files
the old default value has been specified explicitly.

**A.2** The input variables rf1atpol, rf1dir, rf1elfd and rf1phon (and similar input variables for perturbations 2 and 3) have been suppressed.
They were used in the case of Raman calculations (one derivative with respect to an atomic displacement, and two displacements with respect to an electric field),
but have been superceded by rf2_XXX input variables a long time ago.

* * *

### **B.** Most noticeable achievements

**B.1** Low-scaling GW and RPA implementations

A cubic scaling real-space imaginary-time algorithm for GW and RPA has been implemented.
See the theory in [[cite:Liu2016]] and related references. It relies on the minimax time-frequency grids
available in the GreenX library [[cite:Azizi2023]].
At present, only norm-conserving pseudopotentials can be used. The implementation is restricted to non-magnetic materials,
and without spin-orbit coupling.
Still, different types of flows and algorithms (including self-consistency) are available, see [[gwr_task]].

Activate it using [[optdriver]]=6, and specify [[gwr_task]].

New input variables: [[gwr_task]], [[gwr_chi_algo]], [[gwr_sigma_algo]],
[[gwr_np_kgts]], [[gwr_ucsc_batch]], [[gwr_ntau]],
[[gwr_boxcutmin]], [[gwr_max_hwtene]], [[gwr_rpa_ncut]], [[gwr_nstep]]. Also, [[gwr_tolqpe]] replaces the obsolete gw_toldfeig input variable.

A introductory tutorial (work in progress) is available, see [[tutorial:gwr_intro]].
Tests are provided in the newly created subdirectory tests/gwr_suite, see test:gwr_suite__01 to test:gwr_suite_08
See also test:paral_78 and test:paral_79

By M. Giantomassi (MR875, MR907, MR964)

**B.2** GPU porting of ABINIT

While there was already an old CUDA-based port of ABINIT on GPU, two new implementations are provided in ABINITv10 .
One implementation is based on OpenMP, the other one is based on KOKKOS+CUDA.
At present, this is only available for ground-state calculations [[optdriver]]=0.
See the description of the GPU possibilities of ABINIT in the documentation, input variable [[gpu_option]].

The
OpenMP implementation works on NVidia accelerators, if ABINIT has been
compiled with a [[CUDA]] compatible compiler and linked with NVidia FFT/linear algebra
libraries ([cuFFT](https://docs.nvidia.com/cuda/cufft),
[cuBLAS](https://docs.nvidia.com/cuda/cublas) and
[cuSOLVER](https://docs.nvidia.com/cuda/cusolvermp)).
It also works on AMD accelerators (EXPERIMENTAL),
if ABINIT has been compiled with a AMD compatible compiler and linked with NVidia
FFT/linear algebra libraries ([ROCm](https://www.amd.com/fr/graphics/servers-solutions-rocm)
or [HIP](https://github.com/ROCm/HIP)).

The [[KOKKOS]]+[[CUDA]] implementation -- at present -- is only compatible with
NVidia accelerators. It requires that ABINIT has been linked to the
[Kokkos](https://github.com/kokkos/kokkos) and [YAKL](https://github.com/mrnorman/YAKL)
performance libraries. It also uses NVidia FFT/linear algebra libraries
([cuFFT](https://docs.nvidia.com/cuda/cufft), [cuBLAS](https://docs.nvidia.com/cuda/cublas)).
The [[KOKKOS]] GPU implementation can be used in conjuction with openMP threads
on CPU (see [[gpu_kokkos_nthrd]]).

For an expert user of ABINIT on [[GPU]], some additional keywords can be used. See [[gpu_nl_distrib]], [[gpu_nl_splitsize]].

Note: the input variable use_gpu_nvtx has been suppressed and replaced by an option to
be set at configure step: --enable-gpu-nvtx.

Several GPU devices can be detected and used on a node.

GPU regression tests are present in new directories,
[[test:gpu_omp_01]], [[test:gpu_omp_02]], [[test:gpu_omp_03]],
[[test:hpc_gpu_omp_01]] to [[test:hpc_gpu_omp_13]],
as well as in the previously existing directory gpu.

The ABINIT test farm has been upgraded to support those tests. The NVHPC compiler is supported. Also, ABINIT works on Adastra cluster with newest Cray version (23.12).

FOR THE TIME BEING (APRIL 2024), THE DOCUMENTATION TO BUILD THE GPU VERSION OF ABINIT REMAINS TO BE WRITTEN.
Still, there are three examples *ac9 files in the directory doc/build of the package, that are also available from
[https://github.com/abinit/abinit/tree/master/doc/build](https://github.com/abinit/abinit/tree/master/doc/build) .
Moreover, if you want to be a beta tester and collaborate on this topic, please contact Marc Torrent.

By P. Kesterner, M. Sarraute, J.-M. Beuken, L. Baguet and M. Torrent
(MR942, 943, 951, 954, 955, 961, 965, 966, 967, 968, 969, 974, 978)


**B.3** CMake build of ABINIT

ABINIT can now be build using CMake instead of the standard configure+make. This is needed to build the GPU version
of ABINIT relying on KOKKOS. Try:

mkdir build; cd build; cmake ..

More information is available in the ABINIT [[help:../installation|installation guide]].
Still, the usual build procedure, using autotools is to be preferred by non-experts. CMake is to be considered experimental.

Also, the version number of ABINIT is now generated automatically from the git tag information.

By P. Kestener with M. Torrent (MR944, 979)


**B.4** Computation of phonon angular momentum

It is now possible to compute the angular momentum of phonons in anaddb, following the formula provided in [[cite:Zhang2014]]. A new type of file (*_PHANGMOM) has been defined. The phonon angular momentum is also written to *_PHBST.nc. See [[test:rf2_5]].

By M. Mignolet (MR921)


**B.5** New input "supra"variable [[write_files]]

It is now possible to govern the printing of files thanks to the "supra"variable [[write_files]],
instead of using the numerous prt* input variables.

This supravariable is introduced while maintaining the underlying logic of the prt file options inside abinit.
(Experienced) Users can now trigger the presence or absence of a file in a calculation via a string using

[[write_files]] ``<list_of_file_suffixes_with_options\>"

When rationalizing the set of prt* variables and their behavior, new ones were introduced: [[prtevk]], [[prthist]] and [[prtddb]].

See [[test:v9_150]].

By B. Guster (MR920)


**B.6** New capabilities of mrgddb, and netcdf reading/writing for pseudopotential and crystal types.

A new DDB IO interface has been coded.
The netcdf format can be activated with [[iomode]].
The postprocessor mrgddb can convert a single DDB file between text format and netcdf format
Also the IO capabilities have been extended for pseudopotential_type and crystal_t
See [[test:v9_160]] and [[test:v9_161]]

By G. Antonius (MR 927)


**B.7** Improvements for metaGGA.

There are several new input variables related to metaGGA.
The new variable [[xc_taupos]] allows one to control positivity of kinetic energy density, regardless the one for the density.
The new variables [[irdkden]] and [[getkden]] allow one not to read KDEN when reading DEN (and other possible applications).
There are several other improvements for automatic settings (forces, stresses, self-consistent cycle).
The computation of c parameter of TB09 functional, compatible with NCPP and PAW, is now implemented in a specific routine.

By the way, mBJ/TB09 is not adapted to forces/stresses computation because it does not solve a variational problem for the energy. This means that it is probably not suitable for DFPT.

See the new tests [[test:libxc_23]], [[test:libxc_24]]

By M. Torrent (MR 938)

**B.8** Spin-orbit for orbital magnetism.

The spin-orbit coupling has been added to the computation of orbital magnetism.
See [[test:v9_141]].

By J. Zwanziger (MR980)

**B.9** Interface to coupled-cluster CC4S calculations.

The writing of the file needed as input for computations with the CC4S package, <https://manuals.cc4s.org/user-manual>,
allowing to perform coupled-cluster calculations (e.g. CCSD), perturbative triples (and more), can be activated using [[optdriver]]=6 and [[gwr_task]]="CC4S".
See test test:gwr_07

By M. Giantomassi (MR 875, 907)

* * *

### **C.** Changes for the developers (including information about compilers)

**C.1** A new bot called EOS has been included in the test farm, in order to test ABINIT on GPUs.
Several flavors of the NVHPC compiler are available.

By J.-M. Beuken, M. Torrent, M. Sarraute, P. Kesterner (MR949)


**C.2** New bots have been introduced to replace obsolete ones:
bob_gnu_13.2_openmp, higgs_gnu_12.3_cov, scope_gnu_13.2_dep.

By J.-M. Beuken (MR960)


**C.3** Support for Py3.12 has been added to ./runtest.py .
Use importlib if py>3.12, this fixes SyntaxWarning due to invalid escape sequence.

By M. Giantomassi (MR958)

**C.4**
Fixes to make Abinit compile when fft_flavor=fftw3-threads and openMP.

By M. Torrent (MR970)

* * *

### **D.**  Other changes (or on-going developments, not yet finalized, as well as miscellaneous bug fixes)


**D.1**
There are several improvements related to [[tolwfr]].

Now one can couple [[tolwfr]] with other tolerances.
In that case, the SCF loop stops if both criteria are satisfied. The documentation has been upgraded accordingly. See [[test:v9_200]].

The input variable [[tolwfr_diago]] has been added, to distinguish the criterion used for SCF loop
and the one used inside diagonalization algorithms to skip lines.
[[tolwfr_diago]] is set to [[tolwfr]] by default. Note that one can define [[tolwfr_diago]] while [[tolwfr]] is not defined.

There is a new use of [[nbdbuf]] (=-101), which is kind of an automatic buffering.
In that case, the maximum of residuals is computed as max(res*occ) instead of max(res).
Used for both [[tolwfr]] and [[tolwfr_diago]]. This is experimental, so not documented yet.

Apart from LOBPCG ([[wfoptalg]]==114), these developments are new features without altering the previous code behavior, as shown in the testsuite.
See [[test:mpiio_26]], [[test:mpiio_27]], [[test:mpiio_51]],
[[test:paral_31]],[[test:paral_63]],
[[test:paral_66]],[[test:paral_86]]
[[test:psic_03]], [[test:v9_202]].

By L. Baguet (MR947, MR1070))

**D.2**
Introduced the new [[nblock_lobpcg]] variable (default=1), and set [[bandpp]] accordingly (taking into account npband).
If [[bandpp]] is used, set [[nblock_lobpcg]] accordingly (taking into account npband).
[[bandpp]] and [[nblock_lobpcg]] cannot be used simultaneously.
[[bandpp]] should be considered as an internal variable in the future.
[[autoparal]] behavior is not changed, so it sets [[nblock_lobpcg]] according to [[npband]] and [[bandpp]].

A preferable behavior would be to look for the best values of [[npband]] with a fixed value of [[nblock_lobpcg]], and set [[bandpp]] accordingly.
But this would be a big modification of [[autoparal]] and would require more tests.
Ground state parallelization tutorial should be updated too.

See [[test:v9_201]], [[test:v9_202]], [[test:v9_205]].
From L. Baguet (MR957)


**D.3**
Miscellaneous developments for ground state CPU calculations,
including refactoring (e.g. create independent Rayleigh-Ritz and Ortho routines,
originally linked to lobpcg or chebfi modules, and merge of the two RR routines).
Timing is improved thanks to these modifications, also related to GPU development.

From L. Baguet (MR973 and MR981)



**D.4**
Add post comparison of initial [[spgroup]] with final [[spgroup]], with related bug fixes and improvements. Also, information about non-primitive cells.

After echoing all the variables as usual, ABINIT performs a post-analysis of the symmetries and associated [[spgroup]] (and [[spgroupma]]),
with comparison with the initial assessment. In case the post analysis and the initial one give differences,
a comment is issued in the output file, and a more detailed analysis is produced in the log file.
While implementing this feature, several bug fixes and improvements have been done.
Several errors in the (magnetic) space group generator were uncovered.
The documentation of [[spgroupma]] has been fixed.
Unneeded redundant write of [[spgroup]] in the log file is avoided.

The information about the primitive cell when the input contains a non-primitive cell is now echoed in the log file.
See [[test:v9_189]].

By X. Gonze (MR932, MR976)

**D.5**
ANADDB can now read and write DDB files with more than 1000 atoms.

By J. Bouquiaux (MR956)

**D.6**
The new input variable [[eph_ahc_type]] has been introduced to allow computing of adiabatic AHC ZPR (previously, only the non-adiabatic one).
See the new test case [[test:v8_44]], where the new dataset 7 has been added for that purpose.

By S. Ponce (MR962)

**D.7** DMFT susceptibility.

Add local charge and magnetic susceptibility for DMFT in CTQMC
Add two tests, [[test:paral_100]] and [[test:paral_101]] for testing local susceptibility in DMFT.

By F. Gendron (MR963, MR983)

**D.8**
Check on consistency of parallelism input variables [[autoparal]] and [[paral_kgb]] with [[optdriver]].
Introduce [[chkparal]]. Can be disabled by non-zero [[expert_user]] input variable (already existed, from 9.2).
See [[test:v67mbpt_36]].

By X. Gonze (MR982)

**D.9**
New input variable [[vloc_rcut]].
This variable defines the cutoff for the radial mesh used to compute `epsatm`
(the alpha term in the total energy due to the pseudos) and the Bessel transform
for the local part in the case of NC pseudos given in UPF2 format.

This parameter can be used to cut off the numerical noise arising from the large-r tail when integrating V_loc(r) - Z_v/r.
In Quantum Espresso, [[vloc_rcut]] is harcoded to 10 Bohr. However, numerical experiments showed that such value leads to oscillations
in the second order derivatives of the vloc form factors.
For this reason, the default value in Abinit is set to 6.0.

By M. Giantomassi (6 April 2023)

**D.10**
New input variable [[eph_frohl_ntheta]].
Only relevant for [[optdriver]] = 7 and [[eph_task]] = 4, i.e. computation of the electron-phonon self-energy.
This variable defines the angular mesh for the spherical integration of the Frohlich divergence
in the microzone around the Gamma point to accelerate the convergence with the number of q-points.

Numerous tests: [[test:v8_57]], [[test:v8_60]], [[test:v9_60]], [[test:v9_66]] and several tests in the [[tutorial:eph4zpr]]

By M. Giantomassi

**D.11** Progress in the implementation of the strong coupling variational polaron equations.

Introduced a new frohlich_t datatype (replaces the src/78_eph/m_frohlichmodel.f90 module).
The src/78_eph/m_frohlichmodel.f90 module is replaced by the src/78_eph/m_frohlich.f90, containing the frohlich_t datatype.
The previous module provided routines that did not store any data and only produced output to the main output file.
The new datatype is more modular and allows for reusing the Fröhlich model data for multiple calculations.

By V. Vasilchenko (MR971)

**D.12**
Developments in the aTDEP post-processing application.
See e.g. [[test:atdep_38]].

By F. Bottin and J. Bouchet (MR496)

**D.13**
There are significant documentation updates to [[topic:EFG|EFG]] and [[topic:NMR|NMR]] topics, as well as for the [[tutorial:rf2|rf2 tutorial]].

By J. Zwanziger (MR980)


**D.14**
Introduced a warning in m_respfn_driver.f90 regarding non-colinear dfpt in metals for norm-conserving psps. I compared results obtained by dfpt to results obtained via finite differences on Fe bcc. Everything seems alright.
There is no paw implementation yet for [[nspden]]=4.

By M. Mignolet (MR922)

**D.15**
Fixed a parser issue: the input file was not be parsed correctly when more than one environment variable was present in the input file.

By M. Mignolet (MR928)

**D.16**
Fixed d2eig parallel writing.
Reactivated the tests [[test:v6_37]] to [[test:v6_37]],
[[test:v6_50]] to [[test:v6_53]], [[test:v7_50]] to [[test:v6_54]].
Now the communicator is properly passed to the ddb when writing d2eig.

By G. Antonius (MR931)


**D.17**
Fixed a typo in 'bs_nstates' docs where direct diago was mapped to 'bs_algorithm=2' while it is 1.

By F. Goudreault (MR945)

**D.18**
Fix DFT+U + SOC + [[nspden]]=1 case.
Now local magnetic moment is correctly forced to be zero.

By M Torrent (MR948)

**D.19**
Fix in posdoppler routine.

From A. Donkov, through M. Torrent (MR950)

**D.20**
There was a Bug when doing AHC computations with dipoles+quadrupoles activated.
In that case the DDB block dimension is bigger, from $(3*mpert)^2$ to $(3*mpert)^3$
and a reshaping is needed in ddb_get_dielt_zeff.
This has been fixed.
There was also a bug for phonons in supercells, for
the reconstruction of the dynamical matrices using crystal symmetries.
The dimensionless real-space and reciprocal-space primitive translations were used instead of their dimensional
counterpart.

From S. Ponce (MR952 and MR985)

**D.21**
Get rid of all cp added/cp modified lines.

From M. Torrent (MR959)


**D.22**
Introduced new input variable [[invovl_blksliced]].

From M. Torrent

**D.23**
Handle the Debye-Waller when only VB in the active space.
There was a bug in the case where highest
energy state was degenerate and when the interpolation
of the potential was used. In that case the Sternheimer
contribution was omitted.
Other bug fix: missing i eta in the denominator of the Debye-Waller term on the active space.
Such factor is important to exactly cancel the Fan term on the active
space. This missing factor was only present when us-
ing the Fourier interpolation of the perturbed potential.
From S. Ponce (MR972 and commit a84206b4)

**D.24**
Bug fixes for the LRUJ utility, documentation, and subroutines.
Polynomial regression subroutine written from scratch to replace that subroutine.
LRUJ tutorial documentation updated.
LRUJ post-processor now calculates HP error correctly.
Reference files redone and updated.

From L. MacEnulty (MR975)

**D.25** Molecular Berry curvature

DDB support has been added for the molecular Berry curvature. A new block type has been introduced inside anaddb (number 85). It represents the molecular Berry curvature and is of kind d2E.

When computing the molecular Berry curvature (eph_task=14), it is now written to a ddb file (*_BERRY_DBB).

From M. Mignolet (MR984)

**D.26** Bug fixes for Gstore, Molecular Berry curvature and VdW DFT-D3

Gstore module: incorrect parallelization over spins would make abinit crash when computing the GSTORE.nc with nsppol=2.
Molecular Berry curvature module: incorrect initialization of the DDB would result in a crash when using kptopt!=3.
VdW DFT-D3 module: the max/min values of the C6 coefficient would not be displayed correctly if larger than 100.
Partially fix issue #160, removed tests libxc/t34,35 .
Multidataset parsing ignoring ? wildcard when the first index is above 10.

From M. Mignolet (MR999, MR1011, MR1015)

**D.27**

Add new DFT-D3 compatible XC functionals + slight changes in doc. Typo in CPU name corrected.
Suppress -mtune when -march is specified (march=>mtune); this is problematic with clang v<15.
Several improvements of the build system (including better openMP detection).
Autodetect ZDOT bugfix (autotools & cmake) + autodetect some important MPI features (cmake).

From M. Torrent (MR1000, MR1001, MR1008, MR1031, MR1037)

**D.28**

Fix usecase of istwfk==2 when using CHEBFI as diago algorithm .
Add tests checking Kokkos code path.
Fixes for OpenMP GPU (LOBPCG multi-GPU and istwfk==2).
Fix race condition in inwffil (cg_envlop kernel).
Various fixes related to OpenMP GPU code.

From M. Sarraute (MR1002, MR1003, MR1018, MR1020, MR1023, MR1053)

**D.29**

Chebfi2 norm-conserving.
Fix istwfk2.
Fix openblas.

From L. Baguet (MR1021, MR1039, MR1040)


**D.30**

Fix restartxf -1 in Molecular Dynamics. See tests v10 81 & 82.

From F. Brieuc (MR1054)

* * *


## v9.10

Version 9.10, released on June 24, 2023.
List of changes with respect to version 9.8.
<!-- Release notes updated on July 11, 2023. -->

Many thanks to the contributors to the ABINIT project between
September 2022 and April 2023 (with some late contributions until June 2023).
These release notes
are relative to modifications/improvements of ABINIT v9.10 with respect to v9.8.
<!-- Merge requests up to and including MR919. Also, MR923 to MR925 and MR928 are taken into account. -->

The list of contributors includes:
J. Abreu, F. Akhmetov (Radioteddy on github), J.-M. Beuken, A. Blanchet, F. Bruneval, M. Cote, M. Giantomassi, X. Gonze, B. Guster, P. Kesterner,
L. MacEnulty, M. Mignolet, D.D. O'Regan, S. Rostami,
M. Royo, A. Sasani, M. Stengel, M. Torrent, M. Verstraete, A. Zabalo, J. Zwanziger.

It is worthwhile to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases.
This might take some time ...

Xavier

* * *

### **A.** Important remarks and warnings.

**A.1** The names of several tutorial files have been changed to make them easier to understand. "gspw" is now "paral_bandpw", "ucrpa" is now "ucalc_crpa",
"depes" is now "eph_tdep_legacy", "eph" is now "eph_legacy", "ffield" is now "polarization". The tutorial "udet" relying on an old utility
is now superceded by the tutorial "lruj", see section [B.3](#v9.10.B.3). The name of the topic CRPA has been changed to [[topic:CalcUJ]].

**A.2** The default value for [[dosdeltae@anaddb]] has been changed from 1 cm-1 to 0.2 cm-1,
and the default value for [[dossmear@anaddb]] has been changed from 5 cm-1 to 1 cm-1.
Also, the default values for [[dosdeltae@atdep]] has been changed from 4.5d-6 to 0.2 cm-1.
This is to allow default calculations of thermal expansion using abipy to be more stable numerically.

By S. Rostami and X. Gonze (commit 8b7697502c)

* * *

### **B.** Most noticeable achievements

**B.1** Orbital magnetization

The computation of the orbital magnetization and chemical shielding
(in the converse method, that is, with a nuclear dipole moment added)
has been implemented, as described and tested in
[[cite:Zwanziger2023]].  This implementation works for insulators and
metals, with [[nspinor]]=1 and 2. However it works with PAW only, not
with NC pseudopotentials.  Lamb shielding is treated. The PAW atomic
dataset generator "Atompaw" has been updated accordingly to compute
and output the Lamb shielding in xml files.

See [[test:v9_44]], [[test:v9_140]], [[test:v9_141]], [[test:v9_142]], [[test:v9_143]], and [[test:nuc_4]],
with input variables [[orbmag]], [[nucdipmom]], [[lambsig]].

The [[tutorial:nuc| tutorial on properties at nuclei]] has been modified to present such computations.

By J. Zwanziger, M. Torrent and X. Gonze (MR895, 904, 917).


**B.2** Natural optical activity tensor, linear response to a vector potential (orbital magnetic field linear response) and other modifications of the longwave driver.

Computation of the natural optical activity tensor can be performed via the input variable [[lw_natopt]].
This is demonstrated with the test [[test:lw_8]].
The topic [[topic:longwave]] has been upgraded.

In addition, the computation of the linear response to a vector potential in the long-wavelength limit
has been implemented via minimal modifications of the routines that calculate second derivatives
of wavefunctions with respect to wavevector as explained in [[cite:Zabalo2022]].
See test [[test:v9_146]], and input variable [[rf2_dkdk]] with value 2.
Related test [[test:v9_147]].

These are by-products of large modifications of the longwave driver, mimicking the structure of the nonlinear one.
Other improvements related to the large modifications:
i) the number of source code lines has been reduced.
ii) the number of I/O operations has been reduced.
iii) symmetries are now used in order to calculate only linearly independent components of the tensors. This also reduces the number of linear-response functions to precalculate.
iv) the whole structure is now more general, thus facilitating the implementation of future spatial dispersion quantities.

Also, the new input variable [[ffnl_lw]] has been introduced. It allows to reduce memory footprint at the expense of CPU time.

Finally, a bug has been removed (when all KB energies are negative), and a test introduced [[test:v9_145]].

By Miquel Royo, Asier Zabalo and Massimiliano Stengel (MR913).

<a name="v9.10.B.3"></a>
**B.3** Linear response computation of the U and J parameters

The old utility "ujdet" to compute the U and J parameters in DFT+U with the linear response method [[cite:Cococcioni2005]]
has been replaced by the new "lruj" utility. The workflow is different.
The [[tutorial:lruj]] has been written, with three corresponding tests, [[test:lruj_1]], [[test:lruj_2]], [[test:lruj_3]].
See also the input variables [[pawujv]]. The tutorial and corresponding tests "ujdet" have been suppressed.
The tests v5_38, v5_39, v5_40, v6_41 have been suppressed, and replaced by [[test:v9_105]],
[[test:v9_106]], [[test:v9_107]], [[test:v9_108]], [[test:v9_109]]. This is also documented in [[topic:CalcUJ]].

By Lórien MacEnulty with help from David D. O'Regan (MR905, 912).

**B.4** Cumulant method for spectral function

The computation of the electronic spectral function with electron-phonon coupling included is now enabled using the cumulant method.
See [[cite:Nery2018]] and [[cite:Abreu2022]] and other related publications.

Activate it using [[eph_task]]=9 when [[optdriver]] == 7.
New input variable: [[tolcum]].
See the [[test:v9_60]].

By J. Abreu with help from M. Giantomassi (MR 907)


**B.5** Support for norm-conserving pseudopotentials in UPF2 format

Abinit can now read NC pseudos in UPF2 format (both scalar and relativistic version) thanks to the UPF parser imported from quantum espresso and an additional
routine used to convert FR pseudos from (j,kappa) to scalar + SOC term taken from oncvpsp.

The total energy computed with UPF pseudos does not perfectly agree with the one obtained with the corresponding psp8 pseudos.
Most of the difference originates from the value of epsatm $$ \int r^2 (V(r)+\frac{Zv}{r}) dr$$ as the local part in the UPF file is tabulated
on a much larger radial mesh. This should not represent a serious issue as long as total energy calculations are performed with the same pseudopotentials.
Forces, stress tensor and KS eigenvalues are in much better agreement in the systems investigated so far.

See the [[test:v9_130]] and [[test:v9_131]]

By M. Giantomassi (MR896)

**B.6** Initialization of the wavefunctions using atomic pseudo-orbitals

For pseudopotentials that contain the information about atomic local pseudo-orbitals, like the UPF2 format, the wavefunctions
inside ABINIT can be initialized from the Hilbert space spanned by such set of functions,
using the input variable [[wfinit]]=2.

See test [[test:v9_130]].

By M. Giantomassi (MR896)

**B.7** Temperature-dependent XC functionals (free energy), using libXC

The following temperature-dependent XC functionals from libXC are now available:
LDA T-dependent functionals from [[cite:Karasiev2014]], with [[ixc]]=-269, from [[cite:Karasiev2018]],
with [[ixc]]=-318, and from [[cite:Groth2017]], with [[ixc]]=-577.
Previously, the IIT temperature-dependent Free Energy functional of [[cite:Ichimaru1987]], with [[ixc]]=50
had been coded, but not documented. Documentation is delivered in the present release.

See [[test:libxc_22]].

By M. Torrent (MR901)

**B.8** Atomic orbital magnetic moment inside PAW spheres

Implementation of atom-by-atom orbital magnetization integration inside the PAW spheres.
x, y and z components are printed. Decomposition on p, d and f orbitals is also done.
Works only for PAW+U+SOC (and nspden=4).
Works also for orbitals where no U is specified.
New input flag [[prt_lorbmag]]. See test [[test:v9_112]]. Mentioned in [[topic:AtomCentered]].

By A. Sasani & E. Bousquet (MR915).


**B.9** High-temperature DFT: Improvements of the Extended First-Principles Molecular Dynamic (ExtFPMD) calculations.

[[useextfpmd]]=1 now computes contributions using Fermi gas DOS (which was found to be more stable for pure and mixtures).
Old [[useextfpmd]]=1 is now [[useextfpmd]]=4 (tests were changed accordingly)

[[useextfpmd]]=10 controls the Hybrid Thomas-Fermi / Kohn-Sham scheme see [[cite:Hollebon2022]].
This model aims at reducing further the needed number of bands at high temperature compared to Extended FPMD.
No documentation is available since it is work in progress.

Added a new input variable [[extfpmd_nbdbuf]] which specifies the number of bands to use for the buffer if [[useextfpmd]] /= 0.
Among the total number of bands, last [[extfpmd_nbdbuf]] bands
occupation will be set to 0, and ExtFPMD model will take charge of computing
electronic contributions starting from [[nband]] - [[extfpmd_nbdbuf]].
In some cases, setting this input variable to a positive number can solve
convergency problems due to high variations of electron density within the SCF cycle.
Moreover, setting [[extfpmd_nbdbuf]] = [[nband]] should theoretically give
access to Fermi gas orbital free calculations (not tested yet).
This fix has been proposed on the forum by Thomas Gawne (University of Oxford, UK)

By A. Blanchet (MR883 and MR916)


<!--
**B.11** Interface to coupled-cluster CC4S calculations.

The writing of the file needed as input for computations with the CC4S package, <https://manuals.cc4s.org/user-manual>,
allowing to perform coupled-cluster calculations (e.g. CCSD), perturbative triples (and more), can be activated using [[optdriver]]=6 and [[gwr_task]]="CC4S".
See test test:gwr_07 (not yet activated)..

By M. Giantomassi (MR 875, 907)
-->

* * *

### **C.** Changes for the developers (including information about compilers)

**C.1** Improvement of gfortran handling.

Previously Abinit did not add the specific flags for gfortran 12, similar to gfortran 11, because this new version was not anticipated in config/hints/ .
The logic has been changed in config/hints, the default being the newest version.
The old gfortran versions 7, 8 and 9 are treated as specific cases.
So, ABINIT is working with gfortran 12 & 13.

By F Bruneval, with adaptation to gfortran 12 by JM Beuken and to gfortran 13 by M Torrent (MR914, MR919),

**C.2** A new bot for memory profiling has been added to the test farm, named scope_gnu_12.2_mpich .

By JM Beuken (MR888)

* * *

### **D.**  Other changes (or on-going developments, not yet finalized)

**D.1** Low-scaling GW and RPA implementations
The implementation of a cubic scaling real space imaginary time algorithm for GW and RPA is on-going.
By M. Giantomassi (MR875, MR907)

**D.2** Interface to coupled-cluster CC4S calculations.

The implementation of writing a file needed as input for computations with the CC4S package, <https://manuals.cc4s.org/user-manual>,
is working, but not yet in production.
By M. Giantomassi (MR 875, 907)

**D.2** Coulomb interaction with 2D cut-off is now working for the total energy and forces.
This has been tested agains Quantum Espresso implementation. However, stresses are still incorrect.
New test [[test:v9_132]].
By B. Guster, with help from X. Gonze (MR908).

**D.3** Implemented forces and stresses using "gemm" programming model.
This will later allow computing forces/stresses on GPU.
By M. Torrent (MR886 and 887).

**D.4**
Added [[prteig]]=2, in order to print EIG.nc file at each timestep (using already used for other variables TIMx suffix).
(No test - should be added).
By A. Blanchet (MR916).

**D.5**
Allow for band-by-band decomposition of the STM density, using negative values of [[prtstm]].
See test [[test:v4_46]].
By X. Gonze (MR880).

**D.6**
Several improvements for recognition of parallel netcdf for macbookpro and the Zenobe Belgian supercomputer.
By M. Verstraete (MR876).

**D.7**
Improve the initialization of paral_kgb/wfoptalg/istwfk.
By M. Torrent (MR918).

**D.8**
New units are recognized by the input file parser: "meV" (for millielectron-volt);  "S", "Sec" or "Second"; "Kelvin".
By M. Giantomassi (commit 692a4ee0c6) and X. Gonze (commit 39801af30).

**D.9**
New tests of the band parallelism in DFPT: [[test:paral_65]] and [[test:paral_66]].
By M. Giantomassi (commit 31e8aa66d8).

**D.10** Improved developer documentation, section .
[How to add a new test](https://docs.abinit.org/developers/developers_howto/#how-to-add-a-new-test-in-the-test-suite).
By X. Gonze (commit dabc1b905).

**D.11** Fixed typo in CITATION.cff.
By P. Kestener (MR910).

**D.12** New topic [[topic:AtomCentered]] created.
By X. Gonze (commit 425e8c)

**D.13** Improvements of tutorials base3.md and gw1.md in order to better avoid students to commit mistakes.
By X. Gonze (commits 1d56e983f and dfd207458)

**D.14** In m_phgamma.F90, spin-resolved calculations for the case prteliash==3 did not work correctly
since phonon linewidths were calculated for spin=1 only. spin>=2 values are filled by NaNs. This small addition fixes the issue.
Also, a minor format fix.
By F. Akhmetov (Radioteddy on Github). commit bd76768 on abinit github, but directly ported to the trunk/release-9.10 branch.

**D.15** Fix parser problem. The input would not be parsed correctly when more than one environment variable is present in the input file.
By M. Mignolet (MR 928 backported to ABINITv9.10)

**D.16** Improvements of documentation for spinmagntarget and occopt, in the case of ferromagnetic insulators.
Improvements of documentation for rfasr, asr and chneut. Fix timing issues for Fock.
Update doc about the change from npkpt to np_spkpt. Cross refer between U(J) tutorials
By X. Gonze (several commits)

* * *


## v9.8

Version 9.8, released on December 23, 2022.
List of changes with respect to version 9.6.
<!-- Release notes updated on March 15, 2023. -->

Many thanks to the contributors to the ABINIT project between
October 2021 and August 2022, and some late contributions up to April 2023 ! These release notes
are relative to modifications/improvements of ABINIT v9.8 with respect to v9.6.
<!-- Merge requests up to and including MR874. Also, MR881, 882, 885, 891, 892, 894, 897, 898, 899, 900, 902, 903, 911, are taken into account. -->

The list of contributors includes:
B. Amadon, G. Antonius, L. Baguet, S. Bandyopadhyay, L. Bastogne, J.-M. Beuken, J. Bieder, A. Blanchet,
F. Bottin, J. Bouchet, E. Bousquet, F. Brieuc, V. Brousseau-Couture, N. Brouwer, F. Bruneval, M. Cote,
C. Espejo, Ph. Ghosez, M. Giantomassi, O. Gingras, X. Gonze, B. Guster, P. Kesterner,
R. Outerovich, Ch. Paillard, M. Royo, A. Sasani, B. Sataric, M. Schmitt, F. Soubiran,
M. Torrent, M. Verstraete, He Xu, J. Zwanziger.

It is worth to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases.
This might take some time ...

Xavier

* * *

### **A.** Important remarks and warnings. Also, hotfixes for v9.8.3 (A.4 to A.10).

**A.1** Warning: the input variables prtefg and prtfc have been renamed [[nucefg]] and [[nucfc]].

By J. Zwanziger (MR850)

**A.2**
Within PAW in ABINIT, the ZORA relativistic factor for all PAW+spin-orbit coupling has been corrected:
A ^2 was missing for the 1/(1-(v/c)^2) factor.
This changes slightly all the PAW+SOC results (because the factor is very small), without any influence on physical results.
However, previous results cannot be anymore obtained.

By M. Torrent (MR849)

**A.3**
The default value for [[rfdir]] is now (1 1 1), instead of (0 0 0).
The one of [[rfatpol]] is now (1 [[natom]]), instead of (1 1).

By X. Gonze (MR852)

**A.4** The default values of dossmear (@anaddb) and dosdeltae (@anaddb and @atdep) have changed, to smaller values, to improve default accuracy.

By X. Gonze

**A.5**
Correction and cleaning of DFT+U with magnetism ([[nspden]]=4), following the merge 881.
Introduce [[optdcmagpawu]] input to control different choices of DC term (for tests and code comparisons, not useful for production).
Versions before 9.8 is equivalent to [[optdcmagpawu]]=1 (no magnetism in the DC term). Now the default is 3 (magnetism in the DC term). Some refs are changed accordingly.
Change of the tests v9/76,77 and 78 for better precision. Add the tests v9/88 and 89 to tests the combination of PAW+U with different [[pawxcdev]].

By L. Baguet (MR898)

**A.6**
Some users have reported that with the new [[optcell]]=4,5 or 6 option, the calculation
continues when the criterion should be already reached. For example, with [[optcell]]=5 (relaxing b only), the
xx component of the stress could still be large but shouldn't be checked. This has been fixed.

By Xu He (MR897)

**A.7**
More fix on the generation of symmetry-adapted terms in Multibinit. A more complete set of terms are now generated, with improved time and memory performance.

By Xu He and A. Sasani (MR899)

**A.8**
m4 detection netcdf(-f) //: AC_RUN_IFELSE hangs with some compilers then switch to AC_LINK_IFELSE

By JM Beuken (MR900)

**A.9**
Several improvements for external libs: FFTW3 thread safety, openBLAS multithreading, netCDF Fortran parallel

By M. Torrent, with patch from P. Kestener (MR902)

**A.10**
Several fixes, including documentation for icutcoul and related input variables

By. X. Gonze (MR903)

* * *

### **B.** Most noticeable achievements

**B.1** Lattice Wannier functions can be computed, using the SCDM-k algorithm, or the projected Wannier function algorithm.

The [[tutorial:lattice_wannier]] has been written to teach how to construct such lattice Wannier functions.
They are used in localized bases for atomic distortions.
One typical use case is to build an effective Hamiltonian of collective, localized, atomic displacements (see the next achievement).

A script (compare_phbands.py) to compare phonon/LWF band structures is added to the scripts/post_processing directory.

Related input variables:
[[lwfflag@anaddb]],
[[lwf_anchor_ibands@anaddb]],
[[lwf_anchor_proj@anaddb]],
[[lwf_anchor_qpt@anaddb]],
[[lwf_disentangle@anaddb]],
[[lwf_mu@anaddb]],
[[lwf_ngqpt@anaddb]],
[[lwf_nwann@anaddb]],
[[lwf_projector@anaddb]],
[[lwf_sigma@anaddb]].

Related topic [[topic:LWFModel]].
See tests [[test:lwf_1]], [[test:v9_110]], and [[test:v9_111]].

This feature is still under heavy development. The current version should be regarded as a "technology preview".

By He Xu (MR844)

**B.2** Lattice Wannier function dynamics is available inside the second-principle engine MULTIBINIT.

Once Lattice Wannier functions are available, they can be used with different dynamical algorithm to deduce heat capacity, susceptibility, structural phase transition, critical temperature, instead of doing standard molecular dynamics with all atoms.
Related input variables:
[[lwf_constraint@multibinit]],
[[lwf_dt@multibinit]],
[[lwf_dynamics@multibinit]],
[[lwf_init_state@multibinit]],
[[lwf_init_hist_fname@multibinit]],
lwf_mc_avg_amp@multibinit,
[[lwf_nctime@multibinit]],
[[lwf_ntime@multibinit]],
[[lwf_pot_fname@multibinit]],
[[lwf_taut@multibinit]],
lwf_temperature@multibinit,
[[lwf_temperature_end@multibinit]],
[[lwf_temperature_nstep@multibinit]],
[[lwf_temperature_start@multibinit]],
[[lwf_var_temperature@multibinit]].

See the tutorial [[tutorial:lwf_model]] and related tests.

This feature is still under heavy development. The current version should be regarded as a "technology preview".

By He Xu (MR851)

**B.3** Numerous miscellaneous improvements of the second-principle engine MULTIBINIT.

Hybrid Monte-Carlo Mover with NPT ensemble ([[ionmov]]=25).

MULTIBINIT can now be used without 'files' file, as the main ABINIT or ANADDB.
Related input variables [[latt_pot_fname@multibinit]], [[latt_harm_pot_fname@multibinit]],
[[latt_anharm_pot_fname@multibinit]], [[latt_training_set_fname@multibinit]],
[[latt_test_set_fname@multibinit]], [[spin_pot_fname@multibinit]], [[spin_init_hist_fname@multibinit]], [[slc_pot_fname@multibinit]],
[[outdata_prefix@multibinit]].

The MULTIBINIT tutorial has been improved. In particular, it is now starting with a global introduction ([[tutorial:multibinit]]).

The following MULTIBINIT features and input variables have been introduced:
define oblique supercells ([[ncellmat@multibinit]]), print Goal-Function values in CSV format ([[prt_GF_csv@multibinit]]),
impose the seed for MULTIBINIT spin/LWF dynamics to obtain reproducible results ([[randomseed@multibinit]]),
specify three weights for Energy, Forces and Stresses in the calculation of the Goal Function ([[bound_factors@multibinit]]),
specify three weights for Energy, Forces and Stresses in the calculation of the Goal Function during the fit process ([[fit_factors@multibinit]]),
specify three weights for Energy, Forces and Stresses in the calculation of the Goal Function during the optimization process ([[opt_factors@multibinit]]),
specify the relative penalty for the determination of bounding coefficient values ([[bound_penalty@multibinit]]),
activate the generation of pure displacement coefficients ([[fit_dispterms@multibinit]]),
specify the number of anharmonic coefficients per symmetric irreducible atom to add during fit process ([[fit_ncoeff_per_iatom@multibinit]]),
specify the number of coefficients imposed with fixed value as in the input xml during the fit process for the mode
([[fit_nimposecoeff@multibinit]]),
specify the indices of the imposed coefficients with fixed coefficient value during the fit process for the model ([[fit_imposecoeff@multibinit]]),

See the tests in which these input variables are used.

By He Xu, M. Schmitt, A. Sasani, L. Bastogne, S. Bandyopadhyay, and P. Ghosez (MR812, 851, 868, 894)


**B.4** The TDEP formalism implemented in ABINIT (aTDEP), allowing to compute temperature-dependent phonon band structures,
 has been made more robust, and tested extensively. One tutorial ([[tutorial:atdep1]]) is now available.
See the thirty-seven tests from [[test:atdep_01]] to [[test:atdep_37]] and the tests mentioned in the
tutorial [[tutorial:atdep1]]. Related publication, see [[cite:Bottin2020]].

By F. Bottin, J. Bieder and J. Bouchet (MR 836).


**B.5**
The conducti utility can treat Spin-Orbit Coupling for transport properties (within PAW): conductivity, XANES, transport coefficients, ...
See [[cite:Brouwer2021]].
The writing of documentation and tutorial is in progress.

By M. Torrent and N. Brouwer (MR849).

* * *

### **C.** Changes for the developers (including information about compilers)

**C.1** Added support for NAG 7.1

From J.-M. Beuken (MR830).

**C.2** Update build system to allow the use of NVTX library, providing profiling annotations (only when gpu is enabled).
This makes more readable profiling and tracing information when viewed with nsys-ui.
If abinit is built without gpu, annotations completely vanish at compile time.

From P. Kestener (MR843)

**C.3** Improve detection of inlined macros. The developers should now use

    ABI_SFREE(allocatable_array)

instead of

    if(allocated(allocatable_array)) ABI_FREE(allocatable_array)

From M. Giantomassi (MR859)

**C.4** The DDB IO routines have been refactored throughout the different main codes
(respfn, gstate, anaddb, nonlineal, longwave, ddb_interpolate, gruns_new, thmeig).
The ddb merging routines have been refactored (merge_ddb replaces mblktyp1 and mblktyp5).
New MPI treatment of ddb reading (no more free-for-all reading).
Removed several arguments of ddb_from_file: natom, natifc and atifc.
Moved routine dfptnl_doutput from m_ddb to m_nonlinear. Remove dfpt_lw_doutput .

From G. Antonius (MR872)

**C.5** Suppressed CHILDREN and PARENTS sections in ROBODOC header.
They were not maintained automatically anymore, so many had become misleading.

From X. Gonze (MR874)

**C.6**
A new type of file 'GSTORE' has been introduced,  for the electron-phonon matrix element storage, where it is taken advantage
of filters (on energy or band or wavevectors) to reduce the size compared to the more usual EIG1 file type.
It contains also the related metadata.
Related input variables: [[getgstore_filepath]], [[gstore_cplex]], [[gstore_with_vk]], [[gstore_kzone]], [[gstore_qzone]],
[[gstore_kfilter]], [[gstore_brange]], [[gstore_erange]].

By M. Giantomassi (MR870)


* * *

### **D.**  Other changes (or on-going developments, not yet finalized)

**D.1** New tutorial on the use of the Z2pack postprocessor ([[tutorial:z2pack]]).
In this tutorial, the topological transition from Z2 trivial to Z2 non-trivial under pressure is reproduced.
A test has been added but one needs to install z2pack to test it. See doc/tests/tutoplugs/z2.py, and [[tutorial:z2pack]].

From O. Gingras, V. Brousseau-Couture and M. Cote (MR 871)

**D.2** Improved the Wannier90 tutorial ([[tutorial:wannier90]]).
There is a new part, that describes projecting the band structure of La$_2$CuO$_4$ to a single orbital
($d_{x^2-y^2}$) Hamiltonian.
This is the first step in order the study the Mott transition in La$_2$CuO$_4$ using DFT+DMFT.
From this step, one can use TRIQS to study the Mott transition.

From O. Gingras (MR 871)

**D.3** Make abinit compile and run with libxc v6

From M. Torrent (MR885)

**D.4** Modified [[optcell]]=4,5,6 to allow for relaxation of a vectors length and angle without constraining it to be orthogonal to the three others

From C. Paillard (MR858)

**D.5** New input variable for Optic, [[nband_sum@optic]], allowing to select the maximal number of bands included. This is convenient for convergence studies.

From M. Giantomassi (MR817)

**D.6** Implementation of Chebyshev filtering algorithm, version 2.

From B. Sataric, J. Bieder, and M. Torrent (MR 826)

**D.7** Miscellaneous improvements of the iterative Boltzmann transport equation coding. Input variable [[ibte_alpha_mix]].

From M. Giantomassi (MR821)

**D.8** Improvements related to electron-phonon interaction.

Miscellaneous improvement of the electron-phonon part of ABINIT (documentation, bug fixes, improved parallelism)

Work in progress: a new tutorial [[tutorial:eph4isotc]]
to demonstrate the computation of superconducting properties within the isotropic Eliashberg formalism.
See tests in the [[tutorial:eph4isotc]].

By M. Giantomassi (MR870)

**D.9** Improvements of the cRPA determination of the U and J parameters (default keywords, tests and tutos).

From R. Outerovich (MR835)

**D.10** Work on Real-time Time-Dependent Density Functional Theory implementation within ABINIT.

From F. Brieuc (MR853)

**D.11** Test of the Zero-point renormalization for hexagonal systems (Frohlich generalized model). See [[test:v8_60]].

From B. Guster (MR815)

**D.12** Implement [[chksymtnons]]=3: FFT grid and spatial symmetry operations are coherent.

From X. Gonze (MR828)

**D.13** Improvement of VdW-DF.
Corrections to the implementation of vdW-DF non-local functional found in files 56_xc_/m_xc_vdw.F90 and 98_main/vdw_kernelgen have been done, also several debugging lines have been included. From a careful comparison of the computed quantities from a pre-calculated density and gradient obtained with Siesta, it has been observed that Abinit is getting the correct numbers up to the final 3D FFTs which still present differences that lead to incorrect values for the vdW-DF correlation energy.

From C. Espejo (MR829)

**D.14**
Make use of gvnlxc optional in getghc, which allows memory savings in lobpcg2 and chebfi2. Could be useful in other parts of the code.

From L. Baguet (MR831)

**D.15**
Miscellaneous changes in DMFT: keyword for Wannier orthonormalisation,
implementation of calculation of the weight of configuration in CTMQC,
double counting for charge only DFT+DMFT ([[dmft_dc]]=6), work in progress concerning alternate calculation of electronic entropy in DMFT.

From B. Amadon and R. Outerovich (MR833)

**D.16**
Move DDK reading outside of loop for non var matrix element calculations. Should be much more efficient IO.

From M. Verstraete (MR840)

**D.17**
Use wfdgw_t subclass in GW/BSE code. This is needed so that the wfd in EPH does not allocate bks_tab whose size scales badly with nprocs and nkibz.

From M. Giantomassi (MR841)

**D.18**
Add spinat to GSR.nc

From M. Giantomassi (MR842)

**D.19** Orbital magnetism progress

From J. Zwanziger (MR847)

**D.20** Fix [[prtwf]] and [[prtpot]] in DFPT

From M. Giantomassi (MR848).

**D.21** Improved ABINIT+TRIQS python invocation interface.

From O. Gingras (MR851)

**D.22** Bug fixes for PAW+Hybrid

From F. Bruneval and F. Soubiran (MR854)

**D.23** Improvements of the Frohlich model implementation (e.g. dielectric average decomposition).

From B. Guster (MR860)

**D.24** Bug fix of DFT+U in the non-collinear case

From E. Bousquet (MR881)

**D.25** GPU coding: inverse overlap matrix, and non-local operator

From P. Kestener (MR843 and 869)

**D.26** Restructuring of the tutorial index page, doci/tutorial/index.md ..

From X. Gonze

**D.27** A CITATION.cff file has been created. Also, a LICENCE file (pointing toward COPYING).

From X. Gonze

**D.28** Updated the tarball for the fallbacks

From JM Beuken (MR911)

**D.29** Improve the procedure to examine the convergence with respect to ecut and pawecutdg in tutorial PAW1.

From X. Gonze (20230424)

**D.30** Miscellaneous additional bug fixes, typos fixes, or upgrade of build system.

By F. Goudreault (MR816), M. Giantomassi (MR821 and 845), P. Kestener (MR827 and 843),
A. Blanchet (MR832), C. Paillard (MR834), M. Verstraete (MR837),
M. Torrent (MR838 and 873), B. Seddon and X. Gonze (MR839 and 855), L. Baguet (MR857),
J.-M. Beuken (MR882).

* * *

## v9.6

Version 9.6, released on October 4, 2021.
List of changes with respect to version 9.4.
<!-- Release notes updated on November 9, 2021. -->

Many thanks to the contributors to the ABINIT project between
February 2021 and September 2021. These release notes
are relative to modifications/improvements of ABINIT v9.6 with respect to v9.4.
<!-- Merge requests up to and including MR814 except MR812, then also MR818, 819, 820 and 822 are taken into account. -->

The list of contributors includes:
L. Baguet, J.-M. Beuken, J. Bieder, A. Blanchet,
J. Clerouin, C. Espejo, M. Giantomassi, O. Gingras, X. Gonze, F. Goudreault,
B. Guster, Ch. Paillard,
Y. Pouillon, M. Rodriguez-Mayorga, M. Royo, F. Soubiran,
M. Torrent, M. Verstraete, J. Zwanziger.

It is worth to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases.
This might take some time ...

Xavier

* * *

### **A.** Important remarks and warnings.

(nothing to mention for this v9.6)

* * *

### **B.** Most noticeable achievements

**B.1** Band-parallel implementation of DFPT: the memory footprint is now distributed over different processors.
Previously, the memory was distributed only for k-point parallelism.
This is automatically managed, no user action is to be taken to activate this memory saving.

See test [[test:dfpt_04]].

By M. Verstraete (MR784, 803)

**B.2** The Iterative Boltzmann Transport Equation (IBTE) to compute the electric conductivity has been implemented.
To activate the IBTE, use [[ibte_prep]] = 1 with [[eph_task]] -4.
The IBTE solver can also be invoked in standalone mode by providing a SIGEPH file with [[eph_task]] = 8.
Related input variables: [[ibte_niter]], [[ibte_abs_tol]] and [[ibte_alpha_mix]].
See test [[test:v9_65]].

By M. Giantomassi (MR794)

**B.3** The computation of dynamical quadrupoles and flexoelectricity is now available within the GGA.
Test for GGA + longwaves [[test:v9_46]].

Also, the usage of the quadrupoles has been rationalized (and made easier) in anaddb
as the default value of [[dipquad@anaddb]] and [[quadquad@anaddb]] has been changed to 1.
This means that dipole-quadrupole and quadrupole-quadrupole contributions are always included
in the Fourier interpolation of the dynamical matrix provided the DDB provides these terms.

This "default behaviour" is similar to the one used for the dipole-dipole treatment.
Indeed, the default value of [[dipdip@anaddb]] is 1 hence the dipolar term is automatically included
if the DDB contains the Born effective charges and the electronic dielectric tensor.
Still, the user can deactivate the inclusion of the different terms by setting the corresponding
variable to zero for testing purposes.
See the ANADDB input variables and test [[test:lw_6]]

By M. Royo with contribution by M. Giantomassi for the change of default (MR795)


**B.4** Stresses are available within cDFT (constrained DFT).
See tests [[test:v9_01]], [[test:v9_02]] and [[test:v9_03]].

By X. Gonze (MR802)


**B.5** The computation of effective mass renormalization due to electron-phonon coupling, treated in the generalized Frohlich model,
is now available, for cubic materials. An article has been submitted, see <https:arxiv.org/abs/2109.12594>.
Activate it using [[eph_task]]=10.

See test [[test:v9_66]].

By B. Guster (MR800)

**B.6** Important speed-up of the PAW calculations is allowed thanks to the storage of "cprj" coefficients.
See the input variable [[cprj_update_lvl]]. However, at present this is only possible for ground-state
calculations, with several restrictions, spelled in [[cprj_update_lvl]]. So, this is not activated by default.
There is also an internal variable [[cprj_in_memory]] exposed in the documentation.
Other input variables have been introduced in the development process: [[fft_count]] and [[nonlop_ylm_count]].
They allow one to monitor better the number of FFTs and non-local operator applications.

See tests [[test:v9_71]], [[test:v9_72]], [[test:v9_73]] and [[test:v9_74]].

By L. Baguet (MR793).

**B.7** The Extended First-Principles Molecular Dynamics has been implemented.
This method allows one to drastically reduce the needed number of bands for high temperature simulations,
using pure single plane waves description based on the Fermi gas model beyond explicitly computed bands.
The implementation and usage are described in <https://doi.org/10.1016/j.cpc.2021.108215> (Authors: *A. Blanchet, J. Clérouin, M. Torrent, F. Soubiran*).

See [[topic:ExtFPMD]], as well as the input variables [[useextfpmd]] and [[extfpmd_nbcut]],
and test [[test:v9_92]].

By A. Blanchet, J. Clérouin, M. Torrent, F. Soubiran. (MR788).


* * *

### **C.** Changes for the developers (including information about compilers)

**C.1** Supported compilers

* gfort (GNU) compiler: v11 newly supported.
* ifort (Intel) compiler: v21.4 newly supported.
Two new bots introduced in the test farm: alps_intel_21.4_elpa and graphene_gnu_11.2_macports .

By JM Beuken

* * *

### **D.**  Other changes (or on-going developments, not yet finalized)

**D.1** New input variable for "optic": [[prtlincompmatrixelements@optic]].

Added this flag in order to make it possible to print the different elements that are used to build the susceptibility of the linear component of the dielectric tensor.
These elements are namely: the matrix elements, the renormalized electronic eigenvalues, the occupations and the kpt weights.
Everything is dumped into the _OPTIC.nc file by the main process. Thus optimization could be done memory-wise and speed wise if MPI-IO is implemented for this nc file.

[[test:v9_49]] was created which is the same as [[test:v9_48]] except with the aforementioned flag set to 1.
This test checks that everything works well even though we print the matrix elements.
It does not test that matrix elements are well printed because that would require testing of the OPTIC.nc file.
Although it is possible to check that it works well using a simple python script (see the figure in the merge request on Gitlab).
(Note that the tests v9_13 and v9_14 have been moved to v9_47 and v9_48 in this change).

By F. Goudreault (MR776)

**D.2** New radial sine transform for the vdW-DF kernel.

By C. Espejo (MR 797)


**D.3** New test of orbital magnetism, [[test:v9_37]].
Also, on-going work on orbital magnetism, including use with DDK wavefunctions.

By J. Zwanziger (MR767, MR775, MR779 and MR787)

**D.4** Migration to mkdocs==1.1.2 and  mkdocs-material==7.0.6.
mksite.py now requires python >= 3.6 .
Activated search capabilities, available in the new mkdocs version.

By M. Giantomassi (MR774)

**D.5** Fixed bug in make_efg_onsite for [[nspden]]=2 case.

By J. Zwanziger (MR783)

**D.6** Correction of tutorials Rf1 and Rf2 for version 9

By O. Gingras (MR785)

**D.7** Fixed errors and bugs detected by using -ftrapuv intel option

By M. Giantomassi (MR789)

**D.8** Bug fix in [[nspden]]=4 DFPT for Fe

By M. Verstraete (MR790)

**D.9** Fixed a spurious test line 711 of m_occ.F90 that caused abinit to abort in the case of [[occopt]]=9,
if the number of conduction bands was enough to accommodate nqFD but not enough to accommodate nelect.

By Ch. Paillard (MR791)

**D.10** GW methodology with Kohn-Sham density matrix.
Solving a bug producing a segmentation fault when using [[bdgw]] and [[gw1rdm]].
New test [[test:v9_37]].

By M. Rodriguez-Mayorga (MR792)

**D.11** Introduced new input variable use_oldchi.
This input variable is temporary, for testing purposes. It is documented, but not tested.

By Wei Chen (modified line 743 in src/95_drive/screening.F90 on 23 April 2021).

**D.12** The input variable [[rfstrs_ref]] has been introduced, but not yet documented and tested, as this is on-going work.

By M. Royo


**D.13** Miscellaneous additional bug fixes, or upgrade of build system.
in the upgrade of tutorials).
By J. Bieder, M. Giantomassi, Y. Pouillon, M. Torrent, J. Zwanziger.

* * *

## v9.4

Version 9.4, released on February 25, 2021.
List of changes with respect to version 9.2.
<!-- Release notes updated on April 30, 2021. -->

Many thanks to the contributors to the ABINIT project between
November 2020 and April 2021. These release notes
are relative to modifications/improvements of ABINIT v9.4 with respect to v9.2.
<!-- Merge requests up to and including MR766 are taken into account, also MR768 (backported) up to MR772 and MR780, 781, 782. -->

The list of contributors includes:
B. Amadon, L. Baguet, J.-M. Beuken, J. Bieder, E. Bousquet, V. Brousseau, F. Bruneval,
W. Chen, M. Cote, M. Giantomassi, O. Gingras, X. Gonze, F. Goudreault,
B. Guster, T. Karatsu, A. H. Larsen, O. Nadeau, R. Outerovich, Ch. Paillard, G. Petretto,
S. Ponce, Y. Pouillon, G.-M. Rignanese, M. Rodriguez-Mayorga, M. Schmitt,
M. Torrent, M. Verstraete, He Xu, J. Zwanziger.

It is worth to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases.
This might take some time ...

Xavier

* * *

### **A.** Important remarks and warnings.

**A.1** The [[charge]] variable is obsolete, and has been replaced by [[cellcharge]]. Indeed, [[charge]] was quite ambiguous,
and the string `charge` present in some other input variables. For the time being, ABINIT still recognizes [[charge]]
in the input file, but this might not last longer than in ABINITv9.

**A.2** There is a new check, governed by the input variable [[chksymtnons]], to examine
whether the [[tnons]] of all symmetry operations
is zero or a rational number with small denominator, which is required for GW calculations as implemented in ABINIT.
It is always possible to choose
the origin of coordinate with such properties, and ABINIT can gives suggestions. If you do not want to change
your coordinate origin (e.g. you have no intention to perform a GW calculation), set [[chksymtnons]] to zero,
or (but this is more dangerous) set the meta-variable [[expert_user]] to one to disable several checks of input variables at once.

By X. Gonze (MR712)

**A.3** The input variable npkpt has been changed to np_spkpt . Indeed the parallelism governed by npkpt was about spin and k points,
not only k points. For the time being npkpt is still admitted, but will become obsolete at the next major version change.

By X. Gonze

**A.4** When [[nimage]]>1, the default value of [[prtgsr]] is now 0, like for several prt* variables.

By X. Gonze

**A.5** The code does not stop anymore at the first occurence of overlap between PAW spheres being larger than [[pawovlp]]
in case of [[ionmov]]/=0 or [[imgmov]]/=0, but only at the second occurrence per dataset. Indeed, such trespassing might only be transient.
See the description of [[pawovlp]].

By X. Gonze

* * *

### **B.** Most noticeable achievements

**B.1** The RMM-DIIS algorithm has been implemented.
This SCF (ground-state) algorithm is faster, but potentially more unstable,
than the CG or LOBPCG algorithms, for medium to large size systems,
as it has less cubic scaling steps (e.g. orthogonalisation). Typically used for molecular dynamics
or structural relaxations as the restart from the previous time step gives RMM-DIIS less opportunities to fail.
See input variables [[rmm_diis]] and [[rmm_diis_savemem]].
Several tests exist ([[test:paral_32]], [[test:paral_63]], [[test:paral_64]], [[test:v9_29]], [[test:v9_30]])
covering many cases, including NC/PAW and spin-orbit.
Note that the PAW version of RMM-DIIS is more unstable than the NC one and extra operations are needed to make it convergence.
So the speedup for PAW calculations is not as good as the one observed for NC.

By M. Giantomassi (MR757, MR719, MR718)

**B.2** The treatment of quasi-Fermi energies in the valence and conduction band,
with populations of electrons (in the conduction bands) and holes (in the valence bands)
has been implemented (gapped materials only, of course).
This has been used e.g. in [[cite:Paillard2019]].
See the variables [[occopt]]=9, and [[nqfd]].
See also the related input variable: [[ivalence]].
Internal variables ne_qFD and nh_qFD are presently initialized to [[nqfd]], which is NOT INTERNAL.
See test [[test:v9_91]].

By Ch. Paillard (MR755).

**B.3** All the tutorials have been carefully reexamined and improved, when appropriate.
In particular:

- the old pseudopotentials have been replaced by new ones from e.g. [pseudodojo](http://www.pseudo-dojo.org/) or [JTH](https://www.abinit.org/psp-tables);
- the text of the tutorial uses the new convention for launching abinit (e.g. `abinit input_file` instead of `abinit < files_file`) and the pseudopotentials are mentioned in the input file;
- the input file suffix has been changed from `.in` to `.abi`and the output file suffix has been changed from `.out`to `.abo`;
- the input files have been cleaned when adequate, and many have been restructured using a template;
with populations of electrons (in the conduction bands) and holes (in the valence bands)

<!--
Also, a new tutorial, [[tutorial:eph4isotc]], is available (with tests [[test:eph4isotc_1]] to [[test:eph4isotc_4]],
-->
Also, [[tutorial:nlo]] and [[tutorial:eph4zpr]] have been enlarged to new developments
(see [[test:eph4zpr_8]] and [[test:nlo_6]]).

By B. Amadon, L. Baguet, J. Bieder, E. Bousquet, F. Bruneval,
W. Chen, M. Cote, O. Gingras, M. Giantomassi, X. Gonze, F. Goudreault,
B. Guster, O. Nadeau, R. Outerovich, S. Ponce, M. Schmitt,
M. Torrent, M. Verstraete, He Xu, J. Zwanziger (numerous MRs)

**B.4** The GW 1-body reduced density matrix (1RDM) from the linearized Dyson equation has been implemented.
Its effect on the Hartree-Fock expectation values and therefore on the GW quasiparticle energies can be evaluated.
The resulting total energy parts, kinetic energy (including correlation), electron-nucleus, Hartree, Exchange, can be calculated.
Together with the Galitskii-Migdal correlation, it gives a new approximation the self-consistent GW total energy.
See input variables [[gw1rdm]], [[x1rdm]], also [[irdchkprdm]], [[prtchkprdm]] and [[gwgmcorr]].
See tests [[test:v9_33]] to [[test:v9_36]].

Also, some missing tests have been added:
- GW calculations based on Hartree-Fock wavefunctions can use mini Brillouin Zone integration technique, see [[test:v9_31]].
- A new test [[test:v9_40]] has been provided for the computation of the susceptibility matrix
$\chi_0$ with [[inclvkb]].

By Mauricio Rodriguez-Mayorga and F. Bruneval (MR722).

**B.5** The pSIC (polaron self-interaction corrected) methodology has been implemented.
See [[cite:Sadigh2015]] and [[cite:Sadigh2015a]]. This is based on the `images` capability
of ABINIT, that has been extended to different values of the input variable [[cellcharge]]
for different images, and also parallelized. To activate pSIC, use [[imgmov]]=6 with the proper occupation numbers.
See the test examples [[test:v9_22]], [[test:psic_01]], [[test:psic_02]] and [[test:psic_03]].

By X. Gonze (initial test from C. Tantardini) (MR770).

**B.6** The computation of the electric conductivity has been implemented for metals
in the relaxation-time approximation with transport lifetimes computed from the imaginary part of the Fan-Migdal self-energy.
See tests [[test:v9_62]] to [[test:v9_65]].

By O. Nadeau (MR756, MR716)

**B.7** Implementation of the  i-pi client-server protocol as described in [[cite:Kapil2019]].
This option requires [[ionmov]] 28 and the specification of the socket via command line options.
For UNIX socket, use: --ipi {unixsocket}:UNIX .
For INET socket, use  --ipi {host}:{port} .
Usage example:

     abinit run.abi --ipi {unixsocket}:UNIX > run.log

Note that, at present, this feature is mainly used to interface ABINIT
with the ASE optimization routines. Moreover the user is responsible for creating an input
file with tuned tolerances to prevent Abinit from exiting when internal convergence is reached.
See examples available in the [ASE documentation](https://wiki.fysik.dtu.dk/ase/dev/ase/calculators/socketio/socketio.html)

By M. Giantomassi and A. H. Larsen.


* * *

### **C.** Changes for the developers (including information about compilers)

**C.2** Test farm: new and obsolete bots

* Bots introduced in the test farm: scope_gnu_10.2_paral
* Bots removed: cronos2 (replaced by scope_gnu_10.2_paral)

Bigdft tests have been activated on ALPS.

By JM Beuken

**C.3** Supported compilers

* gfort (GNU) compiler: v10 newly supported, v5 obsolete
* ifort (INTEL) compiler: v15 obsolete (but to be reintroduced for licence reasons)

Support for AOCC has been added in the build system.

By JM Beuken

* * *

<a name="v9.4.D.1"></a>
### **D.**  Other changes (or on-going developments, not yet finalized)

**D.1** Calculation of Luttinger parameters (in the Frohlich model) and echo (see [[test:v8_57]]).
By V. Brousseau (MR736)

**D.2** New test [[test:v9_90]] of the treatment of the Coulomb interaction for low dimensional materials (0D and 2D).
By B. Guster.

**D.3** DFPT (including ddk perturbation) can now be done in the presence of [[nucdipmom]].
By J. Zwanziger (MR749)

**D.4** Bug fix and new test [[test:v9_43]] for the use of ANADDB in the presence of a large [[tolsym]] value.
By G. Petretto

**D.5** Several bug fixes related to the treatment of inaccurate atomic positions (and large tolsym).
Several test have been created test:v9_17 to test:v9_20 (NOTE: all these tests are now v9_180 to v9_199).

**D.6** AiiDA+ABINIT developments

  - A AiiDA plugin for Abinit has been developed: <https://github.com/sponce24/aiida-abinit>, also indexed at <https://aiidateam.github.io/aiida-registry/>.
  - AiiDA is now supporting psp8 type pseudopotentials, <https://github.com/aiidateam/aiida-pseudo>. This also implied modifications of the
    [pseudodojo](http://www.pseudo-dojo.org), that now includes .djrepo for each type of pseudopotential.
  - Work on a [common relaxation workflow](https://github.com/aiidateam/aiida-common-workflows) for a dozen of codes, including ABINIT.

By S. Ponce, also with G.-M. Rignanese, G. Petretto, M. Giantomassi.

**D.7** The new input variable dmft_wanorthnorm has been introduced, see [[test:v6_07]] and [[test:v6_46]]. However, it should still be documented.
By B. Amadon.

**D.8** Increase stack size limit inside xmpi_init using POSIX C-API
By M. Giantomassi. MR 770.

**D.9** Document i-pi interface with links to ASE docs.
By M. Giantomassi. MR 770.

**D.10** Correction (from a message on the forum) related to the forces in electron-positron mode.
By M. Torrent. MR 780

**D.11** Miscellaneous additional bug fixes, improvements of documentation including for the build system (many other were made
in the upgrade of tutorials)..
By B. Amadon, L. Baguet, F. Bruneval, T. Karatsu, G. Petretto, Y. Pouillon, M. Torrent, J. Zwanziger.

* * *

## v9.2

Version 9.2, released on September 30, 2020.
List of changes with respect to version 8.10.
Release notes updated on November 10, 2020.

Many thanks to the contributors to the ABINIT project between
October 2018 and November 2020. These release notes
are relative to modifications/improvements of ABINIT v9.2 with respect to v8.10.
Merge requests up to and including MR692 are taken into account, also MR 697-702, 705, 707-710, 712, 715.

The list of contributors includes:
B. Amadon, G. Antonius, L. Baguet, J.-M. Beuken, J. Bieder, J. Bouchet, E. Bousquet, F. Bruneval, G. Brunin, Wei Chen,
J.-B. Charraud, Ph. Ghosez, M. Giantomassi, O. Gingras, X. Gonze, F. Goudreault,
B. Guster, G. Hautier, Xu He, N. Helbig, F. Jollet,
H. Miranda, F. Naccarato, R. Outerovitch, G. Petretto, N. Pike, Y. Pouillon, F. Ricci, M. Royo,
M. Schmitt, M. Stengel, M. Torrent, J. Van Bever, M. Verstraete, J. Zwanziger.

It is worth to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases.
This might take some time ...

Xavier

* * *

### **A.** Important remarks and warnings.

**A.1** At the occasion of the switch from ABINITv8 to ABINITv9, many improvements of the formats and content of files written
    by ABINIT have been made, so the backward compatibility of ABINITv9 may be broken.

In particular:

1. The build system relies on new `.ac9` files (see [B.6](#v9.2.B.6)), superceeding the v8 `.ac` files.
   A bash script (`upgrade-build-config-file.sh`) located in the top level directory of the package can be used
   to convert from the old `.ac`format to `.ac9`.
2. The build system of ABINITv9 does not build anymore the hard dependencies (Linalg, NetCDF4, HDF5, LibXC, ...),
   as this was not sustainable (see [B.6](#v9.2.B.6)) and nowadays most users install prerequisite libraries themselves.
   See also the specialized INSTALL notes for
   [CentOS](/INSTALL_CentOS), [EasyBuild](/INSTALL_EasyBuild), [MacOS](/INSTALL_MacOS), and [Ubuntu](/INSTALL_Ubuntu).
3. The main ABINIT output file now contains sections written in YAML (sometimes replacing text sections, sometimes adding information).
   This means that some user-developed parsing tools might not work anymore,
   and should be adapted to the new ABINITv9 output file (see [B.9](#v9.2.B.9)).
   Note that the YAML output is still under development and modifications may appear in the next versions.
   A python API to extract the results of the calculation will be provided when the implementation is finalized.
4. Several default values have been changed, see [A.3](#v9.2.A.3).


**A.2**
A new account of the ABINIT effort has been published in Computer Phys. Comm. [[cite:Gonze2020]].
It provides description of several new features.
A version of this paper that is not formatted
for Computer Phys. Comm. [is also available](https://www.abinit.org/sites/default/files/ABINIT20.pdf).
The licence allows the authors to put it on the Web.

A second new account of the ABINIT effort has been published in J. Chem. Phys. [[cite:Romero2020]].
The scope of this second paper is different from the first one. It is more a survey of ABINIT,
focusing on its specific capabilities. Still, it contains also some description of some new features.
A version of this paper that is not formatted for J. Chem. Phys.
[is also available](https://www.abinit.org/sites/default/files/ABINIT20_JPC.pdf).
The licence allows the authors to put it on the Web.

Other specific publications are mentioned in the [Suggested acknowledgment page](/theory/acknowledgments).

<a name="v9.2.A.3"></a>
**A.3**  The default values of the following ABINIT input variables have been changed:
    [[ixcrot]], [[chneut]], [[ntime]], [[prtkden]], [[symsigma]] and [[tolsym]]. In particular the new default value
    of [[tolsym]], 1e-5, is more in line with the tolerances of other codes, so that for users of such
    codes, one barrier to the use of ABINIT is removed. By the same token, some bug in the recognition of symmetries
    has been fixed, when [[tolsym]] is close to the default, see the new tests [[test:v9_190]] and [[test:v9_191]].
    The new input variable [[chksymtnons]] has been introduced, to govern the possible automatic alignment
    of the [[tnons]] with the FFT grid (actually needed for GW calculations).
    By X. Gonze (MR 689 and others)

**A.4** The initialization of the wavefunctions when [[paral_kgb]]=1 and [[nspinor]]=2 has been changed, since the previous one could prevent the code to converge.
    By M Torrent (MR 562).

**A.5** The input variable xangst has been disabled. Use [[xcart]] instead, and specify the unit, namely Angstrom.

**A.6** Work is on-going concerning the Coulomb singularity treatment, see [D.32](#v9.2.D.32). The usage of the input variable [[icutcoul]] is changing.
For the time being, use [[gw_icutcoul]] instead.

**A.7** The name of the t-DEP main executable has been changed from `tdep` to `atdep`, in line with [[cite:Romero2020]].
By J. Bieder (MR 642, 641).

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
with many tests, and three tutorials provided.

List of tests: [[test:v8_44]], [[test:v9_50]], [[test:v9_53]], [[test:v9_56]], [[test:v9_60]], [[test:v9_61]],
[[test:eph4mob_1]], [[test:eph4mob_2]], [[test:eph4mob_3]],
[[test:eph4mob_4]],
[[test:eph4mob_5]],
[[test:eph4mob_6]],
[[test:eph4mob_7]],
[[test:eph4zpr_1]],
[[test:eph4zpr_2]],
[[test:eph4zpr_3]],
[[test:eph4zpr_4]],
[[test:eph4zpr_5]],
[[test:eph4zpr_6]],
[[test:eph4zpr_7]].

New input variables: [[brav]], [[dvdb_add_lr]], [[dvdb_qdamp]],
[[dvdb_rspace_cell]], [[eph_doping]],
[[eph_phrange]], [[eph_tols_idelta]], [[eph_ecutosc]], [[eph_restart]],
[[eph_stern]], [[eph_use_ftinterp]], [[eph_phwinfact]],
[[getdvdb]], [[getdvdb_filepath]], [[getkerange_filepath]],
[[getsigeph_filepath]],
[[irddvdb]], [[prteliash]], [[rifcsph]], [[sigma_bsum_range]], [[sigma_erange]],
[[sigma_ngkpt]], [[sigma_nshiftk]], [[sigma_shiftk]], [[symv1scf]].

Note that the new EPH processing unit of ABINIT [[optdriver]]=7 has a different implementation than the one implemented in anaddb.
Three new tutorials are availables, [[tutorial:eph_intro]], [[tutorial:eph4mob]] and [[tutorial:eph4zpr]], and supercede the legacy tutorials
[[tutorial:eph_legacy]] and [[tutorial:eph_tdep_legacy]].
For further details about the implementation and usage, please consult [[cite:Brunin2020b]].

By G. Brunin, H. Miranda, M. Giantomassi, G.-M. Rignanese, G. Hautier.

**B.2** Flexoelectricity and dynamical quadrupoles

A new driver has been included in abinit that allows one to compute
4 spatial dispersion tensorial quantities: the clamped-ion flexoelectric tensor, the dynamical quadrupoles,
the first moment of IFC matrix and the first moment of the piezoelectric force response tensor.
Precalculation of ground state, first and second (d2_dkdk) order response functions is required.
After execution, the driver creates a 3rd order energy derivative database file
that is used by anaddb to compute the mixed and lattice-mediated flexoelectric tensors
or to include the dipole-quadrupole and quadrupole-quadrupole electrostatic interactions
in the calculation of the dynamical matrix.

See the complementary description
in the Sec. V. D of [[cite:Romero2020]], with underlying theory and test calculations
presented in [[cite:Royo2019]].
At the practical level, see [[cite:Romero2020]]:

>   In this way, both perturbations are generalized to finite q, as
>   is already the case for atomic displacements. This enables us
>   to carry out an analytical third order derivative of the energy
>   with respect to two of the standard perturbations, and to the
>   momentum q, which directly provides the sought-after spatial
>   dispersion tensors. Remarkably, by virtue of the 2n+1 theorem,
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
>   Ref. 34 have been adapted to incorporate the dipole-quadrupole
>   and quadrupole-quadrupole electrostatic interactions
>   derived in Ref. 102. This new functionality results in a
>   faster convergence of the phonon bands calculation with respect
>   to the density of q points and, in some materials, represents
>   the only route to obtain the correct sound velocities.

>   Currently, the implementation is restricted to the use of
>   norm-conserving pseudopotentials without non-linear core corrections, and the LDA
>   functional.

A tutorial is in preparation, with tests [[test:lw_1]] to [[test:lw_7]].

See the [[topic:longwave]]. The relevant input variable is [[optdriver]]==10.
New input variables have been defined: [[lw_flexo]], [[lw_qdrpl]], [[prepalw]], [[flexoflag@anaddb]],
[[dipquad@anaddb]], [[quadquad@anaddb]].

This capability is still under development and not completely stable.
Interested users are strongly recommended to contact Miquel Royo (mroyo@icmab.es)
or Massimiliano Stengel (mstengel@icmab.es) before start using it.

By M. Royo, M. Stengel


**B.3** DFT+DMFT

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


**B.4** Spin model within Multibinit

The new capabilities of Multibinit within ABINITv9 are described
fully in the Sec. 4.1 of [[cite:Gonze2020]]. See also Sec. [D.1](#v9.2.D.1).
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

Many new input variables are present, tested and documented:
[[slc_coupling@multibinit|slc_coupling]],
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

List of tests in addition to those of the tutorial: [[test:v8_16]], [[test:v8_23]], [[test:v9_81]], [[test:v9_82]],
[[test:v9_86]], [[test:v9_87]].

By Xu He, N. Helbig, J. Bieder, E. Bousquet, Ph. Ghosez, M. Verstraete

<a name="v9.2.B.5"></a>
**B.5** Constrained DFT

Constrained Density-Functional Theory (see [[topic:ConstrainedDFT]]) is available,
with a new algorithm allowing to impose the constraints to arbitrary precision,
whether it relates to the charge, magnetization, magnetization direction, or magnetisation size,
or a combination thereof for different atoms. The constraints are smeared spherical integrals
with ajustable sphere radius, centered on atoms. The algorithms has been demonstrated for norm-conserving pseudopotentials
as well as PAW. Forces and derivatives with respect to the constraints
are available (i.e. magnetic torque for the non-collinear spin case).
Stresses are still to be coded, will be available in ABINITv9.4 or ABINITv9.6.

New tests: v8#24-29, v8#95-97 and v9#1-3.
New input variables: [[chrgat]], [[constraint_kind]], [[ratsm]].

By X. Gonze.

<a name="v9.2.B.6"></a>
**B.6** Large modifications of the build system

The build system relies on new <hostname>.ac9 files, superceeding the v8 <hostname>.ac files.
Fully documented example files can be found in doc/build/config-examples.
A bash script (`upgrade-build-config-file.sh`) located in the top level directory of the package can be used
   to convert from the old `.ac`format to `.ac9`.

The build system of ABINITv9 does not build anymore the (hard and soft) dependencies (Linalg, NetCDF4, HDF, LibXC, Wannier90, ...), as this was not sustainable.
Three libraries are now mandatory: linalg, NetCDF4/HDF5 and LibXC. Failing to link to them will prevent building ABINIT.
The other libraries are optional, there will only be a warning if they are not available.
If the user does not provide the path to these libraries,
the build system will try to find them in the "usual" directories, and inform the user that it has done so.
The build system also can make suggestions to the user, to complete its *.ac9 file.

Specialized INSTALL notes are available to help the user for
[CentOS](/INSTALL_CentOS), [EasyBuild](/INSTALL_EasyBuild), [MacOS](/INSTALL_MacOS), and [Ubuntu](/INSTALL_Ubuntu).

By Y. Pouillon and JM Beuken

**B.7** New command line interface

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
test v7_94 for the (old implementation) of electron-phonon calculations in anaddb.

!!! important

    The old "files file" interface is still operational although deprecated and will be **REMOVED** in Abinit v10.


By M. Giantomassi (MR 586).


**B.8** Reading strings from the input file

A new mechanism to read strings enclosed between **double quotation marks** from the input file has been activated.
So, many new input keywords are reading strings as data, and, often, can be used alternatively to similar input keywords
that were expecting numerical values such as the `get*` and `ird*` variables.
The goal is to encourage a new approach for performing ABINIT calculations in which multiple datasets
and `get*` variables are replaced by independent input files that are connected together via file paths.

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
- [[pp_dirpath]]: Used in all the input files of the Test Suite
- [[output_file@anaddb]], [[ddb_filepath@anaddb]], [[gkk_filepath@anaddb]], [[eph_prefix@anaddb]].
  See tests [[test:v8_52]] for a standard analysis of the DDB file and
  test v7_94 for the (old implementation) of electron-phonon calculations in anaddb.

By M. Giantomassi

<a name="v9.2.B.9"></a>
**B.9** YAML sections in the output file

YAML sections are now generated in the output file, sometimes replacing text sections, sometime providing new information.
At present there is a YAML section for the components of the total energy, the GS results including forces and stresses as well as a YAML section for GW calculations, and some YAML sections giving information about the iteration status.
<!--
Example of tests: paral#86, v67mbpt#2. See the input variable use_yaml (TO BE DOCUMENTED).
-->
At the occasion of the development of this capability, and its adaptation to the test farm, the
perl script fldiff.pl has been replaced by a Python version.
See related information in Sec. 5.5 of [[cite:Gonze2020]].

By T. Cavignac, M. Giantomassi, GM Rignanese, X Gonze.


**B.10** New approach to define crystalline structures in the Abinit input

The new variable [[structure]] can be used to initialize the lattice vectors
and the atomic positions **from an external file**.
Variables such as [[natom]], [[ntypat]], [[typat]] and [[znucl]] are automatically initialized
and need not to be specified in the ABINIT input.
At present, the code can read netcdf files produced by ABINIT (`GSR.nc`, `WFK.nc`, `DEN.nc`, `HIST.nc`)
and POSCAR files in VASP-5 format.
See the documentation for the syntax and limitations.

By M. Giantomassi.


**B.11** New capabilities of abipy and abiflows

The abipy and abiflows projects have been significantly extended.
See Sec. 6 of [[cite:Gonze2020]], as well as the [gallery of plotting scripts](http://abinit.github.io/abipy/gallery/index.html) &nbsp;
and the [gallery of abipy workflows](http://abinit.github.io/abipy/flow_gallery/index.html) &nbsp;.

By M. Giantomassi, G. Petretto, F. Naccarato.


* * *

### **C.** Changes for the developers (including information about compilers)

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

     added_in_version="9.2.0"

or, for variables introduced prior to v9,

     added_in_version="before_v9"

**C.2** Test farm: new and obsolete bots

* Bots introduced in the test farm: alps_nag_7.0_openmpi, atlas_gnu_9.1_openmpi, buda2_gnu_8.2_cuda, cronos2_gnu_7.4_paral,
    scope_gnu_10.2_mpich, scope_gnu_7.5_dep.
* Bots upgraded: abiref_gnu_5.3_* to abiref_gnu_9.2_*; bob_gnu_5.3_openmp to bob_gnu_7.5_openmp; buda2_gnu_8.1_mpich3
    to buda2_gnu_8.2_mpich3; graphene_gnu_6.4_macports to graphene_gnu_9.2_macports; max2_gnu_5.3_* to max2_gnu_6.5_*; ubu_gnu_5.3_openmpi to ubu_gnu_9.2_openmpi.
* Bots removed: abiref_nag_6.2_openmpi (superceded by alps_nag_7.0_openmpi), atlas_gnu_7.2_fb.ac (no nightly test of the fallbacks anymore), cronos_gnu_5.3_paral(replaced by cronos2), inca_gnu_6.3_py3k (inca too old), tikal_gnu_* (tikal too old).
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

**C.9** Reorganisation of the tests/Psps_for_tests directory, in preparation of the beautification.
Clarification on which pseudopotentials are of recent format, and which are legacy pseudopotentials.
Preparation of beautification.

By X. Gonze (MR 683)

* * *

<a name="v9.2.D.1"></a>
### **D.**  Other changes (or on-going developments, not yet finalized)

**D.1**  Miscellaneous improvements of Multibinit (lattice part)

Miscellaneous improvements have been made to the lattice part of Multibinit.
See the new input variables below, also see the Sec. 4.1.1 of [[cite:Gonze2020]].

New tests: [[test:v8_38]], [[test:v8_94]], [[test:v8_98]], [[test:v8_99]],
[[test:v8_101]],
[[test:v8_102]],
[[test:v8_103]],
[[test:v9_83]],  [[test:v9_84]], [[test:v9_85]].
New input variables are listed below.
Not all these new input variables are present in automatic tests, though.

- [[analyze_anh_pot@multibinit|analyze_anh_pot]]
- [[dyn_chksym@multibinit|dyn_chksym]]
- [[dyn_tolsym@multibinit|dyn_tolsym]]
- [[fit_efs@multibinit|fit_EFS]]
- [[fit_iatom@multibinit|fit_iatom]]
- [[fit_spc_maxs@multibinit|fit_SPC_maxS]]
- [[latt_friction@multibinit|latt_friction]]
- [[latt_taup@multibinit|latt_taup]] NOT TESTED
- [[latt_taut@multibinit|latt_taut]]
- [[opt_coeff@multibinit|opt_coeff]]
- [[opt_effpot@multibinit|opt_effpot]]
- [[opt_ncoeff@multibinit|opt_ncoeff]]
- [[sel_efs@multibinit|sel_EFS]]
- [[strfact@multibinit|strfact]]
- [[test_effpot@multibinit|test_effpot]]
- [[test_prt_ph@multibinit|test_prt_ph]]
- [[ts_option@multibinit|ts_option]] NOT TESTED

Finally, several input variables of the main ABINIT application are also reused for Multibinit, without modification
of meaning, like [[iatfix]], [[natfix]] (and related similar input variables), and also [[tolmxf]].

By M. Schmitt, Xu He, F. Ricci, M. Verstraete, Ph. Ghosez


**D.2** Miscellaneous improvements in the Chern number and orbital magnetization calculations,
including parallelization over k points of the Chern number calculation.

By J. Zwanziger (MR 469, 500, 545, 588)

**D.3** Calculation of Debye-Waller tensor. New test [[test:v8_58]].

By M. Giantomassi

**D.4** Speed-up of susceptibility matrix calculations and GW analytic continuation calculations.
See the new input variable [[gwaclowrank]] and new test [[test:v9_32]].

By F. Bruneval (MR 687).

**D.5** NCPP Wavefunction mixing with Variational Energy
and minor improvements to prepare PAW+Hybrid variational energy.
New test [[test:v7_73]], simple system for testing Hartree-Fock and the SCF algorithms.

By X. Gonze (MR 434, 444, 445).

**D.6** New weight distribution of the Fourier components of the force constants.
Test tolerance in the new integration weights, tests [[test:v8_52]], [[test:v8_53]], [[test:v8_54]].

By H. Miranda and M. Giantomassi

**D.7** Test calculation of velocity matrix elements (DDK) with
 [[optdriver]] 8 and [[wfk_task]] "wfk_ddk”, see [[test:v8_59]]. By the way, the
 other capabilities linked to [[wfk_task]] ("wfk_fullbz", "wfk_einterp", "wfk_optics_fullbz", "wfk_kpts_erange") seem
 not to have been properly advertised.

By M. Giantomassi

**D.8** Upgraded [[tutorial:paral_bandpw]], new version of auto paral (with threads)

By M. Torrent (MR502).

**D.9** Test wannier90 interface with [[nsppol]]=2 and [[nspden]]=2, [[test:wannier90_04]].

By Xu He

**D.10** Mixed precision for FFT transforms. New input variable [[mixprec]]
see [[test:v8_44]], [[test:v9_57]], [[test:v9_60]], and [[test:v9_61]].

From M. Giantomassi (MR491).

**D.11** Multibinit has been interfaced with [scale-up](https://www.secondprinciples.unican.es).

By Marcus Schmitt, Jordan Bieder, Matthieu Verstraete and Philippe Ghosez

**D.12** The following units are now also allowed in input files:

- S Sec Second  for the ABINIT input file;
- nm (for nanometer)  for the ABINIT and ANADDB input files.

**D.13** aTDEP utility:
added [[pdf:aTDEP_Guide| aTDEP guide]],
[[topic:aTDEP|aTDEP topic]], and corresponding input variable documentation.
References: [[pdf:aTDEP_Paper|aTDEP paper]].
Also, see Sec. 4.2 of [[cite:Gonze2020]].

By F. Bottin, J. Bouchet, J. Bieder (MR491,422).

**D.14** Improvements of NLO calculations.
Optimize memory allocation. In particular for [[usepead]] = 1.
Write Raman susceptibilities to netcdf in anaddb.

By G. Petretto  (MR 599).

**D.15** MetaGGA + PAW in progress.
All internal and unit tests are OK. The implementation seems close to the end. Still need to perform physical tests.

By M. Torrent and J.-B. Charraud (MR 587, 558, 625).

**D.16** Implementation of an alternate version of MPI reductions valid when the number of data exceeds 2^32.

By M. Torrent (MR 587)

**D.17** Build system and src support for llvm (clang/flang) and ARM (armflang/armclang/armpl).

By M Torrent (MR 571)

**D.18** Improvement of core WF reading (optics + electron-positron)

By M Torrent (MR 557)

**D.19** Fix bootstrap kernel convergence issues

By Wei Chen (MR 546)

**D.20** Upgrade versions of libraries used with ABINIT.
Upgrade atompaw to 4.1.0.6. Upgrade Libxc to 4+. Prepare the interface to LibXC 5 and a bit of LibXC 6.

By M. Torrent, JM Beuken (MR 649, 532, 470, 465, 441)

**D.21** Write yaml file for fatbands (phonopy format) with aTDEP

By J. Bieder (MR510)

**D.22** Write polarization vectors in GSR.nc.

By H. Miranda (MR 462).

**D.23** Updated user interface of Raman_spec.py . Added one section to tutorial NLO to describe use of Raman_spec.py.

By N. Pike (MR 460, MR 581)

**D.24** XML format for core wavefunctions

By F. Jollet (MR 423)

**D.25** Wavefunction prediction for molecular dynamics.

By F. Jollet (MR 412)

**D.26** Added a preview for the toptic_4.abi file in the optic tutorial.

By F. Goudreault (MR 408)

**D.27** Improve ELPA detection; make abinit compatible with ELPA2019

By M. Torrent (MR 626)

**D.28** Upgrade of [[tutorial:base3]] and [[tutorial:basepar]].

By X. Gonze (MR628)

**D.29** New input variable [[prtprocar]], see test [[test:v9_108]].

By M. Verstraete (MR630)

**D.30** The ABINIT input variable [[supercell_latt]] is documented.
See also [[test:v8_94]].

By X. Gonze

**D.31** New implementation of the cRPA (to compute U and J in DFT+U or DFT+DMFT)
that can consider any subset of orbitals in the system.
Compilation was a bit long on certain compilers, so new keywords have been added to the build system:
enable_crpa_optim and enable_crpa_no_optim

By R. Outerov and B. Amadon (MR622).

<a name="v9.2.D.32"></a>
**D.32** On-going work on refactoring the Coulomb interaction part of ABINIT.

New input variables [[fock_icutcoul]], and [[gw_icutcoul]], that should superceed [[icutcoul]].
New test added for the mini-Brillouin Zone integration, [[gw_icutcoul]]=14, 15, 16, see [[test:v9_181]].

By B. Guster, M. Giantomassi, F. Bruneval and X. Gonze (MR 627, 633, 673, 679, 686).

**D.33** New TB2J python script to compute the superexchange (J) and the Dzyaloshinskii-Moriya (DMI) interactions. The script can be found in [http://gitlab.abinit.org/xuhe/TB2J](http://gitlab.abinit.org/xuhe/TB2J) with doc and tutorials. The script is interfaced with wannier90 and use the w90 output files. The J calculation works in production, the DMI is much more sensitive to disentanglement noise and have to be use with care. An article is under construction to describe the method and its implementation. The script can deliver input data file for the spin model of Multibinit.

By He Xu, M. Verstraete and E. Bousquet (MR 639).

**D.34** Test linear electro-optical coefficient, tutorespfn [[test:optic_5]].

By N. Pike (MR 575).


**D.35** Optic: print spin-decomposition when applicable.
Fixed reflectivity screwed results. Test more thoroughly optics (incl. interfacing using NetCDF),
see new tests [[test:v9_05]] to [[test:v9_12]], also [[test:v9_47]] and [[test:v9_48]]

By X. Gonze (MR 654, 674).

**D.36** Fixed DFPT+PAW+GGA+usexcnhat=1+q<> bug, and also add new related tests
[[test:v9_41]] and [[test:v9_42]].

By M. Torrent and X. Gonze (MR 657).

**D.37** Update PAW tutorial

By F. Jollet (MR 643, 645, 652, 653).

**D.38** Nonlinear xc preliminary

Preliminary work for some changes in the exchange correlation terms in nonlinear (3rd order DFPT).
One term was implemented in both pead and dfptnl routines. Now it is merged in one routine (dfptnl_exc3).
Due to a subtle reordering of nonlinear core correction terms, some pead refs are changed, but the final result ("First order change in electronic dielectric susceptibility tensor") remain the same.

By L. Baguet (MR 650).

**D.39** Anaddb output data prefix

The former writing of anaddb.nc (fixed path and fixed name) has been made more flexible, by
introducing a prefix

By J. Bieder and He Xu (MR702)

**D.40** New "macro" input variable [[expert_user]]

When non-zero [[expert_user]] automatically switch off all checks done by [[chkprim]], [[chkdilatmx]], [[chksymbreak]] and [[chksymtnons]].

By X. Gonze (MR715)


**D.41** Miscellaneous additional bug fixes, improvements of documentation including for the build system.
G. Antonius, L. Baguet, JM Beuken, J. Bieder, E. Bousquet, F. Bruneval, T. Cavignac, M. Giantomassi, X. Gonze,
F. Jollet, R. Outerovitch, N. Pike, Y Pouillon, M. Royo, M. Torrent, J. Van Bever, M. Verstraete, Xu He.

* * *

## v9.0

Version 9.0, released on March 29, 2020.
List of changes with respect to version 8.10.

Many thanks to the contributors to the ABINIT project between
October 2018 and March 2020. These release notes
are relative to modifications/improvements of ABINIT v9.0 with respect to v8.10
(merge requests up to, and including, MR636 are taken into account)

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

* * *

### **A.** Important remarks and warnings.

**A.1** At the occasion of the switch from ABINITv8 to ABINITv9, many improvements of the formats and content of files written
    by ABINIT have been made, so the backward compatibility of ABINITv9 may be broken.
    The present ABINITv9.0 is NOT to be considered a production version. It is a beta release, allowing developers to get feedback
    from the users. Many features will work correctly, of course. Still, beginners are advised
    to stick to ABINITv8.10.3 except if ABINITv8.10.3 is not appropriate (or not working) for them.

In particular:

1. The build system relies on new `.ac9` files (see [B.6](#v9.0.B.6)), superceeding the v8 `.ac` files.
   A bash script (`upgrade-build-config-file.sh`) located in the top level directory of the package can be used
   to convert from the old `.ac`format to `.ac9`.
2. The build system of ABINITv9 does not build anymore the hard dependencies (Linalg, NetCDF4, HDF5, LibXC, ...),
as this was not sustainable (see [B.6](#v9.0.B.6)) and nowadays most users install prerequisite libraries themselves.
3. The main ABINIT output file now contains sections written in YAML (sometimes replacing text sections, sometimes adding information).
    This means that some user-developed parsing tools might not work anymore, and should be adapted to the new ABINITv9 output file (see [B.9](#v9.0.B.9)). Note that the YAML output is still under development and modifications may appear in the next versions. A python API to extract the results of the calculation will be provided when the implementation is finalized.
4. Several default values have been changed, see [A.3](#v9.0.A.3).


**A.2**
A new account of the ABINIT effort has been published in Computer Phys. Comm. [[cite:Gonze2020]]
It provides description of several new features.
A version of this paper that is not formatted for Computer Phys. Comm.
[is also available](https://www.abinit.org/sites/default/files/ABINIT20.pdf).
The licence allows the authors to put it on the Web.

A second new account of the ABINIT effort has been published in J. Chem. Phys. [[cite:Romero2020]].
The scope of this second paper is different from the first one. It is more a survey of ABINIT,
focusing on its specific capabilities. Still, it contains also some description of some new features.
A version of this paper that is not formatted for J. Chem. Phys.
[is also available](https://www.abinit.org/sites/default/files/ABINIT20_JPC.pdf).
The licence allows the authors to put it on the Web.

Other specific publications are mentioned in the [Suggested acknowledgment page](/theory/acknowledgments).

<a name="v9.0.A.3"></a>
**A.3**  The default values of the following ABINIT input variables have been changed:
    [[ixcrot]], [[chneut]], [[ntime]], [[symsigma]], [[prtkden]].

**A.4** The initialization of the wavefunctions when [[paral_kgb]]=1 and [[nspinor]]=2 has been changed, since the previous one could prevent the code to converge.
    By M Torrent (MR 562).

**A.5** The input variable xangst has been disabled. Use xcart instead, and specify the unit, namely Angstrom.

**A.6** The name of the t-DEP main executable has been changed from `tdep" to `atdep`, in line with [[cite:Romero2020]].

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

New input variables:
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

A new driver has been included in abinit that allows one to compute
4 spatial dispersion tensorial quantities: the clamped-ion flexoelectric tensor, the dynamical quadrupoles,
the first moment of IFC matrix and the first moment of the piezoelectric force response tensor.
Precalculation of ground state, first and second (d2_dkdk) order response functions is required.
After execution, the driver creates a 3rd order energy derivative database file
that is used by anaddb to compute the mixed and lattice-mediated flexoelectric tensors
or to include the dipole-quadrupole and quadrupole-quadrupole electrostatic interactions
in the calculation of the dynamical matrix.

See the complementary description
in the Sec. V. D of [[cite:Romero2020]], with underlying theory and test calculations
presented in [[cite:Royo2019]].

At the practical level, see [[cite:Romero2020]]:

>   In this way, both perturbations are generalized to finite q, as
>   is already the case for atomic displacements. This enables us
>   to carry out an analytical third order derivative of the energy
>   with respect to two of the standard perturbations, and to the
>   momentum q, which directly provides the sought-after spatial
>   dispersion tensors. Remarkably, by virtue of the 2n+1 theorem,
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
>   Ref. 34 have been adapted to incorporate the dipole-quadrupole
>   and quadrupole-quadrupole electrostatic interactions
>   derived in Ref. 102. This new functionality results in a
>   faster convergence of the phonon bands calculation with respect
>   to the density of q points and, in some materials, represents
>   the only route to obtain the correct sound velocities.

>   Currently, the implementation is restricted to the use of
>   norm-conserving pseudopotentials without non-linear core corrections, and the LDA
>   functional.

A tutorial is in preparation, with tests [[test:lw_1]] to [[test:lw_7]],
as well as a specific topic.

New input variables have been defined: [[lw_flexo]], [[lw_qdrpl]], [[prepalw]], [[flexoflag@anaddb]],
[[dipquad@anaddb]], [[quadquad@anaddb]].

This capability is still under development and not completely stable.
Interested users are strongly recommended to contact Miquel Royo (mroyo@icmab.es)
or Massimiliano Stengel (mstengel@icmab.es) before start using it.

By M. Royo, M. Stengel, M. Giantomassi.


**B.3** DFT+DMFT

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


**B.4** Spin model within Multibinit

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
[[slc_coupling@multibinit|slc_coupling]] (NT),
spin_calc_correlation_obs (NT and not documented),
spin_calc_traj_obs (NT and not documented),
[[spin_projection_qpoint@multibinit|spin_projection_qpoint]] (NT).

By Xu He, N. Helbig, J. Bieder, E. Bousquet, Ph. Ghosez, M. Verstraete


**B.5** Constrained DFT

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

<a name="v9.0.B.6"></a>
**B.6** Large modifications of the build system

The build system relies on new <hostname>.ac9 files, superceeding the v8 <hostname>.ac files.
Fully documented example files can be found in doc/build/config-examples.
A bash script (`upgrade-build-config-file.sh`) located in the top level directory of the package can be used
   to convert from the old `.ac`format to `.ac9`.

The build system of ABINITv9 does not build anymore the (hard and soft) dependencies (Linalg, NetCDF4, HDF, LibXC, Wannier90, ...), as this was not sustainable.
Three libraries are now mandatory: linalg, NetCDF4/HDF5 and LibXC. Failing to link to them will prevent building ABINIT.
The other libraries are optional, there will only be a warning if they are not available.
If the user does not provide the path to these libraries,
the build system will try to find them in the "usual" directories, and inform the user that it has done so.
The build system also can make suggestions to the user, to complete its *.ac9 file.

By Y. Pouillon and JM Beuken

**B.7** New command line interface

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
test v7_94 for the (old implementation) of electron-phonon calculations in anaddb.

!!! important

    The old "files file" interface is still operational although deprecated and will be **REMOVED** in Abinit v10.


By M. Giantomassi (MR 586).


**B.8** Reading strings from the input file

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
  test v7_94 for the (old implementation) of electron-phonon calculations in anaddb.

By M. Giantomassi

<a name="v9.0.B.9"></a>
**B.9** YAML sections in the output file

YAML sections are now generated in the output file, sometimes replacing text sections, sometime providing new information.
At present there is a YAML section for the components of the total energy, the GS results including forces and stresses as well as a YAML section for GW calculations, and some YAML sections giving information about the iteration status.
<!--
Example of tests: paral#86, v67mbpt#2. See the input variable use_yaml (TO BE DOCUMENTED).
-->
At the occasion of the development of this capability, and its adaptation to the test farm, the
perl script fldiff.pl has been replaced by a Python version.
See related information in Sec. 5.5 of [[cite:Gonze2020]].

By T. Cavignac, M. Giantomassi, GM Rignanese, X Gonze.


**B.10** New approach to define crystalline structures in the Abinit input

The new variable [[structure]] can be used to initialize the lattice vectors
and the atomic positions **from an external file**.
Variables such as [[natom]], [[ntypat]], [[typat]] and [[znucl]] are automatically initialized
and need not to be specified in the ABINIT input.
At present, the code can read netcdf files produced by ABINIT (`GSR.nc`, `WFK.nc`, `DEN.nc`, `HIST.nc`)
and POSCAR files in VASP-5 format.
See the documentation for the syntax and limitations.

By M. Giantomassi.


**B.11** New capabilities of abipy and abiflows

The abipy and abiflows projects have been significantly extended.
See Sec. 6 of [[cite:Gonze2020]], as well as the [gallery of plotting scripts](http://abinit.github.io/abipy/gallery/index.html) &nbsp;
and the [gallery of abipy workflows](http://abinit.github.io/abipy/flow_gallery/index.html) &nbsp;.

By M. Giantomassi, G. Petretto, F. Naccarato.


* * *

### **C.** Changes for the developers (including information about compilers)

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
* Bots upgraded: abiref_gnu_5.3_* to abiref_gnu_9.2_*; abiref_nag_6.2_openmpi to
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

By J. Zwanziger (MR 469, 500, 545, 588)

**D.3** Calculation of Debye-Waller tensor. New test [[test:v8_58]].

By M. Giantomassi


**D.4** Test linear electro-optical coefficient, tutorespfn [[test:optic_5]].

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

**D.8** Upgraded [[tutorial:paral_bandpw]], new version of auto paral (with threads)

By M. Torrent (MR502).

**D.9** Test wannier90 interface with [[nsppol]]=2 and [[nspden]]=2, [[test:wannier90_04]].

By Xu He

**D.10** Mixed precision for FFT transforms. New input variable [[mixprec]]
see [[test:v8_44]], [[test:v9_57]], [[test:v9_60]], and [[test:v9_61]].

From M. Giantomassi (MR491).

**D.11** Multibinit has been interfaced with scale-up,
[[https://www.secondprinciples.unican.es]]

By Marcus Schmitt, Jordan Bieder, Matthieu Verstraete and Philippe Ghosez

**D.12** The following units are now also allowed in input files:

- S Sec Second  for the ABINIT input file;
- nm (for nanometer)  for the ABINIT and ANADDB input files.

**D.13** aTDEP utility:
added [[pdf:aTDEP_Guide| aTDEP guide]],
[[topic:aTDEP|aTDEP topic]], and corresponding input variable documentation.
References: [[pdf:aTDEP_Paper|aTDEP paper]].
Also, see Sec. 4.2 of [[cite:Gonze2020]].

By F. Bottin, J. Bouchet, J. Bieder (MR491,422).

**D.14** Improvements of NLO calculations.
Optimize memory allocation. In particular for [[usepead]] = 1.
Write Raman susceptibilities to netcdf in anaddb.

By G. Petretto  (MR 599).

**D.15** MetaGGA + PAW in progress.
All internal and unit tests are OK. The implementation seems close to the end. Still need to perform physical tests.

By M. Torrent and J.-B. Charraud (MR 587, 558, 625).

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

**D.21** Write yaml file for fatbands (phonopy format) with aTDEP

By J. Bieder (MR510)

**D.22** Write polarization vectors in GSR.nc.

By H. Miranda (MR 462).

**D.23** Updated user interface of Raman_spec.py . Added one section to tutorial NLO to describe use of Raman_spec.py.

By N. Pike (MR 460, MR 581)

**D.24** XML format for core wavefunctions

By F. Jollet (MR 423)

**D.25** Wavefunction prediction for molecular dynamics.

By F. Jollet (MR 412)

**D.26** Added a preview for the toptic_4.abi file in the optic tutorial.

By F. Goudreault (MR 408)

**D.27** Improve ELPA detection; make abinit compatible with ELPA2019

By M. Torrent (MR 626)

**D.28** Upgrade of [[tutorial:base3]] and [[tutorial:basepar]].

By X. Gonze (MR628)

**D.29** New input variable [[prtprocar]], see test [[test:v9_108]].

By M. Verstraete (MR630)

**D.30** The ABINIT input variable [[supercell_latt]] is documented.
See also [[test:v8_94]].

By X. Gonze

**D.31** New implementation of the cRPA (to compute U and J in DFT+U or DFT+DMFT)
that can consider any subset of orbitals in the system.
Compilation was a bit long on certain compilers, so new keywords have been added to the build system:
enable_crpa_optim and enable_crpa_no_optim

By R. Outerov and B. Amadon (MR622).

**D.32** Work on refactoring the Coulomb interaction part of ABINIT.

By B. Guster, M. Giantomassi and X. Gonze (MR 627&633).

**D.33** New TB2J python script to compute the superexchange (J) and the Dzyaloshinskii-Moriya (DMI) interactions. The script can be found in [http://gitlab.abinit.org/xuhe/TB2J](http://gitlab.abinit.org/xuhe/TB2J) with doc and tutorials. The script is interfaced with wannier90 and use the w90 output files. The J calculation works in production, the DMI is much more sensitive to disentanglement noise and have to be use with care. An article is under construction to describe the method and its implementation. The script can deliver input data file for the spin model of Multibinit.

By He Xu, M. Verstraete and E. Bousquet

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

* * *

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

C.2 New [Howto for developers](/developers/developers_howto) (variables, mkparents, robodoc, test_suite).
    By M. Giantomassi.

C.3 New [Howto for the test suite](/developers/testsuite_howto).
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

* * *

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
    Entry point: see the new header of any ABINIT documentation file e.g. the [new user's guide](/guide/new_user).
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
    see test:v7_93, although there is no real reason for using it.
    The new algorithm is very general, respect better the symmetries, and should even supercede
    the use of other values of [[brav@anaddb]].
    By G. Petretto following discussions with GM Rignanese and XGonze, and tests by Henrique Pereira Miranda.

B.6 The Chern number can be computed, in the norm-conserving case as well as in the PAW case.
    See the theory in [[cite:Ceresoli2006]].
    Associated input variable: [[orbmag]].
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

D.3 The [[tutorial:eph_tdep_legacy|tutorial on temperature-dependence of the electronic structure]] has been upgraded, and carefully tested.
    See all tests in `tutorespfn/teph_tdep_legacy*`.
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

D.10 On-going development: main executable `atdep`, for the TDEP algorithm, by Hellman and coworkers.
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
    Entry point: see the new header of any ABINIT documentation file e.g. the [new user's guide](/guide/new_user)
    By F. Jollet and X. Gonze (also tests/fixes by B. Amadon, M. Torrent).

B.2 A central [[theory:bibliography]] database abiref.bib has been created, and linked to the
    above-mentioned topics (B.1) but also to other parts of the ABINIT documentation (e.g. input variable list,
    the tutorials, the theory documents, the acknowledgments).
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

D.9 Test of non-magnetic DFT+U and DFT+U+SO.
    See the new test v5#16
    By M. Torrent.

D.10 Make DFT+U and local EX-exchange compatible with nspden=1/nspinor=2
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
    [[tests/libxc/Input/t41.abi]], [[test:v8_04]].
    Work by M. Giantomassi

B.2 The Fock mixing factor for the HSE hybrid functional can be tuned thanks to the input variable gwfockmix.
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
