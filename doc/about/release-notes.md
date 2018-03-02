
## v8.4

Many thanks to the contributors to the ABINIT project between
January 2017 and May 2017. 
These release notes are relative to modifications/improvements of ABINITv8.4 with respect to v8.2.

The list of contributors includes :
F. Altvater, G. Antonius, L. Baguet, J.-M. Beuken, J. Bieder, E. Bousquet, 
W. Chen, G. Geneste, M. Giantomassi, Y. Gillet, X. Gonze, F. Jollet, A. Martin, 
F. Naccarato, G. Petretto, S. Prokhorenko, F. Ricci, M. Torrent , M. Verstraete, J. Zwanziger

It is worth to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases ...
This might take some time ...

Xavier

* * *

Version 8.4, released on June 2, 2017.

Changes with respect to version 8.2 :

A. Warnings and important remarks 
B. Most noticeable achievements (for users)
C. Changes in the package, for developers
D. Other changes (or on-going developments, not finalized)

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
    old syntax: `Tlibxc#42, Tv8#04` replaced by `[[tests/libxc/Input/t41.in]], [[test:v8_04]]`
    [[tests/libxc/Input/t41.in]], [[test:v8_04]] 
    [[ac:abiref_gnu_5.3_debug.ac]]
    Work by M. Giantomassi

B.2 The Fock mixing factor for the HSE hybrid functional can be tuned thanks to the input variable gwfockmix  .
    (Warning : this is for the GW-type approach to electronic structure only, not for total energies)
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
    WARNING : This capability of ABINIT has not been fully tested. However, the
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

B.8 New algorithms for the displacement of nuclei ([[ionmov]]) :
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

C.1 Management of the test farm : the new bot ubu_intel_17_openmpi 
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

D.7 Ongoing work : Raman intensities, in the PAW case, using DFPT. 
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
     inside ANADDB, work by G Petretto :
     - Updates in the calculation of sound velocity, 
     - Small modifications for [[anaddb:nlflag]] == 3 and added some quantities to the anaddb netcdf file,
     - Fix a bug in the implementation of the new weights

D.16 Concerning Path Integral Molecular Dynamics with Quantum Thermal Bath :
     allow restart from history file.
     By M. Torrent

D.17 Refactored the computation of the electric dipole moment.
     By M. Torrent

D.18 The energy width for the bands in the Boltztrap intrans file is now automatically set.
     Previously a constant 0.4 Ryd, and the user should have checked by hand if
     the value was sufficient. Should be considered a bug fix.
     By M Verstraete

D.19 Numerous miscellaneous additional bug fixes and improvements of documentation by :
     F. Altvater, G. Antonius, J. Bieder, M. Giantomassi, F. Jollet, 
     G. Petretto, M. Verstraete, M. Torrent, J. Zwanziger. 

* * *

## v8.2
    
Many thanks to the contributors to the ABINIT project between
June 2016 and January 2017. These release notes
are relative to modifications/improvements of ABINITv8.2 with respect to v8.0.

Moreover, most of them are also described in the Computer Physics Communications 2016 ABINIT paper, 
doi:10.1016/j.cpc.2016.04.003
    
The list of contributors includes :
B. Amadon, G. Antonius, L. Baguet, J.-M. Beuken, J. Bieder, E. Bousquet, F. Bruneval,
W. Chen, M. Giantomassi, Y. Gillet, X. Gonze, G. Petretto, F. Jollet, A. Martin,
V. Planes, Y. Pouillon, T. Rangel, F. Ricci, M. Torrent , M. Verstraete

It is worth to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases ...
This might take some time ...

Xavier

* * *

Version 8.2, released on February 16, 2017.

Changes with respect to version 8.0 :

A. WARNINGS.
B. Most noticeable achievements (for users)
C. Changes in the package, for developers
D. Other changes (or on-going developments, not finalized)

* * *

A.  WARNINGS AND IMPORTANT REMARKS

A.0 The 2016 article by the ABINIT group is now mentioned in the acknowledgments :
    "Recent developments in the ABINIT software package. 
    Computer. Phys. Communications 205, 106 (2016)".
    See http://www.abinit.org/doc/helpfiles/for-v8.2/users/acknowledgments.html , as well as
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
    Test case : v8#02 .
    The working routines were based on the original implementation of J. Nocera 
    available on netlib.org.  They have been reshaped and translated into modern fortran, 
    then interfaced to ABINIT by F. Bruneval (sources in 45_geomoptim/m_lbfgs.F90).

B.2 A new tutorial is available : a lesson on the calculation of the effective interactions U and J 
    using constrained Random Phase Approximation (cRPA) for DFT+DMFT (or DFT+U) calculations.
    See doc/tutorial/lesson_ucalc_crpa.html as well as the automatic tests tutorial/tucrpa#1-5 .
    This lesson was prepared by B. Amadon.

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

C.1 The version control system that is used for the development of ABINIT has been changed : 
    the whole ABINIT project
    has been ported from bzr to git.
    Work by J.-M. Beuken, Y. Pouillon, M. Giantomassi, X. Gonze, 
    with discussions with many developers.

C.2 New versions of Fortran compilers have been integrated in the test farm:
    - intel 16.0
    - gnu 6.1 and 6.2
    - IBM xlf compiler 14.1
    - NAG 5.3
    Corresponding examples are available in doc/config/build-examples.
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
    New input variables : ddb_shiftq, eph_task, eph_transport, prtphdos, prtphsurf.
    See tests v7#88 and 89.
    Work by M. Giantomassi and G. Antonius.

D.3 The generation of k-point meshes with kptrlatt and shiftk is now tested.
    See test v8#03 .
    Work by M. Giantomassi

D.4 As a follow-up of the Achievement B3 in the release notes of ABINITv8.0 (DMFT + TRIQS),
    new input variables have been defined for DMFT : dmft_tolfreq and dmftctqmc_triqs_nleg.
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

D.9 Unit tests (fftprof) have been set up for the use of the MKL-DFTI routines : 
    unitary#tfftmkl_03 and 04.
    Work by M. Giantomassi

D.10 Ongoing work : Raman intensities, in the PAW case, using DFPT. 
     By L. Baguet and M. Torrent

D.11 prtvcbm working with all parallelizations.
     Tests mpiio 26:28 have been reactivated on 4 processors.
     By M. Torrent

D.12 On going work related to non-collinear DFPT.
     By F. Ricci, S. Prokhorenko, and E. Bousquet

D.13 MBPT : support for the commutator [Vnl r] in the case of NC pseudos 
     with more than one projector per l-channel has been added. 
     Tests in v67mbpt[40] (GW run with psp8 files).
     SOC is not yet available, though.
     The KSS file continues to use the old implementation to maintain backward compatibility 
     hence outkss won't produce the file if multiple-projectors are detected
     inclvkb 1 has been removed. Now the possible values are either 0 or 2
     By M Giantomassi

D.14 For Hirshfeld and Bader : doc and warning.
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
     and improvements of documentation by :
     G. Antonius, L. Baguet, J. Bieder, F. Bruneval,
     M. Giantomassi, Y. Gillet, G. Petretto, Y. Pouillon,
     M. Verstraete, M. Torrent (in particular, for DFPT+PAW).
