## v9.0
  
Version 9.0, released on March 15, 2020.
List of changes with respect to version 8.10.

Many thanks to the contributors to the ABINIT project between
October 2018 and March 2020. These release notes
are relative to modifications/improvements of ABINIT v9.0 with respect to v8.10.
<TBU>
The merge request #408 is the first MR not reported in these release notes. Then, #410-411, #413-416 have also been included.
<ETBU>

The list of contributors includes:
<TO BE UPDATED>
B. Amadon, G. Antonius, L. Baguet, J.-M. Beuken, J. Bieder, E. Bousquet, F. Bruneval, Wei Chen, M. Cote,
J. Denier, G. Geneste, Ph. Ghosez, M. Giantomassi, O. Gingras, X. Gonze, F. Goudreault, Xu He, Y. Jia, F. Jollet,
A. Lherbier, A. Martin, H. Miranda, F. Naccarato, G. Petretto, N. Pike,
S. Ponce, Y. Pouillon, S. Prokhorenko, F. Ricci, M. Torrent, M. van Setten, B. Van Troeye, M. Verstraete, J. Zwanziger.
<END OF TO BE UPDATED>

It is worth to read carefully all the modifications that are mentioned in the present file,
and examine the links to help files or test cases.
This might take some time ...

Xavier

### A. Important remark and warnings.

A.1 At the occasion of the switch from ABINITv8 to ABINITv9, many improvements of the formats and content of files
    have been made, so the backward compatibility of ABINITv9 is often broken. 
    The present ABINITv9.0 is NOT to be considered a production version, but more a beta release, allowing to get feedback
    from the users. Many features will work correctly, of course. Still, beginners are advised
    to stick to ABINITv8.10.3 except if ABINITv8.10.3 is not appropriate (or not working) for them.
    In particular, the build system relies on new <hostname>.ac9 files (see XX), superceeding the v8 <hostname>.ac files.
    ABINITv9 does not build the dependencies (Linalg, NetCDF, LibXC, ...) anymore, as this was not sustainable (see XX).
    The output file now contains sections written in YAML (sometimes replacing text sections, sometimes adding information),
    which means that some user-developed parsing tools might not work anymore, they have to be adapted to the new ABINITv9 output file. (see XX).

* * *

### B. Most noticeable achievements

* * *

### C. Changes for the developers (also compilers)

* * *

### D.  Other changes (or on-going developments, not yet finalized)

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
J. Denier, G. Geneste, Ph. Ghosez, M. Giantomassi, O. Gingras, X. Gonze, F. Goudreault, Xu He, Y. Jia, F. Jollet, 
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
G.-M. Rignanese, M. Torrent , M. Verstraete, J. Zwanziger

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

A.6 Definition of a maximal value for dilatmx, at 1.15 , than can be bypassed by setting chkdilatmx=0.
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

The list of contributors includes :
B. Amadon, G. Antonius, L. Baguet, J.-M. Beuken, J. Bieder, F. Bottin, Y. Bouchet, E. Bousquet,
M. Giantomassi, O. Gingras, Ph. Ghosez, M. Giantomassi, X. Gonze, F. Jollet, J. Junquera, A. Martin,
F. Naccarato, G. Petretto, N. Pike, Y. Pouillon, S. Prokhorenko, M. Torrent , M. Verstraete, J. Wiktor, J. Zwanziger

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

A.2 Some changes of names :

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
    Those topic web pages usually have :
    - a brief introduction;
    - the list of related tutorials -if any-;
    - the list of related input variables (ordered according to their importance for the topics -compulsory, basic, useful or expert-);
    - possibly example input files;
    - list of references.
    Entry point : see the new header of any ABINIT documentation file (e.g. the [new user's guide](..) )
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

D.6 On-going work on hybrid functionals : speed-up of the SCF loop, computation of stresses,
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

D.14 Ongoing work : Raman intensities, in the PAW case, using DFPT.
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

D.26 Numerous miscellaneous additional bug fixes and improvements of documentation by :
     G. Antonius, J. Bieder, M. Giantomassi, F. Jollet,
     G. Petretto, N. Pike, Y. Pouillon, M. Verstraete, M. Torrent.

* * *

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

List of changes with respect to version 8.0 .

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

* * *
