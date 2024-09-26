"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
]

#: List of keywords that are automatically added to all the tests of this suite.
keywords = [
]

#: List of input files
inp_files = [
    # Constrained DFT
    "t01.abi",
    "t02.abi",
    "t03.abi",
    # Structure variable
    "t04.abi",
    # Optic
    "t05.abi",  # GGA-PBE nsppol=1
    "t06.abi",  # optic
    "t07.abi",  # GGA-PBE nsppol=2
    "t08.abi",  # optic
    "t09.abi",  # GGA-PBE PAW nsppol=1
    "t10.abi",  # optic
    "t11.abi",  # GGA-PBE PAW nsppol=2
    "t12.abi",  # optic
    "t13.abi",  # GGA-PBE NC BCC Iron ferromagnetic
    "t14.abi",  # GGA-PBE PAW BCC Iron ferromagnetic

    "t19.abi",  # disable all checks thanks to expert_user
    "t22.abi",  # test different cellcharge for different images, algo pSIC
    "t29.abi",  # RMM-DIIS eigsolver for NC.
    "t30.abi",  # RMM-DIIS eigsolver for PAW.
    # GW/BSE
    "t31.abi",  # HF exchange checking q->0 terms
    "t32.abi",  # AC GW
    "t33.abi",  # GW 1RDM and related quantities
    "t34.abi",  # Same as t33.in but reading checkpoints
    "t35.abi",  # GW 1RDM and related quantities (using only Sigma_x)
    "t36.abi",  # GW 1RDM and related quantities but using Silicon
    "t37.abi",  # GW 1RDM and related quantities but using Silicon with diff. bdgw values
    "t40.abi",  # chi0 with inclvkb=2
    # DFPT
    "t41.abi",
    "t42.abi",
    "t43.abi",
    "t44.abi",  # test orbital magnetism with DDK wavefunctions on Ne atom
    "t46.abi",  # longwave GGA
    "t47.abi",  # metallic iron GGA-PBE PAW nsppol=2
    "t48.abi",  # optic
    "t49.abi",  # optic (same as t14 but prtlincompmatrixelements = 1)
    # Electron-phonon
    "t50.abi",
    "t51.abi",
    "t52.abi",
    "t53.abi",
    "t54.abi",
    "t55.abi",
    "t56.abi",
    "t57.abi",
    "t58.abi",
    "t59.abi",
    "t60.abi",
    "t61.abi",
    "t62.abi",
    "t63.abi",
    "t64.abi",
    "t65.abi",
    "t66.abi",
    "t67.abi",

    # More ground state
    "t71.abi",  # test cprj_update_lvl and nloc_alg options, istwfk=1
    "t72.abi",  # test cprj_update_lvl and nloc_alg options, istwfk>=2
    "t73.abi",  # test cprj_update_lvl options, forces and stress computed at the end of the run
    "t74.abi",  # test cprj_update_lvl options, forces computed during SCF iterations
    "t75.abi",  # test useylm=1 for NCPP with all KB energies being negative
    # test usepawu options (including negative ones), nsppol=nspinor=nspden=1
    "t76.abi",
    # test usepawu options (including negative ones), nsppol=2,nspden=2
    "t77.abi",
    # test usepawu options (including negative ones), nspinor=2,nspden=4
    "t78.abi",
    # test usepawu options (including negative ones), nspinor=2,nspden=1
    "t79.abi",

    # Spin dynamics in multibinit
    "t81.abi",  # set initial spin using rotation q-modulation
    "t82.abi",  # damping
    # New lattice mover in multibinit
    "t83.abi",  # Langevin NVT
    "t84.abi",  # Berendsen NVT
    "t85.abi",  # Velocity Verlet NVE
    "t86.abi",  # Spin lattice coupling, Only Oiju term activated.
    "t87.abi",  # Spin lattice coupling, Oiju and Tijuv.

    # test usepawu options with nspinor=2,nspden=4 and different pawxcdev
    "t88.abi", # LDA
    "t89.abi", # GGA

    # GS Coulomb cut-off
    "t90.abi",  # checkG Coulomb cut-off, large tolerance a.t.m.
    "t91.abi",  # occopt 9 tests on Si
    "t92.abi",  # check extended fpmd routines with low number of bands
    "t93.abi",  # energy, forces for PAW non-collinear, with usexcnhat=0
    "t94.abi",  # energy, stress for PAW non-collinear, with usexcnhat=0
    "t95.abi",  # energy, convergence for NC mGGA input with UPFv2.0 format psp

    # Optics with spin-orbit coupling
    "t100.abi", # Optical conductivity with spin-orbit coupling - ABINIT step
    "t101.abi", # Optical conductivity with spin-orbit coupling - CONDUCTI step
    "t102.abi", # X-ray core spectrocopy with spin-orbit coupling - ABINIT step
    "t103.abi", # X-ray core spectrocopy with spin-orbit coupling - CONDUCTI step
    # Optcell test
    "t104.abi", # Testing optcell 6 to relax 3rd vector not orthogonal to 1st and 2nd vectors

    # PAW+U and lruj/ujdet utilities
    "t105.abi", # Preliminary step for tests 106-109; generate WFK and DEN
    "t106.abi", # Hubbard U (macro_uj=1) test of ujdet subroutines
    "t107.abi", # Hund's J (macro_uj=4); preliminary step for t109
    "t108.abi", # Test of irdden and prtdosm input keywords
    "t109.abi", # Test of lruj post-processing w/ four LRUJ.nc files from t107

    # Lattice Wannier function
    # scdm-k method for lattice wannier function in anaddb (disentangle option 2)
    "t110.abi",
    # projWF method for lattice wannier function in anaddb (disentangle option 2)
    "t111.abi",
    
    # local orbital magnetic moments calculation.
    "t112.abi", # test prt_lorbmag   printing local orbital magnetic moments in atomic spheres

    # GS PAW Hybrid functionals
    "t120.abi",  # test PBE0 and related functionals with PAW

    # UPF2 format for norm-conserving pseudopotentials
    "t130.abi",  # UPF2
    "t131.abi",  # UPF2 with SOC
    "t132.abi",  # Forces using Beigi 2D cut-off

    # more DFPT
    "t140.abi",  # test orbital magnetism with DDK wavefunctions on AlP solid
    "t141.abi",  # test orbital magnetism with DDK wavefunctions on AlP solid with nspinor 2
    "t142.abi",  # test orbital magnetism with DDK wavefunctions on AlP solid with nsppol 2
    "t143.abi",  # test orbital magnetism with DDK wavefunctions and metallic sodium
    "t144.abi",  # test orbital magnetism with discretized DDK wavefunctions on AlP solid
    "t145.abi",  # test quadrupole calculation with all negative KB energies As PSP
    "t146.abi",  # test antisymmetric vector potential response and quadrupoles of HF molecule
    "t147.abi",  # test total vector potential response and quadrupoles of HF molecule
    "t148.abi",  # test flexoelectric tensor calculation with ffnl_lw=1

    # write_files behaviour
    "t150.abi",  # test the behaviour of supravariable write_files

    # DDB file conversion
    "t160.abi",  # convert DDB from text to netcdf
    "t161.abi",  # convert DDB from netcdf to text

    # Symmetry recognition and possible corrections, including inaccurate POSCAR files
    "t180.abi",  # test treatment of inaccurate POSCAR file
    "t181.abi",  # test treatment of inaccurate POSCAR file
    "t183.abi",  # test treatment of inaccurate POSCAR file
    "t184.abi",  # test treatment of inaccurate POSCAR file
    "t185.abi",  # test treatment of inaccurate POSCAR file
    "t186.abi",  # test treatment of inaccurate POSCAR file
    "t187.abi",  # test treatment of inaccurate POSCAR file
    "t188.abi",  # test treatment of inaccurate POSCAR file
    "t189.abi",  # test post treatment of an increase of symmetry due to geometry optimization
    "t190.abi",  # check slight misalignment of rprim, with large tolsym
    "t191.abi",  # check slightly incorrect xred, with large tolsym
    "t192.abi",  # check slightly incorrect rprim and xred, yielding correction to tnons, although still not tolerated.
    "t193.abi",  # check slightly incorrect rprim and xred, yielding correction to tnons, tolerated.
    "t195.abi",  # test treatment of inaccurate POSCAR file
    "t196.abi",  # test treatment of inaccurate POSCAR file
    "t197.abi",  # test treatment of inaccurate POSCAR file
    "t198.abi",  # test treatment of inaccurate POSCAR file
    "t199.abi",  # test treatment of inaccurate POSCAR file

    # new tolwfr usage
    "t200.abi",  # tolwfr with other tolerances
    "t201.abi",  # tolwfr with tolvrs + nbdbuf with lobpcg
    "t202.abi",  # tolwfr_diago with tolvrs and tolvrs+tolwfr with lobpcg

    # nblock_lobpcg
    "t205.abi",  # test nblock_lobpcg input

    # Chebfi2
    "t206.abi",  # check that various algorithm give the same results (including wfoptalg=111). istwfk=1.
    "t207.abi",  # same as t206 but with nspinor=2
]
