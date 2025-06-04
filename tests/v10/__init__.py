"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
]

#: List of keywords that are automatically added to all the tests of this suite.
keywords = [
]

#: List of input files
inp_files = [
    "t01.abi" ,  # test post treatment of an increase of symmetry due to geometry optimization
    "t02.abi" ,  # check Wyckoff positions of trigonal groups 143-167
    "t03.abi" ,  # same as v10[04] but with norm-conserving pseudos (with DFTI)
    "t04.abi" ,  # same as v9[206] but with istwfk>1 (with DFTI)
    "t05.abi" ,  # same as v10[04] but with norm-conserving pseudos (without DFTI)
    "t06.abi" ,  # same as v9[206] but with istwfk>1 (without DFTI)
    "t07.abi" ,  # test oracle and nbdbuf in chebfi (cprj_in_memory=0 and 1).
    "t08.abi" ,  # same as v10[07] but with istwfk>1 (with less datasets)
    "t09.abi" ,  # same as v10[07] but with nspinor=2 (with less datasets)
    "t10.abi" ,  # compare cprj_in_memory=1 with cprj_in_memory=0. PAW, istwfk=1.
    "t11.abi" ,  # same as v10[10], with nsppol=2.
    "t12.abi" ,  # same as v10[10], with nspinor=2.
    "t13.abi" ,  # same as v10[10], with istwfk>1 (with DFTI).
    "t14.abi" ,  # same as v10[10] but with NC pseudos
#    "t15.abi" ,  # same as v10[10] but with NC pseudos, nspinor=2 ! not working yet
    "t16.abi" ,  # same as v10[10] but with NC pseudos, istwfk>1 (with DFTI)
    "t17.abi" ,  # same as v10[10], with istwfk>1 (without DFTI).
    "t18.abi" ,  # same as v10[10] but with NC pseudos, istwfk>1 (without DFTI)
    "t19.abi" ,  # compare cprj_in_memory=1 with cprj_in_memory=0 for cell optimization.
    "t20.abi" ,  # test nvt_langevin MD algorithm (PIMD implementation)
    "t21.abi" ,  # test npt_langevin MD algorithm (PIMD implementation)
    "t22.abi" ,  # test finite-temperature exchange-correlation functionals, and calculation of Sxc (NC case)
    "t23.abi" ,  # test finite-temperature exchange-correlation functionals, and calculation of Sxc (PAW case)
    "t24.abi" ,  # same as v10[10], with dilatxm>1, istwfk>1 (with DFTI).
    "t25.abi" ,  # same as v10[10], with dilatxm>1, istwfk>1 (without DFTI).
    "t26.abi" ,  # test cprj_in_memory when atoms are not ordered by type.
    "t40.abi" ,  # test orbmag calculation when using spatial symmetries for GS nuclear dipole
    "t41.abi" ,  # test orbmag calculation using R2SCAN mGGA
    "t42.abi" ,  # test orbmag calculation using R2SCAN mGGA, nspinor 2, zora
    "t43.abi" ,  # test quadrupoles calculation with xcnlcc
    "t81.abi" ,  # Short MD to test restart on next test
    "t82.abi" ,  # Test restart of MD from the HIST of previous test using restartxf -1
    "t83.abi" ,  # Test variable cell nudged elastic band method
    "t104.abi" , # Test Si spectral function and mobilities :step 1 WFK
    "t105.abi" , # Test Si spectral function and mobilities :step 2 merge DDB
    "t106.abi" , # Test Si spectral function and mobilities :step 3 merge DVDB
    "t107.abi" , # Test Si spectral function and mobilities :step 4 mobilities cumulant and DM with positive doping correction
]
