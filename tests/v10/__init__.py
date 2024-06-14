"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
]

#: List of keywords that are automatically added to all the tests of this suite.
keywords = [
]

#: List of input files
inp_files = [
    # Fake test t01.abi, to initiate the tests/v10 directory . SHOULD BE REPLACED BY THE TEST OF A REAL NEW DEVELOPMENT.
    "t01.abi" ,  # test post treatment of an increase of symmetry due to geometry optimization
    "t02.abi" ,  # check Wyckoff positions of trigonal groups 143-167
    "t03.abi" ,  # same as v10[04] but with norm-conserving pseudos
    "t04.abi" ,  # same as v9[206] but with istwfk>1
    "t10.abi" ,  # compare cprj_in_memory=1 with cprj_in_memory=0. PAW, istwfk=1.
    "t11.abi" ,  # same as v10[10], with nsppol=2.
    "t12.abi" ,  # same as v10[10], with nspinor=2.
    "t13.abi" ,  # same as v10[10], with istwfk>1.
    "t14.abi" ,  # same as v10[10] but with NC pseudos
#    "t15.abi",  # same as v10[10] but with NC pseudos, nspinor=2 ! not working yet
    "t16.abi" ,  # same as v10[10] but with NC pseudos, istwfk>1
    "t40.abi" ,  # test orbmag calculation when using spatial symmetries for GS nuclear dipole 
]
