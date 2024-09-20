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
    "t03.abi" ,  # same as v10[04] but with norm-conserving pseudos (with DFTI)
    "t04.abi" ,  # same as v9[206] but with istwfk>1 (with DFTI)
    "t05.abi" ,  # same as v10[04] but with norm-conserving pseudos (without DFTI)
    "t06.abi" ,  # same as v9[206] but with istwfk>1 (without DFTI)
    "t40.abi" ,  # test orbmag calculation when using spatial symmetries for GS nuclear dipole 
    "t81.abi" ,  # Short MD to test restart on next test
    "t82.abi" ,  # Test restart of MD from the HIST of previous test using restartxf -1
]
