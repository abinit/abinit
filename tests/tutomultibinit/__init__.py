"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
]

#: List of keywords that are automatically added to all the tests of this suite.
keywords = [
]

#: List of input files
inp_files = [
    "tmulti1_1.abi",
    "tmulti1_2.abi",
    "tmulti1_3.abi",
    "tmulti5_1.abi",   # basic Multibinit spin dynamics run
    "tmulti5_2.abi",   # LaFeO3 MvT
    "tmulti5_3.abi",   # 1D AFM chain with DMI (spin canting example)
]
