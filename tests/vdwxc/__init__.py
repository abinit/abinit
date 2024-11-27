"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
"HAVE_WANNIER90",
]

#: List of keywords that are automatically added to all the tests of this suite. 
keywords = [
]

#: List of input files
inp_files = [
#"-t01.abi",
#"-t04.abi",
"t10.abi",
#"-t11.abi",
]
