"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
"HAVE_WANNIER90",
]

#: List of keywords that are automatically added to all the tests of this suite. 
keywords = [
]

subsuites = [
"w90",
]

#: List of input files
inp_files = [
"tw90_1.abi",
"tw90_2.abi",
"tw90_3.abi",
"tw90_4.abi",
]
