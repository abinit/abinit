"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
"HAVE_MPI_IO",
]

#: List of keywords that are automatically added to all the tests of this suite. 
keywords = [
]

#: This suite contains tests executed with different numbers of MPI processes.
is_multi_parallel = True

#: List of input files
inp_files = [
"t01.abi",
"t21.abi",
"t22.abi",
"t24.abi",
"t25.abi",
"t26.abi",
"t27.abi",
"t28.abi",
"t42.abi",
"t49.abi",
"t51.abi",
"t62.abi",
"t69.abi", 
"-t99.abi", # disabled since it have no meaning with my modifications
]
