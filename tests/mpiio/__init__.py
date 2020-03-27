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
"t01.in",
"t21.in",
"t22.in",
"t24.in",
"t25.in",
"t26.in",
"t27.in",
"t28.in",
"t42.in",
"t49.in",
"t51.in",
"t62.in",
"t69.in", 
"-t99.in", # disabled since it have no meaning with my modifications
]
