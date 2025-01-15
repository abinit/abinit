"""Global variables associated to the test suite."""

#: List of CPP variables that should be defined in config.h in order to enable this suite.
need_cpp_vars = [
"HAVE_GPU",
"HAVE_OPENMP_OFFLOAD",
]

#: List of keywords that are automatically added to all the tests of this suite. 
keywords = [
]

#: This suite contains tests executed with different numbers of MPI processes.
is_multi_parallel = True

#: List of input files
inp_files = [
"t01.abi", # PAW, GS, using CHEBFI and LOBPCG (AlP, nspinor==2, nspden==1)
"t02.abi", # PAW, GS+DFPT, only DDK (AlP, nspinor==2, nspden==1)
]
