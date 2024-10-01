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
# Smoke tests on silicium, no forces nor stress (based on libxc#51)
"t01.abi", # Basic test, reference
"t02.abi", # Use GEMM nonlop
"t03.abi", # Use GPU OpenMP
# "HPC" Perf tests, no forces nor stress
"t04.abi", # Basic test, reference
"t05.abi", # Use GPU OpenMP
# Smoke tests on silicium, with batching, no forces nor stress (based on libxc#51)
"t06.abi", # Basic test, reference
"t07.abi", # Use GEMM nonlop
"t08.abi", # Use GPU OpenMP
"t09.abi", # Use GPU OpenMP (NC), no batching
"t10.abi", # Use GPU OpenMP (NC)
# Smoke tests on silicium, with stress (based on libxc#51)
"t11.abi", # Basic test, reference, PAW
"t12.abi", # With batching
"t13.abi", # With batching, Use GEMM nonlop
"t14.abi", # With batching, Use GPU OpenMP
# Smoke tests on methane, with stress&forces (based on libxc#71)
"t15.abi", # Basic test, reference, PAW
"t16.abi", # With batching
"t17.abi", # With batching, Use GEMM nonlop
"t18.abi", # With batching, Use GPU OpenMP
# Weird single-band test with istwfk==2 on Fock (based on libxc#73)
"t19.abi", # With batching, Use GPU OpenMP
]
