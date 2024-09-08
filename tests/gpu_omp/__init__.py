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
"t01.abi",
"t02.abi",
"t03.abi", # test CHEBFI (istwfk==1, npband==1, paral_kgb==1)
"t04.abi", # test CHEBFI (istwfk==2, npband==1, paral_kgb==1)
"t05.abi", # test LOBPCG (istwfk==1, npband==1, paral_kgb==1)
"t06.abi", # test LOBPCG (istwfk==2, npband==1, paral_kgb==1)
"t07.abi", # test CHEBFI (istwfk==1, npband==2, paral_kgb==1)
"t08.abi", # test CHEBFI (istwfk==2, npband==2, paral_kgb==1)
"t09.abi", # test LOBPCG (istwfk==1, npband==2, paral_kgb==1)
"t10.abi", # test LOBPCG (istwfk==2, npband==2, paral_kgb==1)
"t11.abi", # test CHEBFI (istwfk==1, npband==1, paral_kgb==0)
"t12.abi", # test CHEBFI (istwfk==2, npband==1, paral_kgb==0)
"t13.abi", # test LOBPCG (istwfk==1, npband==1, paral_kgb==0)
"t14.abi", # test LOBPCG (istwfk==2, npband==1, paral_kgb==0)
#Norm-conserving
"t15.abi", # test CHEBFI (istwfk==1, npband==1, paral_kgb==0)
"t16.abi", # test LOBPCG (istwfk==1, npband==1, paral_kgb==0)
#old lobppcg/chebfi (wfoptalg=={4,14}
"t17.abi", # test LOBPCG (istwfk==1, npband==1, paral_kgb==0)
"t18.abi", # test LOBPCG (istwfk==2, npband==1, paral_kgb==0)
"t19.abi", # test LOBPCG (istwfk==1, npband==1, paral_kgb==0)
"t20.abi", # test LOBPCG (istwfk==2, npband==1, paral_kgb==0)
]
