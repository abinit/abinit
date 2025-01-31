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
"t02.abi", # PAW, GS+DFPT, only DDK           (AlP, nspinor==2, nspden==1)
"t03.abi", # PAW, GS, old LOBPCG, (NiO, nspinor==1, paral_kgb==0)
"t04.abi", # PAW, GS, old LOBPCG, (NiO, nspinor==2, nspden==4, paral_kgb==0)
"t05.abi", # PAW, GS, LOBPCG2,    (NiO, nspinor==2, nspden==4, paral_kgb==1)
"t06.abi", # PAW, GS, CGWF,       (Fe,  nspinor==2, nspden==4, paral_kgb==0), with fixed_occ
"t07.abi", # NC,  GS, CGWF, mGGA  (Si,  nspinor==2, nspden==4)
"t08.abi", # NC,  GS, CGWF,       (H,   nspinor==2, nspden==4)
"t09.abi", # PAW, GS, CGWF,       (Fe, nspinor==2, nspden==4, paral_kgb==0)
"t10.abi", # NC,  GS+DFPT, CGWF, phon+elfd+effmas (CaO, nspinor==1)
"t11.abi", # NC,  GS+DFPT, CGWF, phon+elfd+strs   (AlAs, nspinor==1)
"t12.abi", # PAW, GS, LOBPCG2,    (NiO, nspinor==2, nspden==4, paral_kgb==1, npband==2)
]
