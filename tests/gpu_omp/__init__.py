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
#DFPT
"t21.abi", # test respfn for phonons and electric field, no qpt
"t22.abi", # prep for t23, t24
"t23.abi", # test respfn for phonons and electric field, with qpt (GEMM nonlop CPU)
"t24.abi", # test respfn for phonons and electric field, with qpt
"t25.abi", # test respfn for stresses, no qpt
#Fock
"t26.abi", # Silicium, norm-conserving
"t27.abi", # Silicium, PAW
"t28.abi", # Methane, PAW
"t29.abi", # Helium, single-band, istwfk==2
#nspinor==2 (and nspden==4 sometimes)
"t30.abi", # PAW, GS, using CHEBFI and LOBPCG (AlP, nspinor==2, nspden==1)
"t31.abi", # PAW, GS+DFPT, only DDK           (AlP, nspinor==2, nspden==1)
"t32.abi", # PAW, GS, old LOBPCG, (NiO, nspinor==1, paral_kgb==0)
"t33.abi", # PAW, GS, old LOBPCG, (NiO, nspinor==2, nspden==4, paral_kgb==0)
"t34.abi", # PAW, GS, LOBPCG2,    (NiO, nspinor==2, nspden==4, paral_kgb==1)
"t35.abi", # PAW, GS, CGWF,       (Fe,  nspinor==2, nspden==4, paral_kgb==0), with fixed_occ
"t36.abi", # NC,  GS, CGWF, mGGA  (Si,  nspinor==2, nspden==4)
"t37.abi", # NC,  GS, CGWF,       (H,   nspinor==2, nspden==4)
"t38.abi", # PAW, GS, CGWF,       (Fe, nspinor==2, nspden==4, paral_kgb==0)
"t39.abi", # NC,  GS+DFPT, CGWF, phon+elfd+effmas (CaO, nspinor==1)
"t40.abi", # NC,  GS+DFPT, CGWF, phon+elfd+strs   (AlAs, nspinor==1)
"t41.abi", # PAW, GS, LOBPCG2,    (NiO, nspinor==2, nspden==4, paral_kgb==1, npband==2)
#GS with distributed/sliced GEMM nonlop
"t42.abi", # test CHEBFI (istwfk==1, npband==2, paral_kgb==1, gpu_nl_splitsize=2, gpu_nl_distrib=1)
"t43.abi", # test CHEBFI (istwfk==2, npband==2, paral_kgb==1, gpu_nl_splitsize=2, gpu_nl_distrib=1)
"t44.abi", # test CHEBFI (istwfk==1, npband==1, paral_kgb==1, gpu_nl_splitsize=4)
"t45.abi", # test CHEBFI (istwfk==2, npband==1, paral_kgb==1, gpu_nl_splitsize=4)
#DMFT
"t46.abi", # Ni+O, nsppol==2, nspden==2
"t47.abi", # Gd, nsppol==2, nspden==2
"t48.abi", # Gd, nsppol==1, nspden==4, nspinor==2
"t49.abi", # Ce, GS LDA, dtset#1
"t50.abi", # Ce, DMFT CTQMC, dtset#2
#Other
"t51.abi", # test CHEBFI (istwfk==1, npband==1, paral_kgb==1) with gpu_thread_limit set to 64
]
