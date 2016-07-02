program detect_fc

 implicit none

! Compaq Fortran compiler
#if defined __LANGUAGE_FORTRAN_90__
#if defined __alpha
 character(len=*),parameter :: arch = "ALPHA"
#else
 character(len=*),parameter :: arch = "UNKNOWN"
#endif
#if defined __linux__
 character(len=*),parameter :: system = "LINUX"
#elif defined __osf__
 character(len=*),parameter :: system = "OSF"
#else
 character(len=*),parameter :: system = "UNKNOWN"
#endif
 write (*,'("COMPAQ_FORTRAN UNKNOWN ",A,1X,A)') arch,system

! GNU Fortran compiler
#elif defined __GNU_FORTRAN
 write (*,'(A)') "GFORTRAN"

! Intel Fortran compiler
#elif defined __INTEL_COMPILER
 integer,parameter :: version = __INTEL_COMPILER
 integer,parameter :: build   = __INTEL_COMPILER_BUILD_DATE
#if defined __i386__
 character(len=*),parameter :: arch = "IA32"
#elif defined __ia64__
 character(len=*),parameter :: arch = "IA64"
#else
 character(len=*),parameter :: arch = "UNKNOWN"
#endif
#if defined __linux__
 character(len=*),parameter :: system = "LINUX"
#else
 character(len=*),parameter :: system = "UNKNOWN"
#endif
 write (*,'("IFORT ",I1,".",I1,1X,I8,1X,A,1X,A)') version/100, &
  & mod(version,100)/10,build,arch,system
#endif

end program detect_fc
