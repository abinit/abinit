# -*- Autoconf -*-
#
# Copyright (C) 2012-2022 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# OpenMP support
#



# ABI_OMP_CHECK_COLLAPSE()
# ------------------------
#
# Check whether OpenMP has a working implementation of COLLAPSE.
#
AC_DEFUN([ABI_OMP_CHECK_COLLAPSE],[
  dnl Init
  abi_omp_has_collapse="unknown"

  dnl Check whether OpenMP's COLLAPSE is working
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether OpenMP COLLAPSE works])
  AC_RUN_IFELSE([AC_LANG_PROGRAM([],
    [[
      integer :: alfa,i,levels
      levels = 1
      open(unit=10,file="conftest.collapse")
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(alfa,i) DEFAULT(shared)
       do alfa=1,1
         do i=1,levels
           write(10,'(I1)') i*alfa
         end do
       end do
!$OMP END PARALLEL DO
       close(unit=10)
    ]])], [abi_omp_has_collapse="yes"; abi_omp_collapse_result=`cat conftest.collapse 2>/dev/null`; rm -f conftest.collapse], [abi_omp_has_collapse="no"; rm -f conftest.collapse])
  test "${abi_omp_collapse_result}" = "1" || abi_omp_has_collapse="no"
  AC_MSG_RESULT([${abi_omp_has_collapse}])
  AC_LANG_POP([Fortran])

  dnl Propagate result
  if test "${abi_omp_has_collapse}" = "yes"; then
    AC_DEFINE([HAVE_OMP_COLLAPSE],1,[Set to 1 if OpenMP has a working implementation of COLLAPSE.])
  fi
]) # ABI_OMP_CHECK_COLLAPSE

# ABI_OMP_CHECK_GPU_OFFLOAD()
# ------------------------
#
# Check whether OpenMP has a working implementation of TARGET offload directives.
#
AC_DEFUN([ABI_OMP_CHECK_GPU_OFFLOAD],[
  dnl Init
  abi_omp_has_gpu_offload="unknown"

  dnl Check whether OpenMP's TARGET offload directives are recognized
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether OpenMP TARGET offload directives are recognized])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
    [[
      integer :: alfa,i,levels,arr(7)
      levels = 7
      arr = 0

      ! Not an actual GPU offload test, as users may not compile on host with GPU.
      ! Only compiler capabilities are checked.
      if (levels == 8) then
!$OMP TARGET ENTER EXIT MAP(to:arr)
!$OMP TARGET PARALLEL DO COLLAPSE(2) PRIVATE(alfa,i) MAP(to:arr)
      do alfa=1,1
        do i=1,levels
          arr(i)=i
        end do
      end do
!$OMP END TARGET PARALLEL DO

!$OMP TARGET DATA USE_DEVICE_ADDR(arr)
      alfa=1
!$OMP END TARGET DATA
!$OMP TARGET EXIT DATA MAP(from:arr)
      end if

    ]])], [abi_omp_has_gpu_offload="yes"], [abi_omp_has_gpu_offload="no"])
  AC_MSG_RESULT([${abi_omp_has_gpu_offload}])
  AC_LANG_POP([Fortran])

  dnl Propagate result
  if test "${abi_omp_has_gpu_offload}" != "yes"; then
    AC_MSG_ERROR([OpenMP GPU offload is enabled but is not supported by your compiler. Perhaps do you have flags missing ?])
  fi
]) # ABI_OMP_CHECK_GPU_OFFLOAD
