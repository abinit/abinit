# -*- Autoconf -*-
#
# Copyright (C) 2012-2025 ABINIT Group (Yann Pouillon, MTorrent)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# OpenMP support for ABINIT
#


# ABI_OMP_CHECK_COLLAPSE()
# ------------------------
#
# Check whether OpenMP has a working Fortran implementation of COLLAPSE.
#
AC_DEFUN([_ABI_OMP_CHECK_COLLAPSE],[
  dnl Init
  abi_omp_has_collapse="unknown"

  dnl Check whether OpenMP's COLLAPSE is working
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether OpenMP COLLAPSE works with Fortran])
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


                    # ------------------------------------ #


# ABI_OMP_CHECK_GPU_OFFLOAD()
# ------------------------
#
# Check whether OpenMP has a working Fortran implementation of TARGET offload directives.
#
AC_DEFUN([_ABI_OMP_CHECK_GPU_OFFLOAD],[
  dnl Init
  abi_omp_has_gpu_offload="unknown"
  abi_omp_has_gpu_get_mapped_ptr="unknown"
  abi_omp_has_gpu_offload_datastructure="unknown"

  dnl Check whether OpenMP's TARGET offload directives are recognized
  AC_LANG_PUSH([Fortran])
  AC_MSG_CHECKING([whether OpenMP TARGET offload directives are recognized in Fortran])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
    [[
      use omp_lib
      use, intrinsic :: iso_c_binding, only : c_ptr
      integer :: alfa,i,levels,arr(7),device_id,rc,num
      logical :: is_bool
      type(c_ptr) :: ptr

      levels = 7
      arr = 0
      is_bool = omp_is_initial_device()
      device_id = omp_get_initial_device()
      device_id = omp_get_default_device()
      rc = omp_target_is_present(ptr, device_id)
      call omp_set_default_device(device_id)
      num = omp_get_num_devices()

      ! Not an actual GPU offload test, as users may not compile on host with GPU.
      ! Only compiler capabilities are checked.
      if (levels == 8) then
!$OMP TARGET ENTER DATA MAP(to:arr)
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

  if test "${abi_omp_has_gpu_offload}" != "yes"; then
    AC_MSG_ERROR([OpenMP GPU offload is enabled but is not supported by your compiler. Perhaps do you have flags missing ?])
  fi

  dnl Check whether compiler implements omp_get_mapped_ptr
  if test "${abi_omp_has_gpu_offload}" == "yes"; then
    AC_MSG_CHECKING([whether compiler implements omp_get_mapped_ptr])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
      [[
        use omp_lib
        use, intrinsic :: iso_c_binding, only : c_ptr
        integer :: device_id
        type(c_ptr) :: ptr,ptr2

        device_id = omp_get_default_device()

        ptr2 = omp_get_mapped_ptr(ptr, device_id)

      ]])], [abi_omp_has_gpu_get_mapped_ptr="yes"], [abi_omp_has_gpu_get_mapped_ptr="no"])
    AC_MSG_RESULT([${abi_omp_has_gpu_get_mapped_ptr}])
    if test "${abi_omp_has_gpu_get_mapped_ptr}" == "yes"; then
      AC_DEFINE([HAVE_OPENMP_GET_MAPPED_PTR], 1,
        [Define to 1 if compiler support omp_get_mapped_ptr Fortran routine from OpenMP 5.1.])
    fi
  fi

  dnl Check whether compiler correctly handles structured types in OMP TARGET directives
  if test "${abi_omp_has_gpu_offload}" == "yes"; then
    AC_MSG_CHECKING([whether compiler correctly handles structured types in OMP TARGET directives])

    dnl At present (2023) only NVHPC implements correctly structured types in OMP TARGET
    if test "${abi_fc_vendor}" == "nvhpc"; then
      abi_omp_has_gpu_offload_datastructure="yes"
    else
      abi_omp_has_gpu_offload_datastructure="no"
    fi

    AC_MSG_RESULT([${abi_omp_has_gpu_offload_datastructure}])
    if test "${abi_omp_has_gpu_offload_datastructure}" == "yes"; then
      AC_DEFINE([HAVE_OPENMP_OFFLOAD_DATASTRUCTURE], 1,
        [Define to 1 if compiler correctly handles structured types in OMP TARGET directives.])
    fi
  fi

  AC_LANG_POP([Fortran])
]) # ABI_OMP_CHECK_GPU_OFFLOAD


                    # ------------------------------------ #


# ABI_OMP_SET_FCFLAGS()
# ------------------------
#
# Adjust openMP Fortran FLAGS according to configure options
#
AC_DEFUN([_ABI_OMP_SET_FCFLAGS],[

# Empty flags if openMP is not used
if test "${abi_openmp_enable}" != "yes"; then
  FCFLAGS_OPENMP=""
  FCFLAGS_OPENMP_OFFLOAD=""
fi

# Handling OpenMP offload
if test "${abi_openmp_offload_enable}" = "yes"; then

  # Fail if GPU support is not enabled
  if test "${abi_gpu_enable}" = "no"; then
    AC_MSG_ERROR([OpenMP offload is enabled with no GPU support (i.e. CUDA, HIP, ...), please check your config.])
  fi

  # Enable OpenMP if OpenMP offload is requested, to remain consistent
  abi_openmp_enable="yes"

  # Manage ARCH variable
  if test "${GPU_ARCH}" = ""; then
    AC_MSG_ERROR([GPU target architecture is not set, please provide GPU_ARCH variable])
    GPU_ARCH=80
  fi
  # Perform pattern replacement in OpenMP offload flags with requested GPU arch
  FCFLAGS_OPENMP_OFFLOAD=`echo "${FCFLAGS_OPENMP_OFFLOAD}" | sed "s/__GPU_ARCH__/$GPU_ARCH/"`

  #FIXME With LLVM 16 embedded with ROCm 5.6.0, issues occurs if OpenMP offload is enabled everywhere
  #Therefore, we enable it only if subfolders where it's needed and link time.
  if test "${abi_fc_vendor}" == "llvm" -a "${abi_gpu_flavor}" == "hip-double"; then
    amd_openmp_flags="${FCFLAGS_OPENMP_OFFLOAD}"
    FC_LDFLAGS_EXTRA="${FCFLAGS_OPENMP_OFFLOAD} ${FC_LDFLAGS_EXTRA}"
    FCFLAGS_OPENMP_OFFLOAD=""
  fi

  #FIXME With Cray CPE 23.12, it seems that "-lcraymp" isn't always added to LDFLAGS, so we add it here
  if test "${abi_fc_vendor}" == "cray" -a "${abi_gpu_flavor}" == "hip-double"; then
    FC_LDFLAGS_EXTRA="-lcraymp ${FC_LDFLAGS_EXTRA}"
  fi

else
  FCFLAGS_OPENMP_OFFLOAD=""
fi

]) # ABI_OMP_SET_FCFLAGS


                    # ------------------------------------ #


# ABI_OPENMP_INIT()
# --------------
#
# Looks for an implementation of openMP.
# Note 1: this is a convenience feature, purely for comfort.
# Note 2: it should be run as early as possible.
#
AC_DEFUN([ABI_OPENMP_INIT], [

  # Use existing OPENMP env variables (Fortran only)
  OPENMP_FCFLAGS=""
  OPENMP_CFLAGS=""
  OPENMP_CXXFLAGS=""
  OPENMP_LDFLAGS=""
  OPENMP_LIBS=""
  test ! -z "${FCFLAGS_OPENMP}" && OPENMP_FCFLAGS="${OPENMP_FCFLAGS} ${FCFLAGS_OPENMP}"
  test ! -z "${FCFLAGS_OPENMP_OFFLOAD}" && OPENMP_FCFLAGS="${OPENMP_FCFLAGS} ${FCFLAGS_OPENMP_OFFLOAD}"

  # Delegate most of the init stage to Steredeg
  SD_OPENMP_INIT([auto optional warn no-cc no-cxx])

  # Init ABINIT openMP variables
  abi_openmp_enable="${sd_openmp_enable}"
  abi_openmp_init="${sd_openmp_init}"

  # Init ABINIT openMP build flags
  abi_openmp_cppflags=""
  abi_openmp_cflags=""
  abi_openmp_cxxflags=""
  abi_openmp_fcflags=""
  abi_openmp_ldflags=""
  abi_openmp_libs=""

  # Init ABINIT openMP build parameters
  if test "${abi_openmp_enable}" = "yes" -o "${abi_openmp_enable}" = "auto"; then
    abi_openmp_cppflags="${sd_openmp_cppflags}"
    abi_openmp_ldflags="${sd_openmp_ldflags}"
    abi_openmp_libs="${sd_openmp_libs}"
    if test "${sd_openmp_cc_ok}" = "yes"; then
      abi_openmp_cflags="${sd_openmp_cflags}"
    fi
    if test "${sd_openmp_cxx_ok}" = "yes"; then
      abi_openmp_cxxflags="${sd_openmp_cxxflags}"
    fi
    if test "${sd_openmp_fc_ok}" = "yes"; then
      abi_openmp_fcflags="${sd_openmp_fcflags}"
    fi
  else
    if test "${abi_openmp_init}" != "def"; then
      AC_MSG_NOTICE([openMP support disabled from command-line])
    fi
    abi_openmp_enable="no"
    abi_openmp_init="cmd"
    abi_openmp_cppflags=""
    abi_openmp_cflags=""
    abi_openmp_cxxflags=""
    abi_openmp_fcflags=""
    abi_openmp_ldflags=""
    abi_openmp_libs=""
  fi # abi_openmp_enable

  # Enable substitution
  AC_SUBST(abi_openmp_enable)
  AC_SUBST(abi_openmp_flavor)
  AC_SUBST(abi_openmp_cppflags)
  AC_SUBST(abi_openmp_cflags)
  AC_SUBST(abi_openmp_cxxflags)
  AC_SUBST(abi_openmp_fcflags)
  AC_SUBST(abi_openmp_ldflags)
  AC_SUBST(abi_openmp_level)
  AC_SUBST(abi_openmp_incs)
  AC_SUBST(abi_openmp_libs)
]) # ABI_OPENMP_INIT


                    # ------------------------------------ #

# ABI_OPENMP_DETECT()
# ----------------
#
# Tries first to determine whether the openMP implementation is usable,
# then takes appropriate actions.
#
AC_DEFUN([ABI_OPENMP_DETECT], [
  # Delegate the actual detection to Steredeg
  SD_OPENMP_DETECT
  AC_MSG_NOTICE([OpenMP support is enabled in Fortran source code only])

  abi_openmp_ok="${sd_openmp_ok}"

  if test "${sd_openmp_enable}" = "yes"; then

    if test "${sd_openmp_ok}" = "yes"; then

      # Force abi_openmp_enable to "yes", for clarity and to avoid having to
      # further test "auto"
      abi_openmp_enable="yes"

      # Check whether OpenMP has a working implementation of COLLAPSE
      _ABI_OMP_CHECK_COLLAPSE

      # Check if openMP offload is available
      abi_omp_has_gpu_offload="no"
      if test "${abi_openmp_offload_enable}" = "yes"; then
        _ABI_OMP_CHECK_GPU_OFFLOAD
      fi

    fi # sd_openmp_ok

  fi # sd_mpi_enable
]) # ABI_OPENMP_DETECT
