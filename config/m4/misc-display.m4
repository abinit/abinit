# -*- Autoconf -*-
#
# Copyright (C) 2005-2022 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Miscellaneous macros
#



# ABI_MSG_END()
# -------------
#
# Prints a message at the end of the configure process.
#
AC_DEFUN([ABI_MSG_END],
[[
# Set OpenMP status reports
#
if test "${abi_openmp_enable}" = "yes"; then
  tmp_omp_collapse="${abi_omp_has_collapse}"
  tmp_omp_gpu_offload="${abi_omp_has_gpu_offload}"
else
  tmp_omp_collapse="ignored"
  tmp_omp_gpu_offload="ignored"
fi

# Set library-related status reports

# Linalg
#if test "${with_linalg_libs}" = ""; then
if test "${sd_linalg_flavor}" != "${sd_linalg_flavor_req}"; then
  tmp_rep_linalg_libs="auto-detected"
else
  tmp_rep_linalg_libs="user-defined"
fi
test "${sd_linalg_flavor}" = "netlib-fallback" && \
  tmp_rep_linalg_libs="internal"
test "${sd_linalg_flavor}" = "none" -o \
     "${sd_linalg_flavor}" = "netlib-fallback" && \
  tmp_rep_linalg_libs="ignored"

# FFT
if test "${sd_fft_init}" = "def"; then
  tmp_rep_fft_libs="auto-detected"
else
  tmp_rep_fft_libs="user-defined"
fi
test "${sd_fft_flavor}" = "none" && \
  tmp_rep_fft_libs="ignored"

# Display values of important configure options and ending message
cat <<EOF

Core build parameters
---------------------

  * C compiler        : ${abi_cc_vendor} version ${abi_cc_version}
  * Fortran compiler  : ${abi_fc_vendor} version ${abi_fc_version}
  * architecture      : ${abi_cpu_vendor} ${abi_cpu_model} (${abi_cpu_bits} bits)
  * debugging         : ${abi_debug_flavor}
  * optimizations     : ${abi_optim_flavor}

  * OpenMP enabled    : ${abi_openmp_enable} (collapse: ${tmp_omp_collapse}; GPU offload: ${tmp_omp_gpu_offload})
  * MPI    enabled    : ${abi_mpi_enable} (flavor: ${abi_mpi_flavor})
  * MPI    in-place   : ${abi_mpi_inplace_enable}
  * MPI-IO enabled    : ${abi_mpi_io_enable}
  * GPU    enabled    : ${abi_gpu_enable} (flavor: ${abi_gpu_flavor})

  * LibXML2 enabled   : ${abi_libxml2_enable}
  * LibPSML enabled   : ${sd_libpsml_enable}
  * XMLF90  enabled   : ${sd_xmlf90_enable}
  * HDF5 enabled      : ${sd_hdf5_enable} (MPI support: ${sd_hdf5_mpi_ok})
  * NetCDF enabled    : ${sd_netcdf_enable} (MPI support: ${sd_netcdf_mpi_ok})
  * NetCDF-F enabled  : ${sd_netcdf_fortran_enable} (MPI support: ${sd_netcdf_fortran_mpi_ok})

  * FFT flavor        : ${sd_fft_flavor} (libs: ${tmp_rep_fft_libs})
  * LINALG flavor     : ${sd_linalg_flavor} (libs: ${tmp_rep_linalg_libs})
  * SCALAPACK enabled : ${sd_linalg_has_scalapack}
  * ELPA enabled      : ${sd_linalg_has_elpa}
  * MAGMA enabled     : ${sd_linalg_has_magma} (magma version >= 1.5 ? ${sd_linalg_has_magma_15})

  * FCFLAGS           : ${FCFLAGS}
  * NVCC_CFLAGS       : ${NVCC_CFLAGS}
  * CPATH             : ${CPATH}

  * Build workflow    : ${abi_build_steps}

${abi_opt_deprecated_count} deprecated options have been used:${abi_opt_deprecated_used}.

Configuration complete.
You may now type "make" to build Abinit.
(or "make -j<n>", where <n> is the number of available processors)

EOF
]]) # ABI_MSG_END



# ABI_MSG_FC_BUGGY(FC_TYPE)
# -------------------------
#
# Prints a message explaining why a compiler has been wrapped, or giving
# advice for its use.
#
AC_DEFUN([ABI_MSG_FC_BUGGY],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  case "$1" in

    absoft)
      ABI_MSG_NOTICE([fc-absoft],[About the ABSoft Fortran compiler])
      ;;

    ibm)
      ABI_MSG_NOTICE([fc-ibm],[About the IBM XL Fortran compiler])
      ;;

    intel)
      ABI_MSG_NOTICE([fc-intel],[About the Intel Fortran compiler])
      ;;

  esac
]) # ABI_MSG_FC_BUGGY


dnl ABI_MSG_NOTICE(FILE, TITLE)
dnl ---------------------------
dnl
dnl Print a framed message to attract users' attention to something.
dnl
AC_DEFUN([ABI_MSG_NOTICE],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl Init
  abi_msg_file="${abinit_srcdir}/config/messages/$1.msg"
  abi_msg_title="$2"
  test "${abi_msg_title}" = "" && abi_msg_title="IMPORTANT NOTE"

  dnl Format title
  abi_msg_spacer="                                                            "
  abi_msg_tmp1=`echo "${abi_msg_title}" | sed -e 's/./ /g'`
  abi_msg_tmp2=`echo "${abi_msg_tmp1}" | grep "${abi_msg_spacer}"`
  abi_msg_spacer=`echo "${abi_msg_spacer}" | sed -e "s/${abi_msg_tmp1}//"`
  test "${abi_msg_tmp2}" = "" || abi_msg_spacer=""
  abi_msg_title="${abi_msg_title}${abi_msg_spacer}"

  if test -s "${abi_msg_file}"; then

  dnl Print header
  echo ""
  echo "        +--------------------------------------------------------------+"
  echo "        | ${abi_msg_title} |"
  echo "        +--------------------------------------------------------------+"

  dnl Format and write message
  while read abi_msg_line; do
    abi_msg_line=`eval echo ${abi_msg_line}`
    abi_msg_spacer="                                                            "
    abi_msg_tmp1=`echo "${abi_msg_line}" | sed -e 's/./ /g'`
    abi_msg_tmp2=`echo "${abi_msg_tmp1}" | grep "${abi_msg_spacer}"`
    test "${abi_msg_tmp1}" = "" || \
      abi_msg_spacer=`echo "${abi_msg_spacer}" | sed -e "s/${abi_msg_tmp1}//"`
    test "${abi_msg_tmp2}" = "" || abi_msg_spacer=""
    echo "        | ${abi_msg_line}${abi_msg_spacer} |"
  done <"${abi_msg_file}"

  dnl Print footer
  echo "        +--------------------------------------------------------------+"
  echo ""

  else
    AC_MSG_WARN([message file ${abi_msg_file} not found])
  fi
]) dnl ABI_MSG_NOTICE


dnl ABI_MSG_NOTICE_L(FILE, TITLE)
dnl ---------------------------
dnl
dnl Print a framed message to attract users' attention to something.
dnl
AC_DEFUN([ABI_MSG_NOTICE_L],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl Init
  abi_msg_file="${abinit_srcdir}/config/messages/$1.msg"
  abi_msg_title="$2"

  TPUT_BOLD="\$(tput bold)"
  TPUT_OFF="\$(tput sgr0)"

  let spacer=64
  let tput_spacer=0
  while test "`echo $abi_msg_title | grep __`" ; do
     abi_msg_title=$(echo "$abi_msg_title" | sed "s/__/$TPUT_BOLD/;s/__/$TPUT_OFF/")
     let tput_spacer=tput_spacer+10
  done
  abi_msg_title=`eval echo ${abi_msg_title}`

  dnl Format title
  let titlel=${#abi_msg_title}
  let Fill=spacer-titlel+tput_spacer
  abi_msg_title=`printf "${abi_msg_title}";printf ' %.0s' $(seq 1 $Fill)`

  if test -s "${abi_msg_file}"; then

    dnl Print header
    echo ""
    echo "  +------------------------------------------------------------------+"
    echo "  | ${abi_msg_title} |"
    echo "  +------------------------------------------------------------------+"

    dnl Format and write message

    while read abi_msg_line; do

      let tput_spacer=0
      while test "`echo $abi_msg_line | grep __`" ; do
          abi_msg_line=$(echo "$abi_msg_line" | sed "s/__/$TPUT_BOLD/;s/__/$TPUT_OFF/")
          let tput_spacer=tput_spacer+10
      done
      abi_msg_line=`eval echo ${abi_msg_line}`

      dnl Format line
      let linel=${#abi_msg_line}
      let Fill=spacer-linel+tput_spacer
      dnl test "$linel" -gt "64" || echo "too long... : $linel, Fill = $Fill"
      abi_msg_line=`printf "${abi_msg_line}";printf ' %.0s' $(seq 1 $Fill)`

      echo "  | ${abi_msg_line} |"
    done <"${abi_msg_file}"

    dnl Print footer
    echo "  +------------------------------------------------------------------+"
    echo ""

  else
    AC_MSG_WARN([message file ${abi_msg_file} not found])
  fi
]) dnl ABI_MSG_NOTICE_L



dnl ABI_MSG_NOTICE_S(TITLE, SOLUTION)
dnl ---------------------------------
dnl
dnl Print a framed message to attract users' attention to something.
dnl
AC_DEFUN([ABI_MSG_NOTICE_S],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl Init
  abi_msg_title="$1"

  TPUT_BOLD="\$(tput bold)"
  TPUT_OFF="\$(tput sgr0)"

  test "${abi_msg_title}" = "" && abi_msg_title="\$(tput bold)IMPORTANT NOTE\$(tput sgr0)"

  let spacer=64
  let tput_spacer=0
  while test "`echo $abi_msg_title | grep __`" ; do
     abi_msg_title=$(echo "$abi_msg_title" | sed "s/__/$TPUT_BOLD/;s/__/$TPUT_OFF/")
     let tput_spacer=tput_spacer+10
  done
  abi_msg_title=`eval echo ${abi_msg_title}`

  dnl Format title
  titlel=${#abi_msg_title}
  let Fill=spacer-titlel+tput_spacer
  abi_msg_title=`printf "${abi_msg_title}";printf ' %.0s' $(seq 1 $Fill)`

  dnl Print header
  echo ""
  echo "  +------------------------------------------------------------------+"
  echo "  | ${abi_msg_title} |"
  echo "  +------------------------------------------------------------------+"

  dnl Format solution if provided
  abi_msg_line="$2"
  if test "${abi_msg_line}" != ""; then
      let spacer=64
      let tput_spacer=0
      while test "`echo $abi_msg_line | grep __`" ; do
          abi_msg_line=$(echo "$abi_msg_line" | sed "s/__/$TPUT_BOLD/;s/__/$TPUT_OFF/")
          let tput_spacer=tput_spacer+10
      done
      abi_msg_line=`eval echo ${abi_msg_line}`

      dnl Format line
      linel=${#abi_msg_line}
      let Fill=spacer-linel+tput_spacer
      dnl test "$linel" -gt "64" || echo "too long... : $linel, Fill = $Fill"
      abi_msg_line=`printf "${abi_msg_line}";printf ' %.0s' $(seq 1 $Fill)`

      echo "  | ${abi_msg_line} |"
      echo "  +------------------------------------------------------------------+"
  fi

  echo ""

]) dnl ABI_MSG_NOTICE_S



dnl ABI_MSG_NOTICE_X(TITLE, SOLUTION)
dnl ---------------------------------
dnl
dnl Print a framed message to attract users' attention to something.
dnl
AC_DEFUN([ABI_MSG_NOTICE_X],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl Init
  abi_msg_title="$1"
  test "${abi_msg_title}" = "" && abi_msg_title="IMPORTANT NOTE"

  dnl Format title
  abi_msg_spacer="                                                                "
  abi_msg_tmp1=`echo "${abi_msg_title}" | sed -e 's/./ /g'`
  abi_msg_tmp2=`echo "${abi_msg_tmp1}" | grep "${abi_msg_spacer}"`
  abi_msg_spacer=`echo "${abi_msg_spacer}" | sed -e "s/${abi_msg_tmp1}//"`
  test "${abi_msg_tmp2}" = "" || abi_msg_spacer=""
  abi_msg_title="${abi_msg_title}${abi_msg_spacer}"

  dnl Print header
  echo ""
  echo "  +------------------------------------------------------------------+"
  echo "  | ${abi_msg_title} |"
  echo "  +------------------------------------------------------------------+"

  dnl Format solution if provided
  abi_msg_title="$2"
  if test "${abi_msg_title}" != ""; then
      abi_msg_spacer="                                                                "
      abi_msg_tmp1=`echo "${abi_msg_title}" | sed -e 's/./ /g'`
      abi_msg_tmp2=`echo "${abi_msg_tmp1}" | grep "${abi_msg_spacer}"`
      abi_msg_spacer=`echo "${abi_msg_spacer}" | sed -e "s/${abi_msg_tmp1}//"`
      test "${abi_msg_tmp2}" = "" || abi_msg_spacer=""
      abi_msg_title="${abi_msg_title}${abi_msg_spacer}"
      echo "  | ${abi_msg_title} |"
      echo "  +------------------------------------------------------------------+"
  fi

  echo ""

]) dnl ABI_MSG_NOTICE_X



# ABI_MSG_SECTION(TITLE)
# ----------------------
#
# Prints a nice title for each section.
#
AC_DEFUN([ABI_MSG_SECTION],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  abi_sec_title="$1"

  dnl Calculate amount of space chars needed for pretty-printing
  abi_sec_spaces="                                                                      "
  abi_sec_tmp="${abi_sec_title}"
  while test "${abi_sec_tmp}" != ""; do
    abi_sec_spaces=`echo "${abi_sec_spaces}" | sed -e 's/^.//'`
    abi_sec_tmp=`echo "${abi_sec_tmp}" | sed -e 's/^.//'`
  done

  echo ""
  echo " =============================================================================="
  echo " === ${abi_sec_title}${abi_sec_spaces} ==="
  echo " =============================================================================="
  echo ""
]) # ABI_MSG_SECTION
