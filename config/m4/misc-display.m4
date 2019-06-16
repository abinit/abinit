# -*- Autoconf -*-
#
# Copyright (C) 2005-2019 ABINIT Group (Yann Pouillon)
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
if test "${enable_openmp}" = "yes"; then
  tmp_omp_collapse="${abi_omp_has_collapse}"
else
  tmp_omp_collapse="ignored"
fi

# Set library-related status reports
#
# Timer
if test "${with_timer_libs}" = ""; then
  tmp_rep_timer_libs="auto-detected"
else
  tmp_rep_timer_libs="user-defined"
fi
test "${lib_timer_flavor}" = "none" -o "${lib_timer_flavor}" = "abinit" && \
  tmp_rep_timer_libs="ignored"

# Linalg
if test "${with_linalg_libs}" = ""; then
  tmp_rep_linalg_libs="auto-detected"
else
  tmp_rep_linalg_libs="user-defined"
fi
test "${lib_linalg_flavor}" = "netlib-fallback" && \
  tmp_rep_linalg_libs="internal"
test "${lib_linalg_flavor}" = "none" -o \
     "${lib_linalg_flavor}" = "netlib-fallback" && \
  tmp_rep_linalg_libs="ignored"

# Algo
if test "${with_algo_libs}" = ""; then
  tmp_rep_algo_libs="auto-detected"
else
  tmp_rep_algo_libs="user-defined"
fi
test "${lib_algo_flavor}" = "none" && \
  tmp_rep_algo_libs="ignored"

# FFT
if test "${with_fft_libs}" = ""; then
  tmp_rep_fft_libs="auto-detected"
else
  tmp_rep_fft_libs="user-defined"
fi
test "${lib_fft_flavor}" = "none" -o "${lib_fft_flavor}" = "abinit" && \
  tmp_rep_fft_libs="ignored"

# Display values of important configure options and ending message
cat <<EOF

Summary of important options:

  * C compiler      : ${abi_cc_vendor} version ${abi_cc_version}
  * Fortran compiler: ${abi_fc_vendor} version ${abi_fc_version}
  * architecture    : ${abi_cpu_vendor} ${abi_cpu_model} (${abi_cpu_bits} bits)

  * debugging       : ${enable_debug}
  * optimizations   : ${enable_optim}

  * OpenMP enabled  : ${enable_openmp} (collapse: ${tmp_omp_collapse})
  * MPI    enabled  : ${enable_mpi}
  * MPI-IO enabled  : ${enable_mpi_io}
  * GPU    enabled  : ${enable_gpu} (flavor: ${lib_gpu_flavor})
  * XML    enabled  : ${enable_xml}

  * TRIO   flavor = ${lib_trio_flavor}
  * TIMER  flavor = ${lib_timer_flavor} (libs: ${tmp_rep_timer_libs})
  * LINALG flavor = ${lib_linalg_flavor} (libs: ${tmp_rep_linalg_libs})
  * ALGO   flavor = ${lib_algo_flavor} (libs: ${tmp_rep_algo_libs})
  * FFT    flavor = ${lib_fft_flavor} (libs: ${tmp_rep_fft_libs})
  * DFT    flavor = ${lib_dft_flavor}

Configuration complete.
You may now type "make" to build ABINIT.
(or, on a SMP machine, "make mj4", or "make multi multi_nprocs=<n>")

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
