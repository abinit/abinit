# -*- Autoconf -*-
#
# Copyright (C) 2006-2020 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# File I/O
#



# ABI_LOAD_OPTIONS()
# ------------------
#
# Looks for options to the configure script in predefined locations:
#
#   * system-wide : /etc/abinit/build/hostname.ac9
#   * per-user    : ~/.abinit/build/hostname.ac9
#   * local       : <current_directory>/hostname.ac9
#
# and eventually in a command-line-specified file. "hostname" is the
# name of the machine without the domain name. The last valid file is
# the one to be read.
#
# Loads an environment file as well if found (hostname.env).
#
# NOTE
#
#   The name choosen for the per-user the config file allows the peaceful
#   coexistence of several option sets on machines sharing the same home
#   directory (e.g. with NFS).
#
AC_DEFUN([ABI_LOAD_OPTIONS],[
  # Setup file names
  abi_hostname=`hostname | sed -e 's/\..*//'`
  abi_sys_options="/etc/abinit/build/${abi_hostname}.ac9"
  abi_per_options="${HOME}/.abinit/build/${abi_hostname}.ac9"
  abi_src_options="${abinit_srcdir}/${abi_hostname}.ac9"
  abi_loc_options="./${abi_hostname}.ac9"
  abi_cmd_options=`eval echo "${with_config_file}"`
  abi_cfg_enable=""
  abi_cfg_options=""
  abi_ac_distcheck=""

  # Check if the option is yes or no
  if test "${with_config_file}" = "no" -o "${with_config_file}" = "yes"; then
    abi_cfg_enable="${with_config_file}"
    abi_cmd_options=""
  else
    abi_cfg_enable="yes"
  fi

  # Some architectures require "./" for files in current directory
  if test "${abi_cmd_options}" != ""; then
    abi_cmd_has_path=`echo "${abi_cmd_options}" | grep '/'`
    if test "${abi_cmd_has_path}" = ""; then
      abi_cmd_options="./${abi_cmd_options}"
    fi
    unset abi_cmd_has_path
  fi

  # Select and read config file
  if test "${abi_cfg_enable}" = "yes"; then
    if test "${with_config_file}" != "" -a \
            ! -e "${with_config_file}"; then
      AC_MSG_ERROR([config file ${with_config_file} not found])
    fi

    for abi_options in "${abi_sys_options}" "${abi_per_options}" \
                       "${abi_src_options}" "${abi_loc_options}" \
                       "${abi_cmd_options}"; do
      if test -s "${abi_options}"; then
        abi_cfg_options="${abi_options}"
      fi
    done

    # Prevent infinite loops
    if grep "From configure.ac Autotools support for ABINIT" \
      "${abi_cfg_options}" >/dev/null 2>&1; then
      AC_MSG_ERROR([infinite loop detected - aborting!])
    fi

    # Source the file
    if test "${abi_cfg_options}" != ""; then
      AC_MSG_NOTICE([reading options from ${abi_cfg_options}])
      . "${abi_cfg_options}"
    else
      AC_MSG_NOTICE([not loading options (no config file available)])
    fi
  else
    AC_MSG_NOTICE([not loading options (disabled from command line)])
  fi

  # Propagate information to "make distcheck"
  abi_ac_distcheck=`${REALPATH} "${abi_cfg_options}"`
  if test "${abi_ac_distcheck}" != ""; then
    abi_ac_distcheck="--with-config-file=\"${abi_ac_distcheck}\""
  fi

  AC_SUBST(abi_ac_distcheck)
]) # ABI_LOAD_OPTIONS



# ABI_LOAD_DBGFLAGS(SUFFIX, COMPILER, VERSION, ARCHITECTURE)
# ----------------------------------------------------------
#
# Looks for per-directory debug flags for a specified compiler
# and a specified version, running on a specified architecture. Load
# them if found.
#
AC_DEFUN([ABI_LOAD_DBGFLAGS],[
  # Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl
  m4_if([$3], , [AC_FATAL([$0: missing argument 3])])dnl
  m4_if([$4], , [AC_FATAL([$0: missing argument 3])])dnl

  # Init
  abi_result=""
  abi_dbgflags_file=""

  # Explore all the possibilities
  for tmp_dbgflags_file in \
    "${ac_top_srcdir}/config/compilers/$2_$1/all/all.dbg" \
    "${ac_top_srcdir}/config/compilers/$2_$1/all/$4.dbg" \
    "${ac_top_srcdir}/config/compilers/$2_$1/$3/all.dbg" \
    "${ac_top_srcdir}/config/compilers/$2_$1/$3/$4.dbg"; do

    if test -s "${tmp_dbgflags_file}"; then
      abi_dbgflags_file="${tmp_dbgflags_file}"
      abi_result=`echo "${abi_dbgflags_file}" | \
        sed -e 's,.*compilers/,,; s,\.dbg$,,; s,^\([[^_]]*\)_\([[^/]]*\),\2: \1,'`
    fi
  done

  # Source the file
  #AC_MSG_NOTICE([checking ${abi_dbgflags_file}])
  if test "${abi_dbgflags_file}" != ""; then
    AC_MSG_NOTICE([loading debug flags for ${abi_result}])
    . "${abi_dbgflags_file}"
  fi
]) # ABI_LOAD_DBGFLAGS



# ABI_LOAD_DIRFLAGS(SUFFIX, COMPILER, VERSION, ARCHITECTURE)
# ----------------------------------------------------------
#
# Looks for per-directory optimization flags for a specified compiler
# and a specified version, running on a specified architecture. Load
# them if found.
#
AC_DEFUN([ABI_LOAD_DIRFLAGS],[
  # Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl
  m4_if([$3], , [AC_FATAL([$0: missing argument 3])])dnl
  m4_if([$4], , [AC_FATAL([$0: missing argument 3])])dnl

  # Init
  abi_result=""
  abi_dirflags_file=""

  # Explore all the possibilities
  for tmp_dirflags_file in \
    "${ac_top_srcdir}/config/compilers/$2_$1/all/all.dir" \
    "${ac_top_srcdir}/config/compilers/$2_$1/all/$4.dir" \
    "${ac_top_srcdir}/config/compilers/$2_$1/$3/all.dir" \
    "${ac_top_srcdir}/config/compilers/$2_$1/$3/$4.dir"; do

    if test -s "${tmp_dirflags_file}"; then
      abi_dirflags_file="${tmp_dirflags_file}"
      abi_result=`echo "${abi_dirflags_file}" | \
        sed -e 's,.*compilers/,,; s,\.dir$,,; s,^\([[^_]]*\)_\([[^/]]*\),\2: \1,'`
    fi
  done

  # Source the file
  #AC_MSG_NOTICE([checking ${abi_dirflags_file}])
  if test "${abi_dirflags_file}" != ""; then
    AC_MSG_NOTICE([loading customizations for ${abi_result}])
    . "${abi_dirflags_file}"
  fi
]) # ABI_LOAD_DIRFLAGS



# ABI_LOAD_OPTFLAGS(SUFFIX, COMPILER, VERSION, ARCHITECTURE)
# ----------------------------------------------------------
#
# Looks for default optimization flags for a specified compiler and
# a specified version, running on a specified architecture. Load them
# if found.
#
AC_DEFUN([ABI_LOAD_OPTFLAGS],[
  # Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl
  m4_if([$3], , [AC_FATAL([$0: missing argument 3])])dnl
  m4_if([$4], , [AC_FATAL([$0: missing argument 3])])dnl

  # Init
  abi_result=""
  abi_optflags_file=""

  # Explore all the possibilities
  for tmp_optflags_file in \
    "${ac_top_srcdir}/config/compilers/generic_$1/all/all.opt" \
    "${ac_top_srcdir}/config/compilers/$2_$1/all/all.opt" \
    "${ac_top_srcdir}/config/compilers/$2_$1/all/$4.opt" \
    "${ac_top_srcdir}/config/compilers/$2_$1/$3/all.opt" \
    "${ac_top_srcdir}/config/compilers/$2_$1/$3/$4.opt"; do

    if test -s "${tmp_optflags_file}"; then
      abi_optflags_file="${tmp_optflags_file}"
      abi_result=`echo "${abi_optflags_file}" | \
        sed -e 's,.*compilers/,,; s,\.opt$,,; s,^\([[^_]]*\)_\([[^/]]*\),\2: \1,'`
    fi
  done

  # Source the file
  #AC_MSG_NOTICE([checking ${abi_optflags_file}])
  if test "${abi_optflags_file}" != ""; then
    AC_MSG_NOTICE([loading optimizations for ${abi_result}])
    . "${abi_optflags_file}"
  else
    AC_MSG_WARN([could not find suitable optimizations])
  fi
]) # ABI_LOAD_OPTFLAGS



# ABI_PROG_MKDIR_P()
# ------------------
#
# Wrapper for the bugged AC_PROG_MKDIR_P macro.
#
AC_DEFUN([ABI_PROG_MKDIR_P],[
  AC_PROG_MKDIR_P
  abi_tmp_mkdir_p=`echo "${MKDIR_P}" | awk '{print [$]1}'`
  if test "${abi_tmp_mkdir_p}" = "config/gnu/install-sh"; then
    AC_MSG_NOTICE([fixing wrong path to mkdir replacement])
    MKDIR_P="${abinit_srcdir}/${MKDIR_P}"
  fi
  unset abi_tmp_mkdir_p
]) # ABI_PROG_MKDIR_P
