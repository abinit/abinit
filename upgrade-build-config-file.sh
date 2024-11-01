#!/bin/bash

# Banner
cat <<EOF
ABINIT build config file upgrader
=================================

Please read the following instructions carefully.

This script will upgrade your former ABINIT build configuration files in the
best possible way, although it will not be able to reproduce complex settings.
Please note that the new config file will be a customized version of the new
config template, located at ~abinit/doc/build/config-template.ac9. The former
custom per-directory Fortran optimization flags (fcflags_opt_*) will not be
retrieved.

When executed without argument, this script will upgrade your default ABINIT
config file. Any config file can also be specified on the command line, one
file at a time. The resulting upgraded file will be written in the current
directory with the same name, if it does not exist yet. This script will
strictly avoid to overwrite any existing file.

To work properly, this script must be run from the directory where it is
located. Running it from any other location will fail, unless you modify the
script to do so.

A typical limitation of the automatic conversion of the former config file is
that the following will happen. For instance, the configuration of LibXC will
be translated from:

    with_dft_flavor="libxc"
    with_libxc_incs="-I/path/to/libxc/include"
    with_libxc_libs="-L/path/to/libxc/lib -lxcf90 -lxc"

to:

    with_libxc="yes"
    LIBXC_CPPFLAGS="-I/path/to/libxc/include"
    LIBXC_FCFLAGS="-I/path/to/libxc/include"
    LIBXC_LIBS="-L/path/to/libxc/lib -lxcf90 -lxc"

while a more appropriate and compact way to configure it would be:

    with_libxc="/path/to/libxc"

The latter is the recommended way to configure external packages. This makes
it easy to use modern dependency management frameworks with ABINIT, e.g. if
you have installed LibXC with EasyBuild, you just have to set the following:

    with_libxc="\${EBROOTLIBXC}"

In the futture, this step will even become superfluous, since the build system
will be able to automatically detect frameworks like EasyBuild, Spack,
PkgConfig, and the ESL Bundle.

You are more than warmly invited to carefully review the new config file, as
some specificities you may have relied on until now might have slipped through
the cracks of the upgrade procedure. Before doing so, please read carefully
all the warnings that may now come from this script.

EOF

# Select config template
cfg_template="doc/build/config-template.ac9"
if test ! -s "${cfg_template}"; then
  echo "Error: File not found: ${cfg_template}" >&2
  exit 1
fi

# Select old config file
if test "${#}" -ge 1; then
  cfg_old="${1}"
else
  declare -a PathArray=("." "./src" "${HOME}/.abinit/build" "/etc/abinit/build")
  Hostname="`hostname | sed -e 's/\..*//'`"
  for val in ${PathArray[@]}; do
    cfg_old="$val/${Hostname}.ac"
    if test -f $cfg_old; then
       echo "found : $cfg_old"
       break
    fi
  done
fi

# Check that old config file exists
if test ! -s "${cfg_old}"; then
  echo "Error: File not found: ${cfg_old}" >&2
  exit 1
fi

# Select new config file
cfg_new=`basename "${cfg_old}"`
cfg_new="${cfg_new%.ac}.ac9"

# Check that new config file doesn't exist
if test -e "${cfg_new}"; then
  echo "Error: New config file already exists: ${cfg_new}" >&2
  exit 2
fi

# Import old configuration
. "${cfg_old}"

                   # ------------------------------------- #

# Checkpoint
echo "=== Starting conversion of ${cfg_old}"
echo ""

# Warn about removed options
test -z "${enable_64bit_flags}" || \
  echo "Warning: enable_64bit_flags removed, use *FLAGS_EXTRA instead" >&2
test -z "${enable_connectors}" || \
  echo "Warning: enable_connectors removed" >&2
test -z "${enable_fast_check}" || \
  echo "Warning: enable_fast_check removed" >&2
test -z "${with_atompaw_bins}" || \
  echo "Warning: with_atompaw_bins removed" >&2
test -z "${with_atompaw_incs}" || \
  echo "Warning: with_atompaw_incs removed" >&2
test -z "${with_atompaw_libs}" || \
  echo "Warning: with_atompaw_libs removed" >&2
test -z "${with_math_flavor}" || \
  echo "Warning: with_math_flavor removed" >&2
test -z "${with_math_incs}" || \
  echo "Warning: with_math_incs removed" >&2
test -z "${with_math_libs}" || \
  echo "Warning: with_math_libs removed" >&2
test -z "${with_timer_flavor}" || \
  echo "Warning: with_timer_flavor removed" >&2
test -z "${with_wannier90_bins}" || \
  echo "Warning: with_wannier90_bins removed" >&2
test -z "${with_yaml_incs}" || \
  echo "Warning: with_yaml_incs removed" >&2
test -z "${with_yaml_libs}" || \
  echo "Warning: with_yaml_libs removed" >&2

                   # ------------------------------------- #

# Notify about changed options
if test "${with_mpi_level}" = "1"; then
  echo "Warning: MPI level 1 dropped, level 2 will be used instead" >&2
  unset with_mpi_level
fi

                   # ------------------------------------- #

# Translate enable_debug
with_debug_flavor="${enable_debug}"
test "${enable_debug}" = "no" && with_debug_flavor="none"
test "${enable_debug}" = "yes" && with_debug_flavor="custom"

# Translate enable_fallbacks
with_fallbacks="${enable_fallbacks}"

# Translate enable_gpu
with_gpu="${enable_gpu}"

# Translate enable_mpi
with_mpi="${enable_mpi}"

# Translate enable_optim
with_debug_flavor="${enable_optim}"
test "${enable_optim}" = "no" && with_debug_flavor="none"
test "${enable_optim}" = "yes" && with_debug_flavor="custom"

# Translate with_algo_flavor
test "${with_algo_flavor}" = "levmar" && with_levmar="yes"

# Translate with_dft_flavor
dft_pkgs=`echo "${with_dft_flavor}" | sed -e 's/\+/ /g'`
for pkg in ${dft_pkgs}; do
  test "${pkg}" != "atompaw" && eval with_${pkg}="yes"
done
test "${with_psml}" = "yes" && with_libpsml="yes"
tmp_atompaw=`echo "${dft_pkgs}" | grep 'atompaw'`
if test "${tmp_atompaw}" != ""; then
  echo "Warning: build-system support for AtomPAW removed" >&2
fi

# Translate with_trio_flavor
trio_pkgs=`echo "${with_trio_flavor}" | sed -e 's/\+/ /g'`
for pkg in ${trio_pkgs}; do
  test "${pkg}" != "yaml" && eval with_${pkg}="yes"
done

# Translate with_fft_incs
FFT_CPPFLAGS="${with_fft_incs}"
FFT_FCFLAGS="${with_fft_incs}"

# Translate with_fft_libs
FFT_LIBS="${with_fft_libs}"

# Translate with_gpu_cflags
GPU_CFLAGS="${with_gpu_cflags}"

# Translate with_gpu_cppflags
GPU_CPPFLAGS="${with_gpu_cppflags}"

# Translate with_gpu_incs
if test -z "${GPU_CPPFLAGS}"; then
 GPU_CPPFLAGS="${with_gpu_incs}"
else
 GPU_CPPFLAGS="${GPU_CPPFLAGS} ${with_gpu_incs}"
fi

# Translate with_gpu_ldflags
GPU_LDFLAGS="${with_gpu_ldflags}"

# Translate with_gpu_libs
GPU_LIBS="${with_gpu_libs}"

# Translate with_gpu_prefix
test ! -z "${with_gpu_prefix}" && with_gpu="${with_gpu_prefix}"

# Translate with_levmar_incs
LEVMAR_CPPFLAGS="${with_levmar_incs}"
LEVMAR_FCFLAGS="${with_levmar_incs}"

# Translate with_levmar_libs
LEVMAR_LIBS="${with_levmar_libs}"

# Translate with_linalg_incs
LINALG_FCFLAGS="${with_linalg_incs}"

# Translate with_linalg_libs
LINALG_LIBS="${with_linalg_libs}"

# Translate with_libxc_incs
LIBXC_CPPFLAGS="${with_libxc_incs}"
LIBXC_FCFLAGS="${with_libxc_incs}"

# Translate with_libxc_libs
LIBXC_LIBS="${with_libxc_libs}"

# Translate with_mpi_incs
MPI_CPPFLAGS="${with_mpi_incs}"

# Translate with_mpi_libs
MPI_LIBS="${with_mpi_libs}"

# Translate with_mpi_prefix
test ! -z "${with_mpi_prefix}" && with_mpi="${with_mpi_prefix}"

# Translate with_libpsml_incs
LIBPSML_FCFLAGS="${with_psml_incs}"

# Translate with_libpsml_libs
LIBPSML_LIBS="${with_psml_libs}"

# Translate with_papi_incs
PAPI_CPPFLAGS="${with_papi_incs}"
PAPI_FCFLAGS="${with_papi_incs}"

# Translate with_papi_libs
PAPI_LIBS="${with_papi_libs}"

# Translate with_triqs_incs
TRIQS_CPPFLAGS="${with_triqs_incs}"
TRIQS_FCFLAGS="${with_triqs_incs}"

# Translate with_triqs_libs
TRIQS_LIBS="${with_triqs_libs}"

# Translate with_wannier90_incs
WANNIER90_FCFLAGS="${with_wannier90_incs}"

# Translate with_wannier90_libs
WANNIER90_LIBS="${with_wannier90_libs}"

                   # ------------------------------------- #

# Look for a Python interpreter
if test -z "${PYTHON}"; then
  test -x "/usr/bin/python3" && PYTHON="/usr/bin/python3"
  test -z "${PYTHON}" && PYTHON="python"
fi

# Create new config file
${PYTHON} >"${cfg_new}" <<EOF
import re

abinit_options = {
    "AR": "${AR}",
    "ARFLAGS": "${ARFLAGS}",
    "ARFLAGS_DEBUG": "${ARFLAGS_DEBUG}",
    "ARFLAGS_EXTRA": "${ARFLAGS_EXTRA}",
    "ARFLAGS_OPTIM": "${ARFLAGS_OPTIM}",
    "BIGDFT_FCFLAGS": "${BIGDFT_FCFLAGS}",
    "BIGDFT_LIBS": "${BIGDFT_LIBS}",
    "CC": "${CC}",
    "CC_LDFLAGS": "${CC_LDFLAGS}",
    "CC_LDFLAGS_DEBUG": "${CC_LDFLAGS_DEBUG}",
    "CC_LDFLAGS_EXTRA": "${CC_LDFLAGS_EXTRA}",
    "CC_LDFLAGS_OPTIM": "${CC_LDFLAGS_OPTIM}",
    "CC_LIBS": "${CC_LIBS}",
    "CC_LIBS_DEBUG": "${CC_LIBS_DEBUG}",
    "CC_LIBS_EXTRA": "${CC_LIBS_EXTRA}",
    "CC_LIBS_OPTIM": "${CC_LIBS_OPTIM}",
    "CFLAGS": "${CFLAGS}",
    "CFLAGS_DEBUG": "${CFLAGS_DEBUG}",
    "CFLAGS_EXTRA": "${CFLAGS_EXTRA}",
    "CFLAGS_OPTIM": "${CFLAGS_OPTIM}",
    "CPP": "${CPP}",
    "CPPFLAGS": "${CPPFLAGS}",
    "CPPFLAGS_DEBUG": "${CPPFLAGS_DEBUG}",
    "CPPFLAGS_EXTRA": "${CPPFLAGS_EXTRA}",
    "CPPFLAGS_OPTIM": "${CPPFLAGS_OPTIM}",
    "CXX": "${CXX}",
    "CXXFLAGS": "${CXXFLAGS}",
    "CXXFLAGS_DEBUG": "${CXXFLAGS_DEBUG}",
    "CXXFLAGS_EXTRA": "${CXXFLAGS_EXTRA}",
    "CXXFLAGS_OPTIM": "${CXXFLAGS_OPTIM}",
    "CXX_LDFLAGS": "${CXX_LDFLAGS}",
    "CXX_LDFLAGS_DEBUG": "${CXX_LDFLAGS_DEBUG}",
    "CXX_LDFLAGS_EXTRA": "${CXX_LDFLAGS_EXTRA}",
    "CXX_LDFLAGS_OPTIM": "${CXX_LDFLAGS_OPTIM}",
    "CXX_LIBS": "${CXX_LIBS}",
    "CXX_LIBS_DEBUG": "${CXX_LIBS_DEBUG}",
    "CXX_LIBS_EXTRA": "${CXX_LIBS_EXTRA}",
    "CXX_LIBS_OPTIM": "${CXX_LIBS_OPTIM}",
    "F77": "${F77}",
    "FC": "${FC}",
    "FCFLAGS": "${FCFLAGS}",
    "FCFLAGS_DEBUG": "${FCFLAGS_DEBUG}",
    "FCFLAGS_EXTRA": "${FCFLAGS_EXTRA}",
    "FCFLAGS_FIXEDFORM": "${FCFLAGS_FIXEDFORM}",
    "FCFLAGS_FREEFORM": "${FCFLAGS_FREEFORM}",
    "FCFLAGS_HINTS": "${FCFLAGS_HINTS}",
    "FCFLAGS_MODDIR": "${FCFLAGS_MODDIR}",
    "FCFLAGS_OPENMP": "${FCFLAGS_OPENMP}",
    "FCFLAGS_OPTIM": "${FCFLAGS_OPTIM}",
    "FC_LDFLAGS": "${FC_LDFLAGS}",
    "FC_LDFLAGS_DEBUG": "${FC_LDFLAGS_DEBUG}",
    "FC_LDFLAGS_EXTRA": "${FC_LDFLAGS_EXTRA}",
    "FC_LDFLAGS_OPTIM": "${FC_LDFLAGS_OPTIM}",
    "FC_LIBS": "${FC_LIBS}",
    "FC_LIBS_DEBUG": "${FC_LIBS_DEBUG}",
    "FC_LIBS_EXTRA": "${FC_LIBS_EXTRA}",
    "FC_LIBS_OPTIM": "${FC_LIBS_OPTIM}",
    "FFT_CPPFLAGS": "${FFT_CPPFLAGS}",
    "FFT_FCFLAGS": "${FFT_FCFLAGS}",
    "FFT_LIBS": "${FFT_LIBS}",
    "FPP": "${FPP}",
    "FPPFLAGS": "${FPPFLAGS}",
    "FPPFLAGS_DEBUG": "${FPPFLAGS_DEBUG}",
    "FPPFLAGS_EXTRA": "${FPPFLAGS_EXTRA}",
    "FPPFLAGS_OPTIM": "${FPPFLAGS_OPTIM}",
    "GPU_CFLAGS": "${GPU_CFLAGS}",
    "GPU_CPPFLAGS": "${GPU_CPPFLAGS}",
    "GPU_FCFLAGS": "${GPU_FCFLAGS}",
    "GPU_LDFLAGS": "${GPU_LDFLAGS}",
    "GPU_LIBS": "${GPU_LIBS}",
    "LD": "${LD}",
    "LEVMAR_CPPFLAGS": "${LEVMAR_CPPFLAGS}",
    "LEVMAR_FCFLAGS": "${LEVMAR_FCFLAGS}",
    "LEVMAR_LIBS": "${LEVMAR_LIBS}",
    "LIBPSML_FCFLAGS": "${LIBPSML_FCFLAGS}",
    "LIBPSML_LIBS": "${LIBPSML_LIBS}",
    "LIBXC_CPPFLAGS": "${LIBXC_CPPFLAGS}",
    "LIBXC_FCFLAGS": "${LIBXC_FCFLAGS}",
    "LIBXC_LIBS": "${LIBXC_LIBS}",
    "LINALG_FCFLAGS": "${LINALG_FCFLAGS}",
    "LINALG_LIBS": "${LINALG_LIBS}",
    "MPI_CPPFLAGS": "${MPI_CPPFLAGS}",
    "MPI_FCFLAGS": "${MPI_FCFLAGS}",
    "MPI_LIBS": "${MPI_LIBS}",
    "NETCDF_FCFLAGS": "${NETCDF_FCFLAGS}",
    "NETCDF_LIBS": "${NETCDF_LIBS}",
    "NM": "${NM}",
    "NVCC": "${NVCC}",
    "NVCC_CFLAGS": "${NVCC_CFLAGS}",
    "NVCC_CPPFLAGS": "${NVCC_CPPFLAGS}",
    "NVCC_LDFLAGS": "${NVCC_LDFLAGS}",
    "NVCC_LIBS": "${NVCC_LIBS}",
    "PAPI_CPPFLAGS": "${PAPI_CPPFLAGS}",
    "PAPI_FCFLAGS": "${PAPI_FCFLAGS}",
    "PAPI_LIBS": "${PAPI_LIBS}",
    "PYFLAGS": "${PYFLAGS}",
    "PYTHON_CPPFLAGS": "${PYTHON_CPPFLAGS}",
    "RANLIB": "${RANLIB}",
    "TRIQS_CPPFLAGS": "${TRIQS_CPPFLAGS}",
    "TRIQS_FCFLAGS": "${TRIQS_FCFLAGS}",
    "TRIQS_LIBS": "${TRIQS_LIBS}",
    "WANNIER90_FCFLAGS": "${WANNIER90_FCFLAGS}",
    "WANNIER90_LIBS": "${WANNIER90_LIBS}",
    "XPP": "${XPP}",
    "XPPFLAGS": "${XPPFLAGS}",
    "XPPFLAGS_DEBUG": "${XPPFLAGS_DEBUG}",
    "XPPFLAGS_EXTRA": "${XPPFLAGS_EXTRA}",
    "XPPFLAGS_OPTIM": "${XPPFLAGS_OPTIM}",
    "enable_avx_safe_mode": "${enable_avx_safe_mode}",
    "enable_bse_unpacked": "${enable_bse_unpacked}",
    "enable_cclock": "${enable_cclock}",
    "enable_expert": "${enable_expert}",
    "enable_exports": "${enable_exports}",
    "enable_fc_wrapper": "${enable_fc_wrapper}",
    "enable_gw_dpc": "${enable_gw_dpc}",
    "enable_hints": "${enable_hints}",
    "enable_libtetra": "${enable_libtetra}",
    "enable_lotf": "${enable_lotf}",
    "enable_maintainer_checks": "${enable_maintainer_checks}",
    "enable_memory_profiling": "${enable_memory_profiling}",
    "enable_mpi_inplace": "${enable_mpi_inplace}",
    "enable_mpi_io": "${enable_mpi_io}",
    "enable_mpi_io_default": "${enable_mpi_io_default}",
    "enable_netcdf_default": "${enable_netcdf_default}",
    "enable_openmp": "${enable_openmp}",
    "enable_stdin": "${enable_stdin}",
    "enable_timer": "${enable_timer}",
    "enable_triqs_v1_4": "${enable_triqs_v1_4}",
    "enable_triqs_v2_0": "${enable_triqs_v2_0}",
    "enable_wannier90_v1": "${enable_wannier90_v1}",
    "enable_zdot_bugfix": "${enable_zdot_bugfix}",
    "fcflags_opt_95_drive": "${fcflags_opt_95_drive}",
    "prefix": "${prefix}",
    "with_bigdft": "${with_bigdft}",
    "with_debug_flavor": "${with_debug_flavor}",
    "with_fallbacks": "${with_fallbacks}",
    "with_fc_vendor": "${with_fc_vendor}",
    "with_fc_version": "${with_fc_version}",
    "with_fft": "${with_fft}",
    "with_fft_flavor": "${with_fft_flavor}",
    "with_gpu": "${with_gpu}",
    "with_gpu_flavor": "${with_gpu_flavor}",
    "with_levmar": "${with_levmar}",
    "with_libpsml": "${with_libpsml}",
    "with_libxc": "${with_libxc}",
    "with_linalg": "${with_linalg}",
    "with_linalg_flavor": "${with_linalg_flavor}",
    "with_mpi": "${with_mpi}",
    "with_mpi_level": "${with_mpi_level}",
    "with_netcdf": "${with_netcdf}",
    "with_optim_flavor": "${with_optim_flavor}",
    "with_papi": "${with_papi}",
    "with_triqs": "${with_triqs}",
    "with_wannier90": "${with_wannier90}",
}

cfg_text = ""
with open("${cfg_template}", "r") as cfg_file:
    for line in cfg_file.readlines():
        for key, val in abinit_options.items():
            if ( re.match("#{}=".format(key), line) and (val != "") ):
                line = "{}=\"{}\"\n".format(key, val)
                break
        cfg_text += line

print(cfg_text)
EOF

# Checkpoint
cat <<EOF

=== Configuration saved to ${cfg_new}

You can now review the output of this script and the new config file. You may
also want to look at the output of the configure script when running it for
the first time with this new file, in particular the early sections reporting
about configure options and environment.

Enjoy your shiny new ABINIT configuration!

EOF
