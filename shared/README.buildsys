Notes for the split of the source tree
======================================


CPP options aligned with the split
----------------------------------

Only used in shared:

- HAVE_AVX_SAFE_MODE (common, libpaw).
- HAVE_CCLOCK (common).
- HAVE_LIBTETRA_ABINIT (common).
- HAVE_LIBXC (common, libpaw).
- HAVE_MPI2_INPLACE (common, libpaw).
- HAVE_PAPI (common).

Only used in src:

- HAVE_BIGDFT
- HAVE_BSE_UNPACKED
- HAVE_FFTW3
- HAVE_LIBPSML
- HAVE_LOTF
- HAVE_MPI_IO_DEFAULT
- HAVE_TRIQS_v1_4
- HAVE_TRIQS_v2_0
- HAVE_WANNIER90
- HAVE_XML
- USE_MACROAVE


CPP options
-----------

Could be only in src:

- HAVE_GPU_*, if m_abi_linalg is moved upwards.
- HAVE_GW_DPC, if removed from defs_basis and abi_common.h.
- HAVE_LEVMAR, if a low-level C subdirectory is created upwards.
- HAVE_MEM_PROFILING, if removed from m_errors (abinit_doctor) and abi_common.h.
- HAVE_NETCDF_DEFAULT, if removed from m_nctk or the latter moved upwards.
- HAVE_NETCDF_MPI, if m_nctk is moved upwards.


Case of m_cppopts_dumper
------------------------

Could be moved to the upper part of ABINIT if it were not used by m_errors.

