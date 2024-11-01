Notes for the split of the source tree
======================================


Relevance
---------

The source tree of ABINIT has gone from a monolithic structure to a set of 3
autonomous blocks: ABINIT Common, LibPAW, and ABINIT Core. The first 2 blocks
represent the shared part of ABINIT, while the latter contains much more
specific code and is tightly bound to the identity of the package.

Knowing exactly how the use of CPP options is distributed within the source
tree of ABINIT is essential to determine the elements on which each build
system will have to focus. For instance, if something is only used in
low-level routines of ABINIT Common, including it in the build system of
ABINIT Core will be useless, and vice-versa. The same holds for
LibPAW-specific features.

The knowledge extracted from the CPP options will actually define the options
of each configure script and the relationships between the build systems,
block by block. It will also help identify a set of small refactoring
operations with the highest impact on the simplicity and velocity of each
build system. Last but not least, it will pinpoint obsolete designs and make
easier the elimination of dead code.

To give a concrete example, since most of the CPP options related to the
Fortran standard version supported by the compiler are used in ABINIT Common,
it is worth taking a few design decisions in order not to repeat the detection
tests within ABINIT Core, thus saving a significant amount of time when
running the configure scripts, one after the other.

Another motivation to simplify the design of build systems adds up in the case
of ABINIT Common and LibPAW: they will be used by more than one DFT code.
Since they will export Fortran modules -- the least standard software
component that has ever existed -- it is of utmost importance to limit the
number of variants produced during the build of each of these blocks.


Accessing information at different levels
-----------------------------------------

One simple strategy to avoid repeating time-consuming detection tests is to
provide a channel that passes information from one block to the others. ABINIT
would benefit from one of the 3 following implementations:

- define a common space where each block would drop the required information
  to be used by the others, e.g. */INSTALL_PREFIX/share/org.abinit/*;
- have each block produce a config file in its source tree / build tree that
  would be included by the others;
- have ABINIT Common install an executable script, e.g. *abinit-config*, that
  the other blocks would call to retrieve essential information and ensure
  consistency at each build step.

Each of these implementations has pros and cons that should be discussed among
the core developers of ABINIT before selecting the most adequate one.


Classifying CPP options
-----------------------

The Python scripts found in the *split/* subdirectory of the *pouillon*
repository of ABINIT provide both the locations of the CPP options in the
source code and some statistics about them. They extract the following
information relevant to the design of build systems:

- how much each CPP option is used in the source code;
- which CPP options are used in each block;
- which CPP options are exclusively used in one of the blocks;
- which CPP options are used in one block plus only once in another block;
- which CPP options are used in one block plus only twice in another block;
- which CPP are used everywhere and should be treated globally.

CPP options only found in one block are fully aligned with the split of the
source tree and should appear in one of the build systems only. Those which
are mostly used in one block only, plus once or twice in another block, are
the best candidates to explore the impact of refactoring operations on the
design of the individual build systems and their mutual relationships. Those
appearing everywhere are those which will define how the build systems will
communicate.


CPP options aligned with the split
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Only used in ABINIT Common:

- ABINIT_SINGLE
- CUDA_COMMON_H
- DEV_LINALG_TIMING
- FC_GNU
- FC_PGI
- GPU_FOUR_HEADER_H
- HAVE_CCLOCK
- HAVE_ELPA_FORTRAN2008
- HAVE_ERRNO_H
- HAVE_FC_ASYNC
- HAVE_FC_BACKTRACE
- HAVE_FC_CONTIGUOUS
- HAVE_FC_CPUTIME
- HAVE_FC_IEEE_ARITHMETIC
- HAVE_FC_IEEE_EXCEPTIONS
- HAVE_FC_IOMSG
- HAVE_FC_ISO_FORTRAN_2008
- HAVE_FC_MACRO_NEWLINE
- HAVE_FC_PRIVATE
- HAVE_FC_PROTECTED
- HAVE_IBM6
- HAVE_INTTYPES_H
- HAVE_LIBTETRA_ABINIT
- HAVE_LINALG_ASL
- HAVE_LINALG_ELPA_2013
- HAVE_LINALG_ELPA_2014
- HAVE_LINALG_ELPA_2015
- HAVE_LINALG_ELPA_2016
- HAVE_LINALG_ELPA_2017
- HAVE_LINALG_ESSL
- HAVE_LINALG_MAGMA_15
- HAVE_LINALG_MKL_IMATCOPY
- HAVE_LINALG_MKL_OMATCOPY
- HAVE_LINALG_PLASMA
- HAVE_MALLOC_H
- HAVE_MATH_H
- HAVE_MPI_TYPE_CREATE_STRUCT
- HAVE_OMP_COLLAPSE
- HAVE_OS_MACOSX
- HAVE_PAPI
- HAVE_STDARG_H
- HAVE_STDDEF_H
- HAVE_STDINT_H
- HAVE_STDIO_H
- HAVE_STDLIB_H
- HAVE_STRING_H
- HAVE_SYS_STAT_H
- HAVE_SYS_TYPES_H
- HAVE_UNISTD_H
- HAVE_VARMACROS
- _ABINIT_CLIB_H
- _ABINIT_COMMON_H
- _CRAY
- _MD5_H
- _XMALLOC_H_
- __GNUC_MINOR__
- __GNUC__
- __STDC_VERSION__
- __STDC__
- __cplusplus

Only used in LibPAW:

- HAVE_AVX_SAFE_MODE
- HAVE_FOX
- HAVE_LIBPAW
- HAVE_LIBPAW_ABINIT
- HAVE_MPI2_INPLACE
- HAVE_YAML
- LIBPAW_HAVE_FOX
- LIBPAW_HAVE_LIBXC
- LIBPAW_HAVE_NETCDF
- LIBPAW_ISO_C_BINDING

Only used in core:

- CUDA_HEADER_H
- CUDA_REC_HEAD_H
- DEBUG_VERBOSE
- DEV_DEBUG_THIS
- DEV_MG_DEBUG_MODE
- DEV_MG_DEBUG_THIS
- DEV_NEW_CODE
- DEV_NEW_HDR
- DEV_RC_BUG
- DEV_USESPLINE
- DEV_YP_DEBUG_PSP
- DEV_YP_VDWXC
- FFT_PRECISION
- HAVE_BIGDFT
- HAVE_BSE_UNPACKED
- HAVE_DFTI
- HAVE_DFTI_MIXED_PRECISION
- HAVE_FC_COMMAND_ARGUMENT
- HAVE_FFTW3
- HAVE_FFTW3_MPI
- HAVE_FFTW3_THREADS
- HAVE_GPU_CUDA_SP
- HAVE_GPU_CUDA_TM
- HAVE_LEVMAR
- HAVE_LIBPSML
- HAVE_LINALG_MKL_THREADS
- HAVE_LOTF
- HAVE_MPI_IALLTOALL
- HAVE_MPI_IBCAST
- HAVE_MPI_IO_DEFAULT
- HAVE_TEST_TIME_PARTITIONING
- HAVE_TRIQS
- HAVE_TRIQS_v1_4
- HAVE_TRIQS_v2_0
- HAVE_WANNIER90
- HAVE_WANNIER90_V1
- HAVE_XML
- M_AB7_KPOINTS_EXPORT_H
- M_AB7_SYMMETRY_EXPORT_H
- READ_FROM_FILE
- _ABINIT_XC_VDW_H

As can be seen, most of the CPP options related to the capabilities of the
Fortran compiler, as well as the interface between Fortran and other languages
like C and C++, are found in the ABINIT Common block. Most of the intricacies
of linear algebra are also found there. These topics will thus be major focus
areas for the build system of ABINIT Common.

The LibPAW block has few CPP options of its own. A finer analysis that
includes information about the history of ABINIT -- that we will not detail
here -- also shows that most of them are now outdated, such as *HAVE_FOX* or
*HAVE_YAML*. Regarding Fortran standards, *LIBPAW_ISO_C_BINDING* could easily
be removed. LibPAW is also the only place in ABINIT where *HAVE_MPI2_INPLACE*
is used. The implementation of the build system of LibPAW will take place in
parallel with a refactoring of its source code.

ABINIT Core is where most of the optional features are managed, as well as
where the different levels of parallelism are interacting. CPP options related
to FFT are of particular importance at this level. All those associated to
experimental developments are found there too. What this tells us *a
posteriori* is that the split of ABINIT into ABINIT Common, LibPAW and ABINIT
Core is actually consistent with the contents of the source code. It is
definitely the core build system of ABINIT that will manage most of the
external dependencies and related aspects like the fallbacks.


Best candidates for a refactoring
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*14_hidewrite/m_cppopts_dumper.F90* could be moved to the upper part of ABINIT
if it were not used by m_errors. Exclusively because of this module, a lot of
detection code that is only useful for ABINIT Core would have to be duplicated
into ABINIT Common.

CPP options which are used almost uniquely in one block and only once in
another block point to a few easy refactoring operations:

- by stopping to use CPP options related to legacy Fortran standards (Fortran
  90/95) and assuming that Fortran compilers support Fortran 2003, the build
  systems of ABINIT Common and LibPAW could be simplified substantially;
- getting rid of the *HAVE_ETSF_IO* CPP option, which has been replaced by a
  direct implementation of the ETSF File Format;
- moving linear algebra-related Fortran modules from ABINIT Common to ABINIT
  Core, after confirming that LibPAW does not require them.

In a second round, the following CPP options could be used exclusively in
ABINIT Core:

- HAVE_GPU_xxx, if the GPU-related modules and C headers, as well as
  m_abi_linalg are moved to ABINIT Core;
- HAVE_LEVMAR, if a low-level C subdirectory is created upwards (level 40);
- HAVE_MEM_PROFILING, if moved away from m_errors (abinit_doctor) and
  abi_common.h;
- HAVE_NETCDF_DEFAULT, if moved away from m_nctk or the latter moved to ABINIT
  Core;
- HAVE_NETCDF_MPI, if m_nctk is moved to ABINIT Core.

The *HAVE_GW_DPC* option could be used exclusively in ABINIT Core if:

- it is defined somewhere else than *10_defs/defs_basis.F90*;
- *28_numeric_noabirule/m_array.F90* is moved upwards.

All these possible refactoring operations only involve small efforts, which is
why they should be discussed among the core developers of ABINIT before
starting ABINIT 9.


Proposed schedule
-----------------

#. Polish the test farm configuration. (YP+JMB, in progress)
#. Make the new build-system interface more user-friendly. (YP+JMB, in
   progress)
#. Move transient/ to  fallbacks/. (YP)
#. Use a fallbacks tarball. (YP)
#. Install fallbacks within ABINIT build dir by default, with FC vendor and
   version. (YP)
#. Make configure use fallbacks automatically. (YP)
#. Write a proper warning for LibPSML (no dynlibs, design flaw). (YP)

Internal ref: YP/2019/Q4/57+77
