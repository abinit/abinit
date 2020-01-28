!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_xomp
!! NAME
!! m_xomp
!!
!! FUNCTION
!!  Thin wrappers and tools for OpenMP parallelization.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_xomp

 use defs_basis,    only : std_out
#ifdef HAVE_OPENMP
 use omp_lib
#endif

 implicit none

 private

 public :: xomp_show_info
 public :: xomp_get_max_threads
 public :: xomp_get_thread_num
 public :: xomp_get_num_threads
 public :: xomp_set_num_threads
 public :: xomp_in_parallel
 public :: xomp_get_num_cores_node

!----------------------------------------------------------------------

CONTAINS  !=========================================================================================================================

!!***

!----------------------------------------------------------------------

!!****f* m_xomp/xomp_show_info
!! NAME
!!  xomp_show_info
!!
!! FUNCTION
!!  Printout of the most important OMP environment variables.
!!
!! INPUTS
!!  unit=unit number for writing. The named constant dev_null defined in defs_basis can be used to avoid any printing.
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      abinit,fftprof,lapackprof
!!
!! CHILDREN
!!      omp_set_num_threads
!!
!! SOURCE

subroutine xomp_show_info(unit)

!Arguments-------------------------
 integer,optional,intent(in) :: unit

!Local variables-------------------
 integer :: my_unt

! *************************************************************************

 my_unt = std_out; if (PRESENT(unit)) my_unt=unit

#ifdef HAVE_OPENMP
 write(my_unt,'(/,a)')  "  ==== OpenMP parallelism is ON ===="
 write(my_unt,'(a,i0)') "- Max_threads:       ",xomp_get_max_threads()
 write(my_unt,'(a,i0)') "- Num_threads:       ",xomp_get_num_threads(open_parallel=.True.)
 write(my_unt,'(a,i0)') "- Num_procs:         ",omp_get_num_procs()
 write(my_unt,'(a,l1)') "- Dynamic:           ",omp_get_dynamic()
 write(my_unt,'(a,l1)') "- Nested:            ",omp_get_nested()
 !write(my_unt,'(a,i0)')"- Thread_limit:      ",omp_get_thread_limit()
 !write(my_unt,'(a,i0)')"- Max_active_levels: ",omp_get_max_active_levels()
#else
 write(my_unt,'(/,a)')  "  ==== OpenMP parallelism is OFF ===="
#endif

 write(my_unt,*)""

end subroutine xomp_show_info
!!***

!----------------------------------------------------------------------

!!****f* m_xomp/xomp_get_max_threads
!! NAME
!!  xomp_get_max_threads
!!
!! FUNCTION
!!  Wrapper for omp_get_max_threads.
!!
!! OUTPUT
!!  Return the maximum number of threads used for the current parallel region that
!!  does not use the clause num_threads. Return 1 if OMP is disabled.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function xomp_get_max_threads()

!Arguments ------------------------------------
!scalars
 integer :: xomp_get_max_threads

! *************************************************************************

#ifdef HAVE_OPENMP
 xomp_get_max_threads = omp_get_max_threads()
#else
 xomp_get_max_threads = 1
#endif

end function xomp_get_max_threads
!!***

!----------------------------------------------------------------------

!!****f* m_xomp/xomp_get_thread_num
!! NAME
!!  xomp_get_thread_num
!!
!! FUNCTION
!!  Wrapper for omp_get_thread_num
!!  Returns a unique thread identification number within the current team.
!!  In a sequential parts of the program, omp_get_thread_num always returns 0.
!!  In parallel regions the return value varies from 0 to omp_get_num_threads-1 inclusive.
!!  The return value of the master thread of a team is always 0.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function xomp_get_thread_num()

!Arguments ------------------------------------
!scalars
 integer :: xomp_get_thread_num

! *************************************************************************

#ifdef HAVE_OPENMP
 xomp_get_thread_num = omp_get_thread_num()
#else
 xomp_get_thread_num = 0
#endif

end function xomp_get_thread_num
!!***

!----------------------------------------------------------------------

!!****f* m_xomp/xomp_get_num_threads
!! NAME
!!  xomp_get_num_threads
!!
!! FUNCTION
!!  Wrapper for omp_get_num_threads.
!!  The omp_get_num_threads function returns the number of threads in the team currently executing
!!  the parallel region from which it is called. The function binds to the closest enclosing PARALLEL directive.
!!  The omp_set_num_threads subroutine and the OMP_NUM_THREADS environment variable control the number of threads in a team.
!!  If you do not explicitly set the number of threads, the run-time environment will use the number of online processors
!!  on the machine by default. If you call omp_get_num_threads from a serial portion of your program or from a
!!  nested parallel region that is serialized, the function returns 1.
!!
!! INPUTS
!!  [open_parallel]= If .TRUE., a temporary OMP parallel region will be open and omp_get_num_threads
!!                   will be called inside this region.
!!                   Default to .FALSE. so that we have consistent with the OMP API.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function xomp_get_num_threads(open_parallel) result(nthreads)

!Arguments ------------------------------------
!scalars
 logical,optional,intent(in) :: open_parallel
 integer :: nthreads

!Local variables-------------------------------
!scalars
 logical :: do_open

! *************************************************************************

 do_open = .FALSE.; if (PRESENT(open_parallel)) do_open = open_parallel
#ifdef HAVE_OPENMP

 if (do_open .and. .not.xomp_in_parallel()) then
!$OMP PARALLEL
!$OMP SINGLE
  nthreads = omp_get_num_threads()
!$OMP END SINGLE
!$OMP END PARALLEL
 else
   nthreads = omp_get_num_threads()
 end if

#else
 nthreads = 1
#endif

end function xomp_get_num_threads
!!***

!----------------------------------------------------------------------

!!****f* m_xomp/xomp_set_num_threads
!! NAME
!!  xomp_set_num_threads
!!
!! FUNCTION
!!  Specifies the number of threads used by default in subsequent parallel sections,
!!  if those do not specify a num_threads clause. The argument of xomp_set_num_threads shall be a positive integer.
!!
!! INPUTS
!!  nthreads = number of threads
!!
!! SIDE EFFECTS
!!  See description.
!!
!! PARENTS
!!      lapackprof,m_argparse,m_fft,m_fft_prof
!!
!! CHILDREN
!!      omp_set_num_threads
!!
!! SOURCE

subroutine xomp_set_num_threads(nthreads)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nthreads

! *************************************************************************

#ifdef HAVE_OPENMP
 call omp_set_num_threads(nthreads)
#else
 if (.FALSE.) write(std_out,*) nthreads
#endif

end subroutine xomp_set_num_threads
!!***

!----------------------------------------------------------------------

!!****f* m_xomp/xomp_in_parallel
!! NAME
!!  xomp_in_parallel
!!
!! FUNCTION
!!  This function returns true if are currently running in parallel, false otherwise
!!
!! CHILDREN
!!      omp_in_parallel
!!
!! SOURCE

function xomp_in_parallel() result(ans)

!Arguments-------------------------
 logical :: ans

! *************************************************************************

#ifdef HAVE_OPENMP
 ans = omp_in_parallel()
#else
 ans = .FALSE.
#endif

end function xomp_in_parallel
!!***

!----------------------------------------------------------------------

!!****f* m_xomp/xomp_get_num_cores_node
!! NAME
!!  xomp_get_num_cores_node
!!
!! FUNCTION
!!  Wrapper for omp_get_num_procs
!!
!! OUTPUT
!!  Return the maximum number of cores in one shared memory system
!!  Return 0 if OMP is disabled.
!!
!! PARENTS
!!  m_abi_linalg
!!
!! CHILDREN
!!
!! SOURCE

function xomp_get_num_cores_node()

!Arguments ------------------------------------
!scalars
 integer :: xomp_get_num_cores_node

! *************************************************************************

#ifdef HAVE_OPENMP
 xomp_get_num_cores_node=omp_get_thread_limit()
 !We test if thread_limit has been set (if not it should be a large value)
 ! In 2012, 4096 cores is the biggest known shared memory system
 if(xomp_get_num_cores_node > 4096) then
    !so if not set, we used system 'num procs' values which should be the default case
    xomp_get_num_cores_node=omp_get_num_procs()
 end if
#else
 xomp_get_num_cores_node=0
#endif

end function xomp_get_num_cores_node

!----------------------------------------------------------------------

END MODULE m_xomp
!!***
