!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_mpi_global
!! NAME
!! m_mpi_global
!!
!! FUNCTION
!! This module contains global variables for multibinit
!!
!! Datatypes:
!!
!! * mpi_global_t
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .

#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"

  module m_multibinit_global
    use defs_basis
    use m_random_xoroshiro128plus, only: rng_t
    use m_xmpi
    use m_errors
    implicit none
!!****m*
    integer, parameter :: master =0
    logical :: iam_master =.False.
    integer :: my_rank, comm, nproc, ierr

    
    type timer_t
       character(20) :: label="mytimer"
       real(dp) :: t0=0.0_dp, t1=0.0_dp, t=0.0_dp
     contains
       procedure :: initialize => timer_initialize
       procedure :: tic => timer_tic
       procedure :: toc => timer_toc
       procedure :: pause => timer_pause
    end type timer_t

    type(timer_t), save :: timer

  contains
    subroutine init_multibinit_global()
      comm = xmpi_world
      nproc = xmpi_comm_size(comm)
      my_rank = xmpi_comm_rank(comm)
      iam_master = (my_rank == master)
    end subroutine init_multibinit_global

    subroutine debug_mpi(msg)
      character(len=*), intent(in) ::msg
      integer :: i
      i=0
      print *, "ierr", ierr
      print *, "rank: ", my_rank,  msg
      call xmpi_bcast(i, master, comm, ierr)
    end subroutine debug_mpi
    
    subroutine timer_initialize(self, label)
      class(timer_t), intent(inout) :: self
      character(*), intent(in) :: label
      self%label=label
    end subroutine timer_initialize

    subroutine timer_tic(self )
      class(timer_t), intent(inout) :: self
      self%t=0.0_dp
      call cpu_time(self%t0)
    end subroutine timer_tic

    subroutine timer_pause(self)
      class(timer_t), intent(inout) :: self
      call cpu_time(self%t1)
      self%t=self%t+self%t1-self%t0
    end subroutine timer_pause

    subroutine timer_toc(self)
      class(timer_t), intent(inout) :: self
      call cpu_time(self%t1)
      print *, self%t1
      self%t=self%t+self%t1-self%t0
      print *, "Total time of ", self%label, ": ", self%t
    end subroutine timer_toc

  end module m_multibinit_global
