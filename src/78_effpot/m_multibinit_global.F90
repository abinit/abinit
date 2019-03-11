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
    implicit none
!!****m*
    integer, parameter :: master =0
    logical :: iam_master =.False.
    integer :: my_rank, comm, nproc, ierr
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

  end module m_multibinit_global
