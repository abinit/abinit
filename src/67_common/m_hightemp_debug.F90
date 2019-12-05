!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_hightemp_debug
!! NAME
!!  m_hightemp_debug
!!
!! FUNCTION
!!  This module provides routines to debug m_hightemp module
!!
!! COPYRIGHT
!!  Copyright (C) 2018-2019 ABINIT group (A. Blanchet)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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

module m_hightemp_debug
  use defs_abitypes
  use m_pawcprj
  use m_xmpi
  use m_hamiltonian,    only : gs_hamiltonian_type

  implicit none

  public :: hightemp_prt_cprj
contains

  !!****f* ABINIT/m_hightemp/hightemp_prt_cprj
  !! NAME
  !! hightemp_prt_cprj
  !!
  !! FUNCTION
  !! Printing cprj values... INTENT : testing values of cprj, increasing
  !! the number of bands, and the temperature.
  !!
  !! INPUTS
  !!
  !! OUTPUT
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine hightemp_prt_cprj(cprj,gs_hamk,istep,mband,mcprj,mpi_enreg,natom)
    ! Arguments -------------------------------
    ! Scalars
    integer, intent(in) :: istep,mband,mcprj,natom
    type(gs_hamiltonian_type), intent(inout) :: gs_hamk
    type(MPI_type), intent(inout) :: mpi_enreg
    ! Arrays
    type(pawcprj_type),intent(in) :: cprj(natom,mcprj*gs_hamk%usecprj)

    ! Local variables -------------------------
    ! Scalars
    integer :: iband,iibandpp,mpierr
    integer :: ii,jj,kk,nlmn,n1dim,n2dim
    ! Arrays
    real(dp) :: cprj_red(2,mband)

    ! *********************************************************************
    ! NE FONCTIONNE QUE POUR LE PARALLELISME BANDES

    call pawcprj_output(cprj)

    cprj_red(:,:)=zero
    call xmpi_barrier(xmpi_world)

    do iibandpp=1,mpi_enreg%bandpp
      iband=(mpi_enreg%me_band*mpi_enreg%bandpp)+iibandpp
      cprj_red(1,iband)=cprj(1,iibandpp)%cp(1,1)
      cprj_red(2,iband)=cprj(1,iibandpp)%cp(2,1)
      ! if(iband==1)write(0,*) istep,iband,cprj(1,iibandpp)%cp(2,1),cprj_red(iband)
    end do

    call xmpi_sum(cprj_red,xmpi_world,mpierr)
    call xmpi_barrier(xmpi_world)

    if(mpi_enreg%me==0) then
      do iband=1, mband
        write(0,*) istep, iband, cprj_red(1,iband), cprj_red(2,iband), sqrt(cprj_red(1,iband)**2 + cprj_red(2,iband)**2)
      end do
    end if

  end subroutine hightemp_prt_cprj


end module m_hightemp_debug
