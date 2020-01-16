!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_hightemp_top
!! NAME
!!  m_hightemp_top
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

module m_hightemp_top
  use defs_abitypes
  use m_pawcprj
  use m_xmpi
  use m_hamiltonian,    only : gs_hamiltonian_type

  implicit none

  public :: hightemp_prt_cg,hightemp_prt_cprj
contains

  !!****f* ABINIT/m_hightemp/hightemp_prt_cg
  !! NAME
  !! hightemp_prt_cg
  !!
  !! FUNCTION
  !! Printing cg values... INTENT : testing values of cgs, increasing
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
  subroutine hightemp_prt_cg(cprj,eigen,gs_hamk,istep,mband,mcprj,mpi_enreg,natom,nkpt,nsppol,occ)
    ! Arguments -------------------------------
    ! Scalars
    integer, intent(in) :: istep,mband,mcprj,natom,nkpt,nsppol
    type(gs_hamiltonian_type), intent(inout) :: gs_hamk
    type(MPI_type), intent(inout) :: mpi_enreg
    ! Arrays
    real(dp), intent(inout) :: eigen(mband*nkpt*nsppol)
    real(dp), intent(inout) :: occ(mband*nkpt*nsppol)
    type(pawcprj_type),intent(in) :: cprj(natom,mcprj*gs_hamk%usecprj)


  end subroutine hightemp_prt_cg

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
  subroutine hightemp_prt_cprj(cprj,eigen,gs_hamk,istep,mband,mcprj,mpi_enreg,natom,nkpt,nsppol,occ)
    ! Arguments -------------------------------
    ! Scalars
    integer, intent(in) :: istep,mband,mcprj,natom,nkpt,nsppol
    type(gs_hamiltonian_type), intent(inout) :: gs_hamk
    type(MPI_type), intent(inout) :: mpi_enreg
    ! Arrays
    real(dp), intent(inout) :: eigen(mband*nkpt*nsppol)
    real(dp), intent(inout) :: occ(mband*nkpt*nsppol)
    type(pawcprj_type),intent(in) :: cprj(natom,mcprj*gs_hamk%usecprj)

    ! Local variables -------------------------
    ! Scalars
    integer :: iband,iibandpp,mpierr
    integer :: ii,jj,kk,nlmn,n1dim,n2dim
    real(dp) :: mod2cprj
    ! Arrays
    real(dp) :: cprj_red(2,mband)

    ! *********************************************************************
    ! NE FONCTIONNE QUE POUR LE PARALLELISME BANDES

    cprj_red(:,:)=zero
    call xmpi_barrier(xmpi_world)

    do iibandpp=1,mpi_enreg%bandpp
      iband=(mpi_enreg%me_band*mpi_enreg%bandpp)+iibandpp
      cprj_red(1,iband)=cprj(1,iibandpp)%cp(1,1)
      cprj_red(2,iband)=cprj(1,iibandpp)%cp(2,1)
      ! if(iband==1)write(0,*) istep,iband,cprj(1,iibandpp)%cp(2,1),cprj_red(iband)
    end do

    call xmpi_sum(cprj_red,xmpi_world,mpierr)

    if(mpi_enreg%me==0) then
      mod2cprj=zero
      do iband=1, mband
        mod2cprj=cprj_red(1,iband)**2+cprj_red(2,iband)**2
        write(10+istep,*) istep,iband,cprj_red(1,iband),cprj_red(2,iband),&
        & mod2cprj,occ(iband),eigen(iband)
      end do
    end if
  end subroutine hightemp_prt_cprj

end module m_hightemp_top
