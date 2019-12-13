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
  use m_hightemp
  use m_pawcprj
  use m_xmpi
  use m_hamiltonian,    only : gs_hamiltonian_type
  use m_kg,             only : mkkin
  use m_spacepar,       only : meanvalue_g

  implicit none

  public :: hightemp_prt_cprj,hightemp_get_e_shiftfactor
contains

  !!****f* ABINIT/m_hightemp/hightemp_get_e_shiftfactor
  !! NAME
  !! hightemp_get_e_shiftfactor
  !!
  !! FUNCTION
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
  subroutine hightemp_get_e_shiftfactor(cg,ecut,ecutsm,effmass_free,eigen,gmet,hightemp,&
  & istwfk,kg,kptns,mband,mcg,mkmem,mpi_enreg,mpw,my_nspinor,nband,nkpt,nsppol,npwarr)
    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: mband,mcg,mpw,mkmem,my_nspinor,nkpt,nsppol
    real(dp),intent(in) :: ecut,ecutsm,effmass_free
    type(MPI_type),intent(inout) :: mpi_enreg
    type(hightemp_type),pointer,intent(inout) :: hightemp
    ! Arrays
    integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt)
    integer,intent(in) :: kg(3,mpw*mkmem)
    real(dp),intent(in) :: cg(2,mcg),eigen(mband*nkpt*nsppol)
    real(dp),intent(in) :: kptns(3,nkpt),gmet(3,3)

    ! Local variables -------------------------
    ! Scalars
    integer :: bdtot_index,blocksize,iband,iblock,iblocksize
    integer :: ikpt,isppol,mpierr,nband_k,nblockbd,npw_k
    integer :: ii
    real(dp) :: ar
    ! Arrays
    integer,allocatable :: kg_k(:,:)
    real(dp),allocatable :: eknk(:),ek_k(:),kinpw(:)

    ! *********************************************************************

    bdtot_index=0

    ABI_ALLOCATE(eknk,(mband*nkpt*nsppol))
    do isppol=1,nsppol
      do ikpt=1,nkpt
        npw_k = npwarr(ikpt)
        nband_k=nband(ikpt+(isppol-1)*nkpt)

        ABI_ALLOCATE(kg_k,(3,npw_k))
        ABI_ALLOCATE(kinpw,(npw_k))
        ABI_ALLOCATE(ek_k,(nband_k))

        kg_k(:,1:npw_k)=kg(:,1:npw_k)
        nblockbd=nband_k/(mpi_enreg%nproc_band*mpi_enreg%bandpp)
        blocksize=nband_k/nblockbd

        call mkkin(ecut,ecutsm,effmass_free,gmet,kg_k,kinpw,kptns(:,ikpt),npw_k,0,0)

        do iblock=1,nblockbd
          do iblocksize=1,blocksize
            iband=(iblock-1)*blocksize+iblocksize
            call meanvalue_g(ar,kinpw,0,istwfk(ikpt),mpi_enreg,npw_k,my_nspinor,&
            & cg(:,1+(iband-1)*npw_k*my_nspinor:iband*npw_k*my_nspinor),&
            & cg(:,1+(iband-1)*npw_k*my_nspinor:iband*npw_k*my_nspinor),0)

            ek_k(iband)=ar
          end do
        end do

        eknk (1+bdtot_index : nband_k+bdtot_index) = ek_k (:)
        ABI_DEALLOCATE(ek_k)
        ABI_DEALLOCATE(kinpw)
        ABI_DEALLOCATE(kg_k)

        bdtot_index=bdtot_index+nband_k
      end do
    end do

    call hightemp%compute_e_shiftfactor(eigen,eknk,mband,mpi_enreg,&
    & nkpt,nsppol)

    ABI_DEALLOCATE(eknk)
  end subroutine hightemp_get_e_shiftfactor

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
