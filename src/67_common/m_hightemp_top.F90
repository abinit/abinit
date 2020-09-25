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
  use m_geometry,       only : metric
  use m_gsphere,        only : getkpgnorm
  use m_kg,             only : kpgio

  implicit none

  public :: hightemp_prt_cg,hightemp_prt_cprj
  !!***

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
  subroutine hightemp_prt_cg(cg,ckpt,ecut,eig_k,ek_k,exchn2n3d,fnameabo,istwfk,kg_k,kpt,&
  & mcg,mpi_enreg,mpw,nband,nkpt,npw_k,nsppol,rprimd)
    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: ckpt,mcg,mpw,nkpt,npw_k,nsppol,exchn2n3d
    real(dp),intent(in) :: ecut
    type(MPI_type),intent(inout) :: mpi_enreg
    character(len=*),intent(in) :: fnameabo
    ! Arrays
    integer,intent(in) :: istwfk(nkpt),kg_k(3,npw_k),nband(nkpt)
    real(dp),intent(in) :: kpt(3,nkpt),rprimd(3,3)
    real(dp),intent(in) :: cg(2,mcg)
    real(dp),intent(in) :: eig_k(nband(ckpt)),ek_k(nband(ckpt))

    ! Local variables -------------------------
    ! Scalars
    integer :: cgshift,cgshift_tot,iband,iout,ipw,ipw_tot
    integer :: mcg_tot,mpierr,npw_tot
    real(dp) :: ucvol,kpg1_tot,kpg2_tot,kpg3_tot
    character(len=50) :: filenameoutpw
    character(len=4) :: mode_paral
    ! Arrays
    integer :: nppw_tot(mpi_enreg%nproc_band)
    real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
    integer,allocatable :: kg_tot(:,:)
    real(dp),allocatable :: cg_tot(:,:),kpgnorm(:),kpgnorm_tot(:),tempwfk(:,:)

    ! *********************************************************************

    write(std_out,'(a,I5.5,a)') 'Writing plane waves coefs of kpt=',ckpt,'...'
    mode_paral='PERS'
    iout=-1
    call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

    ABI_ALLOCATE(kpgnorm,(npw_k))
    call getkpgnorm(gprimd,kpt(:,ckpt),kg_k,kpgnorm,npw_k)

    mcg_tot=mcg
    nppw_tot(:)=0
    nppw_tot(mpi_enreg%me_band + 1)=npw_k
    call xmpi_sum(nppw_tot,mpi_enreg%comm_world,mpierr)
    call xmpi_sum(mcg_tot,mpi_enreg%comm_world,mpierr)

    npw_tot=sum(nppw_tot)
    ABI_ALLOCATE(cg_tot,(2,mcg_tot))
    ABI_ALLOCATE(kpgnorm_tot,(npw_tot))
    ABI_ALLOCATE(kg_tot,(3,npw_tot))

    cg_tot(:,:)=zero
    kpgnorm_tot(:)=zero
    kg_tot(:,:)=0
    do ipw=1,npw_k
      ipw_tot=sum(nppw_tot(1:mpi_enreg%me_band))+ipw
      kpgnorm_tot(ipw_tot)=kpgnorm(ipw)
      kg_tot(:,ipw_tot)=kg_k(:,ipw)

      do iband=1,nband(ckpt)
        cgshift=(iband-1)*npw_k
        cgshift_tot=(iband-1)*npw_tot
        cg_tot(:,cgshift_tot+ipw_tot)=cg(:,cgshift+ipw)
      end do
    end do

    call xmpi_sum(cg_tot,mpi_enreg%comm_world,mpierr)
    call xmpi_sum(kpgnorm_tot,mpi_enreg%comm_world,mpierr)
    call xmpi_sum(kg_tot,mpi_enreg%comm_world,mpierr)


    if(mpi_enreg%me==mpi_enreg%me_kpt) then

      ! Writting plane waves vector coordinates
      write(filenameoutpw,'(A,I5.5)') '_PW_MESH_k',ckpt

      open(file=trim(fnameabo)//trim(filenameoutpw), unit=23)
      do ipw=1, npw_tot
        kpg1_tot=kpt(1,ckpt)+dble(kg_tot(1,ipw))
        kpg2_tot=kpt(2,ckpt)+dble(kg_tot(2,ipw))
        kpg3_tot=kpt(3,ckpt)+dble(kg_tot(3,ipw))

        write(23,'(i14,ES13.5,ES13.5,ES13.5,ES13.5)')&
        & ipw,&
        & gprimd(1,1)*kpg1_tot+gprimd(1,2)*kpg2_tot+gprimd(1,3)*kpg3_tot,&
        & gprimd(2,1)*kpg1_tot+gprimd(2,2)*kpg2_tot+gprimd(2,3)*kpg3_tot,&
        & gprimd(3,1)*kpg1_tot+gprimd(3,2)*kpg2_tot+gprimd(3,3)*kpg3_tot,&
        & kpgnorm_tot(ipw)
      end do
      close(23)

      ! Writting Eigen energies
      write(filenameoutpw,'(A,I5.5)') '_PW_EIG_k',ckpt
      open(file=trim(fnameabo)//trim(filenameoutpw), unit=23)
      write(23,'(ES12.5)') eig_k
      close(23)

      ! Writting Kinetic energies
      write(filenameoutpw,'(A,I5.5)') '_PW_KIN_k',ckpt
      open(file=trim(fnameabo)//trim(filenameoutpw), unit=23)
      write(23,'(ES12.5)') ek_k
      close(23)

      ABI_ALLOCATE(tempwfk,(npw_tot, 3))
      do iband=1, nband(ckpt)

        tempwfk(:,:) = zero
        cgshift_tot=(iband-1)*npw_tot
        do ipw=1, npw_tot
          tempwfk(ipw,2) = ABS(dcmplx(cg_tot(1,cgshift_tot+ipw), cg_tot(2,cgshift_tot+ipw)))
          ! Hartree to eV * Kinetic energy of the pw
          tempwfk(ipw,1) = 27.2114*(2*pi*kpgnorm_tot(ipw))**2/2.
          tempwfk(ipw,3) = ipw
        end do

        write(filenameoutpw, '(A,I5.5,A,I5.5)') '_PW_k',ckpt,'_b',iband
        open(file=trim(fnameabo)//trim(filenameoutpw), unit=23)
        do ipw=1, npw_tot
          write(23,'(i14,ES14.6,ES14.6,i14)') int(tempwfk(ipw, 3)),tempwfk(ipw, 1),tempwfk(ipw,2)
        end do

        close(23)
      end do
      ABI_DEALLOCATE(tempwfk)
    end if
    ABI_DEALLOCATE(kpgnorm)
    ABI_DEALLOCATE(kg_tot)
    ABI_DEALLOCATE(cg_tot)
    ABI_DEALLOCATE(kpgnorm_tot)
  end subroutine hightemp_prt_cg
  !!***

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
  subroutine hightemp_prt_cprj(cprj,eigen,gs_hamk,istep,mband,&
  & mcprj,mpi_enreg,natom,nkpt,nsppol,occ)
    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: istep,mband,mcprj,natom,nkpt,nsppol
    type(gs_hamiltonian_type),intent(inout) :: gs_hamk
    type(MPI_type),intent(inout) :: mpi_enreg
    ! Arrays
    real(dp),intent(inout) :: eigen(mband*nkpt*nsppol)
    real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
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
  !!***

end module m_hightemp_top
!!***
