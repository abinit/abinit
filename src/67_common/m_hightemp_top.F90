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
  subroutine hightemp_prt_cg(ecut,exchn2n3d,istwfk,kpt,&
  & mpi_enreg,mpw,nband,nkpt,npwarr,nsppol,rprimd)
    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: mpw,nkpt,nsppol,exchn2n3d
    real(dp),intent(in) :: ecut
    type(MPI_type),intent(inout) :: mpi_enreg
    ! Arrays
    integer,intent(in) :: istwfk(nkpt),nband(nkpt),npwarr(nkpt)
    real(dp),intent(in) :: kpt(3,nkpt),rprimd(3,3)

    ! Local variables -------------------------
    ! Scalars
    integer :: cgshift,ckpt,iband,ioffkg,iout,ipw
    integer :: incrdegcg,ikpt,ipwbis,krow,mkmem,npw_k
    real(dp) :: tempcgk,energmin,cgsum,tot_cgsum,ucvol
    character(len=20) :: filenameoutpw
    character(len=4) :: mode_paral
    ! Arrays
    real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
    integer,allocatable :: kg(:,:),kg_k(:,:),npwarr1(:),npwtot1(:)
    real(dp),allocatable :: buf(:),cg_k(:,:),kpgnorm(:),tempwfk(:,:)

    mode_paral='PERS'
    iout=-1
    call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)
    ABI_ALLOCATE(kg,(3,mpw*nkpt))

    do ckpt=1,nkpt
      ! PENSER A LA BOUCLE SUR LES KPOINTS

      npw_k=npwarr(ckpt)
      ABI_ALLOCATE(kg_k,(3,npw_k))
      ABI_ALLOCATE(kpgnorm,(npw_k))
      ABI_ALLOCATE(npwarr1,(nkpt))
      ABI_ALLOCATE(npwtot1,(nkpt))
      mkmem=nkpt

      ! S'OCCUPER DE CA !!!
      !    Create positions index for pw
      call kpgio(ecut,exchn2n3d,gmet,istwfk,kg,kpt,mkmem,nband,nkpt,&
      & mode_paral,mpi_enreg,mpw,npwarr1,npwtot1,nsppol)

      ioffkg=0
      do ikpt=1,ckpt-1
        ioffkg=ioffkg+npwarr1(mkmem)
      end do
      kg_k(:,1:npw_k)=kg(:,1+ioffkg:npw_k+ioffkg)
      call getkpgnorm(gprimd,kpt(:,ckpt),kg_k,kpgnorm,npw_k)


      ABI_ALLOCATE(tempwfk,(mpw, 3))
      ABI_ALLOCATE(buf,(2))
      write(std_out,*) nband(ckpt)
      do iband=1, nband(ckpt)
        tempwfk(:,:) = 0_DP
        cgshift=(iband-1)*npw_k
        do ipw=1, npw_k
          tempwfk(ipw,2) = ABS(dcmplx(cg_k(1,cgshift+ipw), cg_k(2,cgshift+ipw)))
          ! Hartree to eV * Kinetic energy of the pw
          tempwfk(ipw,1) = 27.2114*(2*pi*kpgnorm(ipw))**2/2.
          tempwfk(ipw,3) = ipw
        end do

        ! Sorting in energies
        do ipw=1, npw_k
          krow = minloc(tempwfk(ipw:npw_k, 1), dim=1) + ipw - 1
          buf(:) = tempwfk(ipw, :)
          tempwfk(ipw, :) = tempwfk(krow, :)
          tempwfk(krow, :) = buf(:)
        end do

        ! Normalization constant
        tot_cgsum = sum(tempwfk(:,2))

        write(filenameoutpw, '(A,I4.4,A,I4.4)') 'outpw_k', ckpt, '_b', iband
        open(file=filenameoutpw, unit=23)

        ! do ipw = 1, npw_k
        !    write(23,*) tempwfk(ipw-1, 1), tempwfk(ipw,2), 1, tempwfk(ipw-1, 3)
        ! end do

        ! Degenerescence in energies
        incrdegcg = 1
        energmin = tempwfk(1,1)
        cgsum = tempwfk(1,2)
        do ipw=2, npw_k
          if (abs(tempwfk(ipw,1)-energmin) .LE. 1E-10) then
             cgsum = cgsum + tempwfk(ipw,2)
             incrdegcg = incrdegcg + 1
          else
             write(23,*) tempwfk(ipw-1, 1), cgsum/tot_cgsum, incrdegcg, tempwfk(ipw-1, 3)
             energmin = tempwfk(ipw, 1)
             cgsum = tempwfk(ipw,2)
             incrdegcg = 1
          end if
        end do
        close(23)
      end do
    end do
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

end module m_hightemp_top
