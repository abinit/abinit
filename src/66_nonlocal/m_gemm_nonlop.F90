!!****m* ABINIT/m_gemm_nonlop
!! NAME
!! m_gemm_nonlop
!!
!! FUNCTION
!!  This module provides functions to compute the nonlocal operator by means of the BLAS GEMM
!!  routine. By treating ndat simultaneous wavefunctions, it is able to exploit BLAS3 routines,
!!  which leads to excellent CPU efficiency and OpenMP scalability.
!!
!! COPYRIGHT
!! Copyright (C) 2014-2020 ABINIT group (AL)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

! TODO list :
! Don't allocate the full nkpt structures, only those that are treated by this proc: use same init as in m_bandfft_kpt
! support more options (forces & stresses mostly)
! Support RF/other computations (only GS right now)
! handle the case where nloalg(2) < 0, ie no precomputation of ph3d
! more systematic checking of the workflow (right now, only works if init/make/gemm/destroy, no multiple makes, etc)
! Avoid allocating the complex matrix when istwfk > 1
! Merge with chebfi's invovl


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gemm_nonlop

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_abi_linalg

 use defs_abitypes, only : MPI_type
 use m_opernlc_ylm, only : opernlc_ylm
 use m_pawcprj, only : pawcprj_type

 implicit none

 private

 ! Use these routines in order: first call init, then call make_gemm_nonlop for each k point,
 ! then call gemm_nonlop to do the actual computation, and call destroy when done. See gstate and vtorho.
 public :: init_gemm_nonlop
 public :: make_gemm_nonlop
 public :: gemm_nonlop
 public :: destroy_gemm_nonlop
!!***

!!****t* m_gemm_nonlop/gemm_nonlop_type
!! NAME
!! gemm_nonlop_type
!!
!! FUNCTION
!! Contains information needed to apply the nonlocal operator
!!
!! SOURCE
 type,public :: gemm_nonlop_type

   integer :: nprojs

   real(dp), allocatable :: projs(:, :, :)
   ! (2, npw, nprojs)
   real(dp), allocatable :: projs_r(:, :, :)
   ! (1, npw, nprojs)
   real(dp), allocatable :: projs_i(:, :, :)
   ! (1, npw, nprojs)
 end type gemm_nonlop_type
!!***

 type(gemm_nonlop_type), save, public, allocatable :: gemm_nonlop_kpt(:)
 !(nkpt)

 integer, save, public :: gemm_nonlop_ikpt_this_proc_being_treated
 !! This is oh so very crude, but I can't find any other way to do it without passing ikpt deep down to nonlop
 logical, save, public :: gemm_nonlop_use_gemm = .false.
 ! Public variable indicating whether we should call gemm_nonlop or fall back to the usual nonlop. Set to false
 ! in order not to interfere with non-GS calls to nonlop.

contains

!!****f* m_gemm_nonlop/init_gemm_nonlop
!! NAME
!! init_gemm_nonlop
!!
!! FUNCTION
!! Initalization of the gemm_nonlop_kpt array
!!
!! INPUTS
!! nkpt= number of k-points
!!
!! PARENTS
!!      m_gstate
!!
!! CHILDREN
!!      dgemm,opernlc_ylm,xmpi_sum,zgemm
!!
!! SOURCE
 subroutine init_gemm_nonlop(nkpt)

  integer,intent(in) :: nkpt
  integer :: ikpt

! *************************************************************************

  ! TODO only allocate the number of kpt treated by this proc
  ABI_MALLOC(gemm_nonlop_kpt, (nkpt))
  do ikpt=1,nkpt
    gemm_nonlop_kpt(ikpt)%nprojs = -1
  end do

 end subroutine init_gemm_nonlop
!!***

!!****f* m_gemm_nonlop/destroy_gemm_nonlop
!! NAME
!! destroy_gemm_nonlop
!!
!! FUNCTION
!! Initalization of the gemm_nonlop_kpt array
!!
!! INPUTS
!! nkpt= number of k-points
!!
!! PARENTS
!!      m_gstate
!!
!! CHILDREN
!!      dgemm,opernlc_ylm,xmpi_sum,zgemm
!!
!! SOURCE
 subroutine destroy_gemm_nonlop(nkpt)

  integer,intent(in) :: nkpt
  integer :: ikpt

! *************************************************************************

  ! TODO add cycling if kpt parallelism
  do ikpt = 1,nkpt
    if(gemm_nonlop_kpt(ikpt)%nprojs /= -1) then
      ABI_FREE(gemm_nonlop_kpt(ikpt)%projs)
      ABI_FREE(gemm_nonlop_kpt(ikpt)%projs_r)
      ABI_FREE(gemm_nonlop_kpt(ikpt)%projs_i)
      gemm_nonlop_kpt(ikpt)%nprojs = -1
    end if
  end do

  ABI_FREE(gemm_nonlop_kpt)

 end subroutine destroy_gemm_nonlop
!!***


!!****f* m_gemm_nonlop/make_gemm_nonlop
!! NAME
!! make_gemm_nonlop
!!
!! FUNCTION
!! Build the gemm_nonlop array
!!
!! INPUTS
!!
!! PARENTS
!!      m_dft_energy,m_vtorho
!!
!! CHILDREN
!!      dgemm,opernlc_ylm,xmpi_sum,zgemm
!!
!! SOURCE

 subroutine make_gemm_nonlop(ikpt,npw,lmnmax,ntypat,indlmn,nattyp,istwf_k,ucvol,ffnl_k,ph3d_k)

  integer, intent(in) :: ikpt
  integer, intent(in) :: npw, lmnmax,ntypat
  integer, intent(in) :: indlmn(:,:,:)
  integer, intent(in) :: nattyp(ntypat)
  integer, intent(in) :: istwf_k
  real(dp), intent(in) :: ucvol
  real(dp), intent(in) :: ffnl_k(:,:,:,:)
  real(dp), intent(in) :: ph3d_k(:,:,:)

  integer :: nprojs

  real(dp),allocatable :: atom_projs(:,:,:), temp(:)
  integer :: itypat, ilmn, nlmn, ia, iaph3d, shift
  integer :: il, ipw
  logical :: parity

! *************************************************************************

  iaph3d = 1

  ABI_MALLOC(atom_projs, (2, npw, lmnmax))
  ABI_MALLOC(temp, (npw))

  if(gemm_nonlop_kpt(ikpt)%nprojs /= -1) then
    ! We have been here before, cleanup before remaking
    ABI_FREE(gemm_nonlop_kpt(ikpt)%projs)
    ABI_FREE(gemm_nonlop_kpt(ikpt)%projs_r)
    ABI_FREE(gemm_nonlop_kpt(ikpt)%projs_i)
    gemm_nonlop_kpt(ikpt)%nprojs = -1
  end if

  ! build nprojs
  nprojs = 0
  do itypat=1,ntypat
    nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
  end do

  gemm_nonlop_kpt(ikpt)%nprojs = nprojs

  ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs, (2, npw, gemm_nonlop_kpt(ikpt)%nprojs))
  gemm_nonlop_kpt(ikpt)%projs = zero
  if(istwf_k > 1) then
    ! We still allocate the complex matrix, in case we need it for spinors. TODO could be avoided
    ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_r, (1, npw, gemm_nonlop_kpt(ikpt)%nprojs))
    ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_i, (1, npw, gemm_nonlop_kpt(ikpt)%nprojs))
    gemm_nonlop_kpt(ikpt)%projs_r = zero
    gemm_nonlop_kpt(ikpt)%projs_i = zero
  else
    ! Still allocate so we can deallocate it in destroy_gemm_nonlop
    ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_r, (1, 1, 1))
    ABI_MALLOC(gemm_nonlop_kpt(ikpt)%projs_i, (1, 1, 1))
  end if

  shift = 0
  do itypat = 1, ntypat
    nlmn = count(indlmn(3,:,itypat)>0)

    do ia = 1, nattyp(itypat)

      !! build atom_projs, from opernlb
      !! P = 4pi/sqrt(ucvol)* conj(diag(ph3d)) * ffnl * diag(parity), with parity = (-i)^l
      atom_projs(:,:,:) = zero

      ! start from 4pi/sqrt(ucvol)*ffnl
      ! atom_projs(1, :, 1:nlmn) = four_pi/sqrt(ham%ucvol) * ham%ffnl_k(:, 1, 1:nlmn)
      ! TODO vectorize (DCOPY with stride)
      do ipw=1, npw
        atom_projs(1,ipw, 1:nlmn) = four_pi/sqrt(ucvol) * ffnl_k(ipw, 1, 1:nlmn, itypat)
      end do

      ! multiply by (-i)^l
      do ilmn=1,nlmn
        il=mod(indlmn(1,ilmn, itypat),4);
        parity=(mod(il,2)==0)
        if (il>1) then
          ! multiply by -1
          atom_projs(:,:,ilmn) = -atom_projs(:,:,ilmn)
        end if
        if(.not. parity) then
          ! multiply by -i
          temp = atom_projs(2,:,ilmn)
          atom_projs(2,:,ilmn) = -atom_projs(1,:,ilmn)
          atom_projs(1,:,ilmn) =  temp
        end if
      end do

      ! multiply by conj(ph3d)
      do ilmn=1,nlmn
        temp = atom_projs(1, :, ilmn)
        atom_projs(1, :, ilmn) = atom_projs(1, :, ilmn) * ph3d_k(1, :, iaph3d) &
&                              + atom_projs(2, :, ilmn) * ph3d_k(2, :, iaph3d)
        atom_projs(2, :, ilmn) = atom_projs(2, :, ilmn) * ph3d_k(1, :, iaph3d) &
&                              - temp                   * ph3d_k(2, :, iaph3d)
      end do

      !! atom_projs is built, copy to projs
      gemm_nonlop_kpt(ikpt)%projs(:, :, shift+1:shift+nlmn) = atom_projs(:, :, 1:nlmn)
      if(istwf_k > 1) then
        gemm_nonlop_kpt(ikpt)%projs_r(1, :, shift+1:shift+nlmn) = atom_projs(1, :, 1:nlmn)
        gemm_nonlop_kpt(ikpt)%projs_i(1, :, shift+1:shift+nlmn) = atom_projs(2, :, 1:nlmn)
      end if
      shift = shift + nlmn
      iaph3d = iaph3d + 1
    end do
  end do

  ABI_FREE(atom_projs)
  ABI_FREE(temp)

 end subroutine make_gemm_nonlop
!!***

!!****f* m_gemm_nonlop/gemm_nonlop
!! NAME
!! gemm_nonlop
!!
!! FUNCTION
!! Replacement of nonlop. same prototype as nonlop although not all options are implemented.
!!
!! INPUTS
!!
!! PARENTS
!!      m_nonlop
!!
!! CHILDREN
!!      dgemm,opernlc_ylm,xmpi_sum,zgemm
!!
!! SOURCE
 subroutine gemm_nonlop(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
&                 enl,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&
&                 kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                 mpi_enreg,mpsang,mpssoang,natom,nattyp,ndat,ngfft,nkpgin,nkpgout,nloalg,&
&                 nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,phkxredin,&
&                 phkxredout,ph1d,ph3din,ph3dout,signs,sij,svectout,&
&                 tim_nonlop,ucvol,useylm,vectin,vectout,&
&                 use_gpu_cuda)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir
  integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin
  integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO
  integer,intent(in) :: paw_opt,signs,tim_nonlop,useylm
  integer,optional,intent(in) :: use_gpu_cuda
  real(dp),intent(in) :: lambda(ndat),ucvol
  type(MPI_type),intent(in) :: mpi_enreg
  !arrays
  integer,intent(in) :: atindx1(natom),indlmn(6,lmnmax,ntypat),kgin(3,npwin)
  integer,intent(in) :: kgout(3,npwout),nattyp(ntypat),ngfft(18),nloalg(3)
  real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinortot**2,dimekbq)
  real(dp),intent(in) :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
  real(dp),intent(in) :: ffnlout(npwout,dimffnlout,lmnmax,ntypat),gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3),kpgin(npwin,nkpgin*useylm)
  real(dp),intent(in) :: kpgout(npwout,nkpgout*useylm),kptin(3),kptout(3)
  real(dp),intent(in) :: phkxredin(2,natom),phkxredout(2,natom),ph1d(2,*)
  real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
  real(dp),intent(inout) :: ph3din(2,npwin,matblk),ph3dout(2,npwout,matblk)
  real(dp),intent(inout) :: vectin(2,npwin*nspinor*ndat)
  real(dp),intent(inout) :: enlout(nnlout*ndat)
  real(dp),intent(out) :: svectout(2,npwout*nspinor*(paw_opt/3)*ndat)
  real(dp),intent(inout) :: vectout(2,npwout*nspinor*ndat) !vz_i
  type(pawcprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5)*ndat)

  ! locals
  integer :: idat, nprojs, shift, iatom, nlmn, ierr, ibeg, iend
  integer :: cplex, cplex_enl, cplex_fac
  integer :: iatm, ndgxdt, ndgxdtfac, nd2gxdt, nd2gxdtfac, optder, itypat, ilmn
  integer :: cplex_dgxdt(1), cplex_d2gxdt(1)
  real(dp) :: dgxdt_dum_in(1,1,1,1,1), dgxdt_dum_out(1,1,1,1,1),dgxdt_dum_out2(1,1,1,1,1)
  real(dp) :: d2gxdt_dum_in(1,1,1,1,1), d2gxdt_dum_out(1,1,1,1,1),d2gxdt_dum_out2(1,1,1,1,1)
  real(dp), allocatable :: sij_typ(:)
  real(dp), allocatable :: projections(:,:,:), s_projections(:,:,:), vnl_projections(:,:,:)
  real(dp), allocatable :: temp_realvec(:)

! *************************************************************************

  ! We keep the same interface as nonlop, but we don't use many of those
  ABI_UNUSED((/ffnlin,ffnlout,gmet,gprimd,kpgin,kpgout/))
  ABI_UNUSED((/ph1d(1,1),ph3din,ph3dout/))
  ABI_UNUSED((/phkxredin,phkxredout,ucvol/))
  ABI_UNUSED((/mgfft,mpsang,mpssoang/))
  ABI_UNUSED((/kptin,kptout/))
  ABI_UNUSED((/idir,nloalg,ngfft,kgin,kgout,ngfft,only_SO,tim_nonlop,use_gpu_cuda,signs/))
  ABI_UNUSED((/enlout/))

  cplex=2;if (istwf_k>1) cplex=1
  cplex_enl=1;if (paw_opt>0) cplex_enl=2*dimenl1/(lmnmax*(lmnmax+1)) ! is enl complex?
  cplex_fac=max(cplex,dimekbq)
  if ((nspinortot==2.or.cplex_enl==2).and.paw_opt>0.and.choice/=7) cplex_fac=2 ! is vnl_projections complex?

  nprojs = gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%nprojs

  ! These will store the non-local factors for vectin, svectout and vectout respectively
  ABI_MALLOC(projections,(cplex, nprojs,nspinor*ndat))
  ABI_MALLOC(s_projections,(cplex, nprojs,nspinor*ndat))
  ABI_MALLOC(vnl_projections,(cplex_fac, nprojs,nspinor*ndat))
  projections = zero
  s_projections = zero
  vnl_projections = zero

  if(nprojs == 0) then
    ! TODO check if this is correct
    vectout = zero
    if(paw_opt>0) svectout = vectin
    return
  end if

  if(cpopt >= 2) then
    ! retrieve from cprjin
    do idat=1, ndat*nspinor
      shift = 0
      do iatom = 1, natom
        nlmn = cprjin(iatom, idat)%nlmn
        projections(1:cplex, shift+1:shift+nlmn, idat) = cprjin(iatom, idat)%cp(1:cplex, 1:nlmn)
        shift = shift + nlmn
      end do
    end do
  else
    ! opernla
    if(cplex == 2) then
      call abi_zgemm_2r('C', 'N', nprojs, ndat*nspinor, npwin, cone, &
&                gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwin,&
&                vectin, npwin, czero, projections, nprojs)
    else
       ABI_MALLOC(temp_realvec,(MAX(npwout,npwin)*nspinor*ndat))
      ! only compute real part of projections = P^* psi => projections_r = P_r^T psi_r + P_i^T psi_i
      temp_realvec(1:npwin*nspinor*ndat) = vectin(1,1:npwin*nspinor*ndat)
      if(istwf_k == 2 .and. mpi_enreg%me_g0 == 1) then
        do idat=1, ndat*nspinor
          temp_realvec(1+(idat-1)*npwin) = temp_realvec(1+(idat-1)*npwin)/2
        end do
      end if
      call DGEMM('T', 'N', nprojs, ndat*nspinor, npwin, one, &
&                gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwin, &
&                temp_realvec, npwin, zero, projections, nprojs)
      temp_realvec(1:npwin*nspinor*ndat) = vectin(2,1:npwin*nspinor*ndat)
      if(istwf_k == 2 .and. mpi_enreg%me_g0 == 1) then
        do idat=1, ndat*nspinor
          temp_realvec(1+(idat-1)*npwin) = zero
        end do
      end if
      call DGEMM('T', 'N', nprojs, ndat*nspinor, npwin, one, &
&                gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_i, npwin, &
&                temp_realvec, npwin, one , projections, nprojs)
      projections = projections * 2
      ABI_FREE(temp_realvec)
    end if
    call xmpi_sum(projections,mpi_enreg%comm_fft,ierr)

    if(cpopt >= 0) then
      ! store in cprjin
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = 1, natom
          nlmn = cprjin(iatom, idat)%nlmn
          cprjin(iatom, idat)%cp(1:cplex, 1:nlmn) = projections(1:cplex, shift+1:shift+nlmn, idat)
          shift = shift + nlmn
        end do
      end do
    end if
  end if

  if(choice > 0) then

    if(choice /= 7) then
      ! opernlc
      iatm = 0
      ndgxdt = 0
      ndgxdtfac = 0
      nd2gxdt = 0
      nd2gxdtfac = 0
      optder = 0

      shift = 0
      ABI_MALLOC(sij_typ,(((paw_opt+1)/3)*lmnmax*(lmnmax+1)/2))
      do itypat=1, ntypat
        nlmn=count(indlmn(3,:,itypat)>0)
        if (paw_opt>=2) then
          if (cplex_enl==1) then
            do ilmn=1,nlmn*(nlmn+1)/2
              sij_typ(ilmn)=sij(ilmn,itypat)
            end do
          else
            do ilmn=1,nlmn*(nlmn+1)/2
              sij_typ(ilmn)=sij(2*ilmn-1,itypat)
            end do
          end if
        end if

        ibeg = shift+1
        iend = shift+nattyp(itypat)*nlmn

        do idat = 1,ndat
          call opernlc_ylm(atindx1,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_enl,cplex_fac,dgxdt_dum_in,dgxdt_dum_out,dgxdt_dum_out2,&
&         d2gxdt_dum_in,d2gxdt_dum_out,d2gxdt_dum_out2,dimenl1,dimenl2,dimekbq,enl,&
&         projections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
&         vnl_projections(:, ibeg:iend,1+nspinor*(idat-1):nspinor*idat),&
&         s_projections(:, ibeg:iend,1+nspinor*(idat-1):nspinor*idat),&
&         iatm,indlmn(:,:,itypat),itypat,lambda(idat),mpi_enreg,natom,ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,&
&         nattyp(itypat),nlmn,nspinor,nspinortot,optder,paw_opt,sij_typ)
        end do

        shift = shift + nattyp(itypat)*nlmn
        iatm = iatm+nattyp(itypat)
      end do
      ABI_FREE(sij_typ)
    else
      s_projections = projections
    end if

    ! opernlb
    if(paw_opt == 3 .or. paw_opt == 4) then
      ! Get svectout from s_projections
      if(cplex == 2) then
        call abi_zgemm_2r('N', 'N', npwout, ndat*nspinor, nprojs, cone, &
&                      gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwout, &
&                      s_projections, nprojs, czero, svectout, npwout)
      else
        ABI_MALLOC(temp_realvec,(MAX(npwout,npwin)*nspinor*ndat))
        call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs, one, &
&                  gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwout, &
&                  s_projections, nprojs, zero, temp_realvec, npwout)
        svectout(1,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
        call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs, one, &
&                  gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_i, npwout,&
&                  s_projections, nprojs, zero, temp_realvec, npwout)
        svectout(2,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
        ABI_FREE(temp_realvec)
      end if
      if(choice /= 7) svectout = svectout + vectin ! TODO understand this
    end if
    if(paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4) then
      ! Get vectout from vnl_projections
      if(cplex_fac == 2) then
        call abi_zgemm_2r('N', 'N', npwout, ndat*nspinor, nprojs, cone, &
&                      gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs, npwout, &
&                      vnl_projections, nprojs, czero, vectout, npwout)
      else
        ABI_MALLOC(temp_realvec,(MAX(npwout,npwin)*nspinor*ndat))
        call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs, one, &
&                  gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_r, npwout, &
&                  vnl_projections, nprojs, zero, temp_realvec, npwout)
        vectout(1,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
        call DGEMM('N', 'N', npwout, ndat*nspinor, nprojs, one, &
&                  gemm_nonlop_kpt(gemm_nonlop_ikpt_this_proc_being_treated)%projs_i, npwout, &
&                  vnl_projections, nprojs, zero, temp_realvec, npwout)
        vectout(2,1:npwout*nspinor*ndat) = temp_realvec(1:npwout*nspinor*ndat)
        ABI_FREE(temp_realvec)
      end if
    end if
  end if

  ABI_FREE(projections)
  ABI_FREE(s_projections)
  ABI_FREE(vnl_projections)
 end subroutine gemm_nonlop
!***

end module m_gemm_nonlop
!!***
