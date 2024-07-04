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
!! Copyright (C) 2014-2024 ABINIT group (AL)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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
 use m_xomp
 use m_abi_linalg
 use m_gemm_nonlop_projectors

 use defs_abitypes, only : MPI_type
 use m_opernlc_ylm, only : opernlc_ylm
 use m_opernla_gemm, only : opernla_gemm
 use m_opernlb_gemm, only : opernlb_gemm
 use m_opernld_ylm_allwf, only : opernld_ylm_allwf
 use m_opernld_ylm, only : opernld_ylm
 use m_pawcprj, only : pawcprj_type
 use m_geometry, only : strconv
 use m_kg, only : mkkpg
 use m_hamiltonian, only : KPRIME_H_K, K_H_KPRIME, K_H_K, KPRIME_H_KPRIME


#ifdef HAVE_FC_ISO_C_BINDING
 use, intrinsic :: iso_c_binding, only : c_int32_t, c_int64_t, c_float, c_double, c_size_t, c_loc, c_ptr
#endif


 implicit none

 private

 public :: gemm_nonlop


!!***

!----------------------------------------------------------------------


contains

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop/gemm_nonlop
!! NAME
!! gemm_nonlop
!!
!! FUNCTION
!! Replacement of nonlop. same prototype as nonlop although not all options are implemented.
!!
!! INPUTS
!! [gpu_option] = GPU implementation to use, i.e. cuda, openMP, ... (0=not using GPU)
!!
!! SOURCE
 subroutine gemm_nonlop(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
&                 enl,enl_ndat,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&
&                 kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                 mpi_enreg,natom,nattyp,ndat,ngfft,nkpgin,nkpgout,nloalg,&
&                 nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,&
&                 ph3din,ph3dout,signs,sij,svectout,&
&                 tim_nonlop,ucvol,useylm,vectin,vectout,&
&                 atom_proj_shift,select_k,iatom_only,typat,usepaw,&
&                 vectproj,gpu_option)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir
  integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,natom,ndat,nkpgin
  integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO
  integer,intent(in) :: paw_opt,signs,tim_nonlop,useylm,atom_proj_shift,select_k,iatom_only,usepaw
  integer,optional,intent(in) :: gpu_option
  real(dp),intent(in) :: lambda(ndat),ucvol
  type(MPI_type),intent(in) :: mpi_enreg
  !arrays
  integer,intent(in),target :: atindx1(natom),indlmn(6,lmnmax,ntypat),kgin(3,npwin)
  integer,intent(in),target :: kgout(3,npwout),nattyp(ntypat),ngfft(18),nloalg(3),typat(natom)
  real(dp),intent(in),target :: enl(:,:,:,:),enl_ndat(:,:,:,:,:)
  real(dp),intent(in),target :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
  real(dp),intent(in),target :: ffnlout(npwout,dimffnlout,lmnmax,ntypat),gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3),kptin(3),kptout(3)
  real(dp),intent(in),target :: kpgin(npwin,nkpgin*useylm),kpgout(npwout,nkpgout*useylm)
  real(dp),intent(in),target :: sij(dimenl1,ntypat*((paw_opt+1)/3))
  real(dp),intent(inout),target :: ph3din(2,npwin,matblk),ph3dout(2,npwout,matblk)
  real(dp),intent(inout),target :: vectin(2,npwin*nspinor*ndat)
  real(dp),intent(inout) :: enlout(nnlout*ndat)
  real(dp),intent(out),target :: svectout(:,:)
  real(dp),intent(inout),target :: vectout(:,:)
  real(dp),intent(inout),optional, ABI_CONTIGUOUS target :: vectproj(:,:,:)
  type(pawcprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5)*ndat)

  ! locals
  integer :: ii, idat, igrad, nprojs, ngrads, ngrads2, shift, iatom, nlmn, ierr, ibeg, iend, ikin, ikout
  integer :: cplex, cplex_enl, cplex_fac
  integer :: nnlout_test
  integer :: iatm, ndgxdt, ndgxdtfac, nd2gxdt, nd2gxdtfac, optder, itypat, ilmn
  integer :: cplex_dgxdt(9), cplex_d2gxdt(18)
  logical :: local_vectproj,use_enl_ndat
  real(dp) :: d2gxdt_dum_in(1,1,1,1,1), d2gxdt_dum_out(1,1,1,1,1),d2gxdt_dum_out2(1,1,1,1,1)
  real(dp), allocatable :: sij_typ(:)
  real(dp), ABI_CONTIGUOUS pointer :: projections(:,:,:)
  real(dp), allocatable :: s_projections(:,:,:), vnl_projections(:,:,:)
  real(dp), allocatable :: dprojections(:,:,:), temp_realvec_r(:), temp_realvec_i(:)
  real(dp), allocatable, target :: s_dprojections(:,:,:), vnl_dprojections(:,:,:)
  real(dp), allocatable, target :: d2projections(:,:,:)
  real(dp), allocatable :: enlk(:),fnlk(:,:),ddkk(:,:),strnlk(:,:),gmet2(:,:)
  real(dp), allocatable :: work1(:),work2(:),work3(:,:),work4(:,:),work5(:,:,:),work6(:,:,:),work7(:,:,:)
  integer :: idbeg,idend,idfbeg,idfend,dshift,id2beg,id2end,d2shift,dfshift,enlout_shift
  real(dp) :: work(6)
  integer :: ndgxdt_stored,ishift
  integer :: mu0,ic,nu,mu,jc,mua,mub,nua1,nua2,nub1,nub2
  integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
  integer,parameter :: gamma(3,3)=reshape((/1,6,5,6,2,4,5,4,3/),(/3,3/))
  integer          ::  matblk_,natom_,ntypat_,ispden,dimenl2_,ia_beg,ia_end,dimsij,nkpgin_,nkpgout_
  integer, ABI_CONTIGUOUS pointer :: atindx1_(:),indlmn_(:,:,:),nattyp_(:)
  real(dp),pointer :: ffnlin_(:,:,:,:),ffnlout_(:,:,:,:)
  real(dp),pointer :: ph3din_(:,:,:),ph3dout_(:,:,:)
  real(dp),pointer :: enl_(:,:,:,:),enl_ndat_(:,:,:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: sij_(:,:)
  real(dp), ABI_CONTIGUOUS pointer :: kpgin_(:,:),kpgout_(:,:)

! *************************************************************************

  ! We keep the same interface as nonlop, but we don't use many of those
  ABI_UNUSED((/ffnlin,ffnlout,gmet,kpgin,kpgout/))
  ABI_UNUSED((/ph3din,ph3dout/))
  ABI_UNUSED((/ucvol/))
  ABI_UNUSED((/mgfft/))
  ABI_UNUSED((/kptin,kptout/))
  ABI_UNUSED((/idir,nloalg,ngfft,kgin,kgout,ngfft,only_SO,tim_nonlop,gpu_option/))

  ! Check supported options
  if (.not.gemm_nonlop_use_gemm) then
    ABI_BUG('computation not prepared for gemm_nonlop use!')
  end if
  if ( (choice>3.and.choice/=7.and.choice/=5.and.choice/=51.and.signs==2) .or. &
&      (choice>3.and.choice/=7.and.choice/=23.and.choice/=4.and.choice/=54.and.choice/=55.and.choice/=6.and.signs==1) .or. &
&      (useylm/=1) ) then
    ABI_BUG('gemm_nonlop option not supported!')
  end if
  if (signs==1) then
    nnlout_test=0
    if (choice==1) nnlout_test=1
    if (choice==2) nnlout_test=3*natom
    if (choice==3) nnlout_test=6
    if (choice==23) nnlout_test=6+3*natom
    if (nnlout<nnlout_test) then
      ABI_BUG('wrong nnlout size!')
    end if
  end if

  ikin=1; ikout=1;
  select case (select_k)
  case (K_H_K)
    ikin=1; ikout=1;
  case (K_H_KPRIME)
    ikin=2; ikout=1;
  case (KPRIME_H_K)
    ikin=1; ikout=2;
  case (KPRIME_H_KPRIME)
    ikin=2; ikout=2;
  end select
  cplex=2;if (istwf_k>1) cplex=1
  cplex_enl=1;if (paw_opt>0) cplex_enl=2*dimenl1/(lmnmax*(lmnmax+1)) ! is enl complex?
  cplex_fac=max(cplex,dimekbq)
  if ((nspinortot==2.or.cplex_enl==2).and.paw_opt>0.and.choice/=7) cplex_fac=2 ! is vnl_projections complex?
  use_enl_ndat=.false. ; if (size(enl_ndat)>0) use_enl_ndat=.true.

  ! Processing one atom : set pointers to atom-specific arrays (for opernlc)
  if(iatom_only > 0) then
    iatm=atindx1(iatom_only);itypat=typat(iatom_only)
    ia_beg=iatom_only; ia_end=iatom_only
    natom_=1 ; ntypat_=1 ; dimenl2_=1 ; matblk_=1
    ABI_MALLOC(atindx1_,(1))
    ABI_MALLOC(nattyp_,(1))
    atindx1_(1)=1 ; nattyp_(1)=1
    ABI_MALLOC(ph3din_,(2,npwin,1))
    ABI_MALLOC(ph3dout_,(2,npwout,1))
    ph3din_(:,1:npwin,1)=ph3din(:,1:npwin,iatm)
    ph3dout_(:,1:npwout,1)=ph3dout(:,1:npwout,iatm)
    ABI_MALLOC(ffnlin_,(npwin,dimffnlin,lmnmax,1))
    ABI_MALLOC(ffnlout_,(npwout,dimffnlout,lmnmax,1))
    ffnlin_(:,:,:,1)=ffnlin(:,:,:,itypat)
    ffnlout_(:,:,:,1)=ffnlout(:,:,:,itypat)
    ABI_MALLOC(indlmn_,(6,lmnmax,1))
    indlmn_(:,:,1)=indlmn(:,:,itypat)
    if (size(sij)>0) then
      dimsij=size(sij,1)
      ABI_MALLOC(sij_,(dimsij,1))
      if (size(sij,2)==ntypat) then
        sij_(:,1)=sij(:,itypat)
      else if (size(sij)>0) then
        sij_(:,1)=sij(:,1)
      end if
    end if
    if (size(enl_ndat)>0) then
      ABI_MALLOC(enl_ndat_,(size(enl_ndat,1),1,nspinor**2,size(enl_ndat,4),ndat))
      do idat=1,ndat
        do ii=1,size(enl,4)
          do ispden=1,nspinor**2
            if (dimenl2==natom .and. usepaw==1) then
              enl_ndat_(:,1,ispden,idat,ii)=enl_ndat(:,iatom_only,ispden,idat,ii)
            else if (dimenl2==ntypat) then
              enl_ndat_(:,1,ispden,idat,ii)=enl_ndat(:,itypat,ispden,idat,ii)
            else
              enl_ndat_(:,1,ispden,idat,ii)=enl_ndat(:,1,ispden,idat,ii)
            end if
          end do
        end do
      end do
    else if (size(enl)>0) then
      ABI_MALLOC(enl_,(size(enl,1),1,nspinor**2,size(enl,4)))
      do ii=1,size(enl,4)
        do ispden=1,nspinor**2
          if (dimenl2==natom .and. usepaw==1) then
            enl_(:,1,ispden,ii)=enl(:,iatom_only,ispden,ii)
          else if (dimenl2==ntypat) then
            enl_(:,1,ispden,ii)=enl(:,itypat,ispden,ii)
          else
            enl_(:,1,ispden,ii)=enl(:,1,ispden,ii)
          end if
        end do
      end do
    else
      ABI_MALLOC(enl_,(0,0,0,0))
    end if

  ! Usual case : all atoms are processed
  else
    natom_  =natom; ntypat_=ntypat
    ia_beg=1; ia_end=natom
    dimenl2_=dimenl2   ; matblk_=matblk
    atindx1_    => atindx1
    nattyp_     => nattyp
    ffnlin_     => ffnlin
    ffnlout_    => ffnlout
    enl_        => enl
    enl_ndat_   => enl_ndat
    sij_        => sij
    indlmn_     => indlmn
    ph3din_     => ph3din
    ph3dout_    => ph3dout
  end if

  !Eventually re-compute (k+G) vectors (and related data)
  nkpgin_=0
  if (choice==2.or.choice==54) nkpgin_=3
  if (signs==1) then
    if (choice==4) nkpgin_=9
    if (choice==3.or.choice==23.or.choice==6) nkpgin_=3
    if (choice==55) nkpgin_=3
  end if
  if (nkpgin<nkpgin_) then
    ABI_MALLOC(kpgin_,(npwin,nkpgin_))
    call mkkpg(kgin,kpgin_,kptin,nkpgin_,npwin)
  else
    nkpgin_ = nkpgin
    kpgin_  => kpgin
  end if

  nkpgout_=0
  if ((choice==2.or.choice==3.or.choice==54).and.signs==2) nkpgout_=3
  if (nkpgout<nkpgout_) then
    ABI_MALLOC(kpgout_,(npwout,nkpgout_))
    call mkkpg(kgout,kpgout_,kptout,nkpgout_,npwout)
  else
    nkpgout_ = nkpgout
    kpgout_ => kpgout
  end if

  ! The number of projectors used for computation may vary among
  ! nonlop calls, from computing on all atoms to a select one for
  ! some perturbations.
  ! In such cases, projs arrays must be recomputed
  nprojs=0
  do itypat=1,ntypat_
    nprojs = nprojs + count(indlmn_(3,:,itypat)>0)*nattyp_(itypat)
  end do

  ! If vectproj is provided, use it for further calculations, use allocated array otherwise
  local_vectproj=.false.
  if(PRESENT(vectproj)) then
    if(size(vectproj)>1) local_vectproj=.true.
  end if
  if (local_vectproj) projections => vectproj


  if(nprojs == 0) then
    ! TODO check if this is correct
    if(signs == 1) then
      enlout=zero
      return
    end if
    if(signs == 2) then
      vectout = zero
      if(paw_opt>0) svectout = vectin
      return
    end if
  end if

  if(signs == 1) then
    enlout=zero
    ABI_MALLOC(enlk,(ndat))
    enlk=zero
    ABI_MALLOC(fnlk,(3*natom,ndat))
    fnlk=zero
    ABI_MALLOC(ddkk,(6,ndat))
    ddkk=zero
    ABI_MALLOC(strnlk,(6,ndat))
    strnlk=zero
  end if

  ndgxdt = -1
  nd2gxdt = -1

  ndgxdtfac = 0; nd2gxdtfac = 0
  if (choice==2) then
    if (signs==1) ndgxdt=3
    if (signs==2) ndgxdt=1
    if (signs==2) ndgxdtfac=1
  end if
  if (choice==22) then
    if (signs==2) ndgxdtfac=1
  end if
  if (choice==23) then
    if (signs==1) ndgxdt=9
  end if
  if (choice==3) then
    if (signs==1) ndgxdt=6
    if (signs==2) ndgxdt=1
    if (signs==2) ndgxdtfac=1
  end if
  if (choice==4) then
    if(signs==1) ndgxdt=3
    if(signs==1) ndgxdtfac=3
    if(signs==1) nd2gxdt=6
  end if
  if (choice==5) then
    if(signs==1) ndgxdt=3
    if(signs==2) ndgxdt=1
    if(signs==2) ndgxdtfac=1
  end if
  if (choice==51) then
    if(signs==1) ndgxdt=3
    if(signs==2) ndgxdt=1
    if(signs==2) ndgxdtfac=1
  end if
  if (choice==54) then
    if(signs==1) ndgxdt=6
    if(signs==1) ndgxdtfac=6
    if(signs==1) nd2gxdt=9
    if(signs==2) ndgxdt=1
    if(signs==2) nd2gxdt=1
    if(signs==2) ndgxdtfac=1
    if(signs==2) nd2gxdtfac=1
  end if
  if (choice==55) then
    if(signs==1) ndgxdt=9
    if(signs==1) ndgxdtfac=9
    if(signs==1) nd2gxdt=18
  end if
  if (choice==6) then
    if(signs==1) ndgxdt=9
    if(signs==1) ndgxdtfac=9
    if(signs==1) nd2gxdt=54
  end if
  ngrads=0; ngrads2=0
  if(ndgxdt>0) ngrads=ndgxdt; if(ndgxdt>0) ngrads2=nd2gxdt
  if(ndgxdt>0) then
    ABI_CHECK(ndgxdtfac<=ndgxdt,"BUG: ndgxdtfac>ndgxdt!")
  end if
  optder = 0;if (ndgxdtfac>0) optder = 1
  if (nd2gxdtfac>0) optder=2
  cplex_dgxdt(:) = 1 ; cplex_d2gxdt(:) = 1
  ! When istwf_k > 1, gx derivatives can be real or pure imaginary
  ! cplex_dgxdt(i)  = 1 if dgxdt(1,i,:,:)  is real, 2 if it is pure imaginary
  ! cplex_d2gxdt(i) = 1 if d2gxdt(1,i,:,:) is real, 2 if it is pure imaginary
  if(ndgxdt > 0) then
   if (choice==5.or.choice==51) cplex_dgxdt(:) = 2
   if (choice==54.and.signs==1) cplex_dgxdt(4:6) = 2
   !if (choice==54.and.signs==2) cplex_dgxdt(:)   = 2
   if (choice==55.and.signs==1) cplex_dgxdt(7:9) = 2
  end if
  if(nd2gxdt > 0) then
    if (choice==54) cplex_d2gxdt(:) = 2
    if (choice==55.and.signs==1) cplex_d2gxdt(1:18)= 2
  end if

  ! These will store the non-local factors for vectin, svectout and vectout respectively
  if(.not. local_vectproj) then
    ABI_MALLOC(projections,(cplex, nprojs,nspinor*ndat))
  end if
  ABI_MALLOC(s_projections,(cplex, nprojs,nspinor*ndat))
  ABI_MALLOC(vnl_projections,(cplex_fac, nprojs,nspinor*ndat))

  if(.not. local_vectproj) projections = zero
  s_projections = zero
  vnl_projections = zero

  ! Working buffers for storing derivative
  if (ndgxdt>0) then
    ABI_MALLOC(dprojections,(cplex, ndgxdt*nprojs,nspinor*ndat))
    dprojections(:,:,:) = zero
  else
    ABI_MALLOC(dprojections,(1,1,ndat))
  end if

  if (ndgxdtfac>0) then
    ABI_MALLOC(s_dprojections,(cplex, ndgxdtfac*nprojs,nspinor*ndat))
    ABI_MALLOC(vnl_dprojections,(cplex_fac, ndgxdtfac*nprojs,nspinor*ndat))
    s_dprojections(:,:,:) = zero
    vnl_dprojections(:,:,:) = zero
  else
    ABI_MALLOC(s_dprojections,(1,1,ndat))
    ABI_MALLOC(vnl_dprojections,(1,1,ndat))
  end if

  ! Working buffers for storing 2nd-derivative
  if (nd2gxdt>0) then
    ABI_MALLOC(d2projections,(cplex, nd2gxdt*nprojs, nspinor*ndat))
    if(cpopt < 4) then
      d2projections(:,:,:) = zero
    end if
  else
    ABI_MALLOC(d2projections,(1, 1, ndat))
  end if


  ! determine precisely when temp_realvec_r~i needs to be allocated
  ! to factorize allocate (resp. deallocate) at the begining (resp. at the end) of subroutine
  ! to avoid multiple allocate/deallocate that can be costly
  if (cplex /= 2) then
    if ( (cpopt < 2) .or. &
      &  (paw_opt == 3 .or. paw_opt == 4) .or. &
      &  (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4)) then
       ABI_MALLOC(temp_realvec_r,(MAX(npwout,npwin)*nspinor*ndat))
       ABI_MALLOC(temp_realvec_i,(MAX(npwout,npwin)*nspinor*ndat))
    end if
  end if

  if(cpopt >= 2) then
    ! retrieve from cprjin
    if(.not. local_vectproj .and. cpopt/=3) then
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = ia_beg, ia_end
          nlmn = cprjin(iatom, idat)%nlmn
          projections(1:cplex, shift+1:shift+nlmn, idat) = cprjin(iatom, idat)%cp(1:cplex, 1:nlmn)
          shift = shift + nlmn
        end do
      end do
    end if
    if(cpopt==4.and.allocated(dprojections)) then
      ABI_CHECK(cprjin(1,1)%ncpgr>=ndgxdt,"cprjin%ncpgr not correct! (1)")
      ndgxdt_stored = cprjin(1,1)%ncpgr
      ishift=0
      if (((choice==2).or.(choice==3)).and.(ndgxdt_stored>ndgxdt).and.(signs==2)) ishift=idir-ndgxdt
      if ((choice==2).and.(ndgxdt_stored==9).and.(signs==2)) ishift=ishift+6
      if (choice==2.and.(ndgxdt_stored>ndgxdt).and.(signs==1)) ishift=ndgxdt_stored-ndgxdt
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = ia_beg, ia_end
          nlmn  = cprjin(iatom, idat)%nlmn
          do ilmn=1,nlmn
            do igrad=1,ndgxdt
              dprojections(1:cplex, shift + igrad, idat) = &
                cprjin(iatom, idat)%dcp(1:cplex,igrad+ishift,ilmn)
            end do
            shift = shift + ndgxdt
          end do
        end do
      end do
    end if
  end if ! cpopt

  if(cpopt<=1.or.(cpopt<=3.and.(choice==2.or.choice==3.or.choice==5.or.choice==51.or.choice==23.or.choice==54.or.choice==55.or.choice==4))) then

    call opernla_gemm(choice,cplex,cplex_dgxdt,cplex_d2gxdt,dimffnlin,&
    &       d2projections,dprojections,ffnlin,projections,&
    &       idir,indlmn,istwf_k,kpgin_,matblk,mpi_enreg,nd2gxdt,ndgxdt,nkpgin_,&
    &       npwin,nspinor,ph3din,signs,ucvol,ndat,ntypat,lmnmax,nattyp,(ikin==2),&
    &       iatom_only,atom_proj_shift,cpopt,&
    &       nprojs,&
    &       vectin,&
    &       temp_realvec_r,&
    &       gpu_option,gemm_nonlop_is_distributed)

    if(cpopt >= 0) then
      ! store in cprjin
      if(.not. local_vectproj .and. cpopt/=3) then
        do idat=1, ndat*nspinor
          shift = 0
          do iatom = ia_beg, ia_end
            nlmn = cprjin(iatom, idat)%nlmn
            cprjin(iatom, idat)%cp(1:cplex, 1:nlmn) = projections(1:cplex, shift+1:shift+nlmn, idat)
            shift = shift + nlmn
          end do
        end do
      end if
      if(cpopt==1 .or. cpopt==3) then
        ABI_CHECK(cprjin(1,1)%ncpgr>=ndgxdt,"cprjin%ncpgr not correct! (2)")
        do idat=1, ndat*nspinor
          shift = 0
          do iatom = ia_beg, ia_end
            nlmn = cprjin(iatom, idat)%nlmn
            do ilmn=1,nlmn
              do igrad=1,ndgxdt
                cprjin(iatom, idat)%dcp(1:cplex,igrad,ilmn) = &
                &                   dprojections(1:cplex, shift + igrad, idat)
              end do
              shift = shift + ndgxdt
            end do
          end do
        end do
      end if
    end if ! cpopt >= 0
  end if ! cpopt >= 2

  if(choice > 0) then

    if(choice /= 7) then
      ! opernlc
      iatm = 0
      shift = 0; dshift = 0; dfshift = 0; d2shift = 0
      ABI_MALLOC(sij_typ,(((paw_opt+1)/3)*lmnmax*(lmnmax+1)/2))
      do itypat=1, ntypat_
        nlmn=count(indlmn_(3,:,itypat)>0)
        if (paw_opt>=2) then
          if (cplex_enl==1) then
            do ilmn=1,nlmn*(nlmn+1)/2
              sij_typ(ilmn)=sij_(ilmn,itypat)
            end do
          else
            do ilmn=1,nlmn*(nlmn+1)/2
              sij_typ(ilmn)=sij_(2*ilmn-1,itypat)
            end do
          end if
        end if

        ibeg = shift+1
        iend = shift+nattyp_(itypat)*nlmn

        idbeg = dshift+1
        idend = dshift+nattyp_(itypat)*nlmn*ngrads

        idfbeg = dshift+1
        idfend = dshift+nattyp_(itypat)*nlmn*ndgxdtfac

        id2beg = d2shift+1
        id2end = d2shift+nattyp_(itypat)*nlmn*ngrads2

        do idat = 1,ndat
        if(use_enl_ndat) enl_ => enl_ndat_(:,:,:,idat,:)
          call opernlc_ylm(atindx1_,cplex,cplex_dgxdt,cplex_d2gxdt,&
&         cplex_enl,cplex_fac,&
&         dprojections(:, idbeg:idend, 1+nspinor*(idat-1):nspinor*idat),&
&         vnl_dprojections(:, idfbeg:idfend, 1+nspinor*(idat-1):nspinor*idat),&
&         s_dprojections(:, idfbeg:idfend, 1+nspinor*(idat-1):nspinor*idat),&
&         d2projections(:, id2beg:id2end, 1+nspinor*(idat-1):nspinor*idat),&
&         d2gxdt_dum_out,d2gxdt_dum_out2,&
&         dimenl1,dimenl2_,dimekbq,enl_,&
&         projections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
&         vnl_projections(:, ibeg:iend,1+nspinor*(idat-1):nspinor*idat),&
&         s_projections(:, ibeg:iend,1+nspinor*(idat-1):nspinor*idat),&
&         iatm,indlmn_(:,:,itypat),itypat,lambda(idat),mpi_enreg,natom_,&
&         ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,&
&         nattyp_(itypat),nlmn,nspinor,nspinortot,optder,paw_opt,sij_typ)
        end do

        shift = shift + nattyp_(itypat)*nlmn
        dshift = dshift + nattyp_(itypat)*nlmn*ngrads
        dfshift = dshift + nattyp_(itypat)*nlmn*ndgxdtfac
        d2shift = d2shift + nattyp_(itypat)*nlmn*ngrads2
        iatm = iatm+nattyp_(itypat)
      end do
      ABI_FREE(sij_typ)
    else
      s_projections = projections
    end if ! choice /= 7

    ! opernlb
    if(signs==2) then

      call opernlb_gemm(choice,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_fac,&
      &       d2gxdt_dum_in,d2gxdt_dum_out,&
      &       vnl_dprojections,s_dprojections,dimffnlout,ffnlout,&
      &       vnl_projections,s_projections,&
      &       idir,indlmn,kpgout_,matblk,istwf_k,&
      &       nd2gxdt,nd2gxdtfac,ndgxdt,ndgxdtfac,&
      &       nkpgout_,npwout,nspinor,signs,ucvol,ndat,&
      &       ntypat,lmnmax,nattyp,(ikout==2),iatom_only,atom_proj_shift,&
      &       paw_opt,ph3dout,&
      &       nprojs,&
      &       vectin,vectout,svectout,&
      &       temp_realvec_r,temp_realvec_i,&
      &       gpu_option,gemm_nonlop_is_distributed)
    end if

    ! opernld
    if(signs==1) then
      if(choice==2 .or. choice==3 .or. choice==23 .or. choice==4 .or. choice==54 .or. choice==55 .or. choice==6) then
        call opernld_ylm_allwf(choice,cplex,cplex_fac,ddkk,&
        &       dprojections,vnl_dprojections,s_dprojections,d2projections,&
        &       enlk,enlout,fnlk,projections,vnl_projections,s_projections,&
        &       natom,ndat,nd2gxdt,ndgxdt,&
        &       ndgxdtfac,indlmn_,ntypat_,lmnmax,nprojs,nnlout,nspinor,paw_opt,&
        &       strnlk,nattyp_,gpu_option)
      else
        shift=0; dshift=0; dfshift = 0; d2shift = 0; iatm=1
        do itypat=1, ntypat_
          nlmn=count(indlmn_(3,:,itypat)>0)

          ibeg = shift+1
          iend = shift+nattyp_(itypat)*nlmn

          idbeg = dshift+1
          idend = dshift+nattyp_(itypat)*nlmn*ngrads

          idfbeg = dshift+1
          idfend = dshift+nattyp_(itypat)*nlmn*ndgxdtfac

          id2beg = d2shift+1
          id2end = d2shift+nattyp_(itypat)*nlmn*ngrads2

          do idat=1,ndat
            call opernld_ylm             (choice,cplex,cplex_fac,ddkk(:,idat),&
            &       dprojections    (:, idbeg:idend, 1+nspinor*(idat-1):nspinor*idat),&
            &       vnl_dprojections(:, idfbeg:idfend, 1+nspinor*(idat-1):nspinor*idat),&
            &       s_dprojections  (:, idfbeg:idfend, 1+nspinor*(idat-1):nspinor*idat),&
            &       d2projections (:, id2beg:id2end, 1+nspinor*(idat-1):nspinor*idat),&
            &       enlk(idat),enlout(nnlout*(idat-1)+1:nnlout*idat),fnlk(:,idat),&
            &       projections    (:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
            &       vnl_projections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
            &       s_projections  (:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
            &       iatm,natom_,1,nd2gxdt,ndgxdt,ndgxdtfac,&
            &       nattyp_(itypat),nlmn,nnlout,nspinor,paw_opt,strnlk(:,idat))
          end do

          shift = shift + nattyp_(itypat)*nlmn
          dshift = dshift + nattyp_(itypat)*nlmn*ngrads
          dfshift = dshift + nattyp_(itypat)*nlmn*ndgxdtfac
          d2shift = d2shift + nattyp_(itypat)*nlmn*ngrads2
          iatm = iatm+nattyp_(itypat)
        end do
      end if

      ! Reduction in case of parallelism
      if (mpi_enreg%paral_spinor==1) then
        if (size(enlout)>0) then
          call xmpi_sum(enlout,mpi_enreg%comm_spinor,ierr)
        end if
        if (choice==3.or.choice==23) then
          call xmpi_sum(enlk,mpi_enreg%comm_spinor,ierr)
        end if
        if (choice==55) then
          call xmpi_sum(ddkk,mpi_enreg%comm_spinor,ierr)
        end if
      end if

      !Need sometimes gmet
      if ((signs==1.and.paw_opt<=3).and. &
          & (choice==5 .or.choice==51.or.choice==52.or.choice==53.or.&
          & choice==54.or.choice==55)) then
        ABI_MALLOC(gmet2,(3,3))
        gmet2 = MATMUL(TRANSPOSE(gprimd),gprimd)
      end if

      !Coordinate transformations

      ! Derivatives wrt strain
      !  - Convert from reduced to cartesian coordinates
      !  - Substract volume contribution
      if ((choice==3.or.choice==23).and.paw_opt<=3) then
        do idat=1,ndat
          enlout_shift=(idat-1)*nnlout
          call strconv(enlout(enlout_shift+1:enlout_shift+6),gprimd,work)
          enlout(enlout_shift+1:enlout_shift+3)=(work(1:3)-enlk(idat))
          enlout(enlout_shift+4:enlout_shift+6)= work(4:6)
        end do
      end if

      !2nd derivative wrt to k wave vector and atomic position (effective charges):
      ! - convert from cartesian to reduced coordinates
      if (choice==54.and.signs==1.and.paw_opt<=3) then
        ABI_MALLOC(work1,(3))
        ABI_MALLOC(work2,(3))
        do idat=1,ndat
          mu0=0 ! Shift to be applied in enlout array
          enlout_shift=(idat-1)*nnlout
          do mu=1,3*natom
        !   First, real part
            work1(1)=enlout(enlout_shift+mu0+1);work1(2)=enlout(enlout_shift+mu0+3);work1(3)=enlout(enlout_shift+mu0+5)
            work2(:)=gmet2(:,1)*work1(1)+gmet2(:,2)*work1(2)+gmet2(:,3)*work1(3)
            enlout(enlout_shift+mu0+1)=work2(1);enlout(enlout_shift+mu0+3)=work2(2);enlout(enlout_shift+mu0+5)=work2(3)
        !   Then imaginary part
            work1(1)=enlout(enlout_shift+mu0+2);work1(2)=enlout(enlout_shift+mu0+4);work1(3)=enlout(enlout_shift+mu0+6)
            work2(:)=gmet2(:,1)*work1(1)+gmet2(:,2)*work1(2)+gmet2(:,3)*work1(3)
            enlout(enlout_shift+mu0+2)=work2(1);enlout(enlout_shift+mu0+4)=work2(2);enlout(enlout_shift+mu0+6)=work2(3)
            mu0=mu0+6
          end do
        end do !idat
        ABI_FREE(work1)
        ABI_FREE(work2)
      end if

      !2nd derivative wrt to k wave vector and strain (piezoelectric tensor):
      ! - convert from cartesian to reduced coordinates (k point)
      ! - convert from reduced to cartesian coordinates (strain)
      ! - substract volume contribution
      ! - symetrize strain components
      if (choice==55.and.signs==1.and.paw_opt<=3) then
        ABI_MALLOC(work3,(2,3))
        ABI_MALLOC(work4,(2,3))
        ABI_MALLOC(work5,(2,3,6))
        ABI_MALLOC(work7,(2,3,6))
        ABI_MALLOC(work6,(2,3,3))
        do idat=1,ndat
          enlout_shift=(idat-1)*nnlout
          do ic=1,3 ! gamma
            work5=zero
            do jc=1,3 ! nu
              do ii=1,3 ! lambda
                mu=(gamma(jc,ii)-1)*3+1
                work5(1,jc,ii)=gmet2(ic,1)*enlout(enlout_shift+2*mu-1)+gmet2(ic,2)*enlout(enlout_shift+2*mu+1) &
       &         +gmet2(ic,3)*enlout(enlout_shift+2*mu+3)
                work5(2,jc,ii)=gmet2(ic,1)*enlout(enlout_shift+2*mu  )+gmet2(ic,2)*enlout(enlout_shift+2*mu+2) &
       &         +gmet2(ic,3)*enlout(enlout_shift+2*mu+4)
              end do
            end do
            work6=zero
            do jc=1,3 ! nu
              do ii=1,3 ! beta
                work6(1:cplex,ii,jc)=gprimd(ii,1)*work5(1:cplex,jc,1)+gprimd(ii,2)*work5(1:cplex,jc,2) &
       &         +gprimd(ii,3)*work5(1:cplex,jc,3)
              end do
            end do
            do jc=1,3 ! alpha
              do ii=1,3 ! beta
                mu=gamma(jc,ii)
                work7(1:cplex,ic,mu)=gprimd(jc,1)*work6(1:cplex,ii,1)+gprimd(jc,2)*work6(1:cplex,ii,2) &
       &         +gprimd(jc,3)*work6(1:cplex,ii,3)
              end do
            end do
          end do ! gamma

          do ii=1,3 ! alpha
            work3(1,ii)=gprimd(ii,1)*ddkk(2*1-1,idat)+gprimd(ii,2)*ddkk(2*2-1,idat) &
       &     +gprimd(ii,3)*ddkk(2*3-1,idat)
            work3(2,ii)=gprimd(ii,1)*ddkk(2*1  ,idat)+gprimd(ii,2)*ddkk(2*2  ,idat) &
       &     +gprimd(ii,3)*ddkk(2*3  ,idat)
          end do
          do ii=1,3 ! gamma
            work4(1,ii)=gmet2(ii,1)*ddkk(2*1-1,idat)+gmet2(ii,2)*ddkk(2*2-1,idat) &
       &     +gmet2(ii,3)*ddkk(2*3-1,idat)
            work4(2,ii)=gmet2(ii,1)*ddkk(2*1  ,idat)+gmet2(ii,2)*ddkk(2*2  ,idat) &
       &     +gmet2(ii,3)*ddkk(2*3  ,idat)
          end do

          do mu=1,6
            ii=alpha(mu) ! alpha
            ic=beta(mu) ! beta
            do jc=1,3 ! gamma
              work7(1:cplex,jc,mu)=work7(1:cplex,jc,mu)-half &
       &       *(gprimd(ic,jc)*work3(1:cplex,ii)+gprimd(ii,jc)*work3(1:cplex,ic))
              if (ii==ic) work7(1:cplex,jc,mu)=work7(1:cplex,jc,mu)-work4(1:cplex,jc)
            end do
          end do
          do mu=1,6 ! alpha,beta
            do nu=1,3 ! gamma
              mu0=3*(mu-1)+nu
              enlout(enlout_shift+2*mu0-1)=work7(1,nu,mu)
              enlout(enlout_shift+2*mu0  )=work7(2,nu,mu)
            end do
          end do
        end do !idat
        ABI_FREE(work3)
        ABI_FREE(work4)
        ABI_FREE(work5)
        ABI_FREE(work6)
        ABI_FREE(work7)
      end if


     !2nd derivative wrt to 2 strains (elastic tensor):
     ! - convert from reduced to cartesian coordinates
     ! - substract volume contribution
      if (choice==6.and.signs==1.and.paw_opt<=3) then
        ABI_MALLOC(work1,(6))
        ABI_MALLOC(work2,(6))
        ABI_MALLOC(work3,(6+3*natom,6))
        do idat=1,ndat
          mu0=(idat-1)*nnlout ! Shift to be applied in enlout array
          work3(:,:)=reshape(enlout(mu0+1:mu0+6*(6+3*natom)),(/6+3*natom,6/))
          do mu=1,6
            call strconv(work3(1:6,mu),gprimd,work3(1:6,mu))
          end do
          do mu=1,6+3*natom
            work1(1:6)=work3(mu,1:6)
            call strconv(work1,gprimd,work2)
            work3(mu,1:6)=work2(1:6)
          end do
          enlout(mu0+1:mu0+6*(6+3*natom))=reshape(work3(:,:),(/6*(6+3*natom)/))
          call strconv(strnlk(:,idat),gprimd,strnlk(:,idat))
          do mub=1,6
            nub1=alpha(mub);nub2=beta(mub)
            do mua=1,6
              mu=mu0+mua+(3*natom+6)*(mub-1)
              nua1=alpha(mua);nua2=beta(mua)
              if (mua<=3.and.mub<=3) enlout(mu)=enlout(mu)+enlk(idat)
              if (mua<=3) enlout(mu)=enlout(mu)-strnlk(mub,idat)
              if (mub<=3) enlout(mu)=enlout(mu)-strnlk(mua,idat)
              if (nub1==nua2) enlout(mu)=enlout(mu)-0.25d0*strnlk(gamma(nua1,nub2),idat)
              if (nub2==nua2) enlout(mu)=enlout(mu)-0.25d0*strnlk(gamma(nua1,nub1),idat)
              if (nub1==nua1) enlout(mu)=enlout(mu)-0.25d0*strnlk(gamma(nua2,nub2),idat)
              if (nub2==nua1) enlout(mu)=enlout(mu)-0.25d0*strnlk(gamma(nua2,nub1),idat)
            end do
            if (mub<=3) then
              do nua1=1,natom
                nua2=3*(nua1-1);mu=mu0+nua2+6+(3*natom+6)*(mub-1)
                enlout(mu+1:mu+3)=enlout(mu+1:mu+3)-fnlk(nua2+1:nua2+3,idat)
              end do
            end if
          end do
        end do
        ABI_FREE(work1)
        ABI_FREE(work2)
        ABI_FREE(work3)
      end if

    end if !opernld

  end if ! choice>0

! Release memory

  if (iatom_only>0) then
    ABI_FREE(atindx1_)
    ABI_FREE(nattyp_)
    ABI_FREE(ph3din_)
    ABI_FREE(ph3dout_)
    ABI_FREE(ffnlin_)
    ABI_FREE(ffnlout_)
    if(use_enl_ndat) then
      ABI_FREE(enl_ndat_)
    else
      ABI_FREE(enl_)
    end if
    ABI_FREE(indlmn_)
    if (size(sij) > 1) then
      ABI_FREE(sij_)
    end if
  end if

  if (allocated(enlk)) then
    ABI_FREE(enlk)
    ABI_FREE(fnlk)
    ABI_FREE(strnlk)
    ABI_FREE(ddkk)
  end if

  if (nkpgin<nkpgin_) then
    ABI_FREE(kpgin_)
  end if
  if (nkpgout<nkpgout_) then
    ABI_FREE(kpgout_)
  end if

  if (allocated(gmet2)) then
    ABI_FREE(gmet2)
  end if

  if(.not. local_vectproj) then
    ABI_FREE(projections)
  end if
  ABI_FREE(s_projections)
  ABI_FREE(vnl_projections)
  if (allocated(dprojections)) then
    ABI_FREE(dprojections)
  end if
  if (allocated(s_dprojections)) then
    ABI_FREE(s_dprojections)
  end if
  if (allocated(vnl_dprojections)) then
    ABI_FREE(vnl_dprojections)
  end if
  if (allocated(d2projections)) then
    ABI_FREE(d2projections)
  end if
  if (allocated(temp_realvec_r)) then
    ABI_FREE(temp_realvec_r)
    ABI_FREE(temp_realvec_i)
  end if

 end subroutine gemm_nonlop
!***

end module m_gemm_nonlop
!!***
