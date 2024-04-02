!!****m* ABINIT/m_gemm_nonlop_ompgpu
!! NAME
!! m_gemm_nonlop_ompgpu
!!
!! FUNCTION
!!  This module provides functions to compute the nonlocal operator by means of the BLAS GEMM
!!  routine. By treating ndat simultaneous wavefunctions, it is able to exploit BLAS3 routines,
!!  which leads to excellent CPU efficiency and OpenMP scalability.
!!
!! COPYRIGHT
!! Copyright (C) 2014-2024 ABINIT group (MS)
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

module m_gemm_nonlop_ompgpu

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_xomp
 use m_abi_linalg

#ifdef HAVE_FC_ISO_C_BINDING
 use iso_c_binding
#endif

#ifdef HAVE_GPU
  use m_gpu_toolbox
#endif

 use defs_abitypes, only : MPI_type
 use m_opernlc_ylm_ompgpu, only : opernlc_ylm_ompgpu
 use m_opernlc_ylm, only : opernlc_ylm
 use m_opernla_gemm, only : opernla_gemm
 use m_opernlb_gemm, only : opernlb_gemm
 use m_opernld_ylm_allwf_ompgpu, only : opernld_ylm_allwf_ompgpu
 use m_opernld_ylm, only : opernld_ylm
 use m_pawcprj, only : pawcprj_type
 use m_geometry, only : strconv
 use m_kg, only : mkkpg
 use m_gemm_nonlop

#if defined HAVE_MPI2
 use mpi
#endif

 implicit none

 private

 ! then call gemm_nonlop to do the actual computation, and call destroy when done. See gstate and vtorho.
 public :: gemm_nonlop_ompgpu

 ! Those routines are here to assess memory requirements
 public :: gemm_nonlop_ompgpu_work_mem
 public :: gemm_nonlop_ompgpu_static_mem
!!***

!----------------------------------------------------------------------

#ifdef HAVE_OPENMP_OFFLOAD
 integer, save :: gpu_initialised=0
 integer, save :: mod__ndat=0
 integer, save :: mod__nprojs=0
 real(dp), save, allocatable, target :: projections_(:,:,:), s_projections(:,:,:), vnl_projections(:,:,:), dprojections(:,:,:)
 real(dp), save, allocatable :: temp_realvec_r(:),temp_realvec_i(:)
 real(dp), save, allocatable :: sij_typ(:,:)

#endif

contains

 function gemm_nonlop_ompgpu_work_mem(istwfk, ndat, npw, indlmn, nattyp, ntypat, lmnmax) result(req_mem)
   implicit none

   integer, intent(in) :: istwfk, ndat, npw, ntypat, lmnmax
   integer, intent(in) :: indlmn(:,:,:), nattyp(ntypat)

   integer :: nprojs, cplex, itypat
   real(dp) :: req_mem

! *************************************************************************

   cplex=2;if (istwfk>1) cplex=1
   nprojs=0
   do itypat=1,ntypat
     nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
   end do

   req_mem = 0

   if(cplex == 1) then
     req_mem = req_mem + dp * npw * ndat ! temp_realvec_r
     req_mem = req_mem + dp * npw * ndat ! temp_realvec_i
   end if

   req_mem = req_mem + dp * lmnmax * (lmnmax+1)/2 * ntypat  ! sij_typ

   req_mem = req_mem + dp * cplex * int(nprojs, c_size_t) * int(ndat, c_size_t)  ! projections_
   req_mem = req_mem + dp * cplex * int(nprojs, c_size_t) * int(ndat, c_size_t)  ! s_projections
   req_mem = req_mem + dp * cplex * int(nprojs, c_size_t) * int(ndat, c_size_t)  ! vnl_projections

 end function gemm_nonlop_ompgpu_work_mem

!----------------------------------------------------------------------

 function gemm_nonlop_ompgpu_static_mem(npw, indlmn, nattyp, ntypat, nblocks,&
     compute_grad_strain,compute_grad_atom) result(req_mem)
   implicit none

   integer, intent(in) :: npw, ntypat, nblocks
   integer, intent(in) :: indlmn(:,:,:), nattyp(ntypat)
   logical, intent(in), optional :: compute_grad_strain,compute_grad_atom

   integer :: nprojs, ngrads, itypat
   logical :: my_compute_grad_strain,my_compute_grad_atom
   integer(kind=c_size_t) :: req_mem

! *************************************************************************

   my_compute_grad_strain=.false. ; if (present(compute_grad_strain)) my_compute_grad_strain=compute_grad_strain
   my_compute_grad_atom=.false. ; if (present(compute_grad_atom)) my_compute_grad_atom=compute_grad_atom

   nprojs = 0
   do itypat=1,ntypat
     nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
   end do
   nprojs = nprojs / nblocks

   ngrads = 0
   if (my_compute_grad_strain) ngrads = ngrads + 6
   if (my_compute_grad_atom)   ngrads = ngrads + 3

   req_mem = 0

   if(nblocks>1) then
     req_mem = req_mem + dp * 2 * int(npw, c_size_t) * int(nprojs, c_size_t)          !projs_recv
   end if
   ! projs or projs_r + projs_i
   req_mem = req_mem + 2 * dp * int(npw, c_size_t) * int(nprojs, c_size_t)
   if(ngrads>0) then
     ! dprojs or dprojs_r + dprojs_i
     req_mem = req_mem + 2 * dp * int(npw, c_size_t) * int(ngrads, c_size_t) * int(nprojs, c_size_t)
     if(nblocks>1) then
       req_mem = req_mem + dp * 2 * int(npw, c_size_t) * int(ngrads, c_size_t)*int(nprojs, c_size_t)   !dprojs_recv
     end if
   end if

 end function gemm_nonlop_ompgpu_static_mem

!----------------------------------------------------------------------

!Tested usecases :
! - Nvidia GPUs : FC_NVHPC + CUDA
! - AMD GPUs    : FC_LLVM + HIP
! An eventual Intel implementation would use the OneAPI LLVM compiler.
! Homemade CUDA/HIP interfaces would allow the use of GCC.
! But it is likely that OpenMP performance won't be optimal outside GPU vendors compilers.
#ifdef HAVE_OPENMP_OFFLOAD

 subroutine refresh_gemm_nonlop_kpt_ompgpu(ikpt)
  integer,intent(in) :: ikpt

  character(len=500) :: msg

! *************************************************************************

  if(current_ikpt_in_gpu == ikpt) then
    msg="GPU nonlop arrays are already initialized for K-Point. Redundant call"
    ABI_ERROR(msg)
  end if
  if(gemm_nonlop_kpt(ikpt)%nprojs == -1) then
    msg="Was asked to work on a K-point index for which arrays aren't built.\n" // &
&    "Requested K-point index was %d.\n" // &
&    "Please make sure that make_gemm_nonlop is called for this K-point before reaching this routine."
    ABI_ERROR(msg)
  end if

  call free_gemm_nonlop_kpt_ompgpu

  gpu_nonlop_current_ikpt => gemm_nonlop_kpt(ikpt)
  if(allocated(gpu_nonlop_current_ikpt%projs)) then
    current_ikpt_projs   => gpu_nonlop_current_ikpt%projs
    !$OMP TARGET ENTER DATA MAP(to:current_ikpt_projs)
    nullify(current_ikpt_projs_r)
    nullify(current_ikpt_projs_i)
  else
    current_ikpt_projs_r => gpu_nonlop_current_ikpt%projs_r
    current_ikpt_projs_i => gpu_nonlop_current_ikpt%projs_i
    !$OMP TARGET ENTER DATA MAP(to:current_ikpt_projs_r)
    !$OMP TARGET ENTER DATA MAP(to:current_ikpt_projs_i)
    nullify(current_ikpt_projs)
  end if
  if(gpu_nonlop_current_ikpt%ngrads /= -1) then
    if(allocated(gpu_nonlop_current_ikpt%dprojs)) then
      current_ikpt_dprojs   => gpu_nonlop_current_ikpt%dprojs
      !$OMP TARGET ENTER DATA MAP(to:current_ikpt_dprojs)
      nullify(current_ikpt_dprojs_r)
      nullify(current_ikpt_dprojs_i)
    else
      current_ikpt_dprojs_r => gpu_nonlop_current_ikpt%dprojs_r
      current_ikpt_dprojs_i => gpu_nonlop_current_ikpt%dprojs_i
      !$OMP TARGET ENTER DATA MAP(to:current_ikpt_dprojs_i)
      !$OMP TARGET ENTER DATA MAP(to:current_ikpt_dprojs_r)
      nullify(current_ikpt_dprojs)
    end if
  end if

  current_ikpt_in_gpu=ikpt

 end subroutine refresh_gemm_nonlop_kpt_ompgpu

!----------------------------------------------------------------------

 subroutine free_gemm_nonlop_kpt_ompgpu()

  if(current_ikpt_in_gpu == -1) return

  if(associated(current_ikpt_projs)) then
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_projs)
    nullify(current_ikpt_projs)
  end if
  if(associated(current_ikpt_projs_r)) then
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_projs_i)
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_projs_r)
    nullify(current_ikpt_projs_r)
    nullify(current_ikpt_projs_i)
  end if
  if(associated(current_ikpt_dprojs)) then
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_dprojs)
    nullify(current_ikpt_dprojs)
  end if
  if(associated(current_ikpt_dprojs_i)) then
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_dprojs_i)
    !$OMP TARGET EXIT DATA MAP(delete:current_ikpt_dprojs_r)
    nullify(current_ikpt_dprojs_r)
    nullify(current_ikpt_dprojs_i)
  end if

  current_ikpt_in_gpu=-1

 end subroutine free_gemm_nonlop_kpt_ompgpu

!----------------------------------------------------------------------

 subroutine alloc_work_buffers(cplex, cplex_fac, ndat, nprojs, ntypat, lmnmax, npw)

  integer,intent(in) :: cplex, cplex_fac, ndat, nprojs, ntypat, lmnmax, npw

! *************************************************************************

  call free_work_buffers()

  mod__ndat=ndat
  mod__nprojs=nprojs

  if(cplex == 1) then
    ABI_MALLOC(temp_realvec_r,(npw*ndat))
    ABI_MALLOC(temp_realvec_i,(npw*ndat))
    !$OMP TARGET ENTER DATA MAP(alloc:temp_realvec_r,temp_realvec_i)
  end if

  ABI_MALLOC(sij_typ,(lmnmax*(lmnmax+1)/2,ntypat))
  !$OMP TARGET ENTER DATA MAP(alloc:sij_typ)

  ABI_MALLOC(projections_,(cplex, nprojs, ndat))
  ABI_MALLOC(s_projections,(cplex, nprojs, ndat))
  ABI_MALLOC(vnl_projections,(cplex_fac, nprojs, ndat))

  !FIXME Smarter buffer management
#ifdef HAVE_GPU_HIP
  !Work buffer allocated once to save time in HIP (alloc costful)
  !$OMP TARGET ENTER DATA MAP(alloc:projections_,s_projections,vnl_projections)
#endif

  gpu_initialised=1

 end subroutine alloc_work_buffers
!!***

!----------------------------------------------------------------------

 subroutine free_work_buffers()


  !$OMP TARGET EXIT DATA MAP(delete:sij_typ)
  if(allocated(sij_typ)) then
    ABI_FREE(sij_typ)
  end if

  if(allocated(temp_realvec_r)) then
    !$OMP TARGET EXIT DATA MAP(delete:temp_realvec_r,temp_realvec_i)
    ABI_FREE(temp_realvec_r)
    ABI_FREE(temp_realvec_i)
  end if

  !FIXME Smarter buffer management
#ifdef HAVE_GPU_HIP
  !$OMP TARGET EXIT DATA MAP(delete:projections_,s_projections,vnl_projections)
#endif

  if(allocated(projections_)) then
    ABI_FREE(projections_)
    ABI_FREE(s_projections)
    ABI_FREE(vnl_projections)
  end if

  mod__ndat=0
  gpu_initialised=0

 end subroutine free_work_buffers

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_ompgpu/destroy_gpu_nonlop_ompgpu
!! NAME
!! destroy_gemm_nonlop_ompgpu
!!
!! FUNCTION
!! Destruction of the gemm_nonlop_kpt array
!!
!! INPUTS
!! nkpt= number of k-points
!!
!! SOURCE
 subroutine destroy_gemm_nonlop_ompgpu()

  integer :: ikpt
  character(len=200) :: msg

! *************************************************************************

  call free_work_buffers()

 end subroutine destroy_gemm_nonlop_ompgpu
!!***

!----------------------------------------------------------------------

!!****f* m_gemm_nonlop_ompgpu/gemm_nonlop_ompgpu
!! NAME
!! gemm_nonlop_ompgpu
!!
!! FUNCTION
!! Replacement of nonlop. same prototype as nonlop although not all options are implemented.
!!
!! INPUTS
!! [gpu_option] = GPU implementation to use, i.e. cuda, openMP, ... (0=not using GPU)
!!
!! SOURCE
 subroutine gemm_nonlop_ompgpu(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
&                 enl,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&
&                 kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                 mpi_enreg,mpsang,mpssoang,natom,nattyp,ndat,ngfft,nkpgin,nkpgout,nloalg,&
&                 nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,phkxredin,&
&                 phkxredout,ph1d,ph3din,ph3dout,signs,sij,svectout,&
&                 tim_nonlop,ucvol,useylm,vectin,vectout,atom_proj_shift,&
&                 vectproj,gpu_option)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir
  integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin
  integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO
  integer,intent(in) :: paw_opt,signs,tim_nonlop,useylm,atom_proj_shift
  integer,optional,intent(in) :: gpu_option
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
  real(dp),intent(inout),target :: vectin(2,npwin*nspinor*ndat)
  real(dp),intent(inout) :: enlout(nnlout*ndat)
  real(dp),intent(out),target :: svectout(:,:)
  real(dp),intent(inout),target :: vectout(:,:)
  real(dp),intent(inout),optional, ABI_CONTIGUOUS target :: vectproj(:,:,:)
  type(pawcprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5)*ndat)

  ! locals
  complex(dpc), parameter :: cminusone  = (-1._dp,0._dp)
  integer :: ii, ia, idat, igrad, nprojs, ngrads, shift, iatom, nlmn, ierr, ibeg, iend, ikpt
  integer :: cplex, cplex_enl, cplex_fac, proj_shift, grad_shift
  integer :: projs_beg,projs_end,dprojs_beg,dprojs_end
  integer :: nnlout_test
  integer :: iatm, ndgxdt, ndgxdtfac, nd2gxdt, nd2gxdtfac, optder, itypat, ilmn
  integer :: cplex_dgxdt(1), cplex_d2gxdt(1)
  logical :: local_vectproj
  real(dp) :: dgxdt_dum_in(1,1,1,1,1), dgxdt_dum_out(1,1,1,1,1),dgxdt_dum_out2(1,1,1,1,1)
  real(dp) :: d2gxdt_dum_in(1,1,1,1,1), d2gxdt_dum_out(1,1,1,1,1),d2gxdt_dum_out2(1,1,1,1,1)

  real(dp), ABI_CONTIGUOUS pointer :: projections(:,:,:)
  real(dp), allocatable, target :: s_dprojections(:,:,:), vnl_dprojections(:,:,:)
  integer :: ipw, iproj, iblock, nprojs_blk, i1, i2, i
  integer :: nprojs_my_blk
  integer :: rank, nprocs
  logical :: is_last
  real(dp), allocatable :: projs_recv(:,:,:)
  real(dp), ABI_CONTIGUOUS pointer :: projs_(:,:,:),dprojs_(:,:,:)
  integer :: ngrads_tmp
  real(dp), allocatable :: enlk(:),fnlk(:,:),ddkk(:,:),strnlk(:,:)
  integer :: idbeg,idend,dshift,enlout_shift
  real(dp) :: work(6)
  logical :: nld_on_gpu

  logical :: transfer_vectin,transfer_vectout,transfer_svectout
  integer(C_SIZE_T) :: byte_count
  real(dp), ABI_CONTIGUOUS pointer :: vectin_(:,:),vectout_(:,:),svectout_(:,:)

! *************************************************************************

  ! We keep the same interface as nonlop, but we don't use many of those
  ABI_UNUSED((/ffnlin,ffnlout,gmet,kpgin,kpgout/))
  ABI_UNUSED((/ph1d(1,1),ph3din,ph3dout/))
  ABI_UNUSED((/phkxredin,phkxredout,ucvol/))
  ABI_UNUSED((/mgfft,mpsang,mpssoang/))
  ABI_UNUSED((/kptin,kptout/))
  ABI_UNUSED((/idir,nloalg,ngfft,kgin,kgout,ngfft,only_SO,tim_nonlop,gpu_option/))

  ! Check supported options
  if (.not.gemm_nonlop_use_gemm) then
    ABI_BUG('computation not prepared for gemm_nonlop use!')
  end if
  if ( (choice>3.and.choice/=7.and.choice/=5.and.choice/=51.and.signs==2) .or. &
&      (choice>3.and.choice/=7.and.choice/=23.and.signs==1) .or. &
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

  ikpt=gemm_nonlop_ikpt_this_proc_being_treated
  cplex=2;if (istwf_k>1) cplex=1
  cplex_enl=1;if (paw_opt>0) cplex_enl=2*dimenl1/(lmnmax*(lmnmax+1)) ! is enl complex?
  cplex_fac=max(cplex,dimekbq)
  if ((nspinortot==2.or.cplex_enl==2).and.paw_opt>0.and.choice/=7) cplex_fac=2 ! is vnl_projections complex?

  ngrads = gemm_nonlop_kpt(ikpt)%ngrads
  nprojs_my_blk = gemm_nonlop_kpt(ikpt)%nprojs_blk

  ! The number of projectors used for computation may vary among
  ! nonlop calls, from computing on all atoms to a select one for
  ! some perturbations.
  ! In such cases, projs arrays must be recomputed
  nprojs=0
  do itypat=1,ntypat
    nprojs = nprojs + count(indlmn(3,:,itypat)>0)*nattyp(itypat)
  end do
  projs_beg=1; projs_end=nprojs;
  dprojs_beg=1; dprojs_end=nprojs*ngrads;
  if((choice==2 .and. signs==2)) then
    projs_beg=atom_proj_shift+1
    projs_end=projs_beg+nprojs-1
    dprojs_beg=atom_proj_shift*ngrads+1
    dprojs_end=dprojs_beg+nprojs*ngrads-1
  end if

  if(gemm_nonlop_is_distributed) then
    rank = xmpi_comm_rank(gemm_nonlop_block_comm); nprocs = xmpi_comm_size(gemm_nonlop_block_comm)
    is_last = (rank==nprocs-1)
    if(is_last) nprojs_my_blk = gemm_nonlop_kpt(ikpt)%nprojs_last_blk
  end if

  if(choice==1 .or. choice==0 .or. choice==7) ngrads=0
  if(signs==1) then
    ABI_CHECK(ngrads>=3.or.choice/=2 ,"Bad allocation in gemm_nonlop (2)!")
    ABI_CHECK(ngrads>=6.or.choice/=3 ,"Bad allocation in gemm_nonlop (3)!")
    ABI_CHECK(ngrads>=9.or.choice/=23,"Bad allocation in gemm_nonlop (23)!")
  else if(signs==2) then
    ABI_CHECK(ngrads==1.or.(choice==1.or.choice==0.or.choice==7) ,"Bad allocation in gemm_nonlop !")
  end if

  ! Allocate and copy GPU buffers if user doesn't manage them
  transfer_vectin=.not. xomp_target_is_present(c_loc(vectin)) &
      .and. ((cpopt < 2 .and. choice < 2) .or. (cpopt <= 3 .and. choice >= 2) &
      .or. (choice/=7 .and.paw_opt >=3))
  transfer_vectout=.not. xomp_target_is_present(c_loc(vectout)) &
      .and. (signs==2 .and. (paw_opt == 0 .or. paw_opt == 1 .or. paw_opt == 4))
  transfer_svectout=.not. xomp_target_is_present(c_loc(svectout)) &
      .and. (signs==2 .and. (paw_opt == 3 .or. paw_opt == 4))
  !$OMP TARGET ENTER DATA MAP(to:vectin)      IF(transfer_vectin)
  !$OMP TARGET ENTER DATA MAP(alloc:vectout)  IF(transfer_vectout)
  !$OMP TARGET ENTER DATA MAP(alloc:svectout) IF(transfer_svectout)

  !FIXME These seemingly useless pointers are used in BLAS operations for
  ! working around a bug in AOMP LLVM misreading mapped device pointers.
  ! I chose to generalise the workaround to avoid dupplicating each
  ! BLAS call specifically for handling AOMP LLVM.
  vectin_   => vectin
  vectout_  => vectout
  svectout_ => svectout

  !$OMP TARGET ENTER DATA MAP(to:atindx1,indlmn,enl)
  if(gpu_initialised == 0 .or. mod__ndat /= ndat*nspinor .or. nprojs /= mod__nprojs) then
    call alloc_work_buffers(cplex, cplex_fac,&
&        nspinor*ndat, nprojs, ntypat, lmnmax, MAX(npwin,npwout))
  end if

  if(paw_opt>=2 .and. choice > 0 .and. choice /= 7) then
    if (cplex_enl==1) then
      do itypat=1, ntypat
        nlmn=count(indlmn(3,:,itypat)>0)
        do ilmn=1,nlmn*(nlmn+1)/2
          sij_typ(ilmn,itypat)=sij(ilmn,itypat)
        end do
      end do
    else
      do itypat=1, ntypat
        nlmn=count(indlmn(3,:,itypat)>0)
        do ilmn=1,nlmn*(nlmn+1)/2
          sij_typ(ilmn,itypat)=sij(2*ilmn-1,itypat)
        end do
      end do
    end if
  end if
  !$OMP TARGET UPDATE TO(sij_typ)

  if(current_ikpt_in_gpu /= ikpt) then
    call refresh_gemm_nonlop_kpt_ompgpu(ikpt)
  end if

  projs_ => current_ikpt_projs(:,:,projs_beg:projs_end)
  dprojs_ => current_ikpt_dprojs(:,:,dprojs_beg:dprojs_end)

  ! If vectproj is provided, use it for further calculations, use static array otherwise
  projections => projections_
  local_vectproj=.false.
  if(PRESENT(vectproj)) then
    if(size(vectproj)>1) local_vectproj=.true.
  end if
  if (local_vectproj) projections => vectproj

#ifdef HAVE_GPU_CUDA
  !Work buffers allocated at each call to save memory in CUDA
  !$OMP TARGET ENTER DATA MAP(alloc:s_projections,vnl_projections)
  if(.not. local_vectproj) then
    !$OMP TARGET ENTER DATA MAP(alloc:projections_)
  end if
#endif

  if(gemm_nonlop_is_distributed) then
    ABI_MALLOC(projs_recv, (cplex, npwin, MAX(ngrads,1)*gemm_nonlop_kpt(ikpt)%nprojs_last_blk))
    !$OMP TARGET ENTER DATA MAP(alloc:projs_recv)
  end if

  ! These will store the non-local factors for vectin, svectout and vectout respectively
  if(cpopt < 2) then
    projections=zero
    !$OMP TARGET UPDATE TO(projections)
    !call gpu_set_to_zero(projections,   int(cplex,c_size_t)*nprojs*ndat*nspinor)
  end if
  call gpu_set_to_zero(s_projections,   int(cplex,c_size_t)*nprojs*ndat*nspinor)
  call gpu_set_to_zero(vnl_projections, int(cplex_fac,c_size_t)*nprojs*ndat*nspinor)

  if (ngrads>0) then
    ABI_MALLOC(dprojections,(cplex, ngrads*nprojs, nspinor*ndat))
    !$OMP TARGET ENTER DATA MAP(alloc:dprojections)
    ! Working buffers for storing derivative (used for response function at least)
    if (choice > 1) then
      ABI_MALLOC(s_dprojections,(cplex, ngrads*nprojs,nspinor*ndat))
      ABI_MALLOC(vnl_dprojections,(cplex_fac, ngrads*nprojs,nspinor*ndat))
      !$OMP TARGET ENTER DATA MAP(alloc:s_dprojections,vnl_dprojections)
    end if
    if(cpopt < 4) then
      call gpu_set_to_zero(dprojections, int(cplex,c_size_t)*ngrads*nprojs*ndat*nspinor)
    end if
    if(allocated(s_dprojections)) then
      s_dprojections=zero
      !$OMP TARGET UPDATE TO(s_dprojections)
      !call gpu_set_to_zero(s_dprojections,   int(cplex,c_size_t)*ngrads*nprojs*ndat*nspinor)
    end if
    if(allocated(vnl_dprojections)) then
      vnl_dprojections=zero
      !$OMP TARGET UPDATE TO(vnl_dprojections)
      !call gpu_set_to_zero(vnl_dprojections, int(cplex_fac,c_size_t)*ngrads*nprojs*ndat*nspinor)
    end if
  end if

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

  if(signs == 1 .and. choice > 0) then
    enlout=zero
    !$OMP TARGET ENTER DATA MAP(to:enlout)
    ABI_MALLOC(enlk,(ndat))
    enlk=zero
    ABI_MALLOC(fnlk,(3*natom,ndat))
    ABI_MALLOC(ddkk,(6,ndat))
    ABI_MALLOC(strnlk,(6,ndat))
    !$OMP TARGET ENTER DATA MAP(to:enlk)
  end if

  !$OMP TASKWAIT
  if(cpopt >= 2) then
    ! retrieve from cprjin
    if(.not. local_vectproj .and. cpopt/=3) then
      !$OMP PARALLEL DO PRIVATE(shift,idat,iatom,nlmn)
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = 1, natom
          nlmn = cprjin(iatom, idat)%nlmn
          projections(1:cplex, shift+1:shift+nlmn, idat) = cprjin(iatom, idat)%cp(1:cplex, 1:nlmn)
          shift = shift + nlmn
        end do
      end do
      !$OMP TARGET UPDATE TO(projections)
    end if
    if(cpopt==4.and.allocated(dprojections)) then
      ABI_CHECK(cprjin(1,1)%ncpgr>=ngrads,"cprjin%ncpgr not correct! (1)")
      !$OMP PARALLEL DO PRIVATE(shift,idat,iatom,igrad,nlmn)
      do idat=1, ndat*nspinor
        shift = 0
        do iatom = 1, natom
          nlmn  = cprjin(iatom, idat)%nlmn
          do ilmn=1,nlmn
            do igrad=1,ngrads
              dprojections(1:cplex, shift + igrad, idat) = &
                cprjin(iatom, idat)%dcp(1:cplex,igrad,ilmn)
            end do
            shift = shift + ngrads
          end do
        end do
      end do
      !$OMP TARGET UPDATE TO(dprojections)
    end if
  end if ! cpopt

  if(cpopt<=1.or.(cpopt<=3.and.(choice==2.or.choice==3.or.choice==5.or.choice==51.or.choice==23.or.choice==54.or.choice==55.or.choice==4))) then

    call opernla_gemm(choice,cplex,cplex_dgxdt,cplex_d2gxdt,d2gxdt_dum_in,dprojections,projections,&
    &       idir,istwf_k,mpi_enreg,nd2gxdt,ngrads,&
    &       npwin,nspinor,signs,ndat,rank,&
    &       cpopt,nprocs,&
    &       nprojs,gemm_nonlop_kpt(ikpt)%nprojs_blk,nprojs_my_blk,gemm_nonlop_kpt(ikpt)%nprojs_last_blk,&
    &       vectin,projs_,dprojs_,&
    &       gemm_nonlop_kpt(ikpt)%projs_r,gemm_nonlop_kpt(ikpt)%projs_i,&
    &       gemm_nonlop_kpt(ikpt)%dprojs_r,gemm_nonlop_kpt(ikpt)%dprojs_i,&
    &       temp_realvec_r,&
    &       projs_recv,projs_recv,&
    &       gpu_option,gemm_nonlop_is_distributed)

    if(cpopt >= 0) then
      ! store in cprjin
      if(.not. local_vectproj .and. cpopt/=3) then
        !$OMP TARGET UPDATE FROM(projections)
        !$OMP PARALLEL DO PRIVATE(shift,idat,iatom,nlmn)
        do idat=1, ndat*nspinor
          shift = 0
          do iatom = 1, natom
            nlmn = cprjin(iatom, idat)%nlmn
            cprjin(iatom, idat)%cp(1:cplex, 1:nlmn) = projections(1:cplex, shift+1:shift+nlmn, idat)
            shift = shift + nlmn
          end do
        end do
      end if
      if(cpopt==1 .or. cpopt==3) then
        ABI_CHECK(cprjin(1,1)%ncpgr>=ngrads,"cprjin%ncpgr not correct! (2)")
        !$OMP TARGET UPDATE FROM(dprojections)
        !$OMP PARALLEL DO PRIVATE(shift,idat,iatom,igrad,nlmn)
        do idat=1, ndat*nspinor
          shift = 0
          do iatom = 1, natom
            nlmn = cprjin(iatom, idat)%nlmn
            do ilmn=1,nlmn
              do igrad=1,ngrads
                cprjin(iatom, idat)%dcp(1:cplex,igrad,ilmn) = &
                &                   dprojections(1:cplex, shift + igrad, idat)
              end do
              shift = shift + ngrads
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
      ndgxdt = ngrads
      ndgxdtfac = ngrads
      nd2gxdt = 0
      nd2gxdtfac = 0
      optder = 0;if (ndgxdtfac>0 .and. signs == 2) optder=1
      cplex_dgxdt(:) = 1;
      if(ndgxdt > 0) then
       if (choice==5.or.choice==51) cplex_dgxdt(:) = 2
      end if

      shift = 0
      do itypat=1, ntypat
        nlmn=count(indlmn(3,:,itypat)>0)

        ibeg = shift+1
        iend = shift+nattyp(itypat)*nlmn

        call opernlc_ylm_ompgpu(atindx1,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_enl,cplex_fac,&
  &         dprojections,&
  &         vnl_dprojections,&
  &         s_dprojections,&
  &         d2gxdt_dum_in,d2gxdt_dum_out,d2gxdt_dum_out2,dimenl1,dimenl2,dimekbq,enl,&
  &         projections,&
  &         vnl_projections,&
  &         s_projections,&
  &         iatm,indlmn,itypat,lambda,mpi_enreg,natom,ndgxdt,ndgxdtfac,nd2gxdt,nd2gxdtfac,&
  &         nattyp(itypat),nlmn,nspinor,nspinortot,optder,paw_opt,sij_typ,ndat,ibeg-1,iend,nprojs,ntypat)

        shift = shift + nattyp(itypat)*nlmn
        iatm = iatm+nattyp(itypat)
      end do
    else
      !$OMP TARGET DATA USE_DEVICE_PTR(s_projections,projections)
      call copy_gpu_to_gpu(c_loc(s_projections), &
           &               c_loc(projections), &
           &               INT(cplex, c_size_t) * nprojs * nspinor * ndat * dp)
      !$OMP END TARGET DATA
    end if

    ! opernlb
    if(signs==2) then
      call opernlb_gemm(choice,cplex,cplex_dgxdt,cplex_d2gxdt,cplex_fac,&
      &       d2gxdt_dum_out,d2gxdt_dum_out,&
      &       vnl_dprojections,s_dprojections,&
      &       vnl_projections,s_projections,&
      &       idir,istwf_k,mpi_enreg,nd2gxdt,nd2gxdtfac,ndgxdt,ndgxdtfac,&
      &       npwout,nspinor,signs,ndat,rank,&
      &       cpopt,nprocs,paw_opt,&
      &       nprojs,gemm_nonlop_kpt(ikpt)%nprojs_blk,nprojs_my_blk,gemm_nonlop_kpt(ikpt)%nprojs_last_blk,&
      &       vectin_,vectout_,svectout_,&
      &       projs_,&
      &       dprojs_,&
      &       gemm_nonlop_kpt(ikpt)%projs_r,gemm_nonlop_kpt(ikpt)%projs_i,&
      &       gemm_nonlop_kpt(ikpt)%dprojs_r,gemm_nonlop_kpt(ikpt)%dprojs_i,&
      &       temp_realvec_r,temp_realvec_i,&
      &       projs_recv,projs_recv,&
      &       gpu_option,gemm_nonlop_is_distributed)
    end if

    ! opernld
    if(signs==1) then
#if 1
      nld_on_gpu = .true.
      call opernld_ylm_allwf_ompgpu(choice,cplex,cplex_fac,&
      &       dprojections,vnl_dprojections,s_dprojections,d2gxdt_dum_in,&
      &       enlk,enlout,projections,vnl_projections,s_projections,&
      &       ndat,nd2gxdt,ndgxdt,&
      &       ndgxdtfac,indlmn,ntypat,lmnmax,nprojs,nnlout,nspinor,paw_opt,&
      &       nattyp)
#else
      shift=0; dshift=0; iatm=1
      nld_on_gpu = .false.
      !$OMP TARGET UPDATE FROM(dprojections,projections,vnl_projections)
      do itypat=1, ntypat
        nlmn=count(indlmn(3,:,itypat)>0)

        ibeg = shift+1
        iend = shift+nattyp(itypat)*nlmn

        idbeg = dshift+1
        idend = dshift+nattyp(itypat)*nlmn*ngrads

        do idat=1,ndat
          call opernld_ylm             (choice,cplex,cplex_fac,ddkk(:,idat),&
          &       dprojections    (:, idbeg:idend, 1+nspinor*(idat-1):nspinor*idat),&
          &       vnl_dprojections(:, idbeg:idend, 1+nspinor*(idat-1):nspinor*idat),&
          &       s_dprojections  (:, idbeg:idend, 1+nspinor*(idat-1):nspinor*idat),&
          &       d2gxdt_dum_in,&
          &       enlk(idat),enlout(nnlout*(idat-1)+1:nnlout*idat),fnlk(:,idat),&
          &       projections    (:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
          &       vnl_projections(:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
          &       s_projections  (:, ibeg:iend, 1+nspinor*(idat-1):nspinor*idat),&
          &       iatm,natom,1,nd2gxdt,ndgxdt,ndgxdtfac,&
          &       nattyp(itypat),nlmn,nnlout,nspinor,paw_opt,strnlk(:,idat))
        end do

        shift = shift + nattyp(itypat)*nlmn
        dshift = dshift + nattyp(itypat)*nlmn*ngrads
        iatm = iatm+nattyp(itypat)
      end do
#endif

      ! Reduction in case of parallelism
      if (mpi_enreg%paral_spinor==1) then
        if (size(enlout)>0) then
          !$OMP TARGET UPDATE FROM(enlout) if(nld_on_gpu)
          call xmpi_sum(enlout,mpi_enreg%comm_spinor,ierr)
          !$OMP TARGET UPDATE TO(enlout) if(nld_on_gpu)
        end if
        if (choice==3.or.choice==23) then
          !$OMP TARGET UPDATE FROM(enlk) if(nld_on_gpu)
          call xmpi_sum(enlk,mpi_enreg%comm_spinor,ierr)
        end if
      end if

      ! Derivatives wrt strain
      !  - Convert from reduced to cartesian coordinates
      !  - Substract volume contribution
      if ((choice==3.or.choice==23).and.paw_opt<=3) then
        !$OMP TARGET UPDATE FROM(enlout,enlk) if(nld_on_gpu)
        do idat=1,ndat
          enlout_shift=(idat-1)*nnlout
          call strconv(enlout(enlout_shift+1:enlout_shift+6),gprimd,work)
          enlout(enlout_shift+1:enlout_shift+3)=(work(1:3)-enlk(idat))
          enlout(enlout_shift+4:enlout_shift+6)= work(4:6)
        end do
        !$OMP TARGET UPDATE TO(enlout) if(nld_on_gpu)
      end if

    end if !opernld

  end if ! choice>0

  ! Retrieve and release allocated buffers
  !$OMP TARGET EXIT DATA MAP(delete:vectin) IF(transfer_vectin)
  !$OMP TARGET EXIT DATA MAP(from:vectout)   IF(transfer_vectout)
  !$OMP TARGET EXIT DATA MAP(from:svectout)  IF(transfer_svectout)

#ifdef HAVE_GPU_CUDA
  !$OMP TARGET EXIT DATA MAP(delete:s_projections,vnl_projections)
  if(.not. local_vectproj) then
    !$OMP TARGET EXIT DATA MAP(delete:projections_)
  end if
#endif

! Release memory
  if(signs == 1 .and. choice > 0) then
    !$OMP TARGET EXIT DATA MAP(delete:enlk)
    ABI_FREE(enlk)
    ABI_FREE(fnlk)
    ABI_FREE(strnlk)
    ABI_FREE(ddkk)
    !$OMP TARGET UPDATE FROM(enlout) if(nld_on_gpu)
    !$OMP TARGET EXIT DATA MAP(delete:enlout)
  end if

  if(gemm_nonlop_is_distributed) then
    !$OMP TARGET EXIT DATA MAP(delete:projs_recv)
    ABI_FREE(projs_recv)
  end if
  !$OMP TARGET EXIT DATA MAP(delete:atindx1,indlmn,enl)
  if (allocated(dprojections)) then
    !$OMP TARGET EXIT DATA MAP(delete:dprojections)
    ABI_FREE(dprojections)
  end if
  if (allocated(s_dprojections)) then
    !$OMP TARGET EXIT DATA MAP(release:s_dprojections)
    ABI_FREE(s_dprojections)
  end if
  if (allocated(vnl_dprojections)) then
    !$OMP TARGET EXIT DATA MAP(release:vnl_dprojections)
    ABI_FREE(vnl_dprojections)
  end if

 end subroutine gemm_nonlop_ompgpu
!***


#else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! stubs for compiling with OpenMP GPU offload disabled.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine gemm_nonlop_ompgpu(atindx1,choice,cpopt,cprjin,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,&
&                 enl,enlout,ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,&
&                 kgin,kgout,kpgin,kpgout,kptin,kptout,lambda,lmnmax,matblk,mgfft,&
&                 mpi_enreg,mpsang,mpssoang,natom,nattyp,ndat,ngfft,nkpgin,nkpgout,nloalg,&
&                 nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO,paw_opt,phkxredin,&
&                 phkxredout,ph1d,ph3din,ph3dout,signs,sij,svectout,&
&                 tim_nonlop,ucvol,useylm,vectin,vectout,atom_proj_shift,&
&                 vectproj,gpu_option)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir
  integer,intent(in) :: istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin
  integer,intent(in) :: nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO
  integer,intent(in) :: paw_opt,signs,tim_nonlop,useylm,atom_proj_shift
  integer,optional,intent(in) :: gpu_option
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
  real(dp),intent(inout),target :: vectin(2,npwin*nspinor*ndat)
  real(dp),intent(inout) :: enlout(nnlout*ndat)
  real(dp),intent(out),target :: svectout(2,npwout*nspinor*(paw_opt/3)*ndat)
  real(dp),intent(inout),target :: vectout(2,npwout*nspinor*ndat) !vz_i
  real(dp),intent(inout),optional,target :: vectproj(:,:,:)
  type(pawcprj_type),intent(inout) :: cprjin(natom,nspinor*((cpopt+5)/5)*ndat)

  ABI_UNUSED((/choice,cpopt,dimenl1,dimenl2,dimekbq,dimffnlin,dimffnlout,idir/))
  ABI_UNUSED((/istwf_k,lmnmax,matblk,mgfft,mpsang,mpssoang,natom,ndat,nkpgin/))
  ABI_UNUSED((/nkpgout,nnlout,npwin,npwout,nspinor,nspinortot,ntypat,only_SO/))
  ABI_UNUSED((/paw_opt,signs,tim_nonlop,useylm,gpu_option/))
  ABI_UNUSED((/atindx1,indlmn,kgin,kgout,nattyp,ngfft,nloalg/))
  ABI_UNUSED((/enl,ffnlin,ffnlout,gmet,gprimd,kpgin,kpgout,kptin,kptout,phkxredin,phkxredout/))
  ABI_UNUSED((/ucvol,lambda,sij,ph1d(1,1),ph3din,ph3dout,vectin,enlout,svectout,vectout,vectproj/))
  ABI_UNUSED_A(cprjin)
  ABI_UNUSED_A(mpi_enreg)
  ABI_BUG("Unhandled configuration for OpenMP GPU immplementation")

 end subroutine gemm_nonlop_ompgpu

#endif

end module m_gemm_nonlop_ompgpu
!!***
