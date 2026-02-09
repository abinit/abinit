!!****f* ABINIT/m_slicewf
!! NAME
!! m_slicewf
!!
!! FUNCTION
!! This module contains a routine updating the whole wave functions at a given k-point,
!! using the Spectrum Slicing method (2021 implementation using xG abstraction layer)
!! for a given spin-polarization, from a fixed Hamiltonian but might also simply compute 
!! eigenvectors and eigenvalues at this k point. it will also update the matrix elements 
!! of the Hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2025 ABINIT group (BS, IML)
!! This file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

! nvtx related macro definition
#include "nvtx_macros.h"

module m_slicewf

 use defs_abitypes
 use defs_basis
 use m_abicore
 use m_errors
 use m_fstrings
 use m_time

 use m_slice
 !use m_chebfi
 !use m_chebfi2
 use m_invovl

 use m_cgtools,     only : dotprod_g
 use m_dtset,       only : dataset_type

 use m_hamiltonian, only : gs_hamiltonian_type
 use m_pawcprj,     only : pawcprj_type
 use m_nonlop,      only : nonlop
 use m_prep_kgb,    only : prep_getghc, prep_nonlop
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_getghc,      only : multithreaded_getghc
 use m_gemm_nonlop_projectors , only : gemm_nonlop_use_gemm

 use m_xg
 use m_xgTransposer

#if defined(HAVE_GPU_MARKERS)
 use m_nvtx_data
#endif

#if defined(HAVE_GPU)
 use m_gpu_toolbox
#endif

#if defined(HAVE_YAKL)
 use gator_mod
#endif

 use, intrinsic :: iso_c_binding, only: c_associated,c_loc,c_ptr,c_f_pointer,c_double,c_size_t

 use m_xmpi
 use m_xomp
#ifdef HAVE_OPENMP
 use omp_lib
#endif

 implicit none

 private

 integer, parameter :: l_tim_getghc=7
 real(dp), parameter :: inv_sqrt2 = 1/sqrt2

! For use in getghc_gsc1
 integer, save :: l_cpopt
 integer, save :: l_icplx
 integer, save :: l_npw
 integer, save :: l_nband_filter
 integer, save :: l_nspinor
 logical, save :: l_paw
 integer, save :: l_prtvol
 integer, save :: l_sij_opt
 integer, save :: l_paral_kgb
 integer, save :: l_useria
 integer, save :: l_block_sliced

 type(mpi_type),pointer,save :: l_mpi_enreg
 type(gs_hamiltonian_type),pointer,save :: l_gs_hamk

 integer, parameter :: DEBUG_ROWS = 5
 integer, parameter :: DEBUG_COLUMNS = 5

 public :: slicewf

 CONTAINS  !========================================================================================
!!***

!!****f* m_slicewf/slicewf
!! NAME
!! slicewf
!!
!! FUNCTION
!! This routine updates the whole wave functions set at a given k-point,
!! using the Spectrum Slicing method (2025 version using xG abstraction layer)
!!
!! INPUTS
!!  dtset= input variables for this dataset
!!  mpi_enreg= MPI-parallelisation information
!!  nband= number of bands at this k point
!!  npw= number of plane waves at this k point
!!  nspinor= number of spinorial components of the wavefunctions
!!  prtvol= control print volume and debugging
!!
!! OUTPUT
!!  eig(nband)= eigenvalues (hartree) for all bands
!!  enl_out(nband)= contribution of each band to the nl part of energy
!!  resid(nband)= residuals for each band
!!
!! SIDE EFFECTS
!!  cg(2,npw*nspinor*nband)= planewave coefficients of wavefunctions
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!
!! SOURCE

subroutine slicewf(cg,dtset,eig,occ,enl_out,gs_hamk,mpi_enreg,&
&                    nband,npw,nspinor,prtvol,resid)

 implicit none

 ! Arguments ------------------------------------
 integer,intent(in) :: nband,npw,prtvol,nspinor
 type(mpi_type),target,intent(in) :: mpi_enreg
 real(dp),target,intent(inout) :: cg(2,npw*nspinor*nband)
 real(dp),target,intent(out) :: resid(nband)
 real(dp),intent(out) :: enl_out(nband)
 real(dp),target,intent(out) :: eig(nband)
 real(dp),target,intent(in) :: occ(nband)
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),target,intent(inout) :: gs_hamk

 ! Local variables-------------------------------
 ! scalars
 integer, parameter :: tim_slicewf = 2160
 integer, parameter :: tim_nonlop = 1753
 integer :: iband,shift,space,blockdim
 integer :: spacedim,spacecom,gpu_option
 integer :: me_g0,me_g0_fft
 !integer(kind=c_size_t) :: localMem
 type(slice_t) :: slice
 type(xgBlock_t) :: xgx0,xgeigen,xgresidu
 ! arrays
! real(dp) :: cg_temp(2,npw*nspinor*nband)
 real(dp) :: tsec(2)
 !integer(kind=c_size_t) :: sliceMem(2) 
 real(dp), allocatable :: l_gvnlxc(:,:)
 
 ! Parameters for nonlop call in NC
 integer,parameter :: choice=1, paw_opt=0, signs=1
 real(dp) :: gsc_dummy(1,1)
 type(pawcprj_type) :: cprj_dum(gs_hamk%natom,1)

! *********************************************************************

!################ INITIALIZATION  #####################################
!######################################################################

  call timab(tim_slicewf,1,tsec)

!Set module variables
 l_paw = (gs_hamk%usepaw==1)
 l_cpopt=-1;l_sij_opt=0;if (l_paw) l_sij_opt=1
 l_npw = npw
 l_nspinor = nspinor
 l_prtvol = prtvol
 l_mpi_enreg => mpi_enreg
 l_gs_hamk => gs_hamk
 l_nband_filter = nband
 l_paral_kgb = dtset%paral_kgb
 l_block_sliced = dtset%invovl_blksliced

!Variables
 spacedim = l_npw*l_nspinor
 spacecom = l_mpi_enreg%comm_bandspinorfft
 gpu_option = dtset%gpu_option
 blockdim=l_mpi_enreg%nproc_band*l_mpi_enreg%bandpp
 !for debug
 l_useria=dtset%useria

!Depends on istwfk
 if ( gs_hamk%istwf_k > 1 ) then ! Real only
   space = SPACE_CR
   l_icplx = 2
 else ! complex
   space = SPACE_C
   l_icplx = 1
 end if

 me_g0 = -1
 me_g0_fft = -1
 if (space==SPACE_CR) then
   me_g0 = 0
   me_g0_fft = 0
   if (gs_hamk%istwf_k == 2) then
     if (l_mpi_enreg%me_g0 == 1) me_g0 = 1
     if (l_mpi_enreg%me_g0_fft == 1) me_g0_fft = 1
   end if
 end if
 
 ! Temporary cg object used in Slicing
! cg_temp(:,:) = cg(:,:)

!#ifdef HAVE_OPENMP_OFFLOAD
! !$OMP TARGET ENTER DATA MAP(to:cg,cg_temp,eig,resid) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
!#endif
#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET ENTER DATA MAP(to:cg,eig,resid) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif

 !call xgBlock_map(xgx0,cg_temp,space,spacedim,nband,comm=spacecom,me_g0=me_g0,gpu_option=gpu_option)
 call xgBlock_map(xgx0,cg,space,spacedim,nband,comm=spacecom,me_g0=me_g0,gpu_option=gpu_option)

 call xgBlock_map_1d(xgeigen,eig,SPACE_R,nband,gpu_option=gpu_option)

 call xgBlock_map_1d(xgresidu,resid,SPACE_R,nband,gpu_option=gpu_option)

 write(std_out,*) 'calling slice_init'
 write(std_out,*) 'dtset%paral_kgb=', dtset%paral_kgb
 write(std_out,*) 'spacecom=', spacecom, xmpi_comm_size(spacecom)
 write(std_out,*) 'spacedim=', spacedim

 call slice_init(slice,dtset%nslice,nband,spacedim,dtset%tolwfr_diago,dtset%paral_kgb,&
        dtset%paral_slice,dtset%mdeg_filter,dtset%tolfilter,dtset%ecut,l_mpi_enreg%bandpp,&
        space,spacecom,me_g0,me_g0_fft,l_paw,l_mpi_enreg%comm_spinorfft,l_mpi_enreg%comm_band,&
        dtset%spectral_cut,l_gs_hamk%gpu_option,gpu_kokkos_nthrd=dtset%gpu_kokkos_nthrd,&
        gpu_thread_limit=dtset%gpu_thread_limit)
 
 call slice_allschedule(slice, xgx0, getghc_gsc1, getBm1X, xgeigen, nspinor)
 
 if (dtset%nslice>1) then
    ! Release collective cg_temp memory from GPU, will only use active task memory
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET UPDATE FROM(cg) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
    !$OMP TARGET EXIT DATA MAP(delete:cg) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif
!#ifdef HAVE_OPENMP_OFFLOAD
!    !$OMP TARGET UPDATE FROM(cg_temp) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
!    !$OMP TARGET EXIT DATA MAP(delete:cg_temp) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
!#endif
 end if

!################    RUUUUUUUN    #####################################
!######################################################################

 call slice_run(slice, getghc_gsc1, getBm1X, xgeigen, xgresidu, nspinor)

 if (dtset%nslice>1) then
    ! Retransfer collective cg_temp memory on GPU
!#ifdef HAVE_OPENMP_OFFLOAD
!   !$OMP TARGET ENTER DATA MAP(to:cg_temp) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
!#endif
!   call xgBlock_map(xgx0,cg_temp,space,spacedim,nband,comm=spacecom,me_g0=me_g0,gpu_option=gpu_option)
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET ENTER DATA MAP(to:cg) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif
   call xgBlock_map(xgx0,cg,space,spacedim,nband,comm=spacecom,me_g0=me_g0,gpu_option=gpu_option)
 end if

 call slice_allmerge(slice, xgx0, xgeigen, xgresidu)

 ! Free slice memory
 call slice_free(slice)

 if ( .not. l_paw ) then
   call timab(tim_nonlop,1,tsec)
#ifdef FC_CRAY
   ABI_MALLOC(l_gvnlxc,(1,1))
#else
   ABI_MALLOC(l_gvnlxc,(0,0))
#endif
   !end if

   ABI_NVTX_START_RANGE(NVTX_SLICE_NONLOP)
   !Call nonlop
   if (l_paral_kgb==0) then

     call nonlop(choice,l_cpopt,cprj_dum,enl_out,l_gs_hamk,0,eig,mpi_enreg,nband,1,paw_opt,&
&                signs,gsc_dummy,l_tim_getghc,cg,l_gvnlxc)

   else
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET UPDATE FROM(cg) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif
     do iband=1,nband/blockdim
       shift = (iband-1)*blockdim*l_npw*l_nspinor
       call prep_nonlop(choice,l_cpopt,cprj_dum, &
&        enl_out((iband-1)*blockdim+1:iband*blockdim),l_gs_hamk,0,&
&        eig((iband-1)*blockdim+1:iband*blockdim),blockdim,mpi_enreg,1,paw_opt,signs,&
&        gsc_dummy,l_tim_getghc, &
&        cg(:,shift+1:shift+blockdim*l_npw*l_nspinor),&
!&        l_gvnlxc(:,shift+1:shift+blockdim*l_npw*l_nspinor),&
&        l_gvnlxc(:,:),&
&        already_transposed=.false.)
     end do
   end if
   ABI_NVTX_END_RANGE()
   ABI_FREE(l_gvnlxc)
   call timab(tim_nonlop,2,tsec)
 end if

!Free slice
 call slice_free(slice)

#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET UPDATE FROM(cg,eig,resid) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
 !$OMP TARGET EXIT DATA MAP(delete:cg,eig,resid) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif
!#ifdef HAVE_OPENMP_OFFLOAD
! !$OMP TARGET UPDATE FROM(cg,eig,resid) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
! !$OMP TARGET EXIT DATA MAP(delete:cg_temp,cg,eig,resid) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
!#endif

 call timab(tim_slicewf,2,tsec)

 DBG_EXIT("COLL")

end subroutine slicewf
!!***

!----------------------------------------------------------------------

!!****f* m_slicewf/getghc_gsc1
!! NAME
!! getghc_gsc1
!!
!! FUNCTION
!! This routine computes H|C> and possibly S|C> for a given wave function C.
!!  It acts as a driver for getghc, taken into account parallelism, multithreading, etc.
!!
!! SIDE EFFECTS
!!  X  <type(xgBlock_t)>= memory block containing |C>
!!  AX <type(xgBlock_t)>= memory block containing H|C>
!!  BX <type(xgBlock_t)>= memory block containing S|C>
!!
!! SOURCE

subroutine getghc_gsc1(X,AX,BX)

 implicit none

!Arguments ------------------------------------
 type(xgBlock_t), intent(inout) :: X
 type(xgBlock_t), intent(inout) :: AX
 type(xgBlock_t), intent(inout) :: BX
 integer         :: blockdim
 integer         :: spacedim
 type(pawcprj_type) :: cprj_dum(l_gs_hamk%natom,1)

!Local variables-------------------------------
!scalars
 real(dp) :: eval
!arrays
 real(dp), pointer :: cg(:,:)
 real(dp), pointer :: ghc(:,:)
 real(dp), pointer :: gsc(:,:)
 real(dp)          :: l_gvnlxc(1,1)

! *********************************************************************

 ABI_NVTX_START_RANGE(NVTX_GETGHC)

 call xgBlock_getSize(X,spacedim,blockdim)
 call xgBlock_check(X,AX)
 call xgBlock_check(X,BX)

 write(901,*) 'debug A"', spacedim, blockdim
 write(901,*) 'debug A"', xgBlock_getid(AX)
 flush(901)

 call xgBlock_reverseMap(X,cg,rows=1,cols=spacedim*blockdim)
 call xgBlock_reverseMap(AX,ghc,rows=1,cols=spacedim*blockdim)
 call xgBlock_reverseMap(BX,gsc,rows=1,cols=spacedim*blockdim)

 call multithreaded_getghc(l_cpopt,cg,cprj_dum,ghc,gsc,&
   l_gs_hamk,l_gvnlxc,eval,l_mpi_enreg,blockdim,l_prtvol,l_sij_opt,l_tim_getghc,0)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
 call gpu_device_synchronize()
#endif

 if ( .not. l_paw ) call xgBlock_copy(X,BX)

 write(901,*) 'debug A', spacedim, blockdim
 write(901,*) 'debug A', xgBlock_getid(AX)
 flush(901)

 ABI_NVTX_END_RANGE()

end subroutine getghc_gsc1
!!***

!----------------------------------------------------------------------

!!****f* m_slicewf/getBm1X
!! NAME
!! getBm1X
!!
!! FUNCTION
!! This routine computes S^-1|C> for a given wave function C.
!!  It acts as a driver for apply_invovl.
!!
!! SIDE EFFECTS
!!  X  <type(xgBlock_t)>= memory block containing |C>
!!  Bm1X <type(xgBlock_t)>= memory block containing S^-1|C>
!!
!! SOURCE

subroutine getBm1X(X,Bm1X)

 implicit none

!Arguments ------------------------------------
 type(xgBlock_t), intent(inout) :: X
 type(xgBlock_t), intent(inout) :: Bm1X

!Local variables-------------------------------
!scalars
 integer :: blockdim
 integer :: spacedim
!arrays
 real(dp), pointer :: ghc_filter(:,:)
 real(dp), pointer :: gsm1hc_filter(:,:)
 type(pawcprj_type), allocatable :: cwaveprj_next(:,:) !dummy

! *********************************************************************

 ! working bandpp will be equal to blockdim
 call xgBlock_getSize(X,spacedim,blockdim)

 if(l_paw) then

   call xgBlock_reverseMap(X,ghc_filter,rows=1,cols=spacedim*blockdim)
   call xgBlock_reverseMap(Bm1X,gsm1hc_filter,rows=1,cols=spacedim*blockdim)

   !cwaveprj_next is dummy
   if(gemm_nonlop_use_gemm) then
     ABI_MALLOC(cwaveprj_next, (1,1))
   else
     ABI_MALLOC(cwaveprj_next, (l_gs_hamk%natom,l_nspinor*blockdim))
     call pawcprj_alloc(cwaveprj_next,0,l_gs_hamk%dimcprj)
   end if

   ABI_NVTX_START_RANGE(NVTX_INVOVL)
   call apply_invovl(l_gs_hamk, ghc_filter(:,:), gsm1hc_filter(:,:), cwaveprj_next(:,:), &
       spacedim/l_nspinor, blockdim, l_mpi_enreg, l_nspinor, l_block_sliced)
   ABI_NVTX_END_RANGE()

   call pawcprj_free(cwaveprj_next)
   ABI_FREE(cwaveprj_next)

 else

   call xgBlock_copy(X,Bm1X)

 end if

end subroutine getBm1X
!!***

!----------------------------------------------------------------------

end module m_slicewf
!!***
