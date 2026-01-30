!!****f* ABINIT/m_chebfiwf
!! NAME
!! m_chebfiwf
!!
!! FUNCTION
!! This module contains a routine updating the whole wave functions at a given k-point,
!! using the Chebyshev filtering method (2021 implementation using xG abstraction layer)
!! for a given spin-polarization, from a fixed hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point.
!! it will also update the matrix elements of the hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2026 ABINIT group (BS)
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

module m_chebfiwf

 use defs_abitypes
 use defs_basis
 use m_abicore
 use m_errors
 use m_fstrings
 use m_time

 use m_chebfi
 use m_chebfi2
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

 public :: chebfiwf2

 CONTAINS  !========================================================================================
!!***

!!****f* m_chebfiwf/chebfiwf2
!! NAME
!! chebfiwf2
!!
!! FUNCTION
!! This routine updates the whole wave functions set at a given k-point,
!! using the Chebfi method (2021 version using xG abstraction layer)
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

subroutine chebfiwf2(cg,dtset,eig,occ,enl_out,gs_hamk,mpi_enreg,&
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
 integer, parameter :: tim_chebfiwf2 = 1750
 integer, parameter :: tim_nonlop = 1753
 integer :: iband,shift,space,blockdim,total_spacedim,ierr
 integer :: me_g0,me_g0_fft
 integer(kind=c_size_t) :: localMem
 type(chebfi_t) :: chebfi
 type(xgBlock_t) :: xgx0,xgeigen,xgocc,xgresidu
 ! arrays
 real(dp) :: tsec(2)
 integer(kind=c_size_t) :: chebfiMem(2)
 real(dp), allocatable :: gvnlxc(:,:),occ_tmp(:)

 ! Parameters for nonlop call in NC
 integer,parameter :: choice=1, paw_opt=0, signs=1
 real(dp) :: gsc_dummy(1,1)
 type(pawcprj_type) :: cprj_dum(gs_hamk%natom,1)

! *********************************************************************

!################ INITIALIZATION  #####################################
!######################################################################

  call timab(tim_chebfiwf2,1,tsec)

!Set module variables
 l_paw = (gs_hamk%usepaw==1)
 l_cpopt=-1;l_sij_opt=0;if (l_paw) l_sij_opt=1
 l_nspinor = nspinor
 l_prtvol = prtvol
 l_mpi_enreg => mpi_enreg
 l_gs_hamk => gs_hamk
 l_paral_kgb = dtset%paral_kgb
 l_block_sliced = dtset%invovl_blksliced

!Variables
 blockdim=l_mpi_enreg%nproc_band*l_mpi_enreg%bandpp
 !for debug
 l_useria=dtset%useria

!Depends on istwfk
 if ( gs_hamk%istwf_k > 1 ) then ! Real only
   ! SPACE_CR mean that we have complex numbers but no re*im terms only re*re
   ! and im*im so that a vector of complex is consider as a long vector of real
   ! therefore the number of data is (2*npw*nspinor)*nband
   ! This space is completely equivalent to SPACE_R but will correctly set and
   ! get the array data into the xgBlock
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

!Memory info
 if ( prtvol >= 3 ) then
   if (l_mpi_enreg%paral_kgb == 1) then
     total_spacedim = l_icplx*npw*nspinor
     call xmpi_sum(total_spacedim,l_mpi_enreg%comm_bandspinorfft,ierr)
   else
     total_spacedim = 0
   end if
   chebfiMem = chebfi_memInfo(nband,l_icplx*npw*nspinor,space,l_mpi_enreg%paral_kgb, &
&                             total_spacedim,l_mpi_enreg%bandpp) !blockdim
   localMem = (int(2,c_size_t)*npw*nspinor*nband+3*nband)*kind(1.d0) !blockdim
   write(std_out,'(1x,A,F10.6,1x,A)') "Each MPI process calling chebfi should need around ", &
   (localMem+sum(chebfiMem))/1e9,"GB of peak memory as follows :"
   write(std_out,'(4x,A,F10.6,1x,A)') "Permanent memory in chebfiwf : ",real(localMem)/1e9,"GB"
   write(std_out,'(4x,A,F10.6,1x,A)') "Permanent memory in m_chebfi : ",real(chebfiMem(1))/1e9,"GB"
   write(std_out,'(4x,A,F10.6,1x,A)') "Temporary memory in m_chebfi : ",real(chebfiMem(2))/1e9,"GB"
 end if

#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET ENTER DATA MAP(to:cg,eig,resid,occ) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif

 call xgBlock_map(xgx0,cg,space,npw*nspinor,nband,comm=l_mpi_enreg%comm_bandspinorfft,me_g0=me_g0,&
   & gpu_option=dtset%gpu_option)

 call xgBlock_map_1d(xgeigen,eig,SPACE_R,nband,gpu_option=dtset%gpu_option)

 call xgBlock_map_1d(xgresidu,resid,SPACE_R,nband,gpu_option=dtset%gpu_option)

 ! Occupancies in chebyshev are used for convergence criteria only
 if (dtset%nbdbuf==-101.and.nspinor==1.and.dtset%nsppol==1) then
   ABI_MALLOC(occ_tmp,(nband))
   occ_tmp(:) = half*occ(:)
   call xgBlock_map_1d(xgocc,occ_tmp,SPACE_R,nband,gpu_option=dtset%gpu_option)
 else
   call xgBlock_map_1d(xgocc,occ,SPACE_R,nband,gpu_option=dtset%gpu_option)
 end if

 call timab(tim_chebfiwf2,2,tsec)

 ABI_NVTX_START_RANGE(NVTX_CHEBFI2_INIT)
 call chebfi_init(chebfi,nband,npw*nspinor,dtset%tolwfr_diago,dtset%ecut, &
&                 dtset%paral_kgb,l_mpi_enreg%bandpp, &
&                 dtset%nline, dtset%nbdbuf, space,1, &
&                 l_mpi_enreg%comm_bandspinorfft,me_g0,me_g0_fft,l_paw,&
&                 l_mpi_enreg%comm_spinorfft,l_mpi_enreg%comm_band,&
&                 dtset%chebfi_oracle,dtset%oracle_factor,dtset%oracle_min_occ,&
&                 l_gs_hamk%gpu_option,gpu_kokkos_nthrd=dtset%gpu_kokkos_nthrd,&
&                 gpu_thread_limit=dtset%gpu_thread_limit)
 ABI_NVTX_END_RANGE()

!################    RUUUUUUUN    #####################################
!######################################################################

 call chebfi_run(chebfi,xgx0,getghc_gsc1,getBm1X,xgeigen,xgocc,xgresidu,nspinor)

 if (allocated(occ_tmp)) then
   ABI_FREE(occ_tmp)
 end if

 if ( .not. l_paw ) then
   call timab(tim_nonlop,1,tsec)
#ifdef FC_CRAY
   ABI_MALLOC(gvnlxc,(1,1))
#else
   ABI_MALLOC(gvnlxc,(0,0))
#endif
   !end if

   ABI_NVTX_START_RANGE(NVTX_CHEBFI2_NONLOP)
   !Call nonlop
   if (l_paral_kgb==0) then

     call nonlop(choice,l_cpopt,cprj_dum,enl_out,l_gs_hamk,0,eig,mpi_enreg,nband,1,paw_opt,&
&                signs,gsc_dummy,l_tim_getghc,cg,gvnlxc)

   else
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET UPDATE FROM(cg) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif
     do iband=1,nband/blockdim
       shift = (iband-1)*blockdim*npw*nspinor
       call prep_nonlop(choice,l_cpopt,cprj_dum, &
&        enl_out((iband-1)*blockdim+1:iband*blockdim),l_gs_hamk,0,&
&        eig((iband-1)*blockdim+1:iband*blockdim),blockdim,mpi_enreg,1,paw_opt,signs,&
&        gsc_dummy,l_tim_getghc,cg(:,shift+1:shift+blockdim*npw*nspinor),gvnlxc(:,:),&
&        already_transposed=.false.)
     end do
   end if
   ABI_NVTX_END_RANGE()
   ABI_FREE(gvnlxc)
   call timab(tim_nonlop,2,tsec)
 end if

!Free chebfi
 call chebfi_free(chebfi)

#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET UPDATE FROM(cg,eig,resid) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
 !$OMP TARGET EXIT DATA MAP(delete:cg,eig,resid,occ) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif

 call timab(tim_chebfiwf2,2,tsec)

 DBG_EXIT("COLL")

end subroutine chebfiwf2
!!***

!----------------------------------------------------------------------

!!****f* m_chebfiwf/getghc_gsc1
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
 real(dp)          :: gvnlxc(1,1)

! *********************************************************************

 ABI_NVTX_START_RANGE(NVTX_GETGHC)

 call xgBlock_getSize(X,spacedim,blockdim)
 call xgBlock_check(X,AX)
 call xgBlock_check(X,BX)

 call xgBlock_reverseMap(X,cg,rows=1,cols=spacedim*blockdim)
 call xgBlock_reverseMap(AX,ghc,rows=1,cols=spacedim*blockdim)
 call xgBlock_reverseMap(BX,gsc,rows=1,cols=spacedim*blockdim)

 call multithreaded_getghc(l_cpopt,cg,cprj_dum,ghc,gsc,&
   l_gs_hamk,gvnlxc,eval,l_mpi_enreg,blockdim,l_prtvol,l_sij_opt,l_tim_getghc,0)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
 call gpu_device_synchronize()
#endif

 if ( .not. l_paw ) call xgBlock_copy(X,BX)

 ABI_NVTX_END_RANGE()

end subroutine getghc_gsc1
!!***

!----------------------------------------------------------------------

!!****f* m_chebfiwf/getBm1X
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

end module m_chebfiwf
!!***
