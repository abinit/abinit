!!****f* ABINIT/m_slicewf
!! NAME
!! m_slicewf
!!
!! FUNCTION
!! This module contains a routine updating the whole wave functions at a given k-point,
!! using the Chebyshev filtering method (2021 implementation using xG abstraction layer)
!! for a given spin-polarization, from a fixed hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point.
!! it will also update the matrix elements of the hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2025 ABINIT group (BS, I. Lygatsika)
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

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
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
 integer, parameter :: tim_permute = 2164 
 integer :: tim_sliceX_copy ! slice timers
 integer :: iband,shift,space,blockdim,total_spacedim,ierr
 integer :: islice,nslice,i1,i2,j1,j2,npband
 integer :: nband_slice,nband_merge,merge_option
 integer :: spacedim,spacecom,gpu_option
 integer :: me_g0,me_g0_fft
 integer :: my_rank, my_color ! for MPI
 integer(kind=c_size_t) :: localMem
 type(sliceAll_t) :: sliceAll
 type(slice_t) :: slice
 type(xgBlock_t) :: xgx0,xgeigen,xgresidu
 ! arrays
 real(dp) :: tsec(2)
 integer(kind=c_size_t) :: sliceMem(2)
 real(dp), allocatable :: l_gvnlxc(:,:)
 integer, target, allocatable :: pband(:)
 integer, target, allocatable :: idx(:,:)
 integer, target, allocatable :: idx_ovlp(:,:)
 integer, target, allocatable :: idx_merge(:,:)
 integer, target, allocatable :: ndeg(:)
 integer, target, allocatable :: npbandSlice(:)
 integer, pointer :: pband_ptr(:) => NULL()
 integer, pointer :: idx_ptr(:,:) => NULL()
 integer, pointer :: idx_ovlp_ptr(:,:) => NULL()
 integer, pointer :: idx_merge_ptr(:,:) => NULL()
 integer, pointer :: ndeg_ptr(:) => NULL()
 integer, pointer :: npbandSlice_ptr(:) => NULL()
 real(dp), target, allocatable :: sbound(:,:)
 real(dp), pointer :: sbound_ptr(:,:) => NULL()
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

!Memory info TODO for slice
! if ( prtvol >= 3 ) then
!   if (l_mpi_enreg%paral_kgb == 1) then
!     total_spacedim = l_icplx*l_npw*l_nspinor
!     call xmpi_sum(total_spacedim,l_mpi_enreg%comm_bandspinorfft,ierr)
!   else
!     total_spacedim = 0
!   end if
!   sliceMem = chebfi_memInfo(nband,l_icplx*l_npw*l_nspinor,space,l_mpi_enreg%paral_kgb, &
!&                             total_spacedim,l_mpi_enreg%bandpp) !blockdim
!   localMem = (int(2,c_size_t)*l_npw*l_nspinor*nband+3*nband)*kind(1.d0) !blockdim
!   write(std_out,'(1x,A,F10.6,1x,A)') "Each MPI process calling chebfi should need around ", &
!   (localMem+sum(sliceMem))/1e9,"GB of peak memory as follows :"
!   write(std_out,'(4x,A,F10.6,1x,A)') "Permanent memory in slicewf : ",real(localMem)/1e9,"GB"
!   write(std_out,'(4x,A,F10.6,1x,A)') "Permanent memory in m_slice : ",real(sliceMem(1))/1e9,"GB"
!   write(std_out,'(4x,A,F10.6,1x,A)') "Temporary memory in m_slice : ",real(sliceMem(2))/1e9,"GB"
! end if

 ! Two senarios: (choose one of the two)
 ! 1)Send cg (synchronous mem) to GPU and point xgBlock to cg 
!#ifdef HAVE_OPENMP_OFFLOAD
! !$OMP TARGET ENTER DATA MAP(to:cg) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
!#endif
 !call xgBlock_map(xgx0,cg,space,spacedim,nband,comm=spacecom,me_g0=me_g0,gpu_option=gpu_option)
 ! 2) Point xgBlock to cg, on CPU always
 call xgBlock_map(xgx0,cg,space,spacedim,nband,comm=spacecom,me_g0=me_g0,gpu_option=ABI_GPU_DISABLED)

 !ierr = slice_unitTest(xgx0)

!################### ASYNCHRONOUS MEMORY BUFFER #######################
! cg cannot be modified (read/write) in parallel by multiple MPI
! processes due to converged vectors overlapping between neighbor slices.
! The purpose of the buffer is to safely read/write overlap data without 
! blocking communications.
!######################################################################

 npband = l_mpi_enreg%nproc_band
 nslice = dtset%nslice

 ! Memory allocations of size depending on fixed nslice
 ABI_MALLOC(pband, (nband)); pband_ptr => pband                     ! band permutation
 ABI_MALLOC(idx, (nslice,2)); idx_ptr => idx                        ! idx in overlapping mem
 ABI_MALLOC(idx_ovlp, (nslice,2)); idx_ovlp_ptr => idx_ovlp         ! idx in overlap-free mem
 ABI_MALLOC(idx_merge, (nslice,2)); idx_merge_ptr => idx_merge      ! converged idx in overlap-free mem
 ABI_MALLOC(ndeg, (nslice)); ndeg_ptr => ndeg                       ! filter degree per slice
 ABI_MALLOC(sbound, (nslice,4)); sbound_ptr => sbound               ! eigenvalue bounds per slice
 ABI_MALLOC(npbandSlice, (nslice)); npbandSlice_ptr => npbandSlice  ! number of mpi processes per slice

 ! Initialize values
 pband(:) = (/(iband, iband=1,nband)/)
 if (dtset%paral_slice == 0) then
    npbandSlice(:) = (/(npband, islice=1,nslice)/)                 
 end if

 ! ************ Initialize spectrum slicing workspace: Rayleigh values and residuals
 write(std_out,'(a)') '1) Init sliceAll'
 ABI_NVTX_START_RANGE(NVTX_SLICEALL_INIT)
 call sliceAll_init(sliceAll,nband,spacedim,dtset%tolwfr_diago,dtset%ecut,&
&                   dtset%paral_kgb,l_mpi_enreg%bandpp,dtset%mdeg_filter,space,1,&
&                   l_mpi_enreg%comm_bandspinorfft,me_g0,me_g0_fft,l_paw,&
&                   l_mpi_enreg%comm_spinorfft,l_mpi_enreg%comm_band,&
&                   nslice,npband,dtset%tolfilter,dtset%balfilter,&
&                   dtset%nbdbuf,0,dtset%oracle_factor,dtset%oracle_min_occ,& ! oracle=0
&                   l_gs_hamk%gpu_option,gpu_kokkos_nthrd=dtset%gpu_kokkos_nthrd)
 ABI_NVTX_END_RANGE()

 ! ************ Compute Density Of States (DOS)
 write(std_out,'(a)') '2) Compute DOS'
 ABI_NVTX_START_RANGE(NVTX_SLICEALL_DOS)
 call sliceAll_dos(sliceAll,xgx0,getghc_gsc1,nspinor)
 ABI_NVTX_END_RANGE()

 ! ************ Partition spectrum into slices
 write(std_out,'(a,i0,a)') '3) Partition spectrum into ',nslice,' slices'
 ABI_NVTX_START_RANGE(NVTX_SLICEALL_SPLIT)
 call sliceAll_split(sliceAll,idx_ptr,ndeg_ptr,sbound_ptr,pband_ptr,npbandSlice_ptr)
 ABI_NVTX_END_RANGE()

 ! ************ Distribute MPI procs across slices
 if (dtset%paral_slice == 1) then
    call slice_findOptimalNumMpiProcs(npband,nslice,idx_ptr,ndeg_ptr,npbandSlice_ptr)
    write(std_out,*) ' ==== paral_slice option detected'
    write(std_out,*) 'optimal num mpi procs:'
    write(std_out,*) npbandSlice(:)
    write(std_out,*) ' '
    
    ! reset for the moment because not implemented
    npbandSlice(:) = (/(npband, islice=1,nslice)/)
 else if (dtset%paral_slice == 2) then
    
    ! Equal distribution
    npbandSlice(:) = (/(npband/nslice, islice=1,nslice)/)
    
    ! TODO color mpi rank
    if (xmpi_comm_size(comm_cols)>1) then
        my_rank = xmpi_comm_rank(comm_cols)
        my_color = my_rank / npband_slice
   end if
                   
 end if

 ! data transfer? H2D D2H just to accelerate one AX
 ! TODO: deactivate all that, perform entirely on CPU

 ! cg not needed on GPU for slice_run, copy (update D2H) then delete from GPU to free space
!#ifdef HAVE_OPENMP_OFFLOAD
! !$OMP TARGET UPDATE FROM(cg) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
! !$OMP TARGET EXIT DATA MAP(delete:cg) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
!#endif

 ! Permute eigenvectors in RR quotient-increasing order (implemented on CPU only)
 ABI_NVTX_START_RANGE(NVTX_SLICEALL_PERMUTE_COLS)
 call xgBlock_permuteCols(xgx0,spacedim,nband,pband_ptr)
 ABI_NVTX_END_RANGE()

 !ierr = slice_unitTest(xgx0)

 ! ************ Allocate parallel-safe memory buffer on CPU.
 ABI_NVTX_START_RANGE(NVTX_SLICEALL_INIT_ASYNC_BUFFER)
 call sliceAll_initOverlapFree(sliceAll)
 ABI_NVTX_END_RANGE()

 !ierr = slice_unitTest(sliceALl%xgx0_ovlp)
 
 ! Copy range of cg to range of memory buffer.
 ! Assumes that cg is distributed on MPI rows (so each MPI has all bands).
 ! FIXME if cg is distributed on MPI columns, the problem is that
 ! we have to communicate between MPI columns.
 write(std_out,'(a)') '4) Copy to buffer (parallel safe memory space)'

 ! Define index range in buffer
 idx_ovlp(1:nslice,1) = (/ (1 + idx(islice,2) - idx(islice,1) + 1, islice=1,nslice) /)
 idx_ovlp(1:nslice,2) = idx_ovlp(1:nslice,1) + 1

 ! TODO use pointers.
 ! This command is executed on every MPI process, containing its own rows of xgx0 and all cols.
 ! Since no communication takes place, we can read and write to the MPI part independently of others.
 call xgBlock_setBlock(xgx0,spacedim,nband_slice,fcol=i1)

 if (use_subcomm_) then
    ! each mpi has the bands it has to copy. Can do in parallel
 else
    ! all mpis have all bands
 end if

 j1 = 1
 ! TODO: MPI distribute over columns
 do islice=1,nslice
    ABI_NVTX_START_RANGE(NVTX_SLICE_COPY)
    i1 = idx(islice,1)        ! start read range from cg
    i2 = idx(islice,2)        ! end
    nband_slice = i2 - i1 + 1
    j2 = j1 + nband_slice - 1
    idx_ovlp(islice,1) = j1   ! start write range to buffer
    idx_ovlp(islice,2) = j2   ! end
    call slice_blockCopy(xgx0,sliceAll%xgx0_ovlp,i1,j1,i2,j2)
    j1 = j2 + 1
    ABI_NVTX_END_RANGE()
 end do
 write(std_out,*) 'end 4)'
 
!################    RUUUUUUUN    #####################################
!######################################################################

 ! ************** Diagonalize each slice (sequential or parallel TODO)
 ! Internally treat the case of parallel slices
 ! TODO add member variable paral_slice to sliceAll
 call sliceAll_run(sliceAll,dtset%paral_slice,use_subcomm_)

 ! ================
 ! TODO IL 31/01/2025
 ! * diagnostic de convergence en utilisant residual ratio (r_i/r_i^n > ramp)
 ! * Plot residuals in a slice and see where they are large?
 ! ================

!################ OVERLAP-FREE MEMORY -> CG ###########################
!######################################################################
 
 ! ************** Detect converged eigenvalues from each slice
 ABI_NVTX_START_RANGE(NVTX_SLICEALL_MERGE)
 merge_option = 0 ! FIXME hard-coded 
 call sliceAll_merge(sliceAll,idx_ovlp_ptr,idx_merge_ptr,merge_option)
 ABI_NVTX_END_RANGE()

 ! Point xgeigen and xgresidy to CPU objects
 call xgBlock_map_1d(xgeigen,eig,SPACE_R,nband,gpu_option=ABI_GPU_DISABLED)
 call xgBlock_map_1d(xgresidu,resid,SPACE_R,nband,gpu_option=ABI_GPU_DISABLED)

 ! Copy from PART OF overlap-free to PART OF overlapping mem
 ! using final buffer written only once when all slices has converged
 ! -----------------------------------------> race condition
 write(std_out,'(a)') '6) Copy from safe buffer (overlap-free memory space)'
 call xgBlock_reshape(xgeigen, (/1,nband/))
 call xgBlock_reshape(xgresidu, (/1,nband/))
 j1 = 1                      ! start copy to overlapping mem (cg,eig,resid)
 do islice=1,nslice
    ABI_NVTX_START_RANGE(NVTX_SLICE_COPY)
    i1 = idx_merge(islice,1) ! start read from overlap-free mem
    i2 = idx_merge(islice,2)
    nband_merge = i2 - i1 + 1
    j2 = j1 + nband_merge - 1
    write(std_out,*) 'copy to xgx0'
    call slice_blockCopy(sliceAll%xgx0_ovlp,xgx0,i1,j1,i2,j2)
    write(std_out,*) 'copy to xgeigen'
    call slice_blockCopy(sliceAll%xgeigen_ovlp,xgeigen,i1,j1,i2,j2)
    write(std_out,*) 'copy to xgresidu'
    call slice_blockCopy(sliceAll%xgresidu_ovlp,xgresidu,i1,j1,i2,j2)
    j1 = j2 + 1
    ABI_NVTX_END_RANGE()
 end do
 call xgBlock_reshape(xgeigen, (/nband,1/))
 call xgBlock_reshape(xgresidu, (/nband,1/))

 ! Print final eigenvalues after merge
 !write(std_out,*) 'final eigenvalues='
 !call xgBlock_print(xgeigen,std_out)

 ! Free slice parameters
 if (allocated(pband)) ABI_FREE(pband)
 if (allocated(idx)) ABI_FREE(idx)
 if (allocated(idx_ovlp)) ABI_FREE(idx_ovlp)
 if (allocated(idx_merge)) ABI_FREE(idx_merge)
 if (allocated(ndeg)) ABI_FREE(ndeg)
 if (allocated(sbound)) ABI_FREE(sbound)
 if (allocated(npbandSlice)) ABI_FREE(npbandSlice)

 ! Free spectrum slicing workspace
 write(std_out,'(a)') '7) Free sliceAll and overlap-free workspace'
 ABI_NVTX_START_RANGE(NVTX_SLICEALL_FREE_ASYNC_BUFFER)
 call sliceAll_free(sliceAll)
 ABI_NVTX_END_RANGE()
 
 ! Send cg,eig,resid to GPU (needed for nonlop)
 if ( .not. l_paw .and. l_paral_kgb==1 ) then
#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET ENTER DATA MAP(to:cg,eig,resid) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif
 end if

! =====================================================================================
! spectrum slicing finished

 if ( .not. l_paw ) then
   call timab(tim_nonlop,1,tsec)
#ifdef FC_CRAY
   ABI_MALLOC(l_gvnlxc,(1,1))
#else
   ABI_MALLOC(l_gvnlxc,(0,0))
#endif
   !end if

   ABI_NVTX_START_RANGE(NVTX_SLICES_NONLOP)
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

#ifdef HAVE_OPENMP_OFFLOAD
 !$OMP TARGET UPDATE FROM(cg,eig,resid) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
 !$OMP TARGET EXIT DATA MAP(delete:cg,eig,resid) IF(gs_hamk%gpu_option==ABI_GPU_OPENMP)
#endif

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

 call xgBlock_reverseMap(X,cg,rows=1,cols=spacedim*blockdim)
 call xgBlock_reverseMap(AX,ghc,rows=1,cols=spacedim*blockdim)
 call xgBlock_reverseMap(BX,gsc,rows=1,cols=spacedim*blockdim)

 call multithreaded_getghc(l_cpopt,cg,cprj_dum,ghc,gsc,&
   l_gs_hamk,l_gvnlxc,eval,l_mpi_enreg,blockdim,l_prtvol,l_sij_opt,l_tim_getghc,0)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
 call gpu_device_synchronize()
#endif

 if ( .not. l_paw ) call xgBlock_copy(X,BX)

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

   ! IML 17/10 workaround for spectrum slicing when bandpp (abinit) =/= bandpp_slice
   ! 18/10 delete this if the modifs in 66_wfs/m_prep_kgb works in lines
   ! 703, 707, 741 for bandpp = blocksize, by default =/= mpi_enreg%bandpp
   !l_mpi_enreg%bandpp = blockdim

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
