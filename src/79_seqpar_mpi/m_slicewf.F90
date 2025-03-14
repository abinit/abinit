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

 !FIXME Keep those in these modules or moves them together ?
 use m_invovl,             only : invovl_ompgpu_static_mem,invovl_ompgpu_work_mem
 use m_gemm_nonlop_ompgpu, only : gemm_nonlop_ompgpu_static_mem
 use m_getghc_ompgpu,      only : getghc_ompgpu_work_mem

#if defined(HAVE_GPU)
 use m_gpu_toolbox
#endif

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

 public :: slicewf_blocksize
 public :: slicewf

 CONTAINS  !========================================================================================
!!***

subroutine slicewf_blocksize(gs_hamk,ndat,npw,nband,nspinor,paral_kgb,gpu_option,nblk_gemm_nonlop)
   implicit none

   integer,intent(in) :: ndat,npw,nband,nspinor,paral_kgb,gpu_option
   type(gs_hamiltonian_type),intent(in) :: gs_hamk
   integer, intent(out)  :: nblk_gemm_nonlop

   integer(kind=c_size_t) :: nonlop_smem,invovl_smem,getghc_wmem,invovl_wmem
   integer(kind=c_size_t) :: sum_mem,sum_bandpp_mem,sum_other_mem,free_mem
   integer  :: icplx,space,i,ndat_try,rank,nprocs
   real(dp) :: localMem,sliceMem(2)

! *********************************************************************

   free_mem=256*1e9 ! Dummy value
#ifdef HAVE_GPU
   if(gpu_option /= ABI_GPU_DISABLED) then
     call gpu_get_max_mem(free_mem)
     free_mem = 0.95 * free_mem ! Cutting 5% out to be safe
   end if
#else
   ABI_UNUSED(gpu_option)
#endif
   rank = xmpi_comm_rank(xmpi_world); nprocs = xmpi_comm_size(xmpi_world)
   if ( gs_hamk%istwf_k == 2 ) then ! Real only
     space = SPACE_CR
     icplx = 2
   else ! complex
     space = SPACE_C
     icplx = 1
   end if

   call xmpi_barrier(xmpi_world)
   ndat_try=ndat
   nonlop_smem = gemm_nonlop_ompgpu_static_mem(gs_hamk%npw_fft_k, gs_hamk%indlmn, gs_hamk%nattyp, gs_hamk%ntypat, 1)
   invovl_smem = invovl_ompgpu_static_mem(gs_hamk)
   getghc_wmem = getghc_ompgpu_work_mem(gs_hamk, ndat_try)
   invovl_wmem = invovl_ompgpu_work_mem(gs_hamk, ndat_try)

   sliceMem = chebfi_memInfo(nband,icplx*npw*nspinor,space,paral_kgb,icplx*npw*nspinor,ndat)
   localMem  = (npw+2*npw*nspinor+2*nband)*kind(1.d0) !blockdim

   sum_mem          = nonlop_smem+invovl_smem+getghc_wmem+invovl_wmem+sliceMem(1)+sliceMem(2)+localMem
   sum_bandpp_mem   = getghc_wmem+invovl_wmem
   sum_other_mem    = nonlop_smem+invovl_smem+sliceMem(1)+sliceMem(2)+localMem

   nblk_gemm_nonlop=1

   ! No blocking needed, all good !
   if(sum_mem < free_mem) return

   write(std_out,*) "Setting block size..."
   ! How the number of blocks is decided:
   ! We try to divide bandpp with dividers from 1 to 20
   ! If we fail, that means test case is too fat for given hardware, and that's it
   ! This looks stupid but we don't actually expect to process CHEBFI with 20 blocks.
   do i=1,20

     ! Gemm nonlop static memory requirement is higher, split here
     nblk_gemm_nonlop = nblk_gemm_nonlop + 1
     if(modulo(nprocs,nblk_gemm_nonlop)/=0) cycle

     nonlop_smem = gemm_nonlop_ompgpu_static_mem(gs_hamk%npw_fft_k,gs_hamk%indlmn,gs_hamk%nattyp,gs_hamk%ntypat,nblk_gemm_nonlop)

     ! Bandpp~ndat sized buffer memory requirements are higher, split there
     sum_mem          = nonlop_smem+invovl_smem+getghc_wmem+invovl_wmem+sliceMem(1)+sliceMem(2)+localMem
     sum_bandpp_mem   = getghc_wmem+invovl_wmem
     sum_other_mem    = nonlop_smem+invovl_smem+sliceMem(1)+sliceMem(2)+localMem

     write(std_out,'(A,F10.3,1x,A)') "Free mem                                   : ", real(free_mem)/(1024*1024), "MiB"
     write(std_out,*) "Memory requirements of slicewf per MPI task (OpenMP GPU)"
     write(std_out,*) "---------------------------------------------------------"
     write(std_out,*) "Static buffers, computed once and permanently on card :"
     write(std_out,'(A,F10.3,1x,A)') "   gemm_nonlop_ompgpu (projectors)       : ",  real(nonlop_smem,dp)/(1024*1024), "MiB"
     write(std_out,'(A,F10.3,1x,A)') "   invovl_ompgpu (mkinvovl)              : ",  real(invovl_smem,dp)/(1024*1024), "MiB"
     write(std_out,'(A,F10.3,1x,A)') "   slice                               : ",          sliceMem(1)/(1024*1024), "MiB"
     write(std_out,*) "Work buffers, temporary, bandpp sized  :"
     write(std_out,'(A,F10.3,1x,A)') "   getghc (inc. fourwf+gemm_nonlop)      : ",  real(getghc_wmem,dp)/(1024*1024), "MiB"
     write(std_out,'(A,F10.3,1x,A)') "   invovl                                : ",  real(invovl_wmem,dp)/(1024*1024), "MiB"
     write(std_out,'(A,F10.3,1x,A)') "   slice (RR buffers)                  : ",          sliceMem(2)/(1024*1024), "MiB"
     write(std_out,'(A,F10.3,1x,A)') "   slicewf (cg,resid,eig)               : ",              localMem/(1024*1024), "MiB"
     write(std_out,*) "---------------------------------------------------------"
     write(std_out,'(A,F10.3,1x,A)') "Sum                                      : ", real(sum_mem)/(1024*1024), "MiB"
     flush(std_out)

     if(sum_mem < free_mem) exit
   end do
   if(sum_mem > free_mem) then
     ABI_ERROR("It seems the test case you're trying to run is too big to run with given hardware resources !")
   end if

 end subroutine slicewf_blocksize
!!***

!----------------------------------------------------------------------

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
 integer :: iband,shift,space,blockdim,total_spacedim,ierr
 integer :: islice,nslice,i1,i2,j1,j2,npband
 integer :: nband_slice,nband_merge,merge_option
 integer :: spacedim,spacecom,gpu_option
 integer :: me_g0,me_g0_fft
 integer :: tim_sliceX_copy ! slice timers
 real(dp) :: localmem
 type(sliceAll_t) :: sliceAll
 type(slice_t) :: slice
 type(xgBlock_t) :: xgx0,xgeigen,xgresidu
 ! arrays
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
 real(dp) :: tsec(2),chebfiMem(2)
 real(dp), target, allocatable :: sbound(:,:)
 real(dp), pointer :: sbound_ptr(:,:) => NULL()
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

!Memory info TODO
! if ( prtvol >= 3 ) then
!   if (l_mpi_enreg%paral_kgb == 1) then
!     total_spacedim = l_icplx*l_npw*l_nspinor
!     call xmpi_sum(total_spacedim,l_mpi_enreg%comm_bandspinorfft,ierr)
!   else
!     total_spacedim = 0
!   end if
!   chebfiMem = chebfi_memInfo(nband,l_icplx*l_npw*l_nspinor,space,l_mpi_enreg%paral_kgb, &
!&                             total_spacedim,l_mpi_enreg%bandpp) !blockdim
!   localMem = (l_npw+2*l_npw*l_nspinor+2*nband)*kind(1.d0) !blockdim
!   write(std_out,'(1x,A,F10.6,1x,A)') "Each MPI process calling chebfi should need around ", &
!   (localMem+sum(chebfiMem))/1e9,"GB of peak memory as follows :"
!   write(std_out,'(4x,A,F10.6,1x,A)') "Permanent memory in chebfiwf : ",(localMem)/1e9,"GB"
!   write(std_out,'(4x,A,F10.6,1x,A)') "Permanent memory in m_chebfi : ",(chebfiMem(1))/1e9,"GB"
!   write(std_out,'(4x,A,F10.6,1x,A)') "Temporary memory in m_chebfi : ",(chebfiMem(2))/1e9,"GB"
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

 ! ************ Allocate asynchronous memory buffer (always on CPU)
 ABI_NVTX_START_RANGE(NVTX_SLICEALL_INIT_ASYNC_BUFFER)
 call sliceAll_initOverlapFree(sliceAll)
 ABI_NVTX_END_RANGE()

 !ierr = slice_unitTest(sliceALl%xgx0_ovlp)
 
 ! Copy from PART OF cg to PART OF overlap-free mem 
 ! this buffer allows asyncrhonous read/write between slices
 ! ---------------------------------------------> race condition
 write(std_out,'(a)') '4) Copy to safe buffer (overlap-free memory space)'
 j1 = 1
 ! TODO: MPI distribute over columns
 do islice=1,nslice
    ABI_NVTX_START_RANGE(NVTX_SLICE_COPY)
    i1 = idx(islice,1)        ! start read from cg
    i2 = idx(islice,2)        ! end
    nband_slice = i2 - i1 + 1
    j2 = j1 + nband_slice - 1
    idx_ovlp(islice,1) = j1   ! start copy to overlap-free mem
    idx_ovlp(islice,2) = j2   ! end
    call slice_blockCopy(xgx0,sliceAll%xgx0_ovlp,i1,j1,i2,j2)
    j1 = j2 + 1
    ABI_NVTX_END_RANGE()
 end do
 write(std_out,*) 'end 4)'
 
!################    RUUUUUUUN    #####################################
!######################################################################

 ! ************** Diagonalize each slice (sequential or parallel TODO)
 do islice=1,nslice
    
    write(std_out,'(a,i0)') '5) Diago slice ',islice
    
    ! Allocate slice memory, on GPU
    ABI_NVTX_START_RANGE(NVTX_SLICE_INIT)
    write(std_out,*) 'TRACE slice_init'
    call slice_init(sliceAll,slice,islice)
    ABI_NVTX_END_RANGE()
    
    !write(std_out,*) 'slice%xgx0'
    !ierr = slice_unitTest(slice%xgx0)
    
    ! Asynchronous copy: read from X_safe write to X_slice
    ABI_NVTX_START_RANGE(NVTX_SLICE_COPY)
    j1 = idx_ovlp(islice,1)
    j2 = idx_ovlp(islice,2)
    nband_slice = j2 - j1 + 1
    write(std_out,*) 'TRACE slice_blockCopy'
    call xgBlock_copy_from_gpu(slice%xgx0)
    call slice_blockCopy(sliceAll%xgx0_ovlp,slice%xgx0,j1,1,j2,nband_slice) 
    call xgBlock_copy_to_gpu(slice%xgx0)
    ABI_NVTX_END_RANGE()

    !write(std_out,*) 'sliceAll%xgx0_ovlp'
    !ierr = slice_unitTest(sliceAll%xgx0_ovlp)
    !write(std_out,*) 'slice%xgx0'
    !ierr = slice_unitTest(slice%xgx0)
    
    ! Run
    ABI_NVTX_START_RANGE(NVTX_SLICE_RUN)
    write(std_out,*) 'TRACE slice_run'
    call slice_run(slice,getghc_gsc1,getBm1X,nspinor) 
    ABI_NVTX_END_RANGE() 

    !write(std_out,*) 'slice%xgx0'
    !ierr = slice_unitTest(slice%xgx0)
    
    ! Asynchronous copy: read from X_slice write to X_safe
    ABI_NVTX_START_RANGE(NVTX_SLICE_COPY)
    write(std_out,*) 'TRACE slice_blockCopy'
    call xgBlock_reshape(slice%xgeigen, (/1,nband_slice/))
    call xgBlock_reshape(slice%xgresidu, (/1,nband_slice/))
    call xgBlock_copy_from_gpu(slice%xgx0)
    call xgBlock_copy_from_gpu(slice%xgeigen)
    call xgBlock_copy_from_gpu(slice%xgresidu)
    call slice_blockCopy(slice%xgx0,sliceAll%xgx0_ovlp,1,j1,nband_slice,j2)
    call slice_blockCopy(slice%xgeigen,sliceAll%xgeigen_ovlp,1,j1,nband_slice,j2)
    call slice_blockCopy(slice%xgresidu,sliceAll%xgresidu_ovlp,1,j1,nband_slice,j2)
    ABI_NVTX_END_RANGE()
    
    !write(std_out,*) 'sliceAll%xgx0_ovlp'
    !ierr = slice_unitTest(sliceAll%xgx0_ovlp)
    !write(std_out,*) 'slice%xgx0'
    !ierr = slice_unitTest(slice%xgx0)

    ! Clean slice memory
    ABI_NVTX_START_RANGE(NVTX_SLICE_FREE)
    write(std_out,*) 'TRACE slice_free'
    call slice_free(slice)
    ABI_NVTX_END_RANGE()

 end do

 ! ================
 ! TODO semaine 27-31 janvier
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

 !write(std_out,*) 'TRACE:getBm1X spacedim          =', spacedim
 !write(std_out,*) 'TRACE:getBm1X blockdim          =', blockdim
 !write(std_out,*) 'TRACE:getBm1X l_mpi_enreg%bandpp=', l_mpi_enreg%bandpp
 
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
