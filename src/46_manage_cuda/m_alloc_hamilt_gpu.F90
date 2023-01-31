!!****m* ABINIT/m_alloc_hamilt_gpu
!! NAME
!!  m_alloc_hamilt_gpu
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2022 ABINIT group (MT, FDahm)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_alloc_hamilt_gpu

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_dtset
#if defined HAVE_GPU_CUDA
 use m_gpu_toolbox
#endif

#ifdef HAVE_FC_ISO_C_BINDING
 use, intrinsic :: iso_c_binding
#endif

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes, only : MPI_type

 implicit none

 private
!!***

 public :: alloc_hamilt_gpu
 public :: dealloc_hamilt_gpu
!!***

 private :: alloc_nonlop_kokkos
 private :: dealloc_nonlop_kokkos

 !! data type to store pointers to data used on GPU, mostly in gemm nonlop
 !! opernla/b/c
 type, public :: gemm_nonlop_kokkos_type

   logical     :: allocated
   type(c_ptr) ::     projections_gpu
   type(c_ptr) ::   s_projections_gpu
   type(c_ptr) :: vnl_projections_gpu

   type(c_ptr) ::   vectin_gpu
   type(c_ptr) ::  vectout_gpu
   type(c_ptr) :: svectout_gpu

 end type gemm_nonlop_kokkos_type

 type(gemm_nonlop_kokkos_type), save, public, target :: gemm_nonlop_kokkos

contains
!!***

!!****f* ABINIT/alloc_hamilt_gpu
!! NAME
!! alloc_hamilt_gpu
!!
!! FUNCTION
!! allocate several memory pieces on a GPU device for the application
!! of Hamiltonian using a GPU
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  mpi_enreg=information about MPI parallelization
!!  nattyp(ntypat)= # atoms of each type.
!!  option=0: allocate data for local operator (FFT)
!!         1: allocate data for nonlocal operator
!!         2: allocate both
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  use_gpu_cuda= 0 or 1 to know if we use cuda
!!
!! OUTPUT
!!  (no output - only allocation on GPU)
!!
!! SOURCE

subroutine alloc_hamilt_gpu(atindx1,dtset,gprimd,mpi_enreg,nattyp,npwarr,option,psps,use_gpu_cuda)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option,use_gpu_cuda
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(in) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx1(dtset%natom),nattyp(dtset%ntypat),npwarr(dtset%nkpt)
 real(dp),intent(in) :: gprimd(3,3)

!Local variables-------------------------------
!scalars
#if defined HAVE_GPU_CUDA
 integer :: dimekb1_max,dimekb2_max,dimffnl_max,ierr,ikpt,npw_max_loc,npw_max_nonloc
 integer :: cplex
 integer ::npwarr_tmp(dtset%nkpt)

 integer(kind=c_int32_t) :: proj_dim(3)
#endif

! *************************************************************************

 if (use_gpu_cuda==0) return

#if defined(HAVE_GPU_CUDA) && defined(HAVE_FC_ISO_C_BINDING)
!=== Local Hamiltonian ===
 if (option==0.or.option==2) then
!  Compute max of total planes waves
   npw_max_loc=0
   if(mpi_enreg%paral_kgb==1) then
     npwarr_tmp=npwarr
     call xmpi_sum(npwarr_tmp,mpi_enreg%comm_bandfft,ierr)
     npw_max_loc =maxval(npwarr_tmp)
   else
     npw_max_loc=dtset%mpw
   end if
   ! Initialize gpu data needed in fourwf
   ! ndat=bandpp when paral_kgb=1
   ! no matter paral_kgb=0 or 1, we gathet all bands into a single call gpu_fourwf
   call alloc_gpu_fourwf(dtset%ngfft,dtset%bandpp,npw_max_loc,npw_max_loc)

 end if
!=== Nonlocal Hamiltonian ===
 if (option==1.or.option==2) then
!  Compute max of total planes waves
   npw_max_nonloc=0
   if(mpi_enreg%paral_kgb==1) then
     npwarr_tmp=npwarr
     call xmpi_sum(npwarr_tmp,mpi_enreg%comm_band,ierr)
     npw_max_nonloc =maxval(npwarr_tmp)
   else
     npw_max_nonloc=dtset%mpw
   end if
!  Initialize all gpu data needed in nonlop
   dimffnl_max=4
!  if (abs(dtset%berryopt) == 5) dimffnl_max=4
   dimekb1_max=psps%dimekb
   if (dtset%usepaw==0) dimekb2_max=psps%ntypat
   if (dtset%usepaw==1) dimekb2_max=dtset%natom

   ! TODO (PK) : disable this allocation when Kokkos is available
   ! to save memory on GPU side
   call alloc_nonlop_gpu(npw_max_nonloc, &
     &                   npw_max_nonloc,&
     &                   dtset%nspinor,&
     &                   dtset%natom,&
     &                   dtset%ntypat,&
     &                   psps%lmnmax,&
     &                   psps%indlmn,&
     &                   nattyp,&
     &                   atindx1,&
     &                   gprimd,&
     &                   dimffnl_max,&
     &                   dimffnl_max,&
     &                   dimekb1_max,&
     &                   dimekb2_max)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS)
   if (dtset%use_gemm_nonlop_gpu == 1) then

     call alloc_nonlop_kokkos(dtset, &
       &                      psps, &
       &                      npw_max_nonloc,&
       &                      npw_max_nonloc,&
       &                      atindx1,&
       &                      nattyp,&
       &                      use_gpu_cuda)

   end if

#endif

 end if ! option 1 or 2

 call xmpi_barrier(mpi_enreg%comm_cell)

#else

 if (.false.) then
   write(std_out,*) atindx1(1),dtset%natom,gprimd(1,1),mpi_enreg%me,nattyp(1),option
 end if

#endif

end subroutine alloc_hamilt_gpu
!!***

!!****f* ABINIT/dealloc_hamilt_gpu
!! NAME
!! dealloc_hamilt_gpu
!!
!! FUNCTION
!! deallocate several memory pieces on a GPU device used for the application
!! of Hamiltonian using a GPU
!!
!! INPUTS
!!  option=0: deallocate data for local operator (FFT)
!!         1: deallocate data for nonlocal operator
!!         2: deallocate both
!!  use_gpu_cuda= 0 or 1 to know if we use cuda
!!
!! OUTPUT
!!  (no output - only deallocation on GPU)
!!
!! SOURCE

subroutine dealloc_hamilt_gpu(option,use_gpu_cuda)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option,use_gpu_cuda
!arrays

!Local variables-------------------------------

! *************************************************************************

 if (use_gpu_cuda==0) return

#if defined(HAVE_GPU_CUDA) && defined(HAVE_FC_ISO_C_BINDING)

 if (option==0.or.option==2) then
   call free_gpu_fourwf()
 end if
 if (option==1.or.option==2) then
   call free_nonlop_gpu()

#if defined(HAVE_KOKKOS)
   call dealloc_nonlop_kokkos()
#endif

 end if ! option 1 or 2

#else

 if (.false.) then
   write(std_out,*) option
 end if

#endif

end subroutine dealloc_hamilt_gpu
!!***

!!****f* ABINIT/alloc_nonlop_kokkos
!! NAME
!! alloc_hamilt_gpu
!!
!! FUNCTION
!! allocate several memory pieces on a GPU device for the application
!! of Hamiltonian using a GPU
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  nattyp(ntypat)= # atoms of each type.
!!  npwin is the number of plane waves for vectin
!!  npwout is the number of plane waves for vectout
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  use_gpu_cuda= 0 or 1 to know if we use cuda
!!
!! OUTPUT
!!  (no output - only allocation on GPU, member of gemm_nonlop_kokkos)
!!
!! SOURCE

subroutine alloc_nonlop_kokkos(dtset,&
  &                            psps,&
  &                            npwin,&
  &                            npwout,&
  &                            atindx1,&
  &                            nattyp,&
  &                            use_gpu_cuda)

  !Arguments ------------------------------------
  !scalars
  type(dataset_type),        intent(in) :: dtset
  type(pseudopotential_type),intent(in) :: psps
  integer,                   intent(in) :: npwin
  integer,                   intent(in) :: npwout
  integer,                   intent(in) :: atindx1(dtset%natom)
  !arrays
  integer,                   intent(in) :: nattyp(dtset%ntypat)
  !integer,                   intent(in) :: npwarr(dtset%nkpt)

  ! other
  integer,                   intent(in) :: use_gpu_cuda

  integer :: cplex
  integer :: nprojs
  integer :: itypat

  real(dp) :: allocated_size_bytes
  character(len=500)    :: message

  allocated_size_bytes = 0.

  cplex=2; !if (istwf_k>1) cplex=1

  ! compute nprojs
  nprojs = 0
  do itypat = 1,dtset%ntypat
    nprojs = nprojs + count(psps%indlmn(3,:,itypat)>0) * nattyp(itypat)
  end do


  !! allocate memory on device

  if (gemm_nonlop_kokkos % allocated .eqv. .false.) then
    ! These will store the non-local factors for vectin, svectout and vectout respectively
    ABI_MALLOC_CUDA(gemm_nonlop_kokkos%    projections_gpu, INT(cplex, c_size_t) * nprojs * dtset%nspinor*dtset%bandpp * dp)
    ABI_MALLOC_CUDA(gemm_nonlop_kokkos%  s_projections_gpu, INT(cplex, c_size_t) * nprojs * dtset%nspinor*dtset%bandpp * dp)
    ABI_MALLOC_CUDA(gemm_nonlop_kokkos%vnl_projections_gpu, INT(2    , c_size_t) * nprojs * dtset%nspinor*dtset%bandpp * dp)

    allocated_size_bytes = allocated_size_bytes + (2*cplex+2)*nprojs * dtset%nspinor*dtset%bandpp * dp

    ABI_MALLOC_CUDA(gemm_nonlop_kokkos%  vectin_gpu,  INT(2, c_size_t) * npwin  * dtset%nspinor*dtset%bandpp * dp)
    ABI_MALLOC_CUDA(gemm_nonlop_kokkos%  vectout_gpu, INT(2, c_size_t) * npwout * dtset%nspinor*dtset%bandpp * dp)
    ABI_MALLOC_CUDA(gemm_nonlop_kokkos% svectout_gpu, INT(2, c_size_t) * npwout * dtset%nspinor*dtset%bandpp * dp)

    allocated_size_bytes = allocated_size_bytes + &
      & 2 * (npwin+2*npwout) * dtset%nspinor * dtset%bandpp * dp

    gemm_nonlop_kokkos % allocated = .true.

    write(message,*) '  alloc_nonlop_kokkos allocated ', allocated_size_bytes*1e-6, ' MBytes on device memory'
    call wrtout(std_out,message,'COLL')

  end if

end subroutine alloc_nonlop_kokkos
!!***

!!****f* ABINIT/dealloc_nonlop_kokkos
!! NAME
!! dealloc_nonlop_kokkos
!!
!! FUNCTION
!! deallocate several memory pieces on a GPU device used for the application
!! of nonlop operators using Kokkos implementation
!!
!!
!! OUTPUT
!!  (no output - only deallocation on GPU)
!!
!! SOURCE

subroutine dealloc_nonlop_kokkos()

  if (gemm_nonlop_kokkos % allocated) then
    ABI_FREE_CUDA(gemm_nonlop_kokkos%    projections_gpu)
    ABI_FREE_CUDA(gemm_nonlop_kokkos%  s_projections_gpu)
    ABI_FREE_CUDA(gemm_nonlop_kokkos%vnl_projections_gpu)

    ABI_FREE_CUDA(gemm_nonlop_kokkos% vectin_gpu)
    ABI_FREE_CUDA(gemm_nonlop_kokkos% vectout_gpu)
    ABI_FREE_CUDA(gemm_nonlop_kokkos%svectout_gpu)

    gemm_nonlop_kokkos % allocated = .false.
  end if

end subroutine dealloc_nonlop_kokkos
!!***

end module m_alloc_hamilt_gpu
!!***
