!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_alloc_hamilt_gpu
!! NAME
!!  m_alloc_hamilt_gpu
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2019 ABINIT group (MT, FDahm)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_alloc_hamilt_gpu

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_abicore
 use m_xmpi

#if defined HAVE_GPU_CUDA
 use m_gpu_toolbox
#endif

 implicit none

 private
!!***

 public :: alloc_hamilt_gpu
 public :: dealloc_hamilt_gpu
!!***

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
!! PARENTS
!!      gstate,respfn
!!
!! CHILDREN
!!      free_gpu_fourwf,free_nonlop_gpu
!!
!! SOURCE

subroutine alloc_hamilt_gpu(atindx1,dtset,gprimd,mpi_enreg,nattyp,npwarr,option,psps,use_gpu_cuda)

 implicit none

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
 integer ::npwarr_tmp(dtset%nkpt)
#endif

! *************************************************************************

 if (use_gpu_cuda==0) return

#if defined HAVE_GPU_CUDA
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
!  Initialize gpu data needed in fourwf
!  ndat=bandpp when paral_kgb=1
   if(mpi_enreg%paral_kgb==1) then
     call alloc_gpu_fourwf(dtset%ngfft,dtset%bandpp,npw_max_loc,npw_max_loc)
   else
     call alloc_gpu_fourwf(dtset%ngfft,1,npw_max_loc,npw_max_loc)
   end if
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
   call alloc_nonlop_gpu(npw_max_nonloc,npw_max_nonloc,dtset%nspinor,dtset%natom,&
&   dtset%ntypat,psps%lmnmax,psps%indlmn,nattyp,atindx1,gprimd,&
&   dimffnl_max,dimffnl_max,dimekb1_max,dimekb2_max)
 end if
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
!! PARENTS
!!      gstate,respfn
!!
!! CHILDREN
!!      free_gpu_fourwf,free_nonlop_gpu
!!
!! SOURCE

subroutine dealloc_hamilt_gpu(option,use_gpu_cuda)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option,use_gpu_cuda
!arrays

!Local variables-------------------------------

! *************************************************************************

 if (use_gpu_cuda==0) return

#if defined HAVE_GPU_CUDA
 if (option==0.or.option==2) then
   call free_gpu_fourwf()
 end if
 if (option==1.or.option==2) then
   call free_nonlop_gpu()
 end if
#else
 if (.false.) then
   write(std_out,*) option
 end if
#endif

end subroutine dealloc_hamilt_gpu
!!***

end module m_alloc_hamilt_gpu
!!***
