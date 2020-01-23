!!****m* ABINIT/defs_rectypes
!! NAME
!! defs_rectypes
!!
!! FUNCTION
!! This module contains definitions of all structured datatype for recursion
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module defs_rectypes

 use defs_basis
 use m_abicore
 use defs_abitypes,  only : MPI_type
 use m_pawfgr, only : pawfgr_type

 implicit none
!!***

!!****t* defs_rectypes/vec_int
!! NAME
!! vec_int
!!
!! FUNCTION
!! This structured datatype contains 3d vector of indexes
!!
!! SOURCE

 type,public :: vec_int
  integer  :: x,y,z
 end type vec_int

!!***

!!****t* defs_rectypes/metricrec_type
!! NAME
!! metricrec_type
!!
!! FUNCTION
!! This structured datatype used in recursion_type containing
!! information concerning the infenitesimal metrics
!!
!! SOURCE


 type,public :: metricrec_type
!  Real scalars
  real(dp)  :: ucvol     !--Infinitesimal volume
! Integer arrays
  integer,allocatable  :: gcart(:,:)
! Real (real(dp)) arrays
  real(dp)  :: tr(3)     !--Trace of the metrics: (grid step)
  real(dp)  :: rmet(3,3) !--Infinitesimal metrics in recursion

 end type metricrec_type
!!***

!Structures
!!***

!!****t* defs_rectypes/nlpsprec_type
!! NAME
!! nlpsprec_type
!!
!! FUNCTION
!! This structured datatype used in recursion_type containing
!! information concerning psp
!!
!! SOURCE

 type,public :: nlpsprec_type! Integer scalars

  integer :: lmnmax
  ! in recursion, it contains the max counting-index for lmn, it
  ! corresponds to the ubound of the second index of indlmn

  integer :: intlen !nl_ngfft is the maximum, on all psps, of the
  ! interaction length of non-local potential (in grid step).

  integer :: npsp  ! Number of types of pseudopotentials

  logical :: nlpsp !calculate the non-localc pspeudo part in recursion

! Integer arrays
  integer,allocatable :: pspinfo(:,:) !--for any typat:  (momang,typat)=nproj
  integer,allocatable :: indlmn(:,:,:)
  ! indlmn(6,lmnmax,ntypat)
  ! For each type of psp,
  ! array giving l,m,n,lm,ln,spin for i=lmn (if useylm=1)


! Real (real(dp)) arrays
  real(dp),allocatable :: radii(:,:)   !--Local radii of nl psp
                                  !  part radii(l+1,psp)

  real(dp),allocatable :: projec(:,:,:)         !--p(r,lmn,typat), is the
  !product of Y_lm(r)*nlproj_ln(r) calculated on the grid for any psps.

  real(dp),allocatable :: mat_exp_psp_nl(:,:,:,:) !--Matrix of exponential
                                                !  of non-local psp
                                                ! mat_exp_psp_nl(nproj,nproj,l+1,psp)

  real(dp),allocatable :: eival(:,:,:) !--Matrix of exponential(eigenvalues)
  real(dp),allocatable :: eivec(:,:,:,:) !--Matrix of exponential(eigenvectors)

 end type nlpsprec_type
!!***

!Structures
!!***

!!****t* defs_rectypes/ recparall_type
!! NAME
!! recparall_type
!!
!! FUNCTION
!! This structured datatype for distribution of work in recursion
!!
!!
!! SOURCE

 type,public :: recparall_type

  !--Variables for parallelization in recursion
  integer :: min_pt
  integer :: max_pt
  integer :: ntranche
  integer :: npt !--Number of points to calculate for a proc

! Integer arrays
  integer,allocatable ::  displs(:) !--For mpi
  integer,allocatable ::  vcount(:)

! Structured datatypes
  type(vec_int)        :: pt0    !--Initial point to calculate
  type(vec_int)        :: pt1    !--Final point to calculate

 end type recparall_type

!Structures
!!***

!!****t* defs_rectypes/recGPU_type
!! NAME
!! recGPU_type
!!
!! FUNCTION
!! This structured datatype for recursion with GPU
!!
!!
!! SOURCE

 type,public :: recGPU_type

  integer :: nptrec !--number of points for recursion on GPUon
  !any devices (supposing that all devices are the same type)

  integer,pointer :: map(:)  !--The map of devices GPU with respect to
  !                       cpu (-1 if there are not gpu)

  ! Structured datatypes
  type(recparall_type) :: par    !--For distribution of work on procs

 end type recGPU_type

!Structures
!!***

!!****t* defs_rectypes/recursion_type
!! NAME
!! recursion_type
!!
!! FUNCTION
!! This structured datatype contains some variable used in the recursion
!!
!! SOURCE


 type,public :: recursion_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar


  integer :: quitrec !--To quit by recursion

  integer :: min_nrec !--Minimal number of recursions

  integer :: nfftrec !--Recursion nfft

  integer :: ngpu !--total number of GPUs

  integer :: gpudevice !--it corresponds to the GPU-device setted for this proc
  !if negative: no GPU-device is associated

  integer :: load      !--if egual to 0 (default), load on proc-gpu is
  !homogeneous and %par is egual to %GPU%par (or %GPU is not affected)
  !otherwise load=1
  integer :: tp        !--topology of the machine:
  !0: 1 cpu;
  !1: n cpu;
  !2: 1 cpu 1 gpu;
  !3: n cpu n gpu
  !4: n cpu > m gpu;
  !5: n cpu < m gpu (not implemented)

  logical :: debug     !--Debugging variable
  logical :: tronc     !--True if troncation (ngfftrec/=ngfft) is used


! Integer arrays
  integer :: ngfftrec(18) !--Recursion ngfft
! Real (real(dp)) scalar
  real(dp) :: efermi

  real(dp),allocatable :: ZT_p(:,:)  !--Is the Fourier transform of the
                                      !  green kernel

! Structured datatypes
  type(recparall_type)   :: par    !--For distribution of work on procs
  type(pawfgr_type)      :: pawfgr !--For double grid system
  type(MPI_type),pointer :: mpi    !--For recursion mpi
  type(nlpsprec_type)    :: nl     !--For nlpsp
  type(metricrec_type)   :: inf    !--For metrics
#ifdef HAVE_GPU_CUDA
  type(recGPU_type)      :: GPU  !--information concerning cuda
#endif
 end type recursion_type

end module defs_rectypes
!!***
