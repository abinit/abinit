!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_fft
!! NAME
!!  m_fft
!!
!! FUNCTION
!!  This module provides driver routines for sequential FFTs (OpenMP threads are supported).
!!  It also defines generic interfaces for single or double precision arrays.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2018 ABINIT group (MG, MM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_fft

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xomp
 use m_xmpi
 use m_cplxtools
 use m_cgtools
 use m_sgfft
 use m_sg2002
 use m_fftw3
 use m_dfti
 use iso_c_binding

 use defs_abitypes,   only : MPI_type
 use defs_fftdata,    only : mg
 use m_time,          only : cwtime
 use m_numeric_tools, only : r2c
 use m_fstrings,      only : sjoin, itoa
 use m_blas,          only : xscal
 use m_fftcore,       only : get_cache_kb, kpgsph, get_kg, sphere_fft1, sphere, change_istwfk,&
&                            fftalg_info, fftalg_has_mpi, print_ngfft
 use m_mpinfo,        only : destroy_mpi_enreg
 use m_distribfft,    only : distribfft_type, init_distribfft, destroy_distribfft

 implicit none

 private

 public :: fftbox_execute       ! Driver for FFTs on the full box (complex-to-complex version, operated on complex arrays)
 public :: fft_ug               ! Driver for zero-padded FFTs u(g) --> u(r)
 public :: fft_ur               ! Driver for zero-padded FFTs u(r) --> u(g)
 public :: fftpad               ! Driver for (low-level) zero-padded FFTs, note that fft_ug is the preferred interface.
 public :: fft_poisson          ! Solve the poisson equation in G-space starting from n(r).

 public :: fft_use_lib_threads
 public :: fft_allow_ialltoall  ! Allow the use of non-blocking IALLOTOALL in MPI-FFTs algorithms

! Driver routines for MPI version.
 public :: fourdp_mpi           ! MPI FFT of densities/potentials on the full box.
 public :: fourwf_mpi           ! specialized MPI-FFT for wavefunctions.
 !public :: fftmpi_u

 interface fftbox_execute
   module procedure fftbox_execute_ip_spc
   module procedure fftbox_execute_ip_dpc
   module procedure fftbox_execute_op_spc
   module procedure fftbox_execute_op_dpc
 end interface fftbox_execute

 interface fft_ug
   !module procedure fft_ug_dp   TODO
   module procedure fft_ug_spc
   module procedure fft_ug_dpc
 end interface fft_ug
 !public :: fft_ug_spc, fft_ug_dpc

 interface fft_ur
   !module procedure fft_ur_dp      TODO
   module procedure fft_ur_spc
   module procedure fft_ur_dpc
 end interface fft_ur
 !public :: fft_ur_spc, fft_ur_dpc

 interface fftpad
   !module procedure fftpad_dp
   module procedure fftpad_spc
   module procedure fftpad_dpc
 end interface fftpad
!!***

!----------------------------------------------------------------------

!!****t* m_fft/fftbox_plan3_t
!! NAME
!! fftbox_plan3_t
!!
!! FUNCTION
!!  Stores the options passed to the fftbox_ routines.
!!
!! SOURCE

 type,public :: fftbox_plan3_t
   private
   integer :: fftalg=112      ! Flag defining the library to call.
   integer :: fftcache=16     ! Size of the cache (kB). Only used in SG routines.
   integer :: isign=0         ! Sign of the exponential in the FFT
   integer :: nfft            ! Total number of points in the FFT box.
   integer :: ldxyz=-1        ! Physical dimension of the array to transform
   integer :: ndat=-1         ! Number of FFTs associated to the plan.
   !integer :: nthreads=-1    ! The number of threads associated to the plan.
   integer :: dims(3)=-1      ! The number of FFT divisions.
   integer :: embed(3)=-1     ! Leading dimensions of the input,output arrays.
 end type fftbox_plan3_t

 public :: fftbox_plan3       ! Basic interface to create the plan.
 public :: fftbox_plan3_many  ! Advanced interface
 public :: fftbox_plan3_init  ! Low-level constructor
!!***

!----------------------------------------------------------------------

!unitary tests
 public :: fftbox_utests          ! Unit tests for FFTs on the full box.
 public :: fftu_utests            ! Unit tests for the FFTs of wavefunctions.
 public :: fftbox_mpi_utests      ! Unit tests for the MPI-FFTs on the full box.
 public :: fftu_mpi_utests        ! Unit tests for the MPI-FFTs of the wavefunctions.
!!***

 ! Flag used to enable/disable the use of non-blocking IALLTOALL
#ifdef HAVE_MPI_IALLTOALL
 logical,save,private :: ALLOW_IALLTOALL = .True.
#else
 logical,save,private :: ALLOW_IALLTOALL = .False.
#endif

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fft_allow_ialltoall
!! NAME
!!  fft_allow_ialltoall
!!
!! FUNCTION
!!  Allow the use of non-blocking IALLOTOALL in MPI-FFTs algorithms
!!  Mainly used for profiling purposes.
!!
!! PARENTS
!!      m_argparse
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fft_allow_ialltoall(bool)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft_allow_ialltoall'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 logical,intent(in) :: bool

! *************************************************************************

 ALLOW_IALLTOALL = bool
#ifndef HAVE_MPI_IALLTOALL
 ALLOW_IALLTOALL = .False.
#endif

end subroutine fft_allow_ialltoall
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftbox_plan3
!! NAME
!!  fftbox_plan3
!!
!! FUNCTION
!!  Basic interface to construct fftbox_plan3_t
!!
!! INPUTS
!!
!! PARENTS
!!      debug_tools
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fftbox_plan3(plan,dims,fftalg,isign)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftbox_plan3'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,isign
 type(fftbox_plan3_t),intent(out) :: plan
!arrays
 integer,intent(in) :: dims(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: ndat1=1,fftcache0=0

! *************************************************************************

 call fftbox_plan3_init(plan,ndat1,dims,dims,fftalg,fftcache0,isign)

end subroutine fftbox_plan3
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftbox_plan3_many
!! NAME
!!  fftbox_plan3_many
!!
!! FUNCTION
!!  Advanced interface to construct fftbox_plan3_t
!!
!! INPUTS
!!  See fftbox_plan3_t
!!
!! PARENTS
!!      m_fft,m_fft_prof,m_oscillators,m_pawpwij,m_shirley
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fftbox_plan3_many(plan,ndat,dims,embed,fftalg,isign)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftbox_plan3_many'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,isign,ndat
 type(fftbox_plan3_t),intent(out) :: plan
!arrays
 integer,intent(in) :: dims(3),embed(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: fftcache0=0

! *************************************************************************

 call fftbox_plan3_init(plan,ndat,dims,embed,fftalg,fftcache0,isign)

end subroutine fftbox_plan3_many
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftbox_plan3_init
!! NAME
!!  fftbox_plan3_init
!!
!! FUNCTION
!!  Initialize the plan with the options passed to the fttbox_ routines.
!!  Low-level constructor.
!!
!! INPUTS
!!  See fftbox_plan3_t for the meaning of the different arguments.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fftbox_plan3_init(plan,ndat,dims,embed,fftalg,fftcache,isign)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftbox_plan3_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndat,fftalg,fftcache,isign !,nthreads
 type(fftbox_plan3_t),intent(out) :: plan
!arrays
 integer,intent(in) :: dims(3),embed(3)

! *************************************************************************

 plan%ndat     = ndat
 plan%dims     = dims                       !ngfft(1:3)
 plan%embed    = embed                      !ngfft(4:6)
 plan%fftalg   = fftalg                     !ngfft(7)
 if (fftcache > 0) plan%fftcache = fftcache !ngfft(8)
 plan%isign    = isign
 !plan%nthreads = nthreads

 plan%nfft     = PRODUCT(plan%dims)
 plan%ldxyz    = PRODUCT(plan%embed)

end subroutine fftbox_plan3_init
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftbox_execute_ip_spc
!! NAME
!!  fftbox_execute_ip_spc
!!
!! FUNCTION
!!  In-place FFT transform of complex array.
!!  Call (FFTW3|DFTI) routines if available, otherwise we fallback to SG routines
!!  TARGET: spc arrays
!!
!! INPUTS
!!  plan<fftbox_plan3_t>=Structure with the parameters defining the transform.
!!
!! SIDE EFFECTS
!!  ff(plan%ldxyz*plan%ndat) =
!!    In input: the data to transform.
!!    Changed in output, filled with the FFT results.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fftbox_execute_ip_spc(plan,ff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftbox_execute_ip_spc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(fftbox_plan3_t),intent(in) :: plan
!arrays
 complex(spc),intent(inout) :: ff(plan%ldxyz*plan%ndat)

! *************************************************************************

#include "fftbox_ip_driver.finc"

end subroutine fftbox_execute_ip_spc
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftbox_execute_ip_dpc
!! NAME
!!  fftbox_execute_ip_dpc
!!
!! FUNCTION
!!  In-place FFT transform of complex arrays
!!  Call (FFTW3|DFTI) routines if available, otherwise fallback to SG routines
!!  TARGET: dpc arrays
!!
!! INPUTS
!!  plan<fftbox_plan3_t>=Structure with the parameters defining the transform.
!!
!! SIDE EFFECTS
!!  ff(plan%ldxyz*plan%ndat) =
!!    In input: the data to transform.
!!    Changed in output, filled with the FFT results.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fftbox_execute_ip_dpc(plan,ff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftbox_execute_ip_dpc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(fftbox_plan3_t),intent(in) :: plan
!arrays
 complex(dpc),intent(inout) :: ff(plan%ldxyz*plan%ndat)

! *************************************************************************

#include "fftbox_ip_driver.finc"

end subroutine fftbox_execute_ip_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftbox_execute_op_spc
!! NAME
!!  fftbox_execute_op_spc
!!
!! FUNCTION
!!  Out-of-place FFT transform of complex arrays.
!!  Call (FFTW3|DFTI) routines if available, otherwise fallback to SG routines
!!  TARGET: spc arrays
!!
!! INPUTS
!! plan<fftbox_plan3_t>=Structure with the parameters defining the transform.
!! ff(plan%ldxyz*plan%ndat)=The input array to be transformed.
!!
!! OUTPUT
!!  gg(plan%ldxyz*plan%ndat)= The FFT results.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fftbox_execute_op_spc(plan,ff,gg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftbox_execute_op_spc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(fftbox_plan3_t),intent(in) :: plan
!arrays
 complex(spc),intent(in) :: ff(plan%ldxyz*plan%ndat)
 complex(spc),intent(inout) :: gg(plan%ldxyz*plan%ndat)

! *************************************************************************

#include "fftbox_op_driver.finc"

end subroutine fftbox_execute_op_spc
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftbox_execute_op_dpc
!! NAME
!!  fftbox_execute_op_dpc
!!
!! FUNCTION
!!  Out-of-place FFT transform of complex arrays.
!!  Call (FFTW3|DFTI) routines if available, otherwise fallback to SG routines
!!  TARGET: dpc arrays
!!
!! INPUTS
!! plan<fftbox_plan3_t>=Structure with the parameters defining the transform.
!! ff(plan%ldxyz*plan%ndat)=The input array to be transformed.
!!
!! OUTPUT
!!  gg(plan%ldxyz*plan%ndat)= The FFT results.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fftbox_execute_op_dpc(plan,ff,gg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftbox_execute_op_dpc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(fftbox_plan3_t),intent(in) :: plan
!arrays
 complex(dpc),intent(in) :: ff(plan%ldxyz*plan%ndat)
 complex(dpc),intent(inout) :: gg(plan%ldxyz*plan%ndat)

! *************************************************************************

#include "fftbox_op_driver.finc"

end subroutine fftbox_execute_op_dpc
!!***

!----------------------------------------------------------------------
#if 0

!!****f* m_fft/fft_ug_dp
!! NAME
!! fft_ug_dp
!!
!! FUNCTION
!! Driver routine for G-->R transform of wavefunctions with zero-padded FFT.
!! TARGET: double precision real arrays
!!
!! INPUTS
!! npw_k=number of plane waves for this k-point.
!! nfft=Number of FFT points.
!! nspinor=number of spinorial components
!! ndat=Numer of wavefunctions to transform.
!! mgfft=Max number of FFT divisions
!! ngfft(18)=information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! istwfk=Option describing the storage of the wavefunction. (at present must be 1)
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!! gbound_k_k(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!! ug(npw_k*nspinor*ndat)=wavefunctions in reciprocal space
!!
!! OUTPUT
!!  ur(nfft*nspinor*ndat)=wavefunctions in real space.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fft_ug_dp(npw_k,nfft,nspinor,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,ug,ur)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft_ug_dp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,nfft,nspinor,istwf_k,mgfft,ndat
!arrays
 integer,intent(in) :: ngfft(18),gbound_k(2*mgfft+8,2),kg_k(3,npw_k)
 real(dp),intent(in) :: ug(2*npw_k*nspinor*ndat)
 real(dp),intent(out) :: ur(2*nfft*nspinor*ndat)

! *************************************************************************

 MSG_ERROR("This is a Stub")
!#include "fftug_driver.finc"

end subroutine fft_ug_dp
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_fft/fft_ug_spc
!! NAME
!! fft_ug_spc
!!
!! FUNCTION
!! Driver routine for G-->R transform of wavefunctions with zero-padded FFT.
!! TARGET: single precision arrays
!!
!! INPUTS
!! npw_k=number of plane waves for this k-point.
!! nfft=Number of FFT points.
!! nspinor=number of spinorial components
!! ndat=Numer of wavefunctions to transform.
!! mgfft=Max number of FFT divisions
!! ngfft(18)=information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! istwfk=Option describing the storage of the wavefunction. (at present must be 1)
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!! gbound_k_k(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!! ug(npw_k*nspinor*ndat)=wavefunctions in reciprocal space
!!
!! OUTPUT
!!  ur(nfft*nspinor*ndat)=wavefunctions in real space.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fft_ug_spc(npw_k,nfft,nspinor,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,ug,ur)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft_ug_spc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,nfft,nspinor,istwf_k,mgfft,ndat
!arrays
 integer,intent(in) :: ngfft(18),gbound_k(2*mgfft+8,2),kg_k(3,npw_k)
 complex(spc),intent(in) :: ug(npw_k*nspinor*ndat)
 complex(spc),intent(out) :: ur(nfft*nspinor*ndat)

! *************************************************************************

#include "fftug_driver.finc"

end subroutine fft_ug_spc
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fft_ug_dpc
!! NAME
!! fft_ug_dpc
!!
!! FUNCTION
!! Driver routine for G-->R transform of wavefunctions with zero-padded FFT.
!! TARGET: double precision arrays
!!
!! INPUTS
!! npw_k=number of plane waves for this k-point.
!! nfft=Number of FFT points.
!! nspinor=number of spinorial components
!! ndat=Numer of wavefunctions to transform.
!! mgfft=Max number of FFT divisions
!! ngfft(18)=information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! istwfk=Option describing the storage of the wavefunction. (at present must be 1)
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!! gbound_k_k(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!! ug(npw_k*nspinor*ndat)=wavefunctions in reciprocal space
!!
!! OUTPUT
!!  ur(nfft*nspinor*ndat)=wavefunctions in real space.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fft_ug_dpc(npw_k,nfft,nspinor,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,ug,ur)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft_ug_dpc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,nfft,nspinor,istwf_k,mgfft,ndat
!arrays
 integer,intent(in) :: ngfft(18),gbound_k(2*mgfft+8,2),kg_k(3,npw_k)
 complex(dpc),intent(in) :: ug(npw_k*nspinor*ndat)
 complex(dpc),intent(out) :: ur(nfft*nspinor*ndat)

! *************************************************************************

#include "fftug_driver.finc"

end subroutine fft_ug_dpc
!!***

!----------------------------------------------------------------------
#if 0

!!****f* m_fft/fft_ur_dp
!! NAME
!! fft_ur_dp
!!
!! FUNCTION
!! Compute ndat zero-padded FFTs from R- to G-space .
!! Mainly used for the transform of wavefunctions.
!! TARGET: dp real arrays
!!
!! INPUTS
!! npw_k=number of plane waves for this k-point.
!! nfft=Number of FFT points.
!! nspinor=number of spinorial components
!! ndat=Number of wavefunctions to transform.
!! mgfft=Max number of FFT divisions
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! istwfk=Option describing the storage of the wavefunction. (at present must be 1)
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!! gbound_k(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!!
!! SIDE EFFECTS
!!  ur(nfft*nspinor*ndat)=In input: wavefunctions in real space
!!                        Destroyed in output. Do not use ur anymore!
!!
!! OUTPUT
!!  ug(npw_k*nspinor*ndat)=wavefunctions in reciprocal space given on the G-sphere.
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fft_ur_dp(npw_k,nfft,nspinor,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,ur,ug)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft_ur_dp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,nfft,nspinor,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: ngfft(18),gbound_k(2*mgfft+8,2)
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(inout) :: ur(2*nfft*nspinor*ndat)
 real(dp),intent(out) :: ug(2*npw_k*nspinor*ndat)

! *************************************************************************

 MSG_ERROR("This is a Stub")

!#include "fftur_driver.finc"

end subroutine fft_ur_dp
!!***
#endif

!----------------------------------------------------------------------

!!****f* m_fft/fft_ur_spc
!! NAME
!! fft_ur_spc
!!
!! FUNCTION
!! Compute ndat zero-padded FFTs from R- to G-space .
!! Mainly used for the transform of wavefunctions.
!! TARGET: spc complex arrays
!!
!! INPUTS
!! npw_k=number of plane waves for this k-point.
!! nfft=Number of FFT points.
!! nspinor=number of spinorial components
!! ndat=Number of wavefunctions to transform.
!! mgfft=Max number of FFT divisions
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! istwfk=Option describing the storage of the wavefunction. (at present must be 1)
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!! gbound_k(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!!
!! SIDE EFFECTS
!!  ur(nfft*nspinor*ndat)=In input: wavefunctions in real space
!!                        Destroyed in output. Do not use ur anymore!
!!
!! OUTPUT
!!  ug(npw_k*nspinor*ndat)=wavefunctions in reciprocal space given on the G-sphere.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fft_ur_spc(npw_k,nfft,nspinor,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,ur,ug)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft_ur_spc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,nfft,nspinor,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: ngfft(18),gbound_k(2*mgfft+8,2)
 integer,intent(in) :: kg_k(3,npw_k)
 complex(spc),intent(inout) :: ur(nfft*nspinor*ndat)
 complex(spc),intent(out) :: ug(npw_k*nspinor*ndat)

! *************************************************************************

#include "fftur_driver.finc"

end subroutine fft_ur_spc
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fft_ur_dpc
!! NAME
!! fft_ur_dpc
!!
!! FUNCTION
!! Compute ndat zero-padded FFTs from R- to G-space .
!! Mainly used for the transform of wavefunctions.
!! TARGET: dpc complex arrays
!!
!! INPUTS
!! npw_k=number of plane waves for this k-point.
!! nfft=Number of FFT points.
!! nspinor=number of spinorial components
!! ndat=Number of wavefunctions to transform.
!! mgfft=Max number of FFT divisions
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! istwfk=Option describing the storage of the wavefunction. (at present must be 1)
!! kg_k(3,npw_k)=G-vectors in reduced coordinates
!! gbound_k(2*mgfft+8,2)=Table for padded-FFT. See sphereboundary.
!!
!! SIDE EFFECTS
!!  ur(nfft*nspinor*ndat)=In input: wavefunctions in real space
!!                        Destroyed in output. Do not use ur anymore!
!!
!! OUTPUT
!!  ug(npw_k*nspinor*ndat)=wavefunctions in reciprocal space given on the G-sphere.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fft_ur_dpc(npw_k,nfft,nspinor,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,ur,ug)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft_ur_dpc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k,nfft,nspinor,ndat,istwf_k,mgfft
!arrays
 integer,intent(in) :: ngfft(18),gbound_k(2*mgfft+8,2)
 integer,intent(in) :: kg_k(3,npw_k)
 complex(dpc),intent(inout) :: ur(nfft*nspinor*ndat)
 complex(dpc),intent(out) :: ug(npw_k*nspinor*ndat)

! *************************************************************************

#include "fftur_driver.finc"

end subroutine fft_ur_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftpad_spc
!! NAME
!!  fftpad_spc
!!
!! FUNCTION
!!  Driver routine used to transform COMPLEX arrays using 3D zero-padded FFTs.
!!  TARGET: SPC arrays
!!
!! INPUTS
!!  ngfft(18)=Info on the 3D FFT.
!!  nx,ny,nz=Logical dimensions of the FFT mesh.
!!  ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!!  ndat=Number of FFTs
!!  mgfft=MAX(nx,ny,nz), only used to dimension gbound
!!  isign=The sign of the transform.
!!  gbound(2*mgfft+8,2)= The boundaries of the basis sphere of G vectors at a given k-point. See sphereboundary for more info.
!!
!! SIDE EFFECTS
!!  ff(ldx*ldy*ldz*ndat)=
!!    input: The array with the data to be transformed.
!!    output: The results of the FFT.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fftpad_spc(ff,ngfft,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign,gbound)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftpad_spc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign
!arrays
 integer,intent(in) :: ngfft(18),gbound(2*mgfft+8,2)
 complex(spc),target,intent(inout) :: ff(ldx*ldy*ldz*ndat)

!Local variables-------------------------------
!scalars
 integer :: fftalg,fftalga,fftalgc,ncount,p
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: fofr(:,:),ftarr(:,:)

! *************************************************************************

 fftalg=ngfft(7); fftalga=fftalg/100; fftalgc=MOD(fftalg,10)

 select case (fftalga)

 case (FFT_FFTW3)
   call fftw3_fftpad(ff,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign,gbound)

 case (FFT_DFTI)
   call dfti_fftpad(ff,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign,gbound)

 case (FFT_SG)
   ! Goedecker"s routines.
   ! TODO: sg_fftpad is not the fastest routine, here I should call sg_fftrisc but I need
   ! kg_kin that are not available in rho_tw_g, actually one should pass G-G0 due to
   ! the shift introduced by the umklapp.
   ncount = ldx*ldy*ldz*ndat

   ABI_MALLOC(fofr,(2,ldx*ldy*ldz*ndat))
   ABI_MALLOC(ftarr,(2,ldx*ldy*ldz*ndat))

   do p=1,ldx*ldy*ldz*ndat
     fofr(1,p) = REAL(ff(p))
     fofr(2,p) = AIMAG(ff(p))
   end do

   call sg_fftpad(ngfft(8),mgfft,nx,ny,nz,ldx,ldy,ldz,ndat,gbound,isign,fofr,ftarr)
   !
   ! Copy the results.
   do p=1,ldx*ldy*ldz*ndat
     ff(p) = CMPLX(ftarr(1,p), ftarr(2,p))
   end do

   if (isign==-1) then  ! Here there might be numerical exceptions due to the holes.
     ff = ff/(nx*ny*nz)
   end if

   ABI_FREE(ftarr)
   ABI_FREE(fofr)

 case default
   write(msg,'(a,i0,a)')"fftalga = ", fftalga," not coded "
   MSG_ERROR(msg)
 end select

end subroutine fftpad_spc
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftpad_dpc
!! NAME
!!  fftpad_dpc
!!
!! FUNCTION
!!  Driver routine used to transform COMPLEX arrays using 3D zero-padded FFTs.
!!  TARGET: DPC arrays
!!
!! INPUTS
!!  ngfft(18)=Info on the 3D FFT.
!!  nx,ny,nz=Logical dimensions of the FFT mesh.
!!  ldx,ldy,ldz=Physical dimension of the f array (to avoid cache conflicts).
!!  ndat=Number of FFTs
!!  mgfft=MAX(nx,ny,nz), only used to dimension gbound
!!  isign=The sign of the transform.
!!  gbound(2*mgfft+8,2)= The boundaries of the basis sphere of G vectors at a given k-point. See sphereboundary for more info.
!!
!! SIDE EFFECTS
!!  ff(ldx*ldy*ldz*ndat)=
!!    input: The array with the data to be transformed.
!!    output: The results of the FFT.
!!
!! PARENTS
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fftpad_dpc(ff,ngfft,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign,gbound)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftpad_dpc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign
!arrays
 integer,intent(in) :: ngfft(18),gbound(2*mgfft+8,2)
 complex(dpc),target,intent(inout) :: ff(ldx*ldy*ldz*ndat)

!Local variables-------------------------------
!scalars
 integer :: fftalg,fftalga,fftalgc,ncount
 integer :: ivz !vz_d
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: fofr(:,:,:,:,:)
 real(dp),allocatable :: fofrvz(:,:) !vz_d
 real(dp),ABI_CONTIGUOUS pointer :: fpt_ftarr(:,:,:,:,:)

! *************************************************************************

 fftalg=ngfft(7); fftalga=fftalg/100; fftalgc=MOD(fftalg,10)

 select case (fftalga)

 case (FFT_FFTW3)
   call fftw3_fftpad(ff,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign,gbound)

 case (FFT_DFTI)
   call dfti_fftpad(ff,nx,ny,nz,ldx,ldy,ldz,ndat,mgfft,isign,gbound)

 case (FFT_SG)
   ! Goedecker"s routines.
   ! TODO: sg_fftpad is not the fastest routine, here I should call sg_fftrisc but I need
   ! kg_kin that are not available in rho_tw_g, actually one should pass G-G0 due to
   ! the shift introduced by the umklapp.
   ncount = ldx*ldy*ldz*ndat

   ABI_MALLOC(fofr, (2,ldx,ldy,ldz,ndat))
!  call ZCOPY(ncount,ff,1,fofr,1) !vz_d
!  call DCOPY(2*ncount,ff,1,fofr,1)  ! MG
   ! alternatif of ZCOPY from vz
   ABI_ALLOCATE(fofrvz,(2,ncount))     !vz_d
   do ivz=1,ncount                !vz_d
      fofrvz(1,ivz)= real(ff(ivz))  !vz_d
      fofrvz(2,ivz)=aimag(ff(ivz))  !vz_d
   end do                         !vz_d
   call DCOPY(2*ncount,fofrvz,1,fofr,1) !vz_d
   ABI_DEALLOCATE(fofrvz)             !vz_d

   call C_F_pointer(C_loc(ff),fpt_ftarr, shape=(/2,ldx,ldy,ldz,ndat/))

   !  MT nov. 2012: with xlf, need to put explicit boundaries for the 4th dimension
   !  this looks like a compiler bug...
   !if (size(fpt_ftarr)==2*size(ff)) then
   call sg_fftpad(ngfft(8),mgfft,nx,ny,nz,ldx,ldy,ldz,ndat,gbound,isign,fofr,fpt_ftarr)
   !else
   !  call sg_fftpad(ngfft(8),mgfft,nx,ny,nz,ldx,ldy,ldz,ndat,gbound,isign,fofr,fpt_ftarr(:,:,:,1:ldz,dat))
   !end if

   ABI_FREE(fofr)

   if (isign==-1) then  ! Here there might be numerical exceptions due to the holes.
     call xscal(ncount,one/(nx*ny*nz),ff,1)
   end if

 case default
   write(msg,'(a,i0,a)')"fftalga = ", fftalga," not coded "
   MSG_ERROR(msg)
 end select

end subroutine fftpad_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fft_poisson
!! NAME
!! fft_poisson
!!
!! FUNCTION
!!  Driver routine to solve the Poisson equation in G-space given the density, n(r),
!!  in real space of the FFT box.
!!
!! INPUTS
!! ngfft(18)=Info on the 3D FFT.
!! cplex=1 if fofr is real, 2 if fofr is complex
!! nx,ny,nz=Number of FFT points along the three directions.
!! ldx,ldy,ldz=Leading dimension of the array nr and vg.
!! ndat = Number of densities
!! vg(nx*ny*nz)=Potential in reciprocal space.
!!
!! SIDE EFFECTS
!! nr(cplex*ldx*ldy*ldz*ndat)
!!    input: n(r) (real or complex)
!!    output: the hartree potential in real space
!!
!! NOTES
!!   vg is given on the FFT mesh instead of the augmented mesh [ldx,ldy,ldz]
!!   in order to simplify the interface with the other routines operating of vg
!!
!! PARENTS
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fft_poisson(ngfft,cplex,nx,ny,nz,ldx,ldy,ldz,ndat,vg,nr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft_poisson'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nx,ny,nz,ldx,ldy,ldz,ndat
 integer,intent(in) :: ngfft(18)
!arrays
 real(dp),intent(inout) :: nr(cplex*ldx*ldy*ldz*ndat)
 real(dp),intent(in) :: vg(nx*ny*nz)

!Local variables-------------------------------
!scalars
 integer :: fftalga,fftcache

! *************************************************************************

 fftalga = ngfft(7)/100; fftcache = ngfft(8)

 select case (fftalga)
 case (FFT_SG, FFT_SG2002)
   ! Note: according to my tests fftalg 1xx is always faster than 4xx in sequential
   ! hence it does not make sense to provide a fallback for fftalg 4xx.
   ! Use external FFT libraries if you want to run at top-level speed.
   call sg_poisson(fftcache,cplex,nx,ny,nz,ldx,ldy,ldz,ndat,vg,nr)

 case (FFT_FFTW3)
   call fftw3_poisson(cplex,nx,ny,nz,ldx,ldy,ldz,ndat,vg,nr)

 !case (FFT_DFTI)
 !  call dfti_poisson(cplex,nx,ny,nz,ldx,ldy,ldz,ndat,vg,nr)

 case default
   MSG_BUG(sjoin("Wrong value for fftalga: ",itoa(fftalga)))
 end select

end subroutine fft_poisson
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fft_use_lib_threads
!! NAME
!! fft_use_lib_threads
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      fftprof
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fft_use_lib_threads(logvar)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fft_use_lib_threads'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: logvar

! *************************************************************************

 call dfti_use_lib_threads(logvar)
 call fftw3_use_lib_threads(logvar)

end subroutine fft_use_lib_threads
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftbox_utests
!! NAME
!!  fftbox_utests
!!
!! FUNCTION
!! Driver routine for base unitary tests of the FFT routines (sequential version).
!!
!! INPUTS
!! fftalg =fftalg input variable.
!! ndat = Number of transform to execute
!! nthreads = Number of OpenMP threads.
!! [unit]=Output Unit number (DEFAULT std_out)
!!
!! OUTPUT
!!  nfailed=number of failures.
!!
!! PARENTS
!!      fftprof
!!
!! SOURCE

function fftbox_utests(fftalg,ndat,nthreads,unit) result(nfailed)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftbox_utests'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: fftalg,ndat,nthreads
 integer,optional,intent(in) :: unit
 integer :: nfailed
!arrays

!Local variables-------------------------------
!scalars
 integer,parameter :: NSETS=6
 integer :: ifft,ierr,ldxyz,old_nthreads,ount,cplex
 integer :: iset,nx,ny,nz,ldx,ldy,ldz,fftalga,fftalgc
 !integer :: ix,iy,iz,padat,dat
 real(dp),parameter :: ATOL_SP=tol6,ATOL_DP=tol12 ! Tolerances on the absolute error
 real(dp) :: max_abserr
 character(len=500) :: msg,info,library,cplex_mode,padding_mode
 type(fftbox_plan3_t) :: bw_plan,fw_plan
!arrays
 integer :: pars(6,NSETS)
 real(dp) :: crand(2)
 real(dp),allocatable :: fofg(:),fofr_ref(:),fofr(:)
 complex(dpc),allocatable :: ff(:),ff_ref(:),gg(:)
 complex(spc),allocatable :: ffsp(:),ff_refsp(:),ggsp(:)

! *************************************************************************

 nfailed = 0

 ount = std_out; if (PRESENT(unit)) ount = unit

 if (nthreads>0) then
   old_nthreads = xomp_get_max_threads()
   call xomp_set_num_threads(nthreads)
 end if

 ! These values must be compatible with all the FFT routines.
 ! SG library is the most restrictive (only powers of 2,3,5).
 pars = RESHAPE( (/   &
&  12,18,15,12,18,15, &
&  12,18,15,13,19,16, &
&  12,18,15,13,19,15, &
&  12,18,15,12,18,16, &
&  12,18,15,13,18,15, &
&  12,18,15,15,21,18  &
& /), (/6,NSETS/) )

 fftalga=fftalg/100; fftalgc=mod(fftalg,10)

 call fftalg_info(fftalg,library,cplex_mode,padding_mode)

 do iset=1,SIZE(pars,DIM=2)
   nx =pars(1,iset);  ny=pars(2,iset);  nz=pars(3,iset)
   ldx=pars(4,iset); ldy=pars(5,iset); ldz=pars(6,iset)
   !write(std_out,*)pars(1:6,iset)

   ! Create the FFT plans.
   call fftbox_plan3_many(bw_plan,ndat,pars(1,iset),pars(4,iset),fftalg,+1)
   call fftbox_plan3_many(fw_plan,ndat,pars(1,iset),pars(4,iset),fftalg,-1)

   ldxyz = ldx*ldy*ldz
   !
   ! ================================================================
   ! === TEST the single precision version of dfti_c2c_* routines ===
   ! ================================================================
   ABI_MALLOC(ff_refsp, (ldxyz*ndat))
   ABI_MALLOC(ffsp,     (ldxyz*ndat))
   ABI_MALLOC(ggsp,     (ldxyz*ndat))

   do ifft=1,ldxyz*ndat
     call RANDOM_NUMBER(crand)
     ff_refsp(ifft) = DCMPLX(crand(1), crand(2))
   end do

   ! Set the augmentation region to zero, because FFTW3 wrappers
   ! use zscal to scale the results.
   call cplx_setaug_zero_spc(nx,ny,nz,ldx,ldy,ldz,ndat,ff_refsp)
   ffsp = ff_refsp

   ! in-place version.
   call fftbox_execute(bw_plan,ffsp)
   call fftbox_execute(fw_plan,ffsp)

   ierr = COUNT(ABS(ffsp - ff_refsp) > ATOL_SP)
   nfailed = nfailed + ierr

   info = sjoin(library,"c2c_ip_spc :")
   if (ierr /= 0) then
     max_abserr = MAXVAL(ABS(ffsp - ff_refsp))
     write(msg,"(a,es9.2,a)")" FAILED (max_abserr = ",max_abserr,")"
   else
     write(msg,"(a)")" OK"
   end if
   call wrtout(ount,sjoin(info,msg),"COLL")

   ffsp = ff_refsp
   call fftbox_execute(bw_plan,ffsp,ggsp)
   call fftbox_execute(fw_plan,ggsp,ffsp)

   ierr = COUNT(ABS(ffsp - ff_refsp) > ATOL_SP)
   nfailed = nfailed + ierr

   info = sjoin(library,"c2c_op_spc :")
   if (ierr /= 0) then
     max_abserr = MAXVAL(ABS(ffsp - ff_refsp))
     write(msg,"(a,es9.2,a)")" FAILED (max_abserr = ",max_abserr,")"
   else
     write(msg,"(a)")" OK"
   end if
   call wrtout(ount,sjoin(info,msg),"COLL")

   ABI_FREE(ff_refsp)
   ABI_FREE(ffsp)
   ABI_FREE(ggsp)

   ! ================================================================
   ! === TEST the double precision version of dfti_c2c_* routines ===
   ! ================================================================
   ABI_MALLOC(ff_ref, (ldxyz*ndat))
   ABI_MALLOC(ff,     (ldxyz*ndat))
   ABI_MALLOC(gg,     (ldxyz*ndat))

   do ifft=1,ldxyz*ndat
     call RANDOM_NUMBER(crand)
     ff_ref(ifft) = DCMPLX(crand(1), crand(2))
   end do

   ! Set the augmentation region to zero, because FFTW3 wrappers
   ! use zscal to scale the results.
   call cplx_setaug_zero_dpc(nx,ny,nz,ldx,ldy,ldz,ndat,ff_ref)
   ff = ff_ref

   call fftbox_execute(bw_plan,ff)
   call fftbox_execute(fw_plan,ff)

   ierr = COUNT(ABS(ff - ff_ref) > ATOL_DP)
   nfailed = nfailed + ierr

   info = sjoin(library,"c2c_ip_dpc :")
   if (ierr /= 0) then
     max_abserr = MAXVAL(ABS(ff - ff_ref))
     write(msg,"(a,es9.2,a)")" FAILED (max_abserr = ",max_abserr,")"
   else
     write(msg,"(a)")" OK"
   end if
   call wrtout(ount,sjoin(info,msg),"COLL")

   ff = ff_ref
   call fftbox_execute(bw_plan,ff,gg)
   call fftbox_execute(fw_plan,gg,ff)

   ierr = COUNT(ABS(ff - ff_ref) > ATOL_DP)
   nfailed = nfailed + ierr

   info = sjoin(library,"c2c_op_dpc :")
   if (ierr /= 0) then
     max_abserr = MAXVAL(ABS(ff - ff_ref))
     write(msg,"(a,es9.2,a)")" FAILED (max_abserr = ",max_abserr,")"
   else
     write(msg,"(a)")" OK"
   end if
   call wrtout(ount,sjoin(info,msg),"COLL")

   ABI_FREE(ff_ref)
   ABI_FREE(ff)
   ABI_FREE(gg)

   do cplex=1,2
     !
     !if (fftalga == FFT_FFTW3 .and. ndat > 1 .and. cplex==1) then
     !  call wrtout(ount,"Warning: fourdp with FFTW3-wrappers, cplex=2 and ndat>1, might crash if MKL is used","COLL")
     !  !CYCLE
     !end if

     ABI_MALLOC(fofg,     (2*ldxyz*ndat))
     ABI_MALLOC(fofr_ref, (cplex*ldxyz*ndat))
     ABI_MALLOC(fofr,     (cplex*ldxyz*ndat))

     call RANDOM_NUMBER(fofr_ref)
     !call cg_setaug_zero(cplex,nx,ny,nz,ldx,ldy,ldz,ndat,fofr_ref)
     fofr = fofr_ref

     SELECT CASE (fftalga)
     CASE (FFT_FFTW3)
       call fftw3_seqfourdp(cplex,nx,ny,nz,ldx,ldy,ldz,ndat,-1,fofg,fofr)
       call fftw3_seqfourdp(cplex,nx,ny,nz,ldx,ldy,ldz,ndat,+1,fofg,fofr)

     CASE (FFT_DFTI)
       call dfti_seqfourdp(cplex,nx,ny,nz,ldx,ldy,ldz,ndat,-1,fofg,fofr)
       call dfti_seqfourdp(cplex,nx,ny,nz,ldx,ldy,ldz,ndat,+1,fofg,fofr)

     CASE DEFAULT
       ! TODO
       continue
     END SELECT

     !call cg_setaug_zero(cplex,nx,ny,nz,ldx,ldy,ldz,ndat,fofr)

     ierr = COUNT(ABS(fofr - fofr_ref) > ATOL_DP)
     nfailed = nfailed + ierr

     write(info,"(a,i1,a)")sjoin(library,"fourdp (cplex "),cplex,") :"
     if (ierr /= 0) then
       max_abserr = MAXVAL(ABS(fofr - fofr_ref))
       write(msg,"(a,es9.2,a)")" FAILED (max_abserr = ",max_abserr,")"
     else
       write(msg,"(a)")" OK"
     end if
     call wrtout(ount,sjoin(info,msg),"COLL")

     ABI_FREE(fofg)
     ABI_FREE(fofr_ref)
     ABI_FREE(fofr)
   end do
   !
 end do

 if (nthreads > 0) then
   call xomp_set_num_threads(old_nthreads)
 end if

end function fftbox_utests
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftu_utests
!! NAME
!! fftu_utests
!!
!! FUNCTION
!! Unit tests for the FFTs of wavefunctions (sequential version).
!!
!! INPUTS
!!
!! OUTPUT
!!  nfailed=number of failed tests.
!!
!! PARENTS
!!      fftprof
!!
!! SOURCE

function fftu_utests(ecut,ngfft,rprimd,ndat,nthreads,unit) result(nfailed)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftu_utests'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_51_manage_mpi
 use interfaces_52_fft_mpi_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndat,nthreads
 integer :: nfailed
 integer,optional,intent(in) :: unit
 real(dp),intent(in) :: ecut
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer,parameter :: nspinor1=1,mkmem1=1,exchn2n3d0=0,ikg0=0
 integer :: nx,ny,nz,nxyz,ldx,ldy,ldz,ierr,npw_k,mgfft,istwf_k,ikpt,ldxyz,ipw,old_nthreads,ount
 integer :: fftalg,npw_k_test
 real(dp),parameter :: ATOL_SP=tol6, ATOL_DP=tol12 ! Tolerances on the absolute error
 real(dp) :: max_abserr,ucvol
 character(len=500) :: msg,info,library,cplex_mode,padding_mode
!arrays
 integer :: kg_dum(3,0)
 integer,allocatable :: gbound_k(:,:),kg_k(:,:)
 real(dp) :: kpoint(3),crand(2),kpoints(3,9)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: cg(:,:),cg_ref(:,:),cr(:,:)
 complex(spc),allocatable :: ugsp(:),ug_refsp(:),ursp(:)
 complex(dpc),allocatable :: ug(:),ug_ref(:),ur(:)
 type(MPI_type) :: MPI_enreg_seq

! *************************************************************************

 ount = std_out; if (PRESENT(unit)) ount = unit

 nfailed = 0
 fftalg = ngfft(7)

 if (nthreads > 0) then
   old_nthreads = xomp_get_max_threads()
   call xomp_set_num_threads(nthreads)
 end if

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 nx  = ngfft(1);  ny = ngfft(2);  nz = ngfft(3)
 ldx = ngfft(4); ldy = ngfft(5); ldz = ngfft(6)
 mgfft = MAXVAL(ngfft(1:3))

  nxyz =  nx* ny* nz
 ldxyz = ldx*ldy*ldz

 ABI_CALLOC(cg_ref, (2, ldxyz*ndat))
 ABI_CALLOC(cg,     (2, ldxyz*ndat))
 ABI_CALLOC(cr,     (2, ldxyz*ndat))

 ABI_CALLOC(ug_ref, (ldxyz*ndat))
 ABI_CALLOC(ug,     (ldxyz*ndat))
 ABI_CALLOC(ur,     (ldxyz*ndat))

 ABI_CALLOC(ug_refsp, (ldxyz*ndat))
 ABI_CALLOC(ugsp,     (ldxyz*ndat))
 ABI_CALLOC(ursp,     (ldxyz*ndat))

 kpoints = RESHAPE( (/ &
&  0.1, 0.2, 0.3, &
&  0.0, 0.0, 0.0, &
&  0.5, 0.0, 0.0, &
&  0.0, 0.0, 0.5, &
&  0.5, 0.0, 0.5, &
&  0.0, 0.5, 0.0, &
&  0.5, 0.5, 0.0, &
&  0.0, 0.5, 0.5, &
&  0.5, 0.5, 0.5 /), (/3,9/) )

 call fftalg_info(fftalg,library,cplex_mode,padding_mode)

 call initmpi_seq(MPI_enreg_seq)

 do ikpt=1,SIZE(kpoints,DIM=2)
   kpoint = kpoints(:,ikpt)

   istwf_k = set_istwfk(kpoint)

   ! * Calculate the number of G-vectors for this k-point.
   call kpgsph(ecut,exchn2n3d0,gmet,ikg0,0,istwf_k,kg_dum,kpoint,0,MPI_enreg_seq,0,npw_k)
   !
   ! * Allocate and calculate the set of G-vectors.
   ABI_MALLOC(kg_k,(3,npw_k))
   call kpgsph(ecut,exchn2n3d0,gmet,ikg0,0,istwf_k,kg_k,kpoint,mkmem1,MPI_enreg_seq,npw_k,npw_k_test)

   ABI_MALLOC(gbound_k,(2*mgfft+8,2))
   call sphereboundary(gbound_k,istwf_k,kg_k,mgfft,npw_k)

   !if (istwf_k==2) then
   !  do ipw=1,npw_k
   !    write(std_out,*)ipw, kg_k(:,ipw)
   !  end do
   !  stop
   !end if

#if 0
   !TODO
   ! ================================================
   ! === Test the double precision 2*real version ===
   ! ================================================
   do ipw=1,npw_k*ndat
     call RANDOM_NUMBER(crand)
     cg_ref(:,ipw) = crand(:)
   end do

   if (istwf_k == 2) then
     do ipw=1,npw_k*ndat,npw_k
       cg_ref(2,ipw) = zero
     end do
   end if

   cg = cg_ref

   call fft_ug_dp(npw_k,nxyz,nspinor1,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,cg,cr)
   call fft_ur_dp(npw_k,nxyz,nspinor1,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,cr,cg)

   ierr = COUNT(ABS(cg - cg_ref) > ATOL_DP)
   nfailed = nfailed + ierr

   write(info,"(a,i1,a)")sjoin(library,"fftu_dp, istwfk "),istwf_k," :"
   if (ierr /= 0) then
     max_abserr = MAXVAL(ABS(cg - cg_ref))
     write(msg,"(a,es9.2,a)")" FAILED (max_abserr = ",max_abserr,")"
   else
     write(msg,"(a)")" OK"
   end if
   call wrtout(ount,sjoin(info,msg),"COLL")
#endif

   ! =================================================
   ! === Test the single precision complex version ===
   ! =================================================
   do ipw=1,npw_k*ndat
     call RANDOM_NUMBER(crand)
     ug_refsp(ipw) = CMPLX(crand(1), crand(2))
   end do

   if (istwf_k == 2) then
     do ipw=1,npw_k*ndat,npw_k
       ug_refsp(ipw) = REAL(ug_refsp(ipw))
     end do
   end if

   ugsp = ug_refsp

   call fft_ug(npw_k,nxyz,nspinor1,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,ugsp,ursp)
   call fft_ur(npw_k,nxyz,nspinor1,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,ursp,ugsp)

   ierr = COUNT(ABS(ugsp - ug_refsp) > ATOL_SP)
   nfailed = nfailed + ierr

   write(info,"(a,i1,a)")sjoin(library,"fftu_spc, istwfk "),istwf_k," :"
   if (ierr /= 0) then
     max_abserr = MAXVAL(ABS(ugsp - ug_refsp))
     write(msg,"(a,es9.2,a)")" FAILED (max_abserr = ",max_abserr,")"
   else
     write(msg,"(a)")" OK"
   end if
   call wrtout(ount,sjoin(info,msg),"COLL")

   ! =================================================
   ! === Test the double precision complex version ===
   ! =================================================
   do ipw=1,npw_k*ndat
     call RANDOM_NUMBER(crand)
     ug_ref(ipw) = DCMPLX(crand(1), crand(2))
   end do

   if (istwf_k == 2) then
     do ipw=1,npw_k*ndat,npw_k
       ug_ref(ipw) = REAL(ug_ref(ipw))
     end do
   end if

   ug = ug_ref

   call fft_ug(npw_k,nxyz,nspinor1,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,ug,ur)
   call fft_ur(npw_k,nxyz,nspinor1,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,ur,ug)

   ierr = COUNT(ABS(ug - ug_ref) > ATOL_DP)
   nfailed = nfailed + ierr

   write(info,"(a,i1,a)")sjoin(library,"fftu_dpc, istwfk "),istwf_k," :"
   if (ierr /= 0) then
     max_abserr = MAXVAL(ABS(ug - ug_ref))
     write(msg,"(a,es9.2,a)")" FAILED (max_abserr = ",max_abserr,")"
   else
     write(msg,"(a)")" OK"
   end if
   call wrtout(ount,sjoin(info,msg),"COLL")

   ABI_FREE(kg_k)
   ABI_FREE(gbound_k)
 end do

 ABI_FREE(cg_ref)
 ABI_FREE(cg)
 ABI_FREE(cr)

 ABI_FREE(ug_ref)
 ABI_FREE(ug)
 ABI_FREE(ur)

 ABI_FREE(ug_refsp)
 ABI_FREE(ugsp)
 ABI_FREE(ursp)

 call destroy_mpi_enreg(MPI_enreg_seq)

 if (nthreads > 0) then
   call xomp_set_num_threads(old_nthreads)
 end if

end function fftu_utests
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftbox_mpi_utests
!! NAME
!!  fftbox_mpi_utests
!!
!! FUNCTION
!! Driver routine for unit tests of the MPI FFT routines.
!!
!! INPUTS
!! fftalg =fftalg input variable.
!! cplex=1 for r2c, 2 for c2c transforms.
!! ndat = Number of transform to execute
!! nthreads = Number of OpenMP threads.
!! comm_fft=MPI communicator for the FFT
!! [unit]=Output Unit number (DEFAULT std_out)
!!
!! OUTPUT
!!  nfailed=number of failures.
!!
!! PARENTS
!!      fftprof
!!
!! SOURCE

function fftbox_mpi_utests(fftalg,cplex,ndat,nthreads,comm_fft,unit) result(nfailed)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftbox_mpi_utests'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: fftalg,cplex,ndat,nthreads,comm_fft
 integer,optional,intent(in) :: unit
 integer :: nfailed

!Local variables-------------------------------
!scalars
 integer,parameter :: NSETS=6
 integer :: ierr,old_nthreads,ount,iset,mpierr,nfft,me_fft
 integer :: nproc_fft,fftalga,fftalgc,n1,n2,n3,n4,n5,n6
 real(dp),parameter :: ATOL_DP=tol12
 real(dp) :: max_abserr
 real(dp) ::  ctime,wtime,gflops
 character(len=500) :: msg,info,library,cplex_mode,padding_mode
 type(distribfft_type),target :: fftabs
!arrays
 integer :: pars(6,NSETS),ngfft(18)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp),allocatable :: fofg(:,:),fofr(:),fofr_copy(:)

! *************************************************************************

 ount = std_out; if (PRESENT(unit)) ount = unit
 nfailed = 0

 if (nthreads>0) then
   old_nthreads = xomp_get_max_threads()
   call xomp_set_num_threads(nthreads)
 end if

 ! These values must be compatible with all the FFT routines.
 ! SG library is the most restrictive (only powers of 2,3,5).
 pars = RESHAPE( [    &
&  12,18,15,12,18,15, &
&  12,18,15,13,19,16, &
&  12,18,15,13,19,15, &
&  12,18,15,12,18,16, &
&  12,18,15,13,18,15, &
&  12,18,15,15,21,18  &
& ], [6,NSETS] )

 fftalga=fftalg/100; fftalgc=mod(fftalg,10)

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

 call fftalg_info(fftalg,library,cplex_mode,padding_mode)

 !do iset=1,SIZE(pars,DIM=2)
 do iset=1,1
   n1=pars(1,iset); n2=pars(2,iset); n3=pars(3,iset)
   ! For the time being, ngfft(2) and ngfft(3) must be multiple of nproc_fft
   n2 = n2 * nproc_fft; n3 = n3 * nproc_fft
   !n4=pars(4,iset); n5=pars(5,iset); n6=pars(6,iset)
   n4=n1; n5=n2; n6=n3

   !n4=n1+1; n5=n2; n6=n3
   !n4=n1; n5=n2+1; n6=n3
   !n4=n1; n5=n2; n6=n3+1
   !write(std_out,*)pars(1:6,iset)
   !cplex = 1

   !call getng(boxcutmin,ecut,gmet,kpt,me_fft,mgfft,nfft,ngfft,nproc_fft,nsym,paral_fft,symrel,&
   !&ngfftc,use_gpu_cuda,unit) ! optional

   ! Init ngfft
   ! TODO Propagate more info via ngfft, define helper functions, write routine to get fftcache
   ngfft = 0
   ngfft(1:7) = [n1,n2,n3,n4,n5,n6,fftalg]
   ngfft(8)= get_cache_kb()        ! cache
   ngfft(9)=1                      ! paral_fft_
   ngfft(10)=nproc_fft             ! nproc_fft
   ngfft(11)=xmpi_comm_rank(comm_fft)  ! me_fft
   ngfft(12)=ngfft(2)/nproc_fft    ! n2proc
   ngfft(13)=ngfft(3)/nproc_fft    ! n3proc

   !call print_ngfft(ngfft,"ngfft for MPI-fourdp",unit=std_out,mode_paral="COLL",prtvol=0)

   ! Allocate arrays, fill fofr with random numbers and keep a copy.
   nfft = (n1 * n2 * n3) / nproc_fft
   ABI_MALLOC(fofg, (2,nfft*ndat))
   ABI_MALLOC(fofr, (cplex*nfft*ndat))
   ABI_MALLOC(fofr_copy, (cplex*nfft*ndat))

   call RANDOM_NUMBER(fofr)
   fofr_copy = fofr

   call init_distribfft(fftabs,"c",nproc_fft,n2,n3)
   fftn2_distrib => fftabs%tab_fftdp2_distrib
   ffti2_local => fftabs%tab_fftdp2_local
   fftn3_distrib => fftabs%tab_fftdp3_distrib
   ffti3_local => fftabs%tab_fftdp3_local

   call cwtime(ctime,wtime,gflops,"start")

   select case (fftalga)

   case (FFT_SG2002)
     call sg2002_mpifourdp(cplex,nfft,ngfft,ndat,-1,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)
     call sg2002_mpifourdp(cplex,nfft,ngfft,ndat,+1,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)

   case (FFT_FFTW3)
     call fftw3_mpifourdp(cplex,nfft,ngfft,ndat,-1,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)
     call fftw3_mpifourdp(cplex,nfft,ngfft,ndat,+1,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)

   !case (FFT_DFTI)
   !  call dfti_mpifourdp(cplex,nfft,ngfft,ndat,-1,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)
   !  call dfti_mpifourdp(cplex,nfft,ngfft,ndat,+1,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)

   case default
     write(msg,'(a,i0,a)')"fftalg: ",fftalg," does not support MPI-FFT"
     MSG_BUG(msg)
   end select

   call cwtime(ctime,wtime,gflops,"stop")
   if (me_fft == 0) then
     !write(std_out,'(a,i0,2(a,f10.4))')"fftalg: ",fftalg,", cpu_time: ",ctime,", wall_time: ",wtime
   end if

   ! Check if F^{-1} F = I within ATOL_DP
   ierr = COUNT(ABS(fofr - fofr_copy) > ATOL_DP)
   call xmpi_sum(ierr,comm_fft,mpierr)
   nfailed = nfailed + ierr

   if (cplex == 1) info = sjoin(library,"r2c --> c2r :")
   if (cplex == 2) info = sjoin(library,"c2c :")

   if (ierr /= 0) then
     ! Compute the maximum of the absolute error.
     max_abserr = MAXVAL(ABS(fofr - fofr_copy))
     call xmpi_max(max_abserr,comm_fft,mpierr)
     write(msg,"(a,es9.2,a)")" FAILED (max_abserr = ",max_abserr,")"
   else
     write(msg,"(a)")" OK"
   end if
   call wrtout(ount,sjoin(info,msg),"COLL")

   call destroy_distribfft(fftabs)

   ABI_FREE(fofg)
   ABI_FREE(fofr_copy)
   ABI_FREE(fofr)
 end do

 if (nthreads > 0) then
   call xomp_set_num_threads(old_nthreads)
 end if

end function fftbox_mpi_utests
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftu_mpi_utests
!! NAME
!! fftu_mpi_utests
!!
!! FUNCTION
!! Unit tests for the FFTs of wavefunctions (MPI vesion).
!!
!! INPUTS
!!
!! OUTPUT
!!  nfailed=number of failed tests.
!!
!! PARENTS
!!      fftprof
!!
!! SOURCE

function fftu_mpi_utests(fftalg,ecut,rprimd,ndat,nthreads,comm_fft,paral_kgb,unit) result(nfailed)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftu_mpi_utests'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_52_fft_mpi_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fftalg,ndat,nthreads,comm_fft,paral_kgb
 integer :: nfailed
 integer,optional,intent(in) :: unit
 real(dp),intent(in) :: ecut
!arrays
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer,parameter :: nsym1=1,npw0=0,cplex_one=1,istwfk_one=1
 integer :: n1,n2,n3,idat,n4,n5,n6,ierr,npw_k,full_npw_k,istwf_npw_k,cplexwf
 integer :: mgfft,istwf_k,ikpt,old_nthreads,ount,isign,fftalga,fftalgc
 integer :: ig,i1,i2,i3,i3_glob,i3dat,nd3proc,i3_local,g0sender
 integer :: step,me_g0,me_fft,nproc_fft,mpierr,nfft,cplex
 real(dp),parameter :: boxcutmin2=two,ATOL_DP=tol12,RTOL_DP=tol3 ! Tolerances on the absolute and relative error
 real(dp),parameter :: weight_r=one,weight_i=one
 real(dp) :: max_abserr,max_relerr,ucvol,relerr,den,refden
 real(dp) ::  ctime,wtime,gflops
 character(len=500) :: msg,info,library,cplex_mode,padding_mode
 type(distribfft_type) :: fftabs
!arrays
 integer :: symrel(3,3,nsym1),ngfft(18)
 integer,allocatable :: full_kg_k(:,:),istw_kg_k(:,:)
 integer,allocatable :: gbound_k(:,:),kg_k(:,:)
 real(dp) :: dummy_fofg(0,0) !dummy_denpot(0,0,0)
 real(dp) :: kpoint(3),kpoints(3,2)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),allocatable :: fofg(:,:),ref_fofg(:,:),fofg_out(:,:),fofr(:,:,:,:)
 real(dp),allocatable :: density(:,:,:),pot(:,:,:),invpot(:,:,:)
 real(dp),allocatable :: full_fofg(:,:),istwf_fofg(:,:)

! *************************************************************************

 nfailed = 0

 ount = std_out; if (PRESENT(unit)) ount = unit

 if (nthreads > 0) then
   old_nthreads = xomp_get_max_threads()
   call xomp_set_num_threads(nthreads)
 end if

 if (.not. fftalg_has_mpi(fftalg)) then
   write(msg,'(a,i0,a)')"fftalg: ",fftalg," does not support MPI"
   MSG_ERROR(msg)
 end if

 nproc_fft = xmpi_comm_size(comm_fft); me_fft = xmpi_comm_rank(comm_fft)

 symrel = reshape([1,0,0,0,1,0,0,0,1],[3,3,nsym1])
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 kpoints = RESHAPE( [ &
&  0.1, 0.2, 0.3, &
&  0.0, 0.0, 0.0 ], [3,2] )

 call fftalg_info(fftalg,library,cplex_mode,padding_mode)

 fftalga=fftalg/100; fftalgc=mod(fftalg,10)

 do ikpt=1,size(kpoints, dim=2)
   kpoint = kpoints(:,ikpt)

   ! Get full G-sphere (no MPI distribution, no time-reversal)
   call get_kg(kpoint,istwfk_one,ecut,gmet,full_npw_k,full_kg_k)

   ! Get full G-sphere (no MPI distribution, use time-reversal if possible)
   istwf_k = set_istwfk(kpoint)
   call get_kg(kpoint,istwf_k,ecut,gmet,istwf_npw_k,istw_kg_k)

   ! Get FFT box dims.
   ngfft = -1
   ngfft(7) = fftalg
   ngfft(8) = get_cache_kb()

   call getng(boxcutmin2,ecut,gmet,kpoint,me_fft,mgfft,nfft,ngfft,nproc_fft,nsym1,&
&             paral_kgb,symrel,unit=dev_null)

   n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
   ! Do not use augmentation.
   !ngfft(4:6) = ngfft(1:3)
   n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

   call print_ngfft(ngfft,"ngfft for MPI-fourwf",unit=std_out,mode_paral="COLL",prtvol=0)

   ! Compute FFT distribution tables.
   call init_distribfft(fftabs,"c",nproc_fft,n2,n3)

   ! Set to 1 if this node owns G = 0.
   me_g0 = 0; if (fftabs%tab_fftwf2_distrib(1) == me_fft) me_g0 = 1
   g0sender = fftabs%tab_fftwf2_distrib(1)

   ! Allocate u(g) on the full sphere, initialize with random numbers.
   ! and broadcast full data to the other nodes.
   ABI_MALLOC(full_fofg, (2,full_npw_k*ndat))

   if (me_fft == g0sender) then
     if (istwf_k == 1) then
       call RANDOM_NUMBER(full_fofg)
       !full_fofg = one

     else if (istwf_k == 2) then
       ! Special treatment for real wavefunctions.
       ! Fill the irreducible G-vectors first so that we have the u(g) of a real u(r)
       ! Then get u(g) on the full G-sphere.
       ! TODO: Sequential version OK, SIGSEV if mpi_ncpus > 1!
       ABI_MALLOC(istwf_fofg, (2,istwf_npw_k*ndat))
       call RANDOM_NUMBER(istwf_fofg)
       ! Enforce real u in G-space.
       do idat=1,ndat
         istwf_fofg(2,1+(idat-1)*istwf_npw_k) = zero
       end do

       ! from istwfk 2 to 1.
       call change_istwfk(istwf_npw_k,istw_kg_k,istwf_k,full_npw_k,full_kg_k,istwfk_one,n1,n2,n3,ndat,istwf_fofg,full_fofg)
       ABI_FREE(istwf_fofg)

     else
       MSG_ERROR("istwf_k /= [1,2] not available in MPI-FFT mode")
     end if
   end if

   call xmpi_bcast(full_fofg,g0sender,comm_fft,ierr)

   ! Compute sphere boundaries for zero-padded FFTs
   ABI_MALLOC(gbound_k,(2*mgfft+8,2))
   call sphereboundary(gbound_k,istwfk_one,full_kg_k,mgfft,full_npw_k)

   ! Extract my G-vectors from full_kg_k and store them in kg_k.
   !write(std_out,*)"fftwf2_distrib",fftabs%tab_fftwf2_distrib
   do step=1,2
     if (step == 2) then
       ! Allocate my u(g) and my Gs.
       ! Set fofg to zero to bypass a bug with XLF!
       ABI_MALLOC(kg_k, (3, npw_k))
       ABI_CALLOC(fofg, (2,npw_k*ndat))
     end if

     npw_k = 0
     do ig=1,full_npw_k
       i1=full_kg_k(1,ig); if(i1<0) i1=i1+n1; i1=i1+1
       i2=full_kg_k(2,ig); if(i2<0) i2=i2+n2; i2=i2+1
       i3=full_kg_k(3,ig); if(i3<0) i3=i3+n3; i3=i3+1
       if (fftabs%tab_fftwf2_distrib(i2) == me_fft) then
         npw_k = npw_k + 1
         if (step == 2) then
           kg_k(:,npw_k) = full_kg_k(:,ig)
           fofg(:,npw_k) = full_fofg(:,ig)
         end if
       end if
     end do
   end do ! step

   ABI_FREE(istw_kg_k)

   !write(std_out,*)"dist",trim(itoa(me_fft)),fftabs%tab_fftdp3_distrib(:) !== me_fft)
   !write(std_out,*)"local",trim(itoa(me_fft)),fftabs%tab_fftdp3_local(:)

   ! Allocate my local portion of u(r), and keep a copy my u(g).
   ABI_MALLOC(fofr, (2,n4,n5,n6*ndat))
   ABI_MALLOC(ref_fofg, (2,npw_k*ndat))
   ref_fofg = fofg

   call cwtime(ctime,wtime,gflops,"start")

   ! ----------------------------------------------------------------
   ! Test the the backward and the forward transform of wavefunctions
   ! ----------------------------------------------------------------
   ! c2c or [c2r, r2c]
   cplexwf = 2; if (istwf_k==2) cplexwf = 1
   ! FIXME:
   ! There's a bug somewhere in the MPI routines if cplexwf == 1 is used and ndat > 1
   ! Not a serious problem at present since cplexwf==1 is never used in abinit.
   if (ndat>1) cplexwf = 2
   if (nproc_fft > 1) cplexwf = 2
   cplexwf = 2
   !cplexwf = 2; if (istwf_k==2) cplexwf = 1

   do isign = 1,-1,-2
     call fftmpi_u(npw_k,n4,n5,n6,ndat,mgfft,ngfft,&
&      istwfk_one,gbound_k,kg_k,me_g0,fftabs,isign,fofg,fofr,comm_fft,cplexwf=cplexwf)
   end do

   ! The final interface should be:
   !subroutine fftmpi_u(npw_k,n4,n5,n6,nspinor,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,fftabs,isign,fofg,fofr)

   ! Parallel version (must support augmentation in forf(:,:,:,:), hence we have a different API wrt the seq case!
   !call fftmpi_ug(npw_k,n4,n5,n6,nspinor,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,fftabs,fofg,fofr)
   !call fftmpi_ur(npw_k,n4,n5,n6,nspinor,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,fftabs,fofr,fofgout)

   ! Seq version.
   !call fft_ug(npw_k,nxyz,nspinor,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,ugsp,ursp)
   !call fft_ur(npw_k,nxyz,nspinor,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,ursp,ugsp)

   call cwtime(ctime,wtime,gflops,"stop")
   if (me_fft == 0) then
     !write(std_out,'(a,i0,3(a,f10.4))')"fftalg: ",fftalg,", ecut: ",ecut,", cpu_time: ",ctime,", wall_time: ",wtime
   end if

   ! Check if F^{-1} F = I within ATOL_DP
   ierr = COUNT(ABS(fofg - ref_fofg) > ATOL_DP)
   call xmpi_sum(ierr,comm_fft,mpierr)
   nfailed = nfailed + ierr

   write(info,"(a,i1,a)")sjoin(library,"fftu_mpi, istwfk "),istwf_k," :"

   if (ierr /= 0) then
     max_abserr = MAXVAL(ABS(fofg - ref_fofg))
     call xmpi_max(max_abserr,comm_fft,mpierr)
     write(msg,"(a,es9.2,a)")" FAILED (max_abserr = ",max_abserr,")"
   else
     write(msg,"(a)")" OK"
   end if
   call wrtout(ount,sjoin(info,msg),"COLL")

   ! -------------------------------------------
   ! Test the accumulation of density (option 1)
   ! -------------------------------------------
   ABI_CALLOC(density, (cplex_one*n4,n5,n6))
   !fofg = one

   ! Accumulate density. Does not work if cplexwf==1
   call fourwf_mpi(cplex_one,density,fofg,dummy_fofg,fofr,&
&    gbound_k,gbound_k,istwfk_one,kg_k,kg_k,me_g0,mgfft,ngfft,fftabs,n1,n2,n3,&
&    npw_k,npw_k,n4,n5,n6,ndat,1,weight_r,weight_i,comm_fft,cplexwf=cplexwf)

!   Recompute u(r)
!   call fourwf_mpi(cplex_one,density,fofg,dummy_fofg,fofr,&
!&    gbound_k,gbound_k,istwf_one,kg_k,kg_k,me_g0,mgfft,ngfft,fftabs,n1,n2,n3,&
!&    npw_k,npw_k,n4,n5,n6,ndat,0,weight_r,weight_i,comm_fft,cplexwf=cplexwf)
!   call xmpi_sum(density,comm_fft,ierr)
!   if (me_fft == 0) write(std_out,*)"sum density", sum(density)

   !do i3=1,n6
   !  write(110+me_fft,*)i3,density(:,:,i3)
   !end do

   max_relerr = zero

   ! FIXME: This is wrong if ndat > 1
   ! Must unify the treatment of fofr and rhor on the augmented mesh (n4,n5,n6)
   ! back_wf and for_wf return a MPI distributed arrays but fofr is allocated with
   ! the global dimensions.
   nd3proc=(n3-1)/nproc_fft+1
   do i3=1,n3
   !do i3=1,nd3proc
     if( me_fft == fftabs%tab_fftdp3_distrib(i3) ) then
       i3_local = fftabs%tab_fftdp3_local(i3) !+ nd3proc*(idat-1)
       !i3_local = i3
       i3_glob = i3
       !i3_glob = i3_local
       !i3_glob = i3 + nd3proc * me_fft
       !i3dat = i3  + n6 * (idat-1)
       do i2=1,n2
         do i1=1,n1
            den = density(i1,i2,i3_glob)
            refden = zero
            do idat=1,ndat
              i3dat = i3_local + n6 * (idat-1)
              refden = refden + weight_r*fofr(1,i1,i2,i3dat)**2+ weight_i*fofr(2,i1,i2,i3dat)**2
            end do
            relerr = abs(refden - den)
            if (abs(refden) < tol12) refden = tol12
            relerr = relerr / abs(refden)
            max_relerr = max(max_relerr, relerr)
            !if (relerr > RTOL_DP) write(std_out,*)trim(itoa(me_fft)),i1,i2,i3,idat,den,refden
         end do
       end do
     end if
   end do

   call xmpi_max(max_relerr,comm_fft,mpierr)

   write(info,"(a,i1,a)")sjoin(library,"accrho_mpi, istwfk "),istwf_k," :"
   if (max_relerr > RTOL_DP) then
     write(msg,"(a,es9.2,a)")" FAILED (max_relerr = ",max_relerr,")"
   else
     write(msg,"(a)")" OK"
   end if
   call wrtout(ount,sjoin(info,msg),"COLL")

   ABI_FREE(density)

   ! -------------------------------------------
   ! Test the application of the local potential
   ! -------------------------------------------
   cplex = 1
   ABI_MALLOC(pot, (cplex*n4,n5,n6))
   ABI_MALLOC(invpot, (cplex*n4,n5,n6))

   if (me_fft == g0sender) then
     !pot = fofr(1,:,:,:)
     !call RANDOM_NUMBER(pot)
     !where (abs(pot) < tol4) pot = tol4
     ! Simplest potential ever (G=0 to reduce errors due to aliasing)
     pot = four
     invpot = one/pot
   end if

   call xmpi_bcast(pot,g0sender,comm_fft,ierr)
   call xmpi_bcast(invpot,g0sender,comm_fft,ierr)

   ABI_MALLOC(fofg_out, (2,npw_k*ndat))

   ! Compute fofg_out = <G|pot(r)|fofg>
   call fourwf_mpi(cplex,pot,fofg,fofg_out,fofr,&
&    gbound_k,gbound_k,istwfk_one,kg_k,kg_k,me_g0,mgfft,ngfft,fftabs,n1,n2,n3,&
&    npw_k,npw_k,n4,n5,n6,ndat,2,weight_r,weight_i,comm_fft,cplexwf=cplexwf)

   ! Compute fofg = <G|1/pot(r)|fofg_out>
   call fourwf_mpi(cplex,invpot,fofg_out,fofg,fofr,&
&    gbound_k,gbound_k,istwfk_one,kg_k,kg_k,me_g0,mgfft,ngfft,fftabs,n1,n2,n3,&
&    npw_k,npw_k,n4,n5,n6,ndat,2,weight_r,weight_i,comm_fft,cplexwf=cplexwf)

   ! Check if we got the initial u(g) within ATOL_DP
   ierr = COUNT(ABS(fofg - ref_fofg) > ATOL_DP)
   call xmpi_sum(ierr,comm_fft,mpierr)
   nfailed = nfailed + ierr

   write(info,"(a,i1,a)")sjoin(library,"<G|vloc|u>, istwfk "),istwf_k," :"

   if (ierr /= 0) then
     max_abserr = MAXVAL(ABS(fofg - ref_fofg))
     call xmpi_max(max_abserr,comm_fft,mpierr)
     write(msg,"(a,es9.2,a)")" FAILED (max_abserr = ",max_abserr,")"
     !if (me_fft == 0) write(std_out,*)(fofg(:,ig),ref_fofg(:,ig), ig=1,npw_k*ndat)
   else
     write(msg,"(a)")" OK"
   end if
   call wrtout(ount,sjoin(info,msg),"COLL")

   ABI_FREE(fofg_out)

   ABI_FREE(pot)
   ABI_FREE(invpot)

   ABI_FREE(kg_k)
   ABI_FREE(gbound_k)

   ABI_FREE(fofg)
   ABI_FREE(ref_fofg)
   ABI_FREE(fofr)

   ABI_FREE(full_kg_k)
   ABI_FREE(full_fofg)

   call destroy_distribfft(fftabs)
 end do

 if (nthreads > 0) then
   call xomp_set_num_threads(old_nthreads)
 end if

end function fftu_mpi_utests
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fourdp_mpi
!! NAME
!! fourdp_mpi
!!
!! FUNCTION
!! Conduct Fourier transform of REAL or COMPLEX function f(r)=fofr defined on
!! fft grid in real space, to create complex f(G)=fofg defined on full fft grid
!! in reciprocal space, in full storage mode, or the reverse operation.
!! For the reverse operation, the final data is divided by nfftot.
!! REAL case when cplex=1, COMPLEX case when cplex=2
!! Usually used for density and potentials.
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! ndat=Numbre of FFT transforms
!! isign=sign of Fourier transform exponent: current convention uses
!!    +1 for transforming from G to r
!!    -1 for transforming from r to G.
!! fftn2_distrib(2),ffti2_local(2)
!! fftn3_distrib(3),ffti3_local(3)
!! comm_fft=MPI communicator
!!
!! SIDE EFFECTS
!! Input/Output
!! fofg(2,nfft)=f(G), complex.
!! fofr(cplex*nfft)=input function f(r) (real or complex)
!!
!! TODO
!!  Write simplified API for sequential version.
!!
!! PARENTS
!!      fourdp
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fourdp_mpi(cplex,nfft,ngfft,ndat,isign,&
&  fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fourdp_mpi'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,isign,nfft,ndat,comm_fft
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(in) :: fftn2_distrib(ngfft(2)),ffti2_local(ngfft(2))
 integer,intent(in) :: fftn3_distrib(ngfft(3)),ffti3_local(ngfft(3))
 real(dp),intent(inout) :: fofg(2,nfft*ndat),fofr(cplex*nfft*ndat)

!Local variables-------------------------------
!scalars
 integer :: fftalg,fftalga,fftalgc
 character(len=500) :: msg

! *************************************************************************

 fftalg=ngfft(7); fftalga=fftalg/100 ; fftalgc=mod(fftalg,10)

 select case (fftalga)
 case (FFT_SG2002)
   call sg2002_mpifourdp(cplex,nfft,ngfft,ndat,isign,&
&    fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)

 case (FFT_FFTW3)
   call fftw3_mpifourdp(cplex,nfft,ngfft,ndat,isign,&
&    fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)

 !case (FFT_DFTI)
 !   call dfti_mpifourdp(cplex,nfft,ngfft,ndat,isign,&
 !&    fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local,fofg,fofr,comm_fft)

 case default
   write(msg,"(a,i0)")"Wrong fftalg: ",fftalg
   MSG_BUG(msg)
 end select

end subroutine fourdp_mpi
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fourwf_mpi
!!
!! NAME
!! fourwf_mpi
!!
!! FUNCTION
!! Carry out composite Fourier transforms between real and reciprocal (G) space.
!! Wavefunctions, contained in a sphere in reciprocal space,
!! can be FFT to real space. They can also be FFT from real space
!! to a sphere. Also, the density maybe accumulated, and a local
!! potential can be applied.
!!
!! The different options are :
!! - reciprocal to real space and output the result (option=0),
!! - reciprocal to real space and accumulate the density (option=1)
!! - reciprocal to real space, apply the local potential to the wavefunction
!!    in real space and produce the result in reciprocal space (option=2)
!! - real space to reciprocal space (option=3).
!!
!! Schedule of operations
!!(read first the description of the fftalg input variable in abinit_help)
!! - fftalgc=1 : use separate forward and backward transforms
!!     (7/12 savings in execution time);
!! - fftalgc=2 : in case of option=1 and option=2, use routines for composite operation
!!     even faster than 1x1
!!
!! Also for better speed, it uses no F90 construct, except the allocate command and for zeroing arrays.
!!
!! INPUTS
!! cplex= if 1 , denpot is real, if 2 , denpot is complex
!!    (cplex=2 only allowed for option=2, and istwf_k=1)
!!    not relevant if option=0 or option=3, so cplex=0 can be used to minimize memory
!! fftalgc=1 or 2 => simple or composite FFT applications
!! fofgin(2,npwin)=holds input wavefunction in G vector basis sphere.
!!                 (intent(in) but the routine sphere can modify it for another iflag)
!! gboundin(2*mgfft+8,2)=sphere boundary info for reciprocal to real space
!! gboundout(2*mgfft+8,2)=sphere boundary info for real to reciprocal space
!! istwf_k=option parameter that describes the storage of wfs
!! kg_kin(3,npwin)=reduced planewave coordinates, input
!! kg_kout(3,npwout)=reduced planewave coordinates, output
!! me_g0=1 if this MPI node treats the Gamma, 0 otherwise
!! mgfft=maximum size of 1D FFTs
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! distribfft<distribfft_type>=Tables needed for the FFT parallelism
!! n1,n2,n3=1D FFT sizes
!! npwin=number of elements in fofgin array (for option 0, 1 and 2)
!! npwout=number of elements in fofgout array (for option 2 and 3)
!! n4,n5,n6=dimensions of fofr.
!! ndat=Number of FFTs
!! option= if 0: do direct FFT
!!         if 1: do direct FFT, then sum the density
!!         if 2: do direct FFT, multiply by the potential, then do reverse FFT
!!         if 3: do reverse FFT only
!! weight_r=weight to be used for the accumulation of the density in real space
!!         (needed only when option=1)
!! weight_i=weight to be used for the accumulation of the density in real space
!!         (needed only when option=1)
!! comm_fft=MPI communicator.
!! [cplexwf]= 1 c2r or r2c can be used. Default: 2 i.e. use complex-to-complex FFTs.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! for option==0, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                fofr(2,n4,n5,n6) contains the output Fourier Transform of fofgin;
!!                no use of denpot, fofgout and npwout.
!! for option==1, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                denpot(cplex*n4,n5,n6) contains the input density at input,
!!                and the updated density at output (accumulated);
!!                no use of fofgout and npwout.
!! for option==2, fofgin(2,npwin)=holds input wavefunction in G sphere;
!!                denpot(cplex*n4,n5,n6) contains the input local potential;
!!                fofgout(2,npwout) contains the output function;
!! for option==3, fofr(2,n4,n5,n6) contains the input real space wavefunction;
!!                fofgout(2,npwout) contains its output Fourier transform;
!!                no use of fofgin and npwin.
!!
!! PARENTS
!!      fourwf,m_fft
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

subroutine fourwf_mpi(cplex,denpot,fofgin,fofgout,fofr,&
&  gboundin,gboundout,istwf_k,kg_kin,kg_kout,me_g0,mgfft,ngfft,distribfft,n1,n2,n3,&
&  npwin,npwout,n4,n5,n6,ndat,option,weight_r,weight_i,comm_fft,cplexwf)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fourwf_mpi'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istwf_k,mgfft,n1,n2,n3,n4,n5,n6,npwin,ndat
 integer,intent(in) :: npwout,option,comm_fft,me_g0
 integer,intent(in),optional :: cplexwf
 real(dp),intent(in) :: weight_r,weight_i
 type(distribfft_type),intent(in) :: distribfft
!arrays
 integer,intent(in) :: gboundin(2*mgfft+8,2),gboundout(2*mgfft+8,2),ngfft(18)
 integer,intent(in) :: kg_kin(3,npwin),kg_kout(3,npwout)
 real(dp),intent(inout) :: denpot(cplex*n4,n5,n6),fofgin(2,npwin*ndat)
 real(dp),intent(inout) :: fofr(2,n4,n5,n6*ndat)
 real(dp),intent(out) :: fofgout(2,npwout*ndat)

!Local variables-------------------------------
!scalars
 integer :: fftalg,fftalga,fftalgc,idat
 integer :: cplexwf_,i1,i2,i3,iflag,ig,igdat,me_fft,m1i,m1o,m2i,m2o,m3i
 integer :: m3o,max1i,max1o,max2i,max2i_plus,max2o,max2o_plus
 integer :: max3i,max3o,md1,md1i,md1o,md2,md2i,md2o,md2proc
 integer :: md3,md3i,md3o,min1i,min1o,min2i,min2i_moins,min2o,min2o_moins,min3i,min3o
 integer :: nd3proc,nproc_fft,n6eff,i3_glob,i2_loc,i2dat_loc,i3dat
 integer,save :: nwrites_ialltoall=0
 logical :: use_ialltoall
 real(dp) :: fim,fre,xnorm
 character(len=500) :: msg
!arrays
 integer,parameter :: shiftg0(3)=0
 integer,parameter :: symmE(3,3)=reshape([1,0,0,0,1,0,0,0,1],[3,3])
! real(dp) :: tsec(2)
 real(dp) :: weight_array_r(ndat), weight_array_i(ndat)
 real(dp),allocatable :: workf(:,:,:,:)

! *************************************************************************

 !call timab(540,1,tsec)

 fftalg=ngfft(7); fftalga=fftalg/100 ; fftalgc=mod(fftalg,10)

 me_fft=ngfft(11); nproc_fft=ngfft(10)
 !write(std_out,*)"nproc_fft",nproc_fft
 ! Risky since I don't know if this entry is initialized
 ! Continue to pass comm_fft explicitly.
 !comm_fft = ngfft(14)

 if (fftalgc<1 .or. fftalgc>2) then
   write(msg,'(a,i0,3a)')&
&   'The input algorithm number fftalgc=',fftalgc,' is not allowed with MPI-FFT. Must be 1 or 2',ch10,&
&   'Action: change fftalgc in your input file.'
   MSG_ERROR(msg)
 end if

 if (option<0 .or. option>3) then
   write(msg,'(a,i0,3a)')&
&   'The option number',option,' is not allowed.',ch10,&
&   'Only option=0, 1, 2 or 3 are allowed presently.'
   MSG_ERROR(msg)
 end if

 if (option==1 .and. cplex/=1) then
   write(msg,'(a,i0,a)')&
&   'With the option number 1, cplex must be 1 but it is cplex=',cplex,'.'
   MSG_ERROR(msg)
 end if

 if ( ALL(cplex/=(/1,2/)) .and. ANY(option==(/1,2/)) ) then
   write(msg,'(a,i0,a)')' When option is (1,2) cplex must be 1 or 2, but it is cplex=',cplex,'.'
   MSG_ERROR(msg)
 end if

 !write(std_out,*)"in fourwf_mpi with fftalg: ",fftalg,fftalgc

! We use the non-blocking version if IALLTOALL is available and ndat > 1.
 use_ialltoall = .False.
#ifdef HAVE_MPI_IALLTOALL
 use_ialltoall = (ndat > 1)
#endif
 use_ialltoall = (use_ialltoall .and. ALLOW_IALLTOALL)
 if (use_ialltoall .and. nwrites_ialltoall==0) then
   nwrites_ialltoall = 1
   call wrtout(std_out,"- Will use non-blocking ialltoall for MPI-FFT","COLL")
 end if

 md1i=0; md2i=0; md3i=0; m2i=0
 md1o=0; md2o=0; md3o=0; m2o=0

 if (option/=3) then
   ! Compute the dimensions of the small-box enclosing the input G-sphere
   max1i=gboundin(2,1); min1i=gboundin(1,1)
   max2i=gboundin(4,1); min2i=gboundin(3,1)

   if(istwf_k==2 .or. istwf_k==4 .or. istwf_k==6 .or. istwf_k==8)then
     max1i=max(max1i,-min1i)
     min1i=-max1i
   else if (istwf_k==3 .or. istwf_k==5 .or. istwf_k==7 .or. istwf_k==9) then
     max1i=max(max1i,-min1i-1)
     min1i=-max1i-1
   end if
   if (istwf_k>=2 .and. istwf_k<=5) then
     max2i=max(max2i,-min2i)
     min2i=-max2i
   else if (istwf_k>=6 .and. istwf_k<=9) then
     max2i=max(max2i,-min2i-1)
     min2i=-max2i-1
   end if

   max3i=gboundin(4,2); min3i=gboundin(3,2)

   ! Compute arrays size and leading dimensions to avoid cache trashing
   m1i=max1i-min1i+1; md1i=2*(m1i/2)+1
   m2i=max2i-min2i+1; md2i=2*(m2i/2)+1

   !if (.False.) then
   if (nproc_fft/=1) then
     ! Increase max2i in order to have m2i divisible by nproc_fft
     min2i_moins=(((m2i-1)/nproc_fft+1)*nproc_fft-m2i)/2
     max2i_plus=((m2i-1)/nproc_fft+1)*nproc_fft-m2i-min2i_moins
     ! max2i=max2i+((m2i-1)/nproc_fft+1)*nproc_fft-m2i
     max2i=max2i+max2i_plus
     min2i=min2i-min2i_moins
     ! careful, to be checked and make sure the max and min are smaller than size of box
     m2i=max2i-min2i+1; md2i=2*(m2i/2)+1
   end if
   ABI_CHECK(m2i <= n2, "m2i > n2")

   m3i=max3i-min3i+1; md3i=2*(m3i/2)+1
 end if

 if (option==2 .or. option==3) then
   ! Compute the dimensions of the small-box enclosing the output G-sphere
   max1o=gboundout(2,1); min1o=gboundout(1,1)
   max2o=gboundout(4,1); min2o=gboundout(3,1)

   if (istwf_k==2 .or. istwf_k==4 .or. istwf_k==6 .or. istwf_k==8) then
     max1o=max(max1o,-min1o)
     min1o=-max1o
   else if (istwf_k==3 .or. istwf_k==5 .or. istwf_k==7 .or. istwf_k==9) then
     max1o=max(max1o,-min1o-1)
     min1o=-max1o-1
   end if
   if (istwf_k>=2 .and. istwf_k<=5) then
     max2o=max(max2o,-min2o)
     min2o=-max2o
   else if (istwf_k>=6 .and. istwf_k<=9) then
     max2o=max(max2o,-min2o-1)
     min2o=-max2o-1
   end if

   max3o=gboundout(4,2); min3o=gboundout(3,2)

   ! Compute arrays size and leading dimensions to avoid cache trashing
   m1o=max1o-min1o+1; md1o=2*(m1o/2)+1
   m2o=max2o-min2o+1; md2o=2*(m2o/2)+1

   if (nproc_fft/=1) then
     ! Increase max2o in order to have m2o divisible by nproc_fft
     min2o_moins=(((m2o-1)/nproc_fft+1)*nproc_fft-m2o)/2
     max2o_plus=((m2o-1)/nproc_fft+1)*nproc_fft-m2o-min2o_moins
     ! max2o=max2o+((m2o-1)/nproc_fft+1)*nproc_fft-m2o
     max2o=max2o+max2o_plus
     min2o=min2o-min2o_moins
     ! careful, to be checked and make sure the max and min are smaller than size of box
     m2o=max2o-min2o+1; md2o=2*(m2o/2)+1
   end if
   ABI_CHECK(m2o <= n2, "m2o > n2")

   m3o=max3o-min3o+1; md3o=2*(m3o/2)+1
 end if

 md1=max(md1i,md1o)
 md2=max(md2i,md2o)
 md3=max(md3i,md3o)

 md2proc=(max(m2i,m2o)-1)/nproc_fft+1
 n6eff=(n6-1)/nproc_fft+1
 !write(std_out,*)'fourwf_mpi : max1i,max2i,max3i=',max1i,max2i,max3i
 !write(std_out,*)'fourwf_mpi : min1i,min2i,min3i=',min1i,min2i,min3i
 !write(std_out,*)'fourwf_mpi : m1i,m2i,m3i=',m1i,m2i,m3i

 ! Allocate work array in G-space (note exchange 3 <--> 2)
 ABI_ALLOCATE(workf,(2,md1,md3,md2proc*ndat))

 if (option/=3) then
   ! Insert fofgin into the **small** box (array workf) :
   ! Note the switch of md3 and md2, as they are only needed to dimension workf inside "sphere"

   if (nproc_fft > 1) then
     if (istwf_k/=1 )then
       write(msg,'(a,i0,a)')'The value of istwf_k: ',istwf_k,' is not allowed. Only istwf_k=1 is allowed in MPI-FFT'
       !MSG_WARNING(msg)
       MSG_ERROR(msg)
     end if
     call sphere_fft1(fofgin,ndat,npwin,workf,m1i,m2i,m3i,md1,md3,md2proc,kg_kin,distribfft%tab_fftwf2_local)
   else
     iflag=2
     call sphere(fofgin,ndat,npwin,workf,m1i,m2i,m3i,md1,md3,md2proc,kg_kin,istwf_k,iflag,me_g0,shiftg0,symmE,one)
   end if
 end if

 ! Can we use c2r or r2c?
 cplexwf_=2; if (istwf_k==2) cplexwf_=1
 if (present(cplexwf)) cplexwf_ = cplexwf

 if (option==0 .or. ((option==1.or.option==2) .and. fftalgc==1) .or. option==3) then
   ! Treat non-composite operations

   if (option/=3) then
     ! Fourier transform workf(2,md1,md3,md2proc*ndat) to fofr (reciprocal to real space).
     ! FIXME: This is buggy if cplexwx==1

     select case (fftalga)

     case (FFT_SG2002)

!        do idat=1,ndat
!        call sg2002_mpiback_wf(cplexwf_,1,n1,n2,n3,n4,n5,(n6-1)/nproc_fft+1,&
! &        max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3, &
! &         workf(:,:,:,(idat-1)*md2proc+1:idat*md2proc), &
! &        fofr(:,:,:,(idat-1)*n6eff+1:idat*n6eff),comm_fft)
!        enddo
       call sg2002_mpiback_wf(cplexwf_,ndat,n1,n2,n3,n4,n5,(n6-1)/nproc_fft+1,&
&        max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,workf,fofr,comm_fft)

     case (FFT_FFTW3)

       if (use_ialltoall) then
         call fftw3_mpiback_manywf(cplexwf_,ndat,n1,n2,n3,n4,n5,(n6-1)/nproc_fft+1,&
&          max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,workf,fofr,comm_fft)
       else
         call fftw3_mpiback_wf(cplexwf_,ndat,n1,n2,n3,n4,n5,(n6-1)/nproc_fft+1,&
&          max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,workf,fofr,comm_fft)
       end if

     !case (FFT_DFTI)
     !  call dfti_mpiback_wf(cplexwf_,ndat,n1,n2,n3,n4,n5,(n6-1)/nproc_fft+1,&
     !   &        max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,workf,fofr,comm_fft)
     case default
       MSG_ERROR("fftalga does not provide MPI back_wf")
     end select

   end if ! option

   nd3proc=(n3-1)/nproc_fft+1

   if (option==1) then
     ! Accumulate density
     do idat=1,ndat
       do i3=1,nd3proc
         i3_glob = i3 + nd3proc * me_fft
         i3dat = i3  + n6eff * (idat-1)
         do i2=1,n2
           do i1=1,n1
             denpot(i1,i2,i3_glob) = denpot(i1,i2,i3_glob) &
&              + (weight_r*fofr(1,i1,i2,i3dat)**2+ weight_i*fofr(2,i1,i2,i3dat)**2)
           end do
         end do
       end do
     end do
   end if ! option==1

   if (option==2) then

     !  Apply local potential
     if (cplex==1) then
       do idat=1,ndat
         do i3=1,nd3proc
           i3_glob = i3 + nd3proc * me_fft
           i3dat = i3  + n6eff * (idat-1)
           do i2=1,n2
             do i1=1,n1
               fofr(1,i1,i2,i3dat)=denpot(i1,i2,i3_glob)*fofr(1,i1,i2,i3dat)
               fofr(2,i1,i2,i3dat)=denpot(i1,i2,i3_glob)*fofr(2,i1,i2,i3dat)
             end do
           end do
         end do
       end do

     else if (cplex==2) then
       do idat=1,ndat
         do i3=1,(n3-1)/nproc_fft+1
           i3_glob = i3 + nd3proc * me_fft
           i3dat = i3  + n6eff * (idat-1)
           do i2=1,n2
             do i1=1,n1
               fre=fofr(1,i1,i2,i3dat)
               fim=fofr(2,i1,i2,i3dat)
               fofr(1,i1,i2,i3dat)=denpot(2*i1-1,i2,i3_glob)*fre -denpot(2*i1,i2,i3_glob)*fim
               fofr(2,i1,i2,i3dat)=denpot(2*i1-1,i2,i3_glob)*fim +denpot(2*i1,i2,i3_glob)*fre
             end do
           end do
         end do
       end do
     end if ! cplex

   end if ! option==2

   if (option==2 .or. option==3) then
     ! Fourier transform fofr to workf (real to reciprocal space)
     ! output in workf(2,md1,md3,md2proc*ndat)

     select case (fftalga)
     case (FFT_SG2002)

       call sg2002_mpiforw_wf(cplexwf_,ndat,n1,n2,n3,n4,n5,(n6-1)/nproc_fft+1,&
&        max1o,max2o,max3o,m1o,m2o,m3o,md1,md2proc,md3,fofr,workf,comm_fft)

     case (FFT_FFTW3)

       if (use_ialltoall) then
         call fftw3_mpiforw_manywf(cplexwf_,ndat,n1,n2,n3,n4,n5,(n6-1)/nproc_fft+1,&
&          max1o,max2o,max3o,m1o,m2o,m3o,md1,md2proc,md3,fofr,workf,comm_fft)
       else
         call fftw3_mpiforw_wf(cplexwf_,ndat,n1,n2,n3,n4,n5,(n6-1)/nproc_fft+1,&
&          max1o,max2o,max3o,m1o,m2o,m3o,md1,md2proc,md3,fofr,workf,comm_fft)
       end if

     !case (FFT_DFTI)
     !  call dfti_mpiforw_wf(cplexwf_,ndat,n1,n2,n3,n4,n5,(n6-1)/nproc_fft+1,&
     !  &max1o,max2o,max3o,m1o,m2o,m3o,md1,md2proc,md3,fofr,workf,comm_fft)

     case default
       MSG_ERROR("fftalga does not provide MPI back_wf")
     end select

   end if

 else if (fftalgc==2 .and. (option==1 .or. option==2)) then
   !  Treat composite operations

   select case (option)
   case (1)
     !ABI_CHECK(weight_r == weight_i,"weight_r != weight_i")
     weight_array_r(:)=weight_r
     weight_array_i(:)=weight_i

     select case (fftalga)
     case (FFT_SG2002)
       ! Note that here we don' fill fofr. Don't know if someone in
       ! abinit uses option 1 to get both fofr as well as denpot
       call sg2002_accrho(cplexwf_,ndat,n1,n2,n3,n4,n5,n6,(n6-1)/nproc_fft+1,&
&        max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,comm_fft,nproc_fft,me_fft,&
&        workf,denpot,weight_array_r,weight_array_i)

     case (FFT_FFTW3)
       call fftw3_accrho(cplexwf_,ndat,n1,n2,n3,n4,n5,n6,(n6-1)/nproc_fft+1,&
&        max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,comm_fft,nproc_fft,me_fft,&
&        workf,denpot,weight_array_r, weight_array_i)

     case default
       MSG_ERROR("fftalga does not provide accrho")
     end select

   case (2)
     !write(std_out,*)fftalg,option,cplex
     !ABI_CHECK(cplex==1,"cplex!=2 with fftalg 412 is buggy")

     select case (fftalga)

     case (FFT_SG2002)

       if (use_ialltoall) then
         call sg2002_applypot_many(cplexwf_,cplex,ndat,n1,n2,n3,n4,n5,n6,(n6-1)/nproc_fft+1,&
&          max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,&
&          max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft,denpot,workf)

       else
         call sg2002_applypot(cplexwf_,cplex,ndat,n1,n2,n3,n4,n5,n6,(n6-1)/nproc_fft+1,&
&          max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,&
&          max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft,denpot,workf)
       endif

     case (FFT_FFTW3)

       if (use_ialltoall) then

         call fftw3_applypot_many(cplexwf_,cplex,ndat,n1,n2,n3,n4,n5,n6,(n6-1)/nproc_fft+1,&
&          max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,&
&          max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft,denpot,workf)

       else
         call fftw3_applypot(cplexwf_,cplex,ndat,n1,n2,n3,n4,n5,n6,(n6-1)/nproc_fft+1,&
&          max1i,max2i,max3i,m1i,m2i,m3i,md1,md2proc,md3,&
&          max1o,max2o,max3o,m1o,m2o,m3o,comm_fft,nproc_fft,me_fft,denpot,workf)
       end if

     case default
       MSG_ERROR("fftalga does not provide applypot")
     end select

   case default
     write(msg,"(a,i0,a)")"Option ",option," is not supported when fftalgc == 2"
     MSG_ERROR(msg)
   end select

 end if ! End of composite operations

 if (option==2 .or. option==3) then
   ! From FFT box to kg_kout.
   xnorm=one/dble(n1*n2*n3)

   if (nproc_fft > 1) then
     do idat=1,ndat
       do ig=1,npwout
         i1=kg_kout(1,ig); if(i1<0) i1=i1+m1o; i1=i1+1
         i2=kg_kout(2,ig); if(i2<0) i2=i2+m2o; i2=i2+1
         i3=kg_kout(3,ig); if(i3<0) i3=i3+m3o; i3=i3+1

         igdat = ig + (idat-1) * npwout
         i2_loc = (i2-1)/nproc_fft +1
         i2dat_loc = i2_loc + (idat-1) * md2proc

         fofgout(1,igdat)=workf(1,i1,i3,i2dat_loc)*xnorm
         fofgout(2,igdat)=workf(2,i1,i3,i2dat_loc)*xnorm
       end do
     end do
   else
     ! Warning: This call is buggy if istwfk > 2
     iflag=-2
     call sphere(fofgout,ndat,npwout,workf,m1o,m2o,m3o,md1,md3,md2proc,kg_kout,istwf_k,iflag,&
&      me_g0,shiftg0,symmE,xnorm)
   end if
 end if ! if option==2 or 3

 ABI_DEALLOCATE(workf)

!call timab(540,2,tsec)

end subroutine fourwf_mpi
!!***

!----------------------------------------------------------------------

!!****f* m_fft/fftmpi_u
!!
!! NAME
!! fftmpi_u
!!
!! FUNCTION
!!  **** DO NOT USE IT ****
!!  Wrapper function used to perform the MPI FFT of ndat wavefunctions.
!!  Mainly used for unit tests and prototyping. Better integration
!!  will be provided afterwards.
!!
!! INPUTS
!!  See fourwf_mpi
!!
!! OUTPUT
!!  See fourwf_mpi
!!
!! PARENTS
!!      m_fft
!!
!! CHILDREN
!!      fourwf_mpi
!!
!! SOURCE

! The final interface should be:
!subroutine fftmpi_u(npw_k,n4,n5,n6,nspinor,ndat,mgfft,ngfft,istwf_k,kg_k,gbound_k,fftabs,isign,fofg,fofr)

subroutine fftmpi_u(npw_k,n4,n5,n6,ndat,mgfft,ngfft,&
&  istwf_k,gbound_k,kg_k,me_g0,distribfft,isign,fofg,fofr,comm_fft,cplexwf)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fftmpi_u'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,mgfft,n4,n5,n6,ndat,npw_k
 integer,intent(in) :: isign,comm_fft,me_g0,cplexwf
 type(distribfft_type),intent(in) :: distribfft
!arrays
 integer,intent(in) :: gbound_k(2*mgfft+8,2),ngfft(18)
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(inout) :: fofg(2,npw_k*ndat),fofr(2,n4,n5,n6*ndat)

!Local variables-------------------------------
!scalars
 integer,parameter :: npw0=0,cplex0=0
 integer :: n1,n2,n3
 real(dp),parameter :: weight_r=one,weight_i=one
!arrays
 integer :: dummy_kg(0,0)
 real(dp) :: dummy_denpot(0,0,0),dummy_fofg(0,0)

! *************************************************************************

 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)

 if (isign == 1) then
   ! option 0 G --> R
   call fourwf_mpi(cplex0,dummy_denpot,fofg,dummy_fofg,fofr,&
&    gbound_k,gbound_k,istwf_k,kg_k,dummy_kg,me_g0,mgfft,ngfft,distribfft,n1,n2,n3,&
&    npw_k,npw0,n4,n5,n6,ndat,0,weight_r,weight_i,comm_fft,cplexwf=cplexwf)
 else
   ! option 3 R --> G
   call fourwf_mpi(cplex0,dummy_denpot,dummy_fofg,fofg,fofr,&
&    gbound_k,gbound_k,istwf_k,dummy_kg,kg_k,me_g0,mgfft,ngfft,distribfft,n1,n2,n3,&
&    npw0,npw_k,n4,n5,n6,ndat,3,weight_r,weight_i,comm_fft,cplexwf=cplexwf)
 end if

end subroutine fftmpi_u
!!***

END MODULE m_fft
!!***
