!!****m* ABINIT/m_distribfft
!! NAME
!!  m_distribfft
!!
!! FUNCTION
!!  This module provides the definition of the different arrays
!!  used for FFT parallelization with MPI and n2 plane sharing
!!
!! COPYRIGHT
!! Copyright (C) 2011-2019 ABINIT group (FD,MT)
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

MODULE m_distribfft

 use defs_basis
 use m_abicore
 use m_errors

 implicit none

 private
!!***

!!****t* m_distribfft/distribfft_type
!! NAME
!! distribfft_type
!!
!! FUNCTION
!! The distribfft_type structured datatype gather different information
!! for plane sharing for FFT parallelization
!!
!! TODO
!!   1) One should create two separated tables: one for the wavefunctions and the other one for
!!      fourdp on the dense/coarse mesh.
!!
!!   2) Use shorter names --> fftabs_type
!!
!! SOURCE

 type, public :: distribfft_type

  integer :: n2_coarse=0
  ! Number of points along the y directions for the coarse FFT mesh

  integer :: n2_fine=0
  ! Number of points along the y directions for the dense FFT mesh

  !integer :: me_g0
  ! 1 if this MPI node has G=0, 0 otherwise.
  ! Needed for the FFTs of the wavefunctions.

  integer, allocatable :: tab_fftwf2_distrib(:)
  ! rank of the processors which own fft planes in 2nd dimension for fourwf

  integer, allocatable :: tab_fftdp2_distrib(:)
  ! rank of the processors which own fft planes in 2nd dimension for fourdp

  integer, allocatable :: tab_fftdp3_distrib(:)
  ! rank of the processors which own fft planes in 3rd dimension for fourdp

  integer, allocatable :: tab_fftwf2dg_distrib(:)
  ! rank of the processors which own fft planes in 2nd dimension for fourwf on fine grid

  integer, allocatable :: tab_fftdp2dg_distrib(:)
  ! rank of the processors which own fft planes in 2nd dimension for fourdp on fine grid

  integer, allocatable :: tab_fftdp3dg_distrib(:)
  ! rank of the processors which own fft planes in 3rd dimension for fourdp on fine grid

  integer, allocatable :: tab_fftwf2_local(:)
  ! local i2 indices in fourwf

  integer, allocatable :: tab_fftdp2_local(:)
  ! local i2 indices in fourdp

  integer, allocatable :: tab_fftdp3_local(:)
  ! local i3 indices in fourdp

  integer, allocatable :: tab_fftwf2dg_local(:)
  ! local i2 indices in fourwf on fine grid

  integer, allocatable :: tab_fftdp2dg_local(:)
  ! local i2 indices in fourdp on fine grid

  integer, allocatable :: tab_fftdp3dg_local(:)
  ! local i3 indices in fourdp on fine grid

end type distribfft_type

 public :: init_distribfft         ! Initializes mpi information for FFT distribution.
 public :: init_distribfft_seq     ! Initializes a sequential FFT distribution.
 public :: destroy_distribfft      ! Free dynamic memory.
 public :: copy_distribfft         ! Copy datatype.
!!***

CONTAINS !===========================================================
!!***

!!****f* m_distribfft/init_distribfft
!! NAME
!!  init_distribfft
!!
!! FUNCTION
!!  Initializes mpi information for FFT distribution
!!  Note that we always use cyclic distribution mode for the wavefunctions in G-space.
!!  MPI-FFT routines should always be compatible with this distribution.
!!
!! INPUTS
!! grid_type = 'c' or 'f' for information about coarse or fine fft grid
!! nproc_fft = number of process used to distribute the fft
!! n2,n3     = sizes of second and third fft grid
!!
!! SIDE EFFECTS
!!  distribfft = instance of distribfft_type to initialize
!!  Update of "fft distrib" tabs accordingly to the fft parallelisation
!!
!! PARENTS
!!      m_fft,m_fft_prof,m_qparticles,m_wfd,mpi_setup,scfcv,vtorhorec
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_distribfft(distribfft_arg,grid_type,nproc_fft,n2,n3)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nproc_fft,n2,n3
 character(len=1),intent(in) :: grid_type
 type(distribfft_type), intent(inout) :: distribfft_arg

!Local variables-------------------------------
!scalars
 integer :: i2,i3,n2_local,n3_local
 !character(len=500) :: msg

! ***********************************************************************

 DBG_ENTER("COLL")

 !local sizes
 n2_local = n2 / nproc_fft
 n3_local = n3 / nproc_fft

 select case (grid_type)
 case ('c')
    ! Updating information about coarse fft grid
    if(distribfft_arg%n2_coarse > 0) then
      if(n2 == distribfft_arg%n2_coarse) then
        MSG_WARNING("The distribfft passed was already allocated for coarse grid on the same size")
        return
      else
        MSG_ERROR("The distribfft passed was already allocated for coarse grid")
      endif
    end if
    distribfft_arg%n2_coarse = n2
    ! Initialisation of fft distrib tab
    ABI_ALLOCATE(distribfft_arg%tab_fftwf2_distrib,(n2))
    ABI_ALLOCATE(distribfft_arg%tab_fftwf2_local,(n2))
    ABI_ALLOCATE(distribfft_arg%tab_fftdp2_distrib,(n2))
    ABI_ALLOCATE(distribfft_arg%tab_fftdp2_local,(n2))
    ABI_ALLOCATE(distribfft_arg%tab_fftdp3_distrib,(n3))
    ABI_ALLOCATE(distribfft_arg%tab_fftdp3_local,(n3))
    do i2=1, n2
      ! Cyclic distribution of ig2 planes over fft processors
      distribfft_arg%tab_fftwf2_distrib(i2) = modulo((i2-1),nproc_fft)
      distribfft_arg%tab_fftwf2_local(i2)    = (i2-1)/nproc_fft + 1
      ! Block distribution of i2 planes over fft processors for fourdp
      distribfft_arg%tab_fftdp2_distrib(i2) = (i2-1) /  n2_local
      distribfft_arg%tab_fftdp2_local(i2)   = modulo((i2-1),n2_local) + 1
    end do
    do i3=1, n3
      ! Block distribution of i3 planes over fft processors for fourdp
      distribfft_arg%tab_fftdp3_distrib(i3) = (i3-1) /  n3_local
      distribfft_arg%tab_fftdp3_local(i3)   = modulo((i3-1),n3_local) + 1
    end do

 case ('f')
    if(distribfft_arg%n2_fine > 0) then
      if(n2 == distribfft_arg%n2_fine) then
        MSG_WARNING("The distribfft passed was already allocated for fine grid on the same size")
        return
      else
        MSG_ERROR("The distribfft passed was already allocated for fine grid")
      end if
    endif
    distribfft_arg%n2_fine = n2
    ! Updating information about fine fft grid
    ABI_ALLOCATE(distribfft_arg%tab_fftwf2dg_distrib,(n2))
    ABI_ALLOCATE(distribfft_arg%tab_fftwf2dg_local,(n2))
    ABI_ALLOCATE(distribfft_arg%tab_fftdp2dg_distrib,(n2))
    ABI_ALLOCATE(distribfft_arg%tab_fftdp2dg_local,(n2))
    ABI_ALLOCATE(distribfft_arg%tab_fftdp3dg_distrib,(n3))
    ABI_ALLOCATE(distribfft_arg%tab_fftdp3dg_local,(n3))
    do i2=1, n2
      ! Cyclic distribution of ig2 planes over fft processors
      distribfft_arg%tab_fftwf2dg_distrib(i2) = modulo((i2-1),nproc_fft)
      distribfft_arg%tab_fftwf2dg_local(i2)    = (i2-1)/nproc_fft + 1
      ! Block distribution of i2 planes over fft processors for fourdp on fine grid
      distribfft_arg%tab_fftdp2dg_distrib(i2) = (i2-1) /  n2_local
      distribfft_arg%tab_fftdp2dg_local(i2)   = modulo((i2-1),n2_local) + 1
    end do
    do i3=1, n3
      ! Block distribution of i3 planes over fft processors for fourdp on fine grid
      distribfft_arg%tab_fftdp3dg_distrib(i3) = (i3-1) /  n3_local
      distribfft_arg%tab_fftdp3dg_local(i3)   = modulo((i3-1),n3_local) + 1
    end do

 case default
    MSG_ERROR("Unknown kind of fft grid! Only 'c' for coarse grid and 'f' for fine grid are allowed")
 end select

 ! One needs to know if this node has G=0 when we do the FFTs of the wavefunctions
 !distribfft_arg%me_g0 = 0
 !if (distribfft_arg%tab_fftwf2_distrib(1) == me_fft) distribfft_arg%me_g0 = 1

 DBG_EXIT("COLL")

end subroutine init_distribfft
!!***

!===========================================================

!!****f* m_distribfft/init_distribfft_seq
!! NAME
!!  init_distribfft_seq
!!
!! FUNCTION
!!  Initializes a sequential FFT distribution
!!
!! INPUTS
!! grid_type  = 'c' or 'f' for information about coarse or fine fft grid
!! n2,n3 = sizes of second and third fft grid
!! type_four = 'fourdp' or 'fourwf' or 'all' to prepare a call to fourdp/fourwf
!!
!! SIDE EFFECTS
!!  distribfft = instance of t_distribfft to initialize
!!  Update of "fft distrib" tabs accordingly to the fft parallelisation
!!
!! PARENTS
!!      atm2fft,bethe_salpeter,calc_vhxc_me,cut3d,dfpt_atm2fft,dieltcel,eph
!!      ks_ddiago,m_cut3d,m_dvdb,m_fft_prof,m_gsphere,m_ioarr,m_kxc,m_ppmodel
!!      m_screening,m_wfk,multipoles_fftr,pawgrnl,pawmknhat,pawmknhat_psipsi
!!      pawsushat,scfcv,screening,sigma,suscep_stat,susk,suskmm,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_distribfft_seq(distribfft_arg,grid_type,n2,n3,type_four)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: n2,n3
 character(len=1),intent(in) :: grid_type
 character(len=*),intent(in) :: type_four
 type(distribfft_type), intent(inout) :: distribfft_arg

!Local variables-------------------------------
!scalars
 integer :: ii

! ***********************************************************************

 DBG_ENTER("COLL")

 !distribfft_arg%me_g0 = 1

 select case (grid_type)
 case ('c')
   distribfft_arg%n2_coarse = n2
   if (type_four=='fourwf'.or.type_four(1:3)=='all') then
     if (allocated(distribfft_arg%tab_fftwf2_distrib)) then
       ABI_DEALLOCATE(distribfft_arg%tab_fftwf2_distrib)
     end if
     if (allocated(distribfft_arg%tab_fftwf2_local)) then
       ABI_DEALLOCATE(distribfft_arg%tab_fftwf2_local)
     end if
     ABI_ALLOCATE(distribfft_arg%tab_fftwf2_distrib,(n2))
     ABI_ALLOCATE(distribfft_arg%tab_fftwf2_local,(n2))
     distribfft_arg%tab_fftwf2_distrib=0
     distribfft_arg%tab_fftwf2_local=(/(ii,ii=1,n2)/)
   end if
   if (type_four=='fourdp'.or.type_four(1:3)=='all') then
     if (allocated(distribfft_arg%tab_fftdp2_distrib)) then
       ABI_DEALLOCATE(distribfft_arg%tab_fftdp2_distrib)
     end if
     if (allocated(distribfft_arg%tab_fftdp2_local)) then
       ABI_DEALLOCATE(distribfft_arg%tab_fftdp2_local)
     end if
     if (allocated(distribfft_arg%tab_fftdp3_distrib)) then
       ABI_DEALLOCATE(distribfft_arg%tab_fftdp3_distrib)
     end if
     if (allocated(distribfft_arg%tab_fftdp3_local)) then
       ABI_DEALLOCATE(distribfft_arg%tab_fftdp3_local)
     end if
     ABI_ALLOCATE(distribfft_arg%tab_fftdp2_distrib,(n2))
     ABI_ALLOCATE(distribfft_arg%tab_fftdp2_local,(n2))
     ABI_ALLOCATE(distribfft_arg%tab_fftdp3_distrib,(n3))
     ABI_ALLOCATE(distribfft_arg%tab_fftdp3_local,(n3))
     distribfft_arg%tab_fftdp2_distrib=0
     distribfft_arg%tab_fftdp3_distrib=0
     distribfft_arg%tab_fftdp2_local=(/(ii,ii=1,n2)/)
     distribfft_arg%tab_fftdp3_local=(/(ii,ii=1,n3)/)
   end if

 case ('f')
   distribfft_arg%n2_fine = n2
   if (type_four=='fourwf'.or.type_four(1:3)=='all') then
     if (allocated(distribfft_arg%tab_fftwf2dg_distrib)) then
       ABI_DEALLOCATE(distribfft_arg%tab_fftwf2dg_distrib)
     end if
     if (allocated(distribfft_arg%tab_fftwf2dg_local)) then
       ABI_DEALLOCATE(distribfft_arg%tab_fftwf2dg_local)
     end if
     ABI_ALLOCATE(distribfft_arg%tab_fftwf2dg_distrib,(n2))
     ABI_ALLOCATE(distribfft_arg%tab_fftwf2dg_local,(n2))
     distribfft_arg%tab_fftwf2dg_distrib=0
     distribfft_arg%tab_fftwf2dg_local=(/(ii,ii=1,n2)/)
   end if
   if (type_four=='fourdp'.or.type_four(1:3)=='all') then
     if (allocated(distribfft_arg%tab_fftdp2dg_distrib)) then
       ABI_DEALLOCATE(distribfft_arg%tab_fftdp2dg_distrib)
     end if
     if (allocated(distribfft_arg%tab_fftdp2dg_local)) then
       ABI_DEALLOCATE(distribfft_arg%tab_fftdp2dg_local)
     end if
     if (allocated(distribfft_arg%tab_fftdp3dg_distrib)) then
       ABI_DEALLOCATE(distribfft_arg%tab_fftdp3dg_distrib)
     end if
     if (allocated(distribfft_arg%tab_fftdp3dg_local)) then
       ABI_DEALLOCATE(distribfft_arg%tab_fftdp3dg_local)
     end if
     ABI_ALLOCATE(distribfft_arg%tab_fftdp2dg_distrib,(n2))
     ABI_ALLOCATE(distribfft_arg%tab_fftdp2dg_local,(n2))
     ABI_ALLOCATE(distribfft_arg%tab_fftdp3dg_distrib,(n3))
     ABI_ALLOCATE(distribfft_arg%tab_fftdp3dg_local,(n3))
     distribfft_arg%tab_fftdp2dg_distrib=0
     distribfft_arg%tab_fftdp3dg_distrib=0
     distribfft_arg%tab_fftdp2dg_local=(/(ii,ii=1,n2)/)
     distribfft_arg%tab_fftdp3dg_local=(/(ii,ii=1,n3)/)
   end if

 case default
    MSG_ERROR("Unknown kind of fft grid! Only 'c' for coarse grid and 'f' for fine grid are allowed")
 end select

 DBG_EXIT("COLL")

end subroutine init_distribfft_seq
!!***

!===========================================================

!!****f* m_distribfft/destroy_distribfft
!! NAME
!!  destroy_distribfft
!!
!! FUNCTION
!!  Cleans-up the mpi information for FFT distribution
!!  (mostly deallocate parts distribfft(:) ).
!!
!! SIDE EFFECTS
!!  mpi_enreg=information about MPI parallelization
!!
!! PARENTS
!!      atm2fft,dfpt_atm2fft,m_cut3d,m_fft,m_mpinfo,m_qparticles,m_wfd
!!      multipoles_fftr,pawgrnl,pawmknhat,pawmknhat_psipsi,pawsushat,scfcv
!!      vtorhorec
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_distribfft(distribfft_arg)

!Arguments ------------------------------------
 type(distribfft_type), intent(inout) :: distribfft_arg

! ***********************************************************************

 DBG_ENTER("COLL")

 distribfft_arg%n2_coarse=0
 distribfft_arg%n2_fine  =0

 if (allocated(distribfft_arg%tab_fftwf2_distrib)) then
   ABI_DEALLOCATE(distribfft_arg%tab_fftwf2_distrib)
 end if

 if (allocated(distribfft_arg%tab_fftdp2_distrib)) then
   ABI_DEALLOCATE(distribfft_arg%tab_fftdp2_distrib)
 end if
 if (allocated(distribfft_arg%tab_fftdp3_distrib)) then
   ABI_DEALLOCATE(distribfft_arg%tab_fftdp3_distrib)
 end if
 if (allocated(distribfft_arg%tab_fftwf2dg_distrib)) then
  ABI_DEALLOCATE(distribfft_arg%tab_fftwf2dg_distrib)
 end if
 if (allocated(distribfft_arg%tab_fftdp2dg_distrib)) then
  ABI_DEALLOCATE(distribfft_arg%tab_fftdp2dg_distrib)
 end if
 if (allocated(distribfft_arg%tab_fftdp3dg_distrib)) then
  ABI_DEALLOCATE(distribfft_arg%tab_fftdp3dg_distrib)
 end if

 if (allocated(distribfft_arg%tab_fftwf2_local)) then
  ABI_DEALLOCATE(distribfft_arg%tab_fftwf2_local)
 end if

 if (allocated(distribfft_arg%tab_fftdp2_local)) then
  ABI_DEALLOCATE(distribfft_arg%tab_fftdp2_local)
 end if
 if (allocated(distribfft_arg%tab_fftdp3_local)) then
  ABI_DEALLOCATE(distribfft_arg%tab_fftdp3_local)
 end if
 if (allocated(distribfft_arg%tab_fftwf2dg_local)) then
  ABI_DEALLOCATE(distribfft_arg%tab_fftwf2dg_local)
 end if
 if (allocated(distribfft_arg%tab_fftdp2dg_local)) then
   ABI_DEALLOCATE(distribfft_arg%tab_fftdp2dg_local)
 end if
 if (allocated(distribfft_arg%tab_fftdp3dg_local)) then
   ABI_DEALLOCATE(distribfft_arg%tab_fftdp3dg_local)
 end if

 DBG_EXIT("COLL")

end subroutine destroy_distribfft
!!***

!===========================================================

!!****f* m_distribfft/copy_distribfft
!! NAME
!!  copy_distribfft
!!
!! FUNCTION
!!  Cleans-up the mpi information for FFT distribution
!!  (mostly deallocate parts distribfft(:) ).
!!
!! SIDE EFFECTS
!!  mpi_enreg=information about MPI parallelization
!!
!! PARENTS
!!      m_mpinfo
!!
!! CHILDREN
!!
!! SOURCE

subroutine copy_distribfft(distribfft_src, distribfft_dst)

!Arguments ------------------------------------
 type(distribfft_type),intent(in)   :: distribfft_src
 type(distribfft_type),intent(out) :: distribfft_dst

! ***********************************************************************

 DBG_ENTER("COLL")

 distribfft_dst%n2_coarse=  distribfft_src%n2_coarse
 distribfft_dst%n2_fine  =  distribfft_src%n2_fine

 if (allocated(distribfft_src%tab_fftwf2_distrib)) then
   ABI_ALLOCATE(distribfft_dst%tab_fftwf2_distrib,(size(distribfft_src%tab_fftwf2_distrib)))
   distribfft_dst%tab_fftwf2_distrib=distribfft_src%tab_fftwf2_distrib
 end if

 if (allocated(distribfft_src%tab_fftdp2_distrib)) then
  ABI_ALLOCATE(distribfft_dst%tab_fftdp2_distrib,(size(distribfft_src%tab_fftdp2_distrib)))
  distribfft_dst%tab_fftdp2_distrib=distribfft_src%tab_fftdp2_distrib
 end if

 if (allocated(distribfft_src%tab_fftdp3_distrib)) then
  ABI_ALLOCATE(distribfft_dst%tab_fftdp3_distrib,(size(distribfft_src%tab_fftdp3_distrib)))
  distribfft_dst%tab_fftdp3_distrib=distribfft_src%tab_fftdp3_distrib
 end if

 if (allocated(distribfft_src%tab_fftwf2dg_distrib)) then
   ABI_ALLOCATE(distribfft_dst%tab_fftwf2dg_distrib,(size(distribfft_src%tab_fftwf2dg_distrib)))
   distribfft_dst%tab_fftwf2dg_distrib=distribfft_src%tab_fftwf2dg_distrib
 end if

 if (allocated(distribfft_src%tab_fftdp2dg_distrib)) then
   ABI_ALLOCATE(distribfft_dst%tab_fftdp2dg_distrib,(size(distribfft_src%tab_fftdp2dg_distrib)))
   distribfft_dst%tab_fftdp2dg_distrib=distribfft_src%tab_fftdp2dg_distrib
 end if

 if (allocated(distribfft_src%tab_fftdp3dg_distrib)) then
   ABI_ALLOCATE(distribfft_dst%tab_fftdp3dg_distrib,(size(distribfft_src%tab_fftdp3dg_distrib)))
   distribfft_dst%tab_fftdp3dg_distrib=distribfft_src%tab_fftdp3dg_distrib
 end if

 if (allocated(distribfft_src%tab_fftwf2_local)) then
   ABI_ALLOCATE(distribfft_dst%tab_fftwf2_local,(size(distribfft_src%tab_fftwf2_local)))
   distribfft_dst%tab_fftwf2_local=distribfft_src%tab_fftwf2_local
 end if

 if (allocated(distribfft_src%tab_fftdp2_local)) then
   ABI_ALLOCATE(distribfft_dst%tab_fftdp2_local,(size(distribfft_src%tab_fftdp2_local)))
   distribfft_dst%tab_fftdp2_local=distribfft_src%tab_fftdp2_local
 end if

 if (allocated(distribfft_src%tab_fftdp3_local)) then
   ABI_ALLOCATE(distribfft_dst%tab_fftdp3_local,(size(distribfft_src%tab_fftdp3_local)))
   distribfft_dst%tab_fftdp3_local=distribfft_src%tab_fftdp3_local
 end if

 if (allocated(distribfft_src%tab_fftwf2dg_local)) then
   ABI_ALLOCATE(distribfft_dst%tab_fftwf2dg_local,(size(distribfft_src%tab_fftwf2dg_local)))
   distribfft_dst%tab_fftwf2dg_local=distribfft_src%tab_fftwf2dg_local
 end if

 if (allocated(distribfft_src%tab_fftdp2dg_local)) then
   ABI_ALLOCATE(distribfft_dst%tab_fftdp2dg_local,(size(distribfft_src%tab_fftdp2dg_local)))
   distribfft_dst%tab_fftdp2dg_local=distribfft_src%tab_fftdp2dg_local
 end if

 if (allocated(distribfft_src%tab_fftdp3dg_local)) then
   ABI_ALLOCATE(distribfft_dst%tab_fftdp3dg_local,(size(distribfft_src%tab_fftdp3dg_local)))
   distribfft_dst%tab_fftdp3dg_local=distribfft_src%tab_fftdp3dg_local
 end if

 DBG_EXIT("COLL")

end subroutine copy_distribfft
!!***

END MODULE m_distribfft
!!***
