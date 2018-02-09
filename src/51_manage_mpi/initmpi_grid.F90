!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_grid
!! NAME
!!  initmpi_grid
!!
!! FUNCTION
!!  Initializes the MPI information for the grid:
!!    * 2D if parallization KPT/FFT (!paral_kgb & MPI)
!!    * 3D if parallization KPT/FFT/BAND (paral_kgb & MPI)
!!    * 2D in case of an Hartree-Fock calculation
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2018 ABINIT group
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt.
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! PARENTS
!!      mpi_setup
!!
!! CHILDREN
!!      mpi_cart_coords,mpi_cart_create,mpi_cart_sub,mpi_comm_rank,wrtout
!!      xmpi_abort,xmpi_comm_free
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine initmpi_grid(mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

#if defined HAVE_MPI2 && !defined FC_G95
 use mpi
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initmpi_grid'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

#if defined HAVE_MPI1 || (defined HAVE_MPI && defined FC_G95)
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(MPI_type),intent(inout) :: mpi_enreg
!Local variables-------------------------------
!scalars
 integer :: nproc,nproc_eff,spacecomm
 character(len=500) :: msg
#if defined HAVE_MPI
 integer :: commcart_4d,dimcart,ierr,me_cart_4d
 integer :: commcart_2d,me_cart_2d
 logical :: reorder
#endif
!arrays
#if defined HAVE_MPI
 integer,allocatable :: coords(:),sizecart(:)
 logical,allocatable :: periode(:), keepdim(:)
#endif

! *********************************************************************

 DBG_ENTER("COLL")

!Select the correct "world" communicator"
 nproc=mpi_enreg%nproc_cell
 if(mpi_enreg%paral_pert==1) nproc=mpi_enreg%nproc_cell
 spacecomm=mpi_enreg%comm_cell

!Fake values for null communicator
 if (nproc==0) then
   mpi_enreg%nproc_fft    = 0
   mpi_enreg%nproc_band   = 0
   mpi_enreg%nproc_hf    = 0
   mpi_enreg%nproc_kpt    = 0
   mpi_enreg%nproc_spinor   = 0
   mpi_enreg%comm_fft            = xmpi_comm_null
   mpi_enreg%comm_band           = xmpi_comm_null
   mpi_enreg%comm_hf             = xmpi_comm_null
   mpi_enreg%comm_kpt            = xmpi_comm_null
   mpi_enreg%comm_kptband        = xmpi_comm_null
   mpi_enreg%comm_spinor         = xmpi_comm_null
   mpi_enreg%comm_bandspinor     = xmpi_comm_null
   mpi_enreg%comm_spinorfft      = xmpi_comm_null
   mpi_enreg%comm_bandfft        = xmpi_comm_null
   mpi_enreg%comm_bandspinorfft  = xmpi_comm_null
   mpi_enreg%bandpp       = 1
   return
 end if

#if defined HAVE_MPI

 if (mpi_enreg%paral_hf==0) then 
   ! either the option Fock exchange is not active or there is no parallelization on Fock exchange calculation. 

   if (mpi_enreg%nproc_spinor>1) mpi_enreg%paral_spinor=1

    !Effective number of processors used for the grid
   nproc_eff=mpi_enreg%nproc_fft*mpi_enreg%nproc_band *mpi_enreg%nproc_kpt*mpi_enreg%nproc_spinor
   if(nproc_eff/=nproc) then
     write(msg,'(4a,5(a,i0))') &
&     '  The number of band*FFT*kpt*spinor processors, npband*npfft*npkpt*npspinor should be',ch10,&
&     '  equal to the total number of processors, nproc.',ch10,&
&     '  However, npband   =',mpi_enreg%nproc_band,&
&     '           npfft    =',mpi_enreg%nproc_fft,&
&     '           npkpt    =',mpi_enreg%nproc_kpt,&
&     '           npspinor =',mpi_enreg%nproc_spinor,&
&     '       and nproc    =',nproc
     MSG_WARNING(msg)
   end if

   !Nothing to do if only 1 proc
   if (nproc_eff==1) return

   ! Initialize the communicator for Hartree-Fock to xmpi_comm_self
   mpi_enreg%me_hf =0
   mpi_enreg%comm_hf=xmpi_comm_self

   if(mpi_enreg%paral_kgb==0) then
     mpi_enreg%me_fft =0
     mpi_enreg%me_band=0
     mpi_enreg%me_kpt =mpi_enreg%me_cell
     mpi_enreg%me_spinor=0
     mpi_enreg%comm_fft=xmpi_comm_self
     mpi_enreg%comm_band=xmpi_comm_self
     mpi_enreg%comm_kpt=mpi_enreg%comm_cell
     mpi_enreg%comm_spinor=xmpi_comm_self
     mpi_enreg%comm_bandspinor=xmpi_comm_self
     mpi_enreg%comm_kptband=mpi_enreg%comm_cell
     mpi_enreg%comm_spinorfft=xmpi_comm_self
     mpi_enreg%comm_bandfft=xmpi_comm_self
     mpi_enreg%comm_bandspinorfft=xmpi_comm_self
   else
     !  CREATE THE 4D GRID
     !  ==================================================

     !  Create the global cartesian 4D- communicator
     !  valgrind claims this is not deallocated in test v5/72
     !  Can someone knowledgable check?
     dimcart=4
     ABI_ALLOCATE(sizecart,(dimcart))
     ABI_ALLOCATE(periode,(dimcart))
!    MT 2012-june: change the order of the indexes; not sure this is efficient
!    (not efficient on TGCC-Curie).
     sizecart(1)=mpi_enreg%nproc_kpt  ! mpi_enreg%nproc_kpt
     sizecart(2)=mpi_enreg%nproc_band ! mpi_enreg%nproc_band
     sizecart(3)=mpi_enreg%nproc_spinor ! mpi_enreg%nproc_spinor
     sizecart(4)=mpi_enreg%nproc_fft  ! mpi_enreg%nproc_fft
     periode(:)=.false.;reorder=.false.
     call MPI_CART_CREATE(spacecomm,dimcart,sizecart,periode,reorder,commcart_4d,ierr)
     ABI_DEALLOCATE(periode)
     ABI_DEALLOCATE(sizecart)

!    Find the index and coordinates of the current processor
     call MPI_COMM_RANK(commcart_4d, me_cart_4d, ierr)
     ABI_ALLOCATE(coords,(dimcart))
     call MPI_CART_COORDS(commcart_4d, me_cart_4d,dimcart,coords,ierr)
     mpi_enreg%me_kpt =coords(1)
     mpi_enreg%me_band=coords(2)
     mpi_enreg%me_spinor=coords(3)
     mpi_enreg%me_fft =coords(4)
     ABI_DEALLOCATE(coords)

     ABI_ALLOCATE(keepdim,(dimcart))

!    Create the communicator for fft distribution
     keepdim(1)=.false.
     keepdim(2)=.false.
     keepdim(3)=.false.
     keepdim(4)=.true.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_fft,ierr)

!    Create the communicator for band distribution
     keepdim(1)=.false.
     keepdim(2)=.true.
     keepdim(3)=.false.
     keepdim(4)=.false.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_band,ierr)

!    Create the communicator for kpt distribution
     keepdim(1)=.true.
     keepdim(2)=.false.
     keepdim(3)=.false.
     keepdim(4)=.false.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_kpt,ierr)

!    Create the communicator for spinor distribution
     keepdim(1)=.false.
     keepdim(2)=.false.
     keepdim(3)=.true.
     keepdim(4)=.false.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_spinor,ierr)

!    Create the communicator for band-spinor distribution
     keepdim(1)=.false.
     keepdim(2)=.true.
     keepdim(3)=.true.
     keepdim(4)=.false.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_bandspinor,ierr)
     if (ierr /= MPI_SUCCESS ) then
       call xmpi_abort(mpi_enreg%comm_world,ierr)
     end if

!    Create the communicator for kpt-band distribution
     keepdim(1)=.true.
     keepdim(2)=.true.
     keepdim(3)=.false.
     keepdim(4)=.false.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_kptband,ierr)

!    Create the communicator for fft-spinor distribution
     keepdim(1)=.false.
     keepdim(2)=.false.
     keepdim(3)=.true.
     keepdim(4)=.true.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_spinorfft,ierr)

!    Create the communicator for fft-band distribution
     keepdim(1)=.false.
     keepdim(2)=.true.
     keepdim(3)=.false.
     keepdim(4)=.true.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_bandfft,ierr)

!    Create the communicator for fft-band-spinor distribution
     keepdim(1)=.false.
     keepdim(2)=.true.
     keepdim(3)=.true.
     keepdim(4)=.true.
     call MPI_CART_SUB(commcart_4d, keepdim, mpi_enreg%comm_bandspinorfft,ierr)

     ABI_DEALLOCATE(keepdim)
     call xmpi_comm_free(commcart_4d)
   end if

!  Write some data
   write(msg,'(a,4i5)') 'npfft, npband, npspinor and npkpt: ',&
&   mpi_enreg%nproc_fft,mpi_enreg%nproc_band, &
&   mpi_enreg%nproc_spinor,mpi_enreg%nproc_kpt
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a,4i5)') 'me_fft, me_band, me_spinor , me_kpt: ',&
&   mpi_enreg%me_fft,mpi_enreg%me_band,&
&   mpi_enreg%me_spinor, mpi_enreg%me_kpt

 else ! paral_hf==1
!* Option Hartree-Fock is active and more than 1 processor is dedicated to the parallelization over occupied states.

!* Initialize the values of fft, band and spinor communicators, as in the case paral_kgb==0.
   mpi_enreg%me_fft =0
   mpi_enreg%me_band=0
   mpi_enreg%me_spinor=0
   mpi_enreg%comm_fft=xmpi_comm_self
   mpi_enreg%comm_band=xmpi_comm_self
   mpi_enreg%comm_spinor=xmpi_comm_self
   mpi_enreg%comm_bandspinor=xmpi_comm_self
   mpi_enreg%comm_kptband=mpi_enreg%comm_cell
   mpi_enreg%comm_spinorfft=xmpi_comm_self
   mpi_enreg%comm_bandfft=xmpi_comm_self
   mpi_enreg%comm_bandspinorfft=xmpi_comm_self

!* Create the global cartesian 2D- communicator
   dimcart=2
   ABI_ALLOCATE(sizecart,(dimcart))
   ABI_ALLOCATE(periode,(dimcart))
   sizecart(1)=mpi_enreg%nproc_kpt  ! mpi_enreg%nproc_kpt
   sizecart(2)=mpi_enreg%nproc_hf   ! mpi_enreg%nproc_hf
   periode(:)=.false.;reorder=.false.
   call MPI_CART_CREATE(spacecomm,dimcart,sizecart,periode,reorder,commcart_2d,ierr)
   ABI_DEALLOCATE(periode)
   ABI_DEALLOCATE(sizecart)

!* Find the index and coordinates of the current processor
   call MPI_COMM_RANK(commcart_2d, me_cart_2d, ierr)
   ABI_ALLOCATE(coords,(dimcart))
   call MPI_CART_COORDS(commcart_2d, me_cart_2d,dimcart,coords,ierr)
   mpi_enreg%me_kpt =coords(1)
   mpi_enreg%me_hf=coords(2)
   ABI_DEALLOCATE(coords)

   ABI_ALLOCATE(keepdim,(dimcart))

!* Create the communicator for kpt distribution
   keepdim(1)=.true.
   keepdim(2)=.false.
   call MPI_CART_SUB(commcart_2d, keepdim, mpi_enreg%comm_kpt,ierr)

!* Create the communicator for hf distribution
   keepdim(1)=.false.
   keepdim(2)=.true.
   call MPI_CART_SUB(commcart_2d, keepdim, mpi_enreg%comm_hf,ierr)

   ABI_DEALLOCATE(keepdim)
   call xmpi_comm_free(commcart_2d)

!* Write some data
   write(msg,'(a,2(1x,i0))') 'nphf and npkpt: ',mpi_enreg%nproc_hf, mpi_enreg%nproc_kpt
   call wrtout(std_out,msg,'COLL')
   write(msg,'(a,2(1x,i0))') 'me_hf, me_kpt: ',mpi_enreg%me_hf, mpi_enreg%me_kpt
 end if

#endif

 DBG_EXIT("COLL")

end subroutine initmpi_grid
!!***
