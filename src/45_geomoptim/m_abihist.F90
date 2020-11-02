!!****m* ABINIT/m_abihist
!! NAME
!! m_abihist
!!
!! FUNCTION
!! This module contains definition the type abihist
!! and its related routines
!!
!! Datatypes:
!!
!! * abihist: Historical record of atomic positions forces and cell parameters
!!
!! Subroutines:
!!
!! * abihist_init
!! * abihist_free
!! * abihist_bcast
!! * abihist_compare
!! * hist2var
!! * var2hist
!! * vel2hist
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (XG, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_abihist

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_nctk
#if defined HAVE_NETCDF
 use netcdf
#endif

 use m_geometry,  only : fcart2fred, xred2xcart

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_abihist/abihist
!! NAME
!! abihist
!!
!! FUNCTION
!! This type has several vectors, and index scalars to store
!! a proper history of previous evaluations of forces and
!! stresses,velocities,positions and energies
!!
!! It contains:
!! * mxhist                  : Maximum size of history
!! * ihist                   : index of history
!! * acell(3,mxhist)         : Acell
!! * rprimd(3,3,mxhist)      : Rprimd
!! * xred(3,natom,mxhist)    : Xred
!! * fcart(3,natom,mxhist)   : Fcart
!! * strten(6,mxhist)        : STRten
!! * vel(3,natom,mxhist)     : Velocities of atoms
!! * vel_cell(3,natom,mxhist): Velocities of cell
!! * etot(mxhist)            : Electronic total Energy
!! * ekin(mxhist)            : Ionic Kinetic Energy
!! * entropy(mxhist)         : Entropy
!! * time(mxhist)            : Time (or iteration number for GO)
!!
!! NOTES
!! The vectors are not allocated because in some cases
!! not all the vectors are needed, in particular a history
!! of stresses is only needed if optcell/=0, and a history
!! of velocities is needed for ionmov==1
!!
!! Store acell, rprimd and strten even with optcell/=0
!! represent a waste of 12x (dp)[Usually 8 Bytes] per
!! iteration, the reason to store all the records is
!! because some routines (eg bfgs.F90) uses the metric (gmet)
!! for initialize the hessian and we need rprimd for that.
!!
!! SOURCE

 type, public :: abihist

! scalars
! Index of the last element on all records
    integer :: ihist = 0
! Maximun size of the historical records
    integer :: mxhist = 0
! Booleans to know if some arrays are changing
    logical :: isVused  ! If velocities are changing
    logical :: isARused ! If Acell and Rprimd are changing

! arrays
! Vector of (x,y,z)x(mxhist) values of cell dimensions
    real(dp), allocatable :: acell(:,:)
! Vector of (x,y,z)x(x,y,z)x(mxhist) values of primitive vectors
    real(dp), allocatable :: rprimd(:,:,:)
! Vector of (x,y,z)x(natom)x(mxhist) values of reduced coordinates
    real(dp), allocatable :: xred(:,:,:)
! Vector of (x,y,z)x(natom)x(mxhist) values of cartesian forces
    real(dp), allocatable :: fcart(:,:,:)
! Vector of (6)x(mxhist) values of stress tensor
    real(dp), allocatable :: strten(:,:)
! Vector of (x,y,z)x(natom)x(mxhist) values of atomic velocities
    real(dp), allocatable :: vel(:,:,:)
! Vector of (x,y,z)x(x,y,z)x(mxhist) values of cell velocities
    real(dp), allocatable :: vel_cell(:,:,:)
! Vector of (mxhist) values of electronic total energy
    real(dp), allocatable :: etot(:)
! Vector of (mxhist) values of ionic kinetic energy
    real(dp), allocatable :: ekin(:)
! Vector of (mxhist) values of Entropy
    real(dp), allocatable :: entropy(:)
! Vector of (mxhist) values of time (relevant for MD calculations)
    real(dp), allocatable :: time(:)

 end type abihist

 public :: abihist_init             ! Initialize the object
 public :: abihist_free             ! Destroy the object
 public :: abihist_bcast            ! Broadcast the object
 public :: abihist_copy             ! Copy 2 HIST records
 public :: abihist_compare_and_copy ! Compare 2 HIST records; if similar copy
 public :: hist2var                 ! Get xred, acell and rprimd from the history.
 public :: abihist_findIndex            ! Shift history indexes
 public :: var2hist                 ! Append xred, acell and rprimd
 public :: vel2hist                 ! Append velocities and Kinetic Energy
 public :: write_md_hist            ! Write the history into a netcdf file
 public :: write_md_hist_img        ! Write the history into a netcdf file (with images)
 public :: read_md_hist             ! Read the history from a netcdf file
 public :: read_md_hist_img         ! Read the history from a netcdf file (with images)
 public :: get_dims_hist
 public :: read_csts_hist

 interface abihist_init
   module procedure abihist_init_0D
   module procedure abihist_init_1D
 end interface abihist_init
 interface abihist_free
   module procedure abihist_free_0D
   module procedure abihist_free_1D
 end interface abihist_free
 interface abihist_bcast
   module procedure abihist_bcast_0D
   module procedure abihist_bcast_1D
 end interface abihist_bcast

!!***

!----------------------------------------------------------------------

contains  !=============================================================
!!***

!!****f* m_abihist/abihist_init_0D
!! NAME
!! abihist_init_0D
!!
!! FUNCTION
!! Initialize a hist structure - Target: scalar
!!
!! INPUTS
!!
!!  natom = Number of atoms per unitary cell
!!  mxhist = Maximal number of records to store
!!  isVUsed,isARUsed=flags used to initialize hsit structure
!!
!! OUTPUT
!!  hist <type(abihist)> = The hist to initialize
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_abihist
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine abihist_init_0D(hist,natom,mxhist,isVused,isARused)

!Arguments ------------------------------------
 integer,intent(in) :: natom,mxhist
 logical,intent(in) :: isVUsed,isARused
 type(abihist),intent(inout) :: hist

! ***************************************************************

!Initialize indexes
 hist%ihist=1
 hist%mxhist=mxhist

!Initialize flags
 hist%isVused=isVUsed
 hist%isARused=isARUsed


!Allocate all the histories
 ABI_ALLOCATE(hist%acell,(3,mxhist))
 ABI_ALLOCATE(hist%rprimd,(3,3,mxhist))
 ABI_ALLOCATE(hist%xred,(3,natom,mxhist))
 ABI_ALLOCATE(hist%fcart,(3,natom,mxhist))
 ABI_ALLOCATE(hist%strten,(6,mxhist))
 ABI_ALLOCATE(hist%vel,(3,natom,mxhist))
 ABI_ALLOCATE(hist%vel_cell,(3,3,mxhist))
 ABI_ALLOCATE(hist%etot,(mxhist))
 ABI_ALLOCATE(hist%ekin,(mxhist))
 ABI_ALLOCATE(hist%entropy,(mxhist))
 ABI_ALLOCATE(hist%time,(mxhist))

 hist%etot(1)=zero
 hist%ekin(1)=zero
 hist%entropy(1)=zero
 hist%time(1)=zero

 hist%acell(:,1)=zero
 hist%rprimd(:,:,1)=zero
 hist%xred(:,:,1)=zero
 hist%fcart(:,:,1)=zero
 hist%strten(:,1)=zero
 hist%vel(:,:,1)=zero
 hist%vel_cell(:,:,1)=zero

end subroutine abihist_init_0D
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/abihist_init_1D
!! NAME
!! abihist_init_1D
!!
!! FUNCTION
!! Initialize a hist structure - Target: 1D array
!!
!! INPUTS
!!  natom = Number of atoms per unitary cell
!!  mxhist = Maximal number of records to store
!!  isVUsed,isARUsed=flags used to initialize hsit structure
!!
!! OUTPUT
!!  hist(:) <type(abihist)> = The hist to initialize
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine abihist_init_1D(hist,natom,mxhist,isVUsed,isARUsed)

!Arguments ------------------------------------
 integer,intent(in) :: natom,mxhist
 logical,intent(in) :: isVUsed,isARUsed
 type(abihist),intent(inout) :: hist(:)

!Local variables-------------------------------
 integer :: ii

! ***************************************************************

 do ii=1,size(hist)
   call abihist_init_0D(hist(ii),natom,mxhist,isVUsed,isARUsed)
 end do

end subroutine abihist_init_1D
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/abihist_free_0D
!! NAME
!! abihist_free_0D
!!
!! FUNCTION
!! Deallocate all the pointers in a hist structure - Target: scalar
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  hist <type(abihist)> = The hist to deallocate
!!
!! PARENTS
!!      m_abihist
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine abihist_free_0D(hist)

!Arguments ------------------------------------
 type(abihist),intent(inout) :: hist

! ***************************************************************

!Vector of (mxhist) values of cell dimensions
 if (allocated(hist%acell))  then
   ABI_DEALLOCATE(hist%acell)
 end if
!Vector of (mxhist) values of cell primitive translations
 if (allocated(hist%rprimd))  then
   ABI_DEALLOCATE(hist%rprimd)
 end if
!Vector of (mxhist) values of reduced coordinates
 if (allocated(hist%xred))  then
   ABI_DEALLOCATE(hist%xred)
 end if
!Vector of (mxhist) values of cartesian forces
 if (allocated(hist%fcart))  then
   ABI_DEALLOCATE(hist%fcart)
 end if
!Vector of (mxhist) values of stress tensor
 if (allocated(hist%strten))  then
   ABI_DEALLOCATE(hist%strten)
 end if
! Vector of (mxhist) values of atomic velocities
 if (allocated(hist%vel))  then
   ABI_DEALLOCATE(hist%vel)
 end if
! Vector of (mxhist) values of cell velocities
 if (allocated(hist%vel_cell))  then
   ABI_DEALLOCATE(hist%vel_cell)
 end if
! Vector of (mxhist) values of electronic total energy
 if (allocated(hist%etot))  then
   ABI_DEALLOCATE(hist%etot)
 end if
! Vector of (mxhist) values of ionic kinetic energy
 if (allocated(hist%ekin))  then
   ABI_DEALLOCATE(hist%ekin)
 end if
! Vector of (mxhist) values of Entropy
 if (allocated(hist%entropy))  then
   ABI_DEALLOCATE(hist%entropy)
 end if
! Vector of (mxhist) values of time (relevant for MD calculations)
 if (allocated(hist%time))  then
   ABI_DEALLOCATE(hist%time)
 end if

end subroutine abihist_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/abihist_free_1D
!! NAME
!! abihist_free_1D
!!
!! FUNCTION
!! Deallocate all the pointers in a hist structure - Target: 1D array
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  hist(:) <type(abihist)> = The hist to deallocate
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine abihist_free_1D(hist)

!Arguments ------------------------------------
 type(abihist),intent(inout) :: hist(:)

!Local variables-------------------------------
 integer :: ii

! ***************************************************************

 do ii=1,size(hist)
   call abihist_free_0D(hist(ii))
 end do

end subroutine abihist_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/abihist_bcast_0D
!! NAME
!! abihist_bcast_0D
!!
!! FUNCTION
!! Broadcast a hist datastructure (from a root process to all others) - Target: scalar
!!
!! INPUTS
!!  master=ID of the sending node in comm
!!  comm=MPI Communicator
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  hist <type(abihist)> = The hist to broadcast
!!
!! PARENTS
!!      m_abihist
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine abihist_bcast_0D(hist,master,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: master,comm
 type(abihist),intent(inout) :: hist

!Local variables-------------------------------
!scalars
 integer :: bufsize,ierr,indx,nproc,rank
 integer :: sizeA,sizeA1,sizeA2
 integer :: sizeEt,sizeEk,sizeEnt,sizeT
 integer :: sizeR,sizeR1,sizeR2,sizeR3
 integer :: sizeS,sizeS1,sizeS2
 integer :: sizeV,sizeV1,sizeV2,sizeV3
 integer :: sizeVc,sizeVc1,sizeVc2,sizeVc3
 integer :: sizeX,sizeX1,sizeX2,sizeX3
 integer :: sizeF,sizeF1,sizeF2,sizeF3
!arrays
 integer,allocatable :: buffer_i(:)
 real(dp),allocatable :: buffer_r(:)

! ***************************************************************

 ierr=0
 nproc=xmpi_comm_size(comm)
 if (nproc<=1) return

 rank=xmpi_comm_rank(comm)

!=== Broadcast integers and logicals
 ABI_ALLOCATE(buffer_i,(4))
 if (rank==master) then
   buffer_i(1)=hist%ihist
   buffer_i(2)=hist%mxhist
   buffer_i(3)=0;if (hist%isVused)  buffer_i(3)=1
   buffer_i(4)=0;if (hist%isARused) buffer_i(4)=1
 end if
 call xmpi_bcast(buffer_i,master,comm,ierr)
 if (rank/=master) then
   hist%ihist=buffer_i(1)
   hist%mxhist=buffer_i(2)
   hist%isVused  =(buffer_i(3)==1)
   hist%isARused =(buffer_i(4)==1)
 end if
 ABI_DEALLOCATE(buffer_i)

!If history is empty, return
 if (hist%mxhist==0.or.hist%ihist==0) return

!=== Broadcast sizes of arrays
 ABI_ALLOCATE(buffer_i,(23))
 if (rank==master) then
   sizeA1=size(hist%acell,1);sizeA2=size(hist%acell,2)
   sizeEt=size(hist%etot,1);sizeEk=size(hist%ekin,1);
   sizeEnt=size(hist%entropy,1);sizeT=size(hist%time,1)
   sizeR1=size(hist%rprimd,1);sizeR2=size(hist%rprimd,2);sizeR3=size(hist%rprimd,3)
   sizeS1=size(hist%strten,1);sizeS2=size(hist%strten,2)
   sizeV1=size(hist%vel,1);sizeV2=size(hist%vel,2);sizeV3=size(hist%vel,3)
   sizeVc1=size(hist%vel_cell,1);sizeVc2=size(hist%vel_cell,2);sizeVc3=size(hist%vel_cell,3)
   sizeX1=size(hist%xred,1);sizeX2=size(hist%xred,2);sizeX3=size(hist%xred,3)
   sizeF1=size(hist%fcart,1);sizeF2=size(hist%fcart,2);sizeF3=size(hist%fcart,3)
   buffer_i(1)=sizeA1  ;buffer_i(2)=sizeA2
   buffer_i(3)=sizeEt  ;buffer_i(4)=sizeEk
   buffer_i(5)=sizeEnt ;buffer_i(6)=sizeT
   buffer_i(7)=sizeR1  ;buffer_i(8)=sizeR2
   buffer_i(9)=sizeR3  ;buffer_i(10)=sizeS1
   buffer_i(11)=sizeS2 ;buffer_i(12)=sizeV1
   buffer_i(13)=sizeV2 ;buffer_i(14)=sizeV3
   buffer_i(15)=sizeVc1;buffer_i(16)=sizeVc2
   buffer_i(17)=sizeVc3;buffer_i(18)=sizeX1
   buffer_i(19)=sizeX2 ;buffer_i(20)=sizeX3
   buffer_i(21)=sizeF1 ;buffer_i(22)=sizeF2
   buffer_i(23)=sizeF3
 end if
 call xmpi_bcast(buffer_i,master,comm,ierr)

 if (rank/=master) then
   sizeA1 =buffer_i(1) ;sizeA2 =buffer_i(2)
   sizeEt =buffer_i(3) ;sizeEk =buffer_i(4)
   sizeEnt=buffer_i(5);sizeT   =buffer_i(6)
   sizeR1 =buffer_i(7) ;sizeR2 =buffer_i(8)
   sizeR3 =buffer_i(9) ;sizeS1 =buffer_i(10)
   sizeS2 =buffer_i(11);sizeV1 =buffer_i(12)
   sizeV2 =buffer_i(13);sizeV3 =buffer_i(14)
   sizeVc1=buffer_i(15);sizeVc2=buffer_i(16)
   sizeVc3=buffer_i(17);sizeX1 =buffer_i(18)
   sizeX2 =buffer_i(19);sizeX3 =buffer_i(20)
   sizeF1 =buffer_i(21);sizeF2 =buffer_i(22)
   sizeF3 =buffer_i(23)
 end if
 ABI_DEALLOCATE(buffer_i)

!=== Broadcast reals
 sizeA=sizeA1*sizeA2;sizeR=sizeR1*sizeR2*sizeR3;sizeS=sizeS1*sizeS2
 sizeV=sizeV1*sizeV2*sizeV3;sizeVc=sizeVc1*sizeVc2*sizeVc3;
 sizeX=sizeX1*sizeX2*sizeX3;sizeF=sizeF1*sizeF2*sizeF3
 bufsize=sizeA+sizeEt+sizeEk+sizeEnt+sizeT+sizeR+sizeS+sizeV+sizeVc+sizeX+sizeF
 ABI_ALLOCATE(buffer_r,(bufsize))
 if (rank==master) then
   indx=0
   buffer_r(indx+1:indx+sizeA)=reshape(hist%acell(1:sizeA1,1:sizeA2),(/sizeA/))
   indx=indx+sizeA
   buffer_r(indx+1:indx+sizeEt)=hist%etot(1:sizeEt)
   indx=indx+sizeEt
   buffer_r(indx+1:indx+sizeEk)=hist%ekin(1:sizeEk)
   indx=indx+sizeEk
   buffer_r(indx+1:indx+sizeEnt)=hist%entropy(1:sizeEnt)
   indx=indx+sizeEnt
   buffer_r(indx+1:indx+sizeT)=hist%time(1:sizeT)
   indx=indx+sizeT
   buffer_r(indx+1:indx+sizeR)=reshape(hist%rprimd(1:sizeR1,1:sizeR2,1:sizeR3),(/sizeR/))
   indx=indx+sizeR
   buffer_r(indx+1:indx+sizeS)=reshape(hist%strten(1:sizeS1,1:sizeS2),(/sizeS/))
   indx=indx+sizeS
   buffer_r(indx+1:indx+sizeV)=reshape(hist%vel(1:sizeV1,1:sizeV2,1:sizeV3),(/sizeV/))
   indx=indx+sizeV
   buffer_r(indx+1:indx+sizeVc)=reshape(hist%vel_cell(1:sizeVc1,1:sizeVc2,1:sizeVc3),(/sizeVc/))
   indx=indx+sizeVc
   buffer_r(indx+1:indx+sizeX)=reshape(hist%xred(1:sizeX1,1:sizeX2,1:sizeX3),(/sizeX/))
   indx=indx+sizeX
   buffer_r(indx+1:indx+sizeF)=reshape(hist%fcart(1:sizeF1,1:sizeF2,1:sizeF3),(/sizeF/))
 else
   call abihist_free(hist)
   ABI_ALLOCATE(hist%acell,(sizeA1,sizeA2))
   ABI_ALLOCATE(hist%etot,(sizeEt))
   ABI_ALLOCATE(hist%ekin,(sizeEk))
   ABI_ALLOCATE(hist%entropy,(sizeEnt))
   ABI_ALLOCATE(hist%time,(sizeT))
   ABI_ALLOCATE(hist%rprimd,(sizeR1,sizeR2,sizeR3))
   ABI_ALLOCATE(hist%strten,(sizeS1,sizeS2))
   ABI_ALLOCATE(hist%vel,(sizeV1,sizeV2,sizeV3))
   ABI_ALLOCATE(hist%vel_cell,(sizeVc1,sizeVc2,sizeVc3))
   ABI_ALLOCATE(hist%xred,(sizeX1,sizeX2,sizeX3))
   ABI_ALLOCATE(hist%fcart,(sizeF1,sizeF2,sizeF3))
 end if
 call xmpi_bcast(buffer_r,master,comm,ierr)

 if (rank/=master) then
   indx=0
   hist%acell(1:sizeA1,1:sizeA2)=reshape(buffer_r(indx+1:indx+sizeA), &
&                                       (/sizeA1,sizeA2/))
   indx=indx+sizeA
   hist%etot(1:sizeEt)=buffer_r(indx+1:indx+sizeEt)
   indx=indx+sizeEt
   hist%ekin(1:sizeEk)=buffer_r(indx+1:indx+sizeEk)
   indx=indx+sizeEk
   hist%entropy(1:sizeEnt)=buffer_r(indx+1:indx+sizeEnt)
   indx=indx+sizeEnt
   hist%time(1:sizeT)=buffer_r(indx+1:indx+sizeT)
   indx=indx+sizeT
   hist%rprimd(1:sizeR1,1:sizeR2,1:sizeR3)=reshape(buffer_r(indx+1:indx+sizeR), &
&                                                 (/sizeR1,sizeR2,sizeR3/))
   indx=indx+sizeR
   hist%strten(1:sizeS1,1:sizeS2)=reshape(buffer_r(indx+1:indx+sizeS), &
&                                        (/sizeS1,sizeS2/))
   indx=indx+sizeS
   hist%vel(1:sizeV1,1:sizeV2,1:sizeV3)=reshape(buffer_r(indx+1:indx+sizeV), &
&                                              (/sizeV1,sizeV2,sizeV3/))
   indx=indx+sizeV
   hist%vel_cell(1:sizeVc1,1:sizeVc2,1:sizeVc3)=reshape(buffer_r(indx+1:indx+sizeVc), &
&                                                      (/sizeVc1,sizeVc2,sizeVc3/))
   indx=indx+sizeVc
   hist%xred(1:sizeX1,1:sizeX2,1:sizeX3)=reshape(buffer_r(indx+1:indx+sizeX), &
&                                               (/sizeX1,sizeX2,sizeX3/))
   indx=indx+sizeX
   hist%fcart(1:sizeF1,1:sizeF2,1:sizeF3)=reshape(buffer_r(indx+1:indx+sizeF), &
&                                                (/sizeF1,sizeF2,sizeF3/))
 end if
 ABI_DEALLOCATE(buffer_r)

end subroutine abihist_bcast_0D
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/abihist_bcast_1D
!! NAME
!! abihist_bcast_1D
!!
!! FUNCTION
!! Broadcast a hist datastructure (from a root process to all others) - Target: 1D array
!!
!! INPUTS
!!  master=ID of the sending node in comm
!!  comm=MPI Communicator
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  hist(:) <type(abihist)> = The hist to broadcast
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine abihist_bcast_1D(hist,master,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: master,comm
 type(abihist),intent(inout) :: hist(:)

!Local variables-------------------------------
 integer :: ii

! ***************************************************************

 do ii=1,size(hist)
   call abihist_bcast_0D(hist(ii),master,comm)
 end do

end subroutine abihist_bcast_1D
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/var2hist
!!
!! NAME
!! var2hist
!!
!! FUNCTION
!! Set the values of the history "hist"
!! with the values of xred, acell and rprimd
!!
!! INPUTS
!! natom = number of atoms
!! xred(3,natom) = reduced dimensionless atomic
!!                           coordinates
!! acell(3)    = length scales of primitive translations (bohr)
!! rprimd(3,3) = dimensionlal real space primitive translations
!!               (bohr)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist<type abihist>=Historical record of positions, forces
!!      |                    acell, stresses, and energies,
!!
!! PARENTS
!!      m_fit_polynomial_coeff,m_generate_training_set,m_gstateimg,m_mover
!!      m_pred_bfgs,m_pred_delocint,m_pred_diisrelax,m_pred_fire,m_pred_hmc
!!      m_pred_isokinetic,m_pred_isothermal,m_pred_langevin,m_pred_lotf
!!      m_pred_moldyn,m_pred_nose,m_pred_srkna14,m_pred_steepdesc
!!      m_pred_velverlet,m_pred_verlet
!!
!! CHILDREN
!!
!! SOURCE

subroutine var2hist(acell,hist,natom,rprimd,xred,zDEBUG)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
 type(abihist),intent(inout) :: hist
 logical,intent(in) :: zDEBUG
!arrays
 real(dp),intent(in) :: acell(3)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: kk

! *************************************************************

 hist%xred(:,:,hist%ihist)=xred(:,:)
 hist%rprimd(:,:,hist%ihist)=rprimd(:,:)
 hist%acell(:,hist%ihist)=acell(:)

 if(zDEBUG)then
   write (std_out,*) 'Atom positions and cell parameters '
   write (std_out,*) 'ihist: ',hist%ihist
   write (std_out,*) 'xred:'
   do kk=1,natom
     write (std_out,*) xred(:,kk)
   end do
   write(std_out,*) 'rprimd:'
   do kk=1,3
     write(std_out,*) rprimd(:,kk)
   end do
   write(std_out,*) 'acell:'
   write(std_out,*) acell(:)
 end if

end subroutine var2hist
!!***

!!****f* m_abihist/abihist_fi
!!
!! NAME
!! abihist_findIndex
!!
!! FUNCTION
!!
!! INPUTS
!! hist<type abihist>=Historical record of positions, forces
!!      |                    acell, stresses, and energies
!! step = value of the needed step
!!
!! OUTPUT
!! index = index of the step in the hist file
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function abihist_findIndex(hist,step) result(index)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: step
 integer :: index
!arrays
 type(abihist),intent(in) :: hist
!Local variables-------------------------------
!scalars
 integer :: ii,mxhist
!arrays
 character(len=500) :: msg
! *************************************************************

 mxhist = hist%mxhist

 if ((mxhist ==1.and.step/=+1).or.&
&    (mxhist /=1.and.abs(step) >=mxhist)) then
   write(msg,'(a,I0,2a)')' The requested step must be less than ',mxhist,ch10,&
&                     'Action: increase the number of history stored in the hist'
   MSG_BUG(msg)
 end if

 ii = hist%ihist + step

 do while (ii > mxhist)
   ii = ii - mxhist
 end do
 do while (ii <= 0)
   ii = ii + mxhist
 end do

 index = ii

end function abihist_findIndex
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/hist2var
!!
!! NAME
!! hist2var
!!
!! FUNCTION
!! Set the values of xred, acell and rprimd
!! with the values that comes from the history "hist"
!!
!! INPUTS
!! natom = Number of atoms
!! hist<type abihist>=Historical record of positions, forces
!!                           acell, stresses, and energies,
!! zDebug = If true some output will be printed
!!
!! OUTPUT
!!  xred(3,natom) = reduced dimensionless atomic coordinates
!!  acell(3)    = length scales of primitive translations (bohr)
!!  rprimd(3,3) = dimensional real space primitive translations (bohr)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_gstateimg,m_mover,m_precpred_1geo,m_pred_bfgs,m_pred_delocint
!!      m_pred_diisrelax,m_pred_fire,m_pred_hmc,m_pred_isokinetic
!!      m_pred_isothermal,m_pred_langevin,m_pred_lotf,m_pred_moldyn,m_pred_nose
!!      m_pred_srkna14,m_pred_steepdesc,m_pred_velverlet,m_pred_verlet
!!
!! CHILDREN
!!
!! SOURCE

subroutine hist2var(acell,hist,natom,rprimd,xred,zDEBUG)

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom
type(abihist),intent(in) :: hist
logical,intent(in) :: zDEBUG
!arrays
real(dp),intent(out) :: acell(3)
real(dp),intent(out) :: rprimd(3,3)
real(dp),intent(out) :: xred(3,natom)

!Local variables-------------------------------
!scalars
integer :: kk

! *************************************************************

 xred  (:,:)=hist%xred(:,:,hist%ihist)
 acell (:  )=hist%acell(:,hist%ihist)
 rprimd(:,:)=hist%rprimd(:,:,hist%ihist)

 if(zDEBUG)then
   write (std_out,*) 'Atom positions and cell parameters '
   write (std_out,*) 'ihist: ',hist%ihist
   write (std_out,*) 'xred:'
   do kk=1,natom
     write (std_out,*) xred(:,kk)
   end do
   write(std_out,*) 'rprimd:'
   do kk=1,3
     write(std_out,*) rprimd(:,kk)
   end do
   write(std_out,*) 'acell:'
   write(std_out,*) acell(:)
 end if

end subroutine hist2var
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/vel2hist
!!
!! NAME
!! vel2hist
!!
!! FUNCTION
!! Set the values of the history "hist" related with velocities, ie
!! The array of velocities and Kinetic Energy
!!
!! INPUTS
!! amass(natom) = mass of the atoms
!! vel(3,natom)= Velocities of the atoms
!! vel_cell(3,3)= Velocities of the cell
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist<type abihist>=Historical record of positions, forces
!!                               stresses, cell and energies,
!!
!! PARENTS
!!      m_gstateimg,m_mover
!!
!! CHILDREN
!!
!! SOURCE

subroutine vel2hist(amass,hist,vel,vel_cell)

!Arguments ------------------------------------
!scalars
type(abihist),intent(inout) :: hist
!arrays
real(dp),intent(in) :: amass(:)
real(dp),intent(in) :: vel(:,:)
real(dp),intent(in) :: vel_cell(:,:)

!Local variables-------------------------------
!scalars
integer :: ii,jj,natom
real(dp) :: ekin

! *************************************************************

 natom=size(vel,2)

 if (hist%isVused) then

!  Store the velocities
   hist%vel(:,:,hist%ihist)=vel(:,:)
   hist%vel_cell(:,:,hist%ihist)=vel_cell(:,:)

!  Compute the Ionic Kinetic energy
   ekin=zero
   do ii=1,natom
     do jj=1,3
       ekin=ekin+half*amass(ii)*vel(jj,ii)**2
     end do
   end do

 else
   hist%vel(:,:,hist%ihist)=zero
   hist%vel_cell(:,:,hist%ihist)=zero
   ekin=zero
 end if

!Store the Ionic Kinetic Energy
 hist%ekin(hist%ihist)=ekin

end subroutine vel2hist
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/abihist_copy
!! NAME
!! abihist_copy
!!
!! FUNCTION
!! Copy one HIST record in another
!!
!! INPUTS
!!  hist_in <type(abihist)>
!!
!! OUTPUT
!! SIDE EFFECTS
!!  hist_out <type(abihist)>
!!
!! PARENTS
!!      m_effective_potential_file,m_gstateimg
!!
!! CHILDREN
!!
!! SOURCE

subroutine abihist_copy(hist_in,hist_out)

!Arguments ------------------------------------
!scalars
type(abihist),intent(in) :: hist_in
type(abihist),intent(inout) :: hist_out

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! ***************************************************************

!Check
 if (size(hist_in%xred,2)/=size(hist_out%xred,2)) then
   msg='Incompatible sizes for hist_in and hist_out!'
   MSG_BUG(msg)
 end if

!Copy scalars (except ihist and mxhist)
 hist_out%isVused  =hist_in%isVused
 hist_out%isARused =hist_in%isARused

!Copy arrays
 hist_out%acell(:,hist_out%ihist)     = hist_in%acell(:,hist_in%ihist)
 hist_out%rprimd(:,:,hist_out%ihist)  = hist_in%rprimd(:,:,hist_in%ihist)
 hist_out%xred(:,:,hist_out%ihist)    = hist_in%xred(:,:,hist_in%ihist)
 hist_out%fcart(:,:,hist_out%ihist)   = hist_in%fcart(:,:,hist_in%ihist)
 hist_out%strten(:,hist_out%ihist)    = hist_in%strten(:,hist_in%ihist)
 hist_out%vel(:,:,hist_out%ihist)     = hist_in%vel(:,:,hist_in%ihist)
 hist_out%vel_cell(:,:,hist_out%ihist)= hist_in%vel_cell(:,:,hist_in%ihist)
 hist_out%etot(hist_out%ihist)        = hist_in%etot(hist_in%ihist)
 hist_out%ekin(hist_out%ihist)        = hist_in%ekin(hist_in%ihist)
 hist_out%entropy(hist_out%ihist)     = hist_in%entropy(hist_in%ihist)
 hist_out%time(hist_out%ihist)        = hist_in%time(hist_in%ihist)

end subroutine abihist_copy
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/abihist_compare_and_copy
!! NAME
!! abihist_compare
!!
!! FUNCTION
!! Compare 2 HIST records
!!
!! INPUTS
!!  hist_in <type(abihist)>
!!  tolerance
!!  store_all = flag to know if we need to increment ihist (store all the history)
!!              or just call shift (store just the last step)
!!
!! OUTPUT
!!  similar= 1 the records are consistent
!!           0 the records are not consistent
!!
!! SIDE EFFECTS
!!  hist_out <type(abihist)>
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!
!! SOURCE

subroutine abihist_compare_and_copy(hist_in,hist_out,natom,similar,tolerance,store_all)

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom
integer,intent(out) :: similar
real(dp),intent(in) :: tolerance
type(abihist),intent(in) :: hist_in
type(abihist),intent(inout) :: hist_out
logical,intent(in) :: store_all
!Local variables-------------------------------
!scalars
integer :: kk,jj
real(dp) :: maxdiff,diff
real(dp) :: x,y
!array
character(len= 500) :: msg
! ***************************************************************

 ABI_UNUSED(store_all)

 similar=1

 write(msg,'(a,I0,4a)')  'Using values from history, iteration:',hist_in%ihist,ch10,&
&                     'Differences between present history and values stored',ch10,&
&                     'on the previous history.(Relative difference)'
 call wrtout(std_out,msg,'COLL')

 x=hist_out%xred(1,1,hist_out%ihist)
 y=hist_in%xred(1,1,hist_in%ihist)
 maxdiff=2*abs(x-y)/(abs(x)+abs(y))
 do kk=1,natom
   do jj=1,3
     x=hist_out%xred(jj,kk,hist_out%ihist)
     y=hist_in%xred(jj,kk,hist_in%ihist)
     diff=2*abs(x-y)/(abs(x)+abs(y))
     if (diff>maxdiff) maxdiff=diff
   end do
 end do
 write(msg,'(a,e12.5)') 'xred:     ',maxdiff
 call wrtout(std_out,msg,'COLL')

 if (maxdiff>tolerance) similar=0


 x=hist_out%rprimd(1,1,hist_out%ihist)
 y=hist_in%rprimd(1,1,hist_in%ihist)
 maxdiff=2*abs(x-y)/(abs(x)+abs(y))
 do kk=1,3
   do jj=1,3
     x=hist_out%rprimd(jj,kk,hist_out%ihist)
     y=hist_in%rprimd(jj,kk,hist_in%ihist)
     diff=2*abs(x-y)/(abs(x)+abs(y))
     if (diff>maxdiff) maxdiff=diff
   end do
 end do
 write(msg,'(a,e12.5)') 'rprimd:   ',maxdiff
 call wrtout(std_out,msg,'COLL')
 if (maxdiff>tolerance) similar=0


 x=hist_out%acell(1,hist_out%ihist)
 y=hist_in%acell(1,hist_in%ihist)
 maxdiff=2*abs(x-y)/(abs(x)+abs(y))
 do kk=1,3
   x=hist_out%acell(kk,hist_out%ihist)
   y=hist_in%acell(kk,hist_in%ihist)
   diff=2*abs(x-y)/(abs(x)+abs(y))
   if (diff>maxdiff) maxdiff=diff
 end do
 write(msg,'(a,e12.5)') 'acell:    ',maxdiff
 call wrtout(std_out,msg,'COLL')
 if (maxdiff>tolerance) similar=0

 if (similar==1) then
   hist_out%acell(:,hist_out%ihist)     =hist_in%acell(:,hist_in%ihist)
   hist_out%rprimd(:,:,hist_out%ihist)  =hist_in%rprimd(:,:,hist_in%ihist)
   hist_out%xred(:,:,hist_out%ihist)    =hist_in%xred(:,:,hist_in%ihist)
   hist_out%fcart(:,:,hist_out%ihist)   =hist_in%fcart(:,:,hist_in%ihist)
   hist_out%strten(:,hist_out%ihist)    =hist_in%strten(:,hist_in%ihist)
   hist_out%vel(:,:,hist_out%ihist)     =hist_in%vel(:,:,hist_in%ihist)
   hist_out%vel_cell(:,:,hist_out%ihist)=hist_in%vel_cell(:,:,hist_in%ihist)
   hist_out%etot(hist_out%ihist)        =hist_in%etot(hist_in%ihist)
   hist_out%ekin(hist_out%ihist)        =hist_in%ekin(hist_in%ihist)
   hist_out%entropy(hist_out%ihist)     =hist_in%entropy(hist_in%ihist)
   hist_out%time(hist_out%ihist)        =hist_in%time(hist_in%ihist)
 end if

end subroutine abihist_compare_and_copy
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/write_md_hist
!!
!! NAME
!! write_md_hist
!!
!! FUNCTION
!! Write the history file into a netcdf file
!! This version is not compatible with multiple images of the
!!
!! INPUTS
!!  filname=filename of the file where the history will be stored
!!  hist<type abihist>=Historical record of positions, forces, stresses,
!!                        cell dims and energies,
!!  ifirst=1 if first access to the file
!!  itime = index of the step in the hist file
!!  natom=Number of atoms.
!!  nctime=NetCdf TIME between output of molecular dynamics information
!!  ntypat=Number of type of atoms.
!!  typat(natom)=Type of each natom
!!   amu(ntypat)=mass of the atoms (atomic mass unit)
!!  znucl(:)=Nuclear charge for each type of pseudopotential.
!!           WARNING: alchemical mixing is not supported. We assume npsp == ntypat
!!  dtion=time step for Molecular Dynamics
!!
!! TODO
!!  Have you ever heard about ETSF-IO specifications?
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      m_fit_polynomial_coeff,m_mover
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_md_hist(hist,filename,ifirst,itime,natom,nctime,ntypat,&
&                        typat,amu,znucl,dtion,mdtemp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ifirst,itime,natom,nctime,ntypat
 real(dp),intent(in) :: dtion
 character(len=*),intent(in) :: filename
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: amu(ntypat),znucl(:),mdtemp(2)
 type(abihist),intent(inout),target :: hist

!Local variables-------------------------------
#if defined HAVE_NETCDF
!scalars
 integer :: itime_file,ncerr,ncid,npsp
 integer :: xcart_id,xred_id,fcart_id,fred_id
 integer :: vel_id,vel_cell_id,etotal_id,acell_id,rprimd_id,strten_id
 integer :: ekin_id,entropy_id,mdtime_id
 logical :: has_nimage=.false.,need_to_write
 integer, parameter :: imgmov=0
!arrays
#endif

! *************************************************************************

#if defined HAVE_NETCDF

 need_to_write = .FALSE.
 if(nctime==0 .or. ifirst==1) need_to_write = .TRUE.
 if (itime > nctime .and. nctime /= 0) then
     if (mod(itime,nctime) == 0) need_to_write = .TRUE.
 end if
!Return if we don't need to write the HIST file at this step
 if (.not. need_to_write) return

 if (ifirst==1) then
!##### First access: Create NetCDF file and write defs

   write(std_out,*) 'Write iteration in HIST netCDF file (also create it)'
   npsp=size(znucl)

!  Create netCDF file
   ncerr = nf90_create(path=trim(filename),cmode=NF90_CLOBBER,ncid=ncid)
   NCF_CHECK_MSG(ncerr," create netcdf history file")

!  Define all dims and vars
   call def_file_hist(ncid,natom,1,ntypat,npsp,has_nimage)

!  Write variables that do not change
!  (they are not read in a hist structure).
   call write_csts_hist(ncid,dtion,imgmov,typat,znucl,amu,mdtemp)

!  Compute the itime for the hist file
   itime_file = 1
 else
!##### itime>2 access: just open NetCDF file

   if(need_to_write) then

     write(std_out,*) 'Write iteration in HIST netCDF file'

!    Open netCDF file
     ncerr = nf90_open(path=trim(filename),mode=NF90_WRITE, ncid=ncid)
     NCF_CHECK_MSG(ncerr," open netcdf history file")

!    Compute the itime for the hist file
     itime_file = itime
     if(nctime > 0) itime_file = int(anint(real(itime / nctime,sp)))

   end if
 endif

 if(need_to_write) then
  !##### Write variables into the dataset
  !Get the IDs
   call get_varid_hist(ncid,xcart_id,xred_id,fcart_id,fred_id,vel_id,vel_cell_id,&
&       rprimd_id,acell_id,strten_id,etotal_id,ekin_id,entropy_id,mdtime_id,has_nimage)
!Write
   call write_vars_hist(ncid,hist,natom,has_nimage,1,itime_file,&
&       xcart_id,xred_id,fcart_id,fred_id,vel_id,vel_cell_id,&
&       rprimd_id,acell_id,strten_id,etotal_id,ekin_id,entropy_id,mdtime_id)

!##### Close the file
   ncerr = nf90_close(ncid)
   NCF_CHECK_MSG(ncerr," close netcdf history file")
 end if
#endif

end subroutine write_md_hist
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/write_md_hist_img
!!
!! NAME
!! write_md_hist_img
!!
!! FUNCTION
!! Write the history file into a netcdf file
!! This version is compatible with multiple images of the
!!
!! INPUTS
!!  filname= filename of the file where the history will be stored
!!  hist(:)<type abihist>= Historical record of positions, forces, stresses,
!!                         cell dims and energies,
!!    Size(hist) is equal to a number of images to be written
!!  ifirst= 1 if first access to the file
!!  itime = index of the step in the hist file
!!  natom= Number of atoms.
!!  ntypat= Number of type of atoms.
!!  typat(natom)= Type of each natom
!!  amu(ntypat)= Mass of the atoms (atomic mass unit)
!!  znucl(:)= Nuclear charge for each type of pseudopotential.
!!           WARNING: alchemical mixing is not supported. We assume npsp == ntypat
!!  dtion= Time step for Molecular Dynamics
!!  [nimage]= Total number of images of the cell
!!  [comm_img]= MPI communicator over images of the cell
!!  [imgtab(:)]= In case of multiple images, indexes of images to be read
!!               Default is 1,2,3,...
!!               Size must be equal to size(hist)
!!
!! TODO
!!  Have you ever heard about ETSF-IO specifications?
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      m_gstateimg
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_md_hist_img(hist,filename,ifirst,itime,natom,ntypat,&
&                            typat,amu,znucl,dtion,&
&                            nimage,imgmov,mdtemp,comm_img,imgtab) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ifirst,itime,natom,ntypat
 integer,intent(in),optional :: nimage,imgmov,comm_img
 real(dp),intent(in) :: dtion
 character(len=*),intent(in) :: filename
!arrays
 integer,intent(in) :: typat(natom)
 integer,intent(in),optional :: imgtab(:)
 real(dp),intent(in) :: amu(ntypat),znucl(:),mdtemp(2)
 type(abihist),intent(inout),target :: hist(:)

!Local variables-------------------------------
#if defined HAVE_NETCDF
!scalars
 integer :: ii,iimage,iimg,me_img,my_comm_img,my_nimage,ncerr
 integer :: ncid,nimage_,nproc_img,npsp,imgmov_
 integer :: xcart_id,xred_id,fcart_id,fred_id
 integer :: vel_id,vel_cell_id,etotal_id
 integer :: acell_id,rprimd_id,strten_id
 integer :: ekin_id,entropy_id,mdtime_id
 logical :: has_nimage, has_imgmov
 character(len=500) :: msg
 type(abihist),pointer :: hist_
!arrays
 integer,allocatable :: my_imgtab(:)
#endif

! *************************************************************************

#if defined HAVE_NETCDF

!Manage multiple images of the cell
 has_nimage=present(nimage)
 has_imgmov=present(imgmov)
 nimage_=merge(nimage,1,has_nimage)
 imgmov_=merge(imgmov,0,has_imgmov)
 my_nimage=size(hist) ; if (my_nimage==0) return
 my_comm_img=xmpi_comm_self;if(present(comm_img)) my_comm_img=comm_img
 nproc_img=xmpi_comm_size(my_comm_img)
 me_img=xmpi_comm_rank(my_comm_img)
 ABI_ALLOCATE(my_imgtab,(my_nimage))
 if (present(imgtab)) then
  if (size(my_imgtab)/=my_nimage) then
    msg='Inconsistency between hist and imgtab!'
    MSG_BUG(msg)
  end if
  my_imgtab(:)=imgtab(:)
 else
   my_imgtab(:)=(/(iimage,iimage=1,my_nimage)/)
 end if

!Has to access the HIST file sequentially, proc by proc
 do ii=0,nproc_img-1
   call xmpi_barrier(my_comm_img)
   if (me_img==ii) then

!    ##### First access: Create NetCDF file and write defs
     if (ifirst==1.and.me_img==0) then
       npsp=size(znucl)
!      Create netCDF file
       ncerr = nf90_create(path=trim(filename),cmode=NF90_CLOBBER,ncid=ncid)
       NCF_CHECK_MSG(ncerr," create netcdf history file")
!      Define all dims and vars
       call def_file_hist(ncid,natom,nimage_,ntypat,npsp,has_nimage)
!      Write variables that do not change
!      (they are not read in a hist structure).
       call write_csts_hist(ncid,dtion,imgmov_,typat,znucl,amu,mdtemp)
     end if

!    ##### itime>2 access: just open NetCDF file
     if (ifirst/=1.or.me_img/=0) then
!      Open netCDF file
       ncerr = nf90_open(path=trim(filename),mode=NF90_WRITE, ncid=ncid)
       NCF_CHECK_MSG(ncerr," open netcdf history file")
     end if

     write(std_out,*) 'Write iteration in HIST netCDF file'

!    ##### Write variables into the dataset (loop over images)
!    Get the IDs
     call get_varid_hist(ncid,xcart_id,xred_id,fcart_id,fred_id,vel_id,vel_cell_id,&
&         rprimd_id,acell_id,strten_id,etotal_id,ekin_id,entropy_id,mdtime_id,has_nimage)

!    Write
     do iimage=1,my_nimage
       iimg=my_imgtab(iimage)
       hist_ => hist(iimage)
       call write_vars_hist(ncid,hist_,natom,has_nimage,iimg,itime,&
&           xcart_id,xred_id,fcart_id,fred_id,vel_id,vel_cell_id,&
&           rprimd_id,acell_id,strten_id,etotal_id,ekin_id,entropy_id,mdtime_id)
     end do

!    ##### Close the file
     ncerr = nf90_close(ncid)
     NCF_CHECK_MSG(ncerr," close netcdf history file")
     ABI_DEALLOCATE(my_imgtab)

!  End loop on MPI processes
   end if
 end do

#endif

end subroutine write_md_hist_img
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/read_md_hist
!!
!! NAME
!! read_md_hist
!!
!! FUNCTION
!! Read the history file from a netcdf file and store it into a hist dataset structure
!! This version is not compatible with multiple images of the simulation cell.
!!
!! INPUTS
!!  filename = Filename of the NetCDF to read
!!  isVUsed,isARUsed=flags used to initialize hist structure
!!
!! OUTPUT
!!  hist<type abihist>=Historical record of positions, forces, stresses,
!!                     cell dims and energies,
!!
!! PARENTS
!!      m_effective_potential_file,m_mover,m_mover_effpot,m_tdep_readwrite
!!
!! CHILDREN
!!
!! SOURCE

subroutine read_md_hist(filename,hist,isVUsed,isARUsed,readOnlyLast)

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: isVUsed,isARUsed,readOnlyLast
 character(len=*),intent(in) :: filename
!arrays
 type(abihist),intent(inout),target :: hist

!Local variables-------------------------------
#if defined HAVE_NETCDF
!scalars
 integer :: ncerr,ncid,nimage,natom,time,start_time, ntypat
 integer :: nimage_id,natom_id,xyz_id,time_id,six_id, ntypat_id
 integer :: xcart_id,xred_id,fcart_id,fred_id,ekin_id,entropy_id
 integer :: mdtime_id,vel_id,vel_cell_id,etotal_id
 integer :: acell_id,rprimd_id,strten_id
 logical :: has_nimage
#endif

! *************************************************************************

#if defined HAVE_NETCDF

 hist%ihist=0 ; hist%mxhist=0

!Open netCDF file
 ncerr=nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=ncid)
 if(ncerr /= NF90_NOERR) then
   write(std_out,*) 'Could no open ',trim(filename),', starting from scratch'
   return
 else
   write(std_out,*) 'Succesfully open ',trim(filename),' for reading'
   write(std_out,*) 'Extracting information from NetCDF file...'
 end if

!Inquire dimensions IDs and lengths
 call get_dims_hist(ncid,natom,ntypat,nimage,time,&
&     natom_id,ntypat_id,nimage_id,time_id,xyz_id,six_id,has_nimage)

!If only the last step is needing (restarxf==-3 for example)
 if(readOnlyLast)then
   start_time = time
   time = 1
 else
   start_time = 1
 end if

!Allocate hist structure
 call abihist_init(hist,natom,time,isVused,isARused)

!Get the ID of a variables from their name
 call get_varid_hist(ncid,xcart_id,xred_id,fcart_id,fred_id,vel_id,vel_cell_id,&
&     rprimd_id,acell_id,strten_id,etotal_id,ekin_id,entropy_id,mdtime_id,has_nimage)

!Read variables from the dataset and write them into hist
 call read_vars_hist(ncid,hist,natom,time,has_nimage,1,start_time,&
&     xred_id,fcart_id,vel_id,vel_cell_id,rprimd_id,acell_id,&
&     strten_id,etotal_id,ekin_id,entropy_id,mdtime_id)

!Close NetCDF file
 ncerr = nf90_close(ncid)
 NCF_CHECK_MSG(ncerr," close netcdf history file")

#endif

end subroutine read_md_hist
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/read_md_hist_img
!!
!! NAME
!! read_md_hist_img
!!
!! FUNCTION
!! Read the history file from a netcdf file and store it into a hist dataset structure
!! This version is compatible with multiple images of the simulation cell.
!!
!! INPUTS
!!  filename = Filename of the NetCDF to read
!!  isVUsed,isARUsed=flags used to initialize hist structure
!!  [imgtab(:)]= In case of multiple images,indexes of images to be read
!!               Default is 1,2,3,...
!!               Size must be equal to size(hist)
!!
!! OUTPUT
!!  hist(:)<type abihist>=Historical record of positions, forces, stresses,
!!                        cell dims and energies,
!!    Size(hist) is equal to a number of images to be read
!!
!! PARENTS
!!      m_gstateimg
!!
!! CHILDREN
!!
!! SOURCE

subroutine read_md_hist_img(filename,hist,isVUsed,isARused,imgtab)

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: isVUsed,isARused
 character(len=*),intent(in) :: filename
!arrays
 integer,intent(in),optional :: imgtab(:)
 type(abihist),intent(inout),target :: hist(:)

!Local variables-------------------------------
#if defined HAVE_NETCDF
!scalars
 integer :: iimage,iimg,my_nimage,ncerr,ncid,nimage,natom,time
 integer :: nimage_id,natom_id,xyz_id,time_id,six_id, ntypat, ntypat_id
 integer :: xcart_id,xred_id,fcart_id,fred_id,ekin_id,entropy_id
 integer :: mdtime_id,vel_id,vel_cell_id,etotal_id
 integer :: acell_id,rprimd_id,strten_id
 logical :: has_nimage
 character(len=500) :: msg
 type(abihist),pointer :: hist_
 integer,allocatable :: my_imgtab(:)
#endif

! *************************************************************************

#if defined HAVE_NETCDF

 hist%ihist=0 ; hist%mxhist=0

!Open netCDF file
 ncerr=nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=ncid)
 if(ncerr /= NF90_NOERR) then
   write(std_out,*) 'Could no open ',trim(filename),', starting from scratch'
   return
 else
   write(std_out,*) 'Succesfully open ',trim(filename),' for reading'
   write(std_out,*) 'Extracting information from NetCDF file...'
 end if

 !Manage multiple images of the cell
 my_nimage=size(hist)
 if (my_nimage==0) return
 ABI_ALLOCATE(my_imgtab,(my_nimage))
 if (present(imgtab)) then
  if (size(my_imgtab)/=my_nimage) then
    msg='Inconsistency between hist and imgtab!'
    MSG_BUG(msg)
  end if
  my_imgtab(:)=imgtab(:)
 else
   my_imgtab(:)=(/(iimage,iimage=1,my_nimage)/)
 end if

!Inquire dimensions IDs and lengths
 call get_dims_hist(ncid,natom,ntypat,nimage,time,&
&     natom_id,ntypat_id,nimage_id,time_id,xyz_id,six_id,has_nimage)

 if (nimage<maxval(my_imgtab)) then
   msg='Not enough images in the HIST file!'
   MSG_ERROR(msg)
 end if

!Loop over images
 do iimage=1,my_nimage
   iimg=my_imgtab(iimage)
   hist_ => hist(iimage)

!  Allocate hist structure
   call abihist_init(hist_,natom,time,isVused,isARused)

!  Get the ID of a variables from their name
   call get_varid_hist(ncid,xcart_id,xred_id,fcart_id,fred_id,vel_id,vel_cell_id,&
&       rprimd_id,acell_id,strten_id,etotal_id,ekin_id,entropy_id,mdtime_id,has_nimage)

!  Read variables from the dataset and write them into hist
   call read_vars_hist(ncid,hist_,natom,time,has_nimage,iimg,1,&
&       xred_id,fcart_id,vel_id,vel_cell_id,rprimd_id,acell_id,&
&       strten_id,etotal_id,ekin_id,entropy_id,mdtime_id)

 end do

!Close NetCDF file
 ncerr = nf90_close(ncid)
 NCF_CHECK_MSG(ncerr," close netcdf history file")

 ABI_DEALLOCATE(my_imgtab)

#endif

end subroutine read_md_hist_img
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/def_file_hist
!!
!! NAME
!! def_file_hist
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_abihist
!!
!! CHILDREN
!!
!! SOURCE

subroutine def_file_hist(ncid,natom,nimage,ntypat,npsp,has_nimage)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 integer,intent(in) :: natom,nimage,ntypat,npsp
 logical,intent(in) :: has_nimage

!Local variables-------------------------------
#if defined HAVE_NETCDF
!scalars
 integer :: ncerr
 integer :: natom_id,nimage_id,ntypat_id,npsp_id,time_id,xyz_id,six_id
 integer :: xcart_id,xred_id,fcart_id,fred_id,vel_id,vel_cell_id
 integer :: rprimd_id,acell_id,strten_id
 integer :: etotal_id,ekin_id,entropy_id,mdtime_id
 integer :: typat_id,znucl_id,amu_id,dtion_id,imgmov_id, two_id,mdtemp_id
 !character(len=500) :: msg
!arrays
 integer :: dim0(0),dim1(1),dim2(2),dim3(3),dim4(4)
#endif

! *************************************************************************

#if defined HAVE_NETCDF

!1.Define the dimensions

 if (npsp/=ntypat) then
   MSG_WARNING('HIST file does not support alchemical mixing!')
 end if

 ncerr = nf90_def_dim(ncid,"natom",natom,natom_id)
 NCF_CHECK_MSG(ncerr," define dimension natom")

 ncerr = nf90_def_dim(ncid,"ntypat",ntypat,ntypat_id)
 NCF_CHECK_MSG(ncerr," define dimension ntypat")

 if (has_nimage) then
   ncerr = nf90_def_dim(ncid,"nimage",nimage,nimage_id)
   NCF_CHECK_MSG(ncerr," define dimension nimage")
 end if

 ncerr = nf90_def_dim(ncid,"npsp",npsp,npsp_id)
 NCF_CHECK_MSG(ncerr," define dimension npsp")

 ncerr = nf90_def_dim(ncid,"xyz",3,xyz_id)
 NCF_CHECK_MSG(ncerr," define dimension xyz")

 ncerr = nf90_def_dim(ncid,"six",6,six_id)
 NCF_CHECK_MSG(ncerr," define dimension six")

 ncerr = nf90_def_dim(ncid,"time",NF90_UNLIMITED,time_id)
 NCF_CHECK_MSG(ncerr," define dimension time")

 ncerr = nf90_def_dim(ncid,"two",2,two_id)
 NCF_CHECK_MSG(ncerr," define dimension two")

!2.Define the constant variables

 dim1=(/natom_id/)
 call ab_define_var(ncid,dim1,typat_id,NF90_DOUBLE,&
&  "typat","types of atoms","dimensionless" )

 dim1=(/npsp_id/)
 call ab_define_var(ncid,dim1,znucl_id,NF90_DOUBLE,&
&  "znucl","atomic charges","atomic units" )

 dim1=(/ntypat_id/)
 call ab_define_var(ncid,dim1,amu_id,NF90_DOUBLE,&
&  "amu","atomic masses","atomic units" )

 call ab_define_var(ncid,dim0,dtion_id,NF90_DOUBLE,&
&  "dtion","time step","atomic units" )

!mdtemp
 dim1=(/two_id/)
 call ab_define_var(ncid,dim1,mdtemp_id,NF90_DOUBLE,&
&  "mdtemp","Molecular Dynamics Thermostat Temperatures","Kelvin" )

!mdtime
 dim1=(/time_id/)
 call ab_define_var(ncid,dim1,mdtime_id,NF90_DOUBLE,&
& "mdtime","Molecular Dynamics or Relaxation TIME","hbar/Ha" )

!3.Define the evolving variables

!xcart,xred,fcart,fred,vel
 if (has_nimage) then
   call ab_define_var(ncid,dim0,imgmov_id,NF90_INT,&
  &  "imgmov","Image mover","Not relevant" )
   dim4=(/xyz_id,natom_id,nimage_id,time_id/)
   call ab_define_var(ncid,dim4,xcart_id,NF90_DOUBLE,&
&   "xcart","vectors (X) of atom positions in CARTesian coordinates","bohr" )
   call ab_define_var(ncid,dim4,xred_id,NF90_DOUBLE,&
&   "xred","vectors (X) of atom positions in REDuced coordinates","dimensionless" )
   call ab_define_var(ncid,dim4,fcart_id,NF90_DOUBLE,&
&   "fcart","atom Forces in CARTesian coordinates","Ha/bohr" )
   call ab_define_var(ncid,dim4,fred_id,NF90_DOUBLE,&
&   "fred","atom Forces in REDuced coordinates","dimensionless" )
   call ab_define_var(ncid,dim4,vel_id,NF90_DOUBLE,&
&   "vel","VELocities of atoms","bohr*Ha/hbar" )
 else
   dim3=(/xyz_id,natom_id,time_id/)
   call ab_define_var(ncid,dim3,xcart_id,NF90_DOUBLE,&
&   "xcart","vectors (X) of atom positions in CARTesian coordinates","bohr" )
   call ab_define_var(ncid,dim3,xred_id,NF90_DOUBLE,&
&   "xred","vectors (X) of atom positions in REDuced coordinates","dimensionless" )
   call ab_define_var(ncid,dim3,fcart_id,NF90_DOUBLE,&
&   "fcart","atom Forces in CARTesian coordinates","Ha/bohr" )
   call ab_define_var(ncid,dim3,fred_id,NF90_DOUBLE,&
&   "fred","atom Forces in REDuced coordinates","dimensionless" )
   call ab_define_var(ncid,dim3,vel_id,NF90_DOUBLE,&
&   "vel","VELocities of atoms","bohr*Ha/hbar" )
 end if

!rprimd,vel_cell
 if (has_nimage) then
   dim4=(/xyz_id,xyz_id,nimage_id,time_id/)
   call ab_define_var(ncid,dim4,rprimd_id,NF90_DOUBLE,&
&   "rprimd","Real space PRIMitive translations, Dimensional","bohr" )
   call ab_define_var(ncid,dim4,vel_cell_id,NF90_DOUBLE,&
&   "vel_cell","VELocities of CELl","bohr*Ha/hbar" )
 else
   dim3=(/xyz_id,xyz_id,time_id/)
   call ab_define_var(ncid,dim3,rprimd_id,NF90_DOUBLE,&
&   "rprimd","Real space PRIMitive translations, Dimensional","bohr" )
   call ab_define_var(ncid,dim3,vel_cell_id,NF90_DOUBLE,&
&   "vel_cell","VELocities of cell","bohr*Ha/hbar" )
 end if

!acell
 if (has_nimage) then
   dim3=(/xyz_id,nimage_id,time_id/)
   call ab_define_var(ncid,dim3,acell_id,NF90_DOUBLE,&
&   "acell","CELL lattice vector scaling","bohr" )
 else
   dim2=(/xyz_id,time_id/)
   call ab_define_var(ncid,dim2,acell_id,NF90_DOUBLE,&
&   "acell","CELL lattice vector scaling","bohr" )
 end if

!strten
 if (has_nimage) then
   dim3=(/six_id,nimage_id,time_id/)
   call ab_define_var(ncid,dim3,strten_id,NF90_DOUBLE,&
&   "strten","STRess tensor","Ha/bohr^3" )
 else
   dim2=(/six_id,time_id/)
   call ab_define_var(ncid,dim2,strten_id,NF90_DOUBLE,&
&   "strten","STRess tensor","Ha/bohr^3" )
 end if

!etotal,ekin,entropy
 if (has_nimage) then
   dim2=(/nimage_id,time_id/)
   call ab_define_var(ncid,dim2,etotal_id,NF90_DOUBLE,&
&   "etotal","TOTAL Energy","Ha" )
   call ab_define_var(ncid,dim2,ekin_id,NF90_DOUBLE,&
&   "ekin","Energy KINetic ionic","Ha" )
   call ab_define_var(ncid,dim2,entropy_id,NF90_DOUBLE,&
&   "entropy","Entropy","" )
 else
   dim1=(/time_id/)
   call ab_define_var(ncid,dim1,etotal_id,NF90_DOUBLE,&
&   "etotal","TOTAL Energy","Ha" )
   call ab_define_var(ncid,dim1,ekin_id,NF90_DOUBLE,&
&   "ekin","Energy KINetic ionic","Ha" )
   call ab_define_var(ncid,dim1,entropy_id,NF90_DOUBLE,&
&   "entropy","Entropy","" )
 end if

!4.End define mode

 ncerr = nf90_enddef(ncid)
 NCF_CHECK_MSG(ncerr," end define mode")

#endif

end subroutine def_file_hist
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/get_dims_hist
!!
!! NAME
!! get_dims_hist
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_abihist,m_tdep_readwrite
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_dims_hist(ncid,natom,ntypat,nimage,time,&
&          natom_id,ntypat_id,nimage_id,time_id,xyz_id,six_id,has_nimage)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 integer,intent(out) :: natom,nimage,time,ntypat
 integer,intent(out) :: natom_id,nimage_id,time_id,xyz_id,six_id, ntypat_id
 logical,intent(out) :: has_nimage

!Local variables-------------------------------
#if defined HAVE_NETCDF
!scalars
 integer :: ncerr
 character(len=5) :: char_tmp
#endif

! *************************************************************************

#if defined HAVE_NETCDF
!Inquire dimensions IDs

 ncerr = nf90_inq_dimid(ncid,"natom",natom_id)
 NCF_CHECK_MSG(ncerr," inquire dimension ID for natom")

 ncerr = nf90_inq_dimid(ncid,"npsp",ntypat_id)
 NCF_CHECK_MSG(ncerr," inquire dimension ID for npsp")

 ncerr = nf90_inq_dimid(ncid,"xyz",xyz_id)
 NCF_CHECK_MSG(ncerr," inquire dimension ID for xyz")

 ncerr = nf90_inq_dimid(ncid,"time",time_id)
 NCF_CHECK_MSG(ncerr," inquire dimension ID for time")

 ncerr = nf90_inq_dimid(ncid,"six",six_id)
 NCF_CHECK_MSG(ncerr," inquire dimension ID for six")

 ncerr = nf90_inq_dimid(ncid,"nimage",nimage_id)
 has_nimage=(ncerr==nf90_noerr)

!Inquire dimensions lengths

 if (has_nimage) then
   ncerr = nf90_inquire_dimension(ncid,nimage_id,char_tmp,nimage)
   has_nimage=(ncerr==nf90_noerr)
 end if
 if (.not.has_nimage) nimage=1

 ncerr = nf90_inquire_dimension(ncid,natom_id,char_tmp,natom)
 NCF_CHECK_MSG(ncerr," inquire dimension natom")

 ncerr = nf90_inquire_dimension(ncid,ntypat_id,char_tmp,ntypat)
 NCF_CHECK_MSG(ncerr," inquire dimension ntypat")

 ncerr = nf90_inquire_dimension(ncid,time_id,char_tmp,time)
 NCF_CHECK_MSG(ncerr," inquire dimension time")

#endif

end subroutine get_dims_hist
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/get_varid_hist
!!
!! NAME
!! get_varid_hist
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_abihist
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_varid_hist(ncid,xcart_id,xred_id,fcart_id,fred_id,vel_id,vel_cell_id,&
&          rprimd_id,acell_id,strten_id,etotal_id,ekin_id,entropy_id,mdtime_id,has_nimage)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 integer,intent(out) :: xcart_id,xred_id,fcart_id,fred_id,vel_id
 integer,intent(out) :: vel_cell_id,rprimd_id,acell_id,strten_id
 integer,intent(out) :: etotal_id,ekin_id,entropy_id,mdtime_id
 logical,intent(in)  :: has_nimage
!Local variables-------------------------------
#if defined HAVE_NETCDF
!scalars
 integer :: ncerr
#endif

! *************************************************************************

#if defined HAVE_NETCDF

 ncerr = nf90_inq_varid(ncid, "mdtime", mdtime_id)
 NCF_CHECK_MSG(ncerr," get the id for mdtime")

 ncerr = nf90_inq_varid(ncid, "xcart", xcart_id)
 NCF_CHECK_MSG(ncerr," get the id for xcart")

 ncerr = nf90_inq_varid(ncid, "xred", xred_id)
 NCF_CHECK_MSG(ncerr," get the id for xred")

 ncerr = nf90_inq_varid(ncid, "fcart", fcart_id)
 NCF_CHECK_MSG(ncerr," get the id for fcart")

 ncerr = nf90_inq_varid(ncid, "fred", fred_id)
 NCF_CHECK_MSG(ncerr," get the id for fred")

 ncerr = nf90_inq_varid(ncid, "vel", vel_id)
 NCF_CHECK_MSG(ncerr," get the id for vel")

 ncerr = nf90_inq_varid(ncid, "vel_cell", vel_cell_id)
 if(has_nimage) then
   NCF_CHECK_MSG(ncerr," get the id for vel_cell")
 end if

 ncerr = nf90_inq_varid(ncid, "rprimd", rprimd_id)
 NCF_CHECK_MSG(ncerr," get the id for rprimd")

 ncerr = nf90_inq_varid(ncid, "acell", acell_id)
 NCF_CHECK_MSG(ncerr," get the id for acell")

 ncerr = nf90_inq_varid(ncid, "strten", strten_id)
 NCF_CHECK_MSG(ncerr," get the id for strten")

 ncerr = nf90_inq_varid(ncid, "etotal", etotal_id)
 NCF_CHECK_MSG(ncerr," get the id for etotal")

 ncerr = nf90_inq_varid(ncid, "ekin", ekin_id)
 NCF_CHECK_MSG(ncerr," get the id for ekin")

 ncerr = nf90_inq_varid(ncid, "entropy", entropy_id)
 NCF_CHECK_MSG(ncerr," get the id for entropy")

#endif

end subroutine get_varid_hist
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/read_csts_hist
!!
!! NAME
!! read_csts_hist
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_tdep_readwrite
!!
!! CHILDREN
!!
!! SOURCE

subroutine read_csts_hist(ncid,dtion,typat,znucl,amu)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 real(dp),intent(out) :: dtion
!arrays
 integer,intent(out) :: typat(:)
 real(dp),intent(out) :: amu(:),znucl(:)

!Local variables-------------------------------
#if defined HAVE_NETCDF
!scalars
 integer :: ncerr
 integer :: typat_id,znucl_id,amu_id,dtion_id
#endif

! *************************************************************************

#if defined HAVE_NETCDF

!1.Get the IDs
 ncerr = nf90_inq_varid(ncid, "typat", typat_id)
 NCF_CHECK_MSG(ncerr," get the id for typat")

 ncerr = nf90_inq_varid(ncid, "znucl", znucl_id)
 NCF_CHECK_MSG(ncerr," get the id for znucl")

 ncerr = nf90_inq_varid(ncid, "amu", amu_id)
 NCF_CHECK_MSG(ncerr," get the id for amu")

 ncerr = nf90_inq_varid(ncid, "dtion", dtion_id)
 NCF_CHECK_MSG(ncerr," get the id for dtion")

!2.Write the constants
 ncerr = nf90_get_var(ncid, typat_id, typat)
 NCF_CHECK_MSG(ncerr," get variable typat")

 ncerr = nf90_get_var(ncid, znucl_id, znucl)
 NCF_CHECK_MSG(ncerr," get variable znucl")

 ncerr = nf90_get_var(ncid, amu_id, amu)
 NCF_CHECK_MSG(ncerr," get variable amu")

 ncerr = nf90_get_var(ncid, dtion_id, dtion)
 NCF_CHECK_MSG(ncerr," get variable dtion")

#endif

end subroutine read_csts_hist
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/write_csts_hist
!!
!! NAME
!! write_csts_hist
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_abihist
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_csts_hist(ncid,dtion,imgmov,typat,znucl,amu,mdtemp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 real(dp),intent(in) :: dtion
 integer,intent(in) :: imgmov
!arrays
 integer,intent(in) :: typat(:)
 real(dp),intent(in) :: amu(:),znucl(:), mdtemp(2)

!Local variables-------------------------------
#if defined HAVE_NETCDF
!scalars
 integer :: ncerr
 integer :: typat_id,znucl_id,amu_id,dtion_id, imgmov_id, mdtemp_id
#endif

! *************************************************************************

#if defined HAVE_NETCDF

!1.Get the IDs

 ncerr = nf90_inq_varid(ncid, "typat", typat_id)
 NCF_CHECK_MSG(ncerr," get the id for typat")

 ncerr = nf90_inq_varid(ncid, "znucl", znucl_id)
 NCF_CHECK_MSG(ncerr," get the id for znucl")

 ncerr = nf90_inq_varid(ncid, "amu", amu_id)
 NCF_CHECK_MSG(ncerr," get the id for amu")

 ncerr = nf90_inq_varid(ncid, "dtion", dtion_id)
 NCF_CHECK_MSG(ncerr," get the id for dtion")

 if ( nf90_noerr == nf90_inq_varid(ncid, "imgmov", imgmov_id) ) then
   ncerr = nf90_put_var(ncid, imgmov_id, imgmov)
   NCF_CHECK_MSG(ncerr," write variable imgmov")
 end if

 if ( nf90_noerr == nf90_inq_varid(ncid, "mdtemp", mdtemp_id) ) then
   ncerr = nf90_put_var(ncid, mdtemp_id, mdtemp)
   NCF_CHECK_MSG(ncerr," write variable mdtemp")
 end if

!2.Write the constants

 ncerr = nf90_put_var(ncid, typat_id, typat)
 NCF_CHECK_MSG(ncerr," write variable typat")

 ncerr = nf90_put_var(ncid, znucl_id, znucl)
 NCF_CHECK_MSG(ncerr," write variable znucl")

 ncerr = nf90_put_var(ncid, amu_id, amu)
 NCF_CHECK_MSG(ncerr," write variable amu")

 ncerr = nf90_put_var(ncid, dtion_id, dtion)
 NCF_CHECK_MSG(ncerr," write variable dtion")

#endif

end subroutine write_csts_hist
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/write_vars_hist
!!
!! NAME
!! write_vars_hist
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_abihist
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_vars_hist(ncid,hist,natom,has_nimage,iimg,itime,&
&          xcart_id,xred_id,fcart_id,fred_id,vel_id,vel_cell_id,rprimd_id,&
&          acell_id,strten_id,etotal_id,ekin_id,entropy_id,mdtime_id)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,natom,iimg,itime
 integer,intent(in) :: xcart_id,xred_id,fcart_id,fred_id,vel_id
 integer,intent(in) :: vel_cell_id,rprimd_id,acell_id,strten_id
 integer,intent(in) :: etotal_id,ekin_id,entropy_id,mdtime_id
 logical,intent(in) :: has_nimage
 type(abihist),intent(inout),target :: hist

!Local variables-------------------------------
#if defined HAVE_NETCDF
!scalars
 integer :: ncerr
!arrays
 integer :: count2(2),count3(3),count4(4)
 integer :: start1(1),start2(2),start3(3),start4(4)
 real(dp),allocatable :: conv(:,:)
 real(dp),pointer :: xred(:,:),fcart(:,:),rprimd(:,:),vel(:,:),vel_cell(:,:)
#endif

! *************************************************************************

#if defined HAVE_NETCDF

 xred     => hist%xred(:,:,hist%ihist)
 fcart    => hist%fcart(:,:,hist%ihist)
 vel      => hist%vel(:,:,hist%ihist)
 vel_cell => hist%vel_cell(:,:,hist%ihist)
 rprimd   => hist%rprimd(:,:,hist%ihist)

!Variables not depending on images

!mdtime
 start1=(/itime/)
 ncerr = nf90_put_var(ncid,mdtime_id,hist%time(hist%ihist),start=start1)
 NCF_CHECK_MSG(ncerr," write variable mdtime")

!Variables depending on images

 ABI_ALLOCATE(conv,(3,natom))

!xcart,xred,fcart,fred,vel
 if (has_nimage) then
   start4=(/1,1,iimg,itime/);count4=(/3,natom,1,1/)
   call xred2xcart(natom,rprimd,conv,xred)
   ncerr = nf90_put_var(ncid,xcart_id,conv, start = start4,count = count4)
   NCF_CHECK_MSG(ncerr," write variable xcart")
   ncerr = nf90_put_var(ncid,xred_id,xred, start = start4,count = count4)
   NCF_CHECK_MSG(ncerr," write variable xred")
   ncerr = nf90_put_var(ncid,fcart_id,fcart,start = start4,count = count4)
   NCF_CHECK_MSG(ncerr," write variable fcart")
   call fcart2fred(fcart,conv,rprimd,natom)
   ncerr = nf90_put_var(ncid,fred_id,conv, start = start4,count = count4)
   NCF_CHECK_MSG(ncerr," write variable fred")
   ncerr = nf90_put_var(ncid,vel_id,vel, start = start4,count = count4)
   NCF_CHECK_MSG(ncerr," write variable vel")
 else
   start3=(/1,1,itime/);count3=(/3,natom,1/)
   call xred2xcart(natom,rprimd,conv,xred)
   ncerr = nf90_put_var(ncid,xcart_id,conv, start = start3,count = count3)
   NCF_CHECK_MSG(ncerr," write variable xcart")
   ncerr = nf90_put_var(ncid,xred_id,xred, start = start3,count = count3)
   NCF_CHECK_MSG(ncerr," write variable xred")
   ncerr = nf90_put_var(ncid,fcart_id,fcart,start = start3,count = count3)
   NCF_CHECK_MSG(ncerr," write variable fcart")
   call fcart2fred(fcart,conv,rprimd,natom)
   ncerr = nf90_put_var(ncid,fred_id,conv, start = start3,count = count3)
   NCF_CHECK_MSG(ncerr," write variable fred")
   ncerr = nf90_put_var(ncid,vel_id,vel, start = start3,count = count3)
   NCF_CHECK_MSG(ncerr," write variable vel")
 end if

 ABI_DEALLOCATE(conv)

!rprimd,vel_cell
 if (has_nimage) then
   start4=(/1,1,iimg,itime/);count4=(/3,3,1,1/)
   ncerr = nf90_put_var(ncid,rprimd_id,hist%rprimd(:,:,hist%ihist),&
&                       start = start4,count = count4)
   NCF_CHECK_MSG(ncerr," write variable rprimd")
   ncerr = nf90_put_var(ncid,vel_cell_id,hist%vel_cell(:,:,hist%ihist),&
&                       start = start4,count = count4)
   NCF_CHECK_MSG(ncerr," write variable vel_cell")
 else
   start3=(/1,1,itime/);count3=(/3,3,1/)
   ncerr = nf90_put_var(ncid,rprimd_id,hist%rprimd(:,:,hist%ihist),&
&                       start = start3,count = count3)
   NCF_CHECK_MSG(ncerr," write variable rprimd")
 end if

!acell
 if (has_nimage) then
   start3=(/1,iimg,itime/);count3=(/3,1,1/)
   ncerr = nf90_put_var(ncid,acell_id,hist%acell(:,hist%ihist),&
&                       start = start3,count = count3)
   NCF_CHECK_MSG(ncerr," write variable acell")
 else
   start2=(/1,itime/);count2=(/3,1/)
   ncerr = nf90_put_var(ncid,acell_id,hist%acell(:,hist%ihist),&
&                       start = start2,count = count2)
   NCF_CHECK_MSG(ncerr," write variable acell")
 end if

!strten
 if (has_nimage) then
   start3=(/1,iimg,itime/);count3=(/6,1,1/)
   ncerr = nf90_put_var(ncid,strten_id,hist%strten(:,hist%ihist),&
&                       start = start3,count = count3)
   NCF_CHECK_MSG(ncerr," write variable strten")
 else
   start2=(/1,itime/);count2=(/6,1/)
   ncerr = nf90_put_var(ncid,strten_id,hist%strten(:,hist%ihist),&
&                       start = start2,count = count2)
   NCF_CHECK_MSG(ncerr," write variable strten")
 end if

!etotal,ekin,entropy
 if (has_nimage) then
   start2=(/iimg,itime/)
   ncerr = nf90_put_var(ncid,etotal_id,hist%etot(hist%ihist),start=start2)
   NCF_CHECK_MSG(ncerr," write variable etotal")
   ncerr = nf90_put_var(ncid,ekin_id,hist%ekin(hist%ihist),start=start2)
   NCF_CHECK_MSG(ncerr," write variable ekin")
   ncerr = nf90_put_var(ncid,entropy_id,hist%entropy(hist%ihist),start=start2)
   NCF_CHECK_MSG(ncerr," write variable entropy")
 else
   start1=(/itime/)
   ncerr = nf90_put_var(ncid,etotal_id,hist%etot(hist%ihist),start=start1)
   NCF_CHECK_MSG(ncerr," write variable etotal")
   ncerr = nf90_put_var(ncid,ekin_id,hist%ekin(hist%ihist),start=start1)
   NCF_CHECK_MSG(ncerr," write variable ekin")
   ncerr = nf90_put_var(ncid,entropy_id,hist%entropy(hist%ihist),start=start1)
   NCF_CHECK_MSG(ncerr," write variable entropy")
 end if

#endif

end subroutine write_vars_hist
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/read_vars_hist
!!
!! NAME
!! read_vars_hist
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_abihist
!!
!! CHILDREN
!!
!! SOURCE

subroutine read_vars_hist(ncid,hist,natom,time,has_nimage,iimg,start_time,&
&          xred_id,fcart_id,vel_id,vel_cell_id,rprimd_id,acell_id,&
&          strten_id,etotal_id,ekin_id,entropy_id,mdtime_id)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,natom,time,iimg
 integer,intent(in) :: xred_id,fcart_id,vel_id,vel_cell_id,rprimd_id
 integer,intent(in) :: acell_id,strten_id
 integer,intent(in) :: etotal_id,ekin_id,entropy_id,mdtime_id
 integer,intent(in) :: start_time
 logical,intent(in) :: has_nimage
 type(abihist),intent(inout),target :: hist

!Local variables-------------------------------
#if defined HAVE_NETCDF
!scalars
 integer :: ncerr
!arrays
 integer :: count1(1),count2(2),count3(3),count4(4)
 integer :: start1(1),start2(2),start3(3),start4(4)
#endif

! *************************************************************************

#if defined HAVE_NETCDF

!Variables not depending on imes

!mdtime
 start1=(/start_time/);count1=(/time/)
 ncerr = nf90_get_var(ncid,mdtime_id,hist%time(:),count=count1,start=start1)
 NCF_CHECK_MSG(ncerr," read variable mdtime")

!Variables depending on images

!xred,fcart,vel
 if (has_nimage) then
   start4=(/1,1,iimg,start_time/);count4=(/3,natom,1,time/)
   ncerr = nf90_get_var(ncid,xred_id  ,hist%xred(:,:,:),count=count4,start=start4)
   NCF_CHECK_MSG(ncerr," read variable xred")
   ncerr = nf90_get_var(ncid,fcart_id ,hist%fcart(:,:,:),count=count4,start=start4)
   NCF_CHECK_MSG(ncerr," read variable fcart")
   ncerr = nf90_get_var(ncid,vel_id,hist%vel(:,:,:),count=count4,start=start4)
   NCF_CHECK_MSG(ncerr," read variable vel")
 else
   start3=(/1,1,start_time/);count3=(/3,natom,time/)
   ncerr = nf90_get_var(ncid,xred_id  ,hist%xred(:,:,:),count=count3,start=start3)
   NCF_CHECK_MSG(ncerr," read variable xred")
   ncerr = nf90_get_var(ncid,fcart_id ,hist%fcart(:,:,:),count=count3,start=start3)
   NCF_CHECK_MSG(ncerr," read variable fcart")
   ncerr = nf90_get_var(ncid,vel_id,hist%vel(:,:,:),count=count3,start=start3)
   NCF_CHECK_MSG(ncerr," read variable vel")
 end if

!rprimd,vel_cell
 if (has_nimage) then
   start4=(/1,1,iimg,start_time/);count4=(/3,3,start_time,time/)
   ncerr = nf90_get_var(ncid,rprimd_id,hist%rprimd(:,:,:),count=count4,start=start4)
   NCF_CHECK_MSG(ncerr," read variable rprimd")
   ncerr = nf90_get_var(ncid,vel_cell_id,hist%vel_cell(:,:,:),count=count4,start=start4)
   NCF_CHECK_MSG(ncerr," read variable vel_cell")
 else
   start3=(/1,1,start_time/);count3=(/3,3,time/)
   ncerr = nf90_get_var(ncid,rprimd_id,hist%rprimd(:,:,:),count=count3,start=start3)
   NCF_CHECK_MSG(ncerr," read variable rprimd")
 end if

!acell
 if (has_nimage) then
   start3=(/1,iimg,start_time/);count3=(/3,1,time/)
   ncerr = nf90_get_var(ncid,acell_id,hist%acell(:,:),count=count3,start=start3)
   NCF_CHECK_MSG(ncerr," read variable acell")
 else
   start2=(/1,start_time/);count2=(/3,time/)
   ncerr = nf90_get_var(ncid,acell_id,hist%acell(:,:),count=count2,start=start2)
   NCF_CHECK_MSG(ncerr," read variable acell")
 end if

!strten
 if (has_nimage) then
   start3=(/1,iimg,start_time/);count3=(/6,1,time/)
   ncerr = nf90_get_var(ncid, strten_id,hist%strten(:,:),count=count3,start=start3)
   NCF_CHECK_MSG(ncerr," read variable strten")
 else
   start2=(/1,start_time/);count2=(/6,time/)
   ncerr = nf90_get_var(ncid, strten_id,hist%strten(:,:),count=count2,start=start2)
   NCF_CHECK_MSG(ncerr," read variable strten")
 end if

!etotal,ekin,entropy
 if (has_nimage) then
   start2=(/1,start_time/);count2=(/1,time/)
   ncerr = nf90_get_var(ncid,etotal_id,hist%etot(:),count=count2,start=start2)
   NCF_CHECK_MSG(ncerr," read variable etotal")
   ncerr = nf90_get_var(ncid,ekin_id ,hist%ekin(:),count=count2,start=start2)
   NCF_CHECK_MSG(ncerr," read variable ekin")
   ncerr = nf90_get_var(ncid,entropy_id,hist%entropy(:),count=count2,start=start2)
   NCF_CHECK_MSG(ncerr," read variable entropy")
 else
   start1=(/start_time/);count1=(/time/)
   ncerr = nf90_get_var(ncid,etotal_id,hist%etot(:),count=count1,start=start1)
   NCF_CHECK_MSG(ncerr," read variable etotal")
   ncerr = nf90_get_var(ncid,ekin_id,hist%ekin(:),count=count1,start=start1)
   NCF_CHECK_MSG(ncerr," read variable ekin")
   ncerr = nf90_get_var(ncid,entropy_id,hist%entropy(:),count=count1,start=start1)
   NCF_CHECK_MSG(ncerr," read variable entropy")
 end if

#endif

end subroutine read_vars_hist
!!***

end module m_abihist
