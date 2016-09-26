!{\src2tex{textfont=tt}}
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
!! * abihist_ini
!! * abihist_fin
!! * abihist_bcast
!! * abihist_compare
!! * hist2var
!! * var2hist
!! * vel2hist
!!
!! COPYRIGHT
!! Copyright (C) 2001-2016 ABINIT group (XG, SE)
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
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_nctk
#if defined HAVE_NETCDF
 use netcdf
#endif

 use m_abimover, only : abimover

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
!! * histA(3,mxhist)         : Acell
!! * histE(mxhist)           : Energy
!! * histEk(mxhist)          : Ionic Kinetic Energy
!! * histEnt(mxhist)         : Entropy
!! * histT(mxhist)           : Time (Or iteration number for GO)
!! * histR(3,3,mxhist)       : Rprimd
!! * histS(6,mxhist)         : Strten
!! * histV(3,natom,mxhist)   : Velocity
!! * histXF(3,natom,4,mxhist): Xcart,Xred,Fcart,Fred
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
! Vector of (x,y,z)X(mxhist)
    real(dp), allocatable :: histA(:,:)
! Vector of (mxhist) values of energy
    real(dp), allocatable :: histE(:)
! Vector of (mxhist) values of ionic kinetic energy
    real(dp), allocatable :: histEk(:)
! Vector of (mxhist) values of Entropy
    real(dp), allocatable :: histEnt(:)
! Vector of (mxhist) values of time (relevant for MD calculations)
    real(dp), allocatable :: histT(:)
! Vector of (x,y,z)X(x,y,z)X(mxhist)
    real(dp), allocatable :: histR(:,:,:)
! Vector of (stress [6])X(mxhist)
    real(dp), allocatable :: histS(:,:)
! Vector of (x,y,z)X(natom)X(mxhist) values of velocity
    real(dp), allocatable :: histV(:,:,:)
! Vector of (x,y,z)X(natom)X(xcart,xred,fcart,fred)X(mxhist)
    real(dp), allocatable :: histXF(:,:,:,:)

 end type abihist

 public :: abihist_ini         ! Initialize the object
 public :: abihist_fin         ! Free memory
 public :: abihist_bcast       ! Broadcast the object
 public :: abihist_compare     ! Compare 2 HIST records
 public :: hist2var            ! Get xcart, xred, acell and rprimd from the history.
 public :: var2hist            ! Append xcart, xred, acell and rprimd
 public :: vel2hist            ! Append velocities and Kinetic Energy
 public :: write_md_hist       ! Write the history into a netcdf dataset
 public :: read_md_hist        ! Read the history file from a netcdf file
!!***

!----------------------------------------------------------------------

contains  !=============================================================
!!***

!!****f* m_abihist/abihist_ini
!! NAME
!! abihist_ini
!!
!! FUNCTION
!! Initialize the hist structure for a given number
!! configurations and number of atoms
!!
!! INPUTS
!!  natom = Number of atoms per unitary cell
!!  mxhist = Maximal number of records store  
!!
!! OUTPUT
!!  hist <type(abihist)> = The hist to initialize
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_abiimages,mover
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine abihist_ini(hist,natom,mxhist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abihist_ini'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(abihist) ,intent(inout) :: hist
 integer       ,intent(in)  :: natom
 integer       ,intent(in)  :: mxhist

! ***************************************************************

!Initialize indexes
 hist%ihist=1
 hist%mxhist=mxhist

!Allocate all the histories
 ABI_ALLOCATE(hist%histA,(3,mxhist))
 ABI_ALLOCATE(hist%histE,(mxhist))
 ABI_ALLOCATE(hist%histEk,(mxhist))
 ABI_ALLOCATE(hist%histEnt,(mxhist))
 ABI_ALLOCATE(hist%histT,(mxhist))
 ABI_ALLOCATE(hist%histR,(3,3,mxhist))
 ABI_ALLOCATE(hist%histS,(6,mxhist))
 ABI_ALLOCATE(hist%histV,(3,natom,mxhist))
 ABI_ALLOCATE(hist%histXF,(3,natom,4,mxhist))

end subroutine abihist_ini
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/abihist_fin
!! NAME
!! abihist_fin
!!
!! FUNCTION
!! Deallocate all the pointers in a hist
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  hist <type(abihist)> = The hist to deallocate
!!
!! PARENTS
!!      m_abihist,m_abiimages,mover
!!
!! CHILDREN
!!
!! NOTES
!!  At present 8 variables are present in abihist
!!  if a new variable is added in abihist it should
!!  be added also for deallocate here
!!
!! SOURCE

subroutine abihist_fin(hist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abihist_fin'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(abihist),intent(inout) :: hist

! ***************************************************************

! Vector of (x,y,z)X(mxhist)
 if (allocated(hist%histA))  then
   ABI_DEALLOCATE(hist%histA)
 end if
! Vector of (mxhist) values of energy
 if (allocated(hist%histE))  then
   ABI_DEALLOCATE(hist%histE)
 end if
! Vector of (mxhist) values of ionic kinetic energy
 if (allocated(hist%histEk))  then
   ABI_DEALLOCATE(hist%histEk)
 end if
! Vector of (mxhist) values of Entropy
 if (allocated(hist%histEnt))  then
   ABI_DEALLOCATE(hist%histEnt)
 end if
! Vector of (mxhist) values of time (relevant for MD calculations)
 if (allocated(hist%histT))  then
   ABI_DEALLOCATE(hist%histT)
 end if
! Vector of (x,y,z)X(x,y,z)X(mxhist)
 if (allocated(hist%histR))  then
   ABI_DEALLOCATE(hist%histR)
 end if
! Vector of (stress [6])X(mxhist)
 if (allocated(hist%histS))  then
   ABI_DEALLOCATE(hist%histS)
 end if
! Vector of (x,y,z)X(natom)X(mxhist) values of velocity
 if (allocated(hist%histV))  then
   ABI_DEALLOCATE(hist%histV)
 end if
! Vector of (x,y,z)X(natom)X(xcart,xred,fcart,fred)X(mxhist)
 if (allocated(hist%histXF))  then
   ABI_DEALLOCATE(hist%histXF)
 end if

end subroutine abihist_fin
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/abihist_bcast
!! NAME
!! abihist_bcast
!!
!! FUNCTION
!! Broadcast a hist datastructure (from a root process to all others)
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
!!      mover
!!
!! CHILDREN
!!
!! NOTES
!!  At present 8 variables are present in abihist
!!  if a new variable is added in abihist it should
!!  be added also for broadcast here
!!
!! SOURCE

subroutine abihist_bcast(hist,master,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abihist_bcast'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: master,comm
 type(abihist),intent(inout) :: hist

!Local variables-------------------------------
!scalars
 integer :: bufsize,ierr,indx,nproc,rank
 integer :: sizeA,sizeA1,sizeA2
 integer :: sizeE,sizeEk,sizeEnt,sizeT
 integer :: sizeR,sizeR1,sizeR2,sizeR3
 integer :: sizeS,sizeS1,sizeS2
 integer :: sizeV,sizeV1,sizeV2,sizeV3
 integer :: sizeXF,sizeXF1,sizeXF2,sizeXF3,sizeXF4
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
   hist%isVused =(buffer_i(3)==1)
   hist%isARused=(buffer_i(4)==1)
 end if
 ABI_DEALLOCATE(buffer_i)

!If history is empty, return
 if (hist%mxhist==0.or.hist%ihist==0) return

!=== Broadcast sizes of arrays
ABI_ALLOCATE(buffer_i,(18))
 if (rank==master) then
   sizeA1=size(hist%histA,1);sizeA2=size(hist%histA,2)
   sizeE=size(hist%histE,1);sizeEk=size(hist%histEk,1);
   sizeEnt=size(hist%histEnt,1);sizeT=size(hist%histT,1)
   sizeR1=size(hist%histR,1);sizeR2=size(hist%histR,2);sizeR3=size(hist%histR,3)
   sizeS1=size(hist%histS,1);sizeS2=size(hist%histS,2)
   sizeV1=size(hist%histV,1);sizeV2=size(hist%histV,2);sizeV3=size(hist%histV,3)
   sizeXF1=size(hist%histXF,1);sizeXF2=size(hist%histXF,2);sizeXF3=size(hist%histXF,3)
   sizeXF4=size(hist%histXF,4)
   buffer_i(1)=sizeA1  ;buffer_i(2)=sizeA2
   buffer_i(3)=sizeE   ;buffer_i(4)=sizeEk
   buffer_i(5)=sizeEnt ;buffer_i(6)=sizeT 
   buffer_i(7)=sizeR1  ;buffer_i(8)=sizeR2  
   buffer_i(9)=sizeR3  ;buffer_i(10)=sizeS1 
   buffer_i(11)=sizeS2 ;buffer_i(12)=sizeV1 
   buffer_i(13)=sizeV2 ;buffer_i(14)=sizeV3 
   buffer_i(15)=sizeXF1;buffer_i(16)=sizeXF2
   buffer_i(17)=sizeXF3;buffer_i(18)=sizeXF4
 end if
 call xmpi_bcast(buffer_i,master,comm,ierr)

 if (rank/=master) then
   sizeA1 =buffer_i(1) ;sizeA2 =buffer_i(2)
   sizeE  =buffer_i(3) ;sizeEk =buffer_i(4)
   sizeEnt =buffer_i(5);sizeT  =buffer_i(6) 
   sizeR1 =buffer_i(7) ;sizeR2 =buffer_i(8) 
   sizeR3 =buffer_i(9) ;sizeS1 =buffer_i(10)
   sizeS2 =buffer_i(11);sizeV1 =buffer_i(12)
   sizeV2 =buffer_i(13);sizeV3 =buffer_i(14)
   sizeXF1=buffer_i(15);sizeXF2=buffer_i(16)
   sizeXF3=buffer_i(17);sizeXF4=buffer_i(18)   
   
 end if
 ABI_DEALLOCATE(buffer_i)

!=== Broadcast reals
 sizeA=sizeA1*sizeA2;sizeR=sizeR1*sizeR2*sizeR3;sizeS=sizeS1*sizeS2
 sizeV=sizeV1*sizeV2*sizeV3;sizeXF=sizeXF1*sizeXF2*sizeXF3*sizeXF4
 bufsize=sizeA+sizeE+sizeEk+sizeEnt+sizeT+sizeR+sizeS+sizeV+sizeXF
 ABI_ALLOCATE(buffer_r,(bufsize))
 if (rank==master) then
   indx=0
   buffer_r(indx+1:indx+sizeA)=reshape(hist%histA(1:sizeA1,1:sizeA2),(/sizeA/))
   indx=indx+sizeA
   buffer_r(indx+1:indx+sizeE)=hist%histE(1:sizeE)
   indx=indx+sizeE
   buffer_r(indx+1:indx+sizeEk)=hist%histEk(1:sizeEk)
   indx=indx+sizeEk
   buffer_r(indx+1:indx+sizeEnt)=hist%histEnt(1:sizeEnt)
   indx=indx+sizeEnt
   buffer_r(indx+1:indx+sizeT)=hist%histT(1:sizeT)
   indx=indx+sizeT
   buffer_r(indx+1:indx+sizeR)=reshape(hist%histR(1:sizeR1,1:sizeR2,1:sizeR3),(/sizeR/))
   indx=indx+sizeR
   buffer_r(indx+1:indx+sizeS)=reshape(hist%histS(1:sizeS1,1:sizeS2),(/sizeS/))
   indx=indx+sizeS
   buffer_r(indx+1:indx+sizeV)=reshape(hist%histV(1:sizeV1,1:sizeV2,1:sizeV3),(/sizeV/))
   indx=indx+sizeV
   buffer_r(indx+1:indx+sizeXF)=reshape(hist%histXF(1:sizeXF1,1:sizeXF2,1:sizeXF3,1:sizeXF4),(/sizeXF/))
 else
   call abihist_fin(hist)
   ABI_ALLOCATE(hist%histA,(sizeA1,sizeA2))
   ABI_ALLOCATE(hist%histE,(sizeE))
   ABI_ALLOCATE(hist%histEk,(sizeEk))
   ABI_ALLOCATE(hist%histEnt,(sizeEnt))
   ABI_ALLOCATE(hist%histT,(sizeT))
   ABI_ALLOCATE(hist%histR,(sizeR1,sizeR2,sizeR3))
   ABI_ALLOCATE(hist%histS,(sizeS1,sizeS2))
   ABI_ALLOCATE(hist%histV,(sizeV1,sizeV2,sizeV3))
   ABI_ALLOCATE(hist%histXF,(sizeXF1,sizeXF2,sizeXF3,sizeXF4))
 end if
 call xmpi_bcast(buffer_r,master,comm,ierr)

 if (rank/=master) then
   indx=0
   hist%histA(1:sizeA1,1:sizeA2)=reshape(buffer_r(indx+1:indx+sizeA),(/sizeA1,sizeA2/))
   indx=indx+sizeA
   hist%histE(1:sizeE)=buffer_r(indx+1:indx+sizeE)
   indx=indx+sizeE
   hist%histEk(1:sizeEk)=buffer_r(indx+1:indx+sizeEk)
   indx=indx+sizeEk
   hist%histEnt(1:sizeEnt)=buffer_r(indx+1:indx+sizeEnt)
   indx=indx+sizeEnt
   hist%histT(1:sizeT)=buffer_r(indx+1:indx+sizeT)
   indx=indx+sizeT
   hist%histR(1:sizeR1,1:sizeR2,1:sizeR3)=reshape(buffer_r(indx+1:indx+sizeR),(/sizeR1,sizeR2,sizeR3/))
   indx=indx+sizeR
   hist%histS(1:sizeS1,1:sizeS2)=reshape(buffer_r(indx+1:indx+sizeS),(/sizeS1,sizeS2/))
   indx=indx+sizeS
   hist%histV(1:sizeV1,1:sizeV2,1:sizeV3)=reshape(buffer_r(indx+1:indx+sizeV),(/sizeV1,sizeV2,sizeV3/))
   indx=indx+sizeV
   hist%histXF(1:sizeXF1,1:sizeXF2,1:sizeXF3,1:sizeXF4)=reshape(buffer_r(indx+1:indx+sizeXF), &
&                                                       (/sizeXF1,sizeXF2,sizeXF3,sizeXF4/))
 end if
 ABI_DEALLOCATE(buffer_r)

end subroutine abihist_bcast
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/var2hist
!!
!! NAME
!! var2hist
!!
!! FUNCTION
!! Set the values of the history "hist"
!! with the values of xcart, xred, acell and rprimd
!!
!! INPUTS
!! natom = number of atoms
!! xred(3,natom) = reduced dimensionless atomic
!!                           coordinates
!! xcart(3,natom)= cartesian dimensional atomic
!!                           coordinates (bohr)
!! acell(3)    = length scales of primitive translations (bohr)
!! rprimd(3,3) = dimensionlal real space primitive translations
!!               (bohr)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist<type abihist>=Historical record of positions, forces
!!      |                    acell, stresses, and energies,
!!      |                    contains:
!!      | mxhist:  Maximun number of records
!!      | ihist:   Index of present record
!!      | histA:   Historical record of acell(A) and rprimd(R)
!!      | histE:   Historical record of energy(E)
!!      | histEk:  Historical record of kinetic energy(Ek)
!!      | histEnt: Historical record of Entropy
!!      | histT:   Historical record of time(T)
!!      | histR:   Historical record of rprimd(R)
!!      | histS:   Historical record of strten(S)
!!      | histV:   Historical record of velocity(V)
!!      | histXF:  Historical record of positions(X) and forces(F)
!!
!! PARENTS
!!      m_pred_lotf,mover,pred_bfgs,pred_delocint,pred_diisrelax
!!      pred_isokinetic,pred_isothermal,pred_langevin,pred_moldyn,pred_nose
!!      pred_srkna14,pred_steepdesc,pred_verlet
!!
!! CHILDREN
!!
!! SOURCE

subroutine var2hist(acell,hist,natom,rprim,rprimd,xcart,xred,zDEBUG)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'var2hist'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars

integer,intent(in) :: natom
type(abihist),intent(inout) :: hist
logical,intent(in) :: zDEBUG
!arrays
real(dp),intent(in) :: acell(3)
real(dp),intent(in) :: rprimd(3,3)
real(dp),intent(in) :: rprim(3,3)
real(dp),intent(in) :: xcart(3,natom)
real(dp),intent(in) :: xred(3,natom)

!Local variables-------------------------------
!scalars
integer :: kk

! *************************************************************

 hist%histXF(:,:,1,hist%ihist)=xcart(:,:)
 hist%histXF(:,:,2,hist%ihist)=xred(:,:)
 hist%histA(:,hist%ihist)     =acell(:)
 hist%histR(:,:,hist%ihist)   =rprimd(:,:)

 if(zDEBUG)then
   write (std_out,*) 'Atom positions and cell parameters '
   write (std_out,*) 'ihist: ',hist%ihist
   write (std_out,*) 'xcart:'
   do kk=1,natom
     write (std_out,*) xcart(:,kk)
   end do
   write (std_out,*) 'xred:'
   do kk=1,natom
     write (std_out,*) xred(:,kk)
   end do
   write(std_out,*) 'rprim:'
   do kk=1,3
     write(std_out,*) rprim(:,kk)
   end do
   write(std_out,*) 'rprimd:'
   do kk=1,3
     write(std_out,*) rprimd(:,kk)
   end do
   write(std_out,*) 'acell:'
   write(std_out,*) acell(:)
   call chkrprimd(acell,rprim,rprimd,std_out)
 end if

end subroutine var2hist
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/hist2var
!!
!! NAME
!! hist2var
!!
!! FUNCTION
!! Set the values of xcart, xred, acell and rprimd
!! with the values that comes from the history "hist"
!!
!! INPUTS
!! natom = Number of atoms
!! hist<type abihist>=Historical record of positions, forces
!!                           acell, stresses, and energies,
!! zDebug = If true some output will be printed
!!
!! OUTPUT
!!  xred(3,natom) = reduced dimensionless atomic
!!                           coordinates
!!  xcart(3,natom)= cartesian dimensional atomic
!!                           coordinates (bohr)
!!  acell(3)    = length scales of primitive translations (bohr)
!!  rprim(3,3) = dimensionless real space primitive translations
!!  rprimd(3,3) = dimensional real space primitive translations
!!                (bohr)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_pred_lotf,mover,pred_bfgs,pred_delocint,pred_diisrelax
!!      pred_isokinetic,pred_isothermal,pred_langevin,pred_moldyn,pred_nose
!!      pred_srkna14,pred_steepdesc,pred_verlet
!!
!! CHILDREN
!!
!! SOURCE

subroutine hist2var(acell,hist,natom,rprim,rprimd,xcart,xred,zDEBUG)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hist2var'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom
type(abihist),intent(in) :: hist
logical,intent(in) :: zDEBUG
!arrays
real(dp),intent(out) :: acell(3)
real(dp),intent(out) :: rprimd(3,3)
real(dp),intent(out) :: rprim(3,3)
real(dp),intent(out) :: xcart(3,natom)
real(dp),intent(out) :: xred(3,natom)

!Local variables-------------------------------
!scalars
integer :: jj,kk

! *************************************************************

 xcart(:,:) =hist%histXF(:,:,1,hist%ihist)
 xred(:,:)  =hist%histXF(:,:,2,hist%ihist)
 acell(:)   =hist%histA(:,hist%ihist)
 rprimd(:,:)=hist%histR(:,:,hist%ihist)

!Compute rprim from rprimd and acell
 do kk=1,3
   do jj=1,3
     rprim(jj,kk)=rprimd(jj,kk)/acell(kk)
   end do
 end do

 if(zDEBUG)then
   write (std_out,*) 'Atom positions and cell parameters '
   write (std_out,*) 'ihist: ',hist%ihist
   write (std_out,*) 'xcart:'
   do kk=1,natom
     write (std_out,*) xcart(:,kk)
   end do
   write (std_out,*) 'xred:'
   do kk=1,natom
     write (std_out,*) xred(:,kk)
   end do
   write(std_out,*) 'rprim:'
   do kk=1,3
     write(std_out,*) rprim(:,kk)
   end do
   write(std_out,*) 'rprimd:'
   do kk=1,3
     write(std_out,*) rprimd(:,kk)
   end do
   write(std_out,*) 'acell:'
   write(std_out,*) acell(:)
   call chkrprimd(acell,rprim,rprimd,std_out)
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
!! amass = mass of the atoms
!! natom = number of atoms
!! vel(3,natom)= Velocities of the atoms
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist<type abihist>=Historical record of positions, forces
!!                           acell, stresses, and energies,
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!
!! SOURCE

subroutine vel2hist(amass,hist,natom,vel)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vel2hist'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom
type(abihist),intent(inout) :: hist
!arrays
real(dp),intent(in) :: vel(3,natom)
real(dp),intent(in) :: amass(natom)

!Local variables-------------------------------
!scalars
integer :: ii,jj
real(dp) :: ekin

! *************************************************************

!Store the velocities
 hist%histV(:,:,hist%ihist)=vel(:,:)

!Compute the Ionic Kinetic energy
 ekin=0.0
 do ii=1,natom
   do jj=1,3
     ekin=ekin+0.5_dp*amass(ii)*vel(jj,ii)**2
   end do
 end do

!Store the Ionic Kinetic Energy
 hist%histEk(hist%ihist)=ekin

end subroutine vel2hist
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/abihist_compare
!! NAME
!! abihist_compare
!!
!! FUNCTION
!! Compare 2 HIST records
!!
!! INPUTS
!!  hist_in <type(abihist)>
!!  tolerance
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

subroutine abihist_compare(hist_in,hist_out,natom,similar,tolerance)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abihist_compare'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
integer,intent(in) :: natom
integer,intent(out) :: similar
real(dp),intent(in) :: tolerance
type(abihist),intent(in) :: hist_in
type(abihist),intent(inout) :: hist_out

!Local variables-------------------------------
!scalars
integer :: kk,jj
real(dp) :: maxdiff,diff
real(dp) :: x,y

! ***************************************************************

 similar=1

 write(std_out,*) 'Using values from history, iteration:',hist_in%ihist
 write(std_out,*) 'Differences between present history and values stored'
 write(std_out,*) 'on the previous history.(Relative difference)'

 x=hist_out%histXF(1,1,1,hist_out%ihist)
 y=hist_in%histXF(1,1,1,hist_in%ihist)
 maxdiff=2*abs(x-y)/(abs(x)+abs(y))
 do kk=1,natom
   do jj=1,3
     x=hist_out%histXF(jj,kk,1,hist_out%ihist)
     y=hist_in%histXF(jj,kk,1,hist_in%ihist)
     diff=2*abs(x-y)/(abs(x)+abs(y))
     if (diff>maxdiff) maxdiff=diff
   end do
 end do
 write(std_out,'(a,e12.5)') 'xcart:    ',maxdiff
 if (maxdiff>tolerance) similar=0


 x=hist_out%histXF(1,1,2,hist_out%ihist)
 y=hist_in%histXF(1,1,2,hist_in%ihist)
 maxdiff=2*abs(x-y)/(abs(x)+abs(y))
 do kk=1,natom
   do jj=1,3
     x=hist_out%histXF(jj,kk,2,hist_out%ihist)
     y=hist_in%histXF(jj,kk,2,hist_in%ihist)
     diff=2*abs(x-y)/(abs(x)+abs(y))
     if (diff>maxdiff) maxdiff=diff
   end do
 end do
 write(std_out,'(a,e12.5)') 'xred:     ',maxdiff
 if (maxdiff>tolerance) similar=0


 x=hist_out%histR(1,1,hist_out%ihist)
 y=hist_in%histR(1,1,hist_in%ihist)
 maxdiff=2*abs(x-y)/(abs(x)+abs(y))
 do kk=1,3
   do jj=1,3
     x=hist_out%histR(jj,kk,hist_out%ihist)
     y=hist_in%histR(jj,kk,hist_in%ihist)
     diff=2*abs(x-y)/(abs(x)+abs(y))
     if (diff>maxdiff) maxdiff=diff
   end do
 end do
 write(std_out,'(a,e12.5)') 'rprimd:   ',maxdiff
 if (maxdiff>tolerance) similar=0


 x=hist_out%histA(1,hist_out%ihist)
 y=hist_in%histA(1,hist_in%ihist)
 maxdiff=2*abs(x-y)/(abs(x)+abs(y))
 do kk=1,3
   x=hist_out%histA(kk,hist_out%ihist)
   y=hist_in%histA(kk,hist_in%ihist)
   diff=2*abs(x-y)/(abs(x)+abs(y))
   if (diff>maxdiff) maxdiff=diff
 end do
 write(std_out,'(a,e12.5)') 'acell:    ',maxdiff
 if (maxdiff>tolerance) similar=0

 if (similar==1) then

   hist_out%histV(:,:,hist_out%ihist)=   hist_in%histV(:,:,hist_in%ihist)
   hist_out%histE(hist_out%ihist)=       hist_in%histE(hist_in%ihist)
   hist_out%histEk(hist_out%ihist)=      hist_in%histEk(hist_in%ihist)
   hist_out%histEnt(hist_out%ihist)=     hist_in%histEnt(hist_in%ihist)
   hist_out%histT(hist_out%ihist)=       hist_in%histT(hist_in%ihist)
   hist_out%histXF(:,:,3,hist_out%ihist)=hist_in%histXF(:,:,3,hist_in%ihist)
   hist_out%histXF(:,:,4,hist_out%ihist)=hist_in%histXF(:,:,4,hist_in%ihist)
   hist_out%histS(:,hist_out%ihist)=     hist_in%histS(:,hist_in%ihist)

 end if

end subroutine abihist_compare
!!***

!!****f* m_abihist/write_md_hist
!!
!! NAME
!! write_md_hist
!!
!! FUNCTION
!! Write the history into a netcdf dataset
!!
!! INPUTS
!! filname = Filename of the file where the history will be stored
!! hist<type abihist>=Historical record of positions, forces
!!      |                    acell, stresses, and energies,
!!      |                    contains:
!!      | mxhist:  Maximun number of records
!!      | histA:   Historical record of acell(A) and rprimd(R)
!!      | histE:   Historical record of energy(E)
!!      | histEk:  Historical record of Ionic kinetic energy(Ek)
!!      | histEnt: Historical record of Entropy
!!      | histT:   Historical record of time(T) (For MD or iteration for GO)
!!      | histR:   Historical record of rprimd(R)
!!      | histS:   Historical record of strten(S)
!!      | histV:   Historical record of velocity(V)
!!      | histXF:  Historical record of positions(X) and forces(F)
!! mover<abimover>=Data used by the different predictors to update positions, acell, etc.
!!      | ntypat=Number of type of atoms.
!!      | npsp=Number of pseudopotentials.
!!      | typat(natom)=Type of each natom
!!      | znucl(npsp)=Nuclear charge for each type of pseudopotential.
!!      |             WARNING: alchemical mixing is not supported. We assume npsp == ntypat
!! itime=index of the present iteration 
!! icycle=index of the present cycle.
!!
!! TODO
!!  Have you ever heard about ETSF-IO specifications?
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!
!! SOURCE

subroutine write_md_hist(hist,mover,filename,icycle,itime)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'write_md_hist'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: filename 
 type(abimover),intent(in) :: mover
 integer,intent(in) :: itime,icycle
 type(abihist),intent(in) :: hist
!arrays

!Local variables-------------------------------
!scalars
 integer :: ncerr,ncid,npsp,itypat,iatom
 integer :: xyz_id,natom_id,time_id,six_id
 integer :: xcart_id,xred_id,fcart_id,fred_id
 integer :: vel_id,etotal_id,acell_id,rprimd_id,strten_id
 integer :: ntypat_id,npsp_id,typat_id,znucl_id
 integer :: ekin_id,entr_id,mdtime_id,amu_id,dtion_id !amass_id,
!arrays
 integer :: dimXFids(3),dimAids(2),dimEids(1),dimRids(3),dimSids(2)
 integer :: count1(1),start1(1)
 integer :: count2(2),start2(2)
 integer :: count3(3),start3(3)
 integer :: dimids_0d(0)
 real(dp) :: amu(mover%ntypat)

! *************************************************************************

 write(std_out,*) 'ihist @ write_md_hist',hist%ihist
 write(std_out,*) 'mxhist @ write_md_hist',hist%mxhist

#if defined HAVE_NETCDF

!#####################################################################
!### Creation of NetCDF file

 if (itime==1.and.icycle==1) then
!  1. Create netCDF file
   ncerr = nf90_create(path=trim(filename),cmode=NF90_CLOBBER, ncid=ncid)
   NCF_CHECK_MSG(ncerr,"create netcdf history file")

!  2. Define dimensions
   ncerr = nf90_def_dim(ncid,"natom",mover%natom,natom_id)
   NCF_CHECK_MSG(ncerr," define dimension natom")

   ncerr = nf90_def_dim(ncid,"ntypat",mover%ntypat,ntypat_id)
   NCF_CHECK_MSG(ncerr," define dimension ntypat")

   npsp = size(mover%znucl)
   if (npsp /= mover%ntypat) then
     MSG_WARNING("HIST file does not support alchemical mixing!")
   end if 
   ncerr = nf90_def_dim(ncid,"npsp",npsp,npsp_id)
   NCF_CHECK_MSG(ncerr," define dimension npsp")

   ncerr = nf90_def_dim(ncid,"xyz",3,xyz_id)
   NCF_CHECK_MSG(ncerr," define dimension xyz")

   ncerr = nf90_def_dim(ncid,"time",NF90_UNLIMITED,time_id)
   NCF_CHECK_MSG(ncerr," define dimension time")

   ncerr = nf90_def_dim(ncid,"six",6,six_id)
   NCF_CHECK_MSG(ncerr," define dimension six")

!  Dimensions for xcart,xred,fcart,fred and vel
   dimXFids = (/ xyz_id, natom_id, time_id /)
!  Dimensions for acell
   dimAids = (/ xyz_id, time_id /)
!  Dimensions for etotal
   dimEids = (/ time_id /)
!  Dimensions for rprimd
   dimRids = (/ xyz_id, xyz_id, time_id /)
!  Dimensions for strten
   dimSids = (/ six_id, time_id /)

!  3. Define variables and their attributes (units and mnemonics)

   ! ======================================
   ! Constant variables (written only once)
   ! ======================================
   ncerr = nf90_def_var(ncid, "typat", NF90_INT, natom_id, typat_id)
   NCF_CHECK_MSG(ncerr," define variable typat")

   ncerr = nf90_def_var(ncid, "znucl", NF90_DOUBLE, npsp_id, znucl_id)
   NCF_CHECK_MSG(ncerr," define variable znucl")

!   call ab_define_var(ncid, [natom_id], amass_id, NF90_DOUBLE,&
!&   "amass","Mass of each atom, in unit of electronic mass (=amu*1822...)","bohr" )

   call ab_define_var(ncid, [ntypat_id], amu_id, NF90_DOUBLE,&
&   "amu","Masses of each type of atom in atomic mass units", "" )

   ncerr = nf90_def_var(ncid, "dtion", NF90_DOUBLE, dimids_0d, dtion_id)
   NCF_CHECK_MSG(ncerr," define variable dtion")

   ! ==================
   ! Evolving variables
   ! ==================
   call ab_define_var(ncid, dimXFids, xcart_id, NF90_DOUBLE,&
&   "xcart","vectors (X) of atom positions in CARTesian coordinates","bohr" )

   call ab_define_var(ncid, dimXFids, xred_id, NF90_DOUBLE,&
&   "xred","vectors (X) of atom positions in REDuced coordinates","dimensionless" )

   call ab_define_var(ncid, dimXFids, fcart_id, NF90_DOUBLE,&
&   "fcart","atom Forces in CARTesian coordinates","Ha/bohr" )

   call ab_define_var(ncid, dimXFids, fred_id, NF90_DOUBLE,&
&   "fred","atom Forces in REDuced coordinates","dimensionless" )

   call ab_define_var(ncid, dimXFids, vel_id, NF90_DOUBLE,&
&   "vel","VELocity","bohr*Ha/hbar" )

   call ab_define_var(ncid, dimAids, acell_id, NF90_DOUBLE,&
&   "acell","CELL lattice vector scaling","bohr" )

   call ab_define_var(ncid, dimRids, rprimd_id, NF90_DOUBLE,&
&   "rprimd","Real space PRIMitive translations, Dimensional","bohr" )

   call ab_define_var(ncid, dimEids, etotal_id, NF90_DOUBLE,&
&   "etotal","TOTAL Energy","Ha" )

   call ab_define_var(ncid, dimEids, ekin_id, NF90_DOUBLE,&
&   "ekin","Energy KINetic ionic","Ha" )

   call ab_define_var(ncid, dimEids, entr_id, NF90_DOUBLE,&
&   "entropy","Entropy","" )

   call ab_define_var(ncid, dimEids, mdtime_id, NF90_DOUBLE,&
&   "mdtime","Molecular Dynamics TIME","hbar/Ha" )

   call ab_define_var(ncid, dimSids, strten_id, NF90_DOUBLE,&
&   "strten","STRess tensor","Ha/bohr^3" )

!  4. End define mode
   ncerr = nf90_enddef(ncid)
   NCF_CHECK_MSG(ncerr," end define mode")

   ! Write extra info on the cell needed by post-processing tools e.g. vesta, xcrysden,
   ! Note that these variables do not change and are not read in read_md_hist.
   ncerr = nf90_put_var(ncid, typat_id, mover%typat)
   NCF_CHECK_MSG(ncerr," write variable typat")

   ncerr = nf90_put_var(ncid, znucl_id, mover%znucl)
   NCF_CHECK_MSG(ncerr," write variable znucl")

   ! Get amu from amass.
   do itypat=1,mover%ntypat
     do iatom=1,mover%natom
       if (itypat == mover%typat(iatom)) then
         amu(itypat) = mover%amass(iatom) / amu_emass
         exit
       end if
     end do
   end do

   ncerr = nf90_put_var(ncid, amu_id, amu)
   NCF_CHECK_MSG(ncerr," write variable amu")

   !ncerr = nf90_put_var(ncid, amass_id, mover%amass)
   !NCF_CHECK_MSG(ncerr," write variable amass")

   ncerr = nf90_put_var(ncid, dtion_id, mover%dtion)
   NCF_CHECK_MSG(ncerr," write variable dtion")

 else

!  #####################################################################
!  ### Open the NetCDF file for write new iterations
   write(std_out,*) 'OPEN NETCDF FILE'

!  1. Open netCDF file
   ncerr = nf90_open(path=trim(filename),mode=NF90_WRITE, ncid=ncid)
   NCF_CHECK_MSG(ncerr," open netcdf history file")

!  2. Get the ID of a variables from their name
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

   ncerr = nf90_inq_varid(ncid, "acell", acell_id)
   NCF_CHECK_MSG(ncerr," get the id for acell")

   ncerr = nf90_inq_varid(ncid, "rprimd", rprimd_id)
   NCF_CHECK_MSG(ncerr," get the id for rprimd")

   ncerr = nf90_inq_varid(ncid, "etotal", etotal_id)
   NCF_CHECK_MSG(ncerr," get the id for etotal")

   ncerr = nf90_inq_varid(ncid, "ekin", ekin_id)
   NCF_CHECK_MSG(ncerr," get the id for ekin")

   ncerr = nf90_inq_varid(ncid, "entropy", entr_id)
   NCF_CHECK_MSG(ncerr," get the id for entropy")

   ncerr = nf90_inq_varid(ncid, "mdtime", mdtime_id)
   NCF_CHECK_MSG(ncerr," get the id for mdtime")

   ncerr = nf90_inq_varid(ncid, "strten", strten_id)
   NCF_CHECK_MSG(ncerr," get the id for strten")
 end if

!#####################################################################
!### Write variables into the dataset

 count3 = (/ 3, mover%natom, 1 /)
 start3 = (/ 1, 1, 1 /)
 start3(3)=hist%ihist

 ncerr = nf90_put_var(ncid, xcart_id,&
& hist%histXF(:,:,1,hist%ihist),&
& start = start3,count = count3)
 NCF_CHECK_MSG(ncerr," write variable xcart")

 ncerr = nf90_put_var(ncid, xred_id,&
& hist%histXF(:,:,2,hist%ihist),&
& start = start3,count = count3)
 NCF_CHECK_MSG(ncerr," write variable xred")

 ncerr = nf90_put_var(ncid, fcart_id,&
& hist%histXF(:,:,3,hist%ihist),&
& start = start3,count = count3)
 NCF_CHECK_MSG(ncerr," write variable fcart")

 ncerr = nf90_put_var(ncid, fred_id,&
& hist%histXF(:,:,4,hist%ihist),&
& start = start3,count = count3)
 NCF_CHECK_MSG(ncerr," write variable fred")

 ncerr = nf90_put_var(ncid, vel_id,&
& hist%histV(:,:,hist%ihist),&
& start = start3,count = count3)
 NCF_CHECK_MSG(ncerr," write variable vel")

!RPRIMD
 count3(2)=3
 ncerr = nf90_put_var(ncid, rprimd_id,&
& hist%histR(:,:,hist%ihist),&
& start = start3,count = count3)
 NCF_CHECK_MSG(ncerr," write variable rprimd")

!ACELL
 count2 = (/ 3, 1 /)
 start2 = (/ 1, 1 /)
 start2(2)=hist%ihist
 ncerr = nf90_put_var(ncid, acell_id,&
& hist%histA(:,hist%ihist),&
& start = start2,count = count2)
 NCF_CHECK_MSG(ncerr," write variable acell")

!STRTEN
 count2(1)=6
 ncerr = nf90_put_var(ncid, strten_id,&
& hist%histS(:,hist%ihist),&
& start = start2,count = count2)
 NCF_CHECK_MSG(ncerr," write variable strten")

!ETOTAL
 count1 = (/ 1 /)
 start1 = (/ 1 /)
 start1(1)=hist%ihist
 ncerr = nf90_put_var(ncid, etotal_id,&
& hist%histE(hist%ihist), start = (/ hist%ihist /) )
 NCF_CHECK_MSG(ncerr," write variable etotal")

!Ekin
 count1 = (/ 1 /)
 start1 = (/ 1 /)
 start1(1)=hist%ihist
 ncerr = nf90_put_var(ncid, ekin_id,&
& hist%histEk(hist%ihist), start = (/ hist%ihist /) )
 NCF_CHECK_MSG(ncerr," write variable ekin")

!Entropy
 count1 = (/ 1 /)
 start1 = (/ 1 /)
 start1(1)=hist%ihist
 ncerr = nf90_put_var(ncid, entr_id,&
& hist%histEnt(hist%ihist), start = (/ hist%ihist /) )
 NCF_CHECK_MSG(ncerr," write variable entropy")

!MDTIME
 count1 = (/ 1 /)
 start1 = (/ 1 /)
 start1(1)=hist%ihist
 ncerr = nf90_put_var(ncid, mdtime_id,&
& hist%histT(hist%ihist), start = (/ hist%ihist /) )
 NCF_CHECK_MSG(ncerr," write variable mdtime")

!#####################################################################
!### Close NetCDF file
 ncerr = nf90_close(ncid)
 NCF_CHECK_MSG(ncerr," close netcdf history file")
#endif

end subroutine write_md_hist
!!***

!----------------------------------------------------------------------

!!****f* m_abihist/read_md_hist
!!
!! NAME
!! read_md_hist
!!
!! FUNCTION
!! Read the history file from a netcdf file and store it into a hist dataset structure
!!
!! INPUTS
!!  filename = Filename of the NetCDF to read
!!
!! OUTPUT
!! hist<type abihist>=Historical record of positions, forces
!!      |                    acell, stresses, and energies,
!!      |                    contains:
!!      | mxhist:  Maximun number of records
!!      | ihist:   Index of present record of hist
!!      | histA:   Historical record of acell(A) and rprimd(R)
!!      | histE:   Historical record of energy(E)
!!      | histR:   Historical record of rprimd(R)
!!      | histS:   Historical record of strten(S)
!!      | histV:   Historical record of velocity(V)
!!      | histXF:  Historical record of positions(X) and forces(F)
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!
!! SOURCE

subroutine read_md_hist(filename,hist)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_md_hist'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
 type(abihist),intent(inout) :: hist
 character(len=*),intent(in) :: filename

!Local variables-------------------------------
!scalars
 integer :: ncerr,ncid
 integer :: natom,time
 integer :: xyz_id,natom_id,time_id,six_id
 integer :: xcart_id,xred_id,fcart_id,fred_id,ekin_id,entr_id,mdtime_id
 integer :: vel_id,etotal_id,acell_id,rprimd_id,strten_id
 character(len=5) :: char_tmp
 !character(len=500) :: msg
!arrays

! *************************************************************************

#if defined HAVE_NETCDF

!#####################################################################
!### Reading of NetCDF file

!1. Open netCDF file

 ncerr=nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=ncid)

 if(ncerr /= NF90_NOERR) then
   write(std_out,*) 'Could no open ',trim(filename),', starting from scratch'
   hist%ihist=0
   hist%mxhist=0
   return
 else
   write(std_out,*) 'Succesfully open ',trim(filename),' for reading'
   write(std_out,*) 'Extracting information from NetCDF file...'
   hist%ihist=0
   hist%mxhist=0
 end if

!2. Inquire dimensions IDs

 ncerr = nf90_inq_dimid(ncid,"natom",natom_id)
 NCF_CHECK_MSG(ncerr," inquire dimension ID for natom")

 ncerr = nf90_inq_dimid(ncid,"xyz",xyz_id)
 NCF_CHECK_MSG(ncerr," inquire dimension ID for xyz")

 ncerr = nf90_inq_dimid(ncid,"time",time_id)
 NCF_CHECK_MSG(ncerr," inquire dimension ID for time")

 ncerr = nf90_inq_dimid(ncid,"six",six_id)
 NCF_CHECK_MSG(ncerr," inquire dimension ID for six")

!3. Inquire dimensions lengths

 ncerr = nf90_inquire_dimension(ncid,natom_id,char_tmp,natom)
 NCF_CHECK_MSG(ncerr," inquire dimension ID for natom")

 ncerr = nf90_inquire_dimension(ncid,time_id,char_tmp,time)
 NCF_CHECK_MSG(ncerr," inquire dimension ID for time")

 write(std_out,*) 'Number of atoms readed:',natom
 write(std_out,*) 'Number of iterations recorded:',time

! write(msg, '(a,a,i4,a)' )ch10,&
!&  ' gstate : reading',time,' (x,f) history pairs from input wf file.'
! call wrtout(std_out,msg,'COLL')
! call wrtout(ab_out,msg,'COLL')

!4. Allocate hist structure 

 hist%ihist=1
 hist%mxhist=time
 ABI_ALLOCATE(hist%histA,(3,hist%mxhist))
 ABI_ALLOCATE(hist%histE,(hist%mxhist))
 ABI_ALLOCATE(hist%histEk,(hist%mxhist))
 ABI_ALLOCATE(hist%histEnt,(hist%mxhist))
 ABI_ALLOCATE(hist%histT,(hist%mxhist))
 ABI_ALLOCATE(hist%histR,(3,3,hist%mxhist))
 ABI_ALLOCATE(hist%histS,(6,hist%mxhist))
 ABI_ALLOCATE(hist%histV,(3,natom,hist%mxhist))
 ABI_ALLOCATE(hist%histXF,(3,natom,4,hist%mxhist))

!5. Get the ID of a variables from their name

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

 ncerr = nf90_inq_varid(ncid, "acell", acell_id)
 NCF_CHECK_MSG(ncerr," get the id for acell")

 ncerr = nf90_inq_varid(ncid, "rprimd", rprimd_id)
 NCF_CHECK_MSG(ncerr," get the id for rprimd")

 ncerr = nf90_inq_varid(ncid, "etotal", etotal_id)
 NCF_CHECK_MSG(ncerr," get the id for etotal")

 ncerr = nf90_inq_varid(ncid, "ekin", ekin_id)
 NCF_CHECK_MSG(ncerr," get the id for ekin")

 ncerr = nf90_inq_varid(ncid, "entropy", entr_id)
 NCF_CHECK_MSG(ncerr," get the id for entropy")

 ncerr = nf90_inq_varid(ncid, "mdtime", mdtime_id)
 NCF_CHECK_MSG(ncerr," get the id for mdtime")

 ncerr = nf90_inq_varid(ncid, "strten", strten_id)
 NCF_CHECK_MSG(ncerr," get the id for strten")

!#####################################################################
!### Read variables from the dataset and write them into hist

!XCART,XRED,FCART,FRED,VEL
 ncerr = nf90_get_var(ncid, xcart_id,hist%histXF(:,:,1,:))
 NCF_CHECK_MSG(ncerr," read variable xcart")

 ncerr = nf90_get_var(ncid, xred_id,hist%histXF(:,:,2,:))
 NCF_CHECK_MSG(ncerr," read variable xred")

 ncerr = nf90_get_var(ncid, fcart_id,hist%histXF(:,:,3,:))
 NCF_CHECK_MSG(ncerr," read variable fcart")

 ncerr = nf90_get_var(ncid, fred_id,hist%histXF(:,:,4,:))
 NCF_CHECK_MSG(ncerr," read variable fred")

 ncerr = nf90_get_var(ncid, vel_id,hist%histV(:,:,:))
 NCF_CHECK_MSG(ncerr," read variable vel")

!RPRIMD
 ncerr = nf90_get_var(ncid, rprimd_id,hist%histR(:,:,:))
 NCF_CHECK_MSG(ncerr," read variable rprimd")

!ACELL
 ncerr = nf90_get_var(ncid, acell_id,hist%histA(:,:))
 NCF_CHECK_MSG(ncerr," read variable acell")

!STRTEN
 ncerr = nf90_get_var(ncid, strten_id,hist%histS(:,:))
 NCF_CHECK_MSG(ncerr," read variable strten")

!ETOTAL
 ncerr = nf90_get_var(ncid, etotal_id,hist%histE(:))
 NCF_CHECK_MSG(ncerr," read variable etotal")

!Ekin
 ncerr = nf90_get_var(ncid, ekin_id,hist%histEk(:))
 NCF_CHECK_MSG(ncerr," read variable ekin")

!Entropy
 ncerr = nf90_get_var(ncid, entr_id,hist%histEnt(:))
 NCF_CHECK_MSG(ncerr," read variable entropy")

!MDTime
 ncerr = nf90_get_var(ncid, mdtime_id,hist%histT(:))
 NCF_CHECK_MSG(ncerr," read variable mdtime")

!#####################################################################
!### Close NetCDF file
 ncerr = nf90_close(ncid)
 NCF_CHECK_MSG(ncerr," close netcdf history file")
#endif

end subroutine read_md_hist
!!***

!----------------------------------------------------------------------

end module m_abihist
