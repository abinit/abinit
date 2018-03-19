!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtxfase
!!
!! NAME
!! prtxfase
!!
!! FUNCTION
!! Print the values of xcart (X), forces (F)
!! acell (A), Stresses (S), and energy (E)
!! All values come from the history hist
!! Also compute and print max and rms forces.
!! Also compute absolute and relative differences
!! with previous calculation
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! ab_mover<type abimover>=Subset of dtset only related with
!!          |                 movement of ions and acell, contains:
!!          | dtion:  Time step
!!          ! natom:  Number of atoms
!!          | vis:    viscosity
!!          | iatfix: Index of atoms and directions fixed
!!          | amass:  Mass of ions
!! hist<type abihist>=Historical record of positions, forces,
!!                               stresses, cell and energies,
!! itime= time step
!! iout=unit number for printing
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      gettag,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine prtxfase(ab_mover,hist,itime,iout,pos)

 use defs_basis
 use m_profiling_abi
 use m_abimover
 use m_abihist

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtxfase'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
 type(abimover),intent(in) :: ab_mover
 type(abihist),intent(in),target :: hist
 integer,intent(in) :: itime,iout
 integer,intent(in) :: pos
!arrays

!Local variables-------------------------------
!scalars
 integer :: jj,kk,unfixd,iprt
 real(dp) :: val_max,val_rms,ucvol ! Values maximal and RMS, Volume of Unitary cell
 real(dp) :: dEabs,dErel ! Diff of energy absolute and relative
 real(dp) :: ekin
 real(dp) :: angle(3),rmet(3,3)
!character(len=80*(max(ab_mover%natom,3)+1)) :: message
!MGNAG: This is not very safe. One should use line-based output istead of appending chars
! and then outputting everything! For the time being I use this temporary hack to solve the problem with NAG
 character(len=max(80*(max(ab_mover%natom,3)+1),50000)) :: message
 character(len=18)   :: fmt1
 logical :: prtallatoms
!arrays
 logical :: atlist(ab_mover%natom)
 real(dp),allocatable :: fred(:,:),xcart(:,:)
 real(dp),pointer :: acell(:),fcart(:,:),rprimd(:,:),strten(:),vel(:,:),xred(:,:)

! ***********************************************************

 fmt1='(a,a,1p,3e22.14)'

!##########################################################
!### 1. Organize list of atoms to print

 prtallatoms=.TRUE.
 do kk=1,ab_mover%natom
   if (ab_mover%prtatlist(kk)/=kk) prtallatoms=.FALSE.
 end do

 atlist(:)=.FALSE.
 do iprt=1,ab_mover%natom
   if (ab_mover%prtatlist(iprt)>0.and.ab_mover%prtatlist(iprt)<=ab_mover%natom) atlist(ab_mover%prtatlist(iprt))=.TRUE.
 end do

 acell  => hist%acell(:,hist%ihist)
 rprimd => hist%rprimd(:,:,hist%ihist)
 xred   => hist%xred(:,:,hist%ihist)
 fcart  => hist%fcart(:,:,hist%ihist)
 strten => hist%strten(:,hist%ihist)
 vel    => hist%vel(:,:,hist%ihist)

!###########################################################
!### 1. Positions

 ABI_ALLOCATE(xcart,(3,ab_mover%natom))
 call xred2xcart(ab_mover%natom,rprimd,xcart,xred)

 write(message, '(a,a)' )&
& ch10,' Cartesian coordinates (xcart) [bohr]'
 call prtnatom(atlist,iout,message,ab_mover%natom,prtallatoms,xcart)

 write(message, '(a)' )&
& ' Reduced coordinates (xred)'
 call prtnatom(atlist,iout,message,ab_mover%natom,prtallatoms,xred)

 ABI_DEALLOCATE(xcart)

!###########################################################
!### 2. Forces

 if(pos==mover_AFTER)then

   ABI_ALLOCATE(fred,(3,ab_mover%natom))
   call fcart2fred(fcart,fred,rprimd,ab_mover%natom)

!  Compute max |f| and rms f,
!  EXCLUDING the components determined by iatfix
   val_max=0.0_dp
   val_rms=0.0_dp
   unfixd=0
   do kk=1,ab_mover%natom
     do jj=1,3
       if (ab_mover%iatfix(jj,kk) /= 1) then
         unfixd=unfixd+1
         val_rms=val_rms+fcart(jj,kk)**2
         val_max=max(val_max,abs(fcart(jj,kk)**2))
       end if
     end do
   end do
   if ( unfixd /= 0 ) val_rms=sqrt(val_rms/dble(unfixd))

   write(message, '(a,1p,2e12.5,a)' ) &
&   ' Cartesian forces (fcart) [Ha/bohr]; max,rms=',&
&   sqrt(val_max),val_rms,' (free atoms)'
   call prtnatom(atlist,iout,message,ab_mover%natom,prtallatoms,fcart)

   write(message, '(a)' )&
&   ' Reduced forces (fred)'
   call prtnatom(atlist,iout,message,ab_mover%natom,prtallatoms,fred)

   ABI_DEALLOCATE(fred)

 end if

!###########################################################
!### 3. Velocities

!Only if the velocities are being used
 if (hist%isVused)then
!  Only if velocities are recorded in a history
   if (allocated(hist%vel))then
!    Compute max |v| and rms v,
!    EXCLUDING the components determined by iatfix
     val_max=0.0_dp
     val_rms=0.0_dp
     unfixd=0
     do kk=1,ab_mover%natom
       do jj=1,3
         if (ab_mover%iatfix(jj,kk) /= 1) then
           unfixd=unfixd+1
           val_rms=val_rms+vel(jj,kk)**2
           val_max=max(val_max,abs(vel(jj,kk)**2))
         end if
       end do
     end do
     if ( unfixd /= 0 ) val_rms=sqrt(val_rms/dble(unfixd))

     write(message, '(a,1p,2e12.5,a)' ) &
&     ' Cartesian velocities (vel) [bohr*Ha/hbar]; max,rms=',&
&     sqrt(val_max),val_rms,' (free atoms)'
     call prtnatom(atlist,iout,message,ab_mover%natom,prtallatoms,vel)

!    Compute the ionic kinetic energy (no cell shape kinetic energy yet)
     ekin=0.0_dp
     do kk=1,ab_mover%natom
       do jj=1,3
!        Warning : the fixing of atoms is implemented in reduced
!        coordinates, so that this expression is wrong
         if (ab_mover%iatfix(jj,kk) == 0) then
           ekin=ekin+0.5_dp*ab_mover%amass(kk)*vel(jj,kk)**2
         end if
       end do
     end do
     write(message, '(a,1p,e22.14,a)' )&
&     ' Kinetic energy of ions (ekin) [Ha]=',ekin
     call wrtout(iout,message,'COLL')

   end if
 end if

!###########################################################
!### 3. ACELL

!Only if the acell is being used
 if (hist%isARused)then
!  Only if acell is recorded in a history
   if (allocated(hist%acell))then

     write(message, '(a)' ) &
&     ' Scale of Primitive Cell (acell) [bohr]'
     write(message,fmt1)&
&     TRIM(message),ch10,acell(:)
     call wrtout(iout,message,'COLL')
   end if
 end if

!###########################################################
!### 4. RPRIMD

!Only if the acell is being used
 if (hist%isARused)then
!  Only if rprimd is recorded in a history
   if (allocated(hist%rprimd))then
     write(message, '(a)' ) &
&     ' Real space primitive translations (rprimd) [bohr]'
     do kk=1,3
       write(message,fmt1)&
&       TRIM(message),ch10,&
&       rprimd(:,kk)
     end do
     call wrtout(iout,message,'COLL')
   end if
 end if

!###########################################################
!### 5. Unitary cell volume

 if (ab_mover%optcell/=0)then

   ucvol=&
&   rprimd(1,1)*(rprimd(2,2)*rprimd(3,3)-rprimd(3,2)*rprimd(2,3))+&
&   rprimd(2,1)*(rprimd(3,2)*rprimd(1,3)-rprimd(1,2)*rprimd(3,3))+&
&   rprimd(3,1)*(rprimd(1,2)*rprimd(2,3)-rprimd(2,2)*rprimd(1,3))

   write(message, '(a,1p,e22.14)' )&
&   ' Unitary Cell Volume (ucvol) [Bohr^3]=',&
&   ucvol
   call wrtout(iout,message,'COLL')

!  ###########################################################
!  ### 5. Angles and lengths

!  Compute real space metric.
   rmet = MATMUL(TRANSPOSE(rprimd),rprimd)

   angle(1)=acos(rmet(2,3)/sqrt(rmet(2,2)*rmet(3,3)))/two_pi*360.0d0
   angle(2)=acos(rmet(1,3)/sqrt(rmet(1,1)*rmet(3,3)))/two_pi*360.0d0
   angle(3)=acos(rmet(1,2)/sqrt(rmet(1,1)*rmet(2,2)))/two_pi*360.0d0

   write(message, '(a,a)' ) &
&   ' Angles (23,13,12)= [degrees]'
   write(message,fmt1)&
&   TRIM(message),ch10,&
&   angle(:)
   call wrtout(iout,message,'COLL')

   write(message, '(a,a)' ) &
&   ' Lengths [Bohr]'
   write(message,fmt1)&
&   TRIM(message),ch10,&
&   sqrt(rmet(1,1)),sqrt(rmet(2,2)),sqrt(rmet(3,3))
   call wrtout(iout,message,'COLL')


!  ###########################################################
!  ### 5. Stress Tensor

   if(pos==mover_AFTER)then
!    Only if strten is recorded in a history
     if (allocated(hist%strten))then

       write(message, '(a)' ) &
&       ' Stress tensor in cartesian coordinates (strten) [Ha/bohr^3]'

       write(message,fmt1)&
&       TRIM(message),ch10,&
&       strten(1),strten(6),strten(5)
       write(message,fmt1)&
&       TRIM(message),ch10,&
&       strten(6),strten(2),strten(4)
       write(message,fmt1)&
&       TRIM(message),ch10,&
&       strten(5),strten(4),strten(3)
       call wrtout(iout,message,'COLL')
     end if
   end if
 end if

!###########################################################
!### 6. Energy

 if(pos==mover_AFTER)then
   write(message, '(a,1p,e22.14)' )&
&   ' Total energy (etotal) [Ha]=',&
&   hist%etot(hist%ihist)

   if (itime>1)then
     jj = abihist_findIndex(hist,-1)
     dEabs=hist%etot(hist%ihist)-hist%etot(jj)
     dErel=2*dEabs/(abs(hist%etot(hist%ihist))+&
&     abs(hist%etot(jj)))
     write(message, '(a,a,a,a)' )&
&     TRIM(message),ch10,ch10,&
&     ' Difference of energy with previous step (new-old):'
     write(message, '(a,a,10a,a,1p,e12.5,a,10a,a,1p,e12.5)')&
&     TRIM(message),ch10,&
&     (' ',jj=1,10),' Absolute (Ha)=',dEabs,ch10,&
&     (' ',jj=1,10),' Relative     =',dErel
   end if
   call wrtout(iout,message,'COLL')
 end if

 contains
!!***

!!****f* ABINIT/gettag
!!
!! NAME
!! gettag
!!
!! FUNCTION
!! Set the tag associated to each atom,
!!
!! INPUTS
!! prtallatoms = Logical for PRTint ALL ATOMS
!! atlist      = ATom LIST
!! index       = index for each atom
!! natom       = Number of ATOMs
!!
!! OUTPUT
!!  tag = The string to put for each atom
!!
!! PARENTS
!!      prtxfase
!!
!! CHILDREN
!!      gettag,wrtout
!!
!! SOURCE

subroutine gettag(atlist,index,natom,prtallatoms,tag)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gettag'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
  logical,intent(in) :: prtallatoms
  integer,intent(in) :: natom
  logical,intent(in) :: atlist(natom)
  integer,intent(in) :: index
  character(len=7),intent(out)   :: tag

!Local variables -------------------------

! *********************************************************************
!The numbering will be from (1) to (9999)

   if (prtallatoms)then
     tag=''
   elseif (atlist(index)) then
     if (natom<10) then
       write(tag, '(a,I1.1,a)') ' (',index,')'
     elseif (natom<100) then
       write(tag, '(a,I2.2,a)') ' (',index,')'
     elseif (natom<1000) then
       write(tag, '(a,I3.3,a)') ' (',index,')'
     elseif (natom<10000) then
       write(tag, '(a,I4.4,a)') ' (',index,')'
     end if
   end if

 end subroutine gettag
!!***

!!****f* ABINIT/prtnatom
!!
!! NAME
!! prtnatom
!!
!! FUNCTION
!! Print information for N atoms
!!
!!
!! INPUTS
!! prtallatoms = Logical for PRTint ALL ATOMS
!! atlist      = ATom LIST
!! index       = index for each atom
!! natom       = Number of ATOMs
!!
!! OUTPUT
!!  tag = The string to put for aech atom
!!
!! PARENTS
!!      prtxfase
!!
!! CHILDREN
!!      gettag,wrtout
!!
!! SOURCE


subroutine prtnatom(atlist,iout,message,natom,prtallatoms,thearray)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtnatom'
 use interfaces_14_hidewrite
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
  logical,intent(in) :: prtallatoms
  integer,intent(in) :: natom
  logical,intent(in) :: atlist(natom)
  integer,intent(in) :: iout
  character(len=*),intent(inout) :: message
!arrays
  real(dp) :: thearray(3,natom)

!Local variables-------------------------------
!scalars
  integer :: kk
  character(len=7)   :: tag ! Maximal ' (9999)'
  character(len=18)   :: fmt

! *********************************************************************

   fmt='(a,a,1p,3e22.14,a)'

   do kk=1,natom

     if (atlist(kk)) then
       call gettag(atlist,kk,natom,prtallatoms,tag)
       write(message,fmt)&
&       TRIM(message),ch10,&
&       thearray(:,kk),&
&       tag
     end if

   end do
 !MGNAG
 ! Runtime Error: wrtout_cpp.f90, line 896: Buffer overflow on output
   call wrtout(iout,message,'COLL')

 end subroutine prtnatom
!!***

end subroutine prtxfase
!!***
