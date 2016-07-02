!{\src2tex{textfont=tt}}
!!****f* ABINIT/fconv
!!
!! NAME
!! fconv
!!
!! FUNCTION
!! Check maximal absolute value of force (hartree/bohr) against
!! input tolerance; if below tolerance, return iexit=1.
!! Takes into account the fact that the Broyden (or moldyn) step
!! might be the last one (last itime), to print eventually modified message.
!! Stresses are also included in the check, provided that optcell/=0.
!! If optcell=1, takes only the trace into account
!!    optcell=2, takes all components into account
!!    optcell=3, takes traceless stress into account
!!    optcell=4, takes sigma(1 1) into account
!!    optcell=5, takes sigma(2 2) into account
!!    optcell=6, takes sigma(3 3) into account
!!    optcell=7, takes sigma(2,2),(2,3) and (3 3) into account
!!    optcell=8, takes sigma(1,1),(1,3) and (3 3) into account
!!    optcell=9, takes sigma(1,1),(1,2) and (2 2) into account
!! In the case of stresses, target the tensor strtarget, and
!! take into account the factor strfact
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  fcart(3,natom)=cartesian forces on atoms in hartree/bohr
!!  iatfix(3,natom)=1 for frozen atom, 0 for unfrozen
!!  itime=current number of Broyden/Moldyn iterations
!!  natom=number of atoms in unit cell
!!  ntime=maximum number of Broyden/Moldyn iterations allowed
!!  optcell=option for taking stresses into account (see above)
!!  strfact=factor that multiplies the stresses when they are compared to forces.
!!  strtarget(6)=components of the target stress tensor (hartree/bohr^3)
!!  strten(6)=components of the stress tensor (hartree/bohr^3)
!!  tolmxf=tolerance on maximal absolute value of components of forces
!!
!! OUTPUT
!!  writes to unit std_out and to ab_out, and returns
!!
!! SIDE EFFECTS
!! Input/Output
!!  at input  : iexit=  0 if not the last itime,  1 if the last itime
!!  at output : iexit=  0 if not below tolerance, 1 if below tolerance
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine fconv(fcart,iatfix,iexit,itime,natom,ntime,optcell,strfact,strtarget,strten,tolmxf)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_fstrings,  only : indent

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fconv'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime,natom,ntime,optcell
 integer,intent(inout) :: iexit
 real(dp),intent(in) :: strfact,tolmxf
!arrays
 integer,intent(in) :: iatfix(3,natom)
 real(dp),intent(in) :: fcart(3,natom),strtarget(6),strten(6)

!Local variables-------------------------------
!scalars
 integer :: iatom,idir,istr
 real(dp) :: fmax,strdiag
 character(len=500) :: message
!arrays
 real(dp) :: dstr(6)

! *************************************************************************

!Compute maximal component of forces, EXCLUDING any fixed components
 fmax=0.0_dp
 do iatom=1,natom
   do idir=1,3
     if (iatfix(idir,iatom) /= 1) then
       if( abs(fcart(idir,iatom)) >= fmax ) fmax=abs(fcart(idir,iatom))
     end if
   end do
 end do

 dstr(:)=strten(:)-strtarget(:)

!Eventually take into account the stress
 if(optcell==1)then
   strdiag=(dstr(1)+dstr(2)+dstr(3))/3.0_dp
   if(abs(strdiag)*strfact >= fmax ) fmax=abs(strdiag)*strfact
 else if(optcell==2)then
   do istr=1,6
     if(abs(dstr(istr))*strfact >= fmax ) fmax=abs(dstr(istr))*strfact
   end do
 else if(optcell==3)then
!  Must take away the trace from diagonal elements
   strdiag=(dstr(1)+dstr(2)+dstr(3))/3.0_dp
   do istr=1,3
     if(abs(dstr(istr)-strdiag)*strfact >= fmax ) fmax=abs(dstr(istr)-strdiag)*strfact
   end do
   do istr=4,6
     if(abs(dstr(istr))*strfact >= fmax ) fmax=abs(dstr(istr))*strfact
   end do
 else if(optcell==4 .or. optcell==5 .or. optcell==6)then
   if(abs(dstr(optcell-3))*strfact >= fmax ) fmax=abs(dstr(optcell-3))*strfact
 else if(optcell==7)then
   if(abs(dstr(2))*strfact >= fmax ) fmax=abs(dstr(2))*strfact
   if(abs(dstr(3))*strfact >= fmax ) fmax=abs(dstr(3))*strfact
   if(abs(dstr(4))*strfact >= fmax ) fmax=abs(dstr(4))*strfact
 else if(optcell==8)then
   if(abs(dstr(1))*strfact >= fmax ) fmax=abs(dstr(1))*strfact
   if(abs(dstr(3))*strfact >= fmax ) fmax=abs(dstr(3))*strfact
   if(abs(dstr(5))*strfact >= fmax ) fmax=abs(dstr(5))*strfact
 else if(optcell==9)then
   if(abs(dstr(1))*strfact >= fmax ) fmax=abs(dstr(1))*strfact
   if(abs(dstr(2))*strfact >= fmax ) fmax=abs(dstr(2))*strfact
   if(abs(dstr(6))*strfact >= fmax ) fmax=abs(dstr(6))*strfact
 end if

 if (fmax<tolmxf) then
   write(message, '(a,a,i4,a,a,a,es11.4,a,es11.4,a,a)' ) ch10,&
&   ' At Broyd/MD step',itime,', gradients are converged : ',ch10,&
&   '  max grad (force/stress) =',fmax,' < tolmxf=',tolmxf,' ha/bohr (free atoms)',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
   iexit=1
 else
   if(iexit==1)then
     write(message, '(a,a,a,a,i5,a,a,a,es11.4,a,es11.4,a,a)' ) ch10,&
&     ' fconv : WARNING -',ch10,&
&     '  ntime=',ntime,' was not enough Broyd/MD steps to converge gradients: ',ch10,&
&     '  max grad (force/stress) =',fmax,' > tolmxf=',tolmxf,' ha/bohr (free atoms)',ch10
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')

     write(std_out,"(8a)")ch10,&
&     "--- !RelaxConvergenceWarning",ch10,&
&     "message: | ",ch10,TRIM(indent(message)),ch10,&
&     "..."

   else
     write(message, '(a,i4,a,a,a,es11.4,a,es11.4,a,a)' ) &
&     ' fconv : at Broyd/MD step',itime,', gradients have not converged yet. ',ch10,&
&     '  max grad (force/stress) =',fmax,' > tolmxf=',tolmxf,' ha/bohr (free atoms)',ch10
     call wrtout(std_out,message,'COLL')
   end if
   iexit=0
 end if

end subroutine fconv
!!***
