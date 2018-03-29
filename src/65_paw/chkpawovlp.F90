!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkpawovlp
!! NAME
!! chkpawovlp
!!
!! FUNCTION
!! Verify that the paw spheres are not overlapping
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  natom=number of atoms in cell.
!!  ntypat=number of types of atoms in unit cell.
!!  pawovlp=percentage of voluminal overlap ratio allowed to continue execution
!!          (if negative value, execution always continues)
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data:
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).
!!  typat(natom)=type (integer) for each atom
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (only checking)
!!
!! NOTES
!!
!! PARENTS
!!      bethe_salpeter,respfn,scfcv,screening,sigma,wfk_analyze
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine chkpawovlp(natom,ntypat,pawovlp,pawtab,rmet,typat,xred)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_pawtab, only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkpawovlp'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 real(dp) :: pawovlp
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: rmet(3,3),xred(3,natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: ia,ib,ii,t1,t2,t3
 logical :: stop_on_error
 real(dp) :: dd,dif1,dif2,dif3,ha,hb,norm2
 real(dp) :: ratio_percent,va,vb,vv
 character(len=750) :: message
!arrays
 integer :: iamax(2),ibmax(2),iovl(2)
 real(dp) :: norm2_min(2),r2cut(2),ratio_percent_max(2),rcuta(2),rcutb(2)


! *************************************************************************

 DBG_ENTER("COLL")

 iamax(:)=-1;ibmax(:)=-1
 norm2_min(:)=-1.d0;ratio_percent_max(:)=-1.d0
 iovl(:)=0

!Loop on "overlapping" atoms with the maximum overlap
 do ia=1,natom

   rcuta(1)=pawtab(typat(ia))%rpaw
   rcuta(2)=pawtab(typat(ia))%rshp

   do ib=ia,natom

     rcutb(1)=pawtab(typat(ib))%rpaw
     rcutb(2)=pawtab(typat(ib))%rshp
     r2cut(1)=(rcuta(1)+rcutb(1))**2
     r2cut(2)=(rcuta(2)+rcutb(2))**2

!    Visit the box and its first images:
     do t3=-1,1
       do t2=-1,1
         do t1=-1,1

           dif1=xred(1,ia)-(xred(1,ib)+dble(t1))
           dif2=xred(2,ia)-(xred(2,ib)+dble(t2))
           dif3=xred(3,ia)-(xred(3,ib)+dble(t3))
           norm2=sqnrm_pawovlp(dif1,dif2,dif3)

           do ii=1,2

             if(norm2>tol10.and.norm2<r2cut(ii)) then

               iovl(ii)=iovl(ii)+1

!              Compute the overlap ratio:
               dd=sqrt(norm2)
               va=4._dp/3._dp*pi*rcuta(ii)**3
               vb=4._dp/3._dp*pi*rcutb(ii)**3
               ha=(rcutb(ii)**2-(dd-rcuta(ii))**2)/(two*dd)
               hb=(rcuta(ii)**2-(dd-rcutb(ii))**2)/(two*dd)
               vv=pi/3.d0*(ha**2*(three*rcuta(ii)-ha)+hb**2*(three*rcutb(ii)-hb))
               ratio_percent=100._dp*min(vv/min(va,vb),one)
               if (ratio_percent>ratio_percent_max(ii)) then
                 ratio_percent_max(ii)=ratio_percent
                 norm2_min(ii)=norm2
                 iamax(ii)=ia;ibmax(ii)=ib
               end if

             end if
           end do
         end do
       end do
     end do
   end do
 end do

 stop_on_error=(abs(pawovlp)<=tol6.or.(pawovlp>tol6.and.ratio_percent_max(1)>pawovlp))

!Print adapted message with overlap value
 if (iovl(1)+iovl(2)>0) then

   !ii=1: PAW augmentation regions overlap
   !ii=2: compensation charges overlap
   if (iovl(2)==0) ii=1
   if (iovl(2)> 0) ii=2

   if (iovl(ii)>0) then

     if (ii==1) write(message,' (a)' ) 'PAW SPHERES ARE OVERLAPPING!'
     if (ii==2) write(message, '(2a)' )'PAW COMPENSATION DENSITIES ARE OVERLAPPING !!!!'

     if (iovl(ii)==1) then
       write(message, '(3a)' ) trim(message),ch10,&
&       '   There is one pair of overlapping atoms.'
     else
       write(message, '(3a,i5,a)' ) trim(message),ch10,&
&       '   There are ', iovl(1),' pairs of overlapping atoms.'
     end if
     write(message, '(3a,i3,a,i3,a)' ) trim(message),ch10,&
     '   The maximum overlap percentage is obtained for the atoms ',iamax(ii),' and ',ibmax(ii),'.'
     write(message, '(2a,2(a,i3),a,f9.5,a,2(a,i3,a,f9.5,a),a,f5.2,a)' ) trim(message),ch10,&
&     '    | Distance between atoms ',iamax(ii),' and ',ibmax(ii),' is  : ',sqrt(norm2_min(ii)),ch10,&
&     '    | PAW radius of the sphere around atom ',iamax(ii),' is: ',pawtab(typat(iamax(ii)))%rpaw,ch10,&
&     '    | PAW radius of the sphere around atom ',ibmax(ii),' is: ',pawtab(typat(ibmax(ii)))%rpaw,ch10,&
&     '    | This leads to a (voluminal) overlap ratio of ',ratio_percent_max(ii),' %'
     if (ii==2) then
       write(message, '(3a)' ) trim(message),ch10,&
&       'THIS IS DANGEROUS !, as PAW formalism assumes non-overlapping compensation densities.'
     end if

     if (stop_on_error) then
       MSG_ERROR_NOSTOP(message,ia) !ia is dummy
     else
       MSG_WARNING(message)
     end if

   end if

!  Print advice
   if (stop_on_error) then
     write(message, '(3a)' )&
&     '  Action: 1- decrease cutoff radius of PAW dataset',ch10,&
&     '    OR  2- ajust "pawovlp" input variable to allow overlap (risky)'
     MSG_ERROR(message)
   end if

!  Print last message if execution continues:
   if (pawovlp<=tol6) then
     write(message, '(6a)' ) &
&     '       Results might be approximate,',ch10,&
&     '       and even inaccurate (if overlap is too big) !',ch10,&
&     '       Assume experienced user. Execution will continue.',ch10
     call wrtout(std_out,message,'COLL')
   else if (ratio_percent_max(1)<=pawovlp) then
     write(message, '(8a)' ) &
&     '       Overlap ratio seems to be acceptable (less than value',ch10,&
&     '       of "pawovlp" input parameter): execution will continue.',ch10,&
&     '       But be aware that results might be approximate,',ch10,&
&     '       and even inaccurate (depending on your physical system) !',ch10
     call wrtout(std_out,message,'COLL')
   end if

 end if !iovl>0

 DBG_EXIT("COLL")

 contains

   function sqnrm_pawovlp(u1,u2,u3)
!squared norm of a vector

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sqnrm_pawovlp'
!End of the abilint section

   real(dp) :: sqnrm_pawovlp
   real(dp),intent(in) :: u1,u2,u3

   sqnrm_pawovlp=rmet(1,1)*u1*u1+rmet(2,1)*u2*u1+rmet(3,1)*u3*u1&
&   +rmet(1,2)*u1*u2+rmet(2,2)*u2*u2+rmet(3,2)*u3*u2&
&   +rmet(1,3)*u1*u3+rmet(2,3)*u2*u3+rmet(3,3)*u3*u3

 end function sqnrm_pawovlp
end subroutine chkpawovlp
!!***
