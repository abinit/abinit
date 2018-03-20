!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawxenergy
!! NAME
!! pawxenergy
!!
!! FUNCTION
!! Compute contributions to energy for PAW+ local exact exchange calculations
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (BA,FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors.
!!
!! INPUTS
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrhoij <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab <type(pawtab_type)>=paw tabulated starting data:
!!     %lexexch=l used for local exact-exchange
!!     %vex(2*lexexch+1*4)=screened coulomb matrix
!!
!! SIDE EFFECTS
!!  eexex=energy is updated with the contribution of the cuyrrent atom
!!
!! PARENTS
!!      pawdenpot
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine pawxenergy(eexex,pawprtvol,pawrhoij,pawtab)

 use m_profiling_abi
 use defs_basis
 use m_errors

 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawxenergy'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: pawprtvol
 real(dp),intent(inout) :: eexex
 type(pawrhoij_type),intent(in) :: pawrhoij
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
 integer :: irhoij,irhoij1,ispden,jrhoij,jrhoij1,klmn,klmn1,lexexch,ll,m11,m21,m31,m41,n1
 integer :: n2,n3,n4,nk,nn1,nn2
 real(dp) :: eexextemp
 character(len=500) :: message
!arrays
 integer :: indn(3,3)
 real(dp) :: factnk(6)

! *****************************************************

 DBG_ENTER("COLL")

 lexexch=pawtab%lexexch
 if (pawtab%nproju==1) nk=1
 if (pawtab%nproju==2) nk=6
 factnk(1)=one;factnk(2)=one;factnk(3)=one
 factnk(4)=two;factnk(5)=two;factnk(6)=two
 indn(1,1)=1;indn(1,2)=4;indn(1,3)=5
 indn(2,1)=4;indn(2,2)=2;indn(2,3)=6
 indn(3,1)=5;indn(3,2)=6;indn(3,3)=3

!======================================================
!Compute local exact exchange Energy
!-----------------------------------------------------
 eexextemp=zero

 do ispden=1,pawrhoij%nspden
   jrhoij=1
   do irhoij=1,pawrhoij%nrhoijsel
     klmn=pawrhoij%rhoijselect(irhoij)
     if(pawtab%indklmn(3,klmn)==0.and.pawtab%indklmn(4,klmn)==2*lexexch) then
       m11=pawtab%klmntomn(1,klmn);m21=pawtab%klmntomn(2,klmn)
       n1=pawtab%klmntomn(3,klmn);n2=pawtab%klmntomn(4,klmn)
       nn1=(n1*n2)/2+1
       jrhoij1=1
       do irhoij1=1,pawrhoij%nrhoijsel
         klmn1=pawrhoij%rhoijselect(irhoij1)
         if(pawtab%indklmn(3,klmn1)==0.and.pawtab%indklmn(4,klmn1)==2*lexexch) then
           m31=pawtab%klmntomn(1,klmn1);m41=pawtab%klmntomn(2,klmn1)
           n3=pawtab%klmntomn(3,klmn1);n4=pawtab%klmntomn(4,klmn1)
           nn2=(n3*n4)/2+1
           do ll=1,lexexch+1
!            eexextemp=eexextemp-pawtab%vex(m11,m31,m41,m21,ll)*factnk(indn(nn1,nn2))*pawtab%fk(indn(nn1,nn2),ll)&
!            &        *pawrhoij%rhoijp(jrhoij,ispden)*pawrhoij%rhoijp(jrhoij1,ispden)
             eexextemp=eexextemp-pawtab%vex(m11,m31,m41,m21,ll)*pawtab%dltij(klmn)*pawtab%fk(indn(nn1,nn2),ll)&
&             *pawtab%dltij(klmn1)*pawrhoij%rhoijp(jrhoij,ispden)*pawrhoij%rhoijp(jrhoij1,ispden)
           end do
         end if
         jrhoij1=jrhoij1+pawrhoij%cplex
       end do !irhoij1
     end if
     jrhoij=jrhoij+pawrhoij%cplex
   end do !irhoij
 end do ! ispden
 eexextemp=eexextemp/two
 eexex=eexex+eexextemp*pawtab%exchmix

 if (abs(pawprtvol)>=2) then
   write(message, '(a)' )"   Contributions to the direct expression of energy:"
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a,e20.10,a)') "     HF exchange energy  =",eexextemp,ch10
   call wrtout(std_out,message,'COLL')
 end if

 DBG_EXIT("COLL")

 end subroutine pawxenergy
!!***
