!{\src2tex{textfont=tt}}
!!****f* ABINIT/xfpack_x2vin
!! NAME
!! xfpack_x2vin
!!
!! FUNCTION
!! Old option=1, transfer xred, acell, and rprim to vin
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG,MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell0(3)=reference length scales of primitive translations (bohr), needed
!!   for some values of optcell.
!! natom=number of atoms in cell
!! ndim=dimension of vin arrays
!! nsym=order of group.
!! rprimd0(3,3)=reference real space primitive translations,
!!   needed for some values of optcell.
!! optcell=option for the optimisation of the unit cell. Described in abinit_help.
!!  Depending on its value, different part of acell and rprim
!!  are contained in vin.
!! symrel(3,3,nsym)=symmetry operators in terms of action on primitive translations
!! ucvol=unit cell volume (bohr^3), needed for some values of optcell.
!! ucvol0=reference unit cell volume (bohr^3), needed for some values of optcell.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output variables
!! acell(3)=length scales of primitive translations (bohr)
!! rprim(3,3)=dimensionless real space primitive translations
!! vin(ndim)=vector that contains xred and some quantity derived
!!   from acell and rprim, depending on the value of optcell.
!! xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! PARENTS
!!      pred_bfgs,pred_delocint,pred_lbfgs,pred_verlet,xfh_recover_deloc
!!      xfh_recover_new
!!
!! CHILDREN
!!      matr3inv,metric,mkrdim,strainsym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine xfpack_x2vin(acell,acell0,natom,ndim,nsym,optcell,&
  & rprim,rprimd0,symrel,ucvol,ucvol0,vin,xred)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xfpack_x2vin'
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_45_geomoptim, except_this_one => xfpack_x2vin
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ndim,nsym,optcell
 real(dp),intent(in) :: ucvol0
 real(dp),intent(inout) :: ucvol !vz_i
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: acell0(3),rprimd0(3,3)
 real(dp),intent(in) :: acell(3),rprim(3,3)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(out) :: vin(ndim)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk
 real(dp) :: scale
 character(len=500) :: message
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),gprimd0(3,3),rmet(3,3),rprimd(3,3)
 real(dp) :: rprimd_symm(3,3),scaling(3,3)

! *************************************************************************

!!DEBUG
!write(ab_out,*) ''
!write(ab_out,*) 'xfpack_x2vin'
!write(ab_out,*) 'natom=',natom
!write(ab_out,*) 'ndim=',ndim
!write(ab_out,*) 'nsym=',nsym
!write(ab_out,*) 'optcell=',optcell
!write(ab_out,*) 'ucvol=',ucvol
!write(ab_out,*) 'xred='
!do kk=1,natom
!write(ab_out,*) xred(:,kk)
!end do
!write(ab_out,*) 'VECTOR INPUT (vin) xfpack_x2vin INPUT'
!do ii=1,ndim,3
!if (ii+2<=ndim)then
!write(ab_out,*) ii,vin(ii:ii+2)
!else
!write(ab_out,*) ii,vin(ii:ndim)
!end if
!end do
!!DEBUG


!##########################################################
!### 1. Test for compatible ndim

 if(optcell==0 .and. ndim/=3*natom)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=0, ndim MUST be equal to 3*natom,',ch10,&
&   '  while ndim=',ndim,' and 3*natom=',3*natom,'.'
   MSG_BUG(message)
 end if

 if( (optcell==1 .or. optcell==4 .or. optcell==5 .or. optcell==6) &
& .and. ndim/=3*natom+1)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=1,4,5 or 6, ndim MUST be equal to 3*natom+1,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+1=',3*natom+1,'.'
   MSG_BUG(message)
 end if

 if( (optcell==2 .or. optcell==3) &
& .and. ndim/=3*natom+6)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=2 or 3, ndim MUST be equal to 3*natom+6,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+6=',3*natom+6,'.'
   MSG_BUG(message)
 end if

 if( optcell>=7 .and. ndim/=3*natom+3)then
   write(message,'(a,a,a,i4,a,i4,a)' )&
&   '  When optcell=7,8 or 9, ndim MUST be equal to 3*natom+3,',ch10,&
&   '  while ndim=',ndim,' and 3*natom+3=',3*natom+3,'.'
   MSG_BUG(message)
 end if

!##########################################################
!### 2. option=1, transfer xred, acell, and rprim to vin

!Get vin from xred, acell, and rprim
 vin(1:3*natom)= reshape(xred(:,:), (/3*natom/) )

 if(optcell/=0)then
   call mkrdim(acell,rprim,rprimd)
   call strainsym(nsym,rprimd0,rprimd,rprimd_symm,symrel)
   call metric(gmet,gprimd,-1,rmet,rprimd_symm,ucvol)

   if(optcell==1)then

!    vin(3*natom+1)=ucvol**third
     vin(3*natom+1)=(ucvol/ucvol0)**third

   else if(optcell==2 .or. optcell==3 .or. optcell>=7)then

!    Generates gprimd0
     call matr3inv(rprimd0,gprimd0)
     do ii=1,3
       do jj=1,3
         scaling(ii,jj)=0.0_dp
         do kk=1,3
           scaling(ii,jj)=scaling(ii,jj)+rprimd_symm(ii,kk)*gprimd0(jj,kk)
         end do
       end do
     end do
!    Rescale if the volume must be preserved
     if(optcell==3)then
       scale=(ucvol0/ucvol)**third
       scaling(:,:)=scale*scaling(:,:)
     end if
     if(optcell==2 .or. optcell==3)then
       vin(3*natom+1)=scaling(1,1) ; vin(3*natom+4)=(scaling(2,3)+scaling(3,2))*0.5_dp
       vin(3*natom+2)=scaling(2,2) ; vin(3*natom+5)=(scaling(1,3)+scaling(3,1))*0.5_dp
       vin(3*natom+3)=scaling(3,3) ; vin(3*natom+6)=(scaling(1,2)+scaling(2,1))*0.5_dp
     else if(optcell>=7)then
       vin(3*natom+1)=scaling(1,1)
       vin(3*natom+2)=scaling(2,2)
       vin(3*natom+3)=scaling(3,3)
       if(optcell==7)vin(3*natom+1)=(scaling(2,3)+scaling(3,2))*0.5_dp
       if(optcell==8)vin(3*natom+2)=(scaling(1,3)+scaling(3,1))*0.5_dp
       if(optcell==9)vin(3*natom+3)=(scaling(1,2)+scaling(2,1))*0.5_dp
     end if

   else if(optcell==4 .or. optcell==5 .or. optcell==6)then

     vin(3*natom+1)=acell(optcell-3)/acell0(optcell-3)

   end if

 end if


end subroutine xfpack_x2vin
!!***
