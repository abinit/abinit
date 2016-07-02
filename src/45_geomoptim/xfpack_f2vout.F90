!{\src2tex{textfont=tt}}
!!****f* ABINIT/xfpack_f2vout
!! NAME
!! xfpack_f2vout
!!
!! FUNCTION
!! Old option=3, transfer fred and strten to vout
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
!! natom=number of atoms in cell
!! ndim=dimension of vout arrays
!! optcell=option for the optimisation of the unit cell. Described in abinit_help.
!!  Depending on its value, different part of strten
!!  are contained in vout.
!! strtarget(6)=target stresses ; they will be subtracted from strten when vout
!!  is computed.
!! ucvol=unit cell volume (bohr^3), needed for some values of optcell.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output variables
!! fred(3,natom)=grads of Etot wrt reduced coordinates (hartree)
!! strten(6)=components of the stress tensor (hartree/bohr^3)
!! vout(ndim)=vector that contains fred and some quantity derived from
!!   strten, depending on the value of optcell, and taking care ot strtarget
!!
!! PARENTS
!!      pred_bfgs,pred_delocint,pred_verlet,xfh_recover_deloc,xfh_recover_new
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine xfpack_f2vout(fred,natom,ndim,optcell,strtarget,strten,ucvol,vout)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xfpack_f2vout'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ndim,optcell
 real(dp),intent(in) :: ucvol
!arrays
 real(dp),intent(in) :: strtarget(6)
 real(dp),intent(in) :: fred(3,natom),strten(6)
 real(dp),intent(out) :: vout(ndim)

!Local variables-------------------------------
!scalars
 real(dp) :: strdiag
 character(len=500) :: message
!arrays
 real(dp) :: dstr(6)

! *************************************************************************

!!DEBUG
!write(ab_out,*) ''
!write(ab_out,*) 'xfpack_f2vout'
!write(ab_out,*) 'natom=',natom
!write(ab_out,*) 'ndim=',ndim
!write(ab_out,*) 'optcell=',optcell
!write(ab_out,*) 'ucvol=',ucvol
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

!
!Get vout from fred and strten
!
 vout(1:3*natom)= reshape(fred(:,:), (/3*natom/) )
 dstr(:)=strten(:)-strtarget(:)

 if(optcell==1)then

   vout(3*natom+1)=( dstr(1)+dstr(2)+dstr(3))*ucvol

 else if(optcell==2 .or. optcell==3 .or. optcell>=7)then

!  Eventually take away the trace
   strdiag=0.0_dp
   if(optcell==3) strdiag=(dstr(1)+dstr(2)+dstr(3))/3.0_dp
   if(optcell==2 .or. optcell==3)then
     vout(3*natom+1:3*natom+3)=(dstr(1:3)-strdiag)*ucvol
!    For non-diagonal derivatives, must take into account
!    that eps(i,j) AND eps(j,i) are varied at the same time. Thus, derivative
!    is twice larger
     vout(3*natom+4:3*natom+6)=dstr(4:6)*ucvol*2.0_dp
   else if(optcell==7 .or. optcell==8 .or. optcell==9)then
!    Similar to case optcell==2 or optcell==3, but in 2 dimensions.
     vout(3*natom+1:3*natom+3)=dstr(1:3)*ucvol
     vout(3*natom+optcell-6)  =dstr(optcell-3)*ucvol*2.0_dp
   end if

 else if(optcell==4 .or. optcell==5 .or. optcell==6)then

   vout(3*natom+1)=dstr(optcell-3)*ucvol

 end if

end subroutine xfpack_f2vout
!!***
