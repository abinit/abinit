!{\src2tex{textfont=tt}}
!!****f* ABINIT/chiscwrt 
!! NAME
!!  chiscwrt 
!!
!! FUNCTION
!!  Distributes values of chi_org on chi_sc according to ion-ion distances and multiplicities in shells
!!
!! INPUTS
!!  chi_org  = response matrix in original unit cell
!!  disv_org = distances (multiplied by magntization) in original unit cell
!!  nat_org = number of atoms in original unit cell
!!  sdisv_org = radii of shells in original unit cell
!!  smult_org = multiplicities of shells in original unit cell
!!  nsh_org   = number of shells in original unit cell
!!  disv_sc   = distances (multiplied by magntization) in super-cell
!!  nat_sc  = number of atoms in super-cell
!!  sdisv_sc = radii of shells in super-cell (was unused, so suppressed)
!!  smult_sc = multiplicities of shells in super-cell
!!  nsh_sc = number of shells in super-cell
!!  opt =
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2016 ABINIT group (DJA)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! OUTPUT
!!  chi_sc = first line of reponse matrix in super-cell
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      pawuj_det
!!
!! CHILDREN
!!      prmat,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine chiscwrt(chi_org,disv_org,nat_org,sdisv_org,smult_org,nsh_org,chi_sc,&
& disv_sc,nat_sc,smult_sc,nsh_sc,opt,prtvol) 

 use defs_basis
 use m_profiling_abi

 use m_pptools,    only : prmat

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chiscwrt'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)              :: nat_org,nsh_org,nsh_sc 
 integer,intent(in)              :: nat_sc
 integer,intent(in),optional     :: prtvol                  
 integer,intent(in),optional     :: opt
!arrays
 real(dp),intent(in)             :: chi_org(nat_org),disv_org(nat_org),sdisv_org(nsh_org)
 integer,intent(in)              :: smult_org(nsh_org), smult_sc(nsh_sc)
 real(dp),intent(in)             :: disv_sc(nat_sc)
 real(dp),intent(out)            :: chi_sc(nat_sc)

!Local variables-------------------------------
!scalars
 integer                      :: iatom,jatom,jsh,optt
 character(len=500)           :: message
!arrays
 real(dp)                     :: chi_orgl(nat_org)

! *************************************************************************

 if (present(opt)) then
   optt=opt
 else
   optt=1
 end if


 do iatom=1,nat_org
   do jsh=1,nsh_org
     if (disv_org(iatom)==sdisv_org(jsh)) then
       if (opt==2) then
         chi_orgl(iatom)=chi_org(iatom)
       else if (opt==1) then
         chi_orgl(iatom)=chi_org(iatom)*smult_org(jsh)/smult_sc(jsh)
       end if 
       exit
     end if
   end do !iatom
 end do  !jsh
 
 if (prtvol>1) then
   write(message,fmt='(a,150f10.5)')' chiscwrt: chi at input ',chi_org
   call wrtout(std_out,message,'COLL')
   write(message,fmt='(a,150f10.5)')' chiscwrt: chi after division ',chi_orgl
   call wrtout(std_out,message,'COLL')
 end if

 do iatom=1,nat_sc
   do jatom=1,nat_org
     if (disv_org(jatom)==disv_sc(iatom)) then
       chi_sc(iatom)=chi_orgl(jatom)
       exit
     else if (jatom==nat_org) then
       chi_sc(iatom)=0_dp
     end if
   end do
 end do

 if (prtvol>1) then
   write(message,'(a)')' chiscwrt, chi_sc '
   call wrtout(std_out,message,'COLL')
   call prmat(chi_sc,1,nat_sc,1,std_out)
 end if

end subroutine chiscwrt
!!***
