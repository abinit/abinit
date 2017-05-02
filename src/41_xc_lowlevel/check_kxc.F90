!{\src2tex{textfont=tt}}
!!****f* ABINIT/check_kxc
!! NAME
!! check_kxc
!!
!! FUNCTION
!!  Given a XC functional (defined by ixc), check if Kxc (dVxc/drho) is avalaible.
!!
!! COPYRIGHT
!! Copyright (C) 2012-2017 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ixc = internal code for xc functional
!!  optdriver=type of calculation (ground-state, response function, GW, ...)
!!
!! OUTPUT
!!
!! PARENTS
!!      respfn,scfcv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine check_kxc(ixc,optdriver)

 use defs_basis
 use m_errors
 use libxc_functionals

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_kxc'
!End of the abilint section

 implicit none

!Arguments -------------------------------
 integer, intent(in) :: ixc,optdriver

!Local variables -------------------------
 logical :: kxc_available
 character(len=500) :: msg

! *********************************************************************

 kxc_available=.false.

 if (ixc>=0) then
   kxc_available=(ixc/=16.and.ixc/=17.and.ixc/=26.and.ixc/=27)
   if (.not.kxc_available) then
     write(msg,'(a,i0,3a)') &
&     'The selected XC functional (ixc=',ixc,')',ch10,&
&     'does not provide Kxc (dVxc/drho) !'
   end if
 else if (ixc==-406.or.ixc==-427.or.ixc==-428.or.ixc==-456)then 
   kxc_available=.true.
 else ! ixc<0 and not one of the allowed hybrids
   kxc_available=libxc_functionals_has_kxc()
   if (.not.kxc_available) then
     write(msg,'(a,i0,7a)') &
&     'The selected XC functional (ixc=',ixc,'):',ch10,&
&     '   <<',trim(libxc_functionals_fullname()),'>>',ch10,&
&     'does not provide Kxc (dVxc/drho) !'
   end if
 end if

 if (.not.kxc_available) then
   write(msg,'(7a)') trim(msg),ch10,&
&   'However, with the current input options, ABINIT needs Kxc.',ch10,&
&   '>Possible action:',ch10,&
&   'Change the XC functional in psp file or input file.'
   if (optdriver==0) then
     write(msg,'(13a)') trim(msg),ch10,&
&     '>Possible action (2):',ch10,&
&     'If you are using density mixing for the SCF cycle',ch10,&
&     '(iscf>=10, which is the default for PAW),',ch10,&
&     'change to potential mixing (iscf=7, for instance).',ch10,&
&     '>Possible action (3):',ch10,&
&     'Switch to another value of densfor_pred (=5, for instance).'
   end if
   MSG_ERROR(msg)
 end if

end subroutine check_kxc
!!***
