!{\src2tex{textfont=tt}}
!!****f* ABINIT/defad
!! NAME
!! defad
!!
!! FUNCTION
!! Initialisation of aim input variables to their
!! default values.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (no input : initialisation by default values)
!!
!! OUTPUT
!! aim_dtset = the structured entity containing all input variables
!!
!! PARENTS
!!      aim
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine defad(aim_dtset)

 use defs_basis
 use defs_aimprom
 use defs_abitypes
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'defad'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type),intent(out) :: aim_dtset

!Local variables ------------------------------

! *********************************************************************

 aim_dtset%isurf=0
 aim_dtset%crit=0
 aim_dtset%irsur=0
 aim_dtset%foll=0
 aim_dtset%irho=0
 aim_dtset%ivol=0
 aim_dtset%denout=0
 aim_dtset%lapout=0
 aim_dtset%gpsurf=0
 aim_dtset%plden=0
 aim_dtset%dltyp=0

 aim_dtset%batom=1
 aim_dtset%nsa=3
 aim_dtset%nsb=3
 aim_dtset%nsc=3
 aim_dtset%npt=100
 aim_dtset%nth=32
 aim_dtset%nph=48

 aim_dtset%themax=pi
 aim_dtset%themin=zero
 aim_dtset%phimin=zero
 aim_dtset%phimax=two_pi
 aim_dtset%phi0=zero
 aim_dtset%th0=zero
 aim_dtset%folstp=5.d-2
 aim_dtset%dr0=5.d-2
 aim_dtset%atrad=one
 aim_dtset%rmin=one

 aim_dtset%foldep(:)=zero
 aim_dtset%vpts(:,:)=zero
 aim_dtset%ngrid(:)=30
 aim_dtset%scal(:)=one
 aim_dtset%maxatd=1.d1
 aim_dtset%maxcpd=5.d1

 aim_dtset%dpclim=1.d-2
 aim_dtset%lstep=1.d-10
 aim_dtset%lstep2=1.d-5
 aim_dtset%lgrad=1.d-12
 aim_dtset%lgrad2=1.d-5
 aim_dtset%coff1=0.98_dp
 aim_dtset%coff2=0.95_dp

end subroutine defad
!!***
