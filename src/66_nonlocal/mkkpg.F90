!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkkpg
!! NAME
!! mkkpg
!!
!! FUNCTION
!! Compute all (k+G) vectors (in reduced coordinates) for given k point.
!! Eventually compute related data
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  kg(3,npw)=integer coords of planewaves in basis sphere
!!  kpt(3)=k point in terms of recip. translations
!!  nkpg=second dimension of array kpg
!!  npw=number of plane waves in reciprocal space
!!
!! OUTPUT
!!  kpg(npw,3)= (k+G) components
!!  === if nkpg==9 ===
!!    kpg(npw,4:9)= [(k+G)_a].[(k+G)_b] quantities
!!
!! PARENTS
!!      ctocprj,d2frnl,debug_tools,dfpt_nstpaw,dfpt_nstwf,dfpt_rhofermi
!!      dfptnl_resp,forstrnps,getcprj,getgh1c,ks_ddiago,m_io_kss,m_shirley
!!      m_wfd,nonlop_ylm,prep_bandfft_tabs,vtorho,wfd_vnlpsi
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mkkpg(kg,kpg,kpt,nkpg,npw)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkkpg'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpg,npw
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(out) :: kpg(npw,nkpg)

!Local variables-------------------------------
!scalars
 integer :: ipw,mu,mua,mub
 character(len=500) :: message
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)

! *************************************************************************

 DBG_ENTER("COLL")

 if (nkpg==0) return

!-- Test nkpg --
 if (nkpg/=3.and.nkpg/=9) then
   write(message, '(a,i0)' )' Bad value for nkpg !',nkpg
   MSG_BUG(message)
 end if

!-- Compute (k+G) --
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP PRIVATE(mu,ipw)
 do ipw=1,npw
   do mu=1,3
     kpg(ipw,mu)=kpt(mu)+dble(kg(mu,ipw))
   end do
 end do
!$OMP END PARALLEL DO

!-- Compute [(k+G)_a].[(k+G)_b] --
 if (nkpg==9) then
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP PRIVATE(ipw,mu,mua,mub)
   do ipw=1,npw
     do mu=4,9
       mua=alpha(mu-3);mub=beta(mu-3)
       kpg(ipw,mu)=kpg(ipw,mua)*kpg(ipw,mub)
     end do
   end do
!$OMP END PARALLEL DO
 end if

 DBG_EXIT("COLL")

end subroutine mkkpg
!!***
