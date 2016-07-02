!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkqptequiv
!!
!! NAME
!! mkqptequiv
!!
!! FUNCTION
!! This routine determines the equivalence between 
!!   1) qpoints and fermi surface kpoints
!!   2) qpoints under symmetry operations
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   Cryst<crystal_t>=Info on unit cell and symmetries.
!!   kpt_phon = fermi surface kpoints
!!   nkpt_phon = number of kpoints in the full FS set
!!   nqpt = number of qpoints
!!   qpt_full = qpoint coordinates
!!
!! OUTPUT
!!   FSfullpqtofull = mapping of k + q onto k' for k and k' in full BZ
!!   qpttoqpt(itim,isym,iqpt) = qpoint index which transforms to iqpt under isym and with time reversal itim.
!!
!! NOTES
!!   REMOVED 3/6/2008: much too large matrix, and not used at present
!!       FStoqpt = mapping of kpoint pairs (1 irreducible and 1 full) to qpoints
!!
!! PARENTS
!!      elphon,get_tau_k
!!
!! CHILDREN
!!      destroy_kptrank,get_rank_1kpt,mkkptrank,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkqptequiv(FSfullpqtofull,Cryst,kpt_phon,nkpt_phon,nqpt,qpttoqpt,qpt_full,mqtofull)

 use defs_basis
 use defs_elphon
 use m_kptrank
 use m_errors
 use m_profiling_abi

 use m_crystal,    only : crystal_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkqptequiv'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt_phon,nqpt
 type(crystal_t),intent(in) :: Cryst
!arrays
 integer,intent(out) :: FSfullpqtofull(nkpt_phon,nqpt),qpttoqpt(2,Cryst%nsym,nqpt)
 integer,intent(out),optional :: mqtofull(nqpt)
 real(dp),intent(in) :: kpt_phon(3,nkpt_phon),qpt_full(3,nqpt)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,iFSqpt,iqpt,isym,symrankkpt_phon
 character(len=500) :: message
 type(kptrank_type) :: kptrank_t
!arrays
 real(dp) :: tmpkpt(3),gamma_kpt(3)

! *************************************************************************

 call wrtout(std_out,' mkqptequiv : making rankkpt_phon and invrankkpt_phon',"COLL")

 call mkkptrank (kpt_phon,nkpt_phon,kptrank_t)

 FSfullpqtofull = -999
 gamma_kpt(:) = zero

 do ikpt_phon=1,nkpt_phon
   do iqpt=1,nqpt
!    tmpkpt = jkpt = ikpt + qpt
     tmpkpt(:) = kpt_phon(:,ikpt_phon) + qpt_full(:,iqpt)

!    which kpt is it among the full FS kpts?
     call get_rank_1kpt (tmpkpt,symrankkpt_phon,kptrank_t)

     FSfullpqtofull(ikpt_phon,iqpt) = kptrank_t%invrank(symrankkpt_phon)
     if (FSfullpqtofull(ikpt_phon,iqpt) == -1) then
       MSG_ERROR("looks like no kpoint equiv to k+q !!!")
     end if

   end do
 end do

 if (present(mqtofull)) then
   do iqpt=1,nqpt
     tmpkpt(:) = gamma_kpt(:) - qpt_full(:,iqpt)

!    which kpt is it among the full FS kpts?
     call get_rank_1kpt (tmpkpt,symrankkpt_phon,kptrank_t)

     mqtofull(iqpt) = kptrank_t%invrank(symrankkpt_phon)
     if (mqtofull(iqpt) == -1) then
       MSG_ERROR("looks like no kpoint equiv to -q !!!")
     end if
   end do
 end if

 call destroy_kptrank (kptrank_t)

!start over with q grid
 call wrtout(std_out,' mkqptequiv : FSfullpqtofull made. Do qpttoqpt',"COLL")

 call mkkptrank (qpt_full,nqpt,kptrank_t)

 qpttoqpt(:,:,:) = -1
 do iFSqpt=1,nqpt
   do isym=1,Cryst%nsym
     tmpkpt(:) =  Cryst%symrec(:,1,isym)*qpt_full(1,iFSqpt) &
&     + Cryst%symrec(:,2,isym)*qpt_full(2,iFSqpt) &
&     + Cryst%symrec(:,3,isym)*qpt_full(3,iFSqpt)

     call get_rank_1kpt (tmpkpt,symrankkpt_phon,kptrank_t)
     if (kptrank_t%invrank(symrankkpt_phon) == -1) then
       message = "looks like no kpoint equiv to q by symmetry without time reversal!!!"
       MSG_ERROR(message)
     end if
     qpttoqpt(1,isym,kptrank_t%invrank(symrankkpt_phon)) = iFSqpt

     tmpkpt = -tmpkpt
     call get_rank_1kpt (tmpkpt,symrankkpt_phon,kptrank_t)
     if (kptrank_t%invrank(symrankkpt_phon) == -1) then
       message = ' mkqptequiv : Error : looks like no kpoint equiv to q by symmetry with time reversal!!!'
       MSG_ERROR(message)
     end if
     qpttoqpt(2,isym,kptrank_t%invrank(symrankkpt_phon)) = iFSqpt
   end do
 end do

 call destroy_kptrank (kptrank_t)

end subroutine mkqptequiv
!!***
