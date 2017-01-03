!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_full_kgrid
!! NAME
!! get_full_kgrid
!!
!! FUNCTION
!! create full grid of kpoints and find equivalent
!! irred ones. Duplicates work in getkgrid, but need all outputs
!! of kpt_fullbz, and indkpt
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (MVer,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  kpt(3,nkpt)=irreducible kpoints
!!  kptrlatt(3,3)=lattice vectors for full kpoint grid
!!  nkpt=number of irreducible kpoints
!!  nkpt_fullbz=number of kpoints in full brillouin zone
!!  nshiftk=number of kpoint grid shifts
!!  nsym=number of symmetries
!!  shiftk(3,nshiftk)=kpoint shifts
!!  symrel(3,3,nsym)=symmetry matrices in real space
!!
!! OUTPUT
!!  indkpt(nkpt_fullbz)=non-symmetrized indices of the k-points (see symkpt.f)
!!  kpt_fullbz(3,nkpt_fullbz)=kpoints in full brillouin zone
!!
!! NOTES
!!  MG: The present inplementation always assumes kptopt==1 !!!!
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!      destroy_kptrank,get_kpt_fullbz,get_rank_1kpt,mati3inv,mkkptrank
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine get_full_kgrid(indkpt,kpt,kpt_fullbz,kptrlatt,nkpt,&
& nkpt_fullbz,nshiftk,nsym,shiftk,symrel)

 use defs_basis
 use m_kptrank
 use m_profiling_abi
 use m_errors

 use m_numeric_tools,   only : wrap2_pmhalf

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_full_kgrid'
 use interfaces_32_util
 use interfaces_56_recipspace, except_this_one => get_full_kgrid
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nkpt_fullbz,nshiftk,nsym
!arrays
 integer,intent(in) :: kptrlatt(3,3),symrel(3,3,nsym)
 integer,intent(out) :: indkpt(nkpt_fullbz)
 real(dp),intent(in) :: kpt(3,nkpt),shiftk(3,nshiftk)
 real(dp),intent(out) :: kpt_fullbz(3,nkpt_fullbz)

!Local variables-------------------------------
!scalars
 integer :: ikpt,isym,itim,timrev
 integer :: symrankkpt
 character(len=500) :: message
 type(kptrank_type) :: kptrank_t

!arrays
 integer :: inv_symrel(3,3,nsym)
 real(dp) :: k2(3)

! *********************************************************************

!Invert symrels => gives symrels for kpoints

 do isym=1,nsym
   call mati3inv (symrel(:,:,isym),inv_symrel(:,:,isym))
 end do

 call get_kpt_fullbz(kpt_fullbz,kptrlatt,nkpt_fullbz,nshiftk,shiftk)

!make full k-point rank arrays
 call mkkptrank (kpt,nkpt,kptrank_t)

!
!find equivalence to irred kpoints in kpt
!
 indkpt(:) = 0
 timrev=1 ! includes the time inversion symmetry
 do ikpt=1,nkpt_fullbz
   do isym=1,nsym
     do itim=1,(1-2*timrev),-2

       k2(:) = itim*(inv_symrel(:,1,isym)*kpt_fullbz(1,ikpt) + &
&       inv_symrel(:,2,isym)*kpt_fullbz(2,ikpt) + &
&       inv_symrel(:,3,isym)*kpt_fullbz(3,ikpt))

       call get_rank_1kpt (k2,symrankkpt,kptrank_t)
       if (kptrank_t%invrank(symrankkpt) /= -1) indkpt(ikpt) = kptrank_t%invrank(symrankkpt)

     end do ! loop time reversal symmetry
   end do !  loop sym ops

   if (indkpt(ikpt) == 0) then
     write (message,'(a,i0)')' indkpt(ikpt) is still 0: no irred kpoint is equiv to ikpt ',ikpt
     MSG_BUG(message)
   end if
 end do !  loop full kpts

 call destroy_kptrank (kptrank_t)

end subroutine get_full_kgrid
!!***
