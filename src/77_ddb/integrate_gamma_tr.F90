!{\src2tex{textfont=tt}}
!!****f* ABINIT/integrate_gamma_tr
!!
!! NAME
!! integrate_gamma_tr
!!
!! FUNCTION
!! This routine integrates the TRANSPORT electron phonon coupling matrices
!! over the kpoints on the fermi surface. A dependency on qpoint
!! remains for gamma_qpt_in/out
!! Copied from integrate_gamma
!!
!! COPYRIGHT
!! Copyright (C) 2004-2017 ABINIT group (BXu,MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!      elph_ds%qpt_full = qpoint coordinates
!!   FSfullpqtofull = mapping of k+q to k
!!   veloc_sq1 = mean square electronic velocity on constant energy surface
!!   veloc_sq2 = mean square electronic velocity on constant energy surface
!!
!! OUTPUT
!!   elph_tr_ds%gamma_qpt_tr and created elph_tr_ds%gamma_rpt_tr
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine integrate_gamma_tr(elph_ds,FSfullpqtofull,s1,s2, &
&                             veloc_sq1,veloc_sq2,elph_tr_ds)

 use defs_basis
 use defs_elphon
 use m_errors
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'integrate_gamma_tr'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: s1,s2
 type(elph_tr_type), intent(inout) :: elph_tr_ds
 type(elph_type),intent(in) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 real(dp),intent(in) :: veloc_sq1(3,elph_ds%nsppol), veloc_sq2(3,elph_ds%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ikpt_phonq,ib1,ib2,ibeff,ierr,iqpt,iqpt_fullbz,isppol
 integer :: itensor, icomp, jcomp,comm
 integer :: fib1, fib2
 integer :: ik_this_proc
! integer :: ikpttemp
 character(len=500) :: message
 real(dp) :: wtk, wtkpq, interm
 real(dp) :: veloc1_i, veloc1_j, veloc2_i, veloc2_j
!arrays
 real(dp) :: elvelock(3), elvelockpq(3)
 real(dp) :: velocwtk(3), velocwtkpq(3)
 real(dp) :: vvelocwtk(3,3), vvelocwtkpq(3,3)
 real(dp),allocatable :: tmp_gkk(:,:,:,:)

! *************************************************************************

 comm = xmpi_world

!information
 if (elph_ds%gkqwrite == 0) then
   write (message,'(a)')' integrate_gamma_tr : keeping gamma matrices in memory'
   call wrtout(std_out,message,'COLL')
 else if (elph_ds%gkqwrite == 1) then
   write (message,'(a)')' integrate_gamma_tr : reading gamma matrices from disk'
   call wrtout(std_out,message,'COLL')
 else
   write (message,'(3a,i3)')' integrate_gamma_tr : BUG-',ch10,&
&   ' Wrong value for gkqwrite = ',elph_ds%gkqwrite
   MSG_BUG(message)
 end if

!allocate temp variables
 ABI_STAT_ALLOCATE(tmp_gkk,(2,elph_ds%ngkkband**2,elph_ds%nbranch**2,elph_ds%nsppol), ierr)
 ABI_CHECK(ierr==0, 'trying to allocate array tmp_gkkout')

 do iqpt=1,elph_ds%nqptirred
   iqpt_fullbz = elph_ds%qirredtofull(iqpt)
!  write(std_out,*)'iqpt, iqptfullbz  ',iqpt, iqpt_fullbz

   do ik_this_proc =1,elph_ds%k_phon%my_nkpt
     ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

     if (elph_ds%gkqwrite == 0) then
       tmp_gkk = elph_ds%gkk_qpt(:,:,:,ik_this_proc,:,iqpt)
     else if (elph_ds%gkqwrite == 1) then
       read(elph_ds%unitgkq,REC=((iqpt-1)*elph_ds%k_phon%my_nkpt+ik_this_proc)) tmp_gkk
     end if

     ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)

     do isppol=1,elph_ds%nsppol
       do ib1=1,elph_ds%ngkkband !FS bands
         fib1=ib1+elph_ds%minFSband-1 ! full bands
         elvelock(:)=elph_tr_ds%el_veloc(ikpt_phon,fib1,:,isppol)
         wtk=elph_tr_ds%tmp_gkk_intweight1(ib1,ikpt_phon,isppol)
         velocwtk(:)=elph_tr_ds%tmp_velocwtk1(ib1,ikpt_phon,:,isppol)
         vvelocwtk(:,:)=elph_tr_ds%tmp_vvelocwtk1(ib1,ikpt_phon,:,:,isppol)

         do ib2=1,elph_ds%ngkkband ! FS bands
           ibeff=ib2+(ib1-1)*elph_ds%ngkkband ! full bands
           fib2=ib2+elph_ds%minFSband-1
           elvelockpq(:)= elph_tr_ds%el_veloc(ikpt_phonq,fib2,:,isppol)
           wtkpq=elph_tr_ds%tmp_gkk_intweight2(ib2,ikpt_phonq,isppol)
           velocwtkpq(:)=elph_tr_ds%tmp_velocwtk2(ib2,ikpt_phonq,:,isppol)
           vvelocwtkpq(:,:)=elph_tr_ds%tmp_vvelocwtk2(ib2,ikpt_phonq,:,:,isppol)

!          MJV 31/03/2009: Note that the following is valid for any geometry, not just cubic!
!          see eq 5 and 6 of prb 36 4103 (Al-Lehaibi et al 1987)
!          see also Allen PRB 17 3725
!          generalization to tensorial quantities is simple, by keeping the directional
!          references of velock and velockpq as indices.
           do icomp = 1, 3
             do jcomp = 1, 3
               itensor = (icomp-1)*3+jcomp
!              FIXME: could use symmetry i <-> j

               veloc1_i = sqrt(veloc_sq1(icomp,isppol))
               veloc1_j = sqrt(veloc_sq1(jcomp,isppol))
               veloc2_i = sqrt(veloc_sq2(icomp,isppol))
               veloc2_j = sqrt(veloc_sq2(jcomp,isppol))
               if (elph_ds%use_k_fine == 1) then
                 interm = vvelocwtk(icomp,jcomp)*wtkpq/veloc1_i/veloc1_j + &
&                 s1*s2*vvelocwtkpq(icomp,jcomp)*wtk/veloc2_i/veloc2_j - &
&                 s1*velocwtk(jcomp)*velocwtkpq(icomp)/veloc1_j/veloc2_i - &
&                 s2*velocwtk(icomp)*velocwtkpq(jcomp)/veloc1_i/veloc2_j

                 elph_tr_ds%gamma_qpt_tr(:,itensor,:,isppol,iqpt_fullbz) = &
&                 elph_tr_ds%gamma_qpt_tr(:,itensor,:,isppol,iqpt_fullbz) + &
&                 tmp_gkk(:,ibeff,:,isppol)*interm
               else
                 elph_tr_ds%gamma_qpt_tr(:,itensor,:,isppol,iqpt_fullbz) = &
&                 elph_tr_ds%gamma_qpt_tr(:,itensor,:,isppol,iqpt_fullbz) + &
&                 tmp_gkk(:,ibeff,:,isppol) &
&                 *(elvelock(icomp)/veloc1_i - s1*elvelockpq(icomp)/veloc2_i) &
&                 *(elvelock(jcomp)/veloc1_j - s2*elvelockpq(jcomp)/veloc2_j) &
&                 *wtk*wtkpq
               end if
             end do
           end do

         end do
       end do
     end do ! isppol

   end do ! ik
 end do ! iq

 call xmpi_sum (elph_tr_ds%gamma_qpt_tr, comm, ierr)

 ABI_DEALLOCATE(tmp_gkk)


!need prefactor of 1/nkpt for each integration over 1 kpoint index.
!NOT INCLUDED IN elph_ds%gkk_intweight
!Add a factor of 1/2 for the cross terms of (v-v')(v-v')
 elph_tr_ds%gamma_qpt_tr = elph_tr_ds%gamma_qpt_tr* elph_ds%occ_factor*0.5_dp / elph_ds%k_phon%nkpt

 write (message,'(2a)')' integrate_gamma_tr : transport gamma matrices are calculated ',&
& ' in recip space and for irred qpoints'
!call wrtout(std_out,message,'COLL')

end subroutine integrate_gamma_tr
!!***
