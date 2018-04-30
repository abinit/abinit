!{\src2tex{textfont=tt}}
!!****f* ABINIT/integrate_gamma
!!
!! NAME
!! integrate_gamma
!!
!! FUNCTION
!! This routine integrates the electron phonon coupling matrix
!! over the kpoints on the fermi surface. A dependency on qpoint
!! remains for gamma_qpt
!!
!! COPYRIGHT
!! Copyright (C) 2004-2018 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!      elph_ds%qpt_full = qpoint coordinates
!!      elph_ds%nqptirred = number of irred qpoints
!!      elph_ds%qirredtofull = indexing of the GKK qpoints found
!!   FSfullpqtofull = mapping of k+q to k
!!
!! OUTPUT
!!   elph_ds = modified elph_ds%gamma_qpt and created elph_ds%gamma_rpt
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      get_rank_1kpt,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine integrate_gamma(elph_ds,FSfullpqtofull)

 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_kptrank
 use m_errors
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'integrate_gamma'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)

!Local variables-------------------------------
!scalars
 integer :: comm,ikpt_phon,ikpt_phonq,ib1,ib2,ibeff,iqpt,iqpt_fullbz,isppol,ierr
 integer :: irec, symrankkpt_phon,nbranch,nsppol,ngkkband, ik_this_proc
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 real(dp),allocatable :: tmp_gkk(:,:,:,:)

! *************************************************************************

 comm = xmpi_world

 write (message,'(3a)')ch10,' entering integrate_gamma ',ch10
 call wrtout(std_out,message,'COLL')

 nsppol   = elph_ds%nsppol
 nbranch  = elph_ds%nbranch
 ngkkband = elph_ds%ngkkband

 ABI_ALLOCATE(elph_ds%gamma_qpt,(2,nbranch**2,nsppol,elph_ds%nqpt_full))
 elph_ds%gamma_qpt = zero

 ABI_ALLOCATE(tmp_gkk ,(2,ngkkband**2,nbranch**2,nsppol))

 if (elph_ds%gkqwrite == 0) then
   call wrtout(std_out,' integrate_gamma : keeping gamma matrices in memory','COLL')
 else if (elph_ds%gkqwrite == 1) then
   fname=trim(elph_ds%elph_base_name) // '_GKKQ'
   write (message,'(2a)')' integrate_gamma : reading gamma matrices from file ',trim(fname)
   call wrtout(std_out,message,'COLL')
 else
   write (message,'(a,i0)')' Wrong value for gkqwrite = ',elph_ds%gkqwrite
   MSG_BUG(message)
 end if



 do iqpt=1,elph_ds%nqptirred
   iqpt_fullbz = elph_ds%qirredtofull(iqpt)
   call get_rank_1kpt (elph_ds%k_phon%kpt(:,iqpt_fullbz),symrankkpt_phon, elph_ds%k_phon%kptrank_t)
   write (std_out,*) ' iqpt_fullbz in qpt grid only,  rank ', iqpt_fullbz, symrankkpt_phon

   do ik_this_proc =1,elph_ds%k_phon%my_nkpt
     ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

     if (elph_ds%gkqwrite == 0) then
       tmp_gkk = elph_ds%gkk_qpt(:,:,:,ik_this_proc,:,iqpt)
     else if (elph_ds%gkqwrite == 1) then
       irec = (iqpt-1)*elph_ds%k_phon%my_nkpt+ik_this_proc
       if (ikpt_phon == 1) then
         write (std_out,*) ' integrate_gamma  read record ', irec
       end if
       read (elph_ds%unitgkq,REC=irec) tmp_gkk(:,:,:,:)
     end if

     do isppol=1,nsppol
       ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)
!      
       do ib1=1,ngkkband
         do ib2=1,ngkkband
           ibeff = ib2+(ib1-1)*ngkkband
           elph_ds%gamma_qpt(:,:,isppol,iqpt_fullbz) = elph_ds%gamma_qpt(:,:,isppol,iqpt_fullbz) + &
&           tmp_gkk(:,ibeff,:,isppol)&
&           *elph_ds%gkk_intweight(ib1,ikpt_phon,isppol)*elph_ds%gkk_intweight(ib2,ikpt_phonq,isppol)
!          NOTE: if ngkkband==1 we are using trivial weights since average
!          over bands was done in normsq_gkk (nmsq_gam_sumFS or nmsq_pure_gkk)
         end do ! ib2
       end do ! ib1
     end do ! isppol
   end do ! ikpt_phon
 end do ! iqpt

 call xmpi_sum (elph_ds%gamma_qpt, comm, ierr)

 ABI_DEALLOCATE(tmp_gkk)

!need prefactor of 1/nkpt for each integration over 1 kpoint index. NOT INCLUDED IN elph_ds%gkk_intweight
 do iqpt=1,elph_ds%nqptirred
   iqpt_fullbz = elph_ds%qirredtofull(iqpt)
!  elph_ds%gamma_qpt(:,:,:,iqpt_fullbz) = elph_ds%gamma_qpt(:,:,:,iqpt_fullbz) / elph_ds%k_phon%nkpt / n0(1) / n0(1)
!  elph_ds%gamma_qpt(:,:,:,iqpt_fullbz) = elph_ds%gamma_qpt(:,:,:,iqpt_fullbz) / elph_ds%k_phon%nkpt / elph_ds%k_phon%nkpt
   elph_ds%gamma_qpt(:,:,:,iqpt_fullbz) = elph_ds%gamma_qpt(:,:,:,iqpt_fullbz) * elph_ds%occ_factor / elph_ds%k_phon%nkpt
 end do

 call wrtout(std_out,' integrate_gamma: gamma matrices have been calculated for recip space and irred qpoints ',"COLL")

end subroutine integrate_gamma
!!***
