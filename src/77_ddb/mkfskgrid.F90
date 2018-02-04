!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkfskgrid
!!
!! NAME
!! mkfskgrid
!!
!! FUNCTION
!! This routine sets up the full FS kpt grid by symmetry
!!
!! COPYRIGHT
!! Copyright (C) 2004-2018 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nsym    = number of symmetries for the full system
!!  symrec  = reciprocal space symmetries (those for the kpts)
!!  timrev  = 1 if time reversal symmetry is to be used
!!
!! OUTPUT
!!  elph_k datastructure:
!!  elph_k%nkpt           = full number of kpoints close to the FS
!!  elph_k%kpt            = full set of kpoints close to the FS
!!  elph_k%wtkirr         = weights of the irreducible kpoints
!!  elph_k%kphon_irr2full = indices of irred kpoints in full array
!!
!! NOTES
!!  WARNING: supposes kpt grid has full symmetry!! Not always true!!!
!!    but should be for Monkhorst-Pack, efficient grids.
!!    otherwise you get an error message in interpolate_gkk because
!!    an FS kpt can not be found in the gkk file.
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      destroy_kptrank,get_rank_1kpt,mkkptrank,sort_int,wrap2_pmhalf,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkFSkgrid (elph_k, nsym, symrec, timrev)

 use defs_basis
 use defs_elphon
 use m_kptrank
 use m_profiling_abi
 use m_errors
 use m_sort

 use m_numeric_tools,   only : wrap2_pmhalf

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkFSkgrid'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym,timrev
 type(elph_kgrid_type),intent(inout) :: elph_k
!arrays
 integer,intent(in) :: symrec(3,3,nsym)

!Local variables-------------------------------
!scalars
 integer :: ikpt1,ikpt2,isym,itim,new,symrankkpt
 real(dp) :: timsign, res
 character(len=500) :: message

!arrays
 real(dp) :: kpt(3),redkpt(3)
 integer, allocatable :: sortindexing(:), rankallk(:)

 integer, allocatable :: tmpkphon_full2irr(:,:)
 real(dp), allocatable :: tmpkpt(:,:)

! *************************************************************************

 if(timrev /= 1 .and. timrev /= 0)then
   write (message,'(a,i0)')' timrev must be 1 or 0 but found timrev= ',timrev
   MSG_BUG(message)
 end if

 ABI_ALLOCATE(tmpkphon_full2irr,(3,2*elph_k%nkptirr*nsym))
 tmpkphon_full2irr = -1

 ABI_ALLOCATE(tmpkpt,(3,2*elph_k%nkptirr*nsym))

 ABI_ALLOCATE(elph_k%wtkirr,(elph_k%nkptirr))
 elph_k%wtkirr(:) = zero

!first allocation for irred kpoints - will be destroyed below
 call mkkptrank (elph_k%kptirr,elph_k%nkptirr,elph_k%kptrank_t)
 ABI_ALLOCATE(rankallk,(elph_k%kptrank_t%max_rank))

!elph_k%kptrank_t%invrank is used as a placeholder in the following loop
 rankallk = -1
 elph_k%kptrank_t%invrank = -1

!replicate all irred kpts by symmetry to get the full k grid.
 elph_k%nkpt=0 !zero k-points found so far
 do isym=1,nsym
   do itim=0,1
     timsign = one-two*itim
     do ikpt1=1,elph_k%nkptirr
!      generate symmetrics of kpt ikpt1
       kpt(:) = timsign*(symrec(:,1,isym)*elph_k%kptirr(1,ikpt1) + &
&       symrec(:,2,isym)*elph_k%kptirr(2,ikpt1) + &
&       symrec(:,3,isym)*elph_k%kptirr(3,ikpt1))
       
       call get_rank_1kpt (kpt,symrankkpt,elph_k%kptrank_t)

!      is the kpt on the full grid (may have lower symmetry than full spgroup)
!      is kpt among the full FS kpts found already?
       if (elph_k%kptrank_t%invrank(symrankkpt) == -1) then
         elph_k%wtkirr(ikpt1)=elph_k%wtkirr(ikpt1)+1
         elph_k%nkpt=elph_k%nkpt+1

         call wrap2_pmhalf(kpt(1),redkpt(1),res)
         call wrap2_pmhalf(kpt(2),redkpt(2),res)
         call wrap2_pmhalf(kpt(3),redkpt(3),res)
         tmpkpt(:,elph_k%nkpt) = redkpt
         tmpkphon_full2irr(1,elph_k%nkpt) = ikpt1
!        save sym that sends irred kpt ikpt1 onto full kpt
         tmpkphon_full2irr(2,elph_k%nkpt) = isym
         tmpkphon_full2irr(3,elph_k%nkpt) = itim
         
         elph_k%kptrank_t%invrank(symrankkpt) = elph_k%nkpt
         rankallk(elph_k%nkpt) = symrankkpt
       end if

     end do !end loop over irred k points
   end do !end loop over timrev
 end do !end loop over symmetry

 write(message,'(a,i0)')'mkfskgrid: after first evaluation, elph_k%nkpt= ', elph_k%nkpt
 call wrtout(std_out,message,"COLL")

 elph_k%wtkirr(:) = elph_k%wtkirr(:) / elph_k%nkpt

!copy the kpoints and full --> irred kpt map
!reorder the kpts to get rank increasing monotonically with a sort
!also reorder tmpkphon_full2irr
 ABI_ALLOCATE(elph_k%kpt,(3,elph_k%nkpt))
 ABI_ALLOCATE(elph_k%full2irr,(3,elph_k%nkpt))
 ABI_ALLOCATE(sortindexing,(elph_k%nkpt))

 do ikpt1=1,elph_k%nkpt
   sortindexing(ikpt1)=ikpt1
 end do
 call sort_int(elph_k%nkpt, rankallk, sortindexing)
 do ikpt1=1,elph_k%nkpt
   if (sortindexing(ikpt1) < 1 .or. sortindexing(ikpt1) > elph_k%nkpt) then
     MSG_BUG('sorted k ranks are out of bounds: 1 to nkpt')
   end if 
   elph_k%kpt(:,ikpt1) = tmpkpt(:,sortindexing(ikpt1))
   elph_k%full2irr(:,ikpt1) = tmpkphon_full2irr(:,sortindexing(ikpt1))
 end do

 ABI_DEALLOCATE(sortindexing)
 ABI_DEALLOCATE(rankallk)
 ABI_DEALLOCATE(tmpkphon_full2irr)
 ABI_DEALLOCATE(tmpkpt)
 call destroy_kptrank (elph_k%kptrank_t)


!make proper full rank arrays
 call mkkptrank (elph_k%kpt,elph_k%nkpt,elph_k%kptrank_t)


!find correspondence table between irred FS kpoints and a full one
 ABI_ALLOCATE(elph_k%irr2full,(elph_k%nkptirr))
 elph_k%irr2full(:) = 0

 do ikpt1=1,elph_k%nkptirr
   call get_rank_1kpt (elph_k%kptirr(:,ikpt1),symrankkpt,elph_k%kptrank_t)
   elph_k%irr2full(ikpt1) = elph_k%kptrank_t%invrank(symrankkpt)
 end do

!find correspondence table between FS kpoints under symmetry
 ABI_ALLOCATE(elph_k%full2full,(2,nsym,elph_k%nkpt))
 elph_k%full2full(:,:,:) = -999

 do ikpt1=1,elph_k%nkpt
!  generate symmetrics of kpt ikpt1
   do isym=1,nsym
     do itim=0,timrev
       timsign = one-two*itim
       kpt(:) = timsign*(symrec(:,1,isym)*elph_k%kpt(1,ikpt1) + &
&       symrec(:,2,isym)*elph_k%kpt(2,ikpt1) + &
&       symrec(:,3,isym)*elph_k%kpt(3,ikpt1))

!      which kpt is it among the full FS kpts
       call get_rank_1kpt (kpt,symrankkpt,elph_k%kptrank_t)
       ikpt2 = elph_k%kptrank_t%invrank(symrankkpt)
       new=1
       if (ikpt2 /= -1) then
         elph_k%full2full(itim+1,isym,ikpt2) = ikpt1
         new = 0
       end if

       if (new == 1) then
         write(std_out,*) ' mkfskgrid Error: FS kpt ',ikpt1,' has no symmetric under sym', isym,' with itim ',itim
         write(std_out,*) ' redkpt = ', redkpt
         write(std_out,*) ' symrankkpt,ikpt2 = ', symrankkpt,ikpt2
         MSG_ERROR("Fatal error, cannot continue")
       end if
     end do
   end do
 end do

!got nkpt, tmpkpt, kphon_full2irr, kphon_full2full, and wtkirr

end subroutine mkFSkgrid
!!***
