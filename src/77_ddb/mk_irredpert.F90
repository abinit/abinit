!{\src2tex{textfont=tt}}
!!****f* ABINIT/mk_irredpert
!!
!! NAME
!! mk_irredpert
!!
!! FUNCTION
!! This routine finds the symop combinations needed to
!!   calculate other perturbations based on the present one.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2017 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  indsym = indirect indexing array for symmetries
!!  iqptfull = qpoint index in full grid
!!  natom = number of atoms
!!  nbranch = 3*natom = number of phonon branches
!!  nqpt = total number of qpoints
!!  nsym = number of symops
!!  qpt = qpoint in reduced coordinates
!!  qtimrev = 1 or 0, use time reversal symmetry (for Gamma only) or do not
!!  symq = 1 if symmetry preserves present qpoint. From littlegroup_q
!!  symrel = symmetry operations in recip space (or real)
!!
!! OUTPUT
!!  irredpert(7,nbranch,nqpt) = indices for reconstructing perturbations
!!    explained in NOTES. If irredpert < 0 then there is no information
!!    for reconstructing the element, and it should be read in (it is irreducible).
!!
!! NOTES
!!
!!       irredpert (1,,,) = iatsy1
!!       irredpert (2,,,) = iatsy2
!!       irredpert (3,,,) = isym
!!       irredpert (4,,,) = itim
!!       irredpert (5,,,) = factor for idir2 = 1
!!       irredpert (6,,,) = factor for idir2 = 2
!!       irredpert (7,,,) = factor for idir2 = 3
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mk_irredpert(indsym,iqptfull,irredpert,&
&    natom,nbranch,nqpt,nsym,qpt,qtimrev,symq,symrel)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mk_irredpert'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptfull,natom,nbranch,nqpt,nsym,qtimrev
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),qpt(3),symq(4,2,nsym)
 integer,intent(in) :: symrel(3,3,nsym)
 integer,intent(out) :: irredpert(7,nbranch,nbranch,nqpt)

!Local variables -------------------------
!scalars
 integer :: found,iatom1,iatom2,iatsy1,iatsy2,idir1,idir2,idisy1,idisy2
 integer :: ipert1,ipert2,ipertsy1,ipertsy2,isign,isym,ithree,itim
 integer :: noccur,quit,quit1
 real(dp) :: arg1,arg2,im,re
!arrays
 integer :: indsyminv(4,nsym,natom),sym1(3,3),sym2(3,3)

! *************************************************************************

!-------------------------------------------------------
!for present symmetries, find equivalent perturbations
!-------------------------------------------------------
!irredpert (1,,,) = ipert1
!irredpert (2,,,) = ipert2
!irredpert (3,,,) = isym
!irredpert (4,,,) = itim
!irredpert (5,,,) = factor for idir2 = 1
!irredpert (6,,,) = factor for idir2 = 2
!irredpert (7,,,) = factor for idir2 = 3
!

!Invert the indsym array
 do iatom1=1,natom
   do isym=1,nsym
     indsyminv(1,nsym,indsym(4,isym,iatom1)) = -indsym(1,isym,iatom1)
     indsyminv(2,nsym,indsym(4,isym,iatom1)) = -indsym(2,isym,iatom1)
     indsyminv(3,nsym,indsym(4,isym,iatom1)) = -indsym(3,isym,iatom1)
     indsyminv(4,isym,indsym(4,isym,iatom1)) = iatom1
   end do
 end do

!*********************************************************************

!DEBUG
!write(std_out,*)' mk_irredpert : enter mk_irredpert '
!write(std_out,*)' mk_irredpert : qtimrev=',qtimrev
!ENDDEBUG


!Big Big Loop : symmetrize three times, because
!of some cases in which one element is not yet available
!at the first pass, and even at the second one !
 do ithree=1,3

!  Big loop on all elements
   do iatom1=1,natom+2
     do idir1=1,3
       ipert1=(iatom1-1)*3+idir1
       do iatom2=1,natom+2
         do idir2=1,3
           ipert2=(iatom2-1)*3+idir2

!          Will try to eliminate element (idir1,iatom1,idir2,iatom2)
!          so this element should not have been eliminated yet...
           if(irredpert (1,ipert1,ipert2,iqptfull) > 0) cycle
!          write(std_out,*) 'trying to eliminate ',idir1,iatom1,idir2,iatom2

!          Loop on all symmetries, including time-reversal
           quit1=0
           do isym=1,nsym
             do itim=1,2
               isign=3-2*itim

               if(symq(4,itim,isym)==0) cycle
               found=1

!              Here select the symmetric of iatom1
               if(iatom1<=natom)then
                 iatsy1=indsyminv(4,isym,iatom1)
                 sym1(:,:)=symrel(:,:,isym)
               else
                 found=0
               end if

!              Here select the symmetric of iatom2
               if(iatom2<=natom)then
                 iatsy2=indsyminv(4,isym,iatom2)
!                try this: invert all relations symrec -> symrel
                 sym2(:,:)=symrel(:,:,isym)
               else
                 found=0
               end if

!              Now that a symmetric perturbation has been obtained,
!              including the expression of the symmetry matrix, see
!              if the symmetric values are available
!              first if found==1
               if( found==1 ) then

                 noccur=0
                 quit=0
                 do idisy1=1,3
                   ipertsy1 = (iatsy1-1)*3+idisy1
                   do idisy2=1,3
                     ipertsy2 = (iatsy2-1)*3+idisy2
!                    
!                    NOTE : May be transpose on the next line
!                    
                     if(sym1(idir1,idisy1)/=0 .and. sym2(idir2,idisy2)/=0 )then
                       if(irredpert (1,ipertsy1,ipertsy2,iqptfull) < 0)then
!                        nothing, ipertsy1,ipertsy2 has not been eliminated by symmetry,
!                        just leave found == 1


                       else
!                        Not found: ipert1,   ipert2   can not be reconstituted from
!                        ipertsy1, ipertsy2 using isym
!                        write(std_out,*) 'Not found: ',ipert1,ipert2,&
!                        &                 ' can not be reconstituted from ', ipertsy1,ipertsy2,&
!                        &                 ' using ', isym
                         found=0
                         quit=1
                         exit
                       end if

!                      Here, in case the symmetric of the element
!                      is the element, or the symmetric with
!                      respect to permutation of perturbations
!                      (some more conditions on the time-reversal
!                      symmetry must be fulfilled although)
                       if(  idisy1==idir1 .and. iatsy1==iatom1&
&                       .and. idisy2==idir2 .and. iatsy2==iatom2&
&                       .and.(isign==1 .or. qtimrev==1 &
&                       .or. (idir1==idir2 .and. iatom1==iatom2)))&
&                       then
                         noccur=noccur+sym1(idir1,idisy1)*sym2(idir2,idisy2)
                       else if(  idisy1==idir2 .and. iatsy1==iatom2&
&                         .and. idisy2==idir1 .and. iatsy2==iatom1&
&                         .and.(isign==-1 .or. qtimrev==1&
&                         .or. (idir1==idir2 .and. iatom1==iatom2)))&
&                         then
                         noccur=noccur+sym1(idir1,idisy1)*sym2(idir2,idisy2)
                       end if

                     end if
                   end do
                   if(quit==1)exit
                 end do
               end if
!              End first if found==1

!              SHOULD THIS BE CHECKED IN THE ELPHON CASE TOO?
!              second if found == 1
               if(found==1)then
!                In case of phonons, need to take into account the
!                time-reversal symmetry, and the shift back to the unit cell
!                
                 if(ipert1<=natom .and. ipert2<=natom)then
!                  2) Shift the atoms back to the unit cell.
                   arg1=two_pi*( qpt(1)*indsyminv(1,isym,iatom1)&
&                   +qpt(2)*indsyminv(2,isym,iatom1)&
&                   +qpt(3)*indsyminv(3,isym,iatom1) )
                   arg2=two_pi*( qpt(1)*indsyminv(1,isym,iatom2)&
&                   +qpt(2)*indsyminv(2,isym,iatom2)&
&                   +qpt(3)*indsyminv(3,isym,iatom2) )
                   re=cos(arg1)*cos(arg2)+sin(arg1)*sin(arg2)
!                  XG010117 Must use isign
                   im=isign*(cos(arg2)*sin(arg1)-cos(arg1)*sin(arg2))
                 else
                   re=1.0_dp
                   im=0.0_dp
                 end if

!                Final check, could still fail if the
!                element was its own symmetric
                 if( abs(1.0_dp-re*noccur) < 1.0d-6&
&                 .and.  abs(im*noccur)  < 1.0d-6 )then
                   found=0

!                  DEBUG
!                  write(std_out,*)' element is its own symmetric ...'
!                  ENDDEBUG

                 end if
!                End if element was its own symmetric

               end if
!              End second if found == 1

!              third if found == 1
               if(found==1)then

!                DEBUG
!                write(std_out,*)' all found !  isym, isign= ',isym,isign
!                write(std_out,'(9i4)' )((sym1(ii,jj),ii=1,3),jj=1,3)
!                write(std_out,'(9i4)' )((sym2(ii,jj),ii=1,3),jj=1,3)
!                ENDDEBUG

!                The element has been constructed !
!                ipertsy1,ipertsy2 do not correspond to anything: all
!                pertsy's were scanned above, but all the needed ones
!                are present.
                 irredpert (1,ipert1,ipert2,iqptfull) = iatsy1
                 irredpert (2,ipert1,ipert2,iqptfull) = iatsy2
                 irredpert (3,ipert1,ipert2,iqptfull) = isym
                 irredpert (4,ipert1,ipert2,iqptfull) = itim
                 irredpert (5,ipert1,ipert2,iqptfull) = 0
                 irredpert (6,ipert1,ipert2,iqptfull) = 0
                 irredpert (7,ipert1,ipert2,iqptfull) = 0

!                Exit loop on symmetry operations
                 quit1=1
                 exit

               end if
!              end third if found == 1

!              End loop on all symmetries + time-reversal
             end do
!            end isym do
             if(quit1==1)exit
           end do
!          end itim do

!          End big loop on all elements

         end do
       end do
!      end ipert2 do
     end do
   end do
!  end ipert1 do

!  End Big Big Loop
 end do

!DEBUG
!write(std_out,*) ' mk_irredpert : irredpert for qpt ',iqptfull,' = '
!do ibranch=1,nbranch
!do jbranch=1,nbranch
!write(std_out,'(2i4,2x,7i6)') ibranch,jbranch,irredpert(:,ibranch,jbranch,iqptfull)
!end do
!end do
!ENDDEBUG

end subroutine mk_irredpert
!!***
