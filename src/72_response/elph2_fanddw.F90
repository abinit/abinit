!{\src2tex{textfont=tt}}
!!****f* ABINIT/elph2_fanddw
!! NAME
!! elph2_fanddw
!!
!! FUNCTION
!! This routine calculates the zero-point motion corrections
!! due to the Fan term or to the DDW term..
!!
!! COPYRIGHT
!! Copyright (C) 2011-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors .
!!
!! INPUTS
!!  dim_eig2nkq=1 if eig2nkq is to be computed
!!  displ(2*3*natom*3*natom)=the displacements of atoms in cartesian coordinates.
!!  eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom*dim_eig2nkq)=one half second derivatives of the electronic eigenvalues
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  mband= maximum number of bands
!!  natom= number of atoms in the unit cell
!!  nkpt= number of k-points
!!  nsppol= 1 for unpolarized, 2 for spin-polarized
!!  option 1 for Fan term, 2 for DDW term
!!  phfrq(3*natom)=phonon frequencies
!!  (prtvol > 4) if the mode decomposition is to be printed
!!
!! OUTPUT
!!  eigen_corr(mband*nkpt*nsppol)= T=0 correction to the electronic eigenvalues, due to the Fan term.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine elph2_fanddw(dim_eig2nkq,displ,eig2nkq,eigen_corr,gprimd,mband,natom,nkpt,nsppol,option,phfrq,prtvol)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elph2_fanddw'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dim_eig2nkq,mband,natom,nkpt,nsppol,option,prtvol

!arrays
 real(dp) :: gprimd(3,3)
 real(dp),intent(in) :: displ(2*3*natom*3*natom)
 real(dp),intent(in) :: eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom*dim_eig2nkq)
 real(dp),intent(in) :: phfrq(3*natom)
 real(dp),intent(out) :: eigen_corr(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: neigs_per_line=6
 integer :: iatom1,iatom2,idir1,idir2,iband,ikpt,imode,index,isppol, imin, ii
 real(dp) :: d_at1_dir1_re,d_at1_dir1_im
 real(dp) :: d_at1_dir2_re,d_at1_dir2_im
 real(dp) :: d_at2_dir1_re,d_at2_dir1_im
 real(dp) :: d_at2_dir2_re,d_at2_dir2_im
 real(dp) :: e2_im,e2_re
 real(dp), allocatable :: eigen_corr_mode(:)
 character(len=500) :: message
 character(len=20) :: eig_format, line_format
!arrays
 real(dp) :: displ2cart(2,3,3),displ2red(2,3,3),tmp_displ2(2,3,3)

! *********************************************************************

!DEBUG
!write(std_out,*)' elph2_fanddw : enter '
!write(std_out,*)' option=',option
!ENDDEBUG

 if(option/=1 .and. option/=2)then
   write(message,'(a,i0)')' The argument option should be 1 or 2, while it is found that option=',option
   MSG_BUG(message)
 end if

 !printing options
 eig_format='f16.8'
 write(line_format,'(a,i1,a6,a)') '(',neigs_per_line,eig_format,')'

 if (prtvol > 4) then
   write(message,'(a,a)')ch10,' ================================================================================'
   call wrtout(ab_out_default,message,'COLL')
   if (option==1) then
     write(message,'(a)') ' ---- Begin Fan contributions to eigenvalues renormalization by mode ----'
     call wrtout(ab_out_default,message,'COLL')
   else if (option==2) then
     write(message,'(a)') ' ---- Begin DDW contributions to eigenvalues renormalization by mode ----'
     call wrtout(ab_out_default,message,'COLL')
   end if
 end if

 ABI_ALLOCATE(eigen_corr_mode,(mband*nkpt*nsppol))
 
 eigen_corr(:)=zero
 do imode=1,3*natom
   eigen_corr_mode(:)=zero

!  DEBUG
!  write(std_out,*)' Contribution of mode ',imode,' with frequency=',phfrq(imode),' and displacements :'
!  write(std_out,'(2f14.7)' ) displ(1+2*3*natom*(imode-1):2*3*natom*imode)
!  ENDDEBUG

   if (phfrq(imode)>tol6) then
     do iatom1=1,natom
       do iatom2=1,natom
!        DEBUG
!        write(std_out,*)' iatom1,iatom2=',iatom1,iatom2
!        ENDDEBUG

         do idir1=1,3
           do idir2=1,3
!            Compute the mean cartesian displacements
             d_at1_dir1_re=displ(1 + 2*(idir1-1 +3*(iatom1-1 +natom*(imode-1))))
             d_at1_dir1_im=displ(2 + 2*(idir1-1 +3*(iatom1-1 +natom*(imode-1))))
             d_at2_dir2_re=displ(1 + 2*(idir2-1 +3*(iatom2-1 +natom*(imode-1))))
             d_at2_dir2_im=displ(2 + 2*(idir2-1 +3*(iatom2-1 +natom*(imode-1))))

!            DEBUG
!            write(std_out,*)' idir1,idir2=',iatom1,iatom2,idir1,idir2
!            write(std_out,'(a,4f12.5)' )' d_at1_dir1 re,d_at2_dir2 re=',d_at1_dir1_re,d_at2_dir2_re
!            ENDDEBUG

             if(option==1)then
!              Compute the mean displacement correlation at T=0. 
!              Consistent with Eqs.(7) and (8) of PRB51, 8610 (1995), specialized for the contribution of one q point.
!              but generalized to two different atoms. Note that the complex conjugate is taken on the second direction.
               displ2cart(1,idir1,idir2)=(d_at1_dir1_re*d_at2_dir2_re+ &
&               d_at1_dir1_im*d_at2_dir2_im )/(two*phfrq(imode))
               displ2cart(2,idir1,idir2)=(d_at1_dir1_im*d_at2_dir2_re- &
&               d_at1_dir1_re*d_at2_dir2_im )/(two*phfrq(imode))
             else if(option==2)then
!              Compute the mean square displacement correlation of each atom at T=0, and take mean over iatom1 and iatom2. 
!              See Eqs.(7) and (8) of PRB51, 8610 (1995), specialized for the contribution of one q point.
!              Note that the complex conjugate is taken on the second direction.
!              Also, note the overall negative sign, to make it opposite to the Fan term.
               d_at1_dir2_re=displ(1 + 2*(idir2-1 +3*(iatom1-1 +natom*(imode-1))))
               d_at1_dir2_im=displ(2 + 2*(idir2-1 +3*(iatom1-1 +natom*(imode-1))))
               d_at2_dir1_re=displ(1 + 2*(idir1-1 +3*(iatom2-1 +natom*(imode-1))))
               d_at2_dir1_im=displ(2 + 2*(idir1-1 +3*(iatom2-1 +natom*(imode-1))))
               displ2cart(1,idir1,idir2)=-(d_at1_dir1_re*d_at1_dir2_re+ &
&               d_at1_dir1_im*d_at1_dir2_im+ &
&               d_at2_dir1_re*d_at2_dir2_re+ &
&               d_at2_dir1_im*d_at2_dir2_im )/(four*phfrq(imode))
               displ2cart(2,idir1,idir2)=-(d_at1_dir1_im*d_at1_dir2_re- &
&               d_at1_dir1_re*d_at1_dir2_im+ &
&               d_at2_dir1_im*d_at2_dir2_re- &
&               d_at2_dir1_re*d_at2_dir2_im )/(four*phfrq(imode))
             end if
           end do           
         end do           
!        Switch to reduced coordinates in two steps
         tmp_displ2(:,:,:)=zero
         do idir1=1,3
           do idir2=1,3
             tmp_displ2(:,:,idir1)=tmp_displ2(:,:,idir1)+displ2cart(:,:,idir2)*gprimd(idir2,idir1)
           end do
         end do
         displ2red(:,:,:)=zero
         do idir1=1,3
           do idir2=1,3
             displ2red(:,idir1,:)=displ2red(:,idir1,:)+tmp_displ2(:,idir2,:)*gprimd(idir2,idir1)
           end do
         end do
!        Compute the T=0 shift due to this q point
         do idir1=1,3
           do idir2=1,3
             do ikpt=1,nkpt
               do isppol=1,nsppol
                 do iband=1,mband
                   index=iband+mband*(isppol-1 + nsppol*(ikpt-1))
                   e2_re=eig2nkq(1,iband+mband*(isppol-1),ikpt,idir1,iatom1,idir2,iatom2)
                   e2_im=eig2nkq(2,iband+mband*(isppol-1),ikpt,idir1,iatom1,idir2,iatom2)
                   eigen_corr(index)=eigen_corr(index)+&
&                   e2_re*displ2red(1,idir1,idir2)-e2_im*displ2red(2,idir1,idir2)
                   eigen_corr_mode(index)=eigen_corr_mode(index)+&
&                   e2_re*displ2red(1,idir1,idir2)-e2_im*displ2red(2,idir1,idir2)
                 end do  ! band
               end do  ! spin
             end do  ! kpt
           end do  ! dir2
         end do  ! dir1
       end do  ! atom2
     end do  ! atom1
   end if

   if (prtvol > 4) then
     ! Print the corrections by mode
     write(message,'(a,i1)') ' imode= ',imode
     call wrtout(ab_out_default,message,'COLL')

     do ikpt=1,nkpt
       do isppol=1,nsppol
         write(message,'(a,i4,a,i1)')' ikpt= ',ikpt,' ispin= ',isppol
         call wrtout(ab_out_default,message,'COLL')

         imin = mband * (isppol-1 + nsppol*(ikpt-1))
         do ii=0, (mband-1)/neigs_per_line
           write(message, line_format) (eigen_corr_mode(iband+imin), &
&           iband = 1 + ii * neigs_per_line, min(mband, (ii+1)*neigs_per_line))
           call wrtout(ab_out_default,message,'COLL')
         end do
       end do
     end do
   end if

 end do  ! mode

 if (prtvol > 4) then
   if (option==1) then
     write(message,'(a)') ' ---- End Fan contribution to eigenvalues renormalization by mode ----'
     call wrtout(ab_out_default,message,'COLL')
   else if (option==2) then
     write(message,'(a)') ' ---- End DDW contribution to eigenvalues renormalization by mode ----'
     call wrtout(ab_out_default,message,'COLL')
   end if
   write(message,'(a,a)')' ================================================================================', ch10
   call wrtout(ab_out_default,message,'COLL')
 end if

 ABI_DEALLOCATE(eigen_corr_mode)

!DEBUG
!write(std_out,*)' elph2_fanddw : exit'
!write(std_out,*)' eigen_corr(1)=',eigen_corr(1)
!ENDDEBUG

end subroutine elph2_fanddw
!!***

