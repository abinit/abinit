!{\src2tex{textfont=tt}}
!!****f* ABINIT/ftgkk
!!
!! NAME
!! ftgkk
!!
!! FUNCTION
!! If qtor=1 (q->r):
!! Generates the Fourier transform of the recip space gkk matrices
!! to obtain the real space ones.
!! If qtor=0 (r->q):
!! Generates the Fourier transform of the real space gkk matrices
!! to obtain the reciprocal space ones.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2017 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gkqwrite = flag to write recip space matrix elements to disk
!! gkrwrite = flag to write real space matrix elements to disk
!! gprim(3,3)= Normalized coordinates in reciprocal space
!! ikpt_phon0 = starting kpt number for forward FT.
!! natom= Number of atoms in the unit cell
!! nkpt_phon= Number of kpoints used for the FS
!! ngkkband = number of bands kept in gkq and gkr matrix elements (=1 or nband)
!! nkpt_used= number of FS kpoints used, starting at ikpt_phon0
!! nqpt= Number of q points in the Brillouin zone
!!           if qtor=0 this number is read in the input file
!! nrpt= Number of R points in the Big Box
!! qtor= ( q to r : see above )
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!           These coordinates are normalized (=> * acell(3)!!)
!! qpt_full(3,nqpt)= Reduced coordinates of the q vectors in reciprocal space
!!           if qtor=0 these vectors are read in the input file
!! unit_gkk_rpt = fortran unit for writing real-space matrix elements
!! unitgkq = fortran unit for writing reciprocal-space matrix elements
!! wghatm(natom,natom,nrpt)
!!         = Weights associated to a pair of atoms and to a R vector
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/output
!! gkk_qpt(2,3*natom,nFSband,nFSband,nkpt_used,nqpt)
!!  = gkk matrices in recip space coming from the Derivative Data Base
!! gkk_rpt(2,3*natom,nFSband,nFSband,nkpt_phon,nqpt)
!!  = gkk matrices in real space stored in file unit_gkk_rpt
!!
!! PARENTS
!!      get_all_gkr,interpolate_gkk,test_ftgkk
!!
!! CHILDREN
!!
!! NOTES
!!   copied from ftiaf9.f
!!   recip to real space: real space is forced to disk file unit_gkk_rpt
!!                        recip space depends on gkqwrite and unitgkq
!!   real to recip space: real space is forced to disk file unit_gkk_rpt
!!                        recip space is necessarily in memory in gkk_qpt
!!
!!    real space elements are complex, but could be reduced, as (-r) = (+r)*
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine ftgkk (wghatm,gkk_qpt,gkk_rpt,gkqwrite,gkrwrite,gprim,ikpt_phon0,&
&                  natom,nkpt_phon,ngkkband,nkpt_used,nqpt,nrpt,nsppol,&
&                  qtor,rpt,qpt_full,unit_gkk_rpt,unitgkq)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ftgkk'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: gkqwrite,gkrwrite,ikpt_phon0,nkpt_phon,natom,ngkkband
 integer,intent(in) :: nkpt_used,nqpt,nrpt,nsppol,qtor,unit_gkk_rpt,unitgkq
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),qpt_full(3,nqpt)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 real(dp),intent(inout) :: gkk_qpt(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_used,nsppol,nqpt)
 real(dp),intent(inout) :: gkk_rpt(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_used,nsppol,nrpt)

!Local variables -------------------------
!scalars
 integer :: ikpt_phon,iatom,ib1,ieffkpt_phon,ip,iqpt,irpt,isppol
 integer :: jatom
 real(dp) :: im,kr,re
 character(len=500) :: message
!arrays
 real(dp) :: coskr(nqpt,nrpt),ftwght(2,3*natom*3*natom)
 real(dp) :: gkk_qpt_tmp(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_used,nsppol)
 real(dp) :: gkk_rpt_tmp(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_phon,nsppol)
 real(dp) :: kk(3),sinkr(nqpt,nrpt)

! *********************************************************************

!rewind (unit_gkk_rpt)

!prepare the phase factors
 do iqpt=1,nqpt
!  Calculation of the k coordinates in Normalized Reciprocal
!  coordinates
   kk(1)=   qpt_full(1,iqpt)*gprim(1,1)+&
&   qpt_full(2,iqpt)*gprim(1,2)+&
&   qpt_full(3,iqpt)*gprim(1,3)
   kk(2)=   qpt_full(1,iqpt)*gprim(2,1)+&
&   qpt_full(2,iqpt)*gprim(2,2)+&
&   qpt_full(3,iqpt)*gprim(2,3)
   kk(3)=   qpt_full(1,iqpt)*gprim(3,1)+&
&   qpt_full(2,iqpt)*gprim(3,2)+&
&   qpt_full(3,iqpt)*gprim(3,3)
   do irpt=1,nrpt
!    Product of k and r
     kr =        kk(1)*rpt(1,irpt)+&
&     kk(2)*rpt(2,irpt)+&
&     kk(3)*rpt(3,irpt)
     coskr(iqpt,irpt)=cos(two_pi*kr)
     sinkr(iqpt,irpt)=sin(two_pi*kr)
!    DEBUG
!    if (iqpt < 1000 .and. (irpt == 101 .or. irpt == 901)) then
!    write(std_out,*) iqpt,irpt,kk,rpt(:,irpt),coskr(iqpt,irpt), sinkr(iqpt,irpt)
!    end if
!    ENDDEBUG
   end do
 end do



!Recip to real space
 if (qtor==1) then
!  
   if (nkpt_used /= nkpt_phon) write(std_out,*) 'ftgkk: strange usage of nkpt_used for back FT!'
   do irpt=1,nrpt
!    DEBUG
!    write(std_out,*) ' ftgkk : G->R irpt = ',irpt,' / ',nrpt
!    ENDDEBUG
     gkk_rpt_tmp(:,:,:,:,:) = zero

     do iqpt=1,nqpt

!      write(std_out,*) iqpt

       if (gkqwrite == 0) then
         gkk_qpt_tmp(:,:,:,:,:) = gkk_qpt(:,:,:,:,:,iqpt)
       else
         do ikpt_phon=1, nkpt_phon
           read(unitgkq,REC=((iqpt-1)*nkpt_phon+ikpt_phon)) gkk_qpt_tmp(:,:,:,ikpt_phon,:)
         end do
       end if
!      Get the phase factor with normalization!
       re=coskr(iqpt,irpt)/nqpt
       im=sinkr(iqpt,irpt)/nqpt
       do isppol=1,nsppol
         do ikpt_phon=1,nkpt_used
!          DEBUG
!          write(std_out,*) ' ftgkk : G->R ikpt_phon = ',ikpt_phon,' / ',nkpt_used
!          ENDDEBUG
           do ip=1,3*natom*3*natom
!            Real and imaginary part of the real-space gkk matrices -> exp(-i k.r)
             do ib1=1,ngkkband*ngkkband
               gkk_rpt_tmp(1,ib1,ip,ikpt_phon,isppol) = gkk_rpt_tmp(1,ib1,ip,ikpt_phon,isppol)&
&               +re*gkk_qpt_tmp(1,ib1,ip,ikpt_phon,isppol) &
&               +im*gkk_qpt_tmp(2,ib1,ip,ikpt_phon,isppol)
               gkk_rpt_tmp(2,ib1,ip,ikpt_phon,isppol) = gkk_rpt_tmp(2,ib1,ip,ikpt_phon,isppol)&
&               +re*gkk_qpt_tmp(2,ib1,ip,ikpt_phon,isppol) &
&               -im*gkk_qpt_tmp(1,ib1,ip,ikpt_phon,isppol)
             end do
           end do
         end do
       end do
     end do
     if (gkrwrite == 0) then
       gkk_rpt(:,:,:,:,:,irpt) = gkk_rpt_tmp(:,:,:,:,:)
     else
       write (unit_gkk_rpt,REC=irpt) gkk_rpt_tmp
     end if
   end do

!  Real space to recip space
 else if (qtor==0) then

!  write(std_out,*) 'ftgkk : shape(gkk_qpt) = ', shape(gkk_qpt)
   gkk_qpt(:,:,:,:,:,:)=zero

!  rewind (unit_gkk_rpt)
   do irpt=1,nrpt
     if (gkrwrite == 0) then
       gkk_rpt_tmp(:,:,:,:,:) = gkk_rpt(:,:,:,:,:,irpt)
     else
       read(unit_gkk_rpt,REC=irpt) gkk_rpt_tmp
     end if


     do iqpt=1,nqpt

!      Avoid recalculating weights nkpt_used*9 times
       do iatom=1,natom
         do jatom=1,natom
           ip = 3*((iatom-1)*natom+jatom-1)
!          copy same weight for all 3 directions
           ftwght(1,ip+1:ip+3)=coskr(iqpt,irpt)*wghatm(iatom,jatom,irpt)
           ftwght(2,ip+1:ip+3)=sinkr(iqpt,irpt)*wghatm(iatom,jatom,irpt)
         end do
       end do



       do ip=1,3*natom*3*natom
!        Get phase factor
         re = ftwght(1,ip)
         im = ftwght(2,ip)

         do isppol=1,nsppol
           do ikpt_phon=1,nkpt_used


!            DEBUG
!            write(std_out,*) ' ftgkk : R->G ikpt_phon = ',ikpt_phon,' / ',nkpt_used
!            ENDDEBUG
!            effective FS kpt in real space array is ikpt_phon+ikpt_phon0-1 to allow for offset
             ieffkpt_phon = ikpt_phon+ikpt_phon0-1
!            write(std_out,*) 'ftgkk :ikpt_phon,iqpt,ieffkpt_phon ', ikpt_phon,iqpt,ieffkpt_phon

             do ib1=1,ngkkband*ngkkband
!              Real and imaginary part of the gamma matrices
               gkk_qpt(1,ib1,ip,ikpt_phon,isppol,iqpt)=&
&               gkk_qpt(1,ib1,ip,ikpt_phon,isppol,iqpt)&
&               +re*gkk_rpt_tmp(1,ib1,ip,ieffkpt_phon,isppol)&
&               -im*gkk_rpt_tmp(2,ib1,ip,ieffkpt_phon,isppol)
!              !DEBUG
               gkk_qpt(2,ib1,ip,ikpt_phon,isppol,iqpt)=&
&               gkk_qpt(2,ib1,ip,ikpt_phon,isppol,iqpt)&
&               +im*gkk_rpt_tmp(1,ib1,ip,ieffkpt_phon,isppol)&
&               +re*gkk_rpt_tmp(2,ib1,ip,ieffkpt_phon,isppol)
!              !ENDDEBUG

!              if (iqpt < 100 .and. irpt < 100 .and. &
!              &   tmpgkkrim(irpt)**2+tmpgkkrre(irpt)**2 > tol6) then
!              write(std_out,'(2I4,2E16.8,x,2E16.8)') &
!              &   iqpt,irpt,re,im,tmpgkkrre(irpt),tmpgkkrim(irpt)
!              end if

             end do
           end do
!          end ikpt_phon
         end do
!        end isppol
!        write(std_out,'(a)') ' ftgkk :gkk_qpt :'
!        write(std_out,'(4E16.5)') gkk_qpt(:,1,1,,ikpt_phon,1:nqpt)
       end do
!      end ip
     end do
!    end iqpt
   end do
!  end irpt


!  There is no other space to Fourier transform from ??
 else
   write(message,'(a,a,a,i0,a)' )&
&   'The only allowed values for qtor are 0 or 1, while',ch10,&
&   'qtor=',qtor,' has been required.'
   MSG_BUG(message)
 end if

end subroutine ftgkk
!!***
