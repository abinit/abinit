!{\src2tex{textfont=tt}}
!!****f* ABINIT/smpbz
!!
!! NAME
!! smpbz
!!
!! FUNCTION
!! Generate a set of special k (or q) points which samples in a homogeneous way
!! the entire Brillouin zone of a simple lattice, face-centered cubic,
!! body-centered lattice and hexagonal lattice.
!! If kptrlatt is diagonal, the algorithm used here reduces to the usual
!! Monkhorst-Pack set of k points.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (JCC,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  brav = 1 -> simple lattice; 2 -> face-centered cubic;
!!   3 -> body-centered lattice; 4 -> hexagonal lattice (D6h)
!!  iout = unit number for output
!!  kptrlatt(3,3)=integer coordinates of the primitive vectors of the
!!   lattice reciprocal to the k point lattice to be generated here
!!   If diagonal, the three values are the Monkhorst-Pack usual values, in case of simple cubic.
!!  mkpt = maximum number of k points
!!  nshiftk= number of shift vectors in the repeated cell
!!  option= Flag defining what will be printed of iout: 0 for k points, anything else for q points.
!!    Also, for q points, if the Gamma point is present, place it first in the list.
!!  shiftk(3,nshiftk) = vectors that will be used to determine the shifts from (0. 0. 0.).
!!
!! OUTPUT
!!  nkpt = number of k points
!!  spkpt(3,mkpt) = the nkpt first values contain the special k points
!!   obtained by the Monkhorst & Pack method, in reduced coordinates.
!!   These vectors have to be multiplied by the reciprocal basis vectors
!!   gprimd(3,3) (in cartesian coordinates) to obtain the special k points
!!   set in cartesian coordinates.
!!
!! NOTES
!!  also allows for more than one vector in repeated cell.
!!  this routine should be rewritten, to use the Wigner-Seitz cell,
!!  and thus unify the different treatments.
!!  References :
!!  H.J. Monkhorst and J.D. Pack, Phys. Rev. B 13, 5188 (1976)
!!  J.D. Pack and H.J. Monkhorst, Phys. Rev. B 16, 1748 (1977)
!!  A.H. MacDonald, Phys. Rev. B 18, 5897 (1978)
!!  R.A. Evarestov and V.P. Smirnov, Phys. Stat. Sol. (b) 119, 9 (1983)
!!
!! PARENTS
!!      ep_setupqpt,getkgrid,harmonic_thermo,initberry,invars2,m_fock,m_fstab
!!      m_ifc
!!
!! CHILDREN
!!      matr3inv,wrap2_pmhalf,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine smpbz(brav,iout,kptrlatt,mkpt,nkpt,nshiftk,option,shiftk,spkpt)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_numeric_tools,   only : wrap2_pmhalf

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'smpbz'
 use interfaces_14_hidewrite
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: brav,iout,mkpt,nshiftk,option
 integer,intent(out) :: nkpt
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk)
 real(dp),intent(out) :: spkpt(3,mkpt)

!Local variables -------------------------
!scalars
 integer :: ii,ikshft,jj,kk,nkpout,nkptlatt,nn
 real(dp) :: shift
 character(len=500) :: message
!arrays
 integer :: boundmax(3),boundmin(3),ngkpt(3)
 real(dp) :: k1(3),k2(3),kcar(3),klatt(3,3),ktest(3),rlatt(3,3)

! *********************************************************************

!DEBUG
!write(std_out,*)' smpbz : kptrlatt(:)=',kptrlatt(1,1),kptrlatt(2,2)
!ENDDEBUG

 if(option/=0)then
   call wrtout(iout,'       Homogeneous q point set in the B.Z.  ','COLL')
 end if

 if(brav/=1)then
!  Only generate Monkhorst-Pack lattices
   if(kptrlatt(1,2)/=0 .or. kptrlatt(2,1)/=0 .or. &
&   kptrlatt(1,3)/=0 .or. kptrlatt(3,1)/=0 .or. &
&   kptrlatt(2,3)/=0 .or. kptrlatt(3,2)/=0     ) then
     write(message, '(2a,a,3i4,a,a,3i4,a,a,3i4)' )&
&     '  When brav/=1, kptrlatt must be diagonal, while it is',ch10,&
&     '  kptrlatt(:,1)=',kptrlatt(:,1),ch10,&
&     '  kptrlatt(:,2)=',kptrlatt(:,2),ch10,&
&     '  kptrlatt(:,3)=',kptrlatt(:,3)
     MSG_BUG(message)
   end if

   ngkpt(1)=kptrlatt(1,1)
   ngkpt(2)=kptrlatt(2,2)
   ngkpt(3)=kptrlatt(3,3)
!  
   if( (ngkpt(1)<=0.or.ngkpt(2)<=0.or.ngkpt(3)<=0) .and. (ngkpt(1)/=0.or.ngkpt(2)/=0.or.ngkpt(3)/=0) ) then
     write(message, '(5a,i4,a,a,i4,a,a,i4,a,a)' )&
&     '  All ngkpt (or ngqpt) must be strictly positive',ch10,&
&     '  or all ngk(q)pt must be zero (for Gamma sampling), but :',ch10,&
&     ' ngk(q)pt(1) = ',ngkpt(1),ch10,&
&     ' ngk(q)pt(2) = ',ngkpt(2),ch10,&
&     ' ngk(q)pt(3) = ',ngkpt(3),ch10,&
&     ' Action : correct ngkpt or ngqpt in the input file.'
     MSG_BUG(message)
   end if
 end if

!*********************************************************************

 if(brav==1)then

!  Compute the number of k points in the G-space unit cell
!  (will be multiplied by nshiftk later.
   nkptlatt=kptrlatt(1,1)*kptrlatt(2,2)*kptrlatt(3,3) &
&   +kptrlatt(1,2)*kptrlatt(2,3)*kptrlatt(3,1) &
&   +kptrlatt(1,3)*kptrlatt(2,1)*kptrlatt(3,2) &
&   -kptrlatt(1,2)*kptrlatt(2,1)*kptrlatt(3,3) &
&   -kptrlatt(1,3)*kptrlatt(2,2)*kptrlatt(3,1) &
&   -kptrlatt(1,1)*kptrlatt(2,3)*kptrlatt(3,2)

!  Simple Lattice
   call wrtout(std_out,'       Simple Lattice Grid ','COLL')
   if (mkpt<nkptlatt*nshiftk) then
     write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&     '  The value of mkpt is not large enough. It should be',ch10,&
&     '  at least',nkptlatt*nshiftk,',',ch10,&
&     '  Action : set mkpt to that value in the main routine,',ch10,&
&     '  and recompile the code.'
     MSG_BUG(message)
   end if

!  Build primitive vectors of the k lattice
   rlatt(:,:)=kptrlatt(:,:)
   call matr3inv(rlatt,klatt)

!  Now, klatt contains the three primitive vectors of the k lattice,
!  in reduced coordinates. One builds all k vectors that
!  are contained in the first Brillouin zone, with coordinates
!  in the interval [0,1[ . First generate boundaries of a big box.

   do jj=1,3
!    To accomodate the shifts, boundmin and boundmax don't start from 0

     if ( maxval(shiftk(jj,:)) .eq. zero ) then
       boundmin(jj)=-1
     else
       boundmin(jj)=-ceiling(max(maxval(shiftk(jj,:)),tol14)) !-1
     end if

     if ( minval(shiftk(jj,:)) .eq. zero ) then
       boundmax(jj)=0
     else
       boundmax(jj)=-floor(min(minval(shiftk(jj,:)),tol14)) !0
     end if

     do ii=1,3
       if(kptrlatt(ii,jj)<0)boundmin(jj)=boundmin(jj)+kptrlatt(ii,jj)
       if(kptrlatt(ii,jj)>0)boundmax(jj)=boundmax(jj)+kptrlatt(ii,jj)
     end do
   end do

   nn=1
   do kk=boundmin(3),boundmax(3)
     do jj=boundmin(2),boundmax(2)
       do ii=boundmin(1),boundmax(1)
         do ikshft=1,nshiftk
!          Coordinates of the trial k point with respect to the k primitive lattice
           k1(1)=ii+shiftk(1,ikshft)
           k1(2)=jj+shiftk(2,ikshft)
           k1(3)=kk+shiftk(3,ikshft)
!          Reduced coordinates of the trial k point
           k2(:)=k1(1)*klatt(:,1)+k1(2)*klatt(:,2)+k1(3)*klatt(:,3)
!          Eliminate the point if outside [0,1[
           if(k2(1)<-tol10)cycle ; if(k2(1)>one-tol10)cycle
           if(k2(2)<-tol10)cycle ; if(k2(2)>one-tol10)cycle
           if(k2(3)<-tol10)cycle ; if(k2(3)>one-tol10)cycle
!          Wrap the trial values in the interval ]-1/2,1/2] .
           call wrap2_pmhalf(k2(1),k1(1),shift)
           call wrap2_pmhalf(k2(2),k1(2),shift)
           call wrap2_pmhalf(k2(3),k1(3),shift)
           spkpt(:,nn)=k1(:)
           nn=nn+1
         end do
       end do
     end do
   end do
   nkpt=nn-1

   if(nkpt/=nkptlatt*nshiftk)then
     write(message, '(a,i8,a,a,a,i8,a)' )&
&     '  The number of k points ',nkpt,'  is not equal to',ch10,&
&     '  nkptlatt*nshiftk which is',nkptlatt*nshiftk,'.'
     MSG_BUG(message)
   end if

 else if(brav==2)then

!  Face-Centered Lattice
   call wrtout(std_out,'       Face-Centered Lattice Grid ','COLL')
   if (mkpt<ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk/2) then
     write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&     '  The value of mkpt is not large enough. It should be',ch10,&
&     '  at least',(ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/2,',',ch10,&
&     '  Action : set mkpt to that value in the main routine,',ch10,&
&     '  and recompile the code.'
     MSG_BUG(message)
   end if
   nn=1
   if (ngkpt(1)/=ngkpt(2).or.ngkpt(1)/=ngkpt(3)) then
     write(message, '(4a,3(a,i6,a),a)' )&
&     '  For face-centered lattices, the numbers ngqpt(1:3)',ch10,&
&     '  must be equal, while they are :',ch10,&
&     '  ngqpt(1) = ',ngkpt(1),ch10,&
&     '  ngqpt(2) = ',ngkpt(2),ch10,&
&     '  ngqpt(3) = ',ngkpt(3),ch10,&
&     '  Action : modify ngqpt(1:3) in the input file.'
     MSG_BUG(message)
   end if
   if ((ngkpt(1)*nshiftk)/=(((ngkpt(1)*nshiftk)/2)*2)) then
     write(message, '(4a,3(a,i6,a),a)' )&
&     '  For face-centered lattices, the numbers ngqpt(1:3)*nshiftk',ch10,&
&     '  must be even, while they are :',ch10,&
&     '  ngqpt(1)*nshiftk = ',ngkpt(1)*nshiftk,ch10,&
&     '  ngqpt(2)*nshiftk = ',ngkpt(2)*nshiftk,ch10,&
&     '  ngqpt(3)*nshiftk = ',ngkpt(3)*nshiftk,ch10,&
&     '  Action : modify ngqpt(1:3)*nshiftk in the input file.'
     MSG_ERROR(message)
   end if
   if (ngkpt(1)==0.or.ngkpt(2)==0.or.ngkpt(3)==0) then
     spkpt(1,1)=0.0_dp
     spkpt(2,1)=0.0_dp
     spkpt(3,1)=0.0_dp
     nkpt=1
   else
     do kk=1,ngkpt(3)
       do jj=1,ngkpt(2)
         do ii=1,ngkpt(1)
           do ikshft=1,nshiftk
             k1(1)=(ii-1+shiftk(1,ikshft))/ngkpt(1)
             k1(2)=(jj-1+shiftk(2,ikshft))/ngkpt(2)
             k1(3)=(kk-1+shiftk(3,ikshft))/ngkpt(3)
!            Wrap the trial values in the interval ]-1/2,1/2] .
             call wrap2_pmhalf(k1(1),k2(1),shift)
             call wrap2_pmhalf(k1(2),k2(2),shift)
             call wrap2_pmhalf(k1(3),k2(3),shift)
!            Test whether it is inside the FCC BZ.
             ktest(1)=2*k2(1)-1.0d-10
             ktest(2)=2*k2(2)-2.0d-10
             ktest(3)=2*k2(3)-5.0d-10
             if (abs(ktest(1))+abs(ktest(2))+abs(ktest(3))<1.5_dp) then
               kcar(1)=ktest(1)+1.0d-10
               kcar(2)=ktest(2)+2.0d-10
               kcar(3)=ktest(3)+5.0d-10
               spkpt(1,nn)=0.5_dp*kcar(2)+0.5_dp*kcar(3)
               spkpt(2,nn)=0.5_dp*kcar(1)+0.5_dp*kcar(3)
               spkpt(3,nn)=0.5_dp*kcar(1)+0.5_dp*kcar(2)
               nn=nn+1
             end if
           end do
         end do
       end do
     end do
     nkpt=nn-1
     if(nkpt/=ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk/2)then
       write(message, '(a,i8,a,a,a,i8,a)' )&
&       '  The number of k points ',nkpt,'  is not equal to',ch10,&
&       '  (ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/2 which is',&
&       (ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/2,'.'
       MSG_BUG(message)
     end if
   end if

 else if(brav==3)then

!  Body-Centered Lattice (not mandatory cubic !)
   call wrtout(std_out,'       Body-Centered Lattice Grid ','COLL')
   if (mkpt<ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk/4) then
     write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&     '  The value of mkpt is not large enough. It should be',ch10,&
&     '  at least',(ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/4,',',ch10,&
&     '  Action : set mkpt to that value in the main routine,',ch10,&
&     '  and recompile the code.'
     MSG_BUG(message)
   end if
   nn=1
   if ((ngkpt(1)*nshiftk)/=(((ngkpt(1)*nshiftk)/2)*2) .or.&
&   (ngkpt(2)*nshiftk)/=(((ngkpt(2)*nshiftk)/2)*2) .or.&
&   (ngkpt(3)*nshiftk)/=(((ngkpt(3)*nshiftk)/2)*2) ) then
     write(message, '(4a,3(a,i6,a),a)' )&
&     '  For body-centered lattices, the numbers ngqpt(1:3)',ch10,&
&     '  must be even, while they are :',ch10,&
&     '  ngqpt(1)*nshiftk = ',ngkpt(1)*nshiftk,ch10,&
&     '  ngqpt(2)*nshiftk = ',ngkpt(2)*nshiftk,ch10,&
&     '  ngqpt(3)*nshiftk = ',ngkpt(3)*nshiftk,ch10,&
&     '  Action : modify ngqpt(1:3) in the input file.'
     MSG_ERROR(message)
   end if
   if (ngkpt(1)==0.or.ngkpt(2)==0.or.ngkpt(3)==0) then
     spkpt(1,1)=0.0_dp
     spkpt(2,1)=0.0_dp
     spkpt(3,1)=0.0_dp
     nkpt=1
   else
     do kk=1,ngkpt(3)
       do jj=1,ngkpt(2)
         do ii=1,ngkpt(1)
           do ikshft=1,nshiftk
             k1(1)=(ii-1+shiftk(1,ikshft))/ngkpt(1)
             k1(2)=(jj-1+shiftk(2,ikshft))/ngkpt(2)
             k1(3)=(kk-1+shiftk(3,ikshft))/ngkpt(3)
!            Wrap the trial values in the interval ]-1/2,1/2] .
             call wrap2_pmhalf(k1(1),k2(1),shift)
             call wrap2_pmhalf(k1(2),k2(2),shift)
             call wrap2_pmhalf(k1(3),k2(3),shift)
!            Test whether it is inside the BCC BZ.
             ktest(1)=2*k2(1)-1.0d-10
             ktest(2)=2*k2(2)-2.0d-10
             ktest(3)=2*k2(3)-5.0d-10
             if (abs(ktest(1))+abs(ktest(2))<1._dp) then
               if (abs(ktest(1))+abs(ktest(3))<1._dp) then
                 if (abs(ktest(2))+abs(ktest(3))<1._dp) then
                   kcar(1)=ktest(1)+1.0d-10
                   kcar(2)=ktest(2)+2.0d-10
                   kcar(3)=ktest(3)+5.0d-10
                   spkpt(1,nn)=-0.5*kcar(1)+0.5*kcar(2)+0.5*kcar(3)
                   spkpt(2,nn)=0.5*kcar(1)-0.5*kcar(2)+0.5*kcar(3)
                   spkpt(3,nn)=0.5*kcar(1)+0.5*kcar(2)-0.5*kcar(3)
                   nn=nn+1
                 end if
               end if
             end if
           end do
         end do
       end do
     end do
     nkpt=nn-1
     if(nkpt==0)then
       write(message, '(3a)' )&
&       '  BCC lattice, input ngqpt=0, so no kpt is generated.',ch10,&
&       '  Action : modify ngqpt(1:3) in the input file.'
       MSG_ERROR(message)
     end if
     if(nkpt/=(ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/4)then
       write(message, '(a,i8,a,a,a,i8,a)' )&
&       '  The number of k points ',nkpt,'  is not equal to',ch10,&
&       '  (ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/4 which is',&
&       (ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)/4,'.'
       MSG_BUG(message)
     end if
   end if

 else if(brav==4)then

!  Hexagonal Lattice  (D6h)
   call wrtout(std_out,'       Hexagonal Lattice Grid ','COLL')
   if (mkpt<ngkpt(1)*ngkpt(2)*ngkpt(3)) then
     write(message, '(a,a,a,i8,a,a,a,a,a)' )&
&     '  The value of mkpt is not large enough. It should be',ch10,&
&     '  at least',ngkpt(1)*ngkpt(2)*ngkpt(3),',',ch10,&
&     '  Action : set mkpt to that value in the main routine,',ch10,&
&     '  and recompile the code.'
     MSG_BUG(message) 
   end if
   nn=1
   if (ngkpt(1)/=ngkpt(2)) then
     write(message, '(4a,2(a,i6,a),a)' )&
&     '  For hexagonal lattices, the numbers ngqpt(1:2)',ch10,&
&     '  must be equal, while they are :',ch10,&
&     '  ngqpt(1) = ',ngkpt(1),ch10,&
&     '  ngqpt(2) = ',ngkpt(2),ch10,&
&     '  Action : modify ngqpt(1:3) in the input file.'
     MSG_ERROR(message)
   end if
   if (ngkpt(1)==0.or.ngkpt(2)==0.or.ngkpt(3)==0) then
     write(message, '(3a)' )&
&     '  For hexagonal lattices, ngqpt(1:3)=0 is not permitted',ch10,&
&     '  Action : modify ngqpt(1:3) in the input file.'
     MSG_ERROR(message)
   else
     do kk=1,ngkpt(3)
       do jj=1,ngkpt(2)
         do ii=1,ngkpt(1)
           do ikshft=1,nshiftk
             k1(1)=(ii-1+shiftk(1,ikshft))/ngkpt(1)
             k1(2)=(jj-1+shiftk(2,ikshft))/ngkpt(2)
             k1(3)=(kk-1+shiftk(3,ikshft))/ngkpt(3)
!            Wrap the trial values in the interval ]-1/2,1/2] .
             call wrap2_pmhalf(k1(1),k2(1),shift)
             call wrap2_pmhalf(k1(2),k2(2),shift)
             call wrap2_pmhalf(k1(3),k2(3),shift)
             spkpt(:,nn)=k2(:)
             nn=nn+1
           end do
         end do
       end do
     end do
     nkpt=nn-1
     if(nkpt/=ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk)then
       write(message, '(a,i8,a,a,a,i8,a)' )&
&       '  The number of k points ',nkpt,'  is not equal to',ch10,&
&       '  ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk which is',&
&       ngkpt(1)*ngkpt(2)*ngkpt(3)*nshiftk,'.'
       MSG_BUG(message)
     end if
   end if

 else

   write(message, '(a,i6,a,a,a)' )&
&   '  The calling routine asks brav=',brav,'.',ch10,&
&   '  but only brav=1,2,3 or 4 are allowed.'
   MSG_BUG(message)
 end if

 if(option/=0)then

!  Put the Gamma point first
   if(nkpt>1)then
     do ii=1,nkpt
       if(sum(abs(spkpt(:,ii)))<tol8)then
         spkpt(:,ii)=spkpt(:,1)   
         spkpt(:,1)=zero
         exit
       end if
     end do
   end if

   write(message,'(a,i8)')' Grid q points  : ',nkpt
   call wrtout(iout,message,'COLL')
   nkpout=nkpt
   if(nkpt>80)then
     call wrtout(iout,' greater than 80, so only write 20 of them ','COLL')
     nkpout=20
   end if
   do ii=1,nkpout
     write(message, '(1x,i2,a2,3es16.8)' )ii,') ',spkpt(1,ii),spkpt(2,ii),spkpt(3,ii)
     call wrtout(iout,message,'COLL')
   end do

 end if

end subroutine smpbz
!!***
