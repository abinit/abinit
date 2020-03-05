!!****m* ABINIT/m_getshell
!! NAME
!!  m_getshell
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2020 ABINIT group (MVeithen)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_getshell

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_linalg_interfaces

 use m_kpts,            only : getkgrid

 implicit none

 private
!!***

 public :: getshell
!!***

contains
!!***

!!****f* ABINIT/getshell
!! NAME
!! getshell
!!
!! FUNCTION
!! For each k-point, set up the shells of first neighbours and find
!! the weigths required for the finite difference expression
!! of Marzari and Vanderbilt (see PRB 56, 12847 (1997) [[cite:Marzari1997]]).
!!
!! INPUTS
!! gmet(3,3) = metric tensor of reciprocal space
!! kptopt = option for the generation of k points
!! kptrlatt = k-point lattice specification
!! kpt2(3,nkpt2) = reduced coordinates of the k-points in the
!!                 reduced part of the BZ (see below)
!! mkmem = number of k points which can fit in memory
!! nkpt2 = number of k-points in the reduced BZ
!! nkpt3 = number of k-points in the full BZ
!! nshiftk = number of kpoint grid shifts
!! rmet(3,3) = metric tensor of real space
!! rprimd(3,3) = dimensional primitive translations (bohr)
!! shiftk = shift vectors for k point generation
!! wtk2 = weight assigned to each k point
!! comm=MPI communicator
!!
!! OUTPUT
!! kneigh(30,nkpt2) = for each k-point in the reduced part of the BZ
!!                    kneigh stores the index (ikpt) of the neighbouring
!!                    k-points
!! kg_neigh(30,nkpt2,3) = kg-neigh takes values of -1, 0 or 1,
!!                        and can be non-zero only for a single k-point,
!!                        a line of k-points or a plane of k-points.
!!                        The vector joining the ikpt2-th k-point to its
!!                        ineigh-th nearest neighbour is :
!!                        dk(:)-nint(dk(:))+real(kg_neigh(ineigh,ikpt2,:))
!!                        with dk(:)=kpt2(:,kneigh(ineigh,ikpt2))-kpt2(:,ikpt2)
!! kptindex(2,nkpt3)
!!   kptindex(1,ikpt) = ikpt_rbz
!!     ikpt_rbz = index of the k-point in the reduced BZ
!!     ikpt = index of the k-point in the full BZ
!!   kptindex(2,ikpt) = 1: use time-reversal symmetry to transform the
!!                         wavefunction at ikpt_rbz to the wavefunction at ikpt
!!                      0: ikpt belongs already to the reduced BZ
!!                         (no transformation required)
!! kpt3(3,nkpt3) = reduced coordinates of the k-points in the full BZ
!! mvwtk(30,nkpt2) = weights required to evaluate the finite difference
!!                   formula of Marzari and Vanderbilt, computed for each
!!                   k-point in the reduced part of the BZ
!! mkmem_max = maximal number of k-points on each processor (MPI //)
!! nneigh = total number of neighbours required to evaluate the finite
!!          difference formula
!!
!! COMMENTS
!! The array kpt2 holds the reduced coordinates of the k-points in the
!! reduced part of the BZ. For example, in case time-reversal symmetry is
!! used (kptopt = 2) kpt2 samples half the BZ. Since some of the neighbours
!! of these k-points may lie outside the reduced BZ, getshell also needs the
!! coordinates of the k-points in the full BZ.
!! The coordinates of the k-points in the full BZ are stored in kpt3.
!! The weights mvwtk are computed for the k-points kpt2.
!!
!! In case no symmetry is used to reduce the number of k-points,
!! the arrays kpt2 and kpt3 are equal.
!!
!! PARENTS
!!      nonlinear
!!
!! CHILDREN
!!      dgelss,getkgrid,wrtout,xmpi_max
!!
!! SOURCE

subroutine getshell(gmet,kneigh,kg_neigh,kptindex,kptopt,kptrlatt,kpt2,&
& kpt3,mkmem,mkmem_max,mvwtk,&
& nkpt2,nkpt3,nneigh,nshiftk,rmet,rprimd,shiftk,wtk2, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt,mkmem,nkpt2,nkpt3,comm
 integer,intent(inout) :: nshiftk
 integer,intent(out) :: mkmem_max,nneigh
!arrays
 integer,intent(inout) :: kptrlatt(3,3)
 integer,intent(out) :: kneigh(30,nkpt2),kptindex(2,nkpt3),kg_neigh(30,nkpt2,3)
 real(dp),intent(in) :: gmet(3,3),kpt2(3,nkpt2),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk),wtk2(nkpt2)
 real(dp),intent(out) :: kpt3(3,nkpt3),mvwtk(30,nkpt2)

!Local variables-------------------------------
!scalars
 integer :: bis,flag,ier,ii,ikpt,ikpt2,ikpt3,ineigh,info,irank,is1,ishell
 integer :: jj,kptopt_used,mkmem_cp,nkpt_computed,nshell,nsym1,orig
 integer :: wtkflg, coord1, coord2, coord3
 real(dp) :: dist_,kptrlen,last_dist,max_dist,resdm,s1
 character(len=500) :: message
!arrays
 integer :: neigh(0:6,nkpt2),symafm_dummy(1),vacuum(3)
 integer,allocatable :: symrel1(:,:,:)
 real(dp) :: dist(6),dk(3),dk_(3),mat(6,6),rvec(6),sgval(6)
 real(dp) :: shiftk_(3,MAX_NSHIFTK),work(30)
 real(dp),allocatable :: tnons1(:,:),wtk3(:)

!************************************************************************

!In case of MPI //: compute maximum number of k-points per processor
 if (xmpi_paral == 1) then
   mkmem_cp=mkmem
   call xmpi_max(mkmem_cp,mkmem_max,comm,ier)
 else
   mkmem_max = mkmem
 end if

!------------- In case kptopt = 2 set up the whole k-point grid -------------

!kpt3(3,nkpt3) = reduced coordinates of k-points in the full BZ

 if (kptopt == 3) then

   ABI_ALLOCATE(wtk3,(nkpt3))
   kpt3(:,:) = kpt2(:,:)
   wtk3(:) = wtk2(:)
   do ikpt = 1,nkpt3
     kptindex(1,ikpt) = ikpt
     kptindex(2,ikpt) = 0
   end do

 else if (kptopt == 2) then

   ABI_ALLOCATE(wtk3,(nkpt3))
   ii = 5 ; kptopt_used = 3
   symafm_dummy(1) = 1
   shiftk_(:,:) = 0._dp
   shiftk_(:,1:nshiftk) = shiftk(:,1:nshiftk)

   nsym1 = 1
   ABI_ALLOCATE(symrel1,(3,3,nsym1))
   ABI_ALLOCATE(tnons1,(3,nsym1))
   symrel1(:,:,1) = 0
   symrel1(1,1,1) = 1 ; symrel1(2,2,1) = 1 ; symrel1(3,3,1) = 1
   tnons1(:,:) = 0._dp
   vacuum(:) = 0

   call getkgrid(0,0,ii,kpt3,kptopt_used,kptrlatt,&
&   kptrlen,nsym1,nkpt3,nkpt_computed,nshiftk,nsym1,&
&   rprimd,shiftk_,symafm_dummy,symrel1,&
&   vacuum,wtk3)

   if (nkpt_computed /= nkpt3) then
     write(message,'(a,a,a,a,i4,a,a,i4)')&
&     ' The number of k-points in the whole BZ, nkpt_computed= ',nkpt_computed,ch10,&
&     ' is not twice the number of k-points in half the BZ, nkpt3=',nkpt3
     MSG_BUG(message)
   end if

   kptindex(:,:) = 0
   do ikpt3 = 1, nkpt3

     flag = 1
     do ikpt2 = 1, nkpt2

!      In case the k-points differ only by one reciprocal lattice
!      vector, apply shift of one g-vector to kpt(:,ikpt3)
!     MJV 10/2019: this appears to be using time reversal sym, instead of the G vector...
!      could be equivalent to keep points inside the 1st BZ, but the code below is not consistent
!      with this comment

!
!  here k3 = k2 + G 
!
       dk_(:) = kpt3(:,ikpt3) - kpt2(:,ikpt2)
       dk(:) = dk_(:) - nint(dk_(:))
       if (dk(1)*dk(1) + dk(2)*dk(2) + dk(3)*dk(3) < tol10) then
         do ii = 1, 3
           if ((dk(ii)*dk(ii) < tol10).and.(dk_(ii)*dk_(ii) > tol10)) then
!  transform k3 to -k3
!  TODO: I suspect this should be k3 -= G!!
             kpt3(ii,ikpt3) = -1._dp*kpt3(ii,ikpt3)
           end if
         end do
       end if

!
! here k3 = -k2 + G
!
       dk_(:) = kpt3(:,ikpt3) + kpt2(:,ikpt2)
       dk(:) = dk_(:) - nint(dk_(:))
       if (dk(1)*dk(1) + dk(2)*dk(2) + dk(3)*dk(3) < tol10) then
         do ii = 1, 3
           if ((dk(ii)*dk(ii) < tol10).and.(dk_(ii)*dk_(ii) > tol10)) then
!  transform k3 to -k3
!  TODO: I suspect this should be k3 -= G!!
             kpt3(ii,ikpt3) = -1._dp*kpt3(ii,ikpt3)
           end if
         end do
       end if


!
! here k3 = k2
!
       dk(:) = kpt3(:,ikpt3) - kpt2(:,ikpt2)
       if (dk(1)*dk(1) + dk(2)*dk(2) + dk(3)*dk(3) < tol10) then
         kptindex(1,ikpt3) = ikpt2
         kptindex(2,ikpt3) = 0       ! no use of time-reversal symmetry
         flag = 0
         exit
       end if

!
! here k3 = -k2
!
       dk(:) = kpt3(:,ikpt3) + kpt2(:,ikpt2)
       if (dk(1)*dk(1) + dk(2)*dk(2) + dk(3)*dk(3) < tol10) then
         kptindex(1,ikpt3) = ikpt2
         kptindex(2,ikpt3) = 1       ! use time-reversal symmetry
         flag = 0
         exit
       end if

     end do     ! ikpt2

     if (flag == 1) then
       write(message,'(a,i0)')' Could not find a symmetric k-point for ikpt3=  ',ikpt3
       MSG_BUG(message)
     end if
   end do    ! ikpt3

 else
   message = ' the only values for kptopt that are allowed are 2 and 3 '
   MSG_ERROR(message)
 end if   ! condition on kptopt


!--------- Compute the weights required for the Marzari-Vanderbilt ---------
!--------- finite difference formula ---------------------------------------


!Initialize distance between k-points
!The trace of gmet is an upper limit for its largest eigenvalue. Since the
!components of the distance vectors do not exceed 1, 3. * Tr[gmet] is
!an upper limit for the squared shell radius.
!we take something two times larger to make a bug checking.
 dist_ = 0._dp
 do ii = 1,3
   dist_ = dist_ + gmet(ii,ii)
 end do
 max_dist = 3._dp * dist_ * 2._dp
 write(std_out,*)'max_dist',max_dist

!Calculate an upper limit for the residuum
 resdm = rmet(1,1)*rmet(1,1) + rmet(2,2)*rmet(2,2) + rmet(3,3)*rmet(3,3)&
& + rmet(1,2)*rmet(1,2) + rmet(2,3)*rmet(2,3) + rmet(3,1)*rmet(3,1)

!Initialize shell loop
 ishell = 0
 last_dist = 0._dp
 wtkflg = 0
 kneigh(:,:) = 0
 kg_neigh(:,:,:) = 0
 neigh(:,:) = 0

!Loop over shells until the residuum is zero
 do while ((wtkflg == 0).and.(resdm > tol8))
!  Advance shell counter
   ishell = ishell + 1

!  Initialize shell radius with upper limit
   dist(ishell) = max_dist
!  !!  border_flag = 1

!  !write(std_out,*)'gmet'
!  !do ikpt=1,3
!  !write(std_out,*)gmet(ikpt,:)
!  !enddo
!  !write(std_out,*)kpt3(:,1)

!  Find the (squared) radius of the next shell
   do ikpt = 1,nkpt3
!    !write(std_out,*)ikpt
!    !write(std_out,*)kpt3(:,ikpt)
     dk(:) = kpt3(:,1) - kpt3(:,ikpt)
!    !!dk_(:) = dk(:) - nint(dk(:))
!    !!dist_ = 0._dp
!    !!do ii = 1,3
!    !! do jj = 1,3
!    !!  dist_ = dist_ + dk_(ii)*gmet(ii,jj)*dk_(jj)
!    !! end do
!    !!end do
!    !!write(std_out,*)'dist_1', dist_
!    !!   dist_ = 0._dp
!    !!   do ii = 1,3
!    !!    do jj = 1,3
!    !!     dist_ = dist_ + dk(ii)*gmet(ii,jj)*dk(jj)
!    !!    end do
!    !!   end do
!    !!write(std_out,*)'dist_2',dist_
     do coord1 = 0,1  !three loop to search also on the border of the BZ, ie for the k-points (1,k2,k3) and the likes
       do coord2 = 0,1
         do coord3 = 0,1
!          !!      if ((coord1/=0).or.(coord2/=0).or.(coord3/=0)) then
           dist_ = 0._dp
           dk_(:) = dk(:) - nint(dk(:))
           dk_(1) = dk_(1) + real(coord1,dp)
           dk_(2) = dk_(2) + real(coord2,dp)
           dk_(3) = dk_(3) + real(coord3,dp)
           do ii = 1,3
             do jj = 1,3
               dist_ = dist_ + dk_(ii)*gmet(ii,jj)*dk_(jj)
             end do
           end do
!          Note : for ipkt3 = 1, coord1 = coord2 = coord3 = 0, the distance is 0 ;
!          but the next "if" statement is false with the tol8 criteria and the k-point
!          should be ignored even for ishell = 1 and last_dist= 0.
!          !$write(std_out,*)ikpt,coord1,coord2,coord3
!          !$write(std_out,*)dk_
!          !$write(std_out,*)'dist_2', dist_
!          !!      end if
           if ((dist_ < dist(ishell)).and.(dist_ - last_dist > tol8)) then
             dist(ishell) = dist_
           end if
         end do
       end do
     end do

!    !!   if ((dist_ < dist(ishell)).and.(dist_ - last_dist > tol8)) then
!    !!    dist(ishell) = dist_
!    !!    border_flag = 0
!    !!   end if
   end do

!  !!  if (border_flag==1) then !we haven't found any shell in the interior of the BZ, we need to search on the border
!  !!write(std_out,*)ch10
!  !!write(std_out,*)'search on the border'
!  !!   do ikpt = 1,nkpt3
!  !!    dk(:) = kpt3(:,1) - kpt3(:,ikpt)
!  !!    do coord1 = 0,1
!  !!     do coord2 = 0,1
!  !!      do coord3 = 0,1
!  !!       if ((coord1/=0).or.(coord2/=0).or.(coord3/=0)) then
!  !!        dist_ = 0._dp
!  !!        dk_(:) = dk(:) - nint(dk(:))
!  !!        dk_(1) = dk_(1) + real(coord1,dp)
!  !!        dk_(2) = dk_(2) + real(coord2,dp)
!  !!        dk_(3) = dk_(3) + real(coord3,dp)
!  !!        do ii = 1,3
!  !!         do jj = 1,3
!  !!          dist_ = dist_ + dk_(ii)*gmet(ii,jj)*dk_(jj)
!  !!         end do
!  !!        end do
!  !!write(std_out,*)ikpt,coord1,coord2,coord3
!  !!write(std_out,*)dk_
!  !!write(std_out,*)'dist_2', dist_
!  !!       end if
!  !!       if ((dist_ < dist(ishell)).and.(dist_ - last_dist > tol8)) then
!  !!        dist(ishell) = dist_
!  !!       end if
!  !!      end do
!  !!     end do
!  !!    end do
!  !!   end do
!  !!  endif

!  DEBUG
!  !write(std_out,*)ch10
!  write(std_out,*)'ishell, dist = ',ishell,dist(ishell)
!  ENDDEBUG

   if (max_dist-dist(ishell)<tol8) then
     write(message,'(a,i0)')' Cannot find shell number',ishell
     MSG_BUG(message)
   end if

   last_dist = dist(ishell)

!  For each k-point in half the BZ get the shells of nearest neighbours.
!  These neighbours can be out of the zone sampled by kpt2.
!  !$write(std_out,*)'nkpt2', nkpt2, 'nkpt3', nkpt3
   do ikpt2 = 1, nkpt2              ! k-points in half the BZ
     orig = sum(neigh(0:ishell-1,ikpt2))
!    !write(std_out,*)'ikpt2, orig', ikpt2,orig
!    !write(std_out,*) kpt2(:,ikpt2)
     nneigh = 0
     do ikpt3 = 1, nkpt3             ! whole k-point grid
!      !!    if(border_flag==0)then
       dk(:) = kpt3(:,ikpt3) - kpt2(:,ikpt2)
!      !!     dk_(:) = dk(:) - nint(dk(:))
!      !!     dist_ = 0._dp
       do coord1 = -1,1
         do coord2 = -1,1
           do coord3 = -1,1
!            !!        if ((coord1/=0).or.(coord2/=0).or.(coord3/=0)) then
             dist_ = 0._dp
             dk_(:) = dk(:) - nint(dk(:))
             dk_(1) = dk_(1) + real(coord1,dp)
             dk_(2) = dk_(2) + real(coord2,dp)
             dk_(3) = dk_(3) + real(coord3,dp)
             do ii = 1,3
               do jj = 1,3
                 dist_ = dist_ + dk_(ii)*gmet(ii,jj)*dk_(jj)
               end do
             end do
             if (abs(dist_ - dist(ishell)) < tol8) then
               nneigh = nneigh + 1
               kneigh(orig+nneigh,ikpt2) = ikpt3
               kg_neigh(orig+nneigh,ikpt2,1) = coord1
               kg_neigh(orig+nneigh,ikpt2,2) = coord2
               kg_neigh(orig+nneigh,ikpt2,3) = coord3
             end if
!            !!        end if
           end do
         end do
       end do
!      !write(std_out,*)'ikpt3', ikpt3
!      !write(std_out,*) kpt3(:,ikpt3)
!      write(std_out,*) kpt2(:,ikpt2)
!      !write(std_out,*) dk
!      write(std_out,*) dk_
!      !!     do ii = 1,3
!      !!      do jj = 1,3
!      !!       dist_ = dist_ + dk_(ii)*gmet(ii,jj)*dk_(jj)
!      !!      end do
!      !!     end do
!      !write(std_out,*)'dist_', dist_
!      !!     if (abs(dist_ - dist(ishell)) < tol8) then
!      !!      nneigh = nneigh + 1
!      !!      kneigh(orig+nneigh,ikpt2) = ikpt3
!      !!     end if
!      !!    else !search on the border
!      !!     dk(:) = kpt3(:,ikpt3) - kpt2(:,ikpt2)
!      !!     do coord1 = -1,1
!      !!      do coord2 = -1,1
!      !!       do coord3 = -1,1
!      !!        if ((coord1/=0).or.(coord2/=0).or.(coord3/=0)) then
!      !!         dist_ = 0._dp
!      !!         dk_(:) = dk(:) - nint(dk(:))
!      !!         dk_(1) = dk_(1) + real(coord1,dp)
!      !!         dk_(2) = dk_(2) + real(coord2,dp)
!      !!         dk_(3) = dk_(3) + real(coord3,dp)
!      !!         do ii = 1,3
!      !!          do jj = 1,3
!      !!           dist_ = dist_ + dk_(ii)*gmet(ii,jj)*dk_(jj)
!      !!          end do
!      !!         end do
!      !!         if (abs(dist_ - dist(ishell)) < tol8) then
!      !!          nneigh = nneigh + 1
!      !!          kneigh(orig+nneigh,ikpt2) = ikpt3
!      !!          kneigh_border(orig+nneigh,ikpt2,1) = real(coord1,dp)
!      !!          kneigh_border(orig+nneigh,ikpt2,2) = real(coord2,dp)
!      !!          kneigh_border(orig+nneigh,ikpt2,3) = real(coord3,dp)
!      !!         end if
!      !!        end if
!      !!       end do
!      !!      end do
!      !!     end do
!      !!    end if
     end do
     neigh(ishell,ikpt2) = nneigh
   end do


!  Check if the number of points in shell number ishell
!  is the same for each k-point

   flag = 1
   do ikpt = 1,nkpt2
     if (neigh(ishell,ikpt) /= nneigh) flag = 0
   end do

   if (flag == 0) then
     write(message,'(a,i0,a,a)')&
&     ' The number of points in shell number',ishell,' is not the same',&
&     ' for each k-point.'
     MSG_BUG(message)
   end if

   if (nneigh == 0) then
     write(message,'(a,a,a,a)') ch10,&
&     ' getshell: BUG - ',ch10,&
&     ' Cannot find enough neighbor shells'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     wtkflg = 1
   end if

!  Calculate the total number of neighbors
   nneigh = sum(neigh(1:ishell,1))
!  DEBUG
   write(std_out,*)'ishell = ',ishell,'nneigh = ',nneigh
!  ENDDEBUG

!  Find the weights needed to compute the finite difference expression
!  of the ddk
!  **********************************************************************

!  mvwtk(:,:) = 0._dp

!  The weights are calculated for ikpt=1. The results are copied later
   ikpt = 1

!  Calculate the coefficients of the linear system to be solved
   mat(:,:) = 0._dp
   do is1 = 1, ishell
     orig = sum(neigh(0:is1-1,ikpt))
     bis = orig + neigh(is1,ikpt)
     do ineigh = orig+1, bis
       dk_(:) = kpt3(:,kneigh(ineigh,ikpt)) - kpt2(:,ikpt)
       dk(:) = dk_(:) - nint(dk_(:))
       dk(:) = dk(:) + real(kg_neigh(ineigh,ikpt,:),dp)
       mat(1,is1) = mat(1,is1) + dk(1)*dk(1)
       mat(2,is1) = mat(2,is1) + dk(2)*dk(2)
       mat(3,is1) = mat(3,is1) + dk(3)*dk(3)
       mat(4,is1) = mat(4,is1) + dk(1)*dk(2)
       mat(5,is1) = mat(5,is1) + dk(2)*dk(3)
       mat(6,is1) = mat(6,is1) + dk(3)*dk(1)
     end do
   end do

   rvec(1) = rmet(1,1)
   rvec(2) = rmet(2,2)
   rvec(3) = rmet(3,3)
   rvec(4) = rmet(1,2)
   rvec(5) = rmet(2,3)
   rvec(6) = rmet(3,1)

!  DEBUG
   write(std_out,*) " mat(1:6, 1:ishell) : rmet(1:6) for all 6 products dx^2... dxdy..."
   do ii = 1, 6
     write(std_out,*) mat(ii,1:ishell), ' : ', rvec(ii)
   end do
!  ENDDEBUG

!  Solve the linear least square problem
   call dgelss(6,ishell,1,mat,6,rvec,6,sgval,tol8,irank,work,30,info)

   if( info /= 0 ) then
     write(message,'(3a,i0,a)')&
&     ' Singular-value decomposition of the linear system determining the',ch10,&
&     ' weights failed (info).',info,ch10
     MSG_COMMENT(message)
     wtkflg = 1
   end if

!  Check that the system has maximum rank
   if( irank == ishell ) then
!    System has full rank. Calculate the residuum
     s1 = resdm
     resdm = 0._dp
     do is1 = ishell + 1, 6
       resdm = resdm + rvec(is1) * rvec(is1)
     end do

     if( ishell == 6 .and. resdm > tol8 ) then
       write(message,'(4a)')&
&       ' Linear system determining the weights could not be solved',ch10,&
&       ' This should not happen.',ch10
       MSG_COMMENT(message)
       wtkflg = 1
     end if
   else
!    The system is rank deficient
     ishell = ishell - 1
!    DEBUG
     write(std_out,*) 'Shell not linear independent from previous shells. Skipped.'
!    ENDDEBUG
   end if

!  DEBUG
   write(std_out,*) "ishell, nneigh, irank, resdm ", ishell, nneigh, irank, resdm
!  ENDDEBUG

!  end of loop over shells
 end do

!Copy weights
 ikpt=1
 do is1 = 1, ishell
   orig = sum(neigh(0:is1-1,ikpt))
   bis = orig + neigh(is1,ikpt)
   mvwtk(orig+1:bis,1) = rvec(is1)
 end do
 do ikpt = 2,nkpt2
   mvwtk(1:nneigh,ikpt) = mvwtk(1:nneigh,1)
 end do  ! ikpt

!Report weights
 write(std_out,*) 'Neighbors(1:ishell,1) ', neigh(1:ishell,1)
 write(std_out,*) 'Weights (1:ishell) ', rvec(1:ishell)
 write(std_out,*) mvwtk(1:nneigh,1)

!Check the computed weights
 if (wtkflg == 0) then
   do ikpt = 1, nkpt2
     do ii = 1,3
       do jj = 1,3
         s1 = 0._dp
         do ineigh = 1, nneigh
           dk_(:) = kpt3(:,kneigh(ineigh,ikpt)) - kpt2(:,ikpt)
           dk(:) = dk_(:) - nint(dk_(:))
           dk(:) = dk(:) + real(kg_neigh(ineigh,ikpt,:),dp)
           s1 = s1 + dk(ii)*dk(jj)*mvwtk(ineigh,ikpt)
         end do
         if (abs(s1 - rmet(ii,jj)) > tol6) then
           wtkflg = 1
         end if
       end do
     end do
   end do

   if (wtkflg /= 0) then
     write(message,'(a,a,a,a)') ch10,&
&     ' getshell : BUG -',ch10,&
&     ' The calculated weights do not solve the linear system for all k-points.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if
 end if

 if (wtkflg /= 0) then

   message = ' There is a problem with the finite difference expression of the ddk '//ch10&
&        //' If you are very close to a symmetric structure, you might be confusing the algorithm with'//ch10&
&        //' sets of k-points which are not quite part of the same shell. Try rectifying angles and acell.'
   MSG_BUG(message)

 else

   nshell = ishell

   write(message,'(a,a,a,a,a,a,a,i3,a,a,f16.7)') ch10,&
&   ' getshell : finite difference formula of Marzari and Vanderbilt',ch10,&
&   '            (see Marzari and Vanderbilt, PRB 56, 12847 (1997), Appendix B)',& ! [[cite:Marzari1997]]
&   ch10,ch10,&
&   '            number of first neighbours  : ', neigh(1,1),ch10,&
&   '            weight : ',mvwtk(1,1)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   if (nshell > 1) then
     is1 = neigh(1,1) + 1
     write(message,'(a,a,i3,a,a,f16.7)')ch10,&
&     '            number of second neighbours  : ', neigh(2,1),ch10,&
&     '            weight : ',mvwtk(is1,1)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if

   if (nshell > 2) then
     is1 = sum(neigh(1:2,1)) + 1
     write(message,'(a,a,i3,a,a,f16.7)')ch10,&
&     '            number of third neighbours  : ', neigh(3,1),ch10,&
&     '            weight : ',mvwtk(is1,1)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if

   if (nshell > 3) then
     is1 = sum(neigh(1:3,1)) + 1
     write(message,'(a,a,i3,a,a,f16.7)')ch10,&
&     '            number of fourth neighbours  : ', neigh(4,1),ch10,&
&     '            weight : ',mvwtk(is1,1)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if

   if (nshell > 4) then
     is1 = sum(neigh(1:4,1)) + 1
     write(message,'(a,a,i3,a,a,f16.7)')ch10,&
&     '            number of fifth neighbours  : ', neigh(5,1),ch10,&
&     '            weight : ',mvwtk(is1,1)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if

   if (nshell > 5) then
     is1 = sum(neigh(1:5,1)) + 1
     write(message,'(a,a,i3,a,a,f16.7)')ch10,&
&     '            number of sixth neighbours  : ', neigh(6,1),ch10,&
&     '            weight : ',mvwtk(is1,1)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if

 end if

 if (allocated(tnons1))  then
   ABI_DEALLOCATE(tnons1)
 end if
 if (allocated(symrel1))  then
   ABI_DEALLOCATE(symrel1)
 end if

 ABI_DEALLOCATE(wtk3)

end subroutine getshell
!!***

end module m_getshell
!!***
