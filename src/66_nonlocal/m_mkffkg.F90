!!****m* ABINIT/m_mkffkg
!! NAME
!!  m_mkffkg
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, MT, DRH)
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

module m_mkffkg

 use defs_basis
 use m_abicore

 implicit none

 private
!!***

 public :: mkffkg
 public :: dfpt_mkffkg
!!***

contains
!!***

!!****f* ABINIT/dfpt_mkffkg
!! NAME
!! dfpt_mkffkg
!!
!! FUNCTION
!! Prepare the application of the projectors to the shifted wavefunctions,
!! by precomputing the k+G factors and their product with the form factors
!! Do this on a block of plane waves.
!!
!! INPUTS
!!  choice=governs the combination of k+G vectors to be computed
!!  ffnl(npw,nffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  nffnl=3rd dimension of ffnl(2, conventional, or 3 for 2nd derivatives)
!!  idir=direction of the perturbation (needed if choice==2 and ndgxdt==1,
!!       or if choice==5)
!!  indlmn(6,i,ntypat)=array giving l,m,n,lm,ln,spin for i=ln
!!  ipw1 = index of the first plane wave treated in this block
!!  ispinor=1 or 2, gives the spinorial component of ffnl to be used
!!  itypat = type of atom, needed for ffnl
!!  kg_k(3,npw)=integer coords of planewaves in basis sphere
!!  kpg_k(npw,npkg)= (k+G) components and related data
!!  kpt(3)=real components of k point in terms of recip. translations
!!  lmnmax=max. number of (l,n) components over all type of psps
!!  mblkpw=first dimension of kpgx
!!  ndgxdt=number of components of first order derivative
!!  nffkg=number of products of ffnls with combinations of k+G
!!  nincpw=number of plane waves in the block
!!  nkpg=second size of array kpg_k
!!  nlang = number of angular momenta to be treated = 1 + highest ang. mom.
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npw  = total number of plane waves in reciprocal space
!!  ntens=second dimension of kpgx, number of distinct tensorial products
!!  ntypat = number of type of atoms, dimension needed for ffnl
!!
!! OUTPUT
!!  kpgx(mblkpw,ntens)=different tensorial products of k+G
!!  ffkg(nffkg,mblkpw)=different products of ffnls with k+G
!!  parity(nffkg)=parity of the tensorial product of k+G (2 if even, 1 of odd)
!!
!! NOTES
!!  This routine must be thread-safe as it is called inside loops that are OpenMP parallelized.
!!  Please, do not add variables with the save attribute or SIDE EFFECTS.
!!
!! PARENTS
!!      opernl3,opernl4a,opernl4b
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfpt_mkffkg(choice,ffkg,ffnl,gmet,idir,indlmn,ipw1,ispinor,itypat,&
&                  kg_k,kpg_k,kpgx,kpt,lmnmax,mblkpw,ndgxdt,nffkg,nffnl,nincpw,nkpg,nlang,&
&                  npw,ntens,ntypat,parity)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,idir,ipw1,ispinor,itypat,lmnmax,mblkpw,ndgxdt
 integer,intent(in) :: nffkg,nffnl,nincpw,nkpg,nlang,npw,ntens,ntypat
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),kg_k(3,npw)
 integer,intent(out) :: parity(nffkg)
 real(dp),intent(in) :: ffnl(npw,nffnl,lmnmax,ntypat),gmet(3,3),kpg_k(npw,nkpg)
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(out) :: ffkg(nffkg,mblkpw),kpgx(mblkpw,ntens)

!Local variables-------------------------------
!scalars
 integer :: iffkg,ig,ii,ilang,ilang2,ilangx,ilmn,iln,iln0,iproj,ipw,jj
 integer :: nffkge
 real(dp) :: ffkg_now,kpg_x,kpg_y,kpg_z

! *************************************************************************

 jj=0;ilangx=0

!This will be useless after all the modifications have been done
 do ipw=1,nincpw
   kpgx(ipw,1)=1.0d0
 end do

!Initialize kpgx array related to tensors defined below
 if ( nlang>=2 .or. choice==2 .or. choice==3 .or. choice==4 .or. choice==5&
& .or. choice==6 .or. choice==23) then
   if (nkpg>=3) then
     kpgx(1:nincpw,2)=kpg_k(ipw1+1:ipw1+nincpw,1)
     kpgx(1:nincpw,3)=kpg_k(ipw1+1:ipw1+nincpw,2)
     kpgx(1:nincpw,4)=kpg_k(ipw1+1:ipw1+nincpw,3)
   else
     ig=ipw1
     do ipw=1,nincpw
       kpgx(ipw,2)=kpt(1)+dble(kg_k(1,ig))
       kpgx(ipw,3)=kpt(2)+dble(kg_k(2,ig))
       kpgx(ipw,4)=kpt(3)+dble(kg_k(3,ig))
       ig=ig+1
     end do
   end if
 end if
 if (nlang>=3 .or. choice==3 .or. choice==6 .or. choice==23) then
!  Define (k+G) part of rank 2 symmetric tensor (6 components), l=2
!  Compressed storage is 11 22 33 32 31 21
   if (nkpg>=9) then
     kpgx(1:nincpw,5) =kpg_k(ipw1+1:ipw1+nincpw,4)
     kpgx(1:nincpw,6) =kpg_k(ipw1+1:ipw1+nincpw,5)
     kpgx(1:nincpw,7) =kpg_k(ipw1+1:ipw1+nincpw,6)
     kpgx(1:nincpw,8) =kpg_k(ipw1+1:ipw1+nincpw,7)
     kpgx(1:nincpw,9) =kpg_k(ipw1+1:ipw1+nincpw,8)
     kpgx(1:nincpw,10)=kpg_k(ipw1+1:ipw1+nincpw,9)
   else
     do ipw=1,nincpw
       kpgx(ipw, 5) =      kpgx(ipw, 2)*kpgx(ipw, 2)
       kpgx(ipw, 6) =      kpgx(ipw, 3)*kpgx(ipw, 3)
       kpgx(ipw, 7) =      kpgx(ipw, 4)*kpgx(ipw, 4)
       kpgx(ipw, 8) =      kpgx(ipw, 4)*kpgx(ipw, 3)
       kpgx(ipw, 9) =      kpgx(ipw, 4)*kpgx(ipw, 2)
       kpgx(ipw,10) =      kpgx(ipw, 3)*kpgx(ipw, 2)
     end do
   end if
 end if
 if (nlang>=4 .or. ((choice==3.or.choice==23) .and. nlang>=2) .or. choice==6) then
!  Define (k+G) part of rank 3 symmetric tensor (10 components), l=3
!  Compressed storage is 111 221 331 321 311 211 222 332 322 333
   do ipw=1,nincpw
     kpgx(ipw,11) =     kpgx(ipw, 5)*kpgx(ipw, 2)
     kpgx(ipw,12) =     kpgx(ipw, 6)*kpgx(ipw, 2)
     kpgx(ipw,13) =     kpgx(ipw, 7)*kpgx(ipw, 2)
     kpgx(ipw,14) =     kpgx(ipw, 8)*kpgx(ipw, 2)
     kpgx(ipw,15) =     kpgx(ipw, 9)*kpgx(ipw, 2)
     kpgx(ipw,16) =     kpgx(ipw,10)*kpgx(ipw, 2)
     kpgx(ipw,17) =     kpgx(ipw, 6)*kpgx(ipw, 3)
     kpgx(ipw,18) =     kpgx(ipw, 7)*kpgx(ipw, 3)
     kpgx(ipw,19) =     kpgx(ipw, 8)*kpgx(ipw, 3)
     kpgx(ipw,20) =     kpgx(ipw, 7)*kpgx(ipw, 4)
   end do
 end if
 if (((choice==3.or.choice==23) .and. nlang>=3) .or. choice==6) then
!  Add additional tensors for strain gradients
!  Define (k+G) part of rank 4 symmetric tensor (15 components), l=2
!  Compressed storage is 1111 2211 3311 3211 3111 2111 2221 3321 3221
!  3331 2222 3322 3222 3332 3333
   do ipw=1,nincpw
     kpgx(ipw,21) =     kpgx(ipw, 5)*kpgx(ipw, 5)
     kpgx(ipw,22) =     kpgx(ipw, 6)*kpgx(ipw, 5)
     kpgx(ipw,23) =     kpgx(ipw, 7)*kpgx(ipw, 5)
     kpgx(ipw,24) =     kpgx(ipw, 8)*kpgx(ipw, 5)
     kpgx(ipw,25) =     kpgx(ipw, 9)*kpgx(ipw, 5)
     kpgx(ipw,26) =     kpgx(ipw,10)*kpgx(ipw, 5)
     kpgx(ipw,27) =     kpgx(ipw, 6)*kpgx(ipw,10)
     kpgx(ipw,28) =     kpgx(ipw, 7)*kpgx(ipw,10)
     kpgx(ipw,29) =     kpgx(ipw, 8)*kpgx(ipw,10)
     kpgx(ipw,30) =     kpgx(ipw, 7)*kpgx(ipw, 9)
     kpgx(ipw,31) =     kpgx(ipw, 6)*kpgx(ipw, 6)
     kpgx(ipw,32) =     kpgx(ipw, 7)*kpgx(ipw, 6)
     kpgx(ipw,33) =     kpgx(ipw, 8)*kpgx(ipw, 6)
     kpgx(ipw,34) =     kpgx(ipw, 7)*kpgx(ipw, 8)
     kpgx(ipw,35) =     kpgx(ipw, 7)*kpgx(ipw, 7)
   end do
 end if
 if (((choice==3.or.choice==23) .and. nlang>=4) .or. (choice==6 .and. nlang>=2)) then
!  Define (k+G) part of rank 5 symmetric tensor (21 components), l=3
!  Compressed storage is 11111 22111 33111 32111 31111 21111
!  22211 33211 32211     33311 22221 33221 32221 33321 33331
!  22222 33222 32222     33322 33332 33333
   do ipw=1,nincpw
     kpgx(ipw,36) =     kpgx(ipw,21)*kpgx(ipw, 2)
     kpgx(ipw,37) =     kpgx(ipw,22)*kpgx(ipw, 2)
     kpgx(ipw,38) =     kpgx(ipw,23)*kpgx(ipw, 2)
     kpgx(ipw,39) =     kpgx(ipw,24)*kpgx(ipw, 2)
     kpgx(ipw,40) =     kpgx(ipw,25)*kpgx(ipw, 2)
     kpgx(ipw,41) =     kpgx(ipw,26)*kpgx(ipw, 2)
     kpgx(ipw,42) =     kpgx(ipw,27)*kpgx(ipw, 2)
     kpgx(ipw,43) =     kpgx(ipw,28)*kpgx(ipw, 2)
     kpgx(ipw,44) =     kpgx(ipw,29)*kpgx(ipw, 2)
     kpgx(ipw,45) =     kpgx(ipw,30)*kpgx(ipw, 2)
     kpgx(ipw,46) =     kpgx(ipw,31)*kpgx(ipw, 2)
     kpgx(ipw,47) =     kpgx(ipw,32)*kpgx(ipw, 2)
     kpgx(ipw,48) =     kpgx(ipw,33)*kpgx(ipw, 2)
     kpgx(ipw,49) =     kpgx(ipw,34)*kpgx(ipw, 2)
     kpgx(ipw,50) =     kpgx(ipw,35)*kpgx(ipw, 2)
     kpgx(ipw,51) =     kpgx(ipw,31)*kpgx(ipw, 3)
     kpgx(ipw,52) =     kpgx(ipw,32)*kpgx(ipw, 3)
     kpgx(ipw,53) =     kpgx(ipw,33)*kpgx(ipw, 3)
     kpgx(ipw,54) =     kpgx(ipw,34)*kpgx(ipw, 3)
     kpgx(ipw,55) =     kpgx(ipw,35)*kpgx(ipw, 3)
     kpgx(ipw,56) =     kpgx(ipw,35)*kpgx(ipw, 4)
   end do
 end if
 if (choice==6 .and. nlang>=3) then
!  Define (k+G) part of rank 6 symmetric tensor (28 components)
!  Compressed storage is
!  111111 221111 331111 321111 311111 211111 222111 332111 322111
!  333111 222211 332211 322211 333211 333311 222221 332221 322221
!  333221 333321 333331 222222 332222 322222 333222 333322 333332
!  333333
   do ipw=1,nincpw
     kpgx(ipw,57) =     kpgx(ipw,36)*kpgx(ipw, 2)
     kpgx(ipw,58) =     kpgx(ipw,37)*kpgx(ipw, 2)
     kpgx(ipw,59) =     kpgx(ipw,38)*kpgx(ipw, 2)
     kpgx(ipw,60) =     kpgx(ipw,39)*kpgx(ipw, 2)
     kpgx(ipw,61) =     kpgx(ipw,40)*kpgx(ipw, 2)
     kpgx(ipw,62) =     kpgx(ipw,41)*kpgx(ipw, 2)
     kpgx(ipw,63) =     kpgx(ipw,42)*kpgx(ipw, 2)
     kpgx(ipw,64) =     kpgx(ipw,43)*kpgx(ipw, 2)
     kpgx(ipw,65) =     kpgx(ipw,44)*kpgx(ipw, 2)
     kpgx(ipw,66) =     kpgx(ipw,45)*kpgx(ipw, 2)
     kpgx(ipw,67) =     kpgx(ipw,46)*kpgx(ipw, 2)
     kpgx(ipw,68) =     kpgx(ipw,47)*kpgx(ipw, 2)
     kpgx(ipw,69) =     kpgx(ipw,48)*kpgx(ipw, 2)
     kpgx(ipw,70) =     kpgx(ipw,49)*kpgx(ipw, 2)
     kpgx(ipw,71) =     kpgx(ipw,50)*kpgx(ipw, 2)
     kpgx(ipw,72) =     kpgx(ipw,51)*kpgx(ipw, 2)
     kpgx(ipw,73) =     kpgx(ipw,52)*kpgx(ipw, 2)
     kpgx(ipw,74) =     kpgx(ipw,53)*kpgx(ipw, 2)
     kpgx(ipw,75) =     kpgx(ipw,54)*kpgx(ipw, 2)
     kpgx(ipw,76) =     kpgx(ipw,55)*kpgx(ipw, 2)
     kpgx(ipw,77) =     kpgx(ipw,56)*kpgx(ipw, 2)
     kpgx(ipw,78) =     kpgx(ipw,51)*kpgx(ipw, 3)
     kpgx(ipw,79) =     kpgx(ipw,52)*kpgx(ipw, 3)
     kpgx(ipw,80) =     kpgx(ipw,53)*kpgx(ipw, 3)
     kpgx(ipw,81) =     kpgx(ipw,54)*kpgx(ipw, 3)
     kpgx(ipw,82) =     kpgx(ipw,55)*kpgx(ipw, 3)
     kpgx(ipw,83) =     kpgx(ipw,56)*kpgx(ipw, 3)
     kpgx(ipw,84) =     kpgx(ipw,56)*kpgx(ipw, 4)
   end do
 end if
 if (choice==6 .and. nlang==4) then
!  Define (k+G) part of rank 7 symmetric tensor (36 components)
!  Compressed storage is
!  1111111 2211111 3311111 3211111 3111111 2111111 2221111 3321111 3221111
!  3331111 2222111 3322111 3222111 3332111 3333111 2222211 3322211 3222211
!  3332211 3333211 3333311 2222221 3322221 3222221 3332221 3333221 3333321
!  3333331 2222222 3322222 3222222 3332222 3333222 3333322 3333332 3333333
   do ipw=1,nincpw
     kpgx(ipw,85) =     kpgx(ipw,57)*kpgx(ipw, 2)
     kpgx(ipw,86) =     kpgx(ipw,58)*kpgx(ipw, 2)
     kpgx(ipw,87) =     kpgx(ipw,59)*kpgx(ipw, 2)
     kpgx(ipw,88) =     kpgx(ipw,60)*kpgx(ipw, 2)
     kpgx(ipw,89) =     kpgx(ipw,61)*kpgx(ipw, 2)
     kpgx(ipw,90) =     kpgx(ipw,62)*kpgx(ipw, 2)
     kpgx(ipw,91) =     kpgx(ipw,63)*kpgx(ipw, 2)
     kpgx(ipw,92) =     kpgx(ipw,64)*kpgx(ipw, 2)
     kpgx(ipw,93) =     kpgx(ipw,65)*kpgx(ipw, 2)
     kpgx(ipw,94) =     kpgx(ipw,66)*kpgx(ipw, 2)
     kpgx(ipw,95) =     kpgx(ipw,67)*kpgx(ipw, 2)
     kpgx(ipw,96) =     kpgx(ipw,68)*kpgx(ipw, 2)
     kpgx(ipw,97) =     kpgx(ipw,69)*kpgx(ipw, 2)
     kpgx(ipw,98) =     kpgx(ipw,70)*kpgx(ipw, 2)
     kpgx(ipw,99) =     kpgx(ipw,71)*kpgx(ipw, 2)
     kpgx(ipw,100) =    kpgx(ipw,72)*kpgx(ipw, 2)
     kpgx(ipw,101) =    kpgx(ipw,73)*kpgx(ipw, 2)
     kpgx(ipw,102) =    kpgx(ipw,74)*kpgx(ipw, 2)
     kpgx(ipw,103) =    kpgx(ipw,75)*kpgx(ipw, 2)
     kpgx(ipw,104) =    kpgx(ipw,76)*kpgx(ipw, 2)
     kpgx(ipw,105) =    kpgx(ipw,77)*kpgx(ipw, 2)
     kpgx(ipw,106) =    kpgx(ipw,78)*kpgx(ipw, 2)
     kpgx(ipw,107) =    kpgx(ipw,79)*kpgx(ipw, 2)
     kpgx(ipw,108) =    kpgx(ipw,80)*kpgx(ipw, 2)
     kpgx(ipw,109) =    kpgx(ipw,81)*kpgx(ipw, 2)
     kpgx(ipw,110) =    kpgx(ipw,82)*kpgx(ipw, 2)
     kpgx(ipw,111) =    kpgx(ipw,83)*kpgx(ipw, 2)
     kpgx(ipw,112) =    kpgx(ipw,84)*kpgx(ipw, 2)
     kpgx(ipw,113) =    kpgx(ipw,78)*kpgx(ipw, 3)
     kpgx(ipw,114) =    kpgx(ipw,79)*kpgx(ipw, 3)
     kpgx(ipw,115) =    kpgx(ipw,80)*kpgx(ipw, 3)
     kpgx(ipw,116) =    kpgx(ipw,81)*kpgx(ipw, 3)
     kpgx(ipw,117) =    kpgx(ipw,82)*kpgx(ipw, 3)
     kpgx(ipw,118) =    kpgx(ipw,83)*kpgx(ipw, 3)
     kpgx(ipw,119) =    kpgx(ipw,84)*kpgx(ipw, 3)
     kpgx(ipw,120) =    kpgx(ipw,84)*kpgx(ipw, 4)
   end do
 end if

!*****************************************************************************
!
!Packing of composite projectors in ffkg

 iffkg=0

!Treat composite projectors for the energy
 iln0=0
 do ilmn=1,lmnmax
   iln=indlmn(5,ilmn,itypat)
   if (iln>iln0) then
     iln0=iln
     ilang=1+indlmn(1,ilmn,itypat)
     iproj=indlmn(3,ilmn,itypat)
     if(iproj>0)then
       ilang2=(ilang*(ilang+1))/2

       if(ilang==1)then
!        Treat s-component separately
         ig=ipw1
         iffkg=iffkg+1
         do ipw=1,nincpw
           ffkg(iffkg,ipw)=ffnl(ig,1,ilmn,itypat)
           ig=ig+1
         end do
         parity(iffkg)=2
       else
!        Treat other components (could be made faster by treating explicitely
!        each angular momentum)
         do ii=1,ilang2
!          Get the starting address for the relevant tensor
           jj=ii+((ilang-1)*ilang*(ilang+1))/6
           ig=ipw1
           iffkg=iffkg+1
           do ipw=1,nincpw
             ffkg(iffkg,ipw)=ffnl(ig,1,ilmn,itypat)*kpgx(ipw,jj)
             ig=ig+1
           end do
           if(ilang==2 .or. ilang==4)parity(iffkg)=1
           if(ilang==3)parity(iffkg)=2
         end do
       end if

!      End condition if(iproj>0)
     end if

!    End loop on ilang (ilmn)
   end if
 end do

!This is the number of composite projectors for the energy
 nffkge=iffkg

!Second, treat forces : actually, this part could be rationalized,
!since the outcome is a multiplication by three of the number
!of composite projectors for the energy, while less should be needed
 if((choice==2.or.choice==23) .and. ndgxdt/=1)then
   do ii=1,nffkge
     do ipw=1,nincpw
       ffkg(iffkg+1,ipw)=ffkg(ii,ipw)*kpgx(ipw,2)
       ffkg(iffkg+2,ipw)=ffkg(ii,ipw)*kpgx(ipw,3)
       ffkg(iffkg+3,ipw)=ffkg(ii,ipw)*kpgx(ipw,4)
     end do
     parity(iffkg+1)=3-parity(ii)
     parity(iffkg+2)=parity(iffkg+1)
     parity(iffkg+3)=parity(iffkg+1)
     iffkg=iffkg+3
   end do
 end if
!Note that the additional number of projectors for forces is 3*nffkge

!Third, treat first-derivative of the non-local operator
!with respect to an atomic displacement in one direction :
 if(choice==2 .and. ndgxdt==1)then
   do ii=1,nffkge
     do ipw=1,nincpw
       ffkg(iffkg+1,ipw)=ffkg(ii,ipw)*kpgx(ipw,idir+1)
     end do
     parity(iffkg+1)=3-parity(ii)
     iffkg=iffkg+1
   end do
 end if
!Note that the additional number of projectors for this case is nffkge


!Fourth, treat dynamical matrices : like forces, this part could be rationalized.
 if(choice==4)then
   do ii=1,nffkge
     do ipw=1,nincpw
       kpg_x=kpgx(ipw,2) ; kpg_y=kpgx(ipw,3) ; kpg_z=kpgx(ipw,4)
       ffkg_now=ffkg(ii,ipw)
       ffkg(iffkg+1,ipw)=ffkg_now*kpg_x
       ffkg(iffkg+2,ipw)=ffkg_now*kpg_y
       ffkg(iffkg+3,ipw)=ffkg_now*kpg_z
       ffkg(iffkg+4,ipw)=ffkg_now*kpg_x*kpg_x
       ffkg(iffkg+5,ipw)=ffkg_now*kpg_y*kpg_y
       ffkg(iffkg+6,ipw)=ffkg_now*kpg_z*kpg_z
       ffkg(iffkg+7,ipw)=ffkg_now*kpg_z*kpg_y
       ffkg(iffkg+8,ipw)=ffkg_now*kpg_z*kpg_x
       ffkg(iffkg+9,ipw)=ffkg_now*kpg_y*kpg_x
     end do
     parity(iffkg+1:iffkg+3)=3-parity(ii)
     parity(iffkg+4:iffkg+9)=parity(ii)
     iffkg=iffkg+9
   end do
 end if
!Note that the additional number of projectors for dynamical matrices is 9*nffkge

!Treat composite projectors for the stress or 1st derivative contribution
!to frozen-wavefunction part of elastic tensor
!as well as, for ddk perturbation, the part that depend on ffnl(:,2,..)
 if(choice==3 .or. choice==5 .or. choice==6 .or. choice==23)then

   iln0=0
   do ilmn=1,lmnmax
     if (ispinor==indlmn(6,ilmn,itypat)) then
       iln=indlmn(5,ilmn,itypat)
       if (iln>iln0) then
         iln0=iln
         ilang=1+indlmn(1,ilmn,itypat)
         iproj=indlmn(3,ilmn,itypat)
         if(iproj>0)then
!          number of unique tensor components
           if(choice==3 .or. choice==6 .or. choice==23)ilangx=((ilang+2)*(ilang+3))/2
           if(choice==5)ilangx=(ilang*(ilang+1))/2

           do ii=1,ilangx
!            Get the starting address for the relevant tensor
             if(choice==3 .or. choice==6 .or. choice==23)jj=ii+((ilang+1)*(ilang+2)*(ilang+3))/6
             if(choice==5)jj=ii+((ilang-1)*ilang*(ilang+1))/6
             ig=ipw1
             iffkg=iffkg+1
             if(choice==3 .or. choice==6 .or. choice==23)then
               do ipw=1,nincpw
                 ffkg(iffkg,ipw)=ffnl(ig,2,ilmn,itypat)*kpgx(ipw,jj)
                 ig=ig+1
               end do
             else
               do ipw=1,nincpw
                 ffkg(iffkg,ipw)=ffnl(ig,2,ilmn,itypat)*kpgx(ipw,jj)*&
&                 (kpgx(ipw,2)*gmet(1,idir)+ &
&                 kpgx(ipw,3)*gmet(2,idir)+ &
&                 kpgx(ipw,4)*gmet(3,idir) )
                 ig=ig+1
               end do
             end if
             if(ilang==1 .or. ilang==3)parity(iffkg)=2
             if(ilang==2 .or. ilang==4)parity(iffkg)=1
             if(choice==5)parity(iffkg)=3-parity(iffkg)
           end do

!          End condition if(iproj>0)
         end if
!        End condition on iln
       end if
!      End condition if(ispinor=indlmn(6,...))
     end if
!    End loop on ilmn
   end do

!  End condition of stress
 end if

!Treat composite projectors for the 2nd derivative wrt 2 strains
!and wrt one strain and one atomic displacement (internal strain)
!contributions to frozen-wavefunction part of (generalized) elastic tensor.
!There are 3 sets on terms (in historical order):
!first,  terms with ffnl(:,3,...) and rank+4 tensors.
!second, terms with ffnl(:,1,...) and rank+1 tensors.
!third,  terms with ffnl(:,2,...) and rank+3 tensors.

 if(choice==6)then

   iln0=0
   do ilmn=1,lmnmax
     if (ispinor==indlmn(6,ilmn,itypat)) then
       iln=indlmn(5,ilmn,itypat)
       if (iln>iln0) then
         iln0=iln
         ilang=1+indlmn(1,ilmn,itypat)
         iproj=indlmn(3,ilmn,itypat)

         if(iproj>0)then
!          First set of terms
!          number of unique tensor components
           ilangx=((ilang+4)*(ilang+5))/2

           do ii=1,ilangx
!            Get the starting address for the relevant tensor
             jj=ii+((ilang+3)*(ilang+4)*(ilang+5))/6
             ig=ipw1
             iffkg=iffkg+1
             do ipw=1,nincpw
               ffkg(iffkg,ipw)=ffnl(ig,3,ilmn,itypat)*kpgx(ipw,jj)
               ig=ig+1
             end do
             if(ilang==1 .or. ilang==3)parity(iffkg)=2
             if(ilang==2 .or. ilang==4)parity(iffkg)=1
           end do

!          Second set of terms
!          number of unique tensor components
           ilangx=((ilang+1)*(ilang+2))/2

           do ii=1,ilangx
!            Get the starting address for the relevant tensor
             jj=ii+((ilang)*(ilang+1)*(ilang+2))/6
             ig=ipw1
             iffkg=iffkg+1
             do ipw=1,nincpw
               ffkg(iffkg,ipw)=ffnl(ig,1,ilmn,itypat)*kpgx(ipw,jj)
               ig=ig+1
             end do
             if(ilang==1 .or. ilang==3)parity(iffkg)=1
             if(ilang==2 .or. ilang==4)parity(iffkg)=2
           end do

!          Third set of terms
!          number of unique tensor components
           ilangx=((ilang+3)*(ilang+4))/2

           do ii=1,ilangx
!            Get the starting address for the relevant tensor
             jj=ii+((ilang+2)*(ilang+3)*(ilang+4))/6
             ig=ipw1
             iffkg=iffkg+1
             do ipw=1,nincpw
               ffkg(iffkg,ipw)=ffnl(ig,2,ilmn,itypat)*kpgx(ipw,jj)
               ig=ig+1
             end do
             if(ilang==1 .or. ilang==3)parity(iffkg)=1
             if(ilang==2 .or. ilang==4)parity(iffkg)=2
           end do

!          End condition if(iproj>0)
         end if
!        End condition on iln
       end if
!      End condition if(ispinor=indlmn(6,...))
     end if
!    End loop on ilmn
   end do

!  End condition of 2nd strain derivatives
 end if

!For ddk perturbation, treat the part that depend on ffnl(:,1,..)
!no contribution from s state
 if(nlang>=2 .and. choice==5)then
   iln0=0
   do ilmn=1,lmnmax
     if (ispinor==indlmn(6,ilmn,itypat)) then
       iln=indlmn(5,ilmn,itypat)
       if (iln>iln0) then
         iln0=iln
         ilang=1+indlmn(1,ilmn,itypat)
         if (ilang>=2) then
           iproj=indlmn(3,ilmn,itypat)
           if(iproj>0)then
             ilang2=(ilang*(ilang-1))/2

             do ii=1,ilang2
!              Get the starting address for the relevant tensor
               jj=ii+((ilang-2)*(ilang-1)*ilang)/6
               ig=ipw1
               iffkg=iffkg+1
               do ipw=1,nincpw
                 ffkg(iffkg,ipw)=ffnl(ig,1,ilmn,itypat)*kpgx(ipw,jj)
                 ig=ig+1
               end do
               if(ilang==2 .or. ilang==4)parity(iffkg)=2
               if(ilang==3)parity(iffkg)=1
             end do

!            End condition if(iproj>0)
           end if
!          End condition if(ilang>=2)
         end if
!        End condition if(iln>iln0)
       end if
!      End condition if(ispinor=indlmn(6,...))
     end if
!    End loop on ilmn
   end do
!  End condition of p,d or f state
 end if

!DEBUG
!write(std_out,*)' dfpt_mkffkg : exit '
!ENDDEBUG

end subroutine dfpt_mkffkg
!!***

!!****f* ABINIT/mkffkg
!! NAME
!! mkffkg
!!
!! FUNCTION
!! Prepare the application of the projectors to the shifted wavefunctions,
!! by precomputing the k+G factors and their product with the form factors
!! Do this on a block of plane waves.
!!
!! INPUTS
!!  choice=governs the combination of k+G vectors to be computed
!!  ffnl(npw,nffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  nffnl=3rd dimension of ffnl(2, conventional, or 3 for 2nd derivatives)
!!  idir=direction of the perturbation (needed if choice==2 and ndgxdt==1,
!!       or if choice==5)
!!  indlmn(6,i,ntypat)=array giving l,m,n,lm,ln,spin for i=ln
!!  ipw1 = index of the first plane wave treated in this block
!!  ispinor=1 or 2, gives the spinorial component of ffnl to be used
!!  itypat = type of atom, needed for ffnl
!!  kg_k(3,npw)=integer coords of planewaves in basis sphere
!!  kpg_k(npw,npkg)= (k+G) components and related data
!!  kpt(3)=real components of k point in terms of recip. translations
!!  lmnmax=max. number of (l,n) components over all type of psps
!!  mblkpw=first dimension of kpgx
!!  ndgxdt=number of components of first order derivative
!!  nffkg=number of products of ffnls with combinations of k+G
!!  nincpw=number of plane waves in the block
!!  nkpg=second size of array kpg_k
!!  nlang = number of angular momenta to be treated = 1 + highest ang. mom.
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npw  = total number of plane waves in reciprocal space
!!  ntens=second dimension of kpgx, number of distinct tensorial products
!!  ntypat = number of type of atoms, dimension needed for ffnl
!!
!! OUTPUT
!!  kpgx(mblkpw,ntens)=different tensorial products of k+G
!!  ffkg(mblkpw,nffkg)=different products of ffnls with k+G
!!  parity(nffkg)=parity of the tensorial product of k+G (2 if even, 1 of odd)
!!
!! NOTES
!!  This routine must be thread-safe as it is called inside loops that are OpenMP parallelized.
!!  Please, do not add variables with the save attribute or SIDE EFFECTS.
!!
!! PARENTS
!!      opernl2
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkffkg(choice,ffkg,ffnl,gmet,idir,indlmn,ipw1,ispinor,itypat,&
&                  kg_k,kpg_k,kpgx,kpt,lmnmax,mblkpw,ndgxdt,nffkg,nffnl,nincpw,nkpg,nlang,&
&                  npw,ntens,ntypat,parity)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,idir,ipw1,ispinor,itypat,lmnmax,mblkpw,ndgxdt
 integer,intent(in) :: nffkg,nffnl,nincpw,nkpg,nlang,npw,ntens,ntypat
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),kg_k(3,npw)
 integer,intent(out) :: parity(nffkg)
 real(dp),intent(in) :: ffnl(npw,nffnl,lmnmax,ntypat),gmet(3,3),kpg_k(npw,nkpg)
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(out) :: ffkg(mblkpw,nffkg),kpgx(mblkpw,ntens)

!Local variables-------------------------------
!scalars
 integer :: iffkg,ig,ii,ilang,ilang2,ilangx,ilmn,iln,iln0,iproj,ipw,jj
 integer :: nffkge
 real(dp) :: ffkg_now,kpg_x,kpg_y,kpg_z
!arrays

! *************************************************************************

 jj=0;ilangx=0

!This will be useless after all the modifications have been done
 do ipw=1,nincpw
   kpgx(ipw,1)=1.0d0
 end do

!Initialize kpgx array related to tensors defined below
 if ( nlang>=2 .or. choice==2 .or. choice==3 .or. choice==4 .or. choice==5&
& .or. choice==6 .or. choice==23) then
   if (nkpg>=3) then
     kpgx(1:nincpw,2)=kpg_k(ipw1+1:ipw1+nincpw,1)
     kpgx(1:nincpw,3)=kpg_k(ipw1+1:ipw1+nincpw,2)
     kpgx(1:nincpw,4)=kpg_k(ipw1+1:ipw1+nincpw,3)
   else
     ig=ipw1
     do ipw=1,nincpw
       kpgx(ipw,2)=kpt(1)+dble(kg_k(1,ig))
       kpgx(ipw,3)=kpt(2)+dble(kg_k(2,ig))
       kpgx(ipw,4)=kpt(3)+dble(kg_k(3,ig))
       ig=ig+1
     end do
   end if
 end if
 if (nlang>=3 .or. choice==3 .or. choice==6 .or. choice==23) then
!  Define (k+G) part of rank 2 symmetric tensor (6 components), l=2
!  Compressed storage is 11 22 33 32 31 21
   if (nkpg>=9) then
     kpgx(1:nincpw,5) =kpg_k(ipw1+1:ipw1+nincpw,4)
     kpgx(1:nincpw,6) =kpg_k(ipw1+1:ipw1+nincpw,5)
     kpgx(1:nincpw,7) =kpg_k(ipw1+1:ipw1+nincpw,6)
     kpgx(1:nincpw,8) =kpg_k(ipw1+1:ipw1+nincpw,7)
     kpgx(1:nincpw,9) =kpg_k(ipw1+1:ipw1+nincpw,8)
     kpgx(1:nincpw,10)=kpg_k(ipw1+1:ipw1+nincpw,9)
   else
     do ipw=1,nincpw
       kpgx(ipw, 5) =      kpgx(ipw, 2)*kpgx(ipw, 2)
       kpgx(ipw, 6) =      kpgx(ipw, 3)*kpgx(ipw, 3)
       kpgx(ipw, 7) =      kpgx(ipw, 4)*kpgx(ipw, 4)
       kpgx(ipw, 8) =      kpgx(ipw, 4)*kpgx(ipw, 3)
       kpgx(ipw, 9) =      kpgx(ipw, 4)*kpgx(ipw, 2)
       kpgx(ipw,10) =      kpgx(ipw, 3)*kpgx(ipw, 2)
     end do
   end if
 end if
 if (nlang>=4 .or. ((choice==3.or.choice==23) .and. nlang>=2) .or. choice==6) then
!  Define (k+G) part of rank 3 symmetric tensor (10 components), l=3
!  Compressed storage is 111 221 331 321 311 211 222 332 322 333
   do ipw=1,nincpw
     kpgx(ipw,11) =     kpgx(ipw, 5)*kpgx(ipw, 2)
     kpgx(ipw,12) =     kpgx(ipw, 6)*kpgx(ipw, 2)
     kpgx(ipw,13) =     kpgx(ipw, 7)*kpgx(ipw, 2)
     kpgx(ipw,14) =     kpgx(ipw, 8)*kpgx(ipw, 2)
     kpgx(ipw,15) =     kpgx(ipw, 9)*kpgx(ipw, 2)
     kpgx(ipw,16) =     kpgx(ipw,10)*kpgx(ipw, 2)
     kpgx(ipw,17) =     kpgx(ipw, 6)*kpgx(ipw, 3)
     kpgx(ipw,18) =     kpgx(ipw, 7)*kpgx(ipw, 3)
     kpgx(ipw,19) =     kpgx(ipw, 8)*kpgx(ipw, 3)
     kpgx(ipw,20) =     kpgx(ipw, 7)*kpgx(ipw, 4)
   end do
 end if
 if (((choice==3.or.choice==23) .and. nlang>=3) .or. choice==6) then
!  Add additional tensors for strain gradients
!  Define (k+G) part of rank 4 symmetric tensor (15 components), l=2
!  Compressed storage is 1111 2211 3311 3211 3111 2111 2221 3321 3221
!  3331 2222 3322 3222 3332 3333
   do ipw=1,nincpw
     kpgx(ipw,21) =     kpgx(ipw, 5)*kpgx(ipw, 5)
     kpgx(ipw,22) =     kpgx(ipw, 6)*kpgx(ipw, 5)
     kpgx(ipw,23) =     kpgx(ipw, 7)*kpgx(ipw, 5)
     kpgx(ipw,24) =     kpgx(ipw, 8)*kpgx(ipw, 5)
     kpgx(ipw,25) =     kpgx(ipw, 9)*kpgx(ipw, 5)
     kpgx(ipw,26) =     kpgx(ipw,10)*kpgx(ipw, 5)
     kpgx(ipw,27) =     kpgx(ipw, 6)*kpgx(ipw,10)
     kpgx(ipw,28) =     kpgx(ipw, 7)*kpgx(ipw,10)
     kpgx(ipw,29) =     kpgx(ipw, 8)*kpgx(ipw,10)
     kpgx(ipw,30) =     kpgx(ipw, 7)*kpgx(ipw, 9)
     kpgx(ipw,31) =     kpgx(ipw, 6)*kpgx(ipw, 6)
     kpgx(ipw,32) =     kpgx(ipw, 7)*kpgx(ipw, 6)
     kpgx(ipw,33) =     kpgx(ipw, 8)*kpgx(ipw, 6)
     kpgx(ipw,34) =     kpgx(ipw, 7)*kpgx(ipw, 8)
     kpgx(ipw,35) =     kpgx(ipw, 7)*kpgx(ipw, 7)
   end do
 end if
 if (((choice==3.or.choice==23) .and. nlang>=4) .or. (choice==6 .and. nlang>=2)) then
!  Define (k+G) part of rank 5 symmetric tensor (21 components), l=3
!  Compressed storage is 11111 22111 33111 32111 31111 21111
!  22211 33211 32211     33311 22221 33221 32221 33321 33331
!  22222 33222 32222     33322 33332 33333
   do ipw=1,nincpw
     kpgx(ipw,36) =     kpgx(ipw,21)*kpgx(ipw, 2)
     kpgx(ipw,37) =     kpgx(ipw,22)*kpgx(ipw, 2)
     kpgx(ipw,38) =     kpgx(ipw,23)*kpgx(ipw, 2)
     kpgx(ipw,39) =     kpgx(ipw,24)*kpgx(ipw, 2)
     kpgx(ipw,40) =     kpgx(ipw,25)*kpgx(ipw, 2)
     kpgx(ipw,41) =     kpgx(ipw,26)*kpgx(ipw, 2)
     kpgx(ipw,42) =     kpgx(ipw,27)*kpgx(ipw, 2)
     kpgx(ipw,43) =     kpgx(ipw,28)*kpgx(ipw, 2)
     kpgx(ipw,44) =     kpgx(ipw,29)*kpgx(ipw, 2)
     kpgx(ipw,45) =     kpgx(ipw,30)*kpgx(ipw, 2)
     kpgx(ipw,46) =     kpgx(ipw,31)*kpgx(ipw, 2)
     kpgx(ipw,47) =     kpgx(ipw,32)*kpgx(ipw, 2)
     kpgx(ipw,48) =     kpgx(ipw,33)*kpgx(ipw, 2)
     kpgx(ipw,49) =     kpgx(ipw,34)*kpgx(ipw, 2)
     kpgx(ipw,50) =     kpgx(ipw,35)*kpgx(ipw, 2)
     kpgx(ipw,51) =     kpgx(ipw,31)*kpgx(ipw, 3)
     kpgx(ipw,52) =     kpgx(ipw,32)*kpgx(ipw, 3)
     kpgx(ipw,53) =     kpgx(ipw,33)*kpgx(ipw, 3)
     kpgx(ipw,54) =     kpgx(ipw,34)*kpgx(ipw, 3)
     kpgx(ipw,55) =     kpgx(ipw,35)*kpgx(ipw, 3)
     kpgx(ipw,56) =     kpgx(ipw,35)*kpgx(ipw, 4)
   end do
 end if
 if (choice==6 .and. nlang>=3) then
!  Define (k+G) part of rank 6 symmetric tensor (28 components)
!  Compressed storage is
!  111111 221111 331111 321111 311111 211111 222111 332111 322111
!  333111 222211 332211 322211 333211 333311 222221 332221 322221
!  333221 333321 333331 222222 332222 322222 333222 333322 333332
!  333333
   do ipw=1,nincpw
     kpgx(ipw,57) =     kpgx(ipw,36)*kpgx(ipw, 2)
     kpgx(ipw,58) =     kpgx(ipw,37)*kpgx(ipw, 2)
     kpgx(ipw,59) =     kpgx(ipw,38)*kpgx(ipw, 2)
     kpgx(ipw,60) =     kpgx(ipw,39)*kpgx(ipw, 2)
     kpgx(ipw,61) =     kpgx(ipw,40)*kpgx(ipw, 2)
     kpgx(ipw,62) =     kpgx(ipw,41)*kpgx(ipw, 2)
     kpgx(ipw,63) =     kpgx(ipw,42)*kpgx(ipw, 2)
     kpgx(ipw,64) =     kpgx(ipw,43)*kpgx(ipw, 2)
     kpgx(ipw,65) =     kpgx(ipw,44)*kpgx(ipw, 2)
     kpgx(ipw,66) =     kpgx(ipw,45)*kpgx(ipw, 2)
     kpgx(ipw,67) =     kpgx(ipw,46)*kpgx(ipw, 2)
     kpgx(ipw,68) =     kpgx(ipw,47)*kpgx(ipw, 2)
     kpgx(ipw,69) =     kpgx(ipw,48)*kpgx(ipw, 2)
     kpgx(ipw,70) =     kpgx(ipw,49)*kpgx(ipw, 2)
     kpgx(ipw,71) =     kpgx(ipw,50)*kpgx(ipw, 2)
     kpgx(ipw,72) =     kpgx(ipw,51)*kpgx(ipw, 2)
     kpgx(ipw,73) =     kpgx(ipw,52)*kpgx(ipw, 2)
     kpgx(ipw,74) =     kpgx(ipw,53)*kpgx(ipw, 2)
     kpgx(ipw,75) =     kpgx(ipw,54)*kpgx(ipw, 2)
     kpgx(ipw,76) =     kpgx(ipw,55)*kpgx(ipw, 2)
     kpgx(ipw,77) =     kpgx(ipw,56)*kpgx(ipw, 2)
     kpgx(ipw,78) =     kpgx(ipw,51)*kpgx(ipw, 3)
     kpgx(ipw,79) =     kpgx(ipw,52)*kpgx(ipw, 3)
     kpgx(ipw,80) =     kpgx(ipw,53)*kpgx(ipw, 3)
     kpgx(ipw,81) =     kpgx(ipw,54)*kpgx(ipw, 3)
     kpgx(ipw,82) =     kpgx(ipw,55)*kpgx(ipw, 3)
     kpgx(ipw,83) =     kpgx(ipw,56)*kpgx(ipw, 3)
     kpgx(ipw,84) =     kpgx(ipw,56)*kpgx(ipw, 4)
   end do
 end if
 if (choice==6 .and. nlang==4) then
!  Define (k+G) part of rank 7 symmetric tensor (36 components)
!  Compressed storage is
!  1111111 2211111 3311111 3211111 3111111 2111111 2221111 3321111 3221111
!  3331111 2222111 3322111 3222111 3332111 3333111 2222211 3322211 3222211
!  3332211 3333211 3333311 2222221 3322221 3222221 3332221 3333221 3333321
!  3333331 2222222 3322222 3222222 3332222 3333222 3333322 3333332 3333333
   do ipw=1,nincpw
     kpgx(ipw,85) =     kpgx(ipw,57)*kpgx(ipw, 2)
     kpgx(ipw,86) =     kpgx(ipw,58)*kpgx(ipw, 2)
     kpgx(ipw,87) =     kpgx(ipw,59)*kpgx(ipw, 2)
     kpgx(ipw,88) =     kpgx(ipw,60)*kpgx(ipw, 2)
     kpgx(ipw,89) =     kpgx(ipw,61)*kpgx(ipw, 2)
     kpgx(ipw,90) =     kpgx(ipw,62)*kpgx(ipw, 2)
     kpgx(ipw,91) =     kpgx(ipw,63)*kpgx(ipw, 2)
     kpgx(ipw,92) =     kpgx(ipw,64)*kpgx(ipw, 2)
     kpgx(ipw,93) =     kpgx(ipw,65)*kpgx(ipw, 2)
     kpgx(ipw,94) =     kpgx(ipw,66)*kpgx(ipw, 2)
     kpgx(ipw,95) =     kpgx(ipw,67)*kpgx(ipw, 2)
     kpgx(ipw,96) =     kpgx(ipw,68)*kpgx(ipw, 2)
     kpgx(ipw,97) =     kpgx(ipw,69)*kpgx(ipw, 2)
     kpgx(ipw,98) =     kpgx(ipw,70)*kpgx(ipw, 2)
     kpgx(ipw,99) =     kpgx(ipw,71)*kpgx(ipw, 2)
     kpgx(ipw,100) =    kpgx(ipw,72)*kpgx(ipw, 2)
     kpgx(ipw,101) =    kpgx(ipw,73)*kpgx(ipw, 2)
     kpgx(ipw,102) =    kpgx(ipw,74)*kpgx(ipw, 2)
     kpgx(ipw,103) =    kpgx(ipw,75)*kpgx(ipw, 2)
     kpgx(ipw,104) =    kpgx(ipw,76)*kpgx(ipw, 2)
     kpgx(ipw,105) =    kpgx(ipw,77)*kpgx(ipw, 2)
     kpgx(ipw,106) =    kpgx(ipw,78)*kpgx(ipw, 2)
     kpgx(ipw,107) =    kpgx(ipw,79)*kpgx(ipw, 2)
     kpgx(ipw,108) =    kpgx(ipw,80)*kpgx(ipw, 2)
     kpgx(ipw,109) =    kpgx(ipw,81)*kpgx(ipw, 2)
     kpgx(ipw,110) =    kpgx(ipw,82)*kpgx(ipw, 2)
     kpgx(ipw,111) =    kpgx(ipw,83)*kpgx(ipw, 2)
     kpgx(ipw,112) =    kpgx(ipw,84)*kpgx(ipw, 2)
     kpgx(ipw,113) =    kpgx(ipw,78)*kpgx(ipw, 3)
     kpgx(ipw,114) =    kpgx(ipw,79)*kpgx(ipw, 3)
     kpgx(ipw,115) =    kpgx(ipw,80)*kpgx(ipw, 3)
     kpgx(ipw,116) =    kpgx(ipw,81)*kpgx(ipw, 3)
     kpgx(ipw,117) =    kpgx(ipw,82)*kpgx(ipw, 3)
     kpgx(ipw,118) =    kpgx(ipw,83)*kpgx(ipw, 3)
     kpgx(ipw,119) =    kpgx(ipw,84)*kpgx(ipw, 3)
     kpgx(ipw,120) =    kpgx(ipw,84)*kpgx(ipw, 4)
   end do
 end if

!*****************************************************************************
!
!Packing of composite projectors in ffkg

 iffkg=0

!Treat composite projectors for the energy
 iln0=0
 do ilmn=1,lmnmax
   iln=indlmn(5,ilmn,itypat)
   if (iln>iln0) then
     iln0=iln
     ilang=1+indlmn(1,ilmn,itypat)
     iproj=indlmn(3,ilmn,itypat)
     if(iproj>0)then
       ilang2=(ilang*(ilang+1))/2

       if(ilang==1)then
!        Treat s-component separately
         ig=ipw1
         iffkg=iffkg+1
         do ipw=1,nincpw
           ffkg(ipw,iffkg)=ffnl(ig,1,ilmn,itypat)
           ig=ig+1
         end do
         parity(iffkg)=2
       else
!        Treat other components (could be made faster by treating explicitely
!        each angular momentum)
         do ii=1,ilang2
!          Get the starting address for the relevant tensor
           jj=ii+((ilang-1)*ilang*(ilang+1))/6
           ig=ipw1
           iffkg=iffkg+1
           do ipw=1,nincpw
             ffkg(ipw,iffkg)=ffnl(ig,1,ilmn,itypat)*kpgx(ipw,jj)
             ig=ig+1
           end do
           if(ilang==2 .or. ilang==4)parity(iffkg)=1
           if(ilang==3)parity(iffkg)=2
         end do
       end if

!      End condition if(iproj>0)
     end if
!    End loop on ilang (ilmn)
   end if
 end do

!This is the number of composite projectors for the energy
 nffkge=iffkg

!Second, treat forces : actually, this part could be rationalized,
!since the outcome is a multiplication by three of the number
!of composite projectors for the energy, while less should be needed
 if((choice==2.or.choice==23) .and. ndgxdt/=1)then
   do ii=1,nffkge
     do ipw=1,nincpw
       ffkg(ipw,iffkg+1)=ffkg(ipw,ii)*kpgx(ipw,2)
       ffkg(ipw,iffkg+2)=ffkg(ipw,ii)*kpgx(ipw,3)
       ffkg(ipw,iffkg+3)=ffkg(ipw,ii)*kpgx(ipw,4)
     end do
     parity(iffkg+1)=3-parity(ii)
     parity(iffkg+2)=parity(iffkg+1)
     parity(iffkg+3)=parity(iffkg+1)
     iffkg=iffkg+3
   end do
 end if
!Note that the additional number of projectors for forces is 3*nffkge

!Third, treat first-derivative of the non-local operator
!with respect to an atomic displacement in one direction :
 if(choice==2 .and. ndgxdt==1)then
   do ii=1,nffkge
     do ipw=1,nincpw
       ffkg(ipw,iffkg+1)=ffkg(ipw,ii)*kpgx(ipw,idir+1)
     end do
     parity(iffkg+1)=3-parity(ii)
     iffkg=iffkg+1
   end do
 end if
!Note that the additional number of projectors for this case is nffkge


!Fourth, treat dynamical matrices : like forces, this part could be rationalized.
 if(choice==4)then
   do ii=1,nffkge
     do ipw=1,nincpw
       kpg_x=kpgx(ipw,2) ; kpg_y=kpgx(ipw,3) ; kpg_z=kpgx(ipw,4)
       ffkg_now=ffkg(ipw,ii)
       ffkg(ipw,iffkg+1)=ffkg_now*kpg_x
       ffkg(ipw,iffkg+2)=ffkg_now*kpg_y
       ffkg(ipw,iffkg+3)=ffkg_now*kpg_z
       ffkg(ipw,iffkg+4)=ffkg_now*kpg_x*kpg_x
       ffkg(ipw,iffkg+5)=ffkg_now*kpg_y*kpg_y
       ffkg(ipw,iffkg+6)=ffkg_now*kpg_z*kpg_z
       ffkg(ipw,iffkg+7)=ffkg_now*kpg_z*kpg_y
       ffkg(ipw,iffkg+8)=ffkg_now*kpg_z*kpg_x
       ffkg(ipw,iffkg+9)=ffkg_now*kpg_y*kpg_x
     end do
     parity(iffkg+1:iffkg+3)=3-parity(ii)
     parity(iffkg+4:iffkg+9)=parity(ii)
     iffkg=iffkg+9
   end do
 end if
!Note that the additional number of projectors for dynamical matrices is 9*nffkge

!Treat composite projectors for the stress or 1st derivative contribution
!to frozen-wavefunction part of elastic tensor
!as well as, for ddk perturbation, the part that depend on ffnl(:,2,..)
 if(choice==3 .or. choice==5 .or. choice==6 .or. choice==23)then

   iln0=0
   do ilmn=1,lmnmax
     if (ispinor==indlmn(6,ilmn,itypat)) then
       iln=indlmn(5,ilmn,itypat)
       if (iln>iln0) then
         iln0=iln
         ilang=1+indlmn(1,ilmn,itypat)
         iproj=indlmn(3,ilmn,itypat)
         if(iproj>0)then
!          number of unique tensor components
           if(choice==3 .or. choice==6 .or. choice==23)ilangx=((ilang+2)*(ilang+3))/2
           if(choice==5)ilangx=(ilang*(ilang+1))/2

           do ii=1,ilangx
!            Get the starting address for the relevant tensor
             if(choice==3 .or. choice==6 .or. choice==23)jj=ii+((ilang+1)*(ilang+2)*(ilang+3))/6
             if(choice==5)jj=ii+((ilang-1)*ilang*(ilang+1))/6
             ig=ipw1
             iffkg=iffkg+1
             if(choice==3 .or. choice==6 .or. choice==23)then
               do ipw=1,nincpw
                 ffkg(ipw,iffkg)=ffnl(ig,2,ilmn,itypat)*kpgx(ipw,jj)
                 ig=ig+1
               end do
             else
               do ipw=1,nincpw
                 ffkg(ipw,iffkg)=ffnl(ig,2,ilmn,itypat)*kpgx(ipw,jj)*&
&                 (kpgx(ipw,2)*gmet(1,idir)+ &
&                 kpgx(ipw,3)*gmet(2,idir)+ &
&                 kpgx(ipw,4)*gmet(3,idir) )
                 ig=ig+1
               end do
             end if
             if(ilang==1 .or. ilang==3)parity(iffkg)=2
             if(ilang==2 .or. ilang==4)parity(iffkg)=1
             if(choice==5)parity(iffkg)=3-parity(iffkg)
           end do

!          End condition if(iproj>0)
         end if

!        End loop on ilang (ilmn)
       end if
     end if
   end do

!  End condition of stress
 end if

!Treat composite projectors for the 2nd derivative wrt 2 strains
!and wrt one strain and one atomic displacement (internal strain)
!contributions to frozen-wavefunction part of (generalized) elastic tensor.
!There are 3 sets on terms (in historical order):
!first,  terms with ffnl(:,3,...) and rank+4 tensors.
!second, terms with ffnl(:,1,...) and rank+1 tensors.
!third,  terms with ffnl(:,2,...) and rank+3 tensors.

 if(choice==6)then

   iln0=0
   do ilmn=1,lmnmax
     if (ispinor==indlmn(6,ilmn,itypat)) then
       iln=indlmn(5,ilmn,itypat)
       if (iln>iln0) then
         iln0=iln
         ilang=1+indlmn(1,ilmn,itypat)
         iproj=indlmn(3,ilmn,itypat)
         if(iproj>0)then
!          First set of terms
!          number of unique tensor components
           ilangx=((ilang+4)*(ilang+5))/2

           do ii=1,ilangx
!            Get the starting address for the relevant tensor
             jj=ii+((ilang+3)*(ilang+4)*(ilang+5))/6
             ig=ipw1
             iffkg=iffkg+1
             do ipw=1,nincpw
               ffkg(ipw,iffkg)=ffnl(ig,3,ilmn,itypat)*kpgx(ipw,jj)
               ig=ig+1
             end do
             if(ilang==1 .or. ilang==3)parity(iffkg)=2
             if(ilang==2 .or. ilang==4)parity(iffkg)=1
           end do

!          Second set of terms
!          number of unique tensor components
           ilangx=((ilang+1)*(ilang+2))/2

           do ii=1,ilangx
!            Get the starting address for the relevant tensor
             jj=ii+((ilang)*(ilang+1)*(ilang+2))/6
             ig=ipw1
             iffkg=iffkg+1
             do ipw=1,nincpw
               ffkg(ipw,iffkg)=ffnl(ig,1,ilmn,itypat)*kpgx(ipw,jj)
               ig=ig+1
             end do
             if(ilang==1 .or. ilang==3)parity(iffkg)=1
             if(ilang==2 .or. ilang==4)parity(iffkg)=2
           end do

!          Third set of terms
!          number of unique tensor components
           ilangx=((ilang+3)*(ilang+4))/2

           do ii=1,ilangx
!            Get the starting address for the relevant tensor
             jj=ii+((ilang+2)*(ilang+3)*(ilang+4))/6
             ig=ipw1
             iffkg=iffkg+1
             do ipw=1,nincpw
               ffkg(ipw,iffkg)=ffnl(ig,2,ilmn,itypat)*kpgx(ipw,jj)
               ig=ig+1
             end do
             if(ilang==1 .or. ilang==3)parity(iffkg)=1
             if(ilang==2 .or. ilang==4)parity(iffkg)=2
           end do

!          End condition if(iproj>0)
         end if
!        End loop on ilang (ilmn)
       end if
     end if
   end do

!  End condition of 2nd strain derivatives
 end if

!For ddk perturbation, treat the part that depend on ffnl(:,1,..)
!no contribution from s state
 if(nlang>=2 .and. choice==5)then
   iln0=0
   do ilmn=1,lmnmax
     if (ispinor==indlmn(6,ilmn,itypat)) then
       iln=indlmn(5,ilmn,itypat)
       if (iln>iln0) then
         iln0=iln
         ilang=1+indlmn(1,ilmn,itypat)
         if (ilang>=2) then
           iproj=indlmn(3,ilmn,itypat)
           if(iproj>0)then
             ilang2=(ilang*(ilang-1))/2

             do ii=1,ilang2
!              Get the starting address for the relevant tensor
               jj=ii+((ilang-2)*(ilang-1)*ilang)/6
               ig=ipw1
               iffkg=iffkg+1
               do ipw=1,nincpw
                 ffkg(ipw,iffkg)=ffnl(ig,1,ilmn,itypat)*kpgx(ipw,jj)
                 ig=ig+1
               end do
               if(ilang==2 .or. ilang==4)parity(iffkg)=2
               if(ilang==3)parity(iffkg)=1
             end do

!            End condition if(iproj>0)
           end if
!          End loop on ilang>=2
         end if
       end if
     end if
   end do
!  End condition of p,d or f state
 end if

end subroutine mkffkg
!!***

end module m_mkffkg
!!***
