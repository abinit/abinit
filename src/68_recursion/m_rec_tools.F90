!!****m* ABINIT/m_rec_tools
!! NAME
!!  m_rec_tools
!!
!! FUNCTION
!!  This module provides some functions more or less generic used
!!  in the Recursion Mathod
!!
!! COPYRIGHT
!! Copyright (C) 2002-2020 ABINIT group (MMancini)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! NOTES
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

module m_rec_tools

 use defs_basis
 use m_abicore

 use defs_rectypes, only : recparall_type

 implicit none

 private

 public ::            &
   get_pt0_pt1,       &  !--To get pt0 pt1 from inf,sup
   reshape_pot,       &  !--To rescale the potential between 2 grids
   trottersum            !--To calculate the trotter sum

CONTAINS  !===========================================================
!!***


!!****f* m_rec_tools/get_pt0_pt1
!! NAME
!! get_pt0_pt1
!!
!! FUNCTION
!! utility function to get pt0 and pt1 for given inf and sup
!!
!! INPUTS
!! ngfft(3) = fine grid (corresponds to dtset%ngfft(1:3))
!! inf = inferior point in-line coordinate on the coarse grid
!! sup = superior point in-line coordinate on the coarse grid
!!
!! OUTPUT
!! recpar%pt0<type(vec_int)>=Intial point for this proc in x,y,z
!! recpar%pt1<type(vec_int)>=Final point for this proc in x,y,z
!! recpar%min_pt=inferior point in-line coordinate on the fine grid
!! recpar%max_pt=superior point in-line coordinate on the dine grid
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      first_rec,m_rec
!!
!! CHILDREN
!!
!! SOURCE
subroutine get_pt0_pt1(ngfft,gratio,inf,sup,recpar)

 implicit none

!Arguments ------------------------------------
  integer,intent(in)        :: gratio
  integer,intent(in)        :: inf,sup
  type(recparall_type),intent(inout) :: recpar
  integer,intent(in)        :: ngfft(3)
!Local ---------------------------
  integer :: x,y,z,count
  integer :: boz,boy,pt
! *********************************************************************
  count = 0
  do z = 0,ngfft(3)-1,gratio
    boz = z*ngfft(2)
    do y = 0,ngfft(2)-1,gratio
      boy = (boz+y)*ngfft(1)
      do x = 0,ngfft(1)-1,gratio
        pt = boy+x
        if(count >= inf) then
          if(count == inf) then
            recpar%pt0%x = x; recpar%pt0%y = y; recpar%pt0%z = z
            recpar%min_pt = pt
          end if
          recpar%pt1%x = x; recpar%pt1%y = y; recpar%pt1%z = z
          recpar%max_pt = pt
        endif
        count = count+1
        if(count == sup ) return
      end do
    end do
  end do
end subroutine get_pt0_pt1
!!***

!!****f* m_rec_tools/reshape_pot
!! NAME
!! reshape_pot
!!
!! FUNCTION
!! Reshape array on
!!
!! INPUTS
!! nfft=size of the big grid
!! nfftrec=size of the cut grid
!! trasl(3) center point
!! ngfft(3) dimensions of the big grid
!! ngfftrec(3) dimensions of the cut grid
!! pot(0:nfft-1) 3d-array on the big grid
!!
!! OUTPUT
!! potloc is a cut of pot around trasl
!!
!! PARENTS
!!      first_rec,nlenergyrec,vtorhorec
!!
!! CHILDREN
!!
!! SOURCE

subroutine reshape_pot(trasl,nfft,nfftrec,ngfft,ngfftrec,pot,potloc)

 implicit none

 !Arguments ------------------------------------
 integer,  intent(in) :: nfft,nfftrec
 integer,  intent(in) :: trasl(3)
 integer,  intent(in) :: ngfftrec(3),ngfft(3)
 real(dp), intent(in) :: pot(0:nfft-1)
 real(dp), intent(out):: potloc(0:nfftrec-1)
 !Local ----------------------------------------
 ! scalars
 integer :: zz,yy,xx,parz,pary
 integer :: modi,modj,modk
 !character(len=500) :: msg
 ! *********************************************************************
 do zz = 0,ngfftrec(3)-1
   modk = modulo(zz-trasl(3),ngfft(3))*ngfft(2)
   parz = ngfftrec(2)*zz
   do yy = 0,ngfftrec(2)-1
     modj = (modulo(yy-trasl(2),ngfft(2))+modk)*ngfft(1)
     pary = ngfftrec(1)*(yy+parz)
     do xx = 0,ngfftrec(1)-1
       modi = modulo(xx-trasl(1),ngfft(1))+modj
       potloc(xx+pary) = pot(modi)
     end do
   end do
 end do

end subroutine reshape_pot
!!***

!!****f* m_rec_tools/trottersum
!! NAME
!! trottersum
!!
!! FUNCTION
!! Calculate the contribution to the partial fraction decomposition
!! due to a recursion step.
!!
!! INPUTS
!!  dim_trott=dimension of the partial fraction decomposition (PFD)
!!  pi_on_rtrotter=parameter pi/rtrotter
!!  an,bn2=recursion coefficients at irec
!!  exp1=numerical factor  (see density_rec)
!!  coeef_mu=numerical factor
!!
!! OUTPUT
!!
!! SIZE EFFECTS
!!  D,N=denominator and numerator accumalator of PFD
!!  Dold,Nold=denominator and numerator of PFD (old values)
!!  facrec0=used to select irec=0
!!  error=estimated error of recursion at this step
!!  prod_b2=numerical factor
!!
!! PARENTS
!!      density_rec,recursion,recursion_nl
!!
!! CHILDREN
!!
!! SOURCE

subroutine trottersum(dim_trott,error,&
     &                prod_b2,pi_on_rtrotter,&
     &                facrec0,coeef_mu,exp1,&
     &                an,bn2,&
     &                N,D,Nold,Dold)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,  intent(in) :: dim_trott
 real(dp), intent(in) :: an,bn2,exp1,pi_on_rtrotter
 real(dp), intent(inout) :: error,prod_b2
 complex(dp), intent(in) :: coeef_mu
 complex(dp), intent(inout) :: facrec0
 !arrays
 complex(dpc),intent(inout) :: D(0:dim_trott),Dold(0:dim_trott)
 complex(dpc),intent(inout) :: N(0:dim_trott),Nold(0:dim_trott)
 !Local ----------------------------------------
 ! scalars
 integer :: itrot
 real(dp) :: arg
 complex(dpc) :: Dnew,Nnew,zj
 !character(len=500) :: msg
 ! *********************************************************************

 error = zero
 prod_b2 = prod_b2 * exp1 * bn2

 do itrot=0,dim_trott
   arg = pi_on_rtrotter*(real( itrot,dp) + half )
   zj = cmplx(cos(arg),sin(arg),dp)*coeef_mu

   Nnew = zj*facrec0 + &
        & (zj - cmplx(an  ,zero,dp))*N(itrot) - &
        &       cmplx(bn2 ,zero,dp)*Nold(itrot)
   Dnew = (zj - cmplx(an  ,zero,dp))*D(itrot) - &
        &       cmplx(bn2 ,zero,dp)*Dold(itrot)
   Nold(itrot) = N(itrot)
   Dold(itrot) = D(itrot)
   N(itrot) = Nnew
   D(itrot) = Dnew

   !--Error estimator
   error = error + abs(prod_b2/(D(itrot)*Dold(itrot)))
 end do
 facrec0 = czero

end subroutine trottersum
!!***

end module m_rec_tools
!!***
