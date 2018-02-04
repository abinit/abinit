!{\src2tex{textfont=tt}}
!!****f* ABINIT/pspnl_operat_rec
!! NAME
!! pspnl_operat_rec
!!
!! FUNCTION
!! It calculates the non-local projectors used in recursion for any
!! psp non-local:
!! The nl interaction in recursion is $$exp{-V_{NL}/beta}=\sum_A\sum_{lm}
!! \sum{ij}Y_{lm}(\hat{r-R_A}')f^l_i(r-R_A)D^l_{i,j}Y_{lm}(\hat{r-R_A})f^l_j{r-R_A}$$
!! where $D^_{i,j}$ is a matrix  previously (see pspnl_operat_rec).
!! In this routine  the projectors $Y_{lm}(\hat{r-R_A}')f^l_i(r-R_A)$
!! are calculated. So an array of dimensions
!! rset%nl%projec(nfftrec,lmnmax,nlrec%npsp)
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2018 ABINIT group (the_author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! metrec<metricrec_type>=contains information concerning metric in
!!         recursion: grid_step, metric, infinitesimal volume
!! ngfftrec(18)=is the ngfft grid (truncated if different from ngfft) of recursion 
!! debug=debug variable
!!
!!
!! OUTPUT
!! nlrec<nlpsprec_type>%projec= array containig the projectors on the real grid
!! nlrec<nlpsprec_type>%intlen= integer linear size of the non-local grid
!!
!! SIDE EFFECTS
!! nlrec<nlpsprec_type> data set of non-local pseudo for recursion
!! The better Interaction length (Intlen) is also calculated and printed but
!! recursion use intlen=ngfftrec/2
!!
!! NOTES
!!
!! PARENTS
!!      m_rec
!!
!! CHILDREN
!!      gamma_function,initylmr,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pspnl_operat_rec(nlrec,metrec,ngfftrec,debug)

 use defs_basis
 use defs_rectypes
 use m_profiling_abi

 use m_paw_sphharm, only : initylmr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pspnl_operat_rec'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: debug
 type(metricrec_type),intent(in) ::metrec
 type(nlpsprec_type),intent(inout) :: nlrec
!arrays
 integer,intent(in) :: ngfftrec(18)
!Local variables-------------------------------
 integer :: ii,intlen
 integer :: iangol,ipsp,iproj
 integer :: mpsang,jj,kk,rr
 integer :: nfftrec
 integer :: ilmn,il,ilm,in,lmnmax
 real(dp) :: raggio,rloc,denom,step
 real(dp) :: delta_out,partial,err
 character(len=500) :: msg
 real(dp) :: part_sum(3)
 real(dp) :: ylmr_gr_dm(0,0,0)
 real(dp),allocatable :: ylmr(:,:),proj_arr(:,:,:)
 real(dp),allocatable :: radloc(:,:),nrm(:)

! *************************************************************************

 if(debug)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' pspnl_operat_rec : enter '
   call wrtout(std_out,msg,'PERS')
 end if

!#####################################################################
!--CALCULATE THE (SEMI-)MAXIMUM INTERVAL WHERE ALL THE PROJECTORS ARE
!DIFFERENT TO ZERO.
!--For any pseudo potential:
 delta_out = zero
 step = metrec%tr(1)*half !--Scanning step= grid step/2


 do ipsp = 1, nlrec%npsp !--Loop on the pseudos

!  --For any angular moment:
   do iangol = 0,maxval(nlrec%indlmn(1,:,ipsp)) !--Loop on the angular moment
     rloc = nlrec%radii(iangol+1,ipsp) !--Local radius

!    --For any projector
     do iproj = 1,nlrec%pspinfo(iangol+1,ipsp)
!      --Starting point to searching when the projector goes to zero.
!      this correspond to twice the radius wher the projector has its maximum
       raggio = two*sqrt(real(-2+2*iproj+iangol,dp))*rloc
!      --Caculate the gamma func at the denominator
       call gamma_function(real(iangol+2*iproj,dp)-half,denom)
!      --Find the zero
!      --The following while cycle should be replaced by a bisection
!      --method. Bucause this is calculated only 1 time it is not very
!      important.
       err = one
       ii=0
!      tolloc = 1.d0*abs(minval(nlrec%mat_exp_psp_nl(:nlrec%pspinfo(iangol+1,ipsp),:nlrec%pspinfo(iangol+1,ipsp),iangol+1,ipsp)))
       do while(abs(err)>1.d-2)
         raggio = raggio + step
         err = project_prec(raggio,iproj,iangol,rloc)/sqrt(denom)
         ii = ii+1
       end do
       write(std_out,*)'local delta',raggio,ii
       delta_out=maxval((/ delta_out,raggio /))
     end do !end loop on projectors

   end do !enddo on angular moment
 end do !enddo on pseudos

!--CALCULATE how many grid steps correspond to delta_out
 intlen = int(delta_out/metrec%tr(1))
!--I want that intlen is odd
 if(mod(intlen,2)==0) intlen = intlen+1

 write(msg,'(2a,i3,a)') ch10,' Interac. length of non-local psp(grid steps)=',intlen,ch10
 call wrtout(std_out,msg,'COLL')
!#####################################################################

!--Initialisation
 nfftrec = product(ngfftrec(1:3))
 lmnmax = nlrec%lmnmax
 intlen = ngfftrec(1)/2
 nlrec%intlen = intlen !--Setted in recursion variables

!#####################################################################
!--CALCULATE E(q,q')
!--Cration of the exponential*projectors*ylm matrix
 
!--Initialisation
 ABI_ALLOCATE(nlrec%projec,(nfftrec,lmnmax,nlrec%npsp))
 nlrec%projec = zero
 ABI_ALLOCATE(radloc,(3,nfftrec))
 radloc = zero
 ABI_ALLOCATE(nrm,(nfftrec))
 nrm = zero

!--Loop on pseudo types
 pseudodo: do ipsp = 1, nlrec%npsp 
!  --Control if the psp is non-local, else continue
   if(all(nlrec%pspinfo(:,ipsp)==0)) cycle 
!  --Vector which stores localy the upper part of symmetrical 
!  matrix of the exponential of the non-local operator
   mpsang = maxval(nlrec%indlmn(1,:,ipsp))+1
   ABI_ALLOCATE(proj_arr,(nfftrec,maxval(nlrec%pspinfo(:,ipsp)),mpsang))
   ABI_ALLOCATE(ylmr,(mpsang*mpsang,nfftrec))
   proj_arr = zero
   ylmr = zero

!  !debug
!  write(std_out,*)'mpsang,proj num',mpsang,maxval(nlrec%pspinfo(:,ipsp))
!  !enddebug

!  --Calculate the projctors
   do iangol = 0,mpsang-1
     rloc = nlrec%radii(iangol+1,ipsp)
     do iproj = 1,nlrec%pspinfo(iangol+1,ipsp)
       call gamma_function(real(iangol+2*iproj,dp)-half,denom)
       denom = one/sqrt(denom)
       do ii = 0,ngfftrec(1)-1 !--3-loop on coordinates
         do jj = 0,ngfftrec(2)-1
           do kk = 0,ngfftrec(3)-1
!            --Calculate the radii
             part_sum(:) = real((/ ii,jj,kk /)-intlen,dp)*(metrec%tr)
             rr = 1+ii+(jj+kk*ngfftrec(2))*ngfftrec(3)
             radloc(:,rr) = part_sum
             nrm(rr) = sqrt(sum(part_sum(:)**two))
             partial = project_prec(nrm(rr),iproj,iangol,rloc)*denom
             if(abs(partial)>tol12 ) proj_arr(rr,iproj,iangol+1) = partial
           end do
         end do
       end do  !--End 3-loop on coordinates
     end do
   end do
   

!  -------------------------------------------------------------
!  --Calculate the spherical harmonics (Verified: it works well)
   call initylmr(mpsang,1,nfftrec,nrm(:),1,radloc(:,:),ylmr(:,:),ylmr_gr_dm)
!  -------------------------------------------------------------


   do ilmn = 1,lmnmax
     ilm = nlrec%indlmn(4,ilmn,ipsp)
     il = nlrec%indlmn(1,ilmn,ipsp)+1
     in = nlrec%indlmn(3,ilmn,ipsp)
     write(msg,'(2a,i3,2i2)')ch10,'lm,l,n',ilm,il,in  
     call wrtout(std_out,msg,'COLL')

     nlrec%projec(:,ilmn,ipsp) = ylmr(ilm,:)*proj_arr(:,in,il)
   end do
   
   ABI_DEALLOCATE(ylmr)
   ABI_DEALLOCATE(proj_arr)
 end do pseudodo !--end loop on pseudo types


 ABI_DEALLOCATE(radloc)
 ABI_DEALLOCATE(nrm)

 if(debug)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' pspnl_operat_rec : exit '
   call wrtout(std_out,msg,'PERS')
 end if

 contains

   function project_prec(raggio,iproj,iangol,rloc) 
!--Analytical expression of the projectors in hgh-pspeudopotential
!--The gamma function at denominator is missing

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'project_prec'
!End of the abilint section

   real(dp) :: project_prec 
   integer,intent(in) :: iproj,iangol
   real(dp),intent(in) :: raggio,rloc

   project_prec=sqrt2*(raggio/rloc)**real((iangol+2*(iproj-1)),dp)*&
&   exp(-((raggio/rloc)**two)*half)/rloc**onehalf
 end function project_prec

end subroutine pspnl_operat_rec
!!***
