!{\src2tex{textfont=tt}}
!!****f* ABINIT/critic
!! NAME
!! critic
!!
!! FUNCTION
!!     Search for a critical point starting from point vv
!!
!! COPYRIGHT
!! Copyright (C) 2002-2018 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! INPUTS
!! aim_dtset= the structured entity containing all input variables
!! dmax= maximal step
!! sort= 0(default) general CP searching (Newton-Raphson)
!!                  -1,1,3 searching of specific type CP (Popelier)
!!
!! OUTPUT
!! ev= eigenvalues (ordered) of the Hessian in the final point
!! zz=  eigenvectors of the Hessian in the final point
!! ires= if ires==0 => CP found
!!       if ires==1 => CP not found within the maximum steps
!!
!! SIDE EFFECTS
!! vv(3)= starting point and final point
!!
!! PARENTS
!!      aim_follow,cpdrv,critics
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine critic(aim_dtset,vv,ev,zz,dmax,ires,sort)

 use defs_basis
 use defs_aimprom
 use defs_parameters
 use defs_abitypes
 use m_errors
 use m_profiling_abi

 use m_abilasi,  only : jacobi, lubksb, ludcmp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'critic'
 use interfaces_63_bader, except_this_one => critic
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: sort
 integer,intent(out) :: ires
 real(dp),intent(in) :: dmax
!arrays
 real(dp),intent(inout) :: vv(3)
 real(dp),intent(out) :: ev(3),zz(3,3)
!no_abirules
 type(aim_dataset_type), intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer :: iat,id,ii,info,ipos,istep,jii,jj,nrot
 real(dp),parameter :: evol=1.d-3
 real(dp) :: chg,dg,dltcmax,dv,dvold,rr,ss
 logical :: oscl,outof
!arrays
 integer :: ipiv(3)
 real(dp) :: dc(3),ff(3),grho(3),hrho(3,3),lp(3),vold(3),vt(3),yy(3,3)
 real(dp),allocatable :: lamb(:),pom(:,:),pom2(:,:)

!************************************************************************

!DEBUG
!write(std_out,*)' critic : enter '
!ENDDEBUG
 oscl=.false.
 if (sort==3) then
   ABI_ALLOCATE(pom,(4,4))
   ABI_ALLOCATE(pom2,(4,4))
   ABI_ALLOCATE(lamb,(4))
 elseif (sort/=0) then
   ABI_ALLOCATE(pom,(3,3))
   ABI_ALLOCATE(pom2,(3,3))
   ABI_ALLOCATE(lamb,(3))
 end if


 deb=.false.
 istep=0
 ires=0

!DEBUG
!write(std_out,'(":POSIN ",3F16.8)') vv
!do jj=1,3
!vt(jj)=rprimd(1,jj)*vv(1)+rprimd(2,jj)*vv(2)+rprimd(3,jj)*vv(3)
!end do
!write(std_out,'(":RBPOSIN ",3F16.8)') vt
!ENDDEBUG

 call vgh_rho(vv,chg,grho,hrho,rr,iat,ipos,0)

!write(std_out,'(":GRAD ",3F16.8)') grho
!write(std_out,'(":HESSIAN ",/,3F16.8,/,3F16.8,/,3F16.8)') ((hrho(ii,jj),jj=1,3),ii=1,3)

 dg=1.0_dp
 dv=1.0_dp
 dg = vnorm(grho,0)

 if (chg < aim_rhomin) then
   ires=1
!  DEBUG
!  write(std_out,*)' critic : exit, ires=1'
!  ENDDEBUG
   return
 end if

!main cycle => limits (adhoc):
!aim_dtset%lstep - minimal step
!aim_dtset%lgrad - minimal norm of gradient
!aim_maxstep - max number of steps

 do while ((dv>aim_dtset%lstep).and.(dg>aim_dtset%lgrad).and.(istep<aim_maxstep))
   istep=istep+1
   vold(:)=vv(:)
   dvold=dv
   ev(:)=0._dp
   yy(:,:)=0._dp
   call jacobi(hrho,3,3,ev,yy,nrot)   ! eigenval of Hessian
   call ordr(ev,yy,3,-1)  ! ordering

!  modification of the Newton-Raphson step to searching
!  specific type of CP (Popelier algorithm)

   ff(:)=0._dp
   lp(:)=0._dp
   dc(:)=0._dp
   outof=.false.
   do ii=1,3
     do jj=1,3
       ff(ii)=ff(ii)+yy(jj,ii)*grho(jj)
     end do
   end do
   id=sign(1._dp,ev(1))+sign(1._dp,ev(2))+sign(1._dp,ev(3))
   if (id /= sort) then
     outof=.true.
     select case (sort)
     case (-1)
       lp(3)=0.5_dp*(ev(3)-sqrt(ev(3)*ev(3)+4._dp*ff(3)*ff(3)))
       pom(:,:)=0._dp
       pom2(:,:)=0._dp
       lamb(:)=0._dp
       do ii=1,2
         pom(ii,ii)=ev(ii)
         pom(ii,3)=ff(ii)
         pom(3,ii)=ff(ii)
       end do
       call jacobi(pom,3,3,lamb,pom2,nrot)
       call ordr(lamb,pom2,3,1)
       do ii=1,3
         lp(1)=lamb(ii)
         if (abs(pom2(3,ii))>1.0d-24) exit
       end do
       lp(2)=lp(1)

!        write(std_out,*) (ev(ii),ii=1,3)
!        write(std_out,*) (lamb(ii),ii=1,3)
!        write(std_out,*) ':ID  ',id,lp(1),lp(3)

     case (1)
       lp(1)=0.5_dp*(ev(1)+sqrt(ev(1)*ev(1)+4._dp*ff(1)*ff(1)))
       pom(:,:)=0._dp
       pom2(:,:)=0._dp
       lamb(:)=0._dp
       do ii=2,3
         pom(ii-1,ii-1)=ev(ii)
         pom(ii-1,3)=ff(ii)
         pom(3,ii-1)=ff(ii)
       end do
       call jacobi(pom,3,3,lamb,pom2,nrot)
       call ordr(lamb,pom2,3,1)
       do ii=3,1,-1
         lp(2)=lamb(ii)
         if (abs(pom2(3,ii))>1.0d-24) exit
       end do
       lp(3)=lp(2)

     case (3)
       pom(:,:)=0._dp
       pom2(:,:)=0._dp
       lamb(:)=0._dp
       do ii=1,3
         pom(ii,ii)=ev(ii)
         pom(ii,4)=ff(ii)
         pom(4,ii)=ff(ii)
       end do
       call jacobi(pom,4,4,lamb,pom2,nrot)
       call ordr(lamb,pom2,4,1)
       do ii=4,1,-1
         lp(1)=lamb(ii)
         if (abs(pom2(4,ii))>1.0d-24) exit
       end do
       lp(2)=lp(1); lp(3)=lp(1)
     case default
       lp(:)=0._dp
     end select
   end if

   do ii=1,3
     if (abs(ev(ii)-lp(ii))<1.0d-24) then
       outof=.false.
       exit
     end if
   end do
   do ii=1,3                      ! SEARCHING STEP
     do jj=1,3
       if (outof) then
         dc(ii)=dc(ii)+ff(jj)*yy(ii,jj)/(ev(jj)-lp(jj))
       elseif (abs(ev(jj))>1.0d-24) then
         dc(ii)=dc(ii)+ff(jj)*yy(ii,jj)/ev(jj)
       else
         MSG_ERROR("zero eigval of Hessian")
       end if
     end do
   end do

   dltcmax = vnorm(dc,0)
   if (dltcmax>dmax) then                 ! STEP RESTRICTION
     do ii=1,3
       dc(ii)=dc(ii)*dmax/dltcmax
     end do
   end if                                  ! primitive handling of oscillations
   ss=vnorm(dc,0)                          ! usually not needed
   ss=abs(ss-dv)/ss
   if ((ss < evol).and.(oscl)) then
     dc(:)=dc(:)/2._dp
   end if


   do ii=1,3
     vv(ii) = vv(ii) - dc(ii)
   end do

!  DEBUG
!  write(std_out,'(":POSIN ",3F16.8)') vv
!  ENDDEBUG

   call vgh_rho(vv,chg,grho,hrho,rr,iat,ipos,0)
   dg = vnorm(grho,0)

   if (deb) then                 !  DEBUGG OUTPUT
     write(std_out,'("AFTER STEP ===================================")')
     write(std_out,'(":HESSIAN^(-1) ",/,3F16.8,/,3F16.8,/,3F16.8)') ((yy(ii,jii),jii=1,3),ii=1,3)
     write(std_out,'(":DC ",3F16.8)') dc
     write(std_out,*) 'STEP ',istep
     write(std_out,'(":POS ",3F16.8)') vv
     write(std_out,'(":GRAD ",3F16.8)') grho
     write(std_out,*) ':DGRAD,CHG ',dg,chg
     write(std_out,'(":HESSIAN ",/,3F16.8,/,3F16.8,/,3F16.8)') ((hrho(ii,jii),jii=1,3),ii=1,3)
   end if
   vt(:)=vv(:)-vold(:)
   dv=vnorm(vt,0)
   ss=abs(dvold-dv)/dv
   if (ss < evol) oscl=.true.
 end do

!end of main cycle

!the final output

 write(std_out,*) 'iste:',istep, dv, dg
 if (istep>=aim_maxstep)then
   write(std_out,*) ' istep=MAXSTEP ! Examine lstep2 and lgrad2 .'
   if ( (dv>aim_dtset%lstep2) .and. (dg>aim_dtset%lgrad2 )) then
     write(std_out,'(":POSOUT ",3F16.8)') vv
     ires=1
   end if
 end if

 vt(:)=vv(:)

!write(std_out,'(":POSOUT ",3F16.8)') vv
!write(std_out,'(":RBPOSOUT ",3F16.8)') vt

 call vgh_rho(vv,chg,grho,hrho,rr,iat,ipos,0)

!write(std_out,'(":GRAD ",3F16.8)') grho
!write(std_out,'(":HESSIAN ",/,3F16.8,/,3F16.8,/,3F16.8)')&
!& ((hrho(ii,jii),jii=1,3),ii=1,3)


!FINAL INVERSION OF HESSIAN

 call ludcmp(hrho,3,3,ipiv,id,info)
 if (info /= 0) then
   write(std_out,*) 'Error inverting hrho:'
   do ii=1,3
     write(std_out,*) (hrho(ii,jii),jii=1,3)
   end do
   ires=1
!  DEBUG
!  write(std_out,*)' critic : exit, ires=1'
!  ENDDEBUG
   return
!  stop 'ERROR INVERTING HESSIAN'
 end if
 do ii=1,3
   yy(ii,1:3)=0.
   yy(ii,ii)=1.
 end do
 do jii=1,3
   call lubksb(hrho,3,3,ipiv,yy(1,jii))
 end do


!write(std_out,'(":HESSIAN^(-1) ",/,3F16.8,/,3F16.8,/,3F16.8)') ((y(ii,jii),jii=1,3),ii=1,3)

 call vgh_rho(vv,chg,grho,hrho,rr,iat,ipos,0)

!write(std_out,'("LAPLAC:",F16.8)') hrho(1,1)+hrho(2,2)+hrho(3,3)

 call jacobi(hrho,3,3,ev,yy,nrot)
 call ordr(ev,yy,3,1)
 zz(:,:)=yy(:,:)

!do ii=1,3
!do jii=1,3
!zz(ii,jii)=yy(jii,ii)
!end do
!end do

!write(std_out,'(":AUTOVAL ",3F16.8)') (ev(ii),ii=1,3)
!write(std_out,'(":AUTOVEC ",/,3F16.8,/,3F16.8,/,3F16.8)') ((zz(ii,jii),ii=1,3),jii=1,3)

 if (sort/=0)  then
   ABI_DEALLOCATE(pom)
   ABI_DEALLOCATE(pom2)
   ABI_DEALLOCATE(lamb)
 end if

!DEBUG
!write(std_out,*)' critic : exit, ires= ',ires
!ENDDEBUG
end subroutine critic
!!***


!!****f* ABINIT/ordr
!! NAME
!! ordr
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2007-2018 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (to be filled)
!!
!! OUTPUT
!!  (to be filled)
!!
!! PARENTS
!!      critic
!!
!! CHILDREN
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ordr(aa,dd,nn,cff)

 use m_profiling_abi

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ordr'
!End of the abilint section

 implicit none

!Arguments ----------------------------
!scalars
 integer,intent(in) :: cff,nn
!arrays
 real(dp),intent(inout) :: aa(nn),dd(nn,nn)

!Local variables ----------------------
!scalars
 integer :: ii,jj,kk
 real(dp) :: uu

! *********************************************************************

 do ii=1,nn-1
   kk=ii
   uu=aa(ii)
   do jj=ii+1,nn
     if (cff==1) then
       if (aa(jj) >= uu+tol12) then
         kk=jj
         uu=aa(jj)
       end if
     else
       if (aa(jj) <= uu-tol12) then
         kk=jj
         uu=aa(jj)
       end if
     end if
   end do
   if (kk /= ii) then
     aa(kk)=aa(ii)
     aa(ii)=uu
     do jj=1,nn
       uu=dd(jj,ii)
       dd(jj,ii)=dd(jj,kk)
       dd(jj,kk)=uu
     end do
   end if
 end do
end subroutine ordr

!!***
