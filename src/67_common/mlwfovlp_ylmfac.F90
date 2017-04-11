!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp_ylmfac
!! NAME
!! mlwfovlp_ylmfac
!!
!! FUNCTION
!! Routine that produces a factor by which the initial
!! guess of functions will be multiplied for the Wannier90 interface.
!! It is just used if there are rotations, or if the functions required
!! are linear combinations of the ylm real functions.
!! 
!! Example,
!! For a function G(r)= 1/2 s + 1/3 px - 1/2 pz
!!   it would produce a matrix of the following form:
!!   [1/2,-1/2,1/3,0,0...0]
!!
!! The real spherical harmonics are given as factors of complex spherical harmonics 
!! The real spherical harmonics are given in table 3.1 of Wannier90 user guide.
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2017 ABINIT group (T. Rangel, DRH)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  lmax= maximum l value for spherical harmonics
!!  lmax2=number of ylm functions
!!  mband=maximum number of bands
!!  nwan = number of wannier functions
!!  proj_l(mband)= angular part of the projection function (quantum number l)
!!  proj_m(mband)= angular part of the projection function (quantum number m)
!!  proj_x(3,mband)= x axis for the projection.
!!  proj_z(3,mband)= z axis for the projection.
!!
!! OUTPUT
!!  ylmc_fac(lmax2,nwan)=matrix containig a factor for ylm hybrid orbitals
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      mlwfovlp_proj
!!
!! CHILDREN
!!      rotmat,ylm_cmplx,zgesv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mlwfovlp_ylmfac(ylmc_fac,lmax,lmax2,mband,nwan,proj_l,proj_m,proj_x,proj_z)
    
 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_paw_sphharm, only : ylm_cmplx

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mlwfovlp_ylmfac'
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in):: lmax,lmax2,nwan,mband
! arrays
 integer,intent(in) :: proj_l(mband),proj_m(mband)
 real(dp),intent(in) :: proj_x(3,mband),proj_z(3,mband)
 complex(dp),intent(out)::ylmc_fac(lmax2,nwan)
!
!Local variables-------------------------------
!
 integer :: orb_idx(16)=(/1,3,4,2,7,8,6,9,5,13,14,12,15,11,16,10/) !Tab3.1 Wannier90 user guide
 integer :: idum,ii,info,inversion_flag
 integer :: ir,iwan,jj,ll,lm,lmc,mm,mr
 real(dp):: onem,test
! arrays
 integer:: ipiv(lmax2)
 real(dp)::r(3,lmax2),rp(3,lmax2)
 real(dp)::rs2,rs3,rs6,rs12,umat(3,3)
 complex(dp)::crot(lmax2,lmax2),ctor(lmax2,lmax2),orb_lm(lmax2,-5:3,7)
 complex(dp):: ylmcp(lmax2)
 complex(dp):: ylmc_rr(lmax2,lmax2),ylmc_rr_save(lmax2,lmax2)
 complex(dp):: ylmc_rrinv(lmax2,lmax2),ylmc_rp(lmax2,lmax2)
 complex(dp),parameter :: c0=(0._dp,0._dp),c1=(1._dp,0._dp),ci=(0._dp,1._dp)
 character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!DEBUG
!write(std_out,*)'lmax ',lmax,'lmax2 ',lmax2
!write(std_out,*)'mband ',mband,'nwan ',nwan
!
!do iwan=1,nwan
!write(std_out,*)'iwan,proj_l, proj_m',proj_l(iwan),proj_m(iwan)
!write(std_out,*)'iwan,proj_x, proj_z',iwan,proj_x(:,iwan),proj_z(:,iwan)
!end do
!!END DEBUG

!constants for linear combinations of ylm's
 rs2=1._dp/sqrt(2._dp)
 rs3=1._dp/sqrt(3._dp)
 rs6=1._dp/sqrt(6._dp)
 rs12=1._dp/sqrt(12._dp)

!complex lm coefficients for real spherical harmonics in conventional order
!s, py,pz,px, dxy,dyz,dz2,dxz,dx2-y2, fy(3x2-y2),fxyz,fyz2,fz3,fxz2,
!fz(x2-y2),fx(x2-3y2)
 ctor(:,:)=c0
 do ll=0,lmax
   mm=0
   lm= ll**2+ll+mm+1
   ctor(lm,lm)=c1
   if(ll>0) then
     onem=one
     do mm=1,ll
       onem=-onem !(-1^mm)
       lm= ll**2+ll+mm+1
       lmc=ll**2+ll-mm+1
       ctor(lm ,lm )=rs2*c1
       ctor(lmc,lm )=onem*rs2*c1
       ctor(lm ,lmc)=rs2*ci
       ctor(lmc,lmc)=-onem*rs2*ci
     end do
   end if
 end do

 lm=0
 do ll=0,lmax
   do mm=-ll,ll
     lm=lm+1
     ctor(:,lm)=ctor(:,lm)*conjg(ci)**ll
   end do !mm
 end do !ll


!coefficients for basic wannier orbitals in Table 3.1 order
 orb_lm(:,:,:)=c0
 ii=0
 do ll=0,lmax
   do mr=1,2*ll+1
     ii=ii+1
     orb_lm(:,ll,mr)=ctor(:,orb_idx(ii))
   end do
 end do



!coefficients for linear combinations in table 3.2 order
 if(lmax>=1) then
!  s            px
   orb_lm(:,-1,1)=rs2*ctor(:,1)+rs2*ctor(:,4)
   orb_lm(:,-1,2)=rs2*ctor(:,1)-rs2*ctor(:,4)
!  s            px            py
   orb_lm(:,-2,1)=rs3*ctor(:,1)-rs6*ctor(:,4)+rs2*ctor(:,2)
   orb_lm(:,-2,2)=rs3*ctor(:,1)-rs6*ctor(:,4)-rs2*ctor(:,2)
   orb_lm(:,-2,3)=rs3*ctor(:,1)+2._dp*rs6*ctor(:,4)
!  s        px        py        pz
   orb_lm(:,-3,1)=half*(ctor(:,1)+ctor(:,4)+ctor(:,2)+ctor(:,3))
   orb_lm(:,-3,2)=half*(ctor(:,1)+ctor(:,4)-ctor(:,2)-ctor(:,3))
   orb_lm(:,-3,3)=half*(ctor(:,1)-ctor(:,4)+ctor(:,2)-ctor(:,3))
   orb_lm(:,-3,4)=half*(ctor(:,1)-ctor(:,4)-ctor(:,2)+ctor(:,3))
 end if
 if(lmax>=2) then
!  s            px            py
   orb_lm(:,-4,1)=rs3*ctor(:,1)-rs6*ctor(:,4)+rs2*ctor(:,2)
   orb_lm(:,-4,2)=rs3*ctor(:,1)-rs6*ctor(:,4)-rs2*ctor(:,2)
   orb_lm(:,-4,3)=rs3*ctor(:,1)+2._dp*rs6*ctor(:,4)
!  pz           dz2
   orb_lm(:,-4,4)= rs2*ctor(:,3)+rs2*ctor(:,7)
   orb_lm(:,-4,5)=-rs2*ctor(:,3)+rs2*ctor(:,7)
!  s            px            dz2         dx2-y2
   orb_lm(:,-5,1)=rs6*ctor(:,1)-rs2*ctor(:,4)-rs12*ctor(:,7)+half*ctor(:,9)
   orb_lm(:,-5,2)=rs6*ctor(:,1)+rs2*ctor(:,4)-rs12*ctor(:,7)+half*ctor(:,9)
!  s            py            dz2         dx2-y2
   orb_lm(:,-5,3)=rs6*ctor(:,1)-rs2*ctor(:,2)-rs12*ctor(:,7)-half*ctor(:,9)
   orb_lm(:,-5,4)=rs6*ctor(:,1)+rs2*ctor(:,2)-rs12*ctor(:,7)-half*ctor(:,9)
!  s            pz           dz2
   orb_lm(:,-5,5)=rs6*ctor(:,1)-rs2*ctor(:,3)+rs3*ctor(:,7)
   orb_lm(:,-5,6)=rs6*ctor(:,1)+rs2*ctor(:,3)+rs3*ctor(:,7)
 end if

!stuff complex wannier orbital coefficient array
 do iwan=1,nwan
   ylmc_fac(:,iwan)=orb_lm(:,proj_l(iwan),proj_m(iwan))
 end do


!setup to rotate ylmc_fac to new axes if called for
!skip if only s projectors are used
 if ( lmax>0 ) then
!  generate a set of nr=lmax2 random vectors
!  idum=123456
   do ir=1,lmax2
     do ii=1,3
       r(ii,ir) = uniformrandom(idum)-0.5d0
     end do !ii
     call ylm_cmplx(lmax,ylmcp,r(1,ir),r(2,ir),r(3,ir))
     ylmc_rr(ir,:)=conjg(ylmcp(:))
     ylmc_rr_save(ir,:)=conjg(ylmcp(:))
   end do !ir

   ylmc_rrinv(:,:)=c0
   do ii=1,lmax2
     ylmc_rrinv(ii,ii)=c1
   end do !ii
!  calculate inverse of ylmc(ir,lm) matrix
   call ZGESV(lmax2,lmax2,ylmc_rr,lmax2,ipiv,ylmc_rrinv,lmax2,info)

!  check that r points are independent (ie., that matrix inversion wasn't
!  too close to singular)
   ylmc_rr=matmul(ylmc_rrinv,ylmc_rr_save)
   test=zero
   do ii=1,lmax2
     ylmc_rr(ii,ii)=ylmc_rr(ii,ii)-c1
     do jj=1,lmax2
       test=max(abs(ylmc_rr(ii,jj)),test)
     end do !ii
   end do !jj
   if(test>tol8) then
     write(message, '(5a)' )&
&     '  matrix inversion error for wannier rotations',ch10,&
&     '  random vectors r(j,1:nr) are not all independent !! ',ch10,&
&     '  Action : re-seed uniformrandom or maybe just try again'
     MSG_ERROR(message)
   end if !test>tol8

!  end of the preliminaries, now to the rotations of the wannier orbitals
   do iwan=1,nwan
!    don't bother for s orbitals
     if(proj_l(iwan)==0) cycle
!    check for default axes and cycle if found
     if(proj_z(1,iwan)==zero .and. proj_z(2,iwan)==zero .and.&
&     proj_z(3,iwan)== one .and. proj_x(1,iwan)==one .and.&
&     proj_x(2,iwan)==zero .and. proj_x(3,iwan)==zero) cycle

!    get the u matrix that rotates the reference frame
     call rotmat(proj_x(:,iwan),proj_z(:,iwan),inversion_flag,umat)

!    find rotated r-vectors. Optional inversion
!    operation is an extension of the wannier90 axis-setting options
!    which only allow for proper axis rotations
     if(inversion_flag==1) then
       rp(:,:)= -matmul ( umat(:,:),  r(:,:) )
     else
       rp(:,:) = matmul ( umat(:,:) , r(:,:) )
     end if !inversion_flag

     do ir=1,lmax2
!      get the ylm representation of the rotated vectors
       call ylm_cmplx(lmax,ylmcp,rp(1,ir),rp(2,ir),rp(3,ir))
       ylmc_rp(ir,:)=conjg(ylmcp(:))
     end do !ir
!    the matrix product sum(ir) ylmc_rrinv(lm,ir)*ylmc_rp(ir,lm') gives the
!    the complex lmXlm matrix representation of the coordinate rotation
     crot(:,:)=matmul(ylmc_rrinv(:,:),ylmc_rp(:,:))

!    now rotate the current wannier orbital
     ylmcp(:)=matmul(crot(:,:),ylmc_fac(:,iwan))
     ylmc_fac(:,iwan)=ylmcp(:)

!    write(std_out,*)'ylmc_fac',ylmc_fac(:,iwan)
   end do !iwan
 end if !lmax>0

end subroutine mlwfovlp_ylmfac
!!***
