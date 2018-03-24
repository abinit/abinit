!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp_ylmfar
!! NAME
!! mlwfovlp_ylmfar
!!
!! FUNCTION
!! Routine that produces a fator by which the initial
!! guess of functions will be multiplied for the Wannier90 interface.
!! It is just used if there are rotations, or if the functions required
!! are linear combinations of the ylm real functions.
!!
!! Example,
!! For a function G(r)= 1/2 s + 1/3 px - 1/2 pz
!!   it would produce a matrix of the following form:
!!   [1/2,-1/2,1/3,0,0...0]
!!
!! This function is similar to mlwfovlp_ylmfac, but the factors it uses
!! real spherical harmonics instead of complex
!! spherical harmonics. Remember that real spherical harmonics
!! are linear combinations of complex
!! spherical harmonics
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2018 ABINIT group (T. Rangel)
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
!!
!! NOTES
!!
!! PARENTS
!!      mlwfovlp_projpaw
!!
!! CHILDREN
!!      initylmr,matrginv,rotmat
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mlwfovlp_ylmfar(ylmr_fac,lmax,lmax2,mband,nwan,proj_l,proj_m,proj_x,proj_z)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_geometry,    only : rotmat
 use m_paw_sphharm, only : initylmr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mlwfovlp_ylmfar'
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in):: lmax,lmax2,nwan,mband
! arrays
 integer,intent(in) :: proj_l(mband),proj_m(mband)
 real(dp),intent(in) :: proj_x(3,mband),proj_z(3,mband)
 real(dp),intent(out)::ylmr_fac(lmax2,nwan)
!
!Local variables-------------------------------
!
 integer :: idum,ii,inversion_flag
 integer :: ir,iwan,jj,ll,lm,mm,mr
 real(dp):: onem,test
! arrays
 real(dp),allocatable::dummy(:,:),nrm(:)
 real(dp)::r(3,lmax2),rp(3,lmax2)
 real(dp)::rs2,rs3,rs6,rs12,umat(3,3)
 real(dp)::rot(lmax2,lmax2),tor(lmax2,lmax2),orb_lm(lmax2,-5:3,7)
 real(dp):: ylmrp(lmax2)
 real(dp):: ylmr_rr(lmax2,lmax2),ylmr_rr_save(lmax2,lmax2)
 real(dp):: ylmr_rrinv(lmax2,lmax2),ylmr_rp(lmax2,lmax2)
 character(len=500) :: message                   ! to be uncommented, if needed
!no_abirules
!integer :: orb_idx(16)=(/1,3,4,2,7,8,6,9,5,13,14,12,15,11,16,10/) !Tab3.1 Wannier90 user guide

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

!
!mapping lm coefficients for real spherical harmonics
!table 3.1 of Wannier90 user guide with real spherical harmonics in routine initylmr
!s, py,pz,px, dxy,dyz,dz2,dxz,dx2-y2, fy(3x2-y2),fxyz,fyz2,fz3,fxz2,
!fz(x2-y2),fx(x2-3y2)
!note: check ordering of f orbitals, it might be wrong

 tor(:,:)=0.d0
 lm=0
 do ll=0,lmax
   do mm=-ll,ll
     onem=(-1.d0)**mm
     lm=lm+1
     if(ll == 0) then
       tor(lm,lm)=1.d0
     else
       tor(lm,lm)=onem*1.d0
     end if
   end do !mm
 end do !ll
!do lm=1,16
!write(std_out,*)'tor lm=',lm,tor(:,lm)
!end do

!coefficients for basic wannier orbitals in Table 3.1 order
 orb_lm(:,:,:)=0.d0
 ii=0
 do ll=0,lmax
   do mr=1,2*ll+1
     ii=ii+1
     orb_lm(:,ll,mr)= tor(:,ii)
!    write(std_out,*)'ii',ii,'orb_lm',orb_lm(:,ll,mr)
   end do
 end do



!coefficients for linear combinations in table 3.2 order
 if(lmax>=1) then
!  s            px
   orb_lm(:,-1,1)=rs2*tor(:,1)+rs2*tor(:,4)
   orb_lm(:,-1,2)=rs2*tor(:,1)-rs2*tor(:,4)
!  s            px            py
   orb_lm(:,-2,1)=rs3*tor(:,1)-rs6*tor(:,4)+rs2*tor(:,2)
   orb_lm(:,-2,2)=rs3*tor(:,1)-rs6*tor(:,4)-rs2*tor(:,2)
   orb_lm(:,-2,3)=rs3*tor(:,1)+2._dp*rs6*tor(:,4)
!  s        px        py        pz
   orb_lm(:,-3,1)=half*(tor(:,1)+tor(:,4)+tor(:,2)+tor(:,3))
   orb_lm(:,-3,2)=half*(tor(:,1)+tor(:,4)-tor(:,2)-tor(:,3))
   orb_lm(:,-3,3)=half*(tor(:,1)-tor(:,4)+tor(:,2)-tor(:,3))
   orb_lm(:,-3,4)=half*(tor(:,1)-tor(:,4)-tor(:,2)+tor(:,3))
 end if
 if(lmax>=2) then
!  s            px            py
   orb_lm(:,-4,1)=rs3*tor(:,1)-rs6*tor(:,4)+rs2*tor(:,2)
   orb_lm(:,-4,2)=rs3*tor(:,1)-rs6*tor(:,4)-rs2*tor(:,2)
   orb_lm(:,-4,3)=rs3*tor(:,1)+2._dp*rs6*tor(:,4)
!  pz           dz2
   orb_lm(:,-4,4)= rs2*tor(:,3)+rs2*tor(:,7)
   orb_lm(:,-4,5)=-rs2*tor(:,3)+rs2*tor(:,7)
!  s            px            dz2         dx2-y2
   orb_lm(:,-5,1)=rs6*tor(:,1)-rs2*tor(:,4)-rs12*tor(:,7)+half*tor(:,9)
   orb_lm(:,-5,2)=rs6*tor(:,1)+rs2*tor(:,4)-rs12*tor(:,7)+half*tor(:,9)
!  s            py            dz2         dx2-y2
   orb_lm(:,-5,3)=rs6*tor(:,1)-rs2*tor(:,2)-rs12*tor(:,7)-half*tor(:,9)
   orb_lm(:,-5,4)=rs6*tor(:,1)+rs2*tor(:,2)-rs12*tor(:,7)-half*tor(:,9)
!  s            pz           dz2
   orb_lm(:,-5,5)=rs6*tor(:,1)-rs2*tor(:,3)+rs3*tor(:,7)
   orb_lm(:,-5,6)=rs6*tor(:,1)+rs2*tor(:,3)+rs3*tor(:,7)
 end if

!real wannier orbital coefficient array
 do iwan=1,nwan
   ylmr_fac(:,iwan)=orb_lm(:,proj_l(iwan),proj_m(iwan))
 end do


!setup to rotate ylmr_fac to new axes if called for
!skip if only s projetors are used
 if ( lmax>0 ) then
!  generate a set of nr=lmax2 random vetors
   idum=123456
   do ir=1,lmax2
     do ii=1,3
       r(ii,ir) = uniformrandom(idum)-0.5d0
     end do !ii
   end do !ir
   ABI_ALLOCATE(nrm,(lmax2))
   nrm(:)=sqrt(r(1,:)**2+r(2,:)**2+r(3,:)**2)**0.5
   call initylmr(lmax+1,1,lmax2,nrm,1,r(:,:),ylmr_rr_save(:,:),dummy)
   ylmr_rr(:,:)=ylmr_rr_save(:,:)
   do ir=1,lmax2
     ylmr_rr_save(ir,:)=ylmr_rr(:,ir)
   end do
   ABI_DEALLOCATE(nrm)

   ylmr_rrinv(:,:)=0.d0
   do ii=1,lmax2
     ylmr_rrinv(ii,ii)=1.d0
   end do !ii
!  calculate inverse of ylmr(ir,lm) matrix
   ylmr_rrinv(:,:)=ylmr_rr_save(:,:)
   call matrginv(ylmr_rrinv,lmax2,lmax2)

!  check that r points are independent (ie., that matrix inversion wasn't
!  too close to singular)
   ylmr_rr=matmul(ylmr_rrinv,ylmr_rr_save)
   test=0.d0
   do ii=1,lmax2
     ylmr_rr(ii,ii)=ylmr_rr(ii,ii)-1.d0
     do jj=1,lmax2
       test=max(abs(ylmr_rr(ii,jj)),test)
     end do !ii
   end do !jj
   if(test>tol8) then
     write(message, '(5a)' )&
&     '  matrix inversion error for wannier rotations',ch10,&
&     '  random vetors r(j,1:nr) are not all independent !! ',ch10,&
&     '  Action : re-seed uniformrandom or maybe just try again'
     MSG_ERROR(message)
   end if !test>tol8

!  end of the preliminaries, now to the rotations of the wannier orbitals
   do iwan=1,nwan
!    don't bother for s orbitals
     if(proj_l(iwan)==0) cycle
!    check for default axes and cycle if found
     if(proj_z(1,iwan)==0.d0 .and. proj_z(2,iwan)==0.d0 .and.&
&     proj_z(3,iwan)== 1.d0 .and. proj_x(1,iwan)==1.d0 .and.&
&     proj_x(2,iwan)==0.d0 .and. proj_x(3,iwan)==0.d0) cycle

!    get the u matrix that rotates the reference frame
     call rotmat(proj_x(:,iwan),proj_z(:,iwan),inversion_flag,umat)
!
!    find rotated r-vetors. Optional inversion
!    operation is an extension of the wannier90 axis-setting options
!    which only allow for proper axis rotations
     if(inversion_flag==1) then
       rp(:,:)= -matmul ( umat(:,:),  r(:,:) )
     else
       rp(:,:) = matmul ( umat(:,:) , r(:,:) )
     end if !inversion_flag

!    get the ylm representation of the rotated vetors
     ABI_ALLOCATE(nrm,(lmax2))
     nrm(:)=sqrt(rp(1,:)**2+rp(2,:)**2+rp(3,:)**2)**0.5
     call initylmr(lmax+1,1,lmax2,nrm,1,rp(:,:),ylmr_rp(:,:),dummy)
     ylmr_rr(:,:)=ylmr_rp(:,:)
     do ir=1,lmax2
       ylmr_rp(ir,:)=ylmr_rr(:,ir)
     end do
     ABI_DEALLOCATE(nrm)
!    the matrix product sum(ir) ylmr_rrinv(lm,ir)*ylmr_rp(ir,lm') gives the
!    the  lmXlm matrix representation of the coordinate rotation

     rot(:,:)=matmul(ylmr_rrinv(:,:),ylmr_rp(:,:))
!
!    now rotate the current wannier orbital
     ylmrp(:)=matmul(rot(:,:),ylmr_fac(:,iwan))
     ylmr_fac(:,iwan)=ylmrp(:)
   end do !iwan
 end if !lmax>0

!DEBUG
!write (std_out,*) ' mlwfovlp_ylmfar : exit'
!stop
!ENDDEBUG

end subroutine mlwfovlp_ylmfar
!!***
