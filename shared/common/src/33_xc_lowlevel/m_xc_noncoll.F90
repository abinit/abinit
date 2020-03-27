!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_xc_noncoll
!! NAME
!! m_xc_noncoll
!!
!! FUNCTION
!!  This module provides several routines used in non-collinear XC routines
!!  (rotation of the magnetization in order to align it)
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (EB, MT, FR, SPr)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

MODULE m_xc_noncoll

 use defs_basis
 use m_abicore
 use m_errors

 implicit none

 private

! public procedures
 public :: rotate_mag           ! Rotate a non-collinear density wrt a magnetization
 public :: rotate_back_mag      ! Rotate back a collinear XC potential wrt a magnetization
 public :: rotate_back_mag_dfpt ! Rotate back a collinear 1st-order XC potential wrt a magnetization
 public :: test_rotations       ! test whether methods in rotate_back_mag_dfpt give similar results

!Tolerance on magnetization norm
 real(dp),parameter :: m_norm_min=tol8

!Default rotation method for DFPT
 integer,parameter :: rotation_method_default=3

CONTAINS

!===========================================================
!!***

!!****t* m_xc_noncoll/rotate_mag
!! NAME
!!  rotate_mag
!!
!! FUNCTION
!!  Project (rotate) a non-collinear density (stored as density+magn.)
!!  on a magnetization and give a collinear density (stored as [up,dn] or [up+dn,up]).
!!  Align both z-axis.
!!
!!
!! INPUTS
!!  rho_in(vectsize,4)=input non-collinear density and magnetization (1st or 0th order)
!!  mag(vectsize,3)=gs magnetization used for projection (0th order magnetization)
!!  vectsize=size of vector fields
!!  [mag_norm_in(vectsize)]= --optional-- norm of mag(:) at each point of the grid
!!  [rho_out_format]= 1=rho_out is stored as [up,dn]
!!                    2=rho_out is stored as [up+dn,up]
!!                    Default=1
! OUTPUT
!!  rho_out(vectsize,2)=output (projected, collinear) (1st order if rho_in is 1st order NC density matrix)
!!  [mag_norm_out(vectsize)]= --optional-- norm of mag(:) at each point of the grid
!!
!!     Explicit formulae:
!!     rho_out_format=1
!!       rho_out(1) = half*( rho_in(1) + (mag,rho_in(2:4))/|mag|) // where (*,*) is scalar product
!!       rho_out(2) = half*( rho_in(1) - (mag,rho_in(2:4))/|mag|)
!!
!!     rho_out_format=2
!!       rho_out(1) = rho_in(1)
!!       rho_out(2) = half*( rho_in(1) + (mag,rho_in(2:4))/|mag|)
!!
!!
!! PARENTS
!!      dfpt_mkvxc_noncoll,m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine rotate_mag(rho_in,rho_out,mag,vectsize,cplex,&
&                     mag_norm_in,mag_norm_out,rho_out_format) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: vectsize
 integer,intent(in) :: cplex
 integer,intent(in),optional :: rho_out_format
!arrays
 real(dp),intent(in) ::  rho_in(cplex*vectsize,4),mag(vectsize,3)
 real(dp),intent(out) :: rho_out(cplex*vectsize,2)
 real(dp),intent(in),optional :: mag_norm_in(vectsize)
 real(dp),intent(out),optional :: mag_norm_out(vectsize)

!Local variables-------------------------------
!scalars
 integer :: ipt
 logical :: has_mag_norm,out_mag_norm
 real(dp) :: m_norm,mm,rho_up,rhoin_dot_mag
 real(dp) :: rhoin_dot_mag_re,rhoin_dot_mag_im
!arrays

! *************************************************************************

!DBG_ENTER("COLL")

 has_mag_norm=present(mag_norm_in)
 out_mag_norm=present(mag_norm_out)

 if(cplex==1) then
   do ipt=1,vectsize
     if (has_mag_norm) then
       m_norm=mag_norm_in(ipt)
     else
       m_norm=dsqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
     end if

     rhoin_dot_mag=rho_in(ipt,2)*mag(ipt,1)+rho_in(ipt,3)*mag(ipt,2) &
&                 +rho_in(ipt,4)*mag(ipt,3)

     if(m_norm>m_norm_min)then
       mm=rhoin_dot_mag/m_norm
       rho_out(ipt,1)=half*(rho_in(ipt,1)+mm)
       rho_out(ipt,2)=half*(rho_in(ipt,1)-mm)
     else
       rho_out(ipt,1)=half*rho_in(ipt,1)
       rho_out(ipt,2)=half*rho_in(ipt,1)
     end if

     if (out_mag_norm) then
       if (m_norm >m_norm_min) mag_norm_out(ipt)=m_norm
       if (m_norm<=m_norm_min) mag_norm_out(ipt)=zero
     end if

   end do

 else ! cplex==2

   do ipt=1,vectsize
     if (has_mag_norm) then
       m_norm=mag_norm_in(ipt)
     else
       m_norm=dsqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
     end if

     ! real part of m.m^(1)
     rhoin_dot_mag_re=rho_in(2*ipt-1,2)*mag(ipt,1)+rho_in(2*ipt-1,3)*mag(ipt,2) &
&                    +rho_in(2*ipt-1,4)*mag(ipt,3)
     ! imaginary part of m.m^(1)
     rhoin_dot_mag_im=rho_in(2*ipt  ,2)*mag(ipt,1)+rho_in(2*ipt  ,3)*mag(ipt,2) &
&                    +rho_in(2*ipt  ,4)*mag(ipt,3)
     if(m_norm>m_norm_min)then
       mm=rhoin_dot_mag_re/m_norm
       rho_out(2*ipt-1,1)=half*(rho_in(2*ipt-1,1)+mm)
       rho_out(2*ipt-1,2)=half*(rho_in(2*ipt-1,1)-mm)
       mm=rhoin_dot_mag_im/m_norm
       rho_out(2*ipt  ,1)=half*(rho_in(2*ipt ,1)+mm)
       rho_out(2*ipt  ,2)=half*(rho_in(2*ipt ,1)-mm)
     else
       rho_out(2*ipt-1,1)=half*rho_in(2*ipt-1,1)
       rho_out(2*ipt-1,2)=half*rho_in(2*ipt-1,2)
       rho_out(2*ipt  ,1)=half*rho_in(2*ipt  ,1)
       rho_out(2*ipt  ,2)=half*rho_in(2*ipt  ,2)
     end if

     if (out_mag_norm) then
       if (m_norm >m_norm_min) mag_norm_out(ipt)=m_norm
       if (m_norm<=m_norm_min) mag_norm_out(ipt)=zero
     end if

   end do

 end if

 if (present(rho_out_format)) then
   if (rho_out_format==2) then
     do ipt=1,cplex*vectsize
       rho_up=rho_out(ipt,1)
       rho_out(ipt,1)=rho_up+rho_out(ipt,2)
       rho_out(ipt,2)=rho_up
     end do
   end if
 end if

!DBG_EXIT("COLL")

end subroutine rotate_mag
!!***

!----------------------------------------------------------------------

!!****t* m_xc_noncoll/rotate_back_mag
!! NAME
!!  rotate_back_mag
!!
!! FUNCTION
!!  Rotate back a collinear XC potential (stored as up+dn) with respect to
!!   a magnetization and give a non-collinear XC potential
!!   (stored as up_up, dn_dn, Re[up_dn], Im[up_dn])
!!  Note: works only for cplex=1
!!
!! INPUTS
!!  vxc_in(vectsize,2)=input collinear XC potential
!!  mag(vectsize,3)=gs magnetization used for projection
!!  vectsize=size of vector fields
!!  [mag_norm_in(vectsize)]= --optional-- norm of mag(:) at each point of the grid
!!
!! OUTPUT
!!  vxc_out(vectsize,4)=output non-collinear XC potential
!!
!! PARENTS
!!      dfpt_mkvxc_noncoll,m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine rotate_back_mag(vxc_in,vxc_out,mag,vectsize,&
&                          mag_norm_in) ! optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: vectsize
!arrays
 real(dp),intent(in)  :: vxc_in(vectsize,2),mag(vectsize,3)
 real(dp),intent(out) :: vxc_out(vectsize,4)
 real(dp),intent(in),optional :: mag_norm_in(vectsize)

!Local variables-------------------------------
!scalars
 integer :: ipt
 logical :: has_mag_norm
 real(dp) :: dvdn,dvdz,m_norm
!arrays

! *************************************************************************

!DBG_ENTER("COLL")

 has_mag_norm=present(mag_norm_in)

 do ipt=1,vectsize
   if (has_mag_norm) then
     m_norm=mag_norm_in(ipt)
   else
     m_norm=dsqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
   end if

   dvdn=half*(vxc_in(ipt,1)+vxc_in(ipt,2))

   if (m_norm>m_norm_min) then
     dvdz=half*(vxc_in(ipt,1)-vxc_in(ipt,2))/m_norm
     vxc_out(ipt,1)=dvdn+mag(ipt,3)*dvdz
     vxc_out(ipt,2)=dvdn-mag(ipt,3)*dvdz
     vxc_out(ipt,3)= mag(ipt,1)*dvdz
     vxc_out(ipt,4)=-mag(ipt,2)*dvdz
   else
     vxc_out(ipt,1:2)=dvdn
     vxc_out(ipt,3:4)=zero
   end if
 end do

!DBG_EXIT("COLL")

end subroutine rotate_back_mag
!!***

!!****t* m_xc_noncoll/rotate_back_mag_dfpt
!! NAME
!!  rotate_back_mag_dfpt
!!
!! FUNCTION
!!  Rotate back a 1st-order collinear XC potential (stored as up+dn) with respect to
!!   a magnetization and give a 1st-order non-collinear XC potential
!!   (stored as up_up, dn_dn, Re{up_dn}, Im{up_dn}).
!!
!! INPUTS
!!  mag(vectsize,3)=0-order magnetization used for projection
!!  rho1(vectsize,4)=1st-order non-collinear density and magnetization
!!  vxc(vectsize,4)=0-order non-collinear XC potential
!!  kxc(vectsize,nkxc)=0-order XC kernel (associated to vxc)
!!  vxc1_in(vectsize,2)=input 1st-order collinear XC potential
!!  vectsize=size of vector fields
!!  [mag_norm_in(vectsize)]= --optional-- norm of 0-order mag(:) at each point of the grid
!!  [rot_method]=Select method used to compute rotation matrix (1, 2 or 3)
!!  option=if 0, compute only the U^0 vxc^(1) U^0 part
!!         if 1, full first order xc potential
!!
!! NOTES
!!   cplex=1:
!!     V is stored as : V^11, V^22, Re[V^12], Im[V^12] (complex, hermitian)
!!     N is stored as : n, m_x, m_y, m_z               (real)
!!   cplex=2:
!!     V is stored as : V^11, V^22, V^12, i.V^21 (complex)
!!     N is stored as : n, m_x, m_y, mZ          (complex)
!!
!! OUTPUT
!!  vxc1_out(vectsize,4)=output 1st-order non-collinear XC potential
!!
!! PARENTS
!!      dfpt_mkvxc_noncoll,m_pawxc,m_xc_noncoll
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine rotate_back_mag_dfpt(option,vxc1_in,vxc1_out,vxc,kxc,rho1,mag,vectsize,cplex,&
&                               mag_norm_in,rot_method) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: vectsize
 integer,intent(in) :: cplex
 integer,intent(in) :: option
 integer,intent(in),optional  :: rot_method
!arrays
 real(dp),intent(in)          :: kxc(:,:),mag(vectsize,3),vxc(vectsize,4)
 real(dp),intent(in)          :: rho1(cplex*vectsize,4)
 real(dp),intent(in)          :: vxc1_in(cplex*vectsize,2)
 real(dp),intent(in),optional :: mag_norm_in(vectsize)
 real(dp),intent(out)         :: vxc1_out(cplex*vectsize,4)

!Local variables-------------------------------
!scalars
 integer  :: ipt,rotation_method
 logical  :: has_mag_norm
 logical  :: small_angle
 real(dp) :: bxc_over_m,d1,d2,d3,d4,dvdn,dvdz,fact,m_dot_m1,m_norm
 real(dp) :: dvdn_re,dvdn_im,dvdz_re,dvdz_im
 complex(dpc) :: rho_updn
 real(dp) :: mdirx,mdiry,mdirz,mxy,mx1,my1,mz1,wx,wy,wx1,wy1
 real(dp) :: theta0,theta1,theta1_re,theta1_im
 real(dp) :: wx1_re,wx1_im
 real(dp) :: wy1_re,wy1_im
 real(dp) :: mx1_re,mx1_im,my1_re,my1_im,mz1_re,mz1_im
 real(dp) :: m_dot_m1_re,m_dot_m1_im
 real(dp) :: bxc
!arrays
 real(dp)     :: vxc_diag(2),v21tmp(2)
 complex(dpc) :: r1tmp(2,2),u0(2,2),u0_1(2,2),u0_1r1(2,2),u0v1(2,2)
 complex(dpc) :: rho1_updn(2,2),v1tmp(2,2),vxc1tmp(2,2)
 complex(dpc) :: rho1_offdiag(2)

! *************************************************************************

!DBG_ENTER("COLL")

!Optional arguments
 has_mag_norm=present(mag_norm_in)
 rotation_method=rotation_method_default
 if (present(rot_method)) rotation_method=rot_method

!Check Kxc
 if (size(kxc)>3*vectsize) then
     MSG_ERROR('Cannot use Kxc from GGA!')
 end if

 if((rotation_method==1.or.rotation_method==2).and.cplex==2) then
     MSG_ERROR('rotation_method=1 and 2 are not available for cplex=2 case! use ixcrot=3')
 endif


 select case (rotation_method)

!----------------------------------------
! Taylor expansion of U rotation matrix
!----------------------------------------
 case (1)

   do ipt=1,vectsize

     if (has_mag_norm) then
       m_norm=mag_norm_in(ipt)
     else
       m_norm=sqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
     end if

     if(m_norm>m_norm_min) then

!      Define the U^(0) transformation matrix
       rho_updn=(mag(ipt,1)+(zero,one)*mag(ipt,2))
       d1=sqrt(( m_norm+mag(ipt,3))**2+rho_updn**2)
       d2=sqrt((-m_norm+mag(ipt,3))**2+rho_updn**2)
       d3=sqrt(( m_norm-mag(ipt,3))**2+rho_updn**2)
       d4=sqrt(( m_norm+mag(ipt,3))**2-rho_updn**2)
       u0(1,1)=( m_norm+mag(ipt,3))/d1  ! ( m  + mz)/d1
       u0(2,2)=rho_updn/d2              ! ( mx +imy)/d2
       u0(1,2)=(-m_norm+mag(ipt,3))/d2  ! (-m  + mz)/d2
       u0(2,1)=rho_updn/d1              ! ( mx +imy)/d1

!      Define the inverse of U^(0): U^(0)^-1
       u0_1(1,1)= half*d1/m_norm
       u0_1(2,2)= half*d2*(m_norm+mag(ipt,3))/(m_norm*rho_updn)
       u0_1(1,2)= half*d1*(m_norm-mag(ipt,3))/(m_norm*rho_updn)
       u0_1(2,1)=-half*d2/m_norm

!      Diagonalize the GS Vxc^(0): U^(0)^-1 Vxc^(0) U^(0)
!        (Remember the abinit notation for vxc!)
       vxc_diag(1)=half*(vxc(ipt,1)+vxc(ipt,2) &
&                 -sqrt((vxc(ipt,1)-vxc(ipt,2))**2 &
&                 +four*(vxc(ipt,3)**2+vxc(ipt,4)**2)))
       vxc_diag(2)=half*(vxc(ipt,1)+vxc(ipt,2) &
&                 +sqrt((vxc(ipt,1)-vxc(ipt,2))**2 &
&                 +four*(vxc(ipt,3)**2+vxc(ipt,4)**2)))
       v1tmp(1,1)=cmplx(real(vxc1_in(ipt,1),kind=dp),zero)
       v1tmp(2,2)=cmplx(real(vxc1_in(ipt,2),kind=dp),zero)

       !Tranforming the rhor1 with U0
       rho1_updn(1,1)=half*(rho1(ipt,1)+rho1(ipt,4))
       rho1_updn(2,2)=half*(rho1(ipt,1)-rho1(ipt,4))
       rho1_updn(1,2)=half*(rho1(ipt,2)-(zero,one)*rho1(ipt,3))
       rho1_updn(2,1)=half*(rho1(ipt,2)+(zero,one)*rho1(ipt,3))
       u0_1r1=matmul(u0_1,rho1_updn)
       r1tmp=matmul(u0_1r1,u0)
       rho1_offdiag(1)=r1tmp(1,2) ; rho1_offdiag(2)=r1tmp(2,1)

       if (option==0) then ! for xccc alone
         v1tmp(1,2)=cmplx(zero,zero)
         v1tmp(2,1)=cmplx(zero,zero)
       else
         v1tmp(1,2)=-(rho1_offdiag(1)/m_norm)*(vxc_diag(2)-vxc_diag(1))
         v1tmp(2,1)= (rho1_offdiag(2)/m_norm)*(vxc_diag(1)-vxc_diag(2))
       endif

       !Rotate back the "diagonal" xc computing the term U^(0) Vxc1_^(1) U^(0)^-1
       u0v1=matmul(u0,v1tmp)
       vxc1tmp=matmul(u0v1,u0_1)
       vxc1_out(ipt,1)=real(vxc1tmp(1,1),kind=dp)
       vxc1_out(ipt,2)=real(vxc1tmp(2,2),kind=dp)
       vxc1_out(ipt,3)=real( real(vxc1tmp(1,2)),kind=dp)
       vxc1_out(ipt,4)=real(aimag(vxc1tmp(1,2)),kind=dp)

     else ! Magnetization is zero
       dvdn=(vxc1_in(ipt,1)+vxc1_in(ipt,2))*half
       mx1=rho1(ipt,2) ; my1=rho1(ipt,3) ; mz1=rho1(ipt,4)
!      Compute Bxc/|m| from Kxc (zero limit)
       bxc_over_m = half*(half*(kxc(ipt,1)+kxc(ipt,3))-kxc(ipt,2))
       vxc1_out(ipt,1)= dvdn + bxc_over_m*mz1
       vxc1_out(ipt,2)= dvdn - bxc_over_m*mz1
       vxc1_out(ipt,3)= bxc_over_m*mx1
       vxc1_out(ipt,4)=-bxc_over_m*my1
     end if

   end do ! ipt

!----------------------------------------
! Analytical expression of U rotation matrix
!----------------------------------------
 case (2)
   !Alternative method (explicitely calculated rotation matrices)
   !Vxc^(1) =   phixc^(1).Id +                                               // <= change of "electrostatic" XC potential  (phixc^(1) is denoted dvdn)
   !          + bxc^(1)*( Udag^(0).sigma_z.U^(0) )  +                        // <= this part describes the change of XC magnetic field magnitude bxc^(1)
   !          + bxc^(0)*( Udag^(1).sigma_z.U^(0) + Udag^(0).sigma_z.U^(1) )  // <= remaining terms describe the cost of magnetization rotation

   select case(cplex)
   case(1)
     do ipt=1,vectsize

       if (has_mag_norm) then
         m_norm=mag_norm_in(ipt)
       else
         m_norm=dsqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
       end if


       mx1 =rho1(ipt,2);
       my1 =rho1(ipt,3);
       mz1 =rho1(ipt,4)

       dvdn=(vxc1_in(ipt,1)+vxc1_in(ipt,2))*half  !phixc^(1)
       dvdz=(vxc1_in(ipt,1)-vxc1_in(ipt,2))*half  !bxc^(1)

       if (m_norm>m_norm_min) then

         mxy = dsqrt(mag(ipt,1)**2+mag(ipt,2)**2)
         small_angle=(mxy/m_norm<tol8)            !condition for sin(x)~x to be valid
                                                  !even possible to set to tol6
         mdirx=mag(ipt,1)/m_norm
         mdiry=mag(ipt,2)/m_norm
         mdirz=mag(ipt,3)/m_norm

!        dvdn is phixc^(1) (density only part)
!        dvdz is bxc^(1)   (magnetization magnitude part)

         !U^(0)*.Vxc1.U^(0) part
         vxc1_out(ipt,1)= dvdn+dvdz*mdirz
         vxc1_out(ipt,2)= dvdn-dvdz*mdirz
         vxc1_out(ipt,3)= dvdz*mdirx   ! Real part
         vxc1_out(ipt,4)=-dvdz*mdiry   ! Imaginary part, minus sign comes from sigma_y

         !U^(1)*.Vxc0.U^(0) + U^(0)*.Vxc0.U^(1) part

         !bxc = dsqrt(((vxc(ipt,1)-vxc(ipt,2))*half)**2+vxc(ipt,3)**2+vxc(ipt,4)**2) !bxc^(0)
         bxc = (vxc(ipt,1)-vxc(ipt,2))*half/mag(ipt,3)*m_norm
         if (.not.small_angle) then
           wx     = mag(ipt,2)/mxy
           wy     =-mag(ipt,1)/mxy
           theta0 = dacos(mag(ipt,3)/m_norm)

           theta1 = (mdirz*(mdirx*mx1+mdiry*my1))/mxy - mz1*mxy/m_norm**2
           wx1    = (+mag(ipt,1)**2*my1 - mag(ipt,1)*mag(ipt,2)*mx1)/mxy**2/m_norm  ! wx1 multiplied by sin(theta)=mxy/m_norm
           wy1    = (-mag(ipt,2)**2*mx1 + mag(ipt,1)*mag(ipt,2)*my1)/mxy**2/m_norm  ! wx1 multiplied by sin(theta)=mxy/m_norm

           vxc1_out(ipt,1) = vxc1_out(ipt,1) - bxc*dsin(theta0)*theta1
           vxc1_out(ipt,2) = vxc1_out(ipt,2) + bxc*dsin(theta0)*theta1
           vxc1_out(ipt,3) = vxc1_out(ipt,3) - bxc*(wy1+dcos(theta0)*wy*theta1)
           vxc1_out(ipt,4) = vxc1_out(ipt,4) - bxc*(wx1+dcos(theta0)*wx*theta1)
         else
           !zero order terms O(1)
           vxc1_out(ipt,3) = vxc1_out(ipt,3) + bxc*mx1/abs(mag(ipt,3))
           vxc1_out(ipt,4) = vxc1_out(ipt,4) - bxc*my1/abs(mag(ipt,3))
           !first order terms O(theta)
           fact = bxc/(mag(ipt,3)*abs(mag(ipt,3)))
           vxc1_out(ipt,1) = vxc1_out(ipt,1) - (mag(ipt,1)*mx1+mag(ipt,2)*my1)*fact
           vxc1_out(ipt,2) = vxc1_out(ipt,2) + (mag(ipt,1)*mx1+mag(ipt,2)*my1)*fact
           vxc1_out(ipt,3) = vxc1_out(ipt,3) -  mag(ipt,1)*mz1*fact
           vxc1_out(ipt,4) = vxc1_out(ipt,4) +  mag(ipt,2)*mz1*fact
         endif

       else ! Magnetization is zero
!        Compute Bxc/|m| from Kxc (zero limit)
         bxc_over_m = half*(half*(kxc(ipt,1)+kxc(ipt,3))-kxc(ipt,2))
         vxc1_out(ipt,1)= dvdn + bxc_over_m*mz1
         vxc1_out(ipt,2)= dvdn - bxc_over_m*mz1
         vxc1_out(ipt,3)= bxc_over_m*mx1
         vxc1_out(ipt,4)=-bxc_over_m*my1
       end if
     end do ! ipt

   case(2) !cplex=2

     do ipt=1,vectsize

       if (has_mag_norm) then
         m_norm=mag_norm_in(ipt)
       else
         m_norm=dsqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
       end if

       mx1_re= rho1(2*ipt-1,2); mx1_im= rho1(2*ipt,2)
       my1_re= rho1(2*ipt-1,3); my1_im= rho1(2*ipt,3)
       mz1_re= rho1(2*ipt-1,4); mz1_im= rho1(2*ipt,4)

       dvdn_re=(vxc1_in(2*ipt-1,1)+vxc1_in(2*ipt-1,2))*half
       dvdz_re=(vxc1_in(2*ipt-1,1)-vxc1_in(2*ipt-1,2))*half
       dvdn_im=(vxc1_in(2*ipt  ,1)+vxc1_in(2*ipt  ,2))*half
       dvdz_im=(vxc1_in(2*ipt  ,1)-vxc1_in(2*ipt  ,2))*half

       if (m_norm>m_norm_min) then

         mdirx=mag(ipt,1)/m_norm
         mdiry=mag(ipt,2)/m_norm
         mdirz=mag(ipt,3)/m_norm

         mxy = dsqrt(mag(ipt,1)**2+mag(ipt,2)**2)
         small_angle=(mxy/m_norm<tol8)            !condition for sin(x)~x to be valid
                                                  !
         mdirx=mag(ipt,1)/m_norm
         mdiry=mag(ipt,2)/m_norm
         mdirz=mag(ipt,3)/m_norm

!        dvdn is phixc^(1) (density only part)
!        dvdz is bxc^(1)   (magnetization magnitude part)

         !U^(0)*.Vxc1.U^(0) part
         vxc1_out(2*ipt-1,1)= dvdn_re+dvdz_re*mdirz
         vxc1_out(2*ipt  ,1)= dvdn_im+dvdz_im*mdirz
         vxc1_out(2*ipt-1,2)= dvdn_re-dvdz_re*mdirz
         vxc1_out(2*ipt  ,2)= dvdn_im-dvdz_im*mdirz
         !NOTE: change of definition of the potential matrix components
         !      vxc1_out(:,3) =   V_updn
         !      vxc1_out(:,4) = i.V_updn
         vxc1_out(2*ipt-1,3)= dvdz_re*mdirx + dvdz_im*mdiry   !Re[  V^12]
         vxc1_out(2*ipt  ,3)= dvdz_im*mdirx - dvdz_re*mdiry   !Im[  V^12]

         !U^(1)*.Vxc0.U^(0) + U^(0)*.Vxc0.U^(1) part
         !bxc = dsqrt(((vxc(ipt,1)-vxc(ipt,2))*half)**2+vxc(ipt,3)**2+vxc(ipt,4)**2) !bxc^(0)
         bxc = (vxc(ipt,1)-vxc(ipt,2))*half/mag(ipt,3)*m_norm
         if (.not.small_angle) then
           wx     = mag(ipt,2)/mxy
           wy     =-mag(ipt,1)/mxy
           theta0 = dacos(mag(ipt,3)/m_norm)


           theta1_re = (mdirz*(mdirx*mx1_re+mdiry*my1_re))/mxy - mz1_re*mxy/m_norm**2
           theta1_im = (mdirz*(mdirx*mx1_im+mdiry*my1_im))/mxy - mz1_im*mxy/m_norm**2

           wx1_re = (+mag(ipt,1)**2*my1_re - mag(ipt,1)*mag(ipt,2)*mx1_re)/mxy**2/m_norm  ! wx1 multiplied by sin(theta)=mxy/m_norm
           wx1_im = (+mag(ipt,1)**2*my1_im - mag(ipt,1)*mag(ipt,2)*mx1_im)/mxy**2/m_norm
           wy1_re = (-mag(ipt,2)**2*mx1_re + mag(ipt,1)*mag(ipt,2)*my1_re)/mxy**2/m_norm  ! wy1 multiplied by sin(theta)=mxy/m_norm
           wy1_im = (-mag(ipt,2)**2*mx1_im + mag(ipt,1)*mag(ipt,2)*my1_im)/mxy**2/m_norm

           !U^(1)*.Vxc0.U^(0) + U^(0)*.Vxc0.U^(1)
           vxc1_out(2*ipt-1,1) = vxc1_out(2*ipt-1,1) - bxc*dsin(theta0)*theta1_re
           vxc1_out(2*ipt  ,1) = vxc1_out(2*ipt  ,1) - bxc*dsin(theta0)*theta1_im
           vxc1_out(2*ipt-1,2) = vxc1_out(2*ipt-1,2) + bxc*dsin(theta0)*theta1_re
           vxc1_out(2*ipt  ,2) = vxc1_out(2*ipt  ,2) + bxc*dsin(theta0)*theta1_im
           !cplex=1 part:
           !v12 +=   -(bxc)*(wy1+dcos(theta0)*wy*theta1)-
           !       -i.(bxc)*(wx1+dcos(theta0)*wx*theta1)
           vxc1_out(2*ipt-1,3) = vxc1_out(2*ipt-1,3) - bxc*(wy1_re+dcos(theta0)*wy*theta1_re)
           vxc1_out(2*ipt-1,3) = vxc1_out(2*ipt-1,3) + bxc*(wx1_im+dcos(theta0)*wx*theta1_im)
           vxc1_out(2*ipt  ,3) = vxc1_out(2*ipt  ,3) - bxc*(wy1_im+dcos(theta0)*wy*theta1_im)
           vxc1_out(2*ipt  ,3) = vxc1_out(2*ipt  ,3) - bxc*(wx1_re+dcos(theta0)*wx*theta1_re)
         else
           !small theta case:
           !zero order terms O(1)
           vxc1_out(2*ipt-1,3) = vxc1_out(2*ipt-1,3) + bxc*mx1_re/abs(mag(ipt,3))
           vxc1_out(2*ipt-1,3) = vxc1_out(2*ipt-1,3) + bxc*my1_im/abs(mag(ipt,3))
           vxc1_out(2*ipt  ,3) = vxc1_out(2*ipt  ,3) + bxc*mx1_im/abs(mag(ipt,3))
           vxc1_out(2*ipt  ,3) = vxc1_out(2*ipt  ,3) - bxc*my1_re/abs(mag(ipt,3))
           !first order terms:
           fact = bxc/(mag(ipt,3)*abs(mag(ipt,3)))
           vxc1_out(2*ipt-1,1) = vxc1_out(2*ipt-1,1) - (mag(ipt,1)*mx1_re+mag(ipt,2)*my1_re)*fact
           vxc1_out(2*ipt-1,2) = vxc1_out(2*ipt-1,2) + (mag(ipt,1)*mx1_re+mag(ipt,2)*my1_re)*fact
           vxc1_out(2*ipt  ,1) = vxc1_out(2*ipt  ,1) - (mag(ipt,1)*mx1_im+mag(ipt,2)*my1_im)*fact
           vxc1_out(2*ipt  ,2) = vxc1_out(2*ipt  ,2) + (mag(ipt,1)*mx1_im+mag(ipt,2)*my1_im)*fact

           vxc1_out(2*ipt-1,3) = vxc1_out(2*ipt-1,3) -  mag(ipt,1)*mz1_re*fact
           vxc1_out(2*ipt-1,3) = vxc1_out(2*ipt-1,3) -  mag(ipt,2)*mz1_im*fact
           vxc1_out(2*ipt  ,3) = vxc1_out(2*ipt  ,3) +  mag(ipt,2)*mz1_re*fact
           vxc1_out(2*ipt  ,3) = vxc1_out(2*ipt  ,3) -  mag(ipt,1)*mz1_im*fact
         endif

       else ! Magnetization is practically zero
!        Compute Bxc/|m| from Kxc (zero limit)
         bxc_over_m = half*(half*(kxc(ipt,1)+kxc(ipt,3))-kxc(ipt,2))
         vxc1_out(2*ipt-1,1)= dvdn_re + bxc_over_m*mz1_re
         vxc1_out(2*ipt  ,1)= dvdn_im + bxc_over_m*mz1_im
         vxc1_out(2*ipt-1,2)= dvdn_re - bxc_over_m*mz1_re
         vxc1_out(2*ipt  ,2)= dvdn_im - bxc_over_m*mz1_im
         vxc1_out(2*ipt-1,3)= bxc_over_m*( mx1_re+my1_im)
         vxc1_out(2*ipt  ,3)= bxc_over_m*(-my1_re+mx1_im)
       end if
       !finally reconstruct i.V^12 from V^12
       vxc1_out(2*ipt-1,4) =  vxc1_out(2*ipt  ,3)  ! Re[i.V^21] = Im[V^12]
       vxc1_out(2*ipt  ,4) =  vxc1_out(2*ipt-1,3)  ! Im[i.V^21] = Re[V^12]

     end do ! ipt

   end select

!----------------------------------------
! Explicit derivative of the rotated XC functional
!----------------------------------------
 case (3)
   ! Brute-force derivative of Vxc
   ! Explicit calculation of the rotated xc functional
   ! (derivatives of the analytical expression) (Eq. A)
   ! Vxc^(1) =   phixc^(1).Id +                    // <= change of "electrostatic" XC potential  (phixc^(1) is denoted dvdn)
   !           + bxc^(1)*(sigma,m^(0))/|m^(0)|  +  // <= this term is equivalent to ( Udag^(0).sigma_z.U^(0) ) term in rotation_method=2
   !           + bxc^(0)*(sigma,m^(1)))/|m^(0)| -  // <= the last terms are equivalent to ( Udag^(1).sigma_z.U^(0) + Udag^(0).sigma_z.U^(1) )
   !           - bxc^(0)*(sigma,m^(0))*(m^(1),m^(0))/|m^(0)|**3
   select case(cplex)
   case(1)

     do ipt=1,vectsize

       if (has_mag_norm) then
         m_norm=mag_norm_in(ipt)
       else
         m_norm=sqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
       end if

!      dvdn is phixc^(1) (density only part)
!      dvdz is bxc^(1)   (magnetization magnitude part)
       dvdn=(vxc1_in(ipt,1)+vxc1_in(ipt,2))*half
       dvdz=(vxc1_in(ipt,1)-vxc1_in(ipt,2))*half

       mx1=rho1(ipt,2) ; my1=rho1(ipt,3) ; mz1=rho1(ipt,4)

       if(m_norm>m_norm_min) then

         mdirx=mag(ipt,1)/m_norm; mdiry=mag(ipt,2)/m_norm; mdirz=mag(ipt,3)/m_norm

         !This part describes the change of the magnitude of the xc magnetic field
         !and the change of the scalar part of the xc electrostatic potential, 1st + 2nd term in Eq.A
         !phixc^(1).Id + bxc^(1) (sigma,m^(0))/|m^(0)|
         vxc1_out(ipt,1)= dvdn+dvdz*mdirz
         vxc1_out(ipt,2)= dvdn-dvdz*mdirz
         vxc1_out(ipt,3)= dvdz*mdirx   ! Real part
         vxc1_out(ipt,4)=-dvdz*mdiry   ! Imaginary part, minus sign comes from sigma_y

         if (option/=0) then
           !Add remaining contributions comming from the change of magnetization direction
           !projection of m^(1) on gs magnetization direction
           m_dot_m1=(mdirx*rho1(ipt,2)+mdiry*rho1(ipt,3)+mdirz*rho1(ipt,4))

           bxc_over_m =-dsqrt(((vxc(ipt,1)-vxc(ipt,2))*half)**2+vxc(ipt,3)**2+vxc(ipt,4)**2) !this is bxc^(0)
           bxc_over_m = bxc_over_m/m_norm
           vxc1_out(ipt,1) = vxc1_out(ipt,1) + bxc_over_m*( mz1 - mdirz*m_dot_m1 ) !
           vxc1_out(ipt,2) = vxc1_out(ipt,2) + bxc_over_m*(-mz1 + mdirz*m_dot_m1 ) !
           vxc1_out(ipt,3) = vxc1_out(ipt,3) + bxc_over_m*( mx1 - mdirx*m_dot_m1 ) !
           vxc1_out(ipt,4) = vxc1_out(ipt,4) + bxc_over_m*(-my1 + mdiry*m_dot_m1 ) !
         endif

       else
         if (option/=0) then
           !Compute bxc^(0)/|m| from kxc (|m^(0)| -> zero limit)
           bxc_over_m = half*(half*(kxc(ipt,1)+kxc(ipt,3))-kxc(ipt,2))
           vxc1_out(ipt,1)= dvdn + bxc_over_m*mz1
           vxc1_out(ipt,2)= dvdn - bxc_over_m*mz1
           vxc1_out(ipt,3)= bxc_over_m*mx1
           vxc1_out(ipt,4)=-bxc_over_m*my1
         else
           vxc1_out(ipt,1)= dvdn
           vxc1_out(ipt,2)= dvdn
           vxc1_out(ipt,3)= zero
           vxc1_out(ipt,4)= zero
         endif
       end if

     end do ! ipt

   case(2)
     !cplex=2 case

     do ipt=1,vectsize

       if (has_mag_norm) then
         m_norm=mag_norm_in(ipt)
       else
         m_norm=sqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
       end if

!      see cplex=1 case for details
       dvdn_re=(vxc1_in(2*ipt-1,1)+vxc1_in(2*ipt-1,2))*half
       dvdn_im=(vxc1_in(2*ipt  ,1)+vxc1_in(2*ipt  ,2))*half
       dvdz_re=(vxc1_in(2*ipt-1,1)-vxc1_in(2*ipt-1,2))*half
       dvdz_im=(vxc1_in(2*ipt  ,1)-vxc1_in(2*ipt  ,2))*half

       mx1_re=rho1(2*ipt-1,2); mx1_im=rho1(2*ipt,2)
       my1_re=rho1(2*ipt-1,3); my1_im=rho1(2*ipt,3)
       mz1_re=rho1(2*ipt-1,4); mz1_im=rho1(2*ipt,4)

       if(m_norm>m_norm_min) then

         mdirx=mag(ipt,1)/m_norm; mdiry=mag(ipt,2)/m_norm; mdirz=mag(ipt,3)/m_norm

         !first two terms:
         vxc1_out(2*ipt-1,1)= dvdn_re+dvdz_re*mdirz
         vxc1_out(2*ipt  ,1)= dvdn_im+dvdz_im*mdirz
         vxc1_out(2*ipt-1,2)= dvdn_re-dvdz_re*mdirz
         vxc1_out(2*ipt  ,2)= dvdn_im-dvdz_im*mdirz
         !NOTE: change of definition of the potential matrix components
         !      vxc1_out(:,3) =   V_updn
         !      vxc1_out(:,4) = i.V_dnup

         !  V^12 =   dvdz*mx/|m| - i.dvdz*my/|m| = (Re[dvdz]*mx/|m| + Im[dvdz]*my/|m|) + i.(Im[dvdz]*mx/|m| - Re[dvdz]*my/|m|) => vxc1(:,3)
         !  V^21 =   dvdz*mx/|m| + i.dvdz*my/|m| = (Re[dvdz]*mx/|m| - Im[dvdz]*my/|m|) + i.(Im[dvdz]*mx/|m| + Re[dvdz]*my/|m|)
         vxc1_out(2*ipt-1,3)= dvdz_re*mdirx + dvdz_im*mdiry   !Re[V^12]
         vxc1_out(2*ipt  ,3)= dvdz_im*mdirx - dvdz_re*mdiry   !Im[V^12]
         vxc1_out(2*ipt-1,4)= dvdz_re*mdirx - dvdz_im*mdiry   !Re[V^21]
         vxc1_out(2*ipt  ,4)= dvdz_im*mdirx + dvdz_re*mdiry   !Im[V^21]
         if (option/=0) then

           !remaining contributions:
           m_dot_m1_re= mdirx*mx1_re + mdiry*my1_re + mdirz*mz1_re
           m_dot_m1_im= mdirx*mx1_im + mdiry*my1_im + mdirz*mz1_im

           bxc_over_m =-dsqrt(((vxc(ipt,1)-vxc(ipt,2))*half)**2+vxc(ipt,3)**2+vxc(ipt,4)**2) !this is bxc^(0)
           bxc_over_m = bxc_over_m/m_norm
          !bxc_over_m = (vxc(ipt,1)-vxc(ipt,2))*half/mag(ipt,3)

           vxc1_out(2*ipt-1,1) = vxc1_out(2*ipt-1,1) + bxc_over_m*( mz1_re - mdirz*m_dot_m1_re ) ! Re[V^11]
           vxc1_out(2*ipt  ,1) = vxc1_out(2*ipt  ,1) + bxc_over_m*( mz1_im - mdirz*m_dot_m1_im ) ! Im[V^11]
           vxc1_out(2*ipt-1,2) = vxc1_out(2*ipt-1,2) + bxc_over_m*(-mz1_re + mdirz*m_dot_m1_re ) ! Re[V^22]
           vxc1_out(2*ipt  ,2) = vxc1_out(2*ipt  ,2) + bxc_over_m*(-mz1_im + mdirz*m_dot_m1_im ) ! Im[V^22]

           !    v12  += bxc_over_m*(   (mx1    - mdirx*m_dot_m1   ) - i.( my1    - mdiry*m_dot_m1   )   )  <= see cplex=1
           ! Re[v12] += bxc_over_m*(   (mx1_re - mdirx*m_dot_m1_re) +   ( my1_im - mdiry*m_dot_m1_im)   )
           ! Im[v12] += bxc_over_m*(   (mx1_im - mdirx*m_dot_m1_im) +   (-my1_re + mdiry*m_dot_m1_re)   )
           vxc1_out(2*ipt-1,3) = vxc1_out(2*ipt-1,3) + bxc_over_m*( mx1_re - mdirx*m_dot_m1_re ) ! Re[V^12]
           vxc1_out(2*ipt-1,3) = vxc1_out(2*ipt-1,3) + bxc_over_m*( my1_im - mdiry*m_dot_m1_im ) ! Re[V^12]
           vxc1_out(2*ipt  ,3) = vxc1_out(2*ipt  ,3) + bxc_over_m*( mx1_im - mdirx*m_dot_m1_im ) ! Im[V^12]
           vxc1_out(2*ipt  ,3) = vxc1_out(2*ipt  ,3) + bxc_over_m*(-my1_re + mdiry*m_dot_m1_re ) ! Im[V^12]

           !    v21  += bxc_over_m*(   (mx1    - mdirx*m_dot_m1   ) + i.( my1    - mdiry*m_dot_m1   )   )
           ! Re[v21] += bxc_over_m*(   (mx1_re - mdirx*m_dot_m1_re) +   (-my1_im + mdiry*m_dot_m1_im)   )
           ! Im[v21] += bxc_over_m*(   (mx1_im - mdirx*m_dot_m1_im) +   ( my1_re - mdiry*m_dot_m1_re)   )
           ! the 4th component is actually not v21, but rather i.v21, this will be adjusted later
           vxc1_out(2*ipt-1,4) = vxc1_out(2*ipt-1,4) + bxc_over_m*( mx1_re - mdirx*m_dot_m1_re ) ! Re[V^21]
           vxc1_out(2*ipt-1,4) = vxc1_out(2*ipt-1,4) + bxc_over_m*(-my1_im + mdiry*m_dot_m1_im ) ! Re[V^21]
           vxc1_out(2*ipt  ,4) = vxc1_out(2*ipt  ,4) + bxc_over_m*( mx1_im - mdirx*m_dot_m1_im ) ! Im[V^21]
           vxc1_out(2*ipt  ,4) = vxc1_out(2*ipt  ,4) + bxc_over_m*( my1_re - mdiry*m_dot_m1_re ) ! Im[V^21]
         endif
       else
         if(option/=0) then
           !Compute Bxc/|m| from Kxc (|m^(0)| -> zero limit)
           bxc_over_m = half*(half*(kxc(ipt,1)+kxc(ipt,3))-kxc(ipt,2))
           vxc1_out(2*ipt-1,1)= dvdn_re + bxc_over_m*mz1_re
           vxc1_out(2*ipt-1,2)= dvdn_re - bxc_over_m*mz1_re
           vxc1_out(2*ipt  ,1)= dvdn_im + bxc_over_m*mz1_im
           vxc1_out(2*ipt  ,2)= dvdn_im - bxc_over_m*mz1_im

           vxc1_out(2*ipt-1,3)= bxc_over_m*(mx1_re+my1_im)
           vxc1_out(2*ipt  ,3)= bxc_over_m*(mx1_im-my1_re)
           vxc1_out(2*ipt-1,4)= bxc_over_m*(mx1_re-my1_im)
           vxc1_out(2*ipt  ,4)= bxc_over_m*(mx1_im+my1_re)
         else
           vxc1_out(2*ipt-1,1)= dvdn_re
           vxc1_out(2*ipt-1,2)= dvdn_re
           vxc1_out(2*ipt  ,1)= dvdn_im
           vxc1_out(2*ipt  ,2)= dvdn_im

           vxc1_out(2*ipt-1,3)= zero
           vxc1_out(2*ipt  ,3)= zero
           vxc1_out(2*ipt-1,4)= zero
           vxc1_out(2*ipt  ,4)= zero
         endif
       end if

       !finally reconstruct i.V^21 from V^21
       v21tmp(1) = vxc1_out(2*ipt-1,4) !Re[V^21]
       v21tmp(2) = vxc1_out(2*ipt  ,4) !Im[V^21]

       vxc1_out(2*ipt-1,4) =-v21tmp(2) ! Re[i.V^21]=-Im[V^21]
       vxc1_out(2*ipt  ,4) = v21tmp(1) ! Im[i.V^21]= Re[V^21]

     end do ! ipt

   end select !cplex

 end select ! rotation_method

!DBG_EXIT("COLL")

end subroutine rotate_back_mag_dfpt
!!***


!!****f* ABINIT/m_xc_noncoll/test_rotations
!! NAME
!!  test_rotations
!!
!! FUNCTION
!!  Test three different methods in rotate_back_mag_dfpt
!!
!! INPUTS
!!  option= types of tests to perform
!!          0=> only quick tests
!!          1=> quick and slow tests
!!  cplex = complex or real potential and first order magnetization
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!!  For debug purposes
!!
!! PARENTS
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine test_rotations(option,cplex)

!Arguments ------------------------------------
 integer , intent(in)  :: option
 integer , intent(in)  :: cplex

!Local variables-------------------------------
 real(dp) :: m0(1,3),vxc0(1,4),kxc(1,3)
 real(dp) :: n1(cplex,4),vxc1_in(cplex,4),vxc1_out(cplex,4)
 real(dp) :: delta_23(cplex,4) !,delta_12(cplex,4)
 real(dp) :: m0_norm,dvdn,dvdz,err23 !,wrong_comp!,err12
 real(dp) :: theta0,phi0,theta1,phi1,err,m1_norm
 integer  :: dir0,dir1
! *************************************************************************

 DBG_ENTER("COLL")

! if (option/=1 .and. option/=2 ) then
!  write(msg,'(3a,i0)')&
!&  'The argument option should be 1 or 2,',ch10,&
!&  'however, option=',option
!  MSG_BUG(msg)
! end if
!
! if (sizein<1) then
!  write(msg,'(3a,i0)')&
!&  '  The argument sizein should be a positive number,',ch10,&
!&  '  however, sizein=',sizein
!  MSG_ERROR(msg)
! end if

 DBG_EXIT("COLL")

 !write(*,*)    'VXC_NONCOLL TESTS================================================================'
 if (cplex==1) then
    !write(*,*) '  cplex=1------------------------------------------------------------------------'

    !write(*,*) '    TEST: simple  m* orietnations, bxc^(1) part'
    dvdn=zero;dvdz=1.0!
    err23=zero
    do dir0=1,3
    m0=zero; n1=zero
      do dir1=2,4
        m0(1,dir0)=0.1
        m0_norm=sqrt(m0(1,1)**2+m0(1,2)**2+m0(1,3)**2)
        n1(1,dir1)=0.8     ! any number would do here

        vxc0=zero;     ! no bxc^(0) part at all

        vxc1_in(1,1)= dvdn+dvdz
        vxc1_in(1,2)= dvdn-dvdz

        call rotate_back_mag_dfpt(1,vxc1_in,vxc1_out,vxc0,kxc,n1,m0,1,1,rot_method=1)
        call rotate_back_mag_dfpt(1,vxc1_in,vxc1_out,vxc0,kxc,n1,m0,1,1,rot_method=2)
        delta_23=vxc1_out
        call rotate_back_mag_dfpt(1,vxc1_in,vxc1_out,vxc0,kxc,n1,m0,1,1,rot_method=3)
        delta_23=abs(delta_23-vxc1_out)
        err=max(delta_23(1,1),delta_23(1,2),delta_23(1,3),delta_23(1,4))
        if (err23<err) err23=err;

      enddo
    enddo
    !write(*,*) '    maximum mismatch between methods 2 and 3:',err23

    !write(*,*) '    TEST: simple  m* orietnations, bxc^(0) part'

    err23=zero
    dvdn=zero;dvdz=1.0
    do dir0=1,3
    m0=zero; n1=zero
      do dir1=2,4
        m0(1,dir0)=0.1
        m0_norm=sqrt(m0(1,1)**2+m0(1,2)**2+m0(1,3)**2)
        n1(1,dir1)=0.8     ! =m^1, any number would do here

        vxc0(1,1) = dvdn+dvdz*m0(1,3)/m0_norm
        vxc0(1,2) = dvdn-dvdz*m0(1,3)/m0_norm
        vxc0(1,3) = dvdz*m0(1,1)/m0_norm
        vxc0(1,4) =-dvdz*m0(1,2)/m0_norm

        vxc1_in=zero !vxc^(1) collinear is zero

        call rotate_back_mag_dfpt(1,vxc1_in,vxc1_out,vxc0,kxc,n1,m0,1,1,rot_method=1)
        call rotate_back_mag_dfpt(1,vxc1_in,vxc1_out,vxc0,kxc,n1,m0,1,1,rot_method=2)
        delta_23=vxc1_out
        call rotate_back_mag_dfpt(1,vxc1_in,vxc1_out,vxc0,kxc,n1,m0,1,1,rot_method=3)
        delta_23=abs(delta_23-vxc1_out)
        err=maxval(abs(delta_23(1,:)))
        if (err23<err) err23=err;
      enddo
    enddo
    !write(*,*) '    maximum mismatch between methods 2 and 3:',err23

    !write(*,*) '    TEST: general m0 orietnations, bxc^(0) part'

    theta0=zero
    err23=zero
    m0_norm=0.3
    do while(theta0<=pi)
      phi0=zero
      do while(phi0<=2*pi)
        m0(1,1)=m0_norm*sin(theta0)*cos(phi0)
        m0(1,2)=m0_norm*sin(theta0)*sin(phi0)
        m0(1,3)=m0_norm*cos(theta0)

        do  dir1=2,4
          n1=zero
          n1(1,dir1)=0.8     ! =m^1, any number would do here

          !vxc0=zero;     !
          vxc0(1,1) = dvdn+dvdz*m0(1,3)/m0_norm
          vxc0(1,2) = dvdn-dvdz*m0(1,3)/m0_norm
          vxc0(1,3) = dvdz*m0(1,1)/m0_norm
          vxc0(1,4) =-dvdz*m0(1,2)/m0_norm

          vxc1_in=zero

          !call rotate_back_mag_dfpt(vxc1_in,vxc1_out,vxc0,kxc,n1,m0,1,1,rot_method=1)
          call rotate_back_mag_dfpt(1,vxc1_in,vxc1_out,vxc0,kxc,n1,m0,1,1,rot_method=2)
          delta_23=vxc1_out
          call rotate_back_mag_dfpt(1,vxc1_in,vxc1_out,vxc0,kxc,n1,m0,1,1,rot_method=3)
          delta_23=abs(delta_23-vxc1_out)
          err=maxval(abs(delta_23(1,:)))
          if (err23<err) err23=err;
        enddo
        phi0=phi0+2*pi/100.0
      enddo
      theta0=theta0+pi/100.0
    enddo
    !write(*,*) '    maximum mismatch between methods 2 and 3:',err23

    if(option==2) then
    !write(*,*) '    TEST: general m* orietnations, bxc^(0) part'
    dvdn=zero;dvdz=1.0

    theta0=zero
    err23=zero
    m0_norm=0.3
    m1_norm=10.5
    do while(theta0<=pi) !loops on orientation of m^(0)
      phi0=zero
      do while(phi0<=2*pi)
        m0(1,1)=m0_norm*sin(theta0)*cos(phi0)
        m0(1,2)=m0_norm*sin(theta0)*sin(phi0)
        m0(1,3)=m0_norm*cos(theta0)

        vxc0(1,1) = dvdn+dvdz*m0(1,3)/m0_norm
        vxc0(1,2) = dvdn-dvdz*m0(1,3)/m0_norm
        vxc0(1,3) = dvdz*m0(1,1)/m0_norm
        vxc0(1,4) =-dvdz*m0(1,2)/m0_norm

        vxc1_in=zero

        theta1=zero
        do while(theta1<=pi) !loops on orientation of m^(1)
          phi1=zero
          do while(phi1<=2*pi)
            n1(1,1)=zero
            n1(1,2)=m1_norm*sin(theta1)*cos(phi1)
            n1(1,3)=m1_norm*sin(theta1)*sin(phi1)
            n1(1,4)=m1_norm*cos(theta1)

            !vxc0=zero;     !
            !call rotate_back_mag_dfpt(vxc1_in,vxc1_out,vxc0,kxc,n1,m0,1,1,rot_method=1)
            call rotate_back_mag_dfpt(1,vxc1_in,vxc1_out,vxc0,kxc,n1,m0,1,1,rot_method=2)
            delta_23=vxc1_out
            call rotate_back_mag_dfpt(1,vxc1_in,vxc1_out,vxc0,kxc,n1,m0,1,1,rot_method=3)
            delta_23=abs(delta_23-vxc1_out)
            err=maxval(abs(delta_23(1,:)))
            if (err23<err) err23=err;
            phi1=phi1+2*pi/100.0
          enddo
          theta1=theta1+pi/100.0
        enddo

        phi0=phi0+2*pi/100.0
      enddo
      theta0=theta0+pi/100.0
    enddo
    !write(*,*) '    maximum mismatch between methods 2 and 3:',err23
    endif

 !else !cplex=2

 endif

end subroutine test_rotations
!!***


!----------------------------------------------------------------------

END MODULE m_xc_noncoll
!!***
