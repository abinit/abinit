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
!! Copyright (C) 2001-2017 ABINIT group (EB, MT, FR, SPr)
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
 use m_profiling_abi
 use m_errors

 implicit none

 private

! public procedures
 public :: rotate_mag           ! Rotate a non-collinear density wrt a magnetization
 public :: rotate_back_mag      ! Rotate back a collinear XC potential wrt a magnetization
 public :: rotate_back_mag_dfpt ! Rotate back a collinear 1st-order XC potential wrt a magnetization

!Tolerance on magnetization norm
 real(dp),parameter :: m_norm_min=tol8

!Default rotation method for DFPT
 integer,parameter :: rotation_method_default=1

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
!!
!!     Explicit formulae:
!!     1=>
!!     rho_out(1) = half*( rho_in(1) + (mag,rho_in(2:4))/|mag|) // where (*,*) is scalar product
!!     rho_out(2) = half*( rho_in(1) - (mag,rho_in(2:4))/|mag|)
!!
!!     2=>
!!     rho_out(1) = rho_in(1)
!!     rho_out(2) = half*( rho_in(1) + (mag,rho_in(2:4))/|mag|)
!!
!! OUTPUT
!!  rho_out(vectsize,2)=output (projected, collinear) (1st order if rho_in is 1st order NC density matrix)
!!  [mag_norm_out(vectsize)]= --optional-- norm of mag(:) at each point of the grid
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine rotate_mag(rho_in,rho_out,mag,vectsize,cplex,&
&                     mag_norm_in,mag_norm_out,rho_out_format) ! optional arguments


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rotate_mag'
!End of the abilint section

 implicit none

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
 !SPr TODO: make a loop on cplex inside the ipt loop; for now two separate cplex cases to render code more readable
   do ipt=1,vectsize
     if (has_mag_norm) then
       m_norm=mag_norm_in(ipt)
     else
       m_norm=sqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
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
       m_norm=sqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
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
!!
!! CHILDREN
!!
!! SOURCE

subroutine rotate_back_mag(vxc_in,vxc_out,mag,vectsize,&
&                          mag_norm_in) ! optional argument


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rotate_back_mag'
!End of the abilint section

 implicit none

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
 real(dp) :: dvdn_re,dvdn_im
 real(dp) :: dvdz_re,dvdz_im
!arrays

! *************************************************************************

!DBG_ENTER("COLL")

 has_mag_norm=present(mag_norm_in)

 do ipt=1,vectsize
   if (has_mag_norm) then
     m_norm=mag_norm_in(ipt)
   else
     m_norm=sqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
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
!!
!! CHILDREN
!!
!! SOURCE

subroutine rotate_back_mag_dfpt(vxc1_in,vxc1_out,vxc,kxc,rho1,mag,vectsize,cplex,&
&                               mag_norm_in,rot_method) ! optional arguments


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rotate_back_mag_dfpt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: vectsize
 integer,intent(in) :: cplex
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
 real(dp) :: bxc_over_m,d1,d2,d3,d4,dum,dvdn,dvdz,fact,m_dot_m1,m_norm
 real(dp) :: mdirx,mdiry,mdirz,mxy,mx1,my1,mz1,nx,ny,nz,nx1,ny1
 real(dp) :: dvdn_re,dvdn_im,dvdz_re,dvdz_im
 complex(dpc) :: rho_updn
 real(dp) :: theta0,theta1,theta1_re,theta1_im
 real(dp) :: nx1_re,nx1_im
 real(dp) :: ny1_re,ny1_im
 real(dp) :: nz1_re,nz1_im
 real(dp) :: mx1_re,mx1_im,my1_re,my1_im,mz1_re,mz1_im
 real(dp) :: m_dot_m1_re,m_dot_m1_im
 real(dp) :: fact_re,fact_im
!arrays
 real(dp)     :: vxc_diag(2),rho1_offdiag(2)
 complex(dpc) :: r1tmp(2,2),u0(2,2),u0_1(2,2),u0_1r1(2,2),u0v1(2,2)
 complex(dpc) :: rho1_updn(2,2),v1tmp(2,2),vxc1tmp(2,2)

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

 if(rotation_method==1.and.cplex==2) then
     MSG_ERROR('rotation_method=1 is not available for cplex=2 case!')
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
       rho1_updn(1,1)=rho1(ipt,1)+rho1(ipt,4)
       rho1_updn(2,2)=rho1(ipt,1)-rho1(ipt,4)
       rho1_updn(1,2)=rho1(ipt,2)-(zero,one)*rho1(ipt,3)
       rho1_updn(2,1)=rho1(ipt,2)+(zero,one)*rho1(ipt,3)
       u0_1r1=matmul(u0_1,rho1_updn)
       r1tmp=matmul(u0_1r1,u0)
       rho1_offdiag(1)=r1tmp(1,2) ; rho1_offdiag(2)=r1tmp(2,1)
       v1tmp(1,2)=-(rho1_offdiag(1)/m_norm)*(vxc_diag(1)-vxc_diag(2))
       v1tmp(2,1)= (rho1_offdiag(2)/m_norm)*(vxc_diag(2)-vxc_diag(1))
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
   !SPr Alternative method (explicitely calculated rotation matrices)
   !Vxc^(1) =   phixc^(1).Id +                                               // <= change of "electrostatic" XC potential  (phixc^(1) is denoted dvdn)
   !          + bxc^(1)*( Udag^(0).sigma_z.U^(0) )  +                        // <= this part describes the change of XC magnetic field magnitude bxc^(1)
   !          + bxc^(0)*( Udag^(1).sigma_z.U^(0) + Udag^(0).sigma_z.U^(1) )  // <= remaining terms describe the cost of magnetization rotation

   select case(cplex)
   case(1)
     do ipt=1,vectsize

       if (has_mag_norm) then
         m_norm=mag_norm_in(ipt)
       else
         m_norm=sqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
       end if

       if (m_norm>m_norm_min) then

!        dvdn is phixc^(1) (density only part)
!        dvdz is bxc^(1)   (magnetization magnitude part)
         dvdn=(vxc1_in(ipt,1)+vxc1_in(ipt,2))*half
         dvdz=(vxc1_in(ipt,1)-vxc1_in(ipt,2))*half

         mxy = sqrt(mag(ipt,1)**2+mag(ipt,2)**2)
         nx     =-mag(ipt,2)/mxy
         ny     = mag(ipt,1)/mxy
         nz     = zero

         theta0 = acos(mag(ipt,3)/m_norm)
         fact   = sin(theta0)/mxy

         theta1 =  mag(ipt,3)*(mag(ipt,1)*rho1(ipt,2)+mag(ipt,2)*rho1(ipt,3))
         theta1 =  theta1 - rho1(ipt,4)*mxy**2
         theta1 =  theta1/m_norm**2/mxy

         nx1    = mag(ipt,1)*(mag(ipt,2)*rho1(ipt,2)-mag(ipt,1)*rho1(ipt,3))
         nx1    = nx1/(mag(ipt,1)**2+mag(ipt,2)**2)
         nx1    = nx1*fact ! The factor with sin(theta) to regularize the expression
         ny1    = -(nx/ny)*nx1
         !ny1    = mag(ipt,2)*(mag(ipt,2)*rho1(ipt,2)-mag(ipt,1)*rho1(ipt,3))
         !ny1    = ny1/sqrt(mag(ipt,1)**2+mag(ipt,2)**2)/(mag(ipt,1)**2+mag(ipt,2)**2)

         !U^(0)*.Vxc1.U^(0) part
         fact=dvdz/m_norm ; dum=mag(ipt,3)*fact
         vxc1_out(ipt,1)= dvdn+dum
         vxc1_out(ipt,2)= dvdn-dum
         vxc1_out(ipt,3)= mag(ipt,1)*fact ! Real part
         vxc1_out(ipt,4)=-mag(ipt,2)*fact ! Imaginary part

         !U^(1)*.Vxc0.U^(0) + U^(0)*.Vxc0.U^(1)
         bxc_over_m = (vxc(ipt,1)-vxc(ipt,2))*half/mag(ipt,3)
         vxc1_out(ipt,1) = vxc1_out(ipt,1) - (bxc_over_m*m_norm)*sin(theta0)*theta1
         vxc1_out(ipt,2) = vxc1_out(ipt,2) + (bxc_over_m*m_norm)*sin(theta0)*theta1
         vxc1_out(ipt,3) = vxc1_out(ipt,3) + (bxc_over_m*m_norm)*(ny1+cos(theta0)*ny*theta1)
         vxc1_out(ipt,4) = vxc1_out(ipt,4) + (bxc_over_m*m_norm)*(nx1+cos(theta0)*nx*theta1)

       else ! Magnetization is zero
         dvdn=(vxc1_in(ipt,1)+vxc1_in(ipt,2))*half
         mx1=rho1(ipt,2) ; my1=rho1(ipt,3) ; mz1=rho1(ipt,4)
!        Compute Bxc/|m| from Kxc (zero limit)
         bxc_over_m = half*(half*(kxc(ipt,1)+kxc(ipt,3))-kxc(ipt,2))
         vxc1_out(ipt,1)= dvdn + bxc_over_m*mz1
         vxc1_out(ipt,2)= dvdn - bxc_over_m*mz1
         vxc1_out(ipt,3)= bxc_over_m*mx1
         vxc1_out(ipt,4)=-bxc_over_m*my1
       end if
     end do ! ipt

   case(2) !cplex=2
     ! some blablabla
     ! more blablabla
   end select

!----------------------------------------
! Explicit calculation of rotated XC functional
!----------------------------------------
 case (3)
   ! SPr 2nd method for Vxc potential rotation
   ! Explicit calculation of the rotated xc functional
   ! (derivatives of the analytical expression)
   ! Vxc^(1) =   phixc^(1).Id +                    // <= change of "electrostatic" XC potential  (phixc^(1) is denoted dvdn)
   !           + bxc^(1)*(sigma,m^(0))/|m^(0)|  +  // <= this term is equivalent to ( Udag^(0).sigma_z.U^(0) ) term in rotation_method=2
   !           +(bxc^(0)/|m^(0)|)*(sigma,m^(1)) -  // <= the last terms are equivalent to ( Udag^(1).sigma_z.U^(0) + Udag^(0).sigma_z.U^(1) )
   !           -(bxc^(0)/|m^(0)|)*(sigma,m^(0))*(m^(1),m^(0))/|m^(0)|^2
   select case(cplex)
   case(1)

     do ipt=1,vectsize

       if (has_mag_norm) then
         m_norm=mag_norm_in(ipt)
       else
         m_norm=sqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
       end if

!      dvdn is deltaVxc (density only part)
!      dvdz is deltaBxc (magnetization magnitude part)
       dvdn=(vxc1_in(ipt,1)+vxc1_in(ipt,2))*half
       dvdz=(vxc1_in(ipt,1)-vxc1_in(ipt,2))*half

       mx1=rho1(ipt,2) ; my1=rho1(ipt,3) ; mz1=rho1(ipt,4)

       if(m_norm>m_norm_min) then


         !This part describes the change of the magnitude of the xc magnetic field
         ! and the change of the scalar part of the xc electrostatic potential
         fact=dvdz/m_norm ; dum=mag(ipt,3)*fact
         vxc1_out(ipt,1)= dvdn+dum
         vxc1_out(ipt,2)= dvdn-dum
         vxc1_out(ipt,3)= mag(ipt,1)*fact ! Real part
         vxc1_out(ipt,4)=-mag(ipt,2)*fact ! Imaginary part

         !Add remaining contributions comming from the change of magnetization direction
         m_dot_m1=(mag(ipt,1)*rho1(ipt,2)+mag(ipt,2)*rho1(ipt,3)+mag(ipt,3)*rho1(ipt,4))/m_norm
         mdirx=mag(ipt,1)/m_norm
         mdiry=mag(ipt,2)/m_norm
         mdirz=mag(ipt,3)/m_norm

         bxc_over_m = (vxc(ipt,1)-vxc(ipt,2))*half/mag(ipt,3)
         vxc1_out(ipt,1) = vxc1_out(ipt,1) + bxc_over_m*( mz1 - mdirz*m_dot_m1 ) ! bxc is Bxc^(0)/|m|. In principle,
         vxc1_out(ipt,2) = vxc1_out(ipt,2) + bxc_over_m*(-mz1 + mdirz*m_dot_m1 ) ! bxc = (vxc(ipt,1)-vxc(ipt,2))/m_norm/2.0
         vxc1_out(ipt,3) = vxc1_out(ipt,3) + bxc_over_m*( mx1 - mdirx*m_dot_m1 ) ! but for small magnetization, the correct limit
         vxc1_out(ipt,4) = vxc1_out(ipt,4) + bxc_over_m*(-my1 + mdiry*m_dot_m1 ) ! is computed in rhotoxc.F90

       else
!        Compute Bxc/|m| from Kxc (zero limit)
         bxc_over_m = half*(half*(kxc(ipt,1)+kxc(ipt,3))-kxc(ipt,2))
         vxc1_out(ipt,1)= dvdn + bxc_over_m*mz1
         vxc1_out(ipt,2)= dvdn - bxc_over_m*mz1
         vxc1_out(ipt,3)= bxc_over_m*mx1
         vxc1_out(ipt,4)=-bxc_over_m*my1
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

!      dvdn is deltaVxc (density only part)
!      dvdz is deltaBxc (magnetization magnitude part)
       dvdn=(vxc1_in(ipt,1)+vxc1_in(ipt,2))*half
       dvdz=(vxc1_in(ipt,1)-vxc1_in(ipt,2))*half

       mx1=rho1(ipt,2) ; my1=rho1(ipt,3) ; mz1=rho1(ipt,4)

       if(m_norm>m_norm_min) then


         !This part describes the change of the magnitude of the xc magnetic field
         ! and the change of the scalar part of the xc electrostatic potential
         fact=dvdz/m_norm ; dum=mag(ipt,3)*fact
         vxc1_out(ipt,1)= dvdn+dum
         vxc1_out(ipt,2)= dvdn-dum
         vxc1_out(ipt,3)= mag(ipt,1)*fact ! Real part
         vxc1_out(ipt,4)=-mag(ipt,2)*fact ! Imaginary part

         !Add remaining contributions comming from the change of magnetization direction
         m_dot_m1=(mag(ipt,1)*rho1(ipt,2)+mag(ipt,2)*rho1(ipt,3)+mag(ipt,3)*rho1(ipt,4))/m_norm
         mdirx=mag(ipt,1)/m_norm
         mdiry=mag(ipt,2)/m_norm
         mdirz=mag(ipt,3)/m_norm

         bxc_over_m = (vxc(ipt,1)-vxc(ipt,2))*half/mag(ipt,3)
         vxc1_out(ipt,1) = vxc1_out(ipt,1) + bxc_over_m*( mz1 - mdirz*m_dot_m1 ) ! bxc is Bxc^(0)/|m|. In principle,
         vxc1_out(ipt,2) = vxc1_out(ipt,2) + bxc_over_m*(-mz1 + mdirz*m_dot_m1 ) ! bxc = (vxc(ipt,1)-vxc(ipt,2))/m_norm/2.0
         vxc1_out(ipt,3) = vxc1_out(ipt,3) + bxc_over_m*( mx1 - mdirx*m_dot_m1 ) ! but for small magnetization, the correct limit
         vxc1_out(ipt,4) = vxc1_out(ipt,4) + bxc_over_m*(-my1 + mdiry*m_dot_m1 ) ! is computed in rhotoxc.F90

       else
!        Compute Bxc/|m| from Kxc (zero limit)
         bxc_over_m = half*(half*(kxc(ipt,1)+kxc(ipt,3))-kxc(ipt,2))
         vxc1_out(ipt,1)= dvdn + bxc_over_m*mz1
         vxc1_out(ipt,2)= dvdn - bxc_over_m*mz1
         vxc1_out(ipt,3)= bxc_over_m*mx1
         vxc1_out(ipt,4)=-bxc_over_m*my1
       end if

     end do ! ipt

   end select

 end select

!DBG_EXIT("COLL")

end subroutine rotate_back_mag_dfpt
!!***

!----------------------------------------------------------------------

END MODULE m_xc_noncoll
!!***
