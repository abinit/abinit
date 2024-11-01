!!****m* ABINIT/m_cutoff_cylinder
!! NAME
!!  m_cutoff_cylinder
!!
!! FUNCTION
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_cutoff_cylinder

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_splines
 use m_sort

 use m_fstrings,        only : sjoin, itoa
 use m_geometry,        only : normv, metric
 use m_bessel,          only : CALJY0, CALJY1, CALCK0, CALCK1
 use m_numeric_tools,   only : OPERATOR(.x.), quadrature
 use m_paw_numeric,     only : paw_jbessel

 implicit none

 private
!!***

 public :: cutoff_cylinder, K0cos

 !integer,public,parameter :: CYLINDER_BEIGI = 1
 !integer,public,parameter :: CYLINDER_ROZZI = 2
!!***

! private variables used for the integration needed by the cylindrical case.
 integer,save :: npts_,ntrial_,qopt_
 real(dp),save :: ha_,hb_,r0_
 real(dp),save :: qpg_perp_,qpg_para_,qpgx_,qpgy_
 real(dp),save :: zz_,xx_, rho_
 real(dp),save :: hcyl_,rcut_,accuracy_

CONTAINS
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/cutoff_cylinder
!! NAME
!! cutoff_cylinder
!!
!! FUNCTION
!!  Calculate the Fourier components of an effective Coulomb interaction
!!  zeroed outside a finite cylindrical region. Two methods are implemented:
!!
!!   method==1: The interaction in the (say) x-y plane is truncated outside the Wigner-Seitz
!!              cell centered on the wire in the x-y plane. The interaction has infinite
!!              extent along the z axis and the Fourier transform is singular only at the Gamma point.
!!              Only orthorombic Bravais lattices are supported.
!!   method==2: The interaction is truncated outside a cylinder of radius rcut. The cylinder has finite
!!              extent along z. No singularity occurs.
!!
!! INPUTS
!!  boxcenter(3)= center of the wire in the x-y axis
!!  qpt(3)= q-point
!!  ng=number of G vectors
!!  gvec(3,ng)=G vectors in reduced coordinates
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  method=1 for Beigi approach (infinite cylinder with interaction truncated outside the W-S cell)
!!         2 for Rozzi method (finite cylinder)
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  vc_cut(ng)= Fourier components of the effective Coulomb interaction
!!
!! SOURCE

subroutine cutoff_cylinder(qpt, ng, gvec, rcut, hcyl, pdir, boxcenter, rprimd, vc_cut, method, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng,method,comm
 real(dp),intent(in) :: rcut,hcyl
!arrays
 integer,intent(in) :: gvec(3,ng),pdir(3)
 real(dp),intent(in) :: boxcenter(3),qpt(3),rprimd(3,3)
 real(dp),intent(out) :: vc_cut(ng)

!Local variables-------------------------------
!scalars
 integer,parameter :: N0=1000
 integer :: ig,igs,ierr, my_rank, nproc
 real(dp) :: j0,j1,k0,k1,qpg2,qpg_xy,tmp
 real(dp) :: qpg_z,quad,rcut2,hcyl2,c1,c2,ucvol,SMALL
 logical :: q_is_zero
 character(len=500) :: msg
!arrays
 real(dp) :: qpg(3),b1(3),b2(3),b3(3),gmet(3,3),rmet(3,3),gprimd(3,3),qc(3),gcart(3)
!************************************************************************

 ABI_UNUSED(pdir)
 ABI_UNUSED(boxcenter)

 ! ===================================================
 ! === Setup for the quadrature of matrix elements ===
 ! ===================================================
 qopt_    =6         ! Quadrature method, see quadrature routine.
 ntrial_  =30        ! Max number of attempts.
 accuracy_=0.001     ! Fractional accuracy required.
 npts_    =6         ! Initial number of point (only for Gauss-Legendre method).
 SMALL    =tol4      ! Below this value (q+G)_i is treated as zero.
 rcut_    =rcut      ! Radial cutoff, used only if method==2
 hcyl_    =hcyl      ! Lenght of cylinder along z, only if method==2

 !write(msg,'(3a,2(a,i5,a),a,f8.5)')ch10,&
 ! ' cutoff_cylinder: Info on the quadrature method : ',ch10,&
 ! '  Quadrature scheme      = ',qopt_,ch10,&
 ! '  Max number of attempts = ',ntrial_,ch10,&
 ! '  Fractional accuracy    = ',accuracy_
 !call wrtout(std_out, msg)

 ! From reduced to Cartesian coordinates.
 call metric(gmet, gprimd, -1, rmet, rprimd, ucvol)
 b1(:) =two_pi*gprimd(:,1)
 b2(:) =two_pi*gprimd(:,2)
 b3(:) =two_pi*gprimd(:,3)

 qc = b1(:)*qpt(1) + b2(:)*qpt(2) + b3(:)*qpt(3)

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! ================================================
 ! === Different approaches according to method ===
 ! ================================================
 vc_cut = zero

 select case (method)

 case (1)
   ! Infinite cylinder, interaction is zeroed outside the Wigner-Seitz cell.
   ! NB: Beigi's expression holds only if the BZ is sampled only along z.
   !call wrtout(std_out, 'cutoff_cylinder: Using Beigi''s Infinite cylinder')

   if (ANY(qc(1:2) > SMALL)) then
     write(msg,'(5a)')&
      ' found q-points with non zero components in the X-Y plane. ',ch10,&
      ' This is not allowed, see Notes in cutoff_cylinder.F90. ',ch10,&
      ' ACTION: Modify the q-point sampling. '
     ABI_ERROR(msg)
   end if

   ! Check if Bravais lattice is orthorombic and parallel to the Cartesian versors.
   ! In this case the intersection of the WS cell with the x-y plane is a rectangle with -ha_<=x<=ha_ and -hb_<=y<=hb_
   if ((ANY(ABS(rprimd(2:3,  1)) > tol6)) .or. &
       (ANY(ABS(rprimd(1:3:2,2)) > tol6)) .or. &
       (ANY(ABS(rprimd(1:2,  3)) > tol6))) then
     ABI_ERROR('Bravais lattice should be orthorombic and parallel to the Cartesian versors')
   end if

   ha_ = half*SQRT(DOT_PRODUCT(rprimd(:,1),rprimd(:,1)))
   hb_ = half*SQRT(DOT_PRODUCT(rprimd(:,2),rprimd(:,2)))
   r0_ = MIN(ha_,hb_)/N0

   ! For each (q,G) pair evaluate the integral defining the Coulomb cutoff.
   ! NB: the code assumes that all q-vectors are non zero and q_xy/=0.
   igs=1
   ! Skip singularity at Gamma, it will be treated "by hand" in csigme.
   q_is_zero = (normv(qpt, gmet, 'G') < tol4)

   do ig=igs,ng
     if (mod(ig, nproc) /= my_rank) cycle ! MPI parallelism

     gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
     qpg(:)=qc(:)+gcart(:)
     qpgx_=qpg(1); qpgy_=qpg(2); qpg_para_=ABS(qpg(3))
     !write(std_out,*)"qpgx_=",qpgx_, "qpgy_=",qpgy_, "qpg_para=",qpg_para_

     ! Avoid singularity in K_0{qpg_para_\rho) by using a small q along the periodic dimension.
     if (q_is_zero .and. qpg_para_ < tol6) qpg_para_ = tol6

     ! Calculate $ 2\int_{WS} dxdy K_0{qpg_para_\rho) cos(x.qpg_x + y.qpg_y) $
     ! where WS is the Wigner-Seitz cell.
     tmp=zero

     ! Difficult part, integrate on a small cirle of radius r0 using spherical coordinates
     !call quadrature(K0cos_dth_r0,zero,r0_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
     !ABI_CHECK(ierr == 0, "Accuracy not reached")
     !write(std_out,'(i8,a,es14.6)')ig,' 1 ',quad
     !tmp=tmp+quad
     ! Add region with 0<=x<=r0 and y>=+-(SQRT(r0^2-x^2))since WS is rectangular
     !call quadrature(K0cos_dy_r0,zero,r0_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
     !ABI_CHECK(ierr == 0, "Accuracy not reached")
     !write(std_out,'(i8,a,es14.6)')ig,' 2 ',quad
     !tmp=tmp+quad
     ! Get the in integral in the rectangle with x>=r0, should be the easiest but sometimes has problems to converge
     !call quadrature(K0cos_dy,r0_,ha_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
     !ABI_CHECK(ierr == 0, "Accuracy not reached")
     !write(std_out,'(i8,a,es14.6)')ig,' 3 ',quad
     !
     ! More stable method: midpoint integration with Romberg extrapolation ===
     call quadrature(K0cos_dy,zero,ha_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
     !write(std_out,'(i8,a,es14.6)')ig,' 3 ',quad
     ABI_CHECK(ierr == 0, "Accuracy not reached in quadrature!")

     ! Store final result
     ! Factor two comes from the replacement WS -> (1,4) quadrant thanks to symmetries of the integrad.
     tmp = tmp+quad
     vc_cut(ig) = two*(tmp*two)
   end do ! ig

 case (2)
   ! Finite cylinder of length hcyl from Rozzi et al.
   ! TODO add check on hcyl value that should be smaller that 1/deltaq
   if (hcyl_ < zero) then
     write(msg,'(a,f8.4)')' Negative value for cylinder length hcyl_=',hcyl_
     ABI_BUG(msg)
   end if

   if (ABS(hcyl_) > tol12) then
     !write(std_out,'(2(a,f8.4))')' cutoff_cylinder: using finite cylinder of length= ',hcyl_,' rcut= ',rcut_
     hcyl2=hcyl_**2
     rcut2=rcut_**2

     ! No singularity occurs in finite cylinder, thus start from 1.
     do ig=1,ng
       if (mod(ig, nproc) /= my_rank) cycle ! MPI parallelism

       gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
       qpg(:)=qc(:)+gcart(:)
       qpg_para_=ABS(qpg(3)) ; qpg_perp_=SQRT(qpg(1)**2+qpg(2)**2)

       if (qpg_perp_ /= zero .and. qpg_para_ /= zero) then
         ! $ 4\pi\int_0^{R_c} d\rho\rho j_o(qpg_perp_.\rho)\int_0^hcyl dz\cos(qpg_para_*z)/sqrt(\rho^2+z^2) $
         call quadrature(F2,zero,rcut_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
         ABI_CHECK(ierr == 0, "Accuracy not reached")
         vc_cut(ig) = four_pi*quad

       else if (qpg_perp_ == zero .and. qpg_para_ /= zero) then
         ! $ \int_0^h sin(qpg_para_.z)/\sqrt(rcut^2+z^2)dz $
         call quadrature(F3,zero,hcyl_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
         ABI_CHECK(ierr == 0, "Accuracy not reached")

         c1=one/qpg_para_**2-COS(qpg_para_*hcyl_)/qpg_para_**2-hcyl_*SIN(qpg_para_*hcyl_)/qpg_para_
         c2=SIN(qpg_para_*hcyl_)*SQRT(hcyl2+rcut2)
         vc_cut(ig) = four_pi*c1+four_pi*(c2-quad)/qpg_para_

       else if (qpg_perp_ /= zero .and. qpg_para_ == zero) then
         ! $ 4pi\int_0^rcut d\rho \rho J_o(qpg_perp_.\rho) ln((h+\sqrt(h^2+\rho^2))/\rho) $
         call quadrature(F4,zero,rcut_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
         ABI_CHECK(ierr == 0, "Accuracy not reached")
         vc_cut(ig) = four_pi*quad

       else if (qpg_perp_ == zero .and. qpg_para_ == zero) then
         ! Use lim q+G --> 0
         vc_cut(ig) = two_pi*(-hcyl2+hcyl_*SQRT(hcyl2+rcut2)+rcut2*LOG((hcyl_+SQRT(hcyl_+SQRT(hcyl2+rcut2)))/rcut_))

       else
         ABI_BUG('You should not be here!')
       end if

     end do !ig

   else
     ! Infinite cylinder.
     !call wrtout(std_out, ' cutoff_cylinder: using Rozzi''s method with infinite cylinder ')

     do ig=1,ng
       if (mod(ig, nproc) /= my_rank) cycle ! MPI parallelism

       gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
       qpg(:)=qc(:)+gcart(:)
       qpg2  =DOT_PRODUCT(qpg,qpg)
       qpg_z =ABS(qpg(3)) ; qpg_xy=SQRT(qpg(1)**2+qpg(2)**2)

       if (qpg_z > SMALL) then
         ! Analytic expression.
         call CALCK0(qpg_z *rcut_, k0, 1)
         call CALJY1(qpg_xy*rcut_, j1, 0)
         call CALJY0(qpg_xy*rcut_, j0, 0)
         call CALCK1(qpg_z *rcut_, k1, 1)
         vc_cut(ig) = (four_pi/qpg2)*(one+rcut_*qpg_xy*j1*k0-qpg_z*rcut_*j0*k1)
       else
         if (qpg_xy > SMALL) then
           ! Integrate r*Jo(G_xy r)log(r) from 0 up to rcut_
           call quadrature(F5,zero,rcut_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
           ABI_CHECK(ierr == 0, "Accuracy not reached")
           vc_cut(ig)=-four_pi*quad
         else
           ! Analytic expression
           vc_cut(ig)=-pi*rcut_**2*(two*LOG(rcut_)-one)
         end if
       end if
     end do ! ig
   end if !finite/infinite

 case default
   ABI_BUG(sjoin('Wrong value for method:',itoa(method)))
 end select

 ! Collect vc_cut on each core
 call xmpi_sum(vc_cut, comm, ierr)

end subroutine cutoff_cylinder
!!***

!----------------------------------------------------------------------

real(dp) function F1(rho)

 real(dp),intent(in) :: rho

!Local variables-------------------------------
 integer,parameter :: order=0,ll=0
 real(dp) :: arg,bes,besp,bespp
!************************************************************************

 ! F1(\rho;z)= \rho*j_o(qpg_perp_*\rho)/sqrt(\rho**2+z**2)
 arg=rho*qpg_perp_
 call paw_jbessel(bes,besp,bespp,ll,order,arg)

 if (zz_==zero) then
   F1=bes
 else
   F1=bes*rho/SQRT(rho**2+zz_**2)
 end if

end function F1
!!***

!----------------------------------------------------------------------

function F2(xx)

 real(dp),intent(in) :: xx
 real(dp) :: F2

!Local variables-------------------------------
!scalars
 integer :: ierr
 real(dp) :: intr
!************************************************************************

 zz_=xx
 call quadrature(F1,zero,rcut_,qopt_,intr,ierr,ntrial_,accuracy_,npts_)
 ABI_CHECK(ierr == 0, "Accuracy not reached")

 F2=intr*COS(qpg_para_*xx)

end function F2
!!***

!----------------------------------------------------------------------

real(dp) pure function F3(xx)

 real(dp),intent(in) :: xx

!************************************************************************

 ! F3(z)=z*\sin(qpg_para_*z)/\sqrt(rcut^2+z^2)
 F3=xx*SIN(qpg_para_*xx)/SQRT(rcut_**2+xx**2)

end function F3
!!***

!----------------------------------------------------------------------

real(dp) function F4(rho)

 real(dp),intent(in) :: rho

!Local variables-------------------------------
 integer,parameter :: order=0,ll=0
 real(dp) :: arg,bes,besp,bespp
!************************************************************************

 ! $F4(rho)=\rho*j_o(qpg_perp_.\rho) \ln((hcyl+\sqrt(rho^2+hcyl^2))/\rho)$
 if (ABS(rho)<tol12) then
   F4=zero
 else
   arg=rho*qpg_perp_
   call paw_jbessel(bes,besp,bespp,ll,order,arg)
   F4=bes*rho*LOG((hcyl_+SQRT(rho**2+hcyl_**2))/rho)
 end if

end function F4
!!***

!----------------------------------------------------------------------

real(dp) function F5(rho)

 real(dp),intent(in) :: rho

!Local variables-------------------------------
 integer,parameter :: order = 0, ll = 0
 real(dp) :: arg,bes,besp,bespp
!************************************************************************

 ! $F5(\rho)=\rho*j_o(G_perp\rho)log(\rho)$
 if (rho==0) then
   F5=zero
 else
   arg=rho*qpg_perp_
   call paw_jbessel(bes,besp,bespp,ll,order,arg)
   F5=bes*rho*LOG(rho)
 end if

end function F5
!!***

!----------------------------------------------------------------------

real(dp) function K0cos(yy)

 real(dp),intent(in) :: yy

!Local variables-------------------------------
 real(dp) :: k0,rho,arg
!************************************************************************

 ! K0cos(y)=K0(\rho*|qpg_z|)*COS(x.qpg_x+y*qpg_y)
 rho=SQRT(xx_**2+yy**2) ; arg=qpg_para_*rho
 call CALCK0(arg,k0,1)
 K0cos=k0*COS(qpgx_*xx_+qpgy_*yy)

end function K0cos
!!***

!----------------------------------------------------------------------

real(dp) function K0cos_dy(xx)

 real(dp),intent(in) :: xx

!Local variables-------------------------------
 integer :: ierr
 real(dp) :: quad
!************************************************************************

 !! K0cos_dy(x)=\int_{-b/2}^{b/2} K0(|qpg_z|\rho)cos(x.qpg_x+y.qpg_y)dy$
 xx_=xx
 call quadrature(K0cos,-hb_,+hb_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
 ABI_CHECK(ierr == 0, "Accuracy not reached")

 K0cos_dy=quad

end function K0cos_dy
!!***




!----------------------------------------------------------------------

real(dp) function K0cos_dy_r0(xx)

 real(dp),intent(in) :: xx

!Local variables-------------------------------
!scalars
 integer :: ierr
 real(dp) :: quad,yx
!************************************************************************

 ! $ K0cos_dy_r0(x)= \int_{-b/2}^{-y(x)} K0(|qpg_z|\rho) cos(x.qpg_x+y.qpg_y)dy
 !                  +\int_{y(x)}^{b/2} K0(|qpg_z|\rho)cos(x.qpg_x+y.qpg_y)dy$
 ! where y(x)=SQRT(r0^2-x^2) and x<=r0
 !
 xx_=xx; yx=SQRT(r0_**2-xx**2)
 call quadrature(K0cos,-hb_,-yx,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
 ABI_CHECK(ierr == 0, "Accuracy not reached in quadrature")
 K0cos_dy_r0=quad

 call quadrature(K0cos,+yx,+hb_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
 ABI_CHECK(ierr == 0, "Accuracy not reached in quadrature")

 K0cos_dy_r0=quad+K0cos_dy_r0

end function K0cos_dy_r0
!!***

!----------------------------------------------------------------------

real(dp) function K0cos_dth_r0(rho)

 real(dp),intent(in) :: rho

!Local variables-------------------------------
!scalars
 integer :: ierr
 real(dp) :: quad,arg,k0,tmp

!************************************************************************

 ! $ K0cos_dth_r0(\rho)=
 ! \int_{0}^{2pi)} K0(|qpg_z|\rho)cos(\rho.cos(\theta).qpg_x+\rho.sin(\theta).qpg_y) d\theta $
 !
 ! where y(x)=SQRT(r0^2-x^2) and x<=r0
 !
 rho_=rho
 call quadrature(Fcos_th,zero,two_pi,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
 ABI_CHECK(ierr == 0, "Accuracy not reached in quadrature")

 arg=qpg_para_*rho_
 tmp=zero
 if (arg>tol6) then
   call CALCK0(arg,k0,1)
   tmp=k0*rho_
 end if
 K0cos_dth_r0=quad*tmp

end function K0cos_dth_r0
!!***

!----------------------------------------------------------------------

pure real(dp) function Fcos_th(theta)

 real(dp),intent(in) :: theta

!************************************************************************

 ! $ Fcos_th(\theta)=rho*K0(\rho*|qpg_z|)*COS(\rho.COS(\theta).qpg_x+\rho.SIN/(\theta)*qpg_y) $

 !arg=qpg_para_*rho_
 !call CALCK0(arg,k0,1)
 !tmp=k0*rho_
 Fcos_th=COS(rho_*COS(theta)*qpgx_+rho_*SIN(theta)*qpgy_)

end function Fcos_th
!!***

!----------------------------------------------------------------------

!the following functions should be used to deal with the singularity in the Cylindrical cutoff
!TODO Not yet used and indeed are still private

function K0fit(mq,nn) result(vals)

 integer,intent(in) :: nn
 real(dp),intent(in) :: mq
 real(dp) :: vals(nn)

!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: mqh
!arrays
 real(dp),parameter :: cc(7)=(/-0.57721566,0.42278420,0.23069756, &
                                0.03488590,0.00262698,0.00010750,0.00000740/)
 ! *************************************************************************

 if (nn>8.or.nn<1) then
   ABI_ERROR("nn>8.or.nn<1 not implemented")
 end if

 ! === Eq 9.8.5 in Abramovitz ===
 vals(1)=-LOG(mq*half)*I0(mq)
 mqh=mq*half
 do ii=2,nn
   vals(ii)=cc(ii-1)*mqh**(2*(ii-2))
 end do

end function K0fit

real(dp) function K0fit_int(mq,par,nn) result(integ)

 integer,intent(in) :: nn
 real(dp),intent(in) :: mq
 real(dp),intent(in) :: par(nn)

!Local variables-------------------------------
!scalars
 integer :: ii,aa
 real(dp) :: mqh
!arrays
 real(dp),parameter :: cc(7)=(/-0.57721566,0.42278420,0.23069756,&
&                               0.03488590,0.00262698,0.00010750,0.00000740/)
 ! *************************************************************************

 if (nn>8.or.nn<1) then
   ABI_ERROR("nn>8.or.nn<1 not implemented")
 end if

 mqh=mq*half
 integ=-par(1)*int_I0ln(mqh)
 ! primitive of polynomial \sum_0^{N/2} cc_{2i} (x/2)^{2*i}
 do ii=2,nn
  aa=(2*(ii-1)+1)
  integ=integ+par(ii)*two*cc(ii-1)*(mqh**aa)/aa
 end do

end function K0fit_int

real(dp) function I0(xx)

 real(dp),intent(in) :: xx

!Local variables-------------------------------
 real(dp) :: tt

! *************************************************************************

 ! Eq 9.8.1 of Abramovitz, entering the expansion of K0 -->0
 ! Expansion holds for |x|<3.75, Error<1.6*10D-07
 tt=xx/3.75
 I0=one+3.5156229*tt**2+3.0899424*tt**4 +1.2067492*tt**6 &
       +0.2659732*tt**8+0.0360768*tt**10+0.0045813*tt**12
end function I0

! Primitive of x^m Ln(x) for m/=-1
real(dp) function int_xmln(xx,mm)  result(res)

 integer,intent(in) :: mm
 real(dp),intent(in) :: xx

! *********************************************************************

 if (mm==-1) then
   ABI_BUG('invalid value for mm')
 end if

 if (xx<=zero) then
   ABI_BUG(' invalid value for xx')
 end if

 res= (xx**(mm+1))/(mm+1) * (LOG(xx) - one/(mm+1))

end function int_xmln

! Primitive function of ln(x/2)*I0(x) = sum_0^{N/2} 2^{2s+1} c_{2s} T(x/2,2s)
! where T(x,s)=\int x^s ln(x)dx
real(dp) function int_I0ln(xx) result(res)

!Arguments ------------------------------------
 real(dp),intent(in) :: xx

!Local variables-------------------------------
 real(dp) :: yy
! *********************************************************************

 yy=xx*half
 res =  (       one*2    *int_xmln(yy,0)  &
&        +3.5156229*2**3 *int_xmln(yy,2)  &
&        +3.0899424*2**5 *int_xmln(yy,4)  &
&        +1.2067492*2**7 *int_xmln(yy,6)  &
&        +0.2659732*2**9 *int_xmln(yy,8)  &
&        +0.0360768*2**11*int_xmln(yy,10) &
&        +0.0045813*2**13*int_xmln(yy,12) &
&       )

end function int_I0ln
!!***

end module m_cutoff_cylinder
!!***
