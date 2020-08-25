!!****m* ABINIT/m_vcoul/m_cutoff_cylinder
!! NAME
!!  m_cutoff_cylinder
!!
!! FUNCTION
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

module m_cutoff_cylinder

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_splines
 use m_sort

 use m_fstrings,        only : sjoin, itoa
 use m_geometry,        only : normv, metric
 use m_numeric_tools,   only : quadrature
 use m_bessel,          only : CALJY0, CALJY1, CALCK0, CALCK1
 use m_numeric_tools,   only : arth, geop, imin_loc, llsfit_svd, l2norm, OPERATOR(.x.), quadrature, isdiagmat
 use m_paw_numeric,     only : paw_jbessel

 implicit none

 private
!!***

 public :: cutoff_cylinder, K0cos
!!***

! private variables used for the integration needed by the cylindrical case.
 integer,save :: npts_,ntrial_,qopt_
 real(dp),save :: ha_,hb_,r0_
 real(dp),save :: qpg_perp_,qpg_para_,qpgx_,qpgy_
 real(dp),save :: zz_,xx_
 real(dp),save :: hcyl_,rcut_,accuracy_

CONTAINS  
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/cutoff_cylinder
!! NAME
!! cutoff_cylinder
!!
!! FUNCTION
!!  Calculate the Fourier components of an effective Coulomb interaction
!!  zeroed outside a finite cylindrical region.
!!  Two methods are implemented:
!!   method==1: The interaction in the (say) x-y plane is truncated outside the Wigner-Seitz
!!              cell centered on the wire in the x-y plane. The interaction has infinite
!!              extent along the z axis and the Fourier transform is singular only at the Gamma point.
!!              Only orthorombic Bravais lattice are supported.
!!   method==2: The interaction is truncated outside a cylinder of radius rcut. The cylinder has finite
!!              extent along z. No singularity occurs.
!!
!! INPUTS
!!  boxcenter(3)= center of the wire in the x-y axis
!!  gvec(3,ng)=G vectors in reduced coordinates
!!  ng=number of G vectors
!!  qpt(3,nq)= q-points
!!  nq=number of q-points
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!   where: rprimd(i,j)=rprim(i,j)*acell(j)
!!  method=1 for Beigi approach (infinite cylinder with interaction truncated outside the W-S cell)
!!         2 for Rozzi method (finite cylinder)
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  vccut(ng,nq)= Fourier components of the effective Coulomb interaction
!!
!! PARENTS
!!      m_barevcoul,m_vcoul
!!
!! CHILDREN
!!      calck0,calck1,caljy0,caljy1,metric,paw_jbessel,quadrature,wrtout
!!      xmpi_split_work,xmpi_sum
!!
!! SOURCE

subroutine cutoff_cylinder(nq,qpt,ng,gvec,rcut,hcyl,pdir,boxcenter,rprimd,vccut,method,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng,nq,method,comm
 real(dp),intent(in) :: rcut,hcyl
!arrays
 integer,intent(in) :: gvec(3,ng),pdir(3)
 real(dp),intent(in) :: boxcenter(3),qpt(3,nq),rprimd(3,3)
 real(dp),intent(out) :: vccut(ng,nq)

!Local variables-------------------------------
!scalars
 integer,parameter :: N0=1000
 integer :: icount,ig,igs,iq
 integer :: ntasks,ierr,my_start,my_stop
 real(dp) :: j0,j1,k0,k1,qpg2,qpg_xy,tmp
 real(dp) :: qpg_z,quad,rcut2,hcyl2,c1,c2,ucvol,SMALL
 logical :: q_is_zero
 character(len=500) :: msg
!arrays
 real(dp) :: qpg(3),b1(3),b2(3),b3(3),gmet(3,3),rmet(3,3),gprimd(3,3),qc(3),gcart(3)
 real(dp),allocatable :: qcart(:,:)
!************************************************************************

 ABI_UNUSED(pdir)
 ABI_UNUSED(boxcenter)
 !
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

 write(msg,'(3a,2(a,i5,a),a,f8.5)')ch10,&
&  ' cutoff_cylinder: Info on the quadrature method : ',ch10,&
&  '  Quadrature scheme      = ',qopt_,ch10,&
&  '  Max number of attempts = ',ntrial_,ch10,&
&  '  Fractional accuracy    = ',accuracy_
 call wrtout(std_out,msg,'COLL')
 !
 ! === From reduced to Cartesian coordinates ===
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 b1(:)=two_pi*gprimd(:,1)
 b2(:)=two_pi*gprimd(:,2)
 b3(:)=two_pi*gprimd(:,3)

 ABI_MALLOC(qcart,(3,nq))
 do iq=1,nq
   qcart(:,iq)=b1(:)*qpt(1,iq)+b2(:)*qpt(2,iq)+b3(:)*qpt(3,iq)
 end do

 ntasks=nq*ng
 call xmpi_split_work(ntasks,comm,my_start,my_stop)
 !
 ! ================================================
 ! === Different approaches according to method ===
 ! ================================================
 vccut(:,:)=zero

 SELECT CASE (method)

 CASE (1)
   ! === Infinite cylinder, interaction is zeroed outside the Wigner-Seitz cell ===
   ! * Beigi"s expression holds only if the BZ is sampled only along z.
   write(msg,'(2(a,f8.4))')' cutoff_cylinder: Using Beigi''s Infinite cylinder '
   call wrtout(std_out,msg,'COLL')
   if (ANY(qcart(1:2,:)>SMALL)) then
     write(std_out,*)' qcart = ',qcart(:,:)
     write(msg,'(5a)')&
&      ' found q-points with non zero components in the X-Y plane. ',ch10,&
&      ' This is not allowed, see Notes in cutoff_cylinder.F90. ',ch10,&
&      ' ACTION: Modify the q-point sampling. '
     MSG_ERROR(msg)
   end if
   ! * Check if Bravais lattice is orthorombic and parallel to the Cartesian versors.
   !   In this case the intersection of the W-S cell with the x-y plane is a rectangle with -ha_<=x<=ha_ and -hb_<=y<=hb_
   if ( (ANY(ABS(rprimd(2:3,  1))>tol6)).or.&
&       (ANY(ABS(rprimd(1:3:2,2))>tol6)).or.&
&       (ANY(ABS(rprimd(1:2,  3))>tol6))    &
&     ) then
     msg = ' Bravais lattice should be orthorombic and parallel to the cartesian versors '
     MSG_ERROR(msg)
   end if

   ha_=half*SQRT(DOT_PRODUCT(rprimd(:,1),rprimd(:,1)))
   hb_=half*SQRT(DOT_PRODUCT(rprimd(:,2),rprimd(:,2)))
   r0_=MIN(ha_,hb_)/N0
   !
   ! === For each pair (q,G) evaluate the integral defining the cutoff Coulomb  ===
   ! * Here the code assumes that all q-vectors are non zero and q_xy/=0.
   do iq=1,nq
     igs=1
     ! * Skip singularity at Gamma, it will be treated "by hand" in csigme.
     q_is_zero = (normv(qpt(:,iq),gmet,'G')<tol4)
     !if (q_is_zero) igs=2
     qc(:)=qcart(:,iq)
     write(msg,'(2(a,i3))')' entering loop iq: ',iq,' with igs = ',igs
     call wrtout(std_out,msg,'COLL')

     do ig=igs,ng
       icount=ig+(iq-1)*ng; if (icount<my_start.or.icount>my_stop) CYCLE

       gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
       qpg(:)=qc(:)+gcart(:)
       qpgx_=qpg(1) ; qpgy_=qpg(2) ; qpg_para_=ABS(qpg(3))
       write(std_out,*)"qpgx_=",qpgx_
       write(std_out,*)"qpgy_=",qpgy_
       write(std_out,*)"qpg_para=",qpg_para_

       ! Avoid singularity in K_0{qpg_para_\rho) by using a small q along the periodic dimension.
       if (q_is_zero.and.qpg_para_<tol6) then
         qpg_para_ = tol6
         write(std_out,*)"setting qpg_para to=",qpg_para_
       end if
       !
       ! * Calculate $ 2\int_{WS} dxdy K_0{qpg_para_\rho) cos(x.qpg_x + y.qpg_y) $
       !   where WS is the Wigner-Seitz cell.
       tmp=zero
       ! Difficult part, integrate on a small cirle of radius r0 using spherical coordinates
       !call quadrature(K0cos_dth_r0,zero,r0_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
       !if (ierr/=0) MSG_ERROR("Accuracy not reached")
       !write(std_out,'(i8,a,es14.6)')ig,' 1 ',quad
       !tmp=tmp+quad
       ! Add region with 0<=x<=r0 and y>=+-(SQRT(r0^2-x^2))since WS is rectangular
       !call quadrature(K0cos_dy_r0,zero,r0_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
       !if (ierr/=0) MSG_ERROR("Accuracy not reached")
       !write(std_out,'(i8,a,es14.6)')ig,' 2 ',quad
       !tmp=tmp+quad
       ! Get the in integral in the rectangle with x>=r0, should be the easiest but sometimes has problems to converge
       !call quadrature(K0cos_dy,r0_,ha_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
       !if (ierr/=0) MSG_ERROR("Accuracy not reached")
       !write(std_out,'(i8,a,es14.6)')ig,' 3 ',quad
       !
       ! === More stable method: midpoint integration with Romberg extrapolation ===
       call quadrature(K0cos_dy,zero,ha_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
       !write(std_out,'(i8,a,es14.6)')ig,' 3 ',quad
       if (ierr/=0) then
         MSG_ERROR("Accuracy not reached")
       end if
       ! === Store final result ===
       ! * Factor two comes from the replacement WS -> (1,4) quadrant thanks to symmetries of the integrad.
       tmp=tmp+quad
       vccut(ig,iq)=two*(tmp*two)
     end do !ig
   end do !iq

 CASE (2)
   ! === Finite cylinder of length hcyl, from Rozzi et al ===
   ! TODO add check on hcyl value that should be smaller that 1/deltaq
   if (hcyl_<zero) then
     write(msg,'(a,f8.4)')' Negative value for cylinder length hcyl_=',hcyl_
     MSG_BUG(msg)
   end if

   if (ABS(hcyl_)>tol12) then
     write(msg,'(2(a,f8.4))')' cutoff_cylinder: using finite cylinder of length= ',hcyl_,' rcut= ',rcut_
     call wrtout(std_out,msg,'COLL')
     hcyl2=hcyl_**2
     rcut2=rcut_**2

     do iq=1,nq
       write(msg,'(a,i3)')' entering loop iq: ',iq
       call wrtout(std_out,msg,'COLL')
       qc(:)=qcart(:,iq)
       do ig=1,ng
         ! === No singularity occurs in finite cylinder, thus start from 1 ===
         icount=ig+(iq-1)*ng ; if (icount<my_start.or.icount>my_stop) CYCLE
         gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
         qpg(:)=qc(:)+gcart(:)
         qpg_para_=ABS(qpg(3)) ; qpg_perp_=SQRT(qpg(1)**2+qpg(2)**2)

         if (qpg_perp_/=zero.and.qpg_para_/=zero) then
           ! $ 4\pi\int_0^{R_c} d\rho\rho j_o(qpg_perp_.\rho)\int_0^hcyl dz\cos(qpg_para_*z)/sqrt(\rho^2+z^2) $
           call quadrature(F2,zero,rcut_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
           if (ierr/=0) then
             MSG_ERROR("Accuracy not reached")
           end if

           vccut(ig,iq)=four_pi*quad

         else if (qpg_perp_==zero.and.qpg_para_/=zero) then
           ! $ \int_0^h sin(qpg_para_.z)/\sqrt(rcut^2+z^2)dz $
           call quadrature(F3,zero,hcyl_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
           if (ierr/=0) then
             MSG_ERROR("Accuracy not reached")
           end if

           c1=one/qpg_para_**2-COS(qpg_para_*hcyl_)/qpg_para_**2-hcyl_*SIN(qpg_para_*hcyl_)/qpg_para_
           c2=SIN(qpg_para_*hcyl_)*SQRT(hcyl2+rcut2)
           vccut(ig,iq)=four_pi*c1+four_pi*(c2-quad)/qpg_para_

         else if (qpg_perp_/=zero.and.qpg_para_==zero) then
           ! $ 4pi\int_0^rcut d\rho \rho J_o(qpg_perp_.\rho) ln((h+\sqrt(h^2+\rho^2))/\rho) $
           call quadrature(F4,zero,rcut_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
           if (ierr/=0) then
             MSG_ERROR("Accuracy not reached")
           end if

           vccut(ig,iq)=four_pi*quad

         else if (qpg_perp_==zero.and.qpg_para_==zero) then
           ! Use lim q+G --> 0
           vccut(ig,iq)=two_pi*(-hcyl2+hcyl_*SQRT(hcyl2+rcut2)+rcut2*LOG((hcyl_+SQRT(hcyl_+SQRT(hcyl2+rcut2)))/rcut_))

         else
           MSG_BUG('You should not be here!')
         end if

       end do !ig
     end do !iq

   else
     ! === Infinite cylinder ===
     call wrtout(std_out,' cutoff_cylinder: using Rozzi''s method with infinite cylinder ','COLL')
     do iq=1,nq
       write(msg,'(a,i3)')' entering loop iq: ',iq
       call wrtout(std_out,msg,'COLL')
       qc(:)=qcart(:,iq)
       do ig=1,ng
         icount=ig+(iq-1)*ng ; if (icount<my_start.or.icount>my_stop) CYCLE
         gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
         qpg(:)=qc(:)+gcart(:)
         qpg2  =DOT_PRODUCT(qpg,qpg)
         qpg_z =ABS(qpg(3)) ; qpg_xy=SQRT(qpg(1)**2+qpg(2)**2)
         if (qpg_z>SMALL) then
           ! === Analytic expression ===
           call CALCK0(qpg_z *rcut_,k0,1)
           call CALJY1(qpg_xy*rcut_,j1,0)
           call CALJY0(qpg_xy*rcut_,j0,0)
           call CALCK1(qpg_z *rcut_,k1,1)
           vccut(iq,ig)=(four_pi/qpg2)*(one+rcut_*qpg_xy*j1*k0-qpg_z*rcut_*j0*k1)
         else
           if (qpg_xy>SMALL) then
             ! === Integrate r*Jo(G_xy r)log(r) from 0 up to rcut_  ===
             call quadrature(F5,zero,rcut_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
             if (ierr/=0) then
               MSG_ERROR("Accuracy not reached")
             end if
             vccut(ig,iq)=-four_pi*quad
           else
             ! === Analytic expression ===
             vccut(ig,iq)=-pi*rcut_**2*(two*LOG(rcut_)-one)
           end if
         end if
       end do !ig
     end do !iq
   end if !finite/infinite

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong value for method:',itoa(method)))
 END SELECT
 !
 ! === Collect vccut on each node ===
 call xmpi_sum(vccut,comm,ierr)

 ABI_FREE(qcart)

end subroutine cutoff_cylinder
!!***

!----------------------------------------------------------------------

function F1(rho)

 real(dp),intent(in) :: rho
 real(dp) :: F1

!Local variables-------------------------------
!scalars
 integer,parameter :: order=0,ll=0
 real(dp) :: arg,bes,besp,bespp
!************************************************************************

 !F1(\rho;z)= \rho*j_o(qpg_perp_*\rho)/sqrt(\rho**2+z**2)
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
 if (ierr/=0) then
   MSG_ERROR("Accuracy not reached")
 end if

 F2=intr*COS(qpg_para_*xx)

end function F2
!!***

!----------------------------------------------------------------------

pure function F3(xx)

 real(dp),intent(in) :: xx
 real(dp) :: F3
!************************************************************************

 ! F3(z)=z*\sin(qpg_para_*z)/\sqrt(rcut^2+z^2)
 F3=xx*SIN(qpg_para_*xx)/SQRT(rcut_**2+xx**2)

end function F3
!!***

!----------------------------------------------------------------------

function F4(rho)

 real(dp),intent(in) :: rho
 real(dp) :: F4

!Local variables-------------------------------
!scalars
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

function F5(rho)

 real(dp),intent(in) :: rho
 real(dp) :: F5

!Local variables-------------------------------
!scalars
 integer,parameter :: order=0,ll=0
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

function K0cos(yy)

 real(dp),intent(in) :: yy
 real(dp) :: K0cos

!Local variables-------------------------------
!scalars
 real(dp) :: k0,rho,arg
!************************************************************************

 ! K0cos(y)=K0(\rho*|qpg_z|)*COS(x.qpg_x+y*qpg_y)
 rho=SQRT(xx_**2+yy**2) ; arg=qpg_para_*rho
 call CALCK0(arg,k0,1)
 K0cos=k0*COS(qpgx_*xx_+qpgy_*yy)

end function K0cos
!!***

!----------------------------------------------------------------------

function K0cos_dy(xx)

 real(dp),intent(in) :: xx
 real(dp) :: K0cos_dy
!Local variables-------------------------------
!scalars
 integer :: ierr
 real(dp) :: quad
!************************************************************************

 !! K0cos_dy(x)=\int_{-b/2}^{b/2} K0(|qpg_z|\rho)cos(x.qpg_x+y.qpg_y)dy$
 xx_=xx
 call quadrature(K0cos,-hb_,+hb_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
 if (ierr/=0) then
   MSG_ERROR("Accuracy not reached")
 end if

 K0cos_dy=quad

end function K0cos_dy
!!***

end module m_cutoff_cylinder
!!***
