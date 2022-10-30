!!****m* ABINIT/m_vcoul/m_cutoff_surface
!! NAME
!!  m_cutoff_surface
!!
!! FUNCTION
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_cutoff_surface

 use defs_basis
 use m_abicore
 use m_errors

 use m_fstrings, only : sjoin, itoa

 implicit none

 private
!!***

 public :: cutoff_surface

 !integer,public,parameter :: SURFACE_BEIGI = 1
 !integer,public,parameter :: SURFACE_ROZZI = 2

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/cutoff_surface
!! NAME
!! cutoff_surface
!!
!! FUNCTION
!!  Calculate the Fourier components of an effective Coulomb interaction
!!  within a slab of thickness 2*rcut which is symmetric with respect to the xy plane.
!!  In this implementation rcut=L_z/2 where L_z is the periodicity along z
!!
!! INPUTS
!!  qpt(3)=q-point
!!  ng=Number of G vectors.
!!  gvec(3,ng)=G vectors in reduced coordinates.
!!  gprimd(3,3)=Dimensional primitive translations in reciprocal space ($\textrm{bohr}^{-1}$).
!!  gmet(3,3)=Metric in reciprocal space.
!!
!! OUTPUT
!!  vc_cut(ng)=Fourier components of the effective Coulomb interaction.
!!
!! NOTES
!!  The Fourier expression for an interaction truncated along the z-direction (i.e non-zero only if |z|<R) is:
!!
!!  vc(q.G) = 4pi/|q+G|^2 * [ 1 + e^{-((q+G)_xy)*R} * ( (q_z+G_z)/(q+G)_xy * sin((q_z+G_z)R) -
!!   - cos((q_z+G_Z)R)) ]  (1)
!!
!!  Equation (1) diverges when q_xy+G_xy --> 0 for any non zero q_z+G_z
!!  However if we choose R=L/2, where L defines the periodicity along z,
!!  and we limit ourselves to consider q-points with q_z==0, then sin((q_z+G_z)R)=sin(G_Z 2pi/L)=0 for every G.
!!  Under these assumptions we obtain
!!
!!  v(q,G) = 4pi/|q+G|^2 [1-e^{-(q+G)_xy*L/2}\cos((q_z+G_z)R)]
!!
!!  which is always finite when G_z /=0 while it diverges as 4piR/(q+G)_xy as (q+G)_xy -->0
!!  but only in the x-y plane.
!!
!! SOURCE

subroutine cutoff_surface(qpt, ng, gvec, gprimd, rcut, boxcenter, pdir, alpha, vc_cut, method)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: method,ng
 real(dp),intent(in) :: rcut
!arrays
 integer,intent(in) :: gvec(3,ng),pdir(3)
 real(dp),intent(in) :: alpha(3),boxcenter(3),gprimd(3,3),qpt(3)
 real(dp),intent(out) :: vc_cut(ng)

!Local variables-------------------------------
!scalars
 integer :: ig,igs
 real(dp),parameter :: SMALL=tol4  !@WC: was tol6
 real(dp) :: qpg2,qpg_para,qpg_perp
 character(len=500) :: msg
!arrays
 real(dp) :: b1(3),b2(3),b3(3),gcart(3),qc(3),qpg(3)

! *************************************************************************

 ABI_UNUSED(pdir)
 ABI_UNUSED(boxcenter)

 ! From reduced to cartesian coordinates.
 b1(:)=two_pi*gprimd(:,1)
 b2(:)=two_pi*gprimd(:,2)
 b3(:)=two_pi*gprimd(:,3)

 qc = b1*qpt(1) + b2*qpt(2) + b3*qpt(3)

 ! Different approaches according to method
 vc_cut = zero

 select case (method)

 case (1)
   ! Beigi's expression.
   ! q-points with non-zero component along the z-axis are not allowed if
   ! the simplified Eq.1 for the Coulomb interaction is used.
   if (ANY(ABS(qc) > SMALL)) then
     write(std_out,*)qc
     write(msg,'(5a)')&
      'Found q-points with non-zero component along non-periodic direction ',ch10,&
      'This is not allowed, see Notes in cutoff_surface.F90 ',ch10,&
      'ACTION: Modify the q-point sampling'
     ABI_ERROR(msg)
   end if

   ! Calculate truncated Coulomb interaction for a infinite surface
   ! supposing input q-points are different from zero.
   igs=1; if (SQRT(DOT_PRODUCT(qc,qc))<tol16) igs=2 ! avoid (q=0, G=0)
   do ig=igs,ng
     gcart(:) = b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
     qpg(:) = qc(:) + gcart(:)
     qpg2  = DOT_PRODUCT(qpg(:),qpg(:))
     qpg_para = SQRT(qpg(1)**2+qpg(2)**2) ; qpg_perp=qpg(3)
     vc_cut(ig) = four_pi/qpg2*(one-EXP(-qpg_para*rcut)*COS(qpg_perp*rcut))
   end do

 case (2)
   ! Rozzi's method
   ABI_ERROR("Work in progress")
   ABI_UNUSED(alpha) ! just to keep alpha as an argument
   !alpha=?? ; ap1sqrt=SQRT(one+alpha**2)
   do ig=1,ng
     gcart(:) = b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
     qpg(:) = qc(:) + gcart(:)
     qpg2  =DOT_PRODUCT(qpg(:),qpg(:))
     qpg_para=SQRT(qpg(1)**2+qpg(2)**2) ; qpg_perp =qpg(3)
     if (qpg_para>SMALL) then
      vc_cut(ig)=four_pi/qpg2*(one+EXP(-qpg_para*rcut)*(qpg_perp/qpg_para*SIN(qpg_perp*rcut)-COS(qpg_perp*rcut)))
     else
       if (ABS(qpg_perp)>SMALL) then
         vc_cut(ig)=four_pi/qpg_perp**2*(one-COS(qpg_perp*rcut)-qpg_perp*rcut*SIN(qpg_perp*rcut)) ! &
         ! contribution due to finite surface
         ! + 8*rcut*SIN(qpg_perp*rcut)/qpg_perp*LOG((alpha+ap1sqrt)*(one+ap1sqrt)/alpha)
       else
         vc_cut(ig)=-two_pi*rcut**2
       end if
     end if
   end do !ig

 case default
   ABI_BUG(sjoin('Wrong value for method:', itoa(method)))
 end select

end subroutine cutoff_surface
!!***

end module m_cutoff_surface
!!***
