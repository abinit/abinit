!!****m* ABINIT/m_gtermcutoff
!! NAME
!!  m_gtermcutoff
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group ()
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

module m_gtermcutoff

 use defs_basis
 use m_dtset
 use m_errors
 use m_xmpi
 use m_fstrings,        only : sjoin, itoa
 use m_profiling_abi,   only : abimem_record
 use defs_abitypes,     only : MPI_type
 use m_bessel,          only : CALCK0
 use m_numeric_tools,   only : arth, l2norm, OPERATOR(.x.),quadrature

 use m_geometry,        only : normv, metric

 use m_crystal,         only : crystal_t
 use m_gsphere,         only : gsphere_t
 use m_bz_mesh,         only : kmesh_t,kmesh_init

 implicit none

 private
!!***

 public :: termcutoff
!!***
! private variables used for the integration needed by the cylindrical case.
 integer,save  :: npts_,ntrial_,qopt_
 real(dp),save :: ha_,hb_,r0_
 real(dp),save :: gcart_para_,gcartx_,gcarty_
 real(dp),save :: xx_
 real(dp),save :: accuracy_
  
 
contains
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/termcutoff
!! NAME
!! termcutoff
!!
!! FUNCTION
!! 
!! INPUTS
!! 
!! OUTPUT
!!
!! NOTES
!!  In order to incur minimal changes in some portions of the code 
!!  where a cut-off is needed to be applied, one can work only with 
!! the cut-off part of the Coulomb potential.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine termcutoff(gmet,gprimd,nfft,ngfft,gsqcut,ucvol,gcutoff)
 
!Arguments ------------------------------------
!scalars
 integer,intent(in)    :: nfft,ngfft(18)
 real(dp),intent(in)   :: gsqcut
 real(dp),intent(in)   :: ucvol

!arrays
 real(dp),intent(in):: gmet(3,3),gprimd(3,3)

!Local variables-------------------------------
!scalars
 integer,parameter  :: N0=1000
 integer            :: i1,i2,i23,i3,ierr
 integer            :: ii,ig,ing,n1,n2,n3,id(3)
 integer            :: test,opt_surface !opt_cylinder
 real(dp)           :: cutoff,rcut,check
 real(dp)           :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
 real(dp)           :: gqg2p3,gqgm12,gqgm13,gqgm23,gs2,gs3
 real(dp)           :: gcart2,gcart_para,gcart_perp
 real(dp)           :: quad,tmp
 real(dp)           :: pdir(3),vcutgeo(3),alpha(3),rmet(3,3)
 real(dp),parameter :: tolfix=1.0000001_dp
 character(len=50)  :: mode
 character(len=500) :: msg
 type(dataset_type) :: dtset
 type(kmesh_t)      :: Kmesh 
 type(gsphere_t)    :: Gsph
 type(crystal_t)    :: Cryst
!arrays

 real(dp),allocatable :: gq(:,:),gpq(:),gpq2(:),gcart(:)
 real(dp),allocatable,intent(out) :: gcutoff(:)

!Initialize a few quantities
 cutoff=gsqcut*tolfix
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
! call metric(gmet,gprimd,-1,rmet,Cryst%rprimd,ucvol)
 
 ! Initialize container
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3))) 
 ABI_ALLOCATE(gpq,(n1*n2*n3))
 ABI_ALLOCATE(gpq2,(n1*n2*n3))
 ABI_ALLOCATE(gcart,(n1*n2*n3))  
 ABI_ALLOCATE(gcutoff,(n1*n2*n3))
 gpq = zero ; gpq2 = zero ; gcutoff=zero

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...
 
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     gq(ii,ing)=ing-(ing/id(ii))*ngfft(ii)-1
   end do
 end do

 ! Get the cut-off method info from the input file
 ! Assign method to one of the available cases
 mode='NONE' 
 if (dtset%icutcoul==0) mode='SPHERE'
 if (dtset%icutcoul==1) mode='CYLINDER'
 if (dtset%icutcoul==2) mode='SURFACE'
 if (dtset%icutcoul==3) mode='CRYSTAL'
 if (dtset%icutcoul==4) mode='ERF'
 if (dtset%icutcoul==5) mode='ERFC'

  do i3=1,n3
   ! Precompute some products that do not depend on i2 and i1
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2

   do i2=1,n2
     i23=n1*(i2-1 + n2*(i3-1))
     gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
     gqgm12=gq(2,i2)*gmet(1,2)*2
     gqg2p3=gqgm13+gqgm12
     do i1=1,n1
        ii=i1+i23
        gpq(ii)=gs2+gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
        if(gpq(ii)>=tol4) then 
          gpq2(ii) = piinv/gpq(ii)
        end if 
     end do
   end do
 end do

 !Initialize geomtry type to help select CASE
 vcutgeo=dtset%vcutgeo 
 
 SELECT CASE (TRIM(mode))

   CASE('SPHERE') ! Spencer-Alavi method

     ! Calculate rcut for each method 
     rcut= (three*ucvol/four_pi)**(one/three)
     do ig=1,nfft
       if(abs(gpq(ig))<tol4) then
          gcutoff(ig)=two_pi*rcut 
       else if(gpq(ig)<=cutoff) then
          gcutoff(ig)=one-cos(rcut*sqrt(four_pi/gpq2(ig)))
      end if
     end do

   CASE('CYLINDER')

     test=COUNT(ABS(vcutgeo)>tol6)
     ABI_CHECK(test==1,'Wrong cutgeo for cylinder')   

     !Calculate rcut for each method !
     gcutoff(:)=1 ! Neutral cut-off

     ! * Check if Bravais lattice is orthorombic and parallel to the Cartesian versors.
     !   In this case the intersection of the W-S cell with the x-y plane is a rectangle with -ha_<=x<=ha_ and -hb_<=y<=hb_
     if ( (ANY(ABS(Cryst%rprimd(2:3,  1))>tol6)).or.&
&         (ANY(ABS(Cryst%rprimd(1:3:2,2))>tol6)).or.&
&         (ANY(ABS(Cryst%rprimd(1:2,  3))>tol6))    &
&       ) then
       msg = ' Bravais lattice should be orthorombic and parallel to the cartesian versors '
       MSG_ERROR(msg)
     end if

     ha_=half*SQRT(DOT_PRODUCT(Cryst%rprimd(:,1),Cryst%rprimd(:,1)))
     hb_=half*SQRT(DOT_PRODUCT(Cryst%rprimd(:,2),Cryst%rprimd(:,2)))
     r0_=MIN(ha_,hb_)/N0
     !
     ! ===================================================
     ! === Setup for the quadrature of matrix elements ===
     ! ===================================================
     qopt_    =6         ! Quadrature method, see quadrature routine.
     ntrial_  =30        ! Max number of attempts.
     accuracy_=0.001     ! Fractional accuracy required.
     npts_    =6         ! Initial number of point (only for Gauss-Legendre method).

     do ig=1,nfft

       gcart(:)=b1(:)*Gsph%gvec(1,ig)+b2(:)*Gsph%gvec(2,ig)+b3(:)*Gsph%gvec(3,ig)
       gcartx_=gcart(1) ; gcarty_=gcart(2) ; gcart_para_=ABS(gcart(3))

       tmp=zero

       call quadrature(K0cos_dy,zero,ha_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)

       if (ierr/=0) then
         MSG_ERROR("Accuracy not reached")
       end if
       tmp=tmp+quad
       gcutoff(ig)=two*(tmp*two)
     end do !ig


   CASE('SURFACE')

     test=COUNT(vcutgeo/=zero)
     ABI_CHECK(test==2,"Wrong vcutgeo")

     ! === From reduced to cartesian coordinates ===
     a1=Cryst%rprimd(:,1); b1=two_pi*gprimd(:,1)
     a2=Cryst%rprimd(:,2); b2=two_pi*gprimd(:,2)
     a3=Cryst%rprimd(:,3); b3=two_pi*gprimd(:,3)

     ! Calculate rcut for each method !
     rcut = half*SQRT(DOT_PRODUCT(a3,a3))

     !SURFACE Default - Beigi
     opt_surface=1; alpha(:)=zero
     ! Otherwsise use Rozzi's method
     if (ANY(vcutgeo<zero)) opt_surface=2
     pdir(:)=zero
     do ii=1,3
       check=vcutgeo(ii)
       if (ABS(check)>zero) then ! Use Rozzi"s method with a finite surface along x-y
         pdir(ii)=1
         if (check<zero) alpha(ii)=normv(check*Cryst%rprimd(:,ii),rmet,'R')
       end if
     end do

     SELECT CASE (opt_surface)

     !CASE SURFACE 1 - Beigi
     CASE(1)
   
     do ig=1,nfft
       gcart(:)=b1(:)*Gsph%gvec(1,ig)+b2(:)*Gsph%gvec(2,ig)+b3(:)*Gsph%gvec(3,ig)
       gcart2=DOT_PRODUCT(gcart(:),gcart(:))
       gcart_para=SQRT(gcart(1)**2+gcart(2)**2) ; gcart_perp = gcart(3)  
       gcutoff(ig)=one-EXP(-gcart_para*rcut)*COS(gcart_perp*rcut)
     end do !ig
    
     !CASE SURFACE 2 - Rozzi
     CASE(2)

     !!BG: Trigger needed - use the available input value for this 
     do ig=1,nfft
       gcart(:)=b1(:)*Gsph%gvec(1,ig)+b2(:)*Gsph%gvec(2,ig)+b3(:)*Gsph%gvec(3,ig)
       gcart2=DOT_PRODUCT(gcart(:),gcart(:))
       gcart_para=SQRT(gcart(1)**2+gcart(2)**2) ; gcart_perp = gcart(3)
       if (gcart_para>tol4) then
         gcutoff(ig)=one+EXP(-gcart_para*rcut)*(gcart_perp/gcart_para*SIN(gcart_perp*rcut)-COS(gcart_perp*rcut))
       else
       if (ABS(gcart_perp)>tol4) then
         gcutoff(ig)=one-COS(gcart_perp*rcut)-gcart_perp*rcut*SIN(gcart_perp*rcut)
       else
         gcutoff(ig)=-two_pi*rcut**2
       end if
      end if
     end do !ig
   
     CASE DEFAULT
       write(msg,'(a,i3)')' Wrong value of surface method: ',opt_surface
       MSG_BUG(msg)
     END SELECT

   CASE('ERF')

   ! Calculate rcut for each method ! Same as SPHERE
   rcut= (three*Kmesh%nbz*ucvol/four_pi)**(one/three)
  
     do ig=1,nfft
       if(abs(gpq(ig))<tol4) then
          gcutoff(ig)=zero ! @Gamma: initialize quantity in each requiered routine
       else if(gpq(ig)<=cutoff) then
          gcutoff(ig)=exp(-pi/(gpq2(ig)*rcut**2))
       end if
     end do
 
   CASE('ERFC')

   ! Calculate rcut for each method ! Same as SPHERE
   rcut= (three*Kmesh%nbz*ucvol/four_pi)**(one/three)

     do ig=1,nfft
       if(abs(gpq(ig))<tol4) then
          gcutoff(ig)=zero ! @Gamma: initialize quantity in each requiered routine
       else if(gpq(ig)<=cutoff) then
          gcutoff(ig)=one-exp(-pi/(gpq2(ig)*rcut**2))
       end if
     end do

   CASE('CRYSTAL')
     gcutoff(:)=one ! Neutral cut-off
     write(msg,'(a)')'CRYSTAL method: no cut-off applied to G**2 while CRYSTAL method is implied!'
     MSG_WARNING(msg) 
   CASE DEFAULT
     gcutoff=one ! Neutral cut-off
     write(msg,'(a)')'No cut-off applied to G**2!'
     MSG_WARNING(msg)
 END SELECT

 !write(*,*)'This is mode', mode
 !write(*,*)'This is mode', gcutoff

 ABI_DEALLOCATE(gq) 
 ABI_DEALLOCATE(gpq)
 ABI_DEALLOCATE(gpq2)
 ABI_DEALLOCATE(gcart)
! ABI_DEALLOCATE(gcutoff)
 
end subroutine termcutoff 
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
 rho=SQRT(xx_**2+yy**2) ; arg=gcart_para_*rho
 call CALCK0(arg,k0,1)
 K0cos=k0*COS(gcartx_*xx_+gcarty_*yy)

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

end module m_gtermcutoff
!!***
