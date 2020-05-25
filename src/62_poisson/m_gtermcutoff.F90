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
 use m_abicore
 use m_errors
 use m_xmpi
 use m_fstrings,        only : sjoin, itoa
 use m_profiling_abi,   only : abimem_record
 use defs_abitypes,     only : MPI_type
 use m_bessel,          only : CALCK0
 use m_numeric_tools,   only : arth, l2norm, OPERATOR(.x.),quadrature

 use m_geometry,        only : normv, metric

 use m_crystal,         only : crystal_t
 use m_gsphere,         only : gsph_free,gsphere_t, gsph_init, print_gsphere ! print might be deleted after testing
 use m_bz_mesh,         only : kmesh_t, kmesh_init

 implicit none

 private
!!***

!!****t* m_gtermcutoff/gcut_t
!! NAME
!!  gcut_t
!!
!! FUNCTION
!!
!! SOURCE

! type,public :: gcut_t

!  integer  :: nfft
!  ! Number of points in FFT grid

!  integer  :: ng
!   ! Number of G-vectors

!  real(dp) :: ucvol
!    ! Volume of the unit cell

!   ! character(len=50) :: mode
!   ! String defining the cutoff mode, possible values are: sphere,cylinder,surface,crystal

!   ! integer :: pdir(3)
!   ! 1 if the system is periodic along this direction

!   ! real(dp) :: boxcenter(3)
!   ! 1 if the point in inside the cutoff region 0 otherwise
!   ! Reduced coordinates of the center of the box (input variable)

!  real(dp) :: vcutgeo(3)
!    ! For each reduced direction gives the length of the finite system
!    ! 0 if the system is infinite along that particular direction
!    ! negative value to indicate that a finite size has to be used

!  real(dp) :: rprimd(3,3)
!    ! Lattice vectors in real space.

!    ! real(dp),allocatable :: barev(:)
!    ! gtermcuoff(nfft)
!    ! G cut-off array on the FFT grid

! end type gcut_t

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
!!   Apply a cut-off term to the 1/G**2-like terms that appears throughout
!!   the code at the ground-state level as follows: Ewald, NC-PSP, Hartee.
!!   
!! INPUTS
!!   gsqcut     = cutoff on (k+G)^2 (bohr^-2) (sphere for density and potential) (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2))
!!   icutcoul   = Information about the cut-off
!!   ngfft(18)  = Information on the (fine) FFT grid used for the density.
!!   nkpt       = Number of k-points in the Brillouin zone 
!!   rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!   vcutgeo(3)= Info on the orientation and extension of the cutoff region.
!!    
!! OUTPUT
!!   gcutoff  = Cut-off term applied to 1/G**2 terms
!!
!! NOTES
!!  1. In order to incur minimal changes in some portions of the code 
!!  where a cut-off is needed to be applied, one can work only with 
!!  the cut-off part of the Coulomb potential, unlike what is done
!!  in barevcoul module.
!!  2. Fock term has its own legacy cut-off for the moment.
!!
!! PARENTS
!!     atm2fft,mklocl
!! CHILDREN
!!     calck0,paw_jbessel,quadrature
!!
!! SOURCE

subroutine termcutoff(gcutoff,gsqcut,icutcoul,ngfft,nkpt,rprimd,vcutgeo)
 
!Arguments ------------------------------------
!scalars
 integer,intent(in)   :: icutcoul, nkpt
 real(dp),intent(in)  :: gsqcut

!arrays
 integer,intent(in)    :: ngfft(18)
 real(dp),intent(in)   :: rprimd(3,3),vcutgeo(3)

!Local variables-------------------------------
!scalars
 integer,parameter  :: N0=1000
 integer            :: i1,i2,i23,i3,ierr,id(3),ii,ig,ing
 integer            :: n1,n2,n3,nfft
 integer            :: test,opt_surface !opt_cylinder
 real(dp)           :: cutoff,rcut,check,rmet(3,3)
 real(dp)           :: gvecg2p3,gvecgm12,gvecgm13,gvecgm23,gs2,gs3
 real(dp)           :: gcart_para,gcart_perp
 real(dp)           :: quad,tmp,ucvol
 real(dp)           :: pdir(3),alpha(3)
 real(dp),parameter :: tolfix=1.0000001_dp
 character(len=50)  :: mode
 character(len=500) :: msg
! type(gcut_t)       :: gcut  !

!arrays
 real(dp)             :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
 real(dp)             :: gcart(3),gmet(3,3),gprimd(3,3)
 real(dp),allocatable :: gvec(:,:),gpq(:),gpq2(:)
 real(dp),allocatable :: gcutoff(:)

! === Save dimension and other useful quantities in vcut% ===
! gcut%nfft      = PRODUCT(ngfft(1:3))  ! Number of points in the FFT mesh.
! gcut%ucvol     = ucvol                ! Unit cell volume.
! gcut%rprimd    = rprimd(:,:)    ! Dimensional direct lattice.
! gcut%vcutgeo   = vcutgeo(:)     ! Info on the orientation and extension of the cutoff region.
!
!Initialize a few quantities
 cutoff=gsqcut*tolfix
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 nfft=n1*n2*n3
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 
 ! Initialize container
 ABI_ALLOCATE(gvec,(3,MAX(n1,n2,n3))) 
 ABI_ALLOCATE(gpq,(nfft))
 ABI_ALLOCATE(gpq2,(nfft))
 ABI_ALLOCATE(gcutoff,(nfft))
 gcart(:) = zero ; gpq = zero ; gpq2 = zero ; gcutoff = zero

 !In order to speed the routine, precompute the components of gvectors
 !Also check if the booked space was large enough...
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     gvec(ii,ing)=ing-(ing/id(ii))*ngfft(ii)-1
   end do
 end do

 ! Get the cut-off method info from the input file
 ! Assign method to one of the available cases
 mode='NONE' 

 if (icutcoul==0) mode='SPHERE'
 if (icutcoul==1) mode='CYLINDER'
 if (icutcoul==2) mode='SURFACE'
 if (icutcoul==3) mode='CRYSTAL'
 if (icutcoul==4) mode='ERF'
 if (icutcoul==5) mode='ERFC'
 
 !Print in log info about the cut-off method at every call: 
 write(msg,'(3a)')ch10,' 1/G**2 cut-off applied in the following step : cutoff-mode = ',TRIM(mode)
 call wrtout(std_out,msg)
 !!!

  do i3=1,n3
   ! Precompute some products that do not depend on i2 and i1
   gs3=gvec(3,i3)*gvec(3,i3)*gmet(3,3)
   gvecgm23=gvec(3,i3)*gmet(2,3)*2
   gvecgm13=gvec(3,i3)*gmet(1,3)*2

   do i2=1,n2
     i23=n1*(i2-1 + n2*(i3-1))
     gs2=gs3+ gvec(2,i2)*(gvec(2,i2)*gmet(2,2)+gvecgm23)
     gvecgm12=gvec(2,i2)*gmet(1,2)*2
     gvecg2p3=gvecgm13+gvecgm12
     do i1=1,n1
        ii=i1+i23
        gpq(ii)=gs2+gvec(1,i1)*(gvec(1,i1)*gmet(1,1)+gvecg2p3)
        if(gpq(ii)>=tol4) then 
          gpq2(ii) = piinv/gpq(ii)
        end if 
     end do
   end do
 end do

 SELECT CASE (TRIM(mode))

   CASE('SPHERE') ! Spencer-Alavi method

     ! Calculate rcut for each method 
     rcut= (three*ucvol/four_pi)**(one/three)
     do ig=1,nfft
       if(abs(gpq(ig))<tol4) then
          gcutoff(ig)=zero ! two_pi*rcut - value for the V_coul
       else if(gpq(ig)<=cutoff) then
          gcutoff(ig)=one-cos(rcut*sqrt(four_pi/gpq2(ig)))
      end if
     end do

   CASE('CYLINDER')

     test=COUNT(vcutgeo/=zero)
     ABI_CHECK(test==1,'Wrong cutgeo for cylinder')   
        
     !Calculate rcut for each method !
     !

     ! * Check if Bravais lattice is orthorombic and parallel to the Cartesian versors.
     !   In this case the intersection of the W-S cell with the x-y plane is a rectangle with -ha_<=x<=ha_ and -hb_<=y<=hb_
     if ( (ANY(ABS(rprimd(2:3,  1))>tol6)).or.&
&         (ANY(ABS(rprimd(1:3:2,2))>tol6)).or.&
&         (ANY(ABS(rprimd(1:2,  3))>tol6))    &
&       ) then
       msg = ' Bravais lattice should be orthorombic and parallel to the cartesian versors '
       MSG_ERROR(msg)
     end if

     ha_=half*SQRT(DOT_PRODUCT(rprimd(:,1),rprimd(:,1)))
     hb_=half*SQRT(DOT_PRODUCT(rprimd(:,2),rprimd(:,2)))
     r0_=MIN(ha_,hb_)/N0
     !
     ! ===================================================
     ! === Setup for the quadrature of matrix elements ===
     ! ===================================================
     qopt_    =6         ! Quadrature method, see quadrature routine.
     ntrial_  =30        ! Max number of attempts.
     accuracy_=0.001     ! Fractional accuracy required.
     npts_    =6         ! Initial number of point (only for Gauss-Legendre method).

     do ig=1,n1

       gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
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
     a1=rprimd(:,1); b1=two_pi*gprimd(:,1)
     a2=rprimd(:,2); b2=two_pi*gprimd(:,2)
     a3=rprimd(:,3); b3=two_pi*gprimd(:,3)

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
         if (check<zero) alpha(ii)=normv(check*rprimd(:,ii),rmet,'R')
       end if
     end do

     SELECT CASE (opt_surface)

     !CASE SURFACE 1 - Beigi
       CASE(1)

       do i3=1,n3
        do i2=1,n2
         i23=n1*(i2-1 + n2*(i3-1))
         do i1=1,n1
           ii=i1+i23
           gcart(:)=b1(:)*gvec(1,i1)+b2(:)*gvec(2,i2)+b3(:)*gvec(3,i3)
           gcart_para=SQRT(gcart(1)**2+gcart(2)**2) ; gcart_perp = gcart(3)
           if(ABS(gcart_para)>tol4.and.ABS(gcart_perp)>tol4) then
             gcutoff(ii)=one-EXP(-gcart_para*rcut)*COS(gcart_perp*rcut)
           else 
             gcutoff(ii)=zero
           end if
         end do !i1
        end do !i2
       end do !i3
        
       !CASE SURFACE 2 - Rozzi
       CASE(2)

       !!BG: Trigger needed - use the available input value for this 
       do i3=1,n3
        do i2=1,n2
         i23=n1*(i2-1 + n2*(i3-1))
         do i1=1,n1
           ii=i1+i23
           gcart(:)=b1(:)*gvec(1,i1)+b2(:)*gvec(2,i2)+b3(:)*gvec(3,i3)
           gcart_para=SQRT(gcart(1)**2+gcart(2)**2) ; gcart_perp = gcart(3)
           if(ABS(gcart_para)>tol4) then
             gcutoff(ii)=one+EXP(-gcart_para*rcut)*(gcart_perp/gcart_para*SIN(gcart_perp*rcut)-COS(gcart_perp*rcut))
           else if (ABS(gcart_perp)>tol4) then
             gcutoff(ii)=one-COS(gcart_perp*rcut)-gcart_perp/rcut*SIN(gcart_perp*rcut)
           else
             gcutoff(ii)=zero
           end if
         end do !i1
        end do !i2
       end do !i3
   
       CASE DEFAULT
         write(msg,'(a,i3)')' Wrong value of surface method: ',opt_surface
         MSG_BUG(msg)
       END SELECT

   CASE('ERF')

   ! Calculate rcut for each method ! Same as SPHERE
   rcut= (three*nkpt*ucvol/four_pi)**(one/three)
  
     do ig=1,nfft
       if(abs(gpq(ig))<tol4) then
          gcutoff(ig)=zero ! @Gamma: initialize quantity in each requiered routine
       else if(gpq(ig)<=cutoff) then
          gcutoff(ig)=exp(-pi/(gpq2(ig)*rcut**2))
       end if
     end do  !ig
 
   CASE('ERFC')

   ! Calculate rcut for each method ! Same as SPHERE
   rcut= (three*nkpt*ucvol/four_pi)**(one/three)

     do ig=1,nfft
       if(abs(gpq(ig))<tol4) then
          gcutoff(ig)=zero ! @Gamma: initialize quantity in each requiered routine
       else if(gpq(ig)<=cutoff) then
          gcutoff(ig)=one-exp(-pi/(gpq2(ig)*rcut**2))
       end if
     end do !ig

   CASE('CRYSTAL')
     gcutoff(:)=one ! Neutral cut-off
     write(msg,'(a)')'CRYSTAL method: no cut-off applied to G**2 while CRYSTAL method is implied!'
     MSG_WARNING(msg) 
   CASE DEFAULT
     gcutoff=one ! Neutral cut-off
     write(msg,'(a)')'No cut-off applied to G**2!'
     MSG_WARNING(msg)
 END SELECT

 ABI_DEALLOCATE(gvec) 
 ABI_DEALLOCATE(gpq)
 ABI_DEALLOCATE(gpq2)
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
