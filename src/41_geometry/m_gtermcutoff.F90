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
! use m_xmpi
 use m_fstrings,        only : sjoin, itoa
 use m_profiling_abi,   only : abimem_record
! use defs_abitypes,     only : MPI_type
 use m_bessel,          only : CALJY0, CALJY1, CALCK0, CALCK1
 use m_numeric_tools,   only : arth, l2norm, OPERATOR(.x.),quadrature
 use m_paw_numeric,     only : paw_jbessel
 use m_geometry,        only : normv, metric

 implicit none

 private
!!***

!!****t* m_gtermcutoff/gtermcut_t
!! NAME
!!  gtermcut_t
!!
!! FUNCTION
!!
!! SOURCE

!!! type,public :: gtermcut_t

!!!  integer  :: nfft
!!!   ! Number of points in FFT grid

!!!  integer  :: ng
!!!   ! Number of G-vectors

!!!  real(dp) :: ucvol
!!!  ! Volume of the unit cell

!!!   ! integer :: pdir(3)
!!!   ! 1 if the system is periodic along this direction

!!!   ! real(dp) :: boxcenter(3)
!!!   ! 1 if the point in inside the cutoff region 0 otherwise
!!!   ! Reduced coordinates of the center of the box (input variable)

!!!  real(dp) :: rprimd(3,3)
!!!    ! Lattice vectors in real space.

!!!  real(dp),allocatable :: gtermcuoff(:)
!!!    ! gtermcuoff(nfft)
!!!    ! G cut-off array on the FFT grid

!!! end type gtermcut_t

 public :: termcutoff
!!***
! private variables used for the integration needed by the cylindrical case.
 integer,save  :: npts_,ntrial_,qopt_
 real(dp),save :: ha_,hb_,hcyl_,r0_
 real(dp),save :: gcart_para_,gcart_perp_,gcartx_,gcarty_
 real(dp),save :: xx_,zz_,rcut_
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

subroutine termcutoff(gcutoff,gsqcut,icutcoul,ngfft,nkpt,rcut,rprimd,vcutgeo)
 
!Arguments ------------------------------------
!scalars
 integer,intent(in)   :: icutcoul, nkpt
 real(dp),intent(in)  :: gsqcut,rcut

!arrays
 integer,intent(in)    :: ngfft(18)
 real(dp),intent(in)   :: rprimd(3,3),vcutgeo(3)

!Local variables-------------------------------
!scalars
 integer,parameter  :: N0=1000
 integer            :: i1,i2,i23,i3,ierr,id(3),ii,ig,ing,icount
 integer            :: c1,c2,opt_cylinder
 integer            :: my_start,my_stop
 integer            :: n1,n2,n3,nfft
 integer            :: test,opt_surface !opt_cylinder
 real(dp)           :: alpha_fac, ap1sqrt, log_alpha
 real(dp)           :: cutoff,rcut_loc,rcut2,check,rmet(3,3)
 real(dp)           :: gvecg2p3,gvecgm12,gvecgm13,gvecgm23,gs2,gs3
 real(dp)           :: gcart_para,gcart_perp,gcart_x,gcart_y,gcart_xy,gcart_z
 real(dp)           :: j0,j1,k0,k1
 real(dp)           :: quad,tmp,ucvol
 real(dp)           :: hcyl,hcyl2
 real(dp),parameter :: tolfix=1.0000001_dp
 character(len=50)  :: mode
 character(len=500) :: msg
! type(gcut_t)       :: gcut  !

!arrays
 real(dp)             :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
 real(dp)             :: gcart(3),gmet(3,3),gprimd(3,3)
 real(dp)             :: pdir(3),alpha(3)
 real(dp),allocatable :: gvec(:,:),gpq(:),gpq2(:),xx(:)
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
     if(rcut>tol4) then
         rcut_loc = rcut
     else
         rcut_loc = (three*nkpt*ucvol/four_pi)**(one/three)
     endif

     do ig=1,nfft
       if(abs(gpq(ig))<tol4) then
          gcutoff(ig)=zero
       else
          gcutoff(ig)=one-cos(rcut_loc*sqrt(four_pi/gpq2(ig)))
      end if
     end do

   CASE('CYLINDER')

     test=COUNT(vcutgeo/=zero)
     ABI_CHECK(test==1,'Wrong cutgeo for cylinder')   

     ! === From reduced to Cartesian coordinates ===
     call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

     a1=rprimd(:,1); b1=two_pi*gprimd(:,1)
     a2=rprimd(:,2); b2=two_pi*gprimd(:,2)
     a3=rprimd(:,3); b3=two_pi*gprimd(:,3)
        
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

     ! === Beigi method is the default one, i.e infinite cylinder of radius rcut ===
     ! * Negative values to use Rozzi method with finite cylinder of extent hcyl.
     opt_cylinder=1; hcyl=zero; pdir(:)=0
     do ii=1,3
       check=vcutgeo(ii)
       if (ABS(check)>tol6) then
         pdir(ii)=1
         if (check<zero) then  ! use Rozzi's method.
           hcyl=ABS(check)*SQRT(SUM(rprimd(:,ii)**2))
           opt_cylinder=2
         end if
       end if
     end do

     ! Calculate rcut for each method
     if(rcut>tol4) then
       rcut_loc = rcut
     else
       rcut_loc = half*SQRT(DOT_PRODUCT(a1,a1))
     endif

     if (opt_cylinder==1) then
       ABI_CHECK(ALL(pdir == (/0,0,1/)),"Surface must be in the x-y plane")
     end if

     rcut_= rcut_loc

     ha_=half*SQRT(DOT_PRODUCT(rprimd(:,1),rprimd(:,1)))
     hb_=half*SQRT(DOT_PRODUCT(rprimd(:,2),rprimd(:,2)))
     r0_=MIN(ha_,hb_)/N0
     !
     ! ===================================================
     ! === Setup for the quadrature of matrix elements ===
     ! ===================================================
     qopt_    =6        ! Quadrature method, see quadrature routine.
     ntrial_  =30       ! Max number of attempts.
     accuracy_=0.001    ! Fractional accuracy required.
     npts_    =10       ! Initial number of point (only for Gauss-Legendre method).
     hcyl_    =hcyl     ! Lenght of cylinder along z, only if method==2

     write(msg,'(3a,2(a,i5,a),a,f8.5)')ch10,&
&      ' cutoff_cylinder: Info on the quadrature method : ',ch10,&
&      '  Quadrature scheme      = ',qopt_,ch10,&
&      '  Max number of attempts = ',ntrial_,ch10,&
&      '  Fractional accuracy    = ',accuracy_
     call wrtout(std_out,msg,'COLL')

     SELECT CASE (opt_cylinder)
     
     CASE(1)
   
     do i3=1,n3
      do i2=1,n2
       i23=n1*(i2-1 + n2*(i3-1))
       do i1=1,n1
         ii=i1+i23

         gcart(:)=b1(:)*gvec(1,i1)+b2(:)*gvec(2,i2)+b3(:)*gvec(3,i3)
         gcart_x=gcart(1) ; gcart_y=gcart(2) ; gcart_z=ABS(gcart(3))
         gcart_xy = SQRT(gcart_x**2+gcart_y**2) ;

         if (gcart_z>tol4) then
           ! === Analytic expression ===
           call CALCK0(gcart_z *rcut_loc,k0,1)
           call CALJY1(gcart_xy*rcut_loc,j1,0)
           call CALJY0(gcart_xy*rcut_loc,j0,0)
           call CALCK1(gcart_z *rcut_loc,k1,1)
           gcutoff(ii)=one+rcut_loc*gcart_xy*j1*k0-gcart_z*rcut_loc*j0*k1
         else
           if (gcart_xy>tol4) then
             ! === Integrate r*Jo(G_xy r)log(r) from 0 up to rcut_  ===
             call quadrature(F5,zero,rcut_loc,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
             if (ierr/=0) then
               MSG_ERROR("Accuracy not reached")
             end if
               gcutoff(ii)= -quad*gpq(ii)
           else
               gcutoff(ii)= rcut_loc**2*(two*LOG(rcut_)-one)*gpq(ii)/four
          end if
         end if

       end do !i1
      end do !i2
     end do !i3

     CASE(2)

     ! === Finite cylinder of length hcyl, from Rozzi et al ===
     ! TODO add check on hcyl value that should be smaller that 1/deltaq
     if (hcyl_<zero) then
       write(msg,'(a,f8.4)')' Negative value for cylinder length hcyl=',hcyl_
       MSG_BUG(msg)
     end if
      
     npts_=6
     ABI_MALLOC(xx,(npts_))

     write(msg,'(2(a,f8.4))')' cutoff_cylinder: using finite cylinder of length= ',hcyl,' rcut= ',rcut_loc
     call wrtout(std_out,msg,'COLL')
     hcyl_=hcyl
     hcyl2=hcyl**2
     rcut2=rcut_loc**2

     write(*,*)'This',n1,n2,n3

     do i3=1,n3
      do i2=1,n2
       i23=n1*(i2-1 + n2*(i3-1))
       do i1=1,n1
       ii=i1+i23

         gcart(:)=b1(:)*gvec(1,i1)+b2(:)*gvec(2,i2)+b3(:)*gvec(3,i3)
         gcart_para_=ABS(gcart(3)) ; gcart_perp_=SQRT(gcart(1)**2+gcart(2)**2)

         if (gcart_perp_/=zero.and.gcart_para_/=zero) then
           call quadrature(F2,zero,rcut_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
           if (ierr/=0) then
             MSG_ERROR("Accuracy not reached")
           end if

           gcutoff(ii)=SQRT(1/quad)

         else if (gcart_perp_==zero.and.gcart_para_/=zero) then

           ! $ \int_0^h sin(qpg_para_.z)/\sqrt(rcut^2+z^2)dz $
           call quadrature(F3,zero,hcyl,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
           if (ierr/=0) then
             MSG_ERROR("Accuracy not reached")
           end if

           c1=one/gcart_para_**2-COS(gcart_para_*hcyl_)/gcart_para_**2-hcyl_*SIN(gcart_para_*hcyl_)/gcart_para_
           c2=SIN(gcart_para_*hcyl_)*SQRT(hcyl2+rcut2)
           gcutoff(ii)=SQRT(1/(c1+(c2-quad)/gcart_para_))

         else if (gcart_perp_/=zero.and.gcart_para_==zero) then
           ! $ 4pi\int_0^rcut d\rho \rho J_o(qpg_perp_.\rho) ln((h+\sqrt(h^2+\rho^2))/\rho) $
           call quadrature(F4,zero,rcut_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
           if (ierr/=0) then
             MSG_ERROR("Accuracy not reached")
           end if

           gcutoff(ii)=SQRT(1/quad)

         else if (gcart_perp_==zero.and.gcart_para_==zero) then
           ! Use lim q+G --> 0
           gcutoff(ii)=zero
         else
           MSG_BUG('You should not be here!')
         end if

       end do !i1
      end do !i2
     end do !i3

     ABI_FREE(xx)

     CASE DEFAULT
      MSG_BUG(sjoin('Wrong value for cylinder method:',itoa(opt_cylinder)))
     END SELECT


   CASE('SURFACE')

     test=COUNT(vcutgeo/=zero)
     ABI_CHECK(test==2,"Wrong vcutgeo")

     ! === From reduced to cartesian coordinates ===
     a1=rprimd(:,1); b1=two_pi*gprimd(:,1)
     a2=rprimd(:,2); b2=two_pi*gprimd(:,2)
     a3=rprimd(:,3); b3=two_pi*gprimd(:,3)

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

       ! Calculate rcut for each method !
       if(rcut>tol4) then
          rcut_loc = rcut
       else
          rcut_loc = half*SQRT(DOT_PRODUCT(a3,a3))
       endif
 
       do i3=1,n3
        do i2=1,n2
         i23=n1*(i2-1 + n2*(i3-1))
         do i1=1,n1
           ii=i1+i23
           gcart(:)=b1(:)*gvec(1,i1)+b2(:)*gvec(2,i2)+b3(:)*gvec(3,i3)
           gcart_para=SQRT(gcart(1)**2+gcart(2)**2) ; gcart_perp = gcart(3)
           if(gcart_para<tol4.and.ABS(gcart_perp)<tol4) then
             gcutoff(ii)=zero
           else
             gcutoff(ii)=one-EXP(-gcart_para*rcut_loc)*COS(gcart_perp*rcut_loc)
           end if
         end do !i1
        end do !i2
       end do !i3
        
       !CASE SURFACE 2 - Rozzi
       CASE(2)

       !Set the cut-off radius
       if(rcut>tol4) then
          rcut_loc = rcut
       else
          rcut_loc = half*SQRT(DOT_PRODUCT(a3,a3))
       endif

       !In the case of finite, Rozzi's method provide another parameter
       !for the cut-off: alpha
       !!! ATT: alpha = L_x/L_y --> in-plane geometry dependence
       alpha_fac=SQRT(DOT_PRODUCT(a1,a1))/SQRT(DOT_PRODUCT(a2,a2))
       ap1sqrt=SQRT(one+alpha_fac**2)
       log_alpha=LOG((alpha_fac+ap1sqrt)*(one+ap1sqrt)/alpha_fac)

       do i3=1,n3
        do i2=1,n2
         i23=n1*(i2-1 + n2*(i3-1))
         do i1=1,n1
           ii=i1+i23
           gcart(:)=b1(:)*gvec(1,i1)+b2(:)*gvec(2,i2)+b3(:)*gvec(3,i3)
           gcart_para=SQRT(gcart(1)**2+gcart(2)**2) ; gcart_perp = gcart(3)
           if(gcart_para>tol4) then
             gcutoff(ii)=one+EXP(-gcart_para*rcut_loc)*(gcart_perp/gcart_para*&
&                        SIN(gcart_perp*rcut_loc)-COS(gcart_perp*rcut_loc)) 
           else
             if (ABS(gcart_perp)>tol4) then
               gcutoff(ii)=one-COS(-gcart_perp*rcut_loc)-gcart_perp*rcut_loc*SIN(gcart_perp*rcut_loc)
!               gcutoff(ii)=one-COS(-gcart_perp*rcut_loc)-SIN(gcart_perp*rcut_loc) - Altered Rozzi's
             else
               gcutoff(ii)=zero
             endif
           endif
         end do !i1
        end do !i2
       end do !i3
   
       CASE DEFAULT
         write(msg,'(a,i3)')' Wrong value of surface method: ',opt_surface
         MSG_BUG(msg)
       END SELECT

   CASE('ERF')

    ! Calculate rcut for each method ! Same as SPHERE
    if(rcut>tol4) then
        rcut_loc = rcut
    else
        rcut_loc= (three*nkpt*ucvol/four_pi)**(one/three)
    endif
  
     do ig=1,nfft
       if(abs(gpq(ig))<tol4) then
          gcutoff(ig)=zero ! @Gamma: initialize quantity in each requiered routine
       else !if(gpq(ig)<=cutoff) then
          gcutoff(ig)=exp(-pi/(gpq2(ig)*rcut_loc**2))
       end if
     end do  !ig
 
   CASE('ERFC')

   ! Calculate rcut for each method ! Same as SPHERE
     if(rcut>tol4) then
         rcut_loc = rcut
     else
         rcut_loc= (three*nkpt*ucvol/four_pi)**(one/three)
     endif

     do ig=1,nfft
       if(abs(gpq(ig))<tol4) then
          gcutoff(ig)=zero ! @Gamma: initialize quantity in each requiered routine
       else
          gcutoff(ig)=one-exp(-pi/(gpq2(ig)*rcut_loc**2))
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
 arg=rho*gcart_perp_
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

 F2=intr*COS(gcart_para_*xx)

end function F2
!!***

!----------------------------------------------------------------------

pure function F3(xx)

 real(dp),intent(in) :: xx
 real(dp) :: F3
!************************************************************************

 ! F3(z)=z*\sin(qpg_para_*z)/\sqrt(rcut^2+z^2)
 F3=xx*SIN(gcart_para_*xx)/SQRT(rcut_**2+xx**2)

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
   arg=rho*gcart_perp_
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
   arg=rho*gcart_perp_
   call paw_jbessel(bes,besp,bespp,ll,order,arg)
   F5=bes*rho*LOG(rho)
 end if

end function F5
!!***

end module m_gtermcutoff
!!***
