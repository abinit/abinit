!!****m* ABINIT/m_gammapositron
!! NAME
!!  m_gammapositron
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (MT,GJ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_gammapositron

 use defs_basis
 use m_abicore
 use m_errors
 use m_electronpositron

 use defs_abitypes,     only : MPI_type
 use m_numeric_tools,   only : invcb
 use m_xctk,            only : xcden

 implicit none

 private
!!***

 public :: gammapositron
 public :: gammapositron_fft
!!***

contains
!!***

!!****f* ABINIT/gammapositron
!! NAME
!! gammapositron
!!
!! FUNCTION
!! Compute positron electron-positron enhancement factor (contact density) used to compute positron lifetime.
!! Input is positronic rhop(r) and electronic rhoe(r) at a given set of points.
!!
!! INPUTS
!!  grhocore2(ngr)=square of the gradient of core electronic density rhocore (needed for GGA)
!!  grhoe2(ngr)=square of the gradient of valence electronic density rhoer (needed for GGA)
!!  igamma=type of enhancement factor:
!!     1:  Boronski and Nieminen [2]
!!     2:  Boronski and Nieminen, RPA limit [2]
!!     3:  Sterne and Kaiser [3]
!!     4:  Puska, Seitsonen and Nieminen [4]
!!     See references below
!!  ngr=size of grho2 array (0 if LDA, npt if GGA)
!!  npt=number of real space points on which density is provided
!!  rhocore(npt*usecore)=core electron density (bohr^-3)
!!  rhoer(npt)  =valence electron density (bohr^-3)
!!  rhopr(npt)  =positron density (bohr^-3)
!!  usecore     =1 if core density is not zero
!!
!! OUTPUT
!!  gamma(npt,2)    =electron-positron enhancement factor,
!!                    gamma(:,1): using total   electronic density
!!                    gamma(:,2): using valence electronic density
!!
!! NOTES
!!   References for electron-positron correlation functionals:
!!         [1] J. Arponen and E. Pajanne, Ann. Phys. (N.Y.) 121, 343 (1979) [[cite:Arponen1979a]].
!!         [2] E. Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986) [[cite:Boronski1986]].
!!         [3] P.A. Sterne and J.H. Kaiser, Phys. Rev. B 43, 13892 (1991) [[cite:Sterne1991]].
!!         [4] M.J. Puska, A.P. Seitsonen and R.M. Nieminen, Phys. Rev. B 52, 10947 (1994) [[cite:Puska1994]].
!!         [5] B. Barbiellini, M.J. Puska, T. Torsti and R.M.Nieminen, Phys. Rev. B 51, 7341 (1995) [[cite:Barbiellini1995]]
!!
!! PARENTS
!!      gammapositron_fft,poslifetime,posratecore
!!
!! CHILDREN
!!      invcb
!!
!! SOURCE

subroutine gammapositron(gamma,grhocore2,grhoe2,igamma,ngr,npt,rhocore,rhoer,rhopr,usecore)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: igamma,ngr,npt,usecore
!arrays
 real(dp),intent(in) :: grhocore2(ngr*usecore),grhoe2(ngr),rhocore(npt*usecore),rhoer(npt),rhopr(npt)
 real(dp),intent(out) :: gamma(npt,2)

!Local variables-------------------------------
!scalars
 integer :: iloop,ipt
 logical :: gga
 real(dp),parameter :: alpha_gga=0.22d0,rsfac=0.6203504908994000_dp
 real(dp) :: aa,bb,cc,dg1
 real(dp) :: drs,eps,expgga,g0,g1,g2,gg
 real(dp) :: kf,kk,nqtf2,ratio1,ratio2,ratio3,rho1,rho2,rhoe,rhop,sqrs,rs,rse,rsp
!arrays
 real(dp),allocatable :: grho2(:),rhor(:),rsepts(:),rsppts(:)

! *************************************************************************

 gga=(ngr==npt.and.igamma/=0)

 if (usecore/=0.and.usecore/=1) then
   MSG_ERROR('Wrong value for usecore !')
 end if
 if (igamma/=0.and.igamma/=1.and.igamma/=2.and.igamma/=3.and.igamma/=4) then
   MSG_ERROR('Unknown electron-positron correlation !')
 end if

 ABI_ALLOCATE(rhor,(npt))
 ABI_ALLOCATE(rsepts,(npt))
 if (gga)  then
   ABI_ALLOCATE(grho2,(npt))
 end if

!Eventually compute positronic density radii
 if (igamma==1.or.igamma==4) then
   ABI_ALLOCATE(rsppts,(npt))
   call invcb(rhopr(:),rsppts,npt)
   rsppts(:)=rsfac*rsppts(:)
 end if

!Loop: iloop=1: compute enhancement factor using total   electronic density
!iloop=2: compute enhancement factor using valence electronic density
!===================================================================================
 do iloop=1,1+usecore

!  Compute electronic density radii
   if (iloop==1.and.usecore==1) then
     rhor(1:npt)=rhoer(1:npt)+rhocore(1:npt)
   else
     rhor(1:npt)=rhoer(1:npt)
   end if
   call invcb(rhor(:),rsepts,npt)
   rsepts(:)=rsfac*rsepts(:)

!  Gradients for GGA
   if (gga) then
     if (iloop==1.and.usecore==1) then
       grho2(1:npt)=grhoe2(1:npt)+grhocore2(1:npt)
     else
       grho2(1:npt)=grhoe2(1:npt)
     end if
   end if

!  Loop over grid points
   do ipt=1,npt

     rhoe=rhor(ipt)
     rhop=rhopr(ipt)
     rse =rsepts(ipt)
     gg=zero

!    Testing feature: gamma=1
!    -----------------------------------------------------------------------------------
     if (igamma==0) then

       gg=one

!    Boronski and Nieminen
!    -----------------------------------------------------------------------------------
     else if (igamma==1) then

       rsp =rsppts(ipt)
       if (rhoe>rhop) then
         rho1=rhoe;rho2=rhop;rs=rse
       else
         rho1=rhop;rho2=rhoe;rs=rsp
       end if
       drs=-third*rs/rho1;sqrs=sqrt(rs)
       ratio1=rho2/rho1;ratio2=ratio1*ratio1;ratio3=ratio2*ratio1
       g0=one+1.23_dp*rs+0.8295_dp*sqrs**3-1.26_dp*rs**2+0.3286_dp*sqrs**5+sixth   *rs**3
       g1=one+0.51_dp*rs                  +0.65_dp*rs**2-0.51_dp  *sqrs**5+0.176_dp*rs**3
       g2=one+0.60_dp*rs                  +0.63_dp*rs**2-0.48_dp  *sqrs**5+0.167_dp*rs**3
       dg1=drs*(0.51_dp+two*0.65_dp*rs-2.5_dp*0.51_dp*sqrs**3+three*0.176_dp*rs**2)
       kk=half*rho1*dg1
       aa= two  *kk-six    *g1+eight  *g2-two *g0
       bb=-three*kk+11.0_dp*g1-16.0_dp*g2+five*g0
       cc=       kk-four   *g1+eight  *g2-four*g0
       gg=g0+ratio3*aa+ratio2*bb+ratio1*cc

!      Boronski and Nieminen RPA limit
!      -----------------------------------------------------------------------------------
     else if (igamma==2) then

       rs=rse;sqrs=sqrt(rs)
       gg=one !This is experimental to avoid divergences
       if (rs<=20._dp) gg=gg+1.23_dp*rs+0.8295_dp*sqrs**3-1.26_dp*rs**2+0.3286_dp*sqrs**5+sixth*rs**3

!      Sterne and Kaiser
!      -----------------------------------------------------------------------------------
     else if (igamma==3) then

       rs=rse;sqrs=sqrt(rs)
       gg=one !This is experimental to avoid divergences
       if (rs<=20._dp) gg=gg+0.1512_dp*rs+2.414_dp*sqrs**3-2.01_dp*rs**2+0.4466_dp*sqrs**5+0.1667_dp*rs**3

!      Puska, Seitsonen and Nieminen
!      -----------------------------------------------------------------------------------
     else if (igamma==4) then

       rsp =rsppts(ipt)
       if (rhoe>rhop) then
         rho1=rhoe;rho2=rhop;rs=rse
       else
         rho1=rhop;rho2=rhoe;rs=rsp
       end if
       drs=-third*rs/rho1;sqrs=sqrt(rs)
       ratio1=rho2/rho1;ratio2=ratio1*ratio1;ratio3=ratio2*ratio1
       g0=one+1.2300_dp*rs+0.9889_dp*sqrs**3-1.4820_dp*rs**2+0.3956_dp*sqrs**5+sixth*rs**3
       g1=one+2.0286_dp*rs-3.3892_dp*sqrs**3+3.0547_dp*rs**2-1.0540_dp*sqrs**5+sixth*rs**3
       g2=one+0.2499_dp*rs+0.2949_dp*sqrs**3+0.6944_dp*rs**2-0.5339_dp*sqrs**5+sixth*rs**3
       dg1=drs*(2.0286_dp-1.5_dp*3.3892_dp*sqrs+two*3.0547_dp*rs-2.5_dp*1.0540_dp*sqrs**3+three*sixth*rs**2)
       kk=half*rho1*dg1
       aa= two  *kk-six    *g1+eight  *g2-two *g0
       bb=-three*kk+11.0_dp*g1-16.0_dp*g2+five*g0
       cc=       kk-four   *g1+eight  *g2-four*g0
       gg=g0+ratio3*aa+ratio2*bb+ratio1*cc

     end if ! igamma

     if (gga) then
       kf=(three*pi*pi*rhoe)**third
       nqtf2=(rhoe*sqrt(four*kf/pi))**2
       eps=grho2(ipt)/nqtf2
       if (eps<zero) then
         MSG_ERROR('  problem, negative GGA espilon !')
       end if
       expgga=exp(-alpha_gga*eps*third)
       gg=one+(gg-one)*expgga
     end if

!    Store enhancement factor
     gamma(ipt,iloop)=gg

   end do ! ipt
 end do ! iloop

 ABI_DEALLOCATE(rhor)
 ABI_DEALLOCATE(rsepts)
 if (igamma==1.or.igamma==4)  then
   ABI_DEALLOCATE(rsppts)
 end if
 if (gga)  then
   ABI_DEALLOCATE(grho2)
 end if

!Case usecore=0 (no core density)
 if (usecore==0) gamma(:,2)=gamma(:,1)

end subroutine gammapositron
!!***

!!****f* ABINIT/gammapositron_fft
!! NAME
!! gammapositron_fft
!!
!! FUNCTION
!! Compute positron electron-positron enhancement factor on a real space FFT grid.
!!
!! INPUTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  gprimd(3,3)=dimensional primitive translations in reciprocal space (bohr^-1)
!!  igamma=type of enhancement factor:
!!    -1:  gamma=one (test)
!!     1:  Boronski and Nieminen [1]
!!     2:  Boronski and Nieminen, RPA limit [1]
!!     3:  Sterne and Kaiser [2]
!!     4:  Puska, Seitsonen and Nieminen [3]
!!  mpi_enreg=information about MPI parallelization
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  nfft=number of FFT grid points
!!  ngfft(18)=contain all needed information about 3D FFT
!!  rhor_e(nfft)=real space density electronic density (total density)
!!  rhor_p(nfft)=real space density positronic density
!!  xccc3d(n3xccc)=3D core electron density for XC core correction
!!
!! OUTPUT
!!  gamma(nfft,2)=electron-positron enhancement factor,
!!                gamma(:,1): using total   electronic density
!!                gamma(:,2): using valence electronic density
!!
!! NOTES
!!  The input densities (rhor_e and rhor_p) should be positive everywhere
!!  (call mkdenpos routine before entering this one)
!!
!! PARENTS
!!      posdoppler,poslifetime
!!
!! CHILDREN
!!      gammapositron,xcden
!!
!! SOURCE

subroutine gammapositron_fft(electronpositron,gamma,gprimd,igamma,mpi_enreg,&
&                            n3xccc,nfft,ngfft,rhor_e,rhor_p,xccc3d)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: igamma,n3xccc,nfft
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),rhor_e(nfft),rhor_p(nfft),xccc3d(n3xccc)
 real(dp),intent(out) :: gamma(nfft,2)

!Local variables-------------------------------
!scalars
 integer :: cplex,ishift,ngr,ngrad,nspden_ep,usecore
!arrays
 real(dp),parameter :: qphon(3)=(/zero,zero,zero/)
 real(dp),allocatable :: grhocore2(:),grhoe2(:),rhoc(:,:,:),rhoe(:,:,:)

! *************************************************************************

!Several useful constants
 usecore=n3xccc/nfft
 cplex=1;ishift=0;ngrad=1;nspden_ep=1
 if (electronpositron%ixcpositron==3.or.electronpositron%ixcpositron==31) ngrad=2
 ngr=0;if (ngrad==2) ngr=nfft

!Allocate several arrays
 ABI_ALLOCATE(rhoe,(nfft,nspden_ep,ngrad**2))
 ABI_ALLOCATE(grhoe2,(ngr))
 ABI_ALLOCATE(grhocore2,(ngr*usecore))

!Store electronic density and its gradients
 call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_ep,qphon,rhor_e,rhoe)

!Compute squared gradient of the electronic density
 if (ngrad==2) then
   grhoe2(:)=rhoe(:,1,2)**2+rhoe(:,1,3)**2+rhoe(:,1,4)**2
   if (usecore>0) then
     ABI_ALLOCATE(rhoc,(nfft,1,ngrad**2))
     call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden_ep,qphon,xccc3d,rhoc)
     grhocore2(:)=rhoc(:,1,2)**2+rhoc(:,1,3)**2+rhoc(:,1,4)**2
     ABI_DEALLOCATE(rhoc)
   end if
 end if

!Compute enhancement factor on FFT grid
 call gammapositron(gamma,grhocore2,grhoe2,igamma,ngr,nfft,xccc3d,&
& rhoe(:,1,1),rhor_p,usecore)

!Release temporary memory
 ABI_DEALLOCATE(rhoe)
 ABI_DEALLOCATE(grhoe2)
 ABI_DEALLOCATE(grhocore2)

end subroutine gammapositron_fft
!!***

end module m_gammapositron
!!***
