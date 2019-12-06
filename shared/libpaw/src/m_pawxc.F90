!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pawxc
!! NAME
!!  m_pawxc
!!
!! FUNCTION
!!  XC+PAW related operations
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (MT, FJ, TR, GJ, TD)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

module m_pawxc

 USE_DEFS
 USE_MSG_HANDLING
 USE_MEMORY_PROFILING

#ifdef LIBPAW_ISO_C_BINDING
 use iso_c_binding, only : c_ptr,c_loc,c_f_pointer
#endif

#ifdef HAVE_LIBPAW_ABINIT
 use m_xcpositron,  only : xcpositron
 use m_drivexc,     only : drivexc_main, size_dvxc, xcmult, mkdenpos
 use m_xc_noncoll,  only : rotate_mag,rotate_back_mag,rotate_back_mag_dfpt
#endif

 use m_libpaw_libxc

 use m_pawang,      only : pawang_type
 use m_pawrad,      only : pawrad_type, nderiv_gen, pawrad_deducer0, simp_gen

 implicit none

 private

 public :: pawxc          ! Compute xc correlation potential and energies inside a paw sphere. USE (r,theta,phi)
 public :: pawxcpositron  ! Compute electron-positron correlation potential and energies inside a PAW sphere. USE (r,theta,phi)
 public :: pawxc_dfpt     ! Compute first-order change of XC potential and contribution to
                          !   2nd-order change of XC energy inside a PAW sphere. USE (r,theta,phi)
 public :: pawxcsum       ! Compute useful sums of moments of densities needed to compute on-site contributions to XC energy and potential
 public :: pawxcm         ! Compute xc correlation potential and energies inside a paw sphere. USE (L,M) MOMENTS
 public :: pawxcmpositron ! Compute electron-positron correlation potential and energies inside a PAW sphere. USE (L,M) MOMENTS
 public :: pawxcm_dfpt    ! Compute 1st-order change of XC potential and contrib
                          !   to 2nd-order change of XC ene inside a PAW sphere. USE (L,M) MOMENTS
 public :: pawxc_get_nkxc    ! Compute size of XC kernel (Kxc) according to spin polarization and XC type
 public :: pawxc_get_usekden ! Assess whether kinetic energy density has to be computed

!Private procedures
 private :: pawxcsph                   ! Compute XC energy and potential for a spherical density rho(r) given as (up,dn)
 private :: pawxcsphpositron           ! Compute electron-positron XC energy and potential for spherical densities rho_el(r) rho_pos(r)
 private :: pawxcsph_dfpt              ! Compute XC 1st-order potential for a 1st-order spherical density rho1(r)
 private :: pawxc_rotate_mag           ! Rotate a non-collinear density wrt a magnetization
 private :: pawxc_rotate_back_mag      ! Rotate back a collinear XC potential wrt a magnetization
 private :: pawxc_rotate_back_mag_dfpt ! Rotate back a collinear 1st-order XC potential wrt a magnetization

!Wrappers
 private :: pawxc_drivexc_wrapper    ! wrapper for drivexc_main
 private :: pawxc_mkdenpos_wrapper   ! wrapper for mkdenpos
 private :: pawxc_xcmult_wrapper     ! wrapper for xcmult
 private :: pawxc_size_dvxc_wrapper  ! wrapper for size_dvxc
 private :: pawxc_xcpositron_wrapper ! wrapper for xcpositron

!Zero of density
 real(dp),parameter :: rho_min=tol14
!!***

CONTAINS !===========================================================
!!***

!!****f* m_pawxc/pawxc_xcpositron_wrapper
!! NAME
!! pawxc_xcpositron_wrapper
!!
!! FUNCTION
!! Compute electron-positron correlation potentials and energy density.
!! Used electron-positron correlation functional is controlled by ipawxc_xcpositron_wrapper argument.
!! Returns Fxc, Vxc_pos, Vxc_el from input rhor_pos and rhor_el for positron and electrons.
!!
!! INPUTS
!!  grhoe2(ngr)=square of the gradient of electronic density rhoe (needed for GGA)
!!  ixcpositron=type of electron-positron correlation functional:
!!     1 or -1:  LDA zero positron density limit parametrized by Arponen & Pajanne
!!         and provided by Boronski & Nieminen [1,2]
!!     11: LDA zero positron density limit parametrized by Arponen & Pajanne
!!         and fitted by Sterne & Kaiser [1,3]
!!     2:  LDA electron-positron correlation
!!         provided by Puska, Seitsonen, and Nieminen [1,4]
!!     3:  GGA zero positron density limit parametrized by Arponen & Pajanne
!!         and provided by Boronski & Nieminen [1,2,5]
!!     31: GGA zero positron density limit parametrized by Arponen & Pajanne
!!         and fitted by Sterne & Kaiser [1,3,5]
!!     See references below
!!  ngr=size of grho2 array (0 if LDA, npt if GGA)
!!  npt=number of real space points on which density is provided
!!  posdensity0_limit=True if we are in the zero positron density limit
!!  rhoer(npt)=electron density (bohr^-3)
!!  rhopr(npt)=positron density (bohr^-3)
!!
!! OUTPUT
!!  fnxc(npt)=correlation energy per unit volume fxc
!!  vxce(npt)=correlation potential for electron dfxc/drhoe (hartree)
!!  vxcp(npt)=correlation potential for positron dfxc/drhop (hartree)
!!  vxcegr(ngr)= 1/|gradRhoe| dfxc/d|gradRhoe| (empty if LDA, i.e. ngr=0)
!!  Optional outputs:
!!    dvxce(npt)=partial second derivatives of the xc energy wr to the electronic density
!!               dvxce(:)=dVxce/dRhoe
!!    dvxcp(npt)=partial second derivatives of the xc energy wr to the positronic density
!!               dvxcp(:)=dVxcp/drhop
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
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxc_xcpositron_wrapper(fnxc,grhoe2,ixcpositron,ngr,npt,posdensity0_limit,&
&                                   rhoer,rhopr,vxce,vxcegr,vxcp,&
&                                   dvxce,dvxcp) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixcpositron,ngr,npt
 logical,intent(in) :: posdensity0_limit
!arrays
 real(dp),intent(in) :: grhoe2(ngr),rhoer(npt),rhopr(npt)
 real(dp),intent(out) :: fnxc(npt),vxce(npt),vxcegr(ngr),vxcp(npt)
 real(dp),intent(out),optional :: dvxce(npt),dvxcp(npt)

!Local variables-------------------------------

! *************************************************************************

#if defined HAVE_LIBPAW_ABINIT
 call pawxc_xcpositron_abinit()
#else
 call pawxc_xcpositron_local()
#endif
!!***

contains
!!***

#if defined HAVE_LIBPAW_ABINIT
!!****f* pawxc_xcpositron_wrapper/pawxc_xcpositron_abinit
!! NAME
!!  pawxc_xcpositron_abinit
!!
!! FUNCTION
!!  ABINIT version of electron-positron correlation
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxc_xcpositron_abinit()

! *************************************************************************

 if(present(dvxce) .and. present(dvxcp)) then
  call xcpositron(fnxc,grhoe2,ixcpositron,ngr,npt,posdensity0_limit,rhoer,rhopr,vxce,vxcegr,vxcp,&
&  dvxce=dvxce,dvxcp=dvxcp) ! optional arguments
 elseif( present(dvxce) .and. .not. present(dvxcp)) then
  call xcpositron(fnxc,grhoe2,ixcpositron,ngr,npt,posdensity0_limit,rhoer,rhopr,vxce,vxcegr,vxcp,&
&  dvxce=dvxce) ! optional arguments
 elseif( .not. present(dvxce) .and. present(dvxcp)) then
  call xcpositron(fnxc,grhoe2,ixcpositron,ngr,npt,posdensity0_limit,rhoer,rhopr,vxce,vxcegr,vxcp,&
&  dvxcp=dvxcp) ! optional arguments
 else
  call xcpositron(fnxc,grhoe2,ixcpositron,ngr,npt,posdensity0_limit,rhoer,rhopr,vxce,vxcegr,vxcp)
 end if

end subroutine pawxc_xcpositron_abinit
!!***

#else
!!****f* pawxc_xcpositron_wrapper/pawxc_xcpositron_local
!! NAME
!!  pawxc_xcpositron_local
!!
!! FUNCTION
!!  Local version of electron-positron correlation (to use outside ABINIT)
!!  NOT AVAILABLE
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxc_xcpositron_local()

 character(len=*), parameter :: msg='xcpositron only available in ABINIT!'

! *************************************************************************

 MSG_BUG(msg)

end subroutine pawxc_xcpositron_local
!!***
#endif

end subroutine pawxc_xcpositron_wrapper
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxc_size_dvxc_wrapper
!! NAME
!! pawxc_size_dvxc_wrapper
!!
!! FUNCTION
!! Give the size of the array dvxc(npts,ndvxc) and the second dimension of the d2vxc(npts,nd2vxc)
!! needed for the allocations depending on the routine which is called from the drivexc routine
!!
!! INPUTS
!!  ixc= choice of exchange-correlation scheme
!!  order=gives the maximal derivative of Exc computed.
!!    1=usual value (return exc and vxc)
!!    2=also computes the kernel (return exc,vxc,kxc)
!!   -2=like 2, except (to be described)
!!    3=also computes the derivative of the kernel (return exc,vxc,kxc,k3xc)
!!
!! OUTPUT
!!  ndvxc size of the array dvxc(npts,ndvxc) for allocation
!!  ngr2 size of the array grho2_updn(npts,ngr2) for allocation
!!  nd2vxc size of the array d2vxc(npts,nd2vxc) for allocation
!!  nvxcdgr size of the array dvxcdgr(npts,nvxcdgr) for allocation
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxc_size_dvxc_wrapper(ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order)

!Arguments----------------------
 integer, intent(in) :: ixc,nspden,order
 integer, intent(out) :: ndvxc,nd2vxc,ngr2,nvxcdgr

! *************************************************************************

#if defined HAVE_LIBPAW_ABINIT
 call size_dvxc(ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order)
#else
 call pawxc_size_dvxc_local()
#endif
!!***

#if ! defined HAVE_LIBPAW_ABINIT
contains
!!***

!!****f* pawxc_size_dvxc_wrapper/pawxc_size_dvxc_local
!! NAME
!!  pawxc_size_dvxc_local
!!
!! FUNCTION
!!  Local version of size_dvxc routine (to use outside ABINIT)
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxc_size_dvxc_local()

! *************************************************************************

 ngr2=0;nvxcdgr=0;ndvxc=0;nd2vxc=0

!Dimension for the gradient of the density (only allocated for GGA or mGGA)
 if ((ixc>=11.and.ixc<=17).or.(ixc>=23.and.ixc<=24).or.ixc==26.or.ixc==27.or. &
& (ixc>=31.and.ixc<=34)) ngr2=2*min(nspden,2)-1
 if (ixc<0.and.(libxc_functionals_isgga().or.libxc_functionals_ismgga())) &
&  ngr2=2*min(nspden,2)-1

!A-Only Exc and Vxc
 if (order**2 <= 1) then
   if (((ixc>=11 .and. ixc<=15) .or. (ixc>=23 .and. ixc<=24)) .and. ixc/=13) nvxcdgr=3
   if (ixc==16.or.ixc==17.or.ixc==26.or.ixc==27) nvxcdgr=2
   if (ixc<0) nvxcdgr=3
   if (ixc>=31 .and. ixc<=34) nvxcdgr=3 !Native fake metaGGA functionals (for testing purpose only)
 else

!B- Exc+Vxc and other derivatives
!  Definition of ndvxc and nvxcdgr, 2nd dimension of the arrays of 2nd-order derivatives
   if (ixc==1 .or. ixc==21 .or. ixc==22 .or. (ixc>=7 .and. ixc<=10) .or. ixc==13) then
!    Routine xcspol: new Teter fit (4/93) to Ceperley-Alder data, with spin-pol option routine xcspol
!    Routine xcpbe, with different options (optpbe) and orders (order)
     ndvxc=min(nspden,2)+1
   else if (ixc>=2 .and. ixc<=6) then
!    Perdew-Zunger fit to Ceperly-Alder data (no spin-pol)     !routine xcpzca
!    Teter fit (4/91) to Ceperley-Alder values (no spin-pol)   !routine xctetr
!    Wigner xc (no spin-pol)                                   !routine xcwign
!    Hedin-Lundqvist xc (no spin-pol)                          !routine xchelu
!    X-alpha (no spin-pol)                                     !routine xcxalp
     ndvxc=1
   else if (ixc==12 .or. ixc==24) then
!    Routine xcpbe, with optpbe=-2 and different orders (order)
     ndvxc=8
     nvxcdgr=3
   else if ((ixc>=11 .and. ixc<=15 .and. ixc/=13) .or. (ixc==23)) then
!    Routine xcpbe, with different options (optpbe) and orders (order)
     ndvxc=15
     nvxcdgr=3
   else if(ixc==16 .or. ixc==17 .or. ixc==26 .or. ixc==27 ) then
     ndvxc=0
     nvxcdgr=2
   else if (ixc<0) then
     if(libxc_functionals_isgga().or.libxc_functionals_ismgga()) then
       ndvxc=15
     else
       ndvxc=3
     end if
     nvxcdgr=3
   end if

!  Definition of nd2vxc, 2nd dimension of the array of 3rd-order derivatives
   if (order==3) then
     if (ixc==3) nd2vxc=1 ! Non spin polarized LDA case
     if ((ixc>=7 .and. ixc<=10) .or. (ixc==13)) nd2vxc=3*min(nspden,2)-2
!    Following line to be corrected when the calculation of d2vxcar is implemented for these functionals
     if ((ixc>=11 .and. ixc<=15 .and. ixc/=13) .or. (ixc==23.and.ixc<=24)) nd2vxc=1
     if ((ixc<0.and.(.not.(libxc_functionals_isgga().or. &
&                          libxc_functionals_ismgga())))) nd2vxc=3*min(nspden,2)-2
   end if

 end if

end subroutine pawxc_size_dvxc_local
!!***
#endif

end subroutine pawxc_size_dvxc_wrapper
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxc_xcmult_wrapper
!! NAME
!! pawxc_xcmult_wrapper
!!
!! FUNCTION
!! In the case of GGA, multiply the different gradient of spin-density
!! by the derivative of the XC functional with respect
!! to the norm of the gradient, then divide it by the
!! norm of the gradient
!!
!! INPUTS
!!  depsxc(nfft,nspgrad)=derivative of Exc with respect to the (spin-)density,
!!    or to the norm of the gradient of the (spin-)density,
!!    further divided by the norm of the gradient of the (spin-)density
!!   The different components of depsxc will be
!!   for nspden=1,         depsxc(:,1)=d(rho.exc)/d(rho)
!!         and if ngrad=2, depsxc(:,2)=1/2*1/|grad rho_up|*d(rho.exc)/d(|grad rho_up|)
!!                                      +   1/|grad rho|*d(rho.exc)/d(|grad rho|)
!!         (do not forget : |grad rho| /= |grad rho_up| + |grad rho_down|
!!   for nspden=2,         depsxc(:,1)=d(rho.exc)/d(rho_up)
!!                         depsxc(:,2)=d(rho.exc)/d(rho_down)
!!         and if ngrad=2, depsxc(:,3)=1/|grad rho_up|*d(rho.exc)/d(|grad rho_up|)
!!                         depsxc(:,4)=1/|grad rho_down|*d(rho.exc)/d(|grad rho_down|)
!!                         depsxc(:,5)=1/|grad rho|*d(rho.exc)/d(|grad rho|)
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngrad = must be 2
!!  nspden=number of spin-density components
!!  nspgrad=number of spin-density and spin-density-gradient components
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  rhonow(nfft,nspden,ngrad*ngrad)=
!!   at input :
!!    electron (spin)-density in real space and its gradient,
!!    either on the unshifted grid (if ishift==0,
!!      then equal to rhor), or on the shifted grid
!!     rhonow(:,:,1)=electron density in electrons/bohr**3
!!     rhonow(:,:,2:4)=gradient of electron density in el./bohr**4
!!   at output :
!!    rhonow(:,:,2:4) has been multiplied by the proper factor,
!!    described above.
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxc_xcmult_wrapper(depsxc,nfft,ngrad,nspden,nspgrad,rhonow)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ngrad,nspden,nspgrad
!arrays
 real(dp),intent(in) :: depsxc(nfft,nspgrad)
 real(dp),intent(inout) :: rhonow(nfft,nspden,ngrad*ngrad)

! *************************************************************************

#if defined HAVE_LIBPAW_ABINIT
 call xcmult(depsxc,nfft,ngrad,nspden,nspgrad,rhonow)
#else
 call pawxc_xcmult_local()
#endif
!!***

#if ! defined HAVE_LIBPAW_ABINIT
contains
!!***

!!****f* pawxc_xcmult_wrapper/pawxc_xcmult_local
!! NAME
!!  pawxc_xcmult_local
!!
!! FUNCTION
!!  Local version of xcmult routine (to use outside ABINIT)
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxc_xcmult_local()

!Local variables-------------------------------
!scalars
 integer :: idir,ifft
 real(dp) :: rho_tot,rho_up

! *************************************************************************

 do idir=1,3

   if(nspden==1)then
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(depsxc,idir,nfft,rhonow)
     do ifft=1,nfft
       rhonow(ifft,1,1+idir)=rhonow(ifft,1,1+idir)*depsxc(ifft,2)
     end do
   else
!    In the spin-polarized case, there are more factors to take into account
!$OMP PARALLEL DO PRIVATE(ifft,rho_tot,rho_up) SHARED(depsxc,idir,nfft,rhonow)
     do ifft=1,nfft
       rho_tot=rhonow(ifft,1,1+idir)
       rho_up =rhonow(ifft,2,1+idir)
       rhonow(ifft,1,1+idir)=rho_up *depsxc(ifft,3)         + rho_tot*depsxc(ifft,5)
       rhonow(ifft,2,1+idir)=(rho_tot-rho_up)*depsxc(ifft,4)+ rho_tot*depsxc(ifft,5)
     end do
   end if ! nspden==1
 end do ! End loop on directions

end subroutine pawxc_xcmult_local
!!***
#endif

end subroutine pawxc_xcmult_wrapper
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxc_mkdenpos_wrapper
!! NAME
!! pawxc_mkdenpos_wrapper
!!
!! FUNCTION
!! Make a ground-state density positive everywhere :
!! when the density (or spin-density) is smaller than xc_denpos,
!! set it to the value of xc_denpos
!!
!! INPUTS
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components (max. 2)
!!  option=0 if density rhonow is stored as (up,dn)
!!         1 if density rhonow is stored as (up+dn,up)
!!         Active only when nspden=2
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  Input/output
!!  iwarn=At input: iwarn=0 a warning will be printed when rho is negative
!!                  iwarn>0 no warning will be printed out
!!        At output: iwarn is increased by 1
!!  rhonow(nfft,nspden)=electron (spin)-density in real space,
!!     either on the unshifted grid (if ishift==0,
!!     then equal to rhor),or on the shifted grid
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxc_mkdenpos_wrapper(iwarn,nfft,nspden,option,rhonow,xc_denpos)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nspden,option
 integer,intent(inout) :: iwarn
 real(dp),intent(in) :: xc_denpos
!arrays
 real(dp),intent(inout) :: rhonow(nfft,nspden)

! *************************************************************************

#if defined HAVE_LIBPAW_ABINIT
 call mkdenpos(iwarn,nfft,nspden,option,rhonow,xc_denpos)
#else
 call pawxc_mkdenpos_local()
#endif
!!***

#if ! defined HAVE_LIBPAW_ABINIT
contains
!!***

!!****f* pawxc_mkdenpos_wrapper/pawxc_mkdenpos_local
!! NAME
!!  pawxc_mkdenpos_local
!!
!! FUNCTION
!!  Local version of mkdenpos routine (to use outside ABINIT)
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxc_mkdenpos_local()

!Local variables-------------------------------
!scalars
 integer :: ifft,ispden,numneg
 real(dp) :: rhotmp,worst
 character(len=500) :: msg
!arrays
 real(dp) :: rho(2)

! *************************************************************************

 numneg=0;worst=zero

 if(nspden==1)then
!  Non spin-polarized
!$OMP PARALLEL DO PRIVATE(ifft,rhotmp) REDUCTION(MIN:worst) REDUCTION(+:numneg) SHARED(nfft,rhonow)
   do ifft=1,nfft
     rhotmp=rhonow(ifft,1)
     if(rhotmp<xc_denpos)then
       if(rhotmp<-xc_denpos)then
!        This case is probably beyond machine precision considerations
         worst=min(worst,rhotmp)
         numneg=numneg+1
       end if
       rhonow(ifft,1)=xc_denpos
     end if
   end do

 else if (nspden==2) then
!  Spin-polarized

!  rhonow is stored as (up,dn)
   if (option==0) then
!$OMP PARALLEL DO PRIVATE(ifft,ispden,rho,rhotmp) REDUCTION(MIN:worst) REDUCTION(+:numneg) &
!$OMP&SHARED(nfft,nspden,rhonow)
     do ifft=1,nfft
!      For polarized case, rho(1) is spin-up density, rho(2) is spin-down density
       rho(1)=rhonow(ifft,1)
       rho(2)=rhonow(ifft,2)
       do ispden=1,nspden
         if (rho(ispden)<xc_denpos) then
           if (rho(ispden)<-xc_denpos) then
!            This case is probably beyond machine precision considerations
             worst=min(worst,rho(ispden))
             numneg=numneg+1
           end if
           rhonow(ifft,ispden)=xc_denpos
         end if
       end do
     end do

!  rhonow is stored as (up+dn,up)
   else if (option==1) then
!$OMP PARALLEL DO PRIVATE(ifft,ispden,rho,rhotmp) &
!$OMP&REDUCTION(MIN:worst) REDUCTION(+:numneg) &
!$OMP&SHARED(nfft,nspden,rhonow)
     do ifft=1,nfft
!      For polarized case, rho(1) is spin-up density, rho(2) is spin-down density
       rho(1)=rhonow(ifft,2)
       rho(2)=rhonow(ifft,1)-rho(1)
       do ispden=1,nspden
         if (rho(ispden)<xc_denpos) then
           if (rho(ispden)<-xc_denpos) then
!            This case is probably beyond machine precision considerations
             worst=min(worst,rho(ispden))
             numneg=numneg+1
           end if
           rho(ispden)=xc_denpos
           rhonow(ifft,1)=rho(1)+rho(2)
           rhonow(ifft,2)=rho(1)
         end if
       end do
     end do

   end if  ! option
 else
   msg='nspden>2 not allowed !'
   MSG_BUG(msg)
 end if ! End choice between non-spin polarized and spin-polarized.

 if (numneg>0) then
   if (iwarn==0) then
     write(msg,'(a,i10,a,a,a,es10.2,a,e10.2,a,a,a,a)')&
&     'Density went too small (lower than xc_denpos) at',numneg,' points',ch10,&
&     'and was set to xc_denpos=',xc_denpos,'.  Lowest was ',worst,'.',ch10,&
&     'Likely due to too low boxcut or too low ecut for','pseudopotential core charge.'
     MSG_WARNING(msg)
   end if
   iwarn=iwarn+1
 end if

end subroutine pawxc_mkdenpos_local
!!***
#endif

end subroutine pawxc_mkdenpos_wrapper
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxc_get_usekden
!! NAME
!!  pawxc_get_usekden
!!
!! FUNCTION
!!  Check if kinetic energy density has to be computed
!!
!! INPUTS
!!  ixc= choice of exchange-correlation scheme
!!
!! SOURCE

function pawxc_get_usekden()
!Arguments ------------------------------------
  integer :: pawxc_get_usekden

! *************************************************************************

  pawxc_get_usekden=0
  if (libxc_functionals_ismgga()) pawxc_get_usekden=1

end function pawxc_get_usekden
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxc
!! NAME
!! pawxc
!!
!! FUNCTION
!! Start from the density or spin-density, and compute xc correlation
!! potential and energies inside a paw sphere.
!! USE THE DENSITY OVER A WHOLE SPHERICAL GRID (r,theta,phi)
!! Driver of XC functionals.
!!
!! INPUTS
!!  corexc(nrad)=core density on radial grid
!!  ixc= choice of exchange-correlation scheme
!!  lm_size=size of density array rhor (see below)
!!  lmselect(lm_size)=select the non-zero LM-moments of input density rhor
!!  nhat(nrad,lm_size,nspden)=compensation density
!!                                        (total in 1st half and spin-up in 2nd half if nspden=2)
!!  nkxc=second dimension of the kxc array. If /=0, the exchange-correlation kernel must be computed
!!  non_magnetic_xc= if true, handle density/potential as non-magnetic (even if it is)
!!  nrad=size of radial mesh for densities/potentials (might be different from pawrad%mesh_size)
!!  nspden=number of spin-density components
!!  option=0  compute both XC energies (direct+double-counting) and potential
!!         1  compute only XC potential
!!         2  compute only XC energies (direct+double-counting)
!!         3  compute only XC energy by direct scheme
!!         4  compute only XC energy by direct scheme for spherical part of the density
!!         5  compute only XC potential for spherical part of the density
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  rhor(nrad,lm_size,nspden)=electron density in real space in electrons/bohr**3
!!                                       (total in 1st half and spin-up in 2nd half if nspden=2)
!!  coretau(nrad*usekden) = core kinetic energy density
!!  taur(nrad,lm_size,nspden*usekden) = kinetic energy density on radial mesh
!!  usecore= 1 if core density has to be used in Exc/Vxc ; 0 otherwise
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in double counting energy term only
!!             2 if compensation density (nhat) has to be used in Exc/Vxc and double counting energy term
!!  xclevel= XC functional level
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!
!! OUTPUT
!!  == if option=0, 2, 3, or 4 ==
!!    enxc=returned exchange and correlation energy (hartree)
!!  == if option=0 or 2 ==
!!    enxcdc=returned exchange-cor. contribution to double-counting energy
!!  == if option=0, 1 or 5 ==
!!    vxc(nrad,pawang%angl_size,nspden)=xc potential
!!       (spin up in 1st half and spin-down in 2nd half if nspden=2)
!!   == if option=0, 1 or 5 and usekden=1 ==
!!    vxctau(nrad,pawang%angl_size,nspden*usekden)=MetaGGA xc potential
!!       (spin up in 1st half and spin-down in 2nd half if nspden=2)
!!  == if nkxc>0 ==
!!    kxc(nrad,pawang%angl_size,nkxc)=xc kernel
!!        (see notes below for nkxc)
!!  == if nk3xc>0 ==
!!    k3xc(nrad,pawang%angl_size,nk3xc)= derivative of xc kernel
!!        (see notes below for nk3xc)
!!
!! NOTES
!!  Content of Kxc array:
!!   ===== if LDA
!!    if nspden==1: kxc(:,1)= d2Exc/drho2
!!                 (kxc(:,2)= d2Exc/drho_up drho_dn)
!!    if nspden>=2: kxc(:,1)= d2Exc/drho_up drho_up
!!                  kxc(:,2)= d2Exc/drho_up drho_dn
!!                  kxc(:,3)= d2Exc/drho_dn drho_dn
!!    if nspden==4: kxc(:,4:6)= (m_x, m_y, m_z) (magnetization)
!!   ===== if GGA
!!    if nspden==1:
!!       kxc(:,1)= d2Exc/drho2
!!       kxc(:,2)= 1/|grad(rho)| dExc/d|grad(rho)|
!!       kxc(:,3)= 1/|grad(rho)| d2Exc/d|grad(rho)| drho
!!       kxc(:,4)= 1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dExc/d|grad(rho)| )
!!       kxc(:,5)= gradx(rho)
!!       kxc(:,6)= grady(rho)
!!       kxc(:,7)= gradz(rho)
!!    if nspden>=2:
!!       kxc(:,1)= d2Exc/drho_up drho_up
!!       kxc(:,2)= d2Exc/drho_up drho_dn
!!       kxc(:,3)= d2Exc/drho_dn drho_dn
!!       kxc(:,4)= 1/|grad(rho_up)| dEx/d|grad(rho_up)|
!!       kxc(:,5)= 1/|grad(rho_dn)| dEx/d|grad(rho_dn)|
!!       kxc(:,6)= 1/|grad(rho_up)| d2Ex/d|grad(rho_up)| drho_up
!!       kxc(:,7)= 1/|grad(rho_dn)| d2Ex/d|grad(rho_dn)| drho_dn
!!       kxc(:,8)= 1/|grad(rho_up)| * d/d|grad(rho_up)| ( 1/|grad(rho_up)| dEx/d|grad(rho_up)| )
!!       kxc(:,9)= 1/|grad(rho_dn)| * d/d|grad(rho_dn)| ( 1/|grad(rho_dn)| dEx/d|grad(rho_dn)| )
!!       kxc(:,10)=1/|grad(rho)| dEc/d|grad(rho)|
!!       kxc(:,11)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_up
!!       kxc(:,12)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_dn
!!       kxc(:,13)=1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dEc/d|grad(rho)| )
!!       kxc(:,14)=gradx(rho_up)
!!       kxc(:,15)=gradx(rho_dn)
!!       kxc(:,16)=grady(rho_up)
!!       kxc(:,17)=grady(rho_dn)
!!       kxc(:,18)=gradz(rho_up)
!!       kxc(:,19)=gradz(rho_dn)
!!    if nspden==4:
!!       kxc(:,20:22)= (m_x, m_y, m_z) (magnetization)
!!  Dimension of K3xc:
!!   ===== if LDA (xclevel=1) :
!!    if nspden==1: return  k3xc(:,1)=d3Exc/drho3
!!    if nspden>=2, return  k3xc(:,1)=d3Exc/drho_up drho_up drho_up
!!                          k3xc(:,2)=d3Exc/drho_up drho_up drho_dn
!!                          k3xc(:,3)=d3Exc/drho_up drho_dn drho_dn
!!                          k3xc(:,4)=d3Exc/drho_dn drho_dn drho_dn
!!
!! PARENTS
!!      m_pawpsp,pawdenpot
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE
subroutine pawxc(corexc,enxc,enxcdc,ixc,kxc,k3xc,lm_size,lmselect,nhat,nkxc,nk3xc,non_magnetic_xc,&
&                nrad,nspden,option,pawang,pawrad,rhor,usecore,usexcnhat,vxc,xclevel,xc_denpos,&
&                coretau,taur,vxctau) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,lm_size,nkxc,nk3xc,nrad,nspden,option,usecore,usexcnhat,xclevel
 logical,intent(in) :: non_magnetic_xc
 real(dp),intent(in) :: xc_denpos
 real(dp),intent(out) :: enxc,enxcdc
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size)
 real(dp),intent(in) :: corexc(nrad)
 real(dp),intent(in) :: nhat(nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in),target :: rhor(nrad,lm_size,nspden)
 real(dp),intent(in),target,optional:: coretau(:),taur(:,:,:)
 real(dp),intent(out) :: kxc(nrad,pawang%angl_size,nkxc)
 real(dp),intent(out) :: k3xc(nrad,pawang%angl_size,nk3xc)
 real(dp),intent(out),target :: vxc(nrad,pawang%angl_size,nspden)
 real(dp),intent(out),target,optional :: vxctau(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,ilm,ipts,ir,ispden,iwarn,lm_size_eff,mgga,ndvxc,nd2vxc,ngr2,ngrad
 integer :: nkxc_updn,npts,nspden_eff,nspden_updn,nspgrad,nvxcdgr,order
 integer :: usecoretau,usekden,use_laplacian
 logical :: ismgga,with_vxctau
 real(dp) :: enxcr,factor,vxcrho,factor2
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: dgxc(:),dnexcdn(:,:),drho(:),d2rho(:),drhocore(:),dvxcdgr(:,:)
 real(dp),allocatable :: dvxcdlr(:,:),dvxci(:,:),d2vxci(:,:),dylmdr(:,:,:),d2ylmdr(:,:,:)
 real(dp),allocatable :: exci(:),ff(:),grho2_updn(:,:),gxc(:,:,:,:)
 real(dp),allocatable :: rhoarr(:,:),rho_updn(:,:),lrho_updn(:,:),tauarr(:,:)
 real(dp),allocatable :: tau_updn(:,:),vxci(:,:)
 real(dp),allocatable,target :: mag(:,:,:),rhohat(:,:,:),rhonow(:,:,:),taunow(:,:)
 real(dp), pointer :: mag_(:,:),rho_(:,:,:),tau_(:,:,:)
 real(dp), LIBPAW_CONTIGUOUS pointer :: vxc_diag(:,:),vxc_nc(:,:),vxc_updn(:,:,:)
#ifdef LIBPAW_ISO_C_BINDING
 type(C_PTR) :: cptr
#endif

! *************************************************************************

!----------------------------------------------------------------------
!----- Check options
!----------------------------------------------------------------------

!Some flags
 nkxc_updn=merge(nkxc-3,nkxc,nkxc==6.or.nkxc==22)
 ismgga=libxc_functionals_ismgga()
 mgga=merge(1,0,ismgga)
 usekden=pawxc_get_usekden()
 use_laplacian=mgga ! TEMPORARY
 usecoretau=0 ; with_vxctau=.false.

!Compatibility tests
!if (mgga==1.and.usekden==0) then
!   msg='Kinetic energy density needs to be computed'
!   MSG_ERROR(msg)
! end if
 if (usekden==1) then
   if (.not.present(taur)) then
     msg='taur needs to be present!'
     MSG_BUG(msg)
   else if (size(taur)/=nrad*lm_size*nspden) then
     msg='wrong size for taur!'
     MSG_BUG(msg)
   end if
   if (present(coretau)) then
     usecoretau=1
     if (size(coretau)/=nrad) then
       msg='wrong size for coretau!'
       MSG_BUG(msg)
     end if
   end if
   if ((option==0.or.option==1.or.option==5).and.present(vxctau)) then
     with_vxctau=.true.
     if (size(vxctau)/=nrad*pawang%angl_size*nspden) then
       msg='wrong size for vxctau!'
       MSG_BUG(msg)
     end if
   end if
   !?????For call in m_pawpsp/pawpsp_calc do not compute density laplacian ?????
   if (option==4) use_laplacian=0
 end if
 if(nspden==4.and.nk3xc>0) then
   msg='K3xc for nspden=4 not implemented!'
   MSG_ERROR(msg)
 end if
 if(nk3xc>0.and.nkxc_updn==0) then
   msg='nkxc must be non-zero if nk3xc is!'
   MSG_ERROR(msg)
 end if
 if(nspden==4.and.xclevel==2) then
   msg='GGA for nspden=4 not implemented!'
   MSG_ERROR(msg)
 end if
 if(pawang%angl_size==0) then
   msg='pawang%angl_size=0!'
   MSG_BUG(msg)
 end if
 if(.not.allocated(pawang%ylmr)) then
   msg='pawang%ylmr must be allocated!'
   MSG_BUG(msg)
 end if
 if(xclevel==2.and.(.not.allocated(pawang%ylmrgr))) then
   msg='pawang%ylmrgr must be allocated!'
   MSG_BUG(msg)
 end if
 if(option==4.or.option==5) then
   if (pawang%angl_size/=1) then
     msg='When option=4 or 5, pawang%angl_size must be 1!'
     MSG_BUG(msg)
   end if
   if (pawang%ylm_size/=1) then
     msg='When option=4 or 5, pawang%ylm_size must be 1!'
     MSG_BUG(msg)
   end if
   if (abs(pawang%anginit(1,1)-one)>tol12.or.abs(pawang%anginit(2,1))>tol12.or. &
&   abs(pawang%anginit(3,1))>tol12) then
     msg='When option=4 or 5, pawang%anginit must be (1 0 0)!'
     MSG_BUG(msg)
   end if
 end if
 if (option/=1.and.option/=5) then
   if (nrad<pawrad%int_meshsz) then
     msg='When option=0,2,3,4, nrad must be greater than pawrad%int_meshsz!'
     MSG_BUG(msg)
   end if
 end if

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------
 iwarn=0
 nspden_updn=min(nspden,2)
 nspden_eff=nspden_updn;if (nspden==4.and.xclevel==2) nspden_eff=4
 npts=pawang%angl_size
 lm_size_eff=min(lm_size,pawang%ylm_size)
 ngrad=1;if(xclevel==2)ngrad=2
 nspgrad=0;if (xclevel==2) nspgrad=3*nspden_updn-1
 if (option/=1.and.option/=5) enxc=zero
 if (option==0.or.option==2) enxcdc=zero
 if (option/=3.and.option/=4) vxc(:,:,:)=zero
 if (nkxc>0) kxc(:,:,:)=zero
 if (nk3xc>0) k3xc(:,:,:)=zero
 order=1;if (nkxc_updn>0) order=2;if (nk3xc>0) order=3 ! to which der. of the energy the computation must be done

 if (xclevel==0.or.ixc==0) then
   msg='Note that no xc is applied (ixc=0).'
   MSG_WARNING(msg)

 else

!  Allocation of temporary memory space
   LIBPAW_ALLOCATE(rhonow,(nrad,nspden,ngrad*ngrad+1))
   LIBPAW_ALLOCATE(rhoarr,(nrad,nspden))
   if (usexcnhat>0) then
     LIBPAW_ALLOCATE(rhohat,(nrad,lm_size,nspden))
     rhohat(:,:,:)=rhor(:,:,:)+nhat(:,:,:)
   end if
   if (xclevel==2.and.usecore==1) then
     LIBPAW_ALLOCATE(drhocore,(nrad))
     call nderiv_gen(drhocore,corexc,pawrad)
   end if
   if (option/=3.and.option/=4) then
     if (nspden/=4) then
       vxc_updn => vxc
     else
       LIBPAW_POINTER_ALLOCATE(vxc_updn,(nrad,npts,nspden_updn))
       LIBPAW_ALLOCATE(mag,(nrad,npts,3))
     end if
   end if

!  Allocation of mandatory arguments of drivexc
   LIBPAW_ALLOCATE(exci,(nrad))
   LIBPAW_ALLOCATE(vxci,(nrad,nspden_updn))
   LIBPAW_ALLOCATE(rho_updn,(nrad,nspden_updn))
   LIBPAW_ALLOCATE(lrho_updn,(nrad,nspden_updn))
!  Allocation of optional arguments of drivexc
   call pawxc_size_dvxc_wrapper(ixc,ndvxc,ngr2,nd2vxc,nspden_updn,nvxcdgr,order)
   LIBPAW_ALLOCATE(dvxci,(nrad,ndvxc))
   LIBPAW_ALLOCATE(d2vxci,(nrad,nd2vxc))
   LIBPAW_ALLOCATE(dvxcdgr,(nrad,nvxcdgr))
   LIBPAW_ALLOCATE(dvxcdlr,(nrad,nvxcdgr))
   LIBPAW_ALLOCATE(grho2_updn,(nrad,ngr2))
   LIBPAW_ALLOCATE(dnexcdn,(nrad,nspgrad))

!  GGA: convert Ylm derivatives from normalized to standard cartesian coordinates
!  dYlm/dr_i = { dYlm/dr_i^hat - Sum_j[ dYlm/dr_j^hat (r_j/r)] } * (1/r)
!  Meta GGA: convert Ylm second derivatives from normalized to standard cartesian coordinates
!  d2Ylm/dr_i = { d2Ylm/d2r_i^hat - Sum_j[ d2Ylm/d2r_j^hat (r_j/r)] } * (1/r)
   if (xclevel==2) then
     LIBPAW_ALLOCATE(dylmdr,(3,npts,pawang%ylm_size))
     if (mgga==1) then
       !mgga=1
       LIBPAW_ALLOCATE(tauarr,(nrad,nspden))
       LIBPAW_ALLOCATE(taunow,(nrad,nspden))
       LIBPAW_ALLOCATE(tau_updn,(nrad,nspden_updn))
       if (use_laplacian==1) then
         LIBPAW_ALLOCATE(d2ylmdr,(3,npts,pawang%ylm_size))
       end if
     end if
     do ilm=1,pawang%ylm_size
       do ipts=1,npts
         factor=sum(pawang%ylmrgr(1:3,ilm,ipts)*pawang%anginit(1:3,ipts))
         dylmdr(1:3,ipts,ilm)=pawang%ylmrgr(1:3,ilm,ipts)-factor*pawang%anginit(1:3,ipts)
         if (use_laplacian==1) then
           factor2=pawang%ylmrgr(4,ilm,ipts)*pawang%anginit(1,ipts)+pawang%ylmrgr(5,ilm,ipts)*pawang%anginit(2,ipts)&
&                 +pawang%ylmrgr(6,ilm,ipts)*pawang%anginit(3,ipts)
           d2ylmdr(1,ipts,ilm)=pawang%ylmrgr(4,ilm,ipts)-factor2*pawang%anginit(1,ipts)
           d2ylmdr(2,ipts,ilm)=pawang%ylmrgr(5,ilm,ipts)-factor2*pawang%anginit(2,ipts)
           d2ylmdr(3,ipts,ilm)=pawang%ylmrgr(6,ilm,ipts)-factor2*pawang%anginit(3,ipts)
         end if
       end do
     end do
     LIBPAW_ALLOCATE(gxc,(nrad,3,pawang%ylm_size,nspden_updn))
     gxc=zero
   end if

!  ----------------------------------------------------------------------
!  ----- Loop on the angular part and inits
!  ----------------------------------------------------------------------

!  Do loop on the angular part
   do ipts=1,npts

!    Copy the input density for this (theta,phi)
     rhoarr(:,:)=zero
     if (mgga==1) then
       tauarr(:,:)=zero
     end if
     if (usexcnhat< 2) rho_=> rhor
     if (usexcnhat==2) rho_=> rhohat
     if (mgga==1) tau_=> taur
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect(ilm)) then
           rhoarr(1:nrad,ispden)=rhoarr(1:nrad,ispden) &
&           +rho_(1:nrad,ilm,ispden)*pawang%ylmr(ilm,ipts)
           if (mgga==1) then
             tauarr(1:nrad,ispden)=tauarr(1:nrad,ispden) &
&           +tau_(1:nrad,ilm,ispden)*pawang%ylmr(ilm,ipts)
           end if
         end if
       end do
     end do
     if (usecore==1) then
       rhoarr(1:nrad,1)=rhoarr(1:nrad,1)+corexc(1:nrad)
       if (mgga==1) then
         tauarr(1:nrad,1)=tauarr(1:nrad,1)+coretau(1:nrad)
       end if
       if (nspden==2) rhoarr(1:nrad,2)=rhoarr(1:nrad,2)+half*corexc(1:nrad);
       if (nspden==2.and.mgga==1) tauarr(1:nrad,2)=tauarr(1:nrad,2)+half*coretau(1:nrad);
     end if

!    Optionally suppress magnetic part.
     if (non_magnetic_xc) then
       if(nspden==2) rhoarr(:,2)=rhoarr(:,1)*half
       if(nspden==4) rhoarr(:,2:4)=zero
     endif

     rhonow(1:nrad,1:nspden,1)=rhoarr(1:nrad,1:nspden)

!    GGA: compute gradient of density
     if (xclevel==2) then
       rhonow(:,:,2:5)=zero
       LIBPAW_ALLOCATE(drho,(nrad))
       LIBPAW_ALLOCATE(ff,(nrad))
       do ispden=1,nspden
         do ilm=1,lm_size_eff
           if (lmselect(ilm)) then
             ff(1:nrad)=rho_(1:nrad,ilm,ispden)
             call nderiv_gen(drho,ff,pawrad)
             ff(2:nrad)=ff(2:nrad)/pawrad%rad(2:nrad)
             call pawrad_deducer0(ff,nrad,pawrad)
             do ii=1,3
               rhonow(1:nrad,ispden,1+ii)=rhonow(1:nrad,ispden,1+ii) &
&               +drho(1:nrad)*pawang%ylmr(ilm,ipts)*pawang%anginit(ii,ipts) &
&               +ff(1:nrad)*dylmdr(ii,ipts,ilm)
             end do
           end if
         end do
       end do
       if (non_magnetic_xc) then
         do ii=1,3
           if(nspden==2) rhonow(1:nrad,2,1+ii)=rhonow(1:nrad,1,1+ii)*half
           if(nspden==4) rhonow(1:nrad,2:4,1+ii)=zero
         end do
       end if

!      Compute the laplacian of the density
       if (use_laplacian==1) then
         LIBPAW_ALLOCATE(d2rho,(nrad))
         call nderiv_gen(d2rho,drho,pawrad)
         do ispden=1,nspden
           do ilm=1,lm_size_eff
             ff(1:nrad)=rho_(1:nrad,ilm,ispden)
             call nderiv_gen(drho,ff,pawrad,d2rho)
             ff(2:nrad)=ff(2:nrad)/pawrad%rad(2:nrad)
             call pawrad_deducer0(ff,nrad,pawrad)
             rhonow(1:nrad,ispden,5)=rhonow(1:nrad,ispden,5) &
&                                    +d2rho(1:nrad)*pawang%ylmr(ilm,ipts)
             do ii=1,3
               rhonow(1:nrad,ispden,5)=rhonow(1:nrad,ispden,5) &
&                                      +drho(1:nrad)*pawang%anginit(ii,ipts)*dylmdr(ii,ipts,ilm)&
&                                      +ff(1:nrad)*d2ylmdr(ii,ipts,ilm)
              end do
           end do
         end do
         LIBPAW_DEALLOCATE(d2rho)
       end if
       LIBPAW_DEALLOCATE(drho)
       LIBPAW_DEALLOCATE(ff)
       if (usecore==1) then
         do ii=1,3
           rhonow(1:nrad,1,1+ii)=rhonow(1:nrad,1,1+ii) &
&           +drhocore(1:nrad)*pawang%anginit(ii,ipts)
         end do
         if (nspden==2) then
           do ii=1,3
             rhonow(1:nrad,2,1+ii)=rhonow(1:nrad,2,1+ii) &
&             +half*drhocore(1:nrad)*pawang%anginit(ii,ipts)
           end do
         end if
       end if
     end if

!    Storage of density (and gradient) in (up,dn) format
     if (nspden==1) then
       rho_updn(1:nrad,1)=rhonow(1:nrad,1,1)*half
       if (xclevel==2) then
         grho2_updn(1:nrad,1)=quarter*(rhonow(1:nrad,1,2)**2+rhonow(1:nrad,1,3)**2+rhonow(1:nrad,1,4)**2)
         if (mgga==1) then
           if (use_laplacian==1) then
             lrho_updn(1:nrad,1)=rhonow(1:nrad,1,5)*half
           end if
           tau_updn(1:nrad,1)=tauarr(1:nrad,1)*half
         end if
       end if
     else if (nspden==2) then
       rho_updn(1:nrad,1)=rhonow(1:nrad,2,1)
       rho_updn(1:nrad,2)=rhonow(1:nrad,1,1)-rhonow(1:nrad,2,1)
       if (xclevel==2) then
         grho2_updn(1:nrad,1)=rhonow(1:nrad,2,2)**2+rhonow(1:nrad,2,3)**2+rhonow(1:nrad,2,4)**2
         grho2_updn(1:nrad,2)=(rhonow(1:nrad,1,2)-rhonow(1:nrad,2,2))**2 +   &
&         (rhonow(1:nrad,1,3)-rhonow(1:nrad,2,3))**2 +   &
&         (rhonow(1:nrad,1,4)-rhonow(1:nrad,2,4))**2
         grho2_updn(1:nrad,3)=rhonow(1:nrad,1,2)**2+rhonow(1:nrad,1,3)**2+rhonow(1:nrad,1,4)**2
         if (mgga==1) then
           if (use_laplacian==1) then
             lrho_updn(1:nrad,1)=rhonow(1:nrad,2,5)
             lrho_updn(1:nrad,2)=rhonow(1:nrad,1,5)-rhonow(1:nrad,2,5)
           end if
           tau_updn(1:nrad,1)=tauarr(1:nrad,2)
           tau_updn(1:nrad,2)=tauarr(1:nrad,1)-tauarr(1:nrad,2)
         end if
       end if
     else if (nspden==4) then
       mag_ => rhonow(1:nrad,2:4,1)
       mag(1:nrad,ipts,1:3)=mag_(1:nrad,1:3)
       call pawxc_rotate_mag(rhonow(:,:,1),rho_updn,mag_,nrad)
     end if

!    Make the density positive everywhere (but do not care about gradients)
     call pawxc_mkdenpos_wrapper(iwarn,nrad,nspden_updn,0,rho_updn,xc_denpos)

!    Call to main XC driver

!?????? utiliser with_vxctau
     call pawxc_drivexc_wrapper(exci,ixc,mgga,ndvxc,nd2vxc,ngr2,nrad,nspden_updn,nvxcdgr,&
&                               order,rho_updn,use_laplacian,vxci,xclevel, &
&                               dvxc=dvxci,d2vxc=d2vxci,grho2=grho2_updn,lrho=lrho_updn, &
&                               tau=tau_updn,vxctau=vxctau,vxcgrho=dvxcdgr,vxclrho=dvxcdlr)

!    ----------------------------------------------------------------------
!    ----- Accumulate and store XC kernel and its derivative
!    ----------------------------------------------------------------------
     if (nkxc_updn>0.and.ndvxc>0) then
       if (nkxc_updn==1.and.ndvxc==15) then
         kxc(1:nrad,ipts,1)=half*(dvxci(1:nrad,1)+dvxci(1:nrad,9)+dvxci(1:nrad,10))
       else if (nkxc_updn==3.and.ndvxc==15) then
         kxc(1:nrad,ipts,1)=dvxci(1:nrad,1)+dvxci(1:nrad,9)
         kxc(1:nrad,ipts,2)=dvxci(1:nrad,10)
         kxc(1:nrad,ipts,3)=dvxci(1:nrad,2)+dvxci(1:nrad,11)
       else if (nkxc_updn==7.and.ndvxc==8) then
         kxc(1:nrad,ipts,1)=half*dvxci(1:nrad,1)
         kxc(1:nrad,ipts,2)=half*dvxci(1:nrad,3)
         kxc(1:nrad,ipts,3)=quarter*dvxci(1:nrad,5)
         kxc(1:nrad,ipts,4)=eighth*dvxci(1:nrad,7)
       else if (nkxc_updn==7.and.ndvxc==15) then
         kxc(1:nrad,ipts,1)=half*(dvxci(1:nrad,1)+dvxci(1:nrad,9)+dvxci(1:nrad,10))
         kxc(1:nrad,ipts,2)=half*dvxci(1:nrad,3)+dvxci(1:nrad,12)
         kxc(1:nrad,ipts,3)=quarter*dvxci(1:nrad,5)+dvxci(1:nrad,13)
         kxc(1:nrad,ipts,4)=eighth*dvxci(1:nrad,7)+dvxci(1:nrad,15)
       else if (nkxc_updn==19.and.ndvxc==15) then
         kxc(1:nrad,ipts,1)=dvxci(1:nrad,1)+dvxci(1:nrad,9)
         kxc(1:nrad,ipts,2)=dvxci(1:nrad,10)
         kxc(1:nrad,ipts,3)=dvxci(1:nrad,2)+dvxci(1:nrad,11)
         kxc(1:nrad,ipts,4)=dvxci(1:nrad,3)
         kxc(1:nrad,ipts,5)=dvxci(1:nrad,4)
         kxc(1:nrad,ipts,6)=dvxci(1:nrad,5)
         kxc(1:nrad,ipts,7)=dvxci(1:nrad,6)
         kxc(1:nrad,ipts,8)=dvxci(1:nrad,7)
         kxc(1:nrad,ipts,9)=dvxci(1:nrad,8)
         kxc(1:nrad,ipts,10)=dvxci(1:nrad,12)
         kxc(1:nrad,ipts,11)=dvxci(1:nrad,13)
         kxc(1:nrad,ipts,12)=dvxci(1:nrad,14)
         kxc(1:nrad,ipts,13)=dvxci(1:nrad,15)
       else ! Other cases
         kxc(1:nrad,ipts,1:nkxc)=zero
         kxc(1:nrad,ipts,1:min(nkxc,ndvxc))=dvxci(1:nrad,1:min(nkxc,ndvxc))
       end if
       if (nkxc_updn==7) then
         kxc(1:nrad,ipts,5)=rhonow(1:nrad,1,2)
         kxc(1:nrad,ipts,6)=rhonow(1:nrad,1,3)
         kxc(1:nrad,ipts,7)=rhonow(1:nrad,1,4)
       else if (nkxc_updn==19) then
         kxc(1:nrad,ipts,14)=rhonow(1:nrad,1,2)
         kxc(1:nrad,ipts,15)=rhonow(1:nrad,2,2)
         kxc(1:nrad,ipts,16)=rhonow(1:nrad,1,3)
         kxc(1:nrad,ipts,17)=rhonow(1:nrad,2,3)
         kxc(1:nrad,ipts,18)=rhonow(1:nrad,1,4)
         kxc(1:nrad,ipts,19)=rhonow(1:nrad,2,4)
       end if
     end if
     if (nkxc>=nkxc_updn+3) then
       kxc(1:nrad,ipts,nkxc_updn+1)=rhonow(1:nrad,2,1)
       kxc(1:nrad,ipts,nkxc_updn+2)=rhonow(1:nrad,3,1)
       kxc(1:nrad,ipts,nkxc_updn+3)=rhonow(1:nrad,4,1)
     end if

!    kernel derivative :
     if (nk3xc>0.and.nd2vxc>0) then
       k3xc(1:nrad,ipts,1:min(nk3xc,nd2vxc))=d2vxci(1:nrad,1:min(nk3xc,nd2vxc))
     end if

!    ----------------------------------------------------------------------
!    ----- Accumulate and store XC potential
!    ----------------------------------------------------------------------

     if (option/=3.and.option/=4) then

       do ispden=1,nspden_updn
         vxc_updn(1:nrad,ipts,ispden)=vxci(1:nrad,ispden)
       end do

!      For GGAs, additional terms appear
       if (xclevel==2.and.ixc/=13)then
         dnexcdn(1:nrad,1:nspden_updn)=vxci(1:nrad,1:nspden_updn)
!        Treat explicitely spin up, spin down and total spin for spin-polarized
         do ii=1,3
           if(nspden_updn==1.and.ii>=2)exit !exit when ii=1 is finished if non-spin-polarized
           do ir=1,nrad
!            If the norm of the gradient vanishes, then the different terms vanishes
             if(grho2_updn(ir,ii)<1.0d-24) then
               dnexcdn(ir,ii+nspden_updn)=zero;cycle
             end if
!            Compute the derivative of n.e_xc wrt spin up, spin down, or total density
             if(nspden_updn==1)then
               dnexcdn(ir,ii+nspden_updn)=half*dvxcdgr(ir,1) !Definition of dvxcdgr changed in v3.3
               if (nvxcdgr==3) dnexcdn(ir,ii+nspden_updn)=dnexcdn(ir,ii+nspden_updn)+dvxcdgr(ir,3)
             else if(nspden_updn==2)then
               if (nvxcdgr==3) then
                 dnexcdn(ir,ii+nspden_updn)=dvxcdgr(ir,ii)
               else if (ii/=3) then
                 dnexcdn(ir,ii+nspden_updn)=dvxcdgr(ir,ii)
               else if (ii==3) then
                 dnexcdn(ir,ii+nspden_updn)=zero
               end if
             end if
           end do
         end do
         call pawxc_xcmult_wrapper(dnexcdn,nrad,ngrad,nspden_eff,nspgrad,rhonow)
         factor=one;if (nspden_updn==1) factor=half
         if (option/=4.and.option/=5) then
           factor=factor*four_pi
!          Accumulate moments of gxc
           do ispden=1,nspden_updn
             do ilm=1,pawang%ylm_size
               do ii=1,3
                 gxc(1:nrad,ii,ilm,ispden)=gxc(1:nrad,ii,ilm,ispden)+rhonow(1:nrad,ispden,1+ii) &
&                 *pawang%ylmr(ilm,ipts)*pawang%angwgth(ipts)*factor
               end do
             end do
           end do
         else
           do ispden=1,nspden_updn
             gxc(1:nrad,1,1,ispden)=factor*rhonow(1:nrad,ispden,2)
           end do
         end if
       end if

     end if !option

!    ----------------------------------------------------------------------
!    ----- Accumulate and store XC energy
!    ----------------------------------------------------------------------
     if (option/=1.and.option/=5) then
       LIBPAW_ALLOCATE(ff,(nrad))
       ff(1:nrad)=rhoarr(1:nrad,1)*exci(1:nrad)*pawrad%rad(1:nrad)**2
       call simp_gen(enxcr,ff,pawrad)
       if (option/=4) enxc=enxc+enxcr*pawang%angwgth(ipts)
       if (option==4) enxc=enxc+enxcr
       LIBPAW_DEALLOCATE(ff)
     end if

!    ----------------------------------------------------------------------
!    ----- End of the loop on npts (angular part)
!    ----------------------------------------------------------------------
   end do

!  Deallocate temporary memory space
   LIBPAW_DEALLOCATE(exci)
   LIBPAW_DEALLOCATE(vxci)
   LIBPAW_DEALLOCATE(rho_updn)
   LIBPAW_DEALLOCATE(lrho_updn)
   LIBPAW_DEALLOCATE(dvxci)
   LIBPAW_DEALLOCATE(d2vxci)
   LIBPAW_DEALLOCATE(dvxcdgr)
   LIBPAW_DEALLOCATE(dvxcdlr)
   LIBPAW_DEALLOCATE(grho2_updn)
   LIBPAW_DEALLOCATE(dnexcdn)
   if (xclevel==2) then
     if (mgga==1) then
       LIBPAW_DEALLOCATE(tau_updn)
     end if
     if (usecore==1) then
       LIBPAW_DEALLOCATE(drhocore)
     end if
   end if
   LIBPAW_DEALLOCATE(rhonow)

!  ----------------------------------------------------------------------
!  ----- If GGA, modify potential with term from density gradient
!  ----------------------------------------------------------------------
   if (option/=3.and.option/=4.and.xclevel==2.and.ixc/=13) then
!    Compute divergence of gxc and substract it from Vxc
     LIBPAW_ALLOCATE(dgxc,(nrad))
!    Need to multiply gxc by 2 in the non-polarised case
     factor=one;if (nspden_updn==1) factor=two
     if (option/=4.and.option/=5) then
       LIBPAW_ALLOCATE(ff,(nrad))
       do ispden=1,nspden_updn
         do ilm=1,pawang%ylm_size
           do ii=1,3
             ff(1:nrad)=gxc(1:nrad,ii,ilm,ispden)
             call nderiv_gen(dgxc,ff,pawrad)
             ff(2:nrad)=ff(2:nrad)/pawrad%rad(2:nrad)
             call pawrad_deducer0(ff,nrad,pawrad)
             do ipts=1,npts
               vxc_updn(1:nrad,ipts,ispden)=vxc_updn(1:nrad,ipts,ispden) &
&               -factor*(dgxc(1:nrad)*pawang%anginit(ii,ipts)*pawang%ylmr(ilm,ipts) &
&               +ff(1:nrad)*dylmdr(ii,ipts,ilm))
             end do
           end do
         end do
       end do
       LIBPAW_DEALLOCATE(ff)
     else ! option==4 or option==5
       do ispden=1,nspden_updn
         call nderiv_gen(dgxc,gxc(:,1,1,ispden),pawrad)
         vxc_updn(2:nrad,1,ispden)=vxc_updn(2:nrad,1,ispden) &
&         -factor*(dgxc(2:nrad)+two*gxc(2:nrad,1,1,ispden)/pawrad%rad(2:nrad))
         call pawrad_deducer0(vxc(:,1,ispden),nrad,pawrad)
       end do
     end if
     LIBPAW_DEALLOCATE(dgxc)
   end if ! GGA

!  ----------------------------------------------------------------------
!  ----- If non-collinear, rotate back potential according to magnetization
!  ----------------------------------------------------------------------
   if (option/=3.and.option/=4.and.nspden==4) then
     ! Use of C pointers to avoid copies (when ISO C bindings are available)
     ! %@1$ xlf v15 compiler requires a auxilliary cptr variable
#ifdef LIBPAW_ISO_C_BINDING
     cptr=c_loc(vxc_updn(1,1,1))
     call c_f_pointer(cptr,vxc_diag,shape=[nrad*npts,nspden_updn])
     cptr=c_loc(vxc(1,1,1))
     call c_f_pointer(cptr,vxc_nc,shape=[nrad*npts,nspden])
     cptr=c_loc(mag(1,1,1))
     call c_f_pointer(cptr,mag_,shape=[nrad*npts,3])
#else
     LIBPAW_ALLOCATE(vxc_diag,(nrad*npts,nspden_updn))
     LIBPAW_ALLOCATE(vxc_nc,(nrad*npts,nspden))
     LIBPAW_ALLOCATE(mag_,(nrad*npts,3))
     vxc_diag=reshape(vxc_updn,[nrad*npts,nspden_updn])
     mag_=reshape(mag,[nrad*npts,3])
#endif
     call pawxc_rotate_back_mag(vxc_diag,vxc_nc,mag_,nrad*npts)
#ifndef LIBPAW_ISO_C_BINDING
     vxc=reshape(vxc_nc,[nrad,npts,nspden])
     LIBPAW_DEALLOCATE(vxc_diag)
     LIBPAW_DEALLOCATE(mag_)
     LIBPAW_DEALLOCATE(vxc_nc)
#endif
     LIBPAW_POINTER_DEALLOCATE(vxc_updn)
     LIBPAW_DEALLOCATE(mag)
   end if

!  ----------------------------------------------------------------------
!  ----- Accumulate and store XC double-counting energy
!  ----------------------------------------------------------------------
   if (option==0.or.option==2) then
     LIBPAW_ALLOCATE(ff,(nrad))
     do ipts=1,npts !  Do loop on the angular part
!      Compute density for this (theta,phi)
       rhoarr(:,:)=zero
       if (usexcnhat==0) rho_=>rhor
       if (usexcnhat/=0) rho_=>rhohat
       do ispden=1,nspden
         do ilm=1,lm_size_eff
           if (lmselect(ilm)) then
             rhoarr(1:nrad,ispden)=rhoarr(1:nrad,ispden)+rho_(1:nrad,ilm,ispden)*pawang%ylmr(ilm,ipts)
           end if
         end do
       end do
!      Compute integral of Vxc*rho
       if (nspden/=4) then
         ff(:)=vxc(:,ipts,1)*rhoarr(:,nspden)
         if (nspden==2) ff(:)=ff(:)+vxc(:,ipts,2)*(rhoarr(:,1)-rhoarr(:,2))
       else
         ff(:)=half*(vxc(:,ipts,1)*(rhoarr(:,1)+rhoarr(:,4)) &
                    +vxc(:,ipts,2)*(rhoarr(:,1)-rhoarr(:,4))) &
&                   +vxc(:,ipts,3)*rhoarr(:,2)-vxc(:,ipts,4)*rhoarr(:,3)
       end if
       ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
       call simp_gen(vxcrho,ff,pawrad)
       enxcdc=enxcdc+vxcrho*pawang%angwgth(ipts)
     end do ! End of the loop on npts (angular part)
     LIBPAW_DEALLOCATE(ff)
   end if ! option

!  ----------------------------------------------------------------------
!  ----- End
!  ----------------------------------------------------------------------
!  Add the four*pi factor of the Exc and Excdc angular integration
   if (option/=1.and.option/=5) enxc=enxc*four_pi
   if (option==0.or.option==2) enxcdc=enxcdc*four_pi

!  Final memory deallocation
   nullify(rho_)
   LIBPAW_DEALLOCATE(rhoarr)
   if (usexcnhat>0)  then
     LIBPAW_DEALLOCATE(rhohat)
   end if
   if (xclevel==2) then
     LIBPAW_DEALLOCATE(gxc)
     LIBPAW_DEALLOCATE(dylmdr)
     if (mgga==1) then
       LIBPAW_DEALLOCATE(tauarr)
       if (use_laplacian==1) then
         LIBPAW_DEALLOCATE(d2ylmdr)
       end if
     end if
   end if

!  ------------------------------------
!  End IF a xc part has to be computed
 end if

end subroutine pawxc
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxcpositron
!! NAME
!! pawxcpositron
!!
!! FUNCTION
!! Compute electron-positron correlation potential and energies inside a PAW sphere
!! LDA ONLY - USE THE DENSITY OVER A WHOLE SPHERICAL GRID (r,theta,phi)
!! Driver of XC functionals.
!!
!! INPUTS
!!  calctype=type of electronpositron calculation:
!!           calctype=1 : positron in electronic density
!!           calctype=2 : electrons in positronic density
!!  corexc(nrad)=electron core density on radial grid
!!  ixcpositron=choice of electron-positron XC scheme
!!  lm_size=size of density array rhor (see below)
!!  lmselect   (lm_size)=select the non-zero LM-moments of input density rhor    (see below)
!!  lmselect_ep(lm_size)=select the non-zero LM-moments of input density rhor_ep (see below)
!!  nhat   (nrad,lm_size,nspden)=compensation density corresponding to rhor
!!  nhat_ep(nrad,lm_size,nspden)=compensation density corresponding to rhor_ep
!!  nrad=size of radial mesh for densities/potentials (might be different from pawrad%mesh_size)
!!  nspden=number of spin-density components
!!  option=0 compute both XC energies (direct+double-counting) and potential
!!         1 compute only XC potential
!!         2 compute only XC energies (direct+double-counting)
!!         3 compute only XC energy by direct scheme
!!         4 compute only XC energy by direct scheme for spherical part of the density
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  posdensity0_limit=True if we are in the zero positron density limit
!!  rhor(nrad,lm_size,nspden)=electron (or positron) density in real space
!!                             (total in 1st half and spin-up in 2nd half if nspden=2)
!!                             Contents depends on calctype value:
!!                             calctype=1: rhor is the positronic density
!!                             calctype=2: rhor is the electronic density
!!  rhor_ep(nrad,lm_size,nspden)=electron (or positron) density in real space
!!                             (total in 1st half and spin-up in 2nd half if nspden=2)
!!                             Contents depends on calctype value:
!!                             calctype=1: rhor_ep is the electronic density
!!                             calctype=2: rhor_ep is the positronic density
!!  usecore= 1 if core density has to be used in Exc/Vxc for the electronic density ; 0 otherwise
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in double counting energy term only
!!             2 if compensation density (nhat) has to be used in Exc/Vxc and double counting energy term
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!
!! OUTPUT
!!  == if option==0, 2, 3, or 4 ==
!!    enxc=returned exchange and correlation energy (hartree)
!!  == if option==0 or 2 ==
!!    enxcdc=returned exchange-cor. contribution to double-counting energy
!!  == if option==0 or 1 ==
!!    vxc(nrad,pawang%angl_size,nspden)=xc potential
!!       (spin up in 1st half and spin-down in 2nd half if nspden=2)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!
!! PARENTS
!!      pawdenpot
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxcpositron(calctype,corexc,enxc,enxcdc,ixcpositron,lm_size,lmselect,lmselect_ep,&
&                        nhat,nhat_ep,nrad,nspden,option,pawang,pawrad,posdensity0_limit,&
&                        rhor,rhor_ep,usecore,usexcnhat,vxc,xc_denpos)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: calctype,ixcpositron,lm_size,nrad,nspden,option,usecore,usexcnhat
 logical,intent(in) :: posdensity0_limit
 real(dp),intent(in) :: xc_denpos
 real(dp),intent(out) :: enxc,enxcdc
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size),lmselect_ep(lm_size)
 real(dp),intent(in) :: corexc(nrad)
 real(dp),intent(in) :: nhat(nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: nhat_ep(nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: rhor(nrad,lm_size,nspden)
 real(dp),intent(in) :: rhor_ep(nrad,lm_size,nspden)
 real(dp),intent(out) :: vxc(nrad,pawang%angl_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ilm,ipts,iwarn,iwarnp,ngr,ngrad,npts,order
 real(dp) :: enxcr,vxcrho
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: ff(:),fxci(:),grho2(:),rhoarr(:),rhoarr_ep(:),rhoarrdc(:),vxci(:),vxci_ep(:),vxcigr(:)

! *************************************************************************

!----- Check options
 if(ixcpositron==3.or.ixcpositron==31) then
   msg='GGA is not implemented (use pawxcdev/=0)!'
   MSG_ERROR(msg)
 end if
 if(calctype/=1.and.calctype/=2) then
   msg='Invalid value for calctype!'
   MSG_BUG(msg)
 end if
 if(pawang%angl_size==0) then
   msg='pawang%angl_size=0!'
   MSG_BUG(msg)
 end if
 if(.not.allocated(pawang%ylmr)) then
   msg='pawang%ylmr must be allocated!'
   MSG_BUG(msg)
 end if
 if (option/=1) then
   if (nrad<pawrad%int_meshsz) then
     msg='When option=0,2,3,4, nrad must be greater than pawrad%int_meshsz!'
     MSG_BUG(msg)
   end if
 end if

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------

!Initialization and constants
 iwarn=0;iwarnp=1
 npts=pawang%angl_size
 order=1;ngr=0;ngrad=1 ! only LDA here !

!Initializations of output arrays
 if (option/=1) enxc=zero
 if (option==0.or.option==2) enxcdc=zero
 if (option<3) vxc(:,:,:)=zero

 if (ixcpositron==0) then ! No xc at all is applied (usually for testing)
   msg = 'Note that no xc is applied (ixcpositron=0). Returning'
   MSG_WARNING(msg)
   return
 end if

!Allocations
 LIBPAW_ALLOCATE(fxci,(nrad))
 LIBPAW_ALLOCATE(vxci,(nrad))
 LIBPAW_ALLOCATE(rhoarr,(nrad))
 LIBPAW_ALLOCATE(rhoarr_ep,(nrad))
 if (option==0.or.option==2)  then
   LIBPAW_ALLOCATE(rhoarrdc,(nrad))
 end if

!----------------------------------------------------------------------
!----- Loop on the angular part
 do ipts=1,npts

!  ----------------------------------------------------------------------
!  ----- Build several densities
!  ----------------------------------------------------------------------

!  Eventually add compensation density to input density
   rhoarr=zero;rhoarr_ep=zero
   if (usexcnhat==2) then
     do ilm=1,lm_size
       if (lmselect(ilm)) &
&       rhoarr(:)=rhoarr(:)+(rhor(:,ilm,1)+nhat(:,ilm,1))*pawang%ylmr(ilm,ipts)
     end do
     do ilm=1,lm_size
       if (lmselect_ep(ilm)) &
&       rhoarr_ep(:)=rhoarr_ep(:)+(rhor_ep(:,ilm,1)+nhat_ep(:,ilm,1))*pawang%ylmr(ilm,ipts)
     end do
   else
     do ilm=1,lm_size
       if (lmselect(ilm)) rhoarr(:)=rhoarr(:)+rhor(:,ilm,1)*pawang%ylmr(ilm,ipts)
     end do
     do ilm=1,lm_size
       if (lmselect_ep(ilm)) rhoarr_ep(:)=rhoarr_ep(:)+rhor_ep(:,ilm,1)*pawang%ylmr(ilm,ipts)
     end do
   end if

!  Store density for use in double-counting term
   if (option==0.or.option==2) rhoarrdc(:)=rhoarr(:)

!  Eventually add core density
   if (usecore==1) then
     if (calctype==1) rhoarr_ep(:)=rhoarr_ep(:)+corexc(:)
     if (calctype==2) rhoarr   (:)=rhoarr   (:)+corexc(:)
   end if

!  Make the densities positive
   if (calctype==1) then
     if (.not.posdensity0_limit) then
       call pawxc_mkdenpos_wrapper(iwarnp,nrad,1,1,rhoarr,xc_denpos)
     end if
     call pawxc_mkdenpos_wrapper(iwarn ,nrad,1,1,rhoarr_ep,xc_denpos)
   else if (calctype==2) then
     call pawxc_mkdenpos_wrapper(iwarn ,nrad,1,1,rhoarr,xc_denpos)
     if (.not.posdensity0_limit) then
       call pawxc_mkdenpos_wrapper(iwarnp,nrad,1,1,rhoarr_ep,xc_denpos)
     end if
   end if

!  ----------------------------------------------------------------------
!  ----- Compute XC data
!  ----------------------------------------------------------------------

!  electron-positron correlation for the positron
   LIBPAW_ALLOCATE(vxci_ep,(nrad))
   LIBPAW_ALLOCATE(vxcigr,(ngr))
   LIBPAW_ALLOCATE(grho2,(ngr))
   if (calctype==1) then
     call pawxc_xcpositron_wrapper(fxci,grho2,ixcpositron,ngr,nrad,posdensity0_limit,rhoarr_ep,rhoarr,vxci_ep,vxcigr,vxci)
   else if (calctype==2) then
     call pawxc_xcpositron_wrapper(fxci,grho2,ixcpositron,ngr,nrad,posdensity0_limit,rhoarr,rhoarr_ep,vxci,vxcigr,vxci_ep)
   end if
   LIBPAW_DEALLOCATE(vxci_ep)
   LIBPAW_DEALLOCATE(vxcigr)
   LIBPAW_DEALLOCATE(grho2)

!  ----------------------------------------------------------------------
!  ----- Accumulate and store XC potential
!  ----------------------------------------------------------------------
   if (option<3) then
     vxc(:,ipts,1)=vxci(:)
     if (nspden>=2) vxc(:,ipts,2)=vxci(:)
     if (nspden==4) vxc(:,ipts,3:4)=zero
   end if

!  ----------------------------------------------------------------------
!  ----- Accumulate and store XC energies
!  ----------------------------------------------------------------------

!  ----- Calculate Exc term
   if (option/=1) then
     LIBPAW_ALLOCATE(ff,(nrad))
     ff(1:nrad)=fxci(1:nrad)*pawrad%rad(1:nrad)**2
     call simp_gen(enxcr,ff,pawrad)
     LIBPAW_DEALLOCATE(ff)
     if (option/=4) enxc=enxc+enxcr*pawang%angwgth(ipts)
     if (option==4) enxc=enxc+enxcr
   end if

!  ----- Calculate Excdc double counting term
   if (option==0.or.option==2) then
     if (usexcnhat==1) then
       do ilm=1,lm_size
         if (lmselect(ilm)) then
           rhoarrdc(:)=rhoarrdc(:)+nhat(:,ilm,1)*pawang%ylmr(ilm,ipts)
         end if
       end do
     end if
     LIBPAW_ALLOCATE(ff,(nrad))
     ff(1:nrad)=vxci(1:nrad)*rhoarrdc(1:nrad)*pawrad%rad(1:nrad)**2
     call simp_gen(vxcrho,ff,pawrad)
     LIBPAW_DEALLOCATE(ff)
     enxcdc=enxcdc+vxcrho*pawang%angwgth(ipts)
   end if

!  ---------------------------------------------------
!  ----- End of the loop on npts (angular part)
 end do

!Add the four*pi factor of the angular integration
 if (option/=1) enxc=enxc*four_pi
 if (option==0.or.option==2) enxcdc=enxcdc*four_pi

!Deallocations
 LIBPAW_DEALLOCATE(fxci)
 LIBPAW_DEALLOCATE(vxci)
 LIBPAW_DEALLOCATE(rhoarr)
 LIBPAW_DEALLOCATE(rhoarr_ep)
 if (option==0.or.option==2)  then
   LIBPAW_DEALLOCATE(rhoarrdc)
 end if

end subroutine pawxcpositron
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxc_dfpt
!! NAME
!! pawxc_dfpt
!!
!! FUNCTION
!! Compute first-order change of XC potential and contribution to
!! 2nd-order change of XC energy inside a PAW sphere.
!! LDA+GGA - USE THE DENSITY OVER A WHOLE SPHERICAL GRID (r,theta,phi)
!!
!! INPUTS
!!  corexc1(cplex_den*nrad)=first-order change of core density on radial grid
!!  cplex_den= if 1, 1st-order densities are REAL, if 2, COMPLEX
!!  cplex_vxc= if 1, 1st-order XC potential is complex, if 2, COMPLEX
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nrad,pawang%angl_size,nkxc)=GS xc kernel
!!  lm_size=size of density array rhor (see below)
!!  lmselect(lm_size)=select the non-zero LM-moments of input density rhor1
!!  nhat1(cplex_den*nrad,lm_size,nspden)=first-order change of compensation density
!!                                        (total in 1st half and spin-up in 2nd half if nspden=2)
!!  nkxc=second dimension of the kxc array
!!  non_magnetic_xc= if true, handle density/potential as non-magnetic (even if it is)
!!  nrad=size of radial mesh for densities/potentials (might be different from pawrad%mesh_size)
!!  nspden=number of spin-density components
!!  option=0  compute both 2nd-order XC energy and 1st-order potential
!!         1  compute only 1st-order XC potential
!!         2  compute only 2nd-order XC energy, XC potential is temporary computed here
!!         3  compute only 2nd-order XC energy, XC potential is input in vxc1(:)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  rhor1(cplex_den*nrad,lm_size,nspden)=first-order change of density
!!  usecore= 1 if core density has to be used in Exc/Vxc ; 0 otherwise
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in d2Exc only
!!             2 if compensation density (nhat) has to be used in d2Exc and Vxc1
!!  vxc(nrad,pawang%angl_size,nspden)=GS xc potential
!!  xclevel= XC functional level
!!
!! OUTPUT
!!  == if option=0 or 2 or 3 ==
!!    d2enxc   =returned exchange-cor. contribution to 2nd-order XC energy
!!    d2enxc_im=returned IMAGINARY PART of exchange-cor. contribution to 2nd-order XC energy
!!              (optional argument)
!!
!! SIDE EFFECTS
!!    vxc1(cplex_vxc*nrad,pawang%angl_size,nspden)=1st-order XC potential
!!      Output if option==0 or 1
!!      Unused if option==2
!!      Input  if option==3
!!
!! NOTES
!!  Content of Kxc array:
!!   ===== if LDA
!!    if nspden==1: kxc(:,1)= d2Exc/drho2
!!                 (kxc(:,2)= d2Exc/drho_up drho_dn)
!!    if nspden>=2: kxc(:,1)= d2Exc/drho_up drho_up
!!                  kxc(:,2)= d2Exc/drho_up drho_dn
!!                  kxc(:,3)= d2Exc/drho_dn drho_dn
!!    if nspden==4: kxc(:,4:6)= (m_x, m_y, m_z) (magnetization)
!!   ===== if GGA
!!    if nspden==1:
!!       kxc(:,1)= d2Exc/drho2
!!       kxc(:,2)= 1/|grad(rho)| dExc/d|grad(rho)|
!!       kxc(:,3)= 1/|grad(rho)| d2Exc/d|grad(rho)| drho
!!       kxc(:,4)= 1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dExc/d|grad(rho)| )
!!       kxc(:,5)= gradx(rho)
!!       kxc(:,6)= grady(rho)
!!       kxc(:,7)= gradz(rho)
!!    if nspden>=2:
!!       kxc(:,1)= d2Exc/drho_up drho_up
!!       kxc(:,2)= d2Exc/drho_up drho_dn
!!       kxc(:,3)= d2Exc/drho_dn drho_dn
!!       kxc(:,4)= 1/|grad(rho_up)| dEx/d|grad(rho_up)|
!!       kxc(:,5)= 1/|grad(rho_dn)| dEx/d|grad(rho_dn)|
!!       kxc(:,6)= 1/|grad(rho_up)| d2Ex/d|grad(rho_up)| drho_up
!!       kxc(:,7)= 1/|grad(rho_dn)| d2Ex/d|grad(rho_dn)| drho_dn
!!       kxc(:,8)= 1/|grad(rho_up)| * d/d|grad(rho_up)| ( 1/|grad(rho_up)| dEx/d|grad(rho_up)| )
!!       kxc(:,9)= 1/|grad(rho_dn)| * d/d|grad(rho_dn)| ( 1/|grad(rho_dn)| dEx/d|grad(rho_dn)| )
!!       kxc(:,10)=1/|grad(rho)| dEc/d|grad(rho)|
!!       kxc(:,11)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_up
!!       kxc(:,12)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_dn
!!       kxc(:,13)=1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dEc/d|grad(rho)| )
!!       kxc(:,14)=gradx(rho_up)
!!       kxc(:,15)=gradx(rho_dn)
!!       kxc(:,16)=grady(rho_up)
!!       kxc(:,17)=grady(rho_dn)
!!       kxc(:,18)=gradz(rho_up)
!!       kxc(:,19)=gradz(rho_dn)
!!    if nspden==4:
!!       kxc(:,20:22)= (m_x, m_y, m_z) (magnetization)
!!
!! PARENTS
!!      pawdenpot,pawdfptenergy
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxc_dfpt(corexc1,cplex_den,cplex_vxc,d2enxc,ixc,kxc,lm_size,lmselect,nhat1,&
&                 nkxc,non_magnetic_xc,nrad,nspden,option,pawang,pawrad,rhor1,&
&                 usecore,usexcnhat,vxc,vxc1,xclevel,&
&                 d2enxc_im) ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex_den,cplex_vxc,ixc,lm_size,nkxc,nrad,nspden,option
 integer,intent(in) :: usecore,usexcnhat,xclevel
 logical,intent(in) :: non_magnetic_xc
 real(dp),intent(out) :: d2enxc
 real(dp),intent(out),optional :: d2enxc_im
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size)
 real(dp),intent(in) :: corexc1(cplex_den*nrad)
 real(dp),intent(in) :: nhat1(cplex_den*nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in),target :: kxc(nrad,pawang%angl_size,nkxc)
 real(dp),intent(in),target :: vxc(nrad,pawang%angl_size,nspden)
 real(dp),intent(in),target :: rhor1(cplex_den*nrad,lm_size,nspden)
 real(dp),intent(inout),target :: vxc1(cplex_vxc*nrad,pawang%angl_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ilm,ipts,ir,ispden,jr,kr,lm_size_eff,nkxc_cur,npts,nspden_updn
 logical :: need_impart
 real(dp),parameter :: tol24=tol12*tol12
 real(dp) :: coeff_grho,coeff_grho_corr,coeff_grho_dn,coeff_grho_up
 real(dp) :: coeff_grhoim,coeff_grhoim_corr,coeff_grhoim_dn,coeff_grhoim_up
 real(dp) :: dylmdr_ii,factor,factor_ang_intg,ylm_ii
 real(dp) :: grho_grho1,grho_grho1_up,grho_grho1_dn
 real(dp) :: grho_grho1im,grho_grho1im_up,grho_grho1im_dn
 real(dp) :: rho1_dn,rho1_up,rho1im_dn,rho1im_up
 real(dp) :: ro11i,ro11r,ro12i,ro12r,ro21i,ro21r,ro22i,ro22r
 real(dp) :: v11i,v11r,v12i,v12r,v21i,v21r,v22i,v22r,vxcrho
 character(len=500) :: msg
!arrays
 real(dp) :: g0(3),g0_dn(3),g0_up(3),g1(3),g1_dn(3),g1_up(3)
 real(dp) :: g1im(3),g1im_dn(3),g1im_up(3)
 real(dp) :: gxc1i(3,2),gxc1r(3,2)
 real(dp),allocatable :: dgxc1(:),drho1(:,:),drho1core(:,:),dylmdr(:,:,:)
 real(dp),allocatable :: ff(:),gg(:),grho1arr(:,:,:),gxc1(:,:,:,:)
 real(dp),allocatable,target :: rhohat1(:,:,:),rho1arr(:,:)
 real(dp), LIBPAW_CONTIGUOUS pointer :: kxc_(:,:),mag(:,:)
 real(dp), LIBPAW_CONTIGUOUS pointer :: rho1_(:,:,:),rho1_nc(:,:),rho1_updn(:,:)
 real(dp), LIBPAW_CONTIGUOUS pointer :: vxc_(:,:),vxc1_(:,:,:),vxc1_diag(:,:)
 real(dp), LIBPAW_CONTIGUOUS pointer :: vxc1_nc(:,:),vxc1_updn(:,:,:)
#ifdef LIBPAW_ISO_C_BINDING
 type(C_PTR) :: cptr
#endif

! *************************************************************************

!----------------------------------------------------------------------
!----- Check options
!----------------------------------------------------------------------

 if(option<0.or.option>3) then
   msg='wrong option!'
   MSG_BUG(msg)
 end if
 if(option/=3) then
   call pawxc_get_nkxc(nkxc_cur,nspden,xclevel)
   if (nkxc/=nkxc_cur) then
     msg='Wrong dimension for array kxc!'
     MSG_BUG(msg)
   end if
   if(xclevel==2.and.nspden==4) then
     msg='PAW non-collinear magnetism not compatible with GGA!'
     MSG_ERROR(msg)
   end if
 end if
 if(pawang%angl_size==0) then
   msg='pawang%angl_size=0!'
   MSG_BUG(msg)
 end if
 if(.not.allocated(pawang%ylmr)) then
   msg='pawang%ylmr must be allocated!'
   MSG_BUG(msg)
 end if
 if(xclevel==2.and.(.not.allocated(pawang%ylmrgr))) then
   msg='pawang%ylmrgr must be allocated!'
   MSG_BUG(msg)
 end if
 if (option/=1) then
   if (nrad<pawrad%int_meshsz) then
     msg='When option=0,2, nrad must be greater than pawrad%int_meshsz!'
     MSG_BUG(msg)
   end if
 end if

!----------------------------------------------------------------------
!----- Initializations / allocations
!----------------------------------------------------------------------

 npts=pawang%angl_size
 lm_size_eff=min(lm_size,pawang%ylm_size)
 nspden_updn=min(nspden,2)

 need_impart=present(d2enxc_im)
 if (option/=1) then
   d2enxc=zero
   if (need_impart) d2enxc_im=zero
 end if
 if (option<=1) vxc1(:,:,:)=zero

!Special case: no XC applied
 if (ixc==0.or.(nkxc==0.and.option/=3)) then
   msg='Note that no xc is applied (ixc=0). Returning'
   MSG_WARNING(msg)
   return
 end if

 LIBPAW_ALLOCATE(rho1arr,(cplex_den*nrad,nspden))
 if (usexcnhat>0) then
   LIBPAW_ALLOCATE(rhohat1,(cplex_den*nrad,lm_size,nspden))
   rhohat1(:,:,:)=rhor1(:,:,:)+nhat1(:,:,:)
 end if

 if (option==2) then
   LIBPAW_POINTER_ALLOCATE(vxc1_,(cplex_vxc*nrad,npts,nspden))
 else
   vxc1_ => vxc1
 end if

!Need gradients and additional allocations in case of GGA
 if (xclevel==2.and.option/=3) then
   LIBPAW_ALLOCATE(gxc1,(cplex_vxc*nrad,3,pawang%ylm_size,nspden))
   gxc1=zero
   if (usecore==1) then
     LIBPAW_ALLOCATE(drho1core,(nrad,cplex_den))
     if (cplex_den==1)  then
       call nderiv_gen(drho1core(:,1),corexc1,pawrad)
     else
       LIBPAW_ALLOCATE(ff,(nrad))
       LIBPAW_ALLOCATE(gg,(nrad))
       do ir=1,nrad
         ff(ir)=corexc1(2*ir-1)
         gg(ir)=corexc1(2*ir  )
       end do
       call nderiv_gen(drho1core(:,1),ff,pawrad)
       call nderiv_gen(drho1core(:,2),gg,pawrad)
       LIBPAW_DEALLOCATE(ff)
       LIBPAW_DEALLOCATE(gg)
     end if
   end if
!  Convert Ylm derivatives from normalized to standard cartesian coordinates
!  dYlm/dr_i = { dYlm/dr_i^hat - Sum_j[ dYlm/dr_j^hat (r_j/r)] } * (1/r)
   LIBPAW_ALLOCATE(dylmdr,(3,npts,pawang%ylm_size))
   do ilm=1,pawang%ylm_size
     do ipts=1,npts
       factor=sum(pawang%ylmrgr(1:3,ilm,ipts)*pawang%anginit(1:3,ipts))
       dylmdr(1:3,ipts,ilm)=pawang%ylmrgr(1:3,ilm,ipts)-factor*pawang%anginit(1:3,ipts)
     end do
   end do
 end if

!----------------------------------------------------------------------
!----- Accumulate and store 1st-order change of XC potential
!----------------------------------------------------------------------

 if (option/=3) then

   if (nspden/=4) then
     rho1_updn => rho1arr
     vxc1_updn => vxc1_
   else
     LIBPAW_POINTER_ALLOCATE(rho1_updn,(cplex_den*nrad,nspden_updn))
     LIBPAW_POINTER_ALLOCATE(vxc1_updn,(cplex_vxc*nrad,npts,nspden_updn))
     LIBPAW_POINTER_ALLOCATE(rho1_nc,(cplex_den*nrad*npts,nspden))
     LIBPAW_POINTER_ALLOCATE(mag,(nrad,3))
   end if

!  Do loop on the angular part (theta,phi)
   do ipts=1,npts

!    Copy the input 1st-order density for this (theta,phi)
     rho1arr(:,:)=zero
     if (usexcnhat< 2) rho1_=>rhor1
     if (usexcnhat==2) rho1_=>rhohat1
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect(ilm)) rho1arr(:,ispden)=rho1arr(:,ispden) &
&       +rho1_(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
     end do
     if (usecore==1) then
       rho1arr(:,1)=rho1arr(:,1)+corexc1(:)
       if (nspden==2) rho1arr(:,2)=rho1arr(:,2)+half*corexc1(:)
     end if

!    Optionally suppress magnetic part
     if(non_magnetic_xc) then
       if(nspden==2) rho1arr(:,2)=rho1arr(:,1)*half
       if(nspden==4) rho1arr(:,2:4)=zero
     endif

!    Non-collinear magnetism: rotate magnetization and get a collinear density
     if (nspden==4) then
       !Store non rotated rho^(1) for future use
       ii=(ipts-1)*cplex_den*nrad
       do ispden=1,nspden
         rho1_nc(ii+1:ii+cplex_den*nrad,ispden)=rho1arr(1:cplex_den*nrad,ispden)
       end do
       !Extract magnetization from kxc
       do ii=1,3
         mag(1:nrad,ii)=kxc(:,ipts,ii)
       end do
       !Rotate rhoarr1 -> rhoarr1_
       !Should use cplex_den
       call pawxc_rotate_mag(rho1arr,rho1_updn,mag,nrad,rho_out_format=2)
     end if

!    =======================================================================
!    ======================= LDA ===========================================
!    =======================================================================
     if (xclevel==1.or.ixc==13) then

!      Non-spin-polarized
       if (nspden_updn==1) then
         if (cplex_vxc==1) then
           if (cplex_den==1) then  ! cplex_vxc==1 and cplex_den==1
             vxc1_updn(1:nrad,ipts,1)=kxc(1:nrad,ipts,1)*rho1_updn(1:nrad,1)
           else                    ! cplex_vxc==1 and cplex_den==2
             do ir=1,nrad
               vxc1_updn(ir,ipts,1)=kxc(ir,ipts,1)*rho1_updn(2*ir-1,1)
             end do
           end if
         else
           if (cplex_den==1) then  ! cplex_vxc==2 and cplex_den==1
             do ir=1,nrad
               vxc1_updn(2*ir-1,ipts,1)=kxc(ir,ipts,1)*rho1_updn(ir,1)
               vxc1_updn(2*ir  ,ipts,1)=zero
             end do
           else                    ! cplex_vxc==2 and cplex_den==2
             do ir=1,nrad
               vxc1_updn(2*ir-1,ipts,1)=kxc(ir,ipts,1)*rho1_updn(2*ir-1,1)
               vxc1_updn(2*ir  ,ipts,1)=kxc(ir,ipts,1)*rho1_updn(2*ir  ,1)
             end do
           end if
         end if

!        Spin-polarized
       else
         if (cplex_vxc==1) then
           if (cplex_den==1) then  ! cplex_vxc==1 and cplex_den==1
             do ir=1,nrad
               rho1_up=rho1_updn(ir,2);rho1_dn=rho1_updn(ir,1)-rho1_up
               vxc1_updn(ir,ipts,1)=kxc(ir,ipts,1)*rho1_up+kxc(ir,ipts,2)*rho1_dn
               vxc1_updn(ir,ipts,2)=kxc(ir,ipts,2)*rho1_up+kxc(ir,ipts,3)*rho1_dn
             end do
           else                    ! cplex_vxc==1 and cplex_den==2
             do ir=1,nrad
               jr=2*ir-1
               rho1_up=rho1_updn(jr,2);rho1_dn=rho1_updn(jr,1)-rho1_up
               vxc1_updn(ir,ipts,1)=kxc(ir,ipts,1)*rho1_up+kxc(ir,ipts,2)*rho1_dn
               vxc1_updn(ir,ipts,2)=kxc(ir,ipts,2)*rho1_up+kxc(ir,ipts,3)*rho1_dn
             end do
           end if
         else
           if (cplex_den==1) then  ! cplex_vxc==2 and cplex_den==1
             do ir=1,nrad
               jr=2*ir-1
               rho1_up=rho1_updn(ir,2);rho1_dn=rho1_updn(ir,1)-rho1_up
               vxc1_updn(jr,ipts,1)=kxc(ir,ipts,1)*rho1_up+kxc(ir,ipts,2)*rho1_dn
               vxc1_updn(jr,ipts,2)=kxc(ir,ipts,2)*rho1_up+kxc(ir,ipts,3)*rho1_dn
             end do
           else                    ! cplex_vxc==2 and cplex_den==2
             do ir=1,nrad
               jr=2*ir
               rho1_up  =rho1_updn(jr-1,2);rho1_dn  =rho1_updn(jr-1,1)-rho1_up
               rho1im_up=rho1_updn(jr  ,2);rho1im_dn=rho1_updn(jr  ,1)-rho1im_up
               vxc1_updn(jr-1,ipts,1)=kxc(ir,ipts,1)*rho1_up  +kxc(ir,ipts,2)*rho1_dn
               vxc1_updn(jr  ,ipts,1)=kxc(ir,ipts,1)*rho1im_up+kxc(ir,ipts,2)*rho1im_dn
               vxc1_updn(jr-1,ipts,2)=kxc(ir,ipts,2)*rho1_up  +kxc(ir,ipts,3)*rho1_dn
               vxc1_updn(jr  ,ipts,2)=kxc(ir,ipts,2)*rho1im_up+kxc(ir,ipts,3)*rho1im_dn
             end do
           end if
         end if
       end if

     else
!      =======================================================================
!      ======================= GGA ===========================================
!      =======================================================================

!      Compute the gradient of the first-order density
       LIBPAW_ALLOCATE(drho1,(nrad,cplex_den))
       LIBPAW_ALLOCATE(grho1arr,(cplex_den*nrad,nspden,3))
       grho1arr(:,:,1:3)=zero
       if (cplex_den==1) then
         LIBPAW_ALLOCATE(ff,(nrad))
         do ispden=1,nspden_updn
           do ilm=1,lm_size_eff
             if (lmselect(ilm)) then
               ff(1:nrad)=rho1_(1:nrad,ilm,ispden)
               call nderiv_gen(drho1(:,1),ff,pawrad)
               ff(2:nrad)=ff(2:nrad)/pawrad%rad(2:nrad)
               call pawrad_deducer0(ff,nrad,pawrad)
               do ii=1,3
                 ylm_ii=pawang%ylmr(ilm,ipts)*pawang%anginit(ii,ipts)
                 dylmdr_ii=dylmdr(ii,ipts,ilm)
                 grho1arr(1:nrad,ispden,ii)=grho1arr(1:nrad,ispden,ii) &
&                 +drho1(1:nrad,1)*ylm_ii+ff(1:nrad)*dylmdr_ii
               end do
             end if
           end do
         end do
         LIBPAW_DEALLOCATE(ff)
       else
         LIBPAW_ALLOCATE(ff,(nrad))
         LIBPAW_ALLOCATE(gg,(nrad))
         do ispden=1,nspden_updn
           do ilm=1,lm_size_eff
             if (lmselect(ilm)) then
               do ir=1,nrad
                 ff(ir)=rho1_(2*ir-1,ilm,ispden)
                 gg(ir)=rho1_(2*ir  ,ilm,ispden)
               end do
               call nderiv_gen(drho1(:,1),ff,pawrad)
               call nderiv_gen(drho1(:,2),gg,pawrad)
               ff(2:nrad)=ff(2:nrad)/pawrad%rad(2:nrad)
               gg(2:nrad)=gg(2:nrad)/pawrad%rad(2:nrad)
               call pawrad_deducer0(ff,nrad,pawrad)
               call pawrad_deducer0(gg,nrad,pawrad)
               do ii=1,3
                 ylm_ii=pawang%ylmr(ilm,ipts)*pawang%anginit(ii,ipts)
                 dylmdr_ii=dylmdr(ii,ipts,ilm)
                 do ir=2,nrad
                   jr=2*ir
                   grho1arr(jr-1,ispden,ii)=grho1arr(jr-1,ispden,ii) &
&                   +drho1(ir,1)*ylm_ii+ff(ir)*dylmdr_ii
                   grho1arr(jr  ,ispden,ii)=grho1arr(jr  ,ispden,ii) &
&                   +drho1(ir,2)*ylm_ii+gg(ir)*dylmdr_ii
                 end do
               end do
             end if
           end do
         end do
         LIBPAW_DEALLOCATE(ff)
         LIBPAW_DEALLOCATE(gg)
       end if
       if (usecore==1) then
         factor=one;if (nspden_updn==2) factor=half
         if (cplex_den==1) then
           do ispden=1,nspden_updn
             do ii=1,3
               grho1arr(1:nrad,ispden,ii)=grho1arr(1:nrad,ispden,ii) &
&               +factor*drho1core(1:nrad,1)*pawang%anginit(ii,ipts)
             end do
           end do
         else
           do ispden=1,nspden_updn
             do ii=1,3
               do ir=1,nrad
                 jr=2*ir
                 grho1arr(jr-1,ispden,ii)=grho1arr(jr-1,ispden,ii) &
&                 +factor*drho1core(ir,1)*pawang%anginit(ii,ipts)
                 grho1arr(jr  ,ispden,ii)=grho1arr(jr  ,ispden,ii) &
&                 +factor*drho1core(ir,2)*pawang%anginit(ii,ipts)
               end do
             end do
           end do
         end if
       end if
       LIBPAW_DEALLOCATE(drho1)

!      Optionally suppress magnetic part
       if(non_magnetic_xc) then
         do ii=1,3
           if(nspden==2) grho1arr(:,2,ii)=grho1arr(:,1,ii)*half
           if(nspden==4) grho1arr(:,2:4,ii)=zero
         end do
       endif

!      Apply XC kernel
!      Will compute Vxc^(1) as: vxc1 - Nabla .dot. gxc1

!      Scaling factor for angular integrals: four_pi x spin_factor
       factor_ang_intg=four_pi;if (nspden_updn==1) factor_ang_intg=two_pi

!      A- NON POLARIZED SYSTEMS
       if (nspden_updn==1) then

         do ir=1,nrad
           jr=cplex_den*(ir-1)+1 ; kr=cplex_vxc*(ir-1)+1

           g0(:)=kxc(ir,ipts,5:7) ; g1(:)=grho1arr(jr,1,1:3)
           grho_grho1=dot_product(g0,g1)
           coeff_grho=kxc(ir,ipts,3)*rho1_updn(jr,1)+kxc(ir,ipts,4)*grho_grho1
           vxc1_updn(kr,ipts,1)=kxc(ir,ipts,1)*rho1_updn(jr,1)+kxc(ir,ipts,3)*grho_grho1
           gxc1r(:,1)=g1(:)*kxc(ir,ipts,2)+g0(:)*coeff_grho
           !Accumulate gxc1_lm moments as Intg[gxc1(omega).Ylm(omega).d_omega]
           do ilm=1,pawang%ylm_size
             ylm_ii=pawang%ylmr(ilm,ipts)*pawang%angwgth(ipts)*factor_ang_intg
             do ii=1,3
               gxc1(kr,ii,ilm,1)=gxc1(ir,ii,ilm,1)+gxc1r(ii,1)*ylm_ii
             end do
           end do
           if (cplex_vxc==2) then
             if (cplex_den==2) then
               g1im(:)=grho1arr(jr+1,1,1:3)
               grho_grho1im=dot_product(g0,g1im)
               coeff_grhoim=kxc(ir,ipts,3)*rho1_updn(jr+1,1)+kxc(ir,ipts,4)*grho_grho1im
               vxc1_updn(kr+1,ipts,1)=kxc(ir,ipts,1)*rho1_updn(jr+1,1)+kxc(ir,ipts,3)*grho_grho1im
               gxc1i(:,1)=g1im(:)*kxc(ir,ipts,2)+g0(:)*coeff_grhoim
               !Accumulate gxc1_lm moments as Intg[gxc1(omega).Ylm(omega).d_omega]
               do ilm=1,pawang%ylm_size
                 ylm_ii=pawang%ylmr(ilm,ipts)*pawang%angwgth(ipts)*factor_ang_intg
                 do ii=1,3
                   gxc1(kr+1,ii,ilm,1)=gxc1(kr+1,ii,ilm,1)+gxc1i(ii,1)*ylm_ii
                 end do
               end do
             else
               vxc1_updn(kr+1,ipts,1)=zero ; gxc1i(:,1)=zero
             end if
           end if
         end do ! ir

!      B- POLARIZED SYSTEMS (COLLINEAR)
       else ! nspden_updn==2

         do ir=1,nrad
           jr=cplex_den*(ir-1)+1 ; kr=cplex_vxc*(ir-1)+1

           rho1_up=rho1_updn(jr,2);rho1_dn=rho1_updn(jr,1)-rho1_up
           g0_up(1)=kxc(ir,ipts,15);g0_dn(1)=kxc(ir,ipts,14)-kxc(ir,ipts,15)
           g0_up(2)=kxc(ir,ipts,17);g0_dn(2)=kxc(ir,ipts,16)-kxc(ir,ipts,17)
           g0_up(3)=kxc(ir,ipts,19);g0_dn(3)=kxc(ir,ipts,18)-kxc(ir,ipts,19)
           g1_up(:)=grho1arr(jr,2,:);g1_dn(:)=grho1arr(jr,1,:)-grho1arr(jr,2,:)
           g0(:)=g0_up(:)+g0_dn(:);g1(:)=g1_up(:)+g1_dn(:)
           grho_grho1_up=dot_product(g0_up,g1_up)
           grho_grho1_dn=dot_product(g0_dn,g1_dn)
           grho_grho1   =dot_product(g0,g1)
           coeff_grho_corr=kxc(ir,ipts,11)*rho1_up &
&                         +kxc(ir,ipts,12)*rho1_dn &
&                         +kxc(ir,ipts,13)*grho_grho1
           coeff_grho_up=kxc(ir,ipts,6)*rho1_up &
&                       +kxc(ir,ipts,8)*grho_grho1_up
           coeff_grho_dn=kxc(ir,ipts,7)*rho1_dn &
&                       +kxc(ir,ipts,9)*grho_grho1_dn
           vxc1_updn(kr,ipts,1)=kxc(ir,ipts, 1)*rho1_up &
&                          +kxc(ir,ipts, 2)*rho1_dn &
&                          +kxc(ir,ipts, 6)*grho_grho1_up &
&                          +kxc(ir,ipts,11)*grho_grho1
           vxc1_updn(kr,ipts,2)=kxc(ir,ipts, 3)*rho1_dn &
&                          +kxc(ir,ipts, 2)*rho1_up &
&                          +kxc(ir,ipts, 7)*grho_grho1_dn &
&                          +kxc(ir,ipts,12)*grho_grho1
           gxc1r(:,1)=(kxc(ir,ipts,4)+kxc(ir,ipts,10))*g1_up(:) &
&                    +kxc(ir,ipts,10)                 *g1_dn(:) &
&                    +coeff_grho_up                   *g0_up(:) &
&                    +coeff_grho_corr                 *g0(:)
           gxc1r(:,2)=(kxc(ir,ipts,5)+kxc(ir,ipts,10))*g1_dn(:) &
&                    +kxc(ir,ipts,10)                 *g1_up(:) &
&                    +coeff_grho_dn                   *g0_dn(:) &
&                    +coeff_grho_corr                 *g0(:)
           !Accumulate gxc1_lm moments as Intg[gxc1(omega).Ylm(omega).d_omega]
           do ispden=1,nspden_updn
             do ilm=1,pawang%ylm_size
               ylm_ii=pawang%ylmr(ilm,ipts)*pawang%angwgth(ipts)*factor_ang_intg
               do ii=1,3
                 gxc1(kr,ii,ilm,ispden)=gxc1(kr,ii,ilm,ispden)+gxc1r(ii,ispden)*ylm_ii
               end do
             end do
           end do

           if (cplex_vxc==2) then
             if (cplex_den==2) then
               rho1im_up=rho1_updn(jr+1,2);rho1im_dn=rho1_updn(jr+1,1)-rho1im_up
               g1im_up(:)=grho1arr(jr+1,2,:);g1im_dn(:)=grho1arr(jr+1,1,:)-grho1arr(jr+1,2,:)
               g1im(:)=g1im_up(:)+g1im_dn(:)
               grho_grho1im_up=dot_product(g0_up,g1im_up)
               grho_grho1im_dn=dot_product(g0_dn,g1im_dn)
               grho_grho1im   =dot_product(g0,g1im)
               coeff_grhoim_corr=kxc(ir,ipts,11)*rho1im_up &
&                               +kxc(ir,ipts,12)*rho1im_dn &
&                               +kxc(ir,ipts,13)*grho_grho1im
               coeff_grhoim_up=kxc(ir,ipts,6)*rho1im_up &
&                             +kxc(ir,ipts,8)*grho_grho1im_up
               coeff_grhoim_dn=kxc(ir,ipts,7)*rho1im_dn &
&                             +kxc(ir,ipts,9)*grho_grho1im_dn
               vxc1_updn(kr+1,ipts,1)=kxc(ir,ipts, 1)*rho1im_up &
&                                +kxc(ir,ipts, 2)*rho1im_dn &
&                                +kxc(ir,ipts, 6)*grho_grho1im_up   &
&                                +kxc(ir,ipts,11)*grho_grho1im
               vxc1_updn(kr+1,ipts,2)=kxc(ir,ipts, 3)*rho1im_dn &
&                                +kxc(ir,ipts, 2)*rho1im_up &
&                                +kxc(ir,ipts, 7)*grho_grho1im_dn   &
&                                +kxc(ir,ipts,12)*grho_grho1im
               gxc1i(:,1)=(kxc(ir,ipts,4)+kxc(ir,ipts,10))*g1im_up(:) &
&                        +kxc(ir,ipts,10)                 *g1im_dn(:) &
&                        +coeff_grhoim_up                 *g0_up(:)   &
&                        +coeff_grhoim_corr               *g0(:)
               gxc1i(:,2)=(kxc(ir,ipts,5)+kxc(ir,ipts,10))*g1im_dn(:) &
&                        +kxc(ir,ipts,10)                 *g1im_up(:) &
&                        +coeff_grhoim_dn                 *g0_dn(:)   &
&                        +coeff_grhoim_corr               *g0(:)
               !Accumulate gxc1_lm moments as Intg[gxc1(omega).Ylm(omega).d_omega]
               do ispden=1,nspden_updn
                 do ilm=1,pawang%ylm_size
                   ylm_ii=pawang%ylmr(ilm,ipts)*pawang%angwgth(ipts)*factor_ang_intg
                   do ii=1,3
                     gxc1(kr+1,ii,ilm,ispden)=gxc1(kr+1,ii,ilm,ispden)+gxc1i(ii,ispden)*ylm_ii
                   end do
                 end do
               end do
             else
               vxc1_updn(kr+1,ipts,1:2)=zero ; gxc1i(:,1:2)=zero
             end if
           end if

         end do ! ir

       end if ! nspden_updn

       LIBPAW_DEALLOCATE(grho1arr)

     end if ! LDA or GGA

!  ----- End of the loop on npts (angular part)
   end do

!  Deallocate memory
   if (xclevel==2.and.usecore==1)  then
     LIBPAW_DEALLOCATE(drho1core)
   end if
   if (nspden==4) then
     LIBPAW_POINTER_DEALLOCATE(rho1_updn)
     LIBPAW_POINTER_DEALLOCATE(mag)
   end if

 end if ! option/=3

!----------------------------------------------------------------------
!----- If GGA, modify 1st-order potential with term from density gradient
!----------------------------------------------------------------------
 if (xclevel==2.and.ixc/=13.and.option/=3) then
!  Compute divergence of gxc1 and substract it from Vxc1

!  Need to multiply gxc1 by 2 in the non-polarised case
   factor=one;if (nspden_updn==1) factor=two

   LIBPAW_ALLOCATE(dgxc1,(nrad))
   LIBPAW_ALLOCATE(gg,(nrad))
   do ispden=1,nspden_updn
     do ilm=1,pawang%ylm_size
       do ii=1,3
         do ir=1,nrad
           jr=cplex_vxc*(ir-1)+1
           gg(ir)=gxc1(jr,ii,ilm,ispden)
         end do
         call nderiv_gen(dgxc1,gg,pawrad)
         gg(2:nrad)=gg(2:nrad)/pawrad%rad(2:nrad)
         call pawrad_deducer0(gg,nrad,pawrad)
         do ipts=1,npts
           ylm_ii=pawang%ylmr(ilm,ipts)*pawang%anginit(ii,ipts)
           dylmdr_ii=dylmdr(ii,ipts,ilm)
           do ir=1,nrad
             jr=cplex_vxc*(ir-1)+1
             vxc1_(jr,ipts,ispden)=vxc1_(jr,ipts,ispden) &
&               -factor*(dgxc1(ir)*ylm_ii+gg(ir)*dylmdr_ii)
           end do
         end do ! ipts
       end do ! ii
     end do ! ilm
   end do ! ispden
   if (cplex_vxc==2) then
     do ispden=1,nspden_updn
       do ilm=1,pawang%ylm_size
         do ii=1,3
           do ir=1,nrad
             gg(ir)=gxc1(2*ir,ii,ilm,ispden)
           end do
           call nderiv_gen(dgxc1,gg,pawrad)
           gg(2:nrad)=gg(2:nrad)/pawrad%rad(2:nrad)
           call pawrad_deducer0(gg,nrad,pawrad)
           do ipts=1,npts
             ylm_ii=pawang%ylmr(ilm,ipts)*pawang%anginit(ii,ipts)
             dylmdr_ii=dylmdr(ii,ipts,ilm)
             do ir=1,nrad
               vxc1_(2*ir,ipts,ispden)=vxc1_(2*ir,ipts,ispden) &
  &               -factor*(dgxc1(ir)*ylm_ii+gg(ir)*dylmdr_ii)
             end do
           end do ! ipts
         end do ! ii
       end do ! ilm
     end do ! ispden
   end if ! cplex_vxc
   LIBPAW_DEALLOCATE(dgxc1)
   LIBPAW_DEALLOCATE(gg)

 end if ! GGA

!  ----------------------------------------------------------------------
!  ----- If non-collinear, rotate back potential according to magnetization
!  ----------------------------------------------------------------------
 if (option/=3.and.nspden==4) then
    ! Use of C pointers to avoid copies (when ISO C bindings are available)
    ! %@1$ xlf v15 compiler requires a auxilliary cptr variable
#ifdef LIBPAW_ISO_C_BINDING
   cptr=c_loc(vxc1_updn(1,1,1))
   call c_f_pointer(cptr,vxc1_diag,shape=[cplex_vxc*nrad*npts,nspden_updn])
   cptr=c_loc(vxc1_(1,1,1))
   call c_f_pointer(cptr,vxc1_nc,shape=[cplex_vxc*nrad*npts,nspden])
   cptr=c_loc(vxc(1,1,1))
   call c_f_pointer(cptr,vxc_,shape=[nrad*npts,nspden])
   cptr=c_loc(kxc(1,1,1))
   call c_f_pointer(cptr,kxc_,shape=[nrad*npts,3])
   cptr=c_loc(kxc(1,1,nkxc-2))
   call c_f_pointer(cptr,mag,shape=[nrad*npts,3])
#else
   LIBPAW_ALLOCATE(vxc1_diag,(cplex_vxc*nrad*npts,nspden_updn))
   LIBPAW_ALLOCATE(vxc1_nc,(cplex_vxc*nrad*npts,nspden))
   LIBPAW_ALLOCATE(vxc_,(nrad*npts,nspden))
   LIBPAW_ALLOCATE(kxc_,(nrad*npts,3))
   LIBPAW_ALLOCATE(mag,(nrad*npts,3))
   vxc1_diag=reshape(vxc1_updn,[cplex_vxc*nrad*npts,nspden_updn])
   vxc_=reshape(vxc(1:cplex_vxc*nrad,1:npts,1:nspden),[cplex_vxc*nrad*npts,nspden])
   kxc_=reshape(kxc(1:nrad,1:npts,1:3),[nrad*npts,3])
   mag=reshape(kxc(1:nrad,1:npts,nkxc-2:nkxc),[nrad*npts,3])
#endif
   !Should use cplex_den and cplex_vxc
   call pawxc_rotate_back_mag_dfpt(vxc1_diag,vxc1_nc,vxc_,kxc_,rho1_nc,mag,nrad*npts)
#ifndef LIBPAW_ISO_C_BINDING
   vxc1_=reshape(vxc1_nc,[cplex_vxc*nrad,npts,nspden])
   LIBPAW_DEALLOCATE(vxc1_diag)
   LIBPAW_DEALLOCATE(vxc1_nc)
   LIBPAW_DEALLOCATE(vxc_)
   LIBPAW_DEALLOCATE(kxc_)
   LIBPAW_DEALLOCATE(mag)
#endif
   LIBPAW_POINTER_DEALLOCATE(rho1_nc)
   LIBPAW_POINTER_DEALLOCATE(vxc1_updn)
 end if

!----------------------------------------------------------------------
!----- Accumulate and store 2nd-order change of XC energy
!----------------------------------------------------------------------
 if (option/=1) then

!  Do loop on the angular part (theta,phi)
   do ipts=1,npts

!    Copy the input 1st-order density for this (theta,phi)
     rho1arr(:,:)=zero
     if (usexcnhat< 1) rho1_=>rhor1
     if (usexcnhat>=1) rho1_=>rhohat1
     do ispden=1,nspden
       do ilm=1,lm_size_eff
         if (lmselect(ilm)) rho1arr(:,ispden)=rho1arr(:,ispden) &
  &       +rho1_(:,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
     end do
     if (usecore==1) then
       rho1arr(:,1)=rho1arr(:,1)+corexc1(:)
       if (nspden==2) rho1arr(:,2)=rho1arr(:,2)+half*corexc1(:)
     end if

!    ----- Calculate d2Exc=Int[Vxc^(1)^*(r).n^(1)(r).dr]
     LIBPAW_ALLOCATE(ff,(nrad))
     if (need_impart) then
       LIBPAW_ALLOCATE(gg,(nrad))
     end if

!    COLLINEAR MAGNETISM
     if (nspden/=4) then
       if (cplex_vxc==1.and.cplex_den==1) then       ! cplex_vxc==1 and cplex_den==1
         ff(:)=vxc1_(:,ipts,1)*rho1arr(:,nspden)
         if (nspden==2) ff(:)=ff(:)+vxc1_(:,ipts,2)*(rho1arr(:,1)-rho1arr(:,2))
         if (need_impart) gg(:)=zero
       else if (cplex_vxc==2.and.cplex_den==2) then  ! cplex_vxc==2 and cplex_den==2
         if (.not.need_impart) then      ! Real part only
           do ir=1,nrad
             jr=2*ir;v11r=vxc1_(jr-1,ipts,1);v11i=vxc1_(jr,ipts,1)
             ro11r=rho1arr(jr-1,nspden);ro11i=rho1arr(jr,nspden)
             ff(ir)=v11r*ro11r+v11i*ro11i
           end do
           if (nspden==2) then
             do ir=1,nrad
               jr=2*ir;v22r=vxc1_(jr-1,ipts,2);v22i=vxc1_(jr,ipts,2)
               ro22r=rho1arr(jr-1,1)-rho1arr(jr-1,2)
               ro22i=rho1arr(jr  ,1)-rho1arr(jr  ,2)
               ff(ir)=ff(ir)+v22r*ro22r+v22i*ro22i
             end do
           end if
         else
           do ir=1,nrad                  ! Real and imaginary parts
             jr=2*ir;v11r=vxc1_(jr-1,ipts,1);v11i=vxc1_(jr,ipts,1)
             ro11r=rho1arr(jr-1,nspden);ro11i=rho1arr(jr,nspden)
             ff(ir)=v11r*ro11r+v11i*ro11i
             gg(ir)=v11r*ro11i-v11i*ro11r
           end do
           if (nspden==2) then
             do ir=1,nrad
               jr=2*ir;v22r=vxc1_(jr-1,ipts,2);v22i=vxc1_(jr,ipts,2)
               ro22r=rho1arr(jr-1,1)-rho1arr(jr-1,2)
               ro22i=rho1arr(jr  ,1)-rho1arr(jr  ,2)
               ff(ir)=ff(ir)+v22r*ro22r+v22i*ro22i
               gg(ir)=gg(ir)+v22r*ro22i-v22i*ro22r
             end do
           end if
         end if
       else                                          ! other cases for cplex_vxc and cplex_den
         v11i=zero;ro11i=zero
         do ir=1,nrad
           jr=cplex_den*(ir-1)+1 ; kr=cplex_vxc*(ir-1)+1
           ro11r=rho1arr(jr,nspden);if (cplex_den==2) ro11i=rho1arr(jr+1,nspden)
           v11r=vxc1_(kr,ipts,1);if (cplex_vxc==2) v11i=vxc1_(kr+1,ipts,1)
           ff(ir)=v11r*ro11r+v11i*ro11i
           if (need_impart) gg(ir)=v11r*ro11i-v11i*ro11r
         end do
         if (nspden==2) then
           v22i=zero;ro22i=zero
           do ir=1,nrad
             jr=cplex_den*(ir-1)+1 ; kr=cplex_vxc*(ir-1)+1
             ro22r=rho1arr(jr,1)-rho1arr(jr,2)
             if (cplex_den==2) ro22i=rho1arr(jr+1,1)-rho1arr(jr+1,2)
             v22r=vxc1_(kr,ipts,2);if (cplex_vxc==2) v22i=vxc1_(kr+1,ipts,2)
             ff(ir)=ff(ir)+v22r*ro22r+v22i*ro22i
             gg(ir)=gg(ir)+v22r*ro22i-v22i*ro22r
           end do
         end if
       end if ! cplex_vxc and cplex_den

!      NON-COLLINEAR MAGNETISM
     else
       if (cplex_vxc==1.and.cplex_den==1) then   ! cplex_vxc==1 and cplex_den==1
         ff(:)=half*(vxc1_(:,ipts,1)*(rho1arr(:,1)+rho1arr(:,4)) &
&         +vxc1_(:,ipts,2)*(rho1arr(:,1)-rho1arr(:,4))) &
&         +vxc1_(:,ipts,3)*rho1arr(:,2) &
&         -vxc1_(:,ipts,4)*rho1arr(:,3)
         if (need_impart) gg(:)=zero
       else                                      ! other cases for cplex_vxc and cplex_den

!        V is stored as : v^11, v^22, V^12, i.V^21 (each are complex)
!        N is stored as : n, m_x, m_y, mZ          (each are complex)
         do ir=1,nrad
           jr=cplex_den*(ir-1)+1 ; kr=cplex_vxc*(ir-1)+1
           ro11r= rho1arr(jr,1)+rho1arr(jr,4)
           ro22r= rho1arr(jr,1)-rho1arr(jr,4)
           ro12r= rho1arr(jr,2);ro12i=-rho1arr(jr,3)
           ro21r= rho1arr(jr,2);ro21i= rho1arr(jr,3)
           if (cplex_den==2) then
             ro11i=rho1arr(jr+1,1)+rho1arr(jr+1,4)
             ro22i=rho1arr(jr+1,1)-rho1arr(jr+1,4)
             ro12r=ro12r+rho1arr(jr+1,3);ro12i=ro12i+rho1arr(jr+1,2)
             ro21r=ro21r-rho1arr(jr+1,3);ro21i=ro21i+rho1arr(jr+1,2)
           else
             ro11i=zero;ro22i=zero
           end if
           v11r= vxc1_(kr,ipts,1);v22r= vxc1_(kr,ipts,2)
           v12r= vxc1_(kr,ipts,3);v21i=-vxc1_(kr,ipts,1)
           if (cplex_vxc==2) then
             v11i= vxc1_(kr+1,ipts,1);v22i= vxc1_(kr+1,ipts,2)
             v12i= vxc1_(kr+1,ipts,3);v21r= vxc1_(kr+1,ipts,1)
           else
             v11i=zero;v22i=zero
             v12i=zero;v21i=zero
           end if
!          Real part
           ff(ir)=half*(v11r*ro11r+v11i*ro11i+v22r*ro22r+v22i*ro22i &
&                      +v12r*ro12r+v12i*ro12i+v21r*ro21r+v21i*ro21i)
!          Imaginary part
           if (need_impart) &
&            gg(ir)=half*(v11r*ro11i-v11i*ro11r+v22r*ro22i-v22i*ro22r &
&                        +v12r*ro12i-v12i*ro12r+v21r*ro21i-v21i*ro21r)
         end do
       end if ! cplex_vxc and cplex_den
     end if ! nspden

     ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
     call simp_gen(vxcrho,ff,pawrad)
     d2enxc=d2enxc+vxcrho*pawang%angwgth(ipts)
     LIBPAW_DEALLOCATE(ff)

     if (need_impart) then
       gg(1:nrad)=gg(1:nrad)*pawrad%rad(1:nrad)**2
       call simp_gen(vxcrho,gg,pawrad)
       d2enxc_im=d2enxc_im+vxcrho*pawang%angwgth(ipts)
       LIBPAW_DEALLOCATE(gg)
     end if

!    ----- End of the loop on npts (angular part)
   end do

 end if  ! option/=1

!Add the four*pi factor of the angular integration
 if (option/=1) then
   d2enxc=d2enxc*four_pi
   if (need_impart) d2enxc_im=d2enxc_im*four_pi
 end if

!Free memory
 if (usexcnhat>0)  then
   LIBPAW_DEALLOCATE(rhohat1)
 end if
 LIBPAW_DEALLOCATE(rho1arr)
 if (option==2) then
   LIBPAW_POINTER_DEALLOCATE(vxc1_)
 end if
 if (xclevel==2.and.option/=3) then
   LIBPAW_DEALLOCATE(gxc1)
   LIBPAW_DEALLOCATE(dylmdr)
 end if

end subroutine pawxc_dfpt
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxcsph
!! NAME
!! pawxcsph
!!
!! FUNCTION
!! Compute XC energy and potential for a spherical density rho(r) given as (up,dn)
!! Driver of XC functionals. Only treat collinear spins. LDA and GGA
!!
!! INPUTS
!!  exexch= choice of <<<local>>> exact exchange. Active if exexch>0 (only for GGA)
!!  ixc= choice of exchange-correlation scheme (see above and below)
!!  nkxc= size of kxc(nrad,nkxc) (XC kernel)
!!  nrad= dimension of the radial mesh
!!  nspden=number of spin-density components
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  rho_updn(nrad,lm_size,nspden)=electron density in real space
!!             up (ispden=1) and down (ispden=2) parts
!!             If nspden=1, rho_updn(:,:,1) contains (1/2).rho_total
!!  xclevel= XC functional level
!!
!! OUTPUT
!!  exc(nrad)= XC energy density
!!  vxc((nrad,nspden)= XC potential
!!  === Only if nkxc>0 ===
!!  kxc(nrad,nkxc)=exchange and correlation kernel (returned only if nkxc/=0)
!!   Content of Kxc array:
!!   ===== if LDA
!!    if nspden==1: kxc(:,1)= d2Exc/drho2
!!                 (kxc(:,2)= d2Exc/drho_up drho_dn)
!!    if nspden>=2: kxc(:,1)=d2Exc/drho_up drho_up
!!                  kxc(:,2)=d2Exc/drho_up drho_dn
!!                  kxc(:,3)=d2Exc/drho_dn drho_dn
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

 subroutine pawxcsph(exc,exexch,ixc,kxc,nkxc,nrad,nspden,pawrad,rho_updn,vxc,xclevel)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exexch,ixc,nkxc,nrad,nspden,xclevel
 type(pawrad_type),intent(in) :: pawrad
!arrays
 real(dp),intent(in) :: rho_updn(nrad,nspden)
 real(dp),intent(out) :: exc(nrad),kxc(nrad,nkxc),vxc(nrad,nspden)

!Local variables-------------------------------
!scalars
 integer :: ir,ispden,mgga,ndvxc,nd2vxc,ngr2,nspgrad,nvxcdgr,order,use_laplacian
 real(dp),parameter :: tol24=tol12*tol12
 real(dp) :: coeff,grho_tot,grho_up,fact
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: dff(:),dnexcdn(:,:),dvxcdgr(:,:),dvxci(:,:)
 real(dp),allocatable :: grho2(:,:),grho_updn(:,:)

! *************************************************************************

 if(nspden>2)then
   write(msg, '(a,a,a,i0)' )&
&   'Only non-spin-polarised or collinear spin-densities are allowed,',ch10,&
&   'while the argument nspden=',nspden
   MSG_BUG(msg)
 end if
 if(nkxc>3)then
   msg='nkxc>3 not allowed (GGA)!'
   MSG_ERROR(msg)
 end if
 if(nrad>pawrad%mesh_size)then
   msg='nrad > mesh size!'
   MSG_BUG(msg)
 end if

!Compute sizes of arrays and flags
 order=1;if (nkxc>0) order=2
 nspgrad=0;if (xclevel==2) nspgrad=3*nspden-1
 call pawxc_size_dvxc_wrapper(ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order)
 mgga=0 ; use_laplacian=0 !metaGGA contributions are not taken into account here


!--------------------------------------------------------------------------
!-------------- GGA: computation of the gradient of the density
!--------------------------------------------------------------------------

 LIBPAW_ALLOCATE(grho2,(nrad,ngr2))
 if (xclevel==2) then

!  grho_updn contains the gradient of the radial part
!  grho2(:,1:3) contains the squared norm of this gradient (up, dn and total)
   LIBPAW_ALLOCATE(grho_updn,(nrad,nspden))

!  Gradient of radial part of density
   LIBPAW_ALLOCATE(dff,(nrad))
   do ispden=1,nspden
     call nderiv_gen(dff,rho_updn(:,ispden),pawrad)
     grho_updn(:,ispden)=dff(:)
   end do
   LIBPAW_DEALLOCATE(dff)

!  Squared norm of the gradient
   grho2(:,1)=grho_updn(:,1)**2
   if (nspden==2) then
     grho2(:,2)=grho_updn(:,2)**2
     grho2(:,3)=(grho_updn(:,1)+grho_updn(:,2))**2
   end if

 end if

!--------------------------------------------------------------------------
!-------------- Computation of Exc, Vxc (and Kxc)
!--------------------------------------------------------------------------

!Allocate arrays
 LIBPAW_ALLOCATE(dvxci,(nrad,ndvxc))
 LIBPAW_ALLOCATE(dvxcdgr,(nrad,nvxcdgr))

!Call to main XC driver
 call pawxc_drivexc_wrapper(exc,ixc,mgga,ndvxc,nd2vxc,ngr2,nrad,nspden,nvxcdgr, &
&                           order,rho_updn,use_laplacian,vxc,xclevel, &
&                           dvxc=dvxci,exexch=exexch,grho2=grho2,vxcgrho=dvxcdgr)

!Transfer the XC kernel
 if (nkxc>0.and.ndvxc>0) then
   if (nkxc==1.and.ndvxc==15) then
     kxc(1:nrad,1)=half*(dvxci(1:nrad,1)+dvxci(1:nrad,9)+dvxci(1:nrad,10))
   else if (nkxc==3.and.ndvxc==15) then
     kxc(1:nrad,1)=dvxci(1:nrad,1)+dvxci(1:nrad,9)
     kxc(1:nrad,2)=dvxci(1:nrad,10)
     kxc(1:nrad,3)=dvxci(1:nrad,2)+dvxci(1:nrad,11)
   else if (nkxc==7.and.ndvxc==8) then
     kxc(1:nrad,1)=half*dvxci(1:nrad,1)
     kxc(1:nrad,2)=half*dvxci(1:nrad,3)
     kxc(1:nrad,3)=quarter*dvxci(1:nrad,5)
     kxc(1:nrad,4)=eighth*dvxci(1:nrad,7)
   else if (nkxc==7.and.ndvxc==15) then
     kxc(1:nrad,1)=half*(dvxci(1:nrad,1)+dvxci(1:nrad,9)+dvxci(1:nrad,10))
     kxc(1:nrad,2)=half*dvxci(1:nrad,3)+dvxci(1:nrad,12)
     kxc(1:nrad,3)=quarter*dvxci(1:nrad,5)+dvxci(1:nrad,13)
     kxc(1:nrad,4)=eighth*dvxci(1:nrad,7)+dvxci(1:nrad,15)
   else if (nkxc==19.and.ndvxc==15) then
     kxc(1:nrad,1)=dvxci(1:nrad,1)+dvxci(1:nrad,9)
     kxc(1:nrad,2)=dvxci(1:nrad,10)
     kxc(1:nrad,3)=dvxci(1:nrad,2)+dvxci(1:nrad,11)
     kxc(1:nrad,4)=dvxci(1:nrad,3)
     kxc(1:nrad,5)=dvxci(1:nrad,4)
     kxc(1:nrad,6)=dvxci(1:nrad,5)
     kxc(1:nrad,7)=dvxci(1:nrad,6)
     kxc(1:nrad,8)=dvxci(1:nrad,7)
     kxc(1:nrad,9)=dvxci(1:nrad,8)
     kxc(1:nrad,10)=dvxci(1:nrad,12)
     kxc(1:nrad,11)=dvxci(1:nrad,13)
     kxc(1:nrad,12)=dvxci(1:nrad,14)
     kxc(1:nrad,13)=dvxci(1:nrad,15)
   else ! Other cases
     kxc(1:nrad,1:nkxc)=zero
     kxc(1:nrad,1:min(nkxc,ndvxc))=dvxci(1:nrad,1:min(nkxc,ndvxc))
   end if
   if (nkxc==7) then
     kxc(1:nrad,5)=grho_updn(1:nrad,1)  ! Not correct
     kxc(1:nrad,6)=grho_updn(1:nrad,1)  ! Not correct
     kxc(1:nrad,7)=grho_updn(1:nrad,1)  ! Not correct
   else if (nkxc==19) then
     kxc(1:nrad,14)=grho_updn(1:nrad,1) ! Not correct
     kxc(1:nrad,15)=grho_updn(1:nrad,2) ! Not correct
     kxc(1:nrad,16)=grho_updn(1:nrad,1) ! Not correct
     kxc(1:nrad,17)=grho_updn(1:nrad,2) ! Not correct
     kxc(1:nrad,18)=grho_updn(1:nrad,1) ! Not correct
     kxc(1:nrad,19)=grho_updn(1:nrad,2) ! Not correct
   end if
 end if
 LIBPAW_DEALLOCATE(dvxci)

!--------------------------------------------------------------------------
!-------------- GGA: gardient corrections
!--------------------------------------------------------------------------

 if (xclevel==2.and.ixc/=13) then

!  Compute the derivative of Exc with respect to the (spin-)density,
!  or to the norm of the gradient of the (spin-)density,
!  Further divided by the norm of the gradient of the (spin-)density
!  The different components of dnexcdn will be
!  for nspden=1,         dnexcdn(:,1)=d(n.exc)/d(n)
!  and if xclevel=2, dnexcdn(:,2)=1/2*1/|grad n_up|*d(n.exc)/d(|grad n_up|)
!  +   1/|grad n|*d(n.exc)/d(|grad n|)
!  (do not forget : |grad n| /= |grad n_up| + |grad n_down|
!  for nspden=2,         dnexcdn(:,1)=d(n.exc)/d(n_up)
!  dnexcdn(:,2)=d(n.exc)/d(n_down)
!  and if xclevel=2, dnexcdn(:,3)=1/|grad n_up|*d(n.exc)/d(|grad n_up|)
!  dnexcdn(:,4)=1/|grad n_down|*d(n.exc)/d(|grad n_down|)
!  dnexcdn(:,5)=1/|grad n|*d(n.exc)/d(|grad n|)
   LIBPAW_ALLOCATE(dnexcdn,(nrad,nspgrad))
!  LDA term
   dnexcdn(:,1:nspden)=vxc(:,1:nspden)
!  Additional GGA terms
   do ir=1,nrad
     do ispden=1,3  ! spin_up, spin_down and total spin density
       if (nspden==1.and.ispden>=2) exit
!      If the norm of the gradient vanishes, then the different terms
!      vanishes, but the inverse of the gradient diverges,
!      so skip the update.
       if(grho2(ir,ispden)<tol24) then
         dnexcdn(ir,ispden+nspden)=zero;cycle
       end if
!      Compute the derivative of n.e_xc wrt the spin up, spin down,
!      or total density. In the non-spin-polarized case take the coeff.
!      that will be multiplied by the gradient of the total density.
       if (nvxcdgr/=0) then
         if (nspden==1) then
!          Definition of dvxcdgr changed in v3.3
           if (nvxcdgr==3) then
             coeff=half*dvxcdgr(ir,1)+dvxcdgr(ir,3)
           else
             coeff=half*dvxcdgr(ir,1)
           end if
         else if (nspden==2)then
           if (nvxcdgr==3) then
             coeff=dvxcdgr(ir,ispden)
           else if (ispden/=3) then
             coeff=dvxcdgr(ir,ispden)
           else if (ispden==3) then
             coeff=zero
           end if
         end if
       end if
       dnexcdn(ir,ispden+nspden)=coeff
     end do
   end do

!  Calculate grad(rho)*dnexcdn and put it in grho_updn(:,:)
   if (nvxcdgr/=0) then
     if(nspden==1)then
       grho_updn(:,1)=grho_updn(:,1)*dnexcdn(:,2)
     else
       do ir=1,nrad
         grho_up=grho_updn(ir,1);grho_tot=grho_up+grho_updn(ir,2)
         grho_updn(ir,1)=grho_up*dnexcdn(ir,3)+grho_tot*dnexcdn(ir,5)
         grho_updn(ir,2)=(grho_tot-grho_up)*dnexcdn(ir,4)+grho_tot*dnexcdn(ir,5)
       end do
     end if
   end if
   LIBPAW_DEALLOCATE(dnexcdn)

!  Compute Vxc
   LIBPAW_ALLOCATE(dff,(nrad))
   fact=one;if (nspden==1) fact=two
   do ispden=1,nspden
     call nderiv_gen(dff,grho_updn(:,ispden),pawrad)
     vxc(2:nrad,ispden)=vxc(2:nrad,ispden)-fact*(dff(2:nrad)+two*grho_updn(2:nrad,ispden)/pawrad%rad(2:nrad))
     call pawrad_deducer0(vxc(:,ispden),nrad,pawrad)
   end do
   LIBPAW_DEALLOCATE(dff)

 end if ! xclevel==2

!--------------------------------------------------------------------------
!-------------- Deallocations
!--------------------------------------------------------------------------

 LIBPAW_DEALLOCATE(grho2)
 LIBPAW_DEALLOCATE(dvxcdgr)
 if (xclevel==2)  then
   LIBPAW_DEALLOCATE(grho_updn)
 end if

end subroutine pawxcsph
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxcsph_dfpt
!! NAME
!! pawxcsph_dfpt
!!
!! FUNCTION
!! Compute XC 1st-order potential for a 1st-order spherical density rho1(r)
!! associated to a spherical density, both given as (up,dn)
!! Driver of XC functionals. Only treat collinear spins. LDA and GGA
!!
!! INPUTS
!!  cplex_den= if 1, 1st-order densities are REAL, if 2, COMPLEX
!!  cplex_vxc= if 1, 1st-order XC potential is complex, if 2, COMPLEX
!!  ixc= choice of exchange-correlation scheme (see above and below)
!!  nrad= dimension of the radial mesh
!!  nspden=number of spin-density components
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  rho_updn(nrad,lm_size,nspden)=electron density in real space
!!             up (ispden=1) and down (ispden=2) parts
!!             If nspden=1, rho_updn(:,:,1) contains (1/2).rho_total
!!  rho1_updn(nrad,lm_size,nspden)=electron 1st-order density in real space
!!             up (ispden=1) and down (ispden=2) parts
!!             If nspden=1, rho_updn(:,:,1) contains (1/2).rho1_total
!!  xclevel= XC functional level
!!
!! OUTPUT
!!  vxc1((nrad,nspden)= XC 1st-order potential
!!
!! PARENTS
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE


subroutine pawxcsph_dfpt(cplex_den,cplex_vxc,ixc,nrad,nspden,pawrad,rho_updn,rho1_updn,vxc1,xclevel)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex_den,cplex_vxc,ixc,nrad,nspden,xclevel
 type(pawrad_type),intent(in) :: pawrad
!arrays
 real(dp),intent(in) :: rho_updn(nrad,nspden),rho1_updn(cplex_den*nrad,nspden)
 real(dp),intent(out) :: vxc1(cplex_vxc*nrad,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ir,ispden,ivxc,jr,kr,mgga,ndvxc,nd2vxc,ngr2,ngrad,nkxc,nvxcdgr,order,use_laplacian
 real(dp),parameter :: tol24=tol12*tol12
!real(dp) :: coeff_grho_corr,coeff_grho_dn,coeff_grho_up,fact
!real(dp) :: grho_grho1,grho_grho1_dn,grho_grho1_up
 character(len=500) :: msg
!arrays
 integer,parameter :: ikxc(4)=(/1,2,2,3/),irho(4)=(/1,2,1,2/)
 real(dp),allocatable :: dff(:),dgg(:),dvxcdgr(:,:),dvxc(:,:),exc(:),ff(:),gg(:)
 real(dp),allocatable :: grho_updn(:,:),grho1_updn(:,:),grho2(:,:)
 real(dp),allocatable :: kxc(:,:),vxc(:,:)
!real(dp),allocatable :: gxc1i(:,:),gxc1r(:,:),vxc1i(:,:),vxc1r(:,:)

! *************************************************************************

 if(nspden>2)then
   write(msg, '(a,a,a,i0)' )&
&   'Only non-spin-polarised or collinear spin-densities are allowed,',ch10,&
&   'while the argument nspden=',nspden
   MSG_BUG(msg)
 end if
 if(nrad>pawrad%mesh_size)then
   msg='nrad > mesh size!'
   MSG_BUG(msg)
 end if

!Compute sizes of arrays and flags
 order=2 ! We need Kxc
 ngrad=1;if (xclevel==2) ngrad=2 ! ngrad=1 is for LDAs or LSDs; ngrad=2 is for GGAs
 call pawxc_size_dvxc_wrapper(ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order)
 nkxc=2*nspden-1;if (xclevel==2) nkxc=15 ! Not correct for nspden=1
 mgga=0 ; use_laplacian=0 !metaGGA contributions are not taken into account here

!--------------------------------------------------------------------------
!-------------- GGA: computation of the gradients of the densities
!--------------------------------------------------------------------------

 LIBPAW_ALLOCATE(grho2,(nrad,ngr2))
 if (ngrad==2) then

   LIBPAW_ALLOCATE(grho_updn,(nrad,nspden))
   LIBPAW_ALLOCATE(grho1_updn,(cplex_den*nrad,nspden))

!  Gradient of density
   LIBPAW_ALLOCATE(dff,(nrad))
   do ispden=1,nspden
     call nderiv_gen(dff,rho_updn(:,ispden),pawrad)
     grho_updn(:,ispden)=dff(:)
   end do
!  Gradient of 1st-order density
   if (cplex_den==1) then
     do ispden=1,nspden
       call nderiv_gen(dff,rho1_updn(:,ispden),pawrad)
       grho1_updn(:,ispden)=dff(:)
     end do
   else
     LIBPAW_ALLOCATE(ff,(nrad))
     LIBPAW_ALLOCATE(gg,(nrad))
     LIBPAW_ALLOCATE(dgg,(nrad))
     do ispden=1,nspden
       do ir=1,nrad
         ff(ir)=rho1_updn(2*ir-1,ispden)
         gg(ir)=rho1_updn(2*ir  ,ispden)
       end do
       call nderiv_gen(dff,ff,pawrad)
       call nderiv_gen(dgg,gg,pawrad)
       do ir=1,nrad
         grho1_updn(2*ir-1,ispden)=dff(ir)
         grho1_updn(2*ir  ,ispden)=dgg(ir)
       end do
     end do
     LIBPAW_DEALLOCATE(ff)
     LIBPAW_DEALLOCATE(gg)
     LIBPAW_DEALLOCATE(dgg)
   end if
   LIBPAW_DEALLOCATE(dff)

!  Squared norm of the gradient
   grho2(:,1)=grho_updn(:,1)**2
   if (nspden==2) then
     grho2(:,2)=grho_updn(:,2)**2
     grho2(:,3)=(grho_updn(:,1)+grho_updn(:,2))**2
   end if

 end if

!--------------------------------------------------------------------------
!-------------- Computation of Kxc (and Exc, Vxc)
!--------------------------------------------------------------------------

 LIBPAW_ALLOCATE(exc,(nrad))
 LIBPAW_ALLOCATE(vxc,(nrad,nspden))
 LIBPAW_ALLOCATE(dvxc,(nrad,ndvxc))
 LIBPAW_ALLOCATE(dvxcdgr,(nrad,nvxcdgr))

!Call to main XC driver
 call pawxc_drivexc_wrapper(exc,ixc,mgga,ndvxc,nd2vxc,ngr2,nrad,nspden,nvxcdgr, &
&                           order,rho_updn,use_laplacian,vxc,xclevel, &
&                           dvxc=dvxc,grho2=grho2,vxcgrho=dvxcdgr)

!Transfer the XC kernel
 LIBPAW_ALLOCATE(kxc,(nrad,nkxc))
 if (nkxc>0.and.ndvxc>0) then
   if (nkxc==1.and.ndvxc==15) then
     kxc(1:nrad,1)=half*(dvxc(1:nrad,1)+dvxc(1:nrad,9)+dvxc(1:nrad,10))
   else if (nkxc==3.and.ndvxc==15) then
     kxc(1:nrad,1)=dvxc(1:nrad,1)+dvxc(1:nrad,9)
     kxc(1:nrad,2)=dvxc(1:nrad,10)
     kxc(1:nrad,3)=dvxc(1:nrad,2)+dvxc(1:nrad,11)
   else if (nkxc==7.and.ndvxc==8) then
     kxc(1:nrad,1)=half*dvxc(1:nrad,1)
     kxc(1:nrad,2)=half*dvxc(1:nrad,3)
     kxc(1:nrad,3)=quarter*dvxc(1:nrad,5)
     kxc(1:nrad,4)=eighth*dvxc(1:nrad,7)
   else if (nkxc==7.and.ndvxc==15) then
     kxc(1:nrad,1)=half*(dvxc(1:nrad,1)+dvxc(1:nrad,9)+dvxc(1:nrad,10))
     kxc(1:nrad,2)=half*dvxc(1:nrad,3)+dvxc(1:nrad,12)
     kxc(1:nrad,3)=quarter*dvxc(1:nrad,5)+dvxc(1:nrad,13)
     kxc(1:nrad,4)=eighth*dvxc(1:nrad,7)+dvxc(1:nrad,15)
   else if (nkxc==19.and.ndvxc==15) then
     kxc(1:nrad,1)=dvxc(1:nrad,1)+dvxc(1:nrad,9)
     kxc(1:nrad,2)=dvxc(1:nrad,10)
     kxc(1:nrad,3)=dvxc(1:nrad,2)+dvxc(1:nrad,11)
     kxc(1:nrad,4)=dvxc(1:nrad,3)
     kxc(1:nrad,5)=dvxc(1:nrad,4)
     kxc(1:nrad,6)=dvxc(1:nrad,5)
     kxc(1:nrad,7)=dvxc(1:nrad,6)
     kxc(1:nrad,8)=dvxc(1:nrad,7)
     kxc(1:nrad,9)=dvxc(1:nrad,8)
     kxc(1:nrad,10)=dvxc(1:nrad,12)
     kxc(1:nrad,11)=dvxc(1:nrad,13)
     kxc(1:nrad,12)=dvxc(1:nrad,14)
     kxc(1:nrad,13)=dvxc(1:nrad,15)
   else ! Other cases
     kxc(1:nrad,1:nkxc)=zero
     kxc(1:nrad,1:min(nkxc,ndvxc))=dvxc(1:nrad,1:min(nkxc,ndvxc))
   end if
   if (nkxc==7) then
     kxc(1:nrad,5)=zero ! Not correct
     kxc(1:nrad,6)=zero ! Not correct
     kxc(1:nrad,7)=zero ! Not correct
   else if (nkxc==19) then
     kxc(1:nrad,14)=zero ! Not correct
     kxc(1:nrad,15)=zero ! Not correct
     kxc(1:nrad,16)=zero ! Not correct
     kxc(1:nrad,17)=zero ! Not correct
     kxc(1:nrad,18)=zero ! Not correct
     kxc(1:nrad,19)=zero ! Not correct
   end if
 end if

 LIBPAW_DEALLOCATE(exc)
 LIBPAW_DEALLOCATE(vxc)
 LIBPAW_DEALLOCATE(dvxc)
 LIBPAW_DEALLOCATE(dvxcdgr)

!--------------------------------------------------------------------------
!-------------- LDA
!--------------------------------------------------------------------------
 if (ngrad==1.or.ixc==13) then

   do ispden=1,3*nspden-2
     ivxc=1;if (ispden>2) ivxc=2
     if (cplex_vxc==1.and.cplex_den==1) then
       vxc1(:,ivxc)=vxc1(:,ivxc)+kxc(:,ikxc(ii))*rho1_updn(:,irho(ii))
     else
       do ir=1,nrad
         jr=cplex_den*(ir-1);kr=cplex_vxc*(ir-1)
         do ii=1,1+(cplex_den*cplex_vxc)/4
           jr=jr+1;kr=kr+1
           vxc1(kr,ivxc)=vxc1(kr,ivxc)+kxc(ir,ikxc(ii))*rho1_updn(jr,irho(ii))
         end do
       end do
     end if
   end do

!  --------------------------------------------------------------------------
!  -------------- GGA
!  --------------------------------------------------------------------------
 else

!  FOR NSPDEN=1, should eliminate computation of gxc1i(...), vxc1i(...)

!    LIBPAW_ALLOCATE(vxc1r,(nrad,2))
!    LIBPAW_ALLOCATE(vxc1i,(nrad,2))
!    LIBPAW_ALLOCATE(gxc1r,(nrad,2))
!    LIBPAW_ALLOCATE(gxc1i,(nrad,2))
!    do ir=1,nrad
!      if (cplex_vxc==1) then  ! cplex_vxc==1 and (cplex_den==1 or cplex_den=2)
!        jr=cplex_den*(ir-1)+1
!        grho_grho1_up=grho_updn(ir,1)*grho1_updn(jr,1)
!        grho_grho1_dn=grho_updn(ir,2)*grho1_updn(jr,2)
!        vxc1r(ir,1)=(kxc(ir, 1)+kxc(ir, 9))*rho1_updn(jr,1)+kxc(ir,10)*rho1_updn(jr,2) &
! &       +kxc(ir, 5)*grho_grho1_up+kxc(ir,13)*grho_grho1
!        vxc1r(ir,2)=(kxc(ir, 2)+kxc(ir,11))*rho1_updn(jr,2)+kxc(ir,10)*rho1_updn(jr,1) &
! &       +kxc(ir, 6)*grho_grho1_dn+kxc(ir,14)*grho_grho1
!        coeff_grho_corr=kxc(ir,13)*rho1_updn(jr,1)+kxc(ir,14)*rho1_updn(jr,2)+kxc(ir,15)*grho_grho1
!        coeff_grho_up  =kxc(ir, 5)*rho1_updn(jr,1)+kxc(ir, 7)*grho_grho1_up
!        coeff_grho_dn  =kxc(ir, 6)*rho1_updn(jr,2)+kxc(ir, 8)*grho_grho1_dn
!        gxc1r(ir,1)=(kxc(ir, 3)+kxc(ir,12))*grho1_updn(jr,1)+kxc(ir,12)*grho1_updn(jr,2) &
! &       +coeff_grho_up*grho_updn(jr,1)+coeff_grho_corr*(grho_updn(jr,1)+grho_updn(jr,2))
!        gxc1r(ir,2)=(kxc(ir, 4)+kxc(ir,12))*grho1_updn(jr,2)+kxc(ir,12)*grho1_updn(jr,1) &
! &       +coeff_grho_dn*grho_updn(jr,2)+coeff_grho_corr*(grho_updn(jr,1)+grho_updn(jr,2))
!      end if
!      if (grho2(ir,1)<tol24) gxc1r(ir,:)=zero ! ???
!    end do
!
! !  Apply divergence
!    fact=one;if (nspden==1) fact=two  ! Is it true  ? we force nspden=2 for gxc...
!    if (cplex_vxc==1) then
!      LIBPAW_ALLOCATE(dff,(nrad))
!      do ispden=1,nspden
!        call nderiv_gen(dff,gxc1r(:,ispden),pawrad)
!        vxc1(2:nrad,ispden)=vxc1r(2:nrad,ispden)-fact*(dff(2:nrad)+two*gxc1r(2:nrad,ispden)/pawrad%rad(2:nrad))
!        call pawrad_deducer0(vxc1(:,ispden),nrad,pawrad)
!      end do
!      LIBPAW_DEALLOCATE(dff)
!    else
!      LIBPAW_ALLOCATE(dff,(nrad))
!      LIBPAW_ALLOCATE(dgg,(nrad))
!      LIBPAW_ALLOCATE(ff,(nrad))
!      LIBPAW_ALLOCATE(gg,(nrad))
!      do ispden=1,nspden
!        call nderiv_gen(dff,gxc1r(:,ispden),pawrad)
!        call nderiv_gen(dgg,gxc1i(:,ispden),pawrad)
!        ff(2:nrad)=vxc1r(2:nrad,ispden)-fact*(dff(2:nrad)+two*gxc1r(2:nrad,ispden)/pawrad%rad(2:nrad))
!        gg(2:nrad)=vxc1i(2:nrad,ispden)-fact*(dgg(2:nrad)+two*gxc1i(2:nrad,ispden)/pawrad%rad(2:nrad))
!        call pawrad_deducer0(ff,nrad,pawrad)
!        call pawrad_deducer0(gg,nrad,pawrad)
!        do ir=1,nrad
!          vxc1(2*ir-1,ispden)=ff(ir)
!          vxc1(2*ir  ,ispden)=gg(ir)
!        end do
!      end do
!      LIBPAW_DEALLOCATE(dff)
!      LIBPAW_DEALLOCATE(dgg)
!      LIBPAW_DEALLOCATE(ff)
!      LIBPAW_DEALLOCATE(gg)
!    end if
!
!    LIBPAW_DEALLOCATE(vxc1r)
!    LIBPAW_DEALLOCATE(vxc1i)
!    LIBPAW_DEALLOCATE(gxc1r)
!    LIBPAW_DEALLOCATE(gxc1i)

 end if ! ngrad==2

!--------------------------------------------------------------------------
!-------------- Deallocations
!--------------------------------------------------------------------------

 LIBPAW_DEALLOCATE(grho2)
 LIBPAW_DEALLOCATE(kxc)
 if (ngrad==2) then
   LIBPAW_DEALLOCATE(grho_updn)
   LIBPAW_DEALLOCATE(grho1_updn)
 end if

end subroutine pawxcsph_dfpt
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxcsphpositron
!! NAME
!! pawxcsphpositron
!!
!! FUNCTION
!! Compute electron-positron XC energy and potential for spherical densities rho_el(r) rho_pos(r)
!! Driver of XC functionals. LDA and GGA
!!
!! INPUTS
!!  calctype=type of electron-positron calculation:
!!           calctype=1 : positron in electronic density
!!           calctype=2 : electrons in positronic density
!!  ixcpositron= choice of elctron-positron exchange-correlation scheme
!!  nrad= dimension of the radial mesh
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  posdensity0_limit=True if we are in the zero positron density limit
!!  rho(nrad,lm_size)=electron (or positron) density in real space
!!                    Contents depends on calctype value:
!!                    calctype=1: rho is the positronic density
!!                    calctype=2: rho is the electronic density
!!  rho_ep(nrad,lm_size)=electron (or positron) density in real space
!!                      Contents depends on calctype value:
!!                      calctype=1: rho_ep is the electronic density
!!                      calctype=2: rho_ep is the positronic density
!!
!! OUTPUT
!!  fxc(nrad)= electron-positron XC energy per unit volume
!!  vxce(nrad)= electron-positron XC potential for the electron
!!  vxcp(nrad)= electron-positron XC potential for the positron
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

 subroutine pawxcsphpositron(calctype,fxc,ixcpositron,nrad,pawrad,posdensity0_limit,rho,rho_ep,vxce,vxcp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: calctype,ixcpositron,nrad
 logical,intent(in) :: posdensity0_limit
 type(pawrad_type),intent(in) :: pawrad
!arrays
 real(dp),intent(in) :: rho(nrad),rho_ep(nrad)
 real(dp),intent(out) :: fxc(nrad),vxce(nrad),vxcp(nrad)

!Local variables-------------------------------
!scalars
 integer :: ngr
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: dff(:),rhograd(:),rhograd2(:),vxcegr(:)

! *************************************************************************

 if(nrad>pawrad%mesh_size)then
   msg='nrad > mesh size!'
   MSG_BUG(msg)
 end if

!Need gradient of density for GGA
 ngr=0;if (ixcpositron==3.or.ixcpositron==31) ngr=nrad
 LIBPAW_ALLOCATE(rhograd,(ngr))
 LIBPAW_ALLOCATE(rhograd2,(ngr))
 LIBPAW_ALLOCATE(vxcegr,(ngr))
 if (ngr==nrad) then
   if (calctype==1) then
     call nderiv_gen(rhograd,rho_ep,pawrad)
   else if (calctype==2) then
     call nderiv_gen(rhograd,rho,pawrad)
   end if
   rhograd2(:)=rhograd(:)**2
 end if

!---- Computation of Fxc and Vxc for the positron
!rho    is the positronic density
!rho_ep is the electronic density
 if (calctype==1) then
   call pawxc_xcpositron_wrapper(fxc,rhograd2,ixcpositron,ngr,nrad,posdensity0_limit,rho_ep,rho,vxce,vxcegr,vxcp)

!  ---- Computation of Exc and Vxc for the electron
!  rho    is the electronic density
!  rho_ep is the positronic density
 else if (calctype==2) then
   call pawxc_xcpositron_wrapper(fxc,rhograd2,ixcpositron,ngr,nrad,posdensity0_limit,rho,rho_ep,vxce,vxcegr,vxcp)
 end if

 LIBPAW_DEALLOCATE(rhograd2)

!---- GGA - gradient corrections
 if (ngr==nrad) then
   LIBPAW_ALLOCATE(dff,(nrad))
   vxcegr(1:nrad)=vxcegr(1:nrad)*rhograd(1:nrad)
   call nderiv_gen(dff,vxcegr,pawrad)
   vxcp(2:nrad)=vxcp(2:nrad)-(dff(2:nrad)+two*vxcegr(2:nrad)/pawrad%rad(2:nrad))
   call pawrad_deducer0(vxcp,nrad,pawrad)
   LIBPAW_DEALLOCATE(dff)
 end if

 LIBPAW_DEALLOCATE(vxcegr)
 LIBPAW_DEALLOCATE(rhograd)

end subroutine pawxcsphpositron
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxcsum
!! NAME
!! pawxcsum
!!
!! FUNCTION
!! Compute useful sums of moments of densities needed to compute on-site contributions to XC energy and potential
!!  First order sums:
!!    Sum1(1)=Sum_L{Rho1_L(r)**2}
!!    Sum1(2)=Sum_L{Rho1_L(r)*Rho2_L(r)}
!!    Sum1(3)=Sum_L{Rho2_L(r)**2}
!!    With L>0
!!  Second order sums:
!!    Sum2(L,1)=Sum_L1_L2{Rho1_L1(r)*Rho1_L1(r)*Gaunt_(L,L1,L2)}
!!    Sum2(L,2)=Sum_L1_L2{Rho1_L1(r)*Rho2_L2(r)*Gaunt_(L,L1,L2)}
!!    Sum2(L,3)=Sum_L1_L2{Rho2_L2(r)*Rho2_L2(r)*Gaunt_(L,L1,L2)}
!!    With L1>0, L2>0
!!
!! INPUTS
!!  cplex1=if 1, density Rho1 is REAL, if 2, COMPLEX
!!  cplex2=if 1, density Rho2 is REAL, if 2, COMPLEX
!!  cplexsum=if 1, output sums (Sum1 and Sum2) are REAL, if 2, COMPLEX
!!  lmselect1(lm_size)=select the non-zero LM-moments of input density Rho1
!!  lmselect2(lm_size)=select the non-zero LM-moments of input density Rho2
!!  lm_size=number of moments of the density
!!  nrad=number of radial points
!!  nsums=number of sums to compute:
!!        if nsums=1, computes only
!!                    Sum1(1)=Sum_L{Rho1_L(r)*Rho2_L(r)}
!!                    Sum2(L,1)=Sum_L1_L2{Rho1_L1(r)*Rho2_L2(r)*Gaunt_(L,L1,L2)}
!!        if nsums=3, computes all sums (Sum1(1:3), Sum2(1:3)
!!  option= 1: compute first order sums
!!          2: compute first and second order sums
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  rho1(cplex1*nrad,lm_size)=moments of first density on each radial point
!!  rho2(cplex2*nrad,lm_size)=moments of 2nd density on each radial point
!!
!! OUTPUT
!!  sum1(cplexsum*nrad,nsums)=first order sums
!!  === if option>=2
!!    sum2(cplexsum*nrad,lm_size,nsums)=second order sums
!!
!! PARENTS
!!      m_pawxc,poslifetime,posratecore
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

 subroutine pawxcsum(cplex1,cplex2,cplexsum,lmselect1,lmselect2,lm_size,nrad,nsums,&
&                    option,pawang,rho1,rho2,sum1,sum2)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex1,cplex2,cplexsum,lm_size,nrad,nsums,option
!arrays
 logical,intent(in) :: lmselect1(lm_size),lmselect2(lm_size)
 real(dp),intent(in) :: rho1(cplex1*nrad,lm_size),rho2(cplex2*nrad,lm_size)
 real(dp),intent(out) :: sum1(cplexsum*nrad,nsums),sum2(cplexsum*nrad,lm_size,nsums*(option/2))
 type(pawang_type),intent(in) :: pawang

!Local variables-------------------------------
!scalars
 integer :: ilm,ilm1,ilm2,ir,i1r,i2r,i3r,isel
 real(dp) :: fact,ro1i,ro1r,ro2i,ro2r
 character(len=500) :: msg
!arrays

!************************************************************************

 if(nsums/=1.and.nsums/=3) then
   msg='nsums must be 1 or 3!'
   MSG_BUG(msg)
 end if
 if(pawang%gnt_option==0) then
   msg='pawang%gnt_option=0!'
   MSG_BUG(msg)
 end if

 if (option>=1) then

!  SUM1(r)= Sum_L{Rho1_L(r)*Rho2_L(r)} (L>0)
!  --------------------------------------------------
   sum1=zero

!  ===== All input/output densities are REAL ====
   if (cplex1==1.and.cplex2==1.and.cplexsum==1) then
!    One sum to compute
     if (nsums==1) then
       do ilm=2,lm_size
         if (lmselect1(ilm).and.lmselect2(ilm)) then
           sum1(:,1)=sum1(:,1)+rho1(:,ilm)*rho2(:,ilm)
         end if
       end do
!      Three sums to compute
     else
       do ilm=2,lm_size
         if (lmselect1(ilm)) then
           sum1(:,1)=sum1(:,1)+rho1(:,ilm)**2
           if (lmselect2(ilm)) sum1(:,2)=sum1(:,2)+rho1(:,ilm)*rho2(:,ilm)
         end if
         if (lmselect2(ilm)) sum1(:,3)=sum1(:,3)+rho2(:,ilm)**2
       end do
     end if

!    ===== At least one of Rho1 and Rho2 is COMPLEX ====
   else
!    One sum to compute
     if (nsums==1) then
       do ilm=2,lm_size
         if (lmselect1(ilm).and.lmselect2(ilm)) then
           do ir=1,nrad
             i1r=cplex1*(ir-1)+1;i2r=cplex2*(ir-1)+1;i3r=cplexsum*(ir-1)+1
             ro1r=rho1(i1r,ilm);ro1i=zero;if (cplex1==2) ro1i=rho1(i1r+1,ilm)
             ro2r=rho2(i2r,ilm);ro2i=zero;if (cplex2==2) ro2i=rho2(i2r+1,ilm)
             sum1(i3r,1)=sum1(i3r,1)+ro1r*ro2r-ro1i*ro2i
             if (cplexsum==2) sum1(i3r+1,1)=sum1(i3r+1,1)+ro1r*ro2i+ro1i*ro2r
           end do
         end if
       end do
!      Three sums to compute
     else
       do ilm=2,lm_size
         do ir=1,nrad
           i1r=cplex1*(ir-1)+1;i2r=cplex2*(ir-1)+1;i3r=cplexsum*(ir-1)+1
           ro1r=rho1(i1r,ilm);ro1i=zero;if (cplex1==2) ro1i=rho1(i1r+1,ilm)
           ro2r=rho2(i2r,ilm);ro2i=zero;if (cplex2==2) ro2i=rho2(i2r+1,ilm)
           if (lmselect1(ilm)) then
             sum1(i3r,1)=sum1(i3r,1)+ro1r**2-ro1i**2
             if (lmselect2(ilm)) sum1(i3r,2)=sum1(i3r,2)+ro1r*ro2r-ro1i*ro2i
           end if
           if (lmselect2(ilm)) sum1(i3r,3)=sum1(i3r,3)+ro2r**2-ro2i**2
           if (cplexsum==2) then
             if (lmselect1(ilm)) then
               sum1(i3r+1,1)=sum1(i3r+1,1)+two*ro1r*ro1i
               if (lmselect2(ilm)) sum1(i3r+1,2)=sum1(i3r+1,2)+ro1r*ro2i+ro1i*ro2r
             end if
             if (lmselect2(ilm)) sum1(i3r+1,3)=sum1(i3r+1,3)+two*ro2r*ro2i
           end if
         end do
       end do
     end if ! nsums
   end if  ! cplex

 end if !option

 if (option>=2) then

!  SUM2(r,L)= Sum_L1_L2{Rho1_L1(r)*Rho2_L2(r)*Gaunt_(L,L1,L2)}  (L1>0, L2>0)
!  --------------------------------------------------
   sum2=zero
!  ===== All input/output densities are REAL ====
   if (cplex1==1.and.cplex2==1.and.cplexsum==1) then
!    One sum to compute
     if (nsums==1) then
       do ilm=1,lm_size
         do ilm1=2,lm_size
           if (lmselect1(ilm1)) then
             do ilm2=2,ilm1
               if (lmselect2(ilm2)) then
                 isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
                 if (isel>0) then
                   fact=pawang%realgnt(isel);if (ilm1/=ilm2) fact=two*fact
                   sum2(:,ilm,1)=sum2(:,ilm,1)+fact*rho1(:,ilm1)*rho2(:,ilm2)
                 end if
               end if
             end do
           end if
         end do
       end do
!      Three sums to compute
     else
       do ilm=1,lm_size
         do ilm1=2,lm_size
           if (lmselect1(ilm1)) then
             do ilm2=2,ilm1
               if (lmselect1(ilm2)) then
                 isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
                 if (isel>0) then
                   fact=pawang%realgnt(isel);if (ilm1/=ilm2) fact=two*fact
                   sum2(:,ilm,1)=sum2(:,ilm,1)+fact*rho1(:,ilm1)*rho1(:,ilm2)
                 end if
               end if
             end do
           end if
         end do
         do ilm1=2,lm_size
           if (lmselect2(ilm1)) then
             do ilm2=2,ilm1
               if (lmselect2(ilm2)) then
                 isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
                 if (isel>0) then
                   fact=pawang%realgnt(isel);if (ilm1/=ilm2) fact=two*fact
                   sum2(:,ilm,3)=sum2(:,ilm,3)+fact*rho2(:,ilm1)*rho2(:,ilm2)
                 end if
               end if
             end do
           end if
         end do
         do ilm1=2,lm_size
           if (lmselect1(ilm1)) then
             do ilm2=2,ilm1
               if (lmselect2(ilm2)) then
                 isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
                 if (isel>0) then
                   fact=pawang%realgnt(isel)
                   sum2(:,ilm,2)=sum2(:,ilm,2)+fact*rho1(:,ilm1)*rho2(:,ilm2)
                 end if
               end if
             end do
             if (ilm1<lm_size) then
               do ilm2=ilm1+1,lm_size
                 if (lmselect2(ilm2)) then
                   isel=pawang%gntselect(ilm,ilm1+ilm2*(ilm2-1)/2)
                   if (isel>0) then
                     fact=pawang%realgnt(isel)
                     sum2(:,ilm,2)=sum2(:,ilm,2)+fact*rho1(:,ilm1)*rho2(:,ilm2)
                   end if
                 end if
               end do
             end if
           end if
         end do
       end do
     end if ! nsums

!    ===== At least one of Rho1 and Rho2 is COMPLEX ====
   else
!    One sum to compute
     if (nsums==1) then
       do ilm=1,lm_size
         do ilm1=2,lm_size
           if (lmselect1(ilm1)) then
             do ilm2=2,ilm1
               if (lmselect2(ilm2)) then
                 isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
                 if (isel>0) then
                   fact=pawang%realgnt(isel);if (ilm1/=ilm2) fact=two*fact
                   do ir=1,nrad
                     i1r=cplex1*(ir-1)+1;i2r=cplex2*(ir-1)+1;i3r=cplexsum*(ir-1)+1
                     ro1r=rho1(i1r,ilm1);ro1i=zero;if (cplex1==2) ro1i=rho1(i1r+1,ilm1)
                     ro2r=rho2(i2r,ilm2);ro2i=zero;if (cplex2==2) ro2i=rho2(i2r+1,ilm2)
                     sum2(i3r,ilm,1)=sum2(i3r,ilm,1)+fact*(ro1r*ro2r-ro1i*ro2i)
                     if (cplexsum==2) sum2(i3r+1,ilm,1)=sum2(i3r+1,ilm,1)+fact*(ro1r*ro2i+ro1i*ro2r)
                   end do
                 end if
               end if
             end do
           end if
         end do
       end do
!      Three sums to compute
     else
       do ilm=2,lm_size
         do ir=1,nrad
           i1r=cplex1*(ir-1)+1;i2r=cplex2*(ir-1)+1;i3r=cplexsum*(ir-1)+1
           do ilm1=2,lm_size
             if (lmselect1(ilm1)) then
               ro1r=rho1(i1r,ilm1);ro1i=zero;if (cplex1==2) ro1i=rho1(i1r+1,ilm1)
               do ilm2=2,ilm1
                 if (lmselect1(ilm2)) then
                   isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
                   if (isel>0) then
                     fact=pawang%realgnt(isel);if (ilm1/=ilm2) fact=two*fact
                     ro2r=rho1(i1r,ilm2);ro2i=zero;if (cplex1==2) ro2i=rho1(i1r+1,ilm2)
                     sum2(i3r,ilm,1)=sum2(i3r,ilm,1)+fact*(ro1r*ro2r-ro1i*ro2i)
                     if (cplexsum==2) sum2(i3r+1,ilm,1)=sum2(i3r+1,ilm,1)+fact*(ro1r*ro2i+ro1i*ro2r)
                   end if
                 end if
               end do
             end if
           end do
           do ilm1=2,lm_size
             if (lmselect2(ilm1)) then
               ro1r=rho2(i2r,ilm1);ro1i=zero;if (cplex2==2) ro1i=rho2(i2r+1,ilm1)
               do ilm2=2,ilm1
                 if (lmselect2(ilm2)) then
                   isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
                   if (isel>0) then
                     fact=pawang%realgnt(isel);if (ilm1/=ilm2) fact=two*fact
                     ro2r=rho2(i2r,ilm2);ro2i=zero;if (cplex2==2) ro2i=rho2(i2r+1,ilm2)
                     sum2(i3r,ilm,3)=sum2(i3r,ilm,3)+fact*(ro1r*ro2r-ro1i*ro2i)
                     if (cplexsum==2) sum2(i3r+1,ilm,3)=sum2(i3r+1,ilm,3)+fact*(ro1r*ro2i+ro1i*ro2r)
                   end if
                 end if
               end do
             end if
           end do
           do ilm1=2,lm_size
             if (lmselect1(ilm1)) then
               ro1r=rho1(i1r,ilm1);ro1i=zero;if (cplex1==2) ro1i=rho1(i1r+1,ilm1)
               do ilm2=2,ilm1
                 if (lmselect2(ilm2)) then
                   isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
                   if (isel>0) then
                     fact=pawang%realgnt(isel)
                     ro2r=rho2(i2r,ilm2);ro2i=zero;if (cplex2==2) ro2i=rho2(i2r+1,ilm2)
                     sum2(i3r,ilm,2)=sum2(i3r,ilm,2)+fact*(ro1r*ro2r-ro1i*ro2i)
                     if (cplexsum==2) sum2(i3r+1,ilm,2)=sum2(i3r+1,ilm,2)+fact*(ro1r*ro2i+ro1i*ro2r)
                   end if
                 end if
               end do
               if (ilm1<lm_size) then
                 do ilm2=ilm1+1,lm_size
                   if (lmselect2(ilm2)) then
                     isel=pawang%gntselect(ilm,ilm1+ilm2*(ilm2-1)/2)
                     if (isel>0) then
                       fact=pawang%realgnt(isel)
                       ro2r=rho2(i2r,ilm2);ro2i=zero;if (cplex2==2) ro2i=rho2(i2r+1,ilm2)
                       sum2(i3r,ilm,2)=sum2(i3r,ilm,2)+fact*(ro1r*ro2r-ro1i*ro2i)
                       if (cplexsum==2) sum2(i3r+1,ilm,2)=sum2(i3r+1,ilm,2)+fact*(ro1r*ro2i+ro1i*ro2r)
                     end if
                   end if
                 end do
               end if
             end if
           end do
         end do
       end do
     end if ! nsums

   end if  ! cplex

 end if !option

 end subroutine pawxcsum
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxcm
!! NAME
!! pawxcm
!!
!! FUNCTION
!! Start from the density or spin-density, and compute xc correlation
!! potential and energies inside a paw sphere.
!! LDA+GGA - USE A DEVELOPMENT OF THE DENSITY OVER (L,M) MOMENTS
!! Driver of XC functionals.
!!
!! INPUTS
!!  corexc(nrad)=core density on radial grid
!!  exexch= choice of <<<local>>> exact exchange. Active if exexch=3 (only for PBE)
!!  ixc= choice of exchange-correlation scheme
!!  lm_size=size of density array rhor (see below)
!!  lmselect(lm_size)=select the non-zero LM-moments of input density rhor
!!  nhat(nrad,lm_size,nspden)=compensation density
!!                                        (total in 1st half and spin-up in 2nd half if nspden=2)
!!  nkxc=second dimension of the kxc array. If /=0, the exchange-correlation kernel must be computed
!!  non_magnetic_xc= if true, handle density/potential as non-magnetic (even if it is)
!!  nrad=size of radial mesh for densities/potentials (might be different from pawrad%mesh_size)
!!  nspden=number of spin-density components
!!  option=0 compute both XC energies (direct+double-counting) and potential (and Kernel)
!!         1 compute only XC potential (and Kernel)
!!         2 compute only XC energies (direct+double-counting)
!!         3 compute only XC energy by direct scheme
!!         4 compute only XC energy by direct scheme for spherical part of the density
!!         5 compute only XC potential (and Kernel) for spherical part of the density
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawxcdev=order of Vxc development
!!  rhor(nrad,lm_size,nspden)=electron density in real space in electrons/bohr**3
!!                                       (total in 1st half and spin-up in 2nd half if nspden=2)
!!  usecore= 1 if core density has to be used in Exc/Vxc ; 0 otherwise
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in double counting energy term only
!!             2 if compensation density (nhat) has to be used in Exc/Vxc and double counting energy term
!!  xclevel= XC functional level
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!
!! OUTPUT
!!  == if option==0, 2, 3, or 4 ==
!!    enxc=returned exchange and correlation energy (hartree)
!!  == if option==0 or 2 ==
!!    enxcdc=returned exchange-cor. contribution to double-counting energy
!!  == if option==0 or 1 ==
!!    vxc(nrad,lm_size,nspden)=xc potential
!!       (spin up in 1st half and spin-down in 2nd half if nspden=2)
!!  == if nkxc>0 ==
!!    kxc(nrad,lm_size,nkxc)=xc kernel (see notes below for nkxc)
!!
!! NOTES
!!  Content of Kxc array:
!!   ===== if LDA
!!    if nspden==1: kxc(:,1)= d2Exc/drho2
!!                 (kxc(:,2)= d2Exc/drho_up drho_dn)
!!    if nspden>=2: kxc(:,1)= d2Exc/drho_up drho_up
!!                  kxc(:,2)= d2Exc/drho_up drho_dn
!!                  kxc(:,3)= d2Exc/drho_dn drho_dn
!!    if nspden==4: kxc(:,4:6)= (m_x, m_y, m_z) (magnetization)
!!   ===== if GGA
!!    if nspden==1:
!!       kxc(:,1)= d2Exc/drho2
!!       kxc(:,2)= 1/|grad(rho)| dExc/d|grad(rho)|
!!       kxc(:,3)= 1/|grad(rho)| d2Exc/d|grad(rho)| drho
!!       kxc(:,4)= 1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dExc/d|grad(rho)| )
!!       kxc(:,5)= gradx(rho)
!!       kxc(:,6)= grady(rho)
!!       kxc(:,7)= gradz(rho)
!!    if nspden>=2:
!!       kxc(:,1)= d2Exc/drho_up drho_up
!!       kxc(:,2)= d2Exc/drho_up drho_dn
!!       kxc(:,3)= d2Exc/drho_dn drho_dn
!!       kxc(:,4)= 1/|grad(rho_up)| dEx/d|grad(rho_up)|
!!       kxc(:,5)= 1/|grad(rho_dn)| dEx/d|grad(rho_dn)|
!!       kxc(:,6)= 1/|grad(rho_up)| d2Ex/d|grad(rho_up)| drho_up
!!       kxc(:,7)= 1/|grad(rho_dn)| d2Ex/d|grad(rho_dn)| drho_dn
!!       kxc(:,8)= 1/|grad(rho_up)| * d/d|grad(rho_up)| ( 1/|grad(rho_up)| dEx/d|grad(rho_up)| )
!!       kxc(:,9)= 1/|grad(rho_dn)| * d/d|grad(rho_dn)| ( 1/|grad(rho_dn)| dEx/d|grad(rho_dn)| )
!!       kxc(:,10)=1/|grad(rho)| dEc/d|grad(rho)|
!!       kxc(:,11)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_up
!!       kxc(:,12)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_dn
!!       kxc(:,13)=1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dEc/d|grad(rho)| )
!!       kxc(:,14)=gradx(rho_up)
!!       kxc(:,15)=gradx(rho_dn)
!!       kxc(:,16)=grady(rho_up)
!!       kxc(:,17)=grady(rho_dn)
!!       kxc(:,18)=gradz(rho_up)
!!       kxc(:,19)=gradz(rho_dn)
!!    if nspden==4:
!!       kxc(:,20:22)= (m_x, m_y, m_z) (magnetization)
!!
!! PARENTS
!!      m_pawpsp,pawdenpot
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

 subroutine pawxcm(corexc,enxc,enxcdc,exexch,ixc,kxc,lm_size,lmselect,nhat,nkxc,&
&                  non_magnetic_xc,nrad,nspden,option,pawang,pawrad,pawxcdev,rhor,&
&                  usecore,usexcnhat,vxc,xclevel,xc_denpos)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exexch,ixc,lm_size,nkxc,nrad,nspden,option,pawxcdev,usecore
 integer,intent(in) :: usexcnhat,xclevel
 logical,intent(in) :: non_magnetic_xc
 real(dp),intent(in) :: xc_denpos
 real(dp),intent(out) :: enxc,enxcdc
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size)
 real(dp),intent(in) :: corexc(nrad)
 real(dp),intent(in) :: nhat(nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: rhor(nrad,lm_size,nspden)
 real(dp),intent(out) :: kxc(nrad,lm_size,nkxc)
 real(dp),intent(out) :: vxc(nrad,lm_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ilm,ir,ir1,ir2,ispden,iwarn,jr,nspden_updn,nsums
 real(dp),parameter :: delta=1.d-4
 real(dp) :: dvxc1,dvxc2,dvxc3,dvxc4,dvxca,dvxcb,dvxcc,dvxcd
 real(dp) :: fact,invsqfpi,invsqfpi2,sqfpi,sqfpi2,tol_rho
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: d1kxc(:,:),d2kxc(:,:),d1vxc(:,:),d2vxc(:,:)
 real(dp),allocatable :: exc_(:),exci(:),ff(:),gg(:)
 real(dp),allocatable :: kxc1(:,:),kxc2(:,:),kxcdn1(:,:),kxcdn2(:,:),kxci(:,:)
 real(dp),allocatable :: m_norm_inv(:),rho_(:,:),rhoinv(:,:),rhosph(:,:)
 real(dp),allocatable :: v1sum(:,:),v2sum(:,:,:)
 real(dp),allocatable :: vxc1(:,:),vxc2(:,:),vxcdn1(:,:),vxcdn2(:,:),vxci(:,:)
 real(dp),allocatable,target :: rho_nc(:,:),rho_updn(:,:,:),vxc_diag(:,:),vxc_nc(:,:)
 real(dp), LIBPAW_CONTIGUOUS pointer :: mag_nc(:,:),rho_dn(:,:),rho_up(:,:)

!************************************************************************

 if(nkxc>3) then
   msg='Kxc not implemented for GGA!'
   MSG_ERROR(msg)
 end if
 if(nkxc>0.and.nspden==4) then
   msg='Kxc not implemented for non-collinear magnetism!'
   MSG_ERROR(msg)
 end if
 if (option/=1.and.option/=5) then
   if (nrad<pawrad%int_meshsz) then
     msg='When option=0,2,3,4, nrad must be greater than pawrad%int_meshsz!'
     MSG_BUG(msg)
   end if
 end if

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------

!Arrays dimensions and constants
 iwarn=0
 nspden_updn=min(nspden,2)
 sqfpi=sqrt(four_pi);sqfpi2=half*sqfpi
 invsqfpi=one/sqfpi;invsqfpi2=half*invsqfpi
 nsums=2*nspden_updn-1

!Initializations of output arrays
 if (option/=1.and.option/=5) enxc=zero
 if (option==0.or.option==2) enxcdc=zero
 if (option/=3.and.option/=4) vxc(:,:,:)=zero
 if (nkxc/=0) kxc(:,:,:)=zero

 if (xclevel==0.or.ixc==0) then ! No xc at all is applied (usually for testing)
   msg='Note that no xc is applied (ixc=0). Returning'
   MSG_WARNING(msg)
   return
 end if

!----------------------------------------------------------------------
!----- Build several densities
!----------------------------------------------------------------------

!rho_updn contains the effective density used for XC
!with core density and/or compensation density eventually included
!-----------------------------------------------------------------

 LIBPAW_ALLOCATE(rho_updn,(nrad,lm_size,nspden))
 rho_updn(:,:,:)=rhor(:,:,:)
 if (usexcnhat==2) rho_updn(:,:,:)=rho_updn(:,:,:)+nhat(:,:,:)

!Optionally suppressed magnetic part
 if(non_magnetic_xc) then
   if(nspden==2) rho_updn(:,:,2)=rho_updn(:,:,1)*half
   if(nspden==4) rho_updn(:,:,2:4)=zero
 endif

!Add core density
 if (usecore==1) then
   if (nspden==1.or.nspden==4) then
     rho_updn(:,1,1)=rho_updn(:,1,1)+sqfpi*corexc(:)
   else if (nspden==2) then
     rho_updn(:,1,1)=rho_updn(:,1,1)+sqfpi*corexc(:)
     rho_updn(:,1,2)=rho_updn(:,1,2)+sqfpi2*corexc(:)
   end if
 end if

!In case of collinear magnetism, separate up and down contributions
 if (nspden==2) then
   LIBPAW_ALLOCATE(ff,(nrad))
   do ilm=1,lm_size
     ff(:)=rho_updn(:,ilm,2)
     rho_updn(:,ilm,2)=rho_updn(:,ilm,1)-ff(:)
     rho_updn(:,ilm,1)=ff(:)
   end do
   LIBPAW_DEALLOCATE(ff)
 end if

!Direct links to rho_up and rho_dn
 rho_up => rho_updn(:,:,1)
 rho_dn => rho_updn(:,:,nspden_updn)

!rhoSPH contains the spherical part of effective density
!(including Y00 spherical harmonic)
!-----------------------------------------------------------------
 LIBPAW_ALLOCATE(rhosph,(nrad,nspden_updn))

!  Non-magnetic system: rhoSPH(;,1)=(1/2).rhoSPH_total
 if (nspden==1) then
   rhosph(:,1)=rho_updn(:,1,1)*invsqfpi2

!  Collinear magnetism: rhoSPH = (rhoSPH_up, rhoSPH_dn)
 else if (nspden==2) then
   rhosph(:,1:2)=rho_updn(:,1,1:2)*invsqfpi

!  Non-collinear magnetism: rhoSPH = (rhoSPH_up, rhoSPH_dn)
!    obtained by rotating rho_updn
 else if (nspden==4) then
   LIBPAW_ALLOCATE(m_norm_inv,(nrad))
   LIBPAW_ALLOCATE(rho_nc,(nrad,nspden))
   do ispden=1,nspden
     rho_nc(1:nrad,ispden)=rho_updn(1:nrad,1,ispden)*invsqfpi
   end do
   mag_nc => rho_nc(:,2:4)
   call pawxc_rotate_mag(rho_nc,rhosph,mag_nc,nrad,mag_norm_out=m_norm_inv)
   do ir=1,nrad
     m_norm_inv(ir)=merge(invsqfpi/m_norm_inv(ir),zero,m_norm_inv(ir)>rho_min)
   end do
 end if

!Make spherical density positive
 call pawxc_mkdenpos_wrapper(iwarn,nrad,nspden_updn,0,rhosph,xc_denpos)

!----------------------------------------------------------------------
!----- Compute Exc(rhoSPH) and Vxc(rhoSPH)
!----------------------------------------------------------------------

 LIBPAW_ALLOCATE(exci,(nrad))
 LIBPAW_ALLOCATE(vxci,(nrad,nspden_updn))
 LIBPAW_ALLOCATE(kxci,(nrad,nkxc))
 call pawxcsph(exci,exexch,ixc,kxci,nkxc,nrad,nspden_updn,pawrad,rhosph,vxci,xclevel)

!----------------------------------------------------------------------
!----- Compute numerical derivatives of Vxc,Kxc (by finite diff. scheme)
!----------------------------------------------------------------------

 if (option/=4.and.option/=5) then
   LIBPAW_ALLOCATE(exc_,(nrad))
   LIBPAW_ALLOCATE(rho_,(nrad,nspden_updn))

   if (nspden_updn==2) rho_(:,2)=rhosph(:,2)

!  Compute Exc, Vxc for rho+delta_rho
   LIBPAW_ALLOCATE(vxc1,(nrad,nspden_updn))
   LIBPAW_ALLOCATE(kxc1,(nrad,nkxc))
   rho_(:,1)=(one+delta)*rhosph(:,1)
   call pawxcsph(exc_,exexch,ixc,kxc1,nkxc,nrad,nspden_updn,pawrad,rho_,vxc1,xclevel)

!  Compute Exc, Vxc for rho-delta_rho
   LIBPAW_ALLOCATE(vxc2,(nrad,nspden_updn))
   LIBPAW_ALLOCATE(kxc2,(nrad,nkxc))
   rho_(:,1)=(one-delta)*rhosph(:,1)
   call pawxcsph(exc_,exexch,ixc,kxc2,nkxc,nrad,nspden_updn,pawrad,rho_,vxc2,xclevel)

!  Additional terms for spin-polarized systems
   if (nspden_updn==2) then
     rho_(:,1)=rhosph(:,1)

!    Compute Exc, Vxc for rho+delta_rho_down
     LIBPAW_ALLOCATE(vxcdn1,(nrad,nspden_updn))
     LIBPAW_ALLOCATE(kxcdn1,(nrad,nkxc))
     rho_(:,2)=(one+delta)*rhosph(:,2)
     call pawxcsph(exc_,exexch,ixc,kxcdn1,nkxc,nrad,nspden_updn,pawrad,rho_,vxcdn1,xclevel)

!    Compute Exc, Vxc for rho-delta_rho_down
     LIBPAW_ALLOCATE(vxcdn2,(nrad,nspden_updn))
     LIBPAW_ALLOCATE(kxcdn2,(nrad,nkxc))
     rho_(:,2)=(one-delta)*rhosph(:,2)
     call pawxcsph(exc_,exexch,ixc,kxcdn2,nkxc,nrad,nspden_updn,pawrad,rho_,vxcdn2,xclevel)

   end if !nspden_updn==2
   LIBPAW_DEALLOCATE(exc_)
   LIBPAW_DEALLOCATE(rho_)

!  Store inverse of density finite step
   LIBPAW_ALLOCATE(rhoinv,(nrad,nspden_updn))
   fact=one/delta;if (nspden_updn==1) fact=half*fact
   do ispden=1,nspden_updn
     do ir=1,nrad
       if (rhosph(ir,ispden)>rho_min) then
         rhoinv(ir,ispden)=fact/rhosph(ir,ispden)
       else
         rhoinv(ir,ispden)=zero
       end if
     end do
   end do

!  Compute numerical first derivatives of Vxc (by finite difference scheme)
   LIBPAW_ALLOCATE(d1vxc,(nrad,2*nspden_updn-1))
!  Non-magnetic system: compute dVxc/dn
   if (nspden==1) then
     d1vxc(1:nrad,1)=(vxc1(1:nrad,1)-vxc2(1:nrad,1))*half*rhoinv(1:nrad,1)
!    Collinear magnetism: compute dVxc_up/dn_up,dVxc_dn/dn_up,dVxc_dn/dn_dn
   else if (nspden==2) then
     d1vxc(1:nrad,1)=(vxc1(1:nrad,1)-vxc2(1:nrad,1))*half*rhoinv(1:nrad,1)
     d1vxc(1:nrad,2)=(vxc1(1:nrad,2)-vxc2(1:nrad,2))*half*rhoinv(1:nrad,1)
     d1vxc(1:nrad,3)=(vxcdn1(1:nrad,2)-vxcdn2(1:nrad,2))*half*rhoinv(1:nrad,2)
!    Non-collinear magnetism: compute 1/2 d(Vxc_up+Vxc_dn)/dn,1/2 d(Vxc_up-Vxc_dn)/dn
!    1/2 d(Vxc_up-Vxc_dn)/dm
   else if (nspden==4) then
     do ir=1,nrad
       fact=half*rhoinv(ir,1)
       dvxc1=(vxc1  (ir,1)-vxc2  (ir,1))*fact !dVxc_up/dn_up
       dvxc2=(vxc1  (ir,2)-vxc2  (ir,2))*fact !dVxc_dn/dn_up
       fact=half*rhoinv(ir,2)
       dvxc3=(vxcdn1(ir,2)-vxcdn2(ir,2))*fact !dVxc_dn/dn_dn
       dvxca=dvxc1+dvxc3;dvxcb=dvxc1-dvxc3;dvxcc=two*dvxc2 !Temporary terms
       d1vxc(ir,1)=quarter*(dvxca+dvxcc)  ! 1/2 d(Vxc_up+Vxc_dn)/dn
       d1vxc(ir,2)=quarter* dvxcb         ! 1/2 d(Vxc_up-Vxc_dn)/dn
       d1vxc(ir,3)=quarter*(dvxca-dvxcc)  ! 1/2 d(Vxc_up-Vxc_dn)/dm
     end do
   end if

!  Compute numerical second derivatives of Vxc (by finite difference scheme)
   if (option/=3.or.pawxcdev>=2) then
     LIBPAW_ALLOCATE(d2vxc,(nrad,3*nspden_updn-2))
!    Non-magnetic system: compute d2Vxc/dn2
     if (nspden==1) then
       d2vxc(1:nrad,1)=(vxc1(1:nrad,1)+vxc2(1:nrad,1)-two*vxci(1:nrad,1))*rhoinv(1:nrad,1)**2
!      Collinear magnetism: compute d2Vxc_up/dn_up2,d2Vxc_dn/dn_up2,d2Vxc_up/dn_dn2,d2Vxc_dn/dn_dn2
     else if (nspden==2) then
       d2vxc(1:nrad,1)=(vxc1(1:nrad,1)+vxc2(1:nrad,1)-two*vxci(1:nrad,1))*rhoinv(1:nrad,1)**2
       d2vxc(1:nrad,2)=(vxc1(1:nrad,2)+vxc2(1:nrad,2)-two*vxci(1:nrad,2))*rhoinv(1:nrad,1)**2
       d2vxc(1:nrad,3)=(vxcdn1(1:nrad,1)+vxcdn2(1:nrad,1)-two*vxci(1:nrad,1))*rhoinv(1:nrad,2)**2
       d2vxc(1:nrad,4)=(vxcdn1(1:nrad,2)+vxcdn2(1:nrad,2)-two*vxci(1:nrad,2))*rhoinv(1:nrad,2)**2
!      Non-collinear magnetism: compute 1/2 d2(Vxc_up+Vxc_dn)/dn2,1/2 d2(Vxc_up-Vxc_dn)/dn2
!      1/2 d2(Vxc_up+Vxc_dn)/dm2,1/2 d2(Vxc_up-Vxc_dn)/dm2
     else if (nspden==4) then
       do ir=1,nrad
         fact=rhoinv(ir,1)**2
         dvxc1=(vxc1  (ir,1)+vxc2  (ir,1)-two*vxci(ir,1))*fact !d2Vxc_up/dn_up2
         dvxc2=(vxc1  (ir,2)+vxc2  (ir,2)-two*vxci(ir,2))*fact !d2Vxc_dn/dn_up2
         fact=rhoinv(ir,2)**2
         dvxc3=(vxcdn1(ir,1)+vxcdn2(ir,1)-two*vxci(ir,1))*fact !d2Vxc_up/dn_dn2
         dvxc4=(vxcdn1(ir,2)+vxcdn2(ir,2)-two*vxci(ir,2))*fact !d2Vxc_dn/dn_dn2
         dvxca=dvxc1+dvxc4;dvxcb=dvxc1-dvxc4 !Temporary terms
         dvxcc=dvxc2+dvxc3;dvxcd=dvxc2-dvxc3 !Temporary terms
         d2vxc(ir,1)=(dvxca+three*dvxcc)/8._dp  ! 1/2 d2(Vxc_up+Vxc_dn)/dn2
         d2vxc(ir,2)=(dvxcb+dvxcd)/8._dp        ! 1/2 d2(Vxc_up-Vxc_dn)/dn2
         d2vxc(ir,3)=(dvxca-dvxcc)/8._dp        ! 1/2 d2(Vxc_up+Vxc_dn)/dm2
         d2vxc(ir,4)=(dvxcb-three*dvxcd)/8._dp  ! 1/2 d2(Vxc_up-Vxc_dn)/dm2
       end do
     end if
   end if

!  Compute numerical first and second derivatives of Kxc (by finite difference scheme)
   if (nkxc>0) then
!    Non-magnetic system: compute dKxc/dn, d2Kxc/dn2
     if (nspden==1) then
       LIBPAW_ALLOCATE(d1kxc,(nrad,1))
       LIBPAW_ALLOCATE(d2kxc,(nrad,1))
       d1kxc(1:nrad,1)=(kxc1(1:nrad,1)-kxc2(1:nrad,1))*half*rhoinv(1:nrad,1)
       d2kxc(1:nrad,1)=(kxc1(1:nrad,1)+kxc2(1:nrad,1)-two*kxci(1:nrad,1))*rhoinv(1:nrad,1)**2
!      Collinear magnetism: compute dKxc_upup/dn_up,dKxc_updn/dn_up,dKxc_updn/dn_dn,dKxc_dndn/dn_dn
!      compute d2Kxc_upup/dn_up2,d2Kxc_updn/dn_up2,d2Kxc_upup/dn_dn2,d2Kxc_updn/dn_dn2,d2Kxc_dndn/dn_dn2
     else if (nspden==2) then
       LIBPAW_ALLOCATE(d1kxc,(nrad,4))
       LIBPAW_ALLOCATE(d2kxc,(nrad,5))
       d1kxc(1:nrad,1)=(kxc1(1:nrad,1)-kxc2(1:nrad,1))*half*rhoinv(1:nrad,1)     ! dKxc_upup/dn_up
       d1kxc(1:nrad,2)=(kxc1(1:nrad,2)-kxc2(1:nrad,2))*half*rhoinv(1:nrad,1)     ! dKxc_updn/dn_up
       d1kxc(1:nrad,3)=(kxc1(1:nrad,3)-kxc2(1:nrad,3))*half*rhoinv(1:nrad,1)     ! dKxc_dndn/dn_up
       d1kxc(1:nrad,4)=(kxcdn1(1:nrad,3)-kxcdn2(1:nrad,3))*half*rhoinv(1:nrad,2) ! dKxc_dndn/dn_dn
       d2kxc(1:nrad,1)=(kxc1(1:nrad,1)+kxc2(1:nrad,1)-two*kxci(1:nrad,1))*rhoinv(1:nrad,1)**2      ! d2Kxc_upup/dn_up2
       d2kxc(1:nrad,2)=(kxc1(1:nrad,2)+kxc2(1:nrad,2)-two*kxci(1:nrad,2))*rhoinv(1:nrad,1)**2      ! d2Kxc_updn/dn_up2
       d2kxc(1:nrad,3)=(kxcdn1(1:nrad,1)+kxcdn2(1:nrad,1)-two*kxci(1:nrad,1))*rhoinv(1:nrad,2)**2  ! d2Kxc_upup/dn_dn2
       d2kxc(1:nrad,4)=(kxcdn1(1:nrad,2)+kxcdn2(1:nrad,2)-two*kxci(1:nrad,2))*rhoinv(1:nrad,2)**2  ! d2Kxc_updn/dn_dn2
       d2kxc(1:nrad,5)=(kxcdn1(1:nrad,3)+kxcdn2(1:nrad,3)-two*kxci(1:nrad,3))*rhoinv(1:nrad,2)**2  ! d2Kxc_dndn/dn_dn2
     end if
   end if

   LIBPAW_DEALLOCATE(rhoinv)
   LIBPAW_DEALLOCATE(vxc1)
   LIBPAW_DEALLOCATE(vxc2)
   LIBPAW_DEALLOCATE(kxc1)
   LIBPAW_DEALLOCATE(kxc2)
   if (nspden_updn==2) then
     LIBPAW_DEALLOCATE(vxcdn1)
     LIBPAW_DEALLOCATE(vxcdn2)
     LIBPAW_DEALLOCATE(kxcdn1)
     LIBPAW_DEALLOCATE(kxcdn2)
   end if

 end if ! (option/=4 and option/=5)

 LIBPAW_DEALLOCATE(rhosph)

!If non-collinear magnetism, store 1/2(Vxc_up+Vxc_dn) and 1/2(Vxc_up-Vxc_dn)
 if (nspden==4) then
   vxci(:,1)=half*(vxci(:,1)+vxci(:,2))
   vxci(:,2)=vxci(:,1)-vxci(:,2)
 end if

!----------------------------------------------------------------------
!----- Compute useful sums of densities
!----------------------------------------------------------------------

 if (option/=4.and.option/=5) then

!  Non-collinear magnetism: replace rho_dn by (m_0.dot.m)/|m_0|
   if (nspden==4) then
     LIBPAW_POINTER_ALLOCATE(rho_dn,(nrad,lm_size))
     rho_dn(:,1)=zero
     do ilm=2,lm_size
       rho_dn(1:nrad,ilm)=m_norm_inv(1:nrad) &
&        *(rho_updn(1:nrad,1,2)*rho_updn(1:nrad,ilm,2) &
&         +rho_updn(1:nrad,1,3)*rho_updn(1:nrad,ilm,3) &
&         +rho_updn(1:nrad,1,4)*rho_updn(1:nrad,ilm,4))
     end do
   end if

!  Non-magnetic system:
!  Compute
!  V1SUM1(r)=Sum_L{n_L(r)^2}
!  V2SUM1(r,L)=Sum_L1_L2{n_L1(r)*n_L2(r)*Gaunt_(L,L1,L2)}
!  Collinear magnetism:
!  Compute
!  V1SUM1(r)=Sum_L{n^up_L(r)^2}
!  V1SUM2(r)=Sum_L{n^up_L(r)*n^dn_L(r)}
!  V1SUM3(r)=Sum_L{n^dn_L(r)^2}
!  V2SUM1(r,L)=Sum_L1_L2{n^up_L1(r)*n^up_L2(r)*Gaunt_(L,L1,L2)}
!  V2SUM2(r,L)=Sum_L1_L2{n^up_L1(r)*n^dn_L2(r)*Gaunt_(L,L1,L2)}
!  V2SUM3(r,L)=Sum_L1_L2{n^dn_L1(r)*n^dn_L2(r)*Gaunt_(L,L1,L2)}
!  Non-collinear magnetism:
!  Compute
!  V1SUM1(r)=Sum_L{n_L(r)^2}
!  V1SUM2(r)=Sum_L{n_L(r) (m_0.m_L)}/|m_0|
!  V1SUM3(r)=Sum_L{(m_0.m_L)^2}/|m_0|^2
!  V2SUM1(r,L)=Sum_L1_L2{n_L1(r)*n_L2(r)*Gaunt_(L,L1,L2)}
!  V2SUM2(r,L)=Sum_L1_L2{n_L1(r) (m_0.m_L2)*Gaunt_(L,L1,L2)}/|m_0|
!  V2SUM3(r,L)=Sum_L1_L2{(m_0.m_L1)*(m_0.m_L2)*Gaunt_(L,L1,L2)}/|m_0|^2
   if (pawxcdev>=1)  then
     LIBPAW_ALLOCATE(v1sum,(nrad,nsums))
   else
     LIBPAW_ALLOCATE(v1sum,(0,0))
   end if
   if (pawxcdev>=2)  then
     LIBPAW_ALLOCATE(v2sum,(nrad,lm_size,nsums))
   else
     LIBPAW_ALLOCATE(v2sum,(0,0,0))
   end if
   call pawxcsum(1,1,1,lmselect,lmselect,lm_size,nrad,nsums,pawxcdev,pawang,&
&                rho_up,rho_dn,v1sum,v2sum)

 end if !option

!----------------------------------------------------------------------
!----- Accumulate and store XC potential
!----------------------------------------------------------------------

 if (option/=3.and.option/=4) then

!  === First order development
!  ---------------------------
   if (pawxcdev>=1) then

!    Non-magnetic system
     if (nspden_updn==1) then
       vxc(1:nrad,1,1)=vxci(1:nrad,1)*sqfpi
       if (option/=5) then
         vxc(1:nrad,1,1)=vxc(1:nrad,1,1)+v1sum(1:nrad,1)*d2vxc(1:nrad,1)*invsqfpi2
         do ilm=2,lm_size
           if (lmselect(ilm)) then
             vxc(1:nrad,ilm,1)=d1vxc(1:nrad,1)*rho_up(1:nrad,ilm)
           end if
         end do
       end if

!      Magnetic system (including non-collinear magn.)
     else if (nspden_updn==2) then
       vxc(1:nrad,1,1)=vxci(1:nrad,1)*sqfpi
       vxc(1:nrad,1,2)=vxci(1:nrad,2)*sqfpi
       if (option/=5) then
         vxc(1:nrad,1,1)=vxc(1:nrad,1,1)+invsqfpi2*(v1sum(1:nrad,1)*d2vxc(1:nrad,1) &
&         +two*v1sum(1:nrad,2)*d2vxc(1:nrad,2)+v1sum(1:nrad,3)*d2vxc(1:nrad,3))
         vxc(1:nrad,1,2)=vxc(1:nrad,1,2)+invsqfpi2*(v1sum(1:nrad,1)*d2vxc(1:nrad,2) &
&         +two*v1sum(1:nrad,2)*d2vxc(1:nrad,3)+v1sum(1:nrad,3)*d2vxc(1:nrad,4))
         do ilm=2,lm_size
           if (lmselect(ilm)) then
             vxc(1:nrad,ilm,1)=vxc(1:nrad,ilm,1) &
&             +d1vxc(1:nrad,1)*rho_up(1:nrad,ilm)+d1vxc(1:nrad,2)*rho_dn(1:nrad,ilm)
             vxc(1:nrad,ilm,2)=vxc(1:nrad,ilm,2) &
&             +d1vxc(1:nrad,2)*rho_up(1:nrad,ilm)+d1vxc(1:nrad,3)*rho_dn(1:nrad,ilm)
           end if
         end do
       end if
     end if
   end if ! pawxcdev>=1

!  == 2nd order development
!  ---------------------------
   if (pawxcdev>=2.and.option/=5) then

!    Non-magnetic system
     if (nspden_updn==1) then
       do ilm=2,lm_size
         vxc(1:nrad,ilm,1)=vxc(1:nrad,ilm,1)+half*d2vxc(1:nrad,1)*v2sum(1:nrad,ilm,1)
       end do

!      Magnetic system  (including non-collinear magn.)
     else if (nspden_updn==2) then
       do ilm=2,lm_size
         vxc(1:nrad,ilm,1)=vxc(1:nrad,ilm,1)+d2vxc(1:nrad,2)*v2sum(1:nrad,ilm,2) &
&         +half*(d2vxc(1:nrad,1)*v2sum(1:nrad,ilm,1)+d2vxc(1:nrad,3)*v2sum(1:nrad,ilm,3))
         vxc(1:nrad,ilm,2)=vxc(1:nrad,ilm,2)+d2vxc(1:nrad,3)*v2sum(1:nrad,ilm,2) &
&         +half*(d2vxc(1:nrad,2)*v2sum(1:nrad,ilm,1)+d2vxc(1:nrad,4)*v2sum(1:nrad,ilm,3))
       end do
     end if
   end if !pawxcdev=2

!  === Pathological case: if rho(r) is negative, interpolate Vxc
!  -------------------------------------------------------------
   if (lmselect(1)) then
     tol_rho=xc_denpos*(one+tol6)
     do ispden=1,nspden_updn
       ir1=0;ir2=0
       do ir=1,nrad
         if (rho_updn(ir,1,ispden)<tol_rho) then
           if (ir1==0) ir1=ir-1
           ir2=ir+1
         else if (ir1>0) then
           if (ir1>1.or.ir2<nrad) then
             fact=(vxc(ir2,1,ispden)-vxc(ir1,1,ispden))/(pawrad%rad(ir2)-pawrad%rad(ir1))
             do jr=ir1+1,ir2-1
               vxc(jr,1,ispden)=vxc(ir1,1,ispden)+fact*(pawrad%rad(jr)-pawrad%rad(ir1))
             end do
           end if
           ir1=0;ir2=0
         end if
       end do
     end do
   end if

!  === Non-collinear magnetism: "rotate" back the XC potential
!  ------- ---------------------------------------------------
   if (nspden==4) then
     LIBPAW_ALLOCATE(vxc_diag,(nrad,nspden_updn))
     LIBPAW_ALLOCATE(vxc_nc,(nrad,nspden))
     do ilm=1,lm_size
       vxc_diag(:,1)=vxc(:,ilm,1)+vxc(:,ilm,2) ! Get V from (V_up+V_dn)/2
       vxc_diag(:,2)=vxc(:,ilm,1)-vxc(:,ilm,2) !        and (V_up-V_dn)/2
       call pawxc_rotate_back_mag(vxc_diag,vxc_nc,mag_nc,nrad)
       do ispden=1,nspden
         vxc(1:nrad,ilm,ispden)=vxc_nc(1:nrad,ispden)
       end do
     end do
     LIBPAW_DEALLOCATE(vxc_diag)
     LIBPAW_DEALLOCATE(vxc_nc)
   end if
 end if !option/=3 and option/=4

!----------------------------------------------------------------------
!----- Accumulate and store XC kernel
!----------------------------------------------------------------------

 if (nkxc>0) then

!  === First order development
!  ---------------------------
   if (pawxcdev>=1) then
!    Non-magnetic system:
     if (nspden_updn==1) then
       kxc(1:nrad,1,1)=kxci(1:nrad,1)*sqfpi
       if (option/=5.and.option/=4) then
         kxc(1:nrad,1,1)=kxc(1:nrad,1,1)+invsqfpi2*v1sum(1:nrad,1)*d2kxc(1:nrad,1)
         do ilm=2,lm_size
           if (lmselect(ilm)) then
             kxc(1:nrad,ilm,1)=d1kxc(1:nrad,1)*rho_up(1:nrad,ilm)
           end if
         end do
       end if
!      Magnetic system:
     else if (nspden==2) then
       kxc(1:nrad,1,1)=kxci(1:nrad,1)*sqfpi
       kxc(1:nrad,1,2)=kxci(1:nrad,2)*sqfpi
       kxc(1:nrad,1,3)=kxci(1:nrad,3)*sqfpi
       if (option/=5.and.option/=4) then
         kxc(1:nrad,1,1)=kxc(1:nrad,1,1)+invsqfpi2*(v1sum(1:nrad,1)*d2kxc(1:nrad,1) &
&         +two*v1sum(1:nrad,2)*d2kxc(1:nrad,2)+v1sum(1:nrad,3)*d2kxc(1:nrad,3))
         kxc(1:nrad,1,2)=kxc(1:nrad,1,2)+invsqfpi2*(v1sum(1:nrad,1)*d2kxc(1:nrad,2) &
&         +two*v1sum(1:nrad,2)*d2kxc(1:nrad,3)+v1sum(1:nrad,3)*d2kxc(1:nrad,4))
         kxc(1:nrad,1,3)=kxc(1:nrad,1,3)+invsqfpi2*(v1sum(1:nrad,1)*d2kxc(1:nrad,3) &
&         +two*v1sum(1:nrad,2)*d2kxc(1:nrad,4)+v1sum(1:nrad,3)*d2kxc(1:nrad,5))
         do ilm=2,lm_size
           if (lmselect(ilm)) then
             kxc(1:nrad,ilm,1)=kxc(1:nrad,ilm,1) &
&             +d1kxc(1:nrad,1)*rho_up(1:nrad,ilm)+d1kxc(1:nrad,2)*rho_dn(1:nrad,ilm)
             kxc(1:nrad,ilm,2)=kxc(1:nrad,ilm,2) &
&             +d1kxc(1:nrad,2)*rho_up(1:nrad,ilm)+d1kxc(1:nrad,3)*rho_dn(1:nrad,ilm)
             kxc(1:nrad,ilm,3)=kxc(1:nrad,ilm,3) &
&             +d1kxc(1:nrad,3)*rho_up(1:nrad,ilm)+d1kxc(1:nrad,4)*rho_dn(1:nrad,ilm)
           end if
         end do
       end if
     end if
   end if ! pawxcdev>=1

!  == 2nd order development
!  ---------------------------
   if (pawxcdev>=2.and.option/=4.and.option/=5) then

!    Non-magnetic system:
     if (nspden_updn==1) then
       do ilm=2,lm_size
         kxc(1:nrad,ilm,1)=kxc(1:nrad,ilm,1)+half*d2kxc(1:nrad,1)*v2sum(1:nrad,ilm,1)
       end do
!      Magnetic system:
     else if (nspden==2) then
       do ilm=2,lm_size
         kxc(1:nrad,ilm,1)=kxc(1:nrad,ilm,1)+d2kxc(1:nrad,2)*v2sum(1:nrad,ilm,2) &
&         +half*(d2kxc(1:nrad,1)*v2sum(1:nrad,ilm,1)+d2kxc(1:nrad,3)*v2sum(1:nrad,ilm,3))
         kxc(1:nrad,ilm,2)=kxc(1:nrad,ilm,2)+d2kxc(1:nrad,3)*v2sum(1:nrad,ilm,2) &
&         +half*(d2kxc(1:nrad,2)*v2sum(1:nrad,ilm,1)+d2kxc(1:nrad,4)*v2sum(1:nrad,ilm,3))
         kxc(1:nrad,ilm,3)=kxc(1:nrad,ilm,3)+d2kxc(1:nrad,4)*v2sum(1:nrad,ilm,2) &
&         +half*(d2kxc(1:nrad,3)*v2sum(1:nrad,ilm,1)+d2kxc(1:nrad,5)*v2sum(1:nrad,ilm,3))
       end do
     end if
   end if !pawxcdev=2

!  === Pathological case: if rho(r) is negative, interpolate Kxc
!  -------------------------------------------------------------

!  NOT OK for spin polarized
   if (lmselect(1)) then
     tol_rho=xc_denpos*(one+tol6)
     do ispden=1,nspden_updn
       ir1=0;ir2=0
       do ir=1,nrad
         if (rho_updn(ir,1,ispden)<tol_rho) then
           if (ir1==0) ir1=ir-1
           ir2=ir+1
         else if (ir1>0) then
           if (ir1>1.or.ir2<nrad) then
             fact=(kxc(ir2,1,ispden)-kxc(ir1,1,ispden))/(pawrad%rad(ir2)-pawrad%rad(ir1))
             do jr=ir1+1,ir2-1
               kxc(jr,1,ispden)=kxc(ir1,1,ispden)+fact*(pawrad%rad(jr)-pawrad%rad(ir1))
             end do
           end if
           ir1=0;ir2=0
         end if
       end do
     end do
   end if

!  Non-collinear magnetism: need to store magnetization in kxc
   if (nkxc==6.or.nkxc==22) then
     do ilm=2,lm_size
       kxc(1:nrad,ilm,nkxc-2)=rho_updn(1:nrad,ilm,2)
       kxc(1:nrad,ilm,nkxc-1)=rho_updn(1:nrad,ilm,3)
       kxc(1:nrad,ilm,nkxc  )=rho_updn(1:nrad,ilm,4)
     end do
   end if

 end if ! nkxc>0

 if (nspden==4)  then
   LIBPAW_DEALLOCATE(rho_nc)
   LIBPAW_DEALLOCATE(m_norm_inv)
 end if

 LIBPAW_DEALLOCATE(kxci)
 if (nkxc>0.and.option/=4.and.option/=5) then
   LIBPAW_DEALLOCATE(d1kxc)
   LIBPAW_DEALLOCATE(d2kxc)
 end if

!----------------------------------------------------------------------
!----- Accumulate and store XC energies
!----------------------------------------------------------------------

!----- Calculate Exc (direct scheme) term
!----------------------------------------
 if (option/=1.and.option/=5) then
   LIBPAW_ALLOCATE(ff,(nrad))

!  Contribution from spherical part of rho
   if (nspden==1.or.nspden==4) then
     ff(1:nrad)=rho_updn(1:nrad,1,1)*exci(1:nrad)*sqfpi
   else if (nspden==2) then
     ff(1:nrad)=(rho_updn(1:nrad,1,1)+rho_updn(1:nrad,1,2))*exci(1:nrad)*sqfpi
   end if

!  Contribution from aspherical part of rho
   if (option/=4) then

!    First order development
     if (pawxcdev>=1) then
       if (nspden_updn==1) then
         ff(1:nrad)=ff(1:nrad)+half*v1sum(1:nrad,1)*d1vxc(1:nrad,1)
       else if (nspden_updn==2) then
         ff(1:nrad)=ff(1:nrad)+v1sum(1:nrad,2)*d1vxc(1:nrad,2) &
&         +half*(v1sum(1:nrad,1)*d1vxc(1:nrad,1)+v1sum(1:nrad,3)*d1vxc(1:nrad,3))
       end if
     end if

!    Second order development
     if (pawxcdev>=2) then
       LIBPAW_ALLOCATE(gg,(nrad))

       gg=zero
       do ilm=2,lm_size
         if (lmselect(ilm)) then
           gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,1)*rho_up(1:nrad,ilm)
         end if
       end do
       ff(1:nrad)=ff(1:nrad)+gg(1:nrad)*d2vxc(1:nrad,1)/6._dp

       if (nspden_updn==2) then ! Spin polarized (including non-coll. magn.)
         gg=zero
         do ilm=2,lm_size
           if (lmselect(ilm)) then
             gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,3)*rho_dn(1:nrad,ilm)
           end if
         end do
         ff(1:nrad)=ff(1:nrad)+gg(1:nrad)*d2vxc(1:nrad,4)/6._dp
         gg=zero
         do ilm=2,lm_size
           if (lmselect(ilm)) then
             gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,2)*rho_up(1:nrad,ilm)
           end if
         end do
         ff(1:nrad)=ff(1:nrad)+half*gg(1:nrad)*d2vxc(1:nrad,2)
         gg=zero
         do ilm=2,lm_size
           if (lmselect(ilm)) then
             gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,3)*rho_up(1:nrad,ilm)
           end if
         end do
         ff(1:nrad)=ff(1:nrad)+half*gg(1:nrad)*d2vxc(1:nrad,3)
       end if
       LIBPAW_DEALLOCATE(gg)
     end if

   end if ! option/=4

   ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
   call simp_gen(enxc,ff,pawrad)
   LIBPAW_DEALLOCATE(ff)
 end if ! option/=1 and option/=5

 LIBPAW_DEALLOCATE(exci)
 LIBPAW_DEALLOCATE(vxci)
 if (nspden==4.and.option/=4.and.option/=5)  then
   LIBPAW_POINTER_DEALLOCATE(rho_dn)
 end if
 if (allocated(v1sum))  then
   LIBPAW_DEALLOCATE(v1sum)
 end if
 if (allocated(v2sum))  then
   LIBPAW_DEALLOCATE(v2sum)
 end if
 if (allocated(d1vxc)) then
   LIBPAW_DEALLOCATE(d1vxc)
 end if
 if (allocated(d2vxc)) then
   LIBPAW_DEALLOCATE(d2vxc)
 end if

!----- Calculate Excdc double counting term
!------------------------------------------
 if (option==0.or.option==2) then

   LIBPAW_ALLOCATE(ff,(nrad))

!  Build appropriate density (without core density)
   rho_updn(:,:,:)=rhor(:,:,:)
   if (usexcnhat>0) rho_updn(:,:,:)=rho_updn(:,:,:)+nhat(:,:,:)
   if (nspden==2) then
     do ilm=1,lm_size
       ff(:)=rho_updn(:,ilm,2)
       rho_updn(:,ilm,2)=rho_updn(:,ilm,1)-ff(:)
       rho_updn(:,ilm,1)=ff(:)
     end do
   end if

   ff(1:nrad)=zero

!  Non magnetic or collinear magnetic system:
   if (nspden/=4) then
     do ispden=1,nspden_updn
       do ilm=1,lm_size
         if (lmselect(ilm)) ff(1:nrad)=ff(1:nrad)+vxc(1:nrad,ilm,ispden)*rho_updn(1:nrad,ilm,ispden)
       end do
     end do
   else
!    Non-collinear magnetic system:
     do ilm=1,lm_size
       if (lmselect(ilm)) then
         do ir=1,nrad
           dvxca=vxc(ir,ilm,1)+vxc(ir,ilm,2);dvxcb=vxc(ir,ilm,1)-vxc(ir,ilm,2)
           ff(ir)=ff(ir)+half*(dvxca*rho_updn(ir,ilm,1)+dvxcb*rho_updn(ir,ilm,4)) &
&           +vxc(ir,ilm,3)*rho_updn(ir,ilm,2)-vxc(ir,ilm,4)*rho_updn(ir,ilm,3)
         end do
       end if
     end do
   end if

   ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
   call simp_gen(enxcdc,ff,pawrad)
   LIBPAW_DEALLOCATE(ff)

 end if ! option

 LIBPAW_DEALLOCATE(rho_updn)

 end subroutine pawxcm
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxcm_dfpt
!! NAME
!! pawxcm_dfpt
!!
!! FUNCTION
!! Compute first-order change of XC potential and contribution to
!! 2nd-order change of XC energy inside a PAW sphere.
!! LDA+GGA - USE A DEVELOPMENT OF THE DENSITY OVER (L,M) MOMENTS
!!
!! INPUTS
!!  corexc1(cplex_den*nrad)=first-order change of core density on radial grid
!!  cplex_den= if 1, 1st-order densities are REAL, if 2, COMPLEX
!!  cplex_vxc= if 1, 1st-order XC potential is complex, if 2, COMPLEX
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nrad,lm_size,nkxc)=GS xc kernel
!!  lm_size=size of density array rhor (see below)
!!  lmselect(lm_size)=select the non-zero LM-moments of input density rhor1
!!  nhat1(cplex_den*nrad,lm_size,nspden)=first-order change of compensation density
!!                                        (total in 1st half and spin-up in 2nd half if nspden=2)
!!  nkxc=second dimension of the kxc array
!!  non_magnetic_xc= if true, handle density/potential as non-magnetic (even if it is)
!!  nrad=size of radial mesh for densities/potentials (might be different from pawrad%mesh_size)
!!  nspden=number of spin-density components
!!  option=0  compute both 2nd-order XC energy and 1st-order potential
!!         1  compute only 1st-order XC potential
!!         2  compute only 2nd-order XC energy, XC potential is temporary computed here
!!         3  compute only 2nd-order XC energy, XC potential is input in vxc1(:)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  rhor1(cplex_den*nrad,lm_size,nspden)=first-order change of density
!!  usecore= 1 if core density has to be used in Exc/Vxc ; 0 otherwise
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in d2Exc only
!!             2 if compensation density (nhat) has to be used in d2Exc and Vxc1
!!  xclevel= XC functional level
!!
!! OUTPUT
!!  == if option=0 or 2 or 3 ==rho1_updn
!!    d2enxc=returned exchange-cor. contribution to 2nd-order XC energy
!!
!! SIDE EFFECTS
!!    vxc1(cplex_vxc*nrad,pawang%angl_size,nspden)=1st-order XC potential
!!      Output if option==0 or 1
!!      Unused if option==2
!!      Input  if option==3
!!
!! PARENTS
!!      pawdenpot,pawdfptenergy
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

 subroutine pawxcm_dfpt(corexc1,cplex_den,cplex_vxc,d2enxc,ixc,kxc,lm_size,lmselect,nhat1,&
&                   nkxc,non_magnetic_xc,nrad,nspden,option,pawang,pawrad,rhor1,usecore,&
&                   usexcnhat,vxc1,xclevel,&
&                   d2enxc_im) ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex_den,cplex_vxc,ixc,lm_size,nkxc,nrad,nspden,option
 integer,intent(in) :: usecore,usexcnhat,xclevel
 logical,intent(in) :: non_magnetic_xc
 real(dp),intent(out) :: d2enxc
 real(dp),intent(out),optional :: d2enxc_im
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size)
 real(dp),intent(in) :: corexc1(cplex_den*nrad)
 real(dp),intent(in) :: kxc(nrad,lm_size,nkxc)
 real(dp),intent(in) :: nhat1(cplex_den*nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: rhor1(cplex_den*nrad,lm_size,nspden)
 real(dp),intent(inout),target :: vxc1(cplex_vxc*nrad,lm_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ilm,iplex,ir,ivxc,jr,kr,nkxc_cur
 logical :: need_impart
 real(dp) :: invsqfpi,ro1i,ro1r,sqfpi,sqfpi2,v1i,v1r,vxcrho
 character(len=500) :: msg
!arrays
 integer,parameter :: ikxc(4)=(/1,2,2,3/),irho(4)=(/1,2,1,2/)
! real(dp) :: tsec(2)
 real(dp),allocatable :: ff(:),gg(:),rho1_updn(:,:,:)
 real(dp),allocatable :: v1sum(:),v2sum(:,:)
 real(dp),pointer :: vxc1_(:,:,:)

!************************************************************************

!NOTE (MT)
!lmselect and lm_size are not necessarily the same for densities, kxc and vxc1
!This is not taken into account for the moment, but has to be programmed...

!----------------------------------------------------------------------
!----- Check options
!----------------------------------------------------------------------

 if(option<0.or.option>3) then
   msg='wrong option!'
   MSG_BUG(msg)
 end if
 if(option/=3) then
   call pawxc_get_nkxc(nkxc_cur,nspden,xclevel)
   if(nkxc/=nkxc_cur) then
     msg='Wrong size for kxc array!'
     MSG_BUG(msg)
   end if
 end if
 if(nspden==4.and.option/=3) then
   msg='nspden=4 not implemented (for vxc)!'
   MSG_ERROR(msg)
 end if
 if (option/=1) then
   if (nrad<pawrad%int_meshsz) then
     msg='When option=0,2,3, nrad must be greater than pawrad%int_meshsz!'
     MSG_BUG(msg)
   end if
 end if

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------

!Arrays dimensions and constants
 need_impart=present(d2enxc_im)
 sqfpi=sqrt(four_pi);sqfpi2=half*sqfpi;invsqfpi=one/sqfpi

!Initializations of outputs
 if (option/=1) then
   d2enxc=zero
   if (need_impart) d2enxc_im=zero
 end if
 if (option<=1) vxc1(:,:,:)=zero

!Special case: no XC applied
 if (ixc==0.or.(nkxc==0.and.option/=3)) then
   msg='Note that no xc is applied (ixc=0). Returning'
   MSG_WARNING(msg)
   return
 end if

!----------------------------------------------------------------------
!----- Build several densities
!----------------------------------------------------------------------

!rho1_updn contains the effective 1st-order density used for XC
!with 1st-order core density and/or 1st-order compensation density eventually included
!-----------------------------------------------------------------
 LIBPAW_ALLOCATE(rho1_updn,(cplex_den*nrad,lm_size,nspden))
 rho1_updn(:,:,:)=rhor1(:,:,:)
 if (usexcnhat==2) rho1_updn(:,:,:)=rho1_updn(:,:,:)+nhat1(:,:,:)
 if (usecore==1) then
   if (nspden==1.or.nspden==4) then
     rho1_updn(:,1,1)=rho1_updn(:,1,1)+sqfpi*corexc1(:)
   else if (nspden==2) then
     rho1_updn(:,1,1)=rho1_updn(:,1,1)+sqfpi*corexc1(:)
     rho1_updn(:,1,2)=rho1_updn(:,1,2)+sqfpi2*corexc1(:)
   end if
 end if

!Optionally suppressed magnetic part
 if(non_magnetic_xc) then
   if(nspden==2) rho1_updn(:,:,2)=rho1_updn(:,:,1)*half
   if(nspden==4) rho1_updn(:,:,2:4)=zero
 endif

!In case of collinear magnetism, separate up and down contributions
 if (nspden==2) then
   LIBPAW_ALLOCATE(ff,(cplex_den*nrad))
   do ilm=1,lm_size
     ff(:)=rho1_updn(:,ilm,2)
     rho1_updn(:,ilm,2)=rho1_updn(:,ilm,1)-ff(:)
     rho1_updn(:,ilm,1)=ff(:)
   end do
   LIBPAW_DEALLOCATE(ff)
 end if

!
!----------------------------------------------------------------------
!----- Accumulate and store 1st-order change of XC potential
!----------------------------------------------------------------------

 if (option==2) then
   LIBPAW_POINTER_ALLOCATE(vxc1_,(cplex_vxc*nrad,lm_size,nspden))
 else
   vxc1_ => vxc1
 end if

 if (option/=3) then

   vxc1_=zero
   LIBPAW_ALLOCATE(v1sum,(cplex_vxc*nrad))
   LIBPAW_ALLOCATE(v2sum,(cplex_vxc*nrad,lm_size))

   do ii=1,3*nspden-2
     ivxc=1;if (ii>2) ivxc=2

!    === Vxc1 and Rho1 are REAL
     if (cplex_vxc==1.and.cplex_den==1) then  ! cplex_vxc==1 and cplex_den==1
       call pawxcsum(1,1,1,lmselect,lmselect,lm_size,nrad,1,2,pawang,&
&       kxc(:,:,ikxc(ii)),rho1_updn(:,:,irho(ii)),v1sum,v2sum)
       vxc1_(:,1,ivxc)=vxc1_(:,1,ivxc)+invsqfpi*(v1sum(:)+kxc(:,1,ikxc(ii))*rho1_updn(:,1,irho(ii)))
       do ilm=2,lm_size
         vxc1_(:,ilm,ivxc)=vxc1_(:,ilm,ivxc)+v2sum(:,ilm) &
&         +invsqfpi*(kxc(:,ilm,ikxc(ii))*rho1_updn(:,1  ,irho(ii)) &
&         +kxc(:,1  ,ikxc(ii))*rho1_updn(:,ilm,irho(ii)))
       end do

!    === At least one of Vxc1 or Rho1 is COMPLEX
     else
       call pawxcsum(1,cplex_den,cplex_vxc,lmselect,lmselect,lm_size,nrad,1,2,pawang,&
&       kxc(:,:,ikxc(ii)),rho1_updn(:,:,irho(ii)),v1sum,v2sum)
       do ir=1,nrad
         jr=cplex_den*(ir-1);kr=cplex_vxc*(ir-1)
         do iplex=1,1+(cplex_den*cplex_vxc)/4
           jr=jr+1;kr=kr+1
           vxc1_(kr,1,ivxc)=vxc1_(kr,1,ivxc)+invsqfpi*(v1sum(kr)+kxc(ir,1,ikxc(ii))*rho1_updn(jr,1,irho(ii)))
           do ilm=2,lm_size
             vxc1_(kr,ilm,ivxc)=vxc1_(kr,ilm,ivxc)+v2sum(kr,ilm) &
&             +invsqfpi*(kxc(ir,ilm,ikxc(ii))*rho1_updn(jr,1  ,irho(ii)) &
&             +kxc(ir,1  ,ikxc(ii))*rho1_updn(jr,ilm,irho(ii)))
           end do
         end do
       end do

     end if ! cplex_den and vxc_den
   end do ! ii=1,3*nspden-2

   LIBPAW_DEALLOCATE(v1sum)
   LIBPAW_DEALLOCATE(v2sum)

 end if

!----------------------------------------------------------------------
!----- Accumulate and store 2nd-order change of XC energy
!----------------------------------------------------------------------
 if (option/=1) then

   if (.not.non_magnetic_xc) then
!    For usexnhat=1 particular case, add now compensation density
     if (usexcnhat==1) then
       rho1_updn(:,:,1)=rho1_updn(:,:,1)+nhat1(:,:,nspden)
       if (nspden==2) rho1_updn(:,:,2)=rho1_updn(:,:,2)+nhat1(:,:,1)-nhat1(:,:,2)
     end if
   else
!    Has to be magnetic here
     rho1_updn(:,:,:)=rhor1(:,:,:)
     if (usexcnhat>0) rho1_updn(:,:,:)=rho1_updn(:,:,:)+nhat1(:,:,:)
     if (usecore==1) then
       if (nspden==1.or.nspden==4) then
         rho1_updn(:,1,1)=rho1_updn(:,1,1)+sqfpi*corexc1(:)
       else if (nspden==2) then
         rho1_updn(:,1,1)=rho1_updn(:,1,1)+sqfpi*corexc1(:)
         rho1_updn(:,1,2)=rho1_updn(:,1,2)+sqfpi2*corexc1(:)
       end if
     end if
   end if

   LIBPAW_ALLOCATE(ff,(nrad))
   ff=zero
   if (need_impart) then
     LIBPAW_ALLOCATE(gg,(nrad))
     gg=zero
   end if

!  ----- Calculate d2Exc=Int[Vxc^(1)^*(r).n^(1)(r).dr]
   do ii=1,nspden
!    === Vxc1 and Rho1 are REAL
     if (cplex_vxc==1.and.cplex_den==1) then
       do ilm=1,lm_size
         if (lmselect(ilm)) ff(:)=ff(:)+vxc1_(:,ilm,ii)*rho1_updn(:,ilm,ii)
       end do
!      === Vxc1 and Rho1 are COMPLEX
     else if (cplex_vxc==2.and.cplex_den==2) then  ! cplex_vxc==2 and cplex_den==2
       if (.not.need_impart) then      ! Real part only
         do ilm=1,lm_size
           if (lmselect(ilm)) then
             do ir=1,nrad
               jr=2*ir;v1r=vxc1_(jr-1,ilm,ii);v1i=vxc1_(jr,ilm,ii)
               ro1r=rho1_updn(jr-1,ilm,ii);ro1i=rho1_updn(jr,ilm,ii)
               ff(ir)=ff(ir)+v1r*ro1r+v1i*ro1i
             end do
           end if
         end do
       else                            ! Real and imaginary parts
         do ilm=1,lm_size
           if (lmselect(ilm)) then
             do ir=1,nrad
               jr=2*ir;v1r=vxc1_(jr-1,ilm,ii);v1i=vxc1_(jr,ilm,ii)
               ro1r=rho1_updn(jr-1,ilm,ii);ro1i=rho1_updn(jr,ilm,ii)
               ff(ir)=ff(ir)+v1r*ro1r+v1i*ro1i
               gg(ir)=gg(ir)+v1r*ro1i-v1i*ro1r
             end do
           end if
         end do
       end if ! need_impart
!      === Vxc1 and Rho1 are REAL and COMPLEX
     else
       v1i=zero;ro1i=zero
       do ilm=1,lm_size
         if (lmselect(ilm)) then
           do ir=1,nrad
             jr=cplex_vxc*(ir-1)+1;v1r=vxc1_(jr,ilm,ii);;if(cplex_vxc==2)v1i=vxc1_(jr+1,ilm,ii)
             jr=cplex_den*(ir-1)+1;ro1r=rho1_updn(jr,ilm,ii);if(cplex_den==2)ro1i=rho1_updn(jr+1,ilm,ii)
             ff(ir)=ff(ir)+v1r*ro1r+v1i*ro1i
             if (need_impart) gg(ir)=gg(ir)+v1r*ro1i-v1i*ro1r
           end do
         end if
       end do
     end if ! cplex_vxc and cplex_den
   end do ! ii=1,nspden

   ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
   call simp_gen(vxcrho,ff,pawrad)
   d2enxc=d2enxc+vxcrho
   LIBPAW_DEALLOCATE(ff)

   if (need_impart) then
     gg(1:nrad)=gg(1:nrad)*pawrad%rad(1:nrad)**2
     call simp_gen(vxcrho,gg,pawrad)
     d2enxc_im=d2enxc_im+vxcrho
     LIBPAW_DEALLOCATE(gg)
   end if

 end if

 LIBPAW_DEALLOCATE(rho1_updn)
 if (option==2) then
   LIBPAW_POINTER_DEALLOCATE(vxc1_)
 end if

 end subroutine pawxcm_dfpt
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxcmpositron
!! NAME
!! pawxcmpositron
!!
!! FUNCTION
!! Compute electron-positron correlation potential and energies inside a PAW sphere
!! LDA+GGA - USE A DEVELOPMENT OF THE DENSITY OVER (L,M) MOMENTS
!! Driver of XC functionals.
!!
!! INPUTS
!!  calctype=type of electron-positron calculation:
!!           calctype=1 : positron in electronic density
!!           calctype=2 : electrons in positronic density
!!  corexc(nrad)=electron core density on radial grid
!!  ixcpositron=choice of electron-positron XC scheme
!!  lm_size=size of density array rhor (see below)
!!  lmselect   (lm_size)=select the non-zero LM-moments of input density rhor    (see below)
!!  lmselect_ep(lm_size)=select the non-zero LM-moments of input density rhor_ep (see below)
!!  nhat   (nrad,lm_size,nspden)=compensation density corresponding to rhor
!!  nhat_ep(nrad,lm_size,nspden)=compensation density corresponding to rhor_ep
!!  nrad=size of radial mesh for densities/potentials (might be different from pawrad%mesh_size)
!!  nspden=number of spin-density components
!!  option=0 compute both XC energies (direct+double-counting) and potential
!!         1 compute only XC potential
!!         2 compute only XC energies (direct+double-counting)
!!         3 compute only XC energy by direct scheme
!!         4 compute only XC energy by direct scheme for spherical part of the density
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawxcdev=order of Vxc development
!!  posdensity0_limit=True if we are in the zero positron density limit
!!  rhor(nrad,lm_size,nspden)=electron (or positron) density in real space
!!                             (total in 1st half and spin-up in 2nd half if nspden=2)
!!                             Contents depends on calctype value:
!!                             calctype=1: rhor is the positronic density
!!                             calctype=2: rhor is the electronic density
!!  rhor_ep(nrad,lm_size,nspden)=electron (or positron) density in real space
!!                             (total in 1st half and spin-up in 2nd half if nspden=2)
!!                             Contents depends on calctype value:
!!                             calctype=1: rhor_ep is the electronic density
!!                             calctype=2: rhor_ep is the positronic density
!!  usecore= 1 if core density has to be used in Exc/Vxc for the electronic density ; 0 otherwise
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in double counting energy term only
!!             2 if compensation density (nhat) has to be used in Exc/Vxc and double counting energy term
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!
!! OUTPUT
!!  == if option==0, 2, 3, or 4 ==
!!    enxc=returned exchange and correlation energy (hartree)
!!  == if option==0 or 2 ==
!!    enxcdc=returned exchange-cor. contribution to double-counting energy
!!  == if option==0 or 1 ==
!!    vxc(nrad,lm_size,nspden)=xc potential
!!       (spin up in 1st half and spin-down in 2nd half if nspden=2)
!!
!! NOTES
!!
!! PARENTS
!!      pawdenpot
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxcmpositron(calctype,corexc,enxc,enxcdc,ixcpositron,lm_size,lmselect,lmselect_ep,&
&                         nhat,nhat_ep,nrad,nspden,option,pawang,pawrad,pawxcdev,posdensity0_limit,&
&                         rhor,rhor_ep,usecore,usexcnhat,vxc,xc_denpos)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: calctype,ixcpositron,lm_size,nrad,nspden,option,pawxcdev,usecore
 integer,intent(in) :: usexcnhat
 logical,intent(in) :: posdensity0_limit
 real(dp),intent(in) :: xc_denpos
 real(dp),intent(out) :: enxc,enxcdc
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size),lmselect_ep(lm_size)
 real(dp),intent(in) :: corexc(nrad)
 real(dp),intent(in) :: nhat   (nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: nhat_ep(nrad,lm_size,nspden*((usexcnhat+1)/2))
 real(dp),intent(in) :: rhor   (nrad,lm_size,nspden)
 real(dp),intent(in) :: rhor_ep(nrad,lm_size,nspden)
 real(dp),intent(out) :: vxc(nrad,lm_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ilm,ir,ir1,ir2,iwarn,iwarnp,jr
 real(dp),parameter :: delta=1.d-4
 real(dp) :: fact,invsqfpi,sqfpi,rhomin
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: d1vxc(:,:),d2vxc(:,:),fxc_(:),ff(:),fxci(:),gg(:)
 real(dp),allocatable :: rho_(:),rhotot(:,:),rhotot_ep(:,:),rhoinv(:),rhoinv_ep(:)
 real(dp),allocatable :: rhosph(:),rhosph_ep(:),v1sum(:,:),v2sum(:,:,:)
 real(dp),allocatable :: vxce1(:),vxce1_ep(:),vxce2(:),vxce2_ep(:)
 real(dp),allocatable :: vxcp1(:),vxcp1_ep(:),vxcp2(:),vxcp2_ep(:)
 real(dp),allocatable :: vxcei(:),vxcpi(:)

!************************************************************************

!----- Check options
 if(calctype/=1.and.calctype/=2) then
   msg='Invalid value for calctype'
   MSG_BUG(msg)
 end if
 if (option/=1) then
   if (nrad<pawrad%int_meshsz) then
     msg='When option=0,2,3,4, nrad must be greater than pawrad%int_meshsz!'
     MSG_BUG(msg)
   end if
 end if

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------

!Initializations and constants
 iwarn=0;iwarnp=1
 sqfpi=sqrt(four_pi)
 invsqfpi=one/sqfpi

!Initializations of output arrays
 if (option/=1) enxc=zero
 if (option==0.or.option==2) enxcdc=zero
 if (option<3) vxc(:,:,:)=zero

 if (ixcpositron==0) then ! No xc at all is applied (usually for testing)
   msg='Note that no xc is applied (ixc=0). Returning'
   MSG_WARNING(msg)
   return
 end if

!----------------------------------------------------------------------
!----- Build several densities
!----------------------------------------------------------------------

!rhotot/rhotot_ep contain the effective total densities used for XC
!with core density and/or compensation density eventually included
!-----------------------------------------------------------------
!Input density
 LIBPAW_ALLOCATE(rhotot,(nrad,lm_size))
 LIBPAW_ALLOCATE(rhotot_ep,(nrad,lm_size))
 rhotot   (:,:)=rhor   (:,:,1)
 rhotot_ep(:,:)=rhor_ep(:,:,1)
!Eventually add compensation density
 if (usexcnhat==2) then
   rhotot   (:,:)=rhotot   (:,:)+nhat   (:,:,1)
   rhotot_ep(:,:)=rhotot_ep(:,:)+nhat_ep(:,:,1)
 end if
!Eventually add core density
 if (usecore==1) then
   if (calctype==1) rhotot_ep(:,1)=rhotot_ep(:,1)+sqfpi*corexc(:)
   if (calctype==2) rhotot   (:,1)=rhotot   (:,1)+sqfpi*corexc(:)
 end if

!rhoSPH/rhoSPH_ep contain the spherical part of effective densities
!(including Y00 spherical harmonic)
!-----------------------------------------------------------------
 LIBPAW_ALLOCATE(rhosph,(nrad))
 LIBPAW_ALLOCATE(rhosph_ep,(nrad))

 rhosph   (:)=rhotot   (:,1)*invsqfpi
 rhosph_ep(:)=rhotot_ep(:,1)*invsqfpi

!Make spherical densities positive
 if (calctype==1) then
   if (.not.posdensity0_limit) then
     call pawxc_mkdenpos_wrapper(iwarnp,nrad,1,1,rhosph,xc_denpos)
   end if
   call pawxc_mkdenpos_wrapper(iwarn ,nrad,1,1,rhosph_ep,xc_denpos)
 else if (calctype==2) then
   call pawxc_mkdenpos_wrapper(iwarn ,nrad,1,1,rhosph,xc_denpos)
   if (.not.posdensity0_limit) then
     call pawxc_mkdenpos_wrapper(iwarnp,nrad,1,1,rhosph_ep,xc_denpos)
   end if
 end if

!----------------------------------------------------------------------
!----- Compute Exc(rhoSPH,rhoSPH_ep) and Vxc(rhoSPH,rhoSPH_ep)
!----------------------------------------------------------------------

 LIBPAW_ALLOCATE(fxci,(nrad))
 LIBPAW_ALLOCATE(vxcei,(nrad))
 LIBPAW_ALLOCATE(vxcpi,(nrad))
 call pawxcsphpositron(calctype,fxci,ixcpositron,nrad,pawrad,posdensity0_limit,rhosph,rhosph_ep,vxcei,vxcpi)

!----------------------------------------------------------------------
!----- Compute numerical derivatives of Vxc (by finite diff. scheme)
!----------------------------------------------------------------------

 if (option/=4) then

   LIBPAW_ALLOCATE(fxc_,(nrad))
   LIBPAW_ALLOCATE(rho_,(nrad))

!  Compute Vxc for (rho+delta_rho,rho_ep)
   LIBPAW_ALLOCATE(vxce1,(nrad))
   LIBPAW_ALLOCATE(vxcp1,(nrad))
   rho_(:)=(one+delta)*rhosph(:)
   call pawxcsphpositron(calctype,fxc_,ixcpositron,nrad,pawrad,posdensity0_limit,rho_,rhosph_ep,vxce1,vxcp1)

!  Compute Vxc for(rho-delta_rho,rho_ep)
   LIBPAW_ALLOCATE(vxce2,(nrad))
   LIBPAW_ALLOCATE(vxcp2,(nrad))
   rho_(:)=(one-delta)*rhosph(:)
   call pawxcsphpositron(calctype,fxc_,ixcpositron,nrad,pawrad,posdensity0_limit,rho_,rhosph_ep,vxce2,vxcp2)

!  Compute Vxc for (rho,rho_ep+delta_rho_ep)
   LIBPAW_ALLOCATE(vxce1_ep,(nrad))
   LIBPAW_ALLOCATE(vxcp1_ep,(nrad))
   rho_(:)=(one+delta)*rhosph_ep(:)
   call pawxcsphpositron(calctype,fxc_,ixcpositron,nrad,pawrad,posdensity0_limit,rhosph,rho_,vxce1_ep,vxcp1_ep)

!  Compute Vxc for (rho,rho_ep-delta_rho_ep)
   LIBPAW_ALLOCATE(vxce2_ep,(nrad))
   LIBPAW_ALLOCATE(vxcp2_ep,(nrad))
   rho_(:)=(one-delta)*rhosph_ep(:)
   call pawxcsphpositron(calctype,fxc_,ixcpositron,nrad,pawrad,posdensity0_limit,rhosph,rho_,vxce2_ep,vxcp2_ep)

   LIBPAW_DEALLOCATE(fxc_)
   LIBPAW_DEALLOCATE(rho_)

!  Store inverse of density finite step
   LIBPAW_ALLOCATE(rhoinv,(nrad))
   LIBPAW_ALLOCATE(rhoinv_ep,(nrad))
   fact=one/delta
   do ir=1,nrad
     if (rhosph(ir)>rho_min) then
       rhoinv(ir)=fact/rhosph(ir)
     else
       rhoinv(ir)=zero
     end if
     if (rhosph_ep(ir)>rho_min) then
       rhoinv_ep(ir)=fact/rhosph_ep(ir)
     else
       rhoinv_ep(ir)=zero
     end if
   end do

!  Compute numerical first derivatives of Vxc (by finite difference scheme)
   LIBPAW_ALLOCATE(d1vxc,(nrad,3))
   if (calctype==1) then
     d1vxc(:,1)=(vxcp1   (:)-vxcp2   (:))*half*rhoinv   (:)  ! dVxc+/drho+
     d1vxc(:,2)=(vxcp1_ep(:)-vxcp2_ep(:))*half*rhoinv_ep(:)  ! dVxc+/drho-
     d1vxc(:,3)=(vxce1_ep(:)-vxce2_ep(:))*half*rhoinv_ep(:)  ! dVxc-/drho-
   else if (calctype==2) then
     d1vxc(:,1)=(vxce1   (:)-vxce2   (:))*half*rhoinv   (:)  ! dVxc-/drho-
     d1vxc(:,2)=(vxcp1   (:)-vxcp2   (:))*half*rhoinv   (:)  ! dVxc+/drho-
!    d1vxc(:,2)=(vxce1_ep(:)-vxce2_ep(:))*half*rhoinv_ep(:)  ! dVxc-/drho+
     d1vxc(:,3)=(vxcp1_ep(:)-vxcp2_ep(:))*half*rhoinv_ep(:)  ! dVxc+/drho+
   end if

!  Compute numerical second derivatives of Vxc (by finite difference scheme)
   if (option<3.or.pawxcdev>1) then
     LIBPAW_ALLOCATE(d2vxc,(nrad,4))
     if (calctype==1) then
       d2vxc(:,1)=(vxcp1   (:)+vxcp2   (:)-two*vxcpi(:))*rhoinv   (:)**2  ! d2Vxc+/drho+_drho+
       d2vxc(:,2)=(vxce1   (:)+vxce2   (:)-two*vxcei(:))*rhoinv   (:)**2  ! d2Vxc-/drho+_drho+
       d2vxc(:,3)=(vxcp1_ep(:)+vxcp2_ep(:)-two*vxcpi(:))*rhoinv_ep(:)**2  ! d2Vxc+/drho-_drho-
       d2vxc(:,4)=(vxce1_ep(:)+vxce2_ep(:)-two*vxcei(:))*rhoinv_ep(:)**2  ! d2Vxc-/drho-_drho-
     else if (calctype==2) then
       d2vxc(:,1)=(vxce1   (:)+vxce2   (:)-two*vxcei(:))*rhoinv   (:)**2  ! d2Vxc-/drho-_drho-
       d2vxc(:,2)=(vxcp1   (:)+vxcp2   (:)-two*vxcpi(:))*rhoinv   (:)**2  ! d2Vxc+/drho-_drho-
       d2vxc(:,3)=(vxce1_ep(:)+vxce2_ep(:)-two*vxcei(:))*rhoinv_ep(:)**2  ! d2Vxc-/drho+_drho+
       d2vxc(:,4)=(vxcp1_ep(:)+vxcp2_ep(:)-two*vxcpi(:))*rhoinv_ep(:)**2  ! d2Vxc+/drho+_drho+
     end if
   end if ! option

   LIBPAW_DEALLOCATE(rhoinv)
   LIBPAW_DEALLOCATE(rhoinv_ep)
   LIBPAW_DEALLOCATE(vxce1)
   LIBPAW_DEALLOCATE(vxcp1)
   LIBPAW_DEALLOCATE(vxce2)
   LIBPAW_DEALLOCATE(vxcp2)
   LIBPAW_DEALLOCATE(vxce1_ep)
   LIBPAW_DEALLOCATE(vxcp1_ep)
   LIBPAW_DEALLOCATE(vxce2_ep)
   LIBPAW_DEALLOCATE(vxcp2_ep)

 end if ! option/=4

 LIBPAW_DEALLOCATE(rhosph)
 LIBPAW_DEALLOCATE(rhosph_ep)

!----------------------------------------------------------------------
!----- Compute useful sums of densities
!----------------------------------------------------------------------

 if (option<3.or.option/=1) then

!  Compute V1SUM1(r)=Sum_L{n^el_L(r)^2}
!  V1SUM2(r)=Sum_L{n^el_L(r)*n^pos_L(r)}
!  V1SUM3(r)=Sum_L{n^pos_L(r)^2}
!  V2SUM1(r,L)=Sum_L1_L2{n^el_L1(r)*n^el_L2(r)*Gaunt_(L,L1,L2)}
!  V2SUM2(r,L)=Sum_L1_L2{n^el_L1(r)*n^pos_L2(r)*Gaunt_(L,L1,L2)}
!  V2SUM3(r,L)=Sum_L1_L2{n^pos_L1(r)*n^pos_L2(r)*Gaunt_(L,L1,L2)}
   if (pawxcdev>=1)  then
     LIBPAW_ALLOCATE(v1sum,(nrad,3))
   else
     LIBPAW_ALLOCATE(v1sum,(0,0))
   end if
   if (pawxcdev>=2)  then
     LIBPAW_ALLOCATE(v2sum,(nrad,lm_size,3))
   else
     LIBPAW_ALLOCATE(v2sum,(0,0,0))
   end if
   call pawxcsum(1,1,1,lmselect,lmselect_ep,lm_size,nrad,3,pawxcdev,pawang,rhotot,rhotot_ep,v1sum,v2sum)

 end if !option

!----------------------------------------------------------------------
!----- Accumulate and store XC potential
!----------------------------------------------------------------------

 if (option<3) then

!  if (option==0.or.option==2) allocate(vxc_ep(nrad,lm_size))

!  === First order development
!  ---------------------------
   if (pawxcdev>=1) then
     if (calctype==1) vxc(:,1,1)=vxcpi(:)*sqfpi
     if (calctype==2) vxc(:,1,1)=vxcei(:)*sqfpi
     vxc(:,1,1)=vxc(:,1,1)+invsqfpi*(d2vxc(:,2)*v1sum(:,2) &
&     +half*(d2vxc(:,1)*v1sum(:,1)+d2vxc(:,3)*v1sum(:,3)))
     do ilm=2,lm_size
       if (lmselect(ilm))    vxc(:,ilm,1)=vxc(:,ilm,1)+d1vxc(:,1)*rhotot   (:,ilm)
       if (lmselect_ep(ilm)) vxc(:,ilm,1)=vxc(:,ilm,1)+d1vxc(:,2)*rhotot_ep(:,ilm)
     end do
!    if (option==0.or.option==2) then
!    if (calctype==1) vxc_ep(:,1)=vxcei(:)*sqfpi
!    if (calctype==2) vxc_ep(:,1)=vxcpi(:)*sqfpi
!    vxc_ep(:,1)=vxc_ep(:,1,1)+invsqfpi*(d2vxc(:,3)*v1sum(:,2) &
!    &             +half*(d2vxc(:,2)*v1sum(:,1)+d2vxc(:,4)*v1sum(:,3)))
!    do ilm=2,lm_size
!    if (lmselect(ilm))    vxc_ep(:,ilm)=vxc_ep(:,ilm)+d1vxc(:,2)*rhotot   (:,ilm)
!    if (lmselect_ep(ilm)) vxc_ep(:,ilm)=vxc_ep(:,ilm)+d1vxc(:,3)*rhotot_ep(:,ilm)
!    end do
!    end if
   end if ! pawxcdev>=1

!  == 2nd order development
!  ---------------------------
   if (pawxcdev>=2) then
     do ilm=2,lm_size
       vxc(:,ilm,1)=vxc(:,ilm,1)+d2vxc(:,2)*v2sum(:,ilm,2) &
&       +half*(d2vxc(:,1)*v2sum(:,ilm,1)+d2vxc(:,3)*v2sum(:,ilm,3))
     end do
!    if (option==0.or.option==2) then
!    do ilm=2,lm_size
!    vxc_ep(:,ilm)=vxc_ep(:,ilm)+d2vxc(:,3)*v2sum(:,ilm,2) &
!    &                +half*(d2vxc(:,2)*v2sum(:,ilm,1)+d2vxc(:,4)*v2sum(:,ilm,3))
!    end do
!    end if
   end if !pawxcdev=2

!  === Pathological case: if rho(r) is negative, interpolate Vxc
!  -------------------------------------------------------------
   if (lmselect(1)) then
     rhomin=xc_denpos*(one+tol6)
     ir1=0;ir2=0
     do ir=1,nrad
       if (rhotot(ir,1)<rhomin) then
         if (ir1==0) ir1=ir-1
         ir2=ir+1
       else if (ir1>0) then
         if (ir1>1.or.ir2<nrad) then
           fact=(vxc(ir2,1,1)-vxc(ir1,1,1))/(pawrad%rad(ir2)-pawrad%rad(ir1))
           do jr=ir1+1,ir2-1
             vxc(jr,1,1)=vxc(ir1,1,1)+fact*(pawrad%rad(jr)-pawrad%rad(ir1))
           end do
         end if
         ir1=0;ir2=0
       end if
     end do
   end if
!  if (option==0.or.option==2) then
!  if (lmselect_ep(1)) then
!  ir1=0;ir2=0
!  do ir=1,nrad
!  if (rhotot_ep(ir,1)<rho_min) then
!  if (ir1==0) ir1=ir-1
!  ir2=ir+1
!  else if (ir1>0) then
!  if (ir1>1.or.ir2<nrad) then
!  fact=(vxc_ep(ir2,1)-vxc_ep(ir1,1))/(pawrad%rad(ir2)-pawrad%rad(ir1))
!  do jr=ir1+1,ir2-1
!  vxc_ep(jr,1)=vxc_ep(ir1,1)+fact*(pawrad%rad(jr)-pawrad%rad(ir1))
!  end do
!  end if
!  ir1=0;ir2=0
!  end if
!  end do
!  end if
!  end if

!  When vxc is dimensionned as polarized...
   if (nspden>=2) vxc(:,:,2)=vxc(:,:,1)
   if (nspden==4) vxc(:,:,3:4)=zero

 end if !option<3

 LIBPAW_DEALLOCATE(vxcei)
 LIBPAW_DEALLOCATE(vxcpi)

!----------------------------------------------------------------------
!----- Accumulate and store XC energies
!----------------------------------------------------------------------

!----- Calculate Exc (direct scheme) term
!----------------------------------------

 if (option/=1) then
   LIBPAW_ALLOCATE(ff,(nrad))

!  Contribution from spherical part of rho
   ff(:)=fxci(:)*four_pi

!  Contribution from aspherical part of rho
   if (option/=4) then

!    First order development
     if (pawxcdev>=1) then
       ff(:)=ff(:)+v1sum(:,2)*d1vxc(:,2) &
&       +half*(v1sum(:,1)*d1vxc(:,1)+v1sum(:,3)*d1vxc(:,3))
     end if

!    Second order development
     if (pawxcdev>=2) then
       LIBPAW_ALLOCATE(gg,(nrad))
       gg=zero
       do ilm=2,lm_size
         if (lmselect(ilm))    gg(:)=gg(:)+v2sum(:,ilm,1)*rhotot(:,ilm)
       end do
       ff(:)=ff(:)+gg(:)*d2vxc(:,1)/6._dp
       gg=zero
       do ilm=2,lm_size
         if (lmselect(ilm))    gg(:)=gg(:)+v2sum(:,ilm,2)*rhotot(:,ilm)
       end do
       ff(:)=ff(:) +half*gg(:)*d2vxc(:,2)
       gg=zero
       do ilm=2,lm_size
         if (lmselect(ilm))    gg(:)=gg(:)+v2sum(:,ilm,3)*rhotot(:,ilm)
       end do
       ff(:)=ff(:) +half*gg(:)*d2vxc(:,3)
       gg=zero
       do ilm=2,lm_size
         if (lmselect_ep(ilm)) gg(:)=gg(:)+v2sum(:,ilm,3)*rhotot_ep(:,ilm)
       end do
       ff(:)=ff(:)+gg(:)*d2vxc(:,4)/6._dp
       LIBPAW_DEALLOCATE(gg)
     end if ! pawxcdev>=2

   end if ! option/=4

   ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
   call simp_gen(enxc,ff,pawrad)
   LIBPAW_DEALLOCATE(ff)
 end if ! option/=1

 LIBPAW_DEALLOCATE(fxci)
 if (option<3.or.option/=1)  then
   LIBPAW_DEALLOCATE(v1sum)
   LIBPAW_DEALLOCATE(v2sum)
 end if
 if (option<3.or.(option/=4.and.pawxcdev>1))   then
   LIBPAW_DEALLOCATE(d2vxc)
 end if
 if (option/=4)  then
   LIBPAW_DEALLOCATE(d1vxc)
 end if

!----- Calculate Excdc double counting term
!------------------------------------------
 if (option==0.or.option==2) then

!  Build appropriate density
   if (usexcnhat==1) rhotot(:,:)=rhotot(:,:)+nhat(:,:,1)
   if (usecore==1.and.calctype==2) rhotot(:,1)=rhotot(:,1)-sqfpi*corexc(:)

!  Integrate with potential
   LIBPAW_ALLOCATE(ff,(nrad))
   ff(:)=zero
   do ilm=1,lm_size
     if (lmselect(ilm)) ff(:)=ff(:)+vxc(:,ilm,1)*rhotot(:,ilm)
   end do
   ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
   call simp_gen(enxcdc,ff,pawrad)
   LIBPAW_DEALLOCATE(ff)
 end if ! option

 LIBPAW_DEALLOCATE(rhotot)
 LIBPAW_DEALLOCATE(rhotot_ep)

end subroutine pawxcmpositron
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxc_get_nkxc
!! NAME
!! pawxc_get_nkxc
!!
!! FUNCTION
!! Get size of XC kernel array (Kxc) according to spin polarization and XC type
!!
!! INPUTS
!!  nspden= nmber of density spin components
!!  xclevel= XC type
!!
!! OUTPUT
!!  nkxc= size of XC kernel (kxc array)
!!
!! NOTES
!!  Content of Kxc array:
!!   ===== if LDA
!!    if nspden==1: kxc(:,1)= d2Exc/drho2
!!                 (kxc(:,2)= d2Exc/drho_up drho_dn)
!!    if nspden>=2: kxc(:,1)= d2Exc/drho_up drho_up
!!                  kxc(:,2)= d2Exc/drho_up drho_dn
!!                  kxc(:,3)= d2Exc/drho_dn drho_dn
!!    if nspden==4: kxc(:,4:6)= (m_x, m_y, m_z) (magnetization)
!!   ===== if GGA
!!    if nspden==1:
!!       kxc(:,1)= d2Exc/drho2
!!       kxc(:,2)= 1/|grad(rho)| dExc/d|grad(rho)|
!!       kxc(:,3)= 1/|grad(rho)| d2Exc/d|grad(rho)| drho
!!       kxc(:,4)= 1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dExc/d|grad(rho)| )
!!       kxc(:,5)= gradx(rho)
!!       kxc(:,6)= grady(rho)
!!       kxc(:,7)= gradz(rho)
!!    if nspden>=2:
!!       kxc(:,1)= d2Exc/drho_up drho_up
!!       kxc(:,2)= d2Exc/drho_up drho_dn
!!       kxc(:,3)= d2Exc/drho_dn drho_dn
!!       kxc(:,4)= 1/|grad(rho_up)| dEx/d|grad(rho_up)|
!!       kxc(:,5)= 1/|grad(rho_dn)| dEx/d|grad(rho_dn)|
!!       kxc(:,6)= 1/|grad(rho_up)| d2Ex/d|grad(rho_up)| drho_up
!!       kxc(:,7)= 1/|grad(rho_dn)| d2Ex/d|grad(rho_dn)| drho_dn
!!       kxc(:,8)= 1/|grad(rho_up)| * d/d|grad(rho_up)| ( 1/|grad(rho_up)| dEx/d|grad(rho_up)| )
!!       kxc(:,9)= 1/|grad(rho_dn)| * d/d|grad(rho_dn)| ( 1/|grad(rho_dn)| dEx/d|grad(rho_dn)| )
!!       kxc(:,10)=1/|grad(rho)| dEc/d|grad(rho)|
!!       kxc(:,11)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_up
!!       kxc(:,12)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_dn
!!       kxc(:,13)=1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dEc/d|grad(rho)| )
!!       kxc(:,14)=gradx(rho_up)
!!       kxc(:,15)=gradx(rho_dn)
!!       kxc(:,16)=grady(rho_up)
!!       kxc(:,17)=grady(rho_dn)
!!       kxc(:,18)=gradz(rho_up)
!!       kxc(:,19)=gradz(rho_dn)
!!    if nspden==4:
!!       kxc(:,20:22)= (m_x, m_y, m_z) (magnetization)
!!
!! PARENTS
!!      m_pawxc,respfn,nonlinear
!!
!! CHILDREN
!!
!! SOURCE

 subroutine pawxc_get_nkxc(nkxc,nspden,xclevel)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nspden,xclevel
 integer,intent(out) :: nkxc
!arrays

!Local variables-------------------------------
!scalars
!arrays

!************************************************************************

 nkxc=0

 if (nspden==1) then ! Non polarized

   if (xclevel==1) nkxc=1
   if (xclevel==2) nkxc=7

 else if (nspden==2) then ! Polarized

   if (xclevel==1) nkxc=3
   if (xclevel==2) nkxc=19

 else if (nspden==4) then ! Non-collinear

   ! Store magnetization in the 3 last terms of Kxc
   if (xclevel==1) nkxc=6
   if (xclevel==2) nkxc=22

 end if

 end subroutine pawxc_get_nkxc
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxc_drivexc_wrapper
!! NAME
!! pawxc_drivexc_wrapper
!!
!! FUNCTION
!! PAW only
!! Wrapper for drivexc routines
!!
!! NOTES
!! PENDING. Need to manage properly optional arguments:
!! Check that these are present before calling drivexc
!! Probably use better interfaces of fortran 2003 to avoid
!! numerous if/then sentences.
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

 subroutine pawxc_drivexc_wrapper(exc,ixc,mgga,ndvxc,nd2vxc,ngr2,npts,nspden,nvxcgrho,&
&           order,rho,use_laplacian,vxcrho,xclevel, &
&           dvxc,d2vxc,el_temp,exexch,fxcT,grho2,lrho,tau,vxcgrho,vxclrho,vxctau) ! Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,mgga,ndvxc,nd2vxc,ngr2,npts,nspden,nvxcgrho,order,xclevel,use_laplacian
!arrays
 real(dp),intent(in) :: rho(npts,nspden)
 real(dp),intent(out) :: exc(npts),vxcrho(npts,nspden)
 integer,intent(in),optional :: exexch
 real(dp),intent(in),optional :: el_temp
 real(dp),intent(in),optional:: grho2(npts,ngr2),lrho(npts,nspden*use_laplacian),tau(npts,nspden*mgga)
 real(dp),intent(out),optional:: dvxc(npts,ndvxc),d2vxc(npts,nd2vxc),fxcT(npts),vxcgrho(npts,nvxcgrho)
 real(dp),intent(out),optional:: vxclrho(npts,nspden*use_laplacian),vxctau(npts,nspden*mgga)

!Local variables-------------------------------
 character(len=100) :: msg

! *************************************************************************

!One could add here a section for other codes (i.e. BigDFT, ...)
#if defined HAVE_LIBPAW_ABINIT
 call pawxc_drivexc_abinit()
#elif defined HAVE_LIBXC
 call pawxc_drivexc_libxc()
#else
 write(msg,'(5a)') 'libPAW XC driving routine only implemented in the following cases:',ch10, &
&                  ' - ABINIT',ch10,' - libXC'
 MSG_BUG(msg)
#endif

 if (.false.) write(std_out,*) el_temp,lrho(1,1),tau(1,1)
!!***

contains
!!***

#if defined HAVE_LIBPAW_ABINIT
!!****f* m_pawxc/pawxc_drivexc_abinit
!! NAME
!!  pawxc_drivexc_abinit
!!
!! FUNCTION
!!  ABINIT version of XC driving routine
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxc_drivexc_abinit()

! *************************************************************************

 if ((.not.present(dvxc)).or.(.not.present(grho2)).or.(.not.present(vxcgrho))) then
  msg='dvxc, grho2 and vxcgrho should be present in pawxc_drivexc_wrapper'
  MSG_BUG(msg)
end if

!Call to main XC driver
!PENDING: we cannot handle all optional-variable combinations.
!Hence, only two posibilities are considered here:
!1) Pass dvxc, exexch, grho2 and vxcgrho

 if (present(exexch)) then
   call drivexc_main(exc,ixc,mgga,ndvxc,nd2vxc,ngr2,npts,nspden,nvxcgrho,order,rho,use_laplacian,vxcrho,xclevel,&
&   dvxc=dvxc,d2vxc=d2vxc,exexch=exexch,grho2=grho2,lrho=lrho,tau=tau,vxcgrho=vxcgrho,vxclrho=vxclrho,vxctau=vxctau)
 else
!2) Pass only dvxc, grho2 and vxcgrho
   call drivexc_main(exc,ixc,mgga,ndvxc,nd2vxc,ngr2,npts,nspden,nvxcgrho,order,rho,use_laplacian,vxcrho,xclevel,&
&   dvxc=dvxc,d2vxc=d2vxc,grho2=grho2,lrho=lrho,tau=tau,vxcgrho=vxcgrho,vxclrho=vxclrho,vxctau=vxctau)
 end if

end subroutine pawxc_drivexc_abinit
!!***
#endif

#if defined HAVE_LIBXC
!!****f* m_pawxc/pawxc_drivexc_libxc
!! NAME
!!  pawxc_drivexc_libxc
!!
!! FUNCTION
!!  LibXC version of XC driving routine
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

subroutine pawxc_drivexc_libxc()

! *************************************************************************

!Check the compatibility of input arguments
 if (libxc_functionals_ismgga()) then
   msg='MGGA is not yet coded in pawxc_drivexc_wrapper/LIBXC'
   MSG_ERROR(msg)
 end if
 if (libxc_functionals_needs_laplacian()) then
   msg='Laplacian based XC functionals are not yet coded in pawxc_drivexc_wrapper/LIBXC'
   MSG_ERROR(msg)
 end if
 if (ixc>=0) then
   msg='ixc argument should be negative!'
   MSG_BUG(msg)
 end if
 if (ixc/=libxc_functionals_ixc()) then
   msg='The value of ixc differs from the one used to initialize the functional!'
   MSG_BUG(msg)
 end if
 if ((order<1.and.order/=-2).or.order>4) then
   msg='The only allowed values for order are 1, 2, -2, or 3!'
   MSG_BUG(msg)
 end if
 if ((order**2>1).and.(.not.present(dvxc))) then
   msg='The value of order is not compatible with the presence of the array dvxc!'
   MSG_BUG(msg)
 end if
 if ((order==3).and.(.not.present(d2vxc))) then
   msg='The value of order is not compatible with the presence of the array d2vxc!'
   MSG_BUG(msg)
 end if
 if (libxc_functionals_isgga()) then
   if ((.not.present(grho2)).or.(.not.present(vxcgrho)).or.(nvxcgrho==0))  then
     write(msg,'(3a)') 'At least one of the functionals is a GGA,',ch10, &
&      'but not all the necessary optional arguments are present.'
     MSG_BUG(msg)
   end if
   if (ngr2==0.or.nvxcgrho/=3) then
     msg='The values of nvxcgrho or ngr2 are not compatible with GGA!'
     MSG_BUG(msg)
   end if
 end if

!Call LibXC routines
 if (libxc_functionals_isgga()) then
   if (order**2<=1) then
     call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho,exc,vxcrho,&
&               grho2=grho2,vxcgr=vxcgrho)
   else
     call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho,exc,vxcrho,&
&               grho2=grho2,vxcgr=vxcgrho,dvxc=dvxc)
   end if
 else
   if (order**2<=1) then
     call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho,exc,vxcrho)
   else if (order**2<=4) then
     call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho,exc,vxcrho,&
&                                  dvxc=dvxc)
   else
     call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho,exc,vxcrho,&
&                                  dvxc=dvxc,d2vxc=d2vxc)
   end if
 end if

end subroutine pawxc_drivexc_libxc
!!***
#endif

end subroutine pawxc_drivexc_wrapper
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxc_rotate_mag
!! NAME
!! pawxc_rotate_mag
!!
!! FUNCTION
!!  Project (rotate) a non-collinear density (stored as density+magn.)
!!   on a magnetization and give a collinear density (stored as [up,dn] or [up+dn,up]).
!!
!! INPUTS
!!  rho_in(vectsize,4)=input non-collinear density and magnetization
!!  mag(vectsize,3)=magnetization used for projection
!!  vectsize=size of vector fields
!!  [rho_out_format]= 1=rho_out is stored as [up,dn]
!!                    2=rho_out is stored as [up+dn,up]
!!                    Default=1
!!
!! OUTPUT
!!  rho_out(vectsize,2)=output (projected, collinear) density
!!  [mag_norm_out(vectsize)]= --optional-- norm of mag(:) at each point of the grid
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

 subroutine pawxc_rotate_mag(rho_in,rho_out,mag,vectsize,mag_norm_out,rho_out_format)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: vectsize
 integer,intent(in),optional :: rho_out_format
!arrays
 real(dp),intent(in) :: rho_in(vectsize,4),mag(vectsize,3)
 real(dp),intent(out) :: rho_out(vectsize,2)
 real(dp),intent(out),optional :: mag_norm_out(vectsize)

!Local variables-------------------------------
!scalars
#if ! defined HAVE_LIBPAW_ABINIT
 integer :: ipt
 real(dp),parameter :: m_norm_min=tol8
 real(dp) :: m_norm,rhoin_dot_mag,rho_up
#endif
!arrays

! *************************************************************************

!One could add here a section for other codes (i.e. BigDFT, ...)
#if defined HAVE_LIBPAW_ABINIT
 if (present(rho_out_format).and.present(mag_norm_out)) then
   call rotate_mag(rho_in,rho_out,mag,vectsize,1, &
&          rho_out_format=rho_out_format,mag_norm_out=mag_norm_out)
 else if (present(rho_out_format).and..not.present(mag_norm_out)) then
   call rotate_mag(rho_in,rho_out,mag,vectsize,1,rho_out_format=rho_out_format)
 else if (.not.present(rho_out_format).and.present(mag_norm_out)) then
   call rotate_mag(rho_in,rho_out,mag,vectsize,1,mag_norm_out=mag_norm_out)
 else
   call rotate_mag(rho_in,rho_out,mag,vectsize,1)
 end if
#else
 do ipt=1,vectsize
   m_norm=sqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
   rhoin_dot_mag=rho_in(ipt,2)*mag(ipt,1)+rho_in(ipt,3)*mag(ipt,2) &
&               +rho_in(ipt,4)*mag(ipt,3)
   if(m_norm>m_norm_min)then
     rho_out(ipt,1)=half*(rho_in(ipt,1)+rhoin_dot_mag/m_norm)
     rho_out(ipt,2)=half*(rho_in(ipt,1)-rhoin_dot_mag/m_norm)
   else
     rho_out(ipt,1)=half*rho_in(ipt,1)
     rho_out(ipt,2)=half*rho_in(ipt,1)
   end if
   if (present(mag_norm_out).and.m_norm> m_norm_min) mag_norm_out(ipt)=m_norm
   if (present(mag_norm_out).and.m_norm<=m_norm_min) mag_norm_out(ipt)=zero
 end do
 if (present(rho_out_format)) then
   if (rho_out_format==2) then
     do ipt=1,vectsize
       rho_up=rho_out(ipt,1)
       rho_out(ipt,1)=rho_up+rho_out(ipt,2)
       rho_out(ipt,2)=rho_up
     end do
   end if
 end if
#endif

end subroutine pawxc_rotate_mag
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxc_rotate_back_mag
!! NAME
!! pawxc_rotate_back_mag
!!
!! FUNCTION
!!  Rotate back a collinear XC potential (stored as up+dn) with respect to
!!   a magnetization and give a non-collinear XC potential
!!   (stored as up_up, dn_dn, Re{up_dn}, Im{up_dn}).
!!
!! INPUTS
!!  vxc_in(vectsize,2)=input collinear XC potential
!!  mag(vectsize,3)=magnetization used for projection
!!  vectsize=size of vector fields
!!
!! OUTPUT
!!  vxc_out(vectsize,4)=output non-collinear XC potential
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

 subroutine pawxc_rotate_back_mag(vxc_in,vxc_out,mag,vectsize)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: vectsize
!arrays
 real(dp),intent(in) :: vxc_in(vectsize,2),mag(vectsize,3)
 real(dp),intent(out) :: vxc_out(vectsize,4)

!Local variables-------------------------------
!scalars
#if ! defined HAVE_LIBPAW_ABINIT
 integer :: ipt
 real(dp),parameter :: m_norm_min=tol8
 real(dp) :: dvdn,dvdz,m_norm
#endif
!arrays

! *************************************************************************

!One could add here a section for other codes (i.e. BigDFT, ...)
#if defined HAVE_LIBPAW_ABINIT
 call rotate_back_mag(vxc_in,vxc_out,mag,vectsize)
#else
 do ipt=1,vectsize
   m_norm=sqrt(mag(ipt,1)**2+mag(ipt,2)**2+mag(ipt,3)**2)
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
#endif

end subroutine pawxc_rotate_back_mag
!!***

!----------------------------------------------------------------------

!!****f* m_pawxc/pawxc_rotate_back_mag_dfpt
!! NAME
!! pawxc_rotate_back_mag_dfpt
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
!!
!! OUTPUT
!!  vxc1_out(vectsize,4)=output 1st-order non-collinear XC potential
!!
!! PARENTS
!!      m_pawxc
!!
!! CHILDREN
!!      rotate_back_mag_dfpt
!!
!! SOURCE

 subroutine pawxc_rotate_back_mag_dfpt(vxc1_in,vxc1_out,vxc,kxc,rho1,mag,vectsize)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: vectsize
!arrays
 real(dp),intent(in) :: kxc(:,:),mag(vectsize,3),rho1(vectsize,4)
 real(dp),intent(in) :: vxc(vectsize,4),vxc1_in(vectsize,2)
 real(dp),intent(out) :: vxc1_out(vectsize,4)

!Local variables-------------------------------
!scalars
#if ! defined HAVE_LIBPAW_ABINIT
 character(len=100) :: msg
#endif
!arrays

! *************************************************************************

!One could add here a section for other codes (i.e. BigDFT, ...)
#if defined HAVE_LIBPAW_ABINIT
 call rotate_back_mag_dfpt(1,vxc1_in,vxc1_out,vxc,kxc,rho1,mag,vectsize,1)
#else
 msg='[LIBPAW] Non-collinear DFPT not available (only in ABINIT)!'
 MSG_ERROR(msg)
#endif

end subroutine pawxc_rotate_back_mag_dfpt
!!***

!----------------------------------------------------------------------

end module m_pawxc
!!***
