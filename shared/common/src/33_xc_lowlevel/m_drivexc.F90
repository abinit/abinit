!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_drivexc
!! NAME
!!  m_drivexc
!!
!! FUNCTION
!! Driver of XC functionals. Optionally, deliver the XC kernel, or even the derivative
!! of the XC kernel (the third derivative of the XC energy)
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2019 ABINIT group (MT, MJV, CE, TD, XG)
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

module m_drivexc

 use defs_basis
 use m_abicore
 use m_errors
 use libxc_functionals
 use m_numeric_tools,    only : invcb
 use m_xciit,            only : xciit
 use m_xcpbe,            only : xcpbe
 use m_xchcth,           only : xchcth
 use m_xclda,  only : xcpzca, xcspol, xctetr, xcwign, xchelu, xcxalp, xclb

 implicit none

 private
!!***

 public :: drivexc         ! Driver of XC functionals. Optionally, deliver the XC kernel, or even the derivative
 public :: echo_xc_name    ! Write to log and output the xc functional which will be used for this dataset
 public :: check_kxc       ! Given a XC functional (defined by ixc), check if Kxc (dVxc/drho) is avalaible.
 public :: size_dvxc       ! Give the size of the array dvxc(npts,ndvxc) and the second dimension of the d2vxc(npts,nd2vxc)
 public :: xcmult          ! (GGA) Multiply the different gradient of spin-density by the derivative of the XC functional
                           ! with respect to the norm of the gradient, then divide it by the norm of the gradient
 public :: mkdenpos        ! Make a ground-state density positive everywhere.
!!***

contains
!!***

!!****f* m_drivexc/echo_xc_name
!! NAME
!! echo_xc_name
!!
!! FUNCTION
!!  Write to log and output the xc functional which will be used for this dataset
!!
!! INPUTS
!!  ixc = internal code for xc functional
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine echo_xc_name (ixc)

!Arguments -------------------------------
 integer, intent(in) :: ixc

!Local variables -------------------------
 integer :: l_citation
 character(len=500) :: message, citation

! *********************************************************************

 message =''
 citation =''

!normal case (not libxc)
 if (ixc >= 0) then

   select case (ixc)
   case (0)
     message = 'No xc applied (usually for testing) - ixc=0'
     citation = ''
!      LDA,LSD
   case (1)
     message = 'LDA: new Teter (4/93) with spin-polarized option - ixc=1'
     citation = 'S. Goedecker, M. Teter, J. Huetter, PRB 54, 1703 (1996)' ! [[cite:Goedecker1996]]
   case (2)
     message = 'LDA: Perdew-Zunger-Ceperley-Alder - ixc=2'
     citation = 'J.P.Perdew and A.Zunger, PRB 23, 5048 (1981) ' ! [[cite:Perdew1981]]
   case (3)
     message = 'LDA: old Teter (4/91) fit to Ceperley-Alder data - ixc=3'
     citation = ''
   case (4)
     message = 'LDA: Wigner - ixc=4'
     citation = 'E.P.Wigner, Trans. Faraday Soc. 34, 67 (1938)' ! [[cite:Wigner1938]]
   case (5)
     message = 'LDA: Hedin-Lundqvist - ixc=5'
     citation = 'L.Hedin and B.I.Lundqvist, J. Phys. C4, 2064 (1971)' ! [[cite:Hedin1971]]
   case (6)
     message = 'LDA: "X-alpha" xc - ixc=6'
     citation = 'Slater J. C., Phys. Rev. 81, 385 (1951)' ! [[cite:Slater1951]]
   case (7)
     message = 'LDA: Perdew-Wang 92 LSD fit to Ceperley-Alder data - ixc=7'
     citation = 'J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)' ! [[cite:Perdew1992a]]
   case (8)
     message = 'LDA: Perdew-Wang 92 LSD , exchange-only - ixc=8'
     citation = 'J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)' ! [[cite:Perdew1992a]]
   case (9)
     message = 'LDA: Perdew-Wang 92 Ex+Ec_RPA  energy - ixc=9'
     citation = 'J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)' ! [[cite:Perdew1992]]
   case (10)
     message = 'LDA: RPA LSD energy (only the energy !!) - ixc=10'
     citation = ''
!      GGA
   case (11)
     message = 'GGA: Perdew-Burke-Ernzerhof functional - ixc=11'
     citation = 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' ! [[cite:Perdew1996]]
   case (12)
     message = 'GGA: x-only Perdew-Burke-Ernzerhof functional - ixc=12'
     citation = 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)' ! [[cite:Perdew1996]]
   case (13)
     message = 'GGA: LDA (ixc==7) energy, and the xc _potential_ is given by van Leeuwen-Baerends GGA - ixc=13'
     citation = 'R. van Leeuwen and E. J. Baerends PRA 49, 2421 (1994)' ! [[cite:VanLeeuwen1994]]
   case (14)
     message = 'GGA: revPBE functional - ixc=14'
     citation = 'Zhang and Yang, PRL 80, 890 (1998)' ! [[cite:Zhang1998]]
   case (15)
     message = 'GGA: RPBE functional - ixc=15'
     citation = 'Hammer, L. B. Hansen, and J. K. Norskov, PRB 59, 7413 (1999)' ! [[cite:Hammer1999]]
   case (16)
     message = 'GGA: HCTH93 functional - ixc=16'
     citation = 'F.A. Hamprecht, A.J. Cohen, D.J. Tozer, N.C. Handy, JCP 109, 6264 (1998)' ! [[cite:Hamprecht1998]]
   case (17)
     message = 'GGA: HCTH120 functional - ixc=17'
     citation = 'A.D. Boese, N.L. Doltsinis, N.C. Handy, and M. Sprik, JCP 112, 1670 (2000)' ! [[cite:Boese2000]]
   case (23)
     message = 'GGA: Wu Cohen functional - ixc=23'
     citation = 'Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)' ! [[cite:Wu2006]]
   case (24)
     message = 'GGA: C09x exchange functional - ixc=24'
     citation = 'Valentino R. Cooper, PRB 81, 161104(R) (2010)' ! [[cite:Cooper2010]]
   case (26)
     message = 'GGA: HCTH147 functional - ixc=26'
     citation = 'A.D. Boese, N.L. Doltsinis, N.C. Handy, and M. Sprik, JCP 112, 1670 (2000)' ! [[cite:Boese2000]]
   case (27)
     message = 'GGA: HCTH407 functional - ixc=27'
     citation = 'A.D. Boese, and N.C. Handy, JCP 114, 5497 (2001)' ! [[cite:Boese2001]]
!      Fermi-Amaldi
   case (20)
     message = 'Fermi-Amaldi correction - ixc=20'
     citation = ''
   case (21)
     message = 'Fermi-Amaldi correction with LDA(ixc=1) kernel - ixc=21'
     citation = ''
   case (22)
     message = 'Fermi-Amaldi correction with hybrid BPG kernel - ixc=22'
     citation = ''
   case (31)
     message = 'Meta-GGA fake1 - ixc=31'
     citation = ''
   case (32)
     message = 'Meta-GGA fake2 - ixc=32'
     citation = ''
   case (33)
     message = 'Meta-GGA fake3 - ixc=33'
     citation = ''
   case (34)
     message = 'Meta-GGA fake4 - ixc=34'
     citation = ''
   case (40)
     message = 'Hartree-Fock with mixing coefficient alpha=1'
     citation = ''
   case (41)
     message = 'PBE0 with alpha=0.25'
     citation = ''
   case (42)
     message = 'modified PBE0 with alpha=0.33'
     citation = ''
   case (50)
     message = 'LDA at finite T Ichimaru-Iyetomy-Tanaka - ixc=50'
     citation = 'Ichimaru S., Iyetomi H., Tanaka S., Phys. Rep. 149, 91-205 (1987) ' ! [[cite:Ichimaru1987]]
   case default
     write(message,'(a,i0)')" echo_xc_name does not know how to handle ixc = ",ixc
     MSG_WARNING(message)
   end select

   message = " Exchange-correlation functional for the present dataset will be:" // ch10 &
&   // "  " // trim(message)

   l_citation=len_trim(citation)
   citation = " Citation for XC functional:" // ch10 // "  " // trim(citation)

   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   if(l_citation/=0)then
     call wrtout(ab_out,citation,'COLL')
     call wrtout(std_out,citation,'COLL')
   end if

   message =' '
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

 end if ! end libxc if

end subroutine echo_xc_name
!!***

!!****f* m_drivexc/check_kxc
!! NAME
!! check_kxc
!!
!! FUNCTION
!!  Given a XC functional (defined by ixc), check if Kxc (dVxc/drho) is avalaible.
!!
!! INPUTS
!!  ixc = internal code for xc functional
!!  optdriver=type of calculation (ground-state, response function, GW, ...)
!!
!! OUTPUT
!!
!! PARENTS
!!      respfn,scfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine check_kxc(ixc,optdriver)

!Arguments -------------------------------
 integer, intent(in) :: ixc,optdriver

!Local variables -------------------------
 logical :: kxc_available
 character(len=500) :: msg

! *********************************************************************

 kxc_available=.false.

 if (ixc>=0) then
   kxc_available=(ixc/=16.and.ixc/=17.and.ixc/=26.and.ixc/=27)
   if (.not.kxc_available) then
     write(msg,'(a,i0,3a)') &
&     'The selected XC functional (ixc=',ixc,')',ch10,&
&     'does not provide Kxc (dVxc/drho) !'
   end if
 else if (ixc==-406.or.ixc==-427.or.ixc==-428.or.ixc==-456)then
   kxc_available=.true.
 else ! ixc<0 and not one of the allowed hybrids
   kxc_available=libxc_functionals_has_kxc()
   if (.not.kxc_available) then
     write(msg,'(a,i0,7a)') &
&     'The selected XC functional (ixc=',ixc,'):',ch10,&
&     '   <<',trim(libxc_functionals_fullname()),'>>',ch10,&
&     'does not provide Kxc (dVxc/drho) !'
   end if
 end if

 if (.not.kxc_available) then
   write(msg,'(7a)') trim(msg),ch10,&
&   'However, with the current input options, ABINIT needs Kxc.',ch10,&
&   '>Possible action:',ch10,&
&   'Change the XC functional in psp file or input file.'
   if (optdriver==0) then
     write(msg,'(13a)') trim(msg),ch10,&
&     '>Possible action (2):',ch10,&
&     'If you are using density mixing for the SCF cycle',ch10,&
&     '(iscf>=10, which is the default for PAW),',ch10,&
&     'change to potential mixing (iscf=7, for instance).',ch10,&
&     '>Possible action (3):',ch10,&
&     'Switch to another value of densfor_pred (=5, for instance).'
   end if
   MSG_ERROR(msg)
 end if

end subroutine check_kxc
!!***

!!****f* m_drivexc/size_dvxc
!! NAME
!! size_dvxc
!!
!! FUNCTION
!! Give the sizes of the several arrays involved in exchange-correlation calculation
!! needed to allocated them for the drivexc routine
!!
!! INPUTS
!!  ixc= choice of exchange-correlation scheme
!!  order= gives the maximal derivative of Exc computed.
!!    1=usual value (return exc and vxc)
!!    2=also computes the kernel (return exc,vxc,kxc)
!!   -2=like 2, except (to be described)
!!    3=also computes the derivative of the kernel (return exc,vxc,kxc,k3xc)
!!  nspden= number of spin components
!!  [xc_funcs(2)]= <type(libxc_functional_type)>
!!  [add_tfw]= optional flag controling the addition of Weiszacker gradient correction to Thomas-Fermi XC energy
!!
!! OUTPUT
!!  --- All optionals
!!  [usegradient]= [flag] 1 if the XC functional needs the gradient of the density (grho2_updn)
!!  [uselaplacian]= [flag] 1 if the XC functional needs the laplacian of the density (lrho_updn)
!!  [usekden]= [flag] 1 if the XC functional needs the kinetic energy density (lrho_updn)
!!  [nvxcgrho]= size of the array dvxcdgr(npts,nvxcgrho) (derivative of Exc wrt to gradient)
!!  [nvxclrho]= size of the array dvxclpl(npts,nvxclrho) (derivative of Exc wrt to laplacian)
!!  [nvxctau]= size of the array dvxctau(npts,nvxctau) (derivative of Exc wrt to kin. ener. density)
!!  [ndvxc]= size of the array dvxc(npts,ndvxc) (second derivatives of Exc wrt to density and gradient)
!!  [nd2vxc]= size of the array d2vxc(npts,nd2vxc) (third derivatives of Exc wrt density)
!!
!! PARENTS
!!      m_pawxc,rhotoxc
!!
!! CHILDREN
!!
!! SOURCE

subroutine size_dvxc(ixc,order,nspden,&
&          usegradient,uselaplacian,usekden,&
&          nvxcgrho,nvxclrho,nvxctau,ndvxc,nd2vxc,&
&          add_tfw,xc_funcs) ! Optional

!Arguments----------------------
 integer,intent(in) :: ixc,nspden,order
 integer,intent(out),optional :: nvxcgrho,nvxclrho,nvxctau,ndvxc,nd2vxc
 integer,intent(out),optional :: usegradient,uselaplacian,usekden
 logical, intent(in),optional :: add_tfw
 type(libxc_functional_type),intent(in),optional :: xc_funcs(2)

!Local variables----------------
 logical :: libxc_isgga,libxc_ismgga,libxc_ishybrid,my_add_tfw
 logical :: need_gradient,need_laplacian,need_kden

! *************************************************************************

!Several flags
 my_add_tfw=.false.;if (present(add_tfw)) my_add_tfw=add_tfw
 libxc_isgga=.false. ; libxc_ismgga=.false. ; libxc_ishybrid=.false.
 if(ixc<0)then
   if(present(xc_funcs))then
     libxc_isgga=libxc_functionals_isgga(xc_functionals=xc_funcs)
     libxc_ismgga=libxc_functionals_ismgga(xc_functionals=xc_funcs)
     libxc_ishybrid=libxc_functionals_is_hybrid(xc_functionals=xc_funcs)
   else
     libxc_isgga=libxc_functionals_isgga()
     libxc_ismgga=libxc_functionals_ismgga()
     libxc_ishybrid=libxc_functionals_is_hybrid()
   end if
 end if

!Do we use the gradient?
 need_gradient=((ixc>=11.and.ixc<=17).or.(ixc==23.or.ixc==24).or. &
&               (ixc==26.or.ixc==27).or.(ixc>=31.and.ixc<=34).or. &
&               (ixc==41.or.ixc==42).or.ixc==1402000)
 if (ixc<0.and.(libxc_isgga.or.libxc_ismgga.or.libxc_ishybrid)) need_gradient=.true.
 if (my_add_tfw) need_gradient=.true.
 if (present(usegradient)) usegradient=merge(1,0,need_gradient)

!Do we use the laplacian?
 need_laplacian=(ixc==32)
 if (ixc<0) then
   if(present(xc_funcs)) need_laplacian=libxc_functionals_needs_laplacian(xc_functionals=xc_funcs)
   if(.not.present(xc_funcs)) need_laplacian=libxc_functionals_needs_laplacian()
 end if
 if (present(uselaplacian)) uselaplacian=merge(1,0,need_laplacian)

!Do we use the kinetic energy density?
 need_kden=(ixc==31.or.ixc==34)
 if (ixc<0) need_kden=libxc_ismgga
 if (present(usekden)) usekden=merge(1,0,need_kden)

!First derivative(s) of XC functional wrt gradient of density
 if (present(nvxcgrho)) then
   nvxcgrho=0
   if (abs(order)>=1) then
     if (need_gradient) nvxcgrho=3
     if (ixc==13) nvxcgrho=0
     if (ixc==16.or.ixc==17.or.ixc==26.or.ixc==27) nvxcgrho=2
   end if
 end if

!First derivative(s) of XC functional wrt laplacian of density
 if (present(nvxclrho)) then
   nvxclrho=0
   if (abs(order)>=1) then
     if (need_laplacian) nvxclrho=min(nspden,2)
   end if
 end if

!First derivative(s) of XC functional wrt kinetic energy density
 if (present(nvxctau)) then
   nvxctau=0
   if (abs(order)>=1) then
     if (need_kden) nvxctau=min(nspden,2)
   end if
 end if

!Second derivative(s) of XC functional wrt density
 if (present(ndvxc)) then
   ndvxc=0
   if (abs(order)>=2) then
     if (ixc==1.or.ixc==13.or.ixc==21.or.ixc==22.or.(ixc>=7.and.ixc<=10)) then
       ndvxc=min(nspden,2)+1
     else if ((ixc>=2.and.ixc<=6).or.ixc==50) then
       ndvxc=1
     else if (ixc==12.or.ixc==24) then
       ndvxc=8
     else if (ixc==11.or.ixc==12.or.ixc==14.or.ixc==15.or. &
&             ixc==23.or.ixc==41.or.ixc==42.or.ixc==1402000) then
       ndvxc=15
     else if (ixc<0) then
       ndvxc=3 ; if (need_gradient) ndvxc=15
     end if
   end if
 end if

!Third derivative(s) of XC functional wrt density
 if (present(nd2vxc)) then
   nd2vxc=0
   if (abs(order)>=3) then
     if (ixc==3.or.(ixc>=11.and.ixc<=15.and.ixc/=13).or. &
&        ixc==23.or.ixc==24.or.ixc==41.or.ixc==42) then
       nd2vxc=1
     else if ((ixc>=7.and.ixc<=10).or.ixc==13.or.ixc==1402000) then
       nd2vxc=3*min(nspden,2)-2
     else if (ixc<0) then
       if (.not.need_gradient) nd2vxc=3*min(nspden,2)-2
     end if
   end if
 end if

end subroutine size_dvxc
!!***

!!****f* m_drivexc/xcmult
!! NAME
!! xcmult
!!
!! FUNCTION
!! In the case of GGA, multiply the different gradient of spin-density
!! by the derivative of the XC functional with respect
!! to the norm of the gradient, then divide it by the norm of the gradient
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
!!      m_pawxc,rhotoxc
!!
!! CHILDREN
!!
!! SOURCE

subroutine xcmult (depsxc,nfft,ngrad,nspden,nspgrad,rhonow)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ngrad,nspden,nspgrad
!arrays
 real(dp),intent(in) :: depsxc(nfft,nspgrad)
 real(dp),intent(inout) :: rhonow(nfft,nspden,ngrad*ngrad)

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

end subroutine xcmult
!!***

!!****f* m_drivexc/mkdenpos
!! NAME
!! mkdenpos
!!
!! FUNCTION
!! Make a ground-state density positive everywhere:
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
!!  (see side effects)
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
!! NOTES
!!  At this stage, rhonow(:,1:nspden) contains the density in real space,
!!  on the unshifted or shifted grid. Now test for negative densities
!!  Note that, ignoring model core charge, as long as boxcut>=2
!!  the shifted density is derivable from the square of a Fourier
!!  interpolated charge density => CANNOT go < 0.
!!  However, actually can go < 0 to within machine precision;
!!  do not print useless warnings in this case, just fix it.
!!  Fourier interpolated core charge can go < 0 due to Gibbs
!!  oscillations; could avoid this by recomputing the model core
!!  charge at the new real space grid points (future work).
!!
!! PARENTS
!!      bethe_salpeter,m_pawxc,mkcore_wvl,posdoppler,poslifetime,posratecore
!!      psolver_rhohxc,rhohxcpositron,rhotoxc,wvl_initro
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkdenpos(iwarn,nfft,nspden,option,rhonow,xc_denpos)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nspden,option
 integer,intent(inout) :: iwarn
 real(dp),intent(in) :: xc_denpos
!arrays
 real(dp),intent(inout) :: rhonow(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ifft,ispden,numneg
 real(dp) :: rhotmp,worst
 character(len=500) :: message
!arrays
 real(dp) :: rho(2)

! *************************************************************************

 numneg=0
 worst=zero

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

!    rhonow is stored as (up+dn,up)
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
   MSG_BUG('nspden>2 not allowed !')
 end if ! End choice between non-spin polarized and spin-polarized.

 if (numneg>0) then
   if (iwarn==0) then
     write(message,'(a,i0,a,a,a,es10.2,a,e10.2,a,a,a,a)')&
&     'Density went too small (lower than xc_denpos) at ',numneg,' points',ch10,&
&     'and was set to xc_denpos = ',xc_denpos,'. Lowest was ',worst,'.',ch10,&
&     'Likely due to too low boxcut or too low ecut for',' pseudopotential core charge.'
     MSG_WARNING(message)
   end if
   iwarn=iwarn+1
 end if

end subroutine mkdenpos
!!***

!!****f* m_drivexc/drivexc
!! NAME
!! drivexc
!!
!! FUNCTION
!! Driver of XC functionals. Treat spin-polarized as well as non-spin-polarized.
!! Treat local approximations, GGAs, meta-GGAs or hybrid functionals.
!! Optionally, deliver the XC kernel, or even the derivative
!! of the XC kernel (the third derivative of the XC energy)
!!
!! INPUTS
!!  ixc=index of the XC functional
!!  xclevel=XC functional level (lda, gga, etc...)
!!  usegradient=[flag] 1 if the XC functional depends on density gradient (grho2_updn)
!!  uselaplacian=[flag] 1 if the XC functional depends on density laplacian (lrho_updn)
!!  usekden=[flag] 1 if the XC functional depends on kinetic energy density (tau_updn)
!!  order=gives the maximal derivative of Exc computed.
!!    1=usual value (return exc and vxc)
!!    2=also computes the kernel (return exc,vxc,kxc)
!!   -2=like 2, except (to be described)
!!    3=also computes the derivative of the kernel (return exc,vxc,kxc,k3xc)
!!  npts=number of real space points on which the density is provided
!!  nspden=number of spin-density components (1 or 2)
!!  nvxcgrho=number of components of 1st-derivative of Exc wrt density gradient (nvxcgrho)
!!  nvxclrho=number of components of 1st-derivative of Exc wrt density laplacian (nvxclrho)
!!  nvxctau=number of components of 1st-derivative of Exc wrt kinetic energy density (nvxctau)
!!  ndvxc=number of components of  1st-derivative of Vxc (dvxc)
!!  nd2vxc=number of components of  2nd-derivative of Vxc (d2vxc)
!!  rho_updn(npts,nspden)=spin-up and spin-down densities
!!    In the calling routine, spin-down density must be equal to spin-up density.
!!    If nspden=1, only spin-up density must be given (half the total density).
!!    If nspden=2, spin-up and spin-down densities must be given.
!!  === Optional input arguments ===
!!  [grho2_updn(npts,(2*nspden-1)*usegradient)]=the square of the gradients
!!    of spin-up, spin-down, and total density.
!!    If nspden=1, only the square of the gradient of the spin-up density must be given.
!!     In the calling routine, the square of the gradient of the spin-down density must be equal
!!     to the square of the gradient of the spin-up density, and both must be equal to
!!     one-quarter of the square of the gradient of the total density.
!!    If nspden=2, the square of the gradients of spin-up, spin-down, and total density must be given.
!!     Note that the square of the gradient of the total density is usually NOT related to
!!     the square of the gradient of the spin-up and spin-down densities, because the gradients
!!     are not usually aligned. This is not the case when nspden=1.
!!  [lrho_updn(npts,nspden*uselaplacian)]=the Laplacian of spin-up and spin-down densities.
!!    If nspden=1, only the spin-up Laplacian density must be given and must
!!     be equal to the spin-up Laplacian density.
!!    If nspden=2, the Laplacian of spin-up and spin-down densities must be given.
!!  [tau_updn(npts,nspden*usekden)]=the spin-up and spin-down kinetic energy densities.
!!    If nspden=1, only the spin-up kinetic energy density must be given and must
!!     be equal to the half the total kinetic energy density.
!!    If nspden=2, the spin-up and spin-down kinetic energy densities must be given.
!!  [exexch]=choice of <<<local>>> exact exchange. Active if exexch=3 (only for GGA, and NOT for libxc)
!!  [el_temp]= electronic temperature (to be used for finite temperature XC functionals)
!!  [hyb_mixing]= mixing parameter for the native PBEx functionals (ixc=41 and 42)
!!  [xc_funcs(2)]= <type(libxc_functional_type)>: libxc XC functionals.
!!
!! OUTPUT
!!  exc(npts)=exchange-correlation energy density (hartree)
!!  vxcrho(npts,nspden)= (d($\rho$*exc)/d($\rho_up$)) (hartree)
!!                  and  (d($\rho$*exc)/d($\rho_down$)) (hartree)
!!  === Optional output arguments ===
!!  [vxcgrho(npts,nvxcgrho)]=1st-derivative of the xc energy wrt density gradient>
!!                          = 1/$|grad \rho_up|$ (d($\rho$*exc)/d($|grad \rho_up|$))
!!                            1/$|grad \rho_dn|$ (d($\rho$*exc)/d($|grad \rho_dn|$))
!!                            1/$|grad \rho|$ (d($\rho$*exc)/d($|grad \rho|$))
!!  [vxclrho(npts,nvxclrho)]=1st-derivative of the xc energy wrt density laplacian.
!!                          = d($\rho$*exc)/d($\lrho_up$)
!!                            d($\rho$*exc)/d($\lrho_down$)
!!  [vxctau(npts,nvxctau)]=1st-derivative of the xc energy wrt kinetic energy density.
!!                        = d($\rho$*exc)/d($\tau_up$)
!!                          d($\rho$*exc)/d($\tau_down$)
!!  [dvxc(npts,ndvxc)]=partial second derivatives of the XC energy
!!   === Only if abs(order)>1 ===
!!   In case of local energy functional (option=1,-1 or 3):
!!    dvxc(npts,1+nspden)=
!!     if(nspden=1 .and. order==2): dvxci(:,1)=dvxc/d$\rho$ , dvxc(:,2) empty
!!     if(nspden=1 .and. order==-2): also compute dvxci(:,2)=dvxc($\uparrow$)/d$\rho(\downarrow)$
!!     if(nspden=2): dvxc(:,1)=dvxc($\uparrow$)/d$\rho(\downarrow)$,
!!                   dvxc(:,2)=dvxc($\uparrow$)/d$\rho(\downarrow)$,
!!                   dvxc(:,3)=dvxc($\downarrow$)/d$\rho(\downarrow)$
!!   In case of gradient corrected functional (option=2,-2, 4, -4, 5, 6, 7):
!!    dvxc(npts,15)=
!!     dvxc(:,1)= d2Ex/drho_up drho_up
!!     dvxc(:,2)= d2Ex/drho_dn drho_dn
!!     dvxc(:,3)= dEx/d(abs(grad(rho_up))) / abs(grad(rho_up))
!!     dvxc(:,4)= dEx/d(abs(grad(rho_dn))) / abs(grad(rho_dn))
!!     dvxc(:,5)= d2Ex/d(abs(grad(rho_up))) drho_up / abs(grad(rho_up))
!!     dvxc(:,6)= d2Ex/d(abs(grad(rho_dn))) drho_dn / abs(grad(rho_dn))
!!     dvxc(:,7)= 1/abs(grad(rho_up)) * d/d(abs(grad(rho_up)) (dEx/d(abs(grad(rho_up))) /abs(grad(rho_up)))
!!     dvxc(:,8)= 1/abs(grad(rho_dn)) * d/d(abs(grad(rho_dn)) (dEx/d(abs(grad(rho_dn))) /abs(grad(rho_dn)))
!!     dvxc(:,9)= d2Ec/drho_up drho_up
!!     dvxc(:,10)=d2Ec/drho_up drho_dn
!!     dvxc(:,11)=d2Ec/drho_dn drho_dn
!!     dvxc(:,12)=dEc/d(abs(grad(rho))) / abs(grad(rho))
!!     dvxc(:,13)=d2Ec/d(abs(grad(rho))) drho_up / abs(grad(rho))
!!     dvxc(:,14)=d2Ec/d(abs(grad(rho))) drho_dn / abs(grad(rho))
!!     dvxc(:,15)=1/abs(grad(rho)) * d/d(abs(grad(rho)) (dEc/d(abs(grad(rho))) /abs(grad(rho)))
!!  [d2vxc(npts,nd2vxc)]=second derivative of the XC potential=3rd order derivative of XC energy
!!   === Only if abs(order)>1 ===
!!   === At present only available for LDA ===
!!    if nspden=1 d2vxc(npts,1)=second derivative of the XC potential=3rd order derivative of energy
!!    if nspden=2 d2vxc(npts,1), d2vxc(npts,2), d2vxc(npts,3), d2vxc(npts,4) (3rd derivative of energy)
!!  [fxcT(npts)]=XC free energy of the electron gaz at finite temperature (to be used for plasma systems)
!!
!! PARENTS
!!      drivexc_main
!!
!! CHILDREN
!!      invcb,libxc_functionals_end,libxc_functionals_getvxc
!!      libxc_functionals_init,xchcth,xchelu,xciit,xclb,xcpbe,xcpzca,xcspol
!!      xctetr,xcwign,xcxalp
!!
!! SOURCE

subroutine drivexc(ixc,order,npts,nspden,usegradient,uselaplacian,usekden,&
&          rho_updn,exc,vxcrho,nvxcgrho,nvxclrho,nvxctau,ndvxc,nd2vxc, &            ! mandatory arguments
&          grho2_updn,vxcgrho,lrho_updn,vxclrho,tau_updn,vxctau,dvxc,d2vxc, &  ! optional arguments
&          exexch,el_temp,fxcT,hyb_mixing,xc_funcs)                            ! optional parameters

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,npts,nspden
 integer,intent(in) :: ndvxc,nd2vxc,nvxcgrho,nvxclrho,nvxctau,order
 integer,intent(in) :: usegradient,uselaplacian,usekden
 integer,intent(in),optional :: exexch
 real(dp),intent(in),optional :: el_temp,hyb_mixing
!arrays
 real(dp),intent(in) :: rho_updn(npts,nspden)
 real(dp),intent(in),optional :: grho2_updn(npts,(2*nspden-1)*usegradient)
 real(dp),intent(in),optional :: lrho_updn(npts,nspden*uselaplacian),tau_updn(npts,nspden*usekden)
 real(dp),intent(out) :: exc(npts),vxcrho(npts,nspden)
 real(dp),intent(out),optional :: dvxc(npts,ndvxc),d2vxc(npts,nd2vxc),fxcT(:)
 real(dp),intent(out),optional :: vxcgrho(npts,nvxcgrho),vxclrho(npts,nvxclrho),vxctau(npts,nvxctau)
 type(libxc_functional_type),intent(inout),optional :: xc_funcs(2)

!Local variables-------------------------------
!scalars
 integer :: ispden,ixc_from_lib,ixc1,ixc2,ndvxc_x
 integer :: my_exexch,need_ndvxc,need_nd2vxc,need_nvxcgrho,need_nvxclrho,need_nvxctau
 integer :: need_gradient,need_laplacian,need_kden,optpbe
 logical :: has_gradient,has_laplacian,has_kden,libxc_test
 real(dp) :: alpha,my_hyb_mixing
 real(dp),parameter :: rsfac=0.6203504908994000e0_dp
 character(len=500) :: message
!arrays
 real(dp),allocatable :: exci_rpa(:),rhotot(:),rspts(:),vxci_rpa(:,:),zeta(:)
 real(dp),allocatable :: exc_c(:),exc_x(:),vxcrho_c(:,:),vxcrho_x(:,:)
 real(dp),allocatable :: d2vxc_c(:,:),d2vxc_x(:,:),dvxc_c(:,:),dvxc_x(:,:)
 real(dp),allocatable :: vxcgrho_x(:,:)
 type(libxc_functional_type) :: xc_funcs_vwn3(2),xc_funcs_lyp(2)

! *************************************************************************

!optional arguments
 my_exexch=0;if(present(exexch)) my_exexch=exexch
 my_hyb_mixing=0
 if (ixc==41) my_hyb_mixing=quarter
 if (ixc==42) my_hyb_mixing=third
 if (present(hyb_mixing)) my_hyb_mixing=hyb_mixing

! =================================================
! ==         Compatibility tests                 ==
! =================================================

!Check libXC initialization
 if (ixc<0 .or. ixc==1402) then
   libxc_test=libxc_functionals_check(stop_if_error=.true.)
 end if

! Check libXC consistency between ixc passed in input
!  and the one used to initialize the libXC library
 if (ixc<0) then
   ixc_from_lib=libxc_functionals_ixc()
   if (present(xc_funcs)) then
     ixc_from_lib=libxc_functionals_ixc(xc_functionals=xc_funcs)
   else
     ixc_from_lib=libxc_functionals_ixc()
   end if
   if (ixc/=ixc_from_lib) then
     write(message, '(a,i0,2a,i0)')&
&     'The value of ixc specified in input, ixc = ',ixc,ch10,&
&     'differs from the one used to initialize the functional ',ixc_from_lib
     MSG_BUG(message)
   end if
 end if

!Check value of order
 if( (order<1.and.order/=-2).or.order>4)then
   write(message, '(a,i0)' )&
&   'The only allowed values for order are 1,2,-2 or 3, while it is found to be ',order
   MSG_BUG(message)
 end if

!Determine quantities available in input arguments
 has_gradient=.false.;has_laplacian=.false.;has_kden=.false.
 if (usegradient==1) then
   if (.not.present(grho2_updn)) then
     MSG_BUG('missing grho2_updn argument!')
   end if
   if (nvxcgrho>0) then
     if (.not.present(vxcgrho)) then
       MSG_BUG('missing vxcgrho argument!')
     end if
     has_gradient=.true.
   end if
 else if (nvxcgrho>0) then
   MSG_BUG('nvxcgrho>0 and usegradient=0!')
 end if
 if (uselaplacian==1) then
   if (.not.present(lrho_updn)) then
     MSG_BUG('missing lrho_updn argument!')
   end if
   if (nvxclrho>0) then
     if (.not.present(vxclrho)) then
       MSG_BUG('missing vxclrho argument!')
     end if
     has_laplacian=.true.
   end if
 else if (nvxclrho>0) then
   MSG_BUG('nvxclrho>0 and uselaplacian=0!')
 end if
 if (usekden==1) then
   if (.not.present(tau_updn)) then
     MSG_BUG('missing tau_updn argument!')
   end if
   if (nvxctau>0) then
     if (.not.present(vxctau)) then
       MSG_BUG('missing vxctau argument!')
     end if
     has_kden=.true.
   end if
 else if (nvxctau>0) then
   MSG_BUG('nvxctau>0 and usekden=0!')
 end if
 if (abs(order)>=2) then
   if (.not.present(dvxc)) then
     message='order>=2 needs argument dvxc!'
     MSG_BUG(message)
   else if (ndvxc==0) then
     message='order>=2 needs ndvxc>0!'
     MSG_BUG(message)
   end if
 end if
 if (abs(order)>=3) then
   if (.not.present(d2vxc)) then
     message='order>=3 needs argument d2vxc!'
     MSG_BUG(message)
   else if (nd2vxc==0) then
     message='order>=3 needs nd2vxc>0!'
     MSG_BUG(message)
   end if
 end if

!Determine quantities needed by XC functional
 if (present(xc_funcs)) then
   call size_dvxc(ixc,order,nspden,usegradient=need_gradient,&
&     uselaplacian=need_laplacian,usekden=need_kden,&
&     nvxcgrho=need_nvxcgrho,nvxclrho=need_nvxclrho,&
&     nvxctau=need_nvxctau,ndvxc=need_ndvxc,nd2vxc=need_nd2vxc,&
&     xc_funcs=xc_funcs)
 else
   call size_dvxc(ixc,order,nspden,usegradient=need_gradient,&
&     uselaplacian=need_laplacian,usekden=need_kden,&
&     nvxcgrho=need_nvxcgrho,nvxclrho=need_nvxclrho,&
&     nvxctau=need_nvxctau,ndvxc=need_ndvxc,nd2vxc=need_nd2vxc)
 end if
 if ((has_gradient.and.need_gradient>usegradient).or.&
&    (has_laplacian.and.need_laplacian>uselaplacian).or.&
&    (has_kden.and.need_kden>usekden)) then
   write(message, '(3a)' )&
&    'one of the arguments usegradient/uselaplacian/usesekden',ch10,&
&    'doesnt match the requirements of the XC functional!'
   MSG_BUG(message)
 end if
 if ((has_gradient.and.need_nvxcgrho>nvxcgrho).or.&
&    (has_laplacian.and.need_nvxclrho>nvxclrho).or.&
&    (has_kden.and.need_nvxctau>nvxctau).or.&
&    need_ndvxc>ndvxc.or.need_nd2vxc>nd2vxc) then
   write(message, '(3a)' )&
&    'one of the arguments nvxcgrho/nvxclrho/nvxctau/ndvxc/nd2vxc',ch10,&
&    'doesnt match the requirements of the XC functional!'
   MSG_BUG(message)
 end if
 if (abs(order)>1.and.(need_laplacian==1.or.need_kden==1)) then
   message='Derivatives of XC potential are not available in mGGA!'
   MSG_BUG(message)
 end if

!Check other optional arguments
 if (my_exexch/=0.and.usegradient==0) then
   message='exexch argument only valid for GGA!'
   MSG_BUG(message)
 end if
 if(ixc==50) then
   if(.not.(present(el_temp)).or.(.not.present(fxcT)))then
     message = 'el_temp or fxcT is not present but are needed for IIT XC functional.'
     MSG_BUG(message)
   end if
   if (size(fxcT)/=npts) then
     MSG_BUG('fxcT size must be npts!')
   end if
 end if

! =================================================
! ==  Intermediate quantities computation        ==
! =================================================

!If needed, compute rhotot and rs
 if (ixc==1.or.ixc==2.or.ixc==3.or.ixc==4.or.ixc==5.or.&
&    ixc==6.or.ixc==21.or.ixc==22.or.ixc==50) then
   ABI_ALLOCATE(rhotot,(npts))
   ABI_ALLOCATE(rspts,(npts))
   if(nspden==1)then
     rhotot(:)=two*rho_updn(:,1)
   else
     rhotot(:)=rho_updn(:,1)+rho_updn(:,2)
   end if
   call invcb(rhotot,rspts,npts)
   rspts(:)=rsfac*rspts(:)
 end if

!If needed, compute zeta
 if (ixc==1.or.ixc==21.or.ixc==22) then
   ABI_ALLOCATE(zeta,(npts))
   if(nspden==1)then
     zeta(:)=zero
   else
     zeta(:)=two*rho_updn(:,1)/rhotot(:)-one
   end if
 end if

! =================================================
! ==  XC energy, potentiel, ... computation      ==
! =================================================

!>>>>> No exchange-correlation
 if (ixc==0.or.ixc==40) then
   exc=zero ; vxcrho=zero
   if (present(dvxc).and.ndvxc>0) dvxc(:,:)=zero
   if (present(d2vxc).and.nd2vxc>0) d2vxc(:,:)=zero
   if (present(vxcgrho).and.nvxcgrho>0) vxcgrho(:,:)=zero
   if (present(vxclrho).and.nvxclrho>0) vxclrho(:,:)=zero
   if (present(vxctau).and.nvxctau>0) vxctau(:,:)=zero

!>>>>> New Teter fit (4/93) to Ceperley-Alder data, with spin-pol option
 else if (ixc==1 .or. ixc==21 .or. ixc==22) then
!  new Teter fit (4/93) to Ceperley-Alder data, with spin-pol option
   if (order**2 <= 1) then
     call xcspol(exc,npts,nspden,order,rspts,vxcrho,zeta,ndvxc)
   else
     call xcspol(exc,npts,nspden,order,rspts,vxcrho,zeta,ndvxc,dvxc)
   end if

!>>>>> Perdew-Zunger fit to Ceperly-Alder data (no spin-pol)
 else if (ixc==2) then
   if (order**2 <= 1) then
     call xcpzca(exc,npts,order,rhotot,rspts,vxcrho(:,1))
   else
     call xcpzca(exc,npts,order,rhotot,rspts,vxcrho(:,1),dvxc)
   end if

!>>>>> Teter fit (4/91) to Ceperley-Alder values (no spin-pol)
 else if (ixc==3) then
   if (order**2 <= 1) then
     call xctetr(exc,npts,order,rhotot,rspts,vxcrho(:,1))
   else if (order == 2) then
     call xctetr(exc,npts,order,rhotot,rspts,vxcrho(:,1),dvxc=dvxc)
   else if (order == 3) then
     call xctetr(exc,npts,order,rhotot,rspts,vxcrho(:,1),d2vxc=d2vxc,dvxc=dvxc)
   end if

!>>>>> Wigner xc (no spin-pol)
 else if (ixc==4) then
   if (order**2 <= 1) then
     call xcwign(exc,npts,order,rspts,vxcrho(:,1))
   else
     call xcwign(exc,npts,order,rspts,vxcrho(:,1),dvxc)
   end if

!>>>>>  Hedin-Lundqvist xc (no spin-pol)
 else if (ixc==5) then
   if (order**2 <= 1) then
     call xchelu(exc,npts,order,rspts,vxcrho(:,1))
   else
     call xchelu(exc,npts,order,rspts,vxcrho(:,1),dvxc)
   end if

!>>>>> X-alpha (no spin-pol)
 else if (ixc==6) then
   if (order**2 <= 1) then
     call xcxalp(exc,npts,order,rspts,vxcrho(:,1))
   else
     call xcxalp(exc,npts,order,rspts,vxcrho(:,1),dvxc)
   end if

!>>>>> PBE and alternatives
 else if (((ixc>=7.and.ixc<=15).or.(ixc>=23.and.ixc<=24)).and.ixc/=10.and.ixc/=13) then
!  Perdew-Wang LSD is coded in Perdew-Burke-Ernzerhof GGA, with optpbe=1
   if(ixc==7)optpbe=1
!  x-only part of Perdew-Wang
   if(ixc==8)optpbe=-1
!  Exchange + RPA correlation from Perdew-Wang
   if(ixc==9)optpbe=3
!  Perdew-Burke-Ernzerhof GGA
   if(ixc==11)optpbe=2
!  x-only part of PBE
   if(ixc==12)optpbe=-2
!  C09x exchange of V. R. Cooper
   if(ixc==24)optpbe=-4
!  revPBE of Zhang and Yang
   if(ixc==14)optpbe=5
!  RPBE of Hammer, Hansen and Norskov
   if(ixc==15)optpbe=6
!  Wu and Cohen
   if(ixc==23)optpbe=7
   if (ixc >=7.and.ixc<=9) then
     if (order**2 <= 1) then
       call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc)
     else if (order /=3) then
       call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,dvxci=dvxc)
     else if (order ==3) then
       call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,d2vxci=d2vxc,dvxci=dvxc)
     end if
   else if ((ixc >= 11 .and. ixc <= 15) .or. (ixc>=23 .and. ixc<=24)) then
     if (order**2 <= 1) then
       call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,&
&       dvxcdgr=vxcgrho,exexch=my_exexch,grho2_updn=grho2_updn)
     else if (order /=3) then
       if(ixc==12 .or. ixc==24)then
         call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,&
&         dvxcdgr=vxcgrho,dvxci=dvxc,grho2_updn=grho2_updn)
       else if(ixc/=12 .or. ixc/=24) then
         call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,&
&         dvxcdgr=vxcgrho,dvxci=dvxc,grho2_updn=grho2_updn)
       end if
     else if (order ==3) then
       if(ixc==12 .or. ixc==24)then
         call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,&
&         d2vxci=d2vxc,dvxcdgr=vxcgrho,dvxci=dvxc,grho2_updn=grho2_updn)
       else if(ixc/=12 .or. ixc/=24) then
         call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,&
&         d2vxci=d2vxc,dvxcdgr=vxcgrho,dvxci=dvxc,grho2_updn=grho2_updn)
       end if
     end if
   end if

!>>>>> RPA correlation from Perdew-Wang
 else if (ixc==10) then
   if (order**2 <= 1) then
     ABI_ALLOCATE(exci_rpa,(npts))
     ABI_ALLOCATE(vxci_rpa,(npts,2))
     optpbe=3
     call xcpbe(exci_rpa,npts,nspden,optpbe,order,rho_updn,vxci_rpa,ndvxc,nd2vxc)
     optpbe=1
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc)
     exc(:)=exc(:)-exci_rpa(:)
!    PMA: second index of vxcrho is nspden while that of rpa is 2 they can mismatch
     vxcrho(:,1:min(nspden,2))=vxcrho(:,1:min(nspden,2))-vxci_rpa(:,1:min(nspden,2))
     ABI_DEALLOCATE(exci_rpa)
     ABI_DEALLOCATE(vxci_rpa)
   else if (order /=3) then
     ABI_ALLOCATE(exci_rpa,(npts))
     ABI_ALLOCATE(vxci_rpa,(npts,2))
     optpbe=3
     call xcpbe(exci_rpa,npts,nspden,optpbe,order,rho_updn,vxci_rpa,ndvxc,nd2vxc,dvxci=dvxc)
     optpbe=1
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,dvxci=dvxc)
     exc(:)=exc(:)-exci_rpa(:)
     vxcrho(:,:)=vxcrho(:,:)-vxci_rpa(:,:)
     ABI_DEALLOCATE(exci_rpa)
     ABI_DEALLOCATE(vxci_rpa)
   else if (order ==3) then
     ABI_ALLOCATE(exci_rpa,(npts))
     ABI_ALLOCATE(vxci_rpa,(npts,2))
     optpbe=3
     call xcpbe(exci_rpa,npts,nspden,optpbe,order,rho_updn,vxci_rpa,ndvxc,nd2vxc,&
&     d2vxci=d2vxc,dvxci=dvxc)
     optpbe=1
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,&
&     d2vxci=d2vxc,dvxci=dvxc)
     exc(:)=exc(:)-exci_rpa(:)
     vxcrho(:,:)=vxcrho(:,:)-vxci_rpa(:,:)
     ABI_DEALLOCATE(exci_rpa)
     ABI_DEALLOCATE(vxci_rpa)
   end if

!>>>>> LDA xc energy like ixc==7, and Leeuwen-Baerends GGA xc potential
 else if(ixc==13) then
   if (order**2 <= 1) then
     optpbe=1
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc)
     call xclb(grho2_updn,npts,nspden,rho_updn,vxcrho)
   else if (order /=3) then
     optpbe=1
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,dvxci=dvxc)
     call xclb(grho2_updn,npts,nspden,rho_updn,vxcrho)
   else if (order ==3) then
     optpbe=1
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,d2vxci=d2vxc,dvxci=dvxc)
     call xclb(grho2_updn,npts,nspden,rho_updn,vxcrho)
   end if

!>>>>> HTCH93, HTCH120, HTCH107, HTCH147
 else if(ixc==16 .or. ixc==17 .or. ixc==26 .or. ixc==27) then
   call xchcth(vxcgrho,exc,grho2_updn,ixc,npts,nspden,order,rho_updn,vxcrho)

!>>>>> Only for test purpose (test various part of MGGA implementation)
 else if(ixc==31 .or. ixc==32 .or. ixc==33 .or. ixc==34) then
   exc(:)=zero
   vxcrho(:,:)=zero
   vxcgrho(:,:)=zero
   vxclrho(:,:)=zero
   vxctau(:,:)=zero

!>>>>> Perdew-Wang LSD is coded in Perdew-Burke-Ernzerhof GGA, with optpbe=1
   optpbe=1
   select case(ixc)
   case (31)
     alpha=1.00d0-(1.00d0/1.01d0)
!      Compute first LDA XC (exc,vxc) and then add fake MGGA XC (exc,vxc)
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc)
     if (nspden==1) then
!        it should be : exc_tot= exc_spin up + exc_spin down = 2*exc_spin up but this applies to tau and rho (so it cancels)
       exc(:)=exc(:)+alpha*tau_updn(:,1)/rho_updn(:,1)
     else
       do ispden=1,nspden
         exc(:)=exc(:)+alpha*tau_updn(:,ispden)/(rho_updn(:,1)+rho_updn(:,2))
       end do
     end if
     vxctau(:,:)=alpha
   case (32)
     alpha=0.01d0
!      Compute first LDA XC (exc,vxc) and then add fake MGGA XC (exc,vxc)
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc)
     if (nspden==1) then
       exc(:)=exc(:)+2.0d0*alpha*lrho_updn(:,1)
       vxcrho(:,1) =vxcrho(:,1)+2.0d0*alpha*lrho_updn(:,1)
       vxclrho(:,1)=alpha*2.0d0*rho_updn(:,1)
     else
       do ispden=1,nspden
         exc(:)=exc(:)+alpha*lrho_updn(:,ispden)
         vxcrho(:,ispden) =vxcrho(:,ispden)+alpha*(lrho_updn(:,1)+lrho_updn(:,2))
         vxclrho(:,ispden)=alpha*(rho_updn(:,1)+rho_updn(:,2))
       end do
     end if
   case (33)
     alpha=-0.010d0
!      Compute first LDA XC (exc,vxc) and then add fake MGGA XC (exc,vxc)
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc)
     if (nspden==1) then
!        it should be : exc_tot= exc_spin up + exc_spin down = 2*exc_spin up but this applies to grho2 and rho
!        (for grho2 it is a factor 4 to have total energy and for rho it is just a factor 2. So we end with factor 2 only)
       exc(:)=exc(:)+alpha*2.0d0*grho2_updn(:,1)/rho_updn(:,1)
       if(nvxcgrho==2)vxcgrho(:,1:2)=2.0d0*alpha
       if(nvxcgrho==3)vxcgrho(:,3)=2.0d0*alpha
     else
       exc(:)=exc(:)+alpha*grho2_updn(:,3)/(rho_updn(:,1)+rho_updn(:,2))
       if(nvxcgrho==2)vxcgrho(:,1:2)=2.0d0*alpha
       if(nvxcgrho==3)vxcgrho(:,3)=2.0d0*alpha
     end if
   case (34)
     alpha=-0.010d0
!      Compute first LDA XC (exc,vxc) and then add fake MGGA XC (exc,vxc)
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc)
     if (nspden==1) then
       exc(:)=exc(:)+16.0d0*alpha*tau_updn(:,1)
       vxcrho(:,1)=vxcrho(:,1)+16.0d0*alpha*tau_updn(:,1)
       vxctau(:,1)=16.0d0*alpha*rho_updn(:,1)
     else
       do ispden=1,nspden
         exc(:)=exc(:)+8.0d0*alpha*tau_updn(:,ispden)
         vxcrho(:,ispden)=vxcrho(:,ispden)+8.0d0*alpha*(tau_updn(:,1)+tau_updn(:,2))
         vxctau(:,ispden)=8.0d0*alpha*(rho_updn(:,1)+rho_updn(:,2))
       end do
     end if
   end select

!>>>>> Hybrid PBE0 (1/4 and 1/3)
 else if(ixc>=41.and.ixc<=42) then
!  Requires to evaluate exchange-correlation with PBE (optpbe=2)
!  minus hyb_mixing*exchange with PBE (optpbe=-2)
   ndvxc_x=8
   ABI_ALLOCATE(exc_x,(npts))
   ABI_ALLOCATE(vxcrho_x,(npts,nspden))
   ABI_ALLOCATE(vxcgrho_x,(npts,nvxcgrho))
   exc_x=zero;vxcrho_x=zero;vxcgrho_x=zero
   if (order**2 <= 1) then
     optpbe=2 !PBE exchange correlation
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,&
&     dvxcdgr=vxcgrho,exexch=my_exexch,grho2_updn=grho2_updn)
     optpbe=-2 !PBE exchange-only
     call xcpbe(exc_x,npts,nspden,optpbe,order,rho_updn,vxcrho_x,ndvxc,nd2vxc,&
&     dvxcdgr=vxcgrho_x,exexch=my_exexch,grho2_updn=grho2_updn)
     exc=exc-exc_x*my_hyb_mixing
     vxcrho=vxcrho-vxcrho_x*my_hyb_mixing
     vxcgrho=vxcgrho-vxcgrho_x*my_hyb_mixing
   else if (order /=3) then
     ABI_ALLOCATE(dvxc_x,(npts,ndvxc_x))
     optpbe=2 !PBE exchange correlation
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,&
     dvxcdgr=vxcgrho,dvxci=dvxc,grho2_updn=grho2_updn)
     optpbe=-2 !PBE exchange-only
     call xcpbe(exc_x,npts,nspden,optpbe,order,rho_updn,vxcrho_x,ndvxc_x,nd2vxc,&
&     dvxcdgr=vxcgrho_x,dvxci=dvxc_x,grho2_updn=grho2_updn)
     exc=exc-exc_x*my_hyb_mixing
     vxcrho=vxcrho-vxcrho_x*my_hyb_mixing
     vxcgrho=vxcgrho-vxcgrho_x*my_hyb_mixing
     dvxc(:,1:ndvxc_x)=dvxc(:,1:ndvxc_x)-dvxc_x(:,1:ndvxc_x)*my_hyb_mixing
     ABI_DEALLOCATE(dvxc_x)
   else if (order ==3) then
!    The size of exchange-correlation with PBE (optpbe=2)
!    is the one which defines the size for ndvxc.
     ABI_ALLOCATE(dvxc_x,(npts,ndvxc_x))
     ABI_ALLOCATE(d2vxc_x,(npts,nd2vxc))
     optpbe=2 !PBE exchange correlation
     call xcpbe(exc,npts,nspden,optpbe,order,rho_updn,vxcrho,ndvxc,nd2vxc,&
&     d2vxci=d2vxc,dvxcdgr=vxcgrho,dvxci=dvxc,grho2_updn=grho2_updn)
     optpbe=-2 !PBE exchange-only
     call xcpbe(exc_x,npts,nspden,optpbe,order,rho_updn,vxcrho_x,ndvxc_x,nd2vxc,&
&     d2vxci=d2vxc_x,dvxcdgr=vxcgrho_x,dvxci=dvxc_x,grho2_updn=grho2_updn)
     exc=exc-exc_x*my_hyb_mixing
     vxcrho=vxcrho-vxcrho_x*my_hyb_mixing
     vxcgrho=vxcgrho-vxcgrho_x*my_hyb_mixing
     d2vxc=d2vxc-d2vxc_x*my_hyb_mixing
     dvxc(:,1:ndvxc_x)=dvxc(:,1:ndvxc_x)-dvxc_x(:,1:ndvxc_x)*my_hyb_mixing
     ABI_DEALLOCATE(dvxc_x)
     ABI_DEALLOCATE(d2vxc_x)
   end if
   ABI_DEALLOCATE(exc_x)
   ABI_DEALLOCATE(vxcrho_x)
   ABI_DEALLOCATE(vxcgrho_x)

!>>>>> Ichimaru,Iyetomi,Tanaka,  XC at finite temp (e- gaz)
 else if (ixc==50) then
   if (order**2 <= 1) then
     call xciit(exc,fxcT,npts,order,rspts,el_temp,vxcrho(:,1))
   else
     call xciit(exc,fxcT,npts,order,rspts,el_temp,vxcrho(:,1),dvxc)
   end if

!>>>>> GGA counterpart of the B3LYP functional
 else if(ixc==1402000) then
!  Requires to evaluate exchange-correlation
!  with 5/4 B3LYP - 1/4 B3LYPc, where
!  B3LYPc = (0.19 Ec VWN3 + 0.81 Ec LYP)

!  First evaluate B3LYP.
   if(present(xc_funcs))then
     if (abs(order)==1) then
       call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,&
&       vxcrho,grho2=grho2_updn,vxcgr=vxcgrho,xc_functionals=xc_funcs)
     else if (abs(order)==2) then
       call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,&
&       vxcrho,grho2=grho2_updn,vxcgr=vxcgrho,dvxc=dvxc,xc_functionals=xc_funcs)
     else
       call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,&
&       vxcrho,grho2=grho2_updn,vxcgr=vxcgrho,dvxc=dvxc,d2vxc=d2vxc,xc_functionals=xc_funcs)
     end if
   else
     if (abs(order)==1) then
       call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,&
&       vxcrho,grho2=grho2_updn,vxcgr=vxcgrho)
     else if (abs(order)==2) then
       call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,&
&       vxcrho,grho2=grho2_updn,vxcgr=vxcgrho,dvxc=dvxc)
     else
       call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,&
&       vxcrho,grho2=grho2_updn,vxcgr=vxcgrho,dvxc=dvxc,d2vxc=d2vxc)
     end if
   end if

!  Then renormalize B3LYP and subtract VWN3 contribution
   ABI_ALLOCATE(exc_c,(npts))
   ABI_ALLOCATE(vxcrho_c,(npts,nspden))
   if(order**2>1)then
     ABI_ALLOCATE(dvxc_c,(npts,ndvxc))
   end if
   if(order**2>4)then
     ABI_ALLOCATE(d2vxc_c,(npts,nd2vxc))
   end if
   exc_c=zero;vxcrho_c=zero
   call libxc_functionals_init(-30,nspden,xc_functionals=xc_funcs_vwn3)
   if (order**2 <= 1) then
     call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc_c,&
&     vxcrho_c,xc_functionals=xc_funcs_vwn3)
   elseif (order**2 <= 4) then
     dvxc_c=zero
     call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc_c,&
&     vxcrho_c,dvxc=dvxc_c,xc_functionals=xc_funcs_vwn3)
   else
     dvxc_c=zero
     d2vxc_c=zero
     call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc_c,&
&     vxcrho_c,dvxc=dvxc_c,d2vxc=d2vxc,xc_functionals=xc_funcs_vwn3)
   end if
   exc=1.25d0*exc-quarter*0.19d0*exc_c
   vxcrho=1.25d0*vxcrho-quarter*0.19d0*vxcrho_c
   if(order**2>1)dvxc=1.25d0*dvxc-quarter*0.19d0*dvxc_c
   if(order**2>4)d2vxc=1.25d0*d2vxc-quarter*0.19d0*d2vxc_c
   call libxc_functionals_end(xc_functionals=xc_funcs_vwn3)

!  Then subtract LYP contribution
   call libxc_functionals_init(-131,nspden,xc_functionals=xc_funcs_lyp)
   if (order**2 <= 1) then
     call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc_c,&
&     vxcrho_c,grho2=grho2_updn,vxcgr=vxcgrho,xc_functionals=xc_funcs_lyp)
   elseif (order**2 <= 4) then
     dvxc_c=zero
     call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc_c,&
&     vxcrho_c,grho2=grho2_updn,vxcgr=vxcgrho,dvxc=dvxc_c,xc_functionals=xc_funcs_lyp)
   else
     dvxc_c=zero
     d2vxc_c=zero
     call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc_c,&
&     vxcrho_c,grho2=grho2_updn,vxcgr=vxcgrho,dvxc=dvxc_c,d2vxc=d2vxc,xc_functionals=xc_funcs_lyp)
   end if
   exc=exc-quarter*0.81d0*exc_c
   vxcrho=vxcrho-quarter*0.81d0*vxcrho_c
   if(order**2>1)dvxc=dvxc-quarter*0.81d0*dvxc_c
   if(order**2>4)d2vxc=d2vxc-quarter*0.81d0*d2vxc_c
   call libxc_functionals_end(xc_functionals=xc_funcs_lyp)

   ABI_DEALLOCATE(exc_c)
   ABI_DEALLOCATE(vxcrho_c)
   if(allocated(dvxc_c))then
     ABI_DEALLOCATE(dvxc_c)
   end if
   if(allocated(d2vxc_c))then
     ABI_DEALLOCATE(d2vxc_c)
   end if

!>>>>> All libXC functionals
 else if( ixc<0 ) then

!  ===== meta-GGA =====
   if (need_laplacian==1.or.need_kden==1) then
     if (need_laplacian==1.and.need_kden==1) then
       if (present(xc_funcs)) then
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           grho2=grho2_updn,vxcgr=vxcgrho,&
&           lrho=lrho_updn,vxclrho=vxclrho,&
&           tau=tau_updn,vxctau=vxctau,&
&           xc_functionals=xc_funcs)
       else
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           grho2=grho2_updn,vxcgr=vxcgrho,&
&           lrho=lrho_updn,vxclrho=vxclrho,&
&           tau=tau_updn,vxctau=vxctau)
       end if
     else if (need_laplacian==1) then
       if (present(xc_funcs)) then
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           grho2=grho2_updn,vxcgr=vxcgrho,&
&           lrho=lrho_updn,vxclrho=vxclrho,&
&           xc_functionals=xc_funcs)
       else
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           grho2=grho2_updn,vxcgr=vxcgrho,&
&           lrho=lrho_updn,vxclrho=vxclrho)
       end if
     else if (need_kden==1) then
       if (present(xc_funcs)) then
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           grho2=grho2_updn,vxcgr=vxcgrho,&
&           tau=tau_updn,vxctau=vxctau,&
&           xc_functionals=xc_funcs)
       else
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           grho2=grho2_updn,vxcgr=vxcgrho,&
&           tau=tau_updn,vxctau=vxctau)
       end if
     end if
     !Some meta-GGAs can only be used with a LDA correlation (see doc)
     ixc1=(-ixc)/1000;ixc2=(-ixc)-ixc1*1000
     if (ixc1==206 .or. ixc1==207 .or. ixc1==208 .or. ixc1==209 .or. &
&        ixc2==206 .or. ixc2==207 .or. ixc2==208 .or. ixc2==209    )then
       if (present(vxcgrho)) vxcgrho(:,:)=zero
       if (present(vxclrho)) vxclrho(:,:)=zero
       if (present(vxctau)) vxctau(:,:)=zero
     end if

!  ===== GGA =====
   else if (need_gradient==1) then
     if (abs(order)<=1) then
       if (present(xc_funcs)) then
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           grho2=grho2_updn,vxcgr=vxcgrho,&
&           xc_functionals=xc_funcs)
       else
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           grho2=grho2_updn,vxcgr=vxcgrho)
       end if
     else
       if (present(xc_funcs)) then
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           grho2=grho2_updn,vxcgr=vxcgrho,dvxc=dvxc,&
&           xc_functionals=xc_funcs)
       else
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           grho2=grho2_updn,vxcgr=vxcgrho,dvxc=dvxc)
       end if
     end if

!  ===== LDA =====
   else
     if (abs(order)<=1) then
       if (present(xc_funcs)) then
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           xc_functionals=xc_funcs)
       else
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho)
       end if
     else if (abs(order)<=2) then
       if (present(xc_funcs)) then
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           dvxc=dvxc,xc_functionals=xc_funcs)
       else
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           dvxc=dvxc)
       end if
     else if (abs(order)<=3) then
       if (present(xc_funcs)) then
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           dvxc=dvxc,d2vxc=d2vxc,xc_functionals=xc_funcs)
       else
         call libxc_functionals_getvxc(ndvxc,nd2vxc,npts,nspden,order,rho_updn,exc,vxcrho,&
&           dvxc=dvxc,d2vxc=d2vxc)
       end if
     end if

   end if ! mGGA, GGA, LDA
 end if ! libXC

! =================================================
! ==              Finalization                   ==
! =================================================
!Deallocate arrays
 if(allocated(rhotot)) then
   ABI_DEALLOCATE(rhotot)
 end if
 if(allocated(rspts)) then
   ABI_DEALLOCATE(rspts)
 end if
 if(allocated(zeta)) then
   ABI_DEALLOCATE(zeta)
 end if

end subroutine drivexc
!!***

end module m_drivexc
!!***
