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
!!  Copyright (C) 2012-2018 ABINIT group (MT, MJV, CE, TD)
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
 use m_profiling_abi
 use m_errors
 use libxc_functionals

 implicit none

 private
!!***

 public :: drivexc_main    ! Driver of XC functionals. Optionally, deliver the XC kernel, or even the derivative
 public :: echo_xc_name    ! Write to log and output the xc functional which will be used for this dataset
 public :: check_kxc       ! Given a XC functional (defined by ixc), check if Kxc (dVxc/drho) is avalaible.
 public :: size_dvxc       ! Give the size of the array dvxc(npts,ndvxc) and the second dimension of the d2vxc(npts,nd2vxc)
!!***

contains
!!***

!!****f* ABINIT/drivexc_main
!! NAME
!! drivexc_main
!!
!! FUNCTION
!! Driver of XC functionals. Optionally, deliver the XC kernel, or even the derivative
!! of the XC kernel (the third derivative of the XC energy)
!!
!! INPUTS
!!  ixc=index of the XC functional
!!  mgga=1 if functional is metaGGA
!!  ndvxc=size of dvxc(npts,ndvxc)
!!  nd2vxc=size of d2vxc(npts,nd2vxc)
!!  ngr2=size of grho2(npts,ngr2)
!!  npts=number of real space points on which the density (and its gradients) is provided
!!  nspden=number of spin-density components (1 or 2)
!!  nvxcgrho=size of vxcgrho(npts,nvxcgrho)
!!  order=gives the maximal derivative of Exc computed.
!!    1=usual value (return exc and vxc)
!!    2=also computes the kernel (return exc,vxc,kxc)
!!   -2=like 2, except (to be described)
!!    3=also computes the derivative of the kernel (return exc,vxc,kxc,k3xc)
!!  rho(npts,nspden)=the spin-up and spin-down densities
!!    If nspden=1, only the spin-up density must be given.
!!    In the calling routine, the spin-down density must
!!    be equal to the spin-up density5
!!    and both are half the total density.
!!    If nspden=2, the spin-up and spin-down density must be given
!!  xclevel= XC functional level
!!  === Optional input arguments ===
!!  [el_temp]= electronic temperature (to be used for finite temperature XC functionals)
!!  [exexch]=choice of <<<local>>> exact exchange. Active if exexch=3 (only for GGA, and NOT for libxc)
!!  [grho2(npts,ngr2)]=the square of the gradients of
!!    spin-up, spin-down, and total density
!!    If nspden=1, only the square of the gradient of the spin-up density
!!     must be given. In the calling routine, the square of the gradient
!!     of the spin-down density
!!     must be equal to the square of the gradient of the spin-up density,
!!     and both must be equal to one-quarter of the square of the
!!     gradient of the total density.
!!    If nspden=2, the square of the gradients of
!!     spin-up, spin-down, and total density must be given.
!!     Note that the square of the gradient of the total
!!     density is usually NOT related to the square of the
!!     gradient of the spin-up and spin-down densities,
!!     because the gradients are not usually aligned.
!!     This is not the case when nspden=1.
!!  [hyb_mixing]= mixing parameter for the native PBEx functionals (ixc=41 and 42)
!!  [lrho(npts,nspden)]=the Laplacian of spin-up and spin-down densities
!!    If nspden=1, only the spin-up Laplacian density must be given.
!!    In the calling routine, the spin-down Laplacian density must
!!    be equal to the spin-up Laplacian density,
!!    and both are half the total Laplacian density.
!!    If nspden=2, the Laplacian of spin-up and spin-down densities must be given
!!  [tau(npts,nspden)]=the spin-up and spin-down kinetic energy densities
!!    If nspden=1, only the spin-up kinetic energy density must be given.
!!    In the calling routine, the spin-down kinetic energy density must
!!    be equal to the spin-up kinetic energy density,
!!    and both are half the total kinetic energy density.
!!    If nspden=2, the spin-up and spin-down kinetic energy densities must be given
!!  [xc_funcs(2)]= <type(libxc_functional_type)>, optional : libxc XC functionals.
!!  [xc_tb09_c]=c parameter for the mgga TB09 functional, within libxc

!! OUTPUT
!!  exc(npts)=exchange-correlation energy density (hartree)
!!  vxcrho(npts,nspden)= (d($\rho$*exc)/d($\rho_up$)) (hartree)
!!                  and  (d($\rho$*exc)/d($\rho_down$)) (hartree)
!!  === Optional output arguments ===
!!  [dvxc]=partial second derivatives of the xc energy, only if abs(order)>1
!!   In case of local energy functional (option=1,-1 or 3):
!!    dvxc(npts,1+nspden)=              (Hartree*bohr^3)
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
!!  [d2vxc]=second derivative of the XC potential=3rd order derivative of energy, only if abs(order)>2
!!    (only available for LDA and nspden=1)
!!    if nspden=1 d2vxc(npts,1)=second derivative of the XC potential=3rd order derivative of energy
!!    if nspden=2 d2vxc(npts,1), d2vxc(npts,2), d2vxc(npts,3), d2vxc(npts,4) (3rd derivative of energy)
!!  [fxcT(npts)]=XC free energy of the electron gaz at finite temperature (to be used for plasma systems)
!!  [vxcgrho(npts,3)]= 1/$|grad \rho_up|$ (d($\rho$*exc)/d($|grad \rho_up|$)) (hartree)
!!                     1/$|grad \rho_dn|$ (d($\rho$*exc)/d($|grad \rho_dn|$)) (hartree)
!!                and  1/$|grad \rho|$ (d($\rho$*exc)/d($|grad \rho|$))       (hartree)
!!     (will be zero if a LDA functional is used)
!!  [vxclrho(npts,nspden)]=(only for meta-GGA, i.e. optional output)=
!!                       (d($\rho$*exc)/d($\lrho_up$))   (hartree)
!!                  and  (d($\rho$*exc)/d($\lrho_down$)) (hartree)
!!  [vxctau(npts,nspden)]=(only for meta-GGA, i.e. optional output)=
!!    derivative of XC energy density with respect to kinetic energy density (depsxcdtau).
!!                       (d($\rho$*exc)/d($\tau_up$))    (hartree)
!!                  and  (d($\rho$*exc)/d($\tau_down$))  (hartree)
!!
!! PARENTS
!!      m_pawxc,rhotoxc
!!
!! CHILDREN
!!      drivexc
!!
!! SOURCE

subroutine drivexc_main(exc,ixc,mgga,ndvxc,nd2vxc,ngr2,npts,nspden,nvxcgrho,order,rho,vxcrho,xclevel, &
&                       dvxc,d2vxc,el_temp,exexch,fxcT,grho2,& ! Optional arguments
&                       hyb_mixing,lrho,tau,vxcgrho,vxclrho,vxctau,xc_funcs,xc_tb09_c) ! Optional arguments


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'drivexc_main'
 use interfaces_41_xc_lowlevel
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,mgga,ndvxc,nd2vxc,ngr2,npts,nspden,nvxcgrho,order,xclevel
 integer,intent(in),optional :: exexch
 real(dp),intent(in),optional :: el_temp,hyb_mixing,xc_tb09_c
!arrays
 real(dp),intent(in) :: rho(npts,nspden)
 real(dp),intent(in),optional :: grho2(npts,ngr2),lrho(npts,nspden*mgga),tau(npts,nspden*mgga)
 real(dp),intent(out) :: exc(npts),vxcrho(npts,nspden)
 real(dp),intent(out),optional :: dvxc(npts,ndvxc),d2vxc(npts,nd2vxc),fxcT(:),vxcgrho(npts,nvxcgrho)
 real(dp),intent(out),optional :: vxclrho(npts,nspden*mgga),vxctau(npts,nspden*mgga)
 type(libxc_functional_type),intent(inout),optional :: xc_funcs(2)

!Local variables-------------------------------
!scalars
 real(dp) :: hyb_mixing_,xc_tb09_c_
 logical :: is_gga

!  *************************************************************************

!Checks input parameters
 if (mgga==1) then
   if (.not.present(lrho)) then
     MSG_BUG('lrho arg must be present in case of mGGA!')
   end if
   if (.not.present(tau)) then
     MSG_BUG('tau arg must be present in case of mGGA!')
   end if
   if (.not.present(vxclrho)) then
     MSG_BUG('vxclrho arg must be present in case of mGGA!')
   end if
   if (.not.present(vxctau)) then
     MSG_BUG('vxctau arg must be present in case of mGGA!')
   end if
 end if
 if (present(fxcT)) then
   if (.not.present(el_temp)) then
     MSG_BUG('el_temp arg must be present together with fxcT!')
   end if
 end if

 xc_tb09_c_=99.99_dp;if (present(xc_tb09_c)) xc_tb09_c_=xc_tb09_c

 if(ixc==41)hyb_mixing_=quarter
 if(ixc==42)hyb_mixing_=third
 if (present(hyb_mixing)) hyb_mixing_=hyb_mixing

!>>>>> All libXC functionals

 if (ixc<0) then

   if (present(xc_funcs))then
     is_gga=libxc_functionals_isgga(xc_functionals=xc_funcs)
   else
     is_gga=libxc_functionals_isgga()
   end if

   if (mgga==1) then
     if (abs(xc_tb09_c_-99.99_dp)>tol12) then
       if (present(xc_funcs)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         lrho_updn=lrho,vxclrho=vxclrho,tau_updn=tau,vxctau=vxctau, &
&         xc_funcs=xc_funcs,xc_tb09_c=xc_tb09_c_)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho, &
&         lrho_updn=lrho,vxclrho=vxclrho,tau_updn=tau,vxctau=vxctau, &
&         xc_tb09_c=xc_tb09_c_)
       end if
     else
       if (present(xc_funcs)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho, &
&         lrho_updn=lrho,vxclrho=vxclrho,tau_updn=tau,vxctau=vxctau, &
&         xc_funcs=xc_funcs)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho, &
&         lrho_updn=lrho,vxclrho=vxclrho,tau_updn=tau,vxctau=vxctau)
       end if
     end if
   else if (is_gga) then
     if (order**2<=1) then
       if (present(xc_funcs)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho,xc_funcs=xc_funcs)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho)
       end if
     else
       if (present(xc_funcs)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho,dvxc=dvxc,xc_funcs=xc_funcs)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho,dvxc=dvxc)
       end if
     end if
   else
     if (order**2<=1) then
       if (present(xc_funcs)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
         xc_funcs=xc_funcs)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho)
       end if
     else if (order**2<=4) then
       if (present(xc_funcs)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         dvxc=dvxc, xc_funcs=xc_funcs)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         dvxc=dvxc)
       end if
     else
       if (present(xc_funcs)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         dvxc=dvxc,d2vxc=d2vxc, xc_funcs=xc_funcs)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         dvxc=dvxc,d2vxc=d2vxc)
       end if
     end if
   end if

 else

!  Cases with gradient
   if (xclevel==2)then
     if (order**2<=1.or.ixc==16.or.ixc==17.or.ixc==26.or.ixc==27) then
       if (ixc/=13) then
         if (present(exexch)) then
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&           hyb_mixing=hyb_mixing_,grho2_updn=grho2,vxcgrho=vxcgrho, &
&           exexch=exexch)
         else
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&           hyb_mixing=hyb_mixing_,grho2_updn=grho2,vxcgrho=vxcgrho)
         end if
       else
         if (present(exexch)) then
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&           grho2_updn=grho2, &
&           exexch=exexch)
         else
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&           grho2_updn=grho2)
         end if
       end if
     else if (order/=3) then
       if (ixc/=13) then
         if (present(exexch)) then
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&           hyb_mixing=hyb_mixing_,dvxc=dvxc,grho2_updn=grho2,vxcgrho=vxcgrho, &
&           exexch=exexch)
         else
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&           hyb_mixing=hyb_mixing_,dvxc=dvxc,grho2_updn=grho2,vxcgrho=vxcgrho)
         end if
       else
         if (present(exexch)) then
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&           dvxc=dvxc,grho2_updn=grho2, &
&           exexch=exexch)
         else
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&           dvxc=dvxc,grho2_updn=grho2)
         end if
       end if
     else if (order==3) then
       if (ixc/=13) then
         if (present(exexch)) then
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&           dvxc=dvxc,d2vxc=d2vxc,grho2_updn=grho2,vxcgrho=vxcgrho, &
&           hyb_mixing=hyb_mixing_,exexch=exexch)
         else
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&           hyb_mixing=hyb_mixing_,dvxc=dvxc,d2vxc=d2vxc,grho2_updn=grho2,vxcgrho=vxcgrho)
         end if
       else
         if (present(exexch)) then
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&           dvxc=dvxc,d2vxc=d2vxc,grho2_updn=grho2, &
&           exexch=exexch)
         else
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&           dvxc=dvxc,d2vxc=d2vxc,grho2_updn=grho2)
         end if
       end if
     end if


!    Cases without gradient
   else
     if (order**2<=1) then
       if (ixc>=31.and.ixc<=34) then !fake mgga functionals for testing purpose only (based on LDA functional)
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho, &
&         lrho_updn=lrho,vxclrho=vxclrho,tau_updn=tau,vxctau=vxctau)
       else
         if (present(fxcT)) then
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho,fxcT=fxcT,el_temp=el_temp)
         else
           call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho)
         end if
       end if
     else if (order==3.and.(ixc==3.or.ixc>=7.and.ixc<=10)) then
       call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&       dvxc=dvxc,d2vxc=d2vxc)
     else
       if (present(fxcT)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         dvxc=dvxc,fxcT=fxcT,el_temp=el_temp)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         dvxc=dvxc)
       end if
     end if
   end if

 end if

end subroutine drivexc_main
!!***

!!****f* ABINIT/echo_xc_name
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

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'echo_xc_name'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
     citation = 'S. Goedecker, M. Teter, J. Huetter, PRB 54, 1703 (1996)'
   case (2)
     message = 'LDA: Perdew-Zunger-Ceperley-Alder - ixc=2'
     citation = 'J.P.Perdew and A.Zunger, PRB 23, 5048 (1981) '
   case (3)
     message = 'LDA: old Teter (4/91) fit to Ceperley-Alder data - ixc=3'
     citation = ''
   case (4)
     message = 'LDA: Wigner - ixc=4'
     citation = 'E.P.Wigner, Trans. Faraday Soc. 34, 67 (1938)'
   case (5)
     message = 'LDA: Hedin-Lundqvist - ixc=5'
     citation = 'L.Hedin and B.I.Lundqvist, J. Phys. C4, 2064 (1971)'
   case (6)
     message = 'LDA: "X-alpha" xc - ixc=6'
     citation = 'Slater J. C., Phys. Rev. 81, 385 (1951)'
   case (7)
     message = 'LDA: Perdew-Wang 92 LSD fit to Ceperley-Alder data - ixc=7'
     citation = 'J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)'
   case (8)
     message = 'LDA: Perdew-Wang 92 LSD , exchange-only - ixc=8'
     citation = 'J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)'
   case (9)
     message = 'LDA: Perdew-Wang 92 Ex+Ec_RPA  energy - ixc=9'
     citation = 'J.P.Perdew and Y.Wang, PRB 45, 13244 (1992)'
   case (10)
     message = 'LDA: RPA LSD energy (only the energy !!) - ixc=10'
     citation = ''
!      GGA
   case (11)
     message = 'GGA: Perdew-Burke-Ernzerhof functional - ixc=11'
     citation = 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)'
   case (12)
     message = 'GGA: x-only Perdew-Burke-Ernzerhof functional - ixc=12'
     citation = 'J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)'
   case (13)
     message = 'GGA: LDA (ixc==7) energy, and the xc _potential_ is given by van Leeuwen-Baerends GGA - ixc=13'
     citation = 'R. van Leeuwen and E. J. Baerends PRA 49, 2421 (1994)'
   case (14)
     message = 'GGA: revPBE functional - ixc=14'
     citation = 'Zhang and Yang, PRL 80, 890 (1998)'
   case (15)
     message = 'GGA: RPBE functional - ixc=15'
     citation = 'Hammer, L. B. Hansen, and J. K. Norskov, PRB 59, 7413 (1999)'
   case (16)
     message = 'GGA: HCTH93 functional - ixc=16'
     citation = 'F.A. Hamprecht, A.J. Cohen, D.J. Tozer, N.C. Handy, JCP 109, 6264 (1998)'
   case (17)
     message = 'GGA: HCTH120 functional - ixc=17'
     citation = 'A.D. Boese, N.L. Doltsinis, N.C. Handy, and M. Sprik, JCP 112, 1670 (1998)'
   case (23)
     message = 'GGA: Wu Cohen functional - ixc=23'
     citation = 'Z. Wu and R. E. Cohen, PRB 73, 235116 (2006)'
   case (24)
     message = 'GGA: C09x exchange functional - ixc=24'
     citation = 'Valentino R. Cooper, PRB 81, 161104(R) (2010)'
   case (26)
     message = 'GGA: HCTH147 functional - ixc=26'
     citation = 'A.D. Boese, N.L. Doltsinis, N.C. Handy, and M. Sprik, JCP 112, 1670 (1998)'
   case (27)
     message = 'GGA: HCTH407 functional - ixc=27'
     citation = 'A.D. Boese, and N.C. Handy, JCP 114, 5497 (2001)'
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
     citation = 'Ichimaru S., Iyetomi H., Tanaka S., Phys. Rep. 149, 91-205 (1987) '
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

!!****f* ABINIT/check_kxc
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

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'check_kxc'
!End of the abilint section

 implicit none

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

!!****f* ABINIT/size_dvxc
!! NAME
!! size_dvxc
!!
!! FUNCTION
!! Give the size of the array dvxc(npts,ndvxc) and the second dimension of the d2vxc(npts,nd2vxc)
!! needed for the allocations depending on the routine which is called from the drivexc routine
!!
!! INPUTS
!!  [add_tfw]= optional flag controling the addition of Weiszacker gradient correction to Thomas-Fermi XC energy
!!  ixc= choice of exchange-correlation scheme
!!  order=gives the maximal derivative of Exc computed.
!!    1=usual value (return exc and vxc)
!!    2=also computes the kernel (return exc,vxc,kxc)
!!   -2=like 2, except (to be described)
!!    3=also computes the derivative of the kernel (return exc,vxc,kxc,k3xc)
!!  [xc_funcs(2)]= <type(libxc_functional_type)>
!!
!! OUTPUT
!!  ndvxc size of the array dvxc(npts,ndvxc) for allocation
!!  ngr2 size of the array grho2_updn(npts,ngr2) for allocation
!!  nd2vxc size of the array d2vxc(npts,nd2vxc) for allocation
!!  nvxcdgr size of the array dvxcdgr(npts,nvxcdgr) for allocation
!!
!! PARENTS
!!      m_pawxc,rhotoxc
!!
!! CHILDREN
!!
!! SOURCE

subroutine size_dvxc(ixc,ndvxc,ngr2,nd2vxc,nspden,nvxcdgr,order,&
& add_tfw,xc_funcs) ! Optional

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'size_dvxc'
!End of the abilint section

 implicit none

!Arguments----------------------
 integer, intent(in) :: ixc,nspden,order
 integer, intent(out) :: ndvxc,nd2vxc,ngr2,nvxcdgr
 logical, intent(in), optional :: add_tfw
 type(libxc_functional_type),intent(in),optional :: xc_funcs(2)

!Local variables----------------
 logical :: add_tfw_,isgga,ismgga,is_hybrid

! *************************************************************************

 add_tfw_=.false.;if (present(add_tfw)) add_tfw_=add_tfw
 isgga=.false. ; ismgga=.false. ; is_hybrid=.false.
 if(ixc<0)then
   if(present(xc_funcs))then
     isgga=libxc_functionals_isgga(xc_functionals=xc_funcs)
     ismgga=libxc_functionals_ismgga(xc_functionals=xc_funcs)
     is_hybrid=libxc_functionals_is_hybrid(xc_functionals=xc_funcs)
   else
     isgga=libxc_functionals_isgga()
     ismgga=libxc_functionals_ismgga()
     is_hybrid=libxc_functionals_is_hybrid()
   end if
 end if

 ngr2=0
 nvxcdgr=0
 ndvxc=0
 nd2vxc=0

!Dimension for the gradient of the density (only allocated for GGA or mGGA)
 if ((ixc>=11.and.ixc<=17).or.(ixc>=23.and.ixc<=24).or.ixc==26.or.ixc==27.or. &
& (ixc>=31.and.ixc<=34).or.(ixc==41.or.ixc==42).or.ixc==1402000.or.(add_tfw_)) ngr2=2*min(nspden,2)-1
 if (ixc<0.and.isgga.or.ismgga.or.is_hybrid) ngr2=2*min(nspden,2)-1

!A-Only Exc and Vxc
!=======================================================================================
 if (order**2 <= 1) then
   if (((ixc>=11 .and. ixc<=15) .and. ixc/=13) .or. (ixc>=23 .and. ixc<=24) .or. &
&   (ixc==41 .or. ixc==42) .or. ixc==1402000) nvxcdgr=3
   if (ixc==16.or.ixc==17.or.ixc==26.or.ixc==27) nvxcdgr=2
   if (ixc<0) nvxcdgr=3
   if (ixc>=31 .and. ixc<=34) nvxcdgr=3 !Native fake metaGGA functionals (for testing purpose only)
   if (add_tfw_) nvxcdgr=3

!  B- Exc+Vxc and other derivatives
!  =======================================================================================
 else

!  Definition of ndvxc and nvxcdgr, 2nd dimension of the arrays of 2nd-order derivatives
!  -------------------------------------------------------------------------------------
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
   else if ((ixc>=11 .and. ixc<=15 .and. ixc/=13) .or. ixc==23 .or. ixc==41 .or. ixc==42) then
!    Routine xcpbe, with different options (optpbe) and orders (order)
     ndvxc=15
     nvxcdgr=3
   else if(ixc==16 .or. ixc==17 .or. ixc==26 .or. ixc==27 ) then
     nvxcdgr=2
   else if (ixc==50) then
     ndvxc=1 !  IIT xc (no spin-pol)
   else if (ixc==1402000) then
     ndvxc=15
     nvxcdgr=3
   else if (ixc<0) then
     if(libxc_functionals_isgga() .or. libxc_functionals_ismgga() .or. libxc_functionals_is_hybrid()) then
       ndvxc=15
     else
       ndvxc=3
     end if
     nvxcdgr=3
   end if
   if (add_tfw_) nvxcdgr=3

!  Definition of nd2vxc, 2nd dimension of the array of 3rd-order derivatives
!  -------------------------------------------------------------------------------------
   if (order==3) then
     if (ixc==3) nd2vxc=1 ! Non spin polarized LDA case
     if ((ixc>=7 .and. ixc<=10) .or. (ixc==13)) nd2vxc=3*min(nspden,2)-2
!    Following line to be corrected when the calculation of d2vxcar is implemented for these functionals
     if ((ixc>=11 .and. ixc<=15 .and. ixc/=13) .or. (ixc==23.and.ixc<=24) .or. (ixc==41.or.ixc==42)) nd2vxc=1
     if(ixc==1402000)nd2vxc=3*min(nspden,2)-2
     if ((ixc<0.and.(.not.(libxc_functionals_isgga().or.libxc_functionals_ismgga().or.libxc_functionals_is_hybrid() )))) &
&     nd2vxc=3*min(nspden,2)-2
   end if

 end if

end subroutine size_dvxc
!!***

end module m_drivexc
!!***
