!{\src2tex{textfont=tt}}
!!****f* ABINIT/drivexc_main
!! NAME
!! drivexc_main
!!
!! FUNCTION
!! Driver of XC functionals.Optionally, deliver the XC kernel, or even the derivative
!! of the XC kernel (the third derivative of the XC energy)
!!
!! COPYRIGHT
!! Copyright (C) 2012-2017 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!  [exexch]=choice of <<<local>>> exact exchange. Active if exexch=3 (only for GGA)
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
!!  [xc_tb09_c]=c parameter for the TB09 functional
!!  [hyb_mixing]= mixing parameter for the native PBEx functionals (ixc=41 and 42)

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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine drivexc_main(exc,ixc,mgga,ndvxc,nd2vxc,ngr2,npts,nspden,nvxcgrho,order,rho,vxcrho,xclevel, &
&                       dvxc,d2vxc,el_temp,exexch,fxcT,grho2,hyb_mixing,lrho,tau,vxcgrho,vxclrho,vxctau,xc_tb09_c) ! Optional arguments

 use defs_basis
 use m_profiling_abi
 use m_errors
 use libxc_functionals

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'drivexc_main'
 use interfaces_41_xc_lowlevel, except_this_one => drivexc_main
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

!Local variables-------------------------------
!scalars
 real(dp) :: hyb_mixing_,xc_tb09_c_

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

 if (ixc<0) then
   if (mgga==1) then
     if (abs(xc_tb09_c_-99.99_dp)>tol12) then
       if (present(exexch)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho, &
&         lrho_updn=lrho,vxclrho=vxclrho,tau_updn=tau,vxctau=vxctau, &
&         exexch=exexch,xc_tb09_c=xc_tb09_c_)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho, &
&         lrho_updn=lrho,vxclrho=vxclrho,tau_updn=tau,vxctau=vxctau, &
&         xc_tb09_c=xc_tb09_c_)
       end if
     else
       if (present(exexch)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho, &
&         lrho_updn=lrho,vxclrho=vxclrho,tau_updn=tau,vxctau=vxctau, &
&         exexch=exexch)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho, &
&         lrho_updn=lrho,vxclrho=vxclrho,tau_updn=tau,vxctau=vxctau)
       end if
     end if
   else if (libxc_functionals_isgga()) then
     if (order**2<=1) then
       if (present(exexch)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho,exexch=exexch)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho)
       end if
     else
       if (present(exexch)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho,dvxc=dvxc,exexch=exexch)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         grho2_updn=grho2,vxcgrho=vxcgrho,dvxc=dvxc)
       end if
     end if
   else
     if (order**2<=1) then
       if (present(exexch)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
         exexch=exexch)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho)
       end if
     else if (order**2<=4) then
       if (present(exexch)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         dvxc=dvxc, exexch=exexch)
       else
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         dvxc=dvxc)
       end if
     else
       if (present(exexch)) then
         call drivexc(exc,ixc,npts,nspden,order,rho,vxcrho,ndvxc,ngr2,nd2vxc,nvxcgrho, &
&         dvxc=dvxc,d2vxc=d2vxc, exexch=exexch)
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
