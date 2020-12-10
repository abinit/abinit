!!****m* ABINIT/m_psp_hgh
!! NAME
!!  m_psp_hgh
!!
!! FUNCTION
!! Initialize pspcod=2, 3, 10 pseudopotentials (GTH)
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, MT, FD, PT)
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

module m_psp_hgh

 use defs_basis
 use m_splines
 use m_abicore
 use m_errors
 use m_dtset

 use defs_datatypes, only : pseudopotential_type
 use m_special_funcs,  only : abi_derfc
#if defined HAVE_BIGDFT
 use BigDFT_API, only: atomic_info
#endif
 use m_wvl_descr_psp,  only : wvl_descr_psp_fill

 implicit none

 private
!!***

 public :: psp2in
 public :: psp3in
 public :: psp10in
!!***

contains
!!***

!!****f* m_psp_hgh/psp2in
!! NAME
!! psp2in
!!
!! FUNCTION
!! Initialize pspcod=2 pseudopotentials (GTH format):
!! continue to read the file, then compute the corresponding
!! local and non-local potentials.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  ipsp=id in the array of the pseudo-potential.
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  zion=nominal valence of atom as specified in psp file
!!
!! OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy,
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for each (l,n)
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+Zv/r) dr]$ (hartree)
!!  ffspl(mqgrid,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  nproj(mpsang)=number of projection functions for each angular momentum
!!  vlspl(mqgrid_vl,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  dvlspl(mqgrid_vl,2)=dVloc(r)/dr and second derivatives from spline fit (only
!!                      allocated if vlspl_recipSpace is false.
!!
!! SIDE EFFECTS
!!  Input/output
!!  lmax : at input =value of lmax mentioned at the second line of the psp file
!!    at output= 1
!!  psps <type(pseudopotential_type)>=at output, values depending on the read
!!                                    pseudo are set.
!!   | lmnmax(IN)=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!   |           =if useylm=0, max number of (l,n)   comp. over all type of psps
!!   | lnmax(IN)=max. number of (l,n) components over all type of psps
!!   |           angular momentum of nonlocal pseudopotential
!!   | mpsang(IN)= 1+maximum angular momentum for nonlocal pseudopotentials
!!   | mqgrid_ff(IN)=dimension of q (or G) grid for nl form factors (array ffspl)
!!   | mqgrid_vl(IN)=dimension of q (or G) grid or r grid (if vlspl_recipSpace = .false.)
!!   | qgrid_ff(mqgrid_ff)(IN)=values of q on grid from 0 to qmax (bohr^-1) for nl form factors
!!   | qgrid_vl(mqgrid_vl)(IN)=values of q on grid from 0 to qmax (bohr^-1) for Vloc
!!   |                         if vlspl_recipSpace is .true. else values of r on grid from
!!   |                         0 to 2pi / qmax * mqgrid_ff (bohr).
!!   | useylm(IN)=governs the way the nonlocal operator is to be applied:
!!   |            1=using Ylm, 0=using Legendre polynomials
!!   | vlspl_recipSpace(IN)=.true. if pseudo are expressed in reciprocal space.
!!   | gth_params(OUT)=store GTH coefficients and parameters.
!!
!! PARENTS
!!      m_pspini
!!
!! CHILDREN
!!      spline,zhpev
!!
!! SOURCE

subroutine psp2in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps,vlspl,dvlspl,zion)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ipsp,lmax
 real(dp),intent(in) :: zion
 real(dp),intent(out) :: epsatm
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,intent(out) :: indlmn(6,psps%lmnmax),nproj(psps%mpsang)
 real(dp),intent(out) :: dvlspl(psps%mqgrid_vl,2),ekb(psps%lnmax)
 real(dp),intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax) !vz_i
 real(dp),intent(out) :: vlspl(psps%mqgrid_vl,2)

!Local variables-------------------------------
!scalars
 integer :: iln,index,ipsang,kk,ll,mm
 real(dp) :: cc1,cc2,cc3,cc4,h1p,h1s,h2s,rloc,rrp,rrs
 real(dp) :: yp1,ypn
 character(len=500) :: message,errmsg
!arrays
 real(dp),allocatable :: work_space(:),work_spl(:)
 real(dp),allocatable :: dvloc(:)

! ***************************************************************************

!Set various terms to 0 in case not defined below
!GTH values
 rloc=0.d0
 cc1=0.d0
 cc2=0.d0
 cc3=0.d0
 cc4=0.d0
 rrs=0.d0
 h1s=0.d0
 h2s=0.d0
 rrp=0.d0
 h1p=0.d0
 nproj(1:psps%mpsang)=0

!Read and write different lines of the pseudopotential file
 read (tmp_unit,*, err=10, iomsg=errmsg) rloc,cc1,cc2,cc3,cc4
 write(message, '(a,f12.7)' ) ' rloc=',rloc
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')
 write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )'  cc1=',cc1,'; cc2=',cc2,'; cc3=',cc3,'; cc4=',cc4
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 read (tmp_unit,*, err=10, iomsg=errmsg) rrs,h1s,h2s
 write(message, '(a,f12.7,a,f12.7,a,f12.7)' )'  rrs=',rrs,'; h1s=',h1s,'; h2s=',h2s
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 read (tmp_unit,*, err=10, iomsg=errmsg) rrp,h1p
 write(message, '(a,f12.7,a,f12.7)' )'  rrp=',rrp,'; h1p=',h1p
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!Store the coefficients.
 psps%gth_params%set(ipsp)          = .true.
 psps%gth_params%psppar(0, :, ipsp) = (/ rloc, cc1, cc2, cc3, cc4, 0.d0, 0.d0 /)
 psps%gth_params%psppar(1, :, ipsp) = (/ rrs,  h1s, h2s, 0.d0, 0.d0, 0.d0, 0.d0 /)
 psps%gth_params%psppar(2, :, ipsp) = (/ rrp,  h1p, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
 if (dtset%usewvl == 1) then
   call wvl_descr_psp_fill(psps%gth_params, ipsp, 0, int(psps%zionpsp(ipsp)), int(psps%znuclpsp(ipsp)), tmp_unit)
 end if

 if (abs(h1s)>1.d-08) nproj(1)=1
 if (abs(h2s)>1.d-08) nproj(1)=2

 if (abs(h1p)>1.d-08) then
   if(psps%mpsang<2)then
     write(message, '(a,es12.4,a,a,a,i2,a)' )&
&     'With non-zero h1p (=',h1p,'), mpsang should be at least 2,',ch10,&
&     'while mpsang=',psps%mpsang,'.'
     ABI_ERROR(message)
   end if
   nproj(2)=1
   if (lmax<1) then
     write(message, '(a,i5,a,e12.4,a,a,a,a)' )&
&     'Input lmax=',lmax,' disagree with input h1p=',h1p,'.',&
&     'Your pseudopotential is incoherent.',ch10,&
&     'Action: correct your pseudopotential file.'
     ABI_ERROR(message)
   end if
 end if

!Initialize array indlmn array giving l,m,n,lm,ln,s for i=lmn
 index=0;iln=0;indlmn(:,:)=0
 do ipsang=1,lmax+1
   if(nproj(ipsang)>0)then
     ll=ipsang-1
     do kk=1,nproj(ipsang)
       iln=iln+1
       do mm=1,2*ll*psps%useylm+1
         index=index+1
         indlmn(1,index)=ll
         indlmn(2,index)=mm-ll*psps%useylm-1
         indlmn(3,index)=kk
         indlmn(4,index)=ll*ll+(1-psps%useylm)*ll+mm
         indlmn(5,index)=iln
         indlmn(6,index)=1
       end do
     end do
   end if
 end do

!First, the local potential --
!compute q^2V(q) or V(r)
!MJV NOTE: psp2lo should never be called with dvspl unallocated, which
!is possible unless .not.psps%vlspl_recipSpace
 ABI_ALLOCATE(dvloc,(psps%mqgrid_vl))
 call psp2lo(cc1,cc2,cc3,cc4,dvloc,epsatm,psps%mqgrid_vl,psps%qgrid_vl,&
& vlspl(:,1),rloc,psps%vlspl_recipSpace,yp1,ypn,zion)

!Fit spline to (q^2)V(q) or V(r)
 ABI_ALLOCATE(work_space,(psps%mqgrid_vl))
 ABI_ALLOCATE(work_spl,(psps%mqgrid_vl))
 call spline (psps%qgrid_vl,vlspl(:,1),psps%mqgrid_vl,yp1,ypn,work_spl)
 vlspl(:,2)=work_spl(:)
 if (.not.psps%vlspl_recipSpace) then
   dvlspl(:,1) = dvloc
   call spline (psps%qgrid_vl,dvlspl(:,1),psps%mqgrid_vl,yp1,ypn,work_spl)
   dvlspl(:,2)=work_spl(:)
 end if

 ABI_DEALLOCATE(work_space)
 ABI_DEALLOCATE(work_spl)
 ABI_DEALLOCATE(dvloc)


!Second, compute KB energies and form factors and fit splines
 ekb(:)=0.0d0
!First check if any nonlocal projectors are being used
 if (maxval(nproj(1:lmax+1))>0) then
   call psp2nl(ekb,ffspl,h1p,h1s,h2s,psps%lnmax,psps%mqgrid_ff,psps%qgrid_ff,rrp,rrs)
 end if

 return

 ! Handle IO error
 10 continue
 ABI_ERROR(errmsg)

end subroutine psp2in
!!***

!!****f* m_psp_hgh/psp2nl
!! NAME
!! psp2nl
!!
!! FUNCTION
!! Goedecker-Teter-Hutter nonlocal pseudopotential (from preprint of 1996).
!! Uses Gaussians for fully nonlocal form, analytic expressions.
!!
!! INPUTS
!!  h1p=factor defining strength of 1st projector for l=1 channel
!!  h1s=factor defining strength of 1st projector for l=0 channel
!!  h2s=factor defining strength of 2nd projector for l=0 channel
!!  lnmax=max. number of (l,n) components over all type of psps
!!  mqgrid=number of grid points for qgrid
!!  qgrid(mqgrid)=array of |G| values
!!  rrp=core radius for p channel (bohr)
!!  rrs=core radius for s channel (bohr)
!!
!! OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy
!!  ffspl(mqgrid,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum
!!   and each projector
!!
!! PARENTS
!!      m_psp_hgh
!!
!! CHILDREN
!!      spline,zhpev
!!
!! SOURCE

subroutine psp2nl(ekb,ffspl,h1p,h1s,h2s,lnmax,mqgrid,qgrid,rrp,rrs)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lnmax,mqgrid
 real(dp),intent(in) :: h1p,h1s,h2s,rrp,rrs
!arrays
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(inout) :: ekb(lnmax),ffspl(mqgrid,2,lnmax) !vz_i

!Local variables-------------------------------
!scalars
 integer :: iln,iqgrid
 real(dp) :: qmax,yp1,ypn
!arrays
 real(dp),allocatable :: work(:)

! *************************************************************************

 ABI_ALLOCATE(work,(mqgrid))

!Kleinman-Bylander energies ekb were set to zero in calling program

!Compute KB energies
 iln=0
 if (abs(h1s)>1.d-12) then
   iln=iln+1
   ekb(iln)=h1s*32.d0*rrs**3*(pi**(2.5d0)/(4.d0*pi)**2)
 end if
 if (abs(h2s)>1.d-12) then
   iln=iln+1
   ekb(iln) =h2s*(128.d0/15.d0)*rrs**3*(pi**(2.5d0)/(4.d0*pi)**2)
 end if
 if (abs(h1p)>1.d-12) then
   iln=iln+1
   ekb(iln)=h1p*(64.d0/3.d0)*rrp**5*(pi**(2.5d0)/(4.d0*pi)**2)
 end if

!Compute KB form factor
 iln=0

!l=0 first projector
 if (abs(h1s)>1.d-12) then
   iln=iln+1
   do iqgrid=1,mqgrid
     ffspl(iqgrid,1,iln)=exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrs)**2)
   end do
!  Compute yp1,ypn=derivatives of f(q) at q=0, q=qgrid(mqgrid)
   yp1=0.d0
   qmax=qgrid(mqgrid)
   ypn=-4.d0*pi**2*qmax*rrs**2*exp(-0.5d0*(two_pi*qmax*rrs)**2)
!  Fit spline to get second derivatives by spline fit
   call spline(qgrid,ffspl(:,1,iln),mqgrid,yp1,ypn,ffspl(:,2,iln))
!  else
!  or else put first projector nonlocal correction at l=0 to 0
!  ffspl(:,:,iln)=0.0d0
 end if

!l=0 second projector
 if (abs(h2s)>1.d-12) then
   iln=iln+1
   do iqgrid=1,mqgrid
     ffspl(iqgrid,1,iln)=exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrs)**2) * &
&     (3.d0-(two_pi*qgrid(iqgrid)*rrs)**2)
   end do
!  Compute yp1,ypn=derivatives of f(q) at q=0, q=qgrid(mqgrid)
   yp1=0.d0
   qmax=qgrid(mqgrid)
   ypn=4.d0*pi**2*qmax*rrs**2*exp(-0.5d0*(two_pi*qmax*rrs)**2) * &
&   (-5.d0+(two_pi*qmax*rrs)**2)
!  Fit spline to get second derivatives by spline fit
   call spline(qgrid,ffspl(:,1,iln),mqgrid,yp1,ypn,ffspl(:,2,iln))
!  else if(mproj>=2)then
!  or else put second projector nonlocal correction at l=0 to 0
!  ffspl(:,:,iln)=0.0d0
 end if

!l=1 first projector
 if (abs(h1p)>1.d-12) then
   iln=iln+1
   do iqgrid=1,mqgrid
     ffspl(iqgrid,1,iln)=exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrp)**2) * &
&     (two_pi*qgrid(iqgrid))
   end do
!  Compute yp1,ypn=derivatives of f(q) at q=0, q=qgrid(mqgrid)
   yp1=two_pi
   qmax=qgrid(mqgrid)
   ypn=-two_pi*((two_pi*qmax*rrp)**2-1.d0) * exp(-0.5d0*(two_pi*qmax*rrp)**2)
!  Fit spline to get second derivatives by spline fit
   call spline(qgrid,ffspl(:,1,iln),mqgrid,yp1,ypn,ffspl(:,2,iln))
!  else if(mpsang>=2)then
!  or else put first projector l=1 nonlocal correction to 0
!  ffspl(:,:,iln)=0.0d0
 end if

 ABI_DEALLOCATE(work)

end subroutine psp2nl
!!***

!!****f* ABINIT/psp2lo
!! NAME
!! psp2lo
!!
!! FUNCTION
!! Treat local part of Goedecker-Teter-Hutter pseudopotentials (pspcod=2),
!! as well as Hartwigsen-Goedecker-Hutter pseudopotentials (pspcod=3)
!!
!! INPUTS
!!  cc1,2,3,4=parameters from analytic pseudopotential form
!!  mqgrid=number of grid points in q from 0 to qmax.
!!  qgrid(mqgrid)=values of q (or G) on grid from 0 to qmax (bohr^-1)
!!                if vlspl_recipSpace is .true. else values of r on grid from
!!                0 to 2pi / qmax * mqgrid_ff (bohr).
!!  rloc=local pseudopotential core radius (bohr)
!!  vlspl_recipSpace= .true. if computation of vlspl is done in reciprocal space
!!  zion=valence charge of atom
!!  parameters for local potential: rloc,c1,c2,c3,c4
!!
!! OUTPUT
!!  dvloc(mqgrid)=dVloc(r)/dr (only allocated if vlspl_recipSpace is false).
!!  epsatm=$4\pi\int[r^2 (v(r)+\frac{Zv}{r} dr]$
!!{{\ \begin{eqnarray}
!!  q2vq(mqgrid)&=&q^2 v(q) \nonumber \\
!!  &=&-Zv/\pi
!!   +q^2 4\pi\int[(\frac{\sin(2\pi qr)}{2\pi qr})(r^2 v(r)+r Zv)dr]\nonumber\\
!!  &=&\exp(-K^2*rloc^2/2) \nonumber \\
!!  &&   *(-\frac{zion}{\pi}+(\frac{K^2*rloc^3}{\sqrt{2*\pi}}*
!!       (c1+c2*(3-(rloc*K)^2) \nonumber \\
!!  &&    +c3*(15-10(rloc*K)^2+(rloc*K)^4) \nonumber \\
!!  &&    +c4*(105-105*(rloc*K)^2+21*(rloc*K)^4-(rloc*K)^6)) \nonumber
!!\end{eqnarray} }}
!! for GTH vloc with $K=(2\pi q)$.
!!  yp1,ypn=derivative of q^2 v(q) wrt q at q=0 and q=qmax
!!   (needed for spline fitter).
!!
!! PARENTS
!!      m_psp_hgh
!!
!! CHILDREN
!!      spline,zhpev
!!
!! SOURCE

subroutine psp2lo(cc1,cc2,cc3,cc4,dvloc,epsatm,mqgrid,qgrid,q2vq,&
&  rloc,vlspl_recipSpace,yp1,ypn,zion)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mqgrid
 real(dp),intent(in) :: cc1,cc2,cc3,cc4,rloc,zion
 real(dp),intent(out) :: epsatm,yp1,ypn
 logical,intent(in) :: vlspl_recipSpace
!arrays
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(out) :: dvloc(mqgrid),q2vq(mqgrid)

!Local variables-------------------------------
!scalars
 integer :: iqgrid
 real(dp) :: erfValue,gaussValue,polyValue,qmax,rq,rq2
 character(len=500) :: message

! *************************************************************************

!Compute epsatm = lim(q->0) [Vloc(q) + zion/(Pi*q^2)]
 epsatm=2.d0*pi*rloc**2*zion+(2.d0*pi)**(1.5d0)*rloc**3*&
& (cc1+3.d0*cc2+15.d0*cc3+105.d0*cc4)

!If vlspl_recipSpace is .true., we compute V(q)*q^2 in reciprocal space,
!else we compute V(r) in real space.
 if (vlspl_recipSpace) then
   write(message, '(a)' ) '-  Local part computed in reciprocal space.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

!  d(q^2*V(q))/d(q) at q=0 and q=qmax
   qmax=qgrid(mqgrid)
   rq2=(2.d0*pi*qmax*rloc)**2
   yp1=0.d0
   ypn= (2.d0*pi*qmax*rloc**2)*exp(-0.5d0*rq2)* &
&   (2.d0*zion + sqrt(2.d0*pi)*rloc*&
&   (cc1*(2.d0-rq2) + cc2*(6.d0-7.d0*rq2+rq2**2) +&
&   cc3*(30.d0-55.d0*rq2+16.d0*rq2**2-rq2**3) +&
&   cc4*(210.d0-525.d0*rq2+231.d0*rq2**2-29.d0*rq2**3+rq2**4)))
!  ypn has been tested against Maple-derived expression.

!  Compute q^2*vloc(q) on uniform grid
   do iqgrid=1,mqgrid
     rq2=(2.d0*pi*qgrid(iqgrid)*rloc)**2
     q2vq(iqgrid)=exp(-0.5d0*rq2)*(-zion/pi+rq2*(rloc/sqrt(2.d0*pi)) *&
&     ( cc1 + cc2*(3.d0-rq2) + cc3*(15.d0-10.d0*rq2+rq2**2) +&
&     cc4*(105.d0-rq2*(105.d0-rq2*(21.d0-rq2)))  ))
   end do
 else
   write(message, '(a)' ) '-  Local part computed in real space.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

!  Compute derivatives for splines computations
   yp1 = 0.d0
   rq2 = (qgrid(mqgrid) / rloc) ** 2
   erfValue = abi_derfc(sqrt(0.5d0 * rq2))
   ypn = - 2.0d0 * zion / sqrt(2.d0 * pi) / qgrid(mqgrid) / rloc
   ypn = ypn - rq2 * (cc1 + cc2 * rq2 + cc3 * rq2 ** 2 + cc4 * rq2 ** 3) / qgrid(mqgrid)
   ypn = ypn + (2.d0 * cc2 * rq2 + 4.d0 * cc3 * rq2 ** 2 + 6.d0 * cc4 * rq2 ** 3) / qgrid(mqgrid)
   ypn = ypn * exp(-0.5d0 * rq2)
   ypn = ypn + zion / qgrid(mqgrid) ** 2 * erfValue
!  Note that ypn has been calculated on a full-proof a4 paper sheet.

!  Compute local potential and its first derivatives.
   do iqgrid = 1, mqgrid, 1
     rq2 = (qgrid(iqgrid) / rloc) ** 2
!    Compute erf() part
!    Case r = 0
     gaussValue = exp(-0.5d0 * rq2)
     if (qgrid(iqgrid) == 0.d0) then
       q2vq(iqgrid) = -zion / rloc * sqrt(2.d0 / pi)
       dvloc(iqgrid) = 0.d0
     else
       erfValue = abi_derfc(sqrt(0.5d0 * rq2))
       q2vq(iqgrid) = -zion / qgrid(iqgrid) * (1.0d0 - erfValue)
       dvloc(iqgrid) = - sqrt(2.d0 / pi) * zion * gaussValue / (qgrid(iqgrid) * rloc) - &
&       q2vq(iqgrid) / qgrid(iqgrid)
     end if
!    Add the gaussian part
     polyValue = cc1 + cc2 * rq2 + cc3 * rq2 ** 2 + cc4 * rq2 ** 3
     q2vq(iqgrid) = q2vq(iqgrid) + gaussValue * polyValue
     rq = qgrid(iqgrid) / rloc
     dvloc(iqgrid) = dvloc(iqgrid) - qgrid(iqgrid) / rloc ** 2 * gaussValue * polyValue + &
&     gaussValue * (2.0d0 * cc2 * rq / rloc + 3.0d0 * cc3 * rq ** 3 / rloc + &
&     6.0d0 * cc4 * rq ** 5 / rloc)
   end do

   write(message, '(a,f12.7,a,a,f12.7,a,a,a,f12.7)' ) &
&   '  | dr spline step is : ', qgrid(2), ch10, &
&   '  | r > ', qgrid(mqgrid) ,' is set to 0.', ch10, &
&   '  | last non-nul potential value is : ', q2vq(mqgrid)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

end subroutine psp2lo
!!***

!!****f* ABINIT/psp3in
!! NAME
!! psp3in
!!
!! FUNCTION
!! Initialize pspcod=3 pseudopotentials (HGH psps PRB58,3641(1998) [[cite:Hartwigsen1998]]):
!! continue to read the file, then compute the corresponding
!! local and non-local potentials.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  pspso=spin-orbit characteristics, govern the content of ffspl and ekb
!!   if =0 : this input requires NO spin-orbit characteristics of the psp
!!   if =2 : this input requires HGH characteristics of the psp
!!   if =3 : this input requires HFN characteristics of the psp
!!  ipsp=id in the array of the pseudo-potential.
!!  zion=nominal valence of atom as specified in psp file
!!
!! OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy,
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for each (l,n)
!!             if any, spin-orbit components begin at l=mpsang+1
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$ (hartree)
!!  ffspl(mqgrid_ff,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector; if any, spin-orbit components begin at l=mpsang+1
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  nproj(mpssoang)=number of projection functions for each angular momentum
!!  vlspl(mqgrid_ff,2)=q^2 Vloc(q) and second derivatives from spline fit
!!
!! SIDE EFFECTS
!!  Input/output
!!  lmax : at input =value of lmax mentioned at the second line of the psp file
!!    at output= 1
!!  psps <type(pseudopotential_type)>=at output, values depending on the read
!!                                    pseudo are set.
!!   | lmnmax(IN)=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!   |           =if useylm=0, max number of (l,n)   comp. over all type of psps
!!   | lnmax(IN)=max. number of (l,n) components over all type of psps
!!   |           angular momentum of nonlocal pseudopotential
!!   | mpsang(IN)= 1+maximum angular momentum for nonlocal pseudopotentials
!!   | mpssoang(IN)= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!   | mqgrid_ff(IN)=dimension of q (or G) grid for arrays.
!!   | qgrid_ff(mqgrid_ff)(IN)=values of q on grid from 0 to qmax (bohr^-1) for nl form factors
!!   | useylm(IN)=governs the way the nonlocal operator is to be applied:
!!   |            1=using Ylm, 0=using Legendre polynomials
!!
!! PARENTS
!!      m_pspini
!!
!! CHILDREN
!!      spline,zhpev
!!
!! SOURCE

subroutine psp3in(dtset, ekb, epsatm, ffspl, indlmn, ipsp, lmax, nproj, psps, pspso, vlspl, zion)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ipsp,pspso
 integer,intent(inout) :: lmax
 real(dp),intent(in) :: zion
 real(dp),intent(out) :: epsatm
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,intent(out) :: indlmn(6,psps%lmnmax),nproj(psps%mpssoang)
 real(dp),intent(inout) :: ekb(psps%lnmax),ffspl(psps%mqgrid_ff,2,psps%lnmax)!vz_i
 real(dp),intent(out) :: vlspl(psps%mqgrid_ff,2)

!Local variables-------------------------------
!scalars
 integer :: iln,iln0,index,ipsang,jj,kk,ll,mm,mproj,nn,nso
 real(dp) :: cc1,cc2,cc3,cc4,h11d,h11f,h11p,h11s,h22d,h22p,h22s,h33d,h33p,h33s
 real(dp) :: k11d,k11f,k11p,k22d,k22p,k33d,k33p,rloc
 real(dp) :: rrd,rrf,rrp,rrs,yp1,ypn
 character(len=500) :: message,errmsg
!arrays
 real(dp),allocatable :: dvlspl(:),ekb_so(:,:),ekb_sr(:,:),ffspl_so(:,:,:,:)
 real(dp),allocatable :: ffspl_sr(:,:,:,:),work_space(:),work_spl(:)

! ***************************************************************************

!Set various terms to 0 in case not defined below
!HGH values
 rloc=zero ; rrs=zero  ; h11p=zero ; k33p=zero ; k11d=zero;
 cc1=zero  ; h11s=zero ; h22p=zero ; rrd=zero  ; k22d=zero;
 cc2=zero  ; h22s=zero ; h33p=zero ; h11d=zero ; k33d=zero;
 cc3=zero  ; h33s=zero ; k11p=zero ; h22d=zero ; h11f=zero;
 cc4=zero  ; rrp=zero  ; k22p=zero ; h33d=zero ; k11f=zero;
 rrf=zero
 nproj(1:psps%mpssoang)=0

!Read and write different lines of the pseudopotential file

 read (tmp_unit,*,err=10,iomsg=errmsg) rloc,cc1,cc2,cc3,cc4
 write(message, '(a,f12.7)' ) ' rloc=',rloc
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')
 write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )' cc1 =',cc1,'; cc2 =',cc2,'; cc3 =',cc3,'; cc4 =',cc4
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!For the time being, the s state line must be present and is read,
!even for local pseudopotentials (zero must appear)
 read (tmp_unit,*,err=10,iomsg=errmsg) rrs,h11s,h22s,h33s
 write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )' rrs =',rrs,'; h11s=',h11s,'; h22s=',h22s,'; h33s=',h33s
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 if (lmax > 0) then

   read (tmp_unit,*,err=10,iomsg=errmsg) rrp,h11p,h22p,h33p
   write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )' rrp =',rrp,'; h11p=',h11p,'; h22p=',h22p,'; h33p=',h33p
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   read (tmp_unit,*,err=10,iomsg=errmsg) k11p,k22p,k33p
   write(message, '(a,f12.7,a,f12.7,a,f12.7)' )'                    k11p=',k11p,'; k22p=',k22p,'; k33p=',k33p
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

 end if

 if (lmax > 1) then
   read (tmp_unit,*,err=10,iomsg=errmsg) rrd,h11d,h22d,h33d
   write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )' rrd =',rrd,'; h11d=',h11d,'; h22d=',h22d,'; h33d=',h33d
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   read (tmp_unit,*,err=10,iomsg=errmsg) k11d,k22d,k33d
   write(message, '(a,f12.7,a,f12.7,a,f12.7)' )'                    k11d=',k11d,'; k22d=',k22d,'; k33d=',k33d
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

 if (lmax > 2) then
   read (tmp_unit,*,err=10,iomsg=errmsg) rrf,h11f
   write(message, '(a,f12.7,a,f12.7)' )' rrf =',rrf,'; h11f=',h11f
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   read (tmp_unit,*,err=10,iomsg=errmsg) k11f
   write(message, '(a,f12.7)' )'                    k11f=',k11f
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if

 if (abs(h11s)>1.d-08) nproj(1)=1
 if (abs(h22s)>1.d-08) nproj(1)=2
 if (abs(h33s)>1.d-08) nproj(1)=3

 if (abs(h11p)>1.d-08) then
   nproj(2)=1
   if (lmax<1) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,'  does not agree with input h11p=',h11p,ch10,&
&     '  setting lmax to 1'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=1
   end if
 end if

 if (abs(h22p)>1.d-08) then
   nproj(2)=2
   if (lmax<1) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,' does not agree with input h22p=',h22p,ch10,&
&     '  setting lmax to 1'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=1
   end if
 end if

 if (abs(h33p)>1.d-08) then
   nproj(2)=3
   if (lmax<1) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,' does not agree with input h33p=',h33p,ch10,&
&     '  setting lmax to 1'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=1
   end if
 end if

 if (abs(h11d)>1.d-08) then
   nproj(3)=1
   if (lmax<2) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,'  does not agree with input h11d=',h11d,ch10,&
&     '  setting lmax to 2'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=2
   end if
 end if

 if (abs(h22d)>1.d-08) then
   nproj(3)=2
   if (lmax<2) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,'  does not agree with input h22d=',h22d,ch10,&
&     '  setting lmax to 2'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=2
   end if
 end if

 if (abs(h33d)>1.d-08) then
   nproj(3)=3
   if (lmax<2) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,' does not agree with input h33d=',h33d,ch10,&
&     '  setting lmax to 2'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=2
   end if
 end if

 if (abs(h11f)>1.d-08) then
   nproj(4)=1
   if (lmax<3) then
     write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&     ' psp3in : COMMENT -',ch10,&
&     '  input lmax=',lmax,' does not agree with input h11f=',h11f,ch10,&
&     '  setting lmax to 3'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
     lmax=3
   end if
 end if

 if(pspso/=0) then

   if (abs(k11p)>1.d-08) then
     nproj(psps%mpsang+1)=1
     if (lmax<1) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,'  does not agree with input k11p=',k11p,ch10,&
&       '  setting lmax to 1'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=1
     end if
   end if

   if (abs(k22p)>1.d-08) then
     nproj(psps%mpsang+1)=2
     if (lmax<1) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,' does not agree with input k22p=',k22p,ch10,&
&       '  setting lmax to 1'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=1
     end if
   end if


   if (abs(k33p)>1.d-08) then
     nproj(psps%mpsang+1)=3
     if (lmax<1) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,' does not agree with input k33p=',k33p,ch10,&
&       '  setting lmax to 1'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=1
     end if
   end if

   if (abs(k11d)>1.d-08) then
     nproj(psps%mpsang+2)=1
     if (lmax<2) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,'  does not agree with input k11d=',k11d,ch10,&
&       '  setting lmax to 2'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=2
     end if
   end if

   if (abs(k22d)>1.d-08) then
     nproj(psps%mpsang+2)=2
     if (lmax<2) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,'  does not agree with input k22d=',k22d,ch10,&
&       '  setting lmax to 2'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=2
     end if
   end if

   if (abs(k33d)>1.d-08) then
     nproj(psps%mpsang+2)=3
     if (lmax<2) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,' does not agree with input k33d=',k33d,ch10,&
&       '  setting lmax to 2'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=2
     end if
   end if

   if (abs(k11f)>1.d-08) then
     nproj(psps%mpsang+3)=1
     if (lmax<3) then
       write(message, '(a,a,a,a,i5,a,e12.4,a,a)' ) ch10,&
&       ' psp3in : COMMENT -',ch10,&
&       '  input lmax=',lmax,' does not agree with input k11f=',k11f,ch10,&
&       '  setting lmax to 3'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
       lmax=3
     end if
   end if

 end if

!Store the coefficients.
 psps%gth_params%set(ipsp)          = .true.
 psps%gth_params%psppar(0, :, ipsp) = (/ rloc, cc1, cc2, cc3, cc4, zero, zero /)
 psps%gth_params%psppar(1, :, ipsp) = (/ rrs,  h11s, h22s, h33s, zero, zero, zero /)
 psps%gth_params%psppar(2, :, ipsp) = (/ rrp,  h11p, h22p, h33p, zero, zero, zero /)
 psps%gth_params%psppar(3, :, ipsp) = (/ rrd,  h11d, h22d, h33d, zero, zero, zero /)
 psps%gth_params%psppar(4, :, ipsp) = (/ rrf,  h11f, zero, zero, zero, zero, zero /)

!Store the k coefficients
 psps%gth_params%psp_k_par(1, :, ipsp) = (/ zero, zero, zero /)
 psps%gth_params%psp_k_par(2, :, ipsp) = (/ k11p, k22p, k33p /)
 psps%gth_params%psp_k_par(3, :, ipsp) = (/ k11d, k22d, k33d /)
 psps%gth_params%psp_k_par(4, :, ipsp) = (/ k11f, zero, zero /)

!Additionnal wavelet parameters
 if (dtset%usewvl == 1) then
   call wvl_descr_psp_fill(psps%gth_params, ipsp, 0, int(psps%zionpsp(ipsp)), int(psps%znuclpsp(ipsp)), tmp_unit)
 end if

!Initialize array indlmn array giving l,m,n,ln,lm,s for i=lmn
 nso=1;if(pspso/=0) nso=2
 index=0;iln=0;indlmn(:,:)=0
 do nn=1,nso
   do ipsang=1+(nn-1)*(lmax+1),nn*lmax+1
     if (nproj(ipsang)>0) then
       ll=ipsang-(nn-1)*lmax-1
       do kk=1,nproj(ipsang)
         iln=iln+1
         do mm=1,2*ll*psps%useylm+1
           index=index+1
           indlmn(1,index)=ll
           indlmn(2,index)=mm-ll*psps%useylm-1
           indlmn(3,index)=kk
           indlmn(4,index)=ll*ll+(1-psps%useylm)*ll+mm
           indlmn(5,index)=iln
           indlmn(6,index)=nn
         end do
       end do
     end if
   end do
 end do

 ABI_ALLOCATE(dvlspl,(psps%mqgrid_ff))
!First, the local potential --  compute on q grid and fit spline
 call psp2lo(cc1,cc2,cc3,cc4,dvlspl,epsatm,psps%mqgrid_ff,psps%qgrid_ff,vlspl(:,1),rloc,.true.,yp1,ypn,zion)
 ABI_DEALLOCATE(dvlspl)

!DEBUG
!write(std_out,*)' psp3in : after psp2lo '
!ENDDEBUG

!Fit spline to q^2 V(q) (Numerical Recipes subroutine)
 ABI_ALLOCATE(work_space,(psps%mqgrid_ff))
 ABI_ALLOCATE(work_spl,(psps%mqgrid_ff))
 call spline (psps%qgrid_ff,vlspl(:,1),psps%mqgrid_ff,yp1,ypn,work_spl)
 vlspl(:,2)=work_spl(:)
 ABI_DEALLOCATE(work_space)
 ABI_DEALLOCATE(work_spl)

!Second, compute KB energies and form factors and fit splines
 ekb(:)=zero

!Check if any nonlocal projectors are being used
 mproj=maxval(nproj)
 if (mproj>0) then

   ABI_ALLOCATE(ekb_sr,(psps%mpsang,mproj))
   ABI_ALLOCATE(ffspl_sr,(psps%mqgrid_ff,2,psps%mpsang,mproj))
   ABI_ALLOCATE(ekb_so,(psps%mpsang,mproj))
   ABI_ALLOCATE(ffspl_so,(psps%mqgrid_ff,2,psps%mpsang,mproj))

   call psp3nl(ekb_sr,ffspl_sr,h11s,h22s,h33s,h11p,h22p,h33p,h11d,h22d,&
&   h33d,h11f,mproj,psps%mpsang,psps%mqgrid_ff,psps%qgrid_ff,rrd,rrf,rrp,rrs)
   if(pspso/=0) then
     call psp3nl(ekb_so,ffspl_so,zero,zero,zero,k11p,k22p,k33p,k11d,&
&     k22d,k33d,k11f,mproj,psps%mpsang,psps%mqgrid_ff,psps%qgrid_ff,rrd,rrf,rrp,rrs)
   end if


!  Convert ekb and ffspl
   iln0=0
   do jj=1,psps%lmnmax
     iln=indlmn(5,jj)
     if (iln>iln0) then
       iln0=iln
       if (indlmn(6,jj)<=1) then
         ekb(iln)=ekb_sr(1+indlmn(1,jj),indlmn(3,jj))
         ffspl(:,:,iln)=ffspl_sr(:,:,1+indlmn(1,jj),indlmn(3,jj))
       else
         ekb(iln)=ekb_so(1+indlmn(1,jj),indlmn(3,jj))
         ffspl(:,:,iln)=ffspl_so(:,:,1+indlmn(1,jj),indlmn(3,jj))
       end if
     end if
   end do

   ABI_DEALLOCATE(ekb_sr)
   ABI_DEALLOCATE(ffspl_sr)
   ABI_DEALLOCATE(ekb_so)
   ABI_DEALLOCATE(ffspl_so)
 end if

 return

 ! Handle IO error
 10 continue
 ABI_ERROR(errmsg)

end subroutine psp3in
!!***

!!****f* m_psp_hgh/psp3nl
!! NAME
!! psp3nl
!!
!! FUNCTION
!! Hartwigsen-Goedecker-Hutter nonlocal pseudopotential (from preprint of 1998).
!! Uses Gaussians for fully nonlocal form, analytic expressions.
!!
!! INPUTS
!!  h11s=factor defining strength of 1st projector for l=0 channel
!!  h22s=factor defining strength of 2nd projector for l=0 channel
!!  h33s=factor defining strength of 3rd projector for l=0 channel
!!  h11p=factor defining strength of 1st projector for l=1 channel
!!  h22p=factor defining strength of 2nd projector for l=1 channel
!!  h33p=factor defining strength of 2nd projector for l=1 channel
!!  h11d=factor defining strength of 1st projector for l=2 channel
!!  h22d=factor defining strength of 2nd projector for l=2 channel
!!  h33d=factor defining strength of 2nd projector for l=2 channel
!!  h11f=factor defining strength of 1st projector for l=3 channel
!!  mproj=maximum number of projectors in any channel
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=number of grid points for qgrid
!!  qgrid(mqgrid)=array of |G| values
!!  rrd=core radius for d channel (bohr)
!!  rrf=core radius for f channel (bohr)
!!  rrp=core radius for p channel (bohr)
!!  rrs=core radius for s channel (bohr)
!!
!! OUTPUT
!!  ekb(mpsang,mproj)=Kleinman-Bylander energies
!!  ffspl(mqgrid,2,mpssang,mproj)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projectors
!!
!! PARENTS
!!      m_psp_hgh
!!
!! CHILDREN
!!      spline,zhpev
!!
!! SOURCE

subroutine psp3nl(ekb,ffspl,h11s,h22s,h33s,h11p,h22p,h33p,h11d,h22d,&
&                  h33d,h11f,mproj,mpsang,mqgrid,qgrid,rrd,rrf,rrp,rrs)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mproj,mpsang,mqgrid
 real(dp),intent(in) :: h11d,h11f,h11p,h11s,h22d,h22p,h22s,h33d,h33p,h33s,rrd
 real(dp),intent(in) :: rrf,rrp,rrs
!arrays
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(out) :: ekb(mpsang,mproj),ffspl(mqgrid,2,mpsang,mproj)

!Local variables-------------------------------
!scalars
 integer :: info,iproj,iqgrid,ldz,mu,nproj,nu
 real(dp) :: qmax
 character(len=500) :: message
 character :: jobz,uplo
!arrays
 real(dp) :: ap(2,9),rwork1(9),work1(2,9),ww(3),yp1j(3),ypnj(3)
 real(dp),allocatable :: ppspl(:,:,:,:),uu(:,:),work(:),zz(:,:,:)

! *************************************************************************

 ABI_ALLOCATE(ppspl,(mqgrid,2,mpsang,mproj))
 ABI_ALLOCATE(work,(mqgrid))

 qmax=qgrid(mqgrid)
 jobz='v'
 uplo='u'
 ekb(:,:)=0.0d0

!---------------------------------------------------------------
!Treat s channel

 nproj=0
 ap(:,:)=0.0d0
!If there is at least one s-projector
 if  ( abs(h11s) >= 1.0d-8 ) then
   nproj=1 ; ldz=1 ; ap(1,1)=h11s
 end if
 nproj=1
!If there is a second projector
 if  ( abs(h22s) >= 1.0d-8 ) then
   nproj=2 ; ldz=2 ; ap(1,3)=h22s
   ap(1,2)=-0.5d0*sqrt(0.6d0)*h22s
 end if
!If there is a third projector
 if ( abs(h33s) >= 1.0d-8 ) then
   nproj=3 ; ldz=3 ; ap(1,6)=h33s
   ap(1,4)=0.5d0*sqrt(5.d0/21.d0)*h33s
   ap(1,5)=-0.5d0*sqrt(100.d0/63.d0)*h33s
 end if

 if(nproj/=0)then

   ABI_ALLOCATE(uu,(nproj,nproj))
   ABI_ALLOCATE(zz,(2,nproj,nproj))

   if (nproj > 1) then
     call ZHPEV(jobz,uplo,nproj,ap,ww,zz,ldz,work1,rwork1,info)
     uu(:,:)=zz(1,:,:)
   else
     ww(1)=h11s
     uu(1,1)=1.0d0
   end if

!  Initialization of ekb, and spline fitting
   do iproj=1,nproj
     ekb(1,iproj)=ww(iproj)*32.d0*(rrs**3)*(pi**2.5d0)/(4.d0*pi)**2
     if(iproj==1)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,1,1)=exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrs)**2)
       end do
       yp1j(1)=0.d0
       ypnj(1)=-(two_pi*rrs)**2*qmax*exp(-0.5d0*(two_pi*qmax*rrs)**2)
     else if(iproj==2)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,1,2)=2.0d0/sqrt(15.0d0)     &
&         *exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrs)**2) &
&         *( 3.d0-(two_pi*qgrid(iqgrid)*rrs)**2 )
       end do
       yp1j(2)=0.0d0
       ypnj(2)=2.0d0/sqrt(15.0d0)*(two_pi*rrs)**2*qmax &
&       *exp(-0.5d0*(two_pi*qmax*rrs)**2) * (-5.d0+(two_pi*qmax*rrs)**2)
     else if(iproj==3)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,1,3)=(4.0d0/3.0d0)/sqrt(105.0d0)*&
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrs)**2) * &
&         (15.0d0-10.0d0*(two_pi*qgrid(iqgrid)*rrs)**2 + &
&         (two_pi*qgrid(iqgrid)*rrs)**4)
       end do
       yp1j(3)=0.0d0
       ypnj(3)=(4.0d0/3.0d0)/sqrt(105.0d0)*exp(-0.5d0*(two_pi*qmax*rrs)**2) * &
&       (two_pi*rrs)**2*qmax*(-35.0d0+14d0*(two_pi*qmax*rrs)**2-(two_pi*qmax*rrs)**4)
     end if
     call spline(qgrid,ppspl(:,1,1,iproj),mqgrid,&
&     yp1j(iproj),ypnj(iproj),ppspl(:,2,1,iproj))
   end do

!  Linear combination using the eigenvectors
   ffspl(:,:,1,:)=0.0d0
   do mu=1,nproj
     do nu=1,nproj
       do iqgrid=1,mqgrid
         ffspl(iqgrid,1:2,1,mu)=ffspl(iqgrid,1:2,1,mu) &
&         +uu(nu,mu)*ppspl(iqgrid,1:2,1,nu)
       end do
     end do
   end do

   ABI_DEALLOCATE(uu)
   ABI_DEALLOCATE(zz)
 end if !  End condition on nproj(/=0)

!--------------------------------------------------------------------
!Now treat p channel

 nproj=0
 ap(:,:)=0.0d0
!If there is at least one projector
 if  ( abs(h11p) >= 1.0d-8 ) then
   nproj=1 ; ldz=1 ; ap(1,1)=h11p
 end if
!If there is a second projector
 if  ( abs(h22p) >= 1.0d-8 ) then
   nproj=2 ; ldz=2 ; ap(1,3)=h22p
   ap(1,2)=-0.5d0*sqrt(5.d0/7.d0)*h22p
 end if
!If there is a third projector
 if ( abs(h33p) >= 1.0d-8 ) then
   nproj=3 ; ldz=3 ; ap(1,6)=h33p
   ap(1,4)= (1.d0/6.d0)*sqrt(35.d0/11.d0)*h33p
   ap(1,5)=-(1.d0/6.d0)*(14.d0/sqrt(11.d0))*h33p
 end if

 if(nproj/=0)then

   ABI_ALLOCATE(uu,(nproj,nproj))
   ABI_ALLOCATE(zz,(2,nproj,nproj))

   if (nproj > 1) then
     call ZHPEV(jobz,uplo,nproj,ap,ww,zz,ldz,work1,rwork1,info)
     uu(:,:)=zz(1,:,:)
   else
     ww(1)=h11p
     uu(1,1)=1.0d0
   end if

!  Initialization of ekb, and spline fitting
   do iproj=1,nproj
     ekb(2,iproj)=ww(iproj)*64.d0*(rrp**5)*(pi**2.5d0)/(4.d0*pi)**2
     if(iproj==1)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,2,1)=(1.0d0/sqrt(3.0d0))* &
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrp)**2) * (two_pi*qgrid(iqgrid))
       end do
       yp1j(1)=two_pi*(1.0d0/sqrt(3.0d0))
       ypnj(1)=-two_pi*((two_pi*qmax*rrp)**2-1.d0)*exp(-0.5d0*(two_pi*qmax*rrp)**2)*&
&       (1.0d0/sqrt(3.0d0))
     else if(iproj==2)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,2,2)=(2.0d0/sqrt(105.0d0))* &
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrp)**2) * &
&         (two_pi*qgrid(iqgrid))*(5.0d0-(two_pi*qgrid(iqgrid)*rrp)**2)
       end do
       yp1j(2)=(5.0d0*two_pi)*(2.0d0/sqrt(105.0d0))
       ypnj(2)=(2.0d0/sqrt(105.0d0))*two_pi*exp(-0.5d0*(two_pi*qmax*rrp)**2)* &
&       (-8*(two_pi*qmax*rrp)**2 + (two_pi*qmax*rrp)**4 + 5.0d0)
     else if(iproj==3)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,2,3)=(4.0d0/3.0d0)/sqrt(1155d0)*&
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrp)**2) * &
&         (two_pi*qgrid(iqgrid))*&
&         (35.0d0-14.0d0*(two_pi*qgrid(iqgrid)*rrp)**2+(two_pi*qgrid(iqgrid)*rrp)**4)
       end do
       yp1j(3)=(35.0d0*two_pi)*(4.0d0/3.0d0)/sqrt(1155.0d0)
       ypnj(3)=(4.0d0/3.0d0)/sqrt(1155.0d0)*two_pi*exp(-0.5d0*(two_pi*qmax*rrp)**2)* &
&       (35.0d0-77.0d0*(two_pi*qmax*rrp)**2+19.0d0*(two_pi*qmax*rrp)**4 - &
&       (two_pi*qmax*rrp)**6)
     end if
     call spline(qgrid,ppspl(:,1,2,iproj),mqgrid,&
&     yp1j(iproj),ypnj(iproj),ppspl(:,2,2,iproj))
   end do

!  Linear combination using the eigenvectors
   ffspl(:,:,2,:)=0.0d0
   do mu=1,nproj
     do nu=1,nproj
       do iqgrid=1,mqgrid
         ffspl(iqgrid,1:2,2,mu)=ffspl(iqgrid,1:2,2,mu) &
&         +uu(nu,mu)*ppspl(iqgrid,1:2,2,nu)
       end do
     end do
   end do

   ABI_DEALLOCATE(uu)
   ABI_DEALLOCATE(zz)
 end if !  End condition on nproj(/=0)

!-----------------------------------------------------------------------
!Now treat d channel.

 nproj=0
 ap(:,:)=0.0d0
!If there is at least one projector
 if  ( abs(h11d) >= 1.0d-8 ) then
   nproj=1 ; ldz=1 ; ap(1,1)=h11d
 end if
!If there is a second projector
 if  ( abs(h22d) >= 1.0d-8 ) then
   nproj=2 ; ldz=2 ; ap(1,3)=h22d
   ap(1,2)=-0.5d0*sqrt(7.d0/9.d0)*h22d
 end if
!If there is a third projector. Warning : only two projectors are allowed.
 if ( abs(h33d) >= 1.0d-8 ) then
   write(message, '(a,a,a)' )&
&   '  only two d-projectors are allowed ',ch10,&
&   '  Action: check your pseudopotential file.'
   ABI_ERROR(message)
!  nproj=3 ; ldz=3 ; ap(1,6)=h33d
!  ap(1,4)= 0.5d0*sqrt(63.d0/143.d0)*h33d
!  ap(1,5)= -0.5d0*(18.d0/sqrt(143.d0))*h33d
 end if

 if(nproj/=0)then

   ABI_ALLOCATE(uu,(nproj,nproj))
   ABI_ALLOCATE(zz,(2,nproj,nproj))

   if (nproj > 1) then
     call ZHPEV(jobz,uplo,nproj,ap,ww,zz,ldz,work1,rwork1,info)
     uu(:,:)=zz(1,:,:)
   else
     ww(1)=h11d
     uu(1,1)=1.0d0
   end if

!  Initialization of ekb, and spline fitting
   do iproj=1,nproj
     ekb(3,iproj)=ww(iproj)*128.d0*(rrd**7)*(pi**2.5d0)/(4.d0*pi)**2
     if(iproj==1)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,3,1)=(1.0d0/sqrt(15.0d0))* &
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrd)**2) * (two_pi*qgrid(iqgrid))**2
       end do
       yp1j(1)=0.0d0
       ypnj(1)=(1.0d0/sqrt(15.0d0))*(two_pi**2)*&
&       exp(-0.5d0*(two_pi*qmax*rrd)**2)*qmax*(2d0-(two_pi*qmax*rrd)**2)
     else if(iproj==2)then
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,3,2)=(2.0d0/3.0d0)/sqrt(105.0d0)* &
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrd)**2) * &
&         ((two_pi*qgrid(iqgrid))**2)*(7.0d0-(two_pi*qgrid(iqgrid)*rrd)**2)
       end do
       yp1j(2)=0.0d0
       ypnj(2)=(2.0d0/3.0d0)/sqrt(105.0d0)*exp(-0.5d0*(two_pi*qmax*rrd)**2)* &
&       qmax*(two_pi**2)*( (two_pi*qmax*rrd)**4 - 11.0d0*(two_pi*qmax*rrd)**2 + 14.0d0)
     end if
     call spline(qgrid,ppspl(:,1,3,iproj),mqgrid,&
&     yp1j(iproj),ypnj(iproj),ppspl(:,2,3,iproj))
   end do

!  Linear combination using the eigenvectors
   ffspl(:,:,3,:)=0.0d0
   do mu=1,nproj
     do nu=1,nproj
       do iqgrid=1,mqgrid
         ffspl(iqgrid,1:2,3,mu)=ffspl(iqgrid,1:2,3,mu) &
&         +uu(nu,mu)*ppspl(iqgrid,1:2,3,nu)
       end do
     end do
   end do

   ABI_DEALLOCATE(uu)
   ABI_DEALLOCATE(zz)
 end if !  End condition on nproj(/=0)

!-----------------------------------------------------------------------
!Treat now f channel (max one projector ! - so do not use ppspl)

!l=3 first projector
 if (abs(h11f)>1.d-12) then
   ekb(4,1)=h11f*(256.0d0/105.0d0)*(rrf**9)*(pi**2.5d0)/(4.d0*pi)**2
   do iqgrid=1,mqgrid
     ffspl(iqgrid,1,4,1)=(two_pi*qgrid(iqgrid))**3* &
&     exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrf)**2)
   end do
!  Compute yp1,ypn=derivatives of f(q) at q=0, q=qgrid(mqgrid)
   yp1j(1)=0d0
   ypnj(1)=(two_pi**3)*qmax**2*exp(-0.5d0*(two_pi*qmax*rrf)**2)*&
&   (3.0d0-(two_pi*qmax*rrf)**2)
!  Fit spline to get second derivatives by spline fit
   call spline(qgrid,ffspl(:,1,4,1),mqgrid,&
&   yp1j(1),ypnj(1),ffspl(:,2,4,1))
 end if

!-----------------------------------------------------------------------

 ABI_DEALLOCATE(ppspl)
 ABI_DEALLOCATE(work)

end subroutine psp3nl
!!***

!!****f* m_psp_hgh/psp10in
!! NAME
!! psp10in
!!
!! FUNCTION
!! Initialize pspcod=10 pseudopotentials (formalism is the same as in HGH psps
!! PRB58,3641(1998) [[cite:Hartwigsen1998]], but the full h and k matrices are read, allowing for using
!! also subsequent developments such as Theor. Chem. Acc. 114, 145 (2005) [[cite:Dolg2005]]:
!! continue to read the file, then compute the corresponding
!! local and non-local potentials.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  pspso=spin-orbit characteristics, govern the content of ffspl and ekb
!!   if =0 : this input requires NO spin-orbit characteristics of the psp
!!   if =2 : this input requires HGH characteristics of the psp
!!   if =3 : this input requires HFN characteristics of the psp
!!  ipsp=id in the array of the pseudo-potential.
!!  zion=nominal valence of atom as specified in psp file
!!
!! OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy,
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for each (l,n)
!!             if any, spin-orbit components begin at l=mpsang+1
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$ (hartree)
!!  ffspl(mqgrid_ff,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector; if any, spin-orbit components begin at l=mpsang+1
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  nproj(mpssoang)=number of projection functions for each angular momentum
!!  vlspl(mqgrid_ff,2)=q^2 Vloc(q) and second derivatives from spline fit
!!
!! SIDE EFFECTS
!!  Input/output
!!  lmax : at input =value of lmax mentioned at the second line of the psp file
!!    at output= 1
!!  psps <type(pseudopotential_type)>=at output, values depending on the read
!!                                    pseudo are set.
!!   | lmnmax(IN)=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!   |           =if useylm=0, max number of (l,n)   comp. over all type of psps
!!   | lnmax(IN)=max. number of (l,n) components over all type of psps
!!   |           angular momentum of nonlocal pseudopotential
!!   | mpsang(IN)= 1+maximum angular momentum for nonlocal pseudopotentials
!!   | mpssoang(IN)= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!   | mqgrid_ff(IN)=dimension of q (or G) grid for arrays.
!!   | qgrid_ff(mqgrid_ff)(IN)=values of q on grid from 0 to qmax (bohr^-1) for nl form factors
!!   | useylm(IN)=governs the way the nonlocal operator is to be applied:
!!   |            1=using Ylm, 0=using Legendre polynomials
!!
!! PARENTS
!!      m_pspini
!!
!! CHILDREN
!!      spline,zhpev
!!
!! SOURCE

subroutine psp10in(dtset, ekb, epsatm, ffspl, indlmn, ipsp, lmax, nproj, psps, pspso, vlspl, zion)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ipsp,pspso
 integer,intent(inout) :: lmax
 real(dp),intent(in) :: zion
 real(dp),intent(out) :: epsatm
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,intent(out) :: indlmn(6,psps%lmnmax),nproj(psps%mpssoang)
 real(dp),intent(out) :: ekb(psps%lnmax) !vz_i
 real(dp),intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax) !vz_i
 real(dp),intent(out) :: vlspl(psps%mqgrid_ff,2)

!Local variables-------------------------------
!scalars
 integer :: ii,iln,iln0,index,ipsang,jj,kk,ll,mm,mproj,nn,nnonloc,nprl,nso
 real(dp) :: rloc,yp1,ypn
 character(len=500) :: message,errmsg
!arrays
 integer,allocatable :: dummy_nproj(:)
 real(dp) :: cc(4)
 real(dp),allocatable :: dvlspl(:),ekb_so(:,:),ekb_sr(:,:),ffspl_so(:,:,:,:)
 real(dp),allocatable :: ffspl_sr(:,:,:,:),hij(:,:,:),kij(:,:,:),rr(:)
 real(dp),allocatable :: work_space(:),work_spl(:)

! ***************************************************************************

!Set various terms to 0 in case not defined below
!HGH values
 rloc=zero ; cc(:)=zero
 nproj(1:psps%mpssoang)=0

!Read and write different lines of the pseudopotential file

 read (tmp_unit,*,err=10,iomsg=errmsg) rloc,nn,(cc(jj),jj=1,nn)
 write(message, '(a,f12.7)' ) ' rloc=',rloc
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')
 write(message, '(a,i1,a,4f12.7)' )' cc(1:',nn,')=',(cc(jj),jj=1,nn)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!Read the number of the non-local projectors
 read (tmp_unit,*,err=10,iomsg=errmsg) nnonloc
 if (nnonloc/=lmax+1) then
   write(message, '(a,a,a,a,i5,a,i5,a,a,a,a,i5)' ) ch10,&
&   ' psp10in : COMMENT -',ch10,&
&   '  input lmax=',lmax,'  does not agree with input nnonloc=',nnonloc,ch10,&
&   '  which has to be lmax+1.',ch10,&
&   '  Setting lmax to ',nnonloc-1
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   lmax=1
 end if
 ABI_ALLOCATE(rr,(0:lmax))
 ABI_ALLOCATE(hij,(0:lmax,3,3))
 ABI_ALLOCATE(kij,(0:lmax,3,3))
 rr(:)=zero; hij(:,:,:)=zero; kij(:,:,:)=zero

!Read and echo the coefficients of non-local projectors
 prjloop: do ll=0,lmax
   read (tmp_unit,*,err=10,iomsg=errmsg) rr(ll),nprl,(hij(ll,1,jj),jj=1,nprl)
   do ii=2,nprl
     read (tmp_unit,*,err=10,iomsg=errmsg) (hij(ll,ii,jj),jj=ii,nprl)
   end do
   nproj(ll+1)=nprl
   write(message, '(a,i3,a,f12.7,2a,3f12.7,2a,12x,2f12.7,2a,24x,f12.7)' )&
&   ' for angular momentum l =',ll,' r(l) =',rr(ll),ch10,&
&   '   h11, h12, h13 =', (hij(ll,1,jj),jj=1,3),ch10,&
&   '        h22, h23 =', (hij(ll,2,jj),jj=2,3),ch10,&
&   '             h33 =', (hij(ll,3,jj),jj=3,3)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   if (ll==0) cycle
   do ii=1,nprl
     read (tmp_unit,*,err=10,iomsg=errmsg) (kij(ll,ii,jj),jj=ii,nprl)
   end do
   write(message, '(a,3f12.7,2a,12x,2f12.7,2a,24x,f12.7)' )&
&   '   k11, k12, k13 =', (kij(ll,1,jj),jj=1,3),ch10,&
&   '        k22, k23 =', (kij(ll,2,jj),jj=2,3),ch10,&
&   '             k33 =', (kij(ll,3,jj),jj=3,3)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end do prjloop

 if(pspso/=0) then

!  MJV 10/2008: is this correct? For the normal HGH psp there are cases
!  where there are more SO projectors than SR ones! e.g. Pb with 12 electrons.
   do ll=1,lmax
     nproj(psps%mpsang+ll)=nproj(ll+1)
   end do

 end if

!Store the coefficients.
 psps%gth_params%set(ipsp)          = .true.
 psps%gth_params%psppar(0, :, ipsp) = (/ rloc, cc(1), cc(2), cc(3), cc(4), zero, zero /)
 do ii=1,4
   ll=ii-1
   if (ll>lmax) then
     psps%gth_params%psppar(ii,:,ipsp) = (/ zero, zero, zero, zero, zero, zero, zero /)
   else
     psps%gth_params%psppar(ii,:,ipsp) =&
&     (/ rr(ll), hij(ll,1,1), hij(ll,2,2), hij(ll,3,3), hij(ll,1,2), hij(ll,1,3), hij(ll,2,3) /)
   end if
 end do

!Additionnal wavelet parameters
 if (dtset%usewvl == 1) then
   call wvl_descr_psp_fill(psps%gth_params, ipsp, 0, int(psps%zionpsp(ipsp)), int(psps%znuclpsp(ipsp)), tmp_unit)
 end if

!Initialize array indlmn array giving l,m,n,ln,lm,s for i=lmn
 nso=1;if(pspso/=0) nso=2
 index=0;iln=0;indlmn(:,:)=0
 do nn=1,nso
   do ipsang=1+(nn-1)*(lmax+1),nn*lmax+1
     if (nproj(ipsang)>0) then
       ll=ipsang-(nn-1)*lmax-1
       do kk=1,nproj(ipsang)
         iln=iln+1
         do mm=1,2*ll*psps%useylm+1
           index=index+1
           indlmn(1,index)=ll
           indlmn(2,index)=mm-ll*psps%useylm-1
           indlmn(3,index)=kk
           indlmn(4,index)=ll*ll+(1-psps%useylm)*ll+mm
           indlmn(5,index)=iln
           indlmn(6,index)=nn
         end do
       end do
     end if
   end do
 end do

 ABI_ALLOCATE(dvlspl,(psps%mqgrid_ff))
!First, the local potential --  compute on q grid and fit spline
 call psp2lo(cc(1),cc(2),cc(3),cc(4),dvlspl,epsatm,psps%mqgrid_ff,psps%qgrid_ff,&
& vlspl(:,1),rloc,.true.,yp1,ypn,zion)
 ABI_DEALLOCATE(dvlspl)

!Fit spline to q^2 V(q) (Numerical Recipes subroutine)
 ABI_ALLOCATE(work_space,(psps%mqgrid_ff))
 ABI_ALLOCATE(work_spl,(psps%mqgrid_ff))
 call spline (psps%qgrid_ff,vlspl(:,1),psps%mqgrid_ff,yp1,ypn,work_spl)
 vlspl(:,2)=work_spl(:)
 ABI_DEALLOCATE(work_space)
 ABI_DEALLOCATE(work_spl)

!Second, compute KB energies and form factors and fit splines
 ekb(:)=zero

!Check if any nonlocal projectors are being used
 mproj=maxval(nproj)

 if (mproj>0) then

   ABI_ALLOCATE(ekb_sr,(psps%mpsang,mproj))
   ABI_ALLOCATE(ffspl_sr,(psps%mqgrid_ff,2,psps%mpsang,mproj))
   ABI_ALLOCATE(ekb_so,(psps%mpsang,mproj))
   ABI_ALLOCATE(ffspl_so,(psps%mqgrid_ff,2,psps%mpsang,mproj))

   call psp10nl(ekb_sr,ffspl_sr,hij,lmax,mproj,psps%mpsang,psps%mqgrid_ff,&
&   nproj,psps%qgrid_ff,rr)
   if(pspso/=0) then
     ABI_ALLOCATE(dummy_nproj,(psps%mpsang))
     dummy_nproj(1)=0
     do ll=1,lmax
       dummy_nproj(ll+1)=nproj(psps%mpsang+ll)
     end do
     call psp10nl(ekb_so,ffspl_so,kij,lmax,mproj,psps%mpsang,psps%mqgrid_ff,&
&     dummy_nproj,psps%qgrid_ff,rr)
     ABI_DEALLOCATE(dummy_nproj)
   end if

!  Convert ekb and ffspl
   iln0=0
   do jj=1,psps%lmnmax
     iln=indlmn(5,jj)
     if (iln>iln0) then
       iln0=iln
       if (indlmn(6,jj)<=1) then
         ekb(iln)=ekb_sr(1+indlmn(1,jj),indlmn(3,jj))
         ffspl(:,:,iln)=ffspl_sr(:,:,1+indlmn(1,jj),indlmn(3,jj))
       else
         ekb(iln)=ekb_so(1+indlmn(1,jj),indlmn(3,jj))
         ffspl(:,:,iln)=ffspl_so(:,:,1+indlmn(1,jj),indlmn(3,jj))
       end if
     end if
   end do

   ABI_DEALLOCATE(ekb_sr)
   ABI_DEALLOCATE(ffspl_sr)
   ABI_DEALLOCATE(ekb_so)
   ABI_DEALLOCATE(ffspl_so)
 end if

 ABI_DEALLOCATE(rr)
 ABI_DEALLOCATE(hij)
 ABI_DEALLOCATE(kij)

 return

 ! Handle IO error
 10 continue
 ABI_ERROR(errmsg)

end subroutine psp10in
!!***

!!****f* m_psp_hgh/psp10nl
!! NAME
!! psp10nl
!!
!! FUNCTION
!! Hartwigsen-Goedecker-Hutter nonlocal pseudopotential (from preprint of 1998).
!! Uses Gaussians for fully nonlocal form, analytic expressions.
!!
!! INPUTS
!!  hij(0:lmax,3,3)=factor defining strength of (max 3) projectors for each
!!   angular momentum channel l among 0, 1, ..., lmax
!!  lmax=maximum angular momentum
!!  mproj=maximum number of projectors in any channel
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=number of grid points for qgrid
!!  nproj(1:lmax+1)=number of projectors in any channel
!!  qgrid(mqgrid)=array of |G| values
!!  rr(0:lmax)=core radius for each 0<l<lmax channel (bohr)
!!
!! OUTPUT
!!  ekb(mpsang,mproj)=Kleinman-Bylander energies
!!  ffspl(mqgrid,2,mpssang,mproj)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projectors
!!
!! PARENTS
!!      m_psp_hgh
!!
!! CHILDREN
!!      spline,zhpev
!!
!! SOURCE

subroutine psp10nl(ekb,ffspl,hij,lmax,mproj,mpsang,mqgrid,nproj,qgrid,rr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmax,mproj,mpsang,mqgrid
!arrays
 integer,intent(in) :: nproj(mpsang)
 real(dp),intent(in) :: hij(0:lmax,3,3),qgrid(mqgrid),rr(0:lmax)
 real(dp),intent(out) :: ekb(mpsang,mproj),ffspl(mqgrid,2,mpsang,mproj)

!Local variables-------------------------------
!scalars
 integer :: info,ipack,iproj,iqgrid,jproj,ll,numproj
 real(dp) :: qmax,rrl
 character(len=500) :: message
 character :: jobz,uplo
!arrays
 real(dp) :: ap(2,9),rwork1(9),work1(2,9),ww(3),yp1j(3),ypnj(3)
 real(dp),allocatable :: ppspl(:,:,:,:),uu(:,:),work(:),zz(:,:,:)

! *************************************************************************

 ABI_ALLOCATE(ppspl,(mqgrid,2,mpsang,mproj))
 ABI_ALLOCATE(work,(mqgrid))

 qmax=qgrid(mqgrid)
 jobz='v'
 uplo='u'
 ekb(:,:)=zero

 lloop: do ll=0,lmax
   ap(:,:)=zero
   numproj=nproj(ll+1)

!  Fill up the matrix in packed storage
   prjloop: do jproj=1,numproj
     priloop: do iproj=1,jproj
       ipack=iproj+(jproj-1)*jproj/2
       if(mod((jproj-1)*jproj,2)/=0) then
         ABI_ERROR("odd")
       end if
       ap(1,ipack)=hij(ll,iproj,jproj)
     end do priloop
   end do prjloop

   if(numproj/=0)then

     ABI_ALLOCATE(uu,(numproj,numproj))
     ABI_ALLOCATE(zz,(2,numproj,numproj))

     if (numproj > 1) then
       call ZHPEV(jobz,uplo,numproj,ap,ww,zz,numproj,work1,rwork1,info)
       uu(:,:)=zz(1,:,:)
     else
       ww(1)=hij(ll,1,1)
       uu(1,1)=one
     end if

!    Initialization of ekb, and spline fitting

     if (ll==0) then ! s channel

       rrl=rr(0)
       do iproj=1,numproj
         ekb(1,iproj)=ww(iproj)*32.d0*(rrl**3)*(pi**2.5d0)/(4.d0*pi)**2
         if(iproj==1)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,1,1)=exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2)
           end do
           yp1j(1)=zero
           ypnj(1)=-(two_pi*rrl)**2*qmax*exp(-0.5d0*(two_pi*qmax*rrl)**2)
         else if(iproj==2)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,1,2)=2.0d0/sqrt(15.0d0)     &
&             *exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) &
&             *( 3.d0-(two_pi*qgrid(iqgrid)*rrl)**2 )
           end do
           yp1j(2)=zero
           ypnj(2)=2.0d0/sqrt(15.0d0)*(two_pi*rrl)**2*qmax &
&           *exp(-0.5d0*(two_pi*qmax*rrl)**2) * (-5.d0+(two_pi*qmax*rrl)**2)
         else if(iproj==3)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,1,3)=(4.0d0/3.0d0)/sqrt(105.0d0)*&
&             exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) * &
&             (15.0d0-10.0d0*(two_pi*qgrid(iqgrid)*rrl)**2 + &
&             (two_pi*qgrid(iqgrid)*rrl)**4)
           end do
           yp1j(3)=zero
           ypnj(3)=(4.0d0/3.0d0)/sqrt(105.0d0)*exp(-0.5d0*(two_pi*qmax*rrl)**2) * &
&           (two_pi*rrl)**2*qmax*(-35.0d0+14d0*(two_pi*qmax*rrl)**2-(two_pi*qmax*rrl)**4)
         end if
         call spline(qgrid,ppspl(:,1,1,iproj),mqgrid,&
&         yp1j(iproj),ypnj(iproj),ppspl(:,2,1,iproj))
       end do

     else if (ll==1) then ! p channel

       rrl=rr(1)
       do iproj=1,numproj
         ekb(2,iproj)=ww(iproj)*64.d0*(rrl**5)*(pi**2.5d0)/(4.d0*pi)**2
         if(iproj==1)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,2,1)=(1.0d0/sqrt(3.0d0))* &
&             exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) * (two_pi*qgrid(iqgrid))
           end do
           yp1j(1)=two_pi*(1.0d0/sqrt(3.0d0))
           ypnj(1)=-two_pi*((two_pi*qmax*rrl)**2-1.d0)*exp(-0.5d0*(two_pi*qmax*rrl)**2)*&
&           (1.0d0/sqrt(3.0d0))
         else if(iproj==2)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,2,2)=(2.0d0/sqrt(105.0d0))* &
&             exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) * &
&             (two_pi*qgrid(iqgrid))*(5.0d0-(two_pi*qgrid(iqgrid)*rrl)**2)
           end do
           yp1j(2)=(5.0d0*two_pi)*(2.0d0/sqrt(105.0d0))
           ypnj(2)=(2.0d0/sqrt(105.0d0))*two_pi*exp(-0.5d0*(two_pi*qmax*rrl)**2)* &
&           (-8*(two_pi*qmax*rrl)**2 + (two_pi*qmax*rrl)**4 + 5.0d0)
         else if(iproj==3)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,2,3)=(4.0d0/3.0d0)/sqrt(1155d0)*&
&             exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) * &
&             (two_pi*qgrid(iqgrid))*&
&             (35.0d0-14.0d0*(two_pi*qgrid(iqgrid)*rrl)**2+(two_pi*qgrid(iqgrid)*rrl)**4)
           end do
           yp1j(3)=(35.0d0*two_pi)*(4.0d0/3.0d0)/sqrt(1155.0d0)
           ypnj(3)=(4.0d0/3.0d0)/sqrt(1155.0d0)*two_pi*exp(-0.5d0*(two_pi*qmax*rrl)**2)* &
&           (35.0d0-77.0d0*(two_pi*qmax*rrl)**2+19.0d0*(two_pi*qmax*rrl)**4 - &
&           (two_pi*qmax*rrl)**6)
         end if
         call spline(qgrid,ppspl(:,1,2,iproj),mqgrid,&
&         yp1j(iproj),ypnj(iproj),ppspl(:,2,2,iproj))
       end do

     else if (ll==2) then ! d channel

!      If there is a third projector. Warning : only two projectors are allowed.
       if ( numproj>2 ) then
         write(message, '(3a)' )&
&         ' only two d-projectors are allowed ',ch10,&
&         ' Action: check your pseudopotential file.'
         ABI_ERROR(message)
       end if

       rrl=rr(2)
       do iproj=1,numproj
         ekb(3,iproj)=ww(iproj)*128.d0*(rrl**7)*(pi**2.5d0)/(4.d0*pi)**2
         if(iproj==1)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,3,1)=(1.0d0/sqrt(15.0d0))* &
&             exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) * (two_pi*qgrid(iqgrid))**2
           end do
           yp1j(1)=zero
           ypnj(1)=(1.0d0/sqrt(15.0d0))*(two_pi**2)*&
&           exp(-0.5d0*(two_pi*qmax*rrl)**2)*qmax*(2d0-(two_pi*qmax*rrl)**2)
         else if(iproj==2)then
           do iqgrid=1,mqgrid
             ppspl(iqgrid,1,3,2)=(2.0d0/3.0d0)/sqrt(105.0d0)* &
&             exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2) * &
&             ((two_pi*qgrid(iqgrid))**2)*(7.0d0-(two_pi*qgrid(iqgrid)*rrl)**2)
           end do
           yp1j(2)=zero
           ypnj(2)=(2.0d0/3.0d0)/sqrt(105.0d0)*exp(-0.5d0*(two_pi*qmax*rrl)**2)* &
&           qmax*(two_pi**2)*( (two_pi*qmax*rrl)**4 - 11.0d0*(two_pi*qmax*rrl)**2 + 14.0d0)
         end if
         call spline(qgrid,ppspl(:,1,3,iproj),mqgrid,&
&         yp1j(iproj),ypnj(iproj),ppspl(:,2,3,iproj))
       end do

     else if (ll==3) then ! f channel

!      If there is a second projector. Warning : only one projector is allowed.
       if ( numproj>1 ) then
         write(message, '(a,a,a)' )&
&         'only one f-projector is allowed ',ch10,&
&         'Action: check your pseudopotential file.'
         ABI_ERROR(message)
       end if

       rrl=rr(3)
       ekb(4,1)=ww(1)*(256.0d0/105.0d0)*(rrl**9)*(pi**2.5d0)/(4.d0*pi)**2
       do iqgrid=1,mqgrid
         ppspl(iqgrid,1,4,1)=(two_pi*qgrid(iqgrid))**3* &
&         exp(-0.5d0*(two_pi*qgrid(iqgrid)*rrl)**2)
       end do
!      Compute yp1,ypn=derivatives of f(q) at q=0, q=qgrid(mqgrid)
       yp1j(1)=zero
       ypnj(1)=(two_pi**3)*qmax**2*exp(-0.5d0*(two_pi*qmax*rrl)**2)*&
&       (3.0d0-(two_pi*qmax*rrl)**2)
!      Fit spline to get second derivatives by spline fit
       call spline(qgrid,ppspl(:,1,4,1),mqgrid,&
&       yp1j(1),ypnj(1),ppspl(:,2,4,1))

     else
       ABI_ERROR("lmax>3?")
     end if

!    Linear combination using the eigenvectors
     ffspl(:,:,ll+1,:)=zero
     do jproj=1,numproj
       do iproj=1,numproj
         do iqgrid=1,mqgrid
           ffspl(iqgrid,1:2,ll+1,jproj)=ffspl(iqgrid,1:2,ll+1,jproj) &
&           +uu(iproj,jproj)*ppspl(iqgrid,1:2,ll+1,iproj)
         end do
       end do
     end do

     ABI_DEALLOCATE(uu)
     ABI_DEALLOCATE(zz)

!    End condition on numproj(/=0)
   end if

 end do lloop

 ABI_DEALLOCATE(ppspl)
 ABI_DEALLOCATE(work)

end subroutine psp10nl
!!***

end module m_psp_hgh
!!***
