!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp2in
!! NAME
!! psp2in
!!
!! FUNCTION
!! Initialize pspcod=2 pseudopotentials (GTH format):
!! continue to read the file, then compute the corresponding
!! local and non-local potentials.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!      pspatm
!!
!! CHILDREN
!!      psp2lo,psp2nl,spline,wrtout,wvl_descr_psp_fill
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp2in(dtset,ekb,epsatm,ffspl,indlmn,ipsp,lmax,nproj,psps,vlspl,dvlspl,zion)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_splines
 use m_profiling_abi
 use m_errors
#if defined HAVE_BIGDFT
 use BigDFT_API, only: atomic_info
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp2in'
 use interfaces_14_hidewrite
 use interfaces_43_wvl_wrappers
 use interfaces_64_psp, except_this_one => psp2in
!End of the abilint section

 implicit none

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
 write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )&
& '  cc1=',cc1,'; cc2=',cc2,'; cc3=',cc3,'; cc4=',cc4
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 read (tmp_unit,*, err=10, iomsg=errmsg) rrs,h1s,h2s
 write(message, '(a,f12.7,a,f12.7,a,f12.7)' )&
& '  rrs=',rrs,'; h1s=',h1s,'; h2s=',h2s
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 read (tmp_unit,*, err=10, iomsg=errmsg) rrp,h1p
 write(message, '(a,f12.7,a,f12.7)' )&
& '  rrp=',rrp,'; h1p=',h1p
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
     MSG_ERROR(message)
   end if
   nproj(2)=1
   if (lmax<1) then
     write(message, '(a,i5,a,e12.4,a,a,a,a)' )&
&     'Input lmax=',lmax,' disagree with input h1p=',h1p,'.',&
&     'Your pseudopotential is incoherent.',ch10,&
&     'Action : correct your pseudopotential file.'
     MSG_ERROR(message)
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

!DEBUG
!write(std_out,*)' psp2in : after psp2lo '
!stop
!ENDDEBUG

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
 MSG_ERROR(errmsg)

end subroutine psp2in
!!***
