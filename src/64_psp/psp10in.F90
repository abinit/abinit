!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp10in
!! NAME
!! psp10in
!!
!! FUNCTION
!! Initialize pspcod=10 pseudopotentials (formalism is the same as in HGH psps
!! PRB58,3641(1998), but the full h and k matrices are read, allowing for using
!! also subsequent developments such as Theor. Chem. Acc. 114, 145 (2005)):
!! continue to read the file, then compute the corresponding
!! local and non-local potentials.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, FD, SC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!      pspatm
!!
!! CHILDREN
!!      psp10nl,psp2lo,spline,wrtout,wvl_descr_psp_fill
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp10in(dtset, ekb, epsatm, ffspl, indlmn, ipsp, lmax, nproj, psps, pspso, vlspl, zion)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_splines
 use m_errors
#if defined HAVE_BIGDFT
  use BigDFT_API, only: atomic_info
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp10in'
 use interfaces_14_hidewrite
 use interfaces_43_wvl_wrappers
 use interfaces_64_psp, except_this_one => psp10in
!End of the abilint section

 implicit none

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

!DEBUG
!write(std_out,*) ' psp10in : enter '
!ENDDEBUG

!Read and write different lines of the pseudopotential file

 read (tmp_unit,*,err=10,iomsg=errmsg) rloc,nn,(cc(jj),jj=1,nn)
 write(message, '(a,f12.7)' ) ' rloc=',rloc
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')
 write(message, '(a,i1,a,4f12.7)' )&
& ' cc(1:',nn,')=',(cc(jj),jj=1,nn)
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

!DEBUG
!write(std_out,*)' psp10in : after psp2lo '
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
 MSG_ERROR(errmsg)

end subroutine psp10in
!!***
