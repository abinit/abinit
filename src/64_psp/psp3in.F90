!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp3in
!! NAME
!! psp3in
!!
!! FUNCTION
!! Initialize pspcod=3 pseudopotentials (HGH psps PRB58,3641(1998)):
!! continue to read the file, then compute the corresponding
!! local and non-local potentials.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR, FD, PT)
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
!!      psp2lo,psp3nl,spline,wrtout,wvl_descr_psp_fill
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp3in(dtset, ekb, epsatm, ffspl, indlmn, ipsp, lmax, nproj, psps, pspso, vlspl, zion)

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
#define ABI_FUNC 'psp3in'
 use interfaces_14_hidewrite
 use interfaces_43_wvl_wrappers
 use interfaces_64_psp, except_this_one => psp3in
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

!DEBUG
!write(std_out,*) ' psp3in : enter '
!ENDDEBUG

!Read and write different lines of the pseudopotential file

 read (tmp_unit,*,err=10,iomsg=errmsg) rloc,cc1,cc2,cc3,cc4
 write(message, '(a,f12.7)' ) ' rloc=',rloc
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')
 write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )&
& ' cc1 =',cc1,'; cc2 =',cc2,'; cc3 =',cc3,'; cc4 =',cc4
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!For the time being, the s state line must be present and is read,
!even for local pseudopotentials (zero must appear)
 read (tmp_unit,*,err=10,iomsg=errmsg) rrs,h11s,h22s,h33s
 write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )&
& ' rrs =',rrs,'; h11s=',h11s,'; h22s=',h22s,'; h33s=',h33s
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 if (lmax > 0) then

   read (tmp_unit,*,err=10,iomsg=errmsg) rrp,h11p,h22p,h33p
   write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )&
&   ' rrp =',rrp,'; h11p=',h11p,'; h22p=',h22p,'; h33p=',h33p
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   read (tmp_unit,*,err=10,iomsg=errmsg) k11p,k22p,k33p
   write(message, '(a,f12.7,a,f12.7,a,f12.7)' )&
&   '                    k11p=',k11p,'; k22p=',k22p,'; k33p=',k33p
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

 end if

 if (lmax > 1) then

   read (tmp_unit,*,err=10,iomsg=errmsg) rrd,h11d,h22d,h33d
   write(message, '(a,f12.7,a,f12.7,a,f12.7,a,f12.7)' )&
&   ' rrd =',rrd,'; h11d=',h11d,'; h22d=',h22d,'; h33d=',h33d
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

   read (tmp_unit,*,err=10,iomsg=errmsg) k11d,k22d,k33d
   write(message, '(a,f12.7,a,f12.7,a,f12.7)' )&
&   '                    k11d=',k11d,'; k22d=',k22d,'; k33d=',k33d
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
 call psp2lo(cc1,cc2,cc3,cc4,dvlspl,epsatm,psps%mqgrid_ff,psps%qgrid_ff,&
& vlspl(:,1),rloc,.true.,yp1,ypn,zion)
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
 MSG_ERROR(errmsg)

end subroutine psp3in
!!***

!!****f* ABINIT/psp3nl
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
!!      psp3in
!!
!! CHILDREN
!!      spline,zhpev
!!
!! SOURCE

subroutine psp3nl(ekb,ffspl,h11s,h22s,h33s,h11p,h22p,h33p,h11d,h22d,&
&                  h33d,h11f,mproj,mpsang,mqgrid,qgrid,rrd,rrf,rrp,rrs)

 use defs_basis
 use m_splines
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp3nl'
!End of the abilint section

 implicit none

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

!DEBUG
!write(std_out,*)' psp3nl : enter '
!stop
!ENDDEBUG

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

!  End condition on nproj(/=0)
 end if

!DEBUG
!write(std_out,*)' psp3nl : after s channel '
!stop
!ENDDEBUG

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

!  End condition on nproj(/=0)
 end if

!DEBUG
!write(std_out,*)' psp3nl : after p channel '
!stop
!ENDDEBUG

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
&   '  Action : check your pseudopotential file.'
   MSG_ERROR(message)
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

!  End condition on nproj(/=0)
 end if

!DEBUG
!write(std_out,*)' psp3nl : after d channel '
!stop
!ENDDEBUG

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

!DEBUG
!write(std_out,*)' psp3nl : exit '
!stop
!ENDDEBUG

end subroutine psp3nl
!!***
