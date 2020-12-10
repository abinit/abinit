!!****m* ABINIT/m_psp6
!! NAME
!!  m_psp6
!!
!! FUNCTION
!! Initialize pspcod=6 (Pseudopotentials from the fhi98pp code):
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2020 ABINIT group (XG, AF, GJ,FJ,MT, DRH)
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

module m_psp6

 use defs_basis
 use m_splines
 use m_errors
 use m_abicore

 use m_numeric_tools,  only : smooth, ctrap
 use m_psptk,          only : psp5lo, psp5nl, cc_derivatives

 implicit none

 private
!!***

 public :: psp6in
!!***

contains
!!***

!!****f* m_psp6/psp6in
!! NAME
!! psp6in
!!
!! FUNCTION
!! Initialize pspcod=6 (Pseudopotentials from the fhi98pp code):
!! continue to read the corresponding file, then compute the
!! local and non-local potentials.
!!
!! INPUTS
!!  lloc=angular momentum choice of local pseudopotential
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  lnmax=max. number of (l,n) components over all type of psps
!!  mmax=maximum number of points in real space grid in the psp file
!!   angular momentum of nonlocal pseudopotential
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=dimension of q (or G) grid for arrays.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  optnlxccc=option for nl XC core correction (input variable)
!!  positron=0 if electron GS calculation
!!          1 if positron GS calculation
!!          2 if electron GS calculation in presence of the positron
!!  qgrid(mqgrid)=values of q (or |G|) on grid from 0 to qmax
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  zion=nominal valence of atom as specified in psp file
!!  znucl=nuclear number of atom as specified in psp file
!!
!! OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy,
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for each (l,n)
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r} dr]$ (hartree)
!!  ffspl(mqgrid,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  nproj(mpsang)=number of projection functions for each angular momentum
!!  qchrg is not used, and could be suppressed later
!!  vlspl(mqgrid,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  xcccrc=XC core correction cutoff radius (bohr)
!!  xccc1d(n1xccc,6)=1D core charge function and five derivatives, from psp file
!!
!! PARENTS
!!      m_pspini
!!
!! CHILDREN
!!      cc_derivatives
!!
!! SOURCE

subroutine psp6in(ekb,epsatm,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&
&                  mmax,mpsang,mqgrid,nproj,n1xccc,optnlxccc,positron,qchrg,qgrid,&
&                  useylm,vlspl,xcccrc,xccc1d,zion,znucl)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lloc,lmax,lmnmax,lnmax,mmax,mpsang,mqgrid,n1xccc
 integer,intent(in) :: optnlxccc,positron,useylm
 real(dp),intent(in) :: zion,znucl
 real(dp),intent(out) :: epsatm,qchrg,xcccrc
!arrays
 integer,intent(out) :: indlmn(6,lmnmax),nproj(mpsang)
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(out) :: ekb(lnmax),vlspl(mqgrid,2) !vz_i
 real(dp),intent(inout) :: ffspl(mqgrid,2,lnmax) !vz_i
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i

!Local variables-------------------------------
!scalars
 integer :: ii,index,ipsang,irad,jj,jpsang,mm,mmax2
 real(dp) :: al,al_announced,amesh,fchrg,ratio,rchrg,yp1,ypn
 character(len=3) :: testxc
 character(len=500) :: message,errmsg
!arrays
 real(dp),allocatable :: ekb_tmp(:),ffspl_tmp(:,:,:),rad(:),vloc(:)
!real(dp),allocatable :: radbis
 real(dp),allocatable :: vpspll(:,:),wfll(:,:),work_space(:),work_spl(:)

! ***************************************************************************

!File format of formatted fhi psp input, as adapted for use
!by the ABINIT code (the 3 first lines
!have already been read in calling -pspatm- routine) :

!(1) title (character) line
!(2) znucl,zion,pspdat
!(3) pspcod,pspxc,lmax,lloc,mmax,r2well
!(4) rchrg,fchrg,qchrg
!Note : prior to version 2.2, this 4th line started with  4--  ,
!and no core-correction was available.
!(5)-(18) -empty-
!(19) mmax, amesh ( mesh increment r(m+1)/r(m) )
!Then, for ll=0,lmax :
!for  irad=1,mmax  : irad, r(irad), upsp(irad,ll), vpsp(irad,ll)

 read (tmp_unit, '(a3)', err=10, iomsg=errmsg) testxc
 if(testxc/='4--')then
   backspace(tmp_unit)
   read (tmp_unit,*, err=10, iomsg=errmsg) rchrg,fchrg,qchrg
   write(message, '(3f20.14,t64,a)' ) rchrg,fchrg,qchrg,'rchrg,fchrg,qchrg'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 else
   write(message, '(a)' ) '  No XC core correction.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   rchrg=zero ; fchrg=zero ; qchrg=zero
 end if
 do ii=5,18
   read(tmp_unit,*, err=10, iomsg=errmsg)
 end do

 if (positron==1.and.abs(fchrg)<=tol14) then
   write(message,'(5a)')&
&   'You can only perform positronic ground-state calculations (positron=1)',ch10,&
&   'using fhi pseudopotentials with a core density (fchrg>0)',ch10,&
&   'Action: change your psp file (add fchrg>0).'
   ABI_ERROR(message)
 end if
!--------------------------------------------------------------------
!Will now proceed at the reading of pots and wfs

!rad(:)=radial grid r(i)
!vpspll(:,1),...,vpspll(:,4)=nonlocal pseudopotentials
!wfll(:,1),...,wfll(:,4)=reference config. wavefunctions

 ABI_MALLOC(rad,(mmax))
 ABI_MALLOC(vpspll,(mmax,mpsang))
 ABI_MALLOC(wfll,(mmax,mpsang))

!Read atomic pseudopotential for each l, filling up arrays vpspll
!and wfll. Also set up rad array (actually read more than once)
!Note: put each l into vpspll(:,l+1)
 do ipsang=1,lmax+1
   nproj(ipsang)=1
   read(tmp_unit,*, err=10, iomsg=errmsg)mmax2,amesh
   if(ipsang==1)then
     write(message, '(f10.6,t20,a)' ) amesh,' amesh (Hamman grid)'
     al_announced=log(amesh)
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if
   do irad=1,mmax
     read(tmp_unit,*, err=10, iomsg=errmsg)jj,rad(irad),wfll(irad,ipsang),vpspll(irad,ipsang)
   end do
 end do


!Generate core charge function and derivatives, if needed
 if(fchrg>tol14)then

   if (positron==1) then
     call psp6cc(mmax,n1xccc,rchrg,xccc1d,znucl,vh_tnzc=vpspll(:,lloc+1))
   else if(optnlxccc==1)then
     call psp6cc(mmax,n1xccc,rchrg,xccc1d,znucl)
   else if(optnlxccc==2)then
     call psp6cc_drh(mmax,n1xccc,rchrg,xccc1d)
   end if

! The core charge function for pspcod=6 becomes zero beyond rchrg.
! Thus xcccrc must be set equal to rchrg.
   xcccrc=rchrg
 else
   xccc1d(:,:)=zero
   xcccrc=zero
 end if

!Compute in real(dp) al : the announced amesh is inaccurate.
 ratio=rad(mmax)/rad(1)
 al=log(ratio)/dble(mmax-1)

!vloc(:)=Vlocal(r), lloc=0, 1, or 2 or -1 for avg.
 ABI_MALLOC(vloc,(mmax))
!Copy appropriate nonlocal psp for use as local one
 vloc( 1:mmax ) = vpspll( 1:mmax , lloc+1 )

!--------------------------------------------------------------------
!Carry out calculations for local (lloc) pseudopotential.
!Obtain Fourier transform (1-d sine transform) to get q^2 V(q).

 call psp5lo(al,epsatm,mmax,mqgrid,qgrid,&
& vlspl(:,1),rad,vloc,yp1,ypn,zion)

!Fit spline to q^2 V(q) (Numerical Recipes subroutine)
 ABI_MALLOC(work_space,(mqgrid))
 ABI_MALLOC(work_spl,(mqgrid))
 call spline (qgrid,vlspl(:,1),mqgrid,yp1,ypn,work_spl)
 vlspl(:,2)=work_spl(:)
 ABI_FREE(work_space)
 ABI_FREE(work_spl)

!--------------------------------------------------------------------
!Take care of non-local part

 ABI_MALLOC(ekb_tmp,(mpsang))
 ABI_MALLOC(ffspl_tmp,(mqgrid,2,mpsang))

!Zero out all Kleinman-Bylander energies to initialize
 ekb_tmp(:)=zero
 ekb(:)=zero

!Allow for option of no nonlocal corrections (lloc=lmax=0)
 if (lloc==0.and.lmax==0) then
   write(message, '(a,f5.1)' ) ' Note: local psp for atom with Z=',znucl
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 else

!  ----------------------------------------------------------------------
!  Compute KB form factors and fit splines

   call psp5nl(al,ekb_tmp,ffspl_tmp,lmax,mmax,mpsang,mqgrid,qgrid,rad,vloc, vpspll,wfll)

 end if

 jj=0;index=0;indlmn(:,:)=0
 do ipsang=1,lmax+1
!  nproj had been set at 1, by default
   if(abs(ekb_tmp(ipsang))<tol10)then
     nproj(ipsang)=0
   end if
!  Possible values for nproj in this routine : 0 or 1.
   if(nproj(ipsang)==1)then
     if (useylm==1) then
       jj=jj+1
       do mm=1,2*ipsang-1
         index=index+1
         indlmn(1,index)=ipsang-1
         indlmn(2,index)=mm-ipsang
         indlmn(3,index)=1
         indlmn(4,index)=mm+(ipsang-1)*(ipsang-1)
         indlmn(5,index)=jj
         indlmn(6,index)=1
       end do
     else
       jj=jj+1
       index=index+1
       indlmn(1,index)=ipsang-1
       indlmn(2,index)=0
       indlmn(3,index)=1
       indlmn(4,index)=ipsang+(ipsang-1)*(ipsang-1)
       indlmn(5,index)=jj
       indlmn(6,index)=1
     end if
   end if
 end do
!Transfer ekb and ffspl to their definitive location
 jpsang=1
 do ipsang=1,lmax+1
   if(nproj(ipsang)/=0)then
     ekb(jpsang)=ekb_tmp(ipsang)
     ffspl(:,:,jpsang)=ffspl_tmp(:,:,ipsang)
     jpsang=jpsang+1
     if(jpsang>lnmax)then
       write(message,'(3a,2i6)')&
&       'Problem with the dimension of the ekb and ffspl arrays.',ch10,&
&       'ipsang,lnmax=',ipsang,lnmax
     end if
   end if
 end do

 ABI_FREE(ekb_tmp)
 ABI_FREE(ffspl_tmp)
 ABI_FREE(vpspll)
 ABI_FREE(rad)
 ABI_FREE(vloc)
 ABI_FREE(wfll)

 return

 ! Handle IO error
 10 continue
 ABI_ERROR(errmsg)

end subroutine psp6in
!!***

!!****f* m_psp6/psp6cc
!! NAME
!! psp6cc
!!
!! FUNCTION
!! Compute the core charge density, for use in the XC core
!! correction, following the function definition valid
!! for the format 6 of pseudopotentials.
!!
!! INPUTS
!!  mmax=maximum number of points in real space grid in the psp file
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  rchrg=cut-off radius for the core density
!!  znucl=nuclear number of atom as specified in psp file
!!
!! OUTPUT
!!  xccc1d(n1xccc,6)= 1D core charge function and its five first derivatives
!!  Optional output:
!!    vh_tnzc(mmax) = Hartree potential induced by density tild_[n_Z+n_core]
!!                    (pseudized [n_Z+n_core], where n_Z=ions, n_core=core electrons)
!!                    using a simple pseudization scheme
!!
!! PARENTS
!!      m_psp6
!!
!! CHILDREN
!!      cc_derivatives
!!
!! SOURCE

subroutine psp6cc(mmax,n1xccc,rchrg,xccc1d,znucl,&
&                 vh_tnzc) ! optional argument

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mmax,n1xccc
 real(dp),intent(in) :: rchrg,znucl
!arrays
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i
 real(dp),intent(out),optional :: vh_tnzc(mmax)

!Local variables-------------------------------
!scalars
 integer :: i1xccc,irad
 real(dp) :: der1,dern
 character(len=500) :: errmsg
!arrays
 real(dp),allocatable :: ff(:),ff1(:),ff2(:),ff3(:),gg(:),gg1(:),gg2(:),gg3(:)
 real(dp),allocatable :: gg4(:),nc(:),rad(:),work(:),xx(:)

!**********************************************************************

 ABI_MALLOC(ff,(mmax))
 ABI_MALLOC(ff1,(mmax))
 ABI_MALLOC(ff2,(mmax))
 ABI_MALLOC(ff3,(mmax))
 ABI_MALLOC(rad,(mmax))
 ABI_MALLOC(gg,(n1xccc))
 ABI_MALLOC(gg1,(n1xccc))
 ABI_MALLOC(gg2,(n1xccc))
 ABI_MALLOC(gg3,(n1xccc))
 ABI_MALLOC(gg4,(n1xccc))
 ABI_MALLOC(work,(n1xccc))
 ABI_MALLOC(xx,(n1xccc))

!read from pp file the model core charge (ff) and first (ff1) and
!second (ff2) derivative on logarithmic mesh mmax; rad is the radial grid
!the input functions contain the 4pi factor, it must be rescaled.

 do irad=1,mmax
   read(tmp_unit,*, err=10, iomsg=errmsg) rad(irad),ff(irad),ff1(irad),ff2(irad)
   ff(irad)=ff(irad)/four_pi
   ff1(irad)=ff1(irad)/four_pi
   ff2(irad)=ff2(irad)/four_pi
 end do

!Optional output: VHartree(tild_[n_Z+n_core])
 if (present(vh_tnzc)) then
   ABI_MALLOC(nc,(mmax))
   nc=ff ! n_core
   call psden(1,ff,mmax,nc,rchrg,rad,ff1=ff1,ff2=ff2)
   call vhtnzc(ff,rchrg,vh_tnzc,mmax,rad,znucl)
   ABI_FREE(nc)
 end if

 rad(1)=zero

!calculate third derivative ff3 on logarithmic grid
 der1=ff2(1)
 dern=ff2(mmax)
 call spline(rad,ff1,mmax,der1,dern,ff3)

!generate uniform mesh xx in the box cut by rchrg:

 do i1xccc=1,n1xccc
   xx(i1xccc)=(i1xccc-1)* rchrg/dble(n1xccc-1)
 end do

!now interpolate core charge and derivatives on the uniform grid
!core charge, input=ff,  output=gg
 call splint(mmax,rad,ff,ff2,n1xccc,xx,gg)

!first derivative input=ff1, output=gg1
 call splint(mmax,rad,ff1,ff3,n1xccc,xx,gg1)

!normalize gg1
 gg1(:)=gg1(:)*rchrg

!now calculate second to fourth derivative by forward differences
!to avoid numerical noise uses a smoothing function

 call smooth(gg1,n1xccc,10)

 gg2(n1xccc)=zero
 do i1xccc=1,n1xccc-1
   gg2(i1xccc)=(gg1(i1xccc+1)-gg1(i1xccc))*dble(n1xccc-1)
 end do

 call smooth(gg2,n1xccc,10)

 gg3(n1xccc)=zero
 do i1xccc=1,n1xccc-1
   gg3(i1xccc)=(gg2(i1xccc+1)-gg2(i1xccc))*dble(n1xccc-1)
 end do

 call smooth(gg3,n1xccc,10)

 gg4(n1xccc)=zero
 do i1xccc=1,n1xccc-1
   gg4(i1xccc)=(gg3(i1xccc+1)-gg3(i1xccc))*dble(n1xccc-1)
 end do

 call smooth(gg4,n1xccc,10)

!write on xcc1d
 xccc1d(:,1)=gg(:)
 xccc1d(:,2)=gg1(:)
 xccc1d(:,3)=gg2(:)
 xccc1d(:,4)=gg3(:)
 xccc1d(:,5)=gg4(:)

!WARNING : fifth derivative not yet computed
 xccc1d(:,6)=zero

!note: the normalization condition is the following:
!4pi rchrg /dble(n1xccc-1) sum xx^2 xccc1d(:,1) = qchrg
!
!norm=zero
!do i1xccc=1,n1xccc
!norm = norm + four_pi*rchrg/dble(n1xccc-1)*&
!&             xx(i1xccc)**2*xccc1d(i1xccc,1)
!end do
!write(std_out,*) ' norm=',norm
!
!write(std_out,*)' psp1cc : output of core charge density and derivatives '
!write(std_out,*)'   xx          gg           gg1  '
!do i1xccc=1,n1xccc
!write(10, '(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,1),xccc1d(i1xccc,2)
!end do
!write(std_out,*)'   xx          gg2          gg3  '
!do i1xccc=1,n1xccc
!write(11, '(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,3),xccc1d(i1xccc,4)
!end do
!write(std_out,*)'   xx          gg4          gg5  '
!do i1xccc=1,n1xccc
!write(12, '(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,5),xccc1d(i1xccc,6)
!end do
!write(std_out,*)' psp1cc : debug done, stop '
!stop
!ENDDEBUG

 ABI_FREE(ff)
 ABI_FREE(ff1)
 ABI_FREE(ff2)
 ABI_FREE(ff3)
 ABI_FREE(gg)
 ABI_FREE(gg1)
 ABI_FREE(gg2)
 ABI_FREE(gg3)
 ABI_FREE(gg4)
 ABI_FREE(rad)
 ABI_FREE(work)
 ABI_FREE(xx)

 return

 ! Handle IO error
 10 continue
 ABI_ERROR(errmsg)

end subroutine psp6cc
!!***

!!****f* m_psp6/psden
!! NAME
!! psden
!!
!! FUNCTION
!! Calculate a pseudo-density from an original density on a radial grid (regular or logarithmic)
!!
!! INPUTS
!!  ilog=1 if grid is logarithmic, else 0
!!  mesh= dimension of nc
!!  nc(mesh)= density to be pseudized
!!  rc= cut-off radius
!!  rad(mesh) = radial mesh
!!
!! OUTPUT
!!  ff(mesh)= pseudized density
!!
!!SIDE EFFECTS
!!  Optional:
!!    ff1(mesh)= 1st derivative of pseudo density (only r<rc modified)
!!    ff2(mesh)= 2nd derivative of pseudo density (only r<rc modified)
!!
!! NOTES
!!    ff=exp(-(a+b.r^2+c.r^4))
!!
!! PARENTS
!!      m_psp6
!!
!! CHILDREN
!!      cc_derivatives
!!
!! SOURCE

subroutine psden(ilog,ff,mesh,nc,rc,rad,ff1,ff2)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ilog,mesh
 real(dp),intent(in) :: rc
!arrays
 real(dp),intent(in) :: nc(mesh),rad(mesh)
 real(dp),intent(out) :: ff(mesh)
 real(dp),intent(inout),optional :: ff1(mesh),ff2(mesh)

!Local variables-------------------------------
!scalars
 integer :: ii,nc1
 real(dp) :: aa,aa1,aa2,bb,cc,c1,c3,f0,f0p,norm1,norm2,rc1,step
!arrays
 real(dp),allocatable :: fpir(:),gg(:)

! *************************************************************************

 rc1=rc/four

 ABI_MALLOC(fpir,(mesh))
 fpir(1:mesh)=four_pi*rad(1:mesh)**2
 if (ilog==1) fpir(1:mesh)=fpir(1:mesh)*rad(1:mesh)

 if (ilog==0) then
   step=rad(2)-rad(1)
   nc1=int(rc1/step)+1
   rc1=(nc1-1)*step
 else if (ilog==1) then
   step=log(rad(2)/rad(1))
   nc1=int(log(rc1/rad(1))/step)+1
   rc1=rad(nc1)
 end if
 ff(1:nc1)=nc(1:nc1)*fpir(1:nc1)
 call ctrap(nc1,ff(1:nc1),step,c3)
 if (ilog==1) c3=c3+half*ff(1)
 f0=nc(nc1);c1=-log(f0)
 f0p=half*(nc(nc1+1)-nc(nc1-1))/step

 ii=0;aa1=zero;norm1=c3+one
 do while (norm1>c3.and.ii<100)
   ii=ii+1;aa1=aa1+one
   aa=c1-aa1*rc1**4+rc1*(f0p/f0+four*aa1*rc1**3)*half
   bb=-half*(f0p/f0+four*aa1*rc1**3)/rc1
   ff(1:nc1)=fpir(1:nc1)*exp(-aa-bb*rad(1:nc1)**2-aa1*rad(1:nc1)**4)
   call ctrap(nc1,ff(1:nc1),step,norm1)
   if (ilog==1) norm1=norm1+half*ff(1)
 end do
 if (ii==100) then
   ABI_ERROR('Big pb 1 in psden !')
 end if

 ii=0;aa2=zero;norm2=c3-one
 do while (norm2<c3.and.ii<100)
   ii=ii+1;aa2=aa2-one
   aa=c1-aa2*rc1**4+rc1*(f0p/f0+four*aa2*rc1**3)*half
   bb=-half*(f0p/f0+four*aa2*rc1**3)/rc1
   ff(1:nc1)=fpir(1:nc1)*exp(-aa-bb*rad(1:nc1)**2-aa2*rad(1:nc1)**4)
   call ctrap(nc1,ff(1:nc1),step,norm2)
   if (ilog==1) norm2=norm2+half*ff(1)
 end do
 if (ii==100) then
   ABI_ERROR('Big pb 2 in psden !')
 end if

 do while (abs(norm2-c3)>tol10)

   cc=(aa1+aa2)*half
   aa=c1-cc*rc1**4+rc1*(f0p/f0+four*cc*rc1**3)*half
   bb=-half*(f0p/f0+four*cc*rc1**3)/rc1
   ff(1:nc1)=fpir(1:nc1)*exp(-aa-bb*rad(1:nc1)**2-cc*rad(1:nc1)**4)
   call ctrap (nc1,ff(1:nc1),step,norm2)
   if (ilog==1) norm2=norm2+half*ff(1)
   if ((norm1-c3)*(norm2-c3)>zero) then
     aa1=cc
     norm1=norm2
   else
     aa2=cc
   end if

 end do ! while

 ff(1)=exp(-aa);if (ilog==1) ff(1)=ff(1)*exp(-bb*rad(1)**2-cc*rad(1)**4)
 ff(2:nc1)=ff(2:nc1)/fpir(2:nc1)
 if (nc1<mesh) ff(nc1+1:mesh)=nc(nc1+1:mesh)
 if (present(ff1)) ff1(1:nc1)=-(two*bb*rad(1:nc1)+four*cc*rad(1:nc1)**3)*ff(1:nc1)
 if (present(ff2)) ff2(1:nc1)=-(two*bb+12.0_dp*cc*rad(1:nc1)**2)*ff(1:nc1) &
& +(two*bb*rad(1:nc1)+four*cc*rad(1:nc1)**3)**2*ff(1:nc1)

 ABI_MALLOC(gg,(mesh))
 gg(1:mesh)=fpir(1:mesh)*ff(1:mesh)
 call ctrap(mesh,gg(1:mesh),step,norm1)
 if (ilog==1) norm1=norm1+half*gg(1)
 !write(std_out,*) 'psden: tild_nc integral= ',norm1
 ABI_FREE(gg)

 ABI_FREE(fpir)

end subroutine psden
!!***

!!****f* m_psp6/vhtnzc
!! NAME
!! vhtnzc
!!
!! FUNCTION
!! Compute VHartree(tild[n_Z+n_core]) from input ncore
!!
!! INPUTS
!!  mesh=dimension of radial mesh
!!  nc= core density (to be pseudized)
!!  rad(mesh)=radial mesh
!!  rc=cut-off radius
!!  znucl=nuclear number of atom as specified in psp file
!!
!! OUTPUT
!!  vhtnzc(mesh) = hartree potential induced by density tild[n_Z+n_core] (pseudo core density + nucleus)
!!
!! PARENTS
!!      m_psp6
!!
!! CHILDREN
!!      cc_derivatives
!!
!! SOURCE

subroutine vhtnzc(nc,rc,vh_tnzc,mesh,rad,znucl)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mesh
 real(dp),intent(in) :: znucl
 real(dp),intent(in) :: rc
!arrays
 real(dp),intent(in) :: nc(mesh),rad(mesh)
 real(dp),intent(out) :: vh_tnzc(mesh)

!Local variables-------------------------------
!scalars
 integer :: ir,nc1
 real(dp) :: gnorm,rc1,step,yp1,yp2,yp3
!arrays
 real(dp),allocatable :: den1(:),den2(:),den3(:),den4(:),nzc(:),rvhn(:),shapefunc(:)

! *************************************************************************

 rc1=rc/four

 step=log(rad(2)/rad(1))
 nc1=int(log(rc1/rad(1))/step)+1
 rc1=rad(nc1)

 ABI_MALLOC(shapefunc,(mesh))
 shapefunc(1)=one
 shapefunc(2:nc1)=(sin(pi*rad(2:nc1)/rc1)/(pi*rad(2:nc1)/rc1))**2
 if (nc1<mesh) shapefunc(nc1+1:mesh)=zero

 ABI_MALLOC(den1,(mesh))
 den1(1:mesh)=four_pi*shapefunc(1:mesh)*rad(1:mesh)**3
 call ctrap(mesh,den1,step,gnorm)
 gnorm =one/gnorm
 ABI_FREE(den1)

 ABI_MALLOC(nzc,(mesh))
 nzc(1:mesh)=four*pi*nc(1:mesh)*rad(1:mesh)**2-four_pi*shapefunc(1:mesh)*rad(1:mesh)**2*znucl*gnorm
 ABI_FREE(shapefunc)

 ABI_MALLOC(rvhn,(mesh))
 rvhn(1)=zero

 ABI_MALLOC(den1,(mesh))
 ABI_MALLOC(den2,(mesh))
 ABI_MALLOC(den3,(mesh))
 ABI_MALLOC(den4,(mesh))

 den1(1)=zero;den2(1)=zero
 do ir=2,mesh
   den1(ir)= rad(ir)*nzc(ir)
   den2(ir)= den1(ir)/rad(ir)
 end do

!For first few points do stupid integral
 den3(1)=zero;den4(1)=zero
 do ir=2,mesh
   call ctrap(ir,den1(1:ir),step,den3(ir))
   call ctrap(ir,den2(1:ir),step,den4(ir))
 end do

 do ir=1,mesh
   rvhn(ir)=den3(ir)+rad(ir)*(den4(mesh)-den4(ir))
 end do

 ABI_FREE(den1)
 ABI_FREE(den2)
 ABI_FREE(den3)
 ABI_FREE(den4)

 vh_tnzc(2:mesh)=rvhn(2:mesh)/rad(2:mesh)
 yp2=(vh_tnzc(3)-vh_tnzc(2))/(rad(3)-rad(2))
 yp3=(vh_tnzc(4)-vh_tnzc(3))/(rad(4)-rad(3))
 yp1=yp2+(yp2-yp3)*rad(2)/(rad(3)-rad(2))
 vh_tnzc(1)=vh_tnzc(2)-(yp1+yp2)*rad(2)

 ABI_FREE(nzc)
 ABI_FREE(rvhn)

end subroutine vhtnzc
!!***

!!****f* m_psp6/psp6cc_drh
!! NAME
!! psp6cc_drh
!!
!! FUNCTION
!! Compute the core charge density, for use in the XC core
!! correction, following the function definition valid
!! for the format 6 of pseudopotentials.
!! Version modified by DHamann, with consistent treatment
!! of the derivatives in this routine and the remaining of the code.
!!
!! INPUTS
!!  mmax=maximum number of points in real space grid in the psp file
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  rchrg=cut-off radius for the core density
!!
!! OUTPUT
!!  xccc1d(n1xccc,6)= 1D core charge function and its five first derivatives
!!
!! PARENTS
!!      m_psp6
!!
!! CHILDREN
!!      cc_derivatives
!!
!! NOTES
!! Test version by DRH - requires very smooth model core charge
!!
!! SOURCE

subroutine psp6cc_drh(mmax,n1xccc,rchrg,xccc1d)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mmax,n1xccc
 real(dp),intent(in) :: rchrg
!arrays
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i

!Local variables-------------------------------
!scalars
 integer :: irad
 character(len=500) :: errmsg
!arrays
 real(dp),allocatable :: ff(:),ff1(:),ff2(:),rad(:)

!**********************************************************************

 ABI_MALLOC(ff,(mmax))
 ABI_MALLOC(ff1,(mmax))
 ABI_MALLOC(ff2,(mmax))
 ABI_MALLOC(rad,(mmax))

!
!read from pp file the model core charge (ff) and first (ff1) and
!second (ff2) derivative on logarithmic mesh mmax; rad is the radial grid
!the input functions contain the 4pi factor, it must be rescaled.

 !write(std_out,'(a,2i6)') 'drh:psp6cc_drh - mmax,n1xccc',mmax,n1xccc
 do irad=1,mmax
   read(tmp_unit,*,err=10,iomsg=errmsg) rad(irad),ff(irad),ff1(irad),ff2(irad)
   ff(irad)=ff(irad)/4.d0/pi
   ff1(irad)=ff1(irad)/4.d0/pi
   ff2(irad)=ff2(irad)/4.d0/pi
 end do
 rad(1)=0.d0

 call cc_derivatives(rad,ff,ff1,ff2,mmax,n1xccc,rchrg,xccc1d)

 ABI_FREE(ff)
 ABI_FREE(ff1)
 ABI_FREE(ff2)
 ABI_FREE(rad)

 return

 ! Handle IO error
 10 continue
 ABI_ERROR(errmsg)

end subroutine psp6cc_drh
!!***

end module m_psp6
!!***
