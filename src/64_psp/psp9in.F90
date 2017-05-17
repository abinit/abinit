!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp9in
!! NAME
!! psp9in
!!
!! FUNCTION
!! Initialize pspcod=9 (Pseudopotentials from the PSML XML format):
!! continue to read the corresponding file, then compute the
!! local and non-local potentials.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (JJ, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!  mpssoang= 2*maximum angular momentum for nonlocal pseudopotentials - 1
!!  mqgrid=dimension of q (or G) grid for arrays.
!!  pspso=spin-orbit characteristics, govern the content of ffspl and ekb
!!   if =0 : this input requires NO spin-orbit characteristics of the psp
!!   if =2 : this input requires HGH or psp8 characteristics of the psp
!!   if =3 : this input requires HFN characteristics of the psp
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  optnlxccc=option for nl XC core correction (input variable)
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
!!  nproj(mpssoang)=number of projection functions for each angular momentum
!!  qchrg is not used, and could be suppressed later
!!  vlspl(mqgrid,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  xcccrc=XC core correction cutoff radius (bohr)
!!  xccc1d(n1xccc,6)=1D core charge function and five derivatives, from psp file
!!
!! PARENTS
!!      pspatm
!!
!! CHILDREN
!!      cc_derivatives,ps_destroy,psml_reader,psp8lo,psp8nl,spline,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine psp9in(filpsp,ekb,epsatm,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&
&                  mmax,mpsang,mpssoang,mqgrid,nproj,n1xccc,pspso,qchrg,qgrid,&
&                  useylm,vlspl,xcccrc,xccc1d,zion,znucl)


 use defs_basis
 use m_splines
 use m_errors
 use m_profiling_abi
#if defined HAVE_PSML
 use m_psml
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp9in'
 use interfaces_14_hidewrite
 use interfaces_64_psp, except_this_one => psp9in
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lloc,lmax,lmnmax,lnmax,mpsang,mpssoang,mqgrid,n1xccc
 integer,intent(in) :: pspso,useylm
 integer,intent(out) :: mmax
 real(dp),intent(in) :: zion,znucl
 real(dp),intent(out) :: epsatm,qchrg,xcccrc
 character(len=fnlen),intent(in) :: filpsp
!arrays
 integer,intent(out) :: indlmn(6,lmnmax),nproj(mpssoang)
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(out) :: ekb(lnmax),ffspl(mqgrid,2,lnmax),vlspl(mqgrid,2)
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i

!Local variables-------------------------------
!no_abirules
 integer :: ipsang,ir,jj
 integer :: iproj
 integer :: il,index,irelt,ii
 integer :: iln, kk, ll, nn, mm, nso 
 real(dp) :: amesh,fchrg,zval
 real(dp) :: rmax,yp1,ypn,z,chgvps
 real(dp) :: eps
 character(len=500) :: message

 integer, allocatable :: idx_sr(:)
 integer, allocatable :: idx_so(:)
 real(dp),allocatable :: funrad(:),rad(:)
 real(dp),allocatable :: vloc(:),vpspll(:,:)
 real(dp),allocatable ::       ff1(:),ff2(:)
 real(dp),allocatable :: work_spl(:)
 parameter( eps = 1.0d-4 )
!----------------------------------------------------------------------------
! xcccrc           : XC core correction cutoff radius (bohr)
!                    It is defined as the radius where the pseudo-core
!                    charge density becomes zero
!                    (here we have set up a tolerance of 1.d-12).

#if defined HAVE_PSML
 type(ps_t) :: psxml
#endif

! ***************************************************************************

 if(.false.)write(std_out,*)filpsp ! Just to keep filpsp when HAVE_PSML is false
 if(.false.)write(std_out,*)lloc   ! Just to keep lloc when HAVE_PSML is false
 if(.false.)write(std_out,*)lmax   ! Just to keep lmax when HAVE_PSML is false
 if(.false.)write(std_out,*)qgrid  ! Just to keep qgrid when HAVE_PSML is false
 if(.false.)write(std_out,*)useylm ! Just to keep useylm when HAVE_PSML is false
 if(.false.)write(std_out,*)zion   ! Just to keep zion when HAVE_PSML is false
 if(.false.)write(std_out,*)znucl  ! Just to keep znucl when HAVE_PSML is false

#if defined HAVE_PSML

 call ps_destroy(psxml)
 call psml_reader(filpsp,psxml,debug=.true.)

!Define the mesh used by the FHI pseudopotential code
!The atomic number of the element is read from the header of the XML file
 z      = ps_AtomicNumber(psxml)
 zval   = ps_Zpseudo(psxml)

!---

! Compute the valence charge of the reference configuration used to 
! generate the pseudopotential
 chgvps = 0.0_dp
 do il = 1, ps_NValenceShells(psxml)
   chgvps = chgvps + ps_ValenceShellOccupation(psxml, il)
   write(std_out,*)' psp9in : occupation = ', ps_ValenceShellOccupation(psxml, il)
 end do

!DEBUG
 write(std_out,*)' psp9in : valence charge of the reference configuration '
 write(std_out,*)' psp9in : chgvps = ', chgvps
 write(std_out,*)' psp9in : nominal valence charge '
 write(std_out,*)' psp9in : zval   = ', zval
!ENDDEBUG

!Feb 2015: shifted to Hamann grid for convenience - libpsml interpolates anyway

! The following lines are taken from the oncvpsp.f90 subroutine of the oncvpsp
! code implemented by D. Hamann
! The atomic number of the element is read from the header of the XML file
! Logarithmic grid defined by Hamann in oncvpsp code
! z    = psxml%header%z
! amesh = 1.012d0
! al    = dlog(amesh)
! rr1   = .0005d0/z
! mmax  = dlog(45.0d0 /rr1)/al
!
! ABI_ALLOCATE( rad,(mmax) )
!
! do ir = 1, mmax
!   rad(ir) = rr1 * dexp(al*(ir-1))
! end do
!Determine the maximum number of points in the grid ---
 rmax = 6.0_dp
 amesh    = 0.01
 mmax     = int(rmax/amesh)
 if(mod(mmax,2) .eq. 0) mmax = mmax + 1
!DEBUG
 write(std_out,*)' psp9in : parameters to define the points of the grid '
 write(std_out,*)' psp9in : amesh = ', amesh
 write(std_out,*)' psp9in : rmax = ', rmax
 write(std_out,*)' psp9in : mmax = ', mmax
!ENDDEBUG
!---

!Feb 2015: shifted to Hamann grid for convenience - libpsml interpolates anyway
 ABI_ALLOCATE( rad,(mmax))
 rad(1) = zero
 do ir = 2, mmax
   rad(ir) = rad(ir-1) + amesh
 end do
!! DEBUG
! do ir = 2, mmax
!   write(std_out,'(i5,f20.12)')ir, rad(ir)
! end do
!! ENDDEBUG
!---
 
! TODO: should be simple to average these and get difference for SREL+SOC,
! but also the Ekb etc...
 if (ps_Number_Of_Projectors(psxml,SET_LJ) > 0) then
   message = 'For the moment LJ format projectors are not supported; SREL + SO is the internal abinit format'
   MSG_BUG(message)
 end if

 if (ps_Number_Of_Projectors(psxml,SET_UP) > 0 .or. ps_Number_Of_Projectors(psxml,SET_DOWN) > 0) then
   write (message,'(3a)') 'For the moment separate spin up and down format projectors are not supported;',ch10,&
&   ' spin average is the internal abinit format'
   MSG_BUG(message)
 end if

 if(pspso==2) then
   nso=2
 else
   nso=1
 end if

!    Find the number of projectors per angular momentum shell
 nproj(:)=0
 if (ps_Number_Of_Projectors(psxml,SET_NONREL) > 0) then
   idx_sr = ps_Projector_Indexes(psxml,SET_NONREL)
   do iproj = 1, ps_Number_Of_Projectors(psxml,SET_NONREL)
     il = ps_Projector_L(psxml, idx_sr(iproj))
     nproj(il+1) = nproj(il+1) + 1
   end do
 else if (ps_Number_Of_Projectors(psxml,SET_SREL) > 0) then
   idx_sr = ps_Projector_Indexes(psxml,SET_SREL)
   do iproj = 1, ps_Number_Of_Projectors(psxml,SET_SREL)
     il = ps_Projector_L(psxml, idx_sr(iproj))
     nproj(il+1) = nproj(il+1) + 1
   end do
 else ! this should not happen
   MSG_BUG('Your psml potential should have either scalar- or non- relativistic projectors')
 end if

 irelt = 0
 if (nso == 2) then
   idx_so = ps_Projector_Indexes(psxml,SET_SO)
   do iproj = 1, ps_Number_Of_Projectors(psxml,SET_SO)
     il = ps_Projector_L(psxml, idx_so(iproj))
     nproj(il+lmax+2) = nproj(il+lmax+2) + 1
     irelt = 1
   end do
 end if

! Determine whether the atomic calculation to generate the pseudopotential
! is relativistic or not

!DEBUG
 write(std_out,*)' psp9in : pseudopotential generation relativity ', ps_Relativity(psxml) 
 write(std_out,*)' psp9in : SOC pseudopotential? (1=yes, 0 =no) '
 write(std_out,*)' psp9in : irelt = ', irelt
 write(ab_out,*)' psp9in : irelt = ', irelt
!ENDDEBUG


!Take care of the non-linear core corrections

 if (ps_HasCoreCorrections(psxml)) then
!    In Abinit, at least for the Troullier-Martins pseudopotential,
!    the pseudocore charge density and its derivatives (xccc1d)
!    are introduced in a linear grid.
!    This grid is normalized, so the radial coordinates run between
!    from 0 and 1 (from 0 to xcccrc, where xcccrc is the radius
!    where the pseudo-core becomes zero).

!    Allocate some of the arrays
   ABI_ALLOCATE( funrad,(mmax))
   ABI_ALLOCATE( ff1,(mmax))
   ABI_ALLOCATE( ff2,(mmax))

!    Store the value of the pseudo-core charge.
   do jj=1,mmax
     funrad(jj) = ps_CoreCharge_Value(psxml,rad(jj))
   end do
   funrad = funrad / four / pi

!    determine xcccrc where the pseudocore becomes 0
   do jj=mmax,1,-1
     if (funrad(jj) > 1.0d-12) then
       xcccrc=rad(jj)
       exit
     end if
   end do
! for the moment use the fixed value common in oncvpsp outputs
!   if (xcccrc > five) then
!     write (*,*) 'xcccrc jj funrad(jj) ', xcccrc, jj, funrad(jj)
!     MSG_WARNING('The core charge radius has been truncated to 5 bohrs')
!   end if
!   xcccrc = five

! first derivatives at extreme points

! OLD low order formulae
!   yp1       = ( funrad(2) - funrad(1) )/( rad(2) - rad(1) )
!
!!    Find the first derivative of the pseudo-core charge
!!    in the grid.
!   ff1(1)=yp1
!   do jj=2,mmax-1
!     ff1(jj)=(funrad(jj+1)-funrad(jj))/(rad(jj+1)-rad(jj))
!   end do
   
   do jj=4,mmax-3
     ff1(jj) = -1._dp/60._dp*funrad(jj-3) +3._dp/20._dp*funrad(jj-2) -3._dp/4._dp*funrad(jj-1)  &
&     +1._dp/60._dp*funrad(jj+3) -3._dp/20._dp*funrad(jj+2) +3._dp/4._dp*funrad(jj+1)
   end do

   do jj=1,3
     ff1(jj)=  -49._dp/20._dp*funrad(jj) +6._dp*funrad(jj+1) -15._dp/2._dp*funrad(jj+2) &
&     +20._dp/3._dp*funrad(jj+3) -15._dp/4._dp*funrad(jj+4) +6._dp/5._dp*funrad(jj+5) -1._dp/6._dp*funrad(jj+6)
   end do

   ff1(mmax-2:mmax)=zero 

!    Find the second derivative of the pseudo-core charge
!    in the grid.
!    Be careful, this is very noisy at r->0
   call spline( rad, funrad, mmax, yp1, ypn, ff2 )

!    call cc_derivatives to get 3rd 4th and 5th derivatives,
!    and everything splined onto regular grid [0:xcccrc]
!    in xccc1d
   call cc_derivatives(rad,funrad,ff1,ff2,mmax,n1xccc,xcccrc,xccc1d)

!    write(std_out,*) '# psp9in NLCC data ', n1xccc
!    do ii = 1, n1xccc
!    write(std_out,'(7e20.8)')xcccrc*(ii-1.d0)/(n1xccc-1.d0),xccc1d(ii,1),&
! &         xccc1d(ii,2),xccc1d(ii,3),xccc1d(ii,4),xccc1d(ii,5),xccc1d(ii,6)
!    enddo
!!    ENDDEBUG

   ABI_DEALLOCATE( funrad )
   ABI_DEALLOCATE( ff1 )
   ABI_DEALLOCATE( ff2 )
 else
   write(message, '(a)' ) '  No XC core correction.'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   xcccrc = 0.0_dp ; fchrg = 0.0_dp ; qchrg = 0.0_dp
 end if

!Define the local component of the pseudopotential
 ABI_ALLOCATE(vloc,(mmax))

 vloc = zero
 do ir = 1, mmax
   vloc(ir) = ps_LocalPotential_Value(psxml, rad(ir))
 end do 

!--------------------------------------------------------------------
!Carry out calculations for local (lloc) pseudopotential.
!Obtain Fourier transform (1-d sine transform)
!to get q^2 V(q).

 call psp8lo(amesh,epsatm,mmax,mqgrid,qgrid,vlspl(:,1),rad,vloc,yp1,ypn,zion)

!Fit spline to q^2 V(q) (Numerical Recipes subroutine)
 ABI_ALLOCATE(work_spl,(mqgrid))
 call spline (qgrid,vlspl(:,1),mqgrid,yp1,ypn,work_spl)
 vlspl(:,2)=work_spl(:)
 ABI_DEALLOCATE(work_spl)

!--------------------------------------------------------------------
!Take care of non-local part


!Initialize array indlmn array giving l,m,n,lm,ln,s for i=lmn
 index=0;iln=0;indlmn(:,:)=0
 do nn=1,nso
   do ipsang=1+(nn-1)*(lmax+1), nn*lmax+1
     ll=ipsang-(nn-1)*lmax-1
     if (nproj(ipsang)>0) then
       do kk=1,nproj(ipsang)
         iln=iln+1
         do mm=1,2*ll*useylm+1
           index=index+1
           indlmn(1,index)=ll                      ! l angular momentum channel
           indlmn(2,index)=mm-ll*useylm-1          ! hash of position in m
           indlmn(3,index)=kk                      ! index of projector
           indlmn(4,index)=ll*ll+(1-useylm)*ll+mm  ! hash of position in l(l+1) array
           indlmn(5,index)=iln                     ! absolute index of l, n disregarding m values
           indlmn(6,index)=nn                      ! spin orbit index!!! NOT the n shell index
         end do
       end do
     end if
   end do
 end do

!Zero out all Kleinman-Bylander energies to initialize
 do ii = 1, lmnmax ! loop over all possible projectors
   if (indlmn(6,ii) == 1) then
     ekb (indlmn(5,ii)) = ps_Projector_Ekb(psxml, idx_sr(indlmn(5,ii)))
   else if (indlmn(6,ii) == 2) then
     ekb (indlmn(5,ii)) = ps_Projector_Ekb(psxml, idx_so(indlmn(5,ii)))
   end if
 end do


 ABI_ALLOCATE(vpspll,(mmax,lnmax))
 vpspll = zero
 do ii = 1, lmnmax
   if (indlmn(6,ii) == 1) then
     do ir = 1, mmax
       vpspll(ir, indlmn(5,ii)) = ps_Projector_Value(psxml, idx_sr(indlmn(5,ii)), rad(ir))
     end do 
   else if (indlmn(6,ii) == 2) then
     do ir = 1, mmax
       vpspll(ir, indlmn(5,ii)) = ps_Projector_Value(psxml, idx_so(indlmn(5,ii)), rad(ir))
     end do 
   end if
 end do 

!Allow for option of no nonlocal corrections (lloc=lmax=0)
 if (lloc==0.and.lmax==0) then
   write(message, '(a,f5.1)' ) ' Note: local psp for atom with Z=',znucl
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 else
   call psp8nl(amesh,ffspl,indlmn,lmax,lmnmax,lnmax,mmax,&
&   mqgrid,qgrid,rad,vpspll)
 end if

!!DEBUG
!write(std_out,*)'# lmax     = ', lmax
!write(std_out,*)'# lhigh    = ', lhigh
!write(std_out,*)'# lmnmax   = ', lmnmax
!write(std_out,*)'# lnmax    = ', lnmax
!write(std_out,*)'# mpsang   = ', mpsang
!write(std_out,*)'# mpssoang = ', mpssoang
!write(std_out,*)'# nproj    = ', nproj(:)
!write(std_out,*)'# pspso    = ', pspso
!write(std_out,*)'# filpsp   = ', filpsp
!!ENDDEBUG



 ABI_DEALLOCATE( vpspll   )
 ABI_DEALLOCATE( rad      )
 ABI_DEALLOCATE( vloc     )
 if (allocated(idx_sr)) then
   ABI_FREE_NOCOUNT(idx_sr)
 end if
 if (allocated(idx_so)) then
   ABI_FREE_NOCOUNT(idx_so)
 end if

 call ps_destroy(psxml)


#else
!Initialize some arguments, for portability at compile time
 indlmn=0 ; mmax=0 ; nproj=0
 ekb=zero ; epsatm=zero ; ffspl=zero ; qchrg=zero ; vlspl=zero ; xcccrc=zero ; xccc1d=zero
#endif

end subroutine psp9in
!!***
