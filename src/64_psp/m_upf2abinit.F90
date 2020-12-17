!!****m* ABINIT/m_upf2abinit
!! NAME
!!  m_upf2abinit
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2020 ABINIT group (MJV)
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

module m_upf2abinit

 use defs_basis
 use m_splines
 use m_abicore
 use m_errors
 use m_atomdata
 use pseudo_pwscf ! pwscf module with all data explicit!
 use m_read_upf_pwscf, only : read_pseudo

 use defs_datatypes,  only : pseudopotential_type
 use m_io_tools,      only : open_file
 use m_numeric_tools, only : smooth, nderiv, ctrap
 use m_pspheads,      only : upfxc2abi
 use m_paw_numeric,   only : jbessel => paw_jbessel
 use m_psptk,         only : cc_derivatives

 implicit none

 private
!!***

 public :: upf2abinit
!!***

contains
!!***

!!****f* ABINIT/upf2abinit
!! NAME
!! upf2abinit
!!
!! FUNCTION
!!  This routine wraps a call to a PWSCF module, which reads in
!!  a UPF (PWSCF / Espresso) format pseudopotential, then transfers
!!  data to abinit internal variables.
!!  "UPF PWSCF format" (pspcod=11)
!!
!! INPUTS
!!  filpsp = name of file with UPF data
!!  psps = sturcture with global dimension data for pseudopotentials, header info ...
!!
!! OUTPUT
!!  pspxc = index of xc functional for this pseudo
!!  lmax_ = maximal angular momentum
!!  lloc = local component chosen for pseudopotential
!!  mmax = maximum number of points in real space radial grid
!!  znucl = charge of species nucleus
!!  zion = valence charge
!!  epsatm = integral of local potential - coulomb potential of zion
!!  xcccrc = radius for non linear core correction
!!  ekb(dimekb)= Kleinman Bylander energies, see pspatm.F90
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$(hartree)
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  vlspl(mqgrid_vl,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  ffspl(mqgrid_ff,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector; if any, spin-orbit components begin at l=mpsang+1
!!  nproj= number of projectors for each channel
!!  xccc1d(n1xccc*(1-usepaw),6)=1D core charge function and five derivatives,
!!                              from psp file (used in NC only)
!!
!! PARENTS
!!      m_pspini
!!
!! CHILDREN
!!      ctrap
!!
!! SOURCE

subroutine upf2abinit (filpsp, znucl, zion, pspxc, lmax_, lloc, mmax, &
&  psps, epsatm, xcccrc, indlmn, ekb, ffspl, nproj_l, vlspl, xccc1d)

  implicit none

!Arguments -------------------------------

  character(len=fnlen), intent(in) :: filpsp
  type(pseudopotential_type),intent(in) :: psps
  ! used contents:
  !   psps%lmnmax
  !   psps%mqgrid_ff
  !   psps%mqgrid_vl
  !   psps%dimekb
  !   psps%n1xccc
  !   psps%qgrid_ff
  !   psps%qgrid_vl

  integer, intent(out) :: pspxc, lmax_, lloc, mmax
  real(dp), intent(out) :: znucl, zion
  real(dp), intent(out) :: epsatm, xcccrc
  !arrays
  integer, intent(out)  :: indlmn(6,psps%lmnmax)
  integer, intent(out)  :: nproj_l(psps%mpssoang)
  real(dp), intent(inout) :: ekb(psps%dimekb) !vz_i
  real(dp), intent(inout) :: ffspl(psps%mqgrid_ff,2,psps%lnmax) !vz_i
  real(dp), intent(out) :: vlspl(psps%mqgrid_vl,2)
  real(dp), intent(inout) :: xccc1d(psps%n1xccc,6) !vz_i

!Local variables -------------------------
  integer :: ir, iproj, ll, iunit
  real(dp) :: yp1, ypn
  character(len=500) :: msg
  type(atomdata_t) :: atom

  logical, allocatable :: found_l(:)
  real(dp), allocatable :: work_space(:),work_spl(:)
  real(dp), allocatable :: ff(:), ff1(:), ff2(:), rad_cc(:), proj(:,:)

  ! ######### in module pseudo: ############
  !
  !  only npsx = 1 is used here
  !  grids are allocated for much larger fixed length (ndm=2000)
  !  number of species (6) and projectors (8) as well...
  !
  !  psd(npsx) = specied string
  !  pseudotype = uspp / nc string
  !  dft(npsx) = exchange correlation string (20 chars)
  !  lmax(npsx) = maximum l channel
  !  mesh(npsx) = number of points for local pot
  !  nbeta(npsx) = number of projectors (beta functions for uspp)
  !  nlcc(npsx) = flag for presence of NL core correction
  !  zp(npsx) = valence ionic charge
  !  r(ndm,npsx) = radial mesh
  !  rab(ndm,npsx) = dr / di for radial mesh
  !  rho_atc(ndm,npsx) = NLCC pseudocharge density
  !  vloc0(ndm,npsx) = local pseudopotential
  !  betar(ndm, nbrx, npsx) = projector functions in real space mesh
  !  lll(nbrx,npsx) = angular momentum channel for each projector
  !  ikk2(nbrx,npsx) = maximum index for each projector function
  !  dion(nbrx,nbrx,npsx) = dij or Kleinman Bylander energies
  !
  !  ########  end description of pseudo module contents ##########

! *********************************************************************

!call pwscf routine for reading in UPF
 if (open_file (filpsp,msg,newunit=iunit,status='old',form='formatted') /= 0) then
   MSG_ERROR(msg)
 end if

!read in psp data to static data in pseudo module, for ipsx == 1
 call read_pseudo(1,iunit)
 close (iunit)

!convert to Ha units
 vloc0  = half * vloc0
!betar = half * betar ! ???
 dion   = half * dion

!if upf file is a USPP one, stop
 if (pseudotype == 'US') then
   MSG_ERROR('upf2abinit: USPP UPF files not supported')
 end if

!copy over to abinit internal arrays and vars
 call upfxc2abi(dft(1), pspxc)
 lmax_ = lmax(1)

!Check if the local component is one of the angular momentum channels
!effectively if one of the ll is absent from the NL projectors
 ABI_ALLOCATE(found_l,(0:lmax_))
 found_l = .true.
 do ll = 0, lmax_
   if (any(lll(1:nbeta(1),1) == ll)) then
     found_l(ll) = .false.
   end if
 end do
 if (count(found_l) /= 1) then
   lloc = -1
 else
   do ll = 0, lmax_
     if (found_l(ll)) then
       lloc = ll
       exit
     end if
   end do
 end if
 ABI_DEALLOCATE(found_l)
!FIXME: do something about lloc == -1

 call atomdata_from_symbol(atom,psd(1))
 znucl = atom%znucl
 zion = zp(1)
 mmax = mesh(1)

 call psp11lo(rab(1:mmax,1),epsatm,mmax,psps%mqgrid_vl,psps%qgrid_vl,&
& vlspl(:,1),r(1:mmax,1),vloc0(1:mmax,1),yp1,ypn,zion)

!Fit spline to q^2 V(q) (Numerical Recipes subroutine)
 ABI_ALLOCATE(work_space,(psps%mqgrid_vl))
 ABI_ALLOCATE(work_spl,(psps%mqgrid_vl))
 call spline (psps%qgrid_vl,vlspl(:,1),psps%mqgrid_vl,yp1,ypn,work_spl)
 vlspl(:,2)=work_spl(:)
 ABI_DEALLOCATE(work_space)
 ABI_DEALLOCATE(work_spl)

!this has to do the FT of the projectors to reciprocal space
! allocate proj to avoid temporary copy.
 ABI_ALLOCATE(proj, (mmax,1:nbeta(1)))
 proj = betar(1:mmax,1:nbeta(1),1)

 call psp11nl(ffspl, indlmn, mmax, psps%lnmax, psps%lmnmax, psps%mqgrid_ff, &
& nbeta(1), proj, lll(1:nbeta(1),1), ikk2(1:nbeta(1),1), &
& psps%qgrid_ff, r(1:mmax,1), rab(1:mmax,1), psps%useylm)

 ABI_FREE(proj)

 nproj_l = 0
 do iproj = 1, nbeta(1)
   ll = lll(iproj,1)
   nproj_l(ll+1) = nproj_l(ll+1) + 1
 end do

!shape = dimekb  vs. shape = n_proj
 do ll = 1, nbeta(1)
   ekb(ll) = dion(ll,ll,1)
 end do

 xcccrc = zero
 xccc1d = zero
!if we find a core density, do something about it
!rho_atc contains the nlcc density
!rho_at contains the total density
 if (nlcc(1)) then
   ABI_ALLOCATE(ff,(mmax))
   ABI_ALLOCATE(ff1,(mmax))
   ABI_ALLOCATE(ff2,(mmax))
   ff(1:mmax) = rho_atc(1:mmax,1) ! model core charge without derivative factor

   ff1 = zero
   call nderiv(one,ff,ff1,mmax,1) ! first derivative
   ff1(1:mmax) = ff1(1:mmax) / rab(1:mmax,1)
   call smooth(ff1, mmax, 15) ! run 15 iterations of smoothing

   ff2 = zero
   call nderiv(one, ff1, ff2, mmax, 1) ! second derivative
   ff2(1:mmax) = ff2(1:mmax) / rab(1:mmax,1)
   call smooth(ff2, mmax, 15) ! run 10 iterations of smoothing?

!  determine a good rchrg = xcccrc
   do ir = mmax, 1, -1
     if (abs(ff(ir)) > 1.e-6) then
       xcccrc = r(ir,1)
       exit
     end if
   end do
   ABI_ALLOCATE(rad_cc,(mmax))
   rad_cc = r(1:mmax,1)
   rad_cc(1) = zero ! force this so that the core charge covers whole spline interval.
   call cc_derivatives(rad_cc,ff,ff1,ff2,mmax,psps%n1xccc,xcccrc,xccc1d)
   ABI_DEALLOCATE(rad_cc)
   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(ff1)
   ABI_DEALLOCATE(ff2)

 end if !if nlcc present

end subroutine upf2abinit
!!***

!!****f* m_upf2abinit/psp11nl
!! NAME
!! psp11nl
!!
!! FUNCTION
!! Fourier transform the real space UPF projector functions to reciprocal space
!!
!! INPUTS
!!  lmax=maximum ang momentum for which nonlocal form factor is desired.
!!   Usually lmax=1, sometimes = 0 (e.g. for oxygen); lmax <= 2 allowed.
!!  mmax=number of radial grid points for atomic grid
!!  lnmax= maximum index for all l channel projectors, dimension of ffspl
!!  lmnmax= maximum index for all projectors, dimension of indlmn
!!  mqgrid=number of grid points for q grid
!!  n_proj = total number of NL projectors read in
!!  proj = projector data times r, on a real space grid
!!  proj_l = ang mom channel for each projector
!!  proj_np = max number of points used for each projector
!!  qgrid(mqgrid)=values at which form factors are returned
!!  r(mmax)=radial grid values
!!  drdi=derivative of grid point wrt index
!!  useylm = input to use m dependency of NL part, or only Legendre polynomials
!!
!! OUTPUT
!!  ffspl(mqgrid,2,mpsang)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum
!!  indlmn = indexing of each projector, for n, l, m, s, ln, lmn (see pspatm.F90)
!!
!! NOTES
!!
!! PARENTS
!!      m_upf2abinit
!!
!! CHILDREN
!!      ctrap
!!
!! SOURCE

subroutine psp11nl(ffspl,indlmn,mmax,lnmax,lmnmax,mqgrid,n_proj,&
&                  proj, proj_l, proj_np, qgrid, r, drdi, useylm)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mmax, lnmax, lmnmax, mqgrid, useylm, n_proj
!arrays
 integer, intent(in) :: proj_l(n_proj)
 integer, intent(in) :: proj_np(n_proj)
 integer, intent(out) :: indlmn(6,lmnmax)
 real(dp),intent(in) :: r(mmax)
 real(dp),intent(in) :: drdi(mmax)
 real(dp),intent(in) :: proj(mmax,n_proj)
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(inout) :: ffspl(mqgrid,2,lnmax) !vz_i

!Local variables-------------------------------
!scalars
 integer :: iproj, np, ll, llold, ipsang, i_indlmn
 integer :: iproj_1l, ir, iq, mm
 integer :: bessorder
 real(dp) :: res, arg, besfact, dummy, dummy2
 real(dp), allocatable :: work(:)
 character(len=500) :: message

!*************************************************************************
 bessorder = 0 ! never calculate derivatives of bessel functions

 ffspl = zero
 indlmn = 0
 i_indlmn = 0
 llold = -1
 iproj_1l = 1
!big loop over all projectors
 do iproj = 1, n_proj

   if (iproj > lmnmax) then
     write(message,'(a,2i0)') ' Too many projectors found. n_proj, lmnmax =  ',n_proj, lmnmax
     MSG_ERROR(message)
   end if

   np = proj_np(iproj)
   ABI_ALLOCATE(work,(np))
   ll = proj_l(iproj)
   if (ll < llold) then
     message = 'psp11nl : Error: UPF projectors are not in order of increasing ll'
     MSG_ERROR(message)
   else if (ll == llold) then
     iproj_1l = iproj_1l + 1
   else
     iproj_1l = 1
     llold = ll
   end if
!  determine indlmn for this projector (keep in UPF order and enforce that they are in
!  increasing ll)
   do mm = 1, 2*ll*useylm+1
     i_indlmn = i_indlmn + 1
     indlmn(1,i_indlmn) = ll
     indlmn(2,i_indlmn) = mm-ll*useylm-1
     indlmn(3,i_indlmn) = iproj_1l
     indlmn(4,i_indlmn) = ll*ll+(1-useylm)*ll+mm
     indlmn(5,i_indlmn) = iproj
     indlmn(6,i_indlmn) = 1 !spin? FIXME: to get j for relativistic cases
   end do

!  FT projectors to reciprocal space q
   do iq = 1, mqgrid
     arg = two_pi*qgrid(iq)

!    FIXME: add semianalytic form for integral from 0 to first point
     do ir = 1, np
       call jbessel(besfact, dummy, dummy2, ll, bessorder, arg*r(ir))
!      besfact = sin(arg*r(ir))
       work(ir) = drdi(ir) * besfact * proj(ir, iproj) * r(ir) !* r(ir)
     end do
     call ctrap (np, work, one, res)

     ffspl(iq, 1, iproj) = res
   end do
   ABI_DEALLOCATE(work)
 end do  ! iproj

!add derivative of ffspl(:,1,:) for spline interpolation later
 ABI_ALLOCATE(work,(mqgrid))
 do ipsang = 1, lnmax
   call spline(qgrid,ffspl(:,1,ipsang),mqgrid,zero,zero,ffspl(:,2,ipsang))
 end do
 ABI_DEALLOCATE(work)

end subroutine psp11nl
!!***

!!****f* ABINIT/psp11lo
!! NAME
!! psp11lo
!!
!! FUNCTION
!! Compute sine transform to transform from V(r) to q^2 V(q).
!! Computes integrals on logarithmic grid using related uniform
!! grid in exponent and corrected trapezoidal integration.
!! Generalized from psp5lo for non log grids using dr / di
!!
!! INPUTS
!!  drdi=derivative of radial grid wrt index
!!  mmax=number of radial r grid points
!!  mqgrid=number of grid points in q from 0 to qmax.
!!  qgrid(mqgrid)=q grid values (bohr**-1).
!!  rad(mmax)=r grid values (bohr).
!!  vloc(mmax)=V(r) on radial grid.
!!  zion=nominal valence charge of atom.
!!
!! OUTPUT
!!  epsatm=$ 4\pi\int[r^2 (V(r)+\frac{Zv}{r}dr]$.
!!{{\\ \begin{equation}
!!  q2vq(mqgrid)
!!   =q^2 V(q)
!!   = -\frac{Zv}{\pi}
!!     + q^2 4\pi\int[(\frac{\sin(2\pi q r)}{2\pi q r})(r^2 V(r)+r Zv)dr].
!!\end{equation} }}
!!  yp1,ypn=derivative of q^2 V(q) wrt q at q=0 and q=qmax
!!   (needed for spline fitter).
!!
!! PARENTS
!!      m_upf2abinit
!!
!! CHILDREN
!!      ctrap
!!
!! SOURCE

subroutine psp11lo(drdi,epsatm,mmax,mqgrid,qgrid,q2vq,rad,vloc,yp1,ypn,zion)

 implicit none

!Arguments----------------------------------------------------------
!scalars
 integer,intent(in) :: mmax,mqgrid
 real(dp),intent(in) :: zion
 real(dp),intent(out) :: epsatm,yp1,ypn
!arrays
 real(dp),intent(in) :: drdi(mmax)
 real(dp),intent(in) :: qgrid(mqgrid),rad(mmax),vloc(mmax)
 real(dp),intent(out) :: q2vq(mqgrid)

!Local variables-------------------------------
!scalars
 integer :: iq,ir
 real(dp),parameter :: scale=10.0d0
 real(dp) :: arg,result_ctrap,test,ztor1
!arrays
 real(dp),allocatable :: work(:)

! *************************************************************************

 ABI_ALLOCATE(work,(mmax))

!Do q=0 separately (compute epsatm)
!Do integral from 0 to r1
 ztor1=(zion/2.0d0+rad(1)*vloc(1)/3.d0)*rad(1)**2

!Set up integrand for q=0: $ \int[r^2 (V(r)+\frac{Zv}{r}) dr]$
!with extra factor of drdi to convert to uniform grid
 do ir = 1, mmax
!  First handle tail region
   test=vloc(ir)+zion/rad(ir)
!  Ignore small contributions, or impose a cut-off in the case
!  the pseudopotential data are in single precision.
!  (it is indeed expected that vloc is very close to zero beyond 20,
!  so a value larger than 2.0d-8 is considered anomalous)
   if (abs(test)<1.0d-20 .or. (rad(ir)>20.0d0 .and. abs(test)>2.0d-8) ) then
     work(ir)=zero
   else
     work(ir)=rad(ir)*(rad(ir)*vloc(ir)+zion)
   end if
   work(ir)=work(ir)*drdi(ir)
 end do

!Do integral from r(1) to r(max)
 call ctrap(mmax,work,one,result_ctrap)
 epsatm=4.d0*pi*(result_ctrap+ztor1)

 q2vq(1)=-zion/pi

!Loop over q values
 do iq=2,mqgrid
   arg=2.d0*pi*qgrid(iq)
!  ztor1=$ -Zv/\pi + 2q \int_0^{r1}[\sin(2\pi q r)(rV(r)+Zv) dr]$
   ztor1=(vloc(1)*sin(arg*rad(1))/arg-(rad(1)*vloc(1)+zion)* &
&   cos(arg*rad(1)) )/pi

!  set up integrand
   do  ir=1,mmax
!    test=vloc(ir)+zion/rad(ir)
!    Ignore contributions within decade of machine precision (suppressed ...)
!    if ((scale+abs(test)).eq.scale) then
!    work(ir)=zero
!    else
     work(ir)=sin(arg*rad(ir))*(rad(ir)*vloc(ir)+zion)
!    end if
     work(ir)=work(ir)*drdi(ir)
   end do
!  do integral from r(1) to r(mmax)
   call ctrap(mmax,work,one,result_ctrap)

!  store q^2 v(q)
!  FIXME: I only see one factor q, not q^2, but the same is done in other pspXlo.F90
   q2vq(iq)=ztor1+2.d0*qgrid(iq)*result_ctrap

 end do

!Compute derivatives of q^2 v(q) at ends of interval
 yp1=0.0d0
!ypn=$ 2\int_0^\infty[(\sin(2\pi qmax r)+(2\pi qmax r)*\cos(2\pi qmax r)(r V(r)+Z) dr]$
!integral from 0 to r1
 arg=2.0d0*pi*qgrid(mqgrid)
 ztor1=zion*rad(1)*sin(arg*rad(1))
 ztor1=ztor1+ 3.d0*rad(1)*vloc(1)*cos(arg*rad(1))/arg + &
& (rad(1)**2-1.0d0/arg**2)*vloc(1)*sin(arg*rad(1))
!integral from r(mmax) to infinity is overkill; ignore
!set up integrand
 do ir=1,mmax
!  test=vloc(ir)+zion/rad(ir)
!  Ignore contributions within decade of machine precision (supressed ...)
!  if ((scale+abs(test)).eq.scale) then
!  work(ir)=0.0d0
!  else
   work(ir)=(sin(arg*rad(ir))+arg*rad(ir)*cos(arg*rad(ir))) * &
&   (rad(ir)*vloc(ir)+zion)
!  end if
   work(ir)=work(ir)*drdi(ir)
 end do
 call ctrap(mmax,work,one,result_ctrap)
 ypn=2.0d0 * (ztor1 + result_ctrap)

 ABI_DEALLOCATE(work)

end subroutine psp11lo
!!***

end module m_upf2abinit
!!***
