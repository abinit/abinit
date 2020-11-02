!!****m* ABINIT/m_pawpsp
!! NAME
!!  m_pawpsp
!!
!! FUNCTION
!!  Module to read PAW atomic data
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2020 ABINIT group (MT, FJ,TR, GJ, FB, FrD, AF, GMR, DRH)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

module m_pawpsp

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 use m_libpaw_libxc
#if defined LIBPAW_HAVE_FOX
 use fox_sax
#endif

 use m_libpaw_tools, only : libpaw_basename, libpaw_get_free_unit

 use m_pawang, only: pawang_type
 use m_pawtab, only: pawtab_type, wvlpaw_type, wvlpaw_allocate, wvlpaw_rholoc_free, &
&                    pawtab_free, wvlpaw_free, wvlpaw_rholoc_nullify, pawtab_bcast, &
&                    pawtab_set_flags, wvlpaw_allocate, wvlpaw_free, wvlpaw_rholoc_nullify, &
&                    wvlpaw_rholoc_free
 use m_pawxmlps, only: rdpawpsxml_core, paw_setup_t, paw_setuploc, paw_setup_free
 use m_pawrad, only: pawrad_type, pawrad_init, pawrad_free, pawrad_copy, &
&      pawrad_bcast, pawrad_ifromr, simp_gen, nderiv_gen, bound_deriv, pawrad_deducer0, poisson
 use m_paw_numeric, only: paw_splint, paw_spline, paw_smooth, paw_jbessel_4spline
 use m_paw_atom, only: atompaw_shapebes, atompaw_vhnzc, atompaw_shpfun, &
&                     atompaw_dij0, atompaw_kij
 use m_pawxc, only: pawxc, pawxcm, pawxc_get_usekden
 use m_paw_gaussfit, only: gaussfit_projector

 implicit none

 private

 public:: pawpsp_calc_d5         !calculate up to the 5th derivative
 public:: pawpsp_main            !main routine to read psp
 public:: pawpsp_nl              !make paw projector form factors f_l(q)
 public:: pawpsp_read            !read psp from file
 public:: pawpsp_read_header     !read header of psp file
 public:: pawpsp_read_corewf     !read core wavefunction
 public:: pawpsp_read_header_2   !reads pspversion, basis_size and lmn_size
 public:: pawpsp_rw_atompaw      !read and writes ATOMPAW psp with gaussian |p>
 public:: pawpsp_wvl             !wavelet and icoulomb>0 related operations
 public:: pawpsp_wvl_calc        !wavelet related operations
 public:: pawpsp_7in             !reads non-XML atomic data
 public:: pawpsp_17in            !reads XML atomic data
 public:: pawpsp_calc            !calculates atomic quantities from psp info
 public:: pawpsp_read_header_xml !read header of psp file for XML
 public:: pawpsp_read_pawheader  !read header variables from XML objects
 public:: pawpsp_bcast           ! broadcast PAW psp data
 public:: pawpsp_cg              !compute sin FFT transform of a density
 public:: pawpsp_lo              !compute sin FFT transform of local potential

! Private procedures
 private:: pawpsp_wvl_sin2gauss  !convert sin/cos to gaussians
!!***

!-------------------------------------------------------------------------

!!****t* m_pawpsp/pawpsp_header_type
!! NAME
!! pawpsp_header_type
!!
!! FUNCTION
!! For PAW, header related data
!!
!! SOURCE

 type, public :: pawpsp_header_type

!Integer scalars
  integer :: basis_size    ! Number of elements of the wf basis ((l,n) quantum numbers)
  integer :: l_size        ! Maximum value of l+1 leading to a non zero Gaunt coefficient
  integer :: lmn_size      ! Number of elements of the paw basis
  integer :: mesh_size     ! Dimension of (main) radial mesh
  integer :: pawver        ! Version number of paw psp format
  integer :: shape_type    ! Type of shape function
  real(dp) :: rpaw         ! Radius for paw spheres
  real(dp) :: rshp         ! Cut-off radius of shape function

 end type pawpsp_header_type
!!***

CONTAINS
!===========================================================
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_nl
!! NAME
!! pawpsp_nl
!!
!! FUNCTION
!! Make paw projector form factors f_l(q) for each l
!!
!! INPUTS
!!  indlmn(6,lmnmax)= array giving l,m,n,lm,ln,s for i=lmn
!!  lmnmax=max number of (l,m,n) components
!!  lnmax=max number of (l,n) components
!!  mqgrid=number of grid points for q grid
!!  qgrid(mqgrid)=values at which form factors are returned
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!  wfll(:,lnmax)=paw projector on radial grid
!!
!! OUTPUT
!!  ffspl(mqgrid,2,lnmax)= form factor f_l(q) and second derivative
!!
!! NOTES
!!  u_l(r) is the paw projector (input as wfll);
!!  j_l(q) is a spherical Bessel function;
!!  f_l(q) = $ \int_0^{rmax}[j_l(2\pi q r) u_l(r)  r dr]$
!!
!! PARENTS
!!      m_paw_init,m_pawpsp
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_nl(ffspl,indlmn,lmnmax,lnmax,mqgrid,qgrid,radmesh,wfll)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmnmax,lnmax,mqgrid
 type(pawrad_type),intent(in) :: radmesh
!arrays
 integer,intent(in) :: indlmn(6,lmnmax)
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(in) ::  wfll(:,:)
 real(dp),intent(inout) :: ffspl(mqgrid,2,lnmax)

!Local variables-------------------------------
!scalars
 integer :: ilmn,iln,iln0,iq,ir,ll,meshsz,mmax
 real(dp),parameter :: eps=tol14**4,TOLJ=0.001_dp
 real(dp) :: arg,argn,bes
 real(dp) :: besp,qr
 real(dp) :: yp1,ypn
 character(len=100) :: msg
 type(pawrad_type) :: tmpmesh
!arrays
 real(dp),allocatable :: ff(:),gg(:),rr(:),rr2(:),rr2wf(:),rrwf(:),work(:)

!*************************************************************************

!Is mesh beginning with r=0 ?
 if (radmesh%rad(1)>tol10) then
   msg='Radial mesh cannot begin with r<>0!'
   MSG_BUG(msg)
 end if

 meshsz=size(wfll,1)
 if (meshsz>radmesh%mesh_size) then
   msg='wrong size for wfll!'
   MSG_BUG(msg)
 end if

!Init. temporary arrays and variables
 LIBPAW_ALLOCATE(ff,(meshsz))
 LIBPAW_ALLOCATE(gg,(meshsz))
 LIBPAW_ALLOCATE(rr,(meshsz))
 LIBPAW_ALLOCATE(rr2,(meshsz))
 LIBPAW_ALLOCATE(rrwf,(meshsz))
 LIBPAW_ALLOCATE(rr2wf,(meshsz))
 LIBPAW_ALLOCATE(work,(mqgrid))
 rr(1:meshsz) =radmesh%rad(1:meshsz)
 rr2(1:meshsz)=two_pi*rr(1:meshsz)*rr(1:meshsz)
 argn=two_pi*qgrid(mqgrid)
 mmax=meshsz

!Loop on (l,n) projectors
 iln0=0
 do ilmn=1,lmnmax
   iln=indlmn(5,ilmn)
   if(iln>iln0) then
     iln0=iln;ll=indlmn(1,ilmn)

     ir=meshsz
     do while (abs(wfll(ir,iln))<eps)
       ir=ir-1
     end do
     mmax=min(ir+1,meshsz)
     if (mmax/=radmesh%int_meshsz) then
       call pawrad_init(tmpmesh,mesh_size=meshsz,mesh_type=radmesh%mesh_type, &
&       rstep=radmesh%rstep,lstep=radmesh%lstep,r_for_intg=rr(mmax))
     else
       call pawrad_copy(radmesh,tmpmesh)
     end if

     rrwf(:) =rr (:)*wfll(:,iln)
     rr2wf(:)=rr2(:)*wfll(:,iln)

!    1-Compute f_l(0<q<qmax)
     if (mqgrid>2) then
       do iq=2,mqgrid-1
         arg=two_pi*qgrid(iq)
         do ir=1,mmax
           qr=arg*rr(ir)
           call paw_jbessel_4spline(bes,besp,ll,0,qr,TOLJ)
           ff(ir)=bes*rrwf(ir)
         end do
         call simp_gen(ffspl(iq,1,iln),ff,tmpmesh)
       end do
     end if

!    2-Compute f_l(q=0) and first derivative
     ffspl(1,1,iln)=zero;yp1=zero
     if (ll==0) then
       call simp_gen(ffspl(1,1,iln),rrwf,tmpmesh)
     end if
     if (ll==1) then
       call simp_gen(yp1,rr2wf,tmpmesh)
       yp1=yp1*third
     end if

!    3-Compute f_l(q=qmax) and first derivative
     if (mqgrid>1) then
!      if (ll==0.or.ll==1) then
       do ir=1,mmax
         qr=argn*rr(ir)
         call paw_jbessel_4spline(bes,besp,ll,1,qr,TOLJ)
         ff(ir)=bes*rrwf(ir)
         gg(ir)=besp*rr2wf(ir)
       end do
       call simp_gen(ffspl(mqgrid,1,iln),ff,tmpmesh)
       call simp_gen(ypn,gg,tmpmesh)
     else
       ypn=yp1
     end if

!    4-Compute second derivative of f_l(q)
     call paw_spline(qgrid,ffspl(:,1,iln),mqgrid,yp1,ypn,ffspl(:,2,iln))

     call pawrad_free(tmpmesh)

!    End loop on (l,n) projectors
   end if
 end do

 LIBPAW_DEALLOCATE(ff)
 LIBPAW_DEALLOCATE(gg)
 LIBPAW_DEALLOCATE(rr)
 LIBPAW_DEALLOCATE(rr2)
 LIBPAW_DEALLOCATE(rrwf)
 LIBPAW_DEALLOCATE(rr2wf)
 LIBPAW_DEALLOCATE(work)

end subroutine pawpsp_nl
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_lo
!! NAME
!! pawpsp_lo
!!
!! FUNCTION
!! Compute sine transform to transform from V(r) to q^2 V(q).
!! Computes integrals on (generalized) grid using corrected trapezoidal integration.
!!
!! INPUTS
!!  mqgrid=number of grid points in q from 0 to qmax.
!!  qgrid(mqgrid)=q grid values (bohr**-1).
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!  vloc(:)=V(r) on radial grid.
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
!!  yp1,ypn=derivatives of q^2 V(q) wrt q at q=0 and q=qmax (needed for spline fitter).
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_lo(epsatm,mqgrid,qgrid,q2vq,radmesh,vloc,yp1,ypn,zion)

!Arguments----------------------------------------------------------
!scalars
 integer,intent(in) :: mqgrid
 real(dp),intent(in) :: zion
 real(dp),intent(out) :: epsatm,yp1,ypn
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(in) :: vloc(:)
 real(dp),intent(out) :: q2vq(mqgrid)

!Local variables ------------------------------
!scalars
 integer :: iq,ir,irmax,mesh_size
 real(dp) :: arg,r0tor1,r1torm,rmtoin
 logical :: begin_r0
!arrays
 real(dp),allocatable :: ff(:),rvpz(:)

!************************************************************************

 mesh_size=size(vloc)
 irmax=pawrad_ifromr(radmesh,min(20._dp,radmesh%rmax))
 irmax=min(irmax,mesh_size)

!Particular case of a zero potential
 if (maxval(abs(vloc(1:irmax)))<=1.e-20_dp) then
   q2vq=zero;yp1=zero;ypn=zero;epsatm=zero
   return
 end if

 LIBPAW_ALLOCATE(ff,(mesh_size))
 LIBPAW_ALLOCATE(rvpz,(mesh_size))
 ff=zero;rvpz=zero

!Is mesh beginning with r=0 ?
 begin_r0=(radmesh%rad(1)<1.e-20_dp)

!Store r.V+Z
 do ir=1,irmax
   rvpz(ir)=radmesh%rad(ir)*vloc(ir)+zion
 end do

!===========================================
!=== Compute q^2 v(q) for q=0 separately
!===========================================

!Integral from 0 to r1 (only if r1<>0)
 r0tor1=zero;if (.not.begin_r0) &
& r0tor1=(zion*0.5_dp+radmesh%rad(1)*vloc(1)/3._dp)*radmesh%rad(1)**2

!Integral from r1 to rmax
 do ir=1,irmax
   if (abs(rvpz(ir))>1.e-20_dp) then
     ff(ir)=radmesh%rad(ir)*rvpz(ir)
   end if
 end do

 call simp_gen(r1torm,ff,radmesh)

!Integral from rmax to infinity
!This part is neglected... might be improved.
 rmtoin=zero

!Some of the three parts
 epsatm=four_pi*(r0tor1+r1torm+rmtoin)

 q2vq(1)=-zion/pi

!===========================================
!=== Compute q^2 v(q) for other q''s
!===========================================

!Loop over q values
 do iq=2,mqgrid
   arg=two_pi*qgrid(iq)

!  Integral from 0 to r1 (only if r1<>0)
   r0tor1=zero;if (.not.begin_r0) &
&   r0tor1=( vloc(1)/arg*sin(arg*radmesh%rad(1)) &
&   -rvpz(1)    *cos(arg*radmesh%rad(1)) +zion )/pi

!  Integral from r1 to rmax
   do ir=1,irmax
     if (abs(rvpz(ir))>1.e-20_dp) ff(ir)=sin(arg*radmesh%rad(ir))*rvpz(ir)
   end do
   call simp_gen(r1torm,ff,radmesh)

!  Integral from rmax to infinity
!  This part is neglected... might be improved.
   rmtoin=zero

!  Store q^2 v(q)
   q2vq(iq)=-zion/pi + two*qgrid(iq)*(r0tor1+r1torm+rmtoin)
 end do

!===========================================
!=== Compute derivatives of q^2 v(q)
!=== at ends of interval
!===========================================

!yp(0)=zero
 yp1=zero

!yp(qmax)=$ 2\int_0^\infty[(\sin(2\pi qmax r)+(2\pi qmax r)*\cos(2\pi qmax r)(r V(r)+Z) dr]$
 arg=two_pi*qgrid(mqgrid)

!Integral from 0 to r1 (only if r1<>0)
 r0tor1=zero;if (.not.begin_r0) &
& r0tor1=zion*radmesh%rad(1)                  *sin(arg*radmesh%rad(1)) &
& +three*radmesh%rad(1)*vloc(1)/arg         *cos(arg*radmesh%rad(1)) &
& +(radmesh%rad(1)**2-one/arg**2)*vloc(1)*sin(arg*radmesh%rad(1))

!Integral from r1 to rmax
 do ir=1,irmax
   if (abs(rvpz(ir))>1.e-20_dp) ff(ir)=( arg*radmesh%rad(ir)*cos(arg*radmesh%rad(ir)) &
&   +                    sin(arg*radmesh%rad(ir))) *rvpz(ir)
 end do
 call simp_gen(r1torm,ff,radmesh)

!Integral from rmax to infinity
!This part is neglected... might be improved.
 rmtoin=zero

!Some of the three parts
 ypn=two*(r0tor1+r1torm+rmtoin)

 LIBPAW_DEALLOCATE(ff)
 LIBPAW_DEALLOCATE(rvpz)

end subroutine pawpsp_lo
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_cg
!! NAME
!! pawpsp_cg
!!
!! FUNCTION
!! Compute sine transform to transform from n(r) to n(q).
!! Computes integrals on (generalized) grid using corrected trapezoidal integration.
!!
!! INPUTS
!!  mqgrid=number of grid points in q from 0 to qmax.
!!  qgrid(mqgrid)=q grid values (bohr**-1).
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!  nr(:)=n(r) on radial grid.
!!
!! OUTPUT
!!  dnqdq0= 1/q dn(q)/dq for q=0
!!  d2nqdq0 = Gives contribution of d2(tNcore(q))/d2q for q=0
!!            compute \int{(16/15)*pi^5*n(r)*r^6* dr}
!!{{\\ \begin{equation}
!!  nq(mqgrid)= n(q)
!!            = 4\pi\int[(\frac{\sin(2\pi q r)}{2\pi q r})(r^2 n(r))dr].
!!\end{equation} }}
!!  yp1,ypn=derivatives of n(q) wrt q at q=0 and q=qmax (needed for spline fitter).
!!
!! PARENTS
!!      m_dfpt_elt,m_pawpsp,m_psps
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_cg(dnqdq0,d2nqdq0,mqgrid,qgrid,nq,radmesh,nr,yp1,ypn)

!Arguments----------------------------------------------------------
!scalars
 integer,intent(in) :: mqgrid
 real(dp),intent(out) :: dnqdq0,d2nqdq0,yp1,ypn
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: nr(:)
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(out) :: nq(mqgrid)

!Local variables-------------------------------
!scalars
 integer :: iq,ir,mesh_size
 real(dp) :: aexp,arg,bexp,dn,r0tor1,r1torm,rm,rmtoin
 logical :: begin_r0
 !character(len=500) :: msg
!arrays
 real(dp),allocatable :: ff(:),rnr(:)

! *************************************************************************

 mesh_size=min(size(nr),radmesh%mesh_size)
 LIBPAW_ALLOCATE(ff,(mesh_size))
 LIBPAW_ALLOCATE(rnr,(mesh_size))
 ff=zero;rnr=zero

 do ir=1,mesh_size
   rnr(ir)=radmesh%rad(ir)*nr(ir)
 end do

!Is mesh beginning with r=0 ?
 begin_r0=(radmesh%rad(1)<1.d-20)

!Adjustment of an exponentional at r_max (n_exp(r)=aexp*Exp[-bexp*r])
 rm=radmesh%rad(mesh_size)
 dn=one/(12._dp*radmesh%stepint*radmesh%radfact(mesh_size)) &
& *( 3._dp*nr(mesh_size-4) &
&  -16._dp*nr(mesh_size-3) &
&  +36._dp*nr(mesh_size-2) &
&  -48._dp*nr(mesh_size-1) &
&  +25._dp*nr(mesh_size))
 if (dn<0._dp.and. &
& abs(radmesh%rad(mesh_size)*nr(mesh_size))>1.d-20) then
   bexp=-dn/nr(mesh_size)
   if (bexp * rm > 50._dp) then
     ! This solves the problem with the weird core charge used in v4[62] in which bexp x rm ~= 10^3
     !write(msg,"(a,es16.8)")"Tooooo large bexp * rm: ", bexp*rm, ", setting aexp to 0"
     !MSG_WARNING(msg)
     bexp=0.001_dp;aexp=zero
   else
     aexp=nr(mesh_size)*exp(bexp*rm)
     if (abs(aexp)<1.d-20) then
       bexp=0.001_dp;aexp=zero
     end if
   end if
 else
   bexp=0.001_dp;aexp=zero
 end if

!===========================================
!=== Compute n(q) for q=0 separately
!===========================================

!Integral from 0 to r1 (only if r1<>0)
 r0tor1=zero
 if (.not.begin_r0) r0tor1=(rnr(1)*radmesh%rad(1)**2)/3.d0

!Integral from r1 to rmax
 do ir=1,mesh_size
   if (abs(rnr(ir))>1.d-20) ff(ir)=rnr(ir)*radmesh%rad(ir)
 end do
 call simp_gen(r1torm,ff,radmesh)

!Integral from rmax to infinity
!This part is approximated using an exponential density aexp*Exp[-bexp*r]
!(formulae obtained with mathematica)
 rmtoin=aexp*exp(-bexp*rm)/bexp**3*(two+two*bexp*rm+bexp*bexp*rm*rm)

!Some of the three parts
 nq(1)=four_pi*(r0tor1+r1torm+rmtoin)

!===========================================
!=== Compute n(q) for other q''s
!===========================================

!Loop over q values
 do iq=2,mqgrid
   arg=two_pi*qgrid(iq)

!  Integral from 0 to r1 (only if r1<>0)
   r0tor1=zero;if (.not.begin_r0) &
&   r0tor1=nr(1)*(sin(arg*radmesh%rad(1))/arg/arg&
&   -radmesh%rad(1)*cos(arg*radmesh%rad(1))/arg)

!  Integral from r1 to rmax
   do ir=1,mesh_size
     if (abs(rnr(ir))>1.d-20) ff(ir)=sin(arg*radmesh%rad(ir))*rnr(ir)
   end do
   call simp_gen(r1torm,ff,radmesh)

!  Integral from rmax to infinity
!  This part is approximated using an exponential density aexp*Exp[-bexp*r]
!  (formulae obtained with mathematica)
   rmtoin=aexp*exp(-bexp*rm)/(arg**2+bexp**2)**2 &
&   *(arg*(two*bexp+arg**2*rm+bexp**2*rm)*cos(arg*rm) &
&   +(arg**2*(bexp*rm-one)+bexp**2*(bexp*rm+one))*sin(arg*rm))

!  Store q^2 v(q)
   nq(iq)=two/qgrid(iq)*(r0tor1+r1torm+rmtoin)
 end do

!===========================================
!=== Compute derivatives of n(q)
!=== at ends of interval
!===========================================

!yp(0)=zero
 yp1=zero

!yp(qmax)=$ 2\int_0^\infty[(-\sin(2\pi qmax r)+(2\pi qmax r)*\cos(2\pi qmax r) r n(r) dr]$
 arg=two_pi*qgrid(mqgrid)

!Integral from 0 to r1 (only if r1<>0)
 r0tor1=zero;if (.not.begin_r0) &
& r0tor1=two_pi*nr(1)*(3.d0*radmesh%rad(1)/arg /arg*cos(arg*radmesh%rad(1))+ &
& (radmesh%rad(1)**2/arg-3.0d0/arg**3)*sin(arg*radmesh%rad(1)))

!Integral from r1 to rmax
 do ir=1,mesh_size
   if (abs(rnr(ir))>1.d-20) ff(ir)=(two_pi*radmesh%rad(ir)*cos(arg*radmesh%rad(ir)) &
&   - sin(arg*radmesh%rad(ir))/qgrid(mqgrid)) *rnr(ir)
 end do
 call simp_gen(r1torm,ff,radmesh)

!Integral from rmax to infinity
!This part is approximated using an exponential density aexp*Exp[-bexp*r]
!(formulae obtained with mathematica)
 rmtoin=-one/(qgrid(mqgrid)*(arg**2+bexp**2)**3) &
& *aexp*exp(-bexp*rm) &
& *((arg**5*rm-two_pi*arg**4*qgrid(mqgrid)*rm*(bexp*rm-two) &
& +two*arg**3*bexp*(bexp*rm+one)+arg*bexp**3*(bexp*rm+two) &
& -four_pi*arg**2*bexp*qgrid(mqgrid)*(bexp**2*rm**2-three) &
& -two_pi*bexp**3*qgrid(mqgrid)*(bexp**2*rm**2+two*bexp*rm+two))*cos(arg*rm) &
& +(two*arg**2*bexp**3*rm+two_pi*arg**5*qgrid(mqgrid)*rm**2 &
& +arg**4*(bexp*rm-one)+bexp**4*(bexp*rm+one) &
& +four_pi*arg**3*qgrid(mqgrid)*(bexp**2*rm**2+two*bexp*rm-one) &
& +two_pi*arg*bexp**2*qgrid(mqgrid)*(bexp**2*rm**2+four*bexp*rm+6._dp))*sin(arg*rm))

!Some of the three parts
 ypn=two/qgrid(mqgrid)*(r0tor1+r1torm+rmtoin)

!===========================================
!=== Compute 1/q dn(q)/dq at q=0
!===========================================

!Integral from 0 to r1 (only if r1<>0)
 r0tor1=zero
 if (.not.begin_r0) r0tor1=(rnr(1)*radmesh%rad(1)**4)/5.d0

!Integral from r1 to rmax
 do ir=1,mesh_size
   if (abs(rnr(ir))>1.d-20) ff(ir)=rnr(ir)*radmesh%rad(ir)**3
 end do
 call simp_gen(r1torm,ff,radmesh)

!Integral from rmax to infinity
!This part is approximated using an exponential density aexp*Exp[-bexp*r]
!(formulae obtained with mathematica)
 rmtoin=aexp*exp(-bexp*rm)/bexp**5 &
& *(24._dp+24._dp*bexp*rm+12._dp*bexp**2*rm**2+four*bexp**3*rm**3+bexp**4*rm**4)

!Some of the three parts
 dnqdq0=-(2.d0/3.d0)*two_pi**3*(r0tor1+r1torm+rmtoin)

 LIBPAW_DEALLOCATE(ff)
 LIBPAW_DEALLOCATE(rnr)

 d2nqdq0 = 1_dp

end subroutine pawpsp_cg
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_read
!! NAME
!!  pawpsp_read
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!   File format of formatted PAW psp input (the 3 first lines
!!   have already been read in calling -pspatm- routine) :
!!   (1) title (character) line
!!   (2) psps%znuclpsp(ipsp), zion, pspdat
!!   (3) pspcod, pspxc, lmax, lloc, mmax, r2well
!!   (4) psp_version, creatorID
!!   (5) basis_size, lmn_size
!!   (6) orbitals (for l=1 to basis_size)
!!   (7) number_of_meshes
!!   For imsh=1 to number_of_meshes
!!   (8)  mesh_index, mesh_type ,mesh_size, rad_step[, log_step]
!!   (9) r_cut(SPH)
!!   (10) shape_type, r_shape[, shapefunction arguments]
!!   For iln=1 to basis_size
!!   (11) comment(character)
!!   (12) radial mesh index for phi
!!   (13) phi(r) (for ir=1 to phi_meshsz)
!!   For iln=1 to basis_size
!!   (14) comment(character)
!!   (15) radial mesh index for tphi
!!   (16) tphi(r) (for ir=1 to phi_mesh_size)
!!   For iln=1 to basis_size
!!   (17) comment(character)
!!   (18) radial mesh index for tproj
!!   (19) tproj(r) (for ir=1 to proj_mesh_size)
!!   (20) comment(character)
!!   (21) radial mesh index for core_density
!!   (22) core_density (for ir=1 to core_mesh_size)
!!   (23) comment(character)
!!   (24) radial mesh index for pseudo_core_density
!!   (25) tcore_density (for ir=1 to core_mesh_size)
!!   (26) comment(character)
!!   (27) Dij0 (for ij=1 to lmn_size*(lmn_size+1)/2)
!!   (28) comment(character)
!!   (29) Rhoij0 (for ij=1 to lmn_size*(lmn_size+1)/2)
!!   (30) comment(character)
!!   (31) radial mesh index for Vloc, format of Vloc (0=Vbare, 1=VH(tnzc), 2=VH(tnzc) without nhat in XC)
!!   (32) Vloc(r) (for ir=1 to vloc_mesh_size)
!!   ===== Following lines only if shape_type=-1 =====
!!   For il=1 to 2*max(orbitals)+1
!!   (33) comment(character)
!!   (34) radial mesh index for shapefunc
!!   (35) shapefunc(r)*gnorm(l)*r**l (for ir=1 to shape_mesh_size)
!!   (36) comment(character)
!!   (37) radial mesh index for pseudo_valence_density
!!   (38) tvale(r) (for ir=1 to vale_mesh_size)
!!
!!   Comments:
!!   * psp_version= ID of PAW_psp version
!!   4 characters string of the form 'pawn' (with n varying)
!!   * creatorID= ID of psp generator
!!   creatorid=1xyz : psp generated from Holzwarth AtomPAW generator version x.yz
!!   creatorid=2xyz : psp generated from Vanderbilt ultra-soft generator version x.yz
!!   creatorid=-1: psp for tests (for developpers only)
!!   * mesh_type= type of radial mesh
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!   * radial shapefunction type
!!   shape_type=-1 ; gl(r)=numeric (read from psp file)
!!   shape_type= 1 ; gl(r)=k(r).r^l; k(r)=exp[-(r/sigma)**lambda]
!!   shape_type= 2 ; gl(r)=k(r).r^l; k(r)=[sin(pi*r/rshp)/(pi*r/rshp)]**2 if r<=rshp
!!   shape_type= 3 ; gl(r)=Alpha(1,l)*jl(q(1,l)*r)+Alpha(2,l)*jl(q(2,l)*r) for each l
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_read(core_mesh,funit,imainmesh,lmax,&
& ncore,nmesh,pawrad,pawtab,pspversion,radmesh,save_core_msz,&
& tcoretau,tncore,tnvale,tproj,tproj_mesh,usexcnhat_in,usexcnhat_out,vale_mesh,&
& vlocopt,vlocr,vloc_mesh,znucl)

!Arguments ------------------------------------
 integer,intent(in):: funit,lmax,usexcnhat_in
 integer,intent(out) :: imainmesh,pspversion,usexcnhat_out,vlocopt
 logical,intent(in) :: save_core_msz
 real(dp),intent(in):: znucl
!arrays
 real(dp),pointer :: ncore(:),tcoretau(:),tncore(:),tnvale(:),tproj(:,:),vlocr(:)
 type(pawrad_type),intent(inout) :: pawrad
 type(pawrad_type),intent(out)::core_mesh,tproj_mesh,vale_mesh,vloc_mesh
 type(pawrad_type),pointer :: radmesh(:)
 type(pawtab_type),intent(inout) :: pawtab
 integer,intent(out)::nmesh

!Local variables-------------------------------
 integer :: creatorid,imsh
 integer :: icoremesh,ishpfmesh,ivalemesh,ivlocmesh
 integer :: ib,il,ilm,ilmn,iln,iprojmesh
 integer :: ii,ir,iread1,iread2,jj
 integer :: msz,pngau_,ptotgau_
 real(dp):: rc,rread1,rread2
 real(dp) :: yp1,ypn
!arrays
 integer,allocatable :: nprj(:)
 real(dp),allocatable :: shpf(:,:),val(:),vhnzc(:)
 real(dp),allocatable :: work1(:),work2(:),work3(:),work4(:)
 character :: blank=' ',numb=' '
 character(len=80) :: pspline
 character(len=500) :: msg,submsg
 logical :: read_gauss=.false.
 type(pawrad_type)::shpf_mesh

! *************************************************************************

!==========================================================
!Read lines 4 to 11 of the header

!This is important for BigDFT in standalone mode
 call pawpsp_read_header_2(funit,pspversion,pawtab%basis_size,pawtab%lmn_size)

!Check pspversion for wvl-paw
 if(pspversion<4 .and. pawtab%has_wvl>0)then
   write(msg, '(a,i2,a,a)' )&
&   'In reading atomic psp file, finds pspversion=',pspversion,ch10,&
&   'For WVL-PAW, pspversion >= 4 is required.'
   MSG_BUG(msg)
 end if


!Have to maintain compatibility with Abinit v4.2.x
 if (pspversion==1) then
   LIBPAW_ALLOCATE(pawtab%orbitals,(pawtab%basis_size))
   read(funit,*) (pawtab%orbitals(ib), ib=1,pawtab%basis_size)
   pawtab%l_size=2*maxval(pawtab%orbitals)+1
   nmesh=3
   LIBPAW_DATATYPE_ALLOCATE(radmesh,(nmesh))
   read(funit,'(a80)') pspline
   radmesh(1)%lstep=zero
   read(unit=pspline,fmt=*,err=10,end=10) radmesh(1)%mesh_type,&
&   radmesh(1)%rstep,radmesh(1)%lstep
   10 read(funit,*) pawtab%rpaw
   read(funit,*) radmesh(1)%mesh_size,radmesh(2)%mesh_size,&
&   radmesh(3)%mesh_size
   read(funit,'(a80)') pspline
   pawtab%shape_lambda=-1;pawtab%shape_sigma=1.d99
   read(unit=pspline,fmt=*,err=11,end=11) pawtab%shape_type,&
&   pawtab%shape_lambda,pawtab%shape_sigma
   11 read(funit,*) creatorid
   if (pawtab%shape_type==3) pawtab%shape_type=-1
   radmesh(2)%mesh_type=radmesh(1)%mesh_type
   radmesh(3)%mesh_type=radmesh(1)%mesh_type
   radmesh(2)%rstep=radmesh(1)%rstep
   radmesh(3)%rstep=radmesh(1)%rstep
   radmesh(2)%lstep=radmesh(1)%lstep
   radmesh(3)%lstep=radmesh(1)%lstep
 else

!  Here psp file for Abinit 4.3+
   LIBPAW_ALLOCATE(pawtab%orbitals,(pawtab%basis_size))
   read(funit,*) (pawtab%orbitals(ib), ib=1,pawtab%basis_size)
   pawtab%l_size=2*maxval(pawtab%orbitals)+1
   read(funit,*) nmesh
   LIBPAW_DATATYPE_ALLOCATE(radmesh,(nmesh))
   do imsh=1,nmesh
     rread2=zero
     read(funit,'(a80)') pspline
     read(unit=pspline,fmt=*,err=20,end=20) ii,iread1,iread2,rread1,rread2
     20 continue
     if (ii<=nmesh) then
       radmesh(ii)%mesh_type=iread1
       radmesh(ii)%mesh_size=iread2
       radmesh(ii)%rstep=rread1
       radmesh(ii)%lstep=rread2
     else
       write(msg, '(3a)' )&
&       'Index of mesh out of range !',ch10,&
&       'Action : check your pseudopotential file.'
       MSG_ERROR(msg)
     end if
   end do
   read(funit,*) pawtab%rpaw
   read(funit,'(a80)') pspline
   read(unit=pspline,fmt=*) pawtab%shape_type
   pawtab%shape_lambda=-1;pawtab%shape_sigma=1.d99
 end if

!Initialize radial meshes
 do imsh=1,nmesh
   call pawrad_init(radmesh(imsh))
 end do

!==========================================================
!Initialize various dims and indexes

 pawtab%l_size=2*maxval(pawtab%orbitals)+1
 pawtab%lmn2_size=pawtab%lmn_size*(pawtab%lmn_size+1)/2
 pawtab%ij_size=pawtab%basis_size*(pawtab%basis_size+1)/2
 pawtab%usexcnhat=usexcnhat_in

!indlmn calculation (indices for (l,m,n) basis)
 if (allocated(pawtab%indlmn)) then
   LIBPAW_DEALLOCATE(pawtab%indlmn)
 end if
 LIBPAW_ALLOCATE(pawtab%indlmn,(6,pawtab%lmn_size))
 LIBPAW_BOUND1_ALLOCATE(nprj,BOUNDS(0,maxval(pawtab%orbitals)))
 pawtab%indlmn(:,:)=0
 ilmn=0;iln=0;nprj=0
 do ib=1,pawtab%basis_size
   il=pawtab%orbitals(ib)
   nprj(il)=nprj(il)+1
   iln=iln+1
   do ilm=1,2*il+1
     pawtab%indlmn(1,ilmn+ilm)=il
     pawtab%indlmn(2,ilmn+ilm)=ilm-(il+1)
     pawtab%indlmn(3,ilmn+ilm)=nprj(il)
     pawtab%indlmn(4,ilmn+ilm)=il*il+ilm
     pawtab%indlmn(5,ilmn+ilm)=iln
     pawtab%indlmn(6,ilmn+ilm)=1
   end do
   ilmn=ilmn+2*il+1
 end do
 LIBPAW_DEALLOCATE(nprj)
!Are ilmn (found here) and pawtab%lmn_size compatibles ?
 if (ilmn/=pawtab%lmn_size) then
   write(msg, '(a,a,a,a,a)' )&
&   'Calculated lmn size differs from',ch10,&
&   'lmn_size read from pseudo !',ch10,&
&   'Action: check your pseudopotential file.'
   MSG_ERROR(msg)
 end if

!==========================================================
!Here reading shapefunction parameters

!Shapefunction parameters for Abinit 4.3...4.5
 if (pspversion==2) then
   if (pawtab%shape_type==1) read(unit=pspline,fmt=*) ii,pawtab%shape_lambda,pawtab%shape_sigma
   if (pawtab%shape_type==3) pawtab%shape_type=-1
   pawtab%rshp=zero

!Shapefunction parameters for Abinit 4.6+
 else if (pspversion>=3) then
   pawtab%rshp=zero
   if (pawtab%shape_type==-1) read(unit=pspline,fmt=*,err=21,end=21) ii,pawtab%rshp
   if (pawtab%shape_type== 1) read(unit=pspline,fmt=*,err=21,end=21) ii,pawtab%rshp, &
&   pawtab%shape_lambda,pawtab%shape_sigma
   if (pawtab%shape_type== 2) read(unit=pspline,fmt=*,err=21,end=21) ii,pawtab%rshp
   if (pawtab%shape_type== 3) read(unit=pspline,fmt=*,err=21,end=21) ii,pawtab%rshp
 end if
 21 continue
!If shapefunction type is gaussian, check exponent
 if (pawtab%shape_type==1) then
   if (pawtab%shape_lambda<2) then
     write(msg, '(3a)' )&
&     'For a gaussian shape function, exponent lambda must be >1 !',ch10,&
&     'Action: check your psp file.'
     MSG_ERROR(msg)
   end if
 end if
!If shapefunction type is Bessel, deduce here its parameters from rc
 if (pawtab%shape_type==3) then
   LIBPAW_ALLOCATE(pawtab%shape_alpha,(2,pawtab%l_size))
   LIBPAW_ALLOCATE(pawtab%shape_q,(2,pawtab%l_size))
   rc=pawtab%rshp;if (rc<1.d-8) rc=pawtab%rpaw
   do il=1,pawtab%l_size
     call atompaw_shapebes(pawtab%shape_alpha(1:2,il),pawtab%shape_q(1:2,il),il-1,rc)
   end do
 end if

!==========================================================
!Mirror pseudopotential parameters to the output and log files

 write(msg,'(a,i1)')' Pseudopotential format is: paw',pspversion
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')
 write(msg,'(2(a,i3),a,64i4)') &
& ' basis_size (lnmax)=',pawtab%basis_size,' (lmn_size=',&
& pawtab%lmn_size,'), orbitals=',pawtab%orbitals(1:pawtab%basis_size)
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')
 write(msg,'(a,f11.8)')' Spheres core radius: rc_sph=',pawtab%rpaw
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')
 write(msg,'(a,i1,a)')' ',nmesh,' radial meshes are used:'
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')
 do imsh=1,nmesh
   if (radmesh(imsh)%mesh_type==1) &
&   write(msg,'(a,i1,a,i4,a,g12.5)') &
&   '  - mesh ',imsh,': r(i)=step*(i-1), size=',radmesh(imsh)%mesh_size,&
&   ' , step=',radmesh(imsh)%rstep
   if (radmesh(imsh)%mesh_type==2) &
&   write(msg,'(a,i1,a,i4,2(a,g12.5))') &
&   '  - mesh ',imsh,': r(i)=AA*[exp(BB*(i-1))-1], size=',radmesh(imsh)%mesh_size,&
&   ' , AA=',radmesh(imsh)%rstep,' BB=',radmesh(imsh)%lstep
   if (radmesh(imsh)%mesh_type==3) &
&   write(msg,'(a,i1,a,i4,2(a,g12.5))') &
&   '  - mesh ',imsh,': r(i)=AA*exp(BB*(i-2)), size=',radmesh(imsh)%mesh_size,&
&   ' , AA=',radmesh(imsh)%rstep,' BB=',radmesh(imsh)%lstep
   if (radmesh(imsh)%mesh_type==4) &
&   write(msg,'(a,i1,a,i4,a,g12.5)') &
&   '  - mesh ',imsh,': r(i)=-AA*ln(1-(i-1)/n), n=size=',radmesh(imsh)%mesh_size,&
&   ' , AA=',radmesh(imsh)%rstep
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 end do
 if (pawtab%shape_type==-1) then
   write(msg,'(a)')&
   ' Shapefunction is NUMERIC type: directly read from atomic data file'
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 end if
 if (pawtab%shape_type==1) then
   write(msg,'(2a,a,f6.3,a,i3)')&
&   ' Shapefunction is EXP type: shapef(r)=exp(-(r/sigma)**lambda)',ch10,&
&   '                            with sigma=',pawtab%shape_sigma,' and lambda=',pawtab%shape_lambda
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 end if
 if (pawtab%shape_type==2) then
   write(msg,'(a)')&
   ' Shapefunction is SIN type: shapef(r)=[sin(pi*r/rshp)/(pi*r/rshp)]**2'
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 end if
 if (pawtab%shape_type==3) then
   write(msg,'(a)')&
&   ' Shapefunction is BESSEL type: shapef(r,l)=aa(1,l)*jl(q(1,l)*r)+aa(2,l)*jl(q(2,l)*r)'
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 end if
 if (pawtab%rshp<1.d-8) then
   write(msg,'(a)') ' Radius for shape functions = sphere core radius'
 else
   write(msg,'(a,f11.8)') ' Radius for shape functions = ',pawtab%rshp
 end if
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

!==========================================================
!Perfom tests

!Are lmax and orbitals compatibles ?
 if (lmax/=maxval(pawtab%orbitals)) then
   write(msg, '(a,a,a)' )&
&   'lmax /= MAX(orbitals) !',ch10,&
&   'Action: check your pseudopotential file.'
   MSG_ERROR(msg)
 end if

!Only mesh_type=1,2, 3 or 4 allowed
 do imsh=1,nmesh
   if (radmesh(imsh)%mesh_type>4) then
     write(msg, '(a,a,a)' )&
&     'Only mesh types 1,2, 3 or 4 allowed !',ch10,&
&     'Action : check your pseudopotential or input file.'
     MSG_ERROR(msg)
   end if
 end do

!==========================================================
!Read tabulated atomic data

!---------------------------------
!Read wave-functions (phi)
 do ib=1,pawtab%basis_size
   read (funit,*)
   if (pspversion==1) iread1=1
   if (pspversion>1) read (funit,*) iread1
   if (ib==1) then
     call pawrad_free(pawrad)
     call pawrad_init(pawrad,mesh_size=radmesh(iread1)%mesh_size,mesh_type=radmesh(iread1)%mesh_type,&
&     rstep=radmesh(iread1)%rstep,lstep=radmesh(iread1)%lstep,r_for_intg=pawtab%rpaw)
     pawtab%partialwave_mesh_size=pawrad%mesh_size
     pawtab%mesh_size=pawrad_ifromr(pawrad,pawtab%rpaw)+5
     pawtab%mesh_size=min(pawtab%mesh_size,pawrad%mesh_size)
     if (pawtab%mesh_size>pawrad%mesh_size-2) pawtab%mesh_size=pawrad%mesh_size
     imainmesh=iread1
     LIBPAW_ALLOCATE(pawtab%phi,(pawtab%partialwave_mesh_size,pawtab%basis_size))
   else if (iread1/=imainmesh) then
     write(msg, '(a,a,a)' )&
&     'All Phi and tPhi must be given on the same radial mesh !',ch10,&
&     'Action: check your pseudopotential file.'
     MSG_ERROR(msg)
   end if
   read (funit,*) (pawtab%phi(ir,ib),ir=1,pawtab%partialwave_mesh_size)
 end do

!---------------------------------
!Read pseudo wave-functions (tphi)
 LIBPAW_ALLOCATE(pawtab%tphi,(pawtab%partialwave_mesh_size,pawtab%basis_size))
 do ib=1,pawtab%basis_size
   read (funit,*)
   if (pspversion==1) iread1=1
   if (pspversion>1) read (funit,*) iread1
   if (iread1/=imainmesh) then
     write(msg, '(a,a,a)' )&
&     'All Phi and tPhi must be given on the same radial mesh !',ch10,&
&     'Action: check your pseudopotential file.'
     MSG_ERROR(msg)
   end if
   read (funit,*) (pawtab%tphi(ir,ib),ir=1,pawtab%partialwave_mesh_size)
 end do
 write(msg,'(a,i1)') &
& ' Radial grid used for partial waves is grid ',imainmesh
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

!---------------------------------
!Read projectors (tproj)
 do ib=1,pawtab%basis_size
   read (funit,*)
   if (pspversion==1) iread1=2
   if (pspversion>1) read (funit,*) iread1
   if (ib==1) then
     iprojmesh=iread1
     call pawrad_copy(radmesh(iprojmesh),tproj_mesh)
     LIBPAW_POINTER_ALLOCATE(tproj,(tproj_mesh%mesh_size,pawtab%basis_size))
   else if (iread1/=iprojmesh) then
     write(msg, '(a,a,a)' )&
&     'All tprojectors must be given on the same radial mesh !',ch10,&
&     'Action: check your pseudopotential file.'
     MSG_ERROR(msg)
   end if
!  read projectors from a mesh
   read (funit,*) (tproj(ir,ib),ir=1,tproj_mesh%mesh_size)
 end do
 write(msg,'(a,i2)') &
& ' Radial grid used for projectors is grid ',iprojmesh
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

!---------------------------------
!Read gaussian projectors for wavelets
!  -- only if pawtab%has_wvl flag is on
!  -- if not, we skip the lines
 read(funit,'(a80)') pspline
 if(index(trim(pspline),'GAUSSIAN')/=0) read_gauss=.true.
 if (read_gauss) then
   if (pawtab%has_wvl>0) then
     call wvlpaw_allocate(pawtab%wvl)
     jj=0
     do ib=1,pawtab%basis_size
       if(ib/=1) read(funit,*) pspline
!      read Gaussian coefficients
       read(funit,*) pngau_, ptotgau_ !total number of gaussians
       if(ib==1) then
         pawtab%wvl%ptotgau=ptotgau_
         LIBPAW_ALLOCATE(pawtab%wvl%pngau,(pawtab%basis_size))
         LIBPAW_ALLOCATE(pawtab%wvl%parg,(2,pawtab%wvl%ptotgau))
         LIBPAW_ALLOCATE(pawtab%wvl%pfac,(2,pawtab%wvl%ptotgau))
       else
         if(pawtab%wvl%ptotgau/=ptotgau_) then
           write(msg,'(3a)')&
&           'Total number of gaussians, should be the same for all projectors !',ch10,&
&           'Action: check your pseudopotential file.'
           MSG_ERROR(msg)
         end if
       end if !ib==1
       read(funit,*)(pawtab%wvl%parg(:,ii),ii=jj+1,jj+pngau_)
       read(funit,*)(pawtab%wvl%pfac(:,ii),ii=jj+1,jj+pngau_)
       pawtab%wvl%pngau(ib)=pngau_
       jj=jj+pngau_
     end do
     pawtab%has_wvl=2
   else
!    If pawtab%has_wvl=0, we skip the lines
     do ib=1,pawtab%basis_size
       if(ib/=1) read(funit,*)
       read(funit,*) pngau_, ptotgau_
       LIBPAW_ALLOCATE(val, (pngau_  *2))
       read(funit,*) val
       read(funit,*) val
       LIBPAW_DEALLOCATE(val)
     end do
   end if
 end if

!---------------------------------
!Read core density (coredens)
 if(read_gauss) read (funit,*) !if not read_gauss, this line was already read
 if (pspversion==1) iread1=1
 if (pspversion>1) read (funit,*) iread1
 icoremesh=iread1
 call pawrad_copy(radmesh(icoremesh),core_mesh)
 if ((radmesh(icoremesh)%mesh_type/=pawrad%mesh_type).or.&
& (radmesh(icoremesh)%rstep    /=pawrad%rstep)    .or.&
& (radmesh(icoremesh)%lstep    /=pawrad%lstep)) then
   write(msg, '(a,a,a,a,a)' )&
&   'Ncore must be given on a radial mesh with the same',ch10,&
&   'type and step(s) than the main radial mesh (mesh for Phi) !',ch10,&
&   'Action: check your pseudopotential file.'
   MSG_ERROR(msg)
 end if
 LIBPAW_POINTER_ALLOCATE(ncore,(core_mesh%mesh_size))
 read (funit,*) (ncore(ir),ir=1,core_mesh%mesh_size)

!Construct and save VH[z_NC] if requested
 if (pawtab%has_vhnzc==1) then
   LIBPAW_ALLOCATE(pawtab%VHnZC,(pawtab%mesh_size))
   LIBPAW_ALLOCATE(vhnzc,(core_mesh%mesh_size))
   call atompaw_vhnzc(ncore,core_mesh,vhnzc,znucl)
   pawtab%VHnZC(1:pawtab%mesh_size)=vhnzc(1:pawtab%mesh_size)
   pawtab%has_vhnzc=2
   LIBPAW_DEALLOCATE(vhnzc)
 end if

 pawtab%core_mesh_size=pawrad%mesh_size
 if(save_core_msz) pawtab%core_mesh_size=core_mesh%mesh_size
 LIBPAW_ALLOCATE(pawtab%coredens,(pawtab%core_mesh_size))
 pawtab%rcore=core_mesh%rad(pawtab%core_mesh_size)
 pawtab%coredens(1:pawtab%core_mesh_size)=ncore(1:pawtab%core_mesh_size)

!---------------------------------
!Read pseudo core density (tcoredens)
 if(save_core_msz)  then
   LIBPAW_ALLOCATE(pawtab%tcoredens,(pawtab%core_mesh_size,6))
 else
   LIBPAW_ALLOCATE(pawtab%tcoredens,(pawtab%core_mesh_size,1))
 end if
 pawtab%tcoredens=zero
 read (funit,*)
 if (pspversion==1) iread1=1
 if (pspversion>1) read (funit,*) iread1
 if (iread1/=icoremesh) then
   write(msg, '(a,a,a,a,a,a,a,a)' )&
&   'Pseudized core density (tNcore) must be given',ch10,&
&   'on the same radial mesh as core density (Ncore) !',ch10,&
&   'Action: check your pseudopotential file.'
   MSG_ERROR(msg)
 end if
 LIBPAW_POINTER_ALLOCATE(tncore,(core_mesh%mesh_size))
 read (funit,*) (tncore(ir),ir=1,core_mesh%mesh_size)
 if (maxval(abs(tncore(:)))<tol6) then
   pawtab%usetcore=0
 else
   pawtab%usetcore=1
   pawtab%tcoredens(1:pawtab%core_mesh_size,1)=tncore(1:pawtab%core_mesh_size)
 end if
 write(msg,'(a,i1)') &
& ' Radial grid used for (t)core density is grid ',icoremesh
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

!---------------------------------
!Read frozen part of Dij terms (dij0)
 LIBPAW_ALLOCATE(pawtab%dij0,(pawtab%lmn2_size))
 read (funit,*)
 read (funit,*) (pawtab%dij0(ib),ib=1,pawtab%lmn2_size)

!---------------------------------
!Read initial guess of rhoij (rhoij0)
 LIBPAW_ALLOCATE(pawtab%rhoij0,(pawtab%lmn2_size))
 read (funit,*)
 read (funit,*) (pawtab%rhoij0(ib),ib=1,pawtab%lmn2_size)

!---------------------------------
!Read local pseudopotential=Vh(tn_zc) or Vbare
 read (funit,*)
 if (pspversion==1) ivlocmesh=3
 vlocopt=1
 if (pspversion==2) then
   read (funit,*) ivlocmesh
 else if (pspversion>2) then
!  read (funit,fmt=*,err=30,end=30) ivlocmesh,vlocopt
   msg=blank
   read (funit,fmt='(a)') msg
   read (msg,fmt=*) ivlocmesh
   write(numb,'(i1)')ivlocmesh
   ii=index(msg,numb)
   if(len_trim(trim(msg(ii+1:)))/=0)then
     submsg=trim(msg(ii+1:))
     if(len_trim(submsg)/=0)then
       do ii=1,len_trim(submsg)
         numb=submsg(ii:ii)
         if(numb==blank)cycle
         jj=index('0123456789',numb)
         if(jj<1 .or. jj>10)exit
         vlocopt=jj-1
       end do
     end if
   end if
 end if
 usexcnhat_out=0;if (vlocopt==1) usexcnhat_out=1
 call pawrad_copy(radmesh(ivlocmesh),vloc_mesh)
 LIBPAW_POINTER_ALLOCATE(vlocr,(vloc_mesh%mesh_size))
 read (funit,*) (vlocr(ir),ir=1,vloc_mesh%mesh_size)
 write(msg,'(a,i1)') &
& ' Radial grid used for Vloc is grid ',ivlocmesh
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

!---------------------------------
!Eventually read "numeric" shapefunctions (if shape_type=-1)
 if (pawtab%shape_type==-1) then
   LIBPAW_ALLOCATE(pawtab%shapefunc,(pawtab%mesh_size,pawtab%l_size))
   do il=1,pawtab%l_size
     read (funit,*)
     if (pspversion==1) iread1=1
     if (pspversion>1) read (funit,*) iread1
     if (il==1) then
       call pawrad_copy(radmesh(iread1),shpf_mesh)
       ishpfmesh=iread1
       LIBPAW_ALLOCATE(shpf,(shpf_mesh%mesh_size,pawtab%l_size))
     else if (iread1/=ishpfmesh) then
       write(msg, '(a,a,a)' )&
&       'All shape functions must be given on the same radial mesh !',ch10,&
&       'Action: check your pseudopotential file.'
       MSG_ERROR(msg)
     end if
     read (funit,*) (shpf(ir,il),ir=1,shpf_mesh%mesh_size)
   end do
   write(msg,'(a,i1)') &
&   ' Radial grid used for shape functions is grid ',iread1
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')

!  Has to spline shape functions if mesh is not the "main" mesh
   if (ishpfmesh/=imainmesh) then
     msz=shpf_mesh%mesh_size
     LIBPAW_ALLOCATE(work1,(msz))
     LIBPAW_ALLOCATE(work2,(msz))
     LIBPAW_ALLOCATE(work3,(msz))
     LIBPAW_ALLOCATE(work4,(pawtab%mesh_size))
     work3(1:pawtab%mesh_size)=shpf_mesh%rad(1:pawtab%mesh_size)
     work4(1:pawtab%mesh_size)=pawrad%rad(1:pawtab%mesh_size)
     do il=1,pawtab%l_size
       call bound_deriv(shpf(1:msz,il),shpf_mesh,msz,yp1,ypn)
       call paw_spline(work3,shpf(:,il),msz,yp1,ypn,work1)
       call paw_splint(msz,work3,shpf(:,il),work1,pawtab%mesh_size,work4,pawtab%shapefunc(:,il))
     end do
     LIBPAW_DEALLOCATE(work1)
     LIBPAW_DEALLOCATE(work2)
     LIBPAW_DEALLOCATE(work3)
     LIBPAW_DEALLOCATE(work4)
   else
     pawtab%shapefunc(:,:)=shpf(:,:)
   end if
   LIBPAW_DEALLOCATE(shpf)
   call pawrad_free(shpf_mesh)
 end if

!---------------------------------
!Read pseudo valence density (if psp version >=4)
 if (pspversion>=4) then
   read (funit,*)
   read (funit,*) iread1
   ivalemesh=iread1
   call pawrad_copy(radmesh(iread1),vale_mesh)
   LIBPAW_POINTER_ALLOCATE(tnvale,(vale_mesh%mesh_size))
   read (funit,*) (tnvale(ir),ir=1,vale_mesh%mesh_size)
   pawtab%has_tvale=1
   write(msg,'(a,i1)') &
&   ' Radial grid used for pseudo valence density is grid ',ivalemesh
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 else
   pawtab%has_tvale=0
   LIBPAW_POINTER_ALLOCATE(tnvale,(0))
 end if

!---------------------------------
!Initialize (to zero) kinetic energy densities (for testing purpose)
 if (pawtab%has_coretau>0) then
   write(msg,'(5a)' )&
&   'Kinetic energy density is requested but the core kinetic energy density',ch10,&
&   'is not present in the pseudopotential file!',ch10,&
&   'We assume that it is zero (for testing purpose).'
   MSG_WARNING(msg)
   pawtab%coretau_mesh_size=pawtab%mesh_size
   if(save_core_msz) pawtab%coretau_mesh_size=core_mesh%mesh_size
   LIBPAW_ALLOCATE(pawtab%coretau,(pawtab%coretau_mesh_size))
   LIBPAW_ALLOCATE(pawtab%tcoretau,(pawtab%coretau_mesh_size))
   LIBPAW_POINTER_ALLOCATE(tcoretau,(core_mesh%mesh_size))
   pawtab%rcoretau=core_mesh%rad(pawtab%coretau_mesh_size)
   pawtab%coretau=zero ; pawtab%tcoretau=zero ; tcoretau=zero
 endif

end subroutine pawpsp_read
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_read_corewf
!! NAME
!!  pawpsp_read_corewf
!!
!! FUNCTION
!!
!! INPUTS
!! [filename]= (optional) core WF file name
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_paw_optics,m_positron
!!
!! CHILDREN
!!
!! SOURCE
subroutine pawpsp_read_corewf(energy_cor,indlmn_core,lcor,lmncmax,ncor,nphicor,radmesh,phi_cor,&
&                             filename) ! optional argument

!Arguments ------------------------------------
 integer,intent(out) :: lmncmax,nphicor
 character(len=*),optional :: filename
!arrays
 integer,allocatable,intent(inout) :: indlmn_core(:,:),lcor(:),ncor(:)
 real(dp),allocatable,intent(inout) :: phi_cor(:,:),energy_cor(:)
 type(pawrad_type),intent(in) :: radmesh

!Local variables-------------------------------
 integer :: ib,i1,i2,il,ilm,ilmn,iln,ios,jln,nmesh,npts,unt
 real(dp) :: noccor,r1,r2
 logical :: ex,oldformat,usexml
 character(len=8) :: dum,dum1,dum2,dum3,dum4
 character(len=80) :: fline
 character(len=500) :: msg
 character(len=fnlen) :: filename_

!arrays
 integer,allocatable :: meshtp(:),meshsz(:)
 real(dp),allocatable :: rad(:),radstp(:),work(:)
 real(dp),allocatable :: logstp(:),phitmp(:)
 type(pawrad_type) :: tmpmesh

! ************************************************************************

!Check for core WF file existence and XML format
 usexml=.false.;oldformat=.false.
 if (present(filename)) then
!  Core WF file given as optional argument
   filename_=filename;ex=.false.
   inquire(file=trim(filename_),iostat=ios,exist=ex)
   if (ios/=0) then
     write(msg,'(2a)') 'INQUIRE returns an error for file ',trim(filename_)
     MSG_ERROR(msg)
   end if
   if (.not.ex) then
     write(msg,'(3a)') 'This file does not exist: ',trim(filename_),'!'
     MSG_ERROR(msg)
   end if
   unt = libpaw_get_free_unit()
   open(unit=unt,file=trim(filename_),form='formatted',status='old', action="read")
   read(unt,*) fline
   close(unt)
   usexml=(fline(1:5)=='<?xml')
 else
!  Core WF file: new format
   filename_='corewf.abinit';ex=.false.
   inquire(file=trim(filename_),iostat=ios,exist=ex)
   if (ios/=0) then
     write(msg,'(3a)') 'INQUIRE returns an error for file ',trim(filename_),'!'
     MSG_ERROR(msg)
   end if
   if (.not.ex) then
!    Core WF file: new format XML
     filename_='corewf.xml';ex=.false.
     inquire(file=trim(filename_),iostat=ios,exist=ex)
     if (ios/=0) then
       write(msg,'(3a)') 'INQUIRE returns an error for file ',trim(filename_),'!'
       MSG_ERROR(msg)
     end if
     usexml=ex
     if (.not.ex) then
!      Core WF file: old format
       filename_='corewf.dat';ex=.false.
       inquire(file=trim(filename_),iostat=ios,exist=ex)
       if (ios/=0) then
         write(msg,'(3a)') 'INQUIRE returns an error for file ',trim(filename_),'!'
         MSG_ERROR(msg)
       end if
       oldformat=ex
       if (.not.ex) then
!        No core WF file found
         write(msg, '(3a)' )&
&         'Checks for existence of files corewf.abinit[.xml] or corewf.dat',ch10,&
&         'but INQUIRE finds file does not exist!'
         MSG_ERROR(msg)
       end if
     end if
   end if
 end if

!Core WF file is in new XML format
 if ((.not.oldformat).and.(usexml)) then
   call rdpawpsxml_core(energy_cor,filename_,lcor,ncor,nphicor,radmesh,phi_cor)
 endif

!Core WF file is in new (proprietary) format
 if ((.not.oldformat).and.(.not.usexml)) then
   unt = libpaw_get_free_unit()
   open(unt,file=trim(filename_),form='formatted',action="read")
   read(unt,*) ! skip title
   read(unt,*) ! skip relativism,method,nspinor,nsppol
   read(unt,*) ! skip zatom,zcore,pspdat
   read(unt,*) ! skip pspcod,pspxc,lmax
   read(unt,*) ! skip pspfmt,creatorID
   read(unt,*) nphicor
   read(unt,*) ! skip orbitals
   read(unt,*) nmesh
   LIBPAW_ALLOCATE(meshsz,(nmesh))
   LIBPAW_ALLOCATE(meshtp,(nmesh))
   LIBPAW_ALLOCATE(radstp,(nmesh))
   LIBPAW_ALLOCATE(logstp,(nmesh))
   do iln=1,nmesh
     r2=zero;read(unt,'(a80)') fline
     read(unit=fline,fmt=*,err=20,end=20) ib,i1,i2,r1,r2
     20 continue
     if (ib<=nmesh) then
       meshtp(ib)=i1;meshsz(ib)=i2
       radstp(ib)=r1;logstp(ib)=r2
     end if
   end do
   read(unt,*) ! skip rmax(core)
   LIBPAW_ALLOCATE(ncor,(nphicor))
   LIBPAW_ALLOCATE(lcor,(nphicor))
   LIBPAW_ALLOCATE(energy_cor,(nphicor))
   LIBPAW_ALLOCATE(phi_cor,(radmesh%mesh_size,nphicor))
   do iln=1,nphicor
     read(unt,*) ! skip comment
     read(unt,*) i1
     read(unt,*) ncor(iln),lcor(iln)
     read(unt,*) energy_cor(iln)
     energy_cor(iln)=energy_cor(iln)*half ! For consistency reasons (in the legacy coreWF format, energies are in Ry)
     LIBPAW_ALLOCATE(phitmp,(meshsz(i1)))
     read(unt,*) phitmp
     if ((radmesh%mesh_type/=meshtp(i1)) &
&     .or.(radmesh%rstep/=radstp(i1)) &
&     .or.(radmesh%lstep/=logstp(i1))) then
       call pawrad_init(tmpmesh,mesh_size=meshsz(i1),mesh_type=meshtp(i1),rstep=radstp(i1),lstep=logstp(i1))
       npts=radmesh%mesh_size
       if (tmpmesh%rmax<radmesh%rmax+tol8) npts=pawrad_ifromr(radmesh,tmpmesh%rmax)-1
       LIBPAW_ALLOCATE(work,(meshsz(i1)))
       call bound_deriv(phitmp,tmpmesh,meshsz(i1),r1,r2)
       call paw_spline(tmpmesh%rad,phitmp,meshsz(i1),r1,r2,work)
       call paw_splint(meshsz(i1),tmpmesh%rad,phitmp,work,npts,radmesh%rad(1:npts),phi_cor(1:npts,iln))
       if (npts<radmesh%mesh_size) phi_cor(npts+1:radmesh%mesh_size,iln)=zero
       LIBPAW_DEALLOCATE(work)
       call pawrad_free(tmpmesh)
     else
       npts=min(meshsz(i1),radmesh%mesh_size)
       phi_cor(1:npts,iln)=phitmp(1:npts)
       if (npts<radmesh%mesh_size) phi_cor(npts+1:radmesh%mesh_size,iln)=zero
     end if
     LIBPAW_DEALLOCATE(phitmp)
   end do
   LIBPAW_DEALLOCATE(meshsz)
   LIBPAW_DEALLOCATE(meshtp)
   LIBPAW_DEALLOCATE(radstp)
   LIBPAW_DEALLOCATE(logstp)
   close(unt)
 end if

!Core WF file is in old (proprietary) format
 if ((oldformat).and.(.not.usexml)) then
   unt = libpaw_get_free_unit()
   open(unt,file=trim(filename_),form='formatted',action="read")
   do while (dum/='atompaw ')
     read(unt,'(a8)') dum
   end do
   read(unt,'(2i4)') npts,nphicor
   LIBPAW_ALLOCATE(ncor,(nphicor))
   LIBPAW_ALLOCATE(lcor,(nphicor))
   LIBPAW_ALLOCATE(energy_cor,(nphicor))
   LIBPAW_ALLOCATE(phi_cor,(npts,nphicor))
   LIBPAW_ALLOCATE(rad,(npts))
   do iln=1,nphicor
     read(unt,'(a4,i4,a3,i4,a6,f15.7,a8,f15.7)') &
&     dum1,ncor(iln),dum2,lcor(iln),dum3,noccor,dum4,energy_cor(iln)
     energy_cor(iln)=energy_cor(iln)*half ! For consistency reasons (in the legacy coreWF format, energies are in Ry)

     do jln=1,npts
       read(unt,*) rad(jln),phi_cor(jln,iln)
     end do
     read(unt,*)
   end do
   LIBPAW_DEALLOCATE(rad)
   close(unt)
 end if

!Set an array 'a la' indlmn
 lmncmax=0
 do ib=1,nphicor
   il=lcor(ib)
   lmncmax=lmncmax+2*il+1
 end do
 LIBPAW_ALLOCATE(indlmn_core,(6,lmncmax))
 indlmn_core=0;ilmn=0;iln=0
 do ib=1,nphicor
   il=lcor(ib)
   iln=iln+1
   do ilm=1,2*il+1
     indlmn_core(1,ilmn+ilm)=il
     indlmn_core(2,ilmn+ilm)=ilm-(il+1)
     indlmn_core(3,ilmn+ilm)=1
     indlmn_core(4,ilmn+ilm)=il*il+ilm
     indlmn_core(5,ilmn+ilm)=iln
     indlmn_core(6,ilmn+ilm)=1
   end do
   ilmn=ilmn+2*il+1
 end do

end subroutine pawpsp_read_corewf
!!***


!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_rw_atompaw
!! NAME
!!  pawpsp_rw_atompaw
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_rw_atompaw(basis_size,filpsp,wvl)

!Arguments ------------------------------------
 integer,intent(in):: basis_size
 type(wvlpaw_type),intent(in)::wvl
 character(len=fnlen),intent(in)::filpsp
!arrays
 character(strlen) :: pspline
 character(len=fnlen)::fname

!Local variables-------------------------------
 integer :: ib,ii,ios,jj,step,iunt,ount
!arrays

! *************************************************************************
 iunt = libpaw_get_free_unit()
 ount = libpaw_get_free_unit()

 step=0
! Open psp file for reading
  open(unit=iunt,file=trim(filpsp),form='formatted',status='old',action="read")
! Open the file for writing
  write(fname,'(2a)') libpaw_basename(trim(filpsp)),".wvl"
  open(unit=ount,file=fname,form='formatted',status='unknown',action="write")

  read_loop: do
    if(step==0) then
      read(iunt,'(a)',IOSTAT=ios) pspline
      if ( ios /= 0 ) exit read_loop
      if(index(trim(pspline),'CORE_DENSITY')/=0 .and. &
&        index(trim(pspline),'PSEUDO_CORE_DENSITY')==0 .and. &
&        index(trim(pspline),'TCORE_DENSITY')==0 ) then
        step=1
      else
        write(ount,'(a)') trim(pspline)
      end if
    elseif(step==1) then
!Write Gaussian projectors:
      jj=0
      do ib=1,basis_size
        write(ount,'(a,i1,a)') "===== GAUSSIAN_TPROJECTOR ",ib,&
&       " =====   "
        write(ount,'(i5,1x,i5,1x,a)')wvl%pngau(ib),wvl%ptotgau, ":ngauss, total ngauss"
        write(ount,'(3(1x,es23.16))')(wvl%parg(:,ii),&
&        ii=jj+1,jj+wvl%pngau(ib))
        write(ount,'(3(1x,es23.16))')(wvl%pfac(:,ii),&
&        ii=jj+1,jj+wvl%pngau(ib))
        jj=jj+wvl%pngau(ib)
      end do
      write(ount,'(a)') trim(pspline)
      step=0
    end if
  end do read_loop

  close(iunt)
  close(ount)

end subroutine pawpsp_rw_atompaw
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_calc
!! NAME
!! pawpsp_calc
!!
!! FUNCTION
!! Performs tests and compute data related to pspcod=7 or 17 ("PAW pseudopotentials")
!!
!! INPUTS
!!  core_mesh<type(pawrad_type)>= radial mesh for the core density
!!  [coretau_mesh<type(pawrad_type)>]=radial mesh for the core kinetic energy density
!!  imainmesh= serial number of the main mesh
!!  ixc=exchange-correlation choice from main routine data file
!!  lnmax=max. number of (l,n) components over all type of psps
!!            angular momentum of nonlocal pseudopotential
!!  mqgrid_ff=dimension of q (or G) grid for nl form factors (array ffspl)
!!  mqgrid_vl=dimension of q (or G) grid for Vloc (array vlspl)
!!  ncore(core_mesh%mesh_size)= core density
!!  nmesh= number of radial meshes
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawxcdev=choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  pspversion= version of the atompaw code used to generate paw data.
!!  qgrid_ff(mqgrid_ff)=values of q on grid from 0 to qmax (bohr^-1) for nl form factors
!!  qgrid_vl(mqgrid_vl)=values of q on grid from 0 to qmax (bohr^-1) for Vloc
!!  radmesh(nmesh)<type(pawrad_type)>=paw radial meshes and related data
!!  tncore(core_mesh%mesh_size)= pseudo core density
!!  [tcoretau(coretau_mesh%mesh_size)]= pseudo core kinetic energy density
!!  tproj(tproj_mesh%mesh_size)= non-local projectors in real space
!!  tproj_mesh<type(pawrad_type)>= radial mesh for the projectors
!!  usexcnhat=0 if compensation charge density is not included in XC terms
!!            1 if compensation charge density is included in XC terms
!!  vale_mesh<type(pawrad_type)>= radial mesh for the valence density
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!  vlocopt= option for the local potential.(0=Vbare, 1=VH(tnzc) with hat in XC, 2=VH(tnzc) w/o hat in XC)
!!  vlocr(vloc_mesh%mesh_size)= local potential according to vlocopt.
!!  xclevel= XC functional level
!!  zion=nominal valence of atom as specified in psp file
!!  znucl=atomic number of atom as specified in input file to main routine
!!
!! OUTPUT
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$(hartree)
!!  ffspl(mqgrid_ff,2,lnmax)=form factor f_l(q) and second derivative
!!   from spline fit for each angular momentum and each projector;
!!  vlspl(mqgrid_vl,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!  xcccrc=XC core correction cutoff radius (bohr) from psp file
!!
!! SIDE EFFECTS
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  tnvale(vale_mesh%mesh_size)= pseudo valence density (+ nhat in output)
!!  vloc_mesh<type(pawrad_type)>= radial mesh for the local potential
!!
!! NOTES
!!
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_calc(core_mesh,epsatm,ffspl,imainmesh,ixc,lnmax,&
&          mmax,mqgrid_ff,mqgrid_vl,ncore,nmesh,pawrad,pawtab,pawxcdev,pspversion,&
&          qgrid_ff,qgrid_vl,radmesh,tncore,tnvale,tproj,tproj_mesh,usexcnhat,vale_mesh,&
&          vloc_mesh,vlocopt,vlocr,vlspl,xcccrc,xclevel,xc_denpos,zion,znucl,&
&          tcoretau,coretau_mesh) !optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: imainmesh,ixc,lnmax,mqgrid_ff,mqgrid_vl
 integer,intent(in) :: nmesh,pawxcdev,pspversion,usexcnhat,vlocopt
 integer,intent(in) ::mmax
 integer,intent(in) :: xclevel
 real(dp),intent(in) :: xc_denpos,zion,znucl
 real(dp),intent(out) :: epsatm,xcccrc
 type(pawrad_type),intent(in) :: core_mesh,tproj_mesh,vale_mesh
 type(pawrad_type),intent(in),optional :: coretau_mesh
 type(pawrad_type),intent(inout) ::pawrad,vloc_mesh
 type(pawtab_type),intent(inout) :: pawtab
!arrays
 real(dp),intent(in) :: ncore(core_mesh%mesh_size),tncore(core_mesh%mesh_size)
 real(dp),intent(in),optional :: tcoretau(:)
 real(dp),intent(in) :: qgrid_vl(mqgrid_vl),qgrid_ff(mqgrid_ff)
 real(dp),intent(inout) :: ffspl(mqgrid_ff,2,lnmax)
 real(dp),intent(out) :: vlspl(mqgrid_vl,2)
 real(dp),intent(inout) :: tnvale(vale_mesh%mesh_size*pawtab%has_tvale)
 real(dp),intent(inout) :: tproj(tproj_mesh%mesh_size,pawtab%basis_size)
 real(dp),intent(inout) :: vlocr(vloc_mesh%mesh_size)
 type(pawrad_type),intent(in) :: radmesh(nmesh)

!Local variables ------------------------------
!scalars
 integer,parameter :: reduced_mshsz=2501
 integer :: ib,il,ilm,ilmn,iln,ir,isnotzero,itest
 integer :: j0lmn,jlm,jlmn,jln,klmn,msz,msz1,msz_tmp,mst_tmp,nspden,usekden
 logical :: has_dij0,non_magnetic_xc,reduced_ncor,reduced_taucor,reduced_nval,reduced_vloc,testval
 real(dp),parameter :: reduced_rstep=0.00025_dp,rm_vloc=20.0_dp
 real(dp) :: d2nvdq0,intg,intvh,lstep_tmp,qcore,qq,rstep_tmp,yp1,ypn
 character(len=500) :: msg
 type(pawang_type) :: pawang_tmp
 type(pawrad_type) :: rcore_mesh,rcoretau_mesh,rvale_mesh,rvloc_mesh,tproj_mesh_new
!arrays
 real(dp) :: tmp_qgrid(1),tmp_q2vq(1)
 real(dp),allocatable :: ncorwk(:),nhat(:),nhatwk(:),nwk(:),r2k(:)
 real(dp),allocatable :: rtncor(:),rttaucor(:),rtnval(:),rvlocr(:)
 real(dp),allocatable :: vbare(:),vh(:),vhnzc(:),vxc1(:),vxc2(:)
 real(dp),allocatable,target :: work1(:),work2(:),work3(:)
 real(dp),pointer :: tmp1(:),tmp2(:)
 logical :: tmp_lmselect(1)

! *************************************************************************

!==========================================================
!Perfom tests on meshes

!Are radial meshes for Phi and Vloc compatibles ?
! if (vloc_mesh%rmax<pawrad%rmax) then
!   write(msg, '(a,a,a)' )&
!&   'Rmax for Vloc < Rmax !',ch10,&
!&   'Action : check your pseudopotential (increase Vloc meshSize).'
!   MSG_ERROR(msg)
! end if

!Check optional arguments
 usekden=merge(0,1,pawtab%has_coretau==0)
 if (present(tcoretau)) then
   if (usekden>=1) then
     if (.not.(present(coretau_mesh))) then
       msg='tcoretau present but not coretau_mesh!'
       MSG_BUG(msg)
     end if
     if (size(tcoretau)>coretau_mesh%mesh_size) then
       msg='wrong size for tcoretau!'
       MSG_BUG(msg)
     end if
     if (coretau_mesh%mesh_size<pawtab%mesh_size) then
       write(msg, '(a,a,a,a,a)' )&
&       'Mesh size for core kinetic energy density must be equal or larger',ch10,&
&       'than mesh size for PAW augmentation regions !',ch10,&
&       'Action : check your pseudopotential (increase TAUcore meshSize).'
       MSG_ERROR(msg)
     end if
   end if
 end if

!Are mmax and mesh_size for partial waves compatibles ?
 if (mmax/=pawtab%partialwave_mesh_size) then
   write(msg, '(a,a,a)' )&
&   'mmax /= phi_mesh_size in psp file !',ch10,&
&   'Action: check your pseudopotential file.'
   MSG_ERROR(msg)
 end if

!Are radial meshes for (t)Ncore / tTAU and Phi compatibles ?
 if (usekden>=1.and.core_mesh%mesh_size<pawtab%mesh_size) then
   write(msg, '(a,a,a,a,a)' )&
&   'Mesh size for core density must be equal or larger',ch10,&
&   'than mesh size for PAW augmentation regions !',ch10,&
&   'Action : check your pseudopotential (increase Ncore meshSize).'
   MSG_ERROR(msg)
 end if

!Are radial meshes for (t)Nvale and Phi compatibles ?
 if ((pawtab%has_tvale==1).and.(vale_mesh%rmax<pawrad%rmax)) then
   write(msg, '(a,a,a)' )&
&   'Rmax for tNvale < Rmax for Phi !',ch10,&
&   'Action : check your pseudopotential (increase tNvale meshSize).'
   MSG_ERROR(msg)
 end if

!Is PAW radius included inside radial mesh ?
 if (pawtab%rpaw>pawrad%rmax+tol8) then
   write(msg, '(a,a,a)' )&
&   'Radius of PAW sphere is outside the radial mesh !',ch10,&
&   'Action: check your pseudopotential file.'
   MSG_ERROR(msg)
 end if

!Max. radius of mesh for Vloc has to be "small" in order to avoid numeric noise ?
 if (vloc_mesh%rmax>rm_vloc) then
   msz_tmp=pawrad_ifromr(vloc_mesh,rm_vloc);mst_tmp=vloc_mesh%mesh_type
   rstep_tmp=vloc_mesh%rstep;lstep_tmp=vloc_mesh%lstep
   call pawrad_free(vloc_mesh)
   call pawrad_init(vloc_mesh,mesh_size=msz_tmp,mesh_type=mst_tmp,&
&   rstep=rstep_tmp,lstep=lstep_tmp,r_for_intg=rm_vloc)
   write(msg, '(a,i4,a)' ) ' Mesh size for Vloc has been set to ', &
&    vloc_mesh%mesh_size,' to avoid numerical noise.'
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

!This test has been disable... MT 2006-25-10
!For Simpson rule, it is better to have odd mesh sizes
!itest=0
!do imsh=1,nmesh
!if (mod(radmesh(imsh)%mesh_size,2)==0.and.radmesh(imsh)%mesh_type==1) itest=1
!end do
!if (itest==1) then
!  write(msg, '(5a)' ) &
!&   'Regular radial meshes should have odd number of points ',ch10,&
!&   'for better accuracy of integration sheme (Simpson rule).',ch10,&
!&   'Althought it''s not compulsory, you should change mesh sizes in psp file.'
!  MSG_WARNING(msg)
!end if

!Test the compatibilty between Rpaw and mesh for (t)Phi
 if (pspversion>=3) then
   itest=pawrad_ifromr(radmesh(imainmesh),pawtab%rpaw)
!  This test has been disable... MT 2015-02-12
!   if (itest+2>radmesh(imainmesh)%mesh_size) then
!     write(msg, '(9a)' ) &
!&     'Atomic data could produce inaccurate results:',ch10,&
!&     'Wavefunctions and pseudo-wavefunctions should',ch10,&
!&     'be given on a radial mesh larger than the PAW',ch10,&
!&     'spheres (at least 2 additional points) !',ch10,&
!&     'Action: check your pseudopotential file.'
!     MSG_WARNING(msg)
!   end if
   if (abs(pawtab%rpaw-radmesh(imainmesh)%rad(itest))<tol8) itest=itest-1
   ib=0;isnotzero=0
   do while ((isnotzero==0).and.(ib<pawtab%basis_size))
     ib=ib+1;ir=itest
     do while ((isnotzero==0).and.(ir<pawtab%mesh_size))
       ir=ir+1;if (abs(pawtab%phi(ir,ib)-pawtab%tphi(ir,ib))>tol8) isnotzero=1
     end do
   end do
   if (isnotzero>0) then
     write(msg, '(7a)' )&
&     'Atomic data are inconsistent:',ch10,&
&     'For r>=r_paw, pseudo wavefunctions are not',ch10,&
&     'equal to wave functions (Phi(r)/=tPhi(r)) !',ch10,&
&     'Action: check your pseudopotential file.'
     MSG_ERROR(msg)
   end if
 else
!  For compatibility reasons set PAW radius at the end of mesh (older versions)
   if (pawtab%rpaw/=pawrad%rmax) then
     msz_tmp=pawrad%mesh_size;mst_tmp=pawrad%mesh_type
     rstep_tmp=pawrad%rstep;lstep_tmp=pawrad%lstep
     call pawrad_free(pawrad)
     call pawrad_init(pawrad,mesh_size=msz_tmp,mesh_type=mst_tmp,rstep=rstep_tmp,lstep=lstep_tmp)
     pawtab%rpaw=pawrad%rmax
   end if
 end if
!If Vloc is a "Vbare" potential, it has to be localized inside PAW spheres
 if (vlocopt==0.and.(vloc_mesh%rmax>pawtab%rpaw+tol10)) then
   if(vlocr(pawrad_ifromr(vloc_mesh,pawtab%rpaw))>tol10) then
     write(msg, '(7a)' )&
&     'Atomic data are inconsistent:',ch10,&
&     'Local potential is a "Vbare" potential',ch10,&
&     'and is not localized inside PAW sphere !',ch10,&
&     'Vbare is set to zero if r>rpaw.'
     MSG_WARNING(msg)
     do ir=pawrad_ifromr(vloc_mesh,pawtab%rpaw),vloc_mesh%mesh_size
       vlocr(ir)=zero
     end do
   end if
 end if

!==========================================================
!Initializations

 has_dij0=(allocated(pawtab%dij0))

!Allocate/initialize some dummy variables
 tmp_lmselect(1)=.true.
 non_magnetic_xc=.false.
 if (pawxcdev==0) then
   pawang_tmp%l_size_max=1;pawang_tmp%angl_size=1;pawang_tmp%ylm_size=1
   pawang_tmp%use_ls_ylm=0;pawang_tmp%gnt_option=0;pawang_tmp%ngnt=0;pawang_tmp%nsym=0
   LIBPAW_ALLOCATE(pawang_tmp%angwgth,(1))
   pawang_tmp%angwgth(1)=one
   LIBPAW_ALLOCATE(pawang_tmp%anginit,(3,1))
   pawang_tmp%anginit(1,1)=one
   pawang_tmp%anginit(2:3,1)=zero
   LIBPAW_ALLOCATE(pawang_tmp%ylmr,(1,1))
   pawang_tmp%ylmr(1,1)=1._dp/sqrt(four_pi)
   LIBPAW_ALLOCATE(pawang_tmp%ylmrgr,(9,1,1))
   pawang_tmp%ylmrgr(1:9,1,1)=zero
 end if

!==========================================================
!Compute ffspl(q) (and derivatives)

 ffspl=zero
 if (mqgrid_ff>0) then
   call pawpsp_nl(ffspl,pawtab%indlmn,pawtab%lmn_size,lnmax,mqgrid_ff,qgrid_ff,&
&                 tproj_mesh,tproj)
 end if

!==========================================================
!Compute eventually compensation charge radius (i.e. radius for shape functions)

 if (pawtab%shape_type>0.and.pawtab%rshp<1.d-8) then
   pawtab%rshp=pawtab%rpaw
 else if (pawtab%shape_type==-1) then
   ir=pawrad_ifromr(radmesh(imainmesh),pawtab%rpaw)+1;isnotzero=0
   do while ((isnotzero==0).and.(ir>1))
     ir=ir-1;il=0
     do while ((isnotzero==0).and.(il<pawtab%l_size))
       il=il+1;if (pawtab%shapefunc(ir,il)>tol16) isnotzero=1
     end do
   end do
   ir=min(ir+1,pawrad_ifromr(radmesh(imainmesh),pawtab%rpaw))
   pawtab%rshp=radmesh(imainmesh)%rad(ir)
   do il=1,pawtab%l_size
     if (pawtab%shapefunc(ir,il)>tol6) then
       write(msg, '(a,a,a)' )&
&       'Shape function is not zero at PAW radius !',ch10,&
&       'Action: check your pseudopotential file.'
       MSG_ERROR(msg)
     end if
   end do
 end if

!==========================================================
!Compute compensation charge density (nhat)
!Add it to pseudo valence density

 if (pawtab%has_tvale==1) then
   msz=vale_mesh%mesh_size
   LIBPAW_ALLOCATE(nhat,(msz))
!  A-Has to compute norm of nhat (Int[n-tild_n])
   testval=(abs(tnvale(msz))<tol9)
!  A1-If tnvale is not given with enough points,
!  try to compute it from rhoij0 and tphi
   if (.not.testval) then
     msz1=pawtab%mesh_size
!    Compute n and tild_n from phi and tphi
     LIBPAW_ALLOCATE(work1,(msz1))
     LIBPAW_ALLOCATE(work2,(msz1))
     work1=zero
     work2=zero
     do jlmn=1,pawtab%lmn_size
       j0lmn=jlmn*(jlmn-1)/2;jln=pawtab%indlmn(5,jlmn)
       do ilmn=1,jlmn
         klmn=j0lmn+ilmn;iln=pawtab%indlmn(5,ilmn)
         yp1=two;if (ilmn==jlmn) yp1=one
         work1(1:msz1)=work1(1:msz1)+yp1*pawtab%rhoij0(klmn) &
&         *pawtab% phi(1:msz1,iln)*pawtab% phi(1:msz1,jln)
         work2(1:msz1)=work2(1:msz1)+yp1*pawtab%rhoij0(klmn) &
&         *pawtab%tphi(1:msz1,iln)*pawtab%tphi(1:msz1,jln)
       end do
     end do
!    Spline tnvale onto pawrad if needed
     LIBPAW_ALLOCATE(nwk,(msz1))
     if ((vale_mesh%mesh_type/=pawrad%mesh_type).or.(vale_mesh%rstep/=pawrad%rstep).or.&
&     (vale_mesh%lstep/=pawrad%lstep)) then
       LIBPAW_ALLOCATE(work3,(vale_mesh%mesh_size))
       call bound_deriv(tnvale(1:vale_mesh%mesh_size),vale_mesh,vale_mesh%mesh_size,yp1,ypn)
       call paw_spline(vale_mesh%rad,tnvale,vale_mesh%mesh_size,yp1,ypn,work3)
       call paw_splint(vale_mesh%mesh_size,vale_mesh%rad,tnvale,work3,msz1,pawrad%rad(1:msz1),nwk(1:msz1))
       LIBPAW_DEALLOCATE(work3)
     else
       nwk(1:msz1)=tnvale(1:msz1)
     end if
!    Compare tild_n and tnvale (inside aug. region)
     if (maxval(abs((nwk(1:msz1)*four_pi*pawrad%rad(1:msz1)**2)-work2(1:msz1)))<tol6) then
!      If equality then compute Int[n-tild_n]
       work1=work1-work2
       call simp_gen(qq,work1,pawrad)
       qq=qq/four_pi
     else
!      If not equality, will use tnvale
       testval=.true.
       write(msg, '(3a)' ) &
&       'Valence density is not given with enough points',ch10,&
&       'in psp file. Some charge estimations will be coarse.'
       MSG_WARNING(msg)
     end if
     LIBPAW_DEALLOCATE(nwk)
     LIBPAW_DEALLOCATE(work1)
     LIBPAW_DEALLOCATE(work2)
   end if
!  A2-If tnvale is given with enough points, use it
   if (testval) then
     nhat(1:msz)=tnvale(1:msz)*vale_mesh%rad(1:msz)**2
     call simp_gen(qq,nhat,vale_mesh)
     qq=zion/four_pi-qq
   end if
!  B-Compute nhat and add it to pseudo valence density
   call atompaw_shpfun(0,vale_mesh,intg,pawtab,nhat)
   nhat(1:msz)=qq*nhat(1:msz)
   tnvale(1:msz)=tnvale(1:msz)+nhat(1:msz)
 end if

!==========================================================
!If Vloc potential is in "Vbare" format, translate it into VH(tnzc) format

 if (vlocopt==0) then
   write(msg,'(a)') ' Local potential is in "Vbare" format... '
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
   msz=core_mesh%mesh_size
   LIBPAW_ALLOCATE(r2k,(msz))
   call atompaw_shpfun(0,core_mesh,intg,pawtab,r2k)
   r2k(1:msz)=r2k(1:msz)*core_mesh%rad(1:msz)**2
!  Compute VH[4pi.r2.n(r)=4pi.r2.tncore(r)+(Qcore-Z).r2.k(r)]
   LIBPAW_ALLOCATE(nwk,(core_mesh%mesh_size))
   LIBPAW_ALLOCATE(vh,(core_mesh%mesh_size))
   if (core_mesh%mesh_type==5) then
     nwk(1:msz)=tncore(1:msz)*four_pi*core_mesh%rad(1:msz)**2
     call simp_gen(qcore,nwk,core_mesh)
     qcore=znucl-zion-qcore
   else
     nwk(1:msz)=(ncore(1:msz)-tncore(1:msz))*four_pi*core_mesh%rad(1:msz)**2
     ib=1
     do ir=msz,2,-1
       if(abs(nwk(ir))<tol14)ib=ir
     end do
     call simp_gen(qcore,nwk,core_mesh,r_for_intg=core_mesh%rad(ib))
     nwk(1:msz)=tncore(1:msz)*four_pi*core_mesh%rad(1:msz)**2
   end if
   nwk(1:msz)=nwk(1:msz)+r2k(1:msz)*(qcore-znucl)
   call poisson(nwk,0,core_mesh,vh)
   vh(2:msz)=vh(2:msz)/core_mesh%rad(2:msz)
   call pawrad_deducer0(vh,msz,core_mesh)

   LIBPAW_DEALLOCATE(nwk)
!  Eventually spline Vbare
   LIBPAW_ALLOCATE(vbare,(core_mesh%mesh_size))
   if ((core_mesh%mesh_type/=vloc_mesh%mesh_type).or.&
&   (core_mesh%rstep    /=vloc_mesh%rstep)    .or.&
&   (core_mesh%lstep    /=vloc_mesh%lstep)) then
     msz=core_mesh%mesh_size;if (vloc_mesh%rmax<core_mesh%rmax) msz=pawrad_ifromr(core_mesh,vloc_mesh%rmax)
     call bound_deriv(vlocr(1:vloc_mesh%mesh_size),vloc_mesh,vloc_mesh%mesh_size,yp1,ypn)
     LIBPAW_ALLOCATE(work1,(vloc_mesh%mesh_size))
     LIBPAW_ALLOCATE(work2,(vloc_mesh%mesh_size))
     call paw_spline(vloc_mesh%rad,vlocr,vloc_mesh%mesh_size,yp1,ypn,work1)
     call paw_splint(vloc_mesh%mesh_size,vloc_mesh%rad,vlocr,work1,msz,core_mesh%rad(1:msz),vbare)
     LIBPAW_DEALLOCATE(work1)
     LIBPAW_DEALLOCATE(work2)
   else
     msz=min(core_mesh%mesh_size,vloc_mesh%mesh_size)
     vbare(1:msz)=vlocr(1:msz)
   end if
!  Build VH(tnzc) from Vbare
   vlocr(1:msz)=vbare(1:msz)+vh(1:msz)
   if(vloc_mesh%mesh_size>msz)then
     vlocr(msz+1:vloc_mesh%mesh_size)=vh(msz)*vloc_mesh%rad(msz)/vloc_mesh%rad(msz+1:vloc_mesh%mesh_size)
   end if
   LIBPAW_DEALLOCATE(vbare)
   LIBPAW_DEALLOCATE(vh)

!  Compute <tPhi_i|VH(tnzc)|tPhi_j> and int[VH(tnzc)*Qijhat(r)dr] parts of Dij0
!  Note: it is possible as core_mesh and radmesh(imainmesh) have the same steps
   if (has_dij0) then
     msz=radmesh(imainmesh)%mesh_size
     LIBPAW_ALLOCATE(work1,(msz))
     work1(1:msz)=vlocr(1:msz)*r2k(1:msz)
     call simp_gen(intvh,work1,radmesh(imainmesh))
     do jlmn=1,pawtab%lmn_size
       j0lmn=jlmn*(jlmn-1)/2;jlm=pawtab%indlmn(4,jlmn);jln=pawtab%indlmn(5,jlmn)
       do ilmn=1,jlmn
         klmn=j0lmn+ilmn;ilm=pawtab%indlmn(4,ilmn);iln=pawtab%indlmn(5,ilmn)
         if (jlm==ilm) then
           work1(1:msz)=pawtab%tphi(1:msz,iln)*pawtab%tphi(1:msz,jln)*(vlocr(1:msz)-intvh) &
&           -pawtab%phi (1:msz,iln)*pawtab%phi (1:msz,jln)*intvh
           call simp_gen(intg,work1,radmesh(imainmesh))
           pawtab%dij0(klmn)=pawtab%dij0(klmn)+intg
         end if
       end do
     end do
     LIBPAW_DEALLOCATE(work1)
   end if
   LIBPAW_DEALLOCATE(r2k)
 end if

!==========================================================
!If usexcnhat in psp file is different from usexcnhat chosen
!by user, convert VH(tnzc) and Dij0

 if (pawtab%usexcnhat==-1) then
   pawtab%usexcnhat=usexcnhat
 else if (usexcnhat/=pawtab%usexcnhat) then
   if (pawtab%has_tvale==0) then
     write(msg, '(5a)' ) &
&     'It is only possible to modify the use of compensation charge density',ch10,&
&     'for a file format containing the pseudo valence density (format>=paw4 or XML)!',ch10,&
&     'Action: use usexcnhat=-1 in input file or change psp file format.'
     MSG_WARNING(msg)
   else if (usekden>=1) then
     write(msg, '(5a)' ) &
&     'It is not possible to modify the use of compensation charge density',ch10,&
&     'within the metaGGA XC functional (need valence kinetic density)!',ch10,&
&     'Action: use usexcnhat=-1 in input file or change psp file format.'
     MSG_WARNING(msg)
   else
     msz=vloc_mesh%mesh_size
!    Retrieve tvale and nhat onto vloc mesh
     LIBPAW_ALLOCATE(nwk,(msz))
     LIBPAW_ALLOCATE(ncorwk,(msz))
     LIBPAW_ALLOCATE(nhatwk,(msz))
     nwk=zero;ncorwk=zero;nhatwk=zero
     if ((core_mesh%mesh_type/=vloc_mesh%mesh_type).or.&
&        (core_mesh%rstep    /=vloc_mesh%rstep)    .or.&
&        (core_mesh%lstep    /=vloc_mesh%lstep)) then
       LIBPAW_ALLOCATE(work1,(core_mesh%mesh_size))
       msz1=msz;if (core_mesh%rmax<vloc_mesh%rmax) msz1=pawrad_ifromr(vloc_mesh,core_mesh%rmax)
       call bound_deriv(tncore(1:core_mesh%mesh_size),core_mesh,core_mesh%mesh_size,yp1,ypn)
       call paw_spline(core_mesh%rad,tncore,core_mesh%mesh_size,yp1,ypn,work1)
       call paw_splint(core_mesh%mesh_size,core_mesh%rad,tncore,work1,msz1,vloc_mesh%rad(1:msz1),ncorwk(1:msz1))
       LIBPAW_DEALLOCATE(work1)
     else
       msz1=min(core_mesh%mesh_size,msz)
       ncorwk(1:msz1)=tncore(1:msz1)
     end if
     if ((vale_mesh%mesh_type/=vloc_mesh%mesh_type).or.&
&        (vale_mesh%rstep    /=vloc_mesh%rstep)    .or.&
&        (vale_mesh%lstep    /=vloc_mesh%lstep)) then
       LIBPAW_ALLOCATE(work1,(vale_mesh%mesh_size))
       msz1=msz;if (vale_mesh%rmax<vloc_mesh%rmax) msz1=pawrad_ifromr(vloc_mesh,vale_mesh%rmax)
       call bound_deriv(tnvale(1:vale_mesh%mesh_size),vale_mesh,vale_mesh%mesh_size,yp1,ypn)
       call paw_spline(vale_mesh%rad,tnvale,vale_mesh%mesh_size,yp1,ypn,work1)
       call paw_splint(vale_mesh%mesh_size,vale_mesh%rad,tnvale,work1,msz1,vloc_mesh%rad(1:msz1),nwk(1:msz1))
       call bound_deriv(nhat(1:vale_mesh%mesh_size),vale_mesh,vale_mesh%mesh_size,yp1,ypn)
       call paw_spline(vale_mesh%rad,nhat,vale_mesh%mesh_size,yp1,ypn,work1)
       call paw_splint(vale_mesh%mesh_size,vale_mesh%rad,nhat,work1,msz1,vloc_mesh%rad(1:msz1),nhatwk(1:msz1))
       LIBPAW_DEALLOCATE(work1)
     else
       msz1=min(vale_mesh%mesh_size,msz)
       nwk   (1:msz1)=tnvale(1:msz1)
       nhatwk(1:msz1)=nhat  (1:msz1)
     end if

     nwk=nwk-nhatwk
     nwk=sqrt(four_pi)*nwk;nhatwk=sqrt(four_pi)*nhatwk ! 0th-order moment of densities

!    Compute Vxc without nhat (vxc1) and with nhat (vxc2)
     nspden=1
#if defined LIBPAW_HAVE_LIBXC
     if (ixc<0) nspden=libxc_functionals_nspin()
#endif
     if (ixc<0) then
       LIBPAW_ALLOCATE(vxc1,(msz*nspden))
       LIBPAW_ALLOCATE(vxc2,(msz*nspden))
       LIBPAW_ALLOCATE(work1,(msz))
       LIBPAW_ALLOCATE(work2,(msz*nspden))
       LIBPAW_ALLOCATE(work3,(msz*nspden))
       tmp1 => work1
       work2(1:msz)=nwk
       work3(1:msz)=nhatwk
       if (nspden==2) work2(msz+1:2*msz)=half*nwk
       if (nspden==2) work3(msz+1:2*msz)=half*nhatwk
       if (pawxcdev/=0) then
         call pawxcm(ncorwk,yp1,ypn,0,ixc,work1,1,tmp_lmselect,work3,0,non_magnetic_xc,msz,nspden,5,&
&         pawang_tmp,vloc_mesh,pawxcdev,work2,pawtab%usetcore,0,vxc1,xclevel,xc_denpos)
         call pawxcm(ncorwk,yp1,ypn,0,ixc,work1,1,tmp_lmselect,work3,0,non_magnetic_xc,msz,nspden,5,&
&         pawang_tmp,vloc_mesh,pawxcdev,work2,pawtab%usetcore,2,vxc2,xclevel,xc_denpos)
         vxc1=vxc1/sqrt(four_pi);vxc2=vxc2/sqrt(four_pi) ! Deduce Vxc from its first moment
       else
         call pawxc(ncorwk,yp1,ypn,ixc,work1,tmp1,1,tmp_lmselect,work3,0,0,non_magnetic_xc,msz,nspden,5,&
&         pawang_tmp,vloc_mesh,work2,pawtab%usetcore,0,vxc1,xclevel,xc_denpos)
         call pawxc(ncorwk,yp1,ypn,ixc,work1,tmp1,1,tmp_lmselect,work3,0,0,non_magnetic_xc,msz,nspden,5,&
&         pawang_tmp,vloc_mesh,work2,pawtab%usetcore,2,vxc2,xclevel,xc_denpos)
       end if
       LIBPAW_DEALLOCATE(nwk)
       LIBPAW_DEALLOCATE(ncorwk)
       LIBPAW_DEALLOCATE(nhatwk)
       LIBPAW_DEALLOCATE(work1)
       LIBPAW_DEALLOCATE(work2)
       LIBPAW_DEALLOCATE(work3)
     else
       LIBPAW_ALLOCATE(vxc1,(msz))
       LIBPAW_ALLOCATE(vxc2,(msz))
       LIBPAW_ALLOCATE(work1,(msz))
       tmp1 => work1
       if (pawxcdev/=0) then
         call pawxcm(ncorwk,yp1,ypn,0,ixc,work1,1,tmp_lmselect,nhatwk,0,non_magnetic_xc,msz,1,5,&
&         pawang_tmp,vloc_mesh,pawxcdev,nwk,pawtab%usetcore,0,vxc1,xclevel,xc_denpos)
         call pawxcm(ncorwk,yp1,ypn,0,ixc,work1,1,tmp_lmselect,nhatwk,0,non_magnetic_xc,msz,1,5,&
&         pawang_tmp,vloc_mesh,pawxcdev,nwk,pawtab%usetcore,2,vxc2,xclevel,xc_denpos)
         vxc1=vxc1/sqrt(four_pi);vxc2=vxc2/sqrt(four_pi) ! Deduce Vxc from its first moment
       else
         call pawxc(ncorwk,yp1,ypn,ixc,work1,tmp1,1,tmp_lmselect,nhatwk,0,0,non_magnetic_xc,msz,1,5,&
&         pawang_tmp,vloc_mesh,nwk,pawtab%usetcore,0,vxc1,xclevel,xc_denpos)
         call pawxc(ncorwk,yp1,ypn,ixc,work1,tmp1,1,tmp_lmselect,nhatwk,0,0,non_magnetic_xc,msz,1,5,&
&         pawang_tmp,vloc_mesh,nwk,pawtab%usetcore,2,vxc2,xclevel,xc_denpos)
       end if
       LIBPAW_DEALLOCATE(nwk)
       LIBPAW_DEALLOCATE(ncorwk)
       LIBPAW_DEALLOCATE(nhatwk)
       LIBPAW_DEALLOCATE(work1)
     endif
!    Compute difference of XC potentials
     if (usexcnhat==0.and.pawtab%usexcnhat/=0)  vxc1(1:msz)=vxc2(1:msz)-vxc1(1:msz)
     if (usexcnhat/=0.and.pawtab%usexcnhat==0)  vxc1(1:msz)=vxc1(1:msz)-vxc2(1:msz)
!    Modify VH(tnzc)
     vlocr(1:msz)=vlocr(1:msz)-vxc1(1:msz)
     if (has_dij0) then
!      Modify  Dij0
       LIBPAW_ALLOCATE(work2,(pawtab%lmn2_size))
       call atompaw_kij(pawtab%indlmn,work2,pawtab%lmn_size,ncore,0,0,pawtab,pawrad,&
&                       core_mesh,vloc_mesh,vxc1(1:msz),znucl)
       pawtab%dij0=work2
       LIBPAW_DEALLOCATE(work2)
     end if
     LIBPAW_DEALLOCATE(vxc1)
     LIBPAW_DEALLOCATE(vxc2)
   end if ! has_tvale/=0
 end if
 if (pawtab%usexcnhat==0) then
   write(msg,'(a)') &
&   ' Compensation charge density is not taken into account in XC energy/potential'
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 end if
 if (pawtab%usexcnhat==1) then
   write(msg,'(a)') &
&   ' Compensation charge density is taken into account in XC energy/potential'
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 end if

!==========================================================
! Calculate the coefficient beta = \int { vH[nZc](r) - vloc(r) } 4pi r^2 dr
!
 LIBPAW_ALLOCATE(vhnzc,(core_mesh%mesh_size))
 LIBPAW_ALLOCATE(nwk,(core_mesh%mesh_size))
! get vH[nZc]
 call atompaw_vhnzc(ncore,core_mesh,vhnzc,znucl)

!Transpose vlocr mesh into core mesh
 nwk(:)=zero
 if ((core_mesh%mesh_type/=vloc_mesh%mesh_type).or.&
& (core_mesh%rstep    /=vloc_mesh%rstep)    .or.&
& (core_mesh%lstep    /=vloc_mesh%lstep)) then
   msz=core_mesh%mesh_size;if (vloc_mesh%rmax<core_mesh%rmax) msz=pawrad_ifromr(core_mesh,vloc_mesh%rmax)
   call bound_deriv(vlocr(1:vloc_mesh%mesh_size),vloc_mesh,vloc_mesh%mesh_size,yp1,ypn)
   LIBPAW_ALLOCATE(work1,(vloc_mesh%mesh_size))
   LIBPAW_ALLOCATE(work2,(vloc_mesh%mesh_size))
   call paw_spline(vloc_mesh%rad,vlocr,vloc_mesh%mesh_size,yp1,ypn,work1)
   call paw_splint(vloc_mesh%mesh_size,vloc_mesh%rad,vlocr,work1,msz,core_mesh%rad(1:msz),nwk)
   LIBPAW_DEALLOCATE(work1)
   LIBPAW_DEALLOCATE(work2)
 else
   msz=min(core_mesh%mesh_size,vloc_mesh%mesh_size)
   nwk(1:msz)=vlocr(1:msz)
 end if

!Difference
 nwk(1:msz)=vhnzc(1:msz)-nwk(1:msz)
 if (msz<core_mesh%mesh_size) nwk(msz+1:core_mesh%mesh_size)=zero

!Perform the spherical integration
 nwk(1:msz)=nwk(1:msz)*four_pi*core_mesh%rad(1:msz)**2

 call simp_gen(pawtab%beta,nwk,core_mesh)

 LIBPAW_DEALLOCATE(vhnzc)
 LIBPAW_DEALLOCATE(nwk)

 write(msg,'(a,e18.6)') &
&  ' beta integral value: ',pawtab%beta
 call wrtout(std_out,msg,'COLL')


!==========================================================
!Try to optimize CPU time:
!If Vloc mesh size is big, spline Vloc into a smaller log. mesh

 reduced_vloc=(vloc_mesh%mesh_size>int(reduced_mshsz))
 if (reduced_vloc) then
   msz=vloc_mesh%mesh_size
   lstep_tmp=log(0.9999999_dp*vloc_mesh%rmax/reduced_rstep)/dble(reduced_mshsz-2)
   call pawrad_init(rvloc_mesh,mesh_size=reduced_mshsz,mesh_type=3,&
&   rstep=reduced_rstep,lstep=lstep_tmp)
   LIBPAW_ALLOCATE(rvlocr,(reduced_mshsz))
   call bound_deriv(vlocr(1:msz),vloc_mesh,msz,yp1,ypn)
   LIBPAW_ALLOCATE(work1,(msz))
   LIBPAW_ALLOCATE(work2,(msz))
   LIBPAW_ALLOCATE(work3,(msz))
   work3(1:msz)=vloc_mesh%rad(1:msz)
   call paw_spline(work3,vlocr,msz,yp1,ypn,work1)
   call paw_splint(msz,work3,vlocr,work1,reduced_mshsz,rvloc_mesh%rad,rvlocr)
   LIBPAW_DEALLOCATE(work1)
   LIBPAW_DEALLOCATE(work2)
   LIBPAW_DEALLOCATE(work3)
 end if

!Keep VH(tnZc) eventually in memory
 if (pawtab%has_vhtnzc==1) then
   LIBPAW_ALLOCATE(pawtab%vhtnzc,(pawtab%mesh_size))
   if ((reduced_vloc).and.(rvloc_mesh%mesh_type==pawrad%mesh_type)&
&   .and.(rvloc_mesh%rstep==pawrad%rstep).and.(rvloc_mesh%lstep==pawrad%lstep)) then
     pawtab%vhtnzc(1:pawtab%mesh_size)=rvlocr(1:pawtab%mesh_size)
     pawtab%has_vhtnzc=2
   else if ((vloc_mesh%mesh_type==pawrad%mesh_type)&
&     .and.(vloc_mesh%rstep==pawrad%rstep).and.(vloc_mesh%lstep==pawrad%lstep)) then
     pawtab%vhtnzc(1:pawtab%mesh_size)=vlocr(1:pawtab%mesh_size)
     pawtab%has_vhtnzc=2
   else
     msg = 'Vloc mesh is not right !'
     MSG_ERROR(msg)
   end if
 end if

!==========================================================
!Try to optimize CPU time:
!If ncore mesh size is big, spline tncore into a smaller log. mesh

 reduced_ncor=(core_mesh%mesh_size>int(reduced_mshsz)).and.(pawtab%usetcore/=0)
 if (reduced_ncor) then
   msz=core_mesh%mesh_size
   lstep_tmp=log(0.9999999_dp*core_mesh%rmax/reduced_rstep)/dble(reduced_mshsz-2)
   call pawrad_init(rcore_mesh,mesh_size=reduced_mshsz,mesh_type=3,&
&   rstep=reduced_rstep,lstep=lstep_tmp)
   LIBPAW_ALLOCATE(rtncor,(reduced_mshsz))
   call bound_deriv(tncore(1:msz),core_mesh,msz,yp1,ypn)
   LIBPAW_ALLOCATE(work1,(msz))
   LIBPAW_ALLOCATE(work2,(msz))
   LIBPAW_ALLOCATE(work3,(msz))
   work3(1:msz)=core_mesh%rad(1:msz)
   call paw_spline(work3,tncore,msz,yp1,ypn,work1)
   call paw_splint(msz,work3,tncore,work1,reduced_mshsz,rcore_mesh%rad,rtncor)
   LIBPAW_DEALLOCATE(work1)
   LIBPAW_DEALLOCATE(work2)
   LIBPAW_DEALLOCATE(work3)
 end if

!==========================================================
!Try to optimize CPU time:
!If coretau mesh size is big, spline tcoretau into a smaller log. mesh

 reduced_taucor=.false.
 if (usekden>=1.and.present(tcoretau)) then
   reduced_taucor=(coretau_mesh%mesh_size>int(reduced_mshsz)).and.(pawtab%usetcore/=0)
   if (reduced_taucor) then
     msz=coretau_mesh%mesh_size
     lstep_tmp=log(0.9999999_dp*coretau_mesh%rmax/reduced_rstep)/dble(reduced_mshsz-2)
     call pawrad_init(rcoretau_mesh,mesh_size=reduced_mshsz,mesh_type=3,&
&    rstep=reduced_rstep,lstep=lstep_tmp)
     LIBPAW_ALLOCATE(rttaucor,(reduced_mshsz))
     call bound_deriv(tcoretau(1:msz),coretau_mesh,msz,yp1,ypn)
     LIBPAW_ALLOCATE(work1,(msz))
     LIBPAW_ALLOCATE(work2,(msz))
     LIBPAW_ALLOCATE(work3,(msz))
     work3(1:msz)=coretau_mesh%rad(1:msz)
     call paw_spline(work3,tcoretau,msz,yp1,ypn,work1)
     call paw_splint(msz,work3,tcoretau,work1,reduced_mshsz,rcoretau_mesh%rad,rttaucor)
     LIBPAW_DEALLOCATE(work1)
     LIBPAW_DEALLOCATE(work2)
     LIBPAW_DEALLOCATE(work3)
   end if
 end if

!==========================================================
!Try to optimize CPU time:
!If vale mesh size is big, spline tnvale into a smaller log. mesh

 if (pawtab%has_tvale==1) then
   reduced_nval=(vale_mesh%mesh_size>int(reduced_mshsz))
   if (reduced_nval) then
     msz=vale_mesh%mesh_size
     lstep_tmp=log(0.9999999_dp*vale_mesh%rmax/reduced_rstep)/dble(reduced_mshsz-2)
     call pawrad_init(rvale_mesh,mesh_size=reduced_mshsz,mesh_type=3,&
&     rstep=reduced_rstep,lstep=lstep_tmp)
     LIBPAW_ALLOCATE(rtnval,(reduced_mshsz))
     call bound_deriv(tnvale(1:msz),vale_mesh,msz,yp1,ypn)
     LIBPAW_ALLOCATE(work1,(msz))
     LIBPAW_ALLOCATE(work2,(msz))
     LIBPAW_ALLOCATE(work3,(msz))
     work3(1:msz)=vale_mesh%rad(1:msz)
     call paw_spline(work3,tnvale,msz,yp1,ypn,work1)
     call paw_splint(msz,work3,tnvale,work1,reduced_mshsz,rvale_mesh%rad,rtnval)
     LIBPAW_DEALLOCATE(work1)
     LIBPAW_DEALLOCATE(work2)
     LIBPAW_DEALLOCATE(work3)
   end if
 else
   reduced_nval=.false.
 end if
!==========================================================
!Compute Vlspl(q) (and second derivative) from Vloc(r)

!Compute Vlspl(q)=q^2.Vloc(q) from vloc(r)
 if(mqgrid_vl>0) then
   if (reduced_vloc) then
     call pawpsp_lo(epsatm,mqgrid_vl,qgrid_vl,vlspl(:,1),rvloc_mesh,rvlocr,yp1,ypn,zion)
   else
     call pawpsp_lo(epsatm,mqgrid_vl,qgrid_vl,vlspl(:,1),vloc_mesh,vlocr,yp1,ypn,zion)
   end if
!  Compute second derivative of Vlspl(q)
   call paw_spline(qgrid_vl,vlspl(:,1),mqgrid_vl,yp1,ypn,vlspl(:,2))
 else
   ! Only to compute epsatm
   epsatm=zero
   if (reduced_vloc) then
     call pawpsp_lo(epsatm,1,tmp_qgrid,tmp_q2vq,rvloc_mesh,rvlocr,yp1,ypn,zion)
   else
     call pawpsp_lo(epsatm,1,tmp_qgrid,tmp_q2vq,vloc_mesh,vlocr,yp1,ypn,zion)
   end if
 end if
!==========================================================
!Compute tcorespl(q) (and second derivative) from tNcore(r)

 pawtab%mqgrid=mqgrid_vl
 xcccrc=core_mesh%rmax
 LIBPAW_ALLOCATE(pawtab%tcorespl,(pawtab%mqgrid,2))

 if(mqgrid_vl>0.and.pawtab%usetcore/=0) then
!  Compute tcorespl(q)=tNc(q) from tNcore(r)
   if (reduced_ncor) then
     call pawpsp_cg(pawtab%dncdq0,pawtab%d2ncdq0,mqgrid_vl,qgrid_vl,pawtab%tcorespl(:,1),rcore_mesh,rtncor,yp1,ypn)
   else
     call pawpsp_cg(pawtab%dncdq0,pawtab%d2ncdq0,mqgrid_vl,qgrid_vl,pawtab%tcorespl(:,1),core_mesh,tncore,yp1,ypn)
   end if
!  Compute second derivative of tcorespl(q)
   call paw_spline(qgrid_vl,pawtab%tcorespl(:,1),mqgrid_vl,yp1,ypn,pawtab%tcorespl(:,2))
 else
   pawtab%tcorespl=zero
   pawtab%dncdq0=zero
   pawtab%d2ncdq0=zero
 end if

!==========================================================
!Compute tcoretauspl(q)

 if (present(tcoretau)) then
   LIBPAW_ALLOCATE(pawtab%tcoretauspl,(pawtab%mqgrid,2*usekden))
   if (usekden==1) then
     if (coretau_mesh%rmax/=xcccrc) then
       write(msg, '(a,a,a)' )&
&       'Core density and core kinetic density should be given on the same grid!',ch10,&
&       'Action : check your pseudopotential (increase tNvale meshSize).'
       MSG_ERROR(msg)
     end if
     if(mqgrid_vl>0) then
!      Compute tcorespl(q)=tNc(q) from tNcore(r)
       if (reduced_taucor) then
         call pawpsp_cg(pawtab%dtaucdq0,qq,mqgrid_vl,qgrid_vl,pawtab%tcoretauspl(:,1),rcoretau_mesh,rttaucor,yp1,ypn)
       else
         call pawpsp_cg(pawtab%dtaucdq0,qq,mqgrid_vl,qgrid_vl,pawtab%tcoretauspl(:,1),coretau_mesh,tcoretau,yp1,ypn)
       end if
!      Compute second derivative of tcorespl(q)
       call paw_spline(qgrid_vl,pawtab%tcoretauspl(:,1),mqgrid_vl,yp1,ypn,pawtab%tcoretauspl(:,2))
     else
       pawtab%tcoretauspl=zero
       pawtab%dtaucdq0=zero
     end if
   end if
 end if

!==========================================================
!Compute tvalespl(q) (and second derivative) from tNvale(r)

 if (pawtab%has_tvale/=0.and.mqgrid_vl>0) then
   LIBPAW_ALLOCATE(pawtab%tvalespl,(pawtab%mqgrid,2))
   if (reduced_nval) then
     call pawpsp_cg(pawtab%dnvdq0,d2nvdq0,mqgrid_vl,qgrid_vl,pawtab%tvalespl(:,1),rvale_mesh,rtnval,yp1,ypn)
     pawtab%tnvale_mesh_size=rvale_mesh%mesh_size
   else
     call pawpsp_cg(pawtab%dnvdq0,d2nvdq0,mqgrid_vl,qgrid_vl,pawtab%tvalespl(:,1),vale_mesh,tnvale,yp1,ypn)
     pawtab%tnvale_mesh_size=vale_mesh%mesh_size
   end if
!  Compute second derivative of tvalespl(q)
   call paw_spline(qgrid_vl,pawtab%tvalespl(:,1),mqgrid_vl,yp1,ypn,pawtab%tvalespl(:,2))
 else
   pawtab%dnvdq0=zero
   pawtab%tnvale_mesh_size=0
 end if

!==================================================
!Compute Ex-correlation energy for the core density

 nspden=1
#if defined LIBPAW_HAVE_LIBXC
 if (ixc<0) nspden=libxc_functionals_nspin()
#endif

 LIBPAW_ALLOCATE(work1,(core_mesh%mesh_size*nspden))
 LIBPAW_ALLOCATE(work2,(core_mesh%mesh_size))
 LIBPAW_ALLOCATE(work3,(1))
 work1(:)=zero;work2(:)=zero;work3(:)=zero
 tmp1 => work1 ; tmp2 => work1

 if (pawxcdev/=0) then
   call pawxcm(ncore,pawtab%exccore,yp1,0,ixc,work2,1,tmp_lmselect,work3,0,non_magnetic_xc,core_mesh%mesh_size,&
&   nspden,4,pawang_tmp,core_mesh,pawxcdev,work1,1,0,tmp1,xclevel,xc_denpos)
 else
   if (present(tcoretau)) then
     call pawxc(ncore,pawtab%exccore,yp1,ixc,work2,work1,1,tmp_lmselect,work3,0,0,non_magnetic_xc,core_mesh%mesh_size,&
&     nspden,4,pawang_tmp,core_mesh,tmp1,1,0,tmp2,xclevel,xc_denpos,coretau=tcoretau)
   else
     call pawxc(ncore,pawtab%exccore,yp1,ixc,work2,work1,1,tmp_lmselect,work3,0,0,non_magnetic_xc,core_mesh%mesh_size,&
&     nspden,4,pawang_tmp,core_mesh,tmp1,1,0,tmp2,xclevel,xc_denpos)
   end if
 end if

 LIBPAW_DEALLOCATE(work1)
 LIBPAW_DEALLOCATE(work2)
 LIBPAW_DEALLOCATE(work3)

!==================================================
!Compute atomic contribution to Dij (Dij0)
!if not already in memory
 if ((.not.has_dij0).and.(pawtab%has_kij==2.or.pawtab%has_kij==-1)) then
   LIBPAW_ALLOCATE(pawtab%dij0,(pawtab%lmn2_size))
   if (reduced_vloc) then
     call atompaw_dij0(pawtab%indlmn,pawtab%kij,pawtab%lmn_size,ncore,0,pawtab,pawrad,core_mesh,&
&                      rvloc_mesh,rvlocr,znucl)
   else
     call atompaw_dij0(pawtab%indlmn,pawtab%kij,pawtab%lmn_size,ncore,0,pawtab,pawrad,core_mesh,&
&                      vloc_mesh,vlocr,znucl)
   end if
   has_dij0=.true.
 end if
!==================================================
!Compute kinetic operator contribution to Dij

 if (pawtab%has_kij==1.and.has_dij0) then
   LIBPAW_ALLOCATE(pawtab%kij,(pawtab%lmn2_size))
   call atompaw_kij(pawtab%indlmn,pawtab%kij,pawtab%lmn_size,ncore,0,1,pawtab,pawrad,core_mesh,&
&                   vloc_mesh,vlocr,znucl)
   pawtab%has_kij=2
 end if

!pawtab%has_kij=-1 means that kij does not have to be kept in memory
 if (pawtab%has_kij==-1) then
   LIBPAW_DEALLOCATE(pawtab%kij)
   pawtab%has_kij=0
 end if

!==========================================================
!If projectors have to be kept in memory, we need
!them on the main radial mesh (so, spline them if necessary)

 if (pawtab%has_tproj>0) then
   if ((tproj_mesh%mesh_type/=pawrad%mesh_type).or.&
&      (tproj_mesh%rstep    /=pawrad%rstep).or.&
&      (tproj_mesh%lstep    /=pawrad%lstep)) then
     ir=pawrad_ifromr(pawrad,tproj_mesh%rmax)
     call pawrad_init(tproj_mesh_new,mesh_size=ir,mesh_type=pawrad%mesh_type,&
&                     rstep=pawrad%rstep,lstep=pawrad%lstep)
     LIBPAW_ALLOCATE(pawtab%tproj,(tproj_mesh_new%mesh_size,pawtab%basis_size))
     LIBPAW_ALLOCATE(work1,(tproj_mesh%mesh_size))
     do ib=1,pawtab%basis_size
       call bound_deriv(tproj(:,ib),tproj_mesh,tproj_mesh%mesh_size,yp1,ypn)
       call paw_spline(tproj_mesh%rad,tproj(:,ib),tproj_mesh%mesh_size,yp1,ypn,work1)
       call paw_splint(tproj_mesh%mesh_size,tproj_mesh%rad,tproj(:,ib),work1,&
&           tproj_mesh_new%mesh_size,tproj_mesh_new%rad,pawtab%tproj(:,ib))
     end do
     LIBPAW_DEALLOCATE(work1)
     call pawrad_free(tproj_mesh_new)
   else
     LIBPAW_ALLOCATE(pawtab%tproj,(tproj_mesh%mesh_size,pawtab%basis_size))
     pawtab%tproj(:,:)=tproj(:,:)
   end if
   pawtab%has_tproj=2
 end if

!==========================================================
!Free temporary allocated space

 if (pawtab%has_tvale==1)  then
   LIBPAW_DEALLOCATE(nhat)
 end if
 if (reduced_vloc) then
   call pawrad_free(rvloc_mesh)
   LIBPAW_DEALLOCATE(rvlocr)
 end if
 if (reduced_ncor)  then
   call pawrad_free(rcore_mesh)
   LIBPAW_DEALLOCATE(rtncor)
 end if
  if (reduced_taucor)  then
   call pawrad_free(rcoretau_mesh)
   LIBPAW_DEALLOCATE(rttaucor)
 end if
 if (reduced_nval)  then
   call pawrad_free(rvale_mesh)
   LIBPAW_DEALLOCATE(rtnval)
 end if
 if (pawxcdev==0)  then
   LIBPAW_DEALLOCATE(pawang_tmp%angwgth)
   LIBPAW_DEALLOCATE(pawang_tmp%anginit)
   LIBPAW_DEALLOCATE(pawang_tmp%ylmr)
   LIBPAW_DEALLOCATE(pawang_tmp%ylmrgr)
 end if

end subroutine pawpsp_calc
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_calc_d5
!! NAME
!!  pawpsp_calc_d5
!!
!! FUNCTION
!!  Compute the first to the 5th derivatives of
!!  a given function in a pawrad mesh
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_calc_d5(mesh,mesh_size,tcoredens)

!Arguments ------------------------------------
 integer,intent(in) :: mesh_size
 type(pawrad_type),intent(in) :: mesh
 real(dp),intent(inout) :: tcoredens(mesh_size,6)

!Local variables-------------------------------
 integer,parameter :: it=1 !number of steps for smoothing function
 logical,parameter :: use_smooth=.true.

! *************************************************************************

!calculate first derivative from density,
!and store it
 call nderiv_gen(tcoredens(:,2),tcoredens(:,1),mesh)

!get second derivative from density, and store it
 call paw_spline(mesh%rad,tcoredens(:,1),mesh_size,&
&                zero,zero,tcoredens(:,3))

!smooth functions, to avoid numerical instabilities
 if(use_smooth) then
   call paw_smooth(tcoredens(:,2),mesh_size,it)
   call paw_smooth(tcoredens(:,3),mesh_size,it)
 end if

!get third derivative from first derivative:
 call paw_spline(mesh%rad,tcoredens(:,2),mesh_size,&
&                zero,zero,tcoredens(:,4))

!get fourth derivative from second derivative:
 call paw_spline(mesh%rad,tcoredens(:,3),mesh_size,&
&                zero,zero,tcoredens(:,5))

!smooth 3rd and 4th order derivatives
 if(use_smooth) then
   call paw_smooth(tcoredens(:,4),mesh_size,it)
   call paw_smooth(tcoredens(:,5),mesh_size,it)
 end if

!get fifth derivative from third derivative:
 call paw_spline(mesh%rad,tcoredens(:,4),mesh_size,&
&                zero,zero,tcoredens(:,6))

!smooth 5th order derivative
 if(use_smooth) then
   call paw_smooth(tcoredens(:,6),mesh_size,it)
 end if

end subroutine pawpsp_calc_d5
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_vhar2rho
!! NAME
!!  pawpsp_vhar2rho
!!
!! FUNCTION
!!  gets rho(r) from v(r), solving the Poisson equation
!!  \lap v(r) =  4 \pi rho(r)
!!
!! INPUTS
!! radmesh = radial grid (datastructure)
!! vv(:)= potential
!!
!! OUTPUT
!!  rho(:)= density
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_vhar2rho(radmesh,rho,vv)

!Arguments ------------------------------------
 type(pawrad_type),intent(in) :: radmesh
 real(dp), intent(in) :: vv(:)
 real(dp), intent(out):: rho(:)

!Local variables-------------------------------
 integer :: nr
 real(dp) :: dfdr(radmesh%mesh_size),d2fdr(radmesh%mesh_size)

! *************************************************************************

 nr=size(vv)
 if (nr/=size(rho)) then
   MSG_BUG('wrong sizes!')
 end if

!Laplacian =
!\frac{\partial^2}{\partial r^2} + 2/r \frac{\partial}{\partial r}

!Calculate derivatives
 call nderiv_gen(dfdr(1:nr),vv,radmesh,der2=d2fdr(1:nr))

 rho(2:nr)=d2fdr(2:nr) + 2._dp*dfdr(2:nr)/radmesh%rad(2:nr)
 call pawrad_deducer0(rho,nr,radmesh)

 rho(1:nr)=-rho(1:nr)/(4._dp*pi)

end subroutine pawpsp_vhar2rho
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_wvl_calc
!! NAME
!! pawpsp_wvl_calc
!!
!! FUNCTION
!! Performs tests and compute data related to pspcod=7 or 17 ("PAW pseudopotentials")
!!
!! INPUTS
!!  tnvale(vale_mesh%mesh_size)= pseudo valence density (+ nhat in output)
!!  usewvl= flag for wavelets method
!!  vale_mesh<type(pawrad_type)>= radial mesh for the valence density
!!  vloc_mesh<type(pawrad_type)>= radial mesh for the local potential
!!  vlocr(vloc_mesh%mesh_size)= local potential according to vlocopt.
!!
!! OUTPUT
!!  Sets pawtab%rholoc
!!
!! SIDE EFFECTS
!!  pawtab <type(pawtab_type)>= objects are modified
!!
!! NOTES
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_wvl_calc(pawtab,tnvale,usewvl,vale_mesh,vloc_mesh,vlocr)

!Arguments ------------------------------------
!scalars
 integer,intent(in)::usewvl
 type(pawrad_type),intent(in) :: vale_mesh
 type(pawtab_type),intent(inout) :: pawtab
 type(pawrad_type),intent(in) ::vloc_mesh

!arrays
 real(dp),intent(in) :: tnvale(vale_mesh%mesh_size*pawtab%has_tvale)
 real(dp),intent(in) :: vlocr(vloc_mesh%mesh_size)


!Local variables ------------------------------
!scalars
 integer :: msz
 character(len=500) :: msg
!arrays

! *************************************************************************

!If usewvl flag is on, we must have the pawtab%wvl pointer allocated
 if (pawtab%has_wvl==0) then
   msg='pawtab%has_wvl flag should be on o entry'
   MSG_BUG(msg)
 end if
 call wvlpaw_allocate(pawtab%wvl)

!==========================================================
!Change mesh_size of tvalespl
!Compute second derivative from tNvale(r)

 if (pawtab%has_tvale/=0) then
   if(usewvl==1) then
     if(allocated(pawtab%tvalespl)) then
       LIBPAW_DEALLOCATE(pawtab%tvalespl)
     end if
     LIBPAW_ALLOCATE(pawtab%tvalespl,(vale_mesh%mesh_size,2))
     pawtab%tnvale_mesh_size=vale_mesh%mesh_size
     pawtab%tvalespl(:,1)=tnvale
!    Compute second derivative of tvalespl(r)
     call paw_spline(vale_mesh%rad,pawtab%tvalespl(:,1),vale_mesh%mesh_size,zero,zero,pawtab%tvalespl(:,2))
   end if
 else
   pawtab%dnvdq0=zero
   pawtab%tnvale_mesh_size=0
 end if

!==========================================================
!Save rholoc:
!Get local density from local potential
!use the poisson eq.
 msz=vloc_mesh%mesh_size
 call wvlpaw_rholoc_free(pawtab%wvl%rholoc)
 LIBPAW_ALLOCATE(pawtab%wvl%rholoc%d,(msz,4))
 LIBPAW_ALLOCATE(pawtab%wvl%rholoc%rad,(msz))
 pawtab%wvl%rholoc%msz=msz
 pawtab%wvl%rholoc%rad(1:msz)=vloc_mesh%rad(1:msz)

!get rho from v:
 call pawpsp_vhar2rho(vloc_mesh,pawtab%wvl%rholoc%d(:,1),vlocr)
!
!get second derivative, and store it
 call paw_spline(pawtab%wvl%rholoc%rad,pawtab%wvl%rholoc%d(:,1),pawtab%wvl%rholoc%msz,&
& zero,zero,pawtab%wvl%rholoc%d(:,2))

!save also vlocr:
 pawtab%wvl%rholoc%d(:,3)=vlocr

!get second derivative, and store it
 call paw_spline(pawtab%wvl%rholoc%rad,vlocr,pawtab%wvl%rholoc%msz,&
& zero,zero,pawtab%wvl%rholoc%d(:,4))

!Test
!do ii=1,pawtab%wvl%rholoc%msz
!write(503,'(3(f16.10,x))')pawtab%wvl%rholoc%rad(ii),pawtab%wvl%rholoc%d(ii,1),pawtab%wvl%rholoc%d(ii,3)
!end do
!
!Do splint
!
!nmesh=4000
!rread1= (9.9979999d0/real(nmesh-1,dp)) ! 0.0025001d0  !step
!allocate(raux1(nmesh),raux2(nmesh))
!do ii=1,nmesh
!raux1(ii)=rread1*real(ii-1,dp)  !mesh
!end do
!call splint(pawtab%wvl%rholoc%msz,pawtab%wvl%rholoc%rad,pawtab%wvl%rholoc%d(:,1),pawtab%wvl%rholoc%d(:,2),&
!&  nmesh,raux1,raux2,ierr)
!do ii=1,nmesh
!write(401,'(10(f20.7,x))')raux1(ii),raux2(ii),raux2(ii)*raux1(ii)**2
!end do
!deallocate(raux1,raux2)

end subroutine pawpsp_wvl_calc
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_17in
!! NAME
!! pawpsp_17in
!!
!! FUNCTION
!! Initialize pspcod=17 ("PAW  XML pseudopotentials"):
!! continue to read the corresponding file and compute the form factors
!!
!! INPUTS
!!  ipsp= id in the array of the currently read pseudo.
!!  ixc=exchange-correlation choice from main routine data file
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  lnmax=max. number of (l,n) components over all type of psps
!!            angular momentum of nonlocal pseudopotential
!!  mmax=max number of pts in real space grid (already read in the psp file header)
!!  mqgrid_ff=dimension of q (or G) grid for nl form factors (array ffspl)
!!  mqgrid_vl=dimension of q (or G) grid for Vloc (array vlspl)
!!  pawxcdev=choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  pspheads= header of the current pseudopotential
!!  qgrid_ff(psps%mqgrid_ff)=values of q on grid from 0 to qmax (bohr^-1) for nl form factors
!!  qgrid_vl(psps%mqgrid_vl)=values of q on grid from 0 to qmax (bohr^-1) for Vloc
!!  xclevel= XC functional level
!!  zion=nominal valence of atom as specified in psp file
!!  znucl=atomic number of atom as specified in input file to main routine
!!
!! OUTPUT
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$(hartree)
!!  ffspl(psps%mqgrid_ff,2,psps%lnmax)=form factor f_l(q) and second derivative
!!   from spline fit for each angular momentum and each projector;
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  vlspl(psps%mqgrid_vl,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  wvl_crmult,wvl_frmult= variables definining the fine and coarse grids in a wavelets calculation
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!  xcccrc=XC core correction cutoff radius (bohr) from psp file
!!
!! NOTES
!!  Spin-orbit not yet implemented (to be done)
!!  Comments:
!!  * mesh_type= type of radial mesh
!!  mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!  mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!  mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!  mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!  * radial shapefunction type
!!  shape_type=-1 ; gl(r)=numeric (read from psp file)
!!  shape_type= 1 ; gl(r)=k(r).r^l; k(r)=exp[-(r/sigma)**lambda]
!!  shape_type= 2 ; gl(r)=k(r).r^l; k(r)=[sin(pi*r/rshp)/(pi*r/rshp)]**2 if r<=rshp
!!  shape_type= 3 ; gl(r)=Alpha(1,l)*jl(q(1,l)*r)+Alpha(2,l)*jl(q(2,l)*r) for each l
!!
!! PARENTS
!!      m_pawpsp,m_pspini
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_17in(epsatm,ffspl,icoulomb,ipsp,ixc,lmax,&
& lnmax,mmax,mqgrid_ff,mqgrid_vl,pawpsp_header,pawrad,pawtab,&
& pawxcdev, qgrid_ff,qgrid_vl,usewvl,usexcnhat_in,vlspl,xcccrc,&
& xclevel,xc_denpos,zion,znucl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ipsp,ixc,lmax,lnmax,mqgrid_ff,mqgrid_vl,pawxcdev,usexcnhat_in
 integer,intent(inout) ::mmax
 integer,intent(in) :: xclevel,icoulomb,usewvl
 real(dp),intent(in) :: xc_denpos,zion,znucl
 real(dp),intent(out) :: epsatm,xcccrc
 type(pawpsp_header_type),intent(in) :: pawpsp_header
 type(pawrad_type),intent(inout) :: pawrad
 type(pawtab_type),intent(inout) :: pawtab
!arrays
 real(dp),intent(in) :: qgrid_ff(mqgrid_ff),qgrid_vl(mqgrid_vl)
 real(dp),intent(inout) :: ffspl(mqgrid_ff,2,lnmax)
 real(dp),intent(out) :: vlspl(mqgrid_vl,2)

!Local variables ------------------------------
!scalars
 integer :: has_v_minushalf,ib,icoremesh,icoretaumesh,il,ilm,ilmn,ilmn0,iln,imainmesh,imsh,iprojmesh
 integer :: ir,iread1,ishpfmesh,ivalemesh,ivlocmesh,j0lmn,jlm,pngau
 integer :: jlmn,jln,klmn,msz,nmesh,nval,pspversion,shft,sz10,usexcnhat,vlocopt
 real(dp), parameter :: rmax_vloc=10.0_dp
 real(dp) :: fourpi,occ,rc,yp1,ypn
 logical :: save_core_msz
 character(len=500) :: msg
 type(pawrad_type) :: core_mesh,coretau_mesh,shpf_mesh,tproj_mesh,vale_mesh,vloc_mesh
!arrays
 integer,allocatable :: mesh_shift(:),nprj(:)
 real(dp),allocatable :: kij(:),ncore(:),shpf(:,:),tncore(:),coretau(:)
 real(dp),allocatable :: tcoretau(:),tnvale(:),tproj(:,:),vhnzc(:),vlocr(:)
 real(dp),allocatable :: work1(:),work2(:),work3(:),work4(:)
 type(pawrad_type),allocatable :: radmesh(:)

!************************************************************************

 if (.False.) write(std_out,*) ipsp

!==========================================================
!Destroy everything in pawtab but optional flags
 call pawtab_free(pawtab)
!Destroy everything in pawrad
 call pawrad_free(pawrad)

!==========================================================
!Initialize useful data

 pawtab%usexcnhat=usexcnhat_in
 fourpi=4*acos(-1.d0)
 pspversion=pawpsp_header%pawver
 save_core_msz=(usewvl==1 .or. icoulomb .ne. 0)
 imainmesh=-1;icoremesh=-1;icoretaumesh=-1;iprojmesh=-1
 ishpfmesh=-1;ivalemesh=-1;ivlocmesh=-1

!==========================================================
!Initialize partial waves quantum numbers

 pawtab%basis_size=pawpsp_header%basis_size
 LIBPAW_ALLOCATE(pawtab%orbitals,(pawtab%basis_size))
 do ib=1,pawtab%basis_size
   pawtab%orbitals(ib)=paw_setuploc%valence_states%state(ib)%ll
 end do

!==========================================================
!Initialize various dims and indexes

 pawtab%lmn_size=pawpsp_header%lmn_size
 pawtab%lmn2_size=pawtab%lmn_size*(pawtab%lmn_size+1)/2
 pawtab%l_size=2*maxval(pawtab%orbitals)+1
 pawtab%ij_size=pawtab%basis_size*(pawtab%basis_size+1)/2

!indlmn calculation (indices for (l,m,n) basis)
 if (allocated(pawtab%indlmn)) then
   LIBPAW_DEALLOCATE(pawtab%indlmn)
 end if
 LIBPAW_ALLOCATE(pawtab%indlmn,(6,pawtab%lmn_size))
 pawtab%indlmn(:,:)=0
 LIBPAW_BOUND1_ALLOCATE(nprj,BOUNDS(0,maxval(pawtab%orbitals)))
 ilmn=0;iln=0;nprj=0
 do ib=1,pawtab%basis_size
   il=pawtab%orbitals(ib)
   nprj(il)=nprj(il)+1
   iln=iln+1
   do ilm=1,2*il+1
     pawtab%indlmn(1,ilmn+ilm)=il
     pawtab%indlmn(2,ilmn+ilm)=ilm-(il+1)
     pawtab%indlmn(3,ilmn+ilm)=nprj(il)
     pawtab%indlmn(4,ilmn+ilm)=il*il+ilm
     pawtab%indlmn(5,ilmn+ilm)=iln
     pawtab%indlmn(6,ilmn+ilm)=1
   end do
   ilmn=ilmn+2*il+1
 end do
 LIBPAW_DEALLOCATE(nprj)
!Are ilmn (found here) and pawtab%lmn_size compatibles ?
 if (ilmn/=pawtab%lmn_size) then
   write(msg, '(a,a,a,a,a)' )&
&   'Calculated lmn size differs from',ch10,&
&   'lmn_size read from pseudo !',ch10,&
&   'Action: check your pseudopotential file.'
   MSG_ERROR(msg)
 end if

!==========================================================
!Read and initialize radial meshes

 nmesh=paw_setuploc%ngrid
 LIBPAW_DATATYPE_ALLOCATE(radmesh,(nmesh))
 LIBPAW_ALLOCATE(mesh_shift,(nmesh))
 do imsh=1,nmesh
   radmesh(imsh)%mesh_type=-1
   radmesh(imsh)%rstep=zero
   radmesh(imsh)%lstep=zero
   mesh_shift(imsh)=0
   select case(trim(paw_setuploc%radial_grid(imsh)%eq))
     case("r=a*exp(d*i)")
       mesh_shift(imsh)=1
       radmesh(imsh)%mesh_type=3
       radmesh(imsh)%mesh_size=paw_setuploc%radial_grid(imsh)%iend &
&                             -paw_setuploc%radial_grid(imsh)%istart+1+mesh_shift(imsh)
       radmesh(imsh)%rstep=paw_setuploc%radial_grid(imsh)%aa
       radmesh(imsh)%lstep=paw_setuploc%radial_grid(imsh)%dd
     case("r=a*i/(1-b*i)")
       write(msg, '(3a)' )&
&       'The grid r=a*i/(1-b*i) is not implemented in ABINIT !',ch10,&
&       'Action: check your psp file.'
       MSG_ERROR(msg)
     case("r=a*i/(n-i)")
       mesh_shift(imsh)=0
       radmesh(imsh)%mesh_type=5
       radmesh(imsh)%mesh_size=paw_setuploc%radial_grid(imsh)%iend &
&                             -paw_setuploc%radial_grid(imsh)%istart+1+mesh_shift(imsh)
       radmesh(imsh)%rstep=paw_setuploc%radial_grid(imsh)%aa
       radmesh(imsh)%lstep=dble(paw_setuploc%radial_grid(imsh)%nn)
     case("r=a*(exp(d*i)-1)")
       mesh_shift(imsh)=0
       radmesh(imsh)%mesh_type=2
       radmesh(imsh)%mesh_size=paw_setuploc%radial_grid(imsh)%iend &
&                             -paw_setuploc%radial_grid(imsh)%istart+1+mesh_shift(imsh)
       if(paw_setuploc%radial_grid(imsh)%istart==1)radmesh(imsh)%mesh_size=radmesh(imsh)%mesh_size+1
       radmesh(imsh)%rstep=paw_setuploc%radial_grid(imsh)%aa
       radmesh(imsh)%lstep=paw_setuploc%radial_grid(imsh)%dd
     case("r=d*i")
       mesh_shift(imsh)=0
       radmesh(imsh)%mesh_type=1
       radmesh(imsh)%mesh_size=paw_setuploc%radial_grid(imsh)%iend &
&                             -paw_setuploc%radial_grid(imsh)%istart+1+mesh_shift(imsh)
       if(paw_setuploc%radial_grid(imsh)%istart==1)radmesh(imsh)%mesh_size=radmesh(imsh)%mesh_size+1
       radmesh(imsh)%rstep=paw_setuploc%radial_grid(imsh)%dd
     case("r=(i/n+a)^5/a-a^4")
       write(msg, '(3a)' )&
&       'The grid r=(i/n+a)^5/a-a^4 is not implemented in ABINIT !',ch10,&
&       'Action: check your psp file.'
       MSG_ERROR(msg)
   end select
 end do

!Initialize radial meshes
 do imsh=1,nmesh
   call pawrad_init(radmesh(imsh))
 end do

 pawtab%rpaw=pawpsp_header%rpaw

!==========================================================
!Here reading shapefunction parameters

 pawtab%shape_type=pawpsp_header%shape_type
 pawtab%shape_lambda=-1;pawtab%shape_sigma=1.d99
 pawtab%rshp=pawpsp_header%rshp
 pawtab%shape_lambda=paw_setuploc%shape_function%lamb
 if(trim(paw_setuploc%shape_function%gtype)=="gauss")pawtab%shape_lambda=2
 pawtab%shape_sigma=paw_setuploc%shape_function%rc
!If shapefunction type is gaussian, check exponent
 if (pawtab%shape_type==1) then
   if (pawtab%shape_lambda<2) then
     write(msg, '(3a)' )&
&     'For a gaussian shape function, exponent lambda must be >1 !',ch10,&
&     'Action: check your psp file.'
     MSG_ERROR(msg)
   end if
 end if

!If shapefunction type is Bessel, deduce here its parameters from rc
 if (pawtab%shape_type==3) then
   LIBPAW_ALLOCATE(pawtab%shape_alpha,(2,pawtab%l_size))
   LIBPAW_ALLOCATE(pawtab%shape_q,(2,pawtab%l_size))
   rc=pawtab%rshp;if (rc<1.d-8) rc=pawtab%rpaw
   do il=1,pawtab%l_size
     call atompaw_shapebes(pawtab%shape_alpha(1:2,il),pawtab%shape_q(1:2,il),il-1,rc)
   end do
 end if

!==========================================================
!Mirror pseudopotential parameters to the output and log files

 write(msg,'(a,i2)')' Pseudopotential format is: paw',pspversion
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')
 write(msg,'(2(a,i3),a,64i4)') &
& ' basis_size (lnmax)=',pawtab%basis_size,' (lmn_size=',&
& pawtab%lmn_size,'), orbitals=',pawtab%orbitals(1:pawtab%basis_size)
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')
 write(msg,'(a,f11.8)')' Spheres core radius: rc_sph=',pawtab%rpaw
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')
 write(msg,'(a,i1,a)')' ',nmesh,' radial meshes are used:'
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

 do imsh=1,nmesh
   if (radmesh(imsh)%mesh_type==1) &
&   write(msg,'(a,i1,a,i4,a,g12.5)') &
&   '  - mesh ',imsh,': r(i)=step*(i-1), size=',radmesh(imsh)%mesh_size,&
&   ' , step=',radmesh(imsh)%rstep
   if (radmesh(imsh)%mesh_type==2) &
&   write(msg,'(a,i1,a,i4,2(a,g12.5))') &
&   '  - mesh ',imsh,': r(i)=AA*[exp(BB*(i-1))-1], size=',radmesh(imsh)%mesh_size,&
&   ' , AA=',radmesh(imsh)%rstep,' BB=',radmesh(imsh)%lstep
   if (radmesh(imsh)%mesh_type==3) &
&   write(msg,'(a,i1,a,i4,2(a,g12.5))') &
&   '  - mesh ',imsh,': r(i)=AA*exp(BB*(i-2)), size=',radmesh(imsh)%mesh_size,&
&   ' , AA=',radmesh(imsh)%rstep,' BB=',radmesh(imsh)%lstep
   if (radmesh(imsh)%mesh_type==4) &
&   write(msg,'(a,i1,a,i4,a,g12.5)') &
&   '  - mesh ',imsh,': r(i)=-AA*ln(1-(i-1)/n), n=size=',radmesh(imsh)%mesh_size,&
&   ' , AA=',radmesh(imsh)%rstep
   if (radmesh(imsh)%mesh_type==5) &
&   write(msg,'(a,i1,a,i4,2(a,g12.5))') &
&   '  - mesh ',imsh,': r(i)=-AA*i/(NN-i)), n=size=',radmesh(imsh)%mesh_size,&
&   ' , AA=',radmesh(imsh)%rstep,' NN=',radmesh(imsh)%lstep
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 end do
 if (pawtab%shape_type==-1) then
   write(msg,'(a)')&
   ' Shapefunction is NUMERIC type: directly read from atomic data file'
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 end if
 if (pawtab%shape_type==1) then
   write(msg,'(2a,a,f6.3,a,i3)')&
&   ' Shapefunction is EXP type: shapef(r)=exp(-(r/sigma)**lambda)',ch10,&
&   '                            with sigma=',pawtab%shape_sigma,' and lambda=',pawtab%shape_lambda
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 end if
 if (pawtab%shape_type==2) then
   write(msg,'(a)')&
   ' Shapefunction is SIN type: shapef(r)=[sin(pi*r/rshp)/(pi*r/rshp)]**2'
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 end if
 if (pawtab%shape_type==3) then
   write(msg,'(a)')&
&   ' Shapefunction is BESSEL type: shapef(r,l)=aa(1,l)*jl(q(1,l)*r)+aa(2,l)*jl(q(2,l)*r)'
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 end if
 if (pawtab%rshp<1.d-8) then
   write(msg,'(a)') ' Radius for shape functions = sphere core radius'
 else
   write(msg,'(a,f11.8)') ' Radius for shape functions = ',pawtab%rshp
 end if
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

!==========================================================
!Perfom tests

!Are lmax and orbitals compatibles ?
 if (lmax/=maxval(pawtab%orbitals)) then
   write(msg, '(a,a,a)' )&
&   'lmax /= MAX(orbitals) !',ch10,&
&   'Action: check your pseudopotential file.'
   MSG_ERROR(msg)
 end if

!Only mesh_type=1,2, 3 or 5 allowed
 do imsh=1,nmesh
   if (radmesh(imsh)%mesh_type>5) then
     write(msg, '(a,a,a)' )&
&     'Only mesh types 1,2,3 or 5 allowed !',ch10,&
&     'Action : check your pseudopotential or input file.'
     MSG_ERROR(msg)
   end if
 end do

!==========================================================
!Read tabulated atomic data

!---------------------------------
!Read wave-functions (phi)

 do ib=1,pawtab%basis_size
   if (ib==1) then
     do imsh=1,nmesh
       if(trim(paw_setuploc%ae_partial_wave(1)%grid)==trim(paw_setuploc%radial_grid(imsh)%id)) then
         mmax=radmesh(imsh)%mesh_size
         call pawrad_init(pawrad,mesh_size=mmax,mesh_type=radmesh(imsh)%mesh_type, &
&         rstep=radmesh(imsh)%rstep,lstep=radmesh(imsh)%lstep,r_for_intg=pawtab%rpaw)
         pawtab%partialwave_mesh_size=pawrad%mesh_size
         pawtab%mesh_size=pawrad_ifromr(pawrad,pawtab%rpaw)+5
         pawtab%mesh_size=min(pawtab%mesh_size,pawrad%mesh_size)
         if (pawtab%mesh_size>pawrad%mesh_size-2) pawtab%mesh_size=pawrad%mesh_size
         imainmesh=imsh
         exit
       end if
     end do
     LIBPAW_ALLOCATE(pawtab%phi,(pawtab%partialwave_mesh_size,pawtab%basis_size))
   else if (trim(paw_setuploc%ae_partial_wave(ib)%grid)/=trim(paw_setuploc%radial_grid(imainmesh)%id)) then
     write(msg, '(a,a,a)' )&
&     'All Phi and tPhi must be given on the same radial mesh !',ch10,&
&     'Action: check your pseudopotential file.'
     MSG_ERROR(msg)
   end if
   shft=mesh_shift(imainmesh)
   pawtab%phi(1+shft:pawtab%partialwave_mesh_size,ib)= &
&         paw_setuploc%ae_partial_wave(ib)%data(1:pawtab%partialwave_mesh_size-shft) &
&        *pawrad%rad(1+shft:pawtab%partialwave_mesh_size)
   if (shft==1) pawtab%phi(1,ib)=zero
 end do
 write(msg,'(a,i4)') ' mmax= ',mmax
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

!---------------------------------
!Read pseudo wave-functions (tphi)

 LIBPAW_ALLOCATE(pawtab%tphi,(pawtab%partialwave_mesh_size,pawtab%basis_size))
 do ib=1,pawtab%basis_size

   if(trim(paw_setuploc%pseudo_partial_wave(ib)%grid)/=trim(paw_setuploc%radial_grid(imainmesh)%id)) then
     write(msg, '(a,a,a)' )&
&     'All Phi and tPhi must be given on the same radial mesh !',ch10,&
&     'Action: check your pseudopotential file.'
      MSG_ERROR(msg)
   end if
   shft=mesh_shift(imainmesh)
   pawtab%tphi(1+shft:pawtab%partialwave_mesh_size,ib)=&
&         paw_setuploc%pseudo_partial_wave(ib)%data(1:pawtab%partialwave_mesh_size-shft) &
&        *pawrad%rad(1+shft:pawtab%partialwave_mesh_size)
   if (shft==1) pawtab%tphi(1,ib)=zero
 end do
 write(msg,'(a,i1)') &
& ' Radial grid used for partial waves is grid ',imainmesh
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

!---------------------------------
!Read projectors (tproj)

 if (allocated(paw_setuploc%projector_fit)) then
   call wvlpaw_allocate(pawtab%wvl)
   LIBPAW_ALLOCATE(pawtab%wvl%pngau,(pawtab%basis_size))
   do ib=1,pawtab%basis_size
     pawtab%wvl%pngau(ib) = paw_setuploc%projector_fit(ib)%ngauss
   end do
   pawtab%wvl%ptotgau = sum(pawtab%wvl%pngau) * 2
   LIBPAW_ALLOCATE(pawtab%wvl%parg,(2,pawtab%wvl%ptotgau))
   LIBPAW_ALLOCATE(pawtab%wvl%pfac,(2,pawtab%wvl%ptotgau))
   pngau = 1
   do ib=1,pawtab%basis_size
     ! Complex gaussian
     pawtab%wvl%parg(:,pngau:pngau + pawtab%wvl%pngau(ib) - 1) = &
     & paw_setuploc%projector_fit(ib)%expos(:,1:pawtab%wvl%pngau(ib))
     pawtab%wvl%pfac(:,pngau:pngau + pawtab%wvl%pngau(ib) - 1) = &
     & paw_setuploc%projector_fit(ib)%factors(:,1:pawtab%wvl%pngau(ib))
     pngau = pngau + pawtab%wvl%pngau(ib)
     ! Conjugate gaussian
     pawtab%wvl%parg(1,pngau:pngau + pawtab%wvl%pngau(ib) - 1) = &
     & paw_setuploc%projector_fit(ib)%expos(1,1:pawtab%wvl%pngau(ib))
     pawtab%wvl%parg(2,pngau:pngau + pawtab%wvl%pngau(ib) - 1) = &
     & -paw_setuploc%projector_fit(ib)%expos(2,1:pawtab%wvl%pngau(ib))
     pawtab%wvl%pfac(1,pngau:pngau + pawtab%wvl%pngau(ib) - 1) = &
     & paw_setuploc%projector_fit(ib)%factors(1,1:pawtab%wvl%pngau(ib))
     pawtab%wvl%pfac(2,pngau:pngau + pawtab%wvl%pngau(ib) - 1) = &
     & -paw_setuploc%projector_fit(ib)%factors(2,1:pawtab%wvl%pngau(ib))
     pngau = pngau + pawtab%wvl%pngau(ib)
     pawtab%wvl%pngau(ib) = pawtab%wvl%pngau(ib) * 2
   end do
   pawtab%has_wvl=2
 else
   !Nullify wavelet objects for safety:
   pawtab%has_wvl=0
   call wvlpaw_free(pawtab%wvl)
 end if
 do ib=1,pawtab%basis_size
   if (ib==1) then
     do imsh=1,nmesh
       if(trim(paw_setuploc%projector_function(1)%grid)==trim(paw_setuploc%radial_grid(imsh)%id)) then
         iprojmesh=imsh
         exit
       end if
     end do
     call pawrad_copy(radmesh(iprojmesh),tproj_mesh)
     LIBPAW_ALLOCATE(tproj,(tproj_mesh%mesh_size,pawtab%basis_size))
   else if (trim(paw_setuploc%projector_function(ib)%grid)/=trim(paw_setuploc%radial_grid(iprojmesh)%id)) then
     write(msg, '(a,a,a)' )&
&     'All tprojectors must be given on the same radial mesh !',ch10,&
&     'Action: check your pseudopotential file.'
     MSG_ERROR(msg)
   end if
   shft=mesh_shift(iprojmesh)
   tproj(1+shft:tproj_mesh%mesh_size,ib)=paw_setuploc%projector_function(ib)%data(1:tproj_mesh%mesh_size-shft)&
&   *tproj_mesh%rad(1+shft:tproj_mesh%mesh_size)
   if (shft==1) tproj(1,ib)=zero
 end do
 write(msg,'(a,i1)') &
& ' Radial grid used for projectors is grid ',iprojmesh
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

!---------------------------------
!Read core density (coredens)

 do imsh=1,nmesh
   if(trim(paw_setuploc%ae_core_density%grid)==trim(paw_setuploc%radial_grid(imsh)%id)) then
     icoremesh=imsh
     exit
   end if
 end do
 call pawrad_copy(radmesh(icoremesh),core_mesh)
 if ((radmesh(icoremesh)%mesh_type/=pawrad%mesh_type).or.&
& (radmesh(icoremesh)%rstep    /=pawrad%rstep)    .or.&
& (radmesh(icoremesh)%lstep    /=pawrad%lstep)) then
   write(msg, '(a,a,a,a,a)' )&
&   'Ncore must be given on a radial mesh with the same',ch10,&
&   'type and step(s) than the main radial mesh (mesh for Phi) !',ch10,&
&   'Action: check your pseudopotential file.'
   MSG_ERROR(msg)
 end if
 LIBPAW_ALLOCATE(ncore,(core_mesh%mesh_size))
 shft=mesh_shift(icoremesh)
 ncore(1+shft:core_mesh%mesh_size)=paw_setuploc%ae_core_density%data(1:core_mesh%mesh_size-shft)/sqrt(fourpi)
 if (shft==1) call pawrad_deducer0(ncore,core_mesh%mesh_size,core_mesh)

!Construct and save VH[z_NC] if requested
 if (pawtab%has_vhnzc==1) then
   LIBPAW_ALLOCATE(pawtab%VHnZC,(pawtab%mesh_size))
   LIBPAW_ALLOCATE(vhnzc,(core_mesh%mesh_size))
   call atompaw_vhnzc(ncore,core_mesh,vhnzc,znucl)
   pawtab%VHnZC(1:pawtab%mesh_size)=vhnzc(1:pawtab%mesh_size)
   pawtab%has_vhnzc=2
   LIBPAW_DEALLOCATE(vhnzc)
 end if

 pawtab%core_mesh_size=pawtab%mesh_size
 if(save_core_msz) pawtab%core_mesh_size=core_mesh%mesh_size
 LIBPAW_ALLOCATE(pawtab%coredens,(pawtab%core_mesh_size))
 pawtab%rcore=core_mesh%rad(pawtab%core_mesh_size)
 pawtab%coredens(1:pawtab%core_mesh_size)=ncore(1:pawtab%core_mesh_size)

!---------------------------------
!Read pseudo core density (tcoredens)

 do imsh=1,nmesh
   if(trim(paw_setuploc%pseudo_core_density%grid)==trim(paw_setuploc%radial_grid(imsh)%id)) then
     iread1=imsh
     exit
   end if
 end do
 if (iread1/=icoremesh) then
   write(msg, '(a,a,a,a,a,a,a,a)' )&
&   'Pseudized core density (tNcore) must be given',ch10,&
&   'on the same radial mesh as core density (Ncore) !',ch10,&
&   'Action: check your pseudopotential file.'
   MSG_ERROR(msg)
 end if
 LIBPAW_ALLOCATE(tncore,(core_mesh%mesh_size))
 shft=mesh_shift(icoremesh)
 tncore(1+shft:core_mesh%mesh_size)=paw_setuploc%pseudo_core_density%data(1:core_mesh%mesh_size-shft)/sqrt(fourpi)
 if (shft==1) call pawrad_deducer0(tncore,core_mesh%mesh_size,core_mesh)
 if(save_core_msz)  then
   LIBPAW_ALLOCATE(pawtab%tcoredens,(pawtab%core_mesh_size,6))
 else
   LIBPAW_ALLOCATE(pawtab%tcoredens,(pawtab%core_mesh_size,1))
 end if
 if (maxval(abs(tncore(:)))<tol6) then
   pawtab%usetcore=0
   pawtab%tcoredens(1:pawtab%core_mesh_size,:)=zero
 else
   pawtab%usetcore=1
   pawtab%tcoredens(1:pawtab%core_mesh_size,1)=tncore(1:pawtab%core_mesh_size)
 end if
 write(msg,'(a,i1)') &
& ' Radial grid used for (t)core density is grid ',icoremesh
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

!---------------------------------
!Read core kinetic density (coretau)

 if (paw_setuploc%ae_core_kinetic_energy_density%tread.and.pawtab%has_coretau>=1) then
   do imsh=1,nmesh
     if(trim(paw_setuploc%ae_core_kinetic_energy_density%grid)==trim(paw_setuploc%radial_grid(imsh)%id)) then
       icoretaumesh=imsh
       exit
     end if
   end do
   call pawrad_copy(radmesh(icoretaumesh),coretau_mesh)
   if (icoretaumesh/=icoremesh) then
     write(msg, '(5a)' )&
&     'Core kinetic density (TAUcore) must be given',ch10,&
&     'on the same radial mesh as core density (Ncore) !',ch10,&
&     'Action: check your pseudopotential file.'
     MSG_ERROR(msg)
   end if
   LIBPAW_ALLOCATE(coretau,(coretau_mesh%mesh_size))
   shft=mesh_shift(icoretaumesh)
   coretau(1+shft:coretau_mesh%mesh_size)= &
&   paw_setuploc%ae_core_kinetic_energy_density%data(1:coretau_mesh%mesh_size-shft)/sqrt(fourpi)
   if (shft==1) call pawrad_deducer0(coretau,coretau_mesh%mesh_size,coretau_mesh)
   pawtab%coretau_mesh_size=pawtab%mesh_size
   if(save_core_msz) pawtab%coretau_mesh_size=coretau_mesh%mesh_size
   LIBPAW_ALLOCATE(pawtab%coretau,(pawtab%coretau_mesh_size))
   pawtab%rcoretau=coretau_mesh%rad(pawtab%coretau_mesh_size)
   pawtab%coretau(1:pawtab%coretau_mesh_size)=coretau(1:pawtab%coretau_mesh_size)
 else if (pawtab%has_coretau>=1) then
   write(msg, '(5a)' )&
&   'metaGGA exchange-correlation is requested but the core kinetic energy density',ch10,&
&   'is not present in the pseudopotential file!',ch10,&
&   'Action: check your pseudopotential file.'
   MSG_ERROR(msg)
 end if

!---------------------------------
!Read pseudo core kinetic energy density (tcoretau)

 if (paw_setuploc%pseudo_core_kinetic_energy_density%tread.and.pawtab%has_coretau>=1) then
   do imsh=1,nmesh
     if(trim(paw_setuploc%pseudo_core_kinetic_energy_density%grid)==trim(paw_setuploc%radial_grid(imsh)%id)) then
       iread1=imsh
       exit
     end if
   end do
   if (iread1/=icoretaumesh) then
     write(msg, '(5a)' )&
&     'Pseudized core kinetic energy density (tTAUcore) must be given',ch10,&
&     'on the same radial mesh as core kinetic density (TAUcore) !',ch10,&
&     'Action: check your pseudopotential file.'
     MSG_ERROR(msg)
   end if
   LIBPAW_ALLOCATE(tcoretau,(coretau_mesh%mesh_size))
   shft=mesh_shift(icoretaumesh)
   tcoretau(1+shft:coretau_mesh%mesh_size)= &
&   paw_setuploc%pseudo_core_kinetic_energy_density%data(1:coretau_mesh%mesh_size-shft)/sqrt(fourpi)
   if (shft==1) call pawrad_deducer0(tcoretau,coretau_mesh%mesh_size,coretau_mesh)
   LIBPAW_ALLOCATE(pawtab%tcoretau,(pawtab%coretau_mesh_size))
   pawtab%tcoretau(1:pawtab%coretau_mesh_size)=tcoretau(1:pawtab%coretau_mesh_size)
   pawtab%has_coretau=2
   write(msg,'(a,i1)') &
&   ' Radial grid used for (t)coretau kinetic density is grid ',icoretaumesh
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 else if (pawtab%has_coretau>=1) then
   write(msg, '(5a)' )&
&   'metaGGA exchange-correlation is requested but the pseudo core kinetic energy density',ch10,&
&   'is not present in the pseudopotential file!',ch10,&
&   'Action: check your pseudopotential file.'
   MSG_ERROR(msg)
 end if

!---------------------------------
!Read local pseudopotential=Vh(tn_zc) or Vbare

 if ((paw_setuploc%blochl_local_ionic_potential%tread).and.&
& (pawtab%usexcnhat==-1.or.pawtab%usexcnhat==0.or.(pawtab%usexcnhat==1.and.&
& ((.not.paw_setuploc%zero_potential%tread).or.(.not.paw_setuploc%kresse_joubert_local_ionic_potential%tread))))) then
   usexcnhat=0;vlocopt=2
   do imsh=1,nmesh
     if(trim(paw_setuploc%blochl_local_ionic_potential%grid)==trim(paw_setuploc%radial_grid(imsh)%id)) then
       iread1=imsh
       exit
     end if
   end do
   ivlocmesh=iread1
   call pawrad_copy(radmesh(ivlocmesh),vloc_mesh)
   LIBPAW_ALLOCATE(vlocr,(vloc_mesh%mesh_size))
   shft=mesh_shift(ivlocmesh)
   vlocr(1+shft:vloc_mesh%mesh_size)=paw_setuploc%blochl_local_ionic_potential%data(1:vloc_mesh%mesh_size-shft)/sqrt(fourpi)
   if (shft==1) call pawrad_deducer0(vlocr,vloc_mesh%mesh_size,vloc_mesh)
 else if((paw_setuploc%kresse_joubert_local_ionic_potential%tread).and.&
&   (pawtab%usexcnhat==-1.or.pawtab%usexcnhat==1.or.(pawtab%usexcnhat==0.and.&
&   (.not.paw_setuploc%zero_potential%tread)))) then
   usexcnhat=1;vlocopt=1
   do imsh=1,nmesh
     if(trim(paw_setuploc%kresse_joubert_local_ionic_potential%grid)==trim(paw_setuploc%radial_grid(imsh)%id)) then
       iread1=imsh
       exit
     end if
   end do
   ivlocmesh=iread1
   call pawrad_copy(radmesh(ivlocmesh),vloc_mesh)
   LIBPAW_ALLOCATE(vlocr,(vloc_mesh%mesh_size))
   shft=mesh_shift(ivlocmesh)
   vlocr(1+shft:vloc_mesh%mesh_size)= &
&   paw_setuploc%kresse_joubert_local_ionic_potential%data(1:vloc_mesh%mesh_size-shft)/sqrt(fourpi)
   if (shft==1) call pawrad_deducer0(vlocr,vloc_mesh%mesh_size,vloc_mesh)
 else if(paw_setuploc%zero_potential%tread) then
   usexcnhat=0;vlocopt=0
   do imsh=1,nmesh
     if(trim(paw_setuploc%zero_potential%grid)==trim(paw_setuploc%radial_grid(imsh)%id)) then
       iread1=imsh
       exit
     end if
   end do
   ivlocmesh=iread1
!   vloc_mesh%mesh_type=radmesh(ivlocmesh)%mesh_type
!   vloc_mesh%rstep=radmesh(ivlocmesh)%rstep
!   vloc_mesh%lstep=radmesh(ivlocmesh)%lstep
!   vloc_mesh%mesh_size=radmesh(ivlocmesh)%mesh_size
!   vloc_mesh%mesh_size=pawrad_ifromr(radmesh(ivlocmesh),rmax_vloc)
   call pawrad_copy(radmesh(ivlocmesh),vloc_mesh)
   LIBPAW_ALLOCATE(vlocr,(vloc_mesh%mesh_size))
   vlocr=zero
   shft=mesh_shift(ivlocmesh)
   vlocr(1+shft:vloc_mesh%mesh_size)=paw_setuploc%zero_potential%data(1:vloc_mesh%mesh_size-shft)/sqrt(fourpi)
   if (shft==1) call pawrad_deducer0(vlocr,vloc_mesh%mesh_size,vloc_mesh)
 else
   write(msg, '(a,a,a,a,a)' )&
&   'At least one local potential must be given',ch10,&
&   'Action: check your pseudopotential file.'
   MSG_ERROR(msg)
 end if

 write(msg,'(a,i1)') &
& ' Radial grid used for Vloc is grid ',ivlocmesh
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

!-------------------------------------------------
!Read LDA-1/2 potential

 if (paw_setuploc%LDA_minus_half_potential%tread) then
   do imsh=1,nmesh
     if(trim(paw_setuploc%LDA_minus_half_potential%grid)==trim(paw_setuploc%radial_grid(imsh)%id)) then
       iread1=imsh
       exit
     end if
   end do
   if(iread1/=ivlocmesh) then
     write(msg, '(a)' )&
&     'The LDA-1/2 potential must be given on the same grid as the local potential.'
     MSG_ERROR(msg)
   end if
   has_v_minushalf=1
   LIBPAW_ALLOCATE(pawtab%vminushalf,(vloc_mesh%mesh_size))
   shft=mesh_shift(ivlocmesh)
   pawtab%vminus_mesh_size=vloc_mesh%mesh_size
   pawtab%vminushalf(1+shft:vloc_mesh%mesh_size)= &
&   paw_setuploc%LDA_minus_half_potential%data(1:vloc_mesh%mesh_size-shft)/sqrt(fourpi)
   if (shft==1) call pawrad_deducer0(pawtab%vminushalf,vloc_mesh%mesh_size,vloc_mesh)
   write(msg,'(a,i1)') &
&   ' Radial grid used for LDA-1/2 potential is grid ',ivlocmesh
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 else
   has_v_minushalf=0
 end if
 if(has_v_minushalf==0.and.pawtab%has_vminushalf==1) then
   write(msg, '(a)' )&
&     'The LDA-1/2 potential must be given in the XML PAW datafile.'
   MSG_ERROR(msg)
 end if

!---------------------------------
!Eventually read "numeric" shapefunctions (if shape_type=-1)

 if (pawtab%shape_type==-1) then
   LIBPAW_ALLOCATE(pawtab%shapefunc,(pawtab%mesh_size,pawtab%l_size))
   do imsh=1,nmesh
     if(trim(paw_setuploc%shape_function%grid)==trim(paw_setuploc%radial_grid(imsh)%id)) then
       iread1=imsh
       exit
     end if
   end do
   call pawrad_copy(radmesh(iread1),shpf_mesh)
   ishpfmesh=iread1
   LIBPAW_ALLOCATE(shpf,(shpf_mesh%mesh_size,pawtab%l_size))
   shft=mesh_shift(ishpfmesh)
   shpf(1,1)=one
   do ir=2,shpf_mesh%mesh_size
     shpf(ir,1)=paw_setuploc%shape_function%data(ir-shft,1)
   end do
   sz10=size(paw_setuploc%shape_function%data,2)
   if(sz10>=2) then
     do il=2,pawtab%l_size
       shpf(1,il)=zero
       do ir=2,shpf_mesh%mesh_size
         shpf(ir,il)=paw_setuploc%shape_function%data(ir-shft,il)
       end do
     end do
   else
     do il=2,pawtab%l_size
       shpf(1,il)=zero
       do ir=2,shpf_mesh%mesh_size
         shpf(ir,il)=paw_setuploc%shape_function%data(ir-shft,1)*shpf_mesh%rad(ir)**(il-1)
       end do
     end do
   end if
   write(msg,'(a,i1)') &
&   ' Radial grid used for shape functions is grid ',iread1
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')

!  Has to spline shape functions if mesh is not the "main" mesh
   if (ishpfmesh/=imainmesh) then
     msz=shpf_mesh%mesh_size
     LIBPAW_ALLOCATE(work1,(msz))
     LIBPAW_ALLOCATE(work2,(msz))
     LIBPAW_ALLOCATE(work3,(msz))
     LIBPAW_ALLOCATE(work4,(pawrad%mesh_size))
     work3(1:msz)=shpf_mesh%rad(1:msz)
     work4(1:pawrad%mesh_size)=pawrad%rad(1:pawrad%mesh_size)
     do il=1,pawtab%l_size
       call bound_deriv(shpf(1:msz,il),shpf_mesh,msz,yp1,ypn)
       call paw_spline(work3,shpf(:,il),msz,yp1,ypn,work1)
       call paw_splint(msz,work3,shpf(:,il),work1,pawrad%mesh_size,work4,pawtab%shapefunc(:,il))
     end do
     LIBPAW_DEALLOCATE(work1)
     LIBPAW_DEALLOCATE(work2)
     LIBPAW_DEALLOCATE(work3)
     LIBPAW_DEALLOCATE(work4)
   else
     pawtab%shapefunc(:,:)=shpf(:,:)
   end if
   LIBPAW_DEALLOCATE(shpf)
 end if

!---------------------------------
!Read pseudo valence density

 if (paw_setuploc%pseudo_valence_density%tread) then
   do imsh=1,nmesh
     if(trim(paw_setuploc%pseudo_valence_density%grid)==trim(paw_setuploc%radial_grid(imsh)%id)) then
       iread1=imsh
       exit
     end if
   end do
   ivalemesh=iread1
   call pawrad_copy(radmesh(iread1),vale_mesh)
   LIBPAW_ALLOCATE(tnvale,(vale_mesh%mesh_size))
   shft=mesh_shift(ivalemesh)
   tnvale(1+shft:vale_mesh%mesh_size)=paw_setuploc%pseudo_valence_density%data(1:vale_mesh%mesh_size-shft)/sqrt(fourpi)
   if (shft==1) call pawrad_deducer0(tnvale,vale_mesh%mesh_size,vale_mesh)
   pawtab%has_tvale=1
   write(msg,'(a,i1)') &
&   ' Radial grid used for pseudo valence density is grid ',ivalemesh
   call wrtout(ab_out,msg,'COLL')
   call wrtout(std_out,  msg,'COLL')
 else
   pawtab%has_tvale=0
   LIBPAW_ALLOCATE(tnvale,(0))
 end if

!---------------------------------
!Read initial guess of rhoij (rhoij0)

 LIBPAW_ALLOCATE(pawtab%rhoij0,(pawtab%lmn2_size))
 pawtab%rhoij0=zero
 ilmn0=0
 do ib=1,pawtab%basis_size
   il=2*pawtab%orbitals(ib)+1
   occ=paw_setuploc%valence_states%state(ib)%ff
   if (occ<zero)occ=zero
   do ilmn=ilmn0+1,ilmn0+il
     pawtab%rhoij0(ilmn*(ilmn+1)/2)=occ/dble(il)
   end do
   ilmn0=ilmn0+il
 end do

!---------------------------------
!Read Kij terms (kij0) and deduce eventually Dij0

 LIBPAW_ALLOCATE(kij,(pawtab%lmn2_size))
 kij=zero
 nval=paw_setuploc%valence_states%nval
 do jlmn=1,pawtab%lmn_size
   j0lmn=jlmn*(jlmn-1)/2
   jlm=pawtab%indlmn(4,jlmn);jln=pawtab%indlmn(5,jlmn)
   do ilmn=1,jlmn
     klmn=j0lmn+ilmn
     ilm=pawtab%indlmn(4,ilmn);iln=pawtab%indlmn(5,ilmn)
     if (ilm==jlm) kij(klmn)=paw_setuploc%kinetic_energy_differences%data(jln+(iln-1)*nval)
   end do
 end do
 if (vlocopt>0) then
   LIBPAW_ALLOCATE(pawtab%dij0,(pawtab%lmn2_size))
   if (size(pawtab%vminushalf)>0.and.pawtab%has_vminushalf==1) then
     vlocr(1:vloc_mesh%mesh_size)=vlocr(1:vloc_mesh%mesh_size)+pawtab%vminushalf(1:vloc_mesh%mesh_size)
   end if
   call atompaw_dij0(pawtab%indlmn,kij,pawtab%lmn_size,ncore,0,pawtab,pawrad,core_mesh,&
&                    vloc_mesh,vlocr,znucl)
 end if

!Keep eventualy Kij in memory
 if (pawtab%has_kij==1.or.vlocopt==0) then
   LIBPAW_ALLOCATE(pawtab%kij,(pawtab%lmn2_size))
   pawtab%kij(:)=kij(:)
   if (vlocopt> 0) pawtab%has_kij=2
!  This -1 means that pawtab%kij will be freed later
   if (vlocopt==0) pawtab%has_kij=-1
 end if

 LIBPAW_DEALLOCATE(kij)

!---------------------------------
!Read exact-exchange Fock terms for core-valence interactions (ex_cvij)

 if (paw_setuploc%exact_exchange_matrix%tread.eqv..true.) then
   pawtab%has_fock=2
   LIBPAW_ALLOCATE(pawtab%ex_cvij,(pawtab%lmn2_size))
   pawtab%ex_cvij=zero
   nval=paw_setuploc%valence_states%nval
   do jlmn=1,pawtab%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     jlm=pawtab%indlmn(4,jlmn);jln=pawtab%indlmn(5,jlmn)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       ilm=pawtab%indlmn(4,ilmn);iln=pawtab%indlmn(5,ilmn)
       if (ilm==jlm) pawtab%ex_cvij(klmn)=paw_setuploc%exact_exchange_matrix%data(jln+(iln-1)*nval)
     end do
   end do
   pawtab%ex_cc=paw_setuploc%ex_cc
 end if

!==========================================================
!Compute additional atomic data only depending on present DATASET

 call pawpsp_calc(core_mesh,epsatm,ffspl,imainmesh,ixc,lnmax,&
&     mmax,mqgrid_ff,mqgrid_vl,ncore,nmesh,pawrad,pawtab,pawxcdev,pspversion,&
&     qgrid_ff,qgrid_vl,radmesh,tncore,tnvale,tproj,tproj_mesh,usexcnhat,vale_mesh,&
&     vloc_mesh,vlocopt,vlocr,vlspl,xcccrc,xclevel,xc_denpos,zion,znucl,&
&     tcoretau=tcoretau,coretau_mesh=coretau_mesh)

 if(usewvl==1 .or. icoulomb > 0) then
!  Calculate up to the 5th derivative of tcoredens
   call pawpsp_calc_d5(core_mesh,pawtab%core_mesh_size,pawtab%tcoredens)
!  Other wvl related operations
   call pawpsp_wvl_calc(pawtab,tnvale,usewvl,vale_mesh,vloc_mesh,vlocr)
 else if (pawtab%has_wvl>0) then
   call wvlpaw_rholoc_nullify(pawtab%wvl%rholoc)
 end if

!==========================================================
!Free temporary allocated space

 call pawrad_free(radmesh)
 LIBPAW_DATATYPE_DEALLOCATE(radmesh)
 LIBPAW_DEALLOCATE(mesh_shift)

 call pawrad_free(tproj_mesh)
 call pawrad_free(core_mesh)
 call pawrad_free(vloc_mesh)

 if (allocated(vlocr)) then
   LIBPAW_DEALLOCATE(vlocr)
 end if
 if (allocated(ncore)) then
   LIBPAW_DEALLOCATE(ncore)
 end if
 if (allocated(tncore)) then
   LIBPAW_DEALLOCATE(tncore)
 end if
 if (allocated(coretau)) then
   LIBPAW_DEALLOCATE(coretau)
 end if
 if (allocated(tcoretau)) then
   LIBPAW_DEALLOCATE(tcoretau)
 end if
 if (allocated(tproj)) then
   LIBPAW_DEALLOCATE(tproj)
 end if

 if(pawtab%shape_type==-1) then
   call pawrad_free(shpf_mesh)
 end if
 if (paw_setuploc%pseudo_valence_density%tread) then
   call pawrad_free(vale_mesh)
 end if
 if (allocated(tnvale)) then
   LIBPAW_DEALLOCATE(tnvale)
 end if


end subroutine pawpsp_17in
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_7in
!! NAME
!! pawpsp_7in
!!
!! FUNCTION
!! Initialize pspcod=7 ("PAW pseudopotentials"):
!! continue to read the corresponding file and compute the form factors
!!
!! INPUTS
!!  icoulomb==0 : usual reciprocal space computation
!!           =1 : free boundary conditions are used
!!  ipsp=id in the array of the currently read pseudo.
!!  ixc=exchange-correlation choice from main routine data file
!!  lloc=angular momentum choice of local pseudopotential
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  pawxcdev=choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  xclevel= XC functional level
!!  zion=nominal valence of atom as specified in psp file
!!
!! OUTPUT
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$(hartree)
!!  ffspl(psps%mqgrid_ff,2,psps%lnmax)=form factor f_l(q) and second derivative
!!   from spline fit for each angular momentum and each projector;
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  vlspl(mqgrid_vl,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  xc_denpos= lowest allowed density (usually for the computation of the XC functionals)
!!  xcccrc=XC core correction cutoff radius (bohr) from psp file
!!
!! NOTES
!!  Spin-orbit not yet implemented (to be done)
!!
!! PARENTS
!!      m_pawpsp,m_pspini
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_7in(epsatm,ffspl,icoulomb,ixc,&
& lmax,lnmax,mmax,mqgrid_ff,mqgrid_vl,&
& pawrad,pawtab,pawxcdev,qgrid_ff,qgrid_vl,&
& usewvl,usexcnhat_in,vlspl,xcccrc,xclevel,xc_denpos,zion,znucl)

!Arguments ------------------------------------
!scalars
 integer, intent(in):: icoulomb,ixc
 integer, intent(in):: lmax,lnmax,mmax
 integer, intent(in):: mqgrid_ff,mqgrid_vl,pawxcdev
 integer, intent(in):: usewvl,usexcnhat_in,xclevel
 real(dp), intent(in):: xc_denpos,zion,znucl
 real(dp), intent(out):: epsatm,xcccrc
 type(pawrad_type), intent(inout):: pawrad
 type(pawtab_type), intent(inout) :: pawtab
!arrays
 real(dp),intent(in):: qgrid_ff(mqgrid_ff),qgrid_vl(mqgrid_vl)
 real(dp),intent(inout) :: ffspl(mqgrid_ff,2,lnmax)
 real(dp),intent(out) :: vlspl(mqgrid_vl,2)

!Local variables ------------------------------
!scalars
 integer :: imainmesh,nmesh
 integer :: pspversion,usexcnhat,vlocopt
 logical :: save_core_msz
 type(pawrad_type) :: core_mesh,tproj_mesh,vale_mesh,vloc_mesh
!arrays
 real(dp),pointer :: ncore(:),tncore(:),tcoretau(:),tnvale(:),tproj(:,:),vlocr(:)
 type(pawrad_type),pointer :: radmesh(:)

!************************************************************************

!Destroy everything in pawtab but optional flags
 call pawtab_free(pawtab)
!Destroy everything in pawrad
 call pawrad_free(pawrad)

 save_core_msz=(usewvl==1 .or. icoulomb .ne. 0)
 nullify(ncore);nullify(tncore);nullify(tcoretau);nullify(tnvale)
 nullify(tproj);nullify(vlocr)
 nullify(radmesh)

 call pawpsp_read(core_mesh,tmp_unit,imainmesh,lmax,&
&  ncore,nmesh,pawrad,pawtab,pspversion,radmesh,save_core_msz,&
&  tcoretau,tncore,tnvale,tproj,tproj_mesh,usexcnhat_in,usexcnhat,&
&  vale_mesh,vlocopt,vlocr,vloc_mesh,znucl)

 call pawpsp_calc(core_mesh,epsatm,ffspl,imainmesh,ixc,lnmax,&
&     mmax,mqgrid_ff,mqgrid_vl,ncore,nmesh,pawrad,pawtab,pawxcdev,pspversion,&
&     qgrid_ff,qgrid_vl,radmesh,tncore,tnvale,tproj,tproj_mesh,usexcnhat,vale_mesh,&
&     vloc_mesh,vlocopt,vlocr,vlspl,xcccrc,xclevel,xc_denpos,zion,znucl,&
&     tcoretau=tcoretau,coretau_mesh=core_mesh)

 if(usewvl==1 .or. icoulomb > 0) then
!  Calculate up to the 5th derivative of tcoredens
   call pawpsp_calc_d5(core_mesh,pawtab%core_mesh_size,pawtab%tcoredens)
!  Other wvl related operations
   call pawpsp_wvl_calc(pawtab,tnvale,usewvl,vale_mesh,vloc_mesh,vlocr)
 else if (pawtab%has_wvl>0) then
   call wvlpaw_rholoc_nullify(pawtab%wvl%rholoc)
 end if

!==========================================================
!Free temporary allocated space
 call pawrad_free(radmesh)
 call pawrad_free(tproj_mesh)
 call pawrad_free(core_mesh)
 call pawrad_free(vloc_mesh)
 LIBPAW_DATATYPE_DEALLOCATE(radmesh)
 if (associated(vlocr)) then
   LIBPAW_POINTER_DEALLOCATE(vlocr)
 end if
 if (associated(ncore)) then
   LIBPAW_POINTER_DEALLOCATE(ncore)
 end if
 if (associated(tncore)) then
   LIBPAW_POINTER_DEALLOCATE(tncore)
 end if
 if (associated(tnvale)) then
   LIBPAW_POINTER_DEALLOCATE(tnvale)
 end if
 if (associated(tcoretau)) then
   LIBPAW_POINTER_DEALLOCATE(tcoretau)
 end if
 if (associated(tproj)) then
   LIBPAW_POINTER_DEALLOCATE(tproj)
 end if
 if (pspversion>=4)  then
   call pawrad_free(vale_mesh)
 end if

end subroutine pawpsp_7in
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_wvl_sin2gauss
!! NAME
!!  pawpsp_wvl_sin2gauss
!!
!! FUNCTION
!!  Converts a f(x)=sum_i^N_i a_i sin(b_i x)+ c_i cos( d_i x) to
!!    f(x)=sum_j e_j exp(f_j x), where e and f are complex numbers.
!!
!! INPUTS
!! basis_size =  size of the lmn basis
!! mparam = number of terms in the summatory (N_i, see the expression above)
!! nparam = Array containing the parameters (a_i, b_i,c_i,d_i)
!! wvl = wavelets data type
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! On output wvl%pfac and wvl%parg are filled with complex parameters (e_i, f_i)
!!
!! NOTES
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!
!! SOURCE

 subroutine pawpsp_wvl_sin2gauss(basis_size,mparam,nparam,&
& param,wvl)

!Arguments ------------------------------------
  integer,intent(in) :: mparam,basis_size
  integer,intent(in) :: nparam(basis_size)
  real(dp),intent(in) :: param(mparam,basis_size)
  type(wvlpaw_type),intent(inout):: wvl

!Local variables ------------------------------
  integer :: i,ii,ib,ngauss,nterm
  real(dp) :: sep
  real(dp) :: a1(mparam),a2(mparam),a3(mparam),a4(mparam),a5(mparam)
  real(dp) :: b1r(mparam),b2r(mparam),b1i(mparam),b2i(mparam)
 character(len=500) :: message
  !
  !extra variables, use to debug
  !
  !integer::igau,nr,unitp
  !real(dp)::step,rmax
  !real(dp),allocatable::r(:), y(:)
  !complex::fac,arg
  !complex(dp),allocatable::f(:)
!************************************************************************

!  Convert from \sum(sin+cos) expressions to sums of complex gaussians
!  (only works for option=4, see fit_gen)

!  get number of coefficients:
   ii=0
   do ib=1,basis_size
     nterm=nparam(ib)/4  !option=4, there are 4 parameters for each term
     ii=ii+nterm*2 !two gaussians for each term
   end do
!
!  Allocate objects
!
   ngauss=ii
   wvl%ptotgau=ngauss !total number of complex gaussians
   LIBPAW_ALLOCATE(wvl%pfac,(2,ngauss))
   LIBPAW_ALLOCATE(wvl%parg,(2,ngauss))
   LIBPAW_ALLOCATE(wvl%pngau,(basis_size))
   wvl%pngau(1:basis_size)=nparam(1:basis_size)/2 !option=4
!
   ii=0
   do ib=1,basis_size
!
!    Get parameters in sin+cos expansion:
!    Option4: \sum a1 exp(-a2 x^2) ( a3 sin(k x^2) + a4 cos(k x^2))
!
     nterm=nparam(ib)/4  !option=4
!
     a1(1:nterm)=param(1:nterm,ib)
     a2(1:nterm)=param(nterm+1:nterm*2,ib)
     a3(1:nterm)=param(nterm*2+1:nterm*3,ib)
     a4(1:nterm)=param(nterm*3+1:nterm*4,ib)
     sep=1.1d0
     do i=1,nterm
       a5(i)=sep**(i)
     end do

!    First check that "a2" is a positive number (it is multiplied by -1, so
!    that gaussians decay to zero:
     if( any(a2(1:nterm) < tol12) ) then
       message = 'Real part of Gaussians should be a negative number (they should go to zero at infty)'
       MSG_ERROR(message)
     end if

!
!    Now translate them to a sum of complex gaussians:
!    pngau(ib)=nterm*2
!    Two gaussians by term:
!
!    First gaussian
     b1r(1:nterm)= a1(1:nterm)*a4(1:nterm)/2.d0 !coefficient, real
     b1i(1:nterm)=-a1(1:nterm)*a3(1:nterm)/2.d0 !coefficient, imag
     b2r(1:nterm)=-a2(1:nterm)   !exponential, real
     b2i(1:nterm)= a5(1:nterm)   !exponential, imag
!
     wvl%pfac(1,ii+1:ii+nterm)=b1r(1:nterm)
     wvl%pfac(2,ii+1:ii+nterm)=b1i(1:nterm)
     wvl%parg(1,ii+1:ii+nterm)=b2r(1:nterm)
     wvl%parg(2,ii+1:ii+nterm)=b2i(1:nterm)
!    Second gaussian
     wvl%pfac(1,ii+nterm+1:ii+nterm*2)= b1r(1:nterm)
     wvl%pfac(2,ii+nterm+1:ii+nterm*2)=-b1i(1:nterm)
     wvl%parg(1,ii+nterm+1:ii+nterm*2)= b2r(1:nterm)
     wvl%parg(2,ii+nterm+1:ii+nterm*2)=-b2i(1:nterm)
!
     ii=ii+nterm*2
   end do

!  begin debug
!  write(*,*)'pawpsp_wvl_sin2gauss, comment me'
!  nr=3000
!  rmax=10.d0
!  LIBPAW_ALLOCATE(r,(nr))
!  LIBPAW_ALLOCATE(f,(nr))
!  LIBPAW_ALLOCATE(y,(nr))
!  step=rmax/real(nr-1,dp)
!  do ir=1,nr
!  r(ir)=real(ir-1,dp)*step
!  end do
!  !
!  ii=0
!  do ib=1,basis_size
!  unitp=500+ib
!  f(:)=czero
!  !
!  do igau=1,wvl%pngau(ib)
!  ii=ii+1
!  arg=cmplx(wvl%parg(1,ii),wvl%parg(2,ii))
!  fac=cmplx(wvl%pfac(1,ii),wvl%pfac(2,ii))
!  f(:)=f(:)+fac*exp(arg*r(:)**2)
!  end do
!  do ir=1,nr
!  write(unitp,'(3f16.7)')r(ir),real(f(ir))!,y(ir)
!  end do
!  end do
!  LIBPAW_DEALLOCATE(r)
!  LIBPAW_DEALLOCATE(f)
!  LIBPAW_DEALLOCATE(y)
!  end debug

 end subroutine pawpsp_wvl_sin2gauss
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_read_header
!! NAME
!!  pawpsp_read_header
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!
!! SOURCE


subroutine pawpsp_read_header(funit,lloc,lmax,mmax,pspcod,pspxc,r2well,zion,znucl)

!Arguments ------------------------------------
!scalars
 integer,intent(in):: funit
 integer,intent(out):: lloc,lmax,mmax,pspcod,pspxc
 real(dp),intent(out):: r2well,zion,znucl
!Local variables-------------------------------
 integer:: pspdat
 character(len=fnlen):: title
 character(len=500) :: msg

! *************************************************************************

!Read and write some description of file from first line (character data)
 read (funit,'(a)') title
 write(msg, '(a,a)' ) '- ',trim(title)
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

!Read and write more data describing psp parameters
 read (funit,*) znucl,zion,pspdat
 write(msg, '(a,f9.5,f10.5,2x,i8,t47,a)' ) &
& '-',znucl,zion,pspdat,'znucl, zion, pspdat'
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

 read (funit,*) pspcod,pspxc,lmax,lloc,mmax,r2well
 if(pspxc<0) then
   write(msg, '(i5,i8,2i5,i10,f10.5,t47,a)' ) &
&   pspcod,pspxc,lmax,lloc,mmax,r2well,&
&   'pspcod,pspxc,lmax,lloc,mmax,r2well'
 else
   write(msg, '(4i5,i10,f10.5,t47,a)' ) &
&   pspcod,pspxc,lmax,lloc,mmax,r2well,&
&   'pspcod,pspxc,lmax,lloc,mmax,r2well'
 end if
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

end subroutine pawpsp_read_header
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_read_header_2
!! NAME
!!  pawpsp_read_header_2
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Reads pspversion, basis_size and lmn_size
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!
!! SOURCE


subroutine pawpsp_read_header_2(funit,pspversion,basis_size,lmn_size)

!Arguments ------------------------------------
!scalars
 integer,intent(in):: funit
 integer,intent(out) :: pspversion,basis_size,lmn_size

!Local variables-------------------------------
 integer :: creatorid
 character(len=80) :: pspline
 character(len=500) :: msg

! *************************************************************************

!Read psp version in line 4 of the header
 pspversion=1
 read (funit,'(a80)') pspline;pspline=adjustl(pspline)
 if (pspline(1:3)=="paw".or.pspline(1:3)=="PAW") &
& read(unit=pspline(4:80),fmt=*) pspversion
 if (pspversion<1.or.pspversion>5) then
   write(msg, '(a,i2,a,a,a)' )&
&   'This version of PAW psp file (',pspversion,') is not compatible with',ch10,&
&   'current version of Abinit.'
   MSG_ERROR(msg)
 end if

 if (pspversion==1) then
   read (unit=pspline,fmt=*) basis_size,lmn_size
 else
!  Here psp file for Abinit 4.3+
   read (unit=pspline(5:80),fmt=*) creatorid
   read (funit,*) basis_size,lmn_size
 end if

end subroutine pawpsp_read_header_2
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_wvl
!! NAME
!!  pawpsp_wvl
!!
!! FUNCTION
!! WVL+PAW related operations
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_pawpsp,m_pspini
!!
!! CHILDREN
!!
!! SOURCE


subroutine pawpsp_wvl(filpsp,pawrad, pawtab,usewvl, wvl_ngauss, comm_mpi)

!Arguments------------------------------------
!scalars
 integer, optional,intent(in):: comm_mpi
 integer, intent(in):: usewvl, wvl_ngauss(2)
 character(len=fnlen),intent(in)::filpsp
 type(pawrad_type),intent(in) :: pawrad
 type(pawtab_type),intent(inout):: pawtab
!arrays

!Local variables-------------------------------
!scalars
 integer:: ii, me, mparam, nterm_bounds(2)
 type(pawrad_type)::tproj_mesh
 character(len=500) :: msg
!arrays
 integer,allocatable:: ngauss_param(:)
 real(dp),allocatable:: gauss_param(:,:)

! *************************************************************************

 me=0; if (present(comm_mpi))me=xmpi_comm_rank(comm_mpi)

!If usewvl flag is on, we must have the pawtab%wvl pointer allocated
 if (usewvl==1.and.pawtab%has_wvl==0) then
   call wvlpaw_allocate(pawtab%wvl)
   pawtab%has_wvl=1
 end if

!Fit projectors to a sum of Gaussians:
 if (usewvl ==1 .and. pawtab%wvl%ptotgau==0 ) then

   if (pawtab%has_tproj==0) then
     msg='pawtab%tproj must be allocated'
     MSG_BUG(msg)
   end if

!  1) fit projectors to gaussians
   write(msg,'(a,a)')ch10,'Fitting tproj to Gaussians'
   call wrtout(std_out,msg,'COLL')

!  See fit_gen (option==4):
   do ii=1,2
     nterm_bounds(ii)=ceiling(wvl_ngauss(ii)/2.0)
   end do
   mparam=nterm_bounds(2)*4
   LIBPAW_ALLOCATE(gauss_param,(mparam,pawtab%basis_size))
   LIBPAW_ALLOCATE(ngauss_param,(pawtab%basis_size))
!  compute tproj_mesh
   call pawrad_init(tproj_mesh,mesh_size=size(pawtab%tproj,1),&
&    mesh_type=pawrad%mesh_type,rstep=pawrad%rstep, lstep=pawrad%lstep)

   if(present(comm_mpi)) then
     call gaussfit_projector(pawtab%basis_size,mparam,&
&     ngauss_param,nterm_bounds,pawtab%orbitals,&
&     gauss_param,tproj_mesh,&
&     pawtab%rpaw,pawtab%tproj,comm_mpi)
   else
     call gaussfit_projector(pawtab%basis_size,mparam,&
&     ngauss_param,nterm_bounds,pawtab%orbitals,&
&     gauss_param,tproj_mesh,&
&     pawtab%rpaw,pawtab%tproj)
   end if
!  tproj is now as a sum of sin+cos functions,
!  convert it to a sum of complex gaussians and fill %wvl object:
   call pawpsp_wvl_sin2gauss(pawtab%basis_size,mparam,&
&   ngauss_param,gauss_param,pawtab%wvl)
   LIBPAW_DEALLOCATE(gauss_param)
   LIBPAW_DEALLOCATE(ngauss_param)

   if(me==0) then
     call pawpsp_rw_atompaw(pawtab%basis_size,filpsp,pawtab%wvl)
   end if

   pawtab%has_wvl=2

 end if

!Projectors in real space are no more needed
 call pawrad_free(tproj_mesh)
 if(allocated(pawtab%tproj)) then
   LIBPAW_DEALLOCATE(pawtab%tproj)
   pawtab%has_tproj=0
 end if

end subroutine pawpsp_wvl
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_read_header_xml
!! NAME
!!  pawpsp_read_header_xml
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!! This is done instead of: call pawpsxml2ab( psxml, pspheads,1)
!! since pspheads does not exist in PAW library.
!! should we include it to avoid the following code replica?
!! check pspheads commented out in pawpsp_17in, and routine pawpsp_read_xml_2
!!
!! PARENTS
!!      m_pawpsp,m_pspheads
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_read_header_xml(lloc,lmax,pspcod,pspxc,&
& psxml,r2well,zion,znucl)

!Arguments ------------------------------------
!scalars
 type(paw_setup_t),intent(in) :: psxml
 integer,intent(out):: lloc,lmax,pspcod,pspxc
 real(dp),intent(out):: r2well,zion,znucl
!Local variables-------------------------------
 integer :: il
#if defined LIBPAW_HAVE_LIBXC
 integer :: ii
#endif
 character(len=100) :: xclibxc
 character(len=500) :: msg
!arrays

! *************************************************************************

 lloc   = 0
 r2well = 0
 pspcod=17
 znucl=psxml%atom%znucl
 zion =psxml%atom%zval

!lmax:
 lmax = 0
 do il=1,psxml%valence_states%nval
   if(psxml%valence_states%state(il)%ll>lmax) lmax=psxml%valence_states%state(il)%ll
 end do
!pspxc
 select case(trim(psxml%xc_functional%name))
   case('PZ')
     pspxc = 2
#if defined LIBPAW_HAVE_LIBXC
     pspxc = -(libxc_functionals_getid('XC_LDA_X')*1000 &
&             +libxc_functionals_getid('XC_LDA_C_PZ'))
#endif
   case('W')
     pspxc = 4
#if defined LIBPAW_HAVE_LIBXC
     pspxc = -(libxc_functionals_getid('XC_LDA_X')*1000 &
&             +libxc_functionals_getid('XC_LDA_C_WIGNER'))
#endif
   case('HL')
     pspxc = 5
#if defined LIBPAW_HAVE_LIBXC
     pspxc = -(libxc_functionals_getid('XC_LDA_X')*1000 &
&             +libxc_functionals_getid('XC_LDA_C_HL'))
#endif
   case('GL')
#if defined LIBPAW_HAVE_LIBXC
     pspxc = -(libxc_functionals_getid('XC_LDA_X')*1000 &
&             +libxc_functionals_getid('XC_LDA_C_GL'))
#else
     write(msg, '(7a)' )&
&     'The exchange and correlation functional by Gunnarson-Lundqvist', ch10,&
&     'is not implemented in Abinit.',ch10,&
&     'Action : choose another XC functional in the pseudopotential',ch10, &
&     '         generation or compile ABINIT with the libXC library.'
     MSG_ERROR(msg)
#endif
   case('VWN')
#if defined LIBPAW_HAVE_LIBXC
     pspxc = -(libxc_functionals_getid('XC_LDA_X')*1000 &
&             +libxc_functionals_getid('XC_LDA_C_VWN'))
#else
     write(msg, '(7a)' )&
&     'The exchange and correlation functional by Vosko,Wilk and Nusair', ch10,&
&     'is not implemented in Abinit.',ch10,&
&     'Action : choose another XC functional in the pseudopotential',ch10, &
&     '         generation or compile ABINIT with the libXC library.'
     MSG_ERROR(msg)
#endif
   case('PW')
     pspxc = 7
#if defined LIBPAW_HAVE_LIBXC
     pspxc = -(libxc_functionals_getid('XC_LDA_X')*1000 &
&             +libxc_functionals_getid('XC_LDA_C_PW'))
#endif
   case('PBE')
     pspxc = 11
#if defined LIBPAW_HAVE_LIBXC
     pspxc = -(libxc_functionals_getid('XC_GGA_X_PBE')*1000 &
&             +libxc_functionals_getid('XC_GGA_C_PBE'))
#endif
   case('revPBE')
     pspxc = 14
#if defined LIBPAW_HAVE_LIBXC
     pspxc = -(libxc_functionals_getid('XC_GGA_X_PBE_R')*1000 &
&             +libxc_functionals_getid('XC_GGA_C_PBE'))
#endif
   case('RPBE')
     pspxc = 15
#if defined LIBPAW_HAVE_LIBXC
     pspxc = -(libxc_functionals_getid('XC_GGA_X_RPBE')*1000 &
&             +libxc_functionals_getid('XC_GGA_C_PBE'))
#endif
   case('PW91')
#if defined LIBPAW_HAVE_LIBXC
     pspxc = -(libxc_functionals_getid('XC_GGA_X_PW91')*1000 &
&             +libxc_functionals_getid('XC_GGA_C_PW91'))
#else
     write(msg, '(7a)' )&
&     'The exchange and correlation functional by Perdew and Wang 91', ch10,&
&     'is not implemented in Abinit.',ch10,&
&     'Action : choose another XC functional in the pseudopotential',ch10, &
&     '         generation or compile ABINIT with the libXC library.'
     MSG_ERROR(msg)
#endif
   case('BLYP')
#if defined LIBPAW_HAVE_LIBXC
     pspxc = -(libxc_functionals_getid('XC_GGA_X_B88')*1000 &
&             +libxc_functionals_getid('XC_GGA_C_LYP'))
#else
     write(msg, '(7a)' )&
&     'The exchange and correlation functional BLYP', ch10,&
&     'is not implemented in Abinit.',ch10,&
&     'Action : choose another XC functional in the pseudopotential',ch10, &
&     '         generation or compile ABINIT with the libXC library.'
     MSG_ERROR(msg)
#endif
   case DEFAULT
     xclibxc=trim(psxml%xc_functional%name)
     if (xclibxc(1:3)=='XC_'  .or.xclibxc(1:3)=='xc_'  .or. &
&        xclibxc(1:5)=='LDA_X'.or.xclibxc(1:5)=='LDA_C'.or. &
&        xclibxc(1:5)=='lda_x'.or.xclibxc(1:5)=='lda_c'.or. &
&        xclibxc(1:5)=='GGA_X'.or.xclibxc(1:5)=='GGA_C'.or. &
&        xclibxc(1:5)=='gga_x'.or.xclibxc(1:5)=='gga_c') then
#if defined LIBPAW_HAVE_LIBXC
       ii=index(xclibxc,'+')
       if (ii>0) then
         pspxc=-(libxc_functionals_getid(xclibxc(1:ii-1))*1000 &
&               +libxc_functionals_getid(xclibxc(ii+1:)))
       else
         pspxc=-libxc_functionals_getid(xclibxc)
       end if
#else
       msg='Cannot use LibXC functional because ABINIT is not compiled with LibXC !'
       MSG_ERROR(msg)
#endif
!      To be eliminated later (temporary)
     else if(trim(psxml%xc_functional%functionaltype)=='LIBXC')then
#if defined LIBPAW_HAVE_LIBXC
       xclibxc=trim(psxml%xc_functional%name)
       read(unit=xclibxc,fmt=*) pspxc
       pspxc=-pspxc
#else
       msg='Cannot use LibXC functional because ABINIT is not compiled with LibXC !'
       MSG_ERROR(msg)
#endif
     else
       write(msg, '(3a)') 'Unknown XC functional in psp file: ',trim(xclibxc),' !'
       MSG_ERROR(msg)
     end if
 end select

end subroutine pawpsp_read_header_xml
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_read_pawheader
!! NAME
!!  pawpsp_read_pawheader
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_pawpsp,m_pspheads
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_read_pawheader(basis_size,lmax,lmn_size,&
& l_size,mesh_size,pspversion,psxml,rpaw,rshp,shape_type)

!Arguments ------------------------------------
!scalars
 integer,intent(in):: lmax
 integer,intent(out):: basis_size,mesh_size,lmn_size,l_size
 integer,intent(out):: pspversion,shape_type
 real(dp),intent(out)::rpaw,rshp
 type(paw_setup_t),intent(in) :: psxml
!Local variables-------------------------------
 integer::il

! *************************************************************************

!All of this was moved from pawpsxml2ab,
!basis_size
 basis_size=psxml%valence_states%nval
!mesh_size
 do il=1,psxml%ngrid
   if(psxml%radial_grid(il)%id==psxml%idgrid) &
&   mesh_size=psxml%radial_grid(il)%iend-psxml%radial_grid(il)%istart+1
 end do
!lmn_size:
 lmn_size=0
 do il=1,psxml%valence_states%nval
   lmn_size=lmn_size+2*psxml%valence_states%state(il)%ll+1
 end do
!lsize
 l_size=2*lmax+1
!pspversion
 pspversion=10
!rpaw:
 rpaw=0.d0
 if (psxml%rpaw<0.d0) then
   do il=1,psxml%valence_states%nval
     if(psxml%valence_states%state(il)%rc>rpaw) rpaw=psxml%valence_states%state(il)%rc
   end do
 else
   rpaw=psxml%rpaw
 end if
!shape_type, rshp:
 select case(trim(psxml%shape_function%gtype))
   case('gauss')
     shape_type=1
     rshp=rpaw
   case('bessel')
     shape_type=3
     rshp=psxml%shape_function%rc
   case('sinc')
     shape_type=2
     rshp=psxml%shape_function%rc
   case('exp')
     shape_type=1
     rshp=rpaw
   case('num')
     shape_type=-1
     rshp=rpaw
 end select

end subroutine pawpsp_read_pawheader
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_bcast
!! NAME
!! pawpsp_bcast
!!
!! FUNCTION
!! Communicate paw data to all processors
!!
!! INPUTS
!! comm_mpi= communicator used to broadcast data
!! lnmax= Max. number of (l,n) components over all type of psps
!! mqgrid_ff= dimension of ffspl
!! mqgrid_vl= dimension of vlspl
!!
!! OUTPUT
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r}) dr]$(hartree)
!!  ffspl(mqgrid_ff,2,lnmax)=Kleinman-Bylander form factor f_l(q) and derivative
!!  pawrad=<type pawrad_type>
!!  pawtab=<type pawtab_type>
!!  vlspl(mqgrid_vl,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  xcccrc=XC core correction cutoff radius (bohr) from psp file
!!
!! PARENTS
!!      m_pawpsp,m_pspini
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_bcast(comm_mpi,epsatm,ffspl,pawrad,pawtab,vlspl,xcccrc)

!Arguments ------------------------------------
 integer,intent(in) :: comm_mpi
 real(dp),intent(inout) :: epsatm,xcccrc
 real(dp),intent(inout) :: ffspl(:,:,:),vlspl(:,:)
 type(pawrad_type),intent(inout) :: pawrad
 type(pawtab_type),intent(inout) :: pawtab

!Local variables-------------------------------
 integer :: ierr,ii,me,nn_dpr
 integer :: siz_ffspl,siz1_ffspl,siz2_ffspl,siz3_ffspl,siz_vlspl,siz1_vlspl,siz2_vlspl
 integer,allocatable :: list_int(:)
 real(dp),allocatable :: list_dpr(:)

!*************************************************************************

 me=xmpi_comm_rank(comm_mpi)

!Broadcast pawrad
 call pawrad_bcast(pawrad,comm_mpi)

!Broadcast pawtab (only data read from file)
 call pawtab_bcast(pawtab,comm_mpi,only_from_file=.true.)

!Broadcast the sizes of the arrays
 LIBPAW_ALLOCATE(list_int,(5))
 if (me==0) then
   siz1_vlspl=size(vlspl,1); list_int(1)=siz1_vlspl
   siz2_vlspl=size(vlspl,2); list_int(2)=siz2_vlspl
   siz1_ffspl=size(ffspl,1); list_int(3)=siz1_ffspl
   siz2_ffspl=size(ffspl,2); list_int(4)=siz2_ffspl
   siz3_ffspl=size(ffspl,3); list_int(5)=siz3_ffspl
 end if
 call xmpi_bcast(list_int,0,comm_mpi,ierr)
 if (me/=0) then
   siz1_vlspl=list_int(1)
   siz2_vlspl=list_int(2)
   siz1_ffspl=list_int(3)
   siz2_ffspl=list_int(4)
   siz3_ffspl=list_int(5)
 end if
 siz_vlspl=siz1_vlspl*siz2_vlspl
 siz_ffspl=siz1_ffspl*siz2_ffspl*siz3_ffspl
 LIBPAW_DEALLOCATE(list_int)

!Broadcast the reals
 nn_dpr=2+siz_vlspl+siz_ffspl
 LIBPAW_ALLOCATE(list_dpr,(nn_dpr))
 if (me==0) then
   ii=1
   list_dpr(ii)=epsatm ;ii=ii+1
   list_dpr(ii)=xcccrc ;ii=ii+1
   list_dpr(ii:ii+siz_vlspl-1)=reshape(vlspl,(/siz_vlspl/)) ;ii=ii+siz_vlspl
   list_dpr(ii:ii+siz_ffspl-1)=reshape(ffspl,(/siz_ffspl/)) ;ii=ii+siz_ffspl
 end if
 call xmpi_bcast(list_dpr,0,comm_mpi,ierr)
 if (me/=0) then
   ii=1
   epsatm=list_dpr(ii) ;ii=ii+1
   xcccrc=list_dpr(ii) ;ii=ii+1
   vlspl=reshape(list_dpr(ii:ii+siz_vlspl-1),(/siz1_vlspl,siz2_vlspl/))
   ii=ii+siz_vlspl
   ffspl=reshape(list_dpr(ii:ii+siz_ffspl-1),(/siz1_ffspl,siz2_ffspl,siz3_ffspl/))
   ii=ii+siz_ffspl
 end if
 LIBPAW_DEALLOCATE(list_dpr)

end subroutine pawpsp_bcast
!!***

!-------------------------------------------------------------------------

!!****f* m_pawpsp/pawpsp_main
!! NAME
!! pawpsp_main
!!
!! FUNCTION
!! Reads a PAW dataset (atomic data)
!!
!! INPUTS
!!  filpsp=name of the file containing the PAW dataset
!!  usewvl=1 if we use a wavelet basis, 0 other wise (plane waves)
!!  icoulomb=1 if we use a Poisson routine with wavelets, 0 otherwise
!!  ixc=index of the XC correlation functional
!!  xclevel=type of XC functional (1=LDA, 2=GGA, ...)
!!  pawxcdev=order of the developement of the PAW on-site terms
!!           (0: full calculation, 1: order 1, 2:order 2)
!!  usexcnhat=flag controlling the use of compensation charge (nhat) in XC potential
!!  qgrid_ff=size of the mesh for the sin FFT transform of the non-local projectors (form factors)
!!           (plane waves only, 0 otherwise)
!!  qgrid_vl=size of the mesh for the sin FFT transform of the local potential
!!           (plane waves only, 0 otherwise)
!!  ffspl=sin FFT transform of the non-local projectors (form factors) (plane waves only)
!!  vlspl=sin FFT transform of the local potential (plane waves only)
!!  epsatm=$ 4\pi\int[r^2 (V(r)+\frac{Zv}{r}dr]$.
!!  xcccrc=XC core correction cutoff radius (bohr)
!!  zionpsp=valence of atom as specified in input file
!!  znuclpsp=atomic number of atom as specified in input file
!!  ===== Optional arguments for wvl =====
!!    wvl_ngauss
!!  ===== Other optional arguments =====
!!    psxml=datastructure containing a XMP PAW dataset
!!    comm_mpi=MPI communicator
!!    xc_denpos=tolerance on density/potential for the calculation of XC potential
!!              (if density<xc_denpos, density=zero)
!!
!! OUTPUT
!!  pawrad <type(pawrad_type)>=data containing PAW radial grid information
!!  pawtab <type(pawtab_type)>=data containing the PAW dataset (partial waves...)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawpsp_main( &
& pawrad,pawtab,&
& filpsp,usewvl,icoulomb,ixc,xclevel,pawxcdev,usexcnhat,&
& qgrid_ff,qgrid_vl,ffspl,vlspl,epsatm,xcccrc,zionpsp,znuclpsp,&
& wvl_ngauss,psxml,comm_mpi,xc_denpos)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icoulomb,ixc
 integer,intent(in) :: pawxcdev,usewvl,usexcnhat,xclevel
 integer,optional,intent(in) :: comm_mpi
 real(dp),intent(in):: zionpsp,znuclpsp
 real(dp),optional,intent(in) :: xc_denpos
 real(dp),intent(out) :: epsatm,xcccrc
 character(len=fnlen),intent(in):: filpsp   ! name of the psp file
 type(pawrad_type),intent(inout) :: pawrad
 type(pawtab_type),intent(inout) :: pawtab
 type(paw_setup_t),optional,intent(in) :: psxml
!arrays
 integer,optional,intent(in) :: wvl_ngauss(2)
 real(dp),intent(in) :: qgrid_ff(:),qgrid_vl(:)
 real(dp),intent(inout) :: ffspl(:,:,:)
 real(dp),intent(out) :: vlspl(:,:)

!Local variables-------------------------------
 integer :: has_coretau,has_tproj,has_wvl,ipsp,lmax,lloc,lnmax,mmax,me,mqgrid_ff,mqgrid_vl
 integer :: pspcod,pspxc,usexml
 real(dp),parameter :: xc_denpos_default=tol14
 real(dp) :: my_xc_denpos,r2well,zion,znucl
 character(len=500) :: msg
 type(pawpsp_header_type) :: pawpsp_header
!arrays

! *************************************************************************

!Check consistency of parameters
 if (icoulomb/= 0.or.usewvl==1) then
   if (.not.present(wvl_ngauss)) then
     msg='usewvl==1 or icoulomb/=0: a mandatory argument is missing!'
     MSG_BUG(msg)
   end if
 end if

 mqgrid_ff=size(qgrid_ff)
 mqgrid_vl=size(qgrid_vl)
 lnmax=size(ffspl,3)
 if (size(ffspl,1)/=mqgrid_ff.or.size(ffspl,2)/=2) then
   msg='invalid sizes for ffspl!'
   MSG_BUG(msg)
 end if
 if (size(vlspl,1)/=mqgrid_vl.or.size(vlspl,2)/=2) then
   msg='invalid sizes for vlspl!'
   MSG_BUG(msg)
 end if

 my_xc_denpos=xc_denpos_default;if (present(xc_denpos)) my_xc_denpos=xc_denpos
 pawtab%usexcnhat=usexcnhat
 me=0;if (present(comm_mpi))me=xmpi_comm_rank(comm_mpi)

 has_wvl=0; if (usewvl==1.or.icoulomb/=0) has_wvl=1
 has_tproj=0; if (usewvl==1) has_tproj=1
 has_coretau=0 ; if (pawxc_get_usekden(ixc)>=1) has_coretau=1
 call pawtab_set_flags(pawtab,has_coretau=has_coretau,has_tvale=1,has_wvl=has_wvl,has_tproj=has_tproj)

 if(me==0) then
   write(msg, '(a,t38,a)' )'- pspatm: opening atomic psp file',trim(filpsp)
   call wrtout(ab_out,  msg,'COLL')
   call wrtout(std_out,  msg,'COLL')

!  This checks if file is xml or UPF
!  It sets usexml as well
   call pawpsp_check_xml_upf(filpsp)

!  ----------------------------------------------------------------------------
   if (usexml /= 1) then
!    Open the atomic data file, and read the three first lines
     open (unit=tmp_unit,file=filpsp,form='formatted',status='old')
     rewind (unit=tmp_unit)
!    Read first 3 lines of psp file:
     call pawpsp_read_header(tmp_unit,lloc,lmax,mmax,pspcod,&
&     pspxc,r2well,zion,znucl)

   else if (usexml == 1 .and. present(psxml)) then
     write(msg,'(a,a)')  &
&     '- pawpsp : Reading pseudopotential header in XML form from ', trim(filpsp)
     call wrtout(ab_out,msg,'COLL')
     call wrtout(std_out,  msg,'COLL')

!    Return header information
     call pawpsp_read_header_xml(lloc,lmax,pspcod,&
&     pspxc,psxml,r2well,zion,znucl)
!    Fill in pawpsp_header object:
     call pawpsp_read_pawheader(pawpsp_header%basis_size,&
&   lmax,pawpsp_header%lmn_size,&
&   pawpsp_header%l_size,pawpsp_header%mesh_size,&
&   pawpsp_header%pawver,psxml,&
&   pawpsp_header%rpaw,pawpsp_header%rshp,pawpsp_header%shape_type)
   end if

!  Check data for consistency against main routine input
   call pawpsp_consistency()

!  Read rest of the PSP file
   if (pspcod==7) then
!    ABINIT proprietary format
     call pawpsp_7in(epsatm,ffspl,icoulomb,ixc,&
&     lmax,lnmax,mmax,mqgrid_ff,mqgrid_vl,&
&     pawrad,pawtab,pawxcdev,qgrid_ff,qgrid_vl,&
&     usewvl,usexcnhat,vlspl,xcccrc,xclevel,my_xc_denpos,zion,znucl)

   else if (pspcod==17)then
!    XML format
     ipsp=1
     call pawpsp_17in(epsatm,ffspl,icoulomb,ipsp,ixc,lmax,&
&     lnmax,mmax,mqgrid_ff,mqgrid_vl,pawpsp_header,pawrad,pawtab,&
&     pawxcdev,qgrid_ff,qgrid_vl,usewvl,usexcnhat,vlspl,xcccrc,&
&     xclevel,my_xc_denpos,zion,znucl)

   end if
 end if!me==0

 close(unit=tmp_unit)

 write(msg,'(3a)') ' pawpsp: atomic psp has been read ',&
& ' and splines computed',ch10
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,  msg,'COLL')

!Communicate PAW objects
 if(present(comm_mpi)) then
   if(xmpi_comm_size(comm_mpi)>1) then
     call pawpsp_bcast(comm_mpi,epsatm,ffspl,pawrad,pawtab,vlspl,xcccrc)
   end if
 end if

!WVL+PAW:
 if(icoulomb/=0.or.usewvl==1) then
   if(present(comm_mpi))then
    call pawpsp_wvl(filpsp,pawrad,pawtab,usewvl,wvl_ngauss,comm_mpi)
   else
    call pawpsp_wvl(filpsp,pawrad,pawtab,usewvl,wvl_ngauss)
   end if
 end if

contains
!!***

!-------------------------------------------------------------------------

!!****f* pawpsp_main/pawpsp_check_xml_upf
!! NAME
!!  pawpsp_main_checks
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE


subroutine pawpsp_check_xml_upf(filpsp)

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in):: filpsp   ! name of the psp file

!Local variables-------------------------------
 integer :: unt
 character(len=70):: testxml

! *************************************************************************

!  Check if the file pseudopotential file is written in XML
   usexml = 0
   unt = libpaw_get_free_unit()
   open (unit=unt,file=filpsp,form='formatted',status='old',action="read")
   rewind (unit=unt)
   read(unt,*) testxml
   if(testxml(1:5)=='<?xml')then
     usexml = 1
     read(unt,*) testxml
     if(testxml(1:4)/='<paw')then
       msg='Reading a NC pseudopotential for a PAW calculation?'
       MSG_BUG(msg)
     end if
   else
     usexml = 0
   end if
   close (unit=unt)

!  Check if pseudopotential file is a Q-espresso UPF file
   unt = libpaw_get_free_unit()
   open (unit=unt,file=filpsp,form='formatted',status='old',action="read")
   rewind (unit=unt)
   read(unt,*) testxml ! just a string, no relation to xml.
   if(testxml(1:9)=='<PP_INFO>')then
     msg='UPF format not allowed with PAW (USPP part not read yet)!'
     MSG_ERROR(msg)
   end if
   close (unit=unt)

end subroutine pawpsp_check_xml_upf
!!***

!-------------------------------------------------------------------------

!!****f* pawpsp_main/pawpsp_consistency
!! NAME
!!  pawpsp_consistency
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!
!! SOURCE


subroutine pawpsp_consistency()

! *************************************************************************

!Check pspcod=7 or 17
 if(pspcod/=7 .and. pspcod/=17)then
   write(msg, '(a,i2,a,a)' )&
&   'In reading atomic psp file, finds pspcod=',pspcod,ch10,&
&   'This is not an allowed value within PAW.'
   MSG_BUG(msg)
 end if

!Does nuclear charge znuclpsp agree with psp input znucl
 if (abs(znuclpsp-znucl)>tol8) then
   write(msg, '(a,f10.5,2a,f10.5,5a)' )&
&   'Pseudopotential file znucl=',znucl,ch10,&
&   'does not equal input znuclpsp=',znuclpsp,' better than 1e-08 .',ch10,&
&   'znucl is read from the psp file in pspatm_abinit, while',ch10,&
&   'znuclpsp is read in iofn2.'
   MSG_BUG(msg)
 end if

!Does nuclear charge zionpsp agree with psp input zion
 if (abs(zionpsp-zion)>tol8) then
   write(msg, '(a,f10.5,2a,f10.5,5a)' )&
&   'Pseudopotential file zion=',zion,ch10,&
&   'does not equal input zionpsp=',zionpsp,' better than 1e-08 .',ch10,&
&   'zion is read from the psp file in pawpsp_main, while',ch10,&
&   'zionpsp is read in iofn2.'
   MSG_BUG(msg)
 end if

!Check several choices for ixc against pspxc
!ixc is from ABINIT code; pspxc is from atomic psp file
 if (ixc==0) then
   msg='Note that input ixc=0 => no xc is being used.'
   MSG_WARNING(msg)
 else if(ixc/=pspxc) then
   write(msg, '(a,i8,a,a,a,i8,a,a,a,a,a,a,a,a,a,a)' ) &
&   'Pseudopotential file pspxc=',pspxc,',',ch10,&
&   'not equal to input ixc=',ixc,'.',ch10,&
&   'These parameters must agree to get the same xc ',ch10,&
&   'in ABINIT code as in psp construction.',ch10,&
&   'Action : check psp design or input file.',ch10,&
&   'Assume experienced user. Execution will continue.',ch10
   MSG_WARNING(msg)
 end if

 if (lloc>lmax ) then
   write(msg, '(a,2i12,a,a,a,a)' )&
&   'lloc,lmax=',lloc,lmax,ch10,&
&   'chosen l of local psp exceeds range from input data.',ch10,&
&   'Action : check pseudopotential input file.'
   MSG_ERROR(msg)
 end if

end subroutine pawpsp_consistency
!!***

end subroutine pawpsp_main
!!***

!-------------------------------------------------------------------------

end module m_pawpsp
!!***
