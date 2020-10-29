!!****m* ABINIT/m_fft_mesh
!! NAME
!!  m_fft_mesh
!!
!! FUNCTION
!!  This module contains routines and helper functions to perform the setup of the FFT mesh
!!  It also provides a set of tools to test the grid, rotate the mesh according to the symmetry
!!  operations of the space group etc.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG, XG, GMR, VO, LR, RWG, YMN, RS, TR, DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_fft_mesh

 use defs_basis
 use m_errors
 use m_abicore
 use m_hide_blas

 use defs_fftdata,     only : size_goed_fft
 use m_fstrings,       only : sjoin, itoa
 use m_numeric_tools,  only : denominator, mincm, iseven, pfactorize
 use m_symtk,          only : mati3inv
 use m_geometry,       only : xred2xcart
 use m_crystal,        only : crystal_t

 implicit none

 private

 public :: setmesh             ! Perform the setup of the FFT mesh for the GW oscillator strengths.
 public :: check_rot_fft       ! Test whether the mesh is compatible with the rotational part of the space group.
 public :: fft_check_rotrans   ! Test whether the mesh is compatible with the symmetries of the space group.
 public :: rotate_fft_mesh     ! Calculate the FFT index of the rotated mesh.
 public :: denpot_project      ! Compute n(r) + n( $R^{-1}(r-\tau)$ in) / 2
                               ! Mainly used with R = inversion to select the even/odd part under inversion
 public :: cigfft              ! Calculate the FFT index of G-G0.
 public :: ig2gfft             ! Returns the component of a G in the FFT Box from its sequential index.
 public :: g2ifft              ! Returns the index of the G in the FFT box from its reduced coordinates.
 public :: get_gftt            ! Calculate the G"s in the FFT box from ngfft
 public :: calc_ceigr          ! e^{iG.r} on the FFT mesh (complex valued).
 public :: calc_eigr           ! e^{iG.r} on the FFT mesh (version for real array with RE,IM).
 public :: calc_ceikr          ! e^{ik.r} on the FFT mesh (complex valued).
 public :: times_eigr          ! Multiply an array on the real-space mesh by e^{iG0.r}
 public :: times_eikr          ! Multiply an array on the real-space mesh by e^{ik.r}
 public :: phase               ! Compute ph(ig)=$\exp(\pi\ i \ n/ngfft)$ for n=0,...,ngfft/2,-ngfft/2+1,...,-1
 public :: mkgrid_fft          ! Sets the grid of fft (or real space) points to be treated.

 interface calc_ceigr
   module procedure calc_ceigr_spc
   module procedure calc_ceigr_dpc
 end interface calc_ceigr
!!***

!----------------------------------------------------------------------

!!****t* m_fft_mesh/zpad_t
!! NAME
!!  zpad_t
!!
!! FUNCTION
!!   Store tables used for zero-padded FFTs.
!!
!! SOURCE

 type,public :: zpad_t

   integer :: nlinex
   ! Total number of 1D transforms.

   integer :: n_zplanes
   ! Number of z-planes intersecting the sphere.

   integer,allocatable :: zplane(:,:)
   ! zplane(3,n_zplanes)
   ! zplane(1,zpl) : mapping z-plane index -> FFT index_z
   ! zplane(2,zpl) : mapping z-plane index -> igb index in array gbound

   integer,allocatable :: linex2ifft_yz(:,:)
   ! linex2ifft_yz(2,nlinex)
   ! mapping 1D-FFT -> (FFT_index_y, FFT index_z)

 end type zpad_t

 public :: zpad_init
 public :: zpad_free
!!***

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/zpad_init
!! NAME
!!  zpad_init
!!
!! FUNCTION
!!  Creation method
!!
!! INPUTS
!!   mgfft=MAX(nx,ny,nz), only used to dimension gbound
!!   gbound(2*mgfft+8,2)= The boundaries of the basis sphere of G vectors at a given k-point.
!!     See sphereboundary for more info.
!!
!! OUTPUT
!!  zpad<type(zpad_t)>
!!
!! PARENTS
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine zpad_init(zpad,nx,ny,nz,ldx,ldy,ldz,mgfft,gbound)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nx,ny,nz,ldx,ldy,ldz,mgfft
 type(zpad_t),intent(out) :: zpad
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2)

!Local variables-------------------------------
!scalars
 integer :: jj,g3_max,g3_min,gg3,ifft_g3,igb,g2min,g2max,nlinex

! *************************************************************************

 g3_min = gbound(3,2)
 g3_max = gbound(4,2)

 zpad%n_zplanes = g3_max - g3_min + 1

 ABI_MALLOC(zpad%zplane,      (2,nz))
 ABI_MALLOC(zpad%linex2ifft_yz, (2,nx*ny*nz))
 !
 ! Loop over the z-planes intersecting the G-sphere.
 nlinex = 0
 do gg3=1,zpad%n_zplanes
   !
   if (gg3<=g3_max+1) then
     ifft_g3 = gg3
   else
     ifft_g3 = gg3 + nz - zpad%n_zplanes ! Wrap around for negative gg3.
   end if
   !
   ! Select the set of y for this z-plane.
   igb=2*gg3+3
   g2min = gbound(igb  ,2)
   g2max = gbound(igb+1,2)

   zpad%zplane(1,gg3) = ifft_g3
   zpad%zplane(2,gg3) = igb

   !(1:g2max+1,ifft_g3)     ! Positive g_y.
   !(g2min+ny+1:ny,ifft_g3) ! Negative g_y.

   do jj=1,g2max+1
     nlinex = nlinex + 1
     zpad%linex2ifft_yz(1,nlinex) = jj
     zpad%linex2ifft_yz(2,nlinex) = ifft_g3
   end do

   do jj=g2min+ny+1,ny
     nlinex = nlinex + 1
     zpad%linex2ifft_yz(1,nlinex) = jj
     zpad%linex2ifft_yz(2,nlinex) = ifft_g3
   end do
 end do

 zpad%nlinex = nlinex

 RETURN
 ABI_UNUSED((/ldx,ldy,ldz/))

end subroutine zpad_init
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/zpad_free
!! NAME
!!  zpad_free
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine zpad_free(zpad)

!Arguments ------------------------------------
!scalars
 type(zpad_t),intent(inout) :: zpad

! *************************************************************************

 ABI_SFREE(zpad%zplane)
 ABI_SFREE(zpad%linex2ifft_yz)

end subroutine zpad_free
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/setmesh
!!
!! NAME
!! setmesh
!!
!! FUNCTION
!! Calculate the size of the FFT grid for the GW calculation.
!!
!! INPUTS
!!  gmet(3,3)=Reciprocal space metric.
!!  gvec(3,npwvec)=G-vectors in reduced coordinates.
!!  npwvec=Number of G vectors in the array gvec max(npwwfn,npwsigx)
!!  npwsigx=Size of the dielectric or self-energy matrix.
!!  npwwfn=Number of G-vectors in the wavefunctions.
!!  method=Integer flag for FFT grid (see below)
!!  mG0=Number of shells that must be added to take into account umklapp processes.
!!  Cryst<crystal_t>=Data type gathering information on unit cell and symmetries
!!    %nsym=Number of symmetry operations in the SG.
!!    %symrel(3,3,nsym)=Symmetry operations in real space.
!!    %tnons(3,nsym)=Fractional translations.
!!  enforce_sym=Flag to enforce a FFT which fulfils all symmetry operations, both the
!!   rotational part and fractional translations.
!!  [unit]=Output unit, defaults to std_out
!!
!! OUTPUT
!! ngfft(18)=contain all needed information about 3D FFT,
!!  see also ~abinit/doc/variables/vargs.htm#ngfft
!! nfftot= ngfft(1)*ngfft(2)*ngfft(3)=Total number of points in the FFT grid.
!!
!! NOTES
!! Four methods are implemented for the calculation of the mesh:
!!  method=0 --> FFT mesh defined by the user, useful for debugging.
!!  method=1     Roughly takes the FFT box which encloses the larger of the two spheres of radius
!!               aliasing_factor*rwfn and rsigx, where rwfn and rsigx are the radius of the spheres
!!               with npwwfn and npwsigx planewaves respectively. The default aliasing_factor is 1.
!!  method=2 --> Calculates the optimal FFT grid which allows aliasing only outside the sphere of the
!!               npwsigx planewaves (finer than method=1 with aliasing_factor=1).
!!  method=3 --> Calculates the FFT grid needed to expand the density.
!!               (even finer than method=2, roughly corresponds to method=1 with aliasing_factor=2).
!!
!!  See defs_fftdata for a list of allowed sizes of FFT.
!!
!! PARENTS
!!      m_bethe_salpeter,m_screening_driver,m_sigma_driver
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine setmesh(gmet,gvec,ngfft,npwvec,npwsigx,npwwfn,nfftot,method,mG0,Cryst,enforce_sym,unit)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enforce_sym,method,npwsigx,npwvec,npwwfn
 integer,intent(out) :: nfftot
 integer,optional,intent(in) :: unit
 type(crystal_t),target,intent(in) :: Cryst
!arrays
 integer,intent(in) :: gvec(3,npwvec),mG0(3)
 integer,intent(inout) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3)

!Local variables ------------------------------
!scalars
 integer :: aliasing_factor,fftalg,fftalga,fftalgc,ig,ig1,ig1max,ig2,ig2max,ig3,ig3max,ii,idx,ierr
 integer :: is,m1,m2,m3,mm1,mm2,mm3,n1,n2,n3,nsym,nt,ount
 real(dp) :: ecuteff,ecutsigx,ecutwfn,g1,g2,g3,gsq,gsqmax,reff,rsigx,rwfn
 logical :: fft_ok
 character(len=500) :: msg, tnons_warn
!arrays
 integer :: fftnons(3),fftsym(3),mdum(3)
 !integer,allocatable :: pfactors(:),powers(:)
 integer,pointer :: symrel(:,:,:)
 real(dp),pointer :: tnons(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 if (ANY(mg0<0)) then
   write(msg,'(a,3(i0,1x))')' called with wrong value of mG0 = ',mG0
   MSG_BUG(msg)
 end if

 tnons_warn = "Check your fractional translations tnons. "//ch10//&
   "Components should be a rational fraction in 1/8th nor in 1/12th."//ch10//&
   "You may need to polish the structure by running AbiPy `abistruct.py abisanitize` on the input file."//ch10//&
   "to get rid of spurious tnons"

 ount = std_out; if (present(unit)) ount = unit

 nsym   =  Cryst%nsym
 symrel => Cryst%symrel
 tnons  => Cryst%tnons
 !
 ! Calculate the limits of the sphere of npwwfn G-vectors in each direction.
 m1=MAXVAL(ABS(gvec(1,1:npwwfn)))
 m2=MAXVAL(ABS(gvec(2,1:npwwfn)))
 m3=MAXVAL(ABS(gvec(3,1:npwwfn)))
 !
 ! Calculate the limits of the sphere of npsigx G-vectors in each direction.
 ! Ensure that G+G0 will fit into the FFT grid, where G is any of the npwsigx/npweps vectors
 ! and G0 is (i,j,k) [-nG0shell<i,j,k<nG0shell]. This is required when npwsigx>npwwfn since
 ! we have to take into account umklapp G0 vectors to evaluate the oscillator matrix elements
 ! (see rho_tw_g) or to symmetrize these quantities (see also cigfft).
 mm1=MAXVAL(ABS(gvec(1,1:npwsigx)))
 mm2=MAXVAL(ABS(gvec(2,1:npwsigx)))
 mm3=MAXVAL(ABS(gvec(3,1:npwsigx)))

 mm1=mm1+mG0(1)
 mm2=mm2+mG0(2)
 mm3=mm3+mG0(3)

 ! To avoid possible wrap-around errors in cigfft, it is safe to start
 ! with odd divisions so that the FFT box is centered on Gamma
 ! This holds only if npwsigx > npwwfn.
 !if (iseven(mm1)) mm1=mm1+1
 !if (iseven(mm2)) mm2=mm2+1
 !if (iseven(mm3)) mm3=mm3+1

 write(msg,'(2(2a,i8,a,3i6),2a,3i3)')ch10,&
  ' setmesh: npwwfn        = ',npwwfn, '; Max (m1,m2,m3)   = ',m1,m2,m3,ch10,&
  '          npweps/npwsigx= ',npwsigx,'; Max (mm1,mm2,mm3)= ',mm1,mm2,mm3,ch10,&
  '          mG0 added     = ',mG0(:)
 call wrtout(ount,msg,'COLL')
 !
 ! === Different FFT grids according to method ==
 select case (method)

 case (0)
   ! * FFT mesh defined by user, useful for testing.
   n1=ngfft(1)
   n2=ngfft(2)
   n3=ngfft(3)
   write(msg,'(3(a,i3))')' Mesh size enforced by user = ',n1,'x',n2,'x',n3
   MSG_COMMENT(msg)

   ngfft(1)=n1
   ngfft(2)=n2
   ngfft(3)=n3
   ngfft(4)=2*(ngfft(1)/2)+1
   ngfft(5)=2*(ngfft(2)/2)+1
   ngfft(6)=   ngfft(3)
   !ngfft(4:6)=ngfft(1:3)
   nfftot=n1*n2*n3
   RETURN

 case (1)
   aliasing_factor=1
   write(msg,'(2a,i3)')ch10,' using method 1 with aliasing_factor = ',aliasing_factor
   call wrtout(ount,msg,'COLL')
   m1=m1*aliasing_factor
   m2=m2*aliasing_factor
   m3=m3*aliasing_factor

 case (2,3)

   ecutwfn=-one  ! Calculate the radius of the sphere of npwwfn G-vectors.
   do ig=1,npwwfn
     g1=REAL(gvec(1,ig))
     g2=REAL(gvec(2,ig))
     g3=REAL(gvec(3,ig))
     gsq=       gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
          two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
     ecutwfn=MAX(ecutwfn,gsq)
   end do
   rwfn=SQRT(ecutwfn); ecutwfn=two*ecutwfn*pi**2

   ! * Calculate the radius of the sphere of (npwsigx|npweps) G-vectors.
   ecutsigx=-one
   do ig=1,npwsigx
     g1=REAL(gvec(1,ig))
     g2=REAL(gvec(2,ig))
     g3=REAL(gvec(3,ig))
     gsq=      gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
         two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
     ecutsigx=MAX(ecutsigx,gsq)
   end do
   rsigx=SQRT(ecutsigx); ecutsigx=two*ecutsigx*pi**2

   write(msg,'(a,f7.3,3a,f7.3,a)')&
    ' calculated ecutwfn          = ',ecutwfn, ' [Ha] ',ch10,&
    ' calculated ecutsigx/ecuteps = ',ecutsigx,' [Ha]'
   call wrtout(ount,msg,'COLL')
   !
   ! In the calculation of the GW self-energy or of the RPA dielectric matrix,
   ! we have products $ \rho_{12}(r)=u_1*(r) u_2(r) $ of wavefunctions whose Fourier
   ! coefficients lie in the sphere of radius rwfn. Such products will have non
   ! vanishing Fourier coefficients in the whole sphere of radius 2*rwfn since:
   !  $ rho_{12}(G) = \sum_T u_1*(T) u_2(T+G) $.
   ! However, we only need the Fourier coefficients of $rho_{12}$ that lie in the sphere
   ! of radius rsigx. We can thus allow aliasing outside that sphere, so that the FFT box
   ! will only enclose a sphere of radius reff given by:

   reff=rsigx+rwfn
   if (method==3) reff=two*rwfn ! Yields back the GS FFT grid if full wavefunctions are considered.
   ecuteff=two*(pi*reff)**2
   gsqmax=reff**2

   write(msg,'(a,i2,a,f7.3,a)')' using method = ',method,' with ecuteff = ',ecuteff,' [Ha]'
   call wrtout(ount,msg,'COLL')
   !
   ! === Search the limits of the reff sphere in each direction ===
   !ig1max=2*m1+1
   !ig2max=2*m2+1
   !ig3max=2*m3+1
   if (method==2) then
     ig1max=mm1+m1+1
     ig2max=mm2+m2+1
     ig3max=mm3+m3+1
   else if (method==3) then
     ig1max=MAX(2*m1+1,2*mm1+1,mm1+m1+1)
     ig2max=MAX(2*m2+1,2*mm2+1,mm2+m2+1)
     ig3max=MAX(2*m3+1,2*mm3+1,mm3+m3+1)
   else
     MSG_BUG(sjoin("Wrong method:", itoa(method)))
   end if

   m1=-1; m2=-1; m3=-1
   do ig1=0,ig1max
     do ig2=0,ig2max
       do ig3=0,ig3max
         g1=REAL(ig1)
         g2=REAL(ig2)
         g3=REAL(ig3)
         gsq=     gmet(1,1)*g1**2+gmet(2,2)*g2**2+gmet(3,3)*g3**2+ &
             two*(gmet(1,2)*g1*g2+gmet(1,3)*g1*g3+gmet(2,3)*g2*g3)
         if (gsq>gsqmax+tol6) CYCLE ! tol6 to improve portability
         m1=MAX(m1,ig1)
         m2=MAX(m2,ig2)
         m3=MAX(m3,ig3)
       end do
     end do
   end do

 case default
   MSG_BUG(sjoin('Method > 3 or < 0 not allowed in setmesh while method:', itoa(method)))
 end select
 !
 ! * Warning if low npwwfn.
 if (m1<mm1 .or. m2<mm2 .or. m3<mm3) then
   write(msg,'(5a)')&
    'Note that npwwfn is small with respect to npweps or with respect to npwsigx. ',ch10,&
    'Such a small npwwfn is a waste: ',ch10,&
    'You could raise npwwfn without loss in cpu time. '
   MSG_COMMENT(msg)
 end if
 !
 ! Keep the largest of the m/mm and and find the FFT grid which is compatible
 ! with the library and, if required, with the symmetry operations.
 m1=MAX(m1,mm1)
 m2=MAX(m2,mm2)
 m3=MAX(m3,mm3)

 if (enforce_sym==0) then
   ! === Determine the best size for the FFT grid *without* considering the symm ops ===
   ! * Ideally n=2*m+1 but this  could not be allowed by the FFT library.
   call size_goed_fft(m1,n1,ierr)
   call size_goed_fft(m2,n2,ierr)
   call size_goed_fft(m3,n3,ierr)
   ABI_CHECK(ierr == 0, sjoin("size_goed_fft failed", ch10, tnons_warn))
   nfftot=n1*n2*n3

   ! * Check if the FFT is compatible, write ONLY a warning if it breaks the symmetry
   fftnons(1)=n1
   fftnons(2)=n2
   fftnons(3)=n3
   fft_ok=.TRUE.
   rd: do ii=1,3
     do is=1,nsym
       nt=denominator(tnons(ii,is), ierr)
       if (((fftnons(ii)/nt)*nt) /= fftnons(ii)) then
         fft_ok=.FALSE.; EXIT rd
       end if
     end do
   end do rd
   !
   ! * Warn if not compatibile with tnons or rotational part.
   if (.not.fft_ok) then
    MSG_WARNING('FFT mesh is not compatible with non-symmorphic translations')
   end if
   if (.not.(check_rot_fft(nsym,symrel,n1,n2,n3))) then
     MSG_WARNING('FFT mesh is not compatible with rotations')
   end if

 else
   ! === Determine the best size for the FFT grid considering symm ops ===
   ! * Ideally n=2*m+1 but this could not be allowed by the FFT library (at present only Goedecker)
   call wrtout(ount,' Finding a FFT mesh compatible with all the symmetries','COLL')

   ! 1) Find a FFT mesh compatible with the non-symmorphic operations
   fftnons(:)=1
   do ii=1,3
     fftnons(ii)=1
     do is=1,nsym
       nt=denominator(tnons(ii,is), ierr)
       if (((fftnons(ii)/nt)*nt)/=fftnons(ii)) fftnons(ii)=mincm(fftnons(ii),nt)
     end do
   end do
   write(msg,'(a,3(i0,1x))')' setmesh: divisor mesh ',fftnons(:)
   call wrtout(ount,msg,'COLL')
   !
   ! 2) Check if also rotations preserve the grid.
   ! * Use previous m values as Initial guess.
   call size_goed_fft(m1,fftsym(1),ierr)
   call size_goed_fft(m2,fftsym(2),ierr)
   call size_goed_fft(m3,fftsym(3),ierr)
   ABI_CHECK(ierr == 0, sjoin("size_goed_fft failed", ch10, tnons_warn))
   mdum(1)=m1
   mdum(2)=m2
   mdum(3)=m3

   idx=0
   do ! If a FFT division gets too large the code stops in size_goed_fft.
     if ( check_rot_fft(nsym,symrel,fftsym(1),fftsym(2),fftsym(3)) .and. &
         (MOD(fftsym(1),fftnons(1))==0) .and.                           &
         (MOD(fftsym(2),fftnons(2))==0) .and.                           &
         (MOD(fftsym(3),fftnons(3))==0)                                 &
     ) EXIT
     ii=MOD(idx,3)+1
     mdum(ii)=mdum(ii)+1
     call size_goed_fft(mdum(ii),fftsym(ii),ierr)
     ABI_CHECK(ierr == 0, sjoin("size_goed_fft failed", ch10, tnons_warn))
     idx=idx+1
   end do
   !
   ! Got a good FFT grid, Calculate the number of FFT grid points
   n1=fftsym(1)
   n2=fftsym(2)
   n3=fftsym(3); nfftot=n1*n2*n3

   if (.not.( check_rot_fft(nsym,symrel,n1,n2,n3)) &
       .or.( MOD(fftsym(1),fftnons(1))/=0) .and.  &
           ( MOD(fftsym(2),fftnons(2))/=0) .and.  &
           ( MOD(fftsym(3),fftnons(3))/=0)        &
     ) then
     MSG_BUG('Not able to generate a symmetric FFT')
   end if
 end if ! enforce_sym

 write(msg,'(3(a,i5),2a,i12,a)')&
  ' setmesh: FFT mesh size selected  = ',n1,'x',n2,'x',n3,ch10,&
  '          total number of points  = ',nfftot,ch10
 call wrtout(ount,msg,'COLL')
 if (ount /= dev_null) call wrtout(ab_out,msg,'COLL')

 ngfft(1)=n1
 ngfft(2)=n2
 ngfft(3)=n3
 ngfft(4)=2*(ngfft(1)/2)+1
 ngfft(5)=2*(ngfft(2)/2)+1
 ngfft(6)=   ngfft(3)
 !ngfft(4:6) = ngfft(1:3)
 !
 ! === Check the value of fftalg i.e ngfft(7) ===
 ! * Presently only Goedecker"s library or FFTW3 are allowed, see size_goed_fft.F90
 fftalg=ngfft(7); fftalga=fftalg/100; fftalgc=MOD(fftalg,10)

 if ( ALL(fftalga /= (/FFT_SG,FFT_FFTW3, FFT_DFTI/)) ) then
   write(msg,'(6a)')ch10,&
    "Only Goedecker's routines with fftalg=1xx or FFTW3/DFTI routines are allowed in GW calculations. ",ch10,&
    "Action : check the value of fftalg in your input file, ",ch10,&
    "or modify setmesh.F90 to make sure the FFT mesh is compatible with the FFT library. "
   MSG_ERROR(msg)
 end if

! TODO Had to change setmesh to avoid bad values for FFTW3
! if (fftalga==3) then ! check whether mesh is optimal for FFTW3
!   ABI_MALLOC(pfactors,(5))
!   ABI_MALLOC(powers,(6))
!   pfactors = (/2, 3, 5, 7, 11/)
!   do ii=1,3
!     call pfactorize(ngfft(ii),5,pfactors,powers)
!     if (powers(6)/=1 .or. powers(4)/=0 .or. powers(5)/=0) then
!       write(msg,'(a,i0,a)')&
!&        "ngfft(ii) ",ngfft(ii)," contains powers of 7-11 or greater; FFTW3 is not optimal "
!       MSG_WARNING(msg)
!     end if
!   end do
!   ABI_FREE(pfactors)
!   ABI_FREE(powers)
! end if

 DBG_EXIT("COLL")

end subroutine setmesh
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/check_rot_fft
!! NAME
!!   check_rot_fft
!!
!! FUNCTION
!!  Return .TRUE. if the given grid in real space is compatible
!!  with the rotational part of the space group symmetries.
!!
!! INPUTS
!!  nsym=Number of symmetry operations
!!  symrel(3,3,nsym)=Symmetry operations in real space.
!!  nr1,nr2,nr3=FFT divisions.
!!
!! SOURCE

pure function check_rot_fft(nsym,symrel,nr1,nr2,nr3)

!Arguments
!Scalar
 integer,intent(in) :: nr1,nr2,nr3,nsym
 logical :: check_rot_fft
!Arrays
 integer,intent(in) :: symrel(3,3,nsym)

!local variables
 integer :: is

!************************************************************************

 ! The grid is compatible with the symmetries (only rotational part) if
 ! for each symmetry, each n_i and n_j ==> $n_i*R_{ij}/n_j$ is an integer
 check_rot_fft=.TRUE.
 do is=1,nsym
   if ( MOD(symrel(2,1,is)*nr2, nr1) /=0 .or. &
&       MOD(symrel(3,1,is)*nr3, nr1) /=0 .or. &
&       MOD(symrel(1,2,is)*nr1, nr2) /=0 .or. &
&       MOD(symrel(3,2,is)*nr3, nr2) /=0 .or. &
&       MOD(symrel(1,3,is)*nr1, nr3) /=0 .or. &
&       MOD(symrel(2,3,is)*nr2, nr3) /=0      &
&     ) then
     check_rot_fft=.FALSE.; EXIT
   end if
 end do

end function check_rot_fft
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/fft_check_rotrans
!! NAME
!! fft_check_rotrans
!!
!! FUNCTION
!!  Checks if the real space FFT mesh is compatible both with the rotational
!!  and the translational part of space group of the crystal.
!!
!! INPUTS
!!  nsym=Number of symmetries.
!!  symrel(3,3,nsym)=Symmetries in real space in reduced coordinates.
!!  tnons(3,nsym)=Fractional translations.
!!  ngfft(18)=Information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!
!! OUTPUT
!!  err(3,nsym)=The max error for each symmetry. (given in terms of the FFT vectors)
!!  isok=.FALSE. if the FFT mesh does not fulfil all symmetry properties of the crystal.
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!
!! SOURCE

function fft_check_rotrans(nsym,symrel,tnons,ngfft,err) result(isok)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 logical :: isok
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: tnons(3,nsym)
 real(dp),intent(out) :: err(3,nsym)

!Local variables-------------------------------
!scalars
 integer :: isym,ix,iy,iz,jx,jy,jz,ngfft1,ngfft2,ngfft3
 !character(len=500) :: msg
!arrays
 integer :: Rm1(3,3,nsym),r1_FFT(3),red2fft(3,3)
 real(dp) :: Rm1_FFT(3,3,nsym),fft2red(3,3),r2_FFT(3),tnons_FFT(3,nsym)

! *************************************************************************

 ! === Precalculate R^-1 and fractional translations in FFT coordinates ===
 ngfft1=ngfft(1)
 ngfft2=ngfft(2)
 ngfft3=ngfft(3)

 red2fft=RESHAPE((/ngfft1,0,0,0,ngfft2,0,0,0,ngfft3/),(/3,3/))
 fft2red=RESHAPE((/(one/ngfft1),zero,zero,zero,(one/ngfft2),zero,zero,zero,(one/ngfft3)/),(/3,3/))
 !
 ! === For a fully compatible mesh, each Rm1_FFT should be integer ===
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),Rm1(:,:,isym))
   Rm1(:,:,isym)=TRANSPOSE(Rm1(:,:,isym))
   Rm1_FFT(:,:,isym)=MATMUL(MATMUL(red2fft,Rm1(:,:,isym)),fft2red)
   tnons_FFT(:,isym)=MATMUL(red2fft,tnons(:,isym))
 end do

 err(:,:)=smallest_real
 do iz=0,ngfft3-1
   R1_FFT(3)=DBLE(iz)
   do iy=0,ngfft2-1
     R1_FFT(2)=DBLE(iy)
     do ix=0,ngfft1-1
       R1_FFT(1)=DBLE(ix)
       do isym=1,nsym  ! Form R^-1 (r-\tau) in the FFT basis ===
         R2_FFT(:)=MATMUL(Rm1_FFT(:,:,isym),R1_FFT(:)-tnons_FFT(:,isym))
         jx=NINT(R2_FFT(1)); err(1,isym)=MAX(err(1,isym),ABS(R2_FFT(1)-jx)/ngfft1)
         jy=NINT(R2_FFT(2)); err(2,isym)=MAX(err(2,isym),ABS(R2_FFT(2)-jy)/ngfft2)
         jz=NINT(R2_FFT(3)); err(3,isym)=MAX(err(3,isym),ABS(R2_FFT(3)-jz)/ngfft3)
       end do
     end do
   end do
 end do

 isok=.TRUE.
 do isym=1,nsym
   if (ANY(err(:,isym)>tol6)) then
     isok=.FALSE.
     !write(msg,'(a,i3,a,3es14.6)')' symmetry ',isym,') not compatible with FFT grid, error ',err(:,isym)
     !MSG_WARNING(msg)
   end if
 end do

end function fft_check_rotrans
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/rotate_FFT_mesh
!! NAME
!! rotate_FFT_mesh
!!
!! FUNCTION
!!  Find the FFT index of $ R{-1}(r-\tau) $ for each point in the FFT box.
!!  $R$ is a symmetry operation in real space, $\tau$ is the associated
!!  fractional translation.
!!
!! INPUTS
!!  nsym=Number of symmetries.
!!  symrel(3,3,nsym)=Symmetries in real space in reduced coordinates.
!!  tnons(3,nsym)=Fractional translations.
!!  ngfft(18)=Information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!
!! OUTPUT
!!  irottb(ngfftot,nsym)=Indeces of $R^{-1}(r-\tau)$ in the FFT box.
!!  preserve=.FALSE. if the FFT mesh does not fulfil all symmetry properties of the crystal.
!!
!! NOTES
!!  The evaluation of the rotated point $R^{-1}(r-\tau)$ is done using real arithmetic.
!!  As a consequence, if the FFT mesh does not fulfil the symmetry properties
!!  of the crystal, the array irottb will contain the index of the FFT point which
!!  is the closest one to $R^{-1}(r-\tau)$. This might lead to inaccuracies in the
!!  final results, in particular in the description of degenerate states.
!!
!! PARENTS
!!      m_bethe_salpeter,m_chi0,m_classify_bands,m_cohsex,m_dvdb,m_fft_mesh
!!      m_prep_calc_ucrpa,m_screening_driver,m_sigc,m_sigx,m_wfd
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine rotate_fft_mesh(nsym,symrel,tnons,ngfft,irottb,preserve)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 logical,intent(out) :: preserve
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 integer,intent(in) :: ngfft(18)
 integer,intent(out) :: irottb(ngfft(1)*ngfft(2)*ngfft(3),nsym)
 real(dp),intent(in) :: tnons(3,nsym)

!Local variables-------------------------------
!scalars
 integer :: ir1,isym,ix,iy,iz,jx,jy,jz,ngfft1,ngfft2,ngfft3
 !character(len=500) :: msg
!arrays
 integer :: Rm1(3,3,nsym),r1_FFT(3),red2fft(3,3)
 real(dp) :: Rm1_FFT(3,3,nsym),err(3,nsym),fft2red(3,3),r2_FFT(3)
 real(dp) :: tnons_FFT(3,nsym)

! *************************************************************************

 ! === Precalculate R^-1 and fractional translations in FFT coordinates ===
 ngfft1=ngfft(1)
 ngfft2=ngfft(2)
 ngfft3=ngfft(3)

 red2fft=RESHAPE((/ngfft1,0,0,0,ngfft2,0,0,0,ngfft3/),(/3,3/))
 fft2red=RESHAPE((/(one/ngfft1),zero,zero,zero,(one/ngfft2),zero,zero,zero,(one/ngfft3)/),(/3,3/))
 !
 ! === For a fully compatible mesh, each Rm1_FFT should be integer ===
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),Rm1(:,:,isym))
   Rm1(:,:,isym)=TRANSPOSE(Rm1(:,:,isym))
   Rm1_FFT(:,:,isym)=MATMUL(MATMUL(red2fft,Rm1(:,:,isym)),fft2red)
   tnons_FFT(:,isym)=MATMUL(red2fft,tnons(:,isym))
 end do

 err(:,:)=zero

!$OMP PARALLEL DO PRIVATE(R1_FFT,ir1,R2_FFT,jx,jy,jz) reduction(MAX:err)
 do iz=0,ngfft3-1
   R1_FFT(3)=DBLE(iz)
   do iy=0,ngfft2-1
     R1_FFT(2)=DBLE(iy)
     do ix=0,ngfft1-1
       R1_FFT(1)=DBLE(ix)
       ir1=1+ix+iy*ngfft1+iz*ngfft1*ngfft2
       do isym=1,nsym
         ! === Form R^-1 (r-\tau) in the FFT basis ===
         R2_FFT(:)=MATMUL(Rm1_FFT(:,:,isym),R1_FFT(:)-tnons_FFT(:,isym))
         jx=NINT(R2_FFT(1)); err(1,isym)=MAX(err(1,isym),ABS(R2_FFT(1)-jx)/ngfft1)
         jy=NINT(R2_FFT(2)); err(2,isym)=MAX(err(2,isym),ABS(R2_FFT(2)-jy)/ngfft2)
         jz=NINT(R2_FFT(3)); err(3,isym)=MAX(err(3,isym),ABS(R2_FFT(3)-jz)/ngfft3)
         jx=MODULO(jx,ngfft1)
         jy=MODULO(jy,ngfft2)
         jz=MODULO(jz,ngfft3)
         irottb(ir1,isym)=1+jx+jy*ngfft1+jz*ngfft1*ngfft2
       end do
     end do
   end do
 end do

 preserve=.TRUE.
 do isym=1,nsym
   if (ANY(err(:,isym)>tol6)) then
     preserve=.FALSE.
     !write(msg,'(a,i3,a,3es14.6)')' symmetry ',isym,') not compatible with FFT grid, error ',err(:,isym)
     !MSG_WARNING(msg)
   end if
 end do

end subroutine rotate_fft_mesh
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/denpot_project
!! NAME
!!
!! FUNCTION
!!  Compute n(r) + n( $R^{-1}(r-\tau)$ in) / 2
!!  Mainly used with R = inversion to select the even/odd part under inversion
!!
!! INPUTS
!!  cplex=1 for real, 2 for complex data.
!!  ngfft(3)=Mesh divisions of input array
!!  nspden=Number of density components.
!!  in_rhor(cplex * nfftot * nspden)=Input array
!!  one_symrel(3,3)= R operation
!!  tau(3)=Fractional translation.
!!
!! OUTPUT
!!  out_rhor(cplex * nfftot * nspden)=Output array
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine denpot_project(cplex,  ngfft, nspden, in_rhor, one_symrel, one_tnons, out_rhor)

!Arguments-------------------------------------------------------------
!scalars
 integer,intent(in) :: cplex, nspden
!arrays
 integer,intent(in) :: ngfft(18), one_symrel(3,3)
 real(dp),intent(in) :: in_rhor(cplex, product(ngfft(1:3)), nspden)
 real(dp),intent(in) :: one_tnons(3)
 real(dp),intent(out) :: out_rhor(cplex, product(ngfft(1:3)), nspden)

!Local variables--------------------------------------------------------
!scalars
 integer,parameter :: nsym1 = 1, isgn = 1
 integer :: ispden, ii, ifft, ifft_rot, nfft
 logical :: preserve
!arrays
 integer,allocatable :: irottb(:)

! *************************************************************************

 nfft = product(ngfft(1:3))
 ABI_MALLOC(irottb, (nfft))

 call rotate_fft_mesh(nsym1, one_symrel, one_tnons, ngfft, irottb, preserve)
 ABI_CHECK(preserve, "FFT mesh is not compatible with {R, tau}")

 do ispden=1,nspden
   do ifft=1,nfft
     ifft_rot = irottb(ifft)
     do ii=1,cplex
       out_rhor(cplex, ifft, ispden) = (in_rhor(cplex, ifft, ispden) + isgn * in_rhor(cplex, ifft_rot, ispden)) * half
     end do
   end do
 end do

 ABI_FREE(irottb)

end subroutine denpot_project
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/cigfft
!! NAME
!! cigfft
!!
!! FUNCTION
!! For each of the (2*nG0sh+1)**3 vectors G0 around the origin,
!! calculate G-G0 and its FFT index number for all the NPWVEC vectors G.
!!
!! INPUTS
!! mG0(3)= For each reduced direction gives the max G0 component to account for umklapp processes.
!! npwvec=Number of plane waves
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! gvec(3,npwvec)=Reduced coordinates of G vectors.
!!
!! OUTPUT
!! igfft(npwvec,2*mG0(1)+1,2*mG0(2)+1,2*mG0(3)+1)=For each G, and each G0 vector,
!!  it gives the FFT grid index of the G-G0 vector.
!! ierr=Number of G-G0 vectors falling outside the inout FFT box.
!!
!! PARENTS
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine cigfft(mG0,npwvec,ngfft,gvec,igfft,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwvec
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: gvec(3,npwvec)
 integer,intent(in) :: mg0(3),ngfft(18)
 integer,intent(out) :: igfft(npwvec,2*mg0(1)+1,2*mg0(2)+1,2*mg0(3)+1)

!Local variables ------------------------------
!scalars
 integer :: gmg01,gmg02,gmg03,ig,ig01,ig02,ig03,n1,n2,n3
 character(len=500) :: msg
!arrays
 integer :: gmg0(3)
!************************************************************************

 DBG_ENTER("COLL")

 if (ANY(mg0<0)) then
   write(msg,'(a,3i4)')' Found negative value of mg0= ',mg0
   MSG_BUG(msg)
 end if

 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 ierr=0

 do ig=1,npwvec
   do ig01=-mg0(1),mG0(1)
     gmg0(1) = gvec(1,ig)-ig01
     do ig02=-mg0(2),mg0(2)
       gmg0(2) = gvec(2,ig)-ig02
       do ig03=-mg0(3),mg0(3)
         gmg0(3) = gvec(3,ig)-ig03
         ! === Calculate FFT index of G-G0 ===
         ! * Consider possible wrap around errors.
         gmg01=MODULO(gmg0(1),n1)
         gmg02=MODULO(gmg0(2),n2)
         gmg03=MODULO(gmg0(3),n3)
         igfft(ig,ig01+mg0(1)+1,ig02+mg0(2)+1,ig03+mg0(3)+1) = 1+gmg01+gmg02*n1+gmg03*n1*n2
         if ( ANY(gmg0>ngfft(1:3)/2) .or. ANY(gmg0<-(ngfft(1:3)-1)/2) ) then
           igfft(ig,ig01+mg0(1)+1,ig02+mg0(2)+1,ig03+mg0(3)+1) = 0
           ierr=ierr+1
         end if
       end do
     end do
   end do
 end do !ig

 if (ierr/=0) then
   write(msg,'(a,i0,3a)')&
    'Found ',ierr,' G-G0 vectors falling outside the FFT box. ',ch10,&
    'igfft will be set to zero for these particular G-G0 '
   MSG_WARNING(msg)
 end if

 DBG_EXIT("COLL")

end subroutine cigfft
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/ig2gfft
!! NAME
!!  ig2gfft
!!
!! FUNCTION
!!  Return the reduced component of a G-vector in the FFT mesh starting from is index.
!!
!! INPUTS
!!  ig = The index >=1, <=ng
!!  ng = The number of FFT points along this direction.
!!
!! OUTPUT
!!  gc = The reduced component
!!
!! SOURCE

elemental integer function ig2gfft(ig,ng) result (gc)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ig,ng

!************************************************************************

 ! Use the following indexing (N means ngfft of the adequate direction)
 ! 0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= gc
 ! 1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index ig
 !
 if (ig<=0 .or. ig > ng) then
   ! Wrong ig, returns huge. Parent code will likely crash with SIGSEV.
   gc = huge(1)
   return
 end if

 if ( ig  > ng/2 + 1) then
   gc = ig - ng -1
 else
   gc = ig -1
 end if

end function ig2gfft
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/g2ifft
!! NAME
!!  g2ifft
!!
!! FUNCTION
!! Returns the index of G in the FFT box from its reduced coordinates. 0 if not in the BOX.
!!
!! INPUTS
!!  gg(3)=Reduced coordinated of the G vector.
!!  ngfft(18) = Info on the FFT box.
!!
!! OUTPUT
!!  gidx=Index in the FFT box. 0 if G is outside the box.
!!
!! SOURCE

pure integer function g2ifft(gg,ngfft) result (gidx)

!Arguments ------------------------------------
 integer,intent(in) :: gg(3),ngfft(3)

!Local variables-------------------------------
!scalars
 integer :: n1,n2,n3,ig1,ig2,ig3

!************************************************************************

 ! Use the following indexing (N means ngfft of the adequate direction)
 ! 0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= gg
 ! 1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index
 !
 if ( ANY(gg>ngfft(1:3)/2) .or. ANY(gg<-(ngfft(1:3)-1)/2) ) then ! out of the box.
   gidx=0
 else
   n1=ngfft(1)
   n2=ngfft(2)
   n3=ngfft(3)

   ig1=MODULO(gg(1),n1)
   ig2=MODULO(gg(2),n2)
   ig3=MODULO(gg(3),n3)

   gidx = 1 + ig1 + n1*(ig2+ig3*n2)
 end if

end function g2ifft
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/get_gftt
!! NAME
!!  get_gftt
!!
!! FUNCTION
!!  Returns the set of G-vectors in the FFT mesh and the maximal kinetic energy of k+G.
!!
!! INPUTS
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! kpt(3)=input k vector (reduced coordinates --in terms of reciprocal lattice primitive translations)
!! gmet(3,3)=reciprocal space metric (bohr^-2)
!!
!! OUTPUT
!!  gsq_max=Max value of (k+G)^2 for G in the FFT box
!!  gfft(3,nfft_tot) = The reduced components of the G in the FFT mesh (nfft_tot=PRODUCT(ngfft(1:3))
!!
!! PARENTS
!!      bethe_salpeter,calc_sigc_me,cchi0,cchi0q0,cohsex_me,screening,sigma
!!      m_dvdb
!!
!! CHILDREN
!!      xcopy
!!
!! SOURCE

pure subroutine get_gftt(ngfft, kpt, gmet, gsq_max, gfft)

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: gsq_max
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(out) :: gfft(3,ngfft(1)*ngfft(2)*ngfft(3))
 real(dp),intent(in) :: kpt(3),gmet(3,3)

!Local variables-------------------------------
!scalars
 integer :: ifft,g1,g2,g3,i1,i2,i3
 real(dp) :: dsq

!************************************************************************

 ifft=0; gsq_max=smallest_real
 do i3=1,ngfft(3)
   g3 = ig2gfft(i3,ngfft(3))
   do i2=1,ngfft(2)
     g2 = ig2gfft(i2,ngfft(2))
     do i1=1,ngfft(1)
       g1 = ig2gfft(i1,ngfft(1))
       ifft = ifft+1
       gfft(1,ifft) = g1
       gfft(2,ifft) = g2
       gfft(3,ifft) = g3
       dsq=gmet(1,1)*(kpt(1)+dble(i1))**2 &
        +gmet(2,2)*(kpt(2)+dble(i2))**2 &
        +gmet(3,3)*(kpt(3)+dble(i3))**2 &
        +2._dp*(gmet(1,2)*(kpt(1)+dble(i1))*(kpt(2)+dble(i2)) &
        +gmet(2,3)*(kpt(2)+dble(i2))*(kpt(3)+dble(i3)) &
        +gmet(3,1)*(kpt(3)+dble(i3))*(kpt(1)+dble(i1)))
       gsq_max = MAX(dsq,gsq_max)
     end do
   end do
 end do

end subroutine get_gftt
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/calc_ceigr_spc
!! NAME
!! calc_ceigr_spc
!!
!! FUNCTION
!!  Helper function to calculate e^{iG.r} on the FFT mesh.
!!
!! INPUTS
!!  gg(3)=G vector in reduced coordinates.
!!  nfft=Total number of points in the FFT mesh.
!!  nspinor=Number of spinors
!!  ngfft(18)=information about 3D FFT,
!!
!! OUTPUT
!!  ceigr(nfft*nspinor)=e^{ik.r} on the FFT mesh.
!!
!! PARENTS
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine calc_ceigr_spc(gg,nfft,nspinor,ngfft,ceigr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nspinor
!arrays
 integer,intent(in) :: gg(3)
 integer,intent(in) :: ngfft(18)
 complex(spc),intent(out) :: ceigr(nfft*nspinor)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,fft_idx,base,isp
 real(dp) :: gdotr

! *************************************************************************

 if (ALL(gg==0)) then
   ceigr=(1._sp,0._sp)
   RETURN
 end if

 fft_idx=0
 do iz=0,ngfft(3)-1
   do iy=0,ngfft(2)-1
     do ix=0,ngfft(1)-1
       gdotr= two_pi*( gg(1)*(ix/DBLE(ngfft(1))) &
&                     +gg(2)*(iy/DBLE(ngfft(2))) &
&                     +gg(3)*(iz/DBLE(ngfft(3))) )
       fft_idx = fft_idx+1
       ceigr(fft_idx)=CMPLX(DCOS(gdotr),DSIN(gdotr), KIND=spc)
     end do
   end do
 end do

 if (nspinor > 1) then
   do isp=2,nspinor
     base = 1 + (isp-1)*nfft
     call xcopy(nfft,ceigr,1,ceigr(base:),1)
   end do
 end if

end subroutine calc_ceigr_spc
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/calc_ceigr_dpc
!! NAME
!! calc_ceigr_dpc
!!
!! FUNCTION
!!  Helper function to calculate e^{iG.r} on the FFT mesh.
!!
!! INPUTS
!!  gg(3)=G vector in reduced coordinates.
!!  nfft=Total number of points in the FFT mesh.
!!  nspinor=Number of spinors
!!  ngfft(18)=information about 3D FFT,
!!
!! OUTPUT
!!  ceigr(nfft*nspinor)=e^{ik.r} on the FFT mesh.
!!
!! PARENTS
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine calc_ceigr_dpc(gg,nfft,nspinor,ngfft,ceigr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nspinor
!arrays
 integer,intent(in) :: gg(3)
 integer,intent(in) :: ngfft(18)
 complex(dpc),intent(out) :: ceigr(nfft*nspinor)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,fft_idx,base,isp
 real(dp) :: gdotr

! *************************************************************************

 if (ALL(gg==0)) then
   ceigr=cone; RETURN
 end if

 fft_idx=0
 do iz=0,ngfft(3)-1
   do iy=0,ngfft(2)-1
     do ix=0,ngfft(1)-1
       gdotr= two_pi*( gg(1)*(ix/DBLE(ngfft(1))) &
&                     +gg(2)*(iy/DBLE(ngfft(2))) &
&                     +gg(3)*(iz/DBLE(ngfft(3))) )
       fft_idx = fft_idx+1
       ceigr(fft_idx)=DCMPLX(DCOS(gdotr),DSIN(gdotr))
     end do
   end do
 end do

 if (nspinor > 1) then
   do isp=2,nspinor
     base = 1 + (isp-1)*nfft
     call xcopy(nfft,ceigr,1,ceigr(base:),1)
   end do
 end if

end subroutine calc_ceigr_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/calc_eigr
!! NAME
!! calc_eigr
!!
!! FUNCTION
!!  Helper function to calculate e^{iG.r} on the FFT mesh.
!!
!! INPUTS
!!  gg(3)=G vector in reduced coordinates.
!!  nfft=Total number of points in the FFT mesh.
!!  ngfft(18)=information about 3D FFT,
!!
!! OUTPUT
!!  eigr(2*nfft)=e^{ig.r} on the FFT mesh.
!!
!! PARENTS
!!      m_fft_prof
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine calc_eigr(gg,nfft,ngfft,eigr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
!arrays
 integer,intent(in) :: gg(3)
 integer,intent(in) :: ngfft(18)
 real(dp),intent(out) :: eigr(2*nfft)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,fft_idx
 real(dp) :: gdotr

! *************************************************************************

 if (ALL(gg==0)) then
   eigr(1:2*nfft:2)=one
   eigr(2:2*nfft:2)=zero
   RETURN
 end if

 fft_idx=1
 do iz=0,ngfft(3)-1
   do iy=0,ngfft(2)-1
     do ix=0,ngfft(1)-1
       gdotr= two_pi*( gg(1)*(ix/DBLE(ngfft(1))) &
&                     +gg(2)*(iy/DBLE(ngfft(2))) &
&                     +gg(3)*(iz/DBLE(ngfft(3))) )
       eigr(fft_idx  )=DCOS(gdotr)
       eigr(fft_idx+1)=DSIN(gdotr)
       fft_idx = fft_idx+2
     end do
   end do
 end do

end subroutine calc_eigr
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/calc_ceikr
!! NAME
!! calc_ceikr
!!
!! FUNCTION
!!  Helper function to calculate e^{ik.r} on the FFT mesh.
!!
!! INPUTS
!!  kk(3)=k-point in reduced coordinates.
!!  nfft=Total number of points in the FFT mesh.
!!  ngfft(18)=information about 3D FFT,
!!
!! OUTPUT
!!  ceikr(nfft)=e^{ik.r} on the FFT mesh.
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine calc_ceikr(kk,nfft,ngfft,ceikr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
!arrays
 real(dp),intent(in) :: kk(3)
 integer,intent(in) :: ngfft(18)
 complex(dpc),intent(out) :: ceikr(nfft)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,fft_idx
 real(dp) :: kdotr

! *************************************************************************

 !if (ALL(ABS(kk<tol12)) then
 !  ceikr=cone; RETURN
 !end if

 fft_idx=0
 do iz=0,ngfft(3)-1
   do iy=0,ngfft(2)-1
     do ix=0,ngfft(1)-1
       kdotr= two_pi*( kk(1)*(ix/DBLE(ngfft(1))) &
&                     +kk(2)*(iy/DBLE(ngfft(2))) &
&                     +kk(3)*(iz/DBLE(ngfft(3))) )
       fft_idx = fft_idx+1
       ceikr(fft_idx)=DCMPLX(DCOS(kdotr),DSIN(kdotr))
     end do
   end do
 end do

end subroutine calc_ceikr
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/times_eigr
!! NAME
!! times_eigr
!!
!! FUNCTION
!!  Multiply an array on the real-space mesh by e^{iG0.r} where G0 is a reciprocal lattice vector.
!!
!! INPUTS
!!  gg(3)=G vector in reduced coordinates.
!!  ngfft(18)=information about 3D FFT,
!!  nfft=Number of points in the FFT mesh.
!!  ndat=Number of arrays
!!
!! SIDE EFFECTS
!!  ur(2,nfft,ndat)= contains u(r) in input. output: u(r) e^{ig.r} on the real-space FFT mesh.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine times_eigr(gg,ngfft,nfft,ndat,ur)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ndat
!arrays
 integer,intent(in) :: gg(3)
 integer,intent(in) :: ngfft(18)
 real(dp),intent(inout) :: ur(2,nfft,ndat)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,ifft,idat
 real(dp) :: gr
!arrays
 real(dp) :: ph(2),val(2)

! *************************************************************************

 if (all(gg==0)) return

 do idat=1,ndat
   ifft = 0
   do iz=0,ngfft(3)-1
     do iy=0,ngfft(2)-1
       do ix=0,ngfft(1)-1
         ifft = ifft + 1
         gr = two_pi*(gg(1)*(ix/dble(ngfft(1))) &
                     +gg(2)*(iy/dble(ngfft(2))) &
                     +gg(3)*(iz/dble(ngfft(3))) )
         ph(1) = cos(gr); ph(2) = sin(gr)
         val(1) = ur(1,ifft,idat); val(2) = ur(2,ifft,idat)

         ur(1,ifft,idat) = ph(1) * val(1) - ph(2) * val(2)
         ur(2,ifft,idat) = ph(1) * val(2) + ph(2) * val(1)
       end do
     end do
   end do
 end do ! idat

end subroutine times_eigr
!!***

!----------------------------------------------------------------------

!!****f* m_fft_mesh/times_eikr
!! NAME
!! times_eikr
!!
!! FUNCTION
!!  Multiply an array on the real-space mesh by e^{ik.r} where k is a real(dp) vector in reduced coordinates
!!
!! INPUTS
!!  kk(3)=k-vector in reduced coordinates.
!!  ngfft(18)=information about 3D FFT,
!!  nfft=Number of points in the FFT mesh.
!!  ndat=Number of arrays to transform
!!
!! SIDE EFFECTS
!!  ur(2,nfft)= contains u(r) in input. output: u(r) e^{ig.r} on the real-space FFT mesh.
!!
!! PARENTS
!!  m_dvdb
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine times_eikr(kk,ngfft,nfft,ndat,ur)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ndat
!arrays
 real(dp),intent(in) :: kk(3)
 integer,intent(in) :: ngfft(18)
 real(dp),intent(inout) :: ur(2,nfft,ndat)

!Local variables-------------------------------
!scalars
 integer :: ix,iy,iz,ifft,idat
 real(dp) :: kr
!arrays
 real(dp) :: ph(2),val(2)

! *************************************************************************

 if (all(abs(kk) < tol12)) return

 do idat=1,ndat
   ifft = 0
   do iz=0,ngfft(3)-1
     do iy=0,ngfft(2)-1
       do ix=0,ngfft(1)-1
         ifft = ifft + 1
         kr = two_pi*(kk(1)*(ix/dble(ngfft(1))) &
                     +kk(2)*(iy/dble(ngfft(2))) &
                     +kk(3)*(iz/dble(ngfft(3))) )
         ph(1) = cos(kr); ph(2) = sin(kr)
         val(1) = ur(1,ifft,idat); val(2) = ur(2,ifft,idat)

         ur(1,ifft,idat) = ph(1) * val(1) - ph(2) * val(2)
         ur(2,ifft,idat) = ph(1) * val(2) + ph(2) * val(1)
       end do
     end do
   end do
 end do

end subroutine times_eikr
!!***

!!****f* m_fft_mesh/phase
!! NAME
!! phase
!!
!! FUNCTION
!! Compute ph(ig)=$\exp(\pi\ i \ n/ngfft)$ for n=0,...,ngfft/2,-ngfft/2+1,...,-1
!! while ig runs from 1 to ngfft.
!!
!! INPUTS
!!  ngfft=number of points
!!
!! OUTPUT
!!  ph(2*ngfft)=phase array (complex)
!!
!! NOTES
!! XG 990504 : changed the formulation, in order to preserve
!! the invariance between n and -n, that was broken for n=ngfft/2 if ngfft even.
!! Simply suppresses the corresponding sine.
!!
!! PARENTS
!!      m_xctk
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine phase(ngfft,ph)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ngfft
!arrays
 real(dp),intent(out) :: ph(2*ngfft)

!Local variables-------------------------------
!scalars
 integer :: id,ig,nn
 real(dp) :: arg,fac

! *************************************************************************

 id=ngfft/2+2
 fac=pi/dble(ngfft)
 do ig=1,ngfft
   nn=ig-1-(ig/id)*ngfft
   arg=fac*dble(nn)
   ph(2*ig-1)=cos(arg)
   ph(2*ig)  =sin(arg)

 end do
!XG 990504 Here zero the corresponding sine
 if((ngfft/2)*2==ngfft) ph(2*(id-1))=zero

end subroutine phase
!!***

!!****f* ABINIT/mkgrid_fft
!! NAME
!!  mkgrid_fft
!!
!! FUNCTION
!!  Sets the grid of fft (or real space) points to be treated.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_mklocl_realspace
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

subroutine mkgrid_fft(ffti3_local,fftn3_distrib,gridcart,nfft,ngfft,rprimd)

!Arguments ------------------------------------
 integer, intent(in) :: nfft
 integer,intent(in) :: ngfft(18)
 integer, dimension(*), intent(in) :: ffti3_local,fftn3_distrib
 real(dp), dimension(3,nfft), intent(out) :: gridcart
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
 integer :: ind,i1,i2,i3,i3loc,me,nproc
 integer :: n1,n2,n3
 real(dp), dimension(3) :: coord
 real(dp), dimension(3,nfft) :: gridred

! *************************************************************************

 n1    = ngfft(1)
 n2    = ngfft(2)
 n3    = ngfft(3)
 nproc = ngfft(10)
 me    = ngfft(11)

 do i3 = 1, n3, 1
   if(fftn3_distrib(i3) == me) then !MPI
     i3loc=ffti3_local(i3)
     coord(3) = real(i3 - 1, dp) / real(n3, dp)
     do i2 = 1, n2, 1
       coord(2) = real(i2 - 1, dp) / real(n2, dp)
       do i1 = 1, n1, 1
         ind=i1+(i2-1)*n1+(i3loc-1)*n1*n2
         coord(1) = real(i1 - 1, dp) / real(n1, dp)
         gridred(:, ind) = coord(:)
       end do
     end do
   end if
 end do
 call xred2xcart(nfft, rprimd, gridcart, gridred)

end subroutine mkgrid_fft
!!***

END MODULE m_fft_mesh
!!***
