!!****m* ABINIT/m_barevcoul
!! NAME
!!  m_barevcoul
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1999-2025 ABINIT group ()
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_barevcoul

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi

 use m_fstrings,        only : sjoin
 use defs_abitypes,     only : MPI_type
 use m_numeric_tools,   only : arth, l2norm, OPERATOR(.x.)
 use m_geometry,        only : normv
 use m_crystal,         only : crystal_t
 use m_fft,             only : zerosym
 !use m_gsphere,         only : gsphere_t

! Cut-off methods modules
 !use m_cutoff_sphere,   only : cutoff_sphere
 !use m_cutoff_slab,     only : cutoff_slab
 !use m_cutoff_cylinder, only : cutoff_cylinder

 implicit none

 private
!!***

!!****t* m_barevcoul/vcut_t
!! NAME
!!  vcoul_t
!!
!! FUNCTION
!!
!! SOURCE

 type,public :: vcut_t

  integer  :: nfft
  ! Number of points in FFT grid

  integer  :: ng
   ! Number of G-vectors

  real(dp) :: alpha(3)
   ! Lenght of the finite slab

  real(dp) :: rcut
   ! Cutoff radius

  real(dp) :: i_sz
   ! Value of the integration of the Coulomb singularity 4\pi/V_BZ \int_BZ d^3q 1/q^2

  real(dp) :: hcyl
   ! Length of the finite cylinder along the periodic dimension

  real(dp) :: ucvol
    ! Volume of the unit cell

  character(len=50) :: mode
   ! String defining the cutoff mode, possible values are: sphere,cylinder,slab,crystal

  integer :: pdir(3)
   ! 1 if the system is periodic along this direction

  real(dp) :: boxcenter(3)
   ! 1 if the point in inside the cutoff region 0 otherwise
   ! Reduced coordinates of the center of the box (input variable)

  real(dp) :: vcutgeo(3)
    ! For each reduced direction gives the length of the finite system
    ! 0 if the system is infinite along that particular direction
    ! negative value to indicate that a finite size has to be used

  real(dp) :: rprimd(3,3)
    ! Lattice vectors in real space.

  real(dp),allocatable :: qibz(:,:)
   ! qibz(3,nqibz)
   ! q-points in the IBZ.

  real(dp),allocatable :: barev(:)
    ! barev(nfft)
    ! Bare Coulomb potential on the FFT grid
    ! A cut might be applied.

 end type vcut_t

 public :: barevcoul
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/barevcoul
!! NAME
!! barevcoul
!!
!! FUNCTION
!! Compute bare coulomb term in G-space on the FFT mesh i.e. 4pi/(G+q)**2 for a specified q-point
!!
!! INPUTS
!!  qpoint(3)=reduced coordinates for the phonon wavelength
!!  gsqcut=cutoff value on G**2 for sphere inside fft box. (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2))
!!  icutcoul=Option for the Coulomb potential cutoff technique
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  izero=if 1, unbalanced components of V(q,g) are set to zero # Used by the PAW library
!!  nfft=Total number of FFT grid points.
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  barev(nfft)=4pi/(G+q)**2, q+G=0 component is set carefully
!!
!! NOTES
!!  This routine operates on the full FFT mesh. DO NOT PASS MPI_TYPE
!!  One can easily implemente MPI-FFT by just calling this routine and then
!!  extracting the G-vectors treated by the node.
!!
!! SOURCE

subroutine barevcoul(rcut,icutcoul,qpoint,gsqcut,gmet,nfft,nkpt_bz,ngfft,ucvol,izero,barev,shortrange)

!Arguments ------------------------------------
!scalars
 integer,intent(in)         :: icutcoul,nfft,nkpt_bz,izero
 real(dp),intent(in)        :: rcut,gsqcut,ucvol
 logical,intent(in),optional:: shortrange
!arrays
 integer,intent(in)         :: ngfft(18)
 integer                    :: ng!!!!
 real(dp),intent(in)        :: qpoint(3)
 real(dp),intent(inout)     :: gmet(3,3)
 real(dp),intent(inout)     :: barev(nfft)
 !real(dp)                   :: a1(3),a2(3),a3(3)
 real(dp)                   :: b1(3),b2(3),b3(3),rmet(3,3) !,gprimd(3,3),
 type(MPI_type)             :: mpi_enreg   !!!!
 type(crystal_t)            :: Cryst       !!!!
 !type(gsphere_t)            :: Gsph
!Local variables-------------------------------
!scalars
 integer,parameter    :: empty(3,3)=zero
 integer,parameter    :: cplex1=1
 integer              :: comm
 integer              :: ii1,i1,i2,i23,i3,id1,id2,id3,icutcoul_local
 integer              :: ig,ig1,ig2,ig3,ig1min,ig1max,ig2min,ig2max,ig3min,ig3max
 integer              :: ii,ing,n1,n2,n3,npar,npt
 integer              :: opt_cylinder,opt_slab,test
 integer              :: qeq0,qeq05
 real(dp),parameter   :: tolfix=1.000000001e0_dp ! Same value as the one used in hartre
 real(dp)             :: check,step
 real(dp)             :: cutoff,gqg2p3,gqgm12,gqgm13,gqgm23,gs,gs2,gs3,rcut0
 real(dp)             :: bz_plane,dx,integ,q0_vol,q0_volsph
 character(len=500)   :: msg
!arrays
 integer              :: id(3), gamma_pt(3,1)
 real(dp),allocatable :: gq(:,:),gpq(:),gpq2(:)
 real(dp),allocatable :: vcfit(:,:),xx(:),yy(:)
 real(dp),allocatable :: cov(:,:),par(:),qfit(:,:),sigma(:),var(:),qcart(:,:)
 type(vcut_t)               :: vcut        !!!!
!
 comm=mpi_enreg%comm_world
!
! === Save dimension and other useful quantities in vcut% ===
 vcut%nfft      = PRODUCT(ngfft(1:3))  ! Number of points in the FFT mesh.
 ! ng and gvec are not used yet we don't want to use them without being defined
 ng = -1
 vcut%ng        = ng                   ! Number of G-vectors in the Coulomb matrix elements.
 vcut%rcut      = rcut                 ! Cutoff radius for cylinder.
 vcut%hcyl      = zero                 ! Length of finite cylinder (Rozzi"s method, default is Beigi).
 vcut%ucvol     = ucvol                ! Unit cell volume.

 !FBruneval: comment the definitions below since Cryst and dtset have never been initialized!
 !vcut%rprimd    = Cryst%rprimd(:,:)    ! Dimensional direct lattice.
 !vcut%boxcenter = dtset%boxcenter      ! boxcenter at the moment is supposed to be at the origin.
 !vcut%vcutgeo   = dtset%vcutgeo(:)     ! Info on the orientation and extension of the cutoff region.
!
! === Define geometry and cutoff radius (if used) ===
 vcut%mode='NONE'
 icutcoul_local=icutcoul

 ! for short-range exchange (e.g. HSE06), enforce ERFC
 if (PRESENT(shortrange)) then
   if (shortrange) then
      icutcoul_local=5
   end if
 end if
! -------------------------------------

 if (icutcoul_local == 0) vcut%mode = 'SPHERE'
 if (icutcoul_local == 1) vcut%mode = 'CYLINDER'
 if (icutcoul_local == 2) vcut%mode = 'SLAB'
 if (icutcoul_local == 3) vcut%mode = 'CRYSTAL'
 if (icutcoul_local == 4) vcut%mode = 'ERF'
 if (icutcoul_local == 5) vcut%mode = 'ERFC'
 if (icutcoul_local == 6) vcut%mode = 'AUXILIARY_FUNCTION'
 if (icutcoul_local == 7) vcut%mode = 'AUX_GB'

!Initialize a few quantities
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 cutoff = gsqcut * tolfix
 barev=zero

!Some peculiar values of q: q=0 or q on the BZ edge
 qeq0=0; if (qpoint(1)**2+qpoint(2)**2+qpoint(3)**2<1.d-15) qeq0=1
 qeq05=0
 if (qeq0==0) then
   if (abs(abs(qpoint(1))-half)<tol12.or.abs(abs(qpoint(2))-half)<tol12.or. &
&   abs(abs(qpoint(3))-half)<tol12) qeq05=1
 end if

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...

 ABI_MALLOC(gq,(3,max(n1,n2,n3)))
 ABI_MALLOC(gpq,(nfft))
 ABI_MALLOC(gpq2,(nfft))

 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig+qpoint(ii)
   end do
 end do
 ig1max=-1;ig2max=-1;ig3max=-1
 ig1min=n1;ig2min=n2;ig3min=n3

 id1=n1/2+2;id2=n2/2+2;id3=n3/2+2

 ! Triple loop on each dimension
 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
   ! Precompute some products that do not depend on i2 and i1
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2
   do i2=1,n2
     ig2=i2-(i2/id2)*n2-1
     i23=n1*(i2-1 +(n2)*(i3-1))
     gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
     gqgm12=gq(2,i2)*gmet(1,2)*2
     gqg2p3=gqgm13+gqgm12

     do i1=1,n1
        ii=i1+i23
        gpq(ii)= gs2 + gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
        if (gpq(ii)>=tol4) then
          ! gpq2 contains 4*pi / |q+G|**2
          gpq2(ii) = piinv/gpq(ii)
        end if
     end do

     !
     ! Next part looks for ig1min,ig1max that are needed by zerosym
     ! Do the test that eliminates the Gamma point outside of the inner loop
     ii1=1
     if (i23==0 .and. qeq0==1  .and. ig2==0 .and. ig3==0) then
       ii1=2
     end if

     ! Final inner loop on the first dimension (note the lower limit)
     do i1=ii1,n1
       gs = gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)

       !ii=i1+i23

       if (gs<=cutoff) then
         ! Identify min/max indexes (to cancel unbalanced contributions later)
         ! Count (q+g)-vectors with similar norm
         if ((qeq05==1).and.(izero==1)) then
           ig1=i1-(i1/id1)*n1-1
           ig1max=max(ig1max,ig1); ig1min=min(ig1min,ig1)
           ig2max=max(ig2max,ig2); ig2min=min(ig2min,ig2)
           ig3max=max(ig3max,ig3); ig3min=min(ig3min,ig3)
         end if

       end if ! Cut-off
     end do ! End loop on i1
   end do
 end do


 barev(:)=zero

 ! MG: This triggers SIGFPE as cryst is not initialized
 !a1=Cryst%rprimd(:,1); b1=two_pi*gprimd(:,1)
 !a2=Cryst%rprimd(:,2); b2=two_pi*gprimd(:,2)
 !a3=Cryst%rprimd(:,3); b3=two_pi*gprimd(:,3)

 select case(TRIM(vcut%mode))
 case ('CRYSTAL', 'AUXILIARY_FUNCTION', "AUX_GB")
   if (vcut%mode == "CRYSTAL") then
     ! Analytic integration of 4pi/q^2 over the volume element:
     ! $4pi/V \int_V d^3q 1/q^2 =4pi bz_geometric_factor V^(-2/3)$
     ! i_sz=4*pi*bz_geometry_factor*q0_vol**(-two_thirds) where q0_vol= V_BZ/N_k
     ! bz_geometry_factor: sphere=7.79, fcc=7.44, sc=6.188, bcc=6.946, wz=5.255 (see gwa.pdf, appendix A.4)
     q0_vol = two_pi**3 / (nkpt_bz * ucvol)
     vcut%i_sz = four_pi*7.44*q0_vol**(-two_thirds)

   !TODO FBruneval: cryst is not available here, find a workaround!
   !else if (vcut%mode == "AUXILIARY_FUNCTION") then
   !  ! Numerical integration of the exact-exchange divergence through the
   !  ! auxiliary function of Carrier et al. PRB 75, 205126 (2007) [[cite:Carrier2007]].
   !  vcut%i_sz = carrier_isz(cryst, 1, qpoint, rcut, comm)

   !else if (vcut%mode == "AUX_GB") then
   !  ! We use the auxiliary function of a Gygi-Baldereschi variant [[cite:Gigy1986]]
   !  vcut%i_sz = gygi_baldereschi_isz(cryst, 1, qpoint, vc_ecut, ng, gvec)

   else
     ABI_ERROR(sjoin("Need treatment of 1/q^2 singularity! for mode", vcut%mode))
   end if

   do ig=1,nfft
     if (abs(gpq(ig))<tol4) then
       barev(ig) = vcut%i_sz
     else if (gpq(ig)<=cutoff) then
       barev(ig) = gpq2(ig)
     end if
   end do
   

 case('SPHERE') ! Spherical cutoff

   !
   ! Treatment of the divergence at q+g=zero
   !
   ! rcut is not set (rcut<=0), use the default Spencer-Alavi definition: 
   if ( rcut < tol8 ) then
     rcut0= (three*nkpt_bz*ucvol/four_pi)**(one/three)
   else
     rcut0 = rcut
   end if

   do ig=1,nfft
     if (abs(gpq(ig))<tol4) then
       barev(ig) = two_pi*rcut0**two
     else if (gpq(ig)<=cutoff) then
       barev(ig) = gpq2(ig) * (one - cos( rcut0*sqrt(four_pi/gpq2(ig)) ) )
     end if
   end do

 case('CYLINDER')
   !FBruneval: not working. For instance, Cryst is never initialized
   ABI_BUG("Cylinder cutoff coding is not finalized")

   test=COUNT(ABS(vcut%vcutgeo)>tol6)
   ABI_CHECK(test==1,'Wrong cutgeo for cylinder')

   ! === Beigi method is the default one, i.e infinite cylinder of radius rcut ===
   ! * Negative values to use Rozzi method with finite cylinder of extent hcyl.
   opt_cylinder=1; vcut%hcyl=zero; vcut%pdir(:)=0
   do ii=1,3
     check=vcut%vcutgeo(ii)
     if (ABS(check)>tol6) then
       vcut%pdir(ii)=1
       if (check<zero) then  ! use Rozzi's method.
         vcut%hcyl=ABS(check)*NORM2(Cryst%rprimd(:,ii))
         opt_cylinder=2
       end if
     end if
   end do

   test=COUNT(vcut%pdir==1)
   ABI_CHECK((test==1),'Wrong pdir for cylinder')
   if (vcut%pdir(3)/=1) then
     ABI_ERROR("The cylinder must be along the z-axis")
   end if

   ABI_BUG("cutoff cylinder API has changed!")

!   call cutoff_cylinder(nfft,gq,ng,Gsph%gvec,vcut%rcut,vcut%hcyl,vcut%pdir,&
!&                       vcut%boxcenter,Cryst%rprimd,barev,opt_cylinder,comm)

   ! === If Beigi, treat the limit q--> 0 ===
   if (opt_cylinder==1) then
     npar=8; npt=100 ; gamma_pt=RESHAPE((/0,0,0/),(/3,1/))
     ABI_MALLOC(qfit,(3,npt))
     ABI_MALLOC(vcfit,(1,npt))
     if (nfft==1) then
       ABI_ERROR("nfft == 1 not supported when Beigi's method is used")
     endif
     qfit(:,:)=zero
     step=half/(npt*(nfft-1))              ; qfit(3,:)=arth(tol6,step,npt)

     !call cutoff_cylinder(npt,qfit,1,gamma_pt,vcut%rcut,vcut%hcyl,vcut%pdir,&
     !                    vcut%boxcenter,Cryst%rprimd,vcfit,opt_cylinder,comm)

     ABI_MALLOC(xx,(npt))
     ABI_MALLOC(yy,(npt))
     ABI_MALLOC(sigma,(npt))
     ABI_MALLOC(par,(npar))
     ABI_MALLOC(var,(npar))
     ABI_MALLOC(cov,(npar,npar))
     do ii=1,npt
      xx(ii)=normv(qfit(:,ii),gmet,'G')
     end do
     ABI_FREE(qfit)
     sigma=one ; yy(:)=vcfit(1,:)
     ABI_FREE(vcfit)

     bz_plane=l2norm(b1.x.b2)
     dx=(xx(2)-xx(1))
     integ=yy(2)*dx*3.0/2.0
     integ=integ + SUM(yy(3:npt-2))*dx
     integ=integ+yy(npt-1)*dx*3.0/2.0
     write(std_out,*)' simple integral',integ
     q0_volsph=(two_pi)**3/(nkpt_bz*ucvol)
     q0_vol=bz_plane*two*xx(npt)
     write(std_out,*)' q0 sphere : ',q0_volsph,' q0_vol cyl ',q0_vol
     vcut%i_sz=bz_plane*two*integ/q0_vol
     write(std_out,*)' spherical approximation ',four_pi*7.44*q0_volsph**(-two_thirds)
     write(std_out,*)' Cylindrical cutoff value ',vcut%i_sz

     ABI_FREE(xx)
     ABI_FREE(yy)
     ABI_FREE(sigma)
     ABI_FREE(par)
     ABI_FREE(var)
     ABI_FREE(cov)

   else
     ! In Rozzi"s method the lim q+G --> 0 is finite.
     vcut%i_sz=barev(1)
   end if

 CASE('SLAB')
   !FBruneval: not working. For instance, Cryst is never initialized
   ABI_BUG("Slab cutoff coding is not finalized")

   test=COUNT(vcut%vcutgeo/=zero)
   ABI_CHECK(test==2,"Wrong vcutgeo")
   !
   ! Two methods available
   !
   ! === Default is Beigi"s method ===
   opt_slab=1; vcut%alpha(:)=zero
   if (ANY(vcut%vcutgeo<zero)) opt_slab=2
   vcut%pdir(:)=zero
   do ii=1,3
     check=vcut%vcutgeo(ii)
     if (ABS(check)>zero) then ! Use Rozzi"s method with a finite slab along x-y
       vcut%pdir(ii)=1
       if (check<zero) vcut%alpha(ii)=normv(check*Cryst%rprimd(:,ii),rmet,'R')
     end if
   end do

   ! Beigi"s method: the slab must be along x-y and R must be L_Z/2.
   if (opt_slab==1) then
     ABI_CHECK(ALL(vcut%pdir == (/1,1,0/)),"Surface must be in the x-y plane")
     !vcut%rcut = half*SQRT(DOT_PRODUCT(a3,a3))
   end if

   ABI_BUG("cutoff surface API has changed!")
   !call cutoff_slab(nfft,gq,ng,Gsph%gvec,gprimd,vcut%rcut,&
   !   vcut%boxcenter,vcut%pdir,vcut%alpha,barev,opt_slab)

   !
   ! === If Beigi, treat the limit q--> 0 ===
   if (opt_slab==1) then
     ! Integrate numerically in the plane close to 0
     npt=100 ! Number of points in 1D
     gamma_pt=RESHAPE((/0,0,0/),(/3,1/)) ! Gamma point
     ABI_MALLOC(qfit,(3,npt))
     ABI_MALLOC(qcart,(3,npt))
     ABI_MALLOC(vcfit,(1,npt))
     if (nfft==1) then
       ABI_ERROR("nfft == 1 not supported when Beigi's method is used")
     endif
     qfit(:,:)=zero
     qcart(:,:)=zero
     ! Size of the third vector
     bz_plane=l2norm(b3)
     q0_volsph=(two_pi)**3/(nkpt_bz*ucvol)
     ! radius that gives the same volume as q0_volsph
     ! Let's assume that c is perpendicular to the plane
     ! We also assume isotropic BZ around gamma
     step=sqrt((q0_volsph/bz_plane)/pi)/npt

     !step=half/(npt*(Qmesh%nibz-1))
     ! Let's take qpoints along 1 line, the vcut does depend only on the norm
     qcart(1,:)=arth(tol6,step,npt)

     do ii = 1,npt
       qfit(:,ii) = MATMUL(TRANSPOSE(Cryst%rprimd),qcart(:,ii))/(2*pi)
     end do

     ABI_BUG("cutoff surface API has changed!")

!     call cutoff_slab(npt,qfit,1,gamma_pt,gprimd,vcut%rcut,&
!       vcut%boxcenter,vcut%pdir,vcut%alpha,vcfit,opt_slab)

     ABI_MALLOC(xx,(npt))
     ABI_MALLOC(yy,(npt))
     ABI_MALLOC(sigma,(npt))
     do ii=1,npt
      !xx(ii)=qfit(1,:)
      xx(ii)=normv(qfit(:,ii),gmet,'G')
     end do
     ABI_FREE(qfit)
     sigma=one
     yy(:)=vcfit(1,:)
     !yy(:)=one
     ABI_FREE(vcfit)
     dx=(xx(2)-xx(1))
     ! integ = \int dr r f(r)
     integ=xx(2)*yy(2)*dx*3.0/2.0
     integ=integ + DOT_PRODUCT(xx(3:npt-2),yy(3:npt-2))*dx
     integ=integ+xx(npt-1)*yy(npt-1)*dx*3.0/2.0
     write(std_out,*)' simple integral',integ
     q0_vol=bz_plane*pi*xx(npt)**2
     write(std_out,*)' q0 sphere : ',q0_volsph,' q0_vol cyl ',q0_vol
     vcut%i_sz=bz_plane*2*pi*integ/q0_vol
     write(std_out,*)' spherical approximation ',four_pi*7.44*q0_volsph**(-two_thirds)
     write(std_out,*)' Cylindrical cutoff value ',vcut%i_sz

     ABI_FREE(xx)
     ABI_FREE(yy)
   else
     ! In Rozzi"s method the lim q+G --> 0 is finite.
     vcut%i_sz=barev(1)
   end if

 CASE('ERF')

   do ig=1,nfft
     if (abs(gpq(ig))<tol4) then
        !FIXME FBruneval check this value, ERFC value was wrong, so why not this one.
        barev(ig) = zero ! Stupid definition to remember something should be done here.
     else if (gpq(ig)<=cutoff) then
       !FIXME FBruneval shortrange does not make sense here (it is an optional argument that may not be present)
       ! and ERF is long range any way
       if (shortrange) then
         barev(ig) = + gpq2(ig) * exp( -pi * rcut**2 /gpq2(ig) )
       end if
    end if
   end do

 CASE('ERFC')

   do ig=1,nfft
     if (abs(gpq(ig))<tol4) then
        !FBruneval there was a wrong value here
        barev(ig) = pi * rcut**2
     else if (gpq(ig)<=cutoff) then
       ! gpq2 is 4 pi / (q+G)**2
       ! 4 pi / (q+G)**2 * [ 1 - exp( -1/4 * Rc**2 * (q+G)**2 ) ]
       barev(ig) = gpq2(ig) * ( one - exp( -pi * rcut**2 / gpq2(ig) ) )
    end if
   end do

 case default
   write(msg,'(3a)')'No cut-off applied to the Coulomb Potential.', ch10, &
                    'Either icutcoul value not allowed or not defined.'
   ABI_WARNING(msg)
 end select

 if (izero==1) then
   ! Set contribution of unbalanced components to zero
   if (qeq0==1) then !q=0
     call zerosym(barev,cplex1,n1,n2,n3)
   else if (qeq05==1) then
     !q=1/2; this doesn't work in parallel
     ig1=-1;if (mod(n1,2)==0) ig1=1+n1/2
     ig2=-1;if (mod(n2,2)==0) ig2=1+n2/2
     ig3=-1;if (mod(n3,2)==0) ig3=1+n3/2
     if (abs(abs(qpoint(1))-half)<tol12) then
       if (abs(ig1min)<abs(ig1max)) ig1=abs(ig1max)
       if (abs(ig1min)>abs(ig1max)) ig1=n1-abs(ig1min)
     end if
     if (abs(abs(qpoint(2))-half)<tol12) then
       if (abs(ig2min)<abs(ig2max)) ig2=abs(ig2max)
       if (abs(ig2min)>abs(ig2max)) ig2=n2-abs(ig2min)
     end if
     if (abs(abs(qpoint(3))-half)<tol12) then
       if (abs(ig3min)<abs(ig3max)) ig3=abs(ig3max)
       if (abs(ig3min)>abs(ig3max)) ig3=n3-abs(ig3min)
     end if
     call zerosym(barev,cplex1,n1,n2,n3,ig1=ig1,ig2=ig2,ig3=ig3)
   end if
 end if

 ABI_FREE(gq)
 ABI_FREE(gpq)
 ABI_FREE(gpq2)

end subroutine barevcoul
!!***

end module m_barevcoul
!!***
