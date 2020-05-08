!!****m* ABINIT/m_barevcoul
!! NAME
!!  m_barevcoul
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group ()
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

module m_barevcoul

 use defs_basis
 use m_dtset
 use m_errors
 use m_xmpi
 use m_fstrings,        only : sjoin, itoa
 use m_profiling_abi,   only : abimem_record
 use defs_abitypes,     only : MPI_type
 use m_numeric_tools,   only : arth, l2norm, OPERATOR(.x.),quadrature

 use m_geometry,        only : normv, metric

 use m_crystal,         only : crystal_t
 use m_gsphere,         only : gsphere_t
 use m_bz_mesh,         only : kmesh_t,kmesh_init

! Cut-off methods modules 
 use m_cutoff_sphere,   only : cutoff_sphere
 use m_cutoff_surface,  only : cutoff_surface
 use m_cutoff_cylinder, only : cutoff_cylinder, K0cos

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
   ! Lenght of the finite surface

  real(dp) :: rcut
   ! Cutoff radius

  real(dp) :: i_sz
   ! Value of the integration of the Coulomb singularity 4\pi/V_BZ \int_BZ d^3q 1/q^2

  real(dp) :: hcyl
   ! Length of the finite cylinder along the periodic dimension

  real(dp) :: ucvol
    ! Volume of the unit cell

  character(len=50) :: mode
   ! String defining the cutoff mode, possible values are: sphere,cylinder,surface,crystal

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
! private variables used for the integration needed by the cylindrical case.
 integer,save  :: npts_,ntrial_,qopt_
 real(dp),save :: ha_,hb_,r0_
 real(dp),save :: gcart_para_,gcartx_,gcarty_
 real(dp),save :: xx_
 real(dp),save :: accuracy_
  
 
contains
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/barevcoul
!! NAME
!! barevcoul
!!
!! FUNCTION
!! Compute bare coulomb term in G-space on the FFT mesh i.e. 4pi/(G+q)**2
!!
!! INPUTS
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  gsqcut=cutoff value on G**2 for sphere inside fft box. (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2))
!!  icutcoul=Option for the Coulomb potential cutoff technique
!!  divgq0= value of the integration of the Coulomb singularity 4pi\int_BZ 1/q^2 dq. Used if q = Gamma
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  izero=if 1, unbalanced components of V(q,g) are set to zero # Used by the PAW library
!!  nfft=Total number of FFT grid points.
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  barev(nfft)=4pi/(G+q)**2, G=0 component is set to divgq0/pi if q = Gamma.
!!
!! NOTES
!!  This routine operates on the full FFT mesh. DO NOT PASS MPI_TYPE
!!  One can easily implemente MPI-FFT by just calling this routine and then
!!  extracting the G-vectors treated by the node.
!!
!! PARENTS
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine barevcoul(rcut,qphon,gsqcut,gmet,nfft,nkpt_bz,ngfft,ucvol,barev,shortrange)

!Arguments ------------------------------------
!scalars
 integer,intent(in)         :: nfft,nkpt_bz
 real(dp),intent(in)        :: rcut,gsqcut,ucvol
 logical,intent(in),optional:: shortrange
!arrays
 integer,intent(in)         :: ngfft(18)
 integer                    :: ng!!!!
 real(dp),intent(in)        :: qphon(3)
 real(dp),intent(inout)     :: gmet(3,3)
 real(dp),intent(inout)     :: barev(nfft)
 real(dp)                   :: a1(3),a2(3),a3(3)
 real(dp)                   :: b1(3),b2(3),b3(3),gprimd(3,3),rmet(3,3)
 type(dataset_type)         :: dtset
 type(MPI_type)             :: mpi_enreg   !!!!
 type(crystal_t)            :: Cryst       !!!!
 type(gsphere_t)            :: Gsph
 type(vcut_t)               :: vcut        !!!!
!Local variables-------------------------------
!scalars
 integer,parameter    :: empty(3,3)=zero
 integer              :: comm
 integer              :: i1,i2,i23,i3,id1,id2,id3,icutcoul_local
 integer              :: ig,ig1min,ig1max,ig2min,ig2max,ig3min,ig3max
 integer              :: ii,ing,n1,n2,n3,npar,npt
 integer              :: opt_cylinder,opt_surface,test
 real(dp),parameter   :: tolfix=1.000000001e0_dp ! Same value as the one used in hartre
 real(dp)             :: check,step
 real(dp)             :: cutoff,gqg2p3,gqgm12,gqgm13,gqgm23,gs2,gs3,divgq0,rcut0
 real(dp)             :: bz_plane,dx,integ,q0_vol,q0_volsph
 character(len=500)   :: msg
!arrays
 integer              :: id(3), gamma_pt(3,1)
 real(dp),allocatable :: gq(:,:),gpq(:),gpq2(:)
 real(dp),allocatable :: vcfit(:,:),xx(:),yy(:)
 real(dp),allocatable :: cov(:,:),par(:),qfit(:,:),sigma(:),var(:),qcart(:,:)
!
 comm=mpi_enreg%comm_world
!
! === Save dimension and other useful quantities in barev% ===
 vcut%nfft      = PRODUCT(ngfft(1:3))  ! Number of points in the FFT mesh.
 vcut%ng        = ng                   ! Number of G-vectors in the Coulomb matrix elements.
 vcut%rcut      = rcut                 ! Cutoff radius for cylinder.
 vcut%hcyl      = zero                 ! Length of finite cylinder (Rozzi"s method, default is Beigi).
 vcut%ucvol     = ucvol                ! Unit cell volume.

 vcut%rprimd    = Cryst%rprimd(:,:)    ! Dimensional direct lattice.
 vcut%boxcenter = dtset%boxcenter      ! boxcenter at the moment is supposed to be at the origin.
 vcut%vcutgeo   = dtset%vcutgeo(:)     ! Info on the orientation and extension of the cutoff region.
!
! === Define geometry and cutoff radius (if used) ===
 vcut%mode='NONE'
 !icutcoul_local=dtset%icutcoul

! BG: Temporary to circumvent the tests 
 if(shortrange) then 
    icutcoul_local=5
 else 
    icutcoul_local=0
 end if
! -------------------------------------

 if (icutcoul_local==0) vcut%mode='SPHERE'
 if (icutcoul_local==1) vcut%mode='CYLINDER'
 if (icutcoul_local==2) vcut%mode='SURFACE'
 if (icutcoul_local==4) vcut%mode='ERF'
 if (icutcoul_local==5) vcut%mode='ERFC'
!
! Treatment of the divergence at q+g=zero
 rcut0= (three*nkpt_bz*ucvol/four_pi)**(one/three)
 divgq0= two_pi*rcut0**two

!Initialize a few quantities
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 cutoff=gsqcut*tolfix
 barev=zero

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...
 
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
 ABI_ALLOCATE(gpq,(nfft))
 ABI_ALLOCATE(gpq2,(nfft))
 
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig+qphon(ii)
   end do
 end do
 ig1max=-1;ig2max=-1;ig3max=-1
 ig1min=n1;ig2min=n2;ig3min=n3

 id1=n1/2+2;id2=n2/2+2;id3=n3/2+2

 do i3=1,n3
   ! Precompute some products that do not depend on i2 and i1
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2
   do i2=1,n2
     i23=n1*(i2-1 +(n2)*(i3-1))
     gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
     gqgm12=gq(2,i2)*gmet(1,2)*2
     gqg2p3=gqgm13+gqgm12
     do i1=1,n1
        ii=i1+i23
        gpq(ii)=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
        if(gpq(ii)>=tol4) then 
          gpq2(ii) = piinv/gpq(ii)
        end if 
     end do
   end do
 end do

! Old version of the code extracted from m_Fock
! do ig=1,nfft 
!     if(abs(gpq(ig))<tol4) then 
!        barev(ig)=barev(ig)+divgq0
!     else if(gpq(ig)<=cutoff) then
!       if(shortrange) then
!         barev(ig)=barev(ig)+gpq2(ig)*(one-exp(-pi/(gpq2(ig)*rcut**2)))
!       else
!         barev(ig)=barev(ig)+gpq2(ig)*(one-cos(rcut*sqrt(four_pi/gpq2(ig))))
!       end if
!    end if
! end do

 barev(:)=zero

 a1=Cryst%rprimd(:,1); b1=two_pi*gprimd(:,1)
 a2=Cryst%rprimd(:,2); b2=two_pi*gprimd(:,2)
 a3=Cryst%rprimd(:,3); b3=two_pi*gprimd(:,3)

 SELECT CASE (TRIM(vcut%mode))
 CASE('SPHERE') ! Spencer-Alavi method

   do ig=1,nfft
     if(abs(gpq(ig))<tol4) then
        barev(ig)=barev(ig)+divgq0
     else if(gpq(ig)<=cutoff) then
         barev(ig)=barev(ig)+gpq2(ig)*(one-cos(rcut*sqrt(four_pi/gpq2(ig))))
    end if
   end do

 CASE('CYLINDER')

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
         vcut%hcyl=ABS(check)*SQRT(SUM(Cryst%rprimd(:,ii)**2))
         opt_cylinder=2
       end if
     end if
   end do

   test=COUNT(vcut%pdir==1)
   ABI_CHECK((test==1),'Wrong pdir for cylinder')
   if (vcut%pdir(3)/=1) then
     MSG_ERROR("The cylinder must be along the z-axis")
   end if

   call cutoff_cylinder(nfft,gq,ng,Gsph%gvec,vcut%rcut,vcut%hcyl,vcut%pdir,&
&                       vcut%boxcenter,Cryst%rprimd,barev,opt_cylinder,comm)

   ! === If Beigi, treat the limit q--> 0 ===
   if (opt_cylinder==1) then
     npar=8; npt=100 ; gamma_pt=RESHAPE((/0,0,0/),(/3,1/))
     ABI_MALLOC(qfit,(3,npt))
     ABI_MALLOC(vcfit,(1,npt))
     if (nfft==1) then
       MSG_ERROR("nfft == 1 not supported when Beigi's method is used")
     endif
     qfit(:,:)=zero
     step=half/(npt*(nfft-1))              ; qfit(3,:)=arth(tol6,step,npt)

     call cutoff_cylinder(npt,qfit,1,gamma_pt,vcut%rcut,vcut%hcyl,vcut%pdir,&
&                         vcut%boxcenter,Cryst%rprimd,vcfit,opt_cylinder,comm)

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

     do ii=3,npt-2
       integ=integ+yy(ii)*dx
     end do
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

 CASE('SURFACE')

   test=COUNT(vcut%vcutgeo/=zero)
   ABI_CHECK(test==2,"Wrong vcutgeo")
   !
   ! Two methods available
   !
   ! === Default is Beigi"s method ===
   opt_surface=1; vcut%alpha(:)=zero
   if (ANY(vcut%vcutgeo<zero)) opt_surface=2
   vcut%pdir(:)=zero
   do ii=1,3
     check=vcut%vcutgeo(ii)
     if (ABS(check)>zero) then ! Use Rozzi"s method with a finite surface along x-y
       vcut%pdir(ii)=1
       if (check<zero) vcut%alpha(ii)=normv(check*Cryst%rprimd(:,ii),rmet,'R')
     end if
   end do

   ! Beigi"s method: the surface must be along x-y and R must be L_Z/2.
   if (opt_surface==1) then
     ABI_CHECK(ALL(vcut%pdir == (/1,1,0/)),"Surface must be in the x-y plane")
     vcut%rcut = half*SQRT(DOT_PRODUCT(a3,a3))
   end if

   call cutoff_surface(nfft,gq,ng,Gsph%gvec,gprimd,vcut%rcut,&
&    vcut%boxcenter,vcut%pdir,vcut%alpha,barev,opt_surface)

   !
   ! === If Beigi, treat the limit q--> 0 ===
   if (opt_surface==1) then
     ! Integrate numerically in the plane close to 0
     npt=100 ! Number of points in 1D
     gamma_pt=RESHAPE((/0,0,0/),(/3,1/)) ! Gamma point
     ABI_MALLOC(qfit,(3,npt))
     ABI_MALLOC(qcart,(3,npt))
     ABI_MALLOC(vcfit,(1,npt))
     if (nfft==1) then
       MSG_ERROR("nfft == 1 not supported when Beigi's method is used")
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

     call cutoff_surface(npt,qfit,1,gamma_pt,gprimd,vcut%rcut,&
&       vcut%boxcenter,vcut%pdir,vcut%alpha,vcfit,opt_surface)

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
     do ii=3,npt-2
       integ=integ+xx(ii)*yy(ii)*dx
     end do
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
     if(abs(gpq(ig))<tol4) then
        barev(ig)=barev(ig)+divgq0
     else if(gpq(ig)<=cutoff) then
       if(shortrange) then
         barev(ig)=barev(ig)+gpq2(ig)*exp(-pi/(gpq2(ig)*rcut**2))
       end if
    end if
   end do
 
 CASE('ERFC')

   do ig=1,nfft
     if(abs(gpq(ig))<tol4) then
        barev(ig)=barev(ig)+divgq0
     else if(gpq(ig)<=cutoff) then
       if(shortrange) then
         barev(ig)=barev(ig)+gpq2(ig)*(one-exp(-pi/(gpq2(ig)*rcut**2)))
       end if
    end if
   end do

 CASE DEFAULT
   write(msg,'(a,i3)')'No cut-off applied to the Coulomb Potential.' //&
&                     'Either icutcoul value not allowed or not defined.'
   MSG_WARNING(msg)
 END SELECT

 ABI_DEALLOCATE(gq)
 ABI_DEALLOCATE(gpq)
 ABI_DEALLOCATE(gpq2)
  
end subroutine barevcoul
!!***

end module m_barevcoul
!!***
