!!****m* ABINIT/m_vcoul
!! NAME
!!  m_vcoul
!!
!! FUNCTION
!!  This module contains the definition of the vcoul_t as well
!!  as procedures to calculate the Coulomb interaction in reciprocal
!!  space taking into account a possible cutoff in real space.
!!  Procedures to deal with the singularity for q-->0 are also provided.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (MG, FB)
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

MODULE m_vcoul

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_splines
 use m_sort

 use m_fstrings,        only : sjoin, itoa
 use m_special_funcs,   only : abi_derf
 use m_gwdefs,          only : GW_TOLQ0
 use m_io_tools,        only : open_file
 use m_numeric_tools,   only : arth, geop, imin_loc, llsfit_svd, l2norm, OPERATOR(.x.), quadrature, isdiagmat
 use m_bessel,          only : CALJY0, CALJY1, CALCK0, CALCK1
 use m_hide_lapack,     only : matrginv
 use m_geometry,        only : normv, metric
 use m_qplusg,          only : cmod_qpg
 use m_crystal,         only : crystal_t
 use m_bz_mesh,         only : kmesh_t, get_BZ_item
 use m_gsphere,         only : gsphere_t
 use m_paw_numeric,     only : paw_jbessel
 use m_dtfil,           only : isfile

! Cut-off methods modules 
 use m_cutoff_sphere,   only : cutoff_sphere
 use m_cutoff_surface,  only : cutoff_surface
 use m_cutoff_cylinder, only : cutoff_cylinder, K0cos

 implicit none

 private
!!***

!!****t* m_vcoul/vcoul_t
!! NAME
!!  vcoul_t
!!
!! FUNCTION
!!  This data type contains the square root of the Fourier components of the Coulomb interaction
!!  calculated taking into account a possible cutoff. Info on the particular geometry used for the cutoff
!!  as well as quantities required to deal with the Coulomb divergence.
!!
!! SOURCE

 type,public :: vcoul_t

  ! TODO: Remove it
  integer  :: nfft
  ! Number of points in FFT grid

  integer  :: ng
   ! Number of G-vectors

  integer  :: nqibz
   ! Number of irreducible q-points

  integer  :: nqlwl
   ! Number of small q-points around Gamma

  real(dp) :: alpha(3)
   ! Lenght of the finite surface

  real(dp) :: rcut
   ! Cutoff radius

  real(dp) :: i_sz
   ! Value of the integration of the Coulomb singularity 4\pi/V_BZ \int_BZ d^3q 1/q^2

  real(dp) :: i_sz_resid
   ! Residual difference between the i_sz in the sigma self-energy for exchange,
   ! and the i_sz already present in the generalized Kohn-Sham eigenenergies
   ! Initialized to the same value as i_sz

  real(dp) :: hcyl
   ! Length of the finite cylinder along the periodic dimension

  real(dp) :: ucvol
    ! Volume of the unit cell

  character(len=50) :: mode
   ! String defining the cutoff mode, possible values are: sphere,cylinder,surface,crystal

  integer :: pdir(3)
   ! 1 if the system is periodic along this direction

  ! TODO: Remove it
  integer :: ngfft(18)
    ! Information on the FFT grid

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

  real(dp),allocatable :: qlwl(:,:)
   ! qibz(3,nqlwl)
   ! q-points for the treatment of the Coulomb singularity.

  complex(gwpc),allocatable :: vc_sqrt(:,:)
    ! vc_sqrt(ng,nqibz)
    ! Square root of the Coulomb interaction in reciprocal space.
    ! A cut might be applied.

  complex(gwpc),allocatable :: vcqlwl_sqrt(:,:)
    ! vcqs_sqrt(ng,nqlwl)
    ! Square root of the Coulomb term calculated for small q-points

  complex(gwpc),allocatable :: vc_sqrt_resid(:,:)
    ! vc_sqrt_resid(ng,nqibz)
    ! Square root of the residual difference between the Coulomb interaction in the sigma self-energy for exchange,
    ! and the Coulomb interaction already present in the generalized Kohn-Sham eigenenergies (when they come from an hybrid)
    ! Given in reciprocal space. At the call to vcoul_init, it is simply initialized at the value of vc_sqrt(:,:),
    ! and only later modified.
    ! A cut might be applied.

 end type vcoul_t


 public ::  vcoul_init           ! Main creation method.
 public ::  vcoul_plot           ! Plot vc in real and reciprocal space.
 public ::  vcoul_print          ! Report info on the object.
 public ::  vcoul_free           ! Destruction method.
! public ::  cmod_qpg             ! FT of the long ranged Coulomb interaction.
!!***

! private variables used for the integration needed by the cylindrical case.
 integer,save :: npts_,ntrial_,qopt_
 real(dp),save :: hb_,r0_
 real(dp),save :: qpg_para_,qpgx_,qpgy_
 real(dp),save :: xx_,rho_
 real(dp),save :: accuracy_

CONTAINS  !========================================================================================
!!***

!!****f* m_vcoul/vcoul_init
!! NAME
!! vcoul_init
!!
!! FUNCTION
!! Perform general check and initialize the data type containing information on the cutoff technique
!! Note Vcp%vc_sqrt_resid and Vcp%i_sz_resid are simply initialized at the same value as Vcp%vc_sqrt and Vcp%i_sz
!!
!! INPUTS
!!  Qmesh<kmesh_t>=Info on the q-point sampling.
!!  Kmesh<kmesh_t>=Info on the k-point sampling.
!!  Gsph<gsphere_t>=Info of the G sphere.
!!    %gmet(3,3)=Metric in reciprocal space.
!!    %gprimd(3,3)=Dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!    %gvec=G vectors
!!  rcut=Cutoff radius for the cylinder.
!!  gw_icutcoul=Option of the cutoff technique.
!!  vcutgeo(3)= Info on the orientation and extension of the cutoff region.
!!  ng=Number of G-vectors to be used to describe the Coulomb interaction
!!  nqlwl=Number of point around Gamma for treatment of long-wavelength limit
!!  qlwl(3,nqlwl)= The nqlwl "small" q-points
!!  ngfft(18)=Information on the (fine) FFT grid used for the density.
!!  rprimd(3,3)=Direct lattice vectors.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Vcp=<vcoul_t>=Datatype gathering information on the Coulomb interaction.
!!
!! PARENTS
!!      gwls_hamiltonian,mrgscr,setup_bse,setup_bse_interp,setup_screening
!!      setup_sigma
!!
!! CHILDREN
!!      calck0,paw_jbessel,quadrature
!!
!! SOURCE

subroutine vcoul_init(Vcp,Gsph,Cryst,Qmesh,Kmesh,rcut,gw_icutcoul,vcutgeo,ecut,ng,nqlwl,qlwl,ngfft,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng,nqlwl,gw_icutcoul,comm
 real(dp),intent(in) :: rcut, ecut
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(gsphere_t),intent(in) :: Gsph
 type(crystal_t),intent(in) :: Cryst
 type(vcoul_t),intent(out) :: Vcp
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: qlwl(3,nqlwl),vcutgeo(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,ncell=3
 integer :: nmc_max=2500000
 integer :: nmc,nseed
 integer :: i1,i2,i3,ig,imc
 integer :: ii,iqlwl,iq_bz,iq_ibz,npar,npt
 integer :: opt_cylinder,opt_surface,test,rank,nprocs
 integer, allocatable :: seed(:)
 real(dp),parameter :: tolq0=1.d-3
 real(dp) :: b1b1,b2b2,b3b3,b1b2,b2b3,b3b1
 real(dp) :: bz_geometry_factor,bz_plane,check,dx,integ,q0_vol,q0_volsph
 real(dp) :: qbz_norm,step,ucvol,intfauxgb, alfa
 character(len=500) :: msg
!arrays
 integer :: gamma_pt(3,1)
 real(dp) :: a1(3),a2(3),a3(3),bb(3),b1(3),b2(3),b3(3),gmet(3,3),gprimd(3,3)
 real(dp) :: qbz_cart(3),rmet(3,3)
 real(dp),allocatable :: cov(:,:),par(:),qfit(:,:),sigma(:),var(:),qcart(:,:)
 real(dp),allocatable :: vcfit(:,:),vcoul(:,:),vcoul_lwl(:,:),xx(:),yy(:)
 real(dp),allocatable :: qran(:,:)
 real(dp) :: lmin,vlength,qpg2,rcut2
 real(dp) :: q0sph
 real(dp) :: qtmp(3),qmin(3),qmin_cart(3),qpg(3)
 real(dp) :: rprimd_sc(3,3),gprimd_sc(3,3),gmet_sc(3,3),rmet_sc(3,3),ucvol_sc
 real(dp) :: qcart2red(3,3)
 real(dp) :: qqgg(3)

! *************************************************************************

 DBG_ENTER("COLL")
 !
 ! === Test if the first q-point is zero ===
 ! FIXME this wont work if nqptdm/=0
 !if (normv(Qmesh%ibz(:,1),gmet,'G')<GW_TOLQ0)) STOP 'vcoul_init, non zero first point '
 rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 call metric(gmet,gprimd,-1,rmet,Cryst%rprimd,ucvol)
 !
 ! === Save dimension and other useful quantities in Vcp% ===
 Vcp%nfft      = PRODUCT(ngfft(1:3))  ! Number of points in the FFT mesh.
 Vcp%ng        = ng                   ! Number of G-vectors in the Coulomb matrix elements.
 Vcp%nqibz     = Qmesh%nibz           ! Number of irred q-point.
 Vcp%nqlwl     = nqlwl                ! Number of small q-directions to deal with singularity and non Analytic behavior.
 Vcp%rcut      = rcut                 ! Cutoff radius for cylinder.
 Vcp%hcyl      = zero                 ! Length of finite cylinder (Rozzi"s method, default is Beigi).
 Vcp%ucvol     = ucvol                ! Unit cell volume.

 Vcp%rprimd    = Cryst%rprimd(:,:)    ! Dimensional direct lattice.
 Vcp%boxcenter = zero                 ! boxcenter at the moment is supposed to be at the origin.
 Vcp%vcutgeo   = vcutgeo(:)           ! Info on the orientation and extension of the cutoff region.
 Vcp%ngfft(:)  = ngfft(:)             ! Info on the FFT mesh.
 !
 ! === Define geometry and cutoff radius (if used) ===
 Vcp%mode='NONE'
 if (gw_icutcoul==14) Vcp%mode='MINIBZ-ERF'
 if (gw_icutcoul==15) Vcp%mode='MINIBZ-ERFC'
 if (gw_icutcoul==16) Vcp%mode='MINIBZ'
 if (gw_icutcoul==0) Vcp%mode='SPHERE'
 if (gw_icutcoul==1) Vcp%mode='CYLINDER'
 if (gw_icutcoul==2) Vcp%mode='SURFACE'
 if (gw_icutcoul==3) Vcp%mode='CRYSTAL'
 if (gw_icutcoul==4) Vcp%mode='ERF'
 if (gw_icutcoul==5) Vcp%mode='ERFC'
 if (gw_icutcoul==6) Vcp%mode='AUXILIARY_FUNCTION'
 if (gw_icutcoul==7) Vcp%mode='AUX_GB'
 !
 ! === Calculate Fourier coefficient of Coulomb interaction ===
 ! * For the limit q-->0 we consider ng vectors due to a possible anisotropy in case of a cutoff interaction
 ABI_MALLOC(Vcp%qibz,(3,Vcp%nqibz))
 Vcp%qibz = Qmesh%ibz(:,:)
 ABI_MALLOC(Vcp%qlwl,(3,Vcp%nqlwl))
 Vcp%qlwl = qlwl(:,:)

 ! ===============================================
 ! == Calculation of the FT of the Coulomb term ==
 ! ===============================================
 ABI_MALLOC(vcoul    ,(ng,Vcp%nqibz))
 ABI_MALLOC(vcoul_lwl,(ng,Vcp%nqlwl))

 a1=Cryst%rprimd(:,1); b1=two_pi*gprimd(:,1)
 a2=Cryst%rprimd(:,2); b2=two_pi*gprimd(:,2)
 a3=Cryst%rprimd(:,3); b3=two_pi*gprimd(:,3)
 b1b1=dot_product(b1,b1) ; b2b2=dot_product(b2,b2) ; b3b3=dot_product(b3,b3)
 bb(1)=b1b1 ; bb(2)=b2b2 ; bb(3)=b3b3
 b1b2=dot_product(b1,b2) ; b2b3=dot_product(b2,b3) ; b3b1=dot_product(b3,b1)

 SELECT CASE (TRIM(Vcp%mode))

 CASE ('MINIBZ','MINIBZ-ERFC','MINIBZ-ERF')

   ! Mimicking the BerkeleyGW technique
   ! A Monte-Carlo sampling of each miniBZ surrounding each (q+G) point
   ! However - extended to the multiple shifts
   !         - with an adaptative number of MonteCarlo sampling points

   !
   ! Supercell defined by the k-points
   rprimd_sc(:,:) = MATMUL(Cryst%rprimd,Kmesh%kptrlatt)
   call metric(gmet_sc,gprimd_sc,-1,rmet_sc,rprimd_sc,ucvol_sc)

   qcart2red(:,:) = two_pi * Cryst%gprimd(:,:)
   call matrginv(qcart2red,3,3)

   !
   ! Find the largest sphere inside the miniBZ
   ! in order to integrate the divergence analytically
   q0sph = HUGE(1.0_dp)
   do i1 = -ncell+1, ncell
     qtmp(1) = dble(i1) * 0.5_dp
     do i2 = -ncell+1, ncell
       qtmp(2) = dble(i2) * 0.5_dp
       do i3 = -ncell+1, ncell
         qtmp(3) = dble(i3) * 0.5_dp
         if( i1==0 .AND. i2==0 .AND. i3==0 ) cycle
         vlength = normv(qtmp,gmet_sc,'G')
         if (vlength < q0sph) then
           q0sph = vlength
         end if
       enddo
     enddo
   enddo


   !
   ! Setup the random vectors for the Monte Carlo sampling of the miniBZ
   ! at q=0
   ABI_MALLOC(qran,(3,nmc_max))
   call random_seed(size=nseed)
   ABI_MALLOC(seed,(nseed))
   do i1=1,nseed
     seed(i1) = NINT(SQRT(DBLE(i1)*103731))
   end do
   call random_seed(put=seed)
   call random_number(qran)

   ! Overide the first "random vector" with 0
   qran(:,1) = 0.0_dp

   ! Fold qran into the Wignez-Seitz cell around q=0
   do imc=2,nmc_max
     lmin = HUGE(1.0_dp)
     do i1 = -ncell+1, ncell
       qtmp(1) = qran(1,imc) + dble(i1)
       do i2 = -ncell+1, ncell
         qtmp(2) = qran(2,imc) + dble(i2)
         do i3 = -ncell+1, ncell
           qtmp(3) = qran(3,imc) + dble(i3)
           vlength = normv(qtmp,gmet_sc,'G')
           if (vlength < lmin) then
             lmin = vlength
             ! Get the q-vector in cartersian coordinates
             qmin_cart(:) = two_pi * MATMUL( gprimd_sc(:,:) , qtmp )
             ! Transform it back to the reciprocal space
             qmin(:) = MATMUL( qcart2red , qmin_cart )
           end if
         enddo
       enddo
     enddo

     qran(1,imc) = qmin(1)
     qran(2,imc) = qmin(2)
     qran(3,imc) = qmin(3)

   enddo


   rcut2 = Vcp%rcut**2
   vcoul(:,:) = zero

   !
   ! Admittedly, it's a lot of duplication of code hereafter,
   ! but I think it's the clearest and fastest way
   !
   select case( TRIM(Vcp%mode) )
   case('MINIBZ')

     do iq_ibz=1,Vcp%nqibz
       do ig=1,ng
         if( iq_ibz==1 .AND. ig==1 ) cycle

         qpg(:) =  Vcp%qibz(:,iq_ibz) +  Gsph%gvec(:,ig)
         qpg2 = normv(qpg,gmet,'G')**2
         nmc = adapt_nmc(nmc_max,qpg2)
         do imc=1,nmc
           qpg(:) =  Vcp%qibz(:,iq_ibz) +  Gsph%gvec(:,ig) + qran(:,imc)
           qpg2 = normv(qpg,gmet,'G')**2
           vcoul(ig,iq_ibz) = vcoul(ig,iq_ibz) + four_pi / qpg2 / REAL(nmc,dp)
         end do
       end do
     end do

     ! Override q=1, ig=1
     vcoul(1,1) = four_pi**2 * Kmesh%nbz * ucvol / ( 8.0_dp * pi**3 ) * q0sph
     do imc=1,nmc_max
       qpg(:) =  Vcp%qibz(:,1) +  Gsph%gvec(:,1) + qran(:,imc)
       qpg2 = normv(qpg,gmet,'G')**2
       if( qpg2 > q0sph**2 ) then
         vcoul(1,1) = vcoul(1,1) + four_pi / qpg2 / REAL(nmc_max,dp)
       end if
     end do

     vcoul_lwl(:,:) = zero
     do iq_ibz=1,Vcp%nqlwl
       do ig=1,ng
         qpg(:) =  Vcp%qlwl(:,iq_ibz) +  Gsph%gvec(:,ig)
         qpg2 = normv(qpg,gmet,'G')**2
         nmc = adapt_nmc(nmc_max,qpg2)
         do imc=1,nmc
           qpg(:) =  Vcp%qlwl(:,iq_ibz) +  Gsph%gvec(:,ig) + qran(:,imc)
           qpg2 = normv(qpg,gmet,'G')**2
           vcoul_lwl(ig,iq_ibz) = vcoul_lwl(ig,iq_ibz) + four_pi / qpg2 / REAL(nmc,dp)
         end do
       end do
     end do


   case('MINIBZ-ERFC')

     do iq_ibz=1,Vcp%nqibz
       do ig=1,ng
         if( iq_ibz==1 .AND. ig==1 ) cycle

         qpg(:) =  Vcp%qibz(:,iq_ibz) +  Gsph%gvec(:,ig)
         qpg2 = normv(qpg,gmet,'G')**2
         nmc = adapt_nmc(nmc_max,qpg2)
         do imc=1,nmc
           qpg(:) =  Vcp%qibz(:,iq_ibz) +  Gsph%gvec(:,ig) + qran(:,imc)
           qpg2 = normv(qpg,gmet,'G')**2
           vcoul(ig,iq_ibz) = vcoul(ig,iq_ibz) + four_pi / qpg2 / REAL(nmc,dp) &
&                  * (  one - EXP( -0.25d0 * rcut2 * qpg2 ) )
         end do
       end do
     end do

     ! Override q=1, ig=1
     vcoul(1,1) = four_pi**2 * Kmesh%nbz * ucvol / ( 8.0_dp * pi**3 ) &
&        * ( q0sph - SQRT(pi/rcut2) * abi_derf(0.5_dp*SQRT(rcut2)*q0sph) )
     do imc=1,nmc_max
       qpg(:) =  Vcp%qibz(:,1) +  Gsph%gvec(:,1) + qran(:,imc)
       qpg2 = normv(qpg,gmet,'G')**2
       if( qpg2 > q0sph**2 ) then
         vcoul(1,1) = vcoul(1,1) + four_pi / qpg2 / REAL(nmc_max,dp) &
&                      * (  one - EXP( -0.25d0 * rcut2 * qpg2 ) )
       end if
     end do

     vcoul_lwl(:,:) = zero
     do iq_ibz=1,Vcp%nqlwl
       do ig=1,ng
         qpg(:) =  Vcp%qlwl(:,iq_ibz) +  Gsph%gvec(:,ig)
         qpg2 = normv(qpg,gmet,'G')**2
         nmc = adapt_nmc(nmc_max,qpg2)
         do imc=1,nmc
           qpg(:) =  Vcp%qlwl(:,iq_ibz) +  Gsph%gvec(:,ig) + qran(:,imc)
           qpg2 = normv(qpg,gmet,'G')**2
           vcoul_lwl(ig,iq_ibz) = vcoul_lwl(ig,iq_ibz) + four_pi / qpg2 / REAL(nmc,dp) &
&                  * (  one - EXP( -0.25d0 * rcut2 * qpg2 ) )
         end do
       end do
     end do


   case('MINIBZ-ERF')

     do iq_ibz=1,Vcp%nqibz
       do ig=1,ng
         if( iq_ibz==1 .AND. ig==1 ) cycle

         qpg(:) =  Vcp%qibz(:,iq_ibz) +  Gsph%gvec(:,ig)
         qpg2 = normv(qpg,gmet,'G')**2
         nmc = adapt_nmc(nmc_max,qpg2)
         do imc=1,nmc
           qpg(:) =  Vcp%qibz(:,iq_ibz) +  Gsph%gvec(:,ig) + qran(:,imc)
           qpg2 = normv(qpg,gmet,'G')**2
           vcoul(ig,iq_ibz) = vcoul(ig,iq_ibz) + four_pi / qpg2 / REAL(nmc,dp) &
&                  * EXP( -0.25d0 * rcut2 * qpg2 )
         end do
       end do
     end do

     ! Override q=1, ig=1
     vcoul(1,1) = four_pi**2 * Kmesh%nbz * ucvol / ( 8.0_dp * pi**3 ) &
&                  * SQRT(pi/rcut2) * abi_derf(0.5_dp*SQRT(rcut2)*q0sph)

     do imc=1,nmc_max
       qpg(:) =  Vcp%qibz(:,1) +  Gsph%gvec(:,1) + qran(:,imc)
       qpg2 = normv(qpg,gmet,'G')**2
       if( qpg2 > q0sph**2 ) then
         vcoul(1,1) = vcoul(1,1) + four_pi / qpg2 / REAL(nmc_max,dp) *  EXP( -0.25d0 * rcut2 * qpg2 )
       end if
     end do

     vcoul_lwl(:,:) = zero
     do iq_ibz=1,Vcp%nqlwl
       do ig=1,ng
         qpg(:) =  Vcp%qlwl(:,iq_ibz) +  Gsph%gvec(:,ig)
         qpg2 = normv(qpg,gmet,'G')**2
         nmc = adapt_nmc(nmc_max,qpg2)
         do imc=1,nmc
           qpg(:) =  Vcp%qlwl(:,iq_ibz) +  Gsph%gvec(:,ig) + qran(:,imc)
           qpg2 = normv(qpg,gmet,'G')**2
           vcoul_lwl(ig,iq_ibz) = vcoul_lwl(ig,iq_ibz) + four_pi / qpg2 / REAL(nmc,dp) &
&                  * EXP( -0.25d0 * rcut2 * qpg2 )
         end do
       end do
     end do

   end select

   Vcp%i_sz = vcoul(1,1)

   ABI_FREE(qran)
   ABI_FREE(seed)

 CASE ('SPHERE')
   ! TODO check that L-d > R_c > d
   ! * A non-positive value of rcut imposes the recipe of Spencer & Alavi, PRB 77, 193110 (2008) [[cite:Spencer2008]].
   if (Vcp%rcut<tol12) then
     Vcp%rcut = (ucvol*Kmesh%nbz*3.d0/four_pi)**third
     write(msg,'(2a,2x,f8.4,a)')ch10,&
&     ' Using a calculated value for rcut = ',Vcp%rcut, ' to have the same volume as the BvK crystal '
     call wrtout(std_out,msg,'COLL')
   end if

   Vcp%vcutgeo=zero
   call cutoff_sphere(Qmesh%nibz,Qmesh%ibz,ng,Gsph%gvec,gmet,Vcp%rcut,vcoul)

   ! q-points for optical limit.
   call cutoff_sphere(Vcp%nqlwl,Vcp%qlwl,ng,Gsph%gvec,gmet,Vcp%rcut,vcoul_lwl)
   !
   ! === Treat the limit q--> 0 ===
   ! * The small cube is approximated by a sphere, while vc(q=0)=2piR**2.
   ! * if a single q-point is used, the expression for the volume is exact.
   Vcp%i_sz=two_pi*Vcp%rcut**2
   call vcoul_print(Vcp,unit=ab_out)

 CASE ('CYLINDER')
   test=COUNT(ABS(Vcp%vcutgeo)>tol6)
   ABI_CHECK(test==1,'Wrong cutgeo for cylinder')

   ! === Beigi method is the default one, i.e infinite cylinder of radius rcut ===
   ! * Negative values to use Rozzi method with finite cylinder of extent hcyl.
   opt_cylinder=1; Vcp%hcyl=zero; Vcp%pdir(:)=0
   do ii=1,3
     check=Vcp%vcutgeo(ii)
     if (ABS(check)>tol6) then
       Vcp%pdir(ii)=1
       if (check<zero) then  ! use Rozzi's method.
         Vcp%hcyl=ABS(check)*SQRT(SUM(Cryst%rprimd(:,ii)**2))
         opt_cylinder=2
       end if
     end if
   end do

   test=COUNT(Vcp%pdir==1)
   ABI_CHECK((test==1),'Wrong pdir for cylinder')
   if (Vcp%pdir(3)/=1) then
     MSG_ERROR("The cylinder must be along the z-axis")
   end if

   call cutoff_cylinder(Qmesh%nibz,Qmesh%ibz,ng,Gsph%gvec,Vcp%rcut,Vcp%hcyl,Vcp%pdir,&
&    Vcp%boxcenter,Cryst%rprimd,vcoul,opt_cylinder,comm)

   ! q-points for optical limit.
   call cutoff_cylinder(Vcp%nqlwl,Vcp%qlwl,ng,Gsph%gvec,Vcp%rcut,Vcp%hcyl,Vcp%pdir,&
&    Vcp%boxcenter,Cryst%rprimd,vcoul_lwl,opt_cylinder,comm)
   !
   ! === If Beigi, treat the limit q--> 0 ===
   if (opt_cylinder==1) then
     npar=8; npt=100 ; gamma_pt=RESHAPE((/0,0,0/),(/3,1/))
     ABI_MALLOC(qfit,(3,npt))
     ABI_MALLOC(vcfit,(1,npt))
     if (Qmesh%nibz==1) then
       MSG_ERROR("nqibz == 1 not supported when Beigi's method is used")
     endif
     qfit(:,:)=zero
     step=half/(npt*(Qmesh%nibz-1))              ; qfit(3,:)=arth(tol6,step,npt)
     !step=(half/(Qmesh%nibz-1)/tol6)**(one/npt) ; qfit(3,:)=geop(tol6,step,npt)
     call cutoff_cylinder(npt,qfit,1,gamma_pt,Vcp%rcut,Vcp%hcyl,Vcp%pdir,Vcp%boxcenter,Cryst%rprimd,vcfit,opt_cylinder,comm)

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
     !call llsfit_svd(xx,yy,sigma,npar,K0fit,chisq,par,var,cov,info)
     !do ii=1,npt
     ! write(99,*)xx(ii),yy(ii),DOT_PRODUCT(par,K0fit(xx(ii),npar))
     !end do
     bz_plane=l2norm(b1.x.b2)
     !integ=K0fit_int(xx(npt),par,npar)
     !write(std_out,*)' SVD fit : chi-square',chisq
     !write(std_out,*)' fit-parameters : ',par
     !write(std_out,*)' variance ',var
     !write(std_out,*)' bz_plane ',bz_plane
     !write(std_out,*)' SCD integ ',integ
     ! Here Im assuming homogeneous mesh
     dx=(xx(2)-xx(1))
     integ=yy(2)*dx*3.0/2.0
     do ii=3,npt-2
       integ=integ+yy(ii)*dx
     end do
     integ=integ+yy(npt-1)*dx*3.0/2.0
     write(std_out,*)' simple integral',integ
     q0_volsph=(two_pi)**3/(Kmesh%nbz*ucvol)
     q0_vol=bz_plane*two*xx(npt)
     write(std_out,*)' q0 sphere : ',q0_volsph,' q0_vol cyl ',q0_vol
     Vcp%i_sz=bz_plane*two*integ/q0_vol
     write(std_out,*)' spherical approximation ',four_pi*7.44*q0_volsph**(-two_thirds)
     write(std_out,*)' Cylindrical cutoff value ',Vcp%i_sz
     !Vcp%i_sz=four_pi*7.44*q0_vol**(-two_thirds)
     ABI_FREE(xx)
     ABI_FREE(yy)
     ABI_FREE(sigma)
     ABI_FREE(par)
     ABI_FREE(var)
     ABI_FREE(cov)
   else
     ! In Rozzi"s method the lim q+G --> 0 is finite.
     Vcp%i_sz=vcoul(1,1)
   end if

   call vcoul_print(Vcp,unit=ab_out )

 CASE ('SURFACE')

   test=COUNT(Vcp%vcutgeo/=zero)
   ABI_CHECK(test==2,"Wrong vcutgeo")
   !
   ! === Default is Beigi"s method ===
   opt_surface=1; Vcp%alpha(:)=zero
   if (ANY(Vcp%vcutgeo<zero)) opt_surface=2
   Vcp%pdir(:)=zero
   do ii=1,3
     check=Vcp%vcutgeo(ii)
     if (ABS(check)>zero) then ! Use Rozzi"s method with a finite surface along x-y
       Vcp%pdir(ii)=1
       if (check<zero) Vcp%alpha(ii)=normv(check*Cryst%rprimd(:,ii),rmet,'R')
     end if
   end do

   ! Beigi"s method: the surface must be along x-y and R must be L_Z/2.
   if (opt_surface==1) then
     msg="2D geometry, Beigi method, the periodicity must be in the x-y plane. Modify vcutgeo or your geometry."
     ABI_CHECK(ALL(Vcp%pdir == (/1,1,0/)), msg)
     Vcp%rcut = half*SQRT(DOT_PRODUCT(a3,a3))
   end if

   call cutoff_surface(Qmesh%nibz,Qmesh%ibz,ng,Gsph%gvec,gprimd,Vcp%rcut,&
&    Vcp%boxcenter,Vcp%pdir,Vcp%alpha,vcoul,opt_surface)

   ! q-points for optical limit.
   call cutoff_surface(Vcp%nqlwl,Vcp%qlwl,ng,Gsph%gvec,gprimd,Vcp%rcut,&
&    Vcp%boxcenter,Vcp%pdir,Vcp%alpha,vcoul_lwl,opt_surface)
   !
   ! === If Beigi, treat the limit q--> 0 ===
   if (opt_surface==1) then
     ! Integrate numerically in the plane close to 0
     npt=100 ! Number of points in 1D
     gamma_pt=RESHAPE((/0,0,0/),(/3,1/)) ! Gamma point
     ABI_MALLOC(qfit,(3,npt))
     ABI_MALLOC(qcart,(3,npt))
     ABI_MALLOC(vcfit,(1,npt))
     if (Qmesh%nibz==1) then
       MSG_ERROR("nqibz == 1 not supported when Beigi's method is used")
     endif
     qfit(:,:)=zero
     qcart(:,:)=zero
     ! Size of the third vector
     bz_plane=l2norm(b3)
     q0_volsph=(two_pi)**3/(Kmesh%nbz*ucvol)
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

     call cutoff_surface(npt,qfit,1,gamma_pt,gprimd,Vcp%rcut,&
&       Vcp%boxcenter,Vcp%pdir,Vcp%alpha,vcfit,opt_surface)

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
     Vcp%i_sz=bz_plane*2*pi*integ/q0_vol
     write(std_out,*)' spherical approximation ',four_pi*7.44*q0_volsph**(-two_thirds)
     write(std_out,*)' Cylindrical cutoff value ',Vcp%i_sz
     !Vcp%i_sz=four_pi*7.44*q0_vol**(-two_thirds)
     ABI_FREE(xx)
     ABI_FREE(yy)
   else
     ! In Rozzi"s method the lim q+G --> 0 is finite.
     Vcp%i_sz=vcoul(1,1)
   end if

 CASE ('AUXILIARY_FUNCTION')
   !
   ! Numerical integration of the exact-exchange divergence through the
   ! auxiliary function of Carrier et al. PRB 75, 205126 (2007) [[cite:Carrier2007]].
   !
   do iq_ibz=1,Vcp%nqibz
     call cmod_qpg(Vcp%nqibz,iq_ibz,Vcp%qibz,ng,Gsph%gvec,gprimd,vcoul(:,iq_ibz))

     if (iq_ibz==1) then ! The singularity is treated using vcoul_lwl.
       vcoul(1, iq_ibz) = zero
       vcoul(2:,iq_ibz) = four_pi/vcoul(2:,iq_ibz)**2
     else
       vcoul(:,iq_ibz) = four_pi/vcoul(:,iq_ibz)**2
     end if
   end do

   ! q-points for optical limit.
   do iqlwl=1,Vcp%nqlwl
     call cmod_qpg(Vcp%nqlwl,iqlwl,Vcp%qlwl,ng,Gsph%gvec,gprimd,vcoul_lwl(:,iqlwl))
   end do
   vcoul_lwl = four_pi/vcoul_lwl**2
   !
   ! === Integration of 1/q^2 singularity ===
   ! * We use the auxiliary function from PRB 75, 205126 (2007) [[cite:Carrier2007]]
   q0_vol=(two_pi)**3/(Kmesh%nbz*ucvol) ; bz_geometry_factor=zero

   ! It might be useful to perform the integration using the routine quadrature
   ! so that we can control the accuracy of the integral and improve the
   ! numerical stability of the GW results.
   do iq_bz=1,Qmesh%nbz
     qbz_cart(:)=Qmesh%bz(1,iq_bz)*b1(:)+Qmesh%bz(2,iq_bz)*b2(:)+Qmesh%bz(3,iq_bz)*b3(:)
     qbz_norm=SQRT(SUM(qbz_cart(:)**2))
     if (qbz_norm>tolq0) bz_geometry_factor=bz_geometry_factor-faux(Qmesh%bz(:,iq_bz))
   end do

   bz_geometry_factor = bz_geometry_factor + integratefaux()*Qmesh%nbz

   write(msg,'(2a,2x,f8.4)')ch10,&
&    ' integrate q->0 : numerical BZ geometry factor = ',bz_geometry_factor*q0_vol**(2./3.)
   call wrtout(std_out,msg,'COLL')

   Vcp%i_sz=four_pi*bz_geometry_factor  ! Final result stored here

 CASE ('AUX_GB')
   !
   do iq_ibz=1,Vcp%nqibz
     call cmod_qpg(Vcp%nqibz,iq_ibz,Vcp%qibz,ng,Gsph%gvec,gprimd,vcoul(:,iq_ibz))

     if (iq_ibz==1) then ! The singularity is treated using vcoul_lwl.
       vcoul(1, iq_ibz) = zero
       vcoul(2:,iq_ibz) = four_pi/vcoul(2:,iq_ibz)**2
     else
       vcoul(:,iq_ibz) = four_pi/vcoul(:,iq_ibz)**2
     end if
   end do

   ! q-points for optical limit.
   do iqlwl=1,Vcp%nqlwl
     call cmod_qpg(Vcp%nqlwl,iqlwl,Vcp%qlwl,ng,Gsph%gvec,gprimd,vcoul_lwl(:,iqlwl))
   end do
   vcoul_lwl = four_pi/vcoul_lwl**2
   !
   ! === Integration of 1/q^2 singularity ===
   ! * We use the auxiliary function of a Gygi-Baldereschi variant [[cite:Gigy1986]]
   q0_vol=(two_pi)**3/(Kmesh%nbz*ucvol) ; bz_geometry_factor=zero
   ! * the choice of alfa (the width of gaussian) is somehow empirical
   alfa = 150.0/ecut
   !write(msg,'(2a,2x,f12.4,2x,f12.4)')ch10, ' alfa, ecutsigx = ', alfa, ecut
   !call wrtout(std_out,msg,'COLL')

   do iq_bz=1,Qmesh%nbz
     do ig = 1,ng
        qqgg(:) = Qmesh%bz(:,iq_bz) + Gsph%gvec(:,ig)
        qqgg(:) = qqgg(1)*b1(:) + qqgg(2)*b2(:) + qqgg(3)*b3(:)
        qbz_norm=SQRT(SUM(qqgg(:)**2))
     if (qbz_norm>tolq0*0.01) then
        bz_geometry_factor = bz_geometry_factor - EXP(-alfa*qbz_norm**2)/qbz_norm**2
     end if
     end do
   end do

   intfauxgb = ucvol/four_pi/SQRT(0.5*two_pi*alfa)
   !write(msg,'(2a,2x,f12.4,2x,f12.4)')ch10, ' bz_geometry_factor, f(q) integral = ', &
   !&    bz_geometry_factor, intfauxgb
   !call wrtout(std_out,msg,'COLL')

   bz_geometry_factor = bz_geometry_factor + intfauxgb*Qmesh%nbz

   write(msg,'(2a,2x,f8.4)')ch10,&
&    ' integrate q->0 : numerical BZ geometry factor = ',bz_geometry_factor*q0_vol**(2./3.)
   call wrtout(std_out,msg,'COLL')

   Vcp%i_sz=four_pi*bz_geometry_factor

 CASE ('CRYSTAL')

   do iq_ibz=1,Vcp%nqibz
     call cmod_qpg(Vcp%nqibz,iq_ibz,Vcp%qibz,ng,Gsph%gvec,gprimd,vcoul(:,iq_ibz))

     if (iq_ibz==1) then ! The singularity is treated using vcoul_lwl.
       vcoul(1, iq_ibz) = zero
       vcoul(2:,iq_ibz) = four_pi/vcoul(2:,iq_ibz)**2
     else
       vcoul(:,iq_ibz) = four_pi/vcoul(:,iq_ibz)**2
     end if
   end do

   ! q-points for optical limit.
   do iqlwl=1,Vcp%nqlwl
     call cmod_qpg(Vcp%nqlwl,iqlwl,Vcp%qlwl,ng,Gsph%gvec,gprimd,vcoul_lwl(:,iqlwl))
   end do
   vcoul_lwl = four_pi/vcoul_lwl**2
   !
   ! === Integration of 1/q^2 singularity ===
   ! * We use the auxiliary function from PRB 75, 205126 (2007) [[cite:Carrier2007]]
   q0_vol=(two_pi)**3/(Kmesh%nbz*ucvol) ; bz_geometry_factor=zero

   ! $$ MG: this is to restore the previous implementation, it will facilitate the merge
   ! Analytic integration of 4pi/q^2 over the volume element:
   ! $4pi/V \int_V d^3q 1/q^2 =4pi bz_geometric_factor V^(-2/3)$
   ! i_sz=4*pi*bz_geometry_factor*q0_vol**(-two_thirds) where q0_vol= V_BZ/N_k
   ! bz_geometry_factor: sphere=7.79, fcc=7.44, sc=6.188, bcc=6.946, wz=5.255
   ! (see gwa.pdf, appendix A.4)
   Vcp%i_sz=four_pi*7.44*q0_vol**(-two_thirds)

 CASE ('ERF')
   ! Modified long-range only Coulomb interaction thanks to the error function:
   ! * Vc = erf(r/rcut)/r
   ! * The singularity is treated using vcoul_lwl.
   do iq_ibz=1,Qmesh%nibz
     call cmod_qpg(Qmesh%nibz,iq_ibz,Qmesh%ibz,ng,Gsph%gvec,gprimd,vcoul(:,iq_ibz))

     ! * The Fourier transform of the error function reads
     if (iq_ibz==1) then
       vcoul(1, iq_ibz) = zero
       vcoul(2:,iq_ibz) = four_pi/(vcoul(2:,iq_ibz)**2) *  EXP( -0.25d0 * (Vcp%rcut*vcoul(2:,iq_ibz))**2 )
     else
       vcoul(:,iq_ibz)  = four_pi/(vcoul(:, iq_ibz)**2) *  EXP( -0.25d0 * (Vcp%rcut*vcoul(: ,iq_ibz))**2 )
     end if
   end do

   ! q-points for optical limit.
   do iqlwl=1,Vcp%nqlwl
     call cmod_qpg(Vcp%nqlwl,iqlwl,Vcp%qlwl,ng,Gsph%gvec,gprimd,vcoul_lwl(:,iqlwl))
   end do
   vcoul_lwl = four_pi/(vcoul_lwl**2) *  EXP( -0.25d0 * (Vcp%rcut*vcoul_lwl)**2 )

   !
   ! === Treat 1/q^2 singularity ===
   ! * We use the auxiliary function from PRB 75, 205126 (2007) [[cite:Carrier2007]]
   q0_vol=(two_pi)**3/(Kmesh%nbz*ucvol) ; bz_geometry_factor=zero
   do iq_bz=1,Qmesh%nbz
     qbz_cart(:)=Qmesh%bz(1,iq_bz)*b1(:)+Qmesh%bz(2,iq_bz)*b2(:)+Qmesh%bz(3,iq_bz)*b3(:)
     qbz_norm=SQRT(SUM(qbz_cart(:)**2))
     if (qbz_norm>tolq0) bz_geometry_factor=bz_geometry_factor-faux(Qmesh%bz(:,iq_bz))
   end do

   bz_geometry_factor = bz_geometry_factor + integratefaux()*Qmesh%nbz

   write(msg,'(2a,2x,f12.4)')ch10,&
&    ' integrate q->0 : numerical BZ geometry factor = ',bz_geometry_factor*q0_vol**(2./3.)
   call wrtout(std_out,msg,'COLL')
   Vcp%i_sz=four_pi*bz_geometry_factor  ! Final result stored here

 CASE ('ERFC')
   ! * Use a modified short-range only Coulomb interaction thanks to the complementar error function:
   !   $ V_c = [1-erf(r/r_{cut})]/r $
   ! * The Fourier transform of the error function reads
   !   vcoul=four_pi/(vcoul**2) * ( 1.d0 - exp( -0.25d0 * (Vcp%rcut*vcoul)**2 ) )
   do iq_ibz=1,Qmesh%nibz
     call cmod_qpg(Qmesh%nibz,iq_ibz,Qmesh%ibz,ng,Gsph%gvec,gprimd,vcoul(:,iq_ibz))

     if (iq_ibz==1) then
       vcoul(1 ,iq_ibz) = zero
       vcoul(2:,iq_ibz) = four_pi/(vcoul(2:,iq_ibz)**2) * ( one - EXP( -0.25d0 * (Vcp%rcut*vcoul(2:,iq_ibz))**2 ) )
     else
       vcoul(:, iq_ibz) = four_pi/(vcoul(:, iq_ibz)**2) * ( one - EXP( -0.25d0 * (Vcp%rcut*vcoul(:, iq_ibz))**2 ) )
     end if
   end do
   !
   ! q-points for optical limit.
   do iqlwl=1,Vcp%nqlwl
     call cmod_qpg(Vcp%nqlwl,iqlwl,Vcp%qlwl,ng,Gsph%gvec,gprimd,vcoul_lwl(:,iqlwl))
   end do
   vcoul_lwl = four_pi/(vcoul_lwl**2) * ( one - EXP( -0.25d0 * (Vcp%rcut*vcoul_lwl)**2 ) )
   !
   ! === Treat 1/q^2 singularity ===
   ! * There is NO singularity in this case.
   Vcp%i_sz=pi*Vcp%rcut**2 ! Final result stored here

 CASE DEFAULT
   MSG_BUG(sjoin('Unknown cutoff mode:',Vcp%mode))
 END SELECT
 !
 ! === Store final results in complex array ===
 ! * Rozzi"s cutoff can give real negative values

 ABI_MALLOC(Vcp%vc_sqrt,(ng,Vcp%nqibz))
 ABI_MALLOC(Vcp%vc_sqrt_resid,(ng,Vcp%nqibz))
 Vcp%vc_sqrt=CMPLX(vcoul,zero)
 Vcp%vc_sqrt=SQRT(Vcp%vc_sqrt)
 Vcp%vc_sqrt_resid=Vcp%vc_sqrt
 Vcp%i_sz_resid=Vcp%i_sz
 call vcoul_plot(Vcp,Qmesh,Gsph,ng,vcoul,comm)
 ABI_FREE(vcoul)

 ABI_MALLOC(Vcp%vcqlwl_sqrt,(ng,Vcp%nqlwl))
 Vcp%vcqlwl_sqrt=CMPLX(vcoul_lwl,zero)
 Vcp%vcqlwl_sqrt=SQRT(Vcp%vcqlwl_sqrt)
 ABI_FREE(vcoul_lwl)

 call vcoul_print(Vcp,unit=std_out)

 DBG_EXIT("COLL")

contains !===============================================================

 real(dp) function integratefaux()

!Local variables-------------------------------
  integer,parameter :: nref=3
  integer,parameter :: nq=50
  integer :: ierr,iq,iqx1,iqy1,iqz1,iqx2,iqy2,iqz2,miniqy1,maxiqy1,nqhalf
  real(dp) :: invnq,invnq3,qq,weightq,weightxy,weightxyz
  real(dp) :: qq1(3),qq2(3),bb4sinpiqq_2(3,nq),sin2piqq(nq),bb4sinpiqq2_2(3,0:nq),sin2piqq2(3,0:nq)

  ! nq is the number of sampling points along each axis for the numerical integration
  ! nref is the area where the mesh is refined

  integratefaux=zero
  invnq=one/DBLE(nq)
  invnq3=invnq**3
  nqhalf=nq/2

  !In order to speed up the calculation, precompute the sines
  do iq=1,nq
    qq=DBLE(iq)*invnq-half
    bb4sinpiqq_2(:,iq)=bb(:)*four*SIN(pi*qq)**2 ; sin2piqq(iq)=SIN(two_pi*qq)
  end do

  do iqx1=1,nq
    if(modulo(iqx1,nprocs)/=rank) cycle
    qq1(1)=DBLE(iqx1)*invnq-half
    !Here take advantage of the q <=> -q symmetry : arrange the sampling of qx, qy space to avoid duplicating calculations.
    !Need weights to do this ...
    !do iqy1=1,nq
    miniqy1=nqhalf+1 ; maxiqy1=nq
    if(iqx1>=nqhalf)miniqy1=nqhalf
    if(iqx1>nqhalf .and. iqx1<nq)maxiqy1=nq-1
    do iqy1=miniqy1,maxiqy1
      qq1(2)=DBLE(iqy1)*invnq-half
      !By default, the factor of two is for the q <=> -q symmetry
      weightq=invnq3*two
      !But not all qx qy lines have a symmetric one ...
      if( (iqx1==nqhalf .or. iqx1==nq) .and. (iqy1==nqhalf .or. iqy1==nq))weightq=weightq*half

      do iqz1=1,nq
        qq1(3)=DBLE(iqz1)*invnq-half
        !
        ! Refine the mesh for the point close to the origin
        if( abs(iqx1-nqhalf)<=nref .and. abs(iqy1-nqhalf)<=nref .and. abs(iqz1-nqhalf)<=nref ) then
          !Note that the set of point is symmetric around the central point, while weights are taken into account
          do iq=0,nq
            qq2(:)=qq1(:)+ (DBLE(iq)*invnq-half)*invnq
            bb4sinpiqq2_2(:,iq)=bb(:)*four*SIN(pi*qq2(:))**2 ;  sin2piqq2(:,iq)=SIN(two_pi*qq2(:))
          enddo
          do iqx2=0,nq
            qq2(1)=qq1(1) + (DBLE(iqx2)*invnq-half ) *invnq
            do iqy2=0,nq
              qq2(2)=qq1(2) + (DBLE(iqy2)*invnq-half ) *invnq
              weightxy=invnq3*weightq
              if(iqx2==0 .or. iqx2==nq)weightxy=weightxy*half
              if(iqy2==0 .or. iqy2==nq)weightxy=weightxy*half
              do iqz2=0,nq
                qq2(3)=qq1(3) + (DBLE(iqz2)*invnq-half ) *invnq
                weightxyz=weightxy
                if(iqz2==0 .or. iqz2==nq)weightxyz=weightxy*half
                !
                ! Treat the remaining divergence in the origin as if it would be a spherical
                ! integration of 1/q^2
                if( iqx1/=nqhalf .or. iqy1/=nqhalf .or. iqz1/=nqhalf .or. iqx2/=nqhalf .or. iqy2/=nqhalf .or. iqz2/=nqhalf ) then
!                 integratefaux=integratefaux+ faux(qq2) *invnq**6
                  integratefaux=integratefaux+ faux_fast(qq2,bb4sinpiqq2_2(1,iqx2),bb4sinpiqq2_2(2,iqy2),bb4sinpiqq2_2(3,iqz2),&
&                                                        sin2piqq2(1,iqx2),sin2piqq2(2,iqy2),sin2piqq2(3,iqz2)) * weightxyz
                else
                   integratefaux=integratefaux + 7.7955* ( (two_pi)**3/ucvol*invnq3*invnq3 )**(-2./3.) *invnq3*invnq3
                end if
              end do
            end do
          end do
        else
!        integratefaux=integratefaux+faux(qq1)*invnq**3
         integratefaux=integratefaux+ faux_fast(qq1,bb4sinpiqq_2(1,iqx1),bb4sinpiqq_2(2,iqy1),bb4sinpiqq_2(3,iqz1),&
&                                               sin2piqq(iqx1),sin2piqq(iqy1),sin2piqq(iqz1)) * weightq
        end if
      end do
    end do
  end do

  call xmpi_sum(integratefaux,comm,ierr)

 end function integratefaux

 function faux(qq)

  real(dp),intent(in) :: qq(3)
  real(dp) :: faux
  real(dp) :: bb4sinpiqq1_2,bb4sinpiqq2_2,bb4sinpiqq3_2,sin2piqq1,sin2piqq2,sin2piqq3

  bb4sinpiqq1_2=b1b1*four*SIN(pi*qq(1))**2
  bb4sinpiqq2_2=b2b2*four*SIN(pi*qq(2))**2
  bb4sinpiqq3_2=b3b3*four*SIN(pi*qq(3))**2
  sin2piqq1=SIN(two_pi*qq(1))
  sin2piqq2=SIN(two_pi*qq(2))
  sin2piqq3=SIN(two_pi*qq(3))

  faux=faux_fast(qq,bb4sinpiqq1_2,bb4sinpiqq2_2,bb4sinpiqq3_2,sin2piqq1,sin2piqq2,sin2piqq3)

 end function faux

 function faux_fast(qq,bb4sinpiqq1_2,bb4sinpiqq2_2,bb4sinpiqq3_2,sin2piqq1,sin2piqq2,sin2piqq3)

  real(dp),intent(in) :: qq(3)
  real(dp) :: faux_fast
  real(dp) :: bb4sinpiqq1_2,bb4sinpiqq2_2,bb4sinpiqq3_2,sin2piqq1,sin2piqq2,sin2piqq3

  faux_fast= bb4sinpiqq1_2 + bb4sinpiqq2_2 + bb4sinpiqq3_2 &
&       +two*( b1b2 * sin2piqq1*sin2piqq2 &
&             +b2b3 * sin2piqq2*sin2piqq3 &
&             +b3b1 * sin2piqq3*sin2piqq1 &
&            )
  if(Vcp%rcut>tol6)then
    faux_fast=two_pi*two_pi/faux_fast * exp( -0.25d0*Vcp%rcut**2* sum( ( qq(1)*b1(:)+qq(2)*b2(:)+qq(3)*b3(:) )**2 ) )
  else
    faux_fast=two_pi*two_pi/faux_fast
  endif

 end function faux_fast


 function adapt_nmc(nmc_max,qpg2) result(nmc)

 real(dp),intent(in) :: qpg2
 integer,intent(in)  :: nmc_max
 integer :: nmc

 ! Empirical law to decrease the Monte Carlo sampling
 ! for large q+G, for which the accuracy is not an issue
 nmc = NINT( nmc_max / ( 1.0_dp + 1.0_dp * qpg2**6 ) )
 nmc = MIN(nmc_max,nmc)
 nmc = MAX(1,nmc)

 end function adapt_nmc


end subroutine vcoul_init
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/vcoul_plot
!! NAME
!! vcoul_plot
!!
!! FUNCTION
!! Plot vccut(q,G) as a function of |q+G|. Calculate also vc in real space.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_vcoul
!!
!! CHILDREN
!!      calck0,paw_jbessel,quadrature
!!
!! SOURCE

subroutine vcoul_plot(Vcp,Qmesh,Gsph,ng,vc,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng,comm
 type(kmesh_t),intent(in) :: Qmesh
 type(vcoul_t),intent(in) :: Vcp
 type(gsphere_t),intent(in) :: Gsph
!arrays
 real(dp),intent(in) :: vc(ng,Qmesh%nibz)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: icount,idx_Sm1G,ierr,ig,igs,ii,iq_bz,iq_ibz,iqg,ir,isym,itim
 integer :: my_start,my_stop,nqbz,nqibz,nr,ntasks,rank,unt
 real(dp) :: arg,fact,l1,l2,l3,lmax,step,tmp,vcft,vc_bare
 character(len=500) :: msg
 character(len=fnlen) :: filnam
!arrays
 integer,allocatable :: insort(:)
 real(dp) :: b1(3),b2(3),b3(3),gmet(3,3),gprimd(3,3),qbz(3),qpgc(3)
 real(dp),allocatable :: qpg_mod(:),rr(:,:,:),vcr(:,:),vcr_cut(:,:)

!************************************************************************

 if (TRIM(Vcp%mode)/='CYLINDER') RETURN

 rank = xmpi_comm_rank(comm)

 nqibz=Qmesh%nibz; nqbz=Qmesh%nbz
 gmet=Gsph%gmet; gprimd=Gsph%gprimd

 b1(:)=two_pi*gprimd(:,1)
 b2(:)=two_pi*gprimd(:,2)
 b3(:)=two_pi*gprimd(:,3)

 ! === Compare in Fourier space the true Coulomb with the cutted one ===
 if (rank==master) then
   ABI_MALLOC(insort,(nqibz*ng))
   ABI_MALLOC(qpg_mod,(nqibz*ng))
   iqg=1
   do iq_ibz=1,nqibz
     do ig=1,ng
       qpg_mod(iqg)=normv(Qmesh%ibz(:,iq_ibz)+Gsph%gvec(:,ig),gmet,'g')
       insort(iqg)=iqg; iqg=iqg+1
     end do
   end do
   call sort_dp(nqibz*ng,qpg_mod,insort,tol14)

   filnam='_VCoulFT_'
   call isfile(filnam,'new')
   if (open_file(filnam,msg,newunit=unt,status='new',form='formatted') /=0) then
     MSG_ERROR(msg)
   end if
   write(unt,'(a,i3,a,i6,a)')&
&    '#   |q+G|       q-point (Tot no.',nqibz,')        Gvec (',ng,')     vc_bare(q,G)    vc_cutoff(q,G) '

   do iqg=1,nqibz*ng
     iq_ibz=(insort(iqg)-1)/ng +1
     ig=(insort(iqg))-(iq_ibz-1)*ng
     vc_bare=zero
     if (qpg_mod(iqg)>tol16) vc_bare=four_pi/qpg_mod(iqg)**2
     write(unt,'(f12.6,2x,3f8.4,2x,3i6,2x,2es14.6)')&
&      qpg_mod(iqg),Qmesh%ibz(:,iq_ibz),Gsph%gvec(:,ig),vc_bare,vc(ig,iq_ibz)
   end do

   close(unt)
   ABI_FREE(insort)
   ABI_FREE(qpg_mod)
 end if ! rank==master

 ! === Fourier transform back to real space just to check cutoff implementation ===
 ntasks=nqbz*ng
 call xmpi_split_work(ntasks,comm,my_start,my_stop)

 l1=SQRT(SUM(Vcp%rprimd(:,1)**2))
 l2=SQRT(SUM(Vcp%rprimd(:,2)**2))
 l3=SQRT(SUM(Vcp%rprimd(:,3)**2))

 nr=50
 lmax=MAX(l1,l2,l3) ; step=lmax/(nr-1)
 fact=one/(Vcp%ucvol*nqbz)
 !
 ! numb coding
 ABI_MALLOC(rr,(3,nr,3))
 rr=zero
 do ii=1,3
   do ir=1,nr
     rr(ii,ir,ii)=(ir-1)*step
   end do
 end do

 ABI_MALLOC(vcr,(nr,3))
 ABI_MALLOC(vcr_cut,(nr,3))
 vcr=zero; vcr_cut=zero

 do iq_bz=1,nqbz
   call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym,itim)
   if (ABS(qbz(1))<0.01) qbz(1)=zero
   if (ABS(qbz(2))<0.01) qbz(2)=zero
   if (ABS(qbz(3))<0.01) qbz(3)=zero
   igs=1 ; if (ALL(qbz(:)==zero)) igs=2
   do ig=igs,ng
     icount=ig+(iq_bz-1)*ng
     if (icount<my_start.or.icount>my_stop) CYCLE
     idx_Sm1G=Gsph%rottbm1(ig,itim,isym) ! IS{^-1}G
     vcft=vc(idx_Sm1G,iq_ibz)
     qpgc(:)=qbz(:)+Gsph%gvec(:,ig) ; qpgc(:)=b1(:)*qpgc(1)+b2(:)*qpgc(2)+b3(:)*qpgc(3)
     tmp=SQRT(DOT_PRODUCT(qpgc,qpgc)) ; tmp=tmp**2
     do ii=1,3
       do ir=1,nr
         arg=DOT_PRODUCT(rr(:,ir,ii),qpgc)
         vcr_cut(ir,ii)=vcr_cut(ir,ii) + vcft*COS(arg)
         vcr    (ir,ii)=vcr    (ir,ii) + four_pi/tmp*COS(arg)
       end do
     end do
   end do !ig
 end do !iq_ibz

 call xmpi_sum_master(vcr_cut,master,comm,ierr)
 call xmpi_sum_master(vcr    ,master,comm,ierr)

 if (rank==master) then
   filnam='_VCoulR_'
   call isfile(filnam,'new')
   if (open_file(filnam,msg,newunit=unt,status='new',form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if
   write(unt,'(a)')'# length  vc_bare(r)   vc_cut(r) '
   do ir=1,nr
     write(unt,'(7es18.6)')(ir-1)*step,(fact*vcr(ir,ii),fact*vcr_cut(ir,ii),ii=1,3)
   end do
   close(unt)
 end if

 ABI_FREE(rr)
 ABI_FREE(vcr)
 ABI_FREE(vcr_cut)

end subroutine vcoul_plot
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/vcoul_print
!! NAME
!! vcoul_print
!!
!! FUNCTION
!!  Print the content of a Coulomb datatype.
!!
!! INPUTS
!!  Vcp<vcoul_t>=The datatype whose content has to be printed.
!!  [unit]=The unit number for output
!!  [prtvol]=Verbosity level
!!  [mode_paral]=Either "COLL" or "PERS".
!!
!! PARENTS
!!      m_vcoul
!!
!! CHILDREN
!!      calck0,paw_jbessel,quadrature
!!
!! SOURCE

subroutine vcoul_print(Vcp,unit,prtvol,mode_paral)

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral
 type(vcoul_t),intent(in) :: Vcp

!Local variables-------------------------------
!scalars
 integer :: ii,my_unt,my_prtvol,iqlwl
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt    =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_mode   ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral
 my_prtvol=0       ; if (PRESENT(prtvol    )) my_prtvol=prtvol

 SELECT CASE (Vcp%mode)

 CASE ('MINIBZ')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',TRIM(Vcp%mode)
   call wrtout(my_unt,msg,my_mode)

 CASE ('MINIBZ-ERF')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',TRIM(Vcp%mode)
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(5a,f10.4,3a,f10.2,3a,3f10.5,2a)')ch10,&
&    ' === Error function cutoff === ',ch10,ch10,&
&    '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10
   call wrtout(my_unt,msg,my_mode)

 CASE ('MINIBZ-ERFC')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',TRIM(Vcp%mode)
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(5a,f10.4,3a,f10.2,3a,3f10.5,2a)')ch10,&
&    ' === Complement Error function cutoff === ',ch10,ch10,&
&    '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10
   call wrtout(my_unt,msg,my_mode)

 CASE ('SPHERE')
   write(msg,'(5a,f10.4,3a,f10.2,3a,3f10.5,2a)')ch10,&
&    ' === Spherical cutoff === ',ch10,ch10,&
&    '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10,&
&    '  Volume of the sphere .. ',four_pi/three*Vcp%rcut**3,' [Bohr^3] '
!FB: This has no meaning here! &  '  Sphere centered at .... ',Vcp%boxcenter,' (r.l.u) ',ch10
!MG It might be useful if the system is not centered on the origin because in this case the
!   matrix elements of the Coulomb have to be multiplied by a phase depending on boxcenter.
!   I still have to decide if it is useful to code this possibility and which variable use to
!   define the center (boxcenter is used in the tddft part).
   call wrtout(my_unt,msg,my_mode)

 CASE ('CYLINDER')
   ii=imin_loc(ABS(Vcp%pdir-1))
   write(msg,'(5a,f10.4,3a,i2,2a,3f10.2,a)')ch10,&
&    ' === Cylindrical cutoff === ',ch10,ch10,&
&    '  Cutoff radius ............... ',Vcp%rcut,' [Bohr] ',ch10,&
&    '  Axis parallel to direction... ',ii,ch10,&
&    '  Passing through point ....... ',Vcp%boxcenter,' (r.l.u) '
   call wrtout(my_unt,msg,my_mode)

   write(msg,'(2a)')'  Infinite length  ....... ',ch10
   if (Vcp%hcyl/=zero) write(msg,'(a,f8.5,2a)')'  Finite length of ....... ',Vcp%hcyl,' [Bohr] ',ch10
   call wrtout(my_unt,msg,my_mode)

 CASE ('SURFACE')
   write(msg,'(5a,f10.4,3a,3f10.2,2a)')ch10,&
&    ' === Surface cutoff === ',ch10,ch10,&
&    '  Cutoff radius .................... ',Vcp%rcut,' [Bohr] ',ch10,&
&    '  Central plane passing through .... ',Vcp%boxcenter,' (r.l.u) ',ch10
   call wrtout(my_unt,msg,my_mode)
   !write(msg,'(a)')'  Infinite length  .......'
   !if (Vcp%hcyl/=zero) write(msg,'(a,f8.5,a)')'  Finite length of .......',Vcp%hcyl,' [Bohr] '
   !call wrtout(my_unt,msg,my_mode)

 CASE ('AUXILIARY_FUNCTION')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',TRIM(Vcp%mode)
   call wrtout(my_unt,msg,my_mode)

 CASE ('AUX_GB')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',TRIM(Vcp%mode)
   call wrtout(my_unt,msg,my_mode)

 CASE ('CRYSTAL')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',TRIM(Vcp%mode)
   call wrtout(my_unt,msg,my_mode)

 CASE ('ERF')
   write(msg,'(5a,f10.4,3a,f10.2,3a,3f10.5,2a)')ch10,&
&    ' === Error function cutoff === ',ch10,ch10,&
&    '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10
   call wrtout(my_unt,msg,my_mode)

 CASE ('ERFC')
   write(msg,'(5a,f10.4,3a,f10.2,3a,3f10.5,2a)')ch10,&
&    ' === Complement Error function cutoff === ',ch10,ch10,&
&    '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10
   call wrtout(my_unt,msg,my_mode)

 CASE DEFAULT
   MSG_BUG(sjoin('Unknown cutoff mode: ', Vcp%mode))
 END SELECT

 if (Vcp%nqlwl>0) then
   write(msg,'(a,i3)')" q-points for optical limit: ",Vcp%nqlwl
   call wrtout(my_unt,msg,my_mode)
   do iqlwl=1,Vcp%nqlwl
     write(msg,'(1x,i5,a,2x,3f12.6)') iqlwl,')',Vcp%qlwl(:,iqlwl)
     call wrtout(my_unt,msg,my_mode)
   end do
 end if

 !TODO add additional information

end subroutine vcoul_print
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/vcoul_free
!! NAME
!! vcoul_free
!!
!! FUNCTION
!!  Destroy a vcoul_t type
!!
!! SIDE EFFECTS
!!  Vcp<vcoul_t>=the datatype to be destroyed
!!
!! PARENTS
!!      bethe_salpeter,gwls_hamiltonian,mrgscr,screening,setup_sigma,sigma
!!
!! CHILDREN
!!      calck0,paw_jbessel,quadrature
!!
!! SOURCE

subroutine vcoul_free(Vcp)

!Arguments ------------------------------------
!scalars
 type(vcoul_t),intent(inout) :: Vcp

! *************************************************************************

 ABI_SFREE(Vcp%qibz)
 ABI_SFREE(Vcp%qlwl)
 ABI_SFREE(Vcp%vc_sqrt)
 ABI_SFREE(Vcp%vc_sqrt_resid)
 ABI_SFREE(Vcp%vcqlwl_sqrt)

end subroutine vcoul_free
!!***

!----------------------------------------------------------------------

function K0cos_dy_r0(xx)

 real(dp),intent(in) :: xx
 real(dp) :: K0cos_dy_r0

!Local variables-------------------------------
!scalars
 integer :: ierr
 real(dp) :: quad,yx
!************************************************************************

 ! $ K0cos_dy_r0(x)= \int_{-b/2}^{-y(x)} K0(|qpg_z|\rho) cos(x.qpg_x+y.qpg_y)dy
 !                  +\int_{y(x)}^{b/2} K0(|qpg_z|\rho)cos(x.qpg_x+y.qpg_y)dy$
 ! where y(x)=SQRT(r0^2-x^2) and x<=r0
 !
 xx_=xx; yx=SQRT(r0_**2-xx**2)
 call quadrature(K0cos,-hb_,-yx,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
 if (ierr/=0) then
   MSG_ERROR("Accuracy not reached")
 end if
 K0cos_dy_r0=quad

 call quadrature(K0cos,+yx,+hb_,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
 if (ierr/=0) then
   MSG_ERROR("Accuracy not reached")
 end if

 K0cos_dy_r0=quad+K0cos_dy_r0

end function K0cos_dy_r0
!!***

!----------------------------------------------------------------------

function K0cos_dth_r0(rho)

 real(dp),intent(in) :: rho
 real(dp) :: K0cos_dth_r0

!Local variables-------------------------------
!scalars
 integer :: ierr
 real(dp) :: quad,arg,k0,tmp
!************************************************************************

 ! $ K0cos_dth_r0(\rho)=
 ! \int_{0}^{2pi)} K0(|qpg_z|\rho)cos(\rho.cos(\theta).qpg_x+\rho.sin(\theta).qpg_y) d\theta $
 !
 ! where y(x)=SQRT(r0^2-x^2) and x<=r0
 !
 rho_=rho
 call quadrature(Fcos_th,zero,two_pi,qopt_,quad,ierr,ntrial_,accuracy_,npts_)
 if (ierr/=0) then
   MSG_ERROR("Accuracy not reached")
 end if

 arg=qpg_para_*rho_
 tmp=zero
 if (arg>tol6) then
   call CALCK0(arg,k0,1)
   tmp=k0*rho_
 end if
 K0cos_dth_r0=quad*tmp

end function K0cos_dth_r0
!!***

!----------------------------------------------------------------------

pure function Fcos_th(theta)

 real(dp),intent(in) :: theta
 real(dp) :: Fcos_th

!************************************************************************

 ! $ Fcos_th(\theta)=rho*K0(\rho*|qpg_z|)*COS(\rho.COS(\theta).qpg_x+\rho.SIN/(\theta)*qpg_y) $

 !arg=qpg_para_*rho_
 !call CALCK0(arg,k0,1)
 !tmp=k0*rho_
 Fcos_th=COS(rho_*COS(theta)*qpgx_+rho_*SIN(theta)*qpgy_)

end function Fcos_th
!!***

!----------------------------------------------------------------------

!the following functions should be used to deal with the singularity in the Cylindrical cutoff
!TODO Not yet used and indeed are still private

function K0fit(mq,nn) result(vals)

 integer,intent(in) :: nn
 real(dp),intent(in) :: mq
 real(dp) :: vals(nn)
!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: mqh
!arrays
 real(dp),parameter :: cc(7)=(/-0.57721566,0.42278420,0.23069756,&
&                               0.03488590,0.00262698,0.00010750,0.00000740/)
 ! *************************************************************************

 if (nn>8.or.nn<1) then
   MSG_ERROR("nn>8.or.nn<1 not implemented")
 end if

 ! === Eq 9.8.5 in Abramovitz ===
 vals(1)=-LOG(mq*half)*I0(mq)
 mqh=mq*half
 do ii=2,nn
   vals(ii)=cc(ii-1)*mqh**(2*(ii-2))
 end do

end function K0fit

function K0fit_int(mq,par,nn) result(integ)

 integer,intent(in) :: nn
 real(dp),intent(in) :: mq
 real(dp) :: integ
 real(dp),intent(in) :: par(nn)

!Local variables-------------------------------
!scalars
 integer :: ii,aa
 real(dp) :: mqh
!arrays
 real(dp),parameter :: cc(7)=(/-0.57721566,0.42278420,0.23069756,&
&                               0.03488590,0.00262698,0.00010750,0.00000740/)
 ! *************************************************************************

 if (nn>8.or.nn<1) then
   MSG_ERROR("nn>8.or.nn<1 not implemented")
 end if

 mqh=mq*half
 integ=-par(1)*int_I0ln(mqh)
 ! primitive of polynomial \sum_0^{N/2} cc_{2i} (x/2)^{2*i}
 do ii=2,nn
  aa=(2*(ii-1)+1)
  integ=integ+par(ii)*two*cc(ii-1)*(mqh**aa)/aa
 end do

end function K0fit_int

function I0(xx)

 real(dp),intent(in) :: xx
 real(dp) :: I0

!Local variables-------------------------------
 real(dp) :: tt

! *************************************************************************

 ! Eq 9.8.1 of Abramovitz, entering the expansion of K0 -->0
 ! Expansion holds for |x|<3.75, Error<1.6*10D-07
 tt=xx/3.75
 I0=one+3.5156229*tt**2+3.0899424*tt**4 +1.2067492*tt**6 &
       +0.2659732*tt**8+0.0360768*tt**10+0.0045813*tt**12
end function I0

! Primitive of x^m Ln(x) for m/=-1
function int_xmln(xx,mm)  result(res)

 integer,intent(in) :: mm
 real(dp),intent(in) :: xx
 real(dp) :: res

! *********************************************************************

 if (mm==-1) then
   MSG_BUG('invalid value for mm')
 end if

 if (xx<=zero) then
   MSG_BUG(' invalid value for xx')
 end if

 res= (xx**(mm+1))/(mm+1) * (LOG(xx) - one/(mm+1))

end function int_xmln

! Primitive function of ln(x/2)*I0(x) = sum_0^{N/2} 2^{2s+1} c_{2s} T(x/2,2s)
! where T(x,s)=\int x^s ln(x)dx
function int_I0ln(xx) result(res)

!Arguments ------------------------------------
 real(dp),intent(in) :: xx
 real(dp) :: res
!Local variables-------------------------------
 real(dp) :: yy
! *********************************************************************

 yy=xx*half
 res =  (       one*2    *int_xmln(yy,0)  &
&        +3.5156229*2**3 *int_xmln(yy,2)  &
&        +3.0899424*2**5 *int_xmln(yy,4)  &
&        +1.2067492*2**7 *int_xmln(yy,6)  &
&        +0.2659732*2**9 *int_xmln(yy,8)  &
&        +0.0360768*2**11*int_xmln(yy,10) &
&        +0.0045813*2**13*int_xmln(yy,12) &
&       )

end function int_I0ln

END MODULE m_vcoul
!!***
