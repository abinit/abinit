!!****m* ABINIT/m_vcoul
!! NAME
!!  m_vcoul
!!
!! FUNCTION
!!  This module contains the definition of the vcoul_t as well
!!  as procedures to calculate the Coulomb interaction in reciprocal space
!!  taking into account a possible cutoff in real space.
!!  Procedures to deal with the singularity for q --> 0 are also provided.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2022 ABINIT group (MG, FB)
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

module m_vcoul

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_splines
 use m_sort

 use m_fstrings,        only : sjoin, itoa
 use m_special_funcs,   only : abi_derf
 use m_bessel,          only : calck0
 use m_io_tools,        only : open_file
 use m_gwdefs,          only : GW_TOLQ0
 use m_numeric_tools,   only : arth, geop, imin_loc, llsfit_svd, l2norm, OPERATOR(.x.), quadrature, isdiagmat
 use m_hide_lapack,     only : matrginv
 use m_geometry,        only : normv, metric
 use m_qplusg,          only : cmod_qpg
 use m_crystal,         only : crystal_t
 use m_bz_mesh,         only : kmesh_t
 use m_gsphere,         only : gsphere_t
 use m_fftcore,         only : get_kg
 use m_dtfil,           only : isfile

 ! Cut-off methods modules
 use m_cutoff_sphere,   only : cutoff_sphere
 use m_cutoff_surface,  only : cutoff_surface
 use m_cutoff_cylinder, only : cutoff_cylinder

 implicit none

 public :: gw_icutcoul_to_mode
 public :: carrier_isz

 private
!!***

!!****t* m_vcoul/vcoul_t
!! NAME
!!  vcoul_t
!!
!! FUNCTION
!!  This data type contains the square root of the Fourier components of the Coulomb interaction
!!  calculated taking into account a possible cutoff. It also stores info on the particular geometry
!!  used for the cutoff as well as quantities required to deal with the Coulomb divergence for q --> 0.
!!
!! SOURCE

type,public :: vcoul_t

  integer :: ng = -1
   ! Number of G-vectors

  integer :: nqibz = -1
   ! Number of irreducible q-points

  integer :: nqlwl = -1
   ! Number of small q-points around Gamma

  real(dp) :: alpha(3) = -one
   ! Lenght of the finite surface.

  real(dp) :: rcut = -one
   ! Cutoff radius.

  real(dp) :: i_sz = huge(one)
   ! Value of the integration of the Coulomb singularity 4\pi/V_BZ \int_BZ d^3q 1/q^2

  real(dp) :: i_sz_resid = huge(one)
   ! Residual difference between the i_sz in the sigma self-energy for exchange,
   ! and the i_sz already present in the generalized Kohn-Sham eigenenergies
   ! Initialized to the same value as i_sz

  real(dp) :: hcyl = -one
   ! Length of the finite cylinder along the periodic dimension

  real(dp) :: ucvol = -one
    ! Volume of the unit cell

  character(len=50) :: mode
   ! String defining the cutoff mode

  integer :: pdir(3)
   ! 1 if the system is periodic along this direction

  real(dp) :: boxcenter(3) = -1
   ! 1 if the point in inside the cutoff region 0 otherwise
   ! Reduced coordinates of the center of the box (input variable)

  real(dp) :: vcutgeo(3) = huge(one)
   ! For each reduced direction gives the length of the finite system
   ! 0 if the system is infinite along this direction.
   ! negative values indicate that a finite size has to be used.

  real(dp) :: rprimd(3,3) = zero
   ! Lattice vectors in real space.

  real(dp),allocatable :: qibz(:,:)
   ! (3, nqibz)
   ! q-points in the IBZ.

  real(dp),allocatable :: qlwl(:,:)
   ! (3, nqlwl)
   ! q-points for the treatment of the Coulomb singularity.

  complex(gwpc),allocatable :: vc_sqrt(:,:)
   ! (ng, nqibz)
   ! Square root of the Coulomb interaction in reciprocal space.
   ! complex-valued to allow for a possible cutoff (Rozzi's method)

  complex(gwpc),allocatable :: vcqlwl_sqrt(:,:)
   ! (ng, nqlwl)
   ! Square root of the Coulomb term calculated for small q-points

  complex(gwpc),allocatable :: vc_sqrt_resid(:,:)
   ! (ng, nqibz)
   ! Square root of the residual difference between the Coulomb interaction in the sigma self-energy for exchange,
   ! and the Coulomb interaction already present in the generalized Kohn-Sham eigenenergies (when they come from an hybrid)
   ! Given in reciprocal space. At the call to vcoul_init, it is simply initialized at the value of vc_sqrt(:,:),
   ! and only later modified. A cutoff might be applied.

contains
   procedure :: init => vcoul_init    ! Main creation method.
   procedure :: plot => vcoul_plot    ! Plot vc in real and reciprocal space.
   procedure :: print => vcoul_print  ! Report info on the object.
   procedure :: free => vcoul_free    ! Free memory

end type vcoul_t
!!***

!!****t* m_vcoul/mc_t
!! NAME
!!  mc_t
!!
!! FUNCTION
!! Mimicking the BerkeleyGW technique
!! A Monte-Carlo sampling of each miniBZ surrounding each (q+G) point
!! However:
!!    - extended to multiple shifts
!!    - with an adaptative number of MonteCarlo sampling points
!!
!! SOURCE

type, public :: mc_t

  integer :: nmc_max = -1

  real(dp) :: q0sph = -one

  real(dp) :: ucvol = -one

  real(dp) :: gmet(3,3) = -one

  real(dp),allocatable :: qran(:,:)
  ! (3, nmc_max)

contains
  procedure :: init => mc_init
  procedure :: integrate => mc_integrate
  procedure :: free => mc_free
end type mc_t
!!***

!!****t* m_vcoul/vcgen_t
!! NAME
!!  vcgen_t
!!
!! FUNCTION
!!
!! SOURCE

type, public :: vcgen_t

  integer :: nkbz = -1
   ! Number of k-points in full BZ.

  integer :: opt_cylinder

  integer :: opt_surface

  real(dp) :: alpha(3) = -one
   ! Lenght of the finite surface.

  real(dp) :: rcut = -one
   ! Cutoff radius.

  real(dp) :: hcyl = -one
   ! Length of the finite cylinder along the periodic dimension

  character(len=50) :: mode
   ! String defining the cutoff mode

  integer :: pdir(3)
   ! 1 if the system is periodic along this direction

  real(dp) :: boxcenter(3) = -1
   ! 1 if the point in inside the cutoff region 0 otherwise
   ! Reduced coordinates of the center of the box (input variable)

  real(dp) :: vcutgeo(3) = huge(one)
   ! For each reduced direction gives the length of the finite system
   ! 0 if the system is infinite along this direction.
   ! negative values indicate that a finite size has to be used.

  real(dp) :: i_sz = huge(one)
   ! Value of the integration of the Coulomb singularity 4\pi/V_BZ \int_BZ d^3q 1/q^2

  type(mc_t) :: mc
   ! Monte carlo integrator.

contains
  procedure :: init => vcgen_init                  ! Initialize the object
  procedure :: get_vc_sqrt => vcgen_get_vc_sqrt    ! Compute sqrt(vc(q,g))
  procedure :: free => vcgen_free                  ! Free
  !procedure :: print => vcgen_print
end type vcgen_t
!!***

 ! private stuff
 real(dp),parameter :: TOLQ0 = 1.d-3

CONTAINS  !========================================================================================
!!***

!!****f* m_vcoul/gw_icutcoul_to_mode
!! NAME
!! gw_icutcoul_to_mode
!!
!! FUNCTION
!! Convert gw_icutcoul_to_mode to mode string.
!!
!! SOURCE

subroutine gw_icutcoul_to_mode(gw_icutcoul, mode)

!Arguments ------------------------------------
 integer,intent(in) :: gw_icutcoul
 character(len=*),intent(out) :: mode

! *************************************************************************

 mode = 'NONE'
 if (gw_icutcoul == 0) mode = 'SPHERE'
 if (gw_icutcoul == 1) mode = 'CYLINDER'
 if (gw_icutcoul == 2) mode = 'SURFACE'
 if (gw_icutcoul == 3) mode = 'CRYSTAL'
 if (gw_icutcoul == 4) mode = 'ERF'
 if (gw_icutcoul == 5) mode = 'ERFC'
 if (gw_icutcoul == 6) mode = 'AUXILIARY_FUNCTION'
 if (gw_icutcoul == 7) mode = 'AUX_GB'
 if (gw_icutcoul == 14) mode = 'MINIBZ-ERF'
 if (gw_icutcoul == 15) mode = 'MINIBZ-ERFC'
 if (gw_icutcoul == 16) mode = 'MINIBZ'

end subroutine gw_icutcoul_to_mode
!!***

!!****f* m_vcoul/vcoul_init
!! NAME
!! vcoul_init
!!
!! FUNCTION
!! Perform general check and initialize the data type containing information on the cutoff technique
!! Note %vc_sqrt_resid and %i_sz_resid are simply initialized at the same value as %vc_sqrt and %i_sz
!!
!! INPUTS
!!  Gsph=Info of the G sphere.
!!  Qmesh=Info on the q-point sampling.
!!  Kmesh=Info on the k-point sampling.
!!  rcut=Cutoff radius for the cylinder.
!!  gw_icutcoul=Option of the cutoff technique.
!!  vcutgeo(3)= Info on the orientation and extension of the cutoff region.
!!  ng=Number of G-vectors to be used to describe the Coulomb interaction
!!  nqlwl=Number of point around Gamma for treatment of long-wavelength limit
!!  qlwl(3,nqlwl)= The nqlwl "small" q-points
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  vcp=Datatype gathering information on the Coulomb interaction.
!!
!! SOURCE

subroutine vcoul_init(vcp, Gsph, Cryst, Qmesh, Kmesh, rcut, gw_icutcoul, vcutgeo, ecut, ng, nqlwl, qlwl, comm)

!Arguments ------------------------------------
!scalars
 class(vcoul_t),intent(out) :: vcp
 integer,intent(in) :: ng,nqlwl, gw_icutcoul, comm
 real(dp),intent(in) :: rcut, ecut
 type(kmesh_t),target,intent(in) :: Kmesh, Qmesh
 type(gsphere_t),target,intent(in) :: Gsph
 type(crystal_t),intent(in) :: Cryst
!arrays
 real(dp),intent(in) :: qlwl(3,nqlwl),vcutgeo(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: nqibz, nqbz, nkbz, iqlwl, iq_ibz
 integer :: opt_cylinder,opt_surface,my_rank,nprocs
 real(dp) :: bz_geometry_factor,q0_vol, rcut2
 character(len=500) :: msg
 type(mc_t) :: mc
!arrays
 integer, contiguous, pointer :: gvec(:,:)
 real(dp) :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
 real(dp),allocatable :: vcoul(:,:),vcoul_lwl(:,:)
 real(dp),contiguous, pointer :: qibz(:,:), qbz(:,:)

! *************************************************************************

 ! === Test if the first q-point is zero ===
 ! FIXME this wont work if nqptdm/=0
 !if (normv(Qmesh%ibz(:,1),gmet,'G') < GW_TOLQ0)) STOP 'vcoul_init, non zero first point '

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 nqibz = qmesh%nibz; nqbz = qmesh%nbz
 qibz => qmesh%ibz; qbz => qmesh%bz
 nkbz = kmesh%nbz

 ! Save dimension and other useful quantities in vcp
 vcp%ng        = ng                   ! Number of G-vectors in the Coulomb matrix elements.
 vcp%nqibz     = nqibz                ! Number of irred q-point.
 vcp%nqlwl     = nqlwl                ! Number of small q-directions to deal with singularity and non Analytic behavior.
 vcp%rcut      = rcut                 ! Cutoff radius for cylinder.
 vcp%hcyl      = zero                 ! Length of finite cylinder (Rozzi"s method, default is Beigi).
 vcp%ucvol     = cryst%ucvol          ! Unit cell volume.
 vcp%rprimd    = Cryst%rprimd(:,:)    ! Dimensional direct lattice.
 vcp%boxcenter = zero                 ! Boxcenter at the moment is supposed to be at the origin.
 vcp%vcutgeo   = vcutgeo(:)           ! Info on the orientation and extension of the cutoff region.

 gvec => Gsph%gvec

 call gw_icutcoul_to_mode(gw_icutcoul, vcp%mode)

 ABI_MALLOC(vcp%qibz, (3, nqibz))
 vcp%qibz = Qmesh%ibz(:,:)
 ABI_MALLOC(vcp%qlwl, (3, nqlwl))
 vcp%qlwl = qlwl(:,:)

 ! ===============================================
 ! == Calculation of the FT of the Coulomb term ==
 ! ===============================================
 a1 = cryst%rprimd(:,1); a2 = cryst%rprimd(:,2); a3 = cryst%rprimd(:,3)
 b1 = two_pi * cryst%gprimd(:,1); b2 = two_pi * cryst%gprimd(:,2); b3 = two_pi * cryst%gprimd(:,3)

 ABI_MALLOC(vcoul    , (ng, nqibz))
 ABI_MALLOC(vcoul_lwl, (ng, nqlwl))

 select case (trim(vcp%mode))
 case ('MINIBZ', 'MINIBZ-ERFC', 'MINIBZ-ERF')
   call mc%init(cryst%rprimd, cryst%ucvol, cryst%gprimd, cryst%gmet, kmesh%kptrlatt)

   rcut2 = vcp%rcut**2
   do iq_ibz=1,nqibz
     call mc%integrate(vcp%mode, qibz(:, iq_ibz), ng, gvec, rcut2, nkbz, vcoul(:, iq_ibz), comm)
   end do

   ! Treat the limit q --> 0
   vcp%i_sz = vcoul(1, 1)

   do iqlwl=1,nqlwl
     call mc%integrate(vcp%mode, qlwl(:, iqlwl), ng, gvec, rcut2, nkbz, vcoul_lwl(:, iqlwl), comm)
   end do

   call mc%free()

 case ('SPHERE')
   ! A non-positive value of rcut activates the recipe of Spencer & Alavi, PRB 77, 193110 (2008) [[cite:Spencer2008]].
   if (vcp%rcut < tol12) then
     vcp%rcut = (cryst%ucvol * nkbz * 3.d0 / four_pi) ** third
     write(msg,'(2a,2x,f8.4,a)')ch10,' Using calculated rcut: ',vcp%rcut,' to have same volume as the BvK crystal'
     call wrtout(std_out, msg)
   end if
   vcp%vcutgeo = zero

   do iq_ibz=1,nqibz
     call cutoff_sphere(qibz(:,iq_ibz), ng, gvec, cryst%gmet, vcp%rcut, vcoul(:,iq_ibz))
   end do

   ! q-points for optical limit.
   do iqlwl=1,nqlwl
     call cutoff_sphere(qlwl(:,iqlwl), ng, gvec, cryst%gmet, vcp%rcut, vcoul_lwl(:,iqlwl))
   end do

   ! Treat the limit q --> 0
   ! The small cube is approximated by a sphere, while vc(q=0) = 2piR**2.
   ! if a single q-point is used, the expression for the volume is exact.
   vcp%i_sz = two_pi * vcp%rcut**2
   call vcp%print(unit=ab_out)

 case ('CYLINDER')
   call cylinder_setup(cryst, vcp%vcutgeo, vcp%hcyl, vcp%pdir, opt_cylinder)

   do iq_ibz=1,nqibz
     call cutoff_cylinder(qibz(:,iq_ibz), ng, gvec, vcp%rcut, vcp%hcyl, vcp%pdir,&
                          vcp%boxcenter, Cryst%rprimd, vcoul(:,iq_ibz), opt_cylinder, comm)
   end do

   ! q-points for optical limit.
   do iqlwl=1,nqlwl
     call cutoff_cylinder(qlwl(:,iqlwl), ng, gvec, vcp%rcut, vcp%hcyl, vcp%pdir,&
                          vcp%boxcenter, Cryst%rprimd, vcoul_lwl(:,iqlwl), opt_cylinder, comm)
   end do

   ! If Beigi, treat the limit q --> 0.
   if (opt_cylinder == 1) then
     call beigi_cylinder_limit(opt_cylinder, cryst, nqibz, nkbz, vcp%rcut, vcp%hcyl, vcp%boxcenter, vcp%pdir, vcp%i_sz)
   else
     ! In Rozzi's method the lim q+G --> 0 is finite.
     vcp%i_sz = vcoul(1,1)
   end if

   call vcp%print(unit=ab_out)

 case ('SURFACE')
   call surface_setup(cryst, vcp%vcutgeo, vcp%alpha, vcp%rcut, vcp%pdir, opt_surface)

   do iq_ibz=1,nqibz
     call cutoff_surface(qibz(:,iq_ibz), ng, gvec, cryst%gprimd, vcp%rcut, &
                         vcp%boxcenter, vcp%pdir, vcp%alpha, vcoul(:,iq_ibz), opt_surface)
   end do

   ! q-points for optical limit.
   do iqlwl=1,nqlwl
     call cutoff_surface(qlwl(:,iq_ibz), ng, gvec, cryst%gprimd, vcp%rcut, &
                         vcp%boxcenter, vcp%pdir, vcp%alpha, vcoul_lwl(:,iqlwl), opt_surface)
   end do

   ! If Beigi, treat the limit q --> 0.
   if (opt_surface == 1) then
     ! Integrate numerically in the plane close to 0
     call beigi_surface_limit(opt_surface, cryst, nqibz, nkbz, vcp%rcut, vcp%alpha, &
                              vcp%boxcenter, vcp%pdir, vcp%i_sz)
   else
     ! In Rozzi's method the lim q+G --> 0 is finite.
     vcp%i_sz=vcoul(1,1)
   end if

 case ('CRYSTAL', 'AUXILIARY_FUNCTION', "AUX_GB")
   do iq_ibz=1,nqibz
     call cmod_qpg(nqibz, iq_ibz, qibz, ng, gvec, cryst%gprimd, vcoul(:,iq_ibz))

     if (iq_ibz == 1) then
       ! The singularity is treated using vcoul_lwl.
       vcoul(1, iq_ibz) = zero
       vcoul(2:,iq_ibz) = four_pi / vcoul(2:,iq_ibz)**2
     else
       vcoul(:,iq_ibz) = four_pi / vcoul(:,iq_ibz)**2
     end if
   end do ! iq_ibz

   ! q-points for optical limit.
   do iqlwl=1,nqlwl
     call cmod_qpg(nqlwl, iqlwl, qlwl, ng, gvec, cryst%gprimd, vcoul_lwl(:,iqlwl))
   end do
   vcoul_lwl = four_pi/vcoul_lwl**2

   ! Treatment of 1/q^2 singularity

   if (vcp%mode == "CRYSTAL") then
     ! Analytic integration of 4pi/q^2 over the volume element:
     ! $4pi/V \int_V d^3q 1/q^2 =4pi bz_geometric_factor V^(-2/3)$
     ! i_sz=4*pi*bz_geometry_factor*q0_vol**(-two_thirds) where q0_vol= V_BZ/N_k
     ! bz_geometry_factor: sphere=7.79, fcc=7.44, sc=6.188, bcc=6.946, wz=5.255 (see gwa.pdf, appendix A.4)
     q0_vol = (two_pi) **3 / (nkbz*cryst%ucvol); bz_geometry_factor=zero
     vcp%i_sz = four_pi*7.44*q0_vol**(-two_thirds)

   else if (vcp%mode == "AUXILIARY_FUNCTION") then
     ! Numerical integration of the exact-exchange divergence through the
     ! auxiliary function of Carrier et al. PRB 75, 205126 (2007) [[cite:Carrier2007]].
     vcp%i_sz = carrier_isz(cryst, nqbz, qbz, rcut, comm)

   else if (vcp%mode == "AUX_GB") then
     ! We use the auxiliary function of a Gygi-Baldereschi variant [[cite:Gigy1986]]
     vcp%i_sz = gygi_baldereschi_isz(cryst, nqbz, qbz, ecut, ng, gvec)

   else
     ABI_ERROR(sjoin("Need treatment of 1/q^2 singularity! for mode", vcp%mode))
   end if

 case ('ERF')
   ! Modified long-range only Coulomb interaction thanks to the error function:
   ! * Vc = erf(r/rcut)/r
   ! * The singularity is treated using vcoul_lwl.
   do iq_ibz=1,nqibz
     call cmod_qpg(nqibz, iq_ibz, qibz, ng, gvec, cryst%gprimd, vcoul(:,iq_ibz))

     ! The Fourier transform of the error function reads
     if (iq_ibz == 1) then
       vcoul(1, iq_ibz) = zero
       vcoul(2:,iq_ibz) = four_pi/(vcoul(2:,iq_ibz)**2) *  EXP( -0.25d0 * (vcp%rcut*vcoul(2:,iq_ibz))**2 )
     else
       vcoul(:,iq_ibz)  = four_pi/(vcoul(:, iq_ibz)**2) *  EXP( -0.25d0 * (vcp%rcut*vcoul(: ,iq_ibz))**2 )
     end if
   end do

   ! q-points for optical limit.
   do iqlwl=1,nqlwl
     call cmod_qpg(nqlwl, iqlwl, qlwl, ng, gvec, cryst%gprimd, vcoul_lwl(:,iqlwl))
   end do
   vcoul_lwl = four_pi/(vcoul_lwl**2) *  EXP( -0.25d0 * (vcp%rcut*vcoul_lwl)**2 )

   ! === Treat 1/q^2 singularity ===
   ! * We use the auxiliary function from PRB 75, 205126 (2007) [[cite:Carrier2007]]
   vcp%i_sz = carrier_isz(cryst, nqbz, qbz, rcut, comm)

 case ('ERFC')
   ! * Use a modified short-range only Coulomb interaction thanks to the complementary error function:
   !   $ V_c = [1-erf(r/r_{cut})]/r $
   ! * The Fourier transform of the error function reads
   !   vcoul=four_pi/(vcoul**2) * ( 1.d0 - exp( -0.25d0 * (vcp%rcut*vcoul)**2 ) )
   do iq_ibz=1,nqibz
     call cmod_qpg(nqibz, iq_ibz, qibz,ng, gvec, cryst%gprimd, vcoul(:,iq_ibz))

     if (iq_ibz == 1) then
       vcoul(1 ,iq_ibz) = zero
       vcoul(2:,iq_ibz) = four_pi/(vcoul(2:,iq_ibz)**2) * ( one - EXP( -0.25d0 * (vcp%rcut*vcoul(2:,iq_ibz))**2 ) )
     else
       vcoul(:, iq_ibz) = four_pi/(vcoul(:, iq_ibz)**2) * ( one - EXP( -0.25d0 * (vcp%rcut*vcoul(:, iq_ibz))**2 ) )
     end if
   end do ! iq_ibz

   ! q-points for optical limit.
   do iqlwl=1,nqlwl
     call cmod_qpg(nqlwl, iqlwl, qlwl, ng, gvec, cryst%gprimd, vcoul_lwl(:,iqlwl))
   end do
   vcoul_lwl = four_pi/(vcoul_lwl**2) * ( one - EXP( -0.25d0 * (vcp%rcut*vcoul_lwl)**2 ) )

   ! === Treat 1/q^2 singularity ===
   ! * There is NO singularity in this case.
   vcp%i_sz = pi * vcp%rcut**2 ! Final result stored here

 case default
   ABI_BUG(sjoin('Unsupported cutoff mode:', vcp%mode))
 end select

 call wrtout(std_out, sjoin("vcp%i_sz", ftoa(vcp%i_sz)))
 vcp%i_sz_resid = vcp%i_sz

 ! Store final results in complex array as Rozzi's cutoff can give real negative values
 ABI_MALLOC(vcp%vc_sqrt, (ng, nqibz))
 ABI_MALLOC(vcp%vc_sqrt_resid, (ng, nqibz))
 vcp%vc_sqrt = CMPLX(vcoul, zero)
 vcp%vc_sqrt = SQRT(vcp%vc_sqrt)
 vcp%vc_sqrt_resid = vcp%vc_sqrt
 ABI_FREE(vcoul)

 ABI_MALLOC(vcp%vcqlwl_sqrt, (ng, nqlwl))
 vcp%vcqlwl_sqrt = CMPLX(vcoul_lwl, zero)
 vcp%vcqlwl_sqrt = SQRT(vcp%vcqlwl_sqrt)
 ABI_FREE(vcoul_lwl)

 call vcp%print(unit=std_out)

end subroutine vcoul_init
!!***


!!****f* m_vcoul/cylinder_setup
!! NAME
!!  cylinder_setup
!!
!! FUNCTION
!!
!! SOURCE

subroutine cylinder_setup(cryst, vcutgeo, hcyl, pdir, opt_cylinder)

 type(crystal_t),intent(in) :: cryst
 real(dp),intent(in) :: vcutgeo(3)
 real(dp),intent(out) :: hcyl
 integer,intent(out) :: pdir(3), opt_cylinder

!Local variables-------------------------------
 integer :: ii
 real(dp),parameter :: tol999 = 999.0
 real(dp) :: check

! *************************************************************************

 ABI_CHECK(count(abs(vcutgeo) > tol6) == 1, 'Wrong cutgeo for cylinder')

 ! Beigi's method is the default one, i.e infinite cylinder of radius rcut.
 ! Use negative values to use Rozzi's method with finite cylinder of extent hcyl.
 opt_cylinder = 1; hcyl = zero; pdir(:) = 0
 do ii=1,3
   check = vcutgeo(ii)
   if (abs(check) > tol6) then
     pdir(ii) = 1
     if (check < zero) then
       ! use Rozzi's method.
       hcyl = ABS(check) *SQRT(SUM(cryst%rprimd(:,ii)**2))
       opt_cylinder = 2
       ! Check to enter the infinite Rozzi treatment
       if(vcutgeo(3) <= -tol999) hcyl = tol12
     end if
   end if
 end do

 ABI_CHECK((count(pdir == 1) == 1), 'Wrong pdir for cylinder')
 if (pdir(3) /= 1) then
   ABI_ERROR("The cylinder must be along the z-axis")
 end if

end subroutine cylinder_setup
!!***

!!****f* m_vcoul/surface_setup
!! NAME
!!  surface_setup
!!
!! FUNCTION
!!
!! SOURCE

subroutine surface_setup(cryst, vcutgeo, alpha, rcut, pdir, opt_surface)

 type(crystal_t),intent(in) :: cryst
 real(dp),intent(in) :: vcutgeo(3)
 real(dp),intent(out) :: alpha(3)
 real(dp),intent(inout) :: rcut
 integer,intent(out) :: pdir(3), opt_surface

!Local variables-------------------------------
 integer :: ii
 real(dp) :: check
 character(len=500) :: msg

! *************************************************************************

 ABI_CHECK(count(vcutgeo /= zero) == 2, "Wrong vcutgeo")

 ! Default is Beigi's method.
 opt_surface = 1; if (any(vcutgeo < zero)) opt_surface = 2
 pdir(:) = zero; alpha(:)=zero
 do ii=1,3
   check = vcutgeo(ii)
   if (abs(check) > zero) then
     ! Use Rozzi's method with a finite surface along x-y
     pdir(ii) = 1
     if (check < zero) alpha(ii) = normv(check * cryst%rprimd(:,ii), cryst%rmet, 'R')
   end if
 end do

 ! In Beigi's method, the surface must be along x-y and R must be L_Z/2.
 if (opt_surface == 1) then
   msg = "2D Beigi method, the periodicity must be in the x-y plane. Modify vcutgeo and/or your geometry."
   ABI_CHECK(all(pdir == [1, 1, 0]), msg)
   rcut = half*SQRT(DOT_PRODUCT(cryst%rprimd(:,3), cryst%rprimd(:,3)))
 end if

end subroutine surface_setup
!!***

!!****f* m_vcoul/integratefaux
!! NAME
!!  integratefaux
!!
!! FUNCTION
!!
!! SOURCE

real(dp) function integratefaux(rcut, gprimd, ucvol, comm)

 real(dp),intent(in) :: rcut, gprimd(3,3), ucvol
 integer,intent(in) :: comm

!Local variables-------------------------------
 integer,parameter :: nref = 3, nq = 50
 integer :: ierr,iq,iqx1,iqy1,iqz1,iqx2,iqy2,iqz2,miniqy1,maxiqy1,nqhalf, nprocs, my_rank
 real(dp) :: invnq,invnq3,qq,weightq,weightxy,weightxyz
 real(dp) :: qq1(3),qq2(3),bb4sinpiqq_2(3,nq),sin2piqq(nq),bb4sinpiqq2_2(3,0:nq),sin2piqq2(3,0:nq)
 real(dp) :: b1(3), b2(3), b3(3), bb(3)
 real(dp) :: b1b1,b2b2,b3b3,b1b2,b2b3,b3b1

! *************************************************************************

 ! nq is the number of sampling points along each axis for the numerical integration
 ! nref is the area where the mesh is refined

 integratefaux = zero
 invnq = one/DBLE(nq)
 invnq3 = invnq**3
 nqhalf = nq/2

 b1 = two_pi * gprimd(:,1); b2 = two_pi * gprimd(:,2); b3 = two_pi * gprimd(:,3)
 b1b1 = dot_product(b1, b1); b2b2 = dot_product(b2, b2); b3b3 = dot_product(b3, b3)
 bb(1) = b1b1; bb(2) = b2b2; bb(3) = b3b3
 b1b2 = dot_product(b1, b2); b2b3 = dot_product(b2, b3); b3b1 = dot_product(b3, b1)

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! In order to speed up the calculation, precompute the sines
 do iq=1,nq
   qq=DBLE(iq)*invnq-half
   bb4sinpiqq_2(:,iq)=bb(:)*four*SIN(pi*qq)**2 ; sin2piqq(iq)=SIN(two_pi*qq)
 end do

 do iqx1=1,nq
   if (modulo(iqx1, nprocs) /= my_rank) cycle ! MPI parallelism
   qq1(1)=DBLE(iqx1)*invnq-half
   ! Here take advantage of the q <=> -q symmetry:
   ! arrange the sampling of qx, qy space to avoid duplicating calculations. Need weights to do this ...
   !do iqy1=1,nq
   miniqy1 = nqhalf + 1; maxiqy1 = nq
   if (iqx1 >= nqhalf) miniqy1 = nqhalf
   if (iqx1 > nqhalf .and. iqx1 < nq) maxiqy1 = nq - 1

   do iqy1=miniqy1,maxiqy1
     qq1(2) = DBLE(iqy1)*invnq - half
     ! By default, the factor of two is for the q <=> -q symmetry
     weightq = invnq3*two
     ! But not all qx qy lines have a symmetric one ...
     if( (iqx1 == nqhalf .or. iqx1 == nq) .and. (iqy1 == nqhalf .or. iqy1 == nq)) weightq = weightq*half

     do iqz1=1,nq
       qq1(3) = DBLE(iqz1)*invnq - half

       ! Refine the mesh for the point close to the origin
       if( abs(iqx1-nqhalf) <= nref .and. abs(iqy1-nqhalf) <= nref .and. abs(iqz1-nqhalf) <= nref ) then
         ! Note that the set of point is symmetric around the central point, while weights are taken into account
         do iq=0,nq
           qq2(:) = qq1(:)+ (DBLE(iq)*invnq-half)*invnq
           bb4sinpiqq2_2(:,iq) =bb(:)*four*SIN(pi*qq2(:))**2; sin2piqq2(:,iq)=SIN(two_pi*qq2(:))
         end do
         do iqx2=0,nq
           qq2(1)=qq1(1) + (DBLE(iqx2)*invnq-half ) *invnq
           do iqy2=0,nq
             qq2(2)=qq1(2) + (DBLE(iqy2)*invnq-half ) *invnq
             weightxy=invnq3*weightq
             if (iqx2 == 0 .or. iqx2 == nq) weightxy = weightxy*half
             if (iqy2 == 0 .or. iqy2 == nq) weightxy = weightxy*half
             do iqz2=0,nq
               qq2(3) = qq1(3) + (DBLE(iqz2)*invnq - half) * invnq
               weightxyz = weightxy
               if (iqz2 == 0 .or. iqz2 == nq) weightxyz = weightxy*half
               !
               ! Treat the remaining divergence in the origin as if it would be a spherical integration of 1/q^2
               if (iqx1/=nqhalf .or. iqy1/=nqhalf .or. iqz1/=nqhalf .or. &
                   iqx2/=nqhalf .or. iqy2/=nqhalf .or. iqz2/=nqhalf ) then
                 !integratefaux=integratefaux+ faux(qq2, rcut, b1, b2, b3) *invnq**6
                 integratefaux = integratefaux + &
                   faux_fast(qq2, bb4sinpiqq2_2(1,iqx2), bb4sinpiqq2_2(2,iqy2), bb4sinpiqq2_2(3,iqz2), &
                             sin2piqq2(1,iqx2), sin2piqq2(2,iqy2), sin2piqq2(3,iqz2), rcut, b1, b2, b3) * weightxyz
               else
                  integratefaux = integratefaux + 7.7955* ((two_pi)**3/ucvol*invnq3*invnq3 )**(-2./3.) *invnq3*invnq3
               end if
             end do
           end do
         end do
       else
        ! integratefaux=integratefaux+faux(qq1, rcut, b1, b2, b3)*invnq**3
        integratefaux = integratefaux + &
          faux_fast(qq1, bb4sinpiqq_2(1,iqx1), bb4sinpiqq_2(2,iqy1), bb4sinpiqq_2(3,iqz1), &
                    sin2piqq(iqx1), sin2piqq(iqy1), sin2piqq(iqz1), rcut, b1, b2, b3) * weightq
       end if
     end do
   end do
 end do

 call xmpi_sum(integratefaux, comm, ierr)

end function integratefaux
!!***

real(dp) pure function faux(qq, rcut, b1, b2, b3)

!Arguments ------------------------------------
 real(dp),intent(in) :: qq(3)
 real(dp),intent(in) :: rcut
 real(dp),intent(in) :: b1(3), b2(3), b3(3)

!Local variables-------------------------------
 real(dp) :: bb4sinpiqq1_2, bb4sinpiqq2_2, bb4sinpiqq3_2, sin2piqq1, sin2piqq2, sin2piqq3
 real(dp) :: b1b1,b2b2,b3b3

! *************************************************************************

 b1b1 = dot_product(b1, b1); b2b2 = dot_product(b2, b2); b3b3 = dot_product(b3, b3)

 bb4sinpiqq1_2 = b1b1 * four * SIN(pi*qq(1))**2
 bb4sinpiqq2_2 = b2b2 * four * SIN(pi*qq(2))**2
 bb4sinpiqq3_2 = b3b3 * four * SIN(pi*qq(3))**2
 sin2piqq1 = SIN(two_pi*qq(1))
 sin2piqq2 = SIN(two_pi*qq(2))
 sin2piqq3 = SIN(two_pi*qq(3))

 faux = faux_fast(qq, bb4sinpiqq1_2, bb4sinpiqq2_2, bb4sinpiqq3_2, sin2piqq1, sin2piqq2, sin2piqq3, rcut, b1, b2, b3)

end function faux

real(dp) pure function faux_fast(qq, bb4sinpiqq1_2, bb4sinpiqq2_2, bb4sinpiqq3_2, sin2piqq1, sin2piqq2, sin2piqq3, &
                                 rcut, b1, b2, b3)

!Arguments ------------------------------------
 real(dp),intent(in) :: qq(3)
 real(dp),intent(in) :: bb4sinpiqq1_2, bb4sinpiqq2_2, bb4sinpiqq3_2, sin2piqq1, sin2piqq2, sin2piqq3
 real(dp),intent(in) :: rcut
 real(dp),intent(in) :: b1(3), b2(3), b3(3)

!Local variables-------------------------------
 real(dp) :: b1b2,b2b3,b3b1
 b1b2 = dot_product(b1, b2); b2b3 = dot_product(b2, b3); b3b1 = dot_product(b3, b1)

! *************************************************************************

 faux_fast = bb4sinpiqq1_2 + bb4sinpiqq2_2 + bb4sinpiqq3_2 &
      +two*( b1b2 * sin2piqq1*sin2piqq2 &
            +b2b3 * sin2piqq2*sin2piqq3 &
            +b3b1 * sin2piqq3*sin2piqq1 &
           )

 if (rcut > tol6) then
   faux_fast = two_pi*two_pi/faux_fast * exp( -0.25d0*rcut**2* sum( ( qq(1)*b1(:)+qq(2)*b2(:)+qq(3)*b3(:) )**2 ) )
 else
   faux_fast = two_pi*two_pi/faux_fast
 endif

end function faux_fast
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/adapt_nmc
!! NAME
!! adapt_nmc
!!
!! FUNCTION
!! Empirical law to decrease the Monte Carlo sampling
!! for large q+G, for which the accuracy is not an issue

integer pure function adapt_nmc(nmc_max, qpg2) result(nmc)

!Arguments ------------------------------------
 integer,intent(in)  :: nmc_max
 real(dp),intent(in) :: qpg2

! *************************************************************************

 nmc = NINT( nmc_max / ( 1.0_dp + 1.0_dp * qpg2**6 ) )
 nmc = MIN(nmc_max, nmc)
 nmc = MAX(1, nmc)

end function adapt_nmc
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
!! SOURCE

subroutine vcoul_plot(Vcp, Qmesh, Gsph, ng, vc, comm)

!Arguments ------------------------------------
!scalars
 class(vcoul_t),intent(in) :: Vcp
 integer,intent(in) :: ng,comm
 type(kmesh_t),intent(in) :: Qmesh
 type(gsphere_t),intent(in) :: Gsph
!arrays
 real(dp),intent(in) :: vc(ng,Qmesh%nibz)

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: icount,idx_Sm1G,ierr,ig,igs,ii,iq_bz,iq_ibz,iqg,ir,isym,itim
 integer :: my_start,my_stop,nqbz,nqibz,nr,ntasks,my_rank,unt
 real(dp) :: arg,fact,l1,l2,l3,lmax,step,tmp,vcft,vc_bare
 character(len=500) :: msg
 character(len=fnlen) :: filnam
!arrays
 integer,allocatable :: insort(:)
 real(dp) :: b1(3),b2(3),b3(3),gmet(3,3),gprimd(3,3),qbz(3),qpgc(3)
 real(dp),allocatable :: qpg_mod(:),rr(:,:,:),vcr(:,:),vcr_cut(:,:)

!************************************************************************

 if (trim(Vcp%mode) /= 'CYLINDER') RETURN

 my_rank = xmpi_comm_rank(comm)

 nqibz=Vcp%nqibz; nqbz=Qmesh%nbz
 gmet=Gsph%gmet; gprimd=Gsph%gprimd

 b1(:)=two_pi*gprimd(:,1)
 b2(:)=two_pi*gprimd(:,2)
 b3(:)=two_pi*gprimd(:,3)

 ! Compare in Fourier space the true Coulomb with the cutted one.
 if (my_rank == master) then
   ABI_MALLOC(insort, (nqibz * ng))
   ABI_MALLOC(qpg_mod, (nqibz * ng))
   iqg = 1
   do iq_ibz=1,nqibz
     do ig=1,ng
       qpg_mod(iqg) = normv(Qmesh%ibz(:,iq_ibz) + Gsph%gvec(:,ig), gmet,'g')
       insort(iqg) = iqg; iqg = iqg + 1
     end do
   end do
   call sort_dp(nqibz * ng, qpg_mod, insort, tol14)

   filnam='_VCoulFT_'
   call isfile(filnam, 'new')
   if (open_file(filnam, msg, newunit=unt, status='new', form='formatted') /= 0) then
     ABI_ERROR(msg)
   end if
   write(unt,'(a,i3,a,i6,a)')&
    '#   |q+G|       q-point (Tot no.',nqibz,')        Gvec (',ng,')     vc_bare(q,G)    vc_cutoff(q,G) '

   do iqg=1,nqibz*ng
     iq_ibz = (insort(iqg) - 1) / ng + 1
     ig = (insort(iqg)) - (iq_ibz-1) * ng
     vc_bare = zero
     if (qpg_mod(iqg) > tol16) vc_bare = four_pi / qpg_mod(iqg) ** 2
     write(unt,'(f12.6,2x,3f8.4,2x,3i6,2x,2es14.6)')&
       qpg_mod(iqg), Qmesh%ibz(:,iq_ibz), Gsph%gvec(:,ig), vc_bare, vc(ig, iq_ibz)
   end do

   close(unt)
   ABI_FREE(insort)
   ABI_FREE(qpg_mod)
 end if ! my_rank==master

 ! Fourier transform back to real space just to check cutoff implementation.
 ntasks= nqbz * ng
 call xmpi_split_work(ntasks, comm, my_start, my_stop)

 l1 = SQRT(SUM(Vcp%rprimd(:,1)**2))
 l2 = SQRT(SUM(Vcp%rprimd(:,2)**2))
 l3 = SQRT(SUM(Vcp%rprimd(:,3)**2))

 nr = 50
 lmax=MAX(l1,l2,l3) ; step=lmax/(nr-1)
 fact = one / (Vcp%ucvol * nqbz)

 ! numb coding
 ABI_CALLOC(rr, (3, nr, 3))
 do ii=1,3
   do ir=1,nr
     rr(ii,ir,ii)=(ir-1)*step
   end do
 end do

 ABI_CALLOC(vcr, (nr, 3))
 ABI_CALLOC(vcr_cut, (nr, 3))

 do iq_bz=1,nqbz
   call Qmesh%get_BZ_item(iq_bz, qbz, iq_ibz, isym, itim)
   if (ABS(qbz(1))<0.01) qbz(1)=zero
   if (ABS(qbz(2))<0.01) qbz(2)=zero
   if (ABS(qbz(3))<0.01) qbz(3)=zero
   igs=1; if (ALL(qbz(:)==zero)) igs=2
   do ig=igs,ng
     icount=ig+(iq_bz-1)*ng
     if (icount < my_start .or. icount > my_stop) CYCLE
     idx_Sm1G = Gsph%rottbm1(ig,itim,isym) ! IS{^-1}G
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

 if (my_rank == master) then
   filnam='_VCoulR_'
   call isfile(filnam, 'new')
   if (open_file(filnam,msg,newunit=unt,status='new',form='formatted') /= 0) then
     ABI_ERROR(msg)
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
!! SOURCE

subroutine vcoul_print(Vcp, unit, prtvol, mode_paral)

!Arguments ------------------------------------
!scalars
 class(vcoul_t),intent(in) :: Vcp
 integer,intent(in),optional :: prtvol,unit
 character(len=4),intent(in),optional :: mode_paral

!Local variables-------------------------------
!scalars
 integer :: ii,my_unt,my_prtvol,iqlwl
 character(len=4) :: my_mode
 character(len=500) :: msg

! *************************************************************************

 my_unt    =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_mode   ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral
 my_prtvol=0       ; if (PRESENT(prtvol    )) my_prtvol=prtvol

 select case (Vcp%mode)

 case ('MINIBZ')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',trim(Vcp%mode)
   call wrtout(my_unt,msg,my_mode)

 case ('MINIBZ-ERF')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',trim(Vcp%mode)
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(5a,f10.4,3a,f10.2,3a,3f10.5,2a)')ch10,&
     ' === Error function cutoff === ',ch10,ch10,&
     '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10
   call wrtout(my_unt,msg,my_mode)

 case ('MINIBZ-ERFC')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',trim(Vcp%mode)
   call wrtout(my_unt,msg,my_mode)
   write(msg,'(5a,f10.4,3a,f10.2,3a,3f10.5,2a)')ch10,&
     ' === Complement Error function cutoff === ',ch10,ch10,&
     '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10
   call wrtout(my_unt,msg,my_mode)

 case ('SPHERE')
   write(msg,'(5a,f10.4,3a,f10.2,3a,3f10.5,2a)')ch10,&
    ' === Spherical cutoff === ',ch10,ch10,&
    '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10,&
    '  Volume of the sphere .. ',four_pi/three*Vcp%rcut**3,' [Bohr^3] '
    !FB: This has no meaning here! &  '  Sphere centered at .... ',Vcp%boxcenter,' (r.l.u) ',ch10
    !MG It might be useful if the system is not centered on the origin because in this case the
    !   matrix elements of the Coulomb have to be multiplied by a phase depending on boxcenter.
    !   I still have to decide if it is useful to code this possibility and which variable use to
    !   define the center (boxcenter is used in the tddft part).
   call wrtout(my_unt,msg,my_mode)

 case ('CYLINDER')
   ii=imin_loc(ABS(Vcp%pdir-1))
   write(msg,'(5a,f10.4,3a,i2,2a,3f10.2,a)')ch10,&
     ' === Cylindrical cutoff === ',ch10,ch10,&
     '  Cutoff radius ............... ',Vcp%rcut,' [Bohr] ',ch10,&
     '  Axis parallel to direction... ',ii,ch10,&
     '  Passing through point ....... ',Vcp%boxcenter,' (r.l.u) '
   call wrtout(my_unt,msg,my_mode)

   write(msg,'(2a)')'  Infinite length  ....... ',ch10
   if (Vcp%hcyl/=zero) write(msg,'(a,f8.5,2a)')'  Finite length of ....... ',Vcp%hcyl,' [Bohr] ',ch10
   call wrtout(my_unt,msg,my_mode)

 case ('SURFACE')
   write(msg,'(5a,f10.4,3a,3f10.2,2a)')ch10,&
     ' === Surface cutoff === ',ch10,ch10,&
     '  Cutoff radius .................... ',Vcp%rcut,' [Bohr] ',ch10,&
     '  Central plane passing through .... ',Vcp%boxcenter,' (r.l.u) ',ch10
   call wrtout(my_unt,msg,my_mode)
   !write(msg,'(a)')'  Infinite length  .......'
   !if (Vcp%hcyl/=zero) write(msg,'(a,f8.5,a)')'  Finite length of .......',Vcp%hcyl,' [Bohr] '
   !call wrtout(my_unt,msg,my_mode)

 case ('AUXILIARY_FUNCTION')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',trim(Vcp%mode)
   call wrtout(my_unt,msg,my_mode)

 case ('AUX_GB')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',trim(Vcp%mode)
   call wrtout(my_unt,msg,my_mode)

 case ('CRYSTAL')
   write(msg,'(3a)')ch10,' vcoul_init : cutoff-mode = ',trim(Vcp%mode)
   call wrtout(my_unt,msg,my_mode)

 case ('ERF')
   write(msg,'(5a,f10.4,3a,f10.2,3a,3f10.5,2a)')ch10,&
     ' === Error function cutoff === ',ch10,ch10,&
     '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10
   call wrtout(my_unt,msg,my_mode)

 case ('ERFC')
   write(msg,'(5a,f10.4,3a,f10.2,3a,3f10.5,2a)')ch10,&
     ' === Complement Error function cutoff === ',ch10,ch10,&
     '  Cutoff radius ......... ',Vcp%rcut,' [Bohr] ',ch10
   call wrtout(my_unt,msg,my_mode)

 case default
   ABI_BUG(sjoin('Unknown cutoff mode: ', Vcp%mode))
 end select

 if (Vcp%nqlwl > 0) then
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
!!  Free memory
!!
!! SOURCE

subroutine vcoul_free(Vcp)

!Arguments ------------------------------------
 class(vcoul_t),intent(inout) :: Vcp

! *************************************************************************

 ABI_SFREE(Vcp%qibz)
 ABI_SFREE(Vcp%qlwl)
 ABI_SFREE(Vcp%vc_sqrt)
 ABI_SFREE(Vcp%vc_sqrt_resid)
 ABI_SFREE(Vcp%vcqlwl_sqrt)

end subroutine vcoul_free
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/mc_init
!! NAME
!! mc_init
!!
!! FUNCTION
!!
!! SOURCE

subroutine mc_init(mc, rprimd, ucvol, gprimd, gmet, kptrlatt)

!Arguments ------------------------------------
 class(mc_t),intent(out) :: mc
 real(dp),intent(in) :: rprimd(3,3), ucvol, gprimd(3,3), gmet(3,3)
 integer,intent(in) :: kptrlatt(3,3)

!Local variables-------------------------------
 integer,parameter :: ncell=3
 integer :: nseed, i1,i2,i3,imc
 real(dp) :: lmin,vlength, ucvol_sc
 real(dp) :: rprimd_sc(3,3),gprimd_sc(3,3),gmet_sc(3,3),rmet_sc(3,3), qcart2red(3,3)
 real(dp) :: qtmp(3),qmin(3),qmin_cart(3)
 integer, allocatable :: seed(:)

! *************************************************************************

 mc%gmet = gmet
 mc%ucvol = ucvol

 ! Supercell defined by the k-mesh
 rprimd_sc(:,:) = MATMUL(rprimd, kptrlatt)
 call metric(gmet_sc, gprimd_sc, -1, rmet_sc, rprimd_sc, ucvol_sc)

 qcart2red(:,:) = two_pi * gprimd(:,:)
 call matrginv(qcart2red, 3, 3)

 ! Find the largest sphere inside the miniBZ
 ! in order to integrate the divergence analytically
 mc%q0sph = HUGE(one)
 do i1 = -ncell+1, ncell
   qtmp(1) = dble(i1) * 0.5_dp
   do i2 = -ncell+1, ncell
     qtmp(2) = dble(i2) * 0.5_dp
     do i3 = -ncell+1, ncell
       qtmp(3) = dble(i3) * 0.5_dp
       if (i1 == 0 .AND. i2 == 0 .AND. i3 == 0) cycle
       vlength = normv(qtmp, gmet_sc, 'G')
       if (vlength < mc%q0sph) mc%q0sph = vlength
     enddo
   enddo
 enddo

 ! Setup the random vectors for the Monte Carlo sampling of the miniBZ at q = 0
 mc%nmc_max = 2500000
 ABI_MALLOC(mc%qran,(3, mc%nmc_max))
 call random_seed(size=nseed)
 ABI_MALLOC(seed, (nseed))
 do i1=1,nseed
   seed(i1) = NINT(SQRT(DBLE(i1) * 103731))
 end do
 call random_seed(put=seed)
 call random_number(mc%qran)
 ABI_FREE(seed)

 ! Overide the first "random vector" with 0
 mc%qran(:,1) = zero

 ! Fold qran into the Wignez-Seitz cell around q = 0
 do imc=2,mc%nmc_max
   lmin = HUGE(one)
   do i1 = -ncell+1, ncell
     qtmp(1) = mc%qran(1,imc) + dble(i1)
     do i2 = -ncell+1, ncell
       qtmp(2) = mc%qran(2,imc) + dble(i2)
       do i3 = -ncell+1, ncell
         qtmp(3) = mc%qran(3,imc) + dble(i3)
         vlength = normv(qtmp, gmet_sc, 'G')
         if (vlength < lmin) then
           lmin = vlength
           ! Get the q-vector in cartesian coordinates
           qmin_cart(:) = two_pi * MATMUL( gprimd_sc(:,:) , qtmp )
           ! Transform it back to the reciprocal space
           qmin(:) = MATMUL( qcart2red , qmin_cart )
         end if
       enddo
     enddo
   enddo

   mc%qran(:,imc) = qmin(:)
 enddo

end subroutine mc_init
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/mc_integrate
!! NAME
!! mc_integrate
!!
!! FUNCTION
!!
!! SOURCE

subroutine mc_integrate(mc, mode, qibz, ng, gvec, rcut2, nkbz, vcoul, comm)

!Arguments ------------------------------------
 class(mc_t),intent(in) :: mc
 real(dp),intent(in) :: rcut2
 integer,intent(in) :: nkbz, ng, comm
 character(len=*),intent(in) :: mode
 real(dp),intent(in) :: qibz(3)
 integer,intent(in) :: gvec(3, ng)
 real(dp),intent(out) :: vcoul(ng)

!Local variables-------------------------------
 integer,parameter :: master = 0
 integer :: ig, ig0, imc, nmc, my_rank, nprocs, ierr
 logical :: q_is_gamma
 real(dp)  :: qpg2, qpg(3)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 q_is_gamma = all(abs(qibz) < tol16)

 ! Find index of G=0 in gvec.
 ig0 = -1
 do ig=1,ng
   if (all(gvec(:, ig) == 0)) then
      ig0 = ig; exit
   end if
 end do
 ABI_CHECK(ig0 /= -1, "Cannot find G=0 in gvec!")

 vcoul = zero

 select case (trim(mode))
 case('MINIBZ')

   do ig=1,ng
     if (mod(ig, nprocs) /= my_rank) cycle ! MPI parallelism.
     if (q_is_gamma .and. ig == ig0) cycle
     qpg(:) = qibz(:) + gvec(:,ig)
     qpg2 = normv(qpg, mc%gmet, 'G')**2
     nmc = adapt_nmc(mc%nmc_max, qpg2)
     do imc=1,nmc
       qpg(:) = qibz(:) + gvec(:,ig) + mc%qran(:,imc)
       qpg2 = normv(qpg, mc%gmet, 'G')**2
       vcoul(ig) = vcoul(ig) + four_pi / qpg2 / REAL(nmc, dp)
     end do
   end do ! ig

   if (q_is_gamma .and. my_rank == master) then
     ! Override ig0 component
     vcoul(ig0) = four_pi**2 * nkbz * mc%ucvol / ( 8.0_dp * pi**3 ) * mc%q0sph
     do imc=1,mc%nmc_max
       qpg(:) = qibz(:) + gvec(:,ig0) + mc%qran(:,imc)
       qpg2 = normv(qpg, mc%gmet, 'G')**2
       if (qpg2 > mc%q0sph ** 2) vcoul(ig0) = vcoul(ig0) + four_pi / qpg2 / REAL(mc%nmc_max, dp)
     end do
   end if

 case('MINIBZ-ERFC')

   do ig=1,ng
     if (mod(ig, nprocs) /= my_rank) cycle ! MPI parallelism.
     if (q_is_gamma .and. ig == ig0) cycle
     qpg(:) = qibz(:) + gvec(:,ig)
     qpg2 = normv(qpg, mc%gmet, 'G')**2
     nmc = adapt_nmc(mc%nmc_max, qpg2)
     do imc=1,nmc
       qpg(:) = qibz(:) + gvec(:,ig) + mc%qran(:,imc)
       qpg2 = normv(qpg, mc%gmet, 'G')**2
       vcoul(ig) = vcoul(ig) + four_pi / qpg2 / REAL(nmc,dp) * (  one - EXP( -0.25d0 * rcut2 * qpg2 ) )
     end do
   end do ! ig

   if (q_is_gamma .and. my_rank == master) then
     ! Override ig0 component
     vcoul(ig0) = four_pi**2 * nkbz * mc%ucvol / ( 8.0_dp * pi**3 ) &
        * ( mc%q0sph - SQRT(pi/rcut2) * abi_derf(0.5_dp*SQRT(rcut2)*mc%q0sph) )
     do imc=1,mc%nmc_max
       qpg(:) = qibz(:) + gvec(:,ig0) + mc%qran(:,imc)
       qpg2 = normv(qpg, mc%gmet, 'G')**2
       if (qpg2 > mc%q0sph**2) then
         vcoul(ig0) = vcoul(ig0) + four_pi / qpg2 / REAL(mc%nmc_max,dp) * (one - EXP( -0.25d0 * rcut2 * qpg2))
       end if
     end do
   end if

 case('MINIBZ-ERF')

   do ig=1,ng
     if (mod(ig, nprocs) /= my_rank) cycle ! MPI parallelism.
     if (q_is_gamma .and. ig == ig0) cycle
     qpg(:) = qibz(:) + gvec(:,ig)
     qpg2 = normv(qpg, mc%gmet, 'G')**2
     nmc = adapt_nmc(mc%nmc_max, qpg2)
     do imc=1,nmc
       qpg(:) = qibz(:) + gvec(:,ig) + mc%qran(:,imc)
       qpg2 = normv(qpg, mc%gmet, 'G')**2
       vcoul(ig) = vcoul(ig) + four_pi / qpg2 / REAL(nmc,dp) * EXP( -0.25d0 * rcut2 * qpg2 )
     end do
   end do ! ig

   if (q_is_gamma .and. my_rank == master) then
     ! Override ig=ig0 component
     vcoul(ig0) = four_pi**2 * nkbz * mc%ucvol / ( 8.0_dp * pi**3 ) * SQRT(pi/rcut2) * abi_derf(0.5_dp*SQRT(rcut2)*mc%q0sph)

     do imc=1,mc%nmc_max
       qpg(:) = qibz(:) + gvec(:,ig0) + mc%qran(:,imc)
       qpg2 = normv(qpg, mc%gmet, 'G')**2
       if (qpg2 > mc%q0sph**2) then
         vcoul(ig0) = vcoul(ig0) + four_pi / qpg2 / REAL(mc%nmc_max,dp) *  EXP( -0.25d0 * rcut2 * qpg2 )
       end if
     end do
   end if

 case default
   ABI_ERROR(sjoin("Invalid mode:", mode))
 end select

 ! Collect result on each MPI proc.
 call xmpi_sum(vcoul, comm, ierr)

end subroutine mc_integrate
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/mc_free
!! NAME
!! mc_free
!!
!! FUNCTION
!!  Free dynamic memory
!!
!! SOURCE

subroutine mc_free(mc)

!Arguments ------------------------------------
 class(mc_t),intent(inout) :: mc

! *************************************************************************

 ABI_SFREE(mc%qran)

end subroutine mc_free
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/beigi_cylinder_limit
!! NAME
!! beigi_cylinder_limit
!!
!! FUNCTION
!!
!! SOURCE

subroutine beigi_cylinder_limit(opt_cylinder, cryst, nqibz, nkbz, rcut, hcyl, boxcenter, pdir, i_sz)

!Arguments ------------------------------------
 integer,intent(in) :: opt_cylinder, nqibz, nkbz, pdir(3)
 type(crystal_t),intent(in) :: cryst
 real(dp),intent(in) :: rcut, hcyl, boxcenter(3)
 real(dp),intent(out) :: i_sz

!Local variables-------------------------------
 integer :: ii, iq, npar, npt, gamma_pt(3,1)
 real(dp) :: step, bz_plane, dx, integ, q0_vol, q0_volsph, b1(3),b2(3),b3(3)
 real(dp),allocatable :: cov(:,:),par(:),qfit(:,:),sigma(:),var(:), vcfit(:,:),xx(:),yy(:)

! *************************************************************************

 b1 = two_pi * cryst%gprimd(:,1); b2 = two_pi * cryst%gprimd(:,2); b3 = two_pi * cryst%gprimd(:,3)

 npt =100
 npar=8; gamma_pt = RESHAPE(([0, 0, 0]), [3, 1])
 ABI_MALLOC(qfit, (3, npt))
 ABI_MALLOC(vcfit, (1, npt))
 if (nqibz == 1) then
   ABI_ERROR("nqibz == 1 not supported when Beigi's method is used")
 endif
 qfit(:,:)=zero
 step=half/(npt * (nqibz-1))              ; qfit(3,:)=arth(tol6,step,npt)
 !step=(half/(nqibz-1)/tol6)**(one/npt) ; qfit(3,:)=geop(tol6,step,npt)

 do iq=1,npt
   call cutoff_cylinder(qfit(:,iq),1,gamma_pt,rcut,hcyl,pdir,boxcenter,&
                        Cryst%rprimd,vcfit(:,iq),opt_cylinder, xmpi_comm_self)
 end do

 ABI_MALLOC(xx, (npt))
 ABI_MALLOC(yy, (npt))
 ABI_MALLOC(sigma, (npt))
 ABI_MALLOC(par, (npar))
 ABI_MALLOC(var, (npar))
 ABI_MALLOC(cov, (npar, npar))

 do ii=1,npt
   xx(ii) = normv(qfit(:,ii), cryst%gmet, 'G')
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
 !write(std_out,*)' simple integral',integ
 q0_volsph = (two_pi)**3 / (nkbz * cryst%ucvol)
 q0_vol=bz_plane*two*xx(npt)
 !write(std_out,*)' q0 sphere : ',q0_volsph,' q0_vol cyl ',q0_vol
 i_sz = bz_plane * two * integ / q0_vol
 !write(std_out,*)' spherical approximation ',four_pi*7.44*q0_volsph**(-two_thirds)
 !write(std_out,*)' Cylindrical cutoff value ',i_sz
 !i_sz=four_pi*7.44*q0_vol**(-two_thirds)

 ABI_FREE(xx)
 ABI_FREE(yy)
 ABI_FREE(sigma)
 ABI_FREE(par)
 ABI_FREE(var)
 ABI_FREE(cov)

end subroutine beigi_cylinder_limit
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/beigi_surface_limit
!! NAME
!! beigi_surface_limit
!!
!! FUNCTION
!!
!! SOURCE

subroutine beigi_surface_limit(opt_surface, cryst, nqibz, nkbz, rcut, alpha, boxcenter, pdir, i_sz)

!Arguments ------------------------------------
 integer,intent(in) :: opt_surface, nqibz, nkbz, pdir(3)
 type(crystal_t),intent(in) :: cryst
 real(dp),intent(in) :: rcut, alpha(3), boxcenter(3)
 real(dp),intent(out) :: i_sz

!Local variables-------------------------------
 integer :: ii, npt, gamma_pt(3,1)
 real(dp) :: step, bz_plane, dx, integ, q0_vol, q0_volsph, b1(3),b2(3),b3(3)
 real(dp),allocatable :: qfit(:,:),sigma(:),vcfit(:,:),xx(:),yy(:), qcart(:,:)

! *************************************************************************

 b1 = two_pi * cryst%gprimd(:,1); b2 = two_pi * cryst%gprimd(:,2); b3 = two_pi * cryst%gprimd(:,3)

 gamma_pt=RESHAPE([0, 0, 0], [3, 1]) ! Gamma point
 npt=100 ! Number of points in 1D
 ABI_MALLOC(qfit, (3, npt))
 ABI_MALLOC(qcart, (3, npt))
 ABI_MALLOC(vcfit, (1, npt))
 if (nqibz == 1) then
   ABI_ERROR("nqibz == 1 not supported when Beigi's method is used")
 endif
 qfit(:,:)=zero
 qcart(:,:)=zero
 ! Size of the third vector
 bz_plane=l2norm(b3)
 q0_volsph=(two_pi)**3 / (nkbz * cryst%ucvol)
 ! radius that gives the same volume as q0_volsph
 ! Let's assume that c is perpendicular to the plane
 ! We also assume isotropic BZ around gamma
 step=sqrt((q0_volsph/bz_plane)/pi)/npt

 !step=half/(npt*(nqibz-1))
 ! Let's take qpoints along 1 line, the vcut does depend only on the norm
 qcart(1,:) = arth(tol6,step,npt)

 do ii=1,npt
   qfit(:,ii) = MATMUL(TRANSPOSE(Cryst%rprimd),qcart(:,ii)) / (2*pi)
   call cutoff_surface(qfit(:,ii), 1, gamma_pt, cryst%gprimd, rcut, &
                       boxcenter, pdir, alpha, vcfit(:,ii), opt_surface)
 end do

 ABI_MALLOC(xx, (npt))
 ABI_MALLOC(yy, (npt))
 ABI_MALLOC(sigma, (npt))
 do ii=1,npt
   !xx(ii)=qfit(1,:)
   xx(ii) = normv(qfit(:,ii), cryst%gmet, 'G')
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
 !write(std_out,*)' simple integral',integ
 q0_vol=bz_plane*pi*xx(npt)**2
 !write(std_out,*)' q0 sphere : ',q0_volsph,' q0_vol cyl ',q0_vol
 i_sz=bz_plane*2*pi*integ/q0_vol
 !write(std_out,*)' spherical approximation ',four_pi*7.44*q0_volsph**(-two_thirds)
 !write(std_out,*)' Cylindrical cutoff value ',i_sz
 !i_sz=four_pi*7.44*q0_vol**(-two_thirds)
 ABI_FREE(xx)
 ABI_FREE(yy)

end subroutine beigi_surface_limit
!!***

!!****f* m_vcoul/carrier_isz
!! NAME
!!  carrier_isz
!!
!! FUNCTION
!!
!! SOURCE

real(dp) function carrier_isz(cryst, nqbz, qbz, rcut, comm) result(i_sz)

!Arguments ------------------------------------
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: nqbz, comm
 real(dp), intent(in) :: qbz(3, nqbz), rcut

!Local variables-------------------------------
 integer :: iq_bz
 real(dp) :: qbz_norm, bz_geometry_factor, qbz_cart(3), b1(3), b2(3), b3(3)

!************************************************************************

 b1 = two_pi * cryst%gprimd(:,1); b2 = two_pi * cryst%gprimd(:,2); b3 = two_pi * cryst%gprimd(:,3)

 bz_geometry_factor = zero
 do iq_bz=1,nqbz
   qbz_cart(:) = qbz(1,iq_bz)*b1(:) + qbz(2,iq_bz)*b2(:) + qbz(3,iq_bz)*b3(:)
   qbz_norm = SQRT(SUM(qbz_cart(:)**2))
   if (qbz_norm > TOLQ0) bz_geometry_factor = bz_geometry_factor - faux(qbz(:,iq_bz), rcut, b1, b2, b3)
 end do

 bz_geometry_factor = bz_geometry_factor + integratefaux(rcut, cryst%gprimd, cryst%ucvol, comm) * nqbz
 i_sz = four_pi * bz_geometry_factor  ! Final result stored here

end function carrier_isz
!!***

!!****f* m_vcoul/gygi_baldereschi_isz
!! NAME
!!  gygi_baldereschi_isz
!!
!! FUNCTION
!!
!! SOURCE

real(dp) function gygi_baldereschi_isz(cryst, nqbz, qbz, ecut, ng, gvec) result(i_sz)

!Arguments ------------------------------------
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: nqbz, ng
 real(dp), intent(in) :: qbz(3, nqbz), ecut
 integer,intent(in) :: gvec(3,ng)

!Local variables-------------------------------
 integer :: iq_bz, ig
 real(dp) :: bz_geometry_factor, intfauxgb, alfa, qpg2, qpg(3)

!************************************************************************

 ! the choice of alfa (the width of the gaussian) is somehow empirical
 alfa = 150.0 / ecut

 bz_geometry_factor=zero
 do iq_bz=1,nqbz
   do ig = 1,ng
     qpg(:) = qbz(:,iq_bz) + gvec(:,ig)
     qpg2 = normv(qpg, cryst%gmet, 'G')**2
     if (qpg2 > TOLQ0) bz_geometry_factor = bz_geometry_factor - EXP(-alfa*qpg2)/qpg2
   end do
 end do

 intfauxgb = cryst%ucvol/four_pi/SQRT(0.5*two_pi*alfa)
 bz_geometry_factor = bz_geometry_factor + intfauxgb * nqbz

 i_sz = four_pi*bz_geometry_factor

end function gygi_baldereschi_isz
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/vcgen_init
!! NAME
!! vcgen_init
!!
!! FUNCTION
!!
!! SOURCE

subroutine vcgen_init(vcgen, cryst, kptrlatt, nkbz, nqibz, nqbz, qbz, rcut, gw_icutcoul, vcutgeo, ecut, comm)

!Arguments ------------------------------------
 class(vcgen_t),intent(out) :: vcgen
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: kptrlatt(3,3), nkbz, nqibz, nqbz, gw_icutcoul
 real(dp),intent(in) :: qbz(3,nqbz), rcut, ecut, vcutgeo(3)
 integer,intent(in) :: comm

!Local variables-------------------------------
 integer,parameter :: istwfk1 = 1
 integer :: npw_, gvec0(3)
 real(dp) :: q0_vol, bz_geometry_factor, rcut2
 character(len=500) :: msg
 real(dp) :: vcoul0(1), q_gamma(3)
 integer,allocatable :: gvec_(:,:)

! *************************************************************************

 ! Save dimension and other useful quantities in Vcp
 vcgen%rcut      = rcut                 ! Cutoff radius for cylinder.
 vcgen%hcyl      = zero                 ! Length of finite cylinder (Rozzi's method, default is Beigi).
 vcgen%boxcenter = zero                 ! Boxcenter at the moment is supposed to be at the origin.
 vcgen%vcutgeo   = vcutgeo(:)           ! Info on the orientation and extension of the cutoff region.
 vcgen%nkbz      = nkbz

 call gw_icutcoul_to_mode(gw_icutcoul, vcgen%mode)
 q_gamma = zero
 gvec0 = 0

 select case (trim(vcgen%mode))
 case ('MINIBZ', 'MINIBZ-ERFC', 'MINIBZ-ERF')

   call vcgen%mc%init(cryst%rprimd, cryst%ucvol, cryst%gprimd, cryst%gmet, kptrlatt)

   rcut2 = vcgen%rcut**2
   call vcgen%mc%integrate(vcgen%mode, q_gamma, 1, gvec0, rcut2, nkbz, vcoul0, xmpi_comm_self)

   ! Treat the limit q --> 0.
   vcgen%i_sz = vcoul0(1)

 case ('SPHERE')

   ! A non-positive value of rcut activates the recipe of Spencer & Alavi, PRB 77, 193110 (2008) [[cite:Spencer2008]].
   if (vcgen%rcut < tol12) then
     vcgen%rcut = (cryst%ucvol * nkbz * 3.d0 / four_pi) ** third
     write(msg,'(2a,2x,f8.4,a)')ch10,' Using calculated rcut: ',vcgen%rcut,' to have same volume as the BvK crystal'
     call wrtout(std_out, msg)
   end if
   vcgen%vcutgeo = zero

   ! Treat the limit q --> 0
   ! The small cube is approximated by a sphere, while vc(q=0) = 2piR**2.
   ! if a single q-point is used, the expression for the volume is exact.
   vcgen%i_sz = two_pi * vcgen%rcut**2

 case ('CYLINDER')
   call cylinder_setup(cryst, vcgen%vcutgeo, vcgen%hcyl, vcgen%pdir, vcgen%opt_cylinder)

   ! If Beigi, treat the limit q --> 0.
   if (vcgen%opt_cylinder == 1) then
     call beigi_cylinder_limit(vcgen%opt_cylinder, cryst, nqibz, nkbz, &
                               vcgen%rcut, vcgen%hcyl, vcgen%boxcenter, vcgen%pdir, vcgen%i_sz)
   else
     ! In Rozzi's method the lim q+G --> 0 is finite.
     call cutoff_cylinder(q_gamma, 1, gvec0, vcgen%rcut, vcgen%hcyl, vcgen%pdir,&
                          vcgen%boxcenter, cryst%rprimd, vcoul0, vcgen%opt_cylinder, xmpi_comm_self)
     vcgen%i_sz = vcoul0(1)
   end if

 case ('SURFACE')
   call surface_setup(cryst, vcgen%vcutgeo, vcgen%alpha, vcgen%rcut, vcgen%pdir, vcgen%opt_surface)

   ! If Beigi, treat the limit q --> 0.
   if (vcgen%opt_surface == 1) then
     ! Integrate numerically in the plane close to 0
     call beigi_surface_limit(vcgen%opt_surface, cryst, nqibz, nkbz, vcgen%rcut, vcgen%alpha, &
                              vcgen%boxcenter, vcgen%pdir, vcgen%i_sz)
   else
     ! In Rozzi's method the lim q+G --> 0 is finite.
     call cutoff_surface(q_gamma, 1, gvec0, cryst%gprimd, vcgen%rcut, &
                         vcgen%boxcenter, vcgen%pdir, vcgen%alpha, vcoul0, vcgen%opt_surface)
     vcgen%i_sz = vcoul0(1)
   end if

 case ('CRYSTAL', 'AUXILIARY_FUNCTION', "AUX_GB")

   if (vcgen%mode == "CRYSTAL") then
     ! Analytic integration of 4pi/q^2 over the volume element:
     ! $4pi/V \int_V d^3q 1/q^2 =4pi bz_geometric_factor V^(-2/3)$
     ! i_sz=4*pi*bz_geometry_factor*q0_vol**(-two_thirds) where q0_vol= V_BZ/N_k
     ! bz_geometry_factor: sphere=7.79, fcc=7.44, sc=6.188, bcc=6.946, wz=5.255 (see gwa.pdf, appendix A.4)
     q0_vol = (two_pi) **3 / (nkbz*cryst%ucvol); bz_geometry_factor=zero
     vcgen%i_sz = four_pi*7.44*q0_vol**(-two_thirds)

   else if (vcgen%mode == "AUXILIARY_FUNCTION") then
     ! Numerical integration of the exact-exchange divergence through the
     ! auxiliary function of Carrier et al. PRB 75, 205126 (2007) [[cite:Carrier2007]].
     vcgen%i_sz = carrier_isz(cryst, nqbz, qbz, rcut, comm)

   else if (vcgen%mode == "AUX_GB") then
     ! We use the auxiliary function of a Gygi-Baldereschi variant [[cite:Gigy1986]]
     ABI_ERROR("AUX_GB not implemented in vcgen_init")
     !call get_kg(kk_bz, istwfk1, ecut, cryst%gmet, npw_, gvec_)
     !vcgen%i_sz = gygi_baldereschi_isz(cryst, nqbz, qbz, ecut, ng, gvec_)
     !ABI_FREE(gvec_)

   else
     ABI_ERROR(sjoin("Need treatment of 1/q^2 singularity! for mode", vcgen%mode))
   end if

 case ('ERF')
   vcgen%i_sz = carrier_isz(cryst, nqbz, qbz, rcut, xmpi_comm_self)

 case ('ERFC')
   ! === Treat 1/q^2 singularity ===
   ! * There is NO singularity in this case.
   vcgen%i_sz = pi * vcgen%rcut**2 ! Final result stored here

 case default
   ABI_BUG(sjoin('Unsupported cutoff mode:', vcgen%mode))
 end select

 call wrtout(std_out, sjoin("vcgen%i_sz", ftoa(vcgen%i_sz)))

end subroutine vcgen_init
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/vcgen_get_vc_sqrt
!! NAME
!! vcgen_get_vc_sqrt
!!
!! FUNCTION
!!
!! SOURCE

subroutine vcgen_get_vc_sqrt(vcgen, qpt, npw, gvec, q0, cryst, vc_sqrt, comm)

!Arguments ------------------------------------
 class(vcgen_t),intent(in) :: vcgen
 real(dp),intent(in) :: qpt(3), q0(3)
 integer,intent(in) :: npw, gvec(3,npw), comm
 type(crystal_t),intent(in) :: cryst
 complex(gwpc),intent(out) :: vc_sqrt(npw)

!Local variables-------------------------------
 integer :: ig, ig0
 real(dp) :: rcut2
 logical :: q_is_gamma
 real(dp),allocatable :: vcoul(:)

! *************************************************************************

 q_is_gamma = normv(qpt, cryst%gmet, "G") < GW_TOLQ0

 ! Find index of G=0 in gvec.
 ig0 = -1
 do ig=1,npw
   if (all(gvec(:,ig) == 0)) then
    ig0 = ig; exit
   end if
 end do
 ABI_CHECK(ig0 /= -1, "Cannot find G=0 in gvec!")

 ABI_MALLOC(vcoul, (npw))

 select case (trim(vcgen%mode))
 case ('MINIBZ', 'MINIBZ-ERFC', 'MINIBZ-ERF')
   rcut2 = vcgen%rcut**2
   call vcgen%mc%integrate(vcgen%mode, qpt, npw, gvec, rcut2, vcgen%nkbz, vcoul, comm)

 case ('SPHERE')
   call cutoff_sphere(qpt, npw, gvec, cryst%gmet, vcgen%rcut, vcoul)

 case ('CYLINDER')
   call cutoff_cylinder(qpt, npw, gvec, vcgen%rcut, vcgen%hcyl, vcgen%pdir, &
                        vcgen%boxcenter, cryst%rprimd, vcoul, vcgen%opt_cylinder, comm)

 case ('SURFACE')
   call cutoff_surface(qpt, npw, gvec, cryst%gprimd, vcgen%rcut, &
                       vcgen%boxcenter, vcgen%pdir, vcgen%alpha, vcoul, vcgen%opt_surface)

 case ('CRYSTAL', 'AUXILIARY_FUNCTION', "AUX_GB", "ERF", "ERFC")
   ! Compute |q+G| with special treatment of (q=0, g=0).
   do ig=1,npw
     !if (q_is_gamma) then
     if (q_is_gamma .and. ig == ig0) then
       vcoul(ig) = normv(q0 + gvec(:,ig), cryst%gmet, "G")
     else
       vcoul(ig) = normv(qpt + gvec(:,ig), cryst%gmet, "G")
     end if
   end do

   if (vcgen%mode == "ERF") then
     vcoul(:)  = four_pi/(vcoul(:)**2) *  EXP( -0.25d0 * (vcgen%rcut*vcoul(:))**2 )
   else if (vcgen%mode == "ERFC") then
     vcoul(:) = four_pi/(vcoul(:)**2) * ( one - EXP( -0.25d0 * (vcgen%rcut*vcoul(:))**2 ) )
   else
     vcoul = four_pi/vcoul**2
   end if

 case default
   ABI_BUG(sjoin('Unsupported cutoff mode:', vcgen%mode))
 end select

 ! Store final results in complex array as Rozzi's cutoff can give real negative values
 vc_sqrt = SQRT(CMPLX(vcoul, zero))
 ABI_FREE(vcoul)

end subroutine vcgen_get_vc_sqrt
!!***

!----------------------------------------------------------------------

!!****f* m_vcoul/vcgen_free
!! NAME
!! vcgen_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! SOURCE

subroutine vcgen_free(vcgen)

!Arguments ------------------------------------
 class(vcgen_t),intent(inout) :: vcgen

! *************************************************************************

 call vcgen%mc%free()

end subroutine vcgen_free
!!***

!----------------------------------------------------------------------

end module m_vcoul
!!***
