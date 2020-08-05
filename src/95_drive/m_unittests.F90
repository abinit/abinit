!!****m* ABINIT/m_unittests
!! NAME
!! m_unittests
!!
!! FUNCTION
!! Module to implement unit tests
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (HM)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_unittests

 use defs_basis
 use m_errors
 use m_abicore
 use m_crystal
 use m_defs_ptgroups
 use m_ptgroups
 use m_tetrahedron
 use m_htetra
 use m_krank
 use m_symkpt
 use m_sort
 use m_xmpi
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_time,            only : cwtime, cwtime_report
 use m_fstrings,        only : ltoa, itoa, sjoin, strcat
 use m_numeric_tools,   only : linspace, ctrap, simpson_int
 use m_special_funcs,   only : gaussian
 use m_symtk,           only : matr3inv
 use m_io_tools,        only : open_file
 use m_kpts,            only : kpts_ibz_from_kptrlatt, listkk
 use m_geometry,        only : normv

 implicit none

 public :: tetra_unittests
 public :: tetra_zinv_convergence
 public :: kptrank_unittests
!!***

contains

!!****f* m_unittests/rprimd_from_ptgroup
!! NAME
!!  rprimd_from_ptgroup
!!
!! FUNCTION
!!  Generate rprimd compatible with the point group
!!  This is only for testing porposes
!!
function rprimd_from_ptgroup(ptgroup) result(rprim)

 character(len=*),intent(in) :: ptgroup
 real(dp) :: a,b,c, tmp1, tmp2, alpha,beta, gamma, tx,ty,tz
 real(dp) :: rprim(3,3)

 a = one * 5
 b = two * 5
 c = three * 5
 alpha =     pi/two
 beta  =     pi/three
 gamma = two_pi/three

 select case(ptgroup)
 case ('1','-1')
   ! Triclinic
   tmp1 = c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
   tmp2 = c*sqrt( 1 + 2*cos(alpha)*cos(beta)*cos(gamma)-&
            cos(alpha)**2-cos(beta)**2-cos(gamma)**2 )/sin(gamma)
   rprim(:,1) = [a,                  0.0_dp, 0.0_dp]
   rprim(:,2) = [b*cos(gamma), b*sin(gamma), 0.0_dp]
   rprim(:,3) = [c*cos(beta),          tmp1,   tmp2]
 case ('2','m','-2','2/m')
   ! Monoclinic
   rprim(:,1) = [           a,       0.0_dp, 0.0_dp]
   rprim(:,2) = [b*cos(gamma), b*sin(gamma), 0.0_dp]
   rprim(:,3) = [      0.0_dp,       0.0_dp,      c]
 case ('222','mm2','mmm')
   ! Orthorhombic
   rprim(:,1) = [     a,0.0_dp,0.0_dp]
   rprim(:,2) = [0.0_dp,     b,0.0_dp]
   rprim(:,3) = [0.0_dp,0.0_dp,     c]
 case ('4','-4','4/m','422','4mm','-42m','4/mmm')
   ! Tetragonal
   rprim(:,1) = a*[1.0_dp, 0.0_dp, 0.0_dp]
   rprim(:,2) = a*[0.0_dp, 1.0_dp, 0.0_dp]
   rprim(:,3) = a*[0.0_dp, 0.0_dp,    c/a]
 case ('3','-3','32','3m','-3m')
   ! Trigonal
   tx = sqrt((1-c)/2)
   ty = sqrt((1-c)/6)
   tz = sqrt((1+2*c)/3)
   rprim(:,1) = a*[     tx,       -ty, tz]
   rprim(:,2) = a*[ 0.0_dp, 2.0_dp*ty, tz]
   rprim(:,3) = a*[    -tx,       -ty, tz]
 case ('6','-6','6/m','622','6mm','-62m','6/mmm')
   ! Hexagonal
   rprim(:,1) = a*[1.0_dp,        0.0_dp,  0.0_dp]
   rprim(:,2) = a*[ -half,sqrt(3.0_dp)/2,  0.0_dp]
   rprim(:,3) = a*[0.0_dp,        0.0_dp,     c/a]
 case ('23','m-3','432','-43m','m-3m')
   ! cubic
   rprim(:,1) = a*[1.0_dp, 0.0_dp, 0.0_dp]
   rprim(:,2) = a*[0.0_dp, 1.0_dp, 0.0_dp]
   rprim(:,3) = a*[0.0_dp, 0.0_dp, 1.0_dp]
 end select

end function rprimd_from_ptgroup
!!***

!----------------------------------------------------------------------

!!****f* m_unittests/crystal_from_ptgroup
!! NAME
!!  crystal_from_ptgroup
!!
!! FUNCTION
!!  Create a crystal structure from a user defined point group
!!
type(crystal_t) function crystal_from_ptgroup(ptgroup, use_symmetries) result(cryst)

!Arguments -------------------------------
 character(len=*),intent(in) :: ptgroup
 integer,intent(in) :: use_symmetries

!Variables -------------------------------
 type(irrep_t),allocatable :: irr(:)
 logical,parameter :: use_antiferro_true=.true.,remove_inv_false=.false.
 integer,parameter :: npsp1=1,space_group0=0,timrev1=1
 integer :: natom,nclass,ntypat,nsym
 integer :: typat(1)
 real(dp) :: rprimd(3,3)
 real(dp) :: amu(1),xred(3,1),znucl(1),zion(1)
 real(dp),allocatable :: tnons(:,:)
 character(len=5),allocatable :: class_names(:)
 integer,allocatable :: symrel(:,:,:),symafm(:),class_ids(:,:)

 ! Get some lattice from ptgroup
 rprimd = rprimd_from_ptgroup(ptgroup)

 if (use_symmetries == 1) then
   call get_point_group(ptgroup, nsym, nclass, symrel, class_ids, class_names, irr)
   ABI_FREE(class_ids)
   ABI_FREE(class_names)
   call irrep_free(irr)
   ABI_FREE(irr)
 else
   ! Only identity operator.
   nsym = 1
   ABI_MALLOC(symrel, (3, 3, nsym))
   symrel(:,:,1) = identity_3d
 end if

 ABI_MALLOC(symafm, (nsym))
 symafm = 1
 ABI_CALLOC(tnons, (3, nsym))

 amu = one; natom = 1; ntypat = 1; typat = 1; znucl = one; zion = one
 xred(:, 1) = zero

 call crystal_init(amu, cryst, space_group0, natom, npsp1, ntypat, nsym, rprimd, typat, xred, &
                   zion, znucl, timrev1, use_antiferro_true, remove_inv_false, "test", &
                   symrel=symrel, symafm=symafm, tnons=tnons)

 ABI_FREE(symrel)
 ABI_FREE(tnons)
 ABI_FREE(symafm)

end function crystal_from_ptgroup
!!***

!----------------------------------------------------------------------

!!****f* m_unittests/tetra_unittests
!! NAME
!!  tetra_unittests
!!
!! FUNCTION
!!  Unit tests for the tetrahedron routines
!!

subroutine tetra_unittests(ptgroup, ngqpt, use_symmetries, prtvol, comm)

!Arguments -------------------------------
!scalars
 character(len=*),intent(in) :: ptgroup
 integer,intent(in) :: use_symmetries, prtvol, comm
 integer,intent(in) :: ngqpt(3)

!Local variables -------------------------
!scalars
 integer,parameter :: qptopt1 = 1, nqshft1 = 1, bcorr0 = 0, bcorr1 = 1, master = 0
 integer :: nqibz,iq_ibz,nqbz,ierr,nw, my_rank, return_code
 real(dp),parameter :: max_occ1 = one
 real(dp) :: cpu, wall, gflops, dosdeltae, emin, emax, qnorm, int_dos, broad, min_eig, max_eig, mstar
 character(len=80) :: errstr
 logical,parameter :: use_old_tetra = .False.
 type(crystal_t) :: cryst
 type(t_tetrahedron) :: tetraq
 type(htetra_t) :: htetraq
#ifdef HAVE_NETCDF
 integer :: ncid, ncerr
#endif
!arrays
 integer :: in_qptrlatt(3,3), new_qptrlatt(3,3)
 integer,allocatable :: bz2ibz(:,:)
 real(dp) :: qshift(3,nqshft1), rlatt(3,3), qlatt(3,3)
 real(dp),allocatable :: tweight(:,:), dweight(:,:)
 real(dp),allocatable :: qbz(:,:), qibz(:,:), wtq_ibz(:), wdt(:,:)
 real(dp),allocatable :: wmesh(:), eig(:),mat(:), dos(:), idos(:), cauchy_ppart(:)
 complex(dp),allocatable :: zmesh(:), cweight(:,:)

! *********************************************************************

 my_rank = xmpi_comm_rank(comm)
 return_code = 0

 call wrtout(std_out, sjoin(" Tetrahedron unit tests with ptgroup:", ptgroup, ", ngqpt", ltoa(ngqpt)))
 call cwtime(cpu, wall, gflops, "start")

 ! 0. Initialize
 call wrtout(std_out, ' [0] Initializing tetrahedron object...')

 ! Create fake crystal from ptgroup
 cryst = crystal_from_ptgroup(ptgroup, use_symmetries)

 ! Create a regular grid
 in_qptrlatt = 0; in_qptrlatt(1,1) = ngqpt(1); in_qptrlatt(2,2) = ngqpt(2); in_qptrlatt(3,3) = ngqpt(3)
 qshift = zero

 call kpts_ibz_from_kptrlatt(cryst, in_qptrlatt, qptopt1, nqshft1, qshift, &
                             nqibz, qibz, wtq_ibz, nqbz, qbz, new_kptrlatt=new_qptrlatt, bz2ibz=bz2ibz)
 if (prtvol > 0) call cwtime_report(" kpts_ibz_from_kptrlatt", cpu, wall, gflops)

 rlatt = new_qptrlatt; call matr3inv(rlatt, qlatt)

 if (use_old_tetra) then
   ! Initialize old tetrahedra
   call init_tetra(bz2ibz(1,:), cryst%gprimd, qlatt, qbz, nqbz, tetraq, ierr, errstr, comm)
   if (prtvol > 0) call cwtime_report(" init_tetra_old ", cpu, wall, gflops)
 end if

 ! Initialize new tetrahedra
 call htetra_init(htetraq, bz2ibz(1,:), cryst%gprimd, qlatt, qbz, nqbz, qibz, nqibz, ierr, errstr, comm)
 if (prtvol > 0) then
   call htetraq%print(std_out)
   call cwtime_report(" init_htetra", cpu, wall, gflops)
 end if

 ! 1. Compute parabolic band
 call wrtout(std_out, ' [1] Testing tetrahedra with parabolic band dispersion ... ', newlines=1)
 ABI_MALLOC(eig, (nqibz))
 ABI_MALLOC(mat, (nqibz))
 mstar = one
 mstar = ten
 !mstar = 0.05_dp
 !mstar = 0.5_dp

 do iq_ibz=1,nqibz
   qnorm = normv(qibz(:, iq_ibz), cryst%gmet, 'G')
   ! The DOS for this function goes as sqrt(w)
   eig(iq_ibz) = half * qnorm ** 2 / mstar
   !eig(iq_ibz) = cos(twp * qnorm / q0)
   mat(iq_ibz) = one
   !mat(iq_ibz) = abs(one / eig(iq_ibz))
 end do

 ! Prepare DOS calculation
 !emin = minval(eig) - one; emax = maxval(eig) + one
 emin = minval(eig); emax = maxval(eig)
 nw = 300
 dosdeltae = (emax - emin) / (nw - 1)
 !broad = tol2 * eV_Ha
 broad = tol1 * eV_Ha
 min_eig = minval(eig); max_eig = maxval(eig)

 ABI_MALLOC(wmesh, (nw))
 wmesh = linspace(emin, emax, nw)

 if (my_rank == master) then
   write(std_out, *)" min, Max band energy: ", min_eig, max_eig
   write(std_out, *)" energy mesh, Max: ", emin, emax, nw
   write(std_out, *)" Broad: ", broad
   write(std_out, *)" Effective mass: ", mstar
   write(std_out, *)" Use_symmetries: ", use_symmetries
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_create(ncid, "foo_TETRATEST.nc", xmpi_comm_self))
   NCF_CHECK(cryst%ncwrite(ncid))
   ! Add dimensions.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nw", nw), nctkdim_t("nqibz", nqibz) ], defmode=.True.)
   NCF_CHECK(ncerr)

   ! Define arrays with results.
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("eig", "dp", "nqibz"), &
     nctkarr_t("wmesh", "dp", "nw") &
   ])
   NCF_CHECK(ncerr)

   !ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
   !  "eph_task", "symsigma", "nbsum", "bsum_start", "bsum_stop", "symdynmat", &
   !  "imag_only", "symv1scf", "dvdb_add_lr", "mrta"])
   !NCF_CHECK(ncerr)
   !ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
   !  "eta", "wr_step", "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie", "ph_wstep", "ph_smear"])
   !NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eig"), eig))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "wmesh"), wmesh))
#endif
 end if

 ABI_MALLOC(dos, (nw))
 ABI_MALLOC(idos, (nw))
 ABI_MALLOC(cauchy_ppart, (nw))
 ABI_MALLOC(wdt, (nw, 2))

 if (prtvol > 0) call cwtime_report(" init", cpu, wall, gflops)

 ABI_CALLOC(tweight, (nw, nqibz))
 ABI_CALLOC(dweight, (nw, nqibz))

 if (use_old_tetra) then
   ! Compute tetra weights and DOS/IDOS with old implementation
   call tetra_blochl_weights(tetraq, eig, emin, emax, max_occ1, nw, nqibz, bcorr0, tweight, dweight, comm)
   do iq_ibz=1,nqibz
     dweight(:,iq_ibz) = dweight(:,iq_ibz) * mat(iq_ibz)
     tweight(:,iq_ibz) = tweight(:,iq_ibz) * mat(iq_ibz)
   end do
   dos(:)  = sum(dweight, dim=2)
   idos(:) = sum(tweight, dim=2)
   if (my_rank == master) call write_file("parabola", "old_tetra", nw, wmesh, dos, idos=idos)
 end if

 ! Compute DOS/IDOS with finite broadening
 dos = zero; cauchy_ppart = zero
 do iq_ibz=1,nqibz
   dos = dos - wtq_ibz(iq_ibz) * aimag(one / (wmesh - eig(iq_ibz) + j_dpc * broad)) / pi
   cauchy_ppart = cauchy_ppart + wtq_ibz(iq_ibz) * real(one / (wmesh - eig(iq_ibz) + j_dpc * broad))
 end do
 if (my_rank == master) call write_file("parabola", "broad", nw, wmesh, dos, cauchy_ppart=cauchy_ppart)

 ! Compute tetra weights and DOS/IDOS with new implementation by Henrique
 call htetraq%blochl_weights(eig, emin, emax, max_occ1, nw, nqibz, bcorr0, tweight, dweight, comm)

 do iq_ibz=1,nqibz
   dweight(:,iq_ibz) = dweight(:,iq_ibz) * mat(iq_ibz)
   tweight(:,iq_ibz) = tweight(:,iq_ibz) * mat(iq_ibz)
 end do

 dos(:)  = sum(dweight, dim=2)
 idos(:) = sum(tweight, dim=2)

 if (my_rank == master) call write_file("parabola", "htetra_blochl", nw, wmesh, dos, idos=idos)

 ! Compute tetra weights with Blochl corrections and DOS/IDOS with new implementation by Henrique
 call htetraq%blochl_weights(eig, emin, emax, max_occ1, nw, nqibz, bcorr1, tweight, dweight, comm)

 do iq_ibz=1,nqibz
   dweight(:,iq_ibz) = dweight(:,iq_ibz) * mat(iq_ibz)
   tweight(:,iq_ibz) = tweight(:,iq_ibz) * mat(iq_ibz)
 end do
 dos(:)  = sum(dweight, dim=2)
 idos(:) = sum(tweight, dim=2)
 if (my_rank == master) call write_file("parabola", "htetra_blochl_corr", nw, wmesh, dos, idos=idos)

 ! Compute weights using LV integration from TDEP
 call htetraq%blochl_weights(eig, emin, emax, max_occ1, nw, nqibz, 2, tweight, dweight, comm)

 do iq_ibz=1,nqibz
   dweight(:,iq_ibz) = dweight(:,iq_ibz) * mat(iq_ibz)
   tweight(:,iq_ibz) = tweight(:,iq_ibz) * mat(iq_ibz)
 end do
 dos(:)  = sum(dweight, dim=2)
 idos(:) = sum(tweight, dim=2)
 if (my_rank == master) call write_file("parabola", "htetra_lv", nw, wmesh, dos, idos=idos)

 ABI_MALLOC(zmesh, (nw))
 ABI_MALLOC(cweight, (nw, nqibz))
 zmesh = wmesh + j_dpc * broad

 ! Compute weights using SIMTET routines and erange
 call htetraq%weights_wvals_zinv(eig, nw, zmesh, max_occ1, nqibz, 1, cweight, comm, &
                                 erange=[min_eig * (one - tol2), max_eig * (one + tol2)])

 dos(:)  = -sum(aimag(cweight(:,:)), dim=2) / pi
 cauchy_ppart =  sum(real(cweight(:,:)), dim=2)

 if (my_rank == master) call write_file("parabola", "simtet_erange", nw, wmesh, dos, cauchy_ppart=cauchy_ppart)

 ! Compute weights using SIMTET routines (slow but accurate)
 call htetraq%weights_wvals_zinv(eig, nw, zmesh, max_occ1, nqibz, 1, cweight, comm)

 dos(:)  = -sum(aimag(cweight(:,:)), dim=2) / pi
 cauchy_ppart =  sum(real(cweight(:,:)), dim=2)
 if (my_rank == master) call write_file("parabola", "simtet", nw, wmesh, dos, cauchy_ppart=cauchy_ppart)

 ! Use LV integration from TDEP
 call htetraq%weights_wvals_zinv(eig, nw, zmesh, max_occ1, nqibz, 2, cweight, comm)

 dos(:)  = -sum(aimag(cweight(:,:)), dim=2) / pi
 cauchy_ppart =  sum(real(cweight(:,:)), dim=2)
 if (my_rank == master) call write_file("parabola", "zinv_lv", nw, wmesh, dos, cauchy_ppart=cauchy_ppart)

 dos = zero; idos = zero
 do iq_ibz=1,nqibz
   call htetraq%get_onewk_wvals(iq_ibz, bcorr0, nw, wmesh, max_occ1, nqibz, eig, wdt)
   wdt(:,:) = wdt(:,:) * mat(iq_ibz)
   dos(:)  = dos(:)  + wdt(:,1) * wtq_ibz(iq_ibz)
   idos(:) = idos(:) + wdt(:,2) * wtq_ibz(iq_ibz)
 end do
 if (my_rank == master) call write_file("parabola", "htetra_onewk_wvals", nw, wmesh, dos, idos=idos)

 dos = zero; idos = zero
 do iq_ibz=1,nqibz
   call htetraq%get_onewk(iq_ibz, bcorr0, nw, nqibz, eig, emin, emax, max_occ1, wdt)
   wdt(:,:) = wdt(:,:) * mat(iq_ibz)
   dos(:)  = dos(:)  + wdt(:, 1) * wtq_ibz(iq_ibz)
   idos(:) = idos(:) + wdt(:, 2) * wtq_ibz(iq_ibz)
 end do
 if (my_rank == master) call write_file("parabola", "htetra_onewk", nw, wmesh, dos, idos=idos)

 ! 2. Compute energies for a flat band
 call wrtout(std_out, " [2] Testing tetrahedra with flat band at zero ...", pre_newlines=1)
 eig = zero

 if (use_old_tetra) then
   ! Compute DOS using old tetrahedron implementation
   call tetra_blochl_weights(tetraq, eig, emin, emax, max_occ1, nw, nqibz, bcorr0, tweight, dweight, comm)
   dos(:)  = sum(dweight, dim=2)
   idos(:) = sum(tweight, dim=2)
   if (my_rank == master) call write_file("flat", "tetra_old", nw, wmesh, dos, idos=idos)
 end if

 ! Compute DOS using new tetrahedron implementation from Henrique
 call htetraq%blochl_weights(eig, emin, emax, max_occ1, nw, nqibz, bcorr0, tweight, dweight, comm)

 dos(:)  = sum(dweight, dim=2)
 idos(:) = sum(tweight, dim=2)
 if (my_rank == master) call write_file("flat", "htetra", nw, wmesh, dos, idos=idos)

 dos = zero; idos = zero
 do iq_ibz=1,nqibz
   call htetraq%get_onewk_wvals(iq_ibz, bcorr0, nw, wmesh, max_occ1, nqibz, eig, wdt)
   dos(:)  = dos(:)  + wdt(:, 1) * wtq_ibz(iq_ibz)
   idos(:) = idos(:) + wdt(:, 2) * wtq_ibz(iq_ibz)
 end do
 if (my_rank == master) call write_file("flat", "htetra_onewk_wvals", nw, wmesh, dos, idos=idos)

 dos = zero; idos = zero
 do iq_ibz=1,nqibz
   call htetraq%get_onewk(iq_ibz, bcorr0, nw, nqibz, eig, emin, emax, max_occ1, wdt)
   dos(:)  = dos(:)  + wdt(:, 1) * wtq_ibz(iq_ibz)
   idos(:) = idos(:) + wdt(:, 2) * wtq_ibz(iq_ibz)
 end do
 if (my_rank == master) call write_file("flat", "htetra_onewk", nw, wmesh, dos, idos=idos)

 ! Compute weights using SIMTET routines (slow but accurate)
 call htetraq%weights_wvals_zinv(eig, nw, zmesh, max_occ1, nqibz, 1, cweight, comm)

 dos(:)  = -sum(aimag(cweight(:,:)), dim=2) / pi
 cauchy_ppart =  sum(real(cweight(:,:)), dim=2)
 if (my_rank == master) call write_file("flat", "simtet", nw, wmesh, dos, cauchy_ppart=cauchy_ppart)

 ! Use LV integration from TDEP
 call htetraq%weights_wvals_zinv(eig, nw, zmesh, max_occ1, nqibz, 2, cweight, comm)

 dos(:)  = -sum(aimag(cweight(:,:)), dim=2) / pi
 cauchy_ppart =  sum(real(cweight(:,:)), dim=2)
 if (my_rank == master) call write_file("flat", "zinv_lv", nw, wmesh, dos, cauchy_ppart=cauchy_ppart)

 ! 3. Compute tetrahedron for simple TB FCC
 !TODO

 ! Free memory
 ABI_SFREE(wmesh)
 ABI_SFREE(eig)
 ABI_SFREE(mat)
 ABI_SFREE(wdt)
 ABI_SFREE(tweight)
 ABI_SFREE(dweight)
 ABI_SFREE(wtq_ibz)
 ABI_SFREE(dos)
 ABI_SFREE(idos)
 ABI_SFREE(cauchy_ppart)
 ABI_SFREE(qbz)
 ABI_SFREE(qibz)
 ABI_SFREE(bz2ibz)
 ABI_SFREE(zmesh)
 ABI_SFREE(cweight)

 call cryst%free()
 call htetraq%free()
 call destroy_tetra(tetraq)

 if (my_rank == master) then
#ifdef HAVE_NETCDF
   NCF_CHECK(nf90_close(ncid))
#endif
 end if

 !if (return_code /= 0) then
 !  MSG_ERROR("Some of the tetrahedron unit tests failed. See above results")
 !end if

 contains

 subroutine write_file(ekind, algo, nw, wmesh, dos, idos, cauchy_ppart)

  integer,intent(in) :: nw
  character(len=*), intent(in) :: ekind, algo
  real(dp),intent(in) :: wmesh(nw), dos(nw)
  real(dp),optional, intent(in) :: idos(nw)
  real(dp),optional,intent(in) :: cauchy_ppart(nw)

!Local variables -------------------------
  character(len=500) :: msg
  character(len=fnlen) :: fname
  character(len=nctk_slen) :: dos_vname, idos_vname, ppart_vname, key
  integer :: iw,funit
  real(dp) :: my_idos(nw), rerr
! *********************************************************************

  key = strcat(ekind, "_", algo)
  if (prtvol > 0) call cwtime_report("- "//trim(key), cpu, wall, gflops)

  if (.not. present(idos)) then
    call simpson_int(nw, dosdeltae, dos, my_idos)
  else
    my_idos = idos
  end if

  call ctrap(nw, dos, dosdeltae, int_dos)

  rerr = 100 * (int_dos - my_idos(nw)) / my_idos(nw)
  msg = "[OK]:";
  if (abs(rerr) > half .or. abs(my_idos(nw) - one) > two * tol2) then
    msg = "[FAILED]:"
    return_code = return_code + 1
  end if
  write(std_out, "(1x,2(a,1x),/,4x,a,3(f10.5,1x),/)") &
    trim(key), trim(msg), " integral_dos, idos(nw), relative_err: ", int_dos, my_idos(nw), rerr

#ifdef HAVE_NETCDF
  dos_vname = strcat("dos_", key)
  idos_vname = strcat("idos_", key)
  ppart_vname = strcat("cauchy_ppart_", key)

  ! Define arrays with results.
  ncerr = nctk_def_arrays(ncid, [ &
    nctkarr_t(dos_vname, "dp", "nw"), &
    nctkarr_t(idos_vname, "dp", "nw") &
  ])
  NCF_CHECK(ncerr)
  if (present(cauchy_ppart)) then
    NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t(ppart_vname, "dp", "nw")]))
  end if
  NCF_CHECK(nctk_set_datamode(ncid))
  NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, dos_vname), dos))
  NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, idos_vname), my_idos))
  if (present(cauchy_ppart)) then
    NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, ppart_vname), cauchy_ppart))
  end if
#endif

  ! Write results to txt file.
  fname = trim(key)//".dat"
  if (open_file(fname, msg, newunit=funit, action="write") /= 0) then
    MSG_ERROR(msg)
  end if

  if (.not. present(cauchy_ppart)) then
    write(funit, *)"# Energies, DOS, IDOS"
    do iw=1,nw
      write(funit,*) wmesh(iw), dos(iw), my_idos(iw)
    end do
  else
    write(funit, *)"# Energies DOS, IDOS, Cauchy_PrincipalPart"
    do iw=1,nw
      write(funit,*) wmesh(iw), dos(iw), my_idos(iw), cauchy_ppart(iw)
    end do
  end if

  close(funit)

 end subroutine write_file

end subroutine tetra_unittests
!!***

!----------------------------------------------------------------------

!!****f* m_unittests/tetra_zinv_convergence
!! NAME
!!  tetra_zinv_convergence
!!
!! FUNCTION
!!
subroutine tetra_zinv_convergence(ptgroup, use_symmetries, comm)

!Arguments -------------------------------
!scalars
 character(len=*),intent(in) :: ptgroup
 integer,intent(in) :: use_symmetries, comm

!Local variables -------------------------
!scalars
 integer,parameter :: qptopt1 = 1, nqshft1 = 1, master = 0
 integer :: nqibz, iq_ibz, nqbz, ierr, nw, my_rank
 integer :: num_broad, num_meshes, iq_mesh, ibroad
#ifdef HAVE_NETCDF
 integer :: ncid, ncerr
#endif
 real(dp),parameter :: max_occ1 = one
 real(dp) :: cpu, wall, gflops, dosdeltae, emin, emax, qnorm, int_dos, broad, min_eig, max_eig
 character(len=80) :: errstr
 type(crystal_t) :: cryst
 type(htetra_t) :: htetraq
!arrays
 integer :: in_qptrlatt(3,3), new_qptrlatt(3,3), ngqpt(3)
 integer,allocatable :: bz2ibz(:,:), ngqpt_list(:,:)
 real(dp) :: qshift(3,nqshft1), rlatt(3,3), qlatt(3,3)
 real(dp),allocatable :: qbz(:,:), qibz(:,:), wtq_ibz(:), broad_list(:)
 real(dp),allocatable :: wmesh(:), eig(:), mat(:), dos(:), idos(:), cauchy_ppart(:)
 complex(dp),allocatable :: zmesh(:), cweight(:,:)

! *********************************************************************

 my_rank = xmpi_comm_rank(comm)

 call wrtout(std_out, sjoin(" Tetrahedron unit tests with ptgroup:", ptgroup, ", ngqpt", ltoa(ngqpt)))
 call cwtime(cpu, wall, gflops, "start")
 !
 ! 0. Initialize
 call wrtout(std_out, ' 0. Initialize')

 ! Create fake crystal from ptgroup
 cryst = crystal_from_ptgroup(ptgroup, use_symmetries)

 num_broad = 2
 ABI_MALLOC(broad_list, (num_broad))
 broad_list = [tol1, tol2] * eV_Ha

 num_meshes = 3
 ABI_MALLOC(ngqpt_list, (3, num_meshes))
 ngqpt_list(:, 1) = 30
 ngqpt_list(:, 2) = 60
 ngqpt_list(:, 3) = 90

 do iq_mesh=1,num_meshes
   ngqpt = ngqpt_list(:, iq_mesh)

   ! Create a regular grid
   in_qptrlatt = 0; in_qptrlatt(1,1) = ngqpt(1); in_qptrlatt(2,2) = ngqpt(2); in_qptrlatt(3,3) = ngqpt(3)
   qshift = zero

   call kpts_ibz_from_kptrlatt(cryst, in_qptrlatt, qptopt1, nqshft1, qshift, &
                               nqibz, qibz, wtq_ibz, nqbz, qbz, new_kptrlatt=new_qptrlatt, bz2ibz=bz2ibz)
   call cwtime_report(" kpts_ibz_from_kptrlatt", cpu, wall, gflops)

   rlatt = new_qptrlatt; call matr3inv(rlatt, qlatt)

   ! Initialize new tetrahedra
   call htetra_init(htetraq, bz2ibz(1,:), cryst%gprimd, qlatt, qbz, nqbz, qibz, nqibz, ierr, errstr, comm)
   call htetraq%print(std_out)
   call cwtime_report(" init_htetra", cpu, wall, gflops)

   ! 1. Compute parabolic band
   !call wrtout(std_out, ' 1. Begin testing parabolic band dispersion ... ', newlines=1)
   ABI_MALLOC(eig, (nqibz))
   ABI_MALLOC(mat, (nqibz))
   do iq_ibz=1,nqibz
     qnorm = normv(qibz(:, iq_ibz), cryst%gmet, 'G')
     ! The DOS for this function goes as sqrt(w)
     eig(iq_ibz) = half * qnorm ** 2
     !eig(iq_ibz) = cos(qnorm)
     mat(iq_ibz) = one
     !mat(iq_ibz) = abs(one / eig(iq_ibz))
   end do

   ! Prepare DOS calculation
   emin = minval(eig) - one; emax = maxval(eig) + one
   nw = 200
   dosdeltae = (emax - emin) / (nw - 1)
   min_eig = minval(eig); max_eig = maxval(eig)

   if (iq_mesh == 1) then
     ! Use same e-mesh for all q-samplings.
     ABI_MALLOC(dos, (nw))
     ABI_MALLOC(idos, (nw))
     ABI_MALLOC(cauchy_ppart, (nw))
     ABI_MALLOC(wmesh, (nw))
     wmesh = linspace(emin, emax, nw)
     ABI_MALLOC(zmesh, (nw))

     if (my_rank == master) then
       write(std_out, *)" min, Max band energy: ", min_eig, max_eig
       write(std_out, *)" energy mesh, Max: ", emin, emax, nw
       !write(std_out, *)" Broad: ", broad
#ifdef HAVE_NETCDF
       NCF_CHECK(nctk_open_create(ncid, "foo_ZINVCONV.nc", xmpi_comm_self))
       NCF_CHECK(cryst%ncwrite(ncid))

       ! Add dimensions.
       ncerr = nctk_def_dims(ncid, [ &
         nctkdim_t("nw", nw), nctkdim_t("num_meshes", num_meshes), nctkdim_t("num_broad", num_broad) &
         ], defmode=.True.)
       NCF_CHECK(ncerr)

       ! Define arrays with results.
       ncerr = nctk_def_arrays(ncid, [ &
         nctkarr_t("wmesh", "dp", "nw"), &
         nctkarr_t("dos_simple", "dp", "nw, num_broad, num_meshes"), &
         nctkarr_t("cauchy_ppart_simple", "dp", "nw, num_broad, num_meshes"), &
         nctkarr_t("dos_simtet", "dp", "nw, num_broad, num_meshes"), &
         nctkarr_t("cauchy_ppart_simtet", "dp", "nw, num_broad, num_meshes"), &
         nctkarr_t("dos_simtet_erange", "dp", "nw, num_broad, num_meshes"), &
         nctkarr_t("cauchy_ppart_simtet_erange", "dp", "nw, num_broad, num_meshes"), &
         nctkarr_t("ngqpt_list", "int", "three, num_meshes"), &
         nctkarr_t("broad_list", "dp", "num_broad") &
       ])
       NCF_CHECK(ncerr)
       NCF_CHECK(nctk_set_datamode(ncid))
       NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "wmesh"), wmesh))
       NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngqpt_list"), ngqpt_list))
       NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "broad_list"), broad_list))
#endif
     end if
   end if

   ABI_MALLOC(cweight, (nw, nqibz))

   do ibroad=1,num_broad
     broad = broad_list(ibroad)

     call cwtime_report(" init", cpu, wall, gflops)

     dos = zero; cauchy_ppart = zero
     do iq_ibz=1,nqibz
       dos = dos - wtq_ibz(iq_ibz) * aimag(one / (wmesh - eig(iq_ibz) + j_dpc * broad)) / pi
       cauchy_ppart = cauchy_ppart + wtq_ibz(iq_ibz) * real(one / (wmesh - eig(iq_ibz) + j_dpc * broad))
     end do
     call simpson_int(nw, dosdeltae, dos, idos)
     call ctrap(nw, dos, dosdeltae, int_dos)
     call cwtime_report(" finite broadening   ", cpu, wall, gflops)

     if (my_rank == master) then
       write(std_out, "(a, 3(f10.5,1x))") " int_dos, idos, rerr: ", int_dos, idos(nw), 100 * (int_dos - idos(nw)) / idos(nw)
       !call write_file('parabola_gauss.dat', nw, wmesh, dos, idos, cauchy_ppart=cauchy_ppart)
#ifdef HAVE_NETCDF
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "dos_simple"), dos, start=[1, ibroad, iq_mesh])
       NCF_CHECK(ncerr)
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "cauchy_ppart_simple"), cauchy_ppart, start=[1, ibroad, iq_mesh])
       NCF_CHECK(ncerr)
#endif
     end if

     zmesh = wmesh + j_dpc * broad

     ! Compute weights using SIMTET routines and erange
     call htetraq%weights_wvals_zinv(eig, nw, zmesh, max_occ1, nqibz, 1, cweight, comm, &
                                     erange=[min_eig * (one - tol2), max_eig * (one + tol2)])

     dos(:)  = -sum(aimag(cweight(:,:)), dim=2) / pi
     call simpson_int(nw, dosdeltae, dos, idos)
     cauchy_ppart =  sum(real(cweight(:,:)), dim=2)
     call ctrap(nw, dos, dosdeltae, int_dos)
     call cwtime_report(" htetra_weights_wvals_zinv_simtet_erange   ", cpu, wall, gflops)

     ! Compute weights using SIMTET routines (slow but accurate)
     if (my_rank == master) then
       write(std_out, "(a, 3(f10.5,1x))") " int_dos, idos, rerr: ", int_dos, idos(nw), 100 * (int_dos - idos(nw)) / idos(nw)
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "dos_simtet_erange"), dos, start=[1, ibroad, iq_mesh])
       NCF_CHECK(ncerr)
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "cauchy_ppart_simtet_erange"), cauchy_ppart, start=[1, ibroad, iq_mesh])
       NCF_CHECK(ncerr)
     end if

     call htetraq%weights_wvals_zinv(eig, nw, zmesh, max_occ1, nqibz, 1, cweight, comm)

     dos(:)  = -sum(aimag(cweight(:,:)), dim=2) / pi
     call simpson_int(nw, dosdeltae, dos, idos)
     cauchy_ppart =  sum(real(cweight(:,:)), dim=2)
     call ctrap(nw, dos, dosdeltae, int_dos)
     call cwtime_report(" htetra_weights_wvals_zinv_simtet   ", cpu, wall, gflops)

     if (my_rank == master) then
       write(std_out, "(a, 3(f10.5,1x))") " int_dos, idos, rerr: ", int_dos, idos(nw), 100 * (int_dos - idos(nw)) / idos(nw)
#ifdef HAVE_NETCDF
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "dos_simtet"), dos, start=[1, ibroad, iq_mesh])
       NCF_CHECK(ncerr)
       ncerr = nf90_put_var(ncid, nctk_idname(ncid, "cauchy_ppart_simtet"), cauchy_ppart, start=[1, ibroad, iq_mesh])
       NCF_CHECK(ncerr)
#endif
     end if

     ! Use LV integration from TDEP
     call htetraq%weights_wvals_zinv(eig, nw, zmesh, max_occ1, nqibz, 2, cweight, comm)
     dos(:)  = -sum(aimag(cweight(:,:)), dim=2) / pi
     call simpson_int(nw, dosdeltae, dos, idos)
     cauchy_ppart =  sum(real(cweight(:,:)), dim=2)
     call ctrap(nw, dos, dosdeltae, int_dos)
     call cwtime_report(" htetra_weights_wvals_zinv_lv   ", cpu, wall, gflops)

     if (my_rank == master) then
       write(std_out, "(a, 3(f10.5,1x))") " int_dos, idos, rerr: ", int_dos, idos(nw), 100 * (int_dos - idos(nw)) / idos(nw)
     end if

   end do ! ibroad

   ! Free memory
   ABI_SFREE(eig)
   ABI_SFREE(mat)
   ABI_SFREE(wtq_ibz)
   ABI_SFREE(qbz)
   ABI_SFREE(qibz)
   ABI_SFREE(bz2ibz)
   ABI_SFREE(cweight)
   call htetraq%free()
 end do ! iq_mesh

 ABI_SFREE(wmesh)
 ABI_SFREE(zmesh)
 ABI_SFREE(dos)
 ABI_SFREE(idos)
 ABI_SFREE(cauchy_ppart)
 ABI_FREE(broad_list)
 ABI_FREE(ngqpt_list)

 call cryst%free()

 if (my_rank == master) then
#ifdef HAVE_NETCDF
   NCF_CHECK(nf90_close(ncid))
#endif
 end if

end subroutine tetra_zinv_convergence
!!***

!!****f* m_unittests/kptrank_unittests
!! NAME
!!  kptrank_unittests
!!
!! FUNCTION
!!  Test the krank routines
!!
subroutine kptrank_unittests(ptgroup, ngqpt, use_symmetries, comm)

!Arguments -------------------------------
!scalars
 character(len=*),intent(in) :: ptgroup
 integer,intent(in) :: use_symmetries, comm
 integer,intent(in) :: ngqpt(3)

!Local variables -------------------------
!scalars
 integer,parameter :: qptopt1 = 1, nqshft1 = 1, iout0 = 0,chksymbreak0 = 0, sppoldbl1 = 1, master = 0
 integer :: nqibz, iqbz, iq_ibz, iqbz_rank, nqbz, nqibz_symkpt, nqibz_symkpt_new, my_rank, ierr
 integer :: in_qptrlatt(3,3), new_qptrlatt(3,3)
 real(dp) :: cpu, gflops, wall, dksqmax
 real(dp) :: qshift(3, nqshft1)
 character(len=500) :: msg
 type(crystal_t) :: cryst
 type(krank_t) :: krank
!arrays
 integer,allocatable :: bz2ibz(:,:), bz2ibz_symkpt(:,:), bz2ibz_symkpt_new(:,:)
 integer,allocatable :: bz2ibz_listkk(:,:), ibz2bz(:), ibz2bz_new(:)
 real(dp),allocatable :: wtq_fullbz(:), wtq_folded(:), wtq_ibz(:), qbz(:,:),qibz(:,:)

! *********************************************************************

 call wrtout(std_out, sjoin(" kptrank_unittests with ptgroup:", ptgroup, ", and ngqpt:", ltoa(ngqpt)))

 my_rank = xmpi_comm_rank(comm)

 ! Create fake crystal from ptgroup
 cryst = crystal_from_ptgroup(ptgroup, use_symmetries)

 ! Create a regular grid
 in_qptrlatt = 0; in_qptrlatt(1,1) = ngqpt(1); in_qptrlatt(2,2) = ngqpt(2); in_qptrlatt(3,3) = ngqpt(3)
 qshift(:,1) = 0

 call cwtime(cpu, wall, gflops, "start")
 call kpts_ibz_from_kptrlatt(cryst, in_qptrlatt, qptopt1, nqshft1, qshift, &
                             nqibz, qibz, wtq_ibz, nqbz, qbz, new_kptrlatt=new_qptrlatt, bz2ibz=bz2ibz)
 call cwtime_report(" kpts_ibz_from_kptrlatt", cpu, wall, gflops)

 ! Test krank object.
 krank = krank_new(nqbz, qbz)
 do iqbz=1,nqbz
   iqbz_rank = krank%get_index(qbz(:,iqbz))
   ABI_CHECK(iqbz == iqbz_rank, 'wrong q-point')
 end do
 call cwtime_report(" krank basic check", cpu, wall, gflops)

 ABI_MALLOC(wtq_fullbz, (nqbz))
 ABI_MALLOC(wtq_folded, (nqbz))
 ABI_MALLOC(ibz2bz, (nqbz))
 ABI_MALLOC(ibz2bz_new, (nqbz))
 ABI_MALLOC(bz2ibz_symkpt, (6, nqbz))
 ABI_MALLOC(bz2ibz_symkpt_new, (6, nqbz))
 wtq_fullbz = one / nqbz

 ! Test symkpt (note that the above call to kpts_ibz_from_kptrlatt already involves calling this routine)
 call symkpt(chksymbreak0, cryst%gmet, ibz2bz, iout0, qbz, nqbz, &
             nqibz_symkpt, cryst%nsym, cryst%symrec, cryst%timrev, &
             wtq_fullbz, wtq_folded, bz2ibz_symkpt, comm)
 call cwtime_report(" symkpt", cpu, wall, gflops)

 wtq_fullbz = one / nqbz
 call symkpt_new(chksymbreak0, cryst%gmet, ibz2bz_new, iout0, qbz, nqbz, &
                 nqibz_symkpt_new, cryst%nsym, cryst%symrec, cryst%timrev, &
                 bz2ibz_symkpt_new, comm)
 call cwtime_report(" symkpt_new", cpu, wall, gflops)
 ABI_CHECK(nqibz_symkpt == nqibz_symkpt_new, 'Wrong number of qpoints in the IBZ')

 ! Check if ibz is the same
 do iq_ibz=1,nqibz
   if (ibz2bz(iq_ibz) == ibz2bz_new(iq_ibz)) cycle
   MSG_ERROR("The IBZ is different.")
 end do

 ! Check if mapping is the same
 do iqbz=1,nqbz
   if (bz2ibz_symkpt(1, iqbz) == bz2ibz_symkpt_new(1, iqbz)) cycle
   write(msg,*) "Inconsistent mapping:", iqbz, bz2ibz_symkpt(1,iqbz), bz2ibz_symkpt_new(1,iqbz)
   MSG_ERROR(msg)
 end do

 ! Call listkk
 ABI_MALLOC(bz2ibz_listkk,(nqbz, 6))
 call listkk(dksqmax, cryst%gmet, bz2ibz_listkk, qibz, qbz, nqibz, nqbz, cryst%nsym,&
             sppoldbl1, cryst%symafm, cryst%symrec, cryst%timrev, comm, exit_loop=.True., use_symrec=.True.)
 call cwtime_report(" listkk", cpu, wall, gflops)

 ! Check if indkk is the same
 do iqbz=1,nqbz
   if (bz2ibz_listkk(iqbz, 1) /= bz2ibz_symkpt_new(1, iqbz)) then
     ierr = ierr + 1
     if (my_rank == master) write(std_out,*) "Inconsistent ikpt", iqbz, bz2ibz_listkk(iqbz,1), bz2ibz_symkpt_new(1,iqbz)
   end if
   if (bz2ibz_listkk(iqbz, 2) /= bz2ibz_symkpt_new(2, iqbz)) then
     ierr = ierr + 1
     if (my_rank == master) write(std_out,*) "Inconsistent isym", iqbz, bz2ibz_listkk(iqbz,2), bz2ibz_symkpt_new(2,iqbz)
   end if
   if (bz2ibz_listkk(iqbz, 6) /= bz2ibz_symkpt_new(3, iqbz)) then
     ierr = ierr + 1
     if (my_rank == master ) write(std_out,*) "Inconsistent itim", iqbz, bz2ibz_listkk(iqbz,6), bz2ibz_symkpt_new(3,iqbz)
   end if
   if (.not.all(bz2ibz_listkk(iqbz, 3:5) == bz2ibz_symkpt_new(4:, iqbz))) then
     ierr = ierr + 1
     if (my_rank == master) then
       write(std_out,*) "Inconsistent shift:", iqbz
       write(std_out,*) bz2ibz_listkk(iqbz, 3:5)
       write(std_out,*) bz2ibz_symkpt_new(4:, iqbz)
     end if
   end if
 end do

 ABI_SFREE(bz2ibz_symkpt)
 ABI_SFREE(bz2ibz_symkpt_new)
 ABI_SFREE(bz2ibz_listkk)
 ABI_SFREE(wtq_fullbz)
 ABI_SFREE(wtq_folded)
 ABI_SFREE(ibz2bz)
 ABI_SFREE(ibz2bz_new)
 ABI_SFREE(qbz)
 ABI_SFREE(qibz)
 ABI_SFREE(bz2ibz)
 ABI_SFREE(wtq_ibz)

 call cryst%free()
 call krank%free()

end subroutine kptrank_unittests
!!***

end module m_unittests
