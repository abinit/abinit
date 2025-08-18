!!****m*ABINIT/m_anaddb_driver
!! NAME
!!  m_anaddb_driver
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2025 ABINIT group (GA)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_anaddb_driver

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_nctk
 use netcdf
 use m_ifc
 use m_ddb_hdr
 use m_phonons
 use m_supercell
 use m_raman

 use m_fstrings,       only : strcat
 use m_crystal,        only : crystal_t
 use m_anaddb_dataset, only : anaddb_dataset_type
 use m_ddb,            only : ddb_type, asrq0_t
 use m_dynmat,         only : gtdyn9, dfpt_phfrq, dfpt_prtph
 use m_ddb_interpolate, only : ddb_interpolate
 use m_harmonic_thermo, only : harmonic_thermo
 use m_elphon,         only : elphon
 use m_thmeig,         only : thmeig
 use m_relaxpol,       only : relaxpol
 use m_ddb_diel,       only : ddb_diel
 use m_ddb_elast,      only : ddb_elast
 use m_ddb_piezo,      only : ddb_piezo
 use m_ddb_internalstr, only : ddb_internalstr
 use m_gruneisen,      only : gruns_anaddb
 use m_ddb_flexo,      only : ddb_flexo
 use m_lwf,            only : run_lattice_wannier

 implicit none

 private

 public:: anaddb_driver_type

!----------------------------------------------------------------------

!!****t* m_anaddb_driver/anaddb_driver_type
!! NAME
!! anaddb_driver_type
!!
!! FUNCTION
!! The anaddb_driver_type structured datatype
!! Main subroutines run by anaddb.
!!
!! SOURCE

 type anaddb_driver_type

   integer:: natom
   integer:: msize
   integer:: mpert
   integer:: usepaw  ! GA: TODO Remove usepaw

   logical:: do_ifc=.false.
   logical:: do_electric_tensors=.false.
   logical:: do_dielectric_q0=.false.
   logical:: do_dielectric_nonana=.false.
   logical:: do_phonon_dos=.false.
   logical:: do_phonon_bs=.false.

   real(dp):: epsinf(3, 3)
   real(dp):: dchide(3,3,3)
   real(dp):: dielt_rlx(3, 3)
   real(dp):: elast(6, 6)

   real(dp), allocatable:: d2cart(:,:)
   ! d2cart(2,msize)

   real(dp), allocatable:: displ(:)
   ! displ(2*3*natom*3*natom)

   real(dp), allocatable:: phfrq(:)
   ! phfrq(3*natom)

   real(dp), allocatable:: instrain(:,:)
   ! instrain(3*natom,6)

   real(dp), allocatable:: dchidt(:,:,:,:)
   ! dchidt(natom,3,3,3)

   real(dp), allocatable:: fact_oscstr(:,:,:)
   ! fact_oscstr(2,3,3*natom)

   real(dp), allocatable:: zeff(:,:,:)
   ! zeff(3,3,natom)

   real(dp), allocatable:: qdrp_cart(:,:,:,:)
   ! qdrp_cart(3,3,3,natom)

 contains

   procedure :: init => anaddb_driver_init
   ! Initialize object

   procedure :: free => anaddb_driver_free
   ! Free memory

   procedure :: open_write_nc => anaddb_driver_open_write_nc
   ! Open netcdf output and write some info

   procedure :: electric_tensors => anaddb_driver_electric_tensors
   ! Compute dielectric tensor, born effective charges, and quadrupoles

   procedure :: structural_response => anaddb_driver_structural_response
   ! Structural response at fixed polarization

   procedure :: susceptibilities => anaddb_driver_susceptibilities
   ! Compute non-linear optical susceptibilities and first-order derivatives

   procedure :: interatomic_force_constants => anaddb_driver_interatomic_force_constants
   ! Compute the interatomic force constants from a ddb

   procedure :: phdos => anaddb_driver_phdos
   ! Compute phonon density of states

   procedure :: thermal_supercell => anaddb_driver_thermal_supercell
   ! Thermal supercell calculation

   procedure :: harmonic_thermo => anaddb_driver_harmonic_thermo
   ! Phonon density of states and thermodynamical properties

   procedure :: dielectric_q0 => anaddb_driver_dielectric_q0
   ! Dielectric tensor and related properties

   procedure :: nonlinear_response => anaddb_driver_nonlinear_response
   ! Non-linear response: electrooptic and Raman

   procedure :: dielectric_nonana => anaddb_driver_dielectric_nonana
   ! Non-analyticity in the dielectric matrix and raman susceptibility

   procedure :: internal_strain => anaddb_driver_internal_strain
   ! Internal strain

   procedure :: elastic_tensor => anaddb_driver_elastic_tensor
   ! Elastic tensor

   procedure :: piezoelectric_tensor => anaddb_driver_piezoelectric_tensor
   ! Piezoelectric tensor

   procedure :: flexoelectric_tensor => anaddb_driver_flexoelectric_tensor
   ! Flexoelectric tensor

   procedure :: lattice_wannier => anaddb_driver_lattice_wannier
   ! Construct the Lattice Wannier functions

 end type anaddb_driver_type

contains

!!****f* m_anaddb_driver/anaddb_driver_init
!! NAME
!! anaddb_driver_init
!!
!! FUNCTION
!! Initialize object
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_init(driver, dtset)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver
 type(anaddb_dataset_type), intent(in):: dtset

! ************************************************************************

 ! Set control flags
 if (dtset%ifcflag == 1) then
   driver%do_ifc = .true. 
 end if

 if (dtset%ifcflag /= 0 .or. dtset%dieflag /= 0 &
&    .or. dtset%nph2l /= 0 .or. dtset%nlflag /= 0 &
&    .or. dtset%piezoflag /= 0 .or. dtset%flexoflag /= 0 &
&    .or. dtset%polflag /= 0) then
   driver%do_electric_tensors = .true.
 end if

 if ((dtset%dieflag /= 0 .and. dtset%dieflag /= 2) &
&    .or. dtset%nph2l /= 0 .or. dtset%nlflag == 1 &
&    .or. dtset%piezoflag /= 0) then
   driver%do_dielectric_q0 = .true.
 end if

 if (dtset%nph2l /= 0) then
   driver%do_dielectric_nonana = .true. 
 end if

 if (dtset%ifcflag == 1 .and. any(dtset%prtdos==[1, 2])) then
   driver%do_phonon_dos = .true.
 end if

 if (dtset%nph1l /= 0 .or. dtset%nqpath /= 0) then
   driver%do_phonon_bs = .true. 
 end if

 if (dtset%gruns_nddbs /= 0) then
   driver%do_ifc = .false. 
   driver%do_electric_tensors = .false.
   driver%do_dielectric_q0 = .false.
   driver%do_dielectric_nonana = .false. 
   driver%do_phonon_bs = .false. 
   driver%do_phonon_dos = .false.
 end if

 ! Copy dimensions
 driver%natom = dtset%natom
 driver%msize = dtset%msize
 driver%mpert = dtset%mpert
 driver%usepaw = dtset%usepaw

 ! Allocate memory
 ABI_MALLOC(driver%d2cart, (2, driver%msize))
 ABI_MALLOC(driver%displ, (2*3*driver%natom*3*driver%natom))
 ABI_MALLOC(driver%phfrq, (3*driver%natom))
 ABI_MALLOC(driver%instrain, (3*driver%natom, 6))

 ! Electric tensors
 ABI_MALLOC(driver%zeff, (3, 3, driver%natom))
 ABI_MALLOC(driver%qdrp_cart, (3, 3, 3, driver%natom))

 ! oscillator strength and Lyddane-Sachs-Teller relation
 ABI_MALLOC(driver%fact_oscstr, (2, 3, 3*driver%natom))

 ! Susceptibilities
 if (dtset%nlflag > 0) then
   ABI_MALLOC(driver%dchidt, (driver%natom, 3, 3, 3))
 end if

end subroutine anaddb_driver_init
!!***


!!****f* m_anaddb_driver/anaddb_driver_free
!! NAME
!! anaddb_driver_free
!!
!! FUNCTION
!! Deallocate memory
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_free(driver)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver

! ************************************************************************

 ABI_SFREE(driver%zeff)
 ABI_SFREE(driver%qdrp_cart)
 ABI_SFREE(driver%d2cart)
 ABI_SFREE(driver%displ)
 ABI_SFREE(driver%phfrq)
 ABI_SFREE(driver%dchidt)
 ABI_SFREE(driver%instrain)
 ABI_SFREE(driver%fact_oscstr)

end subroutine anaddb_driver_free
!!***

!!****f* m_anaddb_driver/anaddb_driver_open_write_nc
!! NAME
!! anaddb_driver_open_write_nc
!!
!! FUNCTION
!! Open anaddb netcdf output and define dimensions.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_open_write_nc(driver, ana_ncid, dtset, crystal, comm)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(crystal_t),intent(in):: crystal
 integer,intent(in):: comm
 integer, intent(out):: ana_ncid

!Local variables -------------------------------
 integer, parameter:: master = 0
 integer:: natom,lenstr
 integer:: ncerr
 integer:: my_rank

! ************************************************************************

 my_rank = xmpi_comm_rank(comm)

 natom = driver%natom
 lenstr = dtset%lenstr

 ! Open the netcdf file that will contain the anaddb results
 ana_ncid = nctk_noid
 if (my_rank == master) then
   NCF_CHECK_MSG(nctk_open_create(ana_ncid, trim(dtset%prefix_outdata)//"_anaddb.nc", xmpi_comm_self), "Creating anaddb.nc")
   ncerr = nctk_def_dims(ana_ncid, [ &
       nctkdim_t('number_of_atoms', natom), &
       nctkdim_t('natom3', 3*natom), &
       nctkdim_t('number_of_phonon_modes', 3*natom), &
       nctkdim_t('anaddb_input_len', lenstr) &
   ], defmode=.True.)
   NCF_CHECK(ncerr)
   ncerr = nctk_def_arrays(ana_ncid, [ &
     nctkarr_t("anaddb_input_string", "char", "anaddb_input_len") &
   ])
   NCF_CHECK(ncerr)
   !NCF_CHECK(nctk_defnwrite_ivars(ana_ncid, ["anaddb_version"], [1]))

   ncerr = nctk_def_iscalars(ana_ncid, [character(len = nctk_slen) :: &
       "asr", "chneut", "dipdip", "symdynmat", "dipquad", "quadquad"])
   NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ana_ncid))
   ncerr = nf90_put_var(ana_ncid, nctk_idname(ana_ncid, "anaddb_input_string"), dtset%input_string(:lenstr))
   NCF_CHECK(ncerr)
   NCF_CHECK(crystal%ncwrite(ana_ncid))

   ncerr = nctk_write_iscalars(ana_ncid, [character(len = nctk_slen) :: &
     "asr", "chneut", "dipdip", "symdynmat", "dipquad", "quadquad"], &
     [dtset%asr, dtset%chneut, dtset%dipdip, &
      dtset%symdynmat, dtset%dipquad, dtset%quadquad])
   NCF_CHECK(ncerr)
 end if

end subroutine anaddb_driver_open_write_nc
!!***

!!****f* m_anaddb_driver/anaddb_driver_electric_tensors
!! NAME
!! anaddb_driver_electric_tensors
!!
!! FUNCTION
!! get Dielectric tensor, born effective charges, and quadrupole tensor,
!! and write them to netcdf output.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_electric_tensors(driver, dtset, crystal, ddb, ddb_lw, ddb_hdr, ana_ncid, comm)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver
 type(anaddb_dataset_type), intent(inout):: dtset
 type(crystal_t), intent(in):: crystal
 type(ddb_type), intent(inout):: ddb
 type(ddb_type), intent(inout):: ddb_lw
 type(ddb_hdr_type), intent(in):: ddb_hdr
 integer, intent(in):: ana_ncid
 integer, intent(in):: comm

!Local variables -------------------------------
 integer, parameter:: master = 0
 integer:: my_rank
 integer:: ii
 integer:: iblok, iblok_quadrupoles, iblok_epsinf
 integer:: ncerr
 integer:: lwsym
 character(len = 500):: msg
 integer:: units(2)

! ************************************************************************

 units = [std_out, ab_out]

 ! Get Quadrupole tensor
 iblok_quadrupoles = 0
 driver%qdrp_cart = zero
 if (ddb_hdr%has_d3E_lw) then
   write(msg, '(2a, (80a), 2a)') ch10, ('=',ii = 1, 80)
   call wrtout(units, msg)
   lwsym = 1
   iblok_quadrupoles = ddb_lw%get_quadrupoles(ddb_hdr%ddb_version, lwsym, BLKTYP_d3E_lw, driver%qdrp_cart)
 end if

 ! The default value is 1. Here we set the flags to zero if Q*is not available.
 ! GA: FIXME This change of variable value happens after outputting the input.
 !           It should happen before, but I will need to update tests when I change this.
 if (iblok_quadrupoles == 0) then
   dtset%dipquad = 0
   dtset%quadquad = 0
 end if

 ! Get the electronic dielectric tensor (epsinf) and Born effective charges (zeff)
 ! (initialized to one_3D and zero if the derivatives are not available in the DDB file)
 iblok = ddb%get_dielt_zeff(crystal, dtset%rfmeth, dtset%chneut, dtset%selectz, driver%epsinf, driver%zeff)

 ! Try to get epsinf, in case just the DDE are present
 if (iblok == 0) then
   iblok_epsinf = ddb%get_dielt(dtset%rfmeth, driver%epsinf)
 end if

 !if (iblok_epsinf == 0) then
 !GA: Not the greatest way of checking
 if (driver%epsinf(1, 1)==one .and. driver%epsinf(2, 2)==one .and. driver%epsinf(3, 3)==one) then
   write(msg, '(5a)') ch10, &
     ' The DDB file does not contain the derivatives w.r.t. electric field perturbation. ',ch10, &
     ' The program will continue by setting the electronic dielectric tensor to 1. ',ch10
  ! call wrtout([ab_out], msg)
   ABI_WARNING(msg)
 end if

!**********************************************************************
! Write Dielectric tensor, born effective charges, and quadrupoles to netcdf output.

 my_rank = xmpi_comm_rank(comm)
 if (my_rank == master) then
   ncerr = nctk_def_arrays(ana_ncid, [&
   nctkarr_t('emacro_cart', "dp", 'number_of_cartesian_directions, number_of_cartesian_directions'), &
   nctkarr_t('quadrupoles_cart', "dp", 'three, three, three, number_of_atoms'), &
   nctkarr_t('becs_cart', "dp", "number_of_cartesian_directions, number_of_cartesian_directions, number_of_atoms")],&
   defmode=.True.)
   NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ana_ncid))
   NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, 'emacro_cart'), driver%epsinf))
   NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, 'quadrupoles_cart'), driver%qdrp_cart))
   NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, 'becs_cart'), driver%zeff))
 end if

end subroutine anaddb_driver_electric_tensors
!!***

!!****f* m_anaddb_driver/anaddb_driver_structural_response
!! NAME
!! anaddb_driver_structural_response
!!
!! FUNCTION
!! Compute structural response at fixed polarization
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_structural_response(driver, dtset, crystal, ddb)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(crystal_t), intent(in):: crystal
 type(ddb_type), intent(inout):: ddb

!Local variables -------------------------------
 integer:: iblok
 integer:: msize
 real(dp):: etotal
 character(len = 500):: msg
 integer:: rfelfd(4), rfphon(4), rfstrs(4)
 real(dp):: red_ptot(3)
 real(dp):: pel(3)
 real(dp):: strten(6)
 real(dp) :: targetpol(3)
 real(dp):: qphnrm(3), qphon(3, 3)
 integer, allocatable:: d2flg(:)
 real(dp), allocatable:: gred(:,:)

! ************************************************************************

 msize = dtset%msize
 ABI_MALLOC(d2flg, (msize))

 ! Look for the Gamma Block in the DDB
 qphon(:,1)=zero
 qphnrm(1)=zero
 rfphon(1:2)=1
 rfelfd(1:2)=2
 rfstrs(1:2)=0

 !write(std_out,*)"ddb%mpert",ddb%mpert
 call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, dtset%rfmeth)
 !iblok = ddb%get_dielt_zeff(crystal, dtset%rfmeth, dtset%chneut, dtset%selectz, driver%epsinf, driver%zeff)

 if(iblok /= 0)then
   ! Save the second-order derivatives
   driver%d2cart(1:2, 1:msize) = ddb%val(1:2, 1:msize, iblok)
   d2flg(1:msize) = ddb%flg(1:msize, iblok)

 else
   ! the gamma blok has not been found
   if (dtset%relaxat == 0 .and. dtset%relaxstr == 0) then
     ! The gamma blok is not needed
     driver%d2cart(1:2, 1:msize)=zero
     d2flg(1:msize)=1
   else
     ! There is a problem !
     write(msg, '(7a)' )&
      'The dynamical matrix at Gamma is needed, in order to perform ',ch10, &
      "relaxation at constant polarisation (Na Sai's method)",ch10, &
      'However, this was not found in the DDB.',ch10, &
      'Action: complete your DDB with the dynamical matrix at Gamma.'
     ABI_ERROR(msg)
   end if
 end if  ! iblok not found

 ! Extract the block with the total energy
 if (ddb%get_etotal(etotal) == 0) then
   ABI_ERROR("DDB file does not contain GS etotal")
 end if

 ! Extract the polarizability
 iblok = ddb%get_pel(pel, dtset%relaxat, dtset%relaxstr)

 ! Extract the forces
 iblok = ddb%get_gred(gred, dtset%relaxat, dtset%relaxstr)

 ! Extract the stress tensor
 iblok = ddb%get_strten(strten, dtset%relaxat, dtset%relaxstr)

 ! when called from here red_ptot is not set  ! So set it to zero
 red_ptot(:) = zero

 targetpol(:) = dtset%targetpol

 !GA: FIXME: Remove usepaw
 call relaxpol(crystal, d2flg, driver%d2cart, etotal, gred, dtset%iatfix, &
&   ab_out, dtset%istrfix, dtset%mpert, dtset%msize, dtset%natfix, crystal%natom, &
&   dtset%nstrfix, pel, red_ptot, dtset%relaxat, dtset%relaxstr, &
&   strten, targetpol, dtset%usepaw)

 ABI_SFREE(gred)
 ABI_FREE(d2flg)

end subroutine anaddb_driver_structural_response
!!***

!!****f* m_anaddb_driver/anaddb_driver_susceptibilities
!! NAME
!! anaddb_driver_susceptibilities
!!
!! FUNCTION
!! Compute non-linear optical susceptibilities,
!! and if dtset%nlflag < 3, compute first-order change
!! in the linear dielectric susceptibility induced by an atomic displacement.
!! Then write susceptibilites to netcdf output.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_susceptibilities(driver, dtset, ddb, ana_ncid, comm)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(ddb_type), intent(in):: ddb
 integer, intent(in):: ana_ncid
 integer, intent(in):: comm

!Local variables -------------------------------
 integer, parameter:: master = 0
 integer:: my_rank
 integer:: ncerr

! ************************************************************************

 if (ddb%get_dchidet(dtset%ramansr, dtset%nlflag, driver%dchide, driver%dchidt) == 0) then
   ABI_ERROR("Cannot find block corresponding to non-linear optical susceptibilities in DDB file")
 end if

 ! Save to the netcdf
 my_rank = xmpi_comm_rank(comm)
 if (my_rank == master) then
   ncerr = nctk_def_arrays(ana_ncid, [nctkarr_t("dchide", "dp", "three, three, three")], defmode=.True.)
   NCF_CHECK(ncerr)
   NCF_CHECK(nctk_set_datamode(ana_ncid))
   NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, "dchide"), driver%dchide))

   ! dchidt only present if nlflag == 1 or 2
   if (dtset%nlflag < 3) then
     ncerr = nctk_def_arrays(ana_ncid, [nctkarr_t("dchidt", "dp", &
     "number_of_atoms, three, three, three")], defmode=.True.)
     NCF_CHECK(ncerr)
     NCF_CHECK(nctk_set_datamode(ana_ncid))
     NCF_CHECK(nf90_put_var(ana_ncid, nctk_idname(ana_ncid, "dchidt"), driver%dchidt))
   end if
 end if

end subroutine anaddb_driver_susceptibilities
!!***

!!****f* m_anaddb_driver/anaddb_driver_interatomic_force_constants
!! NAME
!! anaddb_driver_interatomic_force_constants
!!
!! FUNCTION
!! Interatomic forces calculation
!! Compute the interatomic force constants from a ddb.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_interatomic_force_constants(driver, ifc, dtset, crystal, ddb, ana_ncid, comm)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver
 type(ifc_type), intent(out):: ifc
 type(anaddb_dataset_type), intent(in):: dtset
 type(crystal_t), intent(in):: crystal
 type(ddb_type), intent(in):: ddb
 integer, intent(in):: ana_ncid
 integer, intent(in):: comm

!Local variables -------------------------------
 integer, parameter:: master = 0
 integer:: my_rank
 integer:: ii
 type(ifc_type):: Ifc_coarse
 integer:: ngqpt_coarse(3)

! ************************************************************************

 if (any(dtset%qrefine(:) > 1)) then
   ! Gaal-Nagy's algorithm in PRB 73 014117 [[cite:GaalNagy2006]]
   ! Build the IFCs using the coarse q-mesh.
   do ii = 1, 3
     ngqpt_coarse(ii) = dtset%ngqpt(ii) / dtset%qrefine(ii)
   end do
   call Ifc_coarse%init(crystal, ddb, &
     dtset%brav, dtset%asr, dtset%symdynmat, dtset%dipdip, dtset%rfmeth, ngqpt_coarse, dtset%nqshft, dtset%q1shft, &
     driver%epsinf, driver%zeff, driver%qdrp_cart, &
     dtset%nsphere, dtset%rifcsph, dtset%prtsrlr, dtset%enunit, comm, &
     dipquad=dtset%dipquad, quadquad=dtset%quadquad)

   ! Now use the coarse q-mesh to fill the entries in dynmat(q)
   ! on the dense q-mesh that cannot be obtained from the DDB file.
   call ifc%init(crystal, ddb, &
    dtset%brav, dtset%asr, dtset%symdynmat, dtset%dipdip, dtset%rfmeth, &
    dtset%ngqpt(1:3), dtset%nqshft, dtset%q1shft, driver%epsinf, driver%zeff, driver%qdrp_cart, &
    dtset%nsphere, dtset%rifcsph, dtset%prtsrlr, dtset%enunit, comm, &
    Ifc_coarse=Ifc_coarse, dipquad=dtset%dipquad, quadquad=dtset%quadquad)
   call Ifc_coarse%free()

 else
   call ifc%init(crystal, ddb, &
     dtset%brav, dtset%asr, dtset%symdynmat, dtset%dipdip, dtset%rfmeth, &
     dtset%ngqpt(1:3), dtset%nqshft, dtset%q1shft, driver%epsinf, driver%zeff, driver%qdrp_cart, &
     dtset%nsphere, dtset%rifcsph, dtset%prtsrlr, dtset%enunit, comm, &
     dipquad=dtset%dipquad, quadquad=dtset%quadquad)
 end if

 call ifc%print(unit=std_out)

 ! Compute speed of sound.
 if (dtset%vs_qrad_tolkms(1) > zero) then
   call ifc%speedofsound(crystal, dtset%vs_qrad_tolkms, ana_ncid, comm)
 end if

 ! Print analysis of the real-space interatomic force constants
 ! TODO: ifc_out should not have side effects
 my_rank = xmpi_comm_rank(comm)
 if (my_rank == master .and. dtset%ifcout /= 0) then
   call ifc%write(dtset%ifcana, dtset%atifcflg, dtset%ifcout, dtset%prt_ifc, ana_ncid, dtset%prefix_outdata)
 end if

end subroutine anaddb_driver_interatomic_force_constants
!!***

!!****f* m_anaddb_driver/anaddb_driver_phdos
!! NAME
!! anaddb_driver_phdos
!!
!! FUNCTION
!! Compute phonon density of states.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_phdos(driver, dtset, crystal, ifc, comm)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(crystal_t), intent(in):: crystal
 type(ifc_type), intent(in):: ifc
 integer, intent(in):: comm

!Local variables -------------------------------
 integer, parameter:: master = 0
 integer:: my_rank
 integer:: ii
 integer:: phdos_ncid, ncerr
 character(len = fnlen):: phibz_prefix
 character(len = 500):: msg
 type(phdos_t):: Phdos

 integer:: units(2)
 integer:: count_wminmax(2)
 real(dp):: wminmax(2)

! ************************************************************************

 my_rank = xmpi_comm_rank(comm)
 units = [std_out, ab_out]

 write(msg, '(a, (80a), 4a)')ch10, ('=',ii = 1, 80), ch10, ch10, ' Calculation of phonon density of states ',ch10
 call wrtout(units, msg)

 ! Only 1 shift in q-mesh
 wminmax = zero
 phibz_prefix = trim(" ")
 !phibz_prefix = "freq_displ" ! Uncomment this line to activate output of PHIBZ
 !                              ^^^^  GA: What the hell? FIXME ^^^^^
 do
   call Phdos%init(crystal, Ifc, dtset%prtdos, dtset%dosdeltae, dtset%dossmear, dtset%ng2qpt, 1, dtset%q2shft, &
                   phibz_prefix, wminmax, count_wminmax, comm, dos_maxmode=dtset%dos_maxmode)
   if (all(count_wminmax == 0)) exit
   wminmax(1) = wminmax(1) - abs(wminmax(1)) * 0.05; wminmax(2) = wminmax(2) + abs(wminmax(2)) * 0.05
   call phdos%free()
   write(msg, "(a, 2f8.5)")"Initial frequency mesh not large enough. Recomputing PHDOS with wmin, wmax: ",wminmax
   call wrtout(std_out, msg)
 end do

 if (my_rank == master) then
   call phdos%print_msqd(dtset%prefix_outdata, dtset%ntemper, dtset%tempermin, dtset%temperinc)
   call phdos%print(strcat(dtset%prefix_outdata, "_PHDOS"))
   call phdos%print_debye(crystal%ucvol)
   call phdos%print_thermo(strcat(dtset%prefix_outdata, "_THERMO"), dtset%ntemper, dtset%tempermin, dtset%temperinc)

   ncerr = nctk_open_create(phdos_ncid, strcat(dtset%prefix_outdata, "_PHDOS.nc"), xmpi_comm_self)
   NCF_CHECK_MSG(ncerr, "Creating PHDOS.nc file")
   NCF_CHECK(crystal%ncwrite(phdos_ncid))
   call phdos%ncwrite(phdos_ncid)
   NCF_CHECK(nf90_close(phdos_ncid))
 end if

 call phdos%free()

end subroutine anaddb_driver_phdos
!!***

!!****f* m_anaddb_driver/anaddb_driver_thermal_supercell
!! NAME
!! anaddb_driver_thermal_supercell
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_thermal_supercell(driver, dtset, crystal, ifc)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(in):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(crystal_t), intent(in):: crystal
 type(ifc_type), intent(in):: ifc

!Local variables -------------------------------
 type(supercell_type), allocatable:: thm_scells(:)

! ************************************************************************

 ABI_MALLOC(thm_scells, (dtset%ntemper))
 call zacharias_supercell_make(crystal, ifc, dtset%ntemper, dtset%thermal_supercell, dtset%tempermin, dtset%temperinc, thm_scells)
 call zacharias_supercell_print(dtset%prefix_outdata, dtset%ntemper, dtset%tempermin, dtset%temperinc, thm_scells)
 call thermal_supercell_free(dtset%ntemper, thm_scells)
 ABI_FREE(thm_scells)

end subroutine anaddb_driver_thermal_supercell
!!***

!!****f* m_anaddb_driver/anaddb_driver_harmonic_thermo
!! NAME
!! anaddb_driver_harmonic_thermo
!!
!! FUNCTION
!! Phonon density of states and thermodynamical properties calculation
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_harmonic_thermo(driver, dtset, crystal, ifc, ddb, comm)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(in):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(crystal_t), intent(in):: crystal
 type(ifc_type), intent(in):: ifc
 type(ddb_type), intent(in):: ddb
 integer, intent(in):: comm

!Local variables -------------------------------
 integer:: ii
 character(len = 500):: msg
 integer:: units(2)

! ************************************************************************

 units = [std_out, ab_out]

 write(msg, '(a, (80a), a, a, a, a, a, a, a, a)' ) ch10, ('=',ii = 1, 80), ch10, ch10, &
  ' Calculation of phonon density of states, ',ch10, &
  '    thermodynamical properties, ',ch10, &
  '    and Debye-Waller factors.',ch10
 call wrtout(units, msg)

 if (dtset%thmflag == 1) then
   call harmonic_thermo(Ifc, crystal, ddb%amu, dtset, ab_out, dtset%prefix_outdata, comm)

 else if (dtset%thmflag == 2) then
   write(msg, '(a, (80a), a, a, a, a)' ) ch10, ('=',ii = 1, 80), ch10, ch10, ' Entering thm9 routine with thmflag = 2 ',ch10
   call wrtout(units, msg)
   call harmonic_thermo(Ifc, crystal, ddb%amu, dtset, ab_out, dtset%prefix_outdata, comm, thmflag=dtset%thmflag)
 end if

end subroutine anaddb_driver_harmonic_thermo
!!***

!!****f* m_anaddb_driver/anaddb_driver_dielectric_q0
!! NAME
!! anaddb_driver_dielectric_q0
!!
!! FUNCTION
!! Compute dielectric tensor at Gamma and related properties:
!! mode effective charges, oscillator strength.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_dielectric_q0(driver, dtset, crystal, ifc, ddb, asrq0, ana_ncid, comm)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(crystal_t), intent(in):: crystal
 type(ifc_type), intent(in):: ifc
 type(ddb_type), intent(in):: ddb
 type(asrq0_t), intent(inout):: asrq0
 integer, intent(in):: ana_ncid
 integer, intent(in):: comm

!Local variables -------------------------------
 integer:: ii, iblok
 integer:: rfelfd(4), rfphon(4), rfstrs(4)
 integer:: units(2)
 character(len = 500):: msg
 real(dp):: qphnrm(3), qphon(3, 3)
 real(dp), allocatable:: eigval(:,:)
 real(dp), allocatable:: eigvec(:,:,:,:,:)
 real(dp), allocatable:: lst(:)

! ************************************************************************

 units = [std_out, ab_out]

 ABI_MALLOC(eigval, (3, driver%natom))
 ABI_MALLOC(eigvec, (2, 3, driver%natom, 3, driver%natom))
 ABI_MALLOC(lst, (dtset%nph2l+1))

 lst = zero

 !***************************************************************
 ! Generates the dynamical matrix at Gamma
 qphon(:,1)=zero; qphnrm(1)=zero
 ! Generation of the dynamical matrix in cartesian coordinates
 if (dtset%ifcflag == 1) then
   ! Get d2cart using the interatomic forces and the
   ! long-range coulomb interaction through Ewald summation
   call gtdyn9(ddb%acell, Ifc%atmfrc, driver%epsinf, Ifc%dipdip, &
     Ifc%dyewq0, driver%d2cart, crystal%gmet, ddb%gprim, dtset%mpert, crystal%natom, &
     Ifc%nrpt, qphnrm(1), qphon, crystal%rmet, ddb%rprim, Ifc%rpt, &
     Ifc%trans, crystal%ucvol, Ifc%wghatm, crystal%xred, driver%zeff, driver%qdrp_cart, &
     Ifc%ewald_option, xmpi_comm_self, &
     dipquad=Ifc%dipquad, quadquad=Ifc%quadquad)

 else if (dtset%ifcflag == 0) then
   ! Look for the information in the DDB
   rfphon(1:2)=1; rfelfd(1:2)=2; rfstrs(1:2)=0
   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, dtset%rfmeth)
   if (iblok == 0) then
     driver%d2cart(:,1:dtset%msize)=zero
     ! GA: I notice this situation happen in test tutorespfn[telast_3]
     !     and I dont understand why the block is not found.
   else
     ! Copy the dynamical matrix in d2cart
     driver%d2cart(:,1:dtset%msize)=ddb%val(:,:,iblok)
     ! Eventually impose the acoustic sum rule
     call asrq0%apply(crystal%natom, dtset%mpert, dtset%msize, crystal%xcart, driver%d2cart)
   end if

 end if  ! end of the generation of the dynamical matrix at gamma.
 !***************************************************************

 ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
 call dfpt_phfrq(ddb%amu, driver%displ, driver%d2cart, eigval, eigvec, crystal%indsym, &
   dtset%mpert, crystal%nsym, crystal%natom, crystal%nsym, crystal%ntypat, driver%phfrq, qphnrm(1), qphon, &
   crystal%rprimd, dtset%symdynmat, crystal%symrel, crystal%symafm, crystal%typat, crystal%ucvol)

 ! calculation of the oscillator strengths, mode effective charge and
 ! dielectric tensor, frequency dependent dielectric tensor (dieflag)
 ! and mode by mode decomposition of epsilon if dieflag == 3
 if (dtset%dieflag /= 0) then
   if (driver%epsinf(1, 1)==one .and. driver%epsinf(2, 2)==one .and. driver%epsinf(3, 3)==one) then
     write(msg, '(7a)') ch10, &
      ' The DDB file does not contain the derivatives w.r.t. electric field perturbation. ',ch10, &
      ' This is mandatory to calculate the dielectric constant, ',ch10, &
      ' Please check your DDB file or use dieflag = 0.',ch10
     ABI_ERROR(msg)
   end if

   write(msg, '(a, (80a), a)' ) ch10, ('=',ii = 1, 80), ch10
   call wrtout(units, msg)

   ! Print the electronic contribution to the dielectric tensor
   ! It can be extracted directly from the DDB if perturbation with E-field is present
   call ddb_diel(crystal, ddb%amu, dtset, driver%dielt_rlx, driver%displ, driver%d2cart, driver%epsinf, driver%fact_oscstr, &
     ab_out, lst, dtset%mpert, crystal%natom, 0, driver%phfrq, comm, ana_ncid)
 end if

 ABI_SFREE(eigval)
 ABI_SFREE(eigvec)
 ABI_SFREE(lst)

end subroutine anaddb_driver_dielectric_q0
!!***

!!****f* m_anaddb_driver/anaddb_driver_nonlinear_response
!! NAME
!! anaddb_driver_nonlinear_response
!!
!! FUNCTION
!! Non-linear response: electrooptic and Raman (q = Gamma, TO modes only)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_nonlinear_response(driver, dtset, crystal, ana_ncid, comm)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(crystal_t), intent(in):: crystal
 integer, intent(in):: ana_ncid
 integer, intent(in):: comm

!Local variables -------------------------------
 integer, parameter:: master = 0
 integer:: my_rank
 real(dp):: qphnrm(3), qphon(3, 3)
 real(dp), allocatable:: rsus(:,:,:)

! ************************************************************************

 my_rank = xmpi_comm_rank(comm)

 ABI_MALLOC(rsus, (3*driver%natom, 3, 3))

 ! Raman susceptibilities for the 1st list (only TO  modes at q = Gamma)
 qphon(:,1)=zero
 qphnrm(1)=zero
 call ramansus(driver%d2cart, driver%dchide, driver%dchidt, driver%displ, dtset%mpert, crystal%natom, driver%phfrq, qphon, qphnrm(1), rsus, crystal%ucvol)

 if (my_rank == master) then
   call defwrite_raman_terms(ana_ncid, crystal%natom, rsus, driver%phfrq)
 end if

 ! EO coef:
 call electrooptic(driver%dchide, dtset%dieflag, driver%epsinf, driver%fact_oscstr, crystal%natom, driver%phfrq, dtset%prtmbm, rsus, crystal%ucvol)

 ABI_SFREE(rsus)

end subroutine anaddb_driver_nonlinear_response
!!***

!!****f* m_anaddb_driver/anaddb_driver_dielectric_nonana
!! NAME
!! anaddb_driver_dielectric_nonana
!!
!! FUNCTION
!! Compute non-analyticity in the dielectric matrix and raman susceptibility
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_dielectric_nonana(driver, dtset, crystal, ddb, ana_ncid, comm)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(crystal_t), intent(in):: crystal
 type(ddb_type), intent(in):: ddb
 integer, intent(in):: ana_ncid
 integer, intent(in):: comm

!Local variables -------------------------------
 integer, parameter:: master = 0
 integer:: ii, iphl2
 integer:: natom, nph2l
 integer:: my_rank
 integer:: units(2)
 character(len = 500):: msg
 real(dp):: qphnrm(3), qphon(3, 3)
 real(dp), allocatable:: eigval(:,:)
 real(dp), allocatable:: eigvec(:,:,:,:,:)
 real(dp), allocatable:: rsus(:,:,:)
 real(dp), allocatable:: lst(:)

! ************************************************************************

 my_rank = xmpi_comm_rank(comm)
 units = [std_out, ab_out]

 natom = dtset%natom
 nph2l = dtset%nph2l
 ABI_MALLOC(eigval, (3, natom))
 ABI_MALLOC(eigvec, (2, 3, natom, 3, natom))
 ABI_MALLOC(rsus, (3*natom, 3, 3))
 ABI_MALLOC(lst, (nph2l+1))
 lst = zero

 write(msg, '(a, (80a), a, a, a, a)' ) ch10, ('=',ii = 1, 80), ch10, ch10, ' Treat the second list of vectors ',ch10
 call wrtout(units, msg)

 if (my_rank == master) then
   iphl2 = 0
   call defwrite_nonana_terms(ana_ncid, iphl2, nph2l, dtset%qph2l, dtset%natom, driver%phfrq, driver%displ, "define")
   if (dtset%nlflag == 1) then
     call defwrite_nonana_raman_terms(ana_ncid, iphl2, nph2l, dtset%natom, rsus, "define")
   end if
 end if

 !Get the log of product of the square of the phonon frequencies without non-analyticities (q = 0)
 !For the Lyddane-Sachs-Teller relation, it is stored in lst(nph2+1)
 if (dtset%dieflag /= 2 .and. dtset%dieflag /= 0) then
   do ii = 4, 3*crystal%natom
     lst(nph2l+1)=lst(nph2l+1)+2*log(driver%phfrq(ii))
   end do
 end if

 ! Examine every wavevector of this list
 do iphl2 = 1, nph2l

   ! Initialisation of the phonon wavevector
   qphon(:,1)=dtset%qph2l(:,iphl2)
   qphnrm(1)=dtset%qnrml2(iphl2)

   !TODO: Quadrupole interactions need to be incorporated here (MR)

   ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
   ! for the second list of wv (can include non-analyticities if q /= 0)
   call dfpt_phfrq(ddb%amu, driver%displ, driver%d2cart, eigval, eigvec, crystal%indsym, &
     dtset%mpert, crystal%nsym, crystal%natom, crystal%nsym, crystal%ntypat, driver%phfrq, qphnrm(1), qphon, &
     crystal%rprimd, dtset%symdynmat, crystal%symrel, crystal%symafm, crystal%typat, crystal%ucvol)

   ! Write the phonon frequencies for the second list of wv (can include non-analyticities if q /= 0)
   call dfpt_prtph(driver%displ, dtset%eivec, dtset%enunit, ab_out, dtset%natom, driver%phfrq, qphnrm(1), qphon)
   ! TODO: Mode effective charge could be printed here for LO modes (EB)

   if (my_rank == master) then
     ! Loop is not MPI-parallelized--> no need for MPI-IO API.
     call defwrite_nonana_terms(ana_ncid, iphl2, nph2l, dtset%qph2l, dtset%natom, driver%phfrq, driver%displ, "write")
   end if

   ! Get the log of product of the square of the phonon frequencies with non-analyticities (q-->0)
   ! for the Lyddane-Sachs-Teller relation
   ! The fourth mode should have positive frequency otherwise there is an instability: LST relationship should not be evaluated
   ! Isn't it tested somewhere else (i.e. stop of the code if there are imaginary freq.)? (EB)
   if (dtset%dieflag /= 2 .and. dtset%dieflag /= 0) then
     do ii = 4, 3*crystal%natom
       lst(iphl2)=lst(iphl2)+2*log(driver%phfrq(ii))
     end do
   end if

   ! Write Raman susceptibilities for the 2nd list (can includes LO modes if q /= 0 0 0)
   if (dtset%nlflag == 1) then
     call ramansus(driver%d2cart, driver%dchide, driver%dchidt, driver%displ, dtset%mpert, crystal%natom, driver%phfrq, qphon, qphnrm(1), rsus, crystal%ucvol)
     if (my_rank == master) then
       call defwrite_nonana_raman_terms(ana_ncid, iphl2, nph2l, dtset%natom, rsus, "write")
     end if
   end if  ! nlflag = 1 (Raman suscep for the 2nd list of wv.)
 end do  ! iphl2

 ! Lyddane-Sachs-Teller relation:
 if (dtset%dieflag /= 2 .and. dtset%dieflag /= 0) then
   call ddb_diel(crystal, ddb%amu, dtset, driver%dielt_rlx, driver%displ, driver%d2cart, driver%epsinf, driver%fact_oscstr, &
     ab_out, lst, dtset%mpert, crystal%natom, nph2l, driver%phfrq, comm, ana_ncid)
 end if

 ABI_SFREE(eigval)
 ABI_SFREE(eigvec)
 ABI_SFREE(rsus)
 ABI_SFREE(lst)

end subroutine anaddb_driver_dielectric_nonana
!!***

!!****f* m_anaddb_driver/anaddb_driver_internal_strain
!! NAME
!! anaddb_driver_internal_strain
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_internal_strain(driver, dtset, ddb, asrq0)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(ddb_type), intent(in):: ddb
 type(asrq0_t), intent(in):: asrq0

!Local variables -------------------------------
 integer:: ii, iblok
 integer:: prt_internalstr
 integer:: units(2)
 integer:: rfelfd(4), rfphon(4), rfstrs(4)
 character(len = 500):: msg
 real(dp):: qphnrm(3), qphon(3, 3)

! ************************************************************************

 units = [std_out, ab_out]

 ! Here treating the internal strain tensors at Gamma point
 write(msg, '(a, a, (80a), a, a, a, a)') ch10, ('=',ii = 1, 80), ch10, ch10, &
  ' Calculation of the internal-strain  tensor',ch10
 call wrtout(units, msg)

 if (dtset%instrflag == 1) then
   call wrtout(std_out, 'instrflag = 1, so extract the internal strain constant from the 2DTE')

   ! looking after the no. of blok that contains the internal strain tensor
   qphon(:,1)=zero; qphnrm(1)=zero
   rfphon(1:2)=0; rfelfd(1:2)=0; rfstrs(1:2)=3

   call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, dtset%rfmeth)
   if (iblok == 0) then
     ABI_ERROR("DDB file must contain both uniaxial and shear strain for piezoelectric, Check your calculations")
   end if

   ! then print the internal stain tensor
   prt_internalstr = 2
   call ddb_internalstr(dtset%asr, ddb%val, asrq0%d2asr, iblok, driver%instrain, ab_out, dtset%mpert, ddb%natom, ddb%nblok, prt_internalstr)
 end if

end subroutine anaddb_driver_internal_strain
!!***

!!****f* m_anaddb_driver/anaddb_driver_elastic_tensor
!! NAME
!! anaddb_driver_elastic_tensor
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_elastic_tensor(driver, dtset, crystal, ddb, asrq0, ana_ncid)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(crystal_t), intent(in):: crystal
 type(ddb_type), intent(in):: ddb
 type(asrq0_t), intent(inout):: asrq0
 integer, intent(in):: ana_ncid

!Local variables -------------------------------
 integer:: ii, iblok, iblok_stress
 integer:: units(2)
 character(len = 500):: msg
 integer:: rfelfd(4), rfphon(4), rfstrs(4)
 real(dp):: qphnrm(3), qphon(3, 3)
 real(dp):: compl(6, 6), compl_clamped(6, 6), compl_stress(6, 6)
 real(dp):: elast_clamped(6, 6), elast_stress(6, 6)

! ************************************************************************

 units = [std_out, ab_out]

 ! here treating the elastic tensors at Gamma Point
 write(msg, '(a, a, (80a), a, a, a, a, a, a)') ch10, ('=',ii = 1, 80), ch10, ch10, &
  ' Calculation of the elastic and compliances tensor (Voigt notation)',ch10
 call wrtout(units, msg)

 call wrtout(std_out, 'so extract the elastic constant from the 2DTE')

 ! look after the blok no. that contains the stress tensor
 qphon(:,1)=zero; qphnrm(1)=zero
 rfphon(1:2)=0; rfelfd(1:2)=0; rfstrs(1:2)=0

 call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, BLKTYP_d1E_xx)
 iblok_stress = iblok

 ! look after the blok no.iblok that contains the elastic tensor
 qphon(:,1)=zero; qphnrm(1)=zero
 rfphon(1:2)=0; rfelfd(1:2)=0; rfstrs(1:2)=3

 ! for both diagonal and shear parts
 call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, dtset%rfmeth)
 if (iblok == 0) then
   ABI_ERROR("DDB file must contain both uniaxial and shear strain when elaflag != 0, Check your calculations")
 end if

 ! print the elastic tensor
 call ddb_elast(dtset, crystal, ddb%val, compl, compl_clamped, compl_stress, asrq0%d2asr, &
   driver%elast, elast_clamped, elast_stress, iblok, iblok_stress, &
   driver%instrain, ab_out, dtset%mpert, crystal%natom, ddb%nblok, ana_ncid)

end subroutine anaddb_driver_elastic_tensor
!!***

!!****f* m_anaddb_driver/anaddb_driver_piezoelectric_tensor
!! NAME
!! anaddb_driver_piezoelectric_tensor
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_piezoelectric_tensor(driver, dtset, crystal, ddb, ana_ncid)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(inout):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(crystal_t), intent(in):: crystal
 type(ddb_type), intent(in):: ddb
 integer, intent(in):: ana_ncid

!Local variables -------------------------------
 integer:: ii, iblok
 integer:: units(2)
 character(len = 500):: msg
 integer:: rfelfd(4), rfphon(4), rfstrs(4)
 real(dp):: qphnrm(3), qphon(3, 3)
 real(dp):: piezo(6, 3)

! ************************************************************************

 units = [std_out, ab_out]

 ! Here treating the piezoelectric tensor at Gamma Point
 write(msg, '(a, a, (80a), a, a, a, a, a)') ch10, ('=',ii = 1, 80), ch10, ch10, &
 ' Calculation of the tensor related to piezoelectric effetc',ch10, &
 '  (Elastic indices in Voigt notation)',ch10
 call wrtout(units, msg)

 call wrtout(std_out, 'extract the piezoelectric constant from the 2DTE')

 ! Looking for the gamma point block
 qphon(:,1)=zero; qphnrm(1)=zero
 rfphon(1:2)=0; rfelfd(1:2)=0; rfstrs(1:2)=3

 ! For both diagonal and shear parts
 call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, dtset%rfmeth)
 if (iblok == 0) then
   ABI_ERROR("DDB file must contain both uniaxial and shear strain for piezoelectric, Check your calculations")
 end if

 ! Then print out the piezoelectric constants
 call ddb_piezo(dtset, ddb%val, driver%dielt_rlx, driver%elast, iblok, &
     & driver%instrain, ab_out, dtset%mpert, crystal%natom, ddb%nblok, piezo, &
     & crystal%ucvol, ana_ncid)

end subroutine anaddb_driver_piezoelectric_tensor
!!***

!!****f* m_anaddb_driver/anaddb_driver_flexoelectric_tensor
!! NAME
!! anaddb_driver_flexoelectric_tensor
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_flexoelectric_tensor(driver, dtset, crystal, ddb, ddb_lw, ddb_hdr, asrq0)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(in):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(crystal_t), intent(in):: crystal
 type(ddb_type), intent(in):: ddb, ddb_lw
 type(ddb_hdr_type), intent(in):: ddb_hdr
 type(asrq0_t), intent(in):: asrq0

!Local variables -------------------------------
 integer:: ii
 integer:: units(2)
 character(len = 500):: msg

! ************************************************************************

 units = [std_out, ab_out]

 ! Here treating the flexoelectric tensor
 write(msg, '(a, a, (80a), a, a, a, a)') ch10, ('=',ii = 1, 80), ch10, ch10, &
 ' Calculation of the tensors related to flexoelectric effect',ch10
 call wrtout(units, msg)

 ! Compute and print the contributions to the flexoelectric tensor
 call ddb_flexo(dtset%asr, asrq0%d2asr, ddb, ddb_lw, ddb_hdr%ddb_version, crystal, &
     & dtset%filename_ddb, dtset%flexoflag, dtset%prtvol, driver%zeff)

end subroutine anaddb_driver_flexoelectric_tensor
!!***

!!****f* m_anaddb_driver/anaddb_driver_lattice_wannier
!! NAME
!! anaddb_driver_lattice_wannier
!!
!! FUNCTION
!! Lattice Wannier function calculation.
!! Compute the Dynamical matrix for a dense Q-mesh
!! Input the eigenvectors and eigenvalues to the Lattcie Wannier module
!! Construct the Lattice Wannier functions
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_driver_lattice_wannier(driver, dtset, crystal, ifc, comm)

!Arguments -------------------------------
 class(anaddb_driver_type), intent(in):: driver
 type(anaddb_dataset_type), intent(in):: dtset
 type(crystal_t),intent(in):: crystal
 type(ifc_type), intent(in):: ifc
 integer,intent(in):: comm

!Local variables -------------------------------
 integer:: ii
 character(len = 500):: msg
 integer:: units(2)

! ************************************************************************

 units = [std_out, ab_out]

 write(msg, '(a, (80a), 4a)')ch10, ('=',ii = 1, 80), ch10, ch10, ' Calculation of lattice Wannier functions ',ch10
 call wrtout(units, msg)
 call run_lattice_wannier(ifc=ifc, crystal=crystal, dtset=dtset, prefix=dtset%prefix_outdata, comm=comm)
 write(msg, '(a, (80a))')ch10, ('=',ii = 1, 80)
 call wrtout(units, msg)

end subroutine anaddb_driver_lattice_wannier
!!***

end module m_anaddb_driver
!!***
