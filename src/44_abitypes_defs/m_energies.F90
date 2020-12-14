!!****m* ABINIT/m_energies
!! NAME
!!  m_energies
!!
!! FUNCTION
!!  This module provides the definition of the energies used
!!  to store energies from GS calculations.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MT, DC)
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

module m_energies

 use defs_basis
 use m_abicore
 use m_errors
 use m_nctk
 use m_dtset
#ifdef HAVE_NETCDF
 use netcdf
#endif

 implicit none

 private

!public parameter
 integer, public, parameter :: n_energies=35
  ! Number of energies stored in energies datastructure

!!***

!!****t* m_energies/energies_type
!! NAME
!! energies_type
!!
!! FUNCTION
!! This structured datatype contains all parts of total energy. Not all
!! attributes may have a value, depending on the scheme used to
!! compute the total energy and several options read from dtset.
!!
!! SOURCE

 type, public :: energies_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  real(dp) :: e_chempot=zero
   ! energy from spatially-varying chemical potential

  real(dp) :: e_constrained_dft=zero
   ! correction to energy from constrained dft, to make it variational

  real(dp) :: e_corepsp=zero
   ! psp core-core energy

  real(dp) :: e_corepspdc=zero
   ! psp core-core energy double-counting

  real(dp) :: e_eigenvalues=zero
   ! Sum of the eigenvalues - Band energy (Hartree)
   ! (valid for double-counting scheme dtset%optene == 1)

  real(dp) :: e_entropy=zero
   ! Entropy energy due to the occupation number smearing (if metal)
   ! Value is multiplied by dtset%tsmear, see %entropy for the entropy alone.
   ! (valid for metals, dtset%occopt>=3 .and. dtset%occopt<=8)

  real(dp) :: entropy=zero

  real(dp) :: e_elecfield=zero
   ! Electric enthalpy, by adding both ionic and electronic contributions

  real(dp) :: e_electronpositron=zero
   ! Electron-positron: electron-positron interaction energy

  real(dp) :: edc_electronpositron=zero
   ! Electron-positron: double-counting electron-positron interaction energy

  real(dp) :: e0_electronpositron=zero
   !  Electron-positron: energy only due to unchanged particles
   !                     (if calctype=1, energy due to electrons only)
   !                     (if calctype=2, energy due to positron only)

  real(dp) :: e_exactX=zero
   ! Fock exact-exchange energy (hartree)

  real(dp) :: e_ewald=zero
   ! Ewald energy (hartree), store also the ion/ion energy for free boundary conditions.

  real(dp) :: e_fermie=zero
   ! Fermie energy

  real(dp) :: e_fock=zero
   ! Fock part of total energy (hartree units)

  real(dp) :: e_fockdc=zero
   ! Fock part of energy double counting (hartree units)

  real(dp) :: e_fock0=zero
   ! Fock part of total energy, evaluated with the frozen Fock operator (usually in the case of the ACE)

  real(dp) :: e_hartree=zero
   ! Hartree part of total energy (hartree units)

  real(dp) :: e_hybcomp_E0=zero
   ! First compensation energy in the case of hybrid functionals, due to the use of two different XC functionals
   ! Term related to energy, at frozen density

  real(dp) :: e_hybcomp_v0=zero
   ! Second compensation energy in the case of hybrid functionals, due to the use of two different XC potentials
   ! Term related to potential, at frozen density

  real(dp) :: e_hybcomp_v=zero
   ! Third compensation energy in the case of hybrid functionals, due to the use of two different XC potentials
   ! Term related to potential, at optimized density

  real(dp) :: e_kinetic=zero
   ! Kinetic energy part of total energy.
   ! (valid for direct scheme, dtset%optene == 0)

  real(dp) :: e_localpsp=zero
   ! Local psp energy (hartree)

  real(dp) :: e_magfield=zero
   ! Orbital magnetic enthalpy, by adding orbital contribution

  real(dp) :: e_monopole=zero
   ! Monopole correction to the total energy for charged supercells

  real(dp) :: e_nlpsp_vfock=zero
   ! Nonlocal pseudopotential part of total energy.

  real(dp) :: e_paw=zero
   ! PAW spherical part energy

  real(dp) :: e_pawdc=zero
   ! PAW spherical part double-counting energy

  real(dp) :: e_sicdc=zero
   ! Self-interaction energy double-counting

  real(dp) :: e_vdw_dftd=zero
   ! Dispersion energy from DFT-D Van der Waals correction (hartree)

  real(dp) :: e_xc=zero
   ! Exchange-correlation energy (hartree)

  real(dp) :: e_xcdc=zero
   ! enxcdc=exchange-correlation double-counting energy (hartree)

  real(dp) :: e_xc_vdw=zero
   ! vdW-DF correction to the XC energy

  real(dp) :: h0=zero
   ! h0=e_kinetic+e_localpsp+e_nlpsp_vfock

  real(dp) :: e_zeeman=zero
   ! Zeeman spin times magnetic field contribution to the XC energy

 end type energies_type

!public procedures.
 public :: energies_init
 public :: energies_copy
 public :: energies_to_array
 public :: energies_eval_eint
 public :: energies_ncwrite
!!***


CONTAINS !===========================================================
!!***

!!****f* m_energies/energies_init
!!
!! NAME
!! energies_init
!!
!! FUNCTION
!! Set zero in all values of a type(energies_type) object
!!
!! INPUTS
!!
!! OUTPUT
!!   energies <type(energies_type)>=values to initialise
!!
!! PARENTS
!!      m_bethe_salpeter,m_electronpositron,m_gstate,m_positron,m_results_gs
!!      m_scfcv_core,m_screening_driver,m_sigma_driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine energies_init(energies)

!Arguments ------------------------------------
!scalars
 type(energies_type),intent(out) :: energies

! *************************************************************************

!@energies_type

 energies%e_chempot     = zero
 energies%e_constrained_dft    = zero
 energies%e_corepsp     = zero
 energies%e_corepspdc   = zero
 energies%e_eigenvalues = zero
 energies%e_elecfield   = zero
 energies%e_electronpositron   = zero
 energies%edc_electronpositron = zero
 energies%e0_electronpositron  = zero
 energies%e_entropy     = zero
 energies%e_exactX      = zero
 energies%entropy       = zero
 energies%e_ewald       = zero
 energies%e_fermie      = zero
 energies%e_fock        = zero
 energies%e_fockdc      = zero
 energies%e_fock0       = zero
 energies%e_hartree     = zero
 energies%e_hybcomp_E0  = zero
 energies%e_hybcomp_v0  = zero
 energies%e_hybcomp_v   = zero
 energies%e_kinetic     = zero
 energies%e_localpsp    = zero
 energies%e_magfield    = zero
 energies%e_monopole    = zero
 energies%e_nlpsp_vfock = zero
 energies%e_paw         = zero
 energies%e_pawdc       = zero
 energies%e_sicdc       = zero
 energies%e_vdw_dftd    = zero
 energies%e_xc          = zero
 energies%e_xcdc        = zero
 energies%e_xc_vdw      = zero
 energies%h0            = zero
 energies%e_zeeman      = zero

end subroutine energies_init
!!***

!----------------------------------------------------------------------

!!****f* m_energies/energies_copy
!!
!! NAME
!! energies_copy
!!
!! FUNCTION
!! Copy a type(energies_type) object into another
!!
!! INPUTS
!!   energies_in <type(energies_type)>=input values (to copy)
!!
!! OUTPUT
!!   energies_out <type(energies_type)>=output values
!!
!! PARENTS
!!      m_afterscfloop,m_electronpositron,m_positron,m_results_gs
!!
!! CHILDREN
!!
!! SOURCE

 subroutine energies_copy(energies_in,energies_out)

!Arguments ------------------------------------
!scalars
 type(energies_type),intent(in)  :: energies_in
 type(energies_type),intent(out) :: energies_out

!*************************************************************************

!@energies_type

 energies_out%e_chempot            = energies_in%e_chempot
 energies_out%e_constrained_dft    = energies_in%e_constrained_dft
 energies_out%e_corepsp            = energies_in%e_corepsp
 energies_out%e_corepspdc          = energies_in%e_corepspdc
 energies_out%e_eigenvalues        = energies_in%e_eigenvalues
 energies_out%e_elecfield          = energies_in%e_elecfield
 energies_out%e_electronpositron   = energies_in%e_electronpositron
 energies_out%edc_electronpositron = energies_in%edc_electronpositron
 energies_out%e0_electronpositron  = energies_in%e0_electronpositron
 energies_out%entropy              = energies_in%entropy
 energies_out%e_entropy            = energies_in%e_entropy
 energies_out%e_ewald              = energies_in%e_ewald
 energies_out%e_exactX             = energies_in%e_exactX
 energies_out%e_fermie             = energies_in%e_fermie
 energies_out%e_fock               = energies_in%e_fock
 energies_out%e_fockdc             = energies_in%e_fockdc
 energies_out%e_fock0              = energies_in%e_fock0
 energies_out%e_hartree            = energies_in%e_hartree
 energies_out%e_hybcomp_E0         = energies_in%e_hybcomp_E0
 energies_out%e_hybcomp_v0         = energies_in%e_hybcomp_v0
 energies_out%e_hybcomp_v          = energies_in%e_hybcomp_v
 energies_out%e_kinetic            = energies_in%e_kinetic
 energies_out%e_localpsp           = energies_in%e_localpsp
 energies_out%e_magfield           = energies_in%e_magfield
 energies_out%e_monopole           = energies_in%e_monopole
 energies_out%e_nlpsp_vfock        = energies_in%e_nlpsp_vfock
 energies_out%e_paw                = energies_in%e_paw
 energies_out%e_pawdc              = energies_in%e_pawdc
 energies_out%e_sicdc              = energies_in%e_sicdc
 energies_out%e_vdw_dftd           = energies_in%e_vdw_dftd
 energies_out%e_xc                 = energies_in%e_xc
 energies_out%e_xcdc               = energies_in%e_xcdc
 energies_out%e_xc_vdw             = energies_in%e_xc_vdw
 energies_out%h0                   = energies_in%h0
 energies_out%e_zeeman             = energies_in%e_zeeman

end subroutine energies_copy
!!***

!----------------------------------------------------------------------

!!****f* m_energies/energies_to_array
!!
!! NAME
!! energies_to_array
!!
!! FUNCTION
!! Transfer a energies datastructure into a single array or
!! transfer an array into a energies datastructure
!!
!! INPUTS
!!   option= 1: copy energies datastructure into an array
!!   option=-1: copy an array into a energies datastructure
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!   energies <type(energies_type)>=energies stored in a datastructure
!!   energies_array=energies stored in a single array
!!
!! PARENTS
!!      m_results_img
!!
!! CHILDREN
!!
!! SOURCE

 subroutine energies_to_array(energies,energies_array,option)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option
!arrays
 real(dp),intent(inout) :: energies_array(n_energies)
 type(energies_type),intent(inout)  :: energies

!*************************************************************************

!@energies_type

 if (option==1) then
   energies_array(1)=energies%e_chempot
   energies_array(2)=energies%e_constrained_dft
   energies_array(3)=energies%e_corepsp
   energies_array(4)=energies%e_corepspdc
   energies_array(5)=energies%e_eigenvalues
   energies_array(6)=energies%e_elecfield
   energies_array(7)=energies%e_electronpositron
   energies_array(8)=energies%edc_electronpositron
   energies_array(9)=energies%e0_electronpositron
   energies_array(10)=energies%entropy
   energies_array(11)=energies%e_entropy
   energies_array(12)=energies%e_ewald
   energies_array(13)=energies%e_exactX
   energies_array(14)=energies%e_fermie
   energies_array(15)=energies%e_fock
   energies_array(16)=energies%e_fockdc
   energies_array(17)=energies%e_fock0
   energies_array(18)=energies%e_hartree
   energies_array(19)=energies%e_hybcomp_E0
   energies_array(20)=energies%e_hybcomp_v0
   energies_array(21)=energies%e_hybcomp_v
   energies_array(22)=energies%e_kinetic
   energies_array(23)=energies%e_localpsp
   energies_array(24)=energies%e_magfield
   energies_array(25)=energies%e_monopole
   energies_array(26)=energies%e_nlpsp_vfock
   energies_array(27)=energies%e_paw
   energies_array(28)=energies%e_pawdc
   energies_array(29)=energies%e_sicdc
   energies_array(30)=energies%e_vdw_dftd
   energies_array(31)=energies%e_xc
   energies_array(32)=energies%e_xcdc
   energies_array(33)=energies%e_xc_vdw
   energies_array(34)=energies%h0
   energies_array(35)=energies%e_zeeman
 end if

 if (option==-1) then
   energies%e_chempot            = energies_array(1)
   energies%e_constrained_dft    = energies_array(2)
   energies%e_corepsp            = energies_array(3)
   energies%e_corepspdc          = energies_array(4)
   energies%e_eigenvalues        = energies_array(5)
   energies%e_elecfield          = energies_array(6)
   energies%e_electronpositron   = energies_array(7)
   energies%edc_electronpositron = energies_array(8)
   energies%e0_electronpositron  = energies_array(9)
   energies%entropy              = energies_array(10)
   energies%e_entropy            = energies_array(11)
   energies%e_ewald              = energies_array(12)
   energies%e_exactX             = energies_array(13)
   energies%e_fermie             = energies_array(14)
   energies%e_fock               = energies_array(15)
   energies%e_fockdc             = energies_array(16)
   energies%e_fock0              = energies_array(17)
   energies%e_hartree            = energies_array(18)
   energies%e_hybcomp_E0         = energies_array(19)
   energies%e_hybcomp_v0         = energies_array(20)
   energies%e_hybcomp_v          = energies_array(21)
   energies%e_kinetic            = energies_array(22)
   energies%e_localpsp           = energies_array(23)
   energies%e_magfield           = energies_array(24)
   energies%e_monopole           = energies_array(25)
   energies%e_nlpsp_vfock        = energies_array(26)
   energies%e_paw                = energies_array(27)
   energies%e_pawdc              = energies_array(28)
   energies%e_sicdc              = energies_array(29)
   energies%e_vdw_dftd           = energies_array(30)
   energies%e_xc                 = energies_array(31)
   energies%e_xcdc               = energies_array(32)
   energies%e_xc_vdw             = energies_array(33)
   energies%h0                   = energies_array(34)
   energies%e_zeeman             = energies_array(35)
 end if

end subroutine energies_to_array
!!***

!----------------------------------------------------------------------

!!****f* m_energies/energies_eval_eint
!!
!! NAME
!! energies_eval_eint
!!
!! FUNCTION
!! Compute the internal energy (Direct and DC as it was in prtene)
!!
!! INPUTS
!!  energies <type(energies_type)>=values of parts of total energy
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryphase
!!   | kptopt
!!   | occopt
!!   | positron=option for electron-positron calculation
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  usewvl= 0 for PW calculation; =1 for WVL calculation
!!
!!
!! OUTPUT
!!  optdc=option for double counting scheme
!!  eint=internal energy with direct scheme
!!  eintdc=internal energy with double counting scheme
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_common,m_entropyDMFT
!!
!! CHILDREN
!!
!! SOURCE

 subroutine energies_eval_eint(energies,dtset,usepaw,optdc,eint,eintdc)

!Arguments ------------------------------------
!scalars
 type(energies_type),intent(in) :: energies
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: usepaw
 integer , intent(out) :: optdc
 real(dp), intent(out) :: eint
 real(dp), intent(out) :: eintdc

!Local variables-------------------------------
! Do not modify the length of this string
!scalars
 logical :: positron
 logical :: wvlbigdft=.false.

! *************************************************************************

!If usewvl: wvlbigdft indicates that the BigDFT workflow will be followed
 wvlbigdft=(dtset%usewvl==1.and.dtset%wvl_bigdft_comp==1)

 optdc=-1;positron=(dtset%positron/=0)
 if (.not.positron) then
   if ((abs(energies%e_xcdc)<1.e-15_dp).and.(abs(energies%e_fockdc)<1.e-15_dp)) optdc=0
   if (abs(energies%e_localpsp)<1.e-15_dp.and.(abs(energies%e_xcdc)>1.e-15_dp.or.abs(energies%e_fockdc)>1.e-15_dp)) optdc=1
   if (abs(energies%e_localpsp)>1.e-15_dp.and.(abs(energies%e_xcdc)>1.e-15_dp.or.abs(energies%e_fockdc)>1.e-15_dp)) optdc=2
   if (wvlbigdft .and. dtset%iscf > 0) optdc=1
 else
   if (abs(energies%edc_electronpositron)<1.e-15_dp) optdc=0
   if (abs(energies%e_electronpositron)<1.e-15_dp.and.abs(energies%edc_electronpositron)>1.e-15_dp) optdc=1
   if (abs(energies%e_electronpositron)>1.e-15_dp.and.abs(energies%edc_electronpositron)>1.e-15_dp) optdc=2
 end if

 eint  = zero
 eintdc = zero

!============= Evaluate some parts of the energy ===========

 if (optdc==0.or.optdc==2) then
   eint = energies%e_kinetic + energies%e_hartree + energies%e_xc+ &
!&  +two*energies%e_fock-energies%e_fock0+&  ! The Fock energy is already included in the non_local part ...
!&  energies%e_nlpsp_vfock - energies%e_fock0+&
&   energies%e_hybcomp_E0 -energies%e_hybcomp_v0 + energies%e_hybcomp_v+&
&   energies%e_localpsp + energies%e_corepsp + energies%e_constrained_dft

!  See similar section in m_scfcv_core.F90
!  XG 20181025 This gives a variational energy in case of NCPP with all bands occupied - not yet for metals.
   if (usepaw==0) eint = eint + energies%e_nlpsp_vfock - energies%e_fock0
!  XG 20181025 I was expecting the following to give also a variational energy in case of PAW, but this is not true.
!  if (usepaw==1) eint = eint + energies%e_paw + energies%e_nlpsp_vfock - energies%e_fock0
!  XG 20181025 So, the following is giving a non-variational expression ...
   if (usepaw==1) eint = eint + energies%e_paw + energies%e_fock

   if (dtset%berryopt==4 .or. dtset%berryopt==6 .or. dtset%berryopt==7 .or.  &
&   dtset%berryopt==14 .or. dtset%berryopt==16 .or. dtset%berryopt==17) eint=eint+energies%e_elecfield    !!HONG
   eint = eint + energies%e_ewald + energies%e_chempot + energies%e_vdw_dftd
   if (positron) eint=eint+energies%e0_electronpositron+energies%e_electronpositron
 end if
 if (optdc>=1) then
   eintdc = energies%e_eigenvalues - energies%e_hartree + energies%e_xc &
&   + energies%e_hybcomp_E0 - energies%e_hybcomp_v0 &
&   - energies%e_xcdc + energies%e_corepsp - energies%e_corepspdc-energies%e_fock0
   if (usepaw==1) eintdc = eintdc + energies%e_pawdc
   if (dtset%berryopt==4 .or. dtset%berryopt==6 .or. dtset%berryopt==7 .or.  &
&   dtset%berryopt==14 .or. dtset%berryopt==16 .or. dtset%berryopt==17) eintdc = eintdc + energies%e_elecfield
   eintdc = eintdc + energies%e_ewald + energies%e_chempot + energies%e_vdw_dftd + energies%e_constrained_dft
   if (positron) eintdc=eintdc-energies%edc_electronpositron &
&   +energies%e0_electronpositron+energies%e_electronpositron
 end if

end subroutine energies_eval_eint
!!***

!----------------------------------------------------------------------

!!****f* m_energies/energies_ncwrite
!! NAME
!! energies_ncwrite
!!
!! FUNCTION
!!  Write the contenc of the datatype in a netcdf file.
!!
!! INPUTS
!!  ncid=NC file handle
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_results_gs
!!
!! CHILDREN
!!
!! SOURCE

subroutine energies_ncwrite(enes, ncid)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 type(energies_type),intent(in) :: enes

!Local variables-------------------------------
!scalars
#ifdef HAVE_NETCDF
 integer :: ncerr

! *************************************************************************

!@energies_type
 ncerr = nctk_defnwrite_dpvars(ncid, [character(len=nctk_slen) :: &
  "e_chempot", "e_constrained_dft", "e_corepsp", "e_corepspdc", "e_eigenvalues", "e_elecfield", &
  "e_electronpositron", "edc_electronpositron", "e0_electronpositron",&
  "e_entropy", "entropy", "e_ewald", &
  "e_exactX","e_fermie", &
  "e_fock", "e_fockdc", "e_fock0", "e_hartree", "e_hybcomp_E0", "e_hybcomp_v0", "e_hybcomp_v", "e_kinetic",&
  "e_localpsp", "e_magfield", "e_monopole", "e_nlpsp_vfock", &
  "e_paw", "e_pawdc", "e_sicdc", "e_vdw_dftd", &
  "e_xc", "e_xcdc", "e_xc_vdw", &
  "h0","e_zeeman"], &
  [enes%e_chempot, enes%e_constrained_dft, enes%e_corepsp, enes%e_corepspdc, enes%e_eigenvalues, enes%e_elecfield, &
   enes%e_electronpositron, enes%edc_electronpositron, enes%e0_electronpositron,&
   enes%e_entropy, enes%entropy, enes%e_ewald, &
   enes%e_exactX, enes%e_fermie, &
   enes%e_fock, enes%e_fockdc,enes%e_fock0,  enes%e_hartree, &
   enes%e_hybcomp_E0, enes%e_hybcomp_v0, enes%e_hybcomp_v, enes%e_kinetic,&
   enes%e_localpsp, enes%e_magfield, enes%e_monopole, enes%e_nlpsp_vfock, &
   enes%e_paw, enes%e_pawdc, enes%e_sicdc, enes%e_vdw_dftd,&
   enes%e_xc, enes%e_xcdc, enes%e_xc_vdw,&
   enes%h0,enes%e_zeeman])

 NCF_CHECK(ncerr)

#else
 MSG_ERROR("ETSF-IO support is not activated.")
#endif

end subroutine energies_ncwrite
!!***

end module m_energies
!!***
