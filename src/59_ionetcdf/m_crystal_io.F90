!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_crystal_io
!! NAME
!! m_crystal_io
!!
!! FUNCTION
!! Module containing the methods used to do IO on crystal objects.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (MG, YP, DC)
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

MODULE m_crystal_io

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_crystal
 use m_atomdata
 use m_nctk
 use iso_c_binding
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_fstrings,       only : sjoin, yesno
 use m_io_tools,       only : file_exists
 use defs_abitypes,    only : hdr_type

 implicit none

 private
!!***

 public :: crystal_from_hdr      ! Initialize the object from the abinit header.
 !public :: crystal_from_dtset    ! Initialize the object from the abinit dataset.
 public :: crystal_ncwrite       ! Dump the object in a netcdf file associated to a ncid.
 public :: crystal_ncwrite_path  ! Dump the object to file.


CONTAINS

!!****f* m_crystal_io/crystal_from_hdr
!! NAME
!!  crystal_from_hdr
!!
!! FUNCTION
!!  Initializes a crystal_t data type starting from the abinit header.
!!
!! INPUTS
!!  hdr<hdr_type>=the abinit header
!!  timrev ==2 => take advantage of time-reversal symmetry
!!         ==1 ==> do not use time-reversal symmetry
!!  remove_inv [optional]= if .TRUE. the inversion symmetry is removed from the set of operations
!!  even if it is present in the header
!!
!! OUTPUT
!!  cryst<crystal_t>= the data type filled with data reported in the abinit header
!!
!! TODO
!!  Add information on the use of time-reversal in the Abinit header.
!!
!! PARENTS
!!      cut3d,eph,fold2Bloch,gstate,m_ddk,m_dvdb,m_ioarr,m_iowf,m_wfd,m_wfk
!!      mlwfovlp_qp,mrgscr,setup_bse,setup_screening,setup_sigma,wfk_analyze
!!
!! CHILDREN
!!      atomdata_from_znucl,crystal_init
!!
!! SOURCE

subroutine crystal_from_hdr(cryst,hdr,timrev,remove_inv)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'crystal_from_hdr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(hdr_type),intent(in) :: hdr
 type(crystal_t),intent(out) :: cryst
 integer,intent(in) :: timrev
 logical,optional,intent(in) :: remove_inv

!Local variables-------------------------------
 integer :: space_group
 logical :: rinv,use_antiferro
! *********************************************************************

 rinv=.FALSE.; if (PRESENT(remove_inv)) rinv=remove_inv
 use_antiferro = (hdr%nspden==2.and.hdr%nsppol==1)

 ! Consistency check
 ABI_CHECK(any(timrev == [1, 2]),"timrev should be in (1|2)")
 if (use_antiferro) then
   ABI_CHECK(ANY(hdr%symafm==-1),"Wrong nspden, nsppol, symafm.")
 end if

 space_group=0 !FIXME not known

 call crystal_init(hdr%amu,cryst,space_group,hdr%natom,hdr%npsp,hdr%ntypat,hdr%nsym,hdr%rprimd,hdr%typat,hdr%xred,&
& hdr%zionpsp,hdr%znuclpsp,timrev,use_antiferro,rinv,hdr%title,&
& symrel=hdr%symrel,tnons=hdr%tnons,symafm=hdr%symafm) ! Optional

end subroutine crystal_from_hdr
!!***

!----------------------------------------------------------------------

!!****f* m_crystal_io/crystal_ncwrite
!! NAME
!! crystal_ncwrite
!!
!! FUNCTION
!! Output system geometry to a file, using the NETCDF file format and ETSF I/O.
!! Data are taken from the crystal_t object.
!!
!! INPUTS
!!  cryst<crystal_t>=Object defining the unit cell and its symmetries.
!!  ncid=NC file handle.
!!
!! OUTPUT
!!  Only writing
!!
!! NOTES
!!  Alchemy not treated, since crystal should be initialized at the beginning of the run.
!!
!! PARENTS
!!      anaddb,eig2tot,exc_spectra,dfpt_looppert,m_haydock,m_phonons,m_shirley
!!      outscfcv,sigma
!!
!! CHILDREN
!!      atomdata_from_znucl
!!
!! SOURCE

integer function crystal_ncwrite(cryst, ncid) result(ncerr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'crystal_ncwrite'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 type(crystal_t),intent(in) :: cryst

#ifdef HAVE_NETCDF
!Local variables-------------------------------
!scalars
 integer :: itypat
 character(len=500) :: msg
 character(len=etsfio_charlen) :: symmorphic
 type(atomdata_t) :: atom
!arrays
 character(len=2) :: symbols(cryst%ntypat)
 character(len=80) :: psp_desc(cryst%ntypat),symbols_long(cryst%ntypat)

! *************************************************************************

 ! @crystal_t

 ! TODO alchemy not treated correctly
 if (isalchemical(cryst)) then
   write(msg,"(3a)")&
    "Alchemical crystals are not fully supported by the netcdf format",ch10,&
    "Important parameters (e.g. znucl, symbols) are not written with the correct value"
   MSG_WARNING(msg)
 end if

 symmorphic = yesno(isymmorphic(cryst))

 ! Define dimensions.
 ncerr = nctk_def_dims(ncid, [ &
   nctkdim_t("complex", 2), nctkdim_t("symbol_length", 2),&
   nctkdim_t("character_string_length", 80), nctkdim_t("number_of_cartesian_directions", 3),&
   nctkdim_t("number_of_reduced_dimensions", 3), nctkdim_t("number_of_vectors", 3),&
   nctkdim_t("number_of_atoms", cryst%natom), nctkdim_t("number_of_atom_species", cryst%ntypat),&
   nctkdim_t("number_of_symmetry_operations", cryst%nsym)], defmode=.True.)
 NCF_CHECK(ncerr)

 ! Define variables
 NCF_CHECK(nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "space_group"]))

 ncerr = nctk_def_arrays(ncid, [ &
  ! Atomic structure and symmetry operations
  nctkarr_t("primitive_vectors", "dp", "number_of_cartesian_directions, number_of_vectors"), &
  nctkarr_t("reduced_symmetry_matrices", "int", &
    "number_of_reduced_dimensions, number_of_reduced_dimensions, number_of_symmetry_operations"), &
  nctkarr_t("reduced_symmetry_translations", "dp", "number_of_reduced_dimensions, number_of_symmetry_operations"), &
  nctkarr_t("atom_species", "int", "number_of_atoms"), &
  nctkarr_t("reduced_atom_positions", "dp", "number_of_reduced_dimensions, number_of_atoms"), &
  nctkarr_t("atomic_numbers", "dp", "number_of_atom_species"), &
  nctkarr_t("atom_species_names", "char", "character_string_length, number_of_atom_species"), &
  nctkarr_t("chemical_symbols", "char", "symbol_length, number_of_atom_species"), &
  nctkarr_t('atomic_mass_units', "dp", "number_of_atom_species"), &
  ! Atomic information.
  nctkarr_t("valence_charges", "dp", "number_of_atom_species"), &  ! NB: This variable is not written if alchemical
  nctkarr_t("pseudopotential_types", "char", "character_string_length, number_of_atom_species") &
 ])
 NCF_CHECK(ncerr)

 ! Some variables require the "symmorphic" attribute.
 NCF_CHECK(nf90_put_att(ncid, vid("reduced_symmetry_matrices"), "symmorphic", symmorphic))
 NCF_CHECK(nf90_put_att(ncid, vid("reduced_symmetry_translations"), "symmorphic", symmorphic))

 ! At this point we have an ETSF-compliant file. Add additional data for internal use in abinit.
 ! TODO add spinat.
 ncerr = nctk_def_arrays(ncid, nctkarr_t('symafm', "int", "number_of_symmetry_operations"))
 NCF_CHECK(ncerr)

 ! Set-up atomic symbols.
 do itypat=1,cryst%ntypat
   call atomdata_from_znucl(atom,cryst%znucl(itypat))
   symbols(itypat) = atom%symbol
   write(symbols_long(itypat),'(a2,a78)') symbols(itypat),REPEAT(CHAR(0),78)
   write(psp_desc(itypat),'(2a)') &
&    cryst%title(itypat)(1:MIN(80,LEN_TRIM(cryst%title(itypat)))),REPEAT(CHAR(0),MAX(0,80-LEN_TRIM(cryst%title(itypat))))
 end do

 ! Write data.
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid("space_group"), cryst%space_group))
 NCF_CHECK(nf90_put_var(ncid, vid("primitive_vectors"), cryst%rprimd))
 NCF_CHECK(nf90_put_var(ncid, vid("reduced_symmetry_matrices"), cryst%symrel))
 NCF_CHECK(nf90_put_var(ncid, vid("reduced_symmetry_translations"), cryst%tnons))
 NCF_CHECK(nf90_put_var(ncid, vid("atom_species"), cryst%typat))
 NCF_CHECK(nf90_put_var(ncid, vid("reduced_atom_positions"), cryst%xred))
 NCF_CHECK(nf90_put_var(ncid, vid("atomic_numbers"), cryst%znucl(1:cryst%ntypat)))
 NCF_CHECK(nf90_put_var(ncid, vid("atom_species_names"), symbols_long))
 NCF_CHECK(nf90_put_var(ncid, vid("chemical_symbols"), symbols))
 NCF_CHECK(nf90_put_var(ncid, vid('atomic_mass_units'), cryst%amu))
 NCF_CHECK(nf90_put_var(ncid, vid("pseudopotential_types"), psp_desc))
 if (cryst%npsp == cryst%ntypat) then
   NCF_CHECK(nf90_put_var(ncid, vid("valence_charges"), cryst%zion))
 end if

 NCF_CHECK(nf90_put_var(ncid, vid("symafm"), cryst%symafm))

#else
 MSG_ERROR("netcdf library not available")
#endif

contains
 integer function vid(vname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vid'
!End of the abilint section

   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end function crystal_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_crystal_io/crystal_ncwrite_path
!! NAME
!! crystal_ncwrite_path
!!
!! FUNCTION
!! Output system geometry to a file, using the NETCDF file format and ETSF I/O.
!!
!! INPUTS
!!  crystal<crystal_t>=Object defining the unit cell and its symmetries.
!!  path=filename
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function crystal_ncwrite_path(crystal, path) result(ncerr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'crystal_ncwrite_path'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 type(crystal_t),intent(in) :: crystal

#ifdef HAVE_NETCDF
!Local variables-------------------------------
!scalars
 integer :: ncid

! *************************************************************************

 ncerr = nf90_noerr
 if (file_exists(path)) then
   NCF_CHECK(nctk_open_modify(ncid, path, xmpi_comm_self))
 else
   ncerr = nctk_open_create(ncid, path, xmpi_comm_self)
   NCF_CHECK_MSG(ncerr, sjoin("creating:", path))
 end if

 NCF_CHECK(crystal_ncwrite(crystal, ncid))
 NCF_CHECK(nf90_close(ncid))
#endif

end function crystal_ncwrite_path
!!***

end module m_crystal_io
!!***
