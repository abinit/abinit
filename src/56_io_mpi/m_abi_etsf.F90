!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abi_etsf
!! NAME
!! m_abi_etsf
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2019 ABINIT group (DCA,YP,MJV,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!  Remove this module. At present it's used only in m_iowf.F90
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_abi_etsf

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
 use m_abicore
 use m_dtset
 use m_xmpi
 use m_errors
 use m_atomdata
 use m_nctk
 use iso_c_binding
#ifdef HAVE_ETSF_IO
 use etsf_io
#endif

 use m_fstrings,   only : endswith

 implicit none

 private

#ifdef HAVE_ETSF_IO
 public :: abi_etsf_dims_init
#endif
 public :: abi_etsf_init

CONTAINS  !===========================================================
!!***

!!****f* m_abi_etsf/etsf_dims_nullify
!! NAME
!! etsf_dims_nullify
!!
!! FUNCTION
!!  Initialize the structure defining ETSF dimensions.
!!
!! INPUTS
!!
!! OUTPUT
!!  dims=structure with ETSF dimensions.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#ifdef HAVE_ETSF_IO

subroutine etsf_dims_nullify(Dims,dimlen)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: dimlen
 type(etsf_dims),intent(inout) :: Dims

!Local variables-------------------------------
 integer :: my_dimlen

! *************************************************************************

  my_dimlen = 3; if (PRESENT(dimlen)) my_dimlen = dimlen

  Dims%max_number_of_angular_momenta                 = etsf_no_dimension
  Dims%max_number_of_basis_grid_points               = etsf_no_dimension
  Dims%max_number_of_coefficients                    = etsf_no_dimension
  Dims%max_number_of_projectors                      = etsf_no_dimension
  Dims%max_number_of_states                          = etsf_no_dimension
  Dims%number_of_atoms                               = etsf_no_dimension
  Dims%number_of_atom_species                        = etsf_no_dimension
  Dims%number_of_cartesian_directions                = my_dimlen
  Dims%number_of_coefficients_dielectric_function    = etsf_no_dimension
  Dims%number_of_components                          = etsf_no_dimension
  Dims%number_of_frequencies_dielectric_function     = etsf_no_dimension
  Dims%number_of_grid_points_vector1                 = etsf_no_dimension
  Dims%number_of_grid_points_vector2                 = etsf_no_dimension
  Dims%number_of_grid_points_vector3                 = etsf_no_dimension
  Dims%number_of_kpoints                             = etsf_no_dimension
  Dims%number_of_localization_regions                = etsf_no_dimension
  Dims%number_of_qpoints_dielectric_function         = etsf_no_dimension
  Dims%number_of_qpoints_gamma_limit                 = etsf_no_dimension
  Dims%number_of_reduced_dimensions                  = my_dimlen
  Dims%number_of_spinor_components                   = etsf_no_dimension
  Dims%number_of_spins                               = etsf_no_dimension
  Dims%number_of_symmetry_operations                 = etsf_no_dimension
  Dims%number_of_vectors                             = my_dimlen
  Dims%real_or_complex_coefficients                  = etsf_no_dimension
  Dims%real_or_complex_density                       = etsf_no_dimension
  Dims%real_or_complex_gw_corrections                = etsf_no_dimension
  Dims%real_or_complex_potential                     = etsf_no_dimension
  Dims%real_or_complex_wavefunctions                 = etsf_no_dimension
  Dims%symbol_length                                 = etsf_chemlen

  !Dimensions for variables that can be splitted.
  Dims%my_max_number_of_coefficients  = etsf_no_dimension
  Dims%my_max_number_of_states        = etsf_no_dimension
  Dims%my_number_of_components        = etsf_no_dimension
  Dims%my_number_of_grid_points_vect1 = etsf_no_dimension
  Dims%my_number_of_grid_points_vect2 = etsf_no_dimension
  Dims%my_number_of_grid_points_vect3 = etsf_no_dimension
  Dims%my_number_of_kpoints           = etsf_no_dimension
  Dims%my_number_of_spins             = etsf_no_dimension

end subroutine etsf_dims_nullify
!!***
#endif

!----------------------------------------------------------------------

!!****f* m_abi_etsf/abi_etsf_dims_init
!! NAME
!! abi_etsf_dims_init
!!
!! FUNCTION
!!  Initialize the structure defining ETSF dimensions.
!!  starting from values stored in the dataset_type, the pseudopotential_type.
!!  and the wave function handler for the BIGDFT part.
!!
!! INPUTS
!!  dtset<type(dataset_type)>=all input variables for this dataset
!!  psps<type(pseudopotential_type)>=variables related to pseudopotentials
!!  wfs<wvl_wf_type>=Object to handle wave functions for the BIGDFT part
!!    Presently not used, likely will be needed when ETSF-IO will be generalized to deal with PAW.
!!  itype = an integer to define what to put in the output file. This can
!!          be one of the following values (maybe a sum latter):
!!          1 for a density file,
!!          2 for a wavefunction file,
!!          4 for a KSS file,
!!          8 for the exchange potential,
!!         16 for the correlation potential.
!!
!! OUTPUT
!!  dims=structure with ETSF dimensions.
!!
!! PARENTS
!!      m_abi_etsf
!!
!! CHILDREN
!!
!! SOURCE

#ifdef HAVE_ETSF_IO

subroutine abi_etsf_dims_init(dims, dtset, itype, psps, wfs)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itype
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_wf_type),intent(in) :: wfs
 type(etsf_dims),intent(inout) :: dims

! *************************************************************************

!Set-up the dimensions
!=====================
 dims%max_number_of_angular_momenta = psps%mpsang
 dims%max_number_of_projectors      = 1

 dims%max_number_of_states   = dtset%mband
 dims%number_of_atoms        = dtset%natom
 dims%number_of_atom_species = dtset%ntypat
 dims%number_of_components   = dtset%nspden

 dims%number_of_kpoints              = dtset%nkpt
 dims%number_of_spinor_components    = dtset%nspinor
 dims%number_of_spins                = dtset%nsppol
 dims%number_of_symmetry_operations  = dtset%nsym

!In the case of BigDFT, the number of coefficients are the number of wavelets.
 if (dtset%usewvl==0) then
   dims%max_number_of_coefficients      = dtset%mpw
   dims%max_number_of_basis_grid_points = etsf_no_dimension
 else
#ifdef HAVE_BIGDFT
   dims%max_number_of_coefficients      = wfs%ks%lzd%Glr%wfd%nvctr_c + 7 * wfs%ks%lzd%Glr%wfd%nvctr_f
   dims%max_number_of_basis_grid_points = wfs%ks%lzd%Glr%wfd%nvctr_c
#else
   BIGDFT_NOTENABLED_ERROR()
   if (.false.) write(std_out,*) wfs%ks
#endif
 end if

 if (dtset%usepaw==1) then
   dims%number_of_grid_points_vector1  = dtset%ngfftdg(1)
   dims%number_of_grid_points_vector2  = dtset%ngfftdg(2)
   dims%number_of_grid_points_vector3  = dtset%ngfftdg(3)
 else if (dtset%usewvl==1) then
#ifdef HAVE_BIGDFT
!In the case of BigDFT, the grid size is not defined by ngfft.
   dims%number_of_grid_points_vector1  = wfs%ks%lzd%Glr%d%n1 * 2
   dims%number_of_grid_points_vector2  = wfs%ks%lzd%Glr%d%n2 * 2
   dims%number_of_grid_points_vector3  = wfs%ks%lzd%Glr%d%n3 * 2
#endif
 else
   dims%number_of_grid_points_vector1  = dtset%ngfft(1)
   dims%number_of_grid_points_vector2  = dtset%ngfft(2)
   dims%number_of_grid_points_vector3  = dtset%ngfft(3)
 end if

!The density real_or_complex.
 dims%real_or_complex_density = etsf_no_dimension
 if (iand(itype, 1) /= 0 .or. dtset%usepaw==1) dims%real_or_complex_density = 1

!The coefficient of wavefunctions real_or_complex.
 dims%real_or_complex_coefficients   = etsf_no_dimension
 if (iand(itype, 2) /= 0 .or. iand(itype, 4) /= 0) then
   if (dtset%usewvl == 0) then
     dims%real_or_complex_coefficients = 2 ! used in plane waves
   else
     dims%real_or_complex_coefficients = 1 ! used in wavelets
   end if
 end if

!The gw corrections real_or_complex.
!Todo: Currently not exported.
!if (.false. .and. iand(itype, 4) /= 0) then
 dims%real_or_complex_gw_corrections = etsf_no_dimension
 if (iand(itype, 4) /= 0) then
   dims%real_or_complex_gw_corrections = 2 ! used in plane waves
 ! dims%real_or_complex_gw_corrections = 1 ! used in plane waves
 end if

!The potential real_or_complex.
 dims%real_or_complex_potential = etsf_no_dimension
 if (iand(itype, 8) /= 0 .or. iand(itype, 16) /= 0) dims%real_or_complex_potential = 1

 dims%real_or_complex_wavefunctions = etsf_no_dimension

end subroutine abi_etsf_dims_init
!!***
#endif

!----------------------------------------------------------------------

!!****f* m_abi_etsf/abi_etsf_init
!! NAME
!! abi_etsf_init
!!
!! FUNCTION
!!  Create a NetCDF file following the ETSF file format specifications.
!!  It declares the dimensions and set-up the different variables, according
!!  to the itype argument.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  filapp = character string giving the root to form the name of the GEO file
!!  itype = an integer to define what to put in the output file. This can
!!          be one of the following values (maybe a sum latter):
!!          1 for a density file,
!!          2 for a wavefunction file,
!!          4 for a KSS file,
!!          8 for the exchange potential,
!!         16 for the correlation potential.
!!  kdep= .true. if the data for the array sizes are dependant on k points.
!!  lmn_size=Number of (l,m,n) elements for the PAW basis set.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!!  Data written in file whose name is filapp//'-etsf.nc'
!!
!! PARENTS
!!      m_iowf
!!
!! CHILDREN
!!
!! SOURCE

subroutine abi_etsf_init(dtset,filapp,itype,kdep,lmn_size,psps,wfs)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itype
 logical,intent(in) :: kdep
 character(len=fnlen),intent(in) :: filapp
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_wf_type),optional,intent(in) :: wfs
!arrays
 integer,intent(in) :: lmn_size(psps%npsp)

!Local variables-------------------------------
#ifdef HAVE_ETSF_IO
!scalars
 integer :: ncid,var_main,usewvl
 logical :: ok
 character(len=80) :: file_title
 character(len=fnlen) :: filetsf
 type(etsf_dims) :: dims
 type(etsf_groups_flags) :: flags
 type(etsf_io_low_error) :: error
#endif

! *************************************************************************

#ifdef HAVE_ETSF_IO
!Initialize the filename
 filetsf = nctk_ncify(filapp)
 call wrtout(std_out,"about to create file "//TRIM(filetsf),'COLL')

 usewvl = dtset%usewvl

!Set-up the dimensions
!=====================
 call abi_etsf_dims_init(dims,dtset,itype,psps,wfs)

!Set-up the variables
!====================
!These mandatory values are always written by the hdr_io_etsf() routine.
 flags%geometry  = etsf_geometry_all
 flags%kpoints   = etsf_kpoints_red_coord_kpt + etsf_kpoints_kpoint_weights
 flags%electrons = etsf_electrons_all - etsf_electrons_x_functional - etsf_electrons_c_functional
 flags%basisdata = etsf_basisdata_basis_set

 if (usewvl==0) then
   flags%basisdata = flags%basisdata + etsf_basisdata_kin_cutoff + etsf_basisdata_n_coeff
 end if

!These variables may be written depending on prt<something> input variables.
 if (itype==1) then
   flags%main = etsf_main_density
   file_title = "Density file"
 else if (itype == 2) then
   if (usewvl==0) then
     flags%basisdata = flags%basisdata + etsf_basisdata_red_coord_pw
   else
     flags%basisdata = flags%basisdata + etsf_basisdata_coord_grid + etsf_basisdata_n_coeff_grid
   end if
   flags%main = etsf_main_wfs_coeff
   file_title = "Wavefunctions file"
 else if (itype==4) then
   if (usewvl==0) then
     flags%basisdata = flags%basisdata + etsf_basisdata_red_coord_pw
   else
     flags%basisdata = flags%basisdata + etsf_basisdata_coord_grid + etsf_basisdata_n_coeff_grid
   end if
   flags%main   = etsf_main_wfs_coeff
   flags%gwdata = etsf_gwdata_all
   file_title = "KSS file"
 else if (itype==8) then
   flags%main = etsf_main_pot_x_only
   file_title = "Exchange potential file"
 else if (itype==16) then
   flags%main = etsf_main_pot_c_only
   file_title = "Correlation potential file"
 else if (itype==24) then
   flags%main = etsf_main_pot_xc
   file_title = "Exchange-correlation potential file"
 end if

!Actually create the file
!========================
!If the group contains main, we remove it for a while to be sure to
!add it at the end, after ABINIT private variables.
 var_main   = flags%main
 flags%main = etsf_main_none

 call etsf_io_data_init(filetsf, flags, dims, file_title, &
& 'File generated by ABINIT with ETSF_IO', ok, error, overwrite = .true., k_dependent = kdep)
 ETSF_CHECK_ERROR(ok, error)

!Add the private ABINIT information when required.
 NCF_CHECK(nctk_open_modify(ncid, filetsf, xmpi_comm_self))

!Add the private data
 call ini_wf_etsf(ncid,usewvl,lmn_size,psps%npsp,psps%ntypat)

! Add the main part as last variables in the ETSF file.
 call etsf_io_main_def(ncid, ok, error, flags = var_main)
 ETSF_CHECK_ERROR(ok, error)

! Close the file.
 NCF_CHECK(nf90_close(ncid))
#else
 ABI_UNUSED(dtset%ntypat)
 ABI_UNUSED(filapp)
 ABI_UNUSED(itype)
 ABI_UNUSED(kdep)
 ABI_UNUSED(lmn_size)
 if(.false. .and. present(wfs))write(std_out,*)' OK'
#endif

end subroutine abi_etsf_init
!!***

!----------------------------------------------------------------------

!!****f* m_abi_etsf/ini_wf_etsf
!! NAME
!! ini_wf_etsf
!!
!! FUNCTION
!! Do initialization of additional dimensions and variables in wavefunction files in ETSF format.
!!
!! INPUTS
!!  usewvl=1 if wavelets are used, 0 otherwise
!!  lmn_size=Number of (l,m,n) elements for the PAW basis set.
!!  npsp=number of pseudopotentials.
!!  ntypat=number of types of atoms in cell.
!!  ncid=the unit of the open NetCDF file.
!!
!! SIDE EFFECTS
!!  New dimensions and variables are added to the initial NetCDF file.
!!
!! PARENTS
!!      m_abi_etsf
!!
!! CHILDREN
!!
!! SOURCE

subroutine ini_wf_etsf(ncid,usewvl,lmn_size,npsp,ntypat)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,usewvl,npsp,ntypat
!arrays
 integer,intent(in) :: lmn_size(npsp)

!Local variables-------------------------------
#ifdef HAVE_ETSF_IO
 integer :: ncerr
#endif

! *************************************************************************

#ifdef HAVE_ETSF_IO
!Define dimensions.
 ncerr = nctk_def_dims(ncid, [nctkdim_t("npsp", npsp), nctkdim_t("codvsnlen", 6),nctkdim_t("psptitlen", 132)])
 NCF_CHECK(ncerr)

 if (usewvl==1) then ! Add the BigDFT private dimensions.
   ncerr = nctk_def_dims(ncid, nctkdim_t("number_of_wavelet_resolutions", 2))
   NCF_CHECK(ncerr)
 end if

!Define scalars.
!If variables already exist, it will check that definitions are coherent.
 ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
&   "date", "ixc", "intxc", "occopt", "pertcase", "headform", "fform", "usepaw", "usewvl"])
 NCF_CHECK(ncerr)

 ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
&  "ecut_eff", "ecutdg", "ecutsm", "etot", "residm", "stmbias", "tphysel", "tsmear"])
 NCF_CHECK(ncerr)

!Multi-dimensional variables.
 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t("istwfk",  "i", "number_of_kpoints"),&
   nctkarr_t("codvsn", "c", "codvsnlen"),&
   nctkarr_t("pspcod",  "i", "npsp"),&
   nctkarr_t("pspdat",  "i", "npsp"),&
   nctkarr_t("pspso",  "i", "npsp"),&
   nctkarr_t("pspxc",  "i", "npsp"),&
   nctkarr_t("qptn",  "dp", "number_of_reduced_dimensions"),&
   nctkarr_t("so_psp",  "i", "npsp"),&
   nctkarr_t("symafm",  "i", "number_of_symmetry_operations"),&
   nctkarr_t("title", "c", "psptitlen, npsp"),&
   nctkarr_t("zionpsp",  "dp", "npsp"),&
   nctkarr_t("znuclpsp",  "dp", "npsp"),&
   nctkarr_t("lmn_size",  "i", "npsp")])
 NCF_CHECK(ncerr)

!Add the BigDFT private variables.
 if (usewvl == 1) then
   ncerr = nctk_def_arrays(ncid, nctkarr_t("number_of_wavelets", "i", "number_of_wavelet_resolutions"))
   NCF_CHECK(ncerr)
 end if
#else

!if ETSF_IO is undefined, do nothing
 ABI_UNUSED(ncid)
 ABI_UNUSED(ntypat)
 ABI_UNUSED(lmn_size)
 ABI_UNUSED(usewvl)
#endif

end subroutine ini_wf_etsf
!!***

!----------------------------------------------------------------------

END MODULE m_abi_etsf
!!***
