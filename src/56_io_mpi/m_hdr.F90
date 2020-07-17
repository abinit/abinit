!!****m* ABINIT/m_hdr
!! NAME
!! m_hdr
!!
!! FUNCTION
!! This module contains the definition of the abinit header and its methods
!! If you have to change the hdr, pay attention to the following subroutines:
!!
!!   hdr_malloc, hdr_init_lowlvl, hdr_free, hdr_bcast and the IO routines
!!   hdr_mpio_skip, hdr_fort_read, hdr_fort_write, hdr_ncread, hdr_ncwrite
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (XG, MB, MT, DC, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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

!#define DEBUG_MODE

! This option enable the output of the new hdr entries in hdr%echo
! Reference files should be updated
!#define DEV_NEW_HDR

module m_hdr

 use defs_basis
 use m_build_info
 use m_xmpi
 use m_abicore
 use m_errors
 use m_crystal
 use m_wffile
 use m_sort
#ifdef HAVE_MPI2
 use mpi
#endif
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk
 use m_dtset

 use m_copy,          only : alloc_copy
 use m_io_tools,      only : flush_unit, isncfile, file_exists, open_file
 use m_fstrings,      only : sjoin, itoa, ftoa, ltoa, replace_ch0, startswith, endswith, ljust, strcat, atoi
 use m_symtk,         only : print_symmetries
 !use m_kpts,          only : kpts_timrev_from_kptopt
 use defs_wvltypes,   only : wvl_internal_type
 use defs_datatypes,  only : ebands_t, pseudopotential_type
 use m_pawtab,        only : pawtab_type
 use m_pawrhoij,      only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, &
                             pawrhoij_io, pawrhoij_inquire_dim

 implicit none

 private
!!***

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!----------------------------------------------------------------------

!!****t* m_hdr/hdr_type
!! NAME
!! hdr_type
!!
!! FUNCTION
!! It contains all the information needed to write a header for a wf, den or pot file.
!! The structure of the header is explained in the abinit_help.html and other associated html files.
!! The datatype is considered as an object, to which are attached a whole
!! set of "methods", actually, different subroutines.
!! A few of these subroutines are: hdr_init, hdr_update, hdr_free, hdr_check, hdr_io, hdr_skip.
!!
!! SOURCE

 type, public :: hdr_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

  integer :: bantot        ! total number of bands (sum of nband on all kpts and spins)
  integer :: date          ! starting date
  integer :: headform      ! format of the header
  integer :: intxc         ! input variable
  integer :: ixc           ! input variable
  integer :: mband         ! maxval(hdr%nband)
  integer :: natom         ! input variable
  integer :: nkpt          ! input variable
  integer :: npsp          ! input variable
  integer :: nspden        ! input variable
  integer :: nspinor       ! input variable
  integer :: nsppol        ! input variable
  integer :: nsym          ! input variable
  integer :: ntypat        ! input variable
  integer :: occopt        ! input variable
  integer :: pertcase      ! the index of the perturbation, 0 if GS calculation
  integer :: usepaw        ! input variable (0=norm-conserving psps, 1=paw)
  integer :: usewvl        ! input variable (0=plane-waves, 1=wavelets)

  integer :: kptopt          ! input variable (defines symmetries used for k-point sampling)
  integer :: pawcpxocc       ! input variable
  integer :: nshiftk_orig=1  ! original number of shifts given in input (changed in inkpts, the actual value is nshiftk)
  integer :: nshiftk=1       ! number of shifts after inkpts.
  integer :: icoulomb        ! input variable.

  real(dp) :: ecut         ! input variable
  real(dp) :: ecutdg       ! input variable (ecut for NC psps, pawecutdg for paw)
  real(dp) :: ecutsm       ! input variable
  real(dp) :: ecut_eff     ! ecut*dilatmx**2 (dilatmx is an input variable)
  real(dp) :: etot         ! EVOLVING variable
  real(dp) :: fermie       ! EVOLVING variable
  real(dp) :: residm       ! EVOLVING variable
  real(dp) :: stmbias      ! input variable
  real(dp) :: tphysel      ! input variable
  real(dp) :: tsmear       ! input variable
  real(dp) :: nelect       ! number of electrons (computed from pseudos and charge)
  real(dp) :: charge       ! input variable

  ! This record is not a part of the hdr_type, although it is present in the
  ! header of the files. This is because it depends on the kind of file
  ! that is written, while all other information does not depend on it.
  ! It was preferred to let it be initialized or defined outside of hdr_type.
  ! integer :: fform         ! file format

  real(dp) :: qptn(3)
  ! the wavevector, in case of a perturbation

  real(dp) :: rprimd(3,3)
  ! EVOLVING variables

  integer :: ngfft(3)
  ! input variable

  integer :: nwvlarr(2)
  ! nwvlarr(2) array holding the number of wavelets for each resolution.

  integer :: kptrlatt_orig(3,3)
  ! Original kptrlatt

  integer :: kptrlatt(3,3)
  ! kptrlatt after inkpts.

  integer, allocatable :: istwfk(:)
  ! input variable istwfk(nkpt)

  integer, allocatable :: lmn_size(:)
  ! lmn_size(npsp) from psps

  integer, allocatable :: nband(:)
  ! input variable nband(nkpt*nsppol)

  integer, allocatable :: npwarr(:)
  ! npwarr(nkpt) array holding npw for each k point

  integer, allocatable :: pspcod(:)
  ! pscod(npsp) from psps

  integer, allocatable :: pspdat(:)
  ! psdat(npsp) from psps

  integer, allocatable :: pspso(:)
  ! pspso(npsp) from psps

  integer, allocatable :: pspxc(:)
  ! pspxc(npsp) from psps

  integer, allocatable :: so_psp(:)
  ! input variable so_psp(npsp)

  integer, allocatable :: symafm(:)
  ! input variable symafm(nsym)

  integer, allocatable :: symrel(:,:,:)
  ! input variable symrel(3,3,nsym)

  integer, allocatable :: typat(:)
  ! input variable typat(natom)

  real(dp), allocatable :: kptns(:,:)
  ! input variable kptns(3,nkpt)

  real(dp), allocatable :: occ(:)
  ! EVOLVING variable occ(bantot)

  real(dp), allocatable :: tnons(:,:)
  ! input variable tnons(3,nsym)

  real(dp), allocatable :: wtk(:)
  ! weight of kpoints wtk(nkpt)

  real(dp),allocatable :: shiftk_orig(:,:)
  ! original shifts given in input (changed in inkpts).

  real(dp),allocatable :: shiftk(:,:)
  ! shiftk(3,nshiftk), shiftks after inkpts

  real(dp),allocatable :: amu(:)
  ! amu(ntypat) ! EVOLVING variable

  real(dp), allocatable :: xred(:,:)
  ! EVOLVING variable xred(3,natom)

  real(dp), allocatable :: zionpsp(:)
  ! zionpsp(npsp) from psps

  real(dp), allocatable :: znuclpsp(:)
  ! znuclpsp(npsp) from psps
  ! Note the difference between (znucl|znucltypat) and znuclpsp !

  real(dp), allocatable :: znucltypat(:)
  ! znucltypat(ntypat) from alchemy

  character(len=8) :: codvsn
  ! version of the code

  character(len=132), allocatable :: title(:)
  ! title(npsp) from psps

  character(len=md5_slen),allocatable :: md5_pseudos(:)
  ! md5pseudos(npsp)
  ! md5 checksums associated to pseudos (read from file)

  ! EVOLVING variable, only for paw
  type(pawrhoij_type), allocatable :: pawrhoij(:)

  contains

  procedure :: free => hdr_free
  ! Deallocates the components of the header.

  procedure :: get_nelect_from_occ => hdr_get_nelect_from_occ
   ! Returns the number of electrons calculated from the occupation factors Hdr%occ

  procedure :: ncwrite => hdr_ncwrite
   ! Writes the header and fform to a Netcdf file.

  procedure :: vs_dtset => hdr_vs_dtset
   ! Check the compatibility of header with dtset.

  procedure :: get_crystal => hdr_get_crystal
   ! Return the crystal structure stored in the header.

  procedure :: bcast => hdr_bcast
   ! Broadcast the header.

  procedure :: compare => hdr_compare
   ! Compare two headers

  procedure :: update => hdr_update
   ! Update the header.

  procedure :: write_to_fname => hdr_write_to_fname
   ! Write the header (requires a string with the file name).

  procedure :: fort_write => hdr_fort_write
   ! Writes the header and fform to unformatted file

  procedure :: backspace => hdr_backspace
   ! Backspace the header (Fortran IO).

  procedure :: echo => hdr_echo
   ! Echo the header.

 end type hdr_type
!!***

 public :: hdr_init                ! Initialize the header and most of its content from dtset and psps.
 public :: hdr_init_lowlvl         ! Low level initialization method for Hdr (no dtset).
 public :: hdr_copy                ! Deep copy of the Header.
 public :: hdr_mpio_skip           ! Skip the abinit header using MPI-IO routines.
                                   ! Return the offset of the first Fortran record after the header.
 public :: hdr_bsize_frecords      ! Compute the size of the Fortran records from the header and formeig.
 public :: hdr_read_from_fname     ! Read the header (requires a string with the file name).
 public :: hdr_skip                ! Skip the header.
 public :: hdr_io                  ! IO of the header.
 public :: hdr_fort_read           ! Reads the header from a logical unit associated to an unformatted file.
 public :: hdr_ncread              ! Reads the header from a Netcdf file.
 public :: hdr_check               ! Compare two headers.
 public :: hdr_compare             ! Check the consistency of two headers.

 public :: abifile_from_varname
 public :: abifile_from_fform
 public :: fform_from_ext          ! Return the value of fform to be used from the file extension.
 public :: varname_from_fname      ! Return the name of the netcdf variable stored in a file from the file extension.

 ! Generic interface of the routines hdr_skip
 interface hdr_skip
   module procedure hdr_skip_int
   module procedure hdr_skip_wfftype
 end interface hdr_skip

 ! Generic interface of the routines hdr_io
 interface hdr_io
   module procedure hdr_io_int
   module procedure hdr_io_wfftype
 end interface hdr_io

 integer,private,parameter :: HDR_KNOWN_HEADFORMS(1) = [80]
 ! The list of headforms used so far.

 integer,private,parameter :: size_hdr_known_headforms = size(HDR_KNOWN_HEADFORMS) ! Need this for Flang
 integer,public,parameter :: HDR_LATEST_HEADFORM = HDR_KNOWN_HEADFORMS(size_hdr_known_headforms)
 ! The latest headform to be used for writing.
!!***

!!****t* m_hdr/abifile_t
!! NAME
!!  abifile_t
!!
!! FUNCTION
!!  Gather information about a binary file with header.
!!  Every file with header must be registered in all_abifiles, see below.
!!
!! SOURCE

 type,public :: abifile_t

   character(len=nctk_slen) :: varname
   ! Name of the netcdf variable associated to the file.
   ! This string is used in fftdatar_write to find the value of fform to be written to file

   integer :: fform
   ! The value of fform associated to this file

   character(len=24) :: ext
   ! Abinit File extension (`.nc` is not included)

   character(len=24) :: class
   ! Each file belongs to a class e.g. wf_planewave, den, pot, data...

   logical :: has_pawrhoij=.True.
   ! True if this file contains pawrhoij when hdr%usepaw == 1.

 end type abifile_t

 ! Notes about abifiles:
 !
 ! *) fform are positive integers >0 and they must be unique inside the list.
 !    One might have used strings xxx.yyy.zzz instead of integers but there's a lot of code
 !    around that relies on integral fforms so we have to live with it.
 !
 ! *) class is used in postprocessing tools e.g. cut3d when we need to know if we are dealing
 !    with wavefunctions or density-like or potential-like data.
 !    Possible values are: "wf_planewave" for wavefunction files, "density" for density-like,
 !    "potential" for potential files, "data" for generic data a.k.a. internal files e.g. GKK matrix elements.
 !
 ! *) varname can appear multiple times, in this case the entries should be ordered chronologically
 !    i.e. the most recent format should follow the older ones. This could be useful if we decide to
 !    remove pawrhoij from a particular file. Let's assume, for example, that we've decided to remove
 !    pawrhoij from the POT file. In this case, abifiles should contain:
 !
 !        abifile_t(varname="potential", fform=102), &                         ! old file with pawrhoij
 !        abifile_t(varname="potential", fform=202, has_pawrhoij=.False.) &    ! new file wo pawrhoij
 !
 ! *) The file extensions is used in fform_from_ext and varname_from_fname.
 !    fform_from_ext returns the most recent fform associated to a file extension.
 !    varname_from_fname is used in post-processing tools e.g. cut3d
 !    to read data from netcdf file without having to prompt the user for the variable name.
 !    In principle, the extension should be unique but there are exceptions e.g. the WFK produced in bigdft mode.
 !    Moreover the files produced by the DFPT code do not have a well-defined extension and, as a consequence,
 !    they require a special treatment. In python I would use regexp but Fortran is not python!

 type(abifile_t),private,parameter :: all_abifiles(50) = [ &

    ! Files with wavefunctions:
    abifile_t(varname="coefficients_of_wavefunctions", fform=2, ext="WFK", class="wf_planewave"), &
    abifile_t(varname="real_space_wavefunctions", fform=200, ext="WFK", class="wf_wavelet"), &    ! Used by wavelets.
    abifile_t(varname="ur_ae", fform=602, ext="PAWAVES", class="wf_rspace"), &                    ! Used in pawmkaewf.
    abifile_t(varname="coefficients_of_wavefunctions", fform=502, ext="KSS", class="wf_planewave"), &

    ! Files with density-like data.
    abifile_t(varname="density", fform=52, ext="DEN", class="density"), &    ! Official
    abifile_t(varname="positron_density", fform=53, ext="POSITRON", class="density"), &
    abifile_t(varname="first_order_density", fform=54, ext="DEN(\d+)", class="density"), &
    abifile_t(varname="pawrhor", fform=55, ext="PAWDEN", class="density"), &
    abifile_t(varname="pawrhor_core", fform=56, ext="ATMDEN_CORE", class="density"), &
    abifile_t(varname="pawrhor_val", fform=57, ext="ATMDEN_VAL", class="density"), &
    abifile_t(varname="pawrhor_full", fform=58, ext="ATMDEN_FULL", class="density"), &
    abifile_t(varname="pawrhor_ntilde_minus_nhat", fform=59, ext="N_TILDE", class="density"), &
    abifile_t(varname="pawrhor_n_one", fform=60, ext="N_ONE", class="density"), &
    abifile_t(varname="pawrhor_nt_one", fform=61, ext="NT_ONE", class="density"), &
    abifile_t(varname="qp_rhor", fform=62, ext="QP_DEN", class="density"), &
    abifile_t(varname="qp_pawrhor", fform=63, ext="QP_PAWDEN", class="density"), &
    abifile_t(varname="grhor_1", fform=67, ext="GDEN1", class="density"), &
    abifile_t(varname="grhor_2", fform=68, ext="GDEN2", class="density"), &
    abifile_t(varname="grhor_3", fform=69, ext="GDEN3", class="density"), &

    !???
    abifile_t(varname="stm", fform=110, ext="STM", class="density"), &
    abifile_t(varname="kinedr", fform=70, ext="KDEN", class="density"), &
    abifile_t(varname="elfr", fform=64, ext="ELF", class="density"), &
    abifile_t(varname="elfr_up", fform=65, ext="ELF_UP", class="density"), &
    abifile_t(varname="elfr_down", fform=66, ext="ELF_DOWN", class="density"), &
    abifile_t(varname="laprhor", fform=71, ext="LDEN", class="density"), &

    ! Files with potentials
    ! Official
    abifile_t(varname="potential", fform=102, ext="POT", class="potential"), &  ! CHECK THESE TWO FILES
    abifile_t(varname="vtrial", fform=103, ext="POT", class="potential"), &
    abifile_t(varname="vhartree", fform=104, ext="VHA", class="potential"), &
    abifile_t(varname="vpsp", fform=105, ext="VPSP", class="potential"), &
    abifile_t(varname="vhartree_vloc", fform=106, ext="VCLMB", class="potential"), &
    abifile_t(varname="vhxc", fform=107, ext="VHXC", class="potential"), &
    abifile_t(varname="exchange_correlation_potential", fform=108, ext="VXC", class="potential"), &

    abifile_t(varname="first_order_potential", fform=109, ext="POT(\d+)", class="potential"), &
    ! fform 111 contains an extra record with rhog1_q(G=0) after the DFPT potential(r).
    abifile_t(varname="first_order_potential", fform=111, ext="POT(\d+)", class="potential"), &

    abifile_t(varname="first_order_vhartree", fform=112, ext="VHA(\d+)", class="potential"), &
    abifile_t(varname="first_order_vpsp", fform=113, ext="VPSP(\d+)", class="potential"), &
    abifile_t(varname="first_order_vxc", fform=114, ext="VXC(\d+)", class="potential"), &

   ! Data used in conducti
    abifile_t(varname="pawnabla", fform=610, ext="OPT1", class="data"), &
    abifile_t(varname="pawnabla_core", fform=611, ext="OPT2", class="data"), &
    abifile_t(varname="pawnabla_loc", fform=612, ext="OPT", class="data"), &

   ! Data used in elphon
    abifile_t(varname="gkk_elements", fform=42, ext="GKK", class="data"), &

   ! DKK matrix elements in netcdf format (optic, eph)
    abifile_t(varname="h1_matrix_elements", fform=43, ext="DKK", class="data"), &

    ! output files that are not supposed to be read by abinit.
    abifile_t(varname="this_file_is_not_read_by_abinit", fform=666, ext="666", class="data"), &

   ! GW files: old 1002, 1102
   !character(len=nctk_slen),public,parameter :: e_ncname="dielectric_function"
   ! FIXME This one should be rewritten
   abifile_t(varname="polarizability", fform=1003, ext="SUS", class="polariz"),  &
   abifile_t(varname="inverse_dielectric_function", fform=1004, ext="SCR", class="epsm1"), &
   !abifile_t(varname="dielectric_function", fform=1002, ext="EPS", class="eps"), &
   !
   ! BSE: TODO. see m_bse_io
   !abifile_t(varname="bse_uresonant_q0", fform=1002, ext="BSR", class="bsreso"), &
   !abifile_t(varname="bse_ucoupling_q0", fform=1002, ext="BSC", class="bscoup"), &

   ! Miscellaneous
   abifile_t(varname="dos_fractions", fform=3000, ext="FATBANDS", class="data"), &
   abifile_t(varname="spectral_weights", fform=5000, ext="FOLD2BLOCH", class="data"), &
   abifile_t(varname="no_fftdatar_write", fform=6000, ext="ABIWAN", class="data"), &
   abifile_t(varname="None", fform=6001, ext="KERANGE", class="data"), &
   abifile_t(varname="None", fform=6002, ext="SIGEPH", class="data") &
  ]

 type(abifile_t),public,parameter :: abifile_none = abifile_t(varname="None", fform=0, ext="None", class="None")
 ! This object is returned when we cannot find the file in abifiles.

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/fform_from_ext
!! NAME
!!  fform_from_ext
!!
!! FUNCTION
!!  Return the value of fform to be used from the file extension. If a file has multiple fforms,
!!  the most recent one is returned. Returns 0 if the extension is not registered.
!!
!! SOURCE

integer function fform_from_ext(abiext) result(fform)

!Arguments ---------------------------------------------
 character(len=*),intent(in) :: abiext

!Local variables-------------------------------
!scalars
 integer :: ii,ind,ierr,pertcase
 character(len=len(abiext)) :: ext

! *********************************************************************
 ! Remove .nc (if any) and work with ext
 ext = abiext
 if (endswith(abiext, ".nc")) then
   ind = index(abiext, ".nc", back=.True.); ext = abiext(:ind-1)
 end if

 fform = 0
 do ii=1,size(all_abifiles)
   if (ext == all_abifiles(ii)%ext) fform = all_abifiles(ii)%fform
 end do
 if (fform /= 0) return
 ! Here we handle special cases.

 ! Handle DEN[pertcase]
 if (startswith(ext, "DEN")) then
   read(ext(4:), *, iostat=ierr) pertcase
   if (ierr == 0) then
     do ii=1,size(all_abifiles)
       if (all_abifiles(ii)%ext == "DEN(\d+)") fform = all_abifiles(ii)%fform
     end do
     return
   end if
 end if

 ! Handle POT[pertcase]
 if (startswith(ext, "POT")) then
   read(ext(4:), *, iostat=ierr) pertcase
   if (ierr == 0) then
     do ii=1,size(all_abifiles)
       if (all_abifiles(ii)%ext == "POT(\d+)") fform = all_abifiles(ii)%fform
     end do
     return
   end if
 end if

 MSG_ERROR(sjoin("Cannot find fform associated to extension:", abiext))

end function fform_from_ext
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/varname_from_fname
!! NAME
!!  varname_from_fname
!!
!! FUNCTION
!!   Return the name of the netcdf variable stored in a file from the file extension.
!!
!! NOTES
!!  The variable names should be consistent with the ones used in outscf.F90
!!
!! SOURCE

character(len=nctk_slen) function varname_from_fname(filename) result(varname)

!Arguments ---------------------------------------------
 character(len=*),intent(in) :: filename

!Local variables-------------------------------
!scalars
 integer :: ind,pertcase,ierr
 logical :: found
 character(len=len(filename)) :: ext

! *********************************************************************

 ! TODO: This should be a recursive function because we have
 ! to scan the string from left to right (extensions could have the same termination)

 ! Find the Abinit file extension. Examples: t43_VHXC.nc
 if (endswith(filename, ".nc")) then
   ind = index(filename, ".nc", back=.True.)
 else
   !MSG_ERROR(sjoin("Don't know how to handle: ", filename))
   ind = len_trim(filename) + 1
 end if

 ext = filename(:ind-1)
 ind = index(ext, "_", back=.True.)
 ABI_CHECK(ind /= 0, "Cannot find `_` in file name!")
 ABI_CHECK(ind /= len_trim(ext), sjoin("Wrong string: ", ext))
 ext = ext(ind+1:)

 found = .True.
 select case (ext)
 case ("DEN")
   varname = "density"
 !case ("DEN1")
 !  varname = "first_order_density"
 case ("POSITRON")
   varname = "positron_density"
 case ("PAWDEN")
   varname = "pawrhor"
   ! TODO: Other paw densities
 case ("ELF")
   varname = "elfr"
 case ("ELF_UP")
   varname = "elfr_up"
 case ("ELF_DOWN")
   varname = "elfr_down"
 case ("GDEN1")
   varname = "grhor_1"
 case ("GDEN2")
   varname = "grhor_2"
 case ("GDEN3")
   varname = "grhor_3"
 case ("KDEN")
   varname = "kinedr"
 case ("LDEN")
   varname = "laprhor"
 case ("POT")
   varname = "vtrial"
 case ("STM")
   varname = "stm"
 case ("VHA")
   varname = "vhartree"
 case ("VPSP")
   varname = "vpsp"
 case ("VHXC")
   varname = "vhxc"
 case ("VXC")
   varname = "exchange_correlation_potential"
 case ("VCLMB")
   varname = "vhartree_vloc"
 case default
   found = .False.
 end select

 if (found) return

 ! Handle DEN[pertcase]
 if (startswith(ext, "DEN")) then
   read(ext(4:), *, iostat=ierr) pertcase
   if (ierr == 0) then
     varname = "first_order_density"; return
   end if
 end if

 ! Handle POT[pertcase]
 if (startswith(ext, "POT")) then
   read(ext(4:), *, iostat=ierr) pertcase
   if (ierr == 0) then
      varname = "first_order_potential"; return
   end if
 end if

 ! Handle VXC[pertcase]
 if (startswith(ext, "VXC")) then
   read(ext(4:), *, iostat=ierr) pertcase
   if (ierr == 0) then
      varname = "first_order_vxc"; return
   end if
 end if

 ! Handle VHA[pertcase]
 if (startswith(ext, "VHA")) then
   read(ext(4:), *, iostat=ierr) pertcase
   if (ierr == 0) then
      varname = "first_order_vhartree"; return
   end if
 end if

 ! Handle VPSP[pertcase]
 if (startswith(ext, "VPSP")) then
   read(ext(4:), *, iostat=ierr) pertcase
   if (ierr == 0) then
      varname = "first_order_vpsp"; return
   end if
 end if

 MSG_ERROR(sjoin("Unknown abinit extension:", ext))

end function varname_from_fname
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/abifile_from_varname
!! NAME
!!  abifile_from_varname
!!
!! FUNCTION
!!  Return the abifile_t object corresponding to the given variable varname
!!  Return abifile_none if not found. This function is used to find the last
!!  value of fform when we write data to file.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(abifile_t) function abifile_from_varname(varname) result(afile)

!Arguments ---------------------------------------------
 character(len=*),intent(in) :: varname

!Local variables-------------------------------
!scalars
 integer :: ii
! *************************************************************************

 afile = abifile_none
 do ii=1,size(all_abifiles)
   if (all_abifiles(ii)%varname == varname) afile = all_abifiles(ii)
 end do

end function abifile_from_varname
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/abifile_from_fform
!! NAME
!!  abifile_from_fform
!!
!! FUNCTION
!!  Return the abifile_t object corresponding to the given fform
!!  Return abifile_none if not found. This function is used to
!!  find the name of the netcdf variable from the fform and
!!  detect whether the file contains pawrhoij.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(abifile_t) function abifile_from_fform(fform) result(afile)

!Arguments ---------------------------------------------
 integer,intent(in) :: fform

!Local variables-------------------------------
!scalars
 integer :: ii
! *************************************************************************

 afile = abifile_none
 do ii=1,size(all_abifiles)
   if (all_abifiles(ii)%fform == fform) afile = all_abifiles(ii)
 end do

end function abifile_from_fform
!!***

!!****f* m_hdr/check_fform
!! NAME
!!  check_fform
!!
!! FUNCTION
!!   This function is used ifdef DEBUG_MODE. It tests whether the value of fform
!!   is registered in all_abifiles.
!!
!! PARENTS
!!      m_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine check_fform(fform)

!Local variables-------------------------------
!scalars
 integer,intent(in) :: fform
#ifdef DEBUG_MODE
 type(abifile_t) :: abifile
 character(len=500) :: msg

! *********************************************************************
 if (fform == 666) return
 abifile = abifile_from_fform(fform)

 if (abifile%fform == 0) then
    MSG_ERROR(sjoin("Cannot find any abifile object associated to fform:", itoa(fform)))
 end if
 if (abifile%fform /= fform) then
    write(msg,"(2a,2(a,i0))") &
      "Input fform does not agree with the one registered in abifile.",ch10,&
      "hdr%fform= ",fform,", abifile%fform= ",abifile%fform
    MSG_ERROR(msg)
 end if

#else
 ABI_UNUSED(fform)
#endif

end subroutine check_fform
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/test_abifiles
!! NAME
!!  test_abifiles
!!
!! FUNCTION
!!  Check the consistency of the internal abifiles table.
!!
!! PARENTS
!!      m_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine test_abifiles()

!Arguments ---------------------------------------------

!Local variables-------------------------------
!scalars
 integer :: ii,nn,ierr
 integer :: all_fforms(size(all_abifiles)),iperm(size(all_abifiles))
! *************************************************************************

 nn = size(all_abifiles)

 do ii=1,nn
   all_fforms(ii) = all_abifiles(ii)%fform
 end do
 iperm = [(ii, ii=1,nn)]
 call sort_int(nn, all_fforms, iperm)

 ierr = 0
 do ii=1,nn-1
   if (all_fforms(ii) == all_fforms(ii+1)) then
     MSG_WARNING(sjoin("fform: ", itoa(all_fforms(ii+1)), "is already in the abifiles list"))
     ierr = ierr + 1
   end if
 end do

 if (ierr /= 0) then
   MSG_ERROR("test_abifiles gave ierr != 0. Aborting now")
 end if

end subroutine test_abifiles
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_malloc
!! NAME
!! hdr_malloc
!!
!! FUNCTION
!!  Allocate memory from dimensions with the exception of pawrhoij.
!!  This is a private routine. Client code should use hdr_init, hdr_fort_read.
!!
!! PARENTS
!!      m_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_malloc(hdr, bantot, nkpt, nsppol, npsp, natom, ntypat, nsym, nshiftk_orig, nshiftk)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: bantot,nkpt,nsppol,npsp,natom,ntypat,nsym,nshiftk,nshiftk_orig
 class(hdr_type),intent(inout) :: hdr
! *************************************************************************

 !@hdt_type
 call hdr%free()

 ABI_MALLOC(hdr%istwfk, (nkpt))
 ABI_MALLOC(hdr%nband, (nkpt*nsppol))
 ABI_MALLOC(hdr%npwarr, (nkpt))
 ABI_MALLOC(hdr%pspcod, (npsp))
 ABI_MALLOC(hdr%pspdat, (npsp))
 ABI_MALLOC(hdr%pspso, (npsp))
 ABI_MALLOC(hdr%pspxc, (npsp))
 ABI_MALLOC(hdr%lmn_size, (npsp))
 ABI_MALLOC(hdr%so_psp, (npsp))
 ABI_MALLOC(hdr%symafm, (nsym))
 ABI_MALLOC(hdr%symrel, (3,3,nsym))
 ABI_MALLOC(hdr%typat, (natom))
 ABI_MALLOC(hdr%kptns, (3,nkpt))
 ABI_MALLOC(hdr%occ, (bantot))
 ABI_MALLOC(hdr%tnons, (3,nsym))
 ABI_MALLOC(hdr%wtk, (nkpt))
 ABI_MALLOC(hdr%xred, (3,natom))
 ABI_MALLOC(hdr%zionpsp, (npsp))
 ABI_MALLOC(hdr%znuclpsp, (npsp))
 ABI_MALLOC(hdr%znucltypat, (ntypat))
 ABI_MALLOC(hdr%title, (npsp))
 ABI_MALLOC(hdr%shiftk, (3,nshiftk))
 ABI_MALLOC(hdr%shiftk_orig, (3,nshiftk_orig))
 ABI_MALLOC(hdr%md5_pseudos, (npsp))
 ABI_MALLOC(hdr%amu, (ntypat))

end subroutine hdr_malloc
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_init
!! NAME
!! hdr_init
!!
!! FUNCTION
!! This subroutine initializes the header structured datatype
!! and most of its content from dtset and psps, and put default values for evolving variables.
!!
!! INPUTS
!! ebands <type(ebands_t)>=band structure information including Brillouin zone description
!! codvsn=code version
!! dtset <type(dataset_type)>=all input variables for this dataset
!! mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!! comm_atom=--optional-- MPI communicator over atoms
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! pertcase=index of the perturbation, or 0 if GS calculation
!! psps <type(pseudopotential_type)>=all the information about psps
!! my_atomtab(:)=Index of the atoms (in global numbering ) treated by current proc (Optional)
!!
!! OUTPUT
!! hdr <type(hdr_type)>=the header, initialized, and for most part of
!!   it, contain its definite values, except for evolving variables
!!
!! PARENTS
!!      dfpt_looppert,gstate,nonlinear,respfn,setup_bse,setup_screening
!!      setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_init(ebands,codvsn,dtset,hdr,pawtab,pertcase,psps,wvl, &
                    mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: pertcase
 integer,intent(in),optional :: comm_atom
 character(len=8),intent(in) :: codvsn
 type(ebands_t),intent(in) :: ebands
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr !vz_i
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type),intent(in) :: wvl
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: image=1
 character(len=500) :: msg

! *************************************************************************

#ifdef DEBUG_MODE
 call test_abifiles()
#endif

 !@hdr_type

! More checking would be needed ...
 if (dtset%ntypat/=psps%ntypat) then
   write(msg,'(a,2(i0,2x))')' dtset%ntypat and psps%ntypat differs. They are: ',dtset%ntypat,psps%ntypat
   MSG_ERROR(msg)
 end if

 if (dtset%npsp/=psps%npsp) then
   write(msg,'(a,2(i0,2x))')' dtset%npsp and psps%npsp differs. They are: ',dtset%npsp,psps%npsp
   MSG_ERROR(msg)
 end if

 ! Note: The structure parameters are taken from the first image!
 if (present(comm_atom)) then
   if (present(mpi_atmtab)) then
     call hdr_init_lowlvl(hdr,ebands,psps,pawtab,wvl,codvsn,pertcase,&
       dtset%natom,dtset%nsym,dtset%nspden,dtset%ecut,dtset%pawecutdg,dtset%ecutsm,dtset%dilatmx,&
       dtset%intxc,dtset%ixc,dtset%stmbias,dtset%usewvl,dtset%pawcpxocc,dtset%pawspnorb,dtset%ngfft,dtset%ngfftdg,&
       dtset%so_psp,dtset%qptn, dtset%rprimd_orig(:,:,image),dtset%xred_orig(:,:,image),&
       dtset%symrel,dtset%tnons,dtset%symafm,dtset%typat,dtset%amu_orig(:,image),dtset%icoulomb,&
       dtset%kptopt,dtset%nelect,dtset%charge,dtset%kptrlatt_orig,dtset%kptrlatt,&
       dtset%nshiftk_orig,dtset%nshiftk,dtset%shiftk_orig,dtset%shiftk,&
       comm_atom=comm_atom,mpi_atmtab=mpi_atmtab)
   else
     call hdr_init_lowlvl(hdr,ebands,psps,pawtab,wvl,codvsn,pertcase,&
       dtset%natom,dtset%nsym,dtset%nspden,dtset%ecut,dtset%pawecutdg,dtset%ecutsm,dtset%dilatmx,&
       dtset%intxc,dtset%ixc,dtset%stmbias,dtset%usewvl,dtset%pawcpxocc,dtset%pawspnorb,dtset%ngfft,dtset%ngfftdg,&
       dtset%so_psp,dtset%qptn, dtset%rprimd_orig(:,:,image),dtset%xred_orig(:,:,image),&
       dtset%symrel,dtset%tnons,dtset%symafm,dtset%typat,dtset%amu_orig(:,image),dtset%icoulomb,&
       dtset%kptopt,dtset%nelect,dtset%charge,dtset%kptrlatt_orig,dtset%kptrlatt,&
       dtset%nshiftk_orig,dtset%nshiftk,dtset%shiftk_orig,dtset%shiftk,&
       comm_atom=comm_atom)
   end if
 else
   call hdr_init_lowlvl(hdr,ebands,psps,pawtab,wvl,codvsn,pertcase,&
     dtset%natom,dtset%nsym,dtset%nspden,dtset%ecut,dtset%pawecutdg,dtset%ecutsm,dtset%dilatmx,&
     dtset%intxc,dtset%ixc,dtset%stmbias,dtset%usewvl,dtset%pawcpxocc,dtset%pawspnorb,dtset%ngfft,dtset%ngfftdg,&
     dtset%so_psp,dtset%qptn, dtset%rprimd_orig(:,:,image),dtset%xred_orig(:,:,image),dtset%symrel,&
     dtset%tnons,dtset%symafm,dtset%typat,dtset%amu_orig(:,image),dtset%icoulomb,&
     dtset%kptopt,dtset%nelect,dtset%charge,dtset%kptrlatt_orig,dtset%kptrlatt,&
     dtset%nshiftk_orig,dtset%nshiftk,dtset%shiftk_orig,dtset%shiftk)
 end if

end subroutine hdr_init
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_free
!! NAME
!! hdr_free
!!
!! FUNCTION
!! This subroutine deallocates the components of the header structured datatype
!!
!! INPUTS
!! hdr <type(hdr_type)>=the header
!!
!! OUTPUT
!!  (only deallocate)
!!
!! PARENTS
!!      bethe_salpeter,conducti_nc,conducti_paw,conducti_paw_core,cut3d
!!      dfpt_looppert,dfptnl_loop,elphon,emispec_paw,eph,finddistrproc,gstate
!!      initaim,inpgkk,inwffil,ioprof,linear_optics_paw,m_bse_io,m_cut3d,m_ddk
!!      m_dvdb,m_hdr,m_io_kss,m_io_screening,m_ioarr,m_wfd,m_wfk,macroave
!!      mrggkk,mrgscr,nonlinear,optic,read_el_veloc,read_gkk,respfn,screening
!!      sigma,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_free(hdr)

!Arguments ------------------------------------
!scalars
 class(hdr_type),intent(inout) :: hdr

! *************************************************************************

 DBG_ENTER("COLL")

 !@hdr_type

 !integer
 ABI_SFREE(hdr%istwfk)
 ABI_SFREE(hdr%lmn_size)
 ABI_SFREE(hdr%nband)
 ABI_SFREE(hdr%npwarr)
 ABI_SFREE(hdr%pspcod)
 ABI_SFREE(hdr%pspdat)
 ABI_SFREE(hdr%pspso)
 ABI_SFREE(hdr%pspxc)
 ABI_SFREE(hdr%so_psp)
 ABI_SFREE(hdr%symafm)
 ABI_SFREE(hdr%symrel)
 ABI_SFREE(hdr%typat)

 !real
 ABI_SFREE(hdr%amu)
 ABI_SFREE(hdr%kptns)
 ABI_SFREE(hdr%occ)
 ABI_SFREE(hdr%tnons)
 ABI_SFREE(hdr%wtk)
 ABI_SFREE(hdr%shiftk)
 ABI_SFREE(hdr%shiftk_orig)
 ABI_SFREE(hdr%xred)
 ABI_SFREE(hdr%zionpsp)
 ABI_SFREE(hdr%znuclpsp)
 ABI_SFREE(hdr%znucltypat)

 !string arrays
 ABI_SFREE(hdr%md5_pseudos)
 ABI_SFREE(hdr%title)

 if (hdr%usepaw==1 .and. allocated(hdr%pawrhoij) ) then
   call pawrhoij_free(hdr%pawrhoij)
   ABI_FREE(hdr%pawrhoij)
 end if

 DBG_EXIT("COLL")

end subroutine hdr_free
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_copy
!! NAME
!! hdr_copy
!!
!! FUNCTION
!! Deep copy of the abinit header.
!!
!! INPUTS
!!  Hdr_in=The header to be copied.
!!
!! OUTPUT
!!  Hdr_cp=The deep copy of Hdr_in.
!!
!! NOTES
!!  The present version deals with versions of the header up to 56.
!!
!! PARENTS
!!      dfpt_looppert,m_io_kss,m_io_screening,m_wfk,optic
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_copy(Hdr_in,Hdr_cp)

!Arguments ------------------------------------
!scalars
 class(hdr_type),intent(in) :: Hdr_in
 type(hdr_type),intent(inout) :: Hdr_cp

!Local variables-------------------------------
!scalars
 integer :: cplex_rhoij,nspden_rhoij,qphase_rhoij

! *************************************************************************

 !@hdr_type

! Integer values
 Hdr_cp%bantot   = Hdr_in%bantot
 Hdr_cp%date     = Hdr_in%date
 Hdr_cp%headform = Hdr_in%headform
 hdr_cp%icoulomb = hdr_in%icoulomb
 Hdr_cp%intxc    = Hdr_in%intxc
 Hdr_cp%ixc      = Hdr_in%ixc
 Hdr_cp%natom    = Hdr_in%natom
 Hdr_cp%nkpt     = Hdr_in%nkpt
 Hdr_cp%npsp     = Hdr_in%npsp
 Hdr_cp%nspden   = Hdr_in%nspden
 Hdr_cp%nspinor  = Hdr_in%nspinor
 Hdr_cp%nsppol   = Hdr_in%nsppol
 Hdr_cp%nsym     = Hdr_in%nsym
 Hdr_cp%ntypat   = Hdr_in%ntypat
 Hdr_cp%occopt   = Hdr_in%occopt
 Hdr_cp%pertcase = Hdr_in%pertcase
 Hdr_cp%usepaw   = Hdr_in%usepaw
 Hdr_cp%usewvl   = Hdr_in%usewvl
 Hdr_cp%mband    = Hdr_in%mband
 ABI_CHECK(hdr_in%mband == maxval(hdr_in%nband), "mband != maxval(hdr_in%nband)")
 hdr_cp%kptopt = hdr_in%kptopt
 hdr_cp%pawcpxocc = hdr_in%pawcpxocc
 hdr_cp%nshiftk_orig = hdr_in%nshiftk_orig
 hdr_cp%nshiftk = hdr_in%nshiftk

 ! Integer arrays
 Hdr_cp%ngfft   = Hdr_in%ngfft
 Hdr_cp%nwvlarr = Hdr_in%nwvlarr
 hdr_cp%kptrlatt = hdr_in%kptrlatt
 hdr_cp%kptrlatt_orig = hdr_in%kptrlatt_orig

! Integer allocatable arrays
 call alloc_copy( Hdr_in%istwfk,  Hdr_cp%istwfk   )
 call alloc_copy( Hdr_in%lmn_size,Hdr_cp%lmn_size )
 call alloc_copy( Hdr_in%nband,   Hdr_cp%nband    )
 call alloc_copy( Hdr_in%npwarr,  Hdr_cp%npwarr   )
 call alloc_copy( Hdr_in%pspcod,  Hdr_cp%pspcod )
 call alloc_copy( Hdr_in%pspdat,  Hdr_cp%pspdat )
 call alloc_copy( Hdr_in%pspso ,  Hdr_cp%pspso  )
 call alloc_copy( Hdr_in%pspxc ,  Hdr_cp%pspxc  )
 call alloc_copy( Hdr_in%so_psp,  Hdr_cp%so_psp )
 call alloc_copy( Hdr_in%symafm,  Hdr_cp%symafm )
 call alloc_copy( Hdr_in%symrel,  Hdr_cp%symrel )
 call alloc_copy( Hdr_in%typat ,  Hdr_cp%typat  )

! Real variables
 Hdr_cp%ecut        = Hdr_in%ecut
 Hdr_cp%ecutdg      = Hdr_in%ecutdg
 Hdr_cp%ecutsm      = Hdr_in%ecutsm
 Hdr_cp%ecut_eff    = Hdr_in%ecut_eff
 Hdr_cp%etot        = Hdr_in%etot
 Hdr_cp%fermie      = Hdr_in%fermie
 Hdr_cp%residm      = Hdr_in%residm
 Hdr_cp%stmbias     = Hdr_in%stmbias
 Hdr_cp%tphysel     = Hdr_in%tphysel
 Hdr_cp%tsmear      = Hdr_in%tsmear
 hdr_cp%nelect      = hdr_in%nelect
 hdr_cp%charge      = hdr_in%charge

 Hdr_cp%qptn(:)     = Hdr_in%qptn(:)
 Hdr_cp%rprimd(:,:) = Hdr_in%rprimd(:,:)

! Real allocatable arrays
 call alloc_copy(Hdr_in%amu, Hdr_cp%amu)
 call alloc_copy( Hdr_in%kptns     ,Hdr_cp%kptns     )
 call alloc_copy( Hdr_in%occ       ,Hdr_cp%occ       )
 call alloc_copy( Hdr_in%tnons     ,Hdr_cp%tnons     )
 call alloc_copy( Hdr_in%wtk       ,Hdr_cp%wtk       )
 call alloc_copy( Hdr_in%xred      ,Hdr_cp%xred      )
 call alloc_copy( Hdr_in%zionpsp   ,Hdr_cp%zionpsp   )
 call alloc_copy( Hdr_in%znuclpsp  ,Hdr_cp%znuclpsp  )
 call alloc_copy( Hdr_in%znucltypat,Hdr_cp%znucltypat)
 call alloc_copy(Hdr_in%shiftk, Hdr_cp%shiftk)
 call alloc_copy(Hdr_in%shiftk_orig, Hdr_cp%shiftk_orig)

! Character arrays
 Hdr_cp%codvsn = Hdr_in%codvsn
! THIS DOES NOT WORK ON XLF: Hdr_cp%title string length becomes huge and segfaults
! call alloc_copy( Hdr_in%title,Hdr_cp%title )
 ABI_MALLOC(Hdr_cp%title,(Hdr_cp%npsp))
 Hdr_cp%title = Hdr_in%title

 ABI_MALLOC(hdr_cp%md5_pseudos, (hdr_cp%npsp))
 hdr_cp%md5_pseudos = hdr_in%md5_pseudos

! For PAW have to copy Pawrhoij ====
! NOTE alchemy requires a different treatment but for the moment it is not available within PAW.
 if (Hdr_in%usepaw==1) then
   cplex_rhoij  = Hdr_in%Pawrhoij(1)%cplex_rhoij
   qphase_rhoij = Hdr_in%Pawrhoij(1)%qphase
   nspden_rhoij = Hdr_in%Pawrhoij(1)%nspden
   ABI_MALLOC(Hdr_cp%Pawrhoij,(Hdr_in%natom))
   call pawrhoij_alloc(Hdr_cp%Pawrhoij,cplex_rhoij,nspden_rhoij,Hdr_in%nspinor,Hdr_in%nsppol,Hdr_in%typat,&
                       lmnsize=Hdr_in%lmn_size(1:Hdr_in%ntypat),qphase=qphase_rhoij)
   call pawrhoij_copy(Hdr_in%Pawrhoij,Hdr_cp%Pawrhoij)
 end if

end subroutine hdr_copy
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_get_nelect_from_occ
!! NAME
!! hdr_get_nelect_from_occ
!!
!! FUNCTION
!!  Return the number of electrons from the occupation numbers
!!  This function is mainly used for debugging purposes, use hdr%nelect and hdr%charge
!!
!! INPUTS
!!  Hdr<hdr_type>
!!
!! OUTPUT
!!  nelect=Number of electrons in the unit cell.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

real(dp) pure function hdr_get_nelect_from_occ(Hdr) result(nelect)

!Arguments ---------------------------------------------
!scalars
 class(hdr_type),intent(in) :: Hdr

!Local variables ---------------------------------------
!scalars
 integer :: idx,isppol,ikibz,nband_k
! *************************************************************************

 ! Cannot use znucl because we might have additional charge or alchemy.
 nelect=zero ; idx=0
 do isppol=1,Hdr%nsppol
   do ikibz=1,Hdr%nkpt
     nband_k=Hdr%nband(ikibz+(isppol-1)*Hdr%nkpt)
     nelect = nelect + Hdr%wtk(ikibz)*SUM(Hdr%occ(idx+1:idx+nband_k))
     idx=idx+nband_k
   end do
 end do

end function hdr_get_nelect_from_occ
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_init_lowlvl
!! NAME
!! hdr_init_lowlvl
!!
!! FUNCTION
!! This subroutine initializes the header structured datatype
!! and most of its content from psps and other input variables that
!! are passed explicitly. It also use default values for evolving variables.
!! Note that Dtset is not required thus rendering the initialization of the header
!! much easier.
!!
!! INPUTS
!! ebands <type(ebands_t)>=band structure information including Brillouin zone description
!! codvsn=code version
!! mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!! comm_atom=--optional-- MPI communicator over atoms
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! pertcase=index of the perturbation, or 0 if GS calculation
!! psps <type(pseudopotential_type)>=all the information about psps
!! For the meaning of the other varialble see the definition of dataset_type.
!!
!! OUTPUT
!! hdr <type(hdr_type)>=the header, initialized, and for most part of
!!   it, contain its definite values, except for evolving variables
!!
!! PARENTS
!!      m_hdr,m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_init_lowlvl(hdr,ebands,psps,pawtab,wvl,&
  codvsn,pertcase,natom,nsym,nspden,ecut,pawecutdg,ecutsm,dilatmx,&
  intxc,ixc,stmbias,usewvl,pawcpxocc,pawspnorb,ngfft,ngfftdg,so_psp,qptn,&
  rprimd,xred,symrel,tnons,symafm,typat,amu,icoulomb,&
  kptopt,nelect,charge,kptrlatt_orig,kptrlatt,&
  nshiftk_orig,nshiftk,shiftk_orig,shiftk,&
  mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym,nspden,intxc,ixc,usewvl,pawcpxocc,pawspnorb,pertcase
 integer,intent(in) :: kptopt,nshiftk_orig,nshiftk,icoulomb
 integer, intent(in),optional :: comm_atom
 real(dp),intent(in) :: ecut,ecutsm,dilatmx,stmbias,pawecutdg,nelect,charge
 character(len=8),intent(in) :: codvsn
 type(ebands_t),intent(in) :: ebands
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type),intent(in) :: wvl
 type(hdr_type),intent(inout) :: hdr
!arrays
 integer,intent(in) :: typat(natom)
 integer,intent(in) :: so_psp(psps%npsp)
 integer,intent(in) :: symrel(3,3,nsym),symafm(nsym)
 integer,intent(in) :: ngfft(18),ngfftdg(18),kptrlatt_orig(3,3),kptrlatt(3,3)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: tnons(3,nsym),amu(psps%ntypat)
 real(dp),intent(in) :: qptn(3) ! the wavevector, in case of a perturbation
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 real(dp),intent(in) :: shiftk_orig(3,nshiftk_orig),shiftk(3,nshiftk)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: bantot,date,nkpt,npsp,ntypat,nsppol,nspinor
 integer :: cplex_rhoij,nspden_rhoij,qphase_rhoij
 integer :: idx,isppol,ikpt,iband,ipsp
 character(len=8) :: date_time

! *************************************************************************

 !@hdr_type
 call date_and_time(date_time)
 read(date_time,'(i8)')date

 npsp   = psps%npsp
 ntypat = psps%ntypat
 nkpt   = ebands%nkpt
 nsppol = ebands%nsppol
 nspinor= ebands%nspinor
 bantot = ebands%bantot

!Transfer dimensions and other scalars to hdr.
 hdr%intxc    =intxc
 hdr%ixc      =ixc
 hdr%natom    =natom
 hdr%npsp     =npsp
 hdr%nspden   =nspden
 hdr%nspinor  =nspinor
 hdr%nsym     =nsym
 hdr%ntypat   =ntypat
 hdr%bantot   =bantot
 hdr%nkpt     =nkpt
 hdr%nshiftk_orig = nshiftk_orig
 hdr%nshiftk = nshiftk
 hdr%nsppol   =nsppol
 hdr%usepaw   =psps%usepaw
 hdr%usewvl   =usewvl !hdr%nwvlarr will be set later since the number !of wavelets have not yet been computed.
 hdr%occopt   =ebands%occopt
 hdr%codvsn   =codvsn
 hdr%date     =date
 hdr%headform =HDR_LATEST_HEADFORM ! Initialize with the latest headform
 hdr%pertcase =pertcase
 hdr%ecut     =ecut
 hdr%ecutsm   =ecutsm
 hdr%ecut_eff =ecut * (dilatmx)**2
 hdr%stmbias  =stmbias
 hdr%tphysel  =ebands%tphysel
 hdr%tsmear   =ebands%tsmear
 hdr%qptn     =qptn
 hdr%rprimd   =rprimd      ! Evolving data

!Default for other data  (all evolving data)
 hdr%etot     =1.0d20
 hdr%fermie   =1.0d20
 hdr%residm   =1.0d20

! Allocate all components of hdr
 call hdr_malloc(hdr, hdr%bantot, hdr%nkpt, hdr%nsppol, hdr%npsp, hdr%natom, hdr%ntypat,&
                 hdr%nsym, hdr%nshiftk_orig, hdr%nshiftk)

!Transfer data from ebands
 hdr%istwfk(1:nkpt) = ebands%istwfk(1:nkpt)
 hdr%kptns(:,:) = ebands%kptns(:,:)
 hdr%nband(1:nkpt*nsppol)=ebands%nband(1:nkpt*nsppol); hdr%mband = maxval(hdr%nband)
 hdr%npwarr(:) = ebands%npwarr(:)
 hdr%wtk(:) = ebands%wtk(:)

!Transfer data from psps
 hdr%pspcod    =psps%pspcod
 hdr%pspdat    =psps%pspdat
 hdr%pspso     =psps%pspso
 hdr%pspxc     =psps%pspxc
 hdr%znuclpsp  =psps%znuclpsp
 hdr%znucltypat=psps%znucltypat
 hdr%zionpsp   =psps%zionpsp
 do ipsp=1,psps%npsp
   write(hdr%title(ipsp), "(A)") psps%title(ipsp)(1:132)
 end do
 hdr%md5_pseudos = psps%md5_pseudos

 hdr%so_psp=so_psp
 hdr%symafm(1:min(size(symafm),size(hdr%symafm)))=symafm(1:min(size(symafm),size(hdr%symafm)))
 hdr%symrel(:,:,1:min(size(symrel,3),size(hdr%symrel,3))) =symrel(:,:,1:min(size(symrel,3),size(hdr%symrel,3)))
 hdr%tnons(:,1:min(size(tnons,2),size(hdr%tnons,2)))=tnons(:,1:min(size(tnons,2),size(hdr%tnons,2)))

 hdr%typat(1:natom) =typat(1:natom)  ! PMA : in tests/v2/t11 size(dtset%typat) is bigger dtset%natom
 hdr%xred(:,1:natom)=xred(:,1:natom) ! Evolving data

 hdr%kptopt = kptopt
 hdr%pawcpxocc = pawcpxocc
 hdr%nelect = nelect
 hdr%charge = charge
 hdr%kptrlatt_orig = kptrlatt_orig
 hdr%kptrlatt = kptrlatt
 hdr%shiftk_orig = shiftk_orig(:, 1:hdr%nshiftk_orig)
 hdr%shiftk = shiftk
 hdr%icoulomb = icoulomb
 hdr%amu = amu

 if (psps%usepaw==1)then
   call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,qphase_rhoij=qphase_rhoij,nspden_rhoij=nspden_rhoij,&
                             nspden=nspden,spnorb=pawspnorb,cpxocc=pawcpxocc,qpt=qptn)
   ABI_MALLOC(hdr%pawrhoij,(natom))
   ! Values of nspden/nspinor/nsppol are dummy ones; they are overwritten later (by hdr_update)
   if (present(comm_atom)) then
     if (present(mpi_atmtab)) then
       call pawrhoij_alloc(hdr%pawrhoij,cplex_rhoij,nspden_rhoij,nspinor,nsppol,typat,qphase=qphase_rhoij,&
                           pawtab=pawtab,comm_atom=comm_atom,mpi_atmtab=mpi_atmtab)
     else
       call pawrhoij_alloc(hdr%pawrhoij,cplex_rhoij,nspden_rhoij,nspinor,nsppol,typat,qphase=qphase_rhoij,&
                           pawtab=pawtab,comm_atom=comm_atom)
     end if
   else
     call pawrhoij_alloc(hdr%pawrhoij,cplex_rhoij,nspden_rhoij,nspinor,nsppol,typat,qphase=qphase_rhoij,&
                         pawtab=pawtab)
   end if
 end if

 if (psps%usepaw==1 .and. usewvl ==0 ) then
   hdr%ngfft(:) =ngfftdg(1:3)
 else if (usewvl==1) then
#if defined HAVE_BIGDFT
   hdr%ngfft(:) = (/ wvl%Glr%d%n1i, wvl%Glr%d%n2i, wvl%Glr%d%n3i /)
#else
 BIGDFT_NOTENABLED_ERROR()
#endif
 else
   hdr%ngfft(:) =ngfft(1:3)
 end if

!Transfer paw data
 if(psps%usepaw==1) then
   hdr%ecutdg   =pawecutdg
   hdr%lmn_size(1:npsp)=pawtab(1:npsp)%lmn_size
 else
   hdr%ecutdg=hdr%ecut
   hdr%lmn_size(:)=psps%lmnmax
 end if

 hdr%occ(:)=zero; idx=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     do iband=1,hdr%nband(ikpt+(isppol-1)*nkpt)
       idx=idx+1
       hdr%occ(idx)=ebands%occ(iband,ikpt,isppol)
     end do
   end do
 end do

 ABI_UNUSED(wvl%h(1))

end subroutine hdr_init_lowlvl
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_read_from_fname
!! NAME
!! hdr_read_from_fname
!!
!! FUNCTION
!!  Read the header from file fname.
!!  Use Fortran IO or Netcdf depending on the extension of the file
!!  Only rank0 process reads the header and then broadcast data to the other
!!  processes inside comm.
!!
!! INPUTS
!!  fname=String with the name of the file.
!!  comm = MPI communicator.
!!
!! OUTPUT
!!  Hdr<hdr_type>=The abinit header.
!!  fform=Kind of the array in the file (0 signals an error)
!!
!! PARENTS
!!      conducti_paw,conducti_paw_core,cut3d,emispec_paw,finddistrproc,ioprof
!!      linear_optics_paw,m_ddk,m_ioarr,m_wfd,m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_read_from_fname(Hdr,fname,fform,comm)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 integer,intent(out) :: fform
 character(len=*),intent(in) :: fname
 type(hdr_type),intent(inout) :: Hdr

!Local variables-------------------------------
 integer,parameter :: rdwr1=1,master=0
 integer :: fh,my_rank,mpierr
 character(len=500) :: msg
 character(len=len(fname)) :: my_fname

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)
 my_fname = fname

 if (nctk_try_fort_or_ncfile(my_fname, msg) /= 0) then
   MSG_ERROR(msg)
 end if

 if (my_rank == master) then
   if (.not.isncfile(my_fname)) then
     ! Use Fortran IO to open the file and read the header.
     if (open_file(my_fname,msg,newunit=fh,form="unformatted", status="old") /= 0) then
       MSG_ERROR(msg)
     end if

     call hdr_fort_read(Hdr,fh,fform,rewind=(rdwr1==1))
     ABI_CHECK(fform /= 0, sjoin("fform == 0 while reading:", my_fname))
     close(fh)

   else
     ! Use Netcdf to open the file and read the header.
#ifdef HAVE_NETCDF
     NCF_CHECK(nctk_open_read(fh, my_fname, xmpi_comm_self))
     call hdr_ncread(Hdr,fh, fform)
     ABI_CHECK(fform /= 0, sjoin("Error while reading:", my_fname))
     NCF_CHECK(nf90_close(fh))
#else
     MSG_ERROR("netcdf support not enabled")
#endif
   end if
 end if

 ! Broadcast fform and the header.
 if (xmpi_comm_size(comm) > 1) then
   call hdr%bcast(master, my_rank, comm)
   call xmpi_bcast(fform, master, comm, mpierr)
 end if

end subroutine hdr_read_from_fname
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_write_to_fname
!! NAME
!! hdr_write_to_fname
!!
!! FUNCTION
!!  Write the header and fform to file fname.
!!  Use Fortran IO or Netcdf depending on the extension of the file
!!
!! INPUTS
!!  fname=String with the name of the file.
!!  fform=Kind of the array in the file
!!  Hdr<hdr_type>=The abinit header.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_ioarr,m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_write_to_fname(Hdr,fname,fform)

!Arguments ------------------------------------
 integer,intent(in) :: fform
 character(len=*),intent(in) :: fname
 class(hdr_type),intent(inout) :: Hdr

!Local variables-------------------------------
 integer :: fh,ierr
 character(len=500) :: msg

! *************************************************************************

 if (.not.isncfile(fname)) then
   ! Use Fortran IO to write the header.
   if (open_file(fname,msg,newunit=fh,form="unformatted", status="unknown") /= 0) then
     MSG_ERROR(msg)
   end if
   call hdr%fort_write(fh, fform, ierr)
   ABI_CHECK(ierr==0, sjoin("Error while writing Abinit header to file:", fname))
   close(fh)

 else
   ! Use Netcdf to open the file and write the header.
#ifdef HAVE_NETCDF
   if (file_exists(fname)) then
     NCF_CHECK(nctk_open_modify(fh,fname, xmpi_comm_self))
   else
     NCF_CHECK_MSG(nctk_open_create(fh, fname, xmpi_comm_self), sjoin("Creating file:",  fname))
   end if

   NCF_CHECK(hdr%ncwrite(fh, fform, nc_define=.True.))
   NCF_CHECK(nf90_close(fh))
#else
   MSG_ERROR("netcdf support not enabled")
#endif
 end if

end subroutine hdr_write_to_fname
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_mpio_skip
!! NAME
!!  hdr_mio_skip
!!
!! FUNCTION
!!  Skip the abinit header in MPI-IO mode. This routine uses local MPI-IO calls hence
!!  it can be safely called by master node only. Note however that in this case the
!!  offset has to be communicated to the other nodes.
!!
!! INPUTS
!!  mpio_fh=MPI-IO file handler
!!
!! TODO
!!  We don't need to read record to skip. We just need to compute the offset from the dimensions.
!!  The algorithm is as follows:
!!
!!   1) master reads and broadcast the header.
!!   2) The offset is computed from the header
!!   3) Open the file with MPI and use the offset to point the data to be read.
!!
!! OUTPUT
!!  fform=kind of the array in the file
!!  offset=The offset of the Fortran record located immediately below the Abinit header.
!!
!! SOURCE

subroutine hdr_mpio_skip(mpio_fh, fform, offset)

!Arguments ------------------------------------
 integer,intent(in) :: mpio_fh
 integer,intent(out) :: fform
 integer(kind=XMPI_OFFSET_KIND),intent(out) :: offset

!Local variables-------------------------------
!scalars
 integer :: bsize_frm,mpi_type_frm
#ifdef HAVE_MPI_IO
 integer :: headform,ierr,mu,usepaw,npsp
!arrays
 integer(kind=MPI_OFFSET_KIND) :: fmarker,positloc
 integer :: statux(MPI_STATUS_SIZE)
#endif
 character(len=500) :: msg

! *************************************************************************

 !@hdr_type
 offset = 0; fform  = 0

 bsize_frm    = xmpio_bsize_frm    ! bsize_frm= Byte length of the Fortran record marker.
 mpi_type_frm = xmpio_mpi_type_frm ! MPI type of the record marker.

#ifdef HAVE_MPI_IO
!Reading the first record of the file -------------------------------------
!read (unitfi)   codvsn,headform,..............
 positloc = bsize_frm + 8*xmpi_bsize_ch
 call MPI_FILE_READ_AT(mpio_fh,positloc,fform,1,MPI_INTEGER,statux,ierr)

 if (ANY(fform == [1,2,51,52,101,102] )) then
   ! This is the old format !read (unitfi) codvsn,fform
   headform=22
   write(msg,'(3a,i0,4a)') &
     "ABINIT version: ",trim(abinit_version)," cannot read old files with headform: ",headform,ch10,&
     "produced by previous versions. Use an old ABINIT version to read this file or ",ch10,&
     "regenerate your files with version >= 8.0."
   MSG_ERROR(msg)

 else
   !read (unitfi)codvsn,headform,fform
   call MPI_FILE_READ_AT(mpio_fh,positloc,headform,1,MPI_INTEGER,statux,ierr)
   positloc = positloc + xmpi_bsize_int
   call MPI_FILE_READ_AT(mpio_fh,positloc,fform,1,MPI_INTEGER,statux,ierr)
 end if

 if (headform < 80) then
   write(msg,'(3a,i0,4a)') &
     "ABINIT version: ",trim(abinit_version)," cannot read old files with headform: ",headform,ch10,&
     "produced by previous versions. Use an old ABINIT version to read this file or ",ch10,&
     "regenerate your files with version >= 8.0."
   MSG_ERROR(msg)
 end if

 ! Skip first record.
 call xmpio_read_frm(mpio_fh,offset,xmpio_single,fmarker,ierr)

!Read npsp and usepaw from the second record and skip it
 positloc  = offset + bsize_frm + xmpi_bsize_int*13
 call MPI_FILE_READ_AT(mpio_fh,positloc,npsp,1,MPI_INTEGER,statux,ierr)
 positloc = positloc +  xmpi_bsize_int*4
 call MPI_FILE_READ_AT(mpio_fh,positloc,usepaw,1,MPI_INTEGER,statux,ierr)
 call xmpio_read_frm(mpio_fh,offset,xmpio_single,fmarker,ierr)

 ! Skip the rest of the file ---------------------------------------------
 do mu=1,3+npsp
   call xmpio_read_frm(mpio_fh,offset,xmpio_single,fmarker,ierr)
 end do

 if (usepaw == 1) then ! skip rhoij records.
   call xmpio_read_frm(mpio_fh,offset,xmpio_single,fmarker,ierr)
   call xmpio_read_frm(mpio_fh,offset,xmpio_single,fmarker,ierr)
 end if

#else
 MSG_ERROR("hdr_mpio_skip cannot be used when MPI-IO is not enabled")
 ABI_UNUSED(mpio_fh)
#endif

end subroutine hdr_mpio_skip
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_bsize_frecords
!! NAME
!!  hdr_bsize_frecords
!!
!! FUNCTION
!!  Compute the size of the Fortran records of the WFK file from the header and formeig.
!!
!! INPUTS
!!  Hdr<hdr_type>=The abinit header.
!!  formeig = 0 for GS WFK, 1 for response function WFK.
!!
!! OUTPUTS
!!  nfrec = Number fof Fortran records
!!  bsize_frecords(nfrec) = Byte size of each records. Allocated inside this routine.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_bsize_frecords(Hdr,formeig,nfrec,bsize_frecords)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: formeig
 integer,intent(out) :: nfrec
 class(hdr_type),intent(in) :: Hdr
!arrays
 integer(XMPI_OFFSET_KIND),allocatable,intent(out) :: bsize_frecords(:)

!Local variables-------------------------------
!scalars
 integer :: max_nfrec,ik_ibz,spin,mband,nband_k,npw_k,band
!arrays
 integer(XMPI_OFFSET_KIND),allocatable :: bsz_frec(:)

!************************************************************************

!@hdr_type
 mband = MAXVAL(Hdr%nband)
 max_nfrec = Hdr%nkpt*Hdr%nsppol * (3 + mband)

 if (formeig==1) max_nfrec = max_nfrec + Hdr%nkpt*Hdr%nsppol*mband
 ABI_MALLOC(bsz_frec, (max_nfrec))

 nfrec = 0
 do spin=1,Hdr%nsppol
   do ik_ibz=1,Hdr%nkpt
     nband_k = Hdr%nband(ik_ibz + (spin-1)*Hdr%nkpt)
     npw_k   = Hdr%npwarr(ik_ibz)

     ! First record: npw, nspinor, nband_disk
     nfrec = nfrec + 1
     bsz_frec(nfrec) = 3*xmpi_bsize_int

     ! Record with kg_k(3,npw_k) vectors
     nfrec = nfrec + 1
     bsz_frec(nfrec) = 3*npw_k*xmpi_bsize_int

     if (formeig==0) then
       ! Record with the eigenvalues
       ! eig_k(nband_k), occ_k(nband_k)
       nfrec = nfrec + 1
       bsz_frec(nfrec) = 2*nband_k*xmpi_bsize_dp

       ! cg_k record
       do band=1,nband_k
         nfrec = nfrec + 1
         bsz_frec(nfrec) = 2*npw_k*Hdr%nspinor*xmpi_bsize_dp
       end do

     else if (formeig==1) then
       do band=1,nband_k
         ! Record with the eigenvalues
         nfrec = nfrec + 1
         bsz_frec(nfrec) = 2*nband_k*xmpi_bsize_dp

         ! cg_k record
         nfrec = nfrec + 1
         bsz_frec(nfrec) = 2*npw_k*Hdr%nspinor*xmpi_bsize_dp
       end do
     else
       MSG_ERROR("Wrong formeig")
     end if

   end do
 end do

 ABI_MALLOC(bsize_frecords, (nfrec))
 bsize_frecords = bsz_frec(1:nfrec)

 ABI_FREE(bsz_frec)

end subroutine hdr_bsize_frecords
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_io_wfftype
!! NAME
!! hdr_io_wfftype
!!
!! FUNCTION
!! This subroutine deals with the I/O of the hdr_type
!! structured variables (read/write/echo).
!! According to the value of rdwr, it reads the header
!! of a file, writes it, or echo the value of the structured
!! variable to a file.
!! Note that, when reading, different records of hdr
!! are allocated here, according to the values of the
!! read variables. Records of hdr should be deallocated
!! correctly by a call to hdr_free when hdr is not used anymore.
!! Two instances of the hdr_io routines are defined :
!!  hdr_io_int to which only the unit number is given
!!  hdr_io_wfftype to which a wffil datatype is given
!!
!! INPUTS
!!  rdwr= if 1, read the hdr structured variable from the header of the file,
!!        if 2, write the header to unformatted file
!!        if 3, echo part of the header to formatted file (records 1 and 2)
!!        if 4, echo the header to formatted file
!!        if 5, read the hdr without rewinding (unformatted)
!!        if 6, write the hdr without rewinding (unformatted)
!!  unitfi=unit number of the file (unformatted if rdwr=1, 2, 5 or 6 formatted if rdwr=3,4)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  The following variables are both input or output :
!!  fform=kind of the array in the file
!!   if rdwr=1,5 : will be output ; if the reading fail, return fform=0
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!  hdr <type(hdr_type)>=the header structured variable
!!   if rdwr=1,5 : will be output
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!
!! NOTES
!! In all cases, the file is supposed to be open already
!! When reading (rdwr=1) or writing (rdwr=2), rewind the file
!! When echoing (rdwr=3) does not rewind the file.
!! When reading (rdwr=5) or writing (rdwr=6), DOES NOT rewind the file
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_io_wfftype(fform,hdr,rdwr,wff)

!Arguments ------------------------------------
 integer,intent(inout) :: fform
 integer,intent(in) :: rdwr
 type(hdr_type),intent(inout) :: hdr
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
#if defined HAVE_MPI
 integer :: ierr
#endif

! *************************************************************************

 DBG_ENTER("COLL")

 if ( wff%iomode==IO_MODE_FORTRAN .or. &
     (wff%iomode==IO_MODE_FORTRAN_MASTER .and.wff%master==wff%me).or. &
     (wff%iomode==IO_MODE_MPI  .and.wff%master==wff%me) ) then
   call hdr_io_int(fform,hdr,rdwr,wff%unwff)
   ! Master node **MUST** flush the output buffer so that the
   ! other nodes can read headform and therefore the Fortran marker length when MPI-IO is used
   if (rdwr == 2) call flush_unit(wff%unwff)
 end if

#if defined HAVE_MPI
!In the parallel case, if the files were not local, need to bcast the data
 if(rdwr==1)then
   if (wff%iomode==IO_MODE_FORTRAN_MASTER .or. wff%iomode==IO_MODE_MPI) then
     if (wff%spaceComm/=MPI_COMM_SELF) then
       call MPI_BCAST(fform,1,MPI_INTEGER,wff%master,wff%spaceComm,ierr)
       call hdr%bcast(wff%master, wff%me, wff%spaceComm)
     end if
     wff%headform=hdr%headform
     if(wff%iomode==IO_MODE_MPI)then
       call hdr_skip_wfftype(wff,ierr)
     end if
   end if
 end if
#if defined HAVE_MPI_IO
 if (rdwr == 2 .and. wff%iomode==IO_MODE_MPI) then
   if (wff%spaceComm/=MPI_COMM_SELF) then
     call xmpi_barrier(wff%spaceComm)
   end if
   wff%headform=hdr%headform
   call hdr_skip_wfftype(wff,ierr)
 end if
#endif
 if (rdwr==5) wff%headform=hdr%headform
#else
 if (rdwr==1.or.rdwr==5) wff%headform=hdr%headform
#endif

 DBG_EXIT("COLL")

end subroutine hdr_io_wfftype
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_io_int
!! NAME
!! hdr_io_int
!!
!! FUNCTION
!! This subroutine deals with the I/O of the hdr_type structured variables (read/write/echo).
!! According to the value of rdwr, it reads the header of a file, writes it, or echo the value of the structured
!! variable to a file. Note that, when reading, different records of hdr are allocated here, according to the values of the
!! read variables. Records of hdr should be deallocated correctly by a call to hdr_free when hdr is not used anymore.
!! Two instances of the hdr_io routines are defined :
!!   hdr_io_int to which only the unit number is given
!!   hdr_io_wfftype to which a wffil datatype is given
!!
!! INPUTS
!!  rdwr= if 1, read the hdr structured variable from the header of the file,
!!        if 2, write the header to unformatted file
!!        if 3, echo part of the header to formatted file (records 1 and 2)
!!        if 4, echo the header to formatted file
!!        if 5, read the hdr without rewinding (unformatted)
!!        if 6, write the hdr without rewinding (unformatted)
!!  unitfi=unit number of the file (unformatted if rdwr=1, 2, 5 or 6 formatted if rdwr=3,4)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  The following variables are both input or output :
!!  fform=kind of the array in the file
!!   if rdwr=1,5 : will be output ; if the reading fail, return fform=0
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!  hdr <type(hdr_type)>=the header structured variable
!!   if rdwr=1,5 : will be output
!!   if rdwr=2,3,4,6 : should be input, will be written or echo to file
!!
!! NOTES
!! In all cases, the file is supposed to be open already
!! When reading (rdwr=1) or writing (rdwr=2), rewind the file
!! When echoing (rdwr=3) does not rewind the file.
!! When reading (rdwr=5) or writing (rdwr=6), DOES NOT rewind the file
!!
!! PARENTS
!!      m_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_io_int(fform,hdr,rdwr,unitfi)

!Arguments ------------------------------------
 integer,intent(inout) :: fform
 integer,intent(in) :: rdwr,unitfi
 type(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------
 integer :: ierr

!*************************************************************************

 DBG_ENTER("COLL")

 select case(rdwr)
 case (1, 5)
   ! Reading the header of an unformatted file
    call hdr_fort_read(Hdr,unitfi,fform,rewind=(rdwr==1))

 case (2, 6)
   ! Writing the header of an unformatted file
   call hdr%fort_write(unitfi, fform, ierr, rewind=(rdwr==2))

 case (3, 4)
   !  Writing the header of a formatted file
   call hdr_echo(Hdr,fform,rdwr,unit=unitfi)
 case default
   MSG_ERROR(sjoin("Wrong value for rdwr: ",itoa(rdwr)))
 end select

 DBG_EXIT("COLL")

end subroutine hdr_io_int
!!***
!----------------------------------------------------------------------

!!****f* m_hdr/hdr_echo
!! NAME
!! hdr_echo
!!
!! FUNCTION
!! Echo the header
!!
!! INPUTS
!!  hdr <type(hdr_type)>=the header structured variable
!!  rdwr= if 3, echo part of the header to formatted file (records 1 and 2)
!!        if 4, echo the header to formatted file
!!  fform=kind of the array in the file
!!  [unit]=unit number of the formatted file [DEFAULT: std_out]
!!  [header]=Optional title.
!!
!! OUTPUT
!!  Only writing
!!
!! TODO
!!   Activate new header, avoid printing tons of lines with occupations.
!!
!! PARENTS
!!      cut3d,initaim,ioprof,m_ddk,m_dvdb,m_hdr,m_wfd,m_wfk,mrggkk,rchkgsheader
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_echo(hdr, fform, rdwr, unit, header)

!Arguments ------------------------------------
 integer,intent(inout) :: fform
 integer,intent(in) :: rdwr
 integer,optional,intent(in) :: unit
 class(hdr_type),intent(inout) :: hdr
 character(len=*),optional,intent(in) :: header

!Local variables-------------------------------
 integer,parameter :: max_ns=6
 integer :: iatom,ii,ikpt,ipsp,isym,ount !,ns
 !character(len=500) :: msg

!*************************************************************************

 ount = std_out; if (present(unit)) ount = unit; if (ount == dev_null) return

 write(ount,'(a)')' ==============================================================================='
 if (present(header)) write(ount, "(a)")ch10//' === '//trim(adjustl(header))//' === '
 if (rdwr==3) write(ount, '(a)' ) ' ECHO of part of the ABINIT file header '
 if (rdwr==4) write(ount, '(a)' ) ' ECHO of the ABINIT file header '
 write(ount, '(a)' ) ' '
 write(ount, '(a)' ) ' First record :'
 write(ount, '(a,a8,2i5)' )  '.codvsn,headform,fform = ',hdr%codvsn, hdr%headform, fform
 write(ount, '(a)' ) ' '
 write(ount, '(a)' ) ' Second record :'
 write(ount, '(a,4i6)') ' bantot,intxc,ixc,natom  =',hdr%bantot, hdr%intxc, hdr%ixc, hdr%natom
 write(ount, '(a,4i6)') ' ngfft(1:3),nkpt         =',hdr%ngfft(1:3), hdr%nkpt
 write(ount, '(a,2i6)') ' nspden,nspinor          =',hdr%nspden, hdr%nspinor
 write(ount, '(a,4i6)' ) ' nsppol,nsym,npsp,ntypat =',hdr%nsppol,hdr%nsym,hdr%npsp,hdr%ntypat
 write(ount, '(a,3i6)' ) ' occopt,pertcase,usepaw  =',hdr%occopt,hdr%pertcase,hdr%usepaw
 write(ount, '(a,3es18.10)') ' ecut,ecutdg,ecutsm      =',hdr%ecut, hdr%ecutdg, hdr%ecutsm
 write(ount, '(a, es18.10)' ) ' ecut_eff                =',hdr%ecut_eff
 write(ount, '(a,3es18.10)') ' qptn(1:3)               =',hdr%qptn(1:3)
 write(ount, '(a,3es18.10)' ) ' rprimd(1:3,1)           =',hdr%rprimd(1:3,1)
 write(ount, '(a,3es18.10)' ) ' rprimd(1:3,2)           =',hdr%rprimd(1:3,2)
 write(ount, '(a,3es18.10)' ) ' rprimd(1:3,3)           =',hdr%rprimd(1:3,3)
 write(ount, '(a,3es18.10)') ' stmbias,tphysel,tsmear  =',hdr%stmbias,hdr%tphysel, hdr%tsmear

#ifdef DEV_NEW_HDR
 write(ount, "(a,2es18.10,i0)") ' nelect,charge,icoulomb  =',hdr%nelect, hdr%charge, hdr%icoulomb
 write(ount, "(a,2i6)")         ' kptopt,pawcpxocc        =',hdr%kptopt, hdr%pawcpxocc
 write(ount, '(a,9(i0,1x))')    ' kptrlatt_orig           = ',hdr%kptrlatt_orig
 write(ount, '(a,9(i0,1x))' )   ' kptrlatt                = ',hdr%kptrlatt

 ns = min(size(hdr%shiftk_orig, dim=2), max_ns)
 write(msg, sjoin("(a,",itoa(3*ns),"(f4.2,1x))")) ' shiftk_orig             = ',hdr%shiftk_orig(:,1:ns)
 if (size(hdr%shiftk_orig, dim=2) > max_ns) msg = sjoin(msg, "...")
 write(ount,"(a)")trim(msg)

 ns = min(size(hdr%shiftk, dim=2), max_ns)
 write(msg, sjoin("(a,",itoa(3*ns),"(f4.2,1x))")) ' shiftk                  = ',hdr%shiftk(:,1:ns)
 if (size(hdr%shiftk, dim=2) > max_ns) msg = sjoin(msg, "...")
 write(ount,"(a)")trim(msg)
#endif

 write(ount, '(a)' )
 if (rdwr==3)then
   write(ount, '(a,i3,a)' ) ' The header contain ',hdr%npsp+2,' additional records.'
 else
   write(ount, '(a)' ) ' Third record :'
   write(ount, '(a,(12i4,8x))') ' istwfk=',hdr%istwfk
   write(ount, '(a,(12i4,8x))') ' nband =',hdr%nband
   write(ount, '(a,(10i5,8x))') ' npwarr=',hdr%npwarr

   write(ount, '(a,(12i4,8x))') ' so_psp=',hdr%so_psp(:)
   !write(ount,'(a,(12f6.2,1x))' )' amu   =',hdr%amu

   write(ount, '(a)') ' symafm='
   write(ount, '(8x,24i3,8x)') hdr%symafm

   write(ount, '(a)' ) ' symrel='
   do isym=1,hdr%nsym/2
     write(ount, '(a,9i4,a,9i4)' ) '        ',hdr%symrel(:,:,2*isym-1),'  ',hdr%symrel(:,:,2*isym)
   end do
   if(2*(hdr%nsym/2)/=hdr%nsym)write(ount, '(a,9i4)' ) '        ',hdr%symrel(:,:,hdr%nsym)

   write(ount, '(a,(12i4,8x))') ' type  =',hdr%typat(:)
   write(ount, '(a)' ) ' kptns =                 (max 50 k-points will be written)'
   do ikpt=1,min(hdr%nkpt,50)
     write(ount, '(a,3es16.6)' ) '        ',hdr%kptns(:,ikpt)
   end do
   write(ount, '(a)' ) ' wtk ='
   do ikpt=1,hdr%nkpt,10
     write(ount, '(a,10f6.2)' ) '        ',hdr%wtk(ikpt:min(hdr%nkpt,ikpt + 10 - 1))
   end do
   write(ount, '(a)' ) '   occ ='
   do ii=1,hdr%bantot,10
     write(ount, '(a,10f6.2)') '        ',hdr%occ(ii:min(hdr%bantot,ii+10-1))
   end do
   write(ount, '(a)' ) ' tnons ='
   do isym=1,hdr%nsym/2
     write(ount, '(a,3f10.6,a,3f10.6)' ) '        ',hdr%tnons(:,2*isym-1),'  ',hdr%tnons(:,2*isym)
   end do
   if(2*(hdr%nsym/2)/=hdr%nsym)write(ount, '(a,3f10.6)' ) '        ',hdr%tnons(:,hdr%nsym)
   write(ount, '(a,(10f6.2,8x))') '  znucl=',hdr%znucltypat(:)
   write(ount,'(a)')

   write(ount, '(a)' ) ' Pseudopotential info :'
   do ipsp=1,hdr%npsp
     write(ount,'(a,a)' ) ' title=',trim(hdr%title(ipsp))
     ! TODO: This part should always be printed.
     !write(ount,'(a,a)' ) '   md5=',trim(hdr%md5_pseudos(ipsp))
     write(ount,'(a,f6.2,a,f6.2,a,i3,a,i6,a,i3,a,i3)' ) &
      '  znuclpsp=',hdr%znuclpsp(ipsp),    ', zionpsp=',  hdr%zionpsp(ipsp),&
      ', pspso=' , hdr%pspso(ipsp),  ', pspdat=',hdr%pspdat(ipsp),          &
      ', pspcod=', hdr%pspcod(ipsp), ', pspxc=', hdr%pspxc(ipsp)

     if(hdr%usepaw==1)then
       write(ount,'(a,i3)' ) '  lmn_size=', hdr%lmn_size(ipsp)
     else
       write(ount,'(a,i3)' ) '  lmnmax  =', hdr%lmn_size(ipsp)
     end if
   end do

   write(ount, '(a)' ) ' '
   write(ount, '(a)' ) ' Last record :'
   write(ount, '(a,es16.6,es22.12,es16.6)' )' residm,etot,fermie=',hdr%residm, hdr%etot, hdr%fermie
   write(ount, '(a)' ) ' xred ='
   do iatom=1,hdr%natom
     write(ount, '(a,3es16.6)' ) '        ',hdr%xred(:,iatom)
   end do

   if (hdr%usepaw==1)then
     call pawrhoij_io(hdr%pawrhoij,ount,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,hdr%headform,"Echo")
   end if

   if (rdwr==3)write(ount, '(a)' ) ' End the ECHO of part of the ABINIT file header '
   if (rdwr==4)write(ount, '(a)' ) ' End the ECHO of the ABINIT file header '
   write(ount,'(a)')' ==============================================================================='
 end if ! rdwr is 3 or 4

 call flush_unit(ount)

end subroutine hdr_echo
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_skip_int
!! NAME
!! hdr_skip_int
!!
!! FUNCTION
!! Skip wavefunction or density file header, after having rewound the file.
!! Two instances of the hdr_skip routines are defined:
!!  hdr_skip_int to which only the unit number is given
!!  hdr_skip_wfftype to which a wffil datatype is given
!!
!! INPUTS
!!  unit = number of unit to be read
!!
!! OUTPUT
!!  ierr = error code returned by the MPI calls
!!
!! SIDE EFFECTS
!!
!! NOTES
!! No checking performed, since hdr_skip is assumed to be used only
!! on temporary wavefunction files.
!! This initialize further reading and checking by rwwf
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_skip_int(unitfi,ierr)

!Arguments ------------------------------------
 integer,intent(in) :: unitfi
 integer,intent(out) :: ierr

!Local variables-------------------------------
 type(wffile_type) :: wff

! *************************************************************************

!Use default values for wff
 wff%unwff=unitfi; wff%iomode=IO_MODE_FORTRAN
 wff%me=0; wff%master=0
!Then, transmit to hdr_skip_wfftype
 call hdr_skip_wfftype(wff,ierr)

end subroutine hdr_skip_int
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_skip_wfftype
!! NAME
!! hdr_skip_wfftype
!!
!! FUNCTION
!! Skip wavefunction or density file header, after having rewound the file.
!! Two instances of the hdr_skip routines are defined :
!!  hdr_skip_int to which only the unit number is given
!!  hdr_skip_wfftype to which a wffil datatype is given
!!
!! INPUTS
!!  unit = number of unit to be read
!!
!! OUTPUT
!!  ierr = error code returned by the MPI calls
!!
!! NOTES
!! No checking performed, since hdr_skip is assumed to be used only
!! on temporary wavefunction files.
!! This initialize further reading and checking by rwwf
!!
!! PARENTS
!!      m_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_skip_wfftype(wff,ierr)

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer, intent(out) :: ierr

!Local variables-------------------------------
 integer :: headform,mu,npsp,unit,usepaw,fform
 integer :: integers(17)
 character(len=8) :: codvsn
 character(len=500) :: msg,errmsg
#if defined HAVE_MPI_IO
 integer(kind=MPI_OFFSET_KIND) :: delim_record,posit,positloc
 integer :: statux(MPI_STATUS_SIZE)
#endif

!*************************************************************************

 !@hdr_type
 unit=wff%unwff; ierr=0

 if( wff%iomode==IO_MODE_FORTRAN .or. (wff%iomode==IO_MODE_FORTRAN_MASTER.and.wff%master==wff%me) ) then

   rewind(unit, err=10, iomsg=errmsg)

   ! Pick off headform from WF file. Support for pre-v9 (length of codvsn was changed from 6 to 8) is implemented.
   ABI_CHECK(read_first_record(unit, codvsn, headform, fform, errmsg) == 0, errmsg)

   if (headform==1   .or. headform==2   .or. &
       headform==51  .or. headform==52  .or. &
       headform==101 .or. headform==102 ) headform=22

   if (headform < 80) then
     write(msg,'(3a,i0,4a)')&
       "ABINIT version: ",trim(abinit_version)," cannot read old files with headform: ",headform,ch10,&
       "produced by previous versions. Use an old ABINIT version to read this file or ",ch10,&
       "regenerate your files with version >= 8.0."
     MSG_ERROR(msg)
   end if

   read (unit, err=10, iomsg=errmsg) integers(1:13),npsp,integers(15:17),usepaw

!  Skip rest of header records
   do mu=1,3+npsp
     read (unit, err=10, iomsg=errmsg)
   end do

   if (usepaw==1) then
     read (unit, err=10, iomsg=errmsg)
     read (unit, err=10, iomsg=errmsg)
   end if

#if defined HAVE_MPI_IO
 else if(wff%iomode==IO_MODE_MPI)then

   headform=wff%headform
   if (headform==1   .or. headform==2   .or. &
      headform==51  .or. headform==52  .or. &
      headform==101 .or. headform==102) headform=22

   if (headform < 80) then
     write(msg,'(3a,i0,4a)')&
       "ABINIT version: ",trim(abinit_version)," cannot read old files with headform: ",headform,ch10,&
       "produced by previous versions. Use an old ABINIT version to read this file or ",ch10,&
       "regenerate your files with version >= 8.0."
     MSG_ERROR(msg)
   end if

!  Causes all previous writes to be transferred to the storage device
   call flush_unit(wff%unwff)
   call MPI_FILE_SYNC(wff%fhwff,ierr)

!  Check FORTRAN record marker length (only at first call)
   if (wff%nbOct_recMarker<=0) then
     call getRecordMarkerLength_wffile(wff)
   end if

   if (wff%master==wff%me) then

!    Reading the first record of the file -------------------------------------
!    read (unitfi)   codvsn,headform,..............
     posit = 0
     call rwRecordMarker(1,posit,delim_record,wff,ierr)

!    Reading the second record of the file ------------------------------------
!    read(unitfi) bantot, hdr%date, hdr%intxc.................
!    Pick off npsp and usepaw from WF file
     positloc  = posit + wff%nbOct_recMarker + wff%nbOct_int*13
     call MPI_FILE_READ_AT(wff%fhwff,positloc,npsp,1,MPI_INTEGER,statux,ierr)

     ! Read usepaw and skip the fortran record
     positloc = positloc +  wff%nbOct_int*4
     call MPI_FILE_READ_AT(wff%fhwff,positloc,usepaw,1,MPI_INTEGER,statux,ierr)
     call rwRecordMarker(1,posit,delim_record,wff,ierr)

     ! Skip the rest of the file ---------------------------------------------
     do mu=1,3+npsp
       call rwRecordMarker(1,posit,delim_record,wff,ierr)
     end do

     if (usepaw==1) then
       call rwRecordMarker(1,posit,delim_record,wff,ierr)
       call rwRecordMarker(1,posit,delim_record,wff,ierr)
     end if

     wff%offwff=posit
   end if

   if (wff%spaceComm/=MPI_COMM_SELF) then
     call MPI_BCAST(wff%offwff,1,wff%offset_mpi_type,wff%master,wff%spaceComm,ierr)
   end if
#endif
 end if

 ! Handle IO-error: write warning and let the caller handle the exception.
 return
10 ierr=1
 MSG_WARNING(errmsg)

end subroutine hdr_skip_wfftype
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_update
!! NAME
!! hdr_update
!!
!! FUNCTION
!! This subroutine update the header structured datatype.
!! Most of its records had been initialized correctly, but some corresponds
!! to evolving variables, or change with the context (like fform),
!! This routine is to be called before writing the header
!! to a file, in order to have up-to-date information.
!!
!! INPUTS
!! bantot=total number of bands
!! etot=total energy (Hartree)
!! fermie=Fermi energy (Hartree)
!! mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!! comm_atom=--optional-- MPI communicator over atoms
!! residm=maximal residual
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! occ(bantot)=occupancies for each band and k point
!! pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!! xred(3,natom)= relative coords of atoms in unit cell (dimensionless)
!! amu(ntypat)=masses in atomic mass units for each kind of atom in cell.
!!
!! OUTPUT
!! hdr <type(hdr_type)>=the header, initialized, and for most part of
!!   it, contain its definite values, except for evolving variables
!!
!! PARENTS
!!      afterscfloop,dfpt_looppert,dfpt_scfcv,gstate,nonlinear,respfn,scfcv
!!      setup_bse,setup_screening,setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_update(hdr,bantot,etot,fermie,residm,rprimd,occ,pawrhoij,xred,amu, &
                      comm_atom,mpi_atmtab) ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bantot
 integer,optional,intent(in) :: comm_atom
 real(dp),intent(in) :: etot,fermie,residm
 class(hdr_type),intent(inout) :: hdr
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: occ(bantot),rprimd(3,3),xred(3,hdr%natom),amu(hdr%ntypat)
 type(pawrhoij_type),intent(inout) :: pawrhoij(:)

! *************************************************************************

 !@hdr_type
!Update of the "evolving" data
 hdr%etot     =etot
 hdr%fermie   =fermie
 hdr%residm   =residm
 hdr%rprimd(:,:)=rprimd(:,:)
 hdr%occ(:)   =occ(:)
 hdr%xred(:,:)=xred(:,:)
 hdr%amu(:) = amu

 if (hdr%usepaw==1) then
   if (present(comm_atom)) then
     if (present(mpi_atmtab)) then
       call pawrhoij_copy(pawrhoij,hdr%pawrhoij,comm_atom=comm_atom,mpi_atmtab=mpi_atmtab)
     else
       call pawrhoij_copy(pawrhoij,hdr%pawrhoij,comm_atom=comm_atom)
     end if
   else
     call pawrhoij_copy(pawrhoij,hdr%pawrhoij)
   end if
 end if

end subroutine hdr_update
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_bcast
!! NAME
!! hdr_bcast
!!
!! FUNCTION
!! This subroutine transmits the header structured datatype
!! initialized on one processor (or a group of processor),
!! to the other processors. It also allocate the needed
!! part of the header.
!!
!! INPUTS
!!  master = id of the master process
!!  me = id of the current process
!!  comm = id of the space communicator handler
!!
!! OUTPUT
!!  (no output)
!!
!! SIDE EFFECTS
!!  hdr <type(hdr_type)>=the header. For the master, it is already
!!   initialized entirely, while for the other procs, everything has
!!   to be transmitted.
!!
!! NOTES
!! This routine is called only in the case of MPI version of the code.
!!
!! PARENTS
!!      elphon,initaim,m_dvdb,m_hdr,m_io_screening,m_ioarr,m_wfk,optic,read_gkk
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_bcast(hdr, master, me, comm)

!Arguments ------------------------------------
 integer, intent(in) :: master,me,comm
 class(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------
!scalars
 integer :: bantot,cplex_rhoij,iatom,ierr,index,index2,ipsp,iq,iq0,ispden,list_size,list_size2
 integer :: lmn2_size,natom,nkpt,npsp,nsel,nspden,nsppol,nsym,nrhoij,ntypat,qphase
 character(len=fnlen) :: list_tmp
!arrays
 integer,allocatable :: list_int(:)
 real(dp),allocatable :: list_dpr(:)
 character(len=fnlen),allocatable :: list_char(:)

! *************************************************************************

 !@hdr_type
 if (xmpi_comm_size(comm) == 1) return ! Nothing to do

 DBG_ENTER("COLL")

!Transmit the integer scalars
 list_size = 43
 ABI_MALLOC(list_int,(list_size))
 if (master==me)then
   list_int(1)=hdr%bantot
   list_int(2)=hdr%date
   list_int(3)=hdr%headform
   list_int(4)=hdr%intxc
   list_int(5)=hdr%ixc
   list_int(6)=hdr%natom
   list_int(7)=hdr%nkpt
   list_int(8)=hdr%npsp
   list_int(9)=hdr%nspden
   list_int(10)=hdr%nspinor
   list_int(11)=hdr%nsppol
   list_int(12)=hdr%nsym
   list_int(13)=hdr%ntypat
   list_int(14)=hdr%occopt
   list_int(15)=hdr%pertcase
   list_int(16)=hdr%usepaw
   list_int(17:19)=hdr%ngfft(1:3)
   list_int(20)=hdr%usewvl
   list_int(21)=hdr%kptopt
   list_int(22)=hdr%pawcpxocc
   list_int(23)=hdr%nshiftk_orig
   list_int(24)=hdr%nshiftk
   list_int(25:33)=reshape(hdr%kptrlatt_orig, [3*3])
   list_int(34:42)=reshape(hdr%kptrlatt, [3*3])
   list_int(43)=hdr%icoulomb
 end if

 call xmpi_bcast(list_int,master,comm,ierr)

 if(master/=me)then
   hdr%bantot  =list_int(1)
   hdr%date    =list_int(2)
   hdr%headform=list_int(3)
   hdr%intxc   =list_int(4)
   hdr%ixc     =list_int(5)
   hdr%natom   =list_int(6)
   hdr%nkpt    =list_int(7)
   hdr%npsp    =list_int(8)
   hdr%nspden  =list_int(9)
   hdr%nspinor =list_int(10)
   hdr%nsppol  =list_int(11)
   hdr%nsym    =list_int(12)
   hdr%ntypat  =list_int(13)
   hdr%occopt  =list_int(14)
   hdr%pertcase=list_int(15)
   hdr%usepaw  =list_int(16)
   hdr%ngfft(1:3)=list_int(17:19)
   hdr%usewvl  =list_int(20)
   hdr%kptopt       = list_int(21)
   hdr%pawcpxocc    = list_int(22)
   hdr%nshiftk_orig = list_int(23)
   hdr%nshiftk      = list_int(24)
   hdr%kptrlatt_orig = reshape(list_int(25:33), [3,3])
   hdr%kptrlatt = reshape(list_int(34:42), [3,3])
   hdr%icoulomb = list_int(43)
 end if
 ABI_FREE(list_int)

 bantot=hdr%bantot
 natom =hdr%natom
 nkpt  =hdr%nkpt
 npsp  =hdr%npsp
 nspden=hdr%nspden
 nsppol=hdr%nsppol
 nsym  =hdr%nsym
 ntypat=hdr%ntypat

 if (master/=me) then
!  Allocate all components of hdr
   call hdr_malloc(hdr, bantot, nkpt, nsppol, npsp, natom, ntypat,&
                   nsym, hdr%nshiftk_orig, hdr%nshiftk)
 end if

!Transmit the integer arrays
 list_size=nkpt*(2+nsppol)+6*npsp+10*nsym+natom
 ABI_MALLOC(list_int,(list_size))
 if (master==me)then
   list_int(1      :nkpt             )=hdr%istwfk ; index=nkpt
   list_int(1+index:nkpt*nsppol+index)=hdr%nband  ; index=index+nkpt*nsppol
   list_int(1+index:nkpt       +index)=hdr%npwarr ; index=index+nkpt
   list_int(1+index:npsp       +index)=hdr%pspcod ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%pspdat ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%pspso  ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%pspxc  ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%lmn_size ; index=index+npsp
   list_int(1+index:npsp       +index)=hdr%so_psp ; index=index+npsp
   list_int(1+index:nsym       +index)=hdr%symafm ; index=index+nsym
   list_int(1+index:nsym*3*3   +index)=reshape(hdr%symrel,(/3*3*nsym/))
   index=index+nsym*3*3
   list_int(1+index:natom      +index)=hdr%typat   ; index=index+natom
 end if

 call xmpi_bcast(list_int,master,comm,ierr)

 if(master/=me)then
   hdr%istwfk=list_int(1      :nkpt             ) ; index=nkpt
   hdr%nband =list_int(1+index:nkpt*nsppol+index) ; index=index+nkpt*nsppol
   hdr%npwarr=list_int(1+index:nkpt       +index) ; index=index+nkpt
   hdr%pspcod=list_int(1+index:npsp       +index) ; index=index+npsp
   hdr%pspdat=list_int(1+index:npsp       +index) ; index=index+npsp
   hdr%pspso =list_int(1+index:npsp       +index) ; index=index+npsp
   hdr%pspxc =list_int(1+index:npsp       +index) ; index=index+npsp
   hdr%lmn_size=list_int(1+index:npsp     +index) ; index=index+npsp
   hdr%so_psp =list_int(1+index:npsp   +index) ; index=index+npsp
   hdr%symafm=list_int(1+index:nsym       +index) ; index=index+nsym
   hdr%symrel=reshape(list_int(1+index:nsym*3*3   +index),(/3,3,nsym/))
   index=index+nsym*3*3
   hdr%typat  =list_int(1+index:natom      +index) ; index=index+natom
 end if
 ABI_FREE(list_int)

!Transmit the double precision scalars and arrays
 list_size = 21+ 3*nkpt+nkpt+bantot + 3*nsym + 3*natom + 2*npsp+ntypat + &
             2 + 3*hdr%nshiftk_orig + 3*hdr%nshiftk + hdr%ntypat
 ABI_MALLOC(list_dpr,(list_size))

 if (master==me)then
   list_dpr(1)=hdr%ecut_eff
   list_dpr(2)=hdr%etot
   list_dpr(3)=hdr%fermie
   list_dpr(4)=hdr%residm
   list_dpr(5:13)=reshape(hdr%rprimd(1:3,1:3),(/9/))
   list_dpr(14)=hdr%ecut
   list_dpr(15)=hdr%ecutdg
   list_dpr(16)=hdr%ecutsm
   list_dpr(17)=hdr%tphysel
   list_dpr(18)=hdr%tsmear
   list_dpr(19:21)=hdr%qptn(1:3)                                 ; index=21
   list_dpr(1+index:3*nkpt +index)=reshape(hdr%kptns,(/3*nkpt/)) ; index=index+3*nkpt
   list_dpr(1+index:nkpt   +index)=hdr%wtk                       ; index=index+nkpt
   list_dpr(1+index:bantot +index)=hdr%occ                       ; index=index+bantot
   list_dpr(1+index:3*nsym +index)=reshape(hdr%tnons,(/3*nsym/)) ; index=index+3*nsym
   list_dpr(1+index:3*natom+index)=reshape(hdr%xred,(/3*natom/)) ; index=index+3*natom
   list_dpr(1+index:npsp   +index)=hdr%zionpsp                   ; index=index+npsp
   list_dpr(1+index:npsp   +index)=hdr%znuclpsp                  ; index=index+npsp
   list_dpr(1+index:ntypat  +index)=hdr%znucltypat               ; index=index+ntypat
   list_dpr(1+index)=hdr%nelect; index=index+1
   list_dpr(1+index)=hdr%charge; index=index+1
   list_dpr(1+index:index+3*hdr%nshiftk_orig) = reshape(hdr%shiftk_orig, [3*hdr%nshiftk_orig])
   index=index+3*hdr%nshiftk_orig
   list_dpr(1+index:index+3*hdr%nshiftk) = reshape(hdr%shiftk, [3*hdr%nshiftk])
   index=index+3*hdr%nshiftk
   list_dpr(1+index:index+hdr%ntypat) = hdr%amu(1:hdr%ntypat)
 end if

 call xmpi_bcast(list_dpr,master,comm,ierr)

 if(master/=me)then
   hdr%ecut_eff=list_dpr(1)
   hdr%etot    =list_dpr(2)
   hdr%fermie  =list_dpr(3)
   hdr%residm  =list_dpr(4)
   hdr%rprimd  =reshape(list_dpr(5:13),(/3,3/))
   hdr%ecut    =list_dpr(14)
   hdr%ecutdg  =list_dpr(15)
   hdr%ecutsm  =list_dpr(16)
   hdr%tphysel =list_dpr(17)
   hdr%tsmear  =list_dpr(18)
   hdr%qptn(1:3)=list_dpr(19:21)                                    ; index=21
   hdr%kptns   =reshape(list_dpr(1+index:3*nkpt +index),(/3,nkpt/)) ; index=index+3*nkpt
   hdr%wtk     =list_dpr(1+index:nkpt   +index)                     ; index=index+nkpt
   hdr%occ     =list_dpr(1+index:bantot +index)                     ; index=index+bantot
   hdr%tnons   =reshape(list_dpr(1+index:3*nsym +index),(/3,nsym/)) ; index=index+3*nsym
   hdr%xred    =reshape(list_dpr(1+index:3*natom+index),(/3,natom/)); index=index+3*natom
   hdr%zionpsp =list_dpr(1+index:npsp   +index)                     ; index=index+npsp
   hdr%znuclpsp=list_dpr(1+index:npsp   +index)                     ; index=index+npsp
   hdr%znucltypat=list_dpr(1+index:ntypat  +index)                  ; index=index+ntypat
   hdr%nelect = list_dpr(1+index); index=index+1
   hdr%charge = list_dpr(1+index); index=index+1
   hdr%shiftk_orig = reshape(list_dpr(1+index:index+3*hdr%nshiftk_orig), [3, hdr%nshiftk_orig])
   index=index+3*hdr%nshiftk_orig
   hdr%shiftk = reshape(list_dpr(1+index:index+3*hdr%nshiftk), [3, hdr%nshiftk])
   index=index+3*hdr%nshiftk
   hdr%amu = list_dpr(1+index:index+hdr%ntypat)
 end if
 ABI_FREE(list_dpr)

!Transmit the characters
 list_size=npsp+1 + npsp
 ABI_MALLOC(list_char,(list_size))
 if (master==me)then
   list_char(1)       =hdr%codvsn  ! Only 8 characters are stored in list_char(1)
   list_char(2:npsp+1)=hdr%title
   list_char(npsp+2:) =hdr%md5_pseudos
 end if

 call xmpi_bcast(list_char,master,comm,ierr)

 if(master/=me)then
   list_tmp=list_char(1)
   hdr%codvsn=list_tmp(1:8)
   do ipsp=2,npsp+1
     list_tmp =list_char(ipsp)
     hdr%title(ipsp-1) =list_tmp(1:min(fnlen,132))
   end do
   do ipsp=npsp+2,2*npsp+1
     hdr%md5_pseudos(ipsp-npsp-1) = list_char(ipsp)(1:md5_slen)
   end do
 end if
 ABI_FREE(list_char)

!Transmit the structured variables in case of PAW
 if (hdr%usepaw==1) then

   nrhoij=0
   if (master==me)then
     cplex_rhoij=hdr%pawrhoij(1)%cplex_rhoij
     qphase=hdr%pawrhoij(1)%qphase
     nspden=hdr%pawrhoij(1)%nspden
     do iatom=1,natom
       nrhoij=nrhoij+hdr%pawrhoij(iatom)%nrhoijsel
     end do
   end if

   call xmpi_bcast(nrhoij,master,comm,ierr)
   call xmpi_bcast(cplex_rhoij,master,comm,ierr)
   call xmpi_bcast(qphase,master,comm,ierr)
   call xmpi_bcast(nspden,master,comm,ierr)

   list_size=natom+nrhoij;list_size2=nspden*nrhoij*cplex_rhoij*qphase
   ABI_MALLOC(list_int,(list_size))
   ABI_MALLOC(list_dpr,(list_size2))
   if (master==me)then
     index=0;index2=0
     do iatom=1,natom
       nsel=hdr%pawrhoij(iatom)%nrhoijsel
       lmn2_size=hdr%pawrhoij(iatom)%lmn2_size
       list_int(1+index)=nsel
       list_int(2+index:1+nsel+index)=hdr%pawrhoij(iatom)%rhoijselect(1:nsel)
       index=index+1+nsel
       do ispden=1,nspden
         do iq=1,qphase
           iq0=merge(0,lmn2_size*cplex_rhoij,iq==1)
           list_dpr(1+index2:nsel*cplex_rhoij+index2)=hdr%pawrhoij(iatom)%rhoijp(iq0+1:iq0+nsel*cplex_rhoij,ispden)
           index2=index2+nsel*cplex_rhoij
         end do
       end do
     end do
   end if

   call xmpi_bcast(list_int,master,comm,ierr)
   call xmpi_bcast(list_dpr,master,comm,ierr)

   if(master/=me)then
     index=0;index2=0
     ABI_MALLOC(hdr%pawrhoij,(natom))
     call pawrhoij_alloc(hdr%pawrhoij,cplex_rhoij,nspden,hdr%nspinor,hdr%nsppol,hdr%typat,&
                         lmnsize=hdr%lmn_size,qphase=qphase)
     do iatom=1,natom
       nsel=list_int(1+index)
       lmn2_size=hdr%pawrhoij(iatom)%lmn2_size
       hdr%pawrhoij(iatom)%nrhoijsel=nsel
       hdr%pawrhoij(iatom)%rhoijselect(1:nsel)=list_int(2+index:1+nsel+index)
       index=index+1+nsel
       do ispden=1,nspden
         do iq=1,qphase
           iq0=merge(0,lmn2_size*cplex_rhoij,iq==1)
           hdr%pawrhoij(iatom)%rhoijp(iq0+1:iq0+nsel*cplex_rhoij,ispden)=list_dpr(1+index2:nsel*cplex_rhoij+index2)
           index2=index2+nsel*cplex_rhoij
         end do
       end do
     end do
   end if
   ABI_FREE(list_int)
   ABI_FREE(list_dpr)
 end if

 hdr%mband = maxval(hdr%nband)

 DBG_EXIT("COLL")

end subroutine hdr_bcast
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/read_first_record
!! NAME
!! read_first_record
!!
!! FUNCTION
!!  Read the first record of the header.
!!  This function is neede to support for pre-Abinitv9 headers:
!!  length of codvsn was changed from 6 to 8 in v9
!!
!! SOURCE

integer function read_first_record(unit, codvsn8, headform, fform, errmsg) result(ierr)

!Arguments ------------------------------------
 integer,intent(in) :: unit
 integer,intent(out) :: headform, fform
 character(len=8),intent(out) :: codvsn8
 character(len=*),intent(out) :: errmsg

!Local variables-------------------------------
 integer :: major, ii
 character(len=6) :: codvsn6

!*************************************************************************

 ! Try pre-v9 first. This read should not fail as we have enough space in the record
 ! Obviously headform and fform are wrong in > Abinit9.
 read(unit, iostat=ierr, iomsg=errmsg) codvsn6, headform, fform
 if (ierr /= 0) then
   call wrtout(std_out, "Fatal error while reading the first record of the Abinit header!")
   return
 end if

 ii = index(codvsn6, ".")
 if (ii == 0 .or. ii == 1) then
   errmsg = sjoin("Cannot find major.minor pattern in codvsn:", codvsn6)
   ierr = 1; return
 end if

 major = atoi(codvsn6(:ii-1))
 !call wrtout(std_out, sjoin("Reading HDR file generated by major version:", itoa(major)))
 if (major > 8) then
   backspace(unit)
   read(unit, iostat=ierr, iomsg=errmsg) codvsn8, headform, fform
   if (ierr /= 0) then
     call wrtout(std_out, "Fatal error while reading the first record of the Abinit header version > 8!")
     return
   end if
 else
   codvsn8 = ""
   codvsn8(1:6) = codvsn6
 end if

end function read_first_record
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_fort_read
!! NAME
!! hdr_fort_read
!!
!! FUNCTION
!! Reads the header from a logical unit associated to a unformatted file.
!! Note that, when reading, different records of hdr are allocated here, according to the values of the
!! read variables. Records of hdr should be deallocated correctly by a call to hdr_free when hdr is not used anymore.
!!
!! INPUTS
!!  unit=unit number of the unformatted file
!!  [rewind]=True to rewind the file. Default: False
!!
!! OUTPUT
!!  Hdr<hdr_type>=The header of the file fully initialized (if fform /=0)
!!  fform=kind of the array in the file.  if the reading fail, return fform=0
!!
!! NOTES
!! The file is supposed to be open already
!!
!! PARENTS
!!      elphon,initaim,inpgkk,m_bse_io,m_cut3d,m_dvdb,m_hdr,m_io_screening
!!      m_ioarr,macroave,mrggkk,rchkgsheader,read_gkk
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_fort_read(Hdr,unit,fform,rewind)

!Arguments ------------------------------------
 integer,intent(out) :: fform
 integer,intent(in) :: unit
 logical,optional,intent(in) :: rewind
 type(hdr_type),intent(out) :: hdr

!Local variables-------------------------------
!integer :: ierr
 integer :: ipsp
 character(len=500) :: msg,errmsg
 real(dp),allocatable :: occ3d(:,:,:)

!*************************************************************************

 !@hdr_type
 DBG_ENTER("COLL")

 if (present(rewind)) then
   if (rewind) rewind(unit, err=10, iomsg=errmsg)
 end if

 ! Reading the first record of the file ------------------------------------
 ! fform is not a record of hdr_type
 ABI_CHECK(read_first_record(unit, hdr%codvsn, hdr%headform, fform, errmsg) == 0, errmsg)

 if (hdr%headform < 80) then
   write(msg,'(3a,i0,4a)') &
     "ABINIT version: ",trim(abinit_version)," cannot read old files with headform: ",hdr%headform,ch10,&
     "produced by previous versions. Use an old ABINIT version to read this file or ",ch10,&
     "regenerate your files with version >= 8.0."
   MSG_ERROR(msg)
 end if

 call check_fform(fform)

!Reading the second record of the file ------------------------------------
 read(unit, err=10, iomsg=errmsg) &
   hdr%bantot, hdr%date, hdr%intxc, hdr%ixc, hdr%natom, hdr%ngfft(1:3),&
   hdr%nkpt, hdr%nspden, hdr%nspinor, hdr%nsppol, hdr%nsym, hdr%npsp, hdr%ntypat, hdr%occopt, hdr%pertcase,&
   hdr%usepaw, hdr%ecut, hdr%ecutdg, hdr%ecutsm, hdr%ecut_eff, hdr%qptn(1:3), hdr%rprimd,&
   hdr%stmbias, hdr%tphysel, hdr%tsmear, hdr%usewvl, hdr%nshiftk_orig, hdr%nshiftk, hdr%mband

 !Allocate all parts of hdr that need to be --------------------------------
 call hdr_malloc(hdr, hdr%bantot, hdr%nkpt, hdr%nsppol, hdr%npsp, hdr%natom, hdr%ntypat,&
                 hdr%nsym, hdr%nshiftk_orig, hdr%nshiftk)

 if (hdr%usepaw==1)  then
   ABI_MALLOC(hdr%pawrhoij,(hdr%natom))
 end if

! Reading the third record of the file ------------------------------------

! Take into account future migration to occ(:,:,:) in the Format
! read 3d matrix with stride and transfer to (stupid) 1d hdr%occ in packed form.
 ABI_MALLOC(occ3d, (hdr%mband,hdr%nkpt,hdr%nsppol))

 read(unit, err=10, iomsg=errmsg) &
   hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:), &
   hdr%so_psp(:), hdr%symafm(:), hdr%symrel(:,:,:), &
   hdr%typat(:), hdr%kptns(:,:), occ3d, &
   hdr%tnons(:,:), hdr%znucltypat(:), hdr%wtk(:)
 ABI_CHECK(hdr%mband == maxval(hdr%nband), "mband != max(hdr%nband). Are you reading an Abinit8 file with Abinit9?")

 call hdr_set_occ(hdr, occ3d)
 ABI_FREE(occ3d)

! Reading the final record of the header  ---------------------------------
 read(unit, err=10, iomsg=errmsg) hdr%residm, hdr%xred(:,:), hdr%etot, hdr%fermie, hdr%amu(:)

 read(unit, err=10, iomsg=errmsg)&
    hdr%kptopt,hdr%pawcpxocc,hdr%nelect,hdr%charge,hdr%icoulomb,&
    hdr%kptrlatt,hdr%kptrlatt_orig, hdr%shiftk_orig,hdr%shiftk

! Reading the records with psp information ---------------------------------
 do ipsp=1,hdr%npsp
   read(unit, err=10, iomsg=errmsg) &
     hdr%title(ipsp), hdr%znuclpsp(ipsp), hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
     hdr%pspcod(ipsp), hdr%pspxc(ipsp), hdr%lmn_size(ipsp), hdr%md5_pseudos(ipsp)
 end do

 if (hdr%usepaw==1) then ! Reading the Rhoij tab if the PAW method was used.
   call pawrhoij_io(hdr%pawrhoij,unit,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,hdr%headform,"Read")
 end if

 DBG_EXIT("COLL")
 return

 ! Handle IO-error: write warning and let the caller handle the exception.
10 fform=0
 MSG_WARNING(errmsg)

end subroutine hdr_fort_read
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_ncread
!! NAME
!! hdr_ncread
!!
!! FUNCTION
!! This subroutine deals with the reading of the hdr_type structured variables
!! It handles variables according to the ETSF format, whenever
!! possible and uses new variables when not available in the ETSF format.
!! Note that, when reading, different records of hdr are allocated here,
!! Records of hdr should be deallocated
!! correctly by a call to hdr_free when hdr is not used anymore.
!!
!! INPUTS
!!  ncid=the unit of the open NetCDF file.
!!
!! OUTPUT
!!  fform=kind of the array in the file. if the reading fails, return fform=0
!!
!! PARENTS
!!      initaim,inwffil,m_ddk,m_dvdb,m_hdr,m_io_screening,m_ioarr,macroave
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_ncread(Hdr, ncid, fform)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 integer,intent(out) :: fform
 type(hdr_type),target,intent(out) :: hdr

#ifdef HAVE_NETCDF
!Local variables-------------------------------
!scalars
 integer :: nresolution, itypat, ii
 character(len=500) :: msg
!arrays
 integer,allocatable :: nband2d(:,:)
 real(dp),allocatable :: occ3d(:,:,:)

! *************************************************************************

 !@hdr_type
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_get_var(ncid, vid("fform"), fform))
 NCF_CHECK(nf90_get_var(ncid, vid("headform"), hdr%headform))

 if (hdr%headform < 80) then
   write(msg,'(3a,i0,4a)')&
     "ABINIT version: ",trim(abinit_version)," cannot read old files with headform: ",hdr%headform,ch10,&
     "produced by previous versions. Use an old ABINIT version to read this file or ",ch10,&
     "regenerate your files with version >= 8.0."
   MSG_ERROR(msg)
 end if

 call check_fform(fform)

 ! First, we read the declaration of code, fform ...
 ! pad the returned string with " " instead of "\0"
 !
 ! Support for pre-v9 (length of codvsn was changed from 6 to 8)
 NCF_CHECK(nctk_get_dim(ncid, "codvsnlen", ii))
 NCF_CHECK(nf90_get_var(ncid, vid("codvsn"), hdr%codvsn(1:ii)))
 call replace_ch0(hdr%codvsn)

 ! Get ETSF dimensions
 NCF_CHECK(nctk_get_dim(ncid, "number_of_atoms", hdr%natom))
 NCF_CHECK(nctk_get_dim(ncid, "number_of_kpoints", hdr%nkpt))
 NCF_CHECK(nctk_get_dim(ncid, "number_of_components", hdr%nspden))
 NCF_CHECK(nctk_get_dim(ncid, "number_of_spinor_components", hdr%nspinor))
 NCF_CHECK(nctk_get_dim(ncid, "number_of_spins", hdr%nsppol))
 NCF_CHECK(nctk_get_dim(ncid, "number_of_symmetry_operations", hdr%nsym))
 NCF_CHECK(nctk_get_dim(ncid, "number_of_atom_species", hdr%ntypat))
 NCF_CHECK(nctk_get_dim(ncid, "number_of_grid_points_vector1", hdr%ngfft(1)))
 NCF_CHECK(nctk_get_dim(ncid, "number_of_grid_points_vector2", hdr%ngfft(2)))
 NCF_CHECK(nctk_get_dim(ncid, "number_of_grid_points_vector3", hdr%ngfft(3)))
 NCF_CHECK(nctk_get_dim(ncid, "max_number_of_states", hdr%mband))
! bantot is used to dimension %occ in hdr_malloc. Note that hdr%bantot != sum(nband) because states
! are packed in hdr%occ and therefore bantot <= hdr%mband * hdr%nkpt * hdr%nsppol
 NCF_CHECK(nctk_get_dim(ncid, "bantot", hdr%bantot))

 ! Read other dimensions, not handled by ETSF format.
 NCF_CHECK(nctk_get_dim(ncid, "npsp", hdr%npsp))
 NCF_CHECK(nctk_get_dim(ncid, "nshiftk_orig", hdr%nshiftk_orig))
 NCF_CHECK(nctk_get_dim(ncid, "nshiftk", hdr%nshiftk))

 ! Read other important scalar variables.
 NCF_CHECK(nf90_get_var(ncid, vid("usepaw"), hdr%usepaw))
 NCF_CHECK(nf90_get_var(ncid, vid("usewvl"), hdr%usewvl))

 nresolution=0
 if (hdr%usewvl == 1) then
   ! This value must be 2...
   NCF_CHECK(nctk_get_dim(ncid, "number_of_wavelet_resolutions", nresolution))
   ! We set the right ngfft, adding the padding space for wavelets.
   hdr%ngfft = hdr%ngfft + 31
 end if

 ! Allocate all parts of hdr that need to be
 call hdr_malloc(hdr, hdr%bantot, hdr%nkpt, hdr%nsppol, hdr%npsp, hdr%natom, hdr%ntypat,&
                 hdr%nsym, hdr%nshiftk_orig, hdr%nshiftk)

 ABI_MALLOC(nband2d, (hdr%nkpt, hdr%nsppol))
 NCF_CHECK(nf90_get_var(ncid, vid("number_of_states"), nband2d))
 hdr%nband(:) = reshape(nband2d, [hdr%nkpt*hdr%nsppol])
 ABI_FREE(nband2d)
 ABI_CHECK(hdr%mband == maxval(hdr%nband), "mband != maxval(hdr%nband)")

 if (hdr%usepaw==1) then
   ABI_MALLOC(hdr%pawrhoij,(hdr%natom))
 end if

!We get then all variables included in ETSF
 if (hdr%usewvl==0) then
   NCF_CHECK(nf90_get_var(ncid, vid("kinetic_energy_cutoff"), hdr%ecut))
   NCF_CHECK(nf90_get_var(ncid, vid("number_of_coefficients"), hdr%npwarr))
 else
   NCF_CHECK(nf90_get_var(ncid, vid("number_of_wavelets"), hdr%nwvlarr))
 end if

! read 3d matrix with stride and transfer to (stupid) 1d hdr%occ in packed form.
 ABI_CALLOC(occ3d, (hdr%mband, hdr%nkpt, hdr%nsppol))
 NCF_CHECK(nf90_get_var(ncid, vid("occupations"), occ3d))
 call hdr_set_occ(hdr, occ3d)
 ABI_FREE(occ3d)

 NCF_CHECK(nf90_get_var(ncid, vid("fermi_energy"), hdr%fermie))
 NCF_CHECK(nf90_get_var(ncid, vid("primitive_vectors"), hdr%rprimd))
 NCF_CHECK(nf90_get_var(ncid, vid("reduced_symmetry_matrices"), hdr%symrel))
 NCF_CHECK(nf90_get_var(ncid, vid("atom_species"), hdr%typat))
 NCF_CHECK(nf90_get_var(ncid, vid("reduced_symmetry_translations"), hdr%tnons))
 NCF_CHECK(nf90_get_var(ncid, vid("reduced_atom_positions"), hdr%xred))
 NCF_CHECK(nf90_get_var(ncid, vid("atomic_numbers"), hdr%znucltypat))
 NCF_CHECK(nf90_get_var(ncid, vid("reduced_coordinates_of_kpoints"), hdr%kptns))
 NCF_CHECK(nf90_get_var(ncid, vid("kpoint_weights"), hdr%wtk))
 NCF_CHECK(nf90_get_var(ncid, vid("date"), hdr%date))
 NCF_CHECK(nf90_get_var(ncid, vid("ecut_eff"), hdr%ecut_eff))
 NCF_CHECK(nf90_get_var(ncid, vid("ecutsm"), hdr%ecutsm))
 NCF_CHECK(nf90_get_var(ncid, vid("etot"), hdr%etot))
 NCF_CHECK(nf90_get_var(ncid, vid("intxc"), hdr%intxc))
 NCF_CHECK(nf90_get_var(ncid, vid("ixc"), hdr%ixc))
 NCF_CHECK(nf90_get_var(ncid, vid("occopt"), hdr%occopt))
 NCF_CHECK(nf90_get_var(ncid, vid("pertcase"), hdr%pertcase))
 NCF_CHECK(nf90_get_var(ncid, vid("qptn"), hdr%qptn))
 NCF_CHECK(nf90_get_var(ncid, vid("residm"), hdr%residm))
 NCF_CHECK(nf90_get_var(ncid, vid("stmbias"), hdr%stmbias))
 NCF_CHECK(nf90_get_var(ncid, vid("tphysel"), hdr%tphysel))
 NCF_CHECK(nf90_get_var(ncid, vid("tsmear"), hdr%tsmear))
 NCF_CHECK(nf90_get_var(ncid, vid("ecutdg"), hdr%ecutdg))

! Multidimensional variables. Be careful with zionpsp if alchemical mixing!
 NCF_CHECK(nf90_get_var(ncid, vid("istwfk"), hdr%istwfk))
 NCF_CHECK(nf90_get_var(ncid, vid("pspcod"), hdr%pspcod))
 NCF_CHECK(nf90_get_var(ncid, vid("pspdat"), hdr%pspdat))
 NCF_CHECK(nf90_get_var(ncid, vid("pspso"), hdr%pspso))
 NCF_CHECK(nf90_get_var(ncid, vid("pspxc"), hdr%pspxc))
 NCF_CHECK(nf90_get_var(ncid, vid("so_psp"), hdr%so_psp))
 NCF_CHECK(nf90_get_var(ncid, vid("symafm"), hdr%symafm))
 NCF_CHECK(nf90_get_var(ncid, vid("zionpsp"), hdr%zionpsp))
 NCF_CHECK(nf90_get_var(ncid, vid("znuclpsp"), hdr%znuclpsp))
 NCF_CHECK(nf90_get_var(ncid, vid("kptopt"), hdr%kptopt))
 NCF_CHECK(nf90_get_var(ncid, vid("pawcpxocc"), hdr%pawcpxocc))
 NCF_CHECK(nf90_get_var(ncid, vid("nelect"), hdr%nelect))
 NCF_CHECK(nf90_get_var(ncid, vid("charge"), hdr%charge))
 NCF_CHECK(nf90_get_var(ncid, vid("kptrlatt_orig"), hdr%kptrlatt_orig))
 NCF_CHECK(nf90_get_var(ncid, vid("kptrlatt"), hdr%kptrlatt))
 NCF_CHECK(nf90_get_var(ncid, vid("shiftk_orig"), hdr%shiftk_orig))
 NCF_CHECK(nf90_get_var(ncid, vid("shiftk"), hdr%shiftk))
 NCF_CHECK(nf90_get_var(ncid, vid("md5_pseudos"), hdr%md5_pseudos))
 NCF_CHECK(nf90_get_var(ncid, vid("amu"), hdr%amu))
 NCF_CHECK(nf90_get_var(ncid, vid("icoulomb"), hdr%icoulomb))
 NCF_CHECK(nf90_get_var(ncid, vid("title"), hdr%title))

 ! Pad the returned string with " " instead of "\0"
 do itypat=1,size(hdr%title)
   call replace_ch0(hdr%title(itypat))
 end do

 NCF_CHECK(nf90_get_var(ncid, vid("lmn_size"), hdr%lmn_size))
 if (hdr%usepaw==1) then
   call pawrhoij_io(hdr%pawrhoij,ncid,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,&
      hdr%headform,"Read",form="netcdf")
 end if

#else
 MSG_ERROR("netcdf support not activated")
#endif

contains
 integer function vid(vname)

   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine hdr_ncread
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_fort_write
!! NAME
!! hdr_fort_write
!!
!! FUNCTION
!!  Writes the header and fform to unformatted file
!!
!! INPUTS
!!  Hdr<hdr_type>=The header of the file.
!!  fform=kind of the array in the file
!!  unit=unit number of the unformatted file
!!  [rewind]=True to rewind the file. Default: False
!!
!! OUTPUT
!!  ierr=Exit status
!!
!! NOTES
!! The file is supposed to be open already
!!
!! PARENTS
!!      m_bse_io,m_dvdb,m_hdr,m_io_kss,m_io_screening,m_ioarr,mrggkk,outgkk
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_fort_write(Hdr,unit,fform,ierr,rewind)

!Arguments ------------------------------------
 integer,intent(out) :: ierr
 integer,intent(in) :: unit,fform
 logical,optional,intent(in) :: rewind
 class(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------
 integer :: headform,ipsp
 character(len=500) :: errmsg
 real(dp),allocatable :: occ3d(:,:,:)

!*************************************************************************

 ! TODO: Change intent to in. Change pawrhoij_io first!
 !@hdr_type
 ierr = 0

 if (present(rewind)) then
   if (rewind) rewind(unit, err=10, iomsg=errmsg)
 end if

 call check_fform(fform)

!Writing always use last format version
 headform = HDR_LATEST_HEADFORM
 write(unit, err=10, iomsg=errmsg) hdr%codvsn, headform, fform

 write(unit, err=10, iomsg=errmsg) &
   hdr%bantot, hdr%date, hdr%intxc, hdr%ixc, hdr%natom, hdr%ngfft(1:3), &
   hdr%nkpt, hdr%nspden, hdr%nspinor, hdr%nsppol, hdr%nsym, hdr%npsp, hdr%ntypat, hdr%occopt, hdr%pertcase,&
   hdr%usepaw, hdr%ecut, hdr%ecutdg, hdr%ecutsm, hdr%ecut_eff, hdr%qptn, hdr%rprimd, &
   hdr%stmbias, hdr%tphysel, hdr%tsmear, hdr%usewvl, hdr%nshiftk_orig, hdr%nshiftk, hdr%mband
 ABI_CHECK(hdr%mband == maxval(hdr%nband), "mband != maxval(hdr%nband)")

 ABI_MALLOC(occ3d, (hdr%mband,hdr%nkpt,hdr%nsppol))
 call hdr_get_occ3d(hdr, occ3d)
 write(unit,err=10, iomsg=errmsg) hdr%istwfk(:), hdr%nband(:), hdr%npwarr(:),&
   hdr%so_psp(:), hdr%symafm(:), hdr%symrel(:,:,:), hdr%typat(:), hdr%kptns(:,:), occ3d, &
   hdr%tnons(:,:), hdr%znucltypat(:), hdr%wtk(:)
 ABI_FREE(occ3d)

 write(unit,err=10, iomsg=errmsg) hdr%residm, hdr%xred(:,:), hdr%etot, hdr%fermie, hdr%amu(:)
 write(unit,err=10, iomsg=errmsg) &
   hdr%kptopt, hdr%pawcpxocc, hdr%nelect, hdr%charge, hdr%icoulomb,&
   hdr%kptrlatt,hdr%kptrlatt_orig, hdr%shiftk_orig(:,1:hdr%nshiftk_orig),hdr%shiftk(:,1:hdr%nshiftk)

 ! Write the records with psp information ---------------------------------
 do ipsp=1,hdr%npsp
   write(unit, err=10, iomsg=errmsg) &
     hdr%title(ipsp), hdr%znuclpsp(ipsp), hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
     hdr%pspcod(ipsp), hdr%pspxc(ipsp), hdr%lmn_size(ipsp), hdr%md5_pseudos(ipsp)
 end do

 if (hdr%usepaw==1) then
   call pawrhoij_io(hdr%pawrhoij,unit,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,headform,"Write")
 end if

 return

 ! Handle IO-error: write warning and let the caller handle the exception.
10 ierr=1
 MSG_WARNING(errmsg)

end subroutine hdr_fort_write
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_backspace
!! NAME
!! hdr_backspace
!!
!! FUNCTION
!!  Backspace the header. Return exit status and error message
!!  The file is supposed to be open already
!!
!! INPUTS
!!  Hdr<hdr_type>=The header of the file.
!!  unit=unit number of the unformatted file
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function hdr_backspace(hdr, unit, msg) result(ierr)

!Arguments ------------------------------------
 class(hdr_type),intent(in) :: hdr
 integer,intent(in) :: unit
 character(len=*),intent(out) :: msg

!Local variables-------------------------------
 integer :: irec

!*************************************************************************

 ierr = 0
 do irec=1,5 + hdr%npsp
   backspace(unit=unit, err=10, iomsg=msg)
 end do

 if (hdr%usepaw == 1) then
   do irec=1,2
     backspace(unit=unit, err=10, iomsg=msg)
   end do
 end if

 return

 ! Handle IO-error
10 ierr = 1

end function hdr_backspace
!!***

!!****f* m_hdr/hdr_ncwrite
!! NAME
!! hdr_ncwrite
!!
!! FUNCTION
!! This subroutine deals with the output of the hdr_type structured variables in ETSF+NETCDF fornat.
!! It handles variables according to the ETSF format, whenever possible and uses new variables
!! when not available in the ETSF format.
!!
!! INPUTS
!!  fform=kind of the array in the file
!!  ncid=the unit of the open NetCDF file.
!!  [nc_define]=Optional flag. If True, the basic dimensions required by the ETSF specification
!!    are written. Default: False.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      ioarr,m_hdr,m_wfk
!!
!! CHILDREN
!!
!! SOURCE

integer function hdr_ncwrite(hdr, ncid, fform, nc_define) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid,fform
 logical,optional,intent(in) :: nc_define
 class(hdr_type),target,intent(in) :: hdr

#ifdef HAVE_NETCDF
!Local variables-------------------------------
!scalars
 logical :: my_define
 character(len=etsfio_charlen) :: basis_set,k_dependent,symmorphic
 !character(len=500) :: msg
!arrays
 integer,allocatable :: arr2d(:,:)
 real(dp),allocatable :: arr3d(:,:,:)
 type(pawrhoij_type),pointer :: rhoij_ptr(:)

! *************************************************************************

 call check_fform(fform)

 !@hdr_type
 my_define = .False.; if (present(nc_define)) my_define = nc_define
 ncerr = nf90_noerr

 k_dependent = "no"; if (any(hdr%nband(1) /= hdr%nband)) k_dependent = "yes"
 symmorphic = "no"; if (all(abs(hdr%tnons) < tol6)) symmorphic = "yes"

 if (my_define) then
   !call wrtout(std_out, "hdr_ncwrite: defining variables")
   NCF_CHECK(nctk_def_basedims(ncid, defmode=.True.))

   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("max_number_of_states", hdr%mband), &
     nctkdim_t("number_of_atoms", hdr%natom), &
     nctkdim_t("number_of_atom_species", hdr%ntypat), &
     nctkdim_t("number_of_components", hdr%nspden), &
     nctkdim_t("number_of_kpoints", hdr%nkpt), &
     nctkdim_t("number_of_spinor_components", hdr%nspinor), &
     nctkdim_t("number_of_spins", hdr%nsppol), &
     nctkdim_t("number_of_symmetry_operations", hdr%nsym) &
   ])
     !nctkdim_t("nshiftk_orig", ebands%nshiftk_orig), &
     !nctkdim_t("nshiftk", ebands%nshiftk)], &
   NCF_CHECK(ncerr)

   ! Define part of geometry section contained in the header.
   ncerr = nctk_def_arrays(ncid, [ &
    ! Atomic structure and symmetry operations
    nctkarr_t("primitive_vectors", "dp", "number_of_cartesian_directions, number_of_vectors"), &
    nctkarr_t("reduced_symmetry_matrices", "int", &
      "number_of_reduced_dimensions, number_of_reduced_dimensions, number_of_symmetry_operations"), &
    nctkarr_t("reduced_symmetry_translations", "dp", "number_of_reduced_dimensions, number_of_symmetry_operations"), &
    nctkarr_t("atom_species", "int", "number_of_atoms"), &
    nctkarr_t("reduced_atom_positions", "dp", "number_of_reduced_dimensions, number_of_atoms"), &
    nctkarr_t("atomic_numbers", "dp", "number_of_atom_species") &
    !nctkarr_t("atom_species_names", "char", "character_string_length, number_of_atom_species"), &
    !nctkarr_t("chemical_symbols", "char", "symbol_length, number_of_atom_species"), &
    ! Atomic information.
    !nctkarr_t("valence_charges", "dp", "number_of_atom_species"), &  ! NB: This variable is not written if alchemical
    !nctkarr_t("pseudopotential_types", "char", "character_string_length, number_of_atom_species") &
   ])
   NCF_CHECK(ncerr)

   ! Some variables require the "symmorphic" attribute.
   NCF_CHECK(nf90_put_att(ncid, vid("reduced_symmetry_matrices"), "symmorphic", symmorphic))
   NCF_CHECK(nf90_put_att(ncid, vid("reduced_symmetry_translations"), "symmorphic", symmorphic))

   ! At this point we have an ETSF-compliant file. Add additional data for internal use in abinit.
   ! TODO add spinat.
   ncerr = nctk_def_arrays(ncid, nctkarr_t('symafm', "int", "number_of_symmetry_operations"))
   NCF_CHECK(ncerr)

   ! Define k-points. Note: monkhorst_pack_folding is replaced by kptrlatt and shiftk
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t("reduced_coordinates_of_kpoints", "dp", "number_of_reduced_dimensions, number_of_kpoints"), &
     nctkarr_t("kpoint_weights", "dp", "number_of_kpoints") &
     !nctkarr_t("monkhorst_pack_folding", "int", "number_of_vectors") &
   ])
   NCF_CHECK(ncerr)

   ! Define states section.
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("number_of_states", "int", "number_of_kpoints, number_of_spins"), &
     nctkarr_t("eigenvalues", "dp", "max_number_of_states, number_of_kpoints, number_of_spins"), &
     nctkarr_t("occupations", "dp", "max_number_of_states, number_of_kpoints, number_of_spins"), &
     nctkarr_t("smearing_scheme", "char", "character_string_length")  &
   ])
   NCF_CHECK(ncerr)

   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "number_of_electrons"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "fermi_energy", "smearing_width"])
   NCF_CHECK(ncerr)
   NCF_CHECK(nctk_set_atomic_units(ncid, "smearing_width"))

   ! Some variables require the specifications of units.
   NCF_CHECK(nctk_set_atomic_units(ncid, "eigenvalues"))
   NCF_CHECK(nctk_set_atomic_units(ncid, "fermi_energy"))
   NCF_CHECK(nf90_put_att(ncid, vid("number_of_states"), "k_dependent", k_dependent))

   ! Define dimensions.
   ncerr = nctk_def_dims(ncid, [&
     nctkdim_t("npsp", hdr%npsp), nctkdim_t("codvsnlen", 8), nctkdim_t("psptitlen", 132)&
   ])
   NCF_CHECK(ncerr)

   if (hdr%usewvl==1) then ! Add the BigDFT private dimensions.
     ncerr = nctk_def_dims(ncid, nctkdim_t("number_of_wavelet_resolutions", 2))
     NCF_CHECK(ncerr)
   end if

   ! Define scalars.
   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
     "date", "ixc", "intxc", "occopt", "pertcase", "headform", "fform", "usepaw", "usewvl"])
   NCF_CHECK(ncerr)

   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
    "ecut_eff", "ecutdg", "ecutsm", "etot", "residm", "stmbias", "tphysel", "tsmear"])
   NCF_CHECK(ncerr)

   ! Multi-dimensional variables.
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t("istwfk", "i", "number_of_kpoints"),&
     nctkarr_t("codvsn", "c", "codvsnlen"),&
     nctkarr_t("pspcod", "i", "npsp"),&
     nctkarr_t("pspdat", "i", "npsp"),&
     nctkarr_t("pspso", "i", "npsp"),&
     nctkarr_t("pspxc", "i", "npsp"),&
     nctkarr_t("qptn", "dp", "number_of_reduced_dimensions"),&
     nctkarr_t("so_psp", "i", "npsp"),&
     nctkarr_t("symafm", "i", "number_of_symmetry_operations"),&
     nctkarr_t("title", "c", "psptitlen, npsp"),&
     nctkarr_t("zionpsp", "dp", "npsp"),&
     nctkarr_t("znuclpsp", "dp", "npsp"),&
     nctkarr_t("lmn_size", "i", "npsp")])
   NCF_CHECK(ncerr)

   ! Add the BigDFT private variables.
   if (hdr%usewvl == 1) then
     ncerr = nctk_def_arrays(ncid, nctkarr_t("number_of_wavelets", "i", "number_of_wavelet_resolutions"))
     NCF_CHECK(ncerr)
   end if

   NCF_CHECK(nctk_def_arrays(ncid, nctkarr_t("basis_set", "char", "character_string_length")))
   if (hdr%usewvl == 0) then
     NCF_CHECK(nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "kinetic_energy_cutoff"]))
     NCF_CHECK(nctk_set_atomic_units(ncid, "kinetic_energy_cutoff"))
     NCF_CHECK(nctk_def_arrays(ncid, nctkarr_t("number_of_coefficients", "int", "number_of_kpoints")))
   end if

   NCF_CHECK(nf90_put_att(ncid, vid("number_of_states"), "k_dependent", k_dependent))

   if (hdr%usewvl == 0) then
     ! Note that here we always use the coarse FFT mesh even if usepaw == 1
     ncerr = nctk_def_dims(ncid, [&
       nctkdim_t("number_of_grid_points_vector1", hdr%ngfft(1)),&
       nctkdim_t("number_of_grid_points_vector2", hdr%ngfft(2)),&
       nctkdim_t("number_of_grid_points_vector3", hdr%ngfft(3))], defmode=.True.)
     NCF_CHECK(ncerr)
   else
     MSG_WARNING("Don't know how to define grid_points for wavelets!")
   end if

   !write(std_out,*)"hdr%nshiftk_orig,hdr%nshiftk",hdr%nshiftk_orig,hdr%nshiftk
   ncerr = nctk_def_dims(ncid, [&
     nctkdim_t("nshiftk_orig", hdr%nshiftk_orig),&
     nctkdim_t("nshiftk", hdr%nshiftk), &
     nctkdim_t("bantot", hdr%bantot), &
     nctkdim_t("md5_slen", md5_slen)], defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "kptopt", "pawcpxocc", "icoulomb"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "nelect", "charge"])
   NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t("kptrlatt_orig", "i", "number_of_reduced_dimensions, number_of_reduced_dimensions"),&
     nctkarr_t("kptrlatt", "i", "number_of_reduced_dimensions, number_of_reduced_dimensions"),&
     nctkarr_t("shiftk_orig", "dp", "number_of_reduced_dimensions, nshiftk_orig"),&
     nctkarr_t("shiftk", "dp", "number_of_reduced_dimensions, nshiftk"), &
     nctkarr_t("amu", "dp", "number_of_atom_species"), &
     nctkarr_t("md5_pseudos", "ch", "md5_slen, npsp") ])
   NCF_CHECK(ncerr)

   !call wrtout(std_out, "hdr_ncwrite completed define mode")
 end if ! my_define

 ! Switch to write mode.
 NCF_CHECK(nctk_set_datamode(ncid))

 ! Associate and write values to ETSF groups.
 if (hdr%usewvl == 0) then
   ! Plane wave case.
   basis_set = "plane_waves"
   NCF_CHECK(nf90_put_var(ncid, vid("basis_set"), basis_set))
   NCF_CHECK(nf90_put_var(ncid, vid("kinetic_energy_cutoff"), hdr%ecut))
   NCF_CHECK(nf90_put_var(ncid, vid("number_of_coefficients"), hdr%npwarr))
 else
   ! Wavelet case.
   basis_set = "daubechies_wavelets"
   NCF_CHECK(nf90_put_var(ncid, vid("basis_set"), basis_set))
   ! Required variable than should enter the standard.
   NCF_CHECK(nf90_put_var(ncid, vid("number_of_wavelets"), hdr%nwvlarr))
 end if

 ! Write electrons
 NCF_CHECK(nf90_put_var(ncid, vid("fermi_energy"), hdr%fermie))
 NCF_CHECK(nf90_put_var(ncid, vid("smearing_width"), hdr%tsmear))
 NCF_CHECK(nf90_put_var(ncid, vid("smearing_scheme"), nctk_string_from_occopt(hdr%occopt)))

 ! transfer data from (stupid) 1d hdr%nband and hdr%occ in packed form to 2d - 3d matrix with stride
 ! native support for array and array syntax is one of the reasons why we still use Fortran
 ! and we program like in C but without the power of C!o

! Also, strange problem with Petrus + Nag5: had to explicitly specify nf90_put_var,
! with explicit definition of start, count and stride .
! Direct calls to NCF_CHECK, see below, were working for selected tests, but not all tests
 ABI_MALLOC(arr2d, (hdr%nkpt, hdr%nsppol))
 arr2d(:,:) = reshape(hdr%nband, [hdr%nkpt, hdr%nsppol])
 ncerr = nf90_put_var(ncid, vid("number_of_states"), arr2d, start=[1,1], count=[hdr%nkpt,hdr%nsppol], stride=[1,1])
 NCF_CHECK(ncerr)
 ABI_FREE(arr2d)

 ABI_MALLOC(arr3d, (hdr%mband, hdr%nkpt, hdr%nsppol))
 call hdr_get_occ3d(hdr, arr3d)
 NCF_CHECK(nf90_put_var(ncid, vid("occupations"), arr3d))
 ABI_FREE(arr3d)

 ! Write geometry
 NCF_CHECK(nf90_put_var(ncid, vid("primitive_vectors"), hdr%rprimd))
 NCF_CHECK(nf90_put_var(ncid, vid("reduced_symmetry_matrices"), hdr%symrel))
 NCF_CHECK(nf90_put_var(ncid, vid("atom_species"), hdr%typat))
 NCF_CHECK(nf90_put_var(ncid, vid("reduced_symmetry_translations"), hdr%tnons))
 NCF_CHECK(nf90_put_var(ncid, vid("reduced_atom_positions"), hdr%xred))
 NCF_CHECK(nf90_put_var(ncid, vid("atomic_numbers"), hdr%znucltypat))

 ! Write k-points.
 NCF_CHECK(nf90_put_var(ncid, vid("reduced_coordinates_of_kpoints"), hdr%kptns))
 NCF_CHECK(nf90_put_var(ncid, vid("kpoint_weights"), hdr%wtk))

 ! Write non-ETSF variables.
 NCF_CHECK(nf90_put_var(ncid, vid("codvsn"), hdr%codvsn))
 NCF_CHECK(nf90_put_var(ncid, vid("title"), hdr%title))

 ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
&  "date", "ixc", "intxc", "occopt", "pertcase", "headform", "fform", "usepaw", "icoulomb"],&
&  [hdr%date, hdr%ixc ,hdr%intxc ,hdr%occopt, hdr%pertcase, HDR_LATEST_HEADFORM, fform, hdr%usepaw, hdr%icoulomb])
 NCF_CHECK(ncerr)

 ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
&  "ecut_eff", "ecutdg", "ecutsm", "etot", "residm", "stmbias", "tphysel", "tsmear"],&
&  [hdr%ecut_eff, hdr%ecutdg, hdr%ecutsm, hdr%etot, hdr%residm, hdr%stmbias, hdr%tphysel, hdr%tsmear])
 NCF_CHECK(ncerr)

!Array variables.

! FIXME Be careful with zionpsp if alchemical mixing!
 NCF_CHECK(nf90_put_var(ncid, vid("istwfk"), hdr%istwfk))
 NCF_CHECK(nf90_put_var(ncid, vid("pspcod"), hdr%pspcod))
 NCF_CHECK(nf90_put_var(ncid, vid("pspdat"), hdr%pspdat))
 NCF_CHECK(nf90_put_var(ncid, vid("pspso"), hdr%pspso))
 NCF_CHECK(nf90_put_var(ncid, vid("pspxc"), hdr%pspxc))
 NCF_CHECK(nf90_put_var(ncid, vid("qptn"), hdr%qptn))
 NCF_CHECK(nf90_put_var(ncid, vid("so_psp"), hdr%so_psp))
 NCF_CHECK(nf90_put_var(ncid, vid("symafm"), hdr%symafm))
 NCF_CHECK(nf90_put_var(ncid, vid("znuclpsp"), hdr%znuclpsp))
 NCF_CHECK(nf90_put_var(ncid, vid("zionpsp"), hdr%zionpsp))
 NCF_CHECK(nf90_put_var(ncid, vid("lmn_size"), hdr%lmn_size))
 NCF_CHECK(nf90_put_var(ncid, vid("usewvl"), hdr%usewvl))

 ! Write hdr%pawrhoij.
 if (hdr%usepaw == 1) then
   ! Dirty trick to bypass check on the intent, but the problem is in the intent(inout) of pawrhoij_io
   rhoij_ptr => hdr%pawrhoij
   call pawrhoij_io(rhoij_ptr,ncid,hdr%nsppol,hdr%nspinor,hdr%nspden,hdr%lmn_size,hdr%typat,&
                    HDR_LATEST_HEADFORM,"Write",form="netcdf")
 end if

 ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
   "kptopt", "pawcpxocc"],[hdr%kptopt, hdr%pawcpxocc])
 NCF_CHECK(ncerr)

 ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
   "nelect", "charge"],[hdr%nelect, hdr%charge])
 NCF_CHECK(ncerr)

 ! NB: In etsf_io the number of electrons is declared as integer.
 ! We use abinit nelect to store the value as real(dp).
 NCF_CHECK(nf90_put_var(ncid, vid("number_of_electrons"), nint(hdr%nelect)))
 NCF_CHECK(nf90_put_var(ncid, vid("kptrlatt_orig"), hdr%kptrlatt_orig))
 NCF_CHECK(nf90_put_var(ncid, vid("kptrlatt"), hdr%kptrlatt))
 NCF_CHECK(nf90_put_var(ncid, vid("shiftk_orig"), hdr%shiftk_orig))
 NCF_CHECK(nf90_put_var(ncid, vid("shiftk"), hdr%shiftk))
 NCF_CHECK(nf90_put_var(ncid, vid("md5_pseudos"), hdr%md5_pseudos))
 NCF_CHECK(nf90_put_var(ncid, vid("amu"), hdr%amu))

#else
 MSG_ERROR("netcdf support not activated")
#endif

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end function hdr_ncwrite
!!***

!!****f* m_hdr/hdr_set_occ
!! NAME
!! hdr_set_occ
!!
!! FUNCTION
!!  Set the occuations hdr%occ(:) from a 3d array with stride.
!!
!! PARENTS
!!      m_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_set_occ(hdr, occ3d)

!Arguments ------------------------------------
 class(hdr_type),intent(inout) :: hdr
 real(dp),intent(in) :: occ3d(hdr%mband,hdr%nkpt,hdr%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ii,band,ikpt,spin

!*************************************************************************

 ii = 0
 do spin=1,hdr%nsppol
   do ikpt=1,hdr%nkpt
     do band=1,hdr%nband(ikpt + (spin-1) * hdr%nkpt)
         ii = ii +1
         hdr%occ(ii) = occ3d(band,ikpt,spin)
     end do
   end do
 end do

end subroutine hdr_set_occ
!!***

!!****f* m_hdr/hdr_get_occ3d
!! NAME
!! hdr_get_occ3d
!!
!! FUNCTION
!!  Return occupations in a 3d array with stride.
!!
!! PARENTS
!!      m_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_get_occ3d(hdr, occ3d)

!Arguments ------------------------------------
 class(hdr_type),intent(in) :: hdr
 real(dp),intent(out) :: occ3d(hdr%mband,hdr%nkpt,hdr%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ii,band,ikpt,spin

!*************************************************************************

 ii = 0; occ3d = huge(one)
 do spin=1,hdr%nsppol
   do ikpt=1,hdr%nkpt
     do band=1,hdr%nband(ikpt + (spin-1) * hdr%nkpt)
         ii = ii +1
         occ3d(band,ikpt,spin) = hdr%occ(ii)
     end do
   end do
 end do

end subroutine hdr_get_occ3d
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_check
!! NAME
!! hdr_check
!!
!! FUNCTION
!! This subroutine compare the header structured variable (hdr)
!! from input data (mostly dtset and psps) with the one (hdr0) of
!! an input data file (e.g. wf, density, potential).
!! Various values are checked for agreement or near agreement in the
!! case of floating point numbers.  The program will exit or produce
!! warning messages when unexpected values are found.
!! A record of the comparison of the headers is written to stdout.
!!
!! Decisions have been taken about whether a restart is allowed.
!! In the self-consistent case, a restart will always be allowed, but
!! one has to distinguish between a direct restart and a restart with
!! translation of wavefunction.
!! In the non-self-consistent case, the conditions below
!! must be fulfilled to allow a restart.
!!
!! INPUTS
!!  fform=integer specification of data type (expected)
!!  fform0=integer specification of data type (from disk file)
!!  mode_paral: COLL or PERS, for all wrtout calls
!!  hdr <type(hdr_type)>=the header structured variable from dtset and psps
!!  hdr0<type(hdr_type)>=the header structured variable from the disk file
!!
!! OUTPUT
!!  restart=1 if direct restart, =2 if translation is needed, =0 if no
!!              restart is possible.
!!  restartpaw= deals with the additional information in the PAW method
!!              =1 if direct restart, =0 if no restart from spherical data is possible.
!!              also 0 if no restart is possible
!!
!! NOTES
!! In the current version of the user interface restarts are allowed from
!! wavefunction files for self-consistent runs and from densities for
!! non-self-consistent runs. The precise conditions under which we will
!! allow a restart in this release are as follows.
!!
!!           self-consistent case : direct restarts
!!           ======================================
!!
!! A direct restart will be allowed provided the following quantities in
!! old and new calculations are the same:
!!
!!   (A) the primitive vectors                             (tprim)
!!   (B) the plane-wave cutoff                             (tecut)
!!   (C) nkpt, kpt(3,nkpt), wtk(nkpt)                      (tkpt)
!!   (D) istwfk(nkpt), the format of wavefunctions         (twfk)
!!   (E) nspinor, the scalar or spinor wf characteristics  (tspinor)
!! For PAW calculations:
!!   (F) the use of PAW method                             (tpaw)
!!   (G) the number of lmn elements for the paw basis      (tlmn)
!!   (H) the energy cutoff for the double (fine) grid      (tdg)
!! For WVL calculations:
!!   (I) the number of wavelets differs                    (twvl)
!!   (J) the space-grid size differs                       (tgrid)
!!
!!            non-self-consistent restarts
!!            ============================
!!
!! A restart will be allowed provided the following quantities in
!! old and new calculation are the same
!!
!!   (A) the primitive vectors                            (tprim)
!!   (B) the number of atoms of each type                 (tatty)
!!   (C) xred(3,natom)                                    (txred)
!!   (D) pseudopotentials (not just pseudocharges)        (tpseu)
!!   (E) the plane-wave cutoff                            (tecut)
!!   (F) ngfft(1:3)                                       (tng)
!! For PAW calculations:
!!   (G) the use of PAW method                            (tpaw)
!!   (H) the number of lmn elements for the paw basis     (tlmn)
!!   (I) the energy cutoff for the double (fine) grid     (tdg)
!!
!! PARENTS
!!      inwffil,m_ddk,m_io_screening,m_ioarr,m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine hdr_check(fform, fform0, hdr, hdr0, mode_paral, restart, restartpaw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fform,fform0
 integer,intent(out) :: restart,restartpaw
 character(len=4),intent(in) :: mode_paral
 type(hdr_type),intent(in) :: hdr,hdr0

!Local variables-------------------------------
 character(len=500) :: bndfmt, occfmt, wtkfmt, zatfmt, typfmt
!scalars
 integer,parameter :: mwarning=5,nkpt_max=5
 integer :: bantot,bantot_eff,ii,ipsp,isppol,istart,istop,isym,itest,iwarning
 integer :: jj,mu,natom,nelm,nkpt,npsp,nsppol,nsym,ntypat,tatty,tband,tdg
 integer :: tecut,tgrid,tkpt,tlmn,tng,tpaw,tprim,tpsch,tpseu,tspinor,tsym,twfk
 integer :: twvl,txred,enough
 real(dp) :: rms
 logical :: tfform2,tfform52
 character(len=500) :: msg
 type(abifile_t) :: abifile,abifile0

! *************************************************************************

 !@hdr_type
 DBG_ENTER("COLL")

 ! We will adopt convention that if things agree between restart
 ! and current calculation then the tflag is 0. Begin by assuming
 ! that there is complete agreement between the files

 tatty = 0; tband = 0; tdg = 0 ; tecut = 0; tkpt = 0;
 tlmn = 0; tng = 0; tpaw = 0; tprim = 0; tpsch = 0; tpseu = 0;
 tspinor=0; tsym = 0; twfk = 0 ; txred = 0 ; twvl = 0 ; tgrid = 0

 ! Write out a header
 write(msg,'(a1,80a,2a1,10x,a,3a1,10x,a,27x,a,a1,10x,19a,27x,12a,a1)' )&
   ch10,('=',ii=1,80),ch10,ch10,&
   '- hdr_check: checking restart file header for consistency -',&
   (ch10,ii=1,3),'current calculation','restart file',ch10,('-',ii=1,19),('-',ii=1,12),ch10
 call wrtout(std_out,msg,mode_paral)

 ! Check validity of fform, and find filetype
 abifile = abifile_from_fform(fform)
 if (abifile%fform == 0) then
    MSG_ERROR(sjoin("Cannot find any abifile object associated to fform:", itoa(fform)))
 end if

 ! Check validity of fform0, and find filetype
 abifile0 = abifile_from_fform(fform0)
 if (abifile0%fform == 0) then
    MSG_ERROR(sjoin("Cannot find any abifile object associated to fform:", itoa(fform0)))
 end if

 write(msg,'(a,a17,3x,2a,a17)') &
  '  calculation expects a ',ljust(abifile%class, 17),'|','  input file contains a ',ljust(abifile0%class, 17)
 call wrtout(std_out,msg,mode_paral)

 write(msg,'(a,a,13x,a,a,a)')&
  '. ABINIT  code version ',hdr%codvsn,'|','  ABINIT  code version ',hdr0%codvsn
 call wrtout(std_out,msg,mode_paral)

 ! Check fform from input, not from header file
 if ( fform /= fform0) then
   write(msg,'(a,i0,a,i0,a)')'input fform=',fform,' differs from disk file fform=',fform0,'.'
   MSG_ERROR(msg)
 end if

 write(msg, '(a,i8,a,i8,a,i4,2x,a,a,i8,a,i8,a,i4)' ) &
  '. date ',hdr %date,' bantot ',hdr %bantot,' natom ',hdr %natom,'|',&
  '  date ',hdr0%date,' bantot ',hdr0%bantot,' natom ',hdr0%natom
 call wrtout(std_out,msg,mode_paral)

 write(msg, '(a,i8,a,i3,3(a,i4),2x,a,a,i8,a,i3,3(a,i4))' )&
  '  nkpt',hdr %nkpt,' nsym',hdr %nsym,' ngfft',hdr %ngfft(1),',',hdr %ngfft(2),',',hdr %ngfft(3),'|',&
  '  nkpt',hdr0%nkpt,' nsym',hdr0%nsym,' ngfft',hdr0%ngfft(1),',',hdr0%ngfft(2),',',hdr0%ngfft(3)
 call wrtout(std_out,msg,mode_paral)

 if (hdr%usewvl == 0) then
   ! Note that the header actually contains ecut_eff=ecut*dilatmx**2
   write(msg,'(a,i3,a,f12.7,12x,a,a,i3,a,f12.7)')&
    '  ntypat',hdr %ntypat,' ecut_eff',hdr %ecut_eff,'|',&
    '  ntypat',hdr0%ntypat,' ecut_eff',hdr0%ecut_eff
   call wrtout(std_out,msg,mode_paral)
 else
   write(msg,'(a,i3,a,f12.7,12x,a,a,i3,a,f12.7)')&
    '  ntypat',hdr %ntypat,' hgrid   ', 2. * hdr %rprimd(1,1) / (hdr %ngfft(1) - 31),'|',&
    '  ntypat',hdr0%ntypat,' hgrid   ', 2. * hdr0%rprimd(1,1) / (hdr0%ngfft(1) - 31)
   call wrtout(std_out,msg,mode_paral)
   ! Check hgrid and rprimd values.
   if (hdr0%rprimd(1,2) /= zero .or. hdr0%rprimd(1,3) /= zero .or. &
    hdr0%rprimd(2,1) /= zero .or. hdr0%rprimd(2,3) /= zero .or. &
    hdr0%rprimd(3,1) /= zero .or. hdr0%rprimd(3,2) /= zero) then
     MSG_ERROR('disk file rprimd is not parallelepipedic.')
   end if
   if (abs(hdr0%rprimd(1,1) / hdr0%ngfft(1) - hdr %rprimd(1,1) / hdr %ngfft(1)) > tol8) then
     write(msg,'(a,F7.4,a,F7.4)')&
      'input wvl_hgrid=', 2. * hdr%rprimd(1,1) / hdr%ngfft(1), &
      'not equal disk file wvl_hgrid=', 2. * hdr0%rprimd(1,1) / hdr0%ngfft(1)
     MSG_COMMENT(msg)
     tgrid = 1
   end if
 end if

 write(msg, '(a,i3,33x,a,a,i3)' )'  usepaw',hdr %usepaw,'|','  usepaw',hdr0%usepaw
 call wrtout(std_out,msg,mode_paral)

 write(msg, '(a,i3,33x,a,a,i3)' )'  usewvl',hdr %usewvl,'|','  usewvl',hdr0%usewvl
 call wrtout(std_out,msg,mode_paral)

 write(msg,'(a,35x,a,a,3(a1,2x,3f12.7,6x,a,2x,3f12.7))')&
  '  rprimd:','|','  rprimd:',ch10,&
  hdr%rprimd(:,1),'|',hdr0%rprimd(:,1),ch10,&
  hdr%rprimd(:,2),'|',hdr0%rprimd(:,2),ch10,&
  hdr%rprimd(:,3),'|',hdr0%rprimd(:,3)
 call wrtout(std_out,msg,mode_paral)

 if (hdr%bantot/=hdr0%bantot) tband=1

 if (hdr%intxc/=hdr0%intxc) then
   write(msg,'(a,i0,a,i0)')'input intxc=',hdr%intxc,' not equal disk file intxc=',hdr0%intxc
   MSG_WARNING(msg)
 end if

 if (hdr%ixc/=hdr0%ixc) then
   write(msg,'(a,i0,a,i0)')'input ixc=',hdr%ixc,' not equal disk file ixc=',hdr0%ixc
   MSG_WARNING(msg)
 end if

 if (hdr%natom/=hdr0%natom) then
   write(msg,'(a,i0,a,i0)')'input natom=',hdr%natom,' not equal disk file natom=',hdr0%natom
   MSG_WARNING(msg)
   tatty=1
 end if

 if ( ANY(hdr%ngfft/=hdr0%ngfft) ) then
   ! For sensible rho(r) or V(r) data, fft grid must be identical
   ! Note, however, that we allow for different FFT meshes and we interpolate the density in the
   ! caller when we are restarting a SCF calculation.
   if (abifile%class == "density" .or. abifile%class == "potential") then
     write(msg, '(10a)' )&
       'FFT grids must be the same to restart from a ',trim(abifile%class),' file.',ch10,&
       "ngfft from file: ", trim(ltoa(hdr0%ngfft(1:3))), ", from input: ", trim(ltoa(hdr%ngfft(1:3))), ch10, &
       'Action: change the FFT grid in the input via ngfft or change the restart file.'
     MSG_ERROR(msg)
   end if
   tng=1
 end if

 if (hdr%nkpt/=hdr0%nkpt) then
   if (abifile%class == "wf_planewave") then
     write(msg,'(a,i0,a,i0)' )'input nkpt=',hdr%nkpt,' not equal disk file nkpt=',hdr0%nkpt
     MSG_COMMENT(msg)
   end if
   tkpt=1; twfk=1
 end if

 if (hdr%nspinor/=hdr0%nspinor) then
   if (abifile%class == "wf_planewave") then
     write(msg,'(a,i0,a,i0)')'input nspinor=',hdr%nspinor,' not equal disk file nspinor=',hdr0%nspinor
     MSG_WARNING(msg)
   end if
   tspinor=1
 end if

 ! No check is present for nspden
 if (hdr%nsppol/=hdr0%nsppol) then
   write(msg,'(a,i0,a,i0)')'input nsppol=',hdr%nsppol,' not equal disk file nsppol=',hdr0%nsppol
   MSG_WARNING(msg)
 end if

 if (hdr%nsym/=hdr0%nsym) then
   write(msg, '(a,i0,a,i0)' )'input nsym=',hdr%nsym,' not equal disk file nsym=',hdr0%nsym
   MSG_WARNING(msg)
   tsym=1
 end if

 if (hdr%ntypat/=hdr0%ntypat) then
   write(msg,'(a,i0,a,i0)')'input ntypat=',hdr%ntypat,' not equal disk file ntypat=',hdr0%ntypat
   call wrtout(std_out,msg,mode_paral)
   MSG_WARNING(msg)
   tatty=1
 end if

 if (hdr%usepaw/=hdr0%usepaw) then
   write(msg,'(a,i0,a,i0)')'input usepaw=',hdr%usepaw,' not equal disk file usepaw=',hdr0%usepaw
   MSG_WARNING(msg)
   tpaw=1
 end if

 if (hdr%usewvl/=hdr0%usewvl) then
   write(msg, '(a,i6,a,i6,a,a)' )&
     'input usewvl=',hdr%usewvl,' not equal disk file usewvl=',hdr0%usewvl, ch10, &
     'Action: change usewvl input variable or your restart file.'
   MSG_ERROR(msg)
 end if

 ! Also examine agreement of floating point data
 if (hdr%usewvl == 0 .and. abs(hdr%ecut_eff-hdr0%ecut_eff)>tol8) then
   write(msg,'(a,f12.6,a,f12.6,a)')'input ecut_eff=',hdr%ecut_eff,' /= disk file ecut_eff=',hdr0%ecut_eff,'.'
   MSG_WARNING(msg)
   tecut=1
 end if

 do ii=1,3
   do jj=1,3
     if (abs(hdr%rprimd(ii,jj)-hdr0%rprimd(ii,jj))>tol6) then
       write(msg, '(a,i1,a,i1,a,1p,e17.9,a,i1,a,i1,a,e17.9)' )&
         'input rprimd(',ii,',',jj,')=',hdr%rprimd(ii,jj),' /= disk file rprimd(',ii,',',jj,')=',hdr0%rprimd(ii,jj)
       MSG_WARNING(msg)
       tprim=1
     end if
   end do
 end do

 ! Below this point many comparisons only make sense if
 ! certain things agree, e.g. nkpt, natom.  Also have to
 ! accomodate different amounts of data in general.

 if (hdr%usepaw==1 .and. hdr0%usepaw==1) then

   ! Compare ecutdg (PAW)
   write(msg, '(a,f12.6,19x,a,a,f12.6)' )'  PAW: ecutdg',hdr %ecutdg,'|','  PAW: ecutdg',hdr0%ecutdg
   call wrtout(std_out,msg,mode_paral)
   if (hdr%ecutdg/=hdr0%ecutdg) then
     write(msg, '(a,f12.6,a,f12.6)' )'input ecutdg=',hdr%ecutdg,'not equal disk file ecutdg=',hdr0%ecutdg
     MSG_WARNING(msg)
     tdg=1
   end if
 end if

 ! Compare nband(nkpt*nsppol) (cannot compare if nkpt and nsppol not same)
 if (hdr%nkpt==hdr0%nkpt .and. hdr%nsppol==hdr0%nsppol) then
   nkpt=hdr%nkpt ; nsppol=hdr%nsppol
   write(msg,'(a,36x,a,a)') '  nband:','|','  nband:'
   call wrtout(std_out,msg,mode_paral)
   do istart = 1,nsppol*nkpt,9
     istop = min(istart + 8,nsppol*nkpt)
     mu = istop - istart + 1
     ! generate a format specifier
     bndfmt = strcat('(2x,',itoa(mu),'i4,t41,a,2x,',itoa(mu),'i4)')
     if (istart<=100) then
       write(msg,fmt=bndfmt) hdr%nband(istart:istop),'    |',hdr0%nband(istart:istop)
       call wrtout(std_out,msg,mode_paral)
       if (istop>100) call wrtout(std_out, '=> stop printing nband after 100 values', mode_paral)
     end if
   end do

   enough = 0
   do isppol=1,nsppol
     do ii=1,nkpt
       if (hdr%nband(ii)/=hdr0%nband(ii)) then
         tband=1
         enough = enough + 1
         if (abifile%class == "wf_planewave") then
           if (enough > 5) then
              write(std_out, "(a)")"Stop writing warnings after 5 values"
              exit
           else
             write(msg,'(a,i0,a,i0,a,i0)' )&
              'kpt num ',ii,' input nband= ',hdr%nband(ii),' not equal disk file nband=',hdr0%nband(ii)
             MSG_WARNING(msg)
           end if
         end if
       end if
     end do
   end do
 end if

 ! Compare the number of wavelets in each resolution.
 if (hdr%usewvl == 1) then
   if (size(hdr%nwvlarr) /= size(hdr0%nwvlarr) .or. size(hdr%nwvlarr) /= 2) then
     write(msg, '(a,i0,a,i0,a,a)' )&
      'input nwvlres= ',size(hdr%nwvlarr),' not equal disk file nwvlres= ',size(hdr0%nwvlarr),' or 2',&
      ' ABINIT is not implemented for wavelet resolutions different from 2.'
     MSG_ERROR(msg)
   end if
 end if

 ! Compare symmetry arrays (integers) symafm(nsym)
 ! only for same number of symmetries nsym
 itest=0
 if (hdr%nsym==hdr0%nsym) then
   nsym=hdr%nsym
   write(msg,'(a,35x,a,a)') '  symafm:','|','  symafm:'
   call wrtout(std_out,msg,mode_paral)
   do istart = 1,nsym,12
     istop=min(istart+11,nsym)
     nelm = istop - istart + 1
     typfmt = strcat('(2x,',itoa(nelm),'i3,t41,a,2x,',itoa(nelm),'i3)')
     write(msg,fmt=typfmt) hdr%symafm(istart:istop),'    |',hdr0%symafm(istart:istop)
     call wrtout(std_out,msg,mode_paral)
   end do
 end if

 if (itest/=0) then
   write(msg,'(a,i0,a)' )'For symmetry number',itest,' input symafm not equal disk file symafm'
   MSG_WARNING(msg)
   tsym=1
 end if

 ! Compare symmetry arrays (integers) symrel(3,3,nsym)
 ! only for same number of symmetries nsym
 itest=0
 if (hdr%nsym==hdr0%nsym) then
   nsym=hdr%nsym
   write(msg,'(a,35x,a,a)') '  symrel:','|','  symrel:'
   call wrtout(std_out,msg,mode_paral)
   do isym=1,nsym
     write(msg,'(2x,9i3,15x,a,2x,9i3)')hdr%symrel(:,:,isym),'|',hdr0%symrel(:,:,isym)
     call wrtout(std_out,msg,mode_paral)
     if(sum(abs(hdr%symrel(:,:,isym)-hdr0%symrel(:,:,isym)))/=0)then
       itest=isym
       exit
     end if
   end do
 end if

 if (itest/=0) then
   write(msg,'(a,i0,a)')'For symmetry number',itest,' input symrel not equal disk file symrel'
   MSG_WARNING(msg)
   tsym=1
 end if

 ! Compare typat(natom)
 if (hdr%natom==hdr0%natom) then
   natom=hdr%natom
   write(msg,'(a,36x,a,a)') '  typat:','|','  typat:'
   call wrtout(std_out,msg,mode_paral)
   do istart = 1,natom,12
     istop=min(istart+11,natom)
     nelm = istop - istart + 1
     typfmt = strcat('(2x,',itoa(nelm),'i3,t41,a,2x,',itoa(nelm),'i3)')
     write(msg,fmt=typfmt) hdr%typat(istart:istop),'    |',hdr0%typat(istart:istop)
     call wrtout(std_out,msg,mode_paral)
   end do
   do ii=1,natom
     if (hdr%typat(ii)/=hdr0%typat(ii)) then
       write(msg, '(a,i0,a,i0,a,i0)' )&
        'For atom number ',ii,' input typat=',hdr%typat(ii),' not equal disk file typat=',hdr0%typat(ii)
       MSG_WARNING(msg)
       tatty=1
     end if
   end do
 end if

 ! Compare so_psp(npsp)
 if (hdr%npsp==hdr0%npsp) then
   npsp=hdr%npsp
   write(msg,'(a,33x,a,a)') '  so_psp  :','|','  so_psp  :'
   call wrtout(std_out,msg,mode_paral)
   do istart = 1,npsp  ,12
     istop=min(istart+11,npsp  )
     nelm = istop - istart + 1
     typfmt = strcat('(2x,',itoa(nelm),'i3,t41,a,2x,',itoa(nelm),'i3)')
     write(msg,fmt=typfmt) hdr%so_psp  (istart:istop),'    |',hdr0%so_psp  (istart:istop)
     call wrtout(std_out,msg,mode_paral)
   end do
   do ii=1,npsp
     if (hdr%so_psp  (ii)/=hdr0%so_psp  (ii)) then
       write(msg,'(a,i0,a,i0,a,i0)')&
         'For pseudopotential number ',ii,' input so_psp =',hdr%so_psp(ii),' not equal disk file so_psp=',hdr0%so_psp(ii)
       MSG_WARNING(msg)
     end if
   end do
 end if

 ! Compare istwfk(nkpt)
 if (hdr%nkpt==hdr0%nkpt) then
   nkpt=hdr%nkpt
   write(msg,'(a,35x,a,a)') '  istwfk:','|','  istwfk:'
   call wrtout(std_out,msg,mode_paral)
   do istart = 1,nkpt,12
     istop=min(istart+11,nkpt)
     nelm = istop - istart + 1
     typfmt = strcat('(2x,',itoa(nelm),'i3,t41,a,2x,',itoa(nelm),'i3)')
     if (istart<=100) then
       write(msg,fmt=typfmt) hdr%istwfk(istart:istop),'    |',hdr0%istwfk(istart:istop)
       call wrtout(std_out,msg,mode_paral)
       if (istop>100) then
         call wrtout(std_out, '=> stop printing istwfk after 100 values' ,mode_paral)
       end if
     end if
   end do
   do ii=1,nkpt
     if (hdr%istwfk(ii)/=hdr0%istwfk(ii)) then
       write(msg, '(a,i0,a,i0,a,i0)' )&
         'For k point number ',ii,' input istwfk=',hdr%istwfk(ii),' not equal disk file istwfk=',hdr0%istwfk(ii)
       MSG_COMMENT(msg)
       twfk=1
     end if
   end do
 end if

!NEW_HDR
 if (any(hdr%kptrlatt /= hdr0%kptrlatt)) then
    write(msg,"(2(a,9(i0,1x)))")"input kptrlatt = ",hdr%kptrlatt," /= disk file kptrlatt = ",hdr0%kptrlatt
    MSG_COMMENT(msg)
 end if
 if (hdr%kptopt /= hdr0%kptopt) then
    MSG_COMMENT(sjoin("input kptopt = ", itoa(hdr%kptopt)," /= disk file kptopt = ", itoa(hdr0%kptopt)))
 end if
 if (hdr%pawcpxocc /= hdr0%pawcpxocc) then
    MSG_WARNING(sjoin("input pawcpxocc = ", itoa(hdr%pawcpxocc)," /= disk file pawcpxocc = ", itoa(hdr0%pawcpxocc)))
 end if
 if (hdr%icoulomb /= hdr0%icoulomb) then
    MSG_WARNING(sjoin("input icoulomb = ", itoa(hdr%icoulomb)," /= disk file icoulomb = ", itoa(hdr0%icoulomb)))
 end if

 if (abs(hdr%nelect - hdr0%nelect) > tol6) then
    MSG_WARNING(sjoin("input nelect = ", ftoa(hdr%nelect)," /= disk file nelect = ",ftoa(hdr0%nelect)))
 end if
 if (abs(hdr%charge - hdr0%charge) > tol6) then
    MSG_WARNING(sjoin("input charge = ", ftoa(hdr%charge)," /= disk file charge = ", ftoa(hdr0%charge)))
 end if

 if (hdr%ntypat==hdr0%ntypat) then
   if (any(abs(hdr%amu - hdr0%amu) > tol6)) then
      MSG_WARNING(sjoin("input amu = ",ltoa(hdr%amu)," /= disk file amu = ",ltoa(hdr0%amu)))
   end if
 end if
!end NEW_HDR

 ! Compare kpt(3,nkpt)
 if (hdr%nkpt==hdr0%nkpt) then
   nkpt=hdr%nkpt
   write(msg,'(a,38x,a,a)') '  kpt:','|','  kpt:'
   call wrtout(std_out,msg,mode_paral)
   do ii = 1,min(nkpt,nkpt_max)
     write(msg,'(2x,3f12.7,2x,a,2x,3f12.7)')hdr%kptns(:,ii),'    |',hdr0%kptns(:,ii)
     call wrtout(std_out,msg,mode_paral)
     if(ii>nkpt_max)then
       call wrtout(std_out,'The number of printed k points is sufficient... stop writing them.',mode_paral)
       exit
     end if
   end do
   iwarning=0
   do ii=1,nkpt
     itest=0
     do mu=1,3
       if(abs( hdr%kptns(mu,ii)-hdr0%kptns(mu,ii) )>tol6)itest=1
     end do
     if (itest==1) then
       write(msg, '(a,i5,a,3es17.7,a,a,3es17.7)' )&
        'kpt num',ii,', input kpt=',hdr%kptns(:,ii),ch10,&
        'not equal  disk file kpt=',hdr0%kptns(:,ii)
       MSG_WARNING(msg)
       tkpt=1 ; iwarning=iwarning+1
       if(iwarning>=mwarning)then
         call wrtout(std_out,'The number of warning messages is sufficient ... stop writing them.',mode_paral)
         exit
       end if
     end if
   end do
 end if

 ! Compare wtk(nkpt)
 if (hdr%nkpt==hdr0%nkpt) then
   nkpt=hdr%nkpt

   write(msg,'(a,38x,a,a)') '  wtk:','|','  wtk:'
   call wrtout(std_out,msg,mode_paral)
   istop = min(nkpt,nkpt_max)
   do ii = 1, istop, 5
     mu = min(5, istop - ii + 1)
     wtkfmt = strcat('(2x,',itoa(mu),'f7.3,t41,a,2x,',itoa(mu),'f7.3)')
     write(msg, wtkfmt)hdr%wtk(ii:min(istop, ii + 5 - 1)),'    |',hdr0%wtk(ii:min(istop, ii + 5 - 1))
     call wrtout(std_out,msg,mode_paral)
   end do
   iwarning=0
   do ii=1,nkpt
     itest=0
     if (abs( hdr%wtk(ii)-hdr0%wtk(ii) )>tol6) then
       write(msg,'(a,i5,a,es17.7,a,a,es17.7)')&
        'kpt num',ii,', input weight=',hdr%wtk(ii),ch10,&
        'not equal to disk file weight=',hdr0%wtk(ii)
       MSG_WARNING(msg)

       tkpt=1 ; iwarning=iwarning+1
       if(iwarning>=mwarning)then
         call wrtout(std_out,'The number of warning messages is sufficient ... stop writing them.',mode_paral)
         exit
       end if
     end if
   end do
 end if

 ! Compare occ(bantot)
 if (hdr%nkpt==hdr0%nkpt.and. hdr%bantot==hdr0%bantot) then
   nkpt=hdr%nkpt
   bantot=hdr%bantot

   write(msg,'(a,38x,a,a)') '  occ:','|','  occ:'
   call wrtout(std_out,msg,mode_paral)
   bantot_eff=min(bantot,9*nkpt_max)
   do istart = 1,bantot_eff,9
     istop = min(istart+8,bantot_eff)
     mu = istop - istart + 1
     occfmt = strcat('(2x,',itoa(mu),'f4.1,t41,a,2x,',itoa(mu),'f4.1)')
     write(msg,fmt=occfmt)hdr%occ(istart:istop),'    |', hdr0%occ(istart:istop)
     call wrtout(std_out,msg,mode_paral)
     if(istart>9*nkpt_max)then
       call wrtout(std_out,'The number of printed occupation numbers is sufficient ... stop writing them.',mode_paral)
       exit
     end if
   end do
   iwarning=0
   do ii=1,bantot
     if (abs( hdr%occ(ii)-hdr0%occ(ii) )>tol6) then
       write(msg,'(a,i0,a,1p,e15.7,a,e15.7)')'band,k: ',ii,', input occ=',hdr%occ(ii),' disk occ=',hdr0%occ(ii)
       MSG_WARNING(msg)
       tband=1 ; iwarning=iwarning+1
       if(iwarning>=mwarning)then
         call wrtout(std_out,'The number of warning msgs is sufficient ... stop writing them.',mode_paral)
         exit
       end if
     end if
   end do
 end if

 ! Compare tnons(3,nsym)
 if (hdr%nsym==hdr0%nsym) then
   nsym=hdr%nsym
   itest=0
   write(msg,'(a,36x,a,a)') '  tnons:','|','  tnons:'
   call wrtout(std_out,msg,mode_paral)
   do isym=1,nsym
     write(msg,'(2x,3f12.7,2x,a,2x,3f12.7)') hdr%tnons(:,isym),'    |',hdr0%tnons(:,isym)
     call wrtout(std_out,msg,mode_paral)
   end do

   do isym=1,nsym
     if( sum(abs(  hdr%tnons(:,isym)-hdr0%tnons(:,isym) )) > tol6) then
       itest=isym
       exit
     end if
   end do
   if (itest/=0) then
     write(msg, '(a,i0,a)' )'For symmetry number ',itest,' input tnons not equal disk file tnons'
     MSG_WARNING(msg)
   end if
 end if

 ! Compare znucltypat(ntypat)
 if (hdr%ntypat==hdr0%ntypat) then
   ntypat=hdr%ntypat

   write(msg,'(a,35x,a,a)') '   znucl:','|','   znucl:'
   call wrtout(std_out,msg,mode_paral)
   do istart = 1,ntypat,6
     istop = min(istart+5,ntypat)
     mu = istop-istart+1
     zatfmt = strcat('(2x,',itoa(mu),'f6.2,t41,a,6x,',itoa(mu),'f6.2)')
     write(msg,fmt=zatfmt) hdr%znucltypat(istart:istop),'    |',hdr0%znucltypat(istart:istop)
     call wrtout(std_out,msg,mode_paral)
   end do

   do ii=1,ntypat
     if (abs(hdr%znucltypat(ii)-hdr0%znucltypat(ii))>tol6) then
       write(msg, '(a,i5,a,f12.6,a,f12.6)' )&
        ' For atom number ',ii,' input znucl=',hdr%znucltypat(ii),' not equal disk file znucl=',hdr0%znucltypat(ii)
       MSG_WARNING(msg)
     end if
   end do
 end if

 ! Should perform some checks related to pertcase and qptn,
 ! that have been introduced in the header in v4.1
 ! Warning: a GS file might be read, while the hdr corresponds
 ! to a RF file (to initialize k+q), and vice-versa (in nonlinear).

 ! Now check agreement of psp headers too
 if (hdr%npsp==hdr0%npsp) then
   npsp=hdr%npsp
   itest=0

   do ipsp=1,npsp
     write(msg,'(a,i3,a,13x,a,a,i3,a)')&
      '  pseudopotential atom type',ipsp,':','|','  pseudopotential atom type',ipsp,':'
     call wrtout(std_out,msg,mode_paral)

     if (hdr%usepaw==1 .and. hdr0%usepaw==1) then
       write(msg,'(a,i3,a,i7,a,i3,5x,a,a,i3,a,i7,a,i3)')&
        '  pspso ',hdr %pspso(ipsp),' pspxc ',hdr %pspxc(ipsp),&
        '  lmn_size ',hdr%lmn_size(ipsp),'|',&
        '  pspso ',hdr0%pspso(ipsp),' pspxc ',hdr0%pspxc(ipsp),&
        '  lmn_size ',hdr0%lmn_size(ipsp)
       call wrtout(std_out,msg,mode_paral)
       if (hdr%lmn_size(ipsp)/=hdr0%lmn_size(ipsp)) then
         write(msg, '(a,i3,a,i3,a,i3)' )&
          'For atom type ',ipsp,' input lmn_size=',hdr%lmn_size(ipsp),&
          'not equal disk file lmn_size=',hdr0%lmn_size(ipsp)
         MSG_WARNING(msg)
         tlmn=1
       end if
     else
       write(msg,'(a,i3,a,i3,23x,a,a,i3,a,i3)')&
        '  pspso ',hdr %pspso(ipsp),' pspxc ',hdr %pspxc(ipsp),'|',&
        '  pspso ',hdr0%pspso(ipsp),' pspxc ',hdr0%pspxc(ipsp)
       call wrtout(std_out,msg,mode_paral)
     end if
     write(msg,'(a,i8,a,i4,a,f5.1,4x,a,a,i8,a,i4,a,f5.1)')&
      '  pspdat ',hdr %pspdat(ipsp),' pspcod ',hdr %pspcod(ipsp),&
      ' zion ',hdr %zionpsp(ipsp),'|',&
      '  pspdat ',hdr0%pspdat(ipsp),' pspcod ',hdr0%pspcod(ipsp),&
      ' zion ',hdr0%zionpsp(ipsp)
     call wrtout(std_out,msg,mode_paral)

     ! Check on md5 values.
     if (hdr%md5_pseudos(ipsp) /= hdr0%md5_pseudos(ipsp)) then
       write(msg, '(a,i0,6a)' )&
       ' Different md5 checksum for pseudo ',ipsp,ch10,&
       ' input md5= ',hdr%md5_pseudos(ipsp),ch10,&
       ' disk  md5= ',hdr0%md5_pseudos(ipsp)
       MSG_WARNING(msg)
       itest=1; tpsch=1
     end if

     ! Second, test
     ! NOTE, XG 000719: should do something about pspso
     ! NOTE, XG 020716: znucl and zion are not written
     if (abs(hdr%znuclpsp(ipsp)-hdr0%znuclpsp(ipsp))>tol6) itest=1
     if (abs(hdr%zionpsp(ipsp)-hdr0%zionpsp(ipsp))>tol6) then
       itest=1; tpsch=1
     end if
     if (hdr%pspdat(ipsp)/= hdr0%pspdat(ipsp)) itest=1
     if (hdr%pspcod(ipsp)/= hdr0%pspcod(ipsp)) itest=1
     if (hdr%pspxc(ipsp) /= hdr0%pspxc(ipsp) )  itest=1
   end do

   if (itest==1) then
     MSG_WARNING('input psp header does not agree perfectly with disk file psp header.')
     tpseu=1
   end if
 end if

 ! Finally, read residm and etotal ("current value" not known), and check xred.
 if (hdr%natom==hdr0%natom) then
   natom=hdr%natom
   write(msg,'(a,37x,a,a)') '  xred:','|','  xred:'
   call wrtout(std_out,msg,mode_paral)
   do ii=1,natom
     write(msg,'(2x,3f12.7,6x,a,2x,3f12.7)') hdr%xred(:,ii),'|',hdr0%xred(:,ii)
     call wrtout(std_out,msg,mode_paral)
   end do

   ! check atom positions one atom at a time and allow possibility
   ! that there is a harmless translation of atoms by a cell vector.
   do ii=1,natom
     rms=0.0_dp
     do jj=1,3
       rms=rms+(hdr%xred(jj,ii)-hdr0%xred(jj,ii) - dble(nint((hdr%xred(jj,ii)-hdr0%xred(jj,ii)))) )**2
     end do
     rms=sqrt(rms/3.0_dp)
     if (rms>tol6) txred=1
   end do
 end if

 ! Run tests here to establish whether this is a valid restart

 ! tfform2 will be true if there is a problem for the wavefunctions
 tfform2 = (hdr%usewvl == 0 .and. &
           (tprim /= 0 .or. tecut /= 0 .or. tkpt /= 0 .or. &
           twfk /=0 .or. tspinor /= 0)) .or. &
           (hdr%usepaw == 1 .and. &
           (tpaw /= 0 .or. tlmn /= 0 .or. tdg /= 0)) .or. &
           (hdr%usewvl == 1 .and. &
           (tatty /= 0 .or. tband /= 0))

 ! tfform52 will be true if there is a problem for the format 52
 tfform52 = tprim /= 0 .or. tatty /= 0 .or. txred /= 0 .or.&
            tpseu /= 0 .or. tecut /= 0 .or. tng /= 0 .or. &
            (hdr%usepaw == 1 .and. (tpaw /= 0 .or. tlmn /= 0 .or. tdg /= 0))

 restart=1; restartpaw=hdr%usepaw

 ! If there is a problem somewhere
 if ( (abifile%class == "wf_planewave"  .and. tfform2  ) .or.  &
      (abifile%class == "density" .and. tfform52 ) .or.  &
      (abifile%class == "wf_wavelet" .and. tfform2 ) ) then

   if (abifile%class == "wf_planewave") then
     restart=2
     MSG_COMMENT('Restart of self-consistent calculation need translated wavefunctions.')
   else if (abifile%class == "density") then
     restart=0
     MSG_WARNING('Illegal restart of non-self-consistent calculation')
   end if

   write(msg,'(a,a1,a)') &
     '  Indeed, critical differences between current calculation and',ch10,&
     '  restart file have been detected in:'
   call wrtout(std_out,msg,mode_paral)

   if ( (abifile%class == "density" .or. abifile%class == "wf_wavelet") .and. tatty /= 0 ) then
     write(msg, '(8x,a)' ) '* the number of atoms of each type'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( abifile%class /= "wf_wavelet" .and. tecut /= 0 ) then
     write(msg, '(8x,a)' ) '* the plane-wave cutoff'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( abifile%class == "wf_wavelent" .and. tband /= 0 ) then
     write(msg, '(8x,a)' ) '* the band and their occupation'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( abifile%class == "wf_planewave" .and. tkpt /= 0 ) then
     write(msg, '(8x,a)' ) '* the number, position, or weight of k-points'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( abifile%class == "wf_planewave" .and. twfk /= 0 ) then
     write(msg, '(8x,a)' ) '* the format of wavefunctions (istwfk)'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( abifile%class == "wf_planewave"  .and. tspinor /= 0 ) then
     write(msg, '(8x,a)' ) '* the scalar/spinor character of the wf (nspinor)'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( abifile%class == "density"  .and. tng /= 0 ) then
     write(msg, '(8x,a)' ) '* the Fourier transform box dimensions'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( tprim /= 0 ) then
     write(msg, '(8x,a)' )'* the vectors defining the unit cell (obtained from rprim and acell)'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( abifile%class == "density"   .and. tpseu /= 0 ) then
     write(msg, '(8x,a)' )'* the pseudopotential files'
     call wrtout(std_out,msg,mode_paral)
   end if
   if ( abifile%class == "density"  .and. txred /= 0 ) then
     write(msg, '(8x,a)' ) '* the positions of the ions in the basis'
     call wrtout(std_out,msg,mode_paral)
   end if

   ! Tests for a restart in the framework of the PAW method
   if (hdr%usepaw/=0 .or. hdr0%usepaw/=0) then
     if (tpaw /= 0 .or. tlmn /= 0) restartpaw=0
     if (restartpaw == 0) then
       write(msg,'(8x,a)') 'Critical differences for a restart within PAW method:'
       call wrtout(std_out,msg,mode_paral)
       if ( tpaw /= 0 ) then
         write(msg, '(8x,a)' ) '* the use of the PAW method'
         call wrtout(std_out,msg,mode_paral)
       else
         if(tlmn/=0)then
           write(msg, '(8x,a)' ) '* the number of lmn elements for the paw basis'
           call wrtout(std_out,msg,mode_paral)
         end if
       end if
     else if (tdg/=0) then
       write(msg,'(a,a,a,a,a,a)') ch10,&
         ' hdr_check: WARNING -',ch10,&
         '  Restart of calculation within PAW may be inconsistent because of:"'
       call wrtout(std_out,msg,mode_paral)
       if(tdg/=0)then
         write(msg, '(8x,a)' )'* the cutoff energy of the paw double (fine) grid'
         call wrtout(std_out,msg,mode_paral)
       end if
     end if
   end if

 else

   if (abifile%class == "wf_planewave" .or. abifile%class == "wf_wavelet") then
     write(msg,'(a,a)') ' hdr_check: ',' Wavefunction file is OK for direct restart of calculation'
     call wrtout(std_out,msg,mode_paral)
   else if (abifile%class == "density") then
     write(msg,'(a,a)') ' hdr_check: ',' Density/Potential file is OK for restart of calculation'
     call wrtout(std_out,msg,mode_paral)
   end if
 end if

 write(msg,'(80a)') ('=',ii=1,80)
 call wrtout(std_out,msg,mode_paral)

end subroutine hdr_check

!!***

!----------------------------------------------------------------------

!!****f* m_wfk/hdr_compare
!! NAME
!!  hdr_compare
!!
!! FUNCTION
!!  Test two hdr_t objects for consistency. Return non-zero value if test fails.
!!
!! INPUTS
!!  hdr1, hdr2 <class(hdr_t)> = hdr handlers to be compared
!!
!! OUTPUT
!!  ierr
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function hdr_compare(hdr1, hdr2) result(ierr)

!Arguments ------------------------------------
!scalars
 class(hdr_type),intent(in) :: hdr1, hdr2

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

!************************************************************************

 ierr = 0

 ! Test basic dimensions
 if (hdr1%nsppol /= hdr2%nsppol) then
   write(msg,'(a,i0,a,i0)')'Different nsppol : ',hdr1%nsppol,' and ',hdr2%nsppol
   ierr = ierr + 1; MSG_WARNING(msg)
 end if
 if (hdr1%nspinor /= hdr2%nspinor) then
   write(msg,'(a,i0,a,i0)')'Different nspinor : ',hdr1%nspinor,' and ',hdr2%nspinor
   ierr = ierr + 1; MSG_WARNING(msg)
 end if
 if (hdr1%nspden /= hdr2%nspden) then
   write(msg,'(a,i0,a,i0)')'Different nspden : ',hdr1%nspden,' and ',hdr2%nspden
   ierr = ierr + 1; MSG_WARNING(msg)
 end if
 if (hdr1%nkpt /= hdr2%nkpt) then
   write(msg,'(a,i0,a,i0)')'Different nkpt : ',hdr1%nkpt,' and ',hdr2%nkpt
   ierr = ierr + 1; MSG_WARNING(msg)
 end if
 if (hdr1%usepaw /= hdr2%usepaw) then
   write(msg,'(a,i0,a,i0)')'Different usepaw : ',hdr1%usepaw,' and ',hdr2%usepaw
   ierr = ierr + 1; MSG_WARNING(msg)
 end if
 if (hdr1%ntypat /= hdr2%ntypat) then
   write(msg,'(a,i0,a,i0)')'Different ntypat : ',hdr1%ntypat,' and ',hdr2%ntypat
   ierr = ierr + 1; MSG_WARNING(msg)
 end if
 if (hdr1%natom /= hdr2%natom) then
   write(msg,'(a,i0,a,i0)')'Different natom  : ',hdr1%natom,' and ',hdr2%natom
   ierr = ierr + 1; MSG_WARNING(msg)
 end if

 ! Return immediately if important dimensions are not equal.
 if (ierr /= 0) return

 ! Test important arrays (rprimd is not tested)
 if (any(hdr1%typat /= hdr2%typat)) then
   write(msg,'(a,i0,a,i0)')'Different ntypat array : ',hdr1%typat(1),' ... and ',hdr2%typat(1)
   ierr = ierr + 1; MSG_WARNING(msg)
 end if
!Should test npwarr, however taking into account differences due to istwfk !
!if (any(hdr1%npwarr /= hdr2%npwarr)) then
!  write(msg,'(a,i0,a,i0)')'Different npwarr array : ',hdr1%npwarr(1),' ... and ',hdr2%npwarr(1)
!  ierr = ierr + 1; MSG_WARNING(msg)
!end if
 if (any(abs(hdr1%kptns - hdr2%kptns) > tol6)) then
   write(msg,'(a,i0,a,i0)')'Different kptns array '
   ierr = ierr + 1; MSG_WARNING(msg)
 end if

end function hdr_compare
!!***

!----------------------------------------------------------------------

!!****f* m_hdr/hdr_vs_dtset
!! NAME
!! hdr_vs_dtset
!!
!! FUNCTION
!!  Check the compatibility of the Abinit header with respect to the
!!  input variables defined in the input file.
!!
!! INPUTS
!!  Dtset<type(dataset_type)>=all input variables for this dataset
!!  Hdr <type(hdr_type)>=the header structured variable
!!
!! OUTPUT
!!  Only check
!!
!! PARENTS
!!      eph,setup_bse,setup_screening,setup_sigma,wfk_analyze
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine hdr_vs_dtset(Hdr,Dtset)

!Arguments ------------------------------------
 class(Hdr_type),intent(in) :: Hdr
 type(Dataset_type),intent(in) :: Dtset

!Local variables-------------------------------
 integer :: ik, jj, ierr
 logical :: test, tsymrel,ttnons, tsymafm
 character(len=5000) :: msg
! *************************************************************************

 ! Check basic dimensions
 ierr = 0
 call compare_int('natom',  Hdr%natom,  Dtset%natom,  ierr)
 call compare_int('nkpt',   Hdr%nkpt,   Dtset%nkpt,   ierr)
 call compare_int('npsp',   Hdr%npsp,   Dtset%npsp,   ierr)
 call compare_int('nspden', Hdr%nspden, Dtset%nspden, ierr)
 call compare_int('nspinor',Hdr%nspinor,Dtset%nspinor,ierr)
 call compare_int('nsppol', Hdr%nsppol, Dtset%nsppol, ierr)
 call compare_int('nsym',   Hdr%nsym,   Dtset%nsym,   ierr)
 call compare_int('ntypat', Hdr%ntypat, Dtset%ntypat, ierr)
 call compare_int('usepaw', Hdr%usepaw, Dtset%usepaw, ierr)
 call compare_int('usewvl', Hdr%usewvl, Dtset%usewvl, ierr)
 call compare_int('kptopt', Hdr%kptopt, Dtset%kptopt, ierr)
 call compare_int('pawcpxocc', Hdr%pawcpxocc, Dtset%pawcpxocc, ierr)
 call compare_int('nshiftk_orig', Hdr%nshiftk_orig, Dtset%nshiftk_orig, ierr)
 call compare_int('nshiftk', Hdr%nshiftk, Dtset%nshiftk, ierr)

 ! The number of fatal errors must be zero.
 if (ierr/=0) then
   write(msg,'(3a)')&
   'Cannot continue, basic dimensions reported in the header do not agree with input file. ',ch10,&
   'Check consistency between the content of the external file and the input file.'
   MSG_ERROR(msg)
 end if

 test=ALL(ABS(Hdr%xred-Dtset%xred_orig(:,1:Dtset%natom,1))<tol6)
 ABI_CHECK(test,'Mismatch in xred')

 test=ALL(Hdr%typat==Dtset%typat(1:Dtset%natom))
 ABI_CHECK(test,'Mismatch in typat')

 ! Check if the lattice from the input file agrees with that read from the KSS file
 if ( (ANY(ABS(Hdr%rprimd - Dtset%rprimd_orig(1:3,1:3,1)) > tol6)) ) then
   write(msg,'(5a,3(3es16.6),3a,3(3es16.6),3a)')ch10,&
   ' real lattice vectors read from Header differ from the values specified in the input file', ch10, &
   ' rprimd from Hdr file   = ',ch10,(Hdr%rprimd(:,jj),jj=1,3),ch10,&
   ' rprimd from input file = ',ch10,(Dtset%rprimd_orig(:,jj,1),jj=1,3),ch10,ch10,&
   ' Modify the lattice vectors in the input file '
   MSG_ERROR(msg)
 end if

 ! Check symmetry operations.
 tsymrel=(ALL(Hdr%symrel==Dtset%symrel(:,:,1:Dtset%nsym)))
 if (.not.tsymrel) then
   write(msg,'(3a)')&
   ' real space symmetries read from Header ',ch10,&
   ' differ from the values inferred from the input file'
   MSG_WARNING(msg)
   tsymrel=.FALSE.
 end if

 ttnons=ALL(ABS(Hdr%tnons-Dtset%tnons(:,1:Dtset%nsym))<tol6)
 if (.not.ttnons) then
   write(msg,'(3a)')&
   ' fractional translations read from Header ',ch10,&
   ' differ from the values inferred from the input file'
   MSG_WARNING(msg)
   ttnons=.FALSE.
 end if

 tsymafm=ALL(Hdr%symafm==Dtset%symafm(1:Dtset%nsym))
 if (.not.tsymafm) then
   write(msg,'(3a)')&
   ' AFM symmetries read from Header ',ch10,&
   ' differ from the values inferred from the input file'
   MSG_WARNING(msg)
   tsymafm=.FALSE.
 end if

 if (.not.(tsymrel.and.ttnons.and.tsymafm)) then
   write(msg,'(a)')' Header '
   call wrtout(std_out,msg,'COLL')
   call print_symmetries(Hdr%nsym,Hdr%symrel,Hdr%tnons,Hdr%symafm)
   write(msg,'(a)')' Dtset  '
   call wrtout(std_out,msg,'COLL')
   call print_symmetries(Dtset%nsym,Dtset%symrel,Dtset%tnons,Dtset%symafm)
   MSG_ERROR('Check symmetry operations')
 end if

 if (abs(Dtset%nelect-hdr%nelect)>tol6) then
   write(msg,'(2(a,f8.2))')"File contains ", hdr%nelect," electrons but nelect initialized from input is ",Dtset%nelect
   MSG_ERROR(msg)
 end if
 if (abs(Dtset%charge-hdr%charge)>tol6) then
   write(msg,'(2(a,f8.2))')"File contains charge ", hdr%charge," but charge from input is ",Dtset%charge
   MSG_ERROR(msg)
 end if

 if (any(hdr%kptrlatt_orig /= dtset%kptrlatt_orig)) then
   write(msg,"(5a)")&
   "hdr%kptrlatt_orig: ",trim(ltoa(reshape(hdr%kptrlatt_orig,[9]))),ch10,&
   "dtset%kptrlatt_orig: ",trim(ltoa(reshape(dtset%kptrlatt_orig, [9])))
   MSG_ERROR(msg)
 end if

 if (any(hdr%kptrlatt /= dtset%kptrlatt)) then
   write(msg,"(5a)")&
   "hdr%kptrlatt: ",trim(ltoa(reshape(hdr%kptrlatt, [9]))),ch10,&
   "dtset%kptrlatt: ",trim(ltoa(reshape(dtset%kptrlatt, [9])))
   MSG_ERROR(msg)
 end if

 if (any(abs(hdr%shiftk_orig - dtset%shiftk_orig(:,1:dtset%nshiftk_orig)) > tol6)) then
   write(msg,"(5a)")&
   "hdr%shiftk_orig: ",trim(ltoa(reshape(hdr%shiftk_orig, [3*hdr%nshiftk_orig]))),ch10,&
   "dtset%shiftk_orig: ",trim(ltoa(reshape(dtset%shiftk_orig, [3*dtset%nshiftk_orig])))
   MSG_ERROR(msg)
 end if

 if (any(abs(hdr%shiftk - dtset%shiftk(:,1:dtset%nshiftk)) > tol6)) then
   write(msg,"(5a)")&
   "hdr%shiftk: ",trim(ltoa(reshape(hdr%shiftk, [3*hdr%nshiftk]))),ch10,&
   "dtset%shiftk: ",trim(ltoa(reshape(dtset%shiftk, [3*dtset%nshiftk])))
   MSG_ERROR(msg)
 end if

 ! Check if the k-points from the input file agrees with that read from the WFK file
 if ((ANY(ABS(Hdr%kptns(:,:) - Dtset%kpt(:,1:Dtset%nkpt)) > tol6))) then
   write(msg,'(9a)')ch10,&
   ' hdr_vs_dtset: ERROR - ',ch10,&
   '  k-points read from Header ',ch10,&
   '  differ from the values specified in the input file',ch10,&
   '  k-points from Hdr file                        | k-points from input file ',ch10
   call wrtout(std_out,msg,'COLL')
   do ik=1,Dtset%nkpt
     if (any(abs(Hdr%kptns(:,ik) - Dtset%kpt(:,ik)) > tol6)) then
       write(msg,'(3(3es16.6,3x))')Hdr%kptns(:,ik),Dtset%kpt(:,ik)
       call wrtout(std_out,msg,'COLL')
     end if
   end do
   MSG_ERROR('Modify the k-mesh in the input file')
 end if

 if (ANY(ABS(Hdr%wtk(:) - Dtset%wtk(1:Dtset%nkpt)) > tol6)) then
   write(msg,'(9a)')ch10,&
   ' hdr_vs_dtset : ERROR - ',ch10,&
   '  k-point weights read from Header ',ch10,&
   '  differ from the values specified in the input file',ch10,&
   '  Hdr file  |  File ',ch10
   call wrtout(std_out,msg,'COLL')
   do ik=1,Dtset%nkpt
     if (abs(Hdr%wtk(ik) - Dtset%wtk(ik)) > tol6) then
       write(msg,'(2(f11.5,1x))')Hdr%wtk(ik),Dtset%wtk(ik)
       call wrtout(std_out,msg,'COLL')
     end if
   end do
   MSG_ERROR('Check the k-mesh and the symmetries of the system. ')
 end if

 ! Check istwfk storage
 if ( (ANY(Hdr%istwfk(:)/=Dtset%istwfk(1:Dtset%nkpt))) ) then
   MSG_COMMENT('istwfk read from Header differs from the values specified in the input file (this is not critical)')
   !call wrtout(std_out, "  Hdr | input ")
   !do ik=1,Dtset%nkpt
   !  write(msg,'(i5,3x,i5)')Hdr%istwfk(ik),Dtset%istwfk(ik)
   !  call wrtout(std_out,msg,'COLL')
   !end do
   !MSG_ERROR('Modify istwfk in the input file.')
 end if

 CONTAINS
!!***

!!****f* hdr_vs_dtset/compare_int
!! NAME
!! compare_int
!!
!! FUNCTION
!!  Compare two int value and may raise an exception on error.
!!
!! INPUTS
!!  vname=Name of the variable
!!  iexp= expected value.
!!  ifound=the actuval value
!!
!! SIDE EFFECTS
!!  ierr=increased by one if values differ
!!
!! PARENTS
!!      hdr_vs_dtset
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

 subroutine compare_int(vname, iexp, ifound, ierr)

!Arguments ------------------------------------
 integer,intent(in) :: iexp,ifound
 integer,intent(inout) :: ierr
 character(len=*),intent(in) :: vname

!Local variables-------------------------------
 character(len=500) :: msg

! *************************************************************************

 if (.not. iexp == ifound) then
   write(msg,'(2a,i0,a,i0)')' Mismatch in '//trim(vname),' Expected = ', iexp, ' Found = ', ifound
   call wrtout(std_out, msg)
   ! Increase ierr to signal we should stop in the caller.
   ierr = ierr + 1
 end if

 end subroutine compare_int
!!***

end subroutine hdr_vs_dtset
!!***

!!****f* m_hdr/hdr_get_crystal
!! NAME
!!  hdr_get_crystal
!!
!! FUNCTION
!!  Initializes a crystal_t data type starting from the abinit header.
!!
!! INPUTS
!!  hdr<hdr_type>=the abinit header
!!  [gw_timrev] ==2 => take advantage of time-reversal symmetry
!!              ==1 ==> do not use time-reversal symmetry
!!    Default: 2
!!    NOTE THAT HERE WE USE THE GW CONVENTIONS  I.E ABINIT_TIMREV + !
!!  [remove_inv] = if .TRUE. the inversion symmetry is removed from the set of operations
!!      even if it is present in the header
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

type(crystal_t) function hdr_get_crystal(hdr, gw_timrev, remove_inv) result(cryst)

!Arguments ------------------------------------
 class(hdr_type),intent(in) :: hdr
 integer,optional,intent(in) :: gw_timrev
 logical,optional,intent(in) :: remove_inv

!Local variables-------------------------------
 integer :: my_timrev, space_group
 logical :: rinv, use_antiferro
! *********************************************************************

 rinv=.FALSE.; if (PRESENT(remove_inv)) rinv=remove_inv
 use_antiferro = hdr%nspden == 2 .and. hdr%nsppol ==1

 if (.not. present(gw_timrev)) then
   ! Get it from kptopt
   !my_timrev = kpts_timrev_from_kptopt(hdr%kptopt) + 1
   my_timrev = 1; if (any(hdr%kptopt == [3, 4])) my_timrev = 0
   my_timrev = my_timrev + 1
 else
   my_timrev = gw_timrev
 end if

 ! Consistency check
 ABI_CHECK(any(my_timrev == [1, 2]), "timrev should be in (1|2)")
 if (use_antiferro) then
   ABI_CHECK(ANY(hdr%symafm == -1), "Wrong nspden, nsppol, symafm.")
 end if

 space_group = 0 ! FIXME not known at this level.

 call crystal_init(hdr%amu,cryst,space_group,hdr%natom,hdr%npsp,hdr%ntypat,hdr%nsym,hdr%rprimd,hdr%typat,hdr%xred,&
   hdr%zionpsp,hdr%znuclpsp,my_timrev,use_antiferro,rinv,hdr%title,&
   symrel=hdr%symrel,tnons=hdr%tnons,symafm=hdr%symafm) ! Optional

end function hdr_get_crystal
!!***

end module m_hdr
!!***
