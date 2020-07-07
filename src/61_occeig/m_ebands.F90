!!****m* ABINIT/m_ebands
!! NAME
!!  m_ebands
!!
!! FUNCTION
!!  This module contains utilities to analyze and retrieve information from the ebands_t.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG, MJV, BXu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! TODO
!! 1) Remove npwarr, istwfk.
!! 2) Use 3d arrays for ebands%nband
!! 3) Solve issue with Hdr dependency
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_ebands

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_htetra
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr
 use m_krank
 use m_skw
 use m_kpts
 use m_sort
 use m_dtset

 use defs_datatypes,   only : ebands_t
 use m_copy,           only : alloc_copy
 use m_io_tools,       only : file_exists, open_file
 use m_time,           only : cwtime, cwtime_report
 use m_fstrings,       only : tolower, itoa, sjoin, ftoa, ltoa, ktoa, strcat, basename, replace
 use m_numeric_tools,  only : arth, imin_loc, imax_loc, bisect, stats_t, stats_eval, simpson_int, wrap2_zero_one, &
                              isdiagmat, interpol3d
 use m_special_funcs,  only : gaussian
 use m_geometry,       only : normv
 use m_cgtools,        only : set_istwfk
 use m_pptools,        only : printbxsf
 use m_occ,            only : getnel, newocc, occ_fd
 use m_nesting,        only : mknesting
 use m_crystal,        only : crystal_t
 use m_bz_mesh,        only : isamek, kpath_t, kpath_new
 use m_fftcore,        only : get_kg

 implicit none

 private

 ! Helper functions
 public :: pack_eneocc             ! Helper function for reshaping (energies|occupancies|derivate of occupancies).
 public :: get_eneocc_vect         ! Reshape (ene|occ|docdde) returning a matrix instead of a vector.
 public :: put_eneocc_vect         ! Put (ene|occ|doccde) in vectorial form into the data type doing a reshape.
 public :: unpack_eneocc           ! Helper function for reshaping (energies|occupancies|derivate of occupancies).

 ! Ebands methods
 public :: ebands_init             ! Main creation method.
 public :: ebands_from_hdr         ! Init object from the abinit header.
 public :: ebands_from_dtset       ! Init object from the abinit dataset.
 public :: ebands_free             ! Destruction method.
 public :: ebands_copy             ! Copy of the ebands_t.
 public :: ebands_move_alloc       ! Transfer allocation.
 public :: ebands_print            ! Printout basic info on the data type.
 public :: ebands_get_bandenergy   ! Returns the band energy of the system.
 public :: ebands_get_valence_idx  ! Gives the index of the (valence|bands at E_f).
 public :: ebands_get_bands_from_erange   ! Return the indices of the mix and max band within an energy window.
 public :: ebands_vcbm_range_from_gaps ! Find band and energy range for states close to the CBM/VBM given input energies.
 public :: ebands_apply_scissors   ! Apply scissors operator (no k-dependency)
 public :: ebands_get_occupied     ! Returns band indeces after wich occupations are less than an input value.
 public :: ebands_enclose_degbands ! Adjust band indeces such that all degenerate states are treated.
 public :: ebands_get_erange       ! Compute the minimum and maximum energy enclosing a list of states.
 public :: ebands_nelect_per_spin  ! Returns number of electrons per spin channel
 public :: ebands_get_minmax       ! Returns min and Max value of (eig|occ|doccde).
 public :: ebands_has_metal_scheme ! .True. if metallic occupation scheme is used.
 public :: ebands_write_bxsf       ! Write 3D energies for Fermi surface visualization (XSF format)
 public :: ebands_update_occ       ! Update the occupation numbers.
 public :: ebands_set_scheme       ! Set the occupation scheme.
 public :: ebands_set_fermie       ! Change the fermi level (assume metallic scheme).
 public :: ebands_set_nelect       ! Change the number of electrons (assume metallic scheme).
 public :: ebands_get_muT_with_fd  ! Change the number of electrons (assume metallic scheme).
 public :: ebands_calc_nelect      ! Compute nelect from Fermi level and Temperature.
 public :: ebands_report_gap       ! Print info on the fundamental and direct gap.
 public :: ebands_ncwrite          ! Write object to NETCDF file (use ncid)
 public :: ebands_ncwrite_path     ! Dump the object into NETCDF file (use filepath)
 public :: ebands_write_nesting    ! Calculate the nesting function and output data to file.
 public :: ebands_expandk          ! Build a new ebands_t in the full BZ.
 public :: ebands_downsample       ! Build a new ebands_t with a downsampled IBZ.
 public :: ebands_chop             ! Build a new ebands_t with selected nbands.
 public :: ebands_get_edos         ! Compute e-DOS from band structure.
 public :: ebands_get_jdos         ! Compute electron joint-DOS from band structure.

 public :: ebands_get_edos_matrix_elements ! Compute e-DOS and other DOS-like quantities involving
                                           ! vectorial or tensorial matrix elements.

 public :: ebands_interp_kmesh     ! Use SWK Interpolate energies on a k-mesh.
 public :: ebands_interp_kpath     ! Interpolate energies on a k-path.
 public :: ebands_interpolate_kpath

 public :: ebands_prtbltztrp          ! Output files for BoltzTraP code.
 public :: ebands_prtbltztrp_tau_out  ! Output files for BoltzTraP code,
 public :: ebands_write               ! Driver routine to write bands in different (txt) formats.
!!***

!----------------------------------------------------------------------

!!****t* m_ebands/edos_t
!! NAME
!! edos_t
!!
!! FUNCTION
!! Store the electronic DOS
!!
!! SOURCE

 type,public :: edos_t

   integer :: nsppol
    ! Number of spins.

   integer :: nkibz
    ! Number of k-points in the IBZ.

   integer :: nw
   ! Number of points in the frequency mesh.

   integer :: ief = 0
   ! Rightmost Index of the energy mesh such as IDOS[mesh[ief]] < nelect.
   ! 0 if Fermi level could not be computed
   ! Note the value of gef stored in edos_t is computed by performing
   ! a linear interpolation between ief and ief + 1

   integer :: intmeth
   ! 1 for gaussian, 2 tetra

   real(dp) :: broad = zero
   ! Gaussian broadening

   real(dp) :: step
   ! Step of the mesh

   real(dp),allocatable :: mesh(:)
   ! mesh(nw)

   real(dp),allocatable :: dos(:,:)
   ! dos(nw, 0:nsppol)
   ! Total DOS, spin up and spin down component.

   real(dp),allocatable :: idos(:,:)
   ! idos(nw, 0:nsppol)
   ! Integrated DOS: (total, spin up, spin down) component.

   real(dp),allocatable :: gef(:)
   ! gef(0:nsppol)
   ! DOS at the Fermi level. Total, spin up, spin down

 contains

   procedure :: free => edos_free
   ! Free memory

   procedure :: write => edos_write
   ! Write results to file (formatted mode)

   procedure :: print => edos_print
   ! Print eDOS info to Fortran unit.

   procedure :: ncwrite => edos_ncwrite
   ! Write eDOS to netcdf file.

 end type edos_t
!!***

!!****t* m_ebands/jdos_t
!! NAME
!! jdos_t
!!
!! FUNCTION
!! Store the electron joint DOS
!!
!! SOURCE

 type,public :: jdos_t

   integer :: nsppol
    ! Number of spins.

   integer :: nkibz
    ! Number of k-points in the IBZ.

   integer :: nw
   ! Number of points in the frequency mesh.

   integer :: intmeth
   ! 1 for gaussian, 2 tetra

   real(dp) :: broad = zero
   ! Gaussian broadening

   real(dp) :: step
   ! Step of the mesh

   real(dp),allocatable :: mesh(:)
   ! mesh(nw)

   real(dp),allocatable :: values(:,:)
   ! dos(nw,0:nsppol)
   ! Total jDOS, spin up and spin down component.

 contains

   procedure :: free => jdos_free
   ! Free memory

   !procedure :: write => jdos_write
   ! Write results to file (formatted mode)

   !procedure :: print => jdos_print
   ! Print jDOS info to Fortran unit.

   procedure :: ncwrite => jdos_ncwrite
   ! Write jDOS to netcdf file.

 end type jdos_t
!!***

!----------------------------------------------------------------------

!!****t* m_ebands/gaps_t
!! NAME
!! gaps_t
!!
!! FUNCTION
!! Structure with information on the fundamental and direct gaps returned by ebands_report_gap.
!!
!! TODO
!! Remove gaps_t, move info about CBM, VBM and Fermi energy (i.e. Fermi level for T --> 0) inside
!! ebands_t and make sure that all the setter methods of ebands_t support the new protocol.
!!
!! SOURCE

 type,public :: gaps_t

   integer :: nsppol
    ! Number of spins.

   integer,allocatable :: fo_kpos(:,:)
    ! fo_kpos(3,nsppol)
    ! fo_kpos(1:2,spin) ==> Indices of the k-points where the homo, lumo states are located (for each spin).
    ! fo_kpos(3,spin)   ==> the index of k-point where the direct gap is located (for each spin).

   real(dp) :: fermie
    ! Fermi energy taken from ebands.

   real(dp) :: nelect
    ! Number of electrons taken from ebands.

   integer,allocatable :: ierr(:)
    ! ierr(nsppol
    !   0 if the gap has been computed.
    !   1 if the system (or spin-channel) is metallic.
    !   2 if gaps were not computed (because there are only valence bands).

   real(dp),allocatable :: fo_values(:,:)
     ! fo_values(2,nsppol)]
     ! Fundamental and direct gaps (in Hartree) for each spin.

   real(dp),allocatable :: vb_max(:), cb_min(:)
     ! vb_max(nsppol)
     ! valence band max and conduction band min for each spin in Ha.
     ! Only for Semiconductors, set to (+, -) huge(one) for metals.

   real(dp),pointer :: kpoints(:,:) => null()
     ! Reference to the k-points of the band structure used to compute the gaps.

   character(len=500),allocatable :: errmsg_spin(:)
     ! errmsg_spin(nsppol)
     ! String with human-readable error messages if ierr(spin) != 0.

 contains

   procedure :: free => gaps_free
   ! Free memory

   procedure :: print => gaps_print
   ! Print info on the gaps

 end type gaps_t

 public :: ebands_get_gaps     ! Build the gaps object from a bandstructure.
 public :: ebands_print_gaps   ! Helper function to print gaps directrly from ebands.
!!***

!!****t* m_ebands/klinterp_t
!! NAME
!! klinterp_t
!!
!! FUNCTION
!!  Linear interpolation of eigenvalue-like quantities (scalars with the same symmetry as the KS eigenvalues)
!!  Used, for instance, to interpolate electron or phonon lifetimes.
!!
!! SOURCE

 type,public :: klinterp_t

   integer :: bsize, nsppol, ndat
   ! Max number of bands, number of independent spin polarization, size of "extra" dimension.

   integer :: nkx, nky, nkz
   ! Number of divisions of the grid enclosing the first unit cell

   real(dp),allocatable :: data_uk_bsd(:,:,:,:)
    ! (nkx*nky*nkz, mband, nsppol, ndat)

 contains

   procedure :: free => klinterp_free
    ! Free dynamic memory

   procedure :: eval_bsd => klinterp_eval_bsd
    ! Interpolate values at an arbitrary k-point.

 end type klinterp_t

 public :: klinterp_new         ! Build interpolator.


!----------------------------------------------------------------------

CONTAINS  !=====================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_print_gaps
!! NAME
!! ebands_print_gaps
!!
!! FUNCTION
!!  Helper function to print gaps directrly from ebands.
!!
!! INPUTS
!!  ebands<ebands_t>=Info on the band structure, the smearing technique and the physical temperature used.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_print_gaps(ebands, unit, header)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(in)  :: ebands
 integer,intent(in) :: unit
 character(len=*),optional,intent(in) :: header

!Local variables-------------------------------
!scalars
 integer :: ierr, spin
 type(gaps_t) :: gaps

! *********************************************************************

 if (unit == dev_null) return

 gaps = ebands_get_gaps(ebands, ierr)
 if (ierr /= 0) then
   do spin=1, ebands%nsppol
     write(unit, "(2a)")"WARNING: " // trim(gaps%errmsg_spin(spin))
   end do
 end if

 if (present(header)) then
   call gaps%print(unit=std_out, header=header)
 else
   call gaps%print(unit=std_out, header=header)
 end if
 call gaps%free()

end subroutine ebands_print_gaps
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_get_gaps
!! NAME
!! ebands_get_gaps
!!
!! FUNCTION
!!  Returns a gaps_t object with info on the fundamental and direct gap.
!!
!! INPUTS
!!  ebands<ebands_t>=Info on the band structure, the smearing technique and the physical temperature used.
!!  [kmask]=Logical mask used to exclude k-points.
!!
!! OUTPUT
!!  ierr=Return code (!=0 signals failure)
!!  gaps<gaps_t>=object with info on the gaps (calleris responsible for freeing the object).
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(gaps_t) function ebands_get_gaps(ebands, ierr, kmask) result(gaps)

!Arguments ------------------------------------
!scalars
 class(ebands_t),target,intent(in)  :: ebands
 integer,intent(out) :: ierr
!arrays
 logical,optional,intent(in) :: kmask(ebands%nkpt)

!Local variables-------------------------------
!scalars
 integer :: ikibz,nband_k,spin,nsppol,ikopt,ivk,ick,ivb,icb
 real(dp),parameter :: tol_fermi=tol6
 real(dp) :: fun_gap,opt_gap
 logical :: ismetal
 !type(ebands_t)  :: tmp_ebands
!arrays
 integer :: val_idx(ebands%nkpt, ebands%nsppol)
 real(dp) :: top_valence(ebands%nkpt), bot_conduct(ebands%nkpt)
 logical :: my_kmask(ebands%nkpt)

! *********************************************************************

 nsppol = ebands%nsppol

 ! Initialize gaps_t
 gaps%nsppol = nsppol
 gaps%nelect = ebands%nelect
 ABI_MALLOC(gaps%fo_kpos, (3, nsppol))
 ABI_MALLOC(gaps%ierr, (nsppol))
 ABI_MALLOC(gaps%fo_values, (2, nsppol))
 ABI_MALLOC(gaps%vb_max, (nsppol))
 ABI_MALLOC(gaps%cb_min, (nsppol))
 ABI_MALLOC(gaps%errmsg_spin, (nsppol))
 gaps%kpoints => ebands%kptns

 gaps%fo_kpos = 0
 gaps%ierr = 0
 gaps%fo_values = zero
 gaps%vb_max = huge(one); gaps%cb_min = -huge(one)
 gaps%errmsg_spin(:) = ""
 gaps%fermie = ebands%fermie

 my_kmask=.TRUE.; if (PRESENT(kmask)) my_kmask=kmask

 val_idx(:,:) = ebands_get_valence_idx(ebands, tol_fermi)

 spin_loop: &
&  do spin=1,nsppol

   ! No output if system is metallic
   ismetal = ANY(val_idx(:,spin) /= val_idx(1,spin))
   if (ismetal) then
     gaps%ierr(spin) = 1
     write(gaps%errmsg_spin(spin), "(a,i0)")" Detected metallic system for spin channel: ", spin
     cycle
   endif

   ivb = val_idx(1, spin)
   icb = ivb + 1

   do ikibz=1,ebands%nkpt
     if (.not. my_kmask(ikibz)) cycle
     nband_k = ebands%nband(ikibz + (spin-1)*ebands%nkpt)
     top_valence(ikibz) = ebands%eig(ivb, ikibz, spin)
     if (icb > nband_k) then
       gaps%ierr(spin) = 2
       gaps%errmsg_spin(spin) = "Not enough states to calculate the band gap."
       cycle spin_loop
     end if
     bot_conduct(ikibz) = ebands%eig(icb, ikibz, spin)
   end do

   ! Minimum of the direct Gaps
   ikopt = imin_loc(bot_conduct - top_valence, MASK=my_kmask)
   opt_gap = bot_conduct(ikopt) - top_valence(ikopt)

   ! Fundamental Gap
   ick = imin_loc(bot_conduct, MASK=my_kmask)
   ivk = imax_loc(top_valence, MASK=my_kmask)

   gaps%vb_max(spin) = ebands%eig(ivb, ivk, spin)
   gaps%cb_min(spin) = ebands%eig(icb, ick, spin)
   fun_gap = ebands%eig(icb, ick, spin) - ebands%eig(ivb, ivk, spin)
   gaps%fo_values(:, spin) = [fun_gap, opt_gap]
   gaps%fo_kpos(:, spin) = [ivk, ick, ikopt]
 end do spin_loop

 ierr = maxval(gaps%ierr)

 ! TODO
 !gaps = ebands_get_gaps(ebands, gap_err)
 !if (ierr /= 0) then
 !  ! In case of error try to enforce semiconductor occupations before calling ebands_get_gaps
 !  ! This might still fail though...
 !  call gaps%free()
 !  call ebands_copy(ebands, tmp_ebands)
 !  call ebands_set_scheme(tmp_ebands, occopt3, dtset%tsmear, dtset%spinmagntarget, dtset%prtvol, update_occ=.False.)
 !  call ebands_set_nelect(tmp_ebands, tmp_ebands%nelect-dtset%eph_extrael, dtset%spinmagntarget, msg)
 !  gaps = ebands_get_gaps(tmp_ebands, ierr)
 !  call ebands_free(tmp_ebands)
 !end if

end function ebands_get_gaps
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/gaps_free
!! NAME
!!  gaps_free
!!
!! FUNCTION
!!  Free the memory allocated in gaps_t
!!
!! PARENTS
!!      m_sigmaph,setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine gaps_free(gaps)

!Arguments ------------------------------------
 class(gaps_t),intent(inout) :: gaps

! *********************************************************************

!integer
 ABI_SFREE(gaps%fo_kpos)
 ABI_SFREE(gaps%ierr)

!real
 ABI_SFREE(gaps%fo_values)
 ABI_SFREE(gaps%vb_max)
 ABI_SFREE(gaps%cb_min)

!chars
 ABI_SFREE(gaps%errmsg_spin)

! nullify pointers
 nullify(gaps%kpoints)

end subroutine gaps_free
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/gaps_print
!! NAME
!! gaps_print
!!
!! FUNCTION
!!  Print info on the fundamental and direct gap.
!!
!! INPUTS
!!  gaps<gaps_t>=Object with info on the gaps.
!!  [header]=Optional title.
!!  [unit]=Optional unit for output (std_out if not specified)
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_sigmaph,setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine gaps_print(gaps, unit, header)

!Arguments ------------------------------------
!scalars
 class(gaps_t),intent(in)  :: gaps
 integer,intent(in),optional :: unit
 character(len=*),intent(in),optional :: header

!Local variables-------------------------------
!scalars
 integer :: spin, ikopt, ivk, ick, my_unt
 real(dp) :: fun_gap, opt_gap
 character(len=500) :: msg

! *********************************************************************

 my_unt =std_out; if (present(unit)) my_unt = unit
 if (my_unt == dev_null) return

 do spin=1,gaps%nsppol
   if (spin == 1) then
     msg = ch10
     if (present(header)) msg = ch10//' === '//trim(adjustl(header))//' === '
     call wrtout(my_unt, msg)
   end if

   if (gaps%ierr(spin) /= 0) then
     call wrtout(my_unt, gaps%errmsg_spin(spin))
     continue
   end if

   ! Get minimum of the direct Gap.
   fun_gap = gaps%fo_values(1, spin)
   opt_gap = gaps%fo_values(2, spin)

   if (any(gaps%fo_kpos(:,spin) == 0)) then
     call wrtout(my_unt, sjoin(" Cannot detect gap for spin: ", itoa(spin)))
     cycle
   end if

   ivk = gaps%fo_kpos(1, spin)
   ick = gaps%fo_kpos(2, spin)
   ikopt = gaps%fo_kpos(3, spin)

   write(msg,'(a,i2,a,2(a,f6.2,a,2a),30x,2a)') &
    '  >>>> For spin ', spin, ch10, &
    '   Minimum direct gap = ',opt_gap*Ha_eV,' (eV), located at k-point     : ', trim(ktoa(gaps%kpoints(:,ikopt))),ch10, &
    '   Fundamental gap    = ',fun_gap*Ha_eV,' (eV), Top of valence bands at: ', trim(ktoa(gaps%kpoints(:,ivk))),ch10, &
                                             '       Bottom of conduction at: ', trim(ktoa(gaps%kpoints(:,ick)))
   call wrtout(my_unt, msg)
   write(msg, "((a,f6.2,2a))")"   Valence Max:    ", gaps%vb_max(spin) * Ha_eV, " (eV) at: ", trim(ktoa(gaps%kpoints(:, ivk)))
   call wrtout(my_unt, msg)
   write(msg, "((a,f6.2,2a))")"   Conduction min: ", gaps%cb_min(spin) * Ha_eV, " (eV) at: ", trim(ktoa(gaps%kpoints(:, ick)))
   call wrtout(my_unt, msg)
 end do ! spin

 if (any(gaps%fo_kpos == 0)) then
   write(msg, "((2(a,f6.2)))")  "   Fermi level:", gaps%fermie * Ha_eV, " (eV) with nelect:", gaps%nelect
   call wrtout(my_unt, msg, newlines=1)
 else
   call wrtout(my_unt, "", newlines=1)
 end if

end subroutine gaps_print
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_init
!! NAME
!! ebands_init
!!
!! FUNCTION
!! This subroutine initializes the ebands_t structured datatype
!!
!! INPUTS
!! bantot=total number of bands (=sum(nband(:))
!! doccde(bantot)=derivative of the occupation numbers with respect to the energy (Ha)
!! eig(bantot)=eigenvalues (hartree)
!! istwfk(nkpt)=parameter that describes the storage of wfs.
!! kptns(3,nkpt)=k points in terms of recip primitive translations
!! nband(nkpt*nsppol)=number of bands
!! nelect=Number of electrons.
!! nkpt=number of k points
!! npwarr(nkpt)=number of planewaves at each k point
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nspinor=Number of spinor components
!! occopt=Occupation options (see input variable)
!! occ(bantot)=occupation numbers
!! tphysel=Physical temperature (input variable)
!! tsmear=Temperature of smearing.
!! wtk(nkpt)=weight assigned to each k point
!! charge=Additional charge added to the unit cell (input variable).
!! kptopt=Option for k-point generation (see input variable)
!! kptrlatt_orig=Original value of kptrlatt given in input
!! nshiftk_orig=Original number of shifts given in input
!! shiftk_orig(3,nshiftk_orig)=Original set of shifts given in input
!! kptrlatt=Value of kptrlatt after inkpts
!! nshiftk=Number of shifts after inkpts
!! shiftk(3,nshiftk)=Set of shifts after inkpts.
!!
!! OUTPUT
!! ebands<ebands_t>=the ebands_t datatype
!!
!! SIDE EFFECTS
!!  %entropy and %fermie initialized to zero.
!!
!! PARENTS
!!      dfpt_looppert,eig2tot,gstate,m_ebands,mlwfovlp_qp,optic,outscfcv
!!      setup_bse,setup_bse_interp,setup_screening,setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_init(bantot,ebands,nelect,doccde,eig,istwfk,kptns,&
  nband,nkpt,npwarr,nsppol,nspinor,tphysel,tsmear,occopt,occ,wtk,&
  charge, kptopt, kptrlatt_orig, nshiftk_orig, shiftk_orig, kptrlatt, nshiftk, shiftk)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bantot,nkpt,nsppol,nspinor,occopt
 real(dp),intent(in) :: nelect,tphysel,tsmear
 type(ebands_t),intent(out) :: ebands
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol),npwarr(nkpt)
 real(dp),intent(in) :: doccde(bantot),eig(bantot),kptns(3,nkpt),occ(bantot)
 real(dp),intent(in) :: wtk(nkpt)
 integer,intent(in) :: kptopt, nshiftk_orig, nshiftk
 real(dp),intent(in) :: charge
 integer,intent(in) :: kptrlatt_orig(3,3),kptrlatt(3,3)
 real(dp),intent(in) :: shiftk_orig(3,nshiftk_orig),shiftk(3,nshiftk)

! *************************************************************************

 ! Copy the scalars
 ! MG TODO here there is a inconsistency in the way occ are treated in the header
 ! (only the states used, bantot. are saved, and the way occ. and energies
 ! are passed to routines (mband,nkpt,nsppol). It might happen that bantot<mband*nktp*nsppol
 ! this should not lead to problems since arrays are passed by reference
 ! anyway the treatment of these arrays have to be rationalized
 ebands%bantot =bantot
 ebands%mband  =MAXVAL(nband(1:nkpt*nsppol))
 ebands%nkpt   =nkpt
 ebands%nspinor=nspinor
 ebands%nsppol =nsppol
 ebands%occopt =occopt

 ebands%entropy=zero  ! Initialize results
 ebands%fermie =zero  ! Initialize results
 ebands%nelect =nelect
 ebands%tphysel=tphysel
 ebands%tsmear =tsmear

 ! Allocate the components
 ABI_MALLOC(ebands%nband,(nkpt*nsppol))
 ABI_MALLOC(ebands%istwfk,(nkpt))
 ABI_MALLOC(ebands%npwarr,(nkpt))
 ABI_MALLOC(ebands%kptns,(3,nkpt))

 ! Copy the arrays
 ebands%nband(1:nkpt*nsppol)=nband(1:nkpt*nsppol)
 ebands%istwfk(1:nkpt)      =istwfk(1:nkpt)
 ebands%npwarr(1:nkpt)      =npwarr(1:nkpt)
 ebands%kptns(1:3,1:nkpt)   =kptns(1:3,1:nkpt)

 ! In ebands, energies and occupations are stored in a matrix (mband,nkpt,nsppol).
 ! put_eneocc_vect is used to reshape the values stored in vectorial form.
 ABI_MALLOC(ebands%eig   ,(ebands%mband,nkpt,nsppol))
 ABI_MALLOC(ebands%occ   ,(ebands%mband,nkpt,nsppol))
 ABI_MALLOC(ebands%doccde,(ebands%mband,nkpt,nsppol))

 call put_eneocc_vect(ebands,'eig',   eig   )
 call put_eneocc_vect(ebands,'occ',   occ   )
 call put_eneocc_vect(ebands,'doccde',doccde)

 ABI_MALLOC(ebands%wtk,(nkpt))
 ebands%wtk(1:nkpt)=wtk(1:nkpt)

 ebands%kptopt = kptopt
 ebands%nshiftk_orig = nshiftk_orig
 ebands%nshiftk = nshiftk
 ebands%charge = charge
 ebands%kptrlatt_orig = kptrlatt_orig
 ebands%kptrlatt = kptrlatt

 call alloc_copy(shiftk_orig, ebands%shiftk_orig)
 call alloc_copy(shiftk, ebands%shiftk)

end subroutine ebands_init
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_from_hdr
!! NAME
!! ebands_from_hdr
!!
!! FUNCTION
!! This subroutine initializes the ebands_t datatype from the abinit header by
!! calling the main creation method.
!!
!! INPUTS
!!  Hdr<hdr_type>=Abinit header.
!!  mband=Maximum number of bands.
!!  ene3d(mband,Hdr%nkpt,Hdr%nsppol)=Energies.
!!  [nelect]=Number of electrons per unit cell.
!!    Optional argument that can be used for performing a ridid shift of the fermi level.
!!    in the case of metallic occupancies.
!!    If not specified, nelect will be initialized from Hdr.
!!
!! OUTPUT
!!  ebands<ebands_t>=The ebands_t datatype completely initialized.
!!
!! PARENTS
!!      elphon,eph,m_iowf,m_wfk,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

type(ebands_t) function ebands_from_hdr(hdr, mband, ene3d, nelect) result(ebands)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband
 type(hdr_type),intent(in) :: hdr
 real(dp),optional,intent(in) :: nelect
!arrays
 real(dp),intent(in) :: ene3d(mband,hdr%nkpt,hdr%nsppol)

!Local variables-------------------------------
!scalars
 real(dp) :: my_nelect
!arrays
 real(dp),allocatable :: ugly_doccde(:),ugly_ene(:)

! *************************************************************************

 my_nelect = hdr%nelect; if (present(nelect)) my_nelect = nelect

 ! Have to use ugly 1d vectors to call ebands_init
 ABI_CALLOC(ugly_doccde, (hdr%bantot))
 ABI_MALLOC(ugly_ene,(hdr%bantot))

 call pack_eneocc(hdr%nkpt,hdr%nsppol,mband,hdr%nband,hdr%bantot,ene3d,ugly_ene)

 call ebands_init(hdr%bantot,ebands,my_nelect,ugly_doccde,ugly_ene,hdr%istwfk,hdr%kptns,hdr%nband,hdr%nkpt,&
   hdr%npwarr,hdr%nsppol,hdr%nspinor,hdr%tphysel,hdr%tsmear,hdr%occopt,hdr%occ,hdr%wtk,&
   hdr%charge, hdr%kptopt, hdr%kptrlatt_orig, hdr%nshiftk_orig, hdr%shiftk_orig, hdr%kptrlatt, hdr%nshiftk, hdr%shiftk)

 ! Copy the fermi level reported in the header
 ebands%fermie = hdr%fermie

 ABI_FREE(ugly_doccde)
 ABI_FREE(ugly_ene)

end function ebands_from_hdr
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_from_dtset
!! NAME
!! ebands_from_dtset
!!
!! FUNCTION
!! Build and return a new ebands_t datatype. Dimensions are taken from the abinit dataset.
!!
!! INPUTS
!!  dtset<dataset_type>=Abinit dataset
!!  npwarr(dtset%nkpt)=Number of G-vectors for each k-point.
!!
!! OUTPUT
!!  ebands<ebands_t>=The ebands_t datatype completely initialized.
!!    The Fermi level and the entropy are set to zero.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(ebands_t) function ebands_from_dtset(dtset, npwarr) result(new)

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: npwarr(dtset%nkpt)

!Local variables-------------------------------
!scalars
 integer :: bantot
!arrays
 real(dp),allocatable :: ugly_doccde(:),ugly_ene(:),ugly_occ(:)
! *************************************************************************

 ! Have to use ugly 1d vectors to call ebands_init
 bantot = sum(dtset%nband)
 ABI_CALLOC(ugly_doccde, (bantot))
 ABI_CALLOC(ugly_ene, (bantot))
 ABI_CALLOC(ugly_occ, (bantot))

 call ebands_init(bantot,new,dtset%nelect,ugly_doccde,ugly_ene,dtset%istwfk,dtset%kptns,dtset%nband,dtset%nkpt,&
  npwarr,dtset%nsppol,dtset%nspinor,dtset%tphysel,dtset%tsmear,dtset%occopt,ugly_occ,dtset%wtk,&
  dtset%charge, dtset%kptopt, dtset%kptrlatt_orig, dtset%nshiftk_orig, dtset%shiftk_orig, &
  dtset%kptrlatt, dtset%nshiftk, dtset%shiftk)

 ABI_FREE(ugly_doccde)
 ABI_FREE(ugly_ene)
 ABI_FREE(ugly_occ)

end function ebands_from_dtset
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_free
!! NAME
!! ebands_free
!!
!! FUNCTION
!! Deallocates the components of the ebands_t structured datatype
!!
!! INPUTS
!!  ebands<ebands_t>=The data type to be deallocated.
!!
!! OUTPUT
!!  Deallocate the dynamic arrays in the ebands_t type.
!!
!! PARENTS
!!      bethe_salpeter,dfpt_looppert,eig2tot,elphon,eph,fold2Bloch,gstate
!!      m_ebands,m_exc_spectra,m_haydock,m_ioarr,m_iowf,m_shirley,m_sigmaph
!!      m_wfk,mlwfovlp_qp,nonlinear,optic,outscfcv,respfn,screening,sigma
!!      wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_free(ebands)

!Arguments ------------------------------------
 class(ebands_t),intent(inout) :: ebands
! *************************************************************************

 ABI_SFREE(ebands%istwfk)
 ABI_SFREE(ebands%nband)
 ABI_SFREE(ebands%npwarr)
 ABI_SFREE(ebands%kptns)
 ABI_SFREE(ebands%eig)
 ABI_SFREE(ebands%linewidth)
 ABI_SFREE(ebands%occ)
 ABI_SFREE(ebands%doccde)
 ABI_SFREE(ebands%wtk)
 ABI_SFREE(ebands%velocity)
 ABI_SFREE(ebands%kTmesh)
 ABI_SFREE(ebands%shiftk_orig)
 ABI_SFREE(ebands%shiftk)

end subroutine ebands_free
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_copy
!! NAME
!!  ebands_copy
!!
!! FUNCTION
!! This subroutine performs a deep copy of an ebands_t datatype.
!!
!! INPUTS
!!  ibands<ebands_t>=The data type to be copied.
!!
!! OUTPUT
!!  obands<ebands_t>=The copy.
!!
!! PARENTS
!!      m_exc_spectra,m_haydock,m_sigmaph,optic,screening,setup_bse
!!      setup_bse_interp,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_copy(ibands, obands)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(in)  :: ibands
 class(ebands_t),intent(out) :: obands

! *********************************************************************

 ! Copy scalars
 obands%bantot  = ibands%bantot
 obands%mband   = ibands%mband
 obands%nkpt    = ibands%nkpt
 obands%nspinor = ibands%nspinor
 obands%nsppol  = ibands%nsppol
 obands%occopt  = ibands%occopt
 obands%kptopt  = ibands%kptopt
 obands%nshiftk_orig = ibands%nshiftk_orig
 obands%nshiftk     = ibands%nshiftk

 obands%charge  = ibands%charge
 obands%entropy = ibands%entropy
 obands%fermie  = ibands%fermie
 obands%nelect  = ibands%nelect
 obands%tphysel = ibands%tphysel
 obands%tsmear  = ibands%tsmear

 obands%kptrlatt_orig = ibands%kptrlatt_orig
 obands%kptrlatt = ibands%kptrlatt

 ! Copy allocatable arrays
 ! integer
 call alloc_copy(ibands%istwfk, obands%istwfk)
 call alloc_copy(ibands%nband , obands%nband )
 call alloc_copy(ibands%npwarr, obands%npwarr)

 ! real
 call alloc_copy(ibands%kptns , obands%kptns )
 call alloc_copy(ibands%eig   , obands%eig   )
 call alloc_copy(ibands%occ   , obands%occ   )
 call alloc_copy(ibands%doccde, obands%doccde)
 call alloc_copy(ibands%wtk   , obands%wtk   )
 call alloc_copy(ibands%shiftk_orig, obands%shiftk_orig)
 call alloc_copy(ibands%shiftk, obands%shiftk)

 if (allocated(ibands%linewidth)) call alloc_copy(ibands%linewidth, obands%linewidth)

end subroutine ebands_copy
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_move_alloc
!! NAME
!!  ebands_move_alloc
!!
!! FUNCTION
!!  Transfer allocate from `from_ebands` to `to_ebands`.
!!  `from_ebands` is destroyed when the routine returns.
!!
!! SOURCE

subroutine ebands_move_alloc(from_ebands, to_ebands)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(inout) :: from_ebands
 class(ebands_t),intent(inout) :: to_ebands

! *********************************************************************

 call ebands_free(to_ebands)
 call ebands_copy(from_ebands, to_ebands)
 call ebands_free(from_ebands)

end subroutine ebands_move_alloc
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_print
!! NAME
!! ebands_print
!!
!! FUNCTION
!! Print the content of the object.
!!
!! INPUTS
!!  ebands<ebands_t>The type containing the data.
!!  [unit]=Unit number (std_out if None)
!!  [header]=title for info
!!  [prtvol]=Verbosity level (0 if None)
!!  [mode_paral]=Either 'COLL' or 'PERS' ('COLL' if None).
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      eph,setup_bse,setup_bse_interp,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_print(ebands, header, unit, prtvol, mode_paral)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=*),optional,intent(in) :: header
 character(len=4),optional,intent(in) :: mode_paral
 class(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
 integer :: spin,ikpt,my_unt,my_prtvol,ii
 character(len=4) :: my_mode
 character(len=500) :: msg
! *************************************************************************

 my_unt   =std_out; if (PRESENT(unit      )) my_unt   =unit
 my_prtvol=0      ; if (PRESENT(prtvol    )) my_prtvol=prtvol
 my_mode  ='COLL' ; if (PRESENT(mode_paral)) my_mode  =mode_paral

 msg=' ==== Info on the ebands_t ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(my_unt,msg,my_mode)

 write(msg,'(6(a,i0,a))')&
   '  Number of spinorial components ...... ',ebands%nspinor,ch10,&
   '  Number of spin polarizations ........ ',ebands%nsppol,ch10,&
   '  Number of k-points in the IBZ ....... ',ebands%nkpt,ch10,&
   '  kptopt .............................. ',ebands%kptopt,ch10,&
   '  Maximum number of bands ............. ',ebands%mband,ch10,&
   '  Occupation option ................... ',ebands%occopt,ch10
 call wrtout(my_unt,msg,my_mode)

 !EBANDS_NEW
 !write(msg,"(a)")&
 !  " kptrlatt ..............................",ltoa(ebands%kptrlatt)

 write(msg,'(2(a,f14.2,a),4(a,f14.6,a))')&
   '  Number of valence electrons ......... ',ebands%nelect,ch10,&
   '  Extra charge ........................ ',ebands%charge,ch10,&
   '  Fermi level  ........................ ',ebands%fermie,ch10,&
   '  Entropy ............................. ',ebands%entropy,ch10,&
   '  Tsmear value ........................ ',ebands%tsmear,ch10,&
   '  Tphysel value ....................... ',ebands%tphysel,ch10
 call wrtout(my_unt,msg,my_mode)

 if (my_prtvol > 10) then
   if (ebands%nsppol == 1)then
     call wrtout(my_unt, sjoin(' New occ. numbers for occopt= ', itoa(ebands%occopt),' , spin-unpolarized case.'), my_mode)
   end if

   do spin=1,ebands%nsppol
     if (ebands%nsppol==2) then
       write(msg,'(a,i4,a,i2)')' New occ. numbers for occopt= ',ebands%occopt,' spin ',spin
       call wrtout(my_unt,msg,my_mode)
     end if

     do ikpt=1,ebands%nkpt
       write(msg,'(2a,i4,3a,f6.3,2a)')ch10,&
         ' k-point number ',ikpt,') ',trim(ktoa(ebands%kptns(:,ikpt))),'; weight: ',ebands%wtk(ikpt), ch10, &
         " eig (Ha), eig (eV), occ, doccde"
       call wrtout(my_unt,msg,my_mode)
       do ii=1,ebands%nband(ikpt+(spin-1)*ebands%nkpt)
         write(msg,'(4(f7.3,1x))')ebands%eig(ii,ikpt,spin), ebands%eig(ii,ikpt,spin) * Ha_eV, &
             ebands%occ(ii,ikpt,spin), ebands%doccde(ii,ikpt,spin)
         call wrtout(my_unt,msg,my_mode)
       end do
     end do !ikpt

   end do !spin

 end if !my_prtvol

end subroutine ebands_print
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/unpack_eneocc
!! NAME
!! unpack_eneocc
!!
!! FUNCTION
!!  Helper function to do a reshape of (energies|occupancies|derivate of occupancies)
!!  initially stored in a vector. Return a 3D array index by (band,ikpt,spin)
!!
!! INPUTS
!!  nkpt=number of k-points
!!  nsppol=number of spin polarizations
!!  mband=Max number of bands over k-points (just to dimension the output)
!!  nbands(nkpt*nsppol)=Number of bands at eack k and spin
!!  vect(:)=The input values to reshape
!!  [val]=Optional value used to initialize the array.
!!
!! OUTPUT
!!  array3d(mband,nkpt,nsppol)=Arrays containing the values of vect.
!!   Note that the first dimension is usually larger than the
!!   number of bands really used for a particular k-point and spin.
!!
!! PARENTS
!!      cchi0q0_intraband,m_ebands,m_ioarr,m_iowf
!!
!! CHILDREN
!!
!! SOURCE

subroutine unpack_eneocc(nkpt,nsppol,mband,nband,vect,array3d,val)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nsppol,mband
 real(dp),optional,intent(in) :: val
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: vect(:)
 real(dp),intent(out) :: array3d(mband,nkpt,nsppol)

!Local variables-------------------------------
 integer :: spin,ikpt,band,idx
! *************************************************************************

 if (present(val)) then
   array3d = val
 else
   array3d = huge(one)
 end if

 idx=0
 ! elements in vect are packed in the first positions.
 do spin=1,nsppol
   do ikpt=1,nkpt
     do band=1,nband(ikpt+(spin-1)*nkpt)
      idx=idx+1
      array3d(band,ikpt,spin)=vect(idx)
     end do
   end do
 end do

end subroutine unpack_eneocc
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/pack_eneocc
!! NAME
!! pack_eneocc
!!
!! FUNCTION
!!  Helper function to do a reshape of (energies|occupancies|derivate of occupancies)
!!  initially stored in a 3D arrays returning a vector.
!!
!! INPUTS
!!  nkpt=number of k-points
!!  nsppol=number of spin polarizations
!!  mband=Max number of bands over k-points (just to dimension the output)
!!  nbands(nkpt*nsppol)=Number of bands at eack k and spin
!!  bantot=Total number of bands
!!  array3d(mband,nkpt,nsppol)=Arrays containing the values to reshape.
!!
!! OUTPUT
!!  vect(bantot)=The input values stored in vector mode. Only the values really
!!   considered at each k-point and spin are copied.
!!
!! PARENTS
!!      cchi0q0_intraband,m_ebands,m_shirley
!!
!! CHILDREN
!!
!! SOURCE

subroutine pack_eneocc(nkpt,nsppol,mband,nband,bantot,array3d,vect)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nsppol,mband,bantot
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: array3d(mband,nkpt,nsppol)
 real(dp),intent(out) :: vect(bantot)

!Local variables-------------------------------
 integer :: spin,ikpt,band,idx

! *************************************************************************

 vect(:)=zero
 idx=0
 do spin=1,nsppol
   do ikpt=1,nkpt
     do band=1,nband(ikpt+(spin-1)*nkpt)
       idx=idx+1
       vect(idx)=array3d(band,ikpt,spin)
     end do
   end do
 end do

end subroutine pack_eneocc
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_eneocc_vect
!! NAME
!! get_eneocc_vect
!!
!! FUNCTION
!!  Retrieve energies or occupations from a ebands_t structure accessing by name.
!!  Results are reported in a vector to facilitate the interface with other abinit routines.
!!
!! INPUTS
!!  ebands<ebands_t>The type containing the data.
!!  arr_name=The name of the quantity to retrieve. Allowed values are
!!   == "eig"    == For the eigenvalues.
!!   == "occ"    == For the occupation numbers.
!!   == "doccde" == For the derivative of the occupancies wrt the energy.
!!
!! OUTPUT
!!  vect(ebands%bantot)=The values required.
!!
!! PARENTS
!!      m_ebands
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_eneocc_vect(ebands,arr_name,vect)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: arr_name
 class(ebands_t),intent(in) :: ebands
 real(dp),intent(out) :: vect(ebands%bantot)

!Local variables-------------------------------
 integer :: nkpt,nsppol,mband,bantot
! *************************************************************************

 mband =ebands%mband; bantot=ebands%bantot; nkpt=ebands%nkpt; nsppol=ebands%nsppol

 select case (arr_name)
 case ('occ')
   call pack_eneocc(nkpt,nsppol,mband,ebands%nband,bantot,ebands%occ,vect)
 case ('eig')
   call pack_eneocc(nkpt,nsppol,mband,ebands%nband,bantot,ebands%eig,vect)
 case ('doccde')
   call pack_eneocc(nkpt,nsppol,mband,ebands%nband,bantot,ebands%doccde,vect)
 case default
   MSG_BUG(sjoin('Wrong arr_name:', arr_name))
 end select

end subroutine get_eneocc_vect
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/put_eneocc_vect
!! NAME
!! put_eneocc_vect
!!
!! FUNCTION
!!  Update the energies or the occupations stored in a ebands_t structure.
!!  The input values are stored in a vector according to the abinit convention
!!  In the data type, on the contrary,  we use 3D arrays (mband,nkpt,nsspol)
!!  which are much easier to use inside loops.
!!
!! INPUTS
!!  vect(ebands%bantot)=The new values to be stored in the structure.
!!  arr_name=The name of the quantity to be saved (CASE insensitive).
!!  Allowed values are
!!   == "eig"    == For the eigenvalues.
!!   == "occ"    == For the occupation numbers.
!!   == "doccde" == For the derivative of the occupancies wrt the energy.
!!
!! OUTPUT
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  ebands<ebands_t>=The object with updated values depending on the value of arr_name
!!
!! PARENTS
!!      dfpt_looppert,m_ebands
!!
!! CHILDREN
!!
!! SOURCE

subroutine put_eneocc_vect(ebands,arr_name,vect)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: arr_name
 class(ebands_t),intent(inout) :: ebands
 real(dp),intent(in) :: vect(ebands%bantot)

!Local variables-------------------------------
 integer :: nkpt,nsppol,mband,bantot
! *************************************************************************

 mband =ebands%mband; bantot=ebands%bantot; nkpt  =ebands%nkpt; nsppol=ebands%nsppol

 select case (tolower(arr_name))
 case ('occ')
   call unpack_eneocc(nkpt,nsppol,mband,ebands%nband,vect,ebands%occ, val=zero)
 case ('eig')
   ! DFPT routines call ebands_init with the wrong bantot. Using maxval(vect) causes SIGFAULT
   ! so I have to recompute the correct bantot here
   !ABI_CHECK(sum(ebands%nband) == ebands%bantot, "bantot and nband are incosistent")
   call unpack_eneocc(nkpt,nsppol,mband,ebands%nband,vect,ebands%eig, val=maxval(vect(1:sum(ebands%nband))))
 case ('doccde')
   call unpack_eneocc(nkpt,nsppol,mband,ebands%nband,vect,ebands%doccde, val=zero)
 case default
   MSG_BUG(sjoin('Wrong arr_name= ', arr_name))
 end select

end subroutine put_eneocc_vect
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_get_bandenergy
!! NAME
!! ebands_get_bandenergy
!!
!! FUNCTION
!!  Return the band energy (weighted sum of occupied eigenvalues)
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!! TODO Likely this expression is not accurate since it is not variatonal
!!  One should use
!!   band_energy = \int e N(e) de   for e<Ef , where N(e) is the e-DOS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure real(dp) function ebands_get_bandenergy(ebands) result(band_energy)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
 integer :: spin,ikibz,nband_k
 real(dp) :: wtk
! *********************************************************************

 band_energy=zero
 do spin=1,ebands%nsppol
   do ikibz=1,ebands%nkpt
     wtk=ebands%wtk(ikibz)
     nband_k=ebands%nband(ikibz+(spin-1)*ebands%nkpt)
     band_energy = band_energy + wtk*SUM( ebands%eig(1:nband_k,ikibz,spin)*ebands%occ(1:nband_k,ikibz,spin) )
   end do
 end do

end function ebands_get_bandenergy
!!***

!!****f* m_ebands/ebands_get_valence_idx
!! NAME
!!  ebands_get_valence_idx
!!
!! FUNCTION
!!  For each k-point and spin polarisation, report:
!!
!!    1) the index of the valence in case of semiconductors at T = 0
!!    2) (band_k - 1) where band_k is the first band whose energy is > Fermi energy + told_fermi
!!
!!  using the value of the Fermi level.
!!
!! INPUTS
!!  ebands<ebands_t>=The object describing the band structure.
!!  tol_fermi[optional]
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function ebands_get_valence_idx(ebands, tol_fermi) result(val_idx)

!Arguments ------------------------------------
!scalars
 real(dp),optional,intent(in) :: tol_fermi
 class(ebands_t),intent(in) :: ebands
!arrays
 integer :: val_idx(ebands%nkpt,ebands%nsppol)

!Local variables-------------------------------
 integer :: band,ikpt,spin,idx,nband_k
 real(dp) :: tol_

! *************************************************************************

 tol_=tol6; if (present(tol_fermi)) tol_ = tol_fermi

 do spin=1,ebands%nsppol
   do ikpt=1,ebands%nkpt
     nband_k = ebands%nband(ikpt+(spin-1)*ebands%nkpt)
     idx = 0
     do band=1,nband_k
       if (ebands%eig(band,ikpt,spin) > ebands%fermie + abs(tol_)) then
         idx=band; exit
       end if
     end do
     val_idx(ikpt,spin) = idx - 1
     if (idx==1) val_idx(ikpt,spin) = idx
     if (idx==0) val_idx(ikpt,spin) = nband_k
   end do
 end do

end function ebands_get_valence_idx
!!***

!!****f* m_ebands/ebands_get_bands_from_erange
!! NAME
!!  ebands_get_bands_from_erange
!!
!! FUNCTION
!! Return the indices of the min and max band index within an energy window.
!!
!! INPUTS
!!  elow, ehigh: Min and max energy
!!
!! OUTPUT
!!  bstart, bstop: Min and max band index. Initialized to bstart = huge(1); bstop = -huge(1)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine ebands_get_bands_from_erange(ebands, elow, ehigh, bstart, bstop)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(in) :: ebands
 real(dp),intent(in) :: elow, ehigh
 integer,intent(out) :: bstart, bstop

!Local variables-------------------------------
 integer :: band, ik, spin

! *************************************************************************

 bstart = huge(1); bstop = -huge(1)
 do spin=1,ebands%nsppol
   do ik=1,ebands%nkpt
     do band=1,ebands%nband(ik+(spin-1)*ebands%nkpt)
       if (ebands%eig(band, ik , spin) >= elow .and. ebands%eig(band, ik , spin) <= ehigh) then
          bstart = min(bstart, band)
          bstop = max(bstop, band)
       end if
     end do
   end do
 end do

end subroutine ebands_get_bands_from_erange
!!***

!!****f* m_ebands/ebands_vcbm_range_from_gaps
!! NAME
!!  ebands_vcbm_range_from_gaps
!!
!! FUNCTION
!! Find band and energy range for states close to the CBM/VBM given input energies in ebands and gaps.
!! Return exit status and error message in msg.
!!
!! INPUTS
!!  gaps<gaps_t>=Object with info on the gaps.
!!  erange(2)=Energy range for holes and electrons. Only those states whose relative position
!!    wrt to the VBM/CBM is <= than erange are icluded. Note that relative positions are always
!!    positive (even for holes). Use a negative value to exclude either holes or electrons.
!!
!! OUTPUT
!!  e_lowhigh(2)=min and Max energy.
!!  band_lowhigh=min and Max band index.
!!  [ks_range]: For each spin and k-point, the min and max band index included in the output set.
!!     if (ik, spin) is not included then ib_work(1, ik, spin) > ib_work(2, ik, spin) = -huge(1)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function ebands_vcbm_range_from_gaps(ebands, gaps, erange, e_lowhigh, band_lowhigh, ks_range, msg) result(ierr)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(in) :: ebands
 class(gaps_t),intent(in) :: gaps
 real(dp),intent(in) :: erange(2)
 real(dp),intent(out) :: e_lowhigh(2)
 integer,intent(out) :: band_lowhigh(2)
 integer,optional,intent(out) :: ks_range(2, ebands%nkpt, ebands%nsppol)
 character(len=*),intent(out) :: msg

!Local variables-------------------------------
 integer :: band, ik, spin, band_low, band_high
 real(dp) :: cmin, vmax, ee, elow, ehigh
 integer,allocatable :: ib_work(:,:,:)

! *************************************************************************

 ABI_MALLOC(ib_work, (2, ebands%nkpt, ebands%nsppol))
 elow = huge(one); ehigh = -huge(one)
 band_low = huge(1); band_high = -huge(1)

 ierr = 1
 do spin=1,ebands%nsppol
   ! Get cmb and vbm with some tolerance
   vmax = gaps%vb_max(spin) + tol2 * eV_Ha
   cmin = gaps%cb_min(spin) - tol2 * eV_Ha
   do ik=1,ebands%nkpt
     ib_work(1, ik, spin) = huge(1)
     ib_work(2, ik, spin) = -huge(1)
     do band=1,ebands%nband(ik+(spin-1)*ebands%nkpt)
        ee = ebands%eig(band, ik, spin)
        if (erange(1) > zero) then
          if (ee <= vmax .and. vmax - ee <= erange(1)) then
            ib_work(1, ik, spin) = min(ib_work(1, ik, spin), band)
            ib_work(2, ik, spin) = max(ib_work(2, ik, spin), band)
            elow = min(elow, ee); ehigh = max(ehigh, ee)
            band_low = min(band_low, band); band_high = max(band_high, band)
            !write(std_out, *), "Adding valence", band
          end if
        end if
        if (erange(2) > zero) then
          if (ee >= cmin .and. ee - cmin <= erange(2)) then
            ib_work(1, ik, spin) = min(ib_work(1, ik, spin), band)
            ib_work(2, ik, spin) = max(ib_work(2, ik, spin), band)
            elow = min(elow, ee); ehigh = max(ehigh, ee)
            band_low = min(band_low, band); band_high = max(band_high, band)
            !write(std_out, *)"Adding conduction", band
          end if
        end if
     end do
   end do
 end do

 e_lowhigh = [elow, ehigh]
 band_lowhigh = [band_low, band_high]

 if (present(ks_range)) ks_range = ib_work
 ABI_FREE(ib_work)

 ! Set exit status and msg. Caller will handle it.
 ierr = 0; msg = ""
 if (elow > ehigh) then
   ierr = 1
   write(msg, *)"Cannot find states close to the band edges with erange: ", erange
 end if

end function ebands_vcbm_range_from_gaps
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_apply_scissors
!! NAME
!!  ebands_apply_scissors
!!
!! FUNCTION
!!  Apply a scissor operator of amplitude scissor_energy.
!!
!! INPUTS
!!  scissor_energy=The energy shift
!!
!! OUTPUT
!!
!! SIDE EFFECT
!!  ebands<ebands_t>=The following quantities are modified:
!!   %eig(mband,nkpt,nsppol)=The band structure after the application of the scissor operator
!!   %fermi_energy
!!
!! PARENTS
!!      screening,setup_bse,setup_bse_interp
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_apply_scissors(ebands, scissor_energy)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: scissor_energy
 class(ebands_t),intent(inout) :: ebands

!Local variables-------------------------------
 integer :: ikpt,spin,ival,nband_k
 real(dp) :: spinmagntarget_
 character(len=500) :: msg
!arrays
 integer :: val_idx(ebands%nkpt,ebands%nsppol)
! *************************************************************************

 ! Get the valence band index for each k and spin
 val_idx(:,:) = ebands_get_valence_idx(ebands)

 do spin=1,ebands%nsppol
   if (any(val_idx(:, spin) /= val_idx(1, spin))) then
     write(msg,'(a,i0,a)')&
      'Trying to apply a scissor operator on a metallic band structure for spin: ',spin,&
      'Assuming you know what you are doing, continuing anyway! '
     MSG_COMMENT(msg)
     !Likely newocc will stop, unless the system is semimetallic ?
   end if
 end do

 ! Apply the scissor
 do spin=1,ebands%nsppol
   do ikpt=1,ebands%nkpt
     nband_k = ebands%nband(ikpt+(spin-1)*ebands%nkpt)
     ival = val_idx(ikpt,spin)

     if (nband_k >= ival+1) then
       ebands%eig(ival+1:,ikpt,spin) = ebands%eig(ival+1:,ikpt,spin)+scissor_energy
     else
       write(msg,'(2a,4(a,i0))')&
        'Not enough bands to apply the scissor operator. ',ch10,&
        'spin= ',spin,' ikpt= ',ikpt,' nband_k= ',nband_k,' but valence index= ',ival
       MSG_COMMENT(msg)
     end if

   end do
 end do

 ! Recalculate the Fermi level and occupation factors.
 ! For Semiconductors only the Fermi level is changed (in the middle of the new gap)
 spinmagntarget_ = -99.99_dp !?; if (PRESENT(spinmagntarget)) spinmagntarget_=spinmagntarget
 call ebands_update_occ(ebands, spinmagntarget_)

end subroutine ebands_apply_scissors
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_get_occupied
!! NAME
!!  ebands_get_occupied
!!
!! FUNCTION
!!  For each k-point and spin polarisation, report the band index
!!  after which the occupation numbers are less than tol_occ.
!!
!! INPUTS
!!  ebands<ebands_t>=The object describing the band structure.
!!  tol_occ[Optional]=Tollerance on the occupation factors.
!!
!! OUTPUT
!!
!! NOTES
!!  We assume that the occupation factors monotonically decrease as a function of energy.
!!  This is not always true for every smearing technique implemented in Abinit.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function ebands_get_occupied(ebands, tol_occ) result(occ_idx)

!Arguments ------------------------------------
!scalars
 real(dp),optional,intent(in) :: tol_occ
 class(ebands_t),intent(in) :: ebands
!arrays
 integer :: occ_idx(ebands%nkpt,ebands%nsppol)

!Local variables-------------------------------
 integer :: band,ikpt,spin,idx,nband_k
 real(dp) :: tol_

! *************************************************************************

 tol_=tol8; if (PRESENT(tol_occ)) tol_=tol_occ

 do spin=1,ebands%nsppol
   do ikpt=1,ebands%nkpt
     nband_k = ebands%nband(ikpt+(spin-1)*ebands%nkpt)

     idx=0
     do band=1,nband_k
       if (ebands%occ(band,ikpt,spin) < ABS(tol_)) then
         idx=band; EXIT
       end if
     end do
     occ_idx(ikpt,spin)=idx-1
     if (idx==1) occ_idx(ikpt,spin)=idx
     if (idx==0) occ_idx(ikpt,spin)=nband_k

   end do
 end do

end function ebands_get_occupied
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_enclose_degbands
!! NAME
!!  ebands_enclose_degbands
!!
!! FUNCTION
!!  Adjust ibmin and ibmax such that all the degenerate states are enclosed
!!  between ibmin and ibmax. The routine works for a given k-point a spin.
!!
!! INPUTS
!!  ebands<ebands_t>=The object describing the band structure.
!!  ikibz=Index of the k-point.
!!  spin=Spin index.
!!  tol_enedif=Tolerance on the energy difference.
!!
!! OUTPUT
!!  changed=.TRUE. if ibmin or ibmax has been changed.
!!  [degblock(2,ndeg)]=Table allocated by the routine containing the index
!!    of the bands in the `ndeg` degenerate sub-sets
!!    degblock(1, ii) = first band index in the ii-th degenerate subset.
!!    degblock(2, ii) = last band index in the ii-th degenerate subset.
!!
!! SIDE EFFECTS
!!  ibmin,ibmax=
!!    Input: initial guess for the indeces
!!    Output: All the denerate states are between ibmin and ibmax
!!
!! PARENTS
!!      m_sigmaph,setup_sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_enclose_degbands(ebands, ikibz, spin, ibmin, ibmax, changed, tol_enedif, degblock)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(in) :: ebands
 integer,intent(in) :: ikibz,spin
 integer,intent(inout) :: ibmin,ibmax
 real(dp),intent(in) :: tol_enedif
 logical,intent(out) :: changed
!arrays
 integer,allocatable,optional,intent(out) :: degblock(:,:)

!Local variables-------------------------------
!scalars
 integer :: ib,ibmin_bkp,ibmax_bkp,ndeg
 real(dp) :: emin,emax

! *************************************************************************

 ibmin_bkp = ibmin; ibmax_bkp = ibmax

 emin = ebands%eig(ibmin,ikibz,spin)
 do ib=ibmin-1,1,-1
   if (ABS(ebands%eig(ib,ikibz,spin) - emin) > tol_enedif) then
     ibmin = ib +1
     EXIT
   else
     ibmin = ib
   end if
 end do

 emax = ebands%eig(ibmax,ikibz,spin)
 do ib=ibmax+1,ebands%nband(ikibz+(spin-1)*ebands%nkpt)
   if ( ABS(ebands%eig(ib,ikibz,spin) - emax) > tol_enedif) then
     ibmax = ib - 1
     EXIT
   else
     ibmax = ib
   end if
 end do

 changed = (ibmin /= ibmin_bkp) .or. (ibmax /= ibmax_bkp)

 ! Compute degeneracy table.
 if (present(degblock)) then
   ! Count number of degeneracies.
   ndeg = 1
   do ib=ibmin+1,ibmax
     if ( abs(ebands%eig(ib,ikibz,spin) - ebands%eig(ib-1,ikibz,spin) ) > tol_enedif) ndeg = ndeg + 1
   end do
   ! Build degblock table.
   ABI_REMALLOC(degblock, (2, ndeg))
   ndeg = 1; degblock(1, 1) = ibmin
   do ib=ibmin+1,ibmax
     if ( abs(ebands%eig(ib,ikibz,spin) - ebands%eig(ib-1,ikibz,spin) ) > tol_enedif) then
       degblock(2, ndeg) = ib - 1
       ndeg = ndeg + 1
       degblock(1, ndeg) = ib
     end if
   end do
   degblock(2, ndeg) = ibmax
 end if

end subroutine ebands_enclose_degbands
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_get_erange
!! NAME
!!  ebands_get_erange
!!
!! FUNCTION
!!  Compute the minimum and maximum energy enclosing a list of states
!!  specified by k-points and band indices.
!!
!! INPUTS
!!  ebands<ebands_t>=The object describing the band structure.
!!  nkpts=Number of k-points
!!  kpoints(3,nkpts)=K-points
!!  band_block(2,nkpts)=Gives for each k-points, the initial and the final band index to include.
!!
!! OUTPUT
!!  emin,emax=min and max energy
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_get_erange(ebands, nkpts, kpoints, band_block, emin, emax)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpts
 real(dp),intent(out) :: emin,emax
 class(ebands_t),intent(in) :: ebands
!arrays
 integer,intent(in) :: band_block(2,nkpts)
 real(dp),intent(in) :: kpoints(3,nkpts)

!Local variables-------------------------------
!scalars
 integer :: spin,ik,ikpt,cnt
 type(krank_t) :: krank

! *************************************************************************

 krank = krank_new(ebands%nkpt, ebands%kptns)
 emin = huge(one); emax = -huge(one); cnt = 0

 do spin=1,ebands%nsppol
   do ik=1,nkpts
     ikpt = krank%get_index(kpoints(:,ik))
     if (ikpt == -1) then
       MSG_WARNING(sjoin("Cannot find k-point:", ktoa(kpoints(:,ik))))
       cycle
     end if
     if (.not. (band_block(1,ik) >= 1 .and. band_block(2,ik) <= ebands%mband)) cycle
     cnt = cnt + 1
     emin = min(emin, minval(ebands%eig(band_block(1,ik):band_block(2,ik), ikpt, spin)))
     emax = max(emax, maxval(ebands%eig(band_block(1,ik):band_block(2,ik), ikpt, spin)))
   end do
 end do

 call krank%free()

 ! This can happen if wrong input.
 if (cnt == 0) then
    MSG_WARNING("None of the k-points/bands provided was found in ebands%")
    emin = minval(ebands%eig); emax = maxval(ebands%eig)
 end if

end subroutine ebands_get_erange
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_nelect_per_spin
!! NAME
!!  ebands_nelect_per_spin
!!
!! FUNCTION
!!   Return number of electrons in each spin channel (computed from occoputation factors if nsppol=2)
!!
!! INPUTS
!!  ebands<ebands_t>=The object describing the band structure.
!!
!! OUTPUT
!!  nelect_per_spin(ebands%nsppol)=For each spin the number of electrons (eventually fractional)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function ebands_nelect_per_spin(ebands) result(nelect_per_spin)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(in) :: ebands
!arrays
 real(dp) :: nelect_per_spin(ebands%nsppol)

!Local variables-------------------------------
!scalars
 integer :: iband,ikpt,spin

! *************************************************************************

 nelect_per_spin = ebands%nelect
 if (ebands%nsppol > 1) then
   nelect_per_spin = zero
   do spin=1,ebands%nsppol
     do ikpt=1,ebands%nkpt
       do iband=1,ebands%nband(ikpt+ebands%nkpt*(spin-1))
         nelect_per_spin(spin) = nelect_per_spin(spin) + ebands%wtk(ikpt)*ebands%occ(iband, ikpt, spin)
       end do
     end do
   end do
 end if

end function ebands_nelect_per_spin
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_get_minmax
!! NAME
!!  ebands_get_minmax
!!
!! FUNCTION
!!  Report the min and max value over k-points and bands of (eig|occ|doccde) for each
!!  spin. Cannot use F90 array syntax due to the internal storage used in abinit.
!!
!! INPUTS
!!  ebands<ebands_t>=The object describing the band structure.
!!  arr_name=The name of the array whose min and Max value has to be calculated.
!!   Possible values: 'occ', 'eig' 'doccde'
!!
!! OUTPUT
!! minmax(2,ebands%nsppol)=For each spin the min and max value of the quantity specified by "arr_name"
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function ebands_get_minmax(ebands, arr_name) result(minmax)

!Arguments ------------------------------------
!scalars
 class(ebands_t),target,intent(in) :: ebands
 character(len=*),intent(in) :: arr_name
!arrays
 real(dp) :: minmax(2,ebands%nsppol)

!Local variables-------------------------------
!scalars
 integer :: band,ikpt,spin,nband_k
 real(dp) :: datum
!arrays
 real(dp), ABI_CONTIGUOUS pointer :: rdata(:,:,:)

! *************************************************************************

 select case (tolower(arr_name))
 case ('occ')
   rdata => ebands%occ
 case ('eig')
   rdata => ebands%eig
 case ('doccde')
   rdata => ebands%doccde
 case default
   MSG_BUG(sjoin('Wrong arr_name:', arr_name))
 end select

 minmax(1,:)=greatest_real
 minmax(2,:)=smallest_real

 do spin=1,ebands%nsppol
   do ikpt=1,ebands%nkpt
     nband_k=ebands%nband(ikpt+(spin-1)*ebands%nkpt)
     do band=1,nband_k
       datum=rdata(band,ikpt,spin)
       minmax(1,spin)=MIN(minmax(1,spin),datum)
       minmax(2,spin)=MAX(minmax(2,spin),datum)
     end do
   end do
 end do

end function ebands_get_minmax
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_has_metal_scheme
!! NAME
!! ebands_metallic_scheme
!!
!! FUNCTION
!! Returns .TRUE. if metallic occupation scheme is used.
!! Note that this does not imply that the system is metallic.
!!
!! INPUTS
!! ebands<ebands_t>=The ebands_t datatype
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure logical function ebands_has_metal_scheme(ebands) result(ans)

!Arguments ------------------------------------
 class(ebands_t),intent(in) :: ebands

! *************************************************************************

 ans = (any(ebands%occopt == [3, 4, 5, 6, 7, 8]))

end function ebands_has_metal_scheme
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_write_bxsf
!! NAME
!!  ebands_write_bxsf
!!
!! FUNCTION
!!  Write 3D energies for Fermi surface visualization (XSF format)
!!
!! INPUTS
!!  ebands<ebands_t>=The object describing the band structure.
!!  crystal<crystal_t>=Info on unit cell and symmetries.
!!  fname=File name for output.
!!
!! OUTPUT
!!  ierr=Status error.
!!
!! SIDE EFFECTS
!!  Produce BXSF file.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function ebands_write_bxsf(ebands, crystal, fname) result(ierr)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: fname
 class(ebands_t),intent(in) :: ebands
 class(crystal_t),intent(in) :: crystal

!Local variables-------------------------------
 logical :: use_timrev

! *************************************************************************

 use_timrev = (crystal%timrev==2)

 call printbxsf(ebands%eig,zero,ebands%fermie,crystal%gprimd,ebands%kptrlatt,ebands%mband,&
   ebands%nkpt,ebands%kptns,crystal%nsym,crystal%use_antiferro,crystal%symrec,crystal%symafm,&
   use_timrev,ebands%nsppol,ebands%shiftk,ebands%nshiftk,fname,ierr)

end function ebands_write_bxsf
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_update_occ
!! NAME
!! ebands_update_occ
!!
!! FUNCTION
!! Calculate new occupation numbers, the Fermi level and the Max occupied band index
!! for each spin channel starting from the the knowledge of eigenvalues.
!!
!! INPUTS
!!  spinmagntarget=if differ from -99.99d0, fix the spin polarization (in Bohr magneton)
!!  [stmbias]=
!!  [prtvol]=Verbosity level (0 for lowest level)
!!  ebands<ebands_t>=Info on the band structure, the smearing technique and the physical temperature used.
!!
!! OUTPUT
!!  see also SIDE EFFECTS.
!!
!! SIDE EFFECTS
!!  === For metallic occupation the following quantites are recalculated ===
!!   %fermie=the new Fermi energy
!!   %entropy=the new entropy associated with the smearing.
!!   %occ(mband,nkpt,nsppol)=occupation numbers
!!   %doccde(mband,nkpt,nsppol)=derivative of occupancies wrt the energy for each band and k point
!!  === In case of semiconductors ===
!!   All the quantitities in ebands are left unchanged with the exception of:
!!   %fermie=Redefined so that it is in the middle of the gap
!!   %entropy=Set to zero
!!
!! PARENTS
!!      bethe_salpeter,elphon,eph,get_nv_fs_temp,get_tau_k,m_ebands,optic
!!      screening,setup_bse,setup_bse_interp,setup_sigma,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_update_occ(ebands, spinmagntarget, stmbias, prtvol)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(inout) :: ebands
 integer,optional,intent(in) :: prtvol
 real(dp),intent(in) :: spinmagntarget
 real(dp),optional,intent(in) :: stmbias

!Local variables-------------------------------
!scalars
 integer :: band,mband,ikibz,nkpt,spin,nsppol,my_prtvol,nband_k
 real(dp) :: entropy,fermie,stmbias_local,ndiff,cbot,vtop,maxocc
 character(len=500) :: msg
!arrays
 real(dp) :: nelect_spin(ebands%nsppol),condbottom(ebands%nsppol),valencetop(ebands%nsppol)
 real(dp),allocatable :: doccde(:),occ(:),eigen(:)

! *************************************************************************

 my_prtvol=0; if (PRESENT(prtvol )) my_prtvol=prtvol
 stmbias_local=zero; if (PRESENT(stmbias)) stmbias_local=stmbias

 if (ebands_has_metal_scheme(ebands)) then
   ! Compute new occupation numbers if metallic occupation.
   if (my_prtvol > 10) then
     call wrtout(std_out, sjoin(' metallic scheme: calling newocc with spinmagntarget:', ftoa(spinmagntarget, fmt="f9.5")))
   end if

   ! newocc assumes eigenvalues and occupations packed in 1d-vector!!
   mband = ebands%mband; nkpt = ebands%nkpt; nsppol = ebands%nsppol

   ABI_MALLOC(eigen, (mband*nkpt*nsppol))
   call get_eneocc_vect(ebands, 'eig', eigen)

   ABI_MALLOC(occ, (mband*nkpt*nsppol))
   ABI_MALLOC(doccde, (mband*nkpt*nsppol))

   call newocc(doccde,eigen,entropy,fermie,spinmagntarget,mband,ebands%nband,&
    ebands%nelect,ebands%nkpt,ebands%nspinor,ebands%nsppol,occ,ebands%occopt,&
    my_prtvol,stmbias_local,ebands%tphysel,ebands%tsmear,ebands%wtk)

   ! Save output in ebands%.
   ebands%entropy = entropy
   ebands%fermie  = fermie
   call put_eneocc_vect(ebands, 'occ', occ)
   call put_eneocc_vect(ebands, 'doccde', doccde)
   ABI_FREE(eigen)
   ABI_FREE(occ)
   ABI_FREE(doccde)

 else
   ! Semiconductor (non magnetic case)
   maxocc = two / (ebands%nsppol*ebands%nspinor)
   !
   ! FIXME here there is an inconsistency btw GW and Abinit
   ! In abinit Fermi is set to HOMO while in GW fermi is in the middle
   ! of Gap. In case of crystal systems, the later convention should be preferable.
   ! Anyway we have to decide and follow a unique convention to avoid problems.
   !
   ! Occupation factors MUST be initialized
   !if (ALL(ABS(ebands%occ) < tol6)) then
   if (ebands%occopt /= 2) then
     ABI_CHECK(ebands%nelect == nint(ebands%nelect), "nelect should be integer")
     mband = nint((ebands%nelect * ebands%nspinor) / 2)
     ebands%occ = zero
     ebands%occ(1:mband,:,:) = maxocc
     !MSG_ERROR("Occupation factors are not initialized, likely due to scf = -2")
   end if

   ! Calculate the valence index for each spin channel.
   do spin=1,ebands%nsppol
     valencetop(spin) = smallest_real
     condbottom(spin) = greatest_real
     do ikibz=1,ebands%nkpt
       nband_k = ebands%nband(ikibz + (spin-1)*ebands%nkpt)
       do band=1,nband_k
         if (ebands%occ(band,ikibz,spin) / maxocc > one-tol6 .and. valencetop(spin) < ebands%eig(band,ikibz,spin)) then
           valencetop(spin) = ebands%eig(band,ikibz,spin)
         end if
         if (ebands%occ(band,ikibz,spin) / maxocc < tol6 .and. condbottom(spin) > ebands%eig(band,ikibz,spin)) then
           condbottom(spin) = ebands%eig(band,ikibz,spin)
         end if
       end do
     end do
   end do

   vtop = MAXVAL(valencetop)
   cbot = MINVAL(condbottom)

   write(msg,'(a,f8.4,3a,f8.4,a)') &
    ' Top of valence: ', vtop * Ha_eV," (eV)", ch10, &
    ' Bottom of conduction: ', cbot * Ha_eV, " (eV)"
   call wrtout(std_out, msg)
   if (ebands%nsppol == 2) then
     if (ABS(vtop - MINVAL(valencetop)) > tol6) then
       call wrtout(std_out, sjoin(' Top of valence is spin: ', itoa(imax_loc(valencetop))))
     end if
     if (ABS(cbot - MAXVAL(condbottom)) > tol6) then
       call wrtout(std_out, ' Bottom of conduction is spin: ', itoa(imin_loc(condbottom)))
     end if
   end if

   ! Save results. Here I dont know if it is better to be consistent with the abinit convention i.e fermi=vtop
   ebands%entropy = zero
   ebands%fermie = (vtop + cbot) / 2
   if (ABS(cbot - vtop) < tol4) ebands%fermie = vtop ! To avoid error on the last digit
 end if

 call wrtout(std_out, sjoin(' Fermi level: ', ftoa(ebands%fermie * Ha_eV, fmt="f8.4"), " (eV)", ch10))

 ! Compute number of electrons for each spin channel.
 nelect_spin(:)=zero
 do spin=1,ebands%nsppol
   do ikibz=1,ebands%nkpt
     nband_k = ebands%nband(ikibz+(spin-1)*ebands%nkpt)
     nelect_spin(spin)= nelect_spin(spin) + ebands%wtk(ikibz) * sum(ebands%occ(1:nband_k,ikibz,spin))
   end do
 end do

 ndiff = ebands%nelect - SUM(nelect_spin)
 if (my_prtvol > 0) then
   write(msg,'(2a,f6.2,2a,f7.4)')ch10,&
    ' Total number of electrons: ', sum(nelect_spin),ch10,&
    ' Input and calculated no. of electrons differ by ',ndiff
   call wrtout(std_out, msg)
 end if

 if (ABS(ndiff) > 5.d-2*ebands%nelect) then
   write(msg,'(2a,2(a,es12.4))')&
    'Too large difference in number of electrons:,',ch10,&
    'Expected = ',ebands%nelect,' Calculated = ',sum(nelect_spin)
   MSG_ERROR(msg)
 end if

end subroutine ebands_update_occ
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_set_scheme
!! NAME
!! ebands_set_scheme
!!
!! FUNCTION
!! Set the occupation scheme and re-calculate new occupation numbers,
!! the Fermi level and the Max occupied band index for each spin channel starting
!! from the the knowledge of eigenvalues. See ebands_update_occ for more info.
!!
!! INPUTS
!! occopt=Occupation options (see input variable)
!! tsmear=Temperature of smearing.
!! spinmagntarget=if differ from -99.99d0, fix the spin polarization (in Bohr magneton)
!! prtvol=Verbosity level (0 for lowest level)
!! [update_occ]=False to avoid recomputing occupation factors (mainly used when a call to set_scheme is followed
!!  by another call to set_nelect (update_occ is expensive for large k-meshes). Default: True.
!!
!! PARENTS
!!      eph,m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_set_scheme(ebands, occopt, tsmear, spinmagntarget, prtvol, update_occ)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(inout) :: ebands
 integer,intent(in) :: occopt
 integer,intent(in) :: prtvol
 real(dp),intent(in) :: tsmear, spinmagntarget
 logical,optional,intent(in) :: update_occ

!Local variables-------------------------------
!scalars
 real(dp),parameter :: stmbias0 = zero
 logical :: my_update_occ

! *************************************************************************

 my_update_occ = .True.; if (present(update_occ)) my_update_occ = update_occ

 if (prtvol > 10) then
   call wrtout(std_out, " Changing occupation scheme in electron bands")
   call wrtout(std_out, sjoin(" occopt:", itoa(ebands%occopt), " ==> ", itoa(occopt)))
   call wrtout(std_out, sjoin(" tsmear:", ftoa(ebands%tsmear), " ==> ", ftoa(tsmear)))
 end if

 ebands%occopt = occopt; ebands%tsmear = tsmear

 if (my_update_occ) then
   call ebands_update_occ(ebands, spinmagntarget, stmbias0, prtvol=prtvol)
   if (prtvol > 10) call wrtout(std_out, sjoin(' Fermi level is now:', ftoa(ebands%fermie)))
 end if

end subroutine ebands_set_scheme
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_set_fermie
!! NAME
!! ebands_set_fermie
!!
!! FUNCTION
!! Set the new Fermi level from eigenenergies eigen and change the number of electrons
!! Compute also new occupation numbers at each k point, from eigenenergies eigen, according to the
!! smearing scheme defined by occopt (and smearing width tsmear or tphysel) as well as
!! entropy and derivative of occupancies wrt the energy for each band and k point.
!!
!! INPUTS
!! fermie=New fermi level
!!
!! OUTPUT
!! msg=String describing the changes in fermie and nelect.
!!
!! NOTES
!! The routine assumes metallic occupation scheme and will abort it this condition is not satisfied.
!! Use ebands_set_scheme before calling this routine, if you have a semiconductor.
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_set_fermie(ebands, fermie, msg)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(inout) :: ebands
 real(dp),intent(in) :: fermie
 character(len=*),intent(out) :: msg

!Local variables-------------------------------
!scalars
 integer,parameter :: option1=1, unitdos0=0
 integer :: mband,nkpt,nsppol
 real(dp),parameter :: dosdeltae0=zero
 real(dp) :: prev_fermie,prev_nelect,maxocc
!arrays
 real(dp),allocatable :: doccde(:),occ(:),eigen(:)

! *************************************************************************

 if (.not. ebands_has_metal_scheme(ebands)) then
   MSG_ERROR("set_fermie assumes a metallic occupation scheme. Use ebands_set_scheme before calling ebands_set_fermie!")
 end if

 prev_fermie = ebands%fermie; prev_nelect = ebands%nelect

 ! newocc assumes eigenvalues and occupations packed in 1d-vector!!
 mband  = ebands%mband
 nkpt   = ebands%nkpt
 nsppol = ebands%nsppol
 maxocc = two/(nsppol*ebands%nspinor)

 ABI_MALLOC(eigen,(mband*nkpt*nsppol))
 call get_eneocc_vect(ebands,'eig',eigen)
 ABI_MALLOC(occ,(mband*nkpt*nsppol))
 ABI_MALLOC(doccde,(mband*nkpt*nsppol))

 ! Get the total number of electrons nelect, given the new fermi energy.
 call getnel(doccde,dosdeltae0,eigen,ebands%entropy,fermie,maxocc,mband,ebands%nband,&
   ebands%nelect,nkpt,nsppol,occ,ebands%occopt,option1,ebands%tphysel,ebands%tsmear,unitdos0,ebands%wtk)

 ! Save changes in ebands%.
 ebands%fermie = fermie
 call put_eneocc_vect(ebands,'occ'   ,occ)
 call put_eneocc_vect(ebands,'doccde',doccde)
 ABI_FREE(eigen)
 ABI_FREE(occ)
 ABI_FREE(doccde)

 write(msg,"(2(a,es16.6),a,2(a,es16.6))") &
   " Old fermi level: ",prev_fermie,", with nelect: ",prev_nelect,ch10, &
   " New fermi level: ",ebands%fermie,", with nelect: ",ebands%nelect

end subroutine ebands_set_fermie
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_set_nelect
!! NAME
!! ebands_set_nelect
!!
!! FUNCTION
!! Set the new number of electrons. recompute Fermi level from eigenenergies
!! and new occupation numbers according to the smearing scheme defined by occopt
!! (and smearing width tsmear or tphysel) as well as
!! entropy and derivative of occupancies wrt the energy for each band and k point.
!!
!! INPUTS
!!  nelect=New number of electrons
!!  spinmagntarget=if differ from -99.99d0, fix the spin polarization (in Bohr magneton)
!!  [prtvol]=Verbosity level
!!
!! OUTPUT
!! msg=String describing the changes in fermie and nelect.
!!
!! NOTES
!! The routine assumes metallic occupation scheme and will abort it this condition is not satisfied.
!! Use ebands_set_scheme before calling this routine, if you have a semiconductor.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_set_nelect(ebands, nelect, spinmagntarget, msg, prtvol)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(inout) :: ebands
 integer,optional,intent(in) :: prtvol
 real(dp),intent(in) :: nelect,spinmagntarget
 character(len=*),intent(out) :: msg

!Local variables-------------------------------
!scalars
 integer :: my_prtvol
 real(dp) :: prev_fermie,prev_nelect

! *************************************************************************

 my_prtvol = 0; if (present(prtvol)) my_prtvol = prtvol

 if (.not. ebands_has_metal_scheme(ebands)) then
   MSG_ERROR("set_nelect assumes a metallic occupation scheme. Use ebands_set_scheme!")
 end if

 prev_fermie = ebands%fermie; prev_nelect = ebands%nelect
 ebands%nelect = nelect
 call ebands_update_occ(ebands, spinmagntarget, prtvol=my_prtvol)

 write(msg,"(2(a,es16.6),a,2(a,es16.6))")&
   " Old fermi level: ",prev_fermie,", with nelect: ",prev_nelect,ch10,&
   " New fermi level: ",ebands%fermie,", with nelect: ",ebands%nelect
 call wrtout(std_out, msg)

end subroutine ebands_set_nelect
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_get_muT_with_fd
!! NAME
!! ebands_get_muT_with_fd
!!
!! FUNCTION
!!  Compute the Fermi level for different temperatures using Fermi-Dirac occupation function (physical T)
!!  Use ebands%nelect provided in input. Does not change input ebands.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_get_muT_with_fd(self, ntemp, kTmesh, spinmagntarget, prtvol, mu_e, comm)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(in) :: self
 integer,intent(in) :: ntemp, prtvol, comm
 real(dp),intent(in) :: spinmagntarget
 real(dp),intent(in) :: kTmesh(ntemp)
 real(dp),intent(out) :: mu_e(ntemp)

!Local variables-------------------------------
!scalars
 integer,parameter :: occopt3 = 3
 integer :: ierr, it, nprocs, my_rank
 real(dp) :: nelect
 real(dp) :: cpu, wall, gflops
 type(ebands_t) :: tmp_ebands
 character(len=500) :: msg

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 call cwtime(cpu, wall, gflops, "start")

 call ebands_copy(self, tmp_ebands)

 mu_e = zero

 do it=1,ntemp
   if (mod(it, nprocs) /= my_rank) cycle ! MPI parallelism inside comm.

   ! Use Fermi-Dirac occopt
   call ebands_set_scheme(tmp_ebands, occopt3, kTmesh(it), spinmagntarget, prtvol)
   mu_e(it) = tmp_ebands%fermie
   !
   ! Check that the total number of electrons is correct
   ! This is to trigger problems as the routines that calculate the occupations in ebands_set_nelect
   ! are different from the occ_fd that will be used in the rest of the code.
   nelect = ebands_calc_nelect(tmp_ebands, kTmesh(it), mu_e(it))

   if (abs(nelect - self%nelect) > tol6) then
     ! For T = 0 the number of occupied states goes in discrete steps (according to the k-point sampling)
     ! for finite doping its hard to find nelect that exactly matches self%nelect.
     ! in this case we print a warning
     write(msg,'(3(a,f10.6))') &
       'Calculated number of electrons nelect: ',nelect, &
       ' does not correspond with ebands%nelect: ',tmp_ebands%nelect,' for kT = ',kTmesh(it)
     if (kTmesh(it) == 0) then
       MSG_WARNING(msg)
     else
       MSG_ERROR(msg)
     end if
   end if
 end do ! it

 call ebands_free(tmp_ebands)
 call xmpi_sum(mu_e, comm, ierr)
 call cwtime_report(" ebands_get_muT_with_fd", cpu, wall, gflops, end_str=ch10)

end subroutine ebands_get_muT_with_fd
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_calc_nelect
!! NAME
!! ebands_calc_nelect
!!
!! FUNCTION
!!  Compute nelect from Fermi level and Temperature.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

real(dp) pure function ebands_calc_nelect(self, kt, fermie) result(nelect)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(in) :: self
 real(dp),intent(in) :: kt, fermie

!Local variables-------------------------------
!scalars
 integer :: spin, ik, ib
 real(dp) :: ofact

! *************************************************************************

 ofact = two / (self%nsppol * self%nspinor)
 nelect = zero
 do spin=1,self%nsppol
   do ik=1,self%nkpt
     do ib=1,self%nband(ik + (spin-1)*self%nkpt)
       nelect = nelect + self%wtk(ik) * occ_fd(self%eig(ib,ik,spin), kt, fermie)
     end do
   end do
 end do

 nelect = ofact * nelect

end function ebands_calc_nelect
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_report_gap
!! NAME
!! ebands_report_gap
!!
!! FUNCTION
!!  Print info on the fundamental and direct gap.
!!
!! INPUTS
!!  ebands<ebands_t>=Info on the band structure, the smearing technique and the physical temperature used.
!!  [header]=Optional title.
!!  [kmask]=Logical mask used to exclude k-points.
!!  [unit]=Optional unit for output (std_out if not specified)
!!  [mode_paral]=Either "COLL" or "PERS", former is default.
!!
!! OUTPUT
!!  writing.
!!  [gaps(3,nsppol)]=Fundamental and direct gaps. The third index corresponds to a "status":
!!      0.0dp if gaps were not computed (because there are only valence bands);
!!     -1.0dp if the system (or spin-channel) is metallic;
!!      1.0dp if the gap has been computed.
!!
!! PARENTS
!!      gstate,m_exc_diago,m_sigmaph,setup_bse,setup_bse_interp,setup_sigma
!!      sigma
!!
!! TODO
!!  Can be replaced by ebands_print_gaps
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_report_gap(ebands, header, kmask, unit, mode_paral, gaps)

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: unit
 character(len=4),intent(in),optional :: mode_paral
 character(len=*),intent(in),optional :: header
 class(ebands_t),intent(in)  :: ebands
!arrays
 real(dp),optional,intent(out) :: gaps(3,ebands%nsppol)
 logical,optional,intent(in) ::  kmask(ebands%nkpt)

!Local variables-------------------------------
!scalars
 integer :: ikibz,nband_k,spin,nsppol,ikopt,ivk,ick,ivb,icb,my_unt,first
 real(dp),parameter :: tol_fermi=tol6
 real(dp) :: fun_gap,opt_gap
 logical :: ismetal
 character(len=4) :: my_mode
 character(len=500) :: msg
!arrays
 integer :: val_idx(ebands%nkpt,ebands%nsppol)
 real(dp) :: top_valence(ebands%nkpt),bot_conduct(ebands%nkpt)
 logical :: my_kmask(ebands%nkpt)

! *********************************************************************

 nsppol = ebands%nsppol

 my_unt =std_out; if (PRESENT(unit      )) my_unt =unit
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral
 my_kmask=.TRUE.; if (PRESENT(kmask     )) my_kmask=kmask

 if (PRESENT(gaps)) gaps=zero

 val_idx(:,:) = ebands_get_valence_idx(ebands, tol_fermi)
 first=0

 ! Initialize the return status for the gaps
 if (PRESENT(gaps)) gaps(1:3,1:nsppol)=zero

 do spin=1,nsppol

   ! No output if system i metallic
   ismetal=ANY(val_idx(:,spin)/=val_idx(1,spin))
   if (ismetal) then
     if (PRESENT(gaps)) gaps(3,nsppol)=-one
     CYCLE
   endif

   first=first+1
   if (first==1) then
     msg=ch10
     if (PRESENT(header)) msg=ch10//' === '//TRIM(ADJUSTL(header))//' === '
     call wrtout(my_unt,msg,my_mode)
   end if

   ivb=val_idx(1,spin)
   icb=ivb+1

   do ikibz=1,ebands%nkpt
     if (.not.my_kmask(ikibz)) CYCLE
     nband_k = ebands%nband(ikibz+(spin-1)*ebands%nkpt)
     top_valence(ikibz) = ebands%eig(ivb,ikibz,spin)
     if (icb>nband_k) then
       GOTO 10 ! Only occupied states are present, no output!
     end if
     bot_conduct(ikibz) = ebands%eig(icb,ikibz,spin)
   end do

   ! Get minimum of the direct Gap
   ikopt= imin_loc(bot_conduct-top_valence,MASK=my_kmask)
   opt_gap=bot_conduct(ikopt)-top_valence(ikopt)

   ! Get fundamental Gap
   ick = imin_loc(bot_conduct,MASK=my_kmask)
   ivk = imax_loc(top_valence,MASK=my_kmask)
   fun_gap = ebands%eig(icb,ick,spin)-ebands%eig(ivb,ivk,spin)

   write(msg,'(a,i2,a,2(a,f8.4,a,3f8.4,a),33x,a,3f8.4)')&
    '  >>>> For spin ',spin,ch10,&
    '   Minimum direct gap = ',opt_gap*Ha_eV,' [eV], located at k-point      : ',ebands%kptns(:,ikopt),ch10,&
    '   Fundamental gap    = ',fun_gap*Ha_eV,' [eV], Top of valence bands at : ',ebands%kptns(:,ivk),ch10,  &
                                              '      Bottom of conduction at : ',ebands%kptns(:,ick)
   call wrtout(my_unt,msg,my_mode)

   if (present(gaps)) gaps(:,spin) = [fun_gap, opt_gap, one]
 end do !spin

 return

 10 continue
 call wrtout(std_out, "Not enough states to calculate the band gap.", "COLL")

end subroutine ebands_report_gap
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_ncwrite
!! NAME
!! ebands_ncwrite
!!
!! FUNCTION
!!  Writes the content of an ebands_t object to a NETCDF file
!!  according to the ETSF-IO specifications. Return nf90_noerr if success.
!!
!! INPUTS
!!  ncid =NC file handle
!!
!! PARENTS
!!      dfpt_looppert,eig2tot,ioarr,m_ebands,m_iowf,m_shirley,pawmkaewf,sigma
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

integer function ebands_ncwrite(ebands, ncid) result(ncerr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 class(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
!scalars
#ifdef HAVE_NETCDF
 integer :: ii,nelect_int
 logical :: write_ngkpt
 character(len=etsfio_charlen) :: smearing,k_dependent
!arrays
 integer :: ngkpt(3)

! *************************************************************************

 smearing = nctk_string_from_occopt(ebands%occopt)

 ! ==============================================
 ! === Write the dimensions specified by ETSF ===
 ! ==============================================
 ncerr = nctk_def_dims(ncid, [ &
   nctkdim_t("max_number_of_states", ebands%mband), &
   nctkdim_t("number_of_spinor_components", ebands%nspinor), &
   nctkdim_t("number_of_spins", ebands%nsppol), &
   nctkdim_t("number_of_kpoints", ebands%nkpt), &
   nctkdim_t("nshiftk_orig", ebands%nshiftk_orig), &
   nctkdim_t("nshiftk", ebands%nshiftk)], &
   defmode=.True.)
 NCF_CHECK(ncerr)

 ! FIXME
 ! Unofficial variables. Notes:
 ! 1) ETSF-IO does not support nshifts > 1
 ! 2) shiftk_orig, nshiftk_orig refers to the values specified in the input (most useful ones).
 ! 3) shiftk, kptrlatt refers to the values computed in inkpts.
 ! 4) Should define a protocol so that abipy understands if we have a path or a mesh.
 !write_kptrlatt = (SUM(ABS(ebands%kptrlatt))/=0)
 !write_kptrlatt = (ebands%kptopt /= 0)

 ngkpt = 0; write_ngkpt = .False.
 if (isdiagmat(ebands%kptrlatt) .and. ebands%nshiftk == 1) then
    write_ngkpt = .True.
    do ii=1,3
      ngkpt(ii) = ebands%kptrlatt(ii, ii)
    end do
    ncerr = nctk_def_dims(ncid, nctkdim_t('ngkpt_nshiftk', ebands%nshiftk_orig))
    NCF_CHECK(ncerr)
 end if

 ! Define k-points
 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t("reduced_coordinates_of_kpoints", "dp", "number_of_reduced_dimensions, number_of_kpoints"), &
   nctkarr_t("kpoint_weights", "dp", "number_of_kpoints"), &
   nctkarr_t("monkhorst_pack_folding", "int", "number_of_vectors") &
 ])
 NCF_CHECK(ncerr)

 ! Define states section.
 ncerr = nctk_def_arrays(ncid, [&
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

 ! Some variables require the specifications of units.
 NCF_CHECK(nctk_set_atomic_units(ncid, "eigenvalues"))
 NCF_CHECK(nctk_set_atomic_units(ncid, "fermi_energy"))

 k_dependent = "no"; if (any(ebands%nband(1) /= ebands%nband)) k_dependent = "yes"
 NCF_CHECK(nf90_put_att(ncid, vid("number_of_states"), "k_dependent", k_dependent))

 ! Write data.
 ! 1) Electrons.
 ! NB: In etsf_io the number of electrons is declared as integer.
 ! We use abinit nelect to store the value as real(dp).
 nelect_int = nint(ebands%nelect)

 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid("fermi_energy"), ebands%fermie))
 NCF_CHECK(nf90_put_var(ncid, vid("number_of_electrons"), nelect_int))
 NCF_CHECK(nf90_put_var(ncid, vid("smearing_width"), ebands%tsmear))
 NCF_CHECK(nf90_put_var(ncid, vid("number_of_states"), ebands%nband, count=[ebands%nkpt, ebands%nsppol]))
 NCF_CHECK(nf90_put_var(ncid, vid("eigenvalues"), ebands%eig))
 NCF_CHECK(nf90_put_var(ncid, vid("occupations"), ebands%occ))
 NCF_CHECK(nf90_put_var(ncid, vid("smearing_scheme"), smearing))

 ! K-points
 NCF_CHECK(nf90_put_var(ncid, vid("reduced_coordinates_of_kpoints"), ebands%kptns))
 NCF_CHECK(nf90_put_var(ncid, vid("kpoint_weights"), ebands%wtk))

 if (write_ngkpt) then
   NCF_CHECK(nf90_put_var(ncid, vid("monkhorst_pack_folding"), ngkpt))
 end if

 ! ===========================================================
 ! === Write abinit-related stuff (not covered by ETSF-IO) ===
 ! ===========================================================
 ! Define variables.
 NCF_CHECK(nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "occopt", "kptopt"], defmode=.True.))
 NCF_CHECK(nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "tphysel", "charge", "nelect"]))

 ncerr = nctk_def_arrays(ncid, nctkarr_t('istwfk', "i", 'number_of_kpoints'))
 NCF_CHECK(ncerr)

 ! Abinit variables defining the K-point sampling.
 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t('kptrlatt_orig', "i", 'number_of_reduced_dimensions, number_of_reduced_dimensions'), &
   nctkarr_t('shiftk_orig',  "dp", 'number_of_reduced_dimensions, nshiftk_orig'), &
   nctkarr_t('kptrlatt', "i", 'number_of_reduced_dimensions, number_of_reduced_dimensions'), &
   nctkarr_t('shiftk',  "dp", 'number_of_reduced_dimensions, nshiftk') &
 ])
 NCF_CHECK(ncerr)

 if (write_ngkpt) then
   ncerr = nctk_def_arrays(ncid, nctkarr_t('ngkpt_shiftk', "dp", "number_of_reduced_dimensions, ngkpt_nshiftk"))
   NCF_CHECK(ncerr)
 end if

 ! Write Abinit variables
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid("tphysel"), ebands%tphysel))
 NCF_CHECK(nf90_put_var(ncid, vid("occopt"), ebands%occopt))
 NCF_CHECK(nf90_put_var(ncid, vid("istwfk"), ebands%istwfk))
 NCF_CHECK(nf90_put_var(ncid, vid("kptopt"), ebands%kptopt))
 NCF_CHECK(nf90_put_var(ncid, vid("charge"), ebands%charge))
 NCF_CHECK(nf90_put_var(ncid, vid("nelect"), ebands%nelect))
 NCF_CHECK(nf90_put_var(ncid, vid('kptrlatt_orig'), ebands%kptrlatt_orig))
 NCF_CHECK(nf90_put_var(ncid, vid('shiftk_orig'), ebands%shiftk_orig))
 NCF_CHECK(nf90_put_var(ncid, vid('kptrlatt'),ebands%kptrlatt))
 NCF_CHECK(nf90_put_var(ncid, vid('shiftk'), ebands%shiftk))

 if (write_ngkpt) then
   NCF_CHECK(nf90_put_var(ncid, vid('ngkpt_shiftk'), ebands%shiftk_orig))
 end if

#else
 MSG_ERROR("netcdf support is not activated. ")
#endif

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end function ebands_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_ncwrite_path
!! NAME
!! ebands_ncwrite_path
!!
!! FUNCTION
!!  Writes the content of an ebands_t object to a NETCDF file
!!
!! INPUTS
!!  cryst=Crystal structure
!!  path=File name
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function ebands_ncwrite_path(ebands, cryst, path) result(ncerr)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: cryst
 character(len=*),intent(in) :: path

!Local variables-------------------------------
!scalars
 integer :: ncid

! *************************************************************************

 ncerr = -1
#ifdef HAVE_NETCDF
 ncerr = nf90_noerr
 if (file_exists(path)) then
    NCF_CHECK(nctk_open_modify(ncid, path, xmpi_comm_self))
 else
   ncerr = nctk_open_create(ncid, path, xmpi_comm_self)
   NCF_CHECK_MSG(ncerr, sjoin("Creating", path))
 end if

 NCF_CHECK(cryst%ncwrite(ncid))
 NCF_CHECK(ebands_ncwrite(ebands, ncid))
 NCF_CHECK(nf90_close(ncid))
#endif

end function ebands_ncwrite_path
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_get_edos
!! NAME
!!  ebands_get_edos
!!
!! FUNCTION
!!  Calculate the electronic density of states from ebands_t
!!
!! INPUTS
!!  ebands<ebands_t>=Band structure object.
!!  cryst<cryst_t>=Info on the crystalline structure.
!!  intmeth= 1 for Gaussian, 2 or -2 for tetrahedra (-2 if Blochl corrections must be included).
!!    If nkpt == 1 (Gamma only), the routine fallbacks to gaussian method.
!!  step=Step on the linear mesh in Ha. If <0, the routine will use the mean of the energy level spacing
!!  broad=Gaussian broadening, If <0, the routine will use a default
!!    value for the broadening computed from the mean of the energy level spacing.
!!    No meaning for tetrahedra
!!  comm=MPI communicator
!!
!! OUTPUT
!!  edos<edos_t>=Electronic DOS and IDOS.
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!
!! SOURCE

type(edos_t) function ebands_get_edos(ebands, cryst, intmeth, step, broad, comm) result(edos)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: intmeth,comm
 real(dp),intent(in) :: step,broad
 class(ebands_t),target,intent(in)  :: ebands
 type(crystal_t),intent(in) :: cryst

!Local variables-------------------------------
!scalars
 integer :: nw,spin,band,ikpt,ief,nproc,my_rank,mpierr,cnt,ierr,bcorr
 real(dp) :: max_ene,min_ene,wtk,max_occ
 character(len=500) :: msg
 type(htetra_t) :: tetra
!arrays
 real(dp) :: eminmax_spin(2,ebands%nsppol)
 real(dp),allocatable :: wme0(:),wdt(:,:),tmp_eigen(:)

! *********************************************************************

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 ierr = 0

 edos%nkibz = ebands%nkpt; edos%nsppol = ebands%nsppol
 edos%intmeth = intmeth

 if (ebands%nkpt == 1) then
   MSG_COMMENT("Cannot use tetrahedra for e-DOS when nkpt == 1. Switching to gaussian method")
   edos%intmeth = 1
 end if

 edos%broad = broad; edos%step = step

 ! Compute the linear mesh so that it encloses all bands.
 eminmax_spin = ebands_get_minmax(ebands, "eig")
 min_ene = minval(eminmax_spin(1,:)); min_ene = min_ene - 0.1_dp * abs(min_ene)
 max_ene = maxval(eminmax_spin(2,:)); max_ene = max_ene + 0.1_dp * abs(max_ene)

 nw = nint((max_ene - min_ene)/edos%step) + 1; edos%nw = nw

 ABI_MALLOC(edos%mesh, (nw))
 edos%mesh = arth(min_ene, edos%step, nw)

 ABI_CALLOC(edos%gef, (0:edos%nsppol))
 ABI_CALLOC(edos%dos,  (nw, 0:edos%nsppol))
 ABI_CALLOC(edos%idos, (nw, 0:edos%nsppol))

 select case (edos%intmeth)
 case (1)
   !call wrtout(std_out, " Computing electron-DOS with Gaussian method")
   !call wrtout(std_out, sjoin(" broadening: ", ftoa(edos%broad * Ha_eV), " (eV), step: ", ftoa(edos%step * Ha_eV), "(eV), npts: ",
   !itoa(nw)))
   ! Gaussian
   ABI_MALLOC(wme0, (nw))
   cnt = 0
   do spin=1,edos%nsppol
     do ikpt=1,ebands%nkpt
       cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle  ! MPI parallelism
       wtk = ebands%wtk(ikpt)
       do band=1,ebands%nband(ikpt+(spin-1)*ebands%nkpt)
          wme0 = edos%mesh - ebands%eig(band, ikpt, spin)
          edos%dos(:, spin) = edos%dos(:, spin) + wtk * gaussian(wme0, edos%broad)
       end do
     end do
   end do
   ABI_FREE(wme0)
   call xmpi_sum(edos%dos, comm, mpierr)

 case (2, -2)
   !call wrtout(std_out, " Computing electron-DOS with tetrahedron method")
   ! Consistency test
   if (any(ebands%nband /= ebands%nband(1)) ) MSG_ERROR('for tetrahedra, nband(:) must be constant')

   ! Build tetra object.
   tetra = tetra_from_kptrlatt(cryst, ebands%kptopt, ebands%kptrlatt, &
     ebands%nshiftk, ebands%shiftk, ebands%nkpt, ebands%kptns, comm, msg, ierr)
   if (ierr /= 0) MSG_ERROR(msg)

   ! For each spin and band, interpolate over kpoints,
   ! calculate integration weights and DOS contribution.
   ABI_MALLOC(tmp_eigen, (ebands%nkpt))
   ABI_MALLOC(wdt, (nw, 2))

   bcorr = 0; if (intmeth == -2) bcorr = 1
   cnt = 0
   do spin=1,ebands%nsppol
     do band=1,ebands%nband(1)
       ! For each band get its contribution
       tmp_eigen = ebands%eig(band,:,spin)
       do ikpt=1,ebands%nkpt
         cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism.

         ! Calculate integration weights at each irred k-point (Blochl et al PRB 49 16223 [[cite:Bloechl1994a]])
         call tetra%get_onewk(ikpt, bcorr, nw, ebands%nkpt, tmp_eigen, min_ene, max_ene, one, wdt)

         edos%dos(:,spin) = edos%dos(:,spin) + wdt(:, 1)*ebands%wtk(ikpt)
         ! IDOS is computed afterwards with simpson
         !edos%idos(:,spin) = edos%idos(:,spin) + wdt(:, 2)
       end do ! ikpt
     end do ! band
   end do ! spin

   call xmpi_sum(edos%dos, comm, mpierr)

   ! Free memory
   ABI_FREE(tmp_eigen)
   ABI_FREE(wdt)
   call tetra%free()

   ! Filter so that dos[i] is always >= 0 and idos is monotonic
   ! IDOS is computed afterwards with simpson
   where (edos%dos(:,1:) <= zero)
     edos%dos(:,1:) = zero
   end where

 case default
   MSG_ERROR(sjoin("Wrong integration method:", itoa(intmeth)))
 end select

 ! Compute total DOS and IDOS
 max_occ = two/(ebands%nspinor*ebands%nsppol)
 edos%dos(:, 0) = max_occ * sum(edos%dos(:,1:), dim=2)

 do spin=1,edos%nsppol
   call simpson_int(nw,edos%step,edos%dos(:,spin),edos%idos(:,spin))
 end do
 edos%idos(:, 0) = max_occ * sum(edos%idos(:,1:), dim=2)

 ! Use bisection to find fermi level.
 ! Warning: this code assumes idos[i+1] >= idos[i]. This condition may not be
 ! fullfilled if we use tetra and this is the reason why we have filtered the DOS.
 ief = bisect(edos%idos(:,0), ebands%nelect)

 ! Handle out of range condition.
 if (ief == 0 .or. ief == nw) then
   write(msg,"(3a)")&
    "Bisection could not find an initial guess for the Fermi level!",ch10,&
    "Possible reasons: not enough bands or wrong number of electrons"
   MSG_WARNING(msg)
   return
 end if

 ! TODO: Use linear interpolation to find an improved estimate of the Fermi level?
 edos%ief = ief
 do spin=0,edos%nsppol
   edos%gef(spin) = edos%dos(ief,spin)
 end do

 if (.False.) then
   write(std_out,*)"fermie from ebands: ",ebands%fermie
   write(std_out,*)"fermie from IDOS: ",edos%mesh(ief)
   write(std_out,*)"gef:from ebands%fermie: " ,edos%dos(bisect(edos%mesh, ebands%fermie), 0)
   write(std_out,*)"gef:from edos: " ,edos%gef(0)
 end if

end function ebands_get_edos
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/edos_free
!! NAME
!!  edos_free
!!
!! FUNCTION
!!  Free the memory allocated in edos_t
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!
!! SOURCE

subroutine edos_free(edos)

!Arguments ------------------------------------
 class(edos_t),intent(inout) :: edos

! *********************************************************************

!real
 ABI_SFREE(edos%mesh)
 ABI_SFREE(edos%dos)
 ABI_SFREE(edos%idos)
 ABI_SFREE(edos%gef)

end subroutine edos_free
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/edos_write
!! NAME
!! edos_write
!!
!! FUNCTION
!! Write results to file.
!!
!! INPUTS
!!  edos<edos_t>=DOS container
!!  path=File name.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!
!! SOURCE

subroutine edos_write(edos, path)

!Arguments ------------------------------------
 character(len=*),intent(in) :: path
 class(edos_t),intent(in) :: edos

!Local variables-------------------------------
 integer :: iw,spin,unt
 real(dp) :: cfact,efermi
 character(len=500) :: msg

! *************************************************************************

 ! Convert everything into eV
 ! I know that Abinit should use Ha but Hartrees are not readable.
 ! Please don't change this code, in case add an optional argument to specify different units.
 cfact = Ha_eV

 if (open_file(path, msg, newunit=unt, form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 ! Write header.
 write(unt,'(a)')'# Electron density of states: Energy in eV, DOS in states/eV per unit cell.'
 write(unt,"(a)")"# The zero of energies corresponds to the Fermi level."

 select case (edos%intmeth)
 case (1)
   write(unt,'(a,es16.8,a,i0)')&
     '# Gaussian method with smearing= ',edos%broad*cfact,' [eV], nkibz= ',edos%nkibz
 case (2)
   write(unt,'(a,i0)')'# Tetrahedron method, nkibz= ',edos%nkibz
 case default
   MSG_ERROR(sjoin("Wrong method:", itoa(edos%intmeth)))
 end select

 if (edos%ief == 0) then
   write(unt,'(a)')'# Fermi level: None'
   efermi = zero
 else
   write(unt,'(a,es16.8,a)')'# Fermi level: ',edos%mesh(edos%ief)*cfact," [eV]"
   efermi = edos%mesh(edos%ief)
 end if

 ! Write data.
 write(unt,"(a)")"# Energy           DOS_TOT          IDOS_TOT         DOS[spin=UP]     IDOS[spin=UP] ..."
 do iw=1,edos%nw
   write(unt,'(es17.8)',advance='no')(edos%mesh(iw) - efermi)*cfact
   do spin=0,edos%nsppol
     write(unt,'(2es17.8)',advance='no')edos%dos(iw,spin)/cfact,edos%idos(iw,spin)
   end do
   write(unt,*)
 end do

 close(unt)

end subroutine edos_write
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/edos_ncwrite
!! NAME
!! edos_ncwrite
!!
!! FUNCTION
!!  Write results to netcdf file.
!!
!! INPUTS
!!  edos<edos_t>=DOS container
!!  ncid=NC file handle.
!!  [prefix]=String prepended to netcdf dimensions/variables (HDF5 poor-man groups)
!!   Empty string if not specified.
!!
!! OUTPUT
!!  ncerr= netcdf exit status.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function edos_ncwrite(edos, ncid, prefix) result(ncerr)

!Arguments ------------------------------------
 integer,intent(in) :: ncid
 class(edos_t),intent(in) :: edos
 character(len=*),optional,intent(in) :: prefix

!Local variables-------------------------------
 character(len=500) :: prefix_

! *************************************************************************

 prefix_ = ""; if (present(prefix)) prefix_ = trim(prefix)

#ifdef HAVE_NETCDF
 ! Define dimensions.
 ncerr = nctk_def_dims(ncid, [ &
   nctkdim_t("nsppol_plus1", edos%nsppol + 1), nctkdim_t("edos_nw", edos%nw)], defmode=.True., prefix=prefix_)
 NCF_CHECK(ncerr)

 ! Define variables
 NCF_CHECK(nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "edos_intmeth", "edos_nkibz", "edos_ief"], prefix=prefix_))
 NCF_CHECK(nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "edos_broad"], prefix=prefix_))

 ncerr = nctk_def_arrays(ncid, [ &
   nctkarr_t("edos_mesh", "dp", "edos_nw"), &
   nctkarr_t("edos_dos", "dp", "edos_nw, nsppol_plus1"), &
   nctkarr_t("edos_idos", "dp", "edos_nw, nsppol_plus1"), &
   nctkarr_t("edos_gef", "dp", "nsppol_plus1") &
 ],  prefix=prefix_)
 NCF_CHECK(ncerr)

 ! Write data.
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("edos_intmeth")), edos%intmeth))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("edos_nkibz")), edos%nkibz))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("edos_ief")), edos%ief))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("edos_broad")), edos%broad))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("edos_mesh")), edos%mesh))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("edos_dos")), edos%dos))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("edos_idos")), edos%idos))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("edos_gef")), edos%gef))

#else
 MSG_ERROR("netcdf library not available")
#endif

contains
  pure function pre(istr) result(ostr)
    character(len=*),intent(in) :: istr
    character(len=len_trim(prefix_) + len_trim(istr)+1) :: ostr
    ostr = trim(prefix_) // trim(istr)
  end function pre

end function edos_ncwrite
!!***

!!****f* m_ebands/edos_print
!! NAME
!! edos_print
!!
!! FUNCTION
!! Print DOS info to Fortran unit.
!!
!! INPUTS
!!  edos<edos_t>=DOS container
!!  [unit]=Unit number for output. Defaults to std_out
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!
!! SOURCE

subroutine edos_print(edos, unit, header)

!Arguments ------------------------------------
 class(edos_t),intent(in) :: edos
 integer,optional,intent(in) :: unit
 character(len=*),optional,intent(in) :: header

!Local variables-------------------------------
 integer :: unt

! *************************************************************************

 unt = std_out; if (present(unit)) unt = unit
 if (unt == dev_null) return

 if (present(header)) then
   write(unt, "(a)") ch10//' === '//trim(adjustl(header))//' === '
 else
   write(unt, "(a)") ch10
 end if

 select case (edos%intmeth)
 case (1)
   write(unt, "(a,f5.1,a)") " Gaussian method with broadening: ", edos%broad * Ha_meV, " (meV)"
 case (2)
   write(unt, "(a)")" Linear tetrahedron method."
 case (-2)
   write(unt, "(a)")" Linear tetrahedron method with Blochl corrections."
 case default
   MSG_ERROR(sjoin("Wrong intmeth:", itoa(edos%intmeth)))
 end select

 write(unt, "(a,f5.1,a, i0)")" Mesh step: ", edos%step * Ha_meV, " (meV) with npts: ", edos%nw
 write(unt, "(2(a,f5.1),a)")" From emin: ", edos%mesh(1) * Ha_eV, " to emax: ", edos%mesh(edos%nw) * Ha_eV, " (eV)"
 write(unt, "(a, i0)")" Number of k-points in the IBZ: ", edos%nkibz

 if (edos%ief == 0) then
   write(unt, "(a, /)")" edos%ief == 0 --> Cannot print quantities at the Fermi level."
   return
 end if

 write(unt,'(a,es16.8,a)')' Fermi level: ',edos%mesh(edos%ief) * Ha_eV, " (eV)"
 write(unt,"(a,es16.8)")" Total electron DOS at Fermi level in states/eV: ", edos%gef(0) / Ha_eV
 if (edos%nsppol == 2) then
   write(unt,"(a,es16.8)")"   g(eF) for spin up:  ", edos%gef(1) / Ha_eV
   write(unt,"(a,es16.8)")"   g(eF) for spin down:", edos%gef(2) / Ha_eV
 end if
 write(unt,"(a,f6.1)")" Total number of electrons at eF: ", edos%idos(edos%ief, 0)
 if (edos%nsppol == 2) then
   write(unt,"(a,es16.8)")"   N(eF) for spin up:  ", edos%idos(edos%ief, 1)
   write(unt,"(a,es16.8)")"   N(eF) for spin down:", edos%idos(edos%ief, 2)
 end if

 write(unt, "(a)")""

end subroutine edos_print
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_write_nesting
!! NAME
!! ebands_write_nesting
!!
!! FUNCTION
!! Calculate the nesting function and output data to file.
!!
!! INPUTS
!!  ebands<ebands_t>=the ebands_t datatype
!!  cryst<crystal_t>=Info on unit cell and symmetries.
!!  filepath=File name for output data.
!!  prtnest = flags governing the format of the output file. see mknesting.
!!  tsmear=Broadening used to approximation the delta function.
!!  fermie_nest
!!  qpath_vertices = vertices of the reciprocal space trajectory
!!
!! OUTPUT
!!  Return non-zero exist status if netsting factor cannot be produced.
!!  The errmsg string gives information on the error.
!!
!! SIDE EFFECTS
!!   Write data to file.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function ebands_write_nesting(ebands,cryst,filepath,prtnest,tsmear,fermie_nest,qpath_vertices,errmsg) result(skip)

!Arguments ------------------------------------
 class(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: prtnest
 real(dp),intent(in) :: tsmear,fermie_nest
 character(len=*),intent(in) :: filepath
 character(len=*),intent(out) :: errmsg
!arrays
 real(dp),intent(in) :: qpath_vertices(:,:)

!Local variables-------------------------------
!scalaras
 integer :: ikpt,spin,iband,nqpath
 real(dp) :: invgauwidth,prefact,fermie
!arrays
 real(dp), allocatable :: fs_weights(:,:,:)

! *********************************************************************

 skip = 0; errmsg = ""
 if (any(ebands%nband /= ebands%nband(1))) then
   errmsg = 'mknesting can not handle variable nband(1:nkpt). Skipped.'//&
     ch10//' Correct input file to get nesting output'
   skip = 1; return
 end if

 if (ebands%nshiftk /= 1) then
   errmsg = 'mknesting does not support nshiftk > 1. Change ngkpt and shiftk to have only one shift after inkpts'
   skip = 1; return
 end if

 ! FIXME: needs to be generalized to complete the k grid for one of the arguments to mknesting
 fermie = ebands%fermie
 nqpath = size(qpath_vertices, dim=2)

 ! Compute weights. Set sigma to 0.1 eV is tsmear is zero
 invgauwidth = one / (0.1_dp * eV_Ha); if (tsmear > tol10) invgauwidth = one / tsmear
 prefact = one / sqrt(pi) * invgauwidth

 ABI_MALLOC(fs_weights,(ebands%nband(1),ebands%nkpt,ebands%nsppol))

 do spin=1,ebands%nsppol
   do ikpt=1,ebands%nkpt
     do iband=1,ebands%nband(1)
       fs_weights(iband, ikpt, spin) = prefact * &
         exp(-(invgauwidth*(ebands%eig(iband,ikpt,spin)-(fermie + fermie_nest)))**2)
     end do
   end do
 end do

 if (any(ebands%kptopt == [3, 4])) then ! no symmetry
   call mknesting(ebands%nkpt,ebands%kptns,ebands%kptrlatt,ebands%nband(1),fs_weights,nqpath,&
     qpath_vertices,1,[zero, zero, zero],filepath,cryst%gprimd,cryst%gmet,prtnest,identity_3d)
 else
   call mknesting(ebands%nkpt,ebands%kptns,ebands%kptrlatt,ebands%nband(1),fs_weights,nqpath,&
     qpath_vertices,1, [zero, zero, zero], filepath,cryst%gprimd,cryst%gmet,prtnest,identity_3d,&
     nsym=cryst%nsym, symrec=cryst%symrec)
 end if

 ABI_FREE(fs_weights)

end function ebands_write_nesting
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_expandk
!! NAME
!! ebands_expandk
!!
!! FUNCTION
!!  Return a new object of type ebands_t corresponding to a list of k-points
!!  specified in input. Symmetry properties of the eigenvectors are used to
!!  symmetrize energies and occupation numbers.
!!
!! INPUTS
!!  inb<ebands_t>=Initial band structure with energies in the IBZ.
!!  ecut_eff=Effective cutoff energy i.e. ecut * dilatmx**2
!!  force_istwfk1=If True, istwfk if forced to 1 for all the k-points in the BZ.
!!
!! OUTPUT
!!  dksqmax=maximal value of the norm**2 of the difference between
!!    a kpt in the BZ and the closest k-point found in the inb%kpts set, using symmetries.
!!  bz2ibz(nkpt2,6)=describe k point number of kpt1 that allows to
!!    generate wavefunctions closest to given kpt2
!!      bz2ibz(:,1)=k point number of kptns1
!!      bz2ibz(:,2)=symmetry operation to be applied to kpt1, to give kpt1a
!!        (if 0, means no symmetry operation, equivalent to identity )
!!      bz2ibz(:,3:5)=shift in reciprocal space to be given to kpt1a,
!!        to give kpt1b, that is the closest to kpt2.
!!      bz2ibz(:,6)=1 if time-reversal was used to generate kpt1a from kpt1, 0 otherwise
!!  outb<ebands_t>=band structure with energies in the BZ.
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_expandk(inb, cryst, ecut_eff, force_istwfk1, dksqmax, bz2ibz, outb)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: ecut_eff
 real(dp),intent(out) :: dksqmax
 logical,intent(in) :: force_istwfk1
 class(ebands_t),intent(in) :: inb
 class(ebands_t),intent(out) :: outb
 class(crystal_t),intent(in) :: cryst
!arrays
 integer,allocatable,intent(out) :: bz2ibz(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: istwfk_1=1,kptopt3=3
 integer :: nkfull,timrev,bantot,sppoldbl,npw_k,nsppol,istw
 integer :: ik_ibz,ikf,isym,itimrev,spin,mband,my_nkibz,comm
 logical :: isirred_k
 !character(len=500) :: msg
!arrays
 integer :: g0(3)
 integer,allocatable :: istwfk(:),nband(:,:),npwarr(:),kg_k(:,:)
 real(dp),allocatable :: kfull(:,:),doccde(:),eig(:),occ(:),wtk(:),my_kibz(:,:)
 real(dp),allocatable :: doccde_3d(:,:,:),eig_3d(:,:,:),occ_3d(:,:,:)

! *********************************************************************

 ABI_CHECK(inb%kptopt /= 0, "ebands_expandk does not support kptopt == 0")

 comm = xmpi_comm_self
 nsppol = inb%nsppol

 ! Note kptopt=3
 call kpts_ibz_from_kptrlatt(cryst, inb%kptrlatt, kptopt3, inb%nshiftk, inb%shiftk, &
   my_nkibz, my_kibz, wtk, nkfull, kfull) ! new_kptrlatt, new_shiftk)

 ABI_FREE(my_kibz)
 ABI_FREE(wtk)

 ! Costruct full BZ and create mapping BZ --> IBZ
 ! Note:
 !   - we don't change the value of nsppol hence sppoldbl is set to 1
 !   - we use symrel so that bz2ibz can be used to reconstruct the wavefunctions.
 !
 sppoldbl = 1 !; if (any(cryst%symafm == -1) .and. inb%nsppol == 1) sppoldbl=2
 ABI_MALLOC(bz2ibz, (nkfull*sppoldbl,6))

 timrev = kpts_timrev_from_kptopt(inb%kptopt)
 call listkk(dksqmax, cryst%gmet, bz2ibz, inb%kptns, kfull, inb%nkpt, nkfull, cryst%nsym, &
   sppoldbl, cryst%symafm, cryst%symrel, timrev, comm, use_symrec=.False.)

 ABI_MALLOC(wtk, (nkfull))
 wtk = one / nkfull ! weights normalized to one

 ABI_MALLOC(istwfk, (nkfull))
 ABI_MALLOC(nband, (nkfull, nsppol))
 ABI_MALLOC(npwarr, (nkfull))

 if (any(cryst%symrel(:,:,1) /= identity_3d) .and. any(abs(cryst%tnons(:,1)) > tol10) ) then
   MSG_ERROR('The first symmetry is not the identity operator!')
 end if

 do ikf=1,nkfull
   ik_ibz = bz2ibz(ikf,1)
   isym = bz2ibz(ikf,2)
   itimrev = bz2ibz(ikf,6)
   g0 = bz2ibz(ikf,3:5)        ! IS(k_ibz) + g0 = k_bz
   isirred_k = (isym == 1 .and. itimrev == 0 .and. all(g0 == 0))

   do spin=1,nsppol
     nband(ikf,spin) = inb%nband(ik_ibz+(spin-1)*inb%nkpt)
   end do

   if (force_istwfk1) then
     call get_kg(kfull(:,ikf),istwfk_1,ecut_eff,cryst%gmet,npw_k,kg_k)
     ABI_FREE(kg_k)
     istwfk(ikf) = 1
     npwarr(ikf) = npw_k
   else
     if (isirred_k) then
       istwfk(ikf) = inb%istwfk(ik_ibz)
       npwarr(ikf) = inb%npwarr(ik_ibz)
     else
       istw = set_istwfk(kfull(:,ikf))
       call get_kg(kfull(:,ikf),istw,ecut_eff,cryst%gmet,npw_k,kg_k)
       ABI_FREE(kg_k)
       istwfk(ikf) = istw
       npwarr(ikf) = npw_k
     end if
   end if
 end do

 ! Recostruct eig, occ and doccde in the BZ.
 bantot = sum(nband); mband = maxval(nband)

 ABI_MALLOC(doccde_3d, (mband, nkfull, nsppol))
 ABI_MALLOC(eig_3d, (mband, nkfull, nsppol))
 ABI_MALLOC(occ_3d, (mband, nkfull, nsppol))

 do spin=1,nsppol
   do ikf=1,nkfull
     ik_ibz = bz2ibz(ikf,1)
     doccde_3d(:,ikf,spin) = inb%doccde(:,ik_ibz,spin)
     eig_3d(:,ikf,spin) = inb%eig(:,ik_ibz,spin)
     occ_3d(:,ikf,spin) = inb%occ(:,ik_ibz,spin)
   end do
 end do

 ! Have to pack data to call ebands_init (I wonder who decided to use vectors!)
 ABI_MALLOC(doccde, (bantot))
 ABI_MALLOC(eig, (bantot))
 ABI_MALLOC(occ, (bantot))

 call pack_eneocc(nkfull,nsppol,mband,nband,bantot,doccde_3d,doccde)
 call pack_eneocc(nkfull,nsppol,mband,nband,bantot,eig_3d,eig)
 call pack_eneocc(nkfull,nsppol,mband,nband,bantot,occ_3d,occ)

 ABI_FREE(doccde_3d)
 ABI_FREE(eig_3d)
 ABI_FREE(occ_3d)

 call ebands_init(bantot,outb,inb%nelect,doccde,eig,istwfk,kfull,&
   nband,nkfull,npwarr,nsppol,inb%nspinor,inb%tphysel,inb%tsmear,inb%occopt,occ,wtk,&
   inb%charge, kptopt3, inb%kptrlatt_orig, inb%nshiftk_orig, inb%shiftk_orig, inb%kptrlatt, inb%nshiftk, inb%shiftk)

 ABI_FREE(istwfk)
 ABI_FREE(nband)
 ABI_FREE(npwarr)
 ABI_FREE(doccde)
 ABI_FREE(eig)
 ABI_FREE(occ)
 ABI_FREE(wtk)
 ABI_FREE(kfull)

end subroutine ebands_expandk
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_downsample
!! NAME
!! ebands_downsample
!!
!! FUNCTION
!!  Return a new ebands_t object of type ebands_t with a coarser IBZ contained in the inititial one.
!!
!! INPUTS
!!  cryst<crystal_t>=Info on unit cell and symmetries.
!!  in_kptrlatt(3,3)=Defines the sampling of the "small" IBZ. Must be submesh of the "fine" mesh.
!!  in_nshiftk= Number of shifts in the coarse k-mesh
!!  in_shiftk(3, in_nshiftk) = Shifts of the coarse k-mesh
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(ebands_t) function ebands_downsample(self, cryst, in_kptrlatt, in_nshiftk, in_shiftk) result(new)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: in_nshiftk
 type(ebands_t),intent(in) :: self
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: in_kptrlatt(3,3)
 real(dp),intent(in) :: in_shiftk(3, in_nshiftk)

!Local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1 = 1
 integer :: new_nkbz , timrev, bantot, new_nkibz, ik_ibz, ikf, spin,mband, comm
 real(dp) :: dksqmax
 character(len=500) :: msg
!arrays
 integer,allocatable :: ibz_c2f(:,:)
 integer :: new_kptrlatt(3,3)
 integer,allocatable :: istwfk(:),nband(:,:),npwarr(:)
 real(dp),allocatable :: new_kbz(:,:), new_wtk(:), new_kibz(:,:), doccde(:), eig(:), occ(:)
 real(dp),allocatable :: doccde_3d(:,:,:), eig_3d(:,:,:), occ_3d(:,:,:), new_shiftk(:,:)

! *********************************************************************

 comm = xmpi_comm_self

 ! Find IBZ associated to the new mesh.
 call kpts_ibz_from_kptrlatt(cryst, in_kptrlatt, self%kptopt, in_nshiftk, in_shiftk, &
   new_nkibz, new_kibz, new_wtk, new_nkbz, new_kbz, new_kptrlatt=new_kptrlatt, new_shiftk=new_shiftk)

 ! Costruct mapping IBZ_coarse --> IBZ_fine
 ! We don't change the value of nsppol hence sppoldbl1 is set to 1
 ABI_MALLOC(ibz_c2f, (new_nkibz*sppoldbl1, 6))

 timrev = kpts_timrev_from_kptopt(self%kptopt)
 call listkk(dksqmax, cryst%gmet, ibz_c2f, self%kptns, new_kibz, self%nkpt, new_nkibz, cryst%nsym, &
   sppoldbl1, cryst%symafm, cryst%symrel, timrev, comm, use_symrec=.False.)

 if (dksqmax > tol12) then
   write(msg, '(a,es16.6,6a)' )&
    "At least one of the k-points could not be generated from a symmetrical one. dksqmax: ",dksqmax, ch10,&
    "kptrlatt of input ebands: ",trim(ltoa(pack(self%kptrlatt, mask=.True.))),ch10, &
    "downsampled K-mesh: ",trim(ltoa(pack(in_kptrlatt, mask=.True.)))
   MSG_ERROR(msg)
 end if

 ABI_MALLOC(istwfk, (new_nkibz))
 ABI_MALLOC(nband, (new_nkibz, self%nsppol))
 ABI_MALLOC(npwarr, (new_nkibz))

 do ik_ibz=1,new_nkibz
   ikf = ibz_c2f(ik_ibz, 1)
   do spin=1,self%nsppol
     nband(ik_ibz, spin) = self%nband(ikf + (spin-1) * self%nkpt)
   end do
   istwfk(ik_ibz) = self%istwfk(ikf)
   npwarr(ik_ibz) = self%npwarr(ikf)
 end do

 ! Recostruct eig, occ and doccde in the new IBZ.
 bantot = sum(nband); mband = maxval(nband)

 ABI_MALLOC(doccde_3d, (mband, new_nkibz, self%nsppol))
 ABI_MALLOC(eig_3d, (mband, new_nkibz, self%nsppol))
 ABI_MALLOC(occ_3d, (mband, new_nkibz, self%nsppol))

 do spin=1,self%nsppol
   do ik_ibz=1,new_nkibz
     ikf = ibz_c2f(ik_ibz, 1)
     doccde_3d(:, ik_ibz, spin) = self%doccde(:, ikf, spin)
     eig_3d(:, ik_ibz, spin) = self%eig(:, ikf, spin)
     occ_3d(:, ik_ibz, spin) = self%occ(:, ikf, spin)
   end do
 end do

 ! Have to pack data to call ebands_init (I wonder who decided to use vectors!)
 ABI_MALLOC(doccde, (bantot))
 ABI_MALLOC(eig, (bantot))
 ABI_MALLOC(occ, (bantot))

 call pack_eneocc(new_nkibz, self%nsppol, mband, nband, bantot, doccde_3d, doccde)
 call pack_eneocc(new_nkibz, self%nsppol, mband, nband, bantot, eig_3d, eig)
 call pack_eneocc(new_nkibz, self%nsppol, mband, nband, bantot, occ_3d, occ)

 ABI_FREE(doccde_3d)
 ABI_FREE(eig_3d)
 ABI_FREE(occ_3d)

 call ebands_init(bantot, new, self%nelect, doccde, eig, istwfk, new_kibz, &
   nband, new_nkibz, npwarr, self%nsppol, self%nspinor, self%tphysel, self%tsmear, self%occopt, occ, new_wtk, &
   self%charge, self%kptopt, in_kptrlatt, in_nshiftk, self%shiftk, new_kptrlatt, size(new_shiftk, dim=2), new_shiftk)

 new%fermie = self%fermie

 ABI_FREE(istwfk)
 ABI_FREE(nband)
 ABI_FREE(npwarr)
 ABI_FREE(doccde)
 ABI_FREE(eig)
 ABI_FREE(occ)
 ABI_FREE(new_kibz)
 ABI_FREE(new_kbz)
 ABI_FREE(new_wtk)
 ABI_FREE(new_shiftk)
 ABI_FREE(ibz_c2f)

end function ebands_downsample
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_chop
!! NAME
!! ebands_chop
!!
!! FUNCTION
!!  Return a new ebands_t object with a selected number of bands between bstart and bstop
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(ebands_t) function ebands_chop(self, bstart, bstop) result(new)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(in)  :: self
 integer,intent(in) :: bstart, bstop

!Local variables ------------------------------
 integer :: mband, nkpt, nsppol

! *********************************************************************

 ABI_CHECK_IRANGE(bstart, 1, self%mband, "Invalid bstart")
 ABI_CHECK_IRANGE(bstop,  1, self%mband, "Invalid bstop")
 ABI_CHECK_ILEQ(bstart, bstop, "bstart should be <= bstop")

 ! First copy the bands
 call ebands_copy(self, new)

 ! Now chop them
 ABI_FREE(new%eig)
 ABI_FREE(new%occ)
 ABI_FREE(new%doccde)

 mband  = bstop - bstart + 1
 nkpt   = self%nkpt
 nsppol = self%nsppol

 ABI_MALLOC(new%eig, (mband, nkpt, nsppol))
 ABI_MALLOC(new%occ, (mband, nkpt, nsppol))
 ABI_MALLOC(new%doccde, (mband, nkpt, nsppol))

 new%mband  = mband
 new%nband  = mband
 new%eig    = self%eig(bstart:bstop,:,:)
 new%occ    = self%occ(bstart:bstop,:,:)
 new%doccde = self%doccde(bstart:bstop,:,:)

 new%bantot = sum(new%nband)

end function ebands_chop
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_sort
!! NAME
!! ebands_sort
!!
!! FUNCTION
!!  Sort eigvalues_k in ascending order and reorder arrays depending on nband_k
!!  Mainly used when interpolating band energies as the interpolator may not produce ordered eigenvalues
!!  and there are routines whose implementation assumes eig(b) <= eig(b+1)
!!
!! SIDE EFFECTS
!!  ebands<ebands_t> = Object with input energies sorted in output.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_sort(self)

!Arguments ------------------------------------
!scalars
 class(ebands_t),intent(inout) :: self

!Local variables-------------------------------
!scalars
 integer :: spin, ik_ibz, band, nband_k
!arrays
 integer :: iperm_k(self%mband)

! *********************************************************************

 do spin=1,self%nsppol
   do ik_ibz=1,self%nkpt
     nband_k = self%nband(ik_ibz + (spin - 1) * self%nkpt)
     iperm_k = [(band, band=1,nband_k)]
     call sort_dp(nband_k, self%eig(:, ik_ibz, spin), iperm_k, tol12)

     ! Shuffle other arrays depending on nband_k
     self%occ(1:nband_k, ik_ibz, spin) = self%occ(iperm_k(1:nband_k), ik_ibz, spin)
     self%doccde(1:nband_k, ik_ibz, spin) = self%doccde(iperm_k(1:nband_k), ik_ibz, spin)
     if (allocated(self%velocity)) then
       self%velocity(:, 1:nband_k, ik_ibz, spin) = self%velocity(:, iperm_k(1:nband_k), ik_ibz, spin)
     end if
   end do
 end do

end subroutine ebands_sort
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_interp_kmesh
!! NAME
!! ebands_interp_kmesh
!!
!! FUNCTION
!!  Interpolate energies on a k-mesh.
!!
!! INPUTS
!!  ebands<ebands_t> = Object with input energies.
!!  cryst<crystal_t> = Crystalline structure.
!!  params(:):
!!     params(0): interpolation type. 1 for star-functions
!!                                    2 for star-functions with group velocities
!!     if star-functions:
!!         params(2): Ratio between star functions and ab-initio k-points.
!!         params(3:4): Activate Fourier filtering (Eq 9 of PhysRevB.61.1639) if params(2) > tol6
!!         params(3)=rcut, params(4) = rsigma
!!  intp_kptrlatt(3,3) = New k-mesh
!!  intp_nshiftk= Number of shifts in new k-mesh.
!!  intp_shiftk(3,intp_nshiftk) = Shifts in new k-mesh.
!!  band_block(2)=Initial and final band index. If [0,0], all bands are used
!!    This is a global variable i.e. all MPI procs must call the routine with the same value.
!!  comm=MPI communicator
!!
!! OUTPUT
!!  New ebands_t object with interpolated energies.
!!
!! NOTES
!!  Fermi level and occupation factors of the interpolate bands are not recomputed by this routine.
!!  This operation is delegated to the caller.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE


type(ebands_t) function ebands_interp_kmesh(ebands, cryst, params, intp_kptrlatt, intp_nshiftk, intp_shiftk, &
        band_block, comm, out_prefix) result(new)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: intp_nshiftk,comm
 type(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: cryst
 character(len=*),optional,intent(in) :: out_prefix
!arrays
 integer,intent(in) :: intp_kptrlatt(3,3),band_block(2)
 real(dp),intent(in) :: params(:)
 real(dp),intent(in) :: intp_shiftk(3,intp_nshiftk)

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: ik_ibz,spin,new_bantot,new_mband,cplex,itype,nb,ib
 integer :: nprocs,my_rank,cnt,ierr,band,new_nkbz,new_nkibz,new_nshiftk
#ifdef HAVE_NETCDF
 integer :: ncid
#endif
 type(skw_t) :: skw
!arrays
 integer :: new_kptrlatt(3,3),my_bblock(2)
 integer,allocatable :: new_istwfk(:),new_nband(:,:),new_npwarr(:)
 real(dp),allocatable :: new_shiftk(:,:),new_kibz(:,:),new_kbz(:,:),new_wtk(:)
 real(dp),allocatable :: new_doccde(:),new_eig(:),new_occ(:)

! *********************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 itype = nint(params(1))
 my_bblock = band_block; if (all(band_block == 0)) my_bblock = [1, ebands%mband]
 nb = my_bblock(2) - my_bblock(1) + 1

 ! Get ibz, new shifts and new kptrlatt.
 call kpts_ibz_from_kptrlatt(cryst, intp_kptrlatt, ebands%kptopt, intp_nshiftk, intp_shiftk, &
   new_nkibz, new_kibz, new_wtk, new_nkbz, new_kbz, new_kptrlatt=new_kptrlatt, new_shiftk=new_shiftk)
 new_nshiftk = size(new_shiftk, dim=2)

 ! Initialize new ebands_t in new IBZ
 ABI_MALLOC(new_istwfk, (new_nkibz))
 new_istwfk = 1
 !do ik_ibz=1,new_nkibz
 !  new_istwfk(ik_ibz) = set_istwfk(new%kptns(:, ik_ibz))
 !end do
 ABI_MALLOC(new_nband, (new_nkibz, ebands%nsppol))
 new_nband = nb
 ABI_MALLOC(new_npwarr, (new_nkibz))
 new_npwarr = maxval(ebands%npwarr)
 new_bantot = sum(new_nband); new_mband = maxval(new_nband)
 ABI_CALLOC(new_doccde, (new_bantot))
 ABI_CALLOC(new_eig, (new_bantot))
 ABI_CALLOC(new_occ, (new_bantot))

 call ebands_init(new_bantot,new,ebands%nelect,new_doccde,new_eig,new_istwfk,new_kibz,&
   new_nband,new_nkibz,new_npwarr,ebands%nsppol,ebands%nspinor,ebands%tphysel,ebands%tsmear,&
   ebands%occopt,new_occ,new_wtk,&
   ebands%charge, ebands%kptopt, intp_kptrlatt, intp_nshiftk, intp_shiftk, new_kptrlatt, new_nshiftk, new_shiftk)

 ! Get fermi level from input ebands.
 new%fermie = ebands%fermie

 ABI_FREE(new_kibz)
 ABI_FREE(new_wtk)
 ABI_FREE(new_shiftk)
 ABI_FREE(new_kbz)
 ABI_FREE(new_istwfk)
 ABI_FREE(new_nband)
 ABI_FREE(new_npwarr)
 ABI_FREE(new_doccde)
 ABI_FREE(new_eig)
 ABI_FREE(new_occ)

 ! Build SKW object for all bands.
 if (itype == 1 .or. itype == 2) then
   cplex = 1; if (kpts_timrev_from_kptopt(ebands%kptopt) == 0) cplex = 2
   skw = skw_new(cryst, params(2:), cplex, ebands%mband, ebands%nkpt, ebands%nsppol, ebands%kptns, ebands%eig, my_bblock, comm)
   if (itype == 2) then
     ABI_CALLOC(new%velocity,(3,new%mband,new%nkpt,new%nsppol))
   end if
 else
   MSG_ERROR(sjoin("Wrong einterp params(1):", itoa(itype)))
 end if

 ! Interpolate eigenvalues and velocities.
 new%eig = zero; cnt = 0
 do spin=1,new%nsppol
   do ik_ibz=1,new%nkpt
     do ib=1,nb
       cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle  ! Mpi parallelism.
       ! Note the difference between band and ib index if band_block.
       band = my_bblock(1) + ib - 1
       select case (itype)
       case (1)
         call skw%eval_bks(band, new%kptns(:,ik_ibz), spin, new%eig(ib,ik_ibz,spin))
       case (2)
         call skw%eval_bks(band, new%kptns(:,ik_ibz), spin, new%eig(ib,ik_ibz,spin), new%velocity(:,ib,ik_ibz,spin))
       case default
         MSG_ERROR(sjoin("Wrong params(1):", itoa(itype)))
       end select
     end do
   end do
 end do
 call xmpi_sum(new%eig, comm, ierr)
 if (itype == 2) call xmpi_sum(new%velocity, comm, ierr)

 ! Sort eigvalues_k in ascending order to be compatible with other ebands routines.
 call ebands_sort(new)
 !call ebands_update_occ(new, dtset%spinmagntarget, prtvol=dtset%prtvol)

 if (my_rank == master .and. itype == 1 .and. present(out_prefix)) then
   ! Write ESKW file with crystal and (interpolated) band structure energies.
   !call wrtout(ab_out, sjoin("- Writing interpolated bands to file:", strcat(prefix, tag)))
#ifdef HAVE_NETCDF
   ! Write crystal and (interpolated) band structure energies.
   NCF_CHECK(nctk_open_create(ncid, strcat(out_prefix, "_ESKW.nc"), xmpi_comm_self))
   NCF_CHECK(cryst%ncwrite(ncid))
   NCF_CHECK(ebands_ncwrite(new, ncid))

   ! TODO
   !NCF_CHECK(skw%ncwrite(ncid))
   ! Define variables specific to SKW algo.
   !ncerr = nctk_def_arrays(ncid, [ &
   ! nctkarr_t("band_block", "int", "two"), &
   ! nctkarr_t("einterp", "dp", "four")], defmode=.True.)
   !NCF_CHECK(ncerr)

   ! Write data.
   !NCF_CHECK(nctk_set_datamode(ncid))
   !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "band_block"), band_block))
   !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "einterp"), params))
   NCF_CHECK(nf90_close(ncid))
#endif
 end if

 call skw%free()

end function ebands_interp_kmesh
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_interp_kpath
!! NAME
!! ebands_interp_kpath
!!
!! FUNCTION
!!  Interpolate energies on a k-path
!!
!! INPUTS
!!  ebands<ebands_t> = Object with input energies.
!!  cryst<crystal_t> = Crystalline structure.
!!  kpath<kpath_t> = Object describing the k-path
!!  params(:):
!!    params(1): 1 for SKW.
!!  band_block(2)=Initial and final band index to be interpolated. [0,0] if all bands are used.
!!    This is a global variable i.e. all MPI procs must call the routine with the same value.
!!  comm=MPI communicator
!!
!! OUTPUT
!!  New ebands_t object with interpolated energies.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(ebands_t) function ebands_interp_kpath(ebands, cryst, kpath, params, band_block, comm) result(new)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: cryst
 type(kpath_t),intent(in) :: kpath
!arrays
 integer,intent(in) :: band_block(2)
 real(dp),intent(in) :: params(:)

!Local variables-------------------------------
!scalars
 integer,parameter :: new_nshiftk=1
 integer :: ik_ibz,spin,new_bantot,new_mband,cplex
 integer :: nprocs,my_rank,cnt,ierr,band,new_nkibz,itype,nb,ib, new_kptopt
 type(skw_t) :: skw
!arrays
 integer,parameter :: new_kptrlatt(3,3)=0
 integer :: my_bblock(2)
 integer,allocatable :: new_istwfk(:),new_nband(:,:),new_npwarr(:)
 real(dp),parameter :: new_shiftk(3,1) = zero
 real(dp),allocatable :: new_wtk(:),new_doccde(:),new_eig(:),new_occ(:)

! *********************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 itype = nint(params(1))
 my_bblock = band_block; if (all(band_block == 0)) my_bblock = [1, ebands%mband]
 nb = my_bblock(2) - my_bblock(1) + 1

 if (ebands%nkpt == 1) then
   MSG_WARNING("Cannot interpolate band energies when nkpt = 1. Returning")
   return
 end if

 ! Initialize new ebands_t.
 new_nkibz = kpath%npts
 ABI_MALLOC(new_istwfk, (new_nkibz))
 new_istwfk = 1
 ABI_MALLOC(new_nband, (new_nkibz, ebands%nsppol))
 new_nband = nb
 ABI_MALLOC(new_npwarr, (new_nkibz))
 new_npwarr = maxval(ebands%npwarr)
 new_bantot = sum(new_nband); new_mband = maxval(new_nband)
 ABI_CALLOC(new_eig, (new_bantot))
 ABI_CALLOC(new_doccde, (new_bantot))
 ABI_CALLOC(new_occ, (new_bantot))
 ABI_CALLOC(new_wtk, (new_nkibz))

 ! Needed by AbiPy to understand that we have a k-path instead of a mesh.
 new_kptopt = -kpath%nbounds

 call ebands_init(new_bantot,new,ebands%nelect,new_doccde,new_eig,new_istwfk,kpath%points,&
   new_nband,new_nkibz,new_npwarr,ebands%nsppol,ebands%nspinor,ebands%tphysel,ebands%tsmear,&
   ebands%occopt,new_occ,new_wtk,&
   ebands%charge, new_kptopt, new_kptrlatt, new_nshiftk, new_shiftk, new_kptrlatt, new_nshiftk, new_shiftk)
 new%fermie = ebands%fermie

 ABI_FREE(new_wtk)
 ABI_FREE(new_istwfk)
 ABI_FREE(new_nband)
 ABI_FREE(new_npwarr)
 ABI_FREE(new_doccde)
 ABI_FREE(new_eig)
 ABI_FREE(new_occ)

 ! Build SKW object for all bands.
 select case (itype)
 case (1)
   cplex = 1; if (kpts_timrev_from_kptopt(ebands%kptopt) == 0) cplex = 2
   skw = skw_new(cryst, params(2:), cplex, ebands%mband, ebands%nkpt, ebands%nsppol, ebands%kptns, ebands%eig, &
                 my_bblock, comm)

 case default
   MSG_ERROR(sjoin("Wrong einterp params(1):", itoa(itype)))
 end select

 ! Interpolate eigenvalues.
 new%eig = zero; cnt = 0
 do spin=1,new%nsppol
   do ik_ibz=1,new%nkpt
     do ib=1,nb
       cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle  ! Mpi parallelism.
       ! Note the difference between band and ib index if band_block.
       band = my_bblock(1) + ib - 1
       select case (itype)
       case (1)
         call skw%eval_bks(band, new%kptns(:,ik_ibz), spin, new%eig(ib,ik_ibz,spin))
       case default
         MSG_ERROR(sjoin("Wrong einterp params(1):", itoa(itype)))
       end select
     end do
   end do
 end do
 call xmpi_sum(new%eig, comm, ierr)

 ! Sort eigvalues_k in ascending order to be compatible with other ebands routines.
 call ebands_sort(new)

 call skw%free()

end function ebands_interp_kpath
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_get_edos_matrix_elements
!! NAME
!!  ebands_get_edos_matrix_elements
!!
!! FUNCTION
!!  Compute electron DOS and weighted e-DOS with weights given by precomputed scalar, vectorial
!!  or tensorial matrix elements.
!!  Weights are provided in input as (..., num_entries, bsize, nkpt, nsppol) tables, see below.
!!
!! INPUTS
!!  ebands<ebands_t>=Band structure object.
!!  cryst<cryst_t>=Info on the crystalline structure.
!!  bsize=Number of bands in bks_vals, bks_vecs and bks_tens
!!    Not necessarily equal to ebands%mband when brange is used.
!!  nvals=Number of scalar entries. Maybe zero
!!  bks_vals=Scalar matrix elements
!!  nvecs=Number of 3d-vectorial entries. Maybe zero
!!  bks_vecs=Vectorial matrix elements in Cartesian Coordinates
!!  ntens=Numer of 3x3 tensorial entries in Cartesian coordinates. Maybe zero
!!  bks_tens= Tensorial matrix elements (3x3) in Cartesian Coordinates
!!  intmeth=
!!    1 for Gaussian,
!!    2 or -2 for tetrahedra (-2 if Blochl corrections must be included).
!!    If nkpt == 1 (Gamma only), the routine fallbacks to the Gaussian method.
!!  step=Step on the linear mesh in Ha. If < 0, the routine will use the mean of the energy level spacing
!!  broad=Gaussian broadening, If <0, the routine will use a default
!!    value for the broadening computed from the mean of the energy level spacing.
!!    No meaning for tetrahedra
!!  comm=MPI communicator
!!  [brange(2)]=Minimum and maximum band index. Default if not present is `full band range`.
!!    If given bsize must be equal: to brange(2) - brange(1) + 1
!!  [erange(2)]=Minimum and maximum energy to be considered. Default if not present is `full energy range`.
!!
!! OUTPUT
!!  out_valsdos: (nw, 2, nvals, nsppol) array with DOS for scalar quantities if nvals > 0
!!  out_vecsdos: (nw, 2, 3, nvecs, nsppol)) array with DOS weighted by vectorial terms if nvecs > 0
!!  out_tensdos: (nw, 2,3, 3, ntens,  nsppol) array with DOS weighted by tensorial terms if ntens > 0
!!
!!   All these arrays are allocated by the routine. The number of points is available in edos%nw.
!!   (nw, 1, ...) stores the weighted DOS (w-DOS)
!!   (nw, 2, ...) stores the integrated w-DOS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(edos_t) function ebands_get_edos_matrix_elements(ebands, cryst, bsize, &
                                                      nvals, bks_vals, nvecs, bks_vecs, ntens, bks_tens, &
                                                      intmeth, step, broad, out_valsdos, out_vecsdos, out_tensdos, comm, &
                                                      brange, erange) result(edos)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bsize, nvals, nvecs, ntens, intmeth, comm
 real(dp),intent(in) :: step, broad
 class(ebands_t),intent(in)  :: ebands
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,optional,intent(in) :: brange(2)
 real(dp),optional,intent(in) :: erange(2)
 real(dp),intent(in) :: bks_vals(nvals, bsize, ebands%nkpt, ebands%nsppol)
 real(dp),intent(in) :: bks_vecs(3, nvecs, bsize, ebands%nkpt, ebands%nsppol)
 real(dp),intent(in) :: bks_tens(3, 3, ntens, bsize, ebands%nkpt, ebands%nsppol)
 real(dp),allocatable,intent(out) :: out_valsdos(:,:,:,:), out_vecsdos(:,:,:,:,:), out_tensdos(:,:,:,:,:,:)

!Local variables-------------------------------
!scalars
 integer :: nproc, my_rank, nw, spin, band, ib, ik_ibz, cnt, idat, ierr, bcorr
 integer :: ii, jj, ief, bmin_, bmax_
 real(dp),parameter :: max_occ1 = one
 real(dp) :: emax, emin, wtk, max_occ
 real(dp) :: cpu, wall, gflops
 logical :: check_erange
 character(len=500) :: msg
 type(htetra_t) :: tetra
!arrays
 real(dp) :: eminmax_spin(2,ebands%nsppol)
 real(dp) :: vsum(3), tsum(3,3)
 real(dp),allocatable :: wme0(:),tmp_eigen(:), weights(:,:)

! *********************************************************************

 call cwtime(cpu, wall, gflops, "start")

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 edos%nkibz = ebands%nkpt; edos%nsppol = ebands%nsppol
 edos%intmeth = intmeth
 if (ebands%nkpt == 1) then
   MSG_COMMENT("Cannot use tetrahedra for e-DOS when nkpt == 1. Switching to gaussian method")
   edos%intmeth = 1
 end if

 edos%broad = broad; edos%step = step

 ! Define band range.
 bmin_ = 1; bmax_ = ebands%mband
 if (present(brange)) then
   bmin_ = brange(1); bmax_ = brange(2)
 end if

 ABI_CHECK_IRANGE(bmin_, 1, ebands%mband, "Wrong bmin:")
 ABI_CHECK_IRANGE(bmax_, bmin_, ebands%mband, "Wrong bmax:")
 ABI_CHECK_IRANGE(bsize, 1, ebands%mband, "Wrong bsize:")

 if (present(erange)) then
   ! use optional args if provided.
   emin = erange(1)
   emax = erange(2)
   check_erange = .True.
 else
   ! Compute the linear mesh so that it encloses all bands.
   eminmax_spin = ebands_get_minmax(ebands, "eig")
   emin = minval(eminmax_spin(1, :)); emin = emin - 0.1_dp * abs(emin)
   emax = maxval(eminmax_spin(2, :)); emax = emax + 0.1_dp * abs(emax)
   check_erange = .False.
 end if

 nw = nint((emax - emin) / edos%step) + 1
 edos%nw = nw

 ABI_MALLOC(edos%mesh, (nw))
 edos%mesh = arth(emin, edos%step, nw)

 ABI_CALLOC(edos%gef, (0:edos%nsppol))
 ABI_CALLOC(edos%dos,  (nw, 0:edos%nsppol))
 ABI_CALLOC(edos%idos, (nw, 0:edos%nsppol))

 ! Allocate output arrays depending on input.
 if (nvals > 0) then
   ABI_CALLOC(out_valsdos, (nw, 2, nvals, ebands%nsppol))
 endif
 if (nvecs > 0) then
   ABI_CALLOC(out_vecsdos, (nw, 2, 3, nvecs, ebands%nsppol))
 end if
 if (ntens > 0) then
   ABI_CALLOC(out_tensdos, (nw, 2, 3, 3, ntens, ebands%nsppol))
 end if

 !call wrtout(std_out, " Computing DOS weighted by matrix elements.")
 select case (intmeth)
 case (1)
   ! Gaussian
   ABI_MALLOC(wme0, (nw))
   cnt = 0
   do spin=1,ebands%nsppol
     do ik_ibz=1,ebands%nkpt
       cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle  ! MPI parallelism
       wtk = ebands%wtk(ik_ibz)
       do band=bmin_,bmax_
         ib = band - bmin_ + 1

         if (check_erange) then
           if (ebands%eig(band, ik_ibz, spin) < emin - five * broad) cycle
           if (ebands%eig(band, ik_ibz, spin) > emax + five * broad) cycle
         end if

         wme0 = edos%mesh - ebands%eig(band, ik_ibz, spin)
         wme0 = gaussian(wme0, broad) * wtk
         edos%dos(:,spin) = edos%dos(:,spin) + wme0(:)

         ! scalars
         do idat=1,nvals
           out_valsdos(:, 1, idat, spin) = out_valsdos(:,1, idat, spin) + wme0(:) * bks_vals(idat, ib, ik_ibz, spin)
           ! FIXME: This is quite inefficient! Integration should be performed outside!
           call simpson_int(nw, step, out_valsdos(:,1, idat, spin), out_valsdos(:,2,idat,spin))
         end do

         ! vectors
         do idat=1,nvecs
           ! get components, symmetrize and accumulate.
           vsum = cryst%symmetrize_cart_vec3(bks_vecs(:, idat, ib, ik_ibz, spin))
           do ii=1,3
             out_vecsdos(:, 1, ii, idat, spin) = out_vecsdos(:, 1, ii, idat, spin) + wme0(:) * vsum(ii)
             call simpson_int(nw, step, out_vecsdos(:,1,ii,idat,spin), out_vecsdos(:,2,ii,idat,spin))
           end do
         end do

         ! tensor
         do idat=1,ntens
           ! get components, symmetrize and accumulate.
           tsum = cryst%symmetrize_cart_tens33(bks_tens(:, :, idat, ib, ik_ibz, spin))
           do ii=1,3
             do jj=1,3
               out_tensdos(:,1,jj,ii,idat,spin) = out_tensdos(:,1,jj,ii,idat,spin) + wme0(:) * tsum(jj,ii)
               call simpson_int(nw, step, out_tensdos(:,1,jj,ii,idat,spin), out_tensdos(:,2,jj,ii,idat,spin))
             end do
           end do
         end do

       end do !band
     end do !ik_ibz
   end do !spin

   ABI_FREE(wme0)
   call xmpi_sum(edos%dos, comm, ierr)
   if (nvals > 0) call xmpi_sum(out_valsdos, comm, ierr)
   if (nvecs > 0) call xmpi_sum(out_vecsdos, comm, ierr)
   if (ntens > 0) call xmpi_sum(out_tensdos, comm, ierr)

 case (2, -2)
   ! Consistency test
   ABI_CHECK(all(ebands%nband == ebands%nband(1)), 'For tetrahedra, nband(:) must be constant')

   ! Build tetra object.
   tetra = tetra_from_kptrlatt(cryst, ebands%kptopt, ebands%kptrlatt, &
     ebands%nshiftk, ebands%shiftk, ebands%nkpt, ebands%kptns, comm, msg, ierr)
   ABI_CHECK(ierr == 0, msg)

   ! For each spin and band, interpolate over kpoints,
   ! calculate integration weights and DOS contribution.
   ABI_MALLOC(tmp_eigen, (ebands%nkpt))
   ABI_MALLOC(weights, (nw, 2))

   ! Blochl's corrections?
   bcorr = 0; if (intmeth == -2) bcorr = 1

   cnt = 0
   do spin=1,ebands%nsppol
     do band=bmin_,bmax_
       ! For each band get its contribution
       tmp_eigen = ebands%eig(band,:,spin)
       ib = band - bmin_ + 1

       if (check_erange) then
         if (all(tmp_eigen < emin)) cycle
         if (all(tmp_eigen > emax)) cycle
       end if

       do ik_ibz=1,ebands%nkpt
         cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism

         call tetra%get_onewk(ik_ibz, bcorr, nw, ebands%nkpt, tmp_eigen, emin, emax, max_occ1, weights)
         weights = weights * ebands%wtk(ik_ibz)

         ! Compute DOS
         edos%dos(:,spin) = edos%dos(:,spin) + weights(:, 1)

         ! scalar
!$OMP PARALLEL DO
         do idat=1,nvals
           out_valsdos(:, :, idat, spin) = out_valsdos(:, :, idat, spin) + weights(:, :) * bks_vals(idat, ib, ik_ibz, spin)
         end do

         ! vector
!$OMP PARALLEL DO PRIVATE(vsum)
         do idat=1,nvecs
           ! get components, symmetrize and accumulate.
           vsum = cryst%symmetrize_cart_vec3(bks_vecs(:, idat, ib, ik_ibz, spin))
           do ii=1,3
             out_vecsdos(:, :, ii, idat, spin) = out_vecsdos(:, :, ii, idat, spin) + weights(:, :) * vsum(ii)
           end do
         end do

         ! tensor
!$OMP PARALLEL DO PRIVATE(tsum)
         do idat=1,ntens
           ! get components, symmetrize and accumulate.
           tsum = cryst%symmetrize_cart_tens33(bks_tens(:, :, idat, ib, ik_ibz, spin))
           do ii=1,3
             do jj=1,3
               out_tensdos(:, :, jj, ii, idat, spin) = out_tensdos(:, :, jj, ii, idat, spin) + weights(:, :) * tsum(jj,ii)
             end do
           end do
         end do

       end do ! ik_ibz
     end do ! band
   end do ! spin

   ! Free memory
   ABI_FREE(weights)
   ABI_FREE(tmp_eigen)
   call tetra%free()

   call xmpi_sum(edos%dos, comm, ierr)
   if (nvals > 0) call xmpi_sum(out_valsdos, comm, ierr)
   if (nvecs > 0) call xmpi_sum(out_vecsdos, comm, ierr)
   if (ntens > 0) call xmpi_sum(out_tensdos, comm, ierr)

 case default
   MSG_ERROR(sjoin("Wrong integration method:", itoa(intmeth)))
 end select

 ! Compute total DOS and IDOS
 max_occ = two / (ebands%nspinor * ebands%nsppol)
 edos%dos(:, 0) = max_occ * sum(edos%dos(:,1:), dim=2)

 do spin=1,edos%nsppol
   call simpson_int(nw, edos%step, edos%dos(:,spin), edos%idos(:,spin))
 end do
 edos%idos(:, 0) = max_occ * sum(edos%idos(:,1:), dim=2)

 ! Use bisection to find the Fermi level.
 ! Warning: this code assumes idos[i+1] >= idos[i]. This condition may not be
 ! fullfilled if we use tetra and this is the reason why we have filtered the DOS.
 ief = bisect(edos%idos(:,0), ebands%nelect)

 ! Handle out of range condition.
 if (ief == 0 .or. ief == nw) then
   write(msg,"(a, f14.2, 4a)") &
    "Bisection could not find an initial guess for the Fermi level with nelect:",ebands%nelect, ch10, &
    "Possible reasons: not enough bands in DOS or wrong number of electrons.", ch10, &
    "Returning from ebands_get_edos_matrix_elements without setting edos%ief !"
   MSG_WARNING(msg)
   return
 end if

 ! TODO: Use linear interpolation to find an improved estimate of the Fermi level?
 edos%ief = ief
 do spin=0,edos%nsppol
   edos%gef(spin) = edos%dos(ief,spin)
 end do

 call cwtime_report(" ebands_get_edos_matrix_elements", cpu, wall, gflops)

end function ebands_get_edos_matrix_elements
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_get_jdos
!! NAME
!! ebands_get_jdos
!!
!! FUNCTION
!!  Compute the joint density of states.
!!
!! INPUTS
!!  ebands<ebands_t>=Band structure object.
!!  cryst<cryst_t>=Info on the crystalline structure.
!!  intmeth= 1 for gaussian, 2 or -2 for tetrahedra (-2 if Blochl corrections must be included).
!!  step=Step on the linear mesh in Ha. If <0, the routine will use the mean of the energy level spacing
!!  broad=Gaussian broadening, If <0, the routine will use a default
!!    value for the broadening computed from the mean of the energy level spacing.
!!    No meaning if tetra method
!!  comm=MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(jdos_t) function ebands_get_jdos(ebands, cryst, intmeth, step, broad, comm, ierr) result (jdos)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: intmeth,comm
 integer,intent(out) :: ierr
 real(dp),intent(in) :: step,broad
 class(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: cryst

!Local variables-------------------------------
!scalars
 integer :: ik_ibz,ibc,ibv,spin,nw,nband_k,nbv,nproc,my_rank,cnt,mpierr,bcorr !iw, unt,
 real(dp) :: wtk,wmax
 type(htetra_t) :: tetra
 character(len=500) :: msg
 !character(len=fnlen) :: path
!arrays
 integer :: val_idx(ebands%nkpt,ebands%nsppol)
 real(dp) :: eminmax(2,ebands%nsppol)
 real(dp),allocatable :: cvmw(:),wdt(:,:)

! *********************************************************************

 ierr = 0
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 jdos%nsppol = ebands%nsppol
 jdos%nkibz = ebands%nkpt

 ! Find the valence band index for each k and spin.
 val_idx = ebands_get_valence_idx(ebands)

 do spin=1,ebands%nsppol
   if (any(val_idx(:,spin) /= val_idx(1,spin))) then
     write(msg,'(a,i0,a)')&
     'Trying to compute JDOS with a metallic band structure for spin: ',spin,&
     'Assuming you know what you are doing, continuing anyway! '
     MSG_COMMENT(msg)
   end if
 end do

 ! Compute the linear mesh so that it encloses all bands.
 !if (.not. present(mesh)) then
 eminmax = ebands_get_minmax(ebands, "eig")
 wmax = maxval(eminmax(2,:) - eminmax(1,:))
 nw = nint(wmax/step) + 1
 ABI_MALLOC(jdos%mesh, (nw))
 jdos%mesh = arth(zero, step, nw)

 jdos%nw = nw
 jdos%intmeth = intmeth
 jdos%broad = broad

 !if (ebands%nkpt == 1) then
 !  MSG_COMMENT("Cannot use tetrahedra for e-DOS when nkpt == 1. Switching to gaussian method")
 !  jdos%intmeth = 1
 !end if

 ABI_CALLOC(jdos%values, (nw, ebands%nsppol))

 select case (intmeth)
 case (1)
   ! Gaussian
   ABI_MALLOC(cvmw, (nw))

   cnt = 0
   do spin=1,ebands%nsppol
     do ik_ibz=1,ebands%nkpt
       wtk = ebands%wtk(ik_ibz)
       nband_k = ebands%nband(ik_ibz + (spin-1)*ebands%nkpt)
       nbv = val_idx(ik_ibz, spin)

       do ibv=1,nbv
         cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
         do ibc=nbv+1,nband_k
           cvmw = ebands%eig(ibc,ik_ibz,spin) - ebands%eig(ibv,ik_ibz,spin) - jdos%mesh
           jdos%values(:, spin) = jdos%values(:, spin) + wtk * gaussian(cvmw, broad)
         end do
       end do

     end do ! ik_ibz
   end do ! spin

   ABI_FREE(cvmw)
   call xmpi_sum(jdos%values, comm, mpierr)

 case (2, -2)
   ! Tetrahedron method
   if (any(ebands%nband /= ebands%nband(1)) ) then
     MSG_WARNING('For tetrahedra, nband(:) must be constant')
     ierr = ierr + 1
   end if
   if (ierr/=0) return

   tetra = tetra_from_kptrlatt(cryst, ebands%kptopt, ebands%kptrlatt, &
     ebands%nshiftk, ebands%shiftk, ebands%nkpt, ebands%kptns, comm, msg, ierr)
   if (ierr /= 0) then
     call tetra%free(); return
   end if

   ! For each spin and band, interpolate over kpoints,
   ! calculate integration weights and DOS contribution.
   ABI_MALLOC(cvmw, (jdos%nkibz))
   ABI_MALLOC(wdt, (nw, 2))

   bcorr = 0; if (intmeth == -2) bcorr = 1
   cnt = 0
   do spin=1,ebands%nsppol
     nbv = val_idx(1, spin)
     do ibv=1,nbv
       do ibc=nbv+1,ebands%mband
         ! For each (c, v) get its contribution
         cvmw = ebands%eig(ibc,:,spin) - ebands%eig(ibv,:,spin)
         do ik_ibz=1,ebands%nkpt
           cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle  ! mpi-parallelism

           ! Calculate integration weights at each irred k-point (Blochl et al PRB 49 16223 [[cite:Bloechl1994a]])
           call tetra%get_onewk(ik_ibz, bcorr, nw, ebands%nkpt, cvmw, jdos%mesh(0), jdos%mesh(nw), one, wdt)
           jdos%values(:,spin) = jdos%values(:,spin) + wdt(:, 1) * ebands%wtk(ik_ibz)
         end do
       end do ! ibc
     end do ! ibv
   end do ! spin

   call xmpi_sum(jdos%values, comm, mpierr)

   ! Free memory
   ABI_FREE(wdt)
   ABI_FREE(cvmw)
   call tetra%free()

 case default
   MSG_ERROR(sjoin("Wrong integration method:", itoa(intmeth)))
 end select

 if (ebands%nsppol == 1) then
   jdos%values(0,:) = two * jdos%values(1,:)
 else
   jdos%values(0,:) = sum(jdos%values(1:2, :), dim=2)
 end if

 ! Write data.
 !if (my_rank == 0) then
 !  path = "jdos_gauss.data"; if (intmeth == 2) path = "jdos_tetra.data"
 !  if (open_file(path, msg, newunit=unt, form="formatted", action="write") /= 0) then
 !    MSG_ERROR(msg)
 !  end if
 !  do iw=1,nw
 !    write(unt,*)jdos%mesh(iw),(jdos(iw,spin), spin=1,ebands%nsppol)
 !  end do
 !  close(unt)
 !end if

end function ebands_get_jdos
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/jdos_ncwrite
!! NAME
!!  jdos_ncwrite
!!
!! FUNCTION
!!  Write JDOS to netcdf file.
!!
!! INPUTS
!!  ncid=NC file handle.
!!  [prefix]=String prepended to netcdf dimensions/variables (HDF5 poor-man groups)
!!   Empty string if not specified.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function jdos_ncwrite(jdos, ncid, prefix) result(ncerr)

!Arguments ------------------------------------
!scalars
 class(jdos_t),intent(inout)  :: jdos
 integer,intent(in) :: ncid
 character(len=*),optional,intent(in) :: prefix

!Local variables-------------------------------
 character(len=500) :: prefix_

! *********************************************************************

 prefix_ = ""; if (present(prefix)) prefix_ = trim(prefix)

#ifdef HAVE_NETCDF
 ! Define dimensions.
 ncerr = nctk_def_dims(ncid, [ &
   nctkdim_t("nsppol_plus1", jdos%nsppol + 1), nctkdim_t("jdos_nw", jdos%nw)], defmode=.True., prefix=prefix_)
 NCF_CHECK(ncerr)

 ! Define variables
 NCF_CHECK(nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "jdos_intmeth", "jdos_nkibz"], prefix=prefix_))
 NCF_CHECK(nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "jdos_broad"], prefix=prefix_))

 ncerr = nctk_def_arrays(ncid, [ &
   nctkarr_t("jdos_mesh", "dp", "jdos_nw"), &
   nctkarr_t("jdos_values", "dp", "jdos_nw, nsppol_plus1") &
 ],  prefix=prefix_)
 NCF_CHECK(ncerr)

 ! Write data.
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("jdos_intmeth")), jdos%intmeth))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("jdos_nkibz")), jdos%nkibz))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("jdos_broad")), jdos%broad))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("jdos_mesh")), jdos%mesh))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, pre("edos_values")), jdos%values))

#else
 MSG_ERROR("netcdf library not available")
#endif

contains
  pure function pre(istr) result(ostr)
    character(len=*),intent(in) :: istr
    character(len=len_trim(prefix_) + len_trim(istr)+1) :: ostr
    ostr = trim(prefix_) // trim(istr)
  end function pre

end function jdos_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/jdos_free
!! NAME
!!  jdos_free
!!
!! FUNCTION
!!  Free memory
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine jdos_free(jdos)

!Arguments ------------------------------------
 class(jdos_t),intent(inout)  :: jdos

! *********************************************************************

 ABI_FREE(jdos%mesh)
 ABI_FREE(jdos%values)

end subroutine jdos_free
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_prtbltztrp
!! NAME
!! ebands_prtbltztrp
!!
!! FUNCTION
!!   Output files for BoltzTraP code, which integrates Boltzmann transport quantities
!!   over the Fermi surface for different T and chemical potentials. Abinit provides
!!   all necessary input files: struct, energy, input file, and def file for the unit
!!   definitions of fortran files in BT.
!!   See http://www.icams.de/content/departments/ams/madsen/boltztrap.html
!!
!! INPUTS
!!  ebands<ebands_t>=Band structure object.
!!  cryst<cryst_t>=Info on the crystalline structure.
!!  fname_radix = radix of file names for output
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! PARENTS
!!      eph,outscfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_prtbltztrp(ebands, crystal, fname_radix, tau_k)

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: crystal
 character(len=fnlen), intent(in) :: fname_radix
!arrays
 real(dp), intent(in), optional :: tau_k(ebands%nsppol,ebands%nkpt,ebands%mband)

!Local variables-------------------------------
!scalars
 integer :: iout, isym, iband, isppol, ikpt, nsppol, nband
 real(dp),parameter :: ha2ryd=two
 real(dp) :: ewindow
 character(len=fnlen) :: filename
 character(len=2) :: so_suffix
 character(len=500) :: msg
!arrays
 real(dp) :: nelec(ebands%nsppol)
 character(len=3) :: spinsuffix(ebands%nsppol)

! *************************************************************************

 !MG FIXME The number of electrons is wrong if the file is produced in a NSCF run.
 ! See https://forum.abinit.org/viewtopic.php?f=19&t=3339

 nelec = ebands_nelect_per_spin(ebands)
 nsppol = ebands%nsppol
 nband = ebands%nband(1)

 so_suffix=""
 if (nsppol > 1 .or. ebands%nspinor > 1) so_suffix="so"

 if (nsppol == 1) then
   spinsuffix(1) = "ns_"
 else
   spinsuffix = ["up_", "dn_"]
 end if

 do isppol=1,nsppol

   !input file for boltztrap: general info, Ef, Nelec, etc...
   filename= trim(fname_radix)//"_"//trim(spinsuffix(isppol))//"BLZTRP.intrans"
   if (open_file(filename, msg, newunit=iout, form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   ewindow = 1.1_dp * (ebands%fermie-minval(ebands%eig(1, :, isppol)))
   write (iout, '(a)') "GENE                      # Format of input: generic format, with Symmetries"
   write (iout, '(a)') "0 0 0 0.0                 # iskip (not presently used) idebug setgap shiftgap"
   write (iout, '(E15.5,a,2F10.4,a)') ebands%fermie*ha2ryd, " 0.0005 ", ewindow*ha2ryd, nelec(isppol), &
&   "  # Fermilevel (Ry), energy grid spacing, energy span around Fermilevel, number of electrons for this spin"
   write (iout, '(a)') "CALC                      # CALC (calculate expansion coeff), NOCALC read from file"
   write (iout, '(a)') "3                         # lpfac, number of latt-points per k-point"
   write (iout, '(a)') "BOLTZ                     # run mode (only BOLTZ is supported)"
   write (iout, '(a)') ".15                       # (efcut) energy range of chemical potential"
   write (iout, '(a)') "300. 10.                  # Tmax, temperature grid spacing"
   write (iout, '(2a)') "-1                        # energyrange of bands given ",&
&   "individual DOS output sig_xxx and dos_xxx (xxx is band number)"
   write (iout, '(a)') "HISTO                     # DOS calculation method. Other possibility is TETRA"
   write (iout, '(a)') "No                        # not using model for relaxation time"
   write (iout, '(a)') "3                         # Number of doping levels coefficients will be output for"
   write (iout, '(a)') "-1.e16 0.0d0 1.e16        # Values of doping levels (in carriers / cm^3"
   close(iout)

!files file, with association of all units for Boltztrap
   filename= trim(fname_radix)//"_"//trim(spinsuffix(isppol))//"BLZTRP.def"
   if (open_file(filename, msg, newunit=iout, form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   write (iout, '(3a)') "5, '", trim(fname_radix)//"_"//trim(spinsuffix(isppol))//"BLZTRP.intrans',      'old',    'formatted',0"
   write (iout, '(3a)') "6, '", trim(fname_radix)//"_BLZTRP", ".outputtrans',      'unknown',    'formatted',0"
   write (iout, '(3a)') "20,'", trim(fname_radix)//"_BLZTRP", ".struct',         'old',    'formatted',0"
   write (iout, '(3a)') "10,'", trim(fname_radix)//"_BLZTRP."//trim(spinsuffix(isppol))//"energy"//trim(so_suffix),&
&   "',         'old',    'formatted',0"
   if (present (tau_k)) then
     write (iout, '(3a)') "11,'", trim(fname_radix)//"_BLZTRP", ".tau_k',         'old',    'formatted',0"
   end if
   write (iout, '(3a)') "48,'", trim(fname_radix)//"_BLZTRP", ".engre',         'unknown',    'unformatted',0"
   write (iout, '(3a)') "49,'", trim(fname_radix)//"_BLZTRP", ".transdos',        'unknown',    'formatted',0"
   write (iout, '(3a)') "50,'", trim(fname_radix)//"_BLZTRP", ".sigxx',        'unknown',    'formatted',0"
   write (iout, '(3a)') "51,'", trim(fname_radix)//"_BLZTRP", ".sigxxx',        'unknown',    'formatted',0"
   write (iout, '(3a)') "21,'", trim(fname_radix)//"_BLZTRP", ".trace',           'unknown',    'formatted',0"
   write (iout, '(3a)') "22,'", trim(fname_radix)//"_BLZTRP", ".condtens',           'unknown',    'formatted',0"
   write (iout, '(3a)') "24,'", trim(fname_radix)//"_BLZTRP", ".halltens',           'unknown',    'formatted',0"
   write (iout, '(3a)') "25,'", trim(fname_radix)//"_BLZTRP", ".trace_fixdoping',     'unknown',    'formatted',0"
   write (iout, '(3a)') "26,'", trim(fname_radix)//"_BLZTRP", ".condtens_fixdoping',           'unknown',    'formatted',0"
   write (iout, '(3a)') "27,'", trim(fname_radix)//"_BLZTRP", ".halltens_fixdoping',           'unknown',    'formatted',0"
   write (iout, '(3a)') "30,'", trim(fname_radix)//"_BLZTRP", "_BZ.dx',           'unknown',    'formatted',0"
   write (iout, '(3a)') "31,'", trim(fname_radix)//"_BLZTRP", "_fermi.dx',           'unknown',    'formatted',0"
   write (iout, '(3a)') "32,'", trim(fname_radix)//"_BLZTRP", "_sigxx.dx',           'unknown',    'formatted',0"
   write (iout, '(3a)') "33,'", trim(fname_radix)//"_BLZTRP", "_sigyy.dx',           'unknown',    'formatted',0"
   write (iout, '(3a)') "34,'", trim(fname_radix)//"_BLZTRP", "_sigzz.dx',           'unknown',    'formatted',0"
   write (iout, '(3a)') "35,'", trim(fname_radix)//"_BLZTRP", "_band.dat',           'unknown',    'formatted',0"
   write (iout, '(3a)') "36,'", trim(fname_radix)//"_BLZTRP", "_band.gpl',           'unknown',    'formatted',0"
   write (iout, '(3a)') "37,'", trim(fname_radix)//"_BLZTRP", "_deriv.dat',           'unknown',    'formatted',0"
   write (iout, '(3a)') "38,'", trim(fname_radix)//"_BLZTRP", "_mass.dat',           'unknown',    'formatted',0"

   close(iout)
 end do !isppol

!file is for geometry symmetries etc
 filename= trim(fname_radix)//"_BLZTRP.struct"
 if (open_file(filename, msg, newunit=iout, form='formatted') /= 0) then
   MSG_ERROR(msg)
 end if

 write (iout, '(a)') "BoltzTraP geometry file generated by ABINIT."

!here we need to print out the unit cell vectors
 write (iout, '(3E20.10)') crystal%rprimd(:,1)
 write (iout, '(3E20.10)') crystal%rprimd(:,2)
 write (iout, '(3E20.10)') crystal%rprimd(:,3)
 write (iout, '(I7)') crystal%nsym

 do isym=1,crystal%nsym
   write (iout,'(3(3I5,2x), a, I5)') &
&   crystal%symrel(1,:,isym), &
&   crystal%symrel(2,:,isym), &
&   crystal%symrel(3,:,isym), &
&   ' ! symmetry rotation matrix isym = ', isym
 end do

 close (iout)

! second file is for eigenvalues
! two file names for each spin, if necessary
 do isppol=1,nsppol
   filename=trim(fname_radix)//"_BLZTRP."//spinsuffix(isppol)//"energy"//trim(so_suffix)

   if (open_file (filename, msg, newunit=iout, form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   write (iout, '(a,I5)') "BoltzTraP eigen-energies file generated by ABINIT. ispin = ", isppol
   write (iout, '(I7, I7, E20.10, a)') &
&    ebands%nkpt, nsppol, ha2ryd*ebands%fermie, '     ! nk, nspin, Fermi level(Ry) : energies below in Ry'

   do ikpt=1,ebands%nkpt
!    these need to be in reduced coordinates
     write (iout, '(3E20.10, I7, a)') &
&      ebands%kptns(1,ikpt), ebands%kptns(2,ikpt), ebands%kptns(3,ikpt), nband, '    ! kpt nband'
     do iband=1,nband
!      output in Rydberg
       write (iout, '(E20.10)') ha2ryd*ebands%eig(iband, ikpt, isppol)
     end do
   end do

   close (iout)
 end do

!this file is for tau_k
 if (present (tau_k)) then
   do isppol = 1, nsppol
     filename= trim(fname_radix)//"_"//spinsuffix(isppol)//"BLZTRP.tau_k"
     if (open_file(filename, msg, newunit=iout, form='formatted') /= 0) then
       MSG_ERROR(msg)
     end if

     write (iout, '(a)') "BoltzTraP tau_k file generated by ANADDB."
     write (iout, '(I7, I7, E20.10, a)')&
&      ebands%nkpt, nsppol, ha2ryd*ebands%fermie, '     ! nk, nspin, Fermi level(Ry) : energies below in Ry'

     do ikpt=1,ebands%nkpt
!      these need to be in reduced coordinates
       write (iout, '(3E20.10, I7, a)') &
&        ebands%kptns(1,ikpt), ebands%kptns(2,ikpt), ebands%kptns(3,ikpt), nband, '    ! kpt nband'
       do iband=1,nband
!        output in eV
         write (iout, '(E20.10)') tau_k(isppol,ikpt,iband)
       end do
     end do
     close (iout)
   end do

 end if

end subroutine ebands_prtbltztrp
!!***

!!****f* m_ebands/ebands_prtbltztrp_tau_out
!! NAME
!! ebands_prtbltztrp_tau_out
!!
!! FUNCTION
!!   output files for BoltzTraP code, which integrates Boltzmann transport quantities
!!   over the Fermi surface for different T and chemical potentials. Abinit provides
!!   all necessary input files: struct, energy, input file, and def file for the unit
!!   definitions of fortran files in BT.
!!   See http://www.icams.de/content/departments/ams/madsen/boltztrap.html
!!   Output T-depedent tau_k, modified from ebands_prtbltztrp
!!
!! INPUTS
!!  eigen(mband*nkpt*nsppol) = array for holding eigenvalues (hartree)
!!  fermie = Fermi level
!!  fname_radix = radix of file names for output
!!  nband = number of bands
!!  nkpt = number of k points.
!!  nsppol = 1 for unpolarized, 2 for spin-polarized
!!  nsym = number of symmetries in space group
!!  rprimd(3,3) = dimensional primitive translations for real space (bohr)
!!  symrel = symmetry operations in reduced coordinates, real space
!!  to be used in future  xred(3,natom) = reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! PARENTS
!!      get_tau_k
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_prtbltztrp_tau_out (eigen, tempermin, temperinc, ntemper, fermie, fname_radix, kpt, &
       nband, nelec, nkpt, nspinor, nsppol, nsym, rprimd, symrel, tau_k)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nsym, nband, nkpt, nsppol, nspinor, ntemper
 real(dp), intent(in) :: tempermin, temperinc
 real(dp), intent(in) :: nelec
 character(len=fnlen), intent(in) :: fname_radix
!arrays
 real(dp), intent(in) :: fermie(ntemper)
 integer, intent(in) :: symrel(3,3,nsym)
 real(dp), intent(in) :: kpt(3,nkpt)
 real(dp), intent(in) :: eigen(nband, nkpt, nsppol)
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), intent(in) :: tau_k(ntemper,nsppol,nkpt,nband)

!Local variables-------------------------------
!scalars
 integer :: iout, isym, iband, isppol, ikpt, itemp
 real(dp) :: Temp
 real(dp),parameter :: ha2ryd = two
 character(len=500) :: msg
 character(len=fnlen) :: filename,appendix

! *************************************************************************

!input file for boltztrap: general info, Ef, Nelec, etc...
 do itemp = 1, ntemper
   write(appendix,"(i0)") itemp
   filename= trim(fname_radix)//"_BLZTRP.intrans_"//trim(appendix)
   if (open_file(filename, msg, newunit=iout, form="formatted", action="write") /= 0) then
     MSG_ERROR(msg)
   end if

   write (iout, '(a)') "GENE                      # Format of input: generic format, with Symmetries"
   write (iout, '(a)') "0 0 0 0.0                 # iskip (not presently used) idebug setgap shiftgap"
   write (iout, '(E15.5,a,F10.4,a)') fermie(itemp)*two, " 0.0005 0.4  ", nelec, &
    "  # Fermilevel (Ry), energy grid spacing, energy span around Fermilevel, number of electrons"
   write (iout, '(a)') "CALC                      # CALC (calculate expansion coeff), NOCALC read from file"
   write (iout, '(a)') "3                         # lpfac, number of latt-points per k-point"
   write (iout, '(a)') "BOLTZ                     # run mode (only BOLTZ is supported)"
   write (iout, '(a)') ".15                       # (efcut) energy range of chemical potential"
   write (iout, '(2f8.2,a)')&
    tempermin+temperinc*dble(itemp),tempermin+temperinc*dble(itemp), "                  # Tmax, temperature grid spacing"
   write (iout, '(2a)') "-1                        # energyrange of bands given ",&
    "individual DOS output sig_xxx and dos_xxx (xxx is band number)"
   write (iout, '(a)') "TETRA                     # DOS calculation method. Other possibility is TETRA"
   write (iout, '(a)') "No                        # not using model for relaxation time"
   write (iout, '(a)') "3                         # Number of doping levels coefficients will be output for"
   write (iout, '(a)') "-1.e16 0.0d0 1.e16        # Values of doping levels (in carriers / cm^3"
   close(iout)
 end do

!files file, with association of all units for Boltztrap
 filename= trim(fname_radix)//"_BLZTRP.def"
 if (open_file(filename, msg, newunit=iout, form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write (iout, '(3a)') "5, '", trim(fname_radix)//"_BLZTRP", ".intrans',      'old',    'formatted',0"
 write (iout, '(3a)') "6, '", trim(fname_radix)//"_BLZTRP", ".outputtrans',      'unknown',    'formatted',0"
 write (iout, '(3a)') "20,'", trim(fname_radix)//"_BLZTRP", ".struct',         'old',    'formatted',0"
 if (nspinor == 1) then
   write (iout, '(3a)') "10,'", trim(fname_radix)//"_BLZTRP", ".energy',         'old',    'formatted',0"
 else if (nspinor == 2) then
   write (iout, '(3a)') "10,'", trim(fname_radix)//"_BLZTRP", ".energyso',         'old',    'formatted',0"
 end if
 write (iout, '(3a)') "10,'", trim(fname_radix)//"_BLZTRP", ".energy',         'old',    'formatted',0"
 write (iout, '(3a)') "11,'", trim(fname_radix)//"_BLZTRP", ".tau_k',         'old',    'formatted',0"
 write (iout, '(3a)') "48,'", trim(fname_radix)//"_BLZTRP", ".engre',         'unknown',    'unformatted',0"
 write (iout, '(3a)') "49,'", trim(fname_radix)//"_BLZTRP", ".transdos',        'unknown',    'formatted',0"
 write (iout, '(3a)') "50,'", trim(fname_radix)//"_BLZTRP", ".sigxx',        'unknown',    'formatted',0"
 write (iout, '(3a)') "51,'", trim(fname_radix)//"_BLZTRP", ".sigxxx',        'unknown',    'formatted',0"
 write (iout, '(3a)') "21,'", trim(fname_radix)//"_BLZTRP", ".trace',           'unknown',    'formatted',0"
 write (iout, '(3a)') "22,'", trim(fname_radix)//"_BLZTRP", ".condtens',           'unknown',    'formatted',0"
 write (iout, '(3a)') "24,'", trim(fname_radix)//"_BLZTRP", ".halltens',           'unknown',    'formatted',0"
 write (iout, '(3a)') "25,'", trim(fname_radix)//"_BLZTRP", ".trace_fixdoping',     'unknown',    'formatted',0"
 write (iout, '(3a)') "26,'", trim(fname_radix)//"_BLZTRP", ".condtens_fixdoping',           'unknown',    'formatted',0"
 write (iout, '(3a)') "27,'", trim(fname_radix)//"_BLZTRP", ".halltens_fixdoping',           'unknown',    'formatted',0"
 write (iout, '(3a)') "30,'", trim(fname_radix)//"_BLZTRP", "_BZ.dx',           'unknown',    'formatted',0"
 write (iout, '(3a)') "31,'", trim(fname_radix)//"_BLZTRP", "_fermi.dx',           'unknown',    'formatted',0"
 write (iout, '(3a)') "32,'", trim(fname_radix)//"_BLZTRP", "_sigxx.dx',           'unknown',    'formatted',0"
 write (iout, '(3a)') "33,'", trim(fname_radix)//"_BLZTRP", "_sigyy.dx',           'unknown',    'formatted',0"
 write (iout, '(3a)') "34,'", trim(fname_radix)//"_BLZTRP", "_sigzz.dx',           'unknown',    'formatted',0"
 write (iout, '(3a)') "35,'", trim(fname_radix)//"_BLZTRP", "_band.dat',           'unknown',    'formatted',0"
 write (iout, '(3a)') "36,'", trim(fname_radix)//"_BLZTRP", "_band.gpl',           'unknown',    'formatted',0"
 write (iout, '(3a)') "37,'", trim(fname_radix)//"_BLZTRP", "_deriv.dat',           'unknown',    'formatted',0"
 write (iout, '(3a)') "38,'", trim(fname_radix)//"_BLZTRP", "_mass.dat',           'unknown',    'formatted',0"
 close(iout)

!file is for geometry symmetries etc
 filename= trim(fname_radix)//"_BLZTRP.struct"
 if (open_file(filename, msg, newunit=iout, form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write (iout, '(a)') "BoltzTraP geometry file generated by ABINIT."

!here we need to print out the unit cell vectors
 write (iout, '(3E20.10)') rprimd(:,1)
 write (iout, '(3E20.10)') rprimd(:,2)
 write (iout, '(3E20.10)') rprimd(:,3)
 write (iout, '(I7)') nsym

 do isym=1, nsym
   write (iout,'(3(3I5,2x), a, I5)') &
    symrel(1,:,isym), symrel(2,:,isym), symrel(3,:,isym), ' ! symmetry rotation matrix isym = ', isym
 end do
 close (iout)

!second file is for eigenvalues
 if (nspinor == 1) then
   filename= trim(fname_radix)//"_BLZTRP.energy"
 else if (nspinor == 2) then
   filename= trim(fname_radix)//"_BLZTRP.energyso"
 end if

 if (open_file(filename, msg, newunit=iout, form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 write (iout, '(a)') "BoltzTraP eigen-energies file generated by ABINIT."
 write (iout, '(I7, I7, E20.10, a)') nkpt, nsppol, ha2ryd*fermie(1), '     ! nk, nspin, Fermi level(Ry) : energies below in Ry'
 do isppol = 1, nsppol
   do ikpt = 1, nkpt
!    these need to be in reduced coordinates
     write (iout, '(3E20.10, I7, a)') kpt(1,ikpt), kpt(2,ikpt), kpt(3,ikpt), nband, '    ! kpt nband'
     do iband = 1, nband
!      output in eV
       write (iout, '(E20.10)') ha2ryd*eigen(iband, ikpt, isppol)
     end do
   end do
 end do
 close (iout)

!this file is for tau_k
 do itemp = 1, ntemper
   Temp=tempermin+temperinc*dble(itemp)

   write(appendix,"(i0)") itemp
   filename= trim(fname_radix)//"_BLZTRP.tau_k_"//trim(appendix)
   if (open_file(filename, msg, newunit=iout, form="formatted", action="write") /= 0) then
     MSG_ERROR(msg)
   end if
   write (iout, '(a,f12.6)') "BoltzTraP tau_k file generated by ANADDB for T= ", Temp
   write (iout, '(I7, I7, E20.10, a)') nkpt, nsppol, ha2ryd*fermie(itemp), &
   '     ! nk, nspin, Fermi level(Ry) : energies below in Ry'
   do isppol = 1, nsppol
     do ikpt = 1, nkpt
!      these need to be in reduced coordinates
       write (iout, '(3E20.10, I7, a)') kpt(1,ikpt), kpt(2,ikpt), kpt(3,ikpt), nband, '    ! kpt nband'
       do iband = 1, nband
!        output in sec
         write (iout, '(E20.10)') tau_k(itemp,isppol,ikpt,iband)
       end do
     end do
   end do
   close (iout)
 end do

end subroutine ebands_prtbltztrp_tau_out
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_write
!! NAME
!! ebands_write
!!
!! FUNCTION
!!  Driver routine to write bands in different (txt) formats.
!!  This routine should be called by a single processor.
!!
!! INPUTS
!!  prtebands=Flag seleecting the output format:
!!    0 --> None
!!    1 --> xmgrace
!!    2 --> gnuplot     (not coded yet)
!!    3 --> EIG format  (not coded yet)
!!  prefix=Prefix for output filename.
!!  [kptbounds(:,:)]=Optional argument giving the extrema of the k-path.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      eph,m_ebands,outscfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_write(ebands, prtebands, prefix, kptbounds)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: prtebands
 type(ebands_t),intent(in) :: ebands
 character(len=*),intent(in) :: prefix
 real(dp),optional,intent(in) :: kptbounds(:,:)

! *********************************************************************

 select case (prtebands)
 case (0)
    return
 case (1)
   !call wrtout(std_out, sjoin(" Writing interpolated bands to:",  path)
   if (present(kptbounds)) then
     call ebands_write_xmgrace(ebands, strcat(prefix, "_EBANDS.agr"), kptbounds=kptbounds)
   else
     call ebands_write_xmgrace(ebands, strcat(prefix, "_EBANDS.agr"))
   end if
 case (2)
   !call wrtout(std_out, sjoin(" Writing interpolated bands to:",  path)
   if (present(kptbounds)) then
     call ebands_write_gnuplot(ebands, prefix, kptbounds=kptbounds)
   else
     call ebands_write_gnuplot(ebands, prefix)
   end if
 case default
   MSG_WARNING(sjoin("Unsupported value for prtebands:", itoa(prtebands)))
 end select

end subroutine ebands_write
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_write_xmgrace
!! NAME
!! ebands_write_xmgrace
!!
!! FUNCTION
!!  Write bands in Xmgrace format. This routine should be called by a single processor.
!!  Use the driver `ebands_write` to support different formats.
!!
!! INPUTS
!!  filename=Filename
!!  [kptbounds(:,:)]=Optional argument giving the extrema of the k-path.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_ebands
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_write_xmgrace(ebands, filename, kptbounds)

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in) :: ebands
 character(len=*),intent(in) :: filename
 real(dp),optional,intent(in) :: kptbounds(:,:)

!Local variables-------------------------------
!scalars
 integer :: unt,ik,spin,band,ii,start,nkbounds
 character(len=500) :: msg
!arrays
 integer :: g0(3)
 integer,allocatable :: bounds2kpt(:)

! *********************************************************************

 nkbounds = 0
 if (present(kptbounds)) then
   if (product(shape(kptbounds)) > 0 ) then
     ! Find correspondence between kptbounds and k-points in ebands.
     nkbounds = size(kptbounds, dim=2)
     ABI_MALLOC(bounds2kpt, (nkbounds))
     bounds2kpt = 1; start = 1
     do ii=1,nkbounds
        do ik=start,ebands%nkpt
          if (isamek(ebands%kptns(:, ik), kptbounds(:, ii), g0)) then
            bounds2kpt(ii) = ik; start = ik + 1; exit
          end if
        end do
     end do
   end if
 end if

 if (open_file(filename, msg, newunit=unt, form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write(unt,'(a)') "# Grace project file"
 write(unt,'(a)') "# Generated by Abinit"
 write(unt,'(4(a,i0))') &
   "# mband: ",ebands%mband,", nkpt: ",ebands%nkpt,", nsppol: ",ebands%nsppol,", nspinor: ",ebands%nspinor
 write(unt,'(a,f8.2,a,i0,2(a,f8.2))') &
   "# nelect: ",ebands%nelect,", occopt: ",ebands%occopt,", tsmear: ",ebands%tsmear,", tphysel: ",ebands%tphysel
 write(unt,'(a,f8.2,a)') "# Energies are in eV. Zero set to efermi, previously it was at: ",ebands%fermie * Ha_eV, " [eV]"
 write(unt,'(a)')"# List of k-points and their index (C notation i.e. count from 0)"
 do ik=1,ebands%nkpt
   write(unt, "(a)")sjoin("#", itoa(ik-1), ktoa(ebands%kptns(:,ik)))
 end do
 write(unt,'(a)') "@page size 792, 612"
 write(unt,'(a)') "@page scroll 5%"
 write(unt,'(a)') "@page inout 5%"
 write(unt,'(a)') "@link page off"
 write(unt,'(a)') "@with g0"
 write(unt,'(a)') "@world xmin 0.00"
 write(unt,'(a,i0)') '@world xmax ',ebands%nkpt
 write(unt,'(a,es16.8)') '@world ymin ',minval((ebands%eig - ebands%fermie) * Ha_eV)
 write(unt,'(a,es16.8)') '@world ymax ',maxval((ebands%eig - ebands%fermie) * Ha_eV)
 write(unt,'(a)') '@default linewidth 1.5'
 write(unt,'(a)') '@xaxis  tick on'
 write(unt,'(a)') '@xaxis  tick major 1'
 write(unt,'(a)') '@xaxis  tick major color 1'
 write(unt,'(a)') '@xaxis  tick major linestyle 3'
 write(unt,'(a)') '@xaxis  tick major grid on'
 write(unt,'(a)') '@xaxis  tick spec type both'
 write(unt,'(a)') '@xaxis  tick major 0, 0'
 if (nkbounds /= 0) then
   write(unt,'(a,i0)') '@xaxis  tick spec ',nkbounds
   do ik=1,nkbounds
     !write(unt,'(a,i0,a,a)') '@xaxis  ticklabel ',ik-1,',', "foo"
     write(unt,'(a,i0,a,i0)') '@xaxis  tick major ',ik-1,' , ',bounds2kpt(ik) - 1
   end do
 end if
 write(unt,'(a)') '@xaxis  ticklabel char size 1.500000'
 write(unt,'(a)') '@yaxis  tick major 10'
 write(unt,'(a)') '@yaxis  label "Band Energy [eV]"'
 write(unt,'(a)') '@yaxis  label char size 1.500000'
 write(unt,'(a)') '@yaxis  ticklabel char size 1.500000'
 ii = -1
 do spin=1,ebands%nsppol
   do band=1,ebands%mband
     ii = ii + 1
     write(unt,'(a,i0,a,i0)') '@    s',ii,' line color ',spin
   end do
 end do
 ii = -1
 do spin=1,ebands%nsppol
   do band=1,ebands%mband
     ii = ii + 1
     write(unt,'(a,i0)') '@target G0.S',ii
     write(unt,'(a)') '@type xy'
     do ik=1,ebands%nkpt
        write(unt,'(i0,1x,es16.8)') ik-1, (ebands%eig(band, ik, spin) - ebands%fermie) * Ha_eV
     end do
     write(unt,'(a)') '&'
   end do
 end do

 close(unt)

 ABI_SFREE(bounds2kpt)

end subroutine ebands_write_xmgrace
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_write_gnuplot
!! NAME
!! ebands_write_gnuplot
!!
!! FUNCTION
!!  Write bands in gnuplot format. This routine should be called by a single processor.
!!  Use the driver `ebands_write` to support different formats.
!!
!! INPUTS
!!  prefix=prefix for files (.data, .gnuplot)
!!  [kptbounds(:,:)]=Optional argument giving the extrema of the k-path.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_ebands
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_write_gnuplot(ebands, prefix, kptbounds)

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in) :: ebands
 character(len=*),intent(in) :: prefix
 real(dp),optional,intent(in) :: kptbounds(:,:)

!Local variables-------------------------------
!scalars
 integer :: unt,gpl_unt,ik,spin,ii,start,nkbounds
 character(len=500) :: msg,fmt
 character(len=fnlen) :: datafile,basefile
!arrays
 integer :: g0(3)
 integer,allocatable :: bounds2kpt(:)

! *********************************************************************

 nkbounds = 0
 if (present(kptbounds)) then
   if (product(shape(kptbounds)) > 0 ) then
     ! Find correspondence between kptbounds and k-points in ebands.
     nkbounds = size(kptbounds, dim=2)
     ABI_MALLOC(bounds2kpt, (nkbounds))
     bounds2kpt = 1; start = 1
     do ii=1,nkbounds
        do ik=start,ebands%nkpt
          if (isamek(ebands%kptns(:, ik), kptbounds(:, ii), g0)) then
            bounds2kpt(ii) = ik; start = ik + 1; exit
          end if
        end do
     end do
   end if
 end if

 datafile = strcat(prefix, "_EBANDS.data")
 if (open_file(datafile, msg, newunit=unt, form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 if (open_file(strcat(prefix, "_EBANDS.gnuplot"), msg, newunit=gpl_unt, form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if
 basefile = basename(datafile)

 write(unt,'(a)') "# Electron band structure data file"
 write(unt,'(a)') "# Generated by Abinit"
 write(unt,'(4(a,i0))') &
   "# mband: ",ebands%mband,", nkpt: ",ebands%nkpt,", nsppol: ",ebands%nsppol,", nspinor: ",ebands%nspinor
 write(unt,'(a,f8.2,a,i0,2(a,f8.2))') &
   "# nelect: ",ebands%nelect,", occopt: ",ebands%occopt,", tsmear: ",ebands%tsmear,", tphysel: ",ebands%tphysel
 write(unt,'(a,f8.2,a)') "# Energies are in eV. Zero set to efermi, Previously it was at: ",ebands%fermie * Ha_eV, " [eV]"
 write(unt,'(a)')"# List of k-points and their index (C notation i.e. count from 0)"
 do ik=1,ebands%nkpt
   write(unt, "(a)")sjoin("#", itoa(ik-1), ktoa(ebands%kptns(:,ik)))
 end do

 fmt = sjoin("(i0,1x,", itoa(ebands%mband), "(es16.8,1x))")
 write(unt,'(a)') ' '
 do spin=1,ebands%nsppol
   write(unt,'(a,i0)') '# [kpt-index, band_1, band_2 ...]  for spin: ',spin
   do ik=1,ebands%nkpt
     write(unt,fmt) ik-1, (ebands%eig(:, ik, spin) - ebands%fermie) * Ha_eV
   end do
   write(unt,'(a)') ' '
 end do

  ! gnuplot script file
  write(gpl_unt,'(a)') '# File to plot phonon bandstructure with gnuplot'
  write(gpl_unt,'(a)') "#set terminal postscript eps enhanced color font 'Times-Roman,26' lw 2"
  write(gpl_unt,'(a)') '#use the next lines to make a nice figure for a paper'
  write(gpl_unt,'(a)') '#set term postscript enhanced eps color lw 0.5 dl 0.5'
  write(gpl_unt,'(a)') '#set pointsize 0.275'
  write(gpl_unt,'(a)') 'set palette defined ( 0 "blue", 3 "green", 6 "yellow", 10 "red" )'
  write(gpl_unt,'(a)') 'unset key'
  write(gpl_unt,'(a)') '# can make pointsize smaller (~0.5). Too small and nothing is printed'
  write(gpl_unt,'(a)') 'set pointsize 0.8'
  write(gpl_unt,'(a)') 'set view 0,0'
  write(gpl_unt,'(a,i0,a)') 'set xrange [0:',ebands%nkpt-1,']'
  write(gpl_unt,'(2(a,es16.8),a)')&
    'set yrange [',minval((ebands%eig - ebands%fermie) * Ha_eV),':',maxval((ebands%eig - ebands%fermie) * Ha_eV),']'
  write(gpl_unt,'(a)') 'set xlabel "Momentum"'
  write(gpl_unt,'(a)') 'set ylabel "Energy [eV]"'
  write(gpl_unt,'(a)') strcat('set title "', replace(basefile, "_", "\\_"), '"')
  if (nkbounds == 0) then
    write(gpl_unt,'(a)') 'set grid xtics'
  else
    write(gpl_unt,"(a)")"# Add vertical lines in correspondence of high-symmetry points."
    write(gpl_unt,'(a)') 'unset xtics'
    do ii=1,nkbounds
      write(gpl_unt,"(a,2(i0,a))") &
        "set arrow from ",bounds2kpt(ii)-1,",graph(0,0) to ",bounds2kpt(ii)-1,",graph(1,1) nohead ls 'dashed'"
      !write(gpl_unt,"(a)")sjoin("set xtics add('kname'", itoa(bounds2kpt(ii)-1), ")")
    end do

  end if
  write(gpl_unt,"(a)")sjoin("mband =", itoa(ebands%mband))
  write(gpl_unt,"(a)")strcat('plot for [i=2:mband] "', basefile, '" u 1:i every :1 with lines linetype -1')
  if (ebands%nsppol == 2) then
    write(gpl_unt,"(a)")strcat('replot for [i=2:mband] "', basefile, '" u 1:i every :2 with lines linetype 4')
  end if
 write(gpl_unt, "(a)")"pause -1"

 close(unt)
 close(gpl_unt)

 ABI_SFREE(bounds2kpt)

end subroutine ebands_write_gnuplot
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_interpolate_kpath
!!
!! NAME
!!  ebands_interpolate_kpath
!!
!! FUNCTION
!!
!! INPUTS
!!  dtset<dataset_type>=Abinit dataset
!!  band_block(2)=Initial and final band index to be interpolated. [0,0] if all bands are used.
!!    This is a global variable i.e. all MPI procs must call the routine with the same value.
!!
!! OUTPUT
!!
!! PARENTS
!!      outscfcv,sigma
!!
!! SOURCE

subroutine ebands_interpolate_kpath(ebands, dtset, cryst, band_block, prefix, comm)

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in) :: ebands
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: comm
 character(len=*),intent(in) :: prefix
!arrays
 integer,intent(in) :: band_block(2)

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0, intp_nshiftk1 = 1
 integer :: my_rank, ndivsm, nbounds, itype !, spin, ik, ib, ii, jj, ierr
 type(ebands_t) :: ebands_kpath
 type(kpath_t) :: kpath
 character(len=500) :: tag !msg
!arrays
 real(dp),allocatable :: bounds(:,:)

! *********************************************************************

 my_rank = xmpi_comm_rank(comm)

 itype = nint(dtset%einterp(1)); tag =  "_SKW"
 tag = "_INTERP"

 ! Generate k-path
 ndivsm = dtset%ndivsm
 if (ndivsm <= 0) then
   MSG_COMMENT("Setting ndivsm to 20 because variable is not given in input file")
   ndivsm = 20
 end if
 nbounds = dtset%nkpath
 if (nbounds <= 0) then
   MSG_COMMENT("Using hard-coded k-path because nkpath not present in input file.")
   nbounds = 5
   ABI_MALLOC(bounds, (3, 5))
   bounds = reshape([zero, zero, zero, half, zero, zero, zero, half, zero, zero, zero, zero, zero, zero, half], [3,5])
 else
   call alloc_copy(dtset%kptbounds, bounds)
 end if

 kpath = kpath_new(bounds, cryst%gprimd, ndivsm)
 call kpath%print(header="Interpolating energies on k-path", unit=std_out)
 ABI_FREE(bounds)

 ! Interpolate bands on k-path.
 ebands_kpath = ebands_interp_kpath(ebands, cryst, kpath, dtset%einterp, band_block, comm)

 if (my_rank == master) then
   call wrtout(ab_out, sjoin("- Writing interpolated bands to file:", strcat(prefix, tag)))
   call ebands_write(ebands_kpath, dtset%prtebands, strcat(prefix, tag), kptbounds=kpath%bounds)
 end if

 call ebands_free(ebands_kpath)
 call kpath%free()

end subroutine ebands_interpolate_kpath
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/klinterp_new
!! NAME
!! klinterp_new
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(klinterp_t) function klinterp_new(cryst, kptrlatt, nshiftk, shiftk, kptopt, kibz, &
                                       bsize, nkibz, nsppol, ndat, values_bksd, comm) result(new)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: nshiftk, kptopt, bsize, nkibz, nsppol, ndat, comm
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: kibz(3, nkibz), shiftk(3,nshiftk), values_bksd(bsize, nkibz, nsppol, ndat)

!Local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1 = 1
 integer :: ierr, nkfull, ikf, spin, band, ik_ibz, timrev, ix, iy, iz, nkx, nky, nkz, idat
 real(dp) :: dksqmax
 character(len=500) :: msg
!arrays
 integer,allocatable :: bz2ibz(:,:)
 real(dp) :: kpt(3)
 real(dp),allocatable :: kfull(:,:)

! *********************************************************************

 ! Check input parameters
 ierr = 0
 if (nkibz == 1) then
   MSG_ERROR_NOSTOP("Cannot interpolate with a single k-point", ierr)
 end if
 if (.not. isdiagmat(kptrlatt)) then
   MSG_ERROR_NOSTOP('kptrlatt is not diagonal. Multiple shifts are not allowed', ierr)
 end if
 if (nshiftk /= 1) then
   MSG_ERROR_NOSTOP('Multiple shifts not allowed', ierr)
 end if
 if (any(abs(shiftk(:, 1)) > tol8)) then
   MSG_ERROR_NOSTOP("shifted k-mesh not implented", ierr)
 end if

 if (ierr /= 0) then
   MSG_ERROR("Linear interpolation cannot be performed. See messages above.")
 end if

 nkx = kptrlatt(1, 1)
 nky = kptrlatt(2, 2)
 nkz = kptrlatt(3, 3)

 new%nkx = nkx; new%nky = nky; new%nkz = nkz
 new%bsize = bsize; new%nsppol = nsppol; new%ndat = ndat

 ! Build list of k-points in the conventional unit cell.
 ! (x,y,z) ordered as required by interpolation routine
 nkfull = nkx * nky * nkz
 ABI_MALLOC(kfull, (3, nkfull))
 ikf = 0
 do iz=1,nkz
   kpt(3) = (iz - 1 + shiftk(3, 1)) / nkz
   do iy=1,nky
     kpt(2) = (iy - 1 + shiftk(2, 1)) / nky
     do ix=1,nkx
       kpt(1) = (ix - 1 + shiftk(1, 1)) / nkx
       ikf = ikf + 1
       kfull(:, ikf) = kpt
     end do
   end do
 end do

 ! Build mapping kfull --> IBZ
 timrev = kpts_timrev_from_kptopt(kptopt)
 ABI_MALLOC(bz2ibz, (nkfull*sppoldbl1, 6))

 call listkk(dksqmax, cryst%gmet, bz2ibz, kibz, kfull, nkibz, nkfull, cryst%nsym,&
   sppoldbl1, cryst%symafm, cryst%symrec, timrev, comm, exit_loop=.True., use_symrec=.True.)

 ABI_FREE(kfull)

 if (dksqmax > tol12) then
   write(msg, '(3a,es16.6,4a)' )&
   'At least one of the k points could not be generated from a symmetrical one.',ch10,&
   'dksqmax: ',dksqmax,ch10,&
   'Action: check k-point input variables',ch10,&
   '        e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
   MSG_ERROR(msg)
 end if

 ABI_CALLOC(new%data_uk_bsd, (nkx*nky*nkz, bsize, nsppol, ndat))

 ! Build array in the full BZ to prepare call to interpol3d.
 ikf = 0
 do iz=1,nkz
   do iy=1,nky
     do ix=1,nkx
       ikf = ikf + 1
       ik_ibz = bz2ibz(ikf, 1)
       new%data_uk_bsd(ikf, 1:bsize, 1:nsppol, 1:ndat) = values_bksd(1:bsize, ik_ibz, 1:nsppol, 1:ndat)
     end do
   end do
 end do

 ABI_FREE(bz2ibz)

end function klinterp_new
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/klinterp_free
!! NAME
!! klinterp_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine klinterp_free(self)

!Arguments ------------------------------------
!scalars
 class(klinterp_t),intent(inout) :: self

! *********************************************************************

 ABI_SFREE(self%data_uk_bsd)

end subroutine klinterp_free
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/klinterp_eval_bsd
!! NAME
!! klinterp_eval_bsd
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine klinterp_eval_bsd(self, kpt, vals_bsd)

!Arguments ------------------------------------
!scalars
 class(klinterp_t),intent(in) :: self
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(out) :: vals_bsd(self%bsize, self%nsppol, self%ndat)

!Local variables-------------------------------
 integer :: spin, idat, band
 real(dp) :: kwrap(3), shift(3)

! *********************************************************************

 call wrap2_zero_one(kpt, kwrap, shift)
 !write(std_out, *)"kwrap:", kwrap
 do idat=1,self%ndat
   do spin=1,self%nsppol
      do band=1,self%bsize
        vals_bsd(band, spin, idat) = interpol3d(kwrap, self%nkx, self%nky, self%nkz, self%data_uk_bsd(:, band, spin, idat))
      end do
   end do
 end do

end subroutine klinterp_eval_bsd
!!***

end module m_ebands
!!***
