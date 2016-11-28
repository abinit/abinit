!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ebands
!! NAME
!!  m_ebands
!!
!! FUNCTION
!!  This module contains utilities to analyze and retrieve information from the ebands_t.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2016 ABINIT group (MG, MJV, BXu)
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
 use m_profiling_abi
 use m_xmpi
 use m_tetrahedron
 use m_bspline
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr

 use defs_datatypes,   only : ebands_t
 use defs_abitypes,    only : hdr_type, dataset_type
 use m_copy,           only : alloc_copy
 use m_io_tools,       only : file_exists, open_file
 use m_fstrings,       only : tolower, itoa, sjoin, ftoa, ltoa, ktoa
 use m_numeric_tools,  only : arth, imin_loc, imax_loc, bisect, stats_t, stats_eval, simpson_int, wrap2_zero_one,&
                              isdiagmat
 use m_special_funcs,  only : dirac_delta
 use m_geometry,       only : normv
 use m_cgtools,        only : set_istwfk
 use m_nesting,        only : mknesting
 use m_crystal,        only : crystal_t
 use m_bz_mesh,        only : kmesh_t, isamek
 use m_fftcore,        only : get_kg

 implicit none

 private

 public :: ebands_init             ! Main creation method.
 public :: ebands_from_hdr         ! Init object from the abinit header.
 public :: ebands_from_dtset       ! Init object from the abinit dataset.
 public :: ebands_free             ! Destruction method.
 public :: ebands_copy             ! Deep copy of the ebands_t.
 public :: ebands_print            ! Printout basic info on the data type.
 public :: unpack_eneocc           ! Helper function for reshaping (energies|occupancies|derivate of occupancies).
 public :: pack_eneocc             ! Helper function for reshaping (energies|occupancies|derivate of occupancies).
 public :: get_eneocc_vect         ! Reshape (ene|occ|docdde) returning a matrix instead of a vector.
 public :: put_eneocc_vect         ! Put (ene|occ|doccde) in vectorial form into the data type doing a reshape.
 public :: get_bandenergy          ! Returns the band energy of the system.
 public :: get_valence_idx         ! Gives the index of the (valence|bands at E_f).
 public :: apply_scissor           ! Apply a scissor operator (no k-dependency)
 public :: get_occupied            ! Returns band indeces after wich occupations are less than an input value.
 public :: enclose_degbands        ! Adjust band indeces such that all degenerate states are treated.
 public :: ebands_nelect_per_spin  ! Returns number of electrons per spin channel
 public :: get_minmax              ! Returns min and Max value of (eig|occ|doccde).
 public :: ebands_edstats          ! Compute statistical parameters of the energy differences e_ks[b+1] - e_ks[b]
 public :: ebands_has_metal_scheme ! .True. if metallic occupation scheme is used.
 public :: ebands_write_bxsf       ! Write 3D energies for Fermi surface visualization (XSF format)
 public :: ebands_update_occ       ! Update the occupation numbers.
 public :: ebands_set_scheme       ! set the occupation scheme.
 public :: ebands_set_fermie       ! Change the fermi level (assume metallic scheme).
 public :: ebands_set_nelect       ! Change the number of electrons (assume metallic scheme).
 public :: ebands_report_gap       ! Print info on the fundamental and optical gap.
 public :: ebands_ncwrite          ! Dump the object into NETCDF file (use ncid)
 public :: ebands_ncwrite_path     ! Dump the object into NETCDF file (use filepath)
 public :: ebands_write_nesting    ! Calculate the nesting function and output data to file.
 public :: ebands_expandk          ! Build a new ebands_t in the full BZ.
 public :: ebands_jdos             ! Compute the joint density of states.
 public :: ebands_bspline

 public :: ebands_prtbltztrp
 public :: ebands_prtbltztrp_tau_out
 public :: ebands_write_xmgrace
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

   integer :: ief=0
   ! Rightmost Index of the energy mesh such as IDOS[mesh[ief]] < nelect.
   ! 0 if Fermi level could not be computed
   ! Note the value of gef stored in edos_t is computed by performing
   ! a linear interpolation between ief and ief+1

   integer :: intmeth
   ! 1 for gaussian, 2 tetra

   real(dp) :: broad=zero
   ! Gaussian broadening

   real(dp) :: step
   ! Step of the mesh

   real(dp),allocatable :: mesh(:)
   ! mesh(nw)

   real(dp),allocatable :: dos(:,:)
   ! dos(nw,0:nsppol)
   ! Total DOS, spin up and spin down component.

   real(dp),allocatable :: idos(:,:)
   ! idos(nw,0:nsppol)
   ! Integrated DOS, spin up and spin down component.

   real(dp),allocatable :: gef(:)
   ! gef(0:nsppol)
   ! DOS at the Fermi level. Total, spin up, spin down

   !type(ebands_t),pointer :: ebands => null()
   ! Reference to the bandstructure.

 end type edos_t

 public :: ebands_get_edos   ! Compute the dos from the band structure.
 public :: edos_free         ! Free memory
 public :: edos_write        ! Write results to file (formatted mode)
 public :: edos_print        ! Print DOS info to Fortran unit.
!!***

!----------------------------------------------------------------------

!!****t* m_ebands/gaps_t
!! NAME
!! gaps_t
!!
!! FUNCTION
!! Structure with information on the fundamental and optical gaps returned by ebands_report_gap.
!!
!! SOURCE

 type,public :: gaps_t

   integer :: nsppol
    ! Number of spins.

   integer,allocatable :: fo_kpos(:,:)
    ! fo_kpos(3,nsppol)
    ! fo_kpos(1:2,spin) ==> Indices of the k-points where the homo, lumo states are located (for each spin).
    ! fo_kpos(3,spin)   ==> the index of k-point where the optical gap is located (for each spin).

   integer,allocatable :: ierr(:)
     ! The third index corresponds to a "status" :
     !   0.0dp if gaps were not computed (because there are only valence bands);
     !  -1.0dp if the system (or spin-channel) is metallic;
     !   1.0dp if the gap was computed

   real(dp),allocatable :: fo_values(:,:)
     ! fo_values(2,nsppol)]
     ! Fundamental and optical gaps (in Hartree) for each spin.

   real(dp),pointer :: kpoints(:,:) => null()
     ! Reference to the k-points of the band structure used to compute the gaps.

   character(len=500),allocatable :: errmsg_spin(:)
     ! errmsg_spin(nsppol)
     ! String with human-readable error messages if ierr(spin) != 0.

 end type gaps_t

 public :: get_gaps      ! Build the object from a bandstructure.
 public :: gaps_free     ! Free the structure.
 public :: gaps_print    ! Print info on the gaps
!!***

!----------------------------------------------------------------------

!!****t* m_ebands/ebspl_t
!! NAME
!! ebspl_t
!!
!! FUNCTION
!!  B-spline interpolation of electronic eigenvalues.
!!
!! SOURCE

 type :: bcoeff_t
   real(dp),allocatable :: vals(:,:,:)
 end type bcoeff_t

 type,public :: ebspl_t

   integer :: nkx,nky,nkz
   ! Number of input data points

   integer :: kxord,kyord,kzord
   ! Order of the spline.

   !real(dp),allocatable :: xvec(:),yvec(:),zvec(:)
   real(dp),allocatable :: xknot(:),yknot(:),zknot(:)
   ! Array of length ndata+korder containing the knot

   type(bcoeff_t),allocatable :: coeff(:,:)
   ! coeff(mband, nsppol)

 end type ebspl_t

 public :: ebspl_new         ! Build B-spline object.
 public :: ebspl_evalk       ! Interpolate eigenvalues at an arbitrary k-point.
 public :: ebspl_free        ! Free memory.


CONTAINS  !=====================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_gaps
!! NAME
!! get_gaps
!!
!! FUNCTION
!!  Returns a structure with info on the fundamental and optical gap.
!!
!! INPUTS
!!  ebands<ebands_t>=Info on the band structure, the smearing technique and the physical temperature used.
!!  [kmask]=Logical mask used to exclude k-points.
!!
!! OUTPUT
!!  retcode=Return code (!=0 signals failure)
!!  gaps<gaps_t>=object with info on the gaps (parent is responsible for freeing the object).
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function get_gaps(ebands,gaps,kmask) result(retcode)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_gaps'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),target,intent(in)  :: ebands
 type(gaps_t),intent(out) :: gaps
!arrays
 logical,optional,intent(in) :: kmask(ebands%nkpt)

!Local variables-------------------------------
!scalars
 integer :: ikibz,nband_k,spin,nsppol,ikopt,ivk,ick,ivb,icb,retcode
 real(dp),parameter :: tol_fermi=tol6
 real(dp) :: fun_gap,opt_gap
 logical :: ismetal
!arrays
 integer :: val_idx(ebands%nkpt,ebands%nsppol)
 real(dp) :: top_valence(ebands%nkpt),bot_conduct(ebands%nkpt)
 logical :: my_kmask(ebands%nkpt)

! *********************************************************************

 nsppol = ebands%nsppol

 ! Initialize gaps_t
 gaps%nsppol = nsppol
 ABI_MALLOC(gaps%fo_kpos, (3,nsppol))
 ABI_MALLOC(gaps%ierr, (nsppol))
 ABI_MALLOC(gaps%fo_values, (2, nsppol))
 ABI_MALLOC(gaps%errmsg_spin, (nsppol))
 gaps%kpoints => ebands%kptns

 gaps%fo_kpos = 0
 gaps%ierr = 0
 gaps%fo_values = zero
 gaps%errmsg_spin(:) = ""

 my_kmask=.TRUE.; if (PRESENT(kmask)) my_kmask=kmask

 val_idx(:,:) = get_valence_idx(ebands,tol_fermi)

 spin_loop: &
&  do spin=1,nsppol

   ! No output if system i metallic
   ismetal=ANY(val_idx(:,spin)/=val_idx(1,spin))
   if (ismetal) then
     gaps%ierr(spin) = 1
     write(gaps%errmsg_spin(spin), "(a,i0)")"Metallic system for spin channel ",spin
     CYCLE
   endif

   ivb=val_idx(1,spin)
   icb=ivb+1

   do ikibz=1,ebands%nkpt
     if (.not.my_kmask(ikibz)) CYCLE
     nband_k=ebands%nband(ikibz+(spin-1)*ebands%nkpt)
     top_valence(ikibz)=ebands%eig(ivb,ikibz,spin)
     if (icb>nband_k) then
       gaps%ierr(spin) = 2
       gaps%errmsg_spin(spin) = "Not enough states to calculate the band gap."
       CYCLE spin_loop
     endif
     bot_conduct(ikibz)=ebands%eig(icb,ikibz,spin)
   end do

   ! Minimum of the optical Gaps
   ikopt= imin_loc(bot_conduct-top_valence,MASK=my_kmask)
   opt_gap=bot_conduct(ikopt)-top_valence(ikopt)

   ! Fundamental Gap ===
   ick = imin_loc(bot_conduct,MASK=my_kmask)
   ivk = imax_loc(top_valence,MASK=my_kmask)
   fun_gap = ebands%eig(icb,ick,spin)-ebands%eig(ivb,ivk,spin)

   gaps%fo_values(:,spin) = [fun_gap, opt_gap]
   gaps%fo_kpos(:,spin) = [ivk, ick, ikopt]
 end do spin_loop

 retcode = MAXVAL(gaps%ierr)

end function get_gaps
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
!!      setup_sigma
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine gaps_free(gaps)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gaps_free'
!End of the abilint section

 implicit none

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gaps_free'
!End of the abilint section

!Arguments ------------------------------------
 type(gaps_t),intent(inout) :: gaps

! *********************************************************************

 !@gaps_t

!integer
 if (allocated(gaps%fo_kpos)) then
   ABI_FREE(gaps%fo_kpos)
 end if
 if (allocated(gaps%ierr)) then
   ABI_FREE(gaps%ierr)
 end if

!real
 if (allocated(gaps%fo_values)) then
   ABI_FREE(gaps%fo_values)
 end if

!chars
 if (allocated(gaps%errmsg_spin)) then
   ABI_FREE(gaps%errmsg_spin)
 end if

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
!!  Print info on the fundamental and optical gap.
!!
!! INPUTS
!!  gaps<gaps_t>=Object with info on the gaps.
!!  [header]=Optional title.
!!  [unit]=Optional unit for output (std_out if not specified)
!!  [mode_paral]=Either "COLL" or "PERS", former is default.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      setup_sigma
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine gaps_print(gaps,header,unit,mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gaps_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: unit
 character(len=4),intent(in),optional :: mode_paral
 character(len=*),intent(in),optional :: header
 type(gaps_t),intent(in)  :: gaps

!Local variables-------------------------------
!scalars
 integer :: spin,ikopt,ivk,ick,my_unt
 real(dp) :: fun_gap,opt_gap
 character(len=4) :: my_mode
 character(len=500) :: msg

! *********************************************************************

 my_unt =std_out; if (PRESENT(unit      )) my_unt =unit
 my_mode='COLL' ; if (PRESENT(mode_paral)) my_mode=mode_paral

 do spin=1,gaps%nsppol

   if (spin==1) then
     msg=ch10
     if (PRESENT(header)) msg=ch10//' === '//TRIM(ADJUSTL(header))//' === '
     call wrtout(my_unt,msg,my_mode)
   end if

   if (gaps%ierr(spin) /= 0) then
     call wrtout(my_unt,gaps%errmsg_spin(spin), my_mode)
     continue
   end if

   ! Get minimum of the optical Gap.
   fun_gap = gaps%fo_values(1,spin)
   opt_gap = gaps%fo_values(2,spin)

   if (any(gaps%fo_kpos(:,spin) == 0)) then
     call wrtout(my_unt,sjoin("Cannot detect gap for spin: ",itoa(spin)),"COLL")
     cycle
   end if

   ivk = gaps%fo_kpos(1,spin)
   ick = gaps%fo_kpos(2,spin)
   ikopt = gaps%fo_kpos(3,spin)

   write(msg,'(a,i2,a,2(a,f8.4,a,3f8.4,a),33x,a,3f8.4)')&
&    '  >>>> For spin ',spin,ch10,&
&    '   Minimum optical gap = ',opt_gap*Ha_eV,' [eV], located at k-point      : ',gaps%kpoints(:,ikopt),ch10,&
&    '   Fundamental gap     = ',fun_gap*Ha_eV,' [eV], Top of valence bands at : ',gaps%kpoints(:,ivk),ch10,  &
&                                              '       Bottom of conduction at : ',gaps%kpoints(:,ick)
   call wrtout(my_unt,msg,my_mode)
 end do !spin

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
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine ebands_init(bantot,ebands,nelect,doccde,eig,istwfk,kptns,&
& nband,nkpt,npwarr,nsppol,nspinor,tphysel,tsmear,occopt,occ,wtk,&
& charge, kptopt, kptrlatt_orig, nshiftk_orig, shiftk_orig, kptrlatt, nshiftk, shiftk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_init'
!End of the abilint section

 implicit none

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_init'
!End of the abilint section

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
 ebands%eig=HUGE(one); ebands%occ=zero; ebands%doccde=zero

 call put_eneocc_vect(ebands,'eig',   eig   )
 call put_eneocc_vect(ebands,'occ',   occ   )
 call put_eneocc_vect(ebands,'doccde',doccde)

 ABI_MALLOC(ebands%wtk,(nkpt))
 ebands%wtk(1:nkpt)=wtk(1:nkpt)

!EBANDS_NEW
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
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

type(ebands_t) function ebands_from_hdr(hdr, mband, ene3d, nelect) result(ebands)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_from_hdr'
!End of the abilint section

 implicit none

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
&  hdr%npwarr,hdr%nsppol,hdr%nspinor,hdr%tphysel,hdr%tsmear,hdr%occopt,hdr%occ,hdr%wtk,&
&  hdr%charge, hdr%kptopt, hdr%kptrlatt_orig, hdr%nshiftk_orig, hdr%shiftk_orig, hdr%kptrlatt, hdr%nshiftk, hdr%shiftk)

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_from_dtset'
!End of the abilint section

 implicit none

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
!!  (only deallocate)
!!
!! PARENTS
!!      bethe_salpeter,dfpt_looppert,eig2tot,elphon,eph,gstate,m_ioarr,m_iowf
!!      m_shirley,m_wfk,mlwfovlp_qp,nonlinear,optic,outscfcv,respfn,screening
!!      sigma,wfk_analyze
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine ebands_free(ebands)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(inout) :: ebands
! *************************************************************************

 DBG_ENTER("COLL")

 if (allocated(ebands%istwfk)) then
   ABI_FREE(ebands%istwfk)
 end if
 if (allocated(ebands%nband)) then
   ABI_FREE(ebands%nband)
 end if
 if (allocated(ebands%npwarr)) then
   ABI_FREE(ebands%npwarr)
 end if
 if (allocated(ebands%kptns)) then
   ABI_FREE(ebands%kptns)
 end if
 if (allocated(ebands%eig)) then
   ABI_FREE(ebands%eig)
 end if
 if (allocated(ebands%occ)) then
   ABI_FREE(ebands%occ)
 end if
 if (allocated(ebands%doccde)) then
   ABI_FREE(ebands%doccde)
 end if
 if (allocated(ebands%wtk)) then
   ABI_FREE(ebands%wtk)
 end if

 if (allocated(ebands%shiftk_orig)) then
   ABI_FREE(ebands%shiftk_orig)
 end if
 if (allocated(ebands%shiftk)) then
   ABI_FREE(ebands%shiftk)
 end if

 DBG_EXIT("COLL")

end subroutine ebands_free
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_copy
!! NAME
!!  ebands_copy
!!
!! FUNCTION
!! This subroutine performs a deep copy of a ebands_t datatype.
!! All the associated pointers in the input object will be copied preserving the shape.
!! If a pointer in ibands happens to be not associated, the corresponding
!! pointer in the copied object will be nullified.
!!
!! INPUTS
!!  ibands<ebands_t>=The data type to be copied.
!!
!! OUTPUT
!!  obands<ebands_t>=The copy.
!!
!! PARENTS
!!      screening,setup_bse,setup_bse_interp,sigma
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine ebands_copy(ibands,obands)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_copy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in)  :: ibands
 type(ebands_t),intent(out) :: obands

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

end subroutine ebands_copy
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
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine ebands_print(ebands,header,unit,prtvol,mode_paral)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_print'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: prtvol,unit
 character(len=*),optional,intent(in) :: header
 character(len=4),optional,intent(in) :: mode_paral
 type(ebands_t),intent(in) :: ebands

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
   if (ebands%nsppol==1)then
     write(msg,'(a,i0,a)')' New occ. numbers for occopt= ',ebands%occopt,' , spin-unpolarized case. '
     call wrtout(my_unt,msg,my_mode)
   end if

   do spin=1,ebands%nsppol
     if (ebands%nsppol==2) then
       write(msg,'(a,i4,a,i2)')' New occ. numbers for occopt= ',ebands%occopt,' spin ',spin
       call wrtout(my_unt,msg,my_mode)
     end if

     do ikpt=1,ebands%nkpt
       write(msg,'(2a,i4,a,3f12.6,a,f6.3)')ch10,&
         ' k-point number ',ikpt,') ',ebands%kptns(:,ikpt),'; weight: ',ebands%wtk(ikpt)
       call wrtout(my_unt,msg,my_mode)
       do ii=1,ebands%nband(ikpt+(spin-1)*ebands%nkpt)
         write(msg,'(3(f7.3,1x))')ebands%eig(ii,ikpt,spin)*Ha_eV,ebands%occ(ii,ikpt,spin),ebands%doccde(ii,ikpt,spin)
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
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine unpack_eneocc(nkpt,nsppol,mband,nband,vect,array3d)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'unpack_eneocc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nsppol,mband
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: vect(:)
 real(dp),intent(out) :: array3d(mband,nkpt,nsppol)

!Local variables-------------------------------
 integer :: spin,ikpt,band,idx
! *************************************************************************

 array3d=HUGE(zero)

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
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine pack_eneocc(nkpt,nsppol,mband,nband,bantot,array3d,vect)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pack_eneocc'
!End of the abilint section

 implicit none

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
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine get_eneocc_vect(ebands,arr_name,vect)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_eneocc_vect'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: arr_name
 type(ebands_t),intent(in) :: ebands
 real(dp),intent(out) :: vect(ebands%bantot)

!Local variables-------------------------------
 integer :: nkpt,nsppol,mband,bantot
! *************************************************************************

 mband =ebands%mband; bantot=ebands%bantot; nkpt=ebands%nkpt; nsppol=ebands%nsppol

 SELECT CASE (arr_name)
 CASE ('occ')
   call pack_eneocc(nkpt,nsppol,mband,ebands%nband,bantot,ebands%occ,vect)
 CASE ('eig')
   call pack_eneocc(nkpt,nsppol,mband,ebands%nband,bantot,ebands%eig,vect)
 CASE ('doccde')
   call pack_eneocc(nkpt,nsppol,mband,ebands%nband,bantot,ebands%doccde,vect)
 CASE DEFAULT
   MSG_BUG(sjoin('Wrong arr_name= ', arr_name))
 END SELECT

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
!!      m_ebands
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine put_eneocc_vect(ebands,arr_name,vect)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'put_eneocc_vect'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: arr_name
 type(ebands_t),intent(inout) :: ebands
 real(dp),intent(in) :: vect(ebands%bantot)

!Local variables-------------------------------
 integer :: nkpt,nsppol,mband,bantot
! *************************************************************************

 mband =ebands%mband
 bantot=ebands%bantot
 nkpt  =ebands%nkpt
 nsppol=ebands%nsppol

 SELECT CASE (tolower(arr_name))
 CASE ('occ')
   call unpack_eneocc(nkpt,nsppol,mband,ebands%nband,vect,ebands%occ)
 CASE ('eig')
   call unpack_eneocc(nkpt,nsppol,mband,ebands%nband,vect,ebands%eig)
 CASE ('doccde')
   call unpack_eneocc(nkpt,nsppol,mband,ebands%nband,vect,ebands%doccde)
 CASE DEFAULT
   MSG_BUG(sjoin('Wrong arr_name= ', arr_name))
 END SELECT

end subroutine put_eneocc_vect
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_bandenergy
!! NAME
!! get_bandenergy
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

pure function get_bandenergy(ebands) result(band_energy)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_bandenergy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in) :: ebands
 real(dp) :: band_energy

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

end function get_bandenergy
!!***

!!****f* m_ebands/get_valence_idx
!! NAME
!!  get_valence_idx
!!
!! FUNCTION
!!  For each k-point and spin polarisation, report:
!!   The index of the valence in case of Semiconductors.
!!   The index of the band at the Fermi energy+toldfe
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

pure function get_valence_idx(ebands,tol_fermi) result(val_idx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_valence_idx'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),optional,intent(in) :: tol_fermi
 type(ebands_t),intent(in) :: ebands
!arrays
 integer :: val_idx(ebands%nkpt,ebands%nsppol)

!Local variables-------------------------------
 integer :: band,ikpt,spin,idx,nband_k
 real(dp) :: tol_

! *************************************************************************

 tol_=tol6; if (PRESENT(tol_fermi)) tol_=tol_fermi

 do spin=1,ebands%nsppol
   do ikpt=1,ebands%nkpt
     nband_k=ebands%nband(ikpt+(spin-1)*ebands%nkpt)

     idx=0
     do band=1,nband_k
       if (ebands%eig(band,ikpt,spin) > ebands%fermie+ABS(tol_)) then
         idx=band; EXIT
       end if
     end do
     val_idx(ikpt,spin)=idx-1
     if (idx==1) val_idx(ikpt,spin)=idx
     if (idx==0) val_idx(ikpt,spin)=nband_k

   end do
 end do

end function get_valence_idx
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/apply_scissor
!! NAME
!!  apply_scissor
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
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine apply_scissor(ebands,scissor_energy)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'apply_scissor'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: scissor_energy
 type(ebands_t),intent(inout) :: ebands

!Local variables-------------------------------
 integer :: ikpt,spin,ival,nband_k
 real(dp) :: spinmagntarget_
 character(len=500) :: msg
!arrays
 integer :: val_idx(ebands%nkpt,ebands%nsppol)
! *************************************************************************

 ! === Get the valence band index for each k and spin ===
 val_idx(:,:) = get_valence_idx(ebands)

 do spin=1,ebands%nsppol
   if (ANY(val_idx(:,spin)/=val_idx(1,spin))) then
     write(msg,'(a,i0,a)')&
      'Trying to apply a scissor operator on a metallic band structure for spin: ',spin,&
      'Assuming you know what you are doing, continuing anyway! '
     MSG_COMMENT(msg)
     !Likely newocc will stop, unless the system is semimetallic ?
   end if
 end do

 ! === Apply the scissor ===
 do spin=1,ebands%nsppol
   do ikpt=1,ebands%nkpt
     nband_k=ebands%nband(ikpt+(spin-1)*ebands%nkpt)
     ival=val_idx(ikpt,spin)

     if (nband_k>=ival+1) then
       ebands%eig(ival+1:,ikpt,spin) = ebands%eig(ival+1:,ikpt,spin)+scissor_energy
     else
       write(msg,'(2a,4(a,i0))')&
        'Not enough bands to apply the scissor operator. ',ch10,&
        'spin= ',spin,' ikpt= ',ikpt,' nband_k= ',nband_k,' but valence index= ',ival
       MSG_COMMENT(msg)
     end if

   end do
 end do

 ! === Recalculate the fermi level and occ. factors ===
 ! * For Semiconductors only the Fermi level is changed (in the middle of the new gap)
 spinmagntarget_=-99.99_dp !?; if (PRESENT(spinmagntarget)) spinmagntarget_=spinmagntarget
 call ebands_update_occ(ebands,spinmagntarget_)

end subroutine apply_scissor
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/get_occupied
!! NAME
!!  get_occupied
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

pure function get_occupied(ebands,tol_occ) result(occ_idx)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_occupied'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),optional,intent(in) :: tol_occ
 type(ebands_t),intent(in) :: ebands
!arrays
 integer :: occ_idx(ebands%nkpt,ebands%nsppol)

!Local variables-------------------------------
 integer :: band,ikpt,spin,idx,nband_k
 real(dp) :: tol_

! *************************************************************************

 tol_=tol8 ; if (PRESENT(tol_occ)) tol_=tol_occ

 do spin=1,ebands%nsppol
   do ikpt=1,ebands%nkpt
     nband_k=ebands%nband(ikpt+(spin-1)*ebands%nkpt)

     idx=0
     do band=1,nband_k
       if (ebands%occ(band,ikpt,spin)<ABS(tol_)) then
         idx=band; EXIT
       end if
     end do
     occ_idx(ikpt,spin)=idx-1
     if (idx==1) occ_idx(ikpt,spin)=idx
     if (idx==0) occ_idx(ikpt,spin)=nband_k

   end do
 end do

end function get_occupied
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/enclose_degbands
!! NAME
!!  enclose_degbands
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
!!
!! SIDE EFFECTS
!!  ibmin,ibmax=
!!    Input: initial guess for the indeces
!!    Output: All the denerate states are between ibmin and ibmax
!!
!! PARENTS
!!      setup_sigma
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine enclose_degbands(ebands,ikibz,spin,ibmin,ibmax,changed,tol_enedif)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'enclose_degbands'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikibz,spin
 integer,intent(inout) :: ibmin,ibmax
 real(dp),intent(in) :: tol_enedif
 logical,intent(out) :: changed
 type(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
!scalars
 integer :: ib,ibmin_bkp,ibmax_bkp
 real(dp) :: emin,emax

! *************************************************************************

 ibmin_bkp = ibmin; ibmax_bkp = ibmax

 emin =  ebands%eig(ibmin,ikibz,spin)
 do ib=ibmin-1,1,-1
   if ( ABS(ebands%eig(ib,ikibz,spin) - emin) > tol_enedif) then
     ibmin = ib +1
     EXIT
   else
     ibmin = ib
   end if
 end do

 emax =  ebands%eig(ibmax,ikibz,spin)
 do ib=ibmax+1,ebands%nband(ikibz+(spin-1)*ebands%nkpt)
   if ( ABS(ebands%eig(ib,ikibz,spin) - emax) > tol_enedif) then
     ibmax = ib - 1
     EXIT
   else
     ibmax = ib
   end if
 end do

 changed = (ibmin /= ibmin_bkp) .or. (ibmax /= ibmax_bkp)

end subroutine enclose_degbands
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_nelect_per_spin
!! NAME
!!  ebands_nelect_per_spin
!!
!! FUNCTION
!!   return number of electrons in each spin channel (computed from occoputation factors if nsppol=2)
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_nelect_per_spin'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in) :: ebands
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

!!****f* m_ebands/get_minmax
!! NAME
!!  get_minmax
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

function get_minmax(ebands,arr_name) result(minmax)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_minmax'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),target,intent(in) :: ebands
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

 SELECT CASE (tolower(arr_name))
 CASE ('occ')
   rdata => ebands%occ
 CASE ('eig')
   rdata => ebands%eig
 CASE ('doccde')
   rdata => ebands%doccde
 CASE DEFAULT
   MSG_BUG(sjoin('Wrong arr_name:', arr_name))
 END SELECT

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

end function get_minmax
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_edstats
!! NAME
!! ebands_edstats
!!
!! FUNCTION
!!  Compute statistical parameters of the energy differences e_ks[b+1] - e_ks[b]
!!  Returns stats_t record with the results (mean, stdev, min, max)
!!
!! INPUTS
!!  ebands<ebands_t>=Band energies.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(stats_t) function ebands_edstats(ebands) result(stats)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_edstats'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
!scalars
 integer :: ikibz,nband_k,spin
!arrays
 real(dp),allocatable :: ediffs(:,:,:),edvec(:)

! *************************************************************************

! Compute energy difference between b+1 and b.
 ABI_CALLOC(ediffs, (ebands%mband-1,ebands%nkpt,ebands%nsppol))

 do spin=1,ebands%nsppol
   do ikibz=1,ebands%nkpt
     nband_k=ebands%nband(ikibz+(spin-1)*ebands%nkpt)
     if (nband_k > 1) then
       ediffs(1:nband_k-1,ikibz,spin) = ebands%eig(2:nband_k,ikibz,spin) - ebands%eig(1:nband_k-1,ikibz,spin)
     end if
   end do
 end do

 ! Calculate the statistical parameters
 ! Not completely correct if nband_k varies but ...
 ABI_MALLOC(edvec, ((ebands%mband-1)*ebands%nkpt*ebands%nsppol))
 edvec = reshape(ediffs, [(ebands%mband-1)*ebands%nkpt*ebands%nsppol])

 stats = stats_eval(edvec)

 ABI_FREE(ediffs)
 ABI_FREE(edvec)

end function ebands_edstats
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_has_metal_scheme'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in) :: ebands

! *************************************************************************

 ans = (any(ebands%occopt == [3,4,5,6,7,8]))

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_write_bxsf'
 use interfaces_61_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: fname
 type(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: crystal

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
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine ebands_update_occ(ebands,spinmagntarget,stmbias,prtvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_update_occ'
 use interfaces_14_hidewrite
 use interfaces_61_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(inout) :: ebands
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
   !  If occupation is metallic have to compute new occupation numbers.
   if (my_prtvol > 10) then
     write(msg,'(a,f9.5)')' metallic scheme, calling newocc with spinmagntarget = ',spinmagntarget
     call wrtout(std_out,msg,'COLL')
   end if

   ! newocc assumes eigenvalues and occupations packed in 1d-vector!!
   mband  = ebands%mband
   nkpt   = ebands%nkpt
   nsppol = ebands%nsppol

   ABI_MALLOC(eigen,(mband*nkpt*nsppol))
   call get_eneocc_vect(ebands,'eig',eigen)

   ABI_MALLOC(occ,(mband*nkpt*nsppol))
   ABI_MALLOC(doccde,(mband*nkpt*nsppol))

   call newocc(doccde,eigen,entropy,fermie,spinmagntarget,mband,ebands%nband,&
&    ebands%nelect,ebands%nkpt,ebands%nspinor,ebands%nsppol,occ,ebands%occopt,&
&    my_prtvol,stmbias_local,ebands%tphysel,ebands%tsmear,ebands%wtk)

   ! Save output in ebands%.
   ebands%entropy = entropy
   ebands%fermie  = fermie
   call put_eneocc_vect(ebands,'occ'   ,occ)
   call put_eneocc_vect(ebands,'doccde',doccde)
   ABI_FREE(eigen)
   ABI_FREE(occ)
   ABI_FREE(doccde)

 else
   !  Semiconductor or Insulator.
   !
   ! FIXME here there is an inconsistency btw GW and Abinit
   ! In abinit Fermi is set to HOMO while in GW fermi is in the middle
   ! of Gap. In case of crystal systems, the later convention should be preferable.
   ! Anyway we have to decide and follow a unique convention to avoid problems.
   !
   ! occupation factors MUST be initialized
   if (ALL(ABS(ebands%occ) < tol6)) then
     msg = "occupation factors are not initialized, likely due to the use of iscf=-2"
     MSG_ERROR(msg)
   end if

   maxocc=two/(ebands%nsppol*ebands%nspinor)

   ! * Calculate the valence index for each spin channel.
   do spin=1,ebands%nsppol
     valencetop(spin)= smallest_real
     condbottom(spin)= greatest_real

     do ikibz=1,ebands%nkpt
       nband_k=ebands%nband(ikibz+(spin-1)*ebands%nkpt)
       do band=1,nband_k
         if (ebands%occ(band,ikibz,spin)/maxocc>one-tol6 .and. valencetop(spin)<ebands%eig(band,ikibz,spin)) then
           valencetop(spin)=ebands%eig(band,ikibz,spin)
         end if
         if (ebands%occ(band,ikibz,spin)/maxocc<tol6 .and. condbottom(spin)>ebands%eig(band,ikibz,spin)) then
           condbottom(spin)=ebands%eig(band,ikibz,spin)
         end if
       end do
     end do

   end do

   vtop=MAXVAL(valencetop)
   cbot=MINVAL(condbottom)

   write(msg,'(a,f6.2,2a,f6.2)')&
&    ' top of valence       [eV] ',vtop*Ha_eV,ch10,&
&    ' bottom of conduction [eV] ',cbot*Ha_eV
   call wrtout(std_out,msg,'COLL')
   if (ebands%nsppol==2) then
     if (ABS(vtop-MINVAL(valencetop))>tol6) then
       write(msg,'(a,i2)')' top of valence is spin ',MAXLOC(valencetop)
       call wrtout(std_out,msg,'COLL')
     end if
     if (ABS(cbot-MAXVAL(condbottom))>tol6) then
       write(msg,'(a,i2)')' bottom of conduction is spin ',MINLOC(condbottom)
       call wrtout(std_out,msg,'COLL')
     end if
   end if

   ! === Save output ===
   ! Here I dont know if it is better to be consistent with the abinit convention i.e fermi=vtop
   ebands%entropy=zero
   ebands%fermie=(vtop+cbot)/2
   if (ABS(cbot-vtop)<1.d-4) ebands%fermie=vtop ! To avoid error on the last digit FIXME is it really needed
 end if

 write(msg,'(a,f6.2,a)')' Fermi energy         [eV] ',ebands%fermie*Ha_eV,ch10
 call wrtout(std_out,msg,'COLL')
 !
 ! === Compute number of electrons for each spin channel ===
 nelect_spin(:)=zero
 do spin=1,ebands%nsppol
   do ikibz=1,ebands%nkpt
     nband_k=ebands%nband(ikibz+(spin-1)*ebands%nkpt)
     nelect_spin(spin)= nelect_spin(spin) + ebands%wtk(ikibz)*SUM(ebands%occ(1:nband_k,ikibz,spin))
   end do
 end do

 ndiff=ebands%nelect-SUM(nelect_spin)
 if (my_prtvol>0) then
   write(msg,'(2a,f6.2,2a,f7.4)')ch10,&
&    ' total number of electrons = ',SUM(nelect_spin),ch10,&
&    ' input and calculated no. of electrons differ by ',ndiff
   call wrtout(std_out,msg,'COLL')
 end if

 if (ABS(ndiff)>5.d-2*ebands%nelect) then
   write(msg,'(2a,2(a,f6.2))')&
&    'Too large difference in no. of electrons:,',ch10,&
&    'Expected= ',ebands%nelect,' Calculated= ',SUM(nelect_spin)
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
!! [prtvol]=Verbosity level (0 for lowest level)
!!
!! SIDE EFFECTS
!! ebands<ebands_t>=Info on the band structure, see above for side effects
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine ebands_set_scheme(ebands,occopt,tsmear,spinmagntarget,prtvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_set_scheme'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(inout) :: ebands
 integer,intent(in) :: occopt
 integer,optional,intent(in) :: prtvol
 real(dp),intent(in) :: tsmear,spinmagntarget

!Local variables-------------------------------
!scalars
 real(dp),parameter :: stmbias0=zero
 integer :: my_prtvol
 character(len=500) :: msg

! *************************************************************************

 my_prtvol = 0; if (present(prtvol)) my_prtvol = prtvol
 ebands%occopt = occopt; ebands%tsmear = tsmear

 if (prtvol > 10) then
   call wrtout(std_out, "Changing occupation scheme in electron bands")
   call wrtout(std_out, sjoin("occopt:", itoa(ebands%occopt), " ==> ", itoa(occopt)))
   call wrtout(std_out, sjoin("tsmear:", ftoa(ebands%tsmear), " ==> ", ftoa(tsmear)))
 end if

 call ebands_update_occ(ebands,spinmagntarget,stmbias0,prtvol=my_prtvol)

 if (prtvol > 10) then
   call wrtout(std_out, sjoin('Fermi level is now:', ftoa(ebands%fermie)))
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
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine ebands_set_fermie(ebands, fermie, msg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_set_fermie'
 use interfaces_61_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(inout) :: ebands
 real(dp),intent(in) :: fermie
 character(len=*),intent(out) :: msg

!Local variables-------------------------------
!scalars
 integer,parameter :: option1=1,unitdos0=0
 integer :: mband,nkpt,nsppol
 real(dp),parameter :: dosdeltae0=zero
 real(dp) :: prev_fermie,prev_nelect,maxocc
!arrays
 real(dp),allocatable :: doccde(:),occ(:),eigen(:)

! *************************************************************************

 if (ebands_has_metal_scheme(ebands)) then
   msg = "set_fermie assumes a metallic occupation scheme. Use ebands_set_scheme before calling ebands_set_ferme!"
   MSG_ERROR(msg)
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

 write(msg,"(2(a,es16.6),a,2(a,es16.6))")&
   " Old fermi level: ",prev_fermie,", with nelect: ",prev_nelect,ch10,&
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
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine ebands_set_nelect(ebands, nelect, spinmagntarget, msg, prtvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_set_nelect'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(inout) :: ebands
 integer,optional,intent(in) :: prtvol
 real(dp),intent(in) :: nelect,spinmagntarget
 character(len=*),intent(out) :: msg

!Local variables-------------------------------
!scalars
 integer :: my_prtvol
 real(dp) :: prev_fermie,prev_nelect

! *************************************************************************

 my_prtvol = 0; if (present(prtvol)) my_prtvol = prtvol

 if (ebands_has_metal_scheme(ebands)) then
   msg = "set_nelect assumes a metallic occupation scheme. Use ebands_set_scheme!"
   MSG_ERROR(msg)
 end if

 prev_fermie = ebands%fermie; prev_nelect = ebands%nelect
 ebands%nelect = nelect
 call ebands_update_occ(ebands,spinmagntarget,prtvol=my_prtvol)

 write(msg,"(2(a,es16.6),a,2(a,es16.6))")&
   "Old fermi level: ",prev_fermie,", with nelect: ",prev_nelect,ch10,&
   "New fermi level: ",ebands%fermie,", with nelect: ",ebands%nelect

end subroutine ebands_set_nelect
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_report_gap
!! NAME
!! ebands_report_gap
!!
!! FUNCTION
!!  Print info on the fundamental and optical gap.
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
!!  [gaps(3,nsppol)]=Fundamental and optical gaps. The third index corresponds to a "status":
!!      0.0dp if gaps were not computed (because there are only valence bands);
!!     -1.0dp if the system (or spin-channel) is metallic;
!!      1.0dp if the gap has been computed.
!!
!! PARENTS
!!      gstate,m_exc_diago,setup_bse,setup_bse_interp,setup_sigma,sigma
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine ebands_report_gap(ebands,header,kmask,unit,mode_paral,gaps)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_report_gap'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: unit
 character(len=4),intent(in),optional :: mode_paral
 character(len=*),intent(in),optional :: header
 type(ebands_t),intent(in)  :: ebands
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

 val_idx(:,:) = get_valence_idx(ebands,tol_fermi)
 first=0

!Initialize the return status for the gaps
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
     endif
     bot_conduct(ikibz) = ebands%eig(icb,ikibz,spin)
   end do

   ! === Get minimum of the optical Gap ===
   ikopt= imin_loc(bot_conduct-top_valence,MASK=my_kmask)
   opt_gap=bot_conduct(ikopt)-top_valence(ikopt)

   ! === Get fundamental Gap ===
   ick = imin_loc(bot_conduct,MASK=my_kmask)
   ivk = imax_loc(top_valence,MASK=my_kmask)
   fun_gap = ebands%eig(icb,ick,spin)-ebands%eig(ivb,ivk,spin)

   write(msg,'(a,i2,a,2(a,f8.4,a,3f8.4,a),33x,a,3f8.4)')&
&    '  >>>> For spin ',spin,ch10,&
&    '   Minimum optical gap = ',opt_gap*Ha_eV,' [eV], located at k-point      : ',ebands%kptns(:,ikopt),ch10,&
&    '   Fundamental gap     = ',fun_gap*Ha_eV,' [eV], Top of valence bands at : ',ebands%kptns(:,ivk),ch10,  &
&                                              '       Bottom of conduction at : ',ebands%kptns(:,ick)
   call wrtout(my_unt,msg,my_mode)

   if (PRESENT(gaps)) then
     gaps(:,spin) = (/fun_gap,opt_gap,one/)
   end if

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

integer function ebands_ncwrite(ebands,ncid) result(ncerr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_ncwrite'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 type(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
!scalars
#ifdef HAVE_NETCDF
 integer :: ii,nelect_int
 logical :: write_ngkpt
 character(len=etsfio_charlen) :: smearing,k_dependent
!arrays
 integer :: ngkpt(3)

! *************************************************************************

 select case (ebands%occopt)
 case (3)
   smearing = "Fermi-Dirac"
 case (4)
   smearing = "cold smearing of N. Marzari with minimization of the bump"
 case (5)
   smearing = "cold smearing of N. Marzari with monotonic function in the tail"
 case (6)
   smearing = "Methfessel and Paxton"
 case (7)
   smearing = "gaussian"
 case (8)
   smearing = "uniform"
 case default
   smearing = "none"
 end select

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
 ! FIXME: in etsf_io the number of electrons is declared as integer!!!
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
 NCF_CHECK(nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "tphysel", "charge"]))

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

 ! Write variables
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid("tphysel"), ebands%tphysel))
 NCF_CHECK(nf90_put_var(ncid, vid("occopt"), ebands%occopt))
 NCF_CHECK(nf90_put_var(ncid, vid("istwfk"), ebands%istwfk))
 NCF_CHECK(nf90_put_var(ncid, vid("kptopt"), ebands%kptopt))
 NCF_CHECK(nf90_put_var(ncid, vid("charge"), ebands%charge))
 NCF_CHECK(nf90_put_var(ncid, vid('kptrlatt_orig'), ebands%kptrlatt_orig))
 NCF_CHECK(nf90_put_var(ncid, vid('shiftk_orig'), ebands%shiftk_orig))
 NCF_CHECK(nf90_put_var(ncid, vid('kptrlatt'),ebands%kptrlatt))
 NCF_CHECK(nf90_put_var(ncid, vid('shiftk'), ebands%shiftk))

 if (write_ngkpt) then
   !write(std_out,*)"nshiftk_orig",nshiftk_orig,"shiftk_orig",shiftk_orig
   NCF_CHECK(nf90_put_var(ncid, vid('ngkpt_shiftk'), ebands%shiftk_orig))
 end if

#else
 MSG_ERROR("netcdf support is not activated. ")
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
!!  path=File name
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function ebands_ncwrite_path(ebands,path) result(ncerr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_ncwrite_path'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: path
 type(ebands_t),intent(in) :: ebands

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
   NCF_CHECK_MSG(ncerr, sjoin("creating", path))
 end if

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
!!  intmeth= 1 for gaussian, 2 for tetra
!!  step=Step on the linear mesh in Ha. If <0, the routine will use the mean of the energy level spacing
!!  broad=Gaussian broadening, If <0, the routine will use a default
!!    value for the broadening computed from the mean of the energy level spacing.
!!    No meaning if method == "tetra"
!!  comm=MPI communicator
!!
!! OUTPUT
!!  edos<edos_t>=Electronic DOS and IDOS.
!!
!! PARENTS
!!      eph
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

type(edos_t) function ebands_get_edos(ebands,cryst,intmeth,step,broad,comm) result(edos)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_get_edos'
 use interfaces_32_util
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: intmeth,comm
 real(dp),intent(in) :: step,broad
 type(ebands_t),target,intent(in)  :: ebands
 type(crystal_t),intent(in) :: cryst

!Local variables-------------------------------
!scalars
 integer,parameter :: bcorr0=0
 integer :: iw,nw,spin,band,ikpt,ief,nbz,nibz,nproc,my_rank,mpierr,cnt,ierr
 real(dp) :: max_ene,min_ene,wtk,max_occ
 character(len=500) :: msg
 character(len=80) :: errstr
 type(stats_t) :: ediffs
 type(t_tetrahedron) :: tetra
!arrays
 real(dp) :: eminmax_spin(2,ebands%nsppol)
 real(dp) :: qlatt(3,3),rlatt(3,3)
 integer,allocatable :: bz2ibz(:)
 real(dp),allocatable :: wme0(:),fullbz(:,:)
 real(dp),allocatable :: tmp_eigen(:),bdelta(:,:),btheta(:,:)

! *********************************************************************

 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 edos%nkibz = ebands%nkpt; edos%intmeth = intmeth; edos%nsppol = ebands%nsppol

 ! Compute the mean value of the energy spacing.
 ediffs = ebands_edstats(ebands)
 edos%broad = broad; if (broad <= tol16) edos%broad = ediffs%mean
 edos%step = step; if (step <= tol16) edos%step = 0.1 * ediffs%mean

 ! Compute the linear mesh so that it encloses all bands.
 eminmax_spin = get_minmax(ebands, "eig")
 min_ene = minval(eminmax_spin(1,:)); min_ene = min_ene - 0.1_dp * abs(min_ene)
 max_ene = maxval(eminmax_spin(2,:)); max_ene = max_ene + 0.1_dp * abs(max_ene)

 nw = nint((max_ene - min_ene)/edos%step) + 1; edos%nw = nw

 ABI_MALLOC(edos%mesh, (nw))
 edos%mesh = arth(min_ene, edos%step, nw)

 ABI_CALLOC(edos%gef, (0:edos%nsppol))
 ABI_CALLOC(edos%dos,  (nw, 0:edos%nsppol))
 ABI_CALLOC(edos%idos, (nw, 0:edos%nsppol))

 select case (intmeth)
 case (1)
   ! Gaussian
   ABI_MALLOC(wme0, (nw))
   cnt = 0
   do spin=1,edos%nsppol
     do ikpt=1,ebands%nkpt
       cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
       wtk = ebands%wtk(ikpt)
       do band=1,ebands%nband(ikpt+(spin-1)*ebands%nkpt)
          wme0 = edos%mesh - ebands%eig(band, ikpt, spin)
          edos%dos(:, spin) = edos%dos(:, spin) + wtk * dirac_delta(wme0, edos%broad)
       end do
     end do
   end do
   ABI_FREE(wme0)
   call xmpi_sum(edos%dos, comm, mpierr)

 case (2)
   ! Consistency test: return immediately if cannot use tetra.
   if (ebands%nkpt<2) then
     MSG_WARNING('at least 2 points are needed for tetrahedrons')
     ierr = ierr + 1
   end if
   if (any(ebands%nband /= ebands%nband(1)) ) then
     MSG_WARNING('for tetrahedrons, nband(:) must be constant')
     ierr = ierr + 1
   end if
   if (ebands%nshiftk > 1) then
     MSG_WARNING(sjoin("for tetrahedrons, nshift must be (0,1) but found:", itoa(ebands%nshiftk)))
     ierr = ierr + 1
   end if
   if (ierr/=0) MSG_ERROR("Aborting now")

   ! convert kptrlatt to double and invert.
   rlatt = dble(ebands%kptrlatt)
   call matr3inv(rlatt,qlatt)  ! qlatt refers to the shortest qpt vectors.

   ! Calculate total number of k-points in the full BZ
   nbz = &
      ebands%kptrlatt(1,1)*ebands%kptrlatt(2,2)*ebands%kptrlatt(3,3) &
     +ebands%kptrlatt(1,2)*ebands%kptrlatt(2,3)*ebands%kptrlatt(3,1) &
     +ebands%kptrlatt(1,3)*ebands%kptrlatt(2,1)*ebands%kptrlatt(3,2) &
     -ebands%kptrlatt(1,2)*ebands%kptrlatt(2,1)*ebands%kptrlatt(3,3) &
     -ebands%kptrlatt(1,3)*ebands%kptrlatt(2,2)*ebands%kptrlatt(3,1) &
     -ebands%kptrlatt(1,1)*ebands%kptrlatt(2,3)*ebands%kptrlatt(3,2)
   nbz = nbz * ebands%nshiftk

   ABI_MALLOC(bz2ibz,(nbz))
   ABI_MALLOC(fullbz,(3, nbz))
   nibz = ebands%nkpt

   ! === Make full kpoint grid and get equivalence to irred kpoints ===
   ! * Note: This routines scales badly wrt Kmesh%nbz
   ! TODO should be rewritten, pass kptopt and test whether listkk is faster.
   !ABI_CHECK(Kmesh%kptopt==1,"get_full_kgrid assumes kptopt==1")
   call get_full_kgrid(bz2ibz,ebands%kptns,fullbz,ebands%kptrlatt,nibz,nbz,ebands%nshiftk,&
     cryst%nsym,ebands%shiftk,cryst%symrel)

   if (ierr == 0) then
     call init_tetra(bz2ibz, cryst%gprimd, qlatt, fullbz, nbz, tetra, ierr, errstr)
     if (ierr/=0) MSG_ERROR(errstr)
   end if
   ABI_FREE(bz2ibz)
   ABI_FREE(fullbz)

   !call tetra_from_kptrlatt(tetra, cryst, ebands%kptopt, ebands%kptrlatt, &
   !  ebands%nshiftk, ebands%shiftk, ebands%nkpt, ebands%kptns)

   ! For each spin and band, interpolate over kpoints,
   ! calculate integration weights and DOS contribution.
   ABI_MALLOC(tmp_eigen, (ebands%nkpt))
   ABI_MALLOC(btheta, (nw, ebands%nkpt))
   ABI_MALLOC(bdelta, (nw, ebands%nkpt))

   cnt = 0
   do spin=1,ebands%nsppol
     do band=1,ebands%nband(1)
       cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
       ! For each band get its contribution
       tmp_eigen = ebands%eig(band,:,spin)

       ! Calculate integration weights at each irred k-point (Blochl et al PRB 49 16223)
       call tetra_blochl_weights(tetra,tmp_eigen,min_ene,max_ene,one,nw,ebands%nkpt,bcorr0,&
         btheta,bdelta,xmpi_comm_self)

       do ikpt=1,ebands%nkpt
         do iw=1,nw
           edos%dos(iw,spin) = edos%dos(iw,spin) + bdelta(iw, ikpt)
           ! IDOS is computed afterwards with simpson
           !edos%idos(iw,spin) = edos%idos(iw,spin) + btheta(iw,ikpt)
         end do
       end do
     end do ! band
   end do !spin

   call xmpi_sum(edos%dos, comm, mpierr)

   ! Free memory
   ABI_FREE(tmp_eigen)
   ABI_FREE(btheta)
   ABI_FREE(bdelta)

   call destroy_tetra(tetra)

   ! Filter so that dos[i] is always >= 0 and idos is monotonic
   ! IDOS is computed afterwards with simpson
   where (edos%dos(:,1:) <= zero) edos%dos(:,1:) = zero

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
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine edos_free(edos)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'edos_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(edos_t),intent(inout) :: edos

! *********************************************************************

 !@edos_t
!real
 if (allocated(edos%mesh)) then
   ABI_FREE(edos%mesh)
 end if
 if (allocated(edos%dos)) then
   ABI_FREE(edos%dos)
 end if
 if (allocated(edos%idos)) then
   ABI_FREE(edos%idos)
 end if
 if (allocated(edos%gef)) then
   ABI_FREE(edos%gef)
 end if

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
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine edos_write(edos, path)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'edos_write'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=*),intent(in) :: path
 type(edos_t),intent(in) :: edos

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

subroutine edos_print(edos, unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'edos_print'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(edos_t),intent(in) :: edos
 integer,optional,intent(in) :: unit

!Local variables-------------------------------
 integer :: unt

! *************************************************************************

 unt = std_out; if (present(unit)) unt = unit

 write(unt,'(a,es16.8,a)')' Fermi level: ',edos%mesh(edos%ief)*Ha_eV," [eV]"
 write(unt,"(a,es16.8)")" Total electron DOS in states/eV : ",edos%gef(0) / Ha_eV
 if (edos%nsppol == 2) then
   write(unt,"(a,es16.8)")"   Spin up:  ",edos%gef(1) / Ha_eV
   write(unt,"(a,es16.8)")"   Spin down:",edos%gef(2) / Ha_eV
 end if

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

integer function ebands_write_nesting(ebands,cryst,filepath,prtnest,tsmear,fermie_nest,&
  qpath_vertices,errmsg) result(skipnest)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_write_nesting'
!End of the abilint section

 implicit none

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_write_nesting'
!End of the abilint section

!Arguments ------------------------------------
 type(ebands_t),intent(in) :: ebands
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

 skipnest = 0; errmsg = ""
 if (any(ebands%nband /= ebands%nband(1))) then
   errmsg = 'mknesting can not handle variable nband(1:nkpt). Skipped.'//&
     ch10//' Correct input file to get nesting output'
   skipnest = 1; return
 end if

 if (ebands%nshiftk /= 1) then
   errmsg = 'mknesting does not support nshiftk > 1. Change ngkpt and shiftk to have only one shift after inkpts'
   skipnest = 1; return
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
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine ebands_expandk(inb, cryst, ecut_eff, force_istwfk1, dksqmax, bz2ibz, outb)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_expandk'
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: ecut_eff
 real(dp),intent(out) :: dksqmax
 logical,intent(in) :: force_istwfk1
 type(ebands_t),intent(in) :: inb
 type(ebands_t),intent(out) :: outb
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,allocatable,intent(out) :: bz2ibz(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: brav1=1,option0=0,istwfk_1=1
 integer :: mkpt,nkfull,timrev,bantot,sppoldbl,npw_k,nsppol,istw
 integer :: ik_ibz,ikf,isym,itimrev,spin,mband
 logical :: isirred_k
 !character(len=500) :: msg
!arrays
 integer :: g0(3)
 integer,allocatable :: istwfk(:),nband(:,:),npwarr(:),kg_k(:,:)
 real(dp),allocatable :: kfull(:,:),doccde(:),eig(:),occ(:),wtk(:)
 real(dp),allocatable :: doccde_3d(:,:,:),eig_3d(:,:,:),occ_3d(:,:,:)

! *********************************************************************

 ABI_CHECK(inb%kptopt /= 0, "ebands_expandk does not support kptopt == 0")

 nsppol = inb%nsppol

 ! Call smpbz to get the full grid of k-points `kfull`
 ! brav1=1 is able to treat all bravais lattices (same option used in getkgrid)
 mkpt= &
    inb%kptrlatt(1,1)*inb%kptrlatt(2,2)*inb%kptrlatt(3,3) &
   +inb%kptrlatt(1,2)*inb%kptrlatt(2,3)*inb%kptrlatt(3,1) &
   +inb%kptrlatt(1,3)*inb%kptrlatt(2,1)*inb%kptrlatt(3,2) &
   -inb%kptrlatt(1,2)*inb%kptrlatt(2,1)*inb%kptrlatt(3,3) &
   -inb%kptrlatt(1,3)*inb%kptrlatt(2,2)*inb%kptrlatt(3,1) &
   -inb%kptrlatt(1,1)*inb%kptrlatt(2,3)*inb%kptrlatt(3,2)
 mkpt = mkpt * inb%nshiftk

 ABI_MALLOC(kfull, (3,mkpt))
 call smpbz(brav1,std_out,inb%kptrlatt,mkpt,nkfull,inb%nshiftk,option0,inb%shiftk,kfull)

 ! Costruct full BZ and create mapping BZ --> IBZ
 ! Note:
 !   - we don't change the value of nsppol hence sppoldbl is set to 1
 !   - we use symrel so that bz2ibz can be used to reconstruct the wavefunctions.
 !
 sppoldbl = 1 !; if (any(cryst%symafm == -1) .and. inb%nsppol == 1) sppoldbl=2
 ABI_MALLOC(bz2ibz, (nkfull*sppoldbl,6))

 timrev = 1; if (any(inb%kptopt == [3, 4])) timrev = 0
 call listkk(dksqmax,cryst%gmet,bz2ibz,inb%kptns,kfull,inb%nkpt,nkfull,cryst%nsym,&
   sppoldbl,cryst%symafm,cryst%symrel,timrev,use_symrec=.False.)
   !sppoldbl,cryst%symafm,cryst%symrec,timrev,use_symrec=.True.)

 ABI_MALLOC(wtk, (nkfull))
 wtk = one/nkfull ! weights normalized to one

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
     eig_3d(:,ikf,spin) = inb%eig   (:,ik_ibz,spin)
     occ_3d(:,ikf,spin) = inb%occ   (:,ik_ibz,spin)
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
   inb%charge, inb%kptopt, inb%kptrlatt_orig, inb%nshiftk_orig, inb%shiftk_orig, inb%kptrlatt, inb%nshiftk, inb%shiftk)

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

!!****f* m_ebands/ebspl_new
!! NAME
!! ebspl_new
!!
!! FUNCTION
!! Build the `ebspl_t` object used to interpolate the band structure.
!!
!! INPUTS
!!  ords(3)=order of the spline for the three directions. ord(1) must be in [0, nkx] where
!!    nkx is the number of points along the x-axis.
!!  band_block(2)=Initial and final band index. If [0,0], all bands are used
!!
!! OUTPUT
!!
!! PARENTS
!!      m_ebands
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

type(ebspl_t) function ebspl_new(ebands, cryst, ords, band_block) result(ebspl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebspl_new'
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: ords(3), band_block(2)

!Local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1=1
 integer :: kxord,kyord,kzord,nxknot,nyknot,nzknot,ierr,nkfull,ikf
 integer :: spin,band,ik_ibz,timrev,ix,iy,iz,nkx,nky,nkz
 real(dp) :: dksqmax
 character(len=500) :: msg
!arrays
 integer :: ngkpt(3)
 integer,allocatable :: bz2ibz(:,:)
 real(dp),allocatable :: xvec(:),yvec(:),zvec(:),xyzdata(:,:,:),kfull(:,:)

! *********************************************************************

 ! Check input parameters
 ierr = 0
 if (ebands%nkpt == 1) then
   MSG_WARNING("Cannot interpolate with a single k-point")
   ierr = ierr + 1
 end if
 if (.not. isdiagmat(ebands%kptrlatt)) then
   MSG_WARNING('kptrlatt is not diagonal. Multiple shifts are not allowed')
   ierr = ierr + 1
 end if
 if (ebands%nshiftk /= 1) then
   MSG_WARNING('Multiple shifts not allowed')
   ierr = ierr + 1
 end if
 if (any(ebands%nband(:) /= ebands%nband(1))) then
   MSG_WARNING("nband must be constant")
   ierr = ierr + 1
 end if
 if (ierr /= 0) then
   MSG_WARNING("bspline interpolation cannot be performed. See warnings above. Returning")
   return
 end if

 ! Build BZ mesh Note that:
 ! 1) k-point coordinates are in [0, 1]
 ! 2) The mesh is closed e.g. (0,0,0) and (1,1,1) are included
 ngkpt(1)=ebands%kptrlatt(1,1)
 ngkpt(2)=ebands%kptrlatt(2,2)
 ngkpt(3)=ebands%kptrlatt(3,3)

 nkx = ngkpt(1)+1; nky = ngkpt(2)+1; nkz = ngkpt(3)+1
 ABI_MALLOC(xvec, (nkx))
 ABI_MALLOC(yvec, (nky))
 ABI_MALLOC(zvec, (nkz))

 ! Multiple shifts are not supported here.
 do ix=1,nkx
   xvec(ix) = (ix-one+ebands%shiftk(1,1)) / ngkpt(1)
 end do
 do iy=1,nky
   yvec(iy) = (iy-one+ebands%shiftk(2,1)) / ngkpt(2)
 end do
 do iz=1,nkz
   zvec(iz) = (iz-one+ebands%shiftk(3,1)) / ngkpt(3)
 end do

 ! Build list of k-points in full BZ (ordered as required by B-spline routines)
 nkfull = nkx*nky*nkz
 ABI_MALLOC(kfull, (3,nkfull))
 ikf = 0
 do iz=1,nkz
   do iy=1,nky
     do ix=1,nkx
       ikf = ikf + 1
       kfull(:,ikf) = [xvec(ix), yvec(iy), zvec(iz)]
     end do
   end do
 end do

 ! Build mapping kfull --> IBZ
 ABI_MALLOC(bz2ibz, (nkfull*sppoldbl1,6))

 timrev = 1; if (any(ebands%kptopt == [3, 4])) timrev = 0
 call listkk(dksqmax,cryst%gmet,bz2ibz,ebands%kptns,kfull,ebands%nkpt,nkfull,cryst%nsym,&
   sppoldbl1,cryst%symafm,cryst%symrec,timrev,use_symrec=.True.)
 ABI_FREE(kfull)

 if (dksqmax > tol12) then
   write(msg, '(3a,es16.6,4a)' )&
   'At least one of the k points could not be generated from a symmetrical one.',ch10,&
   'dksqmax=',dksqmax,ch10,&
   'Action: check k-point input variables',ch10,&
   '        (e.g. kptopt or shiftk might be wrong in the present dataset or the preparatory one.'
   MSG_ERROR(msg)
 end if

 ! Generate knots (ords is input)
 kxord = ords(1); kyord = ords(2); kzord = ords(3)
 nxknot = nkx + kxord
 nyknot = nky + kyord
 nzknot = nkz + kzord

 ebspl%nkx = nkx; ebspl%kxord = kxord
 ebspl%nky = nky; ebspl%kyord = kyord
 ebspl%nkz = nkz; ebspl%kzord = kzord

 ABI_MALLOC(ebspl%xknot,(nxknot))
 ABI_MALLOC(ebspl%yknot,(nyknot))
 ABI_MALLOC(ebspl%zknot,(nzknot))

 call dbsnak(nkx, xvec, kxord, ebspl%xknot)
 call dbsnak(nky, yvec, kyord, ebspl%yknot)
 call dbsnak(nkz, zvec, kzord, ebspl%zknot)

 ABI_MALLOC(xyzdata,(nkx,nky,nkz))
 ABI_DT_MALLOC(ebspl%coeff, (ebands%mband,ebands%nsppol))

 do spin=1,ebands%nsppol
   do band=1,ebands%mband
     if (all(band_block /= 0)) then
       if (band < band_block(1) .or. band > band_block(2)) cycle
     end if

     ABI_MALLOC(ebspl%coeff(band,spin)%vals, (nkx,nky,nkz))

     ! Build array in full bz to prepare call to dbs3in.
     ikf = 0
     do iz=1,nkz
       do iy=1,nky
         do ix=1,nkx
           ikf = ikf + 1
           ik_ibz = bz2ibz(ikf,1)
           xyzdata(ix,iy,iz) = ebands%eig(band,ik_ibz,spin)
         end do
       end do
     end do

     ! Construct 3D tensor for B-spline. Results in coeff(band,spin)%vals
     call dbs3in(nkx,xvec,nky,yvec,nkz,zvec,xyzdata,nkx,nky,kxord,kyord,kzord,ebspl%xknot,ebspl%yknot,ebspl%zknot,&
        ebspl%coeff(band,spin)%vals)
   end do
 end do

 ABI_FREE(xvec)
 ABI_FREE(yvec)
 ABI_FREE(zvec)
 ABI_FREE(bz2ibz)
 ABI_FREE(xyzdata)

end function ebspl_new
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebspl_evalk
!! NAME
!! ebspl_evalk
!!
!! FUNCTION
!!   Interpolate eigenvalues at an arbitrary k-point.
!!
!! INPUTS
!!  band_block(2)=Initial and final band index.
!!  band_block(2)=Index of the first and last band defining the block of states to be interpolated
!!  kpt(3)=K-point in reduced coordinate (will be wrapped in the interval [0,1[
!!  spin=Spin index
!!
!! OUTPUT
!!  oeig=Array of size band_block(2) - band_block(1) + 1 with the interpolated eigenvalues.
!!  [oder1]
!!
!! PARENTS
!!      m_ebands
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine ebspl_evalk(ebspl, band_block, kpt, spin, oeig, oder1)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebspl_evalk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin
 type(ebspl_t),intent(in) :: ebspl
!arrays
 integer,intent(in) :: band_block(2)
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(out) :: oeig(:)
 real(dp),optional,intent(out) :: oder1(:,:)

!Local variables-------------------------------
!scalars
 integer :: band,ib,ii
!arrays
 integer :: iders(3)
 real(dp) :: kred(3),shift(3)

! *********************************************************************

 ABI_CHECK(size(oeig) >= (band_block(2) - band_block(1) + 1), "oeig too small")

 ! Wrap k-point in the interval [0,1[ where 1 is not included (tol12)
 call wrap2_zero_one(kpt, kred, shift)

 ib = 0
 do band=band_block(1),band_block(2)
   ib = ib +1
   ! B-spline interpolation.
   oeig(ib) = dbs3vl(kred(1), kred(2), kred(3), ebspl%kxord, ebspl%kyord, ebspl%kzord,&
                     ebspl%xknot, ebspl%yknot, ebspl%zknot, ebspl%nkx, ebspl%nky, ebspl%nkz,&
                     ebspl%coeff(band,spin)%vals)

   if (present(oder1)) then
     ! Compute firs-order derivatives.
     do ii=1,3
       iders = 0; iders(ii) = 1
       oder1(ii,ib) = dbs3dr(iders(1), iders(2), iders(3), &
                             kred(1), kred(2), kred(3), ebspl%kxord, ebspl%kyord, ebspl%kzord,&
                             ebspl%xknot, ebspl%yknot, ebspl%zknot, ebspl%nkx, ebspl%nky, ebspl%nkz,&
                             ebspl%coeff(band,spin)%vals)
     end do
   end if

 end do

end subroutine ebspl_evalk
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebspl_free
!! NAME
!! ebspl_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! PARENTS
!!      m_ebands
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine ebspl_free(ebspl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebspl_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebspl_t),intent(inout) :: ebspl

!Local variables-------------------------------
!scalars
 integer :: ii,jj

! *********************************************************************

 if (allocated(ebspl%xknot)) then
   ABI_FREE(ebspl%xknot)
 end if
 if (allocated(ebspl%yknot)) then
   ABI_FREE(ebspl%yknot)
 end if
 if (allocated(ebspl%zknot)) then
   ABI_FREE(ebspl%zknot)
 end if

 ! Free B-spline coefficients.
 if (allocated(ebspl%coeff)) then
   do jj=1,size(ebspl%coeff, dim=2)
     do ii=1,size(ebspl%coeff, dim=1)
       if (allocated(ebspl%coeff(ii,jj)%vals)) then
         ABI_FREE(ebspl%coeff(ii,jj)%vals)
       end if
     end do
   end do
   ABI_DT_FREE(ebspl%coeff)
 end if

end subroutine ebspl_free
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_bspline
!! NAME
!! ebands_bspline
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

type(ebands_t) function ebands_bspline(ebands, cryst, ords, new_kptrlatt, new_nshiftk, new_shiftk, comm) result(new)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_bspline'
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 integer,intent(inout) :: new_nshiftk
 type(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: ords(3)
 integer,intent(inout) :: new_kptrlatt(3,3)
 real(dp),intent(inout) :: new_shiftk(3,new_nshiftk)

!Local variables-------------------------------
!scalars
 integer,parameter :: iout0=0,chksymbreak0=0,iscf2=2
 integer :: ik_ibz,spin,new_bantot,new_nkpt,nsppol,new_mband,nkpt_computed,kptopt
 integer :: nprocs,my_rank,cnt,ierr
 real(dp) :: kptrlen
 type(ebspl_t) :: ebspl
!arrays
 integer,parameter :: vacuum0(3)=[0,0,0]
 integer :: kptrlatt_orig(3,3),band_block0(2)
 integer,allocatable :: new_istwfk(:),new_nband(:,:),new_npwarr(:)
 real(dp) :: mynew_shiftk(3,210)
 real(dp),allocatable :: new_kpts(:,:),new_doccde(:),new_eig(:),new_occ(:),new_wtk(:)

! *********************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 nsppol = ebands%nsppol; kptrlatt_orig = ebands%kptrlatt; kptopt = ebands%kptopt
 band_block0 = [1, ebands%mband]

 ! First call to getkgrid to obtain the number of new_kpts.
 ! TODO: write wrapper
 ABI_MALLOC(new_kpts, (3,0))
 ABI_MALLOC(new_wtk, (0))

 ! Be careful as getkgrid expects shiftk(3,8).
 mynew_shiftk = zero; mynew_shiftk(:,1:new_nshiftk) = new_shiftk
 ABI_CHECK(new_nshiftk > 0 .and. new_nshiftk <=210, "new_nshiftk must be in [1,210]")

 call getkgrid(chksymbreak0,iout0,iscf2,new_kpts,kptopt,new_kptrlatt,kptrlen,&
   cryst%nsym,0,new_nkpt,new_nshiftk,cryst%nsym,cryst%rprimd,mynew_shiftk,cryst%symafm,cryst%symrel,vacuum0,new_wtk)

 ABI_FREE(new_kpts)
 ABI_FREE(new_wtk)

 ! Recall getkgrid to get new_kpts and new_wtk.
 ABI_MALLOC(new_kpts,(3,new_nkpt))
 ABI_MALLOC(new_wtk,(new_nkpt))

 call getkgrid(chksymbreak0,iout0,iscf2,new_kpts,kptopt,new_kptrlatt,kptrlen,&
   cryst%nsym,new_nkpt,nkpt_computed,new_nshiftk,cryst%nsym,cryst%rprimd,mynew_shiftk,&
   cryst%symafm,cryst%symrel,vacuum0,new_wtk)
 new_shiftk = mynew_shiftk(:,1:new_nshiftk)

 ! Initialize new ebands_t in new IBZ
 ABI_MALLOC(new_istwfk, (new_nkpt))
 new_istwfk = 1
 ABI_MALLOC(new_nband, (new_nkpt, nsppol))
 new_nband = ebands%mband
 ABI_MALLOC(new_npwarr, (new_nkpt))
 new_npwarr = maxval(ebands%npwarr)
 new_bantot = sum(new_nband); new_mband = maxval(new_nband)
 ABI_MALLOC(new_doccde, (new_bantot))
 ABI_MALLOC(new_eig, (new_bantot))
 ABI_MALLOC(new_occ, (new_bantot))

 call ebands_init(new_bantot,new,ebands%nelect,new_doccde,new_eig,new_istwfk,new_kpts,&
   new_nband,new_nkpt,new_npwarr,ebands%nsppol,ebands%nspinor,ebands%tphysel,ebands%tsmear,&
   ebands%occopt,new_occ,new_wtk,&
   ebands%charge, kptopt, kptrlatt_orig, ebands%nshiftk, ebands%shiftk, new_kptrlatt, new_nshiftk, new_shiftk)

 ABI_FREE(new_kpts)
 ABI_FREE(new_wtk)
 ABI_FREE(new_istwfk)
 ABI_FREE(new_nband)
 ABI_FREE(new_npwarr)
 ABI_FREE(new_doccde)
 ABI_FREE(new_eig)
 ABI_FREE(new_occ)

 ! Build B-spline object.
 ebspl = ebspl_new(ebands, cryst, ords, band_block0)

 ! Spline eigenvalues.
 new%eig = zero; cnt = 0
 do spin=1,new%nsppol
   do ik_ibz=1,new%nkpt
     cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle  ! Mpi parallelism.
     call ebspl_evalk(ebspl, band_block0, new%kptns(:,ik_ibz), spin, new%eig(:,ik_ibz,spin))
   end do
 end do
 call xmpi_sum(new%eig, comm, ierr)

 call ebspl_free(ebspl)

end function ebands_bspline
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_jdos
!! NAME
!! ebands_jdos
!!
!! FUNCTION
!!  Compute the joint density of states.
!!
!! INPUTS
!!  ebands<ebands_t>=Band structure object.
!!  cryst<cryst_t>=Info on the crystalline structure.
!!  intmeth= 1 for gaussian, 2 for tetra
!!  step=Step on the linear mesh in Ha. If <0, the routine will use the mean of the energy level spacing
!!  broad=Gaussian broadening, If <0, the routine will use a default
!!    value for the broadening computed from the mean of the energy level spacing.
!!    No meaning if tetra method
!!  comm=MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      destroy_tetra,get_full_kgrid,init_tetra,matr3inv,tetra_blochl_weights
!!      xmpi_sum
!!
!! SOURCE

subroutine ebands_jdos(ebands, cryst, intmeth, step, broad, comm, ierr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_jdos'
 use interfaces_32_util
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: intmeth,comm
 integer,intent(out) :: ierr
 real(dp),intent(in) :: step,broad
 type(ebands_t),intent(in) :: ebands
 type(crystal_t),intent(in) :: cryst

!Local variables-------------------------------
!scalars
 integer,parameter :: bcorr0=0
 integer :: ik_ibz,ibc,ibv,spin,iw,nw,nband_k,nbv,nproc,my_rank,nbz,nibz,cnt,mpierr,unt
 real(dp) :: wtk,wmax,wstep,wbroad
 type(stats_t) :: ediffs
 type(t_tetrahedron) :: tetra
 character(len=500) :: msg
 character(len=80) :: errstr
 character(len=fnlen) :: path
!arrays
 integer :: val_idx(ebands%nkpt,ebands%nsppol)
 integer,allocatable :: bz2ibz(:)
 real(dp) :: qlatt(3,3),rlatt(3,3)
 real(dp) :: eminmax(2,ebands%nsppol)
 real(dp),allocatable :: jdos(:,:),wmesh(:),cvmw(:),fullbz(:,:)
 real(dp),allocatable :: bdelta(:,:),btheta(:,:)

! *********************************************************************

 ierr = 0
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Find the valence band index for each k and spin ===
 val_idx = get_valence_idx(ebands)

 do spin=1,ebands%nsppol
   if (any(val_idx(:,spin) /= val_idx(1,spin))) then
     write(msg,'(a,i0,a)')&
     'Trying to compute JDOS with a metallic band structure for spin: ',spin,&
     'Assuming you know what you are doing, continuing anyway! '
     MSG_COMMENT(msg)
   end if
 end do

 ! Compute the mean value of the energy spacing.
 ediffs = ebands_edstats(ebands)
 wbroad = broad; if (wbroad <= tol16) wbroad = 0.1 * ediffs%mean
 wstep = step; if (wstep <= tol16) wstep = 0.02 * ediffs%mean

 ! Compute the linear mesh so that it encloses all bands.
 eminmax = get_minmax(ebands, "eig")
 wmax = maxval(eminmax(2,:) - eminmax(1,:))
 nw = nint(wmax/wstep) + 1

 ABI_CALLOC(jdos, (nw, ebands%nsppol))
 ABI_MALLOC(wmesh, (nw))
 wmesh = arth(zero, wstep, nw)

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
           cvmw = ebands%eig(ibc,ik_ibz,spin) - ebands%eig(ibv,ik_ibz,spin) - wmesh
           jdos(:, spin) = jdos(:, spin) + wtk * dirac_delta(cvmw, wbroad)
         end do
       end do

     end do ! ik_ibz
   end do ! spin

   ABI_FREE(cvmw)
   call xmpi_sum(jdos, comm, mpierr)

 case (2)
   ! Tetrahedron method
   ! consistency test: return immediately if cannot use tetra.
   if (ebands%nkpt<2) then
     MSG_WARNING('at least 2 points are needed for tetrahedrons')
     ierr = ierr + 1
   end if
   if (any(ebands%nband /= ebands%nband(1)) ) then
     MSG_WARNING('for tetrahedrons, nband(:) must be constant')
     ierr = ierr + 1
   end if
   if (ebands%nshiftk>1) then
     MSG_WARNING(sjoin("for tetrahedrons, nshiftk must be (0,1) but found: ",itoa(ebands%nshiftk)))
     ierr = ierr + 1
   end if
   if (ierr/=0) return

   ! convert kptrlatt to double and invert.
   rlatt = dble(ebands%kptrlatt)
   call matr3inv(rlatt,qlatt)  ! qlatt refers to the shortest qpt vectors.

   ! Calculate total number of k-points in the full BZ
   nbz = &
      ebands%kptrlatt(1,1)*ebands%kptrlatt(2,2)*ebands%kptrlatt(3,3) &
     +ebands%kptrlatt(1,2)*ebands%kptrlatt(2,3)*ebands%kptrlatt(3,1) &
     +ebands%kptrlatt(1,3)*ebands%kptrlatt(2,1)*ebands%kptrlatt(3,2) &
     -ebands%kptrlatt(1,2)*ebands%kptrlatt(2,1)*ebands%kptrlatt(3,3) &
     -ebands%kptrlatt(1,3)*ebands%kptrlatt(2,2)*ebands%kptrlatt(3,1) &
     -ebands%kptrlatt(1,1)*ebands%kptrlatt(2,3)*ebands%kptrlatt(3,2)
   nbz = nbz * ebands%nshiftk

   ABI_MALLOC(bz2ibz,(nbz))
   ABI_MALLOC(fullbz,(3, nbz))
   nibz = ebands%nkpt

   ! === Make full kpoint grid and get equivalence to irred kpoints ===
   ! * Note: This routines scales badly wrt Kmesh%nbz
   ! TODO should be rewritten, pass kptopt and test whether listkk is faster.
   call get_full_kgrid(bz2ibz,ebands%kptns,fullbz,ebands%kptrlatt,nibz,nbz,ebands%nshiftk,&
     cryst%nsym,ebands%shiftk,cryst%symrel)

   if (ierr == 0) then
     call init_tetra(bz2ibz, cryst%gprimd, qlatt, fullbz, nbz, tetra, ierr, errstr)
     if (ierr/=0) MSG_WARNING(errstr)
   end if
   ABI_FREE(bz2ibz)
   ABI_FREE(fullbz)

   !call tetra_from_kptrlatt(tetra, cryst, ebands%kptopt, ebands%kptrlatt, &
   !  ebands%nshiftk, ebands%shiftk, ebands%nkpt, ebands%kptns)

   if (ierr /= 0) return

   ! For each spin and band, interpolate over kpoints,
   ! calculate integration weights and DOS contribution.
   ABI_MALLOC(cvmw, (nibz))
   ABI_MALLOC(btheta, (nw, nibz))
   ABI_MALLOC(bdelta, (nw, nibz))

   cnt = 0
   do spin=1,ebands%nsppol
     nbv = val_idx(1, spin)
     do ibv=1,nbv
       do ibc=nbv+1,ebands%mband
         cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
         ! For each (c,v) get its contribution
         cvmw = ebands%eig(ibc,:,spin) - ebands%eig(ibv,:,spin)

         ! Calculate integration weights at each irred k-point (Blochl et al PRB 49 16223)
         call tetra_blochl_weights(tetra,cvmw,wmesh(0),wmesh(nw),one,nw,nibz,bcorr0,&
           btheta,bdelta,xmpi_comm_self)

         do ik_ibz=1,nibz
           do iw=1,nw
             jdos(iw,spin) = jdos(iw,spin) + bdelta(iw, ik_ibz)
           end do
         end do

       end do ! ibc
     end do ! ibv
   end do !spin

   call xmpi_sum(jdos, comm, mpierr)

   ! Free memory
   ABI_FREE(btheta)
   ABI_FREE(bdelta)
   ABI_FREE(cvmw)

   call destroy_tetra(tetra)

 case default
   MSG_ERROR(sjoin("Wrong integration method:", itoa(intmeth)))
 end select

 if (ebands%nsppol == 1) jdos = two * jdos

 ! Write data.
 if (my_rank == 0) then
   path = "jdos_gauss.data"; if (intmeth == 2) path = "jdos_tetra.data"
   if (open_file(path, msg, newunit=unt, form="formatted", action="write") /= 0) then
     MSG_ERROR(msg)
   end if
   do iw=1,nw
     write(unt,*)wmesh(iw),(jdos(iw,spin), spin=1,ebands%nsppol)
   end do
   close(unt)
 end if

 ABI_FREE(wmesh)
 ABI_FREE(jdos)

end subroutine ebands_jdos
!!***

!----------------------------------------------------------------------

!!****f* m_ebands/ebands_prtbltztrp
!! NAME
!! ebands_prtbltztrp
!!
!! FUNCTION
!!   output files for BoltzTraP code, which integrates Boltzmann transport quantities
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_prtbltztrp'
!End of the abilint section

 implicit none

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
 character(len=fnlen) :: filename
 character(len=2) :: so_suffix
 character(len=500) :: msg
!arrays
 real(dp) :: nelec(ebands%nsppol)
 character(len=3) :: spinsuffix(ebands%nsppol)

! *************************************************************************

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

   write (iout, '(a)') "GENE                      # Format of input: generic format, with Symmetries"
   write (iout, '(a)') "0 0 0 0.0                 # iskip (not presently used) idebug setgap shiftgap"
   write (iout, '(E15.5,a,F10.4,a)') ebands%fermie*two, " 0.0005 0.4  ", nelec(isppol), &
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
!!  natom = number of atoms in cell.
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
!! SIDE EFFECTS
!!
!! PARENTS
!!      get_tau_k
!!
!! CHILDREN
!!
!! SOURCE

subroutine ebands_prtbltztrp_tau_out (eigen, tempermin, temperinc, ntemper, fermie, fname_radix, kpt, &
&       natom, nband, nelec, nkpt, nspinor, nsppol, nsym, &
&       rprimd, symrel, tau_k)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_prtbltztrp_tau_out'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom, nsym, nband, nkpt, nsppol, nspinor, ntemper
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
&   "  # Fermilevel (Ry), energy grid spacing, energy span around Fermilevel, number of electrons"
   write (iout, '(a)') "CALC                      # CALC (calculate expansion coeff), NOCALC read from file"
   write (iout, '(a)') "3                         # lpfac, number of latt-points per k-point"
   write (iout, '(a)') "BOLTZ                     # run mode (only BOLTZ is supported)"
   write (iout, '(a)') ".15                       # (efcut) energy range of chemical potential"
   write (iout, '(2f8.2,a)')&
&   tempermin+temperinc*dble(itemp),tempermin+temperinc*dble(itemp), "                  # Tmax, temperature grid spacing"
   write (iout, '(2a)') "-1                        # energyrange of bands given ",&
&   "individual DOS output sig_xxx and dos_xxx (xxx is band number)"
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
&   symrel(1,:,isym), &
&   symrel(2,:,isym), &
&   symrel(3,:,isym), &
&   ' ! symmetry rotation matrix isym = ', isym
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

!!****f* m_ebands/ebands_write_xmgrace
!! NAME
!! ebands_write_xmgrace
!!
!! FUNCTION
!!  Write bands in Xmgrace format. This routine should be called by a single processor.
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

subroutine ebands_write_xmgrace(ebands, path, kptbounds)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ebands_write_xmgrace'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ebands_t),intent(in) :: ebands
 character(len=*),intent(in) :: path
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
   ! Find correspondence between kptbounds and k-points in ebands.
   nkbounds = size(kptbounds, dim=2)
   ABI_MALLOC(bounds2kpt, (nkbounds))
   bounds2kpt = 1
   do ii=1,nkbounds
      start = 1
      do ik=start,ebands%nkpt
        if (isamek(ebands%kptns(:, ik), kptbounds(:, ii), g0)) then
          bounds2kpt(ii) = ik; start = ik + 1; exit
        end if
      end do
   end do
 end if

 if (open_file(path, msg, newunit=unt, form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write(unt,'(a)') "# Grace project file"
 write(unt,'(a)') "# Generated by Abinit"
 write(unt,'(4(a,i0))') &
   "# mband: ",ebands%mband,", nkpt: ",ebands%nkpt,", nsppol: ",ebands%nsppol,", nspinor: ",ebands%nspinor
 write(unt,'(a,f8.2,a,i0,2(a,f8.2))') &
   "# nelect: ",ebands%nelect,", occopt: ",ebands%occopt,", tsmear: ",ebands%tsmear,", tphysel: ",ebands%tphysel
 write(unt,'(a,f8.2,a)') "# Energies are in eV. Zero set to efermi: ",ebands%fermie, " [Ha]"
 write(unt,'(a)')"# List of k-points and their index (C notation i.e. count from 0)"
 do ik=1,ebands%nkpt
   write(unt, "(a)")sjoin("#", itoa(ik-1), ktoa(ebands%kptns(:,ik)))
 end do
 !write(unt,'(a)') "@version 50113"
 write(unt,'(a)') "@page size 792, 612"
 write(unt,'(a)') "@page scroll 5%"
 write(unt,'(a)') "@page inout 5%"
 write(unt,'(a)') "@link page off"
 write(unt,'(a)') "@with g0"
 write(unt,'(a)') "@world xmin 0.00"
 write(unt,'(a,i0)') '@world xmax ',ebands%nkpt
 write(unt,'(a,e16.8)') '@world ymin ',minval((ebands%eig - ebands%fermie) * Ha_eV)
 write(unt,'(a,e16.8)') '@world ymax ',maxval((ebands%eig - ebands%fermie) * Ha_eV)
 write(unt,'(a)') '@default linewidth 1.5'
 write(unt,'(a)') '@xaxis  tick on'
 write(unt,'(a)') '@xaxis  tick major 1'
 write(unt,'(a)') '@xaxis  tick major color 1'
 write(unt,'(a)') '@xaxis  tick major linestyle 3'
 write(unt,'(a)') '@xaxis  tick major grid on'
 write(unt,'(a)') '@xaxis  tick spec type both'
 write(unt,'(a)') '@xaxis  tick major 0, 0'
 ! TODO: Test whether this is ok.
 if (nkbounds /= 0) then
   write(unt,'(a,i0)') '@xaxis  tick spec ',nkbounds
   do ik=1,nkbounds
      write(unt,'(a,i0,a,a)') '@xaxis  ticklabel ',ik-1,',', "foo"
      write(unt,'(a,i0,a,i0)') '@xaxis  tick major ',ik-1,' , ',bounds2kpt(ik)
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
        write(unt,'(i0,1x,e16.8)') ik-1, (ebands%eig(band, ik, spin) - ebands%fermie) * Ha_eV
     end do
     write(unt,'(a)') '&'
   end do
 end do

 close(unt)

 if (allocated(bounds2kpt)) then
   ABI_FREE(bounds2kpt)
 end if

end subroutine ebands_write_xmgrace
!!***

end module m_ebands
!!***
