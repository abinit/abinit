!!****m* ABINIT/m_wvl_descr_psp
!! NAME
!!  wvl_descr_psp
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (DC)
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

module m_wvl_descr_psp

 use defs_basis
 use m_errors
 use m_abicore

 use defs_datatypes, only : pseudopotential_type, pseudopotential_gth_type

 implicit none

 private
!!***

 public :: wvl_descr_psp_set
 public :: wvl_descr_psp_fill
 public :: wvl_descr_free
 public :: wvl_descr_atoms_set
 public :: wvl_descr_atoms_set_sym
!!***

contains
!!***

!!****f* defs_wvltypes/wvl_descr_psp_set
!! NAME
!! wvl_descr_psp_set
!!
!! FUNCTION
!! Defines the part of the wvl%atoms%-datastructure which depends on psps
!!
!! INPUTS
!! psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! dtset <type(dataset_type)>=all input variables for this dataset
!!
!! OUTPUT
!! nsppol=number of spin components
!! spinat(3,natom)=spin polarization on each atom
!! wvl <type(wvl_internal_type)> = wavelet type
!!                 | psppar   = The covalence radii for each pseudo
!!                 | pspcod   = the format -or code- of psp generation
!!                 | iasctype = semicore code (see defs_datatype)
!!                 | nzatom   = charge of the nucleus
!!                 | nelpsp   = the ionic pseudo-charge
!!                 | natsc    = number of atoms with semicore
!!
!! PARENTS
!!      m_gstate
!!
!! CHILDREN
!!      astruct_set_symmetries,symmetry_set_n_sym,wrtout
!!
!! SOURCE

subroutine wvl_descr_psp_set(filoccup, nsppol, psps, spinat, wvl)

 use defs_wvltypes
#if defined HAVE_BIGDFT
 use BigDFT_API, only: aoig_set,UNINITIALIZED,dict_init,dict_free,dictionary, &
&                      operator(//), bigdft_mpi, dict_set
 use BigDFT_API, only: psp_data_merge_to_dict, psp_dict_fill_all, atomic_info, &
&                psp_dict_analyse, atomic_data_set_from_dict, merge_input_file_to_dict
#endif
  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in)                    :: nsppol
  type(wvl_internal_type), intent(inout) :: wvl
  type(pseudopotential_type), intent(in) :: psps
  character(len = *), intent(in) :: filoccup
!arrays
  real(dp),intent(in) :: spinat(:,:)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
  integer :: ityp,pspcod
  logical :: exists
  real(dp) :: radii_cf(3)
  character(len=2) :: symbol
  character(len=27) :: filename
  type(dictionary), pointer :: dict
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

!We create the atoms_data structure, the part that is dependent from psp.
 do ityp=1,size(psps%pspcod)
   wvl%atoms%psppar(:,:,ityp)= psps%gth_params%psppar(:,:,ityp)
   wvl%atoms%npspcode(ityp)  = merge(psps%pspcod(ityp),7,psps%pspcod(ityp)/=17)
   wvl%atoms%ixcpsp(ityp)    = psps%pspxc(ityp)
   wvl%atoms%nzatom(ityp)    = int(psps%znucltypat(ityp))
   wvl%atoms%nelpsp(ityp)    = int(psps%ziontypat(ityp))
 end do
 wvl%atoms%donlcc = .false.

 call dict_init(dict)
 radii_cf = UNINITIALIZED(1._dp)
 do ityp = 1, wvl%atoms%astruct%ntypes, 1
   call atomic_info(wvl%atoms%nzatom(ityp), wvl%atoms%nelpsp(ityp), &
&   symbol = symbol)
   write(wvl%atoms%astruct%atomnames(ityp), "(A)") symbol
   filename = 'psppar.' // trim(symbol)
   pspcod=psps%pspcod(ityp);if (pspcod==17) pspcod=7 ! PAW specificity
   call psp_data_merge_to_dict(dict // filename, int(psps%znucltypat(ityp)), &
&   int(psps%ziontypat(ityp)), pspcod, psps%pspxc(ityp), &
&   psps%gth_params%psppar(0:4,0:6,ityp), radii_cf, &
&   UNINITIALIZED(1._dp), UNINITIALIZED(1._dp))
   call psp_dict_fill_all(dict, trim(symbol), psps%pspxc(ityp))
 end do

 call psp_dict_analyse(dict, wvl%atoms)
 inquire(file = filoccup, exist = exists)
 if (exists) then
   call merge_input_file_to_dict(dict//"Atomic occupation",filoccup,bigdft_mpi)
 end if
 ! Need to add spinat in the dictionary, later
 if (spinat(1,1) == zero .and. .false.) then
   call dict_set(dict // "spinat", spinat(1,1))
 end if
 call atomic_data_set_from_dict(dict, "Atomic occupation", wvl%atoms, nsppol)
 call dict_free(dict)

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) nsppol,wvl%h(1),psps%npsp,filoccup,spinat(1,1)
#endif

end subroutine wvl_descr_psp_set
!!***

!!****f* defs_wvltypes/wvl_descr_psp_fill
!! NAME
!! wvl_descr_psp_fill
!!
!! FUNCTION
!! Read the radii from @pspunit, if any and fill values with default otherwise.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_psp_hgh,m_pspini
!!
!! CHILDREN
!!      astruct_set_symmetries,symmetry_set_n_sym,wrtout
!!
!! SOURCE

subroutine wvl_descr_psp_fill(gth_params, ipsp, ixc, nelpsp, nzatom, pspunit)

#if defined HAVE_BIGDFT
  use BigDFT_API, only: atomic_info, UNINITIALIZED, psp_from_data
#endif
  implicit none

!Arguments ------------------------------------
  integer, intent(in) :: ipsp, pspunit, nzatom, nelpsp, ixc
  type(pseudopotential_gth_type), intent(inout) :: gth_params
!Local variables-------------------------------
#if defined HAVE_BIGDFT
  integer :: ios, ii, nzatom_, nelpsp_, npspcode_, ixc_
  real(dp) :: ehomo, radfine
  logical :: exists
  character(len = 2) :: symbol
  character(len=100) :: line
  character(len=500) :: message
#endif

! ***************************************************************************

#if defined HAVE_BIGDFT

  ! check if gth_params%psppar have been set
 if (any(gth_params%psppar == UNINITIALIZED(1._dp))) then
   call atomic_info(nzatom, nelpsp, symbol = symbol)
   ixc_ = ixc
   if (ixc_>0.and.ixc_< 10) ixc_=1
   if (ixc_>0.and.ixc_>=10) ixc_=11
   call psp_from_data(symbol, nzatom_, nelpsp_, npspcode_, ixc_, &
&   gth_params%psppar(:,:,ipsp), exists)
   if(.not. exists) then
     write(message,'(a,a,a,a)')ch10,&
&     "wvl_descr_psp_fill : bug, chemical element not found in BigDFT table",ch10,&
&     "Action: upgrade BigDFT table"
     call wrtout(ab_out,message,'COLL')
     ABI_BUG(message)
   end if
   gth_params%set(ipsp) = .true.
 end if

  ! Try to read radii from pspunit
 if (pspunit /= 0) then
   read (pspunit, '(a100)', iostat = ios) line
   if (ios /= 0) then
     line=''
   end if
     !  We try to read the values from the pseudo.
   read (line, *, iostat = ios) gth_params%radii_cf(ipsp, 1), gth_params%radii_cf(ipsp, 2), &
&   gth_params%radii_cf(ipsp, 3)
   if (ios /= 0 .or. gth_params%radii_cf(ipsp, 3) < zero) then
     read (line, *, iostat = ios) gth_params%radii_cf(ipsp, 1), gth_params%radii_cf(ipsp, 2)
     gth_params%radii_cf(ipsp, 3) = UNINITIALIZED(gth_params%radii_cf(ipsp, 3))
   end if
   if (ios /= 0 .or. gth_params%radii_cf(ipsp, 1) < zero .or. gth_params%radii_cf(ipsp, 2) < zero) then
     gth_params%radii_cf(ipsp, 1) = UNINITIALIZED(gth_params%radii_cf(ipsp, 1))
     gth_params%radii_cf(ipsp, 2) = UNINITIALIZED(gth_params%radii_cf(ipsp, 2))
     write(message, '(a,a,a,a,a,a,a)' ) '-', ch10,&
&     '- wvl_descr_psp_fill : COMMENT -',ch10,&
&     "-  the pseudo-potential does not include geometric information,",ch10,&
&     '-  values have been computed.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,  message,'COLL')
   end if
 end if

  ! Update radii.
 if (gth_params%radii_cf(ipsp, 1) == UNINITIALIZED(gth_params%radii_cf(ipsp, 1))) then
   call atomic_info(nzatom, nelpsp, ehomo = ehomo)
     !assigning the radii by calculating physical parameters
   gth_params%radii_cf(ipsp, 1)=1._dp/sqrt(abs(2._dp*ehomo))
 end if
 if (gth_params%radii_cf(ipsp, 2) == UNINITIALIZED(gth_params%radii_cf(ipsp, 2))) then
   radfine = gth_params%psppar(0, 0, ipsp)
   do ii = 1, 4, 1
     if (gth_params%psppar(ii, 0, ipsp) /= zero) then
       radfine = min(radfine, gth_params%psppar(ii, 0, ipsp))
     end if
   end do
   gth_params%radii_cf(ipsp, 2)=radfine
 end if
 if (gth_params%radii_cf(ipsp, 3) == UNINITIALIZED(gth_params%radii_cf(ipsp, 3))) then
   gth_params%radii_cf(ipsp, 3)=gth_params%radii_cf(ipsp, 2)
 end if

  ! Set flag.
 if (gth_params%radii_cf(ipsp, 1) >= 0.d0 .and. gth_params%radii_cf(ipsp, 2) >= 0.d0) then
   gth_params%hasGeometry(ipsp) = .true.
 else
   gth_params%hasGeometry(ipsp) = .false.
 end if

 write(message, '(a,f12.7,a,f12.7,a,f12.7)' )&
& '  radii_cf(1)=', gth_params%radii_cf(ipsp, 1),&
& '; radii_cf(2)=', gth_params%radii_cf(ipsp, 2),&
& '; rad_cov=', gth_params%radii_cf(ipsp, 3)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

! Some consistency checks on radii, to be moved earlier, but need to update test refs.
! gth_params%radii_cf(ipsp, 3) = min(gth_params%radii_cf(ipsp, 3), gth_params%radii_cf(ipsp, 2))

! This was one before
! maxrad=zero
! do ii=0,2,1
!   if (ii==1) maxrad=zero
!   if (gth_params%psppar(ii,0,ipsp)/=zero) maxrad=max(maxrad,gth_params%psppar(ii,0,ipsp))
! end do
! if (maxrad== zero) then
!   gth_params%radii_cf(ipsp,3)=zero
! else
!   gth_params%radii_cf(ipsp,3)=max( &
!&    min(dtset%wvl_crmult*psps%gth_params%radii_cf(ipsp,1), &
!&    15._dp*maxrad)/dtset%wvl_frmult,gth_params%radii_cf(ipsp,2))
! end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) ipsp,pspunit,nzatom,nelpsp,ixc,gth_params%psppar
#endif

end subroutine wvl_descr_psp_fill
!!***

!!****f* defs_wvltypes/wvl_descr_free
!!
!! NAME
!! wvl_descr_free
!!
!! FUNCTION
!! Free the wvl%atoms% datastructure (deallocate or nullify)
!!
!! INPUTS
!! wvl <type(wvl_internal_type)>=internal variables for wavelets
!!
!! OUTPUT
!! wvl <type(wvl_internal_type)>=internal variables for wavelets
!!
!! PARENTS
!!      m_gstate,m_memeval
!!
!! CHILDREN
!!      astruct_set_symmetries,symmetry_set_n_sym,wrtout
!!
!! SOURCE

subroutine wvl_descr_free(wvl)

 use defs_wvltypes
#if defined HAVE_BIGDFT
 use BigDFT_API, only : deallocate_atoms_data
 use dynamic_memory
#endif
  implicit none

!Arguments ------------------------------------
!scalars
  type(wvl_internal_type), intent(inout) :: wvl
!arrays

!Local variables-------------------------------
!scalars

! *********************************************************************

#if defined HAVE_BIGDFT
!These arrays are pointers on memory handled by ABINIT.
 nullify(wvl%atoms%astruct%sym%irrzon)
 nullify(wvl%atoms%astruct%sym%phnons)
 if (associated(wvl%atoms%nlccpar)) then
   call f_free_ptr(wvl%atoms%nlccpar)
 end if
 call deallocate_atoms_data(wvl%atoms)
#endif
 if(allocated(wvl%npspcode_paw_init_guess)) then
   ABI_DEALLOCATE(wvl%npspcode_paw_init_guess)
 end if
end subroutine wvl_descr_free
!!***

!!****f* defs_wvltypes/wvl_descr_atoms_set
!!
!! NAME
!! wvl_descr_atoms_set
!!
!! FUNCTION
!! Defines wvl%atoms% data structure
!!
!! INPUTS
!! acell(3)=unit cell length scales (bohr)
!! dtset <type(dataset_type)>=all input variables for this dataset
!!
!! OUTPUT
!! wvl <type(wvl_internal_type)>= wavelet type
!!                 | nat      =  number of atoms
!!                 | ntypes   =  number of species
!!                 | alat1    =  acell(1)
!!                 | alat2    =  acell(2)
!!                 | alat3    =  acell(3)
!!                 | iatype   =  types for atoms
!!                 | lfrztyp  =  flag for the movement of atoms.
!!                 | natpol   =  integer related to polarisation at the first step
!!
!! PARENTS
!!      m_gstate,m_memeval
!!
!! CHILDREN
!!      astruct_set_symmetries,symmetry_set_n_sym,wrtout
!!
!! SOURCE

subroutine wvl_descr_atoms_set(acell, icoulomb, natom, ntypat, typat, wvl)

 use defs_wvltypes
#if defined HAVE_BIGDFT
 use BigDFT_API, only: atoms_data_null,f_routine,f_release_routine,&
&                      astruct_set_n_atoms,astruct_set_n_types,&
&                      allocate_atoms_nat,allocate_atoms_ntypes
#endif
  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in)                    :: icoulomb, natom, ntypat
  type(wvl_internal_type), intent(inout) :: wvl
  !arrays
  integer, intent(in)                    :: typat(natom)
  real(dp), intent(in)                   :: acell(3)

!Local variables-------------------------------
!scalars
#if defined HAVE_BIGDFT
  integer :: itype
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 call f_routine('wvl_descr_atoms_set')

 wvl%atoms=atoms_data_null()

!We create the atoms_data structure from this dataset
!to be used later in BigDFT routines.
 if (icoulomb == 0) then
   wvl%atoms%astruct%geocode = 'P'
 else if (icoulomb == 1) then
   wvl%atoms%astruct%geocode = 'F'
 else if (icoulomb == 2) then
   wvl%atoms%astruct%geocode = 'S'
 end if
 write(wvl%atoms%astruct%units, "(A)") "Bohr"

 call astruct_set_n_atoms(wvl%atoms%astruct, natom)
 call astruct_set_n_types(wvl%atoms%astruct, ntypat)

 do itype = 1, ntypat, 1
   write(wvl%atoms%astruct%atomnames(itype), "(A,I2)") "At. type", itype
 end do
 wvl%atoms%astruct%cell_dim(1)   =  acell(1)
 wvl%atoms%astruct%cell_dim(2)   =  acell(2)
 wvl%atoms%astruct%cell_dim(3)   =  acell(3)
 wvl%atoms%astruct%iatype   = typat

 wvl%atoms%astruct%sym%symObj = 0

 call allocate_atoms_nat(wvl%atoms)
 call allocate_atoms_ntypes(wvl%atoms)

 call f_release_routine()

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) icoulomb,natom,ntypat,wvl%h(1),typat(1),acell(1)
#endif

end subroutine wvl_descr_atoms_set
!!***

!!****f* defs_wvltypes/wvl_descr_atoms_set_sym
!!
!! NAME
!! wvl_descr_atoms_set_sym
!!
!! FUNCTION
!! Add symmetry information to  wvl%atoms data structure.
!!
!! INPUTS
!! acell(3)=unit cell length scales (bohr)
!! dtset <type(dataset_type)>=all input variables for this dataset
!!
!! OUTPUT
!! wvl <type(wvl_internal_type)>= wavelet type
!!                 | nat      =  number of atoms
!!                 | ntypes   =  number of species
!!                 | alat1    =  acell(1)
!!                 | alat2    =  acell(2)
!!                 | alat3    =  acell(3)
!!                 | iatype   =  types for atoms
!!                 | lfrztyp  =  flag for the movement of atoms.
!!                 | natpol   =  integer related to polarisation at the first step
!!
!! PARENTS
!!      m_gstate
!!
!! CHILDREN
!!      astruct_set_symmetries,symmetry_set_n_sym,wrtout
!!
!! SOURCE

subroutine wvl_descr_atoms_set_sym(wvl, efield, irrzon, nsppol, nsym, phnons, &
           & symafm, symrel, tnons, tolsym)

 use defs_wvltypes
 use m_ab7_symmetry
#if defined HAVE_BIGDFT
 use BigDFT_API, only: astruct_set_symmetries
#endif
  implicit none

!Arguments ------------------------------------
!scalars
  integer, intent(in)                    :: nsppol,nsym
  real(dp), intent(in)                   :: tolsym
  type(wvl_internal_type), intent(inout) :: wvl
  !arrays
  integer, intent(in)                    :: symafm(3,3,nsym), symrel(3,3,nsym)
  integer, target, intent(in)            :: irrzon(:,:,:)
  real(dp), target, intent(in)           :: phnons(:,:,:)
  real(dp), intent(in)                   :: efield(3), tnons(3,nsym)

!Local variables-------------------------------
!scalars
#if defined HAVE_BIGDFT
  integer :: errno
 character(len=500) :: message
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 write(message, '(a,a)' ) ch10,&
& ' wvl_descr_atoms_set_sym: Create symmetries for the wvl%atoms object.'
 call wrtout(std_out,message,'COLL')

 wvl%atoms%astruct%sym%symObj = -1
 nullify(wvl%atoms%astruct%sym%irrzon)
 nullify(wvl%atoms%astruct%sym%phnons)
 call astruct_set_symmetries(wvl%atoms%astruct, (nsym <= 1), tolsym, efield, nsppol)
 if (nsym > 1) then
   call symmetry_set_n_sym(wvl%atoms%astruct%sym%symObj, nsym, symrel, tnons, symafm, errno)
 end if
 wvl%atoms%astruct%sym%irrzon => irrzon
 wvl%atoms%astruct%sym%phnons => phnons

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) nsppol,nsym,tolsym,wvl%h(1),symafm(1,1,1),symrel(1,1,1),&
& irrzon(1,1,1),phnons(1,1,1),efield(1),tnons(1,1)
#endif

end subroutine wvl_descr_atoms_set_sym
!!***

end module m_wvl_descr_psp
!!***
