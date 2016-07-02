!!****f* defs_wvltypes/wvl_descr_psp_set
!!
!! NAME
!! wvl_descr_psp_set
!!
!! FUNCTION
!! Defines the part of the wvl%atoms%-datastructure which
!! depends on psps
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
!!      gstate
!!
!! CHILDREN
!!      atomic_info,psp_from_data,wrtout
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_descr_psp_set(filoccup, nsppol, psps, spinat, wvl)

 use m_profiling_abi
 use m_errors

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
#if defined HAVE_DFT_BIGDFT
 use BigDFT_API, only: aoig_set,UNINITIALIZED,dict_init,dict_free,dictionary, &
&                      operator(//), bigdft_mpi, dict_set
 use BigDFT_API, only: psp_data_merge_to_dict, psp_dict_fill_all, atomic_info, &
&                psp_dict_analyse, atomic_data_set_from_dict, merge_input_file_to_dict
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_descr_psp_set'
!End of the abilint section

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
#if defined HAVE_DFT_BIGDFT
  integer :: ityp,pspcod
  logical :: exists
  real(dp) :: radii_cf(3)
  character(len=2) :: symbol
  character(len=27) :: filename
  type(dictionary), pointer :: dict
#endif

! *********************************************************************

#if defined HAVE_DFT_BIGDFT

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
!!
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
!!      psp10in,psp2in,psp3in,pspatm
!!
!! CHILDREN
!!      atomic_info,psp_from_data,wrtout
!!
!! SOURCE
subroutine wvl_descr_psp_fill(gth_params, ipsp, ixc, nelpsp, nzatom, pspunit)

  use defs_datatypes
  use m_errors
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only: atomic_info, UNINITIALIZED, psp_from_data
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_descr_psp_fill'
 use interfaces_14_hidewrite
!End of the abilint section

  implicit none

!Arguments ------------------------------------
  integer, intent(in) :: ipsp, pspunit, nzatom, nelpsp, ixc
  type(pseudopotential_gth_type), intent(inout) :: gth_params
!Local variables-------------------------------
#if defined HAVE_DFT_BIGDFT
  integer :: ios, ii, nzatom_, nelpsp_, npspcode_, ixc_
  real(dp) :: ehomo, radfine
  logical :: exists
  character(len = 2) :: symbol
  character(len=100) :: line
  character(len=500) :: message
#endif

! ***************************************************************************

#if defined HAVE_DFT_BIGDFT

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
     MSG_BUG(message)
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
&     "-  the pseudo-potential does not include geometric informations,",ch10,&
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
