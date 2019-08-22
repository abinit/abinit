!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_wvl_wfs
!! NAME
!! m_wvl_wfs
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DC)
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

module m_wvl_wfs

 use defs_basis
 use m_errors
 use m_abicore

 use defs_datatypes, only : pseudopotential_type

 implicit none

 private
!!***

 public :: wvl_wfs_set
 public :: derf_ab
 public :: wvl_wfs_free
 public :: wvl_wfs_lr_copy
!!***

contains
!!***

!!****f* ABINIT/wvl_wfs_set
!! NAME
!! wvl_wfs_set
!!
!! FUNCTION
!! Compute the access keys for the wavefunctions when the positions
!! of the atoms are given.
!!
!! For memory occupation optimisation reasons, the wavefunctions are not allocated
!! here. See the initialisation routines wvl_wfsinp_disk(), wvl_wfsinp_scratch()
!! and wvl_wfsinp_reformat() to do it. After allocation, use wvl_wfs_free()
!! to deallocate all stuff (descriptors and arrays).
!!
!! INPUTS
!!  dtset <type(dataset_type)>=internal variables used by wavelets, describing
!!   | wvl_internal=desciption of the wavelet box.
!!   | natom=number of atoms.
!!  mpi_enreg=informations about MPI parallelization
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!   | keys=its access keys for compact storage.
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_wfs_set(alphadiis, spinmagntarget, kpt, me, natom, nband, nkpt, nproc, nspinor, &
&  nsppol, nwfshist, occ, psps, rprimd, wfs, wtk, wvl, wvl_crmult, wvl_frmult, xred)

 use defs_wvltypes

 use m_geometry, only : xred2xcart
#if defined HAVE_BIGDFT
 use BigDFT_API, only: createWavefunctionsDescriptors, orbitals_descriptors, &
       & orbitals_communicators, allocate_diis_objects, &
       & input_variables, check_linear_and_create_Lzd, check_communications, &
       & INPUT_IG_OFF, nullify_locreg_descriptors
#endif
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom, nkpt, nsppol, nspinor, nband, nwfshist,me,nproc
 real(dp), intent(in) :: spinmagntarget, wvl_crmult, wvl_frmult, alphadiis
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_wf_type),intent(out) :: wfs
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 real(dp), intent(in) :: kpt(3,nkpt)
 real(dp), intent(in) :: wtk(nkpt), occ(:)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
!scalars
 integer :: ii,idata, norb, norbu, norbd
 logical :: parallel
 character(len=500) :: message
 type(input_variables) :: in  ! To be removed, waiting for BigDFT upgrade
!arrays
 real(dp), allocatable :: kpt_(:,:)
 real(dp),allocatable :: xcart(:,:)
#endif

! *********************************************************************

#if defined HAVE_BIGDFT

 parallel = (nproc > 1)

!Consistency checks, are all pseudo true GTH pseudo with geometric informations?
!Skip for PAW case: we do not have GTH parameters
 do idata = 1, psps%npsp, 1
   if (.not. psps%gth_params%set(idata) .and. psps%usepaw==0) then
     write(message, '(a,a,a,a,I0,a,a,a)' ) ch10,&
&     ' wvl_wfs_set:  consistency checks failed,', ch10, &
&     '  no GTH parameters found for type number ', idata, '.', ch10, &
&     '  Check your input pseudo files.'
     MSG_ERROR(message)
   end if
   if (.not. psps%gth_params%hasGeometry(idata)) then
     write(message, '(a,a,a,a,a,a)' ) ch10,&
&     ' wvl_wfs_set:  consistency checks failed,', ch10, &
&     '  the given GTH parameters has no geometry informations.', ch10, &
&     '  Upgrade your input pseudo files to GTH with geometric informations.'
     MSG_ERROR(message)
   end if
 end do

!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, natom))
 call xred2xcart(natom, rprimd, xcart, xred)

!Nullify possibly unset pointers
 nullify(wfs%ks%psi)
 nullify(wfs%ks%hpsi)
 nullify(wfs%ks%psit)

!Static allocations.
 norb = nband / nkpt
 norbu = 0
 norbd = 0
 if (nsppol == 2) then
   if (spinmagntarget < -real(90, dp)) then
     norbu = min(norb / 2, norb)
   else
     norbu = min(norb / 2 + int(spinmagntarget), norb)
   end if
   norbd = norb - norbu
 else
   norbu = norb
   norbd = 0
 end if
 ABI_ALLOCATE(kpt_, (3, nkpt))
 do ii = 1, nkpt
   kpt_(:,ii) = kpt(:,ii) / (/ rprimd(1,1), rprimd(2,2), rprimd(3,3) /) * two_pi
 end do

 call orbitals_descriptors(me, nproc,norb,norbu,norbd,nsppol,nspinor, &
& nkpt,kpt_,wtk,wfs%ks%orbs,.false.)
 ABI_DEALLOCATE(kpt_)
!We copy occ_orig to wfs%ks%orbs%occup
 wfs%ks%orbs%occup(1:norb * nkpt) = occ(1:norb * nkpt)
!We allocate the eigen values storage.
 ABI_ALLOCATE(wfs%ks%orbs%eval,(wfs%ks%orbs%norb * wfs%ks%orbs%nkpts))

 write(message, '(a,a)' ) ch10,&
& ' wvl_wfs_set: Create access keys for wavefunctions.'
 call wrtout(std_out,message,'COLL')

 call nullify_locreg_descriptors(wfs%ks%lzd%Glr)
 wfs%ks%lzd%Glr = wvl%Glr
 call createWavefunctionsDescriptors(me, wvl%h(1), wvl%h(2), wvl%h(3), &
& wvl%atoms, xcart, psps%gth_params%radii_cf, &
& wvl_crmult, wvl_frmult, wfs%ks%lzd%Glr)
!The memory is not allocated there for memory occupation optimisation reasons.

 call orbitals_communicators(me,nproc,wfs%ks%lzd%Glr,wfs%ks%orbs,wfs%ks%comms)

 write(message, '(a,2I8)' ) &
& '  | all orbitals have coarse segments, elements:', &
& wfs%ks%lzd%Glr%wfd%nseg_c, wfs%ks%lzd%Glr%wfd%nvctr_c
 call wrtout(std_out,message,'COLL')
 write(message, '(a,2I8)' ) &
& '  | all orbitals have fine   segments, elements:', &
& wfs%ks%lzd%Glr%wfd%nseg_f, 7 * wfs%ks%lzd%Glr%wfd%nvctr_f
 call wrtout(std_out,message,'COLL')

!allocate arrays necessary for DIIS convergence acceleration
 call allocate_diis_objects(nwfshist,alphadiis,&
& sum(wfs%ks%comms%ncntt(0:nproc-1)), wfs%ks%orbs%nkptsp, wfs%ks%orbs%nspinor, &
& wfs%ks%diis)

 ABI_DATATYPE_ALLOCATE(wfs%ks%confdatarr, (wfs%ks%orbs%norbp))
 call default_confinement_data(wfs%ks%confdatarr,wfs%ks%orbs%norbp)

 call check_linear_and_create_Lzd(me,nproc,INPUT_IG_OFF,wfs%ks%lzd,&
& wvl%atoms,wfs%ks%orbs,nsppol,xcart)
 wfs%ks%lzd%hgrids = wvl%h

!check the communication distribution
 call check_communications(me,nproc,wfs%ks%orbs,wfs%ks%Lzd,wfs%ks%comms)

!Deallocations
 ABI_DEALLOCATE(xcart)

!DEBUG
 write(std_out,*) 'wvl_wfs_set: TODO, update BigDFT sic_input_variables_default()'
!ENDDEBUG
 call sic_input_variables_default(in)

 wfs%ks%SIC                  = in%SIC
 wfs%ks%exctxpar             = "OP2P"
 wfs%ks%c_obj                = 0
 wfs%ks%orthpar%directDiag   = .true.
 wfs%ks%orthpar%norbpInguess = 5
 wfs%ks%orthpar%bsLow        = 300
 wfs%ks%orthpar%bsUp         = 800
 wfs%ks%orthpar%methOrtho    = 0
 wfs%ks%orthpar%iguessTol    = 1.d-4

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) natom,nkpt,nsppol,nspinor,nband,nwfshist,me,nproc,&
& spinmagntarget,wvl_crmult,wvl_frmult,alphadiis,psps%npsp,wfs%ks,wvl%h(1),&
& kpt(1,1),wtk(1),occ(1),rprimd(1,1),xred(1,1)
#endif

end subroutine wvl_wfs_set
!!***

!!****f* ABINIT/derfcf
!! NAME
!! derfcf
!!
!! FUNCTION
!! Some wrappers for BigDFT which uses different names for the same routines.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine derfcf(derfc_yy,yy)

 use m_special_funcs,  only : abi_derfc
 implicit none
!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: yy
 real(dp),intent(out) :: derfc_yy
!Local variables-------------------------------

! *********************************************************************

 derfc_yy = abi_derfc(yy)

end subroutine derfcf
!!***

!!****f* ABINIT/derf_ab
!! NAME
!! derf_ab
!!
!! FUNCTION
!! Some wrappers for BigDFT which uses different names for the same routines.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      mklocl_realspace,mklocl_wavelets
!!
!! CHILDREN
!!
!! SOURCE

subroutine derf_ab(derf_yy,yy)

 use m_special_funcs,  only : abi_derf
 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: yy
 real(dp),intent(out) :: derf_yy

!Local variables-------------------------------

! *********************************************************************

 derf_yy = abi_derf(yy)

end subroutine derf_ab
!!***


!!****f* ABINIT/wvl_wfs_free
!!
!! NAME
!! wvl_wfs_free
!!
!! FUNCTION
!! Freeing routine.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  wfs <type(wvl_wf_type)>=wavefunctions informations in a wavelet basis.
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      deallocate_comms,deallocate_lr,deallocate_lzd_except_glr
!!      deallocate_orbs,f_free_ptr
!!
!! SOURCE


subroutine wvl_wfs_free(wfs)

 use defs_wvltypes
#if defined HAVE_BIGDFT
 use BigDFT_API, only: deallocate_Lzd_except_Glr, deallocate_lr, &
      & deallocate_orbs, deallocate_comms
 use dynamic_memory
#endif
 implicit none

!Arguments ------------------------------------
!scalars
 type(wvl_wf_type),intent(inout) :: wfs

!Local variables -------------------------

! *********************************************************************

#if defined HAVE_BIGDFT
 call deallocate_Lzd_except_Glr(wfs%ks%lzd)
 call deallocate_lr(wfs%ks%lzd%Glr)
 call deallocate_orbs(wfs%ks%orbs)
 call deallocate_comms(wfs%ks%comms)
 if (associated(wfs%ks%orbs%eval))  then
   ABI_DEALLOCATE(wfs%ks%orbs%eval)
 end if
 ABI_DATATYPE_DEALLOCATE(wfs%ks%confdatarr)

 if (associated(wfs%ks%psi)) then
   call f_free_ptr(wfs%ks%psi)
 end if
 if (associated(wfs%ks%hpsi)) then
   call f_free_ptr(wfs%ks%hpsi)
 end if
 if (associated(wfs%ks%psit)) then
   call f_free_ptr(wfs%ks%psit)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) wfs%ks
#endif

end subroutine wvl_wfs_free
!!***


!!****f* defs_wvltypes/wvl_wfs_lr_copy
!!
!! NAME
!! wvl_wfs_lr_copy
!!
!! FUNCTION
!! Copy the wvl%Glr datastructure geometry part to wfs%Glr.
!!
!! INPUTS
!! wvl <type(wvl_internal_type)> = input localisation region
!!
!! OUTPUT
!! wfs <type(wvl_wf_type)> = output localistaion region
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_wfs_lr_copy(wfs, wvl)

 use defs_wvltypes
 implicit none

!Arguments ------------------------------------
!scalars
  type(wvl_internal_type), intent(in)  :: wvl
  type(wvl_wf_type), intent(inout)     :: wfs
!arrays

!Local variables-------------------------------

! *********************************************************************

#if defined HAVE_BIGDFT
!Use global localization region for the moment.
 wfs%ks%lzd%Glr%geocode    = wvl%Glr%geocode
 wfs%ks%lzd%Glr%hybrid_on  = wvl%Glr%hybrid_on
 wfs%ks%lzd%Glr%ns1        = wvl%Glr%ns1
 wfs%ks%lzd%Glr%ns2        = wvl%Glr%ns2
 wfs%ks%lzd%Glr%ns3        = wvl%Glr%ns3
 wfs%ks%lzd%Glr%nsi1       = wvl%Glr%nsi1
 wfs%ks%lzd%Glr%nsi2       = wvl%Glr%nsi2
 wfs%ks%lzd%Glr%nsi3       = wvl%Glr%nsi3
 wfs%ks%lzd%Glr%d%n1       = wvl%Glr%d%n1
 wfs%ks%lzd%Glr%d%n2       = wvl%Glr%d%n2
 wfs%ks%lzd%Glr%d%n3       = wvl%Glr%d%n3
 wfs%ks%lzd%Glr%d%nfl1     = wvl%Glr%d%nfl1
 wfs%ks%lzd%Glr%d%nfu1     = wvl%Glr%d%nfu1
 wfs%ks%lzd%Glr%d%nfl2     = wvl%Glr%d%nfl2
 wfs%ks%lzd%Glr%d%nfu2     = wvl%Glr%d%nfu2
 wfs%ks%lzd%Glr%d%nfl3     = wvl%Glr%d%nfl3
 wfs%ks%lzd%Glr%d%nfu3     = wvl%Glr%d%nfu3
 wfs%ks%lzd%Glr%d%n1i      = wvl%Glr%d%n1i
 wfs%ks%lzd%Glr%d%n2i      = wvl%Glr%d%n2i
 wfs%ks%lzd%Glr%d%n3i      = wvl%Glr%d%n3i
 wfs%ks%lzd%Glr%outofzone  = wvl%Glr%outofzone

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) wvl%h(1),wfs%ks
#endif

end subroutine wvl_wfs_lr_copy
!!***

end module m_wvl_wfs
!!***
