!!****m* ABINIT/m_wvl_projectors
!! NAME
!!  m_wvl_projectors
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

module m_wvl_projectors

 use defs_basis
 use defs_wvltypes
 use m_errors
 use m_abicore
 use m_atomdata

 use defs_datatypes, only : pseudopotential_type
 use m_geometry,     only : xred2xcart

 implicit none

 private
!!***

 public :: wvl_projectors_set
 public :: wvl_projectors_free
!!***

contains
!!***

!!****f* ABINIT/wvl_projectors_set
!!
!! NAME
!! wvl_projectors_set
!!
!! FUNCTION
!! Allocate and compute the access keys for the projectors when the positions
!! of the atoms are given. The array to store projectors
!! is also allocated, use wvl_projectors_free() to free them after use.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=internal variables used by wavelets, describing
!!   | wvl_internal=desciption of the wavelet box.
!!   | natom=number of atoms.
!!  mpi_enreg=information about MPI parallelization
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  proj <type(wvl_projector_type)>=projectors information for wavelets.
!!   | keys=its access keys for compact storage.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_gstate,m_wvl_wfsinp
!!
!! CHILDREN
!!      free_dft_psp_projectors
!!
!! SOURCE

subroutine wvl_projectors_set(me, natom, proj, psps, rprimd, wfs, wvl, wvl_frmult, xred)

#if defined HAVE_BIGDFT
 use BigDFT_API, only: createProjectorsArrays, wvl_timing => timing
#endif
 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom, me
 real(dp), intent(in) :: wvl_frmult
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_projectors_type),intent(inout) :: proj
 type(wvl_wf_type),intent(in) :: wfs
 type(wvl_internal_type), intent(in) :: wvl
!arrays
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)

!Local variables-------------------------------
!scalars
#if defined HAVE_BIGDFT
 integer :: idata
 logical,parameter :: wvl_debug=.false.
 character(len=500) :: message
!arrays
 real(dp),allocatable :: xcart(:,:)
#endif

! *********************************************************************

#if defined HAVE_BIGDFT
!Consistency checks, are all pseudo true GTH pseudo with geometric information?
 do idata = 1, psps%npsp, 1
   if (.not. psps%gth_params%set(idata)) then
     write(message, '(a,a,a,a,I0,a,a,a)' ) ch10,&
&     ' wvl_projectors_set :  consistency checks failed,', ch10, &
&     '  no GTH parameters found for type number ', idata, '.', ch10, &
&     '  Check your input pseudo files.'
     MSG_ERROR(message)
   end if
   if (.not. psps%gth_params%hasGeometry(idata)) then
     write(message, '(a,a,a,a,a,a)' ) ch10,&
&     ' wvl_projectors_set :  consistency checks failed,', ch10, &
&     '  the given GTH parameters has no geometry information.', ch10, &
&     '  Upgrade your input pseudo files to GTH with geometric informatoins.'
     MSG_ERROR(message)
   end if
 end do

 if (wvl_debug) then
   call wvl_timing(me,'CrtProjectors ','ON')
 end if

!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, natom))
 call xred2xcart(natom, rprimd, xcart, xred)
 call createProjectorsArrays(wfs%ks%Lzd%Glr,xcart,wvl%atoms,wfs%ks%orbs,&
 psps%gth_params%radii_cf,wvl_frmult,wvl_frmult,wvl%h(1),wvl%h(2),&
 wvl%h(3),.false.,proj%nlpsp,proj%G)
 write(message, '(a,a,a,a,I0)' ) ch10,&
& ' wvl_projectors_set : allocate projectors data,', ch10, &
& '  size of the compressed array: ', proj%nlpsp%nprojel
 call wrtout(std_out,message,'COLL')

!Deallocations
 ABI_DEALLOCATE(xcart)

 if (wvl_debug) then
   call wvl_timing(me,'CrtProjectors ','OF')
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) natom,me,wvl_frmult,psps%npsp,proj%nlpsp,wfs%ks,wvl%h(1),&
& rprimd(1,1),xred(1,1)
#endif

end subroutine wvl_projectors_set
!!***

!!****f* ABINIT/wvl_projectors_free
!!
!! NAME
!! wvl_projectors_free
!!
!! FUNCTION
!! Freeing routine.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  proj <type(wvl_projectors_type)>=projectors information in a wavelet basis.
!!
!! PARENTS
!!      m_gstate,m_wvl_wfsinp
!!
!! CHILDREN
!!      free_dft_psp_projectors
!!
!! SOURCE

subroutine wvl_projectors_free(proj)

#if defined HAVE_BIGDFT
  use BigDFT_API, only : free_DFT_PSP_projectors,deallocate_gwf
#endif
 implicit none

!Arguments ------------------------------------
!scalars
 type(wvl_projectors_type),intent(inout) :: proj

!Local variables -------------------------
#if defined HAVE_BIGDFT
 integer :: ii
#endif

  ! *********************************************************************

#if defined HAVE_BIGDFT

 call free_DFT_PSP_projectors(proj%nlpsp)

 if (allocated(proj%G)) then
   do ii=1,size(proj%G)
!    MT dec 2014: cannot call bigdft deallocation routine
!    because content of proj%G datastructure was created
!    without f_malloc (without memory profiling).
!    call deallocate_gwf(proj%G(ii))
     if (associated(proj%G(ii)%ndoc)) then
       ABI_DEALLOCATE(proj%G(ii)%ndoc)
     end if
     if (associated(proj%G(ii)%nam)) then
       ABI_DEALLOCATE(proj%G(ii)%nam)
     end if
     if (associated(proj%G(ii)%nshell)) then
       ABI_DEALLOCATE(proj%G(ii)%nshell)
     end if
     if (associated(proj%G(ii)%psiat)) then
       ABI_DEALLOCATE(proj%G(ii)%psiat)
     end if
     if (associated(proj%G(ii)%xp)) then
       ABI_DEALLOCATE(proj%G(ii)%xp)
     end if
   end do
   ABI_DATATYPE_DEALLOCATE(proj%G)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) proj%nlpsp
#endif

end subroutine wvl_projectors_free
!!***

end module m_wvl_projectors
!!***
