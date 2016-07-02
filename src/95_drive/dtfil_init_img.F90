!{\src2tex{textfont=tt}}
!!****f* ABINIT/dtfil_init_img
!! NAME
!! dtfil_init_img
!!
!! FUNCTION
!! Initialize few scalars in the dtfil structured variable
!! when an alogrithm using image of the cell is selected.
!! (initialize index of images from which read files)
!!
!! COPYRIGHT
!! Copyright (C) 2011-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset=<type datasets_type>=input variables for the current dataset
!!  dtsets(0:ndtset_alloc)=<type datasets_type>=input variables for all datasets
!!  idtset=number of the dataset
!!  jdtset(0:ndtset)=actual index of the datasets
!!  ndtset=number of datasets
!!  ndtset_alloc=number of datasets, corrected for allocation of at least one data set
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! dtfil=<type datafiles_type>= only getxxx_from_image flags are modified
!!
!! NOTES
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dtfil_init_img(dtfil,dtset,dtsets,idtset,jdtset,ndtset,ndtset_alloc)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtfil_init_img'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: idtset,ndtset,ndtset_alloc
 type(datafiles_type),intent(out) :: dtfil
 type(dataset_type),intent(in) :: dtset
!arrays
 integer :: jdtset(0:ndtset)
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)

!Local variables -------------------------
!scalars
 integer :: iget
!arrays

! *********************************************************************

 DBG_ENTER("COLL")

!Default values
 dtfil%getwfk_from_image   =0 ! Get standard WFK from previous dataset
 dtfil%getden_from_image   =0 ! Get standard DEN from previous dataset
 dtfil%getpawden_from_image=0 ! Get standard PAWDEN from previous dataset

 if (dtset%optdriver==RUNL_GSTATE.and.dtset%nimage>1) then

!  Define getwfk_from_image
   if (dtset%getwfk/=0.or.dtset%irdwfk/=0) then
     iget=-1
     if(dtset%getwfk<0) iget=jdtset(idtset+dtset%getwfk)
     if(dtset%getwfk>0) iget=dtset%getwfk
     if(dtset%irdwfk>0) iget=0
     if (iget>=0) then
       if (iget==0.or.dtsets(iget)%nimage==dtset%nimage) then
         dtfil%getwfk_from_image=-1     ! Get WFK from the same image of previous dataset
       else if (dtsets(iget)%nimage>1) then
         dtfil%getwfk_from_image=1      ! Get WFK from the first image of previous dataset
       end if
     end if
   end if

!  Define getden_from_image
   if (dtset%getden/=0.or.dtset%irdden/=0) then
     iget=-1
     if(dtset%getden<0) iget=jdtset(idtset+dtset%getden)
     if(dtset%getden>0) iget=dtset%getden
     if(dtset%irdden>0) iget=0
     if (iget>=0) then
       if (iget==0.or.dtsets(iget)%nimage==dtset%nimage) then
         dtfil%getden_from_image=-1     ! Get DEN from the same image of previous dataset
       else if (dtsets(iget)%nimage>1) then
         dtfil%getden_from_image=1      ! Get DEN from the first image of previous dataset
       end if
     end if
   end if

!  Define getpawden_from_image
   if (dtset%getpawden/=0.or.dtset%irdpawden/=0) then
     iget=-1
     if(dtset%getpawden<0) iget=jdtset(idtset+dtset%getpawden)
     if(dtset%getpawden>0) iget=dtset%getpawden
     if(dtset%irdpawden>0) iget=0
     if (iget>=0) then
       if (iget==0.or.dtsets(iget)%nimage==dtset%nimage) then
         dtfil%getpawden_from_image=-1     ! Get PAWDEN from the same image of previous dataset
       else if (dtsets(iget)%nimage>1) then
         dtfil%getpawden_from_image=1      ! Get PAWDEN from the first image of previous dataset
       end if
     end if
   end if
 end if

 DBG_EXIT("COLL")

end subroutine dtfil_init_img
!!***
