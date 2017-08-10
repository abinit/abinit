!{\src2tex{textfont=tt}}
!!****f* ABINIT/macroin2
!! NAME
!! macroin2
!!
!! FUNCTION
!! Treat "macro" input variables, that can :
!! - initialize several other input variables for one given dataset
!! - initialize several other input variables for a set of datasets.
!! Note that the treatment of these different types of macro input variables is different.
!! Documentation of such input variables is very important, including the
!! proper echo, in the output file, of what such input variables have done. 
!!
!! Important information : all the "macro" input variables should be properly
!! identifiable to be so, and it is proposed to make them start with the string "macro".
!!
!! COPYRIGHT
!! Copyright (C) 2009-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ndtset_alloc=number of datasets, corrected for allocation of at
!!               least one data set.
!!
!! OUTPUT
!!  dtsets(0:ndtset_alloc)=<type datafiles_type>contains all input variables,
!!   some of which are given a value here.
!!   The dataset with number 0 should NOT be modified in the present routine.
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine macroin2(dtsets,ndtset_alloc)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'macroin2'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndtset_alloc
!arrays
 type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc) !vz_i

!Local variables -------------------------------
!scalars
 integer :: idtset,pawujat              !,jdtset

!******************************************************************

 do idtset=1,ndtset_alloc
!  Set first PAW+U atom to perform atomic level shift
   if (dtsets(idtset)%typat(1)==0) cycle
   pawujat=dtsets(idtset)%pawujat
   pawujat=pawujat-count(dtsets(idtset)%lpawu( dtsets(idtset)%typat( 1:pawujat ))<0)

   if (dtsets(idtset)%macro_uj>0) then
!    Level shift atom with amplitude pawujv 
     dtsets(idtset)%atvshift(:,:,pawujat)=dtsets(idtset)%pawujv

!    Case level shift only on one spin channel
     if ((dtsets(idtset)%macro_uj==2.or.dtsets(idtset)%macro_uj==3).and.dtsets(idtset)%nsppol==2) then 
       dtsets(idtset)%atvshift(:,2,pawujat)=0_dp 
     end if

   end if ! macro_uj

 end do

end subroutine macroin2
!!***
