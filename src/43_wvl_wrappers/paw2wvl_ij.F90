!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw2wvl_ij
!! NAME
!!  paw2wvl_ij
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2016 ABINIT group (TR, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      nullify_paw_ij_objects
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine paw2wvl_ij(option,paw_ij,wvl)
    
 use defs_basis
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
 use m_paw_ij, only : paw_ij_type

#if defined HAVE_BIGDFT
 use BigDFT_API, only : nullify_paw_ij_objects
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'paw2wvl_ij'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in)::option
 type(wvl_internal_type), intent(inout)::wvl
 type(paw_ij_type),intent(in) :: paw_ij(:)
!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: iatom,iaux,my_natom
 character(len=500) :: message
#endif
 
! *************************************************************************
 
 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT
 my_natom=size(paw_ij)

!Option==1: allocate and copy
 if(option==1) then
   ABI_DATATYPE_ALLOCATE(wvl%paw%paw_ij,(my_natom))
   do iatom=1,my_natom
     call nullify_paw_ij_objects(wvl%paw%paw_ij(iatom))
     wvl%paw%paw_ij(iatom)%cplex          =paw_ij(iatom)%cplex
     wvl%paw%paw_ij(iatom)%cplex_dij      =paw_ij(iatom)%cplex_dij
     wvl%paw%paw_ij(iatom)%has_dij        =paw_ij(iatom)%has_dij
     wvl%paw%paw_ij(iatom)%has_dijfr      =0
     wvl%paw%paw_ij(iatom)%has_dijhartree =0
     wvl%paw%paw_ij(iatom)%has_dijhat     =0
     wvl%paw%paw_ij(iatom)%has_dijso      =0
     wvl%paw%paw_ij(iatom)%has_dijU       =0
     wvl%paw%paw_ij(iatom)%has_dijxc      =0
     wvl%paw%paw_ij(iatom)%has_dijxc_val  =0
     wvl%paw%paw_ij(iatom)%has_exexch_pot =0
     wvl%paw%paw_ij(iatom)%has_pawu_occ   =0
     wvl%paw%paw_ij(iatom)%lmn_size       =paw_ij(iatom)%lmn_size
     wvl%paw%paw_ij(iatom)%lmn2_size      =paw_ij(iatom)%lmn2_size
     wvl%paw%paw_ij(iatom)%ndij           =paw_ij(iatom)%ndij
     wvl%paw%paw_ij(iatom)%nspden         =paw_ij(iatom)%nspden
     wvl%paw%paw_ij(iatom)%nsppol         =paw_ij(iatom)%nsppol
     if (paw_ij(iatom)%has_dij/=0) then
       iaux=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
       ABI_ALLOCATE(wvl%paw%paw_ij(iatom)%dij,(iaux,paw_ij(iatom)%ndij))
       wvl%paw%paw_ij(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)
     end if
   end do

!  Option==2: deallocate
 elseif(option==2) then
   do iatom=1,my_natom
     wvl%paw%paw_ij(iatom)%has_dij=0
     if (associated(wvl%paw%paw_ij(iatom)%dij)) then
       ABI_DEALLOCATE(wvl%paw%paw_ij(iatom)%dij)
     end if
   end do
   ABI_DATATYPE_DEALLOCATE(wvl%paw%paw_ij)

!  Option==3: only copy
 elseif(option==3) then
   do iatom=1,my_natom  
     wvl%paw%paw_ij(iatom)%cplex     =paw_ij(iatom)%cplex
     wvl%paw%paw_ij(iatom)%cplex_dij =paw_ij(iatom)%cplex_dij
     wvl%paw%paw_ij(iatom)%lmn_size  =paw_ij(iatom)%lmn_size
     wvl%paw%paw_ij(iatom)%lmn2_size =paw_ij(iatom)%lmn2_size
     wvl%paw%paw_ij(iatom)%ndij      =paw_ij(iatom)%ndij
     wvl%paw%paw_ij(iatom)%nspden    =paw_ij(iatom)%nspden
     wvl%paw%paw_ij(iatom)%nsppol    =paw_ij(iatom)%nsppol
     wvl%paw%paw_ij(iatom)%dij(:,:)  =paw_ij(iatom)%dij(:,:)
   end do 
   
 else 
   message = 'paw2wvl_ij: option should be equal to 1, 2 or 3'
   MSG_ERROR(message)
 end if

#else
 if (.false.) write(std_out,*) option,wvl%h(1),paw_ij(1)%ndij
#endif

 DBG_EXIT("COLL")

end subroutine paw2wvl_ij
!!***
