!!****m* ABINIT/m_args_gs
!! NAME
!!  m_args_gs
!!
!! FUNCTION
!!  This module provides the definition of the args_gs_type
!!  used to tranfer some arguments to GS calculations,
!!  especially those depending on the image of the cell.
!!
!! COPYRIGHT
!! Copyright (C) 2015-2020 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_args_gs

 use defs_basis
 use m_abicore
 use m_errors

 implicit none

 private
!!***

!!****t* m_args_gs/args_gs_type
!! NAME
!! args_gs_type
!!
!! FUNCTION
!! This structured datatype contains some arguments of a GS calculation,
!! especially the ones depending on the "images".
!!
!! SOURCE

 type, public :: args_gs_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Real (real(dp)) arrays

  real(dp), pointer :: amu(:)
   ! amu(ntypat)
   ! mass of each atom type

  real(dp), pointer :: mixalch(:,:)
   ! mixalch(ntypat)
   ! Mixing coefficients to generate alchemical pseudo atoms

  real(dp), pointer :: dmatpawu(:,:,:,:)
   ! dmatpawu(2*lpawu+1,2*lpawu+1,nspinor*nsppol,natpawu)
   ! Fixed occupation matrix for correlated orbitals (DFT+U or DMFT only)

  real(dp), pointer :: upawu(:)
   ! upawu(ntypat)
   ! Value of U for the DFT+U or DMFT approach

  real(dp), pointer :: jpawu(:)
   ! upawu(ntypat)
   ! Value of J for the DFT+U or DMFT approach

  real(dp), pointer :: rprimd_orig(:,:)
   ! rprimd_orig(3,3)
   ! Original primitive vectors (usually the input variable)

 end type args_gs_type

!public procedures.
 public :: args_gs_init
 public :: args_gs_free
!!***

CONTAINS

!===========================================================
!!***

!!****f* m_args_gs/args_gs_init
!! NAME
!!  args_gs_init
!!args_gs_init
!! FUNCTION
!!  Init a args_gs datastructure
!!
!! INPUTS
!!  amu(:)= mass of each atom type
!!  mixalch(:)= mixing coefficients to generate alchemical pseudo atoms
!!  dmatpawu(:,:,:,:)= fixed occupation matrix for correlated orbitals (DFT+U or DMFT only)
!!  upawu(:) =value of U for the DFT+U or DMFT approach
!!  jpawu(:) =value of J for the DFT+U or DMFT approach
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  args_gs=<type(args_gs_type)>=args_gs datastructure
!!
!! PARENTS
!!      m_gstateimg
!!
!! CHILDREN
!!
!! SOURCE

subroutine args_gs_init(args_gs,amu,mixalch,dmatpawu,upawu,jpawu,rprimd_orig)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in),target :: amu(:),dmatpawu(:,:,:,:),jpawu(:),mixalch(:,:),upawu(:)
 real(dp),intent(in),target :: rprimd_orig(:,:)
 type(args_gs_type),intent(inout) :: args_gs
!Local variables-------------------------------

!************************************************************************

 !@args_gs_type

 args_gs%amu         => amu
 args_gs%mixalch     => mixalch
 args_gs%dmatpawu    => dmatpawu
 args_gs%upawu       => upawu
 args_gs%jpawu       => jpawu
 args_gs%rprimd_orig => rprimd_orig

end subroutine args_gs_init
!!***

!----------------------------------------------------------------------

!!****f* m_args_gs/args_gs_free
!! NAME
!!  args_gs_free
!!
!! FUNCTION
!!  Clean and destroy a args_gs datastructure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  args_gs(:)=<type(args_gs_type)>=args_gs datastructure
!!
!! PARENTS
!!      m_gstateimg
!!
!! CHILDREN
!!
!! SOURCE

subroutine args_gs_free(args_gs)

!Arguments ------------------------------------
!arrays
 type(args_gs_type),intent(inout) :: args_gs
!Local variables-------------------------------

!************************************************************************

 !@args_gs_type

 args_gs%amu         => null()
 args_gs%mixalch     => null()
 args_gs%dmatpawu    => null()
 args_gs%upawu       => null()
 args_gs%jpawu       => null()
 args_gs%rprimd_orig => null()

end subroutine args_gs_free
!!***

!----------------------------------------------------------------------

END MODULE m_args_gs
!!***
