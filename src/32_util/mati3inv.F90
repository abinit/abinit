!{\src2tex{textfont=tt}}
!!****f* ABINIT/mati3inv
!! NAME
!! mati3inv
!!
!! FUNCTION
!! Invert and transpose orthogonal 3x3 matrix of INTEGER elements.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! mm = integer matrix to be inverted
!!
!! OUTPUT
!! mit = inverse of mm input matrix
!!
!! NOTES
!! Used for symmetry operations.
!! This routine applies to ORTHOGONAL matrices only.
!! Since these form a group, inverses are also integer
!! arrays.  Returned array is TRANSPOSE of inverse, as needed.
!! Note use of integer arithmetic.
!! Also: has been designed so that mit can be same storage space as m, in
!! which case m is overwritten by resulting mit.
!!
!! PARENTS
!!      cg_rotate,chkgrp,classify_bands,debug_tools,dfpt_nstdy,get_full_kgrid
!!      get_npert_rbz,getkgrid,ingeo,m_ab7_symmetry,m_crystal,m_ddb,m_dvdb
!!      m_dynmat,m_fft_mesh,m_fock,m_ptgroups,matpointsym,memory_eval,optic
!!      read_gkk,setsym,strainsym,thmeig,wfconv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mati3inv(mm,mit)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mati3inv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 integer,intent(in) :: mm(3,3)
 integer,intent(out) :: mit(3,3)

!Local variables-------------------------------
!scalars
 integer :: dd
 character(len=500) :: message
!arrays
 integer :: tt(3,3)

! *************************************************************************

 tt(1,1) = mm(2,2) * mm(3,3) - mm(3,2) * mm(2,3)
 tt(2,1) = mm(3,2) * mm(1,3) - mm(1,2) * mm(3,3)
 tt(3,1) = mm(1,2) * mm(2,3) - mm(2,2) * mm(1,3)
 tt(1,2) = mm(3,1) * mm(2,3) - mm(2,1) * mm(3,3)
 tt(2,2) = mm(1,1) * mm(3,3) - mm(3,1) * mm(1,3)
 tt(3,2) = mm(2,1) * mm(1,3) - mm(1,1) * mm(2,3)
 tt(1,3) = mm(2,1) * mm(3,2) - mm(3,1) * mm(2,2)
 tt(2,3) = mm(3,1) * mm(1,2) - mm(1,1) * mm(3,2)
 tt(3,3) = mm(1,1) * mm(2,2) - mm(2,1) * mm(1,2)
 dd  = mm(1,1) * tt(1,1) + mm(2,1) * tt(2,1) + mm(3,1) * tt(3,1)

!Make sure matrix is not singular
 if (dd/=0) then
   mit(:,:)=tt(:,:)/dd
 else
   write(message, '(2a,2x,9i5,a)' )&
&   ' Attempting to invert integer array',ch10,&
&   mm(:,:),'   ==> determinant is zero.'
   MSG_BUG(message)
 end if

!If matrix is orthogonal, determinant must be 1 or -1
 if (abs(dd)/=1) then
   write(message, '(2a,2x,9i5,a)' )&
&   ' Absolute value of determinant should be one',ch10,&
&   ' but determinant=',dd
   MSG_BUG(message)
 end if


end subroutine mati3inv
!!***
