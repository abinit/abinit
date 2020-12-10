!!****m* ABINIT/work_var_lotf
!! NAME
!! work_var_lotf
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2005-2020 ABINIT group (MMancini)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module work_var_lotf

 use defs_basis
 use defs_param_lotf
 use m_errors
 use m_abicore

 implicit none

 public ::             &
   work_var_set,       &
   work_var_dealloc,   &
   cutoff_init,        &
   smallfit

 public

 !--Control variables
 ! !--Atomflags variables
 integer,allocatable ::  ifixed(:) !--MMANCINI what is its utility

 !--Quantflags variables
 integer,allocatable :: iq(:)

 !--Cutoff variables
 real(dp) :: rcut,rcut_nbl
 real(dp) :: rcrust

contains

!!***

!!****f* work_var_lotf/work_var_set
!! NAME
!! work_var_set
!!
!! FUNCTION
!!  set some internal variable of lotf
!! INPUTS
!!  natom=number of atoms
!!
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!      dist_pbc
!!
!! SOURCE

 subroutine work_var_set()

  implicit none
  !Local-----------------------------
  integer :: iat

! *************************************************************************

  !--ifixed from ATOMFLAGS is initialized :
   ABI_MALLOC(ifixed,(lotfvar%natom))
   ifixed(:) = 1

  !--FINDS  FITTED ATOMS
  ! ABI_MALLOC(tquant,(lotfvar%natom))
  ! tquant(:) = .true.
  !  nquant = lotfvar%natom
  !nqxx   = lotfvar%natom

   ABI_MALLOC(iq,(lotfvar%natom))
   iq(:)=(/(iat,iat=1,lotfvar%natom)/)

 end subroutine work_var_set
 !!***


!!****f* work_var_lotf/work_var_dealloc
!! NAME
!! work_var_dealloc
!!
!! FUNCTION
!!  deallocate variable
!!
!! INPUTS
!!
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!      dist_pbc
!!
!! SOURCE

 subroutine work_var_dealloc()

! *************************************************************************
   ABI_FREE(iq)
   ABI_FREE(ifixed)
 end subroutine work_var_dealloc
 !!***



!!****f* work_var_lotf/cutoff_init
!! NAME
!!  cutoff_init
!! FUNCTION
!!
!! INPUTS
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!      dist_pbc
!!
!! SOURCE

 subroutine cutoff_init()
  use pbc_lotf,only : pbc_bb_contract
  !Local ---------------------------
  real(dp)  :: bl(3),blmin
  character(len=500) :: message

! *************************************************************************

   if (lotfvar%classic==5 .OR. lotfvar%classic==6) then
     rcut = 4.0d0 / 0.529177 !--check consistency with glue parameters (mind d and its "limits"!)
     rcut_nbl = rcut + rcrust
   end if

   if (lotfvar%me==1.AND.lotfvar%classic==5) then
     write(message,'(2(a,f12.6,a))')&
&     'GLUE radial cutoff used: ',rcut,ch10,&
&     'GLUE NEIGBOURS cutoff used: ',rcut_nbl,ch10
     call wrtout(std_out,message,'COLL')
   end if

  !--Cut-off check respect to cell size :
   bl(:) = pbc_bb_contract()

   if (lotfvar%me==1) then
     write(message,'(3a,3f12.6,a)')&
&     'LENGTH OF REC. CELL VECTORS : ',ch10,&
&     ' bl1, bl2, bl3 : ',   bl(:),ch10
     call wrtout(std_out,message,'COLL')
   end if
   blmin = minval(bl)

   if (rcut_nbl*blmin  > half) then
     write(message,'(2a,2(a,f12.6))')&
&     'LOTF: cut off too large : ',ch10,&
&     ' cut-off (A) is ', rcut_nbl ,  ' min. allowed : ',half/blmin
     ABI_ERROR(message)
   end if
 end subroutine cutoff_init
 !!***


!!****f* work_var_lotf/smallfit
!! NAME
!! smallfit
!!
!! FUNCTION
!!
!! INPUTS
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!      dist_pbc
!!
!! SOURCE

 subroutine smallfit(tau0,ndum)
  use bond_lotf,only : tafit
  USE pbc_lotf,only : dist_pbc,r2
  implicit none

  !Arguments ------------------------
  integer,intent(out) :: ndum
  real(dp),intent(in):: tau0(3,lotfvar%natom)

  !Local ---------------------------
  integer  ::  i, j, jat, nquant
  real(dp) :: rag_fit

! *************************************************************************

   rag_fit = zero

   nquant = lotfvar%natom !--to remember old notation with nquant

   ndum = 0
   do i = 1, lotfvar%natom
     if(.not.tafit(i)) then
       do j=1,nquant
         jat = iq(j)
         call dist_pbc(tau0(:,i),tau0(:,jat))
         if (r2  <  rag_fit) then
           tafit(i) = .true.
           ndum = ndum + 1
           cycle
         end if
       end do !nqtot
     end if !tafit
   end do

 end subroutine smallfit

end module work_var_lotf
!!***
