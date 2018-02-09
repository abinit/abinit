!{\src2tex{textfont=tt}}
!!****f* ABINIT/pspheads_comm
!! NAME
!! pspheads_comm
!!
!! FUNCTION
!! Communicate pspheads to all processors
!!
!! COPYRIGHT
!! Copyright (C) 2009-2018 ABINIT group (DCA, XG, GMR, FrD, AF, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  npsp=number of pseudopotentials
!!  test_paw=0 if no PAW, 1 if PAW
!!
!! SIDE EFFECTS
!!  pspheads(npsp)=<type pspheader_type>=all the important information from the
!!   pseudopotential file headers, as well as the psp file names. On one processor at input,
!!   on all processors at output
!!
!! PARENTS
!!      m_ab7_invars_f90
!!
!! CHILDREN
!!      timab,xmpi_bcast
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pspheads_comm(npsp,pspheads,test_paw)

 use defs_basis
 use defs_datatypes
 use m_profiling_abi
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pspheads_comm'
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: npsp
 integer,intent(inout) :: test_paw
 type(pspheader_type),intent(inout) :: pspheads(npsp)

!Local variables-------------------------------
#if defined HAVE_MPI
!scalars
 integer,parameter :: master=0
 integer :: ierr,comm
!arrays
 integer,allocatable :: list_int(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: list_dpr(:)
 character(len=fnlen),allocatable :: list_char(:)
#endif

!*************************************************************************

#if defined HAVE_MPI
 call timab(48,1,tsec)

 comm = xmpi_world

!Broadcast the characters (file names and titles)
 ABI_ALLOCATE(list_char,(3*npsp))
 list_char(1:npsp)=pspheads(1:npsp)%filpsp
 list_char(npsp+1:2*npsp)=pspheads(1:npsp)%title
 list_char(2*npsp+1:3*npsp)=pspheads(1:npsp)%md5_checksum

 call xmpi_bcast(list_char,master,comm,ierr)

 pspheads(1:npsp)%filpsp=list_char(1:npsp)
 pspheads(1:npsp)%title=list_char(npsp+1:2*npsp)
 pspheads(1:npsp)%md5_checksum=list_char(2*npsp+1:3*npsp)(1:md5_slen)
 ABI_DEALLOCATE(list_char)

!Brodcast the integers
 ABI_ALLOCATE(list_int,(1+13*npsp))
 list_int(1        :   npsp) = pspheads(1:npsp)%nproj(0)
 list_int(1+   npsp: 2*npsp) = pspheads(1:npsp)%nproj(1)
 list_int(1+ 2*npsp: 3*npsp) = pspheads(1:npsp)%nproj(2)
 list_int(1+ 3*npsp: 4*npsp) = pspheads(1:npsp)%nproj(3)
 list_int(1+ 4*npsp: 5*npsp) = pspheads(1:npsp)%lmax
 list_int(1+ 5*npsp: 6*npsp) = pspheads(1:npsp)%xccc
 list_int(1+ 6*npsp: 7*npsp) = pspheads(1:npsp)%pspxc
 list_int(1+ 7*npsp: 8*npsp) = pspheads(1:npsp)%pspdat
 list_int(1+ 8*npsp: 9*npsp) = pspheads(1:npsp)%pspcod
 list_int(1+ 9*npsp:10*npsp) = pspheads(1:npsp)%pspso
 list_int(1+10*npsp:11*npsp) = pspheads(1:npsp)%nprojso(1)
 list_int(1+11*npsp:12*npsp) = pspheads(1:npsp)%nprojso(2)
 list_int(1+12*npsp:13*npsp) = pspheads(1:npsp)%nprojso(3)
 list_int(1+13*npsp)         = test_paw

 call xmpi_bcast(list_int,master,comm,ierr)

 pspheads(1:npsp)%nproj(0)   = list_int(1        :   npsp)
 pspheads(1:npsp)%nproj(1)   = list_int(1+   npsp: 2*npsp)
 pspheads(1:npsp)%nproj(2)   = list_int(1+ 2*npsp: 3*npsp)
 pspheads(1:npsp)%nproj(3)   = list_int(1+ 3*npsp: 4*npsp)
 pspheads(1:npsp)%lmax       = list_int(1+ 4*npsp: 5*npsp)
 pspheads(1:npsp)%xccc       = list_int(1+ 5*npsp: 6*npsp)
 pspheads(1:npsp)%pspxc      = list_int(1+ 6*npsp: 7*npsp)
 pspheads(1:npsp)%pspdat     = list_int(1+ 7*npsp: 8*npsp)
 pspheads(1:npsp)%pspcod     = list_int(1+ 8*npsp: 9*npsp)
 pspheads(1:npsp)%pspso      = list_int(1+ 9*npsp:10*npsp)
 pspheads(1:npsp)%nprojso(1) = list_int(1+10*npsp:11*npsp)
 pspheads(1:npsp)%nprojso(2) = list_int(1+11*npsp:12*npsp)
 pspheads(1:npsp)%nprojso(3) = list_int(1+12*npsp:13*npsp)
 test_paw                    = list_int(1+13*npsp)
 ABI_DEALLOCATE(list_int)

!Unbeliveable, this cannot be sent with the others, for woopy
 ABI_ALLOCATE(list_int,(npsp))
 list_int(1:npsp) = pspheads(1:npsp)%usewvl
 call xmpi_bcast(list_int,master,comm,ierr)
 pspheads(1:npsp)%usewvl     = list_int(1:npsp)
 ABI_DEALLOCATE(list_int)


!Broadcast zionpsp and znuclpsp
 ABI_ALLOCATE(list_dpr,(7*npsp))
 list_dpr(1       :  npsp) = pspheads(1:npsp)%zionpsp
 list_dpr(1+  npsp:2*npsp) = pspheads(1:npsp)%znuclpsp
 list_dpr(1+2*npsp:3*npsp) = pspheads(1:npsp)%GTHradii(0)
 list_dpr(1+3*npsp:4*npsp) = pspheads(1:npsp)%GTHradii(1)
 list_dpr(1+4*npsp:5*npsp) = pspheads(1:npsp)%GTHradii(2)
 list_dpr(1+5*npsp:6*npsp) = pspheads(1:npsp)%GTHradii(3)
 list_dpr(1+6*npsp:7*npsp) = pspheads(1:npsp)%GTHradii(4)
 
 call xmpi_bcast(list_dpr,master,comm,ierr)

 pspheads(1:npsp)%zionpsp     = list_dpr(1       :  npsp)
 pspheads(1:npsp)%znuclpsp    = list_dpr(1+  npsp:2*npsp)
 pspheads(1:npsp)%GTHradii(0) = list_dpr(1+2*npsp:3*npsp)
 pspheads(1:npsp)%GTHradii(1) = list_dpr(1+3*npsp:4*npsp)
 pspheads(1:npsp)%GTHradii(2) = list_dpr(1+4*npsp:5*npsp)
 pspheads(1:npsp)%GTHradii(3) = list_dpr(1+5*npsp:6*npsp)
 pspheads(1:npsp)%GTHradii(4) = list_dpr(1+6*npsp:7*npsp)
 ABI_DEALLOCATE(list_dpr)

!Broadcast additional integers for PAW psps (testpaw was spread, previously)
 if (test_paw==1) then
   ABI_ALLOCATE(list_int,(6*npsp))
   list_int(1       :  npsp)=pspheads(1:npsp)%pawheader%basis_size
   list_int(1+  npsp:2*npsp)=pspheads(1:npsp)%pawheader%l_size
   list_int(1+2*npsp:3*npsp)=pspheads(1:npsp)%pawheader%lmn_size
   list_int(1+3*npsp:4*npsp)=pspheads(1:npsp)%pawheader%mesh_size
   list_int(1+4*npsp:5*npsp)=pspheads(1:npsp)%pawheader%pawver
   list_int(1+5*npsp:6*npsp)=pspheads(1:npsp)%pawheader%shape_type

   call xmpi_bcast(list_int,master,comm,ierr)

   pspheads(1:npsp)%pawheader%basis_size=list_int(1       :  npsp)
   pspheads(1:npsp)%pawheader%l_size    =list_int(1+  npsp:2*npsp)
   pspheads(1:npsp)%pawheader%lmn_size  =list_int(1+2*npsp:3*npsp)
   pspheads(1:npsp)%pawheader%mesh_size =list_int(1+3*npsp:4*npsp)
   pspheads(1:npsp)%pawheader%pawver    =list_int(1+4*npsp:5*npsp)
   pspheads(1:npsp)%pawheader%shape_type=list_int(1+5*npsp:6*npsp)
   ABI_DEALLOCATE(list_int)

!  broadcast rpaw values
   ABI_ALLOCATE(list_dpr,(2*npsp))

   list_dpr(1       :  npsp) = pspheads(1:npsp)%pawheader%rpaw
   list_dpr(1+1*npsp:2*npsp) = pspheads(1:npsp)%pawheader%rshp

   call xmpi_bcast(list_dpr,master,comm,ierr)

   pspheads(1:npsp)%pawheader%rpaw = list_dpr(1       :  npsp)
   pspheads(1:npsp)%pawheader%rshp = list_dpr(1+  npsp:2*npsp)

   ABI_DEALLOCATE(list_dpr)
 end if

 call timab(48,2,tsec)

#else
!Code to use unused dummy arguments
 if(pspheads(1)%lmax == -10) pspheads(1)%lmax=-10
 if(test_paw == -1) test_paw = -1
#endif

end subroutine pspheads_comm
!!***
