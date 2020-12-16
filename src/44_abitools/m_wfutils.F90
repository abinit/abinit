!!****m* ABINIT/m_wfutils
!! NAME
!!  m_wfutils
!!
!! FUNCTION
!!  parameters and function for wave functions copy
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2020 ABINIT group (CS,GZ,FB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_wfutils

 use defs_basis
 use m_abicore
 use m_errors

 use m_time,     only : timab

 implicit none

 private

 public :: setWFParameter ! Definition of wave functions parameters for copy functions.
 public :: wfcopy         ! Copy of wave functions arrays.

! Global variables
 integer,save,public ABI_PROTECTED :: x_cplx   ! fortran data type for wave functions arrays real or complex.
 integer,save,public ABI_PROTECTED :: x_me_g0  ! processor number
 integer,save,public ABI_PROTECTED :: x_npw_k  ! number of plane waves at this k point
 integer,save,public ABI_PROTECTED :: x_nspinor! number of spinorial components of the wavefunctions on current proc
 integer,save,public ABI_PROTECTED :: x_icg    ! shift to be applied on the location of data in the array cg
 integer,save,public ABI_PROTECTED :: x_igsc   ! shift to be applied on the location of data in the array gsc
 integer,save,public ABI_PROTECTED :: x_blocksize

contains
!!***

!!****f* m_wfutils/setWFParameter
!! NAME
!!
!! setWFParameter
!!
!! FUNCTION
!! Initialize wave functions parameters for copy functions
!!
!! INPUTS
!!  cplx=fortran data type for wave functions arrays real or complex
!!  me_g0=1 if this node treats G=0.
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions on current proc
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array gsc
!!  blocksize=size of blocks
!!
!! PARENTS
!!      lobpcgwf
!!
!! CHILDREN
!!
!! SOURCE
!!
subroutine setWFParameter(cplx,me_g0,npw_k,nspinor,icg,igsc,blocksize)

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: cplx,me_g0,npw_k,nspinor
 integer, intent(in) :: icg,igsc,blocksize

! *********************************************************************

 ! Copy values in global variables
 x_cplx=cplx
 x_me_g0=me_g0
 x_npw_k=npw_k
 x_nspinor=nspinor
 x_icg=icg
 x_igsc=igsc
 x_blocksize=blocksize

end subroutine setWFParameter
!!***

! correspondence with abinit. here for real wf
! this is the index of a given band in cg array
integer function x_cgindex(iblocksize)

 implicit none

 integer, intent(in) :: iblocksize

 x_cgindex=x_npw_k*x_nspinor*(iblocksize-1)+x_icg+1

end function x_cgindex

! correspondence with abinit. here for real wf
! this is the index of a given band in gsc array
integer function x_gscindex(iblocksize)

 implicit none

 integer, intent(in) :: iblocksize

 x_gscindex=x_npw_k*x_nspinor*(iblocksize-1)+x_igsc+1

end function x_gscindex

integer function x_windex(iblocksize)

 implicit none

 integer, intent(in) :: iblocksize

 x_windex=x_npw_k*x_nspinor*(iblocksize-1)+1

end function x_windex

integer function wfindex(iblocksize,indtype)

 implicit none

 integer, intent(in) :: iblocksize
 character(len=1), intent(in) :: indtype

 select case(indtype)
 case ('C')
   wfindex=x_cgindex(iblocksize)
 case ('S')
   wfindex=x_gscindex(iblocksize)
 case ('W')
   wfindex=x_windex(iblocksize)
 case default
   MSG_ERROR("Wrong indtype: "//trim(indtype))
 end select

end function wfindex

!!****f* m_wfutils/wfcopy
!! NAME
!!
!! wfcopy
!!
!! FUNCTION
!! copy of wave functions arrays.
!! called from lobpcg in REAL or COMPLEX wave functions.
!!
!! INPUTS
!!   direction=copy direction
!!                'D' direct (global array to local)
!!                'I' indirect (local to global)
!!   size=number of elements
!!   tsrc=source array
!!   incsrc=size increment for tsrc array
!!   tdest=destination array
!!   incdest=increment for tdest array
!!   blockiter=number of block iteration in case REAL
!!   iblock=block index
!!   indtype=indexation type in array
!!   withbbloc=apply block on band for each block
!!
!! TODO
!!  Split the two cases so that we can avoid the array descriptors.
!!
!! PARENTS
!!   lobpcgwf
!!
!! CHILDREN
!!   dcopy, zcopy
!!
!! SOURCE
!!
subroutine wfcopy(direction,size,tsrc,incsrc,tdest,incdest,blockiter,iblock,indtype,&
&                withbbloc,timopt,tim_wfcopy) ! optional arguments

 implicit none

!Arguments ------------------------------------
 character(len=1), intent(in) :: direction
 integer, intent(in) :: size,incsrc,incdest
 integer, intent(in) :: blockiter,iblock
 character(len=1), intent(in) :: indtype
 logical, optional, intent(in) :: withbbloc
 integer, intent(in), optional :: timopt,tim_wfcopy
!arrays
 real(dp), DEV_CONTARRD intent(in) :: tsrc(:,:)
 real(dp), DEV_CONTARRD intent(inout) :: tdest(:,:)

!Local variables ------------------------------------
 integer,parameter :: ndat1=1
 logical :: bblock=.false.
 integer :: lig,g1,g2,vectsize,rvectsize
 integer :: blocksize,bblocksize,iblocksize,iband
 real(dp) :: factor
 real(dp) :: tsec(2)

! *********************************************************************

 if (present(tim_wfcopy).and.present(timopt)) then
   if(abs(timopt)==3) then
     call timab(tim_wfcopy,1,tsec)
   end if
 end if

 if (present(withbbloc)) bblock=withbbloc

 if (indtype == 'C') then
   lig=x_icg
 else if (indtype == 'S') then
   lig=x_igsc
 else
   lig=0
 endif

 rvectsize=x_npw_k*x_nspinor
 if (x_me_g0 == 1) then
   vectsize=2*rvectsize-1
 else
   vectsize=2*rvectsize
 endif
 if (x_cplx == 2) vectsize=x_npw_k*x_nspinor

 blocksize  = x_blocksize
 bblocksize=(iblock-1)*blocksize

 if (direction == 'D') then

   if (x_cplx == 1) then
     ! Pack real and imag part.
     factor=sqrt(two)
     do iblocksize=1,blockiter
       iband=iblocksize
       if (bblock) then
         iband=iblocksize+bblocksize
       endif
       if (x_me_g0 == 1) then
         tdest(1 ,iblocksize)=tsrc(1,wfindex(iband,indtype))
         g1 = wfindex(iband,  indtype)+1
         g2 = wfindex(iband+1,indtype)-1
         call dcopy(rvectsize-1,tsrc(1,g1:g2),1,tdest(2:rvectsize,iblocksize),1)
         tdest(2:rvectsize,iblocksize) = factor * tdest(2:rvectsize,iblocksize)

         call dcopy(rvectsize-1,tsrc(2,g1:g2),1,tdest(rvectsize+1:vectsize,iblocksize),1)
         tdest(rvectsize+1:vectsize,iblocksize) = factor * tdest(rvectsize+1:vectsize,iblocksize)
         ! MG FIXME: Here gfortran4.9 allocates temporary arrays due to factor.
       else
         g1 = wfindex(iband,  indtype)
         g2 = wfindex(iband+1,indtype)-1
         call dcopy(rvectsize,tsrc(1,g1:g2),1,tdest(1:rvectsize,iblocksize),1)
         tdest(1:rvectsize,iblocksize) = factor * tdest(1:rvectsize,iblocksize)

         call dcopy(rvectsize,tsrc(2,g1:g2),1,tdest(rvectsize+1:vectsize,iblocksize),1)
         tdest(rvectsize+1:vectsize,iblocksize) = factor * tdest(rvectsize+1:vectsize,iblocksize)
       end if
     end do
   else
     if (indtype == 'C') then
       if (bblock) then
         g1 = vectsize*((iblock-1)*blocksize)+lig+1
         g2 = vectsize*(iblock*blocksize)+lig
         call zcopy(size,tsrc(:,g1:g2),incsrc,tdest(:,1:blocksize),incdest)
       else
         g1 = lig+1
         g2 = vectsize*((iblock-1)*blocksize)+lig
         call zcopy(size,tsrc(:,g1:g2),incsrc,tdest(:,1:bblocksize),incdest)
       end if
     else if (indtype == 'S') then
       g1 = lig+1
       g2 = vectsize*((iblock-1)*blocksize)+lig
       call zcopy(size,tsrc(:,g1:g2),incsrc,tdest(:,1:bblocksize),incdest)
     else
       call zcopy(size,tsrc,incsrc,tdest,incdest)
     endif
   end if

 else if (direction == 'I') then
   if (x_cplx == 1) then
     ! Unpack real and imag blocks.
     factor=one/sqrt(two)
     do iblocksize=1,blockiter
       iband=iblocksize
       if ( bblock ) then
         iband=iblocksize+(iblock-1)*blocksize
       end if
       if (x_me_g0 == 1) then
         tdest(1,wfindex(iband,indtype))=tsrc(1,iblocksize)
         tdest(2,wfindex(iband,indtype))=zero

         g1 = wfindex(iband,indtype)+1
         g2 = wfindex(iband+1,indtype)-1
         call dcopy(rvectsize-1,tsrc(2:rvectsize,iblocksize),1,tdest(1,g1:g2),1)
         tdest(1,g1:g2) = factor * tdest(1,g1:g2)

         call dcopy(rvectsize-1,tsrc(rvectsize+1:vectsize,iblocksize),1,tdest(2,g1:g2),1)
         tdest(2,g1:g2) = factor * tdest(2,g1:g2)
       else
         g1 = wfindex(iband,indtype)
         g2 = wfindex(iband+1,indtype)-1
         call dcopy(rvectsize,tsrc(1:rvectsize,iblocksize),1,tdest(1,g1:g2),1)
         tdest(1,g1:g2) = factor * tdest(1,g1:g2)
         call dcopy(rvectsize,tsrc(rvectsize+1:vectsize,iblocksize),1,tdest(2,g1:g2),1)
         tdest(2,g1:g2) = factor * tdest(2,g1:g2)
       end if
     end do
   else
     if (indtype == 'C') then
       g1 = vectsize*((iblock-1)*blocksize)+lig+1
       g2 = vectsize*(iblock*blocksize)+lig
       call zcopy(size,tsrc,1,tdest(:,g1:g2),1)
     else if ( indtype == 'S' ) then
       g1 = vectsize*((iblock-1)*blocksize)+lig+1
       g2 = vectsize*(iblock*blocksize)+lig
       call zcopy(size,tsrc,1,tdest(:,g1:g2),1)
     else
       call zcopy(size,tsrc,1,tdest,1)
     end if
   end if
 else
   MSG_ERROR("Wrong direction: "//trim(direction))
 endif

 if (present(tim_wfcopy).and.present(timopt)) then
   if(abs(timopt)==3) then
     call timab(tim_wfcopy,2,tsec)
   end if
 end if

end subroutine wfcopy
!!***

end module m_wfutils
!!***
