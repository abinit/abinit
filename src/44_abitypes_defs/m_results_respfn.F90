!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_results_results
!! NAME
!!  m_results_results
!!
!! FUNCTION
!!  This module provides the definition of the results_respfn_type used
!!  to store results from RESPFN calculations.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (MT,GG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

MODULE m_results_respfn

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_abicore

 implicit none

 private

! public procedures.
 public :: init_results_respfn
 public :: destroy_results_respfn
!!***

!!****t* defs_datatypes/results_respfn_type
!! NAME
!! results_respfn_type
!!
!! FUNCTION
!! This structured datatype contains a subset of the results of a RESPFN
!! calculation, needed to compute the Diagonal Debye-Waller calculation.
!!
!! SOURCE

 type, public :: results_respfn_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar
  integer :: gam_jdtset
  ! jdtset if the results from a q=gamma wavevector calculation, with dataset jdtset have been stored
  ! -jdtset if the results from a q=gamma wavevector calculation, with dataset jdtset should be stored
  ! 0 if no q=gamma wavevector calculation should be stored

  real(dp), allocatable :: gam_eig2nkq(:,:,:,:,:,:,:)  ! one half second derivatives of the electronic eigenvalues
  ! eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom)

 end type results_respfn_type
!!***

CONTAINS

!===========================================================
!!***

!!****f* m_results_out/init_results_respfn
!! NAME
!!  init_results_respfn
!!
!! FUNCTION
!!  Init all scalars and nullify pointers in the structure.
!!
!! INPUTS
!!  ndtset= number of datasets
!!
!! SIDE EFFECTS
!!  results_respfn(:)=<type(results_respfn_type)>=results_respfn datastructure
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_results_respfn(dtsets,ndtset_alloc,results_respfn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndtset_alloc
!arrays
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(results_respfn_type),intent(inout) :: results_respfn
!Local variables-------------------------------
!scalars
 integer :: getgam_eig2nkq,idtset,idtset_gam
 character(len=500) :: message

!************************************************************************

 results_respfn%gam_jdtset=0

!Possibly find the value of gam_jdtset
 do idtset=1,ndtset_alloc
   getgam_eig2nkq=dtsets(idtset)%getgam_eig2nkq
   if(getgam_eig2nkq/=0)then
     do idtset_gam=1,ndtset_alloc
       if(dtsets(idtset_gam)%jdtset==getgam_eig2nkq)exit
     enddo
     if(dtsets(idtset_gam)%optdriver/=RUNL_RESPFN)then
       write(message, '(a,i5,a,i5,2a,i5,3a)' )&
&        'For jdtset=',dtsets(idtset)%jdtset,', getgam_eig2nkq=',getgam_eig2nkq,ch10,&
&        'However this dataset with idtset=',idtset_gam, ' is not a phonon calculation;',ch10,&
&        'Action : correct the value of getgam_eig2nkq for that dataset.'
       MSG_ERROR(message)
     endif
     if(results_respfn%gam_jdtset==0)then
       results_respfn%gam_jdtset=-getgam_eig2nkq ! Store a negative value, indicating that it is expected.
     else if(results_respfn%gam_jdtset/=-getgam_eig2nkq) then
       write(message, '(a,i5,2a,i5,2a)' )&
&        'results_respfn%gam_jdtset=',results_respfn%gam_jdtset,ch10,&
&        'dtsets(idtset)%getgam_eig2nkq=',getgam_eig2nkq,ch10,&
&        'So, it seems that two gamma q point calculations should be stored, while this is not yet allowed.'
       MSG_BUG(message)
     endif
   endif
 enddo

end subroutine init_results_respfn
!!***

!----------------------------------------------------------------------

!!****f* m_results_out/destroy_results_respfn
!! NAME
!!  destroy_results_respfn
!!
!! FUNCTION
!!  Clean and destroy results_respfn datastructure
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  results_respfn(:)=<type(results_respfn_type)>=results_respfn datastructure
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_results_respfn(results_respfn)

!Arguments ------------------------------------
!arrays
 type(results_respfn_type),intent(inout) :: results_respfn

!************************************************************************

 if(results_respfn%gam_jdtset>0)then
!DEBUG
!  write(std_out,*)' results_respfn%gam_jdtset=',results_respfn%gam_jdtset
!  write(std_out,*)' allocated(results_respfn%gam_eig2nkq)=',allocated(results_respfn%gam_eig2nkq)
!  call flush(6)
!ENDDEBUG
   if(allocated(results_respfn%gam_eig2nkq))then
     ABI_DEALLOCATE(results_respfn%gam_eig2nkq)
   endif
 endif

end subroutine destroy_results_respfn

!----------------------------------------------------------------------

END MODULE m_results_respfn
!!***
