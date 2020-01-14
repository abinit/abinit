!!****m* ABINIT/bond_lotf
!! NAME
!! bond_lotf
!!
!! FUNCTION
!!  Define BOND variables and the procedure to
!!  set them.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2019 ABINIT group (MMancini)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module bond_lotf

 use defs_basis
 use m_errors
 use m_abicore

 implicit none

 public

 !--Fittedatoms variables
 integer :: nfit !--dimension needed for fits
 integer :: nfitmax  !--dimension of the fit arrays
 integer,allocatable :: ifit(:)
 logical,allocatable :: tafit(:)

 !--Fittedbonds variables
 integer :: nbondex
 integer :: ibn_tot,ibn_tot2,ibn_tots
 integer,allocatable :: imat(:)
 integer,dimension(:,:),allocatable :: ibnd_mat,ibmat,ibmat_large

 public ::             &
   bond_atom_init,     &
   bond_tafit_init,    &
   bond_matrix_alloc,  &
   bond_matrix_set,    &
   bond_compute,       &
   bond_fit_set,       &
   bond_dealloc


contains
 !!***


!!****f* bond_lotf/bond_tafit_init
!! NAME
!! bond_tafit_init
!!
!! FUNCTION
!!
!! INPUTS
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine bond_tafit_init(nax)

  implicit none

  !Arguments ------------------------
  integer,intent(in) :: nax

! *************************************************************************

   ABI_ALLOCATE(tafit,(nax))

  !--MMANCINI strange!!!!!
  !  tafit(:nax) = tquant(:nax)
   tafit(:nax) = .true.

 end subroutine bond_tafit_init
 !!***

!!****f* bond_lotf/bond_atom_init
!! NAME
!! bond_atom_init
!!
!! FUNCTION
!!
!! INPUTS
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine bond_atom_init(nneigx,nneig,neighl)

  implicit none

  !Arguments ------------------------
  integer,intent(in) :: nneigx
  integer,intent(in) :: nneig(:)
  integer,intent(in) :: neighl(:,:)

  ! nneigx : max number of neighbours
  ! tau0(3,natom) : atomic positions
  ! neighl(nneigx,natom) : list of neighbours
  ! nneig(natom) : number of neighbours
  ! niter  : iteration number (itime)

  !Local ---------------------------
  integer :: i,j,iat,jat, ibn
  integer,allocatable :: ibnd_dum(:,:)
  character(len=500) :: msg

! *************************************************************************

  !--Initialize nbondex
   nbondex = ((nfitmax/2 * nneigx)+1)/2

  !--Now find initial numbers of active/border bonds :
  !  (0)        CLEARS THE BOND(atom) MATRIX :
   ABI_ALLOCATE(ibnd_mat,(2,nbondex))
   ABI_ALLOCATE(ibnd_dum,(2,nfitmax*6))
   ibnd_mat = 0
   ibnd_dum = 0

   ibn_tot  = 0    !-- bonds between fitted atoms
   ibn_tot2 = 0    !--existing but non optimized bonds with border atoms
   ibn_tots = 0    !--total number of bonds in the fit + border zone

   do i =1,nfit
     iat = ifit(i)
     do j = 1,nneig(iat)
       jat = neighl(j,iat)
       if(tafit(jat)) then
         if(jat > iat) then !--jat is a fitted atom
           ibn_tot = ibn_tot + 1
           ibnd_mat(:,ibn_tot) = (/iat,jat/)
         end if
       else    !--jat is a border atom
         ibn_tot2 = ibn_tot2 + 1
         ibnd_dum(:,ibn_tot2) = (/iat,jat/)
       end if
     end do
   end do

   if(ibn_tot2 > (6*nfitmax)) then
     write(msg,'(3a,i8,2a)')&
&     'ERROR: BOND_ATOM_INIT',ch10,&
&     'IBN_TOT2 =  ',ibn_tot2,ch10,&
&     ' ibnd_dum out of bounds, ibn_tot2 too large '
     MSG_ERROR(msg)

   end if

  !--Reorder to keep 'variational' bonds first :
   ibn_tots = ibn_tot + ibn_tot2
   do ibn=ibn_tot+1,ibn_tots
     ibnd_mat(:,ibn) = ibnd_dum(:,ibn-ibn_tot)
   end do


   ABI_DEALLOCATE(ibnd_dum)
 end subroutine bond_atom_init
 !!***



!!****f* bond_lotf/bond_matrix_alloc
!! NAME
!! bond_matrix_alloc
!!
!! FUNCTION
!!  allocate imat,ibmat,ibmat_large
!! INPUTS
!!
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine bond_matrix_alloc(nax,nneigx)

  implicit none

  !Arguments ------------------------
  integer,intent(in) :: nax
  integer,intent(in) :: nneigx

! *************************************************************************

   ABI_ALLOCATE(ibmat,(nneigx,0:nfitmax))
   ABI_ALLOCATE(ibmat_large,(nfitmax,nfitmax))
   ABI_ALLOCATE(imat,(nax))
 end subroutine bond_matrix_alloc
 !!***

!!****f* bond_lotf/bond_matrix_set
!! NAME
!! bond_matrix_set
!!
!! FUNCTION
!!  Set or update bond matrix imat,ibmat,ibmat_large
!!  associates the bond to the atom neighlists
!! INPUTS
!!
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine bond_matrix_set(nneig,neighl)

  implicit none

  !Arguments ------------------------
  integer,intent(in) :: nneig(:)
  integer,intent(in) :: neighl(:,:)
  !Local ---------------------------
  integer :: ibn,ibn2,iat,jat,ii,jj,j2

! *************************************************************************

   ibmat(:,:) = 0
   ibmat_large(:,:) = 0
   imat(:) = 0

   ibn = 0
   do ii =1,nfit
     iat = ifit(ii)
     do jj = 1,nneig(iat)
       jat = neighl(jj,iat)
       if(jat > iat.AND.tafit(jat)) then
         ibn  = ibn  + 1
         ibmat(jj,ii) = ibn
       end if
     end do
     imat(iat) = ii
   end do

  !--ibmat_large: useful for the glue potential
   ibn2 = 0
   do ii =1,nfit
     iat = ifit(ii)
     do jj = 1,nneig(iat)
       jat = neighl(jj,iat)
       if(jat > iat.AND.tafit(jat)) then
         ibn2  = ibn2  + 1
         j2 = imat(jat)
         ibmat_large(j2,ii) = ibn2
         ibmat_large(ii,j2) = ibn2
       end if
     end do
   end do
 end subroutine bond_matrix_set
 !!***


!!****f* bond_lotf/bond_compute
!! NAME
!! bond_compute
!!
!! FUNCTION
!!  Updates bond matrix, associates the bond to the atom neighlists
!!
!! INPUTS
!!  nneig
!!  neighl
!!
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine bond_compute(nneig,neighl)

  implicit none

  !Arguments ------------------------
  integer,intent(in) :: nneig(:)
  integer,intent(in) :: neighl(:,:)
  !Local ---------------------------
  integer :: ii,jj,iat,jat
  integer, allocatable, dimension(:,:) :: ibnd_dum
  character(len=500) :: msg

! *************************************************************************

   ABI_ALLOCATE(ibnd_dum,(2,nfitmax*6))
   ibnd_dum(:,:) = 0
   ibn_tot  = 0    ! bonds between the fitted atoms
   ibn_tot2 = 0    ! existing but non optimized bonds with border atoms
   ibn_tots = 0    ! total number of bonds


   if(nfit <= 0) then
     write(msg,'(a,i8)')' UPDLIS WARNING : nfit <= 0 = ', nfit
     MSG_WARNING(msg)
   end if

   do ii = 1,nfit
     iat = ifit(ii)
     do jj = 1,nneig(iat)
       jat = neighl(jj,iat)
       if(tafit(jat))then
         if(jat > iat) then  ! jat is a fitted atom
           ibn_tot = ibn_tot + 1
           ibnd_mat(:,ibn_tot) = (/ iat, jat /)
         end if
       else  ! jat is a border atom
         ibn_tot2 = ibn_tot2 + 1
         ibnd_dum(:,ibn_tot2) = (/ iat, jat /)
       end if
     end do
   end do

   if(ibn_tot2 > (6*nfitmax)) then
     write(msg,'(3a,i8,2a)')&
&     'ERROR: BOND_ATOM_INIT',ch10,&
&     'IBN_TOT2 =  ',ibn_tot2,ch10,&
&     ' ibnd_dum out of bounds, ibn_tot2 too large '
     MSG_ERROR(msg)
   end if

  !--reorder to keep 'variational' bonds first :
   ibn_tots = ibn_tot + ibn_tot2
   do ii = ibn_tot+1,ibn_tots
     ibnd_mat(:,ii) = ibnd_dum(:,ii-ibn_tot)
   end do

   ABI_DEALLOCATE (ibnd_dum)
 end subroutine bond_compute
 !!***


!!****f* bond_lotf/bond_fit_set
!! NAME
!! bond_fit_set
!!
!! FUNCTION
!!  set nfitmax (or control it), ifit,nfit
!! INPUTS
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine bond_fit_set(nax,nfitdum)

  implicit none

  !Arguments ------------------------
  integer,intent(in) :: nax
  integer,intent(in) :: nfitdum
  !Local -----------------------------
  integer  :: i
  character(len=500) :: msg

! *************************************************************************

   if(.not. allocated(ifit)) then
     nfitmax = nfitdum + 1000
     ABI_ALLOCATE(ifit,(nfitmax))
   elseif(nfitdum > nfitmax) then
     write(msg,'(a)')' BOND_FIT_SET : PROBLEM OF dimensionS !! '
     MSG_ERROR(msg)
   end if

   ifit = 0
   nfit = 0
   do i = 1, nax
     if(tafit(i)) then
       nfit = nfit + 1
       ifit(nfit) = i
     end if
   end do

 end subroutine bond_fit_set
 !!***

!!****f* bond_lotf/bond_dealloc
!! NAME
!! bond_dealloc
!!
!! FUNCTION
!!  deallocate variables
!!
!! INPUTS
!!
!! PARENTS
!!      m_lotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine  bond_dealloc()

! *************************************************************************
   ABI_DEALLOCATE(ibnd_mat)
   ABI_DEALLOCATE(tafit)
   ABI_DEALLOCATE(ifit)
   ABI_DEALLOCATE(ibmat)
   ABI_DEALLOCATE(ibmat_large)
   ABI_DEALLOCATE(imat)
 end subroutine bond_dealloc

end module bond_lotf
!!***
