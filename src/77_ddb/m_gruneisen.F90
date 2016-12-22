!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_gruneisen
!! NAME
!!  m_gruneisen
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  to compute Gruneisen parameters.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2016 ABINIT group (MG)
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

MODULE m_gruneisen

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_ddb
 use m_ifc
 use m_cgtools
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_io_tools,        only : get_unit
 use m_fstrings,        only : sjoin, itoa
 use m_crystal,         only : crystal_t
 use m_kpts,            only : kpts_ibz_from_kptrlatt
 use m_bz_mesh,         only : kpath_t
 use m_anaddb_dataset,  only : anaddb_dataset_type

 implicit none

 private
!!***

!!****t* m_gruneisen/gruns_t
!! NAME
!! gruns_t
!!
!! FUNCTION
!!  Contains the necessary data to interpolate the
!!  phonon bandstructure and eigenvectors in reciprocal space (ie.
!!  interatomic force constants and corresponding real space grid info).
!!
!! SOURCE

 type,public :: gruns_t

   integer :: natom3
    ! 3 * natom

   integer :: nvols
    ! Number of volumes.

   integer :: iv0
    ! Index of the DDB file associated to the equilibrium volume

   real(dp) :: v0
    ! Equilibrium volume.

   real(dp) :: delta_vol
    ! Uniform grid spacing for finite difference

   type(crystal_t),allocatable :: cryst_vol(:)
    ! cryst_vol(nvols)
    ! crystalline structure for the different volumes

   type(ddb_type),allocatable :: ddb_vol(:)
    ! dbb_vol(nvols)
    ! DDB objects for the different volumes

   type(ifc_type),allocatable :: ifc_vol(:)
    ! ifc_vol(nvols)
    ! interatomic force constants for the different volumes

 end type gruns_t

 public :: gruns_new        ! Constructor
 public :: gruns_qpath      ! Compute Grunesein parameters on a q-path
 public :: gruns_qmesh      ! Compute Grunesein parameters on a q-mesh.
 public :: gruns_free       ! Release memory
!!***

contains  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_new
!! NAME
!!  gruns_new
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(gruns_t) function gruns_new(ddb_paths, inp, comm) result(new)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_free'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: comm
 type(anaddb_dataset_type),intent(in) :: inp
!arrays
 character(len=*),intent(in) :: ddb_paths(:)

!Local variables-------------------------------
 integer,parameter :: natifc0=0
 integer :: ivol,iblock,msym,dimekb,lmnmax,mband,mblktyp,natom,nblok,nkpt,ntypat,usepaw,ddbun
!arrays
 integer,allocatable :: atifc0(:)
 real(dp) :: dielt(3,3)
 real(dp),allocatable :: zeff(:,:,:)

! ************************************************************************

 new%nvols = len(ddb_paths)
 ABI_MALLOC(new%cryst_vol, (new%nvols))
 ABI_MALLOC(new%ddb_vol, (new%nvols))
 ABI_MALLOC(new%ifc_vol, (new%nvols))

 ddbun = get_unit()
 do ivol=1,new%nvols
   call ddb_getdims(dimekb,ddb_paths(ivol),lmnmax,mband,mblktyp,msym,natom,nblok,nkpt,ntypat,ddbun,usepaw,DDB_VERSION,comm)
   ABI_MALLOC(atifc0, (natom))
   atifc0 = 0
   call ddb_from_file(new%ddb_vol(ivol), ddb_paths(ivol), inp%brav, natom, natifc0, atifc0, new%cryst_vol(ivol), comm)
   ABI_FREE(atifc0)

   ! Get Dielectric Tensor and Effective Charges
   ! (initialized to one_3D and zero if the derivatives are not available in the DDB file)
   ABI_MALLOC(zeff, (3,3,natom))
   iblock = ddb_get_dielt_zeff(new%ddb_vol(ivol), new%cryst_vol(ivol), inp%rfmeth, inp%chneut, inp%selectz, dielt, zeff)
   if (iblock == 0) then
     call wrtout(ab_out, sjoin("- Cannot find dielectric tensor and Born effective charges in DDB file:", ddb_paths(ivol)))
     call wrtout(ab_out, "Values initialized with zeros")
   else
     call wrtout(ab_out, sjoin("- Found dielectric tensor and Born effective charges in DDB file:", ddb_paths(ivol)))
   end if

   call ifc_init(new%ifc_vol(ivol), new%cryst_vol(ivol), new%ddb_vol(ivol),&
     inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,inp%ngqpt(1:3),inp%nqshft,inp%q1shft,dielt,zeff,&
     inp%nsphere,inp%rifcsph,inp%prtsrlr,inp%enunit)
   ABI_FREE(zeff)
 end do

 new%natom3 = new%cryst_vol(1)%natom * 3
 new%iv0 = 2
 new%v0 = new%cryst_vol(new%iv0)%ucvol
 new%delta_vol = new%cryst_vol(new%iv0+1)%ucvol - new%cryst_vol(new%iv0)%ucvol
 ! TODO
 ! Consistency check

end function gruns_new
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_fourq
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!  qpt(3)=q-point in reduced coordinates.
!!
!! OUTPUT
!!  phfrq(3*natom) = Gruneisen parameters
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gruns_fourq(gruns, qpt, wv0, gvals)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(gruns_t),intent(in) :: gruns
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: wv0(gruns%natom3),gvals(gruns%natom3)

!Local variables-------------------------------
!scalars
 integer :: ivol,ii,natom3,nu
!arrays
 real(dp) :: wvols(gruns%natom3,gruns%nvols),dot(2)
 real(dp) :: eigvec(2,gruns%natom3,gruns%natom3,gruns%nvols),d2cart(2,gruns%natom3,gruns%natom3,gruns%nvols)
 real(dp) :: dddv(2,gruns%natom3,gruns%natom3),displ_cart(2,gruns%natom3)
 real(dp) :: omat(2,gruns%natom3, gruns%natom3)

! ************************************************************************

 natom3 = gruns%natom3
 do ivol=1,gruns%nvols
   call ifc_fourq(gruns%ifc_vol(ivol), gruns%cryst_vol(ivol), qpt, wvols(:,ivol), displ_cart, &
                  out_d2cart=d2cart(:,:,:,ivol), out_eigvec=eigvec(:,:,:,ivol))
   !call massmult_and_breaksym(gruns%cryst_vol(ivol)%natom, gruns%cryst_vol(ivol)%ntypat,
   !  gruns%cryst_vol(ivol)%cryst%typat, gruns%ifc_vol(ivol)%amu, d2cart(:,:,:,ivol))
 end do
 wv0 = wvols(:, gruns%iv0)

 ! Compute dD(q)/dV with finite difference.
 select case (gruns%nvols)
 case (3)
   dddv = -half*d2cart(:,:,:,1) + half*d2cart(:,:,:,3)
 case default
   MSG_ERROR(sjoin("Finite difference with nvols", itoa(gruns%nvols), "not implemented"))
 end select
 dddv = dddv / gruns%delta_vol

 ! Compute -volume/(2w(q)**2) <u(q)|dD(q)/dq|u(q)>
 call zgemm('N','N',natom3,natom3,natom3,cone,dddv,natom3,eigvec(:,:,:,gruns%iv0),natom3,czero,omat,natom3)
 do nu=1,natom3
   if (abs(wv0(nu)) > tol12) then
     dot = cg_zdotc(natom3, eigvec(1,1,nu, gruns%iv0), omat(1,1,nu))
     gvals(nu) = -gruns%v0 * dot(1) / (two * wv0(nu)**2)
     !write(std_out,*)"dot", dot
   else
     gvals(nu) = zero
   end if
 end do

end subroutine gruns_fourq
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_qpath
!! NAME
!!
!! FUNCTION
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

subroutine gruns_qpath(gruns, qpath, comm)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(gruns_t),intent(in) :: gruns
 type(kpath_t),intent(in) :: qpath

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: nprocs,my_rank,iqpt,ierr
!arrays
 !real(dp) :: qpt(3)
 real(dp),allocatable :: gvals(:,:),wv0(:,:)

! ************************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 ABI_CALLOC(wv0, (gruns%natom3, qpath%npts))
 ABI_CALLOC(gvals, (gruns%natom3, qpath%npts))

 do iqpt=1,qpath%npts
   if (mod(iqpt, nprocs) /= my_rank) cycle ! mpi-parallelism
   call gruns_fourq(gruns, qpath%points(:,iqpt), wv0(:,iqpt), gvals(:,iqpt))
 end do

 call xmpi_sum(wv0, comm, ierr)
 call xmpi_sum(gvals, comm, ierr)
 ABI_FREE(wv0)
 ABI_FREE(gvals)

end subroutine gruns_qpath
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_qmesh
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!  ngqpt(3)=q-mesh divisions
!!  nqshift=Number of shifts used to generated the ab-initio q-mesh.
!!  qshift(3,nqshift)=The shifts of the ab-initio q-mesh.
!!  comm=MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gruns_qmesh(gruns, ngqpt, nqshift, qshift, comm)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqshift,comm
 type(gruns_t),intent(in) :: gruns
!arrays
 integer,intent(in) :: ngqpt(3)
 real(dp),intent(in) :: qshift(3,nqshift)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,qptopt1=1
 integer :: nprocs,my_rank,iqibz,nqbz,nqibz,ierr,ii
!arrays
 real(dp),allocatable :: gvals(:,:),wv0(:,:)
 integer :: qptrlatt(3,3)
 real(dp),allocatable :: qibz(:,:),qbz(:,:),wtq(:)

! ************************************************************************

 ! Generate the q-mesh by finding the IBZ and the corresponding weights.
 qptrlatt = 0
 do ii=1,3
   qptrlatt(ii,ii) = ngqpt(ii)
 end do

 ! Get IBZ and BZ.
 call kpts_ibz_from_kptrlatt(gruns%cryst_vol(gruns%iv0), qptrlatt, qptopt1, nqshift, qshift, &
   nqibz, qibz, wtq, nqbz, qbz)

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 ABI_CALLOC(wv0, (gruns%natom3, nqibz))
 ABI_CALLOC(gvals, (gruns%natom3, nqibz))

 do iqibz=1,nqibz
   if (mod(iqibz, nprocs) /= my_rank) cycle ! mpi-parallelism
   call gruns_fourq(gruns, qibz(:,iqibz), wv0(:,iqibz), gvals(:,iqibz))
 end do

 call xmpi_sum(wv0, comm, ierr)
 call xmpi_sum(gvals, comm, ierr)

 ABI_FREE(qibz)
 ABI_FREE(wtq)
 ABI_FREE(qbz)
 ABI_FREE(wv0)
 ABI_FREE(gvals)

end subroutine gruns_qmesh
!!***

!----------------------------------------------------------------------

!!****f* m_gruneisen/gruns_free
!! NAME
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine gruns_free(gruns)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ifc_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(gruns_t),intent(inout) :: gruns

!Local variables-------------------------------
!scalars
 integer :: ii

! ************************************************************************

 if (allocated(gruns%ifc_vol)) then
   do ii=1,size(gruns%cryst_vol)
    call crystal_free(gruns%cryst_vol(ii))
   end do
   ABI_FREE(gruns%cryst_vol)
 end if

 if (allocated(gruns%ddb_vol)) then
   do ii=1,size(gruns%ddb_vol)
    call crystal_free(gruns%ddb_vol(ii))
   end do
   ABI_FREE(gruns%ddb_vol)
 end if

 if (allocated(gruns%ifc_vol)) then
   do ii=1,size(gruns%ifc_vol)
    call ifc_free(gruns%ifc_vol(ii))
   end do
   ABI_FREE(gruns%ifc_vol)
 end if

end subroutine gruns_free
!!***

!----------------------------------------------------------------------

end module m_gruneisen
