!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ephwg
!! NAME
!! m_ephwg
!!
!! FUNCTION
!!  Tools and objects to compute the weights used for the BZ integration of EPH quantities.
!!  More specifically the integration of quantities such as the imaginary part of the self-energy
!!  involving delta functions. Different approaches are available:
!!
!!    1.
!!    2.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_ephwg

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_copy
 use m_tetrahedron
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_crystal
 use m_ifc
 use m_lgroup
 use m_ebands
 use m_eph_double_grid

 use defs_datatypes,    only : ebands_t
 use defs_abitypes,     only : dataset_type
 use m_time,            only : cwtime, sec2str
 use m_symtk,           only : matr3inv
 use m_numeric_tools,   only : arth, inrange, wrap2_pmhalf
 use m_special_funcs,   only : dirac_delta
 use m_fstrings,        only : strcat, ltoa, itoa, ftoa, ktoa, sjoin
 use m_simtet,          only : sim0onei, SIM0TWOI
 use m_kpts,            only : kpts_timrev_from_kptopt, listkk, kpts_ibz_from_kptrlatt
 use m_occ,             only : occ_fd, occ_be

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_ephwg/ephwg_t
!! NAME
!! ephwg_t
!!
!! FUNCTION
!!  Stores electron eigevalues and phonon frequencies in the IBZ (assume same mesh for e and ph).
!!  Provides tools to compute (e_{k+q} - w{q}) in the IBZ(k)
!!  and integrate the delta functions for phonon emission/absorption with the tetrahedron method.
!!
!! SOURCE

type, public :: ephwg_t

  integer :: natom3
  ! 3 * natom

  integer :: nsppol
  ! Number of independent spin polarizations.

  integer :: nbcount
  ! Number of bands treated.

  integer :: bstart
  ! The fist band (global index) starts at bstart.
  ! Used to select bands around the Fermi level.

  integer :: kptopt
  ! Option for k-point generation.

  integer :: timrev
  ! 1 if the use of time-reversal is allowed; 0 otherwise

  integer :: nibz, nbz
  ! Number of q-points in IBZ and full BZ.

  integer :: nq_k
  ! Nunber of points in IBZ(k) i.e. the irreducible wedge
  ! defined by the operations of the little group of k.

  integer, allocatable :: kq2ibz(:)
  ! kq2ibz(nq_k)
  ! Mapping (k + q) --> initial IBZ array

  real(dp),allocatable :: ibz(:,:)
  ! ibz(3, nibz)
  ! The initial IBZ.

  real(dp),allocatable :: bz(:,:)
  ! bz(3, nbz)
  ! points in full BZ.

  real(dp) :: klatt(3, 3)
  ! Reciprocal of lattice vectors for full kpoint grid. Used by init_tetra

  integer,allocatable :: lgk2ibz(:)
  ! lgk2ibz(nq_k)
  ! Mapping Little-group IBZ_k --> initial IBZ

  real(dp),allocatable :: phfrq_ibz(:,:)
  ! (nibz, natom3)
  ! Phonon frequencies in the IBZ

  real(dp),allocatable :: eigkbs_ibz(:, :, :)
  ! (nibz, nbcount, nsppol)
  ! Electron eigenvalues in the IBZ for nbcount states
  ! (not necessarly equal to global nband, see also bstart and bcount)

  type(crystal_t), pointer :: cryst => null()
  ! Pointer to input structure (does not own memory)

  type(lgroup_t) :: lgk
  ! Little group of the k-point

  type(t_tetrahedron) :: tetra_k
  ! Used to evaluate delta(w - e_{k+q} +/- phw_q) with tetrahedron method.

 end type ephwg_t

 public :: ephwg_new             ! Basic Constructor
 public :: ephwg_from_ebands     ! Build object from ebands_t
 public :: ephwg_setup_kpoint    ! Prepare tetrahedron method for given external k-point.
 public :: ephwg_double_grid_setup_kpoint ! Prepare tetrahedron method for given external k-point using double grid routines.
 public :: ephwg_get_deltas      ! Compute weights for $ \int \delta(\omega - \ee_{k+q, b} \pm \omega_{q\nu} $
 public :: ephwg_get_dweights    ! Compute weights for a set of frequencies.
 public :: ephwg_zinv_weights    ! Compute weights for $ \int 1 / (\omega - \ee_{k+q, b} \pm \omega_{q\nu} $
 public :: ephwg_free            ! Free memory
 !public :: ephwg_test

contains
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_new
!! NAME
!! ephwg_new
!!
!! FUNCTION
!!  Initialize the object from the electronic eigenvalues given in the IBZ.
!!
!! INPUTS
!!  cryst<cryst_t>=Crystalline structure.
!!  ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!!  bstart=Index of the first band to be included.
!!  nbcount=Number of bands included
!!  kptopt=Option for the k-point generation.
!!  kptrlatt(3,3)=k-point lattice specification
!!  nshiftk= number of shift vectors.
!!  shiftk(3,nshiftk)=shift vectors for k point generation
!!  nkibz=Number of points in the IBZ
!!  kibz(3,nkibz)=Reduced coordinates of the k-points in the IBZ.
!!  nsppol=Number of independent spin polarizations.
!!  eig_ibz(nbcount, nkibz, nsppol) = Electron eigenvalues for nbcount states
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(ephwg_t) function ephwg_new( &
&  cryst, ifc, bstart, nbcount, kptopt, kptrlatt, nshiftk, shiftk, nkibz, kibz, nsppol, eig_ibz, comm) result(new)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt, nshiftk, nkibz, bstart, nbcount, nsppol, comm
 type(crystal_t),target,intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3, nshiftk), kibz(3, nkibz)
 real(dp),intent(in) :: eig_ibz(nbcount, nkibz, nsppol)

!Local variables-------------------------------
!scalars
 integer :: nprocs, my_rank, ik, ierr, out_nkibz
!arrays
 real(dp) :: rlatt(3,3)
 integer :: out_kptrlatt(3,3)
 real(dp) :: displ_cart(2,3,cryst%natom,3*cryst%natom), phfrq(3*cryst%natom)
 real(dp),allocatable :: out_kibz(:,:), out_wtk(:)

!----------------------------------------------------------------------

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 new%natom3 = ifc%natom * 3
 new%nsppol = nsppol
 new%nbcount = nbcount
 new%bstart = bstart
 new%kptopt = kptopt
 new%timrev = kpts_timrev_from_kptopt(new%kptopt)
 new%nibz = nkibz
 new%cryst => cryst
 call alloc_copy(kibz, new%ibz)

 ! Get full BZ (new%nbz, new%bz) and new kptrlatt for tetra.
 call kpts_ibz_from_kptrlatt(cryst, kptrlatt, kptopt, nshiftk, shiftk, out_nkibz, out_kibz, out_wtk, new%nbz, new%bz, &
   new_kptrlatt=out_kptrlatt)

 rlatt = out_kptrlatt; call matr3inv(rlatt, new%klatt)

 ABI_CHECK(size(out_kibz, dim=2) == new%nibz, "mismatch in nkibz!")
 ABI_FREE(out_kibz)
 ABI_FREE(out_wtk)

 ! Copy eigenvalues in IBZ. Change shape for better performance in other routines.
 ABI_MALLOC(new%eigkbs_ibz, (new%nibz, new%nbcount, new%nsppol))
 do ik=1,new%nibz
   new%eigkbs_ibz(ik, :, :) = eig_ibz(:, ik, :)
 end do

 ! Fourier interpolate frequencies on the same mesh.
 ABI_CALLOC(new%phfrq_ibz, (new%nibz, new%natom3))
 do ik=1,new%nibz
   if (mod(ik, nprocs) /= my_rank) cycle ! mpi-parallelism
   call ifc_fourq(ifc, cryst, new%ibz(:, ik), phfrq, displ_cart)
   new%phfrq_ibz(ik, :) = phfrq
 end do

 call xmpi_sum(new%phfrq_ibz, comm, ierr)

end function ephwg_new
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_from_ebands
!! NAME
!! ephwg_from_ebands
!!
!! FUNCTION
!!  Convenience constructor to initialize the object from an ebands_t object

type(ephwg_t) function ephwg_from_ebands(cryst, ifc, ebands, bstart, nbcount, comm) result(new)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) ::  bstart, nbcount, comm
 type(crystal_t),intent(in) :: cryst
 type(ifc_type),intent(in) :: ifc
 type(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
 real(dp),allocatable :: eig_ibz(:, :, :)

!----------------------------------------------------------------------

 if (bstart == 1 .and. nbcount == ebands%mband) then
   new = ephwg_new(cryst, ifc, bstart, nbcount, ebands%kptopt, ebands%kptrlatt, ebands%nshiftk, ebands%shiftk, ebands%nkpt, &
      ebands%kptns, ebands%nsppol, ebands%eig, comm)
 else
   ABI_CHECK(inrange(bstart, [1, ebands%mband]), "Wrong bstart")
   ABI_CHECK(inrange(bstart + nbcount - 1, [1, ebands%mband]), "Wrong nbcount")
   ! Copy submatrix of eigenvalues
   ABI_MALLOC(eig_ibz, (nbcount, ebands%nkpt, ebands%nsppol))
   eig_ibz = ebands%eig(bstart:bstart+nbcount-1, : , :)
   new = ephwg_new(cryst, ifc, bstart, nbcount, ebands%kptopt, ebands%kptrlatt, ebands%nshiftk, ebands%shiftk, ebands%nkpt, &
      ebands%kptns, ebands%nsppol, eig_ibz, comm)
   ABI_FREE(eig_ibz)
 end if

end function ephwg_from_ebands
!!***

!!****f* m_ephwg/ephwg_setup_kpoint
!! NAME
!! ephwg_setup_kpoint
!!
!! FUNCTION
!!  Set internal tables and object required to compute integration weights for a given k-point.
!!
!! INPUTS
!!  kpoint(3): k-point in reduced coordinates.
!!  prtvol: Verbosity level
!!  comm: MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_setup_kpoint(self, kpoint, prtvol, comm)

 implicit none

!Arguments ------------------------------------
!scalars
 type(ephwg_t),target,intent(inout) :: self
 integer,intent(in) :: prtvol, comm
!arrays
 real(dp),intent(in) :: kpoint(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1=1
 integer :: ierr,ii
 real(dp) :: dksqmax
 character(len=80) :: errorstring
 character(len=500) :: msg
 type(crystal_t),pointer :: cryst
!arrays
 integer,allocatable :: indkk(:,:)

!----------------------------------------------------------------------

 cryst => self%cryst

 ! Get little group of the kpoint.
 call lgroup_free(self%lgk)
 self%lgk = lgroup_new(self%cryst, kpoint, self%timrev, self%nbz, self%bz, self%nibz, self%ibz)
 if (prtvol > 0) call lgroup_print(self%lgk)
 self%nq_k = self%lgk%nibz

 ! Get mapping IBZ_k --> initial IBZ (self%lgrp%ibz --> self%ibz)
 ABI_MALLOC(indkk, (self%nq_k * sppoldbl1, 6))
 call listkk(dksqmax, cryst%gmet, indkk, self%ibz, self%lgk%ibz, self%nibz, self%nq_k, cryst%nsym,&
    sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, comm, use_symrec=.False.)

 if (dksqmax > tol12) then
   write(msg, '(a,es16.6)' ) &
    "At least one of the points in IBZ(k) could not be generated from a symmetrical one. dksqmax: ",dksqmax
   MSG_ERROR(msg)
 end if
 if (allocated(self%lgk2ibz)) then
   ABI_FREE(self%lgk2ibz)
 end if
 call alloc_copy(indkk(:, 1), self%lgk2ibz)
 ABI_FREE(indkk)

 ! Get mapping (k + q) --> initial IBZ.
 do ii=1,self%nq_k
   self%lgk%ibz(:, ii) = self%lgk%ibz(:, ii) + kpoint
 end do
 ABI_MALLOC(indkk, (self%nq_k * sppoldbl1, 6))

 call listkk(dksqmax, cryst%gmet, indkk, self%ibz, self%lgk%ibz, self%nibz, self%nq_k, cryst%nsym,&
    sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, comm, use_symrec=.False.)

 if (dksqmax > tol12) then
   write(msg, '(a,es16.6)' ) &
    "At least one of the points in IBZ(k) + q could not be generated from a symmetrical one. dksqmax: ",dksqmax
   MSG_ERROR(msg)
 end if

 if (allocated(self%kq2ibz)) then
   ABI_FREE(self%kq2ibz)
 end if
 call alloc_copy(indkk(:, 1), self%kq2ibz)
 ABI_FREE(indkk)
 do ii=1,self%nq_k
   self%lgk%ibz(:, ii) = self%lgk%ibz(:, ii) - kpoint
 end do

 ! Get mapping BZ --> IBZ_k (self%bz --> self%lgrp%ibz) required for tetrahedron method
 ABI_MALLOC(indkk, (self%nbz * sppoldbl1, 6))
 call listkk(dksqmax, cryst%gmet, indkk, self%lgk%ibz, self%bz, self%nq_k, self%nbz, cryst%nsym,&
    sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, comm, use_symrec=.False.)

 ! Build tetrahedron object using IBZ(k) as the effective IBZ
 ! This means that input data for tetra routines must be provided in lgk%kibz_q
 call destroy_tetra(self%tetra_k)
 call init_tetra(indkk(:, 1), cryst%gprimd, self%klatt, self%bz, self%nbz, self%tetra_k, ierr, errorstring)
 ABI_CHECK(ierr == 0, errorstring)
 ABI_FREE(indkk)

end subroutine ephwg_setup_kpoint
!!***

!!****f* m_ephwg/ephwg_double_grid_setup_kpoint
!! NAME
!! ephwg_setup_kpoint
!!
!! FUNCTION
!!  Set internal tables and object required to compute integration weights for a given k-point
!!  using the double grid routines to map the different k-points.
!!  This version should be more efficient than its counterpart without the double grid.
!!
!! INPUTS
!!  kpoint(3): k-point in reduced coordinates.
!!  prtvol: Verbosity level
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_double_grid_setup_kpoint(self, eph_doublegrid, kpoint, prtvol)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ephwg_setup_kpoint'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ephwg_t),target,intent(inout) :: self
 type(eph_double_grid_t),intent(inout) :: eph_doublegrid
 integer,intent(in) :: prtvol
!arrays
 real(dp),intent(in) :: kpoint(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: sppoldbl1=1
 integer :: ierr,ii,ik_idx,ik_bz,isym
 real(dp) :: dksqmax
 character(len=80) :: errorstring
 character(len=500) :: msg
 type(crystal_t),pointer :: cryst
!arrays
 integer,allocatable :: indkk(:,:), lgkibz2bz(:), lgkibz2ibz(:), lgkibz2ibzkq(:)
 integer,allocatable :: bz2lgkibz(:), bz2lgkibzkq(:), bz2bz(:), mapping(:,:)
 real(dp) :: kpt_sym(3), kpt(3), wrap_kpt(3), shift
!----------------------------------------------------------------------

 cryst => self%cryst

 ! Get little group of the kpoint.
 call lgroup_free(self%lgk)
 self%lgk = lgroup_new(self%cryst, kpoint, self%timrev, self%nbz, self%bz, self%nibz, self%ibz)
 if (prtvol > 0) call lgroup_print(self%lgk)
 self%nq_k = self%lgk%nibz

 ! get dg%bz --> self%lgrp%ibz
 ABI_SFREE(eph_doublegrid%bz2lgkibz)
 ABI_MALLOC(eph_doublegrid%bz2lgkibz,(eph_doublegrid%dense_nbz))
 call eph_double_grid_bz2ibz(eph_doublegrid, self%lgk%ibz, self%lgk%nibz,&
                             cryst%symrel, cryst%nsym, &
                             eph_doublegrid%bz2lgkibz, has_timrev=1)

#if 0
   ABI_MALLOC(indkk, (eph_doublegrid%dense_nbz * sppoldbl1, 6))
   call listkk(dksqmax, cryst%gmet, indkk, self%lgk%ibz, eph_doublegrid%kpts_dense,& 
               self%lgk%nibz, eph_doublegrid%dense_nbz, cryst%nsym,&
               sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, use_symrec=.False.)

   if (dksqmax > tol12) then
     write(msg, '(a,es16.6)' ) &
      "At least one of the points in IBZ(k) could not be generated from a symmetrical one. dksqmax: ",dksqmax
     MSG_ERROR(msg)
   end if

   !check if same results as listkk
   !eph_doublegrid%bz2lgkibz(:)=indkk(:,1)
   do ii=1,self%nbz
      if (indkk(ii,1).ne.eph_doublegrid%bz2lgkibz(ii)) then
        write(*,*) ii
        write(*,*) eph_doublegrid%kpts_dense(:,ii)
        write(*,*) ii, '-->', indkk(ii,1)
        write(*,*) self%lgk%ibz(:,indkk(ii,1)), indkk(ii,2), indkk(ii,6)
        !check listkk stuff
        kpt = (1-2*indkk(ii,6))*matmul(transpose(cryst%symrel(:,:,indkk(ii,2))),self%lgk%ibz(:,indkk(ii,1)))
        call wrap2_pmhalf(kpt(1),wrap_kpt(1),shift)
        call wrap2_pmhalf(kpt(2),wrap_kpt(2),shift)
        call wrap2_pmhalf(kpt(3),wrap_kpt(3),shift)
        write(*,*) wrap_kpt

        write(*,*) ii, '-->', eph_doublegrid%bz2lgkibz(ii)
        write(*,*) self%lgk%ibz(:,eph_doublegrid%bz2lgkibz(ii)), mapping(ii,1), mapping(ii,2)
        !check double grid stuff
        kpt = (1-2*mapping(ii,2))*matmul(transpose(cryst%symrel(:,:,mapping(ii,1))),self%lgk%ibz(:,eph_doublegrid%bz2lgkibz(ii)))
        call wrap2_pmhalf(kpt(1),wrap_kpt(1),shift)
        call wrap2_pmhalf(kpt(2),wrap_kpt(2),shift)
        call wrap2_pmhalf(kpt(3),wrap_kpt(3),shift)
        write(*,*) wrap_kpt
        write(*,*)
      end if
      ABI_CHECK(indkk(ii,1)==eph_doublegrid%bz2lgkibz(ii),'Unmatching indexes')
   enddo
   ABI_FREE(indkk)
#endif
 
 ! self%lgrp%ibz --> dg%bz
 ABI_CALLOC(lgkibz2bz,(self%lgk%nibz))
 do ii=1,self%nbz
   ik_idx = eph_doublegrid%bz2lgkibz(ii)
   lgkibz2bz(ik_idx) = ii
 enddo

 ! get self%lgrp%ibz --> dg%bz --> self%ibz
 ABI_SFREE(self%lgk2ibz)
 ABI_MALLOC(self%lgk2ibz, (self%nq_k))
 do ii=1,self%nq_k
   ik_idx = lgkibz2bz(ii)
   self%lgk2ibz(ii) = eph_doublegrid%bz2ibz_dense(ik_idx)
 enddo

#if 0
   ! Get mapping IBZ_k --> initial IBZ (self%lgrp%ibz --> self%ibz)
   ABI_MALLOC(indkk, (self%nq_k * sppoldbl1, 6))
   call listkk(dksqmax, cryst%gmet, indkk, self%ibz, self%lgk%ibz, self%nibz, self%nq_k, cryst%nsym,&
      sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, use_symrec=.False.)

   if (dksqmax > tol12) then
     write(msg, '(a,es16.6)' ) &
      "At least one of the points in IBZ(k) could not be generated from a symmetrical one. dksqmax: ",dksqmax
     MSG_ERROR(msg)
   end if

   !check if same results as listkk
   do ii=1,self%nq_k
      ABI_CHECK(self%lgk2ibz(ii)==indkk(ii,1),'Unmatching indexes')
   enddo
   ABI_FREE(indkk)
#endif

 ! calculate k+q
 do ii=1,self%nq_k
   self%lgk%ibz(:, ii) = self%lgk%ibz(:, ii) + kpoint
 end do

 ! get dg%bz --> lgk%ibz (k+q)
 ABI_MALLOC(bz2lgkibzkq, (eph_doublegrid%dense_nbz))
 call eph_double_grid_bz2ibz(eph_doublegrid, self%lgk%ibz, self%lgk%nibz,&
                             cryst%symrel, cryst%nsym, &
                             bz2lgkibzkq, has_timrev=1)

 ! self%lgrp%ibz (k+q) --> dg%bz
 do ii=1,self%nbz
   ik_idx = bz2lgkibzkq(ii)
   lgkibz2bz(ik_idx) = ii
 enddo
 ABI_FREE(bz2lgkibzkq)

 ! get self%lgrp%ibz (k+q) --> dg%bz --> self%ibz
 ABI_SFREE(self%kq2ibz)
 ABI_MALLOC(self%kq2ibz, (self%nq_k))
 do ii=1,self%nq_k
   ik_idx = lgkibz2bz(ii)
   self%kq2ibz(ii) = eph_doublegrid%bz2ibz_dense(ik_idx)
 enddo
 ABI_FREE(lgkibz2bz)

#if 0
   ! Get mapping (k + q) --> initial IBZ.
   ABI_MALLOC(indkk, (self%nq_k * sppoldbl1, 6))
   call listkk(dksqmax, cryst%gmet, indkk, self%ibz, self%lgk%ibz, self%nibz, self%nq_k, cryst%nsym,&
      sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, use_symrec=.False.)

   if (dksqmax > tol12) then
     write(msg, '(a,es16.6)' ) &
      "At least one of the points in IBZ(k) + q could not be generated from a symmetrical one. dksqmax: ",dksqmax
     MSG_ERROR(msg)
   end if

   !check if same results as listkk
   do ii=1,self%nq_k
      ABI_CHECK(self%kq2ibz(ii)==indkk(ii,1),'Unmatching indexes')
   enddo
   ABI_FREE(indkk)
#endif

 ! revert change
 do ii=1,self%nq_k
   self%lgk%ibz(:, ii) = self%lgk%ibz(:, ii) - kpoint
 end do

 ! get self%bz --> dg%bz --> self%lgrp%ibz
 ABI_MALLOC(bz2lgkibz,(self%nbz))

#if 0
   ABI_MALLOC(indkk, (self%nbz * sppoldbl1, 6))
   ABI_MALLOC(bz2bz, (self%nbz * sppoldbl1))
   call listkk(dksqmax, cryst%gmet, indkk, eph_doublegrid%kpts_dense, self%bz,& 
      eph_doublegrid%dense_nbz, self%nbz, cryst%nsym,&
      sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, use_symrec=.False.)

   do ii=1,self%nbz
      ! get self%bz --> self%bz
      bz2bz(ii) = eph_double_grid_get_index(eph_doublegrid,self%bz(:,ii),2)
   end do

   !check if same results as listkk
   do ii=1,self%nbz
      ABI_CHECK(indkk(ii,1)==bz2bz(ii),'Unmatching indexes')
   enddo
   ABI_FREE(indkk)
   ABI_FREE(bz2bz)
#endif

 do ii=1,self%nbz
    ! get self%bz --> dg%bz
    ik_idx = eph_double_grid_get_index(eph_doublegrid,self%bz(:,ii),2)
    ! dg%bz --> self%lgrp%ibz
    bz2lgkibz(ii) = eph_doublegrid%bz2lgkibz(ik_idx)
 end do

#if 0
   ! Get mapping BZ --> IBZ_k (self%bz --> self%lgrp%ibz) required for tetrahedron method
   ABI_MALLOC(indkk, (self%nbz * sppoldbl1, 6))
   call listkk(dksqmax, cryst%gmet, indkk, self%lgk%ibz, self%bz, self%nq_k, self%nbz, cryst%nsym,&
      sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, use_symrec=.False.)

   if (dksqmax > tol12) then
     write(msg, '(a,es16.6)' ) &
      "At least one of the points in IBZ(k) + q could not be generated from a symmetrical one. dksqmax: ",dksqmax
     MSG_ERROR(msg)
   end if

   !check if same results as listkk
   do ii=1,self%nbz
      ABI_CHECK(indkk(ii,1)==bz2lgkibz(ii),'Unmatching indexes')
   enddo
   ABI_FREE(indkk)
#endif

 ! Build tetrahedron object using IBZ(k) as the effective IBZ
 ! This means that input data for tetra routines must be provided in lgk%kibz_q
 call destroy_tetra(self%tetra_k)
 call init_tetra(bz2lgkibz, cryst%gprimd, self%klatt, self%bz, self%nbz, self%tetra_k, ierr, errorstring)
 if (ierr /= 0) then
   MSG_ERROR(errorstring)
 end if
 ABI_FREE(bz2lgkibz)

end subroutine ephwg_double_grid_setup_kpoint
!!***


!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_get_deltas
!! NAME
!! ephwg_get_deltas
!!
!! FUNCTION
!! Compute weights for $ \delta(\omega - \ee_{k+q, b} \pm \omega_{q\nu} $
!! for a given (band, spin) and phonon mode nu.
!!
!! INPUTS
!! band=band index (global index i.e. unshifted)
!! spin=Spin index
!! nu=Phonon branch.
!! nene=number of energies for DOS
!! eminmax=min and  energy in delta (linear mesh)
!! bcorr=1 to include Blochl correction else 0.
!! comm=MPI communicator
!! [broad]=Gaussian broadening
!!
!! OUTPUT
!!  deltaw_pm(nene, nq_k, 2)  (plus, minus) including the weights for BZ integration.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_get_deltas(self, band, spin, nu, nene, eminmax, bcorr, deltaw_pm, comm, &
                            broad)  ! optional

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band, spin, nu, nene, bcorr, comm
 type(ephwg_t),intent(in) :: self
 real(dp),optional,intent(in) :: broad
!arrays
! This arrays have the same order as the little group used in sigmaph.
 real(dp),intent(in) :: eminmax(2)
 real(dp),intent(out) :: deltaw_pm(nene, self%nq_k, 2)

!Local variables-------------------------------
!scalars
 integer :: iq,iq_ibz,ikpq_ibz,ib,ie
 real(dp),parameter :: max_occ1 = one
 real(dp) :: omega_step
!arrays
 real(dp)  :: thetaw(nene, self%nq_k), wme0(nene), pme_k(self%nq_k, 2)

!----------------------------------------------------------------------

 ib = band - self%bstart + 1

 ! fill array for e_{k+q, b} +- w_{q,nu)
 do iq=1,self%nq_k
   iq_ibz = self%lgk2ibz(iq)   ! IBZ_k --> IBZ
   ikpq_ibz = self%kq2ibz(iq)  ! k + q --> IBZ
   pme_k(iq, 1) = self%eigkbs_ibz(ikpq_ibz, ib, spin) - self%phfrq_ibz(iq_ibz, nu)
   pme_k(iq, 2) = self%eigkbs_ibz(ikpq_ibz, ib, spin) + self%phfrq_ibz(iq_ibz, nu)
 end do

 if (present(broad)) then
   omega_step = (eminmax(2) - eminmax(1)) / (nene - 1)
   ! Use thetaw as workspace array
   thetaw(:, 1) = arth(eminmax(1), omega_step, nene)
   do iq=1,self%nq_k
     do ie=1,2
       wme0 = thetaw(:, 1) - pme_k(iq, ie)
       deltaw_pm(:, iq, ie) = dirac_delta(wme0, broad)
     end do
   end do

   ! Multiply by weights
   do ie=1,nene
     deltaw_pm(ie, :, 1) = deltaw_pm(ie, :, 1) * self%lgk%weights
     deltaw_pm(ie, :, 2) = deltaw_pm(ie, :, 2) * self%lgk%weights
   end do

 else
   ! TODO Add routine to compute only delta
   call tetra_blochl_weights(self%tetra_k, pme_k(:,1), eminmax(1), eminmax(2), max_occ1, nene, self%nq_k, &
     bcorr, thetaw, deltaw_pm(:,:,1), comm)
   call tetra_blochl_weights(self%tetra_k, pme_k(:,2), eminmax(1), eminmax(2), max_occ1, nene, self%nq_k, &
     bcorr, thetaw, deltaw_pm(:,:,2), comm)
 end if

end subroutine ephwg_get_deltas
!!***

!!****f* m_ephwg/ephwg_get_dweights
!! NAME
!! ephwg_get_dweights
!!
!! FUNCTION
!! Compute weights for $ \delta(ww - \ee_{k+q, b} \pm \omega_{q\nu} $
!! for a given (qpoint, band, spin) and all phonon modes.
!!
!! INPUTS
!! band=band index (global index i.e. unshifted)
!! spin=Spin index
!! nw=number of energies for DOS
!! bcorr=1 to include Blochl correction else 0.
!! comm=MPI communicator
!! [use_bzsum]= By default the weights are multiplied by the Nstar(q) / Nq where
!!   Nstar(q) is the number of points in the star of the q-point (using the symmetries of the little group of k)
!!   If use_bzsum is set to True, the Nstar(q) coefficient is removed so that the caller can
!!   integrate over the BZ without using symmetries.
!!
!! OUTPUT
!!  deltaw_pm(nw, nq_k, 2)  (plus, minus) including the weights for BZ integration.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_get_dweights(self, qlklist, nqlk, nw, wvals, spin, bcorr, nbsum, deltaw_pm, comm)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: qlklist(nqlk), nqlk, spin, nw, bcorr, nbsum, comm
 type(ephwg_t),intent(in) :: self
!arrays
 real(dp),intent(in) :: wvals(nw)
 real(dp),intent(inout) :: deltaw_pm(2, nw, self%natom3, nbsum, nqlk)

!Local variables-------------------------------
!scalars
 integer,parameter :: nene3=3
 real(dp),parameter :: max_occ1 = one
 integer :: iq, iq_ibz, ikpq_ibz, ib, ib_k, ie, iw, nu, ii, jj, kk
 integer :: ntetra, itetra, iqlk, iqlk_microzone
!arrays
 integer :: ind_ibz(4), counter, mappingsize
 integer,allocatable :: tetra_mask(:), mapping(:)
 real(dp) :: tmp_deltaw_pm(nw,2,4)
 real(dp) :: pme_k(4, 2), weights(nw, 2)
 real(dp) :: enemin, enemax
 real(dp) :: theta_tmp(nene3,4), delta_tmp(nene3,4), eigen_1tetra(4)

!----------------------------------------------------------------------

 !
 ! This routine still has some bug
 ! The results are not consistent with the calculation where all the weights are precomputed.
 ! 

 ntetra = self%tetra_k%ntetra
 deltaw_pm = 0
 mappingsize = maxval(qlklist)
 ABI_CALLOC(tetra_mask,(ntetra))
 ABI_CALLOC(mapping,(mappingsize))
 !first loop to identify the tetrahedron contribution to the qpoints in the microzone
 do jj=1,nqlk
   !get index of this q in the little group of k
   iqlk = qlklist(jj)
   !map iqlk to microzone
   mapping(iqlk) = jj

   ! Get all the tetrahedron corresponding to this q
   do ii=1,self%tetra_k%ibz_tetra_count(iqlk)
     itetra = self%tetra_k%ibz_tetra_mapping(iqlk,ii)
     !this tetrahedra is contributing to the microzone so we will calculate it
     tetra_mask(itetra) = 1
   end do
 end do

 !second loop to accumulate the contributions of the different tetrahedra
 do itetra=1,ntetra
   !this tetrahedra does not contribute to the microzone
   if (tetra_mask(itetra) == 0) cycle

   ! Here we need the original ordering to reference the correct irred kpoints
   ! map from tetra to ibz
   ind_ibz(:) = self%tetra_k%tetra_full(:,1,itetra)

   !loop over nbsum bands
   ! TODO this has to be consistent with nbsum in sigma
   do ib=self%bstart,self%nbcount
     !loop over phonon freqs
     do nu = 1, self%natom3
       ! Fill array for e_{k+q, b} +- w_{q,nu)
       do ii=1,4
         iq = ind_ibz(ii) 
         iq_ibz = self%lgk2ibz(iq)   ! IBZ_k --> IBZ
         ikpq_ibz = self%kq2ibz(iq)  ! k + q --> IBZ
         pme_k(ii, 1) = self%eigkbs_ibz(ikpq_ibz, ib, spin) - self%phfrq_ibz(iq_ibz, nu)
         pme_k(ii, 2) = self%eigkbs_ibz(ikpq_ibz, ib, spin) + self%phfrq_ibz(iq_ibz, nu)
       end do
       
       do ii = 1,2
         !calculate weights of one tetrahedron
         call tetra_get_onetetra_wvals(self%tetra_k, itetra, pme_k(:,ii), bcorr, &
                                       nw, wvals, tmp_deltaw_pm)
         !this tetrahedron give 4 contributions we map these weights to the array
         do jj = 1,4
           !get index of this q
           iqlk = ind_ibz(jj)
           if (iqlk > mappingsize) cycle !this q-point is not in the microzone
           kk = mapping(iqlk)
           if (kk == 0) cycle !this q-point is not in the microzone
           deltaw_pm(ii,:,nu,ib,kk) = deltaw_pm(ii,:,nu,ib,kk) + tmp_deltaw_pm(:,1,jj) / ( self%lgk%weights(iqlk) )
         end do
       end do !pm
     end do !nu
   end do !ib
 end do !itetra

 ABI_FREE(tetra_mask)
 ABI_FREE(mapping)

end subroutine ephwg_get_dweights
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_zinv_weights
!! NAME
!! ephwg_zinv_weights
!!
!! FUNCTION
!! Compute weights for a given (kpoint, qpoint, spin) for all phonon modes.
!!
!! INPUTS
!! qpt(3)
!! band=band index (global index i.e. unshifted)
!! spin=Spin index
!! nu=Phonon branch.
!! nbsigma
!! zvals
!! comm=MPI communicator
!! [use_bzsum]= By default the weights are multiplied by the Nstar(q) / Nq where
!!   Nstar(q) is the number of points in the star of the q-point (using the symmetries of the little group of k)
!!   If use_bzsum is set to True, the Nstar(q) coefficient is removed so that the caller can
!!   integrate over the BZ without using symmetries.
!!
!! OUTPUT
!!  cweights(nz, 2, nbsigma, %natom3)  (plus, minus)
!!  include weights for BZ integration.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_zinv_weights(self, iqlk, nz, nbsigma, zvals, band, spin, cweights, comm, use_bzsum)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqlk, band, spin, nz, nbsigma,comm
 type(ephwg_t),intent(in) :: self
 logical, optional, intent(in) :: use_bzsum
!arraye
 complex(dpc),intent(in) :: zvals(nz, nbsigma)
 complex(dpc),intent(out) :: cweights(nz, 2, nbsigma, self%natom3)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iq_ibz,ikpq_ibz,ib,ii,jj,iz,itetra,nu,iq,nprocs, my_rank, ierr
 real(dp) :: volconst_mult
 logical :: use_bzsum_
!arrays
 real(dp) :: ework(4, 2)
 real(dp),allocatable :: pme_k(:,:,:)
 integer :: ind_ibz(4)
 !complex(dpc) :: SIM0, SIM0I
 complex(dpc) :: VERM(4), VERL(4), VERLI(4),  cint(nz,4)
!----------------------------------------------------------------------

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 use_bzsum_ = .False.; if (present(use_bzsum)) use_bzsum_ = use_bzsum

 ! Allocate array for e_{k+q, b} +- w_{q,nu)
 ABI_MALLOC(pme_k, (self%nq_k, 2, self%natom3))

 ib = band - self%bstart + 1
 do nu=1,self%natom3
   do iq=1,self%nq_k
     iq_ibz = self%lgk2ibz(iq)   ! IBZ_k --> IBZ
     ikpq_ibz = self%kq2ibz(iq)  ! k + q --> IBZ
     pme_k(iq, 1, nu) = self%eigkbs_ibz(ikpq_ibz, ib, spin) - self%phfrq_ibz(iq_ibz, nu)
     pme_k(iq, 2, nu) = self%eigkbs_ibz(ikpq_ibz, ib, spin) + self%phfrq_ibz(iq_ibz, nu)
   end do
 end do

 cweights = zero
 do itetra = 1, self%tetra_k%ntetra
   ! Here we need the original ordering to reference the correct irred kpoints
   ! Cycle if this tetra does not contribute to this k-point.
   ! See also tetra_get_onewk
   if (all(self%tetra_k%tetra_full(:, 1, itetra) /= iqlk)) cycle
   if (mod(itetra, nprocs) /= my_rank) cycle ! MPI parallelism
   ind_ibz = self%tetra_k%tetra_full(:, 1, itetra)
   !volconst_mult = self%tetra_k%vv / 4.d0 !* dble(self%tetra_k%tetra_mult(itetra))
   volconst_mult = self%tetra_k%vv * dble(self%tetra_k%tetra_mult(itetra))

   do nu=1,self%natom3
     ework = pme_k(ind_ibz(:), :, nu)
     do ib=1,nbsigma
       do ii=1,2
         ! Compute weights for nz points.
         do iz=1,nz
           verm = zvals(iz, ib) - ework(:, ii)
           !verm = zvals(iz, ib) - pme_k(ind_ibz(:), ii, nu)
           !call SIM0ONEI(SIM0, SIM0I, VERM)
           !cint(iz,:) = SIM0I / 4.0_dp * volconst_mult
           call SIM0TWOI(VERL, VERLI, VERM)
           cint(iz,:) = verl(:) * volconst_mult
         end do

         ! TODO
         ! Accumulate contributions to ik_ibz (there might be multiple vertexes that map onto ik_ibz)
         do jj=1,4
           if (ind_ibz(jj) == iqlk) then
             !weights(:,1) = weights(:,1) + dtweightde_tmp(:,jj)
             cweights(:, ii, ib, nu) = cweights(:, ii, ib, nu) + cint(:, jj)
           end if
         end do
       end do  ! +-
     end do ! ib
   end do  ! nu

 end do ! itetra

 ! Rescale weights so that the caller can sum over the full BZ.
 if (use_bzsum_) cweights = cweights / ( self%lgk%weights(iqlk) * self%nbz )

 call xmpi_sum(cweights, comm, ierr)

 ! Compare results with naive one-point integration.
 if (.False. .and. my_rank == master) then
   volconst_mult = one
   volconst_mult = self%lgk%weights(iqlk)
   do ib=1,nbsigma
     do nu=1,self%natom3
       write(std_out,*)"# naive vs tetra integration for band, nu", ib - 1 + self%bstart, nu
       do iz=1,nz
         write(std_out,"(5es16.8)") &
           dble(zvals(iz, ib)), one / (zvals(iz, ib) - pme_k(iqlk, 1, nu)) * volconst_mult, cweights(iz, 1, ib, nu)
         write(std_out,"(5es16.8)") &
           dble(zvals(iz, ib)), one / (zvals(iz, ib) - pme_k(iqlk, 2, nu)) * volconst_mult, cweights(iz, 2, ib, nu)
       end do
     end do
   end do
 end if

 ABI_FREE(pme_k)

end subroutine ephwg_zinv_weights
!!***

!!****f* m_ephwg/ephwg_free
!! NAME
!! ephwg_free
!!
!! FUNCTION
!!  Deallocate memory
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

subroutine ephwg_free(self)

 implicit none

!Arguments ------------------------------------
 type(ephwg_t),intent(inout) :: self

!----------------------------------------------------------------------

 ! integer
 ABI_SFREE(self%kq2ibz)

 ! Real
 ABI_SFREE(self%ibz)
 ABI_SFREE(self%bz)
 ABI_SFREE(self%lgk2ibz)
 ABI_SFREE(self%phfrq_ibz)
 ABI_SFREE(self%eigkbs_ibz)

 ! types
 call destroy_tetra(self%tetra_k)
 call lgroup_free(self%lgk)

 ! nullify pointers
 self%cryst => null()

end subroutine ephwg_free
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_test
!! NAME
!! ephwg_test
!!
!! FUNCTION
!!
!! INPUTS
!! dtset<dataset_type>=All input variables for this dataset.
!! ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

!! subroutine ephwg_test(dtset, cryst, ebands, ifc, prefix, comm)
!!
!!  implicit none
!!
!! !Arguments ------------------------------------
!! !scalarster
!!  integer,intent(in) :: comm
!!  type(dataset_type),intent(in) :: dtset
!!  type(crystal_t),intent(in) :: cryst
!!  type(ebands_t),intent(in) :: ebands
!!  type(ifc_type),intent(in) :: ifc
!!  character(len=*),intent(in) :: prefix
!!
!! !Local variables-------------------------------
!! !scalars
!!  integer,parameter :: bcorr0 = 0, intp_nshiftk = 1, master = 0
!!  integer,parameter :: nz = 2, nbsigma= 1
!!  integer :: band, spin, nu, bstart, nbcount, ik_ibz, ii, it, nene, nprocs, my_rank, cnt, ierr, ntemp
!!  integer :: iq, iq_ibz, ikpq_ibz
!! #ifdef HAVE_NETCDF
!!  integer :: ncid, ncerr
!! #endif
!!  real(dp) :: omega_step
!!  real(dp) :: cpu, wall, gflops
!!  character(len=500) :: msg
!!  type(ephwg_t) :: ephwg
!!  type(ebands_t) :: eb_dense
!!  type(edos_t) :: edos
!! !arrays
!!  integer :: intp_kptrlatt(3, 3), band_block(2)
!!  real(dp):: intp_shiftk(3, intp_nshiftk), eminmax(2)
!!  real(dp),allocatable :: deltaw_pm(:,:,:), nq_list(:),fkq_list(:)
!!  real(dp),allocatable :: sdelta(:,:,:,:,:), omega_mesh(:), ktmesh(:)
!!  real(dp),allocatable :: ekq_list(:), phq_list(:)
!!
!!  !complex(dpc) :: zvals(nz, nbsigma)
!!  !complex(dpc) :: cweights(nz, 2, cryst%natom * 3, nbsigma)
!!
!! !----------------------------------------------------------------------
!!
!!  ! Consistency check
!!  ierr = 0
!!  if (nint(dtset%einterp(1)) == 0) then
!!    MSG_ERROR_NOSTOP("einterp must be specified in the input file", ierr)
!!  end if
!!  if (ierr /= 0) then
!!    MSG_ERROR("Fatal error. See above warnings")
!!  end if
!!
!!  nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
!!
!!  ! Interpolate band energies with star-functions
!!  ! FIXME: One should include all degenerate states from the energy window
!!  ! but it seems there's a problem if I don't include all states.
!!  band_block = [1, 7]
!!  !intp_kptrlatt = reshape([8, 0, 0, 0, 8, 0, 0, 0, 8], [3, 3])
!!  intp_kptrlatt = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3]) * dtset%useria
!!  !ABI_CHECK(all(dtset%bs_interp_kmult /= 0), "bs_interp_kmult == 0!")
!!  !do ii = 1, 3
!!  !  intp_kptrlatt(:, ii) = dtset%bs_interp_kmult(ii) * ebands%kptrlatt(:, ii)
!!  !end do
!!  !intp_nshiftk =
!!  intp_shiftk = zero
!!
!!  call ebands_interpolate_kpath(ebands, dtset, cryst, band_block, prefix, comm)
!!
!!  eb_dense = ebands_interp_kmesh(ebands, cryst, dtset%einterp, intp_kptrlatt, intp_nshiftk, intp_shiftk, band_block, comm)
!!  !call ebands_set_scheme(eb_dense, ebands%occopt, ebands%tsmear, dtset%spinmagntarget, dtset%prtvol)
!!
!!  ! Initalize object for eph weights.
!!  bstart = band_block(1); nbcount = band_block(2) - band_block(1) + 1
!!  ephwg = ephwg_new(cryst, ifc, bstart, nbcount, eb_dense%kptopt, intp_kptrlatt, intp_nshiftk, intp_shiftk, &
!!    eb_dense%nkpt, eb_dense%kptns, eb_dense%nsppol, eb_dense%eig, comm)
!!
!!  ! Frequency mesh for delta functions.
!!  eminmax(1) = minval(eb_dense%eig(band_block(1):band_block(2), :, :))
!!  eminmax(2) = maxval(eb_dense%eig(band_block(1):band_block(2), :, :))
!!  omega_step = 0.1_dp * eV_Ha
!!  nene = 1 + nint((eminmax(2) - eminmax(1)) / omega_step)
!!  ABI_MALLOC(omega_mesh, (nene))
!!  omega_mesh = arth(eminmax(1), omega_step, nene)
!!
!!  ! Temperature mesh.
!!  ntemp = nint(dtset%tmesh(3))
!!  ABI_CHECK(ntemp > 0, "ntemp <= 0")
!!  ABI_MALLOC(ktmesh, (ntemp))
!!  ktmesh = arth(dtset%tmesh(1), dtset%tmesh(2), ntemp) * kb_HaK
!!
!!  ! Print important parameters.
!!  if (my_rank == master) then
!!    write(std_out,"(2a)")ch10," === Parameters used to compute EPH weights ==="
!!    write(std_out,"(a)")sjoin("Band block for SKW interpolation:", ltoa(band_block))
!!    write(std_out,"(a)")sjoin("Number of points in the IBZ:", itoa(eb_dense%nkpt))
!!    write(std_out,"(a)")sjoin("Dense k-mesh:", ltoa(pack(transpose(intp_kptrlatt), mask=.True.)))
!!    write(std_out,"(a)")sjoin("Shifts:", ltoa(pack(intp_shiftk, mask=.True.)))
!!    write(std_out,"(a)")sjoin("Number of frequencies in delta functions:", itoa(nene), &
!!      "From:", ftoa(omega_mesh(1) * Ha_eV), "to", ftoa(omega_mesh(nene) * Ha_eV), "[eV]")
!!    write(std_out,"(a)")sjoin("Number of temperatures:", itoa(ntemp), &
!!      "From:", ftoa(kTmesh(1) / kb_HaK), "to", ftoa(kTmesh(ntemp) / kb_HaK), "[K]")
!!    write(std_out,"(a)")ch10
!!  end if
!!
!!  ! Compute electron DOS with tetra.
!!  edos = ebands_get_edos(eb_dense, cryst, 2, omega_step, zero, comm)
!!
!!  ABI_CALLOC(sdelta, (nene, 2, ntemp, eb_dense%nkpt, eb_dense%nsppol))
!!  cnt = 0
!!  do ik_ibz=1,eb_dense%nkpt
!!    cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! mpi-parallelism
!!    call cwtime(cpu, wall, gflops, "start")
!!    call ephwg_setup_kpoint(ephwg, eb_dense%kptns(:, ik_ibz), dtset%prtvol, xmpi_comm_self)
!!
!!    ABI_MALLOC(deltaw_pm, (nene, ephwg%nq_k, 2))
!!    ABI_MALLOC(phq_list, (ephwg%nq_k))
!!    ABI_MALLOC(nq_list, (ephwg%nq_k))
!!    ABI_MALLOC(ekq_list, (ephwg%nq_k))
!!    ABI_MALLOC(fkq_list, (ephwg%nq_k))
!!
!!    do nu = 1, ifc%natom * 3
!!      ! Load phonons in IBZ(k)
!!      do iq=1,ephwg%nq_k
!!        iq_ibz = ephwg%lgk2ibz(iq)  ! IBZ_k --> IBZ
!!        phq_list(iq) = ephwg%phfrq_ibz(iq_ibz, nu)
!!      end do
!!
!!      do band = band_block(1), band_block(2)
!!        do spin = 1, eb_dense%nsppol
!!          call ephwg_get_deltas(ephwg, band, spin, nu, nene, eminmax, bcorr0, deltaw_pm, xmpi_comm_self)
!!          !call ephwg_get_deltas(ephwg, band, spin, nu, nene, eminmax, bcorr0, deltaw_pm, xmpi_comm_self, broad=0.3 * eV_Ha)
!!
!!          ! Load e_{k+q} in IBZ(k)
!!          do iq=1,ephwg%nq_k
!!            ikpq_ibz = ephwg%kq2ibz(iq)  ! k + q --> IBZ
!!            ekq_list(iq) = eb_dense%eig(band, ikpq_ibz, spin)
!!          end do
!!
!!          ! TODO: Check conventions.
!!          ! Note fermi level taken from ebands.
!!          do it=1,ntemp
!!            nq_list = occ_be(phq_list, ktmesh(it), zero)
!!            fkq_list = occ_fd(ekq_list, ktmesh(it), ebands%fermie)
!!            do iq=1,ephwg%nq_k
!!              sdelta(:, 1, it, ik_ibz, spin) = sdelta(:, 1, it, ik_ibz, spin) +  &
!!                (nq_list(iq) + fkq_list(iq)) * deltaw_pm(:, iq, 1)
!!              sdelta(:, 2, it, ik_ibz, spin) = sdelta(:, 2, it, ik_ibz, spin) +  &
!!                (nq_list(iq) + one - fkq_list(iq)) * deltaw_pm(:, iq, 2)
!!            end do
!!          end do
!!
!!        end do ! spin
!!      end do ! band
!!    end do ! nu
!!
!!    ABI_FREE(deltaw_pm)
!!    ABI_FREE(phq_list)
!!    ABI_FREE(nq_list)
!!    ABI_FREE(ekq_list)
!!    ABI_FREE(fkq_list)
!!
!!    call cwtime(cpu, wall, gflops, "stop")
!!    write(msg,'(2(a,i0),2(a,f8.2))')"k-point [", ik_ibz, "/", eb_dense%nkpt, "] completed. cpu:",cpu,", wall:",wall
!!    call wrtout(std_out, msg, do_flush=.True.)
!!  end do
!!
!!  ! Collect data, then master writes data.
!!  call xmpi_sum_master(sdelta, master, comm, ierr)
!!  if (my_rank == master) then
!! #ifdef HAVE_NETCDF
!!    NCF_CHECK(nctk_open_create(ncid, strcat(prefix, "_SDELTA.nc"), xmpi_comm_self))
!!    NCF_CHECK(crystal_ncwrite(cryst, ncid))
!!    NCF_CHECK(ebands_ncwrite(eb_dense, ncid))
!!    NCF_CHECK(edos_ncwrite(edos, ncid))
!!    ! Add dimensions for SDELTA arrays.
!!    NCF_CHECK(nctk_def_dims(ncid, [nctkdim_t("nomega", nene), nctkdim_t("ntemp", ntemp)], defmode=.True.))
!!    ncerr = nctk_def_arrays(ncid, [ &
!!      nctkarr_t("omega_mesh", "dp", "nomega"), &
!!      nctkarr_t("kTmesh", "dp", "ntemp"), &
!!      nctkarr_t("sdelta", "dp", "nomega, two, ntemp, number_of_kpoints, number_of_spins") &
!!    ])
!!    NCF_CHECK(ncerr)
!!    ! Write data
!!    NCF_CHECK(nctk_set_datamode(ncid))
!!    NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "omega_mesh"), omega_mesh))
!!    NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kTmesh"), ktmesh))
!!    NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "sdelta"), sdelta))
!!    NCF_CHECK(nf90_close(ncid))
!! #endif
!!  end if
!!
!!  ABI_FREE(omega_mesh)
!!  ABI_FREE(ktmesh)
!!  ABI_FREE(sdelta)
!!
!!  call ebands_free(eb_dense)
!!  call edos_free(edos)
!!  call ephwg_free(ephwg)
!!
!! end subroutine ephwg_test
!! !!***

end module m_ephwg
!!***
