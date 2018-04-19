!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ephwg
!! NAME
!! m_ephwg
!!
!! FUNCTION
!!  Tools and objects to compute the weights used for the BZ integration of EPH quantities.
!!  More specifically the integration of quantities such as the imaginary part of the self-energy
!!  involving delta functions.
!!  Different approaches are available:
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
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_copy
 use m_tetrahedron
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_crystal
 use m_crystal_io
 use m_ifc
 use m_lgroup
 use m_ebands

 use defs_datatypes,    only : ebands_t
 use defs_abitypes,     only : dataset_type
 use m_time,            only : cwtime, sec2str
 use m_numeric_tools,   only : arth
 use m_special_funcs,   only : dirac_delta
 use m_fstrings,        only : strcat, ltoa, itoa, ftoa, sjoin
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
!!  Stores electron eigevalues and phonon frequencies in the IBZ.
!!  (use same mesh for e and ph)
!!  Provides tools to compute (e_{k+q} - w{q}) in the IBZ(k)
!!  and integrate the delta functions with the tetrahedron method.
!!
!! SOURCE

type, public :: ephwg_t

  integer :: natom3
  ! 3 * natom

  integer :: nsppol
  ! Number of independent spin polarizations.

  integer :: nbcount
  ! Number of bands treated

  integer :: bstart
  ! The fist band starts at bstart. Used to select bands
  ! aroudh the Fermi level.

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
  ! reciprocal of lattice vectors for full kpoint grid.
  ! Used by init_tetra

  integer,allocatable :: lgk2ibz(:)
  ! lgk2ibz(nq_k)
  ! Mapping Litte-group --> initial IBZ

  real(dp),allocatable :: phfrq_ibz(:,:)
  ! (nibz, natom3)
  ! Phonon frequencies in the IBZ

  real(dp),allocatable :: eigkbs_ibz(:, :, :)
  ! (nibz, nbcount, nsppol)
  ! Electron eigenvalues in the IBZ for nbcount states (not necessarly equal to nband)

  type(crystal_t), pointer :: cryst => null()
  ! Pointer to input structure

  type(lgroup_t) :: lgrp_k
  ! Little group of the k-point

  type(t_tetrahedron) :: tetra_k
  ! tetrahedron object used to evaluate delta(w - e_{k+q} - phw_q)

 end type ephwg_t

 public :: ephwg_new             ! Constructor
 public :: ephwg_set_kpoint      ! Prepare tetrahedron method for k-point
 public :: ephwg_delta_weights   ! Compute weights for $\delta(\omega - \ee_{k+q, b} \pm \omega_{q\nu} $
 public :: ephwg_free            ! Free memory

 public :: ephwg_test

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
!!  bstart
!!  nbcount
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ephwg_new'
!End of the abilint section

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

!!****f* m_ephwg/ephwg_set_kpoint
!! NAME
!! ephwg_set_kpoint
!!
!! FUNCTION
!!  Set internal tables and object required to compute the delta functions with a given k-point.
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

subroutine ephwg_set_kpoint(self, kpoint)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ephwg_set_kpoint'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ephwg_t),target,intent(inout) :: self
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
 call lgroup_free(self%lgrp_k)
 self%lgrp_k = lgroup_new(self%cryst, kpoint, self%timrev, self%nbz, self%bz, self%nibz, self%ibz)
 self%nq_k = self%lgrp_k%nibz
 !call lgroup_print(self%lgrp)

 ! Get mapping IBZ_k --> initial IBZ (self%lgrp%ibz --> self%ibz)
 ABI_MALLOC(indkk, (self%nq_k * sppoldbl1, 6))
 call listkk(dksqmax, cryst%gmet, indkk, self%ibz, self%lgrp_k%ibz, self%nibz, self%nq_k, cryst%nsym,&
    sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, use_symrec=.False.)

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
   self%lgrp_k%ibz(:, ii) = self%lgrp_k%ibz(:, ii) + kpoint
 end do
 ABI_MALLOC(indkk, (self%nq_k * sppoldbl1, 6))

 call listkk(dksqmax, cryst%gmet, indkk, self%ibz, self%lgrp_k%ibz, self%nibz, self%nq_k, cryst%nsym,&
    sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, use_symrec=.False.)

 if (dksqmax > tol12) then
   write(msg, '(a,es16.6)' ) &
    "At least one of the points in IBZ(k) + k could not be generated from a symmetrical one. dksqmax: ",dksqmax
   MSG_ERROR(msg)
 end if
 if (allocated(self%kq2ibz)) then
   ABI_FREE(self%kq2ibz)
 end if
 call alloc_copy(indkk(:, 1), self%kq2ibz)
 ABI_FREE(indkk)
 do ii=1,self%nq_k
   self%lgrp_k%ibz(:, ii) = self%lgrp_k%ibz(:, ii) - kpoint
 end do

 ! Get mapping BZ --> IBZ_k (self%bz --> self%lgrp%ibz) required for tetrahedron method
 ABI_MALLOC(indkk, (self%nbz * sppoldbl1, 6))
 call listkk(dksqmax, cryst%gmet, indkk, self%lgrp_k%ibz, self%bz, self%nq_k, self%nbz, cryst%nsym,&
    sppoldbl1, cryst%symafm, cryst%symrel, self%timrev, use_symrec=.False.)

 ! Build tetrahedron object using IBZ(k) as the effective IBZ
 ! This means that input data for tetra routines must be provided in lgrp_k%kibz_q
 call destroy_tetra(self%tetra_k)
 call init_tetra(indkk(:, 1), cryst%gprimd, self%klatt, self%bz, self%nbz, self%tetra_k, ierr, errorstring)
 if (ierr /= 0) then
   MSG_ERROR(errorstring)
 end if
 ABI_FREE(indkk)

end subroutine ephwg_set_kpoint
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_delta_weights
!! NAME
!! ephwg_delta_weights
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
!!  wdt_plus(nene, nq_k, 2, 2)  (plus, minus)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_delta_weights(self, band, spin, nu, nene, eminmax, bcorr, wdt_plus, comm,
                               broad)  ! optional


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ephwg_delta_weights'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band, spin, nu, nene, bcorr, comm
 type(ephwg_t),intent(in) :: self
 real(dp),optional,intent(in) :: broad
!arrays
! This arrays have the same order as the little group used in sigmaph.
 real(dp),intent(in) :: eminmax(2)
 real(dp),intent(out) :: wdt_plus(nene, self%nq_k, 2, 2)

!Local variables-------------------------------
!scalars
 integer :: iq_k,iq_ibz,ikpq_ibz,ib,ie
 real(dp),parameter :: max_occ1 = one
 real(dp) :: omega_step
!arrays
 real(dp)  :: dtweightde_t(nene, self%nq_k),tweight_t(nene, self%nq_k)
 real(dp),allocatable :: pme_k(:,:)

!----------------------------------------------------------------------

 ! Allocate array for e_{k+q, b} +- w_{q,nu)
 ABI_MALLOC(pme_k, (self%nq_k, 2))
 ib = band - self%bstart + 1

 do iq_k=1,self%nq_k
   iq_ibz = self%lgk2ibz(iq_k)   ! IBZ_k --> IBZ
   ikpq_ibz = self%kq2ibz(iq_k)  ! k + q --> IBZ
   pme_k(iq_k, 1) = self%eigkbs_ibz(ikpq_ibz, ib, spin) - self%phfrq_ibz(iq_ibz, nu)
   pme_k(iq_k, 2) = self%eigkbs_ibz(ikpq_ibz, ib, spin) + self%phfrq_ibz(iq_ibz, nu)
 end do

 if present(broad) then
   omega_step = (einmax(2) - eminmax(1)) / (nene - 1)
   ! workspace array
   !eminmax(1) + (ie - 1) * omega_step
   !dtweightde_t(nene, 1)
   !do iq_k=1,self%nq_k
   !wme0 = edos%mesh - pme_k(iq_k, 1)
   !wdt_plus(:, iq_k, 1, 1) = dirac_delta(wme0, broad)
   !end do
 else

   call tetra_blochl_weights(self%tetra_k, pme_k(:,1), eminmax(1), eminmax(2), max_occ1, nene, self%nq_k, &
     bcorr, wdt_plus(:,:,1,2), wdt_plus(:,:,1,1), comm)

   call tetra_blochl_weights(self%tetra_k, pme_k(:,2), eminmax(1), eminmax(2), max_occ1, nene, self%nq_k, &
     bcorr, wdt_plus(:,:,2,2), wdt_plus(:,:,2,1), comm)
 end if

 ! Multiply by weights
 !do ib=1,2
 !  do ie=1,nene
 !    wdt_plus(ie, :, :, ib) = wdt_plus(ie, :, :, ib) * self%lgrp_k%weights
 !  end do
 !end do

 ABI_FREE(pme_k)
 !call tetra_get_onewk(self%tetra_k, ik_ibz, bcorr, nene, nkibz, eig_ibz, eminmax(1), eminmax(2), max_occ1, weights)

end subroutine ephwg_delta_weights
!!***

!----------------------------------------------------------------------

!!****f* m_ephwg/ephwg_gdz_weights
!! NAME
!! ephwg_gdz_weights
!!
!! FUNCTION
!! Compute weights
!! for a given (kpoint, iq_k, spin) for all phonon modes.
!!
!! INPUTS
!! band=band index (global index i.e. unshifted)
!! spin=Spin index
!! nu=Phonon branch.
!! comm=MPI communicator
!!
!! OUTPUT
!!  wdt_plus(nene, nq_k, 2, 2)  (plus, minus)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephwg_gdz_weights(self, iqlk, band, spin, nz, nbsigma, zvals, cweights, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ephwg_gdz_weights'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqlk, band, spin, nz, nbsigma,comm
 type(ephwg_t),intent(in) :: self
!arrays
 complex(dpc),intent(in) :: zvals(nz, nbsigma)
 complex(dpc),intent(out) :: cweights(nz, 2, self%natom3, nbsigma)

!Local variables-------------------------------
!scalars
 integer :: iq_ibz,ikpq_ibz,ib,ii,jj,iz,itetra,nu,iq_k,nprocs, my_rank
 real(dp) :: volconst_mult
!arrays
 real(dp),allocatable :: pme_k(:,:,:)
 integer :: ind_ibz(4)
 !complex(dpc) :: SIM0, SIM0I
 complex(dpc) :: VERM(4), VERL(4), VERLI(4),  cint(nz,4)
!----------------------------------------------------------------------

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Allocate array for e_{k+q, b} +- w_{q,nu)
 ABI_MALLOC(pme_k, (self%nq_k, 2, self%natom3))
 write(std_out,*)"iqlk:", iqlk

 ib = band - self%bstart + 1
 do nu=1,self%natom3
   do iq_k=1,self%nq_k
     iq_ibz = self%lgk2ibz(iq_k)   ! IBZ_k --> IBZ
     ikpq_ibz = self%kq2ibz(iq_k)  ! k + q --> IBZ
     pme_k(iq_k, 1, nu) = self%eigkbs_ibz(ikpq_ibz, ib, spin) - self%phfrq_ibz(iq_ibz, nu)
     pme_k(iq_k, 2, nu) = self%eigkbs_ibz(ikpq_ibz, ib, spin) + self%phfrq_ibz(iq_ibz, nu)
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
     do ib=1,nbsigma
       do ii=1,2
         ! Compute weights for nz points.
         do iz=1,nz
           verm = pme_k(ind_ibz(:), ii, nu) + zvals(iz, ib)
           !call SIM0ONEI(SIM0, SIM0I, VERM)
           !cint(iz,:) = SIM0I / 4.0_dp * volconst_mult
           call SIM0TWOI(VERL, VERLI, VERM)
           cint(iz,:) = verl(:) * volconst_mult
         end do

         ! TODO
         ! Accumulate contributions to ik_ibz
         ! (there might be multiple vertexes that map onto ik_ibz)
         do jj=1,4
           if (ind_ibz(jj) == iqlk) then
             !weights(:,1) = weights(:,1) + dtweightde_tmp(:,jj)
             cweights(:, ii, nu, ib) = cweights(:, ii, nu, ib) + cint(:, jj)
           end if
         end do
       end do  ! +-
     end do ! ib
   end do  ! nu

 end do ! itetra

 call xmpi_sum(cweights, comm, ierr)

 ! Compare results with naive integration.
 if (my_rank == master) then
   volconst_mult = self%lgrp_k%weights(iqlk)
   do ib=1,nbsigma
     do nu=1,self%natom3
       write(std_out,*)"# naive vs tetra integration for band, nu", band, nu
       do iz=1,nz
         write(std_out,*) zvals(iz) one / (zvals(iz, ib) + pme_k(iqlk, 1, nu)) * volconst_mult, cweights(iz, 1, nu, ib)
         write(std_out,*) zvals(iz, one / (zvals(iz, ib) + pme_k(iqlk, 2, nu)) * volconst_mult, cweights(iz, 2, nu, ib)
       end do
     end do
   end do
 end if

 ABI_FREE(pme_k)

end subroutine ephwg_gdz_weights
!!***

!!****f* m_ephwg/ephwg_free
!! NAME
!! ephwg_free
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

subroutine ephwg_free(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ephwg_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(ephwg_t),intent(inout) :: self

!----------------------------------------------------------------------

 ! integer
 if (allocated(self%kq2ibz)) then
   ABI_FREE(self%kq2ibz)
 end if

 ! Real
 if (allocated(self%ibz)) then
   ABI_FREE(self%ibz)
 end if
 if (allocated(self%bz)) then
   ABI_FREE(self%bz)
 end if
 if (allocated(self%lgk2ibz)) then
   ABI_FREE(self%lgk2ibz)
 end if
 if (allocated(self%phfrq_ibz)) then
   ABI_FREE(self%phfrq_ibz)
 end if
 if (allocated(self%eigkbs_ibz)) then
   ABI_FREE(self%eigkbs_ibz)
 end if

 ! types
 call destroy_tetra(self%tetra_k)
 call lgroup_free(self%lgrp_k)

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

subroutine ephwg_test(dtset, cryst, ebands, ifc, prefix, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ephwg_test'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalarster
 integer,intent(in) :: comm
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(ifc_type),intent(in) :: ifc
 character(len=*),intent(in) :: prefix

!Local variables-------------------------------
!scalars
 integer,parameter :: bcorr0 = 0, intp_nshiftk = 1, master = 0
 integer :: band, spin, nu, bstart, nbcount, ik_ibz, ie, it, nene, nprocs, my_rank, cnt, ierr, ntemp
 integer :: iq_k, iq_ibz, ikpq_ibz
#ifdef HAVE_NETCDF
 integer :: ncid, ncerr
#endif
 real(dp) :: omega_step
 real(dp) :: cpu, wall, gflops
 character(len=500) :: msg
 type(ephwg_t) :: ephwg
 type(ebands_t) :: eb_dense
 type(edos_t) :: edos
!arrays
 integer :: intp_kptrlatt(3, 3), band_block(2)
 real(dp):: intp_shiftk(3, intp_nshiftk), eminmax(2)
 real(dp),allocatable :: wdt_plus(:,:,:,:), nq_list(:),fkq_list(:)
 real(dp),allocatable :: sdelta(:,:,:,:,:), omega_mesh(:), ktmesh(:)
 real(dp),allocatable :: ekq_list(:), phq_list(:)
 integer,parameter :: nz = 2, nbsigma= 1
 complex(dpc) :: zvals(nz, nbsigma)
 complex(dpc) :: cweights(nz, 2, cryst%natom * 3, nbsigma)

!----------------------------------------------------------------------

 ! Consistency check
 ierr = 0
 if (dtset%einterp(1) == 0) then
   MSG_ERROR_NOSTOP("einterp must be specified in the input file", ierr)
 end if
 if (ierr /= 0) then
   MSG_ERROR("Fatal error. See above warnings")
 end if

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Interpolate band energies with star-functions
 band_block = [1, 8]
 !intp_kptrlatt = reshape([8, 0, 0, 0, 8, 0, 0, 0, 8], [3, 3])
 intp_kptrlatt = reshape([12, 0, 0, 0, 12, 0, 0, 0, 12], [3, 3])
 !intp_kptrlatt = reshape([16, 0, 0, 0, 16, 0, 0, 0, 16], [3, 3])
 !intp_kptrlatt = reshape([24, 0, 0, 0, 24, 0, 0, 0, 24], [3, 3])
 !intp_kptrlatt = intp_kptrlatt * 2
 intp_shiftk = zero

 call ebands_interpolate_kpath(ebands, dtset, cryst, band_block, prefix, comm)

 eb_dense = ebands_interp_kmesh(ebands, cryst, dtset%einterp, intp_kptrlatt, intp_nshiftk, intp_shiftk, band_block, comm)
 !call ebands_set_scheme(eb_dense, ebands%occopt, ebands%tsmear, dtset%spinmagntarget, dtset%prtvol)

 ! Init object for eph weights.
 bstart = band_block(1); nbcount = band_block(2) - band_block(1) + 1
 ephwg = ephwg_new(cryst, ifc, bstart, nbcount, eb_dense%kptopt, intp_kptrlatt, intp_nshiftk, intp_shiftk, &
   eb_dense%nkpt, eb_dense%kptns, eb_dense%nsppol, eb_dense%eig, comm)

 ! Frequency mesh for delta functions.
 eminmax(1) = minval(eb_dense%eig(band_block(1):band_block(2), :, :))
 eminmax(2) = maxval(eb_dense%eig(band_block(1):band_block(2), :, :))
 omega_step = 0.1_dp * eV_Ha
 nene = 1 + nint((eminmax(2) - eminmax(1)) / omega_step)
 ABI_MALLOC(omega_mesh, (nene))
 do ie = 1, nene
   omega_mesh(ie) = eminmax(1) + (ie - 1) * omega_step
 end do

 ! Temperature mesh.
 ntemp = nint(dtset%tmesh(3))
 ABI_CHECK(ntemp > 0, "ntemp <= 0")
 ABI_MALLOC(ktmesh, (ntemp))
 ktmesh = arth(dtset%tmesh(1), dtset%tmesh(2), ntemp) * kb_HaK

 ! Print important parameters.
 if (my_rank == master) then
   write(std_out,"(a)")sjoin("Band block for SKW interpolation:", ltoa(band_block))
   write(std_out,"(a)")sjoin("Number of points in the IBZ:", itoa(eb_dense%nkpt))
   write(std_out,"(a)")sjoin("Dense k-mesh:", ltoa(pack(transpose(intp_kptrlatt), mask=.True.)))
   write(std_out,"(a)")sjoin("Shifts:", ltoa(pack(intp_shiftk, mask=.True.)))
   write(std_out,"(a)")sjoin("Number of frequencies in delta functions:", itoa(nene), &
     "From:", ftoa(omega_mesh(1) * Ha_eV), "to", ftoa(omega_mesh(nene) * Ha_eV), "[eV]")
   write(std_out,"(a)")sjoin("Number of temperatures:", itoa(ntemp), &
     "From:", ftoa(kTmesh(1) / kb_HaK), "to", ftoa(kTmesh(ntemp) / kb_HaK), "[K]")
 end if

 ! Compute electron DOS with tetra.
 edos = ebands_get_edos(eb_dense, cryst, 2, omega_step, zero, comm)

 ABI_CALLOC(sdelta, (nene, 2, ntemp, eb_dense%nkpt, eb_dense%nsppol))
 cnt = 0
 do ik_ibz=1,eb_dense%nkpt
   cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! mpi-parallelism
   call cwtime(cpu, wall, gflops, "start")
   call ephwg_set_kpoint(ephwg, eb_dense%kptns(:, ik_ibz))

#if 1
  zvals(1, 1) = eb_dense%eig(1, 10, 1) + j_dpc * dtset%zcut !* tol12
  zvals(2, 1) = eb_dense%eig(2, 10, 1) + j_dpc * dtset%zcut !* tol12
  call ephwg_gdz_weights(self=ephwg, iqlk=1, band=5, spin=1, nz=nz, nbsigma=nbsigma, zvals=zvals, cweights=cweights, comm=comm)
  call ephwg_gdz_weights(self=ephwg, iqlk=14, band=5, spin=1, nz=nz, nbsigma=nbsigma, zvals=zvals, cweights=cweights, comm=comm)
  stop "hello"
#endif

   ABI_MALLOC(wdt_plus, (nene, ephwg%nq_k, 2, 2))
   ABI_MALLOC(phq_list, (ephwg%nq_k))
   ABI_MALLOC(nq_list, (ephwg%nq_k))
   ABI_MALLOC(ekq_list, (ephwg%nq_k))
   ABI_MALLOC(fkq_list, (ephwg%nq_k))

   do nu = 1, ifc%natom * 3
     ! Load phonons in IBZ(k)
     do iq_k=1,ephwg%nq_k
       iq_ibz = ephwg%lgk2ibz(iq_k)  ! IBZ_k --> IBZ
       phq_list(iq_k) = ephwg%phfrq_ibz(iq_ibz, nu)
     end do

     do band = band_block(1), band_block(2)
       do spin = 1, eb_dense%nsppol
         call ephwg_delta_weights(ephwg, band, spin, nu, nene, eminmax, bcorr0, wdt_plus, xmpi_comm_self)

         ! Load e_{k+q} in IBZ(k)
         do iq_k=1,ephwg%nq_k
           ikpq_ibz = ephwg%kq2ibz(iq_k)  ! k + q --> IBZ
           ekq_list(iq_k) = eb_dense%eig(band, ikpq_ibz, spin)
         end do

         ! TODO: Check conventions.
         ! Note fermi level taken from ebands.
         do it=1,ntemp
           nq_list = occ_be(phq_list, ktmesh(it), zero)
           fkq_list = occ_fd(ekq_list, ktmesh(it), ebands%fermie)
           do iq_k=1,ephwg%nq_k
             sdelta(:, 1, it, ik_ibz, spin) = sdelta(:, 1, it, ik_ibz, spin) +  &
               (nq_list(iq_k) + fkq_list(iq_k)) * wdt_plus(:, iq_k, 1, 1)
             sdelta(:, 2, it, ik_ibz, spin) = sdelta(:, 2, it, ik_ibz, spin) +  &
               (nq_list(iq_k) + one - fkq_list(iq_k)) * wdt_plus(:, iq_k, 2, 1)
           end do
         end do

       end do ! spin
     end do ! band
   end do ! nu

   ABI_FREE(wdt_plus)
   ABI_FREE(phq_list)
   ABI_FREE(nq_list)
   ABI_FREE(ekq_list)
   ABI_FREE(fkq_list)

   call cwtime(cpu, wall, gflops, "stop")
   write(msg,'(2(a,i0),2(a,f8.2))')"k-point [", ik_ibz, "/", eb_dense%nkpt, "] completed. cpu:",cpu,", wall:",wall
   call wrtout(std_out, msg, do_flush=.True.)
 end do

 ! Collect data, then master writes data.
 call xmpi_sum_master(sdelta, master, comm, ierr)
 if (my_rank == master) then
#ifdef HAVE_NETCDF
   NCF_CHECK(nctk_open_create(ncid, strcat(prefix, "_SDELTA.nc"), xmpi_comm_self))
   NCF_CHECK(crystal_ncwrite(cryst, ncid))
   NCF_CHECK(ebands_ncwrite(eb_dense, ncid))
   NCF_CHECK(edos_ncwrite(edos, ncid))
   ! Add dimensions for SDELTA arrays.
   NCF_CHECK(nctk_def_dims(ncid, [nctkdim_t("nomega", nene), nctkdim_t("ntemp", ntemp)], defmode=.True.))
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("omega_mesh", "dp", "nomega"), &
     nctkarr_t("kTmesh", "dp", "ntemp"), &
     nctkarr_t("sdelta", "dp", "nomega, two, ntemp, number_of_kpoints, number_of_spins") &
   ])
   NCF_CHECK(ncerr)
   ! Write data
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "omega_mesh"), omega_mesh))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kTmesh"), ktmesh))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "sdelta"), sdelta))
   NCF_CHECK(nf90_close(ncid))
#endif
 end if

 ABI_FREE(omega_mesh)
 ABI_FREE(ktmesh)
 ABI_FREE(sdelta)

 call ebands_free(eb_dense)
 call edos_free(edos)
 call ephwg_free(ephwg)

end subroutine ephwg_test
!!***

end module m_ephwg
!!***
