!!****m* ABINIT/m_berry_curvature
!! NAME
!! m_berry_curvature
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2022 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_berry_curvature

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 !use m_krank
 !use m_htetra
 use netcdf
 use m_nctk
 use m_crystal
 use m_dtset
 use m_dtfil

 use m_time,            only : cwtime, cwtime_report, sec2str
 use m_fstrings,        only : strcat, sjoin, ktoa !, tolower, itoa, ftoa, ltoa, strcat
 use defs_datatypes,    only : ebands_t
 !use m_ebands,         only : edos_t, ebands_get_edos
 use m_kpts,            only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map, kpts_sort, kpts_pack_in_stars
 use m_ddb_hdr,         only : ddb_hdr_type
 use m_ddb,             only : ddb_type
 use m_ifc,             only : ifc_type
 use m_gstore,          only : gstore_t, gqk_t

 implicit none

 private

 public :: berry_curvature
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_berry_curvature/berry_curvature
!! NAME
!! berry_curvature
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine berry_curvature(gstore, dtset, dtfil, in_ddb)

!Arguments ------------------------------------
!scalars
 type(gstore_t),target,intent(inout) :: gstore
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(ddb_type),intent(in) :: in_ddb

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: nproc, my_rank, ierr, itemp, ntemp !, ncid
 integer :: my_is, spin, nsppol, natom3, band, ib1, ib2, band1, band2, nb, ebands_timrev
 integer :: ik_ibz, isym_k, trev_k, tsign_k, g0_k(3)
 integer :: ikq_ibz, isym_kq, trev_kq, tsign_kq, g0_kq(3)
 integer :: iq_ibz, isym_q, trev_q, tsign_q, g0_q(3)
 integer :: my_iq, iq_glob, my_ik, ik_glob, my_ip1, my_ip2, ipert1, ipert2, iblok
 !integer :: phmesh_size
 logical :: isirr_k, isirr_q, isirr_kq
 real(dp) :: e_kq_b1, e_k_b2,  f_kq_b1, f_k_b2
 real(dp) :: kt, wmax, cpu, wall, gflops
 !real(dp) :: edos_step, edos_broad !, sigma, ecut, eshift, eig0nk
 !character(len=5000) :: msg
 class(crystal_t),pointer :: cryst
 class(ebands_t),pointer :: ebands
 type(gqk_t),pointer :: gqk
 type(ddb_hdr_type) :: ddb_hdr
 type(ddb_type) :: berry_ddb
 type(ifc_type) :: berry_ifc
!arrays
 integer,allocatable :: my_kqmap(:,:), kmesh_map(:,:)
 real(dp) :: kk_bz(3), kq_bz(3), kk_ibz(3), kq_ibz(3), qq_bz(3), qq_ibz(3), qq(3) !, vk(3)
 !real(dp),allocatable :: ktmesh(:), phmesh(:), eig_k(:,:), eig_kq(:,:), kmesh(:,:), wght_bz(:,:,:)
 complex(dp),allocatable :: gmat(:,:,:)
 real(dp), allocatable :: d2matr(:,:,:,:,:)
 integer, allocatable :: flg(:,:,:,:)

!----------------------------------------------------------------------

 nproc = xmpi_comm_size(gstore%comm); my_rank = xmpi_comm_rank(gstore%comm)

 call wrtout(std_out, " Computing berry curvature", pre_newlines=2)
 call cwtime(cpu, wall, gflops, "start")

 cryst => gstore%cryst; ebands => gstore%ebands
 natom3 = 3 * cryst%natom; nsppol = ebands%nsppol

 ! Consistency check
 ierr = 0
 ABI_CHECK_NOSTOP(gstore%qzone == "ibz", "qzone == 'ibz' is required", ierr)
 ABI_CHECK_NOSTOP(gstore%kzone == "bz", "kzone == 'bz' is required", ierr)
 ABI_CHECK_NOSTOP(gstore%gqk(1)%cplex == 2, "cplex == 2 is required", ierr)
 ABI_CHECK(ierr == 0, "Wrong gstore object for berry_curvature. See messages above")

 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)
 ABI_CALLOC(gmat, (natom3, natom3, gstore%nqibz))

 do spin=1,gstore%nsppol
   my_is = gstore%spin2my_is(spin); if (my_is == 0) cycle
   gqk => gstore%gqk(my_is)
   nb = gqk%nb
   ABI_MALLOC(my_kqmap, (6, gqk%my_nk))

   do my_iq=1,gqk%my_nq
     iq_ibz = gqk%my_q2ibz(1, my_iq); isym_q = gqk%my_q2ibz(2, my_iq)
     trev_q = gqk%my_q2ibz(6, my_iq); g0_q = gqk%my_q2ibz(3:5, my_iq)
     isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
     !isirr_q = (isym_q == 1 .and. trev_q == 0)
     tsign_q = 1; if (trev_q == 1) tsign_q = -1
     qq_ibz = gstore%qibz(:, iq_ibz)
     iq_glob = my_iq + gqk%my_qstart - 1

     ! Find k+q in the IBZ for all my k-points.
     if (kpts_map("symrel", ebands_timrev, cryst, gstore%krank_ibz, gqk%my_nk, gqk%my_kpts, my_kqmap, qpt=qq_ibz) /= 0) then
       ABI_ERROR(sjoin("Cannot map k+q to IBZ with qpt:", ktoa(qq_ibz)))
     end if

     ! Compute gmat: integration over my k-points in the BZ + summation over two band indices.
     do my_ik=1,gqk%my_nk
       ik_glob = my_ik + gqk%my_kstart - 1

       ik_ibz = gqk%my_k2ibz(1, my_ik); isym_k = gqk%my_k2ibz(2, my_ik)
       trev_k = gqk%my_k2ibz(6, my_ik); g0_k = gqk%my_k2ibz(3:5, my_ik)
       isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
       tsign_k = 1; if (trev_k == 1) tsign_k = -1

       ikq_ibz = my_kqmap(1, my_ik); isym_kq = my_kqmap(2, my_ik)
       trev_kq = my_kqmap(6, my_ik); g0_kq = my_kqmap(3:5, my_ik)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
       tsign_kq = 1; if (trev_kq == 1) tsign_kq = -1

       do ib2=1,nb
         band2 = ib2 + gqk%bstart - 1
         e_k_b2 = ebands%eig(band2, ik_ibz, spin)
         f_k_b2 = ebands%occ(band2, ik_ibz, spin)
         !g2 = gaussian(ebands%eig(band2, ik_ibz, spin) - ebands%fermie, sigma)
         do ib1=1,nb
           band1 = ib1 + gqk%bstart - 1
           e_kq_b1 = ebands%eig(band1, ikq_ibz, spin)
           f_kq_b1 = ebands%occ(band1, ikq_ibz, spin)
           !g1 = gaussian(ebands%eig(band1, ikq_ibz, spin) - ebands%fermie, sigma)
           do my_ip1=1,gqk%my_npert
             ipert1 = gqk%my_iperts(my_ip1)
             do my_ip2=1,gqk%my_npert
               ipert2 = gqk%my_iperts(my_ip2)
               ! Need all perts on the same procs as we have to take the outer product (ipert1, ipert2)
               ! (2, my_npert, nb, my_nq, nb, my_nk)
               !gmat(ipert1, ipert2, iq_ibz) = gmat(ipert1, ipert2, iq_bz) + &
               !  gqk%my_g(:,my_ip1,ib1,my_iq,ib2,my_ik) * gqk%my_g(:,my_ip2,ib1,my_iq,ib2,my_ik)

               ! TODO: add finite delta or do Taylor series expansion of numerator + denominator?
             end do
           end do
         end do
       end do
     end do ! my_ik
   end do ! my_iq
   ABI_FREE(my_kqmap)
 end do ! spin

 gmat = j_dpc * gmat / gstore%nkbz
 call xmpi_sum(gmat, gstore%comm, ierr)
 call cwtime_report(" berry_curvature:", cpu, wall, gflops)

 ! Get the header of in_ddb
 !call ddb_hdr_open_read(ddb_hdr, filename, unddb, comm, &
 !                       matom, mtypat, mband, mkpt, msym, dimekb, lmnmax, usepaw, dimonly)
 !call ddb_hdr%free()

 ! Create berry_ddb object to add gmat term so that we can build new ifc including Berry curvature.
 ! Note that ddb%val(2,msize,nblok) are in Cartesian coordinates whereas gmat is in reduced coords.
 call in_ddb%copy(berry_ddb)

 ABI_MALLOC(d2matr, (2,3,in_ddb%mpert,3,in_ddb%mpert))
 ABI_MALLOC(flg, (3,in_ddb%mpert,3,in_ddb%mpert))

 do iblok=1, in_ddb%nblok
   !if (abs(in_ddb%typ(iblok)) == abs(rfmeth)) then
   qq(:) = in_ddb%qpt(:,iblok) / in_ddb%nrm(:,iblok)
   !d2matr = in_ddb%d2matr
   !flg = in_ddb%flg
   !call berry_ddb%set_d2matr(d2matr, flg, iblok)
 end do

 ABI_FREE(d2matr)
 ABI_FREE(flg)

 !if (my_rank == 0) then
 !  call berry_ddb%write(ddb_hdr, filename, fullinit, comm)
 !end if

 !call ifc_init(berry_ifc, cryst, berry_ddb, &
 !  dtset%brav, dtset%asr, dtset%symdynmat, dtset%dipdip, dtset%rfmeth, &
 !  dtset%ddb_ngqpt, ddb_nqshift, ddb_qshifts, dielt, zeff, &
 !  qdrp_cart, nsphere0, dtset%rifcsph, prtsrlr0, dtset%enunit, comm, &
 !  dipquad=dtset%dipquad, quadquad=dtset%quadquad)
 !call berry_ddb%free()
 !call berry_ifc%print(unit=std_out)

 ! Output phonon band structure (requires qpath)
 !if (dtset%prtphbands /= 0) call ifc_mkphbs(ifc, cryst, dtset, dtfil%filnam_ds(4), comm)
 !call berry_ifc%free()
 ABI_FREE(gmat)

end subroutine berry_curvature
!!***

end module m_berry_curvature
!!***
