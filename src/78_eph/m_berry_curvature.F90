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
 use netcdf
 use m_nctk
 use m_crystal
 use m_dtset
 use m_dtfil

 use m_time,            only : cwtime, cwtime_report, sec2str
 use m_fstrings,        only : strcat, sjoin, ktoa
 use defs_datatypes,    only : ebands_t
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

subroutine berry_curvature(gstore, dtset, dtfil, in_ddb, in_ddb_hdr, berry_ddb)

!Arguments ------------------------------------
!scalars
 type(gstore_t),target,intent(inout) :: gstore
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(ddb_type),intent(in) :: in_ddb
 type(ddb_hdr_type),intent(inout) :: in_ddb_hdr
 type(ddb_type),intent(out) :: berry_ddb

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: nproc, my_rank, ierr, itemp, ntemp
 integer :: my_is, spin, nsppol, natom3, band, ib1, ib2, band1, band2, nb, ebands_timrev
 integer :: ik_ibz, isym_k, trev_k, tsign_k, g0_k(3)
 integer :: ikq_ibz, isym_kq, trev_kq, tsign_kq, g0_kq(3)
 integer :: iq_ibz, isym_q, trev_q, tsign_q, g0_q(3)
 integer :: my_iq, iq_glob, my_ik, ik_glob, my_ip1, my_ip2, ipc1, ipc2, ipert1, ipert2, iblok, idir1, idir2
 logical :: isirr_k, isirr_q, isirr_kq
 real(dp) :: e_kq_b1, e_k_b2,  f_kq_b1, f_k_b2, dene, inv_dene2, spin_occ
 real(dp) :: kt, wmax, cpu, wall, gflops
 !character(len=5000) :: msg
 character(len=fnlen) :: ddb_filepath
 class(crystal_t),pointer :: cryst
 class(ebands_t),pointer :: ebands
 type(gqk_t),pointer :: gqk
!arrays
 integer :: units(2)
 integer,allocatable :: my_kqmap(:,:) !, kmesh_map(:,:), flg(:,:,:,:)
 real(dp) :: kk_bz(3), kq_bz(3), kk_ibz(3), kq_ibz(3), qq_bz(3), qq_ibz(3), ddb_qq(3) !, vk(3)
 !real(dp),allocatable :: ktmesh(:), phmesh(:), eig_k(:,:), eig_kq(:,:), kmesh(:,:), wght_bz(:,:,:)
 real(dp), allocatable :: d2matr(:,:,:,:,:)
 complex(dp),allocatable :: gmat(:,:,:)

!----------------------------------------------------------------------

 nproc = xmpi_comm_size(gstore%comm); my_rank = xmpi_comm_rank(gstore%comm)
 units = [std_out, ab_out]

 call wrtout(std_out, " Computing berry curvature", pre_newlines=2)
 call cwtime(cpu, wall, gflops, "start")

 cryst => gstore%cryst; ebands => gstore%ebands
 natom3 = 3 * cryst%natom; nsppol = ebands%nsppol
 spin_occ = one; if (nsppol == 1 .and. dtset%nspinor == 1) spin_occ = two
 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)

 ! Consistency check
 ierr = 0
 ABI_CHECK_NOSTOP(gstore%qzone == "ibz", "qzone == 'ibz' is required", ierr)
 ABI_CHECK_NOSTOP(gstore%kzone == "bz", "kzone == 'bz' is required", ierr)
 ABI_CHECK_NOSTOP(gstore%gqk(1)%cplex == 2, "cplex == 2 is required", ierr)
 ABI_CHECK_NOSTOP(gstore%gqk(1)%pert_comm%nproc == 1, "berry_curvature is not compatible with pert_parallelism", ierr)
 ABI_CHECK(ierr == 0, "Wrong gstore object for berry_curvature. See messages above")

 ABI_CALLOC(gmat, (natom3, natom3, gstore%nqibz))

 ! Loop over collinear spins (if any)
 do spin=1,gstore%nsppol
   my_is = gstore%spin2my_is(spin); if (my_is == 0) cycle
   gqk => gstore%gqk(my_is)
   nb = gqk%nb
   ABI_MALLOC(my_kqmap, (6, gqk%my_nk))

   ! For each q-point in the IBZ treated by me.
   do my_iq=1,gqk%my_nq
     iq_ibz = gqk%my_q2ibz(1, my_iq); isym_q = gqk%my_q2ibz(2, my_iq)
     trev_q = gqk%my_q2ibz(6, my_iq); g0_q = gqk%my_q2ibz(3:5, my_iq)
     isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
     tsign_q = 1; if (trev_q == 1) tsign_q = -1
     qq_ibz = gstore%qibz(:, iq_ibz)
     iq_glob = my_iq + gqk%my_qstart - 1

     ! Find k+q in the IBZ for all my k-points.
     if (kpts_map("symrel", ebands_timrev, cryst, gstore%krank_ibz, gqk%my_nk, gqk%my_kpts, my_kqmap, qpt=qq_ibz) /= 0) then
       ABI_ERROR(sjoin("Cannot map k+q to IBZ with qpt:", ktoa(qq_ibz)))
     end if

     ! Integration over my k-points in the BZ
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

       ! Summation over the two band indices. occupancies f are rescaled in [0, 1] when nsppol 2.
       do ib2=1,nb
         band2 = ib2 + gqk%bstart - 1
         e_k_b2 = ebands%eig(band2, ik_ibz, spin)
         f_k_b2 = ebands%occ(band2, ik_ibz, spin) / spin_occ
         !g2 = gaussian(ebands%eig(band2, ik_ibz, spin) - ebands%fermie, sigma)
         do ib1=1,nb
           band1 = ib1 + gqk%bstart - 1
           e_kq_b1 = ebands%eig(band1, ikq_ibz, spin)
           f_kq_b1 = ebands%occ(band1, ikq_ibz, spin) / spin_occ
           !g1 = gaussian(ebands%eig(band1, ikq_ibz, spin) - ebands%fermie, sigma)
           dene = e_kq_b1 - e_k_b2
           ! TODO: add finite delta or do Taylor series expansion of numerator + denominator?
           !if (abs(dene) < tol12) then ??
           inv_dene2 = one / (dene ** 2)

           ! Loop overt perturbations and accumulate.
           ! Need all perts on the same procs as we have to take the outer product (ipc1, ipc2)
           do my_ip1=1,gqk%my_npert
             ipc1 = gqk%my_iperts(my_ip1)
             do my_ip2=1,gqk%my_npert
               ipc2 = gqk%my_iperts(my_ip2)
               ! my_g(2, my_npert, nb, my_nq, nb, my_nk)
               gmat(ipc1, ipc2, iq_ibz) = gmat(ipc1, ipc2, iq_ibz) &
                 + inv_dene2 * f_k_b2 * (one - f_kq_b1) &
                 ! TODO my_g should be complex
                 !* gqk%my_g(:,my_ip1,ib1,my_iq,ib2,my_ik) * gqk%my_g(:,my_ip2,ib1,my_iq,ib2,my_ik)
                 - inv_dene2 * f_kq_b1 * (one - f_k_b2)
                 !* gqk%my_g(:,my_ip1,ib1,my_iq,ib2,my_ik) * gqk%my_g(:,my_ip2,ib1,my_iq,ib2,my_ik)
             end do
           end do
         end do
       end do

     end do ! my_ik
   end do ! my_iq
   ABI_FREE(my_kqmap)
 end do ! spin

 ! ALL reduce partial terms so that we sum over MPI-distributed dims i.e. spins and k-points in BZ.
 gmat = j_dpc * gmat / gstore%nkbz
 ! Account for spin degeneracy.
 if (nsppol == 1 .and. dtset%nspinor == 1) gmat = two * gmat
 call xmpi_sum(gmat, gstore%comm, ierr)
 call cwtime_report(" berry_curvature:", cpu, wall, gflops)

 if (my_rank == 0) then
   ! Print some results to ab_out for testing purposes.
   do iq_ibz=1,gstore%nqibz
     write(ab_out, "(4a)")ch10,ch10," G(q) matrix for qpt: ", trim(ktoa(gstore%qibz(:,iq_ibz)))
     write(ab_out, "(4(i4,2x), 2(a12,2x))")"idir1", "ipert1", "idir2", "ipert2", "Re(gmat)", "Im(gmat)"
     do ipc2=1,natom3
       idir2 = mod(ipc2-1, 3) + 1; ipert2 = (ipc2 - idir2) / 3 + 1
       do ipc1=1,natom3
         idir1 = mod(ipc1-1, 3) + 1; ipert1 = (ipc1 - idir1) / 3 + 1
         write(ab_out, "(4(i4,2x), 2(es12.5,2x))") &
           idir1, ipert1, idir2, ipert2, real(gmat(ipc1, ipc2, iq_ibz)), aimag(gmat(ipc1, ipc2, iq_ibz))
       end do
     end do
   end do ! iq_ibz
 end if

 ! Create berry_ddb object to add gmat to the input dynamical matrix.
 ! WARNING: gmat has all (3*natom)**2 entries whereas the DDB stores only the irred elements.
 ! on the one hand, one can eploit this to reduce the number of iterations.
 ! on the other hand, one assumes that D(q) and G(q) have the same symmetry.
 ! Note that ddb%val(2,msize,nblok) are in Cartesian coordinates whereas gmat is in reduced coords.
 !
 ! Symmetry properties:
 !
 ! 1) G(-k) = -G(k)^*
 ! 2) In the presence of time reversal symmetry, the Berry curvature would be zero

 call in_ddb%copy(berry_ddb)

 ABI_MALLOC(d2matr, (2,3,in_ddb%mpert,3,in_ddb%mpert))
 !ABI_MALLOC(flg, (3,in_ddb%mpert,3,in_ddb%mpert))
 !ABI_FREE(flg)

 do iblok=1,in_ddb%nblok
   ddb_qq(:) = in_ddb%qpt(:,iblok) / in_ddb%nrm(:,iblok)
   ! In principle the qpts in the DDB and in gstore%qibz should be the same.
   !iq_ibz = gstore%find_qibz(ddb_qq)
   !ABI_CHECK(iq_ibz /= -1, sjoin("Cannot find qpt:", ktoa(ddb_qq), "in gstore%qibz !"))

   !do iq=1,db%nqpt
   !  if (all(abs(db%qpts(:, iq) - qpt) < my_qtol)) then
   !    iqpt = iq; exit
   !  end if
   !end do

   !d2matr = reshape(in_ddb%val(:,:,iblok), [2,3,in_ddb%mpert,3,in_ddb%mpert])
   !gmat(natom3, natom3, iq_ibz))
   !call berry_ddb%set_d2matr(d2matr, in_ddb%flg(:,iblok), iblok)
 end do

 ABI_FREE(d2matr)
 ABI_FREE(gmat)

 ! Write new DDB file with Berry curvature.
 if (my_rank == 0) then
   ddb_filepath = strcat(dtfil%filnam_ds(4), "BERRY_DDB")
   call wrtout(units, "- Writing DDB file with Berry curvature to: ", ddb_filepath)
   call berry_ddb%write_txt(in_ddb_hdr, ddb_filepath) !, fullinit, comm)
 end if

end subroutine berry_curvature
!!***

end module m_berry_curvature
!!***
