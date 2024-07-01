!!****m* ABINIT/m_berry_curvature
!! NAME
!! m_berry_curvature
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2024 ABINIT group (MG, MMignolet)
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
 use m_crystal
 use m_dtset
 use m_dtfil
 use netcdf
 use m_nctk

 use m_time,            only : cwtime, cwtime_report
 use m_fstrings,        only : strcat, sjoin, ktoa
 use defs_datatypes,    only : ebands_t
 use m_kpts,            only : kpts_timrev_from_kptopt, kpts_map, smpbz
 use m_ddb_hdr,         only : ddb_hdr_type, BLKTYP_d2E_mbc
 use m_ddb,             only : ddb_type
 use m_ifc,             only : ifc_type
 use m_gstore,          only : gstore_t, gqk_t
 use m_dynmat,          only : massmult_and_breaksym_cplx

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
!! Computes the molecular Berry curvature and writes the reulst to a new ddb
!! Ref: see D. Saparov PRB 105, 064303 (2022)
!!
!! INPUTS
!! gstore<type(gstore_t)> = gstore object with el-ph matrix element computed
!!   the whole BZ
!! dtset<type(dataset_type)> even though its listed as inout, its not modified
!! dtfil<type(datafiles_type)>
!! ! in_ddb<type(ddb_type)>
!! in_ifc<type(ifc_type)>
!! dielt(3,3)=dielectric tensor.
!! zeff(3,3,natom)=effective charge on each atom, versus electric field and atomic displacement
!! qdrp_cart(3,3,3,natom)=Quadrupole tensor on each atom in cartesian cordinates
!!
!! OUTPUT
!! Generates a new ddb file with the molecular Berry curvature: *_BERRY_DDB
!!
!! SOURCE

subroutine berry_curvature(gstore, dtset, dtfil)

!Arguments ------------------------------------
!scalars
 type(gstore_t),target,intent(inout) :: gstore
 type(dataset_type),intent(inout) :: dtset
 type(datafiles_type),intent(in) :: dtfil

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0, LOG_MODQ = 5
 integer :: nproc, my_rank, ierr, comm, mpert, msize
 integer :: my_is, spin, nsppol, ntypat, natom, natom3, ib1, ib2, band1, band2, nb, ebands_timrev
 integer :: ik_ibz, isym_k, trev_k, tsign_k, g0_k(3)
 integer :: ikq_ibz, isym_kq, trev_kq, tsign_kq, g0_kq(3)
 integer :: iq_ibz, isym_q, trev_q, tsign_q, g0_q(3)
 integer :: my_iq, iq_glob, my_ik, ik_glob
 integer :: my_ip1, my_ip2, ipc1, ipc2, ipert1, ipert2, nblock, idir1, idir2
 logical :: isirr_k, isirr_q, isirr_kq, print_qtime
 real(dp) :: e_b1_k, e_b2_k, e_b1_kq, e_b2_kq, f_b1_k, f_b1_kq, f_b2_k, f_b2_kq, dene, spin_occ, fact(2)
 real(dp) :: cpu_all, wall_all, gflops_all, cpu_q, wall_q, gflops_q
 complex(dp) :: tmp, tmp1, tmp2
 character(len=5000) :: msg
 character(len=fnlen) :: berry_ddb_filepath
 class(crystal_t),pointer :: cryst
 class(ebands_t),pointer :: ebands
 type(gqk_t),pointer :: gqk
 type(crystal_t) :: in_ddb_crystal
 type(ddb_type) :: in_ddb, berry_ddb
 type(ddb_hdr_type) :: in_ddb_hdr, berry_ddb_hdr
!arrays
 integer :: units(2)
 integer,allocatable :: blkflg(:,:,:,:)
 integer,allocatable :: my_kqmap(:,:)
 real(dp) :: qphon(3)
 real(dp) :: qq_ibz(3)
 real(dp),allocatable :: tmp_mat(:,:,:,:,:)
 complex(dp),allocatable :: gmat(:,:,:)
!----------------------------------------------------------------------

 comm = gstore%comm; nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 units = [std_out, ab_out]
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 cryst => gstore%cryst; ebands => gstore%ebands
 natom = cryst%natom; ntypat = cryst%ntypat
 natom3 = 3 * cryst%natom; nsppol = ebands%nsppol
 spin_occ = one; if (nsppol == 1 .and. dtset%nspinor == 1) spin_occ = two
 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)

 if (my_rank == master) then
   call wrtout(std_out, " Computing berry curvature", pre_newlines=2)
   call gstore%print(std_out, header="Gstore", prtvol=dtset%prtvol)
 end if

 ! Consistency check
 ierr = 0
 ABI_CHECK_NOSTOP(gstore%qzone == "ibz", "qzone == 'ibz' is required", ierr)
 ABI_CHECK_NOSTOP(gstore%kzone == "bz", "kzone == 'bz' is required", ierr)
 ABI_CHECK_NOSTOP(gstore%gqk(1)%cplex == 2, "cplex == 2 is required", ierr)
 ! Need all perts on the same procs as we have to take the outer product (ipc1, ipc2)
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
     print_qtime = (my_iq <= LOG_MODQ .or. mod(my_iq, LOG_MODQ) == 0)
     if (print_qtime) call cwtime(cpu_q, wall_q, gflops_q, "start")
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

       ! Summation over the two band indices. NB: occupancies f are rescaled in [0, 1] when nsppol 2.
       ! Here we accumulate:
       !
       ! + <b1,k|D_{-q,p1}|b2,k+q> <b2,k+q|D_{q,p2}|b1,k> / (e_{b1,k} - e_{b2,k+q})^2  <<< term 1
       ! - <b2,k|D_{-q,p1}|b1,k+q> <b1,k+q|D_{q,p2}|b2,k> / (e_{b2,k} - e_{b1,k+q})^2  <<< term 2
       !
       ! where b1 is the initial band and b2 the final band, p1, p2 are atomic perturbations in reduced coords
       ! and we're summing over k in the BZ at fixed q.
       do ib2=1,nb
         band2 = ib2 + gqk%bstart - 1
         e_b2_k = ebands%eig(band2, ik_ibz, spin)
         f_b2_k = ebands%occ(band2, ik_ibz, spin) / spin_occ
         e_b2_kq = ebands%eig(band2, ikq_ibz, spin)
         f_b2_kq = ebands%occ(band2, ikq_ibz, spin) / spin_occ

         do ib1=1,nb
           band1 = ib1 + gqk%bstart - 1
           e_b1_kq = ebands%eig(band1, ikq_ibz, spin)
           f_b1_kq = ebands%occ(band1, ikq_ibz, spin) / spin_occ
           e_b1_k = ebands%eig(band1, ik_ibz, spin)
           f_b1_k = ebands%occ(band1, ik_ibz, spin) / spin_occ

           ! following the formula: m = b1 and mprime = b2
           fact(1) = f_b1_k  * (one - f_b2_kq)
           fact(2) = f_b1_kq * (one - f_b2_k)

           if (all(abs(fact) < tol20)) cycle

           dene = e_b1_k - e_b2_kq
           if (abs(dene) > tol12) then
           ! the tolerance here might need some tweaking
             fact(1) = fact(1) / dene**2
           else
             ! TODO: add finite delta or do Taylor series expansion of numerator + denominator?
             fact(1) = zero
           end if

           dene = e_b2_k - e_b1_kq
           if (abs(dene) > tol12) then
           ! the tolerance here might need some tweaking
             fact(2) = fact(2) / dene**2
           else
             ! TODO: add finite delta or do Taylor series expansion of numerator + denominator?
             fact(2) = zero
           end if

           ! Loop over perturbations and accumulate.
           do my_ip1=1,gqk%my_npert
             ipc1 = gqk%my_iperts(my_ip1)
             do my_ip2=1,gqk%my_npert
               ipc2 = gqk%my_iperts(my_ip2)
               ! my_g(my_npert, nb, my_nq, nb, my_nk)

               ! 1st term
               ! <k, b1| D_{-q,p1}H |k+q, b2> * <k+q, b2| D_{q,p2}H |k, b1>
               tmp1 = fact(1) * conjg(gqk%my_g(my_ip1,ib2,my_iq,ib1,my_ik)) * gqk%my_g(my_ip2,ib2,my_iq,ib1,my_ik)

               ! 2nd term
               ! <k, b2| D_{-q,p1}H |k+q, b1> * <k+q, b1| D_{q,p2}H |k, b2>
               tmp2 = fact(2) * conjg(gqk%my_g(my_ip1,ib1,my_iq,ib2,my_ik)) * gqk%my_g(my_ip2,ib1,my_iq,ib2,my_ik)

               tmp = tmp1 - tmp2
               gmat(ipc1, ipc2, iq_ibz) = gmat(ipc1, ipc2, iq_ibz) + tmp
             end do
           end do

         end do ! ib1
       end do ! ib2
     end do ! my_ik

     if (print_qtime) then
       write(msg,'(4x,3(a,i0),a)')"my_iq [", my_iq, "/", gqk%my_nq, "] (tot: ", gstore%nqibz, ")"
       call cwtime_report(msg, cpu_q, wall_q, gflops_q); if (my_iq == LOG_MODQ) call wrtout(std_out, " ...")
     end if
   end do ! my_iq
   ABI_FREE(my_kqmap)
 end do ! spin

 ! Here we ALL_REDUCE all partial integrals (sum over MPI-distributed dims i.e. spins and k-points in BZ).
 ! Also, account for spin degeneracy as f in [0,1] if collinear.
 gmat = j_dpc * gmat / gstore%nkbz
 if (nsppol == 1 .and. dtset%nspinor == 1) gmat = two * gmat
 call xmpi_sum(gmat, comm, ierr)
 call massmult_and_breaksym_cplx(cryst%natom, cryst%ntypat, cryst%typat, gstore%ifc%amu, gmat, herm_opt=0)
 call cwtime_report(" berry_curvature:", cpu_all, wall_all, gflops_all)

 if (my_rank == master) then
   ! Print some results to ab_out for testing purposes.
   do iq_ibz=1,gstore%nqibz
     write(msg, "(2a,2x,2a)")ch10, ch10, "G(q) matrix for qpt: ", trim(ktoa(gstore%qibz(:,iq_ibz)))
     call wrtout(units, msg)
     write(msg, "(2x,4(a6,2x), 2(a12,2x))")"idir1", "ipert1", "idir2", "ipert2", "Re(gmat)", "Im(gmat)"
     call wrtout(units, msg)
     do ipc2=1,natom3
       idir2 = mod(ipc2-1, 3) + 1; ipert2 = (ipc2 - idir2) / 3 + 1
       do ipc1=1,natom3
         idir1 = mod(ipc1-1, 3) + 1; ipert1 = (ipc1 - idir1) / 3 + 1
         write(msg, "(2x,4(i6,2x), 2(es12.5,2x))") &
           idir1, ipert1, idir2, ipert2, real(gmat(ipc1, ipc2, iq_ibz)), aimag(gmat(ipc1, ipc2, iq_ibz))
         call wrtout(units, msg)
       end do
     end do
   end do ! iq_ibz
 end if
 !stop

 ! Symmetry properties:
 !  1) G(-q) = G(q)^*
 !  2) In the presence of time reversal symmetry, the Berry curvature would be zero
 call wrtout(units, sjoin("- Reading input DDB from file:", dtfil%filddbsin))
 call in_ddb%from_file(dtfil%filddbsin, in_ddb_hdr, in_ddb_crystal, comm, dtset%prtvol)
 call in_ddb_crystal%free()

 ! Initialize ddb header object
 call wrtout(units, "- Initialize berry_ddb_hdr:")
 ! dirty trick otherwise ddb_hdr_init complains about not havinig kpt weights
 ABI_MALLOC(dtset%wtk, (gstore%nkbz))
 dtset%wtk(:) = one / (one * gstore%nkbz)
 call berry_ddb_hdr%init(dtset,in_ddb_hdr%psps, in_ddb_hdr%pawtab, &
                         dscrpt=' Molecular Berry curvature ', &
                         nblok=gstore%nqibz, nkpt=gstore%ebands%nkpt, kpt=gstore%ebands%kptns, &
                         occ=gstore%ebands%occ)
 ABI_FREE(dtset%wtk)
 call in_ddb_hdr%free()
 call in_ddb%free()

 !  Initialize ddb object
 ! nblock = number of block -> number of qpt here
 ! mpert = maximum number of perturbations (atom displacements + electric field + ...)
 ! msize = maximum size of one block of the ddb e.g. 3*mpert * 3*mpert.
 nblock = gstore%nqibz
 mpert = natom+6 ! in principle we only need natom, but some portions of anaddb
                         ! require at least natom+6
 msize = 3*mpert * 3*mpert
 call wrtout(units, "- Initialize berry_ddb:")
 call berry_ddb%init(dtset, nblock, mpert=mpert, with_d2E=.true.)

 ! Set the values for the 2nd order derivatives
 ABI_CALLOC(blkflg, (3, mpert, 3, mpert))
 ABI_CALLOC(tmp_mat, (2, 3, mpert, 3, mpert))
 blkflg(:3,:natom,:3,:natom) = 1 ! all values have been computed
 do iq_ibz=1,gstore%nqibz
   tmp_mat(1,:3,:natom,:3,:natom) = reshape( real(gmat(:,:,iq_ibz)), (/3,natom,3,natom/))
   tmp_mat(2,:3,:natom,:3,:natom) = reshape(aimag(gmat(:,:,iq_ibz)), (/3,natom,3,natom/))
   qphon = gstore%qibz(:,iq_ibz)
   call berry_ddb%set_qpt(iq_ibz, qphon)
   call berry_ddb%set_typ(iq_ibz, BLKTYP_d2E_mbc)
   call berry_ddb%set_d2matr(iq_ibz, tmp_mat, blkflg)
 end do

 ! Output Berry ddb
 if (my_rank == master) then
   berry_ddb_filepath = strcat(dtfil%filnam_ds(4), "_BERRY_DDB")
   call wrtout(units, sjoin("- Writing DDB file with Berry curvature to: ", berry_ddb_filepath), pre_newlines=2)
   call berry_ddb%write(berry_ddb_hdr, berry_ddb_filepath)
 end if

 ! Deallocate ddb object
 call berry_ddb_hdr%free()
 call berry_ddb%free()
 ABI_FREE(gmat)
 ABI_FREE(blkflg)
 ABI_FREE(tmp_mat)

end subroutine berry_curvature
!!***

end module m_berry_curvature
!!***
