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
 use m_crystal
 use m_dtset
 use m_dtfil
 use netcdf
 use m_nctk

 use m_time,            only : cwtime, cwtime_report
 use m_fstrings,        only : strcat, sjoin, ktoa
 use defs_datatypes,    only : ebands_t
 use m_kpts,            only : kpts_timrev_from_kptopt, kpts_map
 use m_ddb_hdr,         only : ddb_hdr_type
 use m_ddb,             only : ddb_type
 use m_ifc,             only : ifc_type
 use m_gstore,          only : gstore_t, gqk_t
 use m_phonons,         only : phdos_t, ifc_mkphbs
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
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine berry_curvature(gstore, dtset, dtfil, in_ddb, in_ifc, dielt, zeff, qdrp_cart)

!Arguments ------------------------------------
!scalars
 type(gstore_t),target,intent(inout) :: gstore
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(ddb_type),intent(in) :: in_ddb
 type(ifc_type),intent(in) :: in_ifc
 real(dp),intent(in) :: dielt(3,3), zeff(3,3,dtset%natom), qdrp_cart(3,3,3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0, nsphere0 = 0, prtsrlr0 = 0, rftyp1 = 1
 integer :: nproc, my_rank, ierr, itemp, ntemp, comm, ddb_nqshift, ncid, mpert, msize
 integer :: my_is, spin, nsppol, natom3, band, ib1, ib2, band1, band2, nb, ebands_timrev, index
 integer :: ik_ibz, isym_k, trev_k, tsign_k, g0_k(3)
 integer :: ikq_ibz, isym_kq, trev_kq, tsign_kq, g0_kq(3)
 integer :: iq_ibz, isym_q, trev_q, tsign_q, g0_q(3)
 integer :: my_iq, iq_glob, my_ik, ik_glob, my_ip1, my_ip2, ipc1, ipc2, ipert1, ipert2, iblok, idir1, idir2
 logical :: isirr_k, isirr_q, isirr_kq
 real(dp) :: e_b1_k, e_b2_k, e_b1_kq, e_b2_kq, f_b1_k, f_b1_kq, f_b2_k, f_b2_kq, dene, spin_occ, fact(2)
 real(dp) :: kt, wmax, cpu, wall, gflops
 character(len=5000) :: msg
 character(len=fnlen) :: berry_ddb_filepath, path
 class(crystal_t),pointer :: cryst
 class(ebands_t),pointer :: ebands
 type(gqk_t),pointer :: gqk
 type(crystal_t) :: ddb_crystal
 type(ddb_type) :: berry_ddb
 type(ddb_hdr_type) :: in_ddb_hdr
 type(ifc_type) :: berry_ifc
 type(phdos_t) :: phdos
!arrays
 integer :: count_wminmax(2), units(2), rfphon(4),rfelfd(4),rfstrs(4)
 integer,allocatable :: my_kqmap(:,:)
 real(dp) :: qphnrm(3), qphon_padded(3,3), g_ri(2), wminmax(2)
 real(dp) :: kk_bz(3), kq_bz(3), kk_ibz(3), kq_ibz(3), qq_bz(3), qq_ibz(3) !, vk(3)
 real(dp),allocatable :: ddb_qshifts(:,:)
 !real(dp),allocatable :: ktmesh(:), phmesh(:), eig_k(:,:), eig_kq(:,:), kmesh(:,:), wght_bz(:,:,:)
 complex(dp),allocatable :: gmat(:,:,:)

!----------------------------------------------------------------------

 comm = gstore%comm; nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
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
       ! + <b1,k|D_{-q,p1}|b2,k+q> <b2,k+q|D_{q,p2}|b2,k> / (e_{b1,k} - e_{b2,k+q})^2  <<< term 1
       ! - <b2,k|D_{-q,p1}|b1,k+q> <b1,k+q|D_{q,p2}|b2,k> / (e_{b2,k} - e_{b1,k+q})^2  <<< term 2
       !
       ! where b2 is the initial band and b1 the final band, p1, p2 are atomic perturbations in reduced coords
       ! and we're summing over k at fixed q.
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

           fact(1) = f_b1_k  * (one - f_b2_kq)
           fact(2) = f_b1_kq * (one - f_b2_k)
           if (all(abs(fact) < tol20)) cycle

           ! TODO: add finite delta or do Taylor series expansion of numerator + denominator?
           dene = e_b1_k - e_b2_kq
           if (abs(dene) > tol8) then
             fact(1) = fact(1) / dene**2
           else
             fact(1) = zero
           end if
           dene = e_b2_k - e_b1_kq
           if (abs(dene) > tol8) then
             fact(2) = fact(2) / dene**2
           else
             fact(2) = zero
           end if

           ! Loop over perturbations and accumulate.
           do my_ip1=1,gqk%my_npert
             ipc1 = gqk%my_iperts(my_ip1)
             do my_ip2=1,gqk%my_npert
               ipc2 = gqk%my_iperts(my_ip2)
               ! my_g(my_npert, nb, my_nq, nb, my_nk)
               gmat(ipc1, ipc2, iq_ibz) = gmat(ipc1, ipc2, iq_ibz) &
                + fact(1) * conjg(gqk%my_g(my_ip1,ib1,my_iq,ib2,my_ik)) * gqk%my_g(my_ip2,ib1,my_iq,ib2,my_ik) &
                - fact(2) * conjg(gqk%my_g(my_ip1,ib2,my_iq,ib1,my_ik)) * gqk%my_g(my_ip2,ib2,my_iq,ib1,my_ik)
             end do
           end do
         end do
       end do

     end do ! my_ik
   end do ! my_iq
   ABI_FREE(my_kqmap)
 end do ! spin

 ! Here we ALL_REDUCE all partial integrals (sum over MPI-distributed dims i.e. spins and k-points in BZ).
 ! Also, account for spin degeneracy as f in [0,1] if collinear.
 gmat = j_dpc * gmat / gstore%nkbz
 if (nsppol == 1 .and. dtset%nspinor == 1) gmat = two * gmat
 call xmpi_sum(gmat, comm, ierr)
 call massmult_and_breaksym_cplx(cryst%natom, cryst%ntypat, cryst%typat, in_ifc%amu, gmat, herm_opt=0)
 call cwtime_report(" berry_curvature:", cpu, wall, gflops)

 if (my_rank == master) then
   ! Print some results to ab_out for testing purposes.
   do iq_ibz=1,gstore%nqibz
     write(msg, "(2a,2x,2a)")ch10,ch10,"G(q) matrix for qpt: ", trim(ktoa(gstore%qibz(:,iq_ibz)))
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
 stop

 ! Create berry_ddb object to add gmat to the input dynamical matrix.
 ! WARNING: gmat has all (3*natom)**2 entries whereas the DDB stores only the irred elements.
 ! on the one hand, one can eploit this to reduce the number of iterations.
 ! on the other hand, one assumes that D(q) and G(q) have the same symmetry.
 ! Note that ddb%val(2,msize,nblok) are in Cartesian coordinates whereas gmat is in reduced coords.
 ! thus we need to use raw=1 to disable any symetrization or transformation to cartesian coordinates.
 !
 ! Symmetry properties:
 !
 ! 1) G(-k) = -G(k)^*
 ! 2) In the presence of time reversal symmetry, the Berry curvature would be zero

 call wrtout(units, sjoin("- Reading input DDB from file:", dtfil%filddbsin))
 call berry_ddb%from_file_txt(dtfil%filddbsin, dtset%brav, in_ddb_hdr, ddb_crystal, comm, prtvol=dtset%prtvol, raw=1)
 call ddb_crystal%free()

 ! mpert = maximum number of perturbations (atom displacements + electric field + ...)
 ! msize = maximum size of one block of the ddb e.g. 3*mpert * 3*mpert.
 mpert = berry_ddb%mpert; msize = berry_ddb%msize
 rfphon(1:2)=1; rfelfd(1:2)=0; rfstrs(1:2)=0; qphnrm = one

 do iq_ibz=1,gstore%nqibz
   ! Find 2d block associated to gstore%qibz(:,iq_ibz)
   qphon_padded = zero; qphon_padded(:,1) = gstore%qibz(:,iq_ibz)
   call berry_ddb%get_block(iblok, qphon_padded, qphnrm, rfphon, rfelfd, rfstrs, rftyp1)
   ABI_CHECK(iblok /= 0, sjoin("Cannot find q-point ", ktoa(gstore%qibz(:,iq_ibz))," in DDB file"))

   do ipc2=1,natom3
     idir2 = mod(ipc2-1, 3) + 1; ipert2 = (ipc2 - idir2) / 3 + 1
     do ipc1=1,natom3
       idir1 = mod(ipc1-1, 3) + 1; ipert1 = (ipc1 - idir1) / 3 + 1
       index = idir1+3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))
       ! TODO: In MGO all flags are set to 1
       print *, "flg:", berry_ddb%flg(index, iblok)
       !if (berry_ddb%flg(index, iblok) == 0) cycle
       g_ri(:) = [real(gmat(ipc1, ipc2, iq_ibz)), aimag(gmat(ipc1, ipc2, iq_ibz))]
       !berry_ddb%val(:,index,iblok) = berry_ddb%val(:,index,iblok) + g_ri
     end do
   end do
 end do

 ABI_FREE(gmat)

 ! Write new DDB file including Berry curvature.
 berry_ddb_filepath = strcat(dtfil%filnam_ds(4), "_BERRY_DDB")
 if (my_rank == master) then
   call wrtout(units, sjoin("- Writing DDB file with Berry curvature to: ", berry_ddb_filepath), pre_newlines=2)
   in_ddb_hdr%dscrpt = "DDB including Berry curvature"
   call berry_ddb%write_txt(in_ddb_hdr, berry_ddb_filepath) !, fullinit, comm)
 end if
 call in_ddb_hdr%free(); call berry_ddb%free()

 ! ===============================
 ! Now reopen berry_ddb with raw=0
 ! ===============================
 call xmpi_barrier(comm)
 call berry_ddb%from_file_txt(berry_ddb_filepath, dtset%brav, in_ddb_hdr, ddb_crystal, comm, &
                              prtvol=dtset%prtvol, raw=0)
 call in_ddb_hdr%free(); call ddb_crystal%free()

 ! Build berry_ifc from berry_ddb.
 ! Set the q-shift for the DDB (well we mainly use gamma-centered q-meshes)
 ddb_nqshift = 1
 ABI_CALLOC(ddb_qshifts, (3, ddb_nqshift))
 ddb_qshifts(:,1) = dtset%ddb_shiftq(:)

 call berry_ifc%init(cryst, berry_ddb, &
   dtset%brav, dtset%asr, dtset%symdynmat, dtset%dipdip, dtset%rfmeth, &
   dtset%ddb_ngqpt, ddb_nqshift, ddb_qshifts, dielt, zeff, &
   qdrp_cart, nsphere0, dtset%rifcsph, prtsrlr0, dtset%enunit, comm, &
   dipquad=dtset%dipquad, quadquad=dtset%quadquad)
 !call berry_ifc%print(unit=std_out)

 ABI_FREE(ddb_qshifts)

 ! Output phonon band structure (requires qpath)
 ! TODO: Change prefix to encode berry curvature?
 if (dtset%prtphbands /= 0) call ifc_mkphbs(berry_ifc, cryst, dtset, dtfil%filnam_ds(4), comm)

 if (dtset%prtphdos == 1) then
   call wrtout(std_out, " Computing Phonon DOS. Use prtphdos 0 to disable this part.")
   wminmax = zero
   do
     call phdos%init(cryst, berry_ifc, dtset%ph_intmeth, dtset%ph_wstep, dtset%ph_smear, dtset%ph_ngqpt, &
                     dtset%ph_nqshift, dtset%ph_qshift, "", wminmax, count_wminmax, comm)
     if (all(count_wminmax == 0)) exit
     wminmax(1) = wminmax(1) - abs(wminmax(1)) * 0.05; wminmax(2) = wminmax(2) + abs(wminmax(2)) * 0.05
     call phdos%free()
     write(msg, "(a, 2f8.5)") "Initial frequency mesh not large enough. Recomputing PHDOS with wmin, wmax: ",wminmax
     call wrtout(std_out, msg)
   end do

   if (my_rank == master) then
     path = strcat(dtfil%filnam_ds(4), "_PHDOS.nc")
     call wrtout(units, sjoin("- Writing phonon DOS to netcdf file:", path))
     NCF_CHECK_MSG(nctk_open_create(ncid, path, xmpi_comm_self), sjoin("Creating PHDOS.nc file:", path))
     NCF_CHECK(cryst%ncwrite(ncid))
     call phdos%ncwrite(ncid)
     NCF_CHECK(nf90_close(ncid))
   end if
   call phdos%free()
 end if ! prtphdos

 call berry_ifc%free(); call berry_ddb%free()

end subroutine berry_curvature
!!***

end module m_berry_curvature
!!***
