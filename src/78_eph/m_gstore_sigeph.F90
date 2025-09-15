!!****m* ABINIT/m_gstore
!! NAME
!! m_gstore
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2025 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gstore_sigeph

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 !use m_htetra
 !use libtetrabz
 use netcdf
 use m_nctk
 !use m_ddb
 !use m_dvdb
 use m_crystal,       only : crystal_t
 !use m_fft
 !use m_hamiltonian
 !use m_pawcprj
 use m_dtset,          only : dataset_type
 use m_dtfil,          only : datafiles_type
 !use m_wfd
 !use m_ephtk
 !use m_mkffnl
 !use m_sigtk

 !use defs_abitypes,    only : mpi_type
 !use m_time,           only : cwtime, cwtime_report, sec2str
 use m_fstrings,       only : tolower, itoa, ftoa, sjoin, ktoa, ltoa, strcat, replace_ch0, yesno, string_in
 !use m_cgtools,        only : cg_zdotc
 !use m_kg,             only : getph
 !use defs_datatypes,   only : pseudopotential_type
 !use m_hdr,            only : hdr_type, fform_from_ext
 use m_ebands, only : ebands_t, gaps_t
 !use m_matrix,         only : matr3inv
 use m_kpts,           only : kpts_timrev_from_kptopt, kpts_map !, kpts_sort, kpts_pack_in_stars, kptrlatt_from_ngkpt
 !use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack
 use m_ifc,            only : ifc_type
 !use m_phonons,        only : pheigvec_rotate
 !use m_pawang,         only : pawang_type
 !use m_pawrad,         only : pawrad_type
 !use m_pawtab,         only : pawtab_type
 !use m_pawfgr,         only : pawfgr_type
 !use m_pstat,          only : pstat_proc
 use m_occ,             only : occ_be, occ_fd
 use m_lgroup,          only : lgroup_t
 use m_gstore,          only : gstore_t !, gqk_t

 implicit none

 private

!!***

!----------------------------------------------------------------------

!!****t* m_gstore_sigeph/seph_t
!! NAME
!! seph_t
!!
!! FUNCTION
!!
!! NOTES
!!
!! SOURCE

type, public :: seph_t

  integer :: nsppol
   ! Number of independent spin polarizations.

!contains

end type seph_t
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_gstore/gstore_sigeph
!! NAME
!!  gstore_sigeph
!!
!! FUNCTION
!!  Compute matrix elements of the e-ph self-energy (Fan Migdal + Debye Waller).
!!  using precomputed e-ph matrix elements.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_sigeph(dtset, dtfil, cryst, ebands, ifc, comm)

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(ifc_type),target,intent(in) :: ifc
 integer,intent(in) :: comm

!Local variables-------------------------------
 integer,parameter :: master = 0, with_cplex1 = 1
 integer :: spin, my_is, my_ik, my_iq, my_ip, in_k, im_kq, ierr, ntemp, gap_err, my_rank
 integer :: it, ik_ibz, ikq_ibz, band_k, band_kq, timrev, ii, glob_ik
 real(dp) :: wqnu, gkq2, weight_qq
 real(dp) :: eig0nk, eig0mk, eig0mkq ! gdw2, gdw2_stern, rtmp !,nqnu,gkq2,gkq2_pf,
 !real(dp) :: cpu, wall, gflops
 logical :: use_lgk
 complex(dpc) :: ieta !, sig_cplx
 type(gaps_t) :: gaps
 type(lgroup_t) :: lg_myk
 type(gstore_t) :: gstore
!arrays
 integer :: units(2), my_kqmap(6)
 real(dp) :: kk(3), qpt(3), kq(3)
 real(dp),allocatable :: mu_e(:), kTmesh(:), nqnu(:), f_mkq(:), f_nk(:)
 complex(dpc),allocatable :: cfact(:)
 !real(dp),allocatable :: g2_pmnk(:,:,:,:)
 complex(dpc),allocatable :: vals_e0ks(:,:,:)
!----------------------------------------------------------------------

 units = [std_out, ab_out]
 my_rank = xmpi_comm_rank(comm)

 call wrtout(std_out, " Computing Fan-Migdal + DW self-energy from GSTORE.nc", pre_newlines=1)

 call gstore%from_ncpath(dtfil%filgstorein, with_cplex1, dtset, cryst, ebands, ifc, comm)
 call gstore%free()

 ABI_CHECK(gstore%qzone == "bz", "gstore_sigeph assumes qzone == `bz`")
 if (gstore%check_cplex_qkzone_gmode(1, "bz", "ibz", "phonon") /= 0) then
   ABI_ERROR("The gstore object is inconsistent with gstore_sigeph. See messages above.")
 end if

 ! Build (linear) mesh of K * temperatures. tsmesh(1:3) = [start, step, num]
 call dtset%get_ktmesh(ntemp, kTmesh)

 ! Compute the chemical potential at the different physical temperatures with Fermi-Dirac.
 ! TODO: One should check that nband is > nbocc to avoid inaccuracies in mu_e.
 ABI_MALLOC(mu_e, (ntemp))
 mu_e(:) = ebands%fermie
 if (dtset%eph_fermie == zero) then
   call ebands%get_muT_with_fd(ntemp, ktmesh, dtset%spinmagntarget, dtset%prtvol, mu_e, gstore%comm)
 end if

 gaps = ebands%get_gaps(gap_err)
 if (gap_err /= 0) then
   ABI_ERROR("Cannot compute fundamental and direct gap (likely metal)")
 end if

 if (my_rank == master) then
   call gaps%print(units, kTmesh=ktmesh, mu_e=mu_e, header="Gaps, band edges and relative position wrt Fermi level")
 end if
 call gaps%free()

 ABI_MALLOC(nqnu, (ntemp))
 ABI_MALLOC(f_nk, (ntemp))
 ABI_MALLOC(f_mkq, (ntemp))
 ABI_MALLOC(cfact, (ntemp))

 ieta = + j_dpc * dtset%zcut
 !use_lgk = dtset%symsigma /= 0
 use_lgk = .True.

 ! Loop over collinear spins.
 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is), cryst => gstore%cryst)
   spin = gstore%my_spins(my_is)

   ! vals_e0ks(ntemp, max_nbcalc, nk_glob, nsppol)
   ! Sigma_eph(omega=eKS, kT, band) for given (ikcalc, spin).
   ! Fan-Migdal + Debye-Waller
   ABI_CALLOC(vals_e0ks, (ntemp, gqk%nb, gqk%glob_nk))  ! nb_k
   ABI_FREE(vals_e0ks)

   ABI_CHECK(allocated(gqk%my_g2), "my_g2 is not allocated")
   ABI_CHECK(allocated(gqk%my_wnuq), "my_wnuq is not allocated")

   ! Loop over k-points in |n,k>
   do my_ik=1,gqk%my_nk
     kk = gqk%my_kpts(:, my_ik); ik_ibz = gqk%my_k2ibz(1, my_ik)
     ! store the results in glob_ik
     glob_ik = gqk%my_k2glob(my_ik)

     ! Compute the little group of the k-point so that we can compute g(k,q) only for q in the IBZ_k
     timrev = kpts_timrev_from_kptopt(ebands%kptopt)
     if (use_lgk) then
       call lg_myk%init(cryst, kk, timrev, gstore%nqbz, gstore%qbz, gstore%nqibz, gstore%qibz, xmpi_comm_self)
     end if

     ! Sum over q-points.
     do my_iq=1,gqk%my_nq
       call gqk%myqpt(my_iq, gstore, weight_qq, qpt)

       ! weight_qq is computed here. It depends whether we are assuming over the full BZ or IBZ_k.
       weight_qq = one / gstore%nqbz
       if (use_lgk) then
         ii = lg_myk%findq_ibzk(qpt)
         if (ii == -1) cycle
         weight_qq = lg_myk%weights(ii)
       end if

       ! Find the image of kq in the IBZ.
       kq = kk + qpt
       if (kpts_map("symrel", ebands%kptopt, cryst, gstore%krank_ibz, 1, kq, my_kqmap) /= 0) then
         ABI_ERROR(sjoin("Cannot map k+q to IBZ with qpt:", ktoa(qpt)))
       end if
       ikq_ibz = my_kqmap(1)

       ! Sum over phonon modes.
       do my_ip=1,gqk%my_npert
         wqnu = gqk%my_wnuq(my_ip, my_iq)
         nqnu(:) = occ_be(wqnu, kTmesh, zero)

         ! Copy data to improve memory access in the loops below.
         ! (my_npert, nb, my_nq, nb, my_nk)
         !g2_pmnk = gqk%my_g2(my_ip,:,my_iq,:,my_ik)

         ! Sum over bands.
         do im_kq=1,gqk%nb ! gqk%nb_kq
           band_kq = im_kq - gqk%bstart + 1
           eig0mkq = ebands%eig(band_kq, ikq_ibz, spin)

           ! Compute electronic occ for all Temps (note mu_e(it) Fermi level)
           do it=1,ntemp
             f_mkq(it) = occ_fd(eig0mkq, kTmesh(it), mu_e(it))
           end do

           ! Loop over the n index in |n,k>.
           do in_k=1,gqk%nb ! gqk%nb_k
             band_k = in_k - gqk%bstart + 1
             eig0nk = ebands%eig(band_k, ik_ibz, spin)

             ! Compute electronic occ for all Temps (note mu_e(it) Fermi level)
             !do it=1,ntemp
             !  f_nk(it) = occ_fd(eig0nk, kTmesh(it), mu_e(it))
             !end do ! it

             if (dtset%eph_ahc_type == 1) then
               cfact(:) =  (nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + ieta) + &
                           (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + ieta)
             else
               cfact(:) =  (two * nqnu + one) / (eig0nk - eig0mkq + ieta)
             end if

             ! Re + Im of Fan-Migdal self-energy
             ! (my_npert, nb_kq, my_nq, nb_k, my_nk)
             gkq2 = gqk%my_g2(my_ip, im_kq, my_iq, in_k, my_ik)
             cfact = cfact * gkq2 * weight_qq
             !se%vals_e0ks(:, in_k, ik_glob) = se%vals_e0ks(:, in_k, ik_glob) + cfact
           end do ! in_k
         end do ! my_ik
       end do ! im_kq
     end do ! my_iq
   end do ! my_ik
   call lg_myk%free()
   !call xmpi_sum(foobar, gstore%comm, ierr)
   end associate
 end do ! my_is

 ABI_FREE(mu_e)
 ABI_FREE(kTmesh)
 ABI_FREE(nqnu)
 ABI_FREE(f_nk)
 ABI_FREE(f_mkq)
 ABI_FREE(cfact)

end subroutine gstore_sigeph
!!***

end module m_gstore_sigeph
!!***
