!!****m* ABINIT/m_gstore
!! NAME
!! m_gstore
!!
!! FUNCTION
!!  Compute matrix elements of the e-ph self-energy (Fan Migdal + Debye Waller).
!!  using precomputed e-ph matrix elements.
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
 use m_crystal,        only : crystal_t
 !use m_hamiltonian
 use m_dtset,          only : dataset_type
 use m_dtfil,          only : datafiles_type
 !use m_wfd
 use m_ephtk
 !use m_mkffnl
 use m_sigtk

 !use defs_abitypes,    only : mpi_type
 !use m_time,           only : cwtime, cwtime_report, sec2str
 use m_fstrings,       only : tolower, itoa, ftoa, sjoin, ktoa, ltoa, strcat, replace_ch0, yesno, string_in
 !use m_cgtools,        only : cg_zdotc
 !use m_kg,             only : getph
 !use defs_datatypes,   only : pseudopotential_type
 !use m_hdr,            only : hdr_type, fform_from_ext
 use m_geometry,       only : phdispl_cart2red_nmodes
 use m_ebands,         only : ebands_t, gaps_t
 use m_kpts,           only : kpts_timrev_from_kptopt, kpts_map !, kpts_sort, kpts_pack_in_stars, kptrlatt_from_ngkpt
 !use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack
 use m_ifc,            only : ifc_type
 !use m_phonons,        only : pheigvec_rotate
 !use m_pawang,         only : pawang_type
 !use m_pawrad,         only : pawrad_type
 !use m_pawtab,         only : pawtab_type
 !use m_pawfgr,         only : pawfgr_type
 !use m_pstat,          only : pstat_proc
 use m_occ,            only : occ_be, occ_fd
 use m_lgroup,         only : lgroup_t
 use m_gstore,         only : gstore_t, gqk_t

 implicit none

 private
 public :: gstore_sigeph
!!***

!----------------------------------------------------------------------

!!****t* m_gstore_sigeph/seph_t
!! NAME
!! seph_t
!!
!! FUNCTION
!!
!! SOURCE

!! type, public :: seph_t
!!
!!   integer :: nsppol
!!    ! Number of independent spin polarizations.
!!
!! !contains
!!
!! end type seph_t
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
 integer,parameter :: master = 0, with_cplex1 = 1, max_ntemp = 50
 integer :: spin, my_is, my_ik, my_iq, my_ip, in_k, im_kq, ierr, ntemp, gap_err, my_rank
 integer :: it, ik_ibz, ikq_ibz, band_k, band_kq, timrev, ii, ikcalc, natom, natom3, nsppol
 real(dp) :: wqnu, gkq2, weight_q, eig0nk, eig0mk, eig0mkq, ediff, gmod2, hmod2, gdw2, gdw2_stern, rtmp !,nqnu,gkq2,gkq2_pf,
 !real(dp) :: cpu, wall, gflops
 logical :: use_lgk, q_is_gamma, intra_band, same_band, imag_only
 complex(dp) :: ieta !, sig_cplx
 type(gaps_t) :: gaps
 type(lgroup_t) :: lg_myk
 type(gstore_t) :: gstore
!arrays
 integer :: units(2), my_kqmap(6)
 real(dp) :: kk(3), qpt(3), kq(3)
 real(dp),allocatable :: gkq_atm(:,:,:),gkq_nu(:,:,:) !,gkq0_atm(:,:,:,:), gaussw_qnu(:)
 real(dp),allocatable :: rfact(:), mu_e(:), kTmesh(:), nqnu(:), f_mkq(:) !, f_nk(:),  g2_pmnk(:,:,:,:)
 real(dp),allocatable :: displ_nu_cart(:,:,:), displ_nu_red(:,:,:)
 complex(dp),allocatable :: cfact(:), vals_e0ks(:,:,:), dvals_de0ks(:,:,:), tpp_red(:,:) !, fmw_frohl_sphcorr(:,:,:,:), cfact_wr(:),
!----------------------------------------------------------------------

 units = [std_out, ab_out]
 my_rank = xmpi_comm_rank(comm)

 call wrtout(std_out, " Computing Fan-Migdal + DW self-energy from GSTORE.nc", pre_newlines=1)
 natom = cryst%natom; natom3 = 3 * cryst%natom; nsppol = gstore%nsppol
 imag_only = .False.

 ! To compute the DW term in the RIA, we need the complex g in the atom representation
 call gstore%from_ncpath(dtfil%filgstorein, with_cplex1, dtset, cryst, ebands, ifc, comm, &
             gvals_vname="gvals", read_dw=.True.)

 ABI_CHECK(gstore%qzone == "bz", "gstore_sigeph assumes qzone == `bz`")
 if (gstore%check_cplex_qkzone_gmode(1, "bz", "ibz", "phonon") /= 0) then
 !if (gstore%check_cplex_qkzone_gmode(with_cplex1, "bz", "ibz", "atom") /= 0) then
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
 !ABI_MALLOC(f_nk, (ntemp))
 ABI_MALLOC(f_mkq, (ntemp))
 ABI_MALLOC(cfact, (ntemp))
 ABI_MALLOC(rfact, (ntemp))

 ieta = + j_dpc * dtset%zcut
 !use_lgk = dtset%symsigma /= 0
 use_lgk = .True.
 use_lgk = .False.

 ! Compute electronic occ for all Temps (note mu_e(it) Fermi level)
 !do it=1,ntemp
 !  f_nk(it) = occ_fd(eig0nk, kTmesh(it), mu_e(it))
 !end do ! it

 ABI_MALLOC(displ_nu_cart, (2, 3, natom))
 ABI_MALLOC(displ_nu_red, (2, 3, natom))
 ABI_MALLOC(tpp_red, (natom3, natom3))

 ! Loop over collinear spins.
 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is), cryst => gstore%cryst)
   spin = gstore%my_spins(my_is)
   ABI_CHECK(allocated(gqk%my_g2), "my_g2 is not allocated")
   ABI_CHECK(allocated(gqk%my_wnuq), "my_wnuq is not allocated")

   ! vals_e0ks(ntemp, max_nbcalc, nk_glob, nsppol)
   ! Sigma_eph(omega=eKS, kT, band) for given (ikcalc, spin).
   ! Fan-Migdal + Debye-Waller
   ABI_CALLOC(vals_e0ks, (ntemp, gqk%nb, gqk%glob_nk))  ! nb_k
   ABI_CALLOC(dvals_de0ks, (ntemp, gqk%nb, gqk%glob_nk))  ! nb_k

   ! Loop over k-points in |n,k>
   do my_ik=1,gqk%my_nk
     kk = gqk%my_kpts(:, my_ik); ik_ibz = gqk%my_k2ibz(1, my_ik)
     ! will store results using glob_ik index.
     ikcalc = gqk%my_k2glob(my_ik)

     ! Compute the little group of the k-point so that we can compute g(k,q) only for q in the IBZ_k.
     if (use_lgk) then
        timrev = kpts_timrev_from_kptopt(ebands%kptopt)
       call lg_myk%init(cryst, kk, timrev, gstore%nqbz, gstore%qbz, gstore%nqibz, gstore%qibz, xmpi_comm_self)
     end if

     ! Sum over q-points.
     do my_iq=1,gqk%my_nq
       call gqk%myqpt(my_iq, gstore, weight_q, qpt)
       q_is_gamma = sum(qpt**2) < tol14

       ! weight_q is computed here. It depends whether we are assuming over the full BZ or IBZ_k.
       weight_q = one / gstore%nqbz
       if (use_lgk) then
         ii = lg_myk%findq_ibzk(qpt); if (ii == -1) cycle; weight_q = lg_myk%weights(ii)
       end if

       ! Find the image of kq in the IBZ.
       kq = kk + qpt
       if (kpts_map("symrel", ebands%kptopt, cryst, gstore%krank_ibz, 1, kq, my_kqmap) /= 0) then
         ABI_ERROR(sjoin("Cannot map k+q to IBZ with k+q:", ktoa(kq)))
       end if
       ikq_ibz = my_kqmap(1)

       ! Sum over phonon modes.
       do my_ip=1,gqk%my_npert
         wqnu = gqk%my_wnuq(my_ip, my_iq)
         nqnu(:) = occ_be(wqnu, kTmesh, zero)
         !print *, "wqnu", wqnu

         displ_nu_cart = gqk%my_displ_cart(:,:,:,my_ip,my_iq)
         call phdispl_cart2red_nmodes(natom, 1, cryst%gprimd, displ_nu_cart, displ_nu_red)

         ! For the DW, we need gkq_atm at q = 0, for all atomic perturbations

         !ABI_MALLOC(gkq_atm, (2, nbcalc_ks, natom3))
         !ABI_MALLOC(gkq_nu, (2, nbcalc_ks, natom3))
         !!gkq_atm = gqk%my_g(my_ip, im_kq, my_iq, in_k, my_ik)
         !call ephtk_gkknu_from_atm(1, nbcalc_ks, 1, natom, gkq_atm, phfrq, displ_red, gkq_nu)
         !ABI_FREE(gkq_atm)
         !ABI_FREE(gkq_nu)

         ! Copy data to improve memory access in the loops below.
         ! (my_npert, nb, my_nq, nb, my_nk)
         !g2_pmnk = gqk%my_g2(my_ip,:,my_iq,:,my_ik)

         ! Compute T_pp'(q,nu) matrix in reduced coordinates.
         if (q_is_gamma) call sigtk_dw_tpp_red(natom, displ_nu_red, tpp_red)

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
             ediff = eig0nk - eig0mk
             !intra_band = q_is_gamma .and. ediff <= TOL_EDIFF
             !same_band = ibsum_kq == band_ks

             if (dtset%eph_ahc_type == 1) then
               cfact(:) =  (nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + ieta) + &
                           (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + ieta)
             else
               cfact(:) =  (two * nqnu + one) / (eig0nk - eig0mkq + ieta)
             end if

             ! Re + Im of Fan-Migdal self-energy
             ! (my_npert, nb_kq, my_nq, nb_k, my_nk)
             gkq2 = weight_q * gqk%my_g2(my_ip, im_kq, my_iq, in_k, my_ik)
             vals_e0ks(:, in_k, ikcalc) = vals_e0ks(:, in_k, ikcalc) + cfact * gkq2
             !print *, "gkq2", gkq2

             ! Derivative of sigma
             ! Accumulate d(Re Sigma) / dw(w=eKS) for state ib_k
             !cfact(x) =  (nqnu + f_mkq      ) / (x - eig0mkq + wqnu + sigma%ieta) + &
             !            (nqnu - f_mkq + one) / (x - eig0mkq - wqnu + sigma%ieta)
             gmod2 = (eig0nk - eig0mkq + wqnu) ** 2
             hmod2 = (eig0nk - eig0mkq - wqnu) ** 2
             rfact(:) = (nqnu + f_mkq      ) * (-gmod2 + aimag(ieta)**2) / (gmod2 + aimag(ieta)**2) ** 2 + &
                        (nqnu - f_mkq + one) * (-hmod2 + aimag(ieta)**2) / (hmod2 + aimag(ieta)**2) ** 2

             dvals_de0ks(:, in_k, ikcalc) = dvals_de0ks(:, in_k, ikcalc) + gkq2 * rfact

             if (q_is_gamma) then
               ! Accumulate DW for each T, add it to Sigma(e0) and Sigma(w) as well
               ! - (2 n_{q\nu} + 1) * gdw2 / (e_nk - e_mk)
               gdw2 = zero
               if (abs(ediff) > EPHTK_WTOL) then
                 cfact(:) = - weight_q * gdw2 * (two * nqnu + one)  / (ediff + ieta)
               else
                 cfact(:) = zero
               endif
               vals_e0ks(:, in_k, ikcalc) = vals_e0ks(:, in_k, ikcalc) + real(cfact)
             end if
           end do ! in_k

         end do ! im_kq
       end do ! my_ip
     end do ! my_iq
     call lg_myk%free()
   end do ! my_ik

   ! TODO: Sternheimer with KS states.

   ! Sum partial terms inside qgk%comm.
   call xmpi_sum(vals_e0ks, gqk%comm%value, ierr)
   call xmpi_sum(dvals_de0ks, gqk%comm%value, ierr)

   ! Average self-energy matrix elements in the degenerate subspace.
   call write_results__(gqk)

   ABI_FREE(vals_e0ks)
   ABI_FREE(dvals_de0ks)
   end associate
 end do ! my_is

 ABI_FREE(kTmesh)
 ABI_FREE(mu_e)
 ABI_FREE(nqnu)
 !ABI_FREE(f_nk)
 ABI_FREE(f_mkq)
 ABI_FREE(cfact)
 ABI_FREE(rfact)
 ABI_FREE(displ_nu_cart)
 ABI_FREE(displ_nu_red)
 ABI_FREE(tpp_red)

 call gstore%free()

contains

subroutine write_results__(gqk)

!Local variables-------------------------------

 type(gqk_t),intent(in) :: gqk
 integer,parameter :: max_ntemp = 50
#if 0
 integer :: ideg,ib,it,ii,iw,nstates,ierr,my_rank,band_ks,ik_ibz,ibc,ib_val,ib_cond,jj
 integer :: nq_ibzk_eff, nelem, imyq, iq_ibz_k, sr_ncid
 logical :: iwrite
 real(dp) :: ravg,kse,kse_prev,dw,fan0,ks_gap,kse_val,kse_cond,qpe_oms,qpe_oms_val,qpe_oms_cond
 real(dp) :: cpu, wall, gflops, invsig2fmts, tau, ravg2
 complex(dpc) :: sig0c,zc,qpe,qpe_prev,qpe_val,qpe_cond,cavg1,cavg2,cavg3,cavg4
 !character(len=5000) :: msg
 integer :: grp_ncid, ncerr
!arrays
 integer, allocatable :: recvcounts(:), displs(:), nq_rank(:), kq_symtab(:,:), my_kq_symtab(:,:)
 integer, ABI_CONTIGUOUS pointer :: bids(:)
 real(dp) :: qp_gaps(self%ntemp),qpoms_gaps(self%ntemp)
 real(dp),allocatable :: aw(:,:,:), a2few_avg(:,:), gather_srate(:,:,:,:), grp_srate(:,:,:,:)
 real(dp) :: ks_enes(self%max_nbcalc), ze0_vals(self%ntemp, self%max_nbcalc)
 real(dp) :: gfw_avg(self%phmesh_size, 3)
 complex(dpc) :: qpoms_enes(self%ntemp, self%max_nbcalc),qp_enes(self%ntemp, self%max_nbcalc)
#endif
! *************************************************************************

 ! Write legend.
 if (spin == 1) then
   write(ab_out,"(a)")repeat("=", 80)
   write(ab_out,"(a)")" Final results in eV."
   write(ab_out,"(a)")" Notations:"
   write(ab_out,"(a)")"     eKS: Kohn-Sham energy. eQP: quasi-particle energy."
   write(ab_out,"(a)")"     eQP - eKS: Difference between the QP and the KS energy."
   write(ab_out,"(a)")"     SE1(eKS): Real part of the self-energy computed at the KS energy, SE2 for imaginary part."
   write(ab_out,"(a)")"     Z(eKS): Renormalization factor."
   write(ab_out,"(a)")"     FAN: Real part of the Fan term at eKS. DW: Debye-Waller term."
   write(ab_out,"(a)")"     DeKS: KS energy difference between this band and band-1, DeQP same meaning but for eQP."
   write(ab_out,"(a)")"     OTMS: On-the-mass-shell approximation with eQP ~= eKS + Sigma(omega=eKS)"
   write(ab_out,"(a)")"     TAU(eKS): Lifetime in femtoseconds computed at the KS energy."
   write(ab_out,"(a)")"     mu_e: Fermi level for given (T, nelect)"
   write(ab_out,"(a)")" "
   write(ab_out,"(a)")" "
 end if

 do ikcalc=1,gqk%glob_nk
   !kcalc =
 do it=1,ntemp
 do in_k=1,gqk%nb ! gqk%nb_k
   print *, "Re SE (eV), Z:", real(vals_e0ks(it, in_k, ikcalc)) * Ha_eV, real(dvals_de0ks(it, in_k, ikcalc))

#if 0
   ! Write header.
   if (it <= max_ntemp) then
     if (nsppol == 1) then
       write(ab_out,"(3a,f6.1,a,f8.3)") &
         "K-point: ", trim(ktoa(kcalc(:,ikcalc))), ", T: ", kTmesh(it) / kb_HaK, &
         " [K], mu_e: ", mu_e(it) * Ha_eV
     else
       write(ab_out,"(3a,i1,a,f6.1,a,f8.3)") &
         "K-point: ", trim(ktoa(kcalc(:,ikcalc))), ", spin: ", spin, ", T: ",kTmesh(it) / kb_HaK, &
         " [K], mu_e: ", mu_e(it) * Ha_eV
     end if
     if (imag_only) then
       write(ab_out,"(a)")"   B    eKS    SE2(eKS)  TAU(eKS)  DeKS"
     else
       write(ab_out,"(a)")"   B    eKS     eQP    eQP-eKS   SE1(eKS)  SE2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP"
     end if
   end if

   do ibc=1,self%nbcalc_ks(ikcalc, spin)
     band_ks = self%bstart_ks(ikcalc, spin) + ibc - 1
     kse = ebands%eig(band_ks, ik_ibz, spin)
     ks_enes(ibc) = kse
     sig0c = self%vals_e0ks(it, ibc)
     dw = self%dw_vals(it, ibc)
     fan0 = real(sig0c) - dw
     ! Compute QP energies with On-the-Mass-Shell approximation and first renormalization i.e. Z(eKS)
     ! TODO: Note that here I use the full Sigma including the imaginary part
     !zc = one / (one - self%dvals_de0ks(it, ibc))
     zc = one / (one - real(self%dvals_de0ks(it, ibc)))
     ze0_vals(it, ibc) = real(zc)
     qpe = kse + real(zc) * real(sig0c)
     qpe_oms = kse + real(sig0c)
     if (ibc == 1) then
       kse_prev = kse; qpe_prev = qpe
     end if
     if (band_ks == ib_val) then
       kse_val = kse; qpe_val = qpe; qpe_oms_val = qpe_oms
     end if
     if (band_ks == ib_cond) then
       kse_cond = kse; qpe_cond = qpe; qpe_oms_cond = qpe_oms
     end if

     if (it <= max_ntemp) then
       if (self%imag_only) then
         ! 1/tau  = 2 Imag(Sigma)
         invsig2fmts = Time_Sec * 1e+15 / two
         tau = 999999.0_dp
         if (abs(aimag(sig0c)) > tol16) tau = invsig2fmts / abs(aimag(sig0c))
         tau = min(tau, 999999.0_dp)
         write(ab_out, "(i4,2(f8.3,1x),f8.1,1x,f8.3)") &
             band_ks, kse * Ha_eV, aimag(sig0c) * Ha_eV, tau, (kse - kse_prev) * Ha_eV
       else
         write(ab_out, "(i4, 10(f8.3,1x))") &
           band_ks, kse * Ha_eV, real(qpe) * Ha_eV, (real(qpe) - kse) * Ha_eV, &
           real(sig0c) * Ha_eV, aimag(sig0c) * Ha_eV, real(zc), &
           fan0 * Ha_eV, dw * Ha_eV, (kse - kse_prev) * Ha_eV, real(qpe - qpe_prev) * Ha_eV
       end if
     end if

     if (ibc > 1) then
       kse_prev = kse; qpe_prev = qpe
     end if
     qpoms_enes(it, ibc) = qpe_oms
     qp_enes(it, ibc) = qpe
     if (kse_val /= huge(one) * tol6 .and. kse_cond /= huge(one) * tol6) then
       ! We have enough states to compute the gap.
       if (it == 1) ks_gap = kse_cond - kse_val
       qpoms_gaps(it) = qpe_oms_cond - qpe_oms_val
       qp_gaps(it) = real(qpe_cond - qpe_val)
     end if
   end do ! ibc

   ! Print KS and QP gaps.
   if (it <= max_ntemp) then
     if (.not. self%imag_only) then
       if (kse_val /= huge(one) * tol6 .and. kse_cond /= huge(one) * tol6) then
         write(ab_out, "(a)")" "
         write(ab_out, "(a,f8.3,1x,2(a,i0),a)")" KS gap: ",ks_gap * Ha_eV, &
           "(assuming bval:", ib_val, " ==> bcond:", ib_cond, ")"
         write(ab_out, "(2(a,f8.3),a)")" QP gap: ",qp_gaps(it) * Ha_eV," (OTMS: ",qpoms_gaps(it) * Ha_eV, ")"
         write(ab_out, "(2(a,f8.3),a)")" QP_gap - KS_gap: ",(qp_gaps(it) - ks_gap) * Ha_eV,&
             " (OTMS: ",(qpoms_gaps(it) - ks_gap) * Ha_eV, ")"
         write(ab_out, "(a)")" "
       end if
     else
       if (kse_val /= huge(one) * tol6 .and. kse_cond /= huge(one) * tol6) then
         write(ab_out, "(a)")" "
         write(ab_out, "(a,f8.3,1x,2(a,i0),a)")" KS gap: ",ks_gap * Ha_eV, "(assuming bval:",ib_val," ==> bcond:",ib_cond,")"
         write(ab_out, "(a)")" "
       end if
     end if

     write(ab_out, "(a)")repeat("=", 92)
   end if

 end do ! it

 if (self%ntemp > max_ntemp .and. (ikcalc == 1 .and. spin == 1)) then
   write(ab_out, "(a,i0,a)")" No more than ", max_ntemp, " temperatures are written to the main output file."
   write(ab_out, "(2a)")" Please use SIGEPH.nc file and AbiPy to analyze the results.",ch10
 end if
#endif

 end do ! in_k
 end do ! it
 end do ! ikcalc
 stop

end subroutine write_results__

end subroutine gstore_sigeph
!!***

end module m_gstore_sigeph
!!***
