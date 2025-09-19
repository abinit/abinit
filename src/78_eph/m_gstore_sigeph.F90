!!****m* ABINIT/m_gstore
!! NAME
!! m_gstore
!!
!! FUNCTION
!!  Compute matrix elements of the e-ph self-energy (Fan Migdal + Debye Waller).
!!  using precomputed e-ph matrix elements.
!!  See also m_sigmaph, for a version in which the g-matrix elements are computed on-the-fly.
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
 integer :: spin, my_is, my_ik, my_iq, my_ip, in_k, im_kq, ierr, ntemp, gap_err, my_rank, ip1, ip2
 integer :: it, ik_bz, ik_ibz, ikq_ibz, band_k, band_kq, timrev, ii, ikcalc, natom, natom3, nsppol
 real(dp) :: wqnu, gkq2, weight_q, eig0nk, eig0mk, eig0mkq, ediff, gmod2, hmod2, gdw2, gdw2_stern, rtmp !,nqnu,gkq2,gkq2_pf,
 !real(dp) :: cpu, wall, gflops
 logical :: use_lgk, q_is_gamma, intra_band, same_band, imag_only
 complex(dp) :: ieta, cfact !, sig_cplx
 type(gaps_t) :: gaps
 type(lgroup_t) :: lg_myk
 type(gstore_t) :: gstore
!arrays
 integer :: units(2), my_kqmap(6)
 integer,allocatable :: phmodes_skip(:)
 real(dp) :: kk(3), qpt(3), kq(3)
 real(dp),allocatable :: gkq_atm(:,:,:),gkq_nu(:,:,:) !,gkq0_atm(:,:,:,:), gaussw_qnu(:)
 real(dp),allocatable :: rfact_t(:), mu_e(:), kTmesh(:), nqnu(:), f_mkq(:) !, f_nk(:),  g2_pmnk(:,:,:,:)
 real(dp),allocatable :: displ_nu_cart(:,:,:), displ_nu_red(:,:,:), dw_vals(:,:,:)
 complex(dp),allocatable :: cfact_t(:), vals_e0ks(:,:,:), dvals_de0ks(:,:,:), tpp_red(:,:) !, fmw_frohl_sphcorr(:,:,:,:), cfact_wr(:),
!----------------------------------------------------------------------

 units = [std_out, ab_out]
 my_rank = xmpi_comm_rank(comm)

 call wrtout(std_out, " Computing Fan-Migdal + DW self-energy from GSTORE.nc", pre_newlines=1)
 natom = cryst%natom; natom3 = 3 * cryst%natom; nsppol = gstore%nsppol
 imag_only = .False.

 ! The Fan-Migdal SE requires |g(k,q)| in the phonon representation but
 ! to compute the DW term in the RIA, we need complex g in the atom representation.
 call gstore%from_ncpath(dtfil%filgstorein, with_cplex1, dtset, cryst, ebands, ifc, comm, &
                         with_gmode="phonon", gvals_name=dtset%gstore_gname, read_dw=.True.)

 ABI_CHECK(gstore%qzone == "bz", "gstore_sigeph assumes qzone == `bz`")

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
 ABI_MALLOC(cfact_t, (ntemp))
 ABI_MALLOC(rfact_t, (ntemp))

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

 ! Setup a mask to skip accumulating the contribution of certain phonon modes.
 !call ephtk_set_phmodes_skip(dtset%natom, dtset%eph_phrange, phmodes_skip)
 !ABI_FREE(phmodes_skip)

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
   ABI_CALLOC(dw_vals, (ntemp, gqk%nb, gqk%glob_nk))  ! nb_k

   ! Loop over k-points in |n,k>
   do my_ik=1,gqk%my_nk
     kk = gqk%my_kpts(:, my_ik); ik_ibz = gqk%my_k2ibz(1, my_ik)
     ! Will store results using glob_ik index.
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
         ! Ignore unstable modes or modes that should be skipped.
         if (wqnu < EPHTK_WTOL) cycle

         nqnu(:) = occ_be(wqnu, kTmesh, zero)
         !print *, "wqnu", wqnu

         ! FIXME: Not sure this is correct!
         displ_nu_cart = gqk%my_displ_cart(:,:,:,my_ip,my_iq)
         call phdispl_cart2red_nmodes(natom, 1, cryst%gprimd, displ_nu_cart, displ_nu_red)

         ! For the DW, we need gkq_atm at q = 0, for all atomic perturbations
         ! Copy data to improve memory access in the loops below.
         ! (my_npert, nb, my_nq, nb, my_nk)
         !g2_pmnk = gqk%my_g2(my_ip,:,my_iq,:,my_ik)

         ! Compute T_pp'(q,nu) matrix in reduced coordinates.
         call sigtk_dw_tpp_red(natom, displ_nu_red, tpp_red)

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
             !same_band = ibsum_kq == band_k

             if (dtset%eph_ahc_type == 1) then
               cfact_t(:) =  (nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + ieta) + &
                             (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + ieta)
             else
               cfact_t(:) =  (two * nqnu + one) / (eig0nk - eig0mkq + ieta)
             end if

             ! Re + Im of Fan-Migdal self-energy
             ! (my_npert, nb_kq, my_nq, nb_k, my_nk)
             gkq2 = weight_q * gqk%my_g2(my_ip, im_kq, my_iq, in_k, my_ik)
             vals_e0ks(:, in_k, ikcalc) = vals_e0ks(:, in_k, ikcalc) + cfact_t * gkq2
             !print *, "gkq2", gkq2

             ! Derivative of sigma
             ! Accumulate d(Re Sigma) / dw(w=eKS) for state in_k
             !cfact(x) =  (nqnu + f_mkq      ) / (x - eig0mkq + wqnu + sigma%ieta) + &
             !            (nqnu - f_mkq + one) / (x - eig0mkq - wqnu + sigma%ieta)
             gmod2 = (eig0nk - eig0mkq + wqnu) ** 2
             hmod2 = (eig0nk - eig0mkq - wqnu) ** 2
             rfact_t(:) = (nqnu + f_mkq      ) * (-gmod2 + aimag(ieta)**2) / (gmod2 + aimag(ieta)**2) ** 2 + &
                        (nqnu - f_mkq + one) * (-hmod2 + aimag(ieta)**2) / (hmod2 + aimag(ieta)**2) ** 2

             dvals_de0ks(:, in_k, ikcalc) = dvals_de0ks(:, in_k, ikcalc) + gkq2 * rfact_t

             ! Compute DW term following XG paper. Check prefactor.
             ! gkq0_atm(2, nbcalc_ks, bsum_start:bsum_stop, natom3)
             ! (nb_k, nb_kq, natom3, my_nk)
             associate (gkq0_atm => gqk%my_gq0nm_atm(:,:,:,my_ik))
             gdw2 = zero
             do ip2=1,natom3
               do ip1=1,natom3
                 cfact = ( &
                   + real(gkq0_atm(in_k, im_kq, ip1)) * real(gkq0_atm(in_k, im_kq, ip2)) &
                   + aimag(gkq0_atm(in_k, im_kq, ip1)) * aimag(gkq0_atm(in_k, im_kq, ip2)) &
                   + real(gkq0_atm(in_k, im_kq, ip2)) * real(gkq0_atm(in_k, im_kq, ip1)) &
                   + aimag(gkq0_atm(in_k, im_kq, ip2)) * aimag(gkq0_atm(in_k, im_kq, ip1)) &
                  !+ gkq0_atm(1, in_k, im_kq, ip1) * gkq0_atm(1, in_k, im_kq, ip2) &
                  !+ gkq0_atm(2, in_k, im_kq, ip1) * gkq0_atm(2, in_k, im_kq, ip2) &
                  !+ gkq0_atm(1, in_k, im_kq, ip2) * gkq0_atm(1, in_k, im_kq, ip1) &
                  !+ gkq0_atm(2, in_k, im_kq, ip2) * gkq0_atm(2, in_k, im_kq, ip1) &
                 )
                 !
                 gdw2 = gdw2 + real(tpp_red(ip1,ip2) * cfact)
               end do
             end do
             end associate
             gdw2 = gdw2 / (four * two * wqnu)
             !print *, "gdw2", gdw2

             ! Accumulate DW for each T, add it to Sigma(e0) and Sigma(w) as well
             ! - (2 n_{q\nu} + 1) * gdw2 / (e_nk - e_mk)
             if (abs(ediff) > EPHTK_WTOL) then
               cfact_t(:) = - weight_q * gdw2 * (two * nqnu + one)  / (ediff + ieta)
             else
               cfact_t(:) = zero
             endif
             dw_vals(:, in_k, ikcalc) = dw_vals(:, in_k, ikcalc) + real(cfact_t)
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
   call xmpi_sum(dw_vals, gqk%comm%value, ierr)

   ! TODO: Average self-energy matrix elements in the degenerate subspace.
   call write_results__(ntemp, gqk)

   ABI_FREE(vals_e0ks)
   ABI_FREE(dvals_de0ks)
   ABI_FREE(dw_vals)
   end associate
 end do ! my_is

 ABI_FREE(kTmesh)
 ABI_FREE(mu_e)
 ABI_FREE(nqnu)
 !ABI_FREE(f_nk)
 ABI_FREE(f_mkq)
 ABI_FREE(cfact_t)
 ABI_FREE(rfact_t)
 ABI_FREE(displ_nu_cart)
 ABI_FREE(displ_nu_red)
 ABI_FREE(tpp_red)

 call gstore%free()

contains

subroutine write_results__(ntemp, gqk)

!Local variables-------------------------------
 integer,intent(in) :: ntemp
 type(gqk_t),intent(in) :: gqk
 integer,parameter :: max_ntemp = 50
 integer :: band_k,ik_ibz,ib_val,ib_cond,jj !,ideg,ib,it,ii,iw,nstates
 !integer :: nq_ibzk_eff, nelem, imyq, iq_ibz_k, sr_ncid
 !logical :: iwrite
 real(dp) :: ravg,kse,kse_prev,dw,fan0,ks_gap,kse_val,kse_cond,qpe_oms,qpe_oms_val,qpe_oms_cond
 !real(dp) :: invsig2fmts, tau, ravg2
 complex(dp) :: sig0c,zc,qpe,qpe_prev,qpe_val,qpe_cond !,cavg1,cavg2,cavg3,cavg4
 !character(len=5000) :: msg
 !integer :: grp_ncid, ncerr
!arrays
 real(dp) :: kcalc(3)
 !integer, allocatable :: recvcounts(:), displs(:), nq_rank(:), kq_symtab(:,:), my_kq_symtab(:,:)
 !integer, ABI_CONTIGUOUS pointer :: bids(:)
 real(dp) :: qp_gaps(ntemp),qpoms_gaps(ntemp)
 !real(dp),allocatable :: aw(:,:,:), a2few_avg(:,:), gather_srate(:,:,:,:), grp_srate(:,:,:,:)
 real(dp) :: ks_enes(gqk%nb), ze0_vals(ntemp, gqk%nb)
 !real(dp) :: gfw_avg(self%phmesh_size, 3)
 complex(dp) :: qpoms_enes(ntemp, gqk%nb),qp_enes(ntemp, gqk%nb) ! nb_k
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

 ! Compute QP energies and Gaps (Note that I'm assuming a non-magnetic semiconductor!)
 ib_val = nint(ebands%nelect / (two / ebands%nspinor)); ib_cond = ib_val + 1

 do ikcalc=1,gqk%glob_nk
   ik_bz = gstore%kglob2bz(ikcalc, spin)
   ik_ibz = gstore%kbz2ibz(1, ik_bz)
   kcalc = gstore%kbz(:, ik_bz)

   kse_val = huge(one) * tol6; kse_cond = huge(one) * tol6
   qp_enes = huge(one) * tol6; qpoms_enes = huge(one) * tol6
   ks_enes = huge(one) * tol6; ze0_vals = huge(one) * tol6
   ks_gap = -one; qpoms_gaps = -one; qp_gaps = -one

   do it=1,ntemp

     ! Write header.
     if (it <= max_ntemp) then
       if (nsppol == 1) then
         write(ab_out,"(3a,f6.1,a,f8.3)") &
           "K-point: ", trim(ktoa(kcalc)), ", T: ", kTmesh(it) / kb_HaK, &
           " [K], mu_e: ", mu_e(it) * Ha_eV
       else
         write(ab_out,"(3a,i1,a,f6.1,a,f8.3)") &
           "K-point: ", trim(ktoa(kcalc)), ", spin: ", spin, ", T: ",kTmesh(it) / kb_HaK, &
           " [K], mu_e: ", mu_e(it) * Ha_eV
       end if
       if (imag_only) then
         write(ab_out,"(a)")"   B    eKS    SE2(eKS)  TAU(eKS)  DeKS"
       else
         write(ab_out,"(a)")"   B    eKS     eQP    eQP-eKS   SE1(eKS)  SE2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP"
       end if
     end if

     ! Loop over bands for this k-point and spin
     do in_k=1,gqk%nb ! gqk%nb_k
       !print *, "Re SE (eV), Z:", real(vals_e0ks(it, in_k, ikcalc)) * Ha_eV, real(dvals_de0ks(it, in_k, ikcalc))
       band_k = in_k - gqk%bstart + 1
       kse = ebands%eig(band_k, ik_ibz, spin)
       ks_enes(in_k) = kse
       sig0c = vals_e0ks(it, in_k, ikcalc)
       dw = dw_vals(it, in_k, ikcalc)
       fan0 = real(sig0c) - dw
       ! Compute QP energies with On-the-Mass-Shell approximation and first renormalization i.e. Z(eKS)
       ! TODO: Note that here I use the full Sigma including the imaginary part
       !zc = one / (one - dvals_de0ks(it, in_k))
       zc = one / (one - real(dvals_de0ks(it, in_k, ikcalc)))
       ze0_vals(it, in_k) = real(zc)
       qpe = kse + real(zc) * real(sig0c)
       qpe_oms = kse + real(sig0c)
       if (in_k == 1) then
         kse_prev = kse; qpe_prev = qpe
       end if
       if (band_k == ib_val) then
         kse_val = kse; qpe_val = qpe; qpe_oms_val = qpe_oms
       end if
       if (band_k == ib_cond) then
         kse_cond = kse; qpe_cond = qpe; qpe_oms_cond = qpe_oms
       end if

       if (it <= max_ntemp) then
         if (imag_only) then
           ! 1/tau  = 2 Imag(Sigma)
           !invsig2fmts = Time_Sec * 1e+15 / two
           !tau = 999999.0_dp
           !if (abs(aimag(sig0c)) > tol16) tau = invsig2fmts / abs(aimag(sig0c))
           !tau = min(tau, 999999.0_dp)
           !write(ab_out, "(i4,2(f8.3,1x),f8.1,1x,f8.3)") &
           !    band_k, kse * Ha_eV, aimag(sig0c) * Ha_eV, tau, (kse - kse_prev) * Ha_eV
         else
           write(ab_out, "(i4, 10(f8.3,1x))") &
             band_k, kse * Ha_eV, real(qpe) * Ha_eV, (real(qpe) - kse) * Ha_eV, &
             real(sig0c) * Ha_eV, aimag(sig0c) * Ha_eV, real(zc), &
             fan0 * Ha_eV, dw * Ha_eV, (kse - kse_prev) * Ha_eV, real(qpe - qpe_prev) * Ha_eV
         end if
       end if

       if (in_k > 1) then
         kse_prev = kse; qpe_prev = qpe
       end if
       qpoms_enes(it, in_k) = qpe_oms
       qp_enes(it, in_k) = qpe
       if (kse_val /= huge(one) * tol6 .and. kse_cond /= huge(one) * tol6) then
         ! We have enough states to compute the gap.
         if (it == 1) ks_gap = kse_cond - kse_val
         qpoms_gaps(it) = qpe_oms_cond - qpe_oms_val
         qp_gaps(it) = real(qpe_cond - qpe_val)
       end if
     end do ! in_k

     ! Print KS and QP gaps.
     if (it <= max_ntemp) then
       if (.not. imag_only) then
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
 end do ! ikcalc

 if (ntemp > max_ntemp) then
   write(ab_out, "(a,i0,a)")" No more than ", max_ntemp, " temperatures are written to the main output file."
   write(ab_out, "(2a)")" Please use SIGEPH.nc file and AbiPy to analyze the results.",ch10
 end if

 call wrtout(std_out, "gstore_sigeph ended OK")
 stop

end subroutine write_results__

end subroutine gstore_sigeph
!!***

end module m_gstore_sigeph
!!***
