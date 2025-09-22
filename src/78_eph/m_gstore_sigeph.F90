!!****m* ABINIT/m_gstore_sigeph
!! NAME
!! m_gstore_sigeph
!!
!! FUNCTION
!!  Compute (diagonal) matrix elements of the e-ph self-energy (Fan Migdal + Debye Waller).
!!  in the KS basis using precomputed e-ph matrix elements.
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
 use netcdf
 use m_nctk
 !use m_dvdb
 use m_crystal,        only : crystal_t
 !use m_hamiltonian
 use m_dtset,          only : dataset_type
 use m_dtfil,          only : datafiles_type
 !use m_wfd
 use m_ephtk
 !use m_mkffnl
 use m_sigtk

 !use m_time,           only : cwtime, cwtime_report, sec2str
 use m_fstrings,       only : tolower, itoa, ftoa, sjoin, ktoa, ltoa, strcat, replace_ch0, yesno, string_in
 !use m_cgtools,        only : cg_zdotc
 !use m_kg,             only : getph
 !use defs_datatypes,   only : pseudopotential_type
 use defs_abitypes,    only : mpi_type
 use m_hdr,            only : hdr_type, fform_from_ext
 use m_geometry,       only : phdispl_cart2red_nmodes
 use m_ebands,         only : ebands_t, gaps_t
 use m_kpts,           only : kpts_timrev_from_kptopt, kpts_map !, kpts_sort, kpts_pack_in_stars, kptrlatt_from_ngkpt
 use m_ioarr,          only : read_rhor
 !use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack
 use m_ifc,            only : ifc_type
 use m_dfpt_cgwf,      only : stern_t
 !use m_pawang,         only : pawang_type
 !use m_pawrad,         only : pawrad_type
 !use m_pawtab,         only : pawtab_type
 !use m_pawfgr,         only : pawfgr_type
 use m_pawrhoij,       only : pawrhoij_type
 !use m_pstat,          only : pstat_proc
 use m_occ,            only : occ_be, occ_fd
 use m_lgroup,         only : lgroup_t
 use m_gstore,         only : gstore_t, gqk_t

 implicit none

 private
 public :: gstore_sigeph

 real(dp),private,parameter :: TOL_EDIFF = 0.001_dp * eV_Ha
!!***

!----------------------------------------------------------------------

!!****t* m_gstore_sigeph/sep_t
!! NAME
!! sep_t
!!
!! FUNCTION
!! Container for the (diagonal) matrix elements of the electron-phonon self-energy
!! in the KS representation i.e. Sigma_eph(omega, T, band, k, spin).
!! Provides methods to compute QP corrections, spectral functions, QP linewidths and
!! save the results to netcdf file.
!!
!! TODO
!!  Fix problem with spin parallelism.
!!
!! SOURCE

 type,public :: sep_t

  integer :: nwr = 0
   ! Number of frequency points along the real axis for Sigma(w) and spectral function A(w)
   ! Odd number so that the mesh is centered on the KS energy.
   ! The spectral function is computed only if nwr > 0 (taken from dtset%nfreqsp)

  integer :: ntemp = 0
  ! Number of temperatures.

  logical :: imag_only

  real(dp) :: wr_step
   ! Step of the linear mesh along the real axis (Ha units).

  complex(dp) :: ieta = zero
   ! Used to shift the poles in the complex plane (Ha units)
   ! Corresponds to `i eta` term in equations.

  real(dp),allocatable :: kTmesh(:)
  ! kTmesh(ntemp)
  ! List of temperatures (kT units).

  real(dp),allocatable :: mu_e(:)
  ! mu_e(ntemp)
  ! chemical potential of electrons for the different temperatures.

  complex(dp),allocatable :: vals_e0ks(:,:,:)
  ! Sigma_eph(omega=eKS, kT, band) for given (ikcalc, spin).
  ! Fan-Migdal + Debye-Waller

  complex(dp),allocatable :: fan_vals(:,:,:)
  ! (ntemp, nb_k, nkcalc)
  ! Fan-Migdal

  !complex(dp),allocatable :: fan_stern_vals(:,:,:)
  ! (ntemp, nb_k, nkcalc)
  ! Fan-Migdal adiabatic Sternheimer part

  complex(dp),allocatable :: dvals_de0ks(:,:,:)
  ! (ntemp, nb_k, nkcalc)
  ! d Re Sigma_eph(omega, kT, band, kcalc) / d omega (omega=eKS)

  real(dp),allocatable :: dw_vals(:,:,:)
  !  dw_vals(ntemp, nb_k, nkcalc) for given (ikcalc, spin)
  !  Debye-Waller term (static).

  !real(dp),allocatable :: dw_stern_vals(:,:,:)
   !  dw_stern_vals(ntemp, nb_k, nkcalc)
   !  Debye-Waller Sternheimer term (static) .

  !complex(dp),allocatable :: vals_wr(:,:,:,:)
   ! vals_wr(nwr, ntemp, nb_k, nkcalc)
   ! Sigma_eph(omega, kT, band)
   ! enk_KS corresponds to nwr/2 + 1.

  !integer :: phmesh_size
   ! Number of phonon frequencies in phonon mesh used for Eliashberg functions and
   ! and other omega-resolved quantities.

  !real(dp),allocatable :: phmesh(:)
   ! phmesh(phmesh_size)
   ! phonon mesh in Ha.

 contains

 procedure :: gather_and_write_results => sep_gather_and_write_results
 ! Write main dimensions and header of sigmaph on a netcdf file.

 procedure :: free => sep_free
 ! Free dynamic memory

 end type sep_t
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_gstore_sigeph/gstore_sigeph
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

subroutine gstore_sigeph(ngfft, ngfftf, dtset, dtfil, cryst, ebands, ifc, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(ifc_type),target,intent(in) :: ifc
 type(mpi_type),intent(inout) :: mpi_enreg
 integer,intent(in) :: comm
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)

!Local variables-------------------------------
 integer,parameter :: master = 0, with_cplex1 = 1, max_ntemp = 50, cplex1 = 1, pawread0 = 0
 integer :: n1, n2, n3, n4, n5, n6, nb_k, nb_kq, bstart_k, bstart_kq
 integer :: spin, my_is, my_ik, my_iq, my_ip, in_k, im_kq, ierr, gap_err, my_rank, ip1, ip2, nu
 integer :: it, ik_bz, ik_ibz, ikq_ibz, band_k, band_kq, timrev_k, ii, ikcalc, natom, natom3, nsppol
 integer :: nfft, nfftf, mgfft, mgfftf !,nkpg,nkpg1,nq,cnt,imyp, q_start, q_stop, restart, enough_stern
 real(dp) :: wqnu, gkq2, weight_q, eig0nk, eig0mk, eig0mkq, ediff, gmod2, hmod2, gdw2, gdw2_stern, rtmp !,nqnu,gkq2,gkq2_pf,
 !real(dp) :: cpu, wall, gflops
 logical :: q_is_gamma, intra_band, same_band
 complex(dp) :: cfact !, sig_cplx
 character(len=5000) :: msg
 type(gaps_t) :: gaps
 type(lgroup_t) :: lg_myk
 type(gstore_t) :: gstore
 type(sep_t) :: sigma
 type(stern_t) :: stern
 type(hdr_type) :: pot_hdr
 type(crystal_t) :: pot_cryst
!arrays
 integer :: units(2), my_kqmap(6)
 integer,allocatable :: phmodes_skip(:)
 real(dp) :: kk(3), qpt(3), kq(3)
 real(dp),allocatable :: vtrial(:,:) !,gvnlx1(:,:),work(:,:,:,:), vcar_ibz(:,:,:,:)
 real(dp),allocatable :: gkq_atm(:,:,:),gkq_nu(:,:,:) !,gkq0_atm(:,:,:,:), gaussw_qnu(:)
 real(dp),allocatable :: rfact_t(:), nqnu(:), f_mkq(:) !, f_nk(:),  g2_pmnk(:,:,:,:)
 real(dp),allocatable :: displ_nu_cart(:,:,:), displ_nu_red(:,:,:)
 complex(dp),allocatable :: cfact_t(:), tpp_red(:,:) !, fmw_frohl_sphcorr(:,:,:,:), cfact_wr(:),
 type(pawrhoij_type),allocatable :: pot_pawrhoij(:)
!----------------------------------------------------------------------

 units = [std_out, ab_out]
 my_rank = xmpi_comm_rank(comm)

 call wrtout(std_out, " Computing Fan-Migdal + DW self-energy from GSTORE.nc", pre_newlines=1)
 natom = cryst%natom; natom3 = 3 * cryst%natom; nsppol = gstore%nsppol

 sigma%imag_only = .False.
 sigma%ieta = + j_dpc * dtset%zcut

 ! The Fan-Migdal SE requires |g(k,q)| in the phonon representation but
 ! to compute the DW term in the RIA, we need complex g in the atom representation.
 call gstore%from_ncpath(dtfil%filgstorein, with_cplex1, dtset, cryst, ebands, ifc, comm, &
                         with_gmode="phonon", gvals_name=dtset%gstore_gname, read_dw=.True.)

 ! Consistency check.
 ierr = 0
 if (gstore%qzone /= "bz") then
   ABI_ERROR_NOSTOP("gstore_sigeph assumes qzone == `bz`", ierr)
 end if
 if (gstore%has_used_lgk /= 0) then
   ABI_ERROR_NOSTOP("gstore_sigeph formalism does not support use_lgk /=0 .", ierr)
 end if
 if (gstore%has_used_lgq /= 0) then
   ABI_ERROR_NOSTOP("gstore_sigeph does not support use_lgq /=0 .", ierr)
 end if
 if (ierr > 1) then
   write(msg,'(a,i0,5a)')&
     'Checking consistency of input data against itself gave ',ierr,' inconsistencies.',ch10,&
     'The details of the problems can be FOUND ABOVE (or in output or log file), in an earlier WARNING.',ch10,&
     'In parallel, the details might not even be printed there. Then, try running in sequential to see the details.'
   ABI_ERROR(msg)
 end if

 ! Check consistency of little group options
 ABI_CHECK(gstore%check_little_group(dtset, msg) == 0, msg)

 ! FFT meshes from input file, not necessarily equal to the ones found in the external files.
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Build (linear) mesh of K * temperatures. tsmesh(1:3) = [start, step, num]
 call dtset%get_ktmesh(sigma%ntemp, sigma%kTmesh)

 ! Compute the chemical potential at the different physical temperatures with Fermi-Dirac.
 ! TODO: One should check that nband is > nbocc to avoid inaccuracies in mu_e.
 ABI_MALLOC(sigma%mu_e, (sigma%ntemp))
 sigma%mu_e(:) = ebands%fermie
 if (dtset%eph_fermie == zero) then
   call ebands%get_muT_with_fd(sigma%ntemp, sigma%ktmesh, dtset%spinmagntarget, dtset%prtvol, sigma%mu_e, gstore%comm)
 end if

 gaps = ebands%get_gaps(gap_err)
 if (gap_err /= 0) then
   ABI_ERROR("Cannot compute fundamental and direct gap (likely metal)")
 end if

 if (my_rank == master) then
   call gaps%print(units, kTmesh=sigma%ktmesh, mu_e=sigma%mu_e, &
                   header="Gaps, band edges and relative position wrt Fermi level")
 end if
 call gaps%free()

 ! Frequency mesh for sigma(w) and spectral functions.
 ! Use GW variables but change default values
 sigma%nwr = dtset%nfreqsp; sigma%wr_step = zero
 if (sigma%nwr > 0) then
   if (mod(sigma%nwr, 2) == 0) sigma%nwr = sigma%nwr + 1
   sigma%wr_step = two * eV_Ha / (sigma%nwr - 1)
   if (dtset%freqspmax /= zero) sigma%wr_step = dtset%freqspmax / (sigma%nwr - 1)
 end if

 ABI_MALLOC(nqnu, (sigma%ntemp))
 !ABI_MALLOC(f_nk, (sigma%ntemp))
 ABI_MALLOC(f_mkq, (sigma%ntemp))
 ABI_MALLOC(cfact_t, (sigma%ntemp))
 ABI_MALLOC(rfact_t, (sigma%ntemp))

 if (dtset%eph_stern /= 0) then
   ! Read the GS potential (vtrial) from input POT file
   ! In principle one may store vtrial in the DVDB but getpot_filepath is simpler to implement.
   call wrtout(units, sjoin(" Reading GS KS potential for Sternheimer from: ", dtfil%filpotin))
   call read_rhor(dtfil%filpotin, cplex1, dtset%nspden, nfftf, ngfftf, pawread0, mpi_enreg, vtrial, pot_hdr, pot_pawrhoij, comm, &
                  allow_interp=.True., want_varname="vtrial")
   pot_cryst = pot_hdr%get_crystal()
   if (gstore%cryst%compare(pot_cryst, header=" Comparing input crystal with POT crystal") /= 0) then
     ABI_ERROR("Crystal structure from WFK and POT do not agree! Check messages above!")
   end if
   call pot_cryst%free(); call pot_hdr%free()
   NOT_IMPLEMENTED_ERROR()
 end if

 ! Compute electronic occ for all Temps (note sigma%mu_e(it) Fermi level)
 !do it=1,sigma%ntemp
 !  f_nk(it) = occ_fd(eig0nk, sigma%kTmesh(it), sigma%mu_e(it))
 !end do ! it

 ABI_MALLOC(displ_nu_cart, (2, 3, natom))
 ABI_MALLOC(displ_nu_red, (2, 3, natom))
 ABI_MALLOC(tpp_red, (natom3, natom3))

 ! Setup a mask to skip accumulating the contribution of certain phonon modes.
 call ephtk_set_phmodes_skip(dtset%natom, dtset%eph_phrange, phmodes_skip)

 ! Loop over collinear spins.
 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is), cryst => gstore%cryst)
   spin = gstore%my_spins(my_is)
   nb_k = gqk%nb_k; nb_kq = gqk%nb_kq
   bstart_k = gqk%bstart_k; bstart_kq = gqk%bstart_kq

   ABI_CHECK(allocated(gqk%my_g2), "my_g2 is not allocated")
   ABI_CHECK(allocated(gqk%my_wnuq), "my_wnuq is not allocated")

   ABI_CALLOC(sigma%vals_e0ks, (sigma%ntemp, nb_k, gqk%glob_nk))
   ABI_CALLOC(sigma%dvals_de0ks, (sigma%ntemp, nb_k, gqk%glob_nk))
   ABI_CALLOC(sigma%fan_vals, (sigma%ntemp, nb_k, gqk%glob_nk))
   ABI_CALLOC(sigma%dw_vals, (sigma%ntemp, nb_k, gqk%glob_nk))

   ! Loop over k-points in |n,k>
   do my_ik=1,gqk%my_nk
     kk = gqk%my_kpts(:, my_ik); ik_ibz = gqk%my_k2ibz(1, my_ik)
     ! Will store results using glob_ik index.
     ikcalc = gqk%my_k2glob(my_ik)

     ! Compute the little group of the k-point so that we can compute g(k,q) only for q in the IBZ_k.
     if (dtset%gstore_use_lgk /= 0) then
       timrev_k = kpts_timrev_from_kptopt(ebands%kptopt)
       call lg_myk%init(cryst, kk, timrev_k, gstore%nqbz, gstore%qbz, gstore%nqibz, gstore%qibz, xmpi_comm_self)
     end if

     ! Sum over q-points.
     do my_iq=1,gqk%my_nq
       call gqk%myqpt(my_iq, gstore, weight_q, qpt)
       q_is_gamma = sum(qpt**2) < tol14

       ! weight_q is computed here. It depends whether we are assuming over the full BZ or IBZ_k.
       weight_q = one / gstore%nqbz
       if (dtset%gstore_use_lgk /= 0) then
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
         nu = my_ip + gqk%my_pert_start - 1

         ! Ignore unstable modes or modes that should be skipped.
         if (ephtk_skip_phmode(nu, wqnu, phmodes_skip, dtset%eph_phrange_w)) cycle

         nqnu(:) = occ_be(wqnu, sigma%kTmesh, zero)
         !print *, "wqnu", wqnu

         ! FIXME: Not sure this is correct!
         displ_nu_cart = gqk%my_displ_cart(:,:,:,my_ip,my_iq)
         call phdispl_cart2red_nmodes(natom, 1, cryst%gprimd, displ_nu_cart, displ_nu_red)

         ! For the DW, we need gkq_atm at q = 0, for all atomic perturbations
         ! Copy data to improve memory access in the loops below.
         ! (my_npert, nb_kq, my_nq, nb_k, my_nk)
         !g2_pmnk = gqk%my_g2(my_ip,:,my_iq,:,my_ik)

         ! Compute T_pp'(q,nu) matrix in reduced coordinates.
         call sigtk_dw_tpp_red(natom, displ_nu_red, tpp_red)

         ! Sum over bands.
         do im_kq=1,nb_kq
           band_kq = im_kq - bstart_kq + 1
           eig0mkq = ebands%eig(band_kq, ikq_ibz, spin)

           ! Compute electronic occ for all Temps (note mu_e(it) Fermi level)
           do it=1,sigma%ntemp
             f_mkq(it) = occ_fd(eig0mkq, sigma%kTmesh(it), sigma%mu_e(it))
           end do

           ! Loop over the n index in |n,k>.
           do in_k=1,nb_k
             band_k = in_k - bstart_k + 1
             eig0nk = ebands%eig(band_k, ik_ibz, spin)
             ediff = eig0nk - eig0mk
             !intra_band = q_is_gamma .and. ediff <= TOL_EDIFF
             !same_band = ibsum_kq == band_k

             if (dtset%eph_ahc_type == 1) then
               cfact_t(:) =  (nqnu + f_mkq      ) / (eig0nk - eig0mkq + wqnu + sigma%ieta) + &
                             (nqnu - f_mkq + one) / (eig0nk - eig0mkq - wqnu + sigma%ieta)
             else
               cfact_t(:) =  (two * nqnu + one) / (eig0nk - eig0mkq + sigma%ieta)
             end if

             ! Re + Im of Fan-Migdal self-energy
             ! (my_npert, nb_kq, my_nq, nb_k, my_nk)
             gkq2 = weight_q * gqk%my_g2(my_ip, im_kq, my_iq, in_k, my_ik)
             !print *, "gkq2", gkq2
             cfact_t = cfact_t * gkq2
             sigma%vals_e0ks(:, in_k, ikcalc) = sigma%vals_e0ks(:, in_k, ikcalc) + cfact_t
             sigma%fan_vals(:, in_k, ikcalc) = sigma%fan_vals(:, in_k, ikcalc) + cfact_t

             ! Derivative of sigma
             ! Accumulate d(Re Sigma) / dw(w=eKS) for state in_k
             !cfact(x) =  (nqnu + f_mkq      ) / (x - eig0mkq + wqnu + sigma%ieta) + &
             !            (nqnu - f_mkq + one) / (x - eig0mkq - wqnu + sigma%ieta)
             gmod2 = (eig0nk - eig0mkq + wqnu) ** 2
             hmod2 = (eig0nk - eig0mkq - wqnu) ** 2
             rfact_t(:) = (nqnu + f_mkq      ) * (-gmod2 + aimag(sigma%ieta)**2) / (gmod2 + aimag(sigma%ieta)**2) ** 2 + &
                          (nqnu - f_mkq + one) * (-hmod2 + aimag(sigma%ieta)**2) / (hmod2 + aimag(sigma%ieta)**2) ** 2

             sigma%dvals_de0ks(:, in_k, ikcalc) = sigma%dvals_de0ks(:, in_k, ikcalc) + gkq2 * rfact_t

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
               cfact_t(:) = - weight_q * gdw2 * (two * nqnu + one)  / (ediff + sigma%ieta)
             else
               cfact_t(:) = zero
             endif
             sigma%dw_vals(:, in_k, ikcalc) = sigma%dw_vals(:, in_k, ikcalc) + real(cfact_t)
             sigma%vals_e0ks(:, in_k, ikcalc) = sigma%vals_e0ks(:, in_k, ikcalc) + real(cfact_t)
           end do ! in_k

         end do ! im_kq
       end do ! my_ip
     end do ! my_iq
     call lg_myk%free()
   end do ! my_ik

   ! TODO: Sternheimer with KS states.

   call sigma%gather_and_write_results(gstore, gqk, dtset, ebands)
   end associate
 end do ! my_is

 call sigma%free()

 ABI_FREE(nqnu)
 !ABI_FREE(f_nk)
 ABI_FREE(f_mkq)
 ABI_FREE(cfact_t)
 ABI_FREE(rfact_t)
 ABI_FREE(displ_nu_cart)
 ABI_FREE(displ_nu_red)
 ABI_FREE(tpp_red)
 ABI_SFREE(vtrial)
 ABI_FREE(phmodes_skip)

 call gstore%free()

end subroutine gstore_sigeph
!!***

!----------------------------------------------------------------------

!!****f* m_gstore_sigeph/sep_gather_and_write_results
!! NAME
!!  sep_gather_and_write_results
!!
!! FUNCTION
!!  Collect results for a given spin, average results in the degenerate subspace
!!  write results to ab_out and netcdf file
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine sep_gather_and_write_results(sigma, gstore, gqk, dtset, ebands)

!Local variables-------------------------------
 class(sep_t),intent(inout) :: sigma
 type(gstore_t),intent(in) :: gstore
 type(gqk_t),intent(in) :: gqk
 type(ebands_t),intent(in) :: ebands
 type(dataset_type),intent(in) :: dtset
 integer,parameter :: max_ntemp = 50
 integer :: it, in_k, ikcalc, ik_bz, spin, ierr, bstart_k, bstop_k, cnt, ndeg
 integer :: band_k,ik_ibz,ib_val,ib_cond,jj,ideg,ii,iw,nstates, nb_k
 !integer :: nq_ibzk_eff, nelem, imyq, iq_ibz_k, sr_ncid
 logical :: changed !, iwrite
 real(dp) :: ravg,kse,kse_prev,dw,fan0,ks_gap,kse_val,kse_cond,qpe_oms,qpe_oms_val,qpe_oms_cond
 real(dp) :: ravg2 ! invsig2fmts, tau
 complex(dp) :: sig0c,zc,qpe,qpe_prev,qpe_val,qpe_cond,cavg1,cavg2,cavg3,cavg4
 character(len=500) :: this_gtype ! msg
 !integer :: grp_ncid, ncerr
 type(degtab_t) :: degtab
!arrays
 real(dp) :: kcalc(3)
 !integer, allocatable :: recvcounts(:), displs(:), nq_rank(:), kq_symtab(:,:), my_kq_symtab(:,:)
 integer,allocatable :: degblock(:,:)
 real(dp) :: qp_gaps(sigma%ntemp),qpoms_gaps(sigma%ntemp)
 !real(dp),allocatable :: aw(:,:,:), a2few_avg(:,:), gather_srate(:,:,:,:), grp_srate(:,:,:,:)
 real(dp) :: ks_enes(gqk%nb_k), ze0_vals(sigma%ntemp, gqk%nb_k)
 !real(dp) :: gfw_avg(sigma%phmesh_size, 3)
 complex(dp) :: qpoms_enes(sigma%ntemp, gqk%nb_k),qp_enes(sigma%ntemp, gqk%nb_k) ! nb_k
!! *************************************************************************

 spin = gqk%spin; nb_k = gqk%nb_k

 ! Sum partial terms inside qgk%comm.
 call xmpi_sum(sigma%vals_e0ks, gqk%comm%value, ierr)
 call xmpi_sum(sigma%dvals_de0ks, gqk%comm%value, ierr)
 call xmpi_sum(sigma%fan_vals, gqk%comm%value, ierr)
 call xmpi_sum(sigma%dw_vals, gqk%comm%value, ierr)

 this_gtype = "KS"
 if (gstore%gtype == "gwpt" .and. dtset%gstore_gname == "gvals") this_gtype = "GWPT"
 if (gstore%gtype == "gwpt" .and. dtset%gstore_gname == "gvals_ks") this_gtype = "KS"

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
   write(ab_out,"(2a)")" Using g(k,q) of type: ", trim(this_gtype)
   write(ab_out,"(a)")" "
   write(ab_out,"(a)")" "
 end if

 ! Compute QP energies and Gaps (Note that I'm assuming a non-magnetic semiconductor!)
 ib_val = nint(ebands%nelect / (two / ebands%nspinor)); ib_cond = ib_val + 1

 do ikcalc=1,gqk%glob_nk
   ik_bz = gstore%kglob2bz(ikcalc, spin)
   ik_ibz = gstore%kbz2ibz(1, ik_bz)
   kcalc = gstore%kbz(:, ik_bz)

   if (abs(dtset%symsigma) == 1) then
     ! Average self-energy matrix elements in the degenerate subspace.
     ! We will have to average the QP corrections over degenerate states if symsigma=1 is used.
     ! Here we make sure that all the degenerate states are included.
     ! Store also band indices of the degenerate sets, used to average final results.
     bstart_k = gqk%bstart_k; bstop_k = gqk%bstop_k
     call ebands%enclose_degbands(ik_ibz, spin, bstart_k, bstop_k, changed, TOL_EDIFF, degblock=degblock)
     !if (changed) then
     !  ABI_WARNING("Changed")
     !end if
     bstart_k = gqk%bstart_k; bstop_k = gqk%bstop_k

     ! Store band indices used for averaging (shifted by bstart_k)
     ndeg = size(degblock, dim=2)
     ABI_MALLOC(degtab%bids, (ndeg))

     do ii=1,ndeg
       ! Make sure boundaries are within the input nk states.
       ! In principle the nk states should be initialized so that all degenerate states are included.
       degblock(1, ii) = max(degblock(1, ii), bstart_k)
       degblock(2, ii) = min(degblock(2, ii), bstop_k)
       cnt = degblock(2, ii) - degblock(1, ii) + 1
       ABI_MALLOC(degtab%bids(ii)%vals, (cnt))
       degtab%bids(ii)%vals = [(jj, jj= &
         degblock(1, ii) - bstart_k + 1, &
         degblock(2, ii) - bstart_k + 1)]
     end do

     ! Average self-energy matrix elements in the degenerate subspace.
     do ideg=1,size(degtab%bids)
       associate (bids => degtab%bids(ideg)%vals)
       nstates = size(bids)
       do it=1,sigma%ntemp
         ! Average QP(T) and Z(T).
         cavg1 = sum(sigma%vals_e0ks(it, bids(:), ikcalc)) / nstates
         cavg2 = sum(sigma%dvals_de0ks(it, bids(:), ikcalc)) / nstates
         cavg3 = sum(sigma%fan_vals(it, bids(:), ikcalc)) / nstates
         !cavg4 = sum(sigma%fan_stern_vals(it, bids(:), ikcalc)) / nstates
         ravg = sum(sigma%dw_vals(it, bids(:), ikcalc)) / nstates
         !ravg2 = sum(sigma%dw_stern_vals(it, bids(:), ikcalc)) / nstates
         do ii=1,nstates
           sigma%vals_e0ks(it, bids(ii), ikcalc) = cavg1
           sigma%dvals_de0ks(it, bids(ii), ikcalc) = cavg2
           sigma%fan_vals(it, bids(ii), ikcalc) = cavg3
           !sigma%fan_stern_vals(it, bids(ii), ikcalc) = cavg4
           sigma%dw_vals(it, bids(ii), ikcalc) = ravg
           !sigma%dw_stern_vals(it, bids(ii), ikcalc) = ravg2
         end do
       end do ! it
       end associate
     end do ! ideg

     call degtab%free()
     ABI_FREE(degblock)
   end if ! symsigma == +1

   kse_val = huge(one) * tol6; kse_cond = huge(one) * tol6
   qp_enes = huge(one) * tol6; qpoms_enes = huge(one) * tol6
   ks_enes = huge(one) * tol6; ze0_vals = huge(one) * tol6
   ks_gap = -one; qpoms_gaps = -one; qp_gaps = -one

   ! Loop over temperatures.
   do it=1,sigma%ntemp
     ! Write header.
     if (it <= max_ntemp) then
       if (ebands%nsppol == 1) then
         write(ab_out,"(3a,f6.1,a,f8.3)") &
           "K-point: ", trim(ktoa(kcalc)), ", T: ", sigma%kTmesh(it) / kb_HaK, &
           " [K], mu_e: ", sigma%mu_e(it) * Ha_eV
       else
         write(ab_out,"(3a,i1,a,f6.1,a,f8.3)") &
           "K-point: ", trim(ktoa(kcalc)), ", spin: ", spin, ", T: ",sigma%kTmesh(it) / kb_HaK, &
           " [K], mu_e: ", sigma%mu_e(it) * Ha_eV
       end if
       if (sigma%imag_only) then
         write(ab_out,"(a)")"   B    eKS    SE2(eKS)  TAU(eKS)  DeKS"
       else
         write(ab_out,"(a)")"   B    eKS     eQP    eQP-eKS   SE1(eKS)  SE2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP"
       end if
     end if

     ! Loop over band n_k for this k-point and spin.
     do in_k=1,nb_k
       band_k = in_k - bstart_k + 1
       kse = ebands%eig(band_k, ik_ibz, spin)
       ks_enes(in_k) = kse
       sig0c = sigma%vals_e0ks(it, in_k, ikcalc)
       dw = sigma%dw_vals(it, in_k, ikcalc)
       fan0 = real(sig0c) - dw
       ! Compute QP energies with On-the-Mass-Shell approximation and first renormalization i.e. Z(eKS)
       ! TODO: Note that here I use the full Sigma including the imaginary part
       !zc = one / (one - sigma%dvals_de0ks(it, in_k))
       zc = one / (one - real(sigma%dvals_de0ks(it, in_k, ikcalc)))
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
         if (sigma%imag_only) then
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
       if (.not. sigma%imag_only) then
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

 if (sigma%ntemp > max_ntemp) then
   write(ab_out, "(a,i0,a)")" No more than ", max_ntemp, " temperatures are written to the main output file."
   write(ab_out, "(2a)")" Please use SIGEPH.nc file and AbiPy to analyze the results.",ch10
 end if

 !stop; call wrtout(std_out, "gstore_sigeph ended OK")

end subroutine sep_gather_and_write_results
!!***

!----------------------------------------------------------------------

!!****f* m_gstore_sigeph/sep_free
!! NAME
!!  sep_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! SOURCE

subroutine sep_free(sigma)

!Arguments ------------------------------------
 class(sep_t),intent(inout) :: sigma
! *********************************************************************

 ABI_SFREE(sigma%kTmesh)
 ABI_SFREE(sigma%mu_e)
 ABI_SFREE(sigma%vals_e0ks)
 ABI_SFREE(sigma%dvals_de0ks)
 ABI_SFREE(sigma%fan_vals)
 ABI_SFREE(sigma%dw_vals)

end subroutine sep_free
!!***

end module m_gstore_sigeph
!!***
 !!***
