!!****m* ABINIT/m_sigx
!! NAME
!!  m_sigx
!!
!! FUNCTION
!!  Calculate diagonal and off-diagonal matrix elements of the exchange part of the self-energy operator.
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2022 ABINIT group (FB, GMR, VO, LR, RWG, MG, RShaltaf)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_sigx

 use defs_basis
 use m_abicore
 use m_gwdefs
 use m_xmpi
 use m_defs_ptgroups
 use m_errors
 use m_time

 use defs_datatypes,  only : pseudopotential_type, ebands_t
 use m_fstrings,      only : itoa, sjoin, ktoa, ltoa
 use m_hide_blas,     only : xdotc, xgemv
 use m_numeric_tools, only : hermitianize
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_t
 use m_fft_mesh,      only : rotate_FFT_mesh, cigfft
 use m_bz_mesh,       only : kmesh_t, findqg0, littlegroup_t
 use m_gsphere,       only : gsphere_t
 use m_vcoul,         only : vcoul_t
 use m_pawpwij,       only : pawpwff_t, pawpwij_t, pawpwij_init, pawpwij_free, paw_rho_tw_g, paw_cross_rho_tw_g
 use m_paw_pwaves_lmn,only : paw_pwaves_lmn_t
 use m_pawang,        only : pawang_type
 use m_pawtab,        only : pawtab_type
 use m_pawfgrtab,     only : pawfgrtab_type
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy, paw_overlap
 use m_paw_nhat,      only : pawmknhat_psipsi
 use m_paw_sym,       only : paw_symcprj
 use m_wfd,           only : wfdgw_t, wave_t
 use m_sigma,         only : sigma_t, sigma_distribute_bks
 use m_oscillators,   only : rho_tw_g
 use m_esymm,         only : esymm_t, esymm_symmetrize_mels, esymm_failed

 implicit none

 private
!!***

 public :: calc_sigx_me
 public :: sigx_symmetrize   ! Symmetrize Sig_x matrix elements
!!***

contains
!!***

!!****f* ABINIT/calc_sigx_me
!! NAME
!! calc_sigx_me
!!
!! FUNCTION
!! Calculate diagonal and off-diagonal matrix elements of the exchange part of the self-energy operator.
!!
!! INPUTS
!! sigmak_ibz=Index of the k-point in the IBZ.
!! bmin, bmax= min and Max band index for GW correction (for this k-point)
!! Gsph_x<gsphere_t>= Info on the G-sphere used for Sigma_x
!! ikcalc=index in the array Sigp%kptgw2bz of the k-point where GW corrections are calculated
!! ltg_k datatype containing information on the little group
!! Kmesh <kmesh_t>
!! x_ngfft(18)=Information about 3D FFT for the oscillator strengths, see ~abinit/doc/variables/vargs.htm#ngfft
!! Vcp <vcoul_t datatype> containing information on the cutoff technique
!! Pawtab(psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Pawang <type(pawang_type)>=paw angular mesh and related data
!! psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! Qmesh <kmesh_t> : datatype gathering information of the q-mesh used
!! Sigp <sigparams_t> (see the definition of this structured datatype)
!! cryst<crystal_t>=Info on unit cell and symmetries
!! qp_ebands<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!! Paw_pwff<pawpwff_t>=Form factor used to calculate the onsite mat. elements of a plane wave.
!! allQP_sym(%nkibz, %nsppol)<esymm_t>=Datatype collecting data on the irreducible representations of the
!!    little group of kcalc in the KS representation as well as the symmetry of the bdgw_k states.
!! prtvol=Flags governing verbosity level.
!!
!! OUTPUT
!!  Sr%x_mat(bmin:bmax,bmin:bmax,%nsppol*Sigp%nsig_ab)=Matrix elements of Sigma_x.
!!
!! NOTES
!!  1) The treatment of the divergence of Gygi+Baldereschi (PRB 1986) [[cite:Gigy1986]] is included.
!!
!!  2) On the symmetrization of Sigma matrix elements
!!     If  Sk = k+G0 then  M_G(k, Sq)= e^{-i( Sq+G).t} M_{ S^-1(G}   (k,q)
!!     If -Sk = k+G0 then  M_G(k,-Sq)= e^{-i(-Sq+G).t} M_{-S^-1(G)}^*(k,q)
!!
!! Notice the absence of G0 in the expression. Moreover, when we sum over the little group, it turns out
!! that there is a cancellation of the phase factor associated to the non-symmorphic operations due to a
!! similar term coming from the symmetrization of \epsilon^{-1}. Mind however that the nonsymmorphic phase
!! has to be considered when epsilon^-1 is reconstructed starting from the q-points in the IBZ.
!!
!!  3) the unitary transformation relating wavefunctions
!!     at symmetric k-points should be taken into account during the symmetrization
!!     of the oscillator matrix elements. In case of G_oW_o and GW_o calculations, however,
!!     it is possible to make an invariant by just including all the degenerate states and
!!     averaging the final results over the degenerate subset. Here we divide the states
!!     where the QP energies are required into complexes. Note however that this approach is not
!!     based on group theory, and it might lead to spurious results in case of accidental degeneracies.
!!

subroutine calc_sigx_me(sigmak_ibz, ikcalc, bmin, bmax, cryst, qp_ebands, Sigp, Sr, Gsph_x, Vcp, Kmesh, Qmesh, &
                        ltg_k, Pawtab, Pawang, Paw_pwff, Pawfgrtab, Paw_onsite, psps, wfd, Wfdf, &
                        allQP_sym, x_ngfft, ngfftf, prtvol, pawcross, tol_empty_in)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: sigmak_ibz,ikcalc,prtvol,bmin,bmax,pawcross
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),target,intent(in) :: qp_ebands
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(vcoul_t),intent(in) :: Vcp
 type(gsphere_t),intent(in) :: Gsph_x
 type(littlegroup_t),intent(in) :: ltg_k
 type(Pseudopotential_type),intent(in) :: psps
 type(sigparams_t),target,intent(in) :: Sigp
 type(sigma_t),intent(inout) :: Sr
 type(pawang_type),intent(in) :: Pawang
 type(wfdgw_t),target,intent(inout) :: wfd,Wfdf
 real(dp),intent(in) :: tol_empty_in
!arrays
 integer,intent(in) :: x_ngfft(18),ngfftf(18)
 type(Pawtab_type),intent(in) :: Pawtab(psps%ntypat)
 type(pawpwff_t),intent(in) :: Paw_pwff(psps%ntypat*psps%usepaw)
 type(esymm_t),target,intent(in) :: allQP_sym(wfd%nkibz, wfd%nsppol)
 type(pawfgrtab_type),intent(inout) :: Pawfgrtab(cryst%natom*psps%usepaw)
 type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(cryst%natom*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: ndat1 = 1, use_pawnhat0 = 0, ider0 = 0
 integer :: gwcalctyp,izero,iab,band_sum,ib,ierr,ig,ig_rot,ii,iik,itim_q,i2
 integer :: ik_bz, ik_ibz, isym_q, iq_bz, iq_ibz, spin, isym, jb, is_idx
 integer :: jik,jk_bz,jk_ibz,kb,nspinor,nsppol,ifft
 integer :: nq_summed,ibsp,dimcprj_gw,dim_rtwg, isym_kgw, isym_ki
 integer :: spad, spadx1, spadx2, irow, npw_k, wtqm, wtqp
 integer :: npwx, x_nfft, x_mgfft, x_fftalga, nsig_ab
 integer :: nfftf, mgfftf, nhat12_grdim, my_nbks, use_padfft, use_padfftf
 real(dp) :: cpu, wall, gflops, fact_spin, theta_mu_minus_esum, theta_mu_minus_esum2, tol_empty
 complex(dpc) :: ctmp,ph_mkgwt,ph_mkt
 complex(gwpc) :: gwpc_sigxme,gwpc_sigxme2,xdot_tmp
 logical :: iscompatibleFFT,q_is_gamma
 character(len=5000) :: msg
 type(wave_t),pointer :: wave_sum, wave_jb
!arrays
 integer :: g0(3), spinor_padx(2,4)
 integer,allocatable :: igfftxg0(:), igfftfxg0(:), x_gbound(:,:), gboundf(:,:)
 integer,allocatable :: ktabr(:,:),irottb(:,:),ktabrf(:,:), proc_distrb(:,:,:)
 real(dp) :: ksum(3), kgw(3), kgw_m_ksum(3), qbz(3), q0(3), spinrot_kbz(4), spinrot_kgw(4), tsec(2)
 real(dp),contiguous, pointer :: qp_ene(:,:,:), qp_occ(:,:,:)
 real(dp),allocatable :: nhat12(:,:,:),grnhat12(:,:,:,:)
 complex(gwpc),allocatable :: vc_sqrt_qbz(:), rhotwg(:), rhotwgp(:), rhotwg_ki(:,:), ur_bdgw(:,:), ur_ibz(:)
 complex(dpc),allocatable  :: sigxcme_tmp(:,:), sigxme_tmp(:,:,:), sigx(:,:,:,:)
 complex(gwpc),allocatable :: ur_ae_sum(:),ur_ae_onsite_sum(:),ur_ps_onsite_sum(:)
 complex(gwpc),allocatable :: ur_ae_bdgw(:,:),ur_ae_onsite_bdgw(:,:),ur_ps_onsite_bdgw(:,:)
 complex(gwpc),contiguous, pointer :: cg_jb(:),cg_sum(:)
 logical :: can_symmetrize(wfd%nsppol)
 logical,allocatable :: bks_mask(:,:,:)
 type(esymm_t),pointer :: QP_sym(:)
 type(sigijtab_t),pointer :: Sigxij_tab(:)
 type(pawcprj_type),allocatable :: Cprj_kgw(:,:),Cprj_ksum(:,:)
 type(pawpwij_t),allocatable :: Pwij_qg(:),Pwij_fft(:)

!************************************************************************

 DBG_ENTER("COLL")

 call timab(430,1,tsec) ! csigme (SigX)
 call cwtime(cpu, wall, gflops, "start")

 ! Initialize some values.
 gwcalctyp = Sigp%gwcalctyp; nspinor = wfd%nspinor; nsppol = wfd%nsppol; npwx = sigp%npwx
 dim_rtwg = 1; if (nspinor == 2) dim_rtwg = 2
 nsig_ab = sigp%nsig_ab; spinor_padx = reshape([0, 0, npwx, npwx, 0, npwx, npwx, 0], [2, 4])
 ABI_CHECK(Sigp%npwx == Gsph_x%ng, "Sigp%npwx != Gsph_x%ng")

 qp_ene => qp_ebands%eig; qp_occ => qp_ebands%occ

 ! Exctract the symmetries of the bands for this k-point
 QP_sym => allQP_sym(sigmak_ibz,1:nsppol)

 ! Index of sigma_k k-point in the BZ array, its image in the IBZ and symmetries
 jk_bz = Sigp%kptgw2bz(ikcalc)
 call kmesh%get_BZ_item(jk_bz, kgw, jk_ibz, isym_kgw, jik, ph_mkgwt)
 spinrot_kgw(:) = cryst%spinrot(:,isym_kgw)

 write(msg,'(6a)') ch10, &
  ' Calculating <nk|Sigma_x|nk> at k: ',trim(ktoa(kgw)), ", for bands: ", trim(ltoa([bmin, bmax])),ch10
 call wrtout(std_out, msg)

 if (any(x_ngfft(1:3) /= wfd%ngfft(1:3)) ) call wfd%change_ngfft(cryst, psps, x_ngfft)
 x_nfft = product(x_ngfft(1:3)); x_mgfft = maxval(x_ngfft(1:3)); x_fftalga = x_ngfft(7) / 100

 if (pawcross==1) mgfftf = MAXVAL(ngfftf(1:3))

 ! Define whether we can use symmetries to sum over the IBZ_kgw
 can_symmetrize = .FALSE.
 if (Sigp%symsigma > 0) then
   can_symmetrize = .TRUE.
   if (gwcalctyp >= 20) then
     do spin=1,nsppol
       can_symmetrize(spin) = .not. esymm_failed(QP_sym(spin))
       if (.not.can_symmetrize(spin)) then
         write(msg,'(a,i0,4a)')&
          "Symmetrization cannot be performed for spin: ",spin,ch10,&
          "band classification encountered the following problem: ",ch10,TRIM(QP_sym(spin)%err_msg)
         ABI_WARNING(msg)
       end if
     end do
   end if
   if (nspinor == 2) then
     ABI_WARNING('Symmetrization with nspinor=2 not implemented')
   end if
 end if

 ! MRM allow lower occ numbers
 ! Normalization of theta_mu_minus_esum. If nsppol==2, qp_occ $\in [0,1]$
 select case (nsppol)
 case (1)
   fact_spin = half; tol_empty = tol_empty_in          ! below this value the state is assumed empty
   if (nspinor == 2) then
     fact_spin = one; tol_empty = half * tol_empty_in  ! below this value the state is assumed empty
   end if
 case (2)
   fact_spin = one; tol_empty = half * tol_empty_in  ! to be consistent and obtain similar results if a metallic
 case default                                        ! spin unpolarized system is treated using nsppol==2
   ABI_BUG(sjoin('Wrong nsppol:', itoa(nsppol)))
 end select

 ! Table for \Sigmax_ij matrix elements.
 Sigxij_tab => Sigp%Sigxij_tab(ikcalc, 1:nsppol)

 ! Remove empty states from the list of states that will be distributed.
 ABI_MALLOC(bks_mask, (wfd%mband, Kmesh%nbz, nsppol))
 bks_mask = .FALSE.

 do spin=1,nsppol
   do ik_bz=1,Kmesh%nbz
     ik_ibz = Kmesh%tab(ik_bz)
     do band_sum=1,Sigp%nbnds
       bks_mask(band_sum, ik_bz, spin) = (abs(qp_occ(band_sum, ik_ibz, spin)) >= tol_empty)  ! MRM allow negative occ
     end do
   end do
 end do

 ! Distribute tasks.
 ABI_MALLOC(proc_distrb, (wfd%mband, Kmesh%nbz, nsppol))
 call sigma_distribute_bks(wfd,Kmesh,ltg_k,Qmesh,nsppol,can_symmetrize,kgw,Sigp%mg0,my_nbks,proc_distrb,bks_mask=bks_mask)
 ABI_FREE(bks_mask)
 call wrtout(std_out, sjoin(" Will sum ", itoa(my_nbks) ," (b, k, s) occupied states in Sigma_x."))

 ! The index of G-G0 in the FFT mesh for the oscillators
 ! Sigp%mG0 gives the MAX G0 component to account for umklapp.
 ABI_MALLOC(igfftxg0, (Gsph_x%ng))

 ! Precompute the FFT index of $ R^{-1}(r-\tau)$
 ! S = \transpose R^{-1} and k_BZ = S k_IBZ
 ! irottb is the FFT index of $R^{-1} (r-\tau)$ used to symmetrize u_Sk.

 ABI_MALLOC(irottb, (x_nfft, cryst%nsym))
 call rotate_FFT_mesh(cryst%nsym, cryst%symrel, cryst%tnons, x_ngfft, irottb, iscompatibleFFT)
 if (.not. iscompatibleFFT) then
   ABI_WARNING("FFT mesh is not compatible with symmetries. Results might be affected by large errors!")
 end if

 ABI_MALLOC(ktabr, (x_nfft, Kmesh%nbz))
 do ik_bz=1,Kmesh%nbz
   isym = Kmesh%tabo(ik_bz)
   do ifft=1,x_nfft
     ktabr(ifft,ik_bz) = irottb(ifft,isym)
   end do
 end do
 ABI_FREE(irottb)

 if (psps%usepaw == 1 .and. pawcross == 1) then
   nfftf = PRODUCT(ngfftf(1:3))
   ABI_MALLOC(irottb, (nfftf, cryst%nsym))
   call rotate_FFT_mesh(cryst%nsym, cryst%symrel, cryst%tnons, ngfftf, irottb, iscompatibleFFT)

   ABI_MALLOC(ktabrf,(nfftf, Kmesh%nbz))
   do ik_bz=1,Kmesh%nbz
     isym=Kmesh%tabo(ik_bz)
     do ifft=1,nfftf
       ktabrf(ifft,ik_bz)=irottb(ifft,isym)
     end do
   end do
   ABI_FREE(irottb)
 end if

 ! Additional allocations for PAW.
 if (psps%usepaw == 1) then
   ABI_MALLOC(Cprj_ksum, (cryst%natom, nspinor))
   call pawcprj_alloc(Cprj_ksum, 0, wfd%nlmn_atm)

   nhat12_grdim = 0
   if (use_pawnhat0 == 1) then
     ! Compensation charge for \phi_a^*\phi_b
     call wrtout(std_out, "Using nhat12")
     ABI_MALLOC(nhat12  ,(2, x_nfft, nspinor**2))
     ABI_MALLOC(grnhat12,(2, x_nfft, nspinor**2, 3*nhat12_grdim))
   end if
 end if

 nq_summed = Kmesh%nbz
 if (Sigp%symsigma > 0) then
   call ltg_k%print(std_out, prtvol, mode_paral='COLL')
   nq_summed = sum(ltg_k%ibzq(:))
 end if ! symsigma

 write(msg,'(2a,i0,a)')ch10,' calc_sigx_me: calculation status (', nq_summed, ' to be completed):'
 call wrtout(std_out, msg)

 ABI_MALLOC(ur_ibz, (x_nfft * nspinor))
 ABI_MALLOC(rhotwg_ki, (npwx * nspinor, bmin:bmax))
 ABI_MALLOC(rhotwg, (npwx * nspinor))
 ABI_MALLOC(rhotwgp, (npwx * nspinor))
 ABI_MALLOC(vc_sqrt_qbz, (npwx))

 ABI_CALLOC(sigxme_tmp, (bmin:bmax, bmin:bmax, nsppol * nsig_ab))
 ABI_CALLOC(sigxcme_tmp, (bmin:bmax, nsppol * nsig_ab))
 ABI_CALLOC(sigx, (2, bmin:bmax, bmin:bmax, nsppol * nsig_ab))

 if (pawcross==1) then
   ABI_MALLOC(ur_ae_sum,(nfftf*nspinor))
   ABI_MALLOC(ur_ae_onsite_sum,(nfftf*nspinor))
   ABI_MALLOC(ur_ps_onsite_sum,(nfftf*nspinor))
 end if

 do spin=1,nsppol
   if (ALL(proc_distrb(:,:,spin) /= wfd%my_rank)) CYCLE ! Spin parallelism.

   ! ===============================================
   ! Load wavefunctions for Sigma_x matrix elements
   ! ===============================================
   ABI_MALLOC_OR_DIE(ur_bdgw, (x_nfft * nspinor, bmin:bmax), ierr)
   call wfd%get_many_ur([(jb, jb=bmin, bmax)], jk_ibz, spin, ur_bdgw)

   if (wfd%usepaw == 1) then
     ! Load cprj for GW states, note the indexing.
     dimcprj_gw = nspinor * (bmax - bmin + 1)
     ABI_MALLOC(Cprj_kgw, (cryst%natom, bmin:bmin+dimcprj_gw-1))
     call pawcprj_alloc(Cprj_kgw, 0, wfd%nlmn_atm)
     ibsp = bmin
     do jb=bmin,bmax
       call wfd%get_cprj(jb, jk_ibz, spin, cryst, Cprj_ksum, sorted=.FALSE.)
       call paw_symcprj(jk_bz, nspinor, 1, cryst, Kmesh, Pawtab, Pawang, Cprj_ksum)
       call pawcprj_copy(Cprj_ksum, Cprj_kgw(:,ibsp:ibsp+(nspinor-1)))
       ibsp = ibsp + nspinor
     end do
     if (pawcross ==1) then
       ABI_MALLOC(ur_ae_bdgw,(nfftf*nspinor,bmin:bmax))
       ABI_MALLOC(ur_ae_onsite_bdgw,(nfftf*nspinor,bmin:bmax))
       ABI_MALLOC(ur_ps_onsite_bdgw,(nfftf*nspinor,bmin:bmax))
       do jb=bmin,bmax
         call wfdf%paw_get_aeur(jb,jk_ibz,spin,cryst,Paw_onsite,psps,Pawtab,Pawfgrtab,&
                                ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum)
         ur_ae_bdgw(:,jb)=ur_ae_sum
         ur_ae_onsite_bdgw(:,jb)=ur_ae_onsite_sum
         ur_ps_onsite_bdgw(:,jb)=ur_ps_onsite_sum
       end do
     end if
   end if

   ! ==============================
   ! ==== Sum over k in the BZ ====
   ! ==============================

   do ik_bz=1,Kmesh%nbz

     ! Parallelization over k-points and spin.
     if (ALL(proc_distrb(:,ik_bz,spin) /= wfd%my_rank)) CYCLE

     ! Find the symmetrical image of ksum in the IBZ
     call kmesh%get_BZ_item(ik_bz, ksum, ik_ibz, isym_ki, iik, ph_mkt)
     spinrot_kbz = cryst%spinrot(:,isym_ki)

     ! Identify q and G0 where q + G0 = k_GW - ksum
     kgw_m_ksum = kgw - ksum
     call findqg0(iq_bz, g0, kgw_m_ksum, Qmesh%nbz, Qmesh%bz, Sigp%mG0)

     ! If symmetries are exploited only q-points in the IBZ_k are computed.
     ! In this case elements are weighted according to wtqp and wtqm. wtqm is for time-reversal.
     wtqp = 1; wtqm = 0
     if (can_symmetrize(spin)) then
       if (ltg_k%ibzq(iq_bz) /= 1) cycle
       wtqp = 0; wtqm = 0
       do isym=1,ltg_k%nsym_sg
         wtqp = wtqp + ltg_k%wtksym(1, isym, iq_bz)
         wtqm = wtqm + ltg_k%wtksym(2, isym, iq_bz)
       end do
     end if

     write(msg,'(2(a,i4),a,i3)')' calc_sigx_me: ik_bz ',ik_bz,'/',Kmesh%nbz,' done by mpi-rank: ',wfd%my_rank
     call wrtout(std_out, msg)

     ! Find the corresponding irreducible q-point.
     ! NB: non-zero umklapp G_o is not allowed. There's a check in setup_sigma
     call qmesh%get_BZ_item(iq_bz, qbz, iq_ibz, isym_q, itim_q)
     q_is_gamma = normv(qbz, cryst%gmet, "G") < GW_TOLQ0

     ! Tables for the FFT of the oscillators.
     !  a) FFT index of G-G0.
     !  b) x_gbound table for the zero-padded FFT performed in rhotwg.
     ABI_MALLOC(x_gbound, (2*x_mgfft+8, 2))
     call Gsph_x%fft_tabs(g0, x_mgfft, x_ngfft, use_padfft, x_gbound, igfftxg0)

     if (any(x_fftalga == [2, 4])) use_padfft = 0 ! Padded-FFT is not coded in rho_tw_g
#ifdef FC_IBM
     ! XLF does not deserve this optimization (problem with [v67mbpt][t03])
     use_padfft = 0
#endif
     if (use_padfft == 0) then
       ABI_FREE(x_gbound)
       ABI_MALLOC(x_gbound, (2*x_mgfft+8, 2*use_padfft))
     end if

     if (pawcross==1) then
       ABI_MALLOC(gboundf,(2*mgfftf+8,2))
       ABI_MALLOC(igfftfxg0,(Gsph_x%ng))
       call Gsph_x%fft_tabs(g0,mgfftf,ngfftf,use_padfftf,gboundf,igfftfxg0)
       if ( ANY(x_fftalga == [2, 4]) ) use_padfftf=0
       if (use_padfftf==0) then
         ABI_FREE(gboundf)
         ABI_MALLOC(gboundf,(2*mgfftf+8,2*use_padfftf))
       end if
     end if

     if (psps%usepaw==1 .and. use_pawnhat0 == 0) then
       ! Evaluate oscillator matrix elements
       ! $ <phj/r|e^{-i(q+G)}|phi/r> - <tphj/r|e^{-i(q+G)}|tphi/r> $ in packed form
       q0 = qbz !;if (q_is_gamma) q0 = (/0.00001_dp,0.00001_dp,0.00001_dp/) ! GW_Q0_DEFAULT
       ABI_MALLOC(Pwij_qg, (psps%ntypat))
       call pawpwij_init(Pwij_qg, npwx, q0, Gsph_x%gvec, cryst%rprimd, psps, Pawtab, Paw_pwff)
     end if

     ! Get Fourier components of the Coulomb interaction in the BZ
     ! In 3D systems, neglecting umklapp,  vc(Sq,sG)=vc(q,G)=4pi/|q+G|
     ! The same relation holds for 0-D systems, but not in 1-D or 2D systems. It depends on S.
     do ig=1,npwx
       ig_rot = Gsph_x%rottb(ig, itim_q, isym_q)
       vc_sqrt_qbz(ig_rot) = Vcp%vc_sqrt_resid(ig, iq_ibz)
     end do

     ! ==========================
     ! Sum over (occupied) bands
     ! ==========================
     do band_sum=1,Sigp%nbnds

       ! Parallelism over bands.
       if (proc_distrb(band_sum, ik_bz, spin) /= wfd%my_rank) CYCLE

       ! Skip empty states. MRM: allow negative occ numbers.
       if (abs(qp_occ(band_sum, ik_ibz, spin)) < tol_empty) CYCLE

       call wfd%get_ur(band_sum, ik_ibz, spin, ur_ibz)

       if (psps%usepaw == 1) then
         ! Load cprj for point ksum, this spin or spinor and *THIS* band.
         ! TODO MG I could avoid doing this but I have to exchange spin and bands ???
         ! For sure there is a better way to do this!
         call wfd%get_cprj(band_sum, ik_ibz, spin, cryst, Cprj_ksum, sorted=.FALSE.)
         call paw_symcprj(ik_bz, nspinor, 1, cryst, Kmesh, Pawtab, Pawang, Cprj_ksum)
         if (pawcross==1) then
           call wfdf%paw_get_aeur(band_sum,ik_ibz,spin,cryst,Paw_onsite,psps,Pawtab,Pawfgrtab,&
                                  ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum)
         end if
       end if

       ! Get all <k-q,band_sum,s|e^{-i(q+G).r}|s,jb,k>
       do jb=bmin,bmax

         if (Psps%usepaw==1 .and. use_pawnhat0 == 1) then
           ABI_ERROR("use_pawnhat is disabled")
           i2=jb; if (nspinor==2) i2=(2*jb-1)
           spad = nspinor - 1

           izero=0
           call pawmknhat_psipsi(Cprj_ksum,Cprj_kgw(:,i2:i2+spad),ider0,izero,cryst%natom,&
                                 cryst%natom,x_nfft,x_ngfft,nhat12_grdim,nspinor,cryst%ntypat,Pawang,Pawfgrtab,&
                                 grnhat12,nhat12,pawtab)

         else
           call rho_tw_g(nspinor,npwx,x_nfft,ndat1,x_ngfft,1,use_padfft,igfftxg0,x_gbound, &
                         ur_ibz        ,iik,ktabr(:,ik_bz),ph_mkt  ,spinrot_kbz, &
                         ur_bdgw(:,jb),jik,ktabr(:,jk_bz),ph_mkgwt,spinrot_kgw, &
                         nspinor,rhotwg_ki(:,jb))

           if (psps%usepaw == 1 .and. use_pawnhat0 == 0) then
             ! Add on-site contribution, projectors are already in BZ.
             i2=jb; if (nspinor==2) i2=(2*jb-1)
             spad = nspinor - 1
             call paw_rho_tw_g(npwx,nspinor,nspinor,cryst%natom,cryst%ntypat,cryst%typat,cryst%xred,Gsph_x%gvec,&
                               Cprj_ksum(:,:),Cprj_kgw(:,i2:i2+spad),Pwij_qg,rhotwg_ki(:,jb))
           end if
           if (psps%usepaw==1.and.pawcross==1) then ! Add paw cross term
             call paw_cross_rho_tw_g(nspinor,npwx,nfftf,ngfftf,1,use_padfftf,igfftfxg0,gboundf,&
               ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum,iik,ktabrf(:,ik_bz),ph_mkt,spinrot_kbz,&
               ur_ae_bdgw(:,jb),ur_ae_onsite_bdgw(:,jb),ur_ps_onsite_bdgw(:,jb),jik,ktabrf(:,jk_bz),ph_mkgwt,spinrot_kgw,&
               nspinor,rhotwg_ki(:,jb))
           end if
         end if

         ! Multiply by the square root of the Coulomb term
         ! In 3-D systems, the factor sqrt(4pi) is included
         do ii=1,nspinor
           spad = (ii-1) * npwx
           rhotwg_ki(spad+1:spad+npwx, jb) = rhotwg_ki(spad+1:spad + npwx, jb) * vc_sqrt_qbz(1:npwx)
         end do

         if (ik_bz == jk_bz) then
           ! Treat analytically the case q --> 0:
           !
           !   * The oscillator is evaluated at q = 0 as it is considered constant in the small cube around Gamma
           !     while the Colulomb term is integrated out.
           !   * If nspinor == 1, we have nonzero contribution only if band_sum == jb
           !   * If nspinor == 2, we evaluate <band_sum,up|jb,up> and <band_sum,dwn|jb,dwn>,
           !     and impose orthonormalization since npwwfn might be < npwvec.
           !   * Note the use of i_sz_resid and not i_sz, to account for the possibility
           !     to have generalized KS basis set from hybrid

           if (nspinor == 1) then
             rhotwg_ki(1, jb) = czero_gw
             if (band_sum == jb) rhotwg_ki(1,jb) = cmplx(sqrt(Vcp%i_sz_resid), 0.0_gwp)
             !rhotwg_ki(1,jb) = czero_gw ! DEBUG

           else
             npw_k = wfd%npwarr(ik_ibz)
             rhotwg_ki(1, jb) = zero; rhotwg_ki(npwx+1, jb) = zero
             if (band_sum == jb) then
               ABI_CHECK(wfd%get_wave_ptr(band_sum, ik_ibz, spin, wave_sum, msg) == 0, msg)
               cg_sum => wave_sum%ug
               ABI_CHECK(wfd%get_wave_ptr(jb, jk_ibz, spin, wave_jb, msg) == 0, msg)
               cg_jb  => wave_jb%ug
               ctmp = xdotc(npw_k, cg_sum(1:), 1, cg_jb(1:), 1)
               rhotwg_ki(1, jb) = cmplx(sqrt(Vcp%i_sz_resid), 0.0_gwp) * real(ctmp)
               ctmp = xdotc(npw_k, cg_sum(npw_k+1:), 1, cg_jb(npw_k+1:), 1)
               rhotwg_ki(npwx+1, jb) = cmplx(sqrt(Vcp%i_sz_resid), 0.0_gwp) * real(ctmp)
             end if
             !rhotwg_ki(1, jb) = zero; rhotwg_ki(npwx+1, jb) = zero
             ! PAW is missing
           end if
         end if
       end do ! jb Got all matrix elements from bmin up to bmax.

       theta_mu_minus_esum  = fact_spin * qp_occ(band_sum, ik_ibz, spin)
       theta_mu_minus_esum2 = sqrt(abs(fact_spin * qp_occ(band_sum, ik_ibz, spin))) ! MBB Nat. orb. funct. approx. sqrt(occ)

       if (abs(theta_mu_minus_esum / fact_spin) >= tol_empty) then     ! MRM: allow negative occ numbers
         do kb=bmin,bmax
           ! Copy the ket Sigma_x |phi_{k,kb}>.
           rhotwgp(:) = rhotwg_ki(:, kb)

           ! Loop over the non-zero row elements of this column.
           ! If gwcalctyp <  20: only diagonal elements since QP == KS.
           ! If gwcalctyp >= 20:
           !      * Only off-diagonal elements connecting states with same character.
           !      * Only the upper triangle if HF, SEX, or COHSEX.

           do irow=1,Sigxij_tab(spin)%col(kb)%size1
             jb = Sigxij_tab(spin)%col(kb)%bidx(irow)
             rhotwg(:) = rhotwg_ki(:,jb)

             ! Calculate bare exchange <phi_jb|Sigma_x|phi_kb>.
             ! Do the scalar product only if band_sum is occupied.
             do iab=1,nsig_ab
               spadx1 = spinor_padx(1, iab); spadx2 = spinor_padx(2, iab)
               xdot_tmp = -XDOTC(npwx, rhotwg(spadx1+1:), 1, rhotwgp(spadx2+1:), 1)
               gwpc_sigxme  = xdot_tmp * theta_mu_minus_esum
               gwpc_sigxme2 = xdot_tmp * theta_mu_minus_esum2

               ! Accumulate and symmetrize Sigma_x matrix elements.
               ! -wtqm comes from time-reversal (exchange of band indeces)
               is_idx = spin; if (nspinor == 2) is_idx = iab
               sigxme_tmp(jb, kb, is_idx) = sigxme_tmp(jb, kb, is_idx) + &
                  (wtqp + wtqm) * DBLE(gwpc_sigxme) + (wtqp - wtqm) * j_gw * AIMAG(gwpc_sigxme)
               if (jb == kb) then
                 sigxcme_tmp(jb, is_idx) = sigxcme_tmp(jb, is_idx) + &
                   (wtqp + wtqm) * DBLE(gwpc_sigxme2) + (wtqp - wtqm) *j_gw * AIMAG(gwpc_sigxme2)
               end if

               sigx(1, jb, kb, is_idx) = sigx(1, jb, kb, is_idx) + wtqp *      gwpc_sigxme
               sigx(2, jb, kb, is_idx) = sigx(2, jb, kb, is_idx) + wtqm *CONJG(gwpc_sigxme)
             end do
           end do ! jb
         end do ! kb
       end if

     end do ! band_sum

     ! Deallocate k-dependent quantities.
     ABI_FREE(x_gbound)
     if (pawcross==1) then
       ABI_FREE(gboundf)
     end if

     if (psps%usepaw==1 .and. use_pawnhat0 == 0) then
       call pawpwij_free(Pwij_qg)
       ABI_FREE(Pwij_qg)
     end if
   end do ! ik_bz Got all diagonal (off-diagonal) matrix elements.

   ABI_FREE(ur_bdgw)
   if (wfd%usepaw == 1) then
     call pawcprj_free(Cprj_kgw)
     ABI_FREE(Cprj_kgw)
     if (pawcross==1) then
       ABI_FREE(ur_ae_bdgw)
       ABI_FREE(ur_ae_onsite_bdgw)
       ABI_FREE(ur_ps_onsite_bdgw)
     end if
   end if
 end do !spin

 ABI_FREE(igfftxg0)
 if (pawcross==1) then
   ABI_FREE(igfftfxg0)
 end if

 ! Gather contributions from all the CPUs.
 call xmpi_sum(sigxme_tmp, wfd%comm, ierr)
 call xmpi_sum(sigxcme_tmp, wfd%comm, ierr)
 call xmpi_sum(sigx, wfd%comm, ierr)

 ! Multiply by constants. For 3D systems sqrt(4pi) is included in vc_sqrt_qbz.
 sigxme_tmp  = (one / (cryst%ucvol * Kmesh%nbz)) * sigxme_tmp  * Sigp%sigma_mixing
 sigxcme_tmp = (one / (cryst%ucvol * Kmesh%nbz)) * sigxcme_tmp * Sigp%sigma_mixing
 sigx        = (one / (cryst%ucvol * Kmesh%nbz)) * sigx        * Sigp%sigma_mixing

 ! If we have summed over the IBZ_q, we have to average over degenerate states.
 ! NOTE: Presently only diagonal terms are considered
 ! TODO QP-SCGW required a more involved approach, there is a check in sigma
 ! TODO it does not work if spinor == 2.

 do spin=1,nsppol
   if (.not. can_symmetrize(spin)) cycle
   call sigx_symmetrize(jk_ibz, spin, bmin, bmax, nsppol, nspinor, nsig_ab, qp_ene, sigx, sigxme_tmp)
 end do

 if (gwcalctyp >= 20) then
   ! Reconstruct the full sigma_x matrix from the upper triangle.
   if (nspinor == 1) then
     do spin=1,nsppol
       call hermitianize(sigxme_tmp(:,:,spin), "Upper")
     end do
   else
     ABI_WARNING("Should hermitianize non-collinear sigma!")
   end if
 end if

 ! Save diagonal elements or ab components of Sigma_x (Hermitian)
 ! TODO It should be hermitian also if nspinor == 2
 do spin=1,nsppol
   do jb=bmin,bmax
     do iab=1,nsig_ab
       is_idx = spin; if (nsig_ab > 1) is_idx = iab
       if (is_idx <= 2) then
         Sr%sigxme(jb,jk_ibz,is_idx)     = DBLE( sigxme_tmp(jb,jb,is_idx))
         Sr%sigxcnofme(jb,jk_ibz,is_idx) = DBLE(sigxcme_tmp(jb,is_idx))
       else
         Sr%sigxme(jb,jk_ibz,is_idx)     =  sigxme_tmp(jb,jb,is_idx)
         Sr%sigxcnofme(jb,jk_ibz,is_idx) = sigxcme_tmp(jb,is_idx)
       end if
     end do
     !if (nsig_ab > 1) then
     !  write(std_out,'(i3,4f8.3,a,f8.3)')jb,Sr%sigxme(jb,jk_ibz,:)*Ha_eV,' Tot ',SUM(Sr%sigxme(jb,jk_ibz,:))*Ha_eV
     !end if
   end do
 end do

 ! Save full exchange matrix in Sr%
 Sr%x_mat(bmin:bmax, bmin:bmax, jk_ibz, :) = sigxme_tmp(bmin:bmax, bmin:bmax,:)
 ABI_FREE(sigxme_tmp)
 ABI_FREE(sigxcme_tmp)

 ! ===========================
 ! ==== Deallocate memory ====
 ! ===========================
 if (psps%usepaw == 1) then
   call pawcprj_free(Cprj_ksum)
   ABI_FREE(Cprj_ksum)
   if (allocated(Pwij_fft)) then
     call pawpwij_free(Pwij_fft)
     ABI_FREE(Pwij_fft)
   end if
   if (use_pawnhat0 == 1) then
     ABI_FREE(nhat12)
     ABI_FREE(grnhat12)
   end if
   if (pawcross == 1) then
     ABI_FREE(ur_ae_sum)
     ABI_FREE(ur_ae_onsite_sum)
     ABI_FREE(ur_ps_onsite_sum)
     ABI_FREE(ktabrf)
   end if
 end if

 ABI_FREE(ur_ibz)
 ABI_FREE(rhotwg_ki)
 ABI_FREE(rhotwg)
 ABI_FREE(rhotwgp)
 ABI_FREE(vc_sqrt_qbz)
 ABI_FREE(ktabr)
 ABI_FREE(sigx)
 ABI_FREE(proc_distrb)

 call timab(430,2,tsec) ! csigme (SigX)
 call cwtime_report(" calc_sigx_me:", cpu, wall, gflops)

 DBG_EXIT("COLL")

end subroutine calc_sigx_me
!!***

!!****f* ABINIT/sigx_symmetrize
!! NAME
!! sigx_symmetrize
!!
!! FUNCTION
!!  Symmetrize Sig_x matrix elements
!!

subroutine sigx_symmetrize(jk_ibz, spin, bmin, bmax, nsppol, nspinor, nsig_ab, qp_ene, sigx, sigxme_tmp)

 integer,intent(in) :: jk_ibz, spin, bmin, bmax, nsppol, nspinor, nsig_ab
 real(dp),intent(in) :: qp_ene(:,:,:)
 complex(dp),intent(in) :: sigx(2, bmin:bmax, bmin:bmax, nsppol * nsig_ab)
 complex(dp),intent(inout) :: sigxme_tmp(bmin:bmax, bmin:bmax, nsppol * nsig_ab)

!Local variables ------------------------------
 integer :: ib, jb, ndegs, ii
 integer,allocatable :: degtab(:,:)
 complex(dpc),allocatable :: sym_sigx(:,:,:)

!************************************************************************

 ! Find number of degenerates subspaces and number of bands in each subspace.
 ! The tolerance is a little bit arbitrary (0.001 eV)
 ! It could be reduced, in particular in case of nearly accidental degeneracies

 ABI_ICALLOC(degtab, (bmin:bmax, bmin:bmax))
 do ib=bmin,bmax
   do jb=bmin,bmax
    if (abs(qp_ene(ib, jk_ibz, spin) - qp_ene(jb, jk_ibz, spin)) < 0.001 / Ha_ev) degtab(ib, jb) = 1
   end do
 end do

 ABI_MALLOC(sym_sigx, (bmin:bmax, bmin:bmax, nsig_ab))
 sym_sigx = czero

 ! Average over degenerate diagonal elements.
 do ib=bmin,bmax
   ndegs=0
   do jb=bmin,bmax
     if (degtab(ib,jb)==1) then
       if (nspinor == 1) then
         sym_sigx(ib, ib, 1) = sym_sigx(ib, ib, 1) + sum(sigx(:,jb,jb,spin))
       else
         do ii=1,nsig_ab
           sym_sigx(ib, ib, ii) = sym_sigx(ib, ib, ii) + sum(sigx(:,jb,jb,ii))
         end do
       end if
     end if
     ndegs = ndegs + degtab(ib,jb)
   end do
   sym_sigx(ib,ib,:) = sym_sigx(ib,ib,:) / ndegs
 end do

 !if (gwcalctyp >= 20) call esymm_symmetrize_mels(QP_sym(spin),bmin,bmax,sigx(:,:,:,spin),sym_sigx(:,:,1))

 ! Copy symmetrized values.
 do ib=bmin,bmax
   do jb=bmin,bmax
     if (nspinor == 1) then
       sigxme_tmp(ib,jb,spin) = sym_sigx(ib,jb,1)
     else
       do ii=1,nsig_ab
         sigxme_tmp(ib,jb,ii) = sym_sigx(ib,jb,ii)
       end do
     end if
   end do
 end do

 ABI_FREE(sym_sigx)
 ABI_FREE(degtab)

end subroutine sigx_symmetrize
!!***

end module m_sigx
!!***
