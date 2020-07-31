!!****m* ABINIT/m_sigc
!! NAME
!!  m_sigc
!!
!! FUNCTION
!!  Compute matrix elements of the correlated part of the e-h self-energy
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

module m_sigc

 use defs_basis
 use m_gwdefs
 use m_abicore
 use m_xmpi
 use m_xomp
 use m_defs_ptgroups
 use m_errors
 use m_time
 use m_splines
 use m_dtset


 use defs_datatypes,  only : pseudopotential_type, ebands_t
 use m_hide_blas,     only : xdotc, xgemv, xgemm
 use m_numeric_tools, only : hermitianize, imin_loc, coeffs_gausslegint
 use m_fstrings,      only : sjoin, itoa
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_t
 use m_bz_mesh,       only : kmesh_t, get_BZ_item, findqg0, littlegroup_t, littlegroup_print
 use m_gsphere,       only : gsphere_t, gsph_fft_tabs
 use m_fft_mesh,      only : get_gftt, rotate_fft_mesh, cigfft
 use m_vcoul,         only : vcoul_t
 use m_wfd,           only : wfd_t, wave_t
 use m_oscillators,   only : rho_tw_g, calc_wfwfg
 use m_screening,     only : epsilonm1_results, epsm1_symmetrizer, epsm1_symmetrizer_inplace, get_epsm1
 use m_ppmodel,       only : setup_ppmodel, ppm_get_qbz, ppmodel_t, calc_sig_ppm
 use m_sigma,         only : sigma_t, sigma_distribute_bks
 use m_esymm,         only : esymm_t, esymm_symmetrize_mels, esymm_failed
 use m_pawang,        only : pawang_type
 use m_pawtab,        only : pawtab_type
 use m_pawfgrtab,     only : pawfgrtab_type
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy, paw_overlap
 use m_pawpwij,       only : pawpwff_t, pawpwij_t, pawpwij_init, pawpwij_free, paw_rho_tw_g, paw_cross_rho_tw_g
 use m_paw_sym,       only : paw_symcprj
 use m_paw_pwaves_lmn,only : paw_pwaves_lmn_t

 implicit none

 private
!!***

 public :: calc_sigc_me
!!***

contains
!!***

!!****f* ABINIT/calc_sigc_me
!! NAME
!! calc_sigc_me
!!
!! FUNCTION
!! Calculate diagonal and off-diagonal matrix elements of the self-energy operator.
!!
!! INPUTS
!! sigmak_ibz=Index of the k-point in the IBZ.
!! minbnd, maxbnd= min and Max band index for GW correction (for this k-point)
!! Dtset <type(dataset_type)>=all input variables in this dataset
!!    %iomode
!!    %paral_kgb
!!    %nspinor=Number of spinorial components.
!!    %gwcomp=If 1 use an extrapolar approximation to accelerate convergence.
!!    %gwencomp=Extrapolar energy.
!! Er <Epsilonm1_results> (see the definition of this structured datatype)
!!    %mqmem=if 0 use out-of-core method in which a single q-slice of espilon is read inside the loop over k
!!      much slower but it requires less memory
!!    %nomega_i=Number of imaginary frequencies.
!!    %nomega_r=Number of real frequencies.
!!    %nomega=Total number of frequencies.
!! Gsph_c<gsphere_t>= info on G-sphere for Sigma_c
!! Gsph_Max<gsphere_t>= info on biggest G-sphere
!!    %nsym=number of symmetry operations
!!    %rottb(ng,timrev,nsym)=index of (IS) G where I is the identity or the inversion
!!      operation and G is one of the ng vectors in reciprocal space
!!    %timrev=2 if time-reversal symmetry is used, 1 otherwise
!!    %gvec(3,ng)=integer coordinates of each plane wave in reciprocal space
!! ikcalc=index in the array Sigp%kptgw2bz of the k-point where GW corrections are calculated
!! Ltg_k datatype containing information on the little group
!! Kmesh <kmesh_t>
!!    %nbz=Number of points in the BZ
!!    %nibz=Number of points in IBZ
!!    %kibz(3,nibz)=k-point coordinates, irreducible Brillouin zone
!!    %kbz(3,nbz)=k-point coordinates, full Brillouin zone
!!    %ktab(nbz)= table giving for each k-point in the BZ (kBZ), the corresponding
!!    %ktabi(nbz)= for each k-point in the BZ defines whether inversion has to be considered
!!    %ktabp(nbz)= phase factor associated to tnons
!! gwc_ngfft(18)=Information about 3D FFT for the oscillator strengths used for the correlation part,
!! Vcp <vcoul_t datatype> containing information on the cutoff technique
!!    %vc_sqrt(npwc,nqibz)= square-root of the coulombian potential for q-points in the IBZ
!! Pawtab(Psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Pawang <type(pawang_type)>=paw angular mesh and related data
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!    %usepaw=1 for PAW, 0 for NC pseudopotentials.
!! Qmesh <kmesh_t> : datatype gathering information of the q-mesh used
!!    %ibz=q points where $\tilde\epsilon^{-1}$ has been computed
!!    %bz(3,nqbz)=coordinates of all q-points in BZ
!! Sigp <sigparams_t> (see the definition of this structured datatype)
!! Cryst<crystal_t>=Info on unit cell and symmetries
!!    %natom=number of atoms in unit cell
!!    %ucvol=unit cell volume
!!    %nsym=number of symmetry operations
!!    %typat(natom)=type of each atom
!! PPm<ppmodel_t>= Datatype gathering information on the Plasmonpole technique (see also ppm_get_qbz).
!!    %model=type of ppmodel
!!    %npwc=number of G-vectors for the correlation part.
!!    %dm2_otq =size of second dimension of otq array
!!    %dm2_bots=size of second dimension of botsq arrays
!!    %dm_eig  =size of second dimension of eig arrays
!! QP_BSt<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!!  Paw_pwff<pawpwff_t>=Form factor used to calculate the onsite mat. elements of a plane wave.
!! allQP_sym(Wfd%nkibz,Wfd%nsppol)<esymm_t>=Datatype collecting data on the irreducible representaions of the
!!    little group of kcalc in the KS representation as well as the symmetry of the bdgw_k states.
!!  Sr=sigma_t (see the definition of this structured datatype)
!!  use_aerhor=1 is aepaw_rhor is used, 0 otherwise.
!!  aepaw_rhor(rho_nfftot,Wfd%nspden*use_aerhor)=AE PAW density used to generate PPmodel paramenters if mqmem==0
!!
!! OUTPUT
!!
!! NOTES
!!  1) The treatment of the divergence of Gygi+Baldereschi (PRB 1986) [[cite:Gigy1986]] is included.
!!  2) The calculation of energy derivative is based on finite elements.
!!  3) On the symmetrization of Sigma matrix elements ***/
!!        If  Sk = k+G0 then  M_G(k, Sq)= e^{-i( Sq+G).t} M_{ S^-1(G}   (k,q)
!!        If -Sk = k+G0 then  M_G(k,-Sq)= e^{-i(-Sq+G).t} M_{-S^-1(G)}^*(k,q)
!!
!!     Notice the absence of G0 in the expression. Moreover, when we sum over the little group, it turns out
!!     that there is a cancellation of the phase factor associated to the non-symmorphic operations due to a
!!     similar term coming from the symmetrization of \epsilon^{-1}. Mind however that the nonsymmorphic phase
!!     has to be considered when epsilon^-1 is reconstructed starting from the q-points in the IBZ.
!!
!!  4) The unitary transformation relating wavefunctions
!!     at symmetric k-points should be taken into account during the symmetrization
!!     of the oscillator matrix elements. In case of G_oW_o and GW_o calculations, however,
!!     it is possible to make an invariant by just including all the degenerate states and
!!     averaging the final results over the degenerate subset.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      calc_coh_comp,calc_sig_ppm,calc_sig_ppm_comp,calc_sigc_cd,calc_wfwfg
!!      coeffs_gausslegint,cwtime,epsm1_symmetrizer,epsm1_symmetrizer_inplace
!!      esymm_symmetrize_mels,findqg0,get_bz_item,get_epsm1,get_gftt
!!      gsph_fft_tabs,littlegroup_print,paw_cross_rho_tw_g,paw_rho_tw_g
!!      paw_symcprj,pawcprj_alloc,pawcprj_copy,pawcprj_free,pawpwij_free
!!      pawpwij_init,ppm_get_qbz,rho_tw_g,rotate_fft_mesh,setup_ppmodel
!!      sigma_distribute_bks,timab,wfd_change_ngfft,wfd_get_cprj
!!      wfd_get_many_ur,wfd_get_ur,wfd_paw_get_aeur,wrtout,xmpi_max,xmpi_sum
!!
!! SOURCE

subroutine calc_sigc_me(sigmak_ibz,ikcalc,nomega_sigc,minbnd,maxbnd,&
& Dtset,Cryst,QP_BSt,Sigp,Sr,Er,Gsph_Max,Gsph_c,Vcp,Kmesh,Qmesh,Ltg_k,&
& PPm,Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,allQP_sym,&
& gwc_ngfft,rho_ngfft,rho_nfftot,rhor,use_aerhor,aepaw_rhor,sigcme_tmp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: sigmak_ibz,ikcalc,rho_nfftot,nomega_sigc,minbnd,maxbnd
 integer,intent(in) :: use_aerhor
 type(crystal_t),intent(in) :: Cryst
 type(ebands_t),target,intent(in) :: QP_BSt
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(vcoul_t),intent(in) :: Vcp
 type(dataset_type),intent(in) :: Dtset
 type(epsilonm1_results),intent(inout) :: Er
 type(gsphere_t),intent(in) :: Gsph_Max,Gsph_c
 type(littlegroup_t),intent(in) :: Ltg_k
 type(ppmodel_t),intent(inout) :: PPm
 type(Pseudopotential_type),intent(in) :: Psps
 type(pawang_type),intent(in) :: pawang
 type(sigparams_t),target,intent(in) :: Sigp
 type(sigma_t),intent(in) :: Sr
 type(wfd_t),target,intent(inout) :: Wfd,Wfdf
!arrays
 integer,intent(in) :: gwc_ngfft(18),rho_ngfft(18)
 real(dp),intent(in) :: rhor(rho_nfftot,Wfd%nspden)
 real(dp),intent(in) :: aepaw_rhor(rho_nfftot,Wfd%nspden*use_aerhor)
 complex(dpc),intent(out) :: sigcme_tmp(nomega_sigc,minbnd:maxbnd,minbnd:maxbnd,Wfd%nsppol*Sigp%nsig_ab)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
 type(esymm_t),target,intent(in) :: allQP_sym(Wfd%nkibz,Wfd%nsppol)
 type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
 type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom*Psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp2=2,ndat1=1
 integer :: npw_k,iab,ib,ib1,ib2,ierr,ig,iggp,igp,ii,iik,itim_q,i1,i2,npls
 integer :: ik_bz,ik_ibz,io,iiw,isym_q,iq_bz,iq_ibz,spin,isym,jb,is_idx
 integer :: band,band1,band2,idle,rank,jik,jk_bz,jk_ibz,kb,nspinor
 integer :: nomega_tot,nq_summed,ispinor,ibsp,dimcprj_gw,npwc
 integer :: spad,spadc,spadc1,spadc2,irow,my_nbks,ndegs,wtqm,wtqp,mod10
 integer :: isym_kgw,isym_ki,gwc_mgfft,use_padfft,gwc_fftalga,gwc_nfftot,nfftf,mgfftf,use_padfftf
 integer :: iwc,ifft
 real(dp) :: cpu_time,wall_time,gflops
 real(dp) :: e0i,fact_sp,theta_mu_minus_e0i,tol_empty,z2,en_high,gw_gsq,w_localmax,w_max
 complex(dpc) :: ctmp,omegame0i2_ac,omegame0i_ac,ph_mkgwt,ph_mkt
 logical :: iscompatibleFFT,q_is_gamma
 character(len=500) :: msg,sigma_type
 type(wave_t),pointer :: wave_sum, wave_jb
 complex(gwpc),allocatable :: botsq(:,:),otq(:,:),eig(:,:)
!arrays
 integer :: g0(3),spinor_padc(2,4),got(Wfd%nproc)
 integer,allocatable :: proc_distrb(:,:,:),extrapolar_distrb(:,:,:,:),degtab(:,:,:)
 integer,allocatable :: igfftcg0(:),gw_gfft(:,:),gw_gbound(:,:),irottb(:,:),ktabr(:,:)
 integer,allocatable :: igfftfcg0(:),gboundf(:,:),ktabrf(:,:),npoles_missing(:)
 real(dp) :: ksum(3),kgw(3),kgw_m_ksum(3),omegap(Er%nomega_i),omegap2(Er%nomega_i),q0(3),tsec(2),qbz(3)
 real(dp) :: gl_knots(Er%nomega_i),gl_wts(Er%nomega_i)
 real(dp) :: spinrot_kbz(4),spinrot_kgw(4)
 real(dp),ABI_CONTIGUOUS pointer :: qp_ene(:,:,:),qp_occ(:,:,:)
 real(dp),allocatable :: omegame0i(:), w_maxval(:)
 complex(gwpc) :: sigcohme(Sigp%nsig_ab)
 complex(gwpc),allocatable :: vc_sqrt_qbz(:),rhotwg(:),rhotwgp(:)
 complex(gwpc),allocatable :: botsq_conjg_transp(:,:),ac_epsm1cqwz2(:,:,:)
 complex(gwpc),allocatable :: epsm1_qbz(:,:,:),epsm1_trcc_qbz(:,:,:), epsm1_tmp(:,:)
 complex(gwpc),allocatable :: sigc_ket(:,:),ket1(:,:),ket2(:,:)
 complex(gwpc),allocatable :: herm_sigc_ket(:,:),aherm_sigc_ket(:,:)
 complex(gwpc),allocatable :: rhotwg_ki(:,:)
 complex(gwpc),allocatable :: sigcme2(:,:),sigcme_3(:),sigcme_new(:),sigctmp(:,:)
 complex(gwpc),allocatable :: wfr_bdgw(:,:),ur_ibz(:),wf1swf2_g(:),usr_bz(:)
 complex(gwpc),allocatable :: ur_ae_sum(:),ur_ae_onsite_sum(:),ur_ps_onsite_sum(:)
 complex(gwpc),allocatable :: ur_ae_bdgw(:,:),ur_ae_onsite_bdgw(:,:),ur_ps_onsite_bdgw(:,:)
 complex(gwpc),allocatable :: otq_transp(:,:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: cg_jb(:),cg_sum(:)
 complex(dpc),allocatable :: sym_cme(:,:,:,:),sigc(:,:,:,:,:)
 logical :: rank_mask(Wfd%nproc),can_symmetrize(Wfd%nsppol)
!logical :: me_calc_poles(Sr%nomega_r+Sr%nomega4sd)
 type(sigijtab_t),pointer :: Sigcij_tab(:)
 type(pawcprj_type),allocatable :: Cprj_kgw(:,:),Cprj_ksum(:,:)
 type(pawpwij_t),allocatable :: Pwij_qg(:),Pwij_fft(:)
 type(esymm_t),pointer :: QP_sym(:)
 integer :: ilwrk
 integer :: neig(Er%nomega_i)
 real(gwp) :: epsm1_ev(Sigp%npwc)
 complex(gwpc),allocatable :: epsm1_sqrt_rhotw(:,:)
 complex(gwpc),allocatable :: rhotw_epsm1_rhotw(:,:,:)

!************************************************************************

 DBG_ENTER("COLL")

 ! Initial check
 ABI_CHECK(Sr%nomega_r==Sigp%nomegasr,"")
 ABI_CHECK(Sr%nomega4sd==Sigp%nomegasrd,"")
 ABI_CHECK(Sigp%npwc==Gsph_c%ng,"")
 ABI_CHECK(Sigp%npwvec==Gsph_Max%ng,"")

 mod10=MOD(Sigp%gwcalctyp,10)

 call timab(424,1,tsec) ! calc_sigc_me
 call timab(431,1,tsec) ! calc_sigc_me
 call timab(432,1,tsec) ! Init
 call cwtime(cpu_time,wall_time,gflops,"start")

 ! Initialize MPI variables
 qp_ene => QP_BSt%eig; qp_occ => QP_BSt%occ

 ! Extract the symmetries of the bands for this k-point
 QP_sym => allQP_sym(sigmak_ibz,1:Wfd%nsppol)

 ! Index of the GW point in the BZ array, its image in IBZ and time-reversal ===
 jk_bz=Sigp%kptgw2bz(ikcalc)
 call get_BZ_item(Kmesh,jk_bz,kgw,jk_ibz,isym_kgw,jik,ph_mkgwt)
 !%call get_IBZ_item(Kmesh,jk_ibz,kibz,wtk)

 ! TODO: the new version based of get_uug won't suppporte kptgw vectors that are not in
 ! the IBZ since one should perform the rotation before entering the band loop
 ! In the old version, the rotation was done in rho_tw_g
 !ABI_CHECK(jik==1,"jik!=1")
 !ABI_CHECK(isym_kgw==1,"isym_kgw!=1")
 !ABI_CHECK((ABS(ph_mkgwt - cone) < tol12),"ph_mkgwt!")

 spinrot_kgw=Cryst%spinrot(:,isym_kgw)
 ib1=minbnd; ib2=maxbnd

 write(msg,'(2a,3f8.3,2a,2(i3,a))')ch10,&
&  ' Calculating <nk|Sigma_c(omega)|nk> at k = ',kgw(:),ch10,&
&  ' bands n = from ',ib1,' to ',ib2,ch10
 call wrtout(std_out,msg,'COLL')

 if ( Dtset%gwaclowrank > 0 ) then
   ! Today we use the same number of eigenvectors irrespective to iw'
   ! Tomorrow we might optimize this further
   neig(:) = MIN(Dtset%gwaclowrank,Sigp%npwc)
 else
   neig(:) = Sigp%npwc
 endif 

 ABI_MALLOC(w_maxval,(minbnd:maxbnd))
 w_maxval = zero

 if ( ANY(gwc_ngfft(1:3) /= Wfd%ngfft(1:3)) ) call Wfd%change_ngfft(Cryst,Psps,gwc_ngfft)
 gwc_mgfft   = MAXVAL(gwc_ngfft(1:3))
 gwc_fftalga = gwc_ngfft(7)/100 !; gwc_fftalgc=MOD(gwc_ngfft(7),10)

 if (Dtset%pawcross==1) mgfftf = MAXVAL(rho_ngfft(1:3))

 can_symmetrize = .FALSE.
 if (Sigp%symsigma>0) then
   can_symmetrize = .TRUE.
   if (Sigp%gwcalctyp >= 20) then
    do spin=1,Wfd%nsppol
      can_symmetrize(spin) = .not. esymm_failed(QP_sym(spin))
      if (.not.can_symmetrize(spin)) then
        write(msg,'(a,i0,4a)')&
&         " Symmetrization cannot be performed for spin: ",spin,ch10,&
&         " band classification encountered the following problem: ",ch10,TRIM(QP_sym(spin)%err_msg)
        MSG_WARNING(msg)
      end if
    end do
   end if
   if (Wfd%nspinor == 2) MSG_WARNING("Symmetrization with nspinor=2 not implemented")
 end if

 ABI_UNUSED(Pawang%l_max)

 ! Print type of calculation.
 sigma_type = sigma_type_from_key(mod10)
 call wrtout(std_out,sigma_type,'COLL')

 ! Set up logical flags for Sigma calculation
 if (mod10==SIG_GW_AC.and.Sigp%gwcalctyp/=1) then
   MSG_ERROR("not implemented")
 end if

 if (mod10==SIG_GW_AC.and.Sigp%gwcomp==1) then
   MSG_ERROR("not implemented")
 end if

 if (mod10==SIG_GW_AC) then
   write(msg,'(3a,i6,a,i6)') ' Using a low-rank formula for AC',&
&           ch10,' Number of epsm1 eigenvectors retained: ',neig(1),' over ',Sigp%npwc
   call wrtout(std_out,msg,'COLL')
 endif


 ! Initialize some values
 nspinor = Wfd%nspinor
 npwc = Sigp%npwc
 spinor_padc(:,:)=RESHAPE([0, 0, npwc, npwc, 0, npwc, npwc, 0], [2, 4])
 ABI_MALLOC(npoles_missing,(minbnd:maxbnd))
 npoles_missing=0

 ! Normalization of theta_mu_minus_e0i
 ! If nsppol==2, qp_occ $\in [0,1]$
 SELECT CASE (Wfd%nsppol)
 CASE (1)
   fact_sp=half; tol_empty=0.01   ! below this value the state is assumed empty
   if (Wfd%nspinor==2) then
    fact_sp=one; tol_empty=0.005  ! below this value the state is assumed empty
   end if
 CASE (2)
   fact_sp=one; tol_empty=0.005   ! to be consistent and obtain similar results if a metallic
 CASE DEFAULT                     ! spin unpolarized system is treated using nsppol==2
   MSG_BUG('Wrong nsppol')
 END SELECT

 ! Allocate arrays used to accumulate the matrix elements of \Sigma_c over
 ! k-points and bands. Note that for AC requires only the imaginary frequencies
 !nomega_sigc=Sr%nomega_r+Sr%nomega4sd
 !if (mod10==SIG_GW_AC) nomega_sigc=Sr%nomega_i
 !
 ! === Define the G-G0 shifts for the FFT of the oscillators ===
 ! * Sigp%mG0 gives the MAX G0 component to account for umklapp.
 ! * Note the size MAX(Sigp%npwx,npwc).
 !
 ! === Precalculate the FFT index of $(R^{-1}(r-\tau))$ ===
 ! * S=\transpose R^{-1} and k_BZ = S k_IBZ
 ! * irottb is the FFT index of $R^{-1} (r-\tau)$ used to symmetrize u_Sk.
 gwc_nfftot = PRODUCT(gwc_ngfft(1:3))
 ABI_MALLOC(irottb,(gwc_nfftot,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,gwc_ngfft,irottb,iscompatibleFFT)
 if (.not.iscompatibleFFT) then
   msg = "FFT mesh is not compatible with symmetries. Results might be affected by large errors!"
   MSG_WARNING(msg)
 end if

 ABI_MALLOC(ktabr,(gwc_nfftot,Kmesh%nbz))
 do ik_bz=1,Kmesh%nbz
   isym=Kmesh%tabo(ik_bz)
   do ifft=1,gwc_nfftot
     ktabr(ifft,ik_bz)=irottb(ifft,isym)
   end do
 end do
 ABI_FREE(irottb)

 if (Psps%usepaw==1 .and. Dtset%pawcross==1) then
   nfftf = PRODUCT(rho_ngfft(1:3))
   ABI_MALLOC(irottb,(nfftf,Cryst%nsym))
   call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,rho_ngfft,irottb,iscompatibleFFT)

   ABI_MALLOC(ktabrf,(nfftf,Kmesh%nbz))
   do ik_bz=1,Kmesh%nbz
     isym=Kmesh%tabo(ik_bz)
     do ifft=1,nfftf
       ktabrf(ifft,ik_bz)=irottb(ifft,isym)
     end do
   end do
   ABI_FREE(irottb)
 end if

 Sigcij_tab => Sigp%Sigcij_tab(ikcalc,1:Wfd%nsppol)

 got=0
 ABI_MALLOC(proc_distrb,(Wfd%mband,Kmesh%nbz,Wfd%nsppol))
 call sigma_distribute_bks(Wfd,Kmesh,Ltg_k,Qmesh,Wfd%nsppol,can_symmetrize,kgw,Sigp%mg0,&
&  my_nbks,proc_distrb,got,global=.TRUE.)

 write(msg,'(a,i0,a)')" Will sum ",my_nbks," (b,k,s) states in Sigma_c."
 call wrtout(std_out,msg,'PERS')

 if (Sigp%gwcomp==1) then
   en_high=MAXVAL(qp_ene(Sigp%nbnds,:,:)) + Sigp%gwencomp
   write(msg,'(6a,e11.4,a)')ch10,&
&    ' Using the extrapolar approximation to accelerate convergence',ch10,&
&    ' with respect to the number of bands included',ch10,&
&    ' with extrapolar energy: ',en_high*Ha_eV,' [eV]'
   call wrtout(std_out,msg,'COLL')
   ABI_MALLOC(wf1swf2_g,(gwc_nfftot*nspinor))
 end if

 if (Sigp%gwcomp == 1) then
   ! Setup of MPI table for extrapolar contributions.
   ABI_MALLOC(extrapolar_distrb,(ib1:ib2,ib1:ib2,Kmesh%nbz,Wfd%nsppol))
   extrapolar_distrb = xmpi_undefined_rank

   do spin=1,Wfd%nsppol
     do ik_bz=1,Kmesh%nbz
        if (ANY(proc_distrb(:,ik_bz,spin) /= xmpi_undefined_rank) ) then ! This BZ point will be calculated.
           rank_mask = .FALSE. ! The set of node that will treat (k,s).
           do band=1,Wfd%mband
             rank = proc_distrb(band,ik_bz,spin)
             if (rank /= xmpi_undefined_rank) rank_mask(rank+1)=.TRUE.
           end do
           do band2=ib1,ib2
             do irow=1,Sigcij_tab(spin)%col(band2)%size1   ! Looping over the non-zero elements of sigma_ij.
               band1 = Sigcij_tab(spin)%col(band2)%bidx(irow)
               idle = imin_loc(got,mask=rank_mask)
               got(idle) = got(idle)+1
               extrapolar_distrb(band1,band2,ik_bz,spin) = idle-1
             end do
           end do
        end if
     end do
   end do

   write(msg,'(a,i0,a)')" Will treat ",COUNT(extrapolar_distrb==Wfd%my_rank)," extrapolar terms."
   call wrtout(std_out,msg,'PERS')
 end if

 ABI_MALLOC(rhotwg_ki, (npwc*nspinor, minbnd:maxbnd))
 rhotwg_ki=czero_gw
 ABI_MALLOC(rhotwg, (npwc*nspinor))
 ABI_MALLOC(rhotwgp, (npwc*nspinor))
 ABI_MALLOC(vc_sqrt_qbz, (npwc))

 if (Er%mqmem==0) then ! Use out-of-core solution for epsilon.
   MSG_COMMENT('Reading q-slices from file. Slower but less memory.')
 end if                                                                                !

 ! Additional allocations for PAW.
 if (Psps%usepaw==1) then
   ABI_MALLOC(Cprj_ksum,(Cryst%natom,nspinor))
   call pawcprj_alloc(Cprj_ksum,0,Wfd%nlmn_atm)
   !
   ! For the extrapolar method we need the onsite terms of the PW in the FT mesh.
   ! * gw_gfft is the set of plane waves in the FFT Box for the oscillators.
   if (Sigp%gwcomp==1) then
     ABI_MALLOC(gw_gfft,(3,gwc_nfftot))
     q0=zero
     call get_gftt(gwc_ngfft,q0,Cryst%gmet,gw_gsq,gw_gfft)
     ABI_MALLOC(Pwij_fft,(Psps%ntypat))
     call pawpwij_init(Pwij_fft,gwc_nfftot,(/zero,zero,zero/),gw_gfft,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
   end if
 end if ! usepaw==1

 if (mod10==SIG_GW_AC) then ! Calculate Gauss-Legendre quadrature knots and weights for analytic continuation

   ABI_MALLOC(rhotw_epsm1_rhotw, (minbnd:maxbnd, minbnd:maxbnd, Er%nomega_i))

   call coeffs_gausslegint(zero,one,gl_knots,gl_wts,Er%nomega_i)

   do io=1,Er%nomega_i ! First frequencies are always real
     if (ABS(AIMAG(one*Er%omega(Er%nomega_r+io))-(one/gl_knots(io)-one)) > 0.0001) then
      write(msg,'(3a)')&
&       ' Frequencies in the SCR file are not compatible with the analytic continuation. ',ch10,&
&       ' Verify the frequencies in the SCR file. '
      MSG_WARNING(msg)
      if (Wfd%my_rank==Wfd%master) write(std_out,*)AIMAG(Er%omega(Er%nomega_r+io)),(one/gl_knots(io)-one)
      MSG_ERROR("Cannot continue!")
     end if
   end do

   ! To calculate \int_0^\infty domegap f(omegap), we calculate \int_0^1 dz f(1/z-1)/z^2.
   omegap(:)=one/gl_knots(:)-one
   omegap2(:)=omegap(:)*omegap(:)
   ABI_MALLOC(ac_epsm1cqwz2, (npwc, npwc, Er%nomega_i))
 end if

 ! Calculate total number of frequencies and allocate related arrays.
 ! sigcme2 is used to accumulate the diagonal matrix elements over k-points and
 ! GW bands, used only in case of ppmodel 3 and 4 (TODO save memory)
 nomega_tot=Sr%nomega_r+Sr%nomega4sd
 ABI_MALLOC(sigcme2,(nomega_tot,ib1:ib2))
 ABI_MALLOC(sigcme_3,(nomega_tot))
 sigcme2=czero_gw; sigcme_3=czero_gw

 ABI_MALLOC(sigctmp,(nomega_sigc,Sigp%nsig_ab))
 sigctmp=czero_gw
 if (mod10/=SIG_GW_AC) then
   ABI_MALLOC(sigc_ket, (npwc*nspinor, nomega_sigc))
 endif

#if 0
 !TODO gmatteo: these arrays are never used in practice. Should we remove them?
 ! Arrays storing the contribution given by the Hermitian/anti-Hermitian part of \Sigma_c
 ABI_MALLOC(aherm_sigc_ket, (npwc*nspinor, nomega_sigc))
 ABI_MALLOC( herm_sigc_ket, (npwc*nspinor, nomega_sigc))
#endif

 sigcme_tmp=czero

 ABI_MALLOC(sigc,(2,nomega_sigc,ib1:ib2,ib1:ib2,Wfd%nsppol*Sigp%nsig_ab))
 sigc=czero

 if( mod10==SIG_QPGW_PPM .or. mod10==SIG_QPGW_CD ) then
   ABI_MALLOC(ket1, (npwc*nspinor, nomega_tot))
   ABI_MALLOC(ket2, (npwc*nspinor, nomega_tot))
 endif
 ABI_MALLOC(omegame0i,(nomega_tot))

 ! Here we divide the states where the QP energies are required into complexes. Note however that this approach is not
 ! based on group theory, and it might lead to spurious results in case of accidental degeneracies.
 nq_summed=Kmesh%nbz
 if (Sigp%symsigma>0) then
   call littlegroup_print(Ltg_k,std_out,Dtset%prtvol,'COLL')
   nq_summed=SUM(Ltg_k%ibzq(:))
   !
   ! Find number of degenerate subspaces and number of bands in each subspace
   ! The tolerance is a little bit arbitrary (0.001 eV)
   ! It could be reduced, in particular in case of nearly accidental degeneracies
   ABI_MALLOC(degtab,(ib1:ib2,ib1:ib2,Wfd%nsppol))
   degtab=0
   do spin=1,Wfd%nsppol
     do ib=ib1,ib2
       do jb=ib1,ib2
        if (ABS(qp_ene(ib,jk_ibz,spin)-qp_ene(jb,jk_ibz,spin))<0.001/Ha_ev) then
          degtab(ib,jb,spin)=1
        end if
       end do
     end do
   end do
 end if !symsigma

 write(msg,'(2a,i6,a)')ch10,' calculation status ( ',nq_summed,' to be completed):'
 call wrtout(std_out,msg,'COLL')

 ! Here we have a problem in case of CD since epsm1q might be huge
 ! TODO if single q (ex molecule) dont allocate epsm1q, avoid waste of memory
 if ( ANY(mod10== [SIG_GW_AC, SIG_GW_CD, SIG_QPGW_CD])) then
   if (.not.(mod10==SIG_GW_CD.and.Er%mqmem==0)) then
     ABI_MALLOC_OR_DIE(epsm1_qbz, (npwc, npwc, Er%nomega), ierr)
   end if
 end if

 ! TODO epsm1_trcc_qbz is needed for SIG_GW_CD with symmetries since
 ! the Hermitian and the anti-Hermitian part have to be symmetrized in a different way.
 !if (mod10==SIG_QPGW_CD) then
 if (mod10==SIG_QPGW_CD) then
   ABI_MALLOC_OR_DIE(epsm1_trcc_qbz, (npwc, npwc, Er%nomega), ierr)
   ABI_MALLOC(epsm1_tmp, (npwc, npwc))
 end if

 ABI_MALLOC(igfftcg0,(Gsph_Max%ng))
 ABI_MALLOC(ur_ibz,(gwc_nfftot*nspinor))
 ABI_MALLOC(usr_bz,(gwc_nfftot*nspinor))

 if (Dtset%pawcross==1) then
   ABI_MALLOC(igfftfcg0,(Gsph_c%ng))
   ABI_MALLOC(ur_ae_sum,(nfftf*nspinor))
   ABI_MALLOC(ur_ae_onsite_sum,(nfftf*nspinor))
   ABI_MALLOC(ur_ps_onsite_sum,(nfftf*nspinor))
 end if
 call timab(432,2,tsec) ! Init

 ! ==========================================
 ! ==== Fat loop over k_i in the full BZ ====
 ! ==========================================
 do spin=1,Wfd%nsppol
   if (ALL(proc_distrb(:,:,spin)/=Wfd%my_rank)) CYCLE
   call timab(433,1,tsec) ! Init spin

   ! Load wavefunctions for GW corrections
   ! TODO: Rotate the functions here instead of calling rho_tw_g
   ABI_MALLOC(wfr_bdgw,(gwc_nfftot*nspinor,ib1:ib2))
   call Wfd%get_many_ur([(jb, jb=ib1,ib2)], jk_ibz, spin, wfr_bdgw)

   if (Wfd%usepaw==1) then
     ! Load cprj for GW states, note the indexing.
     dimcprj_gw=nspinor*(ib2-ib1+1)
     ABI_MALLOC(Cprj_kgw,(Cryst%natom,ib1:ib1+dimcprj_gw-1))
     call pawcprj_alloc(Cprj_kgw,0,Wfd%nlmn_atm)
     ibsp=ib1
     do jb=ib1,ib2
       call Wfd%get_cprj(jb,jk_ibz,spin,Cryst,Cprj_ksum,sorted=.FALSE.)
       call paw_symcprj(jk_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,Cprj_ksum)
       call pawcprj_copy(Cprj_ksum,Cprj_kgw(:,ibsp:ibsp+(nspinor-1)))
       ibsp=ibsp+nspinor
     end do
     if (Dtset%pawcross==1) then
       ABI_MALLOC(ur_ae_bdgw,(nfftf*nspinor,ib1:ib2))
       ABI_MALLOC(ur_ae_onsite_bdgw,(nfftf*nspinor,ib1:ib2))
       ABI_MALLOC(ur_ps_onsite_bdgw,(nfftf*nspinor,ib1:ib2))
       do jb=ib1,ib2
         call Wfdf%paw_get_aeur(jb,jk_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,&
&          ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum)
         ur_ae_bdgw(:,jb)=ur_ae_sum
         ur_ae_onsite_bdgw(:,jb)=ur_ae_onsite_sum
         ur_ps_onsite_bdgw(:,jb)=ur_ps_onsite_sum
       end do
     end if
   end if ! usepaw

   call timab(433,2,tsec) ! Init spin

   do ik_bz=1,Kmesh%nbz
     ! Parallelization over k-points and spin
     ! For the spin there is another check in the inner loop
     if (ALL(proc_distrb(:,ik_bz,spin)/=Wfd%my_rank)) CYCLE

     call timab(434,1,tsec) ! initq

     ! Find the corresponding irreducible k-point
     call get_BZ_item(Kmesh,ik_bz,ksum,ik_ibz,isym_ki,iik,ph_mkt)
     spinrot_kbz(:)=Cryst%spinrot(:,isym_ki)

     ! Identify q and G0 where q + G0 = k_GW - k_i
     kgw_m_ksum=kgw-ksum
     call findqg0(iq_bz,g0,kgw_m_ksum,Qmesh%nbz,Qmesh%bz,Sigp%mG0)

     ! If symsigma, symmetrize the matrix elements.
     ! Sum only q"s in IBZ_k. In this case elements are weighted
     ! according to wtqp and wtqm. wtqm is for time-reversal.
     wtqp=1; wtqm=0
     if (can_symmetrize(spin)) then
       if (Ltg_k%ibzq(iq_bz)/=1) CYCLE
       wtqp=0; wtqm=0
       do isym=1,Ltg_k%nsym_sg
         wtqp=wtqp+Ltg_k%wtksym(1,isym,iq_bz)
         wtqm=wtqm+Ltg_k%wtksym(2,isym,iq_bz)
       end do
     end if

     write(msg,'(3(a,i0),a,i0)')' Sigma_c: ik_bz ',ik_bz,'/',Kmesh%nbz,", spin: ",spin,' done by mpi-rank: ',Wfd%my_rank
     call wrtout(std_out,msg,'PERS')

     ! Find the corresponding irred q-point.
     call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)
     q_is_gamma = (normv(qbz,Cryst%gmet,"G") < GW_TOL_W0)

     !q_is_gamma = (normv(qbz,Cryst%gmet,"G") < 0.7)
     !if (iq_ibz/=2.and.iq_ibz/=1) CYCLE
     !if (ANY(qbz<=-(half-tol16)) .or. ANY(qbz>(half+tol16))) CYCLE
     !if (q_is_gamma) then; write(std_out,*)"skipping q=Gamma"; CYCLE; end if
     !
     ! Tables for the FFT of the oscillators.
     !  a) FFT index of the G-G0.
     !  b) gw_gbound table for the zero-padded FFT performed in rhotwg.
     ABI_MALLOC(gw_gbound,(2*gwc_mgfft+8,2))
     call gsph_fft_tabs(Gsph_c,g0,gwc_mgfft,gwc_ngfft,use_padfft,gw_gbound,igfftcg0)

     if (ANY(gwc_fftalga == [2, 4])) use_padfft=0 ! Pad-FFT is not coded in rho_tw_g
#ifdef FC_IBM
     ! XLF does not deserve this optimization (problem with [v67mbpt][t03])
     use_padfft = 0
#endif
     if (use_padfft==0) then
       ABI_FREE(gw_gbound)
       ABI_MALLOC(gw_gbound,(2*gwc_mgfft+8,2*use_padfft))
     end if

     if (Dtset%pawcross==1) then
       ABI_MALLOC(gboundf,(2*mgfftf+8,2))
       call gsph_fft_tabs(Gsph_c,g0,mgfftf,rho_ngfft,use_padfftf,gboundf,igfftfcg0)
       if ( ANY(gwc_fftalga == (/2,4/)) ) use_padfftf=0
       if (use_padfftf==0) then
         ABI_FREE(gboundf)
         ABI_MALLOC(gboundf,(2*mgfftf+8,2*use_padfftf))
       end if
     end if

     ! Evaluate oscillator matrix elements
     ! $ <phj/r|e^{-i(q+G)}|phi/r> - <tphj/r|e^{-i(q+G)}|tphi/r> $ in packed form.
     if (Psps%usepaw==1) then
       ABI_MALLOC(Pwij_qg,(Psps%ntypat))
       q0 = qbz !;if (q_is_gamma) q0 = (/0.00001_dp,0.00001_dp,0.00001_dp/) ! GW_Q0_DEFAULT
       call pawpwij_init(Pwij_qg,npwc,q0,Gsph_c%gvec,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
     end if

     if (Er%mqmem==0) then
       ! Read q-slice of epsilon^{-1}|chi0 in Er%epsm1(:,:,:,1) (much slower but less memory).
       call get_epsm1(Er,Vcp,0,0,Dtset%iomode,xmpi_comm_self,iqibzA=iq_ibz)
       if (sigma_needs_ppm(Sigp)) then
         if (Wfd%usepaw==1.and.PPm%userho==1) then
           ! Use PAW AE rhor.
           call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,rho_nfftot,Gsph_c%gvec,&
&            rho_ngfft,aepaw_rhor(:,1),iqiA=iq_ibz)
         else
           call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,rho_nfftot,Gsph_c%gvec,&
&            rho_ngfft,rhor(:,1),iqiA=iq_ibz)
         end if
       end if
     end if

     ! Symmetrize PPM parameters and epsm1 (q_IBZ --> q_BZ):
     ! NOTE:
     !    - We are not considering umklapp with G0/=0. In this case,
     !      indeed the equation is different since we have to use G-G0.
     !      A check, however, is performed in sigma.
     !    - If gwcomp==1 and mod10 in [1,2,9], one needs both to set up botsq and epsm1_q
     if (sigma_needs_ppm(Sigp)) then
       ABI_MALLOC(botsq,(PPm%npwc,PPm%dm2_botsq))
       ABI_MALLOC(otq,(PPm%npwc,PPm%dm2_otq))
       ABI_MALLOC(eig,(PPm%dm_eig,PPm%dm_eig))
       call ppm_get_qbz(PPm,Gsph_c,Qmesh,iq_bz,botsq,otq,eig)
     end if

     if ( ANY(mod10 == [SIG_GW_AC,SIG_GW_CD,SIG_QPGW_CD])) then
       ! Numerical integration or model GW with contour deformation or Analytic Continuation
       ! TODO In case of AC we should symmetrize only the imaginary frequencies
       if (mod10==SIG_GW_CD.and.Er%mqmem==0) then
         ! Do in-place symmetrisation
         call Epsm1_symmetrizer_inplace(iq_bz,Er%nomega, npwc,Er,Gsph_c,Qmesh,.TRUE.)
       else
         call Epsm1_symmetrizer(iq_bz,Er%nomega, npwc,Er,Gsph_c,Qmesh,.TRUE.,epsm1_qbz)
       end if

       if (mod10==SIG_GW_AC) then

         call timab(444,1,tsec) ! ac_lrk_diag
         ! Important to set to zero here for all procs since we're going to use a dirty reduction later 
         ! The reduction 'xmpi_sum' does not induce a significant performance loss in the tested systems
         ac_epsm1cqwz2(:,:,:) = czero_gw
         do iiw=1,Er%nomega_i
           ! Use the available MPI tasks to parallelize over iw'
           if ( Dtset%gwpara == 2 .and. MODULO(iiw-1,Wfd%nproc) /= Wfd%my_rank ) CYCLE
         
           ! Prepare the integration weights w_i 1/z_i^2 f(1/z_i-1)..
           ! The first frequencies are always real, skip them.
           z2=gl_knots(iiw)*gl_knots(iiw)
           ac_epsm1cqwz2(:,:,iiw)= gl_wts(iiw)*epsm1_qbz(:,:,Er%nomega_r+iiw)/z2

           ! (epsm1-1) has negative eigenvalues
           ! after the diago, they will be sorted starting from the most negative
           call xheev('V','L',npwc,ac_epsm1cqwz2(:,:,iiw),epsm1_ev)
           neig(iiw) = neig(1)
           do ilwrk=1,neig(iiw)
             ac_epsm1cqwz2(:,ilwrk,iiw) = ac_epsm1cqwz2(:,ilwrk,iiw) * SQRT( -epsm1_ev(ilwrk) )
           end do
         end do
         if ( Dtset%gwpara == 2 ) then
           call xmpi_sum(ac_epsm1cqwz2, Wfd%comm, ierr)
         endif
         call timab(444,2,tsec) ! ac_lrk_diag

       end if

       if (mod10==SIG_QPGW_CD) then
         ! For model GW we need transpose(conjg(epsm1_qbz)) ===
         do io=1,Er%nomega
          epsm1_tmp(:,:) = GWPC_CONJG(epsm1_qbz(:,:,io))
          epsm1_trcc_qbz(:,:,io) = TRANSPOSE(epsm1_tmp)
         end do
       end if
     end if !gwcalctyp

     ! Get Fourier components of the Coulombian interaction in the BZ ===
     ! In 3D systems, neglecting umklapp: vc(Sq,sG) = vc(q,G) = 4pi/|q+G|**2
     ! The same relation holds for 0-D systems, but not in 1-D or 2D systems. It depends on S.
     do ig=1,npwc
       vc_sqrt_qbz(Gsph_c%rottb(ig,itim_q,isym_q))=Vcp%vc_sqrt(ig,iq_ibz)
     end do

     call timab(434,2,tsec) ! initq

     ! Sum over band
     call timab(445,1,tsec) ! loop
     do ib=1,Sigp%nbnds
       call timab(436,1,tsec) ! (1)

       ! Parallelism over spin
       ! This processor has this k-point but what about spin?
       if (proc_distrb(ib,ik_bz,spin)/=Wfd%my_rank) CYCLE

       call Wfd%get_ur(ib,ik_ibz,spin,ur_ibz)

       if (Psps%usepaw==1) then
         ! Load cprj for point ksum, this spin or spinor and *THIS* band.
         ! TODO MG I could avoid doing this but I have to exchange spin and bands ???
         ! For sure there is a better way to do this!
         call Wfd%get_cprj(ib,ik_ibz,spin,Cryst,Cprj_ksum,sorted=.FALSE.)
         call paw_symcprj(ik_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,Cprj_ksum)
         if (Dtset%pawcross==1) then
           call Wfdf%paw_get_aeur(ib,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,&
&              ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum)
         end if
       end if

       call timab(436,2,tsec) ! (1)
       call timab(437,1,tsec) ! rho_tw_g

       ! Get all <k-q,ib,s|e^{-i(q+G).r}|s,jb,k>, at once
       do jb=ib1,ib2

         call rho_tw_g(nspinor,npwc,gwc_nfftot,ndat1,gwc_ngfft,1,use_padfft,igfftcg0,gw_gbound,&
&          ur_ibz        ,iik,ktabr(:,ik_bz),ph_mkt  ,spinrot_kbz,  &
&          wfr_bdgw(:,jb),jik,ktabr(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&          nspinor,rhotwg_ki(:,jb))

         if (Psps%usepaw==1) then
           ! Add on-site contribution, projectors are already in BZ !TODO Recheck this!
           i2=jb; if (nspinor==2) i2=(2*jb-1)
           spad=(nspinor-1)
           call paw_rho_tw_g(npwc,nspinor,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_c%gvec,&
&            Cprj_ksum(:,:),Cprj_kgw(:,i2:i2+spad),Pwij_qg,rhotwg_ki(:,jb))
           if (Dtset%pawcross==1) then ! Add paw cross term
             call paw_cross_rho_tw_g(nspinor,npwc,nfftf,rho_ngfft,1,use_padfftf,igfftfcg0,gboundf,&
&             ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum,iik,ktabrf(:,ik_bz),ph_mkt,spinrot_kbz,&
&             ur_ae_bdgw(:,jb),ur_ae_onsite_bdgw(:,jb),ur_ps_onsite_bdgw(:,jb),jik,ktabrf(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&             nspinor,rhotwg_ki(:,jb))
           end if
         end if

         ! Multiply by the square root of the Coulomb term
         ! In 3-D systems, the factor sqrt(4pi) is included)
         do ii=1,nspinor
           spad=(ii-1)*npwc
           rhotwg_ki(spad+1:spad+npwc,jb) = rhotwg_ki(spad+1:spad+npwc,jb)*vc_sqrt_qbz(1:npwc)
         end do

         ! === Treat analytically the case q --> 0 ===
         ! * The oscillator is evaluated at q=O as it is considered constant in the small cube around Gamma
         !   while the Colulomb term is integrated out.
         ! * In the scalar case we have nonzero contribution only if ib==jb
         ! * For nspinor==2 evalute <ib,up|jb,up> and <ib,dwn|jb,dwn>,
         !   impose orthonormalization since npwwfn might be < npwvec.
         if (ik_bz==jk_bz) then
           if (nspinor==1) then
             rhotwg_ki(1,jb)=czero_gw
             if (ib==jb) rhotwg_ki(1,jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp)

           else
             npw_k = Wfd%npwarr(ik_ibz)
             rhotwg_ki(1, jb) = zero; rhotwg_ki(npwc+1, jb) = zero
             if (ib==jb) then
               ABI_CHECK(Wfd%get_wave_ptr(ib, ik_ibz, spin, wave_sum, msg) == 0, msg)
               cg_sum => wave_sum%ug
               ABI_CHECK(Wfd%get_wave_ptr(jb, jk_ibz, spin, wave_jb, msg) == 0, msg)
               cg_jb  => wave_jb%ug
               ctmp = xdotc(npw_k, cg_sum(1:), 1, cg_jb(1:), 1)
               rhotwg_ki(1,jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp) * real(ctmp)
               ctmp = xdotc(npw_k, cg_sum(npw_k+1:), 1, cg_jb(npw_k+1:), 1)
               rhotwg_ki(npwc+1,jb) = CMPLX(SQRT(Vcp%i_sz),0.0_gwp) * real(ctmp)
               ! PAW is missing
             end if
           end if
         end if
       end do !jb  Got all matrix elements from minbnd up to maxbnd.

       theta_mu_minus_e0i=fact_sp*qp_occ(ib,ik_ibz,spin)

       ! Starting point to evaluate the derivative of Sigma and the Spectral function
       e0i=qp_ene(ib,ik_ibz,spin)

       ! Frequencies for the spectral function, e0i=qp_ene(ib,ik_ibz,spin)
       ! FIXME the interval is not centered on eoi ! WHY?
       if (Sr%nomega_r>0) omegame0i(1:Sr%nomega_r)=DBLE(Sr%omega_r(1:Sr%nomega_r))-e0i

       call timab(437,2,tsec) ! rho_tw_g

       call timab(443,1,tsec) ! ac_lrk_appl
      
       if (mod10==SIG_GW_AC) then
         rhotw_epsm1_rhotw(:,:,:) = czero_gw
         do iiw=1,Er%nomega_i
           ABI_MALLOC(epsm1_sqrt_rhotw, (neig(iiw), minbnd:maxbnd))
           ! epsm1_sqrt_rhotw = SQRT(epsm1) * rho_tw
           call xgemm('C','N',neig(iiw),maxbnd-minbnd+1,npwc,cone_gw,ac_epsm1cqwz2(:,:,iiw),npwc,&
&                    rhotwg_ki,npwc,czero_gw,epsm1_sqrt_rhotw,neig(iiw))
           call xherk('L','C',maxbnd-minbnd+1,neig(iiw),one_gw,epsm1_sqrt_rhotw,neig(iiw),zero_gw,rhotw_epsm1_rhotw(:,:,iiw),maxbnd-minbnd+1)

           ! Get the upper part of rhotw_epsm1_rhotw
           ! that is hermitian by construction
           do jb=minbnd,maxbnd
             do kb=jb+1,maxbnd
               rhotw_epsm1_rhotw(jb,kb,iiw) = CONJG( rhotw_epsm1_rhotw(kb,jb,iiw) )
             end do
           end do
           ABI_FREE(epsm1_sqrt_rhotw)
         end do
         call timab(443,2,tsec) ! ac_lrk_appl
       endif

       do kb=ib1,ib2
         call timab(438,1,tsec) ! (2)

         ! Get frequencies $\omega$-\epsilon_in$ to evaluate  $d\Sigma/dE$, note the spin
         ! subtract e_KS since we have stored e_KS+ Delta \omega in Sr%omega4sd, not required for AC
         do io=Sr%nomega_r+1,nomega_tot
           omegame0i(io)=DBLE(Sr%omega4sd(kb,jk_ibz,io-Sr%nomega_r,spin))-e0i
         end do

         ! Get the ket \Sigma|\phi_{k,kb}> according to the method.
         rhotwgp(:)=rhotwg_ki(:,kb)

         SELECT CASE (mod10)
         CASE (SIG_GW_PPM)
           ! GW WITH Plasmon-Pole Model.
           ! Note that ppmodel 3 or 4 work only in case of standard perturbative approach!
           ! Moreover, for ppmodel 3 and 4, spinorial case is not allowed
           sigc_ket  = czero_gw
           call calc_sig_ppm(PPm,nspinor,npwc,nomega_tot,rhotwgp,botsq,otq,&
&           omegame0i,Sigp%zcut,theta_mu_minus_e0i,eig,npwc,sigc_ket,sigcme_3)

           if (PPm%model==3.or.PPm%model==4) then
             sigcme2(:,kb)=sigcme2(:,kb) + (wtqp+wtqm)*DBLE(sigcme_3(:)) + (wtqp-wtqm)*j_gw*AIMAG(sigcme_3(:))
           end if

         CASE (SIG_GW_AC)
           ! GW with Analytic continuation.

           ! This part is so optimized for AC that there is nothing to do here !

         CASE (SIG_GW_CD)
           ! GW with contour deformation.
           ! Check if pole contributions need to be summed
           ! (this avoids unnecessary splint calls and saves time)
           !me_calc_poles = .TRUE.
           sigc_ket  = czero_gw
           do io=1,nomega_tot
             if (omegame0i(io)>=zero.AND.(ABS(one-theta_mu_minus_e0i)>zero)) then
               !me_calc_poles(io) = .TRUE.
               if ( w_maxval(kb) < ABS(omegame0i(io)) ) w_maxval(kb) = ABS(omegame0i(io))
             else if (omegame0i(io)<zero.AND.(ABS(theta_mu_minus_e0i)>zero)) then
               !me_calc_poles(io) = .TRUE.
               if ( w_maxval(kb) < ABS(omegame0i(io)) ) w_maxval(kb) = ABS(omegame0i(io))
             end if
           end do

           ! Check memory saving
           if (Er%mqmem == 0) then
             call calc_sigc_cd(npwc,npwc,nspinor,nomega_tot,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
&             Er%omega,Er%epsm1(:,:,:,1),omegame0i,theta_mu_minus_e0i,sigc_ket,Dtset%ppmfrq,npoles_missing(kb),&
&              method=Dtset%cd_frqim_method)
           else
             call calc_sigc_cd(npwc,npwc,nspinor,nomega_tot,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
&             Er%omega,epsm1_qbz,omegame0i,theta_mu_minus_e0i,sigc_ket,Dtset%ppmfrq,npoles_missing(kb),&
&              method=Dtset%cd_frqim_method)
           end if

#if 0
           if (wtqm/=0) then
             call calc_sigc_cd(npwc,npwc,nspinor,,nomega_tot,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
&              Er%omega,epsm1_trcc_qbz,omegame0i,theta_mu_minus_e0i,aherm_sigc_ket,Dtset%ppmfrq,npoles_missing(kb),&
&               method=Dtset%cd_frqim_method)

             herm_sigc_ket  = half*(sigc_ket + aherm_sigc_ket)
             aherm_sigc_ket = half*(sigc_ket - aherm_sigc_ket)
           else
             herm_sigc_ket  = sigc_ket
             aherm_sigc_ket = czero_gw
           end if
#endif

         CASE (SIG_QPGW_PPM)
           ! MODEL GW calculation WITH PPm  TODO Spinor not tested.
           ! Calculate \Sigma(E_k) |k> to obtain <j|\Sigma(E_k)|k>
           ABI_MALLOC(sigcme_new,(nomega_tot))
           sigc_ket  = czero_gw
           ket1      = czero_gw
           ket2      = czero_gw

           call calc_sig_ppm(PPm,nspinor,npwc,nomega_tot,rhotwgp,botsq,otq,&
&            omegame0i,Sigp%zcut,theta_mu_minus_e0i,eig,npwc,ket1,sigcme_new)

           if (Sigp%gwcalctyp==28) then
             if (PPm%model/=1.and.PPm%model/=2) then
               ! This is needed to have npwc=PPm%dm2_botsq=PPm%dm2_otq
               write(msg,'(3a)')&
&                'For the time being, gwcalctyp=28 cannot be used with ppmodel=3,4.',ch10,&
&                'Use another Plasmon Pole Model when gwcalctyp=28.'
               MSG_ERROR(msg)
             end if
             ABI_MALLOC(botsq_conjg_transp,(PPm%dm2_botsq,npwc))
             botsq_conjg_transp=TRANSPOSE(botsq) ! Keep these two lines separated, otherwise gfortran messes up
             botsq_conjg_transp=CONJG(botsq_conjg_transp)
             ABI_MALLOC(otq_transp,(PPm%dm2_otq,PPm%npwc))
             otq_transp=TRANSPOSE(otq)

             call calc_sig_ppm(PPm,nspinor,npwc,nomega_tot,rhotwgp,botsq_conjg_transp,otq_transp,&
&              omegame0i,Sigp%zcut,theta_mu_minus_e0i,eig,npwc,ket2,sigcme_3)

             ABI_FREE(botsq_conjg_transp)
             ABI_FREE(otq_transp)
             sigc_ket= half*(ket1+ket2)
           else
             sigc_ket= ket1
           end if

           ABI_FREE(sigcme_new)

         CASE (SIG_QPGW_CD)
           ! MODEL GW with numerical integration.
           ! Check if pole contributions need to be summed
           ! (this avoids unnecessary splint calls and saves time)
           !me_calc_poles = .TRUE.
           sigc_ket  = czero_gw
           ket1      = czero_gw
           ket2      = czero_gw

           do io=1,nomega_tot
             if (omegame0i(io)>=zero.AND.(ABS(one-theta_mu_minus_e0i)>zero)) then
               !me_calc_poles(io) = .TRUE.
               if ( w_maxval(kb) < ABS(omegame0i(io)) ) w_maxval(kb) = ABS(omegame0i(io))
             else if (omegame0i(io)<zero.AND.(ABS(theta_mu_minus_e0i)>zero)) then
               !me_calc_poles(io) = .TRUE.
               if ( w_maxval(kb) < ABS(omegame0i(io)) ) w_maxval(kb) = ABS(omegame0i(io))
             end if
           end do

           ! Calculate \Sigma(E_k)|k> to obtain <j|\Sigma(E_k)|k>
           call calc_sigc_cd(npwc,npwc,nspinor,nomega_tot,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
&            Er%omega,epsm1_qbz,omegame0i,theta_mu_minus_e0i,ket1,Dtset%ppmfrq,npoles_missing(kb),&
&             method=Dtset%cd_frqim_method)

           if (Sigp%gwcalctyp==29) then
             ! Calculate \Sigma^*(E_k)|k> to obtain <k|\Sigma(E_k)|j>^*
             call calc_sigc_cd(npwc,npwc,nspinor,nomega_tot,Er%nomega,Er%nomega_r,Er%nomega_i,rhotwgp,&
&              Er%omega,epsm1_trcc_qbz,omegame0i,theta_mu_minus_e0i,ket2,Dtset%ppmfrq,npoles_missing(kb),&
&               method=Dtset%cd_frqim_method)
             sigc_ket = half*(ket1+ket2)
           else
             sigc_ket = ket1
           end if

         CASE DEFAULT
           MSG_ERROR(sjoin("Unsupported value for mod10:", itoa(mod10)))
         END SELECT

         if (Sigp%gwcomp==1) then
           ! TODO spinor not implemented
           call calc_sig_ppm_comp(npwc,nomega_tot,rhotwgp,botsq,otq,DBLE(Sr%egw(kb,jk_ibz,spin)-en_high),&
&            Sigp%zcut,theta_mu_minus_e0i,sigc_ket,PPm%model,npwc,PPm%dm2_botsq,PPm%dm2_otq)
         end if

         call timab(438,2,tsec) !
         call timab(439,1,tsec) ! sigma_me

         ! Loop over the non-zero row elements of this column.
         ! 1) If gwcalctyp<20 : only diagonal elements since QP==KS.
         ! 2) If gwcalctyp>=20: only off-diagonal elements connecting states with same character.
         do irow=1,Sigcij_tab(spin)%col(kb)%size1
           jb = Sigcij_tab(spin)%col(kb)%bidx(irow)
           rhotwg=rhotwg_ki(:,jb)

           ! Calculate <\phi_j|\Sigma_c|\phi_k>
           ! Different freqs according to method (AC or Perturbative), see nomega_sigc.
           if (mod10==SIG_GW_AC) then
             sigctmp(:,:) = czero_gw
             do iab=1,Sigp%nsig_ab
               do io=1,nomega_sigc
                 omegame0i_ac  = Sr%omega_i(io)-qp_ene(ib,ik_ibz,spin)
                 omegame0i2_ac = omegame0i_ac*omegame0i_ac
                 do iiw=1,Er%nomega_i
                   sigctmp(io,iab) = sigctmp(io,iab) + piinv * rhotw_epsm1_rhotw(jb,kb,iiw) * omegame0i_ac &
&                                                      /(omegame0i2_ac + omegap2(iiw))
                 end do
               end do
             end do
           else
             do iab=1,Sigp%nsig_ab
               spadc1 = spinor_padc(1, iab); spadc2 = spinor_padc(2, iab)
               do io=1,nomega_sigc
                 sigctmp(io,iab) = XDOTC(npwc,rhotwg(spadc1+1:),1,sigc_ket(spadc2+1:,io),1)
               end do
             end do
           end if

           if (Sigp%gwcomp==1) then
             ! Evaluate Extrapolar term TODO this does not work with spinor
             if (extrapolar_distrb(jb,kb,ik_bz,spin) == Wfd%my_rank ) then
               ! Do it once as it does not depend on the ib index summed over.
               extrapolar_distrb(jb,kb,ik_bz,spin) = xmpi_undefined_rank
#if 1
               call calc_wfwfg(ktabr(:,jk_ibz),jik, spinrot_kgw, & ! TODO: why jk_ibz?
                 gwc_nfftot,nspinor,gwc_ngfft,wfr_bdgw(:,jb),wfr_bdgw(:,kb),wf1swf2_g)
#else
               call calc_wfwfg(ktabr(:,jk_bz),jik,spinrot_kgw,&
                 gwc_nfftot,nspinor,gwc_ngfft,wfr_bdgw(:,jb),wfr_bdgw(:,kb),wf1swf2_g)
#endif

               if (Psps%usepaw==1) then
                 i1=jb; i2=kb
                 if (nspinor==2) then
                   i1=(2*jb-1); i2=(2*kb-1)
                 end if
                 spad=(nspinor-1)
                 call paw_rho_tw_g(gwc_nfftot,Sigp%nsig_ab,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,&
&                  gw_gfft,Cprj_kgw(:,i1:i1+spad),Cprj_kgw(:,i2:i2+spad),Pwij_fft,wf1swf2_g)

                 if (Dtset%pawcross==1) then ! Add paw cross term
                   call paw_cross_rho_tw_g(nspinor,npwc,nfftf,rho_ngfft,1,use_padfftf,igfftfcg0,gboundf,&
&                   ur_ae_bdgw(:,jb),ur_ae_onsite_bdgw(:,jb),ur_ps_onsite_bdgw(:,jb),jik,ktabrf(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&                   ur_ae_bdgw(:,kb),ur_ae_onsite_bdgw(:,kb),ur_ps_onsite_bdgw(:,kb),jik,ktabrf(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&                   nspinor,wf1swf2_g)
                 end if
               end if

               ! The static contribution from completeness relation is calculated once ===
               call calc_coh_comp(iq_ibz,Vcp%i_sz,(jb==kb),nspinor,Sigp%nsig_ab,DBLE(Sr%egw(kb,jk_ibz,spin)-en_high),&
&                npwc,Gsph_c%gvec,gwc_ngfft,gwc_nfftot,wf1swf2_g,vc_sqrt_qbz,botsq,otq,sigcohme)

               do io=1,nomega_sigc
                 sigctmp(io,:) = sigctmp(io,:)+sigcohme(:)
               end do
             end if ! gwcomp==1
           end if ! gwcom==1

           ! Accumulate and, in case, symmetrize matrix elements of Sigma_c
           do iab=1,Sigp%nsig_ab
             is_idx=spin; if (nspinor==2) is_idx=iab

             sigcme_tmp(:,jb,kb,is_idx)=sigcme_tmp(:,jb,kb,is_idx) + &
&              (wtqp+wtqm)*DBLE(sigctmp(:,iab)) + (wtqp-wtqm)*j_gw*AIMAG(sigctmp(:,iab))

             sigc(1,:,jb,kb,is_idx)=sigc(1,:,jb,kb,is_idx) + wtqp*      sigctmp(:,iab)
             sigc(2,:,jb,kb,is_idx)=sigc(2,:,jb,kb,is_idx) + wtqm*CONJG(sigctmp(:,iab))
             ! TODO this should be the contribution coming from the anti-hermitian part.
           end do
         end do ! irow used to calculate matrix elements of $\Sigma$

         ! shaltaf (030406): this has to be done in a clean way later.
         ! TODO does not work with spinor.
         if (mod10==SIG_GW_PPM.and.(PPm%model==3.or.PPm%model==4)) then
           sigcme_tmp(:,kb,kb,spin)= sigcme2(:,kb)
           sigc(1,:,kb,kb,spin)= sigcme2(:,kb)
           sigc(2,:,kb,kb,spin)= czero
         end if

         call timab(439,2,tsec) ! csigme(SigC)
       end do !kb to calculate matrix elements of $\Sigma$
     end do !ib

     call timab(445,2,tsec) ! csigme(SigC)

     ! Deallocate k-dependent quantities.
     ABI_FREE(gw_gbound)
     if (Dtset%pawcross==1) then
       ABI_FREE(gboundf)
     end if

     if (sigma_needs_ppm(Sigp)) then
       ABI_FREE(botsq)
       ABI_FREE(otq)
       ABI_FREE(eig)
     end if
     if (Psps%usepaw==1) then
       call pawpwij_free(Pwij_qg)
       ABI_FREE(Pwij_qg)
     end if

   end do ! ik_bz

   ABI_FREE(wfr_bdgw)
   if (Wfd%usepaw==1) then
     call pawcprj_free(Cprj_kgw)
     ABI_FREE(Cprj_kgw)
     if (Dtset%pawcross==1) then
       ABI_FREE(ur_ae_bdgw)
       ABI_FREE(ur_ae_onsite_bdgw)
       ABI_FREE(ur_ps_onsite_bdgw)
     end if
   end if
 end do ! spin

 ABI_FREE(sigcme2)
 ABI_FREE(sigcme_3)
 ABI_FREE(igfftcg0)
 if (Dtset%pawcross==1) then
   ABI_FREE(igfftfcg0)
 end if

 ! Gather contributions from all the CPUs
 call timab(440,1,tsec) ! wfd_barrier
 call timab(440,2,tsec) ! wfd_barrier
 call timab(441,1,tsec) ! xmpi_sum

 call xmpi_sum(sigcme_tmp, Wfd%comm, ierr)
 call xmpi_sum(sigc, Wfd%comm, ierr)
 call timab(441,2,tsec) ! xmpi_sum

 ! Multiply by constants. In 3D systems sqrt(4pi) is included in vc_sqrt_qbz ===
 sigcme_tmp = sigcme_tmp /(Cryst%ucvol*Kmesh%nbz)
 sigc       = sigc       /(Cryst%ucvol*Kmesh%nbz)

 ! If we have summed over the IBZ_q now we have to average over complexes ===
 ! Presently only diagonal terms are considered
 ! TODO QP-SCGW required a more involved approach, there is a check in sigma
 ! TODO it does not work if nspinor==2.
 call timab(442,1,tsec) ! final ops

 do spin=1,Wfd%nsppol
   if (can_symmetrize(spin)) then
     if (mod10==SIG_GW_AC) then ! FIXME here there is a problem in case of AC with symmetries
       ABI_MALLOC(sym_cme, (Sr%nomega_i, ib1:ib2, ib1:ib2, Sigp%nsig_ab))
     else
       ABI_MALLOC(sym_cme, (nomega_tot, ib1:ib2, ib1:ib2, Sigp%nsig_ab))
     end if
     sym_cme=czero

     ! Average over degenerate diagonal elements
     ! NOTE: frequencies for \Sigma_c(\omega) should be equal to avoid spurious results.
     ! another good reason to use a strict criterion for the tolerance on eigenvalues.
     do ib=ib1,ib2
       ndegs=0
       do jb=ib1,ib2
         if (degtab(ib,jb,spin)==1) then
           if (nspinor == 1) then
             sym_cme(:, ib, ib, 1) = sym_cme(:, ib, ib, 1) + SUM(sigc(:,:,jb,jb,spin), DIM=1)
           else
             do ii=1,Sigp%nsig_ab
               sym_cme(:, ib, ib, ii) = sym_cme(:, ib, ib, ii) + SUM(sigc(:,:,jb,jb,ii), dim=1)
             end do
           end if
         end if
         ndegs = ndegs + degtab(ib,jb,spin)
       end do
       sym_cme(:,ib,ib,:) = sym_cme(:,ib,ib,:) / ndegs
     end do

     if (Sigp%gwcalctyp >= 20) then
       do iwc=1,nomega_sigc
         call esymm_symmetrize_mels(QP_sym(spin),ib1,ib2,sigc(:,iwc,:,:,spin),sym_cme(iwc,:,:,1))
       end do
     end if

     ! Copy symmetrized values
     do ib=ib1,ib2
       do jb=ib1,ib2
         !if (mod10==SIG_GW_AC.and.average_real) CYCLE ! this is to check another scheme in case of AC
         if (nspinor == 1) then
           sigcme_tmp(:,ib,jb,spin) = sym_cme(:,ib,jb,1)
         else
           do ii=1,Sigp%nsig_ab
             sigcme_tmp(:,ib,jb,ii) = sym_cme(:,ib,jb,ii)
           end do
         end if
       end do
     end do
     ABI_FREE(sym_cme)
   end if
 end do

 ! Reconstruct the full sigma matrix from the upper triangle (only for HF, SEX and COHSEX)
 !if (Sigp%gwcalctyp>=20 .and. sigma_is_herm(Sigp) ) then
 !  ABI_CHECK(nspinor==1,"cannot hermitianize non-collinear sigma!")
 !  do spin=1,Wfd%nsppol
 !    do io=1,nomega_sigc
 !      call hermitianize(sigcme_tmp(io,:,:,spin),"Upper")
 !    end do
 !  end do
 !end if

 ! GW with contour deformation: check on the number of poles not included.
 if (ANY(mod10 == [SIG_GW_CD,SIG_QPGW_CD])) then
   call xmpi_sum(npoles_missing, Wfd%comm, ierr)
   npls = SUM(npoles_missing)
   if (npls>0) then
     MSG_WARNING(sjoin("Total number of missing poles for contour deformation method:", itoa(npls)))
     do band=minbnd,maxbnd
       npls = npoles_missing(band)
       if (npls>0) then
         write(msg,'(a,2(i0,a))')" For band ",band," there are ",npls," missing poles"
         call wrtout(std_out,msg,"COLL")
       end if
     end do
   end if
   ! Print data on the maximum value needed for the screening along the real axis
   w_localmax = MAXVAL(w_maxval)
   call xmpi_max(w_localmax,w_max, Wfd%comm, ierr)
   write(msg,'(a,f12.5,a)') ' Max omega value used in W(omega): ',w_max*Ha_eV,' [eV]'
   call wrtout(std_out,msg,"COLL")
 end if
 call timab(442,2,tsec) ! final ops

 ! ===========================
 ! ==== Deallocate memory ====
 ! ===========================
 if (Psps%usepaw==1) then
   if (allocated(gw_gfft)) then
     ABI_FREE(gw_gfft)
   end if
   call pawcprj_free(Cprj_ksum)
   ABI_FREE(Cprj_ksum)
   if (allocated(Pwij_fft)) then
     call pawpwij_free(Pwij_fft)
     ABI_FREE(Pwij_fft)
   end if
   if (Dtset%pawcross==1) then
     ABI_FREE(ur_ae_sum)
     ABI_FREE(ur_ae_onsite_sum)
     ABI_FREE(ur_ps_onsite_sum)
     ABI_FREE(ktabrf)
   end if
 end if

 ABI_FREE(npoles_missing)
 ABI_FREE(ur_ibz)
 ABI_FREE(usr_bz)
 ABI_FREE(ktabr)
 ABI_FREE(sigc_ket)
 ABI_FREE(rhotwg_ki)
 ABI_FREE(rhotwg)
 ABI_FREE(rhotwgp)
 ABI_FREE(vc_sqrt_qbz)
 ABI_FREE(omegame0i)
 ABI_FREE(sigctmp)
 ABI_FREE(sigc)
 ABI_FREE(w_maxval)

 if (allocated(ket1)) then
   ABI_FREE(ket1)
 endif
 if (allocated(ket2)) then
   ABI_FREE(ket2)
 endif
 if (allocated(epsm1_qbz)) then
   ABI_FREE(epsm1_qbz)
 end if
 if (allocated(epsm1_trcc_qbz)) then
   ABI_FREE(epsm1_trcc_qbz)
 end if
 if (allocated(epsm1_tmp)) then
   ABI_FREE(epsm1_tmp)
 end if
 if (allocated(degtab)) then
   ABI_FREE(degtab)
 end if
 if (allocated(ac_epsm1cqwz2)) then
   ABI_FREE(ac_epsm1cqwz2)
 end if
 if (allocated(rhotw_epsm1_rhotw)) then
   ABI_FREE(rhotw_epsm1_rhotw)
 end if
 if (allocated(aherm_sigc_ket)) then
   ABI_FREE(aherm_sigc_ket)
 end if
 if (allocated(herm_sigc_ket)) then
   ABI_FREE(herm_sigc_ket)
 end if
 if (Sigp%gwcomp==1) then
   ABI_FREE(wf1swf2_g)
 end if
 if (Sigp%gwcomp==1) then
   ABI_FREE(extrapolar_distrb)
 end if
 ABI_FREE(proc_distrb)

 call timab(431,2,tsec)
 call timab(424,2,tsec) ! calc_sigc_me

 call cwtime(cpu_time,wall_time,gflops,"stop")
 write(std_out,'(2(a,f9.1))')" cpu_time = ",cpu_time,", wall_time = ",wall_time

 DBG_EXIT("COLL")

end subroutine calc_sigc_me
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/calc_coh_comp
!! NAME
!! calc_coh_comp
!!
!! FUNCTION
!!  Calculates the COH-like contribution to the self-energy when
!!  the extrapolar technique and the closure relation is used to
!!  reduce the number of empty states to be summed over in the Green
!!  function entering the definition of the GW self-energy.
!!
!! INPUTS
!! iqibz=index of the irreducible q-point in the array qibz, point which is
!!  related by a symmetry operation to the point q summed over (see csigme).
!!  This index is also used to treat the integrable coulombian singularity at q=0
!! ngfft(18)=contain all needed information about 3D FFT for GW wavefuntions,
!!  see ~abinit/doc/variables/vargs.htm#ngfft
!! nsig_ab=Number of components in the self-energy operator (1 for collinear magnetism)
!! npwc=number of plane waves in $\tilde epsilon^{-1}$
!! nspinor=Number of spinorial components.
!! i_sz=contribution arising from the integrable coulombian singularity at q==0
!! (see csigme for the method used), note that in case of 3-D systems the factor
!! 4pi in the coulombian potential is included in the definition of i_sz
!! gvec(3,npwc)=G vectors in reduced coordinates
!! vc_sqrt(npwc)= square root of the coulombian matrix elements for this q-point
!! botsq = Plasmon-pole parameters
!! otq  = PPm parameters
!!
!! OUTPUT
!! sigcohme=partial contribution to the matrix element of $<jb k|\Sigma_{COH}| kb k>$
!!  coming from this single q-point for completeness trick
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      calc_sigc_me
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine calc_coh_comp(iqibz,i_sz,same_band,nspinor,nsig_ab,ediff,npwc,gvec,&
&  ngfft,nfftot,wfg2_jk,vc_sqrt,botsq,otq,sigcohme)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,npwc,nsig_ab,nspinor,nfftot
 real(dp),intent(in) :: i_sz,ediff
 logical,intent(in) :: same_band
!arrays
 integer,intent(in) :: gvec(3,npwc),ngfft(18)
 complex(gwpc),intent(in) :: botsq(npwc,npwc),otq(npwc,npwc)
 complex(gwpc),intent(in) :: vc_sqrt(npwc)
 complex(gwpc),intent(in) :: wfg2_jk(nfftot*nsig_ab)
 complex(gwpc),intent(out) :: sigcohme(nsig_ab)

!Local variables-------------------------------
!scalars
 integer,save :: enough=0
 integer :: ig,ig4,ig4x,ig4y,ig4z,igp,igmin,ispinor,ngfft1,ngfft2,ngfft3,spad,outofbox
!arrays
 integer :: g2mg1(3)

! *************************************************************************

 DBG_ENTER("COLL")

 ! === Treat the case q --> 0 adequately ===
 ! TODO Better treatment of wings
 igmin=1 ; if (iqibz==1) igmin=2
 !
 ! === Partial contribution to the matrix element of Sigma_c ===
 ! * For nspinor==2, the closure relation reads:
 !  $\sum_s \psi_a^*(1)\psi_b(2) = \delta_{ab} \delta(1-2)$
 !  where a,b are the spinor components. As a consequence, Sigma_{COH} is always
 !  diagonal in spin-space and only diagonal matrix elements have to be calculated.
 ! MG  TODO wfg2_jk should be calculated on an augmented FFT box to avoid spurious wrapping of G1-G2.
 !
 ngfft1 = ngfft(1); ngfft2 = ngfft(2); ngfft3 = ngfft(3)
 sigcohme(:) = czero_gw

 do ispinor=1,nspinor
  spad=(ispinor-1)*nfftot
  outofbox=0

   do igp=igmin,npwc
     do ig=igmin,npwc

      g2mg1 = gvec(:,igp)-gvec(:,ig)
      if (ANY(g2mg1(:)>ngfft(1:3)/2) .or. ANY(g2mg1(:)<-(ngfft(1:3)-1)/2)) then
        outofbox = outofbox+1; CYCLE
      end if

      ig4x=MODULO(g2mg1(1),ngfft1)
      ig4y=MODULO(g2mg1(2),ngfft2)
      ig4z=MODULO(g2mg1(3),ngfft3)
      ig4= 1+ig4x+ig4y*ngfft1+ig4z*ngfft1*ngfft2

      !MG where is neta here, ediff, otq might be close to zero depending on gwecomp
      sigcohme(ispinor) = sigcohme(ispinor) + &
&       half*wfg2_jk(spad+ig4)*vc_sqrt(ig)*vc_sqrt(igp) * botsq(ig,igp) / ( otq(ig,igp) * ( ediff -otq(ig,igp) ) )
     end do
   end do

   if (iqibz==1.and.same_band) then
     sigcohme(ispinor) = sigcohme(ispinor) + half*wfg2_jk(spad+ig4)*i_sz*botsq(1,1) / ( otq(1,1) * (ediff -otq(1,1)) )
   end if
 end do !ispinor

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=50) then
     MSG_WARNING(sjoin('Number of G1-G2 pairs outside the G-sphere for Wfns: ',itoa(outofbox)))
     if (enough==50) then
       call wrtout(std_out,' ========== Stop writing Warnings ==========')
     end if
   end if
 end if

 DBG_EXIT("COLL")

end subroutine calc_coh_comp
!!***

!!****f* ABINIT/calc_sigc_cd
!! NAME
!! calc_sigc_cd
!!
!! FUNCTION
!! Calculate contributions to the self-energy operator with the contour deformation method.
!!
!! INPUTS
!!  nomega=Total number of frequencies where $\Sigma_c$ matrix elements are evaluated.
!!  nomegae=Number of frequencies where $\epsilon^{-1}$ has been evaluated.
!!  nomegaei=Number of imaginary frequencies for $\epsilon^{-1}$ (non zero).
!!  nomegaer=Number of real frequencies for $\epsilon^{-1}$
!!  npwc=Number of G vectors for the correlation part.
!!  npwx=Number of G vectors in rhotwgp for each spinorial component.
!!  nspinor=Number of spinorial components.
!!  theta_mu_minus_e0i=1 if e0i is occupied, 0 otherwise. Fractional occupancy in case of metals.
!!  omegame0i(nomega)=Contains $\omega-\epsilon_{k-q,b1,\sigma}$
!!  epsm1q(npwc,npwc,nomegae)=Symmetrized inverse dielectric matrix (exchange part is subtracted).
!!  omega(nomegae)=Set of frequencies for $\epsilon^{-1}$.
!!  rhotwgp(npwx*nspinor)=Matrix elements: $<k-q,b1,\sigma|e^{-i(q+G)r} |k,b2,\sigma>*vc_sqrt$
!!
!! OUTPUT
!! ket(npwc,nomega)=Contains \Sigma_c(\omega)|\phi> in reciprocal space.
!!
!! SIDE EFFECTS
!! npoles_missing=Incremented with the number of poles whose contribution has not been taken into account due to
!!  limited frequency mesh used for W.
!!
!! PARENTS
!!      calc_sigc_me,m_screen
!!
!! CHILDREN
!!      spline,splint,xgemm,xgemv
!!
!! SOURCE

subroutine calc_sigc_cd(npwc,npwx,nspinor,nomega,nomegae,nomegaer,nomegaei,rhotwgp,&
& omega,epsm1q,omegame0i,theta_mu_minus_e0i,ket,plasmafreq,npoles_missing,calc_poles,method)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,nomegae,nomegaei,nomegaer,npwc,npwx,nspinor
 integer,intent(inout) :: npoles_missing
 real(dp),intent(in) :: theta_mu_minus_e0i,plasmafreq
!arrays
 real(dp),intent(in) :: omegame0i(nomega)
 complex(dpc),intent(in) :: omega(nomegae)
 complex(gwpc),intent(in) :: epsm1q(npwc,npwc,nomegae)
 complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 complex(gwpc),intent(inout) :: ket(nspinor*npwc,nomega)
 logical, intent(in), optional :: calc_poles(nomega)
 integer, intent(in), optional :: method

!Local variables-------------------------------
!scalars
 integer, parameter :: FABIEN=1,TRAPEZOID=2,NSPLINE=3
 integer :: ii,ig,io,ios,ispinor,spadc,spadx,my_err,ierr,GK_LEVEL,INTMETHOD
 integer :: i,j
 real(dp) :: rt_imag,rt_real,local_one,local_zero
 real(dp) :: intsign,temp1,temp2,temp3,temp4
 real(dp) :: alph,inv_alph,beta,alphsq,betasq,inv_beta
 real(dp) :: re_intG,re_intK,im_intG,im_intK,GKttab,tau,ttil
 real(dp) :: ref,imf,r,s,r2,s2
 complex(dpc) :: ct,domegaleft,domegaright
 complex(gwpc) :: fact
!arrays
 real(dp) :: omegame0i_tmp(nomega),tmp_x(2),tmp_y(2)
 real(dp) :: left(nomega),right(nomega)
 real(dp) :: tbeta(nomega),tinv_beta(nomega),tbetasq(nomega)
 real(dp) :: atermr(nomega),aterml(nomega),logup(nomega),logdown(nomega)
 real(dp) :: rtmp_r(nomegaer),rtmp_i(nomegaer)
 real(dp) :: ftab(nomegaei+2),ftab2(nomegaei+2),xtab(nomegaei+2),y(3,nomegaei+2)
 real(dp) :: work(nomegaei+2),work2(nomegaei+2),y2(3,nomegaei+2)
 complex(dpc) :: omega_imag(nomegaei+1)
 complex(gwpc) :: epsrho(npwc,nomegae),epsrho_imag(npwc,nomegaei+1)
 complex(gwpc) :: tfone(npwc,nomegaei+1),tftwo(npwc,nomegaei+1)
 complex(gwpc) :: weight(nomegaei+1,nomega)
 complex(gwpc) :: weight2(nomegaei,nomega)
 logical :: my_calc_poles(nomega)
 real(dp), allocatable :: KronN(:),KronW(:),GaussW(:),fint(:),fint2(:)
!*************************************************************************

 my_calc_poles=.TRUE.; my_err=0

 ! Set integration method for imaginary axis
 INTMETHOD = FABIEN
 if (present(method)) then
   if (method==1) INTMETHOD = FABIEN
   if (method==2) INTMETHOD = TRAPEZOID
   if (method>2) then
     INTMETHOD = NSPLINE
     if (method==3) then
       GK_LEVEL = 15
       ABI_ALLOCATE(KronN,(GK_LEVEL))
       ABI_ALLOCATE(KronW,(GK_LEVEL))
       ABI_ALLOCATE(GaussW,(GK_LEVEL-8))
       ABI_ALLOCATE(fint,(GK_LEVEL))
       ABI_ALLOCATE(fint2,(GK_LEVEL))
       KronN(:) = Kron15N(:); KronW(:) = Kron15W(:); GaussW(:) = Gau7W(:)
     else if (method==4) then
       GK_LEVEL = 23
       ABI_ALLOCATE(KronN,(GK_LEVEL))
       ABI_ALLOCATE(KronW,(GK_LEVEL))
       ABI_ALLOCATE(GaussW,(GK_LEVEL-12))
       ABI_ALLOCATE(fint,(GK_LEVEL))
       ABI_ALLOCATE(fint2,(GK_LEVEL))
       KronN(:) = Kron23N(:); KronW(:) = Kron23W(:); GaussW(:) = Gau11W(:)
     else if (method>4) then
       GK_LEVEL = 31
       ABI_ALLOCATE(KronN,(GK_LEVEL))
       ABI_ALLOCATE(KronW,(GK_LEVEL))
       ABI_ALLOCATE(GaussW,(GK_LEVEL-16))
       ABI_ALLOCATE(fint,(GK_LEVEL))
       ABI_ALLOCATE(fint2,(GK_LEVEL))
       KronN(:) = Kron31N(:); KronW(:) = Kron31W(:); GaussW(:) = Gau15W(:)
     end if
   end if
 end if

 ! Avoid divergences in $\omega - \omega_s$.
 omegame0i_tmp(:)=omegame0i(:)
 do ios=1,nomega
   if (ABS(omegame0i_tmp(ios))<tol6) omegame0i_tmp(ios)=sign(tol6,omegame0i_tmp(ios))
 end do

 do ispinor=1,nspinor
   spadx=(ispinor-1)*npwx; spadc=(ispinor-1)*npwc

   ! Calculate $ \sum_{Gp} (\epsilon^{-1}_{G Gp}(\omega)-\delta_{G Gp}) \rhotwgp(Gp) $
!$omp parallel do
   do io=1,nomegae
     call XGEMV('N',npwc,npwc,cone_gw,epsm1q(:,:,io),npwc,rhotwgp(1+spadx:),1,czero_gw,epsrho(:,io),1)
   end do

   ! Integrand along the imaginary axis.
   epsrho_imag(:,1)=epsrho(:,1)
   epsrho_imag(:,2:nomegaei+1)=epsrho(:,nomegaer+1:nomegae)

   ! Frequency mesh for integral along the imaginary axis.
   omega_imag(1)=omega(1)
   omega_imag(2:nomegaei+1)=omega(nomegaer+1:nomegae)

   ! Original implementation -- saved here for reference during development
   ! === Perform integration along the imaginary axis ===
   !do io=1,nomegaei+1
   !  if (io==1) then
   !    domegaleft  = omega_imag(io)
   !    domegaright =(omega_imag(io+1)-omega_imag(io  ))*half
   !  else if (io==nomegaei+1) then
   !    domegaleft  =(omega_imag(io  )-omega_imag(io-1))*half
   !    domegaright =(omega_imag(io  )-omega_imag(io-1))*half
   !  else
   !    domegaleft  =(omega_imag(io  )-omega_imag(io-1))*half
   !    domegaright =(omega_imag(io+1)-omega_imag(io  ))*half
   !  end if
   !  do ios=1,nomega
   !    omg2 = -AIMAG(omega_imag(io)+domegaright)/REAL(omegame0i_tmp(ios))
   !    omg1 = -AIMAG(omega_imag(io)-domegaleft )/REAL(omegame0i_tmp(ios))
   !    fact = ATAN(omg2)-ATAN(omg1)
   !    ket(spadc+1:spadc+npwc,ios)=ket(spadc+1:spadc+npwc,ios)+epsrho_imag(:,io)*fact
   !  end do
   !end do !io

   !ket(spadc+1:spadc+npwc,:)=ket(spadc+1:spadc+npwc,:)/pi
   ! ---------------- end of original implementation -----------------------

   select case(INTMETHOD)
   case(FABIEN)
     ! Hopefully more effective implementation MS 04.08.2011
     ! Perform integration along imaginary axis using BLAS
     ! First calculate first and last weights
     weight(1,:) = ATAN(-half*AIMAG(omega_imag(2))/REAL(omegame0i_tmp(:)))
     domegaleft  = (three*omega_imag(nomegaei+1)-omega_imag(nomegaei))
     domegaright = (omega_imag(nomegaei+1)+omega_imag(nomegaei))
     right(:)    = -AIMAG(omega_imag(nomegaei+1)-omega_imag(nomegaei))*REAL(omegame0i_tmp(:))
     left(:)     = quarter*AIMAG(domegaleft)*AIMAG(domegaright) &
&                   +REAL(omegame0i_tmp(:))*REAL(omegame0i_tmp(:))
     weight(nomegaei+1,:) = ATAN(right(:)/left(:))
     ! Calculate the rest of the weights
     do io=2,nomegaei
       domegaleft  = (omega_imag(io  )+omega_imag(io-1))
       domegaright = (omega_imag(io+1)+omega_imag(io  ))
       right(:)    = -half*AIMAG(omega_imag(io+1)-omega_imag(io-1))*REAL(omegame0i_tmp(:))
       left(:)     = REAL(omegame0i_tmp(:))*REAL(omegame0i_tmp(:)) &
&       +quarter*AIMAG(domegaleft)*AIMAG(domegaright)
       weight(io,:) = ATAN(right(:)/left(:))
     end do

     ! Use BLAS call to perform matrix-matrix multiplication and accumulation
     fact = CMPLX(piinv,zero)

     call xgemm('N','N',npwc,nomega,nomegaei+1,fact,epsrho_imag,npwc,&
&     weight,nomegaei+1,cone_gw,ket(spadc+1:spadc+npwc,:),npwc)

   case(TRAPEZOID)
!   Trapezoidal rule Transform omega coordinates
     alph     = plasmafreq
     alphsq   = alph*alph
     inv_alph = one/alph

     xtab(1:nomegaei+1) = AIMAG(omega_imag(:))/(AIMAG(omega_imag(:)) + alph)
     xtab(nomegaei+2)   = one

!   Efficient trapezoidal rule with BLAS calls
     tbeta(:)     = REAL(omegame0i_tmp(:))
     tbetasq(:)   = tbeta(:)*tbeta(:)
     tinv_beta(:) = one/tbeta(:)

     do io=1,nomegaei
       atermr(:)    = inv_alph*tinv_beta(:)*((alphsq+tbetasq(:))*xtab(io+1)-tbetasq(:))
       aterml(:)    = inv_alph*tinv_beta(:)*((alphsq+tbetasq(:))*xtab(io  )-tbetasq(:))
       right(:)     = ATAN((atermr(:)-aterml(:))/(one+atermr(:)*aterml(:)))
       logup(:)     = ABS(((alphsq+tbetasq(:))*xtab(io+1)-two*tbetasq(:)) &
&                     *xtab(io+1)+tbetasq(:))
       logdown(:)   = ABS(((alphsq+tbetasq(:))*xtab(io  )-two*tbetasq(:)) &
&                     *xtab(io  )+tbetasq(:))
       ! Trapezoid integration weights
       weight(io,:)  = CMPLX(-(half*alph*tbeta(:)*LOG(logup(:)/logdown(:)) + tbetasq(:) &
&                         *right(:))/(alphsq+tbetasq(:)),zero)
       weight2(io,:) = CMPLX(-right(:),zero)
       ! Linear interpolation coefficients for each section (sum over ig)
       tfone(:,io)   = (epsrho_imag(:,io+1)-epsrho_imag(:,io)) &
&                      /(xtab(io+1)-xtab(io))
       tftwo(:,io)   = epsrho_imag(:,io) - tfone(:,io)*xtab(io)
     end do

     ! Calculate weights for asymptotic behaviour
     atermr(:)   = alph*tinv_beta(:)
     aterml(:)   = inv_alph*tinv_beta(:)*((alphsq+tbetasq(:))*xtab(nomegaei+1)-tbetasq(:))
     logup(:)    = alphsq*xtab(nomegaei+1)*xtab(nomegaei+1)
     logdown(:)  = ABS(((alphsq+tbetasq(:))*xtab(nomegaei+1)-two*tbetasq(:)) &
&                   *xtab(nomegaei+1)+tbetasq(:))
     right(:)     = ATAN((atermr(:)-aterml(:))/(one+atermr(:)*aterml(:)))
     weight (nomegaei+1,:) = CMPLX(-(half*(alphsq*tinv_beta(:)*LOG(logdown(:)/logup(:)) &
&     - tbeta(:)*LOG(xtab(nomegaei+1)*xtab(nomegaei+1))) - alph*right(:)),zero)
     tfone(:,nomegaei+1) = -(zero-epsrho_imag(:,nomegaei+1)*AIMAG(omega_imag(nomegaei+1))) &
                           /(one-xtab(nomegaei+1))

     ! Use BLAS call to perform matrix-matrix multiplication and accumulation
     fact = CMPLX(piinv,zero)

     call xgemm('N','N',npwc,nomega,nomegaei+1,fact,tfone,npwc,&
&     weight ,nomegaei+1,cone_gw,ket(spadc+1:spadc+npwc,:),npwc)
     call xgemm('N','N',npwc,nomega,nomegaei  ,fact,tftwo,npwc,&
&     weight2,nomegaei  ,cone_gw,ket(spadc+1:spadc+npwc,:),npwc)

   case(NSPLINE)
     ! Natural spline followed by Gauss-Kronrod
     ! Transform omega coordinates
     alph     = plasmafreq
     alphsq   = alph*alph
     inv_alph = one/alph

     xtab(1:nomegaei+1) = AIMAG(omega_imag(:))/(AIMAG(omega_imag(:)) + alph)
     xtab(nomegaei+2)   = one

! Gauss-Kronrod integration of spline fit of f(t)/(1-t) in transformed space
! *** OPENMP SECTION *** Added by MS
!!$OMP PARALLEL DO PRIVATE(ig,ftab,ftab2,s,s2,r,r2,y,y2,work,work2,beta,betasq,inv_beta, &
!!$OMP  intsign,io,ii,i,j,re_intG,re_intK,im_intG,im_intK,temp1,temp2,temp3,temp4, &
!!$OMP  ttil,tau,ref,fint,imf,fint2,GKttab)
     do ig=1,npwc
       ! Spline fit
       ftab (1:nomegaei+1) =  REAL(epsrho_imag(ig,1:nomegaei+1))/(one-xtab(1:nomegaei+1))
       ftab2(1:nomegaei+1) = AIMAG(epsrho_imag(ig,1:nomegaei+1))/(one-xtab(1:nomegaei+1))
       ftab (nomegaei+2)   = zero; ftab2(nomegaei+2) = zero
       ! Explicit calculation of spline coefficients
       s  = zero; s2 = zero
       do i = 1, nomegaei+2-1
         r  = ( ftab (i+1) - ftab (i) ) / ( xtab(i+1) - xtab(i) )
         r2 = ( ftab2(i+1) - ftab2(i) ) / ( xtab(i+1) - xtab(i) )
         y (2,i) = r  - s; y2(2,i) = r2 - s2
         s  = r; s2 = r2
       end do
       s = zero; s2 = zero
       r = zero; r2 = zero
       y(2,1) = zero; y2(2,1) = zero
       y(2,nomegaei+2) = zero; y2(2,nomegaei+2) = zero
       do i = 2, nomegaei+2-1
         y (2,i) = y (2,i) + r  * y (2,i-1)
         y2(2,i) = y2(2,i) + r2 * y2(2,i-1)
         work (i) = two * ( xtab(i-1) - xtab(i+1) ) - r  * s
         work2(i) = two * ( xtab(i-1) - xtab(i+1) ) - r2 * s2
         s = xtab(i+1) - xtab(i)
         s2 = s
         r  = s  / work (i)
         r2 = s2 / work2(i)
       end do
       do j = 2, nomegaei+2-1
         i = nomegaei+2+1-j
         y (2,i) = ( ( xtab(i+1) - xtab(i) ) * y (2,i+1) - y (2,i) ) / work (i)
         y2(2,i) = ( ( xtab(i+1) - xtab(i) ) * y2(2,i+1) - y2(2,i) ) / work2(i)
       end do
       do i = 1, nomegaei+2-1
         s = xtab(i+1) - xtab(i); s2 = s;
         r = y(2,i+1) - y(2,i); r2 = y2(2,i+1) - y2(2,i);
         y(3,i) = r / s; y2(3,i) = r2 / s2;
         y(2,i) = three * y(2,i); y2(2,i) = three * y2(2,i);
         y (1,i) = ( ftab (i+1) - ftab (i) ) / s  - ( y (2,i) + r  ) * s
         y2(1,i) = ( ftab2(i+1) - ftab2(i) ) / s2 - ( y2(2,i) + r2 ) * s2
       end do
       ! End of spline interpolation
       do ios=1,nomega
         beta     = REAL(omegame0i_tmp(ios))
         betasq   = beta*beta
         inv_beta = one/beta
         intsign = sign(half*piinv,beta)
         beta = ABS(beta)
         io = 1; re_intG = zero; re_intK = zero; im_intG = zero; im_intK = zero
         do ii=1,GK_LEVEL
           do
             GKttab = two*alph*xtab(io+1)/(beta-(beta-alph)*xtab(io+1))-one
             if (GKttab > KronN(ii)) EXIT
             io = io + 1
           end do
           temp1     = half*(KronN(ii)+one)
           temp2     = temp1 - half
           temp3     = temp2*temp2
           temp4     = half/(temp3 + quarter)
           ttil      = beta*temp1/(alph-(alph-beta)*temp1)
           tau       = ttil - xtab(io)
           ref       = ftab (io) + tau*(y (1,io)+tau*(y (2,io)+tau*y (3,io)))
           fint (ii) = -ref*(one-ttil)*temp4
           imf       = ftab2(io) + tau*(y2(1,io)+tau*(y2(2,io)+tau*y2(3,io)))
           fint2(ii) = -imf*(one-ttil)*temp4
           re_intK   = KronW(ii)*fint (ii)
           im_intK   = KronW(ii)*fint2(ii)
           ket(spadc+ig,ios) = ket(spadc+ig,ios)+intsign*CMPLX(re_intK,im_intK)
           end do ! ii
       end do !ios
     end do !ig
!!$OMP END PARALLEL DO

   end select !intmethod

   local_one = one
   local_zero = zero

   ! ============================================
   ! ==== Add contribution coming from poles ====
   ! ============================================
   ! First see if the contribution has been checked before the routine is entered
   if (present(calc_poles)) then
     my_calc_poles = calc_poles
   else ! Otherwise check locally if there is a contribution
     do ios=1,nomega
       if (omegame0i_tmp(ios)>tol12) then
         if ((local_one-theta_mu_minus_e0i)<tol12) my_calc_poles(ios) = .FALSE.
       end if
       if (omegame0i_tmp(ios)<-tol12) then
         if (theta_mu_minus_e0i<tol12) my_calc_poles(ios) = .FALSE.
       end if
     end do !ios
   end if

   if (ANY(my_calc_poles(:))) then ! Make sure we only enter if necessary
! *** OPENMP SECTION *** Added by MS
!!OMP !write(std_out,'(a,i0)') ' Entering openmp loop. Number of threads: ',xomp_get_num_threads()
!$OMP PARALLEL SHARED(npwc,nomega,nomegaer,theta_mu_minus_e0i,spadc,local_one,local_zero, &
!$OMP                    omega,epsrho,omegame0i_tmp,ket,my_calc_poles) &
!$OMP PRIVATE(ig,ios,rtmp_r,rtmp_i,tmp_x,tmp_y,rt_real,rt_imag,ct,ierr) REDUCTION(+:my_err)
!!OMP $ write(std_out,'(a,i0)') ' Entering openmp loop. Number of threads: ',xomp_get_num_threads()
!$OMP DO
     do ig=1,npwc
       ! Prepare the spline interpolation by filling at once the arrays rtmp_r, rtmp_i
       call spline(DBLE(omega(1:nomegaer)),DBLE(epsrho(ig,1:nomegaer)),nomegaer,local_zero,local_zero,rtmp_r)
       call spline(DBLE(omega(1:nomegaer)),DBLE(AIMAG(epsrho(ig,1:nomegaer))),nomegaer,local_zero,local_zero,rtmp_i)
       ! call spline_complex( DBLE(omega(1:nomegaer)), epsrho(ig,1:nomegaer), nomegaer, zero, zero, rtmp )

       do ios=1,nomega
         if (.NOT.my_calc_poles(ios)) CYCLE

         ! Interpolate real and imaginary part of epsrho at |omegame0i_tmp|.
         tmp_x(1) = ABS(omegame0i_tmp(ios))
         call splint(nomegaer,DBLE(omega(1:nomegaer)),DBLE(epsrho(ig,1:nomegaer)),rtmp_r,1,tmp_x,tmp_y,ierr=ierr)
         if (ig==1.and.ispinor==1) my_err = my_err + ierr
         rt_real = tmp_y(1)

         tmp_x(1) = ABS(omegame0i_tmp(ios))
         call splint(nomegaer,DBLE(omega(1:nomegaer)),DBLE(AIMAG(epsrho(ig,1:nomegaer))),rtmp_i,1,tmp_x,tmp_y)
         rt_imag = tmp_y(1)
         !!call splint_complex(nomegaer,DBLE(omega(1:nomegaer)),epsrho(ig,1:nomegaer),rtmp,1,tmp_x,yfit)

         ct=DCMPLX(rt_real,rt_imag)

         if (omegame0i_tmp(ios)>tol12) then
           ket(spadc+ig,ios)=ket(spadc+ig,ios)+ct*(local_one-theta_mu_minus_e0i)
         end if
         if (omegame0i_tmp(ios)<-tol12) then
           ket(spadc+ig,ios)=ket(spadc+ig,ios)-ct*theta_mu_minus_e0i
         end if

       end do !ios
     end do !ig
!$OMP END DO
!$OMP END PARALLEL
   end if ! ANY(my_calc_poles)
 end do !ispinor

 npoles_missing = npoles_missing + my_err

 if (INTMETHOD>2) then
   ABI_DEALLOCATE(KronN)
   ABI_DEALLOCATE(KronW)
   ABI_DEALLOCATE(GaussW)
   ABI_DEALLOCATE(fint)
   ABI_DEALLOCATE(fint2)
 end if

end subroutine calc_sigc_cd
!!***

!!****f* ABINIT/calc_sig_ppm_comp
!!
!! NAME
!! calc_sig_ppm_comp
!!
!! FUNCTION
!! Calculating contributions to self-energy operator using a plasmon-pole model
!!
!! INPUTS
!!  nomega=number of frequencies to consider
!!  npwc= number of G vectors in the plasmon pole
!!  npwc1= 1 if ppmodel==3, =npwc if ppmodel== 4, 1 for all the other cases
!!  npwc2= 1 if ppmodel==3, =1    if ppmodel== 4, 1 for all the other cases
!!  npwx=number of G vectors in rhotwgp
!!  ppmodel=plasmon pole model
!!  theta_mu_minus_e0i= $\theta(\mu-\epsilon_{k-q,b1,s}), defines if the state is occupied or not
!!  zcut=small imaginary part to avoid the divergence. (see related input variable)
!!  omegame0i(nomega)=frequencies where evaluate \Sigma_c ($\omega$ - $\epsilon_i$
!!  otq(npwc,npwc2)=plasmon pole parameters for this q-point
!!  botsq(npwc,npwc1)=plasmon pole parameters for this q-point
!!  eig(npwc,npwc)=the eigvectors of the symmetrized inverse dielectric matrix for this q point
!!   (first index for G, second index for bands)
!!  rhotwgp(npwx)=oscillator matrix elements divided by |q+G| i.e
!!    $\frac{\langle b1 k-q s | e^{-i(q+G)r | b2 k s \rangle}{|q+G|}$
!!
!! OUTPUT
!!  sigcme(nomega) (to be described), only relevant if ppm3 or ppm4
!!
!!  ket(npwc,nomega):
!!
!!  In case of ppmodel==1,2 it contains
!!
!!   ket(G,omega) = Sum_G2       conjg(rhotw(G)) * Omega(G,G2) * rhotw(G2)
!!                          ---------------------------------------------------
!!                            2 omegatw(G,G2) (omega-E_i + omegatw(G,G2)(2f-1))
!!
!! NOTES
!! Taken from old routine
!!
!! PARENTS
!!      calc_sigc_me
!!
!! CHILDREN
!!
!! SOURCE

subroutine calc_sig_ppm_comp(npwc,nomega,rhotwgp,botsq,otq,omegame0i_io,zcut,theta_mu_minus_e0i,ket,ppmodel,npwx,npwc1,npwc2)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc,npwc1,npwc2,npwx,ppmodel
 real(dp),intent(in) :: omegame0i_io,theta_mu_minus_e0i,zcut
!arrays
 complex(gwpc),intent(in) :: botsq(npwc,npwc1),rhotwgp(npwx),otq(npwc,npwc2)
 complex(gwpc),intent(inout) :: ket(npwc,nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,io
 real(dp) :: den,otw,twofm1_zcut
 complex(gwpc) :: num,rhotwgdp_igp
 logical :: fully_occupied,totally_empty
 character(len=500) :: msg
!arrays
 complex(gwpc),allocatable :: ket_comp(:)

!*************************************************************************

 if (ppmodel/=1.and.ppmodel/=2) then
   write(msg,'(a,i0,a)')' The completeness trick cannot be used when ppmodel is ',ppmodel,' It should be set to 1 or 2. '
   MSG_ERROR(msg)
 end if

 ABI_ALLOCATE(ket_comp,(npwc))
 ket_comp(:)=0.d0

 fully_occupied=(abs(theta_mu_minus_e0i-1.)<0.001)
 totally_empty=(abs(theta_mu_minus_e0i)<0.001)

 if(.not.(totally_empty)) then ! not totally empty
   twofm1_zcut=zcut
   do igp=1,npwc
     rhotwgdp_igp=rhotwgp(igp)
     do ig=1,npwc
       otw=DBLE(otq(ig,igp)) ! in principle otw -> otw - ieta
       num = botsq(ig,igp)*rhotwgdp_igp

       den = omegame0i_io-otw
       if (den**2>zcut**2) then
         ket_comp(ig) = ket_comp(ig) - num/(den*otw)*theta_mu_minus_e0i
       end if
     end do !ig
   end do !igp
 end if ! not totally empty

 if(.not.(fully_occupied)) then ! not fully occupied
   twofm1_zcut=-zcut

   do igp=1,npwc
     rhotwgdp_igp=rhotwgp(igp)
     do ig=1,npwc
       otw=DBLE(otq(ig,igp)) ! in principle otw -> otw - ieta
       num = botsq(ig,igp)*rhotwgdp_igp

       den = omegame0i_io-otw
       if (den**2>zcut**2) then
         ket_comp(ig) = ket_comp(ig) - num/(den*otw)*(1.-theta_mu_minus_e0i)
       end if
     end do !ig
   end do !igp
 end if ! not fully occupied

 do io=1,nomega
   ket(:,io)=ket(:,io)+0.5*ket_comp(:)
 end do

 ABI_DEALLOCATE(ket_comp)

end subroutine calc_sig_ppm_comp
!!***

end module m_sigc
!!***
