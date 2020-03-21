!!****m* ABINIT/m_cohsex
!! NAME
!!  m_cohsex
!!
!! FUNCTION
!! Calculate diagonal and off-diagonal matrix elements of the SEX or COHSEX self-energy operator.
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2020 ABINIT group (FB, GMR, VO, LR, RWG, MG, RShaltaf)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_cohsex

 use defs_basis
 use m_defs_ptgroups
 use m_gwdefs
 use m_xmpi
 use m_errors
 use m_abicore

 use defs_datatypes,  only : pseudopotential_type, ebands_t
 use m_time,          only : timab
 use m_fstrings,      only : sjoin, itoa
 use m_hide_blas,     only : xdotc, xgemv
 use m_numeric_tools, only : hermitianize, imin_loc
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_t
 use m_bz_mesh,       only : kmesh_t, get_BZ_item, findqg0, littlegroup_t, littlegroup_print
 use m_gsphere,       only : gsphere_t, gsph_fft_tabs
 use m_fft_mesh,      only : get_gftt, rotate_fft_mesh, cigfft
 use m_vcoul_dt
 use m_pawpwij,       only : pawpwff_t, pawpwij_t, pawpwij_init, pawpwij_free, paw_rho_tw_g
 use m_wfd,           only : wfd_t
 use m_oscillators,   only : rho_tw_g, calc_wfwfg
 use m_screening,     only : epsm1_symmetrizer, get_epsm1, epsilonm1_results
 use m_esymm,         only : esymm_t, esymm_symmetrize_mels, esymm_failed
 use m_sigma,         only : sigma_t, sigma_distribute_bks
 use m_pawang,        only : pawang_type
 use m_pawtab,        only : pawtab_type
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy, paw_overlap
 use m_paw_sym,       only : paw_symcprj

 implicit none

 private
!!***

 public :: cohsex_me
!!***

contains
!!***

!!****f* ABINIT/cohsex_me
!! NAME
!! cohsex_me
!!
!! FUNCTION
!! Calculate diagonal and off-diagonal matrix elements of the SEX or COHSEX self-energy operator.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (FB, GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! sigmak_ibz=Index of the k-point in the IBZ.
!! minbnd, maxbnd= min and Max band index for GW correction (for this k-point)
!! iomode=Option defining the file format of the SCR file (Fortran, NETCDF)
!! Er <Epsilonm1_results> (see the definition of this structured datatype)
!!    %mqmem=if 0 use out-of-core method in which a single q-slice of espilon is read inside the loop over k
!!    %nomega_i=Number of imaginary frequencies.
!!    %nomega_r=Number of real frequencies.
!!    %nomega=Total number of frequencies.
!! Gsph_c<gsphere_t>= info on the G-sphere for Sigma_x
!!    %nsym=number of symmetry operations
!!    %rottb(ng,timrev,nsym)=index of (IS) G where I is the identity or the inversion
!!      operation and G is one of the ng vectors in reciprocal space
!!    %timrev=2 if time-reversal symmetry is used, 1 otherwise
!!    %gvec(3,Sigp%npwc)=integer coordinates of each plane wave in reciprocal space
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
!!  much slower but it requires less memory
!! QP_BSt<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!!  Paw_pwff<pawpwff_t>=Form factor used to calculate the onsite mat. elements of a plane wave.
!! allQP_sym(%nkibz,%nsppol)<esymm_t>=Datatype collecting data on the irreducible representaions of the
!!  little group of kcalc in the KS representation as well as the symmetry of the bdgw_k states.
!!  Sr=sigma_t (see the definition of this structured datatype)
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
!!      calc_coh,calc_wfwfg,epsm1_symmetrizer,esymm_symmetrize_mels,findqg0
!!      get_bz_item,get_epsm1,get_gftt,gsph_fft_tabs,hermitianize
!!      littlegroup_print,paw_rho_tw_g,paw_symcprj,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawpwij_free,pawpwij_init,rho_tw_g,rotate_fft_mesh
!!      sigma_distribute_bks,timab,wfd_change_ngfft,wfd_get_cprj
!!      wfd_get_many_ur,wfd_get_ur,wrtout,xgemv,xmpi_sum
!!
!! SOURCE

subroutine cohsex_me(sigmak_ibz,ikcalc,nomega_sigc,minbnd,maxbnd,Cryst,QP_BSt,Sigp,Sr,Er,Gsph_c,Vcp,&
& Kmesh,Qmesh,Ltg_k,Pawtab,Pawang,Paw_pwff,Psps,Wfd,allQP_sym,gwc_ngfft,iomode,prtvol,sigcme_tmp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: sigmak_ibz,ikcalc,prtvol,iomode,nomega_sigc,minbnd,maxbnd
 type(crystal_t),intent(in) :: Cryst
 type(ebands_t),target,intent(in) :: QP_BSt
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(vcoul_t),intent(in) :: Vcp
 type(Epsilonm1_results),intent(inout) :: Er
 type(gsphere_t),intent(in) :: Gsph_c
 type(littlegroup_t),intent(in) :: Ltg_k
 type(Pseudopotential_type),intent(in) :: Psps
 type(pawang_type),intent(in) :: pawang
 type(sigparams_t),target,intent(in) :: Sigp
 type(sigma_t),intent(in) :: Sr
 type(wfd_t),target,intent(inout) :: Wfd
!arrays
 integer,intent(in) :: gwc_ngfft(18)
 complex(dpc),intent(out) :: sigcme_tmp(nomega_sigc,minbnd:maxbnd,minbnd:maxbnd,Wfd%nsppol*Sigp%nsig_ab)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
 type(esymm_t),target,intent(in) :: allQP_sym(Wfd%nkibz,Wfd%nsppol)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp=2,ndat1=1
 integer :: iab,ib,ib1,ib2,ierr,ig,ii,iik,itim_q,i1,i2,npwc
 integer :: ik_bz,ik_ibz,io,isym_q,iq_bz,iq_ibz,spin,isym,jb,is_idx
 integer :: band,band1,band2,idle,rank
 integer :: jik,jk_bz,jk_ibz,kb,nspinor,nsppol
 integer :: nomega_tot,nq_summed,ispinor,ibsp,dimcprj_gw
 integer :: spad,spadc,spadc1,spadc2,irow,my_nbks
 integer :: ndegs,wtqm,wtqp,mod10
 integer :: isym_kgw,isym_ki,gwc_mgfft,use_padfft,gwc_fftalga,gwc_nfftot,ifft,npw_k
 real(dp) :: fact_sp,theta_mu_minus_e0i,tol_empty,gw_gsq
 complex(dpc) :: ctmp,ph_mkgwt,ph_mkt
 logical :: iscompatibleFFT,q_is_gamma
 character(len=500) :: msg
!arrays
 integer :: g0(3),spinor_padc(2,4),nbv_ks(Kmesh%nibz,Wfd%nsppol)
 integer,allocatable :: proc_distrb(:,:,:),coh_distrb(:,:,:,:),degtab(:,:,:)
 integer,allocatable :: igfftcg0(:),gw_gfft(:,:),gw_gbound(:,:),irottb(:,:),ktabr(:,:)
 integer :: got(Wfd%nproc)
 real(dp) :: ksum(3),kgw(3),kgw_m_ksum(3),q0(3),tsec(2),qbz(3),spinrot_kbz(4),spinrot_kgw(4)
 real(dp),ABI_CONTIGUOUS pointer :: qp_ene(:,:,:),qp_occ(:,:,:)
 complex(gwpc) :: sigcohme(Sigp%nsig_ab)
 complex(gwpc),allocatable :: vc_sqrt_qbz(:),rhotwg(:),rhotwgp(:),sigsex(:)
 complex(gwpc),allocatable :: epsm1_qbz(:,:,:)
 complex(gwpc),allocatable :: sigc_ket(:,:)
 complex(gwpc),allocatable :: rhotwg_ki(:,:)
 complex(gwpc),allocatable :: sigctmp(:,:)
 complex(gwpc),allocatable :: wfr_bdgw(:,:),ur_sum(:),wf1swf2_g(:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: cg_jb(:),cg_sum(:)
 complex(dpc),allocatable :: sym_cme(:,:,:,:),sigc(:,:,:,:,:)
 logical :: rank_mask(Wfd%nproc),can_symmetrize(Wfd%nsppol)
 logical,allocatable :: bks_mask(:,:,:)
 type(sigijtab_t),pointer :: Sigcij_tab(:)
 type(pawcprj_type),allocatable :: Cprj_kgw(:,:),Cprj_ksum(:,:)
 type(pawpwij_t),allocatable :: Pwij_qg(:),Pwij_fft(:)
 type(esymm_t),pointer :: QP_sym(:)

!************************************************************************

 DBG_ENTER("COLL")

 call timab(423,1,tsec) ! cohsex_me

 ! Initial check
 ABI_CHECK(Sr%nomega_r == Sigp%nomegasr,"")
 ABI_CHECK(Sr%nomega4sd == Sigp%nomegasrd,"")
 !ABI_CHECK(Sigp%npwc==Gsph_c%ng,"")

 ! Initialize some values
 nspinor = Wfd%nspinor; nsppol = Wfd%nsppol
 npwc = sigp%npwc
 spinor_padc = RESHAPE([0, 0, npwc, npwc, 0, npwc, npwc,0], [2, 4])

 qp_ene => QP_BSt%eig; qp_occ => QP_BSt%occ

 ! Extract the symmetries of the bands for this k-point
 QP_sym => allQP_sym(sigmak_ibz,1:nsppol)

 ! Index of the GW point in the BZ array, its image in IBZ and time-reversal ===
 jk_bz=Sigp%kptgw2bz(ikcalc)
 call get_BZ_item(Kmesh,jk_bz,kgw,jk_ibz,isym_kgw,jik,ph_mkgwt)
 !%call get_IBZ_item(Kmesh,jk_ibz,kibz,wtk)
 spinrot_kgw=Cryst%spinrot(:,isym_kgw)
 ib1 = minbnd; ib2 = maxbnd

 write(msg,'(2a,3f8.3,2a,2(i3,a))')ch10,&
&  ' Calculating <nk|Sigma_c(omega)|nk> at k = ',kgw(:),ch10,&
&  ' bands n = from ',ib1,' to ',ib2,ch10
 call wrtout(std_out,msg,'COLL')

 if (ANY(gwc_ngfft(1:3) /= Wfd%ngfft(1:3)) ) call wfd%change_ngfft(Cryst,Psps,gwc_ngfft)
 gwc_mgfft = MAXVAL(gwc_ngfft(1:3))
 gwc_fftalga = gwc_ngfft(7)/100 !; gwc_fftalgc=MOD(gwc_ngfft(7),10)

 can_symmetrize = .FALSE.
 if (Sigp%symsigma>0) then
   can_symmetrize = .TRUE.
   if (Sigp%gwcalctyp >= 20) then
    do spin=1,Wfd%nsppol
      can_symmetrize(spin) = .not. esymm_failed(QP_sym(spin))
      if (.not.can_symmetrize(spin)) then
        write(msg,'(a,i0,4a)')" Symmetrization cannot be performed for spin: ",spin,ch10,&
&         " band classification encountered the following problem: ",ch10,TRIM(QP_sym(spin)%err_msg)
        MSG_WARNING(msg)
      end if
    end do
   end if
   if (nspinor == 2) MSG_WARNING('Symmetrization with nspinor=2 not implemented')
 end if

 mod10=MOD(Sigp%gwcalctyp, 10)

 call timab(491,1,tsec) ! csigme(tot) Overall clock. TODO check this
 call timab(495,1,tsec) ! csigme (SigC)

 ! Normalization of theta_mu_minus_e0i
 ! If nsppol==2, qp_occ $\in [0,1]$
 SELECT CASE (nsppol)
 CASE (1)
   fact_sp=half; tol_empty=0.01   ! below this value the state is assumed empty
   if (nspinor==2) then
    fact_sp=one; tol_empty=0.005  ! below this value the state is assumed empty
   end if
 CASE (2)
   fact_sp=one; tol_empty=0.005   ! to be consistent and obtain similar results if a metallic
 CASE DEFAULT                     ! spin unpolarized system is treated using nsppol==2
   MSG_BUG('Wrong nsppol')
 END SELECT

 call timab(442,1,tsec) ! csigme(init0)

 ! Precalculate the FFT index of $(R^{-1}(r-\tau))$ ===
 ! S=\transpose R^{-1} and k_BZ = S k_IBZ
 ! irottb is the FFT index of $R^{-1} (r-\tau)$ used to symmetrize u_Sk.
 gwc_nfftot = PRODUCT(gwc_ngfft(1:3))
 ABI_MALLOC(irottb,(gwc_nfftot,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,gwc_ngfft,irottb,iscompatibleFFT)
 if (.not.iscompatibleFFT) then
   MSG_WARNING("FFT mesh is not compatible with symmetries. Results might be affected by large errors!")
 end if

 ABI_MALLOC(ktabr,(gwc_nfftot, Kmesh%nbz))
 do ik_bz=1,Kmesh%nbz
   isym=Kmesh%tabo(ik_bz)
   do ifft=1,gwc_nfftot
     ktabr(ifft,ik_bz)=irottb(ifft,isym)
   end do
 end do
 ABI_FREE(irottb)

 ! The number of occupied states for each point in the IBZ and spin.
 do spin=1,nsppol
   do ik_ibz=1,Kmesh%nibz
     nbv_ks(ik_ibz,spin) = COUNT(qp_occ(:,ik_ibz,spin)>=tol_empty)
   end do
 end do

 ! (b,k,s) mask for MPI distribution of the sum over occupied states in the BZ.
 ABI_MALLOC(bks_mask,(Wfd%mband,Kmesh%nbz,nsppol))
 bks_mask=.FALSE.
 do spin=1,nsppol
   do ik_bz=1,Kmesh%nbz
      ik_ibz = Kmesh%tab(ik_bz)
      bks_mask(1:nbv_ks(ik_ibz,spin),ik_bz,spin) = .TRUE.
   end do
 end do

 ! Distribute the individual terms of the sum over the BZ taking into account symmetries and MPI memory distribution.
 ! got is used to optimize the distribution if more than one node can calculate the same (b,k,s) element.
 got=0
 ABI_MALLOC(proc_distrb,(Wfd%mband,Kmesh%nbz,nsppol))
 call sigma_distribute_bks(Wfd,Kmesh,Ltg_k,Qmesh,nsppol,can_symmetrize,kgw,Sigp%mg0,my_nbks,&
&  proc_distrb,got,bks_mask,global=.TRUE.)

 ABI_FREE(bks_mask)

 write(msg,'(a,i0,a)')" Will sum ",my_nbks," (b,k,s) occupied states in (COHSEX|SEX)."
 call wrtout(std_out,msg,'PERS')

 Sigcij_tab => Sigp%Sigcij_tab(ikcalc,1:nsppol)

 if (mod10==SIG_COHSEX) then
   ! Distribute the COHSEX terms, taking into account the symmetries of the Sigma_ij matrix.
   ABI_MALLOC(coh_distrb,(ib1:ib2,ib1:ib2,Kmesh%nbz,nsppol))

   coh_distrb = xmpi_undefined_rank
   do spin=1,nsppol
     do ik_bz=1,Kmesh%nbz
        if (ANY(proc_distrb(:,ik_bz,spin) /= xmpi_undefined_rank) ) then ! This BZ point will be calculated.
           rank_mask = .FALSE. ! To select only those nodes that will treat (k,s).
           do band=1,Wfd%mband
             rank = proc_distrb(band,ik_bz,spin)
             if (rank /= xmpi_undefined_rank) rank_mask(rank+1)=.TRUE.
           end do
           do band2=ib1,ib2
             do irow=1,Sigcij_tab(spin)%col(band2)%size1   ! Looping over the upper triangle of sigma_ij with non-zero elements.
               band1 = Sigcij_tab(spin)%col(band2)%bidx(irow)
               idle = imin_loc(got,mask=rank_mask)
               got(idle) = got(idle)+1
               coh_distrb(band1,band2,ik_bz,spin) = idle-1
             end do
           end do
        end if
     end do
   end do

   write(msg,'(a,i0,a)')" will treat ",COUNT(coh_distrb==Wfd%my_rank)," COH terms."
   call wrtout(std_out,msg,'PERS')
 end if

 ABI_MALLOC(rhotwg_ki, (npwc * nspinor, minbnd:maxbnd))
 rhotwg_ki=czero_gw
 ABI_MALLOC(rhotwg, (npwc * nspinor))
 ABI_MALLOC(rhotwgp  ,(npwc * nspinor))
 ABI_MALLOC(vc_sqrt_qbz, (npwc))

 ! Additional allocations for PAW
 if (Psps%usepaw==1) then
   ABI_DT_MALLOC(Cprj_ksum,(Cryst%natom,nspinor))
   call pawcprj_alloc(Cprj_ksum,0,Wfd%nlmn_atm)

   ! For COHSEX we need the onsite terms of the PW on the FFT mesh.
   ! gw_gfft is the set of plane waves in the FFT Box for the oscillators.
   if (mod10==SIG_COHSEX) then
     ABI_MALLOC(gw_gfft,(3,gwc_nfftot))
     q0=zero
     call get_gftt(gwc_ngfft,q0,Cryst%gmet,gw_gsq,gw_gfft)
     ABI_DT_MALLOC(Pwij_fft,(Psps%ntypat))
     call pawpwij_init(Pwij_fft,gwc_nfftot,(/zero,zero,zero/),gw_gfft,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
   end if
 end if ! usepaw==1

 ! === Calculate total number of frequencies and allocate related arrays ===
 ! sigcme2 is used to accumulate the diagonal matrix elements over k-points and
 ! GW bands, used only in case of ppmodel 3 and 4 (TODO save memory)
 nomega_tot=Sr%nomega_r+Sr%nomega4sd

 ABI_MALLOC(sigctmp, (nomega_sigc,Sigp%nsig_ab))
 sigctmp = czero_gw
 ABI_MALLOC(sigc_ket, (npwc*nspinor, nomega_sigc))

 if (mod10==SIG_COHSEX)  then
   ABI_MALLOC(wf1swf2_g,(gwc_nfftot*nspinor))
 end if

 ! Arrays storing the contribution given by the Hermitian/anti-Hermitian part of \Sigma_c
 !allocate(aherm_sigc_ket(npwc*nspinor,nomega_sigc))
 !allocate( herm_sigc_ket(npwc*nspinor,nomega_sigc))
 ABI_MALLOC(sigsex,(npwc))
 sigcme_tmp=czero

 ABI_MALLOC(sigc,(2,nomega_sigc,ib1:ib2,ib1:ib2,nsppol*Sigp%nsig_ab))
 sigc=czero

 ! Here we divide the states where the QP energies are required into complexes. Note however that this approach is not
 ! based on group theory, and it might lead to spurious results in case of accidental degeneracies.
 nq_summed=Kmesh%nbz
 if (Sigp%symsigma>0) then
   call littlegroup_print(Ltg_k,std_out,prtvol,'COLL')
   nq_summed=SUM(Ltg_k%ibzq(:))
   !
   ! Find number of degenerate states and number of bands in each subspace
   ! The tolerance is a little bit arbitrary (0.001 eV)
   ! It could be reduced, in particular in case of nearly accidental degeneracies
   ABI_MALLOC(degtab,(ib1:ib2,ib1:ib2,nsppol))
   degtab=0
   do spin=1,nsppol
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

 ! TODO if single q (ex molecule) dont allocate epsm1q, avoid waste of memory
 ABI_MALLOC_OR_DIE(epsm1_qbz, (npwc, npwc, 1), ierr)
 ABI_MALLOC(igfftcg0,(Gsph_c%ng))

 ! Out-of-core solution for epsilon.
 if (Er%mqmem==0) then
   MSG_COMMENT('Reading q-slices from file. Slower but less memory.')
 end if

 call timab(442,2,tsec)

 ! ==========================================
 ! ==== Fat loop over k_i in the full BZ ====
 ! ==========================================
 ABI_MALLOC(ur_sum,(gwc_nfftot*nspinor))

 do spin=1,nsppol
   if (ALL(proc_distrb(:,:,spin)/=Wfd%my_rank)) CYCLE

   ABI_MALLOC(wfr_bdgw,(gwc_nfftot*nspinor,ib1:ib2))
   call wfd%get_many_ur([(jb, jb=ib1,ib2)], jk_ibz, spin, wfr_bdgw)

   if (Wfd%usepaw==1) then
     ! Load cprj for GW states, note the indexing.
     dimcprj_gw=nspinor*(ib2-ib1+1)
     ABI_DT_MALLOC(Cprj_kgw,(Cryst%natom,ib1:ib1+dimcprj_gw-1))
     call pawcprj_alloc(Cprj_kgw,0,Wfd%nlmn_atm)
     ibsp=ib1
     do jb=ib1,ib2
       call wfd%get_cprj(jb,jk_ibz,spin,Cryst,Cprj_ksum,sorted=.FALSE.)
       call paw_symcprj(jk_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,Cprj_ksum)
       call pawcprj_copy(Cprj_ksum,Cprj_kgw(:,ibsp:ibsp+(nspinor-1)))
       ibsp=ibsp+nspinor
     end do
   end if

   do ik_bz=1,Kmesh%nbz
     ! Parallelization over k-points and spin
     ! For the spin there is another check in the inner loop.
     if (ALL(proc_distrb(:,ik_bz,spin)/=Wfd%my_rank)) CYCLE

     call timab(443,1,tsec) ! csigme (initq)

     ! Find the corresponding irreducible k-point
     call get_BZ_item(Kmesh,ik_bz,ksum,ik_ibz,isym_ki,iik,ph_mkt)
     spinrot_kbz(:)=Cryst%spinrot(:,isym_ki)

     ! Identify q and G0 where q+G0=k_GW-k_i
     kgw_m_ksum=kgw-ksum
     call findqg0(iq_bz,g0,kgw_m_ksum,Qmesh%nbz,Qmesh%bz,Sigp%mG0)

     ! Symmetrize the matrix elements.
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

     !%write(msg,'(2(a,i4),a,i3)')' csigme : ik_bz ',ik_bz,'/',Kmesh%nbz,' done by processor ',Wfd%my_rank
     !%call wrtout(std_out,msg,'PERS')

     ! Find the corresponding irred q-point.
     call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)
     q_is_gamma = (normv(qbz,Cryst%gmet,"G") < GW_TOL_W0)

     ! Tables for the FFT of the oscillators.
     !  a) FFT index of the G-G0.
     !  b) gw_gbound table for the zero-padded FFT performed in rhotwg.
     ABI_MALLOC(gw_gbound,(2*gwc_mgfft+8,2))
     call gsph_fft_tabs(Gsph_c,g0,gwc_mgfft,gwc_ngfft,use_padfft,gw_gbound,igfftcg0)
     if ( ANY(gwc_fftalga == [2, 4]) ) use_padfft=0 ! Pad-FFT is not coded in rho_tw_g
#ifdef FC_IBM
 ! XLF does not deserve this optimization (problem with [v67mbpt][t03])
 use_padfft = 0
#endif
     if (use_padfft==0) then
       ABI_FREE(gw_gbound)
       ABI_MALLOC(gw_gbound,(2*gwc_mgfft+8,2*use_padfft))
     end if

     if (Psps%usepaw==1) then
       ! Get PAW oscillator matrix elements $ <phj/r|e^{-i(q+G)}|phi/r> - <tphj/r|e^{-i(q+G)}|tphi/r> $ in packed form.
       ABI_DT_MALLOC(Pwij_qg,(Psps%ntypat))
       q0 = qbz !;if (q_is_gamma) q0 = (/0.00001_dp,0.00001_dp,0.00001_dp/) ! GW_Q0_DEFAULT
       call pawpwij_init(Pwij_qg, npwc, q0, Gsph_c%gvec, Cryst%rprimd, Psps, Pawtab, Paw_pwff)
     end if

     if (Er%mqmem==0) then
       ! Read q-slice of epsilon^{-1}|chi0 in Er%epsm1(:,:,:,1) (much slower but less memory).
       call get_epsm1(Er,Vcp,0,0,iomode,xmpi_comm_self,iqibzA=iq_ibz)
     end if

     ! Only omega==0 for SEX or COHSEX
     call Epsm1_symmetrizer(iq_bz, 1, npwc, Er, Gsph_c, Qmesh, .True., epsm1_qbz)

     ! Get Fourier components of the Coulomb interaction in the BZ.
     ! In 3D systems, neglecting umklapp,  vc(Sq,sG)=vc(q,G)=4pi/|q+G|
     ! The same relation holds for 0-D systems, but not in 1-D or 2D systems. It depends on S.
     do ig=1,npwc
       vc_sqrt_qbz(Gsph_c%rottb(ig,itim_q,isym_q)) = Vcp%vc_sqrt(ig,iq_ibz)
     end do

     call timab(443,2,tsec) ! csigme (initq)

     ! Sum over bands.
     do ib=1,Sigp%nbnds
       ! Parallelism over spin
       ! This processor has this k-point but what about spin?
       if (proc_distrb(ib,ik_bz,spin)/=Wfd%my_rank) CYCLE

       ! Skip empty state ib for HF, SEX, and COHSEX.
       if (qp_occ(ib,ik_ibz,spin)<tol_empty) CYCLE

       theta_mu_minus_e0i=fact_sp*qp_occ(ib,ik_ibz,spin)

       call wfd%get_ur(ib,ik_ibz,spin,ur_sum)

       if (Psps%usepaw==1) then
         ! Load cprj for point ksum, this spin or spinor and *THIS* band.
         ! TODO MG I could avoid doing this but I have to exchange spin and bands ???
         ! For sure there is a better way to do this!
         call wfd%get_cprj(ib,ik_ibz,spin,Cryst,Cprj_ksum,sorted=.FALSE.)
         call paw_symcprj(ik_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,Cprj_ksum)
       end if

       do jb=ib1,ib2
         ! Get all <k-q,ib,s|e^{-i(q+G).r}|s,jb,k>, at once.
         call rho_tw_g(nspinor,npwc,gwc_nfftot,ndat1,gwc_ngfft,1,use_padfft,igfftcg0,gw_gbound,&
&          ur_sum       ,iik,ktabr(:,ik_bz),ph_mkt  ,spinrot_kbz,  &
&          wfr_bdgw(:,jb),jik,ktabr(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&          nspinor,rhotwg_ki(:,jb))

         if (Psps%usepaw==1) then
           ! Add on-site contribution, projectors are already in BZ !TODO Recheck this!
           i2=jb; if (nspinor==2) i2=(2*jb-1)
           spad=(nspinor-1)
           call paw_rho_tw_g(npwc,nspinor,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_c%gvec,&
&            Cprj_ksum(:,:),Cprj_kgw(:,i2:i2+spad),Pwij_qg,rhotwg_ki(:,jb))
         end if

         ! Multiply by the square root of the Coulomb term.
         ! In 3-D systems, the factor sqrt(4pi) is included)
         do ii=1,nspinor
           spad = (ii-1) * npwc
           rhotwg_ki(spad+1:spad+npwc,jb) = rhotwg_ki(spad+1:spad+npwc,jb)*vc_sqrt_qbz(1:npwc)
         end do

         ! === Treat analytically the case q --> 0 ===
         ! * The oscillator is evaluated at q=O as it is considered constant in the small cube around Gamma
         !   while the Colulomb term is integrated out
         ! * In the scalar case we have nonzero contribution only if ib==jb
         ! * For nspinor==2 evalute <ib,up|jb,up> and <ib,dwn|jb,dwn>,
         !   impose orthonormalization since npwwfn might be < npwvec.
         if (ik_bz == jk_bz) then
           if (nspinor == 1) then
             rhotwg_ki(1, jb) = czero_gw
             if (ib==jb) rhotwg_ki(1, jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp)
           else
             npw_k = Wfd%npwarr(ik_ibz)
             rhotwg_ki(1, jb) = zero; rhotwg_ki(npwc+1, jb) = zero
             if (ib == jb) then
               cg_sum => Wfd%Wave(ib, ik_ibz, spin)%ug
               cg_jb  => Wfd%Wave(jb, jk_ibz, spin)%ug
               ctmp = xdotc(npw_k, cg_sum(1:), 1, cg_jb(1:), 1)
               rhotwg_ki(1,jb) = CMPLX(SQRT(Vcp%i_sz),0.0_gwp) * real(ctmp)
               ctmp = xdotc(npw_k, cg_sum(npw_k+1:), 1, cg_jb(npw_k+1:), 1)
               rhotwg_ki(npwc+1,jb) = CMPLX(SQRT(Vcp%i_sz),0.0_gwp) * real(ctmp)
               ! PAW is missing

               !rhotwg_ki(1,jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp) * sqrt(half)
               !rhotwg_ki(npwc+1,jb) = CMPLX(SQRT(Vcp%i_sz),0.0_gwp) * sqrt(half)
             end if
           end if
         end if
       end do !jb  Got all matrix elements from minbnd up to maxbnd.

       do kb=ib1,ib2
         ! Get the ket \Sigma|\phi_{k,kb}> according to the method.
         rhotwgp(:) = rhotwg_ki(:,kb)
         sigc_ket = czero_gw

         ! SEX part. TODO add check on theta_mu_minus_e0i
         do ispinor=1,nspinor
           spadc = (ispinor-1) * npwc
           call XGEMV('N',npwc,npwc,cone_gw,epsm1_qbz(:,:,1),npwc,rhotwgp(1+spadc:),1,czero_gw,sigsex,1)

           sigsex(:)= -theta_mu_minus_e0i*sigsex(:)

           do io=1,nomega_tot ! nomega==1 as SEX is energy independent.
             sigc_ket(spadc+1:spadc+npwc,io) = sigsex(:)
           end do
         end do

         ! Loop over the non-zero row elements of this column.
         ! 1) If gwcalctyp<20 : only diagonal elements since QP==KS.
         ! 2) If gwcalctyp>=20:
         !     * Only off-diagonal elements connecting states with same character.
         !     * Only the upper triangle if HF, SEX, or COHSEX.
         do irow=1,Sigcij_tab(spin)%col(kb)%size1
           jb = Sigcij_tab(spin)%col(kb)%bidx(irow)
           rhotwg = rhotwg_ki(:,jb)

           ! Calculate <\phi_j|\Sigma_c|\phi_k>
           ! Different freqs according to method (AC or Perturbative), see nomega_sigc.
           do iab=1,Sigp%nsig_ab
             spadc1=spinor_padc(1,iab); spadc2=spinor_padc(2,iab)
             do io=1,nomega_sigc
               sigctmp(io,iab) = XDOTC(npwc,rhotwg(spadc1+1:),1,sigc_ket(spadc2+1:,io),1)
             end do
           end do

           ! TODO: save wf1swf2_g to avoid having to recalculate it at each q-point.
           if (mod10==SIG_COHSEX) then
             ! Evaluate Static COH. TODO add spinor.
             if (coh_distrb(jb,kb,ik_bz,spin) == Wfd%my_rank) then
               ! COH term is done only once for each k-point.
               ! It does not depend on the index ib summed over.
               coh_distrb(jb,kb,ik_bz,spin) = xmpi_undefined_rank

#if 1
               call calc_wfwfg(ktabr(:,jk_ibz),jik, spinrot_kgw, & ! TODO why jk_ibz?
&                gwc_nfftot,nspinor,gwc_ngfft,wfr_bdgw(:,jb),wfr_bdgw(:,kb),wf1swf2_g)
#else
               ABI_CHECK(jik==1,"jik")
               call calc_wfwfg(ktabr(:,jk_bz),jik, spinrot_kgw, &
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
               end if

               call calc_coh(nspinor,Sigp%nsig_ab,gwc_nfftot,gwc_ngfft,npwc,Gsph_c%gvec,wf1swf2_g,epsm1_qbz(:,:,1),&
&                vc_sqrt_qbz,Vcp%i_sz,iq_ibz,(jb==kb),sigcohme)

               do io=1,nomega_sigc ! Should be 1
                 sigctmp(io,:) = sigctmp(io,:)+sigcohme(:)
               end do

             end if
           end if ! COHSEX

           ! Accumulate and, in case, symmetrize matrix elements of Sigma_c.
           do iab=1,Sigp%nsig_ab
             is_idx = spin; if (nspinor==2) is_idx=iab

             sigcme_tmp(:,jb,kb,is_idx)=sigcme_tmp(:,jb,kb,is_idx) + &
&              (wtqp+wtqm)*DBLE(sigctmp(:,iab)) + (wtqp-wtqm)*j_gw*AIMAG(sigctmp(:,iab))

             sigc(1,:,jb,kb,is_idx)=sigc(1,:,jb,kb,is_idx) + wtqp*      sigctmp(:,iab)
             sigc(2,:,jb,kb,is_idx)=sigc(2,:,jb,kb,is_idx) + wtqm*CONJG(sigctmp(:,iab))
             ! TODO this should be the contribution coming from the anti-hermitian part.
           end do
         end do !jb used to calculate matrix elements of $\Sigma$

       end do !kb to calculate matrix elements of $\Sigma$
     end do !ib

     ! Deallocate k-dependent quantities.
     ABI_FREE(gw_gbound)
     if (Psps%usepaw==1) then
       call pawpwij_free(Pwij_qg)
       ABI_DT_FREE(Pwij_qg)
     end if
   end do !ik_bz

   ABI_FREE(wfr_bdgw)
   if (Wfd%usepaw==1) then
     call pawcprj_free(Cprj_kgw)
     ABI_DT_FREE(Cprj_kgw)
   end if
 end do !spin

 ABI_FREE(igfftcg0)

 ! Gather contributions from all the CPUs.
 call xmpi_sum(sigcme_tmp, wfd%comm, ierr)
 call xmpi_sum(sigc, wfd%comm, ierr)

 ! Multiply by constants
 ! For 3D systems sqrt(4pi) is included in vc_sqrt_qbz.
 sigcme_tmp = sigcme_tmp /(Cryst%ucvol*Kmesh%nbz)
 sigc       = sigc       /(Cryst%ucvol*Kmesh%nbz)

 ! If we have summed over the IBZ_q now we have to average over degenerate states.
 ! Presently only diagonal terms are considered
 ! TODO it does not work if nspinor==2.
 do spin=1,nsppol
   if (can_symmetrize(spin)) then
     ABI_MALLOC(sym_cme, (nomega_tot, ib1:ib2, ib1:ib2, sigp%nsig_ab))
     sym_cme=czero

     ! Average over degenerate diagonal elements.
     ! NOTE: frequencies for \Sigma_c(\omega) should be equal to avoid spurious results.
     ! another good reason to use a strict criterion for the tolerance on eigenvalues.
     do ib=ib1,ib2
       ndegs=0
       do jb=ib1,ib2
         if (degtab(ib,jb,spin)==1) then
           if (nspinor == 1) then
             sym_cme(:, ib, ib, 1) = sym_cme(:, ib, ib, 1) + SUM(sigc(:, :, jb, jb, spin), dim=1)
           else
             do ii=1,sigp%nsig_ab
               sym_cme(:, ib, ib, ii) = sym_cme(:, ib, ib, ii) + SUM(sigc(:, :, jb, jb, ii), dim=1)
             end do
           end if

         end if
         ndegs = ndegs + degtab(ib,jb,spin)
       end do
       sym_cme(:,ib,ib,:) = sym_cme(:,ib,ib,:) / ndegs
     end do

     if (Sigp%gwcalctyp >= 20) then
       call esymm_symmetrize_mels(QP_sym(spin),ib1,ib2,sigc(:,1,:,:,spin),sym_cme(1,:,:,1))
     end if

     ! Copy symmetrized values.
     do ib=ib1,ib2
       do jb=ib1,ib2
         if (nspinor == 1) then
           sigcme_tmp(:,ib,jb,spin) = sym_cme(:,ib,jb,1)
         else
           sigcme_tmp(:,ib,jb,:) = sym_cme(:,ib,jb,:)
         end if
       end do
     end do
     ABI_FREE(sym_cme)
   end if
 end do

 ! Reconstruct the full sigma matrix from the upper triangle (only for HF, SEX and COHSEX)
 if (Sigp%gwcalctyp>=20 .and. sigma_is_herm(Sigp) ) then
   ABI_CHECK(nspinor==1,"cannot hermitianize non-collinear sigma!")
   do spin=1,nsppol
     do io=1,nomega_sigc
       call hermitianize(sigcme_tmp(io,:,:,spin),"Upper")
     end do
   end do
 end if

 ! ===========================
 ! ==== Deallocate memory ====
 ! ===========================
 if (Psps%usepaw==1) then
   if (allocated(gw_gfft)) then
     ABI_FREE(gw_gfft)
   end if
   call pawcprj_free(Cprj_ksum)
   ABI_DT_FREE(Cprj_ksum)
   if (allocated(Pwij_fft)) then
     call pawpwij_free(Pwij_fft)
     ABI_DT_FREE(Pwij_fft)
   end if
 end if

 ABI_FREE(ktabr)
 ABI_FREE(ur_sum)
 ABI_FREE(rhotwg_ki)
 ABI_FREE(rhotwg)
 ABI_FREE(rhotwgp)
 ABI_FREE(vc_sqrt_qbz)
 ABI_FREE(sigc_ket)
 ABI_FREE(epsm1_qbz)
 ABI_FREE(sigctmp)
 ABI_FREE(sigc)
 ABI_FREE(sigsex)
 ABI_FREE(proc_distrb)
 if (mod10==SIG_COHSEX) then
   ABI_FREE(wf1swf2_g)
   ABI_FREE(coh_distrb)
 end if

 if (allocated(degtab)) then
   ABI_FREE(degtab)
 end if

 call timab(495,2,tsec) ! csigme(SigC)
 call timab(491,2,tsec)
 call timab(423,2,tsec) ! cohsex_me

 DBG_EXIT("COLL")

end subroutine cohsex_me
!!***

!!****f* ABINIT/calc_coh
!! NAME
!! calc_coh
!!
!! FUNCTION
!!  Calculates the partial contribution to the COH part of the COHSEX self-energy for a given q-point.
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
!! nfftot=number of points in real space
!! i_sz=contribution arising from the integrable coulombian singularity at q==0
!! (see csigme for the method used), note that in case of 3-D systems the factor
!! 4pi in the coulombian potential is included in the definition of i_sz
!! gvec(3,npwc)=G vectors in reduced coordinates
!! vc_sqrt(npwc)= square root of the coulombian matrix elements for this q-point
!! epsm1q_o(npwc,npwc)= contains $\tilde epsilon^{-1}(q,w=0) - \delta_{G Gp}$ for
!!  the particular q-point considered in the sum
!! wfg2_jk(nsig_ab*nfftot)= Fourier Transform of $\u_{jb k}^*(r) u_{kb k}$
!!  jb,kb=left and righ band indeces definining the left and right states where the
!!  partial contribution to the matrix element of $\Sigma_{COH}$ is evaluated
!!
!! OUTPUT
!! sigcohme=partial contribution to the matrix element of $<jb k \sigma|\Sigma_{COH} | kb k \sigma>$
!!  coming from this single q-point
!!
!! PARENTS
!!      cohsex_me
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine calc_coh(nspinor,nsig_ab,nfftot,ngfft,npwc,gvec,wfg2_jk,epsm1q_o,vc_sqrt,i_sz,iqibz,same_band,sigcohme)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,nfftot,npwc,nsig_ab,nspinor
 real(dp),intent(in) :: i_sz
 logical,intent(in) :: same_band
!arrays
 integer,intent(in) :: gvec(3,npwc),ngfft(18)
 complex(gwpc),intent(in) :: epsm1q_o(npwc,npwc),vc_sqrt(npwc)
 complex(gwpc),intent(in) :: wfg2_jk(nfftot*nsig_ab)
 complex(gwpc),intent(out) :: sigcohme(nsig_ab)

!Local variables-------------------------------
!scalars
 integer,save :: enough=0
 integer :: ig,ig4,ig4x,ig4y,ig4z,igp,igmin,ispinor,spad,outofbox
!arrays
 integer :: g2mg1(3)

! *************************************************************************

 DBG_ENTER("COLL")

 ! === Partial contribution to the matrix element of Sigma_c ===
 ! * For nspinor==2, the closure relation reads:
 !  $\sum_s \psi_a^*(1)\psi_b(2) = \delta_{ab} \delta(1-2)$
 !  where a,b are the spinor components. As a consequence, Sigma_{COH} is always
 !  diagonal in spin-space and only diagonal matrix elements have to be calculated.
 ! MG  TODO wfg2_jk should be calculated on an augmented FFT box to avoid spurious wrapping of G1-G2.
 ! MG: One has to make sure G1-G2 is still in the FFT mesh for each G1 and G2 in chi0 (not always true)
 ! MODULO wraps G1-G2 in the FFT box but the Fourier components are not periodic!

 ! * Treat the case q --> 0 adequately.
 ! TODO Better treatment of wings, check cutoff in the coulombian interaction.
 igmin=1; if (iqibz==1) igmin=2

 sigcohme(:)=czero_gw

 do ispinor=1,nspinor
   spad=(ispinor-1)*nfftot
   outofbox=0

   do igp=igmin,npwc
     do ig=igmin,npwc

      g2mg1 = gvec(:,igp)-gvec(:,ig)
      if (ANY(g2mg1(:)>ngfft(1:3)/2) .or. ANY(g2mg1(:)<-(ngfft(1:3)-1)/2)) then
        outofbox = outofbox+1; CYCLE
      end if

      ig4x=MODULO(g2mg1(1),ngfft(1))
      ig4y=MODULO(g2mg1(2),ngfft(2))
      ig4z=MODULO(g2mg1(3),ngfft(3))
      ig4= 1+ig4x+ig4y*ngfft(1)+ig4z*ngfft(1)*ngfft(2)

      sigcohme(ispinor) = sigcohme(ispinor) + &
&                         half*wfg2_jk(spad+ig4)*epsm1q_o(ig,igp)*vc_sqrt(ig)*vc_sqrt(igp)
     end do !ig
   end do !igp

   if (iqibz==1.and.same_band) then
     sigcohme(ispinor) = sigcohme(ispinor) + half*wfg2_jk(spad+ig4)*epsm1q_o(1,1)*i_sz
   end if
 end do !ispinor

 if (outofbox/=0) then
   enough=enough+1
   if (enough<=50) then
     MSG_WARNING(sjoin(' Number of G1-G2 pairs outside the G-sphere for Wfns:', itoa(outofbox)))
     if (enough==50) then
       call wrtout(std_out,' ========== Stop writing Warnings ==========')
     end if
   end if
 end if

 DBG_EXIT("COLL")

end subroutine calc_coh
!!***

end module m_cohsex
!!***
