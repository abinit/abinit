!{\src2tex{textfont=tt}}
!!****f* ABINIT/cchi0
!! NAME
!! cchi0
!!
!! FUNCTION
!! Main calculation of the independent-particle susceptibility chi0 for qpoint!=0
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! use_tr=If .TRUE. valence states are stored in Wfs_val and only resonant transitions are calculated
!!  (time reversal is assumed)
!! Dtset <type(dataset_type)>=all input variables in this dataset
!! Cryst<crystal_t>= data type gathering info on symmetries and unit cell
!!    %natom=number of atoms
!!    %nsym=number of symmetries
!!    %xred(3,natom)=reduced coordinated of atoms
!!    %typat(natom)=type of each atom
!!    %rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!    %timrev= 2 if time reversal can be used, 1 otherwise
!! qpoint(3)=reciprocal space coordinates of the q wavevector
!! Ep<type(em1params_t_type)>= Parameters related to the calculation of the inverse dielectric matrix.
!!    %nbnds=number of bands summed over
!!    %npwe=number of planewaves for the irreducible polarizability X^0_GGp
!!    %npwvec=maximum number of G vectors
!!     used to define the dimension of some arrays e.g igfft
!!    %nsppol=1 for unpolarized, 2 for spin-polarized
!!    %nomega=total number of frequencies in X^0 (both real and imaginary)
!!    %nomegasf=number of real frequencies used to sample the imaginary part of X^0 (spectral method)
!!    %spmeth=1 if we use the spectral method, 0 for standard Adler-Wiser expression
!!    %spsmear=gaussian broadening used to approximate the delta distribution
!!    %zcut=small imaginary shift to avoid poles in X^0
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! Kmesh <kmesh_t>= datatype gathering parameters related to the k-point sampling
!!    %nibz=number of k-points in the IBZ
!!    %nbz=number of k-points in the BZ
!!    %bz(3,nbz)=reduced coordinates for k-points in the full Brillouin zone
!!    %ibz(3,nibz)=reduced coordinates for k-points in the irreducible wedge
!!    %tab(nbz)=mapping between a kpt in the BZ (array bz) and the irred point in the array ibz
!!    %tabi(nbz)= -1 if inversion is needed to obtain this particular kpt in the BZ, 1 means identity
!!    %tabo(nbz)= for each point in the BZ, the index of the symmetry operation S in reciprocal
!!      space which rotates k_IBZ onto \pm k_BZ (depending on tabi)
!!    %tabp(nbz)= For each k_BZ, it gives the phase factors associated to non-symmorphic operations, i.e
!!      e^{-i 2 \pi k_IBZ \cdot R{^-1}t} == e{-i 2\pi k_BZ cdot t} where :
!!      \transpose R{-1}=S and (S k_IBZ) = \pm k_BZ (depending on ktabi)
!!    %tabr(nfftot,nbz) For each point r on the real mesh and for each k-point in the BZ, tabr
!!      gives the index of (R^-1 (r-t)) in the FFT array where R=\transpose S^{-1} and k_BZ=S k_IBZ.
!!      t is the fractional translation associated to R
!! Gsph_epsG0<gsphere_t data type> The G-sphere used to describe chi0/eps. (including umklapp G0 vectors)
!!    %ng=number of G vectors for chi0
!!    %rottbm1(ng,2,nsym)=contains the index (IS^{-1}) G
!!    %phmGt(Ep%npwe,nsym)=phase factors e^{-iG \cdot t} needed to symmetrize oscillator matrix elements and epsilon
!!    %gprimd(3,3)=dimensional reciprocal space primitive translations (b^-1)
!!    %gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!! nbvw=number of bands in the arrays wfrv
!! ngfft_gw(18)= array containing all the information for 3D FFT for the oscillator strengths (see input variable)
!! nfftot_gw=Total number of points in the GW FFT grid
!! Ltg_q<Little group>=Data type gathering information on the little group of the q-points.
!! Pawtab(Psps%ntypat) <type(pawtab_type)>=paw tabulated starting data
!! Pawang<pawang_type> angular mesh discretization and related data:
!! QP_BSt<ebands_t>=Quasiparticle energies and occupations (for the moment real quantities)
!!   %mband=MAX number of bands over k-points and spin (==Ep%nbnds)
!!   %occ(mband,nkpt,nsppol)=QP occupation numbers, for each k point in IBZ, and each band
!!   %eig(mband,nkpt,nsppol)=GW energies, for self-consistency purposes
!!  Paw_pwff<pawpwff_t>=Form factor used to calculate the onsite mat. elements of a plane wave.
!!  Wfd<wfd_t>=Object used to access the wavefunctions
!!
!! OUTPUT
!!  chi0(Ep%npwe,Ep%npwe,Ep%nomega)=independent-particle susceptibility matrix at wavevector qpoint and
!!   each frequeny defined by Ep%omega and Ep%nomega
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      accumulate_chi0sumrule,approxdelta,assemblychi0_sym,assemblychi0sf
!!      calc_wfwfg,chi0_bbp_mask,completechi0_deltapart,cwtime,flush_unit
!!      get_bz_diff,get_bz_item,get_gftt,get_uug,gsph_fft_tabs,gsph_free
!!      gsph_in_fftbox,hilbert_transform,littlegroup_print,make_transitions
!!      paw_cross_rho_tw_g,paw_rho_tw_g,paw_symcprj,pawcprj_alloc,pawcprj_free
!!      pawpwij_free,pawpwij_init,read_plowannier,rho_tw_g,setup_spectral
!!      symmetrize_afm_chi0,timab,wfd_change_ngfft,wfd_distribute_kb_kpbp
!!      wfd_get_cprj,wfd_get_ur,wfd_paw_get_aeur,wfd_sym_ur,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine cchi0(use_tr,Dtset,Cryst,qpoint,Ep,Psps,Kmesh,QP_BSt,Gsph_epsG0,&
& Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,nbvw,ngfft_gw,nfftot_gw,ngfftf,nfftf_tot,&
& chi0,ktabr,ktabrf,Ltg_q,chi0_sumrule,Wfd,Wfdf)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_xmpi
 use m_blas
 use m_errors
 use m_profiling_abi
 use m_time

 use m_gwdefs,        only : GW_TOL_DOCC, GW_TOL_W0, czero_gw, em1params_t, g0g0w
 use m_numeric_tools, only : imin_loc
 use m_geometry,      only : normv
 use m_io_tools,      only : flush_unit
 use m_crystal,       only : crystal_t
 use m_bz_mesh,       only : kmesh_t, get_BZ_item, get_BZ_diff, littlegroup_t, littlegroup_print
 use m_gsphere,       only : gsphere_t, gsph_fft_tabs,  gsph_free, gsph_in_fftbox
 use m_fft_mesh,      only : get_gftt
 use m_wfd,           only : wfd_get_ur, wfd_t, wfd_distribute_kb_kpbp, wfd_get_cprj, wfd_barrier, wfd_change_ngfft,&
&                            wfd_paw_get_aeur, wfd_sym_ur
 use m_oscillators,   only : rho_tw_g, calc_wfwfg, get_uug
 use m_chi0,          only : hilbert_transform, setup_spectral, assemblychi0_sym, assemblychi0sf, symmetrize_afm_chi0,&
&                            approxdelta
 use m_pawang,        only : pawang_type
 use m_pawtab,        only : pawtab_type
 use m_pawfgrtab,     only : pawfgrtab_type
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_pawpwij,      only : pawpwff_t, pawpwij_t, pawpwij_init, pawpwij_free, paw_rho_tw_g, paw_cross_rho_tw_g
 use m_paw_pwaves_lmn,only : paw_pwaves_lmn_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cchi0'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_65_paw
 use interfaces_70_gw, except_this_one => cchi0
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbvw,nfftot_gw,nfftf_tot
 logical,intent(in) :: use_tr
 type(ebands_t),target,intent(in) :: QP_BSt
 type(kmesh_t),intent(in) :: Kmesh
 type(crystal_t),intent(in) :: Cryst
 type(Dataset_type),intent(in) :: Dtset
 type(em1params_t),intent(in) :: Ep
 type(gsphere_t),intent(in) :: Gsph_epsG0
 type(littlegroup_t),intent(in) :: Ltg_q
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
 type(wfd_t),target,intent(inout) :: Wfd,Wfdf
!arrays
 integer,intent(in) :: ktabr(nfftot_gw,Kmesh%nbz),ktabrf(nfftf_tot*Dtset%pawcross,Kmesh%nbz)
 integer,intent(in) :: ngfft_gw(18),ngfftf(18)
 real(dp),intent(in) :: qpoint(3)
 real(dp),intent(out) :: chi0_sumrule(Ep%npwe)
 complex(gwpc),intent(out) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
 type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp1=1,two_poles=2,one_pole=1,ndat1=1
 integer :: bandinf,bandsup,dim_rtwg,band1,band2,ierr
 integer :: ig1,ig2,iat,ik_bz,ik_ibz,ikmq_bz,ikmq_ibz
 integer :: io,iomegal,iomegar,ispinor1,ispinor2,isym_k,itypatcor,nfft
 integer :: isym_kmq,itim_k,itim_kmq,m1,m2,my_wl,my_wr,size_chi0
 integer :: nfound,nkpt_summed,nspinor,nsppol,mband
 integer :: comm,gw_mgfft,use_padfft,gw_fftalga,lcor,mgfftf,use_padfftf
 integer :: my_nbbp,my_nbbpks,spin
 integer :: nbmax,dummy
 real(dp) :: cpu_time,wall_time,gflops
 real(dp) :: deltaeGW_b1kmq_b2k,deltaeGW_enhigh_b2k,deltaf_b1kmq_b2k
 real(dp) :: e_b1_kmq,en_high,fac,fac2,fac3,f_b1_kmq,factor,max_rest,min_rest,my_max_rest
 real(dp) :: my_min_rest,numerator,spin_fact,weight,wl,wr
 real(dp) :: gw_gsq,memreq
 complex(dpc) :: ph_mkmqt,ph_mkt
 complex(gwpc) :: local_czero_gw
 logical :: qzero,isirred_k,isirred_kmq,luwindow
 character(len=500) :: msg,allup
 type(gsphere_t) :: Gsph_FFT
!arrays
 integer :: G0(3),umklp_k(3),umklp_kmq(3)
 integer :: ucrpa_bands(2)
 integer :: wtk_ltg(Kmesh%nbz),got(Wfd%nproc)
 integer,allocatable :: tabr_k(:),tabr_kmq(:),tabrf_k(:),tabrf_kmq(:)
 integer,allocatable :: igfftepsG0(:),gspfft_igfft(:),igfftepsG0f(:)
 integer,allocatable :: gw_gfft(:,:),gw_gbound(:,:),dummy_gbound(:,:),gboundf(:,:)
 integer,allocatable :: bbp_ks_distrb(:,:,:,:)
 real(dp) :: kbz(3),kmq_bz(3),spinrot_k(4),spinrot_kmq(4),q0(3),tsec(2)
 real(dp),ABI_CONTIGUOUS pointer :: qp_energy(:,:,:),qp_occ(:,:,:)
 real(dp),allocatable :: omegasf(:)
 complex(dpc),allocatable :: green_enhigh_w(:),green_w(:),kkweight(:,:)
 complex(gwpc),allocatable :: sf_chi0(:,:,:),rhotwg(:)
 complex(gwpc),allocatable :: ur1_kmq_ibz(:),ur2_k_ibz(:),wfwfg(:)
 complex(gwpc),allocatable :: usr1_kmq(:),ur2_k(:)
 complex(gwpc),allocatable :: ur_ae1(:),ur_ae_onsite1(:),ur_ps_onsite1(:)
 complex(gwpc),allocatable :: ur_ae2(:),ur_ae_onsite2(:),ur_ps_onsite2(:)
 complex(dpc), allocatable :: coeffW_BZ(:,:,:,:,:,:)
 logical,allocatable :: bbp_mask(:,:)
 type(pawcprj_type),allocatable :: Cprj1_kmq(:,:),Cprj2_k(:,:)
 type(pawpwij_t),allocatable :: Pwij(:),Pwij_fft(:)
!************************************************************************

#define DEV_USE_OLDRHOTWG 1

 DBG_ENTER("COLL")

 call timab(331,1,tsec) ! cchi0
 call cwtime(cpu_time,wall_time,gflops,"start")

 nsppol  = Wfd%nsppol
 nspinor = Wfd%nspinor
 ucrpa_bands(1)=dtset%ucrpa_bands(1)
 ucrpa_bands(2)=dtset%ucrpa_bands(2)
 luwindow=.false.
 if(abs(dtset%ucrpa_window(1)+1_dp)>tol8.and.(abs(dtset%ucrpa_window(2)+1_dp)>tol8)) then
   luwindow=.true.
 endif
! write(6,*)"ucrpa_bands",ucrpa_bands
! write(6,*)"ucrpa_window",dtset%ucrpa_window
! write(6,*)"luwindow",luwindow

!  For cRPA calculation of U: read forlb.ovlp
 if(dtset%ucrpa>=1) then
   call read_plowannier(Cryst,bandinf,bandsup,coeffW_BZ,itypatcor,Kmesh,lcor,luwindow,&
& nspinor,nsppol,pawang,dtset%prtvol,ucrpa_bands)
 endif
! End of reading forlb.ovlp

 if ( ANY(ngfft_gw(1:3) /= Wfd%ngfft(1:3)) ) then
   call wfd_change_ngfft(Wfd,Cryst,Psps,ngfft_gw)
 end if
 gw_mgfft = MAXVAL(ngfft_gw(1:3))
 gw_fftalga = ngfft_gw(7)/100 !; gw_fftalgc=MOD(ngfft_gw(7),10)

 if (Dtset%pawcross==1) mgfftf = MAXVAL(ngfftf(1:3))
 !
 ! === Initialize MPI variables ===
 comm = Wfd%comm
 !
 ! == Copy some values ===
 mband   = Wfd%mband
 nfft    = Wfd%nfft
 ABI_CHECK(Wfd%nfftot==nfftot_gw,"Wrong nfftot_gw")

 dim_rtwg=1; if (nspinor==2) dim_rtwg=4    !can reduce size depending on Ep%nI and Ep%nj
 size_chi0 = Ep%npwe*Ep%nI*Ep%npwe*Ep%nJ*Ep%nomega

 qp_energy => QP_BSt%eig(:,:,:)
 qp_occ    => QP_BSt%occ(:,:,:)
 !
 ! === Initialize the completeness correction  ===
 if (Ep%gwcomp==1) then
   en_high=MAXVAL(qp_energy(Ep%nbnds,:,:)) + Ep%gwencomp
   write(msg,'(a,f8.2,a)')' Using completeness correction with the energy ',en_high*Ha_eV,' [eV]'
   call wrtout(std_out,msg,'COLL')
   !
   ! Allocation of wfwfg and green_enhigh_w moved inside openmp loop
   !
   ! Init the largest G-sphere contained in the FFT box for the wavefunctions.
   call gsph_in_fftbox(Gsph_FFT,Cryst,Wfd%ngfft)

   !call print_gsphere(Gsph_FFT,unit=std_out,prtvol=10)

   ABI_MALLOC(gspfft_igfft,(Gsph_FFT%ng))
   ABI_MALLOC(dummy_gbound,(2*gw_mgfft+8,2))
   !
   ! Mapping between G-sphere and FFT box.
   call gsph_fft_tabs(Gsph_FFT,(/0,0,0/),Wfd%mgfft,Wfd%ngfft,dummy,dummy_gbound,gspfft_igfft)
   ABI_FREE(dummy_gbound)
   !
   if (Psps%usepaw==1) then  ! * Prepare the onsite contributions on the GW FFT mesh.
     ABI_MALLOC(gw_gfft,(3,nfft))
     q0=zero
     call get_gftt(ngfft_gw,q0,Cryst%gmet,gw_gsq,gw_gfft) ! Get the set of plane waves in the FFT Box.
     ABI_DT_MALLOC(Pwij_fft,(Psps%ntypat))
     call pawpwij_init(Pwij_fft,nfft,(/zero,zero,zero/),gw_gfft,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
   end if
   !
 end if
 !
 ! === Setup weights (2 for spin unpolarized sistem, 1 for polarized) ===
 ! * spin_fact is used to normalize the occupation factors to one. Consider also the AFM case.
 SELECT CASE (nsppol)
 CASE (1)
   weight=two/Kmesh%nbz; spin_fact=half
   if (Wfd%nspden==2) then
    weight=one/Kmesh%nbz; spin_fact=half
   end if
   if (nspinor==2) then
    weight=one/Kmesh%nbz; spin_fact=one
   end if
 CASE (2)
   weight=one/Kmesh%nbz; spin_fact=one
 CASE DEFAULT
   MSG_BUG("Wrong nsppol")
 END SELECT
 !
 ! === Weight for points in the IBZ_q ===
 wtk_ltg(:)=1
 if (Ep%symchi==1) then
   do ik_bz=1,Ltg_q%nbz
     wtk_ltg(ik_bz)=0
     if (Ltg_q%ibzq(ik_bz)/=1) CYCLE ! Only k-points in the IBZ_q.
     wtk_ltg(ik_bz)=SUM(Ltg_q%wtksym(:,:,ik_bz))
   end do
 end if

 write(msg,'(a,i2,2a,i2)')&
&  ' Using spectral method for the imaginary part = ',Ep%spmeth,ch10,&
&  ' Using symmetries to sum only over the IBZ_q  = ',Ep%symchi
 call wrtout(std_out,msg,'COLL')

 if (use_tr) then
   ! Special care has to be taken in metals and/or spin dependent systems
   ! as Wfs_val might contain unoccupied states.
   call wrtout(std_out,' Using faster algorithm based on time reversal symmetry.','COLL')
 end if

 ! TODO this table can be calculated for each k-point
 my_nbbpks=0; allup="All"; got=0
 ABI_MALLOC(bbp_ks_distrb,(mband,mband,Kmesh%nbz,nsppol))
 ABI_MALLOC(bbp_mask,(mband,mband))

 do spin=1,nsppol
   do ik_bz=1,Kmesh%nbz
     if (Ep%symchi==1) then
       if (Ltg_q%ibzq(ik_bz)/=1) CYCLE  ! Only IBZ_q
     end if
     !
     ! * Get ik_ibz, non-symmorphic phase, ph_mkt, and symmetries from ik_bz.
     call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k)
     !
     ! * Get index of k-q in the BZ, stop if not found as the weight=one/nkbz is not correct.
     call get_BZ_diff(Kmesh,kbz,qpoint,ikmq_bz,g0,nfound)
     ABI_CHECK(nfound==1,"Check kmesh")
     !
     ! * Get ikmq_ibz, non-symmorphic phase, ph_mkmqt, and symmetries from ikmq_bz.
     call get_BZ_item(Kmesh,ikmq_bz,kmq_bz,ikmq_ibz,isym_kmq,itim_kmq)

     call chi0_bbp_mask(Ep,use_tr,QP_BSt,mband,ikmq_ibz,ik_ibz,spin,spin_fact,bbp_mask)

     call wfd_distribute_kb_kpbp(Wfd,ikmq_ibz,ik_ibz,spin,allup,my_nbbp,bbp_ks_distrb(:,:,ik_bz,spin),got,bbp_mask)
     my_nbbpks = my_nbbpks + my_nbbp
   end do
 end do

 ABI_FREE(bbp_mask)

 write(msg,'(a,i0,a)')" Will sum ",my_nbbpks," (b,b',k,s) states in chi0."
 call wrtout(std_out,msg,'PERS')

 if (Psps%usepaw==1) then
   ABI_DT_MALLOC(Pwij,(Psps%ntypat))
   call pawpwij_init(Pwij,Ep%npwepG0,qpoint,Gsph_epsG0%gvec,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
   ! Allocate statements moved to inside openmp loop
 end if

 SELECT CASE (Ep%spmeth)

 CASE (0)
   call wrtout(std_out,' Calculating chi0(q,omega,G,G")','COLL')
   ! Allocation of green_w moved inside openmp loop

 CASE (1,2)
   call wrtout(std_out,' Calculating Im chi0(q,omega,G,G")','COLL')
   !
   ! Find Max and min resonant transitions for this q, report also treated by this proc.
   call make_transitions(Wfd,1,Ep%nbnds,nbvw,nsppol,Ep%symchi,Cryst%timrev,GW_TOL_DOCC,&
&    max_rest,min_rest,my_max_rest,my_min_rest,Kmesh,Ltg_q,qp_energy,qp_occ,qpoint,bbp_ks_distrb)
   !
   ! Calculate frequency dependent weights for Hilbert transform.
   ABI_MALLOC(omegasf,(Ep%nomegasf))
   ABI_MALLOC(kkweight,(Ep%nomegasf,Ep%nomega))
   !my_wl=1; my_wr=Ep%nomegasf
   call setup_spectral(Ep%nomega,Ep%omega,Ep%nomegasf,omegasf,max_rest,min_rest,my_max_rest,my_min_rest,&
&    0,Ep%zcut,zero,my_wl,my_wr,kkweight)

   if (.not.use_tr) then
     MSG_BUG('spectral method requires time-reversal')
   end if

   memreq = two*gwpc*Ep%npwe**2*(my_wr-my_wl+1)*b2Gb
   write(msg,'(a,f10.4,a)')' memory required per spectral point: ',2.0*gwpc*Ep%npwe**2*b2Mb,' [Mb]'
   call wrtout(std_out,msg,'PERS')
   write(msg,'(a,f10.4,a)')' memory required by sf_chi0: ',memreq,' [Gb]'
   call wrtout(std_out,msg,'PERS')
   if (memreq > two) then
     MSG_WARNING(' Memory required for sf_chi0 is larger than 2.0 Gb!')
   end if
   if ((Dtset%userre > zero).and.(memreq > Dtset%userre)) then
     write(msg,'(a,f8.3,a)')' Memory required by sf_chi0 is larger than userre:',Dtset%userre,'[Gb]'
     MSG_ERROR(msg)
   end if

   ABI_STAT_MALLOC(sf_chi0,(Ep%npwe,Ep%npwe,my_wl:my_wr), ierr)
   ABI_CHECK(ierr==0, 'out-of-memory in sf_chi0')
   sf_chi0=czero_gw

 CASE DEFAULT
   MSG_BUG("Wrong spmeth")
 END SELECT

 nkpt_summed=Kmesh%nbz
 if (Ep%symchi==1) then
   nkpt_summed=Ltg_q%nibz_ltg
   call littlegroup_print(Ltg_q,std_out,Dtset%prtvol,'COLL')
 end if

 write(msg,'(a,i6,a)')' Calculation status : ',nkpt_summed,' to be completed '
 call wrtout(std_out,msg,'COLL')
 !
 ! ============================================
 ! === Begin big fat loop over transitions ===
 ! ============================================
 chi0=czero_gw; chi0_sumrule=zero

 ! === Loop on spin to calculate trace $\chi_{up,up}+\chi_{down,down}$ ===
 ! * Only $\chi_{up,up} for AFM.
 do spin=1,nsppol

   if (ALL(bbp_ks_distrb(:,:,:,spin) /= Wfd%my_rank)) CYCLE

   ! Allocation of arrays that are private to loop
   if (Ep%gwcomp==1)  then
     ABI_MALLOC(wfwfg,(nfft*nspinor**2))
   end if
   if (Ep%gwcomp==1)  then
     ABI_MALLOC(green_enhigh_w,(Ep%nomega))
   end if
   if (Ep%spmeth==0)  then
     ABI_MALLOC(green_w,(Ep%nomega))
   end if
   if (Psps%usepaw==1) then
     ABI_DT_MALLOC(Cprj2_k  ,(Cryst%natom,nspinor))
     call pawcprj_alloc(Cprj2_k,  0,Wfd%nlmn_atm)
     ABI_DT_MALLOC(Cprj1_kmq,(Cryst%natom,nspinor))
     call pawcprj_alloc(Cprj1_kmq,0,Wfd%nlmn_atm)
     if (Dtset%pawcross==1) then
       ABI_MALLOC(ur_ae1,(nfftf_tot*nspinor))
       ABI_MALLOC(ur_ae_onsite1,(nfftf_tot*nspinor))
       ABI_MALLOC(ur_ps_onsite1,(nfftf_tot*nspinor))
       ABI_MALLOC(ur_ae2,(nfftf_tot*nspinor))
       ABI_MALLOC(ur_ae_onsite2,(nfftf_tot*nspinor))
       ABI_MALLOC(ur_ps_onsite2,(nfftf_tot*nspinor))
       ABI_MALLOC(igfftepsG0f,(Ep%npwepG0))
       ABI_MALLOC(tabrf_k,(nfftf_tot))
       ABI_MALLOC(tabrf_kmq,(nfftf_tot))
     end if
   end if

   ABI_MALLOC(rhotwg,(Ep%npwepG0*nspinor**2))
   ABI_MALLOC(tabr_k,(nfft))
   ABI_MALLOC(tabr_kmq,(nfft))

   ABI_MALLOC(ur1_kmq_ibz,(nfft*nspinor))
   ABI_MALLOC(ur2_k_ibz,(nfft*nspinor))

   ABI_MALLOC(usr1_kmq,(nfft*nspinor))
   ABI_MALLOC(ur2_k,   (nfft*nspinor))

   ABI_MALLOC(igfftepsG0,(Ep%npwepG0))

   do ik_bz=1,Kmesh%nbz ! Loop over k-points in the BZ.

     if (Ep%symchi==1) then
       if (Ltg_q%ibzq(ik_bz)/=1) CYCLE  ! Only IBZ_q
     end if

     if (ALL(bbp_ks_distrb(:,:,ik_bz,spin) /= Wfd%my_rank)) CYCLE

     write(msg,'(2(a,i4),a,i2,a,i3)')' ik = ',ik_bz,' / ',Kmesh%nbz,' spin = ',spin,' done by processor ',Wfd%my_rank
     call wrtout(std_out,msg,'PERS')
     !
     ! * Get ik_ibz, non-symmorphic phase, ph_mkt, and symmetries from ik_bz.
     call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k,ph_mkt,umklp_k,isirred_k)

     call get_BZ_diff(Kmesh,kbz,qpoint,ikmq_bz,G0,nfound)
     if (nfound==0) then
       MSG_ERROR("Cannot find kbz - qpoint in Kmesh")
     end if

     ! * Get ikmq_ibz, non-symmorphic phase, ph_mkmqt, and symmetries from ikmq_bz.
     call get_BZ_item(Kmesh,ikmq_bz,kmq_bz,ikmq_ibz,isym_kmq,itim_kmq,ph_mkmqt,umklp_kmq,isirred_kmq)

!BEGIN DEBUG
     if (ANY(umklp_k /=0)) then
       write(msg,'(a,3i2)')" umklp_k /= 0 ",umklp_k
       MSG_ERROR(msg)
     end if
     !if (ANY( g0 /= -umklp_kmq + umklp_k) ) then
     !if (ANY( g0 /= -umklp_kmq ) ) then
     !  write(msg,'(a,3(1x,3i2))')" g0 /= -umklp_kmq + umklp_k ",g0, umklp_kmq, umklp_k
     !  MSG_ERROR(msg)
     !end if
     !g0 = -umklp_k + umklp_kmq
     !g0 = +umklp_k - umklp_kmq
     !if (ANY (ABS(g0) > Ep%mg0) ) then
     !  write(msg,'(a,6(1x,i0))')"  ABS(g0) > Ep%mg0 ",g0,Ep%mg0
     !  MSG_ERROR(msg)
     !end if
!END DEBUG

     ! * Copy tables for rotated FFT points
     tabr_k(:)  =ktabr(:,ik_bz)
     spinrot_k(:)=Cryst%spinrot(:,isym_k)

     tabr_kmq(:)=ktabr(:,ikmq_bz)
     spinrot_kmq(:)=Cryst%spinrot(:,isym_kmq)

     if (Dtset%pawcross==1) then
       tabrf_k(:)  =ktabrf(:,ik_bz)
       tabrf_kmq(:)=ktabrf(:,ikmq_bz)
     end if
     !
     ! Tables for the FFT of the oscillators.
     !  a) FFT index of G-G0.
     !  b) gw_gbound table for the zero-padded FFT performed in rhotwg.
     ABI_MALLOC(gw_gbound,(2*gw_mgfft+8,2))
     call gsph_fft_tabs(Gsph_epsG0,g0,gw_mgfft,ngfft_gw,use_padfft,gw_gbound,igfftepsG0)
     if ( ANY(gw_fftalga == (/2,4/)) ) use_padfft=0 ! Pad-FFT is not coded in rho_tw_g
     if (use_padfft==0) then
       ABI_FREE(gw_gbound)
       ABI_MALLOC(gw_gbound,(2*gw_mgfft+8,2*use_padfft))
     end if

     if (Dtset%pawcross==1) then
        ABI_MALLOC(gboundf,(2*mgfftf+8,2))
       call gsph_fft_tabs(Gsph_epsG0,g0,mgfftf,ngfftf,use_padfftf,gboundf,igfftepsG0f)
       if (ANY(gw_fftalga == (/2,4/))) use_padfftf=0
       if (use_padfftf==0) then
         ABI_FREE(gboundf)
         ABI_MALLOC(gboundf,(2*mgfftf+8,2*use_padfftf))
       end if
     end if

     nbmax=Ep%nbnds
     do band1=1,nbmax ! Loop over "conduction" states.
       if (ALL(bbp_ks_distrb(band1,:,ik_bz,spin) /= Wfd%my_rank)) CYCLE

#if DEV_USE_OLDRHOTWG
       call wfd_get_ur(Wfd,band1,ikmq_ibz,spin,ur1_kmq_ibz)
#else
       call wfd_sym_ur(Wfd,Cryst,Kmesh,band1,ikmq_bz,spin,usr1_kmq,trans="C",with_umklp=.FALSE.)
#endif

       if (Psps%usepaw==1) then
         call wfd_get_cprj(Wfd,band1,ikmq_ibz,spin,Cryst,Cprj1_kmq,sorted=.FALSE.)
         call paw_symcprj(ikmq_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,Cprj1_kmq)
         if (Dtset%pawcross==1) then
           call wfd_paw_get_aeur(Wfdf,band1,ikmq_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae1,ur_ae_onsite1,ur_ps_onsite1)
         end if
       end if

       e_b1_kmq=qp_energy(band1,ikmq_ibz,spin)
       f_b1_kmq=   qp_occ(band1,ikmq_ibz,spin)

       do band2=1,nbmax ! Loop over "valence" states.
!debug         if (.not.luwindow.AND.dtset%ucrpa==1.AND.band1<=ucrpa_bands(2).AND.band1>=ucrpa_bands(1)&
!debug&                                            .AND.band2<=ucrpa_bands(2).AND.band2>=ucrpa_bands(1)) CYClE
       !write(6,*) "ik,band1,band2",ik_bz,band1,band2
         if (luwindow.AND.dtset%ucrpa==1.AND.((QP_Bst%eig(band1,ik_ibz   ,spin)-QP_Bst%fermie)<=dtset%ucrpa_window(2)) &
&                                       .AND.((QP_Bst%eig(band1,ik_ibz   ,spin)-QP_Bst%fermie)>=dtset%ucrpa_window(1)) &
&                                       .AND.((QP_Bst%eig(band2,ikmq_ibz,spin)-QP_Bst%fermie)<=dtset%ucrpa_window(2)) &
&                                       .AND.((QP_Bst%eig(band2,ikmq_ibz,spin)-QP_Bst%fermie)>=dtset%ucrpa_window(1))) CYCLE

         if (bbp_ks_distrb(band1,band2,ik_bz,spin) /= Wfd%my_rank) CYCLE

         deltaf_b1kmq_b2k=spin_fact*(f_b1_kmq-qp_occ(band2,ik_ibz,spin))

         if (Ep%gwcomp==0) then ! Skip negligible transitions.
           if (ABS(deltaf_b1kmq_b2k) < GW_TOL_DOCC) CYCLE
         else
           ! * When the completeness correction is used,
           !   we need to also consider transitions with vanishing deltaf
           !if (qp_occ(band2,ik_ibz,spin) < GW_TOL_DOCC) CYCLE
           !
           ! Rangel This is to compute chi correctly when using the extrapolar method
           if (qp_occ(band2,ik_ibz,spin) < GW_TOL_DOCC .and. (ABS(deltaf_b1kmq_b2k) < GW_TOL_DOCC .or. band1<band2)) CYCLE
         end if

         deltaeGW_b1kmq_b2k=e_b1_kmq-qp_energy(band2,ik_ibz,spin)

#if DEV_USE_OLDRHOTWG
         call wfd_get_ur(Wfd,band2,ik_ibz,spin,ur2_k_ibz)
#else
         call wfd_sym_ur(Wfd,Cryst,Kmesh,band2,ik_bz,spin,ur2_k,with_umklp=.FALSE.,ur_kibz=ur2_k_ibz)
#endif

         if (Psps%usepaw==1) then
           call wfd_get_cprj(Wfd,band2,ik_ibz,spin,Cryst,Cprj2_k,sorted=.FALSE.)
           call paw_symcprj(ik_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,Cprj2_k)
           if (Dtset%pawcross==1) then
             call wfd_paw_get_aeur(Wfdf,band2,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae2,ur_ae_onsite2,ur_ps_onsite2)
           end if
         end if

         SELECT CASE (Ep%spmeth)

         CASE (0) ! Standard Adler-Wiser expression.
          ! * Add the small imaginary of the Time-Ordered RF only for non-zero real omega ! FIXME What about metals?

          if (.not.use_tr) then ! Have to sum over all possible resonant and anti-resonant transitions.
            do io=1,Ep%nomega
              green_w(io) = g0g0w(Ep%omega(io),deltaf_b1kmq_b2k,deltaeGW_b1kmq_b2k,Ep%zcut,GW_TOL_W0,one_pole)
            end do
          else

            if (Ep%gwcomp==0) then ! cannot be completely skipped in case of completeness correction
              if (band1<band2) CYCLE ! Here we GAIN a factor ~2
            end if

            do io=1,Ep%nomega
              !Rangel: In metals, the intra-band transitions term does not contain the antiresonant part
              !green_w(io) = g0g0w(Ep%omega(io),deltaf_b1kmq_b2k,deltaeGW_b1kmq_b2k,Ep%zcut,GW_TOL_W0)
              if (band1==band2) green_w(io) = g0g0w(Ep%omega(io),deltaf_b1kmq_b2k,deltaeGW_b1kmq_b2k,Ep%zcut,GW_TOL_W0,one_pole)
              if (band1/=band2) green_w(io) = g0g0w(Ep%omega(io),deltaf_b1kmq_b2k,deltaeGW_b1kmq_b2k,Ep%zcut,GW_TOL_W0,two_poles)

              if (Ep%gwcomp==1) then ! Calculate the completeness correction
                numerator= -spin_fact*qp_occ(band2,ik_ibz,spin)
                deltaeGW_enhigh_b2k = en_high-qp_energy(band2,ik_ibz,spin)

                if (REAL(Ep%omega(io))<GW_TOL_W0) then ! Completeness correction is NOT valid for real frequencies
                  green_enhigh_w(io) = g0g0w(Ep%omega(io),numerator,deltaeGW_enhigh_b2k,Ep%zcut,GW_TOL_W0,two_poles)
                else
                  green_enhigh_w(io) = local_czero_gw
                end if
                !
                !Rangel Correction for metals
                !if (deltaf_b1kmq_b2k<0.d0) then
                if (band1>=band2 .and. ABS(deltaf_b1kmq_b2k) > GW_TOL_DOCC ) then
                  green_w(io)= green_w(io) - green_enhigh_w(io)
                else ! Disregard green_w, since it is already accounted for through the time-reversal
                  green_w(io)=             - green_enhigh_w(io)
                end if
              end if !gwcomp==1
            end do !io
            !
            if (Ep%gwcomp==1.and.band1==band2) then
              ! Add the "delta part" of the extrapolar method. TODO doesnt work for spinor

              call calc_wfwfg(tabr_k,itim_k,nfft,ngfft_gw,ur2_k_ibz,ur2_k_ibz,wfwfg)

              if (Psps%usepaw==1) then
                call paw_rho_tw_g(nfft,dim_rtwg,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,gw_gfft,&
&                 Cprj2_k,Cprj2_k,Pwij_fft,wfwfg)

                ! Add PAW cross term
                if (Dtset%pawcross==1) then
                  call paw_cross_rho_tw_g(nspinor,Ep%npwepG0,nfftf_tot,ngfftf,1,use_padfftf,igfftepsG0f,gboundf,&
&                  ur_ae2,ur_ae_onsite2,ur_ps_onsite2,itim_kmq,tabrf_kmq,ph_mkmqt,spinrot_kmq,&
&                  ur_ae2,ur_ae_onsite2,ur_ps_onsite2,itim_k  ,tabrf_k  ,ph_mkt  ,spinrot_k,dim_rtwg,wfwfg)
                end if
              end if

              qzero=.FALSE.
              call completechi0_deltapart(ik_bz,qzero,Ep%symchi,Ep%npwe,Gsph_FFT%ng,Ep%nomega,nspinor,&
&               nfft,ngfft_gw,gspfft_igfft,gsph_FFT,Ltg_q,green_enhigh_w,wfwfg,chi0)

            end if
          end if ! use_tr

         CASE (1,2) ! Spectral method, WARNING time-reversal here is always assumed!
           if (deltaeGW_b1kmq_b2k<0) CYCLE
           call approxdelta(Ep%nomegasf,omegasf,deltaeGW_b1kmq_b2k,Ep%spsmear,iomegal,iomegar,wl,wr,Ep%spmeth)
         END SELECT
         !
         ! ==== Form rho-twiddle(r)=u^*_{b1,kmq_bz}(r) u_{b2,kbz}(r) and its FFT transform ====

#if DEV_USE_OLDRHOTWG
         call rho_tw_g(nspinor,Ep%npwepG0,nfft,ndat1,ngfft_gw,1,use_padfft,igfftepsG0,gw_gbound,&
&          ur1_kmq_ibz,itim_kmq,tabr_kmq,ph_mkmqt,spinrot_kmq,&
&          ur2_k_ibz,  itim_k  ,tabr_k  ,ph_mkt  ,spinrot_k,dim_rtwg,rhotwg)
#else
         call get_uug(Ep%npwepG0,nfft,ndat1,ngfft_gw,use_padfft,igfftepsG0,gw_gbound,usr1_kmq,ur2_k,rhotwg)
#endif

         if (Psps%usepaw==1) then! Add PAW on-site contribution, projectors are already in the BZ.
           call paw_rho_tw_g(Ep%npwepG0,dim_rtwg,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_epsG0%gvec,&
&           Cprj1_kmq,Cprj2_k,Pwij,rhotwg)

           ! Add PAW cross term
           if (Dtset%pawcross==1) then
             call paw_cross_rho_tw_g(nspinor,Ep%npwepG0,nfftf_tot,ngfftf,1,use_padfftf,igfftepsG0f,gboundf,&
&             ur_ae1,ur_ae_onsite1,ur_ps_onsite1,itim_kmq,tabrf_kmq,ph_mkmqt,spinrot_kmq,&
&             ur_ae2,ur_ae_onsite2,ur_ps_onsite2,itim_k  ,tabrf_k  ,ph_mkt  ,spinrot_k,dim_rtwg,rhotwg)
           end if
         end if

         SELECT CASE (Ep%spmeth)

         CASE (0) ! Adler-Wiser.
!debug           if(dtset%ucrpa==2)  then
           if(dtset%ucrpa>=1.and..not.luwindow)  then
             fac=one
             fac2=one
             fac3=one
             m1=-1
             m2=-1
             call flush_unit(std_out)
             if(dtset%ucrpa<=2) then
               if (       band1<=ucrpa_bands(2).AND.band1>=ucrpa_bands(1)&
&                    .AND.band2<=ucrpa_bands(2).AND.band2>=ucrpa_bands(1)) then
                 do iat=1, cryst%nattyp(itypatcor)
                   do ispinor1=1,nspinor
                     do ispinor2=1,nspinor
                       do m1=1,2*lcor+1
                         do m2=1,2*lcor+1
                           fac=fac - real(coeffW_BZ(iat,spin,band1,ik_bz,ispinor1,m1)*&
&                                     conjg(coeffW_BZ(iat,spin,band1,ik_bz,ispinor1,m1))* &
&                                     coeffW_BZ(iat,spin,band2,ikmq_bz,ispinor2,m2)*&
&                                     conjg(coeffW_BZ(iat,spin,band2,ikmq_bz,ispinor2,m2)))
                         enddo
                       enddo
                     enddo
                   enddo
                 enddo
                 if(dtset%ucrpa==1) fac=zero
               endif
             else if (dtset%ucrpa>=3) then
               if (band1<=ucrpa_bands(2).AND.band1>=ucrpa_bands(1)) then
                 do iat=1, cryst%nattyp(itypatcor)
                   do ispinor1=1,nspinor
                     do m1=1,2*lcor+1
                       fac2=fac2-real(coeffW_BZ(iat,spin,band1,ik_bz,ispinor1,m1)*&
&                                 conjg(coeffW_BZ(iat,spin,band1,ik_bz,ispinor1,m1)))
                     enddo
                   enddo
                 enddo
                 if(dtset%ucrpa==4) fac2=zero
               endif
               if (band2<=ucrpa_bands(2).AND.band2>=ucrpa_bands(1)) then
                 do iat=1, cryst%nattyp(itypatcor)
                   do ispinor1=1,nspinor
                     do m1=1,2*lcor+1
                       fac3=fac3-real(coeffW_BZ(iat,spin,band2,ikmq_bz,ispinor1,m1)*&
&                                 conjg(coeffW_BZ(iat,spin,band2,ikmq_bz,ispinor1,m1)))
                     enddo
                   enddo
                 enddo
                 if(dtset%ucrpa==4) fac3=zero
               endif
               fac=real(fac2*fac3)
             endif

!             if(dtset%prtvol>=10) write(6,'(6i3,e15.5,a)') ik_bz,ikmq_bz,band1,band2,m1,m2,fac," q=/0"
!             if(dtset%prtvol>=10.and.abs(fac-one)>0.00001) &
!&             write(6,'(a,6i4,e15.5,a)') "*****FAC*********",ik_bz,ikmq_bz,band1,band2,m1,m2,fac," q/=0"
!             if(dtset%prtvol>=10) write(6,'(a,6i4,e15.5,a)') "*****FAC*********",ik_bz,ikmq_bz,band1,band2,m1,m2,fac," q/=0"
             green_w=green_w*fac
           endif

           call assemblychi0_sym(ik_bz,nspinor,Ep,Ltg_q,green_w,Ep%npwepG0,rhotwg,Gsph_epsG0,chi0)

         CASE (1,2) ! Spectral method ! TODO Does not work with spinor
           call assemblychi0sf(ik_bz,nspinor,Ep%symchi,Ltg_q,Ep%npwepG0,Ep%npwe,rhotwg,Gsph_epsG0,&
&            deltaf_b1kmq_b2k,my_wl,iomegal,wl,my_wr,iomegar,wr,Ep%nomegasf,sf_chi0)

         CASE DEFAULT
           MSG_BUG("Wrong spmeth")
         END SELECT
         !
         ! === Accumulating the sum rule on chi0 ===
         ! *Eq. (5.284) in G.D. Mahan Many-Particle Physics 3rd edition.
         ! TODO Does not work with spinor

         factor=spin_fact*qp_occ(band2,ik_ibz,spin)
         call accumulate_chi0sumrule(ik_bz,Ep%symchi,Ep%npwe,factor,deltaeGW_b1kmq_b2k,&
&          Ltg_q,Gsph_epsG0,Ep%npwepG0,rhotwg,chi0_sumrule)
         !
         ! * Include also the completeness correction in the sum rule
         if (Ep%gwcomp==1) then
           factor=-spin_fact*qp_occ(band2,ik_ibz,spin)
           call accumulate_chi0sumrule(ik_bz,Ep%symchi,Ep%npwe,factor,deltaeGW_enhigh_b2k,&
&            Ltg_q,Gsph_epsG0,Ep%npwepG0,rhotwg,chi0_sumrule)
           if (band1==Ep%nbnds) then
             chi0_sumrule(:)=chi0_sumrule(:) + wtk_ltg(ik_bz)*spin_fact*qp_occ(band2,ik_ibz,spin)*deltaeGW_enhigh_b2k
           end if
         end if

       end do !band2
     end do !band1

     ABI_FREE(gw_gbound)
     if (Dtset%pawcross==1) then
       ABI_FREE(gboundf)
     end if
   end do !ik_bz
   !
   ! Deallocation of arrays private to the spin loop.
   ABI_FREE(igfftepsG0)
   ABI_FREE(ur1_kmq_ibz)
   ABI_FREE(ur2_k_ibz)
   ABI_FREE(usr1_kmq)
   ABI_FREE(ur2_k)
   ABI_FREE(rhotwg)
   ABI_FREE(tabr_k)
   ABI_FREE(tabr_kmq)
   if (allocated(green_w)) then
     ABI_FREE(green_w)
   end if
   if (allocated(wfwfg)) then
     ABI_FREE(wfwfg)
   end if
   if (allocated(green_enhigh_w)) then
     ABI_FREE(green_enhigh_w)
   end if
   if (Psps%usepaw==1) then
     call pawcprj_free(Cprj2_k  )
     ABI_DT_FREE(Cprj2_k)
     call pawcprj_free(Cprj1_kmq)
     ABI_DT_FREE(Cprj1_kmq)
     if (Dtset%pawcross==1) then
       ABI_FREE(ur_ae1)
       ABI_FREE(ur_ae_onsite1)
       ABI_FREE(ur_ps_onsite1)
       ABI_FREE(ur_ae2)
       ABI_FREE(ur_ae_onsite2)
       ABI_FREE(ur_ps_onsite2)
       ABI_FREE(tabrf_k)
       ABI_FREE(tabrf_kmq)
       ABI_FREE(gboundf)
       ABI_FREE(igfftepsG0f)
     end if
   end if
 end do !spin
 !
 ! === After big loop over transitions, now MPI ===
 ! * Master took care of the contribution in case of metallic|spin polarized systems.
 SELECT CASE (Ep%spmeth)
 CASE (0) ! Adler-Wiser
   ! * Collective sum of the contributions of each node.
   ! * Looping on frequencies to avoid problems with the size of the MPI packet
   do io=1,Ep%nomega
     call xmpi_sum(chi0(:,:,io),comm,ierr)
   end do
   !
 CASE (1,2) ! Spectral method.

   call hilbert_transform(Ep%npwe,Ep%nomega,Ep%nomegasf,my_wl,my_wr,kkweight,sf_chi0,chi0,Ep%spmeth)
   !
   !  Deallocate here before the xsums
   if (allocated(sf_chi0)) then
     ABI_FREE(sf_chi0)
   end if

   ! === Collective sum of the contributions ===
   ! * Looping over frequencies to avoid problems with the size of the MPI packet
   do io=1,Ep%nomega
     call xmpi_sum(chi0(:,:,io),comm,ierr)
   end do

 CASE DEFAULT
   MSG_BUG("Wrong spmeth")
 END SELECT
 !
 ! Divide by the volume
!$OMP PARALLEL WORKSHARE
   chi0=chi0*weight/Cryst%ucvol
!$OMP END PARALLEL WORKSHARE
 !
 ! === Collect the sum rule ===
 ! * The pi factor comes from Im[1/(x-ieta)] = pi delta(x)
 call xmpi_sum(chi0_sumrule,comm,ierr)
 chi0_sumrule=chi0_sumrule*pi*weight/Cryst%ucvol
 !
 ! *************************************************
 ! **** Now each node has chi0(q,G,Gp,Ep%omega) ****
 ! *************************************************

 ! Impose Hermiticity (valid only for zero or purely imaginary frequencies)
 ! MG what about metals, where we have poles around zero?
 do io=1,Ep%nomega
   if (ABS(REAL(Ep%omega(io)))<0.00001) then
     do ig2=1,Ep%npwe
       do ig1=1,ig2-1
         chi0(ig2,ig1,io) = GWPC_CONJG(chi0(ig1,ig2,io))
       end do
     end do
   end if
 end do
 !
 ! === Symmetrize chi0 in case of AFM system ===
 ! * Reconstruct $chi0{\down,\down}$ from $chi0{\up,\up}$.
 ! * Works only in case of magnetic group Shubnikov type IV.
 if (Cryst%use_antiferro) then
   call symmetrize_afm_chi0(Cryst,Gsph_epsG0,Ltg_q,Ep%npwe,Ep%nomega,chi0)
 end if
 !
 ! =====================
 ! ==== Free memory ====
 ! =====================
 ABI_FREE(bbp_ks_distrb)

 if (allocated(gw_gfft)) then
   ABI_FREE(gw_gfft)
 end if
 if (allocated(kkweight)) then
   ABI_FREE(kkweight)
 end if
 if (allocated(omegasf)) then
   ABI_FREE(omegasf)
 end if

 if (allocated(gspfft_igfft)) then
   ABI_FREE(gspfft_igfft)
 end if

 call gsph_free(Gsph_FFT)

 ! deallocation for PAW.
 if (Psps%usepaw==1) then
   call pawpwij_free(Pwij)
   ABI_DT_FREE(Pwij)
   if (allocated(Pwij_fft)) then
     call pawpwij_free(Pwij_fft)
     ABI_DT_FREE(Pwij_fft)
   end if
 end if

 if(dtset%ucrpa>=1) then
   ABI_DEALLOCATE(coeffW_BZ)
 endif

 call timab(331,2,tsec)
 call cwtime(cpu_time,wall_time,gflops,"stop")
 write(std_out,'(2(a,f9.1))')" cpu_time = ",cpu_time,", wall_time = ",wall_time

 DBG_EXIT("COLL")

end subroutine cchi0
!!***
