!{\src2tex{textfont=tt}}
!!****f* ABINIT/cchi0q0
!! NAME
!! cchi0q0
!!
!! FUNCTION
!! Calculate chi0 in the limit q-->0
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  use_tr=If .TRUE. Wfs_val are allocate and only resonant transitions are evaluated (assumes time reversal symmetry)
!!  Dtset <type(dataset_type)>=all input variables in this dataset
!!  Ep= datatype gathering differening parameters related to the calculation of the inverse dielectric matrix
!!  Gsph_epsG0<gvectors_data_type>: Info on the G-sphere used to describe chi0/espilon (including umklapp)
!!    %ng=number of G vectors
!!    %rottbm1(ng,2,nsym)=contains the index (IS^{-1}) G  in the array gvec
!!    %phmGt(ng,nsym)=phase factor e^{-iG.\tau} needed to symmetrize oscillator matrix elements and chi0
!!    %gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!!    %gprimd(3,3)=dimensional reciprocal space primitive translations (b^-1)
!!  Ep%inclvkb=flag to include (or not) the grad of Vkb
!!  Ltg_q= little group datatype
!!  nbvw=number of bands in the arrays wfrv,wfgv
!!  Kmesh<kmesh_t> The k-point mesh
!!   %kbz(3,nbz)=k-point coordinates, full Brillouin zone
!!   %tab(nbz)= table giving for each k-point in the BZ (kBZ), the corresponding
!!   irreducible point (kIBZ), where kBZ= (IS) kIBZ and I is either the inversion or the identity
!!   %tabi(nbzx)= for each point in the BZ defines whether inversion  has to be
!!   considered in the relation kBZ=(IS) kIBZ (1 => only S; -1 => -S)
!!   %tabo(nbzx)= the symmetry operation S that takes kIBZ to each kBZ
!!   %tabp(nbzx)= phase factor associated to tnons e^{-i 2 \pi k\cdot R{^-1}t}
!!  ktabr(nfftot_gw,Kmesh%nbz) index of R^-(r-t) in the FFT array, where k_BZ = (IS) k_IBZ and S = \transpose R^{-1}
!!  Ep%nbnds=number of bands
!!  ngfft_gw(18)= array containing all the information for 3D FFT for the oscillator strengths.
!!  Ep%nomega=number of frequencies
!!  Cryst<crystal_t>= data type gathering info on symmetries and unit cell
!!   %natom=number of atoms
!!   %nsym=number of symmetry operations
!!   %symrec(3,3,nsym)=symmetry operations in reciprocal space
!!   %typat(natom)=type of each atom
!!   %xred(3,natom)=reduced coordinated of atoms
!!   %rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!   %timrev=2 if time-reversal symmetry can be used, 1 otherwise
!!  Ep%npwe=number of planewaves for sigma exchange (input variable)
!!  nfftot_gw=Total number of points in the GW FFT grid
!!  Ep%omega(Ep%nomega)=frequencies
!!  Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!     %mpsang=1+maximum angular momentum for nonlocal pseudopotential
!!  Pawang<pawang_type> angular mesh discretization and related data:
!!  Pawrad(ntypat*usepaw)<Pawrad_type>=paw radial mesh and related data
!!  Paw_ij(natom*usepaw)<Paw_ij_type)>=paw arrays given on (i,j) channels
!!  QP_BSt<ebands_t>=Quasiparticle energies and occupations (for the moment real quantities)
!!    %mband=MAX number of bands over k-points and spin (==Ep%nbnds)
!!    %occ(mband,nkpt,nsppol)=QP occupation numbers, for each k point in IBZ, and each band
!!    %eig(mband,nkpt,nsppol)=GW energies, for self-consistency purposes
!!  KS_BSt<ebands_t>=KS energies and occupations.
!!    %eig(mband,nkpt,nsppol)=KS energies
!!  Paw_pwff<pawpwff_t>=Form factor used to calculate the onsite mat. elements of a plane wave.
!!  Wfd<wfd_t>=Object used to access the wavefunctions
!!
!! OUTPUT
!!  chi0(Ep%npwe,Ep%npwe,Ep%nomega)=independent-particle susceptibility matrix for wavevector qq,
!!   and frequencies defined by Ep%omega
!!  chi0_lwing(Ep%npwe*Ep%nI,Ep%nomega,3)= Lower wings
!!  chi0_uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)= Upper wings
!!  chi0_head(3,3,Ep%nomega)=Head of chi0.
!!
!! NOTES
!!  *) The terms "head", "wings" and "body" of chi(G,Gp) refer to
!!     G=Gp=0, either G or Gp=0, and neither=0 respectively
!!
!!  *) Symmetry conventions:
!!      1) symmetry in real space is defined as: R_t f(r) = f(R^-1(r-t))
!!      2) S=\transpose R^-1
!!      3) kbz=S kibz
!!
!!  The wavefunctions for the k-point in the BZ are (assuming nondegenerate states):
!!
!!  u(G,b, Sk) = u ( S^-1G,b,k)* e^{-i(Sk+G)*t)
!!  u(G,b,-Sk) = u*(-S^-1G,b,k)* e^{ i(Sk-G)*t)
!!
!!  u(r,b, Sk) = u (R^-1(r-t),b,k) e^{-iSk*t}
!!  u(r,b,-Sk) = u*(R^-1(r-t),b,k) e^{ iSK*t}
!!
!!  The gradient of Vnl(K,Kp) for the k-point in the BZ should be:
!!   gradvnl(SG,SGp,Sk)=S gradvnl(G,Gp,kibz)
!! /***********************************************************************/
!!
!! TODO
!!  Check npwepG0 before Switching on umklapp
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      accumulate_chi0_q0,accumulate_chi0sumrule,accumulate_sfchi0_q0
!!      approxdelta,calc_wfwfg,chi0_bbp_mask,completechi0_deltapart,cwtime
!!      flush_unit,get_bz_item,get_gftt,get_uug,gsph_fft_tabs,gsph_free
!!      gsph_in_fftbox,hilbert_transform,hilbert_transform_headwings
!!      littlegroup_print,make_transitions,paw_cross_ihr_comm
!!      paw_cross_rho_tw_g,paw_rho_tw_g,paw_symcprj,pawcprj_alloc,pawcprj_copy
!!      pawcprj_free,pawhur_free,pawhur_init,pawpwij_free,pawpwij_init
!!      print_gsphere,read_plowannier,rho_tw_g,setup_spectral
!!      symmetrize_afm_chi0,vkbr_free,vkbr_init,wfd_change_ngfft
!!      wfd_distribute_bbp,wfd_get_cprj,wfd_get_ur,wfd_paw_get_aeur,wfd_sym_ur
!!      wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine cchi0q0(use_tr,Dtset,Cryst,Ep,Psps,Kmesh,QP_BSt,KS_BSt,Gsph_epsG0,&
&  Pawang,Pawrad,Pawtab,Paw_ij,Paw_pwff,Pawfgrtab,Paw_onsite,ktabr,ktabrf,nbvw,ngfft_gw,&
&  nfftot_gw,ngfftf,nfftf_tot,chi0,chi0_head,chi0_lwing,chi0_uwing,Ltg_q,chi0_sumrule,Wfd,Wfdf)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_time

 use m_gwdefs,          only : GW_TOL_DOCC, GW_TOL_W0, czero_gw, em1params_t, g0g0w
 use m_numeric_tools,   only : imin_loc
 use m_geometry,        only : normv, vdotw
 use m_crystal,         only : crystal_t
 use m_fft_mesh,        only : get_gftt
 use m_bz_mesh,         only : kmesh_t, get_BZ_item, littlegroup_t, littlegroup_print
 use m_gsphere,         only : gsphere_t, gsph_fft_tabs, gsph_in_fftbox, gsph_free, print_gsphere
 use m_io_tools,        only : flush_unit
 use m_wfd,             only : wfd_get_ur, wfd_t, wfd_distribute_bbp, wfd_get_cprj, &
&                              wfd_barrier, wfd_change_ngfft,wfd_paw_get_aeur, wfd_sym_ur
 use m_oscillators,     only : rho_tw_g, calc_wfwfg
 use m_vkbr,            only : vkbr_t, vkbr_free, vkbr_init, nc_ihr_comm
 use m_chi0,            only : hilbert_transform, setup_spectral, symmetrize_afm_chi0, approxdelta, &
                               accumulate_chi0_q0, accumulate_sfchi0_q0, hilbert_transform_headwings, &
                               completechi0_deltapart, accumulate_chi0sumrule, make_transitions, chi0_bbp_mask
 use m_pawang,          only : pawang_type
 use m_pawrad,          only : pawrad_type
 use m_pawtab,          only : pawtab_type
 use m_paw_ij,          only : paw_ij_type
 use m_pawfgrtab,       only : pawfgrtab_type
 use m_pawcprj,         only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy
 use m_pawpwij,        only : pawpwff_t, pawpwij_t, pawpwij_init, pawpwij_free, paw_rho_tw_g, paw_cross_rho_tw_g
 use m_paw_pwaves_lmn,  only : paw_pwaves_lmn_t
 use m_pawhr,           only : pawhur_t, pawhur_free, pawhur_init, paw_ihr, paw_cross_ihr_comm

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cchi0q0'
 use interfaces_14_hidewrite
 use interfaces_65_paw
 use interfaces_70_gw, except_this_one => cchi0q0
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbvw,nfftot_gw,nfftf_tot
 logical,intent(in) :: use_tr
 type(ebands_t),target,intent(in) :: QP_BSt,KS_BSt
 type(crystal_t),intent(in) :: Cryst
 type(Dataset_type),intent(in) :: Dtset
 type(littlegroup_t),intent(in) :: Ltg_q
 type(em1params_t),intent(in) :: Ep
 type(kmesh_t),intent(in) :: Kmesh
 type(gsphere_t),intent(in) :: Gsph_epsG0
 type(Pseudopotential_type),intent(in) :: Psps
 type(Pawang_type),intent(in) :: Pawang
 type(wfd_t),target,intent(inout) :: Wfd,Wfdf
!arrays
 integer,intent(in) :: ktabr(nfftot_gw,Kmesh%nbz),ktabrf(nfftf_tot*Dtset%pawcross,Kmesh%nbz)
 integer,intent(in) :: ngfft_gw(18),ngfftf(18)
 real(dp),intent(out) :: chi0_sumrule(Ep%npwe)
 complex(gwpc),intent(out) :: chi0(Ep%npwe,Ep%npwe,Ep%nomega)
 complex(dpc),intent(out) :: chi0_lwing(Ep%npwe*Ep%nI,Ep%nomega,3)
 complex(dpc),intent(out) :: chi0_uwing(Ep%npwe*Ep%nJ,Ep%nomega,3)
 complex(dpc),intent(out) :: chi0_head(3,3,Ep%nomega)
 type(Pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 type(Paw_ij_type),intent(in) :: Paw_ij(Cryst%natom*Psps%usepaw)
 type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
 type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp=1,enough=10,two_poles=2,one_pole=1,ndat1=1
 integer :: bandinf,bandsup,lcor,nspinor,npw_k,istwf_k,mband,nfft
 integer :: band1,band2,iat,ig,itim_k,ik_bz,ik_ibz,io,iqlwl,ispinor1,ispinor2,isym_k
 integer :: itypatcor,m1,m2,nkpt_summed,dim_rtwg,use_padfft,gw_fftalga,use_padfftf,mgfftf
 integer :: my_nbbp,my_nbbpks,spin,nsppol !ig1,ig2,
 integer :: comm,ierr,my_wl,my_wr,iomegal,iomegar,gw_mgfft,dummy
 real(dp) :: cpu_time,wall_time,gflops
 real(dp) :: fac,fac1,fac2,fac3,fac4,spin_fact,deltaf_b1b2,weight,factor
 real(dp) :: max_rest,min_rest,my_max_rest,my_min_rest
 real(dp) :: en_high,deltaeGW_enhigh_b2
 real(dp) :: wl,wr,numerator,deltaeGW_b1b2
 real(dp) :: gw_gsq,memreq
 complex(dpc) :: deltaeKS_b1b2
 logical :: qzero,luwindow
 character(len=500) :: msg_tmp,msg,allup
 type(gsphere_t) :: Gsph_FFT
!arrays
 integer,ABI_CONTIGUOUS pointer :: kg_k(:,:)
 integer :: ucrpa_bands(2)
 integer :: wtk_ltg(Kmesh%nbz)
 integer :: got(Wfd%nproc)
 integer,allocatable :: tabr_k(:),tabrf_k(:)
 integer,allocatable :: igffteps0(:),gspfft_igfft(:),igfftepsG0f(:)
 integer,allocatable :: gw_gfft(:,:),gw_gbound(:,:),dummy_gbound(:,:),gboundf(:,:)
 integer,allocatable :: bbp_ks_distrb(:,:,:,:)
 real(dp) :: kbz(3),spinrot_kbz(4),q0(3)
 real(dp),ABI_CONTIGUOUS pointer :: ks_energy(:,:,:),qp_energy(:,:,:),qp_occ(:,:,:)
 real(dp),allocatable :: omegasf(:)
 complex(gwpc) :: rhotwx(3,Wfd%nspinor**2)
 complex(gwpc),allocatable :: rhotwg(:)
 complex(dpc),allocatable :: green_w(:),green_enhigh_w(:)
 complex(dpc),allocatable :: sf_lwing(:,:,:),sf_uwing(:,:,:),sf_head(:,:,:)
 complex(dpc) :: wng(3),chq(3)
 complex(dpc) :: ph_mkt
 complex(gwpc),allocatable :: sf_chi0(:,:,:)
 complex(dpc),allocatable :: kkweight(:,:)
 complex(gwpc),allocatable :: ur1_kibz(:),ur2_kibz(:)
 complex(gwpc),allocatable :: usr1_k(:),ur2_k(:)
 complex(gwpc),allocatable :: wfwfg(:)
 complex(gwpc),allocatable :: ur_ae1(:),ur_ae_onsite1(:),ur_ps_onsite1(:)
 complex(gwpc),allocatable :: ur_ae2(:),ur_ae_onsite2(:),ur_ps_onsite2(:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: ug1(:),ug2(:)
 complex(dpc), allocatable :: coeffW_BZ(:,:,:,:,:,:)
 logical :: gradk_not_done(Kmesh%nibz)
 logical,allocatable :: bbp_mask(:,:)
 type(pawcprj_type),allocatable :: Cprj1_bz(:,:),Cprj2_bz(:,:)
 type(pawcprj_type),allocatable :: Cprj1_ibz(:,:),Cprj2_ibz(:,:)
 type(pawpwij_t),allocatable :: Pwij(:),Pwij_fft(:)
 type(pawhur_t),allocatable :: Hur(:)
 type(vkbr_t),allocatable :: vkbr(:)
!************************************************************************

 DBG_ENTER("COLL")
 call cwtime(cpu_time,wall_time,gflops,"start")

 ! Change FFT mesh if needed
 if (ANY(ngfft_gw(1:3) /= Wfd%ngfft(1:3))) call wfd_change_ngfft(Wfd,Cryst,Psps,ngfft_gw)

 gw_mgfft = MAXVAL(ngfft_gw(1:3)); gw_fftalga = ngfft_gw(7)/100 !; gw_fftalgc=MOD(ngfft_gw(7),10)
 if (Dtset%pawcross==1) mgfftf = MAXVAL(ngfftf(1:3))

 ! Copy important variables.
 comm = Wfd%comm; nsppol = Wfd%nsppol; nspinor = Wfd%nspinor; mband = Wfd%mband
 nfft = Wfd%nfft
 ABI_CHECK(Wfd%nfftot==nfftot_gw,"Wrong nfftot_gw")
 dim_rtwg=1; if (nspinor==2) dim_rtwg=2 ! Can reduce size depending on Ep%nI and Ep%nj

 ucrpa_bands(1)=dtset%ucrpa_bands(1)
 ucrpa_bands(2)=dtset%ucrpa_bands(2)
 luwindow=.false.
 if(abs(dtset%ucrpa_window(1)+1_dp)>tol8.and.(abs(dtset%ucrpa_window(2)+1_dp)>tol8)) then
   luwindow=.true.
 endif

 ! For cRPA calculation of U: read forlb.ovlp
 if(dtset%ucrpa>=1) then
   call read_plowannier(Cryst,bandinf,bandsup,coeffW_BZ,itypatcor,Kmesh,lcor,luwindow,&
& nspinor,nsppol,pawang,dtset%prtvol,ucrpa_bands)
 endif

 ks_energy => KS_BSt%eig
 qp_energy => QP_BSt%eig; qp_occ => QP_BSt%occ

 chi0_lwing = czero; chi0_uwing= czero; chi0_head = czero

 if (Psps%usepaw==0) then
   if (Ep%inclvkb/=0) then
     ! Include the term <n,k|[Vnl,iqr]|n"k>' for q->0.
     ABI_CHECK(nspinor==1,"nspinor+inclvkb not coded")
   else
     MSG_WARNING('Neglecting <n,k|[Vnl,iqr]|m,k>')
   end if
 else
   ! For PAW+LDA+U, precalculate <\phi_i|[Hu,r]|phi_j\>
   ABI_DT_MALLOC(HUr,(Cryst%natom))
   if (Dtset%usepawu/=0) then
     call pawhur_init(hur,nsppol,Dtset%pawprtvol,Cryst,Pawtab,Pawang,Pawrad,Paw_ij)
   end if
 end if

 ! Initialize the completeness correction.
 ABI_MALLOC(green_enhigh_w,(Ep%nomega))
 green_enhigh_w=czero

 if (Ep%gwcomp==1) then
   en_high=MAXVAL(qp_energy(Ep%nbnds,:,:))+Ep%gwencomp
   write(msg,'(a,f8.2,a)')' Using completeness correction with energy ',en_high*Ha_eV,' [eV] '
   call wrtout(std_out,msg,'COLL')
   ABI_MALLOC(wfwfg,(nfft*nspinor**2))

   ! Init the largest G-sphere contained in the FFT box for the wavefunctions.
   call gsph_in_fftbox(Gsph_FFT,Cryst,Wfd%ngfft)
   call print_gsphere(Gsph_FFT,unit=std_out,prtvol=10)

   ABI_MALLOC(gspfft_igfft,(Gsph_FFT%ng))
   ABI_MALLOC(dummy_gbound,(2*gw_mgfft+8,2))

   ! Mapping between G-sphere and FFT box.
   call gsph_fft_tabs(Gsph_FFT, [0, 0, 0],Wfd%mgfft,Wfd%ngfft,dummy,dummy_gbound,gspfft_igfft)
   ABI_FREE(dummy_gbound)

   if (Psps%usepaw==1) then
     ! Prepare the onsite contributions on the GW FFT mesh.
     ABI_MALLOC(gw_gfft,(3,nfft))
     q0=zero
     call get_gftt(ngfft_gw,q0,Cryst%gmet,gw_gsq,gw_gfft) ! The set of plane waves in the FFT Box.
     ABI_DT_MALLOC(Pwij_fft,(Psps%ntypat))
     call pawpwij_init(Pwij_fft,nfft,(/zero,zero,zero/),gw_gfft,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
   end if
 end if

 ! Setup weight (2 for spin unpolarized systems, 1 for polarized).
 ! spin_fact is used to normalize the occupation factors to one.
 ! Consider also the AFM case.
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

 ! Weight for points in the IBZ_q
 wtk_ltg(:)=1
 if (Ep%symchi==1) then
   do ik_bz=1,Ltg_q%nbz
     wtk_ltg(ik_bz)=0
     if (Ltg_q%ibzq(ik_bz)/=1) CYCLE ! Only k in IBZ_q
     wtk_ltg(ik_bz)=SUM(Ltg_q%wtksym(:,:,ik_bz))
   end do
 end if

 write(msg,'(a,i3,a)')' Q-points for long wave-length limit. # ',Ep%nqlwl,ch10
 do iqlwl=1,Ep%nqlwl
   write(msg_tmp,'(1x,i5,a,2x,3f12.6,a)') iqlwl,')',Ep%qlwl(:,iqlwl),ch10
   msg=TRIM(msg)//msg_tmp
 end do
 call wrtout(std_out,msg,'COLL')

 write(msg,'(a,i2,2a,i2)')&
&  ' Using spectral method for the imaginary part = ',Ep%spmeth,ch10,&
&  ' Using symmetries to sum only over the IBZ_q  = ',Ep%symchi
 call wrtout(std_out,msg,'COLL')

 if (use_tr) then
   ! Special care has to be taken in metals and/or spin dependent systems
   ! as Wfs_val might contain unoccupied states.
   call wrtout(std_out,' Using faster algorithm based on time reversal symmetry.')
 else
   call wrtout(std_out,' Using slow algorithm without time reversal symmetry.')
 end if

 ! Evaluate oscillator matrix elements btw partial waves. Note q=Gamma
 if (Psps%usepaw==1) then
   ABI_DT_MALLOC(Pwij,(Psps%ntypat))
   call pawpwij_init(Pwij,Ep%npwepG0,(/zero,zero,zero/),Gsph_epsG0%gvec,Cryst%rprimd,Psps,Pawtab,Paw_pwff)

   ABI_DT_MALLOC(Cprj1_bz,(Cryst%natom,nspinor))
   call pawcprj_alloc(Cprj1_bz,0,Wfd%nlmn_atm)
   ABI_DT_MALLOC(Cprj2_bz,(Cryst%natom,nspinor))
   call pawcprj_alloc(Cprj2_bz,0,Wfd%nlmn_atm)

   ABI_DT_MALLOC(Cprj1_ibz,(Cryst%natom,nspinor))
   call pawcprj_alloc(Cprj1_ibz,0,Wfd%nlmn_atm)
   ABI_DT_MALLOC(Cprj2_ibz,(Cryst%natom,nspinor))
   call pawcprj_alloc(Cprj2_ibz,0,Wfd%nlmn_atm)
   if (Dtset%pawcross==1) then
     ABI_MALLOC(ur_ae1,(nfftf_tot*nspinor))
     ABI_MALLOC(ur_ae_onsite1,(nfftf_tot*nspinor))
     ABI_MALLOC(ur_ps_onsite1,(nfftf_tot*nspinor))
     ABI_MALLOC(ur_ae2,(nfftf_tot*nspinor))
     ABI_MALLOC(ur_ae_onsite2,(nfftf_tot*nspinor))
     ABI_MALLOC(ur_ps_onsite2,(nfftf_tot*nspinor))
     ABI_MALLOC(igfftepsG0f,(Ep%npwepG0))
     ABI_MALLOC(tabrf_k,(nfftf_tot))
   end if
 end if

 ABI_MALLOC(rhotwg,(Ep%npwe*nspinor**2))
 ABI_MALLOC(tabr_k,(nfft))
 ABI_MALLOC(ur1_kibz,(nfft*nspinor))
 ABI_MALLOC(ur2_kibz,(nfft*nspinor))

 ABI_MALLOC(usr1_k,(nfft*nspinor))
 ABI_MALLOC(ur2_k,(nfft*nspinor))
 !
 ! Tables for the FFT of the oscillators.
 !  a) FFT index of the G sphere (only vertical transitions, unlike cchi0, no need to shift the sphere).
 !  b) gw_gbound table for the zero-padded FFT performed in rhotwg.
 ABI_MALLOC(igffteps0,(Gsph_epsG0%ng))
 ABI_MALLOC(gw_gbound,(2*gw_mgfft+8,2))
 call gsph_fft_tabs(Gsph_epsG0,(/0,0,0/),gw_mgfft,ngfft_gw,use_padfft,gw_gbound,igffteps0)
 if ( ANY(gw_fftalga == (/2,4/)) ) use_padfft=0 ! Pad-FFT is not coded in rho_tw_g
#ifdef FC_IBM
 ! XLF does not deserve this optimization (problem with [v67mbpt][t03])
 use_padfft = 0
#endif
 if (use_padfft==0) then
   ABI_FREE(gw_gbound)
   ABI_MALLOC(gw_gbound,(2*gw_mgfft+8,2*use_padfft))
 end if
 if (Dtset%pawcross==1) then
    ABI_MALLOC(gboundf,(2*mgfftf+8,2))
   call gsph_fft_tabs(Gsph_epsG0,(/0,0,0/),mgfftf,ngfftf,use_padfftf,gboundf,igfftepsG0f)
   if ( ANY(gw_fftalga == (/2,4/)) ) use_padfftf=0
   if (use_padfftf==0) then
     ABI_FREE(gboundf)
     ABI_MALLOC(gboundf,(2*mgfftf+8,2*use_padfftf))
   end if
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
     ik_ibz=Kmesh%tab(ik_bz)

     call chi0_bbp_mask(Ep,use_tr,QP_BSt,mband,ik_ibz,ik_ibz,spin,spin_fact,bbp_mask)

     call wfd_distribute_bbp(Wfd,ik_ibz,spin,allup,my_nbbp,bbp_ks_distrb(:,:,ik_bz,spin),got,bbp_mask)
     my_nbbpks = my_nbbpks + my_nbbp
   end do
 end do

 ABI_FREE(bbp_mask)

 write(msg,'(a,i0,a)')" Will sum ",my_nbbpks," (b,b',k,s) states in chi0q0."
 call wrtout(std_out,msg,'PERS')

 SELECT CASE (Ep%spmeth)
 CASE (0)
   call wrtout(std_out,' Calculating chi0(q=(0,0,0),omega,G,G")',"COLL")
   ABI_MALLOC(green_w,(Ep%nomega))

 CASE (1, 2)
   call wrtout(std_out,' Calculating Im chi0(q=(0,0,0),omega,G,G")','COLL')
   !
   ! === Find max and min resonant transitions for this q, report values for this processor ===
   call make_transitions(Wfd,1,Ep%nbnds,nbvw,nsppol,Ep%symchi,Cryst%timrev,GW_TOL_DOCC,&
&    max_rest,min_rest,my_max_rest,my_min_rest,Kmesh,Ltg_q,qp_energy,qp_occ,(/zero,zero,zero/),bbp_ks_distrb)

   !FIXME there is a problem in make_transitions due to MPI_enreg
   !ltest=(MPI_enreg%gwpara==0.or.MPI_enreg%gwpara==2)
   !ABI_CHECK(ltest,"spectral method with gwpara==1 not coded")
   !
   ! === Calculate frequency dependent weights for Kramers Kronig transform ===
   ABI_MALLOC(omegasf,(Ep%nomegasf))
   ABI_MALLOC(kkweight,(Ep%nomegasf,Ep%nomega))
   !my_wl=1; my_wr=Ep%nomegasf
   call setup_spectral(Ep%nomega,Ep%omega,Ep%nomegasf,omegasf,max_rest,min_rest,my_max_rest,my_min_rest,&
&     0,Ep%zcut,zero,my_wl,my_wr,kkweight)

   if (.not.use_tr) then
     MSG_BUG('Hilbert transform requires time-reversal')
   end if

   ! allocate heads and wings of the spectral function.
   ABI_MALLOC(sf_head,(3,3,my_wl:my_wr))
   ABI_MALLOC(sf_lwing,(Ep%npwe,my_wl:my_wr,3))
   ABI_MALLOC(sf_uwing,(Ep%npwe,my_wl:my_wr,3))
   sf_head=czero; sf_lwing=czero; sf_uwing=czero

   memreq = two*gwpc*Ep%npwe**2*(my_wr-my_wl+1)*b2Gb
   write(msg,'(a,f10.4,a)')' memory required per spectral point: ',two*gwpc*Ep%npwe**2*b2Mb,' [Mb]'
   call wrtout(std_out,msg,'PERS')
   write(msg,'(a,f10.4,a)')' memory required by sf_chi0q0:       ',memreq,' [Gb]'
   call wrtout(std_out,msg,'PERS')
   if (memreq > two) then
     MSG_WARNING(' Memory required for sf_chi0q0 is larger than 2.0 Gb!')
   end if
   ABI_STAT_MALLOC(sf_chi0,(Ep%npwe,Ep%npwe,my_wl:my_wr), ierr)
   ABI_CHECK(ierr==0, 'out-of-memory in sf_chi0q0')
   sf_chi0=czero_gw

 CASE DEFAULT
   MSG_BUG("Wrong spmeth")
 END SELECT

 nkpt_summed=Kmesh%nbz
 if (Ep%symchi/=0) then
   nkpt_summed=Ltg_q%nibz_ltg
   call littlegroup_print(Ltg_q,std_out,Dtset%prtvol,'COLL')
 end if

 ABI_DT_MALLOC(vkbr,(Kmesh%nibz))
 gradk_not_done=.TRUE.

 write(msg,'(a,i6,a)')' Calculation status ( ',nkpt_summed,' to be completed):'
 call wrtout(std_out,msg,'COLL')
 !
 ! ============================================
 ! === Begin big fat loop over transitions ====
 ! ============================================
 chi0=czero_gw; chi0_sumrule =zero

 ! Loop on spin to calculate $\chi_{\up,\up} + \chi_{\down,\down}$
 do spin=1,nsppol
   if (ALL(bbp_ks_distrb(:,:,:,spin) /= Wfd%my_rank)) CYCLE

   ! Loop over k-points in the BZ.
   do ik_bz=1,Kmesh%nbz
     if (Ep%symchi==1) then
       if (Ltg_q%ibzq(ik_bz)/=1) CYCLE ! Only IBZ_q
     end if

     if (ALL(bbp_ks_distrb(:,:,ik_bz,spin) /= Wfd%my_rank)) CYCLE

     write(msg,'(2(a,i4),a,i2,a,i3)')' ik= ',ik_bz,'/',Kmesh%nbz,' spin=',spin,' done by mpi rank:',Wfd%my_rank
     call wrtout(std_out,msg,'PERS')

     ! Get ik_ibz, non-symmorphic phase and symmetries from ik_bz.
     call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym_k,itim_k,ph_mkt)
     tabr_k=ktabr(:,ik_bz) ! Table for rotated FFT points
     spinrot_kbz(:)=Cryst%spinrot(:,isym_k)
     if (Dtset%pawcross==1) tabrf_k(:) = ktabrf(:,ik_bz)

     istwf_k =  Wfd%istwfk(ik_ibz)
     npw_k   =  Wfd%npwarr(ik_ibz)
     kg_k    => Wfd%Kdata(ik_ibz)%kg_k

     if (Psps%usepaw==0.and.Ep%inclvkb/=0.and.gradk_not_done(ik_ibz)) then
       ! Include term <n,k|[Vnl,iqr]|n"k>' for q->0.
       call vkbr_init(vkbr(ik_ibz),Cryst,Psps,Ep%inclvkb,istwf_k,npw_k,Kmesh%ibz(:,ik_ibz),kg_k)
       gradk_not_done(ik_ibz)=.FALSE.
     end if

     ! Loop over "conduction" states.
     do band1=1,Ep%nbnds
       if (ALL(bbp_ks_distrb(band1,:,ik_bz,spin) /= Wfd%my_rank)) CYCLE

       ug1 => Wfd%Wave(band1,ik_ibz,spin)%ug
       call wfd_get_ur(Wfd,band1,ik_ibz,spin,ur1_kibz)

       if (Psps%usepaw==1) then
         call wfd_get_cprj(Wfd,band1,ik_ibz,spin,Cryst,Cprj1_ibz,sorted=.FALSE.)
         call pawcprj_copy(Cprj1_ibz,Cprj1_bz)
         call paw_symcprj(ik_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,Cprj1_bz)
         if (Dtset%pawcross==1) then
           call wfd_paw_get_aeur(Wfdf,band1,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae1,ur_ae_onsite1,ur_ps_onsite1)
         end if
       end if

       do band2=1,Ep%nbnds ! Loop over "valence" states.
!       write(6,*) "ik,band1,band2",ik_bz,band1,band2

!         -----------------  cRPA for U
!debug         if (.not.luwindow.AND.dtset%ucrpa==1.AND.band1<=ucrpa_bands(2).AND.band1>=ucrpa_bands(1)&
!debug&                                            .AND.band2<=ucrpa_bands(2).AND.band2>=ucrpa_bands(1)) CYCLE
         if (luwindow.AND.dtset%ucrpa==1.AND.((KS_Bst%eig(band1,ik_ibz,spin)-KS_Bst%fermie)<=dtset%ucrpa_window(2)) &
&                                       .AND.((KS_Bst%eig(band1,ik_ibz,spin)-KS_Bst%fermie)>=dtset%ucrpa_window(1)) &
&                                       .AND.((KS_Bst%eig(band2,ik_ibz,spin)-KS_Bst%fermie)<=dtset%ucrpa_window(2)) &
&                                       .AND.((KS_Bst%eig(band2,ik_ibz,spin)-KS_Bst%fermie)>=dtset%ucrpa_window(1))) CYCLE
!         -----------------  cRPA for U

         if (bbp_ks_distrb(band1,band2,ik_bz,spin) /= Wfd%my_rank) CYCLE

         deltaf_b1b2  =spin_fact*(qp_occ(band1,ik_ibz,spin)-qp_occ(band2,ik_ibz,spin))
         deltaeKS_b1b2= ks_energy(band1,ik_ibz,spin) - ks_energy(band2,ik_ibz,spin)
         deltaeGW_b1b2= qp_energy(band1,ik_ibz,spin) - qp_energy(band2,ik_ibz,spin)

         if (Ep%gwcomp==0) then ! Skip negligible transitions.
           if (ABS(deltaf_b1b2) < GW_TOL_DOCC) CYCLE
         else
           ! when the completeness trick is used, we need to also consider transitions with vanishing deltaf
           !Rangel Correction for metals
           !if (qp_occ(band2,ik_ibz,spin) < GW_TOL_DOCC) CYCLE
           if (qp_occ(band2,ik_ibz,spin) < GW_TOL_DOCC .and. ( ABS(deltaf_b1b2)< GW_TOL_DOCC .or. band1<band2)) CYCLE
         end if

         ug2 => Wfd%Wave(band2,ik_ibz,spin)%ug
         call wfd_get_ur(Wfd,band2,ik_ibz,spin,ur2_kibz)

         if (Psps%usepaw==1) then
           call wfd_get_cprj(Wfd,band2,ik_ibz,spin,Cryst,Cprj2_ibz,sorted=.FALSE.)
           call pawcprj_copy(Cprj2_ibz,Cprj2_bz)
           call paw_symcprj(ik_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,Cprj2_bz)
           if (Dtset%pawcross==1) then
             call wfd_paw_get_aeur(Wfdf,band2,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur_ae2,ur_ae_onsite2,ur_ps_onsite2)
           end if
         end if

         SELECT CASE (Ep%spmeth)
         CASE (0)
           ! Adler-Wiser expression.
           ! Add small imaginary of the Time-Ordered resp function but only for non-zero real omega  FIXME What about metals?

           if (.not.use_tr) then
             ! Adler-Wiser without time-reversal.
             do io=1,Ep%nomega
               green_w(io) = g0g0w(Ep%omega(io),deltaf_b1b2,deltaeGW_b1b2,Ep%zcut,GW_TOL_W0,one_pole)
             end do
           else
             if (Ep%gwcomp==0) then ! cannot be completely skipped in case of completeness correction
               if (band1<band2) CYCLE ! Here we GAIN a factor ~2
             end if

             do io=1,Ep%nomega
               !Rangel: In metals, the intra-band transitions term does not contain the antiresonant part
               !if(abs(deltaeGW_b1b2)>GW_TOL_W0) green_w(io) = g0g0w(Ep%omega(io),deltaf_b1b2,deltaeGW_b1b2,Ep%zcut,GW_TOL_W0)
               if (band1==band2) green_w(io) = g0g0w(Ep%omega(io),deltaf_b1b2,deltaeGW_b1b2,Ep%zcut,GW_TOL_W0,one_pole)
               if (band1/=band2) green_w(io) = g0g0w(Ep%omega(io),deltaf_b1b2,deltaeGW_b1b2,Ep%zcut,GW_TOL_W0,two_poles)

               if (Ep%gwcomp==1) then ! Calculate the completeness correction
                 numerator= -spin_fact*qp_occ(band2,ik_ibz,spin)
                 deltaeGW_enhigh_b2=en_high-qp_energy(band2,ik_ibz,spin)
                 ! Completeness correction is NOT valid for real frequencies
                 if (REAL(Ep%omega(io))<GW_TOL_W0) then
                   green_enhigh_w(io) = g0g0w(Ep%omega(io),numerator,deltaeGW_enhigh_b2,Ep%zcut,GW_TOL_W0,two_poles)
                 else
                   green_enhigh_w(io) = czero_gw
                 endif
                 !
                 !Rangel Correction for metals
                 !if (deltaf_b1b2<0.d0) then
                 if (band1>=band2 .and. abs(deltaf_b1b2) > GW_TOL_DOCC) then
                   green_w(io)= green_w(io) - green_enhigh_w(io)
                 else ! Disregard green_w, since it is already accounted for through the time-reversal
                   green_w(io)=             - green_enhigh_w(io)
                 end if
               end if !gwcomp==1
             end do !io

             if (Ep%gwcomp==1.and.band1==band2) then
               ! Add the "delta part", symmetrization is done inside the routine.
               call calc_wfwfg(tabr_k,itim_k,spinrot_kbz,nfft,nspinor,ngfft_gw,ur2_kibz,ur2_kibz,wfwfg)

               if (Psps%usepaw==1) then
                 call paw_rho_tw_g(nfft,dim_rtwg,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,gw_gfft,&
&                  Cprj2_bz,Cprj2_bz,Pwij_fft,wfwfg)

                ! Add PAW cross term
                if (Dtset%pawcross==1) then
                  call paw_cross_rho_tw_g(nspinor,Ep%npwepG0,nfftf_tot,ngfftf,1,use_padfftf,igfftepsG0f,gboundf,&
&                  ur_ae2,ur_ae_onsite2,ur_ps_onsite2,itim_k,tabrf_k,ph_mkt,spinrot_kbz,&
&                  ur_ae2,ur_ae_onsite2,ur_ps_onsite2,itim_k,tabrf_k,ph_mkt,spinrot_kbz,&
&                  dim_rtwg,wfwfg)
                end if
               end if

               qzero=.TRUE.
               call completechi0_deltapart(ik_bz,qzero,Ep%symchi,Ep%npwe,Gsph_FFT%ng,Ep%nomega,nspinor,&
&                nfft,ngfft_gw,gspfft_igfft,Gsph_FFT,Ltg_q,green_enhigh_w,wfwfg,chi0)
             end if
           end if ! use_tr

         CASE (1, 2)
           ! Spectral method, here time-reversal is always assumed.
           if (deltaeGW_b1b2<0) CYCLE
           call approxdelta(Ep%nomegasf,omegasf,deltaeGW_b1b2,Ep%spsmear,iomegal,iomegar,wl,wr,Ep%spmeth)
         END SELECT

         ! FFT of u^*_{b1,k}(r) u_{b2,k}(r) and (q,G=0) limit using small q and k.p perturbation theory
         call rho_tw_g(nspinor,Ep%npwe,nfft,ndat1,ngfft_gw,1,use_padfft,igffteps0,gw_gbound,&
&          ur1_kibz,itim_k,tabr_k,ph_mkt,spinrot_kbz,&
&          ur2_kibz,itim_k,tabr_k,ph_mkt,spinrot_kbz,&
&          dim_rtwg,rhotwg)

         if (Psps%usepaw==0) then
           ! Matrix elements of i[H,r] for NC pseudopotentials.
           rhotwx = nc_ihr_comm(vkbr(ik_ibz),cryst,psps,npw_k,nspinor,istwf_k,Ep%inclvkb,Kmesh%ibz(:,ik_ibz),ug1,ug2,kg_k)

         else
           ! 1) Add PAW onsite contribution, projectors are already in the BZ.
           call paw_rho_tw_g(Ep%npwe,dim_rtwg,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_epsG0%gvec,&
&            Cprj1_bz,Cprj2_bz,Pwij,rhotwg)

           ! 2) Matrix elements of i[H,r] for PAW.
           rhotwx = paw_ihr(spin,nspinor,npw_k,istwf_k,Kmesh%ibz(:,ik_ibz),Cryst,Pawtab,ug1,ug2,kg_k,Cprj1_ibz,Cprj2_ibz,HUr)

           ! Add PAW cross term
           if (Dtset%pawcross==1) then
             call paw_cross_rho_tw_g(nspinor,Ep%npwepG0,nfftf_tot,ngfftf,1,use_padfftf,igfftepsG0f,gboundf,&
&             ur_ae1,ur_ae_onsite1,ur_ps_onsite1,itim_k,tabrf_k,ph_mkt,spinrot_kbz,&
&             ur_ae2,ur_ae_onsite2,ur_ps_onsite2,itim_k,tabrf_k,ph_mkt,spinrot_kbz,&
&             dim_rtwg,rhotwg)
              ! Add cross-term contribution to the commutator
             if (Dtset%userib/=111) then
               call paw_cross_ihr_comm(rhotwx,nspinor,nfftf_tot,Cryst,Pawfgrtab,Paw_onsite,&
&                   ur_ae1,ur_ae2,ur_ae_onsite1,ur_ae_onsite2,Cprj1_ibz,Cprj2_ibz)
             end if
           end if
         end if

         ! Treat a possible degeneracy between v and c.
         if (ABS(deltaeKS_b1b2)>GW_TOL_W0) then
           rhotwx=-rhotwx/deltaeKS_b1b2
         else
           rhotwx=czero_gw
         end if

         SELECT CASE (Ep%spmeth)
         CASE (0)

!          ---------------- Ucrpa (begin)
           if(dtset%ucrpa>=1.and..not.luwindow)  then
             fac=one
             fac1=zero
             fac2=zero
             fac3=zero
             fac4=one
             m1=-1
             m2=-1
             if(dtset%ucrpa<=2) then
               call flush_unit(std_out)
               if (       band1<=ucrpa_bands(2).AND.band1>=ucrpa_bands(1)&
&                    .AND.band2<=ucrpa_bands(2).AND.band2>=ucrpa_bands(1)) then
!               if(dtset%prtvol>=10)write(6,*)"calculation is in progress",band1,band2,ucrpa_bands(1),ucrpa_bands(2)
                 do iat=1, cryst%nattyp(itypatcor)
                   do ispinor1=1,nspinor
                     do m1=1,2*lcor+1
                       fac1=fac1+ real(coeffW_BZ(iat,spin,band1,ik_bz,ispinor1,m1)*&
&                                 conjg(coeffW_BZ(iat,spin,band1,ik_bz,ispinor1,m1)))
                       fac2=fac2+ real(coeffW_BZ(iat,spin,band2,ik_bz,ispinor1,m1)*&
&                                 conjg(coeffW_BZ(iat,spin,band2,ik_bz,ispinor1,m1)))
                       do ispinor2=1,nspinor
                         do m2=1,2*lcor+1
                           fac=fac - real(coeffW_BZ(iat,spin,band1,ik_bz,ispinor1,m1)*&
&                                    conjg(coeffW_BZ(iat,spin,band1,ik_bz,ispinor1,m1))* &
&                                    coeffW_BZ(iat,spin,band2,ik_bz,ispinor2,m2)*&
&                                    conjg(coeffW_BZ(iat,spin,band2,ik_bz,ispinor2,m2)))
                           fac3=fac3 + real(coeffW_BZ(iat,spin,band1,ik_bz,ispinor1,m1)*&
&                                      conjg(coeffW_BZ(iat,spin,band1,ik_bz,ispinor1,m1))* &
&                                      coeffW_BZ(iat,spin,band2,ik_bz,ispinor2,m2)*&
&                                      conjg(coeffW_BZ(iat,spin,band2,ik_bz,ispinor2,m2)))
!                         if(dtset%prtvol>=10)write(6,*) fac,fac3
                         enddo
                       enddo
!                       if(dtset%prtvol>=10)write(6,*) fac,fac3,fac1,fac2,fac1*fac2
                     enddo
                   enddo
                 enddo
                 fac4=fac
!                 fac=zero
                 if(dtset%ucrpa==1) fac=zero
!                 write(6,'(a,6i4,e15.5,a)') "*****FAC*********",ik_bz,ik_bz,band1,band2,m1,m2,fac," q==0"
               endif
             else if (dtset%ucrpa==3) then
               if        (band1<=ucrpa_bands(2).AND.band1>=ucrpa_bands(1)) then
                 do iat=1, cryst%nattyp(itypatcor)
                   do ispinor1=1,nspinor
                     do m1=1,2*lcor+1
                       fac2=fac2-real(coeffW_BZ(iat,spin,band1,ik_bz,ispinor1,m1)*&
&                                conjg(coeffW_BZ(iat,spin,band1,ik_bz,ispinor1,m1)))
                     enddo
                   enddo
                 enddo
                 if(dtset%ucrpa==4) fac2=zero
               endif
               if        (band2<=ucrpa_bands(2).AND.band2>=ucrpa_bands(1)) then
                 do iat=1, cryst%nattyp(itypatcor)
                   do ispinor1=1,nspinor
                     do m1=1,2*lcor+1
                       fac3=fac3-real(coeffW_BZ(iat,spin,band2,ik_bz,ispinor1,m1)*&
&                                conjg(coeffW_BZ(iat,spin,band2,ik_bz,ispinor1,m1)))
                     enddo
                   enddo
                 enddo
                 if(dtset%ucrpa==4) fac3=zero
               endif
               fac=real(fac2*fac3)
             endif
!             if(dtset%prtvol>=10) write(6,'(6i4,e15.5,a)') ik_bz,ik_bz,band1,band2,m1,m2,fac," q==0"
!             if(abs(fac-one)>0.00001) write(6,'(a,6i4,e15.5,a)') "*****FAC*********",ik_bz,ik_bz,band1,band2,m1,m2,fac," q==0"
!             if(dtset%prtvol>=10) write(6,'(a,6i4,e15.5,a)') "*****FAC*********",ik_bz,ik_bz,band1,band2,m1,m2,fac4," q==0"
             green_w=green_w*fac
           endif
!          ---------------- Ucrpa (end)

           ! Adler-Wiser expression, to be consistent here we use the KS eigenvalues (?)
           call accumulate_chi0_q0(ik_bz,isym_k,itim_k,Ep%gwcomp,nspinor,Ep%npwepG0,Ep,&
&           Cryst,Ltg_q,Gsph_epsG0,chi0,rhotwx,rhotwg,green_w,green_enhigh_w,deltaf_b1b2,chi0_head,chi0_lwing,chi0_uwing)

         CASE (1, 2)
           ! Spectral method, to be consistent here we use the KS eigenvalues.
           call accumulate_sfchi0_q0(ik_bz,isym_k,itim_k,nspinor,Ep%symchi,Ep%npwepG0,Ep%npwe,Cryst,Ltg_q,&
&            Gsph_epsG0,deltaf_b1b2,my_wl,iomegal,wl,my_wr,iomegar,wr,rhotwx,rhotwg,Ep%nomegasf,sf_chi0,sf_head,sf_lwing,sf_uwing)

         CASE DEFAULT
           MSG_BUG("Wrong spmeth")
         END SELECT

         ! Accumulating the sum rule on chi0. Eq. (5.284) in G.D. Mahan Many-Particle Physics 3rd edition.
         factor=spin_fact*qp_occ(band2,ik_ibz,spin)

         call accumulate_chi0sumrule(ik_bz,Ep%symchi,Ep%npwe,factor,deltaeGW_b1b2,&
&          Ltg_q,Gsph_epsG0,Ep%npwepG0,rhotwg,chi0_sumrule)

         if (Ep%gwcomp==1) then
           ! Include also the completeness correction in the sum rule.
           factor=-spin_fact*qp_occ(band2,ik_ibz,spin)
           call accumulate_chi0sumrule(ik_bz,Ep%symchi,Ep%npwe,factor,deltaeGW_enhigh_b2,&
&            Ltg_q,Gsph_epsG0,Ep%npwepG0,rhotwg,chi0_sumrule)
           if (band1==Ep%nbnds) then
             chi0_sumrule(:)=chi0_sumrule(:) + wtk_ltg(ik_bz)*spin_fact*qp_occ(band2,ik_ibz,spin)*deltaeGW_enhigh_b2
           end if
         end if

       end do !band2
     end do !band1

     if (Psps%usepaw==0.and.Ep%inclvkb/=0.and.Ep%symchi==1) then
       call vkbr_free(vkbr(ik_ibz)) ! Not need anymore as we loop only over IBZ.
     end if
   end do !ik_bz
 end do !spin

 ABI_FREE(igffteps0)

 call vkbr_free(vkbr)
 ABI_DT_FREE(vkbr)

 ! === After big fat loop over transitions, now MPI ===
 ! * Master took care of the contribution in case of (metallic|spin) polarized systems.
 SELECT CASE (Ep%spmeth)
 CASE (0)
   ! Adler-Wiser expression
   ! Sum contributions from each proc.
   ! Looping on frequencies to avoid problems with the size of the MPI packet.
   do io=1,Ep%nomega
     call xmpi_sum(chi0(:,:,io),comm,ierr)
   end do

 CASE (1, 2)
   ! Spectral method.
   call hilbert_transform(Ep%npwe,Ep%nomega,Ep%nomegasf,my_wl,my_wr,kkweight,sf_chi0,chi0,Ep%spmeth)

   if (allocated(sf_chi0)) then
     ABI_FREE(sf_chi0)
   end if

   ! Sum contributions from each proc ===
   ! Looping on frequencies to avoid problems with the size of the MPI packet
   do io=1,Ep%nomega
     call xmpi_sum(chi0(:,:,io),comm,ierr)
   end do

   call hilbert_transform_headwings(Ep%npwe,Ep%nomega,Ep%nomegasf,&
&   my_wl,my_wr,kkweight,sf_lwing,sf_uwing,sf_head,chi0_lwing,&
&   chi0_uwing,chi0_head,Ep%spmeth)

 CASE DEFAULT
   MSG_BUG("Wrong spmeth")
 END SELECT

 ! Divide by the volume
!$OMP PARALLEL WORKSHARE
   chi0=chi0*weight/Cryst%ucvol
!$OMP END PARALLEL WORKSHARE

 ! Collect sum rule. pi comes from Im[1/(x-ieta)] = pi delta(x)
 call xmpi_sum(chi0_sumrule,comm,ierr)
 chi0_sumrule=chi0_sumrule * pi * weight / Cryst%ucvol

 ! Collect heads and wings.
 call xmpi_sum(chi0_head,comm,ierr)
 call xmpi_sum(chi0_lwing,comm,ierr)
 call xmpi_sum(chi0_uwing,comm,ierr)

 chi0_head  = chi0_head * weight/Cryst%ucvol
 do io=1,Ep%nomega ! Tensor in the basis of the reciprocal lattice vectors.
   chi0_head(:,:,io) = MATMUL(chi0_head(:,:,io),Cryst%gmet) * (two_pi**2)
 end do
 chi0_lwing = chi0_lwing * weight/Cryst%ucvol
 chi0_uwing = chi0_uwing * weight/Cryst%ucvol

 ! ===============================================
 ! ==== Symmetrize chi0 in case of AFM system ====
 ! ===============================================
 ! Reconstruct $chi0{\down,\down}$ from $chi0{\up,\up}$.
 ! Works only in the case of magnetic group Shubnikov type IV.
 if (Cryst%use_antiferro) then
   call symmetrize_afm_chi0(Cryst,Gsph_epsG0,Ltg_q,Ep%npwe,Ep%nomega,chi0,chi0_head,chi0_lwing,chi0_uwing)
 end if

 ! ===================================================
 ! ==== Construct heads and wings from the tensor ====
 ! ===================================================
 do io=1,Ep%nomega
   do ig=2,Ep%npwe
     wng = chi0_uwing(ig,io,:)
     chi0(1,ig,io) = vdotw(Ep%qlwl(:,1),wng,Cryst%gmet,"G")
     wng = chi0_lwing(ig,io,:)
     chi0(ig,1,io) = vdotw(Ep%qlwl(:,1),wng,Cryst%gmet,"G")
   end do
   chq = MATMUL(chi0_head(:,:,io), Ep%qlwl(:,1))
   chi0(1,1,io) = vdotw(Ep%qlwl(:,1),chq,Cryst%gmet,"G")  ! Use user-defined small q
 end do

 ! Impose Hermiticity (valid only for zero or purely imaginary frequencies)
 ! MG what about metals, where we have poles around zero?
 !  do io=1,Ep%nomega
 !    if (ABS(REAL(Ep%omega(io)))<0.00001) then
 !      do ig2=1,Ep%npwe
 !        do ig1=1,ig2-1
 !         chi0(ig2,ig1,io)=GWPC_CONJG(chi0(ig1,ig2,io))
 !        end do
 !      end do
 !    end if
 !  end do
 !
 ! =====================
 ! ==== Free memory ====
 ! =====================
 ABI_FREE(bbp_ks_distrb)
 ABI_FREE(rhotwg)
 ABI_FREE(tabr_k)
 ABI_FREE(ur1_kibz)
 ABI_FREE(ur2_kibz)
 ABI_FREE(usr1_k)
 ABI_FREE(ur2_k)
 ABI_FREE(gw_gbound)

 if (Dtset%pawcross==1) then
   ABI_FREE(gboundf)
 end if

 if (allocated(green_enhigh_w))  then
   ABI_FREE(green_enhigh_w)
 end if
 if (allocated(gw_gfft))  then
   ABI_FREE(gw_gfft)
 end if
 if (allocated(wfwfg))  then
   ABI_FREE(wfwfg)
 end if
 if (allocated(kkweight))  then
   ABI_FREE(kkweight)
 end if
 if (allocated(omegasf))  then
   ABI_FREE(omegasf)
 end if
 if (allocated(green_w))  then
   ABI_FREE(green_w)
 end if

 if (allocated(sf_head))  then
   ABI_FREE(sf_head)
 end if
 if (allocated(sf_lwing))  then
   ABI_FREE(sf_lwing)
 end if
 if (allocated(sf_uwing))  then
   ABI_FREE(sf_uwing)
 end if
 if (allocated(gspfft_igfft))  then
   ABI_FREE(gspfft_igfft)
 end if

 call gsph_free(Gsph_FFT)

 if (Psps%usepaw==1) then ! deallocation for PAW.
   call pawcprj_free(Cprj1_bz )
   ABI_DT_FREE(Cprj1_bz)
   call pawcprj_free(Cprj2_bz )
   ABI_DT_FREE(Cprj2_bz)
   call pawcprj_free(Cprj1_ibz )
   ABI_DT_FREE(Cprj1_ibz)
   call pawcprj_free(Cprj2_ibz )
   ABI_DT_FREE(Cprj2_ibz)
   call pawpwij_free(Pwij)
   ABI_DT_FREE(Pwij)
   if (allocated(Pwij_fft)) then
     call pawpwij_free(Pwij_fft)
     ABI_DT_FREE(Pwij_fft)
   end if
   call pawhur_free(Hur)
   ABI_DT_FREE(Hur)
   if (Dtset%pawcross==1) then
     ABI_FREE(ur_ae1)
     ABI_FREE(ur_ae_onsite1)
     ABI_FREE(ur_ps_onsite1)
     ABI_FREE(ur_ae2)
     ABI_FREE(ur_ae_onsite2)
     ABI_FREE(ur_ps_onsite2)
     ABI_FREE(tabrf_k)
     ABI_FREE(gboundf)
     ABI_FREE(igfftepsG0f)
   end if
 end if

 if(dtset%ucrpa>=1) then
   ABI_DEALLOCATE(coeffW_BZ)
 endif

 call cwtime(cpu_time,wall_time,gflops,"stop")
 write(std_out,'(2(a,f9.1))')" cpu_time = ",cpu_time,", wall_time = ",wall_time

 DBG_EXIT("COLL")

end subroutine cchi0q0
!!***
