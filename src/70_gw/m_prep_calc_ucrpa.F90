!!****m* ABINIT/m_prep_calc_ucrpa
!! NAME
!!  m_prep_calc_ucrpa
!!
!! FUNCTION
!! Prepare data for the calculation of U with the CRPA method: oscillators strenghs and k-points.
!!
!! COPYRIGHT
!! Copyright (C) 2006-2020 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

MODULE m_prep_calc_ucrpa

 use defs_basis
 use m_abicore
 use m_gwdefs!,        only : czero_gw, cone_gw, j_gw, sigparams_t
 use m_xmpi
 use m_defs_ptgroups
 use m_errors

 use defs_datatypes,  only : pseudopotential_type, ebands_t
 use m_time,          only : timab
 use m_hide_blas,     only : xdotc
 use m_geometry,      only : normv
 use m_crystal,       only : crystal_t
 use m_fft_mesh,      only : rotate_FFT_mesh
 use m_bz_mesh,       only : kmesh_t, get_BZ_item, findqg0, has_IBZ_item
 use m_gsphere,       only : gsphere_t, gsph_fft_tabs
 use m_io_tools,      only : flush_unit, open_file
 use m_vcoul_dt
 use m_pawpwij,       only : pawpwff_t, pawpwij_t, pawpwij_init, pawpwij_free, paw_rho_tw_g, paw_cross_rho_tw_g
 use m_paw_pwaves_lmn,only : paw_pwaves_lmn_t
 use m_pawang,        only : pawang_type
 use m_pawtab,        only : pawtab_type
 use m_pawfgrtab,     only : pawfgrtab_type
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy, paw_overlap
 use m_paw_nhat,      only : pawmknhat_psipsi
 use m_paw_sym,       only : paw_symcprj
 use m_wfd,           only : wfd_t
 use m_oscillators,   only : rho_tw_g
 use m_esymm,         only : esymm_t, esymm_failed
 use m_read_plowannier, only : read_plowannier
 use m_plowannier, only : plowannier_type,operwan_realspace_type

 implicit none

 private

 public :: prep_calc_ucrpa
!!***

contains

!!****f* ABINIT/prep_calc_ucrpa
!! NAME
!! prep_calc_ucrpa
!!
!! FUNCTION
!! Prepare data for the calculation of U with the CRPA method: oscillators strenghs and k-points.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (FB, GMR, VO, LR, RWG, MG, RShaltaf,TApplencourt,BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! sigmak_ibz=Index of the k-point in the IBZ.
!! minbnd, maxbnd= min and Max band index for GW correction (for this k-point)
!! Gsph_x<gsphere_t>= Info on the G-sphere used for Sigma_x
!!    %nsym=number of symmetry operations
!!    %rottb(ng,timrev,nsym)=index of (IS) G where I is the identity or the inversion
!!      operation and G is one of the ng vectors in reciprocal space
!!    %timrev=2 if time-reversal symmetry is used, 1 otherwise
!!    %gvec(3,ng)=integer coordinates of each plane wave in reciprocal space
!! ikcalc=index in the array Sigp%kptgw2bz of the k-point where GW corrections are calculated
!! Kmesh <kmesh_t>
!!    %nbz=Number of points in the BZ
!!    %nibz=Number of points in IBZ
!!    %kibz(3,nibz)=k-point coordinates, irreducible Brillouin zone
!!    %kbz(3,nbz)=k-point coordinates, full Brillouin zone
!!    %ktab(nbz)= table giving for each k-point in the BZ (kBZ), the corresponding
!!    %ktabi(nbz)= for each k-point in the BZ defines whether inversion has to be considered
!!    %ktabp(nbz)= phase factor associated to tnons
!! gwx_ngfft(18)=Information about 3D FFT for the oscillator strengths, see ~abinit/doc/variables/vargs.htm#ngfft
!! gwx_nfftot=number of points of the FFT grid for GW wavefunctions
!! Vcp <vcoul_t datatype> containing information on the cutoff technique
!!    %vc_sqrt(npwx,nqibz)= square-root of the coulombian potential for q-points in the IBZ
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
!!  occ(Sigp%nbnds,Kmesh%nibz,nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!!  Paw_pwff<pawpwff_t>=Form factor used to calculate the onsite mat. elements of a plane wave.
!!  allQP_sym(%nkibz,%nsppol)<esymm_t>=Datatype collecting data on the irreducible representaions of the
!!    little group of kcalc in the KS representation as well as the symmetry of the bdgw_k states.
!! prtvol=Flags governing verbosity level.
!!
!! OUTPUT
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
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      findqg0,flush_unit,get_bz_item,gsph_fft_tabs,paw_cross_rho_tw_g
!!      paw_rho_tw_g,paw_symcprj,pawcprj_alloc,pawcprj_copy,pawcprj_free
!!      pawmknhat_psipsi,pawpwij_free,pawpwij_init,read_plowannier,rho_tw_g
!!      rotate_fft_mesh,timab,wfd_change_ngfft,wfd_get_cprj,wfd_get_ur
!!      wfd_paw_get_aeur,wrtout
!!
!! SOURCE

subroutine prep_calc_ucrpa(sigmak_ibz,ikcalc,itypatcor,minbnd,maxbnd,Cryst,QP_BSt,Sigp,Gsph_x,Vcp,Kmesh,Qmesh,lpawu,&
& M1_q_m,Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,&
& Psps,Wfd,Wfdf,allQP_sym,gwx_ngfft,ngfftf,&
& prtvol,pawcross,plowan_compute,rhot1_q_m,wanbz,rhot1)

#ifndef HAVE_CRPA_OPTIM
#ifdef FC_INTEL
#if  __INTEL_COMPILER<=1700
!DEC$ NOOPTIMIZE
#endif
#endif
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: sigmak_ibz,ikcalc,itypatcor,prtvol,lpawu,minbnd,maxbnd,pawcross,plowan_compute
 type(crystal_t),intent(in) :: Cryst
 type(ebands_t),target,intent(in) :: QP_BSt
 type(kmesh_t),intent(in) :: Kmesh,Qmesh
 type(vcoul_t),intent(in) :: Vcp
 type(gsphere_t),intent(in) :: Gsph_x
! type(littlegroup_t),intent(in) :: Ltg_k
 type(Pseudopotential_type),intent(in) :: Psps
 type(sigparams_t),target,intent(in) :: Sigp
 type(pawang_type),intent(in) :: Pawang
 type(wfd_t),target,intent(inout) :: Wfd,Wfdf
!arrays
 complex(dpc), intent(out) :: rhot1_q_m(cryst%nattyp(itypatcor),Wfd%nspinor,Wfd%nspinor,2*lpawu+1,2*lpawu+1,sigp%npwx,Qmesh%nibz)
 complex(dpc), intent(out) :: M1_q_m(cryst%nattyp(itypatcor),Wfd%nspinor,Wfd%nspinor,2*lpawu+1,2*lpawu+1,sigp%npwx,Qmesh%nibz)
 integer,intent(in) :: gwx_ngfft(18),ngfftf(18)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(pawpwff_t),intent(in) :: Paw_pwff(Psps%ntypat*Psps%usepaw)
 type(esymm_t),target,intent(in) :: allQP_sym(Wfd%nkibz,Wfd%nsppol)
 type(pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom*Psps%usepaw)
 type(paw_pwaves_lmn_t),intent(in) :: Paw_onsite(Cryst%natom)
 type(plowannier_type),intent(in) :: wanbz
 type(operwan_realspace_type),target,intent(inout) :: rhot1(Sigp%npwx,Qmesh%nibz)

!Local variables ------------------------------
!scalars
 integer,parameter :: use_pawnhat=0,ider0=0,ndat1=1
 integer :: bandinf,bandsup
 integer :: gwcalctyp,izero,ib_sum,ib,ib1,ib2,ig,ig_rot,ii,iik,itim_q,i2
 integer :: ik_bz,ik_ibz,isym_q,iq_bz,iq_ibz,spin,isym,itypatcor_read,jb,iat
 integer :: jik,jk_bz,jk_ibz,lcor,m1,m3,nspinor,nsppol,ifft
 integer :: ibsp,dimcprj_gw
 integer :: spad
 integer :: comm
 integer :: ispinor1,ispinor3,isym_kgw,isym_ki,gwx_mgfft,use_padfft,use_padfftf,gwx_fftalga,gwx_fftalgb
 integer :: gwx_nfftot,nfftf,mgfftf,ierr
 integer :: nhat12_grdim
 integer :: iatom1,iatom2,il1,il2,im1,im2,ispinor2,pos1,pos2,wan_jb,wan_ib_sum,pwx
 real(dp) :: fact_sp,theta_mu_minus_esum,tol_empty,norm,weight
 complex(dpc) :: ctmp,scprod,ph_mkgwt,ph_mkt,eikr
 logical :: iscompatibleFFT,q_is_gamma
 character(len=500) :: msg
!arrays
 integer :: g0(3),spinor_padx(2,4)
 integer,pointer :: igfftxg0(:),igfftfxg0(:)
 integer,allocatable :: gwx_gfft(:,:),gwx_gbound(:,:),gboundf(:,:)
 integer,allocatable ::  ktabr(:,:),irottb(:,:),ktabrf(:,:)
 real(dp) :: ksum(3),kgw(3),kgw_m_ksum(3),qbz(3),q0(3),tsec(2)
 real(dp) :: spinrot_kbz(4),spinrot_kgw(4)
 real(dp),pointer :: qp_ene(:,:,:),qp_occ(:,:,:)
 real(dp),allocatable :: nhat12(:,:,:),grnhat12(:,:,:,:)
 complex(gwpc),allocatable :: vc_sqrt_qbz(:)
 complex(gwpc),allocatable :: rhotwg_ki(:,:)
 complex(gwpc),allocatable :: wfr_bdgw(:,:),wfr_sum(:)
 complex(gwpc),allocatable :: ur_ae_sum(:),ur_ae_onsite_sum(:),ur_ps_onsite_sum(:)
 complex(gwpc),allocatable :: ur_ae_bdgw(:,:),ur_ae_onsite_bdgw(:,:),ur_ps_onsite_bdgw(:,:)
 complex(gwpc),pointer :: cg_jb(:),cg_sum(:)
 complex(dpc) :: ovlp(2)
 complex(dpc),allocatable :: coeffW_BZ(:,:,:,:,:,:)
 complex(dpc),pointer :: ptr_rhot(:,:,:,:,:)
 logical :: can_symmetrize(Wfd%nsppol)
 logical,allocatable :: bks_mask(:,:,:)
 type(pawcprj_type),allocatable :: Cprj_kgw(:,:),Cprj_ksum(:,:)
 type(pawpwij_t),allocatable :: Pwij_qg(:),Pwij_fft(:)
 type(esymm_t),pointer :: QP_sym(:)
 !type(plowannier_type) :: wan
 logical     :: ecriture=.FALSE.
 logical     :: l_ucrpa,luwindow
 integer     :: g0_dump(3),iq_ibz_dump,dumint(2)

!************************************************************************

 l_ucrpa=.true.

 DBG_ENTER("COLL")

 !
 ! === Initial check ===
 ABI_CHECK(Sigp%npwx==Gsph_x%ng,'')

 call timab(430,1,tsec) ! csigme (SigX)

 gwcalctyp=Sigp%gwcalctyp
 !
 ! === Initialize MPI variables ===
 comm = Wfd%comm

 !
 ! === Initialize some values ===
 nspinor = Wfd%nspinor
 nsppol  = Wfd%nsppol
 spinor_padx(:,:)=RESHAPE((/0,0,Sigp%npwx,Sigp%npwx,0,Sigp%npwx,Sigp%npwx,0/),(/2,4/))

 qp_ene => QP_BSt%eig(:,:,:)
 qp_occ => QP_BSt%occ(:,:,:)

 ! Exctract the symmetries of the bands for this k-point
 QP_sym => allQP_sym(sigmak_ibz,1:nsppol)

 ib1=minbnd
 ib2=maxbnd

 ! === Read Wannier function coefficients for Ucrpa
 ! === for future computation of rhot_m_q directly in this routine.

 dumint=0
 luwindow=.true.
! write(6,*) "cc",allocated(coeffW_BZ)
 if (plowan_compute <10)then 
   call read_plowannier(Cryst,bandinf,bandsup,coeffW_BZ,itypatcor_read,Kmesh,lcor,luwindow,&
     & nspinor,nsppol,pawang,prtvol,dumint)
   if(lcor/=lpawu) then
     msg = "lcor and lpawu differ in prep_calc_ucrpa"
     MSG_ERROR(msg)
   endif
 endif

 ! === End of read Wannier function coefficients for Ucrpa


 !
 ! === Index of the GW point in the BZ array, its image in IBZ and time-reversal ===
 jk_bz=Sigp%kptgw2bz(ikcalc)
 !write(6,*) "ikcalc,jk_bz",ikcalc,jk_bz
 !write(6,*) "ikcalc",Kmesh%bz(:,ikcalc)
 !write(6,*) "jk_bz",Kmesh%bz(:,jk_bz)
! jk_bz=ikcalc
 call get_BZ_item(Kmesh,jk_bz,kgw,jk_ibz,isym_kgw,jik,ph_mkgwt)
! write(6,*) "jk_ibz",Kmesh%ibz(:,jk_ibz)
! write(6,*) "jk_bz,jk_ibz",jk_bz,jk_ibz,isym_kgw,itim
 !%call get_IBZ_item(Kmesh,jk_ibz,kibz,wtk)
 spinrot_kgw(:)=Cryst%spinrot(:,isym_kgw)
 !
 write(msg,'(2a,3f8.3,a,i4,a,2(i3,a))')ch10,&
&  ' Calculating Oscillator element at k= ',kgw, "k-point number",ikcalc,&
&  ' bands n = from ',ib1,' to ',ib2,ch10
 call wrtout(std_out,msg,'COLL')


 if (ANY(gwx_ngfft(1:3) /= Wfd%ngfft(1:3)) ) then
   call wfd%change_ngfft(Cryst,Psps,gwx_ngfft)
 end if
 gwx_mgfft   = MAXVAL(gwx_ngfft(1:3))
 gwx_fftalga = gwx_ngfft(7)/100
 gwx_fftalgb = MOD(gwx_ngfft(7),100)/10

 if (pawcross==1) then
   mgfftf = MAXVAL(ngfftf(1:3))
 end if

 can_symmetrize = .FALSE.
 if (Sigp%symsigma>0) then
   can_symmetrize = .TRUE.
   if (gwcalctyp >= 20) then
    do spin=1,Wfd%nsppol
      can_symmetrize(spin) = .not.esymm_failed(QP_sym(spin))
      if (.not.can_symmetrize(spin)) then
        write(msg,'(a,i0,4a)')&
&         " Symmetrization cannot be performed for spin: ",spin,ch10,&
&         " band classification encountered the following problem: ",ch10,TRIM(QP_sym(spin)%err_msg)
        MSG_WARNING(msg)
      end if
    end do
   end if
   ABI_CHECK(nspinor==1,'Symmetrization with nspinor=2 not implemented')
 end if

 ABI_ALLOCATE(rhotwg_ki,(Sigp%npwx*nspinor,minbnd:maxbnd))
 rhotwg_ki=czero_gw
 ABI_ALLOCATE(vc_sqrt_qbz,(Sigp%npwx))
 !
 ! === Normalization of theta_mu_minus_esum ===
 ! * If nsppol==2, qp_occ $\in [0,1]$
 SELECT CASE (nsppol)
 CASE (1)
   fact_sp=half; tol_empty=0.01   ! below this value the state is assumed empty
   if (Sigp%nspinor==2) then
    fact_sp=one; tol_empty=0.005  ! below this value the state is assumed empty
   end if
 CASE (2)
   fact_sp=one; tol_empty=0.005 ! to be consistent and obtain similar results if a metallic
 CASE DEFAULT                    ! spin unpolarized system is treated using nsppol==2
   MSG_BUG('Wrong nsppol')
 END SELECT

 ! Remove empty states from the list of states that will be distributed.
 ABI_ALLOCATE(bks_mask,(Wfd%mband,Kmesh%nbz,nsppol))
 bks_mask=.FALSE.
 do spin=1,nsppol
   do ik_bz=1,Kmesh%nbz
     ik_ibz = Kmesh%tab(ik_bz)
     do ib_sum=1,Sigp%nbnds
       bks_mask(ib_sum,ik_bz,spin) = (qp_occ(ib_sum,ik_ibz,spin)>=tol_empty)
     end do
   end do
 end do

! ABI_ALLOCATE(proc_distrb,(Wfd%mband,Kmesh%nbz,nsppol))
! call sigma_distribution(Wfd,Kmesh,Ltg_k,Qmesh,nsppol,can_symmetrize,kgw,Sigp%mg0,my_nbks,proc_distrb,bks_mask=bks_mask)
! call sigma_distribute_bks(Wfd,Kmesh,Ltg_k,Qmesh,nsppol,can_symmetrize,kgw,Sigp%mg0,my_nbks,proc_distrb,bks_mask=bks_mask)
! write(6,*)"lim", ib1,ib2
! do ib_sum= ib1,ib2
!   do ik_bz=1, Kmesh%nbz
!     write(6,*) ib_sum,ik_bz, proc_distrb(ib_sum,ik_bz,1),Wfd%my_rank
!   enddo
! enddo

 ABI_DEALLOCATE(bks_mask)

 write(msg,'(a,i8)')" Will sum all (b,k,s) occupied states in Sigma_x for k-point",ikcalc
 call wrtout(std_out,msg,'PERS')
 !
 ! The index of G-G0 in the FFT mesh the oscillators ===
 ! * Sigp%mG0 gives the MAX G0 component to account for umklapp.
 ! * Note the size MAX(Sigp%npwx,Sigp%npwc).
 ABI_ALLOCATE(igfftxg0,(Gsph_x%ng))
 !
 ! === Precalculate the FFT index of $ R^{-1}(r-\tau) $ ===
 ! * S=\transpose R^{-1} and k_BZ = S k_IBZ
 ! * irottb is the FFT index of $R^{-1} (r-\tau)$ used to symmetrize u_Sk.
 gwx_nfftot = PRODUCT(gwx_ngfft(1:3))
 ABI_ALLOCATE(irottb,(gwx_nfftot,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,gwx_ngfft,irottb,iscompatibleFFT)
 if (.not.iscompatibleFFT) then
   msg = "FFT mesh is not compatible with symmetries. Results might be affected by large errors!"
   MSG_WARNING(msg)
 end if

 ABI_ALLOCATE(ktabr,(gwx_nfftot,Kmesh%nbz))
 do ik_bz=1,Kmesh%nbz
   isym=Kmesh%tabo(ik_bz)
   do ifft=1,gwx_nfftot
     ktabr(ifft,ik_bz)=irottb(ifft,isym)
   end do
 end do
 ABI_DEALLOCATE(irottb)

 if (Psps%usepaw==1 .and. pawcross==1) then
   nfftf = PRODUCT(ngfftf(1:3))
   ABI_ALLOCATE(irottb,(nfftf,Cryst%nsym))
   call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,ngfftf,irottb,iscompatibleFFT)

   ABI_ALLOCATE(ktabrf,(nfftf,Kmesh%nbz))
   do ik_bz=1,Kmesh%nbz
     isym=Kmesh%tabo(ik_bz)
     do ifft=1,nfftf
       ktabrf(ifft,ik_bz)=irottb(ifft,isym)
     end do
   end do
   ABI_DEALLOCATE(irottb)
 end if
 !
 ! === Additional allocations for PAW ===
 if (Psps%usepaw==1) then
   ABI_DATATYPE_ALLOCATE(Cprj_ksum,(Cryst%natom,nspinor))
   call pawcprj_alloc(Cprj_ksum,0,Wfd%nlmn_atm)

   nhat12_grdim=0
   if (use_pawnhat==1) then ! Compensation charge for \phi_a^*\phi_b
     call wrtout(std_out,"Using nhat12","COLL")
     ABI_ALLOCATE(nhat12  ,(2,gwx_nfftot,nspinor**2))
     ABI_ALLOCATE(grnhat12,(2,gwx_nfftot,nspinor**2,3*nhat12_grdim))
   end if
 end if ! usepaw==1
 !

 if (Sigp%symsigma>0) then
   !call littlegroup_print(Ltg_k,std_out,prtvol,'COLL')
   !
   ! === Find number of complexes and number of bands in each complex ===
   ! The tolerance is a little bit arbitrary (0.001 eV)
   ! It could be reduced, in particular in case of nearly accidental degeneracies
!   if (ANY(degtab/=0)) then ! If two states do not belong to the same complex => matrix elements of v_xc differ
!     write(msg,'(a,3f8.3,a)')' Degenerate states at k-point = ( ',kgw(:),' ).'
!     call wrtout(std_out,msg,'COLL')
!     do spin=1,nsppol
!       do ib=ib1,ib2
!         do jb=ib+1,ib2
!           if (degtab(ib,jb,spin)==1) then
!             write(msg,'(a,i2,a,i4,a,i4)')' (spin ',spin,')',ib,' <====> ',jb
!             call wrtout(std_out,msg,'COLL')
!             if (ABS(Sr%vxcme(ib,jk_ibz,spin)-Sr%vxcme(jb,jk_ibz,spin))>ABS(tol6*Sr%vxcme(jb,jk_ibz,spin))) then
!               write(msg,'(7a)')&
!&                ' It seems that an accidental degeneracy is occurring at this k-point ',ch10,&
!&                ' In this case, using symsigma=1 might lead to spurious results as the algorithm ',ch10,&
!&                ' will treat these states as degenerate, and it won''t be able to remove the degeneracy. ',ch10,&
!&                ' In order to avoid this deficiency, run the calculation using symsigma=0'
!               MSG_WARNING(msg)
!             end if
!           end if
!         end do
!       end do
!     end do
!   end if
 end if !symsigma

 ABI_ALLOCATE(wfr_sum,(gwx_nfftot*nspinor))
 if (pawcross==1) then
   ABI_ALLOCATE(ur_ae_sum,(nfftf*nspinor))
   ABI_ALLOCATE(ur_ae_onsite_sum,(nfftf*nspinor))
   ABI_ALLOCATE(ur_ps_onsite_sum,(nfftf*nspinor))
 end if

!!*******************************************
!!  Save ik_bz and Norm of G vectors.
!!   FOR THE UCRPA calculation
!!*******************************************
 if (ikcalc==1) then
   ecriture=.TRUE.
!    open(unit=2011,file='ikbz_COORD',form='formatted',status='unknown')
!    do ik_bz=1,Kmesh%nbz
!          call get_BZ_item(Kmesh,ik_bz,&
!                               kbz_coord,ik_ibz,isym_kgw,iik,ph_mkt)
!          write(2011,*) ik_bz,kbz_coord(:)
!    end do
!    close(2011)

!   if (prtvol>10.and.jk_bz==1) then ! probably just to print one time.
!         !!q=0 Forcement donc divergence pour G=0
!         if (open_file("normeG", msg, unit=2022, form="formatted", status="unknown") /= 0) then
!           MSG_ERROR(msg)
!         end if
!         write(2022,*) 1,real(CMPLX(SQRT(Vcp%i_sz),0.0_gwp)),real((4*3.14159265)**(0.5)/CMPLX(SQRT(Vcp%i_sz),0.0_gwp))
!            write G=0 term for the potential computed elsewhere.
!
!         do ig=2,Sigp%npwx
!               write(2022,*) ig,real(Vcp%vc_sqrt(ig,1)),real((4*3.14159265)**(0.5)/Vcp%vc_sqrt(ig,1))
!         end do
!            write potential for ig, q=0
!         close(2022)
!   end if

!  write header for q point written later.
!   if (Wfd%my_rank==0) then
!     open(unit=2015,file='iqbz_COORD',form='formatted',status='unknown')
!     write (2015,*) "q pour le k", jk_bz, ikcalc
!     close(2015)
!   endif

 else
   ecriture=.FALSE.
 end if

!!*******************************************
!!   End of print if ik_bz for UCRPA calc
!!*******************************************

 !
 ! =======================================
 ! ==== Begin loop over k_i in the BZ ====
 ! =======================================

 do spin=1,nsppol

!   if (ALL(proc_distrb(:,:,spin)/=Wfd%my_rank)) CYCLE
!   write(6,*) "AA",Wfd%my_rank,spin
   !
   ! * Load wavefunctions for GW corrections.
   ABI_STAT_ALLOCATE(wfr_bdgw,(gwx_nfftot*nspinor,ib1:ib2), ierr)
   ABI_CHECK(ierr==0, "out of memory in wfr_bdgw")
   do jb=ib1,ib2
     call wfd%get_ur(jb,jk_ibz,spin,wfr_bdgw(:,jb))
!     write(6,'(a,6i4)')"indforwfd" ,jb,jk_ibz,spin
   end do

   if (Wfd%usepaw==1) then ! * Load cprj for GW states, note the indexing.
     dimcprj_gw=nspinor*(ib2-ib1+1)
     ABI_DATATYPE_ALLOCATE(Cprj_kgw,(Cryst%natom,ib1:ib1+dimcprj_gw-1))
     call pawcprj_alloc(Cprj_kgw,0,Wfd%nlmn_atm)
     ibsp=ib1
     do jb=ib1,ib2
!           write(6,*) "has_cprj",Wfd%Wave(jb,jk_ibz,spin)%has_cprj
           Wfd%Wave(jb,jk_ibz,spin)%has_cprj=1
       call wfd%get_cprj(jb,jk_ibz,spin,Cryst,Cprj_ksum,sorted=.FALSE.)
       call paw_symcprj(jk_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,Cprj_ksum)
       call pawcprj_copy(Cprj_ksum,Cprj_kgw(:,ibsp:ibsp+(nspinor-1)))
       ibsp=ibsp+nspinor
     end do
     if (pawcross==1) then
       ABI_ALLOCATE(ur_ae_bdgw,(nfftf*nspinor,ib1:ib2))
       ABI_ALLOCATE(ur_ae_onsite_bdgw,(nfftf*nspinor,ib1:ib2))
       ABI_ALLOCATE(ur_ps_onsite_bdgw,(nfftf*nspinor,ib1:ib2))
       do jb=ib1,ib2
         call wfdf%paw_get_aeur(jb,jk_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,&
&          ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum)
         ur_ae_bdgw(:,jb)=ur_ae_sum
         ur_ae_onsite_bdgw(:,jb)=ur_ae_onsite_sum
         ur_ps_onsite_bdgw(:,jb)=ur_ps_onsite_sum
       end do
     end if
   end if

   do ik_bz=1,Kmesh%nbz
     !
     ! === Parallelization over k-points and spin ===
     ! * For the spin there is another check in the inner loop
!     if (ALL(proc_distrb(:,ik_bz,spin)/=Wfd%my_rank)) CYCLE
!     write(6,*) "BB",Wfd%my_rank,spin,ik_bz
     !
     ! * Find the corresponding irreducible k-point
     call get_BZ_item(Kmesh,ik_bz,ksum,ik_ibz,isym_ki,iik,ph_mkt)
     spinrot_kbz(:)=Cryst%spinrot(:,isym_ki)
!     write(6,'(a,6i4)')"indices" ,jk_bz,jk_ibz,ik_bz,ik_ibz,spin

     ! * Identify q and G0 where q+G0=k_GW-k_i
     kgw_m_ksum=kgw-ksum
!     write(6,*) "kgw       ",kgw
!     write(6,*) "ksum      ",ksum
     call findqg0(iq_bz,g0,kgw_m_ksum,Qmesh%nbz,Qmesh%bz,Sigp%mG0)

!       if(iq_bz/=1.or.ik_bz/=1) cycle
!     write(6,*) "g0",g0
!     write(6,*) " ik_bz=",ik_bz
!     write(6,*) " iq_bz=",iq_bz

     ! === Symmetrize the matrix elements ===
     ! * Sum only q"s in IBZ_k. In this case elements are weighted
     !   according to wtqp and wtqm. wtqm is for time-reversal.
!     wtqp=1; wtqm=0
!     if (can_symmetrize(spin)) then
!     !  if (Ltg_k%ibzq(iq_bz)/=1) CYCLE
!       wtqp=0; wtqm=0
!       do isym=1,Ltg_k%nsym_sg
!         wtqp=wtqp+Ltg_k%wtksym(1,isym,iq_bz)
!         wtqm=wtqm+Ltg_k%wtksym(2,isym,iq_bz)
!       end do
!     end if


     !
     ! * Find the corresponding irreducible q-point.
     call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)
     q_is_gamma = (normv(qbz,Cryst%gmet,"G") < GW_TOL_W0)

!!*******************************************
!!    Check if qbz belongs to IBZ.
!!    because dielectric matrix is computed in the IBZ.
!!    then output the q point.
!!*******************************************
!      write(6,*) "kkk1",ik_bz,jk_bz,iq_ibz

      if (.NOT.has_IBZ_item(Qmesh,qbz,iq_ibz_dump,g0_dump)) then
        cycle
      end if
      
      write(msg,'(2(a,i4),a,i3)')' prep_calc_ucrpa : ik_bz ',ik_bz,'/',Kmesh%nbz,' done'
      call wrtout(std_out,msg,'PERS')
!      write(6,*) "kkk1p",ik_bz,jk_bz,iq_ibz

!     write(std_out,*)'prep_calc_ucrpa:ik_bz ',ik_bz,'/',Kmesh%nbz,' done by processor ',Wfd%my_rank,"iq_BZ",iq_bz,"iQ_iBZ",iq_ibz
!     write(123,*)'prep_calc_ucrpa:jk_bz ',jk_bz,'ikmq_bz',ik_bz,'iq_ibz',iq_ibz

     !Ecriture du iq_ibz
!     if (ecriture.and.Wfd%my_rank==0) then
!            open(unit=2016,file='iqbz_COORD',form='formatted',status='unknown',position='append')
!            write(2016,*) iq_ibz,qbz(:),Qmesh%wt(iq_ibz)
!          !       write(2011,*) iq_bz,qbz(:)
!            close(2016)
!     end if
!!*******************************************
!!    End of modif for UCRPA.
!!*******************************************

     !
     ! Tables for the FFT of the oscillators.
     !  a) FFT index of the G-G0.
     ABI_ALLOCATE(gwx_gbound,(2*gwx_mgfft+8,2))
     call gsph_fft_tabs(Gsph_x,g0,gwx_mgfft,gwx_ngfft,use_padfft,gwx_gbound,igfftxg0)

     if ( ANY(gwx_fftalga == (/2,4/)) ) use_padfft=0 ! Pad-FFT is not coded in rho_tw_g
#ifdef FC_IBM
 ! XLF does not deserve this optimization (problem with [v67mbpt][t03])
 use_padfft = 0
#endif
     if (use_padfft==0) then
       ABI_DEALLOCATE(gwx_gbound)
       ABI_ALLOCATE(gwx_gbound,(2*gwx_mgfft+8,2*use_padfft))
     end if

     if (pawcross==1) then
       ABI_ALLOCATE(gboundf,(2*mgfftf+8,2))
       ABI_ALLOCATE(igfftfxg0,(Gsph_x%ng))
       call gsph_fft_tabs(Gsph_x,g0,mgfftf,ngfftf,use_padfftf,gboundf,igfftfxg0)
       if ( ANY(gwx_fftalga == (/2,4/)) ) use_padfftf=0
       if (use_padfftf==0) then
         ABI_DEALLOCATE(gboundf)
         ABI_ALLOCATE(gboundf,(2*mgfftf+8,2*use_padfftf))
       end if
     end if
     !
     ! === Evaluate oscillator matrix elements ===
     ! * $ <phj/r|e^{-i(q+G)}|phi/r> - <tphj/r|e^{-i(q+G)}|tphi/r> $ in packed form.
     if (Psps%usepaw==1.and.use_pawnhat==0) then
       q0 = qbz !;if (q_is_gamma) q0 = (/0.00001_dp,0.00001_dp,0.00001_dp/) ! GW_Q0_DEFAULT
       ABI_DATATYPE_ALLOCATE(Pwij_qg,(Psps%ntypat))
       call pawpwij_init(Pwij_qg,Sigp%npwx,q0,Gsph_x%gvec,Cryst%rprimd,Psps,Pawtab,Paw_pwff)
     end if
     !
     ! === Get Fourier components of the Coulomb interaction in the BZ ===
     ! * In 3D systems, neglecting umklapp,  vc(Sq,sG)=vc(q,G)=4pi/|q+G|
     ! * The same relation holds for 0-D systems, but not in 1-D or 2D systems. It depends on S.
     do ig=1,Sigp%npwx
       ig_rot = Gsph_x%rottb(ig,itim_q,isym_q)
       vc_sqrt_qbz(ig_rot)=Vcp%vc_sqrt(ig,iq_ibz)
     end do
     !
!     write(6,*) "kkk2",ik_bz,jk_bz,iq_ibz
     ! === Sum over bands ===
     !do ib_sum=1,Sigp%nbnds
      do ib_sum=ib1,ib2
       !write(6,*) "ib_sum",ib_sum
       !
       ! === Parallelism over spin ===
       ! * This processor has this k-point but what about spin?
!       if (proc_distrb(ib_sum,ik_bz,spin)/=Wfd%my_rank) CYCLE
!       write(6,*) "CC",Wfd%my_rank,spin,ik_bz,ib_sum
       !
       ! * Skip empty states.
       !if (qp_occ(ib_sum,ik_ibz,spin)<tol_empty) CYCLE

       call wfd%get_ur(ib_sum,ik_ibz,spin,wfr_sum)
!       write(6,'(a,3i4)')"indforwfd2" ,ib_sum,ik_ibz,spin
!       write(6,*) wfd_ihave_ug(Wfd,ib_sum,ik_ibz,spin,"Stored"),Wfd%Wave(ib_sum,ik_ibz,spin)%has_ug
!       write(6,*) wfd_ihave_ur(Wfd,ib_sum,ik_ibz,spin,"Stored"),Wfd%Wave(ib_sum,ik_ibz,spin)%has_ur

       if (Psps%usepaw==1) then ! Load cprj for point ksum, this spin or spinor and *THIS* band.
         ! TODO MG I could avoid doing this but I have to exchange spin and bands ???
         ! For sure there is a better way to do this!
         call wfd%get_cprj(ib_sum,ik_ibz,spin,Cryst,Cprj_ksum,sorted=.FALSE.)
         call paw_symcprj(ik_bz,nspinor,1,Cryst,Kmesh,Pawtab,Pawang,Cprj_ksum)
         if (pawcross==1) then
           call wfdf%paw_get_aeur(ib_sum,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,&
&              ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum)
         end if
       end if

       do jb=ib1,ib2 ! Get all <k-q,ib_sum,s|e^{-i(q+G).r}|s,jb,k>
!        if(ib_sum.ne.jb) cycle
!        if(ib_sum.ne.1.and.ib_sum.ne.10) cycle
!        write(6,*) "jb",jb

         if (Psps%usepaw==1.and.use_pawnhat==1) then
           i2=jb; if (nspinor==2) i2=(2*jb-1)
           spad=(nspinor-1)

           izero=0
           call pawmknhat_psipsi(Cprj_ksum,Cprj_kgw(:,i2:i2+spad),ider0,izero,Cryst%natom,&
&            Cryst%natom,gwx_nfftot,gwx_ngfft,nhat12_grdim,nspinor,Cryst%ntypat,Pawang,Pawfgrtab,&
&            grnhat12,nhat12,pawtab)

#if 1
           msg = "reinstate optional Argument in rho_tw_g but mind inca slave!"
           MSG_ERROR(msg)
#else
           call rho_tw_g(nspinor,Sigp%npwx,gwx_nfftot,ndat1,gwx_ngfft,1,use_padfft,igfftxg0,gwx_gbound,&
&            wfr_sum       ,iik,ktabr(:,ik_bz),ph_mkt  ,spinrot_kbz,&
&            wfr_bdgw(:,jb),jik,ktabr(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&            nspinor,rhotwg_ki(:,jb),nhat12=nhat12)
#endif

         else
           call rho_tw_g(nspinor,Sigp%npwx,gwx_nfftot,ndat1,gwx_ngfft,1,use_padfft,igfftxg0,gwx_gbound,&
&            wfr_sum       ,iik,ktabr(:,ik_bz),ph_mkt  ,spinrot_kbz,&
&            wfr_bdgw(:,jb),jik,ktabr(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&            nspinor,rhotwg_ki(:,jb))

           if (Psps%usepaw==1.and.use_pawnhat==0) then ! Add on-site contribution, projectors are already in BZ.
             i2=jb; if (nspinor==2) i2=(2*jb-1)
             spad=(nspinor-1)
             call paw_rho_tw_g(Sigp%npwx,nspinor,nspinor,Cryst%natom,Cryst%ntypat,Cryst%typat,Cryst%xred,Gsph_x%gvec,&
&              Cprj_ksum(:,:),Cprj_kgw(:,i2:i2+spad),Pwij_qg,rhotwg_ki(:,jb))
           end if
           if(iq_bz==1) then
             if((ib_sum/=jb).and.(abs(rhotwg_ki(1,jb))>tol8)) then
               if((ib_sum/=jb).and.(abs(rhotwg_ki(1,jb))>0.01_dp)) then
                 write(std_out,*) "Warning: precision is low, oscillator strengh should be zero and is :",rhotwg_ki(1,jb)
               !else
               !  write(std_out,*) "Warning1: oscillator strengh",rhotwg_ki(1,jb)
               endif
             endif
             if((ib_sum==jb).and.(abs(rhotwg_ki(1,jb)-1_dp)>tol8)) then
               if((ib_sum==jb).and.(abs(rhotwg_ki(1,jb)-1_dp)>0.01_dp))  then
                 write(std_out,*) "Warning: precision is low, oscillator strengh should be one and is :",rhotwg_ki(1,jb)
               !else
               !  write(std_out,*) "Warning1: oscillator strengh",rhotwg_ki(1,jb)
               endif
             endif
           endif
           if (Psps%usepaw==1.and.pawcross==1) then ! Add paw cross term
             call paw_cross_rho_tw_g(nspinor,Sigp%npwx,nfftf,ngfftf,1,use_padfftf,igfftfxg0,gboundf,&
&             ur_ae_sum,ur_ae_onsite_sum,ur_ps_onsite_sum,iik,ktabrf(:,ik_bz),ph_mkt,spinrot_kbz,&
&             ur_ae_bdgw(:,jb),ur_ae_onsite_bdgw(:,jb),ur_ps_onsite_bdgw(:,jb),jik,ktabrf(:,jk_bz),ph_mkgwt,spinrot_kgw,&
&             nspinor,rhotwg_ki(:,jb))
           end if
         end if

!  ************************************8
!    Compute M Matrix in Wannier basis
!  ************************************8

         if (ib_sum.GE.ib1.AND.ib_sum.LE.ib2) then
           call flush_unit(std_out)
           call flush_unit(ab_out)
           if (plowan_compute <10)then
             do iat=1, cryst%nattyp(itypatcor)
               do ispinor1=1,nspinor
                 do ispinor3=1,nspinor
                   do m1=1,2*lcor+1
                     do m3=1,2*lcor+1
                       M1_q_m(iat,ispinor1,ispinor3,m1,m3,:,iq_ibz)=M1_q_m(iat,ispinor1,ispinor3,m1,m3,:,iq_ibz)+&
&                      rhotwg_ki(:,jb)*coeffW_BZ(iat,spin,jb,jk_bz,ispinor3,m3)*conjg(coeffW_BZ(iat,spin,ib_sum,ik_bz,ispinor1,m1))
                     enddo
                   enddo
                 enddo
               enddo
             enddo
           end if
         endif
!  ************************************8
!  ************************************8
         !
         ! === Multiply by the square root of the Coulomb term ===
         ! * In 3-D systems, the factor sqrt(4pi) is included)
         do ii=1,nspinor
           spad=(ii-1)*Sigp%npwx
!!$omp parallel workshare
           rhotwg_ki(spad+1:spad+Sigp%npwx,jb)=rhotwg_ki(spad+1:spad+Sigp%npwx,jb)*vc_sqrt_qbz(1:Sigp%npwx)
!!$omp end parallel workshare
         end do
         !
         ! === Treat analytically the case q --> 0 ===
         ! * The oscillator is evaluated at q=O as it is considered constant in the small cube around Gamma
         !   while the Colulomb term is integrated out.
         ! * In the scalar case we have nonzero contribution only if ib_sum==jb
         ! * For nspinor==2 evalute <ib_sum,up|jb,up> and <ib_sum,dwn|jb,dwn>,
         !   impose orthonormalization since npwwfn might be < npwvec.
         if (ik_bz==jk_bz) then
           if (nspinor==1) then
             rhotwg_ki(1,jb)=czero_gw
             if (ib_sum==jb) rhotwg_ki(1,jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp)
           else
             ! TODO Recheck this!
             cg_sum  => Wfd%Wave(ib,ik_ibz,spin)%ug
             cg_jb   => Wfd%Wave(jb,jk_ibz,spin)%ug

             ctmp = xdotc(Wfd%npwarr(ik_ibz)*Wfd%nspinor,cg_sum,1,cg_jb,1)
             ovlp(1) = REAL(ctmp)
             ovlp(2) = AIMAG(ctmp)

             if (Psps%usepaw==1) then
               i2=(2*jb-1)
               ovlp = ovlp + paw_overlap(Cprj_ksum,Cprj_kgw(:,i2:i2+1),Cryst%typat,Pawtab)
             end if
             !ovlp(2) = -ovlp(1)
             !if (ib_sum==jb) ovlp(2)=cone_gw-ovlp(1)
             if (ib_sum==jb) then
               norm=DBLE(ovlp(1)+ovlp(2))
               ovlp(1)=DBLE(ovlp(1)/norm)
               ovlp(2)=DBLE(ovlp(2)/norm)
             else
               scprod=ovlp(1)+ovlp(2)
               ovlp(1)=ovlp(1)-scprod*half
               ovlp(2)=ovlp(2)-scprod*half
             end if
             rhotwg_ki(1          ,jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp)*ovlp(1)
             rhotwg_ki(Sigp%npwx+1,jb)=CMPLX(SQRT(Vcp%i_sz),0.0_gwp)*ovlp(2)
           end if
         end if

!  ************************************8
!   compute rhotwidle in Wannier basis for UcRPA (summed over k)
!  ************************************8
!         if (cryst%nsym==1) then
           weight=one
!         else
!           weight=Kmesh%wt(jk_bz)
!         endif
         if (ib_sum.GE.ib1.AND.ib_sum.LE.ib2) then
           call flush_unit(std_out)
           call flush_unit(ab_out)
           if (plowan_compute<10)then
             do iat=1, cryst%nattyp(itypatcor)
               do ispinor1=1,nspinor
                 do ispinor3=1,nspinor
                   do m1=1,2*lcor+1
                     do m3=1,2*lcor+1
                       if(m1==2.and.m3==2) then
                       endif
                       rhot1_q_m(iat,ispinor1,ispinor3,m1,m3,:,iq_ibz)=&
                         &rhot1_q_m(iat,ispinor1,ispinor3,m1,m3,:,iq_ibz)+&
                         &rhotwg_ki(:,jb)*coeffW_BZ(iat,spin,jb,jk_bz,ispinor3,m3)&
                         &*conjg(coeffW_BZ(iat,spin,ib_sum,ik_bz,ispinor1,m1))*weight
                     enddo
                   enddo
                 enddo
               enddo
             enddo
           else
             wan_jb=jb-wanbz%bandi_wan+1
             wan_ib_sum=ib_sum-wanbz%bandi_wan+1
             do pwx=1,sigp%npwx
               do iatom1=1,wanbz%natom_wan
               do iatom2=1,wanbz%natom_wan
                 !Loig Vaugier PhD eq. 5.11
                 eikr=exp(- cmplx(0.0,1.0) * two_pi * ( &
   kmesh%bz(1,ik_bz)* ( cryst%xred(1,wanbz%iatom_wan(iatom1)) - cryst%xred(1,wanbz%iatom_wan(iatom2)) )+&
   kmesh%bz(2,ik_bz)* ( cryst%xred(2,wanbz%iatom_wan(iatom1)) - cryst%xred(2,wanbz%iatom_wan(iatom2)) )+&
   kmesh%bz(3,ik_bz)* ( cryst%xred(3,wanbz%iatom_wan(iatom1)) - cryst%xred(3,wanbz%iatom_wan(iatom2)) )))
                 do pos1=1,size(wanbz%nposition(iatom1)%pos,1)
                 do pos2=1,size(wanbz%nposition(iatom2)%pos,1)
                   do il1=1,wanbz%nbl_atom_wan(iatom1)
                   do il2=1,wanbz%nbl_atom_wan(iatom2)
                     ptr_rhot=>rhot1(pwx,iq_ibz)%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl
                     do im1=1,2*wanbz%latom_wan(iatom1)%lcalc(il1)+1
                     do im2=1,2*wanbz%latom_wan(iatom2)%lcalc(il2)+1
                       do ispinor1=1,wanbz%nspinor
                       do ispinor2=1,wanbz%nspinor
      ptr_rhot(im1,im2,spin,ispinor1,ispinor2)=&
      &ptr_rhot(im1,im2,spin,ispinor1,ispinor2)+&
      &rhotwg_ki(pwx,jb)*wanbz%psichi(jk_bz,wan_jb,iatom1)%atom(il1)%matl(im1,spin,ispinor1)*&
      &conjg(wanbz%psichi(ik_bz,wan_ib_sum,iatom2)%atom(il2)%matl(im2,spin,ispinor2))*weight&
      *eikr
                                 enddo!im2
                               enddo!im1
                             enddo!il2
                           enddo!il1
                         enddo!pos2
                       enddo!pos1
                     enddo!iatom2
                   enddo!iatom1
                 enddo!ispinor2
               enddo!ispinor1
             enddo!pwx
           endif!plowan_compute<10
         end if
!  ************************************8
!  ************************************8

       end do !jb Got all matrix elements from minbnd up to maxbnd.

       theta_mu_minus_esum=fact_sp*qp_occ(ib_sum,ik_ibz,spin)

     end do !ib_sum
     !
     ! Deallocate k-dependent quantities.
     ABI_DEALLOCATE(gwx_gbound)
     if (pawcross==1) then
       ABI_DEALLOCATE(gboundf)
     end if

     if (Psps%usepaw==1.and.use_pawnhat==0) then
       call pawpwij_free(Pwij_qg)
       ABI_DATATYPE_DEALLOCATE(Pwij_qg)
     end if

   end do !ik_bz Got all diagonal (off-diagonal) matrix elements.

   ABI_DEALLOCATE(wfr_bdgw)
   if (Wfd%usepaw==1) then
     call pawcprj_free(Cprj_kgw )
     ABI_DATATYPE_DEALLOCATE(Cprj_kgw)
     if (pawcross==1) then
       ABI_DEALLOCATE(ur_ae_bdgw)
       ABI_DEALLOCATE(ur_ae_onsite_bdgw)
       ABI_DEALLOCATE(ur_ps_onsite_bdgw)
     end if
   end if
 end do !spin

 ABI_DEALLOCATE(igfftxg0)
 if (pawcross==1) then
   ABI_DEALLOCATE(igfftfxg0)
 end if
 !
 ! Gather contributions from all the CPUs.
 !
 ! ===========================
 ! ==== Deallocate memory ====
 ! ===========================
 if (Psps%usepaw==1) then
   if (allocated(gwx_gfft))  then
     ABI_DEALLOCATE(gwx_gfft)
   end if
   call pawcprj_free(Cprj_ksum)
   ABI_DATATYPE_DEALLOCATE(Cprj_ksum)
   if (allocated(Pwij_fft)) then
     call pawpwij_free(Pwij_fft)
     ABI_DATATYPE_DEALLOCATE(Pwij_fft)
   end if
   if (use_pawnhat==1) then
     ABI_DEALLOCATE(nhat12)
     ABI_DEALLOCATE(grnhat12)
   end if
   if (pawcross==1) then
     ABI_DEALLOCATE(ur_ae_sum)
     ABI_DEALLOCATE(ur_ae_onsite_sum)
     ABI_DEALLOCATE(ur_ps_onsite_sum)
     ABI_DEALLOCATE(ktabrf)
   end if
 end if

 ABI_DEALLOCATE(wfr_sum)
 ABI_DEALLOCATE(rhotwg_ki)
 ABI_DEALLOCATE(vc_sqrt_qbz)
 ABI_DEALLOCATE(ktabr)
! ABI_DEALLOCATE(proc_distrb)
 if (plowan_compute<10) then
   ABI_DEALLOCATE(coeffW_BZ)
 endif


 call timab(430,2,tsec) ! csigme (SigX)

 DBG_EXIT("COLL")

end subroutine prep_calc_ucrpa
!!***

END MODULE m_prep_calc_ucrpa
!!***
