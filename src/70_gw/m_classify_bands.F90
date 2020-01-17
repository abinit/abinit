!!****m* ABINIT/m_classify_bands
!! NAME
!!  m_classify_bands
!!
!! FUNCTION
!!  Finds the irreducible representation associated to
!!  a set of degenerate bands at a given k-point and spin.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MG)
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

module m_classify_bands

 use defs_basis
 use m_abicore
 use m_esymm
 use m_errors

 use defs_datatypes,   only : pseudopotential_type, ebands_t
 use m_numeric_tools,  only : get_trace
 use m_symtk,          only : mati3inv
 use m_hide_blas,      only : xdotc, xdotu, xcopy
 use m_fft_mesh,       only : rotate_FFT_mesh, calc_ceigr
 use m_crystal,        only : crystal_t
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type, pawtab_get_lsize
 use m_pawfgrtab,      only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_print, pawfgrtab_free
 use m_pawcprj,        only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy
 use m_paw_pwaves_lmn, only : paw_pwaves_lmn_t, paw_pwaves_lmn_init, paw_pwaves_lmn_free
 use m_paw_sphharm,    only : setsym_ylm
 use m_paw_nhat,       only : nhatgrid
 use m_wfd,            only : wfd_t

 implicit none

 private
!!***

 public :: classify_bands
!!***

contains
!!***

!!****f* m_classify_bands/classify_bands
!! NAME
!! classify_bands
!!
!! FUNCTION
!!  This routine finds the irreducible representation associated to
!!  a set of degenerate bands at a given k-point and spin.
!!  The irreducible representation is obtained by rotating the set
!!  of degenerate wavefunctions using the symmetry operations in the little group of k.
!!  Two states are treated as degenerate if their energy differs by less than EDIFF_TOL.
!!
!! INPUTS
!!  Wfd(wfd_t)= structure gathering information on wave functions
!!    %nibz=Number of points in the IBZ.
!!    %nsppol=number of independent spin polarizations
!!    %usepaw=1 if PAW
!!  ik_ibz=The index of the k-point in the IBZ.
!!  spin=The spin index.
!!  ngfft(18)=Info on the FFT mesh to be used for evaluting u(r) and the rotated u(R^{1}(r-t)).
!!    ngfft must be compatible with the symmetries of the crystal and can differ from Wfd%ngfft.
!!    wfd_change_ngfft is called if ANY(Wfd%ngfft(1:3) =/ ngfft).
!!  Cryst<crystal_t>=Type gathering info on the crystal structure.
!!    %nsym=Number of operations in space group
!!    %ntypat=Number of type of atoms (onlu for PAW)
!!    %symrec(3,3,nsym)=Symmetry operations in reciprocal space (reduced coordinates)
!!    %tnons(3,nsym)=Fractional translations
!!    %typat(natom)=Type of each atom
!!  BSt<ebands_t>=Datatype with electronic energies.
!!  Pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data.
!!  Pawang <type(pawang_type)>=paw angular mesh and related data
!!  Psps<pseudopotential_type>
!!    %indlmn(6,lmnmax,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn (for each atom type)
!! Dtfil<datafiles_type>=variables related to files
!!    %unpaw
!! tolsym=Tolerance for the symmetries (input variable)
!!  [EDIFF_TOL]= tolerance on the energy difference of two states (if not specified is set to 0.005 eV)
!!
!! OUTPUT
!!  BSym<Bands_Symmetries>=structure containing info on the little group of the k-point as well
!!    as the character of the representation associated to each set of degenerate states
!!  if BSym%isymmorphic the symmetry analysis cannot be performed, usually it means that
!!   k is at zone border and there are non-symmorphic translations (see Notes)
!!
!! NOTES
!! * Let M(R_t) the irreducible representation associated to the space group symmetry (R_t).
!! * By convention M(R_t) multiplies wave functions as a row vector:
!!
!!    $ R_t \psi_a(r) = \psi_a (R^{-1}(r-\tau)) = \sum_b M(R_t)_{ba} \psi_b $
!!
!!   Therefore, if R_t belongs to the little group of k (i.e. Sk=k+G0), one obtains:
!!
!!    $ M_ab(R_t) = e^{-i(k+G0).\tau} \int e^{iG0.r} u_{ak}(r)^* u_{bk}(R^{-1}(r-\tau)) \,dr $.
!!
!! * The irreducible representation of the small _point_ group of k, M_ab(R), suffices to
!!   classify the degenerate eigenstates provided that particular conditions are fulfilled
!!   (see limitations below). The matrix is indeed given by:
!!
!!    $ M_ab(R) = e^{+ik.\tau} M_ab(R_t) = e^{-iG0.\tau} \int e^{iG0.r} u_{ak}(r)^* u_{bk}(R^{-1}(r-\tau))\,dr $
!!
!!   The phase factor outside the integral should be zero since symmetry analysis at border zone in non-symmorphic
!!   space groups is not available. Anyway it is included in our expressions for the sake of consistency.
!!
!! * For PAW there is an additional onsite terms involving <phi_i|phi_j(R^{-1}(r-\tau)> and
!!   the pseudized version that can be  evaluated using the rotation matrix for
!!    real spherical harmonis, zarot(mp,m,l,R). $ Y_{lm}(Rr)= \sum_{m'} zarot(m',m,ll,R) Y_{lm'}(r) $
!!
!!    $ M^{onsite}_ab(R_t) = sum_{c ij} <\tpsi_a| p_i^c>  <p_j^{c'}|\tpsi_b\> \times
!!       [ <\phi_i^c|\phi_j^{c'}> - <\tphi_i^c|\tphi_j^{c'}> ]. $
!!
!!    $ [ <\phi_i^c|\phi_j^{c'}> - <\tphi_i^c|\tphi_j^{c'}> ] = s_{ij} D_{\mi\mj}^\lj(R^{-1}) $
!!
!!   where c' is the rotated atom i.e c' = R^{-1}( c-\tau) and D is the rotation matrix for
!!   real spherical harmonics.
!!
!!   Remember that zarot(m',m,l,R)=zarot(m,m',l,R^{-1})
!!   and $ Y^l_m(ISG) = sum_{m'} D_{m'm}(S) Y_{m'}^l(G) (-i)^l $
!!       $ D_{m'm}^l (R) = D_{m,m'}^l (R^{-1}) $
!!
!! * LIMITATIONS: The method does not work if k is at zone border and the little group of k
!!                contains a non-symmorphic fractional translation.
!!
!! PARENTS
!!      sigma,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine classify_bands(Wfd,use_paw_aeur,first_band,last_band,ik_ibz,spin,ngfftf,&
& Cryst,BSt,Pawtab,Pawrad,Pawang,Psps,tolsym,BSym,&
& EDIFF_TOL) ! optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,spin,first_band,last_band
 real(dp),intent(in) :: tolsym
 real(dp),intent(in),optional :: EDIFF_TOL
 logical,intent(in) :: use_paw_aeur
 type(crystal_t),intent(in) :: Cryst
 type(pawang_type),intent(in) :: Pawang
 type(pseudopotential_type),intent(in) :: Psps
 type(wfd_t),intent(inout) :: Wfd
 type(ebands_t),target,intent(in) :: BSt
 type(esymm_t),intent(out) :: BSym
!arrays
 integer,intent(in) :: ngfftf(18)
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(Pawrad_type),intent(inout) :: Pawrad(Cryst%ntypat*Wfd%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: nspinor1=1
 integer :: dim_degs,ib1,ib2,ib_stop,ib_start,iclass,idg,sym_idx
 integer :: ir,isym,isym_class,tr_isym,jb1,jb2
 integer :: nr1,nr2,nr3,nsym_class,nfft,cplex !,ifgd,nfgd,ifft_sph
 integer :: ii,jj,lmax
 integer :: optcut,optgr0,optgr1,optgr2,optrad
 real(dp) :: EDIFF_TOL_,arg,fft_fact
 complex(dpc) :: exp_mikg0t,exp_ikg0t,cmat_ab
 logical :: iscompatibleFFT,found,only_trace
 character(len=500) :: msg
!arrays
 integer :: g0(3),toinv(Cryst%nsym),trial(3,3)
 integer,pointer :: Rm1_rmt(:)
 integer,target,allocatable :: irottb(:,:)
 integer,allocatable :: tmp_sym(:,:,:),l_size_atm(:)
 real(dp) :: kpt(3),kpg0(3),omat(2)
 real(dp),pointer :: ene_k(:)
 real(dp),pointer :: zarot(:,:,:,:)
 complex(dpc),allocatable :: eig0r(:,:),tr_emig0r(:,:)
 complex(gwpc),allocatable :: ur1(:),ur2(:),ur2_rot(:)
 type(pawcprj_type),allocatable :: Cprj_b1(:,:),Cprj_b2(:,:),Cprj_b2rot(:,:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)

! *************************************************************************

 DBG_ENTER("COLL")
 !
 ! Consistency check on input.
 ABI_CHECK(Wfd%nspinor==1,'nspinor/=1 not coded')
 !
 ! By default all bands are included
 !first_band=1; last_band=Wfd%nband(ik_ibz,spin)
 ABI_CHECK(first_band==1,"first_band/=1 not coded")
 ABI_CHECK(last_band<=Wfd%nband(ik_ibz,spin),"last_band cannot be > nband_k")

 EDIFF_TOL_=0.005/Ha_eV; if (PRESENT(EDIFF_TOL)) EDIFF_TOL_=ABS(EDIFF_TOL)

 call wfd%change_ngfft(Cryst,Psps,ngfftf)
 !
 ! === Get index of the rotated FFT points ===
 ! * FFT mesh in real space _must_ be compatible with symmetries.
 nr1=Wfd%ngfft(1)
 nr2=Wfd%ngfft(2)
 nr3=Wfd%ngfft(3)
 nfft=Wfd%nfft ! No FFT parallelism

 ABI_MALLOC(irottb,(nfft,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,Wfd%ngfft,irottb,iscompatibleFFT)

 if (.not.iscompatibleFFT) then
   write(msg,'(3a)')&
&    ' For symmetry analysis, the real space FFT mesh must be compatible with the symmetries of the space group',ch10,&
&    ' classify_bands will return. Action: change the input variable ngfftf '
   MSG_WARNING(msg)
   Bsym%err_status=1
   Bsym%err_msg= msg
   RETURN
 end if
 !
 ! only_trace=if .TRUE. only the trace of a single matrix per class is calculated (standard procedure if
 ! only the symmetry of bands is required). If .FALSE. all the matrices for each irreducible representation
 ! are calculated and stored in BSym
 only_trace=.FALSE.
 !
 ! ==========================================
 ! ==== Analyse k-point symmetries first ====
 ! ==========================================
 ! * The analysis is done here so that we already know if there is a problem.
 kpt=Wfd%kibz(:,ik_ibz)
 !
 !----Initialize the Bsym structure for this k-point and spin----!
 ! * NOTE that all the degenerate states should be included! No check is done.

 ene_k => BSt%eig(first_band:,ik_ibz,spin) ! Select a slice of eigenvalues

 call esymm_init(Bsym,kpt,Cryst,only_trace,Wfd%nspinor,first_band,last_band,EDIFF_TOL_,ene_k,tolsym)
 !Bsym%degs_bounds = Bsym%degs_bounds + (first_band -1)

 if (Bsym%err_status/=0) then
   write(msg,'(a,i0,a)')" esymm_init returned err_status= ",Bsym%err_status,&
&    " Band classifications cannot be performed."
   MSG_WARNING(msg)
   RETURN
 end if

 do ii=1,Cryst%nsym
   call mati3inv(Cryst%symrel(:,:,ii),trial)
   trial=transpose(trial)
   found=.FALSE.
   do jj=1,Cryst%nsym
     if (ALL(trial==Cryst%symrel(:,:,jj))) then
       toinv(ii)=jj
       !toinv(jj)=ii
       found=.TRUE.; EXIT
     end if
   end do
   if (.not.found) then
     MSG_ERROR("inverse not found! ")
   end if
 end do

 nullify(zarot)

 if (Wfd%usepaw==1) then ! Allocate cprj_k and cprj_krot to store a set of bands for a single (K,SPIN).
   ABI_DT_MALLOC(Cprj_b1   ,(Cryst%natom,Wfd%nspinor))
   call pawcprj_alloc(Cprj_b1,   0,Wfd%nlmn_atm)
   ABI_DT_MALLOC(Cprj_b2   ,(Cryst%natom,Wfd%nspinor))
   call pawcprj_alloc(Cprj_b2,   0,Wfd%nlmn_atm)
   ABI_DT_MALLOC(Cprj_b2rot,(Cryst%natom,Wfd%nspinor))
   call pawcprj_alloc(Cprj_b2rot,0,Wfd%nlmn_atm)

   !zarot => Pawang%zarot
   lmax = Pawang%l_max-1
   ABI_MALLOC(zarot,(2*lmax+1,2*lmax+1,lmax+1,Cryst%nsym))
   zarot = Pawang%zarot

   ABI_MALLOC(tmp_sym,(3,3,Cryst%nsym))
   do isym=1,Cryst%nsym
     tmp_sym(:,:,isym) = Cryst%symrel(:,:,isym)
     !tmp_sym(:,:,isym) = Cryst%symrel(:,:,toinv(isym))
     !tmp_sym(:,:,isym) = transpose(Cryst%symrel(:,:,isym))
     !tmp_sym(:,:,isym) = Cryst%symrec(:,:,isym)
     !tmp_sym(:,:,isym) = TRANSPOSE(Cryst%symrec(:,:,isym))
   end do
   !% call setsym_ylm(Cryst%rprimd,lmax,Cryst%nsym,3,Cryst%gprimd,tmp_sym,zarot)
   !call setsym_ylm(Cryst%gprimd,lmax,Cryst%nsym,1,Cryst%rprimd,tmp_sym,zarot)
   ABI_FREE(tmp_sym)
   zarot = Pawang%zarot

   cplex=1
   call pawtab_get_lsize(Pawtab,l_size_atm,Cryst%natom,Cryst%typat)
   ABI_DT_MALLOC(Pawfgrtab,(Cryst%natom))
   call pawfgrtab_init(Pawfgrtab,cplex,l_size_atm,Wfd%nspden,Cryst%typat)
   ABI_FREE(l_size_atm)

   optcut=1                     ! use rpaw to construct local_pawfgrtab
   optgr0=0; optgr1=0; optgr2=0 ! dont need gY terms locally
   optrad=1                     ! do store r-R


   call nhatgrid(Cryst%atindx1,Cryst%gmet,Cryst%natom,Cryst%natom,Cryst%nattyp,Wfd%ngfft,Cryst%ntypat,&
&    optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,pawtab,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)

   !call pawfgrtab_print(Pawfgrtab,unit=std_out,Wfd%prtvol=10)

   ABI_DT_MALLOC(Paw_onsite,(Cryst%natom))

   if (use_paw_aeur) then
     MSG_WARNING("Using AE wavefunction for rotation in real space!")
     call paw_pwaves_lmn_init(Paw_onsite,Cryst%natom,Cryst%natom,Cryst%ntypat,&
&                             Cryst%rprimd,Cryst%xcart,Pawtab,Pawrad,Pawfgrtab)
   end if
 end if
 !
 ! ===============================================
 ! ==== Calculate the representation matrices ====
 ! ===============================================
 fft_fact=one/nfft
 ABI_MALLOC(ur1,(nfft))
 ABI_MALLOC(ur2,(nfft))
 ABI_MALLOC(ur2_rot,(nfft))
 !
 ! * Precalculate eig0r = e^{iG0.r} on the FFT mesh.
 ABI_MALLOC(eig0r,(nfft,Bsym%nsym_gk))

 do isym=1,Bsym%nsym_gk
   g0=Bsym%g0(:,isym)
   call calc_ceigr(g0,nfft,nspinor1,Wfd%ngfft,eig0r(:,isym))
 end do

 if (Bsym%can_use_tr) then
   ABI_MALLOC(tr_emig0r,(nfft,Bsym%nsym_trgk))
   do isym=1,Bsym%nsym_trgk
     g0=Bsym%tr_g0(:,isym)
     call calc_ceigr(-g0,nfft,nspinor1,Wfd%ngfft,tr_emig0r(:,isym))
   end do
 end if
 !
 do idg=1,Bsym%ndegs ! Loop over the set of degenerate states.
   ib_start=Bsym%degs_bounds(1,idg)
   ib_stop =Bsym%degs_bounds(2,idg)
   dim_degs=Bsym%degs_dim(idg)

   do ib1=ib_start,ib_stop ! First band index in the degenerate set.
     jb1=ib1-ib_start+1

     ! debugging: use AE wave on dense FFT mesh.
     if (Wfd%usepaw==1..and.use_paw_aeur) then
       call wfd%paw_get_aeur(ib1,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur1)
     else
       call wfd%get_ur(ib1,ik_ibz,spin,ur1)
       if (Wfd%usepaw==1) then
         call wfd%ug2cprj(ib1,ik_ibz,spin,1,0,Cryst%natom,Cryst,Cprj_b1,sorted=.FALSE.)
       end if
     end if

     do ib2=ib_start,ib_stop ! Second band index in the degenerate set.

       if (Bsym%only_trace.and.ib1/=ib2) CYCLE ! Only the diagonal is needed.

       if (ib2==ib1) then
         call xcopy(nfft,ur1,1,ur2,1)
         if (Wfd%usepaw==1) then
           call pawcprj_copy(Cprj_b1,Cprj_b2)
         end if
       else
         !
         ! debugging: use AE wave on dense FFT mesh.
         if (Wfd%usepaw==1.and.use_paw_aeur) then
           call wfd%paw_get_aeur(ib2,ik_ibz,spin,Cryst,Paw_onsite,Psps,Pawtab,Pawfgrtab,ur2)
         else
           call wfd%get_ur(ib2,ik_ibz,spin,ur2)
           if (Wfd%usepaw==1) then
             call wfd%ug2cprj(ib2,ik_ibz,spin,1,0,Cryst%natom,Cryst,Cprj_b2,sorted=.FALSE.)
           end if
         end if
       end if
       !
       ! ===================================================
       ! ==== Loop over the classes of the little group ====
       ! ===================================================
       sym_idx=0
       do iclass=1,Bsym%nclass
         nsym_class=Bsym%nelements(iclass)
         !
         do isym_class=1,nsym_class ! Loop over elements in each class.
           sym_idx=sym_idx+1
           if (Bsym%only_trace.and.isym_class/=1) CYCLE ! Do it once if only the character is required.

           isym=Bsym%sgk2symrec(sym_idx)
           Rm1_rmt => irottb(:,isym)
           !
           ! Classify states according to the irreps of the little group of k.
           kpg0= kpt + Bsym%g0(:,sym_idx)

           arg=-two_pi * DOT_PRODUCT(kpg0,Cryst%tnons(:,isym))

           if (ABS(arg) > tol6) then
             exp_mikg0t=DCMPLX(DCOS(arg),DSIN(arg))
           else
             exp_mikg0t=cone
           end if

           !if (Wfd%usepaw==1) then
           !end if
           !
           ! === Rotate the right wave function and apply the phase ===
           ! * Note that the k-point is the same within a lattice vector.
           do ir=1,nfft
             ur2_rot(ir)=ur2(Rm1_rmt(ir))*eig0r(ir,sym_idx)
           end do

           ! * The matrix element on the FFT mesh.
           cmat_ab = xdotc(nfft,ur1,1,ur2_rot,1)*fft_fact*exp_mikg0t

           if (Wfd%usepaw==1.and..not.use_paw_aeur) then ! Add the on-site contribution.
             call rotate_cprj(kpt,isym,Wfd%nspinor,1,Cryst%natom,Cryst%nsym,Cryst%typat,Cryst%indsym,Cprj_b2,Cprj_b2rot)

             omat = paw_phirotphj(Wfd%nspinor,Cryst%natom,Cryst%typat,&
&              zarot(:,:,:,isym),Pawtab,Psps,Cprj_b1,Cprj_b2rot)

             cmat_ab = cmat_ab + DCMPLX(omat(1),omat(2)) !* exp_mikg0t
           end if
           !
           jb2=ib2-ib_start+1
           Bsym%Calc_irreps(idg)%mat(jb1,jb2,sym_idx)=cmat_ab

         end do !isym_class
       end do !iclass
       !
       ! =========================================================
       ! ==== Loop over the symmetries such that -Sk = k + G0 ====
       ! =========================================================
       ! <-k,a| S |k b>  = e^{i(k+G0).t} \int e^{-ig0.r} u_a u_b(R^{1}(r-t))
       if (Bsym%can_use_tr) then

         do tr_isym=1,Bsym%nsym_trgk

           isym=Bsym%tr_sgk2symrec(tr_isym)
           Rm1_rmt => irottb(:,isym)

           kpg0= kpt + Bsym%tr_g0(:,tr_isym)
           arg= two_pi * DOT_PRODUCT(kpg0,Cryst%tnons(:,isym))

           if (ABS(arg) > tol6) then
             exp_ikg0t=DCMPLX(DCOS(arg),DSIN(arg))
           else
             exp_ikg0t=cone
           end if
           !
           ! === Rotate the right wave function and apply the phase ===
           ! * Note that the k-point is the same within a lattice vector.
           do ir=1,nfft
             ur2_rot(ir)=ur2(Rm1_rmt(ir)) * tr_emig0r(ir,tr_isym)
           end do
           !
           ! * The matrix element on the FFT mesh.
           cmat_ab = xdotu(nfft,ur1,1,ur2_rot,1)*fft_fact*exp_ikg0t

           if (Wfd%usepaw==1.and..not.use_paw_aeur) then ! Add the on-site contribution. ! TODO rechek this part.
               call rotate_cprj(kpt,isym,Wfd%nspinor,1,Cryst%natom,Cryst%nsym,Cryst%typat,Cryst%indsym,Cprj_b2,Cprj_b2rot)
               omat = paw_phirotphj(Wfd%nspinor,Cryst%natom,Cryst%typat,&
&                zarot(:,:,:,isym),Pawtab,Psps,Cprj_b1,Cprj_b2rot,conjg_left=.TRUE.)
             cmat_ab = cmat_ab + DCMPLX(omat(1),omat(2)) !* exp_ikg0t
           end if
           !
           jb2=ib2-ib_start+1
           Bsym%trCalc_irreps(idg)%mat(jb1,jb2,tr_isym)=cmat_ab
         end do ! tr_isym
       end if

     end do !ib2
   end do !ib1
   !
   ! === Calculate the trace for each class ===
   if (Bsym%only_trace) then ! TODO this is valid if only trace.
     MSG_ERROR("Have to reconstruct missing traces")
   else
     do isym=1,Bsym%nsym_gk
       Bsym%Calc_irreps(idg)%trace(isym) = get_trace( Bsym%Calc_irreps(idg)%mat(:,:,isym) )
     end do
     if (Bsym%can_use_tr) then
       do tr_isym=1,Bsym%nsym_trgk
         Bsym%trCalc_irreps(idg)%trace(tr_isym) = get_trace( Bsym%trCalc_irreps(idg)%mat(:,:,tr_isym) )
       end do
     end if
   end if

 end do ! idg

 call esymm_finalize(Bsym,Wfd%prtvol)

 call esymm_print(Bsym,unit=std_out,prtvol=Wfd%prtvol)
 call esymm_print(Bsym,unit=ab_out ,prtvol=Wfd%prtvol)
 !
 ! ===================
 ! === Free memory ===
 ! ===================
 ABI_FREE(irottb)
 ABI_FREE(ur1)
 ABI_FREE(ur2)
 ABI_FREE(ur2_rot)
 ABI_FREE(eig0r)
 if (Bsym%can_use_tr)  then
   ABI_FREE(tr_emig0r)
 end if

 if (Wfd%usepaw==1) then
   call pawcprj_free(Cprj_b1)
   ABI_DT_FREE(Cprj_b1)
   call pawcprj_free(Cprj_b2)
   ABI_DT_FREE(Cprj_b2)
   call pawcprj_free(Cprj_b2rot)
   ABI_DT_FREE(Cprj_b2rot)
   ABI_FREE(zarot)
   call pawfgrtab_free(Pawfgrtab)
   ABI_DT_FREE(Pawfgrtab)
   call paw_pwaves_lmn_free(Paw_onsite)
   ABI_DT_FREE(Paw_onsite)
 end if

 DBG_EXIT("COLL")

end subroutine classify_bands
!!***

!----------------------------------------------------------------------

!!****f* m_classify_bands/rotate_cprj
!! NAME
!! rotate_cprj
!!
!! FUNCTION
!!  Rotate cprj matrix elements by applying the symmetry operation of index isym
!!  that preserves the given k-point within a reciprocal lattice vector.
!!
!! INPUTS
!! isym=index of the symmetry in the symrec arrays that preserves the given k-point within a reciprocal lattice vector
!! ntypat=number of types of atom.
!! natom=number of atoms.
!! Cryst<crystal_t>=Datatype gathering info on the unit cell.
!!   typat(natom)=type of each atom.
!! nbnds=number of bands for this k-point ans spin
!! Cprj_in(natom,nbnds)<type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk>
!!  with all NL projectors at fixed k-point
!!
!! OUTPUT
!! Cprj_out(natom,nbnds) <type(pawcprj_type)>= projection of the smooth PAW wave function onto
!!  projectors centered on equivalent sites of the crystal (non restricted to be in the firs unit cell)
!!  The equivalent site is defined according to the symmetry operation isym. Thus Cprj_out contains
!!
!!  Cprj_out(at,b)=<p_j^{R^{-1}(L_{at}-\tau)} | \tpsi_b> if  R is the isym operation  with fractional translation \tau
!!  L_{at} is the position of the initial atom inside the first unit cell
!!  Note that atom a might be in a cell different from the initial one. No wrapping is done.
!!
!! PARENTS
!!      classify_bands
!!
!! CHILDREN
!!
!! SOURCE

subroutine rotate_cprj(kpoint,isym,nspinor,nbnds,natom,nsym,typat,indsym,Cprj_in,Cprj_out)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbnds,nspinor,natom,isym,nsym
!arrays
 integer,intent(in) :: typat(natom),indsym(4,nsym,natom)
 real(dp),intent(in) :: kpoint(3)
 type(pawcprj_type),intent(in) :: Cprj_in(natom,nspinor*nbnds)
 type(pawcprj_type),intent(out) :: Cprj_out(natom,nspinor*nbnds)

!Local variables-------------------------------
!scalars
 integer :: iat,iband,itypat,iat_sym
 real(dp) :: kdotr0
!arrays
 integer :: r0(3)
 real(dp) :: phase_kr0(2)

! *************************************************************************

 !call pawcprj_copy(cprj_in,cprj_out)
 !RETURN

 do iat=1,natom
   itypat=typat(iat)
   !
   ! The index of the symmetric atom.
   ! * R^{-1} (xred(:,iat)-tnons) = xred(:,iat_sym) + r0.
   ! * phase_kr0 takes into account the case in which rotated atom is in another unit cell.
   iat_sym=indsym(4,isym,iat); r0=indsym(1:3,isym,iat)

   kdotr0 = two_pi*DOT_PRODUCT(kpoint,r0)
   phase_kr0(1) = DCOS(kdotr0)
   phase_kr0(2) = DSIN(kdotr0)

   !phase_kr0 = (/one,zero/)

   do iband=1,nspinor*nbnds
     Cprj_out(iat,iband)%cp(1,:)=  Cprj_in(iat_sym,iband)%cp(1,:)*phase_kr0(1) &
&                                 -Cprj_in(iat_sym,iband)%cp(2,:)*phase_kr0(2)

     Cprj_out(iat,iband)%cp(2,:)=  Cprj_in(iat_sym,iband)%cp(1,:)*phase_kr0(2) &
&                                 +Cprj_in(iat_sym,iband)%cp(2,:)*phase_kr0(1)
   end do
 end do ! iat

end subroutine rotate_cprj
!!***

!----------------------------------------------------------------------

!!****f* m_classify_bands/paw_phirotphj
!! NAME
!! paw_phirotphj
!!
!! FUNCTION
!!  This routine calculates
!!  <\tPsi_1|\tprj_i> <\tprj_j|\tPsi_2> [ <\phi_i|\phi_j(R^{-1}r> - <\tphi_i|\tphi_j(R^{-1}r> ]
!!
!! [ <\phi_i|\phi_j(R^{-1}r> - <\tphi_i|\tphi_j(R^{-1}r> ] = s_ij D_{mi,mi}^{li}(R)
!!
!! INPUTS
!! nspinor=Number of spinorial components.
!! natom=number of atoms
!! typat(natom)=type of eahc atom
!! zarot_isym
!! Pawtab(ntypat)<Pawtab_type>=PAW tabulated starting data
!! Psps<pseudopotential_type>=Info on pseudopotentials.
!! Cprj_b1(natom,nspinor)<type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk>
!!  with all NL projectors at fixed k-point
!! Cprj_b2(natom,nspinor)<type(pawcprj_type)>= projected input wave functions <Proj_i|Cnk>
!!  with all NL projectors at fixed k-point
!! [conjg_left]=.TRUE if the complex conjugate of the left wavefunctions has to be taken. Defaults to .FALSE.
!!
!! OUTPUT
!!  omat(2)=The onsite matrix element.
!!
!! PARENTS
!!
!! SOURCE

function paw_phirotphj(nspinor,natom,typat,zarot_isym,Pawtab,Psps,Cprj_b1,Cprj_b2,conjg_left) result(omat)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nspinor,natom
 logical,optional,intent(in) :: conjg_left
 type(pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: zarot_isym(:,:,:)
 real(dp) :: omat(2)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat)
 type(pawcprj_type),intent(in) :: Cprj_b1(natom,nspinor),Cprj_b2(natom,nspinor)

!Local variables-------------------------------
!scalars
 integer :: iat,il,ilmn,ilpm,im,itypat,jl,jlmn,jlpm,jm,k0lmn,klmn,nlmn
 real(dp) :: dmimj,fij,im_p,re_p,sij
 logical :: do_conjg_left

! *************************************************************************

 do_conjg_left = .FALSE.; if (PRESENT(conjg_left)) do_conjg_left = conjg_left

 if (nspinor/=1) then
   MSG_ERROR("nspinor/=1 not yet coded")
 end if

 ! === Rotate PAW projections ===
 ! * zarot_isym is the rotation matrix of real spherical harmonics associated to symrec(:,:,isym).
 ! * zarot_isym multiply harmonics as row vectors, we need R^{-1} but we read R and invert m,mp in the equation below
 omat=zero

 do iat=1,natom
   itypat=typat(iat)
   nlmn=Pawtab(itypat)%lmn_size

   do jlmn=1,nlmn
     k0lmn=jlmn*(jlmn-1)/2
     jl=Psps%indlmn(1,jlmn,itypat)
     jm=Psps%indlmn(2,jlmn,itypat)
     jlpm=1+jl+jm

     do ilmn=1,jlmn
       il=Psps%indlmn(1,ilmn,itypat)
       im=Psps%indlmn(2,ilmn,itypat)
       if (il/=jl.or.im/=jm) CYCLE ! Selection rule on l and m.
       ilpm=1+il+im

       klmn=k0lmn+ilmn
       sij=Pawtab(itypat)%sij(klmn) !; if (ABS(sij)<tol14) CYCLE

       ! Here we get the matrix associated to R^{-1}.
       dmimj=zarot_isym(ilpm,jlpm,jl+1)

       if (do_conjg_left) then  ! take the complex conjugate of the left cprj.
         re_p=  Cprj_b1(iat,1)%cp(1,ilmn) * Cprj_b2(iat,1)%cp(1,jlmn) &
&              -Cprj_b1(iat,1)%cp(2,ilmn) * Cprj_b2(iat,1)%cp(2,jlmn) &
&              +Cprj_b1(iat,1)%cp(1,jlmn) * Cprj_b2(iat,1)%cp(1,ilmn) &
&              -Cprj_b1(iat,1)%cp(2,jlmn) * Cprj_b2(iat,1)%cp(2,ilmn)

         im_p=  Cprj_b1(iat,1)%cp(1,ilmn) * Cprj_b2(iat,1)%cp(2,jlmn) &
&              +Cprj_b1(iat,1)%cp(2,ilmn) * Cprj_b2(iat,1)%cp(1,jlmn) &
&              -Cprj_b1(iat,1)%cp(1,jlmn) * Cprj_b2(iat,1)%cp(2,ilmn) &
&              -Cprj_b1(iat,1)%cp(2,jlmn) * Cprj_b2(iat,1)%cp(1,ilmn)
       else
         re_p=  Cprj_b1(iat,1)%cp(1,ilmn) * Cprj_b2(iat,1)%cp(1,jlmn) &
&              +Cprj_b1(iat,1)%cp(2,ilmn) * Cprj_b2(iat,1)%cp(2,jlmn) &
&              +Cprj_b1(iat,1)%cp(1,jlmn) * Cprj_b2(iat,1)%cp(1,ilmn) &
&              +Cprj_b1(iat,1)%cp(2,jlmn) * Cprj_b2(iat,1)%cp(2,ilmn)

         im_p=  Cprj_b1(iat,1)%cp(1,ilmn) * Cprj_b2(iat,1)%cp(2,jlmn) &
&              -Cprj_b1(iat,1)%cp(2,ilmn) * Cprj_b2(iat,1)%cp(1,jlmn) &
&              +Cprj_b1(iat,1)%cp(1,jlmn) * Cprj_b2(iat,1)%cp(2,ilmn) &
&              -Cprj_b1(iat,1)%cp(2,jlmn) * Cprj_b2(iat,1)%cp(1,ilmn)
       end if
       !
       ! * Accumulate the atom-centered contributions.
       fij = Pawtab(itypat)%dltij(klmn)/two
       omat(1)= omat(1) + fij*sij*re_p*dmimj
       omat(2)= omat(2) + fij*sij*im_p*dmimj

     end do !ilmn
   end do !jlmn
 end do !iat

end function paw_phirotphj
!!***

end module m_classify_bands
!!***
