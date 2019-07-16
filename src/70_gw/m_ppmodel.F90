!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ppmodel
!! NAME
!! m_ppmodel
!!
!! FUNCTION
!!  Module containing the definition of the ppmodel_t used to deal with
!!  the plasmonpole technique. Methods to operate on the object are also provided.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (MG, GMR, VO, LR, RWG, RS)
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

MODULE m_ppmodel

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_abicore
 use m_array
 use m_linalg_interfaces

 use m_fstrings,       only : sjoin, itoa
 use m_hide_lapack,    only : xhegv
 use m_gwdefs,         only : GW_Q0_DEFAULT, czero_gw
 use m_crystal,        only : crystal_t
 use m_bz_mesh,        only : kmesh_t, get_bz_item
 use m_gsphere,        only : gsphere_t
 use m_vcoul,          only : vcoul_t, cmod_qpg
 use m_fft_mesh,       only : g2ifft
 use m_fft,            only : fourdp
 use m_mpinfo,         only : destroy_mpi_enreg, initmpi_seq

 implicit none

 private
!!***

 integer,public,parameter :: PPM_NONE            = 0
 integer,public,parameter :: PPM_GODBY_NEEDS     = 1
 integer,public,parameter :: PPM_HYBERTSEN_LOUIE = 2
 integer,public,parameter :: PPM_LINDEN_HORSH    = 3
 integer,public,parameter :: PPM_ENGEL_FARID     = 4

 ! Flags giving the status of the pointers defined in ppmodel_t
 integer,private,parameter :: PPM_ISPOINTER   = 1 ! The pointer is used to store the address in memory.
 integer,private,parameter :: PPM_ISALLOCATED = 2 ! The pointer is used as an allocable array.

 ! Flags giving the status of the plasmon-pole tables
 integer,private,parameter :: PPM_NOTAB         = 0
 integer,private,parameter :: PPM_TAB_ALLOCATED = 1
 integer,private,parameter :: PPM_TAB_STORED    = 2

!----------------------------------------------------------------------

!!****t* m_ppmodel/ppmodel_t
!! NAME
!! ppmodel_t
!!
!! FUNCTION
!!  For the GW part of ABINIT, this structure datatype gathers all
!!  the information on the Plasmonpole technique used in the calculations
!!
!! NOTES
!!  If you modify this datatype, please check there there is no creation/destruction/copy routine,
!!  declared in another part of ABINIT, that might need to take into account your modification.
!!  Procedures that should be modified to reflect the changes done in the datatype have the tag @ppmodel_t.
!!
!! SOURCE

 type,public :: ppmodel_t

   !integers
   integer :: dm2_botsq
   ! =npwc if ppmodel=1,2; =1 if ppmodel=3,4

   integer :: dm_eig
   ! =npwc if ppmodel=3;   =0 if ppmodel=1,2,4

   integer :: dm2_otq
   ! =npwc if ppmodel=1,2; =1 if ppmodel=3,4

   integer :: invalid_freq
   ! what to do when PPM frequencies are invalid

   integer :: model
   ! The type of Plasmonpole model

   integer :: mqmem
   ! =nqibz if in-core solutions, =0 for out-of-core for which the last dimension in the PPm arrays has size 1.

   integer :: nqibz
   ! Number of q-points in the IBZ

   integer :: npwc
   ! Number of G vectors in $\tilde \epsilon $

   integer :: userho
   ! 1 if the ppmodel requires rho(G).

   integer :: iq_bz=0
   ! The index of the q-point in the BZ that is referenced by the internal pointer.

   integer :: bigomegatwsq_qbz_stat = PPM_ISPOINTER
   ! Status of bigomegatwsq.

   integer :: omegatw_qbz_stat = PPM_ISPOINTER
   ! Status of omegatw_qbz.

   integer :: eigpot_qbz_stat = PPM_ISPOINTER
   ! Status of new_eigpot_qbz_stat.

   real(dp) :: drude_plsmf
   ! Drude plasma frequency

   real(dp) :: force_plsmf=zero
   ! Force plasma frequency to that set by ppmfreq (not used if set zero).

   ! arrays
   logical,allocatable :: keep_q(:)
   ! .TRUE. if the ppmodel tables for this q in the IBZ are kept in memory.

   integer,allocatable :: has_q(:)
   ! Flag defining the status of the tables for the different q. See the PPM_TAB flags.

   type(array2_gwpc_t),pointer :: bigomegatwsq_qbz   => null()
   ! (Points|Stores) the symmetrized plasmon pole parameters $\tilde\Omega^2_{G Gp}(q_bz)$.

   type(array2_gwpc_t),pointer :: omegatw_qbz  => null()
   ! (Points|Stores) the symmetrized plasmon pole parameters $\tilde\omega_{G Gp}(q_bz)$.

   type(array2_gwpc_t),pointer :: eigpot_qbz   => null()
   ! (Points|Stores) the eigvectors of the symmetrized inverse dielectric matrix

   type(array2_gwpc_t),allocatable :: bigomegatwsq(:)
   ! bigomegatwsq(nqibz)%value(npwc,dm2_botsq)
   ! Plasmon pole parameters $\tilde\Omega^2_{G Gp}(q)$.

   type(array2_gwpc_t),allocatable :: omegatw(:)
   ! omegatw(nqibz)%value(npwc,dm2_otq)
   ! Plasmon pole parameters $\tilde\omega_{G Gp}(q)$.

   type(array2_gwpc_t),allocatable :: eigpot(:)
   ! eigpot(nqibz)%value(dm_eig,dm_eig)
   ! Eigvectors of the symmetrized inverse dielectric matrix

 end type ppmodel_t

 public :: ppm_get_qbz              ! Symmetrize the PPm parameters in the BZ.
 public :: ppm_nullify              ! Nullify all pointers
 public :: ppm_init                 ! Initialize dimensions and pointers
 public :: ppm_free                 ! Destruction method.
 public :: setup_ppmodel            ! Main Driver
 public :: new_setup_ppmodel        ! Main Driver
 public :: getem1_from_PPm          ! Reconstruct e^{-1}(w) from PPm.
 public :: getem1_from_PPm_one_ggp  ! Reconstruct e^{-1}(w) from PPm for one G,G' pair
 public :: get_ppm_eigenvalues
 public :: calc_sig_ppm             ! Matrix elements of the self-energy with ppmodel.
 public :: ppm_times_ket            ! Matrix elements of the self-energy with ppmodel.
 public :: ppm_symmetrizer
 public :: cqratio
!!***

CONTAINS  !==============================================================================
!!***

!!****f* m_ppmodel/ppm_get_qbz
!! NAME
!!  ppm_get_qbz
!!
!! FUNCTION
!!  Calculates the plasmonpole matrix elements in the full BZ zone.
!!
!! INPUTS
!!  PPm<ppmodel_t>=data type containing information on the plasmonpole technique
!!  Gsph<gsphere_t>=data related to the G-sphere
!!    %grottb
!!    %phmSGt
!!  Qmesh<kmesh_t>=Info on the q-mesh
!!    %nbz=number if q-points in the BZ
!!    %tab(nbz)=index of the symmeric q-point in the IBZ, for each point in the BZ
!!    %tabo(nbz)=the operation that rotates q_ibz onto \pm q_bz (depending on tabi)
!!    %tabi(nbz)=-1 if time-reversal has to be considered, 1 otherwise
!!  iq_bz=Index of the q-point in the BZ where PPmodel parameters have to be symmetrized
!!
!! OUTPUT
!!  botsq
!!  otq
!!  eig (only if PPm%ppmodel==3)
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0.
!!  In this case,indeed, the equation is different since we have to consider G-G0.
!!  There is however a check in sigma
!!
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!    If q_bz=Sq_ibz+G0:
!!
!!    $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then :
!!    $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!!
!! * Notice that eig is used only if PPm%model==3
!!
!! PARENTS
!!      calc_sigc_me,m_ppmodel
!!
!! CHILDREN
!!
!! SOURCE

subroutine ppm_get_qbz(PPm,Gsph,Qmesh,iq_bz,botsq,otq,eig)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz
 type(ppmodel_t),target,intent(inout) :: PPm
 type(gsphere_t),target,intent(in) :: Gsph
 type(kmesh_t),intent(in) :: Qmesh
 complex(gwpc),intent(out) :: botsq(:,:),otq(:,:),eig(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,iq_ibz,itim_q,isym_q,iq_curr,isg1,isg2
 !character(len=500) :: msg
!arrays
 integer, ABI_CONTIGUOUS pointer :: grottb(:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: phsgt(:),bigomegatwsq(:,:),omegatw(:,:)
! *********************************************************************

 !@ppmodel_t
 ! Save the index of the q-point for checking purpose.
 PPm%iq_bz = iq_bz

 ! Here there is a problem with the small q, still cannot use BZ methods
 iq_ibz=Qmesh%tab(iq_bz)
 isym_q=Qmesh%tabo(iq_bz)
 itim_q=(3-Qmesh%tabi(iq_bz))/2

 !call get_bz_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q,isirred=q_isirred)

 iq_curr=iq_ibz; if (PPm%mqmem==0) iq_curr=1

 grottb => Gsph%rottb (1:PPm%npwc,itim_q,isym_q)
 phsgt  => Gsph%phmSGt(1:PPm%npwc,isym_q)

 bigomegatwsq => PPm%bigomegatwsq(iq_curr)%vals
 omegatw      => PPm%omegatw(iq_curr)%vals

 ! Symmetrize the PPM parameters
 SELECT CASE (PPm%model)
 CASE (PPM_GODBY_NEEDS, PPM_HYBERTSEN_LOUIE)
   ! Plasmon pole frequencies otq are invariant under symmetry
!$omp parallel do private(isg1, isg2)
   do jj=1,PPm%npwc
     isg2 = grottb(jj)
     do ii=1,PPm%npwc
       isg1 = grottb(ii)
       botsq(isg1,isg2) = bigomegatwsq(ii,jj)*phsgt(ii)*CONJG(phsgt(jj))
       otq  (isg1,isg2) = omegatw(ii,jj)
     end do
   end do

 CASE (PPM_LINDEN_HORSH)
   ! * For notations see pag 22 of Quasiparticle Calculations in solid (Aulbur et al)
   !  If q_bz=Sq_ibz+G0 then:
   !
   ! $\omega^2_{ii}(q_bz) = \omega^2_{ii}(q)$        (otq array)
   ! $\alpha_{ii}(q_bz)   = \alpha_{ii}(q)$          (botq array
   ! $\Phi_{SG-G0}(q_bz)  = \Phi_{G}(q) e^{-iSG.t}$  (eigenvectors of e^{-1}, eig array)
   !
   do ii=1,PPm%npwc ! DM bands index
     botsq(ii,1) = bigomegatwsq(ii,1)
     otq  (ii,1) = omegatw     (ii,1)
     do jj=1,PPm%npwc
       eig(grottb(jj),ii) = PPm%eigpot(iq_curr)%vals(jj,ii) * phsgt(jj)
     end do
   end do
   if (itim_q==2) eig=CONJG(eig) ! Time-reversal

 CASE (PPM_ENGEL_FARID)
   ! * For notations see pag 23 of Quasiparticle Calculations in solid (Aulbur et al)
   ! If q_bz=Sq_ibz+G0 then:
   !
   ! $\omega^2_{ii}(q_bz) = \omega^2_{ii}(q)$        (otq array)
   ! $y_{SG-G0}(q_bz)     = y_{G}(q) e^{-iSG.t}$     (y=Lx)
   !
   do ii=1,PPm%npwc ! DM bands index
     otq(ii,1) = omegatw(ii,1)
     do jj=1,PPm%npwc
       botsq(grottb(jj),ii) = bigomegatwsq(jj,ii)*phsgt(jj)
     end do
   end do

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong PPm%model:',itoa(PPm%model)))
 END SELECT
 !
 ! * Take into account time-reversal symmetry.
 if (itim_q==2) then
!$omp parallel workshare
   botsq=CONJG(botsq)
!$omp end parallel workshare
 end if

end subroutine ppm_get_qbz
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_nullify
!! NAME
!!  ppm_nullify
!!
!! FUNCTION
!!  Nullify dynamic entities in a ppmodel_t object.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_ppmodel,m_screen
!!
!! CHILDREN
!!
!! SOURCE

subroutine ppm_nullify(PPm)

!Arguments ------------------------------------
 type(ppmodel_t),intent(inout) :: PPm
! *********************************************************************

 !@ppmodel_t
 ! types
 nullify(Ppm%bigomegatwsq_qbz)
 nullify(Ppm%omegatw_qbz)
 nullify(Ppm%eigpot_qbz)

end subroutine ppm_nullify
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_free
!! NAME
!!  ppm_free
!!
!! FUNCTION
!!  Deallocate all associated pointers defined in a variable of type ppmodel_t.
!!
!! SIDE EFFECTS
!!  PPm<ppmodel_t>=All dynamic memory is released.
!!
!! PARENTS
!!      m_screen,mrgscr,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine ppm_free(PPm)

!Arguments ------------------------------------
 type(ppmodel_t),intent(inout) :: PPm

!Local variables-------------------------------
!scalars
 integer :: dim_q,iq_ibz

! *********************************************************************

 !@ppmodel_t

 ! Be careful here to avoid dangling pointers.
 if (associated(PPm%bigomegatwsq_qbz)) then
   select case (PPm%bigomegatwsq_qbz_stat)
   case (PPM_ISALLOCATED)
     call array_free(PPm%bigomegatwsq_qbz)
   case (PPM_ISPOINTER)
     nullify(PPm%bigomegatwsq_qbz)
   end select
 end if

 if (associated(PPm%omegatw_qbz)) then
   select case (PPm%omegatw_qbz_stat)
   case (PPM_ISALLOCATED)
     call array_free(PPm%omegatw_qbz)
   case (PPM_ISPOINTER)
     nullify(PPm%omegatw_qbz)
   end select
 end if

 if (associated(PPm%eigpot_qbz)) then
   select case (PPm%eigpot_qbz_stat)
   case (PPM_ISALLOCATED)
     call array_free(PPm%eigpot_qbz)
   case (PPM_ISPOINTER)
     nullify(PPm%eigpot_qbz)
   end select
 end if

 dim_q=PPm%nqibz; if (PPm%mqmem==0) dim_q=1

#if 0
 do iq_ibz=1,dim_q
   call ppm_table_free(PPm,iq_ibz)
 end do
 ABI_FREE(PPm%bigomegatwsq)
 ABI_FREE(PPm%omegatw)
 ABI_FREE(PPm%eigpot)

#else
 if (allocated(PPm%bigomegatwsq)) then
   do iq_ibz=1,dim_q
     call array_free(PPm%bigomegatwsq(iq_ibz))
   end do
   ABI_DT_FREE(PPm%bigomegatwsq)
 end if
 !
 if (allocated(PPm%omegatw)) then
   do iq_ibz=1,dim_q
     call array_free(PPm%omegatw(iq_ibz))
   end do
   ABI_DT_FREE(PPm%omegatw)
 end if
 !
 if (allocated(PPm%eigpot)) then
   do iq_ibz=1,dim_q
     call array_free(PPm%eigpot(iq_ibz))
   end do
   ABI_DT_FREE(PPm%eigpot)
 end if
#endif

 ! logical flags must be deallocated here.
 ABI_SFREE(PPm%keep_q)
 ABI_SFREE(PPm%has_q)

end subroutine ppm_free
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_mallocq
!! NAME
!!  ppm_mallocq
!!
!! FUNCTION
!!  Allocate the ppmodel tables for the selected q-point in the IBZ.
!!
!! INPUT
!!  iq_ibz=Index of the q-point in the IBZ.
!!
!! SIDE EFFECTS
!!  PPm<ppmodel_t>=PPm tables are allocated.
!!
!! PARENTS
!!      m_ppmodel
!!
!! CHILDREN
!!
!! SOURCE

subroutine ppm_mallocq(PPm,iq_ibz)

!Arguments ------------------------------------
 integer,intent(in) :: iq_ibz
 type(ppmodel_t),intent(inout) :: PPm

! *********************************************************************

 !@ppmodel_t
 ABI_MALLOC(PPm%bigomegatwsq(iq_ibz)%vals, (PPm%npwc,PPm%dm2_botsq))

 ABI_MALLOC(PPm%omegatw(iq_ibz)%vals, (PPm%npwc,PPm%dm2_otq))

 ABI_MALLOC(PPm%eigpot(iq_ibz)%vals, (PPm%dm_eig,PPm%dm_eig))

 PPm%has_q(iq_ibz) = PPM_TAB_ALLOCATED

end subroutine ppm_mallocq
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_table_free
!! NAME
!!  ppm_table_free
!!
!! FUNCTION
!!  Free the ppmodel tables for the selected q-point in the IBZ.
!!
!! INPUT
!!  iq_ibz=Index of the q-point in the IBZ.
!!
!! SIDE EFFECTS
!!  PPm<ppmodel_t>=PPm tables are deallocated.
!!
!! PARENTS
!!      m_ppmodel
!!
!! CHILDREN
!!
!! SOURCE

subroutine ppm_table_free(PPm,iq_ibz)

!Arguments ------------------------------------
 integer,intent(in) :: iq_ibz
 type(ppmodel_t),intent(inout) :: PPm

! *********************************************************************

 !@ppmodel_t
 if (allocated(PPm%bigomegatwsq)) call array_free(PPm%bigomegatwsq(iq_ibz))
 if (allocated(PPm%omegatw)) call array_free(PPm%omegatw(iq_ibz))
 if (allocated(PPm%eigpot)) call array_free(PPm%eigpot(iq_ibz))

 PPm%has_q(iq_ibz) = PPM_NOTAB

end subroutine ppm_table_free
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_init
!! NAME
!!  ppm_init
!!
!! FUNCTION
!!  Initialize dimensions and other useful variables related to the PPmodel
!!
!! INPUTS
!! ppmodel=
!! drude_plsmf=
!!
!! SIDE EFFECTS
!!  PPm<ppmodel_t>=Arrays needed to store the ppmodel parameters are allocated with
!!   proper dimensions according to the plasmon-pole model.
!!
!! PARENTS
!!      m_screen,mrgscr,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine ppm_init(PPm,mqmem,nqibz,npwe,ppmodel,drude_plsmf,invalid_freq)

!Arguments ------------------------------------
 integer,intent(in) :: mqmem,nqibz,npwe,ppmodel,invalid_freq
 real(dp),intent(in) :: drude_plsmf
 type(ppmodel_t),intent(out) :: PPm

!Local variables-------------------------------
!scalars
 integer :: dim_q,iq_ibz,ierr
 logical :: ltest
 !character(len=500) :: msg
! *********************************************************************

 DBG_ENTER("COLL")

 !@ppmodel_t
 call ppm_nullify(PPm)

 PPm%invalid_freq=invalid_freq
 PPm%nqibz = nqibz
 PPm%mqmem = mqmem
 ltest = (PPm%mqmem==0.or.PPm%mqmem==PPm%nqibz)
! write(std_out,'(a,I0)') ' PPm%mqmem = ',PPm%mqmem
! write(std_out,'(a,I0)') ' PPm%nqibz = ',PPm%nqibz
! ABI_CHECK(ltest,'Wrong value for mqmem')

 PPm%npwc        = npwe
 PPm%model       = ppmodel
 PPm%drude_plsmf = drude_plsmf
 PPm%userho      = 0
 if (ANY(ppmodel==(/PPM_HYBERTSEN_LOUIE, PPM_LINDEN_HORSH, PPM_ENGEL_FARID/))) PPm%userho=1

 ABI_MALLOC(PPm%keep_q,(nqibz))
 PPm%keep_q=.FALSE.; if (PPm%mqmem>0) PPm%keep_q=.TRUE.

 ABI_MALLOC(PPm%has_q,(nqibz))
 PPm%has_q=PPM_NOTAB
 !
 ! Full q-mesh is stored or out-of-memory solution.
 dim_q=PPm%nqibz; if (PPm%mqmem==0) dim_q=1

 ABI_DT_MALLOC(PPm%bigomegatwsq, (dim_q))
 ABI_DT_MALLOC(PPm%omegatw,(dim_q))
 ABI_DT_MALLOC(PPm%eigpot,(dim_q))

 SELECT CASE (PPm%model)

 CASE (PPM_NONE)
   MSG_WARNING("Called with ppmodel==0")
   PPm%dm2_botsq = 0
   PPm%dm2_otq   = 0
   PPm%dm_eig    = 0
   RETURN

 CASE (PPM_GODBY_NEEDS, PPM_HYBERTSEN_LOUIE)
   PPm%dm2_botsq = PPm%npwc
   PPm%dm2_otq   = PPm%npwc
   PPm%dm_eig    = 1 ! Should be set to 0, but g95 doesnt like zero-sized arrays

 CASE (PPM_LINDEN_HORSH)
   PPm%dm2_botsq = 1
   PPm%dm2_otq   = 1
   PPm%dm_eig    = PPm%npwc

 CASE (PPM_ENGEL_FARID)
   PPm%dm2_botsq = PPm%npwc
   PPm%dm2_otq   = 1
   PPm%dm_eig    = 1 ! Should be set to 0, but g95 doesnt like zero-sized arrays

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong PPm%model:',itoa(PPm%model)))
 END SELECT
 !
 ! Allocate tables depending on the value of keep_q.
 do iq_ibz=1,dim_q
   !%if (keep_q(iq_ibz)) then
#if 1
   ABI_MALLOC_OR_DIE(PPm%bigomegatwsq(iq_ibz)%vals, (PPm%npwc,PPm%dm2_botsq), ierr)
   ABI_MALLOC_OR_DIE(PPm%omegatw(iq_ibz)%vals, (PPm%npwc,PPm%dm2_otq), ierr)
   ABI_MALLOC_OR_DIE(PPm%eigpot(iq_ibz)%vals, (PPm%dm_eig,PPm%dm_eig), ierr)
   !%endif
#else
   call ppm_mallocq(PPm,iq_ibz)
#endif
   !%end if
 end do

 DBG_EXIT("COLL")

end subroutine ppm_init
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/setup_ppmodel
!! NAME
!! setup_ppmodel
!!
!! FUNCTION
!!  Initialize some values of several arrays of the PPm datastructure
!!  that are used in case of plasmonpole calculations
!!  This is a wrapper around different plasmonpole routines.
!!
!! INPUTS
!!  Cryst<crystal_t>=Info on the unit cell and crystal symmetries.
!!  Qmesh<kmesh_t>=the q-mesh used for the inverse dielectric matrix
!!    %nibz=number of irreducible q-points
!!    %ibz(3,%nibz)=the irred q-point
!!  npwe=number of G vectors for the correlation part
!!  nomega=number of frequencies in $\epsilon^{-1}$
!!  omega=frequencies in epsm1
!!  epsm1=the inverse dielctric matrix
!!  ngfftf(18)=contain all needed information about the 3D fine FFT mesh, see ~abinit/doc/variables/vargs.htm#ngfft
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  nfftf=the number of points in the FFT mesh (for this processor)
!!  rhor_tot(nfftf)=the total charge in real space
!!  PPm<ppmodel_t>:
!!    %ppmodel=the type of  plasmonpole model
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  PPm<ppmodel_t>:
!!  == if ppmodel 1 or 2 ==
!!   %omegatw and %bigomegatwsq
!!  == if ppmodel 3 ==
!!   %omegatw, %bigomegatwsq and %eigpot
!!  == if ppmodel 4 ==
!!   %omegatw and %bigomegatwsq
!!
!! NOTES
!! * FFT parallelism not implemented.
!! * TODO: rhor_tot should be replaced by rhog_tot
!!
!! PARENTS
!!      calc_sigc_me,mrgscr,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine setup_ppmodel(PPm,Cryst,Qmesh,npwe,nomega,omega,epsm1,nfftf,gvec,ngfftf,rhor_tot,&
& iqiA) !Optional

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,npwe,nomega
 integer,intent(in),optional :: iqiA
 type(kmesh_t),intent(in) :: Qmesh
 type(crystal_t),intent(in) :: Cryst
 type(ppmodel_t),intent(inout) :: PPm
!arrays
 integer,intent(in) :: gvec(3,npwe),ngfftf(18)
 real(dp),intent(in) :: rhor_tot(nfftf)
 complex(dpc),intent(in) :: omega(nomega)
 complex(gwpc),intent(in) :: epsm1(:,:,:,:)

!Local variables-------------------------------
!scalars
 integer :: nqiA,iq_ibz
 real(dp) :: n_at_G_zero
 logical :: single_q
 character(len=500) :: msg
!scalars
 real(dp) :: qpt(3)
! *************************************************************************

 DBG_ENTER("COLL")

 !@ppmodel_t
 !
 ! === if iqiA is present, then consider only one qpoint to save memory ===
 ! * This means the object has been already initialized
 nqiA=Qmesh%nibz; single_q=.FALSE.
 if (PRESENT(iqiA)) then
   nqiA=1; single_q=.TRUE.
 end if
 !
 ! Allocate plasmonpole parameters
 ! TODO ppmodel==1 by default, should be set to 0 if AC and CD
 SELECT CASE (PPm%model)

 CASE (PPM_NONE)
   MSG_COMMENT(' Skipping Plasmompole model calculation')

 CASE (PPM_GODBY_NEEDS)
   ! * Note: the q-dependency enters only through epsilon^-1.
   do iq_ibz=1,nqiA
     call cppm1par(npwe,nomega,omega,PPm%drude_plsmf,&
&       epsm1(:,:,:,iq_ibz),PPm%omegatw(iq_ibz)%vals,PPm%bigomegatwsq(iq_ibz)%vals)
   end do

 CASE (PPM_HYBERTSEN_LOUIE)
   do iq_ibz=1,nqiA
     qpt = Qmesh%ibz(:,iq_ibz); if (single_q) qpt=Qmesh%ibz(:,iqiA)

     call cppm2par(qpt,npwe,epsm1(:,:,1,iq_ibz),ngfftf,gvec,Cryst%gprimd,rhor_tot,nfftf,Cryst%gmet,&
&      PPm%bigomegatwsq(iq_ibz)%vals,PPm%omegatw(iq_ibz)%vals,PPm%invalid_freq)
   end do
   !
   ! Quick-and-dirty change of the plasma frequency. Never executed in standard runs.
   if (PPm%force_plsmf>tol6) then ! Integrate the real-space density
      n_at_G_zero = SUM(rhor_tot(:))/nfftf
      ! Change the prefactor
      write(msg,'(2(a,es16.8))') 'Forced ppmfreq:',PPm%force_plsmf*Ha_eV,' nelec/ucvol:',n_at_G_zero
      MSG_WARNING(msg)
      PPm%force_plsmf = (PPm%force_plsmf**2)/(four_pi*n_at_G_zero)
      do iq_ibz=1,PPm%nqibz
        PPm%bigomegatwsq(iq_ibz)%vals = PPm%force_plsmf * PPm%bigomegatwsq(iq_ibz)%vals
        PPm%omegatw(iq_ibz)%vals      = PPm%force_plsmf * PPm%omegatw(iq_ibz)%vals
      end do
      write(msg,'(a,es16.8)') 'Plasma frequency forced in HL ppmodel, new prefactor is:',PPm%force_plsmf
      MSG_WARNING(msg)
   end if

 CASE (PPM_LINDEN_HORSH) ! TODO Check better double precision, this routine is in a messy state
   do iq_ibz=1,nqiA
     qpt = Qmesh%ibz(:,iq_ibz); if (single_q) qpt=Qmesh%ibz(:,iqiA)

     call cppm3par(qpt,npwe,epsm1(:,:,1,iq_ibz),ngfftf,gvec,Cryst%gprimd,rhor_tot,nfftf,&
&      PPm%bigomegatwsq(iq_ibz)%vals,PPm%omegatw(iq_ibz)%vals(:,1),PPm%eigpot(iq_ibz)%vals)
   end do

 CASE (PPM_ENGEL_FARID)  ! TODO Check better double precision, this routine is in a messy state
   do iq_ibz=1,nqiA
     qpt = Qmesh%ibz(:,iq_ibz); if (single_q) qpt=Qmesh%ibz(:,iqiA)
     if ((ALL(ABS(qpt)<1.0e-3))) qpt = GW_Q0_DEFAULT ! FIXME

     call cppm4par(qpt,npwe,epsm1(:,:,1,iq_ibz),ngfftf,gvec,Cryst%gprimd,rhor_tot,nfftf,&
&      PPm%bigomegatwsq(iq_ibz)%vals,PPm%omegatw(iq_ibz)%vals(:,1))
   end do

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong PPm%model:',itoa(PPm%model)))
 END SELECT

 DBG_EXIT("COLL")

end subroutine setup_ppmodel
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/getem1_from_ppm
!! NAME
!!  getem1_from_ppm
!!
!! FUNCTION
!!  Calculate the symmetrized inverse dielectric matrix from the parameters of the plasmon-pole model.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_ppmodel
!!
!! CHILDREN
!!
!! SOURCE

subroutine getem1_from_ppm(PPm,mpwc,iqibz,zcut,nomega,omega,Vcp,em1q,only_ig1,only_ig2)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpwc,iqibz,nomega
 type(ppmodel_t),intent(in) :: PPm
 type(vcoul_t),intent(in) :: Vcp
 real(dp),intent(in) :: zcut
 integer,optional,intent(in) :: only_ig1,only_ig2
!arrays
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(out) :: em1q(mpwc,mpwc,nomega)

!Local variables-------------------------------
!scalars
 integer :: ig1,ig2,io,idm,ig1_min,ig2_min,ig1_max,ig2_max
 real(dp) :: den
 complex(dpc) :: qpg1,qpg2,ug1,ug2
 complex(dpc) :: delta,em1ggp,otw,zzpq,yg1,yg2,bot1,bot2,chig1g2
 !character(len=500) :: msg

! *************************************************************************

 ABI_CHECK(PPm%mqmem/=0,'mqmem==0 not implemented')

 !TODO zcut should be an entry in PPm
 delta=CMPLX(zero,zcut)

 ! To save memory, a particular combination of
 ! ig1 and ig2 can be selected
 ig1_min = 1
 ig2_min = 1
 ig1_max = PPm%npwc
 ig2_max = PPm%npwc
 if (present(only_ig1)) then
   ig1_min = only_ig1
   ig1_max = only_ig1
 end if
 if (present(only_ig2)) then
   ig2_min = only_ig2
   ig2_max = only_ig2
 end if

 select case (PPm%model)
 case (PPM_GODBY_NEEDS, PPM_HYBERTSEN_LOUIE)
   do io=1,nomega
     !
     do ig2=ig2_min,ig2_max
       do ig1=ig1_min,ig1_max
        !den = omega(io)**2-REAL(PPm%omegatw(iqibz)%vals(ig1,ig2)**2)
        !if (den**2<zcut**2) den = omega(io)**2-REAL( (PPm%omegatw(iqibz)%vals(ig1,ig2)-delta)**2 )
        den = omega(io)**2 - REAL( (PPm%omegatw(iqibz)%vals(ig1,ig2)-delta)**2 )
        em1ggp = PPm%bigomegatwsq(iqibz)%vals(ig1,ig2)/den
        if (ig1==ig2) em1ggp=em1ggp+one
        em1q(ig1,ig2,io)=em1ggp
        !em1q(ig1,ig2,io)=em1ggp*Vcp%vc_sqrt(ig1,iqibz)*Vcp%vc_sqrt(ig2,iqibz)
       end do
     end do
     !
   end do !io

 case (PPM_LINDEN_HORSH)
   !TODO Check coefficients
   do io=1,nomega
     !
     do ig2=ig2_min,ig2_max
       do ig1=ig1_min,ig1_max
         !
         em1ggp=czero
         do idm=1,PPm%npwc
           !den=omega(io)**2-(PPm%omegatw(iqibz)%vals(ig1,ig2)-delta)**2
           !em1w(io)=em1w(io)+eigvec(ig1,idm,iqibz)*conjg(eigvec(ig2,idm,iqibz))*bigomegatwsq(ig1,ig2,iqibz)/den
           ug1 = PPm%eigpot(iqibz)%vals(ig1,idm)
           ug2 = PPm%eigpot(iqibz)%vals(ig2,idm)
           otw = PPm%bigomegatwsq(iqibz)%vals(idm,1)*PPm%omegatw(iqibz)%vals(idm,1)
           zzpq=PPm%bigomegatwsq(iqibz)%vals(idm,1)
           den=half*REAL(zzpq*otw*( one/(omega(io)-otw+delta) - one/(omega(io)+otw-delta) ))
           em1ggp=em1ggp+ug1*CONJG(ug2)*den
           !eigenvalues(idm,io)=one + half*REAL(zzpq*otw*( one/(omega(io)-otw+delta) - one/(omega(io)+otw-delta) ))
         end do
         if (ig2==ig1) em1ggp=em1ggp+one
         em1q(ig1,ig2,io)=em1ggp
       end do !ig1
     end do !ig2
     !
   end do !iomega

 case (PPM_ENGEL_FARID)
   ! Make e^-1
   do io=1,nomega
     !
     do ig2=ig2_min,ig2_max
       qpg2=one/Vcp%vc_sqrt(ig2,iqibz)
       do ig1=ig1_min,ig1_max
         qpg1=one/Vcp%vc_sqrt(ig1,iqibz)

         chig1g2=czero
         do idm=1,PPm%npwc
           otw =PPm%omegatw(iqibz)%vals(idm,1)
           bot1=PPm%bigomegatwsq(iqibz)%vals(ig1,idm)
           bot2=PPm%bigomegatwsq(iqibz)%vals(ig2,idm)
           yg1=SQRT(otw/four_pi)*qpg1*bot1
           yg2=SQRT(otw/four_pi)*qpg2*bot2
           chig1g2=chig1g2 + yg1*CONJG(yg2)/(omega(io)**2-(otw-delta)**2)
         end do

         em1ggp=four_pi*chig1g2/(qpg1*qpg2)
         if (ig1==ig2) em1ggp=em1ggp+one
         em1q(ig1,ig2,io)=em1ggp !*Vcp%vc_sqrt(ig1,iqibz)*Vcp%vc_sqrt(ig2,iqibz)
       end do !ig1
     end do !ig2
     !
   end do !iomega

 case default
   MSG_BUG(sjoin('Wrong PPm%model:',itoa(PPm%model)))
 end select

end subroutine getem1_from_ppm
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/getem1_from_ppm_one_ggp
!! NAME
!!  getem1_from_ppm
!!
!! FUNCTION
!!  Same as getem1_from_ppm, but does it for a single set of G,G' vectors
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      m_ppmodel,mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine getem1_from_ppm_one_ggp(PPm,iqibz,zcut,nomega,omega,Vcp,em1q,ig1,ig2)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,nomega
 type(ppmodel_t),intent(in) :: PPm
 type(vcoul_t),intent(in) :: Vcp
 real(dp),intent(in) :: zcut
 integer, intent(in) :: ig1,ig2
!arrays
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(out) :: em1q(nomega)

!Local variables-------------------------------
!scalars
 integer :: io,idm !,ig1_min,ig2_min,ig2_max
 real(dp) :: den
 complex(dpc) :: qpg1,qpg2,ug1,ug2
 complex(dpc) :: delta,em1ggp,otw,zzpq,yg1,yg2,bot1,bot2,chig1g2
 !character(len=500) :: msg

! *************************************************************************

 ABI_CHECK(PPm%mqmem/=0,'mqmem==0 not implemented')

 !TODO zcut should be an entry in PPm
 delta=CMPLX(zero,zcut)

 select case (PPm%model)

 case (PPM_GODBY_NEEDS, PPM_HYBERTSEN_LOUIE)
   do io=1,nomega
     !den = omega(io)**2-REAL(PPm%omegatw(iqibz)%vals(ig1,ig2)**2)
     !if (den**2<zcut**2) den = omega(io)**2-REAL( (PPm%omegatw(iqibz)%value(ig1,ig2)-delta)**2 )
     den = omega(io)**2-REAL( (PPm%omegatw(iqibz)%vals(ig1,ig2)-delta)**2 )
     em1ggp = PPm%bigomegatwsq(iqibz)%vals(ig1,ig2)/den
     if (ig1==ig2) em1ggp=em1ggp+one
     em1q(io)=em1ggp
     !em1q(io)=em1ggp*Vcp%vc_sqrt(ig1,iqibz)*Vcp%vc_sqrt(ig2,iqibz)
   end do !io

 case (PPM_LINDEN_HORSH)
   !TODO Check coefficients
   do io=1,nomega
     em1ggp=czero
     do idm=1,PPm%npwc
       !den=omega(io)**2-(PPm%omegatw(iqibz)%vals(ig1,ig2)-delta)**2
       !em1w(io)=em1w(io)+eigvec(ig1,idm,iqibz)*conjg(eigvec(ig2,idm,iqibz))*bigomegatwsq(ig1,ig2,iqibz)/den
       ug1 =PPm%eigpot(iqibz)%vals(ig1,idm)
       ug2 =PPm%eigpot(iqibz)%vals(ig2,idm)
       otw =PPm%bigomegatwsq(iqibz)%vals(idm,1)*PPm%omegatw(iqibz)%vals(idm,1)
       zzpq=PPm%bigomegatwsq(iqibz)%vals(idm,1)
       den=half*REAL(zzpq*otw*( one/(omega(io)-otw+delta) - one/(omega(io)+otw-delta) ))
       em1ggp=em1ggp+ug1*CONJG(ug2)*den
       !eigenvalues(idm,io)=one + half*REAL(zzpq*otw*( one/(omega(io)-otw+delta) - one/(omega(io)+otw-delta) ))
     end do

     if (ig2==ig1) em1ggp=em1ggp+one
     em1q(io)=em1ggp
   end do !iomega

 case (PPM_ENGEL_FARID)
   ! Make e^-1
   do io=1,nomega
     !
     qpg2=one/Vcp%vc_sqrt(ig2,iqibz)
     qpg1=one/Vcp%vc_sqrt(ig1,iqibz)

     chig1g2=czero
     do idm=1,PPm%npwc
       otw =PPm%omegatw(iqibz)%vals(idm,1)
       bot1=PPm%bigomegatwsq(iqibz)%vals(ig1,idm)
       bot2=PPm%bigomegatwsq(iqibz)%vals(ig2,idm)
       yg1=SQRT(otw/four_pi)*qpg1*bot1
       yg2=SQRT(otw/four_pi)*qpg2*bot2
       chig1g2=chig1g2 + yg1*CONJG(yg2)/(omega(io)**2-(otw-delta)**2)
     end do

     em1ggp=four_pi*chig1g2/(qpg1*qpg2)
     if (ig1==ig2) em1ggp=em1ggp+one
     em1q(io)=em1ggp !*Vcp%vc_sqrt(ig1,iqibz)*Vcp%vc_sqrt(ig2,iqibz)

   end do !iomega

 case default
   MSG_BUG(sjoin('Wrong PPm%model:',itoa(PPm%model)))
 end select

end subroutine getem1_from_ppm_one_ggp
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/get_ppm_eigenvalues
!! NAME
!!  get_ppm_eigenvalues
!!
!! FUNCTION
!!  Constructs the inverse dielectri matrix starting from the plasmon-pole
!!  parameters and calculates the frequency-dependent eigenvalues for each
!!  of the nomega frequencies specifies in the array omega.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_ppm_eigenvalues(PPm,iqibz,zcut,nomega,omega,Vcp,eigenvalues)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,nomega
 type(ppmodel_t),intent(in) :: PPm
 type(vcoul_t),intent(in) :: Vcp
 real(dp),intent(in) :: zcut
!arrays
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(out) :: eigenvalues(PPm%npwc,nomega)

!Local variables-------------------------------
!scalars
 integer :: info,lwork,negw,ig1,ig2,idx,sdim,iomega,ierr
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: ww(:),rwork(:)
 complex(dpc),allocatable :: work(:),Adpp(:),eigvec(:,:),wwc(:),vs(:,:),Afull(:,:)
 complex(dpc),allocatable :: em1q(:,:,:)
 logical,allocatable :: bwork(:)
 logical :: sortcplx !BUG in abilint

! *************************************************************************

 ABI_CHECK(PPm%mqmem/=0,'mqmem==0 not implemented')

 ABI_MALLOC(em1q,(PPm%npwc,PPm%npwc,nomega))

 call getem1_from_PPm(PPm,PPm%npwc,iqibz,zcut,nomega,omega,Vcp,em1q)

 do iomega=1,nomega
   !
   if (ABS(REAL(omega(iomega)))>0.00001) then
     !if (.TRUE.) then
     ! === Eigenvalues for a generic complex matrix ===

     lwork=4*2*PPm%npwc
     ABI_MALLOC(wwc,(PPm%npwc))
     ABI_MALLOC(work,(lwork))
     ABI_MALLOC(rwork,(PPm%npwc))
     ABI_MALLOC(bwork,(PPm%npwc))
     ABI_MALLOC(vs,(PPm%npwc,PPm%npwc))
     ABI_MALLOC(Afull,(PPm%npwc,PPm%npwc))
     Afull=em1q(:,:,iomega)

     !for the moment no sort, maybe here I should sort using the real part?
     call ZGEES('V','N',sortcplx,PPm%npwc,Afull,PPm%npwc,sdim,wwc,vs,PPm%npwc,work,lwork,rwork,bwork,info)
     if (info/=0) then
      write(msg,'(2a,i10)')' get_ppm_eigenvalues : Error in ZGEES, diagonalizing complex matrix, info = ',info
      call wrtout(std_out,msg,'COLL')
     end if

     eigenvalues(:,iomega)=wwc(:)

     ABI_FREE(wwc)
     ABI_FREE(work)
     ABI_FREE(rwork)
     ABI_FREE(bwork)
     ABI_FREE(vs)
     ABI_FREE(Afull)

   else
     ! === Hermitian Case ===
     lwork=2*PPm%npwc-1
     ABI_MALLOC(ww,(PPm%npwc))
     ABI_MALLOC(work,(lwork))
     ABI_MALLOC(rwork,(3*PPm%npwc-2))
     ABI_MALLOC(eigvec,(PPm%npwc,PPm%npwc))

     ABI_MALLOC_OR_DIE(Adpp,(PPm%npwc*(PPm%npwc+1)/2), ierr)
     write(std_out,*) 'in hermitian'

     idx=0
     do ig2=1,PPm%npwc
       do ig1=1,ig2
         idx=idx+1
         Adpp(idx)=em1q(ig1,ig2,iomega)
       end do
     end do

     ! For the moment we require also the eigenvectors.
     call ZHPEV('V','U',PPm%npwc,Adpp,ww,eigvec,PPm%npwc,work,rwork,info)

     if (info/=0) then
       write(msg,'(2a,i10)')' get_ppm_eigenvalues : Error diagonalizing matrix, info = ',info
       call wrtout(std_out,msg,'COLL')
     end if
     negw = (COUNT((REAL(ww)<tol6)))
     if (negw/=0) then
       write(msg,'(a,i0,a,i0,a,f8.4)')&
&       'Found negative eigenvalues. No. ',negw,' at iqibz= ',iqibz,' minval= ',MINVAL(REAL(ww))
        MSG_WARNING(msg)
     end if

     eigenvalues(:,iomega)=ww(:)

     ABI_FREE(ww)
     ABI_FREE(work)
     ABI_FREE(rwork)
     ABI_FREE(eigvec)
     ABI_FREE(Adpp)
   end if
   !
 end do !iomega

 ABI_FREE(em1q)

end subroutine get_ppm_eigenvalues
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/cppm1par
!! NAME
!! cppm1par
!!
!! FUNCTION
!! Calculate the plasmon-pole parameters big-omega-twiddle-squared and omega-twiddle from
!! epsilon-twiddle^-1 calculated for nomega (usually 2) frequencies omega=0 and omega=iE0.
!!
!! INPUTS
!!  epsm1(npwc,npwc,nomega)=dielectric matrix at nomega frequencies.
!!  npwc=number of plane waves
!!  nomega=number of frequencies (usually 2)
!!  omega(nomega)=frequencies
!!  omegaplasma=input variable or Drude plasma frequency
!!
!! OUTPUT
!!  bigomegatwsq(npwc,npwc)=parameter of the plasmon-pole model (see gwa.pdf file)
!!  omegatw(npwc,npwc)=parameter of the plasmon-pole model (see gwa.pdf file)
!!
!! TODO
!!  Calculation can be done in place.
!!
!! PARENTS
!!      m_ppmodel
!!
!! CHILDREN
!!
!! SOURCE

subroutine cppm1par(npwc,nomega,omega,omegaplasma,epsm1,omegatw,bigomegatwsq)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc
 real(dp),intent(in) :: omegaplasma
!arrays
 complex(gwpc),intent(in) :: epsm1(npwc,npwc,nomega)
 complex(dpc),intent(in) :: omega(nomega)
 complex(gwpc),intent(out) :: bigomegatwsq(npwc,npwc)
 complex(gwpc),intent(out) :: omegatw(npwc,npwc)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,io,io0,ioe0
 real(dp) :: e0,minomega
 character(len=500) :: msg
 complex(gwpc) :: AA,omegatwsq,diff,ratio
 complex(gwpc) :: epsm1_io0,epsm1_ioe0

! *************************************************************************

 DBG_ENTER("COLL")
 !
 ! === Find omega=0 and omega=imag (closest to omegaplasma) to fit the ppm parameters ===
 minomega=1.0d-3; io0=0
 do io=1,nomega
   if (ABS(omega(io))<minomega) then
     io0=io; minomega=ABS(omega(io))
   end if
 end do
 ABI_CHECK(io0/=0,"omega=0 not found")

 minomega=1.0d-3; e0=200.0; ioe0=0
 do io=1,nomega
   if (REAL(omega(io))<minomega.and.AIMAG(omega(io))>minomega) then
     if (ABS(AIMAG(omega(io))-omegaplasma)<ABS(e0-omegaplasma)) then
       ioe0=io; e0=AIMAG(omega(io))
     end if
   end if
 end do

 write(msg,'(a,f9.4,a)')' Imaginary frequency for fit located at: ',e0*Ha_eV,' [eV] '
 call wrtout(std_out,msg,'COLL')
 ABI_CHECK(ioe0/=0,"Imaginary omega not found")
 !
 ! ================================================================
 ! === Calculate plasmon-pole A parameter A=epsilon^-1(0)-delta ===
 ! ================================================================
 do ig=1,npwc
   do igp=1,npwc
     epsm1_io0  = epsm1(ig,igp,io0)
     epsm1_ioe0 = epsm1(ig,igp,ioe0)

     AA=epsm1_io0
     if (ig==igp) AA=AA-one

     ! === Calculate plasmon-pole omega-twiddle-square parameter ===
     ! XG201009 Strangely, the next formula does not work with gcc43-debug
     ! omegatwsq=(AA/(epsm1_io0-epsm1_ioe0)-one)*e0**2
     ! This seems to be due to precision issue at the level of division by a complex whose norm squared
     ! is below the smallest representable number.
     ! After many trials, I have decided to shift the difference by a small number ... well, not so small ...
     ! for numerical issues
     diff=epsm1_io0-epsm1_ioe0
     diff=diff+cmplx(tol10,tol10)
     ratio=AA/diff
     omegatwsq=(ratio-cone)*e0**2
     !
     ! If omega-twiddle-squared is negative,set omega-twiddle-squared to 1.0 (a reasonable way of treating
     ! such terms, in which epsilon**-1 was originally increasing along this part of the imaginary axis)
     ! (note: originally these terms were ignored in Sigma; this was changed on 6 March 1990.)

     if (REAL(omegatwsq)<=zero) omegatwsq=one
     !
     ! Get omega-twiddle
     ! * Neglect the imag part (if any) in omega-twiddle-squared
     omegatw(ig,igp)=SQRT(REAL(omegatwsq))
     !
     ! Get big-omega-twiddle-squared=-omega-twiddle-squared AA
     bigomegatwsq(ig,igp)=-AA*omegatw(ig,igp)**2

   end do !igp
 end do !ig

 write(msg,'(2a,f15.12,2a,2i5,a)')ch10,&
&  ' cppm1par : omega twiddle minval [eV]  = ',MINVAL(ABS(omegatw))*Ha_eV,ch10,&
&  '            omega twiddle min location = ',MINLOC(ABS(omegatw)),ch10
 call wrtout(std_out,msg,'COLL')

 DBG_EXIT("COLL")

end subroutine cppm1par
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/cppm2par
!! NAME
!! cppm2par
!!
!! FUNCTION
!!  Calculate plasmon-pole parameters of the Hybertsen and Louie model (PRB 34, 5390 (1986) [[cite:Hybertsen1986]])
!!
!! INPUTS
!!  qpt(3)=The coordinates of the q-point in the IBZ.
!!  epsm1(npwc,npwc)=symmetrized inverse dielectric (static limit is used)
!!  gmet(3,3)=metric in reciprocal space
!!  ngfftf(18)=contain all needed information about the 3D fine FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  npwc=number of plane waves in epsm1
!!  rhor(nfftf)=charge density on the real space FFT grid
!!  nfftf= total number of points in the fine FFT mesh  (for this processor)
!!  invalid_freq: what to do when PPM omega is found negative or imaginary
!!     0) drop it (default as specified in Hybersen-Louie original GW paper)
!!     1) set to 1 hartree
!!     2) set to infinity
!!
!! OUTPUT
!!  bigomegatwsq(npwc,npwc)= squared bare plasma frequencies
!!   \Omega^2_{G1 G2}(q) = 4\pi \frac {(q+G1).(q+G2)}/{|q+G1|^2} n(G1-G2)
!!  omegatw(npwc,npwc)= plasmon frequencies \tilde\omega_{G1 G2}(q) where:
!!  \tilde\omega^2_{G1 G2}(q) =
!!    \frac {\Omega^2_{G1 G2}(q)} {\delta_{G1 G2}-\tilde\epsilon^{-1}_{G1 G2} (q, \omega=0)}
!!
!! PARENTS
!!      m_ppmodel
!!
!! CHILDREN
!!
!! SOURCE

subroutine cppm2par(qpt,npwc,epsm1,ngfftf,gvec,gprimd,rhor,nfftf,gmet,bigomegatwsq,omegatw,invalid_freq)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwc,nfftf,invalid_freq
!arrays
 integer,intent(in) :: gvec(3,npwc)
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(in) :: qpt(3),gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: rhor(nfftf)
 complex(gwpc),intent(in) :: epsm1(npwc,npwc)
 complex(gwpc),intent(out) :: bigomegatwsq(npwc,npwc)
 complex(gwpc),intent(out) :: omegatw(npwc,npwc)

!Local variables-------------------------------
!scalars
 integer,parameter :: paral_kgb0=0
 integer :: ig,igp,nimwp,ngfft1,ngfft2,ngfft3,gmgp_idx,ierr
 real(dp) :: lambda,phi,AA
 logical,parameter :: use_symmetrized=.TRUE.,check_imppf=.FALSE.
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 real(dp) :: qlist(3,1)
 real(dp),allocatable :: tmp_rhor(:),qratio(:,:),qplusg(:),rhog_dp(:,:)
 complex(gwpc),allocatable :: omegatwsq(:,:)
 complex(gwpc),allocatable :: rhog(:),rhogg(:,:),temp(:,:)  !MG these should be double precision TODO
!*************************************************************************

 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfftf(2),ngfftf(3),'all')
 !
 ! === Calculate qratio(npwec,npvec) = (q+G).(q+Gp)/|q+G|^2 ===
 ABI_MALLOC_OR_DIE(qratio,(npwc,npwc), ierr)

 call cqratio(npwc,gvec,qpt,gmet,gprimd,qratio)
 !
 ! === Compute the density in G space rhor(R)--> rhog(G) ===
 ABI_MALLOC(rhog_dp,(2,nfftf))
 ABI_MALLOC(rhog,(nfftf))
 ngfft1=ngfftf(1)
 ngfft2=ngfftf(2)
 ngfft3=ngfftf(3)

 ABI_MALLOC(tmp_rhor,(nfftf))
 tmp_rhor=rhor ! To avoid having to use intent(inout).
 call fourdp(1,rhog_dp,tmp_rhor,-1,MPI_enreg_seq,nfftf,1,ngfftf,0)
 ABI_FREE(tmp_rhor)

 rhog(1:nfftf)=CMPLX(rhog_dp(1,1:nfftf),rhog_dp(2,1:nfftf))
 !
 ! Calculate the FFT index of each (G-Gp) vector and assign
 ! the value of the correspondent density simultaneously
 ABI_MALLOC_OR_DIE(rhogg,(npwc,npwc), ierr)

 ierr=0
 do ig=1,npwc
   do igp=1,npwc
     gmgp_idx = g2ifft(gvec(:,ig)-gvec(:,igp),ngfftf)
     if (gmgp_idx/=0) then
       rhogg(ig,igp)=rhog(gmgp_idx)
     else
       ierr=ierr+1
       rhogg(ig,igp)=czero
     end if
   end do
 end do

 if (ierr/=0) then
   write(msg,'(a,i4,3a)')&
&   'Found ',ierr,' G1-G2 vectors falling outside the FFT box. ',ch10,&
&   'Enlarge the FFT mesh to get rid of this problem. '
   MSG_WARNING(msg)
 end if

 rhogg=four_pi*rhogg
 ABI_FREE(rhog_dp)
 ABI_FREE(rhog)
 !
 ! Calculate GPP parameters
 ! unsymmetrized epsm1 -> epsm1=|q+Gp|/|q+G|*epsm1
 ABI_MALLOC(qplusg,(npwc))
 ABI_MALLOC(temp,(npwc,npwc))
 ABI_MALLOC_OR_DIE(omegatwsq,(npwc,npwc), ierr)

 temp = -epsm1(:,:)
 !
 ! RS still not obvious for me whether one shall use the symmetrized inverse DM or the unsymmetrized one
 ! the default here is to use the symmetrized one, I must discuss this with XG
 !
 ! MG it turns out that using the symmetrized inverse DM in the plasmon-pole
 ! equations give the same results for the squared plasmon frequencies omegatwsq while the
 ! squared bare plasma frequencies bigomegatwsq related to the symmetrized dielectric matrix
 ! are obtained multipling by |q+G1|/|q+G2|
 !
 if (.not.use_symmetrized) then
   qlist(:,1) = qpt
   call cmod_qpg(1,1,qlist,npwc,gvec,gprimd,qplusg) !MG TODO here take care of small q
   do ig=1,npwc
     do igp=1,npwc
       temp(ig,igp)=qplusg(igp)/qplusg(ig)*temp(ig,igp)
     end do
   end do
 end if

 nimwp=0
 do ig=1,npwc
   temp(ig,ig)=temp(ig,ig)+one
   do igp=1,npwc
     bigomegatwsq(ig,igp) = rhogg(ig,igp)*qratio(ig,igp)
     omegatwsq(ig,igp)=bigomegatwsq(ig,igp)/temp(ig,igp)
     !
     ! Set to an arbitrary value the omegawsq which become negative or imaginary
     ! in principle these correspond to cases where the imaginary part of epsm1 does not have
     ! a well defined peak. The imaginary part of epsm1 in these cases oscillates  with a small amplitude
     ! since the amplitude A_GGpr=-pi/2*bigomegatwsq/omegatw,
     ! it follows that bigomegatwsq shall be set to zero for these cases
     if ( REAL(omegatwsq(ig,igp))<= tol12 .or. AIMAG(omegatwsq(ig,igp))**2*tol12> REAL(omegatwsq(ig,igp))**2) then
       nimwp=nimwp+1

       if ( invalid_freq == 1 ) then
        ! set omegatwsq to 1 hartree
         omegatwsq(ig,igp)=cone
         AA = epsm1(ig,igp)
         if ( ig == igp ) AA = AA - one
         omegatw(ig,igp)=SQRT(REAL(omegatwsq(ig,igp)))
         bigomegatwsq(ig,igp)=-AA*omegatw(ig,igp)**2
       elseif ( invalid_freq == 2 ) then
         ! set omegatwsq to infinity
         omegatwsq(ig,igp)=cone/tol6
         AA = epsm1(ig,igp)
         if ( ig == igp ) AA = AA - one
         omegatw(ig,igp)=SQRT(REAL(omegatwsq(ig,igp)))
         bigomegatwsq(ig,igp)=-AA*omegatw(ig,igp)**2
       else
         ! simply ignore all cases of omegatw with imaginary values
         bigomegatwsq(ig,igp)=(0.,0.)
         omegatw(ig,igp)=(ten,0.)
       end if
       if (check_imppf) then
         write(msg,'(a,2i8)')' Imaginary plasmon frequency at : ',ig,igp
         call wrtout(std_out,msg,'COLL')
       end if
     else
       ! this part has been added to deal with systems without inversion symmetry
       ! this new implementation gives the same results as the previous one if
       ! omegatwsq is a pure real number and has the advantage of being an improved
       ! approach for systems without an inversion center.
       lambda=ABS(omegatwsq(ig,igp))
       phi=ATAN(AIMAG(omegatwsq(ig,igp))/REAL(omegatwsq(ig,igp)))
       omegatw(ig,igp)=SQRT(lambda/COS(phi))
       bigomegatwsq(ig,igp)=bigomegatwsq(ig,igp)*(1.-(0.,1.)*TAN(phi))
       ! Uncomment the following line and comment the previous to restore the old version.
       !omegatw(ig,igp)=sqrt(real(omegatwsq(ig,igp)))
     end if
   end do
 end do

 write(msg,'(a,3f12.6,a,i8,a,i8)')&
&  ' at q-point : ',qpt,&
&  ' number of imaginary plasmonpole frequencies = ',nimwp,' / ',npwc**2
 call wrtout(std_out,msg,'COLL')

 write(msg,'(2a,f12.8,2a,3i5)')ch10,&
&  ' cppm2par : omega twiddle minval [eV]  = ',MINVAL(ABS(omegatw))*Ha_eV,ch10,&
&  '            omega twiddle min location = ',MINLOC(ABS(omegatw))
 call wrtout(std_out,msg,'COLL')

 call destroy_mpi_enreg(MPI_enreg_seq)

 ABI_FREE(omegatwsq)
 ABI_FREE(rhogg)
 ABI_FREE(temp)
 ABI_FREE(qplusg)
 ABI_FREE(qratio)

end subroutine cppm2par
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/cppm3par
!! NAME
!! cppm3par
!!
!! FUNCTION
!! Calculate the plasmon-pole parameters using the von Linden-Horsh model (PRB 37, 8351, 1988) [[cite:vonderLinden1988]]
!! (see also Pag 22 of Quasiparticle Calculations in Solids [[cite:Aulbur2001]].
!!
!! INPUTS
!! epsm1(npwc,npwc))= symmetrized inverse dielectric
!! ngfftf(18)=contain all needed information about 3D fine FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! npwc=number of plane waves in epsm1
!! qratio=(q+G1).(q+G2)/(|q+G1|.|q+G2|)
!! rhor(nfftf)=charge density on the real space FFT grid
!! nfftf=number of points in the FFT grid (for this processor)
!! gvec(3,npwc)= G vectors in reduced coordinates
!!
!! OUTPUT
!!  omegatw(npwc,npwc)= plasmon pole positions
!!  bigomegatwsq(npwc,npwc)=(E_{q,ii}^{-1}-1)*omegatw
!!   where E^{-1} is the eigenvalue of the inverse dielectric matrix
!!  eigtot(npwc,npwc)=the eigvectors of the symmetrized inverse dielectric matrix
!!   (first index for G, second index for bands)
!!
!! PARENTS
!!      m_ppmodel
!!
!! CHILDREN
!!
!! SOURCE

subroutine cppm3par(qpt,npwc,epsm1,ngfftf,gvec,gprimd,rhor,nfftf,bigomegatwsq,omegatw,eigtot)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,npwc
!arrays
 integer,intent(in) :: gvec(3,npwc),ngfftf(18)
 real(dp),intent(in) :: qpt(3),gprimd(3,3),rhor(nfftf)
 complex(gwpc),intent(in) :: epsm1(npwc,npwc)
 complex(gwpc),intent(out) :: bigomegatwsq(npwc,1),eigtot(npwc,npwc)
 complex(gwpc),intent(out) :: omegatw(npwc)

!Local variables-------------------------------
!TODO these should be dp
!scalars
 integer,parameter :: paral_kgb0=0
 integer :: idx,ierr,ig,igp,ii,jj,ngfft1,ngfft2,ngfft3,gmgp_idx
 real(dp) :: num,qpg_dot_qpgp
 complex(dpc) :: conjg_eig
 logical :: qiszero
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg_seq
!arrays
 real(dp) :: b1(3),b2(3),b3(3),gppq(3),gpq(3),qlist(3,1)
 real(dp),allocatable :: eigval(:),qplusg(:),rhog_dp(:,:),zhpev2(:),tmp_rhor(:)
 complex(dpc),allocatable :: eigvec(:,:),matr(:),mm(:,:),rhog(:),rhogg(:,:)
 complex(dpc),allocatable :: zhpev1(:),zz(:)

!*************************************************************************

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfftf(2),ngfftf(3),'all')

 qiszero = (ALL(ABS(qpt)<1.0e-3))

 b1=two_pi*gprimd(:,1)
 b2=two_pi*gprimd(:,2)
 b3=two_pi*gprimd(:,3)

 ngfft1=ngfftf(1)
 ngfft2=ngfftf(2)
 ngfft3=ngfftf(3)

 ABI_MALLOC(rhog_dp,(2,nfftf))
 ABI_MALLOC(rhog,(nfftf))

 ABI_MALLOC_OR_DIE(rhogg,(npwc,npwc), ierr)
 !
 ! === Compute the density in G space rhog(r)--> rho(G) ===
 ! FIXME this has to be fixed, rho(G) should be passed instead of doing FFT for each q

 ABI_MALLOC(tmp_rhor,(nfftf))
 tmp_rhor=rhor ! To avoid having to use intent(inout).
 call fourdp(1,rhog_dp,tmp_rhor,-1,MPI_enreg_seq,nfftf,1,ngfftf,0)
 ABI_FREE(tmp_rhor)

 rhog(1:nfftf)=CMPLX(rhog_dp(1,1:nfftf),rhog_dp(2,1:nfftf))
 !
 ! Calculate the FFT index of each (G-Gp) vector and assign the value
 ! of the correspondent density simultaneously
 ierr=0
 do ig=1,npwc
   do igp=1,npwc
     gmgp_idx = g2ifft(gvec(:,ig)-gvec(:,igp),ngfftf)
     if (gmgp_idx/=0) then
       rhogg(ig,igp)=rhog(gmgp_idx)
     else
       ierr=ierr+1
       rhogg(ig,igp)=czero
     end if
   end do
 end do

 if (ierr/=0) then
   write(msg,'(a,i4,3a)')&
&   'Found ',ierr,' G1-G2 vectors falling outside the FFT box. ',ch10,&
&   'Enlarge the FFT mesh to get rid of this problem. '
   MSG_WARNING(msg)
 end if
 !
 ! mm(G,Gp) = (q+G) \cdot (q+Gp) n(G-Gp)
 ABI_MALLOC_OR_DIE(mm,(npwc,npwc), ierr)

 do ig=1,npwc
   if (qiszero) then
     ! To be discussed with Riad, here we should use the small q
     ! to be consistent and consider the limit q-->0
     gpq(:)=gvec(:,ig)
   else
     gpq(:)=gvec(:,ig)+qpt
   end if
   do igp=1,npwc
     if (qiszero) then
       gppq(:)=gvec(:,igp)
     else
       gppq(:)=gvec(:,igp)+qpt
     end if
     qpg_dot_qpgp=zero
     do ii=1,3
       qpg_dot_qpgp=qpg_dot_qpgp+&
&        ( gpq(1)*b1(ii) +gpq(2)*b2(ii) +gpq(3)*b3(ii))*&
&        (gppq(1)*b1(ii)+gppq(2)*b2(ii)+gppq(3)*b3(ii))
     end do
     mm(ig,igp)=rhogg(ig,igp)*qpg_dot_qpgp
   end do !igp
 end do !ig

 ABI_FREE(rhog_dp)
 ABI_FREE(rhog)
 ! === Now we have rhogg,rho0 ===
 !
 ! Calculate the dielectric matrix eigenvalues and vectors
 ! Use only the static epsm1 i.e., only the w=0 part (eps(:,:,1,:))
 ABI_MALLOC(eigval,(npwc))

 ABI_MALLOC_OR_DIE(eigvec,(npwc,npwc), ierr)

 ABI_MALLOC(zz,(npwc))
 zz=czero

 ABI_MALLOC(qplusg,(npwc))

 ! Store the susceptibility matrix in upper mode before calling zhpev.
 ABI_MALLOC_OR_DIE(matr,(npwc*(npwc+1)/2), ierr)

 idx=1
 do ii=1,npwc
   do jj=1,ii
     matr(idx)=epsm1(jj,ii); idx=idx+1
   end do
 end do

 ABI_MALLOC(zhpev2,(3*npwc-2))
 ABI_MALLOC(zhpev1,(2*npwc-1))

 call ZHPEV('V','U',npwc,matr,eigval,eigvec,npwc,zhpev1,zhpev2,ierr)
 ABI_FREE(matr)
 ABI_FREE(zhpev2)
 ABI_FREE(zhpev1)

 if (ierr<0) then
   write (msg,'(2a,i4,a)')&
&    ' Failed to calculate the eigenvalues and eigenvectors of the dielectric matrix ',ch10,&
&    ierr*(-1),'-th argument in the matrix has an illegal value. '
   MSG_ERROR(msg)
 end if

 if (ierr>0) then
   write(msg,'(3a,i4,2a)')&
&    ' Failed to calculate the eigenvalues and eigenvectors of the dielectric matrix ',ch10,&
&    ' the algorithm failed to converge; ierr = ', ierr,ch10,&
&    ' off-diagonal elements of an intermediate tridiagonal form did not converge to zero. '
   MSG_ERROR(msg)
 end if
 !
 ! Calculate the PPM parameters and the eigenpotentials needed for
 ! the calculation of the generalized overlap matrix
 ! Note: the eigenpotentials has to be calculated on the FFT (G-Gp) index
 !
 ! Save eigenvectors of \tilde\epsilon^{-1}
 ! MG well it is better to save \Theta otherwise
 ! we have to calculare \Theta for each band, spin, k-point but oh well
 eigtot=eigvec

 qlist(:,1) = qpt
 call cmod_qpg(1,1,qlist,npwc,gvec,gprimd,qplusg) !MG TODO here take care of small q
 !
 ! Basic Equation:
 !
 ! \Theta_{q,ii}(G)=\Psi_{q,ii}(G)/|q+G|
 ! where \Psi_{q,ii}(G) is the eigenvector of \tilde\epsilon^{-1}

 ! \tilde\omega_{ii,q}^2= 4\pi (1-eigenval(ii,q)))
 ! \sum_{G,Gp} \Theta^*_{q,ii}(G) (q+G)\cdot(q+Gp) n(G-Gp) \Theta_{q,ii}(Gp)

 do ii=1,npwc !DM band
   ! Calculate \Theta_{q,ii}(G)
   ! why the first element is not modified? if the problem is the small value of qplusg(1)
   ! we could multiply by sqrt(mod((q+G)(q+G'))) and then add the sing at the end
   if (qiszero)then
     eigvec(2:,ii)=eigvec(2:,ii)/qplusg(2:)
   else
     eigvec(:,ii)=eigvec(:,ii)/qplusg(:)
   end if
   do ig=1,npwc
     conjg_eig=CONJG(eigvec(ig,ii))
     do igp=1,npwc
       if(qiszero .and. ig==1 .and. igp==1)then
         zz(ii)=zz(ii)+conjg_eig*rhogg(ig,igp)*eigvec(igp,ii)
       else
         zz(ii)=zz(ii)+conjg_eig*mm(ig,igp)*eigvec(igp,ii)
       end if
     end do
   end do

   num=one-eigval(ii)
   if (num<=zero) then
!    here I think we should set bigomegatwsq=0 and omegatw to an arbitrary value
!    maybe we can output a warning TO BE discussed with Riad
     if (ABS(num)<1.0d-4) then
       num=1.0d-5
     else
       MSG_ERROR("One or more imaginary plasmon pole energies")
     end if
   end if

   omegatw(ii)=SQRT(4*pi*REAL(zz(ii))/num)
   ! this should be \alpha = 2\pi omegatw * (1-eigenval)
   ! MG check this, in the review I found a factor 2\pi, maybe it is reintroduced later
   bigomegatwsq(ii,1)=num*omegatw(ii)
 end do

 ABI_FREE(rhogg)
 ABI_FREE(mm)
 ABI_FREE(eigval)
 ABI_FREE(zz)
 ABI_FREE(eigvec)
 ABI_FREE(qplusg)

 call destroy_mpi_enreg(MPI_enreg_seq)

 write(msg,'(2a,f12.8,2a,3i5)')ch10,&
&  ' cppm3par : omega twiddle minval [eV]  = ',MINVAL(ABS(omegatw))*Ha_eV,ch10,&
&  '            omega twiddle min location = ',MINLOC(ABS(omegatw))
 call wrtout(std_out,msg,'COLL')

end subroutine cppm3par
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/cppm4par
!! NAME
!! cppm4par
!!
!! FUNCTION
!! Calculate the plasmon-pole parameters using Engel and Farid model (PRB47,15931,1993) [[cite:Engel1993]].
!! See also Quasiparticle Calculations in Solids [[cite:Aulbur2001]] p. 23
!!
!! INPUTS
!!  qpt(3)=Reduced coordinates of the q-point.
!!  epsm1(npwc,npwc)=symmetrized inverse dielectric matrix.
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  ngfftf(18)=contain all needed information about 3D fine FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  npwc=number of plane waves in epsm1
!!  rhor(nfftf)=charge density on the real space FFT grid
!!  gvec(3,npwc)=G vectors in reduced coordinated
!!
!! OUTPUT
!!  bigomegatwsq(npwc,npwc)=plasmon-pole strength
!!  omegatw(npwc)=plasmon-pole frequencies
!!
!! PARENTS
!!      m_ppmodel
!!
!! CHILDREN
!!
!! SOURCE

subroutine cppm4par(qpt,npwc,epsm1,ngfftf,gvec,gprimd,rhor,nfftf,bigomegatwsq,omegatw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,npwc
!arrays
 integer,intent(in) :: gvec(3,npwc),ngfftf(18)
 real(dp),intent(in) :: gprimd(3,3),qpt(3)
 real(dp),intent(in) :: rhor(nfftf)
 complex(gwpc),intent(in) :: epsm1(npwc,npwc)
 complex(gwpc),intent(out) :: bigomegatwsq(npwc,npwc),omegatw(npwc)

!Local variables-------------------------------
!scalars
 integer,parameter :: paral_kgb0=0
 integer :: ierr,ig,igp,ii,ngfft1,ngfft2,ngfft3,gmgp_idx
 real(dp) :: qpg_dot_qpgp
 character(len=500) :: msg
 character(len=80) :: bar
 type(MPI_type) :: MPI_enreg_seq
!arrays
 real(dp) :: b1(3),b2(3),b3(3),gppq(3),gpq(3),qlist(3,1)
 real(dp),allocatable :: eigval(:),qplusg(:),rhog_dp(:,:)
 real(dp),allocatable :: tmp_rhor(:)
 complex(dpc),allocatable :: chi(:,:),chitmp(:,:),chitmps(:,:),eigvec(:,:)
 complex(dpc),allocatable :: mm(:,:),mtemp(:,:),rhog(:)
 complex(dpc),allocatable :: rhogg(:,:),tmp1(:),zz2(:,:)

!*************************************************************************

 DBG_ENTER("COLL")

 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfftf(2),ngfftf(3),'all')

 b1=two_pi*gprimd(:,1)
 b2=two_pi*gprimd(:,2)
 b3=two_pi*gprimd(:,3)
 !
 ! === Calculate density in G space rhog(G) ===
 ABI_MALLOC(rhog_dp,(2,nfftf))
 ABI_MALLOC(rhog,(nfftf))

 ABI_MALLOC_OR_DIE(rhogg,(npwc,npwc), ierr)
 !
 ! Conduct FFT tho(r)-->rhog(G)
 ! FIXME this has to be fixed, rho(G) should be passed instead of doing FFT for each q
 ABI_MALLOC(tmp_rhor,(nfftf))
 tmp_rhor=rhor ! To avoid having to use intent(inout).
 call fourdp(1,rhog_dp,tmp_rhor,-1,MPI_enreg_seq,nfftf,1,ngfftf,0)
 ABI_FREE(tmp_rhor)

 rhog(1:nfftf)=CMPLX(rhog_dp(1,1:nfftf),rhog_dp(2,1:nfftf))
 ABI_FREE(rhog_dp)
 !
 ! Calculate the FFT index of each (G-Gp) vector and assign the value
 ! of the correspondent density simultaneously
 ngfft1=ngfftf(1)
 ngfft2=ngfftf(2)
 ngfft3=ngfftf(3)

 ierr=0
 do ig=1,npwc
   do igp=1,npwc
     gmgp_idx = g2ifft(gvec(:,ig)-gvec(:,igp),ngfftf)
     if (gmgp_idx/=0) then
       rhogg(ig,igp)=rhog(gmgp_idx)
     else
       ierr=ierr+1
       rhogg(ig,igp)=czero
     end if
   end do
 end do

 if (ierr/=0) then
   write(msg,'(a,i0,3a)')&
&   'Found ',ierr,' G1-G2 vectors falling outside the FFT box. ',ch10,&
&   'Enlarge the FFT mesh to get rid of this problem. '
   MSG_WARNING(msg)
 end if

 ABI_FREE(rhog)
 !
 ! Now we have rhogg, calculate the M matrix (q+G1).(q+G2) n(G1-G2)
 ABI_MALLOC_OR_DIE(mm,(npwc,npwc), ierr)

 do ig=1,npwc
   gpq(:)=gvec(:,ig)+qpt
   do igp=1,npwc
     gppq(:)=gvec(:,igp)+qpt
     qpg_dot_qpgp=zero
     do ii=1,3
       qpg_dot_qpgp=qpg_dot_qpgp+&
&        ( gpq(1)*b1(ii) +gpq(2)*b2(ii) +gpq(3)*b3(ii))*&
&        (gppq(1)*b1(ii)+gppq(2)*b2(ii)+gppq(3)*b3(ii))
     end do
     mm(ig,igp)=rhogg(ig,igp)*qpg_dot_qpgp
   end do !igp
 end do !ig

 !MG TODO too much memory in chi, we can do all this stuff inside a loop
 ABI_MALLOC_OR_DIE(chitmp,(npwc,npwc), ierr)

 ABI_MALLOC_OR_DIE(chi,(npwc,npwc), ierr)

 ABI_MALLOC(qplusg,(npwc))
 !
 ! Extract the full polarizability from \tilde \epsilon^{-1}
 ! \tilde\epsilon^{-1}_{G1 G2} = \delta_{G1 G2} + 4\pi \frac{\chi_{G1 G2}}{|q+G1| |q+G2|}
 chitmp(:,:)=epsm1(:,:)

 qlist(:,1) = qpt
 call cmod_qpg(1,1,qlist,npwc,gvec,gprimd,qplusg) !MG TODO here take care of small q

 do ig=1,npwc
   chitmp(ig,ig)=chitmp(ig,ig)-one
 end do
 do ig=1,npwc
   do igp=1,npwc
     chi(ig,igp)=chitmp(ig,igp)*qplusg(ig)*qplusg(igp)/four_pi
   end do
 end do
 ABI_FREE(chitmp)
 !
 ! === Solve chi*X = Lambda M*X where Lambda=-1/em(q)**2 ===
 ABI_MALLOC(eigval,(npwc))

 ABI_MALLOC_OR_DIE(eigvec,(npwc,npwc), ierr)
 ABI_MALLOC_OR_DIE(mtemp,(npwc,npwc), ierr)
 ABI_MALLOC_OR_DIE(chitmps,(npwc,npwc), ierr)
 !
 ! Copy chi and mm into working arrays
 chitmps(:,:)=chi(:,:)
 mtemp(:,:)=mm(:,:)

 call xhegv(1,"Vectors","Upper",npwc,chitmps,mtemp,eigval)
 !
 ! Eigenvectors are normalized as : X_i^* M X_j = \delta_{ij}
 eigvec(:,:)=chitmps(:,:)
 ABI_FREE(mtemp)
 ABI_FREE(chitmps)
 !
 ! === Calculate the plasmon pole parameters ===
 ABI_MALLOC(tmp1,(npwc))

 ABI_MALLOC_OR_DIE(zz2,(npwc,npwc), ierr)
 !
 ! good check:
 ! the lowest plasmon energy on gamma should be
 ! close to experimental plasma energy within an error of 10 %
 ! this error can be reduced further if one includes the non local
 ! commutators in the calculation of the polarizability at q==0
 zz2(:,:)=(0.0,0.0)

 qlist(:,1) = qpt
 call cmod_qpg(1,1,qlist,npwc,gvec,gprimd,qplusg) !MG TODO here take care of small q

 do ii=1,npwc
   !
   ! keeping in mind that the above matrix is negative definite
   ! we might have a small problem with the eigval that correspond to large G vectors
   ! i.e. DM band index, where the eigevalues become very small with
   ! possibility of being small positive numbers (due to numerical problems)
   ! thus as a caution one can use the following condition
   ! this will not affect the result since such a huge plasmon energy give almost zero
   ! contribution to the self correlation energy

   if (eigval(ii)>=zero) then
     eigval(ii) = -1.0d-4
     if (eigval(ii)>1.0d-3) then
       eigval(ii) = -1.0d-22
       write(msg,'(a,i6,a,es16.6)')&
&        ' Imaginary plasmon pole eigenenergy, eigenvector number ',ii,' with eigval',eigval(ii),ch10
       MSG_ERROR(msg)
     end if
   end if
   !
   ! === Save plasmon energies ===
   omegatw(ii)=SQRT(-1/eigval(ii))
   !
   ! Calculate and save scaled plasmon-pole eigenvectors
   ! defined as \sqrt{4\pi} \frac{Mx}{\sqrt{\tilde\omega} |q+G|}
   tmp1(:)=eigvec(:,ii)

   do ig=1,npwc
     do igp=1,npwc
       zz2(ig,ii)=zz2(ig,ii)+mm(ig,igp)*tmp1(igp) ! z--->y
     end do
     bigomegatwsq(ig,ii)=SQRT(four_pi)*zz2(ig,ii)/SQRT(omegatw(ii))
     bigomegatwsq(ig,ii)=bigomegatwsq(ig,ii)/qplusg(ig)
   end do

 end do

 ABI_FREE(tmp1)
 ABI_FREE(eigvec)
 ABI_FREE(eigval)
 ABI_FREE(zz2)

 call destroy_mpi_enreg(MPI_enreg_seq)

 ABI_FREE(qplusg)
 ABI_FREE(chi)
 ABI_FREE(rhogg)
 ABI_FREE(mm)

 bar=REPEAT('-',80)
 write(msg,'(3a)')bar,ch10,' plasmon energies vs q vector shown for lowest 10 bands                 '
 call wrtout(std_out,msg,'COLL')
 write(msg,'(2x,5x,10f7.3)')(REAL(omegatw(ig))*Ha_eV, ig=1,10)
 call wrtout(std_out,msg,'COLL')
 write(msg,'(a)')bar
 call wrtout(std_out,msg,'COLL')

 write(msg,'(2a,f12.8,2a,3i5)')ch10,&
&  ' cppm4par : omega twiddle minval [eV]  = ',MINVAL(ABS(omegatw))*Ha_eV,ch10,&
&  '            omega twiddle min location = ',MINLOC(ABS(omegatw))
 call wrtout(std_out,msg,'COLL')

 DBG_EXIT("COLL")

end subroutine cppm4par
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/cqratio
!! NAME
!! cqratio
!!
!! FUNCTION
!!  Calculate qratio(G,Gp,q)= (q+G)\cdot(q+Gp) / |q+G|^2 needed for Hybertsen-Louie and Plasmonpole model
!!
!! INPUTS
!!  npwc=number of planewaves considered (used for the correlation part)
!!  gvec(3,npwc)=reduced coordinates of the plane waves
!!  q(3)=coordinates of q points
!!  gmet(3,3)=metric in reciprocal space
!!  gprimd(3,3)=reciprocal lattice vectors
!!
!! OUTPUT
!!  qratio(npwc,npwc)=(q+G).(q+Gp)
!!
!! PARENTS
!!      m_ppmodel,mrgscr
!!
!! CHILDREN
!!
!! SOURCE

subroutine cqratio(npwc,gvec,q,gmet,gprimd,qratio)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwc
!arrays
 integer,intent(in) :: gvec(3,npwc)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),q(3)
 real(dp),intent(out) :: qratio(npwc,npwc)

!Local variables ------------------------------
!scalars
 integer :: ig,igp,ii
 real(dp) :: qpg_dot_qpgp
!arrays
 real(dp) :: b1(3),b2(3),b3(3),gppq(3),gpq(3),norm(npwc)

!************************************************************************

 b1=two_pi*gprimd(:,1); b2=two_pi*gprimd(:,2); b3=two_pi*gprimd(:,3)

 norm(:)=zero; qratio=zero

 !FIXME this loops have to be rewritten!!!!
 do ig=1,npwc
   gpq(:)=gvec(:,ig)+q
   norm(ig)=two_pi*SQRT(DOT_PRODUCT(gpq,MATMUL(gmet,gpq)))
!  norm(ig)=normv(gpq,gmet,'g')
 end do
 do ig=1,npwc
   gpq(:)=gvec(:,ig)+q
   do igp=1,npwc
     gppq(:)=gvec(:,igp)+q
     qpg_dot_qpgp=zero
!    qpg_dot_qpgp=vdotw(gpq,gppq,gmet,'g')
     do ii=1,3
       qpg_dot_qpgp=qpg_dot_qpgp+&
&        ( gpq(1)*b1(ii) +  gpq(2)*b2(ii) + gpq(3)*b3(ii))*&
&        (gppq(1)*b1(ii) + gppq(2)*b2(ii) +gppq(3)*b3(ii))
     end do

!    Now calculate qratio = (q+G).(q+Gp)/|q+G|^2
!    when |q+G|^2 and (q+G).(q+Gp) are both zero set (q+G).(q+Gp)/|q+G|^2 = 1
!    when |q+G|^2 is zero and |q+Gp| is not zero set (q+G).(q+Gp)/|q+G|^2 = 0
     if (norm(ig)<0.001) then
       if (norm(igp)<0.001) then     ! Case q=0 and G=Gp=0
         qratio(ig,igp)=one
       else                          ! Case q=0 and G=0 and Gp !=0
         qratio(ig,igp)=zero
       end if
     else if (norm(igp)<0.001) then  ! Case q=0 and G= !0 and Gp=0
       qratio(ig,igp)=zero
     else
       qratio(ig,igp)=qpg_dot_qpgp/norm(ig)**2
     end if

   end do
 end do

end subroutine cqratio
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/calc_sig_ppm
!!
!! NAME
!! calc_sig_ppm
!!
!! FUNCTION
!!  Calculate the contribution to self-energy operator using a plasmon-pole model.
!!
!! INPUTS
!!  nomega=Number of frequencies.
!!  nspinor=Number of spinorial components.
!!  npwc=Number of G vectors in the plasmon pole.
!!  npwx=number of G vectors in rhotwgp, i.e. no. of G-vectors for the exchange part.
!!  theta_mu_minus_e0i= $\theta(\mu-\epsilon_{k-q,b1,s}), defines if the state is occupied or not
!!  zcut=Small imaginary part to avoid the divergence. (see related input variable)
!!  omegame0i(nomega)=Frequencies used to evaluate \Sigma_c ($\omega$ - $\epsilon_i)$
!!  otq(npwc,dm2_otq)=Plasmon pole parameters for this q-point.
!!  PPm<ppmodel_t>=structure gathering info on the Plasmon-pole technique.
!!     %model=type plasmon pole model
!!     %dm2_botsq= 1 if model==3, =npwc if model== 4, 1 for all the other cases
!!     %dm2_otq= 1 if model==3, =1    if model== 4, 1 for all the other cases
!!     %dm_eig=npwc if model=3, 0 otherwise
!!  botsq(npwc,dm2_botsq)=Plasmon pole parameters for this q-point.
!!  eig(dm_eig,dm_eig)=The eigvectors of the symmetrized inverse dielectric matrix for this q point
!!   (first index for G, second index for bands)
!!  rhotwgp(npwx)=oscillator matrix elements divided by |q+G| i.e.
!!    $\frac{\langle b1 k-q s | e^{-i(q+G)r | b2 k s \rangle}{|q+G|}$
!!
!! OUTPUT
!!  sigcme(nomega) (to be described), only relevant if ppm3 or ppm4
!!  ket(npwc,nomega):
!!  === model==1,2 ====
!!    ket(G,omega) = Sum_G2       conjg(rhotw(G)) * Omega(G,G2) * rhotw(G2)
!!                            ---------------------------------------------------
!!                             2 omegatw(G,G2) (omega-E_i + omegatw(G,G2)(2f-1))
!!
!! PARENTS
!!      calc_sigc_me
!!
!! CHILDREN
!!
!! SOURCE

subroutine calc_sig_ppm(PPm,nspinor,npwc,nomega,rhotwgp,botsq,otq,&
& omegame0i,zcut,theta_mu_minus_e0i,eig,npwx,ket,sigcme)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc,npwx,nspinor
 real(dp),intent(in) :: theta_mu_minus_e0i,zcut
 type(ppmodel_t),intent(in) :: PPm
!arrays
 real(dp),intent(in) :: omegame0i(nomega)
 complex(gwpc),intent(in) :: botsq(npwc,PPm%dm2_botsq),eig(PPm%dm_eig,PPm%dm_eig),otq(npwc,PPm%dm2_otq)
 complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 complex(gwpc),intent(inout) :: ket(npwc*nspinor,nomega)
 complex(gwpc),intent(out) :: sigcme(nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,ii,ios,ispinor,spadc,spadx
 real(dp) :: den,ff,inv_den,omegame0i_io,otw,twofm1,twofm1_zcut
 complex(gwpc) :: ct,num,numf,rhotwgdp_igp
 logical :: fully_occupied,totally_empty
 !character(len=500) :: msg
!arrays
 complex(gwpc),allocatable :: rhotwgdpcc(:)

!*************************************************************************

 SELECT CASE (PPm%model)
 CASE (PPM_GODBY_NEEDS, PPM_HYBERTSEN_LOUIE)
   fully_occupied =(ABS(theta_mu_minus_e0i-one)<0.001)
   totally_empty  =(ABS(theta_mu_minus_e0i    )<0.001)

   do ispinor=1,nspinor
     spadx=(ispinor-1)*npwx; spadc=(ispinor-1)*npwc

     if (.not.totally_empty) then
       ! \Bomega^2_{G1G2}/\omegat_{G1G2} M_{G1,G2}. \theta(\mu-e_s) / (\omega+\omegat_{G1G2}-e_s-i\delta)
       twofm1_zcut=zcut
!$omp parallel do private(omegame0i_io, rhotwgdp_igp, otw, num, den)
       do ios=1,nomega
         omegame0i_io=omegame0i(ios)
         do igp=1,npwc
           rhotwgdp_igp=rhotwgp(spadx+igp)
           do ig=1,npwc
             otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
             num = botsq(ig,igp)*rhotwgdp_igp
             den = omegame0i_io+otw
             if (den**2>zcut**2) then
               ket(spadc+ig,ios)=ket(spadc+ig,ios) + num/(den*otw)*theta_mu_minus_e0i
             else
               ket(spadc+ig,ios)=ket(spadc+ig,ios) + &
&                num*CMPLX(den,twofm1_zcut)/((den**2+twofm1_zcut**2)*otw)*theta_mu_minus_e0i
             end if
           end do !ig
         end do !igp
       end do !ios
     end if !not totally empty

     if (.not.(fully_occupied)) then
       ! \Bomega^2_{G1G2}/\omegat_{G1G2} M_{G1,G2}. \theta(e_s-\mu) / (\omega-\omegat_{G1G2}-e_s+i\delta)
       twofm1_zcut=-zcut
!$omp parallel do private(omegame0i_io, rhotwgdp_igp, otw, num, den)
       do ios=1,nomega
         omegame0i_io=omegame0i(ios)
         do igp=1,npwc
           rhotwgdp_igp=rhotwgp(spadx+igp)
           do ig=1,npwc
             otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
             num = botsq(ig,igp)*rhotwgdp_igp
             den=omegame0i_io-otw
             if (den**2>zcut**2) then
               ket(spadc+ig,ios)=ket(spadc+ig,ios) + num/(den*otw)*(one-theta_mu_minus_e0i)
             else
               ket(spadc+ig,ios)=ket(spadc+ig,ios) + &
&                num*CMPLX(den,twofm1_zcut)/((den**2+twofm1_zcut**2)*otw)*(one-theta_mu_minus_e0i)
             end if
           end do !ig
         end do !igp
         !
       end do !ios
     end if !not fully occupied

   end do !ispinor

   ket=ket*half

 CASE (PPM_LINDEN_HORSH,PPM_ENGEL_FARID)
   ABI_CHECK(nspinor == 1, "nspinor/=1 not allowed")

   ! rho-twiddle(G) is formed, introduce rhotwgdpcc, for speed reason
   ABI_MALLOC(rhotwgdpcc,(npwx))

   ff=theta_mu_minus_e0i      ! occupation number f (include poles if ...)
   twofm1=two*ff-one          ! 2f-1
   twofm1_zcut=twofm1*zcut
   rhotwgdpcc(:)=CONJG(rhotwgp(:))

   do ios=1,nomega
     omegame0i_io=omegame0i(ios)
     ct=czero_gw
     do ii=1,npwc ! Loop over the DM bands
       num=czero_gw

       SELECT CASE (PPm%model)
       CASE (PPM_LINDEN_HORSH)
         ! Calculate \beta (eq. 106 pag 47)
         do ig=1,npwc
           num=num+rhotwgdpcc(ig)*eig(ig,ii)
         end do
         numf=num*CONJG(num) !MG this means that we cannot do SCGW
         numf=numf*botsq(ii,1)

       CASE (PPM_ENGEL_FARID)
         do ig=1,npwc
           num=num+rhotwgdpcc(ig)*botsq(ig,ii)
         end do
         numf=num*CONJG(num) !MG this means that we cannot do SCGW

       CASE DEFAULT
         MSG_ERROR("Wrong PPm%model")
       END SELECT

       !numf=num*CONJG(num) !MG this means that we cannot do SCGW
       !if (PPm%model==3) numf=numf*botsq(ii,1)

       otw=DBLE(otq(ii,1)) ! in principle otw -> otw - ieta
       den=omegame0i_io+otw*twofm1

       if (den**2>zcut**2) then
         inv_den=one/den
         ct=ct+numf*inv_den
       else
         inv_den=one/((den**2+twofm1_zcut**2))
         ct=ct+numf*CMPLX(den,twofm1_zcut)*inv_den
       end if

     end do !ii DM bands
     sigcme(ios)=ct*half
   end do !ios
   ABI_FREE(rhotwgdpcc)

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong PPm%model:',itoa(PPm%model)))
 END SELECT

end subroutine calc_sig_ppm
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_symmetrizer
!! NAME
!!  ppm_symmetrizer
!!
!! FUNCTION
!!  Symmetrize the plasmonpole matrix elements in the full BZ zone.
!!
!! INPUTS
!!  iq_bz=Index of the q-point in the BZ where the PPmodel parameters are wanted.
!!  Gsph<gsphere_t>=data related to the G-sphere.
!!  Cryst<crystal_t>=Info on the unit cell and crystal symmetries.
!!  Qmesh<kmesh_t>=the q-mesh used for the inverse dielectric matrix
!!  iq_ibz=Index of the q-point in the BZ.
!!  npwe=number of G vectors for the correlation part
!!  nomega=number of frequencies in $\epsilon^{-1}$
!!  omega=frequencies in epsm1_ggw
!!  epsm1_ggw(npwe,npwe,nomega)=the inverse dielctric matrix
!!  ngfftf(18)=contain all needed information about the 3D fine FFT mesh, see ~abinit/doc/variables/vargs.htm#ngfft
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  nfftf=the number of points in the FFT mesh (for this processor)
!!  rhor_tot(nfftf)=the total charge in real space
!!
!! SIDE EFFECTS
!!  PPm<ppmodel_t>=data type containing information on the plasmonpole technique.
!!  Internal tables are modified so that they (point|store) the plasmon-pole parameters
!!  for the specified q-point in the BZ.
!!
!! PARENTS
!!      m_screen
!!
!! CHILDREN
!!
!! SOURCE

subroutine ppm_symmetrizer(PPm,iq_bz,Cryst,Qmesh,Gsph,npwe,nomega,omega,epsm1_ggw,nfftf,ngfftf,rhor_tot)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,npwe,nomega,iq_bz
 type(crystal_t),intent(in) :: Cryst
 type(gsphere_t),intent(in) :: Gsph
 type(kmesh_t),intent(in) :: Qmesh
 type(ppmodel_t),target,intent(inout) :: PPm
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(in) :: rhor_tot(nfftf)
 complex(dpc),intent(in) :: omega(nomega)
 complex(gwpc),intent(in) :: epsm1_ggw(npwe,npwe,nomega)

!Local variables-------------------------------
!scalars
 integer :: iq_ibz,itim_q,isym_q,iq_curr,ierr
 logical :: q_isirred
 !character(len=500) :: msg
!arrays
 real(dp) :: qbz(3)
! *********************************************************************

 !@ppmodel_t
 !
 ! Save the index of the q-point in the BZ for checking purpose.
 PPm%iq_bz = iq_bz
 !
 ! Here there is a problem with the small q, still cannot use BZ methods
 !iq_ibz=Qmesh%tab(iq_bz)
 !isym_q=Qmesh%tabo(iq_bz)
 !itim_q=(3-Qmesh%tabi(iq_bz))/2

 call get_bz_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q,isirred=q_isirred)
 iq_curr=iq_ibz; if (PPm%mqmem==0) iq_curr=1
 !
 ! =======================================================
 ! ==== Branching for in-core or out-of-core solution ====
 ! =======================================================
 !
 if (PPm%has_q(iq_ibz)==PPM_NOTAB) then
   ! Allocate the tables for this q_ibz
   call ppm_mallocq(PPm,iq_ibz)
 end if

 if (PPm%has_q(iq_ibz)==PPM_TAB_ALLOCATED) then
   ! Calculate the ppmodel tables for this q_ibz
   call new_setup_ppmodel(PPm,iq_ibz,Cryst,Qmesh,npwe,nomega,omega,epsm1_ggw,nfftf,Gsph%gvec,ngfftf,rhor_tot) !Optional
 end if

 if (q_isirred) then
   ! Symmetrization is not needed.
   if (PPm%bigomegatwsq_qbz_stat==PPM_ISALLOCATED)  then
     ABI_FREE(PPm%bigomegatwsq_qbz%vals)
   end if
   if (PPm%omegatw_qbz_stat==PPM_ISALLOCATED)  then
     ABI_FREE(PPm%omegatw_qbz%vals)
   end if
   if (PPm%eigpot_qbz_stat==PPM_ISALLOCATED)  then
     ABI_FREE(PPm%eigpot_qbz%vals)
   end if

   ! Point the data in memory and change the status.
   PPm%bigomegatwsq_qbz => PPm%bigomegatwsq(iq_ibz)
   PPm%bigomegatwsq_qbz_stat = PPM_ISPOINTER

   PPm%omegatw_qbz => PPm%omegatw(iq_ibz)
   PPm%omegatw_qbz_stat = PPM_ISPOINTER

   PPm%eigpot_qbz => PPm%eigpot(iq_ibz)
   PPm%eigpot_qbz_stat = PPM_ISPOINTER
 else
   ! Allocate memory if not done yet.
   if (PPm%bigomegatwsq_qbz_stat==PPM_ISPOINTER) then
      nullify(PPm%bigomegatwsq_qbz)
      ABI_MALLOC_OR_DIE(PPm%bigomegatwsq_qbz%vals, (PPm%npwc,PPm%dm2_botsq), ierr)
      PPm%bigomegatwsq_qbz_stat=PPM_ISALLOCATED
   end if

   if (PPm%omegatw_qbz_stat==PPM_ISPOINTER) then
      nullify(PPm%omegatw_qbz)
      ABI_MALLOC_OR_DIE(PPm%omegatw_qbz%vals,(PPm%npwc,PPm%dm2_otq), ierr)
      PPm%omegatw_qbz_stat=PPM_ISALLOCATED
   end if

   if (PPm%eigpot_qbz_stat==PPM_ISPOINTER) then
     nullify(PPm%eigpot_qbz)
     ABI_MALLOC_OR_DIE(PPm%eigpot_qbz%vals,(PPm%dm_eig,PPm%dm_eig), ierr)
     PPm%eigpot_qbz_stat=PPM_ISALLOCATED
   end if
   !
   ! Calculate new table for this q-point in the BZ.
   ! Beware: Dimensions should not change.
   !botsq => PPm%bigomegatwsq_qbz%vals
   !otq   => PPm%omegatw_qbz%vals
   !eig   => PPm%eigpot_qbz%vals

   call ppm_get_qbz(PPm,Gsph,Qmesh,iq_bz,PPm%bigomegatwsq_qbz%vals,&
&    PPm%omegatw_qbz%vals,PPm%eigpot_qbz%vals)

   ! Release the table in the IBZ if required.
   if (.not.PPm%keep_q(iq_ibz)) call ppm_table_free(PPm,iq_ibz)
 end if

end subroutine ppm_symmetrizer
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/new_setup_ppmodel
!! NAME
!! new_setup_ppmodel
!!
!! FUNCTION
!!  Initialize some values of several arrays of the PPm datastructure
!!  that are used in case of plasmonpole calculations
!!  Just a wrapper around different plasmonpole routines.
!!
!! INPUTS
!!  Cryst<crystal_t>=Info on the unit cell and crystal symmetries.
!!  Qmesh<kmesh_t>=the q-mesh used for the inverse dielectric matrix
!!  iq_ibz=Index of the q-point in the BZ.
!!  npwe=number of G vectors for the correlation part
!!  nomega=number of frequencies in $\epsilon^{-1}$
!!  omega=frequencies in epsm1_ggw
!!  epsm1_ggw(npwe,npwe,nomega)=the inverse dielctric matrix
!!  ngfftf(18)=contain all needed information about the 3D fine FFT mesh, see ~abinit/doc/variables/vargs.htm#ngfft
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  nfftf=the number of points in the FFT mesh (for this processor)
!!  rhor_tot(nfftf)=the total charge in real space
!!  PPm<ppmodel_t>:
!!
!! SIDE EFFECTS
!!  PPm<ppmodel_t>:
!!  == if ppmodel 1 or 2 ==
!!   %omegatw and %bigomegatwsq
!!  == if ppmodel 3 ==
!!   %omegatw, %bigomegatwsq and %eigpot
!!  == if ppmodel 4 ==
!!   %omegatw and %bigomegatwsq
!!
!! NOTES
!! * FFT parallelism not implemented.
!! * TODO: rhor_tot should be replaced by rhog_tot
!!
!! PARENTS
!!      m_ppmodel,m_screen
!!
!! CHILDREN
!!
!! SOURCE

subroutine new_setup_ppmodel(PPm,iq_ibz,Cryst,Qmesh,npwe,nomega,omega,epsm1_ggw,nfftf,gvec,ngfftf,rhor_tot)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,npwe,nomega,iq_ibz
 type(kmesh_t),intent(in) :: Qmesh
 type(crystal_t),intent(in) :: Cryst
 type(ppmodel_t),intent(inout) :: PPm
!arrays
 integer,intent(in) :: gvec(3,npwe),ngfftf(18)
 real(dp),intent(in) :: rhor_tot(nfftf)
 complex(dpc),intent(in) :: omega(nomega)
 complex(gwpc),intent(in) :: epsm1_ggw(npwe,npwe,nomega)

!Local variables-------------------------------
!scalars
 real(dp) :: n_at_G_zero
 character(len=500) :: msg
!scalars
 real(dp) :: qpt(3)
! *************************************************************************

 DBG_ENTER("COLL")

 !@ppmodel_t
 !
 if (PPm%has_q(iq_ibz) /= PPM_TAB_ALLOCATED) then
   MSG_ERROR("ppmodel tables are not allocated")
 end if

 qpt = Qmesh%ibz(:,iq_ibz)
 PPm%has_q(iq_ibz) = PPM_TAB_STORED
 !
 ! Calculate plasmonpole parameters
 ! TODO ppmodel==1 by default, should be set to 0 if AC and CD
 SELECT CASE (PPm%model)

 CASE (PPM_NONE)
   MSG_COMMENT('Skipping Plasmonpole model calculation')

 CASE (PPM_GODBY_NEEDS)
   ! Note: the q-dependency enters only through epsilon^-1.
   call cppm1par(npwe,nomega,omega,PPm%drude_plsmf,epsm1_ggw,&
&     PPm%omegatw(iq_ibz)%vals,PPm%bigomegatwsq(iq_ibz)%vals)

 CASE (PPM_HYBERTSEN_LOUIE)
   call cppm2par(qpt,npwe,epsm1_ggw(:,:,1),ngfftf,gvec,Cryst%gprimd,rhor_tot,nfftf,Cryst%gmet,&
&    PPm%bigomegatwsq(iq_ibz)%vals,PPm%omegatw(iq_ibz)%vals,PPm%invalid_freq)
   !
   ! Quick-and-dirty change of the plasma frequency. Never executed in standard runs.
   if (PPm%force_plsmf>tol6) then ! Integrate the real-space density
      n_at_G_zero = SUM(rhor_tot(:))/nfftf
      ! Change the prefactor
      write(msg,'(2(a,es16.8))') 'Forced ppmfreq:',PPm%force_plsmf*Ha_eV,' nelec/ucvol:',n_at_G_zero
      MSG_WARNING(msg)
      PPm%force_plsmf = (PPm%force_plsmf**2)/(four_pi*n_at_G_zero)
      PPm%bigomegatwsq(iq_ibz)%vals = PPm%force_plsmf * PPm%bigomegatwsq(iq_ibz)%vals
      PPm%omegatw(iq_ibz)%vals      = PPm%force_plsmf * PPm%omegatw(iq_ibz)%vals
      write(msg,'(a,es16.8)') 'Plasma frequency forced in HL ppmodel, new prefactor is:',PPm%force_plsmf
      MSG_WARNING(msg)
   end if

 CASE (PPM_LINDEN_HORSH)
   ! TODO Check better double precision, this routine is in a messy state
   call cppm3par(qpt,npwe,epsm1_ggw(:,:,1),ngfftf,gvec,Cryst%gprimd,rhor_tot,nfftf,&
&                PPm%bigomegatwsq(iq_ibz)%vals,PPm%omegatw(iq_ibz)%vals(:,1),PPm%eigpot(iq_ibz)%vals)

 CASE (PPM_ENGEL_FARID)  ! TODO Check better double precision, this routine is in a messy state
   if ((ALL(ABS(qpt)<1.0e-3))) qpt = GW_Q0_DEFAULT ! FIXME

   call cppm4par(qpt,npwe,epsm1_ggw(:,:,1),ngfftf,gvec,Cryst%gprimd,rhor_tot,nfftf,&
&    PPm%bigomegatwsq(iq_ibz)%vals,PPm%omegatw(iq_ibz)%vals(:,1))

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong PPm%model:',itoa(PPm%model)))
 END SELECT

 DBG_EXIT("COLL")

end subroutine new_setup_ppmodel
!!***

!----------------------------------------------------------------------

!!****f* m_ppmodel/ppm_times_ket
!!
!! NAME
!! ppm_times_ket
!!
!! FUNCTION
!!  Calculate the contribution to self-energy operator using a plasmon-pole model.
!!
!! INPUTS
!!  nomega=Number of frequencies.
!!  nspinor=Number of spinorial components.
!!  npwc=Number of G vectors in the plasmon pole.
!!  npwx=number of G vectors in rhotwgp, i.e. no. of G-vectors for the exchange part.
!!  theta_mu_minus_e0i= $\theta(\mu-\epsilon_{k-q,b1,s}), defines if the state is occupied or not
!!  zcut=Small imaginary part to avoid the divergence. (see related input variable)
!!  omegame0i(nomega)=Frequencies used to evaluate \Sigma_c ($\omega$ - $\epsilon_i)$
!!  PPm<ppmodel_t>=structure gathering info on the Plasmon-pole technique.
!!     %model=type plasmon pole model
!!     %dm2_botsq= 1 if model==3, =npwc if model== 4, 1 for all the other cases
!!     %dm2_otq= 1 if model==3, =1    if model== 4, 1 for all the other cases
!!     %dm_eig=npwc if model=3, 0 otherwise
!!     %botsq(npwc,dm2_botsq)=Plasmon pole parameters for this q-point.
!!     %eig(dm_eig,dm_eig)=The eigvectors of the symmetrized inverse dielectric matrix for this q point
!!       (first index for G, second index for bands)
!!     %otq(npwc,dm2_otq)=Plasmon pole parameters for this q-point.
!!  rhotwgp(npwx)=oscillator matrix elements divided by |q+G| i.e.
!!    $\frac{\langle b1 k-q s | e^{-i(q+G)r | b2 k s \rangle}{|q+G|}$
!!
!! OUTPUT
!!  sigcme(nomega) (to be described), only relevant if ppm3 or ppm4
!!  ket(npwc,nomega):
!!   === model==1,2 ====
!!
!!   ket(G,omega) = Sum_G2       conjg(rhotw(G)) * Omega(G,G2) * rhotw(G2)
!!                          ---------------------------------------------------
!!                            2 omegatw(G,G2) (omega-E_i + omegatw(G,G2)(2f-1))
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ppm_times_ket(PPm,nspinor,npwc,nomega,rhotwgp,omegame0i,zcut,theta_mu_minus_e0i,npwx,ket,sigcme)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc,npwx,nspinor
 real(dp),intent(in) :: theta_mu_minus_e0i,zcut
 type(ppmodel_t),intent(in) :: PPm
!arrays
 real(dp),intent(in) :: omegame0i(nomega)
 complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 complex(gwpc),intent(inout) :: ket(npwc*nspinor,nomega)
 complex(gwpc),intent(out) :: sigcme(nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,ii,ios,ispinor,spadc,spadx
 real(dp) :: den,ff,inv_den,omegame0i_io,otw,twofm1,twofm1_zcut
 complex(gwpc) :: ct,num,numf,rhotwgdp_igp
 logical :: fully_occupied,totally_empty
!arrays
 complex(gwpc),allocatable :: rhotwgdpcc(:)
 complex(gwpc), ABI_CONTIGUOUS pointer :: botsq(:,:),eig(:,:),otq(:,:)

!*************************************************************************

 SELECT CASE (PPm%model)
 CASE (PPM_GODBY_NEEDS, PPM_HYBERTSEN_LOUIE)
   botsq  => PPm%bigomegatwsq_qbz%vals !(1:npwc,1:PPm%dm2_botsq)
   otq    => PPm%omegatw_qbz%vals      !(1:npwc,1:PPm%dm2_otq)
   fully_occupied =(ABS(theta_mu_minus_e0i-one)<0.001)
   totally_empty  =(ABS(theta_mu_minus_e0i    )<0.001)

   do ispinor=1,nspinor
     spadx=(ispinor-1)*npwx; spadc=(ispinor-1)*npwc

     if (.not.totally_empty) then
       ! \Bomega^2_{G1G2}/\omegat_{G1G2} M_{G1,G2}. \theta(\mu-e_s) / (\omega+\omegat_{G1G2}-e_s-i\delta)
       twofm1_zcut=zcut
!$omp parallel do private(omegame0i_io, rhotwgdp_igp, otw, num, den)
       do ios=1,nomega
         omegame0i_io=omegame0i(ios)
         do igp=1,npwc
           rhotwgdp_igp=rhotwgp(spadx+igp)
           do ig=1,npwc
             otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
             num = botsq(ig,igp)*rhotwgdp_igp
             den = omegame0i_io+otw
             if (den**2>zcut**2) then
               ket(spadc+ig,ios)=ket(spadc+ig,ios) + num/(den*otw)*theta_mu_minus_e0i
             else
               ket(spadc+ig,ios)=ket(spadc+ig,ios) + &
&                num*CMPLX(den,twofm1_zcut)/((den**2+twofm1_zcut**2)*otw)*theta_mu_minus_e0i
             end if
           end do !ig
         end do !igp
       end do !ios
     end if !not totally empty

     if (.not.(fully_occupied)) then
       ! \Bomega^2_{G1G2}/\omegat_{G1G2} M_{G1,G2}. \theta(e_s-\mu) / (\omega-\omegat_{G1G2}-e_s+i\delta)
       twofm1_zcut=-zcut
!$omp parallel do private(omegame0i_io, rhotwgdp_igp, otw, num, den)
       do ios=1,nomega
         omegame0i_io=omegame0i(ios)
         do igp=1,npwc
           rhotwgdp_igp=rhotwgp(spadx+igp)
           do ig=1,npwc
             otw=DBLE(otq(ig,igp)) !in principle otw -> otw - ieta
             num = botsq(ig,igp)*rhotwgdp_igp
             den=omegame0i_io-otw
             if (den**2>zcut**2) then
               ket(spadc+ig,ios)=ket(spadc+ig,ios) + num/(den*otw)*(one-theta_mu_minus_e0i)
             else
               ket(spadc+ig,ios)=ket(spadc+ig,ios) + &
&               num*CMPLX(den,twofm1_zcut)/((den**2+twofm1_zcut**2)*otw)*(one-theta_mu_minus_e0i)
             end if
           end do !ig
         end do !igp
       end do !ios
     end if !not fully occupied

   end do !ispinor

!$omp parallel workshare
   ket=ket*half
!$omp end parallel workshare

 CASE (PPM_LINDEN_HORSH,PPM_ENGEL_FARID)
   ABI_CHECK(nspinor == 1, "nspinor/=1 not allowed")

   botsq => PPm%bigomegatwsq_qbz%vals !(1:npwc,1:PPm%dm2_botsq)
   otq   => PPm%omegatw_qbz%vals      !(1:npwc,1:PPm%dm2_otq)
   eig   => PPm%eigpot_qbz%vals       !(1:PPm%dm_eig,1:PPm%dm_eig)

   ! rho-twiddle(G) is formed, introduce rhotwgdpcc, for speed reason
   ABI_MALLOC(rhotwgdpcc,(npwx))

   ff=theta_mu_minus_e0i      ! occupation number f (include poles if ...)
   twofm1=two*ff-one          ! 2f-1
   twofm1_zcut=twofm1*zcut
   rhotwgdpcc(:)=CONJG(rhotwgp(:))

   do ios=1,nomega
     omegame0i_io=omegame0i(ios)
     ct=czero_gw
     do ii=1,npwc ! Loop over the DM bands
       num=czero_gw

       SELECT CASE (PPm%model)
       CASE (PPM_LINDEN_HORSH)
         ! Calculate \beta (eq. 106 pag 47)
         do ig=1,npwc
           num=num+rhotwgdpcc(ig)*eig(ig,ii)
         end do
         numf=num*CONJG(num) !MG this means that we cannot do SCGW
         numf=numf*botsq(ii,1)

       CASE (PPM_ENGEL_FARID)
         do ig=1,npwc
           num=num+rhotwgdpcc(ig)*botsq(ig,ii)
         end do
         numf=num*CONJG(num) !MG this means that we cannot do SCGW

       CASE DEFAULT
         MSG_ERROR("Wrong PPm%model")
       END SELECT

       !numf=num*CONJG(num) !MG this means that we cannot do SCGW
       !if (PPm%model==3) numf=numf*botsq(ii,1)

       otw=DBLE(otq(ii,1)) ! in principle otw -> otw - ieta
       den=omegame0i_io+otw*twofm1

       if (den**2>zcut**2) then
         inv_den=one/den
         ct=ct+numf*inv_den
       else
         inv_den=one/((den**2+twofm1_zcut**2))
         ct=ct+numf*CMPLX(den,twofm1_zcut)*inv_den
       end if

     end do !ii DM bands
     sigcme(ios)=ct*half
   end do !ios
   ABI_FREE(rhotwgdpcc)

 CASE DEFAULT
   MSG_BUG(sjoin('Wrong PPm%model:',itoa(PPm%model)))
 END SELECT

end subroutine ppm_times_ket
!!***

!----------------------------------------------------------------------

END MODULE m_ppmodel
!!***
