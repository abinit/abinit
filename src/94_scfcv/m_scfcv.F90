!!****m* ABINIT/m_scfcv
!! NAME
!!  m_scfcv
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (JB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

module m_scfcv

 use defs_basis
 use defs_wvltypes
 use defs_rectypes
 use m_abicore
 use m_errors
 use m_wffile
 use m_rec
 use m_efield
 use m_entropyDMFT
 use m_hdr
 use m_dtfil

 use defs_datatypes,     only : pseudopotential_type
 use defs_abitypes,      only : MPI_type
 use m_scf_history,      only : scf_history_type
 use m_results_gs ,      only : results_gs_type
 use m_electronpositron, only : electronpositron_type
 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_pawcprj,          only : pawcprj_type
 use m_pawrhoij,         only : pawrhoij_type
 use m_pawfgr,           only : pawfgr_type
 use m_paw_dmft,         only : paw_dmft_type
 use m_paw_uj,           only : macro_uj_type
 use m_orbmag,           only : orbmag_type
 use m_data4entropyDMFT, only : data4entropyDMFT_t, data4entropyDMFT_init, data4entropyDMFT_destroy
 use m_scfcv_core,       only : scfcv_core

 implicit none

 private

 public :: scfcv_init
 public :: scfcv_destroy
 public :: scfcv_run
!!***

! *************************************************************************

!!****t* m_scfcv/scfcv_t
!! NAME
!!  scfcv_t
!!
!! FUNCTION
!!  This structured datatype contains the necessary data
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (JB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

 type, public :: scfcv_t
  !scalars
   integer,pointer :: mcg => null()
   integer,pointer :: mcprj => null()
   integer,pointer :: my_natom => null()
   integer,pointer :: ndtpawuj => null()
   integer,pointer :: pwind_alloc => null()
   integer,pointer :: initialized => null()
   integer,pointer :: nfftf => null()
   real(dp),pointer :: cpus => null()
   real(dp),pointer :: ecore  => null()
   real(dp),pointer :: fatvshift => null()
   type(pawang_type),pointer :: pawang => null()
   type(pseudopotential_type),pointer :: psps => null()
   type(MPI_type),pointer :: mpi_enreg => null()
   type(datafiles_type),pointer :: dtfil => null()
   type(dataset_type),pointer :: dtset => null()
   type(efield_type),pointer :: dtefield => null()
   type(orbmag_type),pointer :: dtorbmag => null()
   type(electronpositron_type),pointer :: electronpositron => null()
   type(hdr_type),pointer :: hdr => null()
   type(pawfgr_type),pointer :: pawfgr => null()
   type(recursion_type),pointer :: rec_set => null()
   type(results_gs_type),pointer :: results_gs => null()
   type(scf_history_type),pointer :: scf_history => null()
   type(wffile_type),pointer :: wffnew => null()
   type(wffile_type),pointer :: wffnow => null()
   type(wvl_data),pointer :: wvl => null()
   type(paw_dmft_type), pointer :: paw_dmft => null()

   !arrays
   integer,pointer :: atindx(:) => null()
   integer,pointer :: atindx1(:) => null()
   integer, pointer :: irrzon(:,:,:) => null()
   integer, pointer :: symrec(:,:,:) => null()
   integer,pointer :: indsym(:,:,:) => null()
   !no_abirules
   integer, pointer :: kg(:,:) => null()
   integer, pointer :: nattyp(:) => null()
   integer, pointer :: npwarr(:) => null()
   integer, pointer :: pwind(:,:,:) => null()
   real(dp), pointer :: dmatpawu(:,:,:,:) => null()
   real(dp), pointer :: phnons(:,:,:) => null()
   real(dp), pointer :: pwnsfac(:,:) => null()
   real(dp), pointer :: ylm(:,:) => null()
   real(dp), pointer :: ylmgr(:,:,:) => null()
   real(dp), pointer :: cg(:,:) => null()
   real(dp), pointer :: eigen(:) => null()
   real(dp), pointer :: occ(:) => null()
   !real(dp), pointer :: rprimd(:,:) => null()
   !real(dp), pointer :: rhog(:,:) => null()
   !real(dp), pointer :: rhor(:,:) => null()
   real(dp), pointer :: taug(:,:) => null()
   real(dp), pointer :: taur(:,:) => null()
   real(dp), pointer :: resid(:) => null()
   type(pawrad_type), pointer :: pawrad(:) => null()
   type(pawtab_type), pointer :: pawtab(:) => null()
   type(macro_uj_type),pointer :: dtpawuj(:) => null()
   type(pawrhoij_type), pointer :: pawrhoij(:) => null()
   type(pawcprj_type),pointer :: cprj(:,:) => null()
   ! PRIVATE ATTRIBUTS
   type(entropyDMFT_t) ABI_PRIVATE :: entropyDMFT
 end type scfcv_t
!!***

! *************************************************************************

contains
!!***


!!****f* ABINIT/m_scfcv/scfcv_init
!! NAME
!!  scfcv_init
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! INPUTS
!!  scfcv=structure of scfcv
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_gstate
!!
!! CHILDREN
!!      scfcv_core
!!
!! SOURCE

subroutine scfcv_init(this,atindx,atindx1,cg,cprj,cpus,&
&  dmatpawu,dtefield,dtfil,dtorbmag,dtpawuj,dtset,ecore,eigen,hdr,&
&  indsym,initialized,irrzon,kg,mcg,mcprj,mpi_enreg,my_natom,nattyp,ndtpawuj,&
&  nfftf,npwarr,occ,pawang,pawfgr,pawrad,pawrhoij,&
&  pawtab,phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,&
&  resid,results_gs,scf_history,fatvshift,&
&  symrec,taug,taur,wvl,ylm,ylmgr,paw_dmft,wffnew,wffnow)


!Arguments ------------------------------------
!scalars
 type(scfcv_t), intent(inout) :: this
 integer,intent(in),target :: mcg,mcprj,my_natom,ndtpawuj,pwind_alloc
 integer,intent(in),target :: initialized,nfftf
 real(dp),intent(in),target :: cpus,ecore
 real(dp),intent(in),target :: fatvshift
 type(MPI_type),intent(in),target :: mpi_enreg
 type(datafiles_type),intent(in),target :: dtfil
 type(dataset_type),intent(in),target :: dtset
 type(efield_type),intent(in),target :: dtefield
 type(orbmag_type),intent(in),target :: dtorbmag
! type(electronpositron_type),pointer :: electronpositron
 type(hdr_type),intent(in),target :: hdr
 type(pawang_type),intent(in),target :: pawang
 type(pawfgr_type),intent(in),target :: pawfgr
 type(pseudopotential_type),intent(in),target :: psps
 type(recursion_type),intent(in),target :: rec_set
 type(results_gs_type),intent(in),target :: results_gs
 type(scf_history_type),intent(in),target :: scf_history
! type(wffile_type),intent(in),target :: wffnew,wffnow
 type(wvl_data),intent(in),target :: wvl
!arrays
 integer,intent(in),target :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(in),target :: indsym(4,dtset%nsym,dtset%natom)
!no_abirules
 integer, intent(in),target :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 integer, intent(in),target :: kg(3,dtset%mpw*dtset%mkmem)
 integer, intent(in),target :: nattyp(psps%ntypat),npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
 integer, intent(in),target :: symrec(3,3,dtset%nsym)
 real(dp), intent(in),target :: cg(2,mcg),dmatpawu(:,:,:,:)
 real(dp), intent(in),target :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in),target :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in),target :: phnons(2,dtset%nfft**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp), intent(in),target :: pwnsfac(2,pwind_alloc)
! real(dp), intent(in),target :: rprimd(3,3)
! real(dp), pointer :: rhog(:,:),rhor(:,:)
 real(dp), pointer :: taug(:,:),taur(:,:)
 real(dp), intent(in),target :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
! real(dp), intent(in),target :: xred(3,dtset%natom),xred_old(3,dtset%natom)
 real(dp), intent(in),target :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in),target :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(macro_uj_type),intent(in),target :: dtpawuj(0:ndtpawuj)
 type(pawrhoij_type), intent(in),target :: pawrhoij(my_natom*psps%usepaw)
 type(pawrad_type), intent(in),target :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type), intent(in),target :: pawtab(psps%ntypat*psps%usepaw)
 !type(dataset_type),intent(in),target :: dtset
! type(electronpositron_type),intent(in),target :: electronpositron
 type(paw_dmft_type), intent(in),target :: paw_dmft
 type(wffile_type),intent(in),target :: wffnew,wffnow
 type(pawcprj_type), allocatable,intent(in),target :: cprj(:,:)
!Local variables -------------------------
!scalars
 logical :: DEBUG=.FALSE.
! *************************************************************************

 DBG_ENTER("COLL")

 if (DEBUG) then
   write(std_out,*) 'INTENT(IN) ARGUMENTS ON SCFCV'
 !  write(std_out,*) 'atindx=',ab_scfcv_in%atindx
 !  write(std_out,*) 'atindx1=',ab_scfcv_in%atindx1
 !  write(std_out,*) 'cpus=',ab_scfcv_in%cpus
 !  write(std_out,*) 'ecore=',ab_scfcv_in%ecore
 !  write(std_out,*) 'fatvshift=',ab_scfcv_in%fatvshift
 !  write(std_out,*) 'indsym=',ab_scfcv_in%indsym
 !  write(std_out,*) 'kg=',ab_scfcv_in%kg
 !  write(std_out,*) 'my_natom=',ab_scfcv_in%my_natom
 !  write(std_out,*) 'nattyp=',ab_scfcv_in%nattyp
 !  write(std_out,*) 'ndtpawuj=',ab_scfcv_in%ndtpawuj
 !  write(std_out,*) 'npwarr=',ab_scfcv_in%npwarr
 !  write(std_out,*) 'phnons=',ab_scfcv_in%phnons
 !  write(std_out,*) 'pwind=',ab_scfcv_in%pwind
 !  write(std_out,*) 'pwind_alloc=',ab_scfcv_in%pwind_alloc
 !  write(std_out,*) 'pwnsfac=',ab_scfcv_in%pwnsfac
 !  write(std_out,*) 'ylm=',ab_scfcv_in%ylm
 !  write(std_out,*) 'ylmgr=',ab_scfcv_in%ylmgr
!!  write(std_out,*) 'pawang=',ab_scfcv_in%pawang
!!  write(std_out,*) 'pawrad=',ab_scfcv_in%pawrad
!!  write(std_out,*) 'pawtab=',ab_scfcv_in%pawtab
!!  write(std_out,*) 'psps=',ab_scfcv_in%psps
 end if

 this%atindx=>atindx
 this%atindx1=>atindx1
 this%cpus=>cpus
 this%ecore=>ecore
 this%fatvshift=>fatvshift
 this%indsym=>indsym
 this%kg=>kg
 this%mcg=>mcg
 this%mcprj=>mcprj
 this%my_natom=>my_natom
 this%nattyp=>nattyp
 this%ndtpawuj=>ndtpawuj
 this%npwarr=>npwarr
 this%pawang=>pawang
 this%pawrad=>pawrad
 this%pawtab=>pawtab
 this%phnons=>phnons
 this%psps=>psps
 this%pwind=>pwind
 this%pwind_alloc=>pwind_alloc
 this%pwnsfac=>pwnsfac
 this%ylm=>ylm
 this%ylmgr=>ylmgr

 this%cg=>cg
 this%cprj=>cprj
 this%dmatpawu=>dmatpawu
 this%dtefield=>dtefield
 this%dtorbmag=>dtorbmag
 this%dtfil=>dtfil
 this%dtpawuj=>dtpawuj
 this%eigen=>eigen
 this%hdr=>hdr
 this%initialized=>initialized
 this%irrzon=>irrzon
 this%mpi_enreg=>mpi_enreg
 this%nfftf=>nfftf
 this%occ=>occ
 this%pawfgr=>pawfgr
 this%pawrhoij=>pawrhoij
 this%pawtab=>pawtab
 this%rec_set=>rec_set
 this%resid=>resid
 this%results_gs=>results_gs
 this%scf_history=>scf_history
 this%symrec=>symrec
 this%taug=>taug
 this%taur=>taur
 this%wvl=>wvl

 this%dtset=>dtset
 !this%electronpositron=>electronpositron
 this%paw_dmft=>paw_dmft
 !this%rhog=>rhog
 !this%rhor=>rhor
 !this%rprimd=>rprimd
 this%wffnew=>wffnew
 this%wffnow=>wffnow
 !this%xred=>xred
 !this%xred_old=>xred_old


 !!!!!!!!! INITIALIZE or REINITIALIZE parallelization here !!
 ! TODO at next step
 !if ( this%dtset%usedmft /= 0 ) then
 !  call data4entropyDMFT_init(this%paw_dmft%forentropyDMFT,&
 !                            this%dtset%natom,&
 !                            this%dtset%typat,&
 !                            this%dtset%lpawu,&
 !                            this%dtset%dmft_t2g==1, &
 !                            this%dtset%upawu,&   !!! Should use this%pawtab%upawu
 !                            this%dtset%jpawu)    !!! Should use this%pawtab%jpawu
 !end if

 !call entropyDMFT_init(this%entropyDMFT,this%dtset,this%pawtab,this%mpi_enreg%comm_cell,this%dtfil%filnam_ds(3),this%dtfil%filnam_ds(4)) ! Do something only if DMFT and dmft_entropy = 1

 DBG_EXIT("COLL")

end subroutine scfcv_init
!!***


!!****f* ABINIT/m_scfcv/scfcv_destroy
!! NAME
!!  scfcv_destroy
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! INPUTS
!!  scfcv=structure of scfcv
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_gstate,m_mover_effpot
!!
!! CHILDREN
!!      scfcv_core
!!
!! SOURCE

subroutine scfcv_destroy(this)

!Arguments ------------------------------------
type(scfcv_t), intent(inout) :: this

!Local variables-------------------------------

! *************************************************************************

 DBG_ENTER("COLL")

 !scalars
 this%mcg => null()
 this%mcprj => null()
 this%my_natom => null()
 this%ndtpawuj => null()
 this%pwind_alloc => null()
 this%initialized => null()
 this%nfftf => null()
 this%cpus => null()
 this%ecore  => null()
 this%fatvshift => null()
 this%pawang => null()
 this%psps => null()
 this%mpi_enreg => null()
 this%dtfil => null()
 this%dtset => null()
 this%dtefield => null()
 this%dtorbmag => null()
 this%electronpositron => null()
 this%hdr => null()
 this%pawfgr => null()
 this%rec_set => null()
 this%results_gs => null()
 this%scf_history => null()
 this%wffnew => null()
 this%wffnow => null()
 this%wvl => null()
 this%paw_dmft => null()

 !arrays
 this%atindx => null()
 this%atindx1 => null()
 this%irrzon => null()
 this%symrec => null()
 this%indsym => null()
 !no_abirules
 this%kg => null()
 this%nattyp => null()
 this%npwarr => null()
 this%pwind => null()
 this%phnons => null()
 this%pwnsfac => null()
 this%ylm => null()
 this%ylmgr => null()
 this%cg => null()
 this%cprj => null()
 this%dmatpawu => null()
 this%eigen => null()
 this%occ => null()
 !this%rprimd
 !this%rhog => null()
 !this%rhor => null()
 this%taug => null()
 this%taur => null()
 this%resid => null()
 this%pawrad => null()
 this%pawtab => null()
 this%dtpawuj => null()
 this%pawrhoij => null()

 ! This call should be done inside destroy_sc_dmft
 !if ( this%dtset%usedmft /= 0 ) then
 !  call data4entropyDMFT_destroy(this%paw_dmft%forentropyDMFT)
 !end if
 !call entropyDMFT_destroy(this%entropyDMFT)

 DBG_EXIT("COLL")

end subroutine scfcv_destroy
!!***

!!****f* ABINIT/m_scfcv/scfcv_run
!! NAME
!!  scfcv_run
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! INPUTS
!!  itime=Relaxation iteration.
!!  scfcv=structure of scfcv
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gstate,m_mover
!!
!! CHILDREN
!!      scfcv_core
!!
!! SOURCE

subroutine scfcv_run(this, itime, electronpositron, rhog, rhor, rprimd, xred, xred_old, conv_retcode)

!Arguments ------------------------------------
 type(scfcv_t), intent(inout) :: this
 integer,intent(in) :: itime
 type(electronpositron_type),pointer:: electronpositron
 real(dp), intent(inout) :: rprimd(3,3)
 real(dp), intent(inout) :: xred(3,this%dtset%natom)
 real(dp), intent(inout) :: xred_old(3,this%dtset%natom)
 real(dp), pointer, intent(inout) :: rhog(:,:)
 real(dp), pointer, intent(inout) :: rhor(:,:)
 integer ,intent(out) :: conv_retcode

!Local variables-------------------------------

! *************************************************************************

 DBG_ENTER("COLL")

!!!  Should be changed if special parallelization.

 ! Moved inside mover.F90 before this call
 !call scfcv_reformatWFK(this,rhog, rhor, rprimd, xred, xred_old)

 ! First initialize the datatype to gather information

 !debug purpose
 !this%electronpositron => electronpositron
 if ( this%dtset%dmft_entropy == 0 ) then
   call scfcv_scfcv(this, itime, electronpositron, rhog, rhor, rprimd, xred, xred_old, conv_retcode)
 elseif ( this%dtset%dmft_entropy >=1 ) then
   call scfcv_runWEntropyDMFT(this, itime, electronpositron,rhog,rhor,rprimd,xred,xred_old,conv_retcode)
 end if

 DBG_EXIT("COLL")

end subroutine scfcv_run
!!***


!!!!****f* ABINIT/m_scfcv/scfcv_reformatWFK
!!!! NAME
!!!!  scfcv_reformatWFK
!!!!
!!!! FUNCTION
!!!!  FIXME: add description.
!!!!
!!!! INPUTS
!!!!  scfcv=structure of scfcv
!!!!  argin(sizein)=description
!!!!
!!!! OUTPUT
!!!!  argout(sizeout)=description
!!!!
!!!! SIDE EFFECTS
!!!!
!!!! NOTES
!!!!
!!!! PARENTS
!!!!
!!!! CHILDREN
!!!!
!!!! SOURCE
!!
!!subroutine scfcv_reformatWFK(this,rhog, rhor, rprimd, xred, xred_old)
!!
!!
!! type(scfcv_t), intent(inout) :: this
!! real(dp), intent(inout) :: xred(3,this%dtset%natom)
!! real(dp), intent(inout) :: xred_old(3,this%dtset%natom)
!! real(dp), intent(inout) :: rprimd(3,3)
!! real(dp), pointer, intent(inout) :: rhog(:,:)
!! real(dp), pointer, intent(inout) :: rhor(:,:)
!!
!!!WVL - reformat the wavefunctions in the case of xred != xred_old
!! if (this%dtset%usewvl == 1 .and. maxval(xred_old - xred) > zero) then
!!!  WVL - Before running scfcv, on non-first geometry step iterations,
!!!  we need to reformat the wavefunctions, taking into acount the new
!!!  coordinates.
!!!  We prepare to change rhog (to be removed) and rhor.
!!   ABI_DEALLOCATE(rhog)
!!   ABI_DEALLOCATE(rhor)
!!
!!   call wvl_wfsinp_reformat(this%dtset, this%mpi_enreg,&
!!&   this%psps, rprimd, this%wvl, xred, xred_old)
!!   this%nfftf = this%dtset%nfft
!!
!!   ABI_ALLOCATE(rhog,(2, this%dtset%nfft))
!!   ABI_ALLOCATE(rhor,(2, this%dtset%nfft))
!!   call wvl_mkrho(this%dtset, this%irrzon, this%mpi_enreg,&
!!&   this%phnons, rhor,this%wvl%wfs,this%wvl%den)
!! end if
!!
!!end subroutine scfcv_reformatWFK
!!!!***

!!****f* ABINIT/m_scfcv/scfcv_runWEntropyDMFT
!! NAME
!!  scfcv_runWEntropyDMFT
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! INPUTS
!!  scfcv=structure of scfcv
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      m_scfcv
!!
!! CHILDREN
!!      scfcv_core
!!
!! SOURCE

subroutine scfcv_runWEntropyDMFT(this,itime,electronpositron,rhog,rhor,rprimd,xred,xred_old,conv_retcode)


!Arguments ------------------------------------
 type(scfcv_t), intent(inout) :: this
 integer,intent(in) :: itime
 type(electronpositron_type),pointer :: electronpositron
 real(dp), intent(inout) :: rprimd(3,3)
 real(dp), intent(inout) :: xred(3,this%dtset%natom)
 real(dp), intent(inout) :: xred_old(3,this%dtset%natom)
 real(dp), pointer, intent(inout) :: rhog(:,:)
 real(dp), pointer, intent(inout) :: rhor(:,:)
 integer , intent(out)   :: conv_retcode

!Local variables-------------------------------

! *************************************************************************

 DBG_ENTER("COLL")

 !if ( this%dtset%usedmft /= 0 ) then
 !  call data4entropyDMFT_init(this%paw_dmft%forentropyDMFT,&
 !                            this%dtset%natom,&
 !                            this%dtset%typat,&
 !                            this%dtset%lpawu,&
 !                            this%dtset%dmft_t2g==1, &
 !                            this%dtset%upawu,&     !!! Should use this%pawtab%upawu
 !                            this%dtset%jpawu)      !!! Should use this%pawtab%jpawu
 !end if

 call entropyDMFT_init(this%entropyDMFT,this%dtset,this%pawtab,this%mpi_enreg%comm_cell,&
&       this%dtfil%filnam_ds(3),this%dtfil%filnam_ds(4)) ! Do something only if DMFT and dmft_entropy = 1

 ! Start loop over all integration points (lambda)
 ! TODO WORK ON PARALLELISATION HERE
 do while (entropyDMFT_nextLambda(this%entropyDMFT,this%dtset,this%pawtab,this%pawang,this%pawrad))

 !-----------------------------------------------------
   call scfcv_scfcv(this, itime, electronpositron,rhog,rhor,rprimd,xred,xred_old,conv_retcode)
 !-----------------------------------------------------
   call entropyDMFT_addIntegrand(this%entropyDMFT,this%dtset, this%results_gs%energies,this%paw_dmft%forentropyDMFT)

 end do !!! End loop for entropy DMFT

 ! GATHER DATA HERE OR INSIDE THE NEXT CALL ?
 call entropyDMFT_computeEntropy(this%entropyDMFT,this%results_gs%energies%entropy)
 !-----------------------------------------------------
 ! This call should be done inside destroy_sc_dmft
 !if ( this%dtset%usedmft /= 0 ) then
 !  call data4entropyDMFT_destroy(this%paw_dmft%forentropyDMFT)
 !end if
 call entropyDMFT_destroy(this%entropyDMFT)

 DBG_EXIT("COLL")

end subroutine scfcv_runWEntropyDMFT
!!***

!!****f* ABINIT/m_scfcv/scfcv_scfcv
!! NAME
!!  scfcv_scfcv
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! INPUTS
!!  scfcv=structure of scfcv
!!  itime: Relaxation step
!!
!! OUTPUT
!!
!! NOTES
!!  Wrapper to scfcv to avoid circular dependencies ...
!!
!! PARENTS
!!      m_scfcv
!!
!! CHILDREN
!!      scfcv_core
!!
!! SOURCE

subroutine scfcv_scfcv(this, itime, electronpositron, rhog, rhor, rprimd, xred, xred_old, conv_retcode)

 type(scfcv_t), intent(inout) :: this
 integer,intent(in) :: itime
 type(electronpositron_type),pointer :: electronpositron
 real(dp), intent(inout) :: rprimd(3,3)
 real(dp), intent(inout) :: xred(3,this%dtset%natom)
 real(dp), intent(inout) :: xred_old(3,this%dtset%natom)
 real(dp), pointer, intent(inout) :: rhog(:,:)
 real(dp), pointer, intent(inout) :: rhor(:,:)
 integer , intent(out)   :: conv_retcode

   call scfcv_core(itime, this%atindx,this%atindx1,this%cg,this%cprj,this%cpus,this%dmatpawu,this%dtefield,this%dtfil,&
    this%dtorbmag,this%dtpawuj,&
    this%dtset,this%ecore,this%eigen,electronpositron,this%fatvshift,this%hdr,this%indsym,&
    this%initialized,this%irrzon,this%kg,this%mcg,this%mcprj,this%mpi_enreg,this%my_natom,this%nattyp,this%ndtpawuj,&
    this%nfftf,this%npwarr,&
    this%occ,this%paw_dmft,this%pawang,this%pawfgr,this%pawrad,this%pawrhoij,this%pawtab,this%phnons,this%psps,this%pwind,&
    this%pwind_alloc,this%pwnsfac,this%rec_set,this%resid,this%results_gs,rhog,rhor,rprimd,&
    this%scf_history,this%symrec,this%taug,this%taur,this%wffnew,this%wvl,xred,xred_old,this%ylm,this%ylmgr,&
    conv_retcode)

end subroutine scfcv_scfcv

end module m_scfcv
!!***
