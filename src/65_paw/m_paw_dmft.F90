!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_dmft
!! NAME
!!  m_paw_dmft
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2006-2019 ABINIT group (BAmadon)
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

MODULE m_paw_dmft

 use defs_basis
 use m_CtqmcInterface
 use m_errors
 use m_abicore
 use m_xmpi
 use m_data4entropyDMFT
 use m_dtset

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes, only : MPI_type
 use m_io_tools,  only : open_file
 use m_pawtab,    only : pawtab_type

 implicit none

 private

 public :: init_dmft
 public :: init_sc_dmft
 public :: construct_nwli_dmft
 public :: destroy_dmft
 public :: destroy_sc_dmft
 public :: print_dmft
 public :: print_sc_dmft
 public :: saveocc_dmft
 public :: readocc_dmft

 private :: init_sc_dmft_paralkgb, destroy_sc_dmft_paralkgb
!!***

!----------------------------------------------------------------------

!!****t* m_paw_dmft/paw_dmft_type
!! NAME
!!  paw_dmft_type
!!
!! FUNCTION
!!  This structured datatype contains the necessary data for the link
!!  between dmft and paw.
!!  occnd(non-diagonal band occupations for self-consistency), band_in
!!  (say which band are taken into account in the calculation), and the
!   dimensions of these arrays.
!!
!! SOURCE

 type, public :: paw_dmft_type

  integer :: dmft_dc
  ! Type of double counting used in DMFT

  integer :: dmft_iter
  ! Nb of iterations for dmft

  integer :: dmft_kspectralfunc
  ! =0 Default
  ! =1 Activate calculation of k-resolved spectral function

  integer :: dmft_solv
  ! choice of solver for DMFT

  integer :: dmftcheck
  ! Check various part of the implementation

  integer :: dmft_log_freq
  ! = 0: do not use log frequencies
  ! = 1: use log frequencies

!  integer :: dmft_mag
!  ! 0 if non magnetic calculation, 1 if magnetic calculation

  integer :: dmft_nwlo
  ! dmft frequencies

  integer :: dmft_nwr
  ! dmft frequencies

  integer :: dmft_nwli
  ! dmft frequencies

  integer :: dmftqmc_l
  ! qmc related input

  integer :: dmftqmc_seed
  ! qmc seed

  integer :: dmftqmc_therm
  ! qmc thermalization

  integer :: dmftctqmc_basis
  ! CTQMC: Basis to perform the CTQMC calculation
  ! for historical reasons and tests
  ! 0 : Slm basis, 1 : diagonalise local Hamiltonian, 2: diago density matrix
  !
  integer :: dmft_blockdiag
  !  Block diagonalize Hamiltonian in the local basis

  integer :: dmftctqmc_check
  ! CTQMC: perform a check on the impurity and/or bath operator
  ! only for debug
  ! 0 : nothing, 1 : impurity, 2 : bath, 3 : both

  integer :: dmftctqmc_correl
  ! CTQMC: Gives analysis for CTQMC
  ! 0 : nothing, 1 : activated Correlations.dat

  integer :: dmftctqmc_gmove
  ! CTQMC: add global move every dmftctqmc_gmove sweeps
  ! >= 0 ; warning inside CT-QMC with warning
  ! == 0 ; no global moves

  integer :: dmftctqmc_grnns
  ! CTQMC: compute green function noise for each imaginary time
  ! 0 : nothing, 1 : activated

  integer :: dmftctqmc_meas
  ! CTQMC: every each dmftctqmc_meas step energy is measured
  ! Speed up caltion about 10% for dmftctqmc_meas = 2

  integer :: dmftctqmc_mrka
  ! CTQMC: Write a temporary file Spectra_RANK.dat with the sweep evolution of
  ! the number of electron for each flavor
  ! The measurement is done every dmftctqmc_meas*dmftctqmc_mrka sweep
  ! e.g. : meas=2 mrka=10 -> every 20 sweeps sum_i c+(ti)c(t'i) is mesured

  integer :: dmftctqmc_mov
  ! CTQMC: Gives movie for CTQMC
  ! 0 : nothing, 1 : 1 file Movie_RANK.tex for each cpu

  integer :: dmftctqmc_order
  ! CTQMC: Gives order in perturbation for CTQMC solver
  ! 0 : nothing, >=1 max order evaluated in Perturbation.dat

  integer :: dmftctqmc_triqs_nleg
  ! CTQMC of TRIQS: Nb of Legendre polynomial used to compute the
  ! Green's function (Phys. Rev. B 84, 075145) [[cite:Boehnke2011]]. Default is 30.

  ! 0 : nothing, >=1 max order evaluated in Perturbation.dat

  real(dp) :: dmftqmc_n
  ! qmc number of sweeps

  integer :: dmftqmc_x2my2d
  ! for doing qmc with x2my2d only (for testing purposes)

  integer :: dmftqmc_t2g
  ! for doing qmc with t2g only (for testing purposes)

  integer :: dmftbandi
  ! Number of bands

  integer :: dmftbandf
  ! Number of bands

  integer :: dmft_read_occnd
  ! Number of bands

  integer :: dmft_rslf
  ! Number of bands

  integer :: dmft_prgn
  ! Precise the way of printing the green function.
  !  =1   print green
  !  =2   print self

  integer :: idmftloop
  ! current iteration in the dmft loop

  integer :: maxlpawu         ! Number of correlated atoms


  integer :: mband
  ! Number of bands

  integer :: mbandc
  ! Total number of bands in the Kohn-Sham Basis for PAW+DMFT

  integer :: natom
  ! Number of atom

  integer :: natpawu         ! Number of correlated atoms

  integer :: nkpt
  ! Number of k-point in the IBZ.

  integer :: nspden

  integer :: nspinor

  integer :: nsppol

  integer :: prtdos

  integer :: prtvol

  integer  :: lpsichiortho

  integer  :: use_fixed_self

  real(dp) :: edmft

  real(dp) :: dmft_charge_prec
  ! Precision on charge required for determination of fermi level (fermi_green) with newton method

  real(dp) :: dmft_fermi_prec
  ! Required precision on Fermi level (fermi_green) during the DMFT SCF cycle, (=> ifermie_cv)
  ! used also for self (new_self)  (=> iself_cv).

  real(dp) :: dmft_mxsf
  ! Mixing coefficient for Self-Energy during the SCF DMFT cycle.

  real(dp) :: dmft_tolfreq
  ! Required precision on local correlated density matrix  (depends on
  ! frequency mesh), used in m_dmft/dmft_solve

  real(dp) :: dmft_lcpr
  ! Required precision on local correlated charge  in order to stop SCF
  ! DMFT cycle (integrate_green) => ichargeloc_cv

  real(dp) :: fermie


  real(dp) :: fermie_lda

  real(dp) :: nelectval

  character(len=fnlen) :: filapp

  character(len=fnlen) :: filnamei

  real(dp) :: temp

  integer, allocatable :: lpawu(:)

  integer, allocatable :: include_bands(:)
  ! for each bands included in the calculation (1..mbandc), include_bands
  ! gives the index in the full band index  (1...mband)

  integer, allocatable :: exclude_bands(:)
  ! gives the bands than are not in the DMFT calculations.

  real(dp), allocatable :: occnd(:,:,:,:,:)
  ! non diagonal band-occupation for each k-point, polarisation.

!  real(dp), allocatable :: phi0phiiint(:)
!  ! non diagonal band-occupation for each k-point, polarisation.

  logical, allocatable :: band_in(:)
  ! true for each band included in the calculation.

  integer, allocatable :: bandc_proc(:)
  ! proc index (on comm_band) for each correlated band in DMFT

  logical, allocatable :: use_bandc(:)
  ! true for each proc wich has at least one band involved in DMFT non diagonal
  ! occupations on band parallelism

  integer :: use_dmft
  ! 1 if non diagonal occupations are used, else 0

  integer :: use_sc_dmft
  ! 1 if calculations have to be carried out self-consistently in the
  ! electronic density.

  complex(dpc), allocatable :: psichi(:,:,:,:,:,:)

  real(dp), allocatable :: eigen_lda(:,:,:)

  real(dp), pointer :: wtk(:) => null()
  real(dp), pointer :: fixed_self(:,:,:,:) => null()
  real(dp), pointer :: omega_lo(:) => null()
  real(dp), pointer :: omega_r(:) => null()
  real(dp), pointer :: wgt_wlo(:) => null()

  type(CtqmcInterface), allocatable :: hybrid(:)
  type(data4entropyDMFT_t) :: forentropyDMFT

  ! MPI realated variables
  integer :: myproc
  integer :: nproc
  integer :: spacecomm

 end type paw_dmft_type
!!***

!----------------------------------------------------------------------

CONTAINS  !========================================================================================
!!***

!!****f* m_paw_dmft/init_sc_dmft
!! NAME
!! init_sc_dmft
!!
!! FUNCTION
!!  Allocate variables used in type paw_dmft_type.
!!
!! INPUTS
!! dmftbandi = lower bound for band states included in DMFT calculation
!! dmftbandf = upper bound for band states included in DMFT calculation
!! mband     = max number of bands
!! nband     = number of bands for each k-point
!! nkpt      = number of k-points
!! nsppol    = number of spin polarisation
!! occ       =  occupations
!! usedmft  = ==1 if dmft is activated
!! use_sc_dmft = for self-consistency in dmft
!!
!! OUTPUTS
!! paw_dmft  = structure of data for dmft
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine init_sc_dmft(bandkss,dmftbandi,dmftbandf,dmft_read_occnd,mband,nband,nkpt,nspden,&
&nspinor,nsppol,occ,usedmft,paw_dmft,use_sc_dmft,dmft_solv,mpi_enreg)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: bandkss,dmft_read_occnd,dmftbandi,dmftbandf,mband,nkpt,nspden,&
& nspinor,nsppol,usedmft,use_sc_dmft,dmft_solv
!type
 type(paw_dmft_type),intent(out) :: paw_dmft
 type(MPI_type), intent(in) :: mpi_enreg
! arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol)
!Local variables ------------------------------------
 integer :: iband,icb,ikpt,isppol,nband_k,bdtot_index
 integer :: myproc,nproc,spacecomm,use_dmft
! integer :: ie,nb_procs
 character(len=500) :: message

!************************************************************************
 use_dmft=abs(usedmft)


! Check processors for DMFT
! Initialise spaceComm, myproc, and nproc

 spacecomm=mpi_enreg%comm_cell
 myproc=mpi_enreg%me_cell
 nproc=mpi_enreg%nproc_cell
 spacecomm=mpi_enreg%comm_world
 myproc=mpi_enreg%me
 nproc=mpi_enreg%nproc
 !print *, " spacecomm,myproc,nproc",spacecomm,myproc,nproc
 paw_dmft%spacecomm=spacecomm
 paw_dmft%myproc=myproc
 paw_dmft%nproc=nproc

 ! Do not comment these lines: it guarantees the parallelism in DMFT/HI or QMC will work.
 if ((use_dmft/=0).and.(xmpi_comm_size(xmpi_world) /= xmpi_comm_size(mpi_enreg%comm_world))) then
   MSG_ERROR("Someone changed the k-point parallelism again")
 end if

 if(use_dmft/=0) then
   write(message,'(2a,i3)') ch10,'-       ( number of procs used in dmft ) =', nproc
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   write(std_out_default,'(2a,i3)') ch10,'       ( current proc is        ) =', myproc
  ! write(ab_out_default,'(2a,i3)') ch10,'       ( current proc is        ) =', myproc
   if(myproc==nproc-1) then
     write(std_out_default,'(2a,i3)') ch10,'      ( last proc            ) =', myproc
  !   write(ab_out_default,'(2a,i3)') ch10,'       ( last proc            ) =', myproc
   endif
 endif
!#ifdef HAVE_MPI
! call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,ie)
! write(6,*) "nprocs,nb_procs",nproc,nb_procs
! if(nb_procs/=nproc)  then
!   message = ' Number of procs used in DMFT is erroneously computed '
!   MSG_ERROR(message)
! endif
!#endif

 paw_dmft%mband       = mband
 paw_dmft%dmftbandf   = dmftbandf
 paw_dmft%dmftbandi   = dmftbandi
 paw_dmft%nkpt        = nkpt

!  Spin related variables and check
 paw_dmft%nsppol      = nsppol
 paw_dmft%nspinor     = nspinor
 paw_dmft%nspden      = nspden
 if(nspinor==2.and.nspden==1.and.use_dmft/=0) then
   message = ' nspinor==2 and nspden =1 and usedmft=1 is not implemented'
   MSG_ERROR(message)
 endif

! if(nspinor==1.and.nspden==1.and.use_dmft/=0) then
!   message = ' nspinor==1 and nspden =1 and usedmft=1 is not implemented'
!   MSG_ERROR(message)
! endif

 paw_dmft%use_dmft    = use_dmft
 if (bandkss/=0) then
   paw_dmft%use_sc_dmft = 0
 else
   paw_dmft%use_sc_dmft = use_sc_dmft
 endif
 paw_dmft%dmft_read_occnd = dmft_read_occnd
 paw_dmft%idmftloop=0
 paw_dmft%mbandc  = 0
 ABI_ALLOCATE(paw_dmft%occnd,(2,mband,mband,nkpt,nsppol*use_dmft))
 ABI_ALLOCATE(paw_dmft%band_in,(mband*use_dmft))
 ABI_ALLOCATE(paw_dmft%include_bands,((dmftbandf-dmftbandi+1)*use_dmft))
 ABI_ALLOCATE(paw_dmft%exclude_bands,(mband*use_dmft))
! allocate(paw_dmft%ph0phiiint()
 paw_dmft%band_in(:)=.false.
 paw_dmft%occnd=zero
 icb=0
 if(use_dmft==1) then
  do iband=1, mband
   if(iband>=paw_dmft%dmftbandi.and.iband<=paw_dmft%dmftbandf) then
    paw_dmft%band_in(iband)=.true.
    paw_dmft%mbandc = paw_dmft%mbandc+1
    paw_dmft%include_bands(paw_dmft%mbandc) = iband
   else
    icb=icb+1
    paw_dmft%exclude_bands(icb)=iband
   endif
  enddo
  bdtot_index=1
  do isppol=1,nsppol
   do ikpt=1,nkpt
    nband_k=nband(ikpt+(isppol-1)*nkpt)
    do iband=1,nband_k
     paw_dmft%occnd(1,iband,iband,ikpt,isppol)=occ(bdtot_index)
     bdtot_index=bdtot_index+1
    end do
   end do
  end do
 else
  paw_dmft%mbandc = 0
 endif

 if(paw_dmft%use_sc_dmft /= 0 .and. mpi_enreg%paral_kgb/=0) then
   call init_sc_dmft_paralkgb(paw_dmft, mpi_enreg)
 end if

 if(paw_dmft%use_dmft > 0 .and. paw_dmft%mbandc /= dmftbandf-dmftbandi+1) then
  write(message, '(3a)' )&
&  ' WARNING init_sc_dmft',ch10,&
&  '  number of bands in dmft is not correctly computed ',ch10, &
&  '  Action : check the code'
  MSG_WARNING(message)
 endif
 if(use_dmft>=1) then
   write(message, '(7a)' ) ch10,ch10," ******************************************", &
&   ch10," DFT+DMFT Method is used", &
&   ch10," ******************************************"
   call wrtout(std_out,  message,'COLL')
   call wrtout(ab_out,  message,'COLL')
 endif

 if(use_dmft>=1) then
   if(dmft_solv==0) then
     write(message, '(a,a)') ch10,' DMFT check: no solver and U=J=0'
   else if(dmft_solv==1) then
     write(message, '(a,a)') ch10,' DMFT check: static solver'
   else if(dmft_solv==-1) then
     write(message, '(a,a)') ch10,' DMFT check: static solver without renormalization of projectors: should recover LDA+U'
   else if(dmft_solv==2) then
     write(message, '(a,a)') ch10,' DMFT uses the Hubbard one solver'
   else if(dmft_solv==4) then
     write(message, '(a,a)') ch10,' DMFT uses the Hirsch Fye solver'
   else if(dmft_solv==5) then
     write(message, '(a,a)') ch10,' DMFT uses the Continuous Time Quantum Monte Carlo solver of ABINIT'
   else if(dmft_solv==6) then
     write(message, '(a,a)') ch10,' DMFT uses the Continuous Time Quantum Monte Carlo solver of TRIQS&
     & (with density density interactions)'
   else if(dmft_solv==7) then
     write(message, '(a,a)') ch10,' DMFT uses the Continuous Time Quantum Monte Carlo solver of TRIQS&
     & (with rotationaly invariant interaction)'
  endif
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')
 endif

end subroutine init_sc_dmft
!!***

!!****f* m_paw_dmft/init_dmft
!! NAME
!! init_dmft
!!
!! FUNCTION
!!  Allocate variables and setup lda hamiltonian and related data
!!  (init_sc_dmft has to been called before)
!!
!! INPUTS
!!  dmatpawu   = fixed occupation matrix of correlated orbitals
!!  eigen      = LDA eigenvalues
!!  fermie_lda = LDA Fermi level
!!  psichi     = <chi|Psi> projection of KS states over atomic !wavefunction
!!  nkpt       = number of k-points
!!  nsppol     = number of spin polarisation
!!  nspinor    = number of spinorial component
!!
!!
!! PARENTS
!!      outscfcv,vtorho
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE
!!
!! NOTE
!! The part of the code which deals
!! with the use of logarithmic frequencies
!! is a modification of the GNU GPL
!! code available on http://dmft.rutgers.edu/ and
!! described in the  RMP paper written by
!! G.Kotliar,  S.Y.Savrasov, K.Haule, V.S.Oudovenko, O.Parcollet, C.A.Marianetti.
!!

subroutine init_dmft(dmatpawu, dtset, fermie_lda, fnametmp_app, fnamei, nspinor, paw_dmft, pawtab, psps, typat)

 use m_splines
 !use m_CtqmcInterface

!Arguments ------------------------------------
!scalars
 integer, intent(in)  :: nspinor
 real(dp), intent(in) :: fermie_lda
!type
 type(pseudopotential_type), intent(in) :: psps
 type(dataset_type),target,intent(in) :: dtset
 type(pawtab_type),intent(in)  :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_dmft_type),intent(inout) :: paw_dmft
 character(len=fnlen), intent(in) :: fnametmp_app
 character(len=fnlen), intent(in) :: fnamei
!arrays
 integer,intent(in) :: typat(dtset%natom)
 real(dp),intent(in),target :: dmatpawu(:,:,:,:)
!Local variables ------------------------------------
 integer :: grid_unt,ikpt,isym,itypat,nsppol,mbandc,maxlpawu, iatom, ifreq, ioerr
 integer :: nflavor,ngrid,iexist2
 character(len=fnlen) :: tmpfil
 real(dp) :: sumwtk
 character(len=500) :: message
 real(dp) :: unit_e,step
 logical :: lexist
! *********************************************************************

 nsppol = dtset%nsppol
 if(dtset%ucrpa==0) then
 write(message,'(6a)') ch10,' ====================================', &
&                      ch10,' =====  Start of DMFT calculation', &
&                      ch10,' ===================================='
 else if(dtset%ucrpa>0) then
 write(message,'(6a)') ch10,' ============================================================', &
&                      ch10,' =====  Initialize construction of Wannier in DMFT routines',&
&                      ch10,' ============================================================'
 endif
 call wrtout(std_out,message,'COLL')

 unit_e=2_dp
!=======================
!==  Check sym
!=======================
 do isym=1,dtset%nsym
   if(dtset%symafm(isym)<0) then
     message = 'symafm negative is not implemented in DMFT '
     MSG_ERROR(message)
   endif
 enddo

! paw_dmft%dmft_mag=0
! do iatom=1,dtset%natom
!   do  ii=1,3
!     if ( dtset(ii,iatom) > 0.001 ) paw_dmft%dmft_mag=1
!   enddo
! enddo

!=======================
!==  Define integers and reals
!=======================
 paw_dmft%fermie_lda=fermie_lda ! in Ha
 paw_dmft%fermie= fermie_lda
 if(nspinor==1) then
   if(paw_dmft%nsppol==2) then
     paw_dmft%nelectval= dtset%nelect-float(paw_dmft%dmftbandi-1)*paw_dmft%nsppol
   else if(paw_dmft%nsppol==1) then
     paw_dmft%nelectval= dtset%nelect-float(paw_dmft%dmftbandi-1)*paw_dmft%nsppol*2
   endif
 else if (nspinor==2) then
   paw_dmft%nelectval= dtset%nelect-float(paw_dmft%dmftbandi-1)*paw_dmft%nsppol
 endif
 paw_dmft%filapp= fnametmp_app
 paw_dmft%filnamei= fnamei
 paw_dmft%natpawu=dtset%natpawu
 paw_dmft%natom=dtset%natom
 paw_dmft%temp=dtset%tsmear!*unit_e
 paw_dmft%dmft_iter=dtset%dmft_iter
 paw_dmft%dmft_kspectralfunc=dtset%dmft_kspectralfunc
 paw_dmft%dmft_dc=dtset%dmft_dc
 !paw_dmft%idmftloop=0
 paw_dmft%prtvol = dtset%prtvol
 paw_dmft%prtdos = dtset%prtdos
 paw_dmft%dmft_tolfreq = dtset%dmft_tolfreq
 paw_dmft%dmft_lcpr = dtset%dmft_tollc
 paw_dmft%dmft_charge_prec = dtset%dmft_charge_prec

!=======================
!==  Fixed self for input
!=======================
 paw_dmft%use_fixed_self=dtset%usedmatpu
 paw_dmft%fixed_self=>dmatpawu

!=======================
!==  Choose solver
!=======================

 paw_dmft%dmft_solv=dtset%dmft_solv
 paw_dmft%dmft_blockdiag=0
 if(paw_dmft%dmft_solv==-2) then
   paw_dmft%dmft_solv=2
   paw_dmft%dmft_blockdiag=1
 endif
!  0: LDA, no solver
!  1: LDA+U
! -1: LDA+U but LDA values are not renormalized !
! if((paw_dmft%dmft_solv==0.and.paw_dmft%prtvol>4).or.&
!&   (paw_dmft%dmft_solv>=-1.and.paw_dmft%dmft_solv<=2)) then
!   call wrtout(std_out,message,'COLL')
!   call wrtout(ab_out,message,'COLL')
! endif

 if(paw_dmft%dmft_solv==0) then
   do itypat=1,psps%ntypat
     if(pawtab(itypat)%lpawu/=-1) then
       if((pawtab(itypat)%upawu)>tol5.or.(pawtab(itypat)%jpawu)>tol5) then
          write(message, '(2a,i5,2a,2e15.6)' )ch10,&
&          ' option dmft_solv=0 requires upaw=jpaw=0 for species',itypat,ch10,&
&          ' Value of upawu and jpawu are here',pawtab(itypat)%upawu,pawtab(itypat)%jpawu
          MSG_ERROR(message)
        endif
     endif
   enddo
 endif

! todo_ab: why upaw and jpawu are not zero (on bigmac) if lpawu==-1 ?
! if(paw_dmft%dmft_solv==0.and.&
!& (maxval(abs(pawtab(:)%upawu))>tol5.or.maxval(abs(pawtab(:)%jpawu))>tol5)) then
!   write(message, '(a,a,2f12.3)' )ch10,&
!&   ' option dmft_solv=0 requires upaw=jpaw=0',maxval(abs(pawtab(:)%upawu)),maxval(abs(pawtab(:)%jpawu))
!    MSG_WARNING(message)
! endif

 paw_dmft%dmftcheck=dtset%dmftcheck

 if(paw_dmft%dmftcheck==-1) then
   message = ' init_dmft: dmftcheck=-1 should not happend here'
   MSG_BUG(message)
 endif
 paw_dmft%dmft_log_freq=1 ! use logarithmic frequencies.
 if(paw_dmft%dmft_solv==6.or.paw_dmft%dmft_solv==7) then
   paw_dmft%dmft_log_freq=0 ! do not use logarithmic frequencies.
 endif
 paw_dmft%dmft_nwli=dtset%dmft_nwli
 if(paw_dmft%dmft_log_freq==1) then
   paw_dmft%dmft_nwlo=dtset%dmft_nwlo
 else
   paw_dmft%dmft_nwlo=dtset%dmft_nwli
 endif
 paw_dmft%dmft_nwr=800
 paw_dmft%dmft_rslf=dtset%dmft_rslf
 paw_dmft%dmft_mxsf=dtset%dmft_mxsf
 paw_dmft%dmftqmc_l=dtset%dmftqmc_l
 paw_dmft%dmftqmc_n=dtset%dmftqmc_n
 paw_dmft%dmftqmc_seed=dtset%dmftqmc_seed
 paw_dmft%dmftqmc_therm=dtset%dmftqmc_therm
 paw_dmft%dmftqmc_t2g=dtset%dmft_t2g
 paw_dmft%dmftqmc_x2my2d=dtset%dmft_x2my2d
 paw_dmft%dmftctqmc_basis =dtset%dmftctqmc_basis
 paw_dmft%dmftctqmc_check =dtset%dmftctqmc_check
 paw_dmft%dmftctqmc_correl=dtset%dmftctqmc_correl
 paw_dmft%dmftctqmc_gmove =dtset%dmftctqmc_gmove
 paw_dmft%dmftctqmc_grnns =dtset%dmftctqmc_grnns
 paw_dmft%dmftctqmc_meas  =dtset%dmftctqmc_meas
 paw_dmft%dmftctqmc_mrka  =dtset%dmftctqmc_mrka
 paw_dmft%dmftctqmc_mov   =dtset%dmftctqmc_mov
 paw_dmft%dmftctqmc_order =dtset%dmftctqmc_order
 paw_dmft%dmftctqmc_triqs_nleg =dtset%dmftctqmc_triqs_nleg

 if ( paw_dmft%dmft_solv >= 4 ) then
 write(message, '(a,a,i6)' )ch10,&
&   '=> Seed for QMC inside DMFT is dmftqmc_seed=',paw_dmft%dmftqmc_seed
   call wrtout(std_out,message,'COLL')
 endif


!=======================
!==  Variables for DMFT itself
!=======================

 mbandc = paw_dmft%mbandc

 ABI_ALLOCATE(paw_dmft%eigen_lda,(paw_dmft%nsppol,paw_dmft%nkpt,paw_dmft%mbandc))
 paw_dmft%eigen_lda=zero

! allocate(paw_dmft%wtk(paw_dmft%nkpt))
 paw_dmft%wtk=>dtset%wtk
 if(dtset%iscf<0) then
   paw_dmft%wtk=one/float(dtset%nkpt)
 endif
 sumwtk=0
 do ikpt=1,paw_dmft%nkpt
   sumwtk=sumwtk+paw_dmft%wtk(ikpt)
 enddo
 if(abs(sumwtk-1_dp)>tol13.and.dtset%iscf>=0) then
   write(message, '(a,f15.10)' )' sum of k-point is incorrect',sumwtk
   MSG_ERROR(message)
 endif
 ABI_ALLOCATE(paw_dmft%lpawu,(paw_dmft%natom))
 do iatom=1,paw_dmft%natom
   paw_dmft%lpawu(iatom)=pawtab(typat(iatom))%lpawu
   if(paw_dmft%dmftqmc_t2g==1.and.paw_dmft%lpawu(iatom)==2) paw_dmft%lpawu(iatom)=1
   if(paw_dmft%dmftqmc_x2my2d==1.and.paw_dmft%lpawu(iatom)==2) paw_dmft%lpawu(iatom)=0
 enddo
 paw_dmft%maxlpawu=maxval(paw_dmft%lpawu(:))
 maxlpawu = paw_dmft%maxlpawu

 ABI_ALLOCATE(paw_dmft%psichi,(nsppol,dtset%nkpt,mbandc,nspinor,dtset%natom,(2*maxlpawu+1)))

 paw_dmft%psichi=cmplx(zero,zero,kind=dp)
 paw_dmft%lpsichiortho=0


!=======================
! Real      frequencies
!=======================

 iexist2=1
 if(dtset%iscf<0.and.(paw_dmft%dmft_solv==5.or.paw_dmft%dmft_solv==8)) then
     tmpfil = trim(paw_dmft%filapp)//'_spectralfunction_realfrequencygrid'
     inquire(file=trim(tmpfil),exist=lexist)!,recl=nrecl)
   !  write(6,*) "inquire",lexist
   grid_unt=2000
     if((.not.lexist)) then
       iexist2=0
       write(message,'(4x,a,i5,3a)') "File number",grid_unt,&
&       " called ",trim(tmpfil)," does not exist"
       call wrtout(std_out,message,'COLL')
       message = "Cannot continue: the missing file coming from Maxent code is needed"
       MSG_WARNING(message)
     endif

     if(iexist2==1) then
#ifdef FC_NAG
       open (unit=grid_unt,file=trim(tmpfil),status='unknown',form='formatted',recl=ABI_RECL)
#else
       open (unit=grid_unt,file=trim(tmpfil),status='unknown',form='formatted')
#endif
       rewind(grid_unt)
      ! if (open_file(tmpfil, message, newunit=grid_unt, status='unknown', form='formatted') /= 0) then
      !   MSG_ERROR(message)
      ! end if
       write(message,'(3a)') ch10,"  == Read  grid frequency in file ",trim(tmpfil)
       call wrtout(std_out,message,'COLL')
       write(message,'(a,a,a,i4)') 'opened file : ', trim(tmpfil), ' unit', grid_unt
       call wrtout(std_out,message,'COLL')
       read(grid_unt,*) ngrid
       ABI_ALLOCATE(paw_dmft%omega_r,(ngrid))
       if(ioerr<0) then
         message = "Error reading grid file"
         MSG_ERROR(message)
       endif
       do ifreq=1,ngrid
         read(grid_unt,*) paw_dmft%omega_r(ifreq)
         paw_dmft%omega_r(ifreq)=paw_dmft%omega_r(ifreq)
       enddo
       if(ioerr<0) then
         message = "Error reading grid file"
         MSG_ERROR(message)
       endif
     endif
 else
  ABI_ALLOCATE(paw_dmft%omega_r,(2*paw_dmft%dmft_nwr))
  ! Set up real frequencies for spectral function in Hubbard one.
   step=0.00005_dp
   paw_dmft%omega_r(2*paw_dmft%dmft_nwr)=pi*step*(two*float(paw_dmft%dmft_nwr-1)+one)
   do ifreq=1,2*paw_dmft%dmft_nwr-1
    paw_dmft%omega_r(ifreq)=pi*step*(two*float(ifreq-1)+one)-paw_dmft%omega_r(2*paw_dmft%dmft_nwr)
  !  write(std_out,*) ifreq,paw_dmft%omega_r(ifreq)
   enddo

 endif
!=======================
! Imaginary frequencies
!=======================
! Set up log frequencies
 if(dtset%ucrpa==0) call construct_nwlo_dmft(paw_dmft)

!=========================================================
!== if we use ctqmc impurity solver
!=========================================================
! IMPORTANT : paw_dmft%hybrid is corrupted somewhere in DMFT routines on
! tikal_psc and max2_open64. Use a local hybrid in qmc_prep even if not optimal.
! Anyway Initializing ctqmc here is not good and produce the same result for
! dmft_iter =1 which speed up the convergency ...
! FIXME : Move this to init_sc_dmft and find bug
 if(paw_dmft%dmft_solv==5) then ! CTQMC initialisation
   write(message,'(a,2x,a,f13.5)') ch10,&
&  " == Initializing CTQMC"
!   call wrtout(std_out,message,'COLL')

   ABI_DATATYPE_ALLOCATE(paw_dmft%hybrid,(paw_dmft%natom))
   do iatom=1,paw_dmft%natom
     if(paw_dmft%lpawu(iatom)/=-1) then
       nflavor=2*(2*paw_dmft%lpawu(iatom)+1)
#ifdef HAVE_MPI
       call CtqmcInterface_init(paw_dmft%hybrid(iatom),paw_dmft%dmftqmc_seed,paw_dmft%dmftqmc_n, &
&       paw_dmft%dmftqmc_therm, paw_dmft%dmftctqmc_meas,nflavor,paw_dmft%dmftqmc_l,one/paw_dmft%temp,zero,&
&       std_out,paw_dmft%spacecomm)
#else
       call CtqmcInterface_init(paw_dmft%hybrid(iatom),paw_dmft%dmftqmc_seed,paw_dmft%dmftqmc_n, &
&       paw_dmft%dmftqmc_therm, paw_dmft%dmftctqmc_meas,nflavor,paw_dmft%dmftqmc_l,one/paw_dmft%temp,zero,&
&       std_out)
#endif
       call CtqmcInterface_setOpts(paw_dmft%hybrid(iatom),&
                                   opt_Fk      =1,&
&                                  opt_order   =paw_dmft%dmftctqmc_order ,&
&                                  opt_movie   =paw_dmft%dmftctqmc_mov   ,&
&                                  opt_analysis=paw_dmft%dmftctqmc_correl,&
&                                  opt_check   =paw_dmft%dmftctqmc_check ,&
&                                  opt_noise   =paw_dmft%dmftctqmc_grnns ,&
&                                  opt_spectra =paw_dmft%dmftctqmc_mrka  ,&
&                                  opt_gmove   =paw_dmft%dmftctqmc_gmove )
     end if
   enddo
   write(message,'(a,2x,a,f13.5)') ch10,&
&  " == Initialization CTQMC done"
   !call wrtout(std_out,message,'COLL')
 endif

 if(paw_dmft%dmftcheck==1.and.paw_dmft%dmft_solv<4) then
   paw_dmft%dmftqmc_l=64
 endif

!************************************************************************
end subroutine init_dmft
!!***

!!****f* m_paw_dmft/construct_nwli_dmft
!! NAME
!! construct_nwli_dmft
!!
!! FUNCTION
!!  Compute linear frequencies
!!
!! INPUTS
!!  paw_dmft=structure for dmft
!!  nwli=number of linear frequencies
!!
!! OUTPUTS
!!  omegali(1:nwli)=computed frequencies
!!
!! PARENTS
!!      m_green,m_paw_dmft
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE
!!

subroutine construct_nwli_dmft(paw_dmft,nwli,omega_li)

 type(paw_dmft_type), intent(in) :: paw_dmft
 integer, intent(in) :: nwli
 real(dp), intent(out) :: omega_li(:)
 !fortran2003 ?
 !real(dp), allocatable, intent(inout) :: omega_li(:)
 integer :: ifreq
 real(dp) :: factor
 character(len=100) :: message

! if (allocated(omega_li)) then
   if (size(omega_li) .ne. nwli) then
     write(message,'(2a,i8,a,i8)') ch10, "Number of linear frequencies asked is", &
       &    nwli, "whereas dimension of array omega_li is", size(omega_li)
     MSG_BUG(message)
!     ABI_DEALLOCATE(omega_li)
!     ABI_ALLOCATE(omega_li,(nwli))
!     write(*,*) "RESIZE"
!     call flush(6)
   endif
!     write(*,*) "NOTHING"
!     call flush(6)
! else
!     write(*,*) "ALLOCATE"
!     call flush(6)
!   ABI_ALLOCATE(omega_li,(nwli))
! endif

! Set up linear frequencies
 factor = pi*paw_dmft%temp
 do ifreq=1,nwli
   omega_li(ifreq)=factor*(real(2*ifreq-1,kind=dp))
   ! (2(ifreq-1)+1 = 2ifreq-1
 enddo
end subroutine construct_nwli_dmft
!!***

!!****f* m_paw_dmft/construct_nwlo_dmft
!! NAME
!! construct_nwlo_dmft
!!
!! FUNCTION
!!  Allocate log frequencies if used and compute them as well as their weight
!!
!! INPUTS
!!  paw_dmft=structure for dmft calculation
!!
!!
!! PARENTS
!!      m_paw_dmft
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE
!!
!! NOTE
!! The part of the code which deals
!! with the use of logarithmic frequencies
!! is a modification of the GNU GPL
!! code available on http://dmft.rutgers.edu/ and
!! described in the  RMP paper written by
!! G.Kotliar,  S.Y.Savrasov, K.Haule, V.S.Oudovenko, O.Parcollet, C.A.Marianetti.
!!

subroutine construct_nwlo_dmft(paw_dmft)
 use m_splines

 type(paw_dmft_type), intent(inout) :: paw_dmft
 integer :: cubic_freq
 integer :: ifreq,ifreq2
 ! for parallel
 integer :: myproc, nproc, spacecomm
 integer :: deltaw, residu, omegaBegin, omegaEnd
 ! end
 character(len=500) :: message
 real(dp) :: deltaomega,expfac,omegamaxmin,prefacexp,AA,BB,CC,nlin,nlog,t1
 real(dp) :: wl
 complex(dpc):: ybcbeg,ybcend
 integer, allocatable :: select_log(:)
 real(dp), allocatable :: omega_lo_tmp(:)
 real(dp), allocatable :: omega_li(:)
 real(dp), allocatable :: wgt_wlo(:)
 complex(dpc), allocatable :: tospline_lo(:), splined_li(:),ysplin2_lo(:)

 ABI_ALLOCATE(omega_lo_tmp,(paw_dmft%dmft_nwlo))
 ABI_ALLOCATE(wgt_wlo,(paw_dmft%dmft_nwlo))

!==  Variables for DMFT related to frequencies
! the part of the code which deals
! with the use of logarithmic frequencies
! is a modification of the GNU GPL
! code available on http://dmft.rutgers.edu/ and
! described in the  RMP paper written by
! G.Kotliar, S.Y.Savrasov, K.Haule, V.S.Oudovenko, O.Parcollet, C.A.Marianetti

!========================================
!== construct log. freq.
 if(paw_dmft%dmft_log_freq==1) then
!=======================================
   cubic_freq=0
   !omegamaxmin=paw_dmft%omega_li(paw_dmft%dmft_nwli)-paw_dmft%omega_li(paw_dmft%dmftqmc_l+1)
   omegamaxmin=pi*paw_dmft%temp*two*real(paw_dmft%dmft_nwli-paw_dmft%dmftqmc_l-1,kind=dp)

   if(cubic_freq==1) then

     if (paw_dmft%dmft_solv .eq. 5 ) then
       write(message, '(2a)') ch10, "Warning : Cubish Mesh not tested with CT-QMC"
       MSG_WARNING(message)
     end if
!  ------------  CUBIC MESH MESH
!    useless
     nlin=dble(paw_dmft%dmft_nwli)
     nlog=dble(paw_dmft%dmft_nwlo)
     AA=((nlin-one)/nlin/(nlog**2-one)-one/(three*nlin))/((nlog**3-one)/(nlog**2-one)-seven/three)
     BB=(one/nlin - seven*AA)/three
     CC=-AA-BB
!    AA=((nlin-one)/nlin/(nlog-one)-one/(nlin))/((nlog**2-one)/(nlog-one)-three)
!    BB=(one/nlin - three*AA)
!    CC=-AA-BB
     write(message, '(a,16x,2(2x,a))') ch10,"  Cubic Mesh Parameters are"
     call wrtout(std_out,message,'COLL')
     write(message, '(3x,a,3(2x,e13.5))') "AA,BB,CC",AA,BB,CC
     call wrtout(std_out,message,'COLL')
     do ifreq=1,paw_dmft%dmft_nwlo
       t1=dble(ifreq)
       !omega_lo_tmp(ifreq)=(AA*t1**3+BB*t1**2+CC)*omegamaxmin+paw_dmft%omega_li(1)
       omega_lo_tmp(ifreq)=(AA*t1**3+BB*t1**2+CC)*omegamaxmin+paw_dmft%temp*pi
!       paw_dmft%omega_lo(ifreq)=(AA*t1**2+BB*t1+CC)*omegamaxmin+paw_dmft%omega_li(1)
!     write(69,*) paw_dmft%omega_lo(ifreq),0.5
     enddo
   else
     if(paw_dmft%dmft_solv<4) paw_dmft%dmftqmc_l=0

!   ------------  LOGARITHMIC MESH
     deltaomega=0.5_dp
     expfac=log(omegamaxmin/deltaomega)/(float(paw_dmft%dmft_nwlo-paw_dmft%dmftqmc_l-1)/two)
     prefacexp=omegamaxmin/(exp(expfac*float(paw_dmft%dmft_nwlo-paw_dmft%dmftqmc_l-1))-one)
     ABI_ALLOCATE(select_log,(paw_dmft%dmft_nwlo))
     select_log=0

!   ------------ IMPOSE LINEAR MESH for w < 2*w_n=(2*l-1)pi/beta
!         Check variables (Already done in chkinp if dmft_solv==5)
     if (paw_dmft%dmftqmc_l .gt. paw_dmft%dmft_nwlo) then
       write(message, '(a,a,i6)' )ch10,&
&       ' ERROR: dmft_nwlo has to be at leat equal to 2xdmftqmc_l :',2*paw_dmft%dmftqmc_l
       MSG_ERROR(message)
     end if
!         End Check

     call construct_nwli_dmft(paw_dmft,paw_dmft%dmftqmc_l,omega_lo_tmp(1:paw_dmft%dmftqmc_l))
     select_log(1:paw_dmft%dmftqmc_l) = (/ (ifreq,ifreq=1,paw_dmft%dmftqmc_l) /)

     !do ifreq=1,paw_dmft%dmftqmc_l
     !  omega_lo_tmp(ifreq)=(two*DBLE(ifreq-1)+one)*pi*paw_dmft%temp
     !  select_log(ifreq)=ifreq
     !enddo

!   ------------ COMPLETE FREQUENCIES WITH LOG MESH
     wl = paw_dmft%temp*pi*real(2*paw_dmft%dmftqmc_l+1,kind=dp)
     do ifreq=1,paw_dmft%dmft_nwlo-paw_dmft%dmftqmc_l
       !omega_lo_tmp(ifreq+paw_dmft%dmftqmc_l)=prefacexp*(exp(expfac*float(ifreq-1))-one)+paw_dmft%omega_li(paw_dmft%dmftqmc_l+1)
       omega_lo_tmp(ifreq+paw_dmft%dmftqmc_l)=prefacexp*(exp(expfac*real(ifreq-1,kind=dp))-one) + wl
!       -------- Impose that the each frequency of the logarithmic mesh is on a Matsubara frequency
! FIXME : This may be done for all solver, not only for QMCs
       if(paw_dmft%dmft_solv>=4) then
         ! Compute the integer "n" of iwn
         ifreq2 = nint((omega_lo_tmp(ifreq+paw_dmft%dmftqmc_l)/(paw_dmft%temp*pi)-one)*half)
         ! compute freq
         omega_lo_tmp(ifreq+paw_dmft%dmftqmc_l)= (dble(ifreq2)*two+one)*pi*paw_dmft%temp

         if((ifreq2+1)>paw_dmft%dmft_nwli) then
           write(message, '(a,a,i8)' )ch10,&
&          ' BUG: init_dmft,   dimension  of array select_log is about to be overflown',&
&          (ifreq2+1)
           MSG_BUG(message)
         endif
         select_log(paw_dmft%dmftqmc_l+ifreq)=ifreq2+1
       endif
     enddo

!       -------- Suppress duplicate frequencies
! FIXME : So this also should be done for all solver and remove useless
! frequencies
     if(paw_dmft%dmft_solv>=4) then
       ifreq2=1
       do ifreq=2,paw_dmft%dmft_nwlo-1
         if(select_log(ifreq2).ne.select_log(ifreq)) then
           ifreq2=ifreq2+1
           omega_lo_tmp(ifreq2)=omega_lo_tmp(ifreq)
           select_log(ifreq2)=select_log(ifreq)
         endif
       enddo
       paw_dmft%dmft_nwlo=ifreq2+1
     endif

     omega_lo_tmp(1)=paw_dmft%temp*pi
     omega_lo_tmp(paw_dmft%dmft_nwlo)=paw_dmft%temp*pi*real(2*paw_dmft%dmft_nwli-1,kind=dp)
  endif

!=======================
!== construct weight for log. freq.
!=======================
  ABI_ALLOCATE(tospline_lo,(paw_dmft%dmft_nwlo))
  ABI_ALLOCATE(splined_li,(paw_dmft%dmft_nwli))
  ABI_ALLOCATE(ysplin2_lo,(paw_dmft%dmft_nwlo))
  if (allocated(omega_li)) then
    ABI_DEALLOCATE(omega_li)
  endif
  ABI_ALLOCATE(omega_li,(1:paw_dmft%dmft_nwli))
  call construct_nwli_dmft(paw_dmft,paw_dmft%dmft_nwli,omega_li)

  !Parallelisation over frequencies!
  ! ============= Set up =============
  myproc = paw_dmft%myproc
  nproc = paw_dmft%nproc
  spacecomm = paw_dmft%spacecomm
  deltaw = paw_dmft%dmft_nwlo / nproc
  residu = paw_dmft%dmft_nwlo - nproc*deltaw
  if ( myproc .LT. nproc - residu ) then
    omegaBegin = 1 + myproc*deltaw
    omegaEnd   = (myproc + 1)*deltaw
  else
    omegaBegin = 1 + myproc*(deltaw + 1) -nproc + residu
    omegaEnd = omegaBegin + deltaw
  end if
  wgt_wlo(:)=zero
  ! ============= END Set up =============
  do ifreq=omegaBegin,omegaEnd
    tospline_lo=cmplx(0_dp,0_dp,kind=dp)
!    do ifreq1=1,paw_dmft%dmft_nwlo
    tospline_lo(ifreq)=cmplx(1_dp,0_dp,kind=dp)
!    tospline_lo(ifreq1)=ifreq1**2-ifreq1
!    enddo
    splined_li=cmplx(0_dp,0_dp,kind=dp)
!    ybcbeg=cmplx(one/tol16**2,zero)
!    ybcend=cmplx(one/tol16**2,zero)
    ybcbeg=czero
    ybcend=czero


!==         spline delta function
    call spline_complex( omega_lo_tmp, tospline_lo, paw_dmft%dmft_nwlo, &
   & ybcbeg, ybcend, ysplin2_lo)
!   do ifreq1=1,paw_dmft%dmft_nwlo
!    write(6588,*) paw_dmft%omega_lo(ifreq1),ysplin2_lo(ifreq1)
!   enddo

    call splint_complex( paw_dmft%dmft_nwlo, omega_lo_tmp, tospline_lo,&
   & ysplin2_lo, paw_dmft%dmft_nwli, omega_li, splined_li)

!==         accumulate weights
    wgt_wlo(ifreq)=zero
    do ifreq2=1,paw_dmft%dmft_nwli
      wgt_wlo(ifreq)=wgt_wlo(ifreq)+real(splined_li(ifreq2),kind=dp)
    enddo
! do ifreq1=1,paw_dmft%dmft_nwlo
!  write(6688,*) paw_dmft%omega_lo(ifreq1),tospline_lo(ifreq1)
! enddo
! do ifreq1=1,paw_dmft%dmft_nwli
!  write(6788,*) paw_dmft%omega_li(ifreq1),splined_li(ifreq1)
  enddo
  ! ============= Gatherall  =============
  call xmpi_sum(wgt_wlo, spacecomm, residu)
  ! ============= END Gatherall ==========
  ! end parallelisation over frequencies

  ABI_DEALLOCATE(tospline_lo)
  ABI_DEALLOCATE(splined_li)
  ABI_DEALLOCATE(ysplin2_lo)
! if(abs(dtset%pawprtvol)>=3) then
    write(message, '(a,18x,2(2x,a21))') ch10,"Log. Freq","weight"
    call wrtout(std_out,message,'COLL')
    do ifreq=1,paw_dmft%dmft_nwlo
      write(message, '(3x,a9,i6,2(2x,e21.14))') "--ifreq--",ifreq,omega_lo_tmp(ifreq),wgt_wlo(ifreq)
      call wrtout(std_out,message,'COLL')
    enddo
    write(message, '(3x,a,i6)') "  Total number of log frequencies is", paw_dmft%dmft_nwlo
    call wrtout(std_out,message,'COLL')
    ifreq2 = 1
    do ifreq=1,min(30,paw_dmft%dmft_nwlo)
      write(message, '(3x,a9,i6,2(2x,e21.14))') "--ifreq--",ifreq,omega_li(ifreq)
      call wrtout(std_out,message,'COLL')
      if ( select_log(ifreq2) .eq. ifreq ) then
        write(message, '(3x,a,i4,2(2x,i5))') "--sel_log",1
        ifreq2 = ifreq+1
      else
        write(message, '(3x,a,i4,2(2x,i5))') "--sel_log",0
      end if
      call wrtout(std_out,message,'COLL')
    enddo
    write(message, '(3x,2a)') "--ifreq--","..."
    call wrtout(std_out,message,'COLL')
    write(message, '(3x,a,i6,2(2x,e13.5))') "--ifreq--",paw_dmft%dmft_nwli,omega_li(paw_dmft%dmft_nwli)
    call wrtout(std_out,message,'COLL')
!   endif
   ABI_DEALLOCATE(select_log)
   ABI_DEALLOCATE(omega_li)

!=========================================================
!== do not construct log. freq. and use linear frequencies
 else
!=========================================================
   write(message, '(a,10x,2(2x,a))') ch10,"   Use of linear frequencies for DMFT calculation"
   call wrtout(std_out,message,'COLL')
   call construct_nwli_dmft(paw_dmft,paw_dmft%dmft_nwli,omega_lo_tmp)
   wgt_wlo=one
 endif

 ! Should be check but since type definition does not initialize pointer with
 ! =>null() (fortran95 and later) it produces conditional jump in valgrind
 !if ( associated(paw_dmft%omega_lo) ) then
 !  ABI_DEALLOCATE(paw_dmft%omega_lo)
 !endif
 !if ( associated(paw_dmft%wgt_wlo) ) then
 !  ABI_DEALLOCATE(paw_dmft%wgt_wlo)
 !endif
 ABI_ALLOCATE(paw_dmft%omega_lo,(paw_dmft%dmft_nwlo))
 ABI_ALLOCATE(paw_dmft%wgt_wlo,(paw_dmft%dmft_nwlo))
 paw_dmft%omega_lo(1:paw_dmft%dmft_nwlo) = omega_lo_tmp(1:paw_dmft%dmft_nwlo)
 paw_dmft%wgt_wlo(1:paw_dmft%dmft_nwlo) = wgt_wlo(1:paw_dmft%dmft_nwlo)
 ABI_DEALLOCATE(omega_lo_tmp)
 ABI_DEALLOCATE(wgt_wlo)
end subroutine construct_nwlo_dmft
!!***

!!****f* m_paw_dmft/destroy_dmft
!! NAME
!! destroy_dmft
!!
!! FUNCTION
!!  deallocate some variables related to paw_dmft
!!
!! INPUTS
!!  paw_dmft
!!
!! OUTPUT
!!
!! PARENTS
!!      outscfcv,vtorho
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine destroy_dmft(paw_dmft)

!Arguments ------------------------------------
!scalars
 type(paw_dmft_type),intent(inout) :: paw_dmft

!Local variables-------------------------------
 integer :: iatom

! *********************************************************************

   if (paw_dmft%dmft_solv == 5 .and. allocated(paw_dmft%hybrid)) then
     do iatom=1, size(paw_dmft%hybrid) !paw_dmft%natom
       !if(paw_dmft%lpawu(iatom)/=-1) then
         call ctqmcinterface_finalize(paw_dmft%hybrid(iatom))
       !endif
     enddo
     ABI_DATATYPE_DEALLOCATE(paw_dmft%hybrid)
   endif
   if (allocated(paw_dmft%psichi))  then
     ABI_DEALLOCATE(paw_dmft%psichi)
   end if
!   paw_dmft%wtk is only an explicit pointer =>dtset%wtk
!   if (associated(paw_dmft%wtk)) deallocate(paw_dmft%wtk)
   paw_dmft%wtk => null()
   paw_dmft%fixed_self => null()
   if (allocated(paw_dmft%eigen_lda))  then
     ABI_DEALLOCATE(paw_dmft%eigen_lda)
   endif
   if (associated(paw_dmft%omega_lo))  then
     ABI_DEALLOCATE(paw_dmft%omega_lo)
   end if
   if (associated(paw_dmft%omega_r))  then
     ABI_DEALLOCATE(paw_dmft%omega_r)
   end if
   if (associated(paw_dmft%wgt_wlo))  then
     ABI_DEALLOCATE(paw_dmft%wgt_wlo)
   end if
   if (allocated(paw_dmft%lpawu))  then
     ABI_DEALLOCATE(paw_dmft%lpawu)
   end if

end subroutine destroy_dmft
!!***

!!****f* m_paw_dmft/destroy_sc_dmft
!! NAME
!! destroy_sc_dmft
!!
!! FUNCTION
!!  deallocate paw_dmft
!!
!! INPUTS
!!  paw_dmft
!!
!! OUTPUT
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine destroy_sc_dmft(paw_dmft)

!Arguments ------------------------------------
!scalars
 type(paw_dmft_type),intent(inout) :: paw_dmft

!Local variables-------------------------------
 character(len=500) :: message

! *********************************************************************

 if (( .not. allocated(paw_dmft%occnd) .or. .not. allocated(paw_dmft%band_in) &
&  .or. .not. allocated(paw_dmft%include_bands) .or. .not. allocated(paw_dmft%exclude_bands)) &
&  .and. paw_dmft%use_dmft == 1 )  then
  write(message, '(a,a,a)' )&
&  '  an array is not allocated and is not deallocated with use_dmft==1 ',ch10, &
&  '  Action : check the code'
  MSG_WARNING(message)
 endif
 if ( allocated(paw_dmft%occnd) )          then
   ABI_DEALLOCATE(paw_dmft%occnd)
 end if
 if ( allocated(paw_dmft%band_in) )        then
   ABI_DEALLOCATE(paw_dmft%band_in)
 end if
 if ( allocated(paw_dmft%include_bands) )  then
   ABI_DEALLOCATE(paw_dmft%include_bands)
 end if
 if ( allocated(paw_dmft%exclude_bands) )  then
   ABI_DEALLOCATE(paw_dmft%exclude_bands)
 end if

 call destroy_sc_dmft_paralkgb(paw_dmft)

end subroutine destroy_sc_dmft
!!***

!!****f* m_paw_dmft/print_dmft
!! NAME
!! print_dmft
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      outscfcv,vtorho
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_dmft(paw_dmft,pawprtvol)

!Arguments ------------------------------------
!type
 type(paw_dmft_type),intent(in) :: paw_dmft
 integer :: pawprtvol

!Local variables-------------------------------
 integer :: ikpt,iband,ifreq,isppol
 character(len=500) :: message
! *********************************************************************

 if( abs(pawprtvol) >= 3 )  then
  write(message,'(4a,3(a,2x,e21.14,a))') &
&   "  -------------------------------------------------",ch10,&
&   "  --- Data for DMFT ",ch10,&
&   "  --- paw_dmft%fermie     = ",paw_dmft%fermie    ,ch10,&
&   "  --- paw_dmft%fermie_lda = ",paw_dmft%fermie_lda,ch10,&
&   "  --- paw_dmft%temp       = ",paw_dmft%temp      ,ch10
  call wrtout(std_out,message,'COLL')
  write(message,'(7(a,15x,i8,a),a,2x,e21.14,2a)') &
&   "  --- paw_dmft%natpawu    = ",paw_dmft%natpawu   ,ch10,&
&   "  --- paw_dmft%dmft_iter  = ",paw_dmft%dmft_iter ,ch10,&
&   "  --- paw_dmft%dmft_solv  = ",paw_dmft%dmft_solv ,ch10,&
&   "  --- paw_dmft%dmft_nwlo  = ",paw_dmft%dmft_nwlo ,ch10,&
&   "  --- paw_dmft%dmft_nwli  = ",paw_dmft%dmft_nwli ,ch10,&
&   "  --- paw_dmft%dmft_dc    = ",paw_dmft%dmft_dc   ,ch10,&
&   "  --- paw_dmft%dmftqmc_l  = ",paw_dmft%dmftqmc_l ,ch10,&
&   "  --- paw_dmft%dmftqmc_n  = ",paw_dmft%dmftqmc_n ,ch10,&
&   "  -------------------------------------------------"
  call wrtout(std_out,message,'COLL')

!  write(message,'(4a,3(a,2x,f8.3,a),8(a,2x,i8,a),a)') "-----------------------------------------------",ch10,&
!&   "--- Data for DMFT ",ch10,&
!&   "--- paw_dmft%fermie     = ",paw_dmft%fermie    ,ch10,&
!&   "--- paw_dmft%fermie_lda = ",paw_dmft%fermie_lda,ch10,&
!&   "--- paw_dmft%temp       = ",paw_dmft%temp      ,ch10,&
!&   "--- paw_dmft%natpawu    = ",paw_dmft%natpawu   ,ch10,&
!&   "--- paw_dmft%dmft_iter  = ",paw_dmft%dmft_iter ,ch10,&
!&   "--- paw_dmft%dmft_solv  = ",paw_dmft%dmft_solv ,ch10,&
!&   "--- paw_dmft%dmft_nwlo  = ",paw_dmft%dmft_nwlo ,ch10,&
!&   "--- paw_dmft%dmft_nwli  = ",paw_dmft%dmft_nwli ,ch10,&
!&   "--- paw_dmft%dmft_dc    = ",paw_dmft%dmft_dc   ,ch10,&
!&   "--- paw_dmft%dmftqmc_l  = ",paw_dmft%dmftqmc_l ,ch10,&
!&   "--- paw_dmft%dmftqmc_n  = ",paw_dmft%dmftqmc_n ,ch10,&
!&   "-----------------------------------------------"
  if(abs(pawprtvol)>10) then
   call wrtout(std_out,message,'COLL')
   write(message, '(a)') " LDA Eigenvalues "
   do isppol=1,paw_dmft%nsppol
    write(message, '(a,i4)') "--isppol--",isppol
    call wrtout(std_out,message,'COLL')
    do ikpt=1,paw_dmft%nkpt
     write(message, '(a,i4,2x,f14.5,a)') "  -k-pt--",ikpt,paw_dmft%wtk(ikpt),"(<-weight(k-pt))"

     call wrtout(std_out,message,'COLL')
     do iband=1,paw_dmft%mbandc
      write(message, '(a,i4,f10.5)') "   -iband--",iband,paw_dmft%eigen_lda(isppol,ikpt,iband)
      call wrtout(std_out,message,'COLL')
     enddo
    enddo
   enddo
   write(message, '(3x,a)') "Log. Freq"
   call wrtout(std_out,message,'COLL')
   do ifreq=1,paw_dmft%dmft_nwlo
    write(message, '(3x,a,i4,2(2x,e13.5))') "--ifreq--",ifreq,paw_dmft%omega_lo(ifreq),paw_dmft%wgt_wlo(ifreq)
    call wrtout(std_out,message,'COLL')
   enddo
  endif
 endif

end subroutine print_dmft
!!***

!!****f* m_paw_dmft/print_sc_dmft
!! NAME
!! print_sc_dmft
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine print_sc_dmft(paw_dmft,pawprtvol)

!Arguments ------------------------------------
!type
 type(paw_dmft_type),intent(in) :: paw_dmft
 integer :: pawprtvol

!Local variables-------------------------------
 integer :: iband
 character(len=500) :: message
! *********************************************************************

 if( abs(pawprtvol) >= 3 )  then
   write(message,'(5a,8(a,2x,i5,a),a)')ch10,"-----------------------------------------------",ch10,&
&    "--- Data for SC DMFT ",ch10,&
&    "--- paw_dmft%mband       = ",paw_dmft%mband,ch10,&
&    "--- paw_dmft%dmftbandf   = ",paw_dmft%dmftbandf,ch10,&
&    "--- paw_dmft%dmftbandi   = ",paw_dmft%dmftbandi,ch10,&
&    "--- paw_dmft%nkpt        = ",paw_dmft%nkpt,ch10,&
&    "--- paw_dmft%nsppol      = ",paw_dmft%nsppol,ch10,&
&    "--- paw_dmft%use_dmft    = ",paw_dmft%use_dmft,ch10,&
&    "--- paw_dmft%use_sc_dmft = ",paw_dmft%use_sc_dmft,ch10,&
&    "--- paw_dmft%mbandc      = ",paw_dmft%mbandc,ch10,&
&    "-----------------------------------------------"
   call wrtout(std_out,message,'COLL')
   write(message, '(a)') " paw_dmft%band_in"
   call wrtout(std_out,message,'COLL')
   write(message, '(100i5)') (iband,iband=1,min(paw_dmft%mband,100))
   call wrtout(std_out,message,'COLL')
   write(message, '(100L3)') (paw_dmft%band_in(iband),iband=1,min(paw_dmft%mband,100))
   call wrtout(std_out,message,'COLL')
   do iband=1,paw_dmft%mbandc
     write(message,*) "include_bands",iband,paw_dmft%include_bands(iband)
     call wrtout(std_out,message,'COLL')
   enddo
   write(message, '(a,a,i4,a)' )ch10,&
&    'The',paw_dmft%mband-paw_dmft%dmftbandf+paw_dmft%dmftbandi-1,&
&    '  Following bands are excluded from the DMFT calculation  '
   call wrtout(std_out,message,'COLL')
   write(message,'(100i5)') (paw_dmft%exclude_bands(iband),iband=1,min(paw_dmft%mband-paw_dmft%dmftbandf+paw_dmft%dmftbandi-1,100))
   call wrtout(std_out,message,'COLL')
 endif

end subroutine print_sc_dmft
!!***

!!****f* m_paw_dmft/saveocc_dmft
!! NAME
!! saveocc_dmft
!!
!! FUNCTION
!!  save occnd on disk
!!
!! INPUTS
!!  paw_dmft
!!
!! OUTPUT
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine saveocc_dmft(paw_dmft)

!Arguments ------------------------------------
 type(paw_dmft_type),intent(inout) :: paw_dmft

!Local variables-------------------------------
!scalars
 character(len=fnlen) :: tmpfil
 integer :: ib,ib1,ikpt,is,unitsaveocc
 character(len=500) :: message
! *********************************************************************
 tmpfil = trim(paw_dmft%filapp)//'_DMFTOCCND'
 if (open_file(tmpfil,message,newunit=unitsaveocc,status='unknown',form='formatted') /= 0) then
   MSG_ERROR(message)
 end if

 rewind(unitsaveocc)
 write(message,'(2a)') ch10,"  == Print DMFT non diagonal occupations on disk"
 call wrtout(std_out,message,'COLL')
 write(message,'(3a,2x,4i5)') "# natom,nsppol,mbandc,nkpt",ch10&
&              ,"####",paw_dmft%natom,paw_dmft%nsppol,paw_dmft%mbandc,paw_dmft%nkpt
 call wrtout(unitsaveocc,message,'COLL')
 do is = 1 , paw_dmft%nsppol
   do ikpt = 1, paw_dmft%nkpt
     do ib = 1, paw_dmft%mbandc
       do ib1 = 1, paw_dmft%mbandc
         write(unitsaveocc,*) is,ikpt,ib,ib1,paw_dmft%occnd(1,paw_dmft%include_bands(ib),paw_dmft%include_bands(ib1),ikpt,is),&
&         paw_dmft%occnd(2,paw_dmft%include_bands(ib),paw_dmft%include_bands(ib1),ikpt,is)
       enddo
     enddo
   enddo
 enddo
 write(message,'(3a)') "# end of record",ch10&
&              ,"####  1234 "
 call wrtout(unitsaveocc,message,'COLL')
 close(unitsaveocc)

end subroutine saveocc_dmft
!!***

!!****f* m_paw_dmft/readocc_dmft
!! NAME
!! readocc_dmft
!!
!! FUNCTION
!!  read occnd on disk
!!
!! INPUTS
!!  paw_dmft   = data structure
!!  filnam_ds3 = root for filname to read (input)
!!  filnam_ds4 = root for filname to read (output)
!!
!! OUTPUT
!!  paw_dmft: occnd
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine readocc_dmft(paw_dmft,filnam_ds3,filnam_ds4)

!Arguments ------------------------------------
 type(paw_dmft_type),intent(inout) :: paw_dmft
 character(len=fnlen) :: filnam_ds3
 character(len=fnlen) :: filnam_ds4

!Local variables-------------------------------
!scalars
 character(len=fnlen) :: tmpfil
 integer :: ib,ib1,ikpt,is,unitsaveocc,dum1,dum2,dum3,dum4,ioerr
 logical :: lexist
 character(len=500) :: message
 character(len=4) :: chtemp
! *********************************************************************
 if(paw_dmft%dmft_read_occnd==0) return
 if(paw_dmft%dmft_read_occnd==1) tmpfil=trim(filnam_ds3)//'_DMFTOCCND'
 if(paw_dmft%dmft_read_occnd==2) tmpfil=trim(filnam_ds4)//'_DMFTOCCND'
 inquire(file=trim(tmpfil),exist=lexist)!,recl=nrecl)
 unitsaveocc=679
 if (lexist) then
   if (open_file(tmpfil,message,unit=unitsaveocc,status='unknown',form='formatted') /= 0) then
     MSG_ERROR(message)
   end if
   rewind(unitsaveocc)
   write(message,'(3a)') ch10,"  == Read DMFT non diagonal occupations on disk"
   call wrtout(std_out,message,'COLL')
   read(unitsaveocc,*)
   read(unitsaveocc,*,iostat=ioerr)&
&              chtemp,dum1,dum2,dum3,dum4
   if(ioerr<0) then
     write(std_out,*) "read",dum1,dum2,dum3,dum4
   endif
   write(message,'(2a,4i4)') ch10,"  == natom, nsppol, nbandc, nkpt  read are",dum1,dum2,dum3,dum4
   call wrtout(std_out,message,'COLL')
   do is = 1 , paw_dmft%nsppol
     do ikpt = 1, paw_dmft%nkpt
       do ib = 1, paw_dmft%mbandc
         do ib1 = 1, paw_dmft%mbandc
           read(unitsaveocc,*) dum1,dum2,dum3,dum4,&
&           paw_dmft%occnd(1,paw_dmft%include_bands(ib),paw_dmft%include_bands(ib1),ikpt,is),&
&           paw_dmft%occnd(2,paw_dmft%include_bands(ib),paw_dmft%include_bands(ib1),ikpt,is)
         enddo
       enddo
     enddo
   enddo
!   write(read,'(3a)') "# end of record",ch10&
!&                ,"####  1234 "
!   call wrtout(unitsaveocc,message,'COLL')
 else
   write(message,'(2a,2x,2a)') ch10,"   File",trim(tmpfil),"is not available"
   call wrtout(std_out,message,'COLL')
   write(message,'(4a)') ch10,"  ==> DMFT Occupations not available for restart", ch10, &
&   "      -> The calculation is started with Fermi Dirac scheme for occupations"
   call wrtout(std_out,message,'COLL')
 endif

end subroutine readocc_dmft
!!***

!!****f* m_paw_dmft/init_sc_dmft_paralkgb
!! NAME
!! init_sc_dmft_paralkgb
!!
!! FUNCTION
!!  Init some values used with KGB parallelism in self consistent DMFT
!!  calculation.
!!
!! INPUTS
!!  paw_dmft   = data structure
!!
!! OUTPUT
!!  paw_dmft: bandc_proc, use_bandc
!!
!! PARENTS
!!      init_sc_dmft
!!
!! CHILDREN
!!
!! SOURCE

subroutine init_sc_dmft_paralkgb(paw_dmft,mpi_enreg)

!Arguments ------------------------------------
 type(paw_dmft_type),intent(inout) :: paw_dmft
 type(MPI_type), intent(in) :: mpi_enreg

!Local variables-------------------------------
!scalars
 integer :: nproc, ib, ibc, proc

! *********************************************************************
 nproc = mpi_enreg%nproc_band

 ABI_ALLOCATE(paw_dmft%bandc_proc,(paw_dmft%mbandc))
 ABI_ALLOCATE(paw_dmft%use_bandc,(nproc))
 paw_dmft%bandc_proc = 0
 paw_dmft%use_bandc = .false.

 do ibc=1,paw_dmft%mbandc
   ib = paw_dmft%include_bands(ibc)
   proc = mod((ib-1)/mpi_enreg%bandpp,nproc)

   paw_dmft%bandc_proc(ibc) = proc

   paw_dmft%use_bandc(proc+1) = .true.
 end do
end subroutine init_sc_dmft_paralkgb

!!***

!!****f* m_paw_dmft/destroy_sc_dmft_paralkgb
!! NAME
!! destroy_sc_dmft_paralkgb
!!
!! FUNCTION
!!   deallocate bandc_proc and use_bandc
!!
!! INPUTS
!!  paw_dmft   = data structure
!!
!! OUTPUT
!!
!! PARENTS
!!      destroy_sc_dmft
!!
!! CHILDREN
!!
!! SOURCE

subroutine destroy_sc_dmft_paralkgb(paw_dmft)

!Arguments ------------------------------------
 type(paw_dmft_type),intent(inout) :: paw_dmft
! *********************************************************************

 if ( allocated(paw_dmft%bandc_proc) )  then
   ABI_DEALLOCATE(paw_dmft%bandc_proc)
 end if

 if ( allocated(paw_dmft%use_bandc) )  then
   ABI_DEALLOCATE(paw_dmft%use_bandc)
 end if

end subroutine destroy_sc_dmft_paralkgb
!!***

END MODULE m_paw_dmft
!!***
