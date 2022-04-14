!!****m* ABINIT/m_rttddft
!! NAME
!!  m_rttddft
!!
!! FUNCTION
!!  Contains various subroutines used in RT-TDDFT
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2022 ABINIT group (FB, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!  m_rttddft_driver
!!  m_rttddft_propagate
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rttddft

 use defs_basis
 use defs_abitypes,      only: MPI_type
 use defs_datatypes,     only: pseudopotential_type

 use m_cgprj,            only: ctocprj
 use m_cgtools,          only: dotprod_g
 use m_dtfil,            only: datafiles_type
 use m_dtset,            only: dataset_type
 use m_efield,           only: efield_type
 use m_energies,         only: energies_type
 use m_getchc,           only: getcsc
 use m_fourier_interpol, only: transgrid
 use m_hamiltonian,      only: init_hamiltonian, gs_hamiltonian_type
 use m_kg,               only: getcut, getph
 use m_mkrho,            only: mkrho
 use m_mpinfo,           only: proc_distrb_cycle
 use m_nonlop,           only: nonlop
 use m_paw_an,           only: paw_an_reset_flags
 use m_paw_correlations, only: setrhoijpbe0
 use m_pawcprj,          only: pawcprj_type
 use m_paw_denpot,       only: pawdenpot
 use m_pawdij,           only: pawdij, symdij
 use m_paw_ij,           only: paw_ij_reset_flags
 use m_paw_mkrho,        only: pawmkrho
 use m_paw_nhat,         only: nhatgrid
 use m_paw_occupancies,  only: pawmkrhoij
 use m_pawrhoij,         only: pawrhoij_type, pawrhoij_free, &
                             & pawrhoij_alloc, pawrhoij_inquire_dim
 use m_paw_tools,        only: chkpawovlp
 use m_rttddft_tdks,     only: tdks_type
 use m_specialmsg,       only: wrtout
 use m_setvtr,           only: setvtr
 use m_xmpi,             only: xmpi_paral, xmpi_sum, &
                             & xmpi_comm_rank !FB-test

 implicit none

 private
!!***

 public :: rttddft_setup_ele_step
 public :: rttddft_init_hamiltonian
 public :: rttddft_calc_density
 public :: rttddft_calc_etot
 public :: rttddft_calc_eig
 public :: rttddft_calc_enl
 public :: rttddft_calc_occ
!!***

contains

!!****f* m_rttddft/rttddft_setup_ele_step
!!
!! NAME
!!  rttddft_setup_ele_step
!!
!! FUNCTION
!!  Init/Update various quantities needed before performing 
!!  propagation of KS orbitals
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  tdks <type(tdks_type)> = Main RT-TDDFT object
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_propagate/rttddft_propagate_ele
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_setup_ele_step(dtset, mpi_enreg, psps, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 type(tdks_type),            intent(inout) :: tdks

 !Local variables-------------------------------
 !scalars
 integer                   :: forces_needed
 integer                   :: my_natom
 integer                   :: optcut, optgr0, optgr1, optgr2, optrad
 integer                   :: stress_needed
 !arrays
 !real(dp),parameter        :: k0(3)=(/zero,zero,zero/)

! ***********************************************************************

 my_natom=mpi_enreg%my_natom

 !** Update various quantities that needs to be changed 
 !** after a change of xred during the nuclear step

 !Compute large sphere G^2 cut-off (gsqcut) and box / sphere ratio
 !FB: @MT Probably not needed? The box didn't change only nuclear pos..
 !if (psps%usepaw==1) then
 !   call getcut(tdks%boxcut,dtset%pawecutdg,tdks%gmet,tdks%gsqcut,dtset%iboxcut, &
 !             & std_out,k0,tdks%pawfgr%ngfft)
 !else
 !   call getcut(tdks%boxcut,dtset%ecut,tdks%gmet,tdks%gsqcut,dtset%iboxcut, &
 !             & std_out,k0,tdks%pawfgr%ngfft)
 !end if
 
 !Compute structure factor phases (exp(2Pi i G.xred)) on coarse and fine grid
 call getph(tdks%atindx,dtset%natom,tdks%pawfgr%ngfftc(1),tdks%pawfgr%ngfftc(2), &
          & tdks%pawfgr%ngfftc(3),tdks%ph1d,tdks%xred)
 if (psps%usepaw==1.and.tdks%pawfgr%usefinegrid==1) then
    call getph(tdks%atindx,dtset%natom,tdks%pawfgr%ngfft(1),tdks%pawfgr%ngfft(2), &
             & tdks%pawfgr%ngfft(3),tdks%ph1df,tdks%xred)
 else
    tdks%ph1df(:,:)=tdks%ph1d(:,:)
 end if
 
 !PAW specific
 if (psps%usepaw==1) then
    !Check for non-overlapping PAW spheres
    call chkpawovlp(dtset%natom,psps%ntypat,dtset%pawovlp,tdks%pawtab,tdks%rmet, &
                  & dtset%typat,tdks%xred)
 
    !Identify parts of the rectangular grid where the density has to be calculated
    !FB: Needed?
    optcut=0;optgr0=dtset%pawstgylm;optgr1=0;optgr2=0;optrad=1-dtset%pawstgylm
    forces_needed=0 !FB TODO needs to be changed if Ehrenfest?
    stress_needed=0
    if ((forces_needed==1).or. &
      & (dtset%xclevel==2.and.dtset%pawnhatxc>0.and.tdks%usexcnhat>0).or. &
      & (dtset%positron/=0.and.forces_needed==2)) then
       optgr1=dtset%pawstgylm; if (stress_needed==1) optrad=1; if (dtset%pawprtwf==1) optrad=1
    end if
    call nhatgrid(tdks%atindx1,tdks%gmet,my_natom,dtset%natom,tdks%nattyp,         &
                & tdks%pawfgr%ngfft,psps%ntypat,optcut,optgr0,optgr1,              &
                & optgr2,optrad,tdks%pawfgrtab,tdks%pawtab,tdks%rprimd,            &
                & dtset%typat, tdks%ucvol,tdks%xred,comm_atom=mpi_enreg%comm_atom, &
                & mpi_atmtab=mpi_enreg%my_atmtab,comm_fft=mpi_enreg%comm_fft,      &
                & distribfft=mpi_enreg%distribfft)
 endif

!!FB: @MT Needed? If yes, then don't forget to put it back in tdks_init/second_setup as well
!!if any nuclear dipoles are nonzero, compute the vector potential in real space (depends on
!!atomic position so should be done for nstep = 1 and for updated ion positions
!if ( any(abs(dtset%nucdipmom(:,:))>tol8) ) then
!   with_vectornd = 1
!else
!   with_vectornd = 0
!end if
!if(allocated(vectornd)) then
!   ABI_FREE(vectornd)
!end if
!ABI_MALLOC(vectornd,(with_vectornd*nfftf,3))
!vectornd=zero
!if(with_vectornd .EQ. 1) then
!   call make_vectornd(1,gsqcut,psps%usepaw,mpi_enreg,dtset%natom,nfftf,ngfftf,dtset%nucdipmom,&
!        & rprimd,vectornd,xred)

end subroutine rttddft_setup_ele_step

!!****f* m_rttddft/rttddft_init_hamiltonian
!!
!! NAME
!!  rttddft_init_hamiltonian
!!
!! FUNCTION
!!  Init/Update various quantities in order to set up
!!  the Hamiltonian
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  gs_hamk <type(gs_hamiltonian_type)> = Hamiltonian object
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  tdks <type(tdks_type)> = Main RT-TDDFT object
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_propagators/rttddft_propagator_er
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_init_hamiltonian(dtset, energies, gs_hamk, istep, mpi_enreg, psps, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                    intent(in)    :: istep
 type(dataset_type),         intent(inout) :: dtset
 type(energies_type),        intent(inout) :: energies
 type(gs_hamiltonian_type),  intent(out)   :: gs_hamk
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 type(tdks_type),            intent(inout) :: tdks

 !Local variables-------------------------------
 !scalars
 character(len=500)        :: msg
 integer                   :: comm
 integer,parameter         :: cplex=1
 integer,parameter         :: ipert=0
 integer                   :: initialized0
 integer                   :: istep_mix
 integer                   :: moved_atm_inside, moved_rhor
 integer                   :: my_natom
 integer                   :: nfftotf
 integer                   :: nzlmopt
 integer                   :: optene
 integer                   :: option
 integer                   :: nkxc, n1xccc, n3xccc
 integer                   :: usecprj_local
 logical                   :: calc_ewald
 logical                   :: tfw_activated
 real(dp)                  :: compch_sph
 real(dp)                  :: vxcavg
 !arrays
 real(dp),allocatable      :: grchempottn(:,:)
 real(dp),allocatable      :: grewtn(:,:)
 real(dp),parameter        :: k0(3)=(/zero,zero,zero/)
 real(dp),allocatable      :: kxc(:,:)
 real(dp)                  :: strsxc(6)
 real(dp)                  :: vpotzero(2)

! ***********************************************************************

 my_natom=mpi_enreg%my_natom

 !** Set up the potential (calls setvtr)
 !**  The following steps have been gathered in the setvtr routine:
 !**  - get Ewald energy and Ewald forces
 !**  - compute local ionic pseudopotential vpsp
 !**  - possibly compute 3D core electron density xccc3d
 !**  - possibly compute 3D core kinetic energy density
 !**  - possibly compute vxc and vhartr
 !**  - set up vtrial
 !** Only the local part of the potential is computed here
 !FB: Are the values of moved_atm_inside and moved_rhor correct?
 optene = 4; nkxc=0; moved_atm_inside=0; moved_rhor=1
 n1xccc=0;if (psps%n1xccc/=0) n1xccc=psps%n1xccc
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=tdks%pawfgr%nfft
 strsxc(:)=zero
 !FB: tfw_activated is a save variable in scfcv, should check where it appears again
 tfw_activated=.false.
 if (dtset%tfkinfunc==12) tfw_activated=.true.
 ABI_MALLOC(grchempottn,(3,dtset%natom))
 ABI_MALLOC(grewtn,(3,dtset%natom))
 ABI_MALLOC(kxc,(tdks%pawfgr%nfft,nkxc))
 calc_ewald = .false.
 if (dtset%ionmov/=0 .or. istep == tdks%first_step) calc_ewald=.true.
 !FB: Should we also add an option to avoid recomputing xccc3d?
 print*, 'FB-test: before setvtr'
 print*, 'FB-test: tdks%pawtab', associated(tdks%pawtab), size(tdks%pawtab)
 print*, 'FB-test: tdks%pawrad', associated(tdks%pawrad), size(tdks%pawrad)
 call setvtr(tdks%atindx1,dtset,energies,tdks%gmet,tdks%gprimd,grchempottn,         &
           & grewtn,tdks%grvdw,tdks%gsqcut,istep,kxc,tdks%pawfgr%mgfft,             &
           & moved_atm_inside,moved_rhor,mpi_enreg,tdks%nattyp,tdks%pawfgr%nfft,    &
           & tdks%pawfgr%ngfft,tdks%ngrvdw,tdks%nhat,tdks%nhatgr,tdks%nhatgrdim,    &
           & nkxc,psps%ntypat,n1xccc,n3xccc,optene,tdks%pawrad,tdks%pawtab,         &
           & tdks%ph1df,psps,tdks%rhog,tdks%rhor,tdks%rmet,tdks%rprimd,strsxc,      &
           & tdks%ucvol,tdks%usexcnhat,tdks%vhartr,tdks%vpsp,tdks%vtrial,tdks%vxc,  &
           & vxcavg,tdks%wvl,tdks%xccc3d,tdks%xred,taur=tdks%taur,                  &
           & vxc_hybcomp=tdks%vxc_hybcomp,vxctau=tdks%vxctau,add_tfw=tfw_activated, &
           & xcctau3d=tdks%xcctau3d,calc_ewald=calc_ewald)
 ABI_FREE(grchempottn)
 ABI_FREE(grewtn)
 ABI_FREE(kxc)

 ! set the zero of the potentials here
 if(dtset%usepotzero==2) tdks%vpsp(:) = tdks%vpsp(:) + tdks%ecore / ( tdks%zion * tdks%ucvol )

 !** Update PAW quantities
 !** Compute energies and potentials in the augmentation regions (spheres)
 !** and pseudopotential strengths (Dij quantities)
 if (psps%usepaw==1)then
   !** Local exact exch.: impose occ. matrix if required
   if (dtset%useexexch/=0) then
      if (xmpi_paral==1.and.mpi_enreg%paral_hf==1) then
         comm=mpi_enreg%comm_kpt
      else
         comm=mpi_enreg%comm_cell
      end if
      istep_mix=1; initialized0=0
      call setrhoijpbe0(dtset,initialized0,istep,istep_mix, &
                     & comm,my_natom,dtset%natom,dtset%ntypat,tdks%pawrhoij,tdks%pawtab, &
                     & dtset%typat,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   end if
   
   !** Computation of on-site densities/potentials/energies
   !** Force the recomputation of on-site potentials and Dij
   call paw_an_reset_flags(tdks%paw_an)
   !FB: @MT Changed self_consistent to false here. Is this right?
   call paw_ij_reset_flags(tdks%paw_ij,self_consistent=.false.)
   option=0; compch_sph=-1.d5; nzlmopt=0
   call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,ipert,         &
                & dtset%ixc,my_natom,dtset%natom,dtset%nspden,psps%ntypat,  &
                & dtset%nucdipmom,nzlmopt,option,tdks%paw_an,tdks%paw_an,   &
                & tdks%paw_ij,tdks%pawang,dtset%pawprtvol,tdks%pawrad,      &
                & tdks%pawrhoij,dtset%pawspnorb,tdks%pawtab,dtset%pawxcdev, &
                & dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,tdks%ucvol, &
                & psps%znuclpsp,comm_atom=mpi_enreg%comm_atom,              &
                & mpi_atmtab=mpi_enreg%my_atmtab,vpotzero=vpotzero)
   !Correct the average potential with the calculated constant vpotzero
   !Correct the total energies accordingly
   !vpotzero(1) = -beta/ucvol
   !vpotzero(2) = -1/ucvol sum_ij rho_ij gamma_ij
   write(msg,'(a,f14.6,2x,f14.6)') &
   & ' average electrostatic smooth potential [Ha] , [eV]', &
   & SUM(vpotzero(:)),SUM(vpotzero(:))*Ha_eV
   call wrtout(std_out,msg,'COLL')
   tdks%vtrial(:,:)=tdks%vtrial(:,:)+SUM(vpotzero(:))
   if(option/=1)then
      !Fix the direct total energy (non-zero only for charged systems)
      energies%e_paw=energies%e_paw-SUM(vpotzero(:))*dtset%cellcharge(1)
      !Fix the double counting total energy accordingly (for both charged AND
      !neutral systems)
      energies%e_pawdc=energies%e_pawdc-SUM(vpotzero(:))*tdks%zion+ &
                          & vpotzero(2)*dtset%cellcharge(1)
   end if
   
   !** Dij computation
   !FB: @MT fatvshift?
   nfftotf=tdks%pawfgr%ngfft(1)*tdks%pawfgr%ngfft(2)*tdks%pawfgr%ngfft(3)
   call pawdij(cplex,dtset%enunit,tdks%gprimd,ipert,my_natom,dtset%natom,            &
             & tdks%pawfgr%nfft,nfftotf,dtset%nspden,psps%ntypat,tdks%paw_an,        &
             & tdks%paw_ij,tdks%pawang,tdks%pawfgrtab,dtset%pawprtvol,tdks%pawrad,   &
             & tdks%pawrhoij,dtset%pawspnorb,tdks%pawtab,dtset%pawxcdev,k0,          &
             & dtset%spnorbscl,tdks%ucvol,dtset%cellcharge(1),tdks%vtrial,           &
             & tdks%vxc,tdks%xred,natvshift=dtset%natvshift,atvshift=dtset%atvshift, &
             & fatvshift=one,comm_atom=mpi_enreg%comm_atom,                          &
             & mpi_atmtab=mpi_enreg%my_atmtab,mpi_comm_grid=mpi_enreg%comm_fft,      &
             & nucdipmom=dtset%nucdipmom)
   
   !Symetrize Dij
   call symdij(tdks%gprimd,tdks%indsym,ipert,my_natom,dtset%natom,dtset%nsym, &
             & psps%ntypat,0,tdks%paw_ij,tdks%pawang,dtset%pawprtvol,         &
             & tdks%pawtab,tdks%rprimd,dtset%symafm,tdks%symrec,              &
             & comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 end if

 !** Initialize most of the Hamiltonian
 !** Allocate all arrays and initialize quantities that do not depend on k and spin.
 !FB: Should recompute cprj if ions have moved right?
 usecprj_local=0; if (psps%usepaw==1 .and. dtset%ionmov==0) usecprj_local=1
 call init_hamiltonian(gs_hamk,psps,tdks%pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,dtset%typat,    &
                     & tdks%xred,dtset%nfft,dtset%mgfft,dtset%ngfft,tdks%rprimd,dtset%nloalg,paw_ij=tdks%paw_ij,    &
                     & ph1d=tdks%ph1d,usecprj=usecprj_local,comm_atom=mpi_enreg%comm_atom,                          &
                     & mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,nucdipmom=dtset%nucdipmom, &
                     & use_gpu_cuda=dtset%use_gpu_cuda)

end subroutine rttddft_init_hamiltonian
 
!!****f* m_rttddft/rttddft_calc_density
!!
!! NAME
!!  rttddft_calc_density
!!
!! FUNCTION
!!  Compute electronic density (in 1/bohr^3) from the WF (cg coefficients)
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  tdks <type(tdks_type)> = Main RT-TDDFT object
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_driver/rttddft
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_calc_density(dtset, mpi_enreg, psps, tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 type(tdks_type),            intent(inout) :: tdks
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps

 !Local variables-------------------------------
 !FB-test integer :: i, me !FB-test
 !scalars
 integer, parameter          :: cplex=1
 integer                     :: cplex_rhoij
 integer                     :: ipert, idir
 integer                     :: my_natom
 integer                     :: nspden_rhoij
 integer                     :: tim_mkrho
 real(dp)                    :: compch_fft
 !arrays
 real(dp)                    :: qpt(3)
 real(dp),allocatable        :: rhowfg(:,:), rhowfr(:,:)
 type(pawrhoij_type),pointer :: pawrhoij_unsym(:)

! ***********************************************************************

 my_natom=mpi_enreg%my_natom

 tim_mkrho=1

 if (psps%usepaw==1) then

   ABI_MALLOC(rhowfg,(2,dtset%nfft))
   ABI_MALLOC(rhowfr,(dtset%nfft,dtset%nspden))

   ! 1-Compute density from WFs (without compensation charge density nhat)
   call mkrho(tdks%cg,dtset,tdks%gprimd,tdks%irrzon,tdks%kg,tdks%mcg,mpi_enreg, &
            & tdks%npwarr,tdks%occ0,tdks%paw_dmft,tdks%phnons,rhowfg,rhowfr,     &
            & tdks%rprimd,tim_mkrho,tdks%ucvol,tdks%wvl%den,tdks%wvl%wfs)

   ! 2-Compute cprj = <\psi_{n,k}|p_{i,j}>
   call ctocprj(tdks%atindx,tdks%cg,1,tdks%cprj,tdks%gmet,tdks%gprimd,0,0,0,           &
              & dtset%istwfk,tdks%kg,dtset%kptns,tdks%mcg,tdks%mcprj,dtset%mgfft,      &
              & dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,tdks%nattyp,   &
              & dtset%nband,dtset%natom,dtset%ngfft,dtset%nkpt,dtset%nloalg,           &
              & tdks%npwarr,dtset%nspinor,dtset%nsppol,psps%ntypat,dtset%paral_kgb,    &
              & tdks%ph1d,psps,tdks%rmet,dtset%typat,tdks%ucvol,tdks%unpaw,tdks%xred,  &
              & tdks%ylm,tdks%ylmgr)

   !FB-test
   !FB-test me=xmpi_comm_rank(mpi_enreg%comm_band)
   !FB-test do i=1, size(tdks%cprj(1,:))
   !FB-test    write(100+me,*) tdks%cprj(1,i)%cp
   !FB-test end do
   !FB-test write(100+me,*) ' '


   !paral atom
   if (my_natom/=dtset%natom) then
     ABI_MALLOC(pawrhoij_unsym,(dtset%natom))
     call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij, &
                             & nspden=dtset%nspden,spnorb=dtset%pawspnorb,        &
                             & cpxocc=dtset%pawcpxocc)
     call pawrhoij_alloc(pawrhoij_unsym,cplex_rhoij,nspden_rhoij,dtset%nspinor, &
                       & dtset%nsppol,dtset%typat,pawtab=tdks%pawtab,use_rhoijp=0)
   else
       pawrhoij_unsym => tdks%pawrhoij
   end if

   ! 3-Compute pawrhoij = \rho_{i,j} = \sum_{n,k}f_{n,k} \tilde{c}^{i,*}_{n,k} \tilde{c}^{j}_{n,k}
   call pawmkrhoij(tdks%atindx,tdks%atindx1,tdks%cprj,tdks%dimcprj,dtset%istwfk,    &
                 & dtset%kptopt,dtset%mband,tdks%mband_cprj,tdks%mcprj,dtset%mkmem, &
                 & mpi_enreg,dtset%natom,dtset%nband,dtset%nkpt,dtset%nspinor,      &
                 & dtset%nsppol,tdks%occ0,dtset%paral_kgb,tdks%paw_dmft,             &
                 & pawrhoij_unsym,tdks%unpaw,dtset%usewvl,dtset%wtk)

   ! 4-Symetrize rhoij, compute nhat and add it to rhor
   ! Note pawrhoij_unsym and pawrhoij are the same, which means that pawrhoij
   ! cannot be distributed over different atomic sites.
   ipert=0; idir=0; qpt(:)=zero; compch_fft=-1e-5_dp
   tdks%nhat = zero
   call pawmkrho(1,compch_fft,cplex,tdks%gprimd,idir,tdks%indsym,ipert,mpi_enreg, &
               & my_natom,dtset%natom,dtset%nspden,dtset%nsym,dtset%ntypat,       &
               & dtset%paral_kgb,tdks%pawang,tdks%pawfgr,tdks%pawfgrtab,          &
               & dtset%pawprtvol,tdks%pawrhoij,pawrhoij_unsym,tdks%pawtab,qpt,    &
               & rhowfg,rhowfr,tdks%rhor,tdks%rprimd,dtset%symafm,tdks%symrec,    &
               & dtset%typat,tdks%ucvol,dtset%usewvl,tdks%xred,pawnhat=tdks%nhat, &
               & rhog=tdks%rhog)

   ! 6-Take care of kinetic energy density
   if(dtset%usekden==1)then
     call mkrho(tdks%cg,dtset,tdks%gprimd,tdks%irrzon,tdks%kg,tdks%mcg,mpi_enreg, &
              & tdks%npwarr,tdks%occ0,tdks%paw_dmft,tdks%phnons,rhowfg,rhowfr,     &
              & tdks%rprimd,tim_mkrho,tdks%ucvol,tdks%wvl%den,tdks%wvl%wfs,option=1)
     !FB: Useful?
     call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,tdks%pawfgr, &
                  & rhowfg,tdks%taug,rhowfr,tdks%taur)
   end if

   ABI_FREE(rhowfg)
   ABI_FREE(rhowfr)

   if (my_natom/=dtset%natom) then
     call pawrhoij_free(pawrhoij_unsym)
     ABI_FREE(pawrhoij_unsym)
   else
      pawrhoij_unsym => NULL()
   end if

 else

   ! 1-Compute density from WFs
   call mkrho(tdks%cg,dtset,tdks%gprimd,tdks%irrzon,tdks%kg,tdks%mcg,mpi_enreg,   &
            & tdks%npwarr,tdks%occ0,tdks%paw_dmft,tdks%phnons,tdks%rhog,tdks%rhor, &
            & tdks%rprimd,tim_mkrho,tdks%ucvol,tdks%wvl%den,tdks%wvl%wfs)
   ! 2-Take care of kinetic energy density
   if(dtset%usekden==1)then
     call mkrho(tdks%cg,dtset,tdks%gprimd,tdks%irrzon,tdks%kg,tdks%mcg,mpi_enreg,   &
              & tdks%npwarr,tdks%occ0,tdks%paw_dmft,tdks%phnons,tdks%taug,tdks%taur, &
              & tdks%rprimd,tim_mkrho,tdks%ucvol,tdks%wvl%den,tdks%wvl%wfs,option=1)
   end if

 endif

 end subroutine rttddft_calc_density

!!****f* m_rttddft/rttddft_calc_energy
!!
!! NAME
!!  rttddft_calc_energy
!!
!! FUNCTION
!!  Computes total energy
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  energies <energies_type> = contains the different contribution ot the total energy
!!
!! OUTPUT
!!  etotal <real(dp)> = the total energy
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_driver/rttddft
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_calc_etot(dtset, energies, etotal)

 implicit none

 !Arguments ------------------------------------
 !scalars
 real(dp),                   intent(out)   :: etotal
 type(dataset_type),         intent(inout) :: dtset
 type(energies_type),        intent(inout) :: energies

 !Local variables-------------------------------
 !scalars
 !arrays
! ***********************************************************************

!  When the finite-temperature VG broadening scheme is used,
!  the total entropy contribution "tsmear*entropy" has a meaning,
!  and gather the two last terms of Eq.8 of VG paper
!  Warning : might have to be changed for fixed moment calculations
 if(dtset%occopt>=3 .and. dtset%occopt<=8) then
   if (abs(dtset%tphysel) < tol10) then
      energies%e_entropy = - dtset%tsmear * energies%entropy
   else
      energies%e_entropy = - dtset%tphysel * energies%entropy
   end if
 else
   energies%e_entropy = zero
 end if

 etotal = energies%e_kinetic     &
      & + energies%e_hartree     &  
      & + energies%e_xc          &
      & + energies%e_localpsp    &   
      & + energies%e_corepsp     &
      & + energies%e_entropy     &
      & + energies%e_ewald       &
      & + energies%e_vdw_dftd    &
      & + energies%e_nlpsp_vfock &
      & + energies%e_paw         
!FB: @MT Should one add the last e_paw contribution or not?
!     & + energies%e_chempot     &
!     & + energies%e_elecfield   &  
!     & + energies%e_magfield    &
!     & + energies%e_nucdip      &
!     & + energies%e_hybcomp_E0  &
!     & - energies%e_hybcomp_v0  &
!     & + energies%e_hybcomp_v   &
!     & + energies%e_constrained_dft

!if (psps%usepaw==0) etotal = etotal + energies%e_nlpsp_vfock - energies%e_fock0
!if (psps%usepaw==1) etotal = etotal + energies%e_paw + energies%e_fock

 end subroutine rttddft_calc_etot

!!****f* m_rttddft/rttddft_calc_eig
!!
!! NAME
!!  rttddft_calc_eig
!!
!! FUNCTION
!!  Computes eigenvalues from cg and ghc = <G|H|C> 
!!  and gsc = <G|S|C> if paw
!!
!! INPUTS
!!  cg <real(2,npw*nspinor*nband)> = the wavefunction coefficients
!!  ghc <real(2,npw*nspinor*nband)> = <G|H|C>
!!  istwfk <integer> = option describing the storage of wfs at k
!!  nband <integer> = number of bands
!!  npw <integer> = number of plane waves
!!  nspinor <integer> = dimension of spinors
!!  me_g0 <integer> = if set to 1 current proc contains G(0,0,0)
!!  comm <integer> = MPI communicator
!!  gsc <real(2,npw*nspinor*nband)> = <G|S|C> (optional - only in PAW)
!!
!! OUTPUT
!!  eig <real(nband)> = the eigenvalues
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_driver/rttddft
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_calc_eig(cg,eig,ghc,istwf_k,nband,npw,nspinor,me_g0,comm,gsc)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,  intent(in)            :: istwf_k
 integer,  intent(in)            :: nband
 integer,  intent(in)            :: npw
 integer,  intent(in)            :: nspinor
 integer,  intent(in)            :: me_g0
 integer,  intent(in)            :: comm
 !arrays
 real(dp), intent(in)            :: cg(2,npw*nspinor*nband)
 real(dp), intent(out)           :: eig(nband)
 real(dp), intent(in)            :: ghc(2,npw*nspinor*nband)
 real(dp), intent(in), optional  :: gsc(2,npw*nspinor*nband)

 !Local variables-------------------------------
 !scalars
 integer   :: iband
 integer   :: shift
 real(dp)  :: dprod_r, dprod_i
 !arrays
! ***********************************************************************

 do iband=1, nband
    shift = npw*nspinor*(iband-1)
    !Compute eigenvalues
    call dotprod_g(eig(iband),dprod_i,istwf_k,npw*nspinor,1, &
                 & ghc(:, shift+1:shift+npw*nspinor),        &
                 & cg(:, shift+1:shift+npw*nspinor),         &
                 & me_g0, comm)
    if (present(gsc)) then
       call dotprod_g(dprod_r,dprod_i,istwf_k,npw*nspinor,1, &
                    & gsc(:, shift+1:shift+npw*nspinor),     &
                    & cg(:, shift+1:shift+npw*nspinor),      &
                    & me_g0, comm)
       eig(iband) = eig(iband)/dprod_r
    end if
 end do

 end subroutine rttddft_calc_eig

!!****f* m_rttddft/rttddft_calc_enl
!!
!! NAME
!!  rttddft_calc_enl
!!
!! FUNCTION
!!  Computes the NL part of energy in NC case
!!
!! INPUTS
!!  cg <real(2,npw*nspinor*nband)> = the wavefunction coefficients
!!  ham_k <gs_hamiltonian_type> = hamiltonian at point k
!!  nband <integer> = number of bands
!!  npw <integer> = number of plane waves
!!  nspinor <integer> = dimension of spinors
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!
!! OUTPUT
!!  enl <real(nband)> = the non local part of the energy in NC case
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_driver/rttddft
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_calc_enl(cg,enl,ham_k,nband,npw,nspinor,mpi_enreg)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                   intent(in)    :: nband
 integer,                   intent(in)    :: npw
 integer,                   intent(in)    :: nspinor
 type(gs_hamiltonian_type), intent(in)    :: ham_k
 type(MPI_type),            intent(in)    :: mpi_enreg
 !arrays
 real(dp),                  intent(inout) :: cg(2,npw*nspinor*nband)
 real(dp),                  intent(out)   :: enl(nband)

 !Local variables-------------------------------
 !scalars
 integer, parameter              :: choice=1
 integer, parameter              :: cpopt=-1
 integer, parameter              :: paw_opt=0
 integer, parameter              :: signs=1
 integer, parameter              :: tim_getghc = 5
 !arrays
 type(pawcprj_type), allocatable :: cprj_dummy(:,:)
 real(dp),           allocatable :: eig_dummy(:)
 real(dp),           allocatable :: gvnlxc_dummy(:,:)
 real(dp),           allocatable :: gsc_dummy(:,:)
! ***********************************************************************

 ABI_MALLOC(cprj_dummy,(ham_k%natom,0))
 ABI_MALLOC(gsc_dummy,(0,0))
 ABI_MALLOC(gvnlxc_dummy, (0, 0))
 ABI_MALLOC(eig_dummy,(nband))
 call nonlop(choice,cpopt,cprj_dummy,enl,ham_k,0,eig_dummy,mpi_enreg,nband, &
            & 1,paw_opt,signs,gsc_dummy,tim_getghc,cg,gvnlxc_dummy)
 ABI_FREE(cprj_dummy)
 ABI_FREE(gsc_dummy)
 ABI_FREE(eig_dummy)
 ABI_FREE(gvnlxc_dummy)

 end subroutine rttddft_calc_enl

!!****f* m_rttddft/rttddft_calc_occ
!!
!! NAME
!!  rttddft_calc_occ
!!
!! FUNCTION
!!  Computes occupations at time t from cg(t), cg0 and occ0
!!  f_{n,k}(t) = \sum_{m} f_{m,k}(0) <\phi_m(0)|\phi_n(t)>   in NC
!!  f_{n,k}(t) = \sum_{m} f_{m,k}(0) <\phi_m(0)|S|\phi_n(t)> in PAW
!!
!! INPUTS
!!  cg <real(2,npw*nspinor*nband)> = the wavefunction coefficients
!!  cg0 <real(2,npw*nspinor*nband)> = the initial wavefunction coefficients
!!  ham_k <gs_hamiltonian_type> = hamiltonian at point k
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  nband <integer> = number of bands
!!  npw <integer> = number of plane waves
!!  nspinor <integer> = dimension of spinors
!!  occ0 <real(nband)> = initial occupations
!!
!! OUTPUT
!!  occ <real(nband)> = the occupations at time t
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_calc_occ(cg,cg0,cprj,cprj0,ham_k,mpi_enreg,nband,npw,nspinor,occ,occ0)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                   intent(in)    :: nband
 integer,                   intent(in)    :: npw
 integer,                   intent(in)    :: nspinor
 type(pawcprj_type),        intent(inout) :: cprj(:,:)
 type(pawcprj_type),        intent(inout) :: cprj0(:,:)
 type(gs_hamiltonian_type), intent(inout) :: ham_k
 type(MPI_type),            intent(inout) :: mpi_enreg
 !arrays
 real(dp),                  intent(inout) :: cg(2,npw*nspinor*nband)
 real(dp),                  intent(inout) :: cg0(2,npw*nspinor*nband)
 real(dp),                  intent(out)   :: occ(nband)
 real(dp),                  intent(in)    :: occ0(nband)

 !Local variables-------------------------------
 !scalars
 integer, parameter :: cpopt = 2 !cprojs already in memory (in cprj) 
 integer            :: iband, jband, ierr, ind
 integer            :: shift
 integer            :: tim_getcsc=6
 !arrays
 real(dp) :: csc(2*nband*nband)
! ***********************************************************************

 do iband = 1, nband
   !** 1 - First calc csc = <cg0|S|cg> (PAW) or <cg0|cg> (NC)
   shift = npw*nspinor*(iband-1)
   print*, 'FB-test: size(cprj,2), nspinor, nband, nspinor*nband', size(cprj,2), nspinor, nband, nspinor*nband
   print*, 'FB-test: size(cprj0,2), nspinor, nband, nspinor*nband', size(cprj0,2), nspinor, nband, nspinor*nband
   call getcsc(csc,cpopt,cg(:,shift+1:shift+npw*nspinor),cg0,cprj(:,iband:iband+(nspinor-1)),cprj0,ham_k,mpi_enreg,nband,tim_getcsc)
   !** 2 - If band parallel then reduce csc
   if (mpi_enreg%nproc_band > 1) then
      call xmpi_sum(csc,mpi_enreg%comm_band,ierr)
   end if
   !** 3 - Calc occupations from csc
   do jband = 1, nband
      ind = 2*jband
      occ(iband) = occ(iband) + occ0(jband)*(csc(ind-1)**2+csc(ind)**2)
   end do
 end do

 end subroutine rttddft_calc_occ

!!!****f* m_rttddft/rttddft_calc_occ
!!!
!!! NAME
!!!  rttddft_calc_occ
!!!
!!! FUNCTION
!!!  Compute occupation numbers at time t
!!!
!!! INPUTS
!!!  dtset <type(dataset_type)> = all input variables for this dataset
!!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!!  tdks <type(tdks_type)> = Main RT-TDDFT object
!!!
!!! OUTPUT
!!!
!!! SIDE EFFECTS
!!!
!!! PARENTS
!!!  m_rttddft_driver/rttddft
!!!
!!! CHILDREN
!!!
!!! SOURCE
!subroutine rttddft_calc_occ(tdks, dtset, mpi_enreg)
!
! implicit none
!
! !Arguments ------------------------------------
! !scalars
! type(dataset_type), intent(inout) :: dtset
! type(MPI_type),     intent(inout) :: mpi_enreg
! type(tdks_type),    intent(inout) :: tdks
!
! !Local variables-------------------------------
! integer :: ii
! real(dp) :: norm
! !scalars
! integer   :: bdtot_index
! integer   :: iband, jband
! integer   :: icg
! integer   :: ikpt, ikpt_loc
! integer   :: istwf_k
! integer   :: isppol
! integer   :: me_distrb
! integer   :: my_ikpt, my_nspinor
! integer   :: nband_k
! integer   :: npw_k
! integer   :: nkpt, nsppol
! integer   :: shift_i, shift_j
! integer   :: spaceComm_distrb
! real(dp)  :: dprod_r, dprod_i, dot_prod
! !arrays
!! ***********************************************************************
!
! !Init MPI
! spaceComm_distrb=mpi_enreg%comm_cell
! if (mpi_enreg%paral_kgb==1) spaceComm_distrb=mpi_enreg%comm_kpt
! me_distrb=xmpi_comm_rank(spaceComm_distrb)
!
! my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
!
! nkpt = dtset%nkpt
! nsppol = dtset%nsppol
!
! tdks%occ = zero
!
! icg=0
! bdtot_index=0
!
! !*** LOOP OVER SPINS
! do isppol=1,nsppol
!
!   ikpt_loc=0
!   !ikg=0
!
!   !*** BIG FAT k POINT LOOP
!   ikpt = 0
!   do while (ikpt_loc < nkpt)
!
!      ikpt_loc = ikpt_loc + 1
!      ikpt = ikpt_loc
!      my_ikpt = mpi_enreg%my_kpttab(ikpt)
!
!      nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
!      istwf_k=dtset%istwfk(ikpt)
!      npw_k=tdks%npwarr(ikpt)
!
!      if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) then
!         bdtot_index=bdtot_index+nband_k
!         cycle
!      end if
!
!      !Calc occupation number at point k
!      !PAW - calc <\phi_n|S|\phi_m> (getgsc)
!      if (dtset%usepaw /=0) then 
!         !do iband = 1, nband_k
!         !   shift_i = icg+npw_k*my_nspinor*(iband-1)
!         !   call getcsc(csc,cpopt,cwavef,cwavef_left,cprj,cprj_left,tdks%gs_hamk,mpi_enreg,ndat,tim_getcsc,)
!         !end do
!      !Norm-conserving - calc <\phi_n|\phi_m> (dotprod_g)
!      else
!         do iband = 1, nband_k
!            shift_i = icg+npw_k*my_nspinor*(iband-1)
!            do jband = 1, nband_k
!               shift_j = icg+npw_k*my_nspinor*(jband-1)
!               call dotprod_g(dprod_r,dprod_i,istwf_k,npw_k*my_nspinor,2,tdks%cg0(:, shift_i+1:shift_i+npw_k*my_nspinor),&
!                            & tdks%cg(:, shift_j+1:shift_j+npw_k*my_nspinor),mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
!               dot_prod = dprod_r**2 + dprod_i**2
!               norm = 0.0_dp
!               do ii = shift_i+1, shift_i+npw_k*my_nspinor
!                  norm = norm + tdks%cg(1,ii)**2+tdks%cg(2,ii)**2
!               end do 
!               !FB-test print*, 'sum(cg(t)^2):', sqrt(norm)
!               norm = 0.0_dp
!               do ii = shift_i+1, shift_i+npw_k*my_nspinor
!                  norm = norm + tdks%cg0(1,ii)**2+tdks%cg0(2,ii)**2
!               end do 
!               !FB-test print*, 'sum(cg(0)^2):', sqrt(norm)
!               tdks%occ(bdtot_index+iband) = tdks%occ(bdtot_index+iband)+tdks%occ0(bdtot_index+jband)*dot_prod
!               !FB-test print*, 'ikpt, iband, jband, bdtot_index, dot_prod :', ikpt, iband, jband, bdtot_index, dot_prod
!               !FB-test print*, 'bdtot_index+iband, tdks%occ(bdtot_index+iband):', bdtot_index+iband, tdks%occ(bdtot_index+iband)
!               !FB-test print*, 'bdtot_index+jband, tdks%occ0(bdtot_index+jband), tdks%occ0(bdtot_index+jband)*dot_prod:', &
!               !FB-test        & bdtot_index+jband, tdks%occ0(bdtot_index+jband), tdks%occ0(bdtot_index+jband)*dot_prod
!            end do
!         end do
!      end if
!
!      !** Also shift array memory if dtset%mkmem/=0
!      if (dtset%mkmem/=0) then
!         icg=icg+npw_k*my_nspinor*nband_k
!         bdtot_index=bdtot_index+nband_k
!      end if
!
!   end do !nkpt
!
! end do !nsppol
!
! print*, 'FB-test - occ', tdks%occ
! print*, 'FB-test - occ0', tdks%occ0
!
! !FB TODO some MPI sum are needed here
!
! end subroutine rttddft_calc_occ

end module m_rttddft
!!***
