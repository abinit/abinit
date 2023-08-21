!!****m* ABINIT/m_rttddft
!! NAME
!!  m_rttddft
!!
!! FUNCTION
!!  Contains various subroutines used in RT-TDDFT
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2022 ABINIT group (FB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

 use m_dtset,            only: dataset_type
 use m_energies,         only: energies_type
 use m_hamiltonian,      only: init_hamiltonian, gs_hamiltonian_type
 use m_kg,               only: getph
 use m_nonlop,           only: nonlop
 use m_paw_an,           only: paw_an_reset_flags
 use m_paw_correlations, only: setrhoijpbe0
 use m_paw_denpot,       only: pawdenpot
 use m_pawdij,           only: pawdij, symdij
 use m_paw_ij,           only: paw_ij_reset_flags
 use m_paw_nhat,         only: nhatgrid
 use m_paw_tools,        only: chkpawovlp
 use m_profiling_abi,    only: abimem_record
 use m_rttddft_tdks,     only: tdks_type
 use m_specialmsg,       only: wrtout
 use m_setvtr,           only: setvtr
 use m_xmpi,             only: xmpi_paral

 implicit none

 private
!!***

 public :: rttddft_setup_ele_step
 public :: rttddft_init_hamiltonian
!!***

contains
!!***

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
!!***

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
!!  energies <energies_type> = contains various contribution to the energy
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  tdks <type(tdks_type)> = Main RT-TDDFT object
!!
!! OUTPUT
!!  gs_hamk <type(gs_hamiltonian_type)> = Hamiltonian object
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
 !FB: @MT Are the values of moved_atm_inside and moved_rhor correct?
 optene = 4; nkxc=0; moved_atm_inside=0; moved_rhor=1
 n1xccc=0;if (psps%n1xccc/=0) n1xccc=psps%n1xccc
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=tdks%pawfgr%nfft
 strsxc(:)=zero
 !FB: tfw_activated is a save variable in scfcv, should maybe check where it appears again
 tfw_activated=.false.
 if (dtset%tfkinfunc==12) tfw_activated=.true.
 ABI_MALLOC(grchempottn,(3,dtset%natom))
 ABI_MALLOC(grewtn,(3,dtset%natom))
 ABI_MALLOC(kxc,(tdks%pawfgr%nfft,nkxc))
 calc_ewald = .false.
 if (dtset%ionmov/=0 .or. istep == tdks%first_step) calc_ewald=.true.
 !FB: Should we also add an option to avoid recomputing xccc3d as for Ewald?
 call setvtr(tdks%atindx1,dtset,energies,tdks%gmet,tdks%gprimd,grchempottn,         &
           & grewtn,tdks%grvdw,tdks%gsqcut,istep,kxc,tdks%pawfgr%mgfft,             &
           & moved_atm_inside,moved_rhor,mpi_enreg,tdks%nattyp,tdks%pawfgr%nfft,    &
           & tdks%pawfgr%ngfft,tdks%ngrvdw,tdks%nhat,tdks%nhatgr,tdks%nhatgrdim,    &
           & nkxc,psps%ntypat,n1xccc,n3xccc,optene,tdks%pawang,tdks%pawrad,         &
           & tdks%pawrhoij,tdks%pawtab,tdks%ph1df,psps,tdks%rhog,tdks%rhor,         &
           & tdks%rmet,tdks%rprimd,strsxc,tdks%ucvol,tdks%usexcnhat,tdks%vhartr,    &
           & tdks%vpsp,tdks%vtrial,tdks%vxc,vxcavg,tdks%wvl,tdks%xccc3d,tdks%xred,  &
           & taur=tdks%taur,vxc_hybcomp=tdks%vxc_hybcomp,vxctau=tdks%vxctau,        &
           & add_tfw=tfw_activated,xcctau3d=tdks%xcctau3d,calc_ewald=calc_ewald)
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
                & dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,            &
                & dtset%xc_taupos,tdks%ucvol,                               &
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
!!***
 
end module m_rttddft
!!***
