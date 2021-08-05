!!****m* ABINIT/m_rttddft
!! NAME
!!  m_rttddft
!!
!! FUNCTION
!!  Contains various subroutines used in RT-TDDFT
!!
!! COPYRIGHT
!!  Copyright (C) 2021 ABINIT group (FB, MT)
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

module m_rttddft

 use defs_basis
 use defs_abitypes,     only: MPI_type
 use defs_datatypes,    only: pseudopotential_type
 
 use m_cgprj,            only: ctocprj
 use m_dtfil,            only: datafiles_type
 use m_dtset,            only: dataset_type
 use m_fourier_interpol, only: transgrid
 use m_kg,               only: getcut, getph
 use m_mkrho,            only: mkrho
 use m_paw_an,           only: paw_an_reset_flags
 use m_paw_correlations, only: setrhoijpbe0
 use m_paw_denpot,       only: pawdenpot
 use m_pawdij,           only: pawdij, symdij
 use m_paw_ij,           only: paw_ij_reset_flags
 use m_paw_mkrho,        only: pawmkrho
 use m_paw_nhat,         only: nhatgrid
 use m_paw_occupancies,  only: pawmkrhoij
 use m_pawrhoij,         only: pawrhoij_type,  pawrhoij_free, &
                               pawrhoij_alloc, pawrhoij_inquire_dim
 use m_paw_tools,        only: chkpawovlp
 use m_rttddft_types,    only: tdks_type 
 use m_specialmsg,       only: wrtout
 use m_setvtr,           only: setvtr
 use m_xmpi,             only: xmpi_paral

 implicit none

 private
!!***

 public :: calc_density
 public :: update_after_nuc_step
 public :: update_paw

contains 

!!****f* m_rttddft/calc_density
!!
!! NAME
!!  tdks_init
!!
!! FUNCTION
!!  Compute electronic density (in 1/bohr^3) from the WF (cg coefficients)
!!
!! INPUTS
!!  tdks <class(tdks_type)> = the tdks object to initialize
!!  dtfil <type datafiles_type> = infos about file names, file unit numbers
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_driver
!!
!! CHILDREN
!!
!! SOURCE
subroutine calc_density(tdks, dtfil, dtset, mpi_enreg, psps)

 implicit none

 !Arguments ------------------------------------
 !scalars
 class(tdks_type),           intent(inout) :: tdks
 type(datafiles_type),       intent(in)    :: dtfil
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 
 !Local variables-------------------------------
 !scalars
 integer                     :: cplex, cplex_rhoij
 integer                     :: ipert, idir
 integer                     :: my_natom
 integer                     :: nspden_rhoij
 integer                     :: tim_mkrho
 real(dp)                    :: compch_fft
 !arrays
 real(dp)                    :: gmet(3,3) 
 real(dp),allocatable        :: nhat(:,:)
 real(dp),allocatable        :: ph1d(:,:)
 real(dp)                    :: qpt(3)
 real(dp)                    :: rmet(3,3) 
 real(dp),allocatable        :: rhowfg(:,:), rhowfr(:,:)
 type(pawrhoij_type),pointer :: pawrhoij_unsym(:)
 
! ***********************************************************************

 my_natom=mpi_enreg%my_natom 

 tim_mkrho=1

 if (psps%usepaw==1) then

   ABI_MALLOC(rhowfg,(2,dtset%nfft))
   ABI_MALLOC(rhowfr,(dtset%nfft,dtset%nspden))
   ABI_MALLOC(ph1d,(2,3*(2*tdks%pawfgr%mgfft+1)*dtset%natom))
   ABI_MALLOC(nhat,(tdks%nfftf,dtset%nspden*psps%usepaw))

   ! 1-Compute structure factor phases for current atomic pos
   call getph(tdks%atindx,dtset%natom,tdks%pawfgr%ngfft(1),tdks%pawfgr%ngfft(2), &
            & tdks%pawfgr%ngfft(3),ph1d,tdks%xred)

   ! 2-Compute density from WFs (without compensation charge density nhat)
   call mkrho(tdks%cg,dtset,tdks%gprimd,tdks%irrzon,tdks%kg,tdks%mcg,mpi_enreg, &
            & tdks%npwarr,tdks%occ,tdks%paw_dmft,tdks%phnons,rhowfg,rhowfr, &
            & tdks%rprimd,tim_mkrho,tdks%ucvol,tdks%wvl%den,tdks%wvl%wfs)
   ! transfer density from the coarse to the fine FFT grid
   call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,tdks%pawfgr, &
                & rhowfg,tdks%rhog,rhowfr,tdks%rhor)

   ! 3-Compute cprj = <\psi_{n,k}|p_{i,j}>
   call ctocprj(tdks%atindx,tdks%cg,1,tdks%cprj,gmet,tdks%gprimd,0,0,0, &
              & dtset%istwfk,tdks%kg,dtset%kptns,tdks%mcg,tdks%mcprj,dtset%mgfft, &
              & dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,tdks%nattyp, &
              & dtset%nband,dtset%natom,dtset%ngfft,dtset%nkpt,dtset%nloalg, &
              & tdks%npwarr,dtset%nspinor,dtset%nsppol,psps%ntypat,dtset%paral_kgb, &
              & ph1d,psps,rmet,dtset%typat,tdks%ucvol,dtfil%unpaw,tdks%xred, &
              & tdks%ylm,tdks%ylmgr)

   !paral atom
   if (my_natom/=dtset%natom) then
     ABI_MALLOC(pawrhoij_unsym,(dtset%natom))
     call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij, &
                             & nspden=dtset%nspden,spnorb=dtset%pawspnorb, &
                             & cpxocc=dtset%pawcpxocc)
     call pawrhoij_alloc(pawrhoij_unsym,cplex_rhoij,nspden_rhoij,dtset%nspinor, &
                       & dtset%nsppol,dtset%typat,pawtab=tdks%pawtab,use_rhoijp=0)
   else
       pawrhoij_unsym => tdks%pawrhoij
   end if

   ! 4-Compute pawrhoij = \rho_{i,j} = \sum_{n,k}f_{n,k} \tilde{c}^{i,*}_{n,k} \tilde{c}^{j}_{n,k} 
   call pawmkrhoij(tdks%atindx,tdks%atindx1,tdks%cprj,tdks%dimcprj,dtset%istwfk, &
                 & dtset%kptopt,dtset%mband,tdks%mband_cprj,tdks%mcprj,dtset%mkmem, &
                 & mpi_enreg,dtset%natom,dtset%nband,dtset%nkpt,dtset%nspinor, &
                 & dtset%nsppol,tdks%occ,dtset%paral_kgb,tdks%paw_dmft, &
                 & pawrhoij_unsym,dtfil%unpaw,dtset%usewvl,dtset%wtk)

   ! 5-Symetrize rhoij, compute nhat and add it to rhor
   ! Note pawrhoij_unsym and pawrhoij are the same, which means that pawrhoij
   ! cannot be distributed over different atomic sites.
   cplex=1; ipert=0; idir=0; qpt(:)=zero; compch_fft=-1e-5_dp
   nhat = zero
   call pawmkrho(1,compch_fft,cplex,tdks%gprimd,idir,tdks%indsym,ipert,mpi_enreg, &
               & my_natom,dtset%natom,dtset%nspden,dtset%nsym,dtset%ntypat, &
               & dtset%paral_kgb,tdks%pawang,tdks%pawfgr,tdks%pawfgrtab, &
               & dtset%pawprtvol,tdks%pawrhoij,pawrhoij_unsym,tdks%pawtab,qpt, &
               & rhowfg,rhowfr,tdks%rhor,tdks%rprimd,dtset%symafm,tdks%symrec, &
               & dtset%typat,tdks%ucvol,dtset%usewvl,tdks%xred,pawnhat=nhat, &
               & rhog=tdks%rhog)

   ! 6-Take care of kinetic energy density
   if(dtset%usekden==1)then
     call mkrho(tdks%cg,dtset,tdks%gprimd,tdks%irrzon,tdks%kg,tdks%mcg,mpi_enreg, &
              & tdks%npwarr,tdks%occ,tdks%paw_dmft,tdks%phnons,rhowfg,rhowfr, &
              & tdks%rprimd,tim_mkrho,tdks%ucvol,tdks%wvl%den,tdks%wvl%wfs,option=1)
     call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,tdks%pawfgr, &
                  & rhowfg,tdks%taug,rhowfr,tdks%taur)
   end if

   ABI_FREE(rhowfg)
   ABI_FREE(rhowfr)
   ABI_FREE(ph1d)
   ABI_FREE(nhat)

   if (my_natom/=dtset%natom) then
     call pawrhoij_free(pawrhoij_unsym)
     ABI_FREE(pawrhoij_unsym)
   else
      pawrhoij_unsym => NULL()
   end if

 else

   ! 1-Compute density from WFs
   call mkrho(tdks%cg,dtset,tdks%gprimd,tdks%irrzon,tdks%kg,tdks%mcg,mpi_enreg, &
            & tdks%npwarr,tdks%occ,tdks%paw_dmft,tdks%phnons,tdks%rhog,tdks%rhor, &
            & tdks%rprimd,tim_mkrho,tdks%ucvol,tdks%wvl%den,tdks%wvl%wfs)
   ! 2-Take care of kinetic energy density
   if(dtset%usekden==1)then
     call mkrho(tdks%cg,dtset,tdks%gprimd,tdks%irrzon,tdks%kg,tdks%mcg,mpi_enreg, &
              & tdks%npwarr,tdks%occ,tdks%paw_dmft,tdks%phnons,tdks%taug,tdks%taur, &
              & tdks%rprimd,tim_mkrho,tdks%ucvol,tdks%wvl%den,tdks%wvl%wfs,option=1)
   end if

 endif

 end subroutine calc_density

!!****f* m_rttddft/update_after_nuc_step
!!
!! NAME
!!  tdks_init
!!
!! FUNCTION
!!  Update various quantities that needs to be recomputed/changed
!!  after a change of xred during the nuc propagation step
!!
!! INPUTS
!!  tdks <class(tdks_type)> = the tdks object to initialize
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_driver
!!
!! CHILDREN
!!
!! SOURCE
subroutine update_after_nuc_step(tdks, dtset, istep, mpi_enreg, psps)

 implicit none

 !Arguments ------------------------------------
 !scalars
 class(tdks_type),           intent(inout) :: tdks
 integer,                    intent(in)    :: istep
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 
 !Local variables-------------------------------
 !scalars
 logical              :: tfw_activated
 integer              :: forces_needed
 integer              :: moved_atm_inside, moved_rhor
 integer              :: my_natom
 integer              :: nkxc, n1xccc, n3xccc
 integer              :: optcut, optgr0, optgr1, optgr2, optrad, optene
 integer              :: stress_needed
 real(dp)             :: vxcavg
 !arrays
 real(dp),parameter   :: k0(3)=(/zero,zero,zero/)
 real(dp),allocatable :: grchempottn(:,:)
 real(dp),allocatable :: grewtn(:,:)
 real(dp),allocatable :: kxc(:,:)
 real(dp)             :: strsxc(6)
 
! ***********************************************************************

 my_natom=mpi_enreg%my_natom 

 !Compute large sphere G^2 cut-off (gsqcut) and box / sphere ratio
 !FB: Needed? Box didn't change only nuclear pos..
 if (psps%usepaw==1) then 
   call getcut(tdks%boxcut,dtset%pawecutdg,tdks%gmet,tdks%gsqcut,dtset%iboxcut, &
             & std_out,k0,tdks%pawfgr%ngfft)
 else
   call getcut(tdks%boxcut,dtset%ecut,tdks%gmet,tdks%gsqcut,dtset%iboxcut, &
             & std_out,k0,tdks%pawfgr%ngfft)
 end if
   
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
    call nhatgrid(tdks%atindx1,tdks%gmet,my_natom,dtset%natom,&
                & tdks%nattyp,tdks%pawfgr%ngfft,psps%ntypat,optcut,optgr0,optgr1, &
                & optgr2,optrad,tdks%pawfgrtab,tdks%pawtab,tdks%rprimd,dtset%typat, &
                & tdks%ucvol,tdks%xred,comm_atom=mpi_enreg%comm_atom, &
                & mpi_atmtab=mpi_enreg%my_atmtab,comm_fft=mpi_enreg%comm_fft, &
                & distribfft=mpi_enreg%distribfft)

 end if

 !!FB: Needed? If yes, then don't forget to put it back in tdks_init/second_setup as well
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
 !endif

 !** Set up the potential (calls setvtr)
 !    The following steps have been gathered in the setvtr routine:
 !    - get Ewald energy and Ewald forces
 !    - compute local ionic pseudopotential vpsp
 !    - possibly compute 3D core electron density xccc3d
 !    - possibly compute 3D core kinetic energy density
 !    - possibly compute vxc and vhartr
 !    - set up vtrial
 optene = 4; nkxc=0; moved_atm_inside=0; moved_rhor=0
 n1xccc=0;if (psps%n1xccc/=0) n1xccc=psps%n1xccc
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=tdks%pawfgr%nfft
 strsxc(:)=zero
 if (dtset%tfkinfunc==12) tfw_activated=.true.
 ABI_MALLOC(grchempottn,(3,dtset%natom))
 ABI_MALLOC(grewtn,(3,dtset%natom))
 ABI_MALLOC(kxc,(tdks%pawfgr%nfft,nkxc))
 call setvtr(tdks%atindx1,dtset,tdks%energies,tdks%gmet,tdks%gprimd,grchempottn,  &
          & grewtn,tdks%grvdw,tdks%gsqcut,istep,kxc,tdks%pawfgr%mgfft,            &
          & moved_atm_inside,moved_rhor,mpi_enreg,tdks%nattyp,tdks%pawfgr%nfft,   &
          & tdks%pawfgr%ngfft,tdks%ngrvdw,tdks%nhat,tdks%nhatgr,tdks%nhatgrdim,   &
          & nkxc,psps%ntypat,psps%n1xccc,n3xccc,optene,tdks%pawrad,tdks%pawtab,   &
          & tdks%ph1df,psps,tdks%rhog,tdks%rhor,tdks%rmet,tdks%rprimd,strsxc,     &
          & tdks%ucvol,tdks%usexcnhat,tdks%vhartr,tdks%vpsp,tdks%vtrial,tdks%vxc, &
          & vxcavg,tdks%wvl,tdks%xccc3d,tdks%xred,taur=tdks%taur,                 &
          & vxc_hybcomp=tdks%vxc_hybcomp,vxctau=tdks%vxctau,add_tfw=tfw_activated,&
          & xcctau3d=tdks%xcctau3d)

   ! set the zero of the potentials here
   if(dtset%usepotzero==2) tdks%vpsp(:) = tdks%vpsp(:) + tdks%ecore / ( tdks%zion * tdks%ucvol )

end subroutine update_after_nuc_step

!!****f* m_rttddft/update_before_ele_step
!!
!! NAME
!!  update_before_ele_step
!!
!! FUNCTION
!!  Update various PAW quantities before performing the RT-TDDFT
!!  propagation step of the KS orbitals. In particular computes
!!  energies and potentials in augmentation spheres ans pseudopotential
!!  strengths Dij
!!
!! INPUTS
!!  tdks <class(tdks_type)> = the tdks object to initialize
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  istep <integer> = step number
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!  m_rttddft_propagate
!!
!! CHILDREN
!!
!! SOURCE
subroutine update_paw(tdks, dtset, istep, mpi_enreg, psps)

 implicit none

 !Arguments ------------------------------------
 !scalars
 class(tdks_type),           intent(inout) :: tdks
 integer,                    intent(in)    :: istep
 type(dataset_type),         intent(inout) :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 
 !Local variables-------------------------------
 !scalars
 character(len=500)   :: msg
 integer              :: comm
 integer,parameter    :: cplex=1
 integer,parameter    :: ipert=0
 integer              :: initialized0
 integer              :: istep_mix
 integer              :: my_natom
 integer              :: nfftotf
 integer              :: nzlmopt
 integer              :: option
 real(dp)             :: compch_sph
 real(dp)             :: hyb_mixing,hyb_mixing_sr
 !arrays
 real(dp),parameter   :: k0(3)=(/zero,zero,zero/)
 real(dp)             :: vpotzero(2)
 
! ***********************************************************************

 my_natom=mpi_enreg%my_natom 

 !** Local exact exch.: impose occ. matrix if required
 if (dtset%useexexch/=0) then
    if (xmpi_paral==1.and.mpi_enreg%paral_hf==1) then
       comm=mpi_enreg%comm_kpt
    else
       comm=mpi_enreg%comm_cell
    end if
    !FB: I changed the values of istep_mix and initialized here for our purpose
    !FB: compared to GS. Hope this is correct..
    istep_mix=1; initialized0=0
    call setrhoijpbe0(dtset,initialized0,istep,istep_mix, &
                   & comm,my_natom,dtset%natom,dtset%ntypat,tdks%pawrhoij,tdks%pawtab, &
                  & dtset%typat,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 end if

 hyb_mixing=zero; hyb_mixing_sr=zero
 if (dtset%usefock==1) then
    hyb_mixing=tdks%fock%fock_common%hyb_mixing
    hyb_mixing_sr=tdks%fock%fock_common%hyb_mixing_sr
 endif

 !** Computation of on-site densities/potentials/energies
 !** Force the recomputation of on-site potentials and Dij
 call paw_an_reset_flags(tdks%paw_an)
 !FB: Changed self_consistent to false here. Is this right?
 call paw_ij_reset_flags(tdks%paw_ij,self_consistent=.false.) 
 !FB: Used option = 0 here so both potentials and energies are recomputed.
 !FB: Would potentials only be sufficient?
 option=0; compch_sph=-1.d5; nzlmopt=0
 vpotzero(:)=zero
 call pawdenpot(compch_sph,tdks%energies%e_paw,tdks%energies%e_pawdc,ipert, &
              & dtset%ixc,my_natom,dtset%natom,dtset%nspden,psps%ntypat,    &
              & dtset%nucdipmom,nzlmopt,option,tdks%paw_an,tdks%paw_an,     &
              & tdks%paw_ij,tdks%pawang,dtset%pawprtvol,tdks%pawrad,        &
              & tdks%pawrhoij,dtset%pawspnorb,tdks%pawtab,dtset%pawxcdev,   &
              & dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,tdks%ucvol,   &
              & psps%znuclpsp,comm_atom=mpi_enreg%comm_atom,                & 
              & mpi_atmtab=mpi_enreg%my_atmtab,hyb_mixing=hyb_mixing,       &
              & hyb_mixing_sr=hyb_mixing_sr,vpotzero=vpotzero)

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
    tdks%energies%e_paw=tdks%energies%e_paw-SUM(vpotzero(:))*dtset%cellcharge(1)
    !Fix the double counting total energy accordingly (for both charged AND
    !neutral systems)
    tdks%energies%e_pawdc=tdks%energies%e_pawdc-SUM(vpotzero(:))*tdks%zion+ &
                        & vpotzero(2)*dtset%cellcharge(1)
 end if

 !** Dij computation
 !FB: fatvshift?
 nfftotf=tdks%pawfgr%ngfft(1)*tdks%pawfgr%ngfft(2)*tdks%pawfgr%ngfft(3)
 call pawdij(cplex,dtset%enunit,tdks%gprimd,ipert,my_natom,dtset%natom,            &
           & tdks%pawfgr%nfft,nfftotf,dtset%nspden,psps%ntypat,tdks%paw_an,        &
           & tdks%paw_ij,tdks%pawang,tdks%pawfgrtab,dtset%pawprtvol,tdks%pawrad,   &
           & tdks%pawrhoij,dtset%pawspnorb,tdks%pawtab,dtset%pawxcdev,k0,          &
           & dtset%spnorbscl,tdks%ucvol,dtset%cellcharge(1),tdks%vtrial,           &
           & tdks%vxc,tdks%xred,natvshift=dtset%natvshift,atvshift=dtset%atvshift, &
           & fatvshift=one,comm_atom=mpi_enreg%comm_atom,                          &
           & mpi_atmtab=mpi_enreg%my_atmtab,mpi_comm_grid=mpi_enreg%comm_fft,      &
           & hyb_mixing=hyb_mixing,hyb_mixing_sr=hyb_mixing_sr,                    &
           & nucdipmom=dtset%nucdipmom)
        
 !Symetrize Dij
 call symdij(tdks%gprimd,tdks%indsym,ipert,my_natom,dtset%natom,dtset%nsym, &
           & psps%ntypat,0,tdks%paw_ij,tdks%pawang,dtset%pawprtvol,         &
           & tdks%pawtab,tdks%rprimd,dtset%symafm,tdks%symrec,              &
           & comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

end subroutine update_paw

end module m_rttddft
!!***
