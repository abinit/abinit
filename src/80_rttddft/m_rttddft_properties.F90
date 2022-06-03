!!****m* ABINIT/m_rttddft_properties
!! NAME
!!  m_rttddft_properties
!!
!! FUNCTION
!!  Contains most of the subroutines to compute
!!  properties (energy, occupations, eigenvalues..)
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2022 ABINIT group (FB)
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

module m_rttddft_properties

 use defs_basis
 use defs_abitypes,      only: MPI_type
 use defs_datatypes,     only: pseudopotential_type

 use m_bandfft_kpt,      only: bandfft_kpt_type
 use m_cgprj,            only: ctocprj
 use m_cgtools,          only: dotprod_g
 use m_dtset,            only: dataset_type
 use m_energies,         only: energies_type
 use m_fourier_interpol, only: transgrid
 use m_hamiltonian,      only: gs_hamiltonian_type
 use m_mkrho,            only: mkrho
 use m_nonlop,           only: nonlop
 use m_pawcprj,          only: pawcprj_type, pawcprj_alloc, pawcprj_get, &
                             & pawcprj_free, pawcprj_mpi_allgather
 use m_paw_mkrho,        only: pawmkrho
 use m_paw_occupancies,  only: pawmkrhoij
 use m_pawrhoij,         only: pawrhoij_type, pawrhoij_free, &
                             & pawrhoij_alloc, pawrhoij_inquire_dim
 use m_profiling_abi,    only: abimem_record !FB: Needed to pass bbot test farm
 use m_rttddft_tdks,     only: tdks_type
 use m_spacepar,         only: meanvalue_g
 use m_xmpi,             only: xmpi_sum, xmpi_comm_rank

 implicit none

 private
!!***

 public :: rttddft_calc_density
 public :: rttddft_calc_etot
 public :: rttddft_calc_eig
 public :: rttddft_calc_enl
 public :: rttddft_calc_kin
 public :: rttddft_calc_occ
!!***

contains
!!***

!!****f* m_rttddft_properties/rttddft_calc_density
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
!! PARENTS
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
            & tdks%npwarr,tdks%occ0,tdks%paw_dmft,tdks%phnons,rhowfg,rhowfr,    &
            & tdks%rprimd,tim_mkrho,tdks%ucvol,tdks%wvl%den,tdks%wvl%wfs)

   ! 2-Compute cprj = <\psi_{n,k}|p_{i,j}>
   call ctocprj(tdks%atindx,tdks%cg,1,tdks%cprj,tdks%gmet,tdks%gprimd,0,0,0,           &
              & dtset%istwfk,tdks%kg,dtset%kptns,tdks%mcg,tdks%mcprj,dtset%mgfft,      &
              & dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,tdks%nattyp,   &
              & dtset%nband,dtset%natom,dtset%ngfft,dtset%nkpt,dtset%nloalg,           &
              & tdks%npwarr,dtset%nspinor,dtset%nsppol,psps%ntypat,dtset%paral_kgb,    &
              & tdks%ph1d,psps,tdks%rmet,dtset%typat,tdks%ucvol,tdks%unpaw,tdks%xred,  &
              & tdks%ylm,tdks%ylmgr)

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
                 & dtset%nsppol,tdks%occ0,dtset%paral_kgb,tdks%paw_dmft,            &
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
!!***

!!****f* m_rttddft_properties/rttddft_calc_energy
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
!!  occ <real(nband*nkpt*nsspol)> = occupation numbers at time t
!!
!! OUTPUT
!!  etotal <real(dp)> = the total energy
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_calc_etot(dtset, energies, etotal, occ)

 implicit none

 !Arguments ------------------------------------
 !scalars
 type(dataset_type),  intent(inout) :: dtset
 type(energies_type), intent(inout) :: energies
 real(dp),            intent(out)   :: etotal
 !arrays
 real(dp),            intent(in)    :: occ(:)

 !Local variables-------------------------------
 !scalars
 real(dp) :: entropy
 !arrays

! ***********************************************************************

 ! Compute electronic entropy
 call rttddft_calc_ent(entropy, dtset, occ)

!  When the finite-temperature VG broadening scheme is used,
!  the total entropy contribution "tsmear*entropy" has a meaning,
!  and gather the two last terms of Eq.8 of VG paper
!  Warning : might have to be changed for fixed moment calculations
 if(dtset%occopt>=3 .and. dtset%occopt<=8) then
   if (abs(dtset%tphysel) < tol10) then
      energies%e_entropy = - dtset%tsmear * entropy
   else
      energies%e_entropy = - dtset%tphysel * entropy
   end if
 else
    !FB - TODO: Might want to clean this ?!
   energies%e_entropy = -entropy
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
!FB: Seeems like all the other contributions are not relevant here @MT?
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
!!***

!!****f* m_rttddft_properties/rttddft_calc_eig
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
!!  istwf_k <integer> = option describing the storage of wfs at k
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
!! PARENTS
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
!!***

!!****f* m_rttddft_properties/rttddft_calc_kin
!!
!! NAME
!!  rttddft_calc_kin
!!
!! FUNCTION
!!  Computes the NL part of energy in NC case
!!
!! INPUTS
!!  cg <real(2,npw*nspinor*nband)> = the wavefunction coefficients
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  ham_k <gs_hamiltonian_type> = hamiltonian at point k
!!  nband <integer> = number of bands
!!  npw <integer> = number of plane waves
!!  nspinor <integer> = dimension of spinors
!!  occ0 <real(nband)> = initial occupations
!!  wk <real> = weight of associated kpt
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  bandfft_kpt <bandfft_kpt_type> = additional info on parallelisation
!!   for the associated kpt
!!
!! OUTPUT
!!  kin <real(nband)> = the non local part of the energy in NC case
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_calc_kin(kin,cg,dtset,ham_k,nband,npw,nspinor,occ0,wk,mpi_enreg,bandfft)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                   intent(in)    :: nband
 integer,                   intent(in)    :: npw
 integer,                   intent(in)    :: nspinor
 real(dp),                  intent(in)    :: wk
 type(dataset_type),        intent(inout) :: dtset
 type(gs_hamiltonian_type), intent(in)    :: ham_k
 type(MPI_type),            intent(in)    :: mpi_enreg
 type(bandfft_kpt_type),    intent(in)    :: bandfft
 !arrays
 real(dp),                  intent(in)    :: cg(2,npw*nspinor*nband)
 real(dp),                  intent(inout) :: kin
 real(dp),                  intent(in)    :: occ0(nband)

 !Local variables-------------------------------
 !scalars
 integer  :: displ
 integer  :: iband, ipw
 integer  :: jpw
 integer  :: me_bandfft
 integer  :: shift
 real(dp) :: ar
 !arrays

! ***********************************************************************

 if (dtset%paral_kgb /= 1) then 
   displ = 0
 else
   me_bandfft = xmpi_comm_rank(mpi_enreg%comm_bandspinorfft) 
   displ = bandfft%rdispls(me_bandfft+1)
 end if

 do iband=1, nband
   if (abs(occ0(iband))>tol8) then 
      shift = npw*nspinor*(iband-1)
      !FB: meanvalue_g does the mpi_sum over the bands inside, that's not very efficient since 
      !FB: we could do it only once at the end
      !FB: From Lucas: meanvalue_g seems slow
      !FB: It maybe useful not to use meanvalue_g at all here
      !call meanvalue_g(ar,ham_k%kinpw_k(1+displ:displ+npw*nspinor),0,1,mpi_enreg,npw,nspinor, &
      !               & cg(:,1+shift:shift+npw*nspinor),cg(:,1+shift:shift+npw*nspinor),0)
      ar = zero
      do ipw = 1, npw
         ar = ar + ham_k%kinpw_k(displ+ipw)*(cg(1,shift+ipw)*cg(1,shift+ipw)+cg(2,shift+ipw)*cg(2,shift+ipw))
      end do
      if(nspinor==2)then
         do ipw = 1+npw, 2*npw
            jpw = ipw - npw
            ar = ar + ham_k%kinpw_k(displ+jpw)*(cg(1,shift+ipw)*cg(1,shift+ipw)+cg(2,shift+ipw)*cg(2,shift+ipw))
         end do
      end if
      kin = kin + wk*occ0(iband)*ar
   end if
 end do

end subroutine rttddft_calc_kin
!!***

!!****f* m_rttddft_properties/rttddft_calc_enl
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
!!***

!!****f* m_rttddft_properties/rttddft_calc_occ
!!
!! NAME
!!  rttddft_calc_occ
!!
!! FUNCTION
!!  Computes occupations at time t from cg(t), cg0 and occ0
!!  In NC:
!!    f_{n,k}(t) = \sum_{m} f_{m,k}(0) <\phi_m(0)|\phi_n(t)>
!!  In PAW:
!!    f_{n,k}(t) = \sum_{m} f_{m,k}(0) <\phi_m(0)|S|\phi_n(t)>
!!               = \sum_{m} f_{m,k}(0) [ <\phi_m(0)|\phi_n(t)> +
!!                 \sum_{i} <\phi_m(0)|p_{i}>\sum_jS_{i,j}<p_{j}|\phi_n(t)> ]
!!
!! INPUTS
!!  cg <real(2,npw*nspinor*nband)> = the wavefunction coefficients
!!  cg0 <real(2,npw*nspinor*nband)> = the initial wavefunction coefficients
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  ham_k <gs_hamiltonian_type> = hamiltonian at point k
!!  ikpt <integer> = indice of the considered k-point
!!  ibg <integer> = indice of the considered k-point for cprj
!!  isppol <integer> = indice of the considered spin-polarization
!!  mpi_enreg <MPI_type> = MPI-parallelisation information
!!  nband_k <integer> = number of bands
!!  npw_k <integer> = number of plane waves
!!  nspinor <integer> = dimension of spinors
!!  occ0 <real(nband)> = initial occupations
!!  tdks <type(tdks_type)> = Main RT-TDDFT object
!!
!! OUTPUT
!!  occ <real(nband)> = the occupations at time t
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_calc_occ(cg,cg0,dtset,ham_k,ikpt,ibg,isppol,mpi_enreg,nband_k,npw_k,nspinor,occ,occ0,tdks)

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,                   intent(in)    :: ikpt
 integer,                   intent(in)    :: ibg
 integer,                   intent(in)    :: isppol
 integer,                   intent(in)    :: nband_k
 integer,                   intent(in)    :: npw_k
 integer,                   intent(in)    :: nspinor
 type(dataset_type),        intent(inout) :: dtset
 type(gs_hamiltonian_type), intent(inout) :: ham_k
 type(MPI_type),            intent(inout) :: mpi_enreg
 type(tdks_type), target,   intent(inout) :: tdks
 !arrays
 real(dp),                  intent(inout) :: cg(2,npw_k*nspinor*nband_k)
 real(dp),                  intent(inout) :: cg0(2,npw_k*nspinor*nband_k)
 real(dp),                  intent(out)   :: occ(nband_k)
 real(dp),                  intent(in)    :: occ0(nband_k)

 !Local variables-------------------------------
 !scalars
 integer                     :: iband, jband, ierr, ind
 integer                     :: nband_cprj_k
 integer                     :: natom
 integer                     :: shift
 !parameters for nonlop
 integer, parameter          :: cpopt = 2
 integer, parameter          :: choice = 1
 integer, parameter          :: idir = 1
 integer, parameter          :: nnlout = 1
 integer, parameter          :: paw_opt = 3
 integer, parameter          :: signs = 1
 integer, parameter          :: tim_nonlop = 16
 logical                     :: cprj_paral_band
 !arrays
 real(dp)                    :: csc(2*nband_k)
 real(dp),allocatable        :: gsc(:,:),gvnlxc(:,:)
 real(dp),allocatable        :: enlout(:),enlout_im(:)
 type(pawcprj_type), pointer :: cprj(:,:), cprj_k(:,:)
 type(pawcprj_type), pointer :: cprj0(:,:), cprj0_k(:,:)

! ***********************************************************************

 !Prepare cprj in PAW case
 cprj_paral_band=.false.
 if (ham_k%usepaw == 1) then
    ! Determine if cprj datastructure is distributed over bands
    cprj_paral_band=(tdks%mband_cprj<dtset%mband)
    nband_cprj_k=nband_k; if (cprj_paral_band) nband_cprj_k=mpi_enreg%bandpp
    natom = dtset%natom
    ! Extract the right cprj for this k-point
    if (dtset%mkmem*dtset%nsppol/=1) then
       ABI_MALLOC(cprj_k,(natom,nspinor*nband_cprj_k))
       ABI_MALLOC(cprj0_k,(natom,nspinor*nband_cprj_k))
       call pawcprj_alloc(cprj_k,0,tdks%dimcprj)
       call pawcprj_alloc(cprj0_k,0,tdks%dimcprj)
       call pawcprj_get(tdks%atindx1,cprj_k,tdks%cprj,natom,1,ibg,ikpt,0,isppol,     &
                      & tdks%mband_cprj,dtset%mkmem,natom,nband_cprj_k,nband_cprj_k, &
                      & nspinor,dtset%nsppol,tdks%unpaw,mpicomm=mpi_enreg%comm_kpt,  &
                      & proc_distrb=mpi_enreg%proc_distrb)
       call pawcprj_get(tdks%atindx1,cprj0_k,tdks%cprj0,natom,1,ibg,ikpt,0,isppol,   &
                      & tdks%mband_cprj,dtset%mkmem,natom,nband_cprj_k,nband_cprj_k, &
                      & nspinor,dtset%nsppol,tdks%unpaw,mpicomm=mpi_enreg%comm_kpt,  &
                      & proc_distrb=mpi_enreg%proc_distrb)
    else
       cprj_k => tdks%cprj
       cprj0_k => tdks%cprj0
    end if
    ! If cprj are distributed over bands, gather them (because we need to mix bands)
    if (cprj_paral_band) then
       ABI_MALLOC(cprj,(natom,nspinor*nband_k))
       ABI_MALLOC(cprj0,(natom,nspinor*nband_k))
       call pawcprj_alloc(cprj,0,tdks%dimcprj)
       call pawcprj_alloc(cprj0,0,tdks%dimcprj)
       call pawcprj_mpi_allgather(cprj_k,cprj,natom,nspinor*nband_cprj_k,mpi_enreg%bandpp, &
                                & tdks%dimcprj,0,mpi_enreg%nproc_band,mpi_enreg%comm_band, &
                                & ierr,rank_ordered=.false.)
       call pawcprj_mpi_allgather(cprj0_k,cprj0,natom,nspinor*nband_cprj_k,mpi_enreg%bandpp, &
                                & tdks%dimcprj,0,mpi_enreg%nproc_band,mpi_enreg%comm_band,   &
                                & ierr,rank_ordered=.false.)
    else
       cprj => cprj_k
       cprj0 => cprj0_k
    end if

   !allocate necessary arrays for nonlop
   ABI_MALLOC(gsc,(0,0))
   ABI_MALLOC(gvnlxc,(0,0))
   ABI_MALLOC(enlout,(nband_k))
   ABI_MALLOC(enlout_im,(nband_k))
 end if

 csc = zero
 do iband = 1, nband_k
   !* 1 - Compute csc = <cg0|cg>
   shift = npw_k*nspinor*(iband-1)
   call zgemv('C',npw_k*nspinor,nband_k,cone,cg0,npw_k*nspinor, &
            & cg(:,shift+1:shift+npw_k*nspinor),1,czero,csc,1)
   !If band parallel then reduce csc
   if (mpi_enreg%nproc_band > 1) then
      call xmpi_sum(csc,mpi_enreg%comm_bandfft,ierr)
   end if
   !* 2 - If PAW, add the additional term \sum_i cprj_i \sum_j S_{i,j} cprj_j
   if (ham_k%usepaw == 1) then
      call nonlop(choice,cpopt,cprj(:,iband:iband+(nspinor-1)),enlout,     &
                & ham_k,idir,(/zero/),mpi_enreg,1,nnlout,paw_opt,signs,    &
                & gsc,tim_nonlop,cg(:,shift+1:shift+npw_k*nspinor),gvnlxc, &
                & cprjin_left=cprj0,enlout_im=enlout_im,ndat_left=nband_k)
      do jband = 1, nband_k
         csc(2*jband-1) = csc(2*jband-1) + enlout(jband)
         csc(2*jband)   = csc(2*jband)   + enlout_im(jband)
      end do
   end if
   !* 3 - Calc occupations from csc and occ0
   do jband = 1, nband_k
      ind = 2*jband
      occ(iband) = occ(iband) + occ0(jband)*(csc(ind-1)**2+csc(ind)**2)
   end do
 end do

 if (ham_k%usepaw == 1) then
   ABI_FREE(gsc)
   ABI_FREE(gvnlxc)
   ABI_FREE(enlout)
   ABI_FREE(enlout_im)
   if (cprj_paral_band) then
      call pawcprj_free(cprj)
      ABI_FREE(cprj)
      call pawcprj_free(cprj0)
      ABI_FREE(cprj0)
   end if
   if (dtset%mkmem*dtset%nsppol/=1) then
      call pawcprj_free(cprj_k)
      ABI_FREE(cprj_k)
      call pawcprj_free(cprj0_k)
      ABI_FREE(cprj0_k)
   end if
 end if

end subroutine rttddft_calc_occ
!!***

!!****f* m_rttddft_properties/rttddft_calc_ent
!!
!! NAME
!!  rttddft_calc_ent
!!
!! FUNCTION
!!  Computes electronic entropy from occupation numbers f_{nk}
!!  S = -2 \sum_k \sum_n w(k) [ f_{nk}*ln(f_{nk}) + (1-f_{nk})*ln(1-f_{nk}) ]
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  occ <real(nband*nkpt*nsspol)> = occupation numbers at time t
!!
!! OUTPUT
!!  entropy <real>
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine rttddft_calc_ent(entropy,dtset,occ)

 implicit none

 !Arguments ------------------------------------
 !scalars
 real(dp),           intent(out)   :: entropy
 type(dataset_type), intent(inout) :: dtset
 !arrays
 real(dp),           intent(in)    :: occ(:)

 !Local variables-------------------------------
 !scalars
 integer  :: band_index
 integer  :: iband, isppol, ikpt
 integer  :: nband_k
 real(dp) :: fnk

! ***********************************************************************

 entropy = zero
 band_index=0

 do isppol = 1, dtset%nsppol
   do ikpt = 1, dtset%nkpt
      nband_k = dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
      do iband = 1, nband_k
         fnk = occ(iband+band_index)*0.5_dp
         if( fnk>tol16 .and. (one-fnk)>tol16 ) then
            entropy = entropy - two*dtset%wtk(ikpt)*(fnk*log(fnk)+(one-fnk)*log(one-fnk))
         end if
      end do
      band_index = band_index + nband_k
   end do
 end do

end subroutine rttddft_calc_ent
!!***

end module m_rttddft_properties
!!***
