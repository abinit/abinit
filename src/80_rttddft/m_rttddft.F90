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
 
 use m_dtfil,           only: datafiles_type
 use m_dtset,           only: dataset_type
 use m_mkrho,           only: mkrho
 use m_paw_mkrho,       only: pawmkrho
 use m_paw_occupancies, only: pawmkrhoij
 use m_rttddft_tdks,    only: tdks_type 

 implicit none

 private
!!***

 public :: calc_density

contains 

!!****f* m_rttddft/calc_density
!!
!! NAME
!! tdks_init
!!
!! FUNCTION
!! Compute electronic density (in 1/bohr^3) from the WF (cg coefficients)
!!
!! INPUTS
!! tdks <class(tdks_type)> = the tdks object to initialize
!! dtfil <type datafiles_type> = infos about file names, file unit numbers
!! dtset <type(dataset_type)> = all input variables for this dataset
!! mpi_enreg <MPI_type> = MPI-parallelisation information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!! m_rttddft_driver
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

end module m_rttddft
!!***
