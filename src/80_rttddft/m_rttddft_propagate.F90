!!****m* ABINIT/m_rttddft_propagate
!! NAME
!!  m_rttddft_propagate
!!
!! FUNCTION
!!  Contains various subroutines to propagate the KS 
!!  orbitals and potentially also the nuclei in RTTDDFT
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

module m_rttddft_propagate

 use defs_basis
 use defs_abitypes,   only: MPI_type
 use defs_datatypes,  only: pseudopotential_type
 
 use m_dtset,         only: dataset_type
 use m_hamiltonian,   only: init_hamiltonian, gs_hamiltonian_type
 use m_rttddft,       only: update_after_nuc_step, update_paw
 use m_rttddft_types, only: tdks_type 
 use m_specialmsg,    only: wrtout
 use m_symtk,         only: symmetrize_xred

 implicit none

 private
!!***

 public :: rttddft_propagate_ele
 public :: rttddft_propagate_nuc

contains 

!!****f* m_rttddft/rttddft_propagate_ele
!!
!! NAME
!!  rttddft_propagate_ele
!!
!! FUNCTION
!!  Main subroutine to propagate time-dependent KS orbitals
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
subroutine rttddft_propagate_ele(tdks, dtset, istep, mpi_enreg, psps)

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
 character(len=500)        :: msg
 integer                   :: bdtot_index
 integer                   :: ibg,icg
 integer                   :: mbdkpsp
 integer                   :: usecprj_local
 type(gs_hamiltonian_type) :: gs_hamk
 !arrays
 real(dp)                  :: grnl(3*dtset%natom)
 real(dp),allocatable      :: eknk(:),eknk_nd(:,:,:,:,:)
 real(dp),allocatable      :: enlxnk(:)
 real(dp),allocatable      :: EigMin(:,:)
 real(dp),allocatable      :: focknk(:),fockfornk(:,:,:)
 real(dp),allocatable      :: grnlnk(:,:)
 
! ***********************************************************************


 write(msg,'(2a,i5,a)') ch10,'--- Iteration',istep,ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 !** Initialize / Update various quantities that needs to be changed 
 !** after the change of xred during the nuclear step
 if (istep == 1 .or. dtset%ionmov /= 0) then
   call update_after_nuc_step(tdks,dtset,istep,mpi_enreg,psps)
 end if

 !** Update PAW quantities
 !** Compute energies and potentials in the augmentation regions (spheres)
 !** and pseudopotential strengths (Dij quantities)
 if (psps%usepaw==1)then
   call update_paw(tdks,dtset,istep,mpi_enreg,psps)
 end if

!tdks%energies%e_eigenvalues = zero
!tdks%energies%e_kinetic     = zero
!tdks%energies%e_nucdip      = zero
!tdks%energies%e_nlpsp_vfock = zero
!if (dtset%usefock/=0) then
!  tdks%energies%e_fock=zero
!  tdks%energies%e_fockdc=zero
!end if
!grnl(:)=zero
!tdks%eigen(:) = zero
!bdtot_index=0
!ibg=0;icg=0
!mbdkpsp=dtset%mband*dtset%nkpt*dtset%nsppol
!ABI_MALLOC(eknk,(mbdkpsp))
!ABI_MALLOC(enlxnk,(mbdkpsp))
!ABI_MALLOC(eknk_nd,(dtset%nsppol,dtset%nkpt,2,dtset%mband,dtset%mband*tdks%paw_dmft%use_dmft))
!ABI_MALLOC(EigMin,(2,dtset%mband))
!ABI_MALLOC(grnlnk,(3*dtset%natom,mbdkpsp*dtset%optforces))
!if (dtset%usefock==0) then
!  ABI_MALLOC(focknk,(mbdkpsp))
!  focknk=zero
!  if (dtset%optforces>0)then
!    ABI_MALLOC(fockfornk,(3,dtset%natom,mbdkpsp))
!    fockfornk=zero
!  end if
!end if
!eknk(:)=zero;enlxnk(:)=zero
!if (dtset%optforces>0) grnlnk(:,:)=zero

!!==== Initialize most of the Hamiltonian ====
!!Allocate all arrays and initialize quantities that do not depend on k and spin.
!usecprj_local=0;if (psps%usepaw==1) usecprj_local=1
!call init_hamiltonian(gs_hamk,psps,tdks%pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
!                    & dtset%typat,tdks%xred,dtset%nfft,dtset%mgfft,dtset%ngfft,tdks%rprimd,dtset%nloalg,&
!                    & paw_ij=tdks%paw_ij,ph1d=tdks%ph1d,usecprj=usecprj_local,fock=tdks%fock,&
!                    & comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
!                    & nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

 !FB: the cprj and rhoij should also be updated once the cg have changed..
 !FB: probably to be done after propagating?

 !FB: Here we should now call the propagator to evolve the cg



 end subroutine rttddft_propagate_ele

!!****f* m_rttddft/rttddft_propagate_nuc
!!
!! NAME
!!  rttddft_propagate_nuc
!!
!! FUNCTION
!!  Main subroutine to propagate nuclei using Ehrenfest dynamics
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
subroutine rttddft_propagate_nuc(tdks, dtset, istep, mpi_enreg, psps)

 implicit none

 !Arguments ------------------------------------
 !scalars
 class(tdks_type),           intent(inout) :: tdks
 integer,                    intent(in)    :: istep
 type(dataset_type),         intent(in)    :: dtset
 type(MPI_type),             intent(inout) :: mpi_enreg
 type(pseudopotential_type), intent(inout) :: psps
 
 !Local variables-------------------------------
 !scalars
 character(len=500)   :: msg
 !arrays
 
! ***********************************************************************

 write(msg,'(2a,i5,a)') ch10,'--- Iteration',istep,ch10
 call wrtout(ab_out,msg)
 if (do_write_log) call wrtout(std_out,msg)

 ! FB: Should we do this? 
 ! Eventually symmetrize atomic coordinates over space group elements:
 call symmetrize_xred(dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,tdks%xred,indsym=tdks%indsym)

 end subroutine rttddft_propagate_nuc

end module m_rttddft_propagate
!!***
