!!****m* ABINIT/m_scf_history
!! NAME
!!  m_scf_history
!!
!! FUNCTION
!!  This module provides the definition of the scf_history_type used to store
!!  various arrays obtained from previous SCF cycles (density, positions, wavefunctions ...),
!!  as needed by the specific SCF algorithm.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2020 ABINIT group (MT)
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

MODULE m_scf_history

 use defs_basis
 use m_abicore
 use m_dtset
 use m_errors

 use defs_abitypes, only : MPI_type
 use m_pawcprj,  only : pawcprj_type, pawcprj_free
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_nullify, pawrhoij_free

 implicit none

 private

! public procedures.
 public :: scf_history_init
 public :: scf_history_free
 public :: scf_history_nullify
!!***

!!****t* m_scf_history/scf_history_type
!! NAME
!! scf_history_type
!!
!! FUNCTION
!! This structured datatype contains various arrays obtained from
!! previous SCF cycles (density, positions...)
!!
!! SOURCE

 type, public :: scf_history_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer scalar

  integer :: history_size
   ! Number of previous SCF cycles stored in history
   ! If history_size<0, scf_history is not used
   ! If history_size=0, scf_history only contains
   !    current values of data (rhor, taur, pawrhoih, xred)
   ! If history_size>0, scf_history contains
   !    current values of data and also
   !    history_size previous values of these data

  integer :: icall
   ! Number of call for the routine extraprho or wf_mixing

  integer :: mcg
   ! Size of cg array

  integer :: mcprj
   ! Size of cprj datsatructure array

  integer :: meigen
   ! Size of eigen array

  integer :: natom
   ! Number of atoms in cell

  integer :: nfft
   ! Size of FFT grid (for density)

  integer :: nspden
   ! Number of independant spin components for density

  integer :: usecg
   ! usecg=0 if the extrapolation/mixing of the density/potential is active but not the one of the wavefunction
   ! usecg=1 if the extrapolation/mixing of the density/potential and wavefunctions is active
   ! usecg=2 if the extrapolation/mixing of the wavefunctions is active but not the one of the density/potential

  integer :: wfmixalg
   ! algorithm used to mix the wavefunctions (in case usecg=2)

  real(dp) :: alpha
   ! alpha mixing coefficient for the prediction of density and wavefunctions
   ! In the case of wavefunction simple mixing, contain wfmix factor

  real(dp) :: beta
   ! beta mixing coefficient for the prediction of density and wavefunctions

! Integer arrays

  integer,allocatable :: hindex(:)
   ! Indexes of SCF cycles in the history
   !
   ! For the density-based schemes (with or without wavefunctions) :
   ! hindex(history_size)
   ! hindex(1) is the newest SCF cycle
   ! hindex(history_size) is the oldest SCF cycle
   !
   ! For wavefunction-based schemes (outer loop of a double loop SCF):
   ! hindex(2*history_size+1)
   ! The odd indices refer to the out wavefunction,
   ! the even indices refer to the in wavefunction (not all such wavefunctions being stored, though).
   ! hindex(1:2) is the newest SCF cycle, hindex(3:4) is the SCF cycle before the newest one ... In case of an
   ! algorithm based on a biorthogonal ensemble of wavefunctions, the reference is stored in hindex(2*history_size+1)
   ! When the index points to a location beyond history_size, the corresponding wavefunction set must be reconstructed
   ! from the existing wavefunctions sets (to be implemented)

! Real (real(dp)) arrays

   real(dp),allocatable :: cg(:,:,:)
    ! cg(2,mcg,history_size)
    ! wavefunction coefficients needed for each SCF cycle of history
    ! Might also contain the wf residuals

   real(dp),allocatable :: deltarhor(:,:,:)
    ! deltarhor(nfft,nspden,history_size)
    ! Difference between electronic density (in real space)
    ! and sum of atomic densities at the end of each SCF cycle of history

   real(dp),allocatable :: eigen(:,:)
    ! eigen(meigen,history_size)
    ! eigenvalues for each SCF cycle of history

   real(dp),allocatable :: atmrho_last(:)
    ! atmrho_last(nfft)
    ! Sum of atomic densities at the end of the LAST SCF cycle

   real(dp),allocatable :: rhor_last(:,:)
    ! rhor_last(nfft,nspden)
    ! Last computed electronic density (in real space)

   real(dp),allocatable :: taur_last(:,:)
    ! taur_last(nfft,nspden*usekden)
    ! Last computed kinetic energy density (in real space)

   real(dp),allocatable :: xreddiff(:,:,:)
    ! xreddiff(3,natom,history_size)
    ! Difference of reduced coordinates of atoms between a
    ! SCF cycle and the previous

   real(dp),allocatable :: xred_last(:,:)
    ! xred_last(3,natom)
    ! Last computed atomic positions (reduced coordinates)

   real(dp),allocatable :: dotprod_sumdiag_cgcprj_ij(:,:,:)
    ! dotprod_sumdiag_cgcprj_mn(2,history_size,history_size)
    ! Container for the scalar products between aligned sets of wavefunctions or their residuals
    ! S_ij=Sum_nk <wf_nk(set i)|wf_nk(set j)> possibly with some weighting factor that might depend on nk.

! Structured datatypes arrays

  type(pawrhoij_type), allocatable :: pawrhoij(:,:)
    ! pawrhoij(natom,history_size)
    ! PAW only: occupancies matrix at the end of each SCF cycle of history

  type(pawrhoij_type), allocatable :: pawrhoij_last(:)
    ! pawrhoij_last(natom)
    ! PAW only: last computed occupancies matrix

  type(pawcprj_type),allocatable :: cprj(:,:,:)
    !cprj(natom,nspinor*mband*mkmem*nsppol,history_size)

 end type scf_history_type
!!***

CONTAINS !===========================================================
!!***

!!****f* m_scf_history/scf_history_init
!! NAME
!!  scf_history_init
!!
!! FUNCTION
!!  Init all scalars and pointers in a scf_history datastructure
!!  according to scf_history%history_size value which has to be
!!  defined before calling this routine
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  mpi_enreg=MPI-parallelisation information
!!  usecg= if ==0 => no handling of wfs (and eigenvalues),
!!         if==1 => handling of density/potential AND wfs and eigen,
!!         if==2 => ONLY handling of wfs and eigen
!!
!! SIDE EFFECTS
!!  scf_history=<type(scf_history_type)>=scf_history datastructure
!!    hindex is always allocated
!!    The density/potential arrays that are possibly allocated are : atmrho_last, deltarhor,
!!      pawrhoij, pawrhoij_last, rhor_last, taur_last, xreddiff, xred_last.
!!    The wfs arrays that are possibly allocated are : cg, cprj and eigen
!!
!! PARENTS
!!      m_gstate,m_scfcv_core
!!
!! CHILDREN
!!
!! SOURCE

subroutine scf_history_init(dtset,mpi_enreg,usecg,scf_history)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: usecg
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 type(scf_history_type),intent(inout) :: scf_history

!Local variables-------------------------------
!scalars
 integer :: jj,mband_cprj,my_natom,my_nspinor,nfft
!arrays

!************************************************************************

 !@scf_history_type

 if (scf_history%history_size<0) then
   call scf_history_nullify(scf_history)
 else

   scf_history%usecg=usecg
   scf_history%wfmixalg=dtset%fockoptmix/100

   nfft=dtset%nfft
   if (dtset%usepaw==1.and.(dtset%pawecutdg>=1.0000001_dp*dtset%ecut)) nfft=dtset%nfftdg
   my_natom=mpi_enreg%my_natom

   if (scf_history%history_size>=0 .and. usecg<2) then
     ABI_MALLOC(scf_history%rhor_last,(nfft,dtset%nspden))
     ABI_MALLOC(scf_history%taur_last,(nfft,dtset%nspden*dtset%usekden))
     ABI_MALLOC(scf_history%xred_last,(3,dtset%natom))
     ABI_MALLOC(scf_history%pawrhoij_last,(my_natom*dtset%usepaw))
     if (dtset%usepaw==1) then
       call pawrhoij_nullify(scf_history%pawrhoij_last)
     end if
   end if

   if (scf_history%history_size>0) then

     scf_history%natom=dtset%natom
     scf_history%nfft=nfft
     scf_history%nspden=dtset%nspden
     scf_history%beta=zero
     scf_history%icall=0

     scf_history%mcg=0
     scf_history%mcprj=0
     scf_history%meigen=0
     if (usecg>0) then
       my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
       scf_history%meigen=dtset%nbandhf*dtset%mkmem*dtset%nsppol
       scf_history%mcg=dtset%mpw*my_nspinor*scf_history%meigen ! This is for scf_history_wf
       if(usecg==1)then
         scf_history%meigen=dtset%mband*dtset%mkmem*dtset%nsppol
         scf_history%mcg=dtset%mpw*my_nspinor*scf_history%meigen ! This is for scf_history (when extrapwf==1)
       endif
       if (dtset%usepaw==1) then
         mband_cprj=dtset%nbandhf
         if(usecg==1)mband_cprj=dtset%mband
         if (dtset%paral_kgb/=0) mband_cprj=mband_cprj/mpi_enreg%nproc_band
         scf_history%mcprj=my_nspinor*mband_cprj*dtset%mkmem*dtset%nsppol
       end if
     end if

     if (usecg<2) then
       ABI_MALLOC(scf_history%hindex,(scf_history%history_size))
       scf_history%alpha=zero
     else
       ABI_MALLOC(scf_history%hindex,(2*scf_history%history_size+1))
       scf_history%alpha=dtset%wfmix
     endif
     scf_history%hindex(:)=0

     if (usecg<2) then
       ABI_MALLOC(scf_history%deltarhor,(nfft,dtset%nspden,scf_history%history_size))
       ABI_MALLOC(scf_history%xreddiff,(3,dtset%natom,scf_history%history_size))
       ABI_MALLOC(scf_history%atmrho_last,(nfft))
       if (dtset%usepaw==1) then
         ABI_MALLOC(scf_history%pawrhoij,(my_natom,scf_history%history_size))
         do jj=1,scf_history%history_size
           call pawrhoij_nullify(scf_history%pawrhoij(:,jj))
         end do
       endif
     end if

     if (scf_history%usecg>0) then
       ABI_MALLOC(scf_history%cg,(2,scf_history%mcg,scf_history%history_size))
       ABI_MALLOC(scf_history%eigen,(scf_history%meigen,scf_history%history_size))
!      Note that the allocation is made even when usepaw==0. Still, scf_history%mcprj=0 ...
       ABI_MALLOC(scf_history%cprj,(dtset%natom,scf_history%mcprj,scf_history%history_size))
     end if

     if (scf_history%usecg==2)then
!      This relatively small matrix is always allocated when usecg==1, even if not used
       ABI_MALLOC(scf_history%dotprod_sumdiag_cgcprj_ij,(2,scf_history%history_size,scf_history%history_size))
     endif

   end if
 end if

end subroutine scf_history_init
!!***

!----------------------------------------------------------------------

!!****f* m_scf_history/scf_history_free
!! NAME
!!  scf_history_free
!!
!! FUNCTION
!!  Clean and destroy a scf_history datastructure
!!
!! SIDE EFFECTS
!!  scf_history(:)=<type(scf_history_type)>=scf_history datastructure
!!
!! PARENTS
!!      m_gstateimg,m_scfcv_core
!!
!! CHILDREN
!!
!! SOURCE

subroutine scf_history_free(scf_history)

!Arguments ------------------------------------
!arrays
 type(scf_history_type),intent(inout) :: scf_history

!Local variables-------------------------------
!scalars
 integer :: jj

!************************************************************************

 !@scf_history_type

 if (allocated(scf_history%pawrhoij_last)) then
   call pawrhoij_free(scf_history%pawrhoij_last)
   ABI_FREE(scf_history%pawrhoij_last)
 end if
 if (allocated(scf_history%pawrhoij)) then
   do jj=1,size(scf_history%pawrhoij,2)
     call pawrhoij_free(scf_history%pawrhoij(:,jj))
   end do
   ABI_FREE(scf_history%pawrhoij)
 end if
 if (allocated(scf_history%cprj)) then
   do jj=1,size(scf_history%cprj,3)
     call pawcprj_free(scf_history%cprj(:,:,jj))
   end do
   ABI_FREE(scf_history%cprj)
 end if

 if (allocated(scf_history%hindex))       then
   ABI_FREE(scf_history%hindex)
 end if
 if (allocated(scf_history%deltarhor))    then
   ABI_FREE(scf_history%deltarhor)
 end if
 if (allocated(scf_history%xreddiff))     then
   ABI_FREE(scf_history%xreddiff)
 end if
 if (allocated(scf_history%atmrho_last))  then
   ABI_FREE(scf_history%atmrho_last)
 end if
 if (allocated(scf_history%xred_last))    then
   ABI_FREE(scf_history%xred_last)
 end if
 if (allocated(scf_history%rhor_last))    then
   ABI_FREE(scf_history%rhor_last)
 end if
 if (allocated(scf_history%taur_last))    then
   ABI_FREE(scf_history%taur_last)
 end if
 if (allocated(scf_history%cg))           then
   ABI_FREE(scf_history%cg)
 end if
 if (allocated(scf_history%eigen))           then
   ABI_FREE(scf_history%eigen)
 end if
 if (allocated(scf_history%dotprod_sumdiag_cgcprj_ij))           then
   ABI_FREE(scf_history%dotprod_sumdiag_cgcprj_ij)
 end if

 scf_history%history_size=-1
 scf_history%usecg=0
 scf_history%icall=0
 scf_history%mcprj=0
 scf_history%mcg=0
 scf_history%meigen=0

end subroutine scf_history_free
!!***

!----------------------------------------------------------------------

!!****f* m_scf_history/scf_history_nullify
!! NAME
!!  scf_history_nullify
!!
!! FUNCTION
!!  Nullify (set to null) an scf_history datastructure
!!
!! SIDE EFFECTS
!!  scf_history(:)=<type(scf_history_type)>=scf_history datastructure
!!
!! PARENTS
!!      m_gstateimg,m_scf_history
!!
!! CHILDREN
!!
!! SOURCE

subroutine scf_history_nullify(scf_history)

!Arguments ------------------------------------
!arrays
 type(scf_history_type),intent(inout) :: scf_history
!Local variables-------------------------------
!scalars

!************************************************************************

 !@scf_history_type
 scf_history%history_size=-1
 scf_history%icall=0
 scf_history%mcprj=0
 scf_history%mcg=0
 scf_history%meigen=0

end subroutine scf_history_nullify
!!***

END MODULE m_scf_history
!!***
