!!****m* ABINIT/m_dfpttd_berrycurv
!! NAME
!!  m_dfpttd_berrycurv
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2023 ABINIT group (FIXME: add author)
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

module m_dfpttd_berrycurv
    
 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_dtset
 use m_dtfil
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_mpinfo,     only : proc_distrb_cycle
 use m_cgtools,    only : dotprod_g

 implicit none

 public :: dfpttd_berrycurv

 private

! *************************************************************************

contains 
!!***

!!****f* ABINIT/m_dfpttd_berrycurv/dfpttd_berrycurv
!! NAME
!!  dfpttd_berrycurv
!!
!! FUNCTION
!! Compute first-order response function contributions to the time-dispersion
!! 3rd order energy derivatives of the dispersion driver.
!! The main inputs are :
!!   - 1st-order WFs for two perturbations i1pert/i1dir,i2pert/i2dir (cg1,cg2)
!!
!! INPUTS
!!  cg1 = first derivative of cg with respect the perturbation i1pert
!!  cg2 = first derivative of cg with respect the perturbation i2pert
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!          if 2, COMPLEX
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gsqcut=large sphere cut-off
!!  mband = maximum number of bands
!!  mk1mem = maximum number of k points for first-order WF
!!           which can fit in core memory
!!  mpert =maximum number of ipert
!!  mpi_enreg=MPI-parallelisation information
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  natom = number of atoms in unit cell
!!  nfft= number of FFT grid points (for this proc) 
!!  ngfft(1:18)=integer array with FFT box dimensions and other 
!!  nkpt = number of k points
!!  nspden = number of spin-density components
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  npwarr(nkpt) = array holding npw for each k point
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  ucvol=volume of the unit cell
!!
!! OUTPUT
!!  d3etot(2,3,mpert,3,mpert,3,mpert) = third derivatives of the energy tensor
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfpttd_berrycurv(cg1,cg2,d3etot_td,dtset,&
 & mband,mk1mem,mpi_enreg,mpw,nkpt,nspinor,nsppol, & 
 & npwarr,occ)
    
 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: mband,mk1mem,mpw
 integer, intent(in) :: nkpt,nspinor,nsppol
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg

!arrays
 integer,intent(in) ::  npwarr(nkpt)
 real(dp),intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: cg2(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(out) :: d3etot_td(2)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer :: bandtot,iband,icg,ierr,ikpt,isppol,istwf_k,me
 integer :: nband_k,npw_k,offset_cgi,size_wf,spaceworld 
 real(dp) :: doti,dotr,wtk_k                                    
!arrays
 real(dp) :: d3etot_k(2)
 real(dp),allocatable :: cwavef1(:,:),cwavef2(:,:)
 real(dp),allocatable :: occ_k(:)
 
! *************************************************************************

 DBG_ENTER("COLL")
 
!Init parallelism
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt 

!Loop over spins
 d3etot_td=zero
 bandtot = 0
 icg=0
 do isppol = 1, nsppol

!  Loop over k-points
   do ikpt = 1, nkpt

     nband_k = dtset%nband(ikpt+(isppol-1)*nkpt)
     npw_k = npwarr(ikpt)
     istwf_k = dtset%istwfk(ikpt)

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,isppol,mpi_enreg%me)) then
       bandtot = bandtot + nband_k
       cycle ! Skip the rest of the k-point loop
     end if

     ABI_MALLOC(occ_k,(nband_k))
     occ_k(:) = occ(1+bandtot:nband_k+bandtot)
     wtk_k    = dtset%wtk(ikpt)

     d3etot_k=zero
     size_wf= dtset%nspinor*npw_k
     ABI_MALLOC(cwavef1,(2,size_wf))
     ABI_MALLOC(cwavef2,(2,size_wf))

     !Loop over bands
     do iband=1,nband_k
    
       if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= mpi_enreg%me_kpt) cycle
       
       !Select bks wf1
       offset_cgi = (iband-1)*size_wf+icg
       cwavef1(:,:)= cg1(:,1+offset_cgi:size_wf+offset_cgi)
       cwavef2(:,:)= cg2(:,1+offset_cgi:size_wf+offset_cgi)

       !Compute the Berry curvature
       !< u_{i,k}^{\lambda1} | u_{i,k}^{\lambda2} >
       call dotprod_g(dotr,doti,istwf_k,size_wf,2,cwavef1,cwavef2,&
      & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

       d3etot_k(1)=d3etot_k(1)+occ_k(iband)*dotr
       d3etot_k(2)=d3etot_k(2)+occ_k(iband)*doti
 
     end do !iband
   
!    Scale d3etot_k contributions by the kpt weight
     d3etot_k(:)=d3etot_k(:)*wtk_k

!    Add the contribution from each k-point. 
     d3etot_td= d3etot_td + d3etot_k

!    Keep track of total number of bands
     bandtot = bandtot + nband_k

!    Shift arrays memory
     icg=icg+npw_k*dtset%nspinor*nband_k

!    Deallocations
     ABI_FREE(cwavef1)
     ABI_FREE(cwavef2)
     ABI_FREE(occ_k)

   end do !ikpt

 end do !isppol

!=== MPI communications ==================
 if (xmpi_paral==1) then
   call xmpi_sum(d3etot_td,spaceworld,ierr)
 end if
 
 DBG_EXIT("COLL")

end subroutine dfpttd_berrycurv
!!***

end module m_dfpttd_berrycurv
!!***
