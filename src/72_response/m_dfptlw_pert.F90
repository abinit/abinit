!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfptlw_pert
!! NAME
!!  m_dfptlw_pert
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2022 ABINIT group (MR)
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

module m_dfptlw_pert
    
 use defs_basis
 use defs_abitypes
 use defs_datatypes
 use m_dtset
 use m_dtfil
 use m_errors
 use m_profiling_abi
 use m_hamiltonian
 use m_cgtools
 use m_pawcprj
 use m_pawfgr
 use m_wfk
 use m_xmpi
 use m_getgh1c
 use m_mklocl

 use m_fstrings, only : itoa, sjoin
 use m_io_tools, only : file_exists
 use m_time, only : cwtime
 use m_kg, only : mkkpg


 implicit none

 public :: dfptlw_pert

 private

! *************************************************************************

contains 
!!***

!!****f* ABINIT/m_dfptlw_pert/dfptlw_pert
!! NAME
!!  dfptlw_pert
!!
!! FUNCTION
!! Compute first-order response function contributions to the spatial-dispersion
!! 3rd order energy derivatives of the longwave driver.
!! The main inputs are :
!!   - GS WFs and Hamiltonian (cg,gs_hamkq)
!!   - 1st-order WFs for two perturbations i1pert/i1dir,i2pert/i2dir (cg1,cg2)
!!   - 1st-order Local+SCF potentials for i1pert and i2pert (vtrial1_i1pert,vtrial1_i2pert)
!!   - 1st-order WFs DDK and 2nd-order WF D2_DKDK (d2_dkdk_f)
!!
!! COPYRIGHT
!! Copyright (C) 2018-2021 ABINIT group (MR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cg(2,mpw*nspinor*mband*mkmem_rbz*nsppol) = array for planewave
!!                                          coefficients of wavefunctions
!!  cg1 = first derivative of cg with respect the perturbation i1pert
!!  cg2 = first derivative of cg with respect the perturbation i2pert
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!          if 2, COMPLEX
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  i1dir,i2dir,i3dir=directions of the corresponding perturbations
!!  i1pert,i2pert,i3pert = type of perturbation that has to be computed
!!  kg(3,mpw*mkmem_rbz)=reduced planewave coordinates
!!  mband = maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem_rbz = maximum number of k points which can fit in core memory
!!  mk1mem = maximum number of k points for first-order WF
!!           which can fit in core memory
!!  mpert =maximum number of ipert
!!  mpi_enreg=MPI-parallelisation information
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  natom = number of atoms in unit cell
!!  nattyp(ntypat)= # atoms of each type.
!!  nfft= number of FFT grid points (for this proc) 
!!  ngfft(1:18)=integer array with FFT box dimensions and other 
!!  nkpt = number of k points
!!  nspden = number of spin-density components
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  npwarr(nkpt) = array holding npw for each k point
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  rho1g1(2,nfft)=G-space RF electron density in electrons/bohr**3 (i1pert)
!!  rho2r1(cplex*nfft,nspden)=RF electron density in electrons/bohr**3 (i2pert)
!!  rprimd(3,3) = dimensional primitive translations (bohr)
!!  ucvol=volume of the unit cell
!!  vtrial1_i1pert(cplex*nfft,nspden)=firs-order local potential
!!  vtrial1_i2pert(cplex*nfft,nspden)=firs-order local potential
!!  ddk_f = wf files
!!  d2_dkdk_f = wf files
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density (dummy) 
!!  xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!  d3etot(2,3,mpert,3,mpert,3,mpert) = third derivatives of the energy tensor
!!
!! SIDE EFFECTS
!!  TO DO!
!!
!! PARENTS
!!      m_dfptlw_loop
!!
!! CHILDREN
!!      dotprod_vn
!!
!! SOURCE

subroutine dfptlw_pert(atindx,cg,cg1,cg2,cplex,dtfil,dtset,d3etot,gs_hamkq,i1dir,i2dir,i3dir,&
& i1pert,i2pert,i3pert,kg,mband,mgfft,mkmem_rbz,mk1mem,mpert,mpi_enreg,mpsang,mpw,natom,nattyp,nfft,ngfft,nkpt,&
& nspden,nspinor,nsppol,npwarr,occ,pawfgr,ph1d,psps,rho1g1,rho2r1,rprimd,&
& ucvol,vtrial1_i1pert,vtrial1_i2pert,ddk_f,d2_dkdk_f,xccc3d1,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,mband,mgfft
 integer,intent(in) :: mk1mem,mkmem_rbz,mpert,mpsang,mpw,natom,nfft,nkpt,nspden
 integer,intent(in) :: nspinor,nsppol
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(pawfgr_type),intent(in) :: pawfgr
 type(wfk_t),intent(inout) :: ddk_f,d2_dkdk_f

!arrays
 integer,intent(in) :: atindx(natom),kg(3,mpw*mkmem_rbz),nattyp(psps%ntypat),ngfft(18),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem_rbz*nsppol)
 real(dp),intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: cg2(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rho1g1(2,nfft),rho2r1(cplex*nfft,dtset%nspden)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: xccc3d1(cplex*nfft),xred(3,natom)
 real(dp),intent(in) :: vtrial1_i1pert(cplex*nfft,nspden)
 real(dp),intent(in) :: vtrial1_i2pert(cplex*nfft,nspden)
 real(dp),intent(inout) :: d3etot(2,3,mpert,3,mpert,3,mpert)

!Variables ------------------------------------
!scalars
 integer :: ii,me,spaceworld
 character(len=1000) :: msg
!arrays
 type(rf_hamiltonian_type) :: rf_hamkq_i1pert, rf_hamkq_i2pert
 
! *************************************************************************

 DBG_ENTER("COLL")

!Anounce start of spatial-dispersion calculation
 write(msg, '(a,80a,a,a,a)' ) ch10,('=',ii=1,80),ch10,&
&   ' ==> Compute spatial-dispersion 3rd-order energy derivatives <== ',ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

!Init parallelism
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt 

!Initialize rf_hamiltonians (the k-dependent part is prepared in getgh1c_setup)
 call init_rf_hamiltonian(cplex,gs_hamkq,i1pert,rf_hamkq_i1pert,& 
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
& mpi_spintab=mpi_enreg%my_isppoltab)

 call init_rf_hamiltonian(cplex,gs_hamkq,i2pert,rf_hamkq_i2pert,& 
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
& mpi_spintab=mpi_enreg%my_isppoltab)




















 call rf_hamkq_i1pert%free()
 call rf_hamkq_i2pert%free()

 DBG_EXIT("COLL")

end subroutine dfptlw_pert
!!***

end module m_dfptlw_pert
!!***
