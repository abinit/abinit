!!****m* ABINIT/m_dfptlw_loop
!! NAME
!!  m_dfptlw_loop
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2022 ABINIT group (MR)
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

module m_dfptlw_loop

 implicit none

 private
!!***

 public :: dfptlw_loop
!!***
! *************************************************************************

contains 
!!***

!!****f* ABINIT/m_dfptlw_loop/dfptlw_loop
!! NAME
!! dfptlw_loop 
!!
!! FUNCTION
!!  Loop over two perturbations j1, j2 and a q gradient
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol) = array for planewave coefficients of wavefunctions
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates
!!  kxc(nfftf,nkxc)=exchange-correlation kernel
!!  mband = maximum number of bands
!!  mgfft = maximum single fft dimension
!!  mkmem = Number of k points treated by this node.
!!  mk1mem = Number of k points for first-order WF treated by this node.
!!  mpert =maximum number of ipert
!!  mpi_enreg=MPI-parallelisation information
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  natom = number of atoms in unit cell
!!  nattyp(ntypat)= # atoms of each type.
!!  nfftf=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!!  ngfftf(1:18)=integer array with FFT box dimensions and other for the "fine" grid (see NOTES in respfn.F90)
!!  nhat=compensation charge density on fine rectangular grid
!!  nkpt  = number of k points
!!  nkxc=second dimension of the array kxc, see rhotoxc.f for a description
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  npwarr(nkpt) = array holding npw for each k point
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data for the GS
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  rfpert(3,mpert,3,mpert,3,mpert) = array defining the type of perturbations
!!       that have to be computed
!!       1   ->   element has to be computed explicitely
!!      -1   ->   use symmetry operations to obtain the corresponding element
!!  rhog(2,nfftf)=array for Fourier transform of GS electron density
!!  rhor(nfftf,nspden)=array for GS electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations (bohr)
!!  ucvol = unit cell volume (bohr^3)
!!  vxc(nfftf,nspden)=Exchange-Correlation GS potential (Hartree)
!!  xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert) = flags for each element of the 3DTE
!!                             (=1 if computed)
!!  d3etot(2,3,mpert,3,mpert,3,mpert) = third derivatives of the energy tensor
!!
!! SIDE EFFECTS
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!
!! PARENTS
!!      m_longwave
!!
!! CHILDREN
!!
!! SOURCE

    
subroutine dfptlw_loop(atindx,blkflg,cg,dtfil,dtset,d3etot,eigen0,gmet,gprimd,&
& hdr,kg,kxc,mband,mgfft,mgfftf,mkmem,mk1mem,&
& mpert,mpi_enreg,mpw,natom,nattyp,ngfftf,nfftf,nhat,nkpt,nkxc,nspinor,nsppol,&
& npwarr,occ,&
& pawfgr,pawrad,pawrhoij,pawtab,&
& psps,rfpert,rhog,rhor,rprimd,ucvol,vxc,xred)

 use defs_basis
 use defs_wvltypes
 use m_errors
 use m_abicore
 use m_hdr
 use m_nctk
 use m_wffile
 use m_wfk
 use m_dtset
 use m_dtfil

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes, only : MPI_type
 use m_time,        only : timab
 use m_io_tools,    only : file_exists
 use m_kg,          only : getph
 use m_inwffil,     only : inwffil
 use m_fft,         only : fourdp
 use m_ioarr,       only : read_rhor
 use m_hamiltonian, only : gs_hamiltonian_type, init_hamiltonian
 use m_pawdij,      only : pawdij, pawdijfr, symdij
 use m_pawfgr,      only : pawfgr_type
 use m_pawfgrtab,   only : pawfgrtab_type
 use m_paw_an,      only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify, paw_an_reset_flags
 use m_paw_ij,      only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify, paw_ij_reset_flags, paw_ij_print
 use m_pawang,      only : pawang_type
 use m_pawrad,      only : pawrad_type
 use m_pawrhoij,    only : pawrhoij_type, pawrhoij_alloc, pawrhoij_free, pawrhoij_nullify, &
&                          pawrhoij_io, pawrhoij_inquire_dim
 use m_paw_nhat,    only : pawmknhat,pawnhatfr
 use m_paw_denpot,  only : pawdenpot
 use m_pawtab,      only : pawtab_type
 use m_rf2,         only : rf2_getidir
 use m_initylmg,    only : initylmg
 use m_atm2fft,     only : dfpt_atm2fft
 use m_dfpt_mkvxc,  only : dfpt_mkvxc
 use m_dfpt_rhotov, only : dfpt_rhotov
 use m_mkcore,      only : dfpt_mkcore
 use m_mklocl,      only : dfpt_vlocal
 use m_dfptnl_pert, only : dfptnl_pert

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfft,mgfftf,mk1mem,mkmem,mpert,mpw,natom,nfftf
 integer,intent(in) :: nkpt,nkxc,nspinor,nsppol
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps

!arrays
 integer,intent(in) :: atindx(natom),kg(3,mk1mem*mpw)
 integer,intent(in) :: nattyp(psps%ntypat),ngfftf(18),npwarr(nkpt)
 integer,intent(in) :: rfpert(3,mpert,3,mpert,3,mpert)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert,3,mpert) 
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),gmet(3,3)
 real(dp),intent(in) :: eigen0(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: gprimd(3,3),kxc(nfftf,nkxc)
 real(dp),intent(in) :: nhat(nfftf,dtset%nspden)
 real(dp),intent(in) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(in) :: vxc(nfftf,dtset%nspden)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
 real(dp),intent(inout) :: d3etot(2,3,mpert,3,mpert,3,mpert) 
 type(pawrhoij_type),intent(in) :: pawrhoij(natom*psps%usepaw)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
! integer ::                                      
! real(dp) ::                                     
!character(len=500) :: msg                   
 
! *************************************************************************

 DBG_ENTER("COLL")
 
! if (option/=1 .and. option/=2 ) then
!  write(msg,'(3a,i0)')&
!&  'The argument option should be 1 or 2,',ch10,&
!&  'however, option=',option
!  MSG_BUG(msg)
! end if
!
! if (sizein<1) then
!  write(msg,'(3a,i0)')&
!&  '  The argument sizein should be a positive number,',ch10,&
!&  '  however, sizein=',sizein
!  MSG_ERROR(msg)
! end if

 DBG_EXIT("COLL")

end subroutine dfptlw_loop
!!***

end module m_dfptlw_loop
!!***
