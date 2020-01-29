!!****m* ABINIT/m_dfptnl_loop
!! NAME
!!  m_dfptnl_loop
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2018-2020 ABINIT group (LB)
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

module m_dfptnl_loop

 implicit none

 private
!!***

 public :: dfptnl_loop
!!***

contains
!!***

!!****f* ABINIT/dfptnl_loop
!! NAME
!! dfptnl_loop
!!
!! FUNCTION
!! Loop over the perturbations j1, j2 and j3
!!
!! COPYRIGHT
!! Copyright (C) 2018-2020 ABINIT group (LB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol) = array for planewave coefficients of wavefunctions
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1)
!!  gsqcut=Fourier cutoff on G^2 for "large sphere" of radius double
!!   that of the basis sphere--appropriate for charge density rho(G),
!!   Hartree potential, and pseudopotentials
!!  kg(3,mpw*mkmem)=reduced planewave coordinates
!!  kxc(nfftf,nkxc)=exchange-correlation kernel
!!  k3xc(nfftf,nk3xc)=third-order exchange-correlation kernel
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
!!  nk3xc=second dimension of the array k3xc
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  npwarr(nkpt) = array holding npw for each k point
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  paw_an0(natom) <type(paw_an_type)>=paw arrays for 0th-order quantities given on angular mesh
!!  paw_ij0(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawang1 <type(pawang_type)>=pawang datastructure containing only the symmetries preserving the perturbation
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid for the GS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data for the GS
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  ph1df(2,3*(2*mgfftf+1)*natom)=one-dimensional structure factor information (fine grid)
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  rfpert(3,mpert,3,mpert,3,mpert) = array defining the type of perturbations
!!       that have to be computed
!!       1   ->   element has to be computed explicitely
!!      -1   ->   use symmetry operations to obtain the corresponding element
!!  rhog(2,nfftf)=array for Fourier transform of GS electron density
!!  rhor(nfftf,nspden)=array for GS electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations (bohr)
!!  ucvol = unit cell volume (bohr^3)
!!  usecprj= 1 if cprj, cprjq, cprj1 arrays are stored in memory
!!  vtrial(nfftf,nspden)=GS Vtrial(r).
!!  vxc(nfftf,nspden)=Exchange-Correlation GS potential (Hartree)
!!  xred(3,natom) = reduced atomic coordinates
!!  nsym1=number of symmetry elements in space group consistent with perturbation
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  symaf1(nsym1)=anti(ferromagnetic) part of symmetry operations
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert) = flags for each element of the 3DTE
!!                             (=1 if computed)
!!  d3etot(2,3,mpert,3,mpert,3,mpert) = third derivatives of the energy tensor
!!                                    = \sum_{i=1}^9 d3etot_i
!!  d3etot_1(2,3,mpert,3,mpert,3,mpert) = 1st term of d3etot
!!  d3etot_2(2,3,mpert,3,mpert,3,mpert) = 2nd term of d3etot
!!  d3etot_3(2,3,mpert,3,mpert,3,mpert) = 3rd term of d3etot
!!  d3etot_4(2,3,mpert,3,mpert,3,mpert) = 4th term of d3etot
!!  d3etot_5(2,3,mpert,3,mpert,3,mpert) = 5th term of d3etot
!!  d3etot_6(2,3,mpert,3,mpert,3,mpert) = 6th term of d3etot
!!  d3etot_7(2,3,mpert,3,mpert,3,mpert) = 7th term of d3etot
!!  d3etot_8(2,3,mpert,3,mpert,3,mpert) = 8th term of d3etot
!!  d3etot_9(2,3,mpert,3,mpert,3,mpert) = 9th term of d3etot
!!
!! SIDE EFFECTS
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!
!! PARENTS
!!      nonlinear
!!
!! CHILDREN
!!      appdig,dfpt_mkcore,dfpt_mkvxc,dfpt_vlocal,dfptnl_mv,dfptnl_resp
!!      dotprod_vn,fourdp,getph,hartre,hdr_free,initylmg,inwffil,read_rhor
!!      status,timab,wffclose,wrtout
!!
!! SOURCE

subroutine dfptnl_loop(atindx,blkflg,cg,dtfil,dtset,d3etot,eigen0,gmet,gprimd,gsqcut, &
& hdr,kg,kxc,k3xc,mband,mgfft,mgfftf,mkmem,mk1mem,&
& mpert,mpi_enreg,mpw,natom,nattyp,ngfftf,nfftf,nhat,nkpt,nkxc,nk3xc,nspinor,nsppol,&
& npwarr,occ,paw_an0,paw_ij0,&
& pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,&
& ph1d,ph1df,psps,rfpert,rhog,rhor,rprimd,ucvol,usecprj,vtrial,vxc,xred,&
& nsym1,indsy1,symaf1,symrc1,&
& d3etot_1,d3etot_2,d3etot_3,d3etot_4,d3etot_5,d3etot_6,d3etot_7,d3etot_8,d3etot_9)

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

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfft,mgfftf,mk1mem,mkmem,mpert,mpw,natom,nfftf
 integer,intent(in) :: nk3xc,nkpt,nkxc,nspinor,nsppol,nsym1,usecprj
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(inout) :: pawang,pawang1
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps

!arrays
 integer,intent(in) :: atindx(natom),kg(3,mk1mem*mpw)
 integer,intent(in) :: nattyp(psps%ntypat),ngfftf(18),npwarr(nkpt)
 integer,intent(in) :: rfpert(3,mpert,3,mpert,3,mpert)
 integer,intent(in) :: indsy1(4,nsym1,dtset%natom),symaf1(nsym1),symrc1(3,3,nsym1)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert,3,mpert) !vz_i
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),gmet(3,3)
 real(dp),intent(in) :: eigen0(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: gprimd(3,3),k3xc(nfftf,nk3xc),kxc(nfftf,nkxc)
 real(dp),intent(in) :: nhat(nfftf,dtset%nspden)
 real(dp),intent(in) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),ph1df(2,3*(2*mgfftf+1)*natom)
 real(dp),intent(in) :: vtrial(nfftf,dtset%nspden),xred(3,natom)
 real(dp),intent(in) :: vxc(nfftf,dtset%nspden)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
 real(dp),intent(inout) :: d3etot(2,3,mpert,3,mpert,3,mpert) !vz_i
 real(dp),intent(inout) :: d3etot_1(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(inout) :: d3etot_2(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(inout) :: d3etot_3(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(inout) :: d3etot_4(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(inout) :: d3etot_5(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(inout) :: d3etot_6(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(inout) :: d3etot_7(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(inout) :: d3etot_8(2,3,mpert,3,mpert,3,mpert)
 real(dp),intent(inout) :: d3etot_9(2,3,mpert,3,mpert,3,mpert)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom*psps%usepaw)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_an_type),intent(in) :: paw_an0(natom*psps%usepaw)
 type(paw_ij_type),intent(in) :: paw_ij0(natom*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=51
 integer :: ask_accurate,comm_cell,counter,cplex,cplex_rhoij,formeig,flag1,flag3
 integer :: has_dijfr,has_diju
 integer :: i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,iatom,idir_dkde,ierr,ii,ireadwf
 integer :: mcg,mpsang,n1,n2,n3,n3xccc,ndir,nfftotf,nhat1grdim,npert_phon,nspden,nspden_rhoij,nwffile
 integer :: option,optene,optfr,optorth,pert1case,pert2case,pert3case
 integer :: qphase_rhoij,rdwrpaw,second_idir,timrev,usexcnhat
 logical :: non_magnetic_xc
 real(dp) :: dummy_real,dummy_real2, dummy_real3, ecut_eff
 character(len=500) :: message
 character(len=fnlen) :: fiden1i,fiwf1i,fiwf2i,fiwf3i,fiwfddk,fnamewff(5)
 type(gs_hamiltonian_type) :: gs_hamkq
 type(wffile_type) :: wff1,wff2,wff3,wfft1,wfft2,wfft3
 type(wfk_t) :: ddk_f(5)
 type(wvl_data) :: wvl
 type(hdr_type) :: hdr_den
!arrays
 integer :: file_index(5)
 real(dp) :: qphon(3),tsec(2)
 real(dp),allocatable :: cg1(:,:),cg2(:,:),cg3(:,:),eigen1(:),eigen2(:),eigen3(:)
 real(dp),allocatable :: nhat1_i1pert(:,:),nhat1_i2pert(:,:),nhat1_i3pert(:,:)
 real(dp),allocatable :: nhat1gr(:,:,:),vresid_dum(:,:)
 real(dp),allocatable :: rho1r1(:,:)
 real(dp),allocatable :: rho2g1(:,:),rho2r1(:,:),rho3r1(:,:),vhartr1_i2pert(:)
 real(dp),allocatable :: vpsp1(:),vxc1_i2pert(:,:),work(:)
 real(dp),allocatable,target :: vtrial1_i2pert(:,:)
 real(dp),pointer :: vtrial1_tmp(:,:)
 real(dp),allocatable :: xccc3d1(:),xccc3d2(:),xccc3d3(:)
 type(pawrhoij_type),allocatable :: pawrhoij1_i1pert(:),pawrhoij1_i2pert(:),pawrhoij1_i3pert(:)
 type(paw_an_type),allocatable :: paw_an1_i2pert(:)
 type(paw_ij_type),allocatable :: paw_ij1_i2pert(:)

! ***********************************************************************

 DBG_ENTER("COLL")

 call timab(503,1,tsec)

 comm_cell = mpi_enreg%comm_cell

 timrev = 1 ! as q=0
 cplex = 2 - timrev
 nspden = dtset%nspden
 ecut_eff = (dtset%ecut)*(dtset%dilatmx)**2
 mpsang = psps%mpsang
 optorth=1;if (psps%usepaw==1) optorth=0

 qphon(:)=zero

 ABI_ALLOCATE(cg1,(2,dtset%mpw*dtset%nspinor*mband*dtset%mk1mem*dtset%nsppol))
 ABI_ALLOCATE(cg2,(2,dtset%mpw*dtset%nspinor*mband*dtset%mk1mem*dtset%nsppol))
 ABI_ALLOCATE(cg3,(2,dtset%mpw*dtset%nspinor*mband*dtset%mk1mem*dtset%nsppol))
 ABI_ALLOCATE(eigen1,(2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_ALLOCATE(eigen2,(2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_ALLOCATE(eigen3,(2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_ALLOCATE(rho1r1,(cplex*nfftf,dtset%nspden))
 ABI_ALLOCATE(rho2r1,(cplex*nfftf,dtset%nspden))
 ABI_ALLOCATE(rho2g1,(2,nfftf))
 ABI_ALLOCATE(rho3r1,(cplex*nfftf,dtset%nspden))

 ask_accurate=1 ; formeig = 1 ; ireadwf = 1
 n1=ngfftf(1) ; n2=ngfftf(2) ; n3=ngfftf(3)
 nfftotf=n1*n2*n3

!==== Initialize most of the Hamiltonian (and derivative) ====
!1) Allocate all arrays and initialize quantities that do not depend on k and spin.
!2) Perform the setup needed for the non-local factors:
!* Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
!* PAW: Initialize the overlap coefficients and allocate the Dij coefficients.
 call init_hamiltonian(gs_hamkq,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,natom,&
& dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,&
& usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda,paw_ij=paw_ij0)

 ABI_ALLOCATE(vpsp1,(cplex*nfftf))
 ABI_ALLOCATE(xccc3d1,(cplex*nfftf))
 ABI_ALLOCATE(xccc3d2,(cplex*nfftf))
 ABI_ALLOCATE(xccc3d3,(cplex*nfftf))
 ABI_ALLOCATE(vhartr1_i2pert,(cplex*nfftf))
 ABI_ALLOCATE(vxc1_i2pert,(cplex*nfftf,dtset%nspden))
 ABI_ALLOCATE(vtrial1_i2pert,(cplex*nfftf,dtset%nspden))

 ABI_ALLOCATE(vresid_dum,(0,0))
! PAW stuff
 usexcnhat = 0
 nhat1grdim=0
 ABI_ALLOCATE(nhat1gr,(0,0,0))
 nhat1gr(:,:,:) = zero
 rdwrpaw=psps%usepaw
!Allocate 1st-order PAW occupancies (rhoij1)
 if (psps%usepaw==1) then
   call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,qphase_rhoij=qphase_rhoij,nspden_rhoij=nspden_rhoij,&
&                        nspden=dtset%nspden,spnorb=dtset%pawspnorb,cplex=cplex,cpxocc=dtset%pawcpxocc)
   ABI_DATATYPE_ALLOCATE(pawrhoij1_i1pert,(natom))
   ABI_DATATYPE_ALLOCATE(pawrhoij1_i2pert,(natom))
   ABI_DATATYPE_ALLOCATE(pawrhoij1_i3pert,(natom))
   call pawrhoij_nullify(pawrhoij1_i1pert)
   call pawrhoij_nullify(pawrhoij1_i2pert)
   call pawrhoij_nullify(pawrhoij1_i3pert)
   call pawrhoij_alloc(pawrhoij1_i1pert,cplex_rhoij,nspden_rhoij,dtset%nspinor,dtset%nsppol,&
&   dtset%typat,qphase=qphase_rhoij,pawtab=pawtab,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call pawrhoij_alloc(pawrhoij1_i2pert,cplex_rhoij,nspden_rhoij,dtset%nspinor,dtset%nsppol,&
&   dtset%typat,qphase=qphase_rhoij,pawtab=pawtab,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call pawrhoij_alloc(pawrhoij1_i3pert,cplex_rhoij,nspden_rhoij,dtset%nspinor,dtset%nsppol,&
&   dtset%typat,qphase=qphase_rhoij,pawtab=pawtab,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 else
   ABI_DATATYPE_ALLOCATE(pawrhoij1_i1pert,(0))
   ABI_DATATYPE_ALLOCATE(pawrhoij1_i2pert,(0))
   ABI_DATATYPE_ALLOCATE(pawrhoij1_i3pert,(0))
 end if

 mcg=mpw*nspinor*mband*mkmem*nsppol

!Allocations/initializations for PAW only
 if(psps%usepaw==1) then
   usexcnhat=maxval(pawtab(:)%usexcnhat)
!  1st-order compensation density
   ABI_ALLOCATE(nhat1_i1pert,(cplex*nfftf,dtset%nspden))
   nhat1_i1pert=zero
   ABI_ALLOCATE(nhat1_i2pert,(cplex*nfftf,dtset%nspden))
   nhat1_i2pert=zero
   ABI_ALLOCATE(nhat1_i3pert,(cplex*nfftf,dtset%nspden))
   nhat1_i3pert=zero

!  1st-order arrays/variables related to the PAW spheres
   ABI_DATATYPE_ALLOCATE(paw_an1_i2pert,(natom))
   ABI_DATATYPE_ALLOCATE(paw_ij1_i2pert,(natom))
   call paw_an_nullify(paw_an1_i2pert)
   call paw_ij_nullify(paw_ij1_i2pert)

   has_dijfr=1
   has_diju=merge(0,1,dtset%usepawu==0)

   call paw_an_init(paw_an1_i2pert,dtset%natom,dtset%ntypat,0,0,dtset%nspden,cplex,dtset%pawxcdev,&
&   dtset%typat,pawang,pawtab,has_vxc=1,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

   call paw_ij_init(paw_ij1_i2pert,cplex,dtset%nspinor,dtset%nsppol,dtset%nspden,0,dtset%natom,&
&   dtset%ntypat,dtset%typat,pawtab,&
&   has_dij=1,has_dijhartree=1,has_dijfr=has_dijfr,has_dijU=has_diju,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
 else
   ABI_ALLOCATE(nhat1_i1pert,(0,0))
   ABI_ALLOCATE(nhat1_i2pert,(0,0))
   ABI_ALLOCATE(nhat1_i3pert,(0,0))
   ABI_DATATYPE_ALLOCATE(paw_an1_i2pert,(0))
   ABI_DATATYPE_ALLOCATE(paw_ij1_i2pert,(0))
 end if ! PAW

 n3xccc=0;if(psps%n1xccc/=0)n3xccc=nfftf
 non_magnetic_xc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)

!Loop over the perturbations j1, j2, j3

 pert1case = 0 ; pert2case = 0 ; pert3case = 0

 do i1pert = 1, mpert
   do i1dir = 1, 3

     if ((maxval(rfpert(i1dir,i1pert,:,:,:,:))==1)) then

       pert1case = i1dir + (i1pert-1)*3
       counter = pert1case
       call appdig(pert1case,dtfil%fnamewff1,fiwf1i)

       call inwffil(ask_accurate,cg1,dtset,dtset%ecut,ecut_eff,eigen1,dtset%exchn2n3d,&
&       formeig,hdr,ireadwf,dtset%istwfk,kg,dtset%kptns,dtset%localrdwf,&
&       dtset%mband,mcg,dtset%mk1mem,mpi_enreg,mpw,&
&       dtset%nband,dtset%ngfft,dtset%nkpt,npwarr,&
&       dtset%nsppol,dtset%nsym,&
&       occ,optorth,dtset%symafm,dtset%symrel,dtset%tnons,&
&       dtfil%unkg1,wff1,wfft1,dtfil%unwff1,fiwf1i,wvl)

       if (ireadwf==1) then
         call WffClose (wff1,ierr)
       end if

       flag1 = 0
       rho1r1(:,:) = zero
       if (dtset%get1den /= 0 .or. dtset%ird1den /= 0) then
         call appdig(pert1case,dtfil%fildens1in,fiden1i)

         call read_rhor(fiden1i, cplex, dtset%nspden, nfftf, ngfftf, rdwrpaw, mpi_enreg, rho1r1, &
         hdr_den, pawrhoij1_i1pert, comm_cell, check_hdr=hdr)
         call hdr_den%free()
       end if

       xccc3d1(:) = zero
       if (psps%usepaw==1 .or. psps%nc_xccc_gspace==1) then
         ndir=1
         call dfpt_atm2fft(atindx,cplex,gmet,gprimd,gsqcut,i1dir,i1pert,&
&         mgfftf,psps%mqgrid_vl,dtset%natom,ndir,nfftf,ngfftf,psps%ntypat,&
&         ph1df,psps%qgrid_vl,dtset%qptn,dtset%typat,ucvol,psps%usepaw,xred,psps,pawtab,&
&         atmrhor1=xccc3d1,optn_in=n3xccc/nfftf,optn2_in=1,optv_in=0,vspl=psps%vlspl)
       else
    !    Norm-conserving psp: compute Vloc(1) in reciprocal sp. and core(1) in real sp.
    !    ------------------------------------------------------------------------------
         if(psps%n1xccc/=0)then
           call dfpt_mkcore(cplex,i1dir,i1pert,dtset%natom,psps%ntypat,n1,psps%n1xccc,&
&           n2,n3,dtset%qptn,rprimd,dtset%typat,ucvol,psps%xcccrc,psps%xccc1d,xccc3d1,xred)
         end if ! psps%n1xccc/=0
       end if ! usepaw

       do i3pert = 1, mpert
         do i3dir = 1, 3

           if ((maxval(rfpert(i1dir,i1pert,:,:,i3dir,i3pert))==1)) then

             pert3case = i3dir + (i3pert-1)*3
             counter = 100*pert3case + pert1case
             call appdig(pert3case,dtfil%fnamewff1,fiwf3i)

             call inwffil(ask_accurate,cg3,dtset,dtset%ecut,ecut_eff,eigen3,dtset%exchn2n3d,&
&             formeig,hdr,ireadwf,dtset%istwfk,kg,dtset%kptns,dtset%localrdwf,&
&             dtset%mband,mcg,dtset%mk1mem,mpi_enreg,mpw,&
&             dtset%nband,dtset%ngfft,dtset%nkpt,npwarr,&
&             dtset%nsppol,dtset%nsym,&
&             occ,optorth,dtset%symafm,dtset%symrel,dtset%tnons,&
&             dtfil%unkg1,wff3,wfft3,dtfil%unwff3,&
&             fiwf3i,wvl)
             if (ireadwf==1) then
               call WffClose (wff3,ierr)
             end if

             flag3 = 0
             rho3r1(:,:) = zero
             if (dtset%get1den /= 0 .or. dtset%ird1den /= 0) then

               call appdig(pert3case,dtfil%fildens1in,fiden1i)

               call read_rhor(fiden1i, cplex, dtset%nspden, nfftf, ngfftf, rdwrpaw, mpi_enreg, rho3r1, &
               hdr_den, pawrhoij1_i3pert, comm_cell, check_hdr=hdr)
               call hdr_den%free()
             end if

             xccc3d3(:) = zero
             if (psps%usepaw==1 .or. psps%nc_xccc_gspace==1) then
               ndir=1
               call dfpt_atm2fft(atindx,cplex,gmet,gprimd,gsqcut,i3dir,i3pert,&
&               mgfftf,psps%mqgrid_vl,dtset%natom,ndir,nfftf,ngfftf,psps%ntypat,&
&               ph1df,psps%qgrid_vl,dtset%qptn,dtset%typat,ucvol,psps%usepaw,xred,psps,pawtab,&
&               atmrhor1=xccc3d3,optn_in=n3xccc/nfftf,optn2_in=1,optv_in=0,vspl=psps%vlspl)
             else
            !    Norm-conserving psp: compute Vloc(1) in reciprocal sp. and core(1) in real sp.
            !    ------------------------------------------------------------------------------
               if(psps%n1xccc/=0)then
                 call dfpt_mkcore(cplex,i3dir,i3pert,dtset%natom,psps%ntypat,n1,psps%n1xccc,&
&                 n2,n3,dtset%qptn,rprimd,dtset%typat,ucvol,psps%xcccrc,psps%xccc1d,xccc3d3,xred)
               end if ! psps%n1xccc/=0
             end if ! usepaw

             do i2pert = 1, mpert
               do i2dir = 1, 3

                 if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then

                   blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1

                   npert_phon = 0
                   if(i1pert<=dtset%natom) npert_phon = npert_phon + 1
                   if(i2pert<=dtset%natom) npert_phon = npert_phon + 1
                   if(i3pert<=dtset%natom) npert_phon = npert_phon + 1
                   if (npert_phon>1) then
                     MSG_ERROR("dfptnl_loop is available with at most one phonon perturbation. Change your input!")
                   end if

                   pert2case = i2dir + (i2pert-1)*3
                   counter = 100*pert2case + pert2case
                   call appdig(pert2case,dtfil%fnamewff1,fiwf2i)

                   call inwffil(ask_accurate,cg2,dtset,dtset%ecut,ecut_eff,eigen2,dtset%exchn2n3d,&
&                   formeig,hdr,ireadwf,dtset%istwfk,kg,dtset%kptns,dtset%localrdwf,&
&                   dtset%mband,mcg,dtset%mk1mem,mpi_enreg,mpw,&
&                   dtset%nband,dtset%ngfft,dtset%nkpt,npwarr,&
&                   dtset%nsppol,dtset%nsym,&
&                   occ,optorth,dtset%symafm,dtset%symrel,dtset%tnons,&
&                   dtfil%unkg1,wff2,wfft2,dtfil%unwff2,&
&                   fiwf2i,wvl)
                   if (ireadwf==1) then
                     call WffClose (wff2,ierr)
                   end if

!                  Read the first-order densities from disk-files
                   rho2r1(:,:) = zero ; rho2g1(:,:) = zero

                   if (dtset%get1den /= 0 .or. dtset%ird1den /= 0) then

                     call appdig(pert2case,dtfil%fildens1in,fiden1i)

                     call read_rhor(fiden1i, cplex, dtset%nspden, nfftf, ngfftf, rdwrpaw, mpi_enreg, rho2r1, &
                     hdr_den, pawrhoij1_i2pert , comm_cell, check_hdr=hdr)
                     call hdr_den%free()

!                    Compute up+down rho1(G) by fft
                     ABI_ALLOCATE(work,(cplex*nfftf))
                     work(:)=rho2r1(:,1)
                     call fourdp(cplex,rho2g1,work,-1,mpi_enreg,nfftf,1,ngfftf,0)
                     ABI_DEALLOCATE(work)

                   end if

                   xccc3d2(:)=zero ; vpsp1(:)=zero
                   !  PAW: compute Vloc(1) and core(1) together in reciprocal space
                   !  --------------------------------------------------------------
                   if (psps%usepaw==1 .or. psps%nc_xccc_gspace==1) then
                     ndir=1
                     call dfpt_atm2fft(atindx,cplex,gmet,gprimd,gsqcut,i2dir,i2pert,&
&                     mgfftf,psps%mqgrid_vl,dtset%natom,ndir,nfftf,ngfftf,psps%ntypat,&
&                     ph1df,psps%qgrid_vl,dtset%qptn,dtset%typat,ucvol,psps%usepaw,xred,psps,pawtab,&
&                     atmrhor1=xccc3d2,atmvlocr1=vpsp1,optn_in=n3xccc/nfftf,optn2_in=1,vspl=psps%vlspl)
                     !    PAW only: we sometimes have to compute 1st-order compensation density
                     !    and eventually add it to density from 1st-order WFs
                     !    ----------------------------------------------------------------------
                     if (psps%usepaw==1) then

                       !Force the computation of nhatfr
                       do iatom=1,dtset%natom
                         pawfgrtab(iatom)%nhatfr_allocated = 0
                         pawfgrtab(iatom)%nhatfr = zero
                       end do

!                      This portion of code works only when npert_phon<=1
                       if (i1pert<=natom.and.usexcnhat==0) then
                         call pawnhatfr(0,i1dir,i1pert,1,dtset%natom,nspden,psps%ntypat,&
&                         pawang,pawfgrtab(i1pert),pawrhoij(i1pert),pawtab,rprimd,&
&                         mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
                       end if
                       if (i2pert<=natom) then
                         call pawnhatfr(0,i2dir,i2pert,1,dtset%natom,nspden,psps%ntypat,&
&                         pawang,pawfgrtab(i2pert),pawrhoij(i2pert),pawtab,rprimd,&
&                         mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
                       end if
                       if (i3pert<=natom.and.usexcnhat==0) then
                         call pawnhatfr(0,i3dir,i3pert,1,dtset%natom,nspden,psps%ntypat,&
&                         pawang,pawfgrtab(i3pert),pawrhoij(i3pert),pawtab,rprimd,&
&                         mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
                       end if

                       if (usexcnhat==0) then

                         call pawmknhat(dummy_real,cplex,0,i1dir,i1pert,0,gprimd,natom,dtset%natom,&
&                         nfftf,ngfftf,nhat1grdim,nspden,psps%ntypat,pawang,pawfgrtab,nhat1gr,nhat1_i1pert,&
&                         pawrhoij1_i1pert,pawrhoij,pawtab,qphon,rprimd,ucvol,dtset%usewvl,xred,&
&                         mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
                         if (flag1==0) then
                           rho1r1(:,:) = rho1r1(:,:) - nhat1_i1pert(:,:)
                           flag1 = 1
                         end if

                         call pawmknhat(dummy_real,cplex,0,i3dir,i3pert,0,gprimd,natom,dtset%natom,&
&                         nfftf,ngfftf,nhat1grdim,nspden,psps%ntypat,pawang,pawfgrtab,nhat1gr,nhat1_i3pert,&
&                         pawrhoij1_i3pert,pawrhoij,pawtab,qphon,rprimd,ucvol,dtset%usewvl,xred,&
&                         mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
                         if (flag3==0) then
                           rho3r1(:,:) = rho3r1(:,:) - nhat1_i3pert(:,:)
                           flag3 = 1
                         end if

                       end if

                       call pawmknhat(dummy_real,cplex,0,i2dir,i2pert,0,gprimd,natom,dtset%natom,&
&                       nfftf,ngfftf,nhat1grdim,nspden,psps%ntypat,pawang,pawfgrtab,nhat1gr,nhat1_i2pert,&
&                       pawrhoij1_i2pert,pawrhoij,pawtab,qphon,rprimd,ucvol,dtset%usewvl,xred,&
&                       mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)

                     end if

                   else

                  !    Norm-conserving psp: compute Vloc(1) in reciprocal sp. and core(1) in real sp.
                  !    ------------------------------------------------------------------------------
                     if(psps%n1xccc/=0)then
                       call dfpt_mkcore(cplex,i2dir,i2pert,dtset%natom,psps%ntypat,n1,psps%n1xccc,&
&                       n2,n3,dtset%qptn,rprimd,dtset%typat,ucvol,psps%xcccrc,psps%xccc1d,xccc3d2,xred)
                     end if ! psps%n1xccc/=0

                     call dfpt_vlocal(atindx,cplex,gmet,gsqcut,i2dir,i2pert,mpi_enreg,psps%mqgrid_vl,dtset%natom,&
&                     nattyp,nfftf,ngfftf,psps%ntypat,n1,n2,n3,ph1df,psps%qgrid_vl,&
&                     dtset%qptn,ucvol,psps%vlspl,vpsp1,xred)

                   end if ! usepaw

                   option=1;optene=0
                   call dfpt_rhotov(cplex,dummy_real,dummy_real,dummy_real,dummy_real,dummy_real,&
&                   gsqcut,i2dir,i2pert,dtset%ixc,kxc,mpi_enreg,dtset%natom,nfftf,ngfftf,nhat,&
&                   nhat1_i2pert,nhat1gr,nhat1grdim,nkxc,nspden,n3xccc,non_magnetic_xc,optene,option,&
&                   dtset%qptn,rhog,rho2g1,rhor,rho2r1,rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1_i2pert,&
&                   vpsp1,vresid_dum,dummy_real,vtrial1_i2pert,vxc,vxc1_i2pert,xccc3d2,dtset%ixcrot)

                   if (psps%usepaw==1.and.usexcnhat==0) then
                     rho2r1(:,:) = rho2r1(:,:) - nhat1_i2pert(:,:)
                   end if

                   if (psps%usepaw==1)then
                     call paw_an_reset_flags(paw_an1_i2pert) ! Force the recomputation of on-site potentials
                     call paw_ij_reset_flags(paw_ij1_i2pert,all=.true.) ! Force the recomputation of Dij
                     optfr=0
                     call pawdijfr(gprimd,i2dir,i2pert,natom,natom,nfftf,ngfftf,nspden,nsppol,&
&                     psps%ntypat,optfr,paw_ij1_i2pert,pawang,pawfgrtab,pawrad,pawtab,cplex,qphon,&
&                     rprimd,ucvol,vpsp1,vtrial,vxc,xred,&
&                     mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)

!                    Computation of "on-site" first-order potentials, first-order densities
                     option=1
                     call pawdenpot(dummy_real,dummy_real2,dummy_real3,i2pert,dtset%ixc,natom,dtset%natom,&
&                     nspden,psps%ntypat,dtset%nucdipmom,&
&                     0,option,paw_an1_i2pert,paw_an0,paw_ij1_i2pert,pawang,&
&                     dtset%pawprtvol,pawrad,pawrhoij1_i2pert,dtset%pawspnorb,pawtab,dtset%pawxcdev,&
&                     dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,ucvol,psps%znuclpsp, &
&                     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
                !    First-order Dij computation
!                     call timab(561,1,tsec)
                     if (has_dijfr>0) then
                       !vpsp1 contribution to Dij already stored in frozen part of Dij
                       ABI_ALLOCATE(vtrial1_tmp,(cplex*nfftf,nspden))
                       vtrial1_tmp=vtrial1_i2pert
                       do ii=1,min(nspden,2)
                         vtrial1_tmp(:,ii)=vtrial1_tmp(:,ii)-vpsp1(:)
                       end do
                     else
                       vtrial1_tmp => vtrial1_i2pert
                     end if
                     call pawdij(cplex,dtset%enunit,gprimd,i2pert,natom,dtset%natom,&
&                     nfftf,nfftotf,dtset%nspden,psps%ntypat,paw_an1_i2pert,paw_ij1_i2pert,pawang,&
&                     pawfgrtab,dtset%pawprtvol,pawrad,pawrhoij1_i2pert,dtset%pawspnorb,pawtab,&
&                     dtset%pawxcdev,qphon,dtset%spnorbscl,ucvol,dtset%charge,vtrial1_tmp,vxc1_i2pert,xred,&
&                     mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
                     if (has_dijfr>0) then
                       ABI_DEALLOCATE(vtrial1_tmp)
                     end if
                     call symdij(gprimd,indsy1,i2pert,natom,dtset%natom,nsym1,psps%ntypat,0,&
&                     paw_ij1_i2pert,pawang1,dtset%pawprtvol,pawtab,rprimd,symaf1,symrc1, &
&                     mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom,&
&                     qphon=qphon)
!                     call timab(561,2,tsec)

                   end if ! end usepaw section

                   nwffile = 1
                   file_index(1) = i2dir + 3*(i2pert-1)
                   fnamewff(1) = dtfil%fnamewff1

                   if (i2pert==natom+2) then

                     nwffile = 3
                     file_index(2) = i2dir+natom*3
                     fnamewff(2) = dtfil%fnamewffddk
!                    As npert_phon<=1 and i2pert==natom+2, i1pert or i3pert is necessarly equal to natom+2
                     if (i3pert==natom+2) then
                       second_idir = i3dir
                     else if (i1pert==natom+2) then
                       second_idir = i1dir
                     else
                       MSG_BUG(" i1pert or i3pert is supposed to be equal to natom+2, which is not the case here.")
                     end if
                     call rf2_getidir(i2dir,second_idir,idir_dkde)
                     file_index(3) = idir_dkde+9+(dtset%natom+6)*3
                     fnamewff(3) = dtfil%fnamewffdkde

                     if (npert_phon==1.and.psps%usepaw==1.and.second_idir/=i2dir) then
                       nwffile = 5
                       file_index(4) = second_idir+natom*3
                       fnamewff(4) = dtfil%fnamewffddk
                       call rf2_getidir(second_idir,i2dir,idir_dkde) ! i2dir and second_idir are reversed
                       file_index(5) = idir_dkde+9+(dtset%natom+6)*3
                       fnamewff(5) = dtfil%fnamewffdkde
                     end if

                   end if

                   do ii=1,nwffile
                     call appdig(file_index(ii),fnamewff(ii),fiwfddk)
                     ! Checking the existence of data file
                     if (.not. file_exists(fiwfddk)) then
                       ! Trick needed to run Abinit test suite in netcdf mode.
                       if (file_exists(nctk_ncify(fiwfddk))) then
                         write(message,"(3a)")"- File: ",trim(fiwfddk),&
                         " does not exist but found netcdf file with similar name."
                         call wrtout(std_out,message,'COLL')
                         fiwfddk = nctk_ncify(fiwfddk)
                       end if
                       if (.not. file_exists(fiwfddk)) then
                         MSG_ERROR('Missing file: '//TRIM(fiwfddk))
                       end if
                     end if
                     write(message,'(2a)')'-dfptnl_loop : read the wavefunctions from file: ',trim(fiwfddk)
                     call wrtout(std_out,message,'COLL')
                     call wrtout(ab_out,message,'COLL')
!                    Note that the unit number for these files is 50,51,52 or 53 (dtfil%unddk=50)
                     call wfk_open_read(ddk_f(ii),fiwfddk,1,dtset%iomode,dtfil%unddk+(ii-1),mpi_enreg%comm_cell)
                   end do

!                  Perform DFPT part of the 3dte calculation
                   call timab(513,1,tsec)
!                  NOTE : eigen2 equals zero here

                   call dfptnl_pert(atindx,cg,cg1,cg2,cg3,cplex,dtfil,dtset,d3etot,eigen0,gs_hamkq,k3xc,indsy1,i1dir,&
&                   i2dir,i3dir,i1pert,i2pert,i3pert,kg,mband,mgfft,mkmem,mk1mem,mpert,mpi_enreg,&
&                   mpsang,mpw,natom,nattyp,nfftf,nfftotf,ngfftf,nkpt,nk3xc,nspden,nspinor,nsppol,nsym1,npwarr,occ,&
&                   pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawtab,pawrhoij,pawrhoij1_i1pert,pawrhoij1_i2pert,pawrhoij1_i3pert,&
&                   paw_an0,paw_an1_i2pert,paw_ij1_i2pert,ph1d,psps,rho1r1,rho2r1,rho3r1,&
&                   rprimd,symaf1,symrc1,ucvol,vtrial,vhartr1_i2pert,vtrial1_i2pert,vxc1_i2pert,&
&                   ddk_f,xccc3d1,xccc3d2,xccc3d3,xred,&
&                   d3etot_1,d3etot_2,d3etot_3,d3etot_4,d3etot_5,d3etot_6,d3etot_7,d3etot_8,d3etot_9)
                   call timab(513,2,tsec)


!                  Eventually close the dot file
                   do ii=1,nwffile
                     call ddk_f(ii)%close()
                   end do

!                   if (psps%usepaw==1) then
!                     do ii=1,natom
!                       pawfgrtab(ii)%nhatfr = zero
!                     end do
!                   end if

                 end if   ! rfpert
               end do    ! i2dir
             end do     ! i2pert

           end if   ! rfpert
         end do    ! i3dir
       end do     ! i3pert

     end if   ! rfpert
   end do    ! i1dir
 end do     ! i1pert

!More memory cleaning
 call gs_hamkq%free()

 ABI_DEALLOCATE(cg1)
 ABI_DEALLOCATE(cg2)
 ABI_DEALLOCATE(cg3)
 ABI_DEALLOCATE(eigen1)
 ABI_DEALLOCATE(eigen2)
 ABI_DEALLOCATE(eigen3)
 ABI_DEALLOCATE(rho1r1)
 ABI_DEALLOCATE(rho2r1)
 ABI_DEALLOCATE(rho2g1)
 ABI_DEALLOCATE(rho3r1)
 ABI_DEALLOCATE(nhat1gr)
 ABI_DEALLOCATE(vresid_dum)
 ABI_DEALLOCATE(vtrial1_i2pert)
 ABI_DEALLOCATE(vxc1_i2pert)
 ABI_DEALLOCATE(vhartr1_i2pert)
 ABI_DEALLOCATE(vpsp1)
 ABI_DEALLOCATE(xccc3d1)
 ABI_DEALLOCATE(xccc3d2)
 ABI_DEALLOCATE(xccc3d3)
 if (psps%usepaw==1) then
   call pawrhoij_free(pawrhoij1_i1pert)
   call pawrhoij_free(pawrhoij1_i2pert)
   call pawrhoij_free(pawrhoij1_i3pert)
   ABI_DEALLOCATE(nhat1_i1pert)
   ABI_DEALLOCATE(nhat1_i2pert)
   ABI_DEALLOCATE(nhat1_i3pert)
   call paw_an_free(paw_an1_i2pert)
   call paw_ij_free(paw_ij1_i2pert)
   ABI_DATATYPE_DEALLOCATE(paw_an1_i2pert)
   ABI_DATATYPE_DEALLOCATE(paw_ij1_i2pert)
 end if
 ABI_DATATYPE_DEALLOCATE(pawrhoij1_i1pert)
 ABI_DATATYPE_DEALLOCATE(pawrhoij1_i2pert)
 ABI_DATATYPE_DEALLOCATE(pawrhoij1_i3pert)

 call timab(503,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfptnl_loop
!!***

end module m_dfptnl_loop
!!***
