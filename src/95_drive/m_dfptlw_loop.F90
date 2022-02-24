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
 use m_kg,          only : getcut,getph
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
!scalars
 integer :: ask_accurate,comm_cell,cplex,formeig
 integer :: i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,ierr,ireadwf,mcg,me,mpsang
 integer :: n1,n2,n3,nhat1grdim,nfftotf,nspden,n3xccc
 integer :: optorth,pawread,pert1case,pert2case,pert3case,timrev,usexcnhat 
 real(dp) :: boxcut,ecut,ecut_eff,gsqcut                                    
 logical :: non_magnetic_xc
 character(len=500) :: msg                   
 character(len=fnlen) :: fiden1i,fiwf1i,fiwf2i,fiwf3i,fiwfddk
 type(gs_hamiltonian_type) :: gs_hamkq
 type(wffile_type) :: wff1,wff2,wff3,wfft1,wfft2,wfft3
 type(wvl_data) :: wvl
 type(hdr_type) :: hdr_den
!arrays
 real(dp),allocatable :: cg1(:,:),cg2(:,:),cg3(:,:),eigen1(:),eigen2(:),eigen3(:)
 real(dp),allocatable :: nhat1(:,:),nhat1gr(:,:,:),ph1d(:,:)
 real(dp),allocatable :: rho1g1(:,:),rho1r1(:,:)
 real(dp),allocatable :: rho2g1(:,:),rho2r1(:,:)
 real(dp),allocatable :: vhartr1(:),vpsp1(:),vresid_dum(:,:)
 real(dp),allocatable,target :: vtrial1_i2pert(:,:),vtrial1_i1pert(:,:)
 real(dp),allocatable :: vxc1(:,:),xccc3d1(:)
 type(pawrhoij_type),allocatable :: pawrhoij_read(:)
 
 
! *************************************************************************

 DBG_ENTER("COLL")
 
!Init parallelism
 comm_cell=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt

!Various initializations
 timrev = 1 ! as q=0
 cplex = 2 - timrev
 nspden = dtset%nspden
 ecut=dtset%ecut
 ecut_eff = ecut*(dtset%dilatmx)**2
 mpsang = psps%mpsang
 optorth=1;if (psps%usepaw==1) optorth=0
 non_magnetic_xc=.true.

 ABI_MALLOC(cg1,(2,dtset%mpw*dtset%nspinor*mband*dtset%mk1mem*dtset%nsppol))
 ABI_MALLOC(cg2,(2,dtset%mpw*dtset%nspinor*mband*dtset%mk1mem*dtset%nsppol))
 ABI_MALLOC(cg3,(2,dtset%mpw*dtset%nspinor*mband*dtset%mk1mem*dtset%nsppol))
 ABI_MALLOC(eigen1,(2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_MALLOC(eigen2,(2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_MALLOC(eigen3,(2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_MALLOC(rho1r1,(cplex*nfftf,dtset%nspden))
 ABI_MALLOC(rho2r1,(cplex*nfftf,dtset%nspden))
 ABI_MALLOC(rho1g1,(2,nfftf))
 ABI_MALLOC(rho2g1,(2,nfftf))

 ask_accurate=1 ; formeig = 1 ; ireadwf = 1
 n1=ngfftf(1) ; n2=ngfftf(2) ; n3=ngfftf(3)
 nfftotf=n1*n2*n3

!Compute large sphere cut-off gsqcut
 call getcut(boxcut,ecut,gmet,gsqcut,dtset%iboxcut,std_out,dtset%qptn,dtset%ngfft)

!Generate the 1-dimensional phases
 ABI_MALLOC(ph1d,(2,3*(2*dtset%mgfft+1)*dtset%natom))
 call getph(atindx,dtset%natom,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),ph1d,xred)

!==== Initialize most of the Hamiltonian (and derivative) ====
!1) Allocate all arrays and initialize quantities that do not depend on k and spin.
!2) Perform the setup needed for the non-local factors:
!3) Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 call init_hamiltonian(gs_hamkq,psps,pawtab,dtset%nspinor,dtset%nsppol,nspden,dtset%natom,&
& dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,ph1d=ph1d,&
& use_gpu_cuda=dtset%use_gpu_cuda)

 ABI_MALLOC(vpsp1,(cplex*nfftf))
 ABI_MALLOC(xccc3d1,(cplex*nfftf))
 ABI_MALLOC(vhartr1,(cplex*nfftf))
 ABI_MALLOC(vxc1,(cplex*nfftf,dtset%nspden))
 ABI_MALLOC(vtrial1_i1pert,(cplex*nfftf,dtset%nspden))
 ABI_MALLOC(vtrial1_i2pert,(cplex*nfftf,dtset%nspden))
 ABI_MALLOC(vresid_dum,(0,0))

!This is necessary to deactivate paw options in the dfpt_rhotov routine
 ABI_MALLOC(pawrhoij_read,(0))
 usexcnhat=0
 n3xccc=0
 pawread=0
 nhat1grdim=0
 ABI_MALLOC(nhat1gr,(0,0,0))
 nhat1gr(:,:,:) = zero
 ABI_MALLOC(nhat1,(cplex*dtset%nfft,nspden))
 nhat1=zero

 mcg=mpw*nspinor*mband*mkmem*nsppol

 pert1case = 0 ; pert2case = 0 ; pert3case = 0

 do i1pert = 1, mpert
   do i1dir = 1, 3

     if ((maxval(rfpert(i1dir,i1pert,:,:,:,:))==1)) then

       pert1case = i1dir + (i1pert-1)*3
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

       rho1r1(:,:) = zero
       if (dtset%get1den /= 0 .or. dtset%ird1den /= 0) then
         call appdig(pert1case,dtfil%fildens1in,fiden1i)

         call read_rhor(fiden1i, cplex, dtset%nspden, nfftf, ngfftf, psps%usepaw, mpi_enreg, rho1r1, &
         hdr_den, pawrhoij_read, comm_cell, check_hdr=hdr)
         call hdr_den%free()
       end if
       xccc3d1(:) = zero

     end if   ! rfpert
   end do    ! i1dir
 end do     ! i1pert

 DBG_EXIT("COLL")

end subroutine dfptlw_loop
!!***

end module m_dfptlw_loop
!!***
