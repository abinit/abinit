!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_rhofermi
!! NAME
!! dfpt_rhofermi
!!
!! FUNCTION
!! This routine computes the fixed contribution to the first-order
!! Fermi energy for metallic occupation and Q=0, as well as the
!! Fermi level charge density needed to compute the remainder of the
!! first-order Fermi energy from the self-consistent local potential
!! at each step in the iteration process.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DRH, DCA, XG, GMR, AR, MB, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!  cgq(2,mpw1*nspinor*mband*mkqmem*nsppol)=pw coefficients of GS wavefunctions at k+q.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL; if 2, COMPLEX
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk>
!!  cprjq(natom,nspinor*mband*mkqmem*nsppol*usecprj)= wave functions at k+q
!!              projected with non-local projectors: cprjq=<p_i|Cnk+q>
!!  doccde_rbz(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy
!!  docckqde(mband*nkpt_rbz*nsppol)=derivative of occkq wrt the energy
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigenq(mband*nkpt_rbz*nsppol)=GS eigenvalues at k+q (hartree)
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  idir=direction of the perturbation
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  ipert=type of the perturbation
!!  irrzon1(nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points.
!!  mband=maximum number of bands
!!  mkmem =number of k points treated by this node (GS data)
!!  mkqmem =number of k+q points treatede by this node (GS data)
!!  mk1mem =number of k points treated by this node.
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  maximum dimension for q points in grids for nonlocal form factors
!!  natom=number of atoms in cell.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs - see comment in respfn.F90)
!!  ngfftf(1:18)=integer array with FFT box dimensions and other for the "fine" grid
!!  nhatfermi(nfft,nspden)=array for fermi-level compensation charge density (PAW only)
!!  nkpt_rbz=number of k points in the IBZ for this perturbation
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym1=number of symmetry elements in space group consistent with
!!    perturbation
!!  occkq(mband*nkpt_rbz*nsppol)=occupation number for each band (often 2)
!!   at each k+q point of the reduced Brillouin zone.
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band and k (usually 2)
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawang1 <type(pawang_type)>=pawang datastr. containing only symmetries preserving the perturbation
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid for the GS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  phnons1(2,nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional real space primitive translations
!!  symaf1(nsym1)=(anti)ferromagnetic part of symmetry operations
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  symrl1(3,3,nsym1)=3x3 matrices of the group symmetries
!!  ucvol=volume of the unit cell
!!  usecprj= 1 if cprj, cprjq, cprj1 arrays are stored in memory
!!  useylmgr1= 1 if ylmgr1 array is allocated
!!  vtrial(nfftf,nspden)=GS potential (Hartree).
!!  vxc(nfftf,nspden)=XC potential (Hartree).
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= spherical harmonics for each G and k+g point
!!  ylmgr1(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics at k+q
!!
!!
!! OUTPUT
!!  eigen1(2*mband*mband*nkpt_rbz*nsppol)=array for holding eigenvalues
!!   (hartree) - only digonal elements computed here
!!  fe1fixed=fixed contribution to the first-order Fermi energy
!!   (nonlocal and kinetic in the case of strain)
!!  nhatfermi(cplex*nfftf,nspden)=fermi-level compensation charge density (PAW only)
!!  rhorfermi(cplex*nfftf,nspden)=fermi-level electronic density
!!
!! NOTES
!!  This routine will NOT work with nspden==4:
!!    at least the use of fftpac should be modified.
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      destroy_hamiltonian,destroy_rf_hamiltonian,dfpt_wfkfermi,fftpac
!!      init_hamiltonian,init_rf_hamiltonian,kpgstr,load_k_hamiltonian
!!      load_k_rf_hamiltonian,load_kprime_hamiltonian,load_spin_hamiltonian
!!      load_spin_rf_hamiltonian,mkffnl,mkkin,mkkpg,occeig,paw_ij_free
!!      paw_ij_init,paw_ij_nullify,pawdijfr,pawmkrho,pawrhoij_alloc
!!      pawrhoij_free,pawrhoij_free_unpacked,pawrhoij_init_unpacked
!!      pawrhoij_mpisum_unpacked,status,symrhg,timab,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_rhofermi(cg,cgq,cplex,cprj,cprjq,&
& doccde_rbz,docckqde,dtfil,dtset,eigenq,eigen0,eigen1,fe1fixed,gmet,gprimd,idir,&
& indsy1,ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,mband,mkmem,mkqmem,mk1mem,mpi_enreg,&
& mpw,mpw1,my_natom,natom,nband_rbz,ncpgr,nfftf,ngfftf,nhatfermi,nkpt_rbz,npwarr,npwar1,nspden,&
& nsppol,nsym1,occkq,occ_rbz,paw_ij,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoijfermi,pawtab,&
& phnons1,ph1d,prtvol,psps,rhorfermi,rmet,rprimd,symaf1,symrc1,symrl1,&
& ucvol,usecprj,useylmgr1,vtrial,vxc,wtk_rbz,xred,ylm,ylm1,ylmgr1)


 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_hamiltonian
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_wfk

 use m_io_tools,    only : get_unit, iomode_from_fname
 use m_pawang,      only : pawang_type
 use m_pawrad,      only : pawrad_type
 use m_pawtab,      only : pawtab_type
 use m_paw_ij,      only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawfgrtab,   only : pawfgrtab_type
 use m_pawrhoij,    only : pawrhoij_type, pawrhoij_init_unpacked, pawrhoij_gather, &
&                          pawrhoij_alloc, pawrhoij_free, pawrhoij_nullify, &
&                          pawrhoij_free_unpacked, pawrhoij_mpisum_unpacked
 use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_get
 use m_pawdij,      only : pawdijfr
 use m_pawfgr,      only : pawfgr_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_rhofermi'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_61_occeig
 use interfaces_65_paw
 use interfaces_66_nonlocal
 use interfaces_67_common
 use interfaces_72_response, except_this_one => dfpt_rhofermi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mband,mk1mem,mkmem,mkqmem
 integer,intent(in) :: mpw,mpw1,my_natom,natom,ncpgr,nfftf,nkpt_rbz,nspden,nsppol,nsym1
 integer,intent(in) :: prtvol,usecprj,useylmgr1
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: fe1fixed
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(in) :: pawang,pawang1
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: indsy1(4,nsym1,natom)
 integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))
 integer,intent(in) :: istwfk_rbz(nkpt_rbz),kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem)
 integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol),ngfftf(18)
 integer,intent(in) :: npwar1(nkpt_rbz),npwarr(nkpt_rbz),symaf1(nsym1)
 integer,intent(in) :: symrc1(3,3,nsym1),symrl1(3,3,nsym1)
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cgq(2,mpw1*dtset%nspinor*mband*mkqmem*nsppol)
 real(dp),intent(in) :: doccde_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: docckqde(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigen0(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigenq(mband*nkpt_rbz*nsppol),gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
 real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol),occkq(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
 real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3),vtrial(nfftf,nspden),vxc(nfftf,nspden),wtk_rbz(nkpt_rbz)
 real(dp),intent(in) :: xred(3,natom),ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 real(dp),intent(out) :: eigen1(2*mband*mband*nkpt_rbz*nsppol)
 real(dp),intent(out) :: nhatfermi(:,:)
 real(dp),intent(out) :: rhorfermi(cplex*nfftf,nspden)
 type(pawcprj_type),intent(in) :: cprj (natom,dtset%nspinor*mband*mkmem *nsppol*usecprj)
 type(pawcprj_type),intent(in) :: cprjq(natom,dtset%nspinor*mband*mkqmem*nsppol*usecprj)
 type(paw_ij_type),intent(in) :: paw_ij(my_natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
 type(pawrhoij_type),target,intent(inout)::pawrhoijfermi(my_natom*psps%usepaw)!vz_i
 type(pawtab_type), intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=17
 integer :: bd2tot_index,bdtot_index,buffer_size,counter,cplex_rhoij
 integer :: dimffnl1,dimffnlk,iatom,iband,ibg,ibgq
 integer :: icg,icgq,ider,idir0,ierr,iexit,ii,ikg,ikg1,ikpt,ilm,ilmn,indx
 integer :: ispden,isppol,istr,istwf_k
 integer :: mbd2kpsp,mcgq,mcgq_disk,mcprjq,mcprjq_disk
 integer :: me,n1,n2,n3,n4,n5,n6,nband_k,nkpg,nkpg1,npw1_k,npw_k,nspden_rhoij
 integer :: optfr,spaceworld
 logical :: paral_atom,qne0
 real(dp) :: arg,fe1norm,invfe1norm,wtk_k
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
!arrays
 integer,allocatable :: kg1_k(:,:),kg_k(:,:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: kpoint(3),kpq(3),tsec(2)
 real(dp) :: ylmgr_dum(1,1,1)
 real(dp),allocatable :: buffer1(:),dkinpw(:),doccde_k(:)
 real(dp),allocatable :: doccde_kq(:),eig0_k(:),eig0_kq(:),eig1_k(:)
 real(dp),allocatable :: fe1fixed_k(:),fe1norm_k(:)
 real(dp),allocatable :: ffnl1(:,:,:,:),ffnlk(:,:,:,:)
 real(dp),allocatable :: kinpw1(:),kpg1_k(:,:),kpg_k(:,:)
 real(dp),allocatable :: occ_k(:),occ_kq(:),ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: rhoaug(:,:,:),rhogfermi(:,:),rhowfr(:,:)
 real(dp),allocatable :: rocceig(:,:),ylm1_k(:,:),ylm_k(:,:),ylmgr1_k(:,:,:)
 type(paw_ij_type),allocatable :: paw_ij1fr(:)
 type(pawrhoij_type),pointer :: pawrhoijfermi_unsym(:)

! *********************************************************************

 DBG_ENTER('COLL')

!Check arguments validity
 if (ipert>natom.and.ipert/=natom+3.and.ipert/=natom+4) then
   MSG_BUG('wrong ipert argument!')
 end if
 if (cplex/=1) then
   MSG_BUG('wrong cplex/=1 argument !')
 end if

!Keep track of total time spent in this routine
 call timab(121,1,tsec)
 call timab(124,1,tsec)

!Retrieve parallelism data
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt
 paral_atom=(my_natom/=dtset%natom)
 my_atmtab=>mpi_enreg%my_atmtab

!Initialize output variables
 fe1fixed=zero
 if (psps%usepaw==0) rhorfermi(:,:)=zero

!Initialisations/allocation of temporary variables
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)
 bdtot_index=0 ; bd2tot_index=0 ; ibg=0 ; ibgq=0 ; icg=0 ; icgq=0
 qne0=(dtset%qptn(1)**2+dtset%qptn(2)**2+dtset%qptn(3)**2>=tol14)
 mbd2kpsp=2*mband**2*nkpt_rbz*nsppol
 fe1norm=zero
 ABI_ALLOCATE(rhoaug,(cplex*n4,n5,n6))
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(kg1_k,(3,mpw1))
 if (psps%usepaw==1) then
   ABI_ALLOCATE(rhowfr,(cplex*dtset%nfft,dtset%nspden))
   rhowfr(:,:)=zero
 end if

 mcgq=mpw1*dtset%nspinor*mband*mkqmem*nsppol;mcgq_disk=0

!Prepare RF PAW files for reading and writing if mkmem, mkqmem or mk1mem==0
 if (psps%usepaw==1) then
   mcprjq=dtset%nspinor*mband*mkqmem*nsppol*usecprj;mcprjq_disk=0
 else
   mcprjq=0;mcprjq_disk=0
 end if

!Initialize most of the Hamiltonian (arrays and quantities that do not depend on k + nl form factors)
 call init_hamiltonian(gs_hamkq,psps,pawtab,dtset%nspinor,nsppol,nspden,natom,&
& dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,&
& paw_ij=paw_ij,usecprj=usecprj,ph1d=ph1d,use_gpu_cuda=dtset%use_gpu_cuda,&
& mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
 call init_rf_hamiltonian(cplex,gs_hamkq,ipert,rf_hamkq,paw_ij1=paw_ij1fr,&
& mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)

!PAW:has to compute frozen part of Dij^(1) (without Vpsp(1) contribution)
 if (psps%usepaw==1) then
   ABI_DATATYPE_ALLOCATE(paw_ij1fr,(my_natom))
   call paw_ij_nullify(paw_ij1fr)
   call paw_ij_init(paw_ij1fr,cplex,dtset%nspinor,dtset%nsppol,dtset%nspden,0,&
&   dtset%natom,dtset%ntypat,dtset%typat,pawtab,has_dijfr=1,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom )
   optfr=1
   ABI_ALLOCATE(buffer1,(0))
   call pawdijfr(cplex,gprimd,idir,ipert,my_natom,natom,nfftf,ngfftf,dtset%nspden,&
&   dtset%ntypat,optfr,paw_ij1fr,pawang,pawfgrtab,pawrad,pawtab,&
&   dtset%qptn,rprimd,ucvol,buffer1,vtrial,vxc,xred,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
   ABI_DEALLOCATE(buffer1)
 end if

!PAW:allocate memory for non-symetrized occupancies matrix at EFermi (pawrhoijfermi)
 pawrhoijfermi_unsym => pawrhoijfermi
 if (psps%usepaw==1) then
   if (paral_atom) then
     ABI_DATATYPE_ALLOCATE(pawrhoijfermi_unsym,(natom))
     cplex_rhoij=max(cplex,dtset%pawcpxocc)
     nspden_rhoij=dtset%nspden;if (dtset%pawspnorb>0.and.dtset%nspinor==2) nspden_rhoij=4
     call pawrhoij_alloc(pawrhoijfermi_unsym,cplex_rhoij,nspden_rhoij,dtset%nspinor,&
&     dtset%nsppol,dtset%typat,pawtab=pawtab,use_rhoijp=0,use_rhoij_=1)
   else
     call pawrhoij_init_unpacked(pawrhoijfermi_unsym)
   end if
 end if

!LOOP OVER SPINS
 do isppol=1,nsppol
   ikg=0;ikg1=0

!  Continue to initialize the Hamiltonian at k+q
   call load_spin_hamiltonian(gs_hamkq,isppol,paw_ij=paw_ij,&
&       mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
   call load_spin_rf_hamiltonian(rf_hamkq,gs_hamkq,isppol,paw_ij1=paw_ij1fr,&
&       mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)

!  Nullify contribution to density at EFermi from this k-point
   rhoaug(:,:,:)=zero

   call timab(125,1,tsec)

!  BIG FAT k POINT LOOP
   do ikpt=1,nkpt_rbz

     counter=100*ikpt+isppol
     call status(counter,dtfil%filstat,iexit,level,'loop ikpt     ')
     nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
     istwf_k=istwfk_rbz(ikpt)
     npw_k=npwarr(ikpt)
     npw1_k=npwar1(ikpt)
     wtk_k=wtk_rbz(ikpt)

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       eigen1(1+bd2tot_index : 2*nband_k**2+bd2tot_index) = zero
       bdtot_index=bdtot_index+nband_k
       bd2tot_index=bd2tot_index+2*nband_k**2
!      Skip the rest of the k-point loop
       cycle
     end if

     ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
     ABI_ALLOCATE(ylm1_k,(npw1_k,psps%mpsang*psps%mpsang*psps%useylm))
     ABI_ALLOCATE(ylmgr1_k,(npw1_k,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))

!    Continue to initialize the Hamiltonian at k+q
     kpoint(:)=kpt_rbz(:,ikpt)
     kpq(:)=kpoint(:)+dtset%qptn(1:3)

     ABI_ALLOCATE(doccde_k,(nband_k))
     ABI_ALLOCATE(doccde_kq,(nband_k))
     ABI_ALLOCATE(eig0_k,(nband_k))
     ABI_ALLOCATE(eig0_kq,(nband_k))
     ABI_ALLOCATE(eig1_k,(2*nband_k**2))
     ABI_ALLOCATE(fe1fixed_k,(nband_k))
     ABI_ALLOCATE(fe1norm_k,(nband_k))
     ABI_ALLOCATE(occ_k,(nband_k))
     ABI_ALLOCATE(occ_kq,(nband_k))
     ABI_ALLOCATE(rocceig,(nband_k,nband_k))

     eig1_k(:)=zero
     eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
     eig0_kq(:)=eigenq(1+bdtot_index:nband_k+bdtot_index)
     occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)
     occ_kq(:)=occkq(1+bdtot_index:nband_k+bdtot_index)
     doccde_k(:)=doccde_rbz(1+bdtot_index:nband_k+bdtot_index)
     doccde_kq(:)=docckqde(1+bdtot_index:nband_k+bdtot_index)

!    For each pair of active bands (m,n), generates the ratios
!    rocceig(m,n)=(occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n))
!    and decide to which band to attribute it.
     call occeig(doccde_k,doccde_kq,eig0_k,eig0_kq,nband_k,&
&     dtset%occopt,occ_k,occ_kq,rocceig)

!    Get plane-wave coeffs and related data at k
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if

!    Get plane-wave coeffs and related data at k+q
     kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm1_k(1:npw1_k,ilm)=ylm1(1+ikg1:npw1_k+ikg1,ilm)
       end do
       if (useylmgr1==1) then
         do ilm=1,psps%mpsang*psps%mpsang
           do ii=1,3
             ylmgr1_k(1:npw1_k,ii,ilm)=ylmgr1(1+ikg1:npw1_k+ikg1,ii,ilm)
           end do
         end do
       end if
     end if

!    Set up the ground-state Hamiltonian, and some parts of the 1st-order Hamiltonian

!    Compute (k+G) vectors
     nkpg=0;if(ipert>=1.and.ipert<=natom) nkpg=3*dtset%nloalg(3)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
     end if

!    Compute (k+q+G) vectors
     nkpg1=0;if(ipert>=1.and.ipert<=natom) nkpg1=3*dtset%nloalg(3)
     ABI_ALLOCATE(kpg1_k,(npw1_k,nkpg1))
     if (nkpg1>0) then
       call mkkpg(kg1_k,kpg1_k,kpq,nkpg1,npw1_k)
     end if

!    ===== Preparation of non-local contributions

     dimffnlk=0;if (ipert<=natom) dimffnlk=1
     ABI_ALLOCATE(ffnlk,(npw_k,dimffnlk,psps%lmnmax,dtset%ntypat))

!    Compute nonlocal form factors ffnlk at (k+G)
     if (ipert<=natom ) then
       ider=0;idir0=0
       call status(counter,dtfil%filstat,iexit,level,'call mkffnl(0)')
       call mkffnl(psps%dimekb,dimffnlk,psps%ekb,ffnlk,psps%ffspl,&
&       gmet,gprimd,ider,idir0,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&       psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,dtset%ntypat,&
&       psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
     end if

!    Compute nonlocal form factors ffnl1 at (k+q+G)
     !-- Atomic displacement perturbation
     if (ipert<=natom) then
       ider=0;idir0=0
     !-- Strain perturbation
     else if (ipert==natom+3.or.ipert==natom+4) then
       if (ipert==natom+3) istr=idir
       if (ipert==natom+4) istr=idir+3
       ider=1;idir0=-istr
     end if
     dimffnl1=1+ider;if (ider==1.and.idir0==0) dimffnl1=dimffnl1+2*psps%useylm
     ABI_ALLOCATE(ffnl1,(npw1_k,dimffnl1,psps%lmnmax,dtset%ntypat))
     call status(counter,dtfil%filstat,iexit,level,'call mkffnl(1)')
     call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1,psps%ffspl,gmet,gprimd,ider,idir0,&
&     psps%indlmn,kg1_k,kpg1_k,kpq,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,&
&     npw1_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1_k,ylmgr1_k)

!    ===== Preparation of kinetic contributions

     ABI_ALLOCATE(dkinpw,(npw_k))
     ABI_ALLOCATE(kinpw1,(npw1_k))

!    Compute the derivative of the kinetic operator vs strain in dkinpw
     if (ipert==natom+3.or.ipert==natom+4) then
       if (ipert==natom+3) istr=idir
       if (ipert==natom+4) istr=idir+3
       call status(counter,dtfil%filstat,iexit,level,'call kpgstr   ')
       call kpgstr(dkinpw,dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,gprimd,istr,&
&       kg_k,kpoint,npw_k)
     end if

!    Compute (1/2) (2 Pi)**2 (k+q+G)**2:
     call status(counter,dtfil%filstat,iexit,level,'call mkkin(1) ')
!     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg1_k,kinpw1,kpq,npw1_k)
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg1_k,kinpw1,kpq,npw1_k,0,0)

!    ===== Load the k/k+q dependent parts of the Hamiltonian

!    Load k-dependent part in the Hamiltonian datastructure
     ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamkq%matblk))
     call load_k_hamiltonian(gs_hamkq,kpt_k=kpoint,npw_k=npw_k,istwf_k=istwf_k,kg_k=kg_k,kpg_k=kpg_k,&
&     ph3d_k=ph3d,compute_ph3d=.true.,compute_gbound=.true.)
     if (size(ffnlk)>0) then
       call load_k_hamiltonian(gs_hamkq,ffnl_k=ffnlk)
     else
       call load_k_hamiltonian(gs_hamkq,ffnl_k=ffnl1)
     end if

!    Load k+q-dependent part in the Hamiltonian datastructure
!        Note: istwf_k is imposed to 1 for RF calculations (should use istwf_kq instead)
     call load_kprime_hamiltonian(gs_hamkq,kpt_kp=kpq,npw_kp=npw1_k,istwf_kp=istwf_k,&
&     kinpw_kp=kinpw1,kg_kp=kg1_k,kpg_kp=kpg1_k,ffnl_kp=ffnl1,&
&     compute_gbound=.true.)
     if (qne0) then
       ABI_ALLOCATE(ph3d1,(2,npw1_k,gs_hamkq%matblk))
       call load_kprime_hamiltonian(gs_hamkq,ph3d_kp=ph3d1,compute_ph3d=.true.)
     end if

!    Load k-dependent part in the 1st-order Hamiltonian datastructure
     call load_k_rf_hamiltonian(rf_hamkq,npw_k=npw_k,dkinpw_k=dkinpw)

!    Compute fixed contributions to 1st-order Fermi energy
!    and Fermi level charge density
     call status(counter,dtfil%filstat,iexit,level,'call dfpt_wfkfermi  ')
     fe1fixed_k(:)=zero ; fe1norm_k(:)=zero

!    Note that dfpt_wfkfermi is called with kpoint, while kpt is used inside dfpt_wfkfermi
     call dfpt_wfkfermi(cg,cgq,cplex,cprj,cprjq,dtfil,eig0_k,eig1_k,fe1fixed_k,&
&     fe1norm_k,gs_hamkq,ibg,ibgq,icg,icgq,idir,ikpt,ipert,isppol,dtset%kptopt,mband,&
&     mcgq,mcprjq,mkmem,mpi_enreg,mpw,nband_k,ncpgr,npw_k,npw1_k,dtset%nspinor,nsppol,occ_k,&
&     pawrhoijfermi_unsym,prtvol,rf_hamkq,rhoaug,rocceig,wtk_k)

!    Free temporary storage
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(kpg1_k)
     ABI_DEALLOCATE(dkinpw)
     ABI_DEALLOCATE(ffnlk)
     ABI_DEALLOCATE(ffnl1)
     ABI_DEALLOCATE(kinpw1)
     ABI_DEALLOCATE(doccde_k)
     ABI_DEALLOCATE(doccde_kq)
     ABI_DEALLOCATE(eig0_k)
     ABI_DEALLOCATE(eig0_kq)
     ABI_DEALLOCATE(occ_kq)
     ABI_DEALLOCATE(rocceig)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylm1_k)
     ABI_DEALLOCATE(ylmgr1_k)
     ABI_DEALLOCATE(ph3d)
     if (allocated(ph3d1)) then
       ABI_DEALLOCATE(ph3d1)
     end if
     call status(counter,dtfil%filstat,iexit,level,'after dfpt_wfkfermi ')

!    Save eigenvalues (hartree)
     eigen1 (1+bd2tot_index : 2*nband_k**2+bd2tot_index) = eig1_k(:)

!    Accumulate sum over k points for 1st-order Fermi energy components
     do iband=1,nband_k
       fe1fixed=fe1fixed+wtk_k*occ_k(iband)*fe1fixed_k(iband)
       fe1norm=fe1norm+wtk_k*occ_k(iband)*fe1norm_k(iband)
     end do

     ABI_DEALLOCATE(eig1_k)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(fe1fixed_k)
     ABI_DEALLOCATE(fe1norm_k)

!    Keep track of total number of bands
!    (all k points so far, even for k points not treated by me)
     bdtot_index=bdtot_index+nband_k
     bd2tot_index=bd2tot_index+2*nband_k**2

!    Shift array memory
     if (mkmem/=0) then
       ibg=ibg+nband_k
       icg=icg+npw_k*dtset%nspinor*nband_k
       ikg=ikg+npw_k
     end if
     if (mkqmem/=0) then
       ibgq=ibgq+dtset%nspinor*nband_k
       icgq=icgq+npw1_k*dtset%nspinor*nband_k
     end if
     if (mk1mem/=0) then
       ikg1=ikg1+npw1_k
     end if

!    End big k point loop
   end do

   call timab(125,2,tsec)

!  Transfer density on augmented fft grid to normal fft grid in real space
!  Also take into account the spin.
   if (psps%usepaw==0) then
     call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,rhorfermi,rhoaug,1)
   else
     call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,rhowfr   ,rhoaug,1)
   end if
 end do ! End loop over spins

 !if(xmpi_paral==1)then
 !  call timab(166,1,tsec)
 !  call wrtout(std_out,'dfpt_rhofermi: loop on k-points and spins done in parallel','COLL')
 !  call xmpi_barrier(spaceworld)
 !  call timab(166,2,tsec)
 !end if

!More memory cleaning
 call destroy_hamiltonian(gs_hamkq)
 call destroy_rf_hamiltonian(rf_hamkq)
 if(psps%usepaw==1) then
   call paw_ij_free(paw_ij1fr)
   ABI_DATATYPE_DEALLOCATE(paw_ij1fr)
 end if
 ABI_DEALLOCATE(rhoaug)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(kg1_k)

 call status(0,dtfil%filstat,iexit,level,'after loops   ')
 call timab(124,2,tsec)


!=== MPI communications ==================
 if(xmpi_paral==1)then
   call timab(129,1,tsec)

!  Identify MPI buffer size
   buffer_size=cplex*dtset%nfft*nspden+2+mbd2kpsp
   ABI_ALLOCATE(buffer1,(buffer_size))

!  Pack rhorfermi, fe1fixed, fe1norm
   indx=cplex*dtset%nfft*nspden
   if (psps%usepaw==0) then
     buffer1(1:indx)=reshape(rhorfermi,(/indx/))
   else
     buffer1(1:indx)=reshape(rhowfr,(/indx/))
   end if
   buffer1(indx+1)=fe1fixed ; buffer1(indx+2)=fe1norm
   indx=indx+2 ; bd2tot_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt_rbz
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       buffer1(indx+1:indx+2*nband_k**2)=eigen1(bd2tot_index+1:bd2tot_index+2*nband_k**2)
       bd2tot_index=bd2tot_index+2*nband_k**2
       indx=indx+2*nband_k**2
     end do
   end do
   if(indx<buffer_size)buffer1(indx+1:buffer_size)=zero

!  Build sum of everything
   call timab(48,1,tsec)
   call xmpi_sum(buffer1,buffer_size,spaceworld,ierr)
   call timab(48,2,tsec)

!  Unpack the final result
   indx=cplex*dtset%nfft*nspden
   if (psps%usepaw==0) then
     rhorfermi(:,:)=reshape(buffer1(1:indx),(/cplex*dtset%nfft,nspden/))
   else
     rhowfr(:,:)=reshape(buffer1(1:indx),(/cplex*dtset%nfft,nspden/))
   end if
   fe1fixed=buffer1(indx+1) ; fe1norm =buffer1(indx+2)
   indx=indx+2 ; bd2tot_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt_rbz
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       eigen1(bd2tot_index+1:bd2tot_index+2*nband_k**2)=buffer1(indx+1:indx+2*nband_k**2)
       bd2tot_index=bd2tot_index+2*nband_k**2
       indx=indx+2*nband_k**2
     end do
   end do
   ABI_DEALLOCATE(buffer1)

!  Accumulate PAW occupancies
   if (psps%usepaw==1) then
     call pawrhoij_mpisum_unpacked(pawrhoijfermi_unsym,spaceworld)
   end if

   call timab(129,2,tsec)
 end if ! if kpt parallel
!=== End communications ==================

 call timab(127,1,tsec)

!Normalize the fixed part of fermie1
 invfe1norm = zero ; if (abs(fe1norm) > tol10) invfe1norm=one/fe1norm
 fe1fixed=fe1fixed*invfe1norm

!Symmetrize the density
 call status(0,dtfil%filstat,iexit,level,'call symrhg   ')
!In order to have the symrhg working in parallel on FFT coefficients, the size
!of irzzon1 and phnons1 should be set to nfftot. Therefore, nsym\=1 does not work.
!We also have the spin-up density, symmetrized, in rhorfermi(:,2).
 ABI_ALLOCATE(rhogfermi,(2,dtset%nfft))
 if (psps%usepaw==0) then
   call symrhg(cplex,gprimd,irrzon1,mpi_enreg,dtset%nfft,dtset%nfft,dtset%ngfft,nspden,&
&   nsppol,nsym1,dtset%paral_kgb,phnons1,rhogfermi,rhorfermi,rprimd,symaf1,symrl1)
 else
   call symrhg(cplex,gprimd,irrzon1,mpi_enreg,dtset%nfft,dtset%nfft,dtset%ngfft,nspden,&
&   nsppol,nsym1,dtset%paral_kgb,phnons1,rhogfermi,rhowfr   ,rprimd,symaf1,symrl1)
 end if

!PAW: Build new rhoij quantities then symetrize them
!Compute and add the compensation density to rhowfr to get the total density
 if (psps%usepaw == 1) then
   if (size(nhatfermi)>0) then
     call pawmkrho(arg,cplex,gprimd,0,indsy1,0,mpi_enreg,&
&     my_natom,natom,nspden,nsym1,dtset%ntypat,dtset%paral_kgb,pawang,pawfgr,&
&     pawfgrtab,-10001,pawrhoijfermi,pawrhoijfermi_unsym,pawtab,dtset%qptn,&
&     rhogfermi,rhowfr,rhorfermi,rprimd,symaf1,symrc1,dtset%typat,ucvol,&
&     dtset%usewvl,xred,pawang_sym=pawang1,pawnhat=nhatfermi)
   else
     call pawmkrho(arg,cplex,gprimd,0,indsy1,0,mpi_enreg,&
&     my_natom,natom,nspden,nsym1,dtset%ntypat,dtset%paral_kgb,pawang,pawfgr,&
&     pawfgrtab,-10001,pawrhoijfermi,pawrhoijfermi_unsym,pawtab,dtset%qptn,&
&     rhogfermi,rhowfr,rhorfermi,rprimd,symaf1,symrc1,dtset%typat,ucvol,&
&     dtset%usewvl,xred,pawang_sym=pawang1)
   end if
   ABI_DEALLOCATE(rhowfr)
   call pawrhoij_free_unpacked(pawrhoijfermi_unsym)
   if (paral_atom) then
     call pawrhoij_free(pawrhoijfermi_unsym)
     ABI_DATATYPE_DEALLOCATE(pawrhoijfermi_unsym)
   end if
 end if
 ABI_DEALLOCATE(rhogfermi)

!Normalize the Fermi level charge density (and associated PAW occupancies)
 rhorfermi(:,:)=invfe1norm*rhorfermi(:,:)
 if (psps%usepaw==1) then
   if (size(nhatfermi)>0) nhatfermi(:,:)=invfe1norm*nhatfermi(:,:)
   do iatom=1,my_natom
     do ispden=1,nspden
       do ilmn=1,pawrhoijfermi(iatom)%nrhoijsel
         pawrhoijfermi(iatom)%rhoijp(ilmn,ispden)=&
&         pawrhoijfermi(iatom)%rhoijp(ilmn,ispden)*invfe1norm
       end do
     end do
   end do
 end if

 call timab(127,2,tsec)
 call timab(121,2,tsec)

 DBG_EXIT('COLL')

end subroutine dfpt_rhofermi
!!***
