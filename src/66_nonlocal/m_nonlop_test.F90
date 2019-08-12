!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_nonlop_test
!! NAME
!!  m_nonlop_test
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2017-2019 ABINIT group (MT)
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

module m_nonlop_test

 implicit none

 private
!!***

 public :: nonlop_test
!!***

contains
!!***

!!****f* ABINIT/nonlop_test
!! NAME
!! nonlop_test
!!
!! FUNCTION
!! This routine is supposed to be used only for testing purpose.
!! It tests the "nonlop" routine application (non-local operator) with respect to Finite Differences.
!! It is not supposed to be used standalone, but via the nonlop_dfpt_test.py script to be found
!! in ~abinit/scripts/post_processing/nonlop_dfpt_test directory. This Python script
!! launches Abinit (several datasets) and analyse the result, in order to compare
!!  <Psi_i|H^(i)|Psi_j> compute with DFPT or Finite Differences.
!! H^(i) is the ith derivative of the Hamiltonian with respect to one or several perturbation(s).
!!
!! INPUTS
!!  cg(2,mcg)=wavefunctions (may be read from disk file)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced coordinates (integers) of G vecs in basis
!!  kpt(3,nkpt)=k points in reduced coordinates
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=number of k points treated by this node.
!!  mpi_enreg=informations about MPI parallelization
!!  mpw= maximum number of plane waves
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nband(nkpt)=number of bands at each k point
!!  nfft=number of FFT grid points
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt=number of k points in Brillouin zone
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt)=number of planewaves in basis and boundary at each k
!!  nspden=Number of spin Density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms
!!  paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  typat(natom)=type of each atom
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!
!! PARENTS
!!      afterscfloop
!!
!! CHILDREN
!!      destroy_hamiltonian,dotprod_g,init_hamiltonian,initylmg
!!      load_k_hamiltonian,load_spin_hamiltonian,mkffnl,mkkpg,nonlop
!!
!! SOURCE

subroutine nonlop_test(cg,eigen,istwfk,kg,kpt,mband,mcg,mgfft,mkmem,mpi_enreg,mpw,my_natom,natom,&
&                      nband,nfft,ngfft,nkpt,nloalg,npwarr,nspden,nspinor,nsppol,ntypat,&
&                       paw_ij,pawtab,ph1d,psps,rprimd,typat,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_abicore
 use m_xmpi
 use m_errors
 use m_hamiltonian
 use m_pawtab
 use m_paw_ij
 use m_pawcprj
 use m_cgtools

 use m_kg,             only : mkkpg
 use m_initylmg,       only : initylmg
 use m_mkffnl,         only : mkffnl
 use m_mpinfo,         only : proc_distrb_cycle
 use m_nonlop,         only : nonlop

!Arguments ------------------------------------
!scalars
 integer :: mband,mcg,mgfft,mkmem,mpw,my_natom,natom,nfft,nkpt,nspden,nspinor,nsppol,ntypat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: istwfk(nkpt),kg(3,mpw*mkmem),nband(nkpt*nsppol),ngfft(18),nloalg(3),npwarr(nkpt),typat(natom)
 real(dp),intent(in) :: cg(2,mcg),eigen(mband*nkpt*nsppol),kpt(3,nkpt), ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
 type(paw_ij_type),intent(in) :: paw_ij(my_natom*psps%usepaw)
!Local variables-------------------------------
!scalars
 integer,parameter :: ndtset_test=6,tim_nonlop=4
 integer,save :: idtset=0
 integer :: bandpp,bdtot_index,blocksize,choice,cplex,cpopt,dimffnl,iatm,iatom,iatom_only
 integer :: iband,iband_last,iband_test,iblock,icg,ider_ffnl,idir,idir_ffnl,idir_nonlop
 integer :: ii,ikg,ikpt,ilm,inlout,isppol,istwf_k,me_distrb,my_nspinor,nband_k,nblockbd
 integer :: nkpg,nnlout,npw_k,paw_opt,signs,spaceComm
 logical :: ex
 character(len=100) :: strg
 real(dp) :: argr,argi
 type(gs_hamiltonian_type) :: gs_hamk
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: kpoint(3),rmet(3,3)
 real(dp),allocatable :: cwavef(:,:),cwavef_out(:,:),enl(:,:,:,:),enlout(:),kpg_k(:,:),lambda(:)
 real(dp),allocatable :: scwavef_out(:,:),ylm(:,:),ylmgr(:,:,:),ylm_k(:,:),ylmgr_k(:,:,:)
 real(dp),allocatable,target :: ffnl(:,:,:,:),ph3d(:,:,:)
 type(pawcprj_type) :: cwaveprj(1,1)

!*************************************************************************

!Increment dataset counter
 idtset=idtset+1
 if (idtset<=2) return

!Data from parallelism
 spaceComm=mpi_enreg%comm_kpt
 me_distrb=mpi_enreg%me_kpt
 my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)

!Initialize Hamiltonian datastructure
 call init_hamiltonian(gs_hamk,psps,pawtab,nspinor,nsppol,nspden,natom,&
& typat,xred,nfft,mgfft,ngfft,rprimd,nloalg,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
& mpi_spintab=mpi_enreg%my_isppoltab,paw_ij=paw_ij,ph1d=ph1d)
 rmet = MATMUL(TRANSPOSE(rprimd),rprimd)

!Check for existence of files in the current directory\
! and set parameters accordingly
 choice=1 ; idir=0 ; signs=1
 if(idtset<ndtset_test)then
   inquire(file='config/signs1',exist=ex) ; if(ex) signs=1
   inquire(file='config/signs2',exist=ex) ; if(ex) signs=2
   do ii=1,100
     if (ii< 10) write(unit=strg,fmt='(a13,i1)') "config/choice",ii
     if (ii>=10) write(unit=strg,fmt='(a13,i2)') "config/choice",ii
     inquire(file=trim(strg),exist=ex)  ; if(ex) choice=ii
   end do
   do ii=1,9
     write(unit=strg,fmt='(a11,i1)') "config/idir",ii
     inquire(file=trim(strg),exist=ex)  ; if(ex) idir=ii
   end do
 else
   inquire(file='config/signsdfpt1',exist=ex)  ; if(ex)signs=1
   inquire(file='config/signsdfpt2',exist=ex)  ; if(ex)signs=2
   do ii=1,100
     if (ii< 10) write(unit=strg,fmt='(a17,i1)') "config/choicedfpt",ii
     if (ii>=10) write(unit=strg,fmt='(a17,i2)') "config/choicedfpt",ii
     inquire(file=trim(strg),exist=ex)  ; if(ex) choice=ii
   end do
   do ii=1,36
     if (ii< 10) write(unit=strg,fmt='(a15,i1)') "config/idirdfpt",ii
     if (ii>=10) write(unit=strg,fmt='(a15,i2)') "config/idirdfpt",ii
     inquire(file=trim(strg),exist=ex) ; if(ex) idir=ii
   end do
 end if
 iatom=1 ; iband_test=-1
 do ii=1,50
   if (ii< 10) write(unit=strg,fmt='(a12,i1)') "config/iatom",ii
   if (ii>=10) write(unit=strg,fmt='(a12,i2)') "config/iatom",ii
   inquire(file=trim(strg),exist=ex)  ; if(ex) iatom=ii
   if (ii< 10) write(unit=strg,fmt='(a12,i1)') "config/iband",ii
   if (ii>=10) write(unit=strg,fmt='(a12,i2)') "config/iband",ii
   inquire(file=trim(strg),exist=ex)  ; if(ex) iband_test=ii
 end do

!Set parameters for the "nonlop" routine according to users choice
 cpopt=-1 ; paw_opt=3*psps%usepaw ; iatm=gs_hamk%atindx(iatom)
 inquire(file='config/dij',exist=ex);if(ex) paw_opt=1*psps%usepaw
 if(signs==1)then
   iatom_only=-1 ; idir_ffnl=0
   idir_nonlop=0 ; cplex=1
   if(choice==1)then
     ider_ffnl=0
     nnlout=1 ; inlout=1
   end if
   if(choice==2)then
     ider_ffnl=0
     nnlout=3*natom ; inlout=3*(iatm-1)+idir ! Atoms are type-sorted in enlout()
   end if
   if(choice==3)then
     ider_ffnl=1
     nnlout=6 ; inlout=idir
   end if
   if(choice==5)then
     ider_ffnl=1
     nnlout=3 ; inlout=idir
   end if
   if(choice==51.or.choice==52)then
     ider_ffnl=1 ; cplex=2
     nnlout=6 ; inlout=2*idir-1
   end if
   if(choice==54)then
     ider_ffnl=2 ; cplex=2
     nnlout=18*natom ; inlout=18*(iatm-1)+2*idir-1 ! Atoms are type-sorted in enlout()
   end if
   if(choice==55)then
     ider_ffnl=2 ; cplex=2
     nnlout=36 ; inlout=2*idir-1
   end if
   if(choice==8)then
     ider_ffnl=2
     nnlout=6 ; inlout=idir
   end if
   if(choice==81)then
     ider_ffnl=2 ; cplex=2
     nnlout=18 ; inlout=2*idir-1
   end if
 else if(signs==2)then
   nnlout=1 ; inlout =1 ; cplex=1
   idir_nonlop=idir ; iatom_only=-1
   if(choice==1)then
     ider_ffnl=0 ; idir_ffnl=0
   end if
   if(choice==2)then
     iatom_only=iatom
     ider_ffnl=0 ; idir_ffnl=0
   end if
   if(choice==3)then
     ider_ffnl=1 ; idir_ffnl=-7
   end if
   if(choice==5)then
     ider_ffnl=1 ; idir_ffnl=4
   end if
   if(choice==51.or.choice==52)then
     ider_ffnl=1 ; idir_ffnl=4 ; cplex=2
   end if
   if(choice==54)then
     iatom_only=iatom
     ider_ffnl=2 ; idir_ffnl=4 ; cplex=2
   end if
   if(choice==8)then
     ider_ffnl=2 ; idir_ffnl=4
   end if
   if(choice==81)then
     ider_ffnl=2 ; idir_ffnl=4 ; cplex=2
   end if
 end if

!Set parameters for the "mkffnl" routine according to users choice
 dimffnl=1+ider_ffnl
 if (ider_ffnl==1.and.(idir_ffnl==0.or.idir_ffnl==4)) dimffnl=2+2*psps%useylm
 if (ider_ffnl==2.and.(idir_ffnl==0.or.idir_ffnl==4)) dimffnl=3+7*psps%useylm
 if (ider_ffnl==1.and.idir_ffnl==-7) dimffnl=2+5*psps%useylm
 if (idir_ffnl>-7.and.idir_ffnl<0) dimffnl=2

!Write recognizable statement in log file
 write(std_out,'(2(a,i2),(a,i1),2(a,i2),(a,i1),(a,i2),(a,i1),2(a,i2))') &
& "TESTDFPT: choice=",choice,", idir(mkffnl)=",idir_ffnl,&
& ", ider(mkffnl)=",ider_ffnl,", dimffnl=",dimffnl,&
& ", idir(nonlop)=",idir_nonlop,", signs=",signs,&
& ", iatom=",iatom_only,", paw_opt=",paw_opt,&
& ", nnlout=",nnlout,", inlout=",inlout

!Compute all spherical harmonics and gradients
 ABI_ALLOCATE(ylm,(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm))
 ABI_ALLOCATE(ylmgr,(mpw*mkmem,9,psps%mpsang*psps%mpsang*psps%useylm))
 if (psps%useylm==1) then
   call initylmg(gs_hamk%gprimd,kg,kpt,mkmem,mpi_enreg,psps%mpsang,mpw,nband,nkpt,&
&   npwarr,nsppol,2,rprimd,ylm,ylmgr)
 else
   ylm=zero ; ylmgr=zero
 end if

!No loop over spins; only do the first one
 bdtot_index=0 ; icg=0 ; isppol=1

!Continue to initialize the Hamiltonian (PAW DIJ coefficients)
 call load_spin_hamiltonian(gs_hamk,isppol,with_nonlocal=.true.)

!No loop over k points; only do the first one
 ikg=0 ; ikpt=1

 nband_k=nband(ikpt+(isppol-1)*nkpt)
 istwf_k=istwfk(ikpt)
 npw_k=npwarr(ikpt)
 kpoint(:)=kpt(:,ikpt)

!My spin/kpoint or not?
 if(.not.proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) then

!  Parallelism over FFT and/or bands: define sizes and tabs
   bandpp=mpi_enreg%bandpp
   nblockbd=nband_k/bandpp
   blocksize=nband_k/nblockbd

!  Several allocations
   ABI_ALLOCATE(lambda,(blocksize))
   ABI_ALLOCATE(enlout,(nnlout*blocksize))
   ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor*blocksize))
   ABI_ALLOCATE(cwavef_out,(2,npw_k))
   if (paw_opt>=3) then
     ABI_ALLOCATE(scwavef_out,(2,npw_k))
     ABI_ALLOCATE(enl,(0,0,0,0))
   else
     ABI_ALLOCATE(scwavef_out,(0,0))
     ABI_ALLOCATE(enl,(gs_hamk%dimekb1,gs_hamk%dimekb2,gs_hamk%nspinor**2,1))
     enl(:,:,:,:)=one
   end if

!  Compute (k+G) vectors and associated spherical harmonics
   nkpg=3*nloalg(3)
   ABI_ALLOCATE(kg_k,(3,mpw))
   ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
   ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
   ABI_ALLOCATE(ylmgr_k,(npw_k,9,psps%mpsang*psps%mpsang*psps%useylm))
   kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
   if (nkpg>0) then
     call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
   end if
   if (psps%useylm==1) then
     do ilm=1,psps%mpsang*psps%mpsang
       ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
     end do
     if (ider_ffnl>=1) then
       do ilm=1,psps%mpsang*psps%mpsang
         do ii=1,3+6*(ider_ffnl/2)
           ylmgr_k(1:npw_k,ii,ilm)=ylmgr(1+ikg:npw_k+ikg,ii,ilm)
         end do
       end do
     end if
   end if

!  Compute non-local form factors
   ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&   gs_hamk%gmet,gs_hamk%gprimd,ider_ffnl,idir_ffnl,psps%indlmn,kg_k,kpg_k,&
&   gs_hamk%kpt_k,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&   npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&   psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

!  Load k-dependent part in the Hamiltonian datastructure
   ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))
   call load_k_hamiltonian(gs_hamk,kpt_k=kpoint,istwf_k=istwf_k,npw_k=npw_k,&
&   kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl,ph3d_k=ph3d,compute_ph3d=.true.)

   do iblock=1,nblockbd

     iband=(iblock-1)*blocksize+1;iband_last=min(iband+blocksize-1,nband_k)
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband_last,isppol,me_distrb)) cycle

!    Select a specific band or all
     if (iband==iband_test.or.iband_test==-1) then

!      Load contribution from block(n,k)
       cwavef(:,1:npw_k*my_nspinor*blocksize)=&
&       cg(:,1+(iblock-1)*npw_k*my_nspinor*blocksize+icg:iblock*npw_k*my_nspinor*blocksize+icg)
       lambda(1:blocksize)= eigen(1+(iblock-1)*blocksize+bdtot_index:iblock*blocksize+bdtot_index)

!      Call NONLOP
       if (paw_opt<3) then
         call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,idir_nonlop,lambda,&
&         mpi_enreg,1,nnlout,paw_opt,signs,scwavef_out,tim_nonlop,cwavef,cwavef_out,&
&         iatom_only=iatom_only,enl=enl)
       else
         call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,idir_nonlop,lambda,&
&         mpi_enreg,1,nnlout,paw_opt,signs,scwavef_out,tim_nonlop,cwavef,cwavef_out,&
&         iatom_only=iatom_only)
       end if

!      Post-processing if nonlop is called with specific options
       if (signs==2) then
         if (paw_opt<3) then
           call dotprod_g(argr,argi,istwf_k,npw_k,cplex,cwavef,cwavef_out,&
&           mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
         else
           call dotprod_g(argr,argi,istwf_k,npw_k,cplex,cwavef,scwavef_out,&
&           mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
         end if
         enlout(inlout)=argr
       end if
       if (signs==1.and.choice==1) then
         call dotprod_g(argr,argi,istwf_k,npw_k,1,cwavef,cwavef,&
&         mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
         enlout(:)=enlout(:)+argr
       end if

!      Write recognizable statements in log file
       if (idtset<ndtset_test) then
         write(std_out,'(a,i3,es24.16)') "TESTDFPT_df:  ",iband,enlout(inlout)
       else
         write(std_out,'(a,i3,es24.16)') "TESTDFPT_dfpt:",iband,enlout(inlout)
       end if

     end if

   end do ! End of loop on block of bands

!  Increment indexes (not used here because only one spin/kpoint)
   icg=icg+npw_k*my_nspinor*nband_k
   ikg=ikg+npw_k

 end if ! Not my spin/kpoint

 bdtot_index=bdtot_index+nband_k

!Memory deallocations
 ABI_DEALLOCATE(enl)
 ABI_DEALLOCATE(enlout)
 ABI_DEALLOCATE(lambda)
 ABI_DEALLOCATE(ph3d)
 ABI_DEALLOCATE(ffnl)
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(cwavef_out)
 ABI_DEALLOCATE(scwavef_out)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(kpg_k)
 ABI_DEALLOCATE(ylm_k)
 ABI_DEALLOCATE(ylmgr_k)
 ABI_DEALLOCATE(ylm)
 ABI_DEALLOCATE(ylmgr)
 call destroy_hamiltonian(gs_hamk)

end subroutine nonlop_test
!!***

end module m_nonlop_test
!!***
