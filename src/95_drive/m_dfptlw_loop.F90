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
!!  d3e_pert1(mpert)=array with the i1pert cases to calculate
!!  d3e_pert2(mpert)=array with the i2pert cases to calculate
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
!!  rmet(3,3)=real space metric tensor in bohr**2
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

    
subroutine dfptlw_loop(atindx,blkflg,cg,d3e_pert1,d3e_pert2,d3etot,dtfil,dtset,eigen0,gmet,gprimd,&
& hdr,kg,kxc,mband,mgfft,mgfftf,mkmem,mk1mem,&
& mpert,mpi_enreg,mpw,natom,nattyp,ngfftf,nfftf,nhat,nkpt,nkxc,nspinor,nsppol,&
& npwarr,occ,&
& pawfgr,pawrad,pawrhoij,pawtab,&
& psps,rfpert,rhog,rhor,rmet,rprimd,ucvol,vxc,xred)

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
 use m_mklocl,      only : dfpt_vlocal, vlocalstr,dfpt_vlocaldq,dfpt_vmetdqdq 
 use m_dfptlw_pert, only : dfptlw_pert
 use m_dynmat,      only : cart39

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
 integer,intent(in) :: atindx(natom),d3e_pert1(mpert),d3e_pert2(mpert)
 integer,intent(in) :: kg(3,mk1mem*mpw)
 integer,intent(in) :: nattyp(psps%ntypat),ngfftf(18),npwarr(nkpt)
 integer,intent(in) :: rfpert(3,mpert,3,mpert,3,mpert)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert,3,mpert) 
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),gmet(3,3)
 real(dp),intent(in) :: eigen0(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: gprimd(3,3),kxc(nfftf,nkxc)
 real(dp),intent(in) :: nhat(nfftf,dtset%nspden)
 real(dp),intent(in) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(in) :: vxc(nfftf,dtset%nspden)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
 real(dp),intent(inout) :: d3etot(2,3,mpert,3,mpert,3,mpert) 
 type(pawrhoij_type),intent(in) :: pawrhoij(natom*psps%usepaw)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: ask_accurate,beta,comm_cell,cplex,delta,dkdk_index,formeig,gamma,g0term
 integer :: i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,idir_dkdk 
 integer :: idq,ierr,ii,ireadwf,istr,mcg,me,mpsang
 integer :: n1,n2,n3,n1dq,n2dq,nhat1grdim,nfftotf,nspden,n3xccc,optene
 integer :: opthartdqdq,optorth,optres,pawread
 integer :: pert1case,pert2case,pert3case,timrev,usexcnhat 
 real(dp) :: boxcut,dum_scl,ecut,ecut_eff,fac,gsqcut   
 logical :: non_magnetic_xc,samepert
 character(len=500) :: message
 character(len=fnlen) :: fiden1i,fiwf1i,fiwf2i,fiwf3i,fiwfddk,fiwfdkdk
 type(gs_hamiltonian_type) :: gs_hamkq
 type(wffile_type) :: wff1,wff2,wff3,wfft1,wfft2,wfft3
 type(wfk_t) :: ddk_f,d2_dkdk_f
 type(wvl_data) :: wvl
 type(hdr_type) :: hdr_den
!arrays
 integer,save :: idx(18)=(/1,1,2,2,3,3,3,2,3,1,2,1,2,3,1,3,1,2/)
 integer :: flg1(3),flg2(3)
 real(dp),allocatable :: cg1(:,:),cg2(:,:),cg3(:,:)
 real(dp),allocatable :: d3etot_t4(:,:),d3etot_t5(:,:),eigen1(:)
 real(dp),allocatable :: nhat1(:,:),nhat1gr(:,:,:),ph1d(:,:)
 real(dp),allocatable :: rho1g1(:,:),rho1r1(:,:)
 real(dp),allocatable :: rho2g1(:,:),rho2r1(:,:)
 real(dp),allocatable :: t4_typeI(:,:,:,:,:,:),t4_typeII(:,:,:,:,:,:,:)
 real(dp),allocatable :: t5_typeI(:,:,:,:,:,:),t5_typeII(:,:,:,:,:,:,:)
 real(dp),allocatable :: vhartr1(:),vhart1dqdq(:),vpsp1(:),vpsp1dqdq(:),vresid_dum(:,:)
 real(dp),allocatable :: vtrial1_i1pert(:,:),vtrial1_i2pert(:,:)
 real(dp),allocatable :: vpsp1_i1pertdq(:,:,:),vpsp1_i2pertdq(:,:,:)
 real(dp),allocatable :: vxc1(:,:),vxc1dqdq(:),work(:),xccc3d1(:)
 real(dp) :: vec1(3),vec2(3)
 type(pawrhoij_type),allocatable :: pawrhoij_read(:)
 
 
! *************************************************************************

 DBG_ENTER("COLL")

!Anounce start of spatial-dispersion calculation
 write(message, '(a,80a,a,a,a)' ) ch10,('=',ii=1,80),ch10,&
&   ' ==> Compute spatial-dispersion 3rd-order energy derivatives <== ',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 
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
 opthartdqdq=1

 ABI_MALLOC(cg1,(2,dtset%mpw*dtset%nspinor*mband*dtset%mk1mem*dtset%nsppol))
 ABI_MALLOC(cg2,(2,dtset%mpw*dtset%nspinor*mband*dtset%mk1mem*dtset%nsppol))
 ABI_MALLOC(eigen1,(2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_MALLOC(rho1r1,(cplex*nfftf,dtset%nspden))
 ABI_MALLOC(rho2r1,(cplex*nfftf,dtset%nspden))
 ABI_MALLOC(rho1g1,(2,nfftf))
 ABI_MALLOC(rho2g1,(2,nfftf))

 ask_accurate=1 ; formeig = 1 ; ireadwf = 1
 n1=ngfftf(1) ; n2=ngfftf(2) ; n3=ngfftf(3)
 nfftotf=n1*n2*n3

!Allocations for type-I terms
 ABI_MALLOC(t4_typeII,(2,3,mpert,3,mpert,3,mpert))
 if (d3e_pert2(natom+3)==1.or.d3e_pert2(natom+4)==1) then
   ABI_MALLOC(t4_typeI,(2,3,mpert,3,3,3))
 end if
 ABI_MALLOC(t5_typeII,(2,3,mpert,3,mpert,3,mpert))
 if (d3e_pert1(natom+3)==1.or.d3e_pert1(natom+4)==1) then
   ABI_MALLOC(t5_typeI,(2,3,mpert,3,3,3))
 end if

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
 xccc3d1(:)=zero

!Specific allocations for strain-gradient perturbation
 if (dtset%lw_flexo==1.or.dtset%lw_flexo==2.or.dtset%lw_flexo==4) then
  ABI_MALLOC(vhart1dqdq,(2*nfftf))
  ABI_MALLOC(vpsp1dqdq,(2*nfftf))
  ABI_MALLOC(vxc1dqdq,(2*nfftf))
 end if

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
       & formeig,hdr,ireadwf,dtset%istwfk,kg,dtset%kptns,dtset%localrdwf,&
       & dtset%mband,mcg,dtset%mk1mem,mpi_enreg,mpw,&
       & dtset%nband,dtset%ngfft,dtset%nkpt,npwarr,&
       & dtset%nsppol,dtset%nsym,&
       & occ,optorth,dtset%symafm,dtset%symrel,dtset%tnons,&
       & dtfil%unkg1,wff1,wfft1,dtfil%unwff1,fiwf1i,wvl)

       if (ireadwf==1) then
         call WffClose (wff1,ierr)
       end if

       rho1r1(:,:) = zero; rho1g1(:,:) = zero
       if (dtset%get1den /= 0 .or. dtset%ird1den /= 0) then
         call appdig(pert1case,dtfil%fildens1in,fiden1i)

         call read_rhor(fiden1i, cplex, dtset%nspden, nfftf, ngfftf, psps%usepaw, mpi_enreg, rho1r1, &
         hdr_den, pawrhoij_read, comm_cell, check_hdr=hdr)
         call hdr_den%free()
       end if

       !Perform FFT rhor1 to rhog1
       ABI_MALLOC(work,(cplex*nfftf))
       work(:)=rho1r1(:,1)
       call fourdp(cplex,rho1g1,work,-1,mpi_enreg,dtset%nfft,1,dtset%ngfft,0)
       ABI_FREE(work)

       !Calculate first-order local and SCF potential
       if (i1pert<=natom) then
         call dfpt_vlocal(atindx,cplex,gmet,gsqcut,i1dir,i1pert,mpi_enreg,psps%mqgrid_vl,dtset%natom,&
       & nattyp,dtset%nfft,dtset%ngfft,psps%ntypat,n1,n2,n3,ph1d,psps%qgrid_vl,&
       & dtset%qptn,ucvol,psps%vlspl,vpsp1,xred)
       !TODO: the strain perturbation is not used as ipert1 for the implemented
       !quantities.
       else if (i1pert==natom+3.or.i1pert==natom+4) then
         istr=i1dir; if (i1pert==natom+4) istr=3+i1dir
         g0term=1
         call vlocalstr(gmet,gprimd,gsqcut,istr,dtset%mgfft,mpi_enreg,&
       & psps%mqgrid_vl,dtset%natom,nattyp,dtset%nfft,dtset%ngfft,dtset%ntypat,ph1d,psps%qgrid_vl,&
       & ucvol,psps%vlspl,vpsp1,g0term=g0term)
       else
         vpsp1(:)=zero
       end if     

       optene=0; optres=1
       call dfpt_rhotov(cplex,dum_scl,dum_scl,dum_scl,dum_scl,dum_scl, &
       & gsqcut,i1dir,i1pert,dtset%ixc,kxc,mpi_enreg,dtset%natom,dtset%nfft,&
       & dtset%ngfft,nhat,nhat1,nhat1gr,nhat1grdim,nkxc,nspden,n3xccc,&
       & non_magnetic_xc,optene,optres,dtset%qptn,rhog,rho1g1,rhor,rho1r1,&
       & rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,vresid_dum,dum_scl,&
       & vtrial1_i1pert,vxc,vxc1,xccc3d1,dtset%ixcrot)

       !Allocate the first-order gradient local potential
       if (i1pert <= natom .or. i1pert <= natom+3) then
         n1dq=1
         ABI_MALLOC(vpsp1_i1pertdq,(2*nfftf,dtset%nspden,n1dq))
       else if (i1pert <= natom+4) then
         n1dq=2
         ABI_MALLOC(vpsp1_i1pertdq,(2*nfftf,dtset%nspden,n1dq))
       else
         n1dq=1
       end if
       ABI_MALLOC(d3etot_t5,(2,n1dq))

       do i2pert = 1, mpert
         do i2dir = 1, 3
       
           if ((maxval(rfpert(i1dir,i1pert,i2dir,i2pert,:,:))==1)) then
       
             pert2case = i2dir + (i2pert-1)*3
             call appdig(pert2case,dtfil%fnamewff1,fiwf2i)

             call inwffil(ask_accurate,cg2,dtset,dtset%ecut,ecut_eff,eigen1,dtset%exchn2n3d,&
             & formeig,hdr,ireadwf,dtset%istwfk,kg,dtset%kptns,dtset%localrdwf,&
             & dtset%mband,mcg,dtset%mk1mem,mpi_enreg,mpw,&
             & dtset%nband,dtset%ngfft,dtset%nkpt,npwarr,&
             & dtset%nsppol,dtset%nsym,&
             & occ,optorth,dtset%symafm,dtset%symrel,dtset%tnons,&
             & dtfil%unkg1,wff2,wfft2,dtfil%unwff2,fiwf2i,wvl)
 
             if (ireadwf==1) then
               call WffClose (wff2,ierr)
             end if

             if (i1pert==i2pert.and.i1dir==i2dir) then
               samepert=.true.
             else
               samepert=.false.
             end if

             rho2r1(:,:) = zero; rho2g1(:,:) = zero
             if (dtset%get1den /= 0 .or. dtset%ird1den /= 0) then
               call appdig(pert2case,dtfil%fildens1in,fiden1i)
      
               call read_rhor(fiden1i, cplex, dtset%nspden, nfftf, ngfftf, psps%usepaw, mpi_enreg, rho2r1, &
               hdr_den, pawrhoij_read, comm_cell, check_hdr=hdr)
               call hdr_den%free()
             end if

             if (.not.samepert) then
               !Perform FFT rhor1 to rhog1
               ABI_MALLOC(work,(cplex*nfftf))
               work(:)=rho2r1(:,1)
               call fourdp(cplex,rho2g1,work,-1,mpi_enreg,dtset%nfft,1,dtset%ngfft,0)
               ABI_FREE(work)
  
               !Calculate first-order local and SCF potential
               if (i2pert<=natom) then
                 call dfpt_vlocal(atindx,cplex,gmet,gsqcut,i2dir,i2pert,mpi_enreg,psps%mqgrid_vl,dtset%natom,&
               & nattyp,dtset%nfft,dtset%ngfft,psps%ntypat,n1,n2,n3,ph1d,psps%qgrid_vl,&
               & dtset%qptn,ucvol,psps%vlspl,vpsp1,xred)
               else if (i2pert==natom+3.or.i2pert==natom+4) then
                 istr=i2dir; if (i2pert==natom+4) istr=3+i2dir
                 g0term=1
                 call vlocalstr(gmet,gprimd,gsqcut,istr,dtset%mgfft,mpi_enreg,&
               & psps%mqgrid_vl,dtset%natom,nattyp,dtset%nfft,dtset%ngfft,dtset%ntypat,ph1d,psps%qgrid_vl,&
               & ucvol,psps%vlspl,vpsp1,g0term=g0term)
               !TODO: the electric field  perturbation is not used as ipert2 
               !for the implemented quantities.
               else
                 vpsp1(:)=zero
               end if     
  
               optene=0; optres=1
               call dfpt_rhotov(cplex,dum_scl,dum_scl,dum_scl,dum_scl,dum_scl, &
               & gsqcut,i2dir,i2pert,dtset%ixc,kxc,mpi_enreg,dtset%natom,dtset%nfft,&
               & dtset%ngfft,nhat,nhat1,nhat1gr,nhat1grdim,nkxc,nspden,n3xccc,&
               & non_magnetic_xc,optene,optres,dtset%qptn,rhog,rho2g1,rhor,rho2r1,&
               & rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,vresid_dum,dum_scl,&
               & vtrial1_i2pert,vxc,vxc1,xccc3d1,dtset%ixcrot)
             end if !samepert  

             !Allocate the first-order gradient local potential
             if (i2pert <= natom .or. i2pert <= natom+3) then
               n2dq=1
               ABI_MALLOC(vpsp1_i2pertdq,(2*nfftf,dtset%nspden,n2dq))
             else if (i2pert <= natom+4) then
               n2dq=2
               ABI_MALLOC(vpsp1_i2pertdq,(2*nfftf,dtset%nspden,n2dq))
             else
               n2dq=1
             end if
             ABI_MALLOC(d3etot_t4,(2,n2dq))

             do i3pert = 1, mpert
               do i3dir = 1, 3

                 if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then

                   blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1

                   !Calculate local potentials for first-order gradient Hamiltonians
                   !gradient of i1pert:
                   if (i1pert<=natom) then
                     !Get q-gradient of first-order local part of the pseudopotential
                     call dfpt_vlocaldq(atindx,2,gmet,gsqcut,i1dir,i1pert,mpi_enreg, &
                     & psps%mqgrid_vl,dtset%natom,nattyp,dtset%nfft,dtset%ngfft,dtset%ntypat,n1,n2,n3, &
                     & ph1d,i3dir,psps%qgrid_vl,dtset%qptn,ucvol,psps%vlspl,vpsp1_i1pertdq(:,1,1))
                   else if (i1pert==natom+3.or.i1pert==natom+4) then
                     istr=i1dir; if (i1pert==natom+4) istr=3+i1dir
                     !Get 2nd q-gradient of first-order local part of the pseudopotential and of the Hartree
                     !(and XC if GGA) contribution from ground state density
                     call dfpt_vmetdqdq(2,gmet,gprimd,gsqcut,istr,i1pert,kxc,mpi_enreg, &
                     & psps%mqgrid_vl,natom,nattyp,dtset%nfft,dtset%ngfft,dtset%ntypat,n1,n2,n3,&
                     & nkxc,nspden,opthartdqdq,ph1d,i3dir,psps%qgrid_vl,&
                     & dtset%qptn,rhog,rhor,ucvol,psps%vlspl,vhart1dqdq,vpsp1dqdq,vxc1dqdq)
                     vpsp1_i1pertdq(:,1,1)=vhart1dqdq(:)+vpsp1dqdq(:)+vxc1dqdq(:)
                     if (i1pert==natom+4) then
                       !Here we need to calculate both extradiagonal shear-strains
                       !because the second gradient of the metric perturbation is
                       !type-I, i.e., non symmetric with respect to the
                       !permutation of the strain indexes. 
                       istr=6+i1dir
                       call dfpt_vmetdqdq(2,gmet,gprimd,gsqcut,istr,i1pert,kxc,mpi_enreg, &
                       & psps%mqgrid_vl,natom,nattyp,dtset%nfft,dtset%ngfft,dtset%ntypat,n1,n2,n3,&
                       & nkxc,nspden,opthartdqdq,ph1d,i3dir,psps%qgrid_vl,&
                       & dtset%qptn,rhog,rhor,ucvol,psps%vlspl,vhart1dqdq,vpsp1dqdq,vxc1dqdq)
                       vpsp1_i1pertdq(:,1,2)=vhart1dqdq(:)+vpsp1dqdq(:)+vxc1dqdq(:)
                     end if
                   end if

                   if (.not.samepert) then
                     !gradient of i2pert:
                     if (i2pert<=natom) then
                       !Get q-gradient of first-order local part of the pseudopotential
                       call dfpt_vlocaldq(atindx,2,gmet,gsqcut,i2dir,i2pert,mpi_enreg, &
                       & psps%mqgrid_vl,dtset%natom,nattyp,dtset%nfft,dtset%ngfft,dtset%ntypat,n1,n2,n3, &
                       & ph1d,i3dir,psps%qgrid_vl,dtset%qptn,ucvol,psps%vlspl,vpsp1_i2pertdq(:,1,1))
                     else if (i2pert==natom+3.or.i2pert==natom+4) then
                       istr=i2dir; if (i2pert==natom+4) istr=3+i2dir
                       !Get 2nd q-gradient of first-order local part of the pseudopotential and of the Hartree
                       !(and XC if GGA) contribution from ground state density
                       call dfpt_vmetdqdq(2,gmet,gprimd,gsqcut,istr,i2pert,kxc,mpi_enreg, &
                       & psps%mqgrid_vl,natom,nattyp,dtset%nfft,dtset%ngfft,dtset%ntypat,n1,n2,n3,&
                       & nkxc,nspden,opthartdqdq,ph1d,i3dir,psps%qgrid_vl,&
                       & dtset%qptn,rhog,rhor,ucvol,psps%vlspl,vhart1dqdq,vpsp1dqdq,vxc1dqdq)
                       vpsp1_i2pertdq(:,1,1)=vhart1dqdq(:)+vpsp1dqdq(:)+vxc1dqdq(:)
                       if (i2pert==natom+4) then
                         !Here we need to calculate both extradiagonal shear-strains
                         !because the second gradient of the metric perturbation is
                         !type-I, i.e., non symmetric with respect to the
                         !permutation of the strain indexes. 
                         istr=6+i2dir
                         call dfpt_vmetdqdq(2,gmet,gprimd,gsqcut,istr,i2pert,kxc,mpi_enreg, &
                         & psps%mqgrid_vl,natom,nattyp,dtset%nfft,dtset%ngfft,dtset%ntypat,n1,n2,n3,&
                         & nkxc,nspden,opthartdqdq,ph1d,i3dir,psps%qgrid_vl,&
                         & dtset%qptn,rhog,rhor,ucvol,psps%vlspl,vhart1dqdq,vpsp1dqdq,vxc1dqdq)
                         vpsp1_i2pertdq(:,1,2)=vhart1dqdq(:)+vpsp1dqdq(:)+vxc1dqdq(:)
                       end if
                     end if
                   end if !samepert

                   !Prepare ddk wf file
                   pert3case = i3dir + natom*3
                   call appdig(pert3case,dtfil%fnamewffddk,fiwfddk)
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
                       ABI_ERROR('Missing file: '//TRIM(fiwfddk))
                     end if
                   end if
                   write(message,'(2a)')'-dfptlw_loop : read the ddk wavefunctions from file: ',trim(fiwfddk)
                   call wrtout(std_out,message,'COLL')
                   !call wrtout(ab_out,message,'COLL')
!                  Note that the unit number for these files is 50,51,52 or 53 (dtfil%unddk=50)
                   call wfk_open_read(ddk_f,fiwfddk,1,dtset%iomode,dtfil%unddk,mpi_enreg%comm_cell)

                   !Prepare d2_dkdk wf file
                   if (i1pert==natom+2) then
                     call rf2_getidir(i1dir,i3dir,idir_dkdk)
                     if (idir_dkdk>6) idir_dkdk=idir_dkdk-3
                     dkdk_index=idir_dkdk+(dtset%natom+6)*3
                     call appdig(dkdk_index,dtfil%fnamewffdkdk,fiwfdkdk)
                     !Check that d2_ddk file exists and open it
                     if (.not. file_exists(fiwfdkdk)) then
                       ! Trick needed to run Abinit test suite in netcdf mode. 
                       if (file_exists(nctk_ncify(fiwfdkdk))) then             
                         write(message,"(3a)")"- File: ",trim(fiwfdkdk),& 
                         " does not exist but found netcdf file with similar name."
                         call wrtout(std_out,message,'COLL')
                         fiwfdkdk = nctk_ncify(fiwfdkdk)
                       end if
                       if (.not. file_exists(fiwfdkdk)) then
                         ABI_ERROR('Missing file: '//TRIM(fiwfdkdk))
                       end if
                     end if
                     write(message,'(2a)')'-dfptlw_loop : read the d2_dkdk wavefunctions from file: ',trim(fiwfdkdk)
                     call wrtout(std_out,message,'COLL')
                     !call wrtout(ab_out,message,'COLL') 
                     call wfk_open_read(d2_dkdk_f,fiwfdkdk,1,dtset%iomode,dtfil%unddk+1,mpi_enreg%comm_cell)

                   end if

                   !Perform the longwave DFPT part of the 3dte calculation
                   call dfptlw_pert(atindx,cg,cg1,cg2,cplex,dtfil,dtset,d3etot,d3etot_t4,d3etot_t5,gs_hamkq,gsqcut,i1dir,&
                   & i2dir,i3dir,i1pert,i2pert,i3pert,kg,kxc,mband,mgfft,mkmem,mk1mem,mpert,mpi_enreg,&
                   & mpsang,mpw,natom,nattyp,n1dq,n2dq,nfftf,ngfftf,nkpt,nkxc,nspden,nspinor,nsppol,npwarr,occ,&
                   & pawfgr,ph1d,psps,rhog,rho1g1,rhor,rho2r1,rmet,rprimd,samepert,ucvol,vpsp1_i1pertdq,vpsp1_i2pertdq,&
                   & vtrial1_i1pert,vtrial1_i2pert,ddk_f,d2_dkdk_f,xccc3d1,xred)

                   !close ddk file
                   call ddk_f%close()

                   !close d2_dkdk file
                   if (i1pert==natom+2) call d2_dkdk_f%close()

                   !Save the type-I terms
                   if (i2pert==natom+3.or.i2pert==natom+4) then
                     gamma=i3dir
                     do idq=1,n2dq
                       if (i2pert==natom+3) then
                         istr=i2dir
                       else 
                         istr=idq*3+i2dir
                       endif 
                       beta=idx(2*istr-1); delta=idx(2*istr)
                       t4_typeI(:,i1dir,i1pert,beta,delta,gamma)=d3etot_t4(:,idq)
                     end do
                   else 
                     t4_typeII(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)=d3etot_t4(:,1)
                   end if
      
                   if (i1pert==natom+3.or.i1pert==natom+4) then
                     gamma=i3dir
                     do idq=1,n1dq
                       if (i1pert==natom+3) then
                         istr=i1dir
                       else 
                         istr=idq*3+i1dir
                       endif 
                       beta=idx(2*istr-1); delta=idx(2*istr)
                       t5_typeI(:,i2dir,i2pert,beta,delta,gamma)=d3etot_t5(:,idq)
                     end do
                   else 
                     t5_typeII(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)=d3etot_t5(:,1)
                   end if

                 end if   ! rfpert
               end do    ! ir3dir
             end do     ! ir3pert
             
             ABI_FREE(vpsp1_i2pertdq)
             ABI_FREE(d3etot_t4)

           end if   ! rfpert
         end do    ! i2dir
       end do     ! i2pert

       ABI_FREE(vpsp1_i1pertdq)
       ABI_FREE(d3etot_t5)

     end if   ! rfpert
   end do    ! i1dir
 end do     ! i1pert

!More memory cleaning
 call gs_hamkq%free()

 ABI_FREE(cg1)
 ABI_FREE(cg2)
 ABI_FREE(eigen1)
 ABI_FREE(rho1r1)
 ABI_FREE(rho2r1)
 ABI_FREE(rho1g1)
 ABI_FREE(rho2g1)
 ABI_FREE(nhat1gr)
 ABI_FREE(nhat1)
 ABI_FREE(vresid_dum)
 ABI_FREE(vtrial1_i1pert)
 ABI_FREE(vtrial1_i2pert)
 ABI_FREE(vxc1)
 ABI_FREE(vhartr1)
 ABI_FREE(vpsp1)
 ABI_FREE(xccc3d1)
 ABI_FREE(pawrhoij_read)
 ABI_FREE(ph1d)

 if (dtset%lw_flexo==1.or.dtset%lw_flexo==2.or.dtset%lw_flexo==4) then
  ABI_FREE(vhart1dqdq)
  ABI_FREE(vpsp1dqdq)
  ABI_FREE(vxc1dqdq)
 end if

!Treatment of T4 and T5 terms that have a q-gradient of a rf Hamiltonian
!they need to be converted to type-II for strain perturbation
 if (d3e_pert2(natom+3)==1.or.d3e_pert2(natom+4)==1) then
   fac=two_pi ** 2
   i3pert= natom+8
   do i1pert = 1, mpert
     do i1dir = 1, 3
       if ((maxval(rfpert(i1dir,i1pert,:,:,:,:))==1)) then
         do i2pert = natom+3, natom+4
           do i2dir = 1, 3
             istr=(i2pert-natom-3)*3+i2dir
             beta=idx(2*istr-1); delta=idx(2*istr)
             do ii=1,2

               !Transform into type-II
               do i3dir=1,3
                 gamma=i3dir
                 t4_typeII(ii,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)= &
               & t4_typeI(ii,i1dir,i1pert,beta,delta,gamma) + &
               & t4_typeI(ii,i1dir,i1pert,delta,gamma,beta) - &
               & t4_typeI(ii,i1dir,i1pert,gamma,beta,delta)
               end do ! i3dir


               !Transform i3dir into reduced coordinates
               do i3dir=1,3
                 vec1(i3dir)=t4_typeII(ii,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
                 flg1(i3dir)=blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) 
               end do
               call cart39(flg1,flg2,transpose(rprimd),natom+2,natom,transpose(gprimd),vec1,vec2)
               do i3dir=1,3
                 t4_typeII(ii,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)=vec2(i3dir)*fac
               end do

             end do ! ii
           end do ! i2dir
         end do ! i2pert
       end if ! rfpert
     end do ! i1dir
   end do ! i1pert
 end if

 if (d3e_pert1(natom+3)==1.or.d3e_pert1(natom+4)==1) then
   fac=two_pi ** 2
   i3pert= natom+8
   do i2pert = 1, mpert
     do i2dir = 1, 3
       if ((maxval(rfpert(:,:,i2dir,i2pert,:,:))==1)) then
         do i1pert = natom+3, natom+4
           do i1dir = 1, 3
             istr=(i1pert-natom-3)*3+i1dir
             beta=idx(2*istr-1); delta=idx(2*istr)
             do ii=1,2

               !Transform into type-II
               do i3dir=1,3
                 gamma=i3dir
                 t5_typeII(ii,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)= &
               & t5_typeI(ii,i2dir,i2pert,beta,delta,gamma) + &
               & t5_typeI(ii,i2dir,i2pert,delta,gamma,beta) - &
               & t5_typeI(ii,i2dir,i2pert,gamma,beta,delta)
               end do ! i3dir


               !Transform i3dir into reduced coordinates
               do i3dir=1,3
                 vec1(i3dir)=t5_typeII(ii,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
                 flg1(i3dir)=blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) 
               end do
               call cart39(flg1,flg2,transpose(rprimd),natom+2,natom,transpose(gprimd),vec1,vec2)
               do i3dir=1,3
                 t5_typeII(ii,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)=vec2(i3dir)*fac
               end do

             end do ! ii
           end do ! i2dir
         end do ! i2pert
       end if ! rfpert
     end do ! i1dir
   end do ! i1pert
 end if

!Incorporate T4 and T5 to d3etot
 d3etot=d3etot+t4_typeII+t5_typeII

!Anounce end of spatial-dispersion calculation
 write(message, '(a,a,a,a)' ) ch10,ch10,&
&   ' -- Spatial-dispersion 3rd-order derivatives completed -- ',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 

 DBG_EXIT("COLL")

end subroutine dfptlw_loop
!!***

end module m_dfptlw_loop
!!***
