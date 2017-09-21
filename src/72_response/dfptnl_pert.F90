!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfptnl_pert
!! NAME
!! dfptnl_pert
!!
!! FUNCTION
!! Compute the linear response part to the 3dte
!!
!! COPYRIGHT
!! Copyright (C) 2016-2016 ABINIT group (LB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol) = array for planewave
!!                                          coefficients of wavefunctions
!!  cg1 = first-order wavefunction relative to the perturbations i1pert
!!  cg3 = first-order wavefunction relative to the perturbations i3pert
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!          if 2, COMPLEX
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  i1dir,i2dir,i3dir=directions of the corresponding perturbations
!!  i1pert,i2pert,i3pert = type of perturbation that has to be computed
!!  kg(3,mpw*mkmem)=reduced planewave coordinates
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  mband = maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem = maximum number of k points which can fit in core memory
!!  mk1mem = maximum number of k points for first-order WF
!!           which can fit in core memory
!!  mpert =maximum number of ipert
!!  mpi_enreg=MPI-parallelisation information
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  natom = number of atoms in unit cell
!!  nfft  = (effective) number of FFT grid points (for this processor)
!!  nkpt  = number of k points
!!  nspden = number of spin-density components
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  npwarr(nkpt) = array holding npw for each k point
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  rprimd(3,3) = dimensional primitive translations (bohr)
!!  vtrial1(cplex*nfft,nspden)=firs-order local potential
!!  xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!  d3etot(2,3,mpert,3,mpert,3,mpert) = matrix of the 3DTEs
!!
!! PARENTS
!!      dfptnl_loop
!!
!! CHILDREN
!!      destroy_hamiltonian,dotprod_g,fftpac,fourwf,init_hamiltonian
!!      load_k_hamiltonian,mkffnl,mkkpg,nonlop,status,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfptnl_pert(atindx,atindx1,cg,cg1,cg2,cg3,cplex,dtfil,dtset,d3etot,eigen0,gs_hamkq,k3xc,indsy1,i1dir,i2dir,i3dir,&
& i1pert,i2pert,i3pert,kg,mband,mgfft,mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,natom,nattyp,nfftf,nfftotf,ngfftf,nkpt,nk3xc,&
& nspden,nspinor,nsppol,nsym1,npwarr,occ,pawang,pawang1,pawfgrtab,pawrad,pawtab,&
& pawrhoij0,pawrhoij1_i1pert,pawrhoij1_i2pert,pawrhoij1_i3pert,&
& paw_an0,paw_an1_i2pert,paw_ij0,paw_ij1_i2pert,pawfgr,ph1d,psps,rho1r1,rho2r1,rho3r1,rprimd,symaf1,symrc1,&
& ucvol,vtrial,vhartr1_i2pert,vtrial1_i2pert,vxc1_i2pert,ddk_f,xccc3d1,xccc3d2,xccc3d3,xred,&
& d3etot_1,d3etot_2,d3etot_3,d3etot_4,d3etot_5,d3etot_6,d3etot_7,d3etot_8,d3etot_9)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_wffile
 use m_wfk
 use m_xmpi
 use m_hamiltonian
 use m_errors
 use m_rf2

 use m_cgtools,    only : dotprod_g,sqnorm_g
 use m_pawang,     only : pawang_type
 use m_pawfgrtab,  only : pawfgrtab_type
 use m_pawrad,     only : pawrad_type
 use m_pawtab,     only : pawtab_type
 use m_pawcprj,    only : pawcprj_type, pawcprj_alloc, pawcprj_free
 use m_paw_ij,     only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify,paw_ij_reset_flags
 use m_pawdij,     only : pawdijfr
 use m_pawfgr,     only : pawfgr_type
 use m_pawrhoij,   only : pawrhoij_type, pawrhoij_alloc , pawrhoij_nullify, pawrhoij_free,&
&                         pawrhoij_init_unpacked, pawrhoij_mpisum_unpacked
 use m_paw_an,     only : paw_an_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptnl_pert'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_53_spacepar
 use interfaces_56_recipspace
 use interfaces_65_paw
 use interfaces_66_nonlocal
 use interfaces_66_wfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,mband,mgfft
 integer,intent(in) :: mk1mem,mkmem,mpert,mpsang,mpw,natom,nfftf,nfftotf,nkpt,nspden,nsym1
 integer,intent(in) :: nk3xc,nspinor,nsppol
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(pawang_type),intent(inout) :: pawang,pawang1
 type(pawfgr_type),intent(in) :: pawfgr
 type(wfk_t),intent(inout) :: ddk_f(5)

!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),kg(3,mpw*mkmem),nattyp(psps%ntypat),ngfftf(18),npwarr(nkpt)
 integer,intent(in) :: indsy1(4,nsym1,dtset%natom),symaf1(nsym1),symrc1(3,3,nsym1)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: cg2(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: cg3(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: eigen0(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: k3xc(nfftf,nk3xc)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rho1r1(cplex*nfftf,dtset%nspden),rho2r1(cplex*nfftf,dtset%nspden)
 real(dp),intent(in) :: rho3r1(cplex*nfftf,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: vtrial(cplex*nfftf,nspden)
 real(dp),intent(in) :: xccc3d1(cplex*nfftf),xccc3d2(cplex*nfftf),xccc3d3(cplex*nfftf),xred(3,natom)
 real(dp),intent(in) :: vxc1_i2pert(cplex*nfftf,nspden),vhartr1_i2pert(cplex*nfftf)
 real(dp),intent(inout) :: vtrial1_i2pert(cplex*nfftf,nspden),d3etot(2,3,mpert,3,mpert,3,mpert)
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
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij0(natom*psps%usepaw)
 type(pawrhoij_type),intent(in),target :: pawrhoij1_i1pert(natom*psps%usepaw)
 type(pawrhoij_type),intent(in)        :: pawrhoij1_i2pert(natom*psps%usepaw)
 type(pawrhoij_type),intent(in),target :: pawrhoij1_i3pert(natom*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_an_type),intent(in) :: paw_an0(natom*psps%usepaw)
 type(paw_an_type),intent(inout) :: paw_an1_i2pert(natom*psps%usepaw)
 type(paw_ij_type),intent(in) :: paw_ij0(natom*psps%usepaw)
 type(paw_ij_type),intent(inout) :: paw_ij1_i2pert(natom*psps%usepaw)

!Local variables-------------------------------
!scalars
 logical :: has_cprj_jband,compute_conjugate,compute_rho21,usetimerev
 integer,parameter :: level=52
 integer :: bandtot,choice,counter,cplex_cprj,cplex_loc,cplex_rhoij,cpopt,dimffnl1,dimffnlk,iband,icg0,ider,ierr,iexit
 integer :: ic,jc,kc,idir0,idir_getgh2c,idir_phon,idir_elfd,ipert_phon,idir_54,ipert_elfd,ispden
 integer :: ia,iatm,ibg,ii,igs,ikg,ikg1,ikpt,ifft,ilm,ipw,isppol,ispinor,istwf_k,jband,jj
 integer :: me,n1,n2,n3,n4,n5,n6,nband_k,ncpgr,nkpg,nkpg1,nnlout,nsp,nspden_rhoij,npert_phon,npw_k,npw1_k,nzlmopt
 integer :: offset_cgi,offset_cgj,offset_eigen,offset_eig0,option,paw_opt,debug_mode,esigns
 integer :: signs,size_wf,size_cprj,spaceComm,tim_fourwf,tim_nonlop,tim_getgh2c,usepaw,useylmgr1
 integer :: sij_opt,usevnl,opt_gvnl2,optnl,optlocal
 real(dp) :: arg,dot1i,dot1r,dot2i,dot2r,doti,dotr,e3tot,lagi,lagi_paw,lagr,lagr_paw
 real(dp) :: rho2r_re,rho2r_im,rho3r_re,rho3r_im
 real(dp) :: sumi,sum_psi1H1psi1,sum_psi1H1psi1_i
 real(dp) :: sum_lambda1psi1psi1,sum_lambda1psi1psi1_i
 real(dp) :: sum_psi0H2psi1a,sum_psi0H2psi1a_i,sum_psi0H2psi1b,sum_psi0H2psi1b_i
 real(dp) :: sum_lambda1psi0S1psi1,sum_lambda1psi0S1psi1_i,tol_test,weight
 character(len=1000) :: msg
!arrays
 integer,allocatable :: kg_k(:,:),kg1_k(:,:)
! real(dp) :: buffer(2)
 real(dp) :: buffer(10),eHxc21_paw(2),eHxc21_nhat(2),exc3(2),exc3_paw(2),enlout(3),kpt(3),eig0_k(mband),dum_svectout(1,1),dum(1)
 real(dp) :: enlout1(2),enlout2(2),enlout_11(2),enlout_12(2),enlout_21(2),enlout_22(2),enlout_31(2),enlout_32(2)
 real(dp) :: rmet(3,3),dum_grad_berry(1,1),wtk_k
 real(dp) :: ylmgr_dum(1,1,1)
 real(dp),allocatable :: cgi(:,:),cgj(:,:),cg_jband(:,:,:),cwavef1(:,:),cwavef2(:,:),cwavef3(:,:),dkinpw(:)
 real(dp),allocatable :: eig1_k_i2pert(:),eig1_k_stored(:)
 real(dp),allocatable :: chi_ij(:,:,:),cwave_right(:,:),cwave_left(:,:),dudk(:,:),dudkde(:,:),dummy_array(:),dummy_array2(:,:)
 real(dp),allocatable :: ffnl1(:,:,:,:),ffnl1_test(:,:,:,:),gh0(:,:),gh1(:,:),gvnl(:,:)
 real(dp),allocatable :: h_cwave(:,:),iddk(:,:),kinpw1(:),kpg_k(:,:),kpg1_k(:,:),nhat21(:,:),nhatfr21(:,:),occ_k(:)
 real(dp),allocatable :: phkxred(:,:),ph3d(:,:,:),ph3d1(:,:,:),rho1r1_tot(:,:),s_cwave(:,:)
 real(dp),allocatable :: vlocal(:,:,:,:),vlocal1_i2pert(:,:,:,:),v_i2pert(:,:),wfraug(:,:,:,:)
 real(dp),allocatable :: ylm(:,:),ylm1(:,:),ylmgr(:,:,:),ylmgr1(:,:,:)
 real(dp),allocatable :: ylm_k(:,:),ylm1_k(:,:),ylmgr1_k(:,:,:)
 real(dp),allocatable :: xc_tmp(:,:)
! type(pawcprj_type),allocatable :: cwaveprj(:,:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:),cwaveprj1(:,:)
 type(pawcprj_type),target :: cprj_empty(0,0)
 type(pawcprj_type),allocatable,target :: cprj_jband(:,:)
 type(pawrhoij_type),allocatable,target  :: pawrhoij21(:)
 type(pawrhoij_type),pointer :: pawrhoij21_unsym(:),pawrhoij11(:)
 type(paw_ij_type),allocatable :: paw_ij_tmp(:)
 type(rf_hamiltonian_type) :: rf_hamkq_i2pert,rf_ham_dum

!***********************************************************************

 DBG_ENTER("COLL")

 me = mpi_enreg%me
 spaceComm=mpi_enreg%comm_cell

 call status(0,dtfil%filstat,iexit,level,'enter         ')

 npert_phon = 0
 if(i1pert<=natom) npert_phon = npert_phon + 1
 if(i2pert<=natom) npert_phon = npert_phon + 1
 if(i3pert<=natom) npert_phon = npert_phon + 1
 if (npert_phon>1) then
   MSG_ERROR("dfptnl_pert is available with at most one phonon perturbation. Change your input!")
 end if

 usepaw = psps%usepaw
 size_cprj = nspinor

 call init_rf_hamiltonian(cplex,gs_hamkq,i2pert,rf_hamkq_i2pert,paw_ij1=paw_ij1_i2pert,has_e1kbsc=.true.)

!Acivate computation of rho^(2:1) and related energy derivatives if needed
 compute_rho21 = .false.
 if (usepaw==1.and.npert_phon==1.and.(i1pert<=natom.or.i3pert<=natom)) then ! so i2pert==natom+2
   compute_rho21 = .true.
   if (i1pert<=natom) then
     ipert_phon  = i1pert
     idir_phon   = i1dir
     ipert_elfd = i3pert
     idir_elfd  = i3dir
     pawrhoij11 => pawrhoij1_i3pert
   else if (i3pert<=natom) then
     ipert_phon  = i3pert
     idir_phon   = i3dir
     ipert_elfd = i1pert
     idir_elfd  = i1dir
     pawrhoij11 => pawrhoij1_i1pert
   end if
   cplex_rhoij=max(cplex,dtset%pawcpxocc);nspden_rhoij=dtset%nspden
   ABI_DATATYPE_ALLOCATE(pawrhoij21,(natom))
   call pawrhoij_nullify(pawrhoij21)
   call pawrhoij_alloc(pawrhoij21,cplex_rhoij,nspden_rhoij,nspinor,dtset%nsppol,&
&     dtset%typat,pawtab=pawtab,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   ABI_DATATYPE_ALLOCATE(cwaveprj0,(natom,size_cprj))
   ABI_DATATYPE_ALLOCATE(cwaveprj1,(natom,size_cprj))
   call pawcprj_alloc(cwaveprj0,1,gs_hamkq%dimcprj)
   call pawcprj_alloc(cwaveprj1,1,gs_hamkq%dimcprj)
!   if (paral_atom) then
!     ABI_DATATYPE_ALLOCATE(pawrhoij1_unsym,(natom))
!     cplex_rhoij=max(cplex,dtset%pawcpxocc);nspden_rhoij=dtset%nspden
!     call pawrhoij_alloc(pawrhoij1_unsym,cplex_rhoij,nspden_rhoij,dtset%nspinor,&
!&     dtset%nsppol,dtset%typat,pawtab=pawtab,use_rhoijp=0,use_rhoij_=1)
!   else
     pawrhoij21_unsym => pawrhoij21
     call pawrhoij_init_unpacked(pawrhoij21_unsym)
!   end if
!  Compute phkxred :
   ABI_ALLOCATE(phkxred,(2,natom))
   do ia=1,natom
     iatm=min(atindx(ia),natom)
     arg=two_pi*(kpt(1)*xred(1,ia)+kpt(2)*xred(2,ia)+kpt(3)*xred(3,ia))
     phkxred(1,iatm)=cos(arg);phkxred(2,iatm)=sin(arg)
   end do
   ABI_DATATYPE_ALLOCATE(paw_ij_tmp,(natom))
   call paw_ij_nullify(paw_ij_tmp)
   cplex_loc=1;nsp=1 ! Force nsppol/nspden to 1 because Dij^(1) due to electric field is spin-independent
   call paw_ij_init(paw_ij_tmp,cplex_loc,dtset%nspinor,nsp,nsp,dtset%pawspnorb,natom,psps%ntypat,&
&   dtset%typat,pawtab,has_dijfr=1,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call paw_ij_reset_flags(paw_ij_tmp,all=.True.)
   call pawdijfr(cplex_loc,gs_hamkq%gprimd,idir_elfd,ipert_elfd,natom,natom,nfftf,ngfftf,nsp,psps%ntypat,&
&     1,paw_ij_tmp,pawang,pawfgrtab,pawrad,pawtab,&
&     (/zero,zero,zero/),rprimd,ucvol,dummy_array2,dummy_array2,dummy_array2,xred,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   ABI_ALLOCATE(chi_ij,(gs_hamkq%dimekb1,gs_hamkq%dimekb2,dtset%nspinor**2))
   if(dtset%nsppol==1) then
     call pawdij2e1kb(paw_ij_tmp,1,mpi_enreg%comm_atom,mpi_enreg%my_atmtab,e1kbfr=chi_ij)
   else
     MSG_ERROR("nonlinear not implemented yet for nsppol>1...")
   end if
   call paw_ij_free(paw_ij_tmp)
   ABI_DATATYPE_DEALLOCATE(paw_ij_tmp)
 else
   ABI_ALLOCATE(phkxred,(0,0))
   ABI_DATATYPE_ALLOCATE(pawrhoij21,(0))
   pawrhoij21_unsym => pawrhoij21
   ABI_DATATYPE_ALLOCATE(cwaveprj0,(0,0))
   ABI_DATATYPE_ALLOCATE(cwaveprj1,(0,0))
 end if

 bandtot = 0
 icg0 = 0
 option = 2
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)

 ABI_ALLOCATE(vlocal,(n4,n5,n6,gs_hamkq%nvloc))
 ABI_ALLOCATE(vlocal1_i2pert,(cplex*n4,n5,n6,gs_hamkq%nvloc))

 ABI_ALLOCATE(wfraug,(2,n4,n5,n6))

 rmet = MATMUL(TRANSPOSE(rprimd),rprimd)

 sumi = zero

!Set up the Ylm for each k point
 ABI_ALLOCATE(ylm,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))
 ABI_ALLOCATE(ylmgr,(dtset%mpw*dtset%mkmem,9,psps%mpsang*psps%mpsang*psps%useylm))
 if (psps%useylm==1) then
   call status(0,dtfil%filstat,iexit,level,'call initylmg ')
   option=2
   call initylmg(gs_hamkq%gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%nband,&
   dtset%nkpt,npwarr,dtset%nsppol,option,rprimd,ylm,ylmgr)
 end if

!Set up the spherical harmonics (Ylm) at k+q
 useylmgr1=0; option=0
 if (psps%useylm==1.and. &
& (i2pert==natom+1.or.i2pert==natom+3.or.i2pert==natom+4.or.(usepaw==1.and.i2pert==natom+2))) then
   useylmgr1=1; option=1
 end if
 ABI_ALLOCATE(ylm1,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))
 ABI_ALLOCATE(ylmgr1,(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))
!To change the following when q/=0
 if (psps%useylm==1) then
   call initylmg(gs_hamkq%gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%nband,&
   dtset%nkpt,npwarr,dtset%nsppol,option,rprimd,ylm1,ylmgr1)
 end if

 debug_mode = 0
 if (dtset%nonlinear_info>3.or.dtset%nonlinear_info==2) debug_mode = 1

!Real parts
 sum_psi1H1psi1 =  zero
 sum_lambda1psi1psi1 = zero
 sum_lambda1psi0S1psi1 = zero
 sum_psi0H2psi1a = zero
 sum_psi0H2psi1b = zero
!Imaginary parts
 sum_psi1H1psi1_i =  zero
 sum_lambda1psi1psi1_i = zero
 sum_lambda1psi0S1psi1_i = zero
 sum_psi0H2psi1a_i = zero
 sum_psi0H2psi1b_i = zero

 compute_conjugate = .false.
!We have to compute < u^(ip1) | H^(ip2) | u^(ip3) >
!For some cases, we want to apply H^(ip2) on < u^(ip1) |, not on | u^(ip3) > (see below)
 if (i3pert<=natom) then ! As npert_phon<=1, this implies that i1pert=natom+2 and i2pert=natom+2
   compute_conjugate = .true.
 end if

!Loop over spins
 do isppol = 1, nsppol

!  Set up local potential vlocal1 with proper dimensioning, from vtrial1
!  Same thing for vlocal from vtrial Also take into account the spin.
   call rf_transgrid_and_pack(isppol,nspden,usepaw,cplex,nfftf,dtset%nfft,dtset%ngfft,&
&   gs_hamkq%nvloc,pawfgr,mpi_enreg,vtrial,vtrial1_i2pert,vlocal,vlocal1_i2pert)

!  Continue to initialize the Hamiltonian
   call load_spin_hamiltonian(gs_hamkq,isppol,vlocal=vlocal,with_nonlocal=.true.)
   call load_spin_rf_hamiltonian(rf_hamkq_i2pert,gs_hamkq,isppol,vlocal1=vlocal1_i2pert,with_nonlocal=.true.)

!  Loop over k-points

   ikg = 0
   ikg1 = 0

   do ikpt = 1, nkpt

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,-1,mpi_enreg%me)) then
       cycle ! Skip the rest of the k-point loop
     end if

     counter = 100*ikpt

     nband_k = dtset%nband(ikpt+(isppol-1)*nkpt)
     npw_k = npwarr(ikpt)
     npw1_k = npw_k ! To change for q/=0
     istwf_k = dtset%istwfk(ikpt)
     ABI_ALLOCATE(occ_k,(nband_k))
     occ_k(:) = occ(1+bandtot:nband_k+bandtot)
     wtk_k    = dtset%wtk(ikpt)

     size_wf = nspinor*npw_k

     kpt(:) = dtset%kptns(:,ikpt)

     ABI_ALLOCATE(cwavef1,(2,npw_k*nspinor))
     ABI_ALLOCATE(cwavef3,(2,npw_k*nspinor))
     if (compute_rho21) then
       ABI_ALLOCATE(cwavef2,(2,npw_k*nspinor))
     end if
     ABI_ALLOCATE(gh0,(2,npw_k*nspinor))
     ABI_ALLOCATE(gvnl,(2,npw_k*nspinor))
     ABI_ALLOCATE(gh1,(2,npw_k*nspinor))

     ABI_ALLOCATE(kg_k,(3,npw_k))
     ABI_ALLOCATE(kg1_k,(3,npw1_k))
     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     ABI_ALLOCATE(ylm1_k,(npw1_k,mpsang*mpsang*psps%useylm))
     ABI_ALLOCATE(ylmgr1_k,(npw1_k,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))

!    Get (k+G) wave vectors and associated spherical harmonics
     kg_k(:,1:npw_k) = kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,mpsang*mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if

!    Get (k+q+G) wave vectors and associated spherical harmonics
     kg1_k(:,1:npw1_k)=kg(:,1+ikg1:npw1_k+ikg1) ! To change for q/=0
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

!!    Set up the ground-state Hamiltonian, and some parts of the 1st-order Hamiltonian
!     call getgh1c_setup(gs_hamkq,rf_hamkq_i2pert,dtset,psps,&                       ! In
!     kpt,kpt,i2dir,i2pert,natom,rmet,gs_hamkq%gprimd,gs_hamkq%gmet,istwf_k,&        ! In
!     npw_k,npw1_k,useylmgr1,kg_k,ylm_k,kg1_k,ylm1_k,ylmgr1_k,&                      ! In
!     dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1,&                 ! Out
!     dummy_array,dummy_array,rf_ham_dum)                                            ! Out

!    Compute (k+G) vectors
     nkpg=0;if(i2pert>=1.and.i2pert<=natom) nkpg=3*dtset%nloalg(3)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpt,nkpg,npw_k)
     end if

!    Compute (k+q+G) vectors
     nkpg1=0;if(i2pert>=1.and.i2pert<=natom) nkpg1=3*dtset%nloalg(3)
     ABI_ALLOCATE(kpg1_k,(npw1_k,nkpg1))
     if (nkpg1>0) then
       call mkkpg(kg1_k,kpg1_k,kpt,nkpg1,npw1_k)
     end if

!    ===== Preparation of the non-local contributions

!    Compute nonlocal form factors ffnl1 at (k+q+G)
     !-- Atomic displacement perturbation
     if (i2pert<=natom) then
       ider=0;idir0=0
!     !-- k-point perturbation (1st-derivative)
!     else if (i2pert==natom+1) then
!       ider=1;idir0=idir
!     !-- k-point perturbation (2nd-derivative)
!     else if (i2pert==natom+10.or.i2pert==natom+11) then
!       ider=2;idir0=4
     !-- Electric field perturbation
     else if (i2pert==natom+2) then
       if (psps%usepaw==1) then
         ider=1;idir0=i2dir
       else
         ider=0;idir0=0
       end if
     !-- Strain perturbation
!     else if (i2pert==natom+3.or.i2pert==natom+4) then
!       if (ipert==natom+3) istr=idir
!       if (ipert==natom+4) istr=idir+3
!       ider=1;idir0=-istr
     end if
     if (compute_rho21) then ! compute_rho21 implies i2pert==natom+2
       ider=1; idir0=4
     end if

!    Compute nonlocal form factors ffnl1 at (k+q+G), for all atoms
     dimffnl1=1+ider
     if (ider==1.and.(idir0==0.or.idir0==4)) dimffnl1=2+2*psps%useylm
!     if (ider==2.and.idir0==4) dimffnl1=3+7*psps%useylm
     ABI_ALLOCATE(ffnl1,(npw1_k,dimffnl1,psps%lmnmax,psps%ntypat))
     call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1,psps%ffspl,gs_hamkq%gmet,gs_hamkq%gprimd,ider,idir0,&
&      psps%indlmn,kg1_k,kpg1_k,kpt,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,&
&      npw1_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1_k,ylmgr1_k)

!    Compute nonlocal form factors ffnl1 at (k+q+G), for all atoms
     if (compute_rho21.and.debug_mode/=0) then
       ABI_ALLOCATE(ffnl1_test,(npw1_k,dimffnl1,psps%lmnmax,psps%ntypat))
       idir0 = 0 ! for nonlop with signs = 1
       call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1_test,psps%ffspl,gs_hamkq%gmet,gs_hamkq%gprimd,ider,idir0,&
  &      psps%indlmn,kg1_k,kpg1_k,kpt,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,&
  &      npw1_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1_k,ylmgr1_k)
     end if

!    ===== Preparation of the kinetic contributions

!    Note that not all these arrays should be allocated in the general case when wtk_k vanishes

!    Compute (1/2) (2 Pi)**2 (k+q+G)**2:
     ABI_ALLOCATE(kinpw1,(npw1_k))
     kinpw1(:)=zero
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gs_hamkq%gmet,kg1_k,kinpw1,kpt,npw1_k,0,0)

     ABI_ALLOCATE(dkinpw,(npw_k)) ! 1st derivative (1st direction)
     dkinpw(:)=zero
!     if(ipert==natom+10 .and. idir>3) then
!       ABI_ALLOCATE(dkinpw2,(npw_k)) ! 1st derivative (2nd directions)
!       dkinpw2(:)=zero
!     end if
!     if(ipert==natom+10) then
!       ABI_ALLOCATE(ddkinpw,(npw_k)) ! 2nd derivative
!       ddkinpw(:)=zero
!     end if

!!-- k-point perturbation (1st-derivative)
! if (ipert==natom+1) then
!!  Compute the derivative of the kinetic operator vs k
!   call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg_k,dkinpw,kpt,npw_k,idir,0) ! 1st derivative
! end if

!  !-- Strain perturbation
! if (ipert==natom+3.or.ipert==natom+4) then
!   if (ipert==natom+3) istr=idir
!   if (ipert==natom+4) istr=idir+3
!!  Compute the derivative of the kinetic operator vs strain
!   call kpgstr(dkinpw,dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,gprimd,istr,kg_k,kpt,npw_k)
! end if

!===== Load the k/k+q dependent parts of the Hamiltonian

!  Load k-dependent part in the Hamiltonian datastructure
   ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamkq%matblk))
   call load_k_hamiltonian(gs_hamkq,kpt_k=kpt,npw_k=npw_k,istwf_k=istwf_k,kg_k=kg_k,kpg_k=kpg_k,&
&   ph3d_k=ph3d,compute_ph3d=.true.,compute_gbound=.true.)
!   if (size(ffnlk)>0) then
!     call load_k_hamiltonian(gs_hamkq,ffnl_k=ffnlk)
!   else
   call load_k_hamiltonian(gs_hamkq,ffnl_k=ffnl1,kpt_k=kpt,npw_k=npw1_k,istwf_k=istwf_k,&
&     kinpw_k=kinpw1,kg_k=kg1_k,kpg_k=kpg1_k,compute_gbound=.true.)
!   end if

!    Load k+q-dependent part in the Hamiltonian datastructure
!      Note: istwf_k is imposed to 1 for RF calculations (should use istwf_kq instead)
!     call load_kprime_hamiltonian(gs_hamkq,kpt_kp=kpt,npw_kp=npw1_k,istwf_kp=istwf_k,&
!&     kinpw_kp=kinpw1,kg_kp=kg1_k,kpg_kp=kpg1_k,ffnl_kp=ffnl1,&
!&     compute_gbound=.true.)
!   if (qne0) then
!     ABI_ALLOCATE(ph3d1,(2,npw1_k,gs_hamkq%matblk))
!     call load_kprime_hamiltonian(gs_hamkq,ph3d_kp=ph3d1,compute_ph3d=.true.)
!   end if

!    Load k-dependent part in the 1st-order Hamiltonian datastructure
     call load_k_rf_hamiltonian(rf_hamkq_i2pert,npw_k=npw_k,dkinpw_k=dkinpw)
!     if (ipert==natom+10) then
!       call load_k_rf_hamiltonian(rf_hamkq,ddkinpw_k=ddkinpw)
!!       if (idir>3) then
!!         call load_k_rf_hamiltonian(rf_hamk_dir2,dkinpw_k=dkinpw2,ddkinpw_k=ddkinpw)
!!       end if
!     end if

!!    Set up the ground-state Hamiltonian, and some parts of the 1st-order Hamiltonian
!     call getgh1c_setup(gs_hamkq,rf_hamkq_i2pert,dtset,psps,&                       ! In
!     kpt,kpt,i2dir,i2pert,natom,rmet,gs_hamkq%gprimd,gs_hamkq%gmet,istwf_k,&        ! In
!     npw_k,npw1_k,useylmgr1,kg_k,ylm_k,kg1_k,ylm1_k,ylmgr1_k,&                      ! In
!     dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1,&                 ! Out
!     dummy_array,dummy_array,rf_ham_dum)                                            ! Out

     ABI_STAT_ALLOCATE(dudk,  (2,nband_k*size_wf), ierr)
     ABI_STAT_ALLOCATE(dudkde,(2,nband_k*size_wf), ierr)
     ABI_STAT_ALLOCATE(eig1_k_i2pert,(2*nband_k), ierr)
     ABI_ALLOCATE(eig1_k_stored,(2*nband_k**2))
     ABI_ALLOCATE(cgi,(2,size_wf))
     ABI_ALLOCATE(cwave_right,(2,size_wf))
     ABI_ALLOCATE(cwave_left,(2,size_wf))

! **************************************************************************************************
!      Read dudk and dudkde
! **************************************************************************************************

     do iband = 1,nband_k

!      Read dude file
       call wfk_read_bks(ddk_f(1), iband, ikpt, isppol, xmpio_single, cg_bks=cwave_right,eig1_bks=eig1_k_i2pert)
!      Copy eig1_k_i2pert in "eig1_k_stored"
       eig1_k_stored(1+(iband-1)*2*nband_k:2*nband_k+(iband-1)*2*nband_k)=eig1_k_i2pert(:)

       if (i2pert==natom+2) then
!        Read dudk file
         call wfk_read_bks(ddk_f(2), iband, ikpt, isppol, xmpio_single, cg_bks=cwave_right,eig1_bks=eig1_k_i2pert)
         offset_cgi = (iband-1)*size_wf+icg0
         cgi(:,:) = cg(:,1+offset_cgi:size_wf+offset_cgi)
!         if (usepaw==0) then
!           call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cgi,cgi,mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
!           if (abs(dotr-1)>tol10.or.abs(doti)>tol10) then
!             print '(2(a,es19.10E3))','       |cgi|^2 = ',dotr,',',doti
!           end if
!           call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cgi,cwave_right,mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
!           if (abs(dotr)>tol10.or.abs(doti)>tol10) then
!             print '(2(a,es19.10E3))',' < cgi | ddk > = ',dotr,',',doti
!           end if
!         end if
!        Copy cwave_right in "dudk"
         dudk(:,1+(iband-1)*size_wf:iband*size_wf)=cwave_right(:,:)

!        Read dudkde file
         call wfk_read_bks(ddk_f(3), iband, ikpt, isppol, xmpio_single, cg_bks=cwave_right,eig1_bks=eig1_k_i2pert)
         offset_cgi = (iband-1)*size_wf+icg0
         cgi(:,:) = cg(:,1+offset_cgi:size_wf+offset_cgi)
!         if (usepaw==0) then
!           call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cgi,cgi,mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
!           if (abs(dotr-1)>tol10.or.abs(doti)>tol10) then
!             print '(2(a,es19.10E3))','       |cgi|^2 = ',dotr,',',doti
!           end if
!         end if
!        Copy cwave_right in "dudkde"
         dudkde(:,1+(iband-1)*size_wf:iband*size_wf)=cwave_right(:,:)
       end if

     end do

     ABI_ALLOCATE(cgj,(2,size_wf))
     ABI_ALLOCATE(iddk,(2,size_wf))

     offset_eig0 = mband*(ikpt-1+nkpt*(isppol-1))
     eig0_k(:) = eigen0(1+offset_eig0:mband+offset_eig0)

     ABI_STAT_ALLOCATE(h_cwave,(2,size_wf), ierr)
     ABI_STAT_ALLOCATE(s_cwave,(2,size_wf), ierr)

!    Allocate work spaces when debug_mode is activated
     has_cprj_jband=.false.
     if (debug_mode/=0) then ! Only for test purposes
       ABI_ALLOCATE(cg_jband,(2,size_wf*nband_k,2))
       cg_jband(:,:,1) = cg(:,1+icg0:size_wf*nband_k+icg0)
       if (i2pert==natom+2) then ! Note the multiplication by "i"
         cg_jband(1,:,2) = -dudk(2,1:size_wf*nband_k)
         cg_jband(2,:,2) =  dudk(1,1:size_wf*nband_k)
       end if
       if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
         ABI_DATATYPE_ALLOCATE(cprj_jband,(natom,size_cprj*nband_k))
         has_cprj_jband=.true.
       else
         ABI_DATATYPE_ALLOCATE(cprj_jband,(natom,0))
       end if
     else
       ABI_ALLOCATE(cg_jband,(2,0,2))
       ABI_DATATYPE_ALLOCATE(cprj_jband,(natom,0))
     end if

!    Loop over bands
     do jband = 1,nband_k
       !  Skip bands not treated by current proc
       if((mpi_enreg%proc_distrb(ikpt,jband,isppol)/=me)) cycle
       if (occ_k(jband)>tol10) then

  !      tol_test = tol8
         offset_cgj = (jband-1)*size_wf+icg0
         cgj(:,:) = cg(:,1+offset_cgj:size_wf+offset_cgj)

! **************************************************************************************************
!        Compute < u^(1) | ( H^(1) - eps^(0) S^(1) ) | u^(1) >
! **************************************************************************************************

         eig1_k_i2pert(:) = eig1_k_stored(1+(jband-1)*2*nband_k:jband*2*nband_k)
         cwavef1(:,:) = cg1(:,1+offset_cgj:size_wf+offset_cgj)
         cwavef3(:,:) = cg3(:,1+offset_cgj:size_wf+offset_cgj)
         if (i2pert==natom+2) then ! Note the multiplication by i
           iddk(1,:) = -dudkde(2,1+(jband-1)*size_wf:jband*size_wf)
           iddk(2,:) =  dudkde(1,1+(jband-1)*size_wf:jband*size_wf)
         else
           iddk(:,:) = zero
         end if
         cwave_right(:,:) = cwavef3(:,:)
         cwave_left(:,:)  = cwavef1(:,:)
         if (compute_conjugate) then
           cwave_right(:,:) = cwavef1(:,:)
           cwave_left(:,:)  = cwavef3(:,:)
         end if

  !      Compute : < u^(ip1) | ( H^(ip2) - eps^(0) S^(ip2) ) | u^(ip3) >
  !           or : < u^(ip3) | ( H^(ip2) - eps^(0) S^(ip2) ) | u^(ip1) >
         call rf2_apply_hamiltonian(cg_jband,cprj_jband,cwave_right,cprj_empty,h_cwave,s_cwave,eig0_k,eig1_k_i2pert,&
  &                                jband,gs_hamkq,iddk,i2dir,i2pert,ikpt,isppol,mkmem,mpi_enreg,nband_k,nsppol,&
!                                   debug_mode,dtset%prtvol,rf_hamkq_i2pert,size_cprj,size_wf,conj=compute_conjugate)
                                   debug_mode,dtset%prtvol,rf_hamkq_i2pert,size_cprj,size_wf)
         call dotprod_g(dotr,doti,gs_hamkq%istwf_k,size_wf,2,cwave_left,h_cwave,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

         if (usepaw==1.and.i2pert/=natom+2) then ! S^(1) is zero for ipert=natom+2
           call dotprod_g(dot2r,dot2i,gs_hamkq%istwf_k,size_wf,2,cwave_left,s_cwave,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
           dotr = dotr - eig0_k(jband)*dot2r
           doti = doti - eig0_k(jband)*dot2i
         end if
         if (compute_conjugate) doti = -doti

! **************************************************************************************************
!        Compute sum_i Lambda_ij^(1) < u_j^(1) | u_i^(1)>
! **************************************************************************************************

!         eig1_k_i2pert(:) = eig1_k_stored(1+(jband-1)*2*nband_k:jband*2*nband_k)
         lagr = zero ; lagi = zero
         lagr_paw = zero ; lagi_paw = zero

         do iband = 1, nband_k
           if(occ_k(iband)>tol10) then

             offset_cgi = (iband-1)*size_wf+icg0
             cwavef3(:,:) = cg3(:,1+offset_cgi:size_wf+offset_cgi)

             if(debug_mode/=0) then
               eig1_k_i2pert(:) = eig1_k_stored(1+(jband-1)*2*nband_k:jband*2*nband_k)
             end if

!            Get Lambda_ij^(1) = < u_i^(0) | H^(1) - eps^(0) S^(1) | u_j^(0) > (see dfpt_cgwf.F90)
             dot1r = eig1_k_i2pert(2*iband-1)
             dot1i = eig1_k_i2pert(2*iband  )

!            Compute < u_j^(1) | S^(0) | u_i^(1) >
             if (usepaw==1) then
               ibg = 0
               call getgsc(cwavef3,cprj_empty,gs_hamkq,s_cwave,ibg,0,0,ikpt,isppol,&
&                    size_wf,size_cprj,size_wf,mpi_enreg,natom,-1,npw_k,nspinor,select_k=KPRIME_H_KPRIME)
             else
               s_cwave(:,:) = cwavef3(:,:)
             end if
             call dotprod_g(dot2r,dot2i,gs_hamkq%istwf_k,size_wf,2,cwavef1,s_cwave,mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
             lagr = lagr + dot1r*dot2r - dot1i*dot2i
             lagi = lagi + dot1r*dot2i + dot1i*dot2r

!            Compute < u_j^(0) | S^(1) | u_i^(1) >
             if (usepaw==1.and.i1pert<=natom) then ! S^(1) is zero for ipert=natom+2

               cpopt=-1+5*gs_hamkq%usecprj ; choice=2 ; signs=2 ; paw_opt=3
               call nonlop(choice,cpopt,cprj_empty,dummy_array,gs_hamkq,i1dir,(/zero/),mpi_enreg,1,nnlout,&
          &     paw_opt,signs,s_cwave,tim_nonlop,cwavef3,dummy_array2,iatom_only=i1pert)
               call dotprod_g(dot2r,dot2i,gs_hamkq%istwf_k,size_wf,2,cgj,s_cwave,mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
               lagr_paw = lagr_paw + dot1r*dot2r - dot1i*dot2i
               lagi_paw = lagi_paw + dot1r*dot2i + dot1i*dot2r

             end if

!            Compute < u_j^(1) | S^(1) | u_i^(0) >
             if (usepaw==1.and.i3pert<=natom) then ! S^(1) is zero for ipert=natom+2

               cgi(:,:) = cg(:,1+offset_cgi:size_wf+offset_cgi)
               cpopt=-1+5*gs_hamkq%usecprj ; choice=2 ; signs=2 ; paw_opt=3
               call nonlop(choice,cpopt,cprj_empty,dummy_array,gs_hamkq,i3dir,(/zero/),mpi_enreg,1,nnlout,&
          &     paw_opt,signs,s_cwave,tim_nonlop,cgi,dummy_array2,iatom_only=i3pert)
               call dotprod_g(dot2r,dot2i,gs_hamkq%istwf_k,size_wf,2,cwavef1,s_cwave,mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
               lagr_paw = lagr_paw + dot1r*dot2r - dot1i*dot2i
               lagi_paw = lagi_paw + dot1r*dot2i + dot1i*dot2r

             end if

           end if
         end do    ! iband

! **************************************************************************************************
!        Sum all band_by_band contributions
! **************************************************************************************************

!        Real part
         sum_psi1H1psi1 = sum_psi1H1psi1 + dtset%wtk(ikpt)*occ_k(jband)*dotr
         sum_lambda1psi1psi1 = sum_lambda1psi1psi1 - dtset%wtk(ikpt)*occ_k(jband)*lagr
         sum_lambda1psi0S1psi1 = sum_lambda1psi0S1psi1 - dtset%wtk(ikpt)*occ_k(jband)*lagr_paw

!        Imaginary part
         sum_psi1H1psi1_i = sum_psi1H1psi1_i + dtset%wtk(ikpt)*occ_k(jband)*doti
         sum_lambda1psi1psi1_i = sum_lambda1psi1psi1_i - dtset%wtk(ikpt)*occ_k(jband)*lagi
         sum_lambda1psi0S1psi1_i = sum_lambda1psi0S1psi1_i - dtset%wtk(ikpt)*occ_k(jband)*lagi_paw

! **************************************************************************************************
!        If compute_rho21 : accumulate rhoij and compute term with H_KV^(2)
! **************************************************************************************************

         if (compute_rho21) then

           if (i1pert<=natom) then ! If true, i3pert==natom+2
             cwave_right = cg3(:,1+offset_cgj:size_wf+offset_cgj)
           else if (i3pert<=natom) then  ! If true, i1pert==natom+2
             cwave_right = cg1(:,1+offset_cgj:size_wf+offset_cgj)
           end if
           choice = 2
           cpopt  = 0
           call getcprj(choice,cpopt,cgj,cwaveprj0,&
&           ffnl1,idir_phon,psps%indlmn,gs_hamkq%istwf_k,kg_k,kpg_k,kpt,psps%lmnmax,&
&           mgfft,mpi_enreg,natom,nattyp,dtset%ngfft,dtset%nloalg,&
&           npw_k,nspinor,psps%ntypat,phkxred,ph1d,ph3d,ucvol,psps%useylm)
           call getcprj(choice,cpopt,cwave_right,cwaveprj1,&
&           ffnl1,idir_phon,psps%indlmn,gs_hamkq%istwf_k,kg_k,kpg_k,kpt,psps%lmnmax,&
&           mgfft,mpi_enreg,natom,nattyp,dtset%ngfft,dtset%nloalg,&
&           npw_k,nspinor,psps%ntypat,phkxred,ph1d,ph3d,ucvol,psps%useylm)

           usetimerev=(dtset%kptopt>0.and.dtset%kptopt<3)
           cplex_cprj=2;if (gs_hamkq%istwf_k>1) cplex_cprj=1
           call paw_dfptnl_accrhoij(atindx,cplex_cprj,cwaveprj0,cwaveprj0,cwaveprj1,cwaveprj1,i1pert,i3pert,isppol,natom,natom,&
&            nspinor,occ_k(jband),pawrhoij21_unsym,usetimerev,wtk_k)
!&           comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)

!          Compute < psi^(0) | H_KV^(pert1pert3) | psi^(pert2) > + < psi^(pert2) | H_KV^(pert1pert3) | psi^(0) >
           cwavef2(:,:) = cg2(:,1+offset_cgj:size_wf+offset_cgj)

!          Read dkk file (for tests only)
           if (debug_mode/=0) then
             if(idir_elfd==i2dir) then
               call wfk_read_bks(ddk_f(2), jband, ikpt, isppol, xmpio_single, cg_bks=iddk)
             else
               call wfk_read_bks(ddk_f(4), jband, ikpt, isppol, xmpio_single, cg_bks=iddk)
             end if
             cg_jband(1,1+size_wf*(jband-1):size_wf*jband,2) = -iddk(2,1:size_wf)
             cg_jband(2,1+size_wf*(jband-1):size_wf*jband,2) =  iddk(1,1:size_wf)
           end if
!          Read dkde file
           if(idir_elfd==i2dir) then
             call wfk_read_bks(ddk_f(3), jband, ikpt, isppol, xmpio_single, cg_bks=iddk)
           else
             call wfk_read_bks(ddk_f(5), jband, ikpt, isppol, xmpio_single, cg_bks=iddk)
           end if
           s_cwave = iddk
           iddk(1,:) = -s_cwave(2,:)
           iddk(2,:) =  s_cwave(1,:)
           call rf2_getidir(idir_phon,idir_elfd,idir_getgh2c)
           call rf2_apply_hamiltonian(cg_jband,cprj_jband,cwavef2,cprj_empty,s_cwave,dummy_array2,eig0_k,eig1_k_i2pert,&
&                                jband,gs_hamkq,iddk,idir_getgh2c,ipert_phon+natom+11,ikpt,isppol,mkmem,mpi_enreg,nband_k,nsppol,&
                                 debug_mode,dtset%prtvol,rf_hamkq_i2pert,size_cprj,size_wf,enl=chi_ij,ffnl1=ffnl1,ffnl1_test=ffnl1_test)
           call dotprod_g(enlout1(1),enlout1(2),gs_hamkq%istwf_k,npw_k*nspinor,2,cgj,s_cwave,&
&                 mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
           sum_psi0H2psi1a   = sum_psi0H2psi1a   + half*dtset%wtk(ikpt)*occ_k(jband)*enlout1(1) ! be careful : factor 0.5
           sum_psi0H2psi1a_i = sum_psi0H2psi1a_i + half*dtset%wtk(ikpt)*occ_k(jband)*enlout1(2) ! be careful : factor 0.5

!          Read ddk file
           if(idir_elfd==i2dir) then
             call wfk_read_bks(ddk_f(2), jband, ikpt, isppol, xmpio_single, cg_bks=iddk)
           else
             call wfk_read_bks(ddk_f(4), jband, ikpt, isppol, xmpio_single, cg_bks=iddk)
           end if
           s_cwave = iddk
           iddk(1,:) = -s_cwave(2,:)
           iddk(2,:) =  s_cwave(1,:)
           call rf2_apply_hamiltonian(cg_jband,cprj_jband,cgj,cprj_empty,s_cwave,dummy_array2,eig0_k,eig1_k_i2pert,&
&                                jband,gs_hamkq,iddk,idir_getgh2c,ipert_phon+natom+11,ikpt,isppol,mkmem,mpi_enreg,nband_k,nsppol,&
                                 debug_mode,dtset%prtvol,rf_hamkq_i2pert,size_cprj,size_wf,enl=chi_ij,ffnl1=ffnl1,ffnl1_test=ffnl1_test)
           call dotprod_g(enlout2(1),enlout2(2),gs_hamkq%istwf_k,npw_k*nspinor,2,cwavef2,s_cwave,&
&                 mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
           sum_psi0H2psi1b   = sum_psi0H2psi1b   + half*dtset%wtk(ikpt)*occ_k(jband)*enlout2(1) ! be careful : factor 0.5
           sum_psi0H2psi1b_i = sum_psi0H2psi1b_i + half*dtset%wtk(ikpt)*occ_k(jband)*enlout2(2) ! be careful : factor 0.5

         end if ! end if compute_rho21

       end if
     end do   ! end loop over bands

! **************************************************************************************************
!    END OF BAND LOOP
! **************************************************************************************************

     ABI_DEALLOCATE(cgi)
     ABI_DEALLOCATE(cgj)
     ABI_DEALLOCATE(iddk)
     ABI_DEALLOCATE(h_cwave)
     ABI_DEALLOCATE(s_cwave)

     ABI_DEALLOCATE(cwave_right)
     ABI_DEALLOCATE(cwave_left)
     ABI_DEALLOCATE(eig1_k_i2pert)
     ABI_DEALLOCATE(eig1_k_stored)

     bandtot = bandtot + nband_k
     icg0 = icg0 + npw_k*nspinor*nband_k
     ikg = ikg + npw_k
     ikg1 = ikg1 + npw1_k

     ABI_DEALLOCATE(cwavef1)
     if (compute_rho21) then
       ABI_DEALLOCATE(cwavef2)
     end if
     ABI_DEALLOCATE(cwavef3)
     ABI_DEALLOCATE(dkinpw)
     ABI_DEALLOCATE(gh0)
     ABI_DEALLOCATE(gh1)
     ABI_DEALLOCATE(gvnl)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(kg1_k)
     ABI_DEALLOCATE(kinpw1)
     ABI_DEALLOCATE(dudk)
     ABI_DEALLOCATE(dudkde)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylm1_k)
     ABI_DEALLOCATE(ylmgr1_k)
     ABI_DEALLOCATE(ffnl1)
     if (compute_rho21.and.debug_mode/=0) then
       ABI_DEALLOCATE(ffnl1_test)
     end if
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(kpg1_k)
     ABI_DEALLOCATE(cg_jband)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(ph3d)
     if (has_cprj_jband) call pawcprj_free(cprj_jband)
     ABI_DATATYPE_DEALLOCATE(cprj_jband)

   end do   ! end loop over k-points
   
 end do   ! end loop over spins

 call destroy_rf_hamiltonian(rf_hamkq_i2pert)

! **************************************************************************************************
!    GATHER BAND-BY-BAND AND XC CONTRIBUTIONS
! **************************************************************************************************

 if (xmpi_paral == 1) then

!  Real parts
   buffer(1)  = sum_psi1H1psi1
   buffer(2)  = sum_lambda1psi1psi1
   buffer(3)  = sum_lambda1psi0S1psi1
   buffer(4)  = sum_psi0H2psi1a
   buffer(5)  = sum_psi0H2psi1b

!  Imaginary parts
   buffer(6)  = sum_psi1H1psi1_i
   buffer(7)  = sum_lambda1psi1psi1_i
   buffer(8)  = sum_lambda1psi0S1psi1_i
   buffer(9)  = sum_psi0H2psi1a_i
   buffer(10) = sum_psi0H2psi1b_i

   call xmpi_sum(buffer,spaceComm,ierr)

!  Real parts
   sum_psi1H1psi1          = buffer(1)
   sum_lambda1psi1psi1     = buffer(2)
   sum_lambda1psi0S1psi1   = buffer(3)
   sum_psi0H2psi1a         = buffer(4)
   sum_psi0H2psi1b         = buffer(5)

!  Imaginary parts
   sum_psi1H1psi1_i        = buffer(6)
   sum_lambda1psi1psi1_i   = buffer(7)
   sum_lambda1psi0S1psi1_i = buffer(8)
   sum_psi0H2psi1a_i       = buffer(9)
   sum_psi0H2psi1b_i       = buffer(10)

!  Accumulate PAW occupancies
   if (compute_rho21) then
     call pawrhoij_mpisum_unpacked(pawrhoij21_unsym,spaceComm)
   end if

 end if

! **************************************************************************************************
!      Compute E_xc^(3) (NOTE : E_H^(3) = 0)
! **************************************************************************************************

!      Compute the third-order xc energy
!      take into account the contribution of the term
!$
!      \frac{d}{d \lambda}
!      \frac{\delta^2 E_{Hxc}}{\delta n(r) \delta n(r\prim)}
!$
!      (seventh term of Eq. (110) of X. Gonze, PRA 52, 1096 (1995)).

!      the following are essentially the 4th and the 3rd terms of PRB 71,125107, but the
!      multiplication for rho1 will be done by dotprod_vn later

!!     in the non spin polarized case xc_tmp has only 1 component

 ABI_ALLOCATE(xc_tmp,(cplex*nfftf,nspden))
 ABI_ALLOCATE(rho1r1_tot,(cplex*nfftf,nspden))
 if (nspden==1)then

!   if (cplex==1) then
!!    This, and the next lines, have to be changed in case cplex=2
!     do ifft=1,nfftf
!       xc_tmp(ifft,1)= k3xc(ifft,1)*(rho2r1(ifft,1)+3*xccc3d2(ifft))*rho3r1(ifft,1)
!     end do
!   else
!     do ifft=1,nfftf   ! 2*ifft-1 denotes the real part, 2*ifft the imaginary part
!       xc_tmp(2*ifft-1,1)= k3xc(ifft,1)*( (rho2r1(2*ifft-1,1)+3*xccc3d2(2*ifft-1))*rho3r1(2*ifft-1,1) &
!&      -( rho2r1(2*ifft,1)+3*xccc3d2(2*ifft))*rho3r1(2*ifft,1))
!       xc_tmp(2*ifft,1)= k3xc(ifft,1)*( (rho2r1(2*ifft-1,1)+3*xccc3d2(2*ifft-1))*rho3r1(2*ifft,1) &
!&      +( rho2r1(2*ifft,1)+3*xccc3d2(2*ifft))*rho3r1(2*ifft-1,1))
!     end do
!   end if

   if (cplex==1) then
     do ifft=1,nfftf
       rho1r1_tot(ifft,1) = rho1r1(ifft,1) + xccc3d1(ifft)
       rho2r_re = rho2r1(ifft,1)+xccc3d2(ifft)
       rho3r_re = rho3r1(ifft,1)+xccc3d3(ifft)
       xc_tmp(ifft,1)= k3xc(ifft,1)*rho2r_re*rho3r_re
     end do
   else
     do ifft=1,nfftf   ! 2*ifft-1 denotes the real part, 2*ifft the imaginary part
       rho1r1_tot(2*ifft-1,1) = rho1r1(2*ifft-1,1) + xccc3d1(2*ifft-1)
       rho1r1_tot(2*ifft  ,1) = rho1r1(2*ifft  ,1) + xccc3d1(2*ifft  )
       rho2r_re = rho2r1(2*ifft-1,1)+xccc3d2(2*ifft-1)
       rho2r_im = rho2r1(2*ifft  ,1)+xccc3d2(2*ifft  )
       rho3r_re = rho3r1(2*ifft-1,1)+xccc3d3(2*ifft-1)
       rho3r_im = rho3r1(2*ifft  ,1)+xccc3d3(2*ifft  )
       xc_tmp(2*ifft-1,1)= k3xc(ifft,1)*(rho2r_re*rho3r_re-rho2r_im*rho3r_im)
       xc_tmp(2*ifft  ,1)= k3xc(ifft,1)*(rho2r_re*rho3r_im+rho2r_im*rho3r_re)
     end do
   end if

 else
   MSG_BUG('NONLINEAR with nspden==2 is not implemented yet')
 end if

!!                  fab: modifications for the spin polarized raman part:
!!                  in the spin polarized case xc_tmp has 2 components
!!                  note that now the non linear core correction is divided by 2
!                   if (nspden==2) then

!                     ABI_ALLOCATE(xc_tmp,(cplex*nfft,2))

!                     if (cplex==1) then
!                       do ifft=1,nfft
!                         xc_tmp(ifft,1)= k3xc(ifft,1)*(rho2r1(ifft,2)+(3._dp/2._dp)*xccc3d2(ifft))*rho3r1(ifft,2)+ &
!&                         k3xc(ifft,2)*(rho2r1(ifft,2)+(3._dp/2._dp)*xccc3d2(ifft))*(rho3r1(ifft,1)-rho3r1(ifft,2))+ &
!&                         k3xc(ifft,2)*((rho2r1(ifft,1)-rho2r1(ifft,2))+(3._dp/2._dp)*xccc3d2(ifft))*rho3r1(ifft,2)+ &
!&                         k3xc(ifft,3)*((rho2r1(ifft,1)-rho2r1(ifft,2))+(3._dp/2._dp)*xccc3d2(ifft))*(rho3r1(ifft,1)-rho3r1(ifft,2))
!                         xc_tmp(ifft,2)= k3xc(ifft,2)*(rho2r1(ifft,2)+(3._dp/2._dp)*xccc3d2(ifft))*rho3r1(ifft,2)+ &
!&                         k3xc(ifft,3)*(rho2r1(ifft,2)+(3._dp/2._dp)*xccc3d2(ifft))*(rho3r1(ifft,1)-rho3r1(ifft,2))+ &
!&                         k3xc(ifft,3)*((rho2r1(ifft,1)-rho2r1(ifft,2))+(3._dp/2._dp)*xccc3d2(ifft))*rho3r1(ifft,2)+ &
!&                         k3xc(ifft,4)*((rho2r1(ifft,1)-rho2r1(ifft,2))+(3._dp/2._dp)*xccc3d2(ifft))*(rho3r1(ifft,1)-rho3r1(ifft,2))
!                       end do

!                     else
!                       do ifft=1,nfft
!!                        These sections should be rewritten, to be easier to read ... (defining intermediate scalars)
!                         xc_tmp(2*ifft-1,1)= k3xc(ifft,1)*&
!&                         ( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*rho3r1(2*ifft-1,2)- &
!&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft,2))+   &
!&                         k3xc(ifft,2)*&
!&                         ( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*(rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2))- &
!&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*(rho3r1(2*ifft,1)-rho3r1(2*ifft,2)))+ &
!&                         k3xc(ifft,2)*&
!&                         ( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*rho3r1(2*ifft-1,2)- &
!&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft,2))+ &
!&                         k3xc(ifft,3)*&
!&                         ( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
!&                         (rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2))- &
!&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*&
!&                         (rho3r1(2*ifft,1)-rho3r1(2*ifft,2)))
!                         xc_tmp(2*ifft,1)=k3xc(ifft,1)*&
!&                         ( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*rho3r1(2*ifft,2)+ &
!&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft-1,2))+   &
!&                         k3xc(ifft,2)*&
!&                         ( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*(rho3r1(2*ifft,1)-rho3r1(2*ifft,2))+ &
!&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*(rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2)))+ &
!&                         k3xc(ifft,2)*&
!&                         ( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*rho3r1(2*ifft,2)+ &
!&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft-1,2))+ &
!&                         k3xc(ifft,3)*&
!&                         ( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
!&                         (rho3r1(2*ifft,1)-rho3r1(2*ifft,2))+ &
!&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*&
!&                         (rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2)))
!!                        fab: now the spin down component
!                         xc_tmp(2*ifft-1,2)= k3xc(ifft,2)*&
!&                         ( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*rho3r1(2*ifft-1,2)- &
!&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft,2))+   &
!&                         k3xc(ifft,3)*( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
!&                         (rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2))- &
!&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*(rho3r1(2*ifft,1)-rho3r1(2*ifft,2)))+ &
!&                         k3xc(ifft,3)*( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
!&                         rho3r1(2*ifft-1,2)- &
!&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft,2))+ &
!&                         k3xc(ifft,4)*( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
!&                         (rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2))- &
!                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*&
!&                         (rho3r1(2*ifft,1)-rho3r1(2*ifft,2)))
!                         xc_tmp(2*ifft,2)=k3xc(ifft,1)*( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
!&                         rho3r1(2*ifft,2)+ &
!&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft-1,2))+   &
!&                         k3xc(ifft,3)*( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
!&                         (rho3r1(2*ifft,1)-rho3r1(2*ifft,2))+ &
!&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*(rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2)))+ &
!&                         k3xc(ifft,3)*( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
!&                         rho3r1(2*ifft,2)+ &
!&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft-1,2))+ &
!&                         k3xc(ifft,4)*( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
!&                         (rho3r1(2*ifft,1)-rho3r1(2*ifft,2))+ &
!&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*&
!&                         (rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2)))
!                       end do

!!                      fab: this is the end if over cplex
!                     end if
!!                    fab: this is the enf if over nspden
!                   end if

 
 call dotprod_vn(cplex,rho1r1_tot,exc3(1),exc3(2),nfftf,nfftotf,nspden,2,xc_tmp,ucvol,mpi_comm_sphgrid=mpi_enreg%comm_fft)
 ABI_DEALLOCATE(xc_tmp)
 ABI_DEALLOCATE(rho1r1_tot)

 exc3_paw = zero
 if (usepaw==1) then

   call paw_dfptnl_energy(exc3_paw,dtset%ixc,natom,natom,psps%ntypat,paw_an0,pawang,dtset%pawprtvol,pawrad,&
&   pawrhoij1_i1pert,pawrhoij1_i2pert,pawrhoij1_i3pert,pawtab,dtset%pawxcdev,mpi_enreg%my_atmtab,mpi_enreg%comm_atom)

 end if

! **************************************************************************************************
!      Compute E_Hxc^(2:1)
! **************************************************************************************************

 eHxc21_paw = zero
 eHxc21_nhat = zero
 if (compute_rho21) then

   if (pawfgr%nfft/=nfftf) then
     write(msg,'(2(a,i10))') 'pawfgr%nfft/=nfftf : pawfgr%nfft=',pawfgr%nfft,' nfftf = ',nfftf
     MSG_ERROR(msg)
   end if
   ABI_ALLOCATE(nhat21,(cplex*nfftf,nspden))
   call pawnhatfr(0,idir_phon,ipert_phon,natom,dtset%natom,nspden,psps%ntypat,&
&   pawang,pawfgrtab,pawrhoij11,pawtab,rprimd)

   call pawmkrho(0,arg,cplex,gs_hamkq%gprimd,idir_phon,indsy1,ipert_phon,mpi_enreg,&
&   natom,natom,nspden,nsym1,psps%ntypat,dtset%paral_kgb,pawang,pawfgr,pawfgrtab,&
&   dtset%pawprtvol,pawrhoij21,pawrhoij21_unsym,pawtab,dtset%qptn,dummy_array2,dummy_array2,&
&   dummy_array2,rprimd,symaf1,symrc1,dtset%typat,ucvol,dtset%usewvl,xred,&
&   pawang_sym=pawang1,pawnhat=nhat21,pawrhoij0=pawrhoij0)

!   if (paral_atom) then
!     call pawrhoij_free(pawrhoij21_unsym)
!     ABI_DATATYPE_DEALLOCATE(pawrhoij21_unsym)
!   end if
   nzlmopt = 0
   call pawdfptenergy(eHxc21_paw,i2pert,ipert_phon,dtset%ixc,natom,dtset%natom,dtset%ntypat,&
&   nzlmopt,nzlmopt,paw_an0,paw_an1_i2pert,paw_ij1_i2pert,pawang,dtset%pawprtvol,&
&   pawrad,pawrhoij1_i2pert,pawrhoij21,pawtab,dtset%pawxcdev,dtset%xclevel)
!&   mpi_atmtab=my_atmtab,comm_atom=my_comm_atom

   ABI_ALLOCATE(v_i2pert,(cplex*nfftf,nspden))
   v_i2pert(:,1) = vhartr1_i2pert(:)
   if(nspden>1) then
     v_i2pert(:,2) = vhartr1_i2pert(:)
   end if
   if (debug_mode/=0) then
     write(msg,'(2(a,i6))') ' DFPTNL_PERT : pawtab(',ipert_phon,')%usexcnhat = ',pawtab(ipert_phon)%usexcnhat
     call wrtout(std_out,msg,'COLL')
   end if
   if (pawtab(ipert_phon)%usexcnhat>0) then
     v_i2pert(:,:) = v_i2pert(:,:) + vxc1_i2pert(:,:)
   end if

   call dotprod_vn(cplex,nhat21,eHxc21_nhat(1),eHxc21_nhat(2),nfftf,nfftotf,nspden,2,v_i2pert,ucvol,mpi_comm_sphgrid=mpi_enreg%comm_fft)

   ABI_DEALLOCATE(v_i2pert)
   ABI_DEALLOCATE(nhat21)

 end if

! **************************************************************************************************
!    ALL TERMS HAVE BEEN COMPUTED
! **************************************************************************************************

!Real part
 e3tot = sum_psi1H1psi1 + sum_lambda1psi1psi1 + sum_lambda1psi0S1psi1 + sum_psi0H2psi1a + sum_psi0H2psi1b
 e3tot = e3tot + half * (eHxc21_paw(1)+eHxc21_nhat(1)) + sixth * (exc3(1) + exc3_paw(1))

!Imaginary part
 sumi = sum_psi1H1psi1_i + sum_lambda1psi1psi1_i + sum_lambda1psi0S1psi1_i + sum_psi0H2psi1a_i + sum_psi0H2psi1b_i
 sumi = sumi   + half * (eHxc21_paw(2)+eHxc21_nhat(2)) + sixth * (exc3(2) + exc3_paw(2))

 write(msg,'(2a,3(a,i2,a,i1))') ch10,'NONLINEAR : ',&
 ' perts : ',i1pert,'.',i1dir,' / ',i2pert,'.',i2dir,' / ',i3pert,'.',i3dir
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 if (dtset%nonlinear_info>0) then
   write(msg,'(10(a,2(a,f18.8)),a)') &
    ch10,'        sum_psi1H1psi1 = ',sum_psi1H1psi1,        ',',sum_psi1H1psi1_i,&
    ch10,'   sum_lambda1psi1psi1 = ',sum_lambda1psi1psi1,   ',',sum_lambda1psi1psi1_i,&
    ch10,' sum_lambda1psi0S1psi1 = ',sum_lambda1psi0S1psi1, ',',sum_lambda1psi0S1psi1_i,&
    ch10,'       sum_psi0H2psi1a = ',sum_psi0H2psi1a,       ',',sum_psi0H2psi1a_i,&
    ch10,'       sum_psi0H2psi1b = ',sum_psi0H2psi1b,       ',',sum_psi0H2psi1b_i,&
    ch10,'                exc3/6 = ',sixth*exc3(1),         ',',sixth*exc3(2),&
    ch10,'            exc3_paw/6 = ',sixth*exc3_paw(1),     ',',sixth*exc3_paw(2),&
    ch10,'          eHxc21_paw/2 = ',half*eHxc21_paw(1),    ',',half*eHxc21_paw(2),&
    ch10,'         eHxc21_nhat/2 = ',half*eHxc21_nhat(1),   ',',half*eHxc21_nhat(2),&
    ch10,' >>>>>>>>>>>>>>> e3tot = ',e3tot,                 ',',sumi,ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
  end if

!Real parts
 d3etot(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)   = e3tot
 d3etot_1(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sum_psi1H1psi1
 d3etot_2(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sum_lambda1psi1psi1
 d3etot_3(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sum_lambda1psi0S1psi1
 d3etot_4(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sum_psi0H2psi1a
 d3etot_5(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sum_psi0H2psi1b
 d3etot_6(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = half * eHxc21_paw(1)
 d3etot_7(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = half * eHxc21_nhat(1)
 d3etot_8(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sixth * exc3(1)
 d3etot_9(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sixth * exc3_paw(1)
!Imaginary parts
 d3etot(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)   = sumi
 d3etot_1(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sum_psi1H1psi1_i
 d3etot_2(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sum_lambda1psi1psi1_i
 d3etot_3(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sum_lambda1psi0S1psi1_i
 d3etot_4(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sum_psi0H2psi1a_i
 d3etot_5(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sum_psi0H2psi1b_i
 d3etot_6(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = half * eHxc21_paw(2)
 d3etot_7(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = half * eHxc21_nhat(2)
 d3etot_8(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sixth * exc3(2)
 d3etot_9(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sixth * exc3_paw(2)

 if (compute_rho21) then
   call pawcprj_free(cwaveprj0)
   call pawcprj_free(cwaveprj1)
   call pawrhoij_free(pawrhoij21)
   ABI_DEALLOCATE(chi_ij)
 end if
 ABI_DEALLOCATE(phkxred)
 ABI_DATATYPE_DEALLOCATE(cwaveprj0)
 ABI_DATATYPE_DEALLOCATE(cwaveprj1)
 ABI_DATATYPE_DEALLOCATE(pawrhoij21)
 ABI_DEALLOCATE(ylm)
 ABI_DEALLOCATE(ylm1)
 ABI_DEALLOCATE(ylmgr)
 ABI_DEALLOCATE(ylmgr1)
 ABI_DEALLOCATE(vlocal1_i2pert)
 ABI_DEALLOCATE(wfraug)

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 DBG_EXIT("COLL")

end subroutine dfptnl_pert
!!***
